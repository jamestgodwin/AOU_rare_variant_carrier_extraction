#!/usr/bin/env python3
"""
Carrier extraction for regenie rare-variant aggregates on All of Us. Run it through
run_carrier_extraction.sh; all settings are in the JSON config (see config.json / README.md).

Steps (each one is skipped when its output already exists, unless listed in "redo"):
  merge    per-chromosome regenie results for one run/ancestry/analysis -> merged .tsv.gz
  index    local plain copy of the merged results + a light lookup index (seconds, once)
  samples  phenotype + sample files (people analysed for this ancestry)
  select   pick aggregates (targets / masks / units / LOG10P) -> aggregates_{PREFIX}.csv
  extract  carriers per aggregate

Why it's fast:
  * selecting aggregates = grep over a slim file of aggregate IDs (+ numpy for LOG10P), no database
  * annotation lookup = zstd -dc | grep -F once per gene/region, then cached
  * genotypes: every variant extracted so far is kept in a local store per chromosome with its
    allele counts precomputed; plink2 only runs for variants the store doesn't have yet, reads the
    .pgen in place through a bucket mount and parses a slimmed .pvar
  * plink 1.9 is replaced by numpy on the packed genotype bytes
  * bucket files are read and written through gcsfuse mounts; gcloud is only the fallback
"""

import argparse
import copy
import fcntl
import glob
import hashlib
import io
import json
import os
import pickle
import re
import shutil
import sqlite3
import subprocess
import sys
import threading
import time
import zlib
from concurrent.futures import ThreadPoolExecutor
from functools import cached_property, lru_cache

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("ERROR: python3 needs numpy and pandas (pip install numpy pandas)")

T0 = time.time()
_UPLOADS = ThreadPoolExecutor(max_workers=2)   # uploads that can run while other work continues
_PENDING = []
TIMES = {}


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def timed(name, t):
    TIMES[name] = TIMES.get(name, 0.0) + time.time() - t


# =====================================================================================
# config
# =====================================================================================
DEFAULTS = {
    "run": None, "ancestry": None, "analysis": "burden",
    "selection": {"name": "top_hits", "targets": [], "targets_file": None, "masks": [], "units": [],
                  "log10p_min": 6.0, "max_aggregates": 500},
    "max_maf": 0.001,
    "inputs": {"results_dir": None, "pheno_file": None, "cov_file": None, "ancestry_file": None,
               "pheno_col": None, "anno_file": None, "mask_file": None, "geno_prefix": None,
               "geno_format": "pgen", "bgen_sample": None, "chr_x": "X", "range_units": ["sliding"]},
    "outputs": {"out_dir": None, "cache_dir": "~/carrier_cache"},
    "performance": {"geno_access": "auto", "use_mounts": True, "parallel_chromosomes": 1,
                    "plink2_threads": 0, "plink2_memory_mb": 0, "genotype_cache_gb": 50,
                    "upload_threads": 16},
    "redo": [],
}


def deep_update(base, new):
    for k, v in new.items():
        if k.startswith("_"):
            continue
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            deep_update(base[k], v)
        else:
            base[k] = v
    return base


class Cfg:
    def __init__(self, path, overrides=None):
        c = deep_update(copy.deepcopy(DEFAULTS), json.load(open(path)))
        if overrides:
            deep_update(c, overrides)
        for key in ("run", "ancestry"):
            if not c[key]:
                sys.exit(f"config: '{key}' is required")
        for key in ("results_dir", "anno_file", "mask_file", "geno_prefix"):
            if not c["inputs"][key]:
                sys.exit(f"config: inputs.{key} is required")
        if not c["outputs"]["out_dir"]:
            sys.exit("config: outputs.out_dir is required")
        self.c = c
        self.sel = c["selection"]
        self.inp = c["inputs"]
        self.perf = c["performance"]
        self.base = f"{c['run']}_{c['ancestry']}_{c['analysis']}"
        self.prefix = f"{self.base}_{self.sel['name']}"
        out = c["outputs"]["out_dir"].rstrip("/")
        self.out_dir = out if out.startswith("gs://") else os.path.expanduser(out)
        self.cache = os.path.expanduser(c["outputs"]["cache_dir"])
        self.merged = f"{self.out_dir}/{self.base}_aggregates_merged.tsv.gz"
        # local results index (see step_index)
        self.results_local = os.path.join(self.cache, "results")
        stem = os.path.join(self.results_local, f"{self.base}_aggregates_merged")
        self.plain = stem + ".tsv"
        self.keys = stem + ".keys"
        self.offsets = stem + ".offsets.npy"
        self.log10p = stem + ".log10p.npy"
        self.chroms = stem + ".chr.npy"
        self.stamp = stem + ".index.json"
        self.old_db = os.path.join(self.results_local, f"{self.base}.sqlite")      # earlier version
        self.meta_out = f"{self.out_dir}/metadata_{c['ancestry']}"
        self.meta_local = os.path.join(self.cache, f"metadata_{c['ancestry']}")
        self.agg_file = f"{self.out_dir}/input/aggregates_{self.prefix}.csv"
        self.carriers_out = f"{self.out_dir}/output/{self.prefix}"
        self.work = os.path.join(self.cache, "work", self.prefix)
        self.redo = set(c["redo"])

    def chr_path(self, template, CHR):
        label = self.inp["chr_x"] if str(CHR) == "X" else str(CHR)
        return os.path.expanduser(template.replace("{CHR}", label))


# =====================================================================================
# storage: gs:// paths are read and written through gcsfuse mounts when possible
# (a mount the Workbench made, or one made here under cache_dir/mnt/); gcloud otherwise
# =====================================================================================
MNT = {"cache": None, "use": True, "out_bucket": None, "threads": 16}
MOUNTS_FILE = os.environ.get("CARRIER_MOUNTS_FILE", "/proc/mounts")
BUCKET_RE = re.compile(r"^[a-z0-9][a-z0-9._-]{1,220}[a-z0-9]$")


def is_gs(p):
    return str(p).startswith("gs://")


def split_gs(p):
    bucket, _, key = p[5:].partition("/")
    return bucket, key.strip("/")


def _unescape(s):
    return re.sub(r"\\([0-7]{3})", lambda m: chr(int(m.group(1), 8)), s)


def gcs_mounts():
    """(bucket, mount point, read-only) for every gcsfuse mount on this machine."""
    out = []
    try:
        with open(MOUNTS_FILE) as f:
            for line in f:
                p = line.split()
                if len(p) >= 4 and "gcsfuse" in p[2]:
                    out.append((p[0], _unescape(p[1]), "ro" in p[3].split(",")))
    except OSError:
        pass
    return out


def own_mount(bucket, rw):
    return os.path.join(MNT["cache"], "mnt", "rw", bucket) if rw else os.path.join(MNT["cache"], "mnt", bucket)


def _alive(path):
    try:
        os.stat(path)
        return True
    except OSError:
        return False


def mounted_roots(bucket):
    """[(mount point, made here, writable)] of live gcsfuse mounts of this bucket; ours first."""
    if not MNT["use"]:
        return []
    own = {own_mount(bucket, False), own_mount(bucket, True)} if MNT["cache"] else set()
    roots = [(mp, mp in own, not ro) for dev, mp, ro in gcs_mounts() if dev == bucket and _alive(mp)]
    return sorted(roots, key=lambda r: not r[1])


_READ_CACHE = {}
_mount_lock = threading.Lock()
_mount_tried = {}


def _workspace_search(bucket, key):
    """Files under ~/workspace mounts that aren't identified by bucket name: match by path tail."""
    others = [mp for dev, mp, _ in gcs_mounts() if dev != bucket and BUCKET_RE.match(dev)]
    parts = key.split("/")
    roots = glob.glob(os.path.expanduser("~/workspace/*")) + glob.glob(os.path.expanduser("~/workspace/*/*"))
    for root in roots:
        if any(root == mp or root.startswith(mp + "/") for mp in others) or not os.path.isdir(root):
            continue
        for k in range(len(parts)):
            cand = os.path.join(root, *parts[k:])
            if os.path.exists(cand):
                return cand
    return None


def read_path(gs_path, auto=False):
    """Local path of a gs:// file through a bucket mount (None if there isn't one).
    auto=True mounts the bucket with gcsfuse if needed (read-only, except the output bucket)."""
    if not MNT["use"]:
        return None
    if gs_path in _READ_CACHE:
        return _READ_CACHE[gs_path]
    bucket, key = split_gs(gs_path)
    found = None
    for root, own, _ in mounted_roots(bucket):
        cand = os.path.join(root, key)
        if os.path.exists(cand):
            found = cand
            break
        if not own:     # a Workbench mount of a folder inside the bucket: match on the path's tail
            parts = key.split("/")
            for k in range(1, len(parts)):
                cand = os.path.join(root, *parts[k:])
                if os.path.exists(cand):
                    found = cand
                    break
            if found:
                break
    if found is None:
        found = _workspace_search(bucket, key)
    if found is None and auto:
        rw = bucket == MNT["out_bucket"]
        if gcsfuse_mount(bucket, rw):
            cand = os.path.join(own_mount(bucket, rw), key)
            found = cand if os.path.exists(cand) else None
    if found:
        _READ_CACHE[gs_path] = found
    return found


def gcsfuse_mount(bucket, rw):
    """Mount gs://bucket under cache_dir/mnt/ (read-only unless rw). The mount stays up for later runs."""
    if not MNT["use"] or not MNT["cache"]:
        return False
    mp = own_mount(bucket, rw)
    with _mount_lock:
        listed = [m for _, m, _ in gcs_mounts() if m == mp]
        if listed and _alive(mp):
            return True
        if (bucket, rw) in _mount_tried:
            return _mount_tried[(bucket, rw)]
        ok = False
        if shutil.which("gcsfuse"):
            if listed:                                           # stale mount left behind
                subprocess.run(["fusermount", "-uz", mp], capture_output=True)
            os.makedirs(mp, exist_ok=True)
            cmd = ["gcsfuse", "--implicit-dirs"] + ([] if rw else ["-o", "ro"]) + [bucket, mp]
            t = time.time()
            try:
                r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                ok = r.returncode == 0
                err = r.stderr.strip()[-300:]
            except subprocess.TimeoutExpired:
                err = "timed out"
            if ok:
                log(f"   mounted gs://{bucket} at {mp} ({'read-write' if rw else 'read-only'}, "
                    f"{time.time() - t:.1f}s; it stays mounted for later runs)")
            else:
                log(f"   gcsfuse mount of {bucket} failed ({err}); using gcloud instead")
        _mount_tried[(bucket, rw)] = ok
        return ok


def write_dir(gs_dir):
    """Local folder that writes into gs://bucket/key (our read-write mount of the bucket), or None."""
    bucket, key = split_gs(gs_dir)
    if not gcsfuse_mount(bucket, rw=True):
        return None
    d = os.path.join(own_mount(bucket, True), key)
    try:
        os.makedirs(d, exist_ok=True)
        return d
    except OSError:
        return None


def gs_ls(path, long=False):
    cmd = ["gcloud", "storage", "ls"] + (["-l"] if long else []) + [path]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.stdout.splitlines() if r.returncode == 0 else []


def exists(path):
    if not is_gs(path):
        return os.path.exists(path)
    bucket, key = split_gs(path)
    for root, own, _ in mounted_roots(bucket):
        if own:                                  # our own mount shows the whole bucket: conclusive
            return os.path.exists(os.path.join(root, key))
    if read_path(path):
        return True
    return bool(gs_ls(path))


def gs_size(path):
    m = read_path(path)
    if m:
        return os.path.getsize(m)
    for line in gs_ls(path, long=True):
        parts = line.split()
        if len(parts) >= 3 and parts[-1] == path and parts[0].isdigit():
            return int(parts[0])
    return None


def copy_files(srcs, dst_dir):
    """Copy files to a folder with gcloud (local or gs://), many at a time."""
    srcs = list(srcs)
    if not srcs:
        return
    if is_gs(dst_dir) or any(is_gs(s) for s in srcs):
        for k in range(0, len(srcs), 400):
            r = subprocess.run(["gcloud", "storage", "cp", "--quiet", *srcs[k:k + 400], dst_dir.rstrip("/") + "/"],
                               capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"gcloud storage cp failed:\n{r.stderr[-2000:]}")
    else:
        os.makedirs(dst_dir, exist_ok=True)
        for s in srcs:
            shutil.copy(s, dst_dir)


def fetch(path, folder, refresh=False):
    """Local copy of a file (gs:// is copied once into folder and reused); local paths returned as is."""
    path = os.path.expanduser(path)
    if not is_gs(path):
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        return path
    local = os.path.join(folder, os.path.basename(path))
    if refresh or not os.path.exists(local):
        os.makedirs(folder, exist_ok=True)
        tmp = local + ".part"
        mounted = read_path(path)
        if mounted:
            shutil.copy2(mounted, tmp)            # keeps the bucket file's time stamp
        else:
            r = subprocess.run(["gcloud", "storage", "cp", "--quiet", path, tmp], capture_output=True, text=True)
            if r.returncode != 0:
                raise FileNotFoundError(f"{path}: {r.stderr.strip()[-500:]}")
        os.replace(tmp, local)
    return local


def publish(local_files, dest_dir):
    """Write outputs to their destination folder: through the bucket mount (many files at once),
    gcloud for anything the mount couldn't take."""
    files = [f for f in local_files if f]
    if not files:
        return
    if not is_gs(dest_dir):
        os.makedirs(dest_dir, exist_ok=True)
        for f in files:
            if os.path.dirname(os.path.abspath(f)) != os.path.abspath(dest_dir):
                shutil.copy(f, dest_dir)
        return
    failed = files
    d = write_dir(dest_dir)
    if d:
        def put(f):
            try:
                shutil.copyfile(f, os.path.join(d, os.path.basename(f)))
                return None
            except OSError:
                return f
        with ThreadPoolExecutor(max_workers=max(1, int(MNT["threads"]))) as ex:
            failed = [f for f in ex.map(put, files) if f]
        if failed:
            log(f"   {len(failed)} files could not be written through the mount; using gcloud for them")
    copy_files(failed, dest_dir)


def decompress_cmd(local_path):
    """argv that prints a (possibly .zst/.gz) local file to stdout."""
    if local_path.endswith(".zst"):
        return ["zstd", "-dc", local_path] if shutil.which("zstd") else ["plink2", "--zst-decompress", local_path]
    if local_path.endswith(".gz"):
        return ["pigz" if shutil.which("pigz") else "gzip", "-dc", local_path]
    return ["cat", local_path]


def decompress_stdin_cmd(name):
    """argv that decompresses stdin (None if that isn't possible for this file type)."""
    if name.endswith(".zst"):
        return ["zstd", "-dc"] if shutil.which("zstd") else None
    if name.endswith(".gz"):
        return ["gzip", "-dc"]
    return ["cat"]


def read_text_table(local_path, **kwargs):
    if local_path.endswith(".zst"):
        p = subprocess.Popen(decompress_cmd(local_path), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, text=True)
        df = pd.read_csv(p.stdout, **kwargs)
        p.stdout.close()
        if p.wait() != 0:
            raise RuntimeError(f"could not decompress {local_path}")
        return df
    return pd.read_csv(local_path, **kwargs)


def chrom_order(c):
    c = str(c)
    return int(c) if c.isdigit() else {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(c, 99)


VARIANT_RE = re.compile(r"^(?:.*[:_])?(?:chr)?([0-9]{1,2}|X|Y|M|MT)[:_-]([0-9]+)[:_-]([ACGTN*]+)[:_-]([ACGTN*]+)$",
                        re.IGNORECASE)


@lru_cache(maxsize=None)
def variant_key(v):
    """'1:10048:C:CT' (or chr1_10048_C_CT, ...) -> (10048, 'C', 'CT'); None if not parseable."""
    m = VARIANT_RE.match(str(v).strip())
    if not m:
        return None
    return int(m.group(2)), m.group(3).upper(), m.group(4).upper()


def save_npy(path, a):
    with open(path + ".part", "wb") as f:
        np.save(f, a)
    os.replace(path + ".part", path)


def save_text(path, text):
    with open(path + ".part", "w") as f:
        f.write(text)
    os.replace(path + ".part", path)


class file_lock:
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        self.f = open(self.path, "w")
        fcntl.flock(self.f, fcntl.LOCK_EX)
        return self

    def __exit__(self, *exc):
        fcntl.flock(self.f, fcntl.LOCK_UN)
        self.f.close()


# =====================================================================================
# step: merge regenie results
# =====================================================================================
OUT_COLS = ["CHR", "POS", "UNIT.MASK.MAX_FREQ", "REF_ALLELE", "MASK.MAX_FREQ", "freq", "sample_size",
            "TEST", "BETA", "SE", "CHISQ", "LOG10P", "DF", "MAC"]


def parse_extra(extra):
    df_val, mac = pd.NA, pd.NA
    if isinstance(extra, str) and extra != "NA":
        for item in re.split(r"[;,]", extra):
            key, _, val = item.partition("=")
            if key == "DF":
                df_val = f"DF={val}"
            elif key == "MAC":
                mac = pd.to_numeric(val, errors="coerce")
    return df_val, mac


def read_regenie(path):
    df = pd.read_csv(path, sep=r"\s+", comment="#", compression="gzip",
                     dtype={"CHROM": str, "ID": str, "ALLELE0": str, "ALLELE1": str, "TEST": str, "EXTRA": str},
                     na_values=["NA"], keep_default_na=False)
    needed = ["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    if "EXTRA" not in df.columns:
        df["EXTRA"] = pd.NA
    extra = df["EXTRA"].map(parse_extra)
    out = pd.DataFrame({
        "CHR": df["CHROM"], "POS": df["GENPOS"], "UNIT.MASK.MAX_FREQ": df["ID"], "REF_ALLELE": df["ALLELE0"],
        "MASK.MAX_FREQ": df["ALLELE1"], "freq": df["A1FREQ"], "sample_size": df["N"].astype("Int64"),
        "TEST": df["TEST"], "BETA": df["BETA"], "SE": df["SE"], "CHISQ": df["CHISQ"], "LOG10P": df["LOG10P"],
        "DF": extra.map(lambda x: x[0]),
        "MAC": pd.to_numeric(extra.map(lambda x: x[1]), errors="coerce").astype("Int64"),
    })
    # regenie only writes MAC= with --htp; for a burden mask 2 * N * A1FREQ = mask allele count
    out["MAC"] = out["MAC"].fillna((2 * out["sample_size"] * out["freq"]).round().astype("Int64"))
    return out[OUT_COLS]


def step_merge(cfg):
    c = cfg.c
    run, anc, analysis = c["run"], c["ancestry"], c["analysis"]
    rdir = os.path.expanduser(cfg.inp["results_dir"]).rstrip("/")
    names = ([l.strip() for l in gs_ls(rdir + "/") if not l.endswith("/")] if is_gs(rdir)
             else sorted(glob.glob(os.path.join(rdir, "*"))))
    pattern = re.compile(rf"^{re.escape(run)}_(?P<anc>.+?)_chr(?P<chr>[0-9]+|X|Y)_(?P<start>[0-9]+)_[0-9]+_"
                         rf".*_(?P<analysis>[^_]+)_{re.escape(run)}\.regenie\.gz$")
    found, ancs, analyses = [], set(), set()
    for p in names:
        m = pattern.match(os.path.basename(p))
        if not m:
            continue
        ancs.add(m.group("anc"))
        if m.group("anc") == anc:
            analyses.add(m.group("analysis"))
            if m.group("analysis") == analysis:
                found.append((chrom_order(m.group("chr")), int(m.group("start")), m.group("chr"), p))
    log(f"merge: run '{run}': ancestries {sorted(ancs)}; analyses for {anc}: {sorted(analyses)}")
    if not found:
        sys.exit(f"no {run}_{anc}_chr*_{analysis}_{run}.regenie.gz files in {rdir}")
    found.sort()
    missing = {str(k) for k in range(1, 23)} - {ch for _, _, ch, _ in found}
    if missing:
        log(f"WARNING: no results for chromosome(s) {sorted(missing, key=int)}")
    files = [p for *_, p in found]
    tmp = os.path.join(cfg.cache, "tmp", "regenie")
    if is_gs(rdir):
        os.makedirs(tmp, exist_ok=True)
        log(f"merge: downloading {len(files)} result files")
        copy_files(files, tmp)
        files = [os.path.join(tmp, os.path.basename(f)) for f in files]
    with ThreadPoolExecutor(max_workers=min(8, os.cpu_count() or 1)) as ex:
        merged = pd.concat(list(ex.map(read_regenie, files)), ignore_index=True)
    merged["_c"] = merged["CHR"].map(chrom_order)
    merged = merged.sort_values(["_c", "POS", "UNIT.MASK.MAX_FREQ", "TEST"], kind="stable").drop(columns="_c")
    merged = merged.drop_duplicates(ignore_index=True)
    n_dup = merged.duplicated(["CHR", "UNIT.MASK.MAX_FREQ", "TEST"]).sum()
    if n_dup:
        log(f"WARNING: {n_dup:,} rows share CHR + aggregate + test but differ in their results")
    os.makedirs(cfg.results_local, exist_ok=True)
    for f in (cfg.stamp, cfg.plain):
        if os.path.exists(f):
            os.remove(f)
    merged.to_csv(cfg.plain + ".part", sep="\t", index=False, na_rep="NA")
    os.replace(cfg.plain + ".part", cfg.plain)
    gz = os.path.join(cfg.results_local, os.path.basename(cfg.merged))
    with open(gz, "wb") as out:      # compressed for the bucket while the index is built
        zipper = subprocess.Popen(["pigz" if shutil.which("pigz") else "gzip", "-3", "-c", cfg.plain], stdout=out)
        build_index(cfg, None)
        if zipper.wait() != 0:
            raise RuntimeError("could not compress the merged results")
    stamp = json.load(open(cfg.stamp))
    stamp["source_size"] = os.path.getsize(gz)
    save_text(cfg.stamp, json.dumps(stamp))
    publish([gz], os.path.dirname(cfg.merged))
    if os.path.dirname(os.path.abspath(gz)) != os.path.abspath(os.path.dirname(cfg.merged)):
        os.remove(gz)
    shutil.rmtree(tmp, ignore_errors=True)
    log(f"merge: {len(merged):,} rows -> {cfg.merged}")


# =====================================================================================
# step: index (plain local copy + grep-able IDs + numpy LOG10P) and select
# =====================================================================================
KEY_SHARDS = 4     # the ID file is split in 4 so it can be grepped on 4 cores at once


def index_ready(cfg):
    return all(os.path.exists(f) for f in (cfg.plain, cfg.offsets, cfg.log10p, cfg.chroms, cfg.stamp)) and \
        all(os.path.exists(f"{cfg.keys}.{k}") for k in range(KEY_SHARDS))


def index_stale(cfg):
    """True if the merged file is visible without network cost and its size changed."""
    src = cfg.merged if not is_gs(cfg.merged) else read_path(cfg.merged)
    if not src or not os.path.exists(src):
        return False
    try:
        return json.load(open(cfg.stamp))["source_size"] != os.path.getsize(src)
    except (OSError, ValueError, KeyError):
        return True


def line_offsets(path):
    """Byte offset of the start of every line, plus the file size at the end."""
    out, pos, last = [np.zeros(1, np.int64)], 0, b"\n"
    with open(path, "rb") as f:
        while True:
            b = f.read(1 << 26)
            if not b:
                break
            a = np.frombuffer(b, np.uint8)
            out.append(np.flatnonzero(a == 10).astype(np.int64) + pos + 1)
            pos += len(b)
            last = b[-1:]
    if last != b"\n":
        out.append(np.array([pos], np.int64))
    return np.concatenate(out)


def gunzip_with_offsets(src, dst):
    """Decompress a .gz file to dst (python zlib: faster than the gzip CLI) and return the byte offset
    of every line start plus the file size, worked out in the same pass."""
    offs, pos, last = [np.zeros(1, np.int64)], 0, b"\n"
    dec, pending = zlib.decompressobj(31), False
    with open(src, "rb") as fin, open(dst, "wb") as fout:
        def emit(data):
            nonlocal pos, last
            if data:
                fout.write(data)
                offs.append(np.flatnonzero(np.frombuffer(data, np.uint8) == 10).astype(np.int64) + pos + 1)
                pos += len(data)
                last = data[-1:]
        while True:
            b = fin.read(1 << 24)
            if not b:
                break
            while b:                                   # a file can hold several gzip members
                emit(dec.decompress(b))
                pending = not dec.eof
                if dec.eof:
                    b = dec.unused_data
                    dec = zlib.decompressobj(31)
                else:
                    b = b""
        emit(dec.flush())
    if pending:
        raise RuntimeError(f"{src} is truncated")
    if last != b"\n":
        offs.append(np.array([pos], np.int64))
    return np.concatenate(offs)


def read_index_columns(path):
    """LOG10P (float) and chromosome code (int8) of every row; pyarrow (multithreaded) if available."""
    if not os.environ.get("CARRIER_NO_PYARROW"):
        try:
            import pyarrow as pa
            import pyarrow.csv as pacsv
            t = pacsv.read_csv(path, read_options=pacsv.ReadOptions(block_size=1 << 26),
                               parse_options=pacsv.ParseOptions(delimiter="\t"),
                               convert_options=pacsv.ConvertOptions(
                                   include_columns=["CHR", "LOG10P"], strings_can_be_null=False,
                                   column_types={"CHR": pa.string(), "LOG10P": pa.float64()},
                                   null_values=["NA", "", "NaN", "nan", "N/A", "NULL", "null"]))
            lp = t.column("LOG10P").to_numpy().astype(np.float64)
            ch = t.column("CHR").combine_chunks().dictionary_encode()
            codes = np.array([chrom_order(c) for c in ch.dictionary.to_pylist()], dtype=np.int8)
            ch = codes[ch.indices.to_numpy()]
            head = pd.read_csv(path, sep="\t", usecols=["CHR", "LOG10P"], nrows=20000, dtype={"CHR": str})
            if not (np.array_equal(head["LOG10P"].to_numpy(np.float64), lp[:len(head)], equal_nan=True) and
                    np.array_equal(np.array([chrom_order(c) for c in head["CHR"]], np.int8), ch[:len(head)])):
                raise ValueError("pyarrow and pandas disagree")
            return lp, ch
        except ImportError:
            pass
        except Exception as e:                          # anything unusual: the pandas parser decides
            log(f"   (pyarrow could not read {os.path.basename(path)}: {str(e).splitlines()[0][:100]}; using pandas)")
    d = pd.read_csv(path, sep="\t", usecols=["CHR", "LOG10P"], dtype={"CHR": "category", "LOG10P": np.float64})
    codes = np.array([chrom_order(c) for c in d["CHR"].cat.categories], dtype=np.int8)
    return d["LOG10P"].values, codes[d["CHR"].cat.codes.values]


def build_index(cfg, source_size, offs=None):
    """From the plain merged file: line offsets (to read any row directly), LOG10P / chromosome as numpy
    arrays, and .keys.0-3 (a tab + the aggregate ID per line, in 4 shards that are grepped in parallel)."""
    t = time.time()
    if offs is None:
        offs = line_offsets(cfg.plain)
    nlines = len(offs) - 1
    cuts = [int(np.searchsorted(offs, offs[-1] * k // KEY_SHARDS)) for k in range(KEY_SHARDS)] + [nlines]
    cuts = [min(c, nlines) for c in cuts]

    def shard(k):
        a, b = int(offs[cuts[k]]), int(offs[cuts[k + 1]])
        with open(f"{cfg.keys}.{k}.part", "wb") as out:
            dd = subprocess.Popen(["dd", f"if={cfg.plain}", "iflag=skip_bytes,count_bytes", f"skip={a}",
                                   f"count={b - a}", "bs=4M", "status=none"], stdout=subprocess.PIPE)
            cut = subprocess.Popen(["cut", "-f3"], stdin=dd.stdout, stdout=subprocess.PIPE)
            dd.stdout.close()
            paste = subprocess.Popen(["paste", "/dev/null", "-"], stdin=cut.stdout, stdout=out)
            cut.stdout.close()
            if paste.wait() != 0 or cut.wait() != 0 or dd.wait() != 0:
                raise RuntimeError(f"could not index {cfg.plain}")
        os.replace(f"{cfg.keys}.{k}.part", f"{cfg.keys}.{k}")

    with ThreadPoolExecutor(max_workers=KEY_SHARDS + 1) as ex:
        fc = ex.submit(read_index_columns, cfg.plain)
        for f in [ex.submit(shard, k) for k in range(KEY_SHARDS)]:
            f.result()
        lp, ch = fc.result()
    if nlines - 1 != len(lp):
        raise RuntimeError(f"{cfg.plain}: {nlines - 1} rows but {len(lp)} LOG10P values")
    save_npy(cfg.offsets, offs)
    save_npy(cfg.log10p, lp)
    save_npy(cfg.chroms, ch)
    save_text(cfg.stamp, json.dumps({"source_size": source_size, "rows": int(len(lp)), "shard_start": cuts[:-1]}))
    for old in (cfg.old_db, cfg.keys):
        if os.path.exists(old):
            os.remove(old)
    log(f"index: {len(lp):,} rows indexed in {time.time() - t:.1f}s")


def step_index(cfg):
    """Plain local copy of the merged results (read straight from the bucket mount if there is one)."""
    t = time.time()
    os.makedirs(cfg.results_local, exist_ok=True)
    old_local = os.path.join(cfg.results_local, os.path.basename(cfg.merged))   # earlier version's copy
    if not is_gs(cfg.merged):
        src = cfg.merged
    else:
        src = read_path(cfg.merged) or (old_local if os.path.exists(old_local) and "index" not in cfg.redo
                                        else fetch(cfg.merged, cfg.results_local, refresh=True))
    if not os.path.exists(src):
        sys.exit(f"index: {cfg.merged} not found (run the merge step)")
    size = os.path.getsize(src)
    if os.path.exists(cfg.stamp):
        os.remove(cfg.stamp)
    offs = None
    if src.endswith(".gz"):
        offs = gunzip_with_offsets(src, cfg.plain + ".part")
    else:
        with open(cfg.plain + ".part", "wb") as out:
            r = subprocess.run(decompress_cmd(src), stdout=out, stderr=subprocess.PIPE)
        if r.returncode != 0:
            raise RuntimeError(f"could not decompress {src}: {r.stderr.decode()[-500:]}")
    os.replace(cfg.plain + ".part", cfg.plain)
    log(f"index: {size / 1e6:.0f} MB merged file decompressed in {time.time() - t:.1f}s")
    build_index(cfg, size, offs)
    if src == old_local and is_gs(cfg.merged):
        os.remove(old_local)


_COORD_TAIL = re.compile(r"_([0-9]+)_([0-9]+)$")


def _derive(agg_id):
    # coding_CHEK2_ENST00000402731.lof.0.001 -> unit, gene, last field (ENST/ENSR/window end), coordinates
    reg = str(agg_id).split(".", 1)[0]
    unit, _, rest = reg.partition("_")
    head, sep, _ = rest.rpartition("_")
    m = _COORD_TAIL.search(reg)
    return (unit, (head if sep else rest).upper(), reg.rpartition("_")[2].upper(),
            int(m.group(1)) if m else None, int(m.group(2)) if m else None)


def derived_columns(d):
    rows = [_derive(x) for x in d["UNIT.MASK.MAX_FREQ"].values]
    out = pd.DataFrame(rows, columns=["unit", "gene", "last", "cstart", "cend"], index=d.index)
    out.insert(1, "mask_name", [str(x).split(".", 1)[0] for x in d["MASK.MAX_FREQ"].values])
    out["cstart"] = pd.to_numeric(out["cstart"], errors="coerce")
    out["cend"] = pd.to_numeric(out["cend"], errors="coerce")
    return out


COORD_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y):([0-9,]+)-([0-9,]+)$", re.IGNORECASE)


def classify_targets(entries):
    t = {"ids": set(), "genes": set(), "ensr": set(), "regions": []}
    for e in entries:
        e = str(e).strip()
        if not e or e.startswith("#"):
            continue
        m = COORD_RE.match(e)
        if m:
            ch = m.group(1).upper()
            t["regions"].append(("23" if ch == "X" else ch, int(m.group(2).replace(",", "")),
                                 int(m.group(3).replace(",", "")), e))
        elif "_" in e and "." in e:
            t["ids"].add(e)
        elif e.upper().startswith("ENSR"):
            t["ensr"].add(e.upper())
        else:
            t["genes"].add(e.upper())
    return t


def _grep_rows(cfg, args):
    """Data-row numbers (0-based) of the .keys lines grep -n matches; the shards are searched at once."""
    shard_start = json.load(open(cfg.stamp))["shard_start"]
    env = {**os.environ, "LC_ALL": "C"}
    procs = [subprocess.Popen(["grep", "-n", *args, f"{cfg.keys}.{k}"], stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, env=env) for k in range(KEY_SHARDS)]
    out = []
    for k, p in enumerate(procs):
        so, se = p.communicate()
        if p.returncode not in (0, 1):
            raise RuntimeError(f"grep failed: {se.decode()[-500:]}")
        lines = np.fromiter((int(l[:l.index(b":")]) for l in so.split(b"\n") if l), dtype=np.int64)
        out.append(lines + shard_start[k] - 2)                # line numbers -> data rows (line 1 = header)
    rows = np.concatenate(out)
    return rows[rows >= 0]


def candidate_rows(cfg, t, targeted):
    """Rows that may match (a superset); the exact filters run on these few rows afterwards."""
    parts = []
    pats = set()
    for i in t["ids"]:
        pats.add("\t" + i)
    for g in t["genes"]:
        pats |= {f"_{g}_", f"_{g}.", f"\t{g}."}
    for e in t["ensr"]:
        pats |= {f"_{e}.", f"\t{e}."}
    if pats:
        os.makedirs(cfg.work, exist_ok=True)
        pf = os.path.join(cfg.work, "select_patterns.txt")
        with open(pf, "w") as f:
            f.write("\n".join(sorted(pats)) + "\n")
        parts.append(_grep_rows(cfg, ["-F", "-i", "-f", pf]))
    if t["regions"]:                                      # every row of those chromosomes (contiguous)
        chrom = np.load(cfg.chroms, mmap_mode="r")
        codes = sorted({chrom_order(ch) for ch, *_ in t["regions"]})
        parts.append(np.flatnonzero(np.isin(chrom, codes)))
    if not targeted:
        lp = np.load(cfg.log10p, mmap_mode="r")
        x = float(cfg.sel["log10p_min"])
        parts.append(np.flatnonzero(lp >= x - 1e-9 * max(1.0, abs(x))))
    if not parts:
        return np.zeros(0, np.int64)
    x = np.sort(np.concatenate(parts), kind="stable")
    return x[np.r_[True, x[1:] != x[:-1]]] if len(x) else x


def read_rows(cfg, rows):
    """These rows of the merged results (runs of consecutive rows are read in one go)."""
    offs = np.load(cfg.offsets, mmap_mode="r")
    rows = np.asarray(rows, dtype=np.int64)
    breaks = np.flatnonzero(np.diff(rows) != 1) + 1
    starts, ends = np.r_[0, breaks], np.r_[breaks, len(rows)]
    with open(cfg.plain, "rb") as f:
        chunks = [f.readline()]
        for s0, e0 in zip(starts, ends):
            if e0 <= s0:
                continue
            a, b = int(offs[rows[s0] + 1]), int(offs[rows[e0 - 1] + 2])
            f.seek(a)
            chunks.append(f.read(b - a))
    return pd.read_csv(io.BytesIO(b"".join(chunks)), sep="\t", dtype={"CHR": str, "DF": str})


def separate(s, sep, n, merge):
    parts = s.str.split(sep, n=(n - 1) if merge else -1, regex=True, expand=True).reindex(columns=range(n))
    return [parts[k].where(parts[k].notna()) for k in range(n)]


def rsplit_last(s, sep):
    parts = s.str.rsplit(sep, n=1, expand=True).reindex(columns=range(2))
    first, last = parts[0], parts[1]
    return first.where(first.notna()), last.where(last.notna())


def tidy_aggregates_table(d):
    # Python version of the R tidy_aggregates_table() (gene/ENST split on the LAST "_")
    d = d.copy()
    unit, rest = separate(d["UNIT.MASK.MAX_FREQ"], "_", 2, merge=True)
    reg, _ = separate(rest, r"\.", 2, merge=True)
    d["unit"], d["reg"] = unit, reg
    is_gene = d["unit"].isin(["coding", "proximal"])
    gene, tmp_region = rsplit_last(d["reg"].where(is_gene), "_")
    nongene = d["reg"].where(~is_gene)
    chrregion = nongene.where(nongene.str.contains("chr|window", na=False))
    _, start, end = separate(chrregion, "_", 3, merge=False)
    region = "chr" + d["CHR"].astype(str) + ":" + start.fillna("NA") + "-" + end.fillna("NA")
    region = region.where(~region.str.contains("NA"))
    region = region.where(~nongene.str.contains("ENSR", na=False), nongene)
    region = region.where(tmp_region.isna(), tmp_region)
    is_rna = d["MASK.MAX_FREQ"].str.contains("rna|RNA", na=False)
    _, rna_region = rsplit_last(d["reg"].where(is_rna), "_")
    region = region.where(~is_rna, rna_region)
    d["gene"], d["region"] = gene, region
    d["singleton"] = d["MASK.MAX_FREQ"].str.contains("singleton", na=False)
    return d


def prep_for_carriers(d):
    # Python version of the R prep_for_carriers()
    d = d.copy()
    annotation, frequency = separate(d["MASK.MAX_FREQ"], r"\.", 2, merge=True)
    tmp99 = frequency.str.contains("99", na=False)
    frequency = frequency.str.replace(r"99.", "", regex=True)
    annotation = annotation.where(~tmp99, annotation + ".99")
    d["annotation"], d["frequency"] = annotation, frequency
    d["CHR"] = d["CHR"].astype(str).str.replace("23", "X", regex=False)
    d["region"] = d["region"].str.replace("chr23", "chrX", regex=False)
    cols = ["CHR", "POS", "unit", "gene", "region", "annotation", "frequency", "MAC", "LOG10P", "BETA", "singleton"]
    return d[cols]


def target_hits(d, t, der=None):
    der = derived_columns(d) if der is None else der
    hit = pd.Series(False, index=d.index)
    if t["ids"]:
        hit |= d["UNIT.MASK.MAX_FREQ"].isin(t["ids"])
    if t["genes"]:
        hit |= der["gene"].isin(t["genes"]) | der["last"].isin(t["genes"])
    if t["ensr"]:
        hit |= der["last"].isin(t["ensr"])
    chrom = d["CHR"].astype(str).replace({"X": "23"})
    for ch, s, e, _ in t["regions"]:
        overlap = ((der["cstart"] <= e) & (der["cend"] >= s)).fillna(False)
        hit |= (chrom == ch) & (overlap | d["POS"].between(s, e))
    return hit


def step_select(cfg):
    t0 = time.time()
    sel = cfg.sel
    targets = list(sel["targets"])
    if sel["targets_file"]:
        tf = fetch(sel["targets_file"], os.path.join(cfg.cache, "tmp"), refresh=True)
        targets += [l.strip() for l in open(tf) if l.strip() and not l.startswith("#")]
    t = classify_targets(targets)
    targeted = any([t["ids"], t["genes"], t["ensr"], t["regions"]])
    if not targeted and sel["log10p_min"] is None:
        sys.exit("selection: give targets / targets_file, or log10p_min")
    log(f"select: {len(t['ids'])} aggregate IDs, {len(t['genes'])} genes/transcripts, {len(t['ensr'])} ENSR, "
        f"{len(t['regions'])} regions; masks {sel['masks'] or 'any'}; units {sel['units'] or 'any'}; "
        f"LOG10P >= {sel['log10p_min']}")
    cand = candidate_rows(cfg, t, targeted)
    if not targeted and not sel["masks"] and not sel["units"] and len(cand) > 20 * sel["max_aggregates"]:
        sys.exit(f"selection: {len(cand):,} result rows have LOG10P >= {sel['log10p_min']} (max_aggregates = "
                 f"{sel['max_aggregates']}); raise log10p_min, add masks/units, or raise selection.max_aggregates")
    rows = read_rows(cfg, cand)
    rows["CHR"] = rows["CHR"].astype(str)
    der = derived_columns(rows)
    keep = target_hits(rows, t, der) if targeted else pd.Series(True, index=rows.index)
    if sel["masks"]:
        keep &= der["mask_name"].isin(sel["masks"]) | rows["MASK.MAX_FREQ"].isin(sel["masks"])
    if sel["units"]:
        keep &= der["unit"].isin(sel["units"])
    if sel["log10p_min"] is not None:
        keep &= rows["LOG10P"] >= float(sel["log10p_min"])
    if t["ids"]:
        keep |= rows["UNIT.MASK.MAX_FREQ"].isin(t["ids"])            # full aggregate IDs are always kept
    rows, der = rows[keep], der[keep]
    rows = rows.sort_values("LOG10P", ascending=False, kind="stable").drop_duplicates(["CHR", "UNIT.MASK.MAX_FREQ"])
    der = der.loc[rows.index]
    if targeted:
        for kind, key in (("aggregate ID", "ids"), ("gene/transcript", "genes"), ("ENSR id", "ensr")):
            for item in sorted(t[key]):
                one = {"ids": set(), "genes": set(), "ensr": set(), "regions": []}
                one[key] = {item}
                if not target_hits(rows, one, der).any():
                    log(f"   no aggregate selected for {kind} {item}")
        for r in t["regions"]:
            if not target_hits(rows, {"ids": set(), "genes": set(), "ensr": set(), "regions": [r]}, der).any():
                log(f"   no aggregate selected for region {r[3]}")
    if rows.empty:
        sys.exit("selection: nothing selected")
    if len(rows) > sel["max_aggregates"]:
        sys.exit(f"selection: {len(rows):,} aggregates; narrow it or raise selection.max_aggregates")
    agg = prep_for_carriers(tidy_aggregates_table(rows))
    agg["ID"] = rows["UNIT.MASK.MAX_FREQ"].values
    bad = agg["region"].isna()
    if bad.any():
        log(f"WARNING: dropping {bad.sum()} aggregates whose region couldn't be parsed: {list(agg.loc[bad, 'ID'][:5])}")
        agg = agg[~bad]
    agg = agg.assign(_c=agg["CHR"].map(chrom_order)).sort_values(["_c", "POS"], kind="stable") \
        .drop(columns="_c").reset_index(drop=True)
    agg["singleton"] = agg["singleton"].map({True: "TRUE", False: "FALSE"})
    agg["MAC"] = pd.to_numeric(agg["MAC"], errors="coerce").round().astype("Int64")
    local = os.path.join(cfg.work, os.path.basename(cfg.agg_file))
    os.makedirs(cfg.work, exist_ok=True)
    agg.to_csv(local, index=False, na_rep="")
    _PENDING.append(_UPLOADS.submit(publish, [local], os.path.dirname(cfg.agg_file)))   # upload in the background
    timed("select", t0)
    log(f"select: {len(agg)} aggregates -> {cfg.agg_file}  ({agg['unit'].value_counts().to_dict()})")
    return agg


# =====================================================================================
# step: phenotype and sample files
# =====================================================================================
def fid_iid(df, file):
    df = df.rename(columns={col: col.lstrip("#") for col in df.columns})
    lower = {col.lower(): col for col in df.columns}
    id_col = next((lower[k] for k in ("iid", "person_id", "research_id", "s") if k in lower), None)
    if id_col is None:
        raise ValueError(f"{file}: no IID / person_id column")
    df = df.rename(columns={id_col: "IID"})
    if "fid" in lower:
        df = df.rename(columns={lower["fid"]: "FID"})
    else:
        df.insert(0, "FID", df["IID"])
    df["FID"], df["IID"] = df["FID"].str.strip(), df["IID"].str.strip()
    return df


def read_any_table(path, folder, **kwargs):
    local = fetch(path, folder)
    first = read_text_table(local, sep="\t", nrows=1, dtype=str)
    sep = "\t" if first.shape[1] > 1 else r"\s+"
    return read_text_table(local, sep=sep, dtype=str, **kwargs)


def step_samples(cfg):
    inp = cfg.inp
    tmp = os.path.join(cfg.cache, "tmp")
    phenos = fid_iid(read_any_table(inp["pheno_file"], tmp), inp["pheno_file"])
    col = inp["pheno_col"] or [x for x in phenos.columns if x not in ("FID", "IID")][0]
    ph = phenos[["FID", "IID", col]].rename(columns={col: "phenotype"})
    ph["phenotype"] = pd.to_numeric(ph["phenotype"], errors="coerce")
    binary = set(ph["phenotype"].dropna().unique()).issubset({0, 1})
    if binary:
        ph["phenotype"] = (ph["phenotype"] + 1).astype("Int64")
    cov = fid_iid(read_any_table(inp["cov_file"], tmp,
                                 usecols=lambda x: x.lstrip("#").lower() in ("fid", "iid", "person_id")), inp["cov_file"])
    samples = cov[["FID", "IID"]].drop_duplicates()
    samples = samples[samples["IID"].isin(set(ph.loc[ph["phenotype"].notna(), "IID"]))]
    if cfg.c["ancestry"] and inp["ancestry_file"]:
        anc = read_any_table(inp["ancestry_file"], tmp,
                             usecols=lambda x: x in ("research_id", "person_id", "s", "ancestry_pred"))
        idc = next(x for x in ("research_id", "person_id", "s") if x in anc.columns)
        groups = sorted(anc["ancestry_pred"].dropna().unique())
        if cfg.c["ancestry"] not in groups:
            sys.exit(f"ancestry '{cfg.c['ancestry']}' not in {inp['ancestry_file']}; available: {groups}")
        samples = samples[samples["IID"].isin(set(anc.loc[anc["ancestry_pred"] == cfg.c["ancestry"], idc].str.strip()))]
    os.makedirs(cfg.meta_local, exist_ok=True)
    pf, sf = os.path.join(cfg.meta_local, "phenotype_file.txt"), os.path.join(cfg.meta_local, "sample_file_no_NA.txt")
    ph.to_csv(pf, sep="\t", index=False, na_rep="NA")
    samples.to_csv(sf, sep="\t", index=False)
    publish([pf, sf], cfg.meta_out)
    log(f"samples: phenotype '{col}' ({'binary 0/1 -> 1/2' if binary else 'quantitative'}); "
        f"{len(samples):,} people in the sample file -> {cfg.meta_out}")


def load_samples(cfg):
    pf = fetch(f"{cfg.meta_out}/phenotype_file.txt", cfg.meta_local)
    sf = fetch(f"{cfg.meta_out}/sample_file_no_NA.txt", cfg.meta_local)
    stamp = [(os.path.getsize(f), os.path.getmtime(f), 3) for f in (pf, sf)]
    os.makedirs(cfg.meta_local, exist_ok=True)
    pk = os.path.join(cfg.meta_local, "samples.pickle")
    try:
        with open(pk, "rb") as f:
            saved = pickle.load(f)
        if saved["stamp"] == stamp:
            return saved["pheno"], saved["keep"], saved["hash"]
    except (OSError, EOFError, KeyError, pickle.UnpicklingError, ValueError):
        pass
    ph = pd.read_csv(pf, sep=r"\s+", dtype=str)
    ph["IID"] = ph["IID"].str.strip()
    ph = ph.drop_duplicates("IID")
    pheno = {"iid": ph["IID"].to_numpy(dtype=object), "str": ph["phenotype"].to_numpy(dtype=object),
             "num": pd.to_numeric(ph["phenotype"], errors="coerce").values.astype(np.float64)}
    keep = set(pd.read_csv(sf, sep=r"\s+", dtype=str)["IID"].str.strip())
    khash = hashlib.sha1("\n".join(sorted(keep)).encode()).hexdigest()[:12]
    try:
        with open(pk + ".part", "wb") as f:
            pickle.dump({"stamp": stamp, "pheno": pheno, "keep": keep, "hash": khash}, f, protocol=4)
        os.replace(pk + ".part", pk)
    except OSError:
        pass
    return pheno, keep, khash


# =====================================================================================
# step: extract -- annotation lookup
# =====================================================================================
def mask_table(cfg, CHR, cache_dir):
    local = fetch(cfg.chr_path(cfg.inp["mask_file"], CHR), cache_dir)
    m = read_text_table(local, sep=r"\s+", header=None, usecols=[0, 1], names=["mask", "annotations"], dtype=str)
    return {r.mask: [a.strip() for a in r.annotations.split(",")] for r in m.itertuples()}


def region_key(row, range_units):
    region = str(row.region)
    if row.unit in range_units:
        return None
    if region.startswith("EN"):
        return region
    _, start, end = region.replace("chr", "").replace(":", "-").split("-")
    return f"{start}_{end}"


def annotation_rows(cfg, CHR, keys, work):
    """Annotation rows for these keys. On a chromosome's first use the annotation file is read
    through the bucket mount and copied to the cache in the same pass."""
    acache = os.path.join(cfg.cache, "annotations")
    path = cfg.chr_path(cfg.inp["anno_file"], CHR)
    if not is_gs(path):
        return grep_annotations(path, keys, work, acache)
    local = os.path.join(acache, os.path.basename(path))
    if os.path.exists(local):
        return grep_annotations(local, keys, work, acache)
    mounted = read_path(path, auto=cfg.perf["geno_access"] != "copy")
    if mounted and decompress_stdin_cmd(path) is not None:
        os.makedirs(acache, exist_ok=True)
        return grep_annotations(mounted, keys, work, acache, tee_to=local)
    return grep_annotations(fetch(path, acache), keys, work, acache)


def grep_annotations(src, keys, work, cache_dir, tee_to=None):
    """Keys looked up before come from a small per-file SQLite cache; the rest are found with one
    pass of zstd -dc | grep -F and added to the cache."""
    cols = ["variant", "set", "annotation"]
    if not keys:
        return pd.DataFrame(columns=cols)
    st = os.stat(src)
    stamp = f"{st.st_size}:{int(st.st_mtime)}"
    os.makedirs(cache_dir, exist_ok=True)
    con = sqlite3.connect(os.path.join(cache_dir, os.path.basename(src) + ".lookup.sqlite"), timeout=60)
    con.execute("CREATE TABLE IF NOT EXISTS meta (stamp TEXT)")
    con.execute("CREATE TABLE IF NOT EXISTS done (k TEXT PRIMARY KEY)")
    con.execute("CREATE TABLE IF NOT EXISTS rows (k TEXT, variant TEXT, \"set\" TEXT, annotation TEXT)")
    con.execute("CREATE INDEX IF NOT EXISTS idx_k ON rows (k)")
    old = con.execute("SELECT stamp FROM meta").fetchone()
    if old is None or old[0] != stamp:             # annotation file changed: start the cache again
        con.execute("DELETE FROM done"); con.execute("DELETE FROM rows"); con.execute("DELETE FROM meta")
        con.execute("INSERT INTO meta VALUES (?)", (stamp,))
    keys = sorted(keys)
    q = ",".join("?" * len(keys))
    done = {r[0] for r in con.execute(f"SELECT k FROM done WHERE k IN ({q})", keys)}
    todo = sorted(set(keys) - done)
    if todo:
        new = _grep(src, todo, work, tee_to)
        found = _assign_keys(new, set(todo))
        con.executemany("INSERT INTO rows VALUES (?, ?, ?, ?)", found)
        con.executemany("INSERT OR IGNORE INTO done VALUES (?)", [(k,) for k in todo])
        con.commit()
    out = pd.read_sql_query(f"SELECT DISTINCT variant, \"set\", annotation FROM rows WHERE k IN ({q})", con,
                            params=keys)
    con.close()
    return out[cols]


def _assign_keys(rows, keys):
    """(key, variant, set, annotation) for each grep hit. Set names end in their key
    (..._ENST00000402731, ..._5687000_5688000), so match on the last one or two '_' fields;
    the few rows that don't match that way are checked by substring."""
    if rows.empty:
        return []
    last1 = rows["set"].str.rsplit("_", n=1).str[-1]
    last2 = rows["set"].str.rsplit("_", n=2).str[-2:].str.join("_")
    out, matched = [], np.zeros(len(rows), dtype=bool)
    for suffix in (last1, last2):
        m = suffix.isin(keys).values
        out += list(zip(suffix[m], rows["variant"][m], rows["set"][m], rows["annotation"][m]))
        matched |= m
    for v, st, an in rows[~matched].itertuples(index=False):
        out += [(k, v, st, an) for k in keys if k in st]
    return out


def _grep(src, keys, work, tee_to=None):
    kf = os.path.join(work, "annotation_keys.txt")
    with open(kf, "w") as f:
        f.write("\n".join(sorted(keys)) + "\n")
    procs = []
    if tee_to:      # bucket file -> tee (local cache copy) -> decompress -> grep, all in one pass
        part = tee_to + ".part"
        fin = open(src, "rb")
        tee = subprocess.Popen(["tee", part], stdin=fin, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        fin.close()
        unzip = subprocess.Popen(decompress_stdin_cmd(src), stdin=tee.stdout, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        tee.stdout.close()
        procs.append(tee)
    else:
        unzip = subprocess.Popen(decompress_cmd(src), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
    procs.append(unzip)
    grep = subprocess.Popen(["grep", "-F", "-f", kf], stdin=unzip.stdout, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, env={**os.environ, "LC_ALL": "C"})
    unzip.stdout.close()
    out, err = grep.communicate()
    for p in procs:
        p.wait()
        if p.returncode != 0:
            raise RuntimeError(f"could not read {src}: {p.stderr.read().decode()[-500:]}")
    if grep.returncode not in (0, 1):
        raise RuntimeError(f"annotation lookup failed: {err.decode()[-500:]}")
    if tee_to:
        st = os.stat(src)
        if os.path.getsize(part) != st.st_size:
            raise RuntimeError(f"incomplete copy of {src}")
        os.utime(part, (st.st_atime, st.st_mtime))      # same time stamp as the bucket file
        os.replace(part, tee_to)
    rows = [l.split() for l in out.decode().splitlines() if l.strip()]
    return pd.DataFrame([x[:3] for x in rows if len(x) >= 3], columns=["variant", "set", "annotation"])


def annotation_index(ann):
    """Annotation rows grouped by set name, so each aggregate only looks at its own sets."""
    if ann is None:
        return None
    sets = ann["set"].to_numpy(dtype=object)
    codes, names = pd.factorize(sets)
    order = np.argsort(codes, kind="stable")
    bounds = np.searchsorted(codes[order], np.arange(len(names) + 1))
    return {"names": list(names), "rows": {nm: order[bounds[k]:bounds[k + 1]] for k, nm in enumerate(names)},
            "variant": ann["variant"].to_numpy(dtype=object), "annotation": ann["annotation"].to_numpy(dtype=object)}


def step_a(row, ai, masks, range_units):
    if row.unit in range_units:
        c, start, end = str(row.region).replace("chr", "").replace(":", "-").split("-")
        return {"mode": "range", "range": (int(start), int(end)), "variants": [], "sets": ""}
    if row.annotation not in masks:
        raise ValueError(f"mask '{row.annotation}' not in the mask file")
    key = region_key(row, range_units)
    matched = [nm for nm in ai["names"] if isinstance(nm, str) and key in nm]
    if len(matched) > 1 and any(nm.startswith(f"{row.unit}_") for nm in matched):
        matched = [nm for nm in matched if nm.startswith(f"{row.unit}_")]
    idx = np.sort(np.concatenate([ai["rows"][nm] for nm in matched])) if matched else np.zeros(0, np.int64)
    ok = pd.Index(masks[row.annotation]).get_indexer(ai["annotation"][idx]) >= 0
    variants = pd.unique(ai["variant"][idx][ok]).tolist()
    return {"mode": "ids", "variants": variants, "sets": ";".join(sorted(matched))}


# =====================================================================================
# step: extract -- genotype files (only touched when the local store lacks a variant)
# =====================================================================================
def small_file(path, gcache):
    """Local copy of a genotype index file (.psam/.fam/.bim/.sample)."""
    if not is_gs(path):
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        return path
    return fetch(path, gcache)


def slim_pvar(cfg, prefix, gcache):
    """Local .pvar with only #CHROM POS ID REF ALT (plink2 parses it much faster than one with INFO)."""
    slim = os.path.join(gcache, os.path.basename(prefix) + ".slim.pvar")
    if os.path.exists(slim):
        return slim
    os.makedirs(gcache, exist_ok=True)
    legacy = os.path.join(gcache, os.path.basename(prefix) + ".pvar")      # full copy from an earlier version
    src = None
    if os.path.exists(legacy):
        src = legacy
    elif is_gs(prefix):
        for ext in (".pvar", ".pvar.zst"):
            src = read_path(prefix + ext, auto=cfg.perf["geno_access"] != "copy")
            if src:
                break
        if src is None:
            ext = ".pvar" if exists(prefix + ".pvar") else ".pvar.zst"
            src = fetch(prefix + ext, gcache)
    else:
        src = prefix + ".pvar" if os.path.exists(prefix + ".pvar") else prefix + ".pvar.zst"
        if not os.path.exists(src):
            raise FileNotFoundError(prefix + ".pvar")
    t = time.time()
    # 1) header: '##' lines (minus INFO/FILTER/FORMAT definitions) and the #CHROM line
    reader = subprocess.Popen(decompress_cmd(src), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE)
    head, k, header_ok = [], 0, False
    for line in reader.stdout:
        if line.startswith(b"#CHROM"):
            k += 1
            f = line.rstrip(b"\n").split(b"\t")
            header_ok = [x.upper() for x in f[:5]] == [b"#CHROM", b"POS", b"ID", b"REF", b"ALT"]
            head.append(b"\t".join(f[:5]) + b"\n" if header_ok else line)
            break
        if not line.startswith(b"##"):           # no #CHROM line: data starts here
            break
        k += 1
        if not line.startswith((b"##INFO", b"##FILTER", b"##FORMAT")):
            head.append(line)
    reader.kill()
    reader.wait()
    # 2) data lines, first five columns only
    with open(slim + ".part", "wb") as out:
        out.write(b"".join(head))
        out.flush()
        p1 = subprocess.Popen(decompress_cmd(src), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["tail", "-n", f"+{k + 1}"], stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        p3 = subprocess.Popen(["cut", "-f1-5"] if header_ok else ["cat"], stdin=p2.stdout, stdout=out)
        p2.stdout.close()
        p3.wait(), p2.wait(), p1.wait()
    if p1.returncode != 0 or p2.returncode != 0 or p3.returncode != 0:
        raise RuntimeError(f"could not read {src}")
    os.replace(slim + ".part", slim)
    log(f"   cached a slim copy of {os.path.basename(src)} ({os.path.getsize(slim) / 1e6:.0f} MB, "
        f"{time.time() - t:.1f}s; once per chromosome)")
    if is_gs(prefix) and os.path.dirname(os.path.abspath(src)) == os.path.abspath(gcache):
        os.remove(src)                           # full local copy no longer needed
    return slim


def genotype_input(cfg, CHR):
    """plink2 input flags. The big genotype file is read in place (bucket mount) when possible; the
    index files (.pvar slimmed, .psam etc.) are cached on local disk, since plink2 reads them in full."""
    fmt = cfg.inp["geno_format"]
    prefix = cfg.chr_path(cfg.inp["geno_prefix"], CHR)
    gcache = os.path.join(cfg.cache, "genotype_index")
    if fmt == "pgen":
        big, small = prefix + ".pgen", [slim_pvar(cfg, prefix, gcache), small_file(prefix + ".psam", gcache)]
    elif fmt == "bed":
        big, small = prefix + ".bed", [small_file(prefix + ".bim", gcache), small_file(prefix + ".fam", gcache)]
    elif fmt == "bgen":
        big, small = prefix + ".bgen", [small_file(cfg.chr_path(cfg.inp["bgen_sample"], CHR), gcache)]
    else:
        raise ValueError("inputs.geno_format must be pgen, bed or bgen")

    access = cfg.perf["geno_access"]
    if is_gs(big):
        big_local = read_path(big, auto=True) if access in ("auto", "mount") else None
        if big_local is None and access == "mount":
            raise RuntimeError(f"{big} is not reachable through a bucket mount")
        if big_local is not None:
            log(f"   chr{CHR}: reading {os.path.basename(big)} through the bucket mount")
        else:   # copy the big file (slow; kept only for this chromosome)
            size = gs_size(big) or 0
            os.makedirs(os.path.join(cfg.work, "geno_copy"), exist_ok=True)
            free = shutil.disk_usage(os.path.join(cfg.work, "geno_copy")).free
            log(f"   chr{CHR}: no bucket mount; copying {size / 1e9:.1f} GB ({free / 1e9:.0f} GB free)")
            if size > free - 5e9:
                raise RuntimeError(f"not enough disk to copy {big}; mount the bucket instead")
            big_local = fetch(big, os.path.join(cfg.work, "geno_copy"))
    else:
        if not os.path.exists(big):
            raise FileNotFoundError(big)
        big_local = big
    if fmt == "pgen":
        return ["--pgen", big_local, "--pvar", small[0], "--psam", small[1]], small[1]
    if fmt == "bed":
        return ["--bed", big_local, "--bim", small[0], "--fam", small[1]], small[1]
    return ["--bgen", big_local, "ref-first", "--sample", small[0]], small[0]


def keep_file_for(sample_file, keep_iids, out):
    """plink2 --keep file using the genotype file's own sample IDs (FID/IID if present)."""
    if sample_file.endswith(".psam"):
        ps = pd.read_csv(sample_file, sep=r"\s+", dtype=str)
        ps.columns = [x.lstrip("#") for x in ps.columns]
        cols = ["FID", "IID"] if "FID" in ps.columns else ["IID"]
    elif sample_file.endswith(".fam"):
        ps = pd.read_csv(sample_file, sep=r"\s+", header=None, usecols=[0, 1], names=["FID", "IID"], dtype=str)
        cols = ["FID", "IID"]
    else:   # .sample: ID_1 ID_2, second line is the type row
        ps = pd.read_csv(sample_file, sep=r"\s+", dtype=str, skiprows=[1]).rename(columns={"ID_1": "FID", "ID_2": "IID"})
        cols = ["FID", "IID"]
    k = ps.loc[ps["IID"].isin(keep_iids), cols]
    with open(out + ".part", "w") as f:
        f.write("#" + "\t".join(cols) + "\n")
        k.to_csv(f, sep="\t", header=False, index=False)
    os.replace(out + ".part", out)
    return len(k)


def plink2_extract(cfg, CHR, positions, windows, ids, geno_flags, keep, work):
    """plink2 -> small .bed of the requested variants (analysed people only). Variants are taken by
    position/window ('pos') and, for IDs that aren't chr:pos:ref:alt, by ID ('id'); both run at once."""
    plink_chr = "X" if str(CHR) == "X" else str(CHR)
    extra = []
    if cfg.perf["plink2_threads"]:
        extra += ["--threads", str(cfg.perf["plink2_threads"])]
    if cfg.perf["plink2_memory_mb"]:
        extra += ["--memory", str(cfg.perf["plink2_memory_mb"])]
    jobs = []
    if len(positions) or windows:
        rf = os.path.join(work, "ranges.txt")
        with open(rf, "w") as f:
            f.write("".join(f"{plink_chr}\t{p}\t{p}\n" for p in positions))
            f.write("".join(f"{plink_chr}\t{s}\t{e}\n" for s, e in windows))
        jobs.append(("pos", ["--extract", "bed1", rf]))
    if ids:
        idf = os.path.join(work, "ids.txt")
        with open(idf, "w") as f:
            f.write("\n".join(ids) + "\n")
        jobs.append(("id", ["--extract", idf]))

    def run(job):
        name, ext = job
        out = os.path.join(work, f"piece_{name}")
        cmd = ["plink2", *geno_flags, "--keep", keep, *ext, "--max-alleles", "2", "--make-bed", "--out", out, *extra]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            if "no variants remaining" in (r.stdout + r.stderr).lower():
                return name, None
            raise RuntimeError(f"plink2 failed (chr{CHR}):\n{(r.stdout + r.stderr)[-1500:]}")
        return name, out

    with ThreadPoolExecutor(max_workers=max(1, len(jobs))) as ex:
        return dict(ex.map(run, jobs))


# ---- small local .bed reader ----------------------------------------------------------
# .bed 2-bit codes for A1 = ALT (plink2 --make-bed): 00 hom A1 -> 2 ALT, 01 missing, 10 het -> 1, 11 hom A2 -> 0
_CODE = np.array([2, -1, 1, 0], dtype=np.int8)
LUT = np.stack([_CODE[(np.arange(256) >> s) & 3] for s in (0, 2, 4, 6)], axis=1)
# per-byte tables (4 people per .bed byte): ALT alleles and called genotypes
_ALT = np.array([2, 0, 1, 0], dtype=np.int16)
_CALLED = np.array([1, 0, 1, 1], dtype=np.int16)
T_ALT = sum(_ALT[(np.arange(256) >> s) & 3] for s in (0, 2, 4, 6)).astype(np.int16)
T_CALLED = sum(_CALLED[(np.arange(256) >> s) & 3] for s in (0, 2, 4, 6)).astype(np.int16)
_SHIFTS = np.array([0, 2, 4, 6], dtype=np.uint8)
_FORCE_DENSE = False   # tests only
MAF_TOL = 1 + 5.684341886080801e-14   # plink 1.9's tolerance on --max-maf


def prepare_rows(Bb, n, csr_maf):
    """For packed .bed rows Bb (variants x bytes, everyone in the file analysed): ALT allele count and
    allele number per variant, and for variants rare enough for any frequency filter up to csr_maf
    (or MAC <= 1) the bytes that aren't 4 x the major homozygote, i.e. carriers and missing calls."""
    nblk, nb = Bb.shape
    pad = nb * 4 - n                          # padding codes in the last byte are 00 (2 ALT, called)
    tol = csr_maf * MAF_TOL
    alt = np.zeros(nblk, np.int64)
    called = np.zeros(nblk, np.int64)
    not_ff = Bb != 0xFF                                  # 0xFF = 4 x REF homozygote
    n_ff = nb - not_ff.sum(axis=1)
    sparse = n_ff * 10 >= nb * 9                         # >= 90% of bytes are 4 x REF homozygote (rare ALT)
    rs = np.flatnonzero(sparse)
    r, c = np.nonzero(not_ff[rs])
    v = Bb[rs[r], c]
    alt[rs] = np.bincount(r, weights=T_ALT[v], minlength=len(rs)).astype(np.int64)
    called[rs] = 4 * n_ff[rs] + np.bincount(r, weights=T_CALLED[v], minlength=len(rs)).astype(np.int64)
    rd = np.flatnonzero(~sparse)
    if len(rd):
        alt[rd] = T_ALT[Bb[rd]].sum(axis=1, dtype=np.int64)
        called[rd] = T_CALLED[Bb[rd]].sum(axis=1, dtype=np.int64)
    alt -= 2 * pad
    an = 2 * (called - pad)
    mia = alt <= an - alt
    mac = np.where(mia, alt, an - alt)
    with np.errstate(divide="ignore", invalid="ignore"):
        maf = np.where(an > 0, mac / an, 0.5)
    rare = (maf <= tol) | (mac <= 1)
    rows_out, cols_out, vals_out = [np.zeros(0, np.int64)], [np.zeros(0, np.int64)], [np.zeros(0, np.uint8)]
    keep_s = rare[rs][r] & mia[rs][r]                    # sparse rows: major is REF (0xFF)
    rows_out.append(rs[r][keep_s]); cols_out.append(c[keep_s]); vals_out.append(v[keep_s])
    for ro in np.flatnonzero(rare & ~(sparse & mia)):    # rare rows not covered above (e.g. ALT major)
        cc = np.flatnonzero(Bb[ro] != (0xFF if mia[ro] else 0x00))
        rows_out.append(np.full(len(cc), ro)); cols_out.append(cc); vals_out.append(Bb[ro, cc])
    ro_all = np.concatenate(rows_out).astype(np.int64)
    order = np.argsort(ro_all, kind="stable")
    counts = np.bincount(ro_all, minlength=nblk)
    return alt, an, counts, np.concatenate(cols_out)[order].astype(np.int32), \
        np.concatenate(vals_out)[order].astype(np.uint8)


def prepare_bed(bed, rows, n, csr_maf, block=None):
    block = block or max(64, min(2000, 20_000_000 // max(bed.shape[1], 1)))     # ~20 MB of .bed per block
    parts = [prepare_rows(np.asarray(bed[rows[s:s + block]]), n, csr_maf) for s in range(0, len(rows), block)]
    if not parts:
        return {"alt": np.zeros(0, np.int64), "an": np.zeros(0, np.int64), "indptr": np.zeros(1, np.int64),
                "nz_byte": np.zeros(0, np.int32), "nz_val": np.zeros(0, np.uint8)}
    counts = np.concatenate([p[2] for p in parts])
    return {"alt": np.concatenate([p[0] for p in parts]), "an": np.concatenate([p[1] for p in parts]),
            "indptr": np.concatenate([[0], np.cumsum(counts)]).astype(np.int64),
            "nz_byte": np.concatenate([p[3] for p in parts]), "nz_val": np.concatenate([p[4] for p in parts])}


class BedSet:
    ARRAYS = ("alt", "an", "indptr", "nz_byte", "nz_val")

    def __init__(self, prefix, prepared=False):
        self.bim = pd.read_csv(prefix + ".bim", sep="\t", header=None,
                               names=["chr", "id", "cm", "pos", "alt", "ref"],
                               dtype={"chr": str, "id": str, "cm": str, "alt": str, "ref": str})
        self.pos = self.bim["pos"].values.astype(np.int64)
        if os.path.exists(prefix + ".fam.pkl"):
            self.fam = pd.read_pickle(prefix + ".fam.pkl")
        else:
            self.fam = pd.read_csv(prefix + ".fam", sep=r"\s+", header=None,
                                   names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"], dtype=str)
        self.n = len(self.fam)
        self.famcol = {c: self.fam[c].to_numpy(dtype=object) for c in ("FID", "IID", "PAT", "MAT", "SEX")}
        self.bim_id, self.bim_alt, self.bim_ref = (self.bim[c].to_numpy(dtype=object) for c in ("id", "alt", "ref"))
        nb = (self.n + 3) // 4
        self.bed = np.memmap(prefix + ".bed", dtype=np.uint8, mode="r", offset=3, shape=(len(self.bim), nb))
        if prepared:
            for name in self.ARRAYS:
                setattr(self, name, np.load(f"{prefix}.{name}.npy", mmap_mode="r"))
            self.minor_is_alt = self.alt <= self.an - self.alt

    @cached_property
    def by_key(self):
        d = {}
        for i, (p, r, a) in enumerate(zip(self.pos, self.bim["ref"].str.upper(), self.bim["alt"].str.upper())):
            d.setdefault((int(p), r, a), i)
        return d

    @cached_property
    def by_id(self):
        return {v: i for i, v in enumerate(self.bim["id"])}

    def prepare(self, max_maf):
        if hasattr(self, "alt"):
            return
        for k, v in prepare_bed(self.bed, np.arange(len(self.bim)), self.n, max_maf).items():
            setattr(self, k, v)
        self.minor_is_alt = self.alt <= self.an - self.alt


STORE_VERSION = 1
BED_MAGIC = bytes([0x6C, 0x1B, 0x01])


class StoreReset(Exception):
    pass


class Store:
    """Every variant extracted so far for one chromosome and sample set, kept as one local .bed with
    allele counts and carrier bytes precomputed ('pos' = taken by position/window, 'id' = by ID).
    Requests it already covers need no plink2 call at all; new variants are appended."""

    def __init__(self, root, kind):
        self.dir = os.path.join(root, kind)
        self.prefix = os.path.join(self.dir, "store")
        self.kind = kind
        self._fresh()
        mf = os.path.join(self.dir, "meta.json")
        if not os.path.exists(mf):
            return
        try:
            m = json.load(open(mf))
            ok = m.get("version") == STORE_VERSION
            if ok and m["nv"]:
                ok = os.path.getsize(self.prefix + ".bed") == 3 + m["nv"] * ((m["n"] + 3) // 4)
            if ok:
                self.meta = m
                if kind == "pos":
                    self.pos = np.load(os.path.join(self.dir, "covered.npy"))
                else:
                    self.ids = set(open(os.path.join(self.dir, "covered.txt")).read().splitlines())
        except (OSError, ValueError, KeyError):
            ok = False
        if not ok:
            log(f"   genotype store {self.dir} was incomplete; starting it again")
            self.reset()

    def _fresh(self):
        self.meta = {"version": STORE_VERSION, "nv": 0, "n": None, "csr_maf": None, "windows": []}
        self.pos = np.zeros(0, np.int64)
        self.ids = set()

    def reset(self):
        shutil.rmtree(self.dir, ignore_errors=True)
        self._fresh()

    def size(self):
        try:
            return os.path.getsize(self.prefix + ".bed")
        except OSError:
            return 0

    def missing(self, positions, windows, ids):
        if self.kind == "id":
            return np.zeros(0, np.int64), [], sorted(set(ids) - self.ids)
        p = np.unique(np.asarray(positions, dtype=np.int64))
        covered = np.isin(p, self.pos)
        cached = [tuple(w) for w in self.meta["windows"]]
        for s, e in cached:
            covered |= (p >= s) & (p <= e)
        new_w = [w for w in windows if not any(s <= w[0] and w[1] <= e for s, e in cached)]
        return p[~covered], new_w, []

    def _arrays(self):
        return {k: np.load(f"{self.prefix}.{k}.npy") for k in BedSet.ARRAYS}

    def _save_arrays(self, arrays):
        for k, v in arrays.items():
            save_npy(f"{self.prefix}.{k}.npy", v)

    def _save_meta(self):
        save_text(os.path.join(self.dir, "meta.json"), json.dumps(self.meta))

    def ensure_csr(self, csr_maf):
        """Carrier bytes must exist for every variant a frequency filter up to csr_maf can keep."""
        if not self.meta["nv"] or (self.meta["csr_maf"] or 0) >= csr_maf:
            return
        nv, n = self.meta["nv"], self.meta["n"]
        bed = np.memmap(self.prefix + ".bed", np.uint8, "r", offset=3, shape=(nv, (n + 3) // 4))
        self._save_arrays(prepare_bed(bed, np.arange(nv), n, csr_maf))
        self.meta["csr_maf"] = csr_maf
        self._save_meta()

    def add(self, piece, new_pos, new_win, new_ids, csr_maf):
        os.makedirs(self.dir, exist_ok=True)
        if piece is not None:
            self._append(piece, csr_maf)
        if self.kind == "pos":
            self.pos = np.union1d(self.pos, np.asarray(new_pos, dtype=np.int64))
            self.meta["windows"] = sorted({tuple(w) for w in self.meta["windows"]} | {tuple(w) for w in new_win})
            save_npy(os.path.join(self.dir, "covered.npy"), self.pos)
        else:
            self.ids |= set(new_ids)
            save_text(os.path.join(self.dir, "covered.txt"), "\n".join(sorted(self.ids)))
        self._save_meta()                       # written last: marks the store as complete

    def _append(self, piece, csr_maf):
        bim_lines = open(piece + ".bim").read().splitlines()
        fam = open(piece + ".fam").read()
        n = fam.count("\n")
        nb = (n + 3) // 4
        if self.meta["nv"] and (self.meta["n"] != n or open(self.prefix + ".fam").read() != fam):
            self.reset()
            raise StoreReset("the sample list changed")
        self.ensure_csr(csr_maf)
        csr = max(csr_maf, self.meta["csr_maf"] or 0)
        key = lambda l: tuple(l.split("\t")[j] for j in (1, 3, 4, 5))          # id, pos, alt, ref
        if self.meta["nv"]:
            with open(self.prefix + ".bim") as f:
                old = {key(l.rstrip("\n")) for l in f}
            take = np.array([i for i, l in enumerate(bim_lines) if key(l) not in old], dtype=np.int64)
        else:
            take = np.arange(len(bim_lines))
        if len(take) == 0:
            return
        piece_bed = np.memmap(piece + ".bed", np.uint8, "r", offset=3, shape=(len(bim_lines), nb))
        first = self.meta["nv"] == 0
        if first:
            with open(self.prefix + ".bed", "wb") as f:
                f.write(BED_MAGIC)
            save_text(self.prefix + ".fam", fam)
            pd.read_csv(self.prefix + ".fam", sep=r"\s+", header=None, dtype=str,
                        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"]).to_pickle(self.prefix + ".fam.pkl")
            open(self.prefix + ".bim", "w").close()
            self.meta["n"] = n
        with open(self.prefix + ".bed", "ab") as f:
            for s in range(0, len(take), 4096):
                f.write(np.ascontiguousarray(piece_bed[take[s:s + 4096]]).tobytes())
        with open(self.prefix + ".bim", "a") as f:
            f.write("".join(bim_lines[i] + "\n" for i in take))
        new = prepare_bed(piece_bed, take, n, csr)
        if not first:
            old_a = self._arrays()
            new["indptr"] = np.concatenate([old_a["indptr"], old_a["indptr"][-1] + new["indptr"][1:]])
            for k in ("alt", "an", "nz_byte", "nz_val"):
                new[k] = np.concatenate([old_a[k], new[k]])
        self._save_arrays(new)
        self.meta["nv"] += len(take)
        self.meta["csr_maf"] = csr

    def bedset(self):
        if not self.meta["nv"]:
            return None
        u = BedSet(self.prefix, prepared=True)
        if len(u.bim) != self.meta["nv"] or len(u.alt) != self.meta["nv"]:
            self.reset()
            raise StoreReset("store files didn't match")
        if not os.path.exists(self.prefix + ".fam.pkl"):
            u.fam.to_pickle(self.prefix + ".fam.pkl")
        return u


def store_root(cfg, CHR, khash):
    inp = cfg.inp
    tag = hashlib.sha1(json.dumps([inp["geno_format"], inp["geno_prefix"], inp["bgen_sample"],
                                   inp["chr_x"]]).encode()).hexdigest()[:10]
    return os.path.join(cfg.cache, "genotypes", f"{tag}_{khash}", f"chr{CHR}")


def genotypes_for(cfg, CHR, ids, windows, keep, khash, work):
    """{'pos': BedSet, 'id': BedSet} holding (at least) every requested variant, from the local
    store; plink2 runs only for positions/windows/IDs the store hasn't seen."""
    root = store_root(cfg, CHR, khash)
    csr = max(float(cfg.c["max_maf"]), 0.01)
    keys = [variant_key(v) for v in ids]
    by_id = sorted({v for v, k in zip(ids, keys) if k is None})
    positions = sorted({k[0] for k in keys if k})
    windows = sorted({tuple(w) for w in windows})
    limit = float(cfg.perf["genotype_cache_gb"]) * 1e9
    with file_lock(os.path.join(root, ".lock")):
        for attempt in range(2):
            stores = {"pos": Store(root, "pos"), "id": Store(root, "id")}
            new = {k: s.missing(positions, windows, by_id) for k, s in stores.items()}
            todo = len(new["pos"][0]) + len(new["pos"][1]) + len(new["id"][2])
            how = "all from the local store"
            try:
                if todo:
                    for k, s in stores.items():
                        if s.size() > limit and (len(new[k][0]) or new[k][1] or new[k][2]):
                            log(f"   genotype store for chr{CHR} passed genotype_cache_gb; starting it again")
                            s.reset()
                            new[k] = s.missing(positions, windows, by_id)
                    t = time.time()
                    geno_flags, sample_file = genotype_input(cfg, CHR)
                    keep_f = os.path.join(root, "keep.txt")
                    if not os.path.exists(keep_f):
                        keep_file_for(sample_file, keep, keep_f)
                    pieces = plink2_extract(cfg, CHR, new["pos"][0], new["pos"][1], new["id"][2],
                                            geno_flags, keep_f, work)
                    if len(new["pos"][0]) or new["pos"][1]:
                        stores["pos"].add(pieces.get("pos"), new["pos"][0], new["pos"][1], [], csr)
                    if new["id"][2]:
                        stores["id"].add(pieces.get("id"), [], [], new["id"][2], csr)
                    how = (f"plink2 for {len(new['pos'][0]):,} new positions, {len(new['pos'][1]):,} windows, "
                           f"{len(new['id'][2]):,} IDs in {time.time() - t:.1f}s; the rest from the local store")
                for s in stores.values():
                    s.ensure_csr(csr)
                unions = {k: s.bedset() for k, s in stores.items() if s.meta["nv"]}
                return unions, how
            except StoreReset as e:
                log(f"   genotype store for chr{CHR}: {e}; rebuilding it")
                if os.path.exists(os.path.join(root, "keep.txt")):
                    os.remove(os.path.join(root, "keep.txt"))
                if attempt:
                    raise


# =====================================================================================
# step: extract -- carriers per aggregate
# =====================================================================================
def carriers_for(row, a, unions, analysed_idx, max_maf):
    """plink 1.9 --max-maf / --max-mac + --recode A + carriers, done on the packed .bed bytes.
    Allele counts come from the store; only bytes that differ from the major homozygote are
    decoded, so the cost depends on the number of carriers, not the number of people."""
    if a["mode"] == "range":
        u = unions.get("pos")
        if u is None:
            return None
        s, e = a["range"]
        rows = np.flatnonzero((u.pos >= s) & (u.pos <= e))
    else:
        keys = [variant_key(v) for v in a["variants"]]
        if any(k is not None for k in keys):        # IDs like 1:10048:C:CT: match on position/REF/ALT
            u = unions.get("pos")
            rows = [u.by_key[k] for k in keys if k is not None and u is not None and k in u.by_key]
        else:                                        # other ID formats: match on the ID itself
            u = unions.get("id")
            rows = [u.by_id[v] for v in a["variants"] if u is not None and v in u.by_id]
        rows = np.array(sorted(set(rows)), dtype=np.int64)
    if u is None or len(rows) == 0:
        return None
    rows = rows[np.lexsort((rows, u.pos[rows]))]     # genomic order (the store appends in batches)
    n_found = len(rows)
    n = u.n
    everyone = len(analysed_idx) == n and not _FORCE_DENSE
    singleton = str(row.singleton).upper() == "TRUE"

    def counted_of(codes, mia):
        g = _CODE[codes]                                          # ALT dosage, -1 = missing
        return np.where(g < 0, -1, np.where(mia, g, 2 - g)).astype(np.int8)

    def frequency_filter(alt, an):
        ref = an - alt
        mia = alt <= ref
        mac = np.where(mia, alt, ref)
        with np.errstate(divide="ignore", invalid="ignore"):
            maf = np.where(an > 0, mac / an, 0.5)
        keep = (mac <= 1) if singleton else (maf <= max_maf * MAF_TOL)
        return mia, keep

    if everyone:
        u.prepare(max_maf)
        minor_is_alt, keep_var = frequency_filter(np.asarray(u.alt[rows]), np.asarray(u.an[rows]))
        rows, minor_is_alt = rows[keep_var], minor_is_alt[keep_var]
        if len(rows) == 0:
            return {"n_found": n_found, "n_after": 0}
        # non-major bytes of the kept variants -> people with a genotype other than major homozygote
        seg = [np.arange(u.indptr[r], u.indptr[r + 1]) for r in rows]
        lens = np.array([len(x) for x in seg])
        pos = np.concatenate(seg) if len(seg) else np.zeros(0, int)
        vrow = np.repeat(np.arange(len(rows)), lens)               # which kept variant each byte belongs to
        c = np.asarray(u.nz_byte[pos]).astype(np.int64)
        codes = (np.asarray(u.nz_val[pos])[:, None] >> _SHIFTS) & 3
        people = c[:, None] * 4 + np.arange(4)
        cnt = counted_of(codes, minor_is_alt[vrow][:, None])
        real = people < n
        hit = (cnt >= 1) & real
        car = np.unique(people[hit])
        counted_car = np.zeros((len(rows), len(car)), dtype=np.int8)   # major homozygote = 0 copies
        sel = real & (cnt != 0)
        pp, vv, cc = people[sel], np.broadcast_to(vrow[:, None], people.shape)[sel], cnt[sel]
        j = np.searchsorted(car, pp)
        inside = (j < len(car)) & (car[np.minimum(j, len(car) - 1)] == pp) if len(car) else np.zeros(len(pp), bool)
        counted_car[vv[inside], j[inside]] = cc[inside]
    else:
        G = LUT[np.asarray(u.bed[rows])].reshape(len(rows), -1)[:, :n][:, analysed_idx]
        called = G >= 0
        alt = np.where(called, G, 0).sum(axis=1).astype(np.int64)
        an = 2 * called.sum(axis=1).astype(np.int64)
        minor_is_alt, keep_var = frequency_filter(alt, an)
        rows, G, minor_is_alt = rows[keep_var], G[keep_var], minor_is_alt[keep_var]
        if len(rows) == 0:
            return {"n_found": n_found, "n_after": 0}
        counted = np.where(G < 0, -1, np.where(minor_is_alt[:, None], G, 2 - G)).astype(np.int8)
        car = np.flatnonzero((counted >= 1).any(axis=0))
        counted_car = counted[:, car]
    vid = u.bim_id[rows]
    a1 = np.where(minor_is_alt, u.bim_alt[rows], u.bim_ref[rows])
    return {"n_found": n_found, "n_after": len(rows), "u": u, "car": car, "counted_car": counted_car,
            "cols": [f"{i}_{x}" for i, x in zip(vid, a1)], "variant_ids": vid, "a1": a1}


_GSTR = np.array(["NA", "0", "1", "2"], dtype=object)


def write_outputs(i, res, u, analysed_idx, pa, out_dir):
    car, counted = res["car"], res["counted_car"]
    # same layout as plink 1.9 --recode A (.raw), carriers only. Genotype cells are single digits
    # (or NA), so each row's genotype part is built as bytes in one numpy step.
    memo = pa["head"]
    new = [k for k in car if memo[k] is None]
    if new:
        w = analysed_idx[new]
        for k, x in zip(new, zip(*(u.famcol[c][w] for c in ("FID", "IID", "PAT", "MAT", "SEX")),
                                 [str(v) for v in pa["str"][new]])):
            memo[k] = "\t".join(x) + "\t"
    head = memo[car]
    M = counted.T
    width = 2 * M.shape[1]
    body = np.full((M.shape[0], width), 9, dtype=np.uint8)              # tabs
    body[:, 0::2] = (M + 48).astype(np.uint8)                            # '0' '1' '2'
    body[:, -1] = 10
    raw = body.tobytes()
    missing = set(np.flatnonzero((M < 0).any(axis=1)).tolist())
    parts = [("\t".join(["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"] + list(res["cols"])) + "\n").encode()]
    for r, h in enumerate(head):
        if r in missing:
            parts.append((h + "\t".join(_GSTR[M[r].astype(np.int16) + 1]) + "\n").encode())
        else:
            parts.append(h.encode() + raw[r * width:(r + 1) * width])
    cf = os.path.join(out_dir, f"carriers_aggregate_{i}.txt")
    with open(cf, "wb") as f:
        f.write(b"".join(parts))

    pheno_num = pa["num"]
    pheno_car = pheno_num[car]
    is_c, is_hom = counted >= 1, counted == 2
    pv = pd.DataFrame({"variant": res["variant_ids"], "counted_allele": res["a1"]})
    if pa["binary"]:
        case_c, ctrl_c = pheno_car == 2, pheno_car == 1
        counts = {"n_cases": pa["n_cases"], "n_controls": pa["n_controls"],
                  "carriers_cases": int(case_c.sum()), "carriers_controls": int(ctrl_c.sum())}
        counts["carrier_pct_cases"] = round(100 * counts["carriers_cases"] / max(counts["n_cases"], 1), 4)
        counts["carrier_pct_controls"] = round(100 * counts["carriers_controls"] / max(counts["n_controls"], 1), 4)
        pv["carriers_cases"] = is_c[:, case_c].sum(axis=1)
        pv["carriers_controls"] = is_c[:, ctrl_c].sum(axis=1)
        pv["hom_cases"] = is_hom[:, case_c].sum(axis=1)
        pv["hom_controls"] = is_hom[:, ctrl_c].sum(axis=1)
    else:
        non = np.ones(len(pheno_num), dtype=bool)
        non[car] = False
        counts = {"n_analysed": len(pheno_num),
                  "mean_pheno_carriers": float(np.nanmean(pheno_car)) if len(car) else np.nan,
                  "mean_pheno_noncarriers": float(np.nanmean(pheno_num[non]))}
        pv["carriers"] = is_c.sum(axis=1)
        pv["hom"] = is_hom.sum(axis=1)
    vf = os.path.join(out_dir, f"variant_counts_aggregate_{i}.tsv")
    pv.to_csv(vf, sep="\t", index=False, na_rep="NA")
    return len(car), counts, [cf, vf]


EMPTY = "FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE\n"


def process_chromosome(cfg, CHR, group, pheno, keep, khash, dry_run):
    t_chr = time.time()
    work = os.path.join(cfg.work, f"chr{CHR}")
    out_local = os.path.join(cfg.work, "output")
    os.makedirs(work, exist_ok=True)
    os.makedirs(out_local, exist_ok=True)
    range_units = set(cfg.inp["range_units"])
    recs, files = [], []

    def record(i, row, **kw):
        r = {"aggregate": i, "ID": row.ID, "CHR": CHR, "unit": row.unit, "region": row.region,
             "annotation": row.annotation, "singleton": row.singleton}
        r.update(kw)
        recs.append(r)
        return r

    # --- step A for every aggregate on this chromosome (one pass over the annotation file, or the cache)
    t = time.time()
    try:
        keys = {k for k in (region_key(r, range_units) for r in group.itertuples()) if k}
        masks = mask_table(cfg, CHR, os.path.join(cfg.cache, "annotations")) if keys else {}
        ann = annotation_rows(cfg, CHR, keys, work) if keys else None
    except Exception as e:
        for i, row in group.iterrows():
            record(i, row, status=f"FAILED annotation/mask files: {str(e).splitlines()[0]}")
        return recs, files
    t_ann = time.time() - t
    timed("annotation", t)
    ai = annotation_index(ann)
    steps = {}
    for i, row in group.iterrows():
        try:
            steps[i] = step_a(row, ai, masks, range_units)
        except Exception as e:
            steps[i] = e
    if dry_run:
        for i, row in group.iterrows():
            a = steps[i]
            if isinstance(a, Exception):
                record(i, row, n_variants_annotated=None, sets="", status=f"FAILED: {a}")
            else:
                record(i, row, n_variants_annotated="range" if a["mode"] == "range" else len(a["variants"]),
                       sets=a["sets"], status="dry run")
        log(f"chr{CHR}: {len(group)} aggregates, annotation lookup {t_ann:.1f}s")
        return recs, files

    # --- genotypes: from the local store; plink2 only for variants it doesn't have yet
    t = time.time()
    try:
        ids = [v for a in steps.values() if isinstance(a, dict) for v in a["variants"]]
        windows = [a["range"] for a in steps.values() if isinstance(a, dict) and a["mode"] == "range"]
        unions, how = genotypes_for(cfg, CHR, ids, windows, keep, khash, work)
    except Exception as e:
        for i, row in group.iterrows():
            record(i, row, status=f"FAILED genotype files: {str(e).splitlines()[0]}")
        log(f"chr{CHR}: ERROR {e}")
        return recs, files
    some = next(iter(unions.values()), None)
    if some is not None:
        iids = some.fam["IID"].values
        j = pd.Index(pheno["iid"]).get_indexer(iids)
        found = j >= 0
        num = np.where(found, pheno["num"][np.maximum(j, 0)], np.nan)
        analysed = found & (pd.Index(list(keep)).get_indexer(iids) >= 0) & ~np.isnan(num) & (num != -9)
        analysed_idx = np.flatnonzero(analysed)
        pa = {"str": pheno["str"][j[analysed_idx]], "num": num[analysed_idx]}
        pa["binary"] = set(np.unique(pa["num"])) <= {1.0, 2.0}
        pa["n_cases"], pa["n_controls"] = int((pa["num"] == 2).sum()), int((pa["num"] == 1).sum())
        pa["head"] = np.full(len(analysed_idx), None, dtype=object)     # FID..PHENOTYPE text per person
    t_geno = time.time() - t
    timed("genotypes", t)
    log(f"chr{CHR}: annotation lookup {t_ann:.1f}s; genotypes {t_geno:.1f}s ({how})"
        + (f"; {len(analysed_idx):,} people analysed" if some is not None else ""))

    # --- per aggregate: numpy filters + carriers
    t_car = time.time()
    for i, row in group.iterrows():
        t = time.time()
        a = steps[i]
        cf = os.path.join(out_local, f"carriers_aggregate_{i}.txt")
        try:
            if isinstance(a, Exception):
                raise a
            n_ann = "range" if a["mode"] == "range" else len(a["variants"])
            res = carriers_for(row, a, unions, analysed_idx, cfg.c["max_maf"]) if some is not None else None
            if a["mode"] == "ids" and not a["variants"]:
                open(cf, "w").write(EMPTY)
                rec = record(i, row, n_variants_annotated=0, n_variants_found=0, n_variants_after_freq_filter=0,
                             n_carriers=0, status="no variants annotated")
            elif res is None:
                open(cf, "w").write(EMPTY)
                rec = record(i, row, n_variants_annotated=n_ann, n_variants_found=0, n_variants_after_freq_filter=0,
                             n_carriers=0, status="no variants in genotype files")
            elif res["n_after"] == 0:
                open(cf, "w").write(EMPTY)
                rec = record(i, row, n_variants_annotated=n_ann, n_variants_found=res["n_found"],
                             n_variants_after_freq_filter=0, n_carriers=0, status="no variants after frequency filter")
            else:
                n_car, counts, fs = write_outputs(i, res, res["u"], analysed_idx, pa, out_local)
                files += fs
                rec = record(i, row, n_variants_annotated=n_ann, n_variants_found=res["n_found"],
                             n_variants_after_freq_filter=res["n_after"], n_carriers=n_car, **counts, status="OK")
            if cf not in files:
                files.append(cf)
        except Exception as e:
            rec = record(i, row, status=f"FAILED: {str(e).splitlines()[0]}")
        rec["seconds"] = round(time.time() - t, 3)
        split = (f" (cases {rec['carriers_cases']}/{rec['n_cases']}, controls {rec['carriers_controls']}/{rec['n_controls']})"
                 if "carriers_cases" in rec else "")
        log(f"   aggregate {i} {row.ID}: {rec['status']}, {rec.get('n_variants_found', 0)} variants, "
            f"{rec.get('n_carriers', 0)} carriers{split} [{rec['seconds']:.2f}s]")
    timed("carriers", t_car)
    shutil.rmtree(os.path.join(cfg.work, "geno_copy"), ignore_errors=True)
    shutil.rmtree(work, ignore_errors=True)
    log(f"chr{CHR}: done in {time.time() - t_chr:.1f}s")
    return recs, files


def step_extract(cfg, aggregates, dry_run=False):
    t = time.time()
    pheno, keep, khash = load_samples(cfg)
    timed("samples", t)
    log(f"extract: {len(aggregates)} aggregates on chromosomes "
        f"{sorted(aggregates['CHR'].astype(str).unique(), key=chrom_order)}; {len(keep):,} people in the sample file")
    shutil.rmtree(os.path.join(cfg.work, "output"), ignore_errors=True)
    chroms = sorted(aggregates["CHR"].astype(str).unique(), key=chrom_order)
    groups = [(ch, aggregates[aggregates["CHR"].astype(str) == ch]) for ch in chroms]
    n_par = max(1, int(cfg.perf["parallel_chromosomes"]))
    if n_par > 1:
        with ThreadPoolExecutor(max_workers=n_par) as ex:
            results = list(ex.map(lambda g: process_chromosome(cfg, g[0], g[1], pheno, keep, khash, dry_run), groups))
    else:
        results = [process_chromosome(cfg, ch, g, pheno, keep, khash, dry_run) for ch, g in groups]
    recs = [r for rs, _ in results for r in rs]
    files = [f for _, fs in results for f in fs]
    summary = pd.DataFrame(recs).sort_values("aggregate")
    summary = summary[[x for x in summary.columns if x != "status"] + ["status"]]
    name = f"dry_run_{cfg.prefix}.tsv" if dry_run else f"extraction_summary_{cfg.prefix}.tsv"
    out_local = os.path.join(cfg.work, "output")
    os.makedirs(out_local, exist_ok=True)
    sf = os.path.join(out_local, name)
    summary.to_csv(sf, sep="\t", index=False, na_rep="NA")
    t = time.time()
    publish(files + [sf], cfg.carriers_out)
    timed("upload", t)
    log(f"extract: wrote {len(files) + 1} files to {cfg.carriers_out} ({time.time() - t:.1f}s)")
    with pd.option_context("display.max_columns", 30, "display.width", 250, "display.max_colwidth", 60):
        print(summary.drop(columns=[x for x in ("sets",) if x in summary.columns and not dry_run]).to_string(index=False))
    if summary["status"].astype(str).str.startswith("FAILED").any():
        return 1
    return 0


# =====================================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("steps", nargs="*", default=["all"],
                    help="all (default) | merge | index | samples | select | extract | dry-run")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=JSON",
                    help='override a config value, e.g. --set selection.name=\'"test"\' --set selection.log10p_min=null')
    args = ap.parse_args()

    overrides = {}
    for kv in args.set:
        key, _, val = kv.partition("=")
        d = overrides
        parts = key.split(".")
        for p in parts[:-1]:
            d = d.setdefault(p, {})
        d[parts[-1]] = json.loads(val)
    cfg = Cfg(args.config, overrides)
    steps = set(args.steps)
    log(f"config {args.config}: {cfg.prefix}; steps {sorted(steps)}")
    os.makedirs(cfg.cache, exist_ok=True)
    MNT.update(cache=cfg.cache, use=bool(cfg.perf["use_mounts"]), threads=cfg.perf["upload_threads"],
               out_bucket=split_gs(cfg.out_dir)[0] if is_gs(cfg.out_dir) else None)
    if "genotypes" in cfg.redo:
        shutil.rmtree(os.path.join(cfg.cache, "genotypes"), ignore_errors=True)
        shutil.rmtree(os.path.join(cfg.cache, "genotype_index"), ignore_errors=True)
    if "annotations" in cfg.redo:
        shutil.rmtree(os.path.join(cfg.cache, "annotations"), ignore_errors=True)

    full = "all" in steps
    selecting = full or bool(steps & {"select", "extract", "dry-run"})
    t = time.time()
    merged_now = False
    need_merge = "merge" in cfg.redo or (not index_ready(cfg) and not exists(cfg.merged))
    if (full and need_merge) or "merge" in steps:
        step_merge(cfg)
        merged_now = True
    if "index" in steps or (selecting and not merged_now and
                            ("index" in cfg.redo or not index_ready(cfg) or index_stale(cfg))):
        step_index(cfg)
    timed("merge/index", t)
    local_meta = [os.path.join(cfg.meta_local, f) for f in ("phenotype_file.txt", "sample_file_no_NA.txt")]
    meta_done = all(os.path.exists(f) for f in local_meta) or \
        (exists(f"{cfg.meta_out}/phenotype_file.txt") and exists(f"{cfg.meta_out}/sample_file_no_NA.txt"))
    if (full and ("samples" in cfg.redo or not meta_done)) or "samples" in steps:
        step_samples(cfg)
    rc = 0
    if selecting:
        agg = step_select(cfg)
        if full or "extract" in steps:
            rc = step_extract(cfg, agg)
        elif "dry-run" in steps:
            rc = step_extract(cfg, agg, dry_run=True)
    t = time.time()
    for f in _PENDING:
        f.result()
    timed("upload", t)
    shown = [f"{k} {v:.1f}s" for k, v in TIMES.items() if v >= 0.05]
    log(f"finished; time by stage: {', '.join(shown) or '-'}")
    sys.exit(rc)


if __name__ == "__main__":
    main()

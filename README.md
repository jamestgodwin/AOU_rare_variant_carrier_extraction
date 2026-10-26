# Fast rare-variant aggregate carrier extraction (All of Us v9)
### James Godwin, Marina Vabistsevits (October 2026, v1)

> *This repository adapts the UK Biobank carrier-extraction pipeline ([UKB_extract_carriers_rare_variant_aggregates](https://github.com/mvab/UKB_extract_carriers_rare_variant_aggregates)) to run on the All of Us Researcher Workbench, and combines its steps into a single, faster command.*

After a regenie rare-variant burden analysis on the All of Us Researcher Workbench, this tool shows **who carries the variants behind an aggregate and what their phenotype is**. For each aggregate it writes:

- the carriers and their genotypes, one column per variant
- carriers per variant, split by cases and controls
- a summary row: variants used, number of carriers, and carrier % in cases vs controls (for quantitative traits, the mean phenotype in carriers vs non-carriers)

---

## Setup

In a JupyterLab terminal on the Workbench:

```bash
git clone https://github.com/<your-account>/rare-variant-carrier-extraction.git ~/carrier_pipeline
cd ~/carrier_pipeline
chmod +x run_carrier_extraction.sh
cp config.example.json config.json
```

Then edit `config.json` as described in the next two sections. You need `python3` (numpy, pandas), `plink2` and `gcloud`, which the Workbench already has.

## 1. Paths you must change (`config.json` → `inputs` / `outputs`)

| key | point it to | example path |
|---|---|---|
| `inputs.results_dir` | folder with your regenie step 2 output files (`*.regenie.gz`) | |
| `inputs.pheno_file` | phenotype file (`FID IID <phenotype>`; binary traits coded 0/1 or 1/2) | |
| `inputs.cov_file` | covariate file. Only its IDs are used, to define who was analysed | |
| `inputs.ancestry_file` | All of Us `ancestry_preds.tsv` (or `null` for no ancestry filter) | `gs://marina-bucket-wb-twinkly-broccoli-3877/PCOS_and_hormones_AllofUs/data/ancestry_preds.tsv` |
| `inputs.anno_file` | the regenie `--anno-file`, with `{CHR}` where the chromosome goes, e.g. `.../annotations_all_chr{CHR}.zst` | `gs://pgen-bucket-wb-cordial-banana-2765/wgs_annotations/v2026.10/regenie_combined/annotations_all_chr{CHR}.zst` |
| `inputs.mask_file` | the regenie `--mask-def` file, e.g. `.../masks_all_chr{CHR}.zst` | `gs://pgen-bucket-wb-cordial-banana-2765/wgs_annotations/v2026.10/regenie_combined/masks_all_chr{CHR}.zst` | 
| `inputs.geno_prefix` | genotype files regenie step 2 used, without extension, e.g. `.../chr{CHR}_genoqc` for `chr{CHR}_genoqc.pgen/.pvar/.psam` | `gs://pgen-bucket-wb-cordial-banana-2765/wgs_v9_pgens_genoqc/chr{CHR}_genoqc_trimmed_merged` |
| `outputs.out_dir` | where the results go: a folder in **your workspace bucket** | |

The regenie result files must be named the way regenie wrote them in the original workflow:

```
{run}_{ancestry}_chr{N}_{start}_{end}_<anything>_{analysis}_{run}.regenie.gz
```

## 2. Parameters you must change

| key | set it to |
|---|---|
| `run` | the run name in your result files, e.g. `MYTRAIT_v9` |
| `ancestry` | the ancestry to analyse, as in the file names and `ancestry_preds.tsv`, e.g. `eur` |
| `analysis` | the analysis label in the file names, e.g. `burden` (or `mac10`) |
| `max_maf` | the MAF limit of the aggregates you extract, e.g. `0.001` for `….0.001` aggregates. It's one value per run, so extract one frequency bin at a time. `singleton` aggregates always use MAC ≤ 1. |
| `selection.name` | a label for this search. Outputs go to `output/{run}_{ancestry}_{analysis}_{name}/`, so change it for each search or the previous one is overwritten |
| `selection.targets` | what to look up (see below) |

Check these too, but only change them if they don't fit:

- `inputs.pheno_col`: the phenotype column, when the phenotype file has several (default: the first one)
- `inputs.geno_format`: `pgen` (default), `bed` or `bgen`
- `inputs.chr_x`: how chromosome X is written in file names
- `inputs.range_units`: the unit names of your sliding-window aggregates (default `["sliding"]`)

## 3. Run

```bash
./run_carrier_extraction.sh config.json dry-run    # preview: which aggregates are selected, how many variants each
./run_carrier_extraction.sh config.json            # extract carriers
```

The first run takes longer (it indexes your results once, and caches files the first time each chromosome is used); after that, searches take seconds. Results appear in `out_dir/output/{run}_{ancestry}_{analysis}_{name}/`:

```
carriers_aggregate_N.txt          carriers of aggregate N: FID IID ... PHENOTYPE + one column per variant (copies of the minor allele)
variant_counts_aggregate_N.tsv    carriers per variant, in cases and controls
extraction_summary_*.tsv          one row per aggregate: variants, carriers, cases/controls, carrier %
```

`N` matches the row of `out_dir/input/aggregates_{...}.csv`, which lists the selected aggregates with their p-values.

---

## Searching genes of interest

Set `selection.targets` in `config.json`, or override it for one run with `--set`. Values are JSON, so strings need the inner quotes.

```bash
# every aggregate of one gene (all masks, all units)
./run_carrier_extraction.sh config.json --set selection.name='"CHEK2"' --set selection.targets='["CHEK2"]'

# a list of genes, loss-of-function masks only
./run_carrier_extraction.sh config.json \
    --set selection.name='"DNA_repair_lof"' \
    --set selection.targets='["CHEK2", "BRCA2", "ATM"]' \
    --set selection.masks='["combined_damaging_lof"]'

# one exact aggregate (as written in the regenie output)
./run_carrier_extraction.sh config.json --set selection.name='"CHEK2_lof"' \
    --set selection.targets='["coding_CHEK2_ENST00000402731.combined_damaging_lof.0.001"]'

# a long gene list from a file (one gene per line)
./run_carrier_extraction.sh config.json --set selection.name='"my_genes"' --set selection.targets_file='"genes.txt"'

# no genes: all aggregates with p < 1e-6
./run_carrier_extraction.sh config.json --set selection.name='"top_hits"' \
    --set selection.targets='[]' --set selection.log10p_min=6
```

What a target can be:

- **gene symbol**, e.g. `CHEK2`, for all of its aggregates (case doesn't matter)
- **transcript**, e.g. `ENST00000402731`
- **full aggregate ID**, which is always kept, whatever the other filters say
- **ENSR ID** or **coordinates** (see the next section)

Narrow the search with these optional filters:

- `selection.masks`: mask names as in your mask file, e.g. `combined_damaging_lof` (any frequency) or `combined_damaging_lof.0.001`
- `selection.units`: e.g. `coding`
- `selection.log10p_min`: keep only aggregates with −log10(p) ≥ this value

Run `dry-run` first to see what a search selects; any target that matches nothing is reported. To find the mask names, dry-run a gene without `masks` and look at the `annotation` column of `out_dir/input/aggregates_*.csv`.

## Non-coding variants

Non-coding aggregates are selected the same way as genes. You name the gene, the regulatory element or a region, and narrow it down with the *unit* (the first part of the aggregate ID). The annotation file used here has three non-coding units:

| unit | what it is |
|---|---|
| `proximal` | promoter / upstream variants of a gene |
| `regulatory` | regulatory elements (ENSR IDs, or regions) |
| `rna` | non-coding RNA genes |

| what you want | `targets` | `units` |
|---|---|---|
| promoter / upstream variants of a gene | `["CHEK2"]` | `["proximal"]` |
| a non-coding RNA gene | `["<RNA gene symbol or transcript>"]` | `["rna"]` |
| one regulatory element | `["ENSR00000012345"]` | – |
| all regulatory elements in a region | `["chr22:28650000-28800000"]` | `["regulatory"]` |
| all non-coding aggregates in a region | `["chr22:28650000-28800000"]` | `["proximal", "regulatory", "rna"]` |
| everything (coding + non-coding) in a region | `["chr22:28650000-28800000"]` | – |

```bash
# all non-coding aggregates around CHEK2 (GRCh38 coordinates)
./run_carrier_extraction.sh config.json \
    --set selection.name='"CHEK2_noncoding"' \
    --set selection.targets='["chr22:28650000-28800000"]' \
    --set selection.units='["proximal", "regulatory", "rna"]'
```

- **Coordinates** (`chr22:start-end` or `22:start-end`, GRCh38) select every aggregate whose region overlaps the range, or whose position lies inside it. This is how to find regulatory elements near a gene, since their IDs don't contain the gene name.
- **Masks:** any non-coding mask defined in your mask file can be added with `selection.masks`.
- **Frequency:** non-coding aggregates are often tested at a higher frequency bin, so set `max_maf` to match (e.g. `0.01`).
- **Check first:** run a `dry-run` to confirm what a search selects before extracting. `rna` aggregates are only recognised when their mask name contains "rna".

## Outputs

All in `out_dir`. `PREFIX` = `{run}_{ancestry}_{analysis}_{selection.name}`, and `N` = the aggregate's row (starting at 0) in the aggregates file.

| file | contents |
|---|---|
| `input/aggregates_{PREFIX}.csv` | the selected aggregates, one row each: chromosome, position, unit, gene, region, mask, frequency bin, MAC, LOG10P, BETA, aggregate ID |
| `output/{PREFIX}/extraction_summary_{PREFIX}.tsv` | one row per aggregate: variants annotated → found in the genotypes → kept after the frequency filter; number of carriers; cases/controls analysed, carriers in each and carrier %; status (`OK` or why not) |
| `output/{PREFIX}/variant_counts_aggregate_N.tsv` | one row per kept variant: the counted (minor) allele, carriers and homozygotes in cases and in controls |
| `output/{PREFIX}/carriers_aggregate_N.txt` | carriers only, in plink `.raw` layout: `FID IID PAT MAT SEX PHENOTYPE` (2 = case, 1 = control), then one column per variant holding copies of the minor allele (0/1/2/NA) |
| `output/{PREFIX}/dry_run_{PREFIX}.tsv` | dry-run only: variants annotated and annotation sets used per aggregate |
| `metadata_{ancestry}/phenotype_file.txt`, `sample_file_no_NA.txt` | the phenotype (binary recoded to 1/2) and the people analysed |
| `{run}_{ancestry}_{analysis}_aggregates_merged.tsv.gz` | all regenie results for this run, ancestry and analysis, merged |

For continuous traits, the summary shows `n_analysed` and the mean phenotype in carriers vs non-carriers instead of case/control columns. The variant counts show `carriers` and `hom`, and `PHENOTYPE` holds the measured value.

Carrier files are participant-level data, and counts of 1–20 must not leave the Workbench.

---

## Tips

- **Data protection:** outputs are participant-level data. Keep them on the Workbench, and never share counts of 1–20 participants outside it (All of Us Data and Statistics Dissemination Policy). `.gitignore` keeps `config.json` and all output files out of this repository.
- **Caching:** everything is cached under `~/carrier_cache`, which is why repeat searches are fast. Delete it to start fresh.
- **Logs:** every run is logged in `logs/`, and the last line shows the time per stage.

#!/bin/bash
# Carrier extraction for regenie rare-variant aggregates (All of Us Researcher Workbench).
#
# Usage:
#   ./run_carrier_extraction.sh config.json                 # everything (steps already done are skipped)
#   ./run_carrier_extraction.sh config.json dry-run         # select + annotation lookup only, no genotypes
#   ./run_carrier_extraction.sh config.json select          # just (re)build the aggregates file
#   ./run_carrier_extraction.sh config.json merge index     # force specific steps
#
# Any config value can be overridden for one run with --set (value in JSON), e.g.
#   ./run_carrier_extraction.sh config.json --set selection.name='"CHEK2"' \
#       --set selection.targets='["CHEK2"]' --set selection.log10p_min=null
#
# Steps: merge -> index -> samples -> select -> extract (see README.md).
# A log of every run is written to logs/ next to this script.

set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "$#" -lt 1 ] || [ ! -f "$1" ]; then
    sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
    if [ "$#" -ge 1 ] && [ ! -f "$1" ] && [ -f "$here/config.example.json" ]; then
        echo "'$1' not found. Start from the example:  cp config.example.json $1  (then fill in your paths)"
    fi
    exit 1
fi
config="$1"; shift

for tool in python3 plink2 gcloud; do
    command -v "$tool" > /dev/null || { echo "ERROR: $tool is not on PATH"; exit 1; }
done
command -v zstd > /dev/null || echo "note: zstd not found; plink2 will decompress the .zst annotation files (slower)"

mkdir -p "$here/logs"
log_file="$here/logs/run_$(date +%Y%m%d_%H%M%S).log"
echo "log: $log_file"

python3 -u "$here/carrier_extraction.py" "$config" "$@" 2>&1 | tee "$log_file"
exit "${PIPESTATUS[0]}"

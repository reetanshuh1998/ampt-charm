#!/usr/bin/env bash
# test_results.sh : regenerate every table and figure of observation/results/ from the shipped analysis outputs
# (observation/results/json/) and check the tables are identical to the shipped ones.  ~30 s, no event generation.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
T=$(mktemp -d "${TMPDIR:-/tmp}/ampt_results_XXXX")
trap 'rm -rf "$T"' EXIT
cd "$ROOT"
python3 -I analysis/audit_claims.py "$T/AUDIT_TABLES.md" > /dev/null
if diff -q "$T/AUDIT_TABLES.md" observation/results/tables/AUDIT_TABLES.md > /dev/null; then
  echo "PASS  AUDIT_TABLES.md (all chi2 of this work, public AMPT and the published models) regenerated identically"
else
  diff "$T/AUDIT_TABLES.md" observation/results/tables/AUDIT_TABLES.md | head -20; echo "FAIL  AUDIT_TABLES.md differs"; exit 1
fi
python3 -I analysis/paper_figs.py "$T/figures" > /dev/null
python3 -I analysis/fig_tilt_law.py "$T/figures/fig_tilt_law" > /dev/null
BULK_V2=sp python3 -I analysis/compare_bulk.py "$T/bulk" observation/results/json/dm25_main.json \
  observation/results/json/dm25_pub.json > /dev/null 2>&1
diff -q "$T/bulk/BULK_FLOW.md" observation/results/bulk/BULK_FLOW.md > /dev/null \
  && echo "PASS  BULK_FLOW.md regenerated identically" || { echo "FAIL  BULK_FLOW.md differs"; exit 1; }
python3 -I analysis/pp_shape_json.py ours=observation/results/json/dm25_ppmain.json \
  public=observation/results/json/dm25_pppub.json > "$T/pp.txt"
diff -q "$T/pp.txt" observation/results/tables/PP_SHAPE.txt > /dev/null \
  && echo "PASS  PP_SHAPE.txt regenerated identically" || { echo "FAIL  PP_SHAPE.txt differs"; exit 1; }
echo "PASS  $(ls "$T/figures"/*.pdf | wc -l) figures regenerated in $T/figures (compare with observation/results/figures)"

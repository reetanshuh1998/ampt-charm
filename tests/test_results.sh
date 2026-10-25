#!/usr/bin/env bash
# test_results.sh : regenerate every table and figure of observation/results/ from the shipped analysis outputs
# (observation/results/json/) and check the tables are identical to the shipped ones.  ~30 s, no event generation.
# Output: regenerated/ (tables, figures; compare the figures with observation/results/figures/).
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
T=${RESULTS_DIR:-$ROOT/regenerated}
rm -rf "$T" && mkdir -p "$T"
cd "$ROOT"
J=observation/results/json
same() {   # same NEW SHIPPED LABEL
  if diff -q "$1" "$2" > /dev/null; then echo "PASS  $3 regenerated identically"
  else diff "$1" "$2" | head -20; echo "FAIL  $3 differs"; exit 1; fi
}
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
python3 -I analysis/tilt_shadow.py "$T/TILT_DM25.md" "dm24 600k=$J/dm24_e_ab.json" "dm25 2.4M=$J/dm25_main.json" \
  "public AMPT=$J/dm25_pub.json" > /dev/null
same "$T/TILT_DM25.md" observation/results/tables/TILT_DM25.md "TILT_DM25.md (charm tilt, v1(c), D0 dv1/dy)"
python3 -I analysis/eps_hessian.py "$T/EPS09S_BAND_DM25.md" $J/dm25_ppmain.json $J/eps/dm25_main__eps*.json.gz > /dev/null
same "$T/EPS09S_BAND_DM25.md" observation/results/tables/EPS09S_BAND_DM25.md "EPS09S_BAND_DM25.md (EPS09s Hessian band, 31 sets)"
python3 -I - "$T/figures" observation/results/figures <<'EOF'
import glob, os, sys
import matplotlib.image as mi
import numpy as np
new, ref = sys.argv[1:3]
bad = [os.path.basename(f) for f in sorted(glob.glob(f'{ref}/*.png'))
       if not os.path.exists(f'{new}/{os.path.basename(f)}')
       or mi.imread(f).shape != mi.imread(f'{new}/{os.path.basename(f)}').shape
       or np.abs(mi.imread(f) - mi.imread(f'{new}/{os.path.basename(f)}')).max() > 0]
n = len(glob.glob(f'{ref}/*.png'))
print(f'PASS  {n} figures regenerated pixel-identical in {new}' if not bad else
      f'WARN  {len(bad)} of {n} figures differ in pixels (fonts or matplotlib version?): {" ".join(bad)}; compare by eye')
EOF

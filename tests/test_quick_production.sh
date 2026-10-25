#!/usr/bin/env bash
# test_quick_production.sh [NAA_JOBS] [NAA_EVENTS] [NPP_JOBS] [NPP_EVENTS]
# A laptop-size production of the reference configuration with NEW seeds (default 4 x 20 Au+Au, 2 x 5000 p+p events,
# ~5 min on 4 cores), then:
#   1. the production analysis (analysis/analyze_arm.py)
#   2. an independent re-implementation (analysis/indep_check.py) on the same events: must agree exactly
#   3. statistical comparison with the published full-statistics results (tests/compare_with_published.py)
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
NA=${1:-4} EA=${2:-20} NP=${3:-2} EP=${4:-5000}
W=${QUICK_DIR:-$ROOT/quick}
rm -rf "$W" && mkdir -p "$W"
echo ">> Au+Au: $NA jobs x $EA events"
"$ROOT/scripts/run_arm.sh" reference_auau quick_aa auau 0 14.5 "$EA" 99000001 17 3 1 "$NA" "" "$W"
echo ">> p+p: $NP jobs x $EP events"
"$ROOT/scripts/run_arm.sh" reference_pp quick_pp pp 0 2.0 "$EP" 99000002 17 5 1 "$NP" "" "$W"
J=$(nproc 2>/dev/null || echo 2)
python3 -I "$ROOT/analysis/analyze_arm.py" "$W/quick_aa" "$W/quick_aa.json" --sys auau -j "$J" > "$W/analyze_aa.log"
python3 -I "$ROOT/analysis/analyze_arm.py" "$W/quick_pp" "$W/quick_pp.json" --sys pp -j "$J" > "$W/analyze_pp.log"
python3 -I "$ROOT/analysis/indep_check.py" "$W/quick_aa" "$W/quick_pp" "$W/indep.json" > "$W/indep.log" 2>&1
python3 -I "$ROOT/analysis/indep_compare.py" "$W/indep.json" "$W/quick_aa.json" "$W/quick_pp.json" "$W/INDEP_CHECK.md" 2>/dev/null \
  | tail -n +5
python3 -I - "$W/INDEP_CHECK.md" <<'EOF'
import sys
bad = [l for l in open(sys.argv[1]) if l.startswith('| ') and not l.startswith('| observable') and '---' not in l
       and l.split('|')[3].strip() not in ('nan', '0.00e+00') and float(l.split('|')[3]) > 1e-12]
print('PASS  independent re-analysis agrees exactly' if not bad else 'FAIL  independent re-analysis differs:\n' + ''.join(bad))
sys.exit(1 if bad else 0)
EOF
python3 -I "$ROOT/tests/compare_with_published.py" "$W/quick_aa.json" "$W/quick_pp.json"

#!/usr/bin/env bash
# eps_band.sh ARMDIR PP.json OUTDIR [NJOBS] : EPS09s NLO Hessian band (90 % CL) of the paper observables for one
# Au+Au sample of the reference configuration (AMPT_HQPDF=1, AMPT_HQSHW=1; README Sec. 2.3).
#   1. analysis/eps09s/eps09s_errw.py: per-event weights of the 31 error sets -> ARMDIR/hqerr_*.txt.gz
#   2. analysis/analyze_arm.py --shw tab:K, K = 1..31 -> OUTDIR/NAME__epsKK.json (same events, re-weighted)
#   3. analysis/eps_hessian.py -> OUTDIR/EPS09S_BAND.md
# PP.json is the analysis of the p+p reference (unweighted, the same for every set).
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
A=${1:?ARMDIR} PP=${2:?PP.json} O=${3:?OUTDIR}
J=${4:-$(nproc 2>/dev/null || echo 2)}
N=$(basename "$A")
mkdir -p "$O"
[ -x "$ROOT/bin/eps09s_rg" ] || make -C "$ROOT" eps-tool > /dev/null
nred=$(ls "$A"/red_*.dat.gz | wc -l); nerr=$(ls "$A"/hqerr_*.txt.gz 2>/dev/null | wc -l)
if [ "$nerr" -eq "$nred" ]; then
  echo "using the existing weights ($nerr hqerr files)"
else
  python3 -I "$ROOT/analysis/eps09s/eps09s_errw.py" "$A" > "$O/epserr_$N.log"
  grep '^# closure' "$O/epserr_$N.log"
fi
for k in $(seq -w 1 31); do
  python3 -I "$ROOT/analysis/analyze_arm.py" "$A" "$O/${N}__eps$k.json" --sys auau -j "$J" --shw "tab:$((10#$k))" \
    > "$O/analyze_${N}__eps$k.log"
done
python3 -I "$ROOT/analysis/eps_hessian.py" "$O/EPS09S_BAND.md" "$PP" "$O/${N}"__eps??.json > /dev/null
echo "wrote $O/EPS09S_BAND.md"

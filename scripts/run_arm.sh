#!/usr/bin/env bash
# run_arm.sh CONFIG ARM SYS BMIN BMAX NEV SEEDBASE SEEDPMUL SEEDPOFF FIRST LAST [NPAR] [OUTDIR]
# Runs jobs FIRST..LAST of one arm in parallel (NPAR at a time, default: number of cores); logs in OUTDIR/logs/ARM/.
# Example (10 jobs of 30 events of the reference Au+Au configuration, seeds of production arm dm24_e_a):
#   scripts/run_arm.sh reference_auau dm24_e_a auau 0 14.5 30 20261107 13 9 1 10
set -euo pipefail
[ $# -ge 11 ] || { sed -n 2,5p "$0"; exit 1; }
ROOT=$(cd "$(dirname "$0")/.." && pwd)
NPAR=${12:-$(nproc 2>/dev/null || echo 2)}
OUTDIR=${13:-$ROOT/prod}
mkdir -p "$OUTDIR/logs/$2"
seq "${10}" "${11}" | xargs -P "$NPAR" -I{} sh -c \
  "'$ROOT/scripts/run_job.sh' $1 $2 $3 $4 $5 $6 $7 $8 $9 {} '$OUTDIR' > '$OUTDIR/logs/$2/{}.log' 2>&1 \
   && echo 'job {} ok' || echo 'job {} FAILED (see $OUTDIR/logs/$2/{}.log)'"

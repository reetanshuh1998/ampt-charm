#!/usr/bin/env bash
# run_job.sh CONFIG ARM SYS BMIN BMAX NEV SEEDBASE SEEDPMUL SEEDPOFF IDX [OUTDIR]
#
# One production job, exactly as on the batch farm (KEKCC lsf_job_dm24/dm25.sh):
#   HIJING seed = SEEDBASE + IDX*7919,  ZPC seed = IDX*SEEDPMUL + SEEDPOFF
#   input card  = inputs/$INPUT of CONFIG with NEVNT, BMIN, BMAX and the two seeds set (lines 9-11, 28-30)
#   environment = $DMENV of CONFIG (the physics switches, docs/CHANGES.pdf)
# Writes OUTDIR/ARM/red_ARM_IIII.dat.gz (reduced ampt.dat, scripts/reduce_dmeson_v2.awk), charm_ARM_IIII.dat.gz
# (ana/charm.dat) and hq_ARM_IIII.log.gz.  SYS = auau | pp.  OUTDIR defaults to prod/.
set -euo pipefail
[ $# -ge 10 ] || { sed -n 2,10p "$0"; exit 1; }
CONFIG=$1 ARM=$2 SYS=$3 BMIN=$4 BMAX=$5 NEV=$6 SEEDBASE=$7 SEEDPMUL=$8 SEEDPOFF=$9 IDX=${10}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUTDIR=${11:-$ROOT/prod}
INPUT='' DMENV=''
# shellcheck disable=SC1090
source "$ROOT/configs/${CONFIG%.env}.env"
AWK=$(command -v gawk || true)
[ -n "$AWK" ] || { echo "ERROR: gawk is required (the light-hadron prescale uses gawk's rand())" >&2; exit 1; }
[ -x "$ROOT/bin/ampt" ] || { echo "ERROR: bin/ampt missing, run scripts/build.sh" >&2; exit 1; }
SEED=$(( SEEDBASE + 10#$IDX * 7919 ))
SEEDP=$(( 10#$IDX * SEEDPMUL + SEEDPOFF ))
TAG=${ARM}_$(printf '%04d' "$((10#$IDX))")
OUT=$OUTDIR/$ARM
mkdir -p "$OUT"
TMP=$(mktemp -d "${TMPDIR:-/tmp}/ampt_${TAG}_XXXX")
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/ana"
cd "$TMP"
ln -sf /dev/null ana/zpc.dat                      # as in production: the large ZPC dump is not kept
ln -sf "$ROOT"/deps/eps09s/eps09sNLO*.dat .
"$AWK" -v n="$NEV" -v sd="$SEED" -v sp="$SEEDP" -v b0="$BMIN" -v b1="$BMAX" -v sys="$SYS" '
  NR==3 && sys=="pp" {print "P               ! PROJ"; next}
  NR==4 && sys=="pp" {print "P               ! TARG"; next}
  NR>=5 && NR<=8 && sys=="pp" {print "1               ! A or Z"; next}
  NR==9  {printf "%-15s ! NEVNT\n", n;  next}
  NR==10 {printf "%-15s ! BMIN\n", b0; next}
  NR==11 {printf "%-15s ! BMAX\n", b1; next}
  NR==28 {printf "%-15s ! ihjsed\n", 0; next}
  NR==29 {printf "%-15s ! HIJING seed\n", sd; next}
  NR==30 {printf "%-15s ! parton cascade seed\n", sp; next}
  {print}' "$ROOT/inputs/$INPUT" > input.ampt
echo "=== $(date '+%F %T') arm=$ARM idx=$IDX sys=$SYS nev=$NEV b=$BMIN-$BMAX input=$INPUT seed=$SEED/$SEEDP config=$CONFIG"
T0=$SECONDS
# shellcheck disable=SC2086
echo 20261017 | env $DMENV LHAPDF_DATA_PATH="$ROOT/deps/pdfsets" OMP_NUM_THREADS=1 "$ROOT/bin/ampt" > ampt.log 2>&1 \
  || { tail -30 ampt.log; echo "ERROR: ampt failed" >&2; exit 1; }
NDONE=$("$AWK" 'NF==11{n++} END{print n+0}' ana/ampt.dat)
[ "$NDONE" -eq "$NEV" ] || { tail -30 ampt.log; echo "ERROR: $NDONE of $NEV events" >&2; exit 1; }
"$AWK" -v seed="$SEED" -f "$ROOT/scripts/reduce_dmeson_v2.awk" ana/ampt.dat | gzip > "$OUT/red_$TAG.dat.gz"
[ -s ana/charm.dat ] && gzip -c ana/charm.dat > "$OUT/charm_$TAG.dat.gz"
grep "^ # hq" ampt.log | gzip > "$OUT/hq_$TAG.log.gz" || true
echo "=== $(date '+%F %T') $TAG done: $NDONE events in $(( SECONDS - T0 )) s"

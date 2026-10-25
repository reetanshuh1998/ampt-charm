#!/usr/bin/env bash
# test_exact_reproduction.sh : regenerate the first events of job 1 of one production arm per configuration used in the
# paper and compare them with the files written on the KEKCC farm (tests/reference/, cut from prod/ARM/*_0001):
#   reference  : dm24_e_a (Au+Au, first 5 events), dm24_epp (p+p, first 300 events)
#   c1 variant : dm25_c1 (Au+Au, 5), dm25_ppc1 (p+p, 300)
#   public AMPT: dm25_pub_a (Au+Au, 5), dm25_pppub (p+p, 300)
#   Duke-Owens : dm23_a (Au+Au, 5; the tilt-law points)
# The seven jobs run in parallel (~3 min with 7 cores, ~6 min on one).
# The production jobs ran 300 (Au+Au) and 30000 (p+p) events; event k depends only on the seeds and events < k, so the
# first events of a shorter run are the same events.
# PASS = charm.dat (every charm quark/hadron record) and the reduced event file identical byte for byte.  If only the
# light-hadron lines of the reduced file differ, your gawk's rand() differs from the farm's (gawk 5.1): the charm
# results are still exactly reproduced, the 25 %/10 % light-hadron prescale is not.
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
REF=$ROOT/tests/reference
OUT=$(mktemp -d "${TMPDIR:-/tmp}/ampt_repro_XXXX")
trap 'rm -rf "$OUT"' EXIT
fail=0
command -v gawk > /dev/null || echo "NOTE  gawk not found: the light-hadron prescale cannot be bit-identical (expect PASS*)"
check() {   # check LABEL ARM REDREF CHARMREF
  local lab=$1 arm=$2
  local red=$OUT/$arm/red_${arm}_0001.dat.gz chm=$OUT/$arm/charm_${arm}_0001.dat.gz
  if cmp -s <(gzip -dc "$chm") "$REF/$4"; then echo "PASS  $lab charm.dat identical ($(wc -l < "$REF/$4") lines)"
  else echo "FAIL  $lab charm.dat differs"; fail=1; fi
  if cmp -s <(gzip -dc "$red") "$REF/$3"; then echo "PASS  $lab reduced event file identical ($(wc -l < "$REF/$3") lines)"
  elif cmp -s <(gzip -dc "$red" | awk 'NF==13 || NF==9 || $1=="Q"') <(awk 'NF==13 || NF==9 || $1=="Q"' "$REF/$3"); then
    echo "PASS* $lab event headers, flow vectors and charm hadrons identical; light-hadron prescale differs (gawk rand)"
  else echo "FAIL  $lab reduced event file differs"; fail=1; fi
}
#    config         arm        sys  bmin bmax nev seed base pmul poff  (docs/CHANGES.pdf, seeds table)
JOBS='
reference_auau dm24_e_a   auau 0 14.5 5   20261107 13 9
reference_pp   dm24_epp   pp   0 2.0  300 20261107 13 9
c1_auau        dm25_c1    auau 0 14.5 5   20261116 31 18
c1_pp          dm25_ppc1  pp   0 2.0  300 20261117 31 19
public_auau    dm25_pub_a auau 0 14.5 5   20261118 31 20
public_pp      dm25_pppub pp   0 2.0  300 20261120 31 22
do_auau        dm23_a     auau 0 14.5 5   20261105 13 7'
while read -r cfg arm sys b0 b1 nev base pm po; do       # all jobs in parallel, one core each
  [ -n "$cfg" ] || continue
  "$ROOT/scripts/run_job.sh" "$cfg" "$arm" "$sys" "$b0" "$b1" "$nev" "$base" "$pm" "$po" 1 "$OUT" > "$OUT/$arm.log" 2>&1 &
done <<< "$JOBS"
wait
while read -r cfg arm sys b0 b1 nev base pm po; do
  [ -n "$cfg" ] || continue
  [ -s "$OUT/$arm/charm_${arm}_0001.dat.gz" ] || { echo "FAIL  $arm: run_job.sh failed"; tail -3 "$OUT/$arm.log"; fail=1; continue; }
  check "$(printf '%-6s %-10s job 1 (%s events)' "$sys" "$arm" "$nev")" "$arm" "red_${arm}_0001_first$nev.dat" \
    "charm_${arm}_0001_first$nev.dat"
done <<< "$JOBS"
[ $fail -eq 0 ] && echo "exact reproduction: OK" || { echo "exact reproduction: FAILED"; exit 1; }

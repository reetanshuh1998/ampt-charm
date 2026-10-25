#!/usr/bin/env bash
# test_exact_reproduction.sh : regenerate the first events of two jobs of the paper production and compare them with the
# files written on the KEKCC farm (tests/reference/, cut from prod/dm24_e_a/*_0001 and prod/dm24_epp/*_0001).
#   Au+Au: arm dm24_e_a, job 1, reference configuration, first 5 events   (~1 min)
#   p+p  : arm dm24_epp, job 1, reference configuration, first 300 events (~15 s)
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
check() {   # check LABEL ARM REDREF CHARMREF
  local lab=$1 arm=$2 red=$OUT/$arm/red_${arm}_0001.dat.gz chm=$OUT/$arm/charm_${arm}_0001.dat.gz
  if cmp -s <(gzip -dc "$chm") "$REF/$4"; then echo "PASS  $lab charm.dat identical ($(wc -l < "$REF/$4") lines)"
  else echo "FAIL  $lab charm.dat differs"; fail=1; fi
  if cmp -s <(gzip -dc "$red") "$REF/$3"; then echo "PASS  $lab reduced event file identical ($(wc -l < "$REF/$3") lines)"
  elif cmp -s <(gzip -dc "$red" | awk 'NF==13 || NF==9 || $1=="Q"') <(awk 'NF==13 || NF==9 || $1=="Q"' "$REF/$3"); then
    echo "PASS* $lab event headers, flow vectors and charm hadrons identical; light-hadron prescale differs (gawk rand)"
  else echo "FAIL  $lab reduced event file differs"; fail=1; fi
}
"$ROOT/scripts/run_job.sh" reference_auau dm24_e_a auau 0 14.5 5   20261107 13 9 1 "$OUT" || fail=1
"$ROOT/scripts/run_job.sh" reference_pp   dm24_epp pp   0 2.0  300 20261107 13 9 1 "$OUT" || fail=1
check "Au+Au dm24_e_a job 1 (5 events)  " dm24_e_a red_dm24_e_a_0001_first5.dat charm_dm24_e_a_0001_first5.dat
check "p+p   dm24_epp job 1 (300 events)" dm24_epp red_dm24_epp_0001_first300.dat charm_dm24_epp_0001_first300.dat
[ $fail -eq 0 ] && echo "exact reproduction: OK" || { echo "exact reproduction: FAILED"; exit 1; }

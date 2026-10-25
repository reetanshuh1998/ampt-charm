#!/usr/bin/env bash
# combine_arms.sh PRODDIR OUT ARM1 ARM2 ... : sample PRODDIR/OUT = symlinks to the red_/charm_/hq_/hqerr_ files of the
# arms, job i of the k-th arm (k = 0, 1, ...) at index 1000 k + i.  This is how the paper samples were combined on the
# farm (README Sec. 2.2): the analysis seeds each file's feed-down decays by its position in the sorted file list, so
# the same arms in the same order give the same sample, file for file.
#   dm25_main   = dm24_e_a dm24_e_b dm25_e_c dm25_e_d dm25_e_e dm25_e_f dm25_e_g dm25_e_h
#   dm25_ppmain = dm24_epp dm25_epp        dm25_pub = dm25_pub_a dm25_pub_b
#   dm24_e_ab   = dm24_e_a dm24_e_b        dm23_ab  = dm23_a dm23_b
set -euo pipefail
P=${1:?PRODDIR} OUT=${2:?OUT}; shift 2
mkdir -p "$P/$OUT"; k=0
for A in "$@"; do
  for f in "$P/$A"/red_"$A"_*.dat.gz; do
    [ -e "$f" ] || { echo "no red_ files in $P/$A" >&2; exit 1; }
    i=${f##*_}; i=${i%.dat.gz}; j=$(printf '%04d' $((10#$i + 1000 * k)))
    for t in red:dat charm:dat hq:log hqerr:txt; do p=${t%%:*}; e=${t##*:}
      if [ -e "$P/$A/${p}_${A}_$i.$e.gz" ]; then ln -sf "../$A/${p}_${A}_$i.$e.gz" "$P/$OUT/${p}_${OUT}_$j.$e.gz"; fi
    done
  done
  k=$((k + 1))
done
echo "combined $(ls "$P/$OUT"/red_* | wc -l) red, $(ls "$P/$OUT"/hq_* 2>/dev/null | wc -l) hq," \
     "$(ls "$P/$OUT"/hqerr_* 2>/dev/null | wc -l) hqerr files into $P/$OUT"

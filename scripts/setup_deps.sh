#!/usr/bin/env bash
# setup_deps.sh : download and verify every external input, build LHAPDF.  Everything goes into deps/.
#   1. official AMPT v1.26t9b/v2.26t9b (Z.-W. Lin), patched by build.sh
#   2. LHAPDF 6.5.1 (the version used for the paper production), built without Python
#   3. CTEQ6.1M grid (LHAPDF set cteq61)
#   4. EPS09s NLO code and tables (Helenius, Eskola, Honkanen, Salgado), used unmodified
# Every archive is checked against the SHA-256 recorded below; the script stops on any mismatch.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DEPS=$ROOT/deps
mkdir -p "$DEPS/dl"
cd "$DEPS/dl"

fetch() {   # fetch NAME URL SHA256
  local name=$1 url=$2 sum=$3
  if [ ! -s "$name" ] || ! echo "$sum  $name" | sha256sum -c --status; then
    echo ">> downloading $name"
    curl -fL --retry 3 -o "$name.part" "$url" && mv "$name.part" "$name"
  fi
  echo "$sum  $name" | sha256sum -c || { echo "ERROR: checksum of $name does not match; remove deps/dl/$name and retry" >&2; exit 1; }
}

fetch ampt-v1.26t9b-v2.26t9b.zip "https://myweb.ecu.edu/linz/ampt/ampt-v1.26t9b-v2.26t9b.zip" \
      9441b5f77c2ab91a57b291abd4afd12de7968f9cbe9f3cc8dbe60fbf5293ed55
fetch LHAPDF-6.5.1.tar.gz "https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.5.1.tar.gz" \
      1256419e2227d1a4f93387fe1da805e648351417d3755e8af5a30a35a6a66751
fetch cteq61.tar.gz "http://lhapdfsets.web.cern.ch/lhapdfsets/current/cteq61.tar.gz" \
      5a6f1393d1ff4d946dd23642628f5ca14a6648ab06bdebcf949a8f550a620b92
fetch eps09s.tgz "https://users.jyu.fi/~kaeskola/EPS09s/eps09s.tgz" \
      294bce5e312a06443bd083c2b681b450b1eb1f2b88a6c0cfc4182e1a85ffdf2c

# ---- LHAPDF 6.5.1 -> deps/lhapdf
if [ ! -x "$DEPS/lhapdf/bin/lhapdf-config" ]; then
  echo ">> building LHAPDF 6.5.1 (a few minutes)"
  rm -rf "$DEPS/build-lhapdf" && mkdir -p "$DEPS/build-lhapdf"
  tar xzf LHAPDF-6.5.1.tar.gz -C "$DEPS/build-lhapdf"
  ( cd "$DEPS/build-lhapdf/LHAPDF-6.5.1" && ./configure --prefix="$DEPS/lhapdf" --disable-python > configure.log 2>&1 \
      && make -j"$(nproc 2>/dev/null || echo 2)" > make.log 2>&1 && make install > install.log 2>&1 ) \
    || { echo "ERROR: LHAPDF build failed, see deps/build-lhapdf/LHAPDF-6.5.1/*.log" >&2; exit 1; }
fi
echo "LHAPDF $("$DEPS/lhapdf/bin/lhapdf-config" --version)"

# ---- CTEQ6.1M grid -> deps/pdfsets/cteq61
mkdir -p "$DEPS/pdfsets"
[ -f "$DEPS/pdfsets/cteq61/cteq61.info" ] || tar xzf cteq61.tar.gz -C "$DEPS/pdfsets"
cp -f "$DEPS/lhapdf/share/LHAPDF/lhapdf.conf" "$DEPS/pdfsets/" 2>/dev/null || true
cp -f "$DEPS/lhapdf/share/LHAPDF/pdfsets.index" "$DEPS/pdfsets/" 2>/dev/null || true

# ---- EPS09s -> deps/eps09s
[ -f "$DEPS/eps09s/eps09s.cpp" ] || tar xzf eps09s.tgz -C "$DEPS"
ls "$DEPS"/eps09s/eps09sNLO*.dat > /dev/null

echo "dependencies ready in $DEPS"

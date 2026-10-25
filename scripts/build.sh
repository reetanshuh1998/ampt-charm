#!/usr/bin/env bash
# build.sh : unpack the official AMPT v1.26t9b/v2.26t9b, apply patches/ampt-v1.26t9b-v2.26t9b_charm.patch, add the new
# files of src_new/ (charm_mod.f, coal_new.f, hqctl.inc, eps09s_f.cpp, Makefile) and compile -> bin/ampt.
# Needs scripts/setup_deps.sh first.  The patched tree stays in build/ampt-src/ for inspection.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DEPS=$ROOT/deps
SRC=$ROOT/build/ampt-src
[ -x "$DEPS/lhapdf/bin/lhapdf-config" ] || { echo "run scripts/setup_deps.sh first" >&2; exit 1; }
rm -rf "$ROOT/build" && mkdir -p "$ROOT/build"
( cd "$ROOT/build" && unzip -q "$DEPS/dl/ampt-v1.26t9b-v2.26t9b.zip" && mv Ampt-v1.26t9b-v2.26t9b ampt-src )
# the patch was made with "diff -ruN a b" on the seven stock files
( cd "$SRC" && patch -p1 --no-backup-if-mismatch < "$ROOT/patches/ampt-v1.26t9b-v2.26t9b_charm.patch" )
cp "$ROOT"/src_new/* "$SRC"/
( cd "$SRC" && make LHAPDF_PREFIX="$DEPS/lhapdf" EPS09S_DIR="$DEPS/eps09s" > make.log 2>&1 ) \
  || { tail -30 "$SRC/make.log"; echo "ERROR: compilation failed (build/ampt-src/make.log)" >&2; exit 1; }
mkdir -p "$ROOT/bin" && cp "$SRC/ampt" "$ROOT/bin/ampt"
echo "built bin/ampt  md5 $(md5sum "$ROOT/bin/ampt" | cut -c1-8)  ($(gfortran --version | head -1))"

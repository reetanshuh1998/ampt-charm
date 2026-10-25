#!/usr/bin/env bash
# check_prereqs.sh : report every missing tool before anything is downloaded or compiled.
miss=0
for t in gfortran g++ make curl unzip patch tar gzip sha256sum md5sum python3; do
  command -v "$t" > /dev/null || { echo "MISSING: $t"; miss=1; }
done
command -v gawk > /dev/null || echo "note: gawk not found; any awk works, but only gawk reproduces the farm's light-hadron prescale bit for bit"
python3 -c "import numpy, matplotlib" 2> /dev/null || { echo "MISSING: python3 packages numpy and matplotlib (pip install -r requirements.txt)"; miss=1; }
if [ $miss -ne 0 ]; then
  echo "Install the missing tools first. Ubuntu/Debian:"
  echo "  sudo apt install build-essential gfortran g++ make gawk curl unzip patch python3-numpy python3-matplotlib"
  echo "Fedora/RHEL: sudo dnf install gcc-gfortran gcc-c++ make gawk curl unzip patch python3-numpy python3-matplotlib"
  echo "macOS (Homebrew): brew install gcc gawk coreutils && pip3 install numpy matplotlib   (untested)"
  echo "Or use the Docker image (README.md)."
  exit 1
fi
echo "prerequisites OK"

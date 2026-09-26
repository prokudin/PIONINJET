#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
: "${FC:=gfortran}"
: "${CXX:=g++}"
: "${LHAPDF_PREFIX:?Set LHAPDF_PREFIX to an existing LHAPDF 6 installation}"
test -f build/kernels.inc || { printf '%s\n' 'Run python3 prepare.py first.' >&2; exit 1; }
make -C build/hoppet -j 4 FC="$FC" F90="$FC" \
    FFLAGS='-O2 -fPIC -std=legacy -fallow-argument-mismatch' \
    F90FLAGS='-O2 -fPIC -fallow-argument-mismatch'
"$FC" -O2 -fPIC -std=legacy -c build/fDSS.f -o build/fDSS.o
"$FC" -O2 -fPIC -std=legacy -c build/DSSV-2008.f -o build/DSSV.o
"$CXX" -O2 -std=c++17 -fPIC -shared bridge.cpp validation/reference.cpp -Ibuild \
    -I"$LHAPDF_PREFIX/include" -L"$LHAPDF_PREFIX/lib" \
    -Wl,-rpath,"$LHAPDF_PREFIX/lib" -o build/libpionjet.so \
    build/fDSS.o build/DSSV.o build/hoppet/libhoppet_v1.a -lLHAPDF -lgfortran

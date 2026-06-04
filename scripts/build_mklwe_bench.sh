#!/usr/bin/env bash
# Configure + build the MKFHE Alg.2 (LWE) noise/timing benchmark
# (src/binfhe/examples/boolean-mklwe-bench.cpp).
#
# Mirrors the README build flags (NATIVE_SIZE=32, WITH_NTL=ON, WITH_NATIVEOPT=ON,
# WITH_OPENMP=OFF). The README pins clang-12; this script defaults to whatever
# cc/c++ are on PATH (g++ works fine). Set CC/CXX to override, e.g.
#   CC=clang-12 CXX=clang++-12 scripts/build_mklwe_bench.sh
#
# Build dependencies: NTL + GMP development headers, e.g. on Debian/Ubuntu:
#   sudo apt-get install -y libntl-dev libgmp-dev
#
# Note: PreLoad.cmake forces the "Unix Makefiles" generator, so do NOT pass -G.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT/build}"

CMAKE_ARGS=(
    -DNATIVE_SIZE=32
    -DWITH_NTL=ON
    -DWITH_NATIVEOPT=ON
    -DWITH_OPENMP=OFF
    -DBUILD_UNITTESTS=OFF
    -DBUILD_BENCHMARKS=OFF
    -DBUILD_EXAMPLES=ON
    -DCMAKE_BUILD_TYPE=Release
)
[[ -n "${CC:-}"  ]] && CMAKE_ARGS+=(-DCMAKE_C_COMPILER="$CC")
[[ -n "${CXX:-}" ]] && CMAKE_ARGS+=(-DCMAKE_CXX_COMPILER="$CXX")

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
cmake "$ROOT" "${CMAKE_ARGS[@]}"
make -j"$(nproc)" boolean-mklwe-bench

echo
echo "Built: $BUILD_DIR/bin/examples/binfhe/boolean-mklwe-bench"

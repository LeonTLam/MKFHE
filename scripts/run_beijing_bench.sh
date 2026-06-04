#!/usr/bin/env bash
# Run the MKFHE Alg.2 (LWE) noise benchmarks for every paper (k, l) config and
# rename the outputs to  beijing_k{k}_l{l}_{growth,noise}.csv.
#
#   k=2  -> l=3 (Bg=512)
#   k=4  -> l=3 (Bg=512)
#   k=8  -> l=3 (Bg=512)
#   k=16 -> l=4 (Bg=128)     # paper's lower-base row
#
# For each config it runs:
#   boolean-mklwe-bench k --noise-growth --reps R   -> beijing_k{k}_l{l}_growth.csv
#   boolean-mklwe-bench k --noise        --reps R   -> beijing_k{k}_l{l}_noise.csv
#
# Then (unless --no-pfail) runs scripts/p_fail_calc.py over the *_growth.csv
# files, emitting a beijing_k{k}_l{l}_growth_pfail.csv next to each.
#
# Usage:
#   scripts/run_beijing_bench.sh [reps] [outdir]
# Examples:
#   scripts/run_beijing_bench.sh                 # 100 reps -> ./beijing_results
#   scripts/run_beijing_bench.sh 100 results     # 100 reps -> ./results
#   REPS_ONLY_K=2 scripts/run_beijing_bench.sh   # just k=2 (quick smoke)
#
# Env overrides:
#   BENCH=/path/to/boolean-mklwe-bench   (default: build/bin/examples/binfhe/...)
#   REPS_ONLY_K=2|4|8|16                 restrict to a single k
#   NO_PFAIL=1                           skip the p_fail post-processing
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REPS="${1:-100}"
OUTDIR="${2:-beijing_results}"
BENCH="${BENCH:-$ROOT/build/bin/examples/binfhe/boolean-mklwe-bench}"
PFAIL="$ROOT/scripts/p_fail_calc.py"

if [[ ! -x "$BENCH" ]]; then
    echo "ERROR: bench binary not found/executable: $BENCH" >&2
    echo "       build it first:  scripts/build_mklwe_bench.sh" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"   # absolutise (the bench writes CSVs in CWD)

ks=(2 4 8 16)
ls_=(3 3 3 4)

for idx in "${!ks[@]}"; do
    k="${ks[$idx]}"
    l="${ls_[$idx]}"
    if [[ -n "${REPS_ONLY_K:-}" && "$k" != "$REPS_ONLY_K" ]]; then
        continue
    fi

    echo "=================================================================="
    echo " k=$k  l=$l   --noise-growth  (reps=$REPS)"
    echo "=================================================================="
    # --csv / --keygen-csv make the bench write the final filenames directly:
    # no fixed-name + rename, and concurrent runs never clobber each other.
    "$BENCH" "$k" --noise-growth --reps "$REPS" \
        --csv        "$OUTDIR/beijing_k${k}_l${l}_growth.csv" \
        --keygen-csv "$OUTDIR/beijing_k${k}_l${l}_growth_keygen.csv"

    echo "=================================================================="
    echo " k=$k  l=$l   --noise          (reps=$REPS)"
    echo "=================================================================="
    "$BENCH" "$k" --noise --reps "$REPS" \
        --csv        "$OUTDIR/beijing_k${k}_l${l}_noise.csv" \
        --keygen-csv "$OUTDIR/beijing_k${k}_l${l}_noise_keygen.csv"
done

echo
echo "Bench CSVs written to $OUTDIR/:"
ls -1 "$OUTDIR"/beijing_*.csv

if [[ -z "${NO_PFAIL:-}" ]]; then
    echo
    echo "=================================================================="
    echo " p_fail post-processing (growth CSVs -> *_pfail.csv)"
    echo "=================================================================="
    # Only the *_growth.csv files carry the 'bootstrap' row p_fail_calc needs;
    # the *_noise.csv component files already include a per-stage p_fail column.
    python3 "$PFAIL" --emit-csv "$OUTDIR"/beijing_*_growth.csv
fi

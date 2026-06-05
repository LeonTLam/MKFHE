#!/usr/bin/env bash
# Run the MKFHE Alg.2 (LWE) noise benchmarks for every paper (k, l) config and
# rename the outputs to  beijing_k{k}_l{l}_{growth,noise}.csv.
#
#   k=2  -> l=3 (Bg=512)
#   k=4  -> l=3 (Bg=512)
#   k=8  -> l=3 (Bg=512)
#   k=16 -> l=4 (Bg=128)     # paper's lower-base row
#
# For each config it runs ONE invocation (one keygen shared by both benches):
#   boolean-mklwe-bench k --noise-all --reps R
#     --growth-csv beijing_k{k}_l{l}_growth.csv   (1->256 add trajectory)
#     --noise-csv  beijing_k{k}_l{l}_noise.csv    (per-stage component noise)
#     --keygen-csv beijing_k{k}_l{l}_keygen.csv
#
# Then (unless NO_PFAIL=1) runs scripts/p_fail_calc.py over the *_growth.csv
# files, emitting a beijing_k{k}_l{l}_growth_pfail.csv next to each.
#
# MEMORY: MKBTKeyGen peaks roughly proportional to k. k=16 (l=4) has been
# observed to need ~16 GB RSS; on smaller hosts it will OOM. A failing k is
# reported and skipped -- the remaining configs and the p_fail step still run.
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
# Note: -e is intentionally OFF so one k's failure (e.g. OOM) does not abort
# the whole sweep.
set -uo pipefail

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
    echo " k=$k  l=$l   --noise-all  (growth + components, reps=$REPS)"
    echo "=================================================================="
    # ONE invocation per k => ONE keygen (MKBTKeyGen) shared by both noise
    # benches. --growth-csv/--noise-csv/--keygen-csv write final filenames
    # directly, so concurrent runs never clobber each other.
    if "$BENCH" "$k" --noise-all --reps "$REPS" \
        --growth-csv "$OUTDIR/beijing_k${k}_l${l}_growth.csv" \
        --noise-csv  "$OUTDIR/beijing_k${k}_l${l}_noise.csv" \
        --keygen-csv "$OUTDIR/beijing_k${k}_l${l}_keygen.csv"; then
        echo "  [ok] k=$k l=$l complete"
    else
        rc=$?
        echo "  [WARN] k=$k l=$l FAILED (exit $rc) -- likely OOM during MKBTKeyGen" >&2
        echo "         (k=16 needs ~16 GB). Skipping; continuing with remaining k." >&2
        # A killed process leaves 0-byte CSVs (header never flushed); drop them.
        for suf in growth noise keygen; do
            f="$OUTDIR/beijing_k${k}_l${l}_${suf}.csv"
            [[ -s "$f" ]] || rm -f "$f"
        done
    fi
done

echo
echo "Bench CSVs written to $OUTDIR/:"
ls -1 "$OUTDIR"/beijing_*.csv 2>/dev/null || echo "  (none)"

if [[ -z "${NO_PFAIL:-}" ]]; then
    echo
    echo "=================================================================="
    echo " p_fail post-processing (growth CSVs -> *_pfail.csv)"
    echo "=================================================================="
    # Only the *_growth.csv files carry the 'bootstrap' row p_fail_calc needs;
    # the *_noise.csv component files already include a per-stage p_fail column.
    python3 "$PFAIL" --emit-csv "$OUTDIR"/beijing_*_growth.csv
fi

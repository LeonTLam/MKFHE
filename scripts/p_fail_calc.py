#!/usr/bin/env python3
"""
p_fail calculator for MKFHE / FINALLY noise_growth CSVs.

Reads one or more `*_growth.csv` (a.k.a. `*_noise.csv` in FINALLY) files
produced by `boolean-mklwe-bench --noise-growth` and estimates the bootstrap
decryption failure probability under three increasingly conservative noise
models. Reports whether each (k, l, B) configuration meets common target
thresholds:
  - 2^-200 : the FINALLY paper's bound at small-B configurations
  - 2^-128 : a common cryptographic-grade threshold
  - 2^-40  : the standard practical p_fail target for FHE schemes

Usage:
    python3 scripts/p_fail_calc.py results/beijing_k2_l3_growth.csv
    python3 scripts/p_fail_calc.py results/beijing_k*_growth.csv      # batch
    python3 scripts/p_fail_calc.py --emit-csv results/*_growth.csv    # + per-file CSV
    python3 scripts/p_fail_calc.py --emit-csv --out-dir pfail results/*_growth.csv

Options:
    --emit-csv         For each input file write <stem>_pfail.csv with the parsed
                       bootstrap distribution + every p_fail estimate + target
                       PASS/FAIL flags (one data row per input).
    --out-dir DIR      Directory for the emitted CSVs (default: next to each
                       input). Implies --emit-csv.
    --quiet            Suppress the per-file text report (useful with --emit-csv).

The CSV must contain a 'bootstrap' row with columns:
    bits_mean, bits_max, bits_stddev, budget_bits

(All other rows are ignored. Files without a 'bootstrap' row -- e.g. the
 per-stage `*_noise.csv` component files, which already carry their own p_fail
 column -- are skipped.)
"""

from __future__ import annotations
import csv
import math
import sys
from pathlib import Path
from dataclasses import dataclass

# Target p_fail thresholds we check against (in log2 units).
TARGETS = [
    ("2^-200 (paper)",  -200),
    ("2^-128 (crypto)", -128),
    ("2^-40 (practical)", -40),
]

# E[log_2 |Z|] for Z ~ N(0,1).  Used to recover linear sigma from
# bits_mean = E[log_2 |noise|].  Derivation:
#   E[log_2 |X|] = log_2(sigma) + E[log_2 |Z|]
#   E[log_2 |Z|] = (ln 2)^-1 * (-gamma_euler/2 - ln(2)/2 + ln(sqrt(2/pi)))
# Empirically this evaluates to ~ -0.6351 in nats, / ln(2) = -0.9163 bits.
LOG2_ABS_GAUSSIAN_MEAN = -0.9163


@dataclass
class BootstrapRow:
    k: int
    N: int
    n: int
    B: int
    l: int
    log2Q: int
    log2q: int
    reps: int
    bits_min: float
    bits_mean: float
    bits_max: float
    bits_stddev: float
    budget_bits: float


def parse_bootstrap_row(path: Path) -> BootstrapRow | None:
    with path.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row.get("step", "").strip().strip('"') != "bootstrap":
                continue
            return BootstrapRow(
                k=int(row["k"]),
                N=int(row["N"]),
                n=int(row["n"]),
                B=int(row["B"]),
                l=int(row["l"]),
                log2Q=int(row["log2Q"]),
                log2q=int(row["log2q"]),
                reps=int(row["reps"]),
                bits_min=float(row["bits_min"]),
                bits_mean=float(row["bits_mean"]),
                bits_max=float(row["bits_max"]),
                bits_stddev=float(row["bits_stddev"]),
                budget_bits=float(row["budget_bits"]),
            )
    return None


def gaussian_pfail_log2(z: float) -> float:
    """log_2 P(|N(0,1)| > z) using a numerically stable approximation.

    For z >= 0:
        P(|N(0,1)| > z) = 2 * Phi(-z) = erfc(z / sqrt(2))
    For large z (where the erfc is below 1e-300), use the asymptotic:
        log P ~ -z^2/2 - 0.5*log(z) + log(2/sqrt(2*pi))
    """
    if z < 0:
        z = 0.0
    if z < 30.0:
        # Direct computation: use erfc which handles up to z ~ 26-27.
        p = math.erfc(z / math.sqrt(2))
        if p > 0.0:
            return math.log2(p)
    # Asymptotic regime (Mills ratio).  Two-sided.
    log_p_nat = (-z * z / 2.0
                 - 0.5 * math.log(z * z + 1.0)
                 + math.log(2.0 / math.sqrt(2.0 * math.pi)))
    return log_p_nat / math.log(2.0)


def estimate_pfail(row: BootstrapRow) -> dict[str, float]:
    """Return three p_fail estimates (log_2 units), each more conservative
    than the previous one.

      1. Empirical log-normal:
            assumes log_2 |noise| ~ N(bits_mean, bits_stddev^2).
            P_fail = P(log_2 |noise| > budget_bits).
         This is the WEAKEST bound -- it treats the log of the noise as
         Gaussian, which has a heavier tail than the true noise distribution.

      2. Gaussian on linear scale:
            recovers linear sigma from bits_mean
              (log_2 sigma = bits_mean + 0.9163)
            then computes P(|N(0, sigma)| > T) where T = 2^(budget_bits + 3).
            Wait -- budget_bits is log_2(q/8), so T = q/8 = 2^budget_bits * 8
            ... no, budget_bits is already log_2(q/8) so T = 2^budget_bits.
            For decryption budget q/4 (binary plaintext, Delta = q/4):
              T_dec = 2 * 2^budget_bits = 2^(budget_bits+1)
            Some implementations use q/8 with a safety factor.
         This is the STANDARD Gaussian bound for FHE correctness analysis.

      3. Subgaussian tightness:
            same as (2) but treats sigma as the subgaussian parameter gamma
            (gamma <= sigma_Gaussian for any distribution).  Identical
            formula, just relabelled, since the tail bound is the same:
              P(|X| > t) <= 2 exp(-t^2 / (2 gamma^2))
            We report this for clarity in the paper's framework.
    """
    out: dict[str, float] = {}

    # 1. Empirical log-normal extrapolation.
    if row.bits_stddev > 0:
        z = (row.budget_bits - row.bits_mean) / row.bits_stddev
        out["empirical_log_normal"] = gaussian_pfail_log2(z)
    else:
        out["empirical_log_normal"] = float("-inf")

    # 2. Gaussian on linear scale, decryption threshold = q/4.
    # log_2 sigma_linear = bits_mean - E[log_2 |Z|]   (subtract because we added)
    # but bits_mean = log_2(sigma) + E[log_2 |Z|] = log_2(sigma) - 0.9163
    # so log_2(sigma) = bits_mean + 0.9163.
    log2_sigma = row.bits_mean - LOG2_ABS_GAUSSIAN_MEAN
    # Decryption budget: q/4 (the Delta encoding).  budget_bits in CSV is log_2(q/8).
    log2_T = row.budget_bits + 1.0   # q/4 = 2 * q/8
    z = 2.0 ** (log2_T - log2_sigma)
    out["gaussian_q_over_4"] = gaussian_pfail_log2(z)

    # 2b. Same Gaussian but with budget = q/8 (conservative).
    z_q8 = 2.0 ** (row.budget_bits - log2_sigma)
    out["gaussian_q_over_8"] = gaussian_pfail_log2(z_q8)

    # 3. Subgaussian tightness bound, same formula as Gaussian but
    # explicitly named.
    out["subgaussian_q_over_4"] = out["gaussian_q_over_4"]

    return out


def check_targets(log2_pfail: float) -> list[tuple[str, bool]]:
    """For each target threshold, return whether the estimate meets it."""
    return [(name, log2_pfail <= log2_target) for name, log2_target in TARGETS]


def print_report(path: Path, row: BootstrapRow) -> None:
    estimates = estimate_pfail(row)

    print(f"\n{'='*78}")
    print(f"  {path.name}")
    print(f"  k={row.k}  N={row.N}  n={row.n}  B={row.B}  l={row.l}")
    print(f"  log2Q={row.log2Q}  log2q={row.log2q}  reps={row.reps}")
    print(f"{'='*78}")

    print(f"\n  Bootstrap noise distribution:")
    print(f"    bits_mean   = {row.bits_mean:7.3f}")
    print(f"    bits_max    = {row.bits_max:7.3f}")
    print(f"    bits_stddev = {row.bits_stddev:7.3f}")
    print(f"    budget_bits = {row.budget_bits:7.3f}  (= log_2(q/8))")
    headroom_mean = row.budget_bits - row.bits_mean
    headroom_max  = row.budget_bits - row.bits_max
    print(f"    headroom    = {headroom_mean:5.2f} bits (mean), "
          f"{headroom_max:5.2f} bits (max over {row.reps} reps)")

    # Theoretical bound (standard FHE analysis).  This is what the paper
    # uses; it assumes Gaussian/subgaussian noise on the LINEAR scale,
    # recovers sigma from bits_mean, and computes the tail probability
    # against the decryption threshold.  This is the bound to compare
    # against the paper's 2^-200 claim.
    # Decode boundary for a boolean message at Delta = q/4 is the MIDPOINT q/8
    # (|e| > q/8 flips the decoded bit), matching FINALLY. This q/8 tail is the
    # correct decryption-failure probability and the one we headline.
    pfail_theoretical = estimates["gaussian_q_over_8"]

    # Empirical "consistency check".  Treats log2|noise| as normal with
    # the measured (mean, stddev).  ALWAYS more pessimistic than the
    # theoretical bound for actually-Gaussian noise, because the log
    # of a half-Gaussian has a thinner upper tail than a log-normal.
    # Useful only as a model-free sanity check ("even under this
    # heavier-tailed model, p_fail is still X").
    pfail_empirical = estimates["empirical_log_normal"]

    # Empirical observable floor: with N reps and zero failures observed,
    # we can only state p_fail <= 1/N at ~63% confidence (one-sided).
    # For a 95% upper bound, p_fail <= 3/N (rule of three).
    empirical_floor_p63 = -math.log2(row.reps) if row.reps > 0 else 0.0
    empirical_floor_p95 = -math.log2(row.reps / 3.0) if row.reps >= 3 else 0.0

    print(f"\n  Theoretical p_fail (standard FHE bound, Gaussian on linear sigma):")
    print(f"    budget = q/8 (decode boundary):  log2(p_fail) = " + _fmt_log2(pfail_theoretical))
    print(f"    budget = q/4 (codeword spacing): log2(p_fail) = " + _fmt_log2(estimates['gaussian_q_over_4']))
    print(f"    -- q/8 is the correct decode boundary (FINALLY); compare to 2^-200 claim.")

    print(f"\n  Empirical consistency check (log-normal extrapolation):")
    print(f"    log2(p_fail) = " + _fmt_log2(pfail_empirical))
    print(f"    -- More pessimistic than the theoretical bound.  PASSING under this")
    print(f"       model means even a heavier-tailed assumption clears the threshold.")
    print(f"       FAILING this does NOT contradict the theoretical bound; it just")
    print(f"       reflects the conservative model.")

    print(f"\n  Empirical observable floor (cannot bound p_fail below this from {row.reps} reps):")
    print(f"    63% confidence:  log2(p_fail) >= {empirical_floor_p63:+7.2f}")
    print(f"    95% confidence:  log2(p_fail) >= {empirical_floor_p95:+7.2f}")
    print(f"    -- The paper's 2^-200 claim is not directly empirically verifiable")
    print(f"       at this rep count; it requires the theoretical analysis above.")

    print(f"\n  Target thresholds (against THEORETICAL bound, Gaussian@q/8):")
    for name, meets in check_targets(pfail_theoretical):
        mark = "PASS" if meets else "FAIL"
        print(f"    [{mark}]  {name}")

    print(f"\n  Target thresholds (against EMPIRICAL log-normal, conservative):")
    for name, meets in check_targets(pfail_empirical):
        mark = "PASS" if meets else "fail*"
        suffix = "" if meets else "  (* heavier-tailed model; theoretical bound may still pass)"
        print(f"    [{mark}]  {name}{suffix}")


def _fmt_log2(val: float) -> str:
    if val < -1e6:
        return f"2^{val:>10.2e}   (essentially zero)"
    if val < -300:
        return f"2^{val:>10.1f}"
    return f"2^{val:>+10.2f}"


# ── CSV emission ──────────────────────────────────────────────────────────────
# Header of the per-input p_fail CSV. One data row per input file: the parsed
# bootstrap distribution, all four p_fail estimates (log2 units), the empirical
# observable floors, and PASS/FAIL flags for each target threshold (evaluated
# against the theoretical Gaussian@q/8 bound (q/8 = decode boundary, FINALLY).
PFAIL_CSV_COLUMNS = [
    "file", "k", "N", "n", "B", "l", "log2Q", "log2q", "reps",
    "bits_min", "bits_mean", "bits_max", "bits_stddev", "budget_bits",
    "headroom_mean_bits", "headroom_max_bits",
    "log2_pfail_gaussian_q4", "log2_pfail_gaussian_q8",
    "log2_pfail_subgaussian_q4", "log2_pfail_empirical_log_normal",
    "log2_floor_p63", "log2_floor_p95",
    "meets_2^-200", "meets_2^-128", "meets_2^-40",
]


def write_pfail_csv(in_path: Path, row: BootstrapRow, out_dir: Path | None) -> Path:
    est = estimate_pfail(row)
    theo = est["gaussian_q_over_8"]   # q/8 = decode boundary (FINALLY)
    meets = {name: ("PASS" if ok else "FAIL") for name, ok in check_targets(theo)}
    floor_p63 = -math.log2(row.reps) if row.reps > 0 else 0.0
    floor_p95 = -math.log2(row.reps / 3.0) if row.reps >= 3 else 0.0

    target_dir = out_dir if out_dir is not None else in_path.parent
    target_dir.mkdir(parents=True, exist_ok=True)
    out_path = target_dir / f"{in_path.stem}_pfail.csv"

    record = {
        "file": in_path.name,
        "k": row.k, "N": row.N, "n": row.n, "B": row.B, "l": row.l,
        "log2Q": row.log2Q, "log2q": row.log2q, "reps": row.reps,
        "bits_min": f"{row.bits_min:.6f}",
        "bits_mean": f"{row.bits_mean:.6f}",
        "bits_max": f"{row.bits_max:.6f}",
        "bits_stddev": f"{row.bits_stddev:.6f}",
        "budget_bits": f"{row.budget_bits:.6f}",
        "headroom_mean_bits": f"{row.budget_bits - row.bits_mean:.6f}",
        "headroom_max_bits": f"{row.budget_bits - row.bits_max:.6f}",
        "log2_pfail_gaussian_q4": f"{theo:.6f}",
        "log2_pfail_gaussian_q8": f"{est['gaussian_q_over_8']:.6f}",
        "log2_pfail_subgaussian_q4": f"{est['subgaussian_q_over_4']:.6f}",
        "log2_pfail_empirical_log_normal": f"{est['empirical_log_normal']:.6f}",
        "log2_floor_p63": f"{floor_p63:.6f}",
        "log2_floor_p95": f"{floor_p95:.6f}",
        "meets_2^-200": meets["2^-200 (paper)"],
        "meets_2^-128": meets["2^-128 (crypto)"],
        "meets_2^-40": meets["2^-40 (practical)"],
    }
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=PFAIL_CSV_COLUMNS)
        writer.writeheader()
        writer.writerow(record)
    return out_path


def main(argv: list[str]) -> int:
    # ── manual flag parsing (paths are the non-flag positional args) ──────────
    emit_csv = False
    quiet = False
    out_dir: Path | None = None
    paths: list[Path] = []
    i = 1
    while i < len(argv):
        a = argv[i]
        if a in ("-h", "--help"):
            print(__doc__)
            return 0
        elif a == "--emit-csv":
            emit_csv = True
        elif a == "--quiet":
            quiet = True
        elif a == "--out-dir":
            if i + 1 >= len(argv):
                print("--out-dir requires a directory argument", file=sys.stderr)
                return 2
            out_dir = Path(argv[i + 1]); emit_csv = True; i += 1
        elif a.startswith("-"):
            print(f"unknown option: {a}", file=sys.stderr)
            return 2
        else:
            paths.append(Path(a))
        i += 1

    if not paths:
        print("usage: p_fail_calc.py [--emit-csv] [--out-dir DIR] [--quiet] "
              "<growth.csv> [more.csv ...]", file=sys.stderr)
        return 2

    missing = [p for p in paths if not p.is_file()]
    if missing:
        print(f"missing files: {missing}", file=sys.stderr)
        return 2

    rows: list[tuple[Path, BootstrapRow]] = []
    for path in paths:
        row = parse_bootstrap_row(path)
        if row is None:
            print(f"  [SKIP]  {path.name} -- no 'bootstrap' row found "
                  f"(component/--noise CSVs already carry per-stage p_fail)",
                  file=sys.stderr)
            continue
        rows.append((path, row))

    written: list[Path] = []
    for path, row in rows:
        if not quiet:
            print_report(path, row)
        if emit_csv:
            out_path = write_pfail_csv(path, row, out_dir)
            written.append(out_path)
            print(f"  [csv]   wrote {out_path}")

    if len(rows) > 1:
        print(f"\n{'='*78}")
        print("  SUMMARY (Gaussian@q/8 decode boundary, log_2 units)")
        print(f"{'='*78}")
        print(f"  {'file':40s}  {'(k,l,B)':18s}  {'log2(p_fail)':>14s}  paper")
        for path, row in rows:
            pfail = estimate_pfail(row)["gaussian_q_over_8"]
            klb = f"({row.k},{row.l},B={row.B})"
            meets_paper = "PASS" if pfail <= -200 else "FAIL"
            pfail_str = (f"2^{pfail:8.2e}" if abs(pfail) > 1e6
                         else f"2^{pfail:+10.2f}")
            print(f"  {path.name:40s}  {klb:18s}  {pfail_str:>14s}  {meets_paper}")

    if written:
        print(f"\n[p_fail] wrote {len(written)} CSV(s):")
        for w in written:
            print(f"  {w}")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

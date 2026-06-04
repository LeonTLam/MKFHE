// boolean-mklwe-bench.cpp
//
// Protocol + timing + noise + keygen-sizing benchmark for the MKFHE
// (SKLC-FHE/MKFHE) NTRU-based MK-FHE, Algorithm 2 (LWE first layer + NTRU
// second layer).
//
// Ports FINALLY's benchmarks/bench_protocol.cpp measurement harness onto the
// MKFHE public API (binfhecontext.h):
//
//   Default mode  (no flag)   — protocol scenario: Setup (timed + sized),
//                               Encryption, NAND+Bootstrap (with 4-input
//                               correctness), level-2 gate, Decryption,
//                               projected summary.
//                               CSV -> bench_protocol_mklwe.csv
//
//   --timing                  — per-op timing: EvalBinGate (full) plus four
//                               instrumented per-component rows (CMux+ExtProduct,
//                               ExternalProduct, KS_NGS_to_LWE, ModSwitch) read
//                               from the library timing profile.
//                               CSV -> bench_timing_mklwe.csv
//
//   --components / --noise    — per-stage LWE noise + erfc decode-error rate for
//                               the public-API-reachable stages (fresh Encrypt,
//                               homomorphic Add, NAND+Bootstrap output). Noise is
//                               computed from the centred phase via public
//                               MKLWECiphertext / MKLWEPrivateKey accessors —
//                               NO library patch required.
//                               CSV -> bench_components_mklwe.csv
//
//   --growth / --noise-growth — LWE noise-budget trajectory: grow noise via
//                               repeated homomorphic Add of Enc(0), then reset
//                               via EvalBinGate(NAND, acc, Enc(1)). Noise in
//                               log2 bits, again from public accessors only.
//                               CSV -> bench_growth_mklwe.csv
//
//   ALL modes ALSO write -> bench_keygen_mklwe.csv (one row per k with
//                            per-primitive timing + per-key material size).
//
// Parameters (128-bit, Alg.2 / LWE), from Table 3:
//   k=2  -> STD128_MKNTRU_LWE     (n=635, q=32749, N=2048, Q~2^27, Bg=512, d=3)
//   k=4  -> STD128_MKNTRU_LWE_2
//   k=8  -> STD128_MKNTRU_LWE_3
//   k=16 -> STD128_MKNTRU_LWE_4   (Bg=128, d=4 — paper's lower-base row)
//
// Key-material size formulas (NATIVE_SIZE=32 → NativeInteger is 4 bytes per
// coefficient):
//
//   LWE SK per party : n bytes  (binary secret stored as one byte per bit)
//   BSK per party    : 2 * (n+1) * (d-1) * N * 4 bytes
//                      [UniEncACCKeyImpl shape [k][1][n+1]; each entry is
//                       UniEncEvalKey(digitsG, 2) with digitsG = d - 1]
//   KSK per party    : 2 * N * (Bks-1) * dks * 4 bytes
//                      [light KSK packs the (Bks-1)*dks scalar LWE keys into
//                       that many RLWE ciphertexts (b, a) of N coeffs each]
//
// --timing per-component coverage depends on the library timing instrumentation
// (binfhe-timing.{h,cpp}, plus ScopedTimer brackets in mk-acc-xzw_B.cpp and the
// MKNTRU_LWE EvalBinGate branch of binfhe-base-scheme.cpp). When the library is
// built without it the four component rows simply carry zero counts and are
// skipped; the "EvalBinGate (full)" row is always emitted.

#include "binfhecontext.h"
#include "binfhe-timing.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace lbcrypto;
using Clock = std::chrono::steady_clock;

static constexpr int N_GATES = 1000;
// MKFHE built with NATIVE_SIZE=32 → NativeInteger is 4 bytes wide.
static constexpr std::size_t BYTES_PER_COEFF = 4;

// ── shared types ──────────────────────────────────────────────────────────────

struct OpResult {
    std::string phase;
    std::string actor;
    std::string op;
    double      ms;
};

struct TimingSummary {
    double min = 0, median = 0, p90 = 0, max = 0, mean = 0, sd = 0;
    std::size_t n = 0;
};
static TimingSummary summarize_ms(std::vector<double> v) {
    TimingSummary s;
    if (v.empty()) return s;
    std::sort(v.begin(), v.end());
    s.n = v.size();
    s.min = v.front();
    s.median = v[v.size() / 2];
    s.p90 = v[std::min(v.size() - 1, (v.size() * 9) / 10)];
    s.max = v.back();
    for (double x : v) s.mean += x;
    s.mean /= static_cast<double>(v.size());
    for (double x : v) s.sd += (x - s.mean) * (x - s.mean);
    s.sd = std::sqrt(s.sd / static_cast<double>(v.size()));
    return s;
}

// ── progress harness ──────────────────────────────────────────────────────────

static int g_step = 0, g_total = 0;
static void set_total(int n) { g_total = n; g_step = 0; }

template <typename Fn>
static OpResult run_op(const std::string& phase, const std::string& actor,
                       const std::string& name, Fn&& fn) {
    ++g_step;
    const int fw = static_cast<int>(std::to_string(g_total).size());
    std::cout << " [" << std::setw(fw) << g_step << "/" << g_total << "] "
              << std::left << std::setw(7) << ("[" + actor + "]") << "  "
              << std::left << std::setw(50) << ("Running: " + name) << " ..."
              << std::flush;
    auto t0 = Clock::now();
    fn();
    double ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    std::cout << "\r [" << std::setw(fw) << g_step << "/" << g_total << "] "
              << std::left << std::setw(7) << ("[" + actor + "]") << "  "
              << std::left << std::setw(6) << "done" << "  "
              << std::left << std::setw(50) << name << std::right << std::fixed
              << std::setprecision(3) << std::setw(10) << ms << " ms\n"
              << std::flush;
    return OpResult{phase, actor, name, ms};
}

static void print_phase_header(const std::string& title) {
    std::cout << "\n=== " << title << " ===\n";
}

static void print_phase_table(const std::string& phase,
                              const std::vector<OpResult>& results) {
    const int WA = 7, WO = 50, WN = 12;
    double total = 0; bool any = false;
    auto rule = [&]() {
        std::cout << "  +" << std::string(WA + 2, '-') << "+"
                  << std::string(WO + 2, '-') << "+" << std::string(WN + 2, '-')
                  << "+\n";
    };
    for (const auto& r : results) {
        if (r.phase != phase) continue;
        if (!any) {
            rule();
            std::cout << "  | " << std::left << std::setw(WA) << "Actor" << " | "
                      << std::left << std::setw(WO) << "Operation" << " | "
                      << std::right << std::setw(WN) << "Time (ms)" << " |\n";
            rule();
            any = true;
        }
        std::cout << "  | " << std::left << std::setw(WA) << r.actor << " | "
                  << std::left << std::setw(WO) << r.op << " | " << std::right
                  << std::fixed << std::setprecision(3) << std::setw(WN) << r.ms
                  << " |\n";
        total += r.ms;
    }
    if (any) {
        rule();
        std::cout << "  Phase total: " << std::fixed << std::setprecision(3)
                  << total << " ms\n";
    }
}

// ── parameter-set table (128-bit Alg.2 LWE) ───────────────────────────────────
// k=2/4/8 use d=3, Bg=512; k=16 drops to d=4, Bg=128 (paper Table 3 row 4).
// Bks/dks are the key-switching base/digits used by the light KSK layout.
struct ParamRow {
    int             k;
    BINFHE_PARAMSET set;
    const char*     name;
    int      N      = 2048;
    int      n      = 635;
    long     B      = 512;     // Bg (gadget base for the 2nd-layer NTRU)
    int      l      = 3;       // gadget digits d = ceil(log2 Q / log2 Bg) = 3
    int      log2Q  = 27;      // Q = PreviousPrime(27, 4096)
    int      log2q  = 15;      // q = 32749
    int      Bks    = 32;      // KS gadget base
    int      dks    = 3;       // KS gadget digits
};
static const std::vector<ParamRow>& param_rows() {
    static const std::vector<ParamRow> rows = {
        {2,  STD128_MKNTRU_LWE,   "STD128_MKNTRU_LWE  (k=2)"},
        {4,  STD128_MKNTRU_LWE_2, "STD128_MKNTRU_LWE_2 (k=4)"},
        {8,  STD128_MKNTRU_LWE_3, "STD128_MKNTRU_LWE_3 (k=8)"},
        {16, STD128_MKNTRU_LWE_4, "STD128_MKNTRU_LWE_4 (k=16)",
            2048, 635, 128, 4, 27, 15, 32, 3},
    };
    return rows;
}

static std::string csv_param_prefix(const ParamRow& r) {
    std::ostringstream os;
    os << r.k << "," << r.N << "," << r.n << "," << r.B << "," << r.l << ","
       << r.log2Q << "," << r.log2q;
    return os.str();
}

// ── key-material size helpers ─────────────────────────────────────────────────

struct KeyMaterialBytes {
    std::size_t lwe_per_party = 0, lwe_total = 0;
    std::size_t bsk_per_party = 0, bsk_total = 0;
    std::size_t ksk_per_party = 0, ksk_total = 0;
};
static KeyMaterialBytes compute_key_bytes(const ParamRow& r) {
    KeyMaterialBytes b;
    const std::size_t N    = static_cast<std::size_t>(r.N);
    const std::size_t n    = static_cast<std::size_t>(r.n);
    const std::size_t d    = static_cast<std::size_t>(r.l);
    const std::size_t Bks  = static_cast<std::size_t>(r.Bks);
    const std::size_t dks  = static_cast<std::size_t>(r.dks);
    const std::size_t K    = static_cast<std::size_t>(r.k);
    const std::size_t bc   = BYTES_PER_COEFF;

    b.lwe_per_party = n;                                      // 1 byte per bit
    b.bsk_per_party = 2 * (n + 1) * (d - 1) * N * bc;         // UniEncACC layout
    b.ksk_per_party = 2 * N * (Bks - 1) * dks * bc;           // light KSK (RLWE packed)

    b.lwe_total = b.lwe_per_party * K;
    b.bsk_total = b.bsk_per_party * K;
    b.ksk_total = b.ksk_per_party * K;
    return b;
}
static double bytes_to_mb(std::size_t b) {
    return static_cast<double>(b) / (1024.0 * 1024.0);
}
static void print_key_size_table(const KeyMaterialBytes& b, int k) {
    std::cout << "\n  Key-material sizes (k=" << k << ", MB):\n";
    auto rule = []() {
        std::cout << "  +" << std::string(26, '-') << "+" << std::string(14, '-')
                  << "+" << std::string(14, '-') << "+\n";
    };
    rule();
    std::cout << "  | " << std::left  << std::setw(24) << "Primitive"
              << " | " << std::right << std::setw(11) << "per-party"
              << " | " << std::setw(11) << "total" << " |\n";
    rule();
    auto row = [&](const char* nm, std::size_t per_p, std::size_t tot) {
        std::cout << "  | " << std::left << std::setw(24) << nm << " | "
                  << std::right << std::fixed << std::setprecision(3)
                  << std::setw(8) << bytes_to_mb(per_p) << " MB"
                  << " | " << std::setw(8) << bytes_to_mb(tot) << " MB |\n";
    };
    row("MKLWE_KeyGen (LWE sk)",  b.lwe_per_party, b.lwe_total);
    row("MKBTKeyGen / BSK",       b.bsk_per_party, b.bsk_total);
    row("MKBTKeyGen / KSK",       b.ksk_per_party, b.ksk_total);
    rule();
    std::cout << "  Total per-party: "
              << std::fixed << std::setprecision(3)
              << bytes_to_mb(b.lwe_per_party + b.bsk_per_party + b.ksk_per_party)
              << " MB     Total (k=" << k << "): "
              << bytes_to_mb(b.lwe_total + b.bsk_total + b.ksk_total) << " MB\n";
}

// One row per k written to a SEPARATE keygen CSV (always emitted in every mode —
// that's where the user-facing keygen times live).
struct KeygenRecord {
    double ctx_ms      = 0;
    double lwe_kg_ms   = 0;
    double bt_kg_ms    = 0;
    KeyMaterialBytes b;
};
static void write_keygen_csv_row(std::ofstream& csv, const ParamRow& row,
                                 const KeygenRecord& kg)
{
    const double total_ms = kg.ctx_ms + kg.lwe_kg_ms + kg.bt_kg_ms;
    const std::size_t total_per_p = kg.b.lwe_per_party + kg.b.bsk_per_party + kg.b.ksk_per_party;
    const std::size_t total_all   = kg.b.lwe_total     + kg.b.bsk_total     + kg.b.ksk_total;
    csv << csv_param_prefix(row) << ","
        << std::fixed << std::setprecision(6)
        << kg.ctx_ms << "," << kg.lwe_kg_ms << "," << kg.bt_kg_ms << "," << total_ms << ","
        << kg.b.lwe_per_party << "," << kg.b.lwe_total << ","
        << kg.b.bsk_per_party << "," << kg.b.bsk_total << ","
        << kg.b.ksk_per_party << "," << kg.b.ksk_total << ","
        << total_per_p << "," << total_all << ","
        << bytes_to_mb(kg.b.lwe_per_party) << "," << bytes_to_mb(kg.b.lwe_total) << ","
        << bytes_to_mb(kg.b.bsk_per_party) << "," << bytes_to_mb(kg.b.bsk_total) << ","
        << bytes_to_mb(kg.b.ksk_per_party) << "," << bytes_to_mb(kg.b.ksk_total) << ","
        << bytes_to_mb(total_per_p)        << "," << bytes_to_mb(total_all)
        << "\n";
}
static void write_keygen_csv_header(std::ofstream& csv) {
    csv << "k,N,n,B,l,log2Q,log2q,"
           "ctx_ms,lwe_keygen_ms,mkbt_keygen_ms,total_ms,"
           "lwe_bytes_per_party,lwe_bytes_total,"
           "bsk_bytes_per_party,bsk_bytes_total,"
           "ksk_bytes_per_party,ksk_bytes_total,"
           "total_bytes_per_party,total_bytes_all,"
           "lwe_mb_per_party,lwe_mb_total,"
           "bsk_mb_per_party,bsk_mb_total,"
           "ksk_mb_per_party,ksk_mb_total,"
           "total_mb_per_party,total_mb_all\n";
}

// Shared setup: runs the three keygen primitives, prints sizes inline (so the
// user sees them as soon as each completes), populates `rec`, and returns sk.
static MKLWEPrivateKey run_setup_with_sizes(BinFHEContext&  cc,
                                            const ParamRow& row,
                                            KeygenRecord&   rec)
{
    rec.b = compute_key_bytes(row);

    {   auto t0 = Clock::now();
        cc.GenerateBinFHEContext(row.set, MKNTRU_LWE);
        rec.ctx_ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
        std::cout << "  GenerateBinFHEContext: " << std::fixed
                  << std::setprecision(3) << rec.ctx_ms << " ms\n";
    }
    MKLWEPrivateKey sk;
    {   auto t0 = Clock::now();
        sk = cc.MKLWE_KeyGen();
        rec.lwe_kg_ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
        std::cout << "  MKLWE_KeyGen:          " << std::fixed
                  << std::setprecision(3) << rec.lwe_kg_ms
                  << " ms     (LWE sk per party = "
                  << std::setprecision(3) << bytes_to_mb(rec.b.lwe_per_party) << " MB"
                  << ", total = " << bytes_to_mb(rec.b.lwe_total) << " MB)\n";
    }
    {   std::cout << "  MKBTKeyGen running (slow)... " << std::flush;
        auto t0 = Clock::now();
        cc.MKBTKeyGen(sk);
        rec.bt_kg_ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
        std::cout << "done in " << std::fixed << std::setprecision(3)
                  << rec.bt_kg_ms << " ms\n";
        std::cout << "    -> BSK per party = " << std::setprecision(3)
                  << bytes_to_mb(rec.b.bsk_per_party) << " MB  (total "
                  << bytes_to_mb(rec.b.bsk_total) << " MB)\n";
        std::cout << "    -> KSK per party = "
                  << bytes_to_mb(rec.b.ksk_per_party) << " MB  (total "
                  << bytes_to_mb(rec.b.ksk_total) << " MB)\n";
    }
    print_key_size_table(rec.b, row.k);
    return sk;
}

// ════════════════════════════════════════════════════════════════════════════
//  LWE noise helpers (public-API only — used by --components and --growth)
// ════════════════════════════════════════════════════════════════════════════
//
// All noise is read straight off the ciphertext with the public accessors
// MKLWECiphertext::GetA()/GetB()/GetModulus()/Getk()/GetLength() and
// MKLWEPrivateKey::GetElement(). No library instrumentation is required.
//
// MKLWE encodes a bit m as  b = <a,s> + e + m*(q/4)  (mklwe-pke.cpp). The
// centred phase v = (b - <a,s>) mod q therefore equals m*(q/4) + e, so for a
// zero message the phase IS the signed noise, and for any message the distance
// to the nearest codeword in {0, q/4} is the noise magnitude.

// Centred phase v = b - sum_i <a_i, s_i>  in (-q/2, q/2], summed over parties.
// Mirrors MKLWEEncryptionScheme::Decrypt's inner-product accumulation exactly.
static long lwe_centered_phase(const MKLWEPrivateKey& sk, const MKLWECiphertext& ct) {
    const NativeInteger q = ct->GetModulus();
    const long          qi = static_cast<long>(q.ConvertToInt());
    const uint32_t      k  = ct->Getk();
    const uint32_t      n  = ct->GetLength();
    const auto&         a  = ct->GetA();

    std::vector<NativeVector> s = sk->GetElement();  // copy so we may SwitchModulus
    const NativeInteger mu = q.ComputeMu();
    NativeInteger inner(0);
    for (uint32_t u = 0; u < k; ++u) {
        s[u].SwitchModulus(q);
        for (uint32_t j = 0; j < n; ++j)
            inner += a[u][j].ModMulFast(s[u][j], q, mu);
        inner.ModEq(q);
    }
    NativeInteger r = ct->GetB();
    r.ModSubFastEq(inner, q);  // (b - <a,s>) mod q
    long v = static_cast<long>(r.ConvertToInt());
    if (v > qi / 2) v -= qi;   // centre into (-q/2, q/2]
    return v;
}

// Noise magnitude = distance from the centred phase to the nearest codeword in
// {0, Delta} (Delta = q/4), i.e. |e| for whichever bit is closest.
static long lwe_noise_magnitude(long v, long qi) {
    const long Delta = qi / 4;
    const long d0    = std::labs(v);
    const long d1    = std::min(std::labs(v - Delta), std::labs(v + Delta));
    return std::min(d0, d1);
}

// Signed noise = phase minus the nearest codeword in {0, +/-Delta}. Used so a
// pool of bootstrap outputs (which may decode to 0 or 1) yields a meaningful
// zero-mean stddev rather than a bimodal spread around 0 and Delta.
static long lwe_signed_noise(long v, long qi) {
    const long Delta   = qi / 4;
    const long cand[3] = {v, v - Delta, v + Delta};
    long best = cand[0];
    for (long c : cand)
        if (std::labs(c) < std::labs(best)) best = c;
    return best;
}

// Manual homomorphic Add (component-wise mod q). MKFHE exposes no public linear
// EvalAdd at the BinFHE layer, so we add via the public accessors. Mirrors
// MKLWEEncryptionScheme::EvalAddEq.
static MKLWECiphertext mk_add(const MKLWECiphertext& c1, const MKLWECiphertext& c2) {
    const NativeInteger q = c1->GetModulus();
    const uint32_t      k = c1->Getk();
    std::vector<NativeVector> a = c1->GetA();  // copy
    const auto& a2 = c2->GetA();
    for (uint32_t u = 0; u < k; ++u)
        a[u].ModAddEq(a2[u]);
    NativeInteger b = c1->GetB();
    b.ModAddFastEq(c2->GetB(), q);
    auto ct = std::make_shared<MKLWECiphertextImpl>(std::move(a), b);
    ct->SetptModulus(c1->GetptModulus());
    return ct;
}

// ── noise statistics (shared by --components) ─────────────────────────────────
// Streaming accumulator: never stores the samples (count can be large), keeps a
// running mean / variance and the worst-case max-abs.
struct NoiseStats {
    double   sum     = 0.0;
    double   sumsq   = 0.0;
    double   max_abs = 0.0;
    uint64_t count   = 0;
    void push(double x) {
        sum   += x;
        sumsq += x * x;
        const double ax = std::fabs(x);
        if (ax > max_abs) max_abs = ax;
        ++count;
    }
    double mean()   const { return count ? sum / static_cast<double>(count) : 0.0; }
    double stddev() const {
        if (count < 2) return 0.0;
        const double m = mean();
        const double v = sumsq / static_cast<double>(count) - m * m;
        return v > 0.0 ? std::sqrt(v) : 0.0;
    }
};

// P[|e| > half_margin] for centred sub-gaussian noise of width sigma:
// erfc(half_margin / (sigma * sqrt2)).
static double erfc_fail_prob(double sigma, double half_margin) {
    if (sigma <= 0.0)       return 0.0;
    if (half_margin <= 0.0) return 1.0;
    return std::erfc(half_margin / (sigma * std::sqrt(2.0)));
}
static double log2_abs(double x) {
    const double ax = std::fabs(x);
    return ax > 0.0 ? std::log2(ax) : 0.0;
}
static void print_noise_row(const std::string& name, const std::string& domain,
                            const NoiseStats& st, double half_margin) {
    const double sigma = st.stddev();
    const double pfail = erfc_fail_prob(sigma, half_margin);
    std::cout << "  | " << std::left  << std::setw(28) << name
              << " | " << std::left  << std::setw(3)  << domain
              << " | " << std::right << std::setw(8)  << std::fixed << std::setprecision(2)
              << log2_abs(st.max_abs)
              << " | " << std::setw(9)  << std::setprecision(2) << log2_abs(sigma)
              << " | " << std::setw(8)  << std::setprecision(2) << log2_abs(half_margin)
              << " | " << std::setw(11) << std::scientific << std::setprecision(2) << pfail
              << std::fixed << " |\n";
}

// ── timing-mode CSV writers (FINALLY schema) ──────────────────────────────────

static void write_timing_row(std::ofstream& csv, const ParamRow& row,
                             const std::string& op_name,
                             std::size_t count_per_bootstrap,
                             const TimingSummary& t)
{
    const double projected = t.median * static_cast<double>(count_per_bootstrap);
    csv << csv_param_prefix(row) << ",\"" << op_name << "\","
        << count_per_bootstrap << "," << t.n << ","
        << std::fixed << std::setprecision(6)
        << t.min << "," << t.median << "," << t.p90 << "," << t.max << ","
        << t.mean << "," << t.sd << "," << projected << "\n";
}

// HOOK POINT — reads the library timing profile (binfhe-timing.{h,cpp}). When
// the library is built with the ScopedTimer brackets these four rows carry real
// per-component numbers; without them the counts are zero and the rows skip.
static bool emit_components_if_instrumented(std::ofstream& csv,
                                            const ParamRow& row,
                                            BinFHEContext&  cc,
                                            int             reps,
                                            MKLWECiphertext& ct1,
                                            MKLWECiphertext& ct2)
{
    (void)cc;
    lbcrypto::ResetTimingProfile();
    for (int i = 0; i < reps; ++i) (void)cc.EvalBinGate(NAND, ct1, ct2);
    const auto tp = lbcrypto::GetTimingProfile();
    std::cerr << "  [debug] counts: cmux=" << tp.cmux_count
              << " ext=" << tp.ext_count
              << " ks="  << tp.ks_count
              << " ms="  << tp.ms_count << "\n";

    auto emit = [&](const char* nm, std::size_t count_per_boot,
                    uint64_t total_ns, uint64_t calls) {
        if (!calls) return;
        const double per_call_ms =
            (static_cast<double>(total_ns) / 1e6) / static_cast<double>(calls);
        TimingSummary t;
        t.n = calls;
        t.min = t.median = t.p90 = t.max = t.mean = per_call_ms;
        t.sd  = 0.0;
        write_timing_row(csv, row, nm, count_per_boot, t);
    };
    const std::size_t n_cmux = static_cast<std::size_t>(row.n) * row.k;
    emit("CMux+ExtProduct", n_cmux, tp.cmux_ns, tp.cmux_count);
    emit("ExternalProduct", n_cmux, tp.ext_ns,  tp.ext_count);
    emit("KS_NGS_to_LWE",   1,      tp.ks_ns,   tp.ks_count);
    emit("ModSwitch",       1,      tp.ms_ns,   tp.ms_count);
    return (tp.cmux_count || tp.ks_count || tp.ms_count);
}

// ── --timing mode ─────────────────────────────────────────────────────────────

static int run_timing(const ParamRow& row, int reps,
                      std::ofstream& timing_csv,
                      std::ofstream& keygen_csv)
{
    std::cout << "\n========================================================\n"
              << "  MKFHE Alg.2 (LWE) - Timing Benchmark  (k=" << row.k << ")\n"
              << "  Parameter set: " << row.name
              << "    reps=" << reps << "\n"
              << "========================================================\n";

    auto cc = BinFHEContext();
    KeygenRecord kg;
    MKLWEPrivateKey sk = run_setup_with_sizes(cc, row, kg);
    write_keygen_csv_row(keygen_csv, row, kg);

    MKLWECiphertext ct1 = cc.Encrypt(sk, 1);
    MKLWECiphertext ct2 = cc.Encrypt(sk, 0);

    // ── EvalBinGate (full) — the only public-API primitive ────────────────────
    {
        (void)cc.EvalBinGate(NAND, ct1, ct2);   // warmup

        std::vector<double> samples;
        samples.reserve(static_cast<std::size_t>(reps));
        std::cout << "\n  timing  " << std::left << std::setw(24) << "EvalBinGate (full)"
                  << " ..." << std::flush;
        for (int r = 0; r < reps; ++r) {
            auto t0 = Clock::now();
            (void)cc.EvalBinGate(NAND, ct1, ct2);
            auto t1 = Clock::now();
            samples.push_back(
                std::chrono::duration<double, std::milli>(t1 - t0).count());
        }
        TimingSummary t = summarize_ms(samples);
        std::cout << "\r  timed   " << std::left << std::setw(24) << "EvalBinGate (full)"
                  << "  min=" << std::fixed << std::setprecision(3) << t.min
                  << "  median=" << t.median << "  p90=" << t.p90
                  << "  max=" << t.max << "  stddev=" << t.sd << " ms\n";
        write_timing_row(timing_csv, row, "EvalBinGate (full)", 1, t);
    }

    const bool emitted = emit_components_if_instrumented(timing_csv, row, cc, reps, ct1, ct2);
    if (!emitted) {
        std::cout << "  [note] MKFHE not instrumented for per-component timing;\n"
                  << "         only the 'EvalBinGate (full)' row was written.\n";
    }
    return 0;
}

// ── --components mode : per-stage LWE noise + erfc decode error ────────────────
//
// FINALLY's --components isolates the four NGS bootstrap primitives because it
// exposes them. MKFHE keeps them inside EvalBinGate, so here we measure the LWE
// noise at the three stages the public API DOES expose:
//
//   Encrypt(0) fresh    : noise a fresh ciphertext carries (LWE, margin q/8).
//   Add (Enc0 + Enc0)   : noise after one homomorphic Add  (LWE, margin q/8).
//   NAND+Bootstrap out  : noise of a bootstrapped gate output — the reset level
//                         (LWE, margin q/8).
//
// margin = q/8 = Delta/2 is the half-window a coefficient may drift before the
// decoded bit flips.
// Core measurement only: keys (cc, sk) are generated ONCE by the caller and
// reused here, so a combined run does not pay MKBTKeyGen twice.
static int run_components(const ParamRow& row, int reps,
                          std::ofstream& comp_csv,
                          BinFHEContext& cc, const MKLWEPrivateKey& sk)
{
    std::cout << "\n  --- Components (per-stage LWE noise) benchmark  (k=" << row.k
              << ", reps=" << reps << ") ---\n";

    const auto& P  = cc.GetParams()->GetMKLWEParams();
    const long  qi = static_cast<long>(P->Getq().ConvertToInt());
    const double lwe_margin = static_cast<double>(qi) / 8.0;   // Delta/2 = q/8

    NoiseStats st_fresh, st_add, st_boot;

    std::cout << "\n  measuring per-stage LWE noise over " << reps << " reps/stage ..."
              << std::flush;
    for (int r = 0; r < reps; ++r) {
        // Fresh Encrypt(0): message 0 => centred phase is the signed noise.
        MKLWECiphertext c0a = cc.Encrypt(sk, 0);
        st_fresh.push(static_cast<double>(lwe_centered_phase(sk, c0a)));

        // Homomorphic Add of two fresh Enc(0): still message 0.
        MKLWECiphertext c0b = cc.Encrypt(sk, 0);
        MKLWECiphertext cadd = mk_add(c0a, c0b);
        st_add.push(static_cast<double>(lwe_centered_phase(sk, cadd)));

        // NAND+Bootstrap output: NAND(1,0)=1 => codeword q/4. Report the signed
        // noise relative to the nearest codeword so stddev is zero-mean.
        MKLWECiphertext a1 = cc.Encrypt(sk, 1);
        MKLWECiphertext b0 = cc.Encrypt(sk, 0);
        MKLWECiphertext out = cc.EvalBinGate(NAND, a1, b0);
        st_boot.push(static_cast<double>(lwe_signed_noise(lwe_centered_phase(sk, out), qi)));
    }
    std::cout << " done\n";

    // ── Table ─────────────────────────────────────────────────────────────────
    std::cout << "\n  Per-stage LWE noise + decode error rate (magnitudes in log2 bits):\n";
    auto rule = []() {
        std::cout << "  +" << std::string(30, '-') << "+" << std::string(5, '-')
                  << "+" << std::string(10, '-') << "+" << std::string(11, '-')
                  << "+" << std::string(10, '-') << "+" << std::string(13, '-') << "+\n";
    };
    rule();
    std::cout << "  | " << std::left << std::setw(28) << "Stage"
              << " | " << std::setw(3) << "dom"
              << " | " << std::right << std::setw(8) << "max|e|"
              << " | " << std::setw(9) << "sigma"
              << " | " << std::setw(8) << "margin"
              << " | " << std::setw(11) << "P[fail]" << " |\n";
    rule();
    print_noise_row("Encrypt(0) fresh",   "LWE", st_fresh, lwe_margin);
    print_noise_row("Add (Enc0 + Enc0)",  "LWE", st_add,   lwe_margin);
    print_noise_row("NAND+Bootstrap out", "LWE", st_boot,  lwe_margin);
    rule();
    std::cout << "  Notes: noise read from the centred phase via public accessors;\n"
                 "         max|e|/sigma in log2 bits; margin = q/8 decode half-window;\n"
                 "         P[fail] = erfc(margin / (sigma*sqrt2)) per-stage tail prob.\n";

    // ── CSV ───────────────────────────────────────────────────────────────────
    auto write_row = [&](const std::string& nm, const std::string& dom,
                         const NoiseStats& st, double margin) {
        comp_csv << csv_param_prefix(row) << ",\"" << nm << "\",\"" << dom << "\","
                 << reps << "," << st.count << ","
                 << std::scientific << std::setprecision(8)
                 << st.mean() << "," << st.stddev() << "," << st.max_abs << ","
                 << std::fixed << std::setprecision(4)
                 << log2_abs(st.stddev()) << "," << log2_abs(margin) << ","
                 << std::scientific << std::setprecision(8)
                 << margin << "," << erfc_fail_prob(st.stddev(), margin) << "\n";
    };
    write_row("Encrypt(0) fresh",   "LWE", st_fresh, lwe_margin);
    write_row("Add (Enc0 + Enc0)",  "LWE", st_add,   lwe_margin);
    write_row("NAND+Bootstrap out", "LWE", st_boot,  lwe_margin);
    return 0;
}

// ── --growth mode : LWE noise-budget trajectory (grow via Add, reset via Boot) ─
//
// Encrypt the bit 0, add fresh Enc(0) ciphertexts (noise grows ~sqrt(count)),
// sampling the LWE noise (log2 bits) at a geometric ladder of add-counts, then
// reset via EvalBinGate(NAND, acc, Enc(1)) — NAND(0,1)=1 — which re-bootstraps
// to a fixed low noise. All operands encrypt 0, so the centred phase IS the
// noise and stays below the decode budget log2(q/8); the chart shows the
// growth->reset sawtooth.
// Core measurement only: keys (cc, sk) are generated ONCE by the caller and
// reused here (no second MKBTKeyGen for a combined run).
static int run_noise_growth(const ParamRow& row, int reps,
                            std::ofstream& growth_csv,
                            BinFHEContext& cc, const MKLWEPrivateKey& sk)
{
    std::cout << "\n  --- Noise-growth (grow via Add, reset via Bootstrap) benchmark  (k="
              << row.k << ", reps=" << reps << ") ---\n";

    const auto& P  = cc.GetParams()->GetMKLWEParams();
    const long  qi = static_cast<long>(P->Getq().ConvertToInt());
    const double budget_bits = std::log2(static_cast<double>(qi) / 8.0);
    const std::vector<std::size_t> add_counts = {1, 2, 4, 8, 16, 32, 64, 128, 256};

    auto noise_bits = [&](const MKLWECiphertext& ct) -> double {
        const long mag = lwe_noise_magnitude(lwe_centered_phase(sk, ct), qi);
        return mag > 0 ? std::log2(static_cast<double>(mag)) : 0.0;
    };

    std::cout << "\n  --growth: " << reps
              << " reps  (LWE noise in log2 bits; decode budget = log2(q/8) = "
              << std::fixed << std::setprecision(1) << budget_bits << " bits)\n";

    std::vector<std::vector<double>> growth(add_counts.size());
    std::vector<double> boot_samples;

    for (int r = 0; r < reps; ++r) {
        MKLWECiphertext acc;
        bool init = false;
        std::size_t ci = 0;
        for (std::size_t count = 1; count <= add_counts.back(); ++count) {
            MKLWECiphertext ct = cc.Encrypt(sk, 0);
            acc = init ? mk_add(acc, ct) : ct;
            init = true;
            if (ci < add_counts.size() && count == add_counts[ci]) {
                growth[ci].push_back(noise_bits(acc));
                ++ci;
            }
        }
        // Bootstrap the (still-decrypts-to-0) accumulator: NAND(0,1)=1 resets noise.
        MKLWECiphertext one = cc.Encrypt(sk, 1);
        MKLWECiphertext out = cc.EvalBinGate(NAND, acc, one);
        boot_samples.push_back(noise_bits(out));
    }

    // ── Table ─────────────────────────────────────────────────────────────────
    std::cout << "\n  Noise growth trajectory (log2 |e| bits):\n";
    auto rule = []() {
        std::cout << "  +" << std::string(26, '-') << "+" << std::string(10, '-')
                  << "+" << std::string(10, '-') << "+" << std::string(10, '-')
                  << "+" << std::string(10, '-') << "+\n";
    };
    rule();
    std::cout << "  | " << std::left << std::setw(24) << "Step"
              << " | " << std::right << std::setw(8) << "min"
              << " | " << std::setw(8) << "mean"
              << " | " << std::setw(8) << "max"
              << " | " << std::setw(8) << "%budget" << " |\n";
    rule();
    auto print_row = [&](const std::string& nm, const TimingSummary& s) {
        std::cout << "  | " << std::left << std::setw(24) << nm
                  << " | " << std::right << std::fixed << std::setprecision(2)
                  << std::setw(8) << s.min
                  << " | " << std::setw(8) << s.mean
                  << " | " << std::setw(8) << s.max
                  << " | " << std::setw(7) << std::setprecision(1)
                  << (budget_bits > 0 ? 100.0 * s.mean / budget_bits : 0.0) << "%" << " |\n";
    };
    std::vector<TimingSummary> gstats(add_counts.size());
    for (std::size_t i = 0; i < add_counts.size(); ++i) {
        gstats[i] = summarize_ms(growth[i]);
        const std::string lbl = (i == 0)
            ? "encrypt (fresh)"
            : "+ Add (" + std::to_string(add_counts[i]) + " cts)";
        print_row(lbl, gstats[i]);
    }
    rule();
    const TimingSummary bstat = summarize_ms(boot_samples);
    print_row("bootstrap (reset)", bstat);
    rule();
    std::cout << "  Add grows noise ~sqrt(count) (~+0.5 bits per doubling); "
                 "Bootstrap resets it.\n";

    // ── CSV ───────────────────────────────────────────────────────────────────
    auto write_row = [&](const std::string& nm, std::size_t n_added,
                         const TimingSummary& s) {
        growth_csv << csv_param_prefix(row) << ",\"" << nm << "\"," << n_added << ","
                   << s.n << "," << std::fixed << std::setprecision(6)
                   << s.min << "," << s.mean << "," << s.max << "," << s.sd << ","
                   << budget_bits << "\n";
    };
    for (std::size_t i = 0; i < add_counts.size(); ++i)
        write_row(i == 0 ? "encrypt" : "add", add_counts[i], gstats[i]);
    write_row("bootstrap", 0, bstat);
    return 0;
}

// ── noise driver : ONE keygen per k, then run growth and/or components ─────────
//
// Generating the bootstrapping keys (MKBTKeyGen) is by far the most expensive
// step, so it is done exactly once here and shared by both noise benches. This
// is why --growth, --noise, and --noise-all all route through this driver: a
// combined run pays a single keygen per k (and writes a single keygen row).
static int run_noise_for_k(const ParamRow& row, int reps,
                           bool do_growth, bool do_components,
                           std::ofstream* growth_csv, std::ofstream* comp_csv,
                           std::ofstream& keygen_csv)
{
    std::cout << "\n========================================================\n"
              << "  MKFHE Alg.2 (LWE) - Noise Benchmark  (k=" << row.k << ")\n"
              << "  Parameter set: " << row.name << "    reps=" << reps << "\n"
              << "  Modes:" << (do_growth ? " noise-growth" : "")
              << (do_components ? " components" : "") << "\n"
              << "========================================================\n";

    auto cc = BinFHEContext();
    KeygenRecord kg;
    MKLWEPrivateKey sk = run_setup_with_sizes(cc, row, kg);  // one keygen per k
    write_keygen_csv_row(keygen_csv, row, kg);               // one keygen row per k

    int rc = 0;
    if (do_growth     && growth_csv) rc |= run_noise_growth(row, reps, *growth_csv, cc, sk);
    if (do_components && comp_csv)    rc |= run_components(row, reps, *comp_csv, cc, sk);
    return rc;
}

// ── default protocol mode ─────────────────────────────────────────────────────

static int run_protocol(const ParamRow& row, int reps,
                        std::ofstream* protocol_csv,
                        std::ofstream& keygen_csv)
{
    std::cout << "\n========================================================\n"
              << "  MKFHE Alg.2 (LWE) - Protocol Benchmark  (k=" << row.k << ")\n"
              << "  Parameter set: " << row.name << "\n"
              << "  Scenario: P1 holds m1, P2 holds m2; "
              << row.k << "-party blind rotation\n"
              << "========================================================\n";

    std::vector<OpResult> results;

    print_phase_header("Phase A - Setup (key + bootstrapping-key generation)");
    auto cc = BinFHEContext();
    KeygenRecord kg;
    kg.b = compute_key_bytes(row);
    MKLWEPrivateKey sk;
    set_total(3);
    results.push_back(run_op("A", "All", "GenerateBinFHEContext", [&]() {
        cc.GenerateBinFHEContext(row.set, MKNTRU_LWE);
    }));
    kg.ctx_ms = results.back().ms;
    results.push_back(run_op("A", "Party", "MKLWE_KeyGen", [&]() {
        sk = cc.MKLWE_KeyGen();
    }));
    kg.lwe_kg_ms = results.back().ms;
    std::cout << "    -> LWE sk per party = " << std::fixed << std::setprecision(3)
              << bytes_to_mb(kg.b.lwe_per_party) << " MB  (total "
              << bytes_to_mb(kg.b.lwe_total) << " MB)\n";
    results.push_back(run_op("A", "Party", "MKBTKeyGen (bootstrapping keys)", [&]() {
        cc.MKBTKeyGen(sk);
    }));
    kg.bt_kg_ms = results.back().ms;
    std::cout << "    -> BSK per party    = "
              << bytes_to_mb(kg.b.bsk_per_party) << " MB  (total "
              << bytes_to_mb(kg.b.bsk_total) << " MB)\n";
    std::cout << "    -> KSK per party    = "
              << bytes_to_mb(kg.b.ksk_per_party) << " MB  (total "
              << bytes_to_mb(kg.b.ksk_total) << " MB)\n";

    print_phase_table("A", results);
    print_key_size_table(kg.b, row.k);
    write_keygen_csv_row(keygen_csv, row, kg);

    print_phase_header("Phase B - Encryption (all four input combinations)");
    const int inputs[4][2] = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    MKLWECiphertext ct_p1[4], ct_p2[4];
    set_total(8);
    for (int r = 0; r < 4; ++r) {
        const int m1 = inputs[r][0], m2 = inputs[r][1];
        results.push_back(run_op("B", "P1",
            "Encrypt P1 m1=" + std::to_string(m1) + "  (row " + std::to_string(r) + ")",
            [&, r, m1]() { ct_p1[r] = cc.Encrypt(sk, m1); }));
        results.push_back(run_op("B", "P2",
            "Encrypt P2 m2=" + std::to_string(m2) + "  (row " + std::to_string(r) + ")",
            [&, r, m2]() { ct_p2[r] = cc.Encrypt(sk, m2); }));
    }
    print_phase_table("B", results);

    print_phase_header("Phase C - Bootstrapped NAND evaluation (server)");
    MKLWECiphertext ct_nand[4];
    set_total(1);
    double nand_boot_ms = 0;
    results.push_back(run_op("C", "Server", "NAND+Bootstrap (row 2: m1=1,m2=0)", [&]() {
        ct_nand[2] = cc.EvalBinGate(NAND, ct_p1[2], ct_p2[2]);
    }));
    nand_boot_ms = results.back().ms;
    for (int r = 0; r < 4; ++r) {
        if (r == 2) continue;
        ct_nand[r] = cc.EvalBinGate(NAND, ct_p1[r], ct_p2[r]);
    }
    print_phase_table("C", results);

    std::cout << "\n  NAND truth table (Bootstrap resets noise each gate):\n";
    std::cout << "    m1  m2 | expected | decrypts | noise | status\n";
    std::cout << "    -------+----------+----------+-------+-------\n";
    int nand_pass = 0;
    {
        const long qi = static_cast<long>(cc.GetParams()->GetMKLWEParams()->Getq().ConvertToInt());
        for (int r = 0; r < 4; ++r) {
            const int m1 = inputs[r][0], m2 = inputs[r][1];
            const int expected = 1 - (m1 & m2);
            MKLWEPlaintext got = 0;
            cc.Decrypt(sk, ct_nand[r], &got);
            const long noise = lwe_noise_magnitude(lwe_centered_phase(sk, ct_nand[r]), qi);
            const bool ok = (static_cast<int>(got) == expected);
            if (ok) ++nand_pass;
            std::cout << "     " << m1 << "   " << m2 << " |    " << expected
                      << "     |    " << static_cast<int>(got) << "     | "
                      << std::setw(5) << noise << " | "
                      << (ok ? "OK" : "FAIL") << "\n";
        }
    }
    std::cout << "    " << nand_pass << "/4 correct\n";

    print_phase_header("Phase D - Level-2 gate + bootstrap timing (server)");
    set_total(1);
    MKLWECiphertext ct_level2;
    {
        std::vector<double> samples;
        samples.reserve(static_cast<std::size_t>(reps));
        ct_level2 = cc.EvalBinGate(NAND, ct_nand[0], ct_nand[3]);
        results.push_back(run_op("D", "Server",
            "NAND(Boot(NAND(0,0)),Boot(NAND(1,1)))+Bootstrap  [min of "
            + std::to_string(reps) + "]",
            [&]() {
            for (int r = 0; r < reps; ++r) {
                auto t0 = Clock::now();
                ct_level2 = cc.EvalBinGate(NAND, ct_nand[0], ct_nand[3]);
                auto t1 = Clock::now();
                samples.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
            }
        }));
        TimingSummary t = summarize_ms(samples);
        std::cout << "    Bootstrap timing over " << reps << " reps:  min="
                  << std::fixed << std::setprecision(2) << t.min
                  << "  median=" << t.median << "  p90=" << t.p90
                  << "  max=" << t.max << "  mean=" << t.mean
                  << "  stddev=" << t.sd << "  ms\n";
        if (protocol_csv) {
            (*protocol_csv) << row.k << "," << t.n << "," << std::fixed << std::setprecision(4)
                            << t.min << "," << t.median << "," << t.p90 << "," << t.max << ","
                            << t.mean << "," << t.sd << "," << nand_boot_ms << "\n";
        }
    }
    print_phase_table("D", results);

    print_phase_header("Phase E - Decryption + correctness (parties)");
    set_total(1);
    MKLWEPlaintext result_timed = 0;
    results.push_back(run_op("E", "Both", "Decrypt NAND(1,0)+Bootstrap  [timed]", [&]() {
        cc.Decrypt(sk, ct_nand[2], &result_timed);
    }));
    print_phase_table("E", results);

    MKLWEPlaintext got_level2 = 0;
    cc.Decrypt(sk, ct_level2, &got_level2);
    const bool level2_ok = (static_cast<int>(got_level2) == 1);
    std::cout << "\n  Correctness summary:\n";
    std::cout << "    NAND+Bootstrap (all 4 inputs): " << nand_pass << "/4 correct\n";
    std::cout << "    Level-2 NAND+Bootstrap:        " << static_cast<int>(got_level2)
              << "  (expected 1)  [" << (level2_ok ? "PASS" : "FAIL") << "]\n";
    const bool scheme_ok = (nand_pass == 4 && level2_ok);
    std::cout << "    Overall: " << (scheme_ok ? "[PASS] scheme correct"
                                              : "[FAIL] failures detected") << "\n";

    auto phase_total = [&](const std::string& ph) {
        double t = 0;
        for (const auto& r : results) if (r.phase == ph) t += r.ms;
        return t;
    };
    const double t_setup = phase_total("A");
    const double t_enc   = phase_total("B");
    const double t_dec   = phase_total("E");
    std::cout << "\n=== Protocol summary (k=" << row.k << ") ===\n\n";
    std::cout << "  k (parties in parameter set):      " << row.k << "\n";
    std::cout << "  One-time setup (KeyGen + MKBTKeyGen): "
              << std::fixed << std::setprecision(1) << t_setup << " ms\n";
    std::cout << "  Encrypt one bit (avg over 8):      " << t_enc / 8.0 << " ms\n";
    std::cout << "  NAND + Bootstrap per gate:         " << std::setprecision(3)
              << nand_boot_ms << " ms\n";
    std::cout << "  Decrypt one bit:                   " << t_dec << " ms\n";
    std::cout << "\n  Projected cost of " << N_GATES << " bootstrapped NAND gates:\n";
    std::cout << "    " << N_GATES << " x " << std::setprecision(1) << nand_boot_ms
              << " ms = " << nand_boot_ms * N_GATES / 1000.0 << " s\n";

    return scheme_ok ? 0 : 1;
}

// ── CLI ───────────────────────────────────────────────────────────────────────

enum class Mode { Protocol, Timing, Components, Growth, NoiseAll };

static void print_usage(const char* prog) {
    std::cout <<
        "MKFHE Alg.2 (LWE) - Protocol/Timing/Noise/Keygen-sizing Benchmark\n\n"
        "Usage: " << prog << " [k] [options]\n\n"
        "Positional:\n"
        "  k             Restrict to a single k in {2, 4, 8, 16} (default: run all)\n\n"
        "Modes (default = full protocol scenario):\n"
        "  --timing                  Per-op timing.   CSV -> bench_timing_mklwe.csv\n"
        "  --components, --noise      Per-stage LWE noise + erfc decode error.\n"
        "                            CSV -> bench_components_mklwe.csv\n"
        "  --growth, --noise-growth   LWE noise-budget trajectory (grow via Add,\n"
        "                            reset via Bootstrap). CSV -> bench_growth_mklwe.csv\n"
        "  --noise-all                Both noise benches sharing ONE keygen per k\n"
        "                            (avoids regenerating bootstrapping keys twice).\n"
        "                            CSVs -> bench_growth_mklwe.csv + bench_components_mklwe.csv\n\n"
        "Options:\n"
        "  --reps N         Sample count (default 20; --timing defaults to 100).\n"
        "  --csv PATH       Write the active mode's CSV to PATH instead of the\n"
        "                   default bench_<mode>_mklwe.csv (single-mode runs only;\n"
        "                   final name written directly -- no rename needed).\n"
        "  --growth-csv P   --noise-all: growth CSV path (default bench_growth_mklwe.csv).\n"
        "  --noise-csv P    --noise-all: components CSV path (default bench_components_mklwe.csv).\n"
        "  --keygen-csv P   Write the keygen CSV to P instead of bench_keygen_mklwe.csv.\n"
        "  -h, --help       Show this help.\n\n"
        "EVERY mode also writes a keygen CSV (one row per k with\n"
        "GenerateBinFHEContext / MKLWE_KeyGen / MKBTKeyGen timings and per-key\n"
        "material sizes in bytes + MB); default name bench_keygen_mklwe.csv.\n\n"
        "Default mode runs the full protocol scenario (Setup, Encryption,\n"
        "NAND+Bootstrap, level-2 gate, Decryption, summary).\n"
        "CSV -> bench_protocol_mklwe.csv\n";
}

int main(int argc, char* argv[]) {
    Mode mode   = Mode::Protocol;
    int  reps   = -1;
    int  only_k = 0;
    std::string csv_path;        // --csv: override the active (single) mode's CSV path
    std::string growth_path;     // --growth-csv: --noise-all growth CSV path
    std::string noise_path;      // --noise-csv:  --noise-all components CSV path
    std::string keygen_path;     // --keygen-csv: override the keygen CSV path

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "-h" || a == "--help") { print_usage(argv[0]); return 0; }
        else if (a == "--timing")                      { mode = Mode::Timing; }
        else if (a == "--components" || a == "--noise"){ mode = Mode::Components; }
        else if (a == "--growth" || a == "--noise-growth") { mode = Mode::Growth; }
        else if (a == "--noise-all")                   { mode = Mode::NoiseAll; }
        else if (a == "--reps" && i + 1 < argc)        { reps = std::atoi(argv[++i]); }
        else if (a == "--csv" && i + 1 < argc)         { csv_path = argv[++i]; }
        else if (a == "--growth-csv" && i + 1 < argc)  { growth_path = argv[++i]; }
        else if (a == "--noise-csv" && i + 1 < argc)   { noise_path = argv[++i]; }
        else if (a == "--keygen-csv" && i + 1 < argc)  { keygen_path = argv[++i]; }
        else if (a.size() > 0 && a[0] != '-')          { only_k = std::atoi(a.c_str()); }
        else {
            std::cerr << "Unknown option: " << a << "\n";
            print_usage(argv[0]); return 1;
        }
    }
    if (reps < 0) reps = (mode == Mode::Timing) ? 100 : 20;
    if (reps < 1) reps = 1;

    // CSV destinations: --csv / --keygen-csv override the defaults so the wrapper
    // can ask for the final filename directly (no fixed-name + rename dance).
    auto pick = [](const std::string& override_path, const char* def) {
        return override_path.empty() ? std::string(def) : override_path;
    };
    const std::string keygen_name = pick(keygen_path, "bench_keygen_mklwe.csv");

    // Keygen CSV is written in EVERY mode.
    std::ofstream keygen_csv(keygen_name);
    write_keygen_csv_header(keygen_csv);

    int rc = 0;
    switch (mode) {
    case Mode::Timing: {
        const std::string name = pick(csv_path, "bench_timing_mklwe.csv");
        std::ofstream timing_csv(name);
        timing_csv << "k,N,n,B,l,log2Q,log2q,operation,count_per_bootstrap,reps,"
                      "min_ms,median_ms,p90_ms,max_ms,mean_ms,stddev_ms,projected_ms\n";
        for (const auto& row : param_rows()) {
            if (only_k != 0 && row.k != only_k) continue;
            rc |= run_timing(row, reps, timing_csv, keygen_csv);
        }
        std::cout << "\n[bench] timing     CSV -> " << name << "\n";
        break;
    }
    case Mode::Components: {
        const std::string name = pick(csv_path, "bench_components_mklwe.csv");
        std::ofstream comp_csv(name);
        comp_csv << "k,N,n,B,l,log2Q,log2q,component,domain,reps,samples,"
                    "mean,stddev,max_abs,log2_stddev,log2_margin,half_margin,p_fail\n";
        for (const auto& row : param_rows()) {
            if (only_k != 0 && row.k != only_k) continue;
            rc |= run_noise_for_k(row, reps, /*growth=*/false, /*components=*/true,
                                  nullptr, &comp_csv, keygen_csv);
        }
        std::cout << "\n[bench] components CSV -> " << name << "\n";
        break;
    }
    case Mode::Growth: {
        const std::string name = pick(csv_path, "bench_growth_mklwe.csv");
        std::ofstream growth_csv(name);
        growth_csv << "k,N,n,B,l,log2Q,log2q,step,n_added,reps,"
                      "bits_min,bits_mean,bits_max,bits_stddev,budget_bits\n";
        for (const auto& row : param_rows()) {
            if (only_k != 0 && row.k != only_k) continue;
            rc |= run_noise_for_k(row, reps, /*growth=*/true, /*components=*/false,
                                  &growth_csv, nullptr, keygen_csv);
        }
        std::cout << "\n[bench] growth     CSV -> " << name << "\n";
        break;
    }
    case Mode::NoiseAll: {
        // Both noise benches, ONE keygen per k (the fix for the double-keygen).
        const std::string gname = pick(growth_path, "bench_growth_mklwe.csv");
        const std::string cname = pick(noise_path,  "bench_components_mklwe.csv");
        std::ofstream growth_csv(gname);
        growth_csv << "k,N,n,B,l,log2Q,log2q,step,n_added,reps,"
                      "bits_min,bits_mean,bits_max,bits_stddev,budget_bits\n";
        std::ofstream comp_csv(cname);
        comp_csv << "k,N,n,B,l,log2Q,log2q,component,domain,reps,samples,"
                    "mean,stddev,max_abs,log2_stddev,log2_margin,half_margin,p_fail\n";
        for (const auto& row : param_rows()) {
            if (only_k != 0 && row.k != only_k) continue;
            rc |= run_noise_for_k(row, reps, /*growth=*/true, /*components=*/true,
                                  &growth_csv, &comp_csv, keygen_csv);
        }
        std::cout << "\n[bench] growth     CSV -> " << gname << "\n";
        std::cout << "[bench] components CSV -> " << cname << "\n";
        break;
    }
    case Mode::Protocol:
    default: {
        const std::string name = pick(csv_path, "bench_protocol_mklwe.csv");
        std::ofstream protocol_csv(name);
        protocol_csv << "k,reps,min_ms,median_ms,p90_ms,max_ms,mean_ms,stddev_ms,cold_nand_boot_ms\n";
        for (const auto& row : param_rows()) {
            if (only_k != 0 && row.k != only_k) continue;
            rc |= run_protocol(row, reps, &protocol_csv, keygen_csv);
        }
        std::cout << "\n[bench] protocol   CSV -> " << name << "\n";
        break;
    }
    }
    std::cout << "[bench] keygen     CSV -> " << keygen_name << "\n";
    return rc;
}

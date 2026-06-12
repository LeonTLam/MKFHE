// noise-regime-test.cpp — TEMPORARY verification probe (not for commit).
// Measures EvalBinGate(NAND,...) output noise under one key for four input
// regimes that differ only in the magnitude of the input a-vectors:
//   A fresh      : a ~ DGG(sigma=2)        (what --noise measures)
//   B sum256     : a ~ N(0, 2*sqrt(256))   (what --growth measures)
//   C chained    : a from KeySwitch output (realistic gate-to-gate)
//   D uniform-a  : a ~ U(Z_q) synthetic    (theoretical worst case)
// If GetMonomial(0)=0 makes c=0 CMux steps noise-free, output noise must rise
// A < B <= C <= D. Also prints sigma(a) per regime to confirm the mechanism.
#include "binfhecontext.h"
#include "math/discreteuniformgenerator.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace lbcrypto;

static long centered(const NativeInteger& v, long q) {
    long x = static_cast<long>(v.ConvertToInt());
    return x > q / 2 ? x - q : x;
}
static long phase(const MKLWEPrivateKey& sk, const MKLWECiphertext& ct) {
    const NativeInteger q = ct->GetModulus();
    const long qi = static_cast<long>(q.ConvertToInt());
    auto s = sk->GetElement();
    NativeInteger mu = q.ComputeMu(), inner(0);
    for (uint32_t u = 0; u < ct->Getk(); ++u) {
        s[u].SwitchModulus(q);
        for (uint32_t j = 0; j < ct->GetLength(); ++j)
            inner += ct->GetA()[u][j].ModMulFast(s[u][j], q, mu);
        inner.ModEq(q);
    }
    NativeInteger r = ct->GetB();
    r.ModSubFastEq(inner, q);
    return centered(r, qi);
}
static long signed_noise(long v, long qi) {
    const long D = qi / 4;
    long best = v;
    for (long c : {v - D, v + D})
        if (std::labs(c) < std::labs(best)) best = c;
    return best;
}
static double sigma_a(const MKLWECiphertext& ct) {
    const long qi = static_cast<long>(ct->GetModulus().ConvertToInt());
    double ss = 0; std::size_t n = 0;
    for (uint32_t u = 0; u < ct->Getk(); ++u)
        for (uint32_t j = 0; j < ct->GetLength(); ++j) {
            double x = static_cast<double>(centered(ct->GetA()[u][j], qi));
            ss += x * x; ++n;
        }
    return std::sqrt(ss / static_cast<double>(n));
}
static MKLWECiphertext add(const MKLWECiphertext& c1, const MKLWECiphertext& c2) {
    auto a = c1->GetA();
    for (uint32_t u = 0; u < c1->Getk(); ++u) a[u].ModAddEq(c2->GetA()[u]);
    NativeInteger b = c1->GetB();
    b.ModAddFastEq(c2->GetB(), c1->GetModulus());
    return std::make_shared<MKLWECiphertextImpl>(std::move(a), b);
}
struct Stats { double sum=0, sumsq=0, mx=0; long n=0;
    void push(double x){sum+=x;sumsq+=x*x;mx=std::max(mx,std::fabs(x));++n;}
    double mean()const{return sum/n;}
    double sd()const{double m=mean();return std::sqrt(sumsq/n-m*m);}
};

int main(int argc, char** argv) {
    const int R = (argc > 1) ? std::atoi(argv[1]) : 40;
    auto cc = BinFHEContext();
    cc.GenerateBinFHEContext(STD128_MKNTRU_LWE, MKNTRU_LWE);
    auto sk = cc.MKLWE_KeyGen();
    std::cout << "keygen..." << std::flush;
    cc.MKBTKeyGen(sk);
    std::cout << " done. reps=" << R << "\n";
    const long qi = static_cast<long>(cc.GetParams()->GetMKLWEParams()->Getq().ConvertToInt());

    Stats sA, sB, sC, sD; double aA=0, aB=0, aC=0, aD=0;
    DiscreteUniformGeneratorImpl<NativeVector> dug;

    for (int r = 0; r < R; ++r) {
        std::cout << "\r rep " << (r+1) << "/" << R << std::flush;
        // A: fresh-fresh (what --noise measures)
        {
            auto c1 = cc.Encrypt(sk, 1), c2 = cc.Encrypt(sk, 0);
            aA += sigma_a(c1);
            sA.push(signed_noise(phase(sk, cc.EvalBinGate(NAND, c1, c2)), qi));
        }
        // B: 256-fold sum NAND fresh (what --growth measures)
        {
            MKLWECiphertext acc = cc.Encrypt(sk, 0);
            for (int i = 1; i < 256; ++i) acc = add(acc, cc.Encrypt(sk, 0));
            auto one = cc.Encrypt(sk, 1);
            aB += sigma_a(acc);
            sB.push(signed_noise(phase(sk, cc.EvalBinGate(NAND, acc, one)), qi));
        }
        // C: chained — inputs are bootstrap outputs (realistic gate-to-gate)
        {
            auto x = cc.EvalBinGate(NAND, cc.Encrypt(sk, 1), cc.Encrypt(sk, 0)); // ->1
            auto y = cc.EvalBinGate(NAND, cc.Encrypt(sk, 0), cc.Encrypt(sk, 0)); // ->1
            aC += sigma_a(x);
            sC.push(signed_noise(phase(sk, cc.EvalBinGate(NAND, x, y)), qi));    // ->0
        }
        // D: synthetic uniform-a inputs with exact phase m*q/4 (worst case)
        {
            auto mk_uniform = [&](int m) {
                auto ct = cc.Encrypt(sk, m);     // template for shape/modulus
                const NativeInteger q = ct->GetModulus();
                NativeInteger mu = q.ComputeMu();
                auto s = sk->GetElement();
                NativeInteger b(static_cast<uint64_t>(m) * (qi / 4));
                dug.SetModulus(q);
                for (uint32_t u = 0; u < ct->Getk(); ++u) {
                    ct->GetA()[u] = dug.GenerateVector(ct->GetLength());
                    s[u].SwitchModulus(q);
                    for (uint32_t j = 0; j < ct->GetLength(); ++j)
                        b.ModAddFastEq(ct->GetA()[u][j].ModMulFast(s[u][j], q, mu), q);
                }
                ct->SetB(b);
                return ct;
            };
            auto c1 = mk_uniform(1), c2 = mk_uniform(0);
            aD += sigma_a(c1);
            sD.push(signed_noise(phase(sk, cc.EvalBinGate(NAND, c1, c2)), qi));
        }
    }
    std::cout << "\n\nregime        sigma(a_in)   out mean   out sigma   out max|e|   E~bits\n";
    auto row = [&](const char* nm, double a, const Stats& s) {
        std::cout << nm << "  " << std::fixed;
        std::cout.width(10); std::cout.precision(1); std::cout << a / R;
        std::cout.width(11); std::cout << s.mean();
        std::cout.width(12); std::cout << s.sd();
        std::cout.width(12); std::cout << s.mx;
        std::cout.width(9);  std::cout.precision(2);
        std::cout << (s.sd() > 0 ? std::log2(s.sd()) : 0) << "\n";
    };
    row("A fresh    ", aA, sA);
    row("B sum256   ", aB, sB);
    row("C chained  ", aC, sC);
    row("D uniform-a", aD, sD);
    std::cout << "\nbudget q/8 = " << qi / 8 << "  (failure if |e| above)\n";
    return 0;
}

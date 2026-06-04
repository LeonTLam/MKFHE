//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

// binfhe-timing.h
//
// Lightweight per-component timing instrumentation for the MKNTRU_LWE
// (Algorithm 2) bootstrapping path. It lets boolean-mklwe-bench.cpp --timing
// attribute a gate's wall-clock to the four internal primitives that the public
// BinFHEContext API otherwise hides:
//
//   CMux+ExtProduct  -> UniEncAccumulatorXZW_B::EvalAcc loop body (n*k iters)
//   ExternalProduct  -> UniEncAccumulatorXZW_B::HbProd (the hybrid product)
//   KS_NGS_to_LWE    -> MKLWEscheme->KeySwitch in the MKNTRU_LWE EvalBinGate
//   ModSwitch        -> MKLWEscheme->ModSwitch (Q->qKS) in that same branch
//
// IMPORTANT (one-definition rule across DSOs): the accumulators g_timing_profile
// MUST live in exactly ONE translation unit (binfhe-timing.cpp). If Reset/Get
// were inline in this header, every .so / executable that includes it would get
// its own copy and the bench would read zero. Hence ResetTimingProfile() and
// GetTimingProfile() are declared here but DEFINED only in binfhe-timing.cpp,
// and g_timing_profile is `extern` here. Verify after build with:
//   nm -D lib/libOPENFHEbinfhe.so | grep g_timing_profile   # exactly one 'B'

#ifndef _BINFHE_TIMING_H_
#define _BINFHE_TIMING_H_

#include <chrono>
#include <cstdint>

namespace lbcrypto {

// Accumulated nanoseconds + call counts for each instrumented primitive.
// All fields are summed across one or more EvalBinGate calls; the bench resets
// them via ResetTimingProfile() before a measured batch and divides by the
// observed call count to get a per-call number.
struct TimingProfile {
    // total nanoseconds spent inside each primitive
    uint64_t cmux_ns = 0;   // EvalAcc loop body (CMux + fused external product)
    uint64_t ext_ns  = 0;   // HbProd (external/hybrid product) only
    uint64_t ks_ns   = 0;   // MKLWE KeySwitch (NGS->LWE key switch)
    uint64_t ms_ns   = 0;   // MKLWE ModSwitch (Q->qKS)

    // number of times each primitive ran
    uint64_t cmux_count = 0;
    uint64_t ext_count  = 0;
    uint64_t ks_count   = 0;
    uint64_t ms_count   = 0;
};

// The single shared profile. DEFINED ONCE in binfhe-timing.cpp.
// thread_local so the (single-threaded under WITH_OPENMP=OFF) bootstrap path
// accumulates without atomics; the bench reads it on the same thread.
extern thread_local TimingProfile g_timing_profile;

// Zero the profile / read a copy of it. DEFINED ONLY in binfhe-timing.cpp so
// there is one canonical accessor regardless of how many DSOs link this header.
void          ResetTimingProfile();
TimingProfile GetTimingProfile();

// RAII timer: on destruction it adds the elapsed nanoseconds into *ns and, if
// provided, bumps *count. Header-inline is fine here because it only ever
// touches the single g_timing_profile via the pointers it is handed -- it holds
// no static state of its own, so there is no per-DSO duplication to worry about.
class ScopedTimer {
public:
    explicit ScopedTimer(uint64_t* ns, uint64_t* count = nullptr)
        : m_start(std::chrono::steady_clock::now()), m_ns(ns), m_count(count) {}

    ~ScopedTimer() {
        const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                 std::chrono::steady_clock::now() - m_start)
                                 .count();
        if (m_ns)
            *m_ns += static_cast<uint64_t>(elapsed);
        if (m_count)
            ++(*m_count);
    }

    ScopedTimer(const ScopedTimer&)            = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;

private:
    std::chrono::steady_clock::time_point m_start;
    uint64_t*                             m_ns;
    uint64_t*                             m_count;
};

}  // namespace lbcrypto

#endif  // _BINFHE_TIMING_H_

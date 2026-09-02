#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <unordered_set>

#include <gtest/gtest.h>

#include "ERF_GTestHashRngCommon.H"

using namespace hash_rng_test;

namespace {

struct SplitMixKnownAnswer
{
    std::uint64_t input;
    std::uint64_t expected;
};

struct CellKeyKnownAnswer
{
    int i;
    int j;
    int k;
    int comp;
    int lev;
    std::uint64_t seed;
    std::uint64_t expected;
};

} // namespace

// Known-answer values use Sebastiano Vigna's public-domain SplitMix64
// reference implementation (https://prng.di.unimi.it/splitmix64.c). The cell
// vectors were obtained by applying that reference finalizer to ERF's documented
// salt sequence with independent unsigned-64-bit wraparound calculations.
TEST(HashRngScalar, KnownAnswerVectors)
{
    constexpr std::array<SplitMixKnownAnswer, 4> splitmix_answers{{
        {UINT64_C(0x0000000000000000), UINT64_C(0xe220a8397b1dcdaf)},
        {UINT64_C(0x0000000000000001), UINT64_C(0x910a2dec89025cc1)},
        {UINT64_C(0xffffffffffffffff), UINT64_C(0xe4d971771b652c20)},
        {UINT64_C(0x0123456789abcdef), UINT64_C(0x157a3807a48faa9d)}
    }};
    for (const auto& answer : splitmix_answers) {
        EXPECT_EQ(erf_hash_rng::splitmix64(answer.input), answer.expected);
    }

    constexpr std::array<CellKeyKnownAnswer, 4> cell_key_answers{{
        {0, 0, 0, 0, 0, UINT64_C(0x0000000000000000), UINT64_C(0xd53226341f0640f4)},
        {1, 2, 3, 1, 0, UINT64_C(0x0123456789abcdef), UINT64_C(0xece8a519c04916b9)},
        {-7, 19, -3, 3, 2, UINT64_C(0xfedcba9876543210), UINT64_C(0xbfe85453268a6b9e)},
        {22, 18, 4, 2, 1, UINT64_C(0x0000000000000400), UINT64_C(0xaa28445367478916)}
    }};
    for (const auto& answer : cell_key_answers) {
        EXPECT_EQ(erf_hash_rng::cell_key(
                      answer.i, answer.j, answer.k, answer.comp, answer.lev, answer.seed),
                  answer.expected);
    }
}

TEST(HashRngScalar, UniformMomentsAndRanges)
{
    constexpr int sample_count = 200000;
    long double sum = 0.0L;
    long double sum_of_squares = 0.0L;
    amrex::Real uniform_min = std::numeric_limits<amrex::Real>::max();
    amrex::Real uniform_max = std::numeric_limits<amrex::Real>::lowest();
    amrex::Real symmetric_min = std::numeric_limits<amrex::Real>::max();
    amrex::Real symmetric_max = std::numeric_limits<amrex::Real>::lowest();

    for (int sample = 0; sample < sample_count; ++sample) {
        const int i = sample - sample_count / 2;
        const int j = (sample * 17) % 1009 - 504;
        const int k = (sample * 31) % 127 - 63;
        const int comp = sample % 4;
        const int lev = (sample / 4) % 3;
        const amrex::Real uniform =
            erf_hash_rng::hash_uniform(i, j, k, comp, lev, kTestSeed);
        const amrex::Real symmetric =
            erf_hash_rng::hash_symmetric(i, j, k, comp, lev, kTestSeed);
        sum += static_cast<long double>(uniform);
        sum_of_squares += static_cast<long double>(uniform) *
                          static_cast<long double>(uniform);
        uniform_min = std::min(uniform_min, uniform);
        uniform_max = std::max(uniform_max, uniform);
        symmetric_min = std::min(symmetric_min, symmetric);
        symmetric_max = std::max(symmetric_max, symmetric);
    }

    const long double mean = sum / sample_count;
    const long double variance = sum_of_squares / sample_count - mean * mean;
    EXPECT_GE(uniform_min, amrex::Real(0.0));
    EXPECT_LT(uniform_max, amrex::Real(1.0));
    EXPECT_GE(symmetric_min, amrex::Real(-1.0));
    EXPECT_LT(symmetric_max, amrex::Real(1.0));
    EXPECT_NEAR(static_cast<double>(mean), 0.5, 0.002);
    EXPECT_NEAR(static_cast<double>(variance), 1.0 / 12.0, 0.001);
}

TEST(HashRngScalar, SplitMixAvalancheFlipsHalfTheOutputBits)
{
    constexpr int sample_count = 1024;
    std::array<long, 64> flipped_by_input_bit{};
    long total_flipped = 0;

    for (int sample = 0; sample < sample_count; ++sample) {
        const std::uint64_t key = erf_hash_rng::cell_key(
            sample, sample * 7 - 3000, 2 * sample + 1, sample % 4,
            sample % 3, kTestSeed);
        const std::uint64_t baseline = erf_hash_rng::splitmix64(key);
        for (int bit = 0; bit < 64; ++bit) {
            const std::uint64_t changed =
                erf_hash_rng::splitmix64(key ^ (UINT64_C(1) << bit));
            const int flipped = std::popcount(baseline ^ changed);
            flipped_by_input_bit[bit] += flipped;
            total_flipped += flipped;
        }
    }

    const double overall_mean =
        static_cast<double>(total_flipped) / static_cast<double>(sample_count * 64);
    EXPECT_GT(overall_mean, 31.5);
    EXPECT_LT(overall_mean, 32.5);
    for (int bit = 0; bit < 64; ++bit) {
        const double bit_mean =
            static_cast<double>(flipped_by_input_bit[bit]) / sample_count;
        EXPECT_GT(bit_mean, 29.0) << "input_bit=" << bit;
        EXPECT_LT(bit_mean, 35.0) << "input_bit=" << bit;
    }
}

TEST(HashRngScalar, SweptCellKeysHaveNoCollisions)
{
    std::unordered_set<std::uint64_t> keys;
    keys.reserve(50000);
    for (int k = -3; k <= 7; ++k) {
        for (int j = -7; j <= 21; ++j) {
            for (int i = -8; i <= 23; ++i) {
                for (int comp = 0; comp < 4; ++comp) {
                    const std::uint64_t key =
                        erf_hash_rng::cell_key(i, j, k, comp, 2, kTestSeed);
                    const auto insertion = keys.insert(key);
                    if (!insertion.second) {
                        ADD_FAILURE() << "collision at i=" << i
                                      << " j=" << j
                                      << " k=" << k
                                      << " comp=" << comp;
                        return;
                    }
                }
            }
        }
    }
    EXPECT_EQ(keys.size(), std::size_t(32 * 29 * 11 * 4));
}

TEST(HashRngScalar, StaggeredVelocityComponentSaltsAreDistinct)
{
    for (int sample = 0; sample < 2048; ++sample) {
        const int i = sample % 47 - 23;
        const int j = (sample * 11) % 43 - 21;
        const int k = (sample * 19) % 17 - 8;
        const std::uint64_t u = erf_hash_rng::cell_key(
            i, j, k, kXVelocityComp, 0, kTestSeed);
        const std::uint64_t v = erf_hash_rng::cell_key(
            i, j, k, kYVelocityComp, 0, kTestSeed);
        const std::uint64_t w = erf_hash_rng::cell_key(
            i, j, k, kZVelocityComp, 0, kTestSeed);
        if (u == v || u == w || v == w) {
            ADD_FAILURE() << "component-key collision at i=" << i
                          << " j=" << j
                          << " k=" << k;
            return;
        }
    }
}

TEST(HashRngScalar, SharedCaseEntryPointIsDeterministic)
{
    for (int case_id = 0; case_id < kHostDeviceCaseCount; ++case_id) {
        const std::uint64_t first = hash_rng_case_word(case_id);
        if (first != hash_rng_case_word(case_id)) {
            ADD_FAILURE() << "case_id=" << case_id;
            return;
        }
    }
}

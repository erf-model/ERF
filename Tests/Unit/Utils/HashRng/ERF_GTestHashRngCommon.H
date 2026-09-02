#ifndef ERF_GTEST_HASH_RNG_COMMON_H_
#define ERF_GTEST_HASH_RNG_COMMON_H_

#include <cstdint>

#include <ERF_HashRNG.H>

namespace hash_rng_test {

constexpr std::uint64_t kTestSeed = UINT64_C(0x0123456789abcdef);
constexpr int kTemperatureComp = 0;
constexpr int kXVelocityComp = 1;
constexpr int kYVelocityComp = 2;
constexpr int kZVelocityComp = 3;
constexpr int kHostDeviceCaseCount = 4096;
constexpr int kRealMantissaBits = (sizeof(amrex::Real) == sizeof(float)) ? 24 : 53;
constexpr std::uint64_t kRealIntegerScale = UINT64_C(1) << kRealMantissaBits;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::uint64_t hash_rng_case_word (const int case_id) noexcept
{
    const int i = (case_id * 37) % 509 - 254;
    const int j = (case_id * 53) % 401 - 200;
    const int k = (case_id * 29) % 97 - 48;
    // comp and lev must not be derived from case_id % 4, which selects the entry
    // point in the switch below: sharing that expression would pin each entry
    // point to a single component salt, so a host/device divergence confined to
    // one (entry point, salt) pair would never be sampled.
    const int comp = (case_id / 4) % 4;
    const int lev = (case_id / 16) % 5;
    const std::uint64_t seed = kTestSeed ^ static_cast<std::uint64_t>(case_id / 20);

    switch (case_id % 4) {
    case 0:
        return erf_hash_rng::cell_key(i, j, k, comp, lev, seed);
    case 1:
        return erf_hash_rng::splitmix64(
            erf_hash_rng::cell_key(i, j, k, comp, lev, seed));
    case 2:
        return static_cast<std::uint64_t>(
            erf_hash_rng::hash_uniform(i, j, k, comp, lev, seed) *
            static_cast<amrex::Real>(kRealIntegerScale));
    default:
        return static_cast<std::uint64_t>(
            (erf_hash_rng::hash_symmetric(i, j, k, comp, lev, seed) +
             amrex::Real(1.0)) *
            static_cast<amrex::Real>(kRealIntegerScale / UINT64_C(2)));
    }
}

} // namespace hash_rng_test

#endif

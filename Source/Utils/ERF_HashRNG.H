#ifndef ERF_HASH_RNG_H_
#define ERF_HASH_RNG_H_

#include <cstdint>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace erf_hash_rng {

/**
 * Apply the standard SplitMix64 integer avalanche to an input word.
 *
 * @param[in]x input word
 * @return   mixed 64-bit word
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::uint64_t
splitmix64 (std::uint64_t x) noexcept
{
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27U)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31U);
}

/**
 * Construct a decomposition-independent key for one indexed component.
 *
 * @param[in]i    global x index
 * @param[in]j    global y index
 * @param[in]k    global z index
 * @param[in]comp component salt
 * @param[in]lev  AMR level
 * @param[in]seed user seed
 * @return      mixed 64-bit cell key
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::uint64_t
cell_key (int i, int j, int k, int comp, int lev, std::uint64_t seed) noexcept
{
    std::uint64_t key = splitmix64(seed ^ UINT64_C(0x243f6a8885a308d3));
    key = splitmix64(key ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(i)) +
                            UINT64_C(0x13198a2e03707344)));
    key = splitmix64(key ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(j)) +
                            UINT64_C(0xa4093822299f31d0)));
    key = splitmix64(key ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(k)) +
                            UINT64_C(0x082efa98ec4e6c89)));
    key = splitmix64(key ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(comp)) +
                            UINT64_C(0x452821e638d01377)));
    return splitmix64(key ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(lev)) +
                             UINT64_C(0xbe5466cf34e90c6c)));
}

/**
 * Map an indexed cell key to a uniformly distributed value in [0, 1).
 *
 * @param[in]i    global x index
 * @param[in]j    global y index
 * @param[in]k    global z index
 * @param[in]comp component salt
 * @param[in]lev  AMR level
 * @param[in]seed user seed
 * @return      uniform value in [0, 1)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
hash_uniform (int i, int j, int k, int comp, int lev, std::uint64_t seed) noexcept
{
    constexpr int mantissa_bits = (sizeof(amrex::Real) == sizeof(float)) ? 24 : 53;
    constexpr std::uint64_t denominator = UINT64_C(1) << mantissa_bits;
    const std::uint64_t numerator = cell_key(i, j, k, comp, lev, seed) >>
                                    (64 - mantissa_bits);
    return static_cast<amrex::Real>(numerator) /
           static_cast<amrex::Real>(denominator);
}

/**
 * Map an indexed cell key to a uniformly distributed value in [-1, 1).
 *
 * @param[in]i    global x index
 * @param[in]j    global y index
 * @param[in]k    global z index
 * @param[in]comp component salt
 * @param[in]lev  AMR level
 * @param[in]seed user seed
 * @return      symmetric value in [-1, 1)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
hash_symmetric (int i, int j, int k, int comp, int lev, std::uint64_t seed) noexcept
{
    return amrex::Real(2.0) * hash_uniform(i, j, k, comp, lev, seed) -
           amrex::Real(1.0);
}

} // namespace erf_hash_rng

#endif

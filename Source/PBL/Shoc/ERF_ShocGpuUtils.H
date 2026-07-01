#ifndef ERF_SHOC_GPU_UTILS_H_
#define ERF_SHOC_GPU_UTILS_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

namespace shoc {

enum class InitRunOn {
    Host,
    Device
};

// SHOC GPU portability rules:
// 1. ShocColumnData, FArrayBox, and MultiFab are host-side ownership objects.
// 2. Device-callable functions must take Array4 views, POD structs, scalars, and options only.
// 3. Device-callable functions must not call .array(), .const_array(), .box(), .nComp(), resize(), or setVal().
// 4. Do not capture this in AMReX GPU lambdas.
// 5. Do not reference anonymous-namespace constexpr Real constants inside GPU lambdas.
// 6. Use SHOC constants through AMREX_GPU_HOST_DEVICE inline accessors.
// 7. Use amrex::min, amrex::max, amrex::Math::abs, and explicit arithmetic in device paths.
// 8. Use a central SHOC FArrayBox initialization helper, not direct FArrayBox::setVal calls.
// 9. Unit tests must allocate column data in backend-safe memory and synchronize before host assertions.

inline InitRunOn
default_init_run_on ()
{
#ifdef AMREX_USE_GPU
    return InitRunOn::Device;
#else
    return InitRunOn::Host;
#endif
}

inline void
set_fab_val (amrex::FArrayBox& fab,
             amrex::Real value,
             InitRunOn run_on)
{
    if (run_on == InitRunOn::Device) {
#ifdef AMREX_USE_GPU
        fab.setVal<amrex::RunOn::Device>(value, fab.box(), 0, fab.nComp());
#else
        fab.setVal<amrex::RunOn::Host>(value, fab.box(), 0, fab.nComp());
#endif
    } else {
        fab.setVal<amrex::RunOn::Host>(value, fab.box(), 0, fab.nComp());
    }
}

inline void
set_fab_val (amrex::FArrayBox& fab, amrex::Real value)
{
    set_fab_val(fab, value, default_init_run_on());
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
height_agl (amrex::Real z, amrex::Real z_sfc) noexcept
{
    return amrex::max(amrex::Real(0.0), z - z_sfc);
}

inline void
sync_if_needed ()
{
    // Test-owned helper for host-side assertions. Production SHOC code should
    // keep synchronization out of the runtime path and use this only through
    // shoc_test::run_and_sync or other explicit test utilities.
#ifdef AMREX_USE_GPU
    amrex::Gpu::streamSynchronize();
#endif
}

} // namespace shoc

#endif

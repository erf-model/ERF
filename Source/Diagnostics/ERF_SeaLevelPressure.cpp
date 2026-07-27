#include "Diagnostics/ERF_SeaLevelPressure.H"

#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>

#include "ERF_IndexDefines.H"
#include "ERF_TerrainMetrics.H"

using namespace amrex;

namespace sea_level_pressure_diagnostics
{

void fill_sea_level_pressure (MultiFab& destination,
                              int destination_component,
                              const MultiFab& conserved_state,
                              const MultiFab& physical_nodal_height,
                              int vapor_component,
                              int klo)
{
    AMREX_ALWAYS_ASSERT(destination_component >= 0);
    AMREX_ALWAYS_ASSERT(destination_component < destination.nComp());
    AMREX_ALWAYS_ASSERT(conserved_state.nComp() > Rho_comp);
    AMREX_ALWAYS_ASSERT(conserved_state.nComp() > RhoTheta_comp);
    AMREX_ALWAYS_ASSERT(vapor_component < 0 || vapor_component < conserved_state.nComp());

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(destination, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const auto destination_arr = destination.array(mfi);
        const auto conserved_arr = conserved_state.const_array(mfi);
        const auto height_arr = physical_nodal_height.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real dry_density = conserved_arr(i, j, klo, Rho_comp);
            const Real dry_density_times_potential_temperature =
                conserved_arr(i, j, klo, RhoTheta_comp);
            const Real vapor_mixing_ratio = vapor_component < 0
                ? zero
                : conserved_arr(i, j, klo, vapor_component) / dry_density;
            const Real height_lowest_level =
                Compute_Z_AtCellCenter(i, j, klo, height_arr);
            const Real height_surface = Compute_Z_AtWFace(i, j, klo, height_arr);

            destination_arr(i, j, k, destination_component) =
                compute_from_erf_lowest_level_state(
                    dry_density, dry_density_times_potential_temperature,
                    vapor_mixing_ratio, height_lowest_level, height_surface);
        });
    }
}

} // namespace sea_level_pressure_diagnostics

#include <ERF_FireSmokeEmission.H>

#ifdef ERF_ENABLE_FIRE

#include <ERF_IndexDefines.H>
#include <AMReX_MFIter.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Print.H>

using namespace amrex;

void inject_smoke_from_fire(
    MultiFab&       cc_source,
    const MultiFab& fire_heat_atm,
    const MultiFab& z_phys_cc,
    const Geometry& geom_atm,
    Real            emission_factor,
    Real            heat_of_combustion,
    int             smoke_comp,
    bool            fire_debug,
    int             step)
{
    if (heat_of_combustion <= 0.0_rt) return;
    if (emission_factor    <= 0.0_rt) return;

    const auto& dx     = geom_atm.CellSizeArray();
    const auto& domain = geom_atm.Domain();
    const int   klo    = domain.smallEnd(2);
    const int   khi    = domain.bigEnd(2);

    for (MFIter mfi(cc_source, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.tilebox();
        // Only inject at k=klo
        const Box& bx0 = makeSlab(bx, 2, klo);
        if (bx0.isEmpty()) continue;

        auto src  = cc_source.array(mfi, smoke_comp);
        auto heat = fire_heat_atm.const_array(mfi);
        auto zcc  = z_phys_cc.const_array(mfi);

        const Real ef  = emission_factor;
        const Real hoc = heat_of_combustion;
        const Real dz  = dx[2];

        ParallelFor(bx0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Cell height — use terrain-aware dz if available
            Real dz_cell = (k < khi)
                ? amrex::max(zcc(i,j,k+1) - zcc(i,j,k), dz * 0.1_rt)
                : dz;
            // smoke_flux [kg/m2/s] = ef * Q [W/m2] / hoc [J/kg]
            Real smoke_flux = ef * amrex::max(heat(i,j,k), 0.0_rt) / hoc;
            // inject as volumetric source [kg/m3/s] = flux / dz
            src(i, j, k) += smoke_flux / dz_cell;
        });
    }

    if (fire_debug) {
        Real src_max  = cc_source.max(smoke_comp);
        Real heat_max = fire_heat_atm.max(0);
        Print() << "[FIRE DEBUG] Phase 4 smoke: step=" << step
                << " fire_heat_atm_max=" << heat_max << " W/m2"
                << " smoke_src_max=" << src_max << " kg/m3/s\n";
    }
}

#endif // ERF_ENABLE_FIRE

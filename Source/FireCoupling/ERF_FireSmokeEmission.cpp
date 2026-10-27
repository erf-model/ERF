#include <ERF_FireSmokeEmission.H>

#ifdef ERF_USE_FIRE

#include <ERF_IndexDefines.H>
#include <AMReX_MFIter.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Print.H>

using namespace amrex;

void inject_smoke_from_fire(
    MultiFab&       cc_source,
    const MultiFab& fire_heat_atm,
    const MultiFab* z_phys_nd,
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
        // Only inject at k = klo, and only from the tile that holds it:
        // makeSlab of a tile above klo is a non-empty box outside the tile,
        // so every z tile of the column added the source once more.
        if (bx.smallEnd(2) > klo || bx.bigEnd(2) < klo) continue;
        const Box bx0 = makeSlab(bx, 2, klo);

        auto src  = cc_source.array(mfi, smoke_comp);
        auto heat = fire_heat_atm.const_array(mfi);
        const bool has_nd = (z_phys_nd != nullptr);
        Array4<const Real> znd;
        if (has_nd) { znd = z_phys_nd->const_array(mfi); }

        const Real ef  = emission_factor;
        const Real hoc = heat_of_combustion;
        const Real dz  = dx[2];

        ParallelFor(bx0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Thickness of the k = 0 cell: the mean of its four node pairs
            // where the nodal heights exist (ERF always passes them), the
            // uniform dz otherwise. The injected mass per column, src dz_cell,
            // then equals the flux on any vertical spacing; only a degenerate
            // cell is floored, as in fire_cell_thickness.
            Real dz_cell = dz;
            if (has_nd) {
                dz_cell = 0.25_rt * ( (znd(i,  j,  k+1) - znd(i,  j,  k))
                                    + (znd(i+1,j,  k+1) - znd(i+1,j,  k))
                                    + (znd(i,  j+1,k+1) - znd(i,  j+1,k))
                                    + (znd(i+1,j+1,k+1) - znd(i+1,j+1,k)) );
                dz_cell = amrex::max(dz_cell, 1.0e-6_rt * dz);
            }
            amrex::ignore_unused(khi);
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
                << " smoke_src_max=" << src_max << " kg/m3/s"
                << " heat_per_kg=" << heat_of_combustion << " J/kg\n";
    }
}

#endif // ERF_USE_FIRE

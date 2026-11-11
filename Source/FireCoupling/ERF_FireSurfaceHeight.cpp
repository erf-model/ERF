#include <ERF_FireSurfaceHeight.H>

using namespace amrex;

void compute_fire_surface_height(
    MultiFab& fire_surface_z,
    const MultiFab* z_phys_nd,
    const Geometry& geom_atm,
    const FireGrid& fg)
{
    // Without nodal heights (callers other than ERF, which always passes
    // them): a flat ground at the floor of the domain.
    if (z_phys_nd == nullptr) {
        fire_surface_z.setVal(geom_atm.ProbLo(2));
        return;
    }

    const int C = fg.C;
    // Atmospheric index of the fire region's lower corner (FireGrid::atm_lo)
    const int atm_nlo_x = fg.atm_lo[0];
    const int atm_nlo_y = fg.atm_lo[1];

    for (MFIter mfi(fire_surface_z, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> zs = fire_surface_z.array(mfi);
        Array4<const Real> z_nd = z_phys_nd->array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            const int i_f = iv[0];
            const int j_f = iv[1];

            // Atmospheric cell holding this fire cell, and its four k = 0 nodes
            const int i_n_lo = i_f / C + atm_nlo_x;
            const int j_n_lo = j_f / C + atm_nlo_y;

            zs(i_f, j_f, 0) = Real(0.25) * (z_nd(i_n_lo,     j_n_lo,     0)
                                    + z_nd(i_n_lo + 1, j_n_lo,     0)
                                    + z_nd(i_n_lo,     j_n_lo + 1, 0)
                                    + z_nd(i_n_lo + 1, j_n_lo + 1, 0));
        });
    }
}

void compute_fire_column_grounds(
    MultiFab& fire_col_ground,
    const MultiFab* z_phys_nd,
    const Geometry& geom_atm,
    const FireGrid& fg)
{
    if (z_phys_nd == nullptr) {
        fire_col_ground.setVal(geom_atm.ProbLo(2));
        return;
    }

    const int C = fg.C;
    // Atmospheric index of the fire region's lower corner (FireGrid::atm_lo)
    const int atm_nlo_x = fg.atm_lo[0];
    const int atm_nlo_y = fg.atm_lo[1];

    for (MFIter mfi(fire_col_ground, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> zg = fire_col_ground.array(mfi);
        Array4<const Real> z_nd = z_phys_nd->array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            const int i_f = iv[0];
            const int j_f = iv[1];

            // Lower-left column of the bilinear stencil, in atmospheric indices
            const Real gx = (Real(i_f) + Real(0.5)) / Real(C) - Real(0.5);
            const Real gy = (Real(j_f) + Real(0.5)) / Real(C) - Real(0.5);
            const int i0 = static_cast<int>(std::floor(gx));
            const int j0 = static_cast<int>(std::floor(gy));

            for (int c = 0; c < 4; ++c) {
                const int ia = i0 + (c & 1);
                const int ja = j0 + ((c >> 1) & 1);

                // Cell (ia,ja) sits between nodes (ia,ja) and (ia+1,ja+1)
                const int i_n = ia + atm_nlo_x;
                const int j_n = ja + atm_nlo_y;

                zg(i_f, j_f, 0, c) = Real(0.25) * (z_nd(i_n,     j_n,     0)
                                           + z_nd(i_n + 1, j_n,     0)
                                           + z_nd(i_n,     j_n + 1, 0)
                                           + z_nd(i_n + 1, j_n + 1, 0));
            }
        });
    }
}

#include <ERF_TerrainSlope.H>
#include <AMReX_LoopConcurrent.H>

using namespace amrex;

void compute_terrain_slopes(
    MultiFab& fire_slopes,
    const MultiFab& z_phys_nd,
    const Geometry& geom_atm,
    const FireGrid& fg)
{
    // Compute dz/dx and dz/dy on fire grid by interpolating from nodal z_phys_nd
    // Fire grid is C times finer than atmospheric grid

    int C = fg.C;
    Real dx_atm = geom_atm.CellSize(0);
    Real dy_atm = geom_atm.CellSize(1);

    // Atmospheric domain
    const Box& domain_atm = geom_atm.Domain();
    int atm_nlo_x = domain_atm.smallEnd(0);
    int atm_nlo_y = domain_atm.smallEnd(1);

    for (MFIter mfi(fire_slopes, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> slopes = fire_slopes.array(mfi);
        Array4<const Real> z_nd = z_phys_nd.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            int i_f = iv[0];
            int j_f = iv[1];

            // Map fire cell to atmospheric cell indices
            // Fire cell (i_f, j_f) is in atmospheric cell (i_a, j_a)
            int i_a = i_f / C;
            int j_a = j_f / C;

            // Nodal indices: z_phys_nd is nodal (one extra in each direction)
            // Node (i_n, j_n) is at corner of cell (i_n-1, j_n-1)
            int i_n_lo = i_a + atm_nlo_x;
            int j_n_lo = j_a + atm_nlo_y;

            // Get the four surrounding nodal values for bilinear interp
            Real z_ll = z_nd(i_n_lo,   j_n_lo,   0);  // (i,   j)
            Real z_lr = z_nd(i_n_lo+1, j_n_lo,   0);  // (i+1, j)
            Real z_ul = z_nd(i_n_lo,   j_n_lo+1, 0);  // (i,   j+1)
            Real z_ur = z_nd(i_n_lo+1, j_n_lo+1, 0);  // (i+1, j+1)

            // Centered finite differences for slopes
            // dz/dx at cell center
            Real dz_dx_left = (z_ul - z_ll) / (2.0 * dy_atm);
            Real dz_dx_right = (z_ur - z_lr) / (2.0 * dy_atm);
            slopes(i_f, j_f, 0, 0) = (z_lr - z_ll) / dx_atm;

            // dz/dy at cell center
            Real dz_dy_bot = (z_lr - z_ll) / (2.0 * dx_atm);
            Real dz_dy_top = (z_ur - z_ul) / (2.0 * dx_atm);
            slopes(i_f, j_f, 0, 1) = (z_ul - z_ll) / dy_atm;
        });
    }
}

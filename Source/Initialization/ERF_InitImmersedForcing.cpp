/**
 * \file ERF_InitImmersedForcing.cpp
 */

#include <ERF.H>

using namespace amrex;

/**
 * Set velocities in cells that are immersed to be 0 (or a very small number)
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_immersed_forcing (int lev)
{
    auto& lev_new = vars_new[lev];
    MultiFab* terrain_blank = terrain_blanking[lev].get();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(lev_new[Vars::cons], TileNoZ()); mfi.isValid(); ++mfi)
    {
        const Box &xbx = mfi.tilebox(IntVect(1,0,0));
        const Box &ybx = mfi.tilebox(IntVect(0,1,0));
        const Box &zbx = mfi.tilebox(IntVect(0,0,1));
        const Real epsilon = 1e-2;

        const Array4<const Real>& t_blank_arr = terrain_blank->const_array(mfi);

        const auto &xvel_arr = lev_new[Vars::xvel].array(mfi);
        const auto &yvel_arr = lev_new[Vars::yvel].array(mfi);
        const auto &zvel_arr = lev_new[Vars::zvel].array(mfi);

        // Set the x,y,z-velocities
        ParallelFor(xbx, ybx, zbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real t_blank = 0.5 * (t_blank_arr(i, j, k) + t_blank_arr(i-1, j, k));
            if (t_blank == 1.0) { xvel_arr(i, j, k) = epsilon; }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real t_blank = 0.5 * (t_blank_arr(i, j, k) + t_blank_arr(i, j-1, k));
            if (t_blank == 1.0) { yvel_arr(i, j, k) = epsilon; }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real t_blank = 0.5 * (t_blank_arr(i, j, k) + t_blank_arr(i, j, k-1));
            if (t_blank == 1.0) { zvel_arr(i, j, k) = epsilon; }
        });
    } //mfi
}
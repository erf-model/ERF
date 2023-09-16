/**
 * \file VelocityToMomentum.cpp
 */
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <Utils.H>

using namespace amrex;

/**
 * Convert velocity to momentum by multiplying by density averaged onto faces.
 * @param[in] xvel_in x-component of velocity
 * @param[in] xvel_ngrow how many cells to grow the tilebox for the x-momentum
 * @param[in] yvel_in y-component of velocity
 * @param[in] yvel_ngrow how many cells to grow the tilebox for the y-momentum
 * @param[in] zvel_in z-component of velocity
 * @param[in] zvel_ngrow how many cells to grow the tilebox for the z-momentum
 * @param[in] density density at cell centers
 * @param[out] xmom x-component of momentum
 * @param[out] ymom y-component of momentum
 * @param[out] zmom z-component of momentum
 * @param[in] l_use_ndiff flag describing whether we will later add explicit numerical diffusion
 */

void VelocityToMomentum( const MultiFab& xvel_in,
                         const IntVect& xvel_ngrow,
                         const MultiFab& yvel_in,
                         const IntVect& yvel_ngrow,
                         const MultiFab& zvel_in,
                         const IntVect& zvel_ngrow,
                         const MultiFab& density,
                         MultiFab& xmom, MultiFab& ymom, MultiFab& zmom,
                         bool l_use_ndiff)
{
    BL_PROFILE_VAR("VelocityToMomentum()",VelocityToMomentum);

    // Loop over boxes = valid boxes grown by ngrow
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // We need momentum in the interior ghost cells (init == real)
        const Box& valid_box = amrex::grow(mfi.validbox(),IntVect(1,1,0));
        Box tbx, tby, tbz;
        if (l_use_ndiff) {
            tbx = amrex::grow(mfi.nodaltilebox(0),xvel_ngrow);
            if (tbx.smallEnd(2) < 0)  tbx.setSmall(2,0);
            tby = amrex::grow(mfi.nodaltilebox(1),yvel_ngrow);
            if (tby.smallEnd(2) < 0)  tby.setSmall(2,0);
            tbz = amrex::grow(mfi.nodaltilebox(2),zvel_ngrow);
            if (tbz.smallEnd(2) < 0)  tbz.setSmall(2,0);
        } else {
            tbx = amrex::grow(mfi.nodaltilebox(0),xvel_ngrow) & surroundingNodes(valid_box,0);
            if (tbx.smallEnd(2) < 0)  tbx.setSmall(2,0);
            tby = amrex::grow(mfi.nodaltilebox(1),yvel_ngrow) & surroundingNodes(valid_box,1);
            if (tby.smallEnd(2) < 0)  tby.setSmall(2,0);
            tbz = amrex::grow(mfi.nodaltilebox(2),zvel_ngrow) & surroundingNodes(valid_box,2);
            if (tbz.smallEnd(2) < 0)  tbz.setSmall(2,0);
        }

        // Conserved/state variables on cell centers -- we use this for density
        const Array4<const Real>& dens_arr = density.array(mfi);

        // Momentum on faces, to be computed
        Array4<Real> const& momx = xmom.array(mfi);
        Array4<Real> const& momy = ymom.array(mfi);
        Array4<Real> const& momz = zmom.array(mfi);

        // Velocity on faces, used in computation
        const Array4<Real const>& velx = xvel_in.array(mfi);
        const Array4<Real const>& vely = yvel_in.array(mfi);
        const Array4<Real const>& velz = zvel_in.array(mfi);

        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momx(i,j,k) = velx(i,j,k) * 0.5 * (dens_arr(i,j,k,Rho_comp) + dens_arr(i-1,j,k,Rho_comp));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momy(i,j,k) = vely(i,j,k) * 0.5 * (dens_arr(i,j,k,Rho_comp) + dens_arr(i,j-1,k,Rho_comp));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momz(i,j,k) = velz(i,j,k) * 0.5 * (dens_arr(i,j,k,Rho_comp) + dens_arr(i,j,k-1,Rho_comp));
        });
    } // end MFIter
}

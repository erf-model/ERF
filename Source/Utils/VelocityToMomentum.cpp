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

void VelocityToMomentum (const MultiFab& xvel_in,
                         const  IntVect& xvel_ngrow,
                         const MultiFab& yvel_in,
                         const  IntVect& yvel_ngrow,
                         const MultiFab& zvel_in,
                         const  IntVect& zvel_ngrow,
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
        Box tbx, tby, tbz;
        if (l_use_ndiff) {
            tbx = mfi.tilebox(IntVect(1,0,0),xvel_ngrow);
            tby = mfi.tilebox(IntVect(0,1,0),yvel_ngrow);
            tbz = mfi.tilebox(IntVect(0,0,1),zvel_ngrow);
        } else {
            tbx = mfi.tilebox(IntVect(1,0,0),IntVect(1,1,0));
            tby = mfi.tilebox(IntVect(0,1,0),IntVect(1,1,0));
            tbz = mfi.tilebox(IntVect(0,0,1),IntVect(1,1,0));
        }

        // Conserved/state variables on cell centers -- we use this for density
        const Array4<const Real>& dens_arr = density.array(mfi);

        // Momentum on faces, to be computed
        Array4<Real> const& momx = xmom.array(mfi);
        Array4<Real> const& momy = ymom.array(mfi);
        Array4<Real> const& momz = zmom.array(mfi);

        // Velocity on faces, used in computation
        const Array4<Real const>& velx = xvel_in.const_array(mfi);
        const Array4<Real const>& vely = yvel_in.const_array(mfi);
        const Array4<Real const>& velz = zvel_in.const_array(mfi);

        ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real rho_x_face = 0.5 * (dens_arr(i,j,k,Rho_comp) + dens_arr(i-1,j,k,Rho_comp));
            momx(i,j,k) = velx(i,j,k) * rho_x_face;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real rho_y_face = 0.5 * (dens_arr(i,j,k,Rho_comp) + dens_arr(i,j-1,k,Rho_comp));
            momy(i,j,k) = vely(i,j,k) * rho_y_face;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real rho_z_face = 0.5 * (dens_arr(i,j,k,Rho_comp) + dens_arr(i,j,k-1,Rho_comp));
            momz(i,j,k) = velz(i,j,k) * rho_z_face;
        });
    } // end MFIter
}

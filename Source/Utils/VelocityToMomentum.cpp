/**
 * \file VelocityToMomentum.cpp
 */
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <Utils.H>
#include <AMReX_MultiFabUtil.H>

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
 * @param[in] domain  Domain at this level
 * @param[in] domain_bcs_type_h   host vector for domain boundary conditions
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
                         const Box& domain,
                         const Vector<BCRec>& domain_bcs_type_h)
{
    BL_PROFILE_VAR("VelocityToMomentum()",VelocityToMomentum);

    const BCRec* bc_ptr_h = domain_bcs_type_h.data();

    amrex::Print() <<" COMING IN WITH NG " << xvel_ngrow << " " << yvel_ngrow << " " << zvel_ngrow << std::endl;

    // Loop over boxes = valid boxes grown by ngrow
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // We need momentum in the interior ghost cells (init == real)
        const Box& bx = mfi.tilebox();
        Box tbx, tby, tbz;

        tbx = mfi.tilebox(IntVect(1,0,0),xvel_ngrow);
        tby = mfi.tilebox(IntVect(0,1,0),yvel_ngrow);
        tbz = mfi.tilebox(IntVect(0,0,1),zvel_ngrow);

#if 0
        if (l_use_ndiff) {
            tbx = mfi.tilebox(IntVect(1,0,0),xvel_ngrow);
            tby = mfi.tilebox(IntVect(0,1,0),yvel_ngrow);
            tbz = mfi.tilebox(IntVect(0,0,1),zvel_ngrow);
        } else {
            tbx = mfi.tilebox(IntVect(1,0,0),IntVect(1,1,1));
            if (tbx.smallEnd(2) < 0) tbx.setSmall(2,0);
            tby = mfi.tilebox(IntVect(0,1,0),IntVect(1,1,1));
            if (tby.smallEnd(2) < 0) tby.setSmall(2,0);
            tbz = mfi.tilebox(IntVect(0,0,1),IntVect(1,1,0));
        }
#endif

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
            momx(i,j,k) = velx(i,j,k) * 0.5 * (dens_arr(i,j,k,Rho_comp) + dens_arr(i-1,j,k,Rho_comp));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momy(i,j,k) = vely(i,j,k) * 0.5 * (dens_arr(i,j,k,Rho_comp) + dens_arr(i,j-1,k,Rho_comp));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momz(i,j,k) = velz(i,j,k) * 0.5 * (dens_arr(i,j,k,Rho_comp) + dens_arr(i,j,k-1,Rho_comp));
        });

        if ( (bx.smallEnd(0) == domain.smallEnd(0)) &&
             (bc_ptr_h[BCVars::cons_bc].lo(0) == ERFBCType::ext_dir) ) {
            ParallelFor(makeSlab(tbx,0,domain.smallEnd(0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                momx(i,j,k) = velx(i,j,k) * dens_arr(i-1,j,k,Rho_comp);
            });
        }
        if ( (bx.bigEnd(0) == domain.bigEnd(0)) &&
             (bc_ptr_h[BCVars::cons_bc].hi(0) == ERFBCType::ext_dir) ) {
            ParallelFor(makeSlab(tbx,0,domain.bigEnd(0)+1), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                momx(i,j,k) = velx(i,j,k) * dens_arr(i,j,k,Rho_comp);
            });
        }
        if ( (bx.smallEnd(1) == domain.smallEnd(1)) &&
             (bc_ptr_h[BCVars::cons_bc].lo(1) == ERFBCType::ext_dir) ) {
            ParallelFor(makeSlab(tby,1,domain.smallEnd(1)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                momy(i,j,k) = vely(i,j,k) * dens_arr(i,j-1,k,Rho_comp);
            });
        }
        if ( (bx.bigEnd(1) == domain.bigEnd(1)) &&
             (bc_ptr_h[BCVars::cons_bc].hi(1) == ERFBCType::ext_dir) ) {
            ParallelFor(makeSlab(tby,1,domain.bigEnd(1)+1), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                momy(i,j,k) = vely(i,j,k) * dens_arr(i,j,k,Rho_comp);
            });
        }
    } // end MFIter
}

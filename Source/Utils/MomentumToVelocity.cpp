/**
 * \file MomentumToVelocity.cpp
 */
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <Utils.H>

using namespace amrex;

/**
 * Convert momentum to velocity by dividing by density averaged onto faces
 *
 * @param[out] xvel x-component of velocity
 * @param[out] yvel y-component of velocity
 * @param[out] zvel z-component of velocity
 * @param[in] density density at cell centers
 * @param[in] xmom_in x-component of momentum
 * @param[in] ymom_in y-component of momentum
 * @param[in] zmom_in z-component of momentum
 * @param[in] domain  Domain at this level
 * @param[in] domain_bcs_type_h   host vector for domain boundary conditions
 */

void
MomentumToVelocity (MultiFab& xvel, MultiFab& yvel, MultiFab& zvel,
                    const MultiFab& density,
                    const MultiFab& xmom_in, const MultiFab& ymom_in, const MultiFab& zmom_in,
                    const Box& domain,
                    const Vector<BCRec>& domain_bcs_type_h)
{
    BL_PROFILE_VAR("MomentumToVelocity()",MomentumToVelocity);

    const BCRec* bc_ptr_h = domain_bcs_type_h.data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // We need velocity in the interior ghost cells (init == real)
        Box bx = mfi.tilebox();

        const Box& tbx = surroundingNodes(bx,0);
        const Box& tby = surroundingNodes(bx,1);
        const Box& tbz = surroundingNodes(bx,2);

        // Conserved variables on cell centers -- we use this for density
        const Array4<const Real>& dens_arr = density.array(mfi);

        // Momentum on faces
        Array4<Real const> const& momx = xmom_in.const_array(mfi);
        Array4<Real const> const& momy = ymom_in.const_array(mfi);
        Array4<Real const> const& momz = zmom_in.const_array(mfi);

        // Velocity on faces
        const Array4<Real>& velx = xvel.array(mfi);
        const Array4<Real>& vely = yvel.array(mfi);
        const Array4<Real>& velz = zvel.array(mfi);

        ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            velx(i,j,k) = momx(i,j,k) * 2.0 / (dens_arr(i,j,k,Rho_comp) + dens_arr(i-1,j,k,Rho_comp));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            vely(i,j,k) = momy(i,j,k) * 2.0 / (dens_arr(i,j,k,Rho_comp) + dens_arr(i,j-1,k,Rho_comp));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            velz(i,j,k) = momz(i,j,k) * 2.0 / (dens_arr(i,j,k,Rho_comp) + dens_arr(i,j,k-1,Rho_comp));
        });

        if ( (bx.smallEnd(0) == domain.smallEnd(0)) &&
             (bc_ptr_h[BCVars::cons_bc].lo(0) == ERFBCType::ext_dir) ) {
            ParallelFor(makeSlab(tbx,0,domain.smallEnd(0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                velx(i,j,k) = momx(i,j,k) / dens_arr(i-1,j,k,Rho_comp);
            });
        }
        if ( (bx.bigEnd(0) == domain.bigEnd(0)) &&
             (bc_ptr_h[BCVars::cons_bc].hi(0) == ERFBCType::ext_dir) ) {
            ParallelFor(makeSlab(tbx,0,domain.bigEnd(0)+1), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                velx(i,j,k) = momx(i,j,k) / dens_arr(i,j,k,Rho_comp);
            });
        }
        if ( (bx.smallEnd(1) == domain.smallEnd(1)) &&
             (bc_ptr_h[BCVars::cons_bc].lo(1) == ERFBCType::ext_dir) ) {
            ParallelFor(makeSlab(tby,1,domain.smallEnd(1)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                vely(i,j,k) = momy(i,j,k) / dens_arr(i,j-1,k,Rho_comp);
            });
        }
        if ( (bx.bigEnd(0) == domain.bigEnd(0)) &&
             (bc_ptr_h[BCVars::cons_bc].hi(0) == ERFBCType::ext_dir) ) {
            ParallelFor(makeSlab(tby,1,domain.bigEnd(1)+1), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                vely(i,j,k) = momy(i,j,k) / dens_arr(i,j,k,Rho_comp);
            });
        }
    } // end MFIter
}

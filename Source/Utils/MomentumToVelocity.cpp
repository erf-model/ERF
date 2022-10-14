/**
 * \file MomentumToVelocity.cpp
 */
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <Utils.H>

using namespace amrex;

/**
 * Convert updated momentum to updated velocity
 */
void
MomentumToVelocity( BoxArray& grids_to_evolve,
                    MultiFab& xvel, MultiFab& yvel, MultiFab& zvel,
                    const MultiFab& cons_in,
                    const MultiFab& xmom_in, const MultiFab& ymom_in, const MultiFab& zmom_in)
{
    BL_PROFILE_VAR("MomentumToVelocity()",MomentumToVelocity);

    // Loop over boxes = valid boxes grown by ngrow
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& valid_bx = grids_to_evolve[mfi.index()];

        // Construct intersection of current tilebox and valid region for updating
        Box bx = mfi.tilebox() & valid_bx;

        const Box& tbx = surroundingNodes(bx,0);
        const Box& tby = surroundingNodes(bx,1);
        const Box& tbz = surroundingNodes(bx,2);

        // Conserved variables on cell centers -- we use this for density
        const Array4<const Real>& cons = cons_in.array(mfi);

        // Momentum on faces
        Array4<Real const> const& momx = xmom_in.array(mfi);
        Array4<Real const> const& momy = ymom_in.array(mfi);
        Array4<Real const> const& momz = zmom_in.array(mfi);

        // Velocity on faces
        const Array4<Real>& velx = xvel.array(mfi);
        const Array4<Real>& vely = yvel.array(mfi);
        const Array4<Real>& velz = zvel.array(mfi);

        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            velx(i,j,k) = momx(i,j,k)/(0.5 * (cons(i,j,k,Rho_comp) + cons(i-1,j,k,Rho_comp)));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            vely(i,j,k) = momy(i,j,k)/(0.5 * (cons(i,j,k,Rho_comp) + cons(i,j-1,k,Rho_comp)));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            velz(i,j,k) = momz(i,j,k)/(0.5 * (cons(i,j,k,Rho_comp) + cons(i,j,k-1,Rho_comp)));
        });
    } // end MFIter
}

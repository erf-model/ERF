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
MomentumToVelocity( MultiFab& xvel,
                    const IntVect& xvel_ngrow,
                    MultiFab& yvel,
                    const IntVect& yvel_ngrow,
                    MultiFab& zvel,
                    const IntVect& zvel_ngrow,
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
        const Box& tbx = amrex::grow(mfi.nodaltilebox(0),xvel_ngrow);
        const Box& tby = amrex::grow(mfi.nodaltilebox(1),yvel_ngrow);
        const Box& tbz = amrex::grow(mfi.nodaltilebox(2),zvel_ngrow);

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

/**
 * \file VelocityToMomentum.cpp
 */
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <Utils.H>

using namespace amrex;

/**
 * Convert velocity to momentum.
 */
void VelocityToMomentum( const MultiFab& xvel_in,
                         const IntVect& xvel_ngrow,
                         const MultiFab& yvel_in,
                         const IntVect& yvel_ngrow,
                         const MultiFab& zvel_in,
                         const IntVect& zvel_ngrow,
                         const MultiFab& cons_in,
                         MultiFab& xmom, MultiFab& ymom, MultiFab& zmom)
{
    BL_PROFILE_VAR("VelocityToMomentum()",VelocityToMomentum);

    // Loop over boxes = valid boxes grown by ngrow
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box tbx = amrex::grow(mfi.nodaltilebox(0),xvel_ngrow); tbx.setSmall(2,0);
        Box tby = amrex::grow(mfi.nodaltilebox(1),yvel_ngrow); tby.setSmall(2,0);
        Box tbz = amrex::grow(mfi.nodaltilebox(2),zvel_ngrow); tbz.setSmall(2,0);

        // Conserved/state variables on cell centers -- we use this for density
        const Array4<const Real>& cons = cons_in.array(mfi);

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
            momx(i,j,k) = velx(i,j,k) * 0.5 * (cons(i,j,k,Rho_comp) + cons(i-1,j,k,Rho_comp));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momy(i,j,k) = vely(i,j,k) * 0.5 * (cons(i,j,k,Rho_comp) + cons(i,j-1,k,Rho_comp));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momz(i,j,k) = velz(i,j,k) * 0.5 * (cons(i,j,k,Rho_comp) + cons(i,j,k-1,Rho_comp));
        });
    } // end MFIter
}

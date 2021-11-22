/**
 * \file MomentumToVelocity.cpp
 */
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <TimeIntegration.H>

using namespace amrex;


/**
 * Convert updated momentum to updated velocity
 */
void MomentumToVelocity( MultiFab& xvel, MultiFab& yvel, MultiFab& zvel,
                         MultiFab& cons_in,
                         MultiFab& xmom_in, MultiFab& ymom_in, MultiFab& zmom_in,
                         int ng,
                         const SolverChoice& solverChoice)
{
    BL_PROFILE_VAR("MomentumToVelocity()",MomentumToVelocity);

    // Loop over boxes
    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& tbx = amrex::grow(mfi.nodaltilebox(0),IntVect(0,ng,ng));
        const Box& tby = amrex::grow(mfi.nodaltilebox(1),IntVect(ng,0,ng));
        const Box& tbz = amrex::grow(mfi.nodaltilebox(2),IntVect(ng,ng,0));

        // Conserved variables on cell centers -- we use this for density
        const Array4<Real>& cons = cons_in.array(mfi);

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
            velx(i,j,k) = momx(i,j,k)/InterpolateDensityFromCellToFace(
                                        i, j, k, cons, NextOrPrev::prev, Coord::x,
                                        solverChoice.spatial_order);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            vely(i,j,k) = momy(i,j,k)/InterpolateDensityFromCellToFace(
                                        i, j, k, cons, NextOrPrev::prev, Coord::y,
                                        solverChoice.spatial_order);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            velz(i,j,k) = momz(i,j,k)/InterpolateDensityFromCellToFace(
                                        i, j, k, cons, NextOrPrev::prev, Coord::z,
                                        solverChoice.spatial_order);
        });
    } // end MFIter
}

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
                         const MultiFab& cons_in,
                         const MultiFab& xmom_in, const MultiFab& ymom_in, const MultiFab& zmom_in,
                         const int l_spatial_order, const int ngrow)
{
    BL_PROFILE_VAR("MomentumToVelocity()",MomentumToVelocity);

    // Loop over boxes = valid boxes grown by ngrow
    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& tbx = amrex::grow(mfi.nodaltilebox(0),IntVect(0,ngrow,ngrow));
        const Box& tby = amrex::grow(mfi.nodaltilebox(1),IntVect(ngrow,0,ngrow));
        const Box& tbz = amrex::grow(mfi.nodaltilebox(2),IntVect(ngrow,ngrow,0));

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
            velx(i,j,k) = momx(i,j,k)/InterpolateFromCellOrFace(
                                        i, j, k, cons, Rho_comp, momx(i,j,k), Coord::x, l_spatial_order);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            vely(i,j,k) = momy(i,j,k)/InterpolateFromCellOrFace(
                                        i, j, k, cons, Rho_comp, momy(i,j,k), Coord::y, l_spatial_order);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            velz(i,j,k) = momz(i,j,k)/InterpolateFromCellOrFace(
                                        i, j, k, cons, Rho_comp, momz(i,j,k), Coord::z, l_spatial_order);
        });
    } // end MFIter
}

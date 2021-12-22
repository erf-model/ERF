/**
 * \file VelocityToMomentum.cpp
 */
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <SpatialStencils.H>

using namespace amrex;

/**
 * Convert velocity to momentum.
 */
void VelocityToMomentum( const MultiFab& xvel_in, const MultiFab& yvel_in, const MultiFab& zvel_in,
                         const MultiFab& cons_in,
                         MultiFab& xmom, MultiFab& ymom, MultiFab& zmom,
                         const int l_spatial_order, const int ngrow)
{
    BL_PROFILE_VAR("VelocityToMomentum()",VelocityToMomentum);

    // Loop over boxes = valid boxes grown by ngrow
    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& tbx = amrex::grow(mfi.nodaltilebox(0),IntVect(0,ngrow,ngrow));
        const Box& tby = amrex::grow(mfi.nodaltilebox(1),IntVect(ngrow,0,ngrow));
        const Box& tbz = amrex::grow(mfi.nodaltilebox(2),IntVect(ngrow,ngrow,0));

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
            momx(i,j,k) = velx(i,j,k)* InterpolateFromCellOrFace(
                                i, j, k, cons, Rho_comp, velx(i,j,k), Coord::x, l_spatial_order);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momy(i,j,k) = vely(i,j,k)* InterpolateFromCellOrFace(
                                i, j, k, cons, Rho_comp, vely(i,j,k), Coord::y, l_spatial_order);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momz(i,j,k) = velz(i,j,k)* InterpolateFromCellOrFace(
                                i, j, k, cons, Rho_comp, velz(i,j,k), Coord::z, l_spatial_order);
        });
    } // end MFIter
}

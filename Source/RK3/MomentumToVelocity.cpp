#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

#include "IndexDefines.H"

using namespace amrex;

void MomentumToVelocity( MultiFab& xvel, MultiFab& yvel, MultiFab& zvel, 
                         MultiFab& cons_in, 
                         MultiFab& xmom_in, MultiFab& ymom_in, MultiFab& zmom_in)
{
    BL_PROFILE_VAR("MomentumToVelocity()",MomentumToVelocity);

    // Loop over boxes
    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) 
    {
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

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
            velx(i,j,k) = 2.*momx(i,j,k)/(cons(i,j,k,Density_comp) + cons(i-1,j,k,Density_comp));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            vely(i,j,k) = 2.*momy(i,j,k)/(cons(i,j,k,Density_comp) + cons(i,j-1,k,Density_comp));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            velz(i,j,k) = 2.*momz(i,j,k)/(cons(i,j,k,Density_comp) + cons(i,j,k-1,Density_comp));
        });
    } // end MFIter
}

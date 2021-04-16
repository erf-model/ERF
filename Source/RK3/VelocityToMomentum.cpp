#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

#include "IndexDefines.H"

using namespace amrex;

void VelocityToMomentum( MultiFab& xvel_in, MultiFab& yvel_in, MultiFab& zvel_in, 
                         MultiFab& cons_in, 
                         MultiFab& xmom, MultiFab& ymom, MultiFab& zmom)
{
    BL_PROFILE_VAR("VelocityToMomentum()",VelocityToMomentum);

    // Loop over boxes
    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) 
    {
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        // Conserved variables on cell centers -- we use this for density
        const Array4<Real>& cons = cons_in.array(mfi);

        // Momentum on faces
        Array4<Real> const& momx = xmom.array(mfi);
        Array4<Real> const& momy = ymom.array(mfi);
        Array4<Real> const& momz = zmom.array(mfi);

        // Velocity on faces
        const Array4<Real const>& velx = xvel_in.array(mfi);
        const Array4<Real const>& vely = yvel_in.array(mfi);
        const Array4<Real const>& velz = zvel_in.array(mfi);

        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momx(i,j,k) = 0.5*velx(i,j,k)*(cons(i,j,k,Density_comp) + cons(i-1,j,k,Density_comp));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momy(i,j,k) = 0.5*vely(i,j,k)*(cons(i,j,k,Density_comp) + cons(i,j-1,k,Density_comp));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momz(i,j,k) = 0.5*velz(i,j,k)*(cons(i,j,k,Density_comp) + cons(i,j,k-1,Density_comp));
        });
    } // end MFIter
}

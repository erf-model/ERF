#include <AMReX_MultiFabUtil.H>

using namespace amrex;

void ConservedToPrimitive(MultiFab& prim_out, 
                          MultiFab& xvel_out, MultiFab& yvel_out, MultiFab& zvel_out,
                          MultiFab& cons_in,  
                          MultiFab& xmom_in , MultiFab& ymom_in, MultiFab& zmom_in)
{
    BL_PROFILE_VAR("ConservedToPrimitive()",ConservedToPrimitive);

    // For now set every component to zero, just in case there are components we're not setting and not using
    prim_out.setVal(0.0);

    // Loop over boxes
    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        // Conserved variables on cell centers
        const Array4<      Real>& cons = cons_in.array(mfi);

        // Primitive variables on cell centers
        const Array4<      Real>& prim = prim_out.array(mfi);

        // Momentum on faces
        Array4<Real const> const& momx = xmom_in.array(mfi);
        Array4<Real const> const& momy = ymom_in.array(mfi);
        Array4<Real const> const& momz = zmom_in.array(mfi);

        // Velocity on faces
        const Array4<Real>& velx = xvel_out.array(mfi);
        const Array4<Real>& vely = yvel_out.array(mfi);
        const Array4<Real>& velz = zvel_out.array(mfi);

        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            velx(i,j,k) = 2.*momx(i,j,k)/(cons(i,j,k,0) + cons(i-1,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            vely(i,j,k) = 2.*momy(i,j,k)/(cons(i,j,k,0) + cons(i,j-1,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            velz(i,j,k) = 2.*momz(i,j,k)/(cons(i,j,k,0) + cons(i,j,k-1,0));
        });

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Density
            prim(i,j,k,0) = cons(i,j,k,0);

            // Velocity
            prim(i,j,k,1) = 0.5*(velx(i,j,k) + velx(i+1,j,k));
            prim(i,j,k,2) = 0.5*(vely(i,j,k) + vely(i,j+1,k));
            prim(i,j,k,3) = 0.5*(velz(i,j,k) + velz(i,j,k+1));

            // Cell-centered momentum
            cons(i,j,k,1) = 0.5*(momx(i,j,k) + momx(i+1,j,k));
            cons(i,j,k,2) = 0.5*(momy(i,j,k) + momy(i,j+1,k));
            cons(i,j,k,3) = 0.5*(momz(i,j,k) + momz(i,j,k+1));

            Real kinenergy = 0.;
            kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
            kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
            kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
            kinenergy *= (0.125/cons(i,j,k,0));

            // internal energy e = E - 1/2 magvel^2
            Real intenergy = (cons(i,j,k,4)-kinenergy)/cons(i,j,k,0);

            // temperature = e / cv
            Real cv_air = 0.718; // (kJ/kg/K)
            prim(i,j,k,4) = intenergy/cv_air;

            // pressure = rho * R * T
            Real Runiv = 8.314; // (J / K / mol)
            prim(i,j,k,5) = cons(i,j,k,0)*Runiv*prim(i,j,k,4);

            // Advected scalar S = (rho S) / rho 
            prim(i,j,k,6) = cons(i,j,k,7) / cons(i,j,k,0);
        });
        
    } // end MFIter
}

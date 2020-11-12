#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

using namespace amrex;

void conservedToPrimitiveStag(MultiFab& prim_in, 
                              MultiFab& u_x, MultiFab& v_y, MultiFab& w_z, 
                              MultiFab& cons_in, 
                              MultiFab& cu_x, MultiFab& cu_y, MultiFab& cu_z)
{
    BL_PROFILE_VAR("conservedToPrimitiveStag()",conservedToPrimitiveStag);

    // Loop over boxes
    for ( MFIter mfi(prim_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<      Real>& cons = cons_in.array(mfi);
        const Array4<      Real>& prim = prim_in.array(mfi);

        // Conserved variables on faces
        Array4<Real const> const& momx = cu_x.array(mfi);
        Array4<Real const> const& momy = cu_y.array(mfi);
        Array4<Real const> const& momz = cu_z.array(mfi);

        // Momentum on faces
        const Array4<Real>& velx = u_x.array(mfi);
        const Array4<Real>& vely = v_y.array(mfi);
        const Array4<Real>& velz = w_z.array(mfi);

        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            velx(i,j,k) = 2*momx(i,j,k)/(cons(i,j,k,0) + cons(i-1,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            vely(i,j,k) = 2*momy(i,j,k)/(cons(i,j,k,0) + cons(i,j-1,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            velz(i,j,k) = 2*momz(i,j,k)/(cons(i,j,k,0) + cons(i,j,k-1,0));
        });
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            prim(i,j,k,0) = cons(i,j,k,0);

            prim(i,j,k,1) = 0.5*(velx(i,j,k) + velx(i+1,j,k));
            prim(i,j,k,2) = 0.5*(vely(i,j,k) + vely(i,j+1,k));
            prim(i,j,k,3) = 0.5*(velz(i,j,k) + velz(i,j,k+1));

            cons(i,j,k,1) = 0.5*(momx(i,j,k) + momx(i+1,j,k));
            cons(i,j,k,2) = 0.5*(momy(i,j,k) + momy(i,j+1,k));
            cons(i,j,k,3) = 0.5*(momz(i,j,k) + momz(i,j,k+1));

            Real kinenergy = 0.;
            kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
            kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
            kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
            kinenergy *= (0.125/cons(i,j,k,0));

            // e = E - 1/2 magvel^2
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

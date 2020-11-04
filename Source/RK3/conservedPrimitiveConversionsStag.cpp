#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

using namespace amrex;

void conservedToPrimitiveStag(MultiFab& prim_in, std::array<MultiFab, AMREX_SPACEDIM>& velStag_in,
                              MultiFab& cons_in, const std::array<MultiFab, AMREX_SPACEDIM>& momStag_in)
{
    BL_PROFILE_VAR("conservedToPrimitiveStag()",conservedToPrimitiveStag);

    // from namelist
    /* 
    // method 1 to create a thread shared array
    // must use if the size of the array is not known at compile time
    // note when passing this into a function, you need to use the type,
    // "Real const * const AMREX_RESTRICT"
    */
    
    // Loop over boxes
    for ( MFIter mfi(prim_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        //const Array4<const Real>& cons = cons_in.array(mfi);
        const Array4<      Real>& cons = cons_in.array(mfi);
        const Array4<      Real>& prim = prim_in.array(mfi);

        AMREX_D_TERM(Array4<Real const> const& momx = momStag_in[0].array(mfi);,
                     Array4<Real const> const& momy = momStag_in[1].array(mfi);,
                     Array4<Real const> const& momz = momStag_in[2].array(mfi););

        AMREX_D_TERM(const Array4<Real>& velx = velStag_in[0].array(mfi);,
                     const Array4<Real>& vely = velStag_in[1].array(mfi);,
                     const Array4<Real>& velz = velStag_in[2].array(mfi););

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
            // method 2 to create a thread private array
            // can use when the size of the array is known at compile-time
            
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

            // Do we need to calculate staggered velocities here as well? (from rho averaged to all faces) -- Ishan

            Real intenergy = (cons(i,j,k,4)-kinenergy)/cons(i,j,k,0);

            // update temperature and pressure in-place using internal energy

            Real cv_air = 0.718; // (kJ/kg/K)
            prim(i,j,k,4) = intenergy/cv_air;

            Real Runiv = 8.314; // (J / K / mol)
            prim(i,j,k,5) = cons(i,j,k,0)*Runiv*prim(i,j,k,4);
        });
        
    } // end MFIter
}

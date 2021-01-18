#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_VisMF.H>

#include <ERF.H>
#include <RK3.H>

using namespace amrex;

void RK3_advance(MultiFab& cons_old,  MultiFab& cons_new,
                 MultiFab& xvel_old, MultiFab& yvel_old, MultiFab& zvel_old, 
                 MultiFab& xvel_new, MultiFab& yvel_new, MultiFab& zvel_new, 
                 MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                 std::array< MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                 const amrex::Geometry geom, const amrex::Real* dxp, const amrex::Real dt)
{
    BL_PROFILE_VAR("RK3stepStag()",RK3stepStag);

    int nvars = cons_old.nComp(); 
    int ngc = 1;

    // ************************************************************************************** 
    // 
    // Allocate temporary MultiFab to hold the primitive variables
    // 
    // ************************************************************************************** 
    MultiFab prim(cons_old.boxArray(),cons_old.DistributionMap(),nvars,2); 

    // ************************************************************************************** 
    // 
    // Allocate temporary MultiFabs to hold the intermediate updates
    // 
    // ************************************************************************************** 

    const BoxArray& ba            = cons_old.boxArray();
    const DistributionMapping& dm = cons_old.DistributionMap();

    // Intermediate solutions -- At cell centers
    MultiFab cons_upd_1(ba,dm,nvars,ngc);
    MultiFab cons_upd_2(ba,dm,nvars,ngc);
    cons_upd_1.setVal(0.0,0,nvars,ngc);
    cons_upd_2.setVal(0.0,0,nvars,ngc);

    // Intermediate solutions -- On faces
    MultiFab xmom_update_1(convert(ba,IntVect(1,0,0)), dm, 1, 1);
    MultiFab ymom_update_1(convert(ba,IntVect(0,1,0)), dm, 1, 1);
    MultiFab zmom_update_1(convert(ba,IntVect(0,0,1)), dm, 1, 1);

    MultiFab xmom_update_2(convert(ba,IntVect(1,0,0)), dm, 1, 1);
    MultiFab ymom_update_2(convert(ba,IntVect(0,1,0)), dm, 1, 1);
    MultiFab zmom_update_2(convert(ba,IntVect(0,0,1)), dm, 1, 1);

    // Initial momentum -- on faces
    MultiFab xmom_old(convert(ba,IntVect(1,0,0)), dm, 1, 1);
    MultiFab ymom_old(convert(ba,IntVect(0,1,0)), dm, 1, 1);
    MultiFab zmom_old(convert(ba,IntVect(0,0,1)), dm, 1, 1);

    // Final momentum -- on faces
    MultiFab xmom_new(convert(ba,IntVect(1,0,0)), dm, 1, 1);
    MultiFab ymom_new(convert(ba,IntVect(0,1,0)), dm, 1, 1);
    MultiFab zmom_new(convert(ba,IntVect(0,0,1)), dm, 1, 1);

    // ************************************************************************************** 
    // 
    // Convert old velocity to old momentum on faces
    // 
    // ************************************************************************************** 
    cons_old.FillBoundary(geom.periodicity());

    VelocityToMomentum(xvel_old, yvel_old, zvel_old, cons_old, xmom_old, ymom_old, zmom_old);

    xmom_old.FillBoundary(geom.periodicity());
    ymom_old.FillBoundary(geom.periodicity());
    zmom_old.FillBoundary(geom.periodicity());

    // ************************************************************************************** 

    Real rho0 = 1.0;
    cons_upd_1.setVal(rho0,0,1,ngc);
    cons_upd_2.setVal(rho0,0,1,ngc);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, 0.0};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // ************************************************************************************** 
    //
    // Return update in the cons_upd_1 and *mom_update_1 MultiFabs
    //
    // ************************************************************************************** 
    RK3_stage(cons_old, cons_upd_1,
              xmom_old, ymom_old, zmom_old, 
              xmom_update_1, ymom_update_1, zmom_update_1, 
              xvel_new, yvel_new, zvel_new, prim,  // These are used as temporary space
              source,
              eta, zeta, kappa,
              faceflux,
              edgeflux_x,
              edgeflux_y,
              edgeflux_z,
              cenflux,
              geom, dxp, dt);

     
    // ************************************************************************************** 
    // 
    // Define updates in the first RK stage
    // 
    // ************************************************************************************** 
    for ( MFIter mfi(cons_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cu_old     = cons_old.array(mfi);
        const Array4<Real> & cu_up1     = cons_upd_1.array(mfi);

        const Array4<Real>& momx_old = xmom_old.array(mfi);
        const Array4<Real>& momy_old = ymom_old.array(mfi);
        const Array4<Real>& momz_old = zmom_old.array(mfi);

        const Array4<Real>& momx_up1 = xmom_update_1.array(mfi);
        const Array4<Real>& momy_up1 = ymom_update_1.array(mfi);
        const Array4<Real>& momz_up1 = zmom_update_1.array(mfi);

        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            // At this point, the update is stored in cu_up1 ... 
            cu_up1(i,j,k,n) += cu_old(i,j,k,n);

            // Now cu_up1 holds the first intermediate solution

        }); // [1:3 indices are not valuable -- momentum flux]

        // momentum flux
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // At this point, the update is stored in momx_up1 ... 
            momx_up1(i,j,k) += momx_old(i,j,k);
            // Now momx_up1 holds the first intermediate solution
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // At this point, the update is stored in momy_up1 ... 
            momy_up1(i,j,k) += momy_old(i,j,k);
            // Now momy_up1 holds the first intermediate solution
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // At this point, the update is stored in momz_up1 ... 
            momz_up1(i,j,k) += momz_old(i,j,k);
            // Now momz_up1 holds the first intermediate solution
        });
    }

    xmom_update_1.FillBoundary(geom.periodicity());
    ymom_update_1.FillBoundary(geom.periodicity());
    zmom_update_1.FillBoundary(geom.periodicity());

    cons_upd_1.FillBoundary(geom.periodicity());
 
    // ************************************************************************************** 
    //
    // Return update in the cons_upd_2 and *mom_update_2 MultiFabs
    //
    // ************************************************************************************** 
    RK3_stage(cons_upd_1, cons_upd_2,
              xmom_update_1, ymom_update_1, zmom_update_1, 
              xmom_update_2, ymom_update_2, zmom_update_2, 
              xvel_new, yvel_new, zvel_new, prim,  // These are used as temporary space
              source,
              eta, zeta, kappa,
              faceflux,
              edgeflux_x,
              edgeflux_y,
              edgeflux_z,
              cenflux,
              geom, dxp, dt);

    // ************************************************************************************** 
    // 
    // Define updates in the second RK stage
    // 
    // ************************************************************************************** 

    for ( MFIter mfi(cons_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cu_old     = cons_old.array(mfi);
        const Array4<Real> & cu_up1 = cons_upd_1.array(mfi);
        const Array4<Real> & cu_up2 = cons_upd_2.array(mfi);

        const Array4<Real>& momx_old = xmom_old.array(mfi);
        const Array4<Real>& momy_old = ymom_old.array(mfi);
        const Array4<Real>& momz_old = zmom_old.array(mfi);

        const Array4<Real>& momx_up1 = xmom_update_1.array(mfi);
        const Array4<Real>& momy_up1 = ymom_update_1.array(mfi);
        const Array4<Real>& momz_up1 = zmom_update_1.array(mfi);

        const Array4<Real>& momx_up2 = xmom_update_2.array(mfi);
        const Array4<Real>& momy_up2 = ymom_update_2.array(mfi);
        const Array4<Real>& momz_up2 = zmom_update_2.array(mfi);

        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            // At this point, the update is stored in cu_up2 ... 
            cu_up2(i,j,k,n) = 0.75*cu_old(i,j,k,n) + 0.25*(cu_up1(i,j,k,n) + cu_up2(i,j,k,n));

            // Now cu_up2 holds the second intermediate solution

        }); // [1:3 indices are not valuable -- momentum flux]

        // momentum flux
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // At this point, the update is stored in momx_up2 ... 
            momx_up2(i,j,k) = 0.75*momx_old(i,j,k) + 0.25*(momx_up1(i,j,k) + momx_up2(i,j,k));
            // Now momx_up2 holds the second intermediate solution
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // At this point, the update is stored in momy_up2 ... 
            momy_up2(i,j,k) = 0.75*momy_old(i,j,k) + 0.25*(momy_up1(i,j,k) + momy_up2(i,j,k));
            // Now momy_up2 holds the second intermediate solution
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // At this point, the update is stored in momz_up2 ... 
            momz_up2(i,j,k) = 0.75*momz_old(i,j,k) + 0.25*(momz_up1(i,j,k) + momz_up2(i,j,k));
            // Now momz_up2 holds the second intermediate solution
        });
    }
        
    xmom_update_2.FillBoundary(geom.periodicity());
    ymom_update_2.FillBoundary(geom.periodicity());
    zmom_update_2.FillBoundary(geom.periodicity());

    cons_upd_2.FillBoundary(geom.periodicity());

    // ************************************************************************************** 
    //
    // Return update in the cons_new and *mom_new MultiFabs
    //
    // ************************************************************************************** 
    RK3_stage(cons_upd_2, cons_new, 
              xmom_update_2, ymom_update_2, zmom_update_2, 
              xmom_new, ymom_new, zmom_new, 
              xvel_new, yvel_new, zvel_new, prim,  // These are used as temporary space
              source,
              eta, zeta, kappa,
              faceflux,   // These are just temporary space
              edgeflux_x, // These are just temporary space
              edgeflux_y, // These are just temporary space
              edgeflux_z, // These are just temporary space
              cenflux,    // These are just temporary space
              geom, dxp, dt);

    for ( MFIter mfi(cons_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cu_old = cons_old.array(mfi);
        const Array4<Real> & cu_up2 = cons_upd_2.array(mfi);
        const Array4<Real> & cu_new = cons_new.array(mfi);

        const Array4<Real>& momx_old = xmom_old.array(mfi);
        const Array4<Real>& momy_old = ymom_old.array(mfi);
        const Array4<Real>& momz_old = zmom_old.array(mfi);

        const Array4<Real>& momx_new = xmom_new.array(mfi);
        const Array4<Real>& momy_new = ymom_new.array(mfi);
        const Array4<Real>& momz_new = zmom_new.array(mfi);

        const Array4<Real>& momx_up2 = xmom_update_2.array(mfi);
        const Array4<Real>& momy_up2 = ymom_update_2.array(mfi);
        const Array4<Real>& momz_up2 = zmom_update_2.array(mfi);

        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            // At this point, the update is stored in cu_new ... 
            cu_new(i,j,k,n) = (1./3.)*cu_old(i,j,k,n) + (2./3.)*(cu_up2(i,j,k,n) + cu_new(i,j,k,n));

            // Now cu_new holds the final solution
            
        }); // [1:3 indices are not valuable -- momentum flux]

        // momentum flux
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // At this point, the update is stored in momx_new ... 
            momx_new(i,j,k) = (1./3.)*momx_old(i,j,k) + (2./3.)*(momx_up2(i,j,k) + momx_new(i,j,k));
            // Now momx_new holds the final solution
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // At this point, the update is stored in momy_new ... 
            momy_new(i,j,k) = (1./3.)*momy_old(i,j,k) + (2./3.)*(momy_up2(i,j,k) + momy_new(i,j,k));
            // Now momy_new holds the final solution
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // At this point, the update is stored in momz_new ... 
            momz_new(i,j,k) = (1./3.)*momz_old(i,j,k) + (2./3.)*(momz_up2(i,j,k) + momz_new(i,j,k));
            // Now momz_new holds the final solution
        });
    }

    // ************************************************************************************** 
    // 
    // Convert new momentum to new velocity on faces
    // 
    // ************************************************************************************** 
    cons_new.FillBoundary(geom.periodicity());
    MomentumToVelocity(xvel_new, yvel_new, zvel_new, cons_new, xmom_new, ymom_new, zmom_new);
}

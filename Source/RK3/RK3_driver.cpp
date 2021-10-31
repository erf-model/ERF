#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>

#include <RK3.H>
#include <utils.H>

using namespace amrex;

void RK3_advance(int level,
                 MultiFab& cons_old,  MultiFab& cons_new,
                 MultiFab& xvel_old, MultiFab& yvel_old, MultiFab& zvel_old,
                 MultiFab& xvel_new, MultiFab& yvel_new, MultiFab& zvel_new,
                 MultiFab& xmom_crse, MultiFab& ymom_crse, MultiFab& zmom_crse,
                 MultiFab& source,
                 std::array< MultiFab, AMREX_SPACEDIM>& faceflux,
                 const amrex::Geometry geom,
                 const amrex::IntVect ref_ratio,
                 const amrex::Real* dxp, const amrex::Real dt,
                 //InterpFaceRegister* ifr,
                 const SolverChoice& solverChoice)
{
    BL_PROFILE_VAR("RK3stepStag()",RK3stepStag);

    int nvars = cons_old.nComp();
    // Determine the number of ghost cells depending on the spatial order
    // **************************************************************************************
    //TODO: Check if this is the only place to specify the number of ghost cells
    //TODO: Also explore how 'ngrow' should be related to the spatial_order
    int ngc = ComputeGhostCells(solverChoice.spatial_order);

    // Allocate temporary MultiFabs to hold the intermediate updates
    // **************************************************************************************
    const BoxArray& ba            = cons_old.boxArray();
    const DistributionMapping& dm = cons_old.DistributionMap();

    // Intermediate solutions (state) -- At cell centers
    MultiFab cons_upd_1(ba,dm,nvars,ngc);
    MultiFab cons_upd_2(ba,dm,nvars,ngc);

    cons_upd_1.setVal(0.0,Rho_comp,nvars,ngc);
    cons_upd_2.setVal(0.0,Rho_comp,nvars,ngc);

    // Intermediate momentum -- On faces
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
    // Prepare for calling the stages of the RK time integration
    // **************************************************************************************

    // Apply BC on old state data at cells before stage 1
    // **************************************************************************************
    cons_old.FillBoundary(geom.periodicity());

    // **************************************************************************************
    // Convert old velocity available on faces to old momentum on faces to be used in time integration
    // **************************************************************************************
    VelocityToMomentum(xvel_old, yvel_old, zvel_old, cons_old, xmom_old, ymom_old, zmom_old, solverChoice);

    // **************************************************************************************
    // HACK HACK HACK -- We copy all of S_old into
    //     cons_upd_1 and cons_upd_2 just so that we can
    //     have values defined in the ghost cells.
    // TODO:  we need to fix this so that cons_upd_1 and cons_upd_2
    //        have interpolated values at the correct time -- not the old time
    // **************************************************************************************
    if (level > 0)
    {
        MultiFab::Copy(cons_upd_1,cons_old,0,0,nvars,ngc);
        MultiFab::Copy(cons_upd_2,cons_old,0,0,nvars,ngc);
        MultiFab::Copy(cons_new  ,cons_old,0,0,nvars,ngc);
    }

    // Apply BC on old momentum data on faces before stage 1
    // **************************************************************************************
    xmom_old.FillBoundary(geom.periodicity());
    ymom_old.FillBoundary(geom.periodicity());
    zmom_old.FillBoundary(geom.periodicity());

#if 0
    if (level > 0)
    {
        amrex::Array<const MultiFab*,3> cmf_const{&xmom_crse, &ymom_crse, &zmom_crse};
        amrex::Array<MultiFab*,3> fmf{&xmom_old , &ymom_old , &zmom_old};

        // Interpolate from coarse faces to fine faces *only* on the coarse-fine boundary
        ifr->interp(fmf,cmf_const,0,1);

        amrex::Array<MultiFab*,3> cmf{&xmom_crse, &ymom_crse, &zmom_crse};

        int nGrow = 1;
        BoxArray fine_grids(cons_old.boxArray());

        // Interpolate from coarse faces on fine faces outside the fine region
        create_umac_grown(level, nGrow, fine_grids, geom, cmf, fmf, ref_ratio);
    }
#endif

    amrex::Vector<MultiFab*> vars_old{&cons_old, &xmom_old, &ymom_old, &zmom_old};

    ERF::applyBCs(geom, vars_old);

    // **************************************************************************************
    // RK3 stage 1: Return update in the cons_upd_1 and [x,y,z]mom_update_1 MultiFabs
    // **************************************************************************************
    RK3_stage(level, cons_old, cons_upd_1,
              xmom_old, ymom_old, zmom_old,
              xmom_update_1, ymom_update_1, zmom_update_1,
              xvel_old, yvel_old, zvel_old,
              source,
              faceflux,
              geom, dxp, dt,
              //ifr,
              solverChoice);

    // **************************************************************************************
    // RK3 stage 1: Define updates in the first RK stage
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
        });

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
    // Apply BC on updated momentum data on faces after stage 1 and before stage 2
    // **************************************************************************************
    xmom_update_1.FillBoundary(geom.periodicity());
    ymom_update_1.FillBoundary(geom.periodicity());
    zmom_update_1.FillBoundary(geom.periodicity());

    // Apply BC on updated state data at cells after stage 1 and before stage 2
    // **************************************************************************************
    cons_upd_1.FillBoundary(geom.periodicity());

    amrex::Vector<MultiFab*> vars1{&cons_upd_1, &xmom_update_1, &ymom_update_1, &zmom_update_1};

    ERF::applyBCs(geom, vars1);

#if 0
    if (level > 0)
    {
        amrex::Array<const MultiFab*,3> cmf_const{&xmom_crse, &ymom_crse, &zmom_crse};
        amrex::Array<MultiFab*,3> fmf{&xmom_update_1 , &ymom_update_1 , &zmom_update_1};

        // Interpolate from coarse faces to fine faces *only* on the coarse-fine boundary
        ifr->interp(fmf,cmf_const,0,1);

        amrex::Array<MultiFab*,3> cmf{&xmom_crse, &ymom_crse, &zmom_crse};

        int nGrow = 1;
        BoxArray fine_grids(cons_old.boxArray());

        // Interpolate from coarse faces on fine faces outside the fine region
        create_umac_grown(level, nGrow, fine_grids, geom, cmf, fmf, ref_ratio);
    }
#endif

    // **************************************************************************************
    // Convert updated momentum to updated velocity on faces after stage 1 and before stage 2
    // Make sure we have filled the ghost cells of cons_upd_1 and *mom_update_1 before this call!
    // **************************************************************************************
    MomentumToVelocity(xvel_new, yvel_new, zvel_new, cons_upd_1, xmom_update_1, ymom_update_1, zmom_update_1, 1, solverChoice);

    // **************************************************************************************
    // RK3 stage 2: Return update in the cons_upd_2 and [x,y,z]mom_update_2 MultiFabs
    // **************************************************************************************
    // TODO: We won't need faceflux, edgeflux, and centflux when using the new code architecture
    RK3_stage(level, cons_upd_1, cons_upd_2,
              xmom_update_1, ymom_update_1, zmom_update_1,
              xmom_update_2, ymom_update_2, zmom_update_2,
              xvel_new, yvel_new, zvel_new,
              source,
              faceflux,
              geom, dxp, dt,
              //ifr,
              solverChoice);

    // **************************************************************************************
    // RK3 stage 2: Define updates in the second RK stage
    // **************************************************************************************
    for ( MFIter mfi(cons_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cu_old     = cons_old.array(mfi);
        const Array4<Real> & cu_up1     = cons_upd_1.array(mfi);
        const Array4<Real> & cu_up2     = cons_upd_2.array(mfi);

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
        });

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
    // TODO: Check if the order of applying BC on cell-centered state or face-centered mom makes any difference
    // Apply BC on updated momentum data on faces after stage 2 and before stage 3
    // **************************************************************************************
    xmom_update_2.FillBoundary(geom.periodicity());
    ymom_update_2.FillBoundary(geom.periodicity());
    zmom_update_2.FillBoundary(geom.periodicity());

    // Apply BC on updated state data at cells after stage 2 and before stage 3
    // **************************************************************************************
    cons_upd_2.FillBoundary(geom.periodicity());

    amrex::Vector<MultiFab*> vars2{&cons_upd_2, &xmom_update_2, &ymom_update_2, &zmom_update_2};

    ERF::applyBCs(geom, vars2);

#if 0
    if (level > 0)
    {
        amrex::Array<const MultiFab*,3> cmf_const{&xmom_crse, &ymom_crse, &zmom_crse};
        amrex::Array<MultiFab*,3> fmf{&xmom_update_2 , &ymom_update_2 , &zmom_update_2};

        // Interpolate from coarse faces to fine faces *only* on the coarse-fine boundary
        ifr->interp(fmf,cmf_const,0,1);

        amrex::Array<MultiFab*,3> cmf{&xmom_crse, &ymom_crse, &zmom_crse};

        int nGrow = 1;
        BoxArray fine_grids(cons_old.boxArray());

        // Interpolate from coarse faces on fine faces outside the fine region
        create_umac_grown(level, nGrow, fine_grids, geom, cmf, fmf, ref_ratio);
    }
#endif

    // **************************************************************************************
    // Convert updated momentum to updated velocity on faces after stage 2 and before stage 3
    // Make sure we have filled the ghost cells of cons_upd_2 and *mom_update_2 before this call!
    // **************************************************************************************
    MomentumToVelocity(xvel_new, yvel_new, zvel_new, cons_upd_2, xmom_update_2, ymom_update_2, zmom_update_2, 1, solverChoice);

    // **************************************************************************************
    // RK3 stage 3: Return update in the cons_new and [x,y,z]mom_new MultiFabs
    // **************************************************************************************
    // TODO: We won't need faceflux, edgeflux, and centflux when using the new code architecture
    RK3_stage(level, cons_upd_2, cons_new,
              xmom_update_2, ymom_update_2, zmom_update_2,
              xmom_new, ymom_new, zmom_new,
              xvel_new, yvel_new, zvel_new,
              source,
              faceflux,
              geom, dxp, dt,
              //ifr,
              solverChoice);
    // **************************************************************************************
    // RK3 stage 3: Define updates in the third RK stage
    // **************************************************************************************
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

        const Array4<Real>& momx_up2 = xmom_update_2.array(mfi);
        const Array4<Real>& momy_up2 = ymom_update_2.array(mfi);
        const Array4<Real>& momz_up2 = zmom_update_2.array(mfi);

        const Array4<Real>& momx_new = xmom_new.array(mfi);
        const Array4<Real>& momy_new = ymom_new.array(mfi);
        const Array4<Real>& momz_new = zmom_new.array(mfi);

        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            // At this point, the update is stored in cu_new ...
            cu_new(i,j,k,n) = (1./3.)*cu_old(i,j,k,n) + (2./3.)*(cu_up2(i,j,k,n) + cu_new(i,j,k,n));
            // Now cu_new holds the final solution
        });

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
    // TODO: Check if the order of applying BC on cell-centered state or face-centered mom makes any difference
    // TODO: Check if we need to apply BC on *mom_new, since we are using them in 'MomentumToVelocity'
    // PKJ comments: Even after applying BC on *mom_new, no noticeable difference in flow data was observed
    // Apply BC on updated momentum data on faces after stage 3 and before going to next time step
    // **************************************************************************************
    ///*
    xmom_new.FillBoundary(geom.periodicity());
    ymom_new.FillBoundary(geom.periodicity());
    zmom_new.FillBoundary(geom.periodicity());
    //*/

    // Apply BC on updated state data at cells after stage 3 and before going to next time step
    // **************************************************************************************
    cons_new.FillBoundary(geom.periodicity());

    amrex::Vector<MultiFab*> vars_new{&cons_new, &xmom_new, &ymom_new, &zmom_new};

    ERF::applyBCs(geom, vars_new);

    // **************************************************************************************
    // Convert updated momentum to updated velocity on faces after stage 3 and before going to next time step
    // **************************************************************************************
    if (level == 1) print_state(xmom_new,IntVect(8,8,8));
    MomentumToVelocity(xvel_new, yvel_new, zvel_new, cons_new, xmom_new, ymom_new, zmom_new, 0, solverChoice);

    // One final application of non-periodic BCs after all RK3 stages have been called, this time on velocities

    xvel_new.FillBoundary(geom.periodicity());
    yvel_new.FillBoundary(geom.periodicity());
    zvel_new.FillBoundary(geom.periodicity());
    amrex::Vector<MultiFab*> vars{&cons_new, &xvel_new, &yvel_new, &zvel_new};

    ERF::applyBCs(geom, vars);

}

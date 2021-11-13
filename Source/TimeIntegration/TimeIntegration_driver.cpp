#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>

#include <TimeIntegration.H>
#include <utils.H>

using namespace amrex;

// TODO: Check if the order of applying BC on cell-centered state or face-centered mom makes any difference

void erf_advance(int level,
                 MultiFab& cons_old,  MultiFab& cons_new,
                 MultiFab& xvel_old, MultiFab& yvel_old, MultiFab& zvel_old,
                 MultiFab& xvel_new, MultiFab& yvel_new, MultiFab& zvel_new,
                 MultiFab& xmom_crse, MultiFab& ymom_crse, MultiFab& zmom_crse,
                 MultiFab& source,
                 std::array< MultiFab, AMREX_SPACEDIM>& faceflux,
                 const amrex::Geometry crse_geom,
                 const amrex::Geometry fine_geom,
                 const amrex::IntVect ref_ratio,
                 const amrex::Real* dxp, const amrex::Real dt,
                 const amrex::Real time,
                 amrex::InterpFaceRegister* ifr,
                 const SolverChoice& solverChoice)
{
    BL_PROFILE_VAR("erf_advance()",erf_advance);

    int nvars = cons_old.nComp();

    // Determine the number of ghost cells depending on the spatial order
    // **************************************************************************************
    //TODO: Check if this is the only place to specify the number of ghost cells
    //TODO: Also explore how 'ngrow' should be related to the spatial_order
    int ngc = ComputeGhostCells(solverChoice.spatial_order);

    const BoxArray& ba            = cons_old.boxArray();
    const DistributionMapping& dm = cons_old.DistributionMap();

    // Initial momentum -- on faces
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > state_old;
    state_old.emplace_back(std::make_unique<amrex::MultiFab>(ba, dm, nvars, cons_old.nGrow())); // cons
    state_old.emplace_back(std::make_unique<amrex::MultiFab>(convert(ba,IntVect(1,0,0)), dm, 1, 1)); // xmom
    state_old.emplace_back(std::make_unique<amrex::MultiFab>(convert(ba,IntVect(0,1,0)), dm, 1, 1)); // ymom
    state_old.emplace_back(std::make_unique<amrex::MultiFab>(convert(ba,IntVect(0,0,1)), dm, 1, 1)); // zmom

    // Final momentum -- on faces
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > state_new;
    state_new.emplace_back(std::make_unique<amrex::MultiFab>(ba, dm, nvars, cons_old.nGrow())); // cons
    state_new.emplace_back(std::make_unique<amrex::MultiFab>(convert(ba,IntVect(1,0,0)), dm, 1, 1)); // xmom
    state_new.emplace_back(std::make_unique<amrex::MultiFab>(convert(ba,IntVect(0,1,0)), dm, 1, 1)); // ymom
    state_new.emplace_back(std::make_unique<amrex::MultiFab>(convert(ba,IntVect(0,0,1)), dm, 1, 1)); // zmom

    // **************************************************************************************
    // Prepare the old-time data for calling the integrator
    // **************************************************************************************

    // Apply BC on old state data at cells
    // **************************************************************************************
    cons_old.FillBoundary(fine_geom.periodicity());

    // We need to apply the boundary conditions here because we are converting from velocity to momentum
    //    which requires having set boundary conditions on density
    amrex::Vector<MultiFab*> vars_orig{&cons_old};
    ERF::applyBCs(fine_geom, vars_orig);

    MultiFab::Copy(*state_old[IntVar::cons], cons_old, 0, 0, cons_old.nComp(), cons_old.nGrow());

    // Convert old velocity available on faces to old momentum on faces to be used in time integration
    // **************************************************************************************
    VelocityToMomentum(xvel_old, yvel_old, zvel_old,
                       *state_old[IntVar::cons],
                       *state_old[IntVar::xmom],
                       *state_old[IntVar::ymom],
                       *state_old[IntVar::zmom],
                       solverChoice);

    // Apply BC on old momentum data on faces before integration
    // **************************************************************************************
    state_old[IntVar::xmom]->FillBoundary(fine_geom.periodicity());
    state_old[IntVar::ymom]->FillBoundary(fine_geom.periodicity());
    state_old[IntVar::zmom]->FillBoundary(fine_geom.periodicity());

    auto interpolate_coarse_fine_faces = [&](Vector<std::unique_ptr<MultiFab> >& S_data) {
        if (level > 0)
        {
            amrex::Array<const MultiFab*,3> cmf_const{&xmom_crse, &ymom_crse, &zmom_crse};
            amrex::Array<MultiFab*,3> fmf{S_data[IntVar::xmom].get(),
                                          S_data[IntVar::ymom].get(),
                                          S_data[IntVar::zmom].get()};

            // Interpolate from coarse faces to fine faces *only* on the coarse-fine boundary
            ifr->interp(fmf,cmf_const,0,1);

            amrex::Array<MultiFab*,3> cmf{&xmom_crse, &ymom_crse, &zmom_crse};

            int nGrow = 1;
            BoxArray fine_grids(cons_old.boxArray());

            // Interpolate from coarse faces on fine faces outside the fine region
            create_umac_grown(level, nGrow, fine_grids, crse_geom, fine_geom, cmf, fmf, ref_ratio);
        }
    };

    auto apply_bcs = [&](Vector<std::unique_ptr<MultiFab> >& S_data) {
        amrex::Vector<MultiFab*> state_p{S_data[IntVar::cons].get(),
                                         S_data[IntVar::xmom].get(),
                                         S_data[IntVar::ymom].get(),
                                         S_data[IntVar::zmom].get()};

        ERF::applyBCs(fine_geom, state_p);
    };

    // are these in the right order?
    interpolate_coarse_fine_faces(state_old);
    apply_bcs(state_old);

    // **************************************************************************************
    // Setup the integrator
    // **************************************************************************************
    TimeIntegrator<amrex::Vector<std::unique_ptr<amrex::MultiFab> > > integrator(state_old);

    auto rhs_fun = [&](Vector<std::unique_ptr<MultiFab> >& S_rhs, const Vector<std::unique_ptr<MultiFab> >& S_data, const Real time) {
        erf_rhs(level, S_rhs, S_data,
                source,
                faceflux,
                fine_geom, dxp, dt,
                ifr,
                solverChoice);
    };

    auto post_update_fun = [&](Vector<std::unique_ptr<MultiFab> >& S_data, const Real time) {
        // Apply BC on updated state and momentum data
        for (auto& mfp : S_data) {
            mfp->FillBoundary(fine_geom.periodicity());
        }

        // are these in the right order?
        apply_bcs(S_data);

        // TODO: we should interpolate coarse data in time first, so that this interplation
        // in space is at the correct time indicated by the `time` function argument.
        interpolate_coarse_fine_faces(S_data);
    };

    // define rhs and 'post update' utility function that is called after calculating
    // any state data (e.g. at RK stages or at the end of a timestep)
    integrator.set_rhs(rhs_fun);
    integrator.set_post_update(post_update_fun);

    // **************************************************************************************
    // Integrate for a single timestep
    // **************************************************************************************
    integrator.advance(state_old, state_new, time, dt);

    // **************************************************************************************
    // Convert updated momentum to updated velocity on faces after we have taken a timestep
    // **************************************************************************************
    MomentumToVelocity(xvel_new, yvel_new, zvel_new,
                       *state_new[IntVar::cons],
                       *state_new[IntVar::xmom],
                       *state_new[IntVar::ymom],
                       *state_new[IntVar::zmom],
                       0, solverChoice);

    // **************************************************************************************
    // Get the final cell centered variables after the step
    // (do this at the very end because its a swap not a copy)
    // **************************************************************************************
    std::swap(cons_new, *state_new[IntVar::cons]);

    // One final application of non-periodic BCs

    xvel_new.FillBoundary(fine_geom.periodicity());
    yvel_new.FillBoundary(fine_geom.periodicity());
    zvel_new.FillBoundary(fine_geom.periodicity());
    amrex::Vector<MultiFab*> vars{&cons_new, &xvel_new, &yvel_new, &zvel_new};

    ERF::applyBCs(fine_geom, vars);

}

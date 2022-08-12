
#include <ERF.H>
#include <SpatialStencils.H>

using namespace amrex;

// Advance a level by dt
// includes a recursive call for finer levels
void
ERF::timeStep (int lev, Real time, int iteration)
{
    if (regrid_int > 0)  // We may need to regrid
    {
        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static Vector<int> last_regrid_step(max_level+1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if
        // it was taken care of during a coarser regrid
        if (lev < max_level && istep[lev] > last_regrid_step[lev])
        {
            if (istep[lev] % regrid_int == 0)
            {
                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
                int old_finest = finest_level;
                regrid(lev, time);

                // mark that we have regridded this level already
                for (int k = lev; k <= finest_level; ++k) {
                    last_regrid_step[k] = istep[k];
                }

                // if there are newly created levels, set the time step
                for (int k = old_finest+1; k <= finest_level; ++k) {
                    dt[k] = dt[k-1] / MaxRefRatio(k-1);
                }
            }
        }
    }

    // Update what we call "old" and "new" time
    t_old[lev] = t_new[lev];
    t_new[lev] += dt[lev];

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
        amrex::Print() << "ADVANCE from time = " << t_old[lev] << " to " << t_new[lev]
                       << " with dt = " << dt[lev] << std::endl;
    }

    // Advance a single level for a single time step, updates flux registers
    Advance(lev, time, dt[lev], iteration, nsubsteps[lev]);

    ++istep[lev];

    if (Verbose())
    {
        amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
    }

    if (lev < finest_level)
    {
        // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev+1]; ++i)
        {
            timeStep(lev+1, time+(i-1)*dt[lev+1], i);
        }

        AverageDownTo(lev); // average lev+1 down to lev
    }
}

// advance a single level for a single time step, updates flux registers
void
ERF::Advance (int lev, Real time, Real dt_lev, int /*iteration*/, int /*ncycle*/)
{
    BL_PROFILE("ERF::Advance()");

    // We must swap the pointers so the previous step's "new" is now this step's "old"
    std::swap(vars_old[lev], vars_new[lev]);

    MultiFab& S_old = vars_old[lev][Vars::cons];
    MultiFab& S_new = vars_new[lev][Vars::cons];

    MultiFab& U_old = vars_old[lev][Vars::xvel];
    MultiFab& V_old = vars_old[lev][Vars::yvel];
    MultiFab& W_old = vars_old[lev][Vars::zvel];

    MultiFab& U_new = vars_new[lev][Vars::xvel];
    MultiFab& V_new = vars_new[lev][Vars::yvel];
    MultiFab& W_new = vars_new[lev][Vars::zvel];

    // We need to set these because otherwise in the first call to erf_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    U_new.setVal(1.e34,U_new.nGrowVect());
    V_new.setVal(1.e34,V_new.nGrowVect());
    W_new.setVal(1.e34,W_new.nGrowVect());

    FillPatch(lev, time, vars_old[lev]);

    MultiFab* S_crse;
    MultiFab rU_crse, rV_crse, rW_crse;

    if (lev > 0)
    {
        S_crse = &vars_old[lev-1][Vars::cons];

        MultiFab& U_crse = vars_old[lev-1][Vars::xvel];
        MultiFab& V_crse = vars_old[lev-1][Vars::yvel];
        MultiFab& W_crse = vars_old[lev-1][Vars::zvel];

        rU_crse.define(U_crse.boxArray(), U_crse.DistributionMap(), 1, U_crse.nGrow());
        rV_crse.define(V_crse.boxArray(), V_crse.DistributionMap(), 1, V_crse.nGrow());
        rW_crse.define(W_crse.boxArray(), W_crse.DistributionMap(), 1, W_crse.nGrow());

        VelocityToMomentum(U_crse, U_crse.nGrowVect(),
                           V_crse, V_crse.nGrowVect(),
                           W_crse, W_crse.nGrowVect(),
                          *S_crse,rU_crse,rV_crse,rW_crse);
    }

    // configure ABLMost params if used MostWall boundary condition
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] != ERF_BC::MOST) {
        if (m_most) setupABLMost(lev);
    }

    // Do an error check
    if (solverChoice.pbl_type == PBLType::MYNN25 &&
        phys_bc_type[Orientation(Direction::z,Orientation::low)] != ERF_BC::MOST) {
        amrex::Error("Must use MOST BC for MYNN2.5 PBL model");
    }

    const auto& local_ref_ratio = (lev > 0) ? ref_ratio[lev-1] : IntVect(1,1,1);

    InterpFaceRegister ifr;
    if (lev > 0)
    {
        ifr.define(S_old.boxArray(), S_old.DistributionMap(), Geom(lev), local_ref_ratio);
    }

    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    int nvars = S_old.nComp();

    // Place-holder for source array -- for now just set to 0
    MultiFab source(ba,dm,nvars,1);
    source.setVal(0.0);

    // These are the actual fluxes we will use to fill the flux registers
    std::array< MultiFab, AMREX_SPACEDIM > flux;

    flux[0].define(convert(ba,IntVect(1,0,0)), dm, nvars, 1);
    flux[1].define(convert(ba,IntVect(0,1,0)), dm, nvars, 1);
    flux[2].define(convert(ba,IntVect(0,0,1)), dm, nvars, 1);

    flux[0].setVal(0.);
    flux[1].setVal(0.);
    flux[2].setVal(0.);

    // We don't need to call FillPatch on cons_mf because we have fillpatch'ed S_old above
    MultiFab cons_mf(ba,dm,nvars,S_old.nGrowVect());
    MultiFab::Copy(cons_mf,S_old,0,0,S_old.nComp(),S_old.nGrowVect());

    MultiFab* r0 = &dens_hse[lev];
    MultiFab* p0 = &pres_hse[lev];

    amrex::Real* dptr_rayleigh_tau      = solverChoice.use_rayleigh_damping ? d_rayleigh_tau[lev].data() : nullptr;
    amrex::Real* dptr_rayleigh_ubar     = solverChoice.use_rayleigh_damping ? d_rayleigh_ubar[lev].data() : nullptr;
    amrex::Real* dptr_rayleigh_vbar     = solverChoice.use_rayleigh_damping ? d_rayleigh_vbar[lev].data() : nullptr;
    amrex::Real* dptr_rayleigh_thetabar = solverChoice.use_rayleigh_damping ? d_rayleigh_thetabar[lev].data() : nullptr;

    // *****************************************************************
    // Update the cell-centered state and face-based velocity using
    // a time integrator.
    // Inputs:
    //          S_old    (state on cell centers)
    //          U_old    (x-velocity on x-faces)
    //          V_old    (y-velocity on y-faces)
    //          W_old    (z-velocity on z-faces)
    //          source   (source term on cell centers)
    // Outputs:
    //          S_new    (state on cell centers)
    //          U_new    (x-velocity on x-faces)
    //          V_new    (y-velocity on y-faces)
    //          W_new    (z-velocity on z-faces)
    // *****************************************************************

    erf_advance(lev,
                cons_mf, S_new,
                U_old, V_old, W_old,
                U_new, V_new, W_new,
                rU_crse, rV_crse, rW_crse,
                source, flux,
                Geom(lev), dt_lev, time, &ifr,
                r0, p0,
                dptr_rayleigh_tau, dptr_rayleigh_ubar,
                dptr_rayleigh_vbar, dptr_rayleigh_thetabar);

    // *****************************************************************
    // Now fill the flux registers with the fluxes so we can reflux later
    // *****************************************************************
    const auto& dx = Geom(lev).CellSize();

    if (finest_level > 0 && do_reflux)
    {
        if (lev < finest_level) {
            flux_registers[lev+1]->setVal(0.0);
            flux_registers[lev+1]->CrseInit(flux[0],0,0,0,nvars,-dx[1]*dx[2],amrex::FluxRegister::ADD);
            flux_registers[lev+1]->CrseInit(flux[1],1,0,0,nvars,-dx[0]*dx[2],amrex::FluxRegister::ADD);
            flux_registers[lev+1]->CrseInit(flux[2],2,0,0,nvars,-dx[0]*dx[1],amrex::FluxRegister::ADD);
        }

        if (lev > 0) {
          flux_registers[lev]->FineAdd(flux[0], 0, 0, 0, nvars, dx[1]*dx[2]);
          flux_registers[lev]->FineAdd(flux[1], 1, 0, 0, nvars, dx[0]*dx[2]);
          flux_registers[lev]->FineAdd(flux[2], 2, 0, 0, nvars, dx[0]*dx[1]);
        }
    }
}

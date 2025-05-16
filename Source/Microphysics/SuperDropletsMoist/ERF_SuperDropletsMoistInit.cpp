#include <sys/time.h>
#include "ERF_SuperDropletsMoist.H"
#include "ERF_MaterialProperties.H"

#ifdef ERF_USE_PARTICLES

/*! Define the super-droplet moisture model parameters from provided inputs */
void SuperDropletsMoist::Define (SolverChoice& a_sc /*!< Solver choices */)
{
    m_fac_cond = lcond / a_sc.c_p;
}

/*! Read inputs from file */
void SuperDropletsMoist::readInputs ()
{
    BL_PROFILE("SuperDropletsMoist::readInputs");
    ParmParse pp(m_name);

    // include phase change in super-droplet dynamics?
    m_flag_phase_change = true; //default
    pp.query("include_phase_change", m_flag_phase_change);
    // include advection in super-droplet dynamics?
    m_flag_advection = true; //default
    pp.query("include_advection", m_flag_advection);
    // include coalescence in super-droplet dynamics?
    m_flag_coalescence = true; //default
    pp.query("include_coalescence", m_flag_coalescence);

    // initial distribution type
    m_init_type = SupDropInit::init_uniform;
    pp.query("initial_distribution_type", m_init_type);

    // minimum radius for rain
    m_r_rain = 4.0e-5; // 40 micrometers
    pp.query("radius_raindrop", m_r_rain);

    // whether to run in kinematic mode
    m_kinematic_mode = false;
    pp.query("kinematic_mode", m_kinematic_mode);

    // simulation dimensionality
    m_dimensionality = SDMSimulationDim::three_d;
    pp.query("dimensionality", m_dimensionality);

    // recycle super-droplets
    m_recycle_particles = false;
    pp.query("recycle_particles", m_recycle_particles);


    // get vapour/condensate species names
    m_species.clear();
    // add water
    m_idx_w = m_species.size();
    m_species.push_back(Species::Name::H2O);
    // add other species
    std::string species_input = "species";
    if (pp.contains(species_input.c_str())) {
        int num_species = pp.countval(species_input.c_str());
        Species::Name sp_name;
        for (int i = 0; i < num_species; i++) {
            pp.get(species_input.c_str(), sp_name, i);
            m_species.push_back(sp_name);
        }
    }
    m_num_species = m_species.size();
    m_qstate_nonmoist_size = (m_num_species-1)*2; // qv, qc for each

    // get aerosol names
    m_aerosols.clear();
    std::string aerosol_input = "aerosols";
    if (pp.contains(aerosol_input.c_str())) {
        int num_aerosols = pp.countval(aerosol_input.c_str());
        Species::Name aero_name;
        for (int i = 0; i < num_aerosols; i++) {
            pp.get(aerosol_input.c_str(), aero_name, i);
            m_aerosols.push_back(aero_name);
        }
    }
    m_num_aerosols = m_aerosols.size();

    // number of time steps between writing distribution  diagnostics to file
    m_diagnostics_iter = 1; //default
    pp.query("diagnostics_interval", m_diagnostics_iter);

    // number of substeps for phase change process
    m_num_substeps_phase_change = 1; //default
    pp.query("num_substeps_phase_change", m_num_substeps_phase_change);

    // let superdroplets relax to a physically correct size at initialization?
    m_init_phase_change = false; //default
    pp.query("initial_phase_change_relaxation", m_init_phase_change);
    // time (in seconds) of initial relaxation
    m_init_phase_change_time = 10.0; //default
    pp.query("initial_phase_change_relaxation_time", m_init_phase_change_time);

    return;
}

/*! Initializes the super-droplet moisture model: allocates the moisture model
    variables Multifabs and the super-droplet particle container. */
void SuperDropletsMoist::Init ( const MultiFab&   a_cons_vars,  /*!< Conserved variables */
                                const BoxArray&,                /*!< Grids */
                                const Geometry&   a_geom,       /*!< Computational domain */
                                const Real&       a_dt,         /*!< Timestep */
                                MFPtr&,
                                MFPtr& )
{
    BL_PROFILE("SuperDropletsMoist::Init()");
    m_dt = a_dt;
    m_geom = a_geom;

    m_mic_var_map = {   MicVar_SD::q_t,
                        MicVar_SD::q_v,
                        MicVar_SD::q_c,
                        MicVar_SD::dqcdt,
                        MicVar_SD::q_r,
                        MicVar_SD::rh,
                        MicVar_SD::rain_accum };
    AMREX_ALWAYS_ASSERT(m_qmoist_size == m_mic_var_map.size());

    /* allocate microphysics multifabs */
    m_mic_fab_vars.resize( MicVar_SD::NumVars
                          +(m_num_species-1)*MicVar_SD_Species::NumVars);
    for (auto i(0); i < m_mic_fab_vars.size(); i++) {
      m_mic_fab_vars[i] = std::make_shared<MultiFab> ( a_cons_vars.boxArray(),
                                                       a_cons_vars.DistributionMap(),
                                                       1,
                                                       a_cons_vars.nGrowVect() );
      m_mic_fab_vars[i]->setVal(0.0);
    }

    /* create the super-droplet particle container */
    m_super_droplets = new SuperDropletPC ( a_geom,
                                            a_cons_vars.DistributionMap(),
                                            a_cons_vars.boxArray(),
                                            m_species,
                                            m_aerosols,
                                            m_name );

    amrex::Print() << "SuperDropletsMoist:\n"
                   << "    diagnostics_interval: " << m_diagnostics_iter << "\n"
                   << "    cloud/rain radius: " << m_r_rain << " [m]\n"
                   << "    kinematic mode: " << (m_kinematic_mode?"true":"false") << "\n"
                   << "    dimensionality: " << amrex::getEnumNameString(m_dimensionality)  << "\n"
                   << "    include phase change: "
                   << (m_flag_phase_change ? "true" : "false") << "\n"
                   << "    include particle advection: "
                   << (m_flag_advection ? "true" : "false") << "\n"
                   << "    include coalescence: "
                   << (m_flag_coalescence ? "true" : "false") << "\n"
                   << "    Recycle particles: " << (m_recycle_particles ? "true" : "false") << "\n"
                   << "    number of substeps (phase change): " << m_num_substeps_phase_change << "\n"
                   << "    initial phase change relaxation: "
                   << (m_init_phase_change ? "true" : "false") << "\n";
    if (m_init_phase_change) {
        amrex::Print()  << "    initial phase change relaxation time: "
                        << m_init_phase_change_time << "\n";
    }

}

/*! Initializes particles in the super-droplet moisture model: if not a restart run, initialize the particles
 *  in the particle container for the super-droplets. */
void SuperDropletsMoist::InitParticles ( MFPtr& a_z_phys_nd /*!< terrain */)
{
    BL_PROFILE("SuperDropletsMoist::InitParticles()");

    if (m_init_type == SuperDropletsMoistInitializations::init_rhoc) {
        /* The conserved variables are not set up yet; the initial condensate
           density is not available. So, just initialize with a uniform distribution
           for now; set the radius and multiplicity from condensate density when
           Update_Micro_Vars() is called for the first time. */
        m_super_droplets->InitializeParticles(a_z_phys_nd);
    } else {
        m_super_droplets->InitializeParticles(a_z_phys_nd);
        amrex::Print() << "Initialized "
                       << m_super_droplets->NumSuperDroplets()
                       << " super-droplets representing "
                       << m_super_droplets->TotalNumberOfParticles()
                       << " particles in super-droplets moisture model.\n";
    }
}

/*! Restarts particles in the super-droplet moisture model */
void SuperDropletsMoist::RestartParticles ( ParGDBBase*, const std::string& a_fname )
{
    BL_PROFILE("SuperDropletsMoist::RestartParticles()");

    amrex::Print() << "Reading in " << m_name << " particle data from restart file.\n";

    struct timeval total_start, total_end;
    gettimeofday(&total_start, NULL);
    m_super_droplets->Restart(a_fname, m_name);
    m_super_droplets->Redistribute();
    gettimeofday(&total_end,NULL);
    long long total_wtime;
    total_wtime = (   (total_end.tv_sec   * 1000000 + total_end.tv_usec  )
                   -  (total_start.tv_sec * 1000000 + total_start.tv_usec) );
    Real total_wtime_sec = (double) total_wtime / 1000000.0;
    ParallelDescriptor::ReduceRealMax( &total_wtime_sec,
                                       1,
                                       ParallelDescriptor::IOProcessorNumber() );

    amrex::Print() << "Restarted "
                   << m_super_droplets->NumSuperDroplets()
                   << " super-droplets representing "
                   << m_super_droplets->TotalNumberOfParticles()
                   << " particles in super-droplets moisture model "
                   << "(" << total_wtime_sec << " seconds).\n";
}

/*! Finish any initializations that depend on the conserved state variables
    that were not available during Init(); also, overwrite the rhoq2 component
    with the quantity computed from the this model. */
void SuperDropletsMoist::FinishInit (const int& /* a_lev */,
                                     MultiFab& a_cons_vars, /*!< Conserved variables */
                                     const Vector<MFPtr>& a_z_phys_nd /*!< terrain */)
{
    m_super_droplets->DensityScaling(*(m_mic_fab_vars[MicVar_SD::rho]));

    if (m_init_type == SuperDropletsMoistInitializations::init_rhoc) {

        /* initial super-droplets attributes computed from condensate mass density */
        MultiFab rho_c ( m_mic_fab_vars[MicVar_SD::q_c]->boxArray(),
                         m_mic_fab_vars[MicVar_SD::q_c]->DistributionMap(),
                         1,
                         m_mic_fab_vars[MicVar_SD::q_c]->nGrowVect() );
        MultiFab::Copy( rho_c, *m_mic_fab_vars[MicVar_SD::q_c], 0, 0, 1, rho_c.nGrowVect() );
        ratioToDensity(rho_c);

        m_super_droplets->SetAttributes(rho_c);
        amrex::Print() << "Initialized "
                       << m_super_droplets->NumSuperDroplets()
                       << " super-droplets representing "
                       << m_super_droplets->TotalNumberOfParticles()
                       << " particles in super-droplets moisture model.\n";

    } else {

        /* call the phase change function so that the super-droplets "relax" to their
         * physical size corresponding to the initial flow */
        if (m_flag_phase_change && m_init_phase_change) {
            phaseChange(m_init_phase_change_time, a_z_phys_nd, true);
        }
    }

    computeQcQrWater();
    computeQtWater();

    for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();
        auto states_arr = a_cons_vars.array(mfi);
        auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->array(mfi);
        auto q_r_arr = m_mic_fab_vars[MicVar_SD::q_r]->array(mfi);
        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            states_arr(i,j,k,RhoQ2_comp) = states_arr(i,j,k,Rho_comp)*q_c_arr(i,j,k);
            states_arr(i,j,k,RhoQ3_comp) = states_arr(i,j,k,Rho_comp)*q_r_arr(i,j,k);
        });
    }

    computeQcSpecies();
    computeQtSpecies();

    for (int is = 1; is < m_num_species; is++) {
        for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {
            const auto& box = mfi.tilebox();
            auto states_arr = a_cons_vars.array(mfi);
            auto q_c_arr = m_mic_fab_vars[s_qc_idx(is)]->array(mfi);
            auto qc_comp = q_qc_idx(is);
            ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                states_arr(i,j,k,qc_comp) = states_arr(i,j,k,Rho_comp)*q_c_arr(i,j,k);
            });
        }
    }

    m_super_droplets->Diagnostics(-1, 0.0, (m_diagnostics_iter>0));

    return;
}

#endif

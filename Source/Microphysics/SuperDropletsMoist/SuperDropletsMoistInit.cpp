#include "SuperDropletsMoist.H"
#include "MaterialProperties.H"

#ifdef ERF_USE_PARTICLES

/*! Define the super-droplet moisture model parameters from provided inputs */
void SuperDropletsMoist::Define (SolverChoice& /* a_sc */)
{
    readInputs();
}

/*! Read inputs from file */
void SuperDropletsMoist::readInputs ()
{
    BL_PROFILE("SuperDropletsMoist::readInputs");
    ParmParse pp(m_name);

    m_flag_phase_change = true;
    m_flag_advection = true;
    m_flag_coalescence = true;
    pp.query("include_phase_change", m_flag_phase_change);
    pp.query("include_advection", m_flag_advection);
    pp.query("include_coalescence", m_flag_coalescence);

    m_init_type = SupDropInit::init_uniform;
    pp.query("initial_distribution_type", m_init_type);

    m_aerosols.clear();
    std::string aerosol_input = "aerosols";
    if (pp.contains(aerosol_input.c_str())) {
        int num_aerosols = pp.countval(aerosol_input.c_str());
        std::string aero_name;
        for (int i = 0; i < num_aerosols; i++) {
            pp.get(aerosol_input.c_str(), aero_name, i);
            m_aerosols.push_back(aero_name);
        }
    }

    m_diagnostics_iter = 1;
    pp.query("diagnostics_interval", m_diagnostics_iter);

    return;
}

/*! Initializes the super-droplet moisture model: allocates the moisture model
    variables Multifabs and the super-droplet particle container. */
void SuperDropletsMoist::Init ( const MultiFab&   a_cons_vars,  /*!< Conserved variables */
                                const BoxArray&,                /*!< Grids */
                                const Geometry&   a_geom,       /*!< Computational domain */
                                const Real&       a_dt,         /*!< Timestep */
                                MFPtr&            a_z_phys_nd,   /*!< terrain */
                                MFPtr& )
{
    BL_PROFILE("SuperDropletsMoist::Init()");
    m_dt = a_dt;
    m_geom = a_geom;

    m_mic_var_map.resize(m_qmoist_size);
    m_mic_var_map = {MicVar_SD::q_t, MicVar_SD::q_v, MicVar_SD::q_c};

    /* allocate microphysics multifabs */
    for (auto i(0); i < MicVar_SD::NumVars; i++) {
      m_mic_fab_vars[i] = std::make_shared<MultiFab> ( a_cons_vars.boxArray(),
                                                       a_cons_vars.DistributionMap(),
                                                       1,
                                                       a_cons_vars.nGrowVect() );
      m_mic_fab_vars[i]->setVal(0.0);
    }

    /* create the super-droplet particle container */
    std::string vapour_mat = MaterialNames::h2o; // water
    m_super_droplets = new SuperDropletPC ( a_geom,
                                            a_cons_vars.DistributionMap(),
                                            a_cons_vars.boxArray(),
                                            vapour_mat, m_aerosols,
                                            m_name );

    if (m_init_type == SuperDropletsMoistInitializations::init_rhoc) {
        /* The conserved variables are not set up yet; the initial condensate
           density is not available. So, just initialize with a uniform distribution
           for now; set the radius and multiplicity from condensate density when
           Update_Micro_Vars() is called for the first time. */
        m_super_droplets->InitializeParticles( SupDropInit::init_uniform,
                                               a_z_phys_nd );
    } else {
        m_super_droplets->InitializeParticles(a_z_phys_nd);
        amrex::Print() << "Initialized "
                       << m_super_droplets->NumSuperDroplets()
                       << " super-droplets representing "
                       << m_super_droplets->TotalNumberOfParticles()
                       << " particles in super-droplets moisture model.\n";
    }

    m_super_droplets->Diagnostics();

    amrex::Print() << "SuperDropletsMoist:\n"
                   << "    diagnostics_interval: " << m_diagnostics_iter << "\n"
                   << "    include phase change: "
                   << (m_flag_phase_change ? "true" : "false") << "\n"
                   << "    include particle advection: "
                   << (m_flag_advection ? "true" : "false") << "\n"
                   << "    include coalescence: "
                   << (m_flag_coalescence ? "true" : "false") << "\n";

}

/*! Finish any initializations that depend on the conserved state variables
    that were not available during Init(); also, overwrite the rhoq2 component
    with the quantity computed from the this model. */
void SuperDropletsMoist::FinishInit (MultiFab& a_cons_vars /*!< Conserved variables */)
{
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
    }

    computeQc();
    computeQt();

    for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();
        auto states_arr = a_cons_vars.array(mfi);
        auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->array(mfi);
        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            states_arr(i,j,k,RhoQ2_comp) = states_arr(i,j,k,Rho_comp)*q_c_arr(i,j,k);
        });
    }

    return;
}


#endif

#include "SuperDropletsMoist.H"

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

    m_init_type = ERFParticleInitializations::init_default;
    pp.query("initial_distribution_type", m_init_type);

    return;
}

/*! Initializes the super-droplet moisture model: allocates the moisture model
    variables Multifabs and the super-droplet particle container. */
void SuperDropletsMoist::Init ( const MultiFab&   a_cons_vars,  /*!< Conserved variables */
                                const BoxArray&,                /*!< Grids */
                                const Geometry&   a_geom,       /*!< Computational domain */
                                const Real&       a_dt          /*!< Timestep */ )
{
    BL_PROFILE("SuperDropletsMoist::Init()");
    m_dt = a_dt;

    m_mic_var_map.resize(m_qmoist_size);
    m_mic_var_map = {MicVar_SD::rho_v, MicVar_SD::rho_c};

    /* allocate vapour and condensate multifabs */
    for (auto i(0); i < MicVar_SD::NumVars; i++) {
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
                                            m_name );

    if (m_init_type == SuperDropletsMoistInitializations::init_rhoc) {
        /* The conserved variables are not set up yet; the initial condensate
           density is not available. So, just initialize with a uniform distribution
           for now; set the radius and multiplicity from condensate density when
           Update_Micro_Vars() is called for the first time. */
        /*TODO: add terrain*/
        m_super_droplets->InitializeParticles(SuperDropletInitializations::init_uniform);
    } else {
        m_super_droplets->InitializeParticles(/*TODO: add terrain*/);
        amrex::Print() << "Initialized "
                       << m_super_droplets->NumSuperDroplets()
                       << " super-droplets representing "
                       << m_super_droplets->TotalNumberOfParticles()
                       << " particles in super-droplets moisture model.\n";
    }

}

/*! Finish any initializations that depend on the conserved state variables
    that were not available during Init() */
void SuperDropletsMoist::FinishInit ()
{
    if (m_init_type == SuperDropletsMoistInitializations::init_rhoc) {
        /* initial super-droplets attributes computed from condensate mass density */
        m_super_droplets->SetAttributes(*(m_mic_fab_vars[MicVar_SD::rho_c]));
        amrex::Print() << "Initialized "
                       << m_super_droplets->NumSuperDroplets()
                       << " super-droplets representing "
                       << m_super_droplets->TotalNumberOfParticles()
                       << " particles in super-droplets moisture model.\n";
    }
    return;
}


#endif

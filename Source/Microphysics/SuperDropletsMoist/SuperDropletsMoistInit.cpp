#include "SuperDropletsMoist.H"

#ifdef ERF_USE_PARTICLES

/*! Define the super-droplet moisture model parameters from provided inputs */
void SuperDropletsMoist::Define (SolverChoice& /* a_sc */)
{
}

/*! Initializes the super-droplet moisture model: allocates the moisture model
    variables Multifabs and the super-droplet particle container. */
void SuperDropletsMoist::Init ( const MultiFab&   a_cons_vars,  /*!< Conserved variables */
                                const BoxArray&,                /*!< Grids */
                                const Geometry&   a_geom,       /*!< Computational domain */
                                const Real&       a_dt          /*!< Timestep */ )
{
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
    m_super_droplets->InitializeParticles(/*TODO: add terrain*/);
    amrex::Print() << "Initialized " << m_super_droplets->TotalNumberOfParticles()
                   << " particles in super-droplets moisture model.\n";
}

#endif

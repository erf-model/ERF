#ifndef _WIN32
#include <sys/time.h>
#endif
#include "ERF_SuperDropletsMoist.H"
#include "ERF_MaterialProperties.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! \brief Define the super-droplet moisture model parameters from provided inputs
 *
 * This function initializes basic model parameters from the solver choices object.
 * Currently it only sets the specific heat capacity (Cp) from the solver choices.
 *
 * \param[in] a_sc Solver choices containing model configuration parameters
 */
void SuperDropletsMoist::Define (SolverChoice& a_sc)
{
    BL_PROFILE("SuperDropletsMoist::Define()");
    m_Cp = a_sc.c_p;
}

/*! \brief Read model configuration parameters from input file
 *
 * This function reads all configuration parameters for the SuperDropletsMoist
 * model from the input file. Parameters include:
 * - Model flags (phase change, advection, coalescence)
 * - Initial distribution type
 * - Rain threshold radius
 * - Model modes (kinematic mode, dimensionality)
 * - Particle recycling options
 * - Species and aerosol definitions
 * - Diagnostic and output options
 * - Phase change parameters
 */
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
    m_init_type = SDMoistInit::uniform;
    pp.query("distribution_type", m_init_type);

    // minimum radius for rain
    m_r_rain = Real(4.0e-5); // 40 micrometers
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
    m_init_phase_change_time = Real(10.0); //default
    pp.query("initial_phase_change_relaxation_time", m_init_phase_change_time);

    return;
}

/*! \brief Initialize the super-droplet moisture model
 *
 * Allocates the moisture model variable MultiFabs and creates the
 * super-droplet particle container. This function sets up:
 * one The mapping between moisture variable indices and internal arrays
 * two MultiFabs for all moisture model variables
 * three The SuperDropletPC particle container
 *
 * After initialization, it prints configuration summary to output.
 *
 * \param[in] a_cons_vars Conserved variables MultiFab
 * \param[in] a_geom Geometry information for computational domain
 * \param[in] a_dt Timestep size
 */
void SuperDropletsMoist::Init ( const MultiFab&   a_cons_vars,
                                const BoxArray&,
                                const Geometry&   a_geom,
                                const Real&       a_dt,
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
    m_mic_fab_vars.resize(  MicVar_SD::NumVars
                          + (m_num_species-1) * MicVar_SD_Species::NumVars
                          + m_num_aerosols * MicVar_SD_Aerosols::NumVars );
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
                                            m_dt,
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

/*! \brief Initialize particles in the super-droplet moisture model
 *
 * If this is not a restart run, this function initializes the particles
 * in the particle container for the super-droplets. The initialization
 * depends on the configured initialization type:
 * - For condensate_density: Initially creates particles with uniform distribution,
 *   which will be updated later during FinishInit() based on condensate density
 * - For other types: Initializes particles directly with the configured distribution
 *
 * \param[in] a_z_phys_nd MultiFab containing terrain height information
 */
void SuperDropletsMoist::InitParticles ( MFPtr& a_z_phys_nd )
{
    BL_PROFILE("SuperDropletsMoist::InitParticles()");

    if (m_init_type == SDMoistInit::condensate_density) {
        /* The conserved variables are not set up yet; the initial condensate
           density is not available. So, just initialize with a uniform distribution
           for now; set the radius and multiplicity from condensate density when
           Update_Micro_Vars() is called for the first time. */
        m_super_droplets->InitializeParticles(zero, a_z_phys_nd);
    } else {
        m_super_droplets->InitializeParticles(zero, a_z_phys_nd);
        amrex::Print() << "Initialized "
                       << m_super_droplets->NumSuperDroplets()
                       << " super-droplets representing "
                       << m_super_droplets->TotalNumberOfParticles()
                       << " particles in super-droplets moisture model.\n";
    }
}

/*! \brief Restart particles in the super-droplet moisture model from checkpoint file
 *
 * This function restarts superdroplet particles from a checkpoint file.
 * It performs the following operations:
 * one Reads particle data from the specified restart file
 * two Redistributes particles to appropriate processors/grids
 * three Measures and reports the time taken to perform the restart
 * Real(4.) Outputs statistics about the restarted particle population
 *
 * \param[in] a_gdb Unused particle grid database pointer
 * \param[in] a_fname File name for the checkpoint file to restart from
 */
void SuperDropletsMoist::RestartParticles ( ParGDBBase* /* a_gdb */, const std::string& a_fname )
{
    BL_PROFILE("SuperDropletsMoist::RestartParticles()");

    amrex::Print() << "Reading in " << m_name << " particle data from restart file.\n";

#ifndef _WIN32
    struct timeval total_start, total_end;
    gettimeofday(&total_start, NULL);
#endif
    m_super_droplets->Restart(a_fname, m_name);
    m_super_droplets->Redistribute();
#ifndef _WIN32
    gettimeofday(&total_end,NULL);
    long long total_wtime;
    total_wtime = (   (total_end.tv_sec   * 1000000 + total_end.tv_usec  )
                   -  (total_start.tv_sec * 1000000 + total_start.tv_usec) );
    Real total_wtime_sec = (double) total_wtime / Real(1000000.0);
    ParallelDescriptor::ReduceRealMax( &total_wtime_sec,
                                       1,
                                       ParallelDescriptor::IOProcessorNumber() );
#else
    Real total_wtime_sec = zero;
#endif

    amrex::Print() << "Restarted "
                   << m_super_droplets->NumSuperDroplets()
                   << " super-droplets representing "
                   << m_super_droplets->TotalNumberOfParticles()
                   << " particles in super-droplets moisture model "
                   << "(" << total_wtime_sec << " seconds).\n";
}

/*! \brief Complete initialization using now-available state variables
 *
 * This function finalizes initialization steps that depend on conserved state
 * variables that were not available during Init(). It performs:
 * one Particle density scaling based on air density
 * two For condensate_density initialization type: sets particle attributes from condensate density
 * three For other initialization types: optionally performs initial phase change relaxation
 * Real(4.) Computes cloud/rain water and total water content for all species
 * Real(5.) Updates the rhoq2 component in conserved variables with computed cloud water
 * Real(6.) Runs initial diagnostics for superdroplets
 *
 * \param[in] a_lev Unused AMR level parameter
 * \param[in,out] a_cons_vars Conserved variables MultiFab to be updated
 * \param[in] a_z_phys_nd Vector of MultiFabs containing terrain information
 */
void SuperDropletsMoist::FinishInit (const int& /* a_lev */,
                                     MultiFab& a_cons_vars,
                                     const Vector<MFPtr>& a_z_phys_nd)
{
    BL_PROFILE("SuperDropletsMoist::FinishInit()");
    m_super_droplets->DensityScaling(*(m_mic_fab_vars[MicVar_SD::rho]));

    if (m_init_type == SDMoistInit::condensate_density) {

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

    m_super_droplets->Diagnostics(-1, zero, (m_diagnostics_iter>0));

    return;
}

#endif

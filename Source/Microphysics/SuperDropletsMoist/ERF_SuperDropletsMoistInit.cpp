#include <algorithm>

#include "ERF_SuperDropletsMoist.H"
#include "ERF_MaterialProperties.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! \brief Define the super-droplet moisture model parameters from provided inputs
 *
 * This function initializes basic model parameters from the solver choices object.
 *
 * \param[in] a_sc Solver choices containing model configuration parameters
 */
void SuperDropletsMoist::Define (SolverChoice& a_sc)
{
    BL_PROFILE("SuperDropletsMoist::Define()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(a_sc.c_p) && a_sc.c_p > Real(0),
                                     "SDM requires finite positive c_p");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_sc.rdOcp == R_d / a_sc.c_p &&
                                     std::isfinite(a_sc.rdOcp) && a_sc.rdOcp > Real(0),
                                     "SDM requires rdOcp = R_d/c_p");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!a_sc.anelastic.empty(),
                                     "SDM requires per-level anelastic configuration");
    for (int mode : a_sc.anelastic) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mode == 0 || mode == 1,
                                         "SDM anelastic selectors must be 0 or 1");
    }
    m_Cp = a_sc.c_p;
    m_rdOcp = a_sc.rdOcp;
    m_anelastic = a_sc.anelastic;
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
    pp.queryAdd("include_phase_change", m_flag_phase_change);
    // include advection in super-droplet dynamics?
    m_flag_advection = true; //default
    pp.queryAdd("include_advection", m_flag_advection);
    // include coalescence in super-droplet dynamics?
    m_flag_coalescence = true; //default
    pp.queryAdd("include_coalescence", m_flag_coalescence);

    // include cold processes?
    m_with_ice = true;
    pp.queryAdd("include_cold_processes", m_with_ice);

    // initial distribution type
    m_init_type = SDMoistInit::uniform;
    pp.queryAdd("distribution_type", m_init_type);

    // minimum radius for rain
    m_r_rain = Real(4.0e-5); // 40 micrometers
    pp.queryAdd("radius_raindrop", m_r_rain);

    m_rime_ratio = 0.3;
    pp.queryAdd("rime_mass_ratio", m_rime_ratio);

    // whether to run in kinematic mode
    m_kinematic_mode = false;
    pp.queryAdd("kinematic_mode", m_kinematic_mode);

    // simulation dimensionality
    m_dimensionality = SDMSimulationDim::three_d;
    pp.queryAdd("dimensionality", m_dimensionality);

    // recycle super-droplets
    m_recycle_particles = false;
    pp.queryAdd("recycle_particles", m_recycle_particles);


    // get vapour/condensate species names
    m_species.clear();
    // add water
    m_idx_w = static_cast<int>(m_species.size());
    m_species.push_back(Species::Name::H2O);
    // add ice
    if (m_with_ice) {
        m_idx_i = m_species.size();
        m_species.push_back(Species::Name::ice);
    }
    // add other species
    m_istart_sp = m_species.size();
    std::string species_input = "species";
    if (pp.contains(species_input.c_str())) {
        int num_species = pp.countval(species_input.c_str());
        Species::Name sp_name;
        for (int i = 0; i < num_species; i++) {
            pp.get(species_input.c_str(), sp_name, i);
            // water and, with cold processes, ice are already present
            if (std::find(m_species.begin(), m_species.end(), sp_name) != m_species.end()) {
                amrex::Abort("SuperDropletsMoist: species "+amrex::getEnumNameString(sp_name)
                             +" is already modeled; remove it from super_droplets_moisture.species");
            }
            m_species.push_back(sp_name);
        }
    }
    m_num_species = static_cast<int>(m_species.size());
    m_num_nonmoist_sp = m_num_species - m_istart_sp;
    m_qstate_nonmoist_size = (m_num_nonmoist_sp)*2; // qv, qc for each

    if (m_num_nonmoist_sp > 5) {
        amrex::Abort("SuperDropletsMoist: cannot run with more than 5 non-moist species (see ERF_IndexDefines.H, definitions of NBCVAR_max and NVAR_max).");
    }

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
    m_num_aerosols = static_cast<int>(m_aerosols.size());

    // number of time steps between writing distribution  diagnostics to file
    m_diagnostics_iter = 1; //default
    pp.queryAdd("diagnostics_interval", m_diagnostics_iter);

    // number of substeps for phase change process
    m_num_substeps_phase_change = 1; //default
    pp.queryAdd("num_substeps_phase_change", m_num_substeps_phase_change);

    // let superdroplets relax to a physically correct size at initialization?
    m_init_phase_change = false; //default
    pp.queryAdd("initial_phase_change_relaxation", m_init_phase_change);
    // time (in seconds) of initial relaxation
    m_init_phase_change_time = Real(10.0); //default
    pp.queryAdd("initial_phase_change_relaxation_time", m_init_phase_change_time);

    return;
}

/*! \brief Initialize the super-droplet moisture model
 *
 * Allocates the moisture model variable MultiFabs and creates the
 * super-droplet particle container. This function sets up:
 * 1. The mapping between moisture variable indices and internal arrays
 * 2. MultiFabs for all moisture model variables
 * 3. The SuperDropletPC particle container
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
                        MicVar_SD::q_i,
                        MicVar_SD::q_g,
                        MicVar_SD::rh_w,
                        MicVar_SD::rain_accum,
                        MicVar_SD::graupel_accum,
                        MicVar_SD::snow_accum,
                        MicVar_SD::precip_accum };
    AMREX_ALWAYS_ASSERT(m_qmoist_size == m_mic_var_map.size());

    /* allocate microphysics multifabs */
    const int num_mic_vars =   MicVar_SD::NumVars
                             + m_num_nonmoist_sp * MicVar_SD_Species::NumVars
                             + m_num_aerosols * MicVar_SD_Aerosols::NumVars;
    m_mic_fab_vars.resize(1);
    m_mic_fab_vars[0].resize(num_mic_vars);
    for (auto i(0); i < num_mic_vars; i++) {
      m_mic_fab_vars[0][i] = std::make_shared<MultiFab> ( a_cons_vars.boxArray(),
                                                          a_cons_vars.DistributionMap(),
                                                          1,
                                                          a_cons_vars.nGrowVect() );
      m_mic_fab_vars[0][i]->setVal(0.0);
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
                   << "    include cold processes: " << (m_with_ice ? "true" : "false") << "\n"
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
void SuperDropletsMoist::InitParticles ( const int a_lev, MFPtr& a_z_phys_nd )
{
    BL_PROFILE("SuperDropletsMoist::InitParticles()");

    m_super_droplets->InitializeParticles(a_lev, zero, a_z_phys_nd);

    if (m_init_type != SDMoistInit::condensate_density) {
        amrex::Print() << "Initialized "
                       << m_super_droplets->NumberOfParticlesAtLevel(a_lev)
                       << " super-droplets at level " << a_lev
                       << " in super-droplets moisture model.\n";
    }
}

/*! \brief Restart particles in the super-droplet moisture model from checkpoint file
 *
 * This function restarts superdroplet particles from a checkpoint file.
 * It performs the following operations:
 * 1. Reads particle data from the specified restart file
 * 2. Redistributes particles to appropriate processors/grids
 * 3. Measures and reports the time taken to perform the restart
 * 4. Outputs statistics about the restarted particle population
 *
 * \param[in] a_gdb Unused particle grid database pointer
 * \param[in] a_fname File name for the checkpoint file to restart from
 */
void SuperDropletsMoist::RestartParticles ( ParGDBBase* /* a_gdb */, const std::string& a_fname )
{
    BL_PROFILE("SuperDropletsMoist::RestartParticles()");

    amrex::Print() << "Reading in " << m_name << " particle data from restart file.\n";

    m_super_droplets->Restart(a_fname, m_name);
    m_super_droplets->Redistribute();
    m_restarted_particles = true;

    amrex::Print() << "Restarted "
                   << m_super_droplets->NumSuperDroplets()
                   << " super-droplets representing "
                   << m_super_droplets->TotalNumberOfParticles()
                   << " particles in super-droplets moisture model.\n";
}

/*! \brief Complete initialization using now-available state variables
 *
 * This function finalizes initialization steps that depend on conserved state
 * variables that were not available during Init(). It performs:
 * 1. Particle density scaling based on air density
 * 2. For condensate_density initialization type: sets particle attributes from condensate density
 * 3. For other initialization types: optionally performs initial phase change relaxation
 * 4. Computes cloud/rain water and total water content for all species
 * 5. Updates the rhoq2 component in conserved variables with computed cloud water
 * 6. Runs initial diagnostics for superdroplets
 *
 * \param[in] a_lev AMR level
 * \param[in,out] a_cons_vars Conserved variables MultiFab to be updated
 * \param[in] a_z_phys_nd Vector of MultiFabs containing terrain information
 */
void SuperDropletsMoist::FinishInit (const int& a_lev,
                                     MultiFab& a_cons_vars,
                                     const Vector<MFPtr>& a_z_phys_nd)
{
    BL_PROFILE("SuperDropletsMoist::FinishInit()");
    const int lev = a_lev;
    if (lev == 0 && !m_restarted_particles) {
        // All native flow fields are available before this call. Finish local
        // particle attributes on every level before redistribution can move them.
        for (int k = 0; k <= m_super_droplets->finestLevel(); ++k) {
            m_super_droplets->DensityScaling(*m_mic_fab_vars[k][MicVar_SD::rho], k);
            if (m_init_type == SDMoistInit::condensate_density) {
                auto& qc = *m_mic_fab_vars[k][MicVar_SD::q_c];
                MultiFab rho_c(qc.boxArray(), qc.DistributionMap(), 1, qc.nGrowVect());
                MultiFab::Copy(rho_c, qc, 0, 0, 1, qc.nGrowVect());
                SetCurrentLevel(k);
                ratioToDensity(rho_c);
                m_super_droplets->SetAttributes(rho_c, k);
            }
        }
        SetCurrentLevel(lev);
        if (m_init_type == SDMoistInit::condensate_density) {
            amrex::Print() << "Initialized "
                           << m_super_droplets->NumSuperDroplets()
                           << " super-droplets representing "
                           << m_super_droplets->TotalNumberOfParticles()
                           << " particles in super-droplets moisture model.\n";
        }
    }
    if (lev == 0) {
        // Restart setup may also change the grids after particles were read.
        // Establish ownership before any projection or initial relaxation.
        m_super_droplets->Redistribute();
        if (!m_restarted_particles) {
            m_super_droplets->SplitMergeAtLevelBoundary();
        }
    }

    const bool relax = !m_restarted_particles &&
                       m_init_type != SDMoistInit::condensate_density &&
                       m_flag_phase_change && m_init_phase_change;
    if (relax) {
        phaseChange(m_init_phase_change_time, a_z_phys_nd, lev);
    }

    computeQcQrWater(*a_z_phys_nd[lev]);
    computeQiQgQsWater(*a_z_phys_nd[lev]);
    computeQtWater();

    // Restart state already includes fine averages in covered coarse cells.
    // A masked local particle projection must not replace that checkpoint data.
    if (!m_restarted_particles) {
        for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {
            const auto& box = mfi.tilebox();
            auto states_arr = a_cons_vars.array(mfi);
            auto q_c_arr = m_mic_fab_vars[lev][MicVar_SD::q_c]->array(mfi);
            auto q_r_arr = m_mic_fab_vars[lev][MicVar_SD::q_r]->array(mfi);
            auto q_i_arr = m_mic_fab_vars[lev][MicVar_SD::q_i]->array(mfi);
            auto q_g_arr = m_mic_fab_vars[lev][MicVar_SD::q_g]->array(mfi);
            auto q_s_arr = m_mic_fab_vars[lev][MicVar_SD::q_s]->array(mfi);
            ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                states_arr(i,j,k,RhoQ2_comp) = states_arr(i,j,k,Rho_comp)*q_c_arr(i,j,k);
                states_arr(i,j,k,RhoQ3_comp) = states_arr(i,j,k,Rho_comp)*q_i_arr(i,j,k);
                states_arr(i,j,k,RhoQ4_comp) = states_arr(i,j,k,Rho_comp)*q_r_arr(i,j,k);
                states_arr(i,j,k,RhoQ5_comp) = states_arr(i,j,k,Rho_comp)*q_s_arr(i,j,k);
                states_arr(i,j,k,RhoQ6_comp) = states_arr(i,j,k,Rho_comp)*q_g_arr(i,j,k);
            });
        }
    }

    // Appended condensate is a projection on every level, including before
    // the first step; fine-level state cannot wait for the first Advance.
    computeQcSpecies(*a_z_phys_nd[lev]);
    computeQtSpecies();

    if (!m_restarted_particles) {
        for (int is = m_istart_sp; is < m_num_species; is++) {
            for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {
                const auto& box = mfi.tilebox();
                auto states_arr = a_cons_vars.array(mfi);
                auto q_c_arr = m_mic_fab_vars[lev][s_qc_idx(is,m_istart_sp)]->array(mfi);
                auto qc_comp = q_qc_idx(is,m_istart_sp);
                ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    states_arr(i,j,k,qc_comp) = states_arr(i,j,k,Rho_comp)*q_c_arr(i,j,k);
                });
            }
        }
    }

    // Relaxation changes particle mass, vapor and latent heat together. Keep
    // that exchange in the initial conserved state, not just its condensate.
    if (relax && !m_kinematic_mode) {
        Copy_Micro_to_State(a_cons_vars);
    }

    if (lev == 0) {
        m_super_droplets->Diagnostics(-1, lev, zero, (m_diagnostics_iter>0));
    }

    return;
}

#endif

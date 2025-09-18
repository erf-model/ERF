#include <random>
#include "ERF_Constants.H"
#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! Add super-droplet method-specific attributes to particles */
void SuperDropletPC::add_superdroplet_attributes()
{
    BL_PROFILE("SuperDropletPC::add_superdroplets_attributes()");
    const bool communicate_this_comp = true;
    int count(0);
    for (int i = 0; i < SuperDropletsRealIdxSoA_RT::ncomps; i++) {
        AddRealComp(communicate_this_comp);
        count++;
    }
    for (int i = 0; i < m_num_aerosols; i++) {
        AddRealComp(communicate_this_comp);
        count++;
    }
    for (int i = 0; i < m_num_species; i++) {
        AddRealComp(communicate_this_comp);
        count++;
    }
    Print() << "SuperDropletPC(" << m_name << "): added " << count << " real-type attibute(s).\n";
    return;
}

/*! Read inputs from file */
void SuperDropletPC::readInputs ()
{
    BL_PROFILE("SuperDropletPC::readInputs");
    ParmParse pp(m_name);

    /* default values */
    m_density_scaling = false;
    m_nucleate_particles = false;
    m_advect_w_flow = true;
    m_advect_w_gravity = true;
    m_prescribed_advection = false;
    m_distribution_grid_size = 100;
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    m_bindist_rmin = 1e-6;
    m_bindist_rmax = 5e-3;
#endif
    m_sigma0 = 0.62;
    m_place_randomly_in_cells = true;
    m_recycle_threshold = 0.01;

    /* Newton solver parameters */
    m_newton_rtol = 1.0e-6;
    m_newton_atol = 1.0e-99;
    m_newton_stol = 1.0e-12;
    m_newton_maxits = 10;

    /* phase change eqn time integration */
    m_mass_change_cfl = 1000.0;
    m_mass_change_ti = SDMassChangeTIMethod::BE; // backward Euler

    /* log file for unconverged particles */
    m_mass_change_logging = false;
    m_mass_change_log_fname = "unconverged_superdroplets.log";

    /* recycled particle position bounds */
    const Geometry& geom = m_gdb->Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    m_recyc_xmin = plo[0];
    m_recyc_xmax = phi[0];
    m_recyc_ymin = plo[1];
    m_recyc_ymax = phi[1];
    m_recyc_zmin = plo[2];
    m_recyc_zmax = phi[2];

    std::string coal_kernel_name = "";
    std::string term_vel_name = "";

    m_coalescence_kernel = SDCoalescenceKernelType::sedimentation;
    coal_kernel_name = "sedimentation";
    m_include_brownian_coalescence = false;
    term_vel_name = "CloudRainShima";

    /* read these parameters if specified */
    pp.query("density_scaling", m_density_scaling);
    pp.query("nucleate_particles", m_nucleate_particles);
    pp.query("advect_with_flow", m_advect_w_flow);
    pp.query("advect_with_gravity", m_advect_w_gravity);
    pp.query("prescribed_advection", m_prescribed_advection);
    pp.query("newton_solver_rtol", m_newton_rtol);
    pp.query("newton_solver_atol", m_newton_atol);
    pp.query("newton_solver_stol", m_newton_stol);
    pp.query("newton_solver_maxits", m_newton_maxits);
    pp.query("mass_change_unconverged_log", m_mass_change_logging);
    pp.query("mass_change_unconverged_log_filename", m_mass_change_log_fname);
    pp.query("distribution_grid_size", m_distribution_grid_size);
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    pp.query("distribution_rmin", m_bindist_rmin);
    pp.query("distribution_rmax", m_bindist_rmax);
#endif
    pp.query("include_brownian_coalescence", m_include_brownian_coalescence);
    pp.query("sigma0", m_sigma0);
    pp.query("place_randomly_in_cells", m_place_randomly_in_cells);

    pp.query("recycle_threshold", m_recycle_threshold);
    pp.query("recycle_xmin", m_recyc_xmin);
    pp.query("recycle_xmax", m_recyc_xmax);
    pp.query("recycle_ymin", m_recyc_ymin);
    pp.query("recycle_ymax", m_recyc_ymax);
    pp.query("recycle_zmin", m_recyc_zmin);
    pp.query("recycle_zmax", m_recyc_zmax);

    std::string ti_name = "backward_euler";
    pp.query("mass_change_cfl", m_mass_change_cfl);
    pp.query("mass_change_ti_method", ti_name);
    if (ti_name == "rk3bs") {
        m_mass_change_ti = SDMassChangeTIMethod::RK3BS;
    } else if (ti_name == "rk4") {
        m_mass_change_ti = SDMassChangeTIMethod::RK4;
    } else if (ti_name == "backward_euler") {
        m_mass_change_ti = SDMassChangeTIMethod::BE;
    } else if (ti_name == "crank_nicolson") {
        m_mass_change_ti = SDMassChangeTIMethod::CN;
    } else if (ti_name == "dirk2") {
        m_mass_change_ti = SDMassChangeTIMethod::DIRK2;
    } else {
        amrex::Abort("Error in SuperDropletPC::readInputs() - invalid choice for mass change time integrator!");
    }

    pp.query("coalescence_kernel", coal_kernel_name);
    if (coal_kernel_name == "golovin") {
        m_coalescence_kernel = SDCoalescenceKernelType::golovin;
    } else if (coal_kernel_name == "sedimentation") {
        m_coalescence_kernel = SDCoalescenceKernelType::sedimentation;
    } else if (coal_kernel_name == "Longs") {
        m_coalescence_kernel = SDCoalescenceKernelType::Longs;
    } else if (coal_kernel_name == "Halls") {
        m_coalescence_kernel = SDCoalescenceKernelType::Halls;
    } else {
        amrex::Abort("Error in SuperDropletPC::readInputs() - invalid kernel choice!");
    }

    pp.query("terminal_velocity_model", term_vel_name);
    if (term_vel_name == "AtlasUlbrich") {
        m_term_vel_type = SDTerminalVelocityType::AtlasUlbrich;
    } else if (term_vel_name == "RogersYau") {
        m_term_vel_type = SDTerminalVelocityType::RogersYau;
    } else if (term_vel_name == "CloudRainShima") {
        m_term_vel_type = SDTerminalVelocityType::CloudRainShima;
    } else {
        amrex::Abort("Error in SuperDropletPC::readInputs() - invalid terminal velocity choice!");
    }

    {
        Vector<int> bin_size = {1,1,1};
        pp.queryarr("coalescence_bin_size", bin_size);
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            m_coalescence_bin_size[i] = bin_size[i];
        }
    }

    pp.query("num_initializations", m_num_initializations);
    m_initializations.resize(m_num_initializations);
    m_num_sd_per_cell = 0;
    const auto dx_h = Geom(m_lev).CellSize();
    const Real cell_volume = dx_h[0]*dx_h[1]*dx_h[2];
    for (int i = 0; i < m_num_initializations; i++) {
        m_initializations[i] = std::make_unique<SDInitialization>();
        m_initializations[i]->setDefaults(Geom(0), m_species_mat,m_aerosol_mat);

        char i_str[12]; sprintf(i_str, "%d", i);
        std::string prefix = m_name + "." + std::string(i_str);
        Print() << "Querying inputs file using the string: '" << m_name << "'!\n";
        m_initializations[i]->readInputs(m_name, Geom(0), m_species_mat, m_aerosol_mat);
        Print() << "Querying inputs file using the string: '" << prefix << "'!\n";
        m_initializations[i]->readInputs(prefix, Geom(0), m_species_mat, m_aerosol_mat);
        m_num_sd_per_cell += m_initializations[i]->numSDPerCell(cell_volume);
    }

    return;
}

/*! define super-droplets */
void SuperDropletPC::define (  const std::vector<Species::Name>& a_species_mat,
                               const std::vector<Species::Name>& a_aerosol_mat,
                               const BoxArray&                   a_ba,
                               const DistributionMapping&        a_dmap )
{
    BL_PROFILE("SuperDropletPC::define()");
    m_num_sd_per_cell = 0;
    m_num_unconverged_particles = 0;

    m_species_mat.clear();
    m_aerosol_mat.clear();

    setSpeciesMaterial( a_species_mat );
    m_num_species = m_species_mat.size();
    setAerosolMaterial( a_aerosol_mat );
    m_num_aerosols = m_aerosol_mat.size();

    AMREX_ASSERT(m_num_species  > 0);
    AMREX_ASSERT(m_num_species  <= SupDropInit::num_species_max);
    AMREX_ASSERT(m_num_aerosols <= SupDropInit::num_aerosols_max);

    add_superdroplet_attributes();
    readInputs();

#ifdef AMREX_USE_GPU
    AMREX_ASSERT(!m_mass_change_logging);
#endif
    if (m_mass_change_logging) {
        m_mass_change_log = fopen(m_mass_change_log_fname.c_str(), "w");
    }

    /* initialize random engine */
    {
        unsigned long int seed;
        int fix_seed = 0;
        ParmParse pp_erf("erf"); pp_erf.query("fix_random_seed", fix_seed);
        if (fix_seed) {
            Print() << "Using fixed seed for SuperDropletPC random engine.\n";
            seed = 1024UL;
        } else {
            std::random_device rd;
            std::uniform_int_distribution<unsigned long int> dist(0, std::numeric_limits<unsigned long int>::max());
            seed = dist(rd);
        }
        m_rndeng.seed(seed);
    }

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    m_mass_ln_R_mf.define(a_ba, a_dmap, m_distribution_grid_size, 0);
    m_num_ln_R_mf.define(a_ba, a_dmap, m_distribution_grid_size, 0);
#else
    amrex::ignore_unused(a_ba);
    amrex::ignore_unused(a_dmap);
#endif
}

/*! Initialize the particles */
void SuperDropletPC::InitializeParticles (const Real a_t, const MFPtr& a_ptr)
{
    amrex::ignore_unused(a_t);
    BL_PROFILE("SuperDropletPC::InitializeParticles()");
    Print() << "SuperDropletPC(" << m_name << "):\n"
            << "    Density scaling: " << (m_density_scaling ? "true" : "false") << "\n"
            << "    Nucleate particles: " << (m_nucleate_particles ? "true" : "false") << "\n"
            << "    Recycling threshold: " << m_recycle_threshold << "\n"
            << "    Recycling bounding box: " <<    "[" << m_recyc_xmin << ", " << m_recyc_xmax << "] "
                                              << " x [" << m_recyc_ymin << ", " << m_recyc_ymax << "] "
                                              << " x [" << m_recyc_zmin << ", " << m_recyc_zmax << "]\n"
            << "    Advect with flow: " << (m_advect_w_flow ? "true" : "false") << "\n"
            << "    Advect with gravity: " << (m_advect_w_gravity ? "true" : "false") << "\n"
            << "    Prescribed advection: " << (m_prescribed_advection ? "true" : "false") << "\n"
            << "    Random initial placement: " << (m_place_randomly_in_cells ? "true" : "false") << "\n"
            << "    Coalescence bin size: " << m_coalescence_bin_size << "\n"
            << "    Include Brownian coaslescence: "
            << (m_include_brownian_coalescence ? "true" : "false") << "\n";
    Print() << "    Coalescence kernel: ";
    if (m_coalescence_kernel == SDCoalescenceKernelType::golovin) {
        Print() << "golovin" << "\n";
    } else if (m_coalescence_kernel == SDCoalescenceKernelType::sedimentation) {
        Print() << "sedimentation" << "\n";
    } else if (m_coalescence_kernel == SDCoalescenceKernelType::Longs) {
        Print() << "Longs" << "\n";
    } else if (m_coalescence_kernel == SDCoalescenceKernelType::Halls) {
        Print() << "Halls" << "\n";
    }
    Print() << "    Mass change time integrator: ";
    if (m_mass_change_ti == SDMassChangeTIMethod::RK3BS) {
        Print() << "rk3bs";
    } else if (m_mass_change_ti == SDMassChangeTIMethod::RK4) {
        Print() << "rk4";
    } else if (m_mass_change_ti == SDMassChangeTIMethod::BE) {
        Print() << "backward_euler";
    } else if (m_mass_change_ti == SDMassChangeTIMethod::CN) {
        Print() << "crank_nicolson";
    } else if (m_mass_change_ti == SDMassChangeTIMethod::DIRK2) {
        Print() << "dirk2";
    }
    Print() << " (cfl = " << m_mass_change_cfl << ")\n";
    Print() << "    Terminal velocity model: ";
    if (m_term_vel_type == SDTerminalVelocityType::AtlasUlbrich) {
        Print() << "AtlasUlbrich" << "\n";
    } else if (m_term_vel_type == SDTerminalVelocityType::RogersYau) {
        Print() << "RogersYau" << "\n";
    } else if (m_term_vel_type ==  SDTerminalVelocityType::CloudRainShima) {
        Print() << "CloudRainShima" << "\n";
    }

    for (int i = 0; i < m_num_initializations; i++) {
        Print() << "SuperDropletPC(" << m_name << ") Initialization";
        if (m_num_initializations > 1) { Print() << " " << i; }
        Print() << ":\n";
        m_initializations[i]->printParameters(m_species_mat, m_aerosol_mat);
        initializeParticles( a_ptr, *(m_initializations[i]) );
        Print() << "    Particle container size: " << NumSuperDroplets() << "\n";
    }
}

/*! Sets the initial number of the super-droplets per cell as a box with a uniform distribution */
void SuperDropletPC::setNumSDBoxDistribution (iMultiFab& a_num_sd, /*!< integer Multifab with number of superdroplets in each grid cell */
                                              const int a_n_per_cell, /*!< number of superdroplets per cell */
                                              const MFPtr& a_height_ptr, /*!< terrain */
                                              const RealBox& a_box /*!< box within which to initialize particles */ )
{
    BL_PROFILE("SuperDropletPC::setNumSDBoxDistribution()");
    a_num_sd.setVal(0);

    const auto dx = Geom(m_lev).CellSizeArray();
    const auto plo = Geom(m_lev).ProbLoArray();

    for(MFIter mfi = MakeMFIter(m_lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();
        auto num_superdroplets_arr = a_num_sd[mfi].array();
        if (a_height_ptr) {
            const auto height_arr = (*a_height_ptr)[mfi].array();
            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x = plo[0] + (i + 0.5)*dx[0];
                Real y = plo[1] + (j + 0.5)*dx[1];
                Real z = 0.125 * (height_arr(i,j  ,k  ) + height_arr(i+1,j  ,k  ) +
                                  height_arr(i,j+1,k  ) + height_arr(i+1,j+1,k  ) +
                                  height_arr(i,j  ,k+1) + height_arr(i+1,j  ,k+1) +
                                  height_arr(i,j+1,k+1) + height_arr(i+1,j+1,k  ) );
                if (a_box.contains(RealVect(x,y,z))) {
                    num_superdroplets_arr(i,j,k) = a_n_per_cell;
                }
            });
        } else {
            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x = plo[0] + (i + 0.5)*dx[0];
                Real y = plo[1] + (j + 0.5)*dx[1];
                Real z = plo[2] + (k + 0.5)*dx[2];
                if (a_box.contains(RealVect(x,y,z))) {
                    num_superdroplets_arr(i,j,k) = a_n_per_cell;
                }
            });
        }
    }

    return;
}

/*! Sets the initial number of the super-droplets per cell as a bubble with a uniform distribution */
void SuperDropletPC::setNumSDBubbleDistribution ( iMultiFab& a_num_sd, /*!< integer Multifab with number of superdroplets in each grid cell */
                                                  const int a_n_per_cell, /*!< number of superdroplets per cell */
                                                  const MFPtr& a_height_ptr, /*!< terrain */
                                                  const RealBox& a_box /*!< box within which to initialize particles */ )
{
    BL_PROFILE("SuperDropletPC::setNumSDBubbleDistribution()");
    a_num_sd.setVal(0);

    const auto dx = Geom(m_lev).CellSizeArray();
    const auto plo = Geom(m_lev).ProbLoArray();

    for(MFIter mfi = MakeMFIter(m_lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();
        auto num_superdroplets_arr = a_num_sd[mfi].array();
        if (a_height_ptr) {
            const auto height_arr = (*a_height_ptr)[mfi].array();
            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x = plo[0] + (i + 0.5)*dx[0];
                Real y = plo[1] + (j + 0.5)*dx[1];
                Real z = 0.125 * (height_arr(i,j  ,k  ) + height_arr(i+1,j  ,k  ) +
                                  height_arr(i,j+1,k  ) + height_arr(i+1,j+1,k  ) +
                                  height_arr(i,j  ,k+1) + height_arr(i+1,j  ,k+1) +
                                  height_arr(i,j+1,k+1) + height_arr(i+1,j+1,k  ) );

                // Extract bubble params from box parameters
                const auto& x_c = a_box.lo();       // center
                const auto& x_r = a_box.hi();       // radius

                Real rad = 0.0;
                if (x_r[0] > 0) rad += std::pow((x - x_c[0])/x_r[0], 2);
                if (x_r[1] > 0) rad += std::pow((y - x_c[1])/x_r[1], 2);
                if (x_r[2] > 0) rad += std::pow((z - x_c[2])/x_r[2], 2);
                rad = std::sqrt(rad);

                if(rad <= 1.0){
                    num_superdroplets_arr(i,j,k) = a_n_per_cell;
                }
            });
        } else {
            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x = plo[0] + (i + 0.5)*dx[0];
                Real y = plo[1] + (j + 0.5)*dx[1];
                Real z = plo[2] + (k + 0.5)*dx[2];

                // Extract bubble params from box parameters
                const auto& x_c = a_box.lo();       // center
                const auto& x_r = a_box.hi();       // radius

                Real rad = 0.0;
                if (x_r[0] > 0) rad += std::pow((x - x_c[0])/x_r[0], 2);
                if (x_r[1] > 0) rad += std::pow((y - x_c[1])/x_r[1], 2);
                if (x_r[2] > 0) rad += std::pow((z - x_c[2])/x_r[2], 2);
                rad = std::sqrt(rad);

                if(rad <= 1.0){
                    num_superdroplets_arr(i,j,k) = a_n_per_cell;
                }
            });
        }
    }

    return;
}

/*! Initialize super-droplets in domain given an initialization type

    + The number of physical particles per cell is computed from the initial number density
      ("initial_number_density"), if specified (if not specified, it is taken to be 1).
    + The number of super-droplets per cell is computed from the initial super-droplet number
      density ("initial_super_droplet_density"), if specified (if not specified, it is set to
      the initial particles per cell ("inital_particles_per_cell"), whose default value is 1).
    + The initial particle mass is set to the sum of the condensate mass and aerosol mass, both
      of whose default values are 0.0.
    + The initial particle radius is the "equivalent radius" for the condensate material given
      the initial mass.
*/
void SuperDropletPC::initializeParticles ( const MFPtr& a_height_ptr, /*!< terrain */
                                           const SDInitialization& a_init /*!< initialization parameters */ )
{
    BL_PROFILE("SuperDropletPC::initializeParticles");

    const auto dx_h = Geom(m_lev).CellSize();
    const Real cell_volume = dx_h[0]*dx_h[1]*dx_h[2];
    const auto dx = Geom(m_lev).CellSizeArray();
    const auto plo = Geom(m_lev).ProbLoArray();

    // number of super-droplets per cell
    int num_sd_per_cell = a_init.numSDPerCell(cell_volume);
    m_num_sd_per_cell += num_sd_per_cell;

    // number of physical particles per cell
    Real num_par_per_cell = a_init.numParticlesPerCell(cell_volume);
    if (!num_par_per_cell) { return; }

    Print() << "    Number of physical particles per cell: " << num_par_per_cell << "\n"
            << "    Number of super droplets per cell: " << num_sd_per_cell << "\n";

    const int num_sp  = m_num_species;
    const int num_ae = m_num_aerosols;
    const Real rho_w = m_species_mat[m_idx_w]->m_density;
    const int idx_w = m_idx_w;
    AMREX_ALWAYS_ASSERT(idx_w >= 0);

    const auto sampled_multiplicity = a_init.sampledMultiplicity();

    iMultiFab num_superdroplets( ParticleBoxArray(m_lev),
                                 ParticleDistributionMap(m_lev),
                                 1, 0 );

    if (a_init.m_type == SupDropInit::init_uniform) {
        setNumSDBoxDistribution( num_superdroplets,
                                 num_sd_per_cell,
                                 a_height_ptr,
                                 a_init.m_init_particle_box );
    } else if (a_init.m_type == SupDropInit::init_bubble) {
        Print() << "Initializing as bubble distribution!";
        setNumSDBubbleDistribution( num_superdroplets,
                                    num_sd_per_cell,
                                    a_height_ptr,
                                    a_init.m_init_particle_box );
    } else if (a_init.m_type == SupDropInit::init_null) {
        num_superdroplets.setVal(0);
    } else {
        amrex::Print() << "Error: " << a_init.m_type
                        << " is not a valid initialization for "
                        << m_name << " particle species.\n";
        amrex::Error("See error message!");
    }

    iMultiFab offsets( ParticleBoxArray(m_lev),
                       ParticleDistributionMap(m_lev),
                       1, 0 );
    offsets.setVal(0);

    for(MFIter mfi = MakeMFIter(m_lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();

        int np = 0;
        {
            int ncell = num_superdroplets[mfi].numPts();
            const int* in = num_superdroplets[mfi].dataPtr();
            int* out = offsets[mfi].dataPtr();
            np = Scan::PrefixSum<int>( ncell,
                                       [=] AMREX_GPU_DEVICE (int i) -> int { return in[i]; },
                                       [=] AMREX_GPU_DEVICE (int i, int const &x) { out[i] = x; },
                                       Scan::Type::exclusive,
                                       Scan::retSum );
        }
        auto offset_arr = offsets[mfi].array();

        auto& particle_tile = DefineAndReturnParticleTile(m_lev, mfi);

        auto my_proc = ParallelDescriptor::MyProc();
        auto nprocs = ParallelDescriptor::NProcs();
        Long pid;
        {
            pid = ParticleType::NextID();
            ParticleType::NextID(pid+np);
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( static_cast<Long>(pid + np) < LastParticleID,
                                          "Error: overflow on particle id numbers!" );

        auto size_old = static_cast<Long>(particle_tile.size());
        auto size_new = size_old + np;
        particle_tile.resize(size_new);
        auto* aos = &particle_tile.GetArrayOfStructs()[0];
        auto& soa = particle_tile.GetStructOfArrays();

        /* SoA attributes */
        auto* vx_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::vx).data() + size_old;
        auto* vy_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::vy).data() + size_old;
        auto* vz_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::vz).data() + size_old;
        auto* mass_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::mass).data() + size_old;

        /* Runtime-added SoA attributes */
        int rt_offset = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::radius).data() + size_old;
        auto* mult_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::multiplicity).data() + size_old;
        auto* vterm_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::term_vel).data() + size_old;
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        auto* condt_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::cond_tendency).data() + size_old;
#endif
        auto* uid_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::uid).data() + size_old;

        SDSpeciesMassArr sp_mass_ptrs;
        for (int i = 0; i < num_sp; i++) {
            sp_mass_ptrs[i] = soa.GetRealData(idx_s(i,num_ae,num_sp)).data() + size_old;
        }

        SDAerosolMassArr ae_mass_ptrs;
        for (int i = 0; i < num_ae; i++) {
            ae_mass_ptrs[i] = soa.GetRealData(idx_a(i,num_ae,num_sp)).data() + size_old;
        }

        Gpu::DeviceVector<Real> species_mass_d(num_sp*np);
        {
            Vector<Real> multiplicity_h(np, 0.0);
            for (int i = 0; i < num_sp; i++) {
                Vector<Real> species_mass_h;
                if (sampled_multiplicity) {
                    a_init.getSpeciesDistribution( species_mass_h,
                                                   multiplicity_h,
                                                   cell_volume,
                                                   i,
                                                   np,
                                                   m_species_mat[i]->m_density,
                                                   m_rndeng );
                } else {
                    a_init.getSpeciesDistribution( species_mass_h,
                                                   i,
                                                   np,
                                                   m_species_mat[i]->m_density,
                                                   m_rndeng );
                }
                Gpu::copy( Gpu::hostToDevice,
                           species_mass_h.begin(),
                           species_mass_h.end(),
                           species_mass_d.begin() + (i*np) );
            }
        }
        Gpu::DeviceVector<Real> aerosol_mass_d(num_ae*np);
        Gpu::DeviceVector<Real> multiplicity_d(np);
        {
            Vector<Real> multiplicity_h(np, 0.0);
            for (int i = 0; i < num_ae; i++) {
                Vector<Real> aerosol_mass_h;
                if (sampled_multiplicity) {
                    a_init.getAerosolDistribution( aerosol_mass_h,
                                                   multiplicity_h,
                                                   cell_volume,
                                                   i,
                                                   np,
                                                   m_aerosol_mat[i]->m_density,
                                                   m_rndeng );
                } else {
                    a_init.getAerosolDistribution( aerosol_mass_h,
                                                   i,
                                                   np,
                                                   m_aerosol_mat[i]->m_density,
                                                   m_rndeng );
                }
                Gpu::copy( Gpu::hostToDevice,
                           aerosol_mass_h.begin(),
                           aerosol_mass_h.end(),
                           aerosol_mass_d.begin() + (i*np) );
            }
            if (sampled_multiplicity) {
                Gpu::copy( Gpu::hostToDevice,
                           multiplicity_h.begin(),
                           multiplicity_h.end(),
                           multiplicity_d.begin() );
            }
        }
        Gpu::synchronize();

        Gpu::DeviceVector<ParticleReal> sp_density(num_sp);
        Gpu::DeviceVector<int> sp_solubility(num_sp);
        {
            Vector<ParticleReal> sp_density_h(num_sp);
            Vector<int> sp_solubility_h(num_sp);
            for (int i = 0; i < num_sp; i++) {
                sp_density_h[i] = m_species_mat[i]->m_density;
                sp_solubility_h[i] = static_cast<int>(m_species_mat[i]->m_is_soluble);
            }
            Gpu::copy(  Gpu::hostToDevice,
                        sp_density_h.begin(),
                        sp_density_h.end(),
                        sp_density.begin() );
            Gpu::copy(  Gpu::hostToDevice,
                        sp_solubility_h.begin(),
                        sp_solubility_h.end(),
                        sp_solubility.begin() );
        }

        Gpu::DeviceVector<ParticleReal> ae_density(num_ae);
        Gpu::DeviceVector<int> ae_solubility(num_ae);
        {
            Vector<ParticleReal> ae_density_h(num_ae);
            Vector<int> ae_solubility_h(num_ae);
            for (int i = 0; i < num_ae; i++) {
                ae_density_h[i] = m_aerosol_mat[i]->m_density;
                ae_solubility_h[i] = static_cast<int>(m_aerosol_mat[i]->m_is_soluble);
            }
            Gpu::copy(  Gpu::hostToDevice,
                        ae_density_h.begin(),
                        ae_density_h.end(),
                        ae_density.begin() );
            Gpu::copy(  Gpu::hostToDevice,
                        ae_solubility_h.begin(),
                        ae_solubility_h.end(),
                        ae_solubility.begin() );
        }

        auto species_mass = species_mass_d.data();
        auto aerosol_mass = aerosol_mass_d.data();
        auto mult_arr = multiplicity_d.data();
        auto sp_rho_arr = sp_density.data();
        auto sp_sol_arr = sp_solubility.data();
        auto ae_rho_arr = ae_density.data();
        auto ae_sol_arr = ae_solubility.data();

        auto num_superdroplets_arr = num_superdroplets[mfi].array();
        auto random_place = m_place_randomly_in_cells;

        ParallelForRNG(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k, const RandomEngine& rnd_engine) noexcept
        {
            int num_sd_this_cell = num_superdroplets_arr(i,j,k);
            Real num_to_add = num_par_per_cell;
            Real n_par_per_supdrop = std::ceil(num_par_per_cell/num_sd_per_cell);

            Real mult_scale = 1.0;
            if (sampled_multiplicity) {
                Real mult_sum = 0.0;
                int start = offset_arr(i,j,k);
                for (int n = start; n < start+num_sd_this_cell; n++) { mult_sum += mult_arr[n]; }
                mult_scale = num_par_per_cell / mult_sum;
                //Print()<<"initializeParticles: Scale mult by mult_scale "<<mult_scale<<" based on mult_sum "<<mult_sum<< "\n";
            }

            int start = offset_arr(i,j,k);
            for (int n = start; n < start+num_sd_this_cell; n++) {

                auto& p = aos[n+size_old];
                p.id()  = pid + n;
                p.cpu() = my_proc;

                if (random_place) {
                    p.pos(0) = plo[0] + (i + Random(rnd_engine))*dx[0];
                    p.pos(1) = plo[1] + (j + Random(rnd_engine))*dx[1];
                    p.pos(2) = plo[2] + (k + Random(rnd_engine))*dx[2];
                } else {
                    p.pos(0) = plo[0] + (i + 0.5)*dx[0];
                    p.pos(1) = plo[1] + (j + 0.5)*dx[1];
                    p.pos(2) = plo[2] + (k + 0.5)*dx[2];
                }

                p.idata(SuperDropletsIntIdxAoS::k) = k;
                vx_ptr[n] = vy_ptr[n] = vz_ptr[n] = 0.0;

                Real mult_this_sd = 0;
                if (sampled_multiplicity) {
                    mult_this_sd = std::ceil(mult_arr[n] * mult_scale);
                } else {
                    mult_this_sd = n_par_per_supdrop;
                }
                if (mult_this_sd < num_to_add) {
                    mult_ptr[n] = mult_this_sd;
                } else {
                    mult_ptr[n] = num_to_add;
                    //Print() << " scaling SD# "<<n<< " : replace mult_ptr with num_to_add "<<num_to_add<<"\n";
                }
                num_to_add -= mult_ptr[n];
                if (mult_ptr[n] == 0) { mult_ptr[n] = 1; }

                for (int ctr = 0; ctr < num_sp; ctr++) {
                    sp_mass_ptrs[ctr][n] = species_mass[ctr*np+n];
                }
                for (int ctr = 0; ctr < num_ae; ctr++) {
                    ae_mass_ptrs[ctr][n] = aerosol_mass[ctr*np+n];
                }

                radius_ptr[n] = SD_effective_radius( n, idx_w,
                                                     rho_w,
                                                     num_sp, num_ae,
                                                     sp_sol_arr, ae_sol_arr,
                                                     sp_mass_ptrs, ae_mass_ptrs,
                                                     sp_rho_arr, ae_rho_arr );
                //Print()<<" SD# "<<n<<" with radius="<<radius_ptr[n]<< " m\n";
                mass_ptr[n] = SD_total_mass( n, num_sp, num_ae, sp_mass_ptrs, ae_mass_ptrs);

                vterm_ptr[n] = 0.0;
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
                condt_ptr[n] = 0.0;
#endif
                uid_ptr[n] = ParticleReal((pid+n-1)*nprocs + my_proc + 1);
            }
        });
        Gpu::synchronize();

        const auto height_arr = (*a_height_ptr)[mfi].array();
        ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int num_sd_this_cell = num_superdroplets_arr(i,j,k);
            int start = offset_arr(i,j,k);
            for (int n = start; n < start+num_sd_this_cell; n++) {
                auto& p = aos[n];
                Real x = p.pos(0);
                Real y = p.pos(1);
                Real z = p.pos(2);
                Real r[3] = { (x-plo[0])/dx[0] - i,
                              (y-plo[1])/dx[1] - j,
                              (z-plo[2])/dx[2] - k };

                Real sx[] = { amrex::Real(1.) - r[0], r[0]};
                Real sy[] = { amrex::Real(1.) - r[1], r[1]};

                Real height_at_pxy_lo = 0.;
                Real height_at_pxy_hi = 0.;
                for (int ii = 0; ii < 2; ++ii) {
                    for (int jj = 0; jj < 2; ++jj) {
                        height_at_pxy_lo += sx[ii] * sy[jj]
                                            * height_arr(i+ii,j+jj,k);
                        height_at_pxy_hi += sx[ii] * sy[jj]
                                            * height_arr(i+ii,j+jj,k+1);
                    }
                }

                p.pos(2) = height_at_pxy_lo
                           + r[2] * (height_at_pxy_hi - height_at_pxy_lo);
           }
        });
        Gpu::synchronize();
    }

    return;
}

/*! Set super-droplets multiplicity and radius in the domain from a given condensate mass density:
    Given the condensate mass density, we compute the mass of condensate per super-droplet from
    the mass per cell and the number of super-droplets per cell. We vary the multiplicity by a
    random amount for each super-droplet. The mass per physical particle is then computed from the
    mass per super-droplet and the multiplicity. The equivalent radius is then comptued from the
    particle mass and the density of condensate. */
void SuperDropletPC::SetAttributes (MultiFab& a_rhoc /*!< mass density of condensate */)
{
    BL_PROFILE("SuperDropletPC::SetAttributes");

    const auto plo = Geom(m_lev).ProbLoArray();
    const auto dx_h = Geom(m_lev).CellSize();
    const Real cell_volume = dx_h[0]*dx_h[1]*dx_h[2];
    const auto dxi = Geom(m_lev).InvCellSizeArray();
    const auto domain = Geom(m_lev).Domain();

    const int num_sd_per_cell = m_num_sd_per_cell;
    const int num_sp  = m_num_species;
    const int num_ae = m_num_aerosols;

    // condensate density
    const Real rho_w = m_species_mat[m_idx_w]->m_density;
    const int idx_w = m_idx_w;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, m_lev); pti.isValid(); ++pti) {
        auto& particle_tile = ParticlesAt(m_lev, pti);
        auto& soa = particle_tile.GetStructOfArrays();
        auto& aos = particle_tile.GetArrayOfStructs();
        auto *p_pbox = aos().data();
        const int n = aos.numParticles();

        /* SoA attributes */
        auto* mass_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::mass).data();

        /* Runtime-added SoA attributes */
        int rt_offset = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::radius).data();
        auto* mult_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::multiplicity).data();

        SDSpeciesMassArr sp_mass_ptrs;
        for (int i = 0; i < num_sp; i++) {
            sp_mass_ptrs[i] = soa.GetRealData(idx_s(i,num_ae,num_sp)).data();
        }

        SDAerosolMassArr ae_mass_ptrs;
        for (int i = 0; i < num_ae; i++) {
            ae_mass_ptrs[i] = soa.GetRealData(idx_a(i,num_ae,num_sp)).data();
        }

        auto condensate_mass_density = a_rhoc[pti.index()].array();

        ParallelForRNG(n, [=] AMREX_GPU_DEVICE (int i, const RandomEngine& rnd_engine)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            auto iv = getParticleCell(p, plo, dxi, domain);

            const Real mass_condensate_cell = condensate_mass_density(iv[0],iv[1],iv[2],0) * cell_volume;
            const Real mass_condensate_sd = mass_condensate_cell / num_sd_per_cell;

            Real mult_rnd = -mult_ptr[i]/3 + 2*mult_ptr[i]/3*Random(rnd_engine);
            mult_ptr[i] += mult_rnd;

            ParticleReal species_mass_total = 0.0;
            for (int ctr = 0; ctr < num_sp; ctr++) {
                if (ctr != idx_w) {
                    species_mass_total += sp_mass_ptrs[ctr][n];
                }
            }

            ParticleReal aerosol_mass_total = 0.0;
            for (int ctr = 0; ctr < num_ae; ctr++) {
                aerosol_mass_total += ae_mass_ptrs[ctr][n];
            }

            const Real mass_particle = mass_condensate_sd / mult_ptr[i] + aerosol_mass_total + species_mass_total;
            mass_ptr[i] = mass_particle;

            Real radius_cubed = mass_particle / ((4.0/3.0)*PI*rho_w);
            Real radius = (radius_cubed == 0.0 ? 0.0 : std::cbrt(radius_cubed));
            radius_ptr[i] = radius;
        });

    }

    return;
}

/*! Scale the multiplicities with density of air */
void SuperDropletPC::DensityScaling (const MultiFab& a_rho /*!< density of air */)
{
    BL_PROFILE("SuperDropletPC::DensityScaling");
    if (!m_density_scaling) { return; }

    const auto plo = Geom(m_lev).ProbLoArray();
    const auto dxi = Geom(m_lev).InvCellSizeArray();
    const auto domain = Geom(m_lev).Domain();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, m_lev); pti.isValid(); ++pti) {
        auto& particle_tile = ParticlesAt(m_lev, pti);
        auto& soa = particle_tile.GetStructOfArrays();
        auto& aos = particle_tile.GetArrayOfStructs();
        auto *p_pbox = aos().data();
        const int n = aos.numParticles();

        /* Runtime-added SoA attributes */
        int rt_offset = SuperDropletsRealIdxSoA::ncomps;
        auto* mult_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::multiplicity).data();

        auto density = a_rho[pti.index()].const_array();

        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            auto iv = getParticleCell(p, plo, dxi, domain);

            auto rho_air = density(iv[0],iv[1],iv[2],0);
            mult_ptr[i] *= rho_air;
        });

    }

    return;
}

#endif

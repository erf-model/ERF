#include <random>
#include "ERF_Constants.H"
#include "SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! Add super-droplet method-specific attributes to particles */
void SuperDropletPC::add_superdroplet_attributes()
{
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
    Print() << "SuperDropletPC(" << m_name << "): added " << count << " real-type attibute(s).\n";
    return;
}

/*! Read inputs from file */
void SuperDropletPC::readInputs ()
{
    ERFPC::readInputs();

    BL_PROFILE("SuperDropletPC::readInputs");
    ParmParse pp(m_name);

    /* default values */
    m_nucleate_particles = false;
    m_max_multiplicity = 1000000;
    m_numdens_init = -1;
    m_numdens_sd_init = m_numdens_init / m_max_multiplicity;
    m_advect_w_flow = true;
    m_advect_w_gravity = true;
    m_distribution_grid_size = 100;
    m_ppc_seed = 0;
    m_seed_mass = 0.0;

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

    std::string coal_kernel_name = "";
    std::string term_vel_name = "";
    if (m_vapour_mat->name() == MaterialNames::h2o) {
        m_mass_condensate_mean = 4.1887902e-42; // => radius = 1e-15 [m]
        m_radius_condensate_min = 1e-15;
        m_radius_condensate_max = 1e-15;
        m_condensate_init_type = SupDropInit::attrib_init_lnr;
        m_coalescence_kernel = SDCoalescenceKernelType::sedimentation;
        coal_kernel_name = "sedimentation";
        m_include_brownian_coalescence = false;
        term_vel_name = "CloudRainShima";
    }

    /* read these parameters if specified */
    pp.query("nucleate_particles", m_nucleate_particles);
    pp.query("initial_number_density", m_numdens_init);
    pp.query("initial_super_droplet_density", m_numdens_sd_init);
    pp.query("maximum_multiplicity", m_max_multiplicity);
    pp.query("initial_condensate_distribution_type", m_condensate_init_type);
    pp.query("initial_condensate_mass", m_mass_condensate_mean);
    pp.query("initial_condensate_min_radius", m_radius_condensate_min);
    pp.query("initial_condensate_max_radius", m_radius_condensate_max);
    pp.query("advect_with_flow", m_advect_w_flow);
    pp.query("advect_with_gravity", m_advect_w_gravity);
    pp.query("newton_solver_rtol", m_newton_rtol);
    pp.query("newton_solver_atol", m_newton_atol);
    pp.query("newton_solver_stol", m_newton_stol);
    pp.query("newton_solver_maxits", m_newton_maxits);
    pp.query("mass_change_unconverged_log", m_mass_change_logging);
    pp.query("mass_change_unconverged_log_filename", m_mass_change_log_fname);
    pp.query("distribution_grid_size", m_distribution_grid_size);
    pp.query("include_brownian_coalescence", m_include_brownian_coalescence);

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
    } else if (term_vel_name == "CloudRainShima") {
        m_term_vel_type = SDTerminalVelocityType::CloudRainShima;
    } else {
        amrex::Abort("Error in SuperDropletPC::readInputs() - invalid terminal velocity choice!");
    }

    pp.query("initial_seeds_per_cell", m_ppc_seed);
    pp.query("seed_condensate_mass", m_seed_mass);

    {
        Vector<int> bin_size = {1,1,1};
        pp.queryarr("coalescence_bin_size", bin_size);
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            m_coalescence_bin_size[i] = bin_size[i];
        }
    }

    {
        Vector<Real> particle_box_lo(AMREX_SPACEDIM);
        Vector<Real> particle_box_hi(AMREX_SPACEDIM);

        // Defaults
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            particle_box_lo[i] = Geom(0).ProbLo(i);
            particle_box_hi[i] = Geom(0).ProbHi(i);
        }

        pp.queryAdd("particle_box_lo", particle_box_lo, AMREX_SPACEDIM);
        AMREX_ASSERT(particle_box_lo.size() == AMREX_SPACEDIM);

        pp.queryAdd("particle_box_hi", particle_box_hi, AMREX_SPACEDIM);
        AMREX_ASSERT(particle_box_hi.size() == AMREX_SPACEDIM);

        m_init_particle_box.setLo(particle_box_lo);
        m_init_particle_box.setHi(particle_box_hi);
    }

    for (int i = 0; i < m_num_aerosols; i++) {
        // default values
        m_aerosol_init_type[i] = SupDropInit::attrib_init_const;
        if (m_aerosol_mat[i]->name() == MaterialNames::nacl) {
            m_mass_aerosol_min[i] = 1.0e-22;
            m_mass_aerosol_mean[i] = 1.0e-19;
            // the following values for radius will result
            // in mean salt mass of ~O(1e-19) kg
            m_radius_aerosol_min[i] = 1.0e-9;
            m_radius_aerosol_max[i] = 5.0e-8;
        }

        {
            std::string key = "initial_aerosol_distribution_type_"+m_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_aerosol_init_type[i]);
        }
        {
            std::string key = "initial_aerosol_min_mass_" + m_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_mass_aerosol_min[i]);
        }
        {
            std::string key = "initial_aerosol_mean_mass_" + m_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_mass_aerosol_mean[i]);
        }
        {
            std::string key = "initial_aerosol_min_radius_" + m_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_radius_aerosol_min[i]);
        }
        {
            std::string key = "initial_aerosol_max_radius_" + m_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_radius_aerosol_max[i]);
        }
    }

    return;
}

/*! Get real-type particle attribute names */
Vector<std::string> SuperDropletPC::varNames () const
{
    BL_PROFILE("ERFPCPC::varNames()");
    amrex::Vector<std::string> retval = {   AMREX_D_DECL("xvel","yvel","zvel"),
                                            "particle_mass",
                                            "radius",
                                            "multiplicity",
                                            "superdroplet_mass",
                                            "terminal_velocity",
                                            "t_coalescence" };
    for (int i = 0; i < m_num_aerosols; i++) {
        retval.push_back(std::string("aerosol_mass_"+m_aerosol_mat[i]->name()));
    }
    return retval;
}

/*! Get Eulerian plot variable names */
Vector<std::string> SuperDropletPC::meshPlotVarNames () const
{
    BL_PROFILE("ERFPCPC::varNames()");
    amrex::Vector<std::string> retval = { AMREX_D_DECL("mass_flux_x",
                                                       "mass_flux_y",
                                                       "mass_flux_z"),
                                          "number_density",
                                          "mass_density",
                                          "radius",
                                          ("mass_density_"+m_vapour_mat->name()),
                                          AMREX_D_DECL (
                                              ("mass_flux_x_"+m_vapour_mat->name()),
                                              ("mass_flux_y_"+m_vapour_mat->name()),
                                              ("mass_flux_z_"+m_vapour_mat->name()) ) };
    for (int i = 0; i < m_num_aerosols; i++) {
        retval.push_back(std::string("aerosol_mass_density_"+m_aerosol_mat[i]->name()));
        retval.push_back(std::string("aerosol_mass_flux_x_"+m_aerosol_mat[i]->name()));
        retval.push_back(std::string("aerosol_mass_flux_y_"+m_aerosol_mat[i]->name()));
        retval.push_back(std::string("aerosol_mass_flux_z_"+m_aerosol_mat[i]->name()));
    }
    return retval;
}

/*! define super-droplets */
void SuperDropletPC::define (  const std::string&              a_vap_mat,
                               const std::vector<std::string>& a_aerosol_mat)
{
    m_num_unconverged_particles = 0;

    m_aerosol_mat.clear();

    setVapourMaterial( a_vap_mat );
    setAerosolMaterial( a_aerosol_mat );
    m_num_aerosols = m_aerosol_mat.size();
    m_mass_aerosol_min.resize(m_num_aerosols);
    m_mass_aerosol_mean.resize(m_num_aerosols);
    m_radius_aerosol_min.resize(m_num_aerosols);
    m_radius_aerosol_max.resize(m_num_aerosols);
    m_aerosol_init_type.resize(m_num_aerosols);
    AMREX_ASSERT(m_num_aerosols <= SupDropInit::num_aerosols_max);

    add_superdroplet_attributes();
    readInputs();

#ifdef AMREX_USE_CUDA
    AMREX_ASSERT(!m_mass_change_logging);
#endif
    if (m_mass_change_logging) {
        m_mass_change_log = fopen(m_mass_change_log_fname.c_str(), "w");
    }
}

/*! Initialize super-droplets in domain given an initialization type */
void SuperDropletPC::InitializeParticles (const std::string& a_initialization_type, /*!< init type */
                                          const std::unique_ptr<amrex::MultiFab>& a_height_ptr /*!< terrain */)
{
    BL_PROFILE("SuperDropletPC::initializeParticles");

    if (a_initialization_type == SupDropInit::init_uniform) {
        initializeParticlesUniformDistribution( a_height_ptr, m_init_particle_box );
    } else if (a_initialization_type == SupDropInit::init_null) {
        initializeParticlesNull( a_height_ptr );
    } else {
        amrex::Print() << "Error: " << a_initialization_type
                        << " is not a valid initialization for "
                        << m_name << " particle species.\n";
        amrex::Error("See error message!");
    }
    return;
}

/*! Uniform distribution: initializes super-droplets uniformly throughout the domain.

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
void SuperDropletPC::initializeParticlesUniformDistribution (const std::unique_ptr<amrex::MultiFab>& a_height_ptr, /*!< terrain */
                                                             const RealBox& a_particle_init_domain /*!< box within which to initialize particles */ )
{
    BL_PROFILE("SuperDropletPC::initializeParticlesUniformDistribution");

    const auto dx = Geom(m_lev).CellSizeArray();
    const auto plo = Geom(m_lev).ProbLoArray();

    const Real cell_volume = dx[0]*dx[1]*dx[2];

    // number of super-droplets per cell
    int num_sd_per_cell = 0;
    if (m_numdens_sd_init >= 0) {
        num_sd_per_cell = static_cast<int>(std::ceil(m_numdens_sd_init*cell_volume));
    } else {
        num_sd_per_cell = m_ppc_init;
    }
    m_num_sd_per_cell = num_sd_per_cell;

    // number of physical particles per cell
    Real num_par_per_cell = 0.0;
    if (m_numdens_init >= 0) {
        num_par_per_cell = std::ceil(m_numdens_init*cell_volume);
        if (!num_par_per_cell) {
            return;
        }
    } else {
        num_par_per_cell = 1;
    }

    const int n_aerosols = m_num_aerosols;
    const int n_aerosols_max = SupDropInit::num_aerosols_max;
    const Real mat_density = m_vapour_mat->density();

    const int n_seeds = m_ppc_seed;
    const Real condensate_mass_seed = m_seed_mass;

    Print() << "SuperDropletPC(" << m_name << "):\n"
            << "    Number of physical particles per cell: " << num_par_per_cell << "\n"
            << "    Number of super droplets per cell: " << num_sd_per_cell << "\n"
            << "    Initial particle box: " << a_particle_init_domain << "\n"
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
    } else if (m_term_vel_type ==  SDTerminalVelocityType::CloudRainShima) {
        Print() << "CloudRainShima" << "\n";
    }

    Print() << "    Vapour material: " << m_vapour_mat->name() << "\n";

    Print() << "    Condensate initial distribution: " << m_condensate_init_type << " (";
    if (m_condensate_init_type == SupDropInit::attrib_init_const) {
        Print() << "value=" << m_mass_condensate_mean;
    } else if (m_condensate_init_type == SupDropInit::attrib_init_exp) {
        Print() << "mean=" << m_mass_condensate_mean;
    } else if (m_condensate_init_type == SupDropInit::attrib_init_lnr) {
        Print() << "min=" << m_radius_condensate_min
                << ", max=" << m_radius_condensate_max;
    }
    Print() << ")\n";

    if (m_aerosol_mat.size() > 0) {
        Print() << "    Aerosol materials:\n";
        for (unsigned long i=0; i < m_aerosol_mat.size(); i++) {
            Print() << "        "
                    << m_aerosol_mat[i]->name()
                    << " (Initial distribution: " << m_aerosol_init_type[i];
            if (m_aerosol_init_type[i] == SupDropInit::attrib_init_const) {
                Print() << ", value=" << m_mass_aerosol_mean[i];
            } else if (m_aerosol_init_type[i] == SupDropInit::attrib_init_exp) {
                Print() << ", min=" << m_mass_aerosol_min[i]
                        << ", mean=" << m_mass_aerosol_mean[i];
            } else if (m_aerosol_init_type[i] == SupDropInit::attrib_init_lnr) {
                Print() << ", min=" << m_radius_aerosol_min[i]
                        << ", max=" << m_radius_aerosol_max[i];
            }
            Print() << ")" << "\n";
        }
    }
    Print() << "    Number of seed particles per cell: " << m_ppc_seed << "\n"
            << "    Seed condensate mass: " << m_seed_mass << "\n";

    iMultiFab num_superdroplets( ParticleBoxArray(m_lev),
                                 ParticleDistributionMap(m_lev),
                                 1, 0 );
    num_superdroplets.setVal(0);
    for(MFIter mfi = MakeMFIter(m_lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();
        auto num_superdroplets_arr = num_superdroplets[mfi].array();
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
                if (a_particle_init_domain.contains(RealVect(x,y,z))) {
                    num_superdroplets_arr(i,j,k) = num_sd_per_cell;
                }
            });
        } else {
            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x = plo[0] + (i + 0.5)*dx[0];
                Real y = plo[1] + (j + 0.5)*dx[1];
                Real z = plo[2] + (k + 0.5)*dx[2];
                if (a_particle_init_domain.contains(RealVect(x,y,z))) {
                    num_superdroplets_arr(i,j,k) = num_sd_per_cell;
                }
            });
        }
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
        particle_tile.resize(np);
        auto* aos = &particle_tile.GetArrayOfStructs()[0];
        auto& soa = particle_tile.GetStructOfArrays();

        /* SoA attributes */
        auto* vx_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::vx).data();
        auto* vy_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::vy).data();
        auto* vz_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::vz).data();
        auto* mass_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::mass).data();

        /* Runtime-added SoA attributes */
        int rt_offset = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::radius).data();
        auto* supdrop_mass_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::sd_mass).data();
        auto* mult_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::multiplicity).data();
        auto* vterm_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::term_vel).data();
        auto* tcoal_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::t_coalescence).data();

        GpuArray<ParticleReal*,n_aerosols_max> aerosol_mass_ptrs;
        for (int i = 0; i < n_aerosols; i++) {
            aerosol_mass_ptrs[i] = soa.GetRealData(   rt_offset
                                                    + SuperDropletsRealIdxSoA_RT::ncomps
                                                    + i ).data();
        }

        auto my_proc = ParallelDescriptor::MyProc();
        Long pid;
        {
            pid = ParticleType::NextID();
            ParticleType::NextID(pid+np);
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( static_cast<Long>(pid + np) < LastParticleID,
                                          "Error: overflow on particle id numbers!" );

        Gpu::DeviceVector<Real> aerosol_mass_d;
        aerosol_mass_d.resize(n_aerosols*np);
        for (int i = 0; i < n_aerosols; i++) {
            auto aerosol_density = m_aerosol_mat[i]->density();
            Vector<Real> aerosol_mass_h;
            aerosol_mass_h.resize(np);
            if (m_aerosol_init_type[i] == SupDropInit::attrib_init_const) {
                for (int n = 0; n < np; n++) {
                    aerosol_mass_h[n] = m_mass_aerosol_mean[i];
                }
            } else if (m_aerosol_init_type[i] == SupDropInit::attrib_init_exp) {
                std::random_device rd;
                std::mt19937 rng(rd());
                auto delta = m_mass_aerosol_mean[i] - m_mass_aerosol_min[i];
                std::exponential_distribution<Real> ed(1.0/delta);
                for (int n = 0; n < np; n++) {
                    aerosol_mass_h[n] = ed(rng) + m_mass_aerosol_min[i];
                }
            } else if (m_aerosol_init_type[i] == SupDropInit::attrib_init_lnr) {
                std::random_device rd;
                std::mt19937 rng(rd());
                std::uniform_real_distribution<> urd(0.0, 1.0);
                Real delta =   std::log(m_radius_aerosol_max[i])
                             - std::log(m_radius_aerosol_min[i]);
                for (int n = 0; n < np; n++) {
                    Real term = std::log(m_radius_aerosol_min[i]) + urd(rng)*delta;
                    Real dry_r = std::exp(term);
                    aerosol_mass_h[n] = (4.0/3.0) * PI
                                        * dry_r * dry_r * dry_r
                                        * aerosol_density;
                }
            } else {
                Abort("Unknown m_aerosol_init_type!");
            }
            Gpu::copy( Gpu::hostToDevice,
                       aerosol_mass_h.begin(),
                       aerosol_mass_h.end(),
                       aerosol_mass_d.begin() + (i*np) );
        }

        Gpu::DeviceVector<Real> condensate_mass_d;
        {
            Vector<Real> condensate_mass_h;
            condensate_mass_h.resize(np);
            condensate_mass_d.resize(np);
            if (m_condensate_init_type == SupDropInit::attrib_init_exp) {
                std::random_device rd;
                std::mt19937 rng(rd());
                std::exponential_distribution<Real> ed(1.0/m_mass_condensate_mean);
                for (int n = 0; n < np; n++) {
                    condensate_mass_h[n] = ed(rng);
                }
            } else if (m_condensate_init_type == SupDropInit::attrib_init_const) {
                for (int n = 0; n < np; n++) {
                    condensate_mass_h[n] = m_mass_condensate_mean;
                }
            } else if (m_condensate_init_type == SupDropInit::attrib_init_lnr) {
                std::random_device rd;
                std::mt19937 rng(rd());
                std::uniform_real_distribution<> urd(0.0, 1.0);
                Real delta =    std::log(m_radius_condensate_max)
                             -  std::log(m_radius_condensate_min);
                for (int n = 0; n < np; n++) {
                    Real term = std::log(m_radius_condensate_min) + urd(rng)*delta;
                    Real radius = std::exp(term);
                    condensate_mass_h[n] =    (4.0/3.0) * PI
                                            * radius * radius * radius
                                            * mat_density;
                }
            } else {
                Abort("Unknown m_condensate_init_type!");
            }
            Gpu::copy( Gpu::hostToDevice,
                       condensate_mass_h.begin(),
                       condensate_mass_h.end(),
                       condensate_mass_d.begin() );
        }
        Gpu::synchronize();

        auto aerosol_mass = aerosol_mass_d.data();
        auto condensate_mass = condensate_mass_d.data();

        auto num_superdroplets_arr = num_superdroplets[mfi].array();

        ParallelForRNG(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                       const RandomEngine& rnd_engine) noexcept
        {
            int num_sd_this_cell = num_superdroplets_arr(i,j,k);
            Real num_to_add = num_par_per_cell;
            Real n_par_per_supdrop = std::ceil(num_par_per_cell/num_sd_per_cell);

            int start = offset_arr(i,j,k);
            for (int n = start; n < start+num_sd_this_cell; n++) {
                Real x = plo[0] + (i + Random(rnd_engine))*dx[0];
                Real y = plo[1] + (j + Random(rnd_engine))*dx[1];
                Real z = plo[2] + (k + Random(rnd_engine))*dx[2];

                Real multiplicity = 0;
                if (num_to_add > n_par_per_supdrop) {
                    multiplicity = n_par_per_supdrop;
                } else {
                    multiplicity = num_to_add;
                }
                num_to_add -= multiplicity;

                auto& p = aos[n];
                p.id()  = pid + n;
                p.cpu() = my_proc;
                p.pos(0) = x;
                p.pos(1) = y;
                p.pos(2) = z;

                p.idata(SuperDropletsIntIdxAoS::k) = k;

                vx_ptr[n] = vy_ptr[n] = vz_ptr[n] = 0.0;
                mult_ptr[n] = multiplicity;

                ParticleReal aerosol_mass_total = 0.0;
                for (int ctr = 0; ctr < n_aerosols; ctr++) {
                    aerosol_mass_ptrs[ctr][n] = aerosol_mass[ctr*np+n];
                    aerosol_mass_total += aerosol_mass[ctr*np+n];
                }

                ParticleReal cond_mass = 0.0, par_radius = 0.0;
                if (condensate_mass[n] > 0.0) {
                    cond_mass = condensate_mass[n];
                    par_radius = std::cbrt(cond_mass/((4.0/3.0)*PI*mat_density));
                } else {
                    par_radius = 1.0e-15;
                    cond_mass = (4.0/3.0)*PI*par_radius*par_radius*par_radius*mat_density;
                }

                mass_ptr[n] = cond_mass + aerosol_mass_total;
                radius_ptr[n] = par_radius;
                supdrop_mass_ptr[n] = mass_ptr[n]*multiplicity;
                vterm_ptr[n] = 0.0;
                tcoal_ptr[n] = 0.0;
            }

            /* Seed particles */
            for (int ns = 0; ns < n_seeds; ns++) {
                int n = start + Random_int(num_sd_this_cell, rnd_engine);

                ParticleReal aerosol_mass_total = 0.0;
                for (int ctr = 0; ctr < n_aerosols; ctr++) {
                    aerosol_mass_total += aerosol_mass_ptrs[ctr][n];
                }
                auto par_mass = condensate_mass_seed + aerosol_mass_total;
                auto par_radius = std::cbrt(par_mass/((4.0/3.0)*PI*mat_density));
                if (par_radius == 0.0) {
                    par_radius = 1.0e-16;
                }

                mult_ptr[n] = std::max( 1.0,
                                        std::floor(supdrop_mass_ptr[n]/par_mass) );

                mass_ptr[n] = par_mass;
                radius_ptr[n] = par_radius;
                supdrop_mass_ptr[n] = par_mass*mult_ptr[n];
            }
        });
        Gpu::synchronize();

        if (a_height_ptr) {

            const auto height_arr = (*a_height_ptr)[mfi].array();
            ParallelForRNG(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                           const RandomEngine& rnd_engine) noexcept
            {
                int num_sd_this_cell = num_superdroplets_arr(i,j,k);
                int start = offset_arr(i,j,k);
                for (int n = start; n < start+num_sd_this_cell; n++) {
                    auto& p = aos[n];
                    Real x = p.pos(0);
                    Real y = p.pos(1);
                    Real r[3] = { (x-plo[0])/dx[0] - i,
                                  (y-plo[1])/dx[1] - j,
                                  Random(rnd_engine) };

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

                    Real z = height_at_pxy_lo
                             + r[2] * (height_at_pxy_hi - height_at_pxy_lo);
                    p.pos(2) = z;
               }
            });
            Gpu::synchronize();

        }
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
    const auto dx = Geom(m_lev).CellSize();
    const auto dxi = Geom(m_lev).InvCellSizeArray();
    const auto domain = Geom(m_lev).Domain();

    const Real cell_volume = dx[0]*dx[1]*dx[2];
    const int num_sd_per_cell = m_num_sd_per_cell;

    const int n_aerosols = m_num_aerosols;
    const int n_aerosols_max = SupDropInit::num_aerosols_max;

    // condensate density
    const Real mat_density = m_vapour_mat->density();

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
        auto* supdrop_mass_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::sd_mass).data();

        GpuArray<ParticleReal*,n_aerosols_max> aerosol_mass_ptrs;
        for (int i = 0; i < n_aerosols; i++) {
            aerosol_mass_ptrs[i] = soa.GetRealData(   rt_offset
                                                    + SuperDropletsRealIdxSoA_RT::ncomps
                                                    + i ).data();
        }

        auto condensate_mass_density = a_rhoc[pti.index()].array();

        ParallelForRNG(n, [=] AMREX_GPU_DEVICE (int i, RandomEngine const& rnd_engine)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            auto iv = getParticleCell(p, plo, dxi, domain);

            const Real mass_condensate_cell = condensate_mass_density(iv[0],iv[1],iv[2],0) * cell_volume;
            const Real mass_condensate_sd = mass_condensate_cell / num_sd_per_cell;

            Real mult_rnd = -mult_ptr[i]/3 + 2*mult_ptr[i]/3*Random(rnd_engine);
            mult_ptr[i] += mult_rnd;

            ParticleReal aerosol_mass_total = 0.0;
            for (int ctr = 0; ctr < n_aerosols; ctr++) {
                aerosol_mass_total += aerosol_mass_ptrs[ctr][n];
            }
            const Real mass_particle = mass_condensate_sd / mult_ptr[i] + aerosol_mass_total;
            mass_ptr[i] = mass_particle;

            Real radius_cubed = mass_particle / ((4.0/3.0)*PI*mat_density);
            Real radius = (radius_cubed == 0.0 ? 0.0 : std::cbrt(radius_cubed));
            radius_ptr[i] = radius;
            supdrop_mass_ptr[i] = mass_ptr[i] * mult_ptr[i];
        });

    }

    return;
}

#endif

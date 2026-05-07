#include <random>
#include "ERF_Constants.H"
#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;
using namespace SDPCDefn;

/*! Add super-droplet method-specific attributes to particles */
void SuperDropletPC::add_superdroplet_attributes()
{
    BL_PROFILE("SuperDropletPC::add_superdroplets_attributes()");
    const bool communicate_this_comp = true;
    for (int i = 0; i < SuperDropletsIntIdxSoA_RT::ncomps; i++) {
        AddIntComp(communicate_this_comp);
    }
    Print() << "SuperDropletPC(" << m_name << "): added " << SuperDropletsIntIdxSoA_RT::ncomps << " int-type attribute(s).\n";
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
    Print() << "SuperDropletPC(" << m_name << "): added " << count << " real-type attribute(s).\n";
    return;
}

/*! Read inputs from file */
void SuperDropletPC::readInputs (const amrex::Real a_dt)
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
    m_sigma0 = Real(0.62);
    m_place_randomly_in_cells = true;
    m_deac_threshold = Real(0.01);
    m_save_inactive = false;

    /* Newton solver parameters */
    m_newton_rtol = Real(1.0e-6);
    m_newton_atol = Real(1.0e-99);
    m_newton_stol = Real(1.0e-12);
    m_newton_maxits = 10;

    /* phase change eqn time integration */
    m_mass_change_cfl = Real(1000.0);
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

    m_coalescence_kernel = SDCoalescenceKernelType::sedimentation;
    coal_kernel_name = "sedimentation";
    m_include_brownian_coalescence = false;
    m_term_vel_type_w = SDTerminalVelocityType::CloudRainShima;

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

    if (pp.contains("recycle_threshold")) {
        char err_msg[100];
        snprintf(err_msg, sizeof(err_msg), "Use \"inactive_threshold\" instead of \"recycle_threshold\"");
        amrex::Abort(err_msg);
    }
    pp.query("inactive_threshold", m_deac_threshold);
    pp.query("write_inactive_plt", m_save_inactive);
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

    {
        std::string inp_string = "terminal_velocity_model";
        if (pp.contains(inp_string.c_str())) {
            pp.get(inp_string.c_str(), m_term_vel_type_w);
        }
    }

    {
        Vector<int> bin_size = {1,1,1};
        pp.queryarr("coalescence_bin_size", bin_size);
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            m_coalescence_bin_size[i] = bin_size[i];
        }
    }

    const auto dx_h = Geom(m_lev).CellSize();
    const Real cell_volume = dx_h[0]*dx_h[1]*dx_h[2];

    pp.query("num_initializations", m_num_initializations);
    m_initializations.resize(m_num_initializations);
    m_num_sd_per_cell = 0;
    for (int i = 0; i < m_num_initializations; i++) {
        m_initializations[i] = std::make_unique<SDInitialization>();
        m_initializations[i]->setDefaults(Geom(0), m_species_mat,m_aerosol_mat);

        char i_str[12]; snprintf(i_str, sizeof(i_str), "%d", i);
        std::string prefix = m_name + "." + std::string(i_str);
        m_initializations[i]->readInputs(m_name, Geom(0), m_species_mat, m_aerosol_mat);
        m_initializations[i]->readInputs(prefix, Geom(0), m_species_mat, m_aerosol_mat);
        m_num_sd_per_cell += m_initializations[i]->numSDPerCell(cell_volume);
    }

    pp.query("num_injections", m_num_injections);
    m_injections.resize(m_num_injections);
    for (int i = 0; i < m_num_injections; i++) {
        m_injections[i] = std::make_unique<SDInjection>();
        m_injections[i]->setDefaults(Geom(0), m_species_mat,m_aerosol_mat);

        char i_str[12]; snprintf(i_str, sizeof(i_str), "%d", i);
        std::string str = m_name + ".injection";
        std::string prefix = str + "." + std::string(i_str);
        m_injections[i]->readInputs(str, Geom(0), m_species_mat, m_aerosol_mat, a_dt);
        m_injections[i]->readInputs(prefix, Geom(0), m_species_mat, m_aerosol_mat, a_dt);
    }

    return;
}

/*! define super-droplets */
void SuperDropletPC::define (  const std::vector<Species::Name>& a_species_mat,
                               const std::vector<Species::Name>& a_aerosol_mat,
                               const BoxArray&                   a_ba,
                               const DistributionMapping&        a_dmap,
                               const amrex::Real                 a_dt )
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

    AMREX_ALWAYS_ASSERT(m_num_species  > 0);
    AMREX_ALWAYS_ASSERT(m_num_species  <= SupDropInit::num_species_max);
    AMREX_ALWAYS_ASSERT(m_num_aerosols <= SupDropInit::num_aerosols_max);
    AMREX_ALWAYS_ASSERT(m_idx_w >= 0);

    add_superdroplet_attributes();
    readInputs(a_dt);

    // Initialize device properties for efficient GPU access
    initializeDeviceProperties();

    // Initialize staggered z-levels for non-uniform vertical grids (from terrain_z_levels)
    {
        ParmParse pp_erf("erf");
        int n_zlevels = pp_erf.countval("terrain_z_levels");
        if (n_zlevels > 0) {
            Vector<Real> zlevels_h(n_zlevels);
            pp_erf.getarr("terrain_z_levels", zlevels_h, 0, n_zlevels);
            m_zlevels_d.resize(n_zlevels);
            Gpu::copy(Gpu::hostToDevice, zlevels_h.begin(), zlevels_h.end(),
                      m_zlevels_d.begin());
        }
    }

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

    if (m_save_inactive) {
        m_mf_buf.define(a_ba, a_dmap, 1, 0);
        m_mf_buf.setVal(0.0);
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
            << "  Density scaling: " << (m_density_scaling ? "true" : "false") << "\n"
            << "  Nucleate particles: " << (m_nucleate_particles ? "true" : "false") << "\n"
            << "  Inactive particles threshold: " << m_deac_threshold << "\n"
            << "  Recycling bounding box: " <<    "[" << m_recyc_xmin << ", " << m_recyc_xmax << "] "
                                              << " x [" << m_recyc_ymin << ", " << m_recyc_ymax << "] "
                                              << " x [" << m_recyc_zmin << ", " << m_recyc_zmax << "]\n"
            << "  Advect with flow: " << (m_advect_w_flow ? "true" : "false") << "\n"
            << "  Advect with gravity: " << (m_advect_w_gravity ? "true" : "false") << "\n"
            << "  Prescribed advection: " << (m_prescribed_advection ? "true" : "false") << "\n"
            << "  Random initial placement: " << (m_place_randomly_in_cells ? "true" : "false") << "\n"
            << "  Coalescence bin size: " << m_coalescence_bin_size << "\n"
            << "  Include Brownian coaslescence: "
            << (m_include_brownian_coalescence ? "true" : "false") << "\n";
    Print() << "  Coalescence kernel: ";
    if (m_coalescence_kernel == SDCoalescenceKernelType::golovin) {
        Print() << "golovin" << "\n";
    } else if (m_coalescence_kernel == SDCoalescenceKernelType::sedimentation) {
        Print() << "sedimentation" << "\n";
    } else if (m_coalescence_kernel == SDCoalescenceKernelType::Longs) {
        Print() << "Longs" << "\n";
    } else if (m_coalescence_kernel == SDCoalescenceKernelType::Halls) {
        Print() << "Halls" << "\n";
    }
    Print() << "  Mass change time integrator: ";
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
    Print() << "    Terminal velocity model: "
            << getEnumNameString(m_term_vel_type_w) << "\n";

    for (int i = 0; i < m_num_initializations; i++) {
        Print() << "  SuperDropletPC(" << m_name << ") Initialization";
        if (m_num_initializations > 1) { Print() << " " << i; }
        Print() << ":\n";
        m_initializations[i]->printParameters(m_species_mat, m_aerosol_mat);
        addParticles( a_ptr, *(m_initializations[i]) );
        Print() << "  Particle container size: " << NumSuperDroplets() << "\n";
    }
    Print() << "  Total number of superdroplets per cell: " << m_num_sd_per_cell << "\n";

    for (int i = 0; i < m_num_injections; i++) {
        Print() << "  SuperDropletPC(" << m_name << ") Injection";
        if (m_num_injections > 1) { Print() << " " << i; }
        Print() << ":\n";
        m_injections[i]->printParameters(m_species_mat, m_aerosol_mat);
    }
}

/*! Inject particles */
void SuperDropletPC::InjectParticles (const Real a_t, const MFPtr& a_ptr, const Real a_dt)
{
    amrex::ignore_unused(a_t);

    for (int i = 0; i < m_num_injections; i++) {
        m_injections[i]->updateDt(a_dt);
        if (    (m_injections[i]->m_inj_rate > 0)
             && (a_t >= m_injections[i]->m_tstart)
             && (a_t <= m_injections[i]->m_tstop) ) {
            Print() << "SuperDropletPC(" << m_name << "): "
                    << " injecting particles (" << i << ").\n";
            addParticles( a_ptr, *(m_injections[i]) );
            Print() << "    Particle container size: " << NumSuperDroplets() << "\n";
        }
    }
}

/*! Set super-droplets multiplicity and radius in the domain from a given condensate mass density:
    Given the condensate mass density, we compute the mass of condensate per super-droplet from
    the mass per cell and the number of super-droplets per cell. We vary the multiplicity by a
    random amount for each super-droplet. The mass per physical particle is then computed from the
    mass per super-droplet and the multiplicity. The equivalent radius is then computed from the
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

    forEachParticleTile(m_lev,
        [&](ParIterType& /*pti*/, int grid, ParticleType* p_pbox,
            const SDProcess::ParticlePointers& ptrs, int np)
    {
        auto condensate_mass_density = a_rhoc[grid].array();

        ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int i, const RandomEngine& rnd_engine)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            auto iv = getParticleCell(p, plo, dxi, domain);

            const Real mass_condensate_cell = condensate_mass_density(iv[0],iv[1],iv[2],0) * cell_volume;
            const Real mass_condensate_sd = mass_condensate_cell / num_sd_per_cell;

            Real mult_rnd = -ptrs.mult_ptr[i]/3 + 2*ptrs.mult_ptr[i]/3*Random(rnd_engine);
            ptrs.mult_ptr[i] += mult_rnd;

            ParticleReal species_mass_total = zero;
            for (int ctr = 0; ctr < num_sp; ctr++) {
                if (ctr != idx_w) {
                    species_mass_total += ptrs.sp_mass_ptrs[ctr][np];
                }
            }

            ParticleReal aerosol_mass_total = zero;
            for (int ctr = 0; ctr < num_ae; ctr++) {
                aerosol_mass_total += ptrs.ae_mass_ptrs[ctr][np];
            }

            const Real mass_particle = mass_condensate_sd / ptrs.mult_ptr[i] + aerosol_mass_total + species_mass_total;
            ptrs.mass_ptr[i] = mass_particle;

            Real radius_cubed = mass_particle / (four_thirds_pi*rho_w);
            Real radius = (radius_cubed == zero ? zero : std::cbrt(radius_cubed));
            ptrs.radius_ptr[i] = radius;
        });
    }); // end forEachParticleTile
}

/*! Scale the multiplicities with density of air */
void SuperDropletPC::DensityScaling (const MultiFab& a_rho /*!< density of air */)
{
    BL_PROFILE("SuperDropletPC::DensityScaling");
    if (!m_density_scaling) { return; }

    const auto plo = Geom(m_lev).ProbLoArray();
    const auto dxi = Geom(m_lev).InvCellSizeArray();
    const auto domain = Geom(m_lev).Domain();

    forEachParticleTile(m_lev,
        [&](ParIterType& /*pti*/, int grid, ParticleType* p_pbox,
            const SDProcess::ParticlePointers& ptrs, int np)
    {
        auto density = a_rho[grid].const_array();

        ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            auto iv = getParticleCell(p, plo, dxi, domain);

            auto rho_air = density(iv[0],iv[1],iv[2],0);
            ptrs.mult_ptr[i] *= rho_air;
        });
    }); // end forEachParticleTile
}

#endif

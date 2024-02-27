#include <ERF_Constants.H>
#include <SuperDropletPC.H>

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
    Print() << "SuperDropletPC: added " << count << " real-type attibute(s).\n";
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
    m_mass_condensate_init = 0.0;

    /* read these parameters if specified */
    pp.query("nucleate_particles", m_nucleate_particles);
    pp.query("initial_number_density", m_numdens_init);
    pp.query("initial_super_droplet_density", m_numdens_sd_init);
    pp.query("maximum_multiplicity", m_max_multiplicity);
    pp.query("initial_condensate_mass", m_mass_condensate_init);

    for (int i = 0; i < m_num_aerosols; i++) {
        m_mass_aerosol_init[i] = 0.0;
        std::string key = "initial_aerosol_mass_" + m_aerosol_mat[i]->name();
        pp.query(key.c_str(), m_mass_aerosol_init[i]);
    }

    return;
}

/*! Initialize super-droplets in domain given an initialization type */
void SuperDropletPC::InitializeParticles (const std::string& a_initialization_type, /*!< init type */
                                          const std::unique_ptr<amrex::MultiFab>& a_height_ptr /*!< terrain */)
{
    BL_PROFILE("SuperDropletPC::initializeParticles");

    if (a_initialization_type == SuperDropletInitializations::init_uniform) {
        initializeParticlesUniformDistribution( a_height_ptr );
    } else if (a_initialization_type == SuperDropletInitializations::init_null) {
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
void SuperDropletPC::initializeParticlesUniformDistribution (const std::unique_ptr<amrex::MultiFab>& a_height_ptr /*!< terrain */)
{
    BL_PROFILE("SuperDropletPC::initializeParticlesUniformDistribution");

    const Real* dx = Geom(m_lev).CellSize();
    const Real* plo = Geom(m_lev).ProbLo();

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

    // initial aerosol mass of each physical particle
    const int num_aerosols = m_num_aerosols;
    Array<Real,SuperDropletInitializations::num_aerosols_max> aerosol_mass;
    Real aerosol_mass_total = 0.0;
    for (int i = 0; i < num_aerosols; i++) {
        aerosol_mass[i] = m_mass_aerosol_init[i];
        aerosol_mass_total += m_mass_aerosol_init[i];
    }
    // initial mass of each physical particle
    const Real par_mass = m_mass_condensate_init + aerosol_mass_total;

    // initial radius
    const Real mat_density = m_vapour_mat->density();
    const Real par_radius = std::exp( std::log(par_mass / ((4.0/3.0)*PI*mat_density))/3.0 );

    Print() << "SuperDropletPC(" << m_name << "):\n"
            << "    Number of physical particles per cell: " << num_par_per_cell << "\n"
            << "    Number of super droplets per cell: " << num_sd_per_cell << "\n"
            << "    Initial radius: " << par_radius << "\n"
            << "    Initial mass: " << par_mass << "\n"
            << "    Vapour material: " << m_vapour_mat->name() << "\n";
    if (m_aerosol_mat.size() > 0) {
        Print() << "    Aerosol materials:\n";
        for (int i=0; i < m_aerosol_mat.size(); i++) {
            Print() << "        " << m_aerosol_mat[i]->name() << "\n";
        }
    }



    iMultiFab num_superdroplets( ParticleBoxArray(m_lev),
                                 ParticleDistributionMap(m_lev),
                                 1, 0 );
    num_superdroplets.setVal(0);
    for(MFIter mfi = MakeMFIter(m_lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();
        auto num_superdroplets_arr = num_superdroplets[mfi].array();
        ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            num_superdroplets_arr(i,j,k) = num_sd_per_cell;
        });
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

        GpuArray<ParticleReal*,SuperDropletInitializations::num_aerosols_max> aerosol_mass_ptrs;
        for (int i = 0; i < num_aerosols; i++) {
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

        ParallelForRNG(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                       const RandomEngine& rnd_engine) noexcept
        {
            Real num_to_add = num_par_per_cell;
            Real n_par_per_supdrop = std::ceil(num_par_per_cell/num_sd_per_cell);

            int start = offset_arr(i,j,k);
            for (int n = start; n < start+num_sd_per_cell; n++) {
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

                vx_ptr[n] = 0.0;
                vy_ptr[n] = 0.0;
                vz_ptr[n] = 0.0;

                mass_ptr[n] = par_mass;

                radius_ptr[n] = par_radius;
                mult_ptr[n] = multiplicity;
                supdrop_mass_ptr[n] = par_mass*multiplicity;

                for (int i = 0; i < num_aerosols; i++) {
                    aerosol_mass_ptrs[i][n] = aerosol_mass[i];
                }
           }
        });

        if (a_height_ptr) {

            const auto height_arr = (*a_height_ptr)[mfi].array();
            ParallelForRNG(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                           const RandomEngine& rnd_engine) noexcept
            {
                int start = offset_arr(i,j,k);
                for (int n = start; n < start+num_sd_per_cell; n++) {
                    Real r = Random(rnd_engine);
                    Real z = height_arr(i,j,k) + r * ( height_arr(i,j,k+1) - height_arr(i,j,k) );
                    auto& p = aos[n];
                    p.pos(2) = z;
               }
            });

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
void SuperDropletPC::SetAttributes (MultiFab& a_mass_density /*!< mass density of condensate */)
{
    BL_PROFILE("SuperDropletPC::SetAttributes");

    const auto plo = Geom(m_lev).ProbLoArray();
    const auto dx = Geom(m_lev).CellSize();
    const auto dxi = Geom(m_lev).InvCellSizeArray();
    const auto domain = Geom(m_lev).Domain();

    const Real cell_volume = dx[0]*dx[1]*dx[2];
    const int num_sd_per_cell = m_num_sd_per_cell;

    Real mass_aerosol = 0.0;
    for (int i = 0; i < m_num_aerosols; i++) {
        mass_aerosol += m_mass_aerosol_init[i];
    }

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

        auto mass_density = a_mass_density[pti.index()].array();

        ParallelForRNG(n, [=] AMREX_GPU_DEVICE (int i, RandomEngine const& rnd_engine)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            auto iv = getParticleCell(p, plo, dxi, domain);

            const Real mass_cell = mass_density(iv[0],iv[1],iv[2],0) * cell_volume;
            const Real mass_sd = mass_cell / num_sd_per_cell;

            Real mult_rnd = -mult_ptr[i]/3 + 2*mult_ptr[i]/3*Random(rnd_engine);
            mult_ptr[i] += mult_rnd;

            const Real mass_particle = mass_sd / mult_ptr[i] + mass_aerosol;
            mass_ptr[i] = mass_particle;

            Real radius_cubed = mass_particle / ((4.0/3.0)*PI*mat_density);
            Real radius = std::exp( std::log(radius_cubed) / 3.0 );
            radius_ptr[i] = radius;
            supdrop_mass_ptr[i] = mass_ptr[i] * mult_ptr[i];
        });

    }

    return;
}

#endif

#include <ERF_Constants.H>
#include <SuperDropletPC.H>

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! Add super-droplet method-specific attributes to particles */
void SuperDropletPC::add_superdroplet_attributes()
{
    const bool communicate_this_comp = true;
    int count(0);
    for (int i = 0; i < SuperDropletsRealIdxSoA::ncomps; i++) {
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
    m_mass_aerosol_init = 0.0;

    /* read these parameters if specified */
    pp.query("nucleate_particles", m_nucleate_particles);
    pp.query("initial_number_density", m_numdens_init);
    pp.query("initial_super_droplet_density", m_numdens_sd_init);
    pp.query("maximum_multiplicity", m_max_multiplicity);
    pp.query("initial_condensate_mass", m_mass_condensate_init);
    pp.query("initial_aerosol_mass", m_mass_aerosol_init);

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
    if (m_numdens_sd_init >= 0) {
        m_num_sd_per_cell = static_cast<int>(std::ceil(m_numdens_sd_init*cell_volume));
    } else {
        m_num_sd_per_cell = m_ppc_init;
    }
    // number of physical particles per cell
    if (m_numdens_init >= 0) {
        m_num_par_per_cell = std::ceil(m_numdens_init*cell_volume);
        if (!m_num_par_per_cell) {
            return;
        }
    } else {
        m_num_par_per_cell = 1;
    }

    // initial mass of each physical particle
    const Real par_mass = m_mass_condensate_init + m_mass_aerosol_init;

    // initial radius
    const Real mat_density = 1.0; /* TODO: material object */
    const Real par_radius = std::exp( std::log(par_mass / ((4.0/3.0)*PI*mat_density))/3.0 );

    Print() << "SuperDropletPC(" << m_name << "):\n"
            << "    Number of physical particles per cell: " << m_num_par_per_cell << "\n"
            << "    Number of super droplets per cell: " << m_num_sd_per_cell << "\n"
            << "    Initial radius: " << par_radius << "\n"
            << "    Initial mass: " << par_mass << "\n";

    for(MFIter mfi = MakeMFIter(m_lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();

        const FArrayBox* height_ptr = nullptr;
        if (a_height_ptr) {
            const auto& height = (*a_height_ptr)[mfi];
#ifdef AMREX_USE_GPU
            std::unique_ptr<FArrayBox> hostfab;
            if (height.arena()->isManaged() || height.arena()->isDevice()) {
                hostfab = std::make_unique<FArrayBox>(height.box(), height.nComp(),
                                                      The_Pinned_Arena());
                Gpu::dtoh_memcpy_async(hostfab->dataPtr(), height.dataPtr(),
                                       height.size()*sizeof(Real));
                Gpu::streamSynchronize();
                height_ptr = hostfab.get();
            }
#else
            height_ptr = &height;
#endif
        }

        auto& particle_tile = DefineAndReturnParticleTile(m_lev, mfi);
        RandomEngine rnd_engine;

        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {

            Real num_to_add = m_num_par_per_cell;
            Real n_par_per_supdrop = std::ceil(m_num_par_per_cell/m_num_sd_per_cell);

            for (int n = 0; n < m_num_sd_per_cell; n++) {
                Real r[3] = {Random(rnd_engine), Random(rnd_engine), Random(rnd_engine)};
                Real v[3] = {0.0, 0.0, 0.0};

                Real x = plo[0] + (iv[0] + r[0])*dx[0];
                Real y = plo[1] + (iv[1] + r[1])*dx[1];
                Real z = plo[2] + (iv[2] + r[2])*dx[2];
                if (height_ptr) {
                    z = (*height_ptr)(iv)
                         + r[2] * (   (*height_ptr)(iv + IntVect(AMREX_D_DECL(0, 0, 1)))
                                    - (*height_ptr)(iv) );
                }

                if (!num_to_add) {
                    break;
                }

                Real multiplicity = 0;
                if (num_to_add > n_par_per_supdrop) {
                    multiplicity = n_par_per_supdrop;
                } else {
                    multiplicity = num_to_add;
                }
                num_to_add -= multiplicity;

                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                p.pos(0) = x;
                p.pos(1) = y;
                p.pos(2) = z;

                p.rdata(SuperDropletsRealIdxAoS::vx) = v[0];
                p.rdata(SuperDropletsRealIdxAoS::vy) = v[1];
                p.rdata(SuperDropletsRealIdxAoS::vz) = v[2];

                p.rdata(SuperDropletsRealIdxAoS::mass) = par_mass;

                p.idata(SuperDropletsIntIdxAoS::i) = iv[0];
                p.idata(SuperDropletsIntIdxAoS::j) = iv[1];
                p.idata(SuperDropletsIntIdxAoS::k) = iv[2];

                particle_tile.push_back(p);

                particle_tile.push_back_real(SuperDropletsRealIdxSoA::radius, par_radius);
                particle_tile.push_back_real(SuperDropletsRealIdxSoA::sol_mass, m_mass_aerosol_init);
                particle_tile.push_back_real(SuperDropletsRealIdxSoA::multiplicity, multiplicity);
           }
           AMREX_ALWAYS_ASSERT(num_to_add < 1);
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

    // condensate density
    const Real mat_density = 1.0; /* TODO: material object */

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, m_lev); pti.isValid(); ++pti) {
        auto& particle_tile = ParticlesAt(m_lev, pti);
        auto& soa = particle_tile.GetStructOfArrays();
        auto& aos = particle_tile.GetArrayOfStructs();
        auto *p_pbox = aos().data();
        const int n = aos.numParticles();

        auto* mult_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::multiplicity).data();
        auto* radius_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::radius).data();

        auto mass_density = a_mass_density[pti.index()].array();

        Real mass_aerosol = m_mass_aerosol_init;
        int num_sd_per_cell = m_num_sd_per_cell;

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
            p.rdata(SuperDropletsRealIdxAoS::mass) = mass_particle;

            Real radius_cubed = mass_particle / ((4.0/3.0)*PI*mat_density);
            Real radius = std::exp( std::log(radius_cubed) / 3.0 );
            radius_ptr[i] = radius;

        });

    }

    return;
}

#endif

#include <ERF_Constants.H>
#include <SuperDropletPC.H>

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! Add super-droplet method-specific attributes to particles */
void SuperDropletPC::add_superdroplet_attributes()
{
    const bool communicate_this_comp = true;

    {
        int count(0);
        for (int i = 0; i < SuperDropletsRealIdxSoA::ncomps; i++) {
            AddRealComp(communicate_this_comp);
            count++;
        }
        Print() << "SuperDropletPC: added " << count << " real-type attibute(s).\n";
    }
    {
        int count(0);
        for (int i = 0; i < SuperDropletsIntIdxSoA::ncomps; i++) {
            AddIntComp(communicate_this_comp);
            count++;
        }
        Print() << "SuperDropletPC: added " << count << " int-type attibute(s).\n";
    }

    return;
}

/*! Read inputs from file */
void SuperDropletPC::readInputs ()
{
    ERFPC::readInputs();

    BL_PROFILE("SuperDropletPC::readInputs");
    ParmParse pp(m_name);

    m_nucleate_particles = false;
    pp.query("nucleate_particles", m_nucleate_particles);

    m_numdens_init = -1;
    pp.query("initial_number_density", m_numdens_init);

    m_radius_init = -1;
    pp.query("initial_radius", m_radius_init);

    m_max_multiplicity = 100;
    pp.query("maximum_multiplicity", m_max_multiplicity);

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

/*! Uniform distribution: the number of super-droplets per grid cell is specified
    by "initial_particles_per_cell" (default is 1), and they are randomly distributed.
    The multiplicity is determined from "initial_number_density" (if specified) or 1.*/
void SuperDropletPC::initializeParticlesUniformDistribution (const std::unique_ptr<amrex::MultiFab>& a_height_ptr /*!< terrain */)
{
    BL_PROFILE("SuperDropletPC::initializeParticlesUniformDistribution");

    const Real* dx = Geom(m_lev).CellSize();
    const Real* plo = Geom(m_lev).ProbLo();

    const Real cell_volume = dx[0]*dx[1]*dx[2];

    // number of super-droplets per cell
    const int num_sd_per_cell = m_ppc_init;
    // number of particles per cell (not super-droplets)
    const int num_par_per_cell = (m_numdens_init < 0 ? 1 : (int) (m_numdens_init * cell_volume));
    if (!num_par_per_cell) {
        return;
    }

    // initial radius
    const Real par_radius = (m_radius_init < 0 ? 1.0e-15 : m_radius_init);

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

            int num_to_add = num_par_per_cell;
            int num_par_per_superdroplet = amrex::max(1,(int)std::ceil(num_par_per_cell/num_sd_per_cell));

            for (int n = 0; n < num_sd_per_cell; n++) {
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
                int multiplicity = std::min( num_to_add, num_par_per_superdroplet );
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

                p.rdata(SuperDropletsRealIdxAoS::mass) = 0.0;

                p.idata(SuperDropletsIntIdxAoS::i) = iv[0];
                p.idata(SuperDropletsIntIdxAoS::j) = iv[1];
                p.idata(SuperDropletsIntIdxAoS::k) = iv[2];

                particle_tile.push_back(p);

                particle_tile.push_back_real(SuperDropletsRealIdxSoA::radius, par_radius);
                particle_tile.push_back_real(SuperDropletsRealIdxSoA::sol_mass, 0.0);
                particle_tile.push_back_int(SuperDropletsIntIdxSoA::multiplicity, multiplicity);
           }
           AMREX_ALWAYS_ASSERT(num_to_add == 0);
        }

    }

    return;
}

/*! Set super-droplets multiplicity and radius in the domain from a given condensate mass density */
void SuperDropletPC::SetAttributes (MultiFab& a_mass_density /*!< mass density of condensate */)
{
    BL_PROFILE("SuperDropletPC::SetAttributes");

    const auto plo = Geom(m_lev).ProbLoArray();
    const auto dx = Geom(m_lev).CellSize();
    const auto dxi = Geom(m_lev).InvCellSizeArray();
    const auto domain = Geom(m_lev).Domain();

    const Real cell_volume = dx[0]*dx[1]*dx[2];

    // number of super-droplets per cell
    const int num_sd_per_cell = m_ppc_init;
    // initial radius
    const Real par_radius = (m_radius_init < 0 ? 1.0e-15 : m_radius_init);
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

        auto* mult_ptr = soa.GetIntData(SuperDropletsIntIdxSoA::multiplicity).data();
        auto* radius_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::radius).data();

        const int max_multiplicity = m_max_multiplicity;
        auto mass_density = a_mass_density[pti.index()].array();

        ParallelForRNG(n, [=] AMREX_GPU_DEVICE (int i, RandomEngine const& rnd_engine)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            auto iv = getParticleCell(p, plo, dxi, domain);

            const Real mass_cell = mass_density(iv[0],iv[1],iv[2],0) * cell_volume;
            const Real mass_per_sd = mass_cell / num_sd_per_cell;

            p.rdata(SuperDropletsRealIdxAoS::mass) = mass_per_sd;

            int multiplicity = Random_int(max_multiplicity, rnd_engine)+1; // 1-100
            Real radius_cubed = mass_per_sd / ((4.0/3.0)*PI*multiplicity*mat_density);
            Real radius = std::exp( std::log(radius_cubed) / 3.0 );

            mult_ptr[i] = multiplicity;
            radius_ptr[i] = radius;

        });

    }

    return;
}

#endif

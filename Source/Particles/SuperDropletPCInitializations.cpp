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

    return;
}

/*! Initialize particles in domain */
void SuperDropletPC::InitializeParticles (const std::unique_ptr<amrex::MultiFab>& a_height_ptr)
{
    BL_PROFILE("SuperDropletPC::initializeParticles");

    if (m_initialization_type == SuperDropletInitializations::init_uniform) {
        initializeParticlesUniformDistribution( a_height_ptr );
    } else if (m_initialization_type == SuperDropletInitializations::init_null) {
        initializeParticlesNull( a_height_ptr );
    } else {
        amrex::Print() << "Error: " << m_initialization_type
                        << " is not a valid initialization for "
                        << m_name << " particle species.\n";
        amrex::Error("See error message!");
    }
    return;
}

/*! Uniform distribution: the number of particles per grid cell is specified
 *  by "initial_particles_per_cell", and they are randomly distributed. */
void SuperDropletPC::initializeParticlesUniformDistribution (const std::unique_ptr<amrex::MultiFab>& a_height_ptr)
{
    BL_PROFILE("SuperDropletPC::initializeParticlesUniformDistribution");

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
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

        auto& particle_tile = DefineAndReturnParticleTile(lev, mfi);

        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
            std::vector<Real> rnd_reals(m_ppc_init*3);
            amrex::FillRandom(rnd_reals.data(), rnd_reals.size());
            for (int n = 0; n < m_ppc_init; n++) {
                Real r[3] = {rnd_reals[3*n], rnd_reals[3*n+1], rnd_reals[3*n+2]};
                Real v[3] = {0.0, 0.0, 0.0};

                Real x = plo[0] + (iv[0] + r[0])*dx[0];
                Real y = plo[1] + (iv[1] + r[1])*dx[1];
                Real z = plo[2] + (iv[2] + r[2])*dx[2];
                if (height_ptr) {
                    z = (*height_ptr)(iv)
                         + r[2] * (   (*height_ptr)(iv + IntVect(AMREX_D_DECL(0, 0, 1)))
                                    - (*height_ptr)(iv) );
                }

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

                particle_tile.push_back_real(SuperDropletsRealIdxSoA::radius, 0.0);
                particle_tile.push_back_int(SuperDropletsIntIdxSoA::multiplicity, 1);
           }
        }

    }

    return;
}

#endif

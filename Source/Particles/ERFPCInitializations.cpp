#include <AMReX_Random.H>
#include <ERFPC.H>

#ifdef ERF_USE_PARTICLES

#include <AMReX_ParmParse.H>

using namespace amrex;

/*! Read inputs from file */
void ERFPC::readInputs()
{
    BL_PROFILE("ERFPC::readInputs");

    ParmParse pp(m_name);

    m_initialization_type = ERFParticleInitializations::init_default;
    pp.query("initial_distribution_type", m_initialization_type);

    m_ppc_init = 1;
    pp.query("inital_particles_per_cell", m_ppc_init);

    m_advect_w_flow = (m_name == ERFParticleNames::tracers ? true : false);
    pp.query("advect_with_flow", m_advect_w_flow);

    m_advect_w_gravity = (m_name == ERFParticleNames::hydro ? true : false);
    pp.query("advect_with_gravity", m_advect_w_gravity);

    return;
}

/*! Initialize particles in domain */
void ERFPC::InitializeParticles(const std::unique_ptr<amrex::MultiFab>& a_height_ptr)
{
    BL_PROFILE("ERFPC::initializeParticles");

    if (m_initialization_type == ERFParticleInitializations::init_default) {
        initializeParticlesDefault( a_height_ptr );
    } else if (m_initialization_type == ERFParticleInitializations::init_uniform) {
        initializeParticlesUniformDistribution( a_height_ptr );
    } else {
        amrex::Print() << "Error: " << m_initialization_type
                        << " is not a valid initialization for "
                        << m_name << " particle species.\n";
        amrex::Error("See error message!");
    }
    return;
}

/*! Default particle initialization */
void ERFPC::initializeParticlesDefault(const std::unique_ptr<amrex::MultiFab>& a_height_ptr)
{
    BL_PROFILE(" ERFPC::initializeParticlesDefault");

    if (m_name == ERFParticleNames::tracers) {
        initializeParticlesDefaultTracersWoA( a_height_ptr );
    } else if (m_name == ERFParticleNames::hydro) {
        initializeParticlesDefaultHydro( a_height_ptr );
    }
}

/*! Default initialization for tracer particles for WoA case (ref: AA) */
void ERFPC::initializeParticlesDefaultTracersWoA(const std::unique_ptr<amrex::MultiFab>& a_height_ptr)
{
    BL_PROFILE("ERFPC::initializeParticlesDefaultTracersWoA");

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        Gpu::HostVector<ParticleType> host_particles;

        if (a_height_ptr) {

            const auto& height = (*a_height_ptr)[mfi];
            const FArrayBox* height_ptr = nullptr;
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
            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
                if (iv[0] == 3) {
                    Real r[3] = {0.5, 0.5, 0.5};  // this means place at cell center
                    Real v[3] = {0.0, 0.0, 0.0};  // with 0 initial velocity

                    Real x = plo[0] + (iv[0] + r[0])*dx[0];
                    Real y = plo[1] + (iv[1] + r[1])*dx[1];
                    Real z = (*height_ptr)(iv)
                              + r[2] * (   (*height_ptr)(iv + IntVect(AMREX_D_DECL(0, 0, 1)))
                                         - (*height_ptr)(iv) );

                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.pos(2) = z;

                    p.rdata(ERFParticlesRealIdx::vx) = v[0];
                    p.rdata(ERFParticlesRealIdx::vy) = v[1];
                    p.rdata(ERFParticlesRealIdx::vz) = v[2];

                    p.idata(ERFParticlesIntIdx::i) = iv[0];
                    p.idata(ERFParticlesIntIdx::j) = iv[1];
                    p.idata(ERFParticlesIntIdx::k) = iv[2];

                    host_particles.push_back(p);
               }
            }

        } else {

            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
                if (iv[0] == 3) {
                    Real r[3] = {0.5, 0.5, 0.5};  // this means place at cell center
                    Real v[3] = {0.0, 0.0, 0.0};  // with 0 initial velocity

                    Real x = plo[0] + (iv[0] + r[0])*dx[0];
                    Real y = plo[1] + (iv[1] + r[1])*dx[1];
                    Real z = plo[2] + (iv[2] + r[2])*dx[2];

                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.pos(2) = z;

                    p.rdata(ERFParticlesRealIdx::vx) = v[0];
                    p.rdata(ERFParticlesRealIdx::vy) = v[1];
                    p.rdata(ERFParticlesRealIdx::vz) = v[2];

                    p.idata(ERFParticlesIntIdx::i) = iv[0];
                    p.idata(ERFParticlesIntIdx::j) = iv[1];
                    p.idata(ERFParticlesIntIdx::k) = iv[2];

                    host_particles.push_back(p);
               }
            }

        }

        auto& particles = GetParticles(lev);
        auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);
    }

    return;
}

/*! Default initialization for hydro particles (ref: AA) */
void ERFPC::initializeParticlesDefaultHydro(const std::unique_ptr<amrex::MultiFab>& a_height_ptr)
{
    BL_PROFILE("ERFPC::initializeParticlesDefaultHydro");

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        Gpu::HostVector<ParticleType> host_particles;

        if (a_height_ptr) {

            const auto& height = (*a_height_ptr)[mfi];
            const FArrayBox* height_ptr = nullptr;
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
            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
                // This is a random choice to put them above the ground and let them fall
                if (iv[2] == 13) {
                    Real r[3] = {0.5, 0.5, 0.5};  // this means place at cell center
                    Real v[3] = {0.0, 0.0, 0.0};  // with 0 initial velocity

                    Real x = plo[0] + (iv[0] + r[0])*dx[0];
                    Real y = plo[1] + (iv[1] + r[1])*dx[1];
                    Real z = (*height_ptr)(iv)
                             + r[2]*(   (*height_ptr)(iv + IntVect(AMREX_D_DECL(0, 0, 1)))
                                      - (*height_ptr)(iv) );

                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.pos(2) = z;

                    p.rdata(ERFParticlesRealIdx::vx) = v[0];
                    p.rdata(ERFParticlesRealIdx::vy) = v[1];
                    p.rdata(ERFParticlesRealIdx::vz) = v[2];

                    p.idata(ERFParticlesIntIdx::i) = iv[0];
                    p.idata(ERFParticlesIntIdx::j) = iv[1];
                    p.idata(ERFParticlesIntIdx::k) = iv[2];

                    host_particles.push_back(p);
               }
            }

        } else {

            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
                // This is a random choice to put them above the ground and let them fall
                if (iv[2] == 23) {
                    Real r[3] = {0.5, 0.5, 0.5};  // this means place at cell center
                    Real v[3] = {0.0, 0.0, 0.0};  // with 0 initial velocity

                    Real x = plo[0] + (iv[0] + r[0])*dx[0];
                    Real y = plo[1] + (iv[1] + r[1])*dx[1];
                    Real z = plo[2] + (iv[2] + r[2])*dx[2];

                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.pos(2) = z;

                    p.rdata(ERFParticlesRealIdx::vx) = v[0];
                    p.rdata(ERFParticlesRealIdx::vy) = v[1];
                    p.rdata(ERFParticlesRealIdx::vz) = v[2];

                    p.idata(ERFParticlesIntIdx::i) = iv[0];
                    p.idata(ERFParticlesIntIdx::j) = iv[1];
                    p.idata(ERFParticlesIntIdx::k) = iv[2];

                    host_particles.push_back(p);
               }
            }

        }

        auto& particles = GetParticles(lev);
        auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);
    }
}

/*! Uniform distribution: the number of particles per grid cell is specified
 *  by "initial_particles_per_cell", and they are randomly distributed. */
void ERFPC::initializeParticlesUniformDistribution(const std::unique_ptr<amrex::MultiFab>& a_height_ptr)
{
    BL_PROFILE("ERFPC::initializeParticlesUniformDistribution");

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        Gpu::HostVector<ParticleType> host_particles;

        if (a_height_ptr) {

            const auto& height = (*a_height_ptr)[mfi];
            const FArrayBox* height_ptr = nullptr;
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
            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
                std::vector<Real> rnd_reals(m_ppc_init*3);
                amrex::FillRandom(rnd_reals.data(), rnd_reals.size());
                for (int n = 0; n < m_ppc_init; n++) {
                    Real r[3] = {rnd_reals[3*n], rnd_reals[3*n+1], rnd_reals[3*n+2]};
                    Real v[3] = {0.0, 0.0, 0.0};

                    Real x = plo[0] + (iv[0] + r[0])*dx[0];
                    Real y = plo[1] + (iv[1] + r[1])*dx[1];
                    Real z = (*height_ptr)(iv)
                              + r[2] * (   (*height_ptr)(iv + IntVect(AMREX_D_DECL(0, 0, 1)))
                                         - (*height_ptr)(iv) );

                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.pos(2) = z;

                    p.rdata(ERFParticlesRealIdx::vx) = v[0];
                    p.rdata(ERFParticlesRealIdx::vy) = v[1];
                    p.rdata(ERFParticlesRealIdx::vz) = v[2];

                    p.idata(ERFParticlesIntIdx::i) = iv[0];
                    p.idata(ERFParticlesIntIdx::j) = iv[1];
                    p.idata(ERFParticlesIntIdx::k) = iv[2];

                    host_particles.push_back(p);
               }
            }

        } else {

            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
                std::vector<Real> rnd_reals(m_ppc_init*3);
                amrex::FillRandom(rnd_reals.data(), rnd_reals.size());
                for (int n = 0; n < m_ppc_init; n++) {
                    Real r[3] = {rnd_reals[3*n], rnd_reals[3*n+1], rnd_reals[3*n+2]};
                    Real v[3] = {0.0, 0.0, 0.0};

                    Real x = plo[0] + (iv[0] + r[0])*dx[0];
                    Real y = plo[1] + (iv[1] + r[1])*dx[1];
                    Real z = plo[2] + (iv[2] + r[2])*dx[2];

                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.pos(2) = z;

                    p.rdata(ERFParticlesRealIdx::vx) = v[0];
                    p.rdata(ERFParticlesRealIdx::vy) = v[1];
                    p.rdata(ERFParticlesRealIdx::vz) = v[2];

                    p.idata(ERFParticlesIntIdx::i) = iv[0];
                    p.idata(ERFParticlesIntIdx::j) = iv[1];
                    p.idata(ERFParticlesIntIdx::k) = iv[2];

                    host_particles.push_back(p);
               }
            }

        }

        auto& particles = GetParticles(lev);
        auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);
    }

    return;
}

#endif

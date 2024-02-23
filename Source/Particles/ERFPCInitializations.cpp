#include <ERFPC.H>

#ifdef ERF_USE_PARTICLES

#include <AMReX_ParmParse.H>

using namespace amrex;

/*! Read inputs from file */
void ERFPC::readInputs ()
{
    BL_PROFILE("ERFPC::readInputs");

    ParmParse pp(m_name);

    m_initialization_type = ERFParticleInitializations::init_default;
    pp.query("initial_distribution_type", m_initialization_type);

    m_ppc_init = 1;
    pp.query("initial_particles_per_cell", m_ppc_init);

    m_advect_w_flow = (m_name == ERFParticleNames::tracers ? true : false);
    pp.query("advect_with_flow", m_advect_w_flow);

    m_advect_w_gravity = (m_name == ERFParticleNames::hydro ? true : false);
    pp.query("advect_with_gravity", m_advect_w_gravity);

    return;
}

/*! Initialize particles in domain */
void ERFPC::InitializeParticles (const std::unique_ptr<amrex::MultiFab>& a_height_ptr)
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
void ERFPC::initializeParticlesDefault (const std::unique_ptr<amrex::MultiFab>& a_height_ptr)
{
    BL_PROFILE(" ERFPC::initializeParticlesDefault");

    if (m_name == ERFParticleNames::tracers) {
        initializeParticlesDefaultTracersWoA( a_height_ptr );
    } else if (m_name == ERFParticleNames::hydro) {
        initializeParticlesDefaultHydro( a_height_ptr );
    }
}

/*! Default initialization for tracer particles for WoA case (ref: AA) */
void ERFPC::initializeParticlesDefaultTracersWoA (const std::unique_ptr<amrex::MultiFab>& a_height_ptr)
{
    BL_PROFILE("ERFPC::initializeParticlesDefaultTracersWoA");

    const int lev = 0;
    const auto dx = Geom(lev).CellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();

    iMultiFab num_particles( ParticleBoxArray(lev),
                             ParticleDistributionMap(lev),
                             1, 0 );
    num_particles.setVal(0);
    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();
        auto num_particles_arr = num_particles[mfi].array();
        ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i == 3) {
                num_particles_arr(i,j,k) = 1;
            }
        });
    }

    iMultiFab offsets( ParticleBoxArray(lev),
                       ParticleDistributionMap(lev),
                       1, 0 );
    offsets.setVal(0);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();

        int np = 0;
        {
            int ncell = num_particles[mfi].numPts();
            const int* in = num_particles[mfi].dataPtr();
            int* out = offsets[mfi].dataPtr();
            np = Scan::PrefixSum<int>( ncell,
                                       [=] AMREX_GPU_DEVICE (int i) -> int { return in[i]; },
                                       [=] AMREX_GPU_DEVICE (int i, int const &x) { out[i] = x; },
                                       Scan::Type::exclusive,
                                       Scan::retSum );
        }
        auto offset_arr = offsets[mfi].array();

        auto& particle_tile = DefineAndReturnParticleTile(lev, mfi);
        particle_tile.resize(np);
        auto aos = &particle_tile.GetArrayOfStructs()[0];
        auto& soa = particle_tile.GetStructOfArrays();
        auto* vx_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::vx).data();
        auto* vy_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::vy).data();
        auto* vz_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::vz).data();
        auto* mass_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::mass).data();

        auto my_proc = ParallelDescriptor::MyProc();
        Long pid;
        {
            pid = ParticleType::NextID();
            ParticleType::NextID(pid+np);
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( static_cast<Long>(pid + np) < LastParticleID,
                                          "Error: overflow on particle id numbers!" );

        if (a_height_ptr) {

            const auto height_arr = (*a_height_ptr)[mfi].array();
            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == 3) {
                    Real r[3] = {0.5, 0.5, 0.5};  // this means place at cell center
                    Real v[3] = {0.0, 0.0, 0.0};  // with 0 initial velocity

                    Real x = plo[0] + (i + r[0])*dx[0];
                    Real y = plo[1] + (j + r[1])*dx[1];
                    Real z = height_arr(i,j,k)
                              + r[2] * (   height_arr(i,j,k+1)
                                         - height_arr(i,j,k) );

                    auto& p = aos[offset_arr(i,j,k)];
                    p.id()  = pid + offset_arr(i,j,k);
                    p.cpu() = my_proc;
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.pos(2) = z;

                    p.idata(ERFParticlesIntIdxAoS::k) = k;

                    vx_ptr[offset_arr(i,j,k)] = v[0];
                    vy_ptr[offset_arr(i,j,k)] = v[1];
                    vz_ptr[offset_arr(i,j,k)] = v[2];

                    mass_ptr[offset_arr(i,j,k)] = 0.0;
               }
            });

        } else {

            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == 3) {
                    Real r[3] = {0.5, 0.5, 0.5};  // this means place at cell center
                    Real v[3] = {0.0, 0.0, 0.0};  // with 0 initial velocity

                    Real x = plo[0] + (i + r[0])*dx[0];
                    Real y = plo[1] + (j + r[1])*dx[1];
                    Real z = plo[2] + (k + r[2])*dx[2];

                    auto& p = aos[offset_arr(i,j,k)];
                    p.id()  = pid + offset_arr(i,j,k);
                    p.cpu() = my_proc;
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.pos(2) = z;

                    p.idata(ERFParticlesIntIdxAoS::k) = k;

                    vx_ptr[offset_arr(i,j,k)] = v[0];
                    vy_ptr[offset_arr(i,j,k)] = v[1];
                    vz_ptr[offset_arr(i,j,k)] = v[2];

                    mass_ptr[offset_arr(i,j,k)] = 0.0;
               }
            });

        }

    }

    return;
}

/*! Default initialization for hydro particles (ref: AA) */
void ERFPC::initializeParticlesDefaultHydro (const std::unique_ptr<amrex::MultiFab>& a_height_ptr)
{
    BL_PROFILE("ERFPC::initializeParticlesDefaultHydro");

    const int lev = 0;
    const auto dx = Geom(lev).CellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();

    iMultiFab num_particles( ParticleBoxArray(lev),
                             ParticleDistributionMap(lev),
                             1, 0 );
    num_particles.setVal(0);
    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();
        auto num_particles_arr = num_particles[mfi].array();
        ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // This is a random choice to put them above the ground and let them fall
            if (k == 13) {
                num_particles_arr(i,j,k) = 1;
            }
        });
    }

    iMultiFab offsets( ParticleBoxArray(lev),
                       ParticleDistributionMap(lev),
                       1, 0 );
    offsets.setVal(0);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();

        int np = 0;
        {
            int ncell = num_particles[mfi].numPts();
            const int* in = num_particles[mfi].dataPtr();
            int* out = offsets[mfi].dataPtr();
            np = Scan::PrefixSum<int>( ncell,
                                       [=] AMREX_GPU_DEVICE (int i) -> int { return in[i]; },
                                       [=] AMREX_GPU_DEVICE (int i, int const &x) { out[i] = x; },
                                       Scan::Type::exclusive,
                                       Scan::retSum );
        }
        auto offset_arr = offsets[mfi].array();

        auto& particle_tile = DefineAndReturnParticleTile(lev, mfi);
        particle_tile.resize(np);
        auto aos = &particle_tile.GetArrayOfStructs()[0];
        auto& soa = particle_tile.GetStructOfArrays();
        auto* vx_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::vx).data();
        auto* vy_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::vy).data();
        auto* vz_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::vz).data();
        auto* mass_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::mass).data();

        auto my_proc = ParallelDescriptor::MyProc();
        Long pid;
        {
            pid = ParticleType::NextID();
            ParticleType::NextID(pid+np);
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( static_cast<Long>(pid + np) < LastParticleID,
                                          "Error: overflow on particle id numbers!" );

        if (a_height_ptr) {

            const auto height_arr = (*a_height_ptr)[mfi].array();
            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k == 13) {
                    Real r[3] = {0.5, 0.5, 0.5};  // this means place at cell center
                    Real v[3] = {0.0, 0.0, 0.0};  // with 0 initial velocity

                    Real x = plo[0] + (i + r[0])*dx[0];
                    Real y = plo[1] + (j + r[1])*dx[1];
                    Real z = height_arr(i,j,k)
                              + r[2] * (   height_arr(i,j,k+1)
                                         - height_arr(i,j,k) );

                    auto& p = aos[offset_arr(i,j,k)];
                    p.id()  = pid + offset_arr(i,j,k);
                    p.cpu() = my_proc;
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.pos(2) = z;

                    p.idata(ERFParticlesIntIdxAoS::k) = k;

                    vx_ptr[offset_arr(i,j,k)] = v[0];
                    vy_ptr[offset_arr(i,j,k)] = v[1];
                    vz_ptr[offset_arr(i,j,k)] = v[2];

                    mass_ptr[offset_arr(i,j,k)] = 0.0;
               }
            });

        } else {

            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k == 13) {
                    Real r[3] = {0.5, 0.5, 0.5};  // this means place at cell center
                    Real v[3] = {0.0, 0.0, 0.0};  // with 0 initial velocity

                    Real x = plo[0] + (i + r[0])*dx[0];
                    Real y = plo[1] + (j + r[1])*dx[1];
                    Real z = plo[2] + (k + r[2])*dx[2];

                    auto& p = aos[offset_arr(i,j,k)];
                    p.id()  = pid + offset_arr(i,j,k);
                    p.cpu() = my_proc;
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.pos(2) = z;

                    p.idata(ERFParticlesIntIdxAoS::k) = k;

                    vx_ptr[offset_arr(i,j,k)] = v[0];
                    vy_ptr[offset_arr(i,j,k)] = v[1];
                    vz_ptr[offset_arr(i,j,k)] = v[2];

                    mass_ptr[offset_arr(i,j,k)] = 0.0;
               }
            });

        }

    }

    return;
}

/*! Uniform distribution: the number of particles per grid cell is specified
 *  by "initial_particles_per_cell", and they are randomly distributed. */
void ERFPC::initializeParticlesUniformDistribution (const std::unique_ptr<amrex::MultiFab>& a_height_ptr)
{
    BL_PROFILE("ERFPC::initializeParticlesUniformDistribution");

    const int lev = 0;
    const auto dx = Geom(lev).CellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();

    int particles_per_cell = m_ppc_init;

    iMultiFab num_particles( ParticleBoxArray(lev),
                             ParticleDistributionMap(lev),
                             1, 0 );
    num_particles.setVal(0);
    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();
        auto num_particles_arr = num_particles[mfi].array();
        ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            num_particles_arr(i,j,k) = particles_per_cell;
        });
    }

    iMultiFab offsets( ParticleBoxArray(lev),
                       ParticleDistributionMap(lev),
                       1, 0 );
    offsets.setVal(0);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();

        int np = 0;
        {
            int ncell = num_particles[mfi].numPts();
            const int* in = num_particles[mfi].dataPtr();
            int* out = offsets[mfi].dataPtr();
            np = Scan::PrefixSum<int>( ncell,
                                       [=] AMREX_GPU_DEVICE (int i) -> int { return in[i]; },
                                       [=] AMREX_GPU_DEVICE (int i, int const &x) { out[i] = x; },
                                       Scan::Type::exclusive,
                                       Scan::retSum );
        }
        auto offset_arr = offsets[mfi].array();

        auto& particle_tile = DefineAndReturnParticleTile(lev, mfi);
        particle_tile.resize(np);
        auto aos = &particle_tile.GetArrayOfStructs()[0];
        auto& soa = particle_tile.GetStructOfArrays();
        auto* vx_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::vx).data();
        auto* vy_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::vy).data();
        auto* vz_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::vz).data();
        auto* mass_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::mass).data();

        auto my_proc = ParallelDescriptor::MyProc();
        Long pid;
        {
            pid = ParticleType::NextID();
            ParticleType::NextID(pid+np);
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( static_cast<Long>(pid + np) < LastParticleID,
                                          "Error: overflow on particle id numbers!" );

        if (a_height_ptr) {

            const auto height_arr = (*a_height_ptr)[mfi].array();
            ParallelForRNG(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                           const RandomEngine& rnd_engine) noexcept
            {
                int start = offset_arr(i,j,k);
                for (int n = start; n < start+particles_per_cell; n++) {
                    Real r[3] = {Random(rnd_engine), Random(rnd_engine), Random(rnd_engine)};
                    Real v[3] = {0.0, 0.0, 0.0};

                    Real x = plo[0] + (i + r[0])*dx[0];
                    Real y = plo[1] + (j + r[1])*dx[1];
                    Real z = height_arr(i,j,k)
                              + r[2] * (   height_arr(i,j,k+1)
                                         - height_arr(i,j,k) );

                    auto& p = aos[n];
                    p.id()  = pid + n;
                    p.cpu() = my_proc;
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.pos(2) = z;

                    p.idata(ERFParticlesIntIdxAoS::k) = k;

                    vx_ptr[n] = v[0];
                    vy_ptr[n] = v[1];
                    vz_ptr[n] = v[2];

                    mass_ptr[n] = 0.0;
               }
            });

        } else {

            ParallelForRNG(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                           const RandomEngine& rnd_engine) noexcept
            {
                int start = offset_arr(i,j,k);
                for (int n = start; n < start+particles_per_cell; n++) {
                    Real r[3] = {Random(rnd_engine), Random(rnd_engine), Random(rnd_engine)};
                    Real v[3] = {0.0, 0.0, 0.0};

                    Real x = plo[0] + (i + r[0])*dx[0];
                    Real y = plo[1] + (j + r[1])*dx[1];
                    Real z = plo[2] + (k + r[2])*dx[2];

                    auto& p = aos[n];
                    p.id()  = pid + n;
                    p.cpu() = my_proc;
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.pos(2) = z;

                    p.idata(ERFParticlesIntIdxAoS::k) = k;

                    vx_ptr[n] = v[0];
                    vy_ptr[n] = v[1];
                    vz_ptr[n] = v[2];

                    mass_ptr[n] = 0.0;
               }
            });
        }
    }

    return;
}

#endif

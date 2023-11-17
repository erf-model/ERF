#include "HydroPC.H"

// #include <AMReX_HydroParticle_mod_K.H>
#include <ERF_Constants.H>

using namespace amrex;

void
HydroPC::
InitParticles ()
{
    BL_PROFILE("HydroPC::InitParticles");

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();

        Gpu::HostVector<ParticleType> host_particles;
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
            if (iv[2] == 23) { // This is a random choice to put them above the ground and let them fall
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

                p.rdata(HydroRealIdx::vx) = v[0];
                p.rdata(HydroRealIdx::vy) = v[1];
                p.rdata(HydroRealIdx::vz) = v[2];

                p.idata(HydroIntIdx::k) = iv[2];  // particles carry their z-index

                host_particles.push_back(p);
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

void
HydroPC::
InitParticles (const MultiFab& a_z_height)
{
    BL_PROFILE("HydroPC::InitParticles");

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const auto& height = a_z_height[mfi];
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
        Gpu::HostVector<ParticleType> host_particles;
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
            if (iv[2] == 13) { // This is a random choice to put them above the ground and let them fall
                Real r[3] = {0.5, 0.5, 0.5};  // this means place at cell center
                Real v[3] = {0.0, 0.0, 0.0};  // with 0 initial velocity

                Real x = plo[0] + (iv[0] + r[0])*dx[0];
                Real y = plo[1] + (iv[1] + r[1])*dx[1];
                Real z = (*height_ptr)(iv) + r[2]*((*height_ptr)(iv + IntVect(AMREX_D_DECL(0, 0, 1))) - (*height_ptr)(iv));

                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                p.pos(0) = x;
                p.pos(1) = y;
                p.pos(2) = z;

                p.rdata(HydroRealIdx::vx) = v[0];
                p.rdata(HydroRealIdx::vy) = v[1];
                p.rdata(HydroRealIdx::vz) = v[2];

                p.idata(HydroIntIdx::k) = iv[2];  // particles carry their z-index

                host_particles.push_back(p);
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

/*
  /brief Uses midpoint method to advance particles using umac.
*/
void
HydroPC::AdvectWithGravity (int lev, Real dt, bool use_terrain, const MultiFab& a_z_height)
{
    BL_PROFILE("HydroPC::AdvectWithGravity()");
    AMREX_ASSERT(lev >= 0 && lev < GetParticles().size());

    const auto      strttime = amrex::second();
    const Geometry& geom = m_gdb->Geom(lev);
    const Box& domain = geom.Domain();
    const auto plo = geom.ProbLoArray();
    const auto dx  = geom.CellSizeArray();
    const auto dxi = geom.InvCellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        int grid    = pti.index();
        auto& ptile = ParticlesAt(lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        const int n = aos.numParticles();
        auto *p_pbox = aos().data();

        auto zheight = use_terrain ? a_z_height[grid].array() : Array4<Real>{};

        amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }

            ParticleReal v = p.rdata(HydroRealIdx::vz);

            // Define acceleration to be (gravity minus drag) where drag is defined
            // such the particles will reach a terminal velocity of 5.0 (totally arbitrary)
            ParticleReal terminal_vel = 5.0;
            ParticleReal grav = CONST_GRAV;
            ParticleReal drag = CONST_GRAV * (v * v) / (terminal_vel*terminal_vel);

            ParticleReal half_dt = 0.5 * dt;

            // Update the particle velocity over first half of step (dt/2)
            p.rdata(HydroRealIdx::vz) -= (grav - drag) * half_dt;

            // Update the particle position over (dt)
            p.pos(2) += static_cast<ParticleReal>(ParticleReal(0.5)*dt*p.rdata(HydroRealIdx::vz));

            // Update the particle velocity over second half of step (dt/2)
            p.rdata(HydroRealIdx::vz) -= (grav - drag) * half_dt;

            // also update z-coordinate here
            IntVect iv(
                AMREX_D_DECL(int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                             int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                             p.idata(0)));
            iv[0] += domain.smallEnd()[0];
            iv[1] += domain.smallEnd()[1];
            ParticleReal zlo, zhi;
            if (use_terrain) {
                zlo = zheight(iv[0], iv[1], iv[2]);
                zhi = zheight(iv[0], iv[1], iv[2]+1);
            } else {
                zlo =  iv[2]    * dx[2];
                zhi = (iv[2]+1) * dx[2];
            }
            if (p.pos(2) > zhi) { // need to be careful here
                p.idata(0) += 1;
            } else if (p.pos(2) <= zlo) {
                p.idata(0) -= 1;
            }
        });
    }

    if (m_verbose > 1)
    {
        auto stoptime = amrex::second() - strttime;

#ifdef AMREX_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
                ParallelReduce::Max(stoptime, ParallelContext::IOProcessorNumberSub(),
                                    ParallelContext::CommunicatorSub());

                amrex::Print() << "HydroParticleContainer::AdvectWithGravity() time: " << stoptime << '\n';
#ifdef AMREX_LAZY
        });
#endif
    }
}

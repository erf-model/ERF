#include <ERFPC.H>

#ifdef ERF_USE_PARTICLES

#include <IndexDefines.H>
#include <ERF_Constants.H>
#include <AMReX_TracerParticle_mod_K.H>

using namespace amrex;

/*! Evolve particles for one time step */
void ERFPC::EvolveParticles ( int                                        a_lev,
                              Real                                       a_dt_lev,
                              Vector<Vector<MultiFab>>&                  a_flow_vars,
                              const Vector<std::unique_ptr<MultiFab>>&   a_z_phys_nd )
{
    BL_PROFILE("ERFPCPC::EvolveParticles()");

    if (m_advect_w_flow) {
        MultiFab* flow_vel( &a_flow_vars[a_lev][Vars::xvel] );
        AdvectWithFlow( flow_vel, a_lev, a_dt_lev, a_z_phys_nd[a_lev] );
    }

    if (m_advect_w_gravity) {
        AdvectWithGravity( a_lev, a_dt_lev, a_z_phys_nd[a_lev] );
    }
    return;
}

/*! Uses midpoint method to advance particles using flow velocity. */
void ERFPC::AdvectWithFlow ( MultiFab*                           a_umac,
                             int                                 a_lev,
                             Real                                a_dt,
                             const std::unique_ptr<MultiFab>&    a_z_height )
{
    BL_PROFILE("ERFPCPC::AdvectWithUmac()");
    AMREX_ASSERT(OK(a_lev, a_lev, a_umac[0].nGrow()-1));
    AMREX_ASSERT(a_lev >= 0 && a_lev < GetParticles().size());

    AMREX_D_TERM(AMREX_ASSERT(a_umac[0].nGrow() >= 1);,
                 AMREX_ASSERT(a_umac[1].nGrow() >= 1);,
                 AMREX_ASSERT(a_umac[2].nGrow() >= 1););

    AMREX_D_TERM(AMREX_ASSERT(!a_umac[0].contains_nan());,
                 AMREX_ASSERT(!a_umac[1].contains_nan());,
                 AMREX_ASSERT(!a_umac[2].contains_nan()););

    const auto      strttime = amrex::second();
    const Geometry& geom = m_gdb->Geom(a_lev);
    const Box& domain = geom.Domain();
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const auto dx  = geom.CellSizeArray();

    Vector<std::unique_ptr<MultiFab> > raii_umac(AMREX_SPACEDIM);
    Vector<MultiFab*> umac_pointer(AMREX_SPACEDIM);
    if (OnSameGrids(a_lev, a_umac[0]))
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            umac_pointer[i] = &a_umac[i];
        }
    }
    else
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
            int ng = a_umac[i].nGrow();
            raii_umac[i] = std::make_unique<MultiFab>
                (amrex::convert(m_gdb->ParticleBoxArray(a_lev), IntVect::TheDimensionVector(i)),
                 m_gdb->ParticleDistributionMap(a_lev), a_umac[i].nComp(), ng);
            umac_pointer[i] = raii_umac[i].get();
            umac_pointer[i]->ParallelCopy(a_umac[i],0,0,a_umac[i].nComp(),ng,ng);
        }
    }

    for (int ipass = 0; ipass < 2; ipass++)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti)
        {
            int grid    = pti.index();
            auto& ptile = ParticlesAt(a_lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            const int n = aos.numParticles();
            auto *p_pbox = aos().data();
            const FArrayBox* fab[AMREX_SPACEDIM] = { AMREX_D_DECL(&((*umac_pointer[0])[grid]),
                                                                  &((*umac_pointer[1])[grid]),
                                                                  &((*umac_pointer[2])[grid])) };

            //array of these pointers to pass to the GPU
            amrex::GpuArray<amrex::Array4<const Real>, AMREX_SPACEDIM>
                const umacarr {{AMREX_D_DECL((*fab[0]).array(),
                                             (*fab[1]).array(),
                                             (*fab[2]).array() )}};

            bool use_terrain = (a_z_height != nullptr);
            auto zheight = use_terrain ? (*a_z_height)[grid].array() : Array4<Real>{};

            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }
                ParticleReal v[AMREX_SPACEDIM];
                mac_interpolate(p, plo, dxi, umacarr, v);
                if (ipass == 0)
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.rdata(dim) = p.pos(dim);
                        p.pos(dim) += static_cast<ParticleReal>(ParticleReal(0.5)*a_dt*v[dim]);
                    }
                }
                else
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.pos(dim) = p.rdata(dim) + static_cast<ParticleReal>(a_dt*v[dim]);
                        p.rdata(dim) = v[dim];
                    }

                    // Update z-coordinate carried by the particle
                    IntVect iv(
                        AMREX_D_DECL(int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                                     int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                                     p.idata(ERFParticlesIntIdx::k)));
                    iv[0] += domain.smallEnd()[0];
                    iv[1] += domain.smallEnd()[1];
                    ParticleReal zlo, zhi;
                    if (use_terrain) {
                        Real lx = (p.pos(0)-plo[0])*dxi[0] - static_cast<Real>(iv[0]-domain.smallEnd()[0]);
                        Real ly = (p.pos(1)-plo[1])*dxi[1] - static_cast<Real>(iv[1]-domain.smallEnd()[1]);
                        zlo =  zheight(iv[0]  ,iv[1]  ,iv[2]  ) * (1.0-lx) * (1.0-ly) +
                               zheight(iv[0]+1,iv[1]  ,iv[2]  ) *      lx  * (1.0-ly) +
                               zheight(iv[0]  ,iv[1]+1,iv[2]  ) * (1.0-lx) * ly +
                               zheight(iv[0]+1,iv[1]+1,iv[2]  ) *      lx  * ly;
                        zhi =  zheight(iv[0]  ,iv[1]  ,iv[2]+1) * (1.0-lx) * (1.0-ly) +
                               zheight(iv[0]+1,iv[1]  ,iv[2]+1) *      lx  * (1.0-ly) +
                               zheight(iv[0]  ,iv[1]+1,iv[2]+1) * (1.0-lx) * ly +
                               zheight(iv[0]+1,iv[1]+1,iv[2]+1) *      lx  * ly;
                    } else {
                        zlo =  iv[2]    * dx[2];
                        zhi = (iv[2]+1) * dx[2];
                    }
                    if (p.pos(2) > zhi) { // need to be careful here
                        p.idata(ERFParticlesIntIdx::k) += 1;
                    } else if (p.pos(2) <= zlo) {
                        p.idata(ERFParticlesIntIdx::k) -= 1;
                    }
                }
            });
        }
    }

    if (m_verbose > 1)
    {
        auto stoptime = amrex::second() - strttime;

#ifdef AMREX_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
                ParallelReduce::Max(stoptime, ParallelContext::IOProcessorNumberSub(),
                                    ParallelContext::CommunicatorSub());

                amrex::Print() << "ERFPC::AdvectWithFlow() time: " << stoptime << '\n';
#ifdef AMREX_LAZY
        });
#endif
    }
}

void ERFPC::AdvectWithGravity (  int                                 a_lev,
                                 Real                                a_dt,
                                 const std::unique_ptr<MultiFab>&    a_z_height )
{
    BL_PROFILE("ERFPC::AdvectWithGravity()");
    AMREX_ASSERT(a_lev >= 0 && a_lev < GetParticles().size());

    const auto      strttime = amrex::second();
    const Geometry& geom = m_gdb->Geom(a_lev);
    const Box& domain = geom.Domain();
    const auto plo = geom.ProbLoArray();
    const auto dx  = geom.CellSizeArray();
    const auto dxi = geom.InvCellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti)
    {
        int grid    = pti.index();
        auto& ptile = ParticlesAt(a_lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        const int n = aos.numParticles();
        auto *p_pbox = aos().data();

        bool use_terrain = (a_z_height != nullptr);
        auto zheight = use_terrain ? (*a_z_height)[grid].array() : Array4<Real>{};

        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }

            ParticleReal v = p.rdata(ERFParticlesRealIdx::vz);

            // Define acceleration to be (gravity minus drag) where drag is defined
            // such the particles will reach a terminal velocity of 5.0 (totally arbitrary)
            ParticleReal terminal_vel = 5.0;
            ParticleReal grav = CONST_GRAV;
            ParticleReal drag = CONST_GRAV * (v * v) / (terminal_vel*terminal_vel);

            ParticleReal half_dt = 0.5 * a_dt;

            // Update the particle velocity over first half of step (a_dt/2)
            p.rdata(ERFParticlesRealIdx::vz) -= (grav - drag) * half_dt;

            // Update the particle position over (a_dt)
            p.pos(2) += static_cast<ParticleReal>(ParticleReal(0.5)*a_dt*p.rdata(ERFParticlesRealIdx::vz));

            // Update the particle velocity over second half of step (a_dt/2)
            p.rdata(ERFParticlesRealIdx::vz) -= (grav - drag) * half_dt;

            // also update z-coordinate here
            IntVect iv(
                AMREX_D_DECL(int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                             int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                             p.idata(ERFParticlesIntIdx::k)));
            iv[0] += domain.smallEnd()[0];
            iv[1] += domain.smallEnd()[1];
            ParticleReal zlo, zhi;
            if (use_terrain) {
                Real lx = (p.pos(0)-plo[0])*dxi[0] - static_cast<Real>(iv[0]-domain.smallEnd()[0]);
                Real ly = (p.pos(1)-plo[1])*dxi[1] - static_cast<Real>(iv[1]-domain.smallEnd()[1]);
                zlo =  zheight(iv[0]  ,iv[1]  ,iv[2]  ) * (1.0-lx) * (1.0-ly) +
                       zheight(iv[0]+1,iv[1]  ,iv[2]  ) *      lx  * (1.0-ly) +
                       zheight(iv[0]  ,iv[1]+1,iv[2]  ) * (1.0-lx) * ly +
                       zheight(iv[0]+1,iv[1]+1,iv[2]  ) *      lx  * ly;
                zhi =  zheight(iv[0]  ,iv[1]  ,iv[2]+1) * (1.0-lx) * (1.0-ly) +
                       zheight(iv[0]+1,iv[1]  ,iv[2]+1) *      lx  * (1.0-ly) +
                       zheight(iv[0]  ,iv[1]+1,iv[2]+1) * (1.0-lx) * ly +
                       zheight(iv[0]+1,iv[1]+1,iv[2]+1) *      lx  * ly;
            } else {
                zlo =  iv[2]    * dx[2];
                zhi = (iv[2]+1) * dx[2];
            }
            if (p.pos(2) > zhi) { // need to be careful here
                p.idata(ERFParticlesIntIdx::k) += 1;
            } else if (p.pos(2) <= zlo) {
                p.idata(ERFParticlesIntIdx::k) -= 1;
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

                amrex::Print() << "ERFPC::AdvectWithGravity() time: " << stoptime << '\n';
#ifdef AMREX_LAZY
        });
#endif
    }
}

#endif

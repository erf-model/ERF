#include <ERFPC.H>

#ifdef ERF_USE_PARTICLES

#include <ERF_IndexDefines.H>
#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include <ERF_TerrainConversion.H>
#include <AMReX_TracerParticle_mod_K.H>

using namespace amrex;

/*! Evolve particles for one time step */
void ERFPC::EvolveParticles ( int                                        a_lev,
                              Real                                       a_dt_lev,
                              Vector<Vector<MultiFab>>&                  a_flow_vars,
                              const Vector<std::unique_ptr<MultiFab>>&   a_z_phys_nd )
{
    BL_PROFILE("ERFPCPC::EvolveParticles()");

    if (m_verbose > 0) {
        Long np_total = 0;
        int finest = m_gdb->finestLevel();
        amrex::Print() << "[" << m_name << "] Evolving particles on level " << a_lev
                       << ": ";
        for (int lev = 0; lev <= finest; lev++) {
            Long np_lev = NumberOfParticlesAtLevel(lev, true, true);
            ParallelDescriptor::ReduceLongSum(np_lev);
            amrex::Print() << "L" << lev << "=" << np_lev;
            if (lev < finest) { amrex::Print() << " "; }
            np_total += np_lev;
        }
        amrex::Print() << " (total=" << np_total << ")\n";
    }

    if (m_advect_w_flow) {
        MultiFab* flow_vel( &a_flow_vars[a_lev][Vars::xvel] );
        AdvectWithFlow( flow_vel, a_lev, a_dt_lev, a_z_phys_nd[a_lev] );
    }

    if (m_advect_w_gravity) {
        AdvectWithGravity( a_lev, a_dt_lev, a_z_phys_nd[a_lev] );
    }

    // Per-level Redistribute keeps each level's particles on that level for
    // the duration of the coarse step, avoiding double-advection of particles
    // that geometrically cross between levels mid-step.  Cross-level moves
    // happen only at regrid time.  For fine levels, particles that leave the
    // level's BA must first be routed back to the coarse level.
    if (a_lev == 0) {
        Redistribute(0, 0);
    } else {
        ExtractAndRouteOORParticles(a_lev);
    }

    ComputeTemperature( a_flow_vars[a_lev][Vars::cons], a_lev, a_dt_lev, a_z_phys_nd[a_lev] );

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
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const Box&  dom = geom.Domain();
    const int   k_max = dom.bigEnd(AMREX_SPACEDIM-1) - dom.smallEnd(AMREX_SPACEDIM-1);

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
            IntVect ng = a_umac[i].nGrowVect();
            raii_umac[i] = std::make_unique<MultiFab>
                (convert(m_gdb->ParticleBoxArray(a_lev), IntVect::TheDimensionVector(i)),
                 m_gdb->ParticleDistributionMap(a_lev), a_umac[i].nComp(), ng);
            umac_pointer[i] = raii_umac[i].get();
            umac_pointer[i]->ParallelCopy(a_umac[i],0,0,a_umac[i].nComp(),ng,ng);
        }
    }

    bool periodic_in_z = (geom.isPeriodic(2));

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
            auto& soa  = ptile.GetStructOfArrays();
            const int n = aos.numParticles();
            auto *p_pbox = aos().data();

            Array<ParticleReal*,AMREX_SPACEDIM> v_ptr;
            v_ptr[0] = soa.GetRealData(ERFParticlesRealIdx::vx).data();
            v_ptr[1] = soa.GetRealData(ERFParticlesRealIdx::vy).data();
            v_ptr[2] = soa.GetRealData(ERFParticlesRealIdx::vz).data();

            const FArrayBox* fab[AMREX_SPACEDIM] = { AMREX_D_DECL(&((*umac_pointer[0])[grid]),
                                                                  &((*umac_pointer[1])[grid]),
                                                                  &((*umac_pointer[2])[grid])) };

            //array of these pointers to pass to the GPU
            GpuArray<Array4<const Real>, AMREX_SPACEDIM>
                const umacarr {{AMREX_D_DECL((*fab[0]).array(),
                                             (*fab[1]).array(),
                                             (*fab[2]).array() )}};

            auto zheight = (*a_z_height)[grid].array();

            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }

                ParticleReal v[AMREX_SPACEDIM];

                int pk = int(amrex::Math::floor((p.pos(AMREX_SPACEDIM-1) - plo[AMREX_SPACEDIM-1])
                                                * dxi[AMREX_SPACEDIM-1]));
                if (pk - 1 < zheight.begin[2] || pk + 2 >= zheight.end[2]) {
                    v[0] = 0; v[1] = 0; v[2] = 0;
                } else {
                    mac_interpolate(p, plo, dxi, umacarr, v);
                    if (amrex::isnan(v[0]) || amrex::isnan(v[1]) || amrex::isnan(v[2])) {
                        v[0] = 0; v[1] = 0; v[2] = 0;
                    }
                }

                // Advance in physical (x, y, z); pos(2) is zeta on disk so go
                // through z_from_zeta / zeta_from_z around the update.
                if (ipass == 0) {
                    const Real x0 = static_cast<Real>(p.pos(0));
                    const Real y0 = static_cast<Real>(p.pos(1));
                    const Real zeta0 = static_cast<Real>(p.pos(AMREX_SPACEDIM-1));
                    const Real z_phys0 = static_cast<Real>(ERF::ParticlePos::z_from_zeta(
                        x0, y0, zeta0, plo, dxi, zheight));
                    const Real x_h = x0 + static_cast<Real>(Real(0.5)*a_dt*v[0]);
                    const Real y_h = y0 + static_cast<Real>(Real(0.5)*a_dt*v[1]);
                    const Real z_h = z_phys0 + static_cast<Real>(Real(0.5)*a_dt*v[2]);
                    v_ptr[0][i] = static_cast<ParticleReal>(x0);
                    v_ptr[1][i] = static_cast<ParticleReal>(y0);
                    v_ptr[2][i] = static_cast<ParticleReal>(zeta0);
                    p.pos(0) = static_cast<ParticleReal>(x_h);
                    p.pos(1) = static_cast<ParticleReal>(y_h);
                    p.pos(AMREX_SPACEDIM-1) = static_cast<ParticleReal>(
                        ERF::ParticlePos::zeta_from_z(x_h, y_h, z_h, plo, dxi, zheight, k_max));
                } else {
                    const Real x0 = static_cast<Real>(v_ptr[0][i]);
                    const Real y0 = static_cast<Real>(v_ptr[1][i]);
                    const Real zeta0 = static_cast<Real>(v_ptr[2][i]);
                    const Real z_phys0 = static_cast<Real>(ERF::ParticlePos::z_from_zeta(
                        x0, y0, zeta0, plo, dxi, zheight));
                    const Real x_n = x0 + static_cast<Real>(a_dt*v[0]);
                    const Real y_n = y0 + static_cast<Real>(a_dt*v[1]);
                    const Real z_n = z_phys0 + static_cast<Real>(a_dt*v[2]);
                    p.pos(0) = static_cast<ParticleReal>(x_n);
                    p.pos(1) = static_cast<ParticleReal>(y_n);
                    p.pos(AMREX_SPACEDIM-1) = static_cast<ParticleReal>(
                        ERF::ParticlePos::zeta_from_z(x_n, y_n, z_n, plo, dxi, zheight, k_max));
                    v_ptr[0][i] = static_cast<ParticleReal>(v[0]);
                    v_ptr[1][i] = static_cast<ParticleReal>(v[1]);
                    v_ptr[2][i] = static_cast<ParticleReal>(v[2]);
                }

                // Particle crossed below the bottom: move to 0.2*dz above floor
                if (!periodic_in_z) {
                    if (p.pos(AMREX_SPACEDIM-1) < plo[AMREX_SPACEDIM-1]) {
                        p.pos(AMREX_SPACEDIM-1) = plo[AMREX_SPACEDIM-1]
                                                + Real(0.2) / dxi[AMREX_SPACEDIM-1];
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

                Print() << "ERFPC::AdvectWithFlow() time: " << stoptime << '\n';
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
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const Box&  dom = geom.Domain();
    const int   k_max = dom.bigEnd(AMREX_SPACEDIM-1) - dom.smallEnd(AMREX_SPACEDIM-1);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti)
    {
        int grid    = pti.index();
        auto& ptile = ParticlesAt(a_lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        auto& soa  = ptile.GetStructOfArrays();
        const int n = aos.numParticles();
        auto *p_pbox = aos().data();

        auto vz_ptr = soa.GetRealData(ERFParticlesRealIdx::vz).data();

        auto zheight = (*a_z_height)[grid].array();

        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }

            ParticleReal v = vz_ptr[i];

            // Define acceleration to be (gravity minus drag) where drag is defined
            // such the particles will reach a terminal velocity of Real(5.0) (totally arbitrary)
            ParticleReal terminal_vel = Real(5.0);
            ParticleReal grav = CONST_GRAV;
            ParticleReal drag = CONST_GRAV * (v * v) / (terminal_vel*terminal_vel);

            ParticleReal myhalf_dt = myhalf * a_dt;

            // Update the particle velocity over first half of step (a_dt/2)
            vz_ptr[i] -= (grav - drag) * myhalf_dt;

            // Advance in physical z, then back to zeta.
            {
                const Real x0 = static_cast<Real>(p.pos(0));
                const Real y0 = static_cast<Real>(p.pos(1));
                const Real z_phys0 = static_cast<Real>(ERF::ParticlePos::z_from_zeta(
                    x0, y0, static_cast<Real>(p.pos(AMREX_SPACEDIM-1)), plo, dxi, zheight));
                const Real z_phys_n = z_phys0 + static_cast<Real>(Real(0.5) * a_dt * vz_ptr[i]);
                p.pos(AMREX_SPACEDIM-1) = static_cast<ParticleReal>(
                    ERF::ParticlePos::zeta_from_z(x0, y0, z_phys_n, plo, dxi, zheight, k_max));
            }

            // Update the particle velocity over second half of step (a_dt/2)
            vz_ptr[i] -= (grav - drag) * myhalf_dt;
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

                Print() << "ERFPC::AdvectWithGravity() time: " << stoptime << '\n';
#ifdef AMREX_LAZY
        });
#endif
    }
}

/*! Uses midpoint method to advance particles using flow velocity. */
void ERFPC::ComputeTemperature (const MultiFab&                     a_ucons,
                                int                                 a_lev,
                                Real                                /*a_dt*/,
                                const std::unique_ptr<MultiFab>&    a_z_height )
{
    BL_PROFILE("ERFPCPC::ComputeTemperature()");
    AMREX_ASSERT(a_lev >= 0 && a_lev < GetParticles().size());

    const auto strttime = amrex::second();
    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    // Compute temperature
    MultiFab T_mf(a_ucons.boxArray(), a_ucons.DistributionMap(), 1, a_ucons.nGrowVect());
    T_mf.setVal(0.0);
    for ( MFIter mfi(a_ucons); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();
        auto states_array = a_ucons.const_array(mfi);
        auto tabs_array  = T_mf.array(mfi);
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            tabs_array(i,j,k)  = getTgivenRandRTh(states_array(i,j,k,Rho_comp),
                                                  states_array(i,j,k,RhoTheta_comp));
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti)
    {
        int grid    = pti.index();
        auto& ptile = ParticlesAt(a_lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        auto& soa  = ptile.GetStructOfArrays();
        const int n = aos.numParticles();
        auto *p_pbox = aos().data();

        auto* T_ptr = soa.GetRealData(ERFParticlesRealIdx::temperature).data();
        auto temperature_arr  = T_mf.array(grid);

        auto zheight = (*a_z_height)[grid].array();

        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }

            ParticleReal temperature = 0;
            int pk = int(amrex::Math::floor((p.pos(AMREX_SPACEDIM-1) - plo[AMREX_SPACEDIM-1])
                                            * dxi[AMREX_SPACEDIM-1]));
            if (pk >= zheight.begin[2] && pk + 2 < zheight.end[2]) {
                cic_interpolate( p, plo, dxi, temperature_arr, &temperature, 1 );
            }
            T_ptr[i] = temperature;
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

                Print() << "ERFPC::ComputeTemperature() time: " << stoptime << '\n';
#ifdef AMREX_LAZY
        });
#endif
    }
}

#endif

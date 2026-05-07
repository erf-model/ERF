#include <ERFPC.H>

#ifdef ERF_USE_PARTICLES

#include <ERF_IndexDefines.H>
#include <ERF_Constants.H>
#include <ERF_EOS.H>
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

    // Redistribute particles.  For fine levels, use
    // ExtractAndRouteOORParticles which handles partial-z refinement
    // (particles escaping the fine level's z-extent) by recomputing
    // k-indices for the target level before per-level Redistribute.
    if (a_lev == 0) {
        Redistribute(0, 0);
    } else {
        ExtractAndRouteOORParticles(a_lev, a_z_phys_nd);
    }

    // After redistribution, recompute k-indices from z-position for all
    // particles using each level's geometry.  ERFParticlesAssignor uses
    // idata(k) for grid assignment, so it must stay in sync with pos(z).
    for (int lev = 0; lev <= a_lev; lev++) {
        const auto& plev = GetParticles();
        if (lev >= static_cast<int>(plev.size())) { continue; }
        if (plev[lev].empty()) { continue; }

        const Geometry& glev = m_gdb->Geom(lev);
        const auto plo_lev = glev.ProbLoArray();
        const auto dxi_lev = glev.InvCellSizeArray();
        const int k_max_lev = glev.Domain().bigEnd(AMREX_SPACEDIM-1);

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            auto& aos = ParticlesAt(lev, pti).GetArrayOfStructs();
            const int np = aos.numParticles();
            auto* p_pbox = aos().data();
            auto zheight = (*a_z_phys_nd[lev])[pti.index()].array();

            ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                ParticleType& p = p_pbox[i];
                if (p.id() > 0) {
                    int k_guess = int(amrex::Math::floor(
                        (Real(p.pos(AMREX_SPACEDIM-1)) - plo_lev[AMREX_SPACEDIM-1])
                        * dxi_lev[AMREX_SPACEDIM-1]));
                    p.idata(ERFParticlesIntIdxAoS::k) = amrex::max(0, amrex::min(k_guess, k_max_lev));
                    update_location_idata(p, plo_lev, dxi_lev, zheight);
                }
            });
        }
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

    auto domlo = lbound(geom.Domain());

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
            v_ptr[0] = soa.GetRealData(ERFParticlesRealIdxSoA::vx).data();
            v_ptr[1] = soa.GetRealData(ERFParticlesRealIdxSoA::vy).data();
            v_ptr[2] = soa.GetRealData(ERFParticlesRealIdxSoA::vz).data();

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

                // Partial-z tiles: zero v for particles outside stencil (k-1..k+2); Redistribute will move them.
                int pk = p.idata(ERFParticlesIntIdxAoS::k);
                if (pk - 1 < zheight.begin[2] || pk + 2 >= zheight.end[2]) {
                    v[0] = 0; v[1] = 0; v[2] = 0;
                } else {
                    mac_interpolate_mapped_z(p, plo, dxi, umacarr, zheight, v);
                    if (amrex::isnan(v[0]) || amrex::isnan(v[1]) || amrex::isnan(v[2])) {
                        v[0] = 0; v[1] = 0; v[2] = 0;
                    }
                }

                if (ipass == 0) {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        v_ptr[dim][i] = p.pos(dim);
                        p.pos(dim) += static_cast<ParticleReal>(ParticleReal(0.5)*a_dt*v[dim]);
                    }
                } else {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.pos(dim) = v_ptr[dim][i] + static_cast<ParticleReal>(a_dt*v[dim]);
                        v_ptr[dim][i] = v[dim];
                    }
                }

                // Update z-coordinate carried by the particle
                update_location_idata(p,plo,dxi,zheight);

                // If the particle crossed below the bottom surface, move it up to 0.2*dz above the surface
                if (!periodic_in_z) {
                    if (p.idata(ERFParticlesIntIdxAoS::k) < 0) {
                        int ii = domlo.x + int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0]));
                        int jj = domlo.y + int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1]));
                        int kk = 0;

                        // Update the stored particle location
                        p.idata(ERFParticlesIntIdxAoS::k) = kk;

                        Real lx = (p.pos(0)-plo[0])*dxi[0] - static_cast<amrex::Real>(ii);
                        Real ly = (p.pos(1)-plo[1])*dxi[1] - static_cast<amrex::Real>(jj);
                        auto zlo = zheight(ii  ,jj  ,kk  ) * (one-lx) * (one-ly) +
                                   zheight(ii+1,jj  ,kk  ) *      lx  * (one-ly) +
                                   zheight(ii  ,jj+1,kk  ) * (one-lx) * ly +
                                   zheight(ii+1,jj+1,kk  ) *      lx  * ly;
                        auto zhi = zheight(ii  ,jj  ,kk+1) * (one-lx) * (one-ly) +
                                   zheight(ii+1,jj  ,kk+1) *      lx  * (one-ly) +
                                   zheight(ii  ,jj+1,kk+1) * (one-lx) * ly +
                                   zheight(ii+1,jj+1,kk+1) *      lx  * ly;
                        p.pos(2) = zlo + Real(0.2) * (zhi - zlo);
                    } // k < 0
                } // !periodic
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

        auto vz_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::vz).data();

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

            // Update the particle position over (a_dt)
            p.pos(2) += static_cast<ParticleReal>(ParticleReal(0.5)*a_dt*vz_ptr[i]);

            // Update the particle velocity over second half of step (a_dt/2)
            vz_ptr[i] -= (grav - drag) * myhalf_dt;

            // also update z-coordinate here
            update_location_idata(p,plo,dxi,zheight);
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

        auto* T_ptr = soa.GetRealData(ERFParticlesRealIdxSoA::temperature).data();
        auto temperature_arr  = T_mf.array(grid);

        auto zheight = (*a_z_height)[grid].array();

        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }

            ParticleReal temperature = 0;
            int pk = p.idata(ERFParticlesIntIdxAoS::k);
            if (pk >= zheight.begin[2] && pk + 2 < zheight.end[2]) {
                cic_interpolate_mapped_z( p, plo, dxi, temperature_arr, zheight, &temperature, 1 );
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

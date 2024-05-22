#include "SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include <AMReX_TracerParticle_mod_K.H>
#include "IndexDefines.H"
#include "TerminalVelocity.H"

using namespace amrex;

/*! Evolve particles for one time step */
void SuperDropletPC::AdvectParticles ( int                   a_lev,
                                       Real                  a_dt,
                                       const MultiFab* const a_flow_vel,
                                       const MultiFab&       a_density,
                                       const MultiFab&       a_pressure,
                                       const MultiFab&       a_temperature,
                                       const Vector<MFPtr>&  a_z_phys_nd )
{
    BL_PROFILE("SuperDropletPC::AdvectParticles()");
    const MFPtr& z_height = a_z_phys_nd[a_lev];

    AMREX_ASSERT(OK(a_lev, a_lev, a_flow_vel[0].nGrow()-1));
    AMREX_ASSERT(a_lev >= 0 && a_lev < GetParticles().size());
    AMREX_D_TERM(AMREX_ASSERT(a_flow_vel[0].nGrow() >= 1);,
                 AMREX_ASSERT(a_flow_vel[1].nGrow() >= 1);,
                 AMREX_ASSERT(a_flow_vel[2].nGrow() >= 1););

    AMREX_D_TERM(AMREX_ASSERT(!a_flow_vel[0].contains_nan());,
                 AMREX_ASSERT(!a_flow_vel[1].contains_nan());,
                 AMREX_ASSERT(!a_flow_vel[2].contains_nan()););

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto domain = geom.Domain();
    const int k_lo = domain.smallEnd(2);

    const bool advect_w_flow = m_advect_w_flow;
    const bool advect_w_gravity = m_advect_w_gravity;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {
        int grid    = pti.index();
        auto& ptile = ParticlesAt(a_lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        auto& soa  = ptile.GetStructOfArrays();
        const int n = aos.numParticles();
        auto *p_pbox = aos().data();

        Array<ParticleReal*,AMREX_SPACEDIM> v_ptr;
        v_ptr[0] = soa.GetRealData(SuperDropletsRealIdxSoA::vx).data();
        v_ptr[1] = soa.GetRealData(SuperDropletsRealIdxSoA::vy).data();
        v_ptr[2] = soa.GetRealData(SuperDropletsRealIdxSoA::vz).data();

        const FArrayBox* fab[AMREX_SPACEDIM] = { AMREX_D_DECL(&(a_flow_vel[0][grid]),
                                                              &(a_flow_vel[1][grid]),
                                                              &(a_flow_vel[2][grid])) };
        amrex::GpuArray<amrex::Array4<const Real>, AMREX_SPACEDIM>
            const umacarr {{AMREX_D_DECL((*fab[0]).array(),
                                         (*fab[1]).array(),
                                         (*fab[2]).array() )}};

        const auto& density_arr = a_density[grid].array();
        const auto& temperature_arr = a_temperature[grid].array();
        const auto& pressure_arr = a_pressure[grid].array();

        bool use_terrain = (z_height != nullptr);
        auto zheight = use_terrain ? (*z_height)[grid].array() : Array4<Real>{};

        int rt_offset = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::radius).data();
        auto* vterm_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::term_vel).data();
        auto* mult_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::multiplicity).data();
        auto* supdrop_mass_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::sd_mass).data();

        TerminalVelocity<ParticleReal> term_vel { m_vapour_mat->density() };
        auto term_vel_type = m_term_vel_type;

        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (mult_ptr[i] == 0) { return; }

            for (int ipass = 0; ipass < 2; ipass++) {

                ParticleReal v[AMREX_SPACEDIM];
                v[0] = v[1] = v[2] = 0.0;

                if (advect_w_flow) {
                    if (use_terrain) {
                        mac_interpolate_mapped_z(p, plo, dxi, umacarr, zheight, v);
                    } else {
                        mac_interpolate(p, plo, dxi, umacarr, v);
                    }
                }

                ParticleReal terminal_vel = 0.0;
                if (advect_w_gravity) {
                    if (term_vel_type == SDTerminalVelocityType::AtlasUlbrich) {
                        terminal_vel = term_vel.AtlasUlbrich( radius_ptr[i] );
                    } else if (term_vel_type == SDTerminalVelocityType::CloudRainShima) {
                        ParticleReal density, pressure, temperature;
                        if (use_terrain) {
                            cic_interpolate_mapped_z( p, plo, dxi, density_arr, zheight, &density, 1 );
                            cic_interpolate_mapped_z( p, plo, dxi, pressure_arr, zheight, &pressure, 1 );
                            cic_interpolate_mapped_z( p, plo, dxi, temperature_arr, zheight, &temperature, 1 );
                        } else {
                            cic_interpolate( p, plo, dxi, density_arr, &density, 1 );
                            cic_interpolate( p, plo, dxi, pressure_arr, &pressure, 1 );
                            cic_interpolate( p, plo, dxi, temperature_arr, &temperature, 1 );
                        }
                        terminal_vel = term_vel.CloudRainShima( radius_ptr[i],
                                                                density,
                                                                pressure,
                                                                temperature );
                    }
                }
                vterm_ptr[i] = terminal_vel;
                v[2] -= terminal_vel;

                if (ipass == 0) {

                    for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                        v_ptr[dim][i] = p.pos(dim);
                        p.pos(dim) += static_cast<ParticleReal>(ParticleReal(0.5)*a_dt*v[dim]);
                    }
                    // Update z-coordinate carried by the particle
                    update_location_idata(p,plo,dxi,zheight);

                } else {

                    for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                        p.pos(dim) = v_ptr[dim][i] + static_cast<ParticleReal>(a_dt*v[dim]);
                        v_ptr[dim][i] = v[dim];
                    }
                    // Update z-coordinate carried by the particle
                    update_location_idata(p,plo,dxi,zheight);

                }
            }

            // check for ground impact
            auto iv = getParticleCell(p, plo, dxi, domain);
            auto z_ground = plo[2];
            if (use_terrain) { z_ground = zheight(iv[0],iv[1],k_lo); }
            if (p.pos(2) < z_ground) {
                p.pos(2) = z_ground - 0.01*dx[2];
                v_ptr[0][i] = v_ptr[1][i] = v_ptr[2][i] = 0.0;
                v_ptr[2][i] = vterm_ptr[i] = 0.0;
                mult_ptr[i] = 0.0;
                supdrop_mass_ptr[i] = 0.0;
            }

        });
        Gpu::synchronize();
    }

    Redistribute();
}

#endif

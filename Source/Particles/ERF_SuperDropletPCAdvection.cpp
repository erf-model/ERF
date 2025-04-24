#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_IndexDefines.H"
#include "ERF_TerminalVelocity.H"

using namespace amrex;

/*! Evolve particles for one time step */
void SuperDropletPC::AdvectParticles ( int                   a_lev,
                                       Real                  /*a_time*/,
                                       Real                  a_dt,
                                       const MultiFab* const a_flow_vel,
                                       const MultiFab&       a_density,
                                       const MultiFab&       a_pressure,
                                       const MultiFab&       a_temperature,
                                       const Vector<MFPtr>&  a_z_phys_nd,
                                       const BCTypeArr&      a_bctypes )
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

    const auto is_periodic = geom.isPeriodic();
    auto is_periodic_z = is_periodic[2];

    const bool advect_w_flow = m_advect_w_flow;
    const bool advect_w_gravity = m_advect_w_gravity;
    const Real rho_w = m_vapour_mat->density();
    const int num_aerosols = m_num_aerosols;

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

        auto zheight = (*z_height)[grid].array();

        int rt_offset = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::radius).data();
        auto* vterm_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::term_vel).data();
        auto* mult_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::multiplicity).data();

        SDAerosolMassArr aerosol_mass_ptrs;
        Gpu::DeviceVector<ParticleReal> aerosol_density(num_aerosols);
        Gpu::DeviceVector<int> aerosol_solubility(num_aerosols);
        Vector<ParticleReal> aerosol_density_h(num_aerosols);
        Vector<int> aerosol_solubility_h(num_aerosols);
        for (int i = 0; i < num_aerosols; i++) {
            aerosol_mass_ptrs[i] = soa.GetRealData(   rt_offset
                                                    + SuperDropletsRealIdxSoA_RT::ncomps
                                                    + i ).data();
            aerosol_density_h[i] = m_aerosol_mat[i]->density();
            aerosol_solubility_h[i] = static_cast<int>(m_aerosol_mat[i]->isSoluble());
        }
        Gpu::copy(  Gpu::hostToDevice,
                    aerosol_density_h.begin(),
                    aerosol_density_h.end(),
                    aerosol_density.begin() );
        Gpu::copy(  Gpu::hostToDevice,
                    aerosol_solubility_h.begin(),
                    aerosol_solubility_h.end(),
                    aerosol_solubility.begin() );

        TerminalVelocity<ParticleReal> term_vel { m_vapour_mat->density() };
        auto term_vel_type = m_term_vel_type;

        auto aero_rho_arr = aerosol_density.data();
        auto aero_sol_arr = aerosol_solubility.data();

        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (mult_ptr[i] == 0) { return; }

            ParticleReal v[AMREX_SPACEDIM];
            v[0] = v[1] = v[2] = 0.0;

            if (is_periodic_z) {
                mac_interpolate(p, plo, dxi, umacarr, v);
            } else {
                mac_interpolate_mapped_z(p, plo, dxi, umacarr, zheight, v);
            }

            // compute effective radius
            ParticleReal r_eff = 0.0;
            {
                ParticleReal m_w = 4.0/3.0 * PI * rho_w
                                   * radius_ptr[i]*radius_ptr[i]*radius_ptr[i];
                ParticleReal m_s = 0.0;
                ParticleReal m_p = 0.0;
                ParticleReal rho_p = 0.0;
                for (int j = 0; j < num_aerosols; j++) {
                    if (aero_sol_arr[j]) {
                        m_s += aerosol_mass_ptrs[j][i];
                    } else {
                        m_p += aerosol_mass_ptrs[j][i];
                        rho_p += aero_rho_arr[j]*aerosol_mass_ptrs[j][i];
                    }
                }
                if (m_p > 0.0) { rho_p /= m_p; }
                else           { rho_p = 1.0; }
                auto m_t = m_w + m_s + (rho_w/rho_p)*m_p;
                r_eff = std::cbrt(m_t / (4.0/3.0*PI*rho_w));
            }

            ParticleReal terminal_vel = 0.0;
            if (term_vel_type == SDTerminalVelocityType::AtlasUlbrich) {
                terminal_vel = term_vel.AtlasUlbrich( r_eff );
            } else if (term_vel_type == SDTerminalVelocityType::RogersYau) {
                terminal_vel = term_vel.RogersYau( r_eff );
            } else if (term_vel_type == SDTerminalVelocityType::CloudRainShima) {
                ParticleReal density, pressure, temperature;
                if (is_periodic_z) {
                    cic_interpolate( p, plo, dxi, density_arr, &density, 1 );
                    cic_interpolate( p, plo, dxi, pressure_arr, &pressure, 1 );
                    cic_interpolate( p, plo, dxi, temperature_arr, &temperature, 1 );
                } else {
                    cic_interpolate_mapped_z( p, plo, dxi, density_arr, zheight, &density, 1 );
                    cic_interpolate_mapped_z( p, plo, dxi, pressure_arr, zheight, &pressure, 1 );
                    cic_interpolate_mapped_z( p, plo, dxi, temperature_arr, zheight, &temperature, 1 );
                }
                terminal_vel = term_vel.CloudRainShima( r_eff,
                                                        density,
                                                        pressure,
                                                        temperature );
            }

            for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                v_ptr[dim][i] = v[dim];
            }
            vterm_ptr[i] = terminal_vel;

            if (advect_w_flow) {
                for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                    p.pos(dim) += static_cast<ParticleReal>(a_dt*v[dim]);
                }
            }
            if (advect_w_gravity) {
                p.pos(AMREX_SPACEDIM-1) -= static_cast<ParticleReal>(a_dt*terminal_vel);
            }

            // Update z-coordinate carried by the particle
            update_location_idata(p,plo,dxi,zheight);

        });
        Gpu::synchronize();
    }

    applyBoundaryTreatment(a_lev, a_z_phys_nd, a_bctypes);
    Redistribute();
}

#endif

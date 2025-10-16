#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_IndexDefines.H"
#include "ERF_TerminalVelocity.H"

using namespace amrex;

/*! Evolve particles for one time step */
void SuperDropletPC::AdvectParticles ( int                   a_lev,
                                       Real                  a_time,
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
    const bool prescribed_advection = m_prescribed_advection;

    const int num_sp = m_num_species;
    const int num_ae = m_num_aerosols;

    const int idx_w = m_idx_w;
    const int idx_i = m_idx_i;

    const auto vterm_type_w = m_term_vel_type_w;
    const auto vterm_type_i = m_term_vel_type_i;

    Real rho_w = m_species_mat[m_idx_w]->m_density;
    Real rho_i = DBL_MAX;
    if (idx_i >= 0) { rho_i = m_species_mat[m_idx_i]->m_density; }

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

        int rtoff_i = SuperDropletsIntIdxSoA::ncomps;
        auto* active_ptr = soa.GetIntData(rtoff_i+SuperDropletsIntIdxSoA_RT::active).data();
        int rtoff_r = SuperDropletsRealIdxSoA::ncomps;
        auto* vterm_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::term_vel).data();
        auto* mult_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::multiplicity).data();

        const ParticleReal *a_ptr(nullptr), *c_ptr(nullptr);
        if (idx_i >= 0) {
            a_ptr = soa.GetRealData(idx_ice_a(num_ae,num_sp)).data();
            c_ptr = soa.GetRealData(idx_ice_c(num_ae,num_sp)).data();
        }

        SDSpeciesMassArr sp_mass_ptrs;
        Gpu::DeviceVector<ParticleReal> sp_density(num_sp);
        Gpu::DeviceVector<int> sp_solubility(num_sp);
        {
            Vector<ParticleReal> sp_density_h(num_sp);
            Vector<int> sp_solubility_h(num_sp);
            for (int i = 0; i < num_sp; i++) {
                sp_mass_ptrs[i] = soa.GetRealData(idx_s(i,num_ae,num_sp)).data();
                sp_density_h[i] = m_species_mat[i]->m_density;
                sp_solubility_h[i] = static_cast<int>(m_species_mat[i]->m_is_soluble);
            }
            Gpu::copy(  Gpu::hostToDevice,
                        sp_density_h.begin(),
                        sp_density_h.end(),
                        sp_density.begin() );
            Gpu::copy(  Gpu::hostToDevice,
                        sp_solubility_h.begin(),
                        sp_solubility_h.end(),
                        sp_solubility.begin() );
        }

        SDAerosolMassArr ae_mass_ptrs;
        Gpu::DeviceVector<ParticleReal> ae_density(num_ae);
        Gpu::DeviceVector<int> ae_solubility(num_ae);
        {
            Vector<ParticleReal> ae_density_h(num_ae);
            Vector<int> ae_solubility_h(num_ae);
            for (int i = 0; i < num_ae; i++) {
                ae_mass_ptrs[i] = soa.GetRealData(idx_a(i,num_ae,num_sp)).data();
                ae_density_h[i] = m_aerosol_mat[i]->m_density;
                ae_solubility_h[i] = static_cast<int>(m_aerosol_mat[i]->m_is_soluble);
            }
            Gpu::copy(  Gpu::hostToDevice,
                        ae_density_h.begin(),
                        ae_density_h.end(),
                        ae_density.begin() );
            Gpu::copy(  Gpu::hostToDevice,
                        ae_solubility_h.begin(),
                        ae_solubility_h.end(),
                        ae_solubility.begin() );
        }

        TerminalVelocity<ParticleReal> term_vel { rho_w, rho_i };

        auto sp_rho_arr = sp_density.data();
        auto sp_sol_arr = sp_solubility.data();
        auto ae_rho_arr = ae_density.data();
        auto ae_sol_arr = ae_solubility.data();

        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (active_ptr[i] == 0) { return; }

            ParticleReal v[AMREX_SPACEDIM];
            v[0] = v[1] = v[2] = 0.0;

            if (is_periodic_z) {
                mac_interpolate(p, plo, dxi, umacarr, v);
            } else {
                mac_interpolate_mapped_z(p, plo, dxi, umacarr, zheight, v);
            }

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

            if (prescribed_advection) {
               if (a_time < 600) {
                   v[2] = 2.0*sin(PI*a_time/600)/density;
               } else {
                   v[2] = 0.0;
               }
            }

            // compute effective radius
            auto r_eff = SD_effective_radius( i, idx_w,
                                              rho_w,
                                              num_sp, num_ae,
                                              sp_sol_arr, ae_sol_arr,
                                              sp_mass_ptrs, ae_mass_ptrs,
                                              sp_rho_arr, ae_rho_arr );

            ParticleReal terminal_vel = 0.0;
            auto par_phase = SD_phase(i, idx_w, idx_i, sp_mass_ptrs);
            if (par_phase == SDPhase::water) {
                if (vterm_type_w == SDTerminalVelocityType::AtlasUlbrich) {
                    terminal_vel = term_vel.AtlasUlbrich( r_eff );
                } else if (vterm_type_w == SDTerminalVelocityType::RogersYau) {
                    terminal_vel = term_vel.RogersYau( r_eff );
                } else if (vterm_type_w == SDTerminalVelocityType::CloudRainShima) {
                    terminal_vel = term_vel.CloudRainShima( r_eff,
                                                            density,
                                                            pressure,
                                                            temperature );
                } else {
                    amrex::Abort("Invalid option for water droplet terminal velocity model");
                }

            } else if (par_phase == SDPhase::ice) {
                AMREX_ALWAYS_ASSERT(idx_i >= 0);
                if (vterm_type_i == SDTerminalVelocityType::AtlasUlbrich) {
                    terminal_vel = term_vel.AtlasUlbrich( r_eff );
                } else if (vterm_type_i == SDTerminalVelocityType::RogersYau) {
                    terminal_vel = term_vel.RogersYau( r_eff );
                } else if (vterm_type_i == SDTerminalVelocityType::IceBohm) {
                    auto m_total = SD_total_mass( i, num_sp, num_ae, sp_mass_ptrs, ae_mass_ptrs);
                    terminal_vel = term_vel.IceBohm( m_total,
                                                     a_ptr[i], c_ptr[i],
                                                     ice_rho(a_ptr[i],c_ptr[i],sp_mass_ptrs[idx_i][i]),
                                                     density, temperature );
                } else {
                    amrex::Abort("Invalid option for ice particle terminal velocity model");
                }
            } else {
                amrex::Abort("Unknown value for particle phase (must be ice or water)");
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

#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_IndexDefines.H"
#include "ERF_TerminalVelocity.H"
#include "ERF_InterpolationUtils.H"

using namespace amrex;
using namespace SDPCDefn;

namespace {
    /*! \brief Traits for particle advection process */
    struct AdvectionTraits : SDProcess::DefaultTraits {
        static constexpr bool needs_velocity  = true;
        static constexpr bool needs_term_vel  = true;
        static constexpr bool needs_ice_axes  = true;  // for IceBohm terminal velocity
    };
}

/*! Evolve particles for one time step */
void SuperDropletPC::AdvectParticles ( int                   a_lev,
                                       Real                  a_time,
                                       Real                  a_dt,
                                       const MultiFab* const a_flow_vel,
                                       const MultiFab&       a_density,
                                       const MultiFab&       a_pressure,
                                       const MultiFab&       a_temperature,
                                       const Vector<MFPtr>&  a_z_phys_nd,
                                       const BCTypeArr&      a_bctypes,
                                       const bool            a_recycle )
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

    // Build process context using helper method
    const auto ctx = buildProcessContext(a_lev);

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto is_periodic_z = geom.isPeriodic(2);

    const bool advect_w_flow = m_advect_w_flow;
    const bool advect_w_gravity = m_advect_w_gravity;
    const bool prescribed_advection = m_prescribed_advection;

    const auto vterm_type_w = m_term_vel_type_w;
    const auto vterm_type_i = m_term_vel_type_i;

    Real rho_i = DBL_MAX;
    if (ctx.idx_ice >= 0) { rho_i = m_species_mat[m_idx_i]->m_density; }

    // Terminal velocity calculator (shared across tiles)
    TerminalVelocity<ParticleReal> term_vel { ctx.rho_water, rho_i };

    forEachParticleTile<AdvectionTraits>(a_lev, ctx,
        [&](ParIterType& pti, int grid, ParticleType* p_pbox,
            const SDProcess::ParticlePointers& ptrs,
            const SDProcess::ProcessContext& ctx)
    {
        auto zheight = (*z_height)[grid].array();

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

        ParallelFor(ptrs.num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (ptrs.active_ptr[i] == 0) { return; }

            ParticleReal v[AMREX_SPACEDIM];
            v[0] = v[1] = v[2] = 0.0;

            if (is_periodic_z) {
                mac_interpolate(p, ctx.plo, ctx.dxi, umacarr, v);
            } else {
                mac_interpolate_mapped_z(p, ctx.plo, ctx.dxi, umacarr, zheight, v);
            }

            // Define field values array to store interpolated results
            ParticleReal field_values[3]; // density, pressure, temperature

            // Define array of field arrays to interpolate from
            const Array4<const Real> field_arrays[3] = {
                density_arr,
                pressure_arr,
                temperature_arr
            };

            // Use the interpolation helper function to interpolate all fields at once
            ERF::Interpolation::interpolateFields(
                p, ctx.plo, ctx.dxi, field_arrays, field_values, 3,
                is_periodic_z ? 1 : 0, is_periodic_z ? nullptr : &zheight
            );

            // Extract interpolated values
            ParticleReal density = field_values[0];
            ParticleReal pressure = field_values[1];
            ParticleReal temperature = field_values[2];

            if (prescribed_advection) {
               if (a_time < 600) {
                   v[2] = 2.0*sin(PI*a_time/600)/density;
               } else {
                   v[2] = 0.0;
               }
            }

            // compute effective radius
            auto r_eff = SD_effective_radius( i, ctx.idx_water,
                                              ctx.rho_water,
                                              ctx.num_species, ctx.num_aerosols,
                                              ptrs.sp_sol_arr, ptrs.ae_sol_arr,
                                              ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs,
                                              ptrs.sp_rho_arr, ptrs.ae_rho_arr );

            ParticleReal terminal_vel = 0.0;
            auto par_phase = SD_phase(i, ctx.idx_water, ctx.idx_ice, ptrs.sp_mass_ptrs);
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
                AMREX_ALWAYS_ASSERT(ctx.idx_ice >= 0);
                if (vterm_type_i == SDTerminalVelocityType::AtlasUlbrich) {
                    terminal_vel = term_vel.AtlasUlbrich( r_eff );
                } else if (vterm_type_i == SDTerminalVelocityType::RogersYau) {
                    terminal_vel = term_vel.RogersYau( r_eff );
                } else if (vterm_type_i == SDTerminalVelocityType::IceBohm) {
                    auto m_total = SD_total_mass( i, ctx.num_species, ctx.num_aerosols,
                                                  ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs);
                    terminal_vel = term_vel.IceBohm( m_total,
                                                     ptrs.a_ptr[i], ptrs.c_ptr[i],
                                                     ice_rho(ptrs.a_ptr[i],ptrs.c_ptr[i],
                                                             ptrs.sp_mass_ptrs[ctx.idx_ice][i]),
                                                     density, temperature );
                } else {
                    amrex::Abort("Invalid option for ice particle terminal velocity model");
                }
            } else {
                amrex::Abort("Unknown value for particle phase (must be ice or water)");
            }

            for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                ptrs.v_ptr[dim][i] = v[dim];
            }
            ptrs.vterm_ptr[i] = terminal_vel;

            if (advect_w_flow) {
                for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                    p.pos(dim) += static_cast<ParticleReal>(a_dt*v[dim]);
                }
            }
            if (advect_w_gravity) {
                p.pos(AMREX_SPACEDIM-1) -= static_cast<ParticleReal>(a_dt*terminal_vel);
            }

            // Update z-coordinate carried by the particle
            update_location_idata(p,ctx.plo,ctx.dxi,zheight);

        });
        Gpu::synchronize();
    }); // end forEachParticleTile

    applyBoundaryTreatment(a_lev, a_z_phys_nd, a_bctypes, a_recycle);
    Redistribute();
}

#endif

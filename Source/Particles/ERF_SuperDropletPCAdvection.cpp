#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_IndexDefines.H"
#include "ERF_TerminalVelocity.H"
#include "ERF_InterpolationUtils.H"

using namespace amrex;
using namespace SDPCDefn;

/*! \brief Field indices for advection interpolation */
AMREX_ENUM(InterpFieldsAdv, density, pressure, temperature, NUM_FIELDS);

/*! \brief Advect superdroplet particles for one time step
 * \param[in] a_lev AMR level
 * \param[in] a_time Current simulation time
 * \param[in] a_dt Timestep size for advection
 * \param[in] a_flow_vel Array of face-based velocities
 * \param[in] a_density Density field
 * \param[in] a_pressure Pressure field
 * \param[in] a_temperature Temperature field
 * \param[in] a_z_phys_nd Array of terrain heights
 * \param[in] a_bctypes Array of boundary condition types
 * \param[in] a_recycle Flag to enable particle recycling
 */
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

    const auto ctx = buildProcessContext(a_lev);
    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto is_periodic_z = geom.isPeriodic(2);

    const bool advect_w_flow = m_advect_w_flow;
    const bool advect_w_gravity = m_advect_w_gravity;
    const bool prescribed_advection = m_prescribed_advection;

    const auto vterm_type_w = m_term_vel_type_w;

    // Terminal velocity calculator (shared across tiles)
    TerminalVelocity<ParticleReal> term_vel { ctx.rho_water };

    forEachParticleTile(a_lev, ctx,
        [&](ParIterType& /*pti*/, int grid, ParticleType* p_pbox,
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
            v[0] = v[1] = v[2] = zero;

            if (is_periodic_z) {
                mac_interpolate(p, ctx.plo, ctx.dxi, umacarr, v);
            } else {
                mac_interpolate_mapped_z(p, ctx.plo, ctx.dxi, umacarr, zheight, v);
            }

            // Interpolate density, pressure, temperature at particle position
            constexpr int nf = static_cast<int>(InterpFieldsAdv::NUM_FIELDS);
            ParticleReal fv[nf];
            const Array4<const Real> fa[nf] = {
                density_arr, pressure_arr, temperature_arr
            };
            ERF::Interpolation::interpolateFields(
                p, ctx.plo, ctx.dxi, fa, fv, nf,
                is_periodic_z ? 1 : 0, is_periodic_z ? nullptr : &zheight
            );
            const auto density     = fv[static_cast<int>(InterpFieldsAdv::density)];
            const auto pressure    = fv[static_cast<int>(InterpFieldsAdv::pressure)];
            const auto temperature = fv[static_cast<int>(InterpFieldsAdv::temperature)];

            if (prescribed_advection) {
               if (a_time < 600) {
                   v[2] = two*sin(PI*a_time/600)/density;
               } else {
                   v[2] = zero;
               }
            }

            // compute effective radius
            auto r_eff = SD_effective_radius( i, ctx.idx_water,
                                              ctx.rho_water,
                                              ctx.num_species, ctx.num_aerosols,
                                              ptrs.sp_sol_arr, ptrs.ae_sol_arr,
                                              ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs,
                                              ptrs.sp_rho_arr, ptrs.ae_rho_arr );

            ParticleReal terminal_vel = zero;
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
                amrex::Abort("Invalid option for terminal velocity model");
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

    // After redistribution, update k-index for particles that moved to new tiles.
    // This is needed because update_location_idata skips particles outside the
    // local tile bounds (to avoid out-of-bounds array access). After Redistribute,
    // those particles are now on the correct tile and can update their k-index.
    forEachParticleTile(a_lev, ctx,
        [&](ParIterType& /*pti*/, int grid, ParticleType* p_pbox,
            const SDProcess::ParticlePointers& ptrs,
            const SDProcess::ProcessContext& ctx)
    {
        auto zheight = (*z_height)[grid].array();

        ParallelFor(ptrs.num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() > 0) {
                update_location_idata(p, ctx.plo, ctx.dxi, zheight);
            }
        });
        Gpu::synchronize();
    });
}

#endif

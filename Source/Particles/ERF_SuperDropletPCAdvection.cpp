#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_IndexDefines.H"
#include "ERF_TerminalVelocity.H"
#include "ERF_InterpolationUtils.H"
#include "ERF_TerrainConversion.H"

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
    AMREX_ALWAYS_ASSERT((a_lev >= 0) && (a_lev <= finestLevel()));

    const MFPtr& z_height = a_z_phys_nd[a_lev];

    AMREX_ASSERT(OK(a_lev, a_lev, a_flow_vel[0].nGrow()-1));
    AMREX_D_TERM(AMREX_ASSERT(a_flow_vel[0].nGrow() >= 1);,
                 AMREX_ASSERT(a_flow_vel[1].nGrow() >= 1);,
                 AMREX_ASSERT(a_flow_vel[2].nGrow() >= 1););
    AMREX_D_TERM(AMREX_ASSERT(!a_flow_vel[0].contains_nan());,
                 AMREX_ASSERT(!a_flow_vel[1].contains_nan());,
                 AMREX_ASSERT(!a_flow_vel[2].contains_nan()););

    const auto ctx = buildProcessContext(a_lev);
    const Geometry& geom = m_gdb->Geom(a_lev);
    const Box& dom = geom.Domain();
    const int  k_max = dom.bigEnd(AMREX_SPACEDIM-1) - dom.smallEnd(AMREX_SPACEDIM-1);

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

            mac_interpolate(p, ctx.plo, ctx.dxi, umacarr, v);
            if (amrex::isnan(v[0]) || amrex::isnan(v[1]) || amrex::isnan(v[2])) {
                v[0] = zero; v[1] = zero; v[2] = zero;
            }

            // Interpolate density, pressure, temperature at particle position
            constexpr int nf = static_cast<int>(InterpFieldsAdv::NUM_FIELDS);
            ParticleReal fv[nf];
            const Array4<const Real> fa[nf] = {
                density_arr, pressure_arr, temperature_arr
            };
            ERF::Interpolation::interpolateFields(
                p, ctx.plo, ctx.dxi, fa, fv, nf
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

            // Advance in physical (x, y, z); pos(2) is zeta on disk so go
            // through z_from_zeta / zeta_from_z around the update.
            if (advect_w_flow || advect_w_gravity) {
                const Real x0 = static_cast<Real>(p.pos(0));
                const Real y0 = static_cast<Real>(p.pos(1));
                const Real z_phys0 = static_cast<Real>(ERF::ParticlePos::z_from_zeta(
                    x0, y0, p.pos(AMREX_SPACEDIM-1), ctx.plo, ctx.dxi, zheight));
                Real x_n = x0;
                Real y_n = y0;
                Real z_n = z_phys0;
                if (advect_w_flow) {
                    x_n += static_cast<Real>(a_dt * v[0]);
                    y_n += static_cast<Real>(a_dt * v[1]);
                    z_n += static_cast<Real>(a_dt * v[AMREX_SPACEDIM-1]);
                }
                if (advect_w_gravity) {
                    z_n -= static_cast<Real>(a_dt * terminal_vel);
                }
                int qi = int(amrex::Math::floor((x_n - ctx.plo[0]) * ctx.dxi[0]));
                int qj = int(amrex::Math::floor((y_n - ctx.plo[1]) * ctx.dxi[1]));
                int qk = int(amrex::Math::floor((z_n - ctx.plo[AMREX_SPACEDIM-1])
                                                * ctx.dxi[AMREX_SPACEDIM-1]));
                if (qi     < zheight.begin[0] || qi + 1 >= zheight.end[0] ||
                    qj     < zheight.begin[1] || qj + 1 >= zheight.end[1] ||
                    qk     < zheight.begin[2] || qk + 1 >= zheight.end[2]) { return; }
                p.pos(0) = static_cast<ParticleReal>(x_n);
                p.pos(1) = static_cast<ParticleReal>(y_n);
                p.pos(AMREX_SPACEDIM-1) = static_cast<ParticleReal>(
                    ERF::ParticlePos::zeta_from_z(x_n, y_n, z_n, ctx.plo, ctx.dxi, zheight, k_max));
            }
        });
        Gpu::synchronize();
    }); // end forEachParticleTile

    applyBoundaryTreatment(a_lev, a_z_phys_nd, a_bctypes, a_recycle);
    Redistribute();
}

#endif

#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_IndexDefines.H"
#include "ERF_TerminalVelocity.H"
#include "ERF_InterpolationUtils.H"

using namespace amrex;
using namespace SDPCDefn;

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

        const ParticleReal *a_ptr(nullptr), *c_ptr(nullptr);
        if (idx_i >= 0) {
            a_ptr = soa.GetRealData(idx_ice_a(num_ae,num_sp)).data();
            c_ptr = soa.GetRealData(idx_ice_c(num_ae,num_sp)).data();
        }

        SDSpeciesMassArr sp_mass_ptrs;
        SDAerosolMassArr ae_mass_ptrs;
        setupMassPointers(soa, sp_mass_ptrs, ae_mass_ptrs);

        // Get pointers to persistent device data
        const ParticleReal* sp_rho_arr = getSpeciesDensitiesDevice();
        const int* sp_sol_arr = getSpeciesSolubilitiesDevice();
        const ParticleReal* ae_rho_arr = getAerosolDensitiesDevice();
        const int* ae_sol_arr = getAerosolSolubilitiesDevice();

        TerminalVelocity<ParticleReal> term_vel { rho_w, rho_i };

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
                p, plo, dxi, field_arrays, field_values, 3,
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

    applyBoundaryTreatment(a_lev, a_z_phys_nd, a_bctypes, a_recycle);
    Redistribute();
}

#endif

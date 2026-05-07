#include "ERF_SuperDropletPC.H"
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuBuffer.H>
#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_Constants.H"
#include "ERF_SuperDropletPCMassChange.H"
#include "ERF_InterpolationUtils.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;
using namespace SDPCDefn;

/*! \brief Field indices for liquid-vapour interpolation */
AMREX_ENUM(InterpFieldsLV, e_sat, sat_ratio, temperature, pressure, NUM_FIELDS);

/*! Compute mass change of particles due to evaporation and condensation */
void SuperDropletPC::MassChange ( int                                         a_lev,
                                  Real                                        a_dt,
                                  const Species::Name&                        a_vap_name,
                                  const MultiFab&                             a_temperature,
                                  const MultiFab&                             a_pressure,
                                  const MultiFab&                             a_sat_pressure,
                                  const MultiFab&                             a_sat_ratio,
                                  const Vector<std::unique_ptr<MultiFab>>&    a_z_phys_nd,
                                  const bool                                  a_is_water )
{
    using namespace SDMassChangeUtils_LV;

    BL_PROFILE("SuperDropletPC::MassChange()");
    AMREX_ASSERT( a_lev == m_lev );

    const auto ctx = buildProcessContext(a_lev);

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto is_periodic_z = geom.isPeriodic(2);

    const std::unique_ptr<MultiFab>& z_height = a_z_phys_nd[a_lev];

    // Find vapour species index
    int idx_vap = -1;
    Real mat_density = -1;
    for (int i = 0; i < m_num_species; i++) {
        if (m_species_mat[i]->m_name == a_vap_name) {
            idx_vap = i;
            mat_density = m_species_mat[i]->m_density;
        }
    }
    AMREX_ALWAYS_ASSERT(idx_vap >= 0);
    AMREX_ALWAYS_ASSERT(mat_density >= 0);
    const MaterialProperties vapour_mat(*(m_species_mat[idx_vap]));
    const MaterialPropertiesCore& vapour_mat_core = vapour_mat;

    const bool log_unconverged = m_mass_change_logging;
    [[maybe_unused]] FILE* file_handle = m_mass_change_log;
    auto cfl = m_mass_change_cfl;
    auto ti_choice = m_mass_change_ti;

    // Solver setup (shared across tiles)
    dRsqdt<ParticleReal> drsqdt{ vapour_mat.m_lat_vap,
                                 therco,
                                 vapour_mat.m_Rv,
                                 mat_density };

    NewtonSolver< dRsqdt<ParticleReal>, ParticleReal > newton_solver{ drsqdt,
                                                                      m_newton_rtol,
                                                                      m_newton_atol,
                                                                      m_newton_stol,
                                                                      m_newton_maxits };

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    dRdt<ParticleReal> drdt{ vapour_mat.m_lat_vap,
                             therco,
                             vapour_mat.m_Rv,
                             mat_density };
    constexpr int rtoff_r = SuperDropletsRealIdxSoA::ncomps;
#endif

    forEachParticleTile(a_lev, ctx,
        [&](ParIterType& pti, int grid, ParticleType* p_pbox,
            const SDProcess::ParticlePointers& ptrs,
            const SDProcess::ProcessContext& ctx)
    {
        amrex::ignore_unused(pti);
        auto zheight = (*z_height)[grid].array();

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        auto& soa = ParticlesAt(a_lev, pti).GetStructOfArrays();
        auto* condt_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::cond_tendency).data();
#endif

        const auto& sat_pressure_arr = a_sat_pressure[grid].array();
        const auto& sat_ratio_arr = a_sat_ratio[grid].array();
        const auto& temperature_arr = a_temperature[grid].array();
        const auto& pressure_arr = a_pressure[grid].array();

        Gpu::Buffer<Long> unconverged_particles({0});
        auto* unconverged_particles_ptr = unconverged_particles.data();

        AMREX_ASSERT_WITH_MESSAGE( ti_choice == SDMassChangeTIMethod::RK4 ||
                                   ti_choice == SDMassChangeTIMethod::RK3BS ||
                                   ti_choice == SDMassChangeTIMethod::BE ||
                                   ti_choice == SDMassChangeTIMethod::CN ||
                                   ti_choice == SDMassChangeTIMethod::DIRK2,
                                   "ERROR: invalid time integrator choice!" );

        ParallelFor(ptrs.num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (ptrs.active_ptr[i] == 0) { return; }

            // Interpolate saturation pressure, saturation ratio, temperature, pressure
            constexpr int nf = static_cast<int>(InterpFieldsLV::NUM_FIELDS);
            ParticleReal fv[nf];
            const Array4<const Real> fa[nf] = {
                sat_pressure_arr, sat_ratio_arr, temperature_arr, pressure_arr
            };
            ERF::Interpolation::interpolateFields(
                p, ctx.plo, ctx.dxi, fa, fv, nf,
                is_periodic_z ? 1 : 0, is_periodic_z ? nullptr : &zheight
            );
            const auto e_sat       = fv[static_cast<int>(InterpFieldsLV::e_sat)];
            const auto sat_ratio   = fv[static_cast<int>(InterpFieldsLV::sat_ratio)];
            const auto temperature = fv[static_cast<int>(InterpFieldsLV::temperature)];
            const auto pressure    = fv[static_cast<int>(InterpFieldsLV::pressure)];

            ParticleReal solute_moles = zero;
            if (a_is_water) {
                for (int j = 0; j < ctx.num_species; j++) {
                    if (j != idx_vap) {
                        solute_moles += (ptrs.sp_mass_ptrs[j][i]*ptrs.sp_ion_arr[j]/ptrs.sp_mw_arr[j]);
                    }
                }
                for (int j = 0; j < ctx.num_aerosols; j++) {
                    solute_moles += (ptrs.ae_mass_ptrs[j][i]*ptrs.ae_ion_arr[j]/ptrs.ae_mw_arr[j]);
                }
            }

            auto coeff_curv = vapour_mat_core.coeffCurv(temperature);
            auto coeff_sol = vapour_mat_core.coeffVPSolute();
            auto coeff_moldiff = vapour_mat_core.coeffMolecularDiffusion(temperature, pressure);

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
            if (a_is_water) {
                condt_ptr[i] = drdt( ptrs.radius_ptr[i],
                                     sat_ratio,
                                     temperature,
                                     e_sat,
                                     coeff_moldiff,
                                     coeff_curv,
                                     coeff_sol,
                                     solute_moles);
            }
#endif

            TI< dRsqdt<ParticleReal>,
                NewtonSolver<dRsqdt<ParticleReal>, ParticleReal>,
                ParticleReal > ti { drsqdt, newton_solver, a_dt, 100,
                                    sat_ratio, temperature, e_sat,
                                    coeff_moldiff, coeff_curv, coeff_sol,
                                    solute_moles,
                                    cfl, 1e-40, 1e-3, 1e-6, false, false };

            auto r_init = SD_effective_radius( i, ctx.idx_water,
                                               ctx.rho_water,
                                               ctx.num_species, ctx.num_aerosols,
                                               ptrs.sp_sol_arr, ptrs.ae_sol_arr,
                                               ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs,
                                               ptrs.sp_rho_arr, ptrs.ae_rho_arr );

            auto r_sq = r_init*r_init;
            bool success = false;
            if (ti_choice == SDMassChangeTIMethod::RK4) {
                ti.rk4(r_sq, success);
            } else if (ti_choice == SDMassChangeTIMethod::RK3BS) {
                ti.rk3bs(r_sq, success);
            } else if (ti_choice == SDMassChangeTIMethod::BE) {
                ti.be(r_sq, success);
            } else if (ti_choice == SDMassChangeTIMethod::CN) {
                ti.cn(r_sq, success);
            } else if (ti_choice == SDMassChangeTIMethod::DIRK2) {
                ti.dirk212(r_sq, success);
            }

            if (!success) {
                if (log_unconverged) {
#ifndef AMREX_USE_GPU
                    fprintf(file_handle,
                            "r=%1.16e, S=%1.16e, T=%1.16e, e=%1.16e, sol_mass=%1.16e\n",
                            ptrs.radius_ptr[i], sat_ratio, temperature, e_sat, solute_moles );
#endif
                }
                Gpu::Atomic::Add(unconverged_particles_ptr, Long(1));
            } else {
                // update particle attributes
                auto r_new = std::sqrt(r_sq);
                auto d_mass = four_thirds_pi*mat_density * (r_new*r_new*r_new - r_init*r_init*r_init);
                ptrs.sp_mass_ptrs[idx_vap][i] += d_mass;
                // don't let it go negative
                ptrs.sp_mass_ptrs[idx_vap][i] = std::max(ptrs.sp_mass_ptrs[idx_vap][i],amrex::Real(0));

                // Update particle attributes (radius and mass)
                SuperDropletPC::updateParticleAttributes(
                    i, ptrs.radius_ptr, ptrs.mass_ptr, ctx.idx_water, ctx.rho_water,
                    ctx.num_species, ctx.num_aerosols, ptrs.sp_sol_arr, ptrs.ae_sol_arr,
                    ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs, ptrs.sp_rho_arr, ptrs.ae_rho_arr);
            }

        });
        Gpu::synchronize();
        m_num_unconverged_particles += *(unconverged_particles.copyToHost());
    }); // end forEachParticleTile

}

#endif

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

namespace {
    /*! \brief Traits for liquid-vapour mass change process
     *  Requires ionization and molecular weight for Köhler equation */
    struct MassChangeLVTraits : SDProcess::DefaultTraits {
        static constexpr bool needs_ionization = true;
        static constexpr bool needs_mol_weight = true;
    };

    /*! \brief Traits for solid-liquid mass change (freezing/melting)
     *  Extends IceFullTraits with multiplicity for particle splitting */
    struct MassChangeSLTraits : SDProcess::IceFullTraits {
        static constexpr bool needs_multiplicity = true;
    };

    /*! \brief Traits for solid-vapour mass change (deposition/sublimation)
     *  Requires ice axes/rime (not Tfz) plus terminal velocity and multiplicity */
    struct MassChangeSVTraits : SDProcess::DefaultTraits {
        static constexpr bool needs_term_vel     = true;
        static constexpr bool needs_multiplicity = true;
        static constexpr bool needs_ice_axes     = true;
        static constexpr bool needs_ice_rime     = true;
    };

    /*! \brief Field indices for liquid-vapour interpolation */
    AMREX_ENUM(InterpFieldsLV, e_sat, sat_ratio, temperature, pressure, NUM_FIELDS);

    /*! \brief Field indices for solid-liquid interpolation */
    AMREX_ENUM(InterpFieldsSL, temperature, sat_ratio, NUM_FIELDS);

    /*! \brief Field indices for solid-vapour interpolation */
    AMREX_ENUM(InterpFieldsSV, e_sat, e_sat_ratio_wi, sat_ratio, density, temperature, pressure, moist_density, NUM_FIELDS);
}

/*! Compute mass change of particles due to evaporation and condensation
 *  (liquid <--> vapour) */
void SuperDropletPC::MassChange_LV (  int                                         a_lev,
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

    BL_PROFILE("SuperDropletPC::MassChange_LV()");
    AMREX_ASSERT( a_lev == m_lev );

    // Build process context using helper method
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

    forEachParticleTile<MassChangeLVTraits>(a_lev, ctx,
        [&](ParIterType& pti, int grid, ParticleType* p_pbox,
            const SDProcess::ParticlePointers& ptrs,
            const SDProcess::ProcessContext& ctx)
    {
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

        ParallelFor(ptrs.num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (ptrs.active_ptr[i] == 0) { return; }

            // skip ice particles
            auto par_phase = SD_phase(i, ctx.idx_water, ctx.idx_ice, ptrs.sp_mass_ptrs);
            if (a_is_water && (par_phase == SDPhase::ice)) { return; }

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

            ParticleReal solute_moles = 0.0;
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

            auto coeff_curv = vapour_mat.coeffCurv(temperature);
            auto coeff_sol = vapour_mat.coeffVPSolute();
            auto coeff_moldiff = vapour_mat.coeffMolecularDiffusion(temperature, pressure);

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
            } else {
                printf("ERROR: invalid time integrator choice!\n");
                return;
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
                auto d_mass = (4.0/3.0)*PI*mat_density * (r_new*r_new*r_new - r_init*r_init*r_init);
                ptrs.sp_mass_ptrs[idx_vap][i] += d_mass;
                // don't let it go negative
                ptrs.sp_mass_ptrs[idx_vap][i] = std::max(ptrs.sp_mass_ptrs[idx_vap][i],0.0);

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

/*! Compute mass change of particles due to freezing/melting
 *  (solid <--> liquid) */
void SuperDropletPC::MassChange_SL (  int                                         a_lev,
                                      Real                                        a_dt,
                                      const MultiFab&                             a_temperature,
                                      const MultiFab&                             a_sat_ratio,
                                      const Vector<std::unique_ptr<MultiFab>>&    a_z_phys_nd )
{
    BL_PROFILE("SuperDropletPC::MassChange_SL()");
    AMREX_ASSERT( a_lev == m_lev );
    amrex::ignore_unused(a_dt);

    if (m_idx_i < 0) { return; } // ice not being modeled

    // Build process context using helper method
    const auto ctx = buildProcessContext(a_lev);

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto is_periodic_z = geom.isPeriodic(2);

    const std::unique_ptr<MultiFab>& z_height = a_z_phys_nd[a_lev];

    AMREX_ALWAYS_ASSERT((ctx.idx_water >= 0) && (ctx.idx_ice >= 0) && (ctx.idx_water != ctx.idx_ice));

    forEachParticleTile<MassChangeSLTraits>(a_lev, ctx,
        [&](ParIterType& pti, int grid, ParticleType* p_pbox,
            const SDProcess::ParticlePointers& ptrs,
            const SDProcess::ProcessContext& ctx)
    {
        auto zheight = (*z_height)[grid].array();

        const auto& sat_ratio_arr = a_sat_ratio[grid].array();
        const auto& temperature_arr = a_temperature[grid].array();

        ParallelFor(ptrs.num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (ptrs.active_ptr[i] == 0) { return; }

            // Interpolate temperature, saturation ratio
            constexpr int nf = static_cast<int>(InterpFieldsSL::NUM_FIELDS);
            ParticleReal fv[nf];
            const Array4<const Real> fa[nf] = {
                temperature_arr, sat_ratio_arr
            };
            ERF::Interpolation::interpolateFields(
                p, ctx.plo, ctx.dxi, fa, fv, nf,
                is_periodic_z ? 1 : 0, is_periodic_z ? nullptr : &zheight
            );
            const auto temperature = fv[static_cast<int>(InterpFieldsSL::temperature)];
            const auto sat_ratio   = fv[static_cast<int>(InterpFieldsSL::sat_ratio)];

            auto par_phase = SD_phase(i, ctx.idx_water, ctx.idx_ice, ptrs.sp_mass_ptrs);
            if (par_phase == SDPhase::water) {
                // SD is water, check for freezing
                if ((temperature <= ptrs.Tfz_ptr[i]) && (sat_ratio > 1.0)) {
                    ptrs.sp_mass_ptrs[ctx.idx_ice][i] = ptrs.sp_mass_ptrs[ctx.idx_water][i];
                    ptrs.sp_mass_ptrs[ctx.idx_water][i] = 0.0;
                    ptrs.a_ptr[i] = ptrs.c_ptr[i] = std::cbrt(ptrs.sp_mass_ptrs[ctx.idx_ice][i]/(4.0/3.0*PI*ctx.rho_ice));
                    ptrs.mrime_ptr[i] = 0.0;
                    ptrs.nmono_ptr[i] = 1.0;
                }

            } else if (par_phase == SDPhase::ice) {

                // SD is ice, check for melting
                if (temperature > tmelt /* from ERF_Constants.H */) {
                    ptrs.sp_mass_ptrs[ctx.idx_water][i] = ptrs.sp_mass_ptrs[ctx.idx_ice][i];
                    ptrs.sp_mass_ptrs[ctx.idx_ice][i] = 0.0;
                    ptrs.a_ptr[i] = ptrs.c_ptr[i] = ptrs.mrime_ptr[i] = ptrs.nmono_ptr[i] = 0.0;
                }

            } else {
                amrex::Abort("Unknown value for particle phase (must be ice or water)");
            }

            // update particle attributes
            SuperDropletPC::updateParticleAttributes(
                i, ptrs.radius_ptr, ptrs.mass_ptr, ctx.idx_water, ctx.rho_water,
                ctx.num_species, ctx.num_aerosols, ptrs.sp_sol_arr, ptrs.ae_sol_arr,
                ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs, ptrs.sp_rho_arr, ptrs.ae_rho_arr);

        });
        Gpu::synchronize();
    }); // end forEachParticleTile

}

/*! Compute mass change of particles due to deposition and sublimation
 *  (solid <--> vapour) */
void SuperDropletPC::MassChange_SV (  int                                      a_lev,
                                      Real                                     a_dt,
                                      const MultiFab&                          a_density,
                                      const MultiFab&                          a_temperature,
                                      const MultiFab&                          a_pressure,
                                      const MultiFab&                          a_moist_density,
                                      const MultiFab&                          a_sat_pressure,
                                      const MultiFab&                          a_sat_pressure_wi,
                                      const MultiFab&                          a_sat_ratio,
                                      const Vector<std::unique_ptr<MultiFab>>& a_z_phys_nd )
{
    using namespace SDMassChangeUtils_SV;

    BL_PROFILE("SuperDropletPC::MassChange_SV()");
    AMREX_ASSERT( a_lev == m_lev );

    // Build process context using helper method
    const auto ctx = buildProcessContext(a_lev);

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto is_periodic_z = geom.isPeriodic(2);

    const std::unique_ptr<MultiFab>& z_height = a_z_phys_nd[a_lev];

    AMREX_ALWAYS_ASSERT((ctx.idx_water >= 0) && (ctx.idx_ice >= 0) && (ctx.idx_water != ctx.idx_ice));
    const auto mat_prop(*(m_species_mat[m_idx_i]));

    // dMdt functor (shared across tiles)
    dMdt<ParticleReal> dmdt{ mat_prop.m_lat_vap,
                             therco,
                             mat_prop.m_Rv,
                             mat_prop.m_density };

    forEachParticleTile<MassChangeSVTraits>(a_lev, ctx,
        [&](ParIterType& pti, int grid, ParticleType* p_pbox,
            const SDProcess::ParticlePointers& ptrs,
            const SDProcess::ProcessContext& ctx)
    {
        auto zheight = (*z_height)[grid].array();

        const auto& sat_pressure_arr = a_sat_pressure[grid].array();
        const auto& sat_pressure_wi_arr = a_sat_pressure_wi[grid].array();
        const auto& sat_ratio_arr = a_sat_ratio[grid].array();
        const auto& density_arr = a_density[grid].array();
        const auto& temperature_arr = a_temperature[grid].array();
        const auto& pressure_arr = a_pressure[grid].array();
        const auto& moist_density_arr = a_moist_density[grid].array();

        ParallelFor(ptrs.num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (ptrs.active_ptr[i] == 0) { return; }

            // skip water particles
            auto par_phase = SD_phase(i, ctx.idx_water, ctx.idx_ice, ptrs.sp_mass_ptrs);
            if (par_phase == SDPhase::water) { return; }

            // Interpolate fields for solid-vapour mass change
            constexpr int nf = static_cast<int>(InterpFieldsSV::NUM_FIELDS);
            ParticleReal fv[nf];
            const Array4<const Real> fa[nf] = {
                sat_pressure_arr, sat_pressure_wi_arr, sat_ratio_arr,
                density_arr, temperature_arr, pressure_arr, moist_density_arr
            };
            ERF::Interpolation::interpolateFields(
                p, ctx.plo, ctx.dxi, fa, fv, nf,
                is_periodic_z ? 1 : 0, is_periodic_z ? nullptr : &zheight
            );
            const auto e_sat          = fv[static_cast<int>(InterpFieldsSV::e_sat)];
            const auto e_sat_ratio_wi = fv[static_cast<int>(InterpFieldsSV::e_sat_ratio_wi)];
            const auto sat_ratio      = fv[static_cast<int>(InterpFieldsSV::sat_ratio)];
            // Note: density is interpolated but unused in this function
            const auto temperature    = fv[static_cast<int>(InterpFieldsSV::temperature)];
            const auto pressure       = fv[static_cast<int>(InterpFieldsSV::pressure)];
            const auto moist_density  = fv[static_cast<int>(InterpFieldsSV::moist_density)];

            auto coeff_moldiff = mat_prop.coeffMolecularDiffusion(temperature, pressure);

            TI< dMdt<ParticleReal>,
                ParticleReal > ti { dmdt, a_dt, a_dt, 100,
                                    ptrs.a_ptr[i], ptrs.c_ptr[i], ptrs.vterm_ptr[i],
                                    sat_ratio, moist_density, temperature, e_sat,
                                    coeff_moldiff,
                                    false };

            auto mass_old = ptrs.sp_mass_ptrs[ctx.idx_ice][i];
            auto vol_old = (4.0/3.0) * PI * ptrs.a_ptr[i]*ptrs.a_ptr[i]*ptrs.c_ptr[i];
            auto rhoi_old = mass_old / vol_old;

            // compute new mass
            bool success = false;
            auto mass_new = mass_old;
            ti.fe(mass_new, success);
            AMREX_ALWAYS_ASSERT(success);
            mass_new = std::max(mass_new, 3.8403e-24); // lower limit: 1nm spherical ice particle
            auto d_mass = mass_new - mass_old;

            // compute new volume
            auto d_vol = dmdt.dVolume(  d_mass, ptrs.a_ptr[i], ptrs.c_ptr[i], rhoi_old,
                                        sat_ratio, temperature,  e_sat, e_sat_ratio_wi,
                                        coeff_moldiff );
            auto vol_new = vol_old + d_vol;
            vol_new = std::max(vol_new, mass_new/mat_prop.m_density);
            d_vol = vol_new - vol_old;
            auto rhoi_new = mass_new / vol_new;

            // compute growth ratio and new radii
            auto gr_star = dmdt.growthRatioStar( d_mass,
                                                 ptrs.a_ptr[i], ptrs.c_ptr[i], ptrs.vterm_ptr[i],
                                                 moist_density, temperature, coeff_moldiff );
            auto d_loga = dmdt.dLogRadius(gr_star, std::log(vol_new/vol_old));
            auto a_new = ptrs.a_ptr[i] * std::exp(d_loga);
            auto c_new = ptrs.c_ptr[i] * std::exp(gr_star*d_loga);

            // Smaller than 1 um --> sphere with true ice density
            if (std::min(a_new,c_new) <= 1.0e-6) {
                rhoi_new = mat_prop.m_density;
                vol_new = mass_new / rhoi_new;
                c_new = a_new = std::cbrt(vol_new/(4.0*PI/3.0));
            }
            // limit particle size to 1 nm
            if (std::min(a_new,c_new) <= 1.0e-9 ) {
                c_new = a_new = 1.0e-9;
                rhoi_new = mat_prop.m_density;
                vol_new = (4.0/3.0)*PI*a_new*a_new*c_new;
                mass_new = vol_new * rhoi_new;
            }

            // rime mass
            auto mrime_new = ptrs.mrime_ptr[i];
            if (d_mass <= 0.0) { mrime_new = ptrs.mrime_ptr[i] * mass_new/mass_old; }

            // update attributes
            ptrs.a_ptr[i] = a_new;
            ptrs.c_ptr[i] = c_new;
            ptrs.mrime_ptr[i] = mrime_new;
            ptrs.sp_mass_ptrs[ctx.idx_ice][i] = mass_new;

            // update particle attributes
            SuperDropletPC::updateParticleAttributes(
                i, ptrs.radius_ptr, ptrs.mass_ptr, ctx.idx_water, ctx.rho_water,
                ctx.num_species, ctx.num_aerosols, ptrs.sp_sol_arr, ptrs.ae_sol_arr,
                ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs, ptrs.sp_rho_arr, ptrs.ae_rho_arr);

        });
        Gpu::synchronize();
    }); // end forEachParticleTile

}

#endif

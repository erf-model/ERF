#include "ERF_SuperDropletPC.H"
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuBuffer.H>
#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_Constants.H"
#include "ERF_SuperDropletPCMassChange.H"
#include "ERF_InterpolationUtils.H"
#include "ERF_MicrophysicsUtils.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;
using namespace SDPCDefn;

/*! \brief Field indices for liquid-vapour interpolation */
AMREX_ENUM(InterpFieldsLV, e_sat, sat_ratio, temperature, pressure, NUM_FIELDS);

/*! \brief Field indices for solid-liquid interpolation */
AMREX_ENUM(InterpFieldsSL, temperature, sat_ratio, moist_density, NUM_FIELDS);

/*! \brief Field indices for solid-vapour interpolation */
AMREX_ENUM(InterpFieldsSV, e_sat, e_sat_ratio_wi, sat_ratio, density, temperature, pressure, moist_density, NUM_FIELDS);

/*! Compute mass change of particles due to evaporation and condensation
 *  (liquid <--> vapour) */
void SuperDropletPC::MassChange_LV (  int                                         a_lev,
                                      double                                      a_dt,
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

    const auto proc_ctx = buildProcessContext(a_lev);

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

    // Ventilation enhancement applies to the water species only
    const bool use_ventilation = m_mass_change_ventilation && a_is_water;
    const auto vent_alpha1 = static_cast<ParticleReal>(m_vent_alpha1);
    const auto vent_beta1  = static_cast<ParticleReal>(m_vent_beta1);
    const auto vent_alpha2 = static_cast<ParticleReal>(m_vent_alpha2);
    const auto vent_beta2  = static_cast<ParticleReal>(m_vent_beta2);
    const auto vent_fcap   = static_cast<ParticleReal>(m_vent_fcap);

    // Solver setup (shared across tiles)
    dRsqdt<ParticleReal> drsqdt{ static_cast<ParticleReal>(vapour_mat.m_lat_vap),
                                 static_cast<ParticleReal>(therco),
                                 static_cast<ParticleReal>(vapour_mat.m_Rv),
                                 static_cast<ParticleReal>(mat_density) };

    NewtonSolver< dRsqdt<ParticleReal>, ParticleReal > newton_solver{ drsqdt,
                                                                      static_cast<ParticleReal>(m_newton_rtol),
                                                                      static_cast<ParticleReal>(m_newton_atol),
                                                                      static_cast<ParticleReal>(m_newton_stol),
                                                                      m_newton_maxits };

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    dRdt<ParticleReal> drdt{ static_cast<ParticleReal>(vapour_mat.m_lat_vap),
                             static_cast<ParticleReal>(therco),
                             static_cast<ParticleReal>(vapour_mat.m_Rv),
                             static_cast<ParticleReal>(mat_density) };
    constexpr int rtoff_r = SuperDropletsRealIdx::ncomps;
#endif

    forEachParticleTile(a_lev, proc_ctx,
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

            // skip particles with an ice core (pure ice or mixed); the melt step owns
            // wet-ice/vapour exchange
            auto par_phase = SD_phase(i, ctx.idx_water, ctx.idx_ice, ptrs.sp_mass_ptrs);
            if (a_is_water && (par_phase != SDPhase::water)) { return; }

            if (ERF::Interpolation::stencilOutOfBoundsZ(p, ctx.plo, ctx.dxi, zheight)) { return; }

            // Interpolate saturation pressure, saturation ratio, temperature, pressure
            constexpr int nf = static_cast<int>(InterpFieldsLV::NUM_FIELDS);
            ParticleReal fv[nf];
            const Array4<const Real> fa[nf] = {
                sat_pressure_arr, sat_ratio_arr, temperature_arr, pressure_arr
            };
            ERF::Interpolation::interpolateFields(
                p, ctx.plo, ctx.dxi, fa, fv, nf
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

            auto coeff_curv = vapour_mat_core.coeffCurv(static_cast<Real>(temperature));
            auto coeff_sol = vapour_mat_core.coeffVPSolute();
            auto coeff_moldiff = vapour_mat_core.coeffMolecularDiffusion(static_cast<Real>(temperature), static_cast<Real>(pressure));

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
            if (a_is_water) {
                ParticleReal f_v_diag = use_ventilation
                    ? ventilationFactor<ParticleReal>( ptrs.radius_ptr[i],
                                                       vent_alpha1, vent_beta1,
                                                       vent_alpha2, vent_beta2,
                                                       vent_fcap )
                    : ParticleReal(1);
                condt_ptr[i] = drdt( ptrs.radius_ptr[i],
                                     sat_ratio,
                                     temperature,
                                     e_sat,
                                     coeff_moldiff,
                                     coeff_curv,
                                     coeff_sol,
                                     solute_moles,
                                     f_v_diag);
            }
#endif

            TI< dRsqdt<ParticleReal>,
                NewtonSolver<dRsqdt<ParticleReal>, ParticleReal>,
                ParticleReal > ti { drsqdt, newton_solver, static_cast<ParticleReal>(a_dt), 100,
                                    sat_ratio, temperature, e_sat,
                                    static_cast<ParticleReal>(coeff_moldiff),
                                    static_cast<ParticleReal>(coeff_curv),
                                    static_cast<ParticleReal>(coeff_sol),
                                    solute_moles,
                                    static_cast<ParticleReal>(cfl),
                                    ParticleReal(1e-40), ParticleReal(1e-3), ParticleReal(1e-6),
                                    false, false,
                                    use_ventilation,
                                    vent_alpha1, vent_beta1,
                                    vent_alpha2, vent_beta2,
                                    vent_fcap };

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
                ptrs.sp_mass_ptrs[idx_vap][i] = std::max(ptrs.sp_mass_ptrs[idx_vap][i],ParticleReal(0));

                // Update particle attributes (radius and mass)
                SuperDropletPC::updateParticleAttributes(
                    i, ptrs.radius_ptr, ptrs.mass_ptr, ctx.idx_water, ctx.rho_water,
                    ctx.num_species, ctx.num_aerosols, ptrs.sp_sol_arr, ptrs.ae_sol_arr,
                    ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs, ptrs.sp_rho_arr, ptrs.ae_rho_arr);
            }

        });
        Gpu::synchronize();
#ifdef AMREX_USE_OMP
#pragma omp atomic
#endif
        m_num_unconverged_particles += static_cast<long>(*(unconverged_particles.copyToHost()));
    }); // end forEachParticleTile

}

/*! Compute mass change of particles due to freezing/melting
 *  (solid <--> liquid) */
void SuperDropletPC::MassChange_SL (  int                                         a_lev,
                                      double                                      a_dt,
                                      const MultiFab&                             a_temperature,
                                      const MultiFab&                             a_sat_ratio,
                                      const MultiFab&                             a_moist_density,
                                      const Vector<std::unique_ptr<MultiFab>>&    a_z_phys_nd )
{
    BL_PROFILE("SuperDropletPC::MassChange_SL()");
    AMREX_ASSERT(a_lev >= 0 && a_lev <= finestLevel());

    if (m_idx_i < 0) { return; } // ice not being modeled

    const auto proc_ctx = buildProcessContext(a_lev);

    const std::unique_ptr<MultiFab>& z_height = a_z_phys_nd[a_lev];

    AMREX_ALWAYS_ASSERT((proc_ctx.idx_water >= 0) && (proc_ctx.idx_ice >= 0) && (proc_ctx.idx_water != proc_ctx.idx_ice));

    // Water material and constants for the gradual melt rate (Seifert-Beheng 2006)
    const MaterialProperties water_mat(*(m_species_mat[m_idx_w]));
    const MaterialPropertiesCore& water_core = water_mat;
    const ParticleReal melt_L_f = static_cast<ParticleReal>(water_core.m_lat_fus);
    const ParticleReal melt_D   = static_cast<ParticleReal>(diffelq); // vapour diffusivity (water)
    const ParticleReal e_sat_w_T0 = static_cast<ParticleReal>(erf_esatw(static_cast<Real>(tmelt))*Real(100));
    SDMassChangeUtils_SV::dMdt<ParticleReal> dmelt{ static_cast<ParticleReal>(water_core.m_lat_vap),
                                                    static_cast<ParticleReal>(therco),
                                                    static_cast<ParticleReal>(water_core.m_Rv),
                                                    static_cast<ParticleReal>(water_core.m_density) };
    const ParticleReal ice_mass_min = static_cast<ParticleReal>(3.8403e-24); // ~1nm ice sphere

    forEachParticleTile(a_lev, proc_ctx,
        [&](ParIterType& /*pti*/, int grid, ParticleType* p_pbox,
            const SDProcess::ParticlePointers& ptrs,
            const SDProcess::ProcessContext& ctx)
    {
        auto zheight = (*z_height)[grid].array();

        const auto& sat_ratio_arr = a_sat_ratio[grid].array();
        const auto& temperature_arr = a_temperature[grid].array();
        const auto& moist_density_arr = a_moist_density[grid].array();

        ParallelFor(ptrs.num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (ptrs.active_ptr[i] == 0) { return; }

            if (ERF::Interpolation::stencilOutOfBoundsZ(p, ctx.plo, ctx.dxi, zheight)) { return; }

            // Interpolate temperature, saturation ratio, moist density
            constexpr int nf = static_cast<int>(InterpFieldsSL::NUM_FIELDS);
            ParticleReal fv[nf];
            const Array4<const Real> fa[nf] = {
                temperature_arr, sat_ratio_arr, moist_density_arr
            };
            ERF::Interpolation::interpolateFields(
                p, ctx.plo, ctx.dxi, fa, fv, nf
            );
            const auto temperature = fv[static_cast<int>(InterpFieldsSL::temperature)];
            const auto sat_ratio   = fv[static_cast<int>(InterpFieldsSL::sat_ratio)];
            const auto moist_density = fv[static_cast<int>(InterpFieldsSL::moist_density)];

            auto par_phase = SD_phase(i, ctx.idx_water, ctx.idx_ice, ptrs.sp_mass_ptrs);
            if (par_phase == SDPhase::water) {
                // SD is liquid water, check for freezing
                if ((temperature <= ptrs.Tfz_ptr[i]) && (sat_ratio > one)) {
                    auto m_water_i = ptrs.sp_mass_ptrs[ctx.idx_water][i];
                    if (m_water_i < ice_mass_min) {
                        // Too little water to make a resolvable ice crystal -- e.g. an
                        // essentially dry aerosol carrying only the trace initial water.
                        // Leave it as a dry aerosol instead of seeding a sub-floor ice
                        // particle (which would settle at a spurious IceBohm speed); it
                        // re-activates by condensation when it reaches supersaturation.
                        ptrs.sp_mass_ptrs[ctx.idx_water][i] = zero;
                    } else {
                        ptrs.sp_mass_ptrs[ctx.idx_ice][i] = m_water_i;
                        ptrs.sp_mass_ptrs[ctx.idx_water][i] = zero;
                        ptrs.a_ptr[i] = ptrs.c_ptr[i] = std::cbrt(m_water_i/(four_thirds_pi*ctx.rho_ice));
                        ptrs.mrime_ptr[i] = zero;
                        ptrs.nmono_ptr[i] = one;
                    }
                }

            } else {

                // SD has an ice core (pure ice or mixed wet ice)
                auto m_ice   = ptrs.sp_mass_ptrs[ctx.idx_ice][i];
                auto m_water = ptrs.sp_mass_ptrs[ctx.idx_water][i];

                if ((temperature > tmelt) && (m_ice > zero)) {

                    // Gradual melting (Seifert-Beheng 2006): meltwater accumulates on the
                    // particle, which becomes a mixed ice-water (wet) super-droplet.
                    auto e_sat_w_T = static_cast<ParticleReal>(erf_esatw(static_cast<Real>(temperature))*Real(100));
                    auto e_inf = sat_ratio * e_sat_w_T;
                    auto mdot = dmelt.meltRate( ptrs.a_ptr[i], ptrs.c_ptr[i], ptrs.vterm_ptr[i],
                                                temperature, e_inf, e_sat_w_T0,
                                                moist_density, melt_D, melt_L_f );

                    if (mdot < zero) { // melt only; sublimational cooling does not refreeze here
                        auto dm = std::max(mdot*static_cast<ParticleReal>(a_dt), -m_ice);
                        auto m_ice_new = m_ice + dm;
                        if (m_ice_new <= ice_mass_min) {
                            // fully melted: pure water droplet
                            ptrs.sp_mass_ptrs[ctx.idx_water][i] = m_water + m_ice;
                            ptrs.sp_mass_ptrs[ctx.idx_ice][i]   = zero;
                            ptrs.a_ptr[i] = ptrs.c_ptr[i] = ptrs.mrime_ptr[i] = zero;
                            ptrs.nmono_ptr[i] = zero;
                        } else {
                            // partial melt: shrink ice core, keep aspect ratio and apparent density
                            auto frac  = m_ice_new/m_ice;
                            auto scale = std::cbrt(frac);
                            ptrs.a_ptr[i] *= scale;
                            ptrs.c_ptr[i] *= scale;
                            ptrs.mrime_ptr[i] *= frac;
                            ptrs.sp_mass_ptrs[ctx.idx_ice][i]   = m_ice_new;
                            ptrs.sp_mass_ptrs[ctx.idx_water][i] = m_water - dm; // dm<0 -> water grows
                        }
                    }

                } else if ((temperature < tmelt) && (m_water > zero) && (m_ice > zero)) {

                    // Mixed particle in subfreezing air: refreeze the liquid film onto the core
                    auto m_ice_new = m_ice + m_water;
                    auto scale = std::cbrt(m_ice_new/m_ice);
                    ptrs.a_ptr[i] *= scale;
                    ptrs.c_ptr[i] *= scale;
                    ptrs.sp_mass_ptrs[ctx.idx_ice][i]   = m_ice_new;
                    ptrs.sp_mass_ptrs[ctx.idx_water][i] = zero;
                }
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
                                      double                                   a_dt,
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
    AMREX_ASSERT(a_lev >= 0 && a_lev <= finestLevel());

    const auto proc_ctx = buildProcessContext(a_lev);

    const std::unique_ptr<MultiFab>& z_height = a_z_phys_nd[a_lev];

    AMREX_ALWAYS_ASSERT((proc_ctx.idx_water >= 0) && (proc_ctx.idx_ice >= 0) && (proc_ctx.idx_water != proc_ctx.idx_ice));
    const MaterialProperties mat_prop(*(m_species_mat[m_idx_i]));
    const MaterialPropertiesCore& mat_prop_core = mat_prop;

    // dMdt functor (shared across tiles)
    dMdt<ParticleReal> dmdt{ static_cast<ParticleReal>(mat_prop_core.m_lat_vap),
                             static_cast<ParticleReal>(therco),
                             static_cast<ParticleReal>(mat_prop_core.m_Rv),
                             static_cast<ParticleReal>(mat_prop_core.m_density) };

    forEachParticleTile(a_lev, proc_ctx,
        [&](ParIterType& /*pti*/, int grid, ParticleType* p_pbox,
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

            // deposition/sublimation acts on a pure-ice surface only; skip water and
            // mixed wet-ice particles (the latter are melting, not depositing)
            auto par_phase = SD_phase(i, ctx.idx_water, ctx.idx_ice, ptrs.sp_mass_ptrs);
            if (par_phase != SDPhase::ice) { return; }

            if (ERF::Interpolation::stencilOutOfBoundsZ(p, ctx.plo, ctx.dxi, zheight)) { return; }

            // Interpolate fields for solid-vapour mass change
            constexpr int nf = static_cast<int>(InterpFieldsSV::NUM_FIELDS);
            ParticleReal fv[nf];
            const Array4<const Real> fa[nf] = {
                sat_pressure_arr, sat_pressure_wi_arr, sat_ratio_arr,
                density_arr, temperature_arr, pressure_arr, moist_density_arr
            };
            ERF::Interpolation::interpolateFields(
                p, ctx.plo, ctx.dxi, fa, fv, nf
            );
            const auto e_sat          = fv[static_cast<int>(InterpFieldsSV::e_sat)];
            const auto e_sat_ratio_wi = fv[static_cast<int>(InterpFieldsSV::e_sat_ratio_wi)];
            const auto sat_ratio      = fv[static_cast<int>(InterpFieldsSV::sat_ratio)];
            // Note: density is interpolated but unused in this function
            const auto temperature    = fv[static_cast<int>(InterpFieldsSV::temperature)];
            const auto pressure       = fv[static_cast<int>(InterpFieldsSV::pressure)];
            const auto moist_density  = fv[static_cast<int>(InterpFieldsSV::moist_density)];

            auto coeff_moldiff = mat_prop_core.coeffMolecularDiffusion(temperature, pressure);

            TI< dMdt<ParticleReal>,
                ParticleReal > ti { dmdt, ParticleReal(a_dt), ParticleReal(a_dt), 100,
                                    ptrs.a_ptr[i], ptrs.c_ptr[i], ptrs.vterm_ptr[i],
                                    sat_ratio, moist_density, temperature, e_sat,
                                    ParticleReal(coeff_moldiff),
                                    false };

            auto mass_old = ptrs.sp_mass_ptrs[ctx.idx_ice][i];
            auto vol_old = four_thirds_pi * ptrs.a_ptr[i]*ptrs.a_ptr[i]*ptrs.c_ptr[i];
            auto rhoi_old = mass_old / vol_old;

            // compute new mass
            bool success = false;
            auto mass_new = mass_old;
            ti.fe(mass_new, success);
            AMREX_ALWAYS_ASSERT(success);

            // Fully sublimated: release the ice back to a (dry) aerosol particle
            // instead of pinning a 1 nm "ice" sphere. A sub-resolution pinned
            // crystal would otherwise be given a spurious IceBohm fall speed (the
            // d = 2a characteristic length -> 0) and settle out unphysically; as a
            // water/aerosol particle its effective radius is set by the aerosol
            // mass, so it advects with the flow and can re-condense or re-freeze
            // when it reaches favourable conditions. Mirrors the full-melt path.
            const auto ice_mass_min = amrex::ParticleReal(3.8403e-24); // ~1 nm ice sphere
            if (mass_new <= ice_mass_min) {
                ptrs.sp_mass_ptrs[ctx.idx_ice][i] = zero;
                ptrs.a_ptr[i] = ptrs.c_ptr[i] = zero;
                ptrs.mrime_ptr[i] = zero;
                ptrs.nmono_ptr[i] = zero;
            } else {
                auto d_mass = mass_new - mass_old;

                // compute new volume
                auto d_vol = dmdt.dVolume(  d_mass, ptrs.a_ptr[i], ptrs.c_ptr[i], rhoi_old,
                                            sat_ratio, temperature,  e_sat, e_sat_ratio_wi,
                                            coeff_moldiff );
                auto vol_new = vol_old + d_vol;
                vol_new = std::max(vol_new, mass_new/mat_prop_core.m_density);
                auto rhoi_new = mass_new / vol_new;

                // compute growth ratio and new radii
                auto gr_star = dmdt.growthRatioStar( d_mass,
                                                     ptrs.a_ptr[i], ptrs.c_ptr[i], ptrs.vterm_ptr[i],
                                                     moist_density, temperature, coeff_moldiff );
                auto d_loga = dmdt.dLogRadius(gr_star, std::log(vol_new/vol_old));
                auto a_new = ptrs.a_ptr[i] * std::exp(d_loga);
                auto c_new = ptrs.c_ptr[i] * std::exp(gr_star*d_loga);

                // Smaller than 1 um --> sphere with true ice density. (Full
                // sublimation, where the size would otherwise fall below 1 nm, is
                // handled above by releasing the particle to an aerosol.)
                if (std::min(a_new,c_new) <= amrex::Real(1.0e-6)) {
                    rhoi_new = mat_prop_core.m_density;
                    vol_new = mass_new / rhoi_new;
                    c_new = a_new = std::cbrt(vol_new/four_thirds_pi);
                }

                // rime mass
                auto mrime_new = ptrs.mrime_ptr[i];
                if (d_mass <= zero) { mrime_new = ptrs.mrime_ptr[i] * mass_new/mass_old; }

                // update attributes
                ptrs.a_ptr[i] = a_new;
                ptrs.c_ptr[i] = c_new;
                ptrs.mrime_ptr[i] = mrime_new;
                ptrs.sp_mass_ptrs[ctx.idx_ice][i] = mass_new;
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

#endif

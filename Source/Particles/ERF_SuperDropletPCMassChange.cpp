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

namespace SDMassChangeUtils_SV {
    AMREX_GPU_CONSTANT amrex::Real s_tb_habit[121] = {
        amrex::Real(1.000000e+00), amrex::Real(9.799412e-01), amrex::Real(9.600636e-01), amrex::Real(9.399397e-01), amrex::Real(9.200258e-01),
        amrex::Real(9.001191e-01), amrex::Real(8.800351e-01), amrex::Real(8.570378e-01), amrex::Real(8.340652e-01), amrex::Real(8.109611e-01),
        amrex::Real(7.830689e-01), amrex::Real(7.550922e-01), amrex::Real(7.030723e-01), amrex::Real(5.370318e-01), amrex::Real(4.677351e-01),
        amrex::Real(5.248075e-01), amrex::Real(6.309573e-01), amrex::Real(8.128305e-01), amrex::Real(1.096478e+00), amrex::Real(1.479108e+00),
        amrex::Real(1.905461e+00), amrex::Real(2.089296e+00), amrex::Real(2.290868e+00), amrex::Real(2.398833e+00), amrex::Real(2.454709e+00),
        amrex::Real(2.426610e+00), amrex::Real(2.371374e+00), amrex::Real(2.290868e+00), amrex::Real(2.137962e+00), amrex::Real(1.995262e+00),
        amrex::Real(1.862087e+00), amrex::Real(1.737801e+00), amrex::Real(1.621810e+00), amrex::Real(1.513561e+00), amrex::Real(1.396368e+00),
        amrex::Real(1.288250e+00), amrex::Real(1.188502e+00), amrex::Real(1.096478e+00), amrex::Real(1.000000e+00), amrex::Real(9.225714e-01),
        amrex::Real(8.511380e-01), amrex::Real(7.852356e-01), amrex::Real(7.244360e-01), amrex::Real(6.683439e-01), amrex::Real(6.165950e-01),
        amrex::Real(5.754399e-01), amrex::Real(5.370318e-01), amrex::Real(5.011872e-01), amrex::Real(4.677351e-01), amrex::Real(4.365158e-01),
        amrex::Real(4.073803e-01), amrex::Real(3.801894e-01), amrex::Real(3.548134e-01), amrex::Real(3.311311e-01), amrex::Real(3.162278e-01),
        amrex::Real(3.019952e-01), amrex::Real(2.917427e-01), amrex::Real(2.851018e-01), amrex::Real(2.818383e-01), amrex::Real(2.786121e-01),
        amrex::Real(2.754229e-01), amrex::Real(2.786121e-01), amrex::Real(2.818383e-01), amrex::Real(2.851018e-01), amrex::Real(2.917427e-01),
        amrex::Real(2.985383e-01), amrex::Real(3.090295e-01), amrex::Real(3.198895e-01), amrex::Real(3.311311e-01), amrex::Real(3.467369e-01),
        amrex::Real(3.672823e-01), amrex::Real(3.935501e-01), amrex::Real(4.265795e-01), amrex::Real(4.570882e-01), amrex::Real(4.897788e-01),
        amrex::Real(5.248075e-01), amrex::Real(5.623413e-01), amrex::Real(6.095369e-01), amrex::Real(6.606934e-01), amrex::Real(7.161434e-01),
        amrex::Real(7.852356e-01), amrex::Real(8.609938e-01), amrex::Real(9.549926e-01), amrex::Real(1.047129e+00), amrex::Real(1.148154e+00),
        amrex::Real(1.258925e+00), amrex::Real(1.380384e+00), amrex::Real(1.496236e+00), amrex::Real(1.603245e+00), amrex::Real(1.698244e+00),
        amrex::Real(1.778279e+00), amrex::Real(1.840772e+00), amrex::Real(1.883649e+00), amrex::Real(1.905461e+00), amrex::Real(1.905461e+00),
        amrex::Real(1.883649e+00), amrex::Real(1.862087e+00), amrex::Real(1.840772e+00), amrex::Real(1.798871e+00), amrex::Real(1.737801e+00),
        amrex::Real(1.698244e+00), amrex::Real(1.640590e+00), amrex::Real(1.584893e+00), amrex::Real(1.549173e+00), amrex::Real(1.513910e+00),
        amrex::Real(1.476046e+00), amrex::Real(1.452112e+00), amrex::Real(1.428894e+00), amrex::Real(1.412863e+00), amrex::Real(1.393157e+00),
        amrex::Real(1.376892e+00), amrex::Real(1.361131e+00), amrex::Real(1.348963e+00), amrex::Real(1.336903e+00), amrex::Real(1.327089e+00),
        amrex::Real(1.317953e+00), amrex::Real(1.308881e+00), amrex::Real(1.302867e+00), amrex::Real(1.293898e+00), amrex::Real(1.287953e+00),
        amrex::Real(1.279087e+00)
    };
}

namespace {
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

    const auto ctx = buildProcessContext(a_lev);

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto is_periodic_z = geom.isPeriodic(2);

    const std::unique_ptr<MultiFab>& z_height = a_z_phys_nd[a_lev];

    AMREX_ALWAYS_ASSERT((ctx.idx_water >= 0) && (ctx.idx_ice >= 0) && (ctx.idx_water != ctx.idx_ice));

    forEachParticleTile(a_lev, ctx,
        [&](ParIterType& /*pti*/, int grid, ParticleType* p_pbox,
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
                if ((temperature <= ptrs.Tfz_ptr[i]) && (sat_ratio > one)) {
                    ptrs.sp_mass_ptrs[ctx.idx_ice][i] = ptrs.sp_mass_ptrs[ctx.idx_water][i];
                    ptrs.sp_mass_ptrs[ctx.idx_water][i] = zero;
                    ptrs.a_ptr[i] = ptrs.c_ptr[i] = std::cbrt(ptrs.sp_mass_ptrs[ctx.idx_ice][i]/(four_thirds_pi*ctx.rho_ice));
                    ptrs.mrime_ptr[i] = zero;
                    ptrs.nmono_ptr[i] = one;
                }

            } else if (par_phase == SDPhase::ice) {

                // SD is ice, check for melting
                if (temperature > tmelt /* from ERF_Constants.H */) {
                    ptrs.sp_mass_ptrs[ctx.idx_water][i] = ptrs.sp_mass_ptrs[ctx.idx_ice][i];
                    ptrs.sp_mass_ptrs[ctx.idx_ice][i] = zero;
                    ptrs.a_ptr[i] = ptrs.c_ptr[i] = ptrs.mrime_ptr[i] = ptrs.nmono_ptr[i] = zero;
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

    const auto ctx = buildProcessContext(a_lev);

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto is_periodic_z = geom.isPeriodic(2);

    const std::unique_ptr<MultiFab>& z_height = a_z_phys_nd[a_lev];

    AMREX_ALWAYS_ASSERT((ctx.idx_water >= 0) && (ctx.idx_ice >= 0) && (ctx.idx_water != ctx.idx_ice));
    const MaterialProperties mat_prop(*(m_species_mat[m_idx_i]));
    const MaterialPropertiesCore& mat_prop_core = mat_prop;

    // dMdt functor (shared across tiles)
    dMdt<ParticleReal> dmdt{ mat_prop_core.m_lat_vap,
                             therco,
                             mat_prop_core.m_Rv,
                             mat_prop_core.m_density };

    forEachParticleTile(a_lev, ctx,
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

            auto coeff_moldiff = mat_prop_core.coeffMolecularDiffusion(temperature, pressure);

            TI< dMdt<ParticleReal>,
                ParticleReal > ti { dmdt, a_dt, a_dt, 100,
                                    ptrs.a_ptr[i], ptrs.c_ptr[i], ptrs.vterm_ptr[i],
                                    sat_ratio, moist_density, temperature, e_sat,
                                    coeff_moldiff,
                                    false };

            auto mass_old = ptrs.sp_mass_ptrs[ctx.idx_ice][i];
            auto vol_old = four_thirds_pi * ptrs.a_ptr[i]*ptrs.a_ptr[i]*ptrs.c_ptr[i];
            auto rhoi_old = mass_old / vol_old;

            // compute new mass
            bool success = false;
            auto mass_new = mass_old;
            ti.fe(mass_new, success);
            AMREX_ALWAYS_ASSERT(success);
            mass_new = std::max(mass_new, amrex::Real(3.8403e-24)); // lower limit: 1nm spherical ice particle
            auto d_mass = mass_new - mass_old;

            // compute new volume
            auto d_vol = dmdt.dVolume(  d_mass, ptrs.a_ptr[i], ptrs.c_ptr[i], rhoi_old,
                                        sat_ratio, temperature,  e_sat, e_sat_ratio_wi,
                                        coeff_moldiff );
            auto vol_new = vol_old + d_vol;
            vol_new = std::max(vol_new, mass_new/mat_prop_core.m_density);
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
            if (std::min(a_new,c_new) <= amrex::Real(1.0e-6)) {
                rhoi_new = mat_prop_core.m_density;
                vol_new = mass_new / rhoi_new;
                c_new = a_new = std::cbrt(vol_new/four_thirds_pi);
            }
            // limit particle size to 1 nm
            if (std::min(a_new,c_new) <= amrex::Real(1.0e-9) ) {
                c_new = a_new = amrex::Real(1.0e-9);
                rhoi_new = mat_prop_core.m_density;
                vol_new = four_thirds_pi*a_new*a_new*c_new;
                mass_new = vol_new * rhoi_new;
            }

            // rime mass
            auto mrime_new = ptrs.mrime_ptr[i];
            if (d_mass <= zero) { mrime_new = ptrs.mrime_ptr[i] * mass_new/mass_old; }

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

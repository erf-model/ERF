#include <ERF.H>
#include <ERF_RadStruct.H>
#include <ERF_RadiationDiagnostics.H>
#include <ERF_TwoStreamColumn.H>
#include <ERF_PrognosticCloudFraction.H>
#include <ERF_AerosolOpticalDepth.H>
#include <ERF_SimplifiedSEB.H>
#include <ERF_SolarGeometry.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Gpu.H>
#include <ERF_IndexDefines.H>
#include <ERF_EOS.H>
#include <cmath>
#include <limits>

using namespace amrex;


namespace {
void fill_or_copy_seb_field(
    MultiFab* seb_mf,
    LandSurface& lsm,
    int lev,
    const char* field_name,
    amrex::Real fallback_value)
{
    if (seb_mf == nullptr) return;

    std::string varname(field_name);
    int lsm_idx = lsm.Get_DataIdx(lev, varname);
    if (lsm_idx >= 0) {
        if (MultiFab* lsm_ptr = lsm.Get_Data_Ptr(lev, lsm_idx)) {
            MultiFab::Copy(*seb_mf, *lsm_ptr, 0, 0, 1, 0);
            return;
        }
    }
    seb_mf->setVal(fallback_value);
}
}

/**
 * @file ERF_AdvanceTwoStreamRadiation.cpp
 * @brief Two-stream radiation driver with per-level heating rates.
 *
 * Computes SW/LW fluxes and per-level heating rates using real per-column
 * vertical sweeps over the atmospheric grid. Reads temperature and density
 * from the state and properly accumulates optical depth through all
 * vertical levels.
 *
 * Capabilities include:
 * - Height-varying optical depth: optional "cloud_layer" tau_profile_type
 *   adds cloud_tau_per_layer on top of the clear-sky background within
 *   [cloud_base_height_m, cloud_top_height_m].
 * - Cloud fraction masking: blends clear-sky and cloudy-column fluxes via
 *   F = (1 - cloud_fraction) * F_clear + cloud_fraction) * F_cloudy.
 * - Diffuse (scattered) SW flux: during the downward SW sweep, each layer's
 *   direct-beam attenuation also generates a diffuse flux contribution
 *   via compute_sw_diffuse_flux() (Meador-Weaver two-stream scattering),
 *   accumulated layer-by-layer alongside the direct beam.
 * - Per-level heating rate output: writes SW/LW heating rates to a
 *   2-component heating-rate MultiFab (component 0 = SW, component 1 = LW),
 *   mirroring the RRTMGP convention.
 *
 * Includes support for:
 * - Height-varying surface properties (albedo, emissivity, temperature)
 * - Dynamic cloud fraction from relative humidity and cloud water
 * - Prescribed bulk aerosol optical depth
 * - Dynamic solar geometry (time-varying solar position)
 * - Surface energy balance diagnostics and prognostic updates
 *
 * Vertical orientation follows ERF: k = kmin is the surface layer and
 * k = kmax the top layer. SW sweeps downward from kmax to kmin; LW sweeps
 * downward (TOA -> surface) and then upward (surface -> TOA). Layer
 * temperature is obtained from rho*theta through the Exner function.
 *
 * Note: CSV diagnostics (SW_surface, heating_rate_max, etc.) are maintained
 * for backward RegTest compatibility. heating_rate_max is the max(|Q_sw|+|Q_lw|)
 * observed during the column evaluation.
 */

void ERF::compute_twostream_radiation_diagnostics(
    int lev,
    int nstep,
    amrex::Real time_step,
    std::string const& call_site
    )
{
    const auto& rad_choice = solverChoice.radChoice;

    // Only proceed if TwoStream radiation is enabled
    if (rad_choice.rad_type != RadType::TwoStream) {
        return;
    }

    // Create RadiationDiagnostics instance for this level with controls
    RadiationDiagnostics rad_diag(rad_choice.verbosity, rad_choice.diag_file, 
                                   rad_choice.diag_enable, rad_choice.diag_stdout_enable,
                                   rad_choice.diag_tagged_enable, rad_choice.diag_regtest_line_enable,
                                   rad_choice.diag_csv_enable, rad_choice.diag_callsite_mode,
                                   rad_choice.diag_dedup_tol);

    // ========================================
    // GPU-Safe ParallelFor Implementation with Cloud Fraction
    // Blending, Diffuse (Scattering) SW Flux, and Per-Level Heating Rate
    // Output (qheating_rates MultiFab)
    // ========================================

    // Initialize global diagnostics
    amrex::Real SW_surface = 0.0;
    amrex::Real SW_TOA = 0.0;
    amrex::Real F_up_surface = 0.0;
    amrex::Real F_down_toa = 0.0;
    amrex::Real heating_rate_max = 0.0;
    amrex::Real seb_residual_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real seb_residual_max  = std::numeric_limits<amrex::Real>::quiet_NaN();
    
    //  Prognostic SEB surface temperature and moisture diagnostics
    amrex::Real t_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real t_s_max  = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real q_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real q_s_max  = std::numeric_limits<amrex::Real>::quiet_NaN();

    // Get state at this level (conservative variables: density, RhoTheta, etc.)
    const auto& state_cons = vars_old[lev][Vars::cons];

    // Only compute radiation if we have valid state data
    if (state_cons.nComp() > 0 ) {

    // Prepare to compute TOA values (used for diagnostics output)
    // Compute dynamic cos_zenith if enabled, otherwise use static value
    amrex::Real cos_zenith;
    if (rad_choice.solar_geometry_dynamic_enable) {
        // Convert absolute simulation time to UTC seconds within the day [0, 86400)
        amrex::Real time_utc_seconds = std::fmod(t_old[lev], 86400.0);
        if (time_utc_seconds < 0.0) time_utc_seconds += 86400.0;
        cos_zenith = compute_cos_zenith_angle(
            time_utc_seconds,
            rad_choice.latitude_deg,
            rad_choice.longitude_deg,
            rad_choice.day_of_year,
            rad_choice.time_zone_offset_hours);
    } else {
        // and earlier: Use static solar zenith angle
        amrex::Real zenith_rad = rad_choice.solar_zenith_deg * M_PI / 180.0;
        cos_zenith = std::cos(zenith_rad);
    }
    SW_TOA = rad_choice.sw_enabled ? (rad_choice.S0 * std::max(0.0, cos_zenith)) : 0.0;

        // Host-side storage for reduction results (will be set by device-side reduction)
        amrex::Real max_heating_global = 0.0;
        amrex::Real sw_surface_sum = 0.0;
        amrex::Real lw_net_sum = 0.0;
        amrex::Long n_columns_total = 0;

        // SEB residual diagnostics
        amrex::Real seb_residual_sum = 0.0;
        amrex::Long n_seb_columns = 0;

        // cloud fraction used to blend clear-sky and cloudy-column results.
        // cloud_fraction == 0.0 (default) means only the clear-sky column is
        // ever evaluated, and the blend below reduces to F = F_clear exactly.
        amrex::Real cloud_fraction = rad_choice.cloud_fraction;

        // Compute UTC seconds within the day for dynamic solar geometry
        amrex::Real time_utc_seconds = 0.0;
        if (rad_choice.solar_geometry_dynamic_enable) {
            time_utc_seconds = std::fmod(t_old[lev], 86400.0);
            if (time_utc_seconds < 0.0) time_utc_seconds += 86400.0;
        }

        // Note: qheating_rates[lev] is expected to be allocated with
        // 2 components by the caller whenever rad_choice.rad_type ==
        // RadType::TwoStream (see Source/ERF_MakeNewArrays.cpp). If not yet
        // allocated (e.g. before initialization), this function still safely
        // computes and logs CSV diagnostics but skips the per-level heating write.
        MultiFab* qheating_mf = qheating_rates[lev].get();

        if (rad_choice.seb_enable) {
            fill_or_copy_seb_field(twostream_alb_sw[lev].get(), lsm, lev, "sfc_alb_dir_vis", rad_choice.surface_albedo_sw);
            fill_or_copy_seb_field(twostream_emiss_lw[lev].get(), lsm, lev, "sfc_emis", rad_choice.surface_emissivity_lw);
             
            // Gate t_sfc fill on prognostic mode: when seb_prognostic_enable is true,
            // t_sfc is owned and evolved by the prognostic update, not reset by fill_or_copy.
            // This prevents silently overwriting the prognostic state before the update reads it.
            //fill_or_copy_seb_field(twostream_t_sfc[lev].get(), lsm, lev, "t_sfc", rad_choice.surface_temp_k);
            if (!rad_choice.seb_prognostic_enable) {
                fill_or_copy_seb_field(twostream_t_sfc[lev].get(), lsm, lev, "t_sfc", rad_choice.surface_temp_k);
            }
             
            fill_or_copy_seb_field(sw_flux_sfc[lev].get(), lsm, lev, "sav", rad_choice.seb_sw_flux_default);
            fill_or_copy_seb_field(lw_flux_sfc[lev].get(), lsm, lev, "fira", rad_choice.seb_lw_flux_default);
            fill_or_copy_seb_field(hfx_sfc[lev].get(), lsm, lev, "hfx", rad_choice.seb_hfx_default);
            fill_or_copy_seb_field(lh_sfc[lev].get(), lsm, lev, "lh", rad_choice.seb_lh_default);
            fill_or_copy_seb_field(grdflx_sfc[lev].get(), lsm, lev, "grdflx", rad_choice.seb_grdflx_default);
             
            // Gate q_sfc fill on prognostic mode: same reasoning as t_sfc.
            if (!rad_choice.seb_prognostic_enable) {
                fill_or_copy_seb_field(q_sfc[lev].get(), lsm, lev, "noahmp_water_vapor_mixing_ratio_2m_vegetated", rad_choice.seb_q_sfc_default);
            }
            fill_or_copy_seb_field(t_deep[lev].get(), lsm, lev, "smstav", rad_choice.seb_t_deep_default);
            fill_or_copy_seb_field(q_deep[lev].get(), lsm, lev, "smstot", rad_choice.seb_q_deep_default);
        }

        // Sequential loop over all boxes (each box handled with GPU-safe ParallelFor)
        for (MFIter mfi(state_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto& state_arr = state_cons.const_array(mfi);
            const Geometry& geom_lev = geom[lev];

            // Get z_phys_cc for nonuniform dz support if available
            Array4<const amrex::Real> z_phys_cc_arr;
            //bool has_z_phys = false;
            if (z_phys_cc[lev] != nullptr) {
                z_phys_cc_arr = z_phys_cc[lev]->const_array(mfi);
                //has_z_phys = true;
            }

            //  Wire LSM surface property fields or use standalone fallback MultiFabs
            // Priority:
            // 1. If LSM is active (lsm.Get_DataIdx() returns >=0), use real LSM fields
            // 2. Otherwise, use standalone fallback MultiFabs (allocated and constant-filled from RadChoice scalars)
            // The resolve_surface_*() helpers implement the full precedence chain with finite guards
            
            // SW albedo: Try LSM field "sfc_alb_dir_vis" (simplified: broadband approx from vis-direct only;
            // future work: full 4-band vis/nir dir/dif support is planned).
            bool has_hetero_alb_sw = false;
            Array4<const amrex::Real> hetero_alb_sw_arr;
            {
                //int lsm_idx = lsm.Get_DataIdx(lev, "sfc_alb_dir_vis");
                std::string varname_alb = "sfc_alb_dir_vis";
                int lsm_idx = lsm.Get_DataIdx(lev, varname_alb);
                if (lsm_idx >= 0) {
                    auto lsm_ptr = lsm.Get_Data_Ptr(lev, lsm_idx);
                    if (lsm_ptr) {
                        hetero_alb_sw_arr = lsm_ptr->const_array(mfi);
                        has_hetero_alb_sw = true;
                    }
                } else if (twostream_alb_sw[lev]) {
                    hetero_alb_sw_arr = twostream_alb_sw[lev]->const_array(mfi);
                    has_hetero_alb_sw = true;
                }
            }

            // LW emissivity: Try LSM field "sfc_emis"
            bool has_hetero_emiss_lw = false;
            Array4<const amrex::Real> hetero_emiss_lw_arr;
            {
                //int lsm_idx = lsm.Get_DataIdx(lev, "sfc_emis");
                std::string varname_emiss = "sfc_emis";
                int lsm_idx = lsm.Get_DataIdx(lev, varname_emiss);
                if (lsm_idx >= 0) {
                    auto lsm_ptr = lsm.Get_Data_Ptr(lev, lsm_idx);
                    if (lsm_ptr) {
                        hetero_emiss_lw_arr = lsm_ptr->const_array(mfi);
                        has_hetero_emiss_lw = true;
                    }
                } else if (twostream_emiss_lw[lev]) {
                    hetero_emiss_lw_arr = twostream_emiss_lw[lev]->const_array(mfi);
                    has_hetero_emiss_lw = true;
                }
            }

            // Surface temperature: Try LSM field "t_sfc"
            bool has_t_sfc_field = false;
            Array4<const amrex::Real> t_sfc_arr;
            {
                //int lsm_idx = lsm.Get_DataIdx(lev, "t_sfc");
                std::string varname_t_sfc = "t_sfc";
                int lsm_idx = lsm.Get_DataIdx(lev, varname_t_sfc);
                if (lsm_idx >= 0) {
                    auto lsm_ptr = lsm.Get_Data_Ptr(lev, lsm_idx);
                    if (lsm_ptr) {
                        t_sfc_arr = lsm_ptr->const_array(mfi);
                        has_t_sfc_field = true;
                    }
                } else if (twostream_t_sfc[lev]) {
                    t_sfc_arr = twostream_t_sfc[lev]->const_array(mfi);
                    has_t_sfc_field = true;
                }
            }


            // Create a 2D box for (i,j) iteration over the horizontal extent
            // One GPU thread per (i,j) column; k-loop is sequential within each thread
            const auto& lo = bx.loVect();
            const auto& hi = bx.hiVect();
            Box xy_box(IntVect(lo[0], lo[1], 0), IntVect(hi[0], hi[1], 0));

            // Count columns in this box for later averaging
            amrex::Long n_cols = static_cast<amrex::Long>(bx.length(0)) *
                                 static_cast<amrex::Long>(bx.length(1));
            n_columns_total += n_cols;

            // Clear-sky-column heating rates are written directly
            // into the real qheating_rates MultiFab when available.
            // Fall back to a throwaway local FArrayBox otherwise (keeps the
            // kernel call GPU-safe even if qheating_rates isn't allocated).
            FArrayBox qheating_fallback_fab;
            Array4<amrex::Real> qheating_clear_arr;
            if (qheating_mf != nullptr) {
                qheating_clear_arr = qheating_mf->array(mfi);
            } else {
                qheating_fallback_fab.resize(bx, 2);
                qheating_clear_arr = qheating_fallback_fab.array();
            }

            // Cloudy-column heating rates always go into a scratch
            // FArrayBox; only used/blended when cloud_fraction > 0.
            FArrayBox qheating_cloudy_fab;
            Array4<amrex::Real> qheating_cloudy_arr;
            if (cloud_fraction > 0.0) {
                qheating_cloudy_fab.resize(bx, 2);
                qheating_cloudy_arr = qheating_cloudy_fab.array();
            }

            // GPU-safe reduction using ReduceOps (per-column results aggregated on device)
            amrex::Real max_heating_box = 0.0;
            amrex::Real sw_sum_box = 0.0;
            amrex::Real lw_sum_box = 0.0;

            // Device-side reduction: compute max heating and sum of surface fluxes
            ReduceOps<ReduceOpMax, ReduceOpSum, ReduceOpSum> reduce_ops;
            ReduceData<amrex::Real, amrex::Real, amrex::Real> reduce_data(reduce_ops);

            using ReduceTuple = typename decltype(reduce_data)::Type;

         // Launch parallel kernel over (i,j) columns
            reduce_ops.eval(xy_box, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int /*k_unused*/) -> ReduceTuple
                {
                    // Clear-sky column (always evaluated; this is the sole
                    // contributor when cloud_fraction == 0.0, matching earlier behavior)
                    amrex::Real max_heating_clear = 0.0;
                    amrex::Real sw_flux_clear = 0.0;
                    amrex::Real lw_net_clear = 0.0;
                    vertical_two_stream_sweep(
                        i, j, bx, geom_lev, state_arr, rad_choice, /*cloudy=*/false,
                        qheating_clear_arr,
                        max_heating_clear, sw_flux_clear, lw_net_clear,
                        z_phys_cc_arr,
                        time_utc_seconds,
                        has_hetero_alb_sw, &hetero_alb_sw_arr,
                        has_hetero_emiss_lw, &hetero_emiss_lw_arr,
                        has_t_sfc_field, &t_sfc_arr);

                    amrex::Real max_heating_col = max_heating_clear;
                    amrex::Real sw_flux_col = sw_flux_clear;
                    amrex::Real lw_net_col = lw_net_clear;

                    // Cloudy column only needs to be evaluated when there is a
                    // nonzero cloud fraction; this keeps the cloud_fraction==0
                    // path numerically and computationally identical.
                    if (cloud_fraction > 0.0) {
                         amrex::Real max_heating_cloudy = 0.0;
                         amrex::Real sw_flux_cloudy = 0.0;
                         amrex::Real lw_net_cloudy = 0.0;
                         vertical_two_stream_sweep(
                            i, j, bx, geom_lev, state_arr, rad_choice, /*cloudy=*/true,
                            qheating_cloudy_arr,
                            max_heating_cloudy, sw_flux_cloudy, lw_net_cloudy,
                             z_phys_cc_arr,
                             time_utc_seconds,
                             has_hetero_alb_sw, &hetero_alb_sw_arr,
                             has_hetero_emiss_lw, &hetero_emiss_lw_arr,
                             has_t_sfc_field, &t_sfc_arr);

                        // Blend clear-sky and cloudy-column results
                        sw_flux_col = (1.0 - cloud_fraction) * sw_flux_clear +
                                      cloud_fraction * sw_flux_cloudy;
                        lw_net_col = (1.0 - cloud_fraction) * lw_net_clear +
                                     cloud_fraction * lw_net_cloudy;
                        max_heating_col = std::max(max_heating_clear, max_heating_cloudy);

                        // Blend per-level heating rates in place
                        // into qheating_clear_arr (which is the real output
                        // MultiFab when qheating_mf != nullptr).
                        int kmin = bx.smallEnd(2);
                        int kmax = bx.bigEnd(2);
                        for (int k = kmin; k <= kmax; ++k) {
                            for (int comp = 0; comp < 2; ++comp) {
                                amrex::Real q_clear_val = qheating_clear_arr(i, j, k, comp);
                                amrex::Real q_cloudy_val = qheating_cloudy_arr(i, j, k, comp);
                                qheating_clear_arr(i, j, k, comp) =
                                    (1.0 - cloud_fraction) * q_clear_val +
                                    cloud_fraction * q_cloudy_val;
                            }
                        }
                    }

                    // Return tuple for reduction
                    return {max_heating_col, sw_flux_col, lw_net_col};
                }
            );

            // Copy results from device to host
            amrex::Gpu::synchronize();
            auto reduce_tuple = reduce_data.value(reduce_ops);
            max_heating_box = amrex::get<0>(reduce_tuple);
            sw_sum_box = amrex::get<1>(reduce_tuple);
            lw_sum_box = amrex::get<2>(reduce_tuple);

            // Accumulate box results into global results
            max_heating_global = std::max(max_heating_global, max_heating_box);
            sw_surface_sum += sw_sum_box;
            lw_net_sum += lw_sum_box;
        }

         // Warn if diagnostic is requested but SEB infrastructure isn't enabled
        if (rad_choice.seb_diagnostic_enable && !rad_choice.seb_enable) {
            static bool warned_seb_misconfig = false;
            if (!warned_seb_misconfig && ParallelDescriptor::IOProcessor()) {
                Print() << "WARNING: erf.radiation.seb_diagnostic_enable=true but "
                           "seb_enable=false; SEB residual diagnostics will report NaN. "
                           "Set erf.radiation.seb_enable=true to enable SEB field "
                           "population.\n";
                warned_seb_misconfig = true;
            }
        }
        // Compute SEB residual diagnostics if enabled
        if (rad_choice.seb_diagnostic_enable && rad_choice.seb_enable) {
            seb_residual_max = 0.0; 
            // Second loop over boxes to compute SEB residual from populated SEB MultiFabs
            for (MFIter mfi(state_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const auto& lo = bx.loVect();
                const auto& hi = bx.hiVect();
                Box xy_box(IntVect(lo[0], lo[1], 0), IntVect(hi[0], hi[1], 0));

                // Get SEB field arrays
                Array4<const amrex::Real> sw_flux_arr = sw_flux_sfc[lev]->const_array(mfi);
                Array4<const amrex::Real> lw_flux_arr = lw_flux_sfc[lev]->const_array(mfi);
                Array4<const amrex::Real> hfx_arr = hfx_sfc[lev]->const_array(mfi);
                Array4<const amrex::Real> lh_arr = lh_sfc[lev]->const_array(mfi);
                Array4<const amrex::Real> grdflx_arr = grdflx_sfc[lev]->const_array(mfi);

                // Count columns and compute residuals
                amrex::Long n_cols_box = static_cast<amrex::Long>(bx.length(0)) *
                                        static_cast<amrex::Long>(bx.length(1));
                amrex::Real residual_sum_box = 0.0;
                amrex::Real residual_max_box = 0.0;

                // GPU-safe reduction for SEB residuals
                ReduceOps<ReduceOpSum, ReduceOpMax> seb_reduce_ops;
                ReduceData<amrex::Real, amrex::Real> seb_reduce_data(seb_reduce_ops);

                using SEBReduceTuple = typename decltype(seb_reduce_data)::Type;

                seb_reduce_ops.eval(xy_box, seb_reduce_data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int /*k_unused*/) -> SEBReduceTuple {
                        amrex::Real sw_net = sw_flux_arr(i, j, 0);
                        amrex::Real lw_net = lw_flux_arr(i, j, 0);
                        amrex::Real hfx = hfx_arr(i, j, 0);
                        amrex::Real lh = lh_arr(i, j, 0);
                        amrex::Real grdflx = grdflx_arr(i, j, 0);

                        // Compute residual using helper function
                        amrex::Real residual = diagnose_seb_residual(sw_net, lw_net, hfx, lh, grdflx);

                        // Return sum and abs(max) of residual
                        return {residual, std::abs(residual)};
                    }
                );

                // Copy results from device to host
                amrex::Gpu::synchronize();
                auto seb_reduce_tuple = seb_reduce_data.value(seb_reduce_ops);
                residual_sum_box = amrex::get<0>(seb_reduce_tuple);
                residual_max_box = amrex::get<1>(seb_reduce_tuple);

                // Accumulate into global results
                seb_residual_sum += residual_sum_box;
                seb_residual_max = std::max(seb_residual_max, residual_max_box);
                n_seb_columns += n_cols_box;
            }
        }

        //  Prognostic SEB surface temperature and moisture evolution
        // Only run if prognostic mode is enabled and Noah-MP is NOT driving LSM at this level
        if (rad_choice.seb_prognostic_enable && rad_choice.seb_enable &&
            call_site == "post_dycore") {
            // Check if Noah-MP is active at this level by attempting to get the LSM t_sfc field
            std::string varname_t_sfc_prog = "t_sfc";
            int lsm_idx_t_sfc = lsm.Get_DataIdx(lev, varname_t_sfc_prog);
            bool noahmp_active = (lsm_idx_t_sfc >= 0);
            
            if (!noahmp_active) {
                // Noah-MP is NOT active; proceed with prognostic update
                
                // Initialize diagnostics for T_s and q_s
                amrex::Real t_s_sum = 0.0;
                amrex::Real t_s_max_val = -std::numeric_limits<amrex::Real>::max();
                amrex::Real q_s_sum = 0.0;
                amrex::Real q_s_max_val = -std::numeric_limits<amrex::Real>::max();
                amrex::Long n_prog_columns = 0;
                
                // Third loop over boxes for prognostic SEB update
                for (MFIter mfi(state_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.tilebox();
                    const auto& lo = bx.loVect();
                    const auto& hi = bx.hiVect();
                    Box xy_box(IntVect(lo[0], lo[1], 0), IntVect(hi[0], hi[1], 0));
                    
                    // Get SEB field arrays (read-only)
                    Array4<const amrex::Real> sw_flux_arr = sw_flux_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> lw_flux_arr = lw_flux_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> hfx_arr = hfx_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> lh_arr = lh_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> grdflx_arr = grdflx_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> t_deep_arr = t_deep[lev]->const_array(mfi);
                    Array4<const amrex::Real> q_deep_arr = q_deep[lev]->const_array(mfi);
                    
                    // Get SEB state arrays (read-write for prognostic update)
                    Array4<amrex::Real> t_s_arr = twostream_t_sfc[lev]->array(mfi);
                    Array4<amrex::Real> q_s_arr = q_sfc[lev]->array(mfi);
                    // Count columns and prepare for reductions
                    amrex::Long n_cols_box = static_cast<amrex::Long>(bx.length(0)) *
                                             static_cast<amrex::Long>(bx.length(1));
                    amrex::Real t_s_sum_box = 0.0;
                    amrex::Real t_s_max_box = -std::numeric_limits<amrex::Real>::max();
                    amrex::Real q_s_sum_box = 0.0;
                    amrex::Real q_s_max_box = -std::numeric_limits<amrex::Real>::max();
                    
                    // GPU-safe update for prognostic T_s and q_s with reductions
                    ReduceOps<ReduceOpSum, ReduceOpMax, ReduceOpSum, ReduceOpMax> prog_reduce_ops;
                    ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real> prog_reduce_data(prog_reduce_ops);
                    
                    using ProgReduceTuple = typename decltype(prog_reduce_data)::Type;
                    
                    prog_reduce_ops.eval(xy_box, prog_reduce_data,
                            [=,  C_s=rad_choice.seb_surface_heat_capacity,
                            tau=rad_choice.seb_restore_timescale_s,
                            d_s=rad_choice.seb_moisture_layer_depth_m,
                            tau_q=rad_choice.seb_moisture_restore_timescale_s,
                            t_min=rad_choice.seb_prognostic_t_min_k,
                            t_max=rad_choice.seb_prognostic_t_max_k,
                            q_min=rad_choice.seb_prognostic_q_min,
                            q_max=rad_choice.seb_prognostic_q_max] 
                            AMREX_GPU_DEVICE (int i, int j, int /*k_unused*/) -> ProgReduceTuple {
                            amrex::Real t_s_old = t_s_arr(i, j, 0);
                            amrex::Real q_s_old = q_s_arr(i, j, 0);
                            
                            // Read forcing data
                            amrex::Real sw_net = sw_flux_arr(i, j, 0);
                            amrex::Real lw_net = lw_flux_arr(i, j, 0);
                            amrex::Real hfx = hfx_arr(i, j, 0);
                            amrex::Real lh = lh_arr(i, j, 0);
                            amrex::Real grdflx = grdflx_arr(i, j, 0);
                            amrex::Real t_deep = t_deep_arr(i, j, 0);
                            amrex::Real q_deep = q_deep_arr(i, j, 0);
                            
                            // Compute SEB residual
                            amrex::Real seb_res = diagnose_seb_residual(sw_net, lw_net, hfx, lh, grdflx);
                            
                            // Compute tendencies
                            amrex::Real dT_s_dt = prognostic_dTs_dt(seb_res, t_s_old, t_deep,
                                                                     C_s, tau);
                            amrex::Real dq_s_dt = prognostic_dqs_dt(lh, q_s_old, q_deep,
                                                                     d_s, tau_q);
                            
                            // Perform Euler update
                            amrex::Real t_s_new = t_s_old + time_step * dT_s_dt;
                            amrex::Real q_s_new = q_s_old + time_step * dq_s_dt;
                            
                            // Clamp to valid ranges
                            t_s_new = amrex::max(t_min, amrex::min(t_max, t_s_new));
                            q_s_new = amrex::max(q_min, amrex::min(q_max, q_s_new));
                            
                            // Write back updated values (this modifies the device array)
                            t_s_arr(i, j, 0) = t_s_new;
                            q_s_arr(i, j, 0) = q_s_new;
                            
                            // Return for reduction: sum T_s, max T_s, sum q_s, max q_s
                            return {t_s_new, std::abs(t_s_new), q_s_new, std::abs(q_s_new)};
                        });
                    
                    // Copy results from device to host
                    amrex::Gpu::synchronize();
                    auto prog_reduce_tuple = prog_reduce_data.value(prog_reduce_ops);
                    t_s_sum_box = amrex::get<0>(prog_reduce_tuple);
                    t_s_max_box = amrex::get<1>(prog_reduce_tuple);
                    q_s_sum_box = amrex::get<2>(prog_reduce_tuple);
                    q_s_max_box = amrex::get<3>(prog_reduce_tuple);
                    
                    // Accumulate into global results
                    t_s_sum += t_s_sum_box;
                    t_s_max_val = std::max(t_s_max_val, t_s_max_box);
                    q_s_sum += q_s_sum_box;
                    q_s_max_val = std::max(q_s_max_val, q_s_max_box);
                    n_prog_columns += n_cols_box;
                }
                
                // Compute mean values from sums
                if (n_prog_columns > 0) {
                    t_s_mean = t_s_sum / static_cast<amrex::Real>(n_prog_columns);
                    t_s_max = t_s_max_val;
                    q_s_mean = q_s_sum / static_cast<amrex::Real>(n_prog_columns);
                    q_s_max = q_s_max_val;
                } else {
                    t_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
                    t_s_max = std::numeric_limits<amrex::Real>::quiet_NaN();
                    q_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
                    q_s_max = std::numeric_limits<amrex::Real>::quiet_NaN();
                }
            } else {
                // Noah-MP is active; skip prognostic update for this level
                // Leave t_s and q_s as populated by LSM passthrough
                t_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
                t_s_max = std::numeric_limits<amrex::Real>::quiet_NaN();
                q_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
                q_s_max = std::numeric_limits<amrex::Real>::quiet_NaN();
            }
        }
        // equivalent to a single-column value for spatially UNIFORM atmospheres
        // (as in the current SW_ClearSky_Analytical / LW_Isothermal RegTests).
        // Cloud layer and scattering tests are ALSO spatially uniform (identical
        // cloud/tau/scattering parameters applied to every column), so the
        // domain-averaged value still equals the true single-column flux there.
        // True horizontal heterogeneity (e.g., patchy clouds varying by column)
        // remains deferred to future work; see RAD_DEVELOPMENT.md.
        if (n_columns_total > 0) {
            SW_surface = sw_surface_sum / static_cast<amrex::Real>(n_columns_total);
            F_up_surface = lw_net_sum / static_cast<amrex::Real>(n_columns_total);
        }
        heating_rate_max = max_heating_global;

        // Compute SEB residual mean from sum
        if (rad_choice.seb_diagnostic_enable && rad_choice.seb_enable && n_seb_columns > 0) {
            seb_residual_mean = seb_residual_sum / static_cast<amrex::Real>(n_seb_columns);
        } else {
            // When feature is disabled, use NaN for backward compatibility
            seb_residual_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
            seb_residual_max = std::numeric_limits<amrex::Real>::quiet_NaN();
        }

        // For LW: in isothermal test, override with analytical value
        if (rad_choice.lw_enabled && rad_choice.isothermal_test) {
            amrex::Real sigma = 5.670374419e-8;
            amrex::Real T = rad_choice.T_iso_K;
            amrex::Real I_thermal = sigma * T * T * T * T;
            F_down_toa = I_thermal;
            F_up_surface = I_thermal;
            heating_rate_max = 0.0;
        }
    }

    // Logging output
    if (rad_choice.verbosity >= 1 && ParallelDescriptor::IOProcessor()) {
        Print() << "Radiation diagnostics at step " << nstep << ":\n"
                << "  SW TOA = " << SW_TOA << " W/m^2\n"
                << "  SW surface = " << SW_surface << " W/m^2\n"
                << "  LW up (surface) = " << F_up_surface << " W/m^2\n"
                << "  LW down (TOA) = " << F_down_toa << " W/m^2\n"
                << "  Max heating rate = " << heating_rate_max << " K/s\n";
        if (rad_choice.seb_diagnostic_enable && std::isfinite(seb_residual_mean)) {
            Print() << "  SEB residual (mean) = " << seb_residual_mean << " W/m^2\n"
                    << "  SEB residual (max) = " << seb_residual_max << " W/m^2\n";
        }
        if (rad_choice.seb_prognostic_enable && std::isfinite(t_s_mean)) {
            Print() << "  Surface temperature (mean) = " << t_s_mean << " K\n"
                    << "  Surface temperature (max) = " << t_s_max << " K\n"
                    << "  Surface moisture (mean) = " << q_s_mean << " kg/kg\n"
                    << "  Surface moisture (max) = " << q_s_max << " kg/kg\n";
        }
    }

    rad_diag.append(nstep, time_step, call_site, SW_surface, SW_TOA,
                    F_up_surface, F_down_toa, heating_rate_max,
                    seb_residual_mean, seb_residual_max,
                    t_s_mean, t_s_max, q_s_mean, q_s_max);
}

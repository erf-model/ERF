/**
 * @file ERF_DustLayer.cpp
 * @brief Implementation of the DustLayer container class.
 */

#include <ERF_DustLayer.H>
#include <ERF_DustPrerequisites.H>
#include <ERF_DustGrid.H>
#include <ERF_DustSurfaceReader.H>
#include <ERF_PhreeqcReader.H>
#include <ERF_DustThreshold.H>
#include <ERF_DustEmission.H>
#include <ERF_DustSuppression.H>
#include <ERF_DustWindExtract.H>
#include <ERF.H>
#include <ERF_SurfaceLayer.H>
#include <AMReX_Print.H>
#include <cmath>

void DustLayer::initialize(const ERF&          erf,
                           const SurfaceLayer* surface_layer,
                           const DustParams&   dust_params)
{
    // Step 1: Call verify_dust_prerequisites
    verify_dust_prerequisites(erf, surface_layer, dust_params);

    // Step 2: Create dust grid
    m_dg = create_dust_grid(erf.boxArray(0), erf.DistributionMap(0),
                            erf.Geom(0), dust_params.grid_ratio);

    if (dust_params.dust_debug) {
        amrex::Box dust_domain = m_dg.ba.minimalBox();
        int dust_nx = dust_domain.length(0);
        int dust_ny = dust_domain.length(1);
        amrex::Print() << "[DUST DEBUG] Created dust grid: " 
                       << dust_nx << "x" << dust_ny << "x1 cells, "
                       << "grid_ratio=" << dust_params.grid_ratio << "\n";
    }

    // Phase 2: grid creation debug output
    if (dust_params.dust_debug) {
        const amrex::Box& dom = m_dg.geom.Domain();
        amrex::Print() << "[DUST DEBUG] Phase 2: grid extent ["
                       << dom.smallEnd() << "," << dom.bigEnd()
                       << "] grid_ratio=" << dust_params.grid_ratio
                       << " boxes=" << m_dg.ba.size() << "\n";
    }

    // Step 3: Store dust_params
    m_params = dust_params;

    // Step 4: Allocate each MultiFab with 1 ghost cell (except z which has 0)
    amrex::IntVect ng(1, 1, 0);

    dust_ustar_t = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
    dust_soil_type = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
    dust_silt_fraction = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
    dust_crust_index = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
    dust_moisture_flag = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
    dust_suppression = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
    dust_emission_flux = std::make_unique<amrex::MultiFab>(
        m_dg.ba, m_dg.dm, dust_params.n_size_bins, ng);

    if (dust_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Allocated MultiFabs: "
                       << "dust_ustar_t (1 comp), "
                       << "dust_soil_type (1 comp), "
                       << "dust_silt_fraction (1 comp), "
                       << "dust_crust_index (1 comp), "
                       << "dust_moisture_flag (1 comp), "
                       << "dust_suppression (1 comp), "
                       << "dust_emission_flux (" << dust_params.n_size_bins << " comp)\n";
    }

    // Step 5: Fill initial values using setVal

    // dust_ustar_t: Bagnold threshold for the coarsest bin (bin 0)
    // u*_t = A * sqrt(rho_p * g * d / rho_a)
    // Bagnold (1941). The Physics of Blown Sand and Desert Dunes, Methuen, London.
    amrex::Real g = 9.81;  // m/s^2
    amrex::Real rho_a = 1.225;  // kg/m^3 (standard air density)
    amrex::Real d_bin0 = dust_params.bin_diameter_um[0] * 1.0e-6;  // convert um to m
    amrex::Real ustar_t = dust_params.threshold_A_coeff
                        * std::sqrt(dust_params.particle_density * g * d_bin0 / rho_a);
    dust_ustar_t->setVal(ustar_t);

    if (dust_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Set dust_ustar_t (Bagnold threshold bin 0): "
                       << "u*_t=" << ustar_t << " m/s, "
                       << "d=" << dust_params.bin_diameter_um[0] << " um, "
                       << "rho_p=" << dust_params.particle_density << " kg/m^3, "
                       << "A=" << dust_params.threshold_A_coeff << "\n";
    }

    // Phase 3: surface maps debug output (before population from files)
    if (dust_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Phase 3: silt_frac min="
                       << dust_params.silt_fraction
                       << " (default, will be updated from file if provided)\n";
    }

    // dust_soil_type: 0.0 (undefined; set by surface reader in Phase 3)
    dust_soil_type->setVal(0.0);

    // dust_silt_fraction: silt_fraction from params
    dust_silt_fraction->setVal(dust_params.silt_fraction);

    // dust_crust_index: crust_index from params
    dust_crust_index->setVal(dust_params.crust_index);

    // dust_moisture_flag: 0.0 (dry surface at initialization)
    dust_moisture_flag->setVal(0.0);

    // dust_suppression: 0.0 (no suppression at initialization)
    dust_suppression->setVal(0.0);

    // dust_emission_flux: 0.0 (zero emission at initialization)
    dust_emission_flux->setVal(0.0);

    if (dust_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Set initial values: "
                       << "dust_soil_type=0.0, "
                       << "dust_silt_fraction=" << dust_params.silt_fraction << ", "
                       << "dust_crust_index=" << dust_params.crust_index << ", "
                       << "dust_moisture_flag=0.0, "
                       << "dust_suppression=0.0, "
                       << "dust_emission_flux=0.0\n";
    }

    // Allocate efflorescence and base u*_t MultiFabs.
    dust_efflor     = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
    dust_ustar_base = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
    dust_efflor->setVal(0.0);
    // Store the Bagnold u*_t computed from DustParams as the base value.
    // update_ustar_t_from_chemistry modifies dust_ustar_t; dust_ustar_base is read-only.
    //dust_ustar_base->Copy(*dust_ustar_t, 0, 0, 1, amrex::IntVect(1,1,0));
    amrex::MultiFab::Copy(*dust_ustar_base, *dust_ustar_t, 0, 0, 1, amrex::IntVect(1,1,0));

    if (dust_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Allocated internal MultiFabs: "
                       << "dust_efflor (1 comp), "
                       << "dust_ustar_base (1 comp)\n";
    }

    // Phase 5: Allocate surface moisture MultiFab (Phase 9 will update it from ERF surface layer).
    dust_surf_moist = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
    dust_surf_moist->setVal(0.0);

    // Phase 6: Allocate u* input MultiFab. Phase 6: filled from DustParams::test_ustar.
    // Phase 9: replaced by MRF surface layer u* extraction.
    dust_ustar_in = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
    dust_ustar_in->setVal(dust_params.test_ustar);

    if (dust_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Phase 6: dust_ustar_in (placeholder) = "
                       << dust_params.test_ustar << " m/s\n";
    }

    // Phase 8: Allocate re-treatment flag. Set to 0 at startup (no re-treatment needed initially).
    dust_retreat_flag = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
    dust_retreat_flag->setVal(0.0);

    if (dust_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Phase 8: dust_retreat_flag allocated, "
                       << "supp_tau_base_s=" << dust_params.supp_tau_base_s
                       << " s, test_surf_temp_K=" << dust_params.test_surf_temp_K
                       << " K, test_wind_speed=" << dust_params.test_wind_speed
                       << " m/s\n";
    }

    // Phase 9: Allocate atmospheric field MultiFabs
    dust_wind_ref = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 2, amrex::IntVect(1,1,0));
    dust_pblh     = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
    dust_tsfc     = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
    dust_wind_ref->setVal(0.0);
    dust_pblh->setVal(0.0);
    dust_tsfc->setVal(dust_params.test_surf_temp_K);

    if (dust_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Phase 9: dust_wind_ref(2), dust_pblh, dust_tsfc allocated\n";
    }

    // Populate surface MultiFabs from external rasters if paths are given in dust_params.
    // MultiFabs retain their setVal defaults above if paths are empty.
    populate_dust_surface_maps(*dust_soil_type, *dust_silt_fraction,
                               *dust_crust_index, *dust_moisture_flag,
                               *dust_suppression, m_dg, dust_params);

    if (dust_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Surface maps populated: "
                       << "soil_type_file=\"" << dust_params.soil_type_file << "\", "
                       << "silt_fraction_file=\"" << dust_params.silt_fraction_file << "\", "
                       << "crust_index_file=\"" << dust_params.crust_index_file << "\", "
                       << "moisture_flag_file=\"" << dust_params.moisture_flag_file << "\", "
                       << "suppression_file=\"" << dust_params.suppression_file << "\"\n";
    }

    // Phase 3: surface maps debug output (after population from files)
    if (dust_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Phase 3: silt_frac min="
                       << dust_silt_fraction->min(0)
                       << " max=" << dust_silt_fraction->max(0)
                       << " soil_type_max=" << dust_soil_type->max(0) << "\n";
    }

    // Load blast schedule if a file is specified.
    // ERF_DustBlastSchedule.H adapts ERF_IgnitionSchedule.H from Source/Fire/.
#ifdef ERF_USE_DUST
    if (!dust_params.blast_schedule_file.empty()) {
        load_blast_schedule(dust_params.blast_schedule_file, m_blast_schedule);
        m_has_blast_schedule = !m_blast_schedule.empty();
        if (m_has_blast_schedule) {
            amrex::Print() << "[DUST] Blast schedule loaded from: "
                           << dust_params.blast_schedule_file << "\n";
        }
        if (dust_params.dust_debug && m_has_blast_schedule) {
            amrex::Print() << "[DUST DEBUG] Phase 7: blast schedule loaded, "
                           << m_blast_schedule.events.size() << " events\n";
        }
    }
#endif

    // Phase 5: Compute initial u*_t. At startup all modifiers are zero so u*_t equals
    // dust_ustar_base where crust_index is also zero, floored to USTAR_T_MIN.
    recompute_dust_ustar_t(*dust_ustar_t, *dust_ustar_base, *dust_crust_index,
                           *dust_efflor, *dust_surf_moist, *dust_suppression,
                           dust_params.alpha_crust, dust_params.alpha_efflor);

    if (dust_params.dust_debug) {
        amrex::Real ut_min = dust_ustar_t->min(0);
        amrex::Real ut_max = dust_ustar_t->max(0);
        amrex::Print() << "[DUST DEBUG] Phase 5: u*_t after init modifiers: "
                       << "min=" << ut_min << " m/s, max=" << ut_max << " m/s "
                       << "(USTAR_T_MIN=" << DustThresholdConst::USTAR_T_MIN << " m/s)\n";
    }

    // Phase 5: initial u*_t debug output
    if (dust_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Phase 5: u*_t init min="
                       << dust_ustar_t->min(0)
                       << " max=" << dust_ustar_t->max(0) << " [m/s]\n";
    }

    // Phase 8: suppression and retreat flag debug output
    if (dust_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Phase 8: suppression_max="
                       << dust_suppression->max(0)
                       << " tau_base=" << dust_params.supp_tau_base_s << " s\n";
    }

    // Step 6: Print status message
    amrex::Box dust_domain = m_dg.ba.minimalBox();
    int dust_nx = dust_domain.length(0);
    int dust_ny = dust_domain.length(1);
    amrex::Print() << "[DUST] DustLayer initialized: grid_ratio="
                   << m_dg.grid_ratio << ", dust cells="
                   << dust_nx << "x" << dust_ny << "x1, "
                   << "n_size_bins=" << dust_params.n_size_bins << ", "
                   << "z0_dust=" << dust_params.z0_dust << " m\n";

    if (dust_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] PHREEQC configuration: "
                       << "phreeqc_output_file=\"" << dust_params.phreeqc_output_file << "\", "
                       << "phreeqc_update_interval_s=" << dust_params.phreeqc_update_interval_s << " s\n";
    }
}


void DustLayer::advance(amrex::Real             dt,
                        const DustParams&       dust_params,
                        const SurfaceLayer*     surface_layer,
                        const amrex::MultiFab*  xvel_mf,
                        const amrex::MultiFab*  yvel_mf,
                        const amrex::MultiFab*  z_phys_cc_mf,
                        int                     nz)
{
    ++m_step;
    m_time += dt;

    // Phase 9: Extract wind and surface fields from atmospheric solver
    bool have_atm = (surface_layer && xvel_mf && yvel_mf && z_phys_cc_mf && nz > 0);

    if (have_atm) {
        // u* from SurfaceLayer. Accessor: get_u_star(0).
        if (surface_layer->get_u_star(0)) {
            fill_dust_ustar_from_surface_layer(*dust_ustar_in,
                                               *surface_layer->get_u_star(0), m_dg);
        }
        // 10m wind by vertical interpolation. Algorithm from ERF_FireWindExtract.cpp.
        fill_dust_wind_from_interpolation(*dust_wind_ref, *xvel_mf, *yvel_mf,
                                          *z_phys_cc_mf, m_dg, m_params.zref, nz);
        // Surface temperature. Accessor: get_t_surf(0), NOT get_t_star.
        if (surface_layer->get_t_surf(0))
            fill_dust_scalar_from_atm(*dust_tsfc, *surface_layer->get_t_surf(0), m_dg);
        // PBL height. Accessor: get_pblh(0).
        if (surface_layer->get_pblh(0))
            fill_dust_scalar_from_atm(*dust_pblh, *surface_layer->get_pblh(0), m_dg);

        if (m_params.dust_debug)
            amrex::Print() << "[DUST DEBUG] Phase 9: step=" << m_step
                           << " u*_max=" << dust_ustar_in->max(0)
                           << " u_10m_max=" << dust_wind_ref->max(0)
                           << " PBLH_max=" << dust_pblh->max(0) << "\n";
    } else {
        // Placeholder path for standalone regression tests (Phases 1-8).
        dust_ustar_in->setVal(m_params.test_ustar);
        dust_tsfc->setVal(m_params.test_surf_temp_K);
        dust_wind_ref->setVal(m_params.test_wind_speed);
        if (m_params.dust_debug)
            amrex::Print() << "[DUST DEBUG] Phase 9: placeholder path"
                           << " test_ustar=" << m_params.test_ustar << "\n";
    }

    // Scalar values for Phase 8 suppression decay. Use max() as representative.
    // Phase 13 will use per-cell extracted fields.
    amrex::Real T_sfc = have_atm ? dust_tsfc->max(0) : m_params.test_surf_temp_K;
    amrex::Real u_10m = have_atm ? dust_wind_ref->max(0) : m_params.test_wind_speed;

    // Physics inserted in Phases 5-13. PHREEQC reader called here (Phase 4).
    // Call PHREEQC reader if update interval has elapsed.
    // The interval is set by dust_params.phreeqc_update_interval_s.
    // File-based coupling is appropriate because geochemical processes
    // evolve on timescales of days to weeks, much longer than the
    // atmospheric timestep.
    bool do_phreeqc = (m_last_phreeqc_update < 0.0) ||
                      (m_time - m_last_phreeqc_update >=
                       dust_params.phreeqc_update_interval_s);

    if (do_phreeqc && !dust_params.phreeqc_output_file.empty()) {
        if (dust_params.dust_debug) {
            amrex::Print() << "[DUST DEBUG] PHREEQC update triggered at step="
                           << m_step << ", time=" << m_time << " s\n";
        }
        update_dust_from_phreeqc(*dust_ustar_t,
                                 *dust_ustar_base,
                                 *dust_crust_index,
                                 *dust_silt_fraction,
                                 *dust_efflor,
                                 *dust_suppression,
                                 *dust_emission_flux,
                                 m_dg,
                                 dust_params);
        m_last_phreeqc_update = m_time;
        if (dust_params.dust_debug) {
            amrex::Print() << "[DUST DEBUG] PHREEQC update completed\n";
        }
    } else if (dust_params.dust_debug && !dust_params.phreeqc_output_file.empty()) {
        amrex::Print() << "[DUST DEBUG] Step=" << m_step << ", time=" << m_time << " s. "
                       << "PHREEQC update not due (next update at "
                       << m_last_phreeqc_update + dust_params.phreeqc_update_interval_s << " s)\n";
    }

    // Phase 5: Copy moisture_flag to surf_moist. Replaced by ERF surface layer fields in Phase 9.
    amrex::MultiFab::Copy(*dust_surf_moist, *dust_moisture_flag, 0, 0, 1, amrex::IntVect(1,1,0));

    // Phase 8: Advance suppression coverage decay and set re-treatment flag.
    // Phase 9 uses extracted T and wind fields instead of test values.
    if (dt > 0.0) {
        advance_dust_suppression(*dust_suppression,
                                 *dust_retreat_flag,
                                 T_sfc,
                                 u_10m,
                                 dt,
                                 m_params.supp_tau_base_s);
    }

    if (m_params.dust_debug) {
        amrex::Real s_max = dust_suppression->max(0);
        amrex::Real r_sum = dust_retreat_flag->sum(0);
        amrex::Print() << "[DUST DEBUG] Phase 8: suppression coverage max="
                       << s_max << ", retreat_flag sum=" << r_sum
                       << " at step=" << m_step << "\n";
    }

    // Phase 5: Recompute u*_t after any PHREEQC or moisture changes this timestep.
    recompute_dust_ustar_t(*dust_ustar_t, *dust_ustar_base, *dust_crust_index,
                           *dust_efflor, *dust_surf_moist, *dust_suppression,
                           m_params.alpha_crust, m_params.alpha_efflor);

    if (m_params.dust_debug) {
        amrex::Real ut_min = dust_ustar_t->min(0);
        amrex::Real ut_max = dust_ustar_t->max(0);
        amrex::Print() << "[DUST DEBUG] Phase 5: u*_t at step=" << m_step
                       << " min=" << ut_min << " max=" << ut_max << " [m/s]\n";
    }

    // Phase 5: u*_t per step debug output
    if (m_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Phase 5: step=" << m_step
                       << " u*_t min=" << dust_ustar_t->min(0)
                       << " max=" << dust_ustar_t->max(0) << "\n";
    }

    // Phase 6: Note: u* is now updated in Phase 9 section above (have_atm path).
    // dust_ustar_in is set either from extracted MRF surface layer or test value.

    // Phase 6: Compute emission flux using current u*, u*_t, and silt fraction.
    // Zero emission when test_ustar = 0 (default until Phase 9).
    compute_dust_emission_flux(*dust_emission_flux,
                               *dust_ustar_t,
                               *dust_ustar_in,
                               *dust_silt_fraction,
                               m_params.n_size_bins,
                               m_params.rho_air);

    if (m_params.dust_debug) {
        amrex::Real f_max = dust_emission_flux->max(0);
        amrex::Real f_sum = dust_emission_flux->sum(0);
        amrex::Print() << "[DUST DEBUG] Phase 6: emission_flux bin0 at step=" << m_step
                       << " max=" << f_max << " sum=" << f_sum
                       << " [kg/m^2/s]\n";
    }

    // Phase 6: emission flux per step debug output
    if (m_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Phase 6: step=" << m_step
                       << " emission_flux_bin0 max=" << dust_emission_flux->max(0) << "\n";
    }

    // Apply any blast events scheduled in this timestep.
    // Blast injection adds to emission flux computed above.
    // ERF_DustBlastSchedule.H documents the CSV format and injection formula.
#ifdef ERF_USE_DUST
    if (m_has_blast_schedule && dt > 0.0) {
        apply_blast_schedule(*dust_emission_flux,
                             m_dg,
                             m_blast_schedule,
                             m_time,
                             m_time - dt,
                             dt,
                             m_params.n_size_bins,
                             m_params.blast_reactivity);

        if (m_params.dust_debug) {
            amrex::Real f_max = dust_emission_flux->max(0);
            amrex::Print() << "[DUST DEBUG] Phase 7: emission_flux bin0 after blast step="
                           << m_step << " max=" << f_max << " [kg/m^2/s]\n";
        }
    }
#endif

    // Phase 8: suppression per step debug output
    if (m_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Phase 8: step=" << m_step
                       << " suppression_max=" << dust_suppression->max(0)
                       << " retreat_sum=" << dust_retreat_flag->sum(0) << "\n";
    }
}
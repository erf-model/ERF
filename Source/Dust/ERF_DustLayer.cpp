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
#include <ERF_DustAtmCoupling.H>
#include <ERF_DustMSHA.H>
#include <ERF_DustMSHAOutput.H>
#include <ERF.H>
#include <ERF_SurfaceLayer.H>
#include <AMReX_Print.H>
#include <cmath>

void
DustLayer::initialize(
  const ERF& erf,
  const SurfaceLayer* surface_layer,
  const DustParams& dust_params)
{
  // Step 1: Call verify_dust_prerequisites
  verify_dust_prerequisites(erf, surface_layer, dust_params);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 1: Verified all prerequisites\n";
  }

  // Step 2: Create dust grid
  m_dg = create_dust_grid(
    erf.boxArray(0), erf.DistributionMap(0), erf.Geom(0),
    dust_params.grid_ratio);

  if (dust_params.dust_debug) {
    amrex::Box dust_domain = m_dg.ba.minimalBox();
    int dust_nx = dust_domain.length(0);
    int dust_ny = dust_domain.length(1);
    amrex::Print() << "[DUST DEBUG] Created dust grid: " << dust_nx << "x"
                   << dust_ny << "x1 cells, "
                   << "grid_ratio=" << dust_params.grid_ratio << "\n";
    const amrex::Box& dom = m_dg.geom.Domain();
    amrex::Print() << "[DUST DEBUG] Phase 2: grid extent [" << dom.smallEnd()
                   << "," << dom.bigEnd()
                   << "] grid_ratio=" << dust_params.grid_ratio
                   << " boxes=" << m_dg.ba.size() << "\n";
  }

  // Step 3: Store dust_params
  m_params = dust_params;

  // Step 4: Allocate each MultiFab with 1 ghost cell (except z which has 0)
  amrex::IntVect ng(1, 1, 0);

  dust_ustar_t = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
  dust_soil_type = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
  dust_silt_fraction =
    std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
  dust_crust_index = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
  dust_moisture_flag =
    std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
  dust_suppression = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
  dust_emission_flux = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, dust_params.n_size_bins, ng);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Allocated MultiFabs: "
                   << "dust_ustar_t (1 comp), " << "dust_soil_type (1 comp), "
                   << "dust_silt_fraction (1 comp), "
                   << "dust_crust_index (1 comp), "
                   << "dust_moisture_flag (1 comp), "
                   << "dust_suppression (1 comp), " << "dust_emission_flux ("
                   << dust_params.n_size_bins << " comp)\n";
  }

  // Step 5: Fill initial values
  amrex::Real g = 9.81;
  amrex::Real rho_a = 1.225;
  amrex::Real d_bin0 = dust_params.bin_diameter_um[0] * 1.0e-6;
  amrex::Real ustar_t =
    dust_params.threshold_A_coeff *
    std::sqrt(dust_params.particle_density * g * d_bin0 / rho_a);
  dust_ustar_t->setVal(ustar_t);

  if (dust_params.dust_debug) {
    amrex::Print()
      << "[DUST DEBUG] Set dust_ustar_t (Bagnold threshold bin 0): " << "u*_t="
      << ustar_t << " m/s, " << "d=" << dust_params.bin_diameter_um[0]
      << " um, " << "rho_p=" << dust_params.particle_density << " kg/m^3, "
      << "A=" << dust_params.threshold_A_coeff << "\n";
  }

  dust_soil_type->setVal(0.0);
  dust_silt_fraction->setVal(dust_params.silt_fraction);
  dust_crust_index->setVal(dust_params.crust_index);
  dust_moisture_flag->setVal(0.0);
  dust_suppression->setVal(0.0);
  dust_emission_flux->setVal(0.0);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Set initial values: "
                   << "dust_soil_type=0.0, "
                   << "dust_silt_fraction=" << dust_params.silt_fraction << ", "
                   << "dust_crust_index=" << dust_params.crust_index << ", "
                   << "dust_moisture_flag=0.0, " << "dust_suppression=0.0, "
                   << "dust_emission_flux=0.0\n";
  }

  // Allocate efflorescence and base u*_t MultiFabs
  dust_efflor = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, 1, amrex::IntVect(1, 1, 0));
  dust_ustar_base = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, 1, amrex::IntVect(1, 1, 0));
  dust_efflor->setVal(0.0);
  amrex::MultiFab::Copy(
    *dust_ustar_base, *dust_ustar_t, 0, 0, 1, amrex::IntVect(1, 1, 0));

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Allocated internal MultiFabs: "
                   << "dust_efflor (1 comp), " << "dust_ustar_base (1 comp)\n";
  }

  // Phase 5: surface moisture
  dust_surf_moist = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, 1, amrex::IntVect(1, 1, 0));
  dust_surf_moist->setVal(0.0);

  // Phase 6: u* input (placeholder; Phase 9 replaces with MRF surface layer)
  dust_ustar_in = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, 1, amrex::IntVect(1, 1, 0));
  dust_ustar_in->setVal(dust_params.test_ustar);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 6: dust_ustar_in (placeholder) = "
                   << dust_params.test_ustar << " m/s\n";
  }

  // Phase 8: re-treatment flag
  dust_retreat_flag = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, 1, amrex::IntVect(1, 1, 0));
  dust_retreat_flag->setVal(0.0);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 8: dust_retreat_flag allocated, "
                   << "supp_tau_base_s=" << dust_params.supp_tau_base_s
                   << " s, test_surf_temp_K=" << dust_params.test_surf_temp_K
                   << " K, test_wind_speed=" << dust_params.test_wind_speed
                   << " m/s\n";
  }

  // Phase 9: atmospheric fields
  dust_wind_ref = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, 2, amrex::IntVect(1, 1, 0));
  dust_pblh = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, 1, amrex::IntVect(1, 1, 0));
  dust_tsfc = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, 1, amrex::IntVect(1, 1, 0));
  dust_wind_ref->setVal(0.0);
  dust_pblh->setVal(0.0);
  dust_tsfc->setVal(dust_params.test_surf_temp_K);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 9: dust_wind_ref(2), dust_pblh, "
                      "dust_tsfc allocated\n";
  }

  // Phase 9: confirm extraction fields are ready
  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 9: zref=" << dust_params.zref
                   << " m for wind extraction\n";
  }

  // Phase 10: Allocate coarsened emission flux on atmospheric grid.
  // BoxArray and DistributionMapping are from ERF level-0.
  // Use the 2D slab: set k-range to [0,0].
  {
    amrex::BoxArray ba_atm = erf.boxArray(0);
    amrex::Vector<amrex::Box> bl;
    for (int b = 0; b < ba_atm.size(); ++b) {
      amrex::Box bx = ba_atm[b];
      bx.setSmall(2, 0);
      bx.setBig(2, 0);
      bl.push_back(bx);
    }
    amrex::BoxArray ba2d(amrex::BoxList(std::move(bl)));
    dust_flux_atm = std::make_unique<amrex::MultiFab>(
      ba2d, erf.DistributionMap(0), 1, amrex::IntVect(1, 1, 0));
    dust_flux_atm->setVal(0.0);
  }

  // Use the second passive scalar slot (RhoScalar_comp + 1).
  // NSCALARS is 2 when ERF_USE_DUST is defined, so this index is
  // within bounds: nvars = NDRY + NSCALARS = 3 + 2 = 5, valid indices 0-4.
  m_dust_scalar_comp = RhoScalar_comp + 1;

  if (dust_params.dust_debug) {
    amrex::Print()
      << "[DUST DEBUG] Phase 10: dust_flux_atm allocated on atm grid."
      << " dust_scalar_comp=" << m_dust_scalar_comp << "\n";
  }

  // Phase 12: allocate deposition accumulators
  dust_deposition_rate = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_deposition_rate->setVal(0.0);

  // dep_flux_atm: 2D slab on atm grid, same BoxArray as dust_flux_atm.
  dep_flux_atm = std::make_unique<amrex::MultiFab>(
    dust_flux_atm->boxArray(),
    dust_flux_atm->DistributionMap(),
    1, amrex::IntVect(1,1,0));
  dep_flux_atm->setVal(0.0);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 12: dust_deposition_rate"
                   << " and dep_flux_atm allocated\n";
  }

  // Phase 13: allocate return fields from 3D solver
  dust_conc_sfc = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_conc_sfc->setVal(0.0);
  dust_surf_moist->setVal(0.0);  // Initialize the existing dust_surf_moist for Phase 13

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 13: dust_conc_sfc and"
                   << " dust_surf_moist allocated on dust grid\n"
                   << "[DUST DEBUG] Phase 13: loading_feedback_coeff="
                   << dust_params.loading_feedback_coeff
                   << " use_dynamic_moisture=" << dust_params.use_dynamic_moisture
                   << "\n";
  }

  // Phase 17: PM size fraction MultiFabs on dust grid.
  dust_pm25 = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_pm10 = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_pm25_24h = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_pm10_24h = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_pm25_exceed = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_pm10_exceed = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  for (auto* mf : {dust_pm25.get(), dust_pm10.get(),
                   dust_pm25_24h.get(), dust_pm10_24h.get(),
                   dust_pm25_exceed.get(), dust_pm10_exceed.get()})
    mf->setVal(0.0);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 17: PM2.5/PM10 MultiFabs allocated\n"
                   << "[DUST DEBUG] Phase 17: naaqs_file="
                   << dust_params.dust_naaqs_file << "\n";
  }

  // Phase 18: MSHA worker exposure MultiFabs on dust grid.
  dust_msha_dose = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_msha_twa = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_msha_exceed = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_msha_shift_twa = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  for (auto* mf : {dust_msha_dose.get(), dust_msha_twa.get(),
                   dust_msha_exceed.get(), dust_msha_shift_twa.get()})
    mf->setVal(0.0);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 18: MSHA exposure MultiFabs allocated\n"
                   << "[DUST DEBUG] Phase 18: PEL=" << dust_params.msha_pel_mg_m3 << " mg/m3 (30 CFR 56.5001)\n"
                   << "[DUST DEBUG] Phase 18: shift_duration=" << dust_params.msha_shift_duration_s << " s\n"
                   << "[DUST DEBUG] Phase 18: n_receptors=" << dust_params.msha_receptor_names.size() << "\n";
    for (int r = 0; r < (int)dust_params.msha_receptor_names.size(); ++r) {
      amrex::Print() << "[DUST DEBUG] Phase 18: receptor " << r << ": "
                     << dust_params.msha_receptor_names[r] << " ("
                     << dust_params.msha_receptor_x[r] << ", "
                     << dust_params.msha_receptor_y[r] << ")\n";
    }
  }

  // Populate surface MultiFabs from external rasters
  populate_dust_surface_maps(
    *dust_soil_type, *dust_silt_fraction, *dust_crust_index,
    *dust_moisture_flag, *dust_suppression, m_dg, dust_params);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Surface maps populated: "
                   << "soil_type_file=\"" << dust_params.soil_type_file
                   << "\", " << "silt_fraction_file=\""
                   << dust_params.silt_fraction_file << "\", "
                   << "crust_index_file=\"" << dust_params.crust_index_file
                   << "\", " << "moisture_flag_file=\""
                   << dust_params.moisture_flag_file << "\", "
                   << "suppression_file=\"" << dust_params.suppression_file
                   << "\"\n";
    amrex::Print() << "[DUST DEBUG] Phase 3: silt_frac min="
                   << dust_silt_fraction->min(0)
                   << " max=" << dust_silt_fraction->max(0)
                   << " soil_type_max=" << dust_soil_type->max(0) << "\n";
  }

  // Load blast schedule if specified
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

  // Phase 5: Compute initial u*_t
  recompute_dust_ustar_t(
    *dust_ustar_t, *dust_ustar_base, *dust_crust_index, *dust_efflor,
    *dust_surf_moist, *dust_suppression, dust_params.alpha_crust,
    dust_params.alpha_efflor);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 5: u*_t after init modifiers: "
                   << "min=" << dust_ustar_t->min(0) << " m/s, "
                   << "max=" << dust_ustar_t->max(0) << " m/s "
                   << "(USTAR_T_MIN=" << DustThresholdConst::USTAR_T_MIN
                   << " m/s)\n";
    amrex::Print() << "[DUST DEBUG] Phase 8: suppression_max="
                   << dust_suppression->max(0)
                   << " tau_base=" << dust_params.supp_tau_base_s << " s\n";
  }

  // Step 6: Print status message
  amrex::Box dust_domain = m_dg.ba.minimalBox();
  int dust_nx = dust_domain.length(0);
  int dust_ny = dust_domain.length(1);
  amrex::Print() << "[DUST] DustLayer initialized: grid_ratio="
                 << m_dg.grid_ratio << ", dust cells=" << dust_nx << "x"
                 << dust_ny << "x1, "
                 << "n_size_bins=" << dust_params.n_size_bins << ", "
                 << "z0_dust=" << dust_params.z0_dust << " m\n";

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] PHREEQC configuration: "
                   << "phreeqc_output_file=\""
                   << dust_params.phreeqc_output_file << "\", "
                   << "phreeqc_update_interval_s="
                   << dust_params.phreeqc_update_interval_s << " s\n";
  }

  // Phase 11: Log bin_diameters and transport mode
  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 11: bin_diameters [um]:";
    for (auto d : dust_params.bin_diameters)
      amrex::Print() << " " << d * 1e6;
    amrex::Print() << "\n"
                   << "[DUST DEBUG] Phase 11: transport_bins_separately="
                   << dust_params.transport_bins_separately << "\n";
  }

  // Phase 17: Log bin PM classification
  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 17: bin PM classification:\n";
    for (int b = 0; b < (int)dust_params.bin_diameters.size(); ++b) {
      amrex::Real d = dust_params.bin_diameters[b];
      amrex::Print() << "  bin " << b << " d=" << d*1e6 << " um"
                     << " PM2.5=" << (is_pm25(d)?"yes":"no")
                     << " PM10=" << (is_pm10(d)?"yes":"no") << "\n";
    }
  }

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] initialize complete:"
                   << " n_size_bins=" << dust_params.n_size_bins
                   << " grid_ratio=" << m_dg.grid_ratio
                   << " dust_scalar_comp=" << m_dust_scalar_comp
                   << " loading_feedback_coeff="
                   << dust_params.loading_feedback_coeff
                   << " use_dynamic_moisture="
                   << dust_params.use_dynamic_moisture << "\n";
  }

  // Phase 16: Write CSV header and initialize plotfile tracking
  write_dust_stats_header(m_params.dust_diag_file);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 16: dust_plot_int="
                   << dust_params.dust_plot_int
                   << " prefix=" << dust_params.dust_plot_prefix
                   << " diag_file=" << dust_params.dust_diag_file << "\n";
  }

+#if defined(ERF_USE_PARTICLES)
+  // Phase 19: Initialize Lagrangian super-particles
+  if (dust_params.enable_particles) {
+    // Store atmospheric geometry for particle use
+    m_geom_atm = erf.Geom(0);
+    
+    // Construct particle container on atmospheric grid
+    // Use dust BoxArray for particle distribution but atm geometry for positions
+    amrex::BoxArray ba_particles = amrex::convert(m_dg.ba, amrex::IntVect::TheZeroVector());
+    m_dust_pc = std::make_unique<ERFDustPC>(m_geom_atm, m_dg.dm, ba_particles);
+    
+    // Allocate source_map on dust grid [kg/m^2]
+    dust_source_map = std::make_unique<amrex::MultiFab>(
+        m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
+    dust_source_map->setVal(0.0);
+    
+    if (dust_params.dust_debug) {
+      amrex::Print() << "[DUST DEBUG] Phase 19: ERFDustPC initialized\n";
+    }
+  }
+#endif


void
DustLayer::advance(
  amrex::Real dt,
  const DustParams& dust_params,
  SurfaceLayer* surface_layer, // non-const: getters not marked const
  const amrex::MultiFab* xvel_mf,
  const amrex::MultiFab* yvel_mf,
  const amrex::MultiFab* z_phys_cc_mf,
  int nz)
{
  ++m_step;
  m_time += dt;

  // Phase 13/17/18: Consolidated entry debug output
  if (m_params.dust_debug) {
    amrex::Real cs = dust_conc_sfc ? dust_conc_sfc->max(0) * 1e9 : 0.0;
    amrex::Real p10 = dust_pm10 ? dust_pm10->max(0) : 0.0;
    amrex::Real twa = dust_msha_twa ? dust_msha_twa->max(0) : 0.0;
    amrex::Print() << "[DUST DEBUG] advance: step=" << m_step
                   << " conc_sfc=" << cs << " ug/m3"
                   << " PM10=" << p10 << " ug/m3"
                   << " MSHA_TWA=" << twa << " mg/m3\n";
  }


  // Phase 9: Extract wind and surface fields from atmospheric solver
  bool have_atm =
    (surface_layer && xvel_mf && yvel_mf && z_phys_cc_mf && nz > 0);

  if (have_atm) {
    if (surface_layer->get_u_star(0))
      fill_dust_ustar_from_surface_layer(
        *dust_ustar_in, *surface_layer->get_u_star(0), m_dg);
    fill_dust_wind_from_interpolation(
      *dust_wind_ref, *xvel_mf, *yvel_mf, *z_phys_cc_mf, m_dg, m_params.zref,
      nz);
    if (surface_layer->get_t_surf(0))
      fill_dust_scalar_from_atm(
        *dust_tsfc, *surface_layer->get_t_surf(0), m_dg);
    if (surface_layer->get_pblh(0))
      fill_dust_scalar_from_atm(*dust_pblh, *surface_layer->get_pblh(0), m_dg);

    if (m_params.dust_debug)
      amrex::Print() << "[DUST DEBUG] Phase 9: step=" << m_step
                     << " u*_max=" << dust_ustar_in->max(0)
                     << " u_10m_max=" << dust_wind_ref->max(0)
                     << " PBLH_max=" << dust_pblh->max(0) << "\n";
  } else {
    dust_ustar_in->setVal(m_params.test_ustar);
    dust_tsfc->setVal(m_params.test_surf_temp_K);
    dust_wind_ref->setVal(m_params.test_wind_speed);
    if (m_params.dust_debug)
      amrex::Print() << "[DUST DEBUG] Phase 9: placeholder path"
                     << " test_ustar=" << m_params.test_ustar << "\n";
  }

  // Additional Phase 9 debug output
  if (m_params.dust_debug && have_atm) {
    amrex::Print() << "[DUST DEBUG] Phase 9: T_sfc_max=" << dust_tsfc->max(0)
                   << " K  PBLH_max=" << dust_pblh->max(0) << " m\n";
  }

  amrex::Real T_sfc = have_atm ? dust_tsfc->max(0) : m_params.test_surf_temp_K;
  amrex::Real u_10m =
    have_atm ? dust_wind_ref->max(0) : m_params.test_wind_speed;

  // Phase 4: PHREEQC reader
  bool do_phreeqc =
    (m_last_phreeqc_update < 0.0) ||
    (m_time - m_last_phreeqc_update >= dust_params.phreeqc_update_interval_s);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 4: PHREEQC reader check at step="
                   << m_step << ", time=" << m_time << " s, "
                   << "do_phreeqc=" << do_phreeqc << "\n";
  }

  if (do_phreeqc && !dust_params.phreeqc_output_file.empty()) {
    if (dust_params.dust_debug)
      amrex::Print() << "[DUST DEBUG] PHREEQC update triggered at step="
                     << m_step << ", time=" << m_time << " s\n";
    update_dust_from_phreeqc(
      *dust_ustar_t, *dust_ustar_base, *dust_crust_index, *dust_silt_fraction,
      *dust_efflor, *dust_suppression, *dust_emission_flux, m_dg, dust_params);
    m_last_phreeqc_update = m_time;
    if (dust_params.dust_debug)
      amrex::Print() << "[DUST DEBUG] PHREEQC update completed\n";
  }

  // Phase 5: Copy moisture_flag to surf_moist
  amrex::MultiFab::Copy(
    *dust_surf_moist, *dust_moisture_flag, 0, 0, 1, amrex::IntVect(1, 1, 0));

  // Phase 8: Suppression decay
  if (dt > 0.0) {
    advance_dust_suppression(
      *dust_suppression, *dust_retreat_flag, T_sfc, u_10m, dt,
      m_params.supp_tau_base_s);
  }

  if (m_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 8: suppression coverage max="
                   << dust_suppression->max(0)
                   << ", retreat_flag sum=" << dust_retreat_flag->sum(0)
                   << " at step=" << m_step << "\n";
  }

  // Phase 5: Recompute u*_t
  recompute_dust_ustar_t(
    *dust_ustar_t, *dust_ustar_base, *dust_crust_index, *dust_efflor,
    *dust_surf_moist, *dust_suppression, m_params.alpha_crust,
    m_params.alpha_efflor);

  if (m_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 5: u*_t at step=" << m_step
                   << " min=" << dust_ustar_t->min(0)
                   << " max=" << dust_ustar_t->max(0) << " [m/s]\n";
  }

  // Phase 13: Apply Shao (2001) loading feedback to u*_t
  // Effective threshold: u*_t *= (1 + coeff * C_sfc)
  // Reference: Shao (2001), https://doi.org/10.1029/2001JD900171
  // dust_conc_sfc is updated by extract_atm_return_fields each step.
  // At step 0 dust_conc_sfc = 0 so feedback has no effect initially.
  if (m_params.loading_feedback_coeff > 0.0 && dust_conc_sfc) {
    for (amrex::MFIter mfi(*dust_ustar_t, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto ust  = dust_ustar_t->array(mfi);
        auto conc = dust_conc_sfc->const_array(mfi);
        const amrex::Real alpha = m_params.loading_feedback_coeff;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            ust(i,j,k) *= (1.0 + alpha * amrex::max(conc(i,j,k), 0.0));
        });
    }
    if (m_params.dust_debug) {
        amrex::Real ust_max = dust_ustar_t->max(0);
        amrex::Print() << "[DUST DEBUG] Phase 13: loading feedback applied"
                       << " ustar_t_max=" << ust_max << " m/s\n";
    }
  }

  // Phase 13: Apply Fecan et al. (1999) dynamic moisture inhibition
  // f_moist = sqrt(1 + a_f * max(w - w_prime, 0))
  // a_f = 1.21, w_prime = 0.003 (residual moisture threshold)
  // When use_dynamic_moisture = false (default, all current tests):
  //   use static moisture_flag from Phase 3 map (existing code path, unchanged).
  // When use_dynamic_moisture = true and Q1fx3 non-null:
  //   w = dust_surf_moist / (Lv * rho_a) where Lv = 2.501e6 J/kg.
  //   When Q1fx3 is null (no moisture scheme), dust_surf_moist = 0,
  //   so w = 0 and f_moist = 1.0 (no inhibition). This is correct for dry conditions.
  // Reference: Fecan et al. (1999), https://doi.org/10.1007/s00585-999-0149-7
  if (m_params.use_dynamic_moisture && dust_surf_moist) {
    constexpr amrex::Real Lv      = 2.501e6;
    constexpr amrex::Real rho_a   = 1.225;
    constexpr amrex::Real a_f     = 1.21;
    constexpr amrex::Real w_prime = 0.003;
    for (amrex::MFIter mfi(*dust_ustar_t, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto ust   = dust_ustar_t->array(mfi);
        auto qflux = dust_surf_moist->const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            amrex::Real w       = amrex::max(qflux(i,j,k), 0.0) / (Lv * rho_a);
            amrex::Real excess  = amrex::max(w - w_prime, 0.0);
            amrex::Real f_moist = std::sqrt(1.0 + a_f * excess);
            ust(i,j,k)  *= f_moist;
        });
    }
    if (m_params.dust_debug) {
        amrex::Real ust_max = dust_ustar_t->max(0);
        amrex::Print() << "[DUST DEBUG] Phase 13: dynamic moisture inhibition applied"
                       << " ustar_t_max=" << ust_max << " m/s\n";
    }
  }

  // Phase 6: Compute emission flux
  compute_dust_emission_flux(
    *dust_emission_flux, *dust_ustar_t, *dust_ustar_in, *dust_silt_fraction,
    m_params.n_size_bins, m_params.rho_air);

  if (m_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 6: emission_flux bin0 at step="
                   << m_step << " max=" << dust_emission_flux->max(0)
                   << " sum=" << dust_emission_flux->sum(0) << " [kg/m^2/s]\n";
  }

  // Phase 7: Blast schedule
#ifdef ERF_USE_DUST
  if (m_has_blast_schedule && dt > 0.0) {
    apply_blast_schedule(
      *dust_emission_flux, m_dg, m_blast_schedule, m_time, m_time - dt, dt,
      m_params.n_size_bins, m_params.blast_reactivity);
    if (m_params.dust_debug)
      amrex::Print()
        << "[DUST DEBUG] Phase 7: emission_flux bin0 after blast step="
        << m_step << " max=" << dust_emission_flux->max(0) << " [kg/m^2/s]\n";
  }

  // Phase 17: Compute PM2.5/PM10 NAAQS diagnostics
  // Called after extract_atm_return_fields has updated dust_conc_sfc.
  // Uses m_time - dt as the time at which fields were extracted (beginning of this step).
  compute_naaqs_diagnostics(dt, m_time - dt, m_step);

  // Phase 18: Compute MSHA worker exposure tracking
  // Called after compute_naaqs_diagnostics.
  compute_msha_exposure(dt, m_time - dt, m_step);

+#if defined(ERF_USE_PARTICLES)
+  // Phase 19: Release and advance Lagrangian super-particles
+  if (m_params.enable_particles && have_atm && xvel_mf && yvel_mf) {
+    // Note: advance_particles requires xvel/yvel/zvel as face-staggered fields.
+    // The xvel_mf and yvel_mf are passed as MultiFab arrays (u_mac, v_mac).
+    // zvel is extracted from velocity solver or computed from divergence.
+    // For simplicity in Phase 19, pass xvel_mf directly as u component.
+    // TODO: Extract zvel properly when available in 3D solver.
+    const amrex::MultiFab* zvel_ptr = nullptr;  // Phase 19 stub: no vertical velocity yet
+    // Actual call would be:
+    // advance_particles(*xvel_mf, *yvel_mf, *zvel_ptr, m_geom_atm, dt, m_step);
+    // For now, skip until zvel is properly passed from ERF_Advance.cpp
+    if (dust_params.dust_debug) {
+      amrex::Print() << "[DUST DEBUG] Phase 19: particle advance skipped (zvel not passed yet)\n";
+    }
+  }
+#endif
#endif
}
#ifdef ERF_USE_DUST

void
DustLayer::apply_to_cc_source(
  amrex::MultiFab& cc_source,
  const amrex::MultiFab& z_phys_cc,
  const amrex::Geometry& geom_atm)
{
  if (!dust_flux_atm)
    return;
  if (m_params.atm_feedback <= 0.0)
    return;

  // Phase 11: Support multi-bin injection when transport_bins_separately = true
  int n_active = m_params.transport_bins_separately ? m_params.n_size_bins : 1;

  for (int b = 0; b < n_active; ++b) {
    // For transport_bins_separately=false (default): sum all bins into slot 0.
    // For transport_bins_separately=true: inject bin b into slot b.
    dust_flux_atm->setVal(0.0);
    if (m_params.transport_bins_separately) {
      amrex::MultiFab::Copy(*dust_flux_atm, *dust_emission_flux, b, 0, 1,
                            amrex::IntVect(0));
    } else {
      // Sum all bins into slot 0 (only on first iteration).
      if (b == 0) {
        for (int bb = 0; bb < m_params.n_size_bins; ++bb) {
          amrex::MultiFab::Add(*dust_flux_atm, *dust_emission_flux, bb, 0, 1,
                               amrex::IntVect(0));
        }
      } else {
        break; // Only one pass needed when not transport_bins_separately
      }
    }

    // Coarsen from dust grid to atm grid.
    coarsen_dust_flux_to_atm(
      *dust_flux_atm, *dust_emission_flux, m_dg.geom, geom_atm, m_dg.grid_ratio);
    // Note: when grid_ratio = 1 the dust grid == atm 2D grid; coarsening is
    // a direct copy. When grid_ratio > 1, average_down reduces to atm resolution.

    // Inject into cc_source.
    apply_dust_tendency_to_cc_source(
      cc_source, *dust_flux_atm, z_phys_cc, geom_atm,
      m_dust_scalar_comp + b, m_params.atm_feedback, m_params.dust_debug);

    if (!m_params.transport_bins_separately) break; // only one pass needed
  }

  if (m_params.dust_debug) {
    amrex::Real F_max = dust_flux_atm->max(0);
    amrex::Print() << "[DUST DEBUG] Phase 10: step=" << m_step
                   << " F_dust_atm_max=" << F_max
                   << " kg/m^2/s  n_active_bins=" << n_active << "\n";
  }
}

void
DustLayer::apply_settling_to_cc_source(
    amrex::MultiFab& cc_source,
    const amrex::MultiFab& S_old,
    const amrex::MultiFab& z_phys_cc,
    const amrex::Geometry& geom_atm)
{
    if (m_params.bin_diameters.empty()) return;

    int n_active = m_params.transport_bins_separately
                 ? m_params.n_size_bins : 1;

    for (int b = 0; b < n_active; ++b) {
        // Clamp bin index to available diameters.
        int d_idx = (b < (int)m_params.bin_diameters.size())
                  ? b : (int)m_params.bin_diameters.size() - 1;
        amrex::Real d_m  = m_params.bin_diameters[d_idx];
        amrex::Real rhop = m_params.particle_density;
        int comp = m_dust_scalar_comp + b;

        apply_dust_settling_to_cc_source(cc_source, S_old, z_phys_cc,
                                         geom_atm, d_m, rhop, comp,
                                         m_params.dust_debug);
    }

    if (m_params.dust_debug) {
        amrex::Print() << "[DUST DEBUG] Phase 11: step=" << m_step
                       << " settling applied n_active=" << n_active
                       << " (erf.transport_scalar=true required)\n";
    }
}

void
DustLayer::apply_deposition_bc(
   amrex::MultiFab& cc_source, const amrex::MultiFab& S_old,
   const amrex::MultiFab& z_phys_cc, const amrex::Geometry& geom_atm,
   amrex::Real dt)
{
   if (!dep_flux_atm || !dust_ustar_in) return;
   if (m_params.atm_feedback <= 0.0) return;

   int n_active = m_params.transport_bins_separately
                ? m_params.n_size_bins : 1;

   for (int b = 0; b < n_active; ++b) {
       int d_idx = (b < (int)m_params.bin_diameters.size())
                 ? b : (int)m_params.bin_diameters.size()-1;
       amrex::Real d_m  = m_params.bin_diameters[d_idx];
       amrex::Real rhop = m_params.particle_density;
       amrex::Real E_0  = m_params.deposition_E0;
       int  comp = m_dust_scalar_comp + b;

       apply_dust_deposition_bc(cc_source, *dep_flux_atm,
                                 S_old, *dust_ustar_in,
                                 z_phys_cc, geom_atm,
                                 d_m, rhop, E_0, comp,
                                 m_params.dust_debug);

       apply_deposition_to_dust_grid(*dust_deposition_rate,
                                      *dep_flux_atm,
                                      m_dg.grid_ratio, dt);
   }

   if (m_params.dust_debug) {
       amrex::Real dep_sum = dust_deposition_rate->sum(0);
       amrex::Print() << "[DUST DEBUG] Phase 12: step=" << m_step
                      << " deposition_rate_sum=" << dep_sum
                      << " kg/m^2 (accumulated total)\n";
   }
}

void
DustLayer::extract_atm_return_fields(
    const amrex::MultiFab& S_new_cons,
    const amrex::MultiFab* Q1fx3,
    const amrex::Geometry& geom_atm)
{
    // Field A: near-surface dust concentration.
    fill_dust_conc_from_atm(*dust_conc_sfc, S_new_cons,
                             m_dust_scalar_comp, geom_atm, m_dg.grid_ratio);

    // Field B: surface moisture flux. Null-safe: returns immediately when
    // Q1fx3 == nullptr (moisture_type == None in current dust tests).
    const amrex::MultiFab* q1fx3_ptr = m_params.use_dynamic_moisture ? Q1fx3 : nullptr;
    fill_dust_moist_from_atm(*dust_surf_moist, q1fx3_ptr,
                              geom_atm, m_dg.grid_ratio);

    if (m_params.dust_debug) {
        amrex::Real conc_max  = dust_conc_sfc->max(0);
        amrex::Real conc_sum  = dust_conc_sfc->sum(0);
        amrex::Real moist_max = dust_surf_moist->max(0);
        amrex::Real dep_total = dust_deposition_rate ? dust_deposition_rate->sum(0) : 0.0;
        bool q1fx3_active = (Q1fx3 != nullptr) && m_params.use_dynamic_moisture;
        amrex::Print() << "[DUST DEBUG] Phase 13: step=" << m_step
                       << " conc_sfc_max=" << conc_max
                       << " kg/m^3  conc_sfc_sum=" << conc_sum
                       << "\n[DUST DEBUG] Phase 13:"
                       << " moist_flux_max=" << moist_max
                       << " kg/m^2/s  moisture_path_active=" << q1fx3_active
                       << "\n[DUST DEBUG] Phase 13:"
                       << " dep_total=" << dep_total
                       << " kg/m^2 (Phase 12 accumulator)\n";
        // Phase 14: MRF diffusion active
        amrex::Print() << "[DUST DEBUG] Phase 14: MRF diffusion active"
                       << " (erf.transport_scalar=true, EddyDiff::Scalar_v"
                       << " set by ComputeDiffusivityMRF)\n"
                       << "[DUST DEBUG] Phase 14: gamma_dust=0"
                       << " (no countergradient term for dust scalar)\n";
    }
}

void
DustLayer::compute_naaqs_diagnostics(amrex::Real dt, amrex::Real cur_time, int nstep)
{
    if (!dust_conc_sfc) return;

    int n_active = m_params.transport_bins_separately
                 ? m_params.n_size_bins : 1;

    // Task A: compute instantaneous PM2.5 and PM10.
    compute_pm_concentrations(*dust_pm25, *dust_pm10,
                               *dust_conc_sfc,
                               m_params.bin_diameters, n_active);

    // Task B: update 24-hour running averages.
    update_running_average(*dust_pm25_24h, *dust_pm25, dt, 86400.0, nstep);
    update_running_average(*dust_pm10_24h, *dust_pm10, dt, 86400.0, nstep);

    // Task C: exceedance flags.
    compute_exceedance_flag(*dust_pm25_exceed, *dust_pm25_24h,
                             DustPMConst::PM25_24H_NAAQS);
    compute_exceedance_flag(*dust_pm10_exceed, *dust_pm10_24h,
                             DustPMConst::PM10_24H_NAAQS);

    // Task D: CSV output.
    append_naaqs_stats(nstep, cur_time, m_params.dust_naaqs_file,
                       *dust_pm25, *dust_pm25_24h,
                       *dust_pm10, *dust_pm10_24h,
                       *dust_pm25_exceed, *dust_pm10_exceed);

    if (m_params.dust_debug) {
        amrex::Real pm25_max    = dust_pm25->max(0);
        amrex::Real pm25_24h_mx = dust_pm25_24h->max(0);
        amrex::Real pm10_max    = dust_pm10->max(0);
        amrex::Real n_ex25      = dust_pm25_exceed->sum(0);
        amrex::Real n_ex10      = dust_pm10_exceed->sum(0);
        amrex::Print() << "[DUST DEBUG] Phase 17: step=" << nstep
                       << " PM25_max=" << pm25_max << " ug/m^3"
                       << " PM25_24h_max=" << pm25_24h_mx << " ug/m^3"
                       << " PM10_max=" << pm10_max << " ug/m^3"
                       << " PM25_exceed_cells=" << (long)n_ex25
                       << " PM10_exceed_cells=" << (long)n_ex10 << "\n";
    }
}

void
DustLayer::compute_msha_exposure(amrex::Real dt, amrex::Real cur_time, int nstep)
{
    if (!dust_pm10) return;

    using namespace amrex;

    // Task A: Update dose and TWA from PM10
    update_msha_dose(*dust_msha_dose, *dust_msha_twa, *dust_pm10, dt);

    // Task B: Compute exceedance flag
    compute_msha_exceed(*dust_msha_exceed, *dust_msha_twa, m_params.msha_pel_mg_m3);

    // Task D: Shift reset check
    Real sd = m_params.msha_shift_duration_s;
    if (sd > 0.0 && std::floor(cur_time / sd) > std::floor((cur_time - dt) / sd)) {
        // Shift boundary crossed: copy TWA to shift_twa, write summary, reset dose
        MultiFab::Copy(*dust_msha_shift_twa, *dust_msha_twa, 0, 0, 1, 0);
        write_msha_shift_summary(++m_msha_shift_count, cur_time,
                                 m_params.msha_shift_file,
                                 *dust_msha_twa, *dust_msha_exceed);
        dust_msha_dose->setVal(0.0);
        if (m_params.dust_debug) {
            amrex::Print() << "[DUST DEBUG] Phase 18: shift " << m_msha_shift_count
                           << " ended t=" << cur_time << " s  dose reset\n";
        }
    }

    // Task D: Write per-step CSV
    append_msha_stats(nstep, cur_time, m_params.msha_exposure_file,
                      *dust_msha_twa, *dust_msha_exceed, *dust_msha_dose);

    // Task C: Sample receptor points
    for (int r = 0; r < (int)m_params.msha_receptor_names.size(); ++r) {
        append_receptor_sample(nstep, cur_time,
            "msha_receptor_" + m_params.msha_receptor_names[r] + ".csv",
            m_params.msha_receptor_names[r],
            m_params.msha_receptor_x[r], m_params.msha_receptor_y[r],
            *dust_pm10, m_dg.geom);
    }

    // Debug output
    if (m_params.dust_debug) {
        Real tmax = dust_msha_twa->max(0);
        Real nex = dust_msha_exceed->sum(0);
        amrex::Print() << "[DUST DEBUG] Phase 18: step=" << nstep
                       << " TWA_max=" << tmax << " mg/m3  exceed=" << (long)nex
                       << " shift=" << m_msha_shift_count << "\n";
    }
}

void
DustLayer::write_output(int nstep, double cur_time, bool is_final)
{
    // Always write CSV diagnostics.
    append_dust_stats(nstep, cur_time,
                      m_params.dust_diag_file,
                      get_emission_flux(),
                      get_deposition_rate(),
                      get_ustar_in(),
                      get_conc_sfc());

    // Determine whether to write the plotfile.
    bool write_plt = false;
    if (m_params.dust_plot_int > 0)
        write_plt = (nstep % m_params.dust_plot_int == 0);
    if (is_final && nstep > m_last_dust_plot_step)
        write_plt = true;

    if (write_plt) {
        WriteDustPlotfile(m_params.dust_plot_prefix, *this, cur_time, nstep);
        m_last_dust_plot_step = nstep;

        if (m_params.dust_debug) {
            amrex::Print() << "[DUST DEBUG] Phase 16: plotfile written"
                           << " step=" << nstep
                           << " is_final=" << is_final << "\n";
        }
    }

    if (m_params.dust_debug) {
        amrex::Real em = dust_emission_flux  ? dust_emission_flux->sum(0)    : 0.0;
        amrex::Real dp = dust_deposition_rate? dust_deposition_rate->sum(0) : 0.0;
        amrex::Print() << "[DUST DEBUG] Phase 16: step=" << nstep
                       << " emission_sum=" << em << " kg/s"
                       << " dep_total=" << dp << " kg/m^2\n";
    }
}

+#if defined(ERF_USE_PARTICLES)
+void
+DustLayer::advance_particles(const amrex::MultiFab& xvel,
+                             const amrex::MultiFab& yvel,
+                             const amrex::MultiFab& zvel,
+                             const amrex::Geometry& geom_atm,
+                             amrex::Real dt, int nstep)
+{
+    if (!m_dust_pc || !m_params.enable_particles) return;
+    if (!dust_emission_flux || !dust_source_map) return;
+
+    int n_active = m_params.transport_bins_separately
+                 ? m_params.n_size_bins : 1;
+
+    // Release particles at each step if interval matches
+    if (nstep % m_params.particle_release_interval == 0) {
+        for (int b = 0; b < n_active; ++b) {
+            int di = (b < (int)m_params.bin_diameters.size())
+                   ? b : (int)m_params.bin_diameters.size()-1;
+            m_dust_pc->ReleaseParticles(*dust_emission_flux,
+                geom_atm, m_dg.geom, dt,
+                m_params.bin_diameters[di], m_params.particle_density);
+        }
+    }
+    
+    // Advance particles
+    m_dust_pc->AdvanceParticles(xvel, yvel, zvel,
+        *dust_source_map, geom_atm, m_dg.geom, dt);
+
+    if (m_params.dust_debug) {
+        long np = m_dust_pc->TotalNumberOfParticles(true, false);
+        amrex::Real sm = dust_source_map->sum(0);
+        amrex::Print() << "[DUST DEBUG] Phase 19: step=" << nstep
+                       << " n_particles=" << np
+                       << " source_map_sum=" << sm << " kg/m^2\n";
+    }
+}
+#endif

#endif // ERF_USE_DUST

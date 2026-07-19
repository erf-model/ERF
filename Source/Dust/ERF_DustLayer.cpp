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
#include <ERF_DustTerrainSlope.H>
#include <ERF_DustEmission.H>
#include <ERF_DustSuppression.H>
#include <ERF_DustWindExtract.H>
#include <ERF_DustAtmCoupling.H>
#include <ERF_DustMSHA.H>
#include <ERF_DustMSHAOutput.H>
#include <ERF_FireWindExtract.H>
#include <ERF.H>
#include <ERF_SurfaceLayer.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFabUtil.H>
#include <cmath>

void
DustLayer::initialize(
  const ERF& erf,
  const SurfaceLayer* surface_layer,
  const amrex::MultiFab& z_phys_nd_atm,
  const DustParams& dust_params)
{
  verify_dust_prerequisites(erf, surface_layer, dust_params);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 1: Verified all prerequisites\n";
  }

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

  m_params = dust_params;

  amrex::IntVect ng(1, 1, 0);
  dust_slopes        = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 2, ng);
  dust_curvature     = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, 0);
  dust_slopes->setVal(0.0);
  dust_curvature->setVal(0.0);
  compute_dust_terrain_slopes(
    *dust_slopes, z_phys_nd_atm, erf.Geom(0), m_dg, dust_params.terrain_file);
  dust_slopes->FillBoundary(m_dg.geom.periodicity());
  compute_terrain_curvature(*dust_curvature, *dust_slopes, m_dg.geom);

  dust_ustar_t       = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
  dust_soil_type     = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
  dust_silt_fraction = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
  dust_crust_index   = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
  dust_moisture_flag = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
  dust_suppression   = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
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

  amrex::Real g    = 9.81;
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

  dust_surf_moist = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, 1, amrex::IntVect(1, 1, 0));
  dust_surf_moist->setVal(0.0);

  dust_ustar_in = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, 1, amrex::IntVect(1, 1, 0));
  dust_ustar_in->setVal(dust_params.test_ustar);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 6: dust_ustar_in (placeholder) = "
                   << dust_params.test_ustar << " m/s\n";
  }

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
    amrex::Print() << "[DUST DEBUG] Phase 9: zref=" << dust_params.zref
                   << " m for wind extraction\n";
  }

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

  m_dust_scalar_comp = RhoScalar_comp + 1;

  if (dust_params.dust_debug) {
    amrex::Print()
      << "[DUST DEBUG] Phase 10: dust_flux_atm allocated on atm grid."
      << " dust_scalar_comp=" << m_dust_scalar_comp << "\n";
  }

  dust_deposition_rate = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_deposition_rate->setVal(0.0);

  dep_flux_atm = std::make_unique<amrex::MultiFab>(
    dust_flux_atm->boxArray(),
    dust_flux_atm->DistributionMap(),
    1, amrex::IntVect(1,1,0));
  dep_flux_atm->setVal(0.0);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 12: dust_deposition_rate"
                   << " and dep_flux_atm allocated\n";
  }

  dust_conc_sfc = std::make_unique<amrex::MultiFab>(
    m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_conc_sfc->setVal(0.0);
  dust_surf_moist->setVal(0.0);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 13: dust_conc_sfc and"
                   << " dust_surf_moist allocated on dust grid\n"
                   << "[DUST DEBUG] Phase 13: loading_feedback_coeff="
                   << dust_params.loading_feedback_coeff
                   << " use_dynamic_moisture=" << dust_params.use_dynamic_moisture
                   << "\n";
  }

  dust_pm25      = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_pm10      = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_pm25_24h  = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_pm10_24h  = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
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

  dust_msha_dose      = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_msha_twa       = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_msha_exceed    = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_msha_shift_twa = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  for (auto* mf : {dust_msha_dose.get(), dust_msha_twa.get(),
                   dust_msha_exceed.get(), dust_msha_shift_twa.get()})
    mf->setVal(0.0);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 18: MSHA exposure MultiFabs allocated\n"
                   << "[DUST DEBUG] Phase 18: PEL=" << dust_params.msha_pel_mg_m3
                   << " mg/m3 (30 CFR 56.5001)\n"
                   << "[DUST DEBUG] Phase 18: shift_duration="
                   << dust_params.msha_shift_duration_s << " s\n"
                   << "[DUST DEBUG] Phase 18: n_receptors="
                   << dust_params.msha_receptor_names.size() << "\n";
    for (int r = 0; r < (int)dust_params.msha_receptor_names.size(); ++r) {
      amrex::Print() << "[DUST DEBUG] Phase 18: receptor " << r << ": "
                     << dust_params.msha_receptor_names[r] << " ("
                     << dust_params.msha_receptor_x[r] << ", "
                     << dust_params.msha_receptor_y[r] << ")\n";
    }
  }

#if defined(ERF_USE_PARTICLES)
  if (dust_params.enable_particles) {
    m_geom_atm = erf.Geom(0);
    m_dust_pc = std::make_unique<ERFDustPC>(
        m_geom_atm, erf.DistributionMap(0), erf.boxArray(0));
    dust_source_map = std::make_unique<amrex::MultiFab>(
        m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
    dust_source_map->setVal(0.0);

    if (dust_params.dust_debug) {
      amrex::Print() << "[DUST DEBUG] Phase 19: ERFDustPC initialized"
                     << " on dust grid " << m_dg.ba.size() << " boxes\n";
    }
  }
#endif

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

  // Phase 22: load haul road vehicle schedule.
  load_road_schedule(dust_params.road_schedule_file, m_road_schedule);
  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 22: road_schedule_file="
                   << dust_params.road_schedule_file
                   << " n_roads=" << m_road_schedule.roads.size() << "\n";
    for (int r = 0; r < (int)m_road_schedule.roads.size(); ++r) {
        const auto& ev = m_road_schedule.roads[r];
        amrex::Print() << "  road " << r << " name=" << ev.name
                       << " bbox=[" << ev.x_lo << "," << ev.y_lo << ","
                       << ev.x_hi << "," << ev.y_hi << "]"
                       << " W=" << ev.vehicle_weight_t << " t"
                       << " silt=" << ev.silt_pct << "%"
                       << " vmt=" << ev.vmt_per_h << "/h"
                       << " t=[" << ev.start_time_s << ","
                       << ev.end_time_s << "] s\n";
    }
  }
#endif

  recompute_dust_ustar_t(
    *dust_ustar_t, *dust_ustar_base, *dust_crust_index, *dust_efflor,
    *dust_surf_moist, *dust_suppression, dust_params.alpha_crust,
    dust_params.alpha_efflor, dust_slopes.get());

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

  amrex::Box dust_domain2 = m_dg.ba.minimalBox();
  int dust_nx2 = dust_domain2.length(0);
  int dust_ny2 = dust_domain2.length(1);
  amrex::Print() << "[DUST] DustLayer initialized: grid_ratio="
                 << m_dg.grid_ratio << ", dust cells=" << dust_nx2 << "x"
                 << dust_ny2 << "x1, "
                 << "n_size_bins=" << dust_params.n_size_bins << ", "
                 << "z0_dust=" << dust_params.z0_dust << " m\n";

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] PHREEQC configuration: "
                   << "phreeqc_output_file=\""
                   << dust_params.phreeqc_output_file << "\", "
                   << "phreeqc_update_interval_s="
                   << dust_params.phreeqc_update_interval_s << " s\n";
  }

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 11: bin_diameters [um]:";
    for (auto d : dust_params.bin_diameters)
      amrex::Print() << " " << d * 1e6;
    amrex::Print() << "\n"
                   << "[DUST DEBUG] Phase 11: transport_bins_separately="
                   << dust_params.transport_bins_separately << "\n";
  }

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

  write_dust_stats_header(m_params.dust_diag_file);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 16: dust_plot_int="
                   << dust_params.dust_plot_int
                   << " prefix=" << dust_params.dust_plot_prefix
                   << " diag_file=" << dust_params.dust_diag_file << "\n";
  }

  // Phase 20: allocate dust_site_id and populate from bounding boxes.
  dust_site_id = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_site_id->setVal(0.0);

  int n_sites = (int)dust_params.site_names.size();
  if (n_sites > 0) {
    populate_dust_site_id(*dust_site_id, m_dg.geom,
                          dust_params.site_x_lo, dust_params.site_y_lo,
                          dust_params.site_x_hi, dust_params.site_y_hi);

    auto counts = count_site_cells(*dust_site_id, n_sites);

    if (dust_params.dust_debug) {
      amrex::Print() << "[DUST DEBUG] Phase 20: " << n_sites
                     << " mine sites registered\n";
      amrex::Print() << "[DUST DEBUG] Phase 20: unassigned cells="
                     << counts[0] << "\n";
      for (int s = 0; s < n_sites; ++s) {
        amrex::Print() << "[DUST DEBUG] Phase 20: site " << (s+1)
                       << " name=" << dust_params.site_names[s]
                       << " file=" << dust_params.site_phreeqc_files[s]
                       << " bbox=["
                       << dust_params.site_x_lo[s] << ","
                       << dust_params.site_y_lo[s] << ","
                       << dust_params.site_x_hi[s] << ","
                       << dust_params.site_y_hi[s] << "]"
                       << " cells=" << counts[s+1] << "\n";
      }
    }

    m_site_ustar_t_factors.resize(n_sites, 1.0);
    for (int s = 0; s < n_sites; ++s) {
      if (!dust_params.site_phreeqc_files[s].empty()) {
        if (dust_params.dust_debug) {
          amrex::Print() << "[DUST DEBUG] Phase 20: site " << (s+1)
                         << " PHREEQC file: "
                         << dust_params.site_phreeqc_files[s] << "\n";
        }
      }
    }
  } else {
    if (dust_params.dust_debug)
      amrex::Print() << "[DUST DEBUG] Phase 20: no site bounding boxes"
                     << " defined; single-site mode (Phase 4 table)\n";
  }

  // Phase 21: PHREEQC deposition feedback file writer
  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 21: phreeqc_feedback_interval_s="
                   << dust_params.phreeqc_feedback_interval_s
                   << " feedback_file=" << dust_params.phreeqc_feedback_file
                   << " site_summary_file="
                   << dust_params.phreeqc_site_summary_file << "\n";
  }

  // Phase 23: critical material flux MultiFab.
  dust_cm_flux = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
  dust_cm_flux->setVal(0.0);

  if (dust_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 23: cm_fractions=";
    for (auto f : dust_params.cm_fractions)
      amrex::Print() << f << " ";
    amrex::Print() << "\n"
                   << "[DUST DEBUG] Phase 23: cm_budget_file="
                   << dust_params.cm_budget_file << "\n";
    if (dust_params.cm_fractions.empty())
      amrex::Print() << "[DUST DEBUG] Phase 23: cm_fractions empty,"
                     << " CM tracking disabled\n";
  }
} // end initialize()

void
DustLayer::advance(
  amrex::Real dt,
  const DustParams& dust_params,
  SurfaceLayer* surface_layer,
  const amrex::MultiFab* xvel_mf,
  const amrex::MultiFab* yvel_mf,
  const amrex::MultiFab* zvel_mf,
  const amrex::MultiFab* z_phys_cc_mf,
  const amrex::Geometry* geom_atm,
  int nz)
{
  ++m_step;
  m_time += dt;

  if (m_params.dust_debug) {
    amrex::Real cs  = dust_conc_sfc  ? dust_conc_sfc->max(0) * 1e9 : 0.0;
    amrex::Real p10 = dust_pm10      ? dust_pm10->max(0) : 0.0;
    amrex::Real twa = dust_msha_twa  ? dust_msha_twa->max(0) : 0.0;
    amrex::Real site_max = dust_site_id ? dust_site_id->max(0) : 0.0;
    amrex::Real ef_max = dust_emission_flux ? dust_emission_flux->max(0) : 0.0;
    int  n_roads_active = 0;
    if (m_road_schedule.loaded) {
        for (const auto& ev : m_road_schedule.roads) {
            bool active = (m_time >= ev.start_time_s) &&
                          (ev.end_time_s < 0.0 || m_time <= ev.end_time_s);
            if (active) ++n_roads_active;
        }
    }
    amrex::Print() << "[DUST DEBUG] advance: step=" << m_step
                   << " emission_flux_max=" << ef_max << " kg/m^2/s"
                   << " active_roads=" << n_roads_active
                   << " conc_sfc=" << cs << " ug/m3"
                   << " PM10=" << p10 << " ug/m3"
                   << " MSHA_TWA=" << twa << " mg/m3"
                   << " site_id_max=" << (int)site_max
                   << " n_sites=" << (int)m_params.site_names.size() << "\n";
  }

  bool have_atm =
    (surface_layer && xvel_mf && yvel_mf && z_phys_cc_mf && nz > 0);

  if (have_atm) {
    if (surface_layer->get_u_star(0))
      fill_dust_ustar_from_surface_layer(
        *dust_ustar_in, *surface_layer->get_u_star(0), m_dg);
    fill_dust_wind_from_interpolation(
      *dust_wind_ref, *xvel_mf, *yvel_mf, *z_phys_cc_mf, m_dg, m_params.zref, nz);
    if (m_params.use_terrain_wind && dust_slopes && dust_curvature) {
      apply_farsite_terrain_wind(
        *dust_wind_ref, *dust_slopes, *dust_curvature, m_params.k_ridge,
        m_params.k_shelter, m_params.k_valley, m_params.k_deflect);
      // Recompute dust_ustar_in from the terrain-corrected wind using log-profile.
      // This ensures the FARSITE terrain wind corrections (ridge speed-up, lee sheltering,
      // valley channeling) are properly reflected in the friction velocity used for
      // dust emission calculations.
      compute_dust_ustar_from_wind(
        *dust_ustar_in, *dust_wind_ref, m_params.zref, m_params.z0_dust);
    }
    if (surface_layer->get_t_surf(0))
      fill_dust_scalar_from_atm(*dust_tsfc, *surface_layer->get_t_surf(0), m_dg);
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

  if (m_params.dust_debug && have_atm) {
    amrex::Print() << "[DUST DEBUG] Phase 9: T_sfc_max=" << dust_tsfc->max(0)
                   << " K  PBLH_max=" << dust_pblh->max(0) << " m\n";
  }

  amrex::Real T_sfc = have_atm ? dust_tsfc->max(0) : m_params.test_surf_temp_K;
  amrex::Real u_10m = have_atm ? dust_wind_ref->max(0) : m_params.test_wind_speed;

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

  amrex::MultiFab::Copy(
    *dust_surf_moist, *dust_moisture_flag, 0, 0, 1, amrex::IntVect(1, 1, 0));

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

  recompute_dust_ustar_t(
    *dust_ustar_t, *dust_ustar_base, *dust_crust_index, *dust_efflor,
    *dust_surf_moist, *dust_suppression, m_params.alpha_crust,
    m_params.alpha_efflor, dust_slopes.get());

  if (m_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 7: crust_index before u*_t computation at step=" << m_step
                   << " min=" << dust_crust_index->min(0)
                   << " max=" << dust_crust_index->max(0) << "\n";
    amrex::Print() << "[DUST DEBUG] Phase 5: u*_t at step=" << m_step
                   << " min=" << dust_ustar_t->min(0)
                   << " max=" << dust_ustar_t->max(0) << " [m/s]\n";
  }

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
        ust(i,j,k) *= f_moist;
      });
    }
    if (m_params.dust_debug) {
      amrex::Real ust_max = dust_ustar_t->max(0);
      amrex::Print() << "[DUST DEBUG] Phase 13: dynamic moisture inhibition applied"
                     << " ustar_t_max=" << ust_max << " m/s\n";
    }
  }

  compute_dust_emission_flux(
    *dust_emission_flux, *dust_ustar_t, *dust_ustar_in, *dust_silt_fraction,
    m_params.n_size_bins, m_params.rho_air);

  if (m_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 6: emission_flux bin0 at step="
                   << m_step << " max=" << dust_emission_flux->max(0)
                   << " sum=" << dust_emission_flux->sum(0) << " [kg/m^2/s]\n";
  }

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

  // Phase 22: apply haul road vehicle resuspension (additive to emission flux).
  apply_road_schedule(*dust_emission_flux, m_dg.geom, m_road_schedule,
                      m_time, dt, m_params.dust_debug,
                      m_params.road_diag_file, m_step);

  // Phase 23: compute critical material flux and write budget.
  if (!m_params.cm_fractions.empty() && dust_emission_flux && dust_cm_flux) {
    int n_active = m_params.transport_bins_separately
                 ? m_params.n_size_bins : 1;
    compute_cm_flux(*dust_cm_flux, *dust_emission_flux,
                    m_params.cm_fractions, n_active);
    append_cm_budget(m_params.cm_budget_file,
                     *dust_cm_flux,
                     dust_site_id.get(),
                     m_dg.geom,
                     m_time, m_step,
                     m_params.site_names);
    if (m_params.dust_debug) {
      amrex::Real cm_max = dust_cm_flux->max(0);
      amrex::Real cm_sum = dust_cm_flux->sum(0);
      amrex::Print() << "[DUST DEBUG] Phase 23: step=" << m_step
                     << " cm_flux_max=" << cm_max << " kg_CM/m^2/s"
                     << " cm_flux_sum=" << cm_sum << " kg_CM/m^2/s\n";
    }
  }

  compute_naaqs_diagnostics(dt, m_time - dt, m_step);
  compute_msha_exposure(dt, m_time - dt, m_step);
#endif

#if defined(ERF_USE_PARTICLES)
  // Phase 19 diagnostic: check why block may not execute
  if (m_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 19 check: step=" << m_step
                   << " enable_particles=" << m_params.enable_particles
                   << " m_dust_pc=" << (m_dust_pc ? "valid" : "null")
                   << " xvel_mf=" << (xvel_mf ? "valid" : "null")
                   << " yvel_mf=" << (yvel_mf ? "valid" : "null")
                   << " zvel_mf=" << (zvel_mf ? "valid" : "null")
                   << " geom_atm=" << (geom_atm ? "valid" : "null")
                   << " interval_check=" << (m_step % m_params.particle_release_interval)
                   << " emission_max=" << dust_emission_flux->max(0) << "\n";
  }
  if (m_params.enable_particles && m_dust_pc && xvel_mf && yvel_mf && zvel_mf
      && geom_atm && (m_step % m_params.particle_release_interval == 0)) {
    amrex::Real d_m   = m_params.bin_diameters.empty() ? 7.0e-6 : m_params.bin_diameters[0];
    amrex::Real rho_p = m_params.particle_density;
    m_dust_pc->ReleaseParticles(*dust_emission_flux, *geom_atm, m_dg.geom,
                                 dt, d_m, rho_p);
    if (dust_source_map) {
      m_dust_pc->AdvanceParticles(*xvel_mf, *yvel_mf, *zvel_mf,
                                   *dust_source_map, *geom_atm, m_dg.geom, dt);
    }
    if (m_params.dust_debug) {
      long np = m_dust_pc->TotalNumberOfParticles();
      amrex::Real sm = dust_source_map ? dust_source_map->sum(0) : 0.0;
      amrex::Print() << "[DUST DEBUG] Phase 19: step=" << m_step
                     << " n_particles=" << np
                     << " source_map_sum=" << sm << " kg/m^2\n";
    }
  }
#endif

} // end advance()

#ifdef ERF_USE_DUST

void
DustLayer::apply_to_cc_source(
  amrex::MultiFab& cc_source,
  const amrex::MultiFab& z_phys_cc,
  const amrex::Geometry& geom_atm)
{
  if (!dust_flux_atm) return;
  if (m_params.atm_feedback <= 0.0) return;

  int n_active = m_params.transport_bins_separately ? m_params.n_size_bins : 1;

  // Temporary 1-component MultiFab on the dust grid for bin summing.
  // Needed when transport_bins_separately=false to sum all bins first
  // before coarsening. Allocated once per call (small 2D slab, cheap).
  amrex::MultiFab dust_bin_tmp(m_dg.ba, m_dg.dm, 1, amrex::IntVect(0));

  for (int b = 0; b < n_active; ++b) {
    dust_bin_tmp.setVal(0.0);

    if (m_params.transport_bins_separately) {
      // Extract bin b from dust_emission_flux into the temporary.
      amrex::MultiFab::Copy(dust_bin_tmp, *dust_emission_flux, b, 0, 1,
                            amrex::IntVect(0));
    } else {
      // Sum all bins into the temporary (b==0 only; loop breaks below).
      for (int bb = 0; bb < m_params.n_size_bins; ++bb) {
        amrex::MultiFab::Add(dust_bin_tmp, *dust_emission_flux, bb, 0, 1,
                             amrex::IntVect(0));
      }
    }

    // FIX for grid_ratio > 1: average_down directly from the dust grid
    // (dust_bin_tmp, m_dg.geom) to the atm 2D slab (dust_flux_atm, geom_atm_2d).
    // Previously this called MultiFab::Copy(*dust_flux_atm, *dust_emission_flux,...)
    // which fails when grid_ratio > 1 because the BoxArrays have different sizes.
    dust_flux_atm->setVal(0.0);
    {
      // Build a 2D geometry for dust_flux_atm (same physical domain, atm resolution).
      // geom_atm is 3D; we need the 2D slab version with z-extent = [0,0].
      amrex::Box atm_domain_2d = geom_atm.Domain();
      atm_domain_2d.setSmall(2, 0);
      atm_domain_2d.setBig(2, 0);
      amrex::RealBox prob_2d = geom_atm.ProbDomain();
      prob_2d.setHi(2, prob_2d.lo(2) + 1.0); // dummy 1 m z-extent
      amrex::Geometry geom_atm_2d(atm_domain_2d, prob_2d,
                                   amrex::CoordSys::cartesian,
                                   {false, false, false});

      amrex::average_down(dust_bin_tmp, *dust_flux_atm,
                          m_dg.geom, geom_atm_2d,
                          0, 1,
                          amrex::IntVect(m_dg.grid_ratio, m_dg.grid_ratio, 1));
    }

    apply_dust_tendency_to_cc_source(
      cc_source, *dust_flux_atm, z_phys_cc, geom_atm,
      m_dust_scalar_comp + b, m_params.atm_feedback, m_params.dust_debug);

    if (!m_params.transport_bins_separately) break;
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

    int n_active = m_params.transport_bins_separately ? m_params.n_size_bins : 1;

    for (int b = 0; b < n_active; ++b) {
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

   int n_active = m_params.transport_bins_separately ? m_params.n_size_bins : 1;

   for (int b = 0; b < n_active; ++b) {
       int d_idx = (b < (int)m_params.bin_diameters.size())
                 ? b : (int)m_params.bin_diameters.size()-1;
       amrex::Real d_m  = m_params.bin_diameters[d_idx];
       amrex::Real rhop = m_params.particle_density;
       amrex::Real E_0  = m_params.deposition_E0;
       int comp = m_dust_scalar_comp + b;

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
    fill_dust_conc_from_atm(*dust_conc_sfc, S_new_cons,
                             m_dust_scalar_comp, geom_atm, m_dg.grid_ratio);

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

    int n_active = m_params.transport_bins_separately ? m_params.n_size_bins : 1;

    compute_pm_concentrations(*dust_pm25, *dust_pm10,
                               *dust_conc_sfc,
                               m_params.bin_diameters, n_active);

    update_running_average(*dust_pm25_24h, *dust_pm25, dt, 86400.0, nstep);
    update_running_average(*dust_pm10_24h, *dust_pm10, dt, 86400.0, nstep);

    compute_exceedance_flag(*dust_pm25_exceed, *dust_pm25_24h,
                             DustPMConst::PM25_24H_NAAQS);
    compute_exceedance_flag(*dust_pm10_exceed, *dust_pm10_24h,
                             DustPMConst::PM10_24H_NAAQS);

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

    update_msha_dose(*dust_msha_dose, *dust_msha_twa, *dust_pm10, dt);
    compute_msha_exceed(*dust_msha_exceed, *dust_msha_twa, m_params.msha_pel_mg_m3);

    Real sd = m_params.msha_shift_duration_s;
    if (sd > 0.0 && cur_time > dt &&
        std::floor(cur_time / sd) > std::floor((cur_time - dt) / sd)) {
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

    append_msha_stats(nstep, cur_time, m_params.msha_exposure_file,
                      *dust_msha_twa, *dust_msha_exceed, *dust_msha_dose);

    for (int r = 0; r < (int)m_params.msha_receptor_names.size(); ++r) {
        append_receptor_sample(nstep, cur_time,
            "msha_receptor_" + m_params.msha_receptor_names[r] + ".csv",
            m_params.msha_receptor_names[r],
            m_params.msha_receptor_x[r], m_params.msha_receptor_y[r],
            *dust_pm10, m_dg.geom);
    }

    if (m_params.dust_debug) {
        Real tmax = dust_msha_twa->max(0);
        Real nex  = dust_msha_exceed->sum(0);
        amrex::Print() << "[DUST DEBUG] Phase 18: step=" << nstep
                       << " TWA_max=" << tmax << " mg/m3  exceed=" << (long)nex
                       << " shift=" << m_msha_shift_count << "\n";
    }
}

void
DustLayer::write_output(int nstep, double cur_time, bool is_final)
{
    append_dust_stats(nstep, cur_time,
                      m_params.dust_diag_file,
                      get_emission_flux(),
                      get_deposition_rate(),
                      get_ustar_in(),
                      get_conc_sfc());

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

    // Phase 21: PHREEQC deposition feedback file writer
    write_phreeqc_feedback(nstep, cur_time, is_final);

    if (m_params.dust_debug) {
        amrex::Real em = dust_emission_flux   ? dust_emission_flux->sum(0)   : 0.0;
        amrex::Real dp = dust_deposition_rate ? dust_deposition_rate->sum(0) : 0.0;
        amrex::Print() << "[DUST DEBUG] Phase 16: step=" << nstep
                       << " emission_sum=" << em << " kg/s"
                       << " dep_total=" << dp << " kg/m^2\n";
    }
}

void
DustLayer::write_phreeqc_feedback(int nstep, amrex::Real cur_time,
                                   bool is_final)
{
    if (!dust_deposition_rate) return;
    const amrex::Real interval = m_params.phreeqc_feedback_interval_s;
    if (interval <= 0.0 && !is_final) return;

    if (nstep == m_last_phreeqc_write_step) return;

    bool do_write = false;
    if (interval > 0.0 &&
        (cur_time - m_last_phreeqc_write_time) >= interval - 0.5*interval*1e-6)
        do_write = true;
    if (is_final && nstep > m_last_phreeqc_write_step)
        do_write = true;

    if (!do_write) return;

    m_last_phreeqc_write_step = nstep;
    m_last_phreeqc_write_time = cur_time;

    write_phreeqc_deposition_file(
        m_params.phreeqc_feedback_file,
        *dust_deposition_rate,
        m_dg.geom,
        cur_time, nstep);

    append_phreeqc_site_summary(
        m_params.phreeqc_site_summary_file,
        *dust_deposition_rate,
        dust_site_id.get(),
        m_dg.geom,
        cur_time,
        m_params.site_names);

    if (m_params.dust_debug) {
        amrex::Real dep_sum = dust_deposition_rate->sum(0);
        amrex::Print() << "[DUST DEBUG] Phase 21: PHREEQC feedback written"
                       << " step=" << nstep
                       << " time=" << cur_time
                       << " dep_sum=" << dep_sum << " kg/m^2"
                       << " is_final=" << is_final << "\n";
    }
}


#endif // ERF_USE_DUST
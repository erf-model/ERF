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

  // Determine scalar component index. Use RhoAdv_comp if available.
  // RhoAdv_comp is the first passive advected scalar in ERF_IndexDefines.H.
  m_dust_scalar_comp = RhoAdv_comp;

  if (dust_params.dust_debug) {
    amrex::Print()
      << "[DUST DEBUG] Phase 10: dust_flux_atm allocated on atm grid."
      << " dust_scalar_comp=" << m_dust_scalar_comp << "\n";
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
}

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

  // Sum emission flux over bins into dust_flux_atm.
  // Phase 11 will handle per-bin injection; here only total mass is used.
  dust_flux_atm->setVal(0.0);
  for (int b = 0; b < m_params.n_size_bins; ++b) {
    amrex::MultiFab::Add(
      *dust_flux_atm, *dust_emission_flux, b, 0, 1, amrex::IntVect(0));
  }

  // Coarsen from dust grid to atm grid.
  coarsen_dust_flux_to_atm(
    *dust_flux_atm, *dust_emission_flux, m_dg.geom, geom_atm, m_dg.grid_ratio);
  // Note: when grid_ratio = 1 the dust grid == atm 2D grid; coarsening is
  // a direct copy. When grid_ratio > 1, average_down reduces to atm resolution.

  // Inject into cc_source.
  apply_dust_tendency_to_cc_source(
    cc_source, *dust_flux_atm, z_phys_cc, geom_atm, m_dust_scalar_comp,
    m_params.atm_feedback, m_params.dust_debug);

  if (m_params.dust_debug) {
    amrex::Print() << "[DUST DEBUG] Phase 10: step=" << m_step
                   << " dust injected into cc_source comp="
                   << m_dust_scalar_comp << "\n";
  }
}

#endif // ERF_USE_DUST

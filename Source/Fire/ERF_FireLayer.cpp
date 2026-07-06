#include <ERF_FireLayer.H>
#include <ERF.H>
#include <ERF_SurfaceLayer.H>
#include <ERF_FirePrerequisites.H>
#include <ERF_FireGrid.H>
#include <ERF_FireWindExtract.H>
#include <ERF_TerrainSlope.H>

using namespace amrex;

void FireLayer::initialize(const ERF& erf,
                            const SurfaceLayer* surface_layer_ptr,
                            const MultiFab& z_phys_nd_atm,
                            const FireParams& fire_params)
{
    m_params = fire_params;
    verify_fire_prerequisites(erf, surface_layer_ptr, fire_params);
    m_fg = create_fire_grid(erf.boxArray(0), erf.DistributionMap(0),
                            erf.Geom(0), fire_params.grid_ratio);
    m_nz = erf.Geom(0).Domain().length(2);

    // Allocate MultiFabs (fire_phi needs 1 ghost cell for gradient stencils)
    fire_phi          = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 1);
    fire_wind_ref     = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_wind_eff     = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_slopes       = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 1);
    fire_curvature    = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_ros          = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_fuel_load    = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_fuel_mc      = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 3, 0);
    fire_mext         = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_heat_flux    = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_spread_vec   = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_arrival_time = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_disp_accum   = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
     
    // Phase 5: Heat flux and diagnostics fields
    fire_fireline_intensity = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_flame_length = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);

    fire_phi->setVal(1.0_rt);
    fire_wind_ref->setVal(0.0_rt);
    fire_wind_eff->setVal(0.0_rt);
    fire_slopes->setVal(0.0_rt);
    fire_curvature->setVal(0.0_rt);
    fire_ros->setVal(0.0_rt);
    fire_arrival_time->setVal(-1.0_rt);  // -1 = unburned
    fire_disp_accum->setVal(0.0_rt);     // zero accumulated displacement

    // Fuel load
    FuelModelParams fp = get_anderson_fuel_params(fire_params.fuel_model_id);
    fire_fuel_load->setVal((fp.w_d1+fp.w_d10+fp.w_d100+fp.w_lh+fp.w_lw)*4.88243);
    m_fuel_load_initial_kg_m2 = (fp.w_d1+fp.w_d10+fp.w_d100+fp.w_lh+fp.w_lw)*4.88243_rt;
    m_fuel_bed_depth_ft = fp.delta;

    // Fuel moisture
    fire_fuel_mc->setVal(0.0_rt);
    for (MFIter mfi(*fire_fuel_mc); mfi.isValid(); ++mfi) {
        Array4<Real> mc = fire_fuel_mc->array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            mc(iv,0) = fire_params.moisture_1hr;
            mc(iv,1) = fire_params.moisture_10hr;
            mc(iv,2) = fire_params.moisture_100hr;
        });
    }

    // Moisture of extinction
    Real dead_load = fp.w_d1 + fp.w_d10 + fp.w_d100;
    Real sigma_weighted = dead_load > 0.0_rt
        ? (fp.w_d1 * fp.sigma_d1) / dead_load : fp.sigma_d1;
    fire_mext->setVal(compute_moisture_of_extinction(sigma_weighted));

    // Terrain slopes and curvature
    compute_terrain_slopes(*fire_slopes, z_phys_nd_atm, erf.Geom(0), m_fg,
                           m_params.terrain_file_name);
    fire_slopes->FillBoundary(m_fg.geom.periodicity());
    compute_terrain_curvature(*fire_curvature, *fire_slopes, m_fg.geom);

    // Ignition
    m_ignition_x = fire_params.ignition_x;
    m_ignition_y = fire_params.ignition_y;
    m_ignition_r = fire_params.ignition_r;
    initialize_ignition(*fire_phi, m_fg.geom, m_ignition_x, m_ignition_y, m_ignition_r);
    fire_phi->FillBoundary(m_fg.geom.periodicity());

    // Mark ignition cells in arrival time (t=0)
    for (MFIter mfi(*fire_phi); mfi.isValid(); ++mfi) {
        Array4<const Real> phi_arr = fire_phi->const_array(mfi);
        Array4<Real>       at_arr  = fire_arrival_time->array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            if (phi_arr(iv) < 0.0_rt) at_arr(iv) = 0.0_rt;
        });
    }

    // FARSITE parameters
    m_fp.phi_threshold   = fire_params.farsite_phi_threshold;
    m_fp.coeff_a         = fire_params.farsite_coeff_a;
    m_fp.coeff_b         = fire_params.farsite_coeff_b;
    m_fp.coeff_c         = fire_params.farsite_coeff_c;
    m_fp.use_anderson_lw = fire_params.farsite_use_anderson_lw;
    m_fp.gaussian_sigma  = fire_params.farsite_gaussian_sigma;
    m_fp.cfl_fire        = fire_params.farsite_cfl_fire;

    // Rothermel coefficients
    m_rc = compute_rothermel_params(fp, fire_params.moisture_1hr,
                                    fire_params.moisture_10hr,
                                    fire_params.moisture_100hr);

    amrex::Print() << "[FIRE] FireLayer initialized: C=" << m_fg.C
                   << ", fuel_model=" << fire_params.fuel_model_id
                   << ", grid=" << m_fg.ba.size() << " boxes" << std::endl;

    if (m_params.fire_debug) {
        IntVect max_extent  = m_fg.geom.Domain().size();
        RealBox prob_domain = m_fg.geom.ProbDomain();
        Real dx_fire = (prob_domain.hi(0) - prob_domain.lo(0)) / max_extent[0];
        Real dy_fire = (prob_domain.hi(1) - prob_domain.lo(1)) / max_extent[1];
        amrex::Print() << "[FIRE DEBUG] Fire grid: dx=" << dx_fire
                       << " m, dy=" << dy_fire << " m, "
                       << max_extent[0] << " x " << max_extent[1]
                       << " cells, grid_ratio=" << m_fg.C << std::endl;
    }
}

void FireLayer::advance(Real time, Real dt, SurfaceLayer& surface_layer,
                        const MultiFab& xvel, const MultiFab& yvel,
                        const MultiFab& z_phys_cc,
                        const MultiFab& T_atm_k0, const MultiFab& RH_atm_k0)
{
    m_current_time = time;
    m_dt_atm       = dt;
    ++m_step;  // Increment step counter

    if (m_params.fire_debug)
        amrex::Print() << "[FIRE DEBUG] Starting fire advance step with dt=" << dt << std::endl;

    // 1. Wind extraction
    fill_fire_wind_from_interpolation(*fire_wind_ref, xvel, yvel, z_phys_cc,
                                      m_fg, m_params.wind_ref_ht, m_nz);
    if (m_params.fire_debug)
        amrex::Print() << "[FIRE DEBUG] Wind extraction completed. Max reference wind: "
                       << fire_wind_ref->max(0) << " m/s" << std::endl;

    // 2. Copy to effective wind
    MultiFab::Copy(*fire_wind_eff, *fire_wind_ref, 0, 0, 2, 0);

    // 3. WAF
    if (m_params.use_waf) {
        if (m_params.fire_debug)
            amrex::Print() << "[FIRE DEBUG] Applying Wind Adjustment Factor (formula: "
                           << m_params.waf_formula << ")" << std::endl;
        apply_waf_to_wind();
    }

    // 4. Terrain wind
    if (m_params.use_terrain_wind) {
        if (m_params.fire_debug)
            amrex::Print() << "[FIRE DEBUG] Applying FARSITE terrain wind corrections" << std::endl;
        apply_farsite_terrain_wind(*fire_wind_eff, *fire_slopes, *fire_curvature,
                                   m_params.k_ridge, m_params.k_shelter,
                                   m_params.k_valley, m_params.k_deflect);
    }

    if (m_params.fire_debug)
        amrex::Print() << "[FIRE DEBUG] Effective wind computed. Max effective wind: "
                       << fire_wind_eff->max(0) << " m/s" << std::endl;

    // 5. Fuel moisture
    if (m_params.fire_debug)
        amrex::Print() << "[FIRE DEBUG] Updating fuel moisture from atmospheric state" << std::endl;
    advance_fuel_moisture(dt, T_atm_k0, RH_atm_k0);
    if (m_params.fire_debug)
        amrex::Print() << "[FIRE DEBUG] Fuel moisture update completed. Max 1-hour moisture: "
                       << fire_fuel_mc->max(0) << std::endl;

    // 5b. Dynamic Rothermel update
    if (m_params.moisture_dynamic) {
        long num_cells   = fire_fuel_mc->boxArray().numPts();
        Real avg_mc_1hr   = (num_cells > 0) ? fire_fuel_mc->sum(0) / Real(num_cells) : m_params.moisture_1hr;
        Real avg_mc_10hr  = (num_cells > 0) ? fire_fuel_mc->sum(1) / Real(num_cells) : m_params.moisture_10hr;
        Real avg_mc_100hr = (num_cells > 0) ? fire_fuel_mc->sum(2) / Real(num_cells) : m_params.moisture_100hr;
        avg_mc_1hr   = amrex::max(0.01_rt, amrex::min(avg_mc_1hr,   0.40_rt));
        avg_mc_10hr  = amrex::max(0.01_rt, amrex::min(avg_mc_10hr,  0.40_rt));
        avg_mc_100hr = amrex::max(0.01_rt, amrex::min(avg_mc_100hr, 0.40_rt));
        FuelModelParams fp_cur = get_anderson_fuel_params(m_params.fuel_model_id);
        m_rc = compute_rothermel_params(fp_cur, avg_mc_1hr, avg_mc_10hr, avg_mc_100hr);
    }

    // 6. ROS
    compute_ros_field(*fire_ros, *fire_wind_eff, *fire_slopes, m_rc);
    if (m_params.fire_debug)
        amrex::Print() << "[FIRE DEBUG] Rate-of-spread computed. Max: " << fire_ros->max(0)
                       << " m/s, Mean: " << fire_ros->sum(0)/fire_ros->boxArray().numPts()
                       << " m/s" << std::endl;

    // 7. FARSITE subcycle
    // The subcycle resets phi to 0 then stamps the propagated front ring.
    // Displacement is accumulated across atmospheric steps in fire_disp_accum.
    // A cell is only stamped when accum magnitude >= 0.5 * dx_fire.
    int n_substeps = advance_fire_subcycle(*fire_phi, *fire_spread_vec,
                                           *fire_disp_accum,
                                           *fire_arrival_time,   // now mutable
                                           *fire_wind_eff, *fire_ros,
                                           m_fg.geom, dt, m_fp);

    if (m_params.fire_debug)
        amrex::Print() << "[FIRE DEBUG] Level-set propagation completed with "
                       << n_substeps << " fire subcycles" << std::endl;

    // 8. Re-apply burned interior from arrival time.
    // The FARSITE subcycle only stamps newly propagated front cells.
    // Restore all previously burned cells from arrival_time here.
    for (MFIter mfi(*fire_phi); mfi.isValid(); ++mfi) {
        Array4<Real>       p  = fire_phi->array(mfi);
        Array4<const Real> at = fire_arrival_time->const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            if (at(iv) >= 0.0_rt) p(iv) = -1.0_rt;
        });
    }

    // 9. Record arrival time for cells newly burned this step.
    Real cur_time = m_current_time + dt;
    for (MFIter mfi(*fire_phi); mfi.isValid(); ++mfi) {
        Array4<const Real> p  = fire_phi->const_array(mfi);
        Array4<Real>       at = fire_arrival_time->array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            if (p(iv) < 0.0_rt && at(iv) < 0.0_rt) at(iv) = cur_time;
        });
    }

    if (m_params.fire_debug) {
        long num_fire_cells = 0;
        for (MFIter mfi(*fire_phi); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<const Real> phi = fire_phi->array(mfi);
            for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k)
                for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
                    for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
                        if (phi(i,j,k) < 0.0_rt) num_fire_cells++;
        }
        amrex::ParallelDescriptor::ReduceLongSum(num_fire_cells);
        amrex::Print() << "[FIRE DEBUG] Number of active fire cells: " << num_fire_cells << std::endl;
    }

    fire_phi->FillBoundary(m_fg.geom.periodicity());

    Real max_ros  = fire_ros->max(0);
    Real mean_ros = fire_ros->sum(0) / fire_ros->boxArray().numPts();
    Real phi_min  = fire_phi->min(0);
    amrex::Print() << "[FIRE] t=" << m_current_time
                   << "  substeps=" << n_substeps
                   << "  phi_min=" << fire_phi->min(0)
                   << "  max_ROS=" << fire_ros->max(0) << " m/s"
                   << "  mean_ROS=" << fire_ros->sum(0)/fire_ros->boxArray().numPts()
                   << " m/s"
                   << "  Q_max=" << fire_heat_flux->max(0) << " W/m2"
                   << "  I_B_max=" << fire_fireline_intensity->max(0) << " kW/m"
                   << "  L_max=" << fire_flame_length->max(0) << " m"
                   << std::endl;
    
    // Write fire statistics to CSV if enabled
    if (m_params.write_fire_stats_csv) {
        static bool csv_header_written = false;
        if (!csv_header_written) {
            write_fire_stats_header(m_params.fire_stats_csv_file);
            csv_header_written = true;
        }
        append_fire_stats(*fire_phi, *fire_arrival_time, m_fg.geom,
                         m_step, m_current_time, m_params.fire_stats_csv_file,
                         fire_ros.get(), fire_heat_flux.get());
    }
}


void FireLayer::apply_waf_to_wind()
{
    Real waf = 0.4_rt;
    if      (m_params.waf_formula == "andrews")      waf = compute_waf_unsheltered(m_fuel_bed_depth_ft);
    else if (m_params.waf_formula == "behaviorplus")  waf = compute_waf_behaviorplus(m_fuel_bed_depth_ft);
    else amrex::Print() << "[FIRE WARNING] Unknown waf_formula='" << m_params.waf_formula
                        << "'. Using default WAF=" << waf << std::endl;
    fire_wind_eff->mult(waf, 0, 2, 0);
}

void FireLayer::advance_fuel_moisture(Real dt_s,
                                      const MultiFab& T_atm_k0,
                                      const MultiFab& RH_atm_k0)
{
    if (!m_params.moisture_dynamic) return;

    Real dt_hours     = dt_s / 3600.0_rt;
    Real precip_mm_hr = m_params.precip_rate_mm_hr;
    int  C            = m_fg.C;

    FuelModelParams fp = get_anderson_fuel_params(m_params.fuel_model_id);

    MultiFab T_fire (fire_fuel_mc->boxArray(), fire_fuel_mc->DistributionMap(), 1, 0);
    MultiFab RH_fire(fire_fuel_mc->boxArray(), fire_fuel_mc->DistributionMap(), 1, 0);

    for (MFIter mfi(T_fire, false); mfi.isValid(); ++mfi) {
        Array4<Real>       T_f  = T_fire.array(mfi);
        Array4<Real>       RH_f = RH_fire.array(mfi);
        Array4<const Real> T_a  = T_atm_k0.const_array(mfi);
        Array4<const Real> RH_a = RH_atm_k0.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv_f) {
            int i_a = iv_f[0] / C;
            int j_a = iv_f[1] / C;
            T_f (iv_f[0], iv_f[1], 0) = T_a (i_a, j_a, 0);
            RH_f(iv_f[0], iv_f[1], 0) = RH_a(i_a, j_a, 0);
        });
    }

    for (MFIter mfi(*fire_fuel_mc); mfi.isValid(); ++mfi) {
        Array4<Real>       mc   = fire_fuel_mc->array(mfi);
        Array4<Real>       mext = fire_mext->array(mfi);
        Array4<const Real> T_f  = T_fire.array(mfi);
        Array4<const Real> RH_f = RH_fire.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv_f) {
            int i_f = iv_f[0], j_f = iv_f[1];
            Real T_C = T_f(i_f, j_f, 0) - 273.15_rt;
            Real RH  = RH_f(i_f, j_f, 0) * 100.0_rt;
            mc(i_f,j_f,0,0) = advance_fuel_moisture_one_class(
                mc(i_f,j_f,0,0), RH, T_C, precip_mm_hr, dt_hours, FuelMoistureConst::TAU_1HR);
            mc(i_f,j_f,0,1) = advance_fuel_moisture_one_class(
                mc(i_f,j_f,0,1), RH, T_C, precip_mm_hr, dt_hours, FuelMoistureConst::TAU_10HR);
            mc(i_f,j_f,0,2) = advance_fuel_moisture_one_class(
                mc(i_f,j_f,0,2), RH, T_C, precip_mm_hr, dt_hours, FuelMoistureConst::TAU_100HR);
            Real dead_load       = fp.w_d1 + fp.w_d10 + fp.w_d100;
            Real sigma_weighted  = dead_load > 0.0_rt
                ? (fp.w_d1 * fp.sigma_d1) / dead_load : fp.sigma_d1;
            mext(i_f,j_f,0) = compute_moisture_of_extinction(sigma_weighted);
        });
    }
}

void FireLayer::compute_heat_flux_and_diagnostics(Real dt_fire_s)
{
    FuelModelParams fp = get_anderson_fuel_params(m_params.fuel_model_id);

    // Compute aggregate dead SAV (load-weighted mean of 1-hr, 10-hr, 100-hr classes)
    Real dead_load = fp.w_d1 + fp.w_d10 + fp.w_d100;
    Real sigma_agg = (dead_load > 1.0e-10_rt)
        ? (fp.w_d1*fp.sigma_d1 + fp.w_d10*FIRE_SIGMA_D10 + fp.w_d100*FIRE_SIGMA_D100) / dead_load
        : fp.sigma_d1;

    // Use user-supplied tau_residence_s if > 0, otherwise derive from fuel SAV
    Real tau_res_s = (m_params.tau_residence_s > 0.0_rt)
        ? m_params.tau_residence_s
        : compute_residence_time_s(sigma_agg, fp.rho_p);

    // Fill heat flux and deplete fuel load
    fill_fire_heat_flux(*fire_heat_flux, *fire_fuel_load,
                        *fire_phi, fp, tau_res_s, dt_fire_s);

    // Compute diagnostics: Byram fireline intensity and Thomas flame length
    Real h_kJ_per_kg = fp.heat_content * 2.326_rt;   // BTU/lb -> kJ/kg
    fill_fire_diagnostics(*fire_fireline_intensity, *fire_flame_length,
                          *fire_phi, *fire_ros, *fire_fuel_load,
                          m_fuel_load_initial_kg_m2, h_kJ_per_kg);
}

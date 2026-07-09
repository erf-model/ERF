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

    fire_phi        = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 3);
    fire_wind_ref   = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_wind_eff   = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_wind_extract_z = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_slopes     = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 1);
    fire_curvature  = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_ros        = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_fuel_load  = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_fuel_mc    = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 5, 0);
    fire_mext       = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_heat_flux  = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_spread_vec = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_arrival_time = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_disp_accum   = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_surface_temp = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_surface_rh   = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);

    // Phase 5: Heat flux and diagnostics fields
    fire_fireline_intensity = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_flame_length = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);

    // Phase 9 diagnostics
    fire_flame_temp = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_crown_fraction_burned = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    if (m_params.compute_flame_tilt) {
        fire_flame_tilt = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    }
    if (m_params.crown.enable) {
        fire_crown_active = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
        fire_crown_load = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    }

    fire_phi->setVal(1.0);
    fire_spread_vec->setVal(0.0_rt);
    fire_wind_ref->setVal(0.0);
    fire_wind_eff->setVal(0.0);
    fire_wind_extract_z->setVal(0.0);
    fire_slopes->setVal(0.0);
    fire_curvature->setVal(0.0);
    fire_ros->setVal(0.0);
    fire_arrival_time->setVal(-1.0_rt);
    fire_disp_accum->setVal(0.0_rt);
    fire_surface_temp->setVal(0.0);
    fire_surface_rh->setVal(0.0);
    fire_fireline_intensity->setVal(0.0_rt);
    fire_flame_length->setVal(0.0_rt);
    fire_flame_temp->setVal(0.0_rt);
    fire_crown_fraction_burned->setVal(0.0_rt);
    if (fire_flame_tilt) {
        fire_flame_tilt->setVal(0.0_rt);
    }
    if (fire_crown_active) {
        fire_crown_active->setVal(0.0_rt);
    }
    if (fire_crown_load) {
        fire_crown_load->setVal(m_params.crown.canopy_bulk_den * m_params.crown.canopy_depth);
    }

    // Phase 6: fire-atmosphere coupling MultiFabs
    const BoxArray& ba_atm = erf.boxArray(0);
    const DistributionMapping& dm_atm = erf.DistributionMap(0);
    BoxArray ba_atm_2d = ba_atm;
    ba_atm_2d.coarsen(IntVect(1, 1, erf.Geom(0).Domain().length(2)));
    m_Q_atm_prev = std::make_unique<MultiFab>(ba_atm_2d, dm_atm, 1, 0);
    m_Q_lat_atm_prev = std::make_unique<MultiFab>(ba_atm_2d, dm_atm, 1, 0);
    m_Q_atm_prev->setVal(0.0_rt);
    m_Q_lat_atm_prev->setVal(0.0_rt);

    fire_latent_flux = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_latent_flux->setVal(0.0_rt);

    // Phase 8: Albini spotting diagnostics
    if (m_params.spotting.enable) {
        fire_albini_data = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 4, 0);
        fire_albini_data->setVal(0.0_rt);
        if (m_params.fire_debug) {
            amrex::Print() << "[FIRE DEBUG] Albini spotting enabled: "
                           << "I_B_min=" << m_params.spotting.I_B_min << " kW/m, "
                           << "P_base=" << m_params.spotting.P_base << "\n";
        }
    }

    // Phase 12: Allocate per-cell acceleration state for temporal model.
    // Only allocated when both enable and use_temporal are true.
    // When disabled or size-based, fire_accel_state stays nullptr.
    if (m_params.accel.enable && m_params.accel.use_temporal) {
        fire_accel_state = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 3, 0);
        fire_accel_state->setVal(0.0_rt);
        if (m_params.fire_debug) {
            amrex::Print() << "[FIRE DEBUG] Fire acceleration enabled (temporal model): "
                           << "A_point=" << m_params.accel.A_point << " 1/min, "
                           << "A_line=" << m_params.accel.A_line << " 1/min\n";
        }
    } else if (m_params.accel.enable) {
        if (m_params.fire_debug) {
            amrex::Print() << "[FIRE DEBUG] Fire acceleration enabled (size-based model): "
                           << "L_acc=" << m_params.accel.L_acc << " m\n";
        }
    }

    FuelModelParams fp = get_anderson_fuel_params(fire_params.fuel_model_id);
    fire_fuel_load->setVal((fp.w_d1+fp.w_d10+fp.w_d100+fp.w_lh+fp.w_lw)*4.88243);
    m_fuel_load_initial_kg_m2 = (fp.w_d1+fp.w_d10+fp.w_d100+fp.w_lh+fp.w_lw)*4.88243_rt;
    m_fuel_bed_depth_ft = fp.delta;

    fire_fuel_mc->setVal(0.0);
    for (MFIter mfi(*fire_fuel_mc); mfi.isValid(); ++mfi) {
        Array4<Real> mc = fire_fuel_mc->array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            mc(iv,0) = fire_params.moisture_1hr;
            mc(iv,1) = fire_params.moisture_10hr;
            mc(iv,2) = fire_params.moisture_100hr;
            mc(iv,3) = fire_params.moisture_live;  // live herbaceous (Phase 15, new component)
            mc(iv,4) = fire_params.moisture_live;  // live woody      (Phase 15, new component)
        });
    }

    Real dead_load = fp.w_d1+fp.w_d10+fp.w_d100;
    Real sigma_weighted = dead_load > 0.0_rt ? (fp.w_d1*fp.sigma_d1)/dead_load : fp.sigma_d1;
    fire_mext->setVal(compute_moisture_of_extinction(sigma_weighted));

    compute_terrain_slopes(*fire_slopes, z_phys_nd_atm, erf.Geom(0), m_fg, m_params.terrain_file_name);
    fire_slopes->FillBoundary(m_fg.geom.periodicity());
    compute_terrain_curvature(*fire_curvature, *fire_slopes, m_fg.geom);

    m_ignition_x = fire_params.ignition_x;
    m_ignition_y = fire_params.ignition_y;
    m_ignition_r = fire_params.ignition_r;

    initialize_ignition(*fire_phi, m_fg.geom, m_ignition_x, m_ignition_y, m_ignition_r);
    fire_phi->FillBoundary(m_fg.geom.periodicity());

    for (MFIter mfi(*fire_phi); mfi.isValid(); ++mfi) {
        auto phi_arr = fire_phi->const_array(mfi);
        auto at_arr  = fire_arrival_time->array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            if (phi_arr(iv) < 0.0_rt) at_arr(iv) = 0.0_rt;
        });
    }

    // Phase 10: Load spatial fuel map from file when specified.
    // File reading is CPU-only on rank 0; broadcast to all ranks; copy to device.
    if (!m_params.fuel_map.fuel_map_file.empty()) {
        const int fire_nx = m_fg.ba.minimalBox().length(0);
        const int fire_ny = m_fg.ba.minimalBox().length(1);
        std::vector<int> h_fuel_codes;
        int nodata_val = -9999;
        bool ok = false;
        if (m_params.fuel_map.fuel_map_format == "lcp") {
            ok = read_lcp_fuel_map(m_params.fuel_map.fuel_map_file,
                                   fire_nx, fire_ny, h_fuel_codes);
        } else {
            ok = read_ascii_fuel_map(m_params.fuel_map.fuel_map_file,
                                     fire_nx, fire_ny, h_fuel_codes, nodata_val);
        }
        if (ok) {
            m_d_fuel_codes.resize(h_fuel_codes.size());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             h_fuel_codes.begin(), h_fuel_codes.end(),
                             m_d_fuel_codes.begin());
            fire_fuel_model = std::make_unique<amrex::MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
            fill_fuel_model_mf(*fire_fuel_model, m_d_fuel_codes.data(),
                               m_fg.geom, fire_nx);
            m_has_spatial_fuel = true;
            if (m_params.fire_debug) {
                amrex::Print() << "[FIRE DEBUG] Loaded spatial fuel map '"
                               << m_params.fuel_map.fuel_map_file << "': "
                               << fire_nx << "x" << fire_ny << " cells\n";
            }
        } else {
            amrex::Print() << "[FIRE] WARNING: Cannot read fuel map '"
                           << m_params.fuel_map.fuel_map_file
                           << "'; using uniform fuel_model_id="
                           << m_params.fuel_model_id << "\n";
        }
    }

    // Phase 10: Apply firebreak barriers after ignition stamp.
    if (!m_params.firebreaks.empty() && fire_phi) {
        apply_firebreaks(*fire_phi, m_params.firebreaks, m_fg.geom);
        if (m_params.fire_debug) {
            amrex::Print() << "[FIRE DEBUG] Applied "
                           << m_params.n_firebreaks << " firebreak(s)\n";
        }
    }

    // Phase 11: Polygon ignition (initial fire perimeter from vertex file).
    // Applied at t=0 as part of initialization, before the schedule.
    if (!m_params.ignition.polygon_file.empty()) {
        std::vector<amrex::Real> xs, ys;
        // Vertex file is read on rank 0 only; broadcast to all ranks inside
        // read_polygon_vertices() before returning.
        read_polygon_vertices(m_params.ignition.polygon_file, xs, ys);
        if (m_params.ignition.polygon_type == "polyline") {
            init_phi_from_polyline(*fire_phi, m_fg.geom, xs, ys,
                                   m_params.ignition.polyline_width);
        } else {
            init_phi_from_polygon(*fire_phi, m_fg.geom, xs, ys);
        }
        fire_phi->FillBoundary(m_fg.geom.periodicity());
        if (m_params.fire_debug) {
            amrex::Print() << "[FIRE DEBUG] Polygon ignition applied from '"
                           << m_params.ignition.polygon_file << "' ("
                           << m_params.ignition.polygon_type << ")\n";
        }
    }

    // Phase 11: Load ignition schedule if specified.
    // File reading and broadcast are handled inside load_ignition_schedule().
    if (!m_params.ignition.ignition_schedule_file.empty()) {
        load_ignition_schedule(m_params.ignition.ignition_schedule_file,
                               m_ignition_schedule);
        m_has_schedule = true;
        if (m_params.fire_debug) {
            amrex::Print() << "[FIRE DEBUG] Loaded ignition schedule: "
                           << m_ignition_schedule.events.size()
                           << " events\n";
        }
    }

    m_fp.phi_threshold      = fire_params.farsite_phi_threshold;
    m_fp.coeff_a            = fire_params.farsite_coeff_a;
    m_fp.coeff_b            = fire_params.farsite_coeff_b;
    m_fp.coeff_c            = fire_params.farsite_coeff_c;
    m_fp.use_anderson_lw    = fire_params.farsite_use_anderson_lw;
    m_fp.gaussian_sigma     = fire_params.farsite_gaussian_sigma;
    m_fp.cfl_fire           = fire_params.farsite_cfl_fire;

    m_rc = compute_rothermel_params(fp, fire_params.moisture_1hr,
                                    fire_params.moisture_10hr,
                                    fire_params.moisture_100hr);

    // Phase 13A: Build per-fuel wind height tables and copy to device.
    // When use_per_fuel_wind_ht = false, all entries equal wind_ref_ht (no-op).
    m_use_per_fuel_wind_ht = m_params.use_per_fuel_wind_ht;
    {
        auto h_fcwh = build_fcwh_table(m_params.wind_ref_ht, m_params.use_per_fuel_wind_ht);
        auto h_fcz0 = build_fcz0_table();
        m_d_fcwh.resize(h_fcwh.size());
        m_d_fcz0.resize(h_fcz0.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_fcwh.begin(), h_fcwh.end(), m_d_fcwh.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_fcz0.begin(), h_fcz0.end(), m_d_fcz0.begin());
        if (m_params.fire_debug && m_params.use_per_fuel_wind_ht) {
            amrex::Print() << "[FIRE DEBUG] Per-fuel wind height enabled. "
                           << "FM1 fcwh=" << h_fcwh[1] << " m, FM4 fcwh=" << h_fcwh[4] << " m\n";
        }
    }

    // Phase 13B: Pre-compute alternative ROS model coefficients.
    {
        FuelModelParams fp_ros = get_anderson_fuel_params(m_params.fuel_model_id);
        if (m_params.ros_model == "balbi") {
            m_bc_default = compute_balbi_params(fp_ros, m_params.balbi);
            // Build per-fuel Balbi table when spatial fuel map is active
            if (m_has_spatial_fuel) {
                // TODO: Implement build_fuel_balbi_table if per-fuel variation is needed
                // For now, just use default for all cells
            }
            if (m_params.fire_debug) {
                amrex::Print() << "[FIRE DEBUG] ROS model: Balbi (2009), A_coeff="
                               << m_bc_default.A_coeff << " m/s, v_b="
                               << m_bc_default.v_b << " m/s\n";
            }
        } else if (m_params.ros_model == "cheney_gould") {
            m_cgc = compute_cheney_gould_params(m_params.cheney_gould);
            if (m_params.fire_debug) {
                amrex::Print() << "[FIRE DEBUG] ROS model: Cheney-Gould (1998), "
                               << "moisture=" << m_params.cheney_gould.moisture
                               << "%, curing=" << m_params.cheney_gould.curing << "\n";
            }
        } else if (m_params.ros_model == "behave") {
            // Phase 15: Pre-compute BEHAVE multi-class coefficients.
            FuelModelParams fp_bh = get_anderson_fuel_params(m_params.fuel_model_id);
            m_bs_default = compute_behave_state(fp_bh,
                                                m_params.moisture_1hr,
                                                m_params.moisture_10hr,
                                                m_params.moisture_100hr,
                                                m_params.moisture_live,
                                                m_params.moisture_live);
            if (m_params.fire_debug) {
                amrex::Print() << "[FIRE DEBUG] ROS model: BEHAVE multi-class Rothermel, "
                               << "R0=" << m_bs_default.r_0 * 0.00508_rt << " m/s\n";
            }
        } else if (m_params.ros_model == "macarthur") {
            if (m_params.fire_debug) {
                amrex::Print() << "[FIRE DEBUG] ROS model: MacArthur (1966) Australian formula\n";
            }
        } else {
            // Default: Rothermel — already initialised above via m_rc
            if (m_params.fire_debug) {
                amrex::Print() << "[FIRE DEBUG] ROS model: Rothermel (1972)\n";
            }
        }
    }

    amrex::Print() << "[FIRE] FireLayer initialized: C=" << m_fg.C
                   << ", fuel_model=" << fire_params.fuel_model_id
                   << ", grid=" << m_fg.ba.size() << " boxes" << std::endl;

    if (m_params.fire_debug) {
        IntVect max_extent = m_fg.geom.Domain().size();
        RealBox prob_domain = m_fg.geom.ProbDomain();
        Real dx_fire = (prob_domain.hi(0)-prob_domain.lo(0))/max_extent[0];
        Real dy_fire = (prob_domain.hi(1)-prob_domain.lo(1))/max_extent[1];
        amrex::Print() << "[FIRE DEBUG] Fire grid: dx=" << dx_fire << " m, dy=" << dy_fire
                       << " m, extent=" << max_extent[0] << "x" << max_extent[1]
                       << ", grid_ratio=" << m_fg.C << std::endl;
        amrex::Print() << "[FIRE DEBUG] Coupling type: " << fire_params.coupling_type
                       << " (passive=" << fire_params.is_passive()
                       << ", lagged=" << fire_params.is_lagged()
                       << ", synchronous=" << fire_params.is_synchronous() << ")" << std::endl;
        amrex::Print() << "[FIRE DEBUG] Fire-atmosphere feedback multiplier: "
                       << fire_params.fire_atm_feedback << std::endl;
        amrex::Print() << "[FIRE DEBUG] Heat flux alfg (e-folding height): "
                       << fire_params.heat_flux_alfg << " m" << std::endl;
        amrex::Print() << "[FIRE DEBUG] Inject latent heat: " << fire_params.inject_latent << std::endl;
    }
}

void FireLayer::advance(Real time, Real dt, SurfaceLayer& surface_layer,
                        const MultiFab& xvel, const MultiFab& yvel,
                        const MultiFab& z_phys_cc,
                        const MultiFab& T_atm_k0, const MultiFab& RH_atm_k0)
{
    m_current_time = time;
    m_dt_atm       = dt;
    ++m_step;

    if (m_params.fire_debug)
        amrex::Print() << "[FIRE DEBUG] Starting fire advance step with dt=" << dt << std::endl;

    fill_fire_wind_from_interpolation(*fire_wind_ref, *fire_wind_extract_z, xvel, yvel, z_phys_cc,
                                      m_fg, m_params.wind_ref_ht, m_nz,
                                      m_use_per_fuel_wind_ht ? fire_fuel_model.get() : nullptr,
                                      m_use_per_fuel_wind_ht ? m_d_fcwh.data() : nullptr,
                                      m_use_per_fuel_wind_ht ? 13 : 0);
    if (m_params.fire_debug) {
        amrex::Print() << "[FIRE DEBUG] Wind extraction completed. Max reference wind: "
                       << fire_wind_ref->max(0) << " m/s" << std::endl;
        if (m_step > 0)
            amrex::Print() << "[FIRE DEBUG] Wind extraction height range: min="
                           << fire_wind_extract_z->min(0) << " m  max="
                           << fire_wind_extract_z->max(0) << " m" << std::endl;
    }

    MultiFab::Copy(*fire_wind_eff, *fire_wind_ref, 0, 0, 2, 0);

    if (m_params.use_waf) {
        if (m_params.fire_debug)
            amrex::Print() << "[FIRE DEBUG] Applying Wind Adjustment Factor (formula: "
                           << m_params.waf_formula << ")" << std::endl;
        apply_waf_to_wind();
    }

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

    if (m_params.fire_debug)
        amrex::Print() << "[FIRE DEBUG] Updating fuel moisture from atmospheric state" << std::endl;
    advance_fuel_moisture(dt, T_atm_k0, RH_atm_k0);
    if (m_params.fire_debug) {
        amrex::Print() << "[FIRE DEBUG] Fuel moisture update completed. Max 1-hour moisture: "
                       << fire_fuel_mc->max(0) << std::endl;
        amrex::Print() << "[FIRE DEBUG] Surface temp range: min="
                       << fire_surface_temp->min(0) << " K  max="
                       << fire_surface_temp->max(0) << " K" << std::endl;
        amrex::Print() << "[FIRE DEBUG] Surface RH range:   min="
                       << fire_surface_rh->min(0) << "    max="
                       << fire_surface_rh->max(0) << std::endl;
    }

    if (m_params.moisture_dynamic) {
        long nc = fire_fuel_mc->boxArray().numPts();
        Real avg1   = (nc>0) ? fire_fuel_mc->sum(0)/Real(nc) : m_params.moisture_1hr;
        Real avg10  = (nc>0) ? fire_fuel_mc->sum(1)/Real(nc) : m_params.moisture_10hr;
        Real avg100 = (nc>0) ? fire_fuel_mc->sum(2)/Real(nc) : m_params.moisture_100hr;
        avg1   = amrex::max(0.01_rt, amrex::min(avg1,   0.40_rt));
        avg10  = amrex::max(0.01_rt, amrex::min(avg10,  0.40_rt));
        avg100 = amrex::max(0.01_rt, amrex::min(avg100, 0.40_rt));
        FuelModelParams fp_cur = get_anderson_fuel_params(m_params.fuel_model_id);
        m_rc = compute_rothermel_params(fp_cur, avg1, avg10, avg100);
        if (m_params.fire_debug)
            amrex::Print() << "[FIRE DEBUG] Updated Rothermel coefficients with avg moisture: "
                           << "M_1hr=" << avg1 << " M_10hr=" << avg10
                           << " M_100hr=" << avg100 << " R0=" << m_rc.R0 << " m/s" << std::endl;
        
        // Phase 13B: Moisture coupling for Balbi and Cheney-Gould models
        if (m_params.moisture_dynamic && m_params.ros_model == "balbi") {
            // Recompute Balbi coefficients with updated moisture
            FuelModelParams fp_balbi = get_anderson_fuel_params(m_params.fuel_model_id);
            fp_balbi.Mx = avg1;  // Use 1-hr moisture as representative
            m_bc_default = compute_balbi_params(fp_balbi, m_params.balbi);
        }
        if (m_params.moisture_dynamic && m_params.ros_model == "cheney_gould") {
            // Update Cheney-Gould with current 1-hr moisture converted to percent
            FireParams::CheneyGouldParams cgp_cur = m_params.cheney_gould;
            cgp_cur.moisture = avg1 * 100.0_rt;  // fraction → percent
            m_cgc = compute_cheney_gould_params(cgp_cur);
        }
        // Phase 15: Update BEHAVE state when dynamic moisture is enabled.
        if (m_params.moisture_dynamic && m_params.ros_model == "behave") {
            FuelModelParams fp_bh = get_anderson_fuel_params(m_params.fuel_model_id);
            // Domain-average live moisture from components 3 and 4
            long nc_live = fire_fuel_mc->boxArray().numPts();
            Real avg_lh  = (nc_live > 0) ? fire_fuel_mc->sum(3) / Real(nc_live) : m_params.moisture_live;
            Real avg_lw  = (nc_live > 0) ? fire_fuel_mc->sum(4) / Real(nc_live) : m_params.moisture_live;
            avg_lh = amrex::max(0.30_rt, amrex::min(avg_lh, 2.50_rt));
            avg_lw = amrex::max(0.30_rt, amrex::min(avg_lw, 2.50_rt));
            m_bs_default = compute_behave_state(fp_bh, avg1, avg10, avg100, avg_lh, avg_lw);
        }
    }

    // Phase 11: Apply any scheduled ignition events due this timestep.
    // Time window: (m_current_time - dt, m_current_time].
    if (m_has_schedule && fire_phi) {
        apply_scheduled_ignitions(*fire_phi, m_fg.geom,
                                  m_ignition_schedule,
                                  m_current_time,
                                  m_current_time - dt);
        // fill_boundary after any phi modification to propagate ghost cells
        //amrex::FillBoundary(*fire_phi, m_fg.geom);
        fire_phi->FillBoundary(m_fg.geom.periodicity());        
    }

    // Phase 13B: ROS model dispatch.
    // All models write into fire_ros [m/s].
    // Rothermel is the default (ros_model = "rothermel" or unrecognised string).
    if (m_params.ros_model == "balbi") {
        const BalbiComputed* d_tbl_ptr = m_d_balbi_table.empty() ? nullptr : m_d_balbi_table.data();
        int tbl_sz = static_cast<int>(m_d_balbi_table.size());
        fill_balbi_ros(*fire_ros, *fire_wind_eff, *fire_slopes,
                       m_bc_default,
                       m_has_spatial_fuel ? fire_fuel_model.get() : nullptr,
                       d_tbl_ptr, tbl_sz);
    } else if (m_params.ros_model == "cheney_gould") {
        fill_cheney_gould_ros(*fire_ros, *fire_wind_eff, m_cgc);
    } else if (m_params.ros_model == "behave") {
        // Phase 15: BEHAVE multi-class Rothermel model
        FuelModelParams fp_behave = get_anderson_fuel_params(m_params.fuel_model_id);
        fill_behave_ros(*fire_ros, *fire_wind_eff, *fire_slopes,
                        fp_behave,
                        m_bs_default,
                        m_params.moisture_dynamic ? fire_fuel_mc.get() : nullptr,
                        m_params.moisture_dynamic);
    } else if (m_params.ros_model == "macarthur") {
        fill_macarthur_ros(*fire_ros, *fire_wind_eff);
    } else {
        // Default: Rothermel (1972) — existing code path, unchanged
        compute_ros_field(*fire_ros, *fire_wind_eff, *fire_slopes, m_rc);
    }
    if (m_params.fire_debug) {
        // Compute masked ROS diagnostics (only for burning cells where phi < 0)
        amrex::Real ros_max_burning = 0.0_rt;
        amrex::Real ros_sum_burning = 0.0_rt;
        long n_burning_cells = 0;

        for (amrex::MFIter mfi(*fire_ros); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto ros_arr = fire_ros->const_array(mfi);
            auto phi_arr = fire_phi->const_array(mfi);
            
            amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
                if (phi_arr(i, j, k, 0) < 0.0_rt) {
                    amrex::Real ros_val = ros_arr(i, j, k, 0);
                    if (ros_val > ros_max_burning) {
                        ros_max_burning = ros_val;
                    }
                    ros_sum_burning += ros_val;
                    ++n_burning_cells;
                }
            });
        }

        amrex::ParallelDescriptor::ReduceRealMax(ros_max_burning);
        amrex::ParallelDescriptor::ReduceRealSum(ros_sum_burning);
        amrex::ParallelDescriptor::ReduceLongSum(n_burning_cells);

        amrex::Real ros_mean_burning = (n_burning_cells > 0) ? ros_sum_burning / static_cast<amrex::Real>(n_burning_cells) : 0.0_rt;
        amrex::Print() << "[FIRE DEBUG] Rate-of-spread computed. Max: " << ros_max_burning
                       << " m/s, Mean: " << ros_mean_burning
                       << " m/s" << std::endl;
    }

    // Phase 10: Fuel boundary blending when spatial fuel map is active.
    if (m_has_spatial_fuel &&
        m_params.fuel_map.blending_fraction > 0.0_rt &&
        fire_fuel_model && fire_ros) {
        apply_fuel_boundary_blending(*fire_ros, *fire_fuel_model,
                                      m_params.fuel_map.blending_fraction);
    }

    fire_phi->FillBoundary(m_fg.geom.periodicity());

    // Phase 12: Apply fire acceleration scaling to ROS.
    // Reduces ROS for small fires not yet at quasi-steady-state.
    // Returns immediately when accel.enable = false (zero cost when disabled).
    if (m_params.accel.enable && fire_ros && fire_phi) {
        apply_fire_acceleration(*fire_ros, *fire_phi, m_fg.geom,
                                m_params.accel, dt,
                                fire_accel_state.get(),
                                m_params.fire_debug);
    }

    int n_substeps = 0;
    
    if (m_params.propagation_method == "levelset") {
        // --- Level-set path ---
        // CFL-based subcycling (same structure as FARSITE)
        amrex::Real time_remaining = dt;
        while (time_remaining > 1.0e-14) {
            amrex::Real max_ros = fire_ros->max(0);
            amrex::Real dt_ls   = (max_ros > 1.0e-10)
                ? m_params.levelset_cfl * std::min(m_fg.geom.CellSize()[0],
                                                   m_fg.geom.CellSize()[1]) / max_ros
                : time_remaining;
            dt_ls = std::min(dt_ls, time_remaining);

            fire_levelset::advect_levelset_weno5z_rk3(*fire_phi, *fire_wind_eff,
                                            *fire_ros, m_fg.geom, dt_ls);
            fire_phi->FillBoundary(m_fg.geom.periodicity());

            ++m_levelset_subcycle_count;
            if (m_levelset_subcycle_count % m_params.levelset_reinit_every == 0) {
                amrex::Real dtau = (m_params.levelset_reinit_dtau > 0.0)
                    ? m_params.levelset_reinit_dtau
                    : 0.5 * std::min(m_fg.geom.CellSize()[0], m_fg.geom.CellSize()[1]);
                fire_levelset::reinitialize_phi(*fire_phi, m_fg.geom,
                                      m_params.levelset_reinit_iters, dtau);
                fire_phi->FillBoundary(m_fg.geom.periodicity());
            }

            // Update arrival time for newly burned cells (phi < 0)
            {
                const amrex::Real t_now = m_current_time + (dt - time_remaining);
                for (amrex::MFIter mfi(*fire_phi); mfi.isValid(); ++mfi) {
                    auto p  = fire_phi->const_array(mfi);
                    auto at = fire_arrival_time->array(mfi);
                    amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const amrex::IntVect& iv) noexcept {
                        if (p(iv) < 0.0_rt && at(iv) < 0.0_rt) at(iv) = t_now;
                    });
                }
            }
            time_remaining -= dt_ls;
        }
        n_substeps = m_levelset_subcycle_count;
    } else {
        // --- Default: FARSITE Lagrangian path (unchanged) ---
        n_substeps = advance_fire_subcycle(*fire_phi, *fire_spread_vec,
                                          *fire_disp_accum,
                                          *fire_arrival_time,
                                          *fire_wind_eff, *fire_ros,
                                          m_fg.geom, dt,
                                          m_current_time,
                                          m_fp);
    }

    if (m_params.fire_debug) {
        amrex::Print() << "[FIRE DEBUG] Fire front propagation completed with "
                       << n_substeps << " fire subcycles" << std::endl;
        long num_fire_cells = 0;
        for (MFIter mfi(*fire_phi); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<const Real> phi_arr = fire_phi->array(mfi);
            for (int k=bx.smallEnd(2); k<=bx.bigEnd(2); ++k)
                for (int j=bx.smallEnd(1); j<=bx.bigEnd(1); ++j)
                    for (int i=bx.smallEnd(0); i<=bx.bigEnd(0); ++i)
                        if (phi_arr(i,j,k) < 0.0_rt) ++num_fire_cells;
        }
        amrex::ParallelDescriptor::ReduceLongSum(num_fire_cells);
        amrex::Print() << "[FIRE DEBUG] Number of active fire cells: " << num_fire_cells << std::endl;
    }

    fire_phi->FillBoundary(m_fg.geom.periodicity());

    compute_heat_flux_and_diagnostics(dt);

    // Phase 8: Albini ember spotting
    // Apply stochastic spotting at the specified interval.
    // fire_wind_eff provides the 2-D wind field for trajectory integration.
    // fire_fuel_load provides residual fuel for re-entry filtering.
    if (m_params.spotting.enable && fire_albini_data && fire_wind_eff) {
        if (m_step % m_params.spotting.spotting_interval == 0) {
            fire_albini_data->setVal(0.0_rt);
            FuelModelParams fp_sp = get_anderson_fuel_params(m_params.fuel_model_id);
            std::string fuel_sys  = m_params.spotting.fuel_system;
            compute_albini_spotting(
                *fire_phi,
                *fire_albini_data,
                *fire_wind_eff,
                *fire_ros,
                m_fg.geom,
                fp_sp,
                m_params.spotting,
                m_step,
                fire_fuel_load.get(),
                &fuel_sys,
                m_params.fuel_model_id,
                m_params.fire_debug,
                fire_fuel_load.get(),
                m_fuel_load_initial_kg_m2);
        }
    }

    if (m_params.fire_debug) {
        amrex::Real phi_min  = fire_phi->min(0, 0);   // nghost=0
        amrex::Real phi_max  = fire_phi->max(0, 0);
        
        // Compute masked ROS diagnostics (only for burning cells where phi < 0)
        amrex::Real ros_max_burning = 0.0_rt;
        amrex::Real ros_sum_burning = 0.0_rt;
        long n_burning_cells = 0;

        for (amrex::MFIter mfi(*fire_ros); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto ros_arr = fire_ros->const_array(mfi);
            auto phi_arr = fire_phi->const_array(mfi);
            
            amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
                if (phi_arr(i, j, k, 0) < 0.0_rt) {
                    amrex::Real ros_val = ros_arr(i, j, k, 0);
                    if (ros_val > ros_max_burning) {
                        ros_max_burning = ros_val;
                    }
                    ros_sum_burning += ros_val;
                    ++n_burning_cells;
                }
            });
        }

        amrex::ParallelDescriptor::ReduceRealMax(ros_max_burning);
        amrex::ParallelDescriptor::ReduceRealSum(ros_sum_burning);
        amrex::ParallelDescriptor::ReduceLongSum(n_burning_cells);

        amrex::Real ros_max = ros_max_burning;
        amrex::Real ros_mean = (n_burning_cells > 0) ? ros_sum_burning / static_cast<amrex::Real>(n_burning_cells) : 0.0_rt;
        
        amrex::Real Q_max    = fire_heat_flux ? fire_heat_flux->max(0) : 0.0;
        amrex::Real I_B_max  = fire_fireline_intensity ? fire_fireline_intensity->max(0) : 0.0;
        amrex::Real L_max    = fire_flame_length ? fire_flame_length->max(0) : 0.0;
        amrex::Print() << "[FIRE] t=" << m_current_time
                       << "  substeps=" << n_substeps
                       << "  phi_min=" << phi_min
                       << "  phi_max=" << phi_max
                       << "  max_ROS=" << ros_max << " m/s"
                       << "  mean_ROS=" << ros_mean << " m/s"
                       << "  Q_max=" << Q_max << " W/m2"
                       << "  I_B_max=" << I_B_max << " kW/m"
                       << "  L_max=" << L_max << " m\n";
    }

    if (m_params.write_fire_stats_csv) {
        static bool csv_header_written = false;
        if (!csv_header_written) {
            write_fire_stats_header(m_params.fire_stats_csv_file);
            csv_header_written = true;
        }
        append_fire_stats(*fire_phi, *fire_arrival_time, m_fg.geom,
                         m_step, m_current_time, m_params.fire_stats_csv_file,
                         fire_ros.get(), fire_heat_flux.get(), fire_albini_data.get());
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
    int C = m_fg.C;

    // Map atmospheric T and RH to fire grid and store persistently.
    // This operation runs every timestep to ensure
    // fire_surface_temp and fire_surface_rh are available for plotfile output.
    for (MFIter mfi(*fire_surface_temp, false); mfi.isValid(); ++mfi) {
        Array4<Real> T_f  = fire_surface_temp->array(mfi);
        Array4<Real> RH_f = fire_surface_rh->array(mfi);
        Array4<const Real> T_atm  = T_atm_k0.const_array(mfi);
        Array4<const Real> RH_atm = RH_atm_k0.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv_f) {
            int ia = iv_f[0]/C, ja = iv_f[1]/C;
            T_f(iv_f[0],iv_f[1],0)  = T_atm(ia,ja,0);
            RH_f(iv_f[0],iv_f[1],0) = RH_atm(ia,ja,0);
        });
    }

    if (!m_params.moisture_dynamic) { return; }
    Real dt_hours = dt_s / 3600.0_rt;
    FuelModelParams fp = get_anderson_fuel_params(m_params.fuel_model_id);
    Real precip_mm_hr = m_params.precip_rate_mm_hr;

    for (MFIter mfi(*fire_fuel_mc); mfi.isValid(); ++mfi) {
        Array4<Real> mc   = fire_fuel_mc->array(mfi);
        Array4<Real> mext = fire_mext->array(mfi);
        Array4<const Real> T_f  = fire_surface_temp->const_array(mfi);
        Array4<const Real> RH_f = fire_surface_rh->const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv_f) {
            int i = iv_f[0], j = iv_f[1];
            Real T_C = T_f(i,j,0) - 273.15_rt;
            Real RH  = RH_f(i,j,0) * 100.0_rt;
            // Existing 3 dead fuel classes (unchanged)
            mc(i,j,0,0) = advance_fuel_moisture_one_class(mc(i,j,0,0),RH,T_C,precip_mm_hr,dt_hours,FuelMoistureConst::TAU_1HR);
            mc(i,j,0,1) = advance_fuel_moisture_one_class(mc(i,j,0,1),RH,T_C,precip_mm_hr,dt_hours,FuelMoistureConst::TAU_10HR);
            mc(i,j,0,2) = advance_fuel_moisture_one_class(mc(i,j,0,2),RH,T_C,precip_mm_hr,dt_hours,FuelMoistureConst::TAU_100HR);
            // Phase 15: live fuel moisture (components 3 and 4)
            // Live fuels respond slowly to atmospheric conditions.
            // Use TAU_100HR as a lower bound; live moisture is bounded [0.30, 2.50].
            if (mc.nComp() >= 5) {
                Real lh_new = advance_fuel_moisture_one_class(mc(i,j,0,3),RH,T_C,0.0_rt,dt_hours,FuelMoistureConst::TAU_100HR);
                Real lw_new = advance_fuel_moisture_one_class(mc(i,j,0,4),RH,T_C,0.0_rt,dt_hours,FuelMoistureConst::TAU_100HR);
                mc(i,j,0,3) = amrex::max(0.30_rt, amrex::min(lh_new, 2.50_rt));  // live herba: 30%–250%
                mc(i,j,0,4) = amrex::max(0.30_rt, amrex::min(lw_new, 2.50_rt));  // live woody: 30%–250%
            }
            Real dead_load = fp.w_d1+fp.w_d10+fp.w_d100;
            Real sw = dead_load>0.0_rt ? (fp.w_d1*fp.sigma_d1)/dead_load : fp.sigma_d1;
            mext(i,j,0) = compute_moisture_of_extinction(sw);
        });
    }
}

void FireLayer::compute_heat_flux_and_diagnostics(Real dt_fire_s)
{
    FuelModelParams fp = get_anderson_fuel_params(m_params.fuel_model_id);

    Real dead_load = fp.w_d1 + fp.w_d10 + fp.w_d100;
    Real sigma_agg = (dead_load > 1.0e-10_rt)
        ? (fp.w_d1*fp.sigma_d1 + fp.w_d10*FIRE_SIGMA_D10 + fp.w_d100*FIRE_SIGMA_D100) / dead_load
        : fp.sigma_d1;

    Real tau_sav = compute_residence_time_s(sigma_agg, fp.rho_p);
    Real tau_sav_floor = (m_params.tau_residence_s > 0.0_rt) ? m_params.tau_residence_s : tau_sav;

    if (m_params.fire_debug) {
        // Compute masked ROS diagnostics (only for burning cells where phi < 0)
        amrex::Real ros_max_burning = 0.0_rt;
        amrex::Real ros_sum_burning = 0.0_rt;
        long n_burning_cells = 0;

        for (amrex::MFIter mfi(*fire_ros); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto ros_arr = fire_ros->const_array(mfi);
            auto phi_arr = fire_phi->const_array(mfi);
            
            amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
                if (phi_arr(i, j, k, 0) < 0.0_rt) {
                    amrex::Real ros_val = ros_arr(i, j, k, 0);
                    if (ros_val > ros_max_burning) {
                        ros_max_burning = ros_val;
                    }
                    ros_sum_burning += ros_val;
                    ++n_burning_cells;
                }
            });
        }

        amrex::ParallelDescriptor::ReduceRealMax(ros_max_burning);
        amrex::ParallelDescriptor::ReduceRealSum(ros_sum_burning);
        amrex::ParallelDescriptor::ReduceLongSum(n_burning_cells);

        amrex::Print() << "[FIRE DEBUG] tau_sav=" << tau_sav
                       << " s  (dx_fire=" << m_fg.geom.CellSize(0)
                       << " m, max_ROS=" << ros_max_burning << " m/s)" << std::endl;
    }

    fill_fire_heat_flux(*fire_heat_flux, *fire_fuel_load,
                        *fire_phi, *fire_ros, fp,
                        m_fg.geom.CellSize(0), tau_sav_floor, dt_fire_s,
                        m_has_spatial_fuel ? fire_fuel_model.get() : nullptr);

    const Real h_kJ_per_kg = fp.heat_content * 2.326_rt;
    const Real h_fuel_Jkg = fp.heat_content * 2326.0_rt;

    if (fire_crown_fraction_burned) {
        fire_crown_fraction_burned->setVal(0.0_rt);
    }

    if (m_params.crown.enable && fire_crown_active && fire_crown_load) {
        MultiFab surface_ros(m_fg.ba, m_fg.dm, 1, 0);
        MultiFab surface_intensity(m_fg.ba, m_fg.dm, 1, 0);
        MultiFab surface_flame_length(m_fg.ba, m_fg.dm, 1, 0);
        MultiFab::Copy(surface_ros, *fire_ros, 0, 0, 1, 0);

        fill_fire_diagnostics(surface_intensity, surface_flame_length,
                              *fire_phi, surface_ros, *fire_fuel_load,
                              m_fuel_load_initial_kg_m2, h_kJ_per_kg);

        const auto& crown = m_params.crown;
        const Real canopy_base_ht = crown.canopy_base_ht;
        const Real canopy_bulk_den = crown.canopy_bulk_den;
        const Real canopy_depth = crown.canopy_depth;
        const Real foliar_moisture = crown.foliar_moisture;
        const Real M_c = crown.M_c;
        const Real default_moisture_10hr = m_params.moisture_10hr;
        const Real I_B_crit = van_wagner_critical_intensity(
            canopy_base_ht, foliar_moisture, M_c);
        const Real fixed_u10_ms = (crown.wind_10m_kmh > 0.0_rt)
            ? crown.wind_10m_kmh / 3.6_rt
            : -1.0_rt;
        const Real h_crown_Jkg = crown.h_crown_BTU_lb * 2326.0_rt;
        const int crown_model_id = (crown.ros_model == "rothermel1991") ? 1
            : (crown.ros_model == "van_wagner_proxy") ? 2 : 0;
        const bool use_dynamic_mc = (m_params.moisture_dynamic && fire_fuel_mc != nullptr);
        const bool use_passive_blend = crown.use_passive_blend;

        for (MFIter mfi(*fire_ros); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const phi_arr = fire_phi->const_array(mfi);
            auto const wind_arr = fire_wind_eff->const_array(mfi);
            auto const surface_ros_arr = surface_ros.const_array(mfi);
            auto const surface_I_B_arr = surface_intensity.const_array(mfi);
            Array4<const Real> mc_arr;
            if (use_dynamic_mc) {
                mc_arr = fire_fuel_mc->const_array(mfi);
            }
            auto ros_arr = fire_ros->array(mfi);
            auto crown_active_arr = fire_crown_active->array(mfi);
            auto crown_load_arr = fire_crown_load->array(mfi);
            auto crown_frac_arr = fire_crown_fraction_burned->array(mfi);
            auto heat_flux_arr = fire_heat_flux->array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (phi_arr(i, j, k) >= 0.0_rt) {
                    ros_arr(i, j, k) = surface_ros_arr(i, j, k);
                    crown_active_arr(i, j, k) = 0.0_rt;
                    crown_frac_arr(i, j, k) = 0.0_rt;
                    return;
                }

                const Real R_surface = amrex::max(surface_ros_arr(i, j, k), 0.0_rt);
                const Real I_surface = amrex::max(surface_I_B_arr(i, j, k), 0.0_rt);
                const Real U10_ms = (fixed_u10_ms > 0.0_rt)
                    ? fixed_u10_ms
                    : std::sqrt(wind_arr(i, j, k, 0) * wind_arr(i, j, k, 0)
                              + wind_arr(i, j, k, 1) * wind_arr(i, j, k, 1));
                const Real moisture_10hr = use_dynamic_mc ? mc_arr(i, j, k, 1) : default_moisture_10hr;

                Real R_active = R_surface;
                if (crown_model_id == 1) {
                    R_active = compute_rothermel_1991_crown_ros(R_surface);
                } else if (crown_model_id == 2) {
                    R_active = compute_van_wagner_proxy_ros(canopy_bulk_den, foliar_moisture);
                } else {
                    R_active = cruz_crown_ros(U10_ms, canopy_bulk_den, moisture_10hr);
                }
                R_active = amrex::max(R_active, R_surface);

                const bool crown_now_active = (crown_active_arr(i, j, k) >= 0.5_rt) || (I_surface >= I_B_crit);
                crown_active_arr(i, j, k) = crown_now_active ? 1.0_rt : 0.0_rt;

                Real R_total = R_surface;
                if (use_passive_blend) {
                    R_total = compute_van_wagner_passive_blend(R_surface, R_active, I_surface, I_B_crit);
                } else if (crown_now_active) {
                    R_total = R_active;
                }
                ros_arr(i, j, k) = amrex::max(R_total, 0.0_rt);
                crown_frac_arr(i, j, k) = compute_crown_fraction_burned(ros_arr(i, j, k), R_surface, R_active);

                if (crown_now_active && crown_load_arr(i, j, k) > 0.0_rt) {
                    const Real tau_crown = amrex::max(canopy_depth / amrex::max(R_active, 1.0e-6_rt), 1.0_rt);
                    heat_flux_arr(i, j, k) += crown_load_arr(i, j, k) * h_crown_Jkg / tau_crown;
                    crown_load_arr(i, j, k) *= std::exp(-dt_fire_s / tau_crown);
                    crown_load_arr(i, j, k) = amrex::max(crown_load_arr(i, j, k), 0.0_rt);
                }
            });
        }
    }

    fill_fire_diagnostics(*fire_fireline_intensity, *fire_flame_length,
                          *fire_phi, *fire_ros, *fire_fuel_load,
                          m_fuel_load_initial_kg_m2, h_kJ_per_kg);

    Real M_f = m_params.moisture_1hr;
    if (m_params.moisture_dynamic && fire_fuel_mc) {
        const long nc = fire_fuel_mc->boxArray().numPts();
        const Real avg1 = (nc > 0) ? fire_fuel_mc->sum(0) / Real(nc) : m_params.moisture_1hr;
        const Real avg10 = (nc > 0) ? fire_fuel_mc->sum(1) / Real(nc) : m_params.moisture_10hr;
        const Real avg100 = (nc > 0) ? fire_fuel_mc->sum(2) / Real(nc) : m_params.moisture_100hr;
        M_f = (dead_load > 1.0e-10_rt)
            ? (fp.w_d1 * avg1 + fp.w_d10 * avg10 + fp.w_d100 * avg100) / dead_load
            : avg1;
    }

    if (fire_flame_temp) {
        fill_flame_temperature(*fire_flame_temp, *fire_fireline_intensity, *fire_phi,
                               m_params.flame_temp_method, h_fuel_Jkg, M_f,
                               m_params.flame_temp_T_amb, m_params.fire_debug);
    }

    if (fire_flame_tilt) {
        fill_flame_tilt_angle(*fire_flame_tilt, *fire_fireline_intensity, *fire_wind_eff,
                              m_params.flame_tilt_rho_air, m_params.flame_tilt_T_amb,
                              m_params.fire_debug);
    }
}

void FireLayer::update_atm_flux_buffer(const amrex::Geometry& geom_atm)
{
    if (!m_params.injects_flux()) {
        if (m_params.fire_debug) {
            amrex::Print() << "[FIRE DEBUG] Skipping flux buffer update: coupling_type is passive (injects_flux=false)" << std::endl;
        }
        return;
    }
    if (!fire_heat_flux || !m_Q_atm_prev) { return; }

    if (m_params.fire_debug) {
        amrex::Print() << "[FIRE DEBUG] Updating atmosphere flux buffer for coupling_type="
                       << m_params.coupling_type << ". Current max heat flux: "
                       << fire_heat_flux->max(0) << " W/m2" << std::endl;
    }

    coarsen_fire_flux_to_atm(*m_Q_atm_prev, *fire_heat_flux,
                             geom_atm, m_fg.geom, m_fg.C);

    if (m_params.inject_latent && m_Q_lat_atm_prev) {
        FuelModelParams fp = get_anderson_fuel_params(m_params.fuel_model_id);
        amrex::Real h_fuel_Jkg = fp.heat_content * 2326.0_rt;

        amrex::Real M_f = m_params.moisture_1hr;
        if (m_params.moisture_dynamic && fire_fuel_mc) {
            long nc = fire_fuel_mc->boxArray().numPts();
            amrex::Real avg1   = (nc > 0) ? fire_fuel_mc->sum(0) / amrex::Real(nc) : m_params.moisture_1hr;
            amrex::Real avg10  = (nc > 0) ? fire_fuel_mc->sum(1) / amrex::Real(nc) : m_params.moisture_10hr;
            amrex::Real avg100 = (nc > 0) ? fire_fuel_mc->sum(2) / amrex::Real(nc) : m_params.moisture_100hr;
            amrex::Real dead_load = fp.w_d1 + fp.w_d10 + fp.w_d100;
            M_f = (dead_load > 1e-10_rt)
                ? (fp.w_d1*avg1 + fp.w_d10*avg10 + fp.w_d100*avg100) / dead_load
                : avg1;
        }

        compute_fire_latent_flux(*fire_latent_flux, *fire_heat_flux, M_f, h_fuel_Jkg);
        coarsen_fire_flux_to_atm(*m_Q_lat_atm_prev, *fire_latent_flux,
                                 geom_atm, m_fg.geom, m_fg.C);
        if (m_params.fire_debug) {
            amrex::Print() << "[FIRE DEBUG] Latent heat flux computed. Max latent flux: "
                           << fire_latent_flux->max(0) << " W/m2, fuel moisture: " << M_f << std::endl;
        }
    } else {
        if (m_params.fire_debug) {
            amrex::Print() << "[FIRE DEBUG] Skipping latent heat injection (inject_latent=false or no moisture)" << std::endl;
        }
        m_Q_lat_atm_prev->setVal(0.0_rt);
    }
}

void FireLayer::apply_fire_coupling_to_cc_source(
    amrex::MultiFab& cc_source,
    const amrex::MultiFab& S_old,
    const amrex::MultiFab& z_phys_cc,
    const amrex::Geometry& geom_atm,
    bool has_moisture)
{
    if (!m_params.injects_flux()) { return; }
    if (!m_Q_atm_prev) { return; }
    if (m_params.fire_atm_feedback <= 0.0_rt) { return; }

    if (m_params.fire_debug) {
        amrex::Print() << "[FIRE DEBUG] Applying fire coupling to atmosphere (coupling_type="
                       << m_params.coupling_type << ", feedback=" << m_params.fire_atm_feedback
                       << ", max_Q_prev=" << m_Q_atm_prev->max(0) << " W/m2)" << std::endl;
    }

    const amrex::MultiFab* Q_lat_ptr = (m_params.inject_latent && has_moisture)
        ? m_Q_lat_atm_prev.get() : nullptr;

    // Pass fire_debug so apply_fire_tendency_to_cc_source prints tendency diagnostics.
    // With fire_debug=true you will see per-RK-stage output:
    //   [FIRE COUPLING] Q_atm_max=...  alfg=...
    //   [FIRE COUPLING] RhoTheta tendency sum=...  max=...  expected_surface_max=...
    apply_fire_tendency_to_cc_source(
        cc_source,
        *m_Q_atm_prev,
        Q_lat_ptr,
        z_phys_cc,
        S_old,
        geom_atm,
        m_params.heat_flux_alfg,
        m_params.fire_atm_feedback,
        has_moisture,
        m_params.fire_debug);
}
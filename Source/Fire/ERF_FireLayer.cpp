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

    // Verify prerequisites
    verify_fire_prerequisites(erf, surface_layer_ptr, fire_params);

    // Create fire grid using public AmrCore accessor (boxArray) instead of
    // protected AmrMesh::grids
    m_fg = create_fire_grid(erf.boxArray(0), erf.DistributionMap(0),
                            erf.Geom(0), fire_params.grid_ratio);

    // Store number of vertical levels
    m_nz = erf.Geom(0).Domain().length(2);

    // Allocate MultiFabs on fire grid (fire_phi with 1 ghost cell for gradients)
    fire_phi        = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 1);
    fire_wind_ref   = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_wind_eff   = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_slopes     = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_curvature  = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_ros        = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_fuel_load  = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_fuel_mc    = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 3, 0);
    fire_mext       = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_heat_flux  = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_spread_vec = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);

    // Initialize values
    fire_phi->setVal(1.0);  // All unburned initially
    fire_wind_ref->setVal(0.0);
    fire_wind_eff->setVal(0.0);
    fire_slopes->setVal(0.0);
    fire_curvature->setVal(0.0);
    fire_ros->setVal(0.0);

    // Set fuel load [kg/m²]
    FuelModelParams fp = get_anderson_fuel_params(fire_params.fuel_model_id);
    Real fuel_load_lb_ft2 = fp.w_d1 + fp.w_d10 + fp.w_d100 + fp.w_lh + fp.w_lw;
    Real fuel_load_kg_m2  = fuel_load_lb_ft2 * 4.88243;  // Convert lb/ft² → kg/m²
    fire_fuel_load->setVal(fuel_load_kg_m2);

    // Store fuel bed depth from fuel model [ft]
    m_fuel_bed_depth_ft = fp.delta;

    // Set fuel moisture [fraction]
    fire_fuel_mc->setVal(0.0);
    for (MFIter mfi(*fire_fuel_mc); mfi.isValid(); ++mfi) {
        Array4<Real> mc = fire_fuel_mc->array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            mc(iv, 0) = fire_params.moisture_1hr;
            mc(iv, 1) = fire_params.moisture_10hr;
            mc(iv, 2) = fire_params.moisture_100hr;
        });
    }

    // Initialize moisture of extinction (Phase 4)
    // Use weighted SAV of dead fuels
    Real dead_load = fp.w_d1 + fp.w_d10 + fp.w_d100;
    Real sigma_weighted = dead_load > 0.0_rt
        ? (fp.w_d1 * fp.sigma_d1) / dead_load
        : fp.sigma_d1;
    Real M_x = compute_moisture_of_extinction(sigma_weighted);
    fire_mext->setVal(M_x);

    // Compute terrain slopes and curvature (static fields, computed once at init)
    compute_terrain_slopes(*fire_slopes, z_phys_nd_atm, erf.Geom(0), m_fg, m_params.terrain_file_name);
    compute_terrain_curvature(*fire_curvature, *fire_slopes, m_fg.geom);

    // Store ignition parameters for phase 3
    m_ignition_x = fire_params.ignition_x;
    m_ignition_y = fire_params.ignition_y;
    m_ignition_r = fire_params.ignition_r;

    // Initialize ignition using phase 3 function
    initialize_ignition(*fire_phi, m_fg.geom, m_ignition_x, m_ignition_y, m_ignition_r);
    fire_phi->FillBoundary(m_fg.geom.periodicity());

    // Initialize FARSITE parameters
    m_fp.phi_threshold      = fire_params.farsite_phi_threshold;
    m_fp.coeff_a            = fire_params.farsite_coeff_a;
    m_fp.coeff_b            = fire_params.farsite_coeff_b;
    m_fp.coeff_c            = fire_params.farsite_coeff_c;
    m_fp.use_anderson_lw    = fire_params.farsite_use_anderson_lw;
    m_fp.gaussian_sigma     = fire_params.farsite_gaussian_sigma;
    m_fp.cfl_fire           = fire_params.farsite_cfl_fire;

    // Precompute Rothermel coefficients
    m_rc = compute_rothermel_params(fp, fire_params.moisture_1hr,
                                    fire_params.moisture_10hr,
                                    fire_params.moisture_100hr);

    amrex::Print() << "[FIRE] FireLayer initialized: "
                   << "C=" << m_fg.C << ", "
                   << "fuel_model=" << fire_params.fuel_model_id << ", "
                   << "grid=" << m_fg.ba.size() << " boxes" << std::endl;
    
    // Debug: Print fire grid mesh resolution
    if (m_params.fire_debug) {
        IntVect max_extent = m_fg.geom.Domain().size();
        RealBox prob_domain = m_fg.geom.ProbDomain();
        Real dx_fire = (prob_domain.hi(0) - prob_domain.lo(0)) / max_extent[0];
        Real dy_fire = (prob_domain.hi(1) - prob_domain.lo(1)) / max_extent[1];
        amrex::Print() << "[FIRE DEBUG] Fire grid mesh resolution: dx=" << dx_fire 
                       << " m, dy=" << dy_fire << " m, "
                       << "extent: " << max_extent[0] << " x " << max_extent[1]
                       << " cells, grid_ratio=" << m_fg.C << std::endl;
    }
}

// surface_layer is non-const because SurfaceLayer's get_u_star / get_z0 /
// get_olen accessors are not marked const.
void FireLayer::advance(Real time, Real dt, SurfaceLayer& surface_layer,
                        const MultiFab& xvel,
                        const MultiFab& yvel,
                        const MultiFab& z_phys_cc,
                        const MultiFab& T_atm_k0,
                        const MultiFab& RH_atm_k0)
{
    m_current_time = time;
    m_dt_atm       = dt;

    if (m_params.fire_debug) {
        amrex::Print() << "[FIRE DEBUG] Starting fire advance step with dt=" << dt << std::endl;
    }

    // 1. Extract wind at reference height using direct vertical interpolation
    fill_fire_wind_from_interpolation(*fire_wind_ref,
                                      xvel, yvel, z_phys_cc,
                                      m_fg, m_params.wind_ref_ht, m_nz);

    if (m_params.fire_debug) {
        Real max_wind_ref = fire_wind_ref->max(0);
        amrex::Print() << "[FIRE DEBUG] Wind extraction completed. Max reference wind: " << max_wind_ref << " m/s" << std::endl;
    }

    // 2. Copy reference wind to effective wind
    MultiFab::Copy(*fire_wind_eff, *fire_wind_ref, 0, 0, 2, 0);

    // 3. Apply Wind Adjustment Factor if enabled
    if (m_params.use_waf) {
        if (m_params.fire_debug) {
            amrex::Print() << "[FIRE DEBUG] Applying Wind Adjustment Factor (formula: " << m_params.waf_formula << ")" << std::endl;
        }
        apply_waf_to_wind();
    }

    // 4. Apply terrain wind corrections (FARSITE) if enabled
    if (m_params.use_terrain_wind) {
        if (m_params.fire_debug) {
            amrex::Print() << "[FIRE DEBUG] Applying FARSITE terrain wind corrections" << std::endl;
        }
        apply_farsite_terrain_wind(*fire_wind_eff, *fire_slopes, *fire_curvature,
                                   m_params.k_ridge, m_params.k_shelter,
                                   m_params.k_valley, m_params.k_deflect);
    }

    if (m_params.fire_debug) {
        Real max_wind_eff = fire_wind_eff->max(0);
        amrex::Print() << "[FIRE DEBUG] Effective wind computed. Max effective wind: " << max_wind_eff << " m/s" << std::endl;
    }

    // 5. Update fuel moisture from atmospheric state (Phase 2)
    // Called before ROS computation so updated moisture affects fire spread
    if (m_params.fire_debug) {
        amrex::Print() << "[FIRE DEBUG] Updating fuel moisture from atmospheric state" << std::endl;
    }
    advance_fuel_moisture(dt, T_atm_k0, RH_atm_k0);

    if (m_params.fire_debug) {
        Real max_mc_1hr = fire_fuel_mc->max(0);
        amrex::Print() << "[FIRE DEBUG] Fuel moisture update completed. Max 1-hour moisture: " << max_mc_1hr << std::endl;
    }

    // 6. Compute rate-of-spread
    compute_ros_field(*fire_ros, *fire_wind_eff, *fire_slopes, m_rc);

    if (m_params.fire_debug) {
        Real max_ros_temp  = fire_ros->max(0);
        Real mean_ros_temp = fire_ros->sum(0) / fire_ros->boxArray().numPts();
        amrex::Print() << "[FIRE DEBUG] Rate-of-spread computed. Max: " << max_ros_temp << " m/s, Mean: " << mean_ros_temp << " m/s" << std::endl;
    }

    // 7. (Phase 3) Advance level-set using FARSITE subcycle
    int n_substeps = advance_fire_subcycle(*fire_phi, *fire_spread_vec,
                                           *fire_wind_eff, *fire_ros,
                                           m_fg.geom, dt, m_fp);

    if (m_params.fire_debug) {
        amrex::Print() << "[FIRE DEBUG] Level-set propagation completed with " << n_substeps << " fire subcycles" << std::endl;
    }

    // Fill ghost cells after propagation so gradient stencils in the next
    // ComputeRothermellSpreadRate call (next atmospheric step) are valid.
    fire_phi->FillBoundary(m_fg.geom.periodicity());

    // Diagnostics
    Real max_ros  = fire_ros->max(0);
    Real mean_ros = fire_ros->sum(0) / fire_ros->boxArray().numPts();
    Real phi_min  = fire_phi->min(0);
    
    amrex::Print() << "[FIRE] t=" << m_current_time
                   << "  substeps=" << n_substeps
                   << "  phi_min=" << phi_min
                   << "  max_ROS=" << max_ros << " m/s"
                   << "  mean_ROS=" << mean_ros << " m/s" << std::endl;

    // (Phase 6) Heat flux:             not implemented in Phase 3
}

void FireLayer::init_ignition(Real center_x, Real center_y, Real radius)
{
    // Signed-distance initialisation of the level-set:
    //   fire_phi < 0  →  burned (inside circle)
    //   fire_phi > 0  →  unburned (outside circle)
    //   fire_phi ≈ 0  →  fire front

    const Real dx = m_fg.geom.CellSize(0);
    const Real dy = m_fg.geom.CellSize(1);

    // Copy ProbLo values into local scalars so the GPU lambda captures them
    // by value without implicitly capturing 'this' (fixes -Wdeprecated-this-capture).
    const Real problo_x = m_fg.geom.ProbLo(0);
    const Real problo_y = m_fg.geom.ProbLo(1);

    for (MFIter mfi(*fire_phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> phi = fire_phi->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            const Real x    = problo_x + (iv[0] + 0.5_rt) * dx;
            const Real y    = problo_y + (iv[1] + 0.5_rt) * dy;
            const Real dist = std::sqrt((x - center_x)*(x - center_x) +
                                                (y - center_y)*(y - center_y));
            phi(iv) = dist - radius;
        });
    }
}

void FireLayer::apply_waf_to_wind()
{
    // Wind Adjustment Factor — uses fuel bed depth from fuel model

    Real waf = 0.4_rt;  // Default: ~40% of open wind (typical closed forest)

    if (m_params.waf_formula == "andrews") {
        // Andrews (1982) unsheltered WAF for given fuel-bed depth
        waf = compute_waf_unsheltered(m_fuel_bed_depth_ft);
    } else if (m_params.waf_formula == "behaviorplus") {
        waf = compute_waf_behaviorplus(m_fuel_bed_depth_ft);
    }

    fire_wind_eff->mult(waf, 0, 2, 0);
}
void FireLayer::advance_fuel_moisture(Real dt_s,
                                      const MultiFab& T_atm_k0,
                                      const MultiFab& RH_atm_k0)
{
    if (!m_params.moisture_dynamic) {
        return;  // Skip if moisture evolution is disabled
    }

    // Convert time step from seconds to hours
    Real dt_hours = dt_s / 3600.0_rt;

    // Get fuel model parameters for weighted SAV computation
    FuelModelParams fp = get_anderson_fuel_params(m_params.fuel_model_id);

    // Precipitation rate (uniform, scalar)
    Real precip_mm_hr = m_params.precip_rate_mm_hr;

    // Refinement factor C
    int C = m_fg.C;

    // Create fire-grid versions of atmospheric fields
    // This ensures indices match between fire grid MFIter and atmospheric data access
    MultiFab T_fire(fire_fuel_mc->boxArray(), fire_fuel_mc->DistributionMap(), 1, 0);
    MultiFab RH_fire(fire_fuel_mc->boxArray(), fire_fuel_mc->DistributionMap(), 1, 0);

    // Fill fire-grid MultiFabs by coarse-to-fine prolongation
    // Use host-side loop to map atmospheric grid to fire grid
    for (MFIter mfi(T_fire, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> T_f = T_fire.array(mfi);
        Array4<Real> RH_f = RH_fire.array(mfi);
        Array4<const Real> T_atm = T_atm_k0.const_array(mfi);
        Array4<const Real> RH_atm = RH_atm_k0.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv_f) {
            int i_f = iv_f[0];
            int j_f = iv_f[1];

            // Map fire grid to atmospheric grid
            int i_a = i_f / C;
            int j_a = j_f / C;

            // Fill with atmospheric values
            T_f(i_f, j_f, 0) = T_atm(i_a, j_a, 0);
            RH_f(i_f, j_f, 0) = RH_atm(i_a, j_a, 0);
        });
    }

    // Update each fire cell using fire-grid indices
    for (MFIter mfi(*fire_fuel_mc); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> mc = fire_fuel_mc->array(mfi);
        Array4<Real> mext = fire_mext->array(mfi);

        // Get atmospheric arrays from fire-grid MultiFabs (co-located indices)
        Array4<const Real> T_f = T_fire.array(mfi);
        Array4<const Real> RH_f = RH_fire.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv_f) {
            int i_f = iv_f[0];
            int j_f = iv_f[1];

            // Read atmospheric state (now directly from fire-grid aligned arrays)
            Real T_K = T_f(i_f, j_f, 0);      // Potential temperature [K]
            Real RH_frac = RH_f(i_f, j_f, 0); // Relative humidity [0-1]

            // Convert to Celsius and percent
            Real T_C = T_K - 273.15_rt;
            Real RH = RH_frac * 100.0_rt;

            // Update moisture for each fuel class
            Real M_1hr_new = advance_fuel_moisture_one_class(
                mc(i_f, j_f, 0, 0), RH, T_C, precip_mm_hr, dt_hours,
                FuelMoistureConst::TAU_1HR);
            Real M_10hr_new = advance_fuel_moisture_one_class(
                mc(i_f, j_f, 0, 1), RH, T_C, precip_mm_hr, dt_hours,
                FuelMoistureConst::TAU_10HR);
            Real M_100hr_new = advance_fuel_moisture_one_class(
                mc(i_f, j_f, 0, 2), RH, T_C, precip_mm_hr, dt_hours,
                FuelMoistureConst::TAU_100HR);

            // Write updated values
            mc(i_f, j_f, 0, 0) = M_1hr_new;
            mc(i_f, j_f, 0, 1) = M_10hr_new;
            mc(i_f, j_f, 0, 2) = M_100hr_new;

            // Recompute moisture of extinction from weighted SAV
            Real dead_load = fp.w_d1 + fp.w_d10 + fp.w_d100;
            Real sigma_weighted = dead_load > 0.0_rt
                ? (fp.w_d1 * fp.sigma_d1) / dead_load
                : fp.sigma_d1;
            Real M_x = compute_moisture_of_extinction(sigma_weighted);
            mext(i_f, j_f, 0) = M_x;
        });
    }
}

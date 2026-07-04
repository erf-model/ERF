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

    // Allocate MultiFabs on fire grid (fire_phi with 1 ghost cell for gradients)
    fire_phi        = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 1);
    fire_wind_ref   = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_wind_eff   = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_slopes     = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_curvature  = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_ros        = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_fuel_load  = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_fuel_mc    = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 3, 0);
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

    // Compute terrain slopes and curvature (static fields, computed once at init)
    compute_terrain_slopes(*fire_slopes, z_phys_nd_atm, erf.Geom(0), m_fg);
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
}

// surface_layer is non-const because SurfaceLayer's get_u_star / get_z0 /
// get_olen accessors are not marked const.
void FireLayer::advance(Real dt, SurfaceLayer& surface_layer)
{
    m_current_time = 0.0_rt;  // Time is tracked externally; this is a placeholder
    m_dt_atm       = dt;

    // 1. Extract MOST wind at reference height
    fill_fire_wind_from_most(*fire_wind_ref,
                             *surface_layer.get_u_star(0),
                             *surface_layer.get_z0(0),
                             *surface_layer.get_olen(0),
                             *surface_layer.get_mac_avg(0, 0),
                             *surface_layer.get_mac_avg(0, 1),
                             m_fg, m_params.wind_ref_ht);

    // 2. Copy reference wind to effective wind
    MultiFab::Copy(*fire_wind_eff, *fire_wind_ref, 0, 0, 2, 0);

    // 3. Apply Wind Adjustment Factor if enabled
    if (m_params.use_waf) {
        apply_waf_to_wind();
    }

    // 4. Apply terrain wind corrections (FARSITE) if enabled
    if (m_params.use_terrain_wind) {
        apply_farsite_terrain_wind(*fire_wind_eff, *fire_slopes, *fire_curvature,
                                   m_params.k_ridge, m_params.k_shelter,
                                   m_params.k_valley, m_params.k_deflect);
    }

    // 5. Compute rate-of-spread
    compute_ros_field(*fire_ros, *fire_wind_eff, *fire_slopes, m_rc);

    // 6. (Phase 3) Advance level-set using FARSITE subcycle
    int n_substeps = advance_fire_subcycle(*fire_phi, *fire_spread_vec,
                                           *fire_wind_eff, *fire_ros,
                                           m_fg.geom, dt, m_fp);

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
    // Wind Adjustment Factor — constant per call for now.
    // A full implementation would compute WAF per cell from canopy properties.

    Real waf = 0.4_rt;  // Default: ~40% of open wind (typical closed forest)

    if (m_params.waf_formula == "andrews") {
        // Andrews (1982) unsheltered WAF for ~2 m fuel-bed depth (6.5 ft)
        waf = compute_waf_unsheltered(6.5_rt);
    } else if (m_params.waf_formula == "behaviorplus") {
        waf = compute_waf_behaviorplus(6.5_rt);
    }

    fire_wind_eff->mult(waf, 0, 2, 0);
}
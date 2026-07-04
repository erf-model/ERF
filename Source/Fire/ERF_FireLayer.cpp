#include <ERF_FireLayer.H>
#include <ERF.H>
#include <ERF_SurfaceLayer.H>
#include <ERF_FirePrerequisites.H>
#include <ERF_FireGrid.H>
#include <ERF_FireWindExtract.H>
#include <ERF_TerrainSlope.H>

using namespace amrex;

void FireLayer::initialize(const ERF& erf, const FireParams& fire_params)
{
    m_params = fire_params;

    // Verify prerequisites
    verify_fire_prerequisites(erf, erf.m_SurfaceLayer.get(), fire_params);

    // Create fire grid
    m_fg = create_fire_grid(erf.grids[0], erf.DistributionMap(0),
                            erf.Geom(0), fire_params.grid_ratio);

    // Get atmospheric terrain
    const MultiFab& z_phys_nd_atm = *erf.z_phys_nd[0];

    // Allocate MultiFabs on fire grid
    fire_phi = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_wind_ref = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_wind_eff = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_slopes = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 2, 0);
    fire_curvature = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_ros = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_fuel_load = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
    fire_fuel_mc = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 3, 0);
    fire_heat_flux = std::make_unique<MultiFab>(m_fg.ba, m_fg.dm, 1, 0);
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
    Real fuel_load_kg_m2 = fuel_load_lb_ft2 * 4.88243;  // Convert to SI
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

    // Compute terrain slopes and curvature (static fields)
    compute_terrain_slopes(*fire_slopes, z_phys_nd_atm, erf.Geom(0), m_fg);
    compute_terrain_curvature(*fire_curvature, *fire_slopes, m_fg.geom);

    // Initialize ignition
    init_ignition(fire_params.ignition_x, fire_params.ignition_y,
                  fire_params.ignition_r);

    // Precompute Rothermel coefficients
    m_rc = compute_rothermel_params(fp, fire_params.moisture_1hr,
                                     fire_params.moisture_10hr,
                                     fire_params.moisture_100hr);

    amrex::Print() << "[FIRE] FireLayer initialized: "
                   << "C=" << m_fg.C << ", "
                   << "fuel_model=" << fire_params.fuel_model_id << ", "
                   << "grid=" << m_fg.ba.size() << " boxes" << std::endl;
}

void FireLayer::advance(Real dt, const SurfaceLayer& surface_layer)
{
    // Compute wind extraction pipeline:
    // 1. MOST wind at reference height
    fill_fire_wind_from_most(*fire_wind_ref,
                             *surface_layer.get_u_star(0),
                             *surface_layer.get_z0(0),
                             *surface_layer.get_olen(0),
                             *surface_layer.get_mac_avg(0, 0),
                             *surface_layer.get_mac_avg(0, 1),
                             m_fg, m_params.wind_ref_ht);

    // 2. Copy to effective wind
    MultiFab::Copy(*fire_wind_eff, *fire_wind_ref, 0, 0, 2, 0);

    // 3. Apply WAF if enabled
    if (m_params.use_waf) {
        apply_waf_to_wind();
    }

    // 4. Apply terrain wind corrections if enabled
    if (m_params.use_terrain_wind) {
        apply_farsite_terrain_wind(*fire_wind_eff, *fire_slopes, *fire_curvature,
                                   m_params.k_ridge, m_params.k_shelter,
                                   m_params.k_valley, m_params.k_deflect);
    }

    // 5. Compute rate-of-spread
    compute_ros_field(*fire_ros, *fire_wind_eff, *fire_slopes, m_rc);

    // Print diagnostics
    Real max_ros = fire_ros->max(0);
    Real mean_ros = fire_ros->sum(0) / fire_ros->boxArray().numPts();
    amrex::Print() << "[FIRE] t= " << std::scientific << std::setprecision(6)
                   << "  max_ROS= " << max_ros << " m/s"
                   << "  mean_ROS= " << mean_ros << " m/s" << std::endl;

    // (Phase 3) Level-set advancement: not implemented in Phase 2
    // (Phase 6) Heat flux: not implemented in Phase 2
}

void FireLayer::init_ignition(Real center_x, Real center_y, Real radius)
{
    // Initialize fire level-set as a circle
    // fire_phi < 0 (burned) inside circle
    // fire_phi > 0 (unburned) outside circle

    Real dx = m_fg.geom.CellSize(0);
    Real dy = m_fg.geom.CellSize(1);

    for (MFIter mfi(*fire_phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> phi = fire_phi->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            // Physical coordinates of cell center
            Real x = m_fg.geom.ProbLo(0) + (iv[0] + 0.5) * dx;
            Real y = m_fg.geom.ProbLo(1) + (iv[1] + 0.5) * dy;

            // Distance from ignition center
            Real dist = std::sqrt((x - center_x)*(x - center_x) +
                                          (y - center_y)*(y - center_y));

            // Level-set: distance function
            // Negative inside, positive outside, zero on boundary
            phi(iv) = dist - radius;
        });
    }
}

void FireLayer::apply_waf_to_wind()
{
    // Apply Wind Adjustment Factor to wind field
    // For now, use a constant WAF value
    // In a full implementation, WAF would depend on canopy properties per cell

    Real waf = 0.4;  // Default WAF for forest (about 40% of open wind)

    if (m_params.waf_formula == "andrews") {
        // Andrews formula for fuel bed depth of 2 m
        // WAF ≈ 0.4 for typical forest
        Real h_ft = 6.5;  // ~2 m fuel bed depth
        waf = compute_waf_unsheltered(h_ft);
    } else if (m_params.waf_formula == "behaviorplus") {
        waf = compute_waf_behaviorplus(6.5);
    }

    // Multiply wind by WAF
    fire_wind_eff->mult(waf, 0, 2, 0);
}

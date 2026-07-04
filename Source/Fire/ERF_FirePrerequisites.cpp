#include <ERF_FirePrerequisites.H>
#include <ERF.H>
#include <ERF_SurfaceLayer.H>
#include <ERF_IndexDefines.H>

using namespace amrex;

void verify_fire_prerequisites(const ERF& erf,
                               const SurfaceLayer* surface_layer,
                               const FireParams& fire_params)
{
    // Check 1: Surface layer BC type is "surface_layer"
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        surface_layer != nullptr,
        "fire module requires SurfaceLayer. "
        "Check: phys_bc_type[zlo] == surface_layer AND erf.most.z0 is set");

    // Check 2: SurfaceLayer pointer not null
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        surface_layer != nullptr,
        "[FIRE] Internal error: m_SurfaceLayer is nullptr");

    // Check 3: u_star is computed
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        surface_layer->get_u_star(0) != nullptr,
        "[FIRE] u_star not computed. Call make_SurfaceLayer_at_level(0,...) first");

    // Check 4: z0 is set
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        surface_layer->get_z0(0) != nullptr,
        "[FIRE] Roughness length z0 not set. Add: erf.most.z0 = 0.1");

    // Check 5: Obukhov length is computed
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        surface_layer->get_olen(0) != nullptr,
        "[FIRE] Obukhov length not computed");

    // Check 6: Wind averages are computed
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        (surface_layer->get_mac_avg(0, 0) != nullptr &&
         surface_layer->get_mac_avg(0, 1) != nullptr),
        "[FIRE] Wind averages not computed. Call SurfaceLayer::update_fluxes() first");

    // Check 7: No z-decomposition (all z-levels on same rank)
    const BoxArray& ba = erf.grids[0];
    const Box& domain = erf.Geom(0).Domain();
    int domain_nz = domain.length(2);

    for (int i = 0; i < ba.size(); ++i) {
        const Box& b = ba[i];
        int box_nz = b.length(2);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            box_nz == domain_nz,
            "[FIRE] Cannot decompose in z direction. "
            "Set: erf.max_grid_size_z = " << domain_nz <<
            " or erf.blocking_factor_z = " << domain_nz);
    }

    // Check 8: grid_ratio >= 1
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        fire_params.grid_ratio >= 1,
        "[FIRE] Fire grid_ratio must be >= 1. Set: erf.fire.grid_ratio >= 1");

    // Check 9: All boxes divisible by grid_ratio in x,y
    int C = fire_params.grid_ratio;
    for (int i = 0; i < ba.size(); ++i) {
        const Box& b = ba[i];
        int nx = b.length(0);
        int ny = b.length(1);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            (nx % C == 0 && ny % C == 0),
            "[FIRE] Box sizes not divisible by grid_ratio. "
            "Adjust erf.max_grid_size so all x,y lengths are divisible by " << C);
    }

    // Check 10: DistributionMapping size matches grid size
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        erf.DistributionMap(0).size() == erf.grids[0].size(),
        "[FIRE] Internal error: dmap size != grids size");

    // Check 11: Domain z-index starts at 0 (no AMR z-indexing)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        erf.Geom(0).Domain().smallEnd(2) == 0,
        "[FIRE] Domain z must start at index 0 (AMR not supported)");

    // Check 12: Domain height > wind_ref_ht
    Real prob_hi_z = erf.Geom(0).ProbHi(2);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        prob_hi_z > fire_params.wind_ref_ht,
        "[FIRE] Domain height " << prob_hi_z << " must exceed wind_ref_ht "
        << fire_params.wind_ref_ht << ". Increase geometry.prob_hi(2)");

    amrex::Print() << "[FIRE] All prerequisites verified successfully" << std::endl;
}

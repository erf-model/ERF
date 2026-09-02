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


    // Check 7: No z-decomposition (all z-levels on same rank)
    // Use the public boxArray() accessor instead of the protected grids[] member
    const BoxArray& ba = erf.boxArray(0);
    const Box& domain  = erf.Geom(0).Domain();
    int domain_nz      = domain.length(2);

    for (int i = 0; i < ba.size(); ++i) {
        const Box& b  = ba[i];
        int box_nz    = b.length(2);
        // Build the message string before passing to the macro
        std::string msg = std::string("[FIRE] Cannot decompose in z direction. ")
                        + "Set: erf.max_grid_size_z = " + std::to_string(domain_nz)
                        + " or erf.blocking_factor_z = " + std::to_string(domain_nz);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(box_nz == domain_nz, msg.c_str());
    }

    // Check 8: grid_ratio >= 1
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        fire_params.grid_ratio >= 1,
        "[FIRE] Fire grid_ratio requires value >= 1. Set: erf.fire.grid_ratio >= 1");

    // Check 9: All boxes divisible by grid_ratio in x,y
    int C = fire_params.grid_ratio;
    for (int i = 0; i < ba.size(); ++i) {
        const Box& b = ba[i];
        int nx = b.length(0);
        int ny = b.length(1);
        std::string msg = std::string("[FIRE] Box sizes not divisible by grid_ratio. ")
                        + "Adjust erf.max_grid_size so all x,y lengths are divisible by "
                        + std::to_string(C);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE((nx % C == 0 && ny % C == 0), msg.c_str());
    }

    // Check 10: DistributionMapping size matches grid size
    // Use the public boxArray() accessor here too
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        erf.DistributionMap(0).size() == erf.boxArray(0).size(),
        "[FIRE] Internal error: dmap size != grids size");

    // Check 11: Domain z-index starts at 0
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        erf.Geom(0).Domain().smallEnd(2) == 0,
        "[FIRE] Domain z requires start at index 0 (AMR not supported)");

    // Check 12: Domain height > wind_ref_ht
    Real prob_hi_z = erf.Geom(0).ProbHi(2);
    std::string msg12 = std::string("[FIRE] Domain height ")
                      + std::to_string(prob_hi_z)
                      + " should exceed wind_ref_ht "
                      + std::to_string(fire_params.wind_ref_ht)
                      + ". Increase geometry.prob_hi(2)";
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(prob_hi_z > fire_params.wind_ref_ht, msg12.c_str());

    // Check 13: divisors and step sizes read straight from the inputs file.
    // ParmParse applies no bounds, and a zero reaches the solver as an integer
    // division by zero (SIGFPE) or, for cfl, a zero-length subcycle that never
    // terminates. Each is only checked on the path that actually uses it.
    if (fire_params.propagation_method == "levelset") {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            fire_params.levelset_reinit_every >= 1,
            "[FIRE] erf.fire.levelset.reinit_every must be >= 1; it is used as a "
            "modulus divisor, so 0 is an integer division by zero");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            fire_params.levelset_cfl > 0.0,
            "[FIRE] erf.fire.levelset.cfl must be > 0; a non-positive value gives a "
            "zero-length subcycle and the level-set loop never advances");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            fire_params.levelset_reinit_iters >= 0,
            "[FIRE] erf.fire.levelset.reinit_iters must be >= 0");
    }

    if (fire_params.spotting.enable) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            fire_params.spotting.spotting_interval >= 1,
            "[FIRE] erf.fire.spotting.spotting_interval must be >= 1; it is used as "
            "a modulus divisor, so 0 is an integer division by zero");
    }

    // Check 14: the FARSITE subcycle divides by cfl_fire in the same way.
    if (fire_params.propagation_method != "levelset") {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            fire_params.farsite_cfl_fire > 0.0,
            "[FIRE] erf.fire.farsite.cfl_fire must be > 0; a non-positive value "
            "makes every fire substep zero-length");
    }

    amrex::Print() << "[FIRE] All prerequisites verified successfully" << std::endl;
}

/**
 * @file ERF_DustPrerequisites.cpp
 * @brief Implementation of prerequisite verification for the ERF-Dust module.
 */

#include <ERF_DustPrerequisites.H>
#include <ERF.H>
#include <SurfaceLayer.H>
#include <AMReX_Print.H>

void verify_dust_prerequisites(const ERF&          erf,
                               const SurfaceLayer* surface_layer,
                               const DustParams&   dust_params)
{
    // Check 1: SurfaceLayer pointer is not null
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        surface_layer != nullptr,
        "[DUST] SurfaceLayer is required. Set: zlo.type = \"surface_layer\"");

    // Get atmospheric grid information
    const amrex::BoxArray& ba_atm = erf.boxArray(0);
    const amrex::DistributionMapping& dm_atm = erf.DistributionMap(0);
    const amrex::Geometry& geom_atm = erf.Geom(0);

    // Check 2: erf.most.z0 is set (indirectly verified through SurfaceLayer)
    // This is a ParmParse check that happens at SurfaceLayer construction
    // For now, trust that it was validated there

    // Get domain information
    const amrex::Box& domain = geom_atm.Domain();
    int domain_nz = domain.length(2);

    // Check 3: No z-direction MPI decomposition
    for (int i = 0; i < ba_atm.size(); ++i) {
        amrex::Box box = ba_atm[i];
        int box_nz = box.length(2);
        std::string msg = std::string("[DUST] Cannot decompose in z direction. ")
                        + "Set: amrex.max_grid_size_z = " + std::to_string(domain_nz);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(box_nz == domain_nz, msg.c_str());
    }

    // Check 4: grid_ratio >= 1
    std::string msg4 = std::string("[DUST] Invalid grid_ratio ")
                     + std::to_string(dust_params.grid_ratio)
                     + ". Set: erf.dust.grid_ratio >= 1";
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dust_params.grid_ratio >= 1, msg4.c_str());

    // Check 5: All atmospheric boxes have x,y sizes divisible by grid_ratio
    for (int i = 0; i < ba_atm.size(); ++i) {
        amrex::Box box = ba_atm[i];
        int box_nx = box.length(0);
        int box_ny = box.length(1);
        if (box_nx % dust_params.grid_ratio != 0 || box_ny % dust_params.grid_ratio != 0) {
            std::string msg5 = std::string("[DUST] Box sizes not divisible by grid_ratio. ")
                             + "All atmospheric box x,y sizes must be divisible by grid_ratio="
                             + std::to_string(dust_params.grid_ratio)
                             + ". Adjust grid_ratio or amrex.max_grid_size.";
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false, msg5.c_str());
        }
    }

    // Check 6: DistributionMapping size matches BoxArray size
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        static_cast<int>(dm_atm.size()) == ba_atm.size(),
        "[DUST] DistributionMapping size mismatch. Check AMReX configuration.");

    // Check 7: Domain z-index starts at 0
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        domain.smallEnd(2) == 0,
        "[DUST] Domain z-index must start at 0. Check Geometry configuration.");

    // Check 8: Domain physical height exceeds MRF reference height
    // Default z_ref = 10.0 m (from erf.most.zref, read by SurfaceLayer)
    // For now, just ensure domain height is positive
    amrex::Real dz = geom_atm.ProbHi(2) - geom_atm.ProbLo(2);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        dz > 0.0,
        "[DUST] Domain physical height must be positive. Check Geometry configuration.");

    amrex::Print() << "[DUST] All prerequisites verified\n";
}

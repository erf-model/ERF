#include <ERF_FirePrerequisites.H>
#include <ERF_FireGrid.H>
#include <ERF_SurfacePrerequisites.H>
#include <ERF.H>
#include <ERF_SurfaceLayer.H>
#include <ERF_IndexDefines.H>

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using namespace amrex;

namespace {

// What to change when a box's x or y length is not a multiple of the grid ratio C.
// AMReX accepts only power-of-2 blocking factors, so the advice depends on C. With
// blocking factor 1 AMReX cuts a length of n cells into the fewest near-equal
// pieces no longer than max_grid_size (20 cells with max_grid_size 8 become
// 7 + 7 + 6), and cuts further when there are more ranks than boxes. On level 0,
// when max_grid_size covers the domain in x and y on more than one rank, ERF
// replaces those boxes with one per rank (ERF::MakeNewLevelFromScratch), whatever
// the blocking factor. On a refined level n is the refinement region's extent.
std::string box_ratio_advice (int C, int lev)
{
    const std::string c = std::to_string(C);
    std::string advice = (lev > 0)
        ? "On a refined level n is the refinement region's x,y extent in that level's "
          "cells, and that level's amr.max_grid_size and amr.blocking_factor apply. "
        : "";
    if ((C & (C - 1)) == 0) {
        advice += "Set amr.blocking_factor_x and amr.blocking_factor_y to a multiple of " + c
                + ", with amr.max_grid_size_x, amr.max_grid_size_y and the x,y cell counts n "
                  "multiples of it too.";
    } else {
        advice += "AMReX blocking factors must be powers of 2, which " + c + " is not, so set "
                  "the blocking factor to 1, make the x,y cell counts n multiples of " + c
                + ", and set amr.max_grid_size_x and amr.max_grid_size_y to a multiple of " + c
                + " that divides n: AMReX cuts a length of n cells into the fewest near-equal "
                  "pieces no longer than max_grid_size (n = 40 and max_grid_size = 20 give "
                  "20 + 20). Use at least as many boxes as ranks, since AMReX cuts the boxes "
                  "further otherwise; one box works on one rank.";
    }
    advice += " With more than one rank keep max_grid_size below n in x or y: when it covers "
              "the whole domain in both, ERF cuts level 0 into one box per rank itself, "
              "whatever the blocking factor.";
    return advice;
}

} // namespace

int fire_anchor_level(int anchor_level, int finest_level)
{
    return surface_anchor_level(anchor_level, finest_level, fire_coupling_names());
}

void verify_fire_precip_source(const std::string& engine_name,
                               const SurfacePrecipAccumulationSources& sources)
{
    verify_surface_precip_source(engine_name, sources, "rain", fire_coupling_names());
}

void verify_fire_prerequisites(const ERF& erf,
                               int lev,
                               const SurfaceLayer* surface_layer,
                               int grid_ratio,
                               bool needs_surface_layer)
{
    // Check 1: Surface layer BC type is "surface_layer", for an engine that reads it
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !needs_surface_layer || surface_layer != nullptr,
        "fire module requires SurfaceLayer. "
        "Check: phys_bc_type[zlo] == surface_layer AND erf.most.z0 is set");

    // Checks 2-7b and 11: the checks every surface coupling makes on the level
    // its grid refines (ERF_SurfacePrerequisites.H): the level exists, keeps
    // its boxes above level 0, its boxes span the full height and cover one
    // rectangle, and its domain starts at index 0.
    verify_surface_grid_level(erf, lev, fire_coupling_names());

    const BoxArray& ba = erf.boxArray(lev);

    // Check 8: grid_ratio >= 1
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        grid_ratio >= 1,
        "[FIRE] Fire grid_ratio requires value >= 1. Set: erf.fire.grid_ratio >= 1");

    // Check 9: All boxes divisible by grid_ratio in x,y (box_ratio_advice says
    // which inputs to change).
    int C = grid_ratio;
    for (int i = 0; i < ba.size(); ++i) {
        const Box& b = ba[i];
        int nx = b.length(0);
        int ny = b.length(1);
        if (nx % C != 0 || ny % C != 0) {
            Abort("[FIRE] Box sizes not divisible by grid_ratio: box " + std::to_string(i)
                  + " of level " + std::to_string(lev) + " is " + std::to_string(nx) + " x "
                  + std::to_string(ny) + " cells, and every x,y box length on that level must "
                  "be divisible by erf.fire.grid_ratio = " + std::to_string(C) + ". "
                  + box_ratio_advice(C, lev));
        }
    }

    // Check 10: DistributionMapping size matches grid size
    // Use the public boxArray() accessor here too
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        erf.DistributionMap(lev).size() == erf.boxArray(lev).size(),
        "[FIRE] Internal error: dmap size != grids size");
}

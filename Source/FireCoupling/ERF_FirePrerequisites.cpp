#include <ERF_FirePrerequisites.H>
#include <ERF.H>
#include <ERF_SurfaceLayer.H>
#include <ERF_IndexDefines.H>
#include <AMReX_ParmParse.H>

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
    const int lev = (anchor_level < 0) ? finest_level : anchor_level;
    if (lev > finest_level) {
        Abort("[FIRE] erf.fire.anchor_level = " + std::to_string(anchor_level)
              + " is above the finest level of this run (" + std::to_string(finest_level)
              + "). Set it to a level from 0 to " + std::to_string(finest_level)
              + ", or leave it unset for the finest level.");
    }
    return lev;
}

void verify_fire_precip_source(const std::string& engine_name,
                               const SurfacePrecipAccumulationSources& sources)
{
    if (surface_precip_has_any_source(sources)) { return; }
    std::string moisture_model = "None";
    ParmParse pp_erf("erf");
    pp_erf.query("moisture_model", moisture_model);
    Abort("[FIRE] The fire engine \"" + engine_name + "\" takes the rain of every column from the "
          "microphysics surface precipitation accumulators, but erf.moisture_model = \"" + moisture_model
          + "\" provides none (Kessler, SAM, Morrison, WSM6 and WDM6 with precipitation do; Kessler_NoRain, "
          "SatAdj, MoistNoCondensation and None do not). Choose one of those, or set the engine's inputs so "
          "it does not read the rain.");
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

    // Check 2: the level the fire grid refines exists. fire_anchor_level()
    // resolves the input and aborts above the finest level; this is the backstop
    // for any other caller.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev >= 0 && lev <= erf.finestLevel(),
        "[FIRE] Internal error: the fire grid's level is outside 0..finest_level");

    // Check 3: static refinement. The fire grid is built once from the boxes of
    // its level, and its heat, moisture and smoke go into that level's source,
    // so the level must keep its boxes for the whole run. Level 0 never regrids.
    if (lev > 0) {
        ParmParse pp_erf("erf");
        int regrid_int = -1;
        pp_erf.query("regrid_int", regrid_int);
        if (regrid_int > 0) {
            Abort("[FIRE] The fire grid on level " + std::to_string(lev)
                  + " is built once from that level's boxes, so the refinement must not change, but "
                  "erf.regrid_int = " + std::to_string(regrid_int) + " regrids it. Use static refinement "
                  "boxes with erf.regrid_int = -1, or set erf.fire.anchor_level = 0.");
        }
    }

    // Check 7: No z-decomposition (all z-levels on same rank)
    // Use the public boxArray() accessor instead of the protected grids[] member
    const BoxArray& ba = erf.boxArray(lev);
    const Box& domain  = erf.Geom(lev).Domain();
    int domain_nz      = domain.length(2);

    for (int i = 0; i < ba.size(); ++i) {
        const Box& b  = ba[i];
        int box_nz    = b.length(2);
        // Build the message string before passing to the macro
        std::string msg = std::string("[FIRE] Cannot decompose in z direction. ")
                        + "Every box on level " + std::to_string(lev)
                        + " must span the level's full height, since the fire model interpolates "
                          "through whole columns: set amr.max_grid_size_z = " + std::to_string(domain_nz)
                        + " (or larger), and on a refined level use a refinement box that reaches the "
                          "domain top (erf.<box>.in_box_hi without a z bound)";
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(box_nz == domain_nz, msg.c_str());
    }

    // Check 7b: the boxes cover one rectangle. The fire grid treats the edges of
    // their bounding rectangle as its domain edges and fills its ghost cells
    // there; a hole or a second patch would leave cells with no neighbour data.
    {
        BoxList bl_2d;
        for (int i = 0; i < ba.size(); ++i) {
            Box b = ba[i];
            b.setSmall(2, 0);
            b.setBig(2, 0);
            bl_2d.push_back(b);
        }
        const BoxArray ba_2d(std::move(bl_2d));
        const Box region = ba_2d.minimalBox();
        if (region.numPts() != ba_2d.numPts()) {
            Abort("[FIRE] The boxes of level " + std::to_string(lev) + " cover "
                  + std::to_string(ba_2d.numPts()) + " of the " + std::to_string(region.numPts())
                  + " columns of their bounding rectangle, but the fire grid needs one rectangle. "
                  "Refine one rectangle on that level (a single refinement box), or set "
                  "erf.fire.anchor_level to a level whose boxes form one.");
        }
    }

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

    // Check 11: Domain z-index starts at 0
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        erf.Geom(lev).Domain().smallEnd(2) == 0,
        "[FIRE] Domain z requires start at index 0");

    // Check 11b: the level's domain starts at index 0 in x and y.
    //
    // A fire cell maps back to its atmospheric cell by i_a = i_f / C + atm_lo,
    // with atm_lo the lower corner of the fire region in the level's index space
    // (FireGrid::atm_lo), in ERF_FireGrid.H, ERF_FireAtmCoupling.H,
    // ERF_FireWindExtract.cpp, ERF_FireSurfaceHeight.cpp and the engine's own maps; create_fire_grid() takes
    // the region's physical position from ProbLo plus its index offset from the
    // domain's corner. Every AMReX level has a zero-based domain, so this
    // asserts the convention rather than relying on it silently.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        erf.Geom(lev).Domain().smallEnd(0) == 0 && erf.Geom(lev).Domain().smallEnd(1) == 0,
        "[FIRE] Domain x and y require start at index 0: the fire-to-atmosphere "
        "index mapping assumes a zero-based domain");
}

#include <ERF_FloodPrerequisites.H>
#include <ERF.H>
#include <AMReX_ParmParse.H>

#include <string>

using namespace amrex;

int
flood_anchor_level (int anchor_level, int finest_level)
{
    const int lev = (anchor_level < 0) ? finest_level : anchor_level;
    if (lev > finest_level) {
        Abort("[FLOOD] erf.flood.anchor_level = " + std::to_string(anchor_level)
              + " is above the finest level of this run (" + std::to_string(finest_level)
              + "). Set it to a level from 0 to " + std::to_string(finest_level)
              + ", or leave it unset for the finest level.");
    }
    return lev;
}

void
verify_flood_precip_source (const std::string& engine_name,
                            const SurfacePrecipAccumulationSources& sources)
{
    if (surface_precip_has_any_source(sources)) { return; }
    std::string moisture_model = "None";
    ParmParse pp_erf("erf");
    pp_erf.query("moisture_model", moisture_model);
    Abort("[FLOOD] The flood engine \"" + engine_name + "\" takes the precipitation of every column from "
          "the microphysics surface precipitation accumulators, but erf.moisture_model = \"" + moisture_model
          + "\" provides none (Kessler, SAM, Morrison, WSM6 and WDM6 with precipitation do; Kessler_NoRain, "
          "SatAdj, MoistNoCondensation and None do not). Choose one of those, or set the engine's inputs so "
          "it does not read the precipitation.");
}

void
verify_flood_prerequisites (const ERF& erf, int lev)
{
    // Check 1: the level the flood grid refines exists. flood_anchor_level()
    // resolves the input and aborts above the finest level; this is the backstop
    // for any other caller.
    if (lev < 0 || lev > erf.finestLevel()) {
        Abort("[FLOOD] Internal error: the flood grid's level " + std::to_string(lev)
              + " is outside 0.." + std::to_string(erf.finestLevel()));
    }

    // Check 2: static refinement. The flood grid is built once from the boxes
    // of its level, so the level must keep its boxes for the whole run. Level 0
    // never regrids.
    if (lev > 0) {
        ParmParse pp_erf("erf");
        int regrid_int = -1;
        pp_erf.query("regrid_int", regrid_int);
        if (regrid_int > 0) {
            Abort("[FLOOD] The flood grid on level " + std::to_string(lev)
                  + " is built once from that level's boxes, so the refinement must not change, but "
                  "erf.regrid_int = " + std::to_string(regrid_int) + " regrids it. Use static refinement "
                  "boxes with erf.regrid_int = -1, or set erf.flood.anchor_level = 0.");
        }
    }

    // Check 3: no z decomposition. The precipitation and the plot fields live
    // on the k = 0 slab of each box of the level, one per box, paired box for
    // box with the flood grid.
    const BoxArray& ba = erf.boxArray(lev);
    const Box& domain  = erf.Geom(lev).Domain();
    const int domain_nz = domain.length(2);
    for (int i = 0; i < ba.size(); ++i) {
        if (ba[i].length(2) != domain_nz) {
            Abort("[FLOOD] Every box on level " + std::to_string(lev)
                  + " must span the level's full height, since the flood grid pairs one surface slab "
                    "with each box: leave amr.no_box_split_dir at its default of 2, or set "
                    "amr.max_grid_size_z to at least " + std::to_string(domain_nz)
                  + ", and on a refined level use a refinement box that reaches the domain top");
        }
    }

    // Check 4: the boxes cover one rectangle. The flood grid treats the edges
    // of their bounding rectangle as its domain edges.
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
            Abort("[FLOOD] The boxes of level " + std::to_string(lev) + " cover "
                  + std::to_string(ba_2d.numPts()) + " of the " + std::to_string(region.numPts())
                  + " columns of their bounding rectangle, but the flood grid needs one rectangle. "
                  "Refine one rectangle on that level (a single refinement box), or set "
                  "erf.flood.anchor_level to a level whose boxes form one.");
        }
    }

    // Check 5: a zero-based domain. A flood cell maps back to its atmospheric
    // cell by i_a = i_f / C + atm_lo (ERF_FloodGrid.H); every AMReX level has a
    // zero-based domain, so this asserts the convention rather than relying on
    // it silently.
    if (domain.smallEnd(0) != 0 || domain.smallEnd(1) != 0 || domain.smallEnd(2) != 0) {
        Abort("[FLOOD] Internal error: the flood coupling needs a domain that starts at index 0");
    }
}

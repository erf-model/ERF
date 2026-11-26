#include <ERF_SurfacePrerequisites.H>
#include <ERF.H>
#include <AMReX_ParmParse.H>

#include <string>

using namespace amrex;

int
surface_anchor_level (int anchor_level, int finest_level, const SurfaceCouplingNames& names)
{
    const int lev = (anchor_level < 0) ? finest_level : anchor_level;
    if (lev > finest_level) {
        Abort(names.tag() + " " + names.key("anchor_level") + " = " + std::to_string(anchor_level)
              + " is above the finest level of this run (" + std::to_string(finest_level)
              + "). Set it to a level from 0 to " + std::to_string(finest_level)
              + ", or leave it unset for the finest level.");
    }
    return lev;
}

void
verify_surface_precip_source (const std::string& engine_name,
                              const SurfacePrecipAccumulationSources& sources,
                              const std::string& what,
                              const SurfaceCouplingNames& names)
{
    if (surface_precip_has_any_source(sources)) { return; }
    std::string moisture_model = "None";
    ParmParse pp_erf("erf");
    pp_erf.query("moisture_model", moisture_model);
    Abort(names.tag() + " The " + names.lower + " engine \"" + engine_name + "\" takes the " + what
          + " of every column from the microphysics surface precipitation accumulators, but "
            "erf.moisture_model = \"" + moisture_model
          + "\" provides none (Kessler, SAM, Morrison, WSM6 and WDM6 with precipitation do; Kessler_NoRain, "
            "SatAdj, MoistNoCondensation and None do not). Choose one of those, or set the engine's inputs so "
            "it does not read the " + what + ".");
}

void
verify_surface_grid_level (const ERF& erf, int lev, const SurfaceCouplingNames& names)
{
    const std::string tag  = names.tag();
    const std::string grid = names.grid();

    // Check 1: the level the grid refines exists. surface_anchor_level()
    // resolves the input and aborts above the finest level; this is the
    // backstop for any other caller.
    if (lev < 0 || lev > erf.finestLevel()) {
        Abort(tag + " Internal error: the " + grid + "'s level " + std::to_string(lev)
              + " is outside 0.." + std::to_string(erf.finestLevel()));
    }

    // Check 2: static refinement, with regrid_int read the way ERF reads it
    verify_surface_static_refinement(lev, erf.pp_prefix, names);

    // Check 3: no z decomposition. The grid pairs one surface slab with each
    // box of the level.
    const BoxArray& ba = erf.boxArray(lev);
    const Box& domain  = erf.Geom(lev).Domain();
    const int domain_nz = domain.length(2);
    for (int i = 0; i < ba.size(); ++i) {
        if (ba[i].length(2) != domain_nz) {
            Abort(tag + " Cannot decompose in z direction. Every box on level " + std::to_string(lev)
                  + " must span the level's full height, since the " + grid + " pairs one surface slab "
                    "with each box: leave amr.no_box_split_dir at its default "
                    "of 2, or set amr.max_grid_size_z to at least " + std::to_string(domain_nz)
                  + ", and on a refined level use a refinement box that reaches the domain top "
                    "(erf.<box>.in_box_hi without a z bound)");
        }
    }

    // Check 4: the boxes cover one rectangle. The grid treats the edges of
    // their bounding rectangle as its domain edges and fills its ghost cells
    // there; a hole or a second patch would leave cells with no neighbor data.
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
            Abort(tag + " The boxes of level " + std::to_string(lev) + " cover "
                  + std::to_string(ba_2d.numPts()) + " of the " + std::to_string(region.numPts())
                  + " columns of their bounding rectangle, but the " + grid + " needs one rectangle. "
                    "Refine one rectangle on that level (a single refinement box), or set "
                  + names.key("anchor_level") + " to a level whose boxes form one.");
        }
    }

    // Check 5: a zero-based domain. A surface cell maps back to its
    // atmospheric cell by i_a = i_s / C + atm_lo (ERF_SurfaceGrid.H), and
    // create_surface_grid() takes the region's physical position from ProbLo
    // plus its index offset from the domain's corner. Every AMReX level has a
    // zero-based domain, so this asserts the convention rather than relying on
    // it silently.
    if (domain.smallEnd(0) != 0 || domain.smallEnd(1) != 0 || domain.smallEnd(2) != 0) {
        Abort(tag + " Internal error: the " + names.lower
              + " coupling needs a level domain that starts at index 0");
    }
}

void
verify_surface_static_refinement (int lev, const std::string& pp_prefix, const SurfaceCouplingNames& names)
{
    // The grid is built once from the boxes of its level, so the level must
    // keep its boxes for the whole run. Level 0 never regrids. ERF reads
    // regrid_int with its own prefix (ERF::ReadParameters), which a MultiBlock
    // block sets to its own name.
    if (lev == 0) { return; }
    ParmParse pp(pp_prefix);
    int regrid_int = -1;
    pp.query("regrid_int", regrid_int);
    if (regrid_int > 0) {
        const std::string key = pp_prefix + ".regrid_int";
        Abort(names.tag() + " The " + names.grid() + " on level " + std::to_string(lev)
              + " is built once from that level's boxes, so the refinement must not change, but "
              + key + " = " + std::to_string(regrid_int) + " regrids it. Use static refinement "
                "boxes with " + key + " = -1, or set " + names.key("anchor_level") + " = 0.");
    }
}

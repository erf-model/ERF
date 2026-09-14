/**
 * \file ERF_EBInit.cpp
 *
 * Construction of the embedded boundary (EB) geometry used by ERF, including
 * the terrain, plane, box and sphere implicit functions as well as the
 * STL-based building geometry.
 */

#include "ERF.H"

#include "AMReX_EB2_IF_Box.H"
#include "AMReX_EB2_IF_Sphere.H"
#include "AMReX_EB2_IF_Plane.H"
#include "AMReX_EB2_IF_Union.H"

#include "ERF_EBIFTerrain.H"
#include "ERF_EBIFBuildings.H"

using namespace amrex;

// Initialize embedded boundary geometries
void
ERF::initializeEB ()
{
    //
    // Construct the EB data structures and store in a separate class
    //
    std::string geometry ="terrain";
    ParmParse pp_eb2("eb2");
    pp_eb2.queryAdd("geometry", geometry);
    if ( solverChoice.terrain_type == TerrainType::EB ||
         solverChoice.terrain_type == TerrainType::ImmersedForcing)
    {
        constexpr int ngrow_for_eb = 4;  // This is the default in amrex but we need to explicitly pass it here since
                               // we want to also pass the build_coarse_level_by_coarsening argument
        const bool build_eb_for_multigrid = (solverChoice.terrain_type == TerrainType::EB &&
                                            ((solverChoice.project_initial_velocity[0] == 1) ||
                                            solverChoice.anelastic[0] == 1));
        // Note this just needs to be an integer > number of V-cycles one might use
        const int max_coarsening_level = (build_eb_for_multigrid) ? 100 : 0;
        const bool build_coarse_level_by_coarsening(false);

        // Define GeometryShop using the implicit function
        if (geometry == "terrain") {
            // Query building STL parameters upfront
            //
            // Note these live under the "erf" prefix, not "eb2", because "eb2" is amrex's own
            // namespace and amrex already reads eb2.stl_* for its own eb2.geometry = stl path.
            ParmParse pp_erf("erf");
            std::string buildings_stl_file;
            bool has_buildings_stl = pp_erf.query("buildings_stl_file", buildings_stl_file);

            Real stl_scale = 1.0;
            Array<Real,3> stl_center = {zero, zero, zero};
            int stl_reverse_normal = 0;

            if (has_buildings_stl) {
                pp_erf.query("buildings_stl_scale", stl_scale);
                pp_erf.query("buildings_stl_center", stl_center);
                pp_erf.query("buildings_stl_reverse_normal", stl_reverse_normal);
            }

            // Decide whether the terrain surface is part of the implicit function.  By default
            // we include it whenever the user has specified terrain by any of the means that
            // init_terrain_surface understands -- not just erf.terrain_file_name -- but the
            // user may also state the choice outright with erf.buildings_only.
            bool buildings_only = (has_buildings_stl && !prob->terrain_is_specified());

            bool buildings_only_in = false;
            if (pp_erf.query("buildings_only", buildings_only_in)) {
                if (buildings_only_in && !has_buildings_stl) {
                    Abort("erf.buildings_only is true but no erf.buildings_stl_file was given");
                }
                buildings_only = buildings_only_in;
            }

            // Determine geometry mode and log configuration
            std::string mode_description;
            if (buildings_only) {
                mode_description = "Buildings-only mode (no terrain)";
            } else if (has_buildings_stl) {
                mode_description = "Terrain + 3D buildings";
            } else {
                mode_description = "Terrain-only mode (no buildings)";
            }

            Print() << "Building EB geometry: " << mode_description << "\n";
            if (has_buildings_stl) {
                Print() << "  STL file: " << buildings_stl_file << "\n";
                Print() << "  STL scale: " << stl_scale << "\n";
                Print() << "  STL center: " << stl_center[0] << " " << stl_center[1] << " " << stl_center[2] << "\n";
                Print() << "  Reverse normals: " << stl_reverse_normal << "\n";
            }

            // Lambda to build EB geometry from shop
            auto build_eb = [&](auto const& gshop) {
                if (build_eb_for_multigrid) {
                    EB2::Build(gshop, geom[max_level], max_level, max_coarsening_level,
                                ngrow_for_eb, build_coarse_level_by_coarsening);
                } else {
                    EB2::Build(gshop, this->Geom(), ngrow_for_eb);
#if USE_FC_FACTORY
                    EB2::BuildFC();
#endif
                }
            };

            // Build the appropriate implicit function and geometry
            if (buildings_only) {
                // Buildings-only
                BuildingsIF buildings_if(buildings_stl_file, stl_scale, stl_center,
                                        stl_reverse_normal, geom[max_level]);
                auto gshop = EB2::makeShop(buildings_if);
                build_eb(gshop);

            } else {
                // Load terrain (from file or custom init_my_custom_terrain)
                Box terrain_bx(surroundingNodes(geom[max_level].Domain()));
                terrain_bx.grow(3);
                FArrayBox terrain_fab(makeSlab(terrain_bx,2,0),1);
                double dummy_time = 0.0;
                prob->init_terrain_surface(geom[max_level], terrain_fab, dummy_time);
                TerrainIF terrain_if(terrain_fab, geom[max_level], stretched_dz_d[max_level]);

                if (has_buildings_stl) {
                    // Terrain + Buildings
                    BuildingsIF buildings_if(buildings_stl_file, stl_scale, stl_center,
                                            stl_reverse_normal, geom[max_level]);
                    auto combined_if = EB2::makeUnion(terrain_if, buildings_if);
                    auto gshop = EB2::makeShop(combined_if);
                    build_eb(gshop);
                } else {
                    // Terrain-only
                    auto gshop = EB2::makeShop(terrain_if);
                    build_eb(gshop);
                }
            }

            Print() << "EB geometry built successfully: " << mode_description << ".\n";
        } else if (geometry == "plane") {
            RealArray plane_point{zero, zero, zero};
            RealArray plane_normal{zero, zero, -one}; // pointing into the solid region
            pp_eb2.queryAdd("plane_point", plane_point);
            pp_eb2.queryAdd("plane_normal", plane_normal);
            EB2::PlaneIF implicit_fun(plane_point, plane_normal, true);
            auto gshop = EB2::makeShop(implicit_fun);
            if (build_eb_for_multigrid) {
                EB2::Build(gshop, geom[max_level], max_level, max_coarsening_level,
                            ngrow_for_eb, build_coarse_level_by_coarsening);
            } else {
                EB2::Build(gshop, this->Geom(), ngrow_for_eb);
#if USE_FC_FACTORY
                EB2::BuildFC();
#endif
            }
        } else if (geometry == "box") {
            RealArray box_lo{zero, zero, zero};
            RealArray box_hi{zero, zero, zero};
            pp_eb2.queryAdd("box_lo", box_lo);
            pp_eb2.queryAdd("box_hi", box_hi);
            EB2::BoxIF implicit_fun(box_lo, box_hi, false);
            auto gshop = EB2::makeShop(implicit_fun);
            if (build_eb_for_multigrid) {
                EB2::Build(gshop, geom[max_level], max_level, max_coarsening_level,
                            ngrow_for_eb, build_coarse_level_by_coarsening);
            } else {
                EB2::Build(gshop, this->Geom(), ngrow_for_eb);
#if USE_FC_FACTORY
                EB2::BuildFC();
#endif
            }
        } else if (geometry == "sphere") {
            auto ProbLoArr = geom[max_level].ProbLoArray();
            auto ProbHiArr = geom[max_level].ProbHiArray();
            const Real xcen = myhalf * (ProbLoArr[0] + ProbHiArr[0]);
            const Real ycen = myhalf * (ProbLoArr[1] + ProbHiArr[1]);
            RealArray sphere_center = {xcen, ycen, zero};
            EB2::SphereIF implicit_fun(myhalf, sphere_center, false);
            auto gshop = EB2::makeShop(implicit_fun);
            if (build_eb_for_multigrid) {
                EB2::Build(gshop, geom[max_level], max_level, max_coarsening_level,
                            ngrow_for_eb, build_coarse_level_by_coarsening);
            } else {
                EB2::Build(gshop, this->Geom(), ngrow_for_eb);
#if USE_FC_FACTORY
                EB2::BuildFC();
#endif
            }
        }
    }

    if ( solverChoice.buildings_type == BuildingsType::ImmersedForcing) {
        constexpr int ngrow_for_eb = 4;
        if (geometry == "terrain") {
            Box buildings_bx(surroundingNodes(geom[max_level].Domain())); buildings_bx.grow(3);
            FArrayBox buildings_fab(makeSlab(buildings_bx,2,0),1);
            double dummy_time = 0.0;
            prob->init_buildings_surface(geom[max_level], buildings_fab, dummy_time);
            TerrainIF implicit_fun(buildings_fab, geom[max_level], stretched_dz_d[max_level]);
            auto gshop = EB2::makeShop(implicit_fun);
            EB2::Build(gshop, this->Geom(), ngrow_for_eb);
#if USE_FC_FACTORY
            EB2::BuildFC();
#endif
        } else if (geometry == "plane") {
            amrex::Abort("plane geometry is not supported with ImmersedForcing for buildings");
        } else if (geometry == "box") {
            RealArray box_lo{zero, zero, zero};
            RealArray box_hi{zero, zero, zero};
            pp_eb2.queryAdd("box_lo", box_lo);
            pp_eb2.queryAdd("box_hi", box_hi);
            EB2::BoxIF implicit_fun(box_lo, box_hi, false);
            auto gshop = EB2::makeShop(implicit_fun);
            EB2::Build(gshop, this->Geom(), ngrow_for_eb);
#if USE_FC_FACTORY
            EB2::BuildFC();
#endif
        } else if (geometry == "sphere") {
            amrex::Abort("sphere geometry is not supported with ImmersedForcing for buildings");
        }
    }
}

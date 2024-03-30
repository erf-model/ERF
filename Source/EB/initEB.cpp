#include <AMReX_ParmParse.H>
#include <AMReX_EB2.H>

#include <algorithm>
#include <ERF.H>
#include <FlowerIF.H>
#include <TerrainIF.H>
#include <prob_common.H>

using namespace amrex;

void ERF::MakeEBGeometry()
{
   /******************************************************************************
   * ERF.geometry=<string> specifies the EB geometry. <string> can be one of    *
   * box, cylinder, flower or terrain */

    ParmParse pp("eb2");

    std::string geom_type;
    pp.query("geometry", geom_type);

    /******************************************************************************
    *                                                                            *
    *  CONSTRUCT EB                                                              *
    *                                                                            *
    ******************************************************************************/

    int max_coarsening_level = 0;

    if(geom_type == "cylinder") {
        amrex::Print() << "\n Building cylinder geometry." << std::endl;
        make_eb_cylinder();

    } else if (geom_type == "terrain") {
        amrex::Print() << "\n Building EB geometry based on idealized terrain." << std::endl;
        Real dummy_time = 0.0;
        Box bx(surroundingNodes(Geom(0).Domain())); bx.grow(2);
        FArrayBox terrain_fab(makeSlab(bx,2,0),1);
        prob->init_custom_terrain(Geom(0), terrain_fab, dummy_time);
        TerrainIF ebterrain(terrain_fab, Geom(0));
        auto gshop = EB2::makeShop(ebterrain);
        EB2::Build(gshop, geom.back(), max_level, max_level+max_coarsening_level);

    } else if (geom_type == "flower") {
        FlowerIF flower(0.2, 0.1, 6, {AMREX_D_DECL(0.5,0.5,0.5)}, false);
        auto gshop = EB2::makeShop(flower);
        EB2::Build(gshop, geom.back(), max_level, max_level+max_coarsening_level);

    } else if(geom_type == "box") {
        amrex::Print() << "\n Building box geometry." << std::endl;
        make_eb_box();

    } else {

        amrex::Print() << "\n No EB geometry declared in inputs => "
                       << " Will build all regular geometry." << std::endl;
        make_eb_regular();
    }
    amrex::Print() << "Done making the geometry ebfactory.\n" << std::endl;
}

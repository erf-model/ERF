#include <AMReX_ParmParse.H>
#include <AMReX_EB2.H>

#include <algorithm>
#include <ERF.H>
#include <ERF_TerrainIF.H>
#include <ERF_ProbCommon.H>

using namespace amrex;

void ERF::MakeEBGeometry()
{
   /******************************************************************************
   * ERF.geometry=<string> specifies the EB geometry. <string> can be either     *
   * box or terrain */

    ParmParse pp("eb2");

    std::string geom_type = "terrain";
    // pp.query("geometry", geom_type);

    /******************************************************************************
    *                                                                            *
    *  CONSTRUCT EB                                                              *
    *                                                                            *
    ******************************************************************************/

    int lev = 0;

    int max_coarsening_level;
    if (solverChoice.anelastic[0] == 1) {
        max_coarsening_level = 100;
    } else {
        max_coarsening_level = 0;
    }

    if (geom_type == "terrain") {
        amrex::Print() << "\n Building EB geometry based on idealized terrain." << std::endl;
        Real dummy_time = 0.0;

        Box bx(surroundingNodes(Geom(lev).Domain())); bx.grow(3);
        FArrayBox terrain_fab(makeSlab(bx,2,0),1);

        prob->init_terrain_surface(Geom(lev), terrain_fab, dummy_time);

        TerrainIF ebterrain(terrain_fab, Geom(lev), stretched_dz_d[lev]);

        auto gshop = EB2::makeShop(ebterrain);

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

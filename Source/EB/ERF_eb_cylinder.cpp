#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <ERF.H>

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Function to create a simple cylinder EB.                                     *
 *                                                                              *
 ********************************************************************************/
void ERF::make_eb_cylinder()
{
    // Initialise cylinder parameters
    bool inside = true;
    Real radius = 0.0002;
    int direction = 0;
    Vector<Real> centervec(3);
    Real rotation  = 0;
    int rotation_axe  = 0;

    // Get cylinder information from inputs file.                               *
    ParmParse pp("cylinder");

    pp.query("internal_flow", inside);
    pp.query("radius", radius);
    pp.query("direction", direction);
    pp.query("rotation",   rotation);
    pp.query("rotation_axe",   rotation_axe);
    pp.getarr("center", centervec, 0, 3);
    Array<Real, AMREX_SPACEDIM> center = {AMREX_D_DECL(centervec[0], centervec[1], centervec[2])};

    rotation = (rotation/180.)*M_PI;

    // Print info about cylinder
    amrex::Print() << " " << std::endl;
    amrex::Print() << " Internal Flow: " << inside << std::endl;
    amrex::Print() << " Radius:    " << radius << std::endl;
    amrex::Print() << " Direction: " << direction << std::endl;
    amrex::Print() << " Rotation angle(rad): " << rotation << std::endl;
    amrex::Print() << " Rotation axe: " << rotation_axe << std::endl;
#if (AMREX_SPACDEIM == 3)
    amrex::Print() << " Center:    " << center[0] << ", " << center[1] << ", " << center[2]
                   << std::endl;
#else
    amrex::Print() << " Center:    " << center[0] << ", " << center[1] << std::endl;
#endif

    // Build the Cylinder implficit function representing the curved walls
    EB2::CylinderIF my_cyl(radius, direction, center, inside);

    auto my_cyl_rot = EB2::rotate(my_cyl, rotation, rotation_axe);

    // Generate GeometryShop
    auto gshop = EB2::makeShop(my_cyl_rot);

    // Build index space
    int max_level_here = 0;
    int max_coarsening_level = 100;
    EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
}

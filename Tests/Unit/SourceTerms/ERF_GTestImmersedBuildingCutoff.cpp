#include <gtest/gtest.h>

#include <AMReX_Geometry.H>

#include "ERF_DataStruct.H"

TEST(ImmersedBuildingCutoff, AllowsNominalFiftyMetreSpacingWithRoundoff)
{
    using namespace amrex;
    SolverChoice choice;
    choice.buildings_type = BuildingsType::ImmersedForcing;

    const Box domain(IntVect(0,0,0), IntVect(3,3,0));
    auto geometry_with_spacing = [&] (Real spacing) {
        const Real extent = Real(4) * spacing;
        const RealBox bounds({Real(0), Real(0), Real(0)},
                             {extent, extent, Real(1)});
        return Geometry(domain, bounds, CoordSys::cartesian, {0,0,0});
    };

    EXPECT_TRUE(choice.use_building_immersed_forcing(
        geometry_with_spacing(Real(50) * (one + Real(5.e-7)))));
    EXPECT_FALSE(choice.use_building_immersed_forcing(
        geometry_with_spacing(Real(50) * (one + Real(3.e-6)))));

    choice.buildings_type = BuildingsType::None;
    EXPECT_FALSE(choice.use_building_immersed_forcing(
        geometry_with_spacing(Real(50))));
}

#include <cstdlib>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

#include "ERF_GTestIceDepositionCommon.H"

using namespace ice_deposition_test;

// Motivation: ERF's ice habit is set by the Chen-Lamb inherent growth ratio
// Gamma(T). Takahashi et al. (1991) show plates (phi=c/a < 1) near -15 C and
// columns (phi > 1) near -6 C. This grows a single crystal for 10 min at
// constant T and water saturation through the production deposition routines and
// checks that the predicted habit follows that signature.
TEST(IceDepositionScalar, HabitTransitionsMatchTakahashi)
{
    const CrystalState plate  = grow_crystal(amrex::Real(273.15 - 15.0));  // TH91: thin plate
    const CrystalState column = grow_crystal(amrex::Real(273.15 -  6.0));  // TH91: column
    const CrystalState cold   = grow_crystal(amrex::Real(273.15 - 23.0));  // TH91: near-equant

    // plate at -15 C: c << a
    EXPECT_LT(aspect_ratio(plate), 1.0) << "phi(-15C)=" << aspect_ratio(plate);
    // column at -6 C: c > a
    EXPECT_GT(aspect_ratio(column), 1.0) << "phi(-6C)=" << aspect_ratio(column);
    // the plate is much flatter than the column is elongated/than the cold crystal
    EXPECT_LT(aspect_ratio(plate), aspect_ratio(cold));
    EXPECT_LT(aspect_ratio(plate), aspect_ratio(column));

    // the crystal actually grew (mass and size increased from the 2 um seed)
    EXPECT_GT(plate.a, amrex::Real(2.0e-6));
    EXPECT_GT(mass_ug(plate), amrex::Real(0.0));
}

// Motivation: the dendrite regime near -15 C is the mass and a-axis maximum in
// TH91; the crystal there should outgrow the colder -23 C crystal.
TEST(IceDepositionScalar, DendriteRegimeGrowsLargest)
{
    const CrystalState d15 = grow_crystal(amrex::Real(273.15 - 15.0));
    const CrystalState d23 = grow_crystal(amrex::Real(273.15 - 23.0));
    EXPECT_GT(d15.a, d23.a);
    EXPECT_GT(mass_ug(d15), mass_ug(d23));
}

// Motivation: deposition and sublimation are the two branches of the same
// production routine (dMdt::rhs_func), whose mass rate is proportional to S-1.
// The crystal grows when the air is supersaturated over ice and sublimates
// (loses mass) when subsaturated, with magnitude scaling with |S-1|.
TEST(IceDepositionScalar, SublimationShrinksAtSubsaturation)
{
    const amrex::Real T = amrex::Real(273.15 - 15.0);
    const amrex::Real dep = deposition_mass_rate(T, amrex::Real(1.05));
    EXPECT_GT(dep, amrex::Real(0.0));                                   // deposition
    EXPECT_LT(deposition_mass_rate(T, amrex::Real(0.95)), amrex::Real(0.0));  // sublimation
    EXPECT_NEAR(deposition_mass_rate(T, amrex::Real(1.0)), amrex::Real(0.0),
                std::abs(dep) * amrex::Real(1.0e-3));                   // no exchange at S = 1
    EXPECT_GT(std::abs(deposition_mass_rate(T, amrex::Real(0.90))),
              std::abs(deposition_mass_rate(T, amrex::Real(0.95))));    // |rate| grows with |S-1|
}

// Motivation: emit mass, a-, c-axis, aspect ratio and apparent density vs T from
// the production routines so they can be plotted against TH91 (Welss Fig.10).
// Skipped unless ERF_ICEDEP_CURVES names an output CSV (no-op in CI).
TEST(IceDepositionScalar, DumpCurvesForPlotting)
{
    const char* out = std::getenv("ERF_ICEDEP_CURVES");
    if (out == nullptr || out[0] == '\0') {
        GTEST_SKIP() << "set ERF_ICEDEP_CURVES=<file.csv> to dump the growth curves";
    }
    std::ofstream f(out);
    ASSERT_TRUE(f.is_open()) << "cannot open " << out;
    f << "T_C,mass_ug,a_um,c_um,phi,rho\n";
    f.precision(8);
    for (int tc = -23; tc <= -3; ++tc) {
        const CrystalState x = grow_crystal(amrex::Real(273.15 + tc));
        const amrex::Real rho = x.mass / (four_thirds_pi * x.a * x.a * x.c);
        f << tc << ',' << mass_ug(x) << ',' << x.a * 1e6 << ',' << x.c * 1e6 << ','
          << aspect_ratio(x) << ',' << rho << '\n';
    }
    f.close();
    EXPECT_TRUE(f.good());
}

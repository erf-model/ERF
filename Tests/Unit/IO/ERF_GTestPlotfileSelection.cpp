#include "ERF_PlotfileSelection.H"

#include <gtest/gtest.h>

using erf_plotfile::Plot3DSelectionCapabilities;

TEST(Plotfile3DSelection, ConservedComponentsFollowActiveState)
{
    const Plot3DSelectionCapabilities kessler{
        7, 3, 0, 1, MoistureType::Kessler};

    EXPECT_EQ(erf_plotfile::plot3d_conserved_component_index("rhoQ1"), 4);
    EXPECT_EQ(erf_plotfile::plot3d_conserved_component_index("rhoQ3"), 6);
    EXPECT_EQ(erf_plotfile::plot3d_conserved_component_index("rhoQ4"), 7);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rhoQ3", kessler));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rhoQ4", kessler));
}

TEST(Plotfile3DSelection, MoistureDiagnosticsMatchCapabilities)
{
    const Plot3DSelectionCapabilities kessler{
        7, 3, 0, 1, MoistureType::Kessler};

    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qv", kessler));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qc", kessler));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qrain", kessler));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qt", kessler));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qp", kessler));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qi", kessler));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qsnow", kessler));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rel_humidity", kessler));

    const Plot3DSelectionCapabilities superdroplets{
        7, 3, 0, 7, MoistureType::SuperDroplets};
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rel_humidity", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("condensation_rate", superdroplets));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("snow_accum", superdroplets));
}

TEST(Plotfile3DSelection, PressureReadersShareAllocationPredicate)
{
    EXPECT_TRUE(erf_plotfile::plot3d_needs_pressure({"temp", "VPD"}));
    EXPECT_TRUE(erf_plotfile::plot3d_needs_pressure({"pressure"}));
    EXPECT_TRUE(erf_plotfile::plot3d_needs_pressure({"eq_pot_temp"}));
    EXPECT_TRUE(erf_plotfile::plot3d_needs_pressure({"qsat"}));
    EXPECT_TRUE(erf_plotfile::plot3d_needs_pressure({"dpdx"}));
#ifdef ERF_COMPUTE_ERROR
    EXPECT_TRUE(erf_plotfile::plot3d_needs_pressure({"pp_err"}));
#endif
}

TEST(Plotfile3DSelection, ParticleCountsIncludeConfiguredUnallocatedContainers)
{
    const amrex::Vector<std::string> requested{
        "aerosols_count", "tracer_particles_count"};
    const amrex::Vector<std::string> configured{
        "tracer_particles", "aerosols"};

    const auto selected = erf_plotfile::plot3d_selected_particle_count_names(requested, configured);
    ASSERT_EQ(selected.size(), 2);
    EXPECT_EQ(selected[0], "aerosols");
    EXPECT_EQ(selected[1], "tracer_particles");
}

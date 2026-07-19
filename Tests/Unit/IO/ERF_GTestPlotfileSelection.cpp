#include "ERF_PlotfileSelection.H"

#include <array>

#include <gtest/gtest.h>

using erf_plotfile::Plot3DSelectionCapabilities;

namespace {

Plot3DSelectionCapabilities
make_capabilities (const MoistureType moisture_type,
                   const int conserved_state_size = 32,
                   const int moist_state_size = 11,
                   const int moist_numconc_size = 0,
                   const int qmoist_size = 7)
{
    Plot3DSelectionCapabilities caps;
    caps.conserved_state_size = conserved_state_size;
    caps.moist_state_size = moist_state_size;
    caps.moist_numconc_size = moist_numconc_size;
    caps.qmoist_size = qmoist_size;
    caps.moisture_type = moisture_type;
    caps.moisture = erf_plotfile::plot3d_moisture_capabilities(moisture_type);
    return caps;
}

} // namespace

TEST(Plotfile3DSelection, MoistureCapabilityTruthTable)
{
    struct Expected {
        MoistureType type;
        bool vapor;
        bool cloud_liquid;
        bool cloud_ice;
        bool rain;
        bool snow;
        bool graupel;
        bool cloud_number;
        bool ice_number;
        bool rain_number;
        bool snow_number;
        bool graupel_number;
        bool rain_accumulation;
        bool snow_accumulation;
        bool graupel_accumulation;
    };

    const std::array<Expected, 12> expected{{
        {MoistureType::None, false, false, false, false, false, false,
         false, false, false, false, false, false, false, false},
        {MoistureType::MoistNoCondensation, true, false, false, false, false, false,
         false, false, false, false, false, false, false, false},
        {MoistureType::SatAdj, true, true, false, false, false, false,
         false, false, false, false, false, false, false, false},
        {MoistureType::Kessler_NoRain, true, true, false, false, false, false,
         false, false, false, false, false, false, false, false},
        {MoistureType::Kessler, true, true, false, true, false, false,
         false, false, false, false, false, true, false, false},
        {MoistureType::SAM_NoPrecip_NoIce, true, true, false, false, false, false,
         false, false, false, false, false, false, false, false},
        {MoistureType::SAM_NoIce, true, true, false, true, false, false,
         false, false, false, false, false, true, false, false},
        {MoistureType::SAM, true, true, true, true, true, true,
         false, false, false, false, false, true, true, true},
        {MoistureType::Morrison_NoIce, true, true, false, true, false, false,
         true, false, true, false, false, true, false, false},
        {MoistureType::Morrison, true, true, true, true, true, true,
         true, true, true, true, true, true, true, true},
        {MoistureType::WSM6, true, true, true, true, true, true,
         false, false, false, false, false, true, true, true},
        {MoistureType::SuperDroplets, true, true, false, false, false, false,
         false, false, false, false, false, true, false, false}
    }};

    for (const auto& item : expected) {
        const auto caps = erf_plotfile::plot3d_moisture_capabilities(item.type);
        EXPECT_EQ(caps.vapor, item.vapor);
        EXPECT_EQ(caps.cloud_liquid, item.cloud_liquid);
        EXPECT_EQ(caps.cloud_ice, item.cloud_ice);
        EXPECT_EQ(caps.rain, item.rain);
        EXPECT_EQ(caps.snow, item.snow);
        EXPECT_EQ(caps.graupel, item.graupel);
        EXPECT_EQ(caps.cloud_number, item.cloud_number);
        EXPECT_EQ(caps.ice_number, item.ice_number);
        EXPECT_EQ(caps.rain_number, item.rain_number);
        EXPECT_EQ(caps.snow_number, item.snow_number);
        EXPECT_EQ(caps.graupel_number, item.graupel_number);
        EXPECT_EQ(caps.rain_accumulation, item.rain_accumulation);
        EXPECT_EQ(caps.snow_accumulation, item.snow_accumulation);
        EXPECT_EQ(caps.graupel_accumulation, item.graupel_accumulation);
    }
}

TEST(Plotfile3DSelection, ConservedComponentsRequireActiveBounds)
{
    const auto caps = make_capabilities(MoistureType::Kessler, 7, 3, 0, 1);

    EXPECT_EQ(erf_plotfile::plot3d_conserved_component_index("rhoQ1"), 4);
    EXPECT_EQ(erf_plotfile::plot3d_conserved_component_index("rhoQ3"), 6);
    EXPECT_EQ(erf_plotfile::plot3d_conserved_component_index("rhoQ4"), 7);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rhoQ3", caps));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rhoQ4", caps));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rhoQ99", caps));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("not_a_fixed_field", caps));
}

TEST(Plotfile3DSelection, SemanticAndStorageChecksAreIndependent)
{
    const auto no_rain = make_capabilities(MoistureType::Kessler_NoRain, 32, 3, 0, 1);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qv", no_rain));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qc", no_rain));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qrain", no_rain));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qp", no_rain));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rain_accum", no_rain));

    const auto no_precip = make_capabilities(MoistureType::SAM_NoPrecip_NoIce, 32, 6, 0, 3);
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qrain", no_precip));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rain_accum", no_precip));

    const auto no_ice = make_capabilities(MoistureType::Morrison_NoIce, 32, 11, 5, 3);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qrain", no_ice));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qi", no_ice));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("nc", no_ice));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("nr", no_ice));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("ni", no_ice));

    const auto too_narrow = make_capabilities(MoistureType::Kessler, 7, 1, 0, 1);
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qrain", too_narrow));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("moist_density", too_narrow));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rain_accum", too_narrow));

    const auto superdroplets = make_capabilities(MoistureType::SuperDroplets, 32, 3, 0, 7);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rel_humidity", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("condensation_rate", superdroplets));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qrain", superdroplets));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qp", superdroplets));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("snow_accum", superdroplets));
}

TEST(Plotfile3DSelection, OptionalStorageGroupsAreExplicit)
{
    auto caps = make_capabilities(MoistureType::None, 4, 0, 0, 0);
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("u_t_avg", caps));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qsrc_sw", caps));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("Kmv", caps));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("diss", caps));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("walldist", caps));

    caps.time_average_storage = true;
    caps.radiation_heating_storage = true;
    caps.eddy_diffusivity_storage = true;
    caps.dissipation_storage = true;
    caps.wall_distance_storage = true;
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("u_t_avg", caps));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qsrc_lw", caps));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("Lturb", caps));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("diss", caps));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("walldist", caps));
}

TEST(Plotfile3DSelection, PressureReadersShareAllocationPredicate)
{
    EXPECT_TRUE(erf_plotfile::plot3d_needs_pressure({"temp", "VPD"}));
    EXPECT_TRUE(erf_plotfile::plot3d_needs_pressure({"pressure"}));
    EXPECT_TRUE(erf_plotfile::plot3d_needs_pressure({"eq_pot_temp"}));
    EXPECT_TRUE(erf_plotfile::plot3d_needs_pressure({"qsat"}));
    EXPECT_TRUE(erf_plotfile::plot3d_needs_pressure({"dpdx"}));
    EXPECT_FALSE(erf_plotfile::plot3d_needs_pressure({"theta"}));
#ifdef ERF_COMPUTE_ERROR
    EXPECT_TRUE(erf_plotfile::plot3d_needs_pressure({"pp_err"}));
#endif
}

TEST(Plotfile3DSelection, ParticleCountsHandleUnknownAndDuplicateNames)
{
    const amrex::Vector<std::string> requested{
        "aerosols_count", "unknown_count", "aerosols_count", "not_a_count"};
    const amrex::Vector<std::string> configured{"tracer_particles", "aerosols"};

    const auto selected = erf_plotfile::plot3d_selected_particle_count_names(requested, configured);
    ASSERT_EQ(selected.size(), 1);
    EXPECT_EQ(selected[0], "aerosols");
}

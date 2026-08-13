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

// Motivation: Keep selector and writer range formulas identical so aggregate
// fields cannot pass selection with an incomplete source range.
TEST(Plotfile3DSelection, AggregateRangesMatchFixedWriterLayouts)
{
    struct Expected {
        int moist_state_size;
        int moist_numconc_size;
        int total_first;
        int total_last;
        int nonprecip_first;
        int nonprecip_last;
        int precip_first;
        int precip_last;
    };

    const std::array<Expected, 7> expected{{
        {2, 0, 1, 2, 1, 2, 0, -1},
        {3, 0, 1, 3, 1, 2, 3, 3},
        {6, 0, 1, 6, 1, 3, 4, 6},
        {11, 5, 1, 6, 1, 3, 4, 6},
        {3, 0, 1, 3, 1, 2, 3, 3},
        {3, 4, 0, -1, 1, 2, 0, -1},
        {0, 0, 0, -1, 0, -1, 0, -1}
    }};

    for (const auto& item : expected) {
        const auto total = erf_plotfile::plot3d_total_mass_q_range(
            item.moist_state_size, item.moist_numconc_size);
        const auto nonprecip = erf_plotfile::plot3d_nonprecipitating_q_range(
            item.moist_state_size);
        const auto precip = erf_plotfile::plot3d_precipitating_q_range(
            item.moist_state_size, item.moist_numconc_size);
        EXPECT_EQ(total.first_q, item.total_first);
        EXPECT_EQ(total.last_q, item.total_last);
        EXPECT_EQ(nonprecip.first_q, item.nonprecip_first);
        EXPECT_EQ(nonprecip.last_q, item.nonprecip_last);
        EXPECT_EQ(precip.first_q, item.precip_first);
        EXPECT_EQ(precip.last_q, item.precip_last);
        EXPECT_EQ(total.valid(), item.total_last >= item.total_first && item.total_first >= 1);
        EXPECT_EQ(nonprecip.valid(), item.nonprecip_last >= item.nonprecip_first && item.nonprecip_first >= 1);
        EXPECT_EQ(precip.valid(), item.precip_last >= item.precip_first && item.precip_first >= 1);
    }
}

// Motivation: Cloud ice is stored in q3; its availability must not depend
// on the unrelated q4 rain component. The synthetic truncated state isolates
// selector source mapping and is not a supported production SAM layout.
TEST(Plotfile3DSelection, CloudIceUsesQ3SourceComponent)
{
    const int through_q3 =
        erf_plotfile::plot3d_q_conserved_component_index(3) + 1;

    const auto cold_through_q3 = make_capabilities(
        MoistureType::SAM,
        through_q3,
        6,
        0,
        3);

    EXPECT_TRUE(
        erf_plotfile::plot3d_source_component_available(
            cold_through_q3, 3));
    EXPECT_FALSE(
        erf_plotfile::plot3d_source_component_available(
            cold_through_q3, 4));

    EXPECT_TRUE(
        erf_plotfile::plot3d_fixed_variable_available(
            "qi", cold_through_q3));
    EXPECT_FALSE(
        erf_plotfile::plot3d_fixed_variable_available(
            "qrain", cold_through_q3));
}

// Motivation: SuperDroplets replaces its constructor sentinel during
// readInputs() and stores qv, qc, and qr in RhoQ1-RhoQ3. This regression
// prevents the temporary sentinel from being mistaken for a missing q3.
TEST(Plotfile3DSelection, SuperDropletsWaterStateIncludesRainComponent)
{
    const int through_q3 =
        erf_plotfile::plot3d_q_conserved_component_index(3) + 1;
    const auto superdroplets = make_capabilities(
        MoistureType::SuperDroplets, through_q3, 3, 0, 7);

    EXPECT_EQ(erf_plotfile::plot3d_q_conserved_component_index(1), 4);
    EXPECT_EQ(erf_plotfile::plot3d_q_conserved_component_index(2), 5);
    EXPECT_EQ(erf_plotfile::plot3d_q_conserved_component_index(3), 6);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qv", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qc", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qrain", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qt", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qn", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qp", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("moist_density", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_q_range_available(
        superdroplets, erf_plotfile::plot3d_total_mass_q_range(3, 0)));
    EXPECT_TRUE(erf_plotfile::plot3d_q_range_available(
        superdroplets, erf_plotfile::plot3d_nonprecipitating_q_range(3)));
    EXPECT_TRUE(erf_plotfile::plot3d_q_range_available(
        superdroplets, erf_plotfile::plot3d_precipitating_q_range(3, 0)));
}

// Motivation: Correcting the production SuperDroplets width must not weaken
// the generic rule that qrain, qt, and qp require a complete q3 source.
TEST(Plotfile3DSelection, TruncatedSuperDropletsStateRejectsQ3Fields)
{
    const int through_q2 =
        erf_plotfile::plot3d_q_conserved_component_index(2) + 1;
    const auto truncated_superdroplets = make_capabilities(
        MoistureType::SuperDroplets, through_q2, 3, 0, 7);

    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qv", truncated_superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qc", truncated_superdroplets));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qrain", truncated_superdroplets));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qt", truncated_superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qn", truncated_superdroplets));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qp", truncated_superdroplets));
}

// Motivation: Bounds hardening must not remove valid aggregate diagnostics
// from ordinary warm-rain and mixed-phase schemes.
TEST(Plotfile3DSelection, StandardAggregateRangesRemainAvailable)
{
    const auto kessler = make_capabilities(MoistureType::Kessler, 7, 3, 0, 1);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qt", kessler));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qn", kessler));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qp", kessler));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("moist_density", kessler));

    const auto morrison = make_capabilities(MoistureType::Morrison, 15, 11, 5, 3);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qt", morrison));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qn", morrison));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qp", morrison));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("moist_density", morrison));

    const auto truncated = make_capabilities(MoistureType::Kessler, 6, 3, 0, 1);
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qt", truncated));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qp", truncated));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qn", truncated));
}

// Motivation: Invalid metadata must fail closed before a writer range is
// constructed from negative or inconsistent sizes.
TEST(Plotfile3DSelection, MalformedAggregateSizesAreRejected)
{
    EXPECT_FALSE(erf_plotfile::plot3d_total_mass_q_range(3, -1).valid());
    EXPECT_FALSE(erf_plotfile::plot3d_total_mass_q_range(3, 4).valid());
    EXPECT_FALSE(erf_plotfile::plot3d_total_mass_q_range(0, 0).valid());
    EXPECT_FALSE(erf_plotfile::plot3d_precipitating_q_range(3, -1).valid());
    EXPECT_FALSE(erf_plotfile::plot3d_precipitating_q_range(3, 4).valid());

    const auto malformed = make_capabilities(MoistureType::Morrison, 15, 3, 4, 3);
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qt", malformed));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qp", malformed));
}

// Motivation: Physical species availability depends on the selected scheme,
// not merely on allocated state width or placeholder components.
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
        {MoistureType::MoistNoCondensation, true, true, false, false, false, false,
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
        {MoistureType::SuperDroplets, true, true, true, true, true, true,
         false, false, false, false, false, true, true, false}
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

// Motivation: NoCondensation still stores vapor and cloud water in RhoQ1 and
// RhoQ2.  A diagnosed zero in the cloud-water field is therefore distinct
// from an unavailable SHOC diagnostic, which is represented by ERF's -999
// plotfile sentinel before the native SHOC lifecycle has run.
TEST(Plotfile3DSelection, NoCondensationCloudWaterIsSelectable)
{
    const auto no_condensation = make_capabilities(
        MoistureType::MoistNoCondensation, 6, 2, 0, 1);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qv", no_condensation));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qc", no_condensation));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qrain", no_condensation));
}

// Motivation: Raw rhoQ requests must never index beyond the active conserved
// state even though the public inventory lists rhoQ1 through rhoQ11.
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

// Motivation: A field is selectable only when both its physical species and
// its concrete source storage exist.
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

    const int through_q3 =
        erf_plotfile::plot3d_q_conserved_component_index(3) + 1;
    const auto superdroplets = make_capabilities(
        MoistureType::SuperDroplets, through_q3, 3, 0, 7);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qt", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qrain", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qp", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rel_humidity", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("condensation_rate", superdroplets));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("snow_accum", superdroplets));
    // graupel accumulation is rejected on semantic grounds alone; the storage
    // is wide enough for it
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("graup_accum", superdroplets));
}

// Motivation: Optional diagnostics must be rejected before the writer can
// dereference unallocated configuration-dependent storage.
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

// Motivation: Every diagnostic that reads pressure must participate in the
// same allocation predicate, including standalone VPD and conditional pp_err.
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

// Motivation: Configured particle containers may be unallocated at output
// time, but unknown and duplicate count requests must not create extra fields.
TEST(Plotfile3DSelection, ParticleCountsHandleUnknownAndDuplicateNames)
{
    const amrex::Vector<std::string> requested{
        "aerosols_count", "unknown_count", "aerosols_count", "not_a_count"};
    const amrex::Vector<std::string> configured{"tracer_particles", "aerosols"};

    const auto selected = erf_plotfile::plot3d_selected_particle_count_names(requested, configured);
    ASSERT_EQ(selected.size(), 1);
    EXPECT_EQ(selected[0], "aerosols");
}

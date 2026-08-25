#include "ERF_PlotfileSelection.H"

#include <array>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using erf_plotfile::Plot3DSelectionCapabilities;

namespace {

// The production conserved-state width for a scheme that carries all eleven
// moist components.  Only the dry names are bounded by this: moisture names are
// selected from the index map instead.
constexpr int full_state_size = NDRY + NSCALARS + 11;

Plot3DSelectionCapabilities
make_capabilities (const MoistureType moisture_type,
                   const int conserved_state_size = full_state_size)
{
    Plot3DSelectionCapabilities caps;
    caps.conserved_state_size = conserved_state_size;
    caps.moisture_indices = MoistureComponentIndices::from_moisture_model(moisture_type);
    return caps;
}

// Every moisture output variable named in ERF.H's derived_names, plus the raw
// moist state names.  Used to sweep both output paths over the same list.
const std::vector<std::string>&
moisture_variable_names ()
{
    static const std::vector<std::string> names{
        "moist_density", "qv", "qc", "qi", "qrain", "qsnow", "qgraup",
        "qt", "qn", "qp", "qsat", "nc", "ni", "nr", "ns", "ng", "nn",
        "rain_accum", "snow_accum", "graup_accum",
        "rel_humidity", "condensation_rate", "precipitable"};
    return names;
}

const std::array<MoistureType, 13>&
all_moisture_types ()
{
    static const std::array<MoistureType, 13> types{
        MoistureType::None, MoistureType::MoistNoCondensation, MoistureType::SatAdj,
        MoistureType::Kessler_NoRain, MoistureType::Kessler,
        MoistureType::SAM_NoPrecip_NoIce, MoistureType::SAM_NoIce, MoistureType::SAM,
        MoistureType::Morrison_NoIce, MoistureType::Morrison,
        MoistureType::WSM6, MoistureType::WDM6, MoistureType::SuperDroplets};
    return types;
}

} // namespace

// Motivation: The index map is now the only description of a scheme's moist
// state, so its conserved components must match the slots each scheme actually
// reads and writes in Copy_State_to_Micro / Copy_Micro_to_State.
TEST(Plotfile3DSelection, MoistureIndexMapMatchesSchemeStateLayouts)
{
    struct Expected {
        MoistureType type;
        int qv; int qc; int qi; int qr; int qs; int qg;
        int nc; int ni; int nr; int ns; int ng; int nn;
    };
    constexpr int no = MoistureComponentIndices::absent;

    const std::array<Expected, 13> expected{{
        {MoistureType::None, no, no, no, no, no, no, no, no, no, no, no, no},
        {MoistureType::MoistNoCondensation, RhoQ1_comp, RhoQ2_comp, no, no, no, no,
         no, no, no, no, no, no},
        {MoistureType::SatAdj, RhoQ1_comp, RhoQ2_comp, no, no, no, no,
         no, no, no, no, no, no},
        {MoistureType::Kessler_NoRain, RhoQ1_comp, RhoQ2_comp, no, no, no, no,
         no, no, no, no, no, no},
        {MoistureType::Kessler, RhoQ1_comp, RhoQ2_comp, no, RhoQ3_comp, no, no,
         no, no, no, no, no, no},
        {MoistureType::SAM_NoPrecip_NoIce, RhoQ1_comp, RhoQ2_comp, no, no, no, no,
         no, no, no, no, no, no},
        // The SAM class always allocates six moist components, so rain stays in
        // the fourth slot even with the frozen species switched off.
        {MoistureType::SAM_NoIce, RhoQ1_comp, RhoQ2_comp, no, RhoQ4_comp, no, no,
         no, no, no, no, no, no},
        {MoistureType::SAM, RhoQ1_comp, RhoQ2_comp, RhoQ3_comp, RhoQ4_comp,
         RhoQ5_comp, RhoQ6_comp, no, no, no, no, no, no},
        // Liquid-only Morrison is two-moment: cloud and rain number are
        // integrated even though the frozen species are not.
        {MoistureType::Morrison_NoIce, RhoQ1_comp, RhoQ2_comp, no, RhoQ4_comp, no, no,
         RhoQ7_comp, no, RhoQ9_comp, no, no, no},
        {MoistureType::Morrison, RhoQ1_comp, RhoQ2_comp, RhoQ3_comp, RhoQ4_comp,
         RhoQ5_comp, RhoQ6_comp, RhoQ7_comp, RhoQ8_comp, RhoQ9_comp,
         RhoQ10_comp, RhoQ11_comp, no},
        {MoistureType::WSM6, RhoQ1_comp, RhoQ2_comp, RhoQ3_comp, RhoQ4_comp,
         RhoQ5_comp, RhoQ6_comp, no, no, no, no, no, no},
        // WDM6 puts the aerosol reservoir in the slot Morrison uses for ice
        // number, and carries no snow or graupel number.
        {MoistureType::WDM6, RhoQ1_comp, RhoQ2_comp, RhoQ3_comp, RhoQ4_comp,
         RhoQ5_comp, RhoQ6_comp, RhoQ7_comp, no, RhoQ9_comp, no, no, RhoQ8_comp},
        {MoistureType::SuperDroplets, RhoQ1_comp, RhoQ2_comp, RhoQ3_comp, RhoQ4_comp,
         RhoQ5_comp, RhoQ6_comp, no, no, no, no, no, no}
    }};

    for (const auto& item : expected) {
        const auto mi = MoistureComponentIndices::from_moisture_model(item.type);
        EXPECT_EQ(mi.qv, item.qv);
        EXPECT_EQ(mi.qc, item.qc);
        EXPECT_EQ(mi.qi, item.qi);
        EXPECT_EQ(mi.qr, item.qr);
        EXPECT_EQ(mi.qs, item.qs);
        EXPECT_EQ(mi.qg, item.qg);
        EXPECT_EQ(mi.nc, item.nc);
        EXPECT_EQ(mi.ni, item.ni);
        EXPECT_EQ(mi.nr, item.nr);
        EXPECT_EQ(mi.ns, item.ns);
        EXPECT_EQ(mi.ng, item.ng);
        EXPECT_EQ(mi.nn, item.nn);
        EXPECT_EQ(mi.has_moisture(), item.type != MoistureType::None);
    }
}

// Motivation: The accumulation and moist-diagnostic slots differ per scheme, so
// a wrong slot silently writes another field's data under the right name.
TEST(Plotfile3DSelection, QmoistDiagnosticSlotsMatchSchemeMicVarMaps)
{
    constexpr int no = MoistureComponentIndices::absent;
    constexpr int derived = MoistureComponentIndices::computed_from_state;

    struct Expected {
        MoistureType type;
        int rain_accum; int snow_accum; int graup_accum; int rel_hum; int cond_rate;
    };

    const std::array<Expected, 13> expected{{
        {MoistureType::None, no, no, no, no, no},
        {MoistureType::MoistNoCondensation, no, no, no, no, no},
        // SatAdj publishes no qmoist arrays, so relative humidity is recovered
        // from the conserved state at output time.
        {MoistureType::SatAdj, no, no, no, derived, no},
        {MoistureType::Kessler_NoRain, no, no, no, no, no},
        {MoistureType::Kessler, 0, no, no, no, no},
        {MoistureType::SAM_NoPrecip_NoIce, no, no, no, no, no},
        {MoistureType::SAM_NoIce, 0, no, no, no, no},
        {MoistureType::SAM, 0, 1, 2, no, no},
        {MoistureType::Morrison_NoIce, 0, no, no, no, no},
        {MoistureType::Morrison, 0, 1, 2, no, no},
        {MoistureType::WSM6, 0, 1, 2, no, no},
        {MoistureType::WDM6, 0, 1, 2, no, no},
        // SuperDroplets has its own qmoist layout, and never fills the graupel
        // accumulation slot it allocates.
        {MoistureType::SuperDroplets, 8, 10, no, 7, 3}
    }};

    for (const auto& item : expected) {
        const auto mi = MoistureComponentIndices::from_moisture_model(item.type);
        EXPECT_EQ(mi.rain_accum, item.rain_accum);
        EXPECT_EQ(mi.snow_accum, item.snow_accum);
        EXPECT_EQ(mi.graup_accum, item.graup_accum);
        EXPECT_EQ(mi.rel_hum, item.rel_hum);
        EXPECT_EQ(mi.cond_rate, item.cond_rate);

        EXPECT_EQ(mi.qmoist_index_for_var("rain_accum"), item.rain_accum);
        EXPECT_EQ(mi.qmoist_index_for_var("snow_accum"), item.snow_accum);
        EXPECT_EQ(mi.qmoist_index_for_var("graup_accum"), item.graup_accum);
        EXPECT_EQ(mi.qmoist_index_for_var("rel_humidity"), item.rel_hum);
        EXPECT_EQ(mi.qmoist_index_for_var("condensation_rate"), item.cond_rate);
    }
}

// Motivation: A diagnostic derived from the state must be distinguishable from
// one read out of a qmoist array, since the writer takes a different path for
// each and the absent case must write nothing at all.
TEST(Plotfile3DSelection, RelativeHumiditySourceIsExplicit)
{
    const auto satadj = make_capabilities(MoistureType::SatAdj);
    EXPECT_EQ(satadj.moisture_indices.rel_hum,
              MoistureComponentIndices::computed_from_state);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rel_humidity", satadj));

    const auto sdm = make_capabilities(MoistureType::SuperDroplets);
    EXPECT_GE(sdm.moisture_indices.rel_hum, 0);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rel_humidity", sdm));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("condensation_rate", sdm));

    const auto sam = make_capabilities(MoistureType::SAM);
    EXPECT_EQ(sam.moisture_indices.rel_hum, MoistureComponentIndices::absent);
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rel_humidity", sam));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("condensation_rate", sam));
}

// Motivation: Aggregates are summed over the component list the map hands out,
// so the list must contain exactly the carried species -- no allocated-but-idle
// slot, and no number concentration mixed into a mass sum.
TEST(Plotfile3DSelection, AggregateCompListsCoverOnlyCarriedSpecies)
{
    const auto morrison = MoistureComponentIndices::from_moisture_model(MoistureType::Morrison);
    const auto total = morrison.total_water_comps();
    ASSERT_EQ(total.size, 6);
    EXPECT_EQ(total.comp[0], RhoQ1_comp);
    EXPECT_EQ(total.comp[5], RhoQ6_comp);
    EXPECT_EQ(morrison.nonprecipitating_comps().size, 3);
    EXPECT_EQ(morrison.precipitating_comps().size, 3);

    // Liquid-only Morrison allocates the frozen species but never integrates
    // them, so they must not enter any sum.
    const auto no_ice = MoistureComponentIndices::from_moisture_model(MoistureType::Morrison_NoIce);
    const auto no_ice_total = no_ice.total_water_comps();
    ASSERT_EQ(no_ice_total.size, 3);
    EXPECT_EQ(no_ice_total.comp[0], RhoQ1_comp);
    EXPECT_EQ(no_ice_total.comp[1], RhoQ2_comp);
    EXPECT_EQ(no_ice_total.comp[2], RhoQ4_comp);
    ASSERT_EQ(no_ice.nonprecipitating_comps().size, 2);
    ASSERT_EQ(no_ice.precipitating_comps().size, 1);
    EXPECT_EQ(no_ice.precipitating_comps().comp[0], RhoQ4_comp);

    // Kessler keeps rain in the third slot
    const auto kessler = MoistureComponentIndices::from_moisture_model(MoistureType::Kessler);
    ASSERT_EQ(kessler.precipitating_comps().size, 1);
    EXPECT_EQ(kessler.precipitating_comps().comp[0], RhoQ3_comp);

    // Schemes without precipitation have no precipitating sum to write
    const auto satadj = MoistureComponentIndices::from_moisture_model(MoistureType::SatAdj);
    EXPECT_TRUE(satadj.precipitating_comps().empty());
    EXPECT_EQ(satadj.total_water_comps().size, 2);

    // A dry run has nothing to aggregate
    const MoistureComponentIndices dry;
    EXPECT_TRUE(dry.total_water_comps().empty());
    EXPECT_TRUE(dry.nonprecipitating_comps().empty());
    EXPECT_TRUE(dry.precipitating_comps().empty());
}

// Motivation: This is the behavior the moisture-index criterion buys us -- a
// component the scheme allocates but never fills is not selectable, so a
// plotfile never publishes untouched memory as data.
TEST(Plotfile3DSelection, AllocatedButUnintegratedComponentsAreNotSelectable)
{
    // Morrison_NoIce runs the Morrison class, which allocates all eleven moist
    // components.  The frozen species and their numbers are never integrated.
    const auto no_ice = make_capabilities(MoistureType::Morrison_NoIce);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qv", no_ice));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qc", no_ice));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qrain", no_ice));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qi", no_ice));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qsnow", no_ice));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qgraup", no_ice));
    // ... but the two-moment warm-rain numbers are integrated and must not be
    // dropped, which is what the index map now records.
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("nc", no_ice));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("nr", no_ice));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("ni", no_ice));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("ns", no_ice));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("ng", no_ice));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("nn", no_ice));

    // SAM_NoPrecip_NoIce likewise allocates six moist components and carries two
    const auto no_precip = make_capabilities(MoistureType::SAM_NoPrecip_NoIce);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qc", no_precip));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qrain", no_precip));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qp", no_precip));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rain_accum", no_precip));
}

// Motivation: nn and ni share a conserved component between schemes, so the two
// must never both be selectable and each must resolve to its own scheme's slot.
TEST(Plotfile3DSelection, AerosolAndIceNumberDoNotAlias)
{
    const auto wdm6 = make_capabilities(MoistureType::WDM6);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("nn", wdm6));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("ni", wdm6));
    EXPECT_EQ(wdm6.moisture_indices.comp_for_var("nn"), RhoQ8_comp);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("nc", wdm6));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("nr", wdm6));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("ns", wdm6));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("ng", wdm6));

    const auto morrison = make_capabilities(MoistureType::Morrison);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("ni", morrison));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("nn", morrison));
    EXPECT_EQ(morrison.moisture_indices.comp_for_var("ni"), RhoQ8_comp);
}

// Motivation: Raw rhoQ requests are gated by the moisture map, while the dry
// components remain bounded by the allocated state width.
TEST(Plotfile3DSelection, RawStateNamesUseTheMoistureMapAboveRhoQ1)
{
    const auto kessler = make_capabilities(MoistureType::Kessler);

    EXPECT_EQ(erf_plotfile::plot3d_conserved_component_index("density"), Rho_comp);
    EXPECT_EQ(erf_plotfile::plot3d_conserved_component_index("rhoQ1"), RhoQ1_comp);
    EXPECT_EQ(erf_plotfile::plot3d_conserved_component_index("rhoQ3"), RhoQ3_comp);
    EXPECT_EQ(erf_plotfile::plot3d_conserved_component_index("rhoQ4"), RhoQ4_comp);
    EXPECT_EQ(erf_plotfile::plot3d_conserved_component_index("rhoQ99"), -1);

    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rhoQ1", kessler));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rhoQ3", kessler));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rhoQ4", kessler));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("not_a_fixed_field", kessler));

    // The frozen slots Morrison_NoIce allocates are not writable as raw state
    const auto no_ice = make_capabilities(MoistureType::Morrison_NoIce);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rhoQ4", no_ice));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rhoQ3", no_ice));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rhoQ7", no_ice));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rhoQ8", no_ice));

    // Dry state names still fail closed on a state narrower than the name list
    auto truncated = make_capabilities(MoistureType::None, RhoScalar_comp);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rhoKE", truncated));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rhoadv_0", truncated));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rhoQ1", truncated));
}

// Motivation: The conserved state does not end at the moist species.
// SuperDroplets places each non-water species' qt/qv above the moist window, at
// RhoQ1_comp + qstate_moist_size + 2*k, and those components are allocated,
// integrated data that the moisture map deliberately does not describe.  Gating
// every "rhoQn" by the map alone silently drops them from both output paths, so
// the split between "map decides" and "allocated width decides" is pinned here.
TEST(Plotfile3DSelection, NonWaterSpeciesAboveTheMoistWindowUseTheAllocatedWidth)
{
    // SuperDroplets: qv, qc, qi, qrain, qsnow, qgraup, then two components per
    // non-water species.
    constexpr int moist_size   = 6;
    constexpr int per_species  = 2;

    const auto sdm_indices =
        MoistureComponentIndices::from_moisture_model(MoistureType::SuperDroplets);

    Plot3DSelectionCapabilities with_species;
    erf_plotfile::plot3d_set_state_capabilities(with_species, sdm_indices,
                                                moist_size, moist_size + per_species);

    // Inside the moist window the map is still the authority
    for (const char* name : {"rhoQ1", "rhoQ2", "rhoQ3", "rhoQ4", "rhoQ5", "rhoQ6"}) {
        EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available(name, with_species))
            << "moist component " << name << " was dropped";
    }

    // The one non-water species occupies the next two slots and is real data
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rhoQ7", with_species));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rhoQ8", with_species));

    // ... and the state stops there
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rhoQ9", with_species));

    // The same scheme without a non-water species stops at the moist window
    Plot3DSelectionCapabilities no_species;
    erf_plotfile::plot3d_set_state_capabilities(no_species, sdm_indices,
                                                moist_size, moist_size);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rhoQ6", no_species));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rhoQ7", no_species));

    // The allocated width must never *widen* the moist window: Morrison_NoIce
    // allocates all eleven moist components but integrates only some of them,
    // and the ones it skips stay unselectable however wide the state is.
    Plot3DSelectionCapabilities no_ice;
    erf_plotfile::plot3d_set_state_capabilities(
        no_ice, MoistureComponentIndices::from_moisture_model(MoistureType::Morrison_NoIce),
        11, 11);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rhoQ4", no_ice));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rhoQ3", no_ice));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rhoQ8", no_ice));

    // A caller that does not describe the layout gets the conservative answer:
    // the map decides about every rhoQn slot, so nothing unintegrated leaks out.
    Plot3DSelectionCapabilities undescribed;
    undescribed.moisture_indices = sdm_indices;
    undescribed.conserved_state_size = NDRY + NSCALARS + moist_size + per_species;
    EXPECT_EQ(undescribed.moist_state_size, NMOIST_max);
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("rhoQ7", undescribed));

    // The dry names are bounded by the allocated width in every case
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("density", with_species));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("rhoadv_0", with_species));
}

// Motivation: A dry run must select no moisture variable at all, whatever the
// allocated state width happens to be.
TEST(Plotfile3DSelection, DryRunSelectsNoMoistureVariables)
{
    const auto dry = make_capabilities(MoistureType::None);
    for (const auto& name : moisture_variable_names()) {
        EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available(name, dry))
            << "dry run selected moisture variable " << name;
    }
    // Dry diagnostics are untouched by the moisture map
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("theta", dry));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("mucape", dry));
}

// Motivation: The 3D plotfile and the subvolume writer must not disagree about
// which moisture variables exist; both now answer from the same index map, and
// this pins that down for every scheme.
TEST(Plotfile3DSelection, PlotfileAndSubvolumeAgreeOnMoistureNames)
{
    for (const auto type : all_moisture_types()) {
        const auto caps = make_capabilities(type);
        const auto& mi = caps.moisture_indices;
        for (const auto& name : moisture_variable_names()) {
            EXPECT_EQ(erf_plotfile::plot3d_fixed_variable_available(name, caps),
                      mi.has_derived_var(name))
                << "output paths disagree on " << name
                << " for moisture model " << amrex::getEnumNameString(type);
            EXPECT_TRUE(mi.query_var(name).governed)
                << name << " is not governed by the moisture map";
        }
    }
}

// Motivation: Names with no moisture dependence must pass through the map
// untouched, so it can be used as a filter over a mixed variable list.
TEST(Plotfile3DSelection, NonMoistureNamesAreNotGovernedByTheMap)
{
    const auto mi = MoistureComponentIndices::from_moisture_model(MoistureType::None);
    for (const char* name : {"theta", "temp", "pressure", "mucape", "qv_hse",
                             "walldist", "density", "not_a_field"}) {
        EXPECT_FALSE(mi.query_var(name).governed) << name;
        EXPECT_TRUE(mi.has_derived_var(name)) << name;
    }
}

// Motivation: NoCondensation still stores vapor and cloud water in RhoQ1 and
// RhoQ2.  A diagnosed zero in the cloud-water field is therefore distinct
// from an unavailable SHOC diagnostic, which is represented by ERF's -999
// plotfile sentinel before the native SHOC lifecycle has run.
TEST(Plotfile3DSelection, NoCondensationCloudWaterIsSelectable)
{
    const auto no_condensation = make_capabilities(MoistureType::MoistNoCondensation);
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qv", no_condensation));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qc", no_condensation));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("qn", no_condensation));
    EXPECT_TRUE(erf_plotfile::plot3d_fixed_variable_available("moist_density", no_condensation));
    EXPECT_FALSE(erf_plotfile::plot3d_fixed_variable_available("qrain", no_condensation));
}

// Motivation: The reflectivity kernels read rain, snow and graupel by hardcoded
// component, so the layout -- not just the presence of the species -- decides
// whether they can run.
TEST(Plotfile3DSelection, ReflectivitySpeciesTestChecksTheLayout)
{
    for (const auto type : {MoistureType::SAM, MoistureType::Morrison,
                            MoistureType::WSM6, MoistureType::WDM6,
                            MoistureType::SuperDroplets}) {
        EXPECT_TRUE(MoistureComponentIndices::from_moisture_model(type)
                        .has_reflectivity_species())
            << amrex::getEnumNameString(type);
    }
    for (const auto type : {MoistureType::None, MoistureType::Kessler,
                            MoistureType::SAM_NoIce, MoistureType::Morrison_NoIce,
                            MoistureType::SatAdj}) {
        EXPECT_FALSE(MoistureComponentIndices::from_moisture_model(type)
                         .has_reflectivity_species())
            << amrex::getEnumNameString(type);
    }
}

// Motivation: Optional diagnostics must be rejected before the writer can
// dereference unallocated configuration-dependent storage.
TEST(Plotfile3DSelection, OptionalStorageGroupsAreExplicit)
{
    auto caps = make_capabilities(MoistureType::None);
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
    EXPECT_TRUE(erf_plotfile::plot3d_needs_pressure({"rel_humidity"}));
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

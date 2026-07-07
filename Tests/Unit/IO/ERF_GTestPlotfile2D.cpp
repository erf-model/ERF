#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_set>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <gtest/gtest.h>

#include "ERF_Plotfile2DCatalog.H"
#include "ERF_Plotfile2DInterpolator.H"
#include "ERF_Plotfile2DFill.H"
#include "ERF_Plotfile2DMetadata.H"
#include "ERF_Plotfile2DSampledField.H"
#include "ERF_Plotfile2DSampledLevel.H"
#include "ERF_Plotfile2DWaterPath.H"
#include "ERF_Plotfile2DUtils.H"
#include "ERF_DataStruct.H"
#include "ERF_IndexDefines.H"
#include "Diagnostics/ERF_SurfaceFluxDiagnostics.H"

using namespace plotfile2d;

namespace
{

bool contains (const std::string& haystack, const std::string& needle)
{
    return haystack.find(needle) != std::string::npos;
}

amrex::Box make_test_domain ()
{
    return amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(1, 1, 2));
}

amrex::BoxArray make_test_boxarray ()
{
    amrex::BoxArray ba(make_test_domain());
    ba.maxSize(8);
    return ba;
}

amrex::BoxArray make_test_slab_boxarray ()
{
    amrex::BoxList bl = make_test_boxarray().boxList();
    for (auto& b : bl) {
        b.setRange(2, 0);
    }
    return amrex::BoxArray(std::move(bl));
}

amrex::BoxArray make_test_nodal_boxarray ()
{
    return amrex::convert(make_test_boxarray(), amrex::IntVect(1, 1, 1));
}

amrex::BoxArray make_test_xface_boxarray ()
{
    auto ba = make_test_boxarray();
    ba.surroundingNodes(0);
    return ba;
}

amrex::BoxArray make_test_yface_boxarray ()
{
    auto ba = make_test_boxarray();
    ba.surroundingNodes(1);
    return ba;
}

amrex::BoxArray make_test_zface_boxarray ()
{
    auto ba = make_test_boxarray();
    ba.surroundingNodes(2);
    return ba;
}

amrex::Geometry make_test_geometry ()
{
    amrex::Geometry geom;
    amrex::Array<int, AMREX_SPACEDIM> is_per{0, 0, 0};
    amrex::Real prob_lo[AMREX_SPACEDIM] = {0.0, 0.0, 0.0};
    amrex::Real prob_hi[AMREX_SPACEDIM] = {2.0, 2.0, 3.0};
    amrex::RealBox rb(prob_lo, prob_hi);
    geom.define(make_test_domain(), rb, 0, is_per);
    return geom;
}

void set_component_value_at_k (amrex::MultiFab& mf, int comp, int k_level, amrex::Real value)
{
    for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto& arr = mf.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (k == k_level) {
                arr(i, j, k, comp) = value;
            }
        });
    }
    amrex::Gpu::streamSynchronize();
}

void add_sampled_level_definition (const std::string& level_name,
                                   const std::string& coordinate,
                                   const std::vector<amrex::Real>& values,
                                   const std::vector<std::string>& fields,
                                   const char* units = nullptr,
                                   const char* interpolation = nullptr,
                                   const char* missing_value = nullptr)
{
    amrex::ParmParse pp("erf");
    const std::string base = "plot2d.level_set." + level_name + ".";
    pp.add(base + "coordinate", coordinate);
    if (units != nullptr) {
        pp.add(base + "units", units);
    }
    if (interpolation != nullptr) {
        pp.add(base + "interpolation", interpolation);
    }
    if (missing_value != nullptr) {
        pp.add(base + "missing_value", missing_value);
    }
    pp.addarr(base + "values", values);
    pp.addarr(base + "fields", fields);
}

} // namespace

// Motivation: The 2D output component order must follow the canonical built-in
// order, not the user input order, so plotfiles stay stable across equivalent
// input decks.
TEST(Plotfile2D, RequestedVariablesUseCanonicalOrder)
{
    const amrex::Vector<std::string> requested{"mapfac", "z_surf"};
    const amrex::Vector<std::string> available{"z_surf", "landmask", "mapfac"};

    const auto selection = select_requested_plot_variables(requested, available);

    const amrex::Vector<std::string> expected{"z_surf", "mapfac"};
    EXPECT_EQ(selection.accepted, expected);
    EXPECT_TRUE(selection.unavailable.empty());
}

// Motivation: Duplicate user requests should not expand the component list or
// create duplicate missing-variable diagnostics.
TEST(Plotfile2D, DuplicateRequestsAreCollapsed)
{
    const amrex::Vector<std::string> requested{"mapfac", "mapfac", "bogus", "bogus"};
    const amrex::Vector<std::string> available{"z_surf", "mapfac"};

    const auto selection = select_requested_plot_variables(requested, available);

    const amrex::Vector<std::string> expected_accepted{"mapfac"};
    const amrex::Vector<std::string> expected_unavailable{"bogus"};
    EXPECT_EQ(selection.accepted, expected_accepted);
    EXPECT_EQ(selection.unavailable, expected_unavailable);
}

// Motivation: Unknown plot names should be surfaced to the caller so users can
// tell whether a misspelling or unsupported diagnostic caused the omission.
TEST(Plotfile2D, UnknownRequestsAreReported)
{
    const amrex::Vector<std::string> requested{"bogus"};
    const amrex::Vector<std::string> available{"z_surf", "mapfac"};

    const auto selection = select_requested_plot_variables(requested, available);

    EXPECT_TRUE(selection.accepted.empty());
    const amrex::Vector<std::string> expected_unavailable{"bogus"};
    EXPECT_EQ(selection.unavailable, expected_unavailable);
}

// Motivation: A mixed request must keep valid variables in canonical order
// while still reporting every unavailable name exactly once.
TEST(Plotfile2D, MixedRequestsPreserveValidCanonicalOrder)
{
    const amrex::Vector<std::string> requested{"bogus", "mapfac", "z_surf", "bogus2"};
    const amrex::Vector<std::string> available{"z_surf", "landmask", "mapfac"};

    const auto selection = select_requested_plot_variables(requested, available);

    const amrex::Vector<std::string> expected_accepted{"z_surf", "mapfac"};
    const amrex::Vector<std::string> expected_unavailable{"bogus", "bogus2"};
    EXPECT_EQ(selection.accepted, expected_accepted);
    EXPECT_EQ(selection.unavailable, expected_unavailable);
}

// Motivation: An empty request should be a no-op so callers can distinguish
// "nothing requested" from "everything filtered away".
TEST(Plotfile2D, EmptyRequestsReturnEmptyLists)
{
    const amrex::Vector<std::string> requested;
    const amrex::Vector<std::string> available{"z_surf", "mapfac"};

    const auto selection = select_requested_plot_variables(requested, available);

    EXPECT_TRUE(selection.accepted.empty());
    EXPECT_TRUE(selection.unavailable.empty());
}

// Motivation: The catalog order defines the canonical 2D plotfile layout, so
// this test catches accidental reordering before it changes output metadata.
TEST(Plotfile2D, CatalogNamesMatchCanonicalOrder)
{
    const amrex::Vector<std::string> expected{
        "z_surf", "landmask", "mapfac", "lat_m", "lon_m",
        "u_star", "w_star", "t_star", "q_star", "Olen", "pblh",
        "t_surf", "q_surf", "z0", "OLR", "sens_flux", "laten_flux",
        "surf_pres", "integrated_qv", "integrated_qc", "integrated_qi",
        "integrated_qr", "integrated_qs", "integrated_qg",
        "surface_diagnostic_source",
        "sensible_heat_flux", "latent_heat_flux",
        "shoc_u_star", "shoc_Olen", "shoc_wthv_sfc"
    };

    EXPECT_EQ(plotfile2d::diagnostic_names(), expected);
}

// Motivation: The condensed water-path catalog must stay ordered after
// integrated_qv so plotfile readers and metadata remain deterministic.
TEST(Plotfile2D, WaterPathDiagnosticsFollowIntegratedQv)
{
    const auto* qv = plotfile2d::find_diagnostic("integrated_qv");
    const auto* qc = plotfile2d::find_diagnostic("integrated_qc");
    const auto* qi = plotfile2d::find_diagnostic("integrated_qi");
    const auto* qr = plotfile2d::find_diagnostic("integrated_qr");
    const auto* qs = plotfile2d::find_diagnostic("integrated_qs");
    const auto* qg = plotfile2d::find_diagnostic("integrated_qg");
    const auto* surface_source = plotfile2d::find_diagnostic("surface_diagnostic_source");

    ASSERT_NE(qv, nullptr);
    ASSERT_NE(qc, nullptr);
    ASSERT_NE(qi, nullptr);
    ASSERT_NE(qr, nullptr);
    ASSERT_NE(qs, nullptr);
    ASSERT_NE(qg, nullptr);
    ASSERT_NE(surface_source, nullptr);

    EXPECT_EQ(qv->id, plotfile2d::DiagnosticID::IntegratedQv);
    EXPECT_EQ(qc->id, plotfile2d::DiagnosticID::IntegratedQc);
    EXPECT_EQ(qi->id, plotfile2d::DiagnosticID::IntegratedQi);
    EXPECT_EQ(qr->id, plotfile2d::DiagnosticID::IntegratedQr);
    EXPECT_EQ(qs->id, plotfile2d::DiagnosticID::IntegratedQs);
    EXPECT_EQ(qg->id, plotfile2d::DiagnosticID::IntegratedQg);
    EXPECT_EQ(surface_source->id, plotfile2d::DiagnosticID::SurfaceDiagnosticSource);

    for (const auto* descriptor : {qc, qi, qr, qs, qg}) {
        EXPECT_EQ(descriptor->category, plotfile2d::DiagnosticCategory::ColumnIntegral);
        EXPECT_STREQ(descriptor->units, "kg/m^2");
        EXPECT_EQ(descriptor->missing_policy, plotfile2d::MissingPolicy::FillZeroWhenUnavailable);
    }
}

// Motivation: Each built-in 2D diagnostic name must be unique so user input
// maps to one output component and one metadata record.
TEST(Plotfile2D, CatalogNamesAreUnique)
{
    const auto names = plotfile2d::diagnostic_names();
    std::unordered_set<std::string> unique_names(names.begin(), names.end());

    EXPECT_EQ(unique_names.size(), names.size());
}

// Motivation: Each built-in 2D diagnostic ID must be unique so the catalog can
// stay stable even if a display name changes later.
TEST(Plotfile2D, CatalogIdsAreUnique)
{
    std::unordered_set<int> unique_ids;
    for (const auto& descriptor : plotfile2d::diagnostic_catalog()) {
        unique_ids.insert(static_cast<int>(descriptor.id));
    }

    EXPECT_EQ(unique_ids.size(), plotfile2d::diagnostic_catalog().size());
}

// Motivation: Catalog metadata feeds documentation and future diagnostics, so
// each built-in entry needs a name, long name, units, category, and policy.
TEST(Plotfile2D, CatalogDescriptorsHaveRequiredMetadata)
{
    for (const auto& descriptor : plotfile2d::diagnostic_catalog()) {
        EXPECT_NE(descriptor.name, nullptr);
        EXPECT_NE(descriptor.long_name, nullptr);
        EXPECT_NE(descriptor.units, nullptr);
        EXPECT_FALSE(std::string(descriptor.name).empty());
        EXPECT_FALSE(std::string(descriptor.long_name).empty());
        EXPECT_FALSE(std::string(descriptor.units).empty());

        bool valid_category = false;
        switch (descriptor.category) {
        case plotfile2d::DiagnosticCategory::Geometry:
        case plotfile2d::DiagnosticCategory::SurfaceLayer:
        case plotfile2d::DiagnosticCategory::Radiation:
        case plotfile2d::DiagnosticCategory::SurfaceFlux:
        case plotfile2d::DiagnosticCategory::PBL:
        case plotfile2d::DiagnosticCategory::SurfaceState:
        case plotfile2d::DiagnosticCategory::ColumnIntegral:
        case plotfile2d::DiagnosticCategory::SampledLevel:
            valid_category = true;
            break;
        }
        EXPECT_TRUE(valid_category);

        bool valid_missing_policy = false;
        switch (descriptor.missing_policy) {
        case plotfile2d::MissingPolicy::AlwaysAvailable:
        case plotfile2d::MissingPolicy::FillZeroWhenUnavailable:
        case plotfile2d::MissingPolicy::FillMinus999WhenUnavailable:
            valid_missing_policy = true;
            break;
        }
        EXPECT_TRUE(valid_missing_policy);
    }
}

// Motivation: The selection helper should be able to resolve a known catalog
// name back to the descriptor that defines its metadata.
TEST(Plotfile2D, FindDiagnosticReturnsDescriptorForKnownName)
{
    const auto* descriptor = plotfile2d::find_diagnostic("sens_flux");

    ASSERT_NE(descriptor, nullptr);
    EXPECT_STREQ(descriptor->name, "sens_flux");
    EXPECT_EQ(descriptor->id, plotfile2d::DiagnosticID::SensFlux);
    EXPECT_EQ(descriptor->category, plotfile2d::DiagnosticCategory::SurfaceFlux);
    EXPECT_EQ(descriptor->missing_policy, plotfile2d::MissingPolicy::FillMinus999WhenUnavailable);
}

// Motivation: pblh now follows the active PBL diagnostic provider, so the
// catalog metadata must describe the provider rather than only SurfaceLayer.
TEST(Plotfile2D, FindDiagnosticReturnsDescriptorForPblh)
{
    const auto* descriptor = plotfile2d::find_diagnostic("pblh");

    ASSERT_NE(descriptor, nullptr);
    EXPECT_STREQ(descriptor->name, "pblh");
    EXPECT_EQ(descriptor->id, plotfile2d::DiagnosticID::Pblh);
    EXPECT_EQ(descriptor->category, plotfile2d::DiagnosticCategory::SurfaceLayer);
    EXPECT_STREQ(descriptor->units, "m");
    EXPECT_TRUE(contains(descriptor->long_name, "active PBL diagnostic provider"));
    EXPECT_EQ(descriptor->missing_policy, plotfile2d::MissingPolicy::FillMinus999WhenUnavailable);
}

// Motivation: The new provenance mask is a public 2D output, so its catalog
// entry must be findable and carry the metadata that documents its output
// convention.
TEST(Plotfile2D, FindDiagnosticReturnsDescriptorForSurfaceDiagnosticSource)
{
    const auto* descriptor = plotfile2d::find_diagnostic("surface_diagnostic_source");

    ASSERT_NE(descriptor, nullptr);
    EXPECT_STREQ(descriptor->name, "surface_diagnostic_source");
    EXPECT_EQ(descriptor->id, plotfile2d::DiagnosticID::SurfaceDiagnosticSource);
    EXPECT_FALSE(std::string(descriptor->long_name).empty());
    EXPECT_FALSE(std::string(descriptor->units).empty());
    EXPECT_EQ(descriptor->missing_policy, plotfile2d::MissingPolicy::FillMinus999WhenUnavailable);
}

// Motivation: Pressure sampled-level definitions must parse with the provided
// units and preserve the requested value list and field order.
TEST(Plotfile2DSampledLevel, ParsesPressureLevelSet)
{
    add_sampled_level_definition("pressure_parse_pressure",
                                 "pressure",
                                 {850.0, 700.0},
                                 {"theta", "temp", "qv", "qc"},
                                 "hPa",
                                 "linear");

    const auto level_set = parse_sampled_level_definition("pressure_parse_pressure", "erf");

    EXPECT_EQ(level_set.coordinate, SampledCoordinate::Pressure);
    EXPECT_EQ(level_set.units, "hPa");
    ASSERT_EQ(level_set.values.size(), 2u);
    EXPECT_DOUBLE_EQ(level_set.values[0], 850.0);
    EXPECT_DOUBLE_EQ(level_set.values[1], 700.0);
    ASSERT_EQ(level_set.fields.size(), 4u);
    EXPECT_EQ(level_set.fields[0], "theta");
    EXPECT_EQ(level_set.fields[1], "temp");
    EXPECT_EQ(level_set.fields[2], "qv");
    EXPECT_EQ(level_set.fields[3], "qc");
    EXPECT_EQ(level_set.interpolation, SampledInterpolation::Linear);
    EXPECT_EQ(validate_sampled_level_definition(level_set), "");
}

// Motivation: Height-agl sampled-level definitions should default to meters
// and linear interpolation when those keys are omitted.
TEST(Plotfile2DSampledLevel, ParsesHeightAGLLevelSet)
{
    add_sampled_level_definition("height_agl_parse_height",
                                 "height_agl",
                                 {100.0, 500.0},
                                 {"theta", "qv", "qc"});

    const auto level_set = parse_sampled_level_definition("height_agl_parse_height", "erf");

    EXPECT_EQ(level_set.coordinate, SampledCoordinate::HeightAGL);
    EXPECT_EQ(level_set.units, "m");
    EXPECT_EQ(level_set.interpolation, SampledInterpolation::Linear);
    ASSERT_EQ(level_set.values.size(), 2u);
    EXPECT_DOUBLE_EQ(level_set.values[0], 100.0);
    EXPECT_DOUBLE_EQ(level_set.values[1], 500.0);
}

// Motivation: Unknown coordinates should fail fast with a clear parser error
// instead of falling through to a later interpolation failure.
TEST(Plotfile2DSampledLevel, RejectsUnknownCoordinate)
{
    // AMReX initialization can leave helper threads alive, so use the
    // threadsafe death-test style instead of fork-based death tests.
    GTEST_FLAG_SET(death_test_style, "threadsafe");

    add_sampled_level_definition("bad_coordinate",
                                 "bogus_coordinate",
                                 {10.0},
                                 {"theta"});

    EXPECT_DEATH(parse_sampled_level_definition("bad_coordinate", "erf"),
                 "Unknown sampled-level coordinate");
}

// Motivation: Unsupported units should be rejected at parse time so sampled
// metadata cannot advertise a unit that the sampler does not canonicalize.
TEST(Plotfile2DSampledLevel, RejectsUnsupportedUnits)
{
    // AMReX initialization can leave helper threads alive, so use the
    // threadsafe death-test style instead of fork-based death tests.
    GTEST_FLAG_SET(death_test_style, "threadsafe");

    add_sampled_level_definition("bad_units",
                                 "pressure",
                                 {850.0},
                                 {"theta"},
                                 "kPa",
                                 "linear");

    EXPECT_DEATH(parse_sampled_level_definition("bad_units", "erf"),
                 "unsupported units");
}

// Motivation: Isentropic output needs a crossing policy, so it must stay
// rejected until that policy exists.
TEST(Plotfile2DSampledLevel, RejectsIsentropicUntilCrossingPolicyExists)
{
    // AMReX initialization can leave helper threads alive, so use the
    // threadsafe death-test style instead of fork-based death tests.
    GTEST_FLAG_SET(death_test_style, "threadsafe");

    add_sampled_level_definition("bad_isentropic",
                                 "isentropic",
                                 {300.0},
                                 {"theta"});

    EXPECT_DEATH(parse_sampled_level_definition("bad_isentropic", "erf"),
                 "crossing policy");
}

// Motivation: Derived moisture aliases are out of scope for this sampled
// output path and must not be accepted as canonical field names.
TEST(Plotfile2DSampledLevel, RejectsDerivedMoistureAliases)
{
    const amrex::Vector<std::string> requested{"qrain", "qsnow", "qgraup"};
    SolverChoice sc;
    const auto selection = select_requested_sampled_fields(requested, sc);

    EXPECT_TRUE(selection.accepted.empty());
    ASSERT_EQ(selection.unavailable.size(), 3u);
    EXPECT_EQ(selection.unavailable[0], "qrain");
    EXPECT_EQ(selection.unavailable[1], "qsnow");
    EXPECT_EQ(selection.unavailable[2], "qgraup");
}

// Motivation: Dry runs should expose the thermodynamic and geometric sampled
// fields but not moisture species.
TEST(Plotfile2DSampledField, DryRunExcludesMoistureSpecies)
{
    SolverChoice sc;
    sc.moisture_type = MoistureType::None;

    const auto available = available_sampled_field_names(sc);

    for (const char* name : {"rho", "theta", "temp", "pressure", "height_msl", "height_agl",
                             "u_east", "v_north", "w", "wind_speed", "wind_dir"}) {
        EXPECT_NE(std::find(available.begin(), available.end(), name), available.end());
    }
    for (const char* name : {"qv", "qc", "qi", "qr", "qs", "qg"}) {
        EXPECT_EQ(std::find(available.begin(), available.end(), name), available.end());
    }
}

// Motivation: Dry-run thermodynamic sampled fields still use qv = 0 in the
// EOS helpers even when the sampled qv field itself is unavailable.
TEST(Plotfile2DSampledField, DryRunUsesZeroVaporForThermodynamicFields)
{
    const auto ba = make_test_boxarray();
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab z_phys_cc(ba, dm, 1, 0);
    amrex::MultiFab z_phys_nd(make_test_nodal_boxarray(), dm, 1, 0);
    amrex::MultiFab dst(ba, dm, 2, 0);

    cons.setVal(amrex::Real(0.0));
    z_phys_cc.setVal(amrex::Real(0.0));
    z_phys_nd.setVal(amrex::Real(0.0));
    dst.setVal(amrex::Real(-7.0));

    set_component_value_at_k(cons, Rho_comp, 0, amrex::Real(2.0));
    set_component_value_at_k(cons, RhoTheta_comp, 0, amrex::Real(500.0));

    SolverChoice sc;
    sc.moisture_type = MoistureType::None;
    const auto& mi = sc.moisture_indices;
    const auto available = available_sampled_field_names(sc);
    EXPECT_EQ(std::find(available.begin(), available.end(), "qv"), available.end());

    Plotfile2DOutputDescriptor pressure_descriptor;
    pressure_descriptor.name = "pressure_k_0";
    pressure_descriptor.long_name = "Pressure sampled on model index levels";
    pressure_descriptor.units = "Pa";
    pressure_descriptor.category = DiagnosticCategory::SampledLevel;
    pressure_descriptor.missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
    pressure_descriptor.missing_value = amrex::Real(-999.0);
    pressure_descriptor.sampled_level = SampledLevelMetadata{
        "dry_run",
        "pressure",
        SampledVerticalCoordinateMetadata{
            "model_index",
            amrex::Real(0.0),
            "1",
            amrex::Real(0.0),
            "1",
            "none"
        }
    };

    Plotfile2DOutputDescriptor temp_descriptor;
    temp_descriptor.name = "temp_k_0";
    temp_descriptor.long_name = "Temperature sampled on model index levels";
    temp_descriptor.units = "K";
    temp_descriptor.category = DiagnosticCategory::SampledLevel;
    temp_descriptor.missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
    temp_descriptor.missing_value = amrex::Real(-999.0);
    temp_descriptor.sampled_level = SampledLevelMetadata{
        "dry_run",
        "temp",
        SampledVerticalCoordinateMetadata{
            "model_index",
            amrex::Real(0.0),
            "1",
            amrex::Real(0.0),
            "1",
            "none"
        }
    };

    plotfile2d::fill_sampled_level_component(dst, 0, pressure_descriptor,
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             mi, 0, 2);
    plotfile2d::fill_sampled_level_component(dst, 1, temp_descriptor,
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             mi, 0, 2);
    amrex::Gpu::streamSynchronize();

    EXPECT_DOUBLE_EQ(dst.min(0), getPgivenRTh(amrex::Real(500.0), amrex::Real(0.0)));
    EXPECT_DOUBLE_EQ(dst.max(0), getPgivenRTh(amrex::Real(500.0), amrex::Real(0.0)));
    EXPECT_DOUBLE_EQ(dst.min(1), getTgivenRandRTh(amrex::Real(2.0), amrex::Real(500.0), amrex::Real(0.0)));
    EXPECT_DOUBLE_EQ(dst.max(1), getTgivenRandRTh(amrex::Real(2.0), amrex::Real(500.0), amrex::Real(0.0)));
}

// Motivation: Wind fields should be first-class sampled-level fields, not
// special-case output names that bypass descriptor metadata.
TEST(Plotfile2DSampledField, WindFieldsAreCataloguedAndAlwaysAvailable)
{
    SolverChoice sc;
    const auto available = available_sampled_field_names(sc);

    for (const char* name : {"u_east", "v_north", "w", "wind_speed", "wind_dir"}) {
        EXPECT_NE(std::find(available.begin(), available.end(), name), available.end());
        const auto* descriptor = find_sampled_field(name);
        ASSERT_NE(descriptor, nullptr);
        EXPECT_STRNE(descriptor->long_name, "");
        if (std::string(name) == "wind_dir") {
            EXPECT_STREQ(descriptor->units, "degrees");
        } else {
            EXPECT_STREQ(descriptor->units, "m/s");
        }
    }
}

// Motivation: Scheme-aware availability should follow the active moisture
// indices rather than a hard-coded microphysics switch.
TEST(Plotfile2DSampledField, KesslerAvailabilityMatchesMoistureIndices)
{
    SolverChoice sc;
    sc.moisture_type = MoistureType::Kessler;
    sc.moisture_indices = MoistureComponentIndices(RhoQ1_comp, RhoQ2_comp, -1, RhoQ3_comp);

    const auto available = available_sampled_field_names(sc);

    for (const char* name : {"qv", "qc", "qr"}) {
        EXPECT_NE(std::find(available.begin(), available.end(), name), available.end());
    }
    for (const char* name : {"qi", "qs", "qg"}) {
        EXPECT_EQ(std::find(available.begin(), available.end(), name), available.end());
    }
}

// Motivation: Full bulk schemes should expose all sampled microphysical mass
// species.
TEST(Plotfile2DSampledField, FullIceAvailabilityIncludesAllMassSpecies)
{
    SolverChoice sc;
    sc.moisture_type = MoistureType::Morrison;
    sc.moisture_indices = MoistureComponentIndices(RhoQ1_comp, RhoQ2_comp,
                                                   RhoQ3_comp, RhoQ4_comp,
                                                   RhoQ5_comp, RhoQ6_comp);

    const auto available = available_sampled_field_names(sc);

    for (const char* name : {"qv", "qc", "qi", "qr", "qs", "qg"}) {
        EXPECT_NE(std::find(available.begin(), available.end(), name), available.end());
    }
}

// Motivation: The bracket rule must work for both increasing and decreasing
// coordinates so sampled pressure and height behave the same way.
TEST(Plotfile2DInterpolator, BracketRuleHandlesIncreasingAndDecreasingCoordinates)
{
    EXPECT_TRUE(sampled_target_is_bracketed(amrex::Real(5.0), amrex::Real(0.0), amrex::Real(10.0)));
    EXPECT_TRUE(sampled_target_is_bracketed(amrex::Real(5.0), amrex::Real(10.0), amrex::Real(0.0)));
    EXPECT_FALSE(sampled_target_is_bracketed(amrex::Real(15.0), amrex::Real(0.0), amrex::Real(10.0)));
}

// Motivation: Linear interpolation is the core numerical step for sampled
// height and pressure outputs.
TEST(Plotfile2DInterpolator, LinearInterpolateMatchesSyntheticValues)
{
    EXPECT_DOUBLE_EQ(linear_interpolate(amrex::Real(10.0), amrex::Real(20.0),
                                        amrex::Real(0.0), amrex::Real(10.0),
                                        amrex::Real(5.0)),
                     amrex::Real(15.0));
}

// Motivation: Model-index sampling should copy the requested level exactly.
TEST(Plotfile2DInterpolator, ModelIndexCopiesExactLevel)
{
    const auto ba = make_test_boxarray();
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab z_phys_cc(ba, dm, 1, 0);
    amrex::MultiFab z_phys_nd(make_test_nodal_boxarray(), dm, 1, 0);
    amrex::MultiFab dst(ba, dm, 1, 0);

    cons.setVal(amrex::Real(0.0));
    z_phys_cc.setVal(amrex::Real(0.0));
    z_phys_nd.setVal(amrex::Real(0.0));
    dst.setVal(amrex::Real(-7.0));

    set_component_value_at_k(cons, Rho_comp, 1, amrex::Real(2.0));
    set_component_value_at_k(cons, RhoTheta_comp, 1, amrex::Real(520.0));

    Plotfile2DOutputDescriptor descriptor;
    descriptor.name = "theta_k_1";
    descriptor.long_name = "Potential temperature sampled on model index levels";
    descriptor.units = "K";
    descriptor.category = DiagnosticCategory::SampledLevel;
    descriptor.missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
    descriptor.missing_value = amrex::Real(-999.0);
    descriptor.sampled_level = SampledLevelMetadata{
        "native_k",
        "theta",
        SampledVerticalCoordinateMetadata{
            "model_index",
            amrex::Real(1.0),
            "1",
            amrex::Real(1.0),
            "1",
            "none"
        }
    };

    plotfile2d::fill_sampled_level_component(dst, 0, descriptor,
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 2);
    amrex::Gpu::streamSynchronize();

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(520.0) / amrex::Real(2.0));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(520.0) / amrex::Real(2.0));
}

// Motivation: Height-msl sampling should interpolate a simple synthetic
// column without changing the sampled-value convention.
TEST(Plotfile2DInterpolator, HeightMSLInterpolatesSyntheticColumn)
{
    const auto ba = make_test_boxarray();
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab z_phys_cc(ba, dm, 1, 0);
    amrex::MultiFab z_phys_nd(make_test_nodal_boxarray(), dm, 1, 0);
    amrex::MultiFab dst(ba, dm, 1, 0);

    cons.setVal(amrex::Real(0.0));
    z_phys_cc.setVal(amrex::Real(0.0));
    z_phys_nd.setVal(amrex::Real(0.0));
    dst.setVal(amrex::Real(-7.0));

    set_component_value_at_k(cons, Rho_comp, 0, amrex::Real(10.0));
    set_component_value_at_k(cons, Rho_comp, 1, amrex::Real(20.0));
    set_component_value_at_k(cons, Rho_comp, 2, amrex::Real(30.0));
    set_component_value_at_k(z_phys_cc, 0, 0, amrex::Real(0.0));
    set_component_value_at_k(z_phys_cc, 0, 1, amrex::Real(1000.0));
    set_component_value_at_k(z_phys_cc, 0, 2, amrex::Real(2000.0));

    Plotfile2DOutputDescriptor descriptor;
    descriptor.name = "rho_z_msl_500m";
    descriptor.long_name = "Density sampled on height-msl levels";
    descriptor.units = "kg/m^3";
    descriptor.category = DiagnosticCategory::SampledLevel;
    descriptor.missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
    descriptor.missing_value = amrex::Real(-999.0);
    descriptor.sampled_level = SampledLevelMetadata{
        "upper_air",
        "rho",
        SampledVerticalCoordinateMetadata{
            "height_msl",
            amrex::Real(500.0),
            "m",
            amrex::Real(500.0),
            "m",
            "linear"
        }
    };

    plotfile2d::fill_sampled_level_component(dst, 0, descriptor,
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 2);
    amrex::Gpu::streamSynchronize();

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(15.0));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(15.0));
}

// Motivation: Height-agl sampling should use the local terrain height rather
// than a global altitude origin.
TEST(Plotfile2DInterpolator, HeightAGLUsesLocalSurfaceHeight)
{
    const auto ba = make_test_boxarray();
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab z_phys_cc(ba, dm, 1, 0);
    amrex::MultiFab z_phys_nd(make_test_nodal_boxarray(), dm, 1, 0);
    amrex::MultiFab dst(ba, dm, 1, 0);

    cons.setVal(amrex::Real(0.0));
    z_phys_cc.setVal(amrex::Real(0.0));
    z_phys_nd.setVal(amrex::Real(0.0));
    dst.setVal(amrex::Real(-7.0));

    set_component_value_at_k(cons, Rho_comp, 0, amrex::Real(10.0));
    set_component_value_at_k(cons, Rho_comp, 1, amrex::Real(20.0));
    set_component_value_at_k(cons, Rho_comp, 2, amrex::Real(30.0));
    set_component_value_at_k(z_phys_nd, 0, 0, amrex::Real(100.0));
    set_component_value_at_k(z_phys_nd, 0, 1, amrex::Real(200.0));
    set_component_value_at_k(z_phys_nd, 0, 2, amrex::Real(300.0));
    set_component_value_at_k(z_phys_nd, 0, 3, amrex::Real(400.0));

    Plotfile2DOutputDescriptor descriptor;
    descriptor.name = "rho_z_agl_50m";
    descriptor.long_name = "Density sampled on height-agl levels";
    descriptor.units = "kg/m^3";
    descriptor.category = DiagnosticCategory::SampledLevel;
    descriptor.missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
    descriptor.missing_value = amrex::Real(-999.0);
    descriptor.sampled_level = SampledLevelMetadata{
        "bl",
        "rho",
        SampledVerticalCoordinateMetadata{
            "height_agl",
            amrex::Real(50.0),
            "m",
            amrex::Real(50.0),
            "m",
            "linear"
        }
    };

    plotfile2d::fill_sampled_level_component(dst, 0, descriptor,
                                             cons, nullptr, z_phys_nd, false,
                                             MoistureComponentIndices(), 0, 2);
    amrex::Gpu::streamSynchronize();

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(10.0));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(10.0));
}

// Motivation: Pressure sampling should interpolate a monotonic column using
// the EOS-derived pressure coordinate.
TEST(Plotfile2DInterpolator, PressureInterpolatesMonotonicColumn)
{
    const auto ba = make_test_boxarray();
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab z_phys_cc(ba, dm, 1, 0);
    amrex::MultiFab z_phys_nd(make_test_nodal_boxarray(), dm, 1, 0);
    amrex::MultiFab dst(ba, dm, 1, 0);

    cons.setVal(amrex::Real(0.0));
    z_phys_cc.setVal(amrex::Real(0.0));
    z_phys_nd.setVal(amrex::Real(0.0));
    dst.setVal(amrex::Real(-7.0));

    const amrex::Real qv = amrex::Real(0.0);
    set_component_value_at_k(cons, Rho_comp, 0, amrex::Real(300.0));
    set_component_value_at_k(cons, Rho_comp, 1, amrex::Real(310.0));
    set_component_value_at_k(cons, Rho_comp, 2, amrex::Real(320.0));
    set_component_value_at_k(cons, RhoTheta_comp, 0, getRhoThetagivenP(amrex::Real(90000.0), qv));
    set_component_value_at_k(cons, RhoTheta_comp, 1, getRhoThetagivenP(amrex::Real(80000.0), qv));
    set_component_value_at_k(cons, RhoTheta_comp, 2, getRhoThetagivenP(amrex::Real(70000.0), qv));

    Plotfile2DOutputDescriptor descriptor;
    descriptor.name = "rho_p_85000Pa";
    descriptor.long_name = "Density sampled on pressure levels";
    descriptor.units = "kg/m^3";
    descriptor.category = DiagnosticCategory::SampledLevel;
    descriptor.missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
    descriptor.missing_value = amrex::Real(-999.0);
    descriptor.sampled_level = SampledLevelMetadata{
        "upper_air",
        "rho",
        SampledVerticalCoordinateMetadata{
            "pressure",
            amrex::Real(85000.0),
            "Pa",
            amrex::Real(85000.0),
            "Pa",
            "linear"
        }
    };

    plotfile2d::fill_sampled_level_component(dst, 0, descriptor,
                                             cons, nullptr, z_phys_nd, false,
                                             MoistureComponentIndices(), 0, 2);
    amrex::Gpu::streamSynchronize();

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(305.0));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(305.0));
}

// Motivation: Targets outside the sampled column must write the configured
// missing value rather than extrapolating.
TEST(Plotfile2DInterpolator, OutOfRangeTargetFillsMissingValue)
{
    const auto ba = make_test_boxarray();
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab z_phys_cc(ba, dm, 1, 0);
    amrex::MultiFab z_phys_nd(make_test_nodal_boxarray(), dm, 1, 0);
    amrex::MultiFab dst(ba, dm, 1, 0);

    cons.setVal(amrex::Real(0.0));
    z_phys_cc.setVal(amrex::Real(0.0));
    z_phys_nd.setVal(amrex::Real(0.0));
    dst.setVal(amrex::Real(-7.0));

    set_component_value_at_k(cons, Rho_comp, 0, amrex::Real(300.0));
    set_component_value_at_k(cons, Rho_comp, 1, amrex::Real(310.0));
    set_component_value_at_k(cons, Rho_comp, 2, amrex::Real(320.0));
    set_component_value_at_k(cons, RhoTheta_comp, 0, getRhoThetagivenP(amrex::Real(90000.0), amrex::Real(0.0)));
    set_component_value_at_k(cons, RhoTheta_comp, 1, getRhoThetagivenP(amrex::Real(80000.0), amrex::Real(0.0)));
    set_component_value_at_k(cons, RhoTheta_comp, 2, getRhoThetagivenP(amrex::Real(70000.0), amrex::Real(0.0)));

    Plotfile2DOutputDescriptor descriptor;
    descriptor.name = "rho_p_65000Pa";
    descriptor.long_name = "Density sampled on pressure levels";
    descriptor.units = "kg/m^3";
    descriptor.category = DiagnosticCategory::SampledLevel;
    descriptor.missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
    descriptor.missing_value = amrex::Real(-999.0);
    descriptor.sampled_level = SampledLevelMetadata{
        "upper_air",
        "rho",
        SampledVerticalCoordinateMetadata{
            "pressure",
            amrex::Real(65000.0),
            "Pa",
            amrex::Real(65000.0),
            "Pa",
            "linear"
        }
    };

    plotfile2d::fill_sampled_level_component(dst, 0, descriptor,
                                             cons, nullptr, z_phys_nd, false,
                                             MoistureComponentIndices(), 0, 2);
    amrex::Gpu::streamSynchronize();

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(-999.0));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(-999.0));
}

// Motivation: Microphysical sampled fields are reported as mixing ratios, so
// the conserved rho*q state must be divided by rho before output.
TEST(Plotfile2DInterpolator, MicrophysicsFieldsUseMixingRatio)
{
    const auto ba = make_test_boxarray();
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab z_phys_cc(ba, dm, 1, 0);
    amrex::MultiFab z_phys_nd(make_test_nodal_boxarray(), dm, 1, 0);
    amrex::MultiFab dst(ba, dm, 1, 0);

    cons.setVal(amrex::Real(0.0));
    z_phys_cc.setVal(amrex::Real(0.0));
    z_phys_nd.setVal(amrex::Real(0.0));
    dst.setVal(amrex::Real(-7.0));

    set_component_value_at_k(cons, Rho_comp, 0, amrex::Real(2.0));
    set_component_value_at_k(cons, RhoQ2_comp, 0, amrex::Real(0.004));

    Plotfile2DOutputDescriptor descriptor;
    descriptor.name = "qc_k_0";
    descriptor.long_name = "Cloud water sampled on model index levels";
    descriptor.units = "kg/kg";
    descriptor.category = DiagnosticCategory::SampledLevel;
    descriptor.missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
    descriptor.missing_value = amrex::Real(-999.0);
    descriptor.sampled_level = SampledLevelMetadata{
        "microphysics",
        "qc",
        SampledVerticalCoordinateMetadata{
            "model_index",
            amrex::Real(0.0),
            "1",
            amrex::Real(0.0),
            "1",
            "none"
        }
    };

    plotfile2d::fill_sampled_level_component(dst, 0, descriptor,
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(RhoQ1_comp, RhoQ2_comp, -1, -1),
                                             0, 2);
    amrex::Gpu::streamSynchronize();

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(0.002));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(0.002));
}

// Motivation: Sampled winds should be cell-centered outputs even though their
// source velocities live on staggered faces.
TEST(Plotfile2DInterpolator, SampledWindsDestaggerWithoutRotation)
{
    const auto cell_ba = make_test_boxarray();
    const auto xba = make_test_xface_boxarray();
    const auto yba = make_test_yface_boxarray();
    const auto zba = make_test_zface_boxarray();
    amrex::DistributionMapping dm(cell_ba);
    amrex::MultiFab cons(cell_ba, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab z_phys_cc(cell_ba, dm, 1, 0);
    amrex::MultiFab z_phys_nd(make_test_nodal_boxarray(), dm, 1, 0);
    amrex::MultiFab xvel(xba, dm, 1, 0);
    amrex::MultiFab yvel(yba, dm, 1, 0);
    amrex::MultiFab zvel(zba, dm, 1, 0);
    amrex::MultiFab dst(cell_ba, dm, 5, 0);

    cons.setVal(amrex::Real(0.0));
    z_phys_cc.setVal(amrex::Real(0.0));
    z_phys_nd.setVal(amrex::Real(0.0));
    xvel.setVal(amrex::Real(1.0));
    yvel.setVal(amrex::Real(2.0));
    zvel.setVal(amrex::Real(3.0));
    dst.setVal(amrex::Real(-7.0));

    const auto wind_sources = plotfile2d::SampledWindSources{&xvel, &yvel, &zvel, nullptr, nullptr};

    auto make_descriptor = [](const char* field_name, const char* units) {
        Plotfile2DOutputDescriptor descriptor;
        descriptor.name = field_name;
        descriptor.long_name = field_name;
        descriptor.units = units;
        descriptor.category = DiagnosticCategory::SampledLevel;
        descriptor.missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
        descriptor.missing_value = amrex::Real(-999.0);
        descriptor.sampled_level = SampledLevelMetadata{
            "upper_air_winds",
            field_name,
            SampledVerticalCoordinateMetadata{
                "model_index",
                amrex::Real(0.0),
                "1",
                amrex::Real(0.0),
                "1",
                "none"
            }
        };
        return descriptor;
    };

    plotfile2d::fill_sampled_level_component(dst, 0, make_descriptor("u_east", "m/s"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 2, wind_sources);
    plotfile2d::fill_sampled_level_component(dst, 1, make_descriptor("v_north", "m/s"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 2, wind_sources);
    plotfile2d::fill_sampled_level_component(dst, 2, make_descriptor("w", "m/s"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 2, wind_sources);
    plotfile2d::fill_sampled_level_component(dst, 3, make_descriptor("wind_speed", "m/s"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 2, wind_sources);
    plotfile2d::fill_sampled_level_component(dst, 4, make_descriptor("wind_dir", "degrees"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 2, wind_sources);
    amrex::Gpu::streamSynchronize();

    const amrex::Real pi = amrex::Real(3.141592653589793238462643383279502884L);
    const amrex::Real expected_dir =
        amrex::Real(270.0) - std::atan2(amrex::Real(2.0), amrex::Real(1.0)) * amrex::Real(180.0) / pi;
    const amrex::Real expected_speed = std::sqrt(amrex::Real(5.0));

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(1.0));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(1.0));
    EXPECT_DOUBLE_EQ(dst.min(1), amrex::Real(2.0));
    EXPECT_DOUBLE_EQ(dst.max(1), amrex::Real(2.0));
    EXPECT_DOUBLE_EQ(dst.min(2), amrex::Real(3.0));
    EXPECT_DOUBLE_EQ(dst.max(2), amrex::Real(3.0));
    EXPECT_DOUBLE_EQ(dst.min(3), expected_speed);
    EXPECT_DOUBLE_EQ(dst.max(3), expected_speed);
    EXPECT_NEAR(dst.min(4), expected_dir, 1.0e-12);
    EXPECT_NEAR(dst.max(4), expected_dir, 1.0e-12);
}

// Motivation: Plotting winds must use earth-relative components when rotation
// coefficients are supplied.
TEST(Plotfile2DInterpolator, SampledWindsApplyRotation)
{
    const auto cell_ba = make_test_boxarray();
    const auto xba = make_test_xface_boxarray();
    const auto yba = make_test_yface_boxarray();
    const auto zba = make_test_zface_boxarray();
    amrex::DistributionMapping dm(cell_ba);
    amrex::MultiFab cons(cell_ba, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab z_phys_cc(cell_ba, dm, 1, 0);
    amrex::MultiFab z_phys_nd(make_test_nodal_boxarray(), dm, 1, 0);
    amrex::MultiFab xvel(xba, dm, 1, 0);
    amrex::MultiFab yvel(yba, dm, 1, 0);
    amrex::MultiFab zvel(zba, dm, 1, 0);
    amrex::MultiFab cos_alpha(cell_ba, dm, 1, 0);
    amrex::MultiFab sin_alpha(cell_ba, dm, 1, 0);
    amrex::MultiFab dst(cell_ba, dm, 4, 0);

    cons.setVal(amrex::Real(0.0));
    z_phys_cc.setVal(amrex::Real(0.0));
    z_phys_nd.setVal(amrex::Real(0.0));
    xvel.setVal(amrex::Real(1.0));
    yvel.setVal(amrex::Real(2.0));
    zvel.setVal(amrex::Real(3.0));
    cos_alpha.setVal(amrex::Real(0.0));
    sin_alpha.setVal(amrex::Real(1.0));
    dst.setVal(amrex::Real(-7.0));

    const auto wind_sources = plotfile2d::SampledWindSources{&xvel, &yvel, &zvel, &cos_alpha, &sin_alpha};

    auto make_descriptor = [](const char* field_name, const char* units) {
        Plotfile2DOutputDescriptor descriptor;
        descriptor.name = field_name;
        descriptor.long_name = field_name;
        descriptor.units = units;
        descriptor.category = DiagnosticCategory::SampledLevel;
        descriptor.missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
        descriptor.missing_value = amrex::Real(-999.0);
        descriptor.sampled_level = SampledLevelMetadata{
            "upper_air_winds",
            field_name,
            SampledVerticalCoordinateMetadata{
                "model_index",
                amrex::Real(0.0),
                "1",
                amrex::Real(0.0),
                "1",
                "none"
            }
        };
        return descriptor;
    };

    plotfile2d::fill_sampled_level_component(dst, 0, make_descriptor("u_east", "m/s"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 2, wind_sources);
    plotfile2d::fill_sampled_level_component(dst, 1, make_descriptor("v_north", "m/s"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 2, wind_sources);
    plotfile2d::fill_sampled_level_component(dst, 2, make_descriptor("wind_speed", "m/s"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 2, wind_sources);
    plotfile2d::fill_sampled_level_component(dst, 3, make_descriptor("wind_dir", "degrees"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 2, wind_sources);
    amrex::Gpu::streamSynchronize();

    const amrex::Real pi = amrex::Real(3.141592653589793238462643383279502884L);
    const amrex::Real expected_dir =
        amrex::Real(270.0) - std::atan2(amrex::Real(1.0), amrex::Real(-2.0)) * amrex::Real(180.0) / pi;

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(-2.0));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(-2.0));
    EXPECT_DOUBLE_EQ(dst.min(1), amrex::Real(1.0));
    EXPECT_DOUBLE_EQ(dst.max(1), amrex::Real(1.0));
    EXPECT_DOUBLE_EQ(dst.min(2), std::sqrt(amrex::Real(5.0)));
    EXPECT_DOUBLE_EQ(dst.max(2), std::sqrt(amrex::Real(5.0)));
    EXPECT_NEAR(dst.min(3), expected_dir, 1.0e-12);
    EXPECT_NEAR(dst.max(3), expected_dir, 1.0e-12);
}

// Motivation: Wind direction is undefined for calm winds, so the output should
// use the configured missing value instead of an arbitrary angle.
TEST(Plotfile2DInterpolator, WindDirReturnsMissingForCalmWinds)
{
    const auto cell_ba = make_test_boxarray();
    const auto xba = make_test_xface_boxarray();
    const auto yba = make_test_yface_boxarray();
    const auto zba = make_test_zface_boxarray();
    amrex::DistributionMapping dm(cell_ba);
    amrex::MultiFab cons(cell_ba, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab z_phys_cc(cell_ba, dm, 1, 0);
    amrex::MultiFab z_phys_nd(make_test_nodal_boxarray(), dm, 1, 0);
    amrex::MultiFab xvel(xba, dm, 1, 0);
    amrex::MultiFab yvel(yba, dm, 1, 0);
    amrex::MultiFab zvel(zba, dm, 1, 0);
    amrex::MultiFab dst(cell_ba, dm, 2, 0);

    cons.setVal(amrex::Real(0.0));
    z_phys_cc.setVal(amrex::Real(0.0));
    z_phys_nd.setVal(amrex::Real(0.0));
    xvel.setVal(amrex::Real(0.0));
    yvel.setVal(amrex::Real(0.0));
    zvel.setVal(amrex::Real(0.0));
    dst.setVal(amrex::Real(-7.0));

    const auto wind_sources = plotfile2d::SampledWindSources{&xvel, &yvel, &zvel, nullptr, nullptr};

    Plotfile2DOutputDescriptor speed_descriptor;
    speed_descriptor.name = "wind_speed_k_0";
    speed_descriptor.long_name = "Horizontal wind speed sampled on model index levels";
    speed_descriptor.units = "m/s";
    speed_descriptor.category = DiagnosticCategory::SampledLevel;
    speed_descriptor.missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
    speed_descriptor.missing_value = amrex::Real(-999.0);
    speed_descriptor.sampled_level = SampledLevelMetadata{
        "upper_air_winds",
        "wind_speed",
        SampledVerticalCoordinateMetadata{
            "model_index",
            amrex::Real(0.0),
            "1",
            amrex::Real(0.0),
            "1",
            "none"
        }
    };

    Plotfile2DOutputDescriptor dir_descriptor = speed_descriptor;
    dir_descriptor.name = "wind_dir_k_0";
    dir_descriptor.long_name = "Meteorological wind direction sampled on model index levels";
    dir_descriptor.units = "degrees";
    dir_descriptor.sampled_level->source_field = "wind_dir";

    plotfile2d::fill_sampled_level_component(dst, 0, speed_descriptor,
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 2, wind_sources);
    plotfile2d::fill_sampled_level_component(dst, 1, dir_descriptor,
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 2, wind_sources);
    amrex::Gpu::streamSynchronize();

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(0.0));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(0.0));
    EXPECT_DOUBLE_EQ(dst.min(1), amrex::Real(-999.0));
    EXPECT_DOUBLE_EQ(dst.max(1), amrex::Real(-999.0));
}

// Motivation: Wind direction is circular, so ERF should interpolate the wind
// vector and derive direction from the interpolated vector.
TEST(Plotfile2DInterpolator, SampledWindsInterpolateVectorComponents)
{
    const auto cell_ba = make_test_boxarray();
    const auto xba = make_test_xface_boxarray();
    const auto yba = make_test_yface_boxarray();
    const auto zba = make_test_zface_boxarray();
    amrex::DistributionMapping dm(cell_ba);
    amrex::MultiFab cons(cell_ba, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab z_phys_cc(cell_ba, dm, 1, 0);
    amrex::MultiFab z_phys_nd(make_test_nodal_boxarray(), dm, 1, 0);
    amrex::MultiFab xvel(xba, dm, 1, 0);
    amrex::MultiFab yvel(yba, dm, 1, 0);
    amrex::MultiFab zvel(zba, dm, 1, 0);
    amrex::MultiFab dst(cell_ba, dm, 5, 0);

    cons.setVal(amrex::Real(0.0));
    z_phys_cc.setVal(amrex::Real(0.0));
    z_phys_nd.setVal(amrex::Real(0.0));
    dst.setVal(amrex::Real(-7.0));

    set_component_value_at_k(xvel, 0, 0, amrex::Real(0.0));
    set_component_value_at_k(xvel, 0, 1, amrex::Real(2.0));
    set_component_value_at_k(yvel, 0, 0, amrex::Real(2.0));
    set_component_value_at_k(yvel, 0, 1, amrex::Real(0.0));
    zvel.setVal(amrex::Real(3.0));
    set_component_value_at_k(z_phys_cc, 0, 0, amrex::Real(0.0));
    set_component_value_at_k(z_phys_cc, 0, 1, amrex::Real(1000.0));

    const auto wind_sources = plotfile2d::SampledWindSources{&xvel, &yvel, &zvel, nullptr, nullptr};

    Plotfile2DOutputDescriptor base_descriptor;
    base_descriptor.category = DiagnosticCategory::SampledLevel;
    base_descriptor.missing_policy = MissingPolicy::FillMinus999WhenUnavailable;
    base_descriptor.missing_value = amrex::Real(-999.0);
    base_descriptor.sampled_level = SampledLevelMetadata{
        "upper_air_winds",
        "u_east",
        SampledVerticalCoordinateMetadata{
            "height_msl",
            amrex::Real(500.0),
            "m",
            amrex::Real(500.0),
            "m",
            "linear"
        }
    };

    auto make_descriptor = [&](const char* field_name, const char* units) {
        Plotfile2DOutputDescriptor descriptor = base_descriptor;
        descriptor.name = std::string(field_name) + "_z_msl_500m";
        descriptor.long_name = field_name;
        descriptor.units = units;
        descriptor.sampled_level->source_field = field_name;
        return descriptor;
    };

    plotfile2d::fill_sampled_level_component(dst, 0, make_descriptor("u_east", "m/s"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 1, wind_sources);
    plotfile2d::fill_sampled_level_component(dst, 1, make_descriptor("v_north", "m/s"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 1, wind_sources);
    plotfile2d::fill_sampled_level_component(dst, 2, make_descriptor("w", "m/s"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 1, wind_sources);
    plotfile2d::fill_sampled_level_component(dst, 3, make_descriptor("wind_speed", "m/s"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 1, wind_sources);
    plotfile2d::fill_sampled_level_component(dst, 4, make_descriptor("wind_dir", "degrees"),
                                             cons, &z_phys_cc, z_phys_nd, true,
                                             MoistureComponentIndices(), 0, 1, wind_sources);
    amrex::Gpu::streamSynchronize();

    const amrex::Real pi = amrex::Real(3.141592653589793238462643383279502884L);
    const amrex::Real expected_dir =
        amrex::Real(270.0) - std::atan2(amrex::Real(1.0), amrex::Real(1.0)) * amrex::Real(180.0) / pi;

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(1.0));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(1.0));
    EXPECT_DOUBLE_EQ(dst.min(1), amrex::Real(1.0));
    EXPECT_DOUBLE_EQ(dst.max(1), amrex::Real(1.0));
    EXPECT_DOUBLE_EQ(dst.min(2), amrex::Real(3.0));
    EXPECT_DOUBLE_EQ(dst.max(2), amrex::Real(3.0));
    EXPECT_DOUBLE_EQ(dst.min(3), std::sqrt(amrex::Real(2.0)));
    EXPECT_DOUBLE_EQ(dst.max(3), std::sqrt(amrex::Real(2.0)));
    EXPECT_NEAR(dst.min(4), expected_dir, 1.0e-12);
    EXPECT_NEAR(dst.max(4), expected_dir, 1.0e-12);
}

// Motivation: Sampled wind fields should use the same deterministic naming and
// descriptor metadata path as thermodynamic sampled fields.
TEST(Plotfile2DMetadata, FormatsSampledWindLevelMetadata)
{
    SampledLevelDefinition upper_air;
    upper_air.name = "upper_air_winds";
    upper_air.coordinate = SampledCoordinate::Pressure;
    upper_air.units = "hPa";
    upper_air.values = {amrex::Real(850.0)};
    upper_air.fields = {"u_east", "v_north", "wind_speed", "wind_dir"};
    upper_air.interpolation = SampledInterpolation::Linear;

    const auto descriptors =
        build_sampled_level_output_descriptors_from_definitions(
            amrex::Vector<SampledLevelDefinition>{upper_air},
            amrex::Vector<std::string>{},
            SolverChoice());

    ASSERT_EQ(descriptors.size(), 4u);
    EXPECT_EQ(descriptors[0].name, "u_east_p_850hPa");
    EXPECT_EQ(descriptors[1].name, "v_north_p_850hPa");
    EXPECT_EQ(descriptors[2].name, "wind_speed_p_850hPa");
    EXPECT_EQ(descriptors[3].name, "wind_dir_p_850hPa");
    EXPECT_EQ(descriptors[0].units, "m/s");
    EXPECT_EQ(descriptors[1].units, "m/s");
    EXPECT_EQ(descriptors[2].units, "m/s");
    EXPECT_EQ(descriptors[3].units, "degrees");
    EXPECT_EQ(descriptors[0].category, DiagnosticCategory::SampledLevel);
    ASSERT_TRUE(descriptors[3].sampled_level.has_value());
    EXPECT_EQ(descriptors[3].sampled_level->source_field, "wind_dir");
    EXPECT_EQ(descriptors[3].sampled_level->vertical_coordinate.canonical_value, amrex::Real(85000.0));
    EXPECT_EQ(descriptors[3].sampled_level->vertical_coordinate.canonical_units, "Pa");
}

// Motivation: Water-path diagnostics should disappear from the available list
// when the active moisture model does not expose the corresponding conserved
// species.
TEST(Plotfile2DWaterPath, NoMoistureExcludesCondensedWaterPaths)
{
    SolverChoice sc;
    sc.moisture_type = MoistureType::None;

    const auto available = plotfile2d::available_diagnostic_names(sc);

    for (const char* name : {"integrated_qc", "integrated_qi",
                             "integrated_qr", "integrated_qs",
                             "integrated_qg"}) {
        EXPECT_EQ(std::find(available.begin(), available.end(), name), available.end());
    }
}

// Motivation: Scheme-aware availability should follow the active moisture
// indices rather than a hard-coded microphysics switch.
TEST(Plotfile2DWaterPath, KesslerAvailabilityMatchesMoistureIndices)
{
    SolverChoice sc;
    sc.moisture_type = MoistureType::Kessler;
    sc.moisture_indices = MoistureComponentIndices(RhoQ1_comp, RhoQ2_comp, -1, RhoQ3_comp);

    const auto available = plotfile2d::available_diagnostic_names(sc);

    EXPECT_NE(std::find(available.begin(), available.end(), "integrated_qc"), available.end());
    EXPECT_NE(std::find(available.begin(), available.end(), "integrated_qr"), available.end());
    EXPECT_EQ(std::find(available.begin(), available.end(), "integrated_qi"), available.end());
    EXPECT_EQ(std::find(available.begin(), available.end(), "integrated_qs"), available.end());
    EXPECT_EQ(std::find(available.begin(), available.end(), "integrated_qg"), available.end());
}

// Motivation: Full bulk schemes with cloud liquid, cloud ice, rain, snow, and
// graupel should expose all five condensed water paths.
TEST(Plotfile2DWaterPath, FullIceAvailabilityIncludesAllCondensedWaterPaths)
{
    SolverChoice sc;
    sc.moisture_type = MoistureType::Morrison;
    sc.moisture_indices = MoistureComponentIndices(RhoQ1_comp, RhoQ2_comp,
                                                   RhoQ3_comp, RhoQ4_comp,
                                                   RhoQ5_comp, RhoQ6_comp);

    const auto available = plotfile2d::available_diagnostic_names(sc);

    for (const char* name : {"integrated_qc", "integrated_qi",
                             "integrated_qr", "integrated_qs",
                             "integrated_qg"}) {
        EXPECT_NE(std::find(available.begin(), available.end(), name), available.end());
    }

    const auto descriptors = plotfile2d::active_condensed_water_path_descriptors(sc);
    ASSERT_EQ(descriptors.size(), 5u);
    EXPECT_EQ(descriptors[0].name, std::string("integrated_qc"));
    EXPECT_EQ(descriptors[0].source_component, RhoQ2_comp);
    EXPECT_EQ(descriptors[1].name, std::string("integrated_qi"));
    EXPECT_EQ(descriptors[1].source_component, RhoQ3_comp);
    EXPECT_EQ(descriptors[2].name, std::string("integrated_qr"));
    EXPECT_EQ(descriptors[2].source_component, RhoQ4_comp);
    EXPECT_EQ(descriptors[3].name, std::string("integrated_qs"));
    EXPECT_EQ(descriptors[3].source_component, RhoQ5_comp);
    EXPECT_EQ(descriptors[4].name, std::string("integrated_qg"));
    EXPECT_EQ(descriptors[4].source_component, RhoQ6_comp);
}

// Motivation: The water-path helper must reduce several selected condensed
// species in one pass without changing unrelated output components.
TEST(Plotfile2DWaterPath, FillCondensedWaterPathsUsesColumnSumAndMetric)
{
    const auto ba3d = make_test_boxarray();
    const auto ba2d = make_test_slab_boxarray();
    amrex::DistributionMapping dm(ba3d);
    // Allocate by named ERF component index, not by a literal count. The RhoQ*
    // offsets depend on the conserved-state layout.
    amrex::MultiFab src(ba3d, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab dst(ba2d, dm, 4, 0);
    amrex::MultiFab detJ(ba3d, dm, 1, 0);
    const amrex::Geometry geom = make_test_geometry();

    dst.setVal(amrex::Real(7.0));
    src.setVal(amrex::Real(0.0));
    detJ.setVal(amrex::Real(2.0));

    set_component_value_at_k(src, RhoQ2_comp, 0, amrex::Real(1.0));
    set_component_value_at_k(src, RhoQ2_comp, 1, amrex::Real(1.0));
    set_component_value_at_k(src, RhoQ2_comp, 2, amrex::Real(1.0));
    set_component_value_at_k(src, RhoQ3_comp, 0, amrex::Real(3.0));
    set_component_value_at_k(src, RhoQ3_comp, 1, amrex::Real(3.0));
    set_component_value_at_k(src, RhoQ3_comp, 2, amrex::Real(3.0));

    SolverChoice sc;
    sc.moisture_type = MoistureType::Kessler;
    sc.moisture_indices = MoistureComponentIndices(RhoQ1_comp, RhoQ2_comp, -1, RhoQ3_comp);

    const amrex::Vector<std::string> plot_var_names{
        "integrated_qv", "integrated_qc", "integrated_qr"
    };
    const auto selected =
        plotfile2d::selected_condensed_water_path_components(plot_var_names, sc);

    ASSERT_EQ(selected.n, 2);
    EXPECT_EQ(selected.dst_comp[0], 1);
    EXPECT_EQ(selected.dst_comp[1], 2);
    EXPECT_EQ(selected.src_comp[0], RhoQ2_comp);
    EXPECT_EQ(selected.src_comp[1], RhoQ3_comp);

    const auto old_mesh_type = SolverChoice::mesh_type;
    SolverChoice::set_mesh_type(MeshType::VariableDz);
    plotfile2d::fill_condensed_water_paths(dst, src, selected, geom, detJ);
    SolverChoice::set_mesh_type(old_mesh_type);

    amrex::Gpu::streamSynchronize();

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(7.0));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(7.0));
    EXPECT_DOUBLE_EQ(dst.min(1), amrex::Real(6.0));
    EXPECT_DOUBLE_EQ(dst.max(1), amrex::Real(6.0));
    EXPECT_DOUBLE_EQ(dst.min(2), amrex::Real(18.0));
    EXPECT_DOUBLE_EQ(dst.max(2), amrex::Real(18.0));
}

// Motivation: Downstream tools rely on the sidecar for units, category, and
// missing-value policy for new water-path diagnostics.
TEST(Plotfile2DMetadata, FormatsWaterPathDiagnostics)
{
    const amrex::Vector<std::string> varnames{"integrated_qc"};
    const std::string json = plotfile2d::format_2d_metadata_json(varnames);

    EXPECT_NE(json.find("\"name\": \"integrated_qc\""), std::string::npos);
    EXPECT_NE(json.find("\"long_name\": \"Column-integrated cloud liquid water\""), std::string::npos);
    EXPECT_NE(json.find("\"units\": \"kg/m^2\""), std::string::npos);
    EXPECT_NE(json.find("\"category\": \"ColumnIntegral\""), std::string::npos);
    EXPECT_NE(json.find("\"missing_policy\": \"FillZeroWhenUnavailable\""), std::string::npos);
    EXPECT_NE(json.find("\"missing_value\": 0"), std::string::npos);
}

// Motivation: The W m^-2 diagnostics are public 2D outputs, so their catalog
// entries must expose the intended units, category, and missing-value policy.
TEST(Plotfile2D, FindDiagnosticReturnsDescriptorForSurfaceFluxComposition)
{
    const auto* sensible = plotfile2d::find_diagnostic("sensible_heat_flux");
    ASSERT_NE(sensible, nullptr);
    EXPECT_STREQ(sensible->name, "sensible_heat_flux");
    EXPECT_EQ(sensible->id, plotfile2d::DiagnosticID::SensibleHeatFlux);
    EXPECT_EQ(sensible->category, plotfile2d::DiagnosticCategory::SurfaceFlux);
    EXPECT_STREQ(sensible->units, "W m^-2");
    EXPECT_EQ(sensible->missing_policy, plotfile2d::MissingPolicy::FillMinus999WhenUnavailable);

    const auto* latent = plotfile2d::find_diagnostic("latent_heat_flux");
    ASSERT_NE(latent, nullptr);
    EXPECT_STREQ(latent->name, "latent_heat_flux");
    EXPECT_EQ(latent->id, plotfile2d::DiagnosticID::LatentHeatFlux);
    EXPECT_EQ(latent->category, plotfile2d::DiagnosticCategory::SurfaceFlux);
    EXPECT_STREQ(latent->units, "W m^-2");
    EXPECT_EQ(latent->missing_policy, plotfile2d::MissingPolicy::FillMinus999WhenUnavailable);
}

// Motivation: Native SHOC state_update preserves consumed flux snapshots, but
// the 2D writer must not turn an absent host latent-flux field into a zero
// SHOC snapshot. Missing flux components should still fill -999.
TEST(Plotfile2D, NativeShocConsumedFluxSelectionRequiresHostFluxField)
{
    EXPECT_TRUE(plotfile2d::use_native_shoc_consumed_flux_source(true, true, true));
    EXPECT_FALSE(plotfile2d::use_native_shoc_consumed_flux_source(true, true, false));
    EXPECT_FALSE(plotfile2d::use_native_shoc_consumed_flux_source(true, false, true));
    EXPECT_FALSE(plotfile2d::use_native_shoc_consumed_flux_source(false, true, true));
    EXPECT_FALSE(plotfile2d::use_native_shoc_consumed_flux_source(false, false, false));
}

// Motivation: Sensible and latent fluxes are selected independently, so one
// missing host field must not force the other component to use the wrong
// source.
TEST(Plotfile2D, NativeShocConsumedFluxSelectionIsComponentSpecific)
{
    const bool native_shoc_owns_scalar_fluxes = true;
    const bool native_shoc_has_consumed_flux_diagnostics = true;
    const bool host_sensible_flux_available = true;
    const bool host_latent_flux_available = false;

    const bool use_sensible =
        plotfile2d::use_native_shoc_consumed_flux_source(
            native_shoc_owns_scalar_fluxes,
            native_shoc_has_consumed_flux_diagnostics,
            host_sensible_flux_available);

    const bool use_latent =
        plotfile2d::use_native_shoc_consumed_flux_source(
            native_shoc_owns_scalar_fluxes,
            native_shoc_has_consumed_flux_diagnostics,
            host_latent_flux_available);

    EXPECT_TRUE(use_sensible);
    EXPECT_FALSE(use_latent);
}

// Motivation: Native SHOC diagnostics are public 2D outputs, so their catalog
// entries must expose the intended category and missing-value policy.
TEST(Plotfile2D, FindDiagnosticReturnsDescriptorForNativeShocSurfaceDiagnostics)
{
    const auto* ustar = plotfile2d::find_diagnostic("shoc_u_star");
    ASSERT_NE(ustar, nullptr);
    EXPECT_STREQ(ustar->name, "shoc_u_star");
    EXPECT_EQ(ustar->id, plotfile2d::DiagnosticID::ShocUStar);
    EXPECT_EQ(ustar->category, plotfile2d::DiagnosticCategory::PBL);
    EXPECT_STREQ(ustar->units, "m/s");
    EXPECT_EQ(ustar->missing_policy, plotfile2d::MissingPolicy::FillMinus999WhenUnavailable);

    const auto* olen = plotfile2d::find_diagnostic("shoc_Olen");
    ASSERT_NE(olen, nullptr);
    EXPECT_STREQ(olen->name, "shoc_Olen");
    EXPECT_EQ(olen->id, plotfile2d::DiagnosticID::ShocOlen);
    EXPECT_EQ(olen->category, plotfile2d::DiagnosticCategory::PBL);
    EXPECT_STREQ(olen->units, "m");
    EXPECT_EQ(olen->missing_policy, plotfile2d::MissingPolicy::FillMinus999WhenUnavailable);

    const auto* wthv = plotfile2d::find_diagnostic("shoc_wthv_sfc");
    ASSERT_NE(wthv, nullptr);
    EXPECT_STREQ(wthv->name, "shoc_wthv_sfc");
    EXPECT_EQ(wthv->id, plotfile2d::DiagnosticID::ShocWthvSfc);
    EXPECT_EQ(wthv->category, plotfile2d::DiagnosticCategory::PBL);
    EXPECT_STREQ(wthv->units, "K m s^-1");
    EXPECT_EQ(wthv->missing_policy, plotfile2d::MissingPolicy::FillMinus999WhenUnavailable);
}

// Motivation: Unknown catalog lookups should fail cleanly so callers can
// distinguish a typo from a valid built-in diagnostic.
TEST(Plotfile2D, FindDiagnosticReturnsNullForUnknownName)
{
    EXPECT_EQ(plotfile2d::find_diagnostic("not_a_real_2d_name"), nullptr);
}

// Motivation: The selection helper must follow the catalog order, not the user
// request order, so equivalent input decks produce the same component layout.
TEST(Plotfile2D, SelectionUsesCatalogCanonicalOrder)
{
    const auto available = plotfile2d::diagnostic_names();
    const amrex::Vector<std::string> requested{"surf_pres", "z0", "z_surf", "lon_m"};

    const auto selection = select_requested_plot_variables(requested, available);

    const amrex::Vector<std::string> expected{"z_surf", "lon_m", "z0", "surf_pres"};
    EXPECT_EQ(selection.accepted, expected);
}

// Motivation: A warning must name both the input parameter and the skipped
// variable so the user can fix the correct namelist entry.
TEST(Plotfile2D, WarningMessageNamesParameterAndVariable)
{
    const amrex::Vector<std::string> available{"z_surf", "landmask", "mapfac"};
    const std::string message =
        format_unavailable_2d_plot_var_warning("erf.plot2d_vars_1", "bogus", available);

    EXPECT_TRUE(contains(message, "erf.plot2d_vars_1"));
    EXPECT_TRUE(contains(message, "bogus"));
    EXPECT_TRUE(contains(message, "2D plot variable"));
    EXPECT_TRUE(contains(message, "skipped"));
    EXPECT_TRUE(contains(message, "z_surf"));
}

// Motivation: Warnings should report the full ParmParse name so users can
// locate the exact plot2d namelist entry in a multi-stream input deck.
TEST(Plotfile2D, ParameterNameIncludesParmParsePrefix)
{
    EXPECT_EQ(format_plot2d_parameter_name("erf", "plot2d_vars_1"), "erf.plot2d_vars_1");
}

// Motivation: Helper formatting should not invent a separator when the prefix
// is empty, which keeps standalone call sites readable.
TEST(Plotfile2D, ParameterNameHandlesEmptyPrefix)
{
    EXPECT_EQ(format_plot2d_parameter_name("", "plot2d_vars_1"), "plot2d_vars_1");
}

// Motivation: The internal component-count guard needs to expose enough
// context to debug a drift between the canonical variable list and the fill
// blocks.
TEST(Plotfile2D, ComponentMismatchMessageReportsCounts)
{
    const std::string message = format_2d_component_count_error(3, 5, 7);

    EXPECT_TRUE(contains(message, "internal error"));
    EXPECT_TRUE(contains(message, "level 3"));
    EXPECT_TRUE(contains(message, "filled 5"));
    EXPECT_TRUE(contains(message, "expected 7"));
    EXPECT_TRUE(contains(message, "inconsistent"));
}

// Motivation: Invalid stream indices should fail with a clear developer-facing
// message instead of silently choosing a file prefix.
TEST(Plotfile2D, InvalidStreamMessageReportsAllowedValues)
{
    const std::string message = format_invalid_2d_stream_error(0);

    EXPECT_TRUE(contains(message, "invalid stream index 0"));
    EXPECT_TRUE(contains(message, "expected 1 or 2"));
}

// Motivation: Missing-value fills should not corrupt neighboring output
// components in the 2D plotfile MultiFab.
TEST(Plotfile2DFill, FillComponentWithValueTouchesOnlyRequestedComponent)
{
    const auto ba = make_test_boxarray();
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab dst(ba, dm, 2, 0);

    dst.setVal(amrex::Real(1.5), 0, 1, 0);
    dst.setVal(amrex::Real(2.5), 1, 1, 0);

    plotfile2d::fill_component_with_value(dst, 1, amrex::Real(-999.0));
    amrex::Gpu::streamSynchronize();

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(1.5));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(1.5));
    EXPECT_DOUBLE_EQ(dst.min(1), amrex::Real(-999.0));
    EXPECT_DOUBLE_EQ(dst.max(1), amrex::Real(-999.0));
}

// Motivation: Surface, top-of-column, and future interpolated diagnostics rely
// on copying the intended source level into the 2D output component.
TEST(Plotfile2DFill, FillComponentFromKLevelUsesRequestedSourceLevelAndComponent)
{
    const auto ba = make_test_boxarray();
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab src(ba, dm, 2, 0);
    amrex::MultiFab dst(ba, dm, 1, 0);

    src.setVal(amrex::Real(-1.0));
    set_component_value_at_k(src, 1, 2, amrex::Real(23.5));

    plotfile2d::fill_component_from_klevel(dst, 0, src, 2, 1);
    amrex::Gpu::streamSynchronize();

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(23.5));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(23.5));
}

// Motivation: Optional diagnostics such as SurfaceLayer, radiation, and SHOC
// fields must write documented missing values when their source is unavailable.
TEST(Plotfile2DFill, FillComponentFromKLevelOrValueUsesMissingValueWhenSourceIsNull)
{
    const auto ba = make_test_boxarray();
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab dst(ba, dm, 1, 0);

    dst.setVal(amrex::Real(4.0));

    plotfile2d::fill_component_from_klevel_or_value(dst, 0, nullptr, 0, amrex::Real(-999.0));
    amrex::Gpu::streamSynchronize();

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(-999.0));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(-999.0));
}

// Motivation: The same helper must handle available and unavailable diagnostic
// sources without changing the component order logic in the writer.
TEST(Plotfile2DFill, FillComponentFromKLevelOrValueCopiesWhenSourceExists)
{
    const auto ba = make_test_boxarray();
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab src(ba, dm, 1, 0);
    amrex::MultiFab dst(ba, dm, 1, 0);

    src.setVal(amrex::Real(-2.0));
    set_component_value_at_k(src, 0, 1, amrex::Real(17.25));

    plotfile2d::fill_component_from_klevel_or_value(dst, 0, &src, 1, amrex::Real(-999.0));
    amrex::Gpu::streamSynchronize();

    EXPECT_DOUBLE_EQ(dst.min(0), amrex::Real(17.25));
    EXPECT_DOUBLE_EQ(dst.max(0), amrex::Real(17.25));
}

// Motivation: The 2D assembly refactor must preserve the W m^-2 conversion
// used by the public sensible_heat_flux diagnostic.
TEST(Plotfile2DFill, SensibleHeatFluxFillUsesSurfaceFluxConversion)
{
    const auto ba = make_test_boxarray();
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab src(ba, dm, 1, 0);
    amrex::MultiFab dst(ba, dm, 1, 0);

    src.setVal(amrex::Real(-7.0));
    set_component_value_at_k(src, 0, 1, amrex::Real(1.75));

    plotfile2d::fill_sensible_heat_flux_from_klevel_or_missing(
        dst, 0, &src, 1, amrex::Real(-999.0));
    amrex::Gpu::streamSynchronize();

    const amrex::Real expected =
        surface_flux_diagnostics::sensible_heat_flux_wm2_from_rhotheta_flux(amrex::Real(1.75));
    EXPECT_DOUBLE_EQ(dst.min(0), expected);
    EXPECT_DOUBLE_EQ(dst.max(0), expected);
}

// Motivation: The 2D assembly refactor must preserve the W m^-2 conversion
// used by the public latent_heat_flux diagnostic.
TEST(Plotfile2DFill, LatentHeatFluxFillUsesSurfaceFluxConversion)
{
    const auto ba = make_test_boxarray();
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab src(ba, dm, 1, 0);
    amrex::MultiFab dst(ba, dm, 1, 0);

    src.setVal(amrex::Real(-7.0));
    set_component_value_at_k(src, 0, 2, amrex::Real(-0.5));

    plotfile2d::fill_latent_heat_flux_from_klevel_or_missing(
        dst, 0, &src, 2, amrex::Real(-999.0));
    amrex::Gpu::streamSynchronize();

    const amrex::Real expected =
        surface_flux_diagnostics::latent_heat_flux_wm2_from_rhoqv_flux(amrex::Real(-0.5));
    EXPECT_DOUBLE_EQ(dst.min(0), expected);
    EXPECT_DOUBLE_EQ(dst.max(0), expected);
}

// Motivation: The JSON sidecar is a public metadata contract, so diagnostic
// categories must serialize to stable, readable strings.
TEST(Plotfile2DMetadata, CategoryStringsMatchPublicMetadataSchema)
{
    EXPECT_STREQ(diagnostic_category_to_string(DiagnosticCategory::Geometry), "Geometry");
    EXPECT_STREQ(diagnostic_category_to_string(DiagnosticCategory::SurfaceLayer), "SurfaceLayer");
    EXPECT_STREQ(diagnostic_category_to_string(DiagnosticCategory::Radiation), "Radiation");
    EXPECT_STREQ(diagnostic_category_to_string(DiagnosticCategory::SurfaceFlux), "SurfaceFlux");
    EXPECT_STREQ(diagnostic_category_to_string(DiagnosticCategory::PBL), "PBL");
    EXPECT_STREQ(diagnostic_category_to_string(DiagnosticCategory::SurfaceState), "SurfaceState");
    EXPECT_STREQ(diagnostic_category_to_string(DiagnosticCategory::ColumnIntegral), "ColumnIntegral");
    EXPECT_STREQ(diagnostic_category_to_string(DiagnosticCategory::SampledLevel), "SampledLevel");
}

// Motivation: Downstream readers need to distinguish always-available fields
// from fields that use zero or -999 as their documented unavailable value.
TEST(Plotfile2DMetadata, MissingPolicyStringsAndValuesMatchPublicMetadataSchema)
{
    EXPECT_STREQ(missing_policy_to_string(MissingPolicy::AlwaysAvailable), "AlwaysAvailable");
    EXPECT_EQ(missing_value_json(MissingPolicy::AlwaysAvailable), "null");

    EXPECT_STREQ(missing_policy_to_string(MissingPolicy::FillZeroWhenUnavailable), "FillZeroWhenUnavailable");
    EXPECT_EQ(missing_value_json(MissingPolicy::FillZeroWhenUnavailable), "0");

    EXPECT_STREQ(missing_policy_to_string(MissingPolicy::FillMinus999WhenUnavailable),
                 "FillMinus999WhenUnavailable");
    EXPECT_EQ(missing_value_json(MissingPolicy::FillMinus999WhenUnavailable), "-999");
}

// Motivation: Metadata strings come from catalog descriptors. The sidecar must
// remain valid JSON if future descriptors contain quotes or control characters.
TEST(Plotfile2DMetadata, EscapesJsonStrings)
{
    const std::string input = "\"quoted\" slash\\ line\n tab\t return\r";
    const std::string expected = "\\\"quoted\\\" slash\\\\ line\\n tab\\t return\\r";

    EXPECT_EQ(escape_json_string(input), expected);
}

// Motivation: The sidecar must describe plotfile components by output index,
// so its variable order must match the selected plotfile component order.
TEST(Plotfile2DMetadata, FormatsSelectedVariablesInOutputOrder)
{
    const amrex::Vector<std::string> selected{
        "z_surf", "pblh", "latent_heat_flux", "shoc_wthv_sfc"
    };

    const std::string json = format_2d_metadata_json(selected);

    EXPECT_TRUE(contains(json, "\"format_version\": 2"));
    EXPECT_TRUE(contains(json, "\"kind\": \"ERF 2D plotfile metadata\""));
    EXPECT_TRUE(contains(json, "\"n_variables\": 4"));
    EXPECT_TRUE(contains(json, "\"component_index\": 0"));
    EXPECT_TRUE(contains(json, "\"component_index\": 1"));
    EXPECT_TRUE(contains(json, "\"component_index\": 2"));
    EXPECT_TRUE(contains(json, "\"component_index\": 3"));
    EXPECT_TRUE(contains(json, "\"name\": \"z_surf\""));
    EXPECT_TRUE(contains(json, "\"name\": \"pblh\""));
    EXPECT_TRUE(contains(json, "\"name\": \"latent_heat_flux\""));
    EXPECT_TRUE(contains(json, "\"name\": \"shoc_wthv_sfc\""));
    EXPECT_TRUE(contains(json, "\"long_name\": \"Surface elevation\""));
    EXPECT_TRUE(contains(json, "\"units\": \"m\""));
    EXPECT_TRUE(contains(json, "\"category\": \"Geometry\""));
    EXPECT_TRUE(contains(json, "\"category\": \"SurfaceLayer\""));
    EXPECT_TRUE(contains(json, "\"category\": \"SurfaceFlux\""));
    EXPECT_TRUE(contains(json, "\"category\": \"PBL\""));
    EXPECT_TRUE(contains(json, "\"missing_policy\": \"AlwaysAvailable\""));
    EXPECT_TRUE(contains(json, "\"missing_policy\": \"FillMinus999WhenUnavailable\""));
    EXPECT_TRUE(contains(json, "\"missing_value\": null"));
    EXPECT_TRUE(contains(json, "\"missing_value\": -999"));

    const auto z_pos = json.find("\"name\": \"z_surf\"");
    const auto pblh_pos = json.find("\"name\": \"pblh\"");
    const auto latent_pos = json.find("\"name\": \"latent_heat_flux\"");
    const auto shoc_wthv_pos = json.find("\"name\": \"shoc_wthv_sfc\"");
    ASSERT_NE(z_pos, std::string::npos);
    ASSERT_NE(pblh_pos, std::string::npos);
    ASSERT_NE(latent_pos, std::string::npos);
    ASSERT_NE(shoc_wthv_pos, std::string::npos);
    EXPECT_LT(z_pos, pblh_pos);
    EXPECT_LT(pblh_pos, latent_pos);
    EXPECT_LT(latent_pos, shoc_wthv_pos);
}

// Motivation: The sidecar must describe the components written to this
// plotfile, not every possible catalog entry.
TEST(Plotfile2DMetadata, FormatsOnlySelectedVariables)
{
    const amrex::Vector<std::string> selected{"surf_pres"};
    const std::string json = format_2d_metadata_json(selected);

    EXPECT_TRUE(contains(json, "\"n_variables\": 1"));
    EXPECT_TRUE(contains(json, "\"name\": \"surf_pres\""));
    EXPECT_FALSE(contains(json, "\"name\": \"z_surf\""));
}

// Motivation: Adding new catalog entries should fail fast if the metadata
// sidecar cannot serialize their category, missing policy, or descriptor
// fields.
TEST(Plotfile2DMetadata, FormatsEveryCatalogDiagnostic)
{
    const auto names = plotfile2d::diagnostic_names();
    const std::string json = format_2d_metadata_json(names);

    EXPECT_TRUE(contains(json, std::string("\"n_variables\": ") + std::to_string(names.size())));
    EXPECT_TRUE(contains(json, "\"name\": \"shoc_u_star\""));
    EXPECT_TRUE(contains(json, "\"name\": \"shoc_Olen\""));
    EXPECT_TRUE(contains(json, "\"name\": \"shoc_wthv_sfc\""));
    EXPECT_TRUE(contains(json, "\"category\": \"PBL\""));
    for (const auto& name : names) {
        EXPECT_TRUE(contains(json, "\"name\": \"" + name + "\""));
    }
}

// Motivation: Sampled-level metadata must preserve the sampled source field
// and vertical-coordinate contract without looking up dynamic names in the
// static catalog.
TEST(Plotfile2DMetadata, FormatsSampledLevelVerticalCoordinate)
{
    Plotfile2DOutputDescriptor descriptor;
    descriptor.name = "theta_p_850hPa";
    descriptor.long_name = "Potential temperature sampled on pressure levels";
    descriptor.units = "K";
    descriptor.category = plotfile2d::DiagnosticCategory::SampledLevel;
    descriptor.missing_policy = plotfile2d::MissingPolicy::FillMinus999WhenUnavailable;
    descriptor.missing_value = amrex::Real(-999.0);
    descriptor.sampled_level = plotfile2d::SampledLevelMetadata{
        "upper_air",
        "theta",
        plotfile2d::SampledVerticalCoordinateMetadata{
            "pressure",
            amrex::Real(850.0),
            "hPa",
            amrex::Real(85000.0),
            "Pa",
            "linear"
        }
    };

    const amrex::Vector<Plotfile2DOutputDescriptor> descriptors{descriptor};
    const std::string json = format_2d_metadata_json(descriptors);

    EXPECT_TRUE(contains(json, "\"format_version\": 2"));
    EXPECT_TRUE(contains(json, "\"name\": \"theta_p_850hPa\""));
    EXPECT_TRUE(contains(json, "\"category\": \"SampledLevel\""));
    EXPECT_TRUE(contains(json, "\"source_field\": \"theta\""));
    EXPECT_TRUE(contains(json, "\"vertical_coordinate\""));
    EXPECT_TRUE(contains(json, "\"type\": \"pressure\""));
    EXPECT_TRUE(contains(json, "\"value\": 850"));
    EXPECT_TRUE(contains(json, "\"units\": \"hPa\""));
    EXPECT_TRUE(contains(json, "\"canonical_value\": 85000"));
    EXPECT_TRUE(contains(json, "\"canonical_units\": \"Pa\""));
    EXPECT_TRUE(contains(json, "\"interpolation\": \"linear\""));
}

// Motivation: Static diagnostics remain on metadata version 2 so existing
// plotfile readers see the same public catalog fields.
TEST(Plotfile2DMetadata, StaticDiagnosticsRemainMetadataVersionTwo)
{
    const std::string json = format_2d_metadata_json(amrex::Vector<std::string>{"z_surf"});

    EXPECT_TRUE(contains(json, "\"format_version\": 2"));
    EXPECT_TRUE(contains(json, "\"name\": \"z_surf\""));
    EXPECT_FALSE(contains(json, "\"source_field\""));
}

// Motivation: Sampled-level metadata must preserve the output component
// order that the writer assembles.
TEST(Plotfile2DMetadata, SampledLevelMetadataPreservesOutputOrder)
{
    SampledLevelDefinition upper_air;
    upper_air.name = "upper_air";
    upper_air.coordinate = SampledCoordinate::Pressure;
    upper_air.units = "hPa";
    upper_air.values = {amrex::Real(850.0), amrex::Real(700.0)};
    upper_air.fields = {"theta", "qv"};
    upper_air.interpolation = SampledInterpolation::Linear;

    SolverChoice sc;
    sc.moisture_type = MoistureType::Kessler;
    sc.moisture_indices = MoistureComponentIndices(RhoQ1_comp, RhoQ2_comp, -1, RhoQ3_comp);

    const auto descriptors =
        build_sampled_level_output_descriptors_from_definitions(
            amrex::Vector<SampledLevelDefinition>{upper_air},
            amrex::Vector<std::string>{"z_surf"},
            sc);

    ASSERT_EQ(descriptors.size(), 5u);
    EXPECT_EQ(descriptors[0].name, "z_surf");
    EXPECT_EQ(descriptors[1].name, "theta_p_850hPa");
    EXPECT_EQ(descriptors[2].name, "qv_p_850hPa");
    EXPECT_EQ(descriptors[3].name, "theta_p_700hPa");
    EXPECT_EQ(descriptors[4].name, "qv_p_700hPa");

    const std::string json = format_2d_metadata_json(descriptors);
    EXPECT_TRUE(contains(json, "\"component_index\": 0"));
    EXPECT_TRUE(contains(json, "\"component_index\": 4"));
    EXPECT_TRUE(contains(json, "\"source_field\": \"theta\""));
    EXPECT_TRUE(contains(json, "\"source_field\": \"qv\""));
}

// Motivation: A sampled-only 2D stream should still build a descriptor list
// with one output per requested field and target value.
TEST(Plotfile2D, SampledOnlyStreamBuildsOutputDescriptorList)
{
    SampledLevelDefinition upper_air;
    upper_air.name = "upper_air";
    upper_air.coordinate = SampledCoordinate::Pressure;
    upper_air.units = "hPa";
    upper_air.values = {amrex::Real(850.0)};
    upper_air.fields = {"theta", "temp"};
    upper_air.interpolation = SampledInterpolation::Linear;

    const auto descriptors =
        build_sampled_level_output_descriptors_from_definitions(
            amrex::Vector<SampledLevelDefinition>{upper_air},
            amrex::Vector<std::string>{},
            SolverChoice());

    ASSERT_EQ(descriptors.size(), 2u);
    EXPECT_EQ(descriptors[0].name, "theta_p_850hPa");
    EXPECT_EQ(descriptors[1].name, "temp_p_850hPa");
    EXPECT_EQ(descriptors[0].category, DiagnosticCategory::SampledLevel);
    EXPECT_EQ(descriptors[0].missing_policy, MissingPolicy::FillMinus999WhenUnavailable);
    ASSERT_TRUE(descriptors[0].sampled_level.has_value());
    EXPECT_EQ(descriptors[0].sampled_level->level_set_name, "upper_air");
    EXPECT_EQ(descriptors[0].sampled_level->source_field, "theta");
}

// Motivation: Static built-in diagnostics must stay ahead of sampled-level
// diagnostics in the assembled output descriptor list.
TEST(Plotfile2D, SampledDescriptorsAppendAfterStaticDiagnostics)
{
    SampledLevelDefinition upper_air;
    upper_air.name = "upper_air";
    upper_air.coordinate = SampledCoordinate::Pressure;
    upper_air.units = "hPa";
    upper_air.values = {amrex::Real(850.0)};
    upper_air.fields = {"theta", "temp"};
    upper_air.interpolation = SampledInterpolation::Linear;

    const auto descriptors =
        build_sampled_level_output_descriptors_from_definitions(
            amrex::Vector<SampledLevelDefinition>{upper_air},
            amrex::Vector<std::string>{"z_surf", "landmask"},
            SolverChoice());

    ASSERT_EQ(descriptors.size(), 4u);
    EXPECT_EQ(descriptors[0].name, "z_surf");
    EXPECT_EQ(descriptors[1].name, "landmask");
    EXPECT_EQ(descriptors[2].name, "theta_p_850hPa");
    EXPECT_EQ(descriptors[3].name, "temp_p_850hPa");
    EXPECT_NE(descriptors[0].static_diagnostic, nullptr);
    EXPECT_NE(descriptors[1].static_diagnostic, nullptr);
}

// Motivation: Sampled descriptors must keep level-set order, target order,
// and field order so the output layout stays deterministic.
TEST(Plotfile2D, SampledDescriptorOrderIsLevelSetTargetField)
{
    SampledLevelDefinition upper_air;
    upper_air.name = "upper_air";
    upper_air.coordinate = SampledCoordinate::Pressure;
    upper_air.units = "hPa";
    upper_air.values = {amrex::Real(850.0), amrex::Real(700.0)};
    upper_air.fields = {"theta", "temp"};
    upper_air.interpolation = SampledInterpolation::Linear;

    SampledLevelDefinition native_k;
    native_k.name = "native_k";
    native_k.coordinate = SampledCoordinate::ModelIndex;
    native_k.values = {amrex::Real(10.0)};
    native_k.fields = {"temp"};
    native_k.interpolation = SampledInterpolation::None;

    const auto descriptors =
        build_sampled_level_output_descriptors_from_definitions(
            amrex::Vector<SampledLevelDefinition>{upper_air, native_k},
            amrex::Vector<std::string>{},
            SolverChoice());

    ASSERT_EQ(descriptors.size(), 5u);
    EXPECT_EQ(descriptors[0].name, "theta_p_850hPa");
    EXPECT_EQ(descriptors[1].name, "temp_p_850hPa");
    EXPECT_EQ(descriptors[2].name, "theta_p_700hPa");
    EXPECT_EQ(descriptors[3].name, "temp_p_700hPa");
    EXPECT_EQ(descriptors[4].name, "temp_k_10");
}

// Motivation: Sampled-level output names are a public plotfile contract, so
// the coordinate and value tags must stay deterministic.
TEST(Plotfile2DSampledLevel, FormatsOutputNames)
{
    EXPECT_EQ(sampled_output_name("theta", SampledCoordinate::Pressure, amrex::Real(850.0), "hPa"),
              "theta_p_850hPa");
    EXPECT_EQ(sampled_output_name("qv", SampledCoordinate::HeightAGL, amrex::Real(500.0), "m"),
              "qv_z_agl_500m");
    EXPECT_EQ(sampled_output_name("qg", SampledCoordinate::ModelIndex, amrex::Real(10.0), "1"),
              "qg_k_10");
    EXPECT_EQ(sampled_output_name("u_east", SampledCoordinate::Pressure, amrex::Real(850.0), "hPa"),
              "u_east_p_850hPa");
    EXPECT_EQ(sampled_output_name("v_north", SampledCoordinate::Pressure, amrex::Real(850.0), "hPa"),
              "v_north_p_850hPa");
    EXPECT_EQ(sampled_output_name("w", SampledCoordinate::HeightAGL, amrex::Real(500.0), "m"),
              "w_z_agl_500m");
    EXPECT_EQ(sampled_output_name("wind_speed", SampledCoordinate::ModelIndex, amrex::Real(10.0), "1"),
              "wind_speed_k_10");
    EXPECT_EQ(sampled_output_name("wind_dir", SampledCoordinate::Pressure, amrex::Real(700.0), "hPa"),
              "wind_dir_p_700hPa");
}

// Motivation: The metadata file should travel with the native AMReX 2D
// plotfile directory and not collide with other output streams.
TEST(Plotfile2DMetadata, MetadataFilenameLivesInsidePlotfileDirectory)
{
    EXPECT_EQ(metadata_json_filename("plt2d_00010"), "plt2d_00010/2DMetadata.json");
}

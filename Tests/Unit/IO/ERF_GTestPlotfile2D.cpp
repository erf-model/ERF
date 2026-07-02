#include <string>
#include <unordered_set>

#include <gtest/gtest.h>

#include "ERF_Plotfile2DCatalog.H"
#include "ERF_Plotfile2DUtils.H"

using namespace plotfile2d;

namespace
{

bool contains (const std::string& haystack, const std::string& needle)
{
    return haystack.find(needle) != std::string::npos;
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
        "surf_pres", "integrated_qv", "surface_diagnostic_source",
        "sensible_heat_flux", "latent_heat_flux",
        "shoc_u_star", "shoc_Olen", "shoc_wthv_sfc"
    };

    EXPECT_EQ(plotfile2d::diagnostic_names(), expected);
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

#include <array>
#include <string>

#include <gtest/gtest.h>

#include "ERF_SurfaceDiagnosticSource.H"

using namespace surface_diagnostics;

// Motivation: surface_diagnostic_source writes numeric codes to plotfiles, so
// the code values must stay stable across releases.
TEST(SurfaceDiagnosticSource, PlotCodesAreStable)
{
    EXPECT_EQ(to_int(SurfaceDiagnosticSource::Missing), 0);
    EXPECT_EQ(to_int(SurfaceDiagnosticSource::SurfaceLayerLand), 1);
    EXPECT_EQ(to_int(SurfaceDiagnosticSource::LSMLand), 2);
    EXPECT_EQ(to_int(SurfaceDiagnosticSource::SurfaceLayerFallback), 3);
    EXPECT_EQ(to_int(SurfaceDiagnosticSource::SurfaceLayerSea), 4);
    EXPECT_EQ(to_int(SurfaceDiagnosticSource::Custom), 5);
    EXPECT_EQ(to_int(SurfaceDiagnosticSource::RICO), 6);
}

// Motivation: source names make docs and debugging messages easier to keep in
// sync with the numeric plotfile convention.
TEST(SurfaceDiagnosticSource, SourceNamesAreNonEmpty)
{
    const std::array<SurfaceDiagnosticSource, 7> values{
        SurfaceDiagnosticSource::Missing,
        SurfaceDiagnosticSource::SurfaceLayerLand,
        SurfaceDiagnosticSource::LSMLand,
        SurfaceDiagnosticSource::SurfaceLayerFallback,
        SurfaceDiagnosticSource::SurfaceLayerSea,
        SurfaceDiagnosticSource::Custom,
        SurfaceDiagnosticSource::RICO
    };

    for (const auto value : values) {
        EXPECT_FALSE(std::string(source_name(value)).empty());
    }
}

// Motivation: scalar source classification must distinguish valid LSM values
// from fallback values without changing the underlying flux calculation.
TEST(SurfaceDiagnosticSource, ClassifiesLsmAndFallbackPaths)
{
    EXPECT_EQ(classify_scalar_source(false, false, true,  true,  true),
              SurfaceDiagnosticSource::LSMLand);
    EXPECT_EQ(classify_scalar_source(true, false, true,  true,  true),
              SurfaceDiagnosticSource::LSMLand);
    EXPECT_EQ(classify_scalar_source(false, true, true,  true,  true),
              SurfaceDiagnosticSource::LSMLand);

    EXPECT_EQ(classify_scalar_source(false, false, true,  true,  false),
              SurfaceDiagnosticSource::SurfaceLayerFallback);
    EXPECT_EQ(classify_scalar_source(true, false, true,  true,  false),
              SurfaceDiagnosticSource::Custom);
    EXPECT_EQ(classify_scalar_source(false, true, true,  true,  false),
              SurfaceDiagnosticSource::RICO);

    EXPECT_EQ(classify_scalar_source(false, false, true,  false, false),
              SurfaceDiagnosticSource::SurfaceLayerLand);

    EXPECT_EQ(classify_scalar_source(false, false, false, false, false),
              SurfaceDiagnosticSource::SurfaceLayerSea);
}

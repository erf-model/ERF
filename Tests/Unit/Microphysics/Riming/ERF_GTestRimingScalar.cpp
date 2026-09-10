#include <gtest/gtest.h>

#include "ERF_GTestRimingCommon.H"

using namespace riming_test;

// Motivation: ERF sets the density of rime accreted onto ice from the
// Heymsfield-Pflaum (1985) parameterization and the ice surface temperature
// (Source/Particles/ERF_SuperDropletPCRiming.H, SDRiming). This checks the rime
// density bounds and impact-speed densification, the deposition-heating sign and
// monotonicity of the ice surface temperature, and the viscosity helper -- all on
// the host.
TEST(RimingScalar, RimePropertiesHoldOnHost)
{
    for (int i = 0; i < kNumMetrics; ++i) {
        EXPECT_LE(riming_metric(i), kTol) << metric_label(i);
    }
}

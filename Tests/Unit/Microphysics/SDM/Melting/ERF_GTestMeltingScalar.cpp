#include <gtest/gtest.h>

#include "ERF_GTestMeltingCommon.H"

using namespace melting_test;

// Motivation: ERF melts ice super-droplets at the Seifert-Beheng (2006, Eq. 72)
// rate (Source/Particles/ERF_SuperDropletPCMassChange.H, dMdt::meltRate). This
// checks the melt sign above 0 C, acceleration with temperature, ventilation
// enhancement, the capacitance size dependence, and the freezing-point threshold
// -- all on the host.
TEST(MeltingScalar, MeltRatePropertiesHoldOnHost)
{
    for (int i = 0; i < kNumMetrics; ++i) {
        EXPECT_LE(melting_metric(i), kTol) << metric_label(i);
    }
}

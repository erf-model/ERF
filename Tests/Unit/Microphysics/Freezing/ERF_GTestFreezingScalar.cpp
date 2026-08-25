#include <gtest/gtest.h>

#include "ERF_GTestFreezingCommon.H"

using namespace freezing_test;

// Motivation: ERF assigns each ice-nucleating super-droplet a freezing
// temperature from the Niemand et al. (2012) INAS density
// (Source/Particles/ERF_SuperDropletPCDefinitions.H, INAS_Niemand2012). This
// checks the n_s and P_freeze monotonicities, the inverse-CDF round-trip, and the
// homogeneous-freezing and clamping edge cases -- all on the host.
TEST(FreezingScalar, INASPropertiesHoldOnHost)
{
    for (int i = 0; i < kNumMetrics; ++i) {
        EXPECT_LE(freezing_metric(i), kTol) << metric_label(i);
    }
}

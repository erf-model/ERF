#include <gtest/gtest.h>

#include "ERF_GTestCondensationCommon.H"

using namespace condensation_test;

// Motivation: ERF's warm-phase super-droplet growth integrates dR^2/dt with the
// Koehler curvature and solute terms (Source/Particles/ERF_SuperDropletPCMassChange.H,
// SDMassChangeUtils_LV). This checks the diffusional r^2-law, the growth/evaporation
// sign, the Kelvin and Raoult limits, Koehler equilibrium and activation, the
// time-integrated growth, and the ventilation factor -- all on the host.
TEST(CondensationScalar, GrowthPropertiesHoldOnHost)
{
    for (int i = 0; i < kNumMetrics; ++i) {
        EXPECT_LE(condensation_metric(i), kTol) << metric_label(i);
    }
}

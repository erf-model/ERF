#include <gtest/gtest.h>

#include "ERF_GTestCoalescenceCommon.H"

using namespace coalescence_test;

// Motivation: ERF's super-droplet collisions use the kernels in
// Source/Particles/ERF_SuperDropletPCCoalescence.H (CollisionKernel). This checks
// the Golovin additive form, the symmetry / non-negativity / zero-at-no-shear /
// monotonicity of the sedimentation, Long, Hall and Brownian kernels, and the
// bounds and cutoffs of the Beard-Grover and Erfani-Mitchell riming efficiencies
// -- all on the host.
TEST(CoalescenceScalar, KernelPropertiesHoldOnHost)
{
    for (int i = 0; i < kNumMetrics; ++i) {
        EXPECT_LE(coalescence_metric(i), kTol) << metric_label(i);
    }
}

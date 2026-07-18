#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "ERF_NOAHMP_ResultPolicy.H"

using amrex::Real;

namespace
{

constexpr Real undefined_value = lsm_undefined;

} // namespace

// Motivation: Noah-MP can return finite-looking defaults for an unprocessed
// cell, so the HFX gate must invalidate every registered provider result.
TEST(NoahMPResultPolicy, UnprocessedCellsRemainUndefined)
{
    EXPECT_FALSE(noahmp_result_policy::result_is_valid(Real(-9999.0)));
    EXPECT_FALSE(noahmp_result_policy::result_is_valid(
        std::numeric_limits<Real>::infinity()));

#define EXPECT_UNPROCESSED(alias, lsm, out)                                  \
    EXPECT_EQ(noahmp_result_policy::select_result(                            \
        false, Real(42.0), undefined_value, NoahmpOutputComp::out),            \
        undefined_value);
    NOAHMP_RESULT_FIELDS(EXPECT_UNPROCESSED)
#undef EXPECT_UNPROCESSED

    for (const int output_comp : {NoahmpOutputComp::hfx,
                                  NoahmpOutputComp::lh,
                                  NoahmpOutputComp::tau_ew,
                                  NoahmpOutputComp::tau_ns}) {
        EXPECT_EQ(noahmp_result_policy::select_result(
                      false, Real(42.0), undefined_value, output_comp),
                  undefined_value);
    }

    EXPECT_EQ(noahmp_result_policy::select_result(
                  false, Real(42.0), undefined_value, -1), undefined_value);
}

// Motivation: HFX equal to zero is a valid processed result, and processed
// fields retain valid zero and signed physical values.
TEST(NoahMPResultPolicy, ProcessedValuesPreservePhysicalSignsAndZeros)
{
    EXPECT_TRUE(noahmp_result_policy::result_is_valid(Real(0.0)));
    EXPECT_EQ(noahmp_result_policy::select_result(
                  true, Real(0.0), undefined_value, NoahmpOutputComp::hfx),
              Real(0.0));
    EXPECT_EQ(noahmp_result_policy::select_result(
                  true, Real(-12.5), undefined_value, NoahmpOutputComp::o_grdflx),
              Real(-12.5));
    EXPECT_EQ(noahmp_result_policy::select_result(
                  true, Real(0.25), undefined_value, NoahmpOutputComp::tsk),
              Real(0.25));
}

// Motivation: Native humidity and vegetation outputs have physical domain
// constraints beyond the generic Noah-MP fill-value check.
TEST(NoahMPResultPolicy, NativeHumidityAndVegetationUseFieldRules)
{
    EXPECT_EQ(noahmp_result_policy::select_result(
                  true, Real(-0.001), undefined_value,
                  NoahmpOutputComp::o_noahmp_q2m_veg), undefined_value);
    EXPECT_EQ(noahmp_result_policy::select_result(
                  true, Real(-0.001), undefined_value,
                  NoahmpOutputComp::o_noahmp_q2m_bare), undefined_value);

    for (const Real value : {Real(0.0), Real(1.0)}) {
        EXPECT_EQ(noahmp_result_policy::select_result(
                      true, value, undefined_value, NoahmpOutputComp::o_noahmp_fveg),
                  value);
    }
    for (const Real value : {Real(-0.001), Real(1.001)}) {
        EXPECT_EQ(noahmp_result_policy::select_result(
                      true, value, undefined_value, NoahmpOutputComp::o_noahmp_fveg),
                  undefined_value);
    }
}

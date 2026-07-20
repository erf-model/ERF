#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "ERF_NOAHMP_ResultPolicy.H"

using amrex::Real;

namespace
{

constexpr Real undefined_value = lsm_undefined;

Real valid_registered_value (int output_comp)
{
    if (output_comp == NoahmpOutputComp::o_noahmp_t2m_veg ||
        output_comp == NoahmpOutputComp::o_noahmp_t2m_bare ||
        output_comp == NoahmpOutputComp::tsk) {
        return Real(290.0);
    }
    if (output_comp == NoahmpOutputComp::o_noahmp_q2m_veg ||
        output_comp == NoahmpOutputComp::o_noahmp_q2m_bare) {
        return Real(0.005);
    }
    if (output_comp == NoahmpOutputComp::o_noahmp_fveg) {
        return Real(0.5);
    }
    return Real(0.3);
}

void expect_registered_results (const noahmp_result_policy::CellPolicy& policy,
                                bool expect_valid)
{
#define CHECK_RESULT(alias, lsm, out)                                      \
    EXPECT_EQ(policy.select(valid_registered_value(NoahmpOutputComp::out),   \
                            undefined_value, NoahmpOutputComp::out),        \
              expect_valid ? valid_registered_value(NoahmpOutputComp::out)   \
                           : undefined_value);
    NOAHMP_RESULT_FIELDS(CHECK_RESULT)
#undef CHECK_RESULT
}

} // namespace

// Motivation: the production HFX-derived policy must invalidate every
// plausible registered result for each provider sentinel/nonfinite gate.
TEST(NoahMPResultPolicy, InvalidHfxBlocksPlausibleRegisteredResults)
{
    for (const Real hfx : {Real(-9999.0), undefined_value,
                           std::numeric_limits<Real>::infinity(),
                           std::numeric_limits<Real>::quiet_NaN()}) {
        const noahmp_result_policy::CellPolicy policy(hfx);
        EXPECT_FALSE(policy.processed);
        expect_registered_results(policy, false);
    }
}

// Motivation: zero HFX is a valid processed result, and ordinary finite HFX
// values preserve valid field-specific representatives.
TEST(NoahMPResultPolicy, ValidHfxPreservesPlausibleRegisteredResults)
{
    for (const Real hfx : {Real(0.0), Real(25.0)}) {
        const noahmp_result_policy::CellPolicy policy(hfx);
        EXPECT_TRUE(policy.processed);
        expect_registered_results(policy, true);
    }
}

// Motivation: processed-cell field rules remain active after the HFX gate.
TEST(NoahMPResultPolicy, ProcessedFieldRulesRemainActive)
{
    const noahmp_result_policy::CellPolicy policy(Real(0.0));

    EXPECT_EQ(policy.select(Real(-0.001), undefined_value,
                            NoahmpOutputComp::o_noahmp_q2m_veg), undefined_value);
    EXPECT_EQ(policy.select(Real(-0.001), undefined_value,
                            NoahmpOutputComp::o_noahmp_q2m_bare), undefined_value);
    EXPECT_EQ(policy.select(Real(-0.001), undefined_value,
                            NoahmpOutputComp::o_noahmp_fveg), undefined_value);
    EXPECT_EQ(policy.select(Real(1.001), undefined_value,
                            NoahmpOutputComp::o_noahmp_fveg), undefined_value);
    EXPECT_EQ(policy.select(Real(0.0), undefined_value,
                            NoahmpOutputComp::o_noahmp_fveg), Real(0.0));
    EXPECT_EQ(policy.select(Real(1.0), undefined_value,
                            NoahmpOutputComp::o_noahmp_fveg), Real(1.0));
    EXPECT_EQ(policy.select(Real(-9999.0), undefined_value,
                            NoahmpOutputComp::o_grdflx), undefined_value);
    EXPECT_EQ(policy.select(Real(-12.5), undefined_value,
                            NoahmpOutputComp::o_grdflx), Real(-12.5));
}

// Motivation: generic soil outputs use the HFX gate plus the generic result
// validity rule, without adding soil-specific physical restrictions.
TEST(NoahMPResultPolicy, SoilSelectionUsesHfxAndGenericRules)
{
    const noahmp_result_policy::CellPolicy invalid_hfx(Real(-9999.0));
    EXPECT_EQ(invalid_hfx.select(Real(0.3), undefined_value, -1), undefined_value);

    const noahmp_result_policy::CellPolicy valid_hfx(Real(0.0));
    EXPECT_EQ(valid_hfx.select(Real(0.3), undefined_value, -1), Real(0.3));
    EXPECT_EQ(valid_hfx.select(Real(-9999.0), undefined_value, -1), undefined_value);
    EXPECT_EQ(valid_hfx.select(undefined_value, undefined_value, -1), undefined_value);
}

// Motivation: heat, moisture, and signed stress fluxes share the production
// HFX gate while retaining the existing kinematic conversion formulas.
TEST(NoahMPResultPolicy, FluxSelectionUsesHfxGate)
{
    constexpr Real density = Real(2.0);
    constexpr Real latent_heat_flux = Real(100.0);
    constexpr Real tau_ew = Real(-4.0);
    constexpr Real tau_ns = Real(6.0);

    const auto invalid = noahmp_result_policy::CellPolicy(Real(-9999.0)).select_fluxes(
        latent_heat_flux, tau_ew, tau_ns, density, Cp_d, L_v, undefined_value);
    EXPECT_EQ(invalid.temperature, undefined_value);
    EXPECT_EQ(invalid.mixing_ratio, undefined_value);
    EXPECT_EQ(invalid.tau_ew, undefined_value);
    EXPECT_EQ(invalid.tau_ns, undefined_value);

    const auto zero_hfx = noahmp_result_policy::CellPolicy(Real(0.0)).select_fluxes(
        latent_heat_flux, tau_ew, tau_ns, density, Cp_d, L_v, undefined_value);
    EXPECT_EQ(zero_hfx.temperature, Real(0.0));
    EXPECT_EQ(zero_hfx.mixing_ratio, latent_heat_flux / (density * L_v));
    EXPECT_EQ(zero_hfx.tau_ew, tau_ew / density);
    EXPECT_EQ(zero_hfx.tau_ns, tau_ns / density);

    const auto nonzero_hfx = noahmp_result_policy::CellPolicy(Real(10.0)).select_fluxes(
        latent_heat_flux, tau_ew, tau_ns, density, Cp_d, L_v, undefined_value);
    EXPECT_EQ(nonzero_hfx.temperature, Real(10.0) / (density * Cp_d));
    EXPECT_EQ(nonzero_hfx.mixing_ratio, latent_heat_flux / (density * L_v));
    EXPECT_EQ(nonzero_hfx.tau_ew, tau_ew / density);
    EXPECT_EQ(nonzero_hfx.tau_ns, tau_ns / density);
}

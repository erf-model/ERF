#include <algorithm>
#include <array>
#include <string>

#include <gtest/gtest.h>

#include <ERF_KesslerUtils.H>

#include "ERF_GTestSatAdjCommon.H"

using namespace satadj_test;

namespace {

struct KesslerNoRainCellResult {
    CellState state;
    amrex::Real phase_change{amrex::Real(0.0)};
};

#ifdef AMREX_USE_FLOAT
constexpr amrex::Real kAsymptoticAbsMixTol = amrex::Real(2.0e-7);
constexpr amrex::Real kAsymptoticAbsThermoTol = amrex::Real(2.0e-4);
constexpr amrex::Real kFiniteGapMixMin = amrex::Real(2.0e-5);
constexpr amrex::Real kFiniteGapThermoMin = amrex::Real(2.0e-3);
#else
constexpr amrex::Real kAsymptoticAbsMixTol = amrex::Real(1.0e-13);
constexpr amrex::Real kAsymptoticAbsThermoTol = amrex::Real(1.0e-10);
constexpr amrex::Real kFiniteGapMixMin = amrex::Real(1.0e-6);
constexpr amrex::Real kFiniteGapThermoMin = amrex::Real(1.0e-4);
#endif

amrex::Real asymptotic_mix_tol (const amrex::Real phase_change_scale)
{
    return std::max(kAsymptoticAbsMixTol,
                    kAsymptoticRelTol * amrex::Math::abs(phase_change_scale));
}

amrex::Real asymptotic_thermo_tol (const CellState& initial,
                                   const amrex::Real phase_change_scale)
{
    const amrex::Real theta_over_t = initial.theta / initial.tabs;
    const amrex::Real thermo_scale = theta_over_t * kFacCond * amrex::Math::abs(phase_change_scale);
    return std::max(kAsymptoticAbsThermoTol,
                    kAsymptoticRelTol * thermo_scale);
}

KesslerNoRainCellResult apply_kessler_norain_scalar (const CellState& initial)
{
    CellState state = initial;
    state.qc = std::max(amrex::Real(0.0), state.qc);

    const amrex::Real qsat_initial = qsat(initial.tabs, initial.pres_mbar);
    const amrex::Real dqsat_initial = dqsat_dt(initial.tabs, initial.pres_mbar);
    const KesslerSaturationAdjustment saturation_adjustment =
        kessler_saturation_adjustment(state.qv,
                                      state.qc,
                                      qsat_initial,
                                      dqsat_initial,
                                      true,
                                      kFacCond);

    const amrex::Real phase_change =
        saturation_adjustment.dq_vapor_to_cloud - saturation_adjustment.dq_cloud_to_vapor;

    state.qv += -saturation_adjustment.dq_vapor_to_cloud
              +  saturation_adjustment.dq_cloud_to_vapor;
    state.qc += phase_change;

    const amrex::Real theta_over_t = state.theta / state.tabs;
    state.theta += theta_over_t * kFacCond * phase_change;

    state.qv = std::max(amrex::Real(0.0), state.qv);
    state.qc = std::max(amrex::Real(0.0), state.qc);
    state.tabs = getTgivenPandTh(amrex::Real(100.0) * initial.pres_mbar,
                                 state.theta,
                                 kRdOcp);

    return KesslerNoRainCellResult{state, phase_change};
}

void expect_asymptotically_close (const CellState& initial,
                                  const CellState& satadj_state,
                                  const KesslerNoRainCellResult& kessler_result)
{
    const amrex::Real phase_change_scale = std::max(
        amrex::Math::abs(satadj_state.qc - initial.qc),
        amrex::Math::abs(kessler_result.phase_change));

    EXPECT_NEAR(kessler_result.state.qv,
                satadj_state.qv,
                asymptotic_mix_tol(phase_change_scale));
    EXPECT_NEAR(kessler_result.state.qc,
                satadj_state.qc,
                asymptotic_mix_tol(phase_change_scale));
    EXPECT_NEAR(kessler_result.state.theta,
                satadj_state.theta,
                asymptotic_thermo_tol(initial, phase_change_scale));
    EXPECT_NEAR(kessler_result.state.tabs,
                satadj_state.tabs,
                asymptotic_thermo_tol(initial, phase_change_scale));
}

} // namespace

// Motivation: Kessler_NoRain linearizes saturation adjustment about the
// initial state, while SatAdj solves the nonlinear fixed-pressure adjustment.
// For sufficiently small supersaturation, the two should agree asymptotically
// even though they are not an exact equality contract at finite amplitude.
TEST(SatAdjKesslerComparison, SmallCondensationAsymptoticallyConsistent)
{
    const amrex::Real tabs = amrex::Real(290.0);
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qsat_initial = qsat(tabs, pres_mbar);
    // Float builds need perturbations above the absolute tolerance floor for this
    // asymptotic comparison to remain diagnostic. Double builds keep smaller
    // perturbations to exercise the limiting behavior closer to roundoff.
#ifdef AMREX_USE_FLOAT
    const std::array<amrex::Real, 3> supersaturation = {
        amrex::Real(1.0e-6), amrex::Real(3.0e-6), amrex::Real(1.0e-5)};
#else
    const std::array<amrex::Real, 3> supersaturation = {
        amrex::Real(1.0e-8), amrex::Real(1.0e-7), amrex::Real(1.0e-6)};
#endif

    for (const amrex::Real delta_qv : supersaturation) {
        SCOPED_TRACE("delta_qv=" + std::to_string(static_cast<double>(delta_qv)));
        const CellState initial = make_cell_state(tabs,
                                                  pres_mbar,
                                                  qsat_initial + delta_qv,
                                                  amrex::Real(0.0));
        CellState satadj_state = initial;
        adjust(satadj_state);
        const KesslerNoRainCellResult kessler_result = apply_kessler_norain_scalar(initial);

        EXPECT_GT(satadj_state.qc, amrex::Real(0.0));
        EXPECT_GT(kessler_result.phase_change, amrex::Real(0.0));
        expect_asymptotically_close(initial, satadj_state, kessler_result);
    }
}

// Motivation: The same asymptotic agreement should hold on the cloudy
// subsaturated side when only a small amount of cloud water evaporates.
TEST(SatAdjKesslerComparison, SmallCloudEvaporationAsymptoticallyConsistent)
{
    const amrex::Real tabs = amrex::Real(290.0);
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qsat_initial = qsat(tabs, pres_mbar);
    // Float builds need perturbations above the absolute tolerance floor for this
    // asymptotic comparison to remain diagnostic. Double builds keep smaller
    // perturbations to exercise the limiting behavior closer to roundoff.
#ifdef AMREX_USE_FLOAT
    const std::array<amrex::Real, 3> subsaturation = {
        amrex::Real(1.0e-6), amrex::Real(3.0e-6), amrex::Real(1.0e-5)};
#else
    const std::array<amrex::Real, 3> subsaturation = {
        amrex::Real(1.0e-8), amrex::Real(1.0e-7), amrex::Real(1.0e-6)};
#endif

    for (const amrex::Real deficit_qv : subsaturation) {
        SCOPED_TRACE("deficit_qv=" + std::to_string(static_cast<double>(deficit_qv)));
        const CellState initial = make_cell_state(tabs,
                                                  pres_mbar,
                                                  qsat_initial - deficit_qv,
                                                  amrex::Real(5.0e-4));
        CellState satadj_state = initial;
        adjust(satadj_state);
        const KesslerNoRainCellResult kessler_result = apply_kessler_norain_scalar(initial);

        EXPECT_LT(kessler_result.phase_change, amrex::Real(0.0));
        EXPECT_LT(satadj_state.qc, initial.qc);
        expect_asymptotically_close(initial, satadj_state, kessler_result);
    }
}

// Motivation: Kessler_NoRain is only a first-order approximation to SatAdj.
// At finite amplitude the two updates should conserve total nonprecipitating
// water individually, but should show a measurable state gap so future tests do
// not accidentally tighten this into a false exact-equality contract.
TEST(SatAdjKesslerComparison, FiniteAdjustmentDifferencesAreDocumented)
{
    struct Case {
        std::string name;
        CellState initial;
    };

    const amrex::Real tabs = amrex::Real(290.0);
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qsat_initial = qsat(tabs, pres_mbar);
    CellState evap_then_recond;
    ASSERT_TRUE(find_evaporation_then_recondensation_state(evap_then_recond));
    const std::array<Case, 2> cases = {{
        {"finite-condensation",
         make_cell_state(tabs,
                         pres_mbar,
                         qsat_initial + amrex::Real(8.0e-4),
                         amrex::Real(4.0e-4))},
        {"finite-evaporation-then-recondensation",
         evap_then_recond}
    }};

    for (const Case& test_case : cases) {
        SCOPED_TRACE(test_case.name);
        CellState satadj_state = test_case.initial;
        adjust(satadj_state);
        const KesslerNoRainCellResult kessler_result = apply_kessler_norain_scalar(test_case.initial);

        const amrex::Real qv_gap = amrex::Math::abs(satadj_state.qv - kessler_result.state.qv);
        const amrex::Real qc_gap = amrex::Math::abs(satadj_state.qc - kessler_result.state.qc);
        const amrex::Real theta_gap = amrex::Math::abs(satadj_state.theta - kessler_result.state.theta);
        const amrex::Real tabs_gap = amrex::Math::abs(satadj_state.tabs - kessler_result.state.tabs);

        const amrex::Real qt_initial = test_case.initial.qv + test_case.initial.qc;
        const amrex::Real qt_satadj = satadj_state.qv + satadj_state.qc;
        const amrex::Real qt_kessler = kessler_result.state.qv + kessler_result.state.qc;

        EXPECT_NEAR(qt_satadj, qt_initial, asymptotic_mix_tol(qt_initial))
            << "case=" << test_case.name
            << " qt_initial=" << qt_initial
            << " qt_satadj=" << qt_satadj;

        EXPECT_NEAR(qt_kessler, qt_initial, asymptotic_mix_tol(qt_initial))
            << "case=" << test_case.name
            << " qt_initial=" << qt_initial
            << " qt_kessler=" << qt_kessler;

        EXPECT_TRUE(std::max(qv_gap, qc_gap) > kFiniteGapMixMin
                    || std::max(theta_gap, tabs_gap) > kFiniteGapThermoMin)
            << "qv_gap=" << qv_gap
            << " qc_gap=" << qc_gap
            << " theta_gap=" << theta_gap
            << " tabs_gap=" << tabs_gap;
    }
}

// Motivation: The asymptotic tests show that Kessler_NoRain approaches SatAdj
// for small adjustments, while the finite-amplitude test documents
// non-equivalence. This characterization bounds the finite gaps relative to the
// phase-change scale so the difference can be interpreted physically without
// turning the schemes into a false equality contract.
TEST(SatAdjKesslerComparison, FiniteAdjustmentGapMagnitudesAreCharacterized)
{
    struct Case {
        std::string name;
        CellState initial;
    };

    const amrex::Real tabs = amrex::Real(290.0);
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qsat_initial = qsat(tabs, pres_mbar);
    CellState evap_then_recond;
    ASSERT_TRUE(find_evaporation_then_recondensation_state(evap_then_recond));

    const std::array<Case, 3> cases = {{
        {"finite-condensation",
         make_cell_state(tabs, pres_mbar, qsat_initial + amrex::Real(8.0e-4), amrex::Real(4.0e-4))},
        {"finite-cloudy-evaporation",
         make_cell_state(tabs, pres_mbar, qsat_initial - amrex::Real(5.0e-4), amrex::Real(1.0e-3))},
        {"finite-evaporation-then-recondensation",
         evap_then_recond}
    }};

    for (const auto& c : cases) {
        SCOPED_TRACE(c.name);
        CellState satadj_state = c.initial;
        adjust(satadj_state);
        const KesslerNoRainCellResult kessler_result = apply_kessler_norain_scalar(c.initial);

        const amrex::Real qv_gap = amrex::Math::abs(satadj_state.qv - kessler_result.state.qv);
        const amrex::Real qc_gap = amrex::Math::abs(satadj_state.qc - kessler_result.state.qc);
        const amrex::Real theta_gap = amrex::Math::abs(satadj_state.theta - kessler_result.state.theta);
        const amrex::Real tabs_gap = amrex::Math::abs(satadj_state.tabs - kessler_result.state.tabs);

        const amrex::Real phase_change_scale = std::max(
            amrex::Math::abs(satadj_state.qc - c.initial.qc),
            amrex::Math::abs(kessler_result.phase_change));

        // Guard: finite cases must produce a nonzero phase-change scale so the
        // normalized-gap divisions below are well-defined.
        ASSERT_GT(phase_change_scale, amrex::Real(0.0))
            << "case=" << c.name
            << " satadj_phase_change=" << satadj_state.qc - c.initial.qc
            << " kessler_phase_change=" << kessler_result.phase_change;

        // latent-temperature scale = (L/cp) * phase_change_scale
        const amrex::Real latent_temp_scale = kFacCond * phase_change_scale;

        // Normalized gaps
        const amrex::Real water_norm = phase_change_scale > amrex::Real(0.0) ?
            std::max(qv_gap, qc_gap) / phase_change_scale : amrex::Real(0.0);
        const amrex::Real thermo_norm = latent_temp_scale > amrex::Real(0.0) ?
            std::max(theta_gap, tabs_gap) / latent_temp_scale : amrex::Real(0.0);

        SCOPED_TRACE("qv_gap=" + std::to_string(static_cast<double>(qv_gap))
                     + " qc_gap=" + std::to_string(static_cast<double>(qc_gap))
                     + " theta_gap=" + std::to_string(static_cast<double>(theta_gap))
                     + " tabs_gap=" + std::to_string(static_cast<double>(tabs_gap))
                     + " phase_change_scale=" + std::to_string(static_cast<double>(phase_change_scale))
                     + " water_norm=" + std::to_string(static_cast<double>(water_norm))
                     + " thermo_norm=" + std::to_string(static_cast<double>(thermo_norm)));

        EXPECT_TRUE(std::isfinite(qv_gap) && std::isfinite(qc_gap)
                    && std::isfinite(theta_gap) && std::isfinite(tabs_gap))
            << "case=" << c.name
            << " qv_gap=" << qv_gap << " qc_gap=" << qc_gap
            << " theta_gap=" << theta_gap << " tabs_gap=" << tabs_gap;
        EXPECT_TRUE(std::isfinite(water_norm) && std::isfinite(thermo_norm))
            << "case=" << c.name
            << " water_norm=" << water_norm << " thermo_norm=" << thermo_norm;

        // Loose sanity upper bounds to catch absurd regressions while avoiding
        // false failures in float/GPU builds.
        constexpr amrex::Real kMaxNormalizedWaterGap = amrex::Real(0.5);
        constexpr amrex::Real kMaxNormalizedThermoGap = amrex::Real(0.5);
        EXPECT_LT(qv_gap / phase_change_scale, kMaxNormalizedWaterGap)
            << "case=" << c.name
            << " qv_gap=" << qv_gap
            << " phase_change_scale=" << phase_change_scale
            << " normalized_qv_gap=" << qv_gap / phase_change_scale
            << " satadj_qv=" << satadj_state.qv
            << " kessler_qv=" << kessler_result.state.qv;
        EXPECT_LT(qc_gap / phase_change_scale, kMaxNormalizedWaterGap)
            << "case=" << c.name
            << " qc_gap=" << qc_gap
            << " phase_change_scale=" << phase_change_scale
            << " normalized_qc_gap=" << qc_gap / phase_change_scale
            << " satadj_qc=" << satadj_state.qc
            << " kessler_qc=" << kessler_result.state.qc;
        if (latent_temp_scale > amrex::Real(0.0)) {
            EXPECT_LT(theta_gap / latent_temp_scale, kMaxNormalizedThermoGap)
                << "case=" << c.name
                << " theta_gap=" << theta_gap
                << " latent_temp_scale=" << latent_temp_scale
                << " normalized_theta_gap=" << theta_gap / latent_temp_scale
                << " satadj_theta=" << satadj_state.theta
                << " kessler_theta=" << kessler_result.state.theta;
            EXPECT_LT(tabs_gap / latent_temp_scale, kMaxNormalizedThermoGap)
                << "case=" << c.name
                << " tabs_gap=" << tabs_gap
                << " latent_temp_scale=" << latent_temp_scale
                << " normalized_tabs_gap=" << tabs_gap / latent_temp_scale
                << " satadj_tabs=" << satadj_state.tabs
                << " kessler_tabs=" << kessler_result.state.tabs;
        }
    }
}
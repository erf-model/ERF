#include "ERF_Constants.H"
#include "ERF_EOS.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_ShocEnergyFixer.H"
#include "ERF_ShocImplicit.H"
#include "ERF_ShocPDF.H"
#include "ERF_ShocStructure.H"
#include "ERF_ShocTKE.H"
#include "ERF_ShocMoments.H"
#include "ERF_ShocTestUtils.H"

#include <gtest/gtest.h>

#include <cmath>

namespace
{
using namespace amrex::literals;

amrex::Real
column_moist_energy (const ShocColumnData& col)
{
    const auto rho = col.rho.const_array();
    const auto dz = col.dz.const_array();
    const auto tabs = col.tabs.const_array();
    const auto zt = col.zt.const_array();
    const auto qv = col.qv.const_array();
    const auto qc = col.qc.const_array();
    const auto qi = col.qi.const_array();
    const auto u = col.u.const_array();
    const auto v = col.v.const_array();
    const auto tke = col.tke.const_array();
    amrex::Real total = 0.0;
    for (int k = 0; k < col.layout.nlev; ++k) {
        const amrex::Real mass = rho(0,k,0) * dz(0,k,0);
        total += mass
               * (Cp_d * tabs(0,k,0)
                  + CONST_GRAV * zt(0,k,0)
                  + amrex::Real(0.5) * (u(0,k,0) * u(0,k,0) + v(0,k,0) * v(0,k,0))
                  + tke(0,k,0)
                  + L_v * (qv(0,k,0) + qc(0,k,0))
                  + (L_v + amrex::Real(3.34e5)) * qi(0,k,0));
    }
    return total;
}

amrex::Real
column_total_water (const ShocColumnData& col)
{
    const auto rho = col.rho.const_array();
    const auto dz = col.dz.const_array();
    const auto qv = col.qv.const_array();
    const auto qc = col.qc.const_array();
    const auto qi = col.qi.const_array();

    amrex::Real total = 0.0_rt;
    for (int k = 0; k < col.layout.nlev; ++k) {
        total += rho(0,k,0) * dz(0,k,0)
               * (qv(0,k,0) + qc(0,k,0) + qi(0,k,0));
    }
    return total;
}

amrex::Real
qsat_liquid_from_thetal_and_pressure (amrex::Real thetal,
                                      amrex::Real p_pa)
{
    const amrex::Real temp = getTgivenPandTh(p_pa, thetal, R_d / Cp_d);
    amrex::Real qs = 0.0_rt;
    erf_qsatw(temp, 0.01_rt * p_pa, qs);
    return qs;
}

amrex::Real
deterministic_pdf_condensate (amrex::Real qw,
                              amrex::Real qs,
                              amrex::Real temp)
{
    const amrex::Real beta =
        (R_d / R_v) * (L_v / (R_d * temp)) * (L_v / (Cp_d * temp));
    return amrex::max(0.0_rt, (qw - qs) / (1.0_rt + beta * qs));
}

void
zero_pdf_moments (ShocColumnData& col)
{
    auto w = col.w.array();
    auto w_sec = col.w_sec.array();
    auto shoc_cldfrac = col.shoc_cldfrac.array();
    auto shoc_ql2 = col.shoc_ql2.array();
    auto wqls_sec = col.wqls_sec.array();
    auto shoc_cond = col.shoc_cond.array();
    auto shoc_evap = col.shoc_evap.array();
    auto wthv_sec = col.wthv_sec.array();
    auto thl_sec = col.thl_sec.array();
    auto qw_sec = col.qw_sec.array();
    auto qwthl_sec = col.qwthl_sec.array();
    auto wthl_sec = col.wthl_sec.array();
    auto wqw_sec = col.wqw_sec.array();
    auto uw_sec = col.uw_sec.array();
    auto vw_sec = col.vw_sec.array();
    auto wtke_sec = col.wtke_sec.array();
    auto w3 = col.w3.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        w(0,k,0) = 0.0_rt;
        w_sec(0,k,0) = 0.0_rt;
        shoc_cldfrac(0,k,0) = 0.0_rt;
        shoc_ql2(0,k,0) = 0.0_rt;
        wqls_sec(0,k,0) = 0.0_rt;
        shoc_cond(0,k,0) = 0.0_rt;
        shoc_evap(0,k,0) = 0.0_rt;
        wthv_sec(0,k,0) = 0.0_rt;
    }

    for (int k = 0; k <= col.layout.nlev; ++k) {
        thl_sec(0,k,0) = 0.0_rt;
        qw_sec(0,k,0) = 0.0_rt;
        qwthl_sec(0,k,0) = 0.0_rt;
        wthl_sec(0,k,0) = 0.0_rt;
        wqw_sec(0,k,0) = 0.0_rt;
        uw_sec(0,k,0) = 0.0_rt;
        vw_sec(0,k,0) = 0.0_rt;
        wtke_sec(0,k,0) = 0.0_rt;
        w3(0,k,0) = 0.0_rt;
    }
}

void
set_uniform_implicit_column (ShocColumnData& col,
                             amrex::Real thetal_val,
                             amrex::Real exner_val,
                             amrex::Real qv_val,
                             amrex::Real qc_val,
                             amrex::Real qi_val,
                             amrex::Real tke_val,
                             amrex::Real u_val,
                             amrex::Real v_val)
{
    const amrex::Real rho_val = col.rho.const_array()(0,0,0);
    auto rho = col.rho.array();
    auto thetal = col.thetal.array();
    auto theta = col.theta.array();
    auto exner = col.exner.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();
    auto tabs = col.tabs.array();
    auto theta_v = col.theta_v.array();
    auto host_dse = col.host_dse.array();
    auto tk = col.tk.array();
    auto tkh = col.tkh.array();
    auto tke = col.tke.array();
    auto u = col.u.array();
    auto v = col.v.array();
    auto shoc_ql = col.shoc_ql.array();
    const auto zt = col.zt.const_array();

    const amrex::Real ql_total = qc_val + qi_val;
    const amrex::Real tabs_val = thetal_val * exner_val
                               + (L_v * qc_val + shoc::latent_sublimation() * qi_val) / Cp_d;
    const amrex::Real theta_val = tabs_val / amrex::max(exner_val, 1.0e-12_rt);
    const amrex::Real qw_val = qv_val + ql_total;

    for (int k = 0; k < col.layout.nlev; ++k) {
        rho(0,k,0) = rho_val;
        thetal(0,k,0) = thetal_val;
        theta(0,k,0) = theta_val;
        exner(0,k,0) = exner_val;
        qv(0,k,0) = qv_val;
        qc(0,k,0) = qc_val;
        qi(0,k,0) = qi_val;
        qw(0,k,0) = qw_val;
        tabs(0,k,0) = tabs_val;
        theta_v(0,k,0) = theta_val * (1.0_rt + 0.61_rt * qv_val - ql_total);
        host_dse(0,k,0) = Cp_d * tabs_val + CONST_GRAV * zt(0,k,0);
        tk(0,k,0) = 0.0_rt;
        tkh(0,k,0) = 0.0_rt;
        tke(0,k,0) = tke_val;
        u(0,k,0) = u_val;
        v(0,k,0) = v_val;
        shoc_ql(0,k,0) = ql_total;
    }
}
}

TEST(ShocPhysical, ColumnHeatBudgetTracksSurfaceFlux)
{
    auto col = shoc_test::make_column(6);
    ShocRuntimeOptions opts;

    auto tk = col.tk.array();
    auto tkh = col.tkh.array();
    auto exner = col.exner.array();
    auto thetal = col.thetal.array();
    auto theta = col.theta.array();
    auto tabs = col.tabs.array();
    auto host_dse = col.host_dse.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();
    for (int k = 0; k < col.layout.nlev; ++k) {
        tk(0,k,0) = 0.0;
        tkh(0,k,0) = 0.0;
        exner(0,k,0) = 1.0;
        thetal(0,k,0) = 300.0;
        theta(0,k,0) = 300.0;
        tabs(0,k,0) = 300.0;
        qv(0,k,0) = 0.008;
        qc(0,k,0) = 0.0;
        qi(0,k,0) = 0.0;
        qw(0,k,0) = qv(0,k,0);
        host_dse(0,k,0) = Cp_d * tabs(0,k,0) + CONST_GRAV * col.zt.const_array()(0,k,0);
        col.tke_tend.array()(0,k,0) = 0.0;
    }
    shoc::set_fab_val(col.surf_sens_flux, 0.02, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.0, shoc::InitRunOn::Host);

    const amrex::Real before = column_moist_energy(col);
    const amrex::Real dt = 10.0;
    const amrex::Real rho_sfc = col.rho.const_array()(0,0,0);

    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(col, opts, dt);
    });

    const amrex::Real after = column_moist_energy(col);
    const amrex::Real expected = dt * rho_sfc * Cp_d * 0.02;
    EXPECT_NEAR(after - before, expected, 5.0e-9 * amrex::max(amrex::Real(1.0), amrex::Math::abs(expected)));
}

TEST(ShocPhysical, StrongerSurfaceHeatingRaisesMeanThetaMore)
{
    auto weak = shoc_test::make_column(6);
    auto strong = shoc_test::make_column(6);
    ShocRuntimeOptions opts;

    for (auto* col : {&weak, &strong}) {
        auto tk = col->tk.array();
        auto tkh = col->tkh.array();
        auto exner = col->exner.array();
        for (int k = 0; k < col->layout.nlev; ++k) {
            tk(0,k,0) = 1.0;
            tkh(0,k,0) = 1.0;
            exner(0,k,0) = 1.0;
            col->tke_tend.array()(0,k,0) = 0.0;
        }
    }

    shoc::set_fab_val(weak.surf_sens_flux, 0.01, shoc::InitRunOn::Host);
    shoc::set_fab_val(strong.surf_sens_flux, 0.03, shoc::InitRunOn::Host);

    const amrex::Real weak_before = weak.thetal.const_array()(0,0,0);
    const amrex::Real strong_before = strong.thetal.const_array()(0,0,0);

    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(weak, opts, 10.0);
        ShocImplicit::update_prognostics(strong, opts, 10.0);
    });

    const amrex::Real weak_delta = weak.thetal.const_array()(0,0,0) - weak_before;
    const amrex::Real strong_delta = strong.thetal.const_array()(0,0,0) - strong_before;
    EXPECT_GT(strong_delta, weak_delta);
}

TEST(ShocPhysical, UnstableSurfaceForcingDeepensPblMoreThanStableForcing)
{
    auto stable = shoc_test::make_column(8);
    auto unstable = shoc_test::make_column(8);

    shoc::set_fab_val(stable.surf_sens_flux, -0.01, shoc::InitRunOn::Host);
    shoc::set_fab_val(unstable.surf_sens_flux, 0.02, shoc::InitRunOn::Host);
    shoc::set_fab_val(stable.surf_lat_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(unstable.surf_lat_flux, 0.0, shoc::InitRunOn::Host);

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(stable);
        ShocStructure::diagnose_pblh(stable);
        ShocStructure::diagnose_surface_layer(unstable);
        ShocStructure::diagnose_pblh(unstable);
    });

    EXPECT_GE(unstable.pblh.const_array()(0,0,0), stable.pblh.const_array()(0,0,0));
}

TEST(ShocPhysical, TwoSmallStepsTrackOneLargeStep)
{
    auto one_step = shoc_test::make_column(5);
    auto two_step = shoc_test::make_column(5);
    ShocRuntimeOptions opts;

    for (auto* col : {&one_step, &two_step}) {
        auto tk = col->tk.array();
        auto tkh = col->tkh.array();
        auto exner = col->exner.array();
        for (int k = 0; k < col->layout.nlev; ++k) {
            tk(0,k,0) = 1.0;
            tkh(0,k,0) = 1.0;
            exner(0,k,0) = 1.0;
            col->tke_tend.array()(0,k,0) = 0.0;
        }
        shoc::set_fab_val(col->surf_sens_flux, 0.02, shoc::InitRunOn::Host);
        shoc::set_fab_val(col->surf_lat_flux, 5.0e-5, shoc::InitRunOn::Host);
    }

    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(one_step, opts, 10.0);
        ShocImplicit::update_prognostics(two_step, opts, 5.0);
        ShocImplicit::update_prognostics(two_step, opts, 5.0);
    });

    const auto th_one = one_step.thetal.const_array();
    const auto th_two = two_step.thetal.const_array();
    const auto qw_one = one_step.qw.const_array();
    const auto qw_two = two_step.qw.const_array();

    for (int k = 0; k < one_step.layout.nlev; ++k) {
        EXPECT_NEAR(th_one(0,k,0), th_two(0,k,0), 5.0e-3);
        EXPECT_NEAR(qw_one(0,k,0), qw_two(0,k,0), 5.0e-5);
    }
}

TEST(ShocPhysical, SubfreezingCondensatePreservesIceSeed)
{
    auto col = shoc_test::make_column(5);
    ShocRuntimeOptions opts;

    auto exner = col.exner.array();
    auto thetal = col.thetal.array();
    auto qw = col.qw.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto tabs = col.tabs.array();
    auto tk = col.tk.array();
    auto tkh = col.tkh.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        exner(0,k,0) = 1.0;
        thetal(0,k,0) = 265.0;
        tabs(0,k,0) = 265.0;
        qw(0,k,0) = 0.004;
        qv(0,k,0) = 0.003;
        qc(0,k,0) = 0.0;
        qi(0,k,0) = 1.0e-3;
        tk(0,k,0) = 0.1;
        tkh(0,k,0) = 0.1;
        col.tke_tend.array()(0,k,0) = 0.0;
    }

    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(col, opts, 10.0);
    });

    const auto qc_new = col.qc.const_array();
    const auto qi_new = col.qi.const_array();
    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_NEAR(qi_new(0,k,0), 1.0e-3_rt, 1.0e-14_rt);
        EXPECT_GE(qc_new(0,k,0), 0.0_rt);
    }
}

TEST(ShocPhysical, TerrainLikeVerticalGridStaysFinite)
{
    auto col = shoc_test::make_column(6);
    ShocRuntimeOptions opts;
    auto zi = col.zi.array();
    auto zt = col.zt.array();
    auto dz = col.dz.array();

    for (int k = 0; k <= col.layout.nlev; ++k) {
        zi(0,k,0) = 15.0 * k + 2.5 * k * k;
        if (k < col.layout.nlev) {
            const amrex::Real zlo = zi(0,k,0);
            const amrex::Real zhi = zi(0,k+1,0);
            zt(0,k,0) = 0.5 * (zlo + zhi);
            dz(0,k,0) = zhi - zlo;
        }
    }

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
        ShocStructure::diagnose_pblh(col);
        ShocStructure::diagnose_length_and_brunt(col, opts, 300.0, 500.0);
        ShocTKE::diagnose_tke_and_diffusivities(col, opts, 300.0);
        ShocMoments::diagnose_moments(col, opts);
        ShocPDF::diagnose_pdf(col, opts, 10.0);
    });

    const auto mix = col.shoc_mix.const_array();
    const auto brunt = col.brunt.const_array();
    const auto cldfrac = col.shoc_cldfrac.const_array();
    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_TRUE(std::isfinite(mix(0,k,0)));
        EXPECT_TRUE(std::isfinite(brunt(0,k,0)));
        EXPECT_TRUE(std::isfinite(cldfrac(0,k,0)));
    }
}

TEST(ShocPhysical, SupersaturationCondensesInDeterministicPdfLimit)
{
    // In the zero-variance PDF limit, SHOC's saturation variable reduces to
    // s = (qw - qsat) / (1 + beta*qsat). Supersaturation should therefore
    // produce cloud liquid, full cloud fraction, positive condensation, and no evaporation.
    auto col = shoc_test::make_column(3);
    ShocRuntimeOptions opts;
    opts.extra_shoc_diags = true;
    const amrex::Real dt = 10.0_rt;
    const amrex::Real thetal_val = 300.0_rt;
    const amrex::Real p_mid_val = 90000.0_rt;
    const amrex::Real qs = qsat_liquid_from_thetal_and_pressure(thetal_val, p_mid_val);
    const amrex::Real qw_val = 1.25_rt * qs;
    const amrex::Real temp = getTgivenPandTh(p_mid_val, thetal_val, R_d / Cp_d);
    const amrex::Real expected_ql = deterministic_pdf_condensate(qw_val, qs, temp);

    ASSERT_GT(qs, 0.0_rt);
    ASSERT_GT(qw_val, 0.0_rt);
    ASSERT_GT(expected_ql, 0.0_rt);
    ASSERT_LT(expected_ql, qw_val);

    auto thetal = col.thetal.array();
    auto p_mid = col.p_mid.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();
    auto shoc_ql = col.shoc_ql.array();

    zero_pdf_moments(col);
    for (int k = 0; k < col.layout.nlev; ++k) {
        thetal(0,k,0) = thetal_val;
        p_mid(0,k,0) = p_mid_val;
        qv(0,k,0) = qw_val;
        qc(0,k,0) = 0.0_rt;
        qi(0,k,0) = 0.0_rt;
        qw(0,k,0) = qw_val;
        shoc_ql(0,k,0) = 0.0_rt;
    }

    shoc_test::run_and_sync([&] {
        ShocPDF::diagnose_pdf(col, opts, dt);
    });

    const auto shoc_ql_out = col.shoc_ql.const_array();
    const auto shoc_cldfrac = col.shoc_cldfrac.const_array();
    const auto shoc_cond = col.shoc_cond.const_array();
    const auto shoc_evap = col.shoc_evap.const_array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_GT(shoc_ql_out(0,k,0), 0.0_rt);
        EXPECT_NEAR(shoc_cldfrac(0,k,0), 1.0_rt, 1.0e-12_rt);
        EXPECT_GT(shoc_cond(0,k,0), 0.0_rt);
        EXPECT_NEAR(shoc_evap(0,k,0), 0.0_rt, 1.0e-14_rt);
        EXPECT_NEAR(shoc_ql_out(0,k,0), expected_ql,
                    1.0e-10_rt * amrex::max(1.0_rt, expected_ql));
    }
}

TEST(ShocPhysical, SubsaturationEvaporatesExistingCondensateInDeterministicPdfLimit)
{
    // In the zero-variance PDF limit, subsaturation gives s < 0, so the new
    // diagnosed PDF condensate should be zero. If old condensate was present,
    // extra SHOC diagnostics should report evaporation.
    auto col = shoc_test::make_column(3);
    ShocRuntimeOptions opts;
    opts.extra_shoc_diags = true;
    const amrex::Real dt = 10.0_rt;
    const amrex::Real thetal_val = 300.0_rt;
    const amrex::Real p_mid_val = 90000.0_rt;
    const amrex::Real qs = qsat_liquid_from_thetal_and_pressure(thetal_val, p_mid_val);
    const amrex::Real old_ql = 1.0e-4_rt;
    const amrex::Real qw_val = 0.75_rt * qs;
    const amrex::Real qv_val = qw_val - old_ql;

    ASSERT_GT(qs, 0.0_rt);
    ASSERT_GT(qw_val, 0.0_rt);
    ASSERT_GT(old_ql, 0.0_rt);
    ASSERT_GT(qv_val, 0.0_rt);
    ASSERT_EQ(deterministic_pdf_condensate(qw_val, qs,
                                           getTgivenPandTh(p_mid_val, thetal_val, R_d / Cp_d)),
              0.0_rt);

    auto thetal = col.thetal.array();
    auto p_mid = col.p_mid.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();
    auto shoc_ql = col.shoc_ql.array();

    zero_pdf_moments(col);
    for (int k = 0; k < col.layout.nlev; ++k) {
        thetal(0,k,0) = thetal_val;
        p_mid(0,k,0) = p_mid_val;
        qv(0,k,0) = qv_val;
        qc(0,k,0) = old_ql;
        qi(0,k,0) = 0.0_rt;
        qw(0,k,0) = qw_val;
        shoc_ql(0,k,0) = old_ql;
    }

    shoc_test::run_and_sync([&] {
        ShocPDF::diagnose_pdf(col, opts, dt);
    });

    const auto shoc_ql_out = col.shoc_ql.const_array();
    const auto shoc_cldfrac = col.shoc_cldfrac.const_array();
    const auto shoc_cond = col.shoc_cond.const_array();
    const auto shoc_evap = col.shoc_evap.const_array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_NEAR(shoc_ql_out(0,k,0), 0.0_rt, 1.0e-14_rt);
        EXPECT_NEAR(shoc_cldfrac(0,k,0), 0.0_rt, 1.0e-14_rt);
        EXPECT_NEAR(shoc_cond(0,k,0), 0.0_rt, 1.0e-14_rt);
        EXPECT_NEAR(shoc_evap(0,k,0), old_ql / dt, 1.0e-12_rt);
    }
}

TEST(ShocPhysical, PdfStateReconstructionConservesWaterAndLiquidWaterPotentialTemperature)
{
    const amrex::Real p = 90000.0_rt;
    const amrex::Real exner = getExnergivenP(p, R_d / Cp_d);

    {
        const amrex::Real thetal = 300.0_rt;
        const amrex::Real qw = 0.012_rt;
        const amrex::Real pdf_ql = 0.002_rt;
        const amrex::Real qi_seed = 0.0_rt;
        const amrex::Real unclamped_tabs = thetal * exner + (L_v / Cp_d) * pdf_ql;
        amrex::Real tabs = 0.0_rt;
        amrex::Real qv = 0.0_rt;
        amrex::Real qc = 0.0_rt;
        amrex::Real qi = 0.0_rt;

        ASSERT_GT(qw, 0.0_rt);
        ASSERT_GT(pdf_ql, 0.0_rt);
        ASSERT_LT(pdf_ql, qw);
        ASSERT_GT(unclamped_tabs, shoc::constants::min_temp());
        ASSERT_GT(unclamped_tabs, shoc::constants::freezing_temp());

        shoc::reconstruct_pdf_state(thetal, qw, exner, qi_seed, pdf_ql,
                                    tabs, qv, qc, qi);

        EXPECT_NEAR(qv + qc + qi, qw, 1.0e-14_rt);

        const amrex::Real theta = tabs / exner;
        const amrex::Real ql_total = qc + qi;
        const amrex::Real reconstructed_thetal =
            theta - (L_v / Cp_d) * ql_total / exner;

        EXPECT_NEAR(reconstructed_thetal, thetal, 1.0e-12_rt);
        EXPECT_GT(qc, 0.0_rt);
        EXPECT_NEAR(qi, 0.0_rt, 1.0e-14_rt);
        EXPECT_GT(qv, 0.0_rt);
    }

    {
        const amrex::Real thetal = 260.0_rt;
        const amrex::Real qw = 0.012_rt;
        const amrex::Real pdf_ql = 0.002_rt;
        const amrex::Real qi_seed = 1.0e-5_rt;
        const amrex::Real unclamped_tabs = thetal * exner + (L_v / Cp_d) * pdf_ql;
        amrex::Real tabs = 0.0_rt;
        amrex::Real qv = 0.0_rt;
        amrex::Real qc = 0.0_rt;
        amrex::Real qi = 0.0_rt;

        ASSERT_GT(qw, 0.0_rt);
        ASSERT_GT(pdf_ql, 0.0_rt);
        ASSERT_LT(pdf_ql, qw);
        ASSERT_GT(unclamped_tabs, shoc::constants::min_temp());
        ASSERT_LT(unclamped_tabs, shoc::constants::freezing_temp());

        shoc::reconstruct_pdf_state(thetal, qw, exner, qi_seed, pdf_ql,
                                    tabs, qv, qc, qi);

        EXPECT_NEAR(qv + qc + qi, qw, 1.0e-14_rt);
        EXPECT_NEAR(qc, pdf_ql, 1.0e-14_rt);
        EXPECT_NEAR(qi, qi_seed, 1.0e-14_rt);
        EXPECT_NEAR(qv, qw - pdf_ql - qi_seed, 1.0e-14_rt);
        EXPECT_NEAR(tabs, shoc::temperature_from_thetal(thetal, pdf_ql, qi_seed, exner),
                    1.0e-12_rt);
    }
}

TEST(ShocPhysical, DryUniformNoFluxColumnIsConserved)
{
    // A dry, horizontally/vertically uniform column with no surface fluxes,
    // no condensate, and no turbulent diffusivity should not create heat,
    // moisture, condensate, or TKE changes.
    auto col = shoc_test::make_column(6);
    ShocRuntimeOptions opts;
    const amrex::Real dt = 10.0_rt;

    zero_pdf_moments(col);
    set_uniform_implicit_column(col, 300.0_rt, 1.0_rt,
                                0.0_rt, 0.0_rt, 0.0_rt,
                                1.0e-3_rt, 2.0_rt, -1.0_rt);
    shoc::set_fab_val(col.surf_sens_flux, 0.0_rt, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 0.0_rt, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, 0.0_rt, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.0_rt, shoc::InitRunOn::Host);

    const amrex::Real energy_before = column_moist_energy(col);
    const amrex::Real water_before = column_total_water(col);

    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(col, opts, dt);
    });

    EXPECT_NEAR(column_moist_energy(col), energy_before,
                1.0e-12_rt * amrex::max(1.0_rt, amrex::Math::abs(energy_before)));
    EXPECT_NEAR(column_total_water(col), water_before, 1.0e-14_rt);

    const auto qv = col.qv.const_array();
    const auto qc = col.qc.const_array();
    const auto qi = col.qi.const_array();
    const auto qw = col.qw.const_array();
    const auto theta_tend = col.theta_tend.const_array();
    const auto qv_tend = col.qv_tend.const_array();
    const auto qc_tend = col.qc_tend.const_array();
    const auto qi_tend = col.qi_tend.const_array();
    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_NEAR(qv(0,k,0), 0.0_rt, 1.0e-14_rt);
        EXPECT_NEAR(qc(0,k,0), 0.0_rt, 1.0e-14_rt);
        EXPECT_NEAR(qi(0,k,0), 0.0_rt, 1.0e-14_rt);
        EXPECT_NEAR(qw(0,k,0), 0.0_rt, 1.0e-14_rt);
        EXPECT_NEAR(theta_tend(0,k,0), 0.0_rt, 1.0e-14_rt);
        EXPECT_NEAR(qv_tend(0,k,0), 0.0_rt, 1.0e-14_rt);
        EXPECT_NEAR(qc_tend(0,k,0), 0.0_rt, 1.0e-14_rt);
        EXPECT_NEAR(qi_tend(0,k,0), 0.0_rt, 1.0e-14_rt);
    }
}

TEST(ShocPhysical, ColumnWaterBudgetTracksSurfaceLatentFlux)
{
    // With a uniform moist column, no initial condensate, no internal gradients,
    // and zero turbulent diffusivity, the column-integrated total water change
    // should equal rho_surface * surface_latent_flux * dt.
    auto col = shoc_test::make_column(6);
    ShocRuntimeOptions opts;
    const amrex::Real dt = 10.0_rt;

    zero_pdf_moments(col);
    set_uniform_implicit_column(col, 300.0_rt, 1.0_rt,
                                0.008_rt, 0.0_rt, 0.0_rt,
                                1.0e-3_rt, 2.0_rt, -1.0_rt);
    shoc::set_fab_val(col.surf_sens_flux, 0.0_rt, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 1.0e-4_rt, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, 0.0_rt, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.0_rt, shoc::InitRunOn::Host);

    const amrex::Real water_before = column_total_water(col);
    const amrex::Real rho_sfc = col.rho.const_array()(0,0,0);
    const amrex::Real expected_delta = dt * rho_sfc * col.surf_lat_flux.const_array()(0,0,0);

    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(col, opts, dt);
    });

    const amrex::Real water_after = column_total_water(col);
    EXPECT_NEAR(water_after - water_before, expected_delta,
                1.0e-9_rt * amrex::max(1.0_rt, amrex::Math::abs(expected_delta)));

    const auto qv = col.qv.const_array();
    const auto qc = col.qc.const_array();
    const auto qi = col.qi.const_array();
    const auto qw = col.qw.const_array();
    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_GE(qv(0,k,0), -1.0e-14_rt);
        EXPECT_GE(qc(0,k,0), -1.0e-14_rt);
        EXPECT_GE(qi(0,k,0), -1.0e-14_rt);
        EXPECT_GE(qw(0,k,0), -1.0e-14_rt);
    }
}

TEST(ShocPhysical, BuoyancyProductionSignFollowsVirtualHeatFlux)
{
    // In the default TKE path, SHOC computes buoyancy production as
    // (g / theta0) * w'theta_v'. Stable virtual heat flux should produce
    // negative buoyancy production and should not increase TKE when shear is zero.
    // Unstable virtual heat flux should produce positive buoyancy production
    // and should increase TKE when it exceeds dissipation.
    auto stable = shoc_test::make_column(4);
    auto unstable = shoc_test::make_column(4);
    ShocRuntimeOptions opts;
    opts.shoc_1p5tke = false;
    opts.signed_tke_production = false;

    for (auto* col : {&stable, &unstable}) {
        auto u = col->u.array();
        auto v = col->v.array();
        auto shoc_mix = col->shoc_mix.array();
        auto tke = col->tke.array();
        auto tabs = col->tabs.array();
        auto p_mid = col->p_mid.array();
        auto brunt = col->brunt.array();
        auto tk = col->tk.array();
        auto tkh = col->tkh.array();
        auto pblh = col->pblh.array();

        pblh(0,0,0) = 1000.0_rt;
        for (int k = 0; k < col->layout.nlev; ++k) {
            u(0,k,0) = 2.0_rt;
            v(0,k,0) = -1.0_rt;
            shoc_mix(0,k,0) = 100.0_rt;
            tke(0,k,0) = 5.0e-4_rt;
            tabs(0,k,0) = 300.0_rt;
            p_mid(0,k,0) = 90000.0_rt;
            brunt(0,k,0) = 0.0_rt;
            tk(0,k,0) = 0.0_rt;
            tkh(0,k,0) = 0.0_rt;
        }
    }

    const amrex::Real stable_old_tke = stable.tke.const_array()(0,0,0);
    const amrex::Real unstable_old_tke = unstable.tke.const_array()(0,0,0);
    for (int k = 0; k < stable.layout.nlev; ++k) {
        stable.wthv_sec.array()(0,k,0) = -0.02_rt;
        unstable.wthv_sec.array()(0,k,0) = 0.05_rt;
    }

    shoc_test::run_and_sync([&] {
        ShocTKE::diagnose_tke_and_diffusivities(stable, opts, 10.0_rt);
        ShocTKE::diagnose_tke_and_diffusivities(unstable, opts, 10.0_rt);
    });

    for (int k = 0; k < stable.layout.nlev; ++k) {
        EXPECT_LT(stable.buoy_prod.const_array()(0,k,0), 0.0_rt);
        EXPECT_GT(unstable.buoy_prod.const_array()(0,k,0), 0.0_rt);
        EXPECT_LE(stable.tke.const_array()(0,k,0), stable_old_tke + 1.0e-14_rt);
        EXPECT_GT(unstable.tke.const_array()(0,k,0), unstable_old_tke);
        EXPECT_NEAR(stable.shear_prod.const_array()(0,k,0), 0.0_rt, 1.0e-14_rt);
        EXPECT_NEAR(unstable.shear_prod.const_array()(0,k,0), 0.0_rt, 1.0e-14_rt);
    }
}

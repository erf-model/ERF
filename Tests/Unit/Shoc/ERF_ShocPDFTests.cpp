#include "ERF_ShocMoments.H"
#include "ERF_ShocPDF.H"
#include "ERF_ShocStructure.H"
#include "ERF_ShocTKE.H"
#include "ERF_ShocTestUtils.H"
#include "ERF_ShocTypes.H"

#include <gtest/gtest.h>

#include <cmath>
#include <random>

namespace
{
void run_pdf_pipeline (ShocColumnData& col,
                       const ShocRuntimeOptions& opts,
                       amrex::Real dt)
{
    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
        ShocStructure::diagnose_pblh(col);
        ShocStructure::diagnose_length_and_brunt(col, opts, 500.0, 500.0);
        ShocTKE::diagnose_tke_and_diffusivities(col, opts, 300.0);
        ShocMoments::diagnose_moments(col, opts);
        ShocPDF::diagnose_pdf(col, opts, dt);
    });
}

void setup_translated_e3sm_pdf_column (ShocColumnData& col, amrex::Real w3_value)
{
    const auto fixture =
        shoc_test::read_named_fixture_vectors("pdf/e3sm_assumed_pdf_nondegenerate_columns.txt");
    const auto& thetal_in = shoc_test::fixture_values(fixture, "thetal");
    const auto& qw_in = shoc_test::fixture_values(fixture, "qw");
    const auto& pres_in = shoc_test::fixture_values(fixture, "pres");
    const auto& zi_in = shoc_test::fixture_values(fixture, "zi_grid");
    const auto& zt_in = shoc_test::fixture_values(fixture, "zt_grid");
    const auto& wsec_in = shoc_test::fixture_values(fixture, "w_sec");
    const auto& thl_sec_in = shoc_test::fixture_values(fixture, "thl_sec");
    const auto& qw_sec_in = shoc_test::fixture_values(fixture, "qw_sec");
    const auto& wthl_sec_in = shoc_test::fixture_values(fixture, "wthl_sec");
    const auto& wqw_sec_in = shoc_test::fixture_values(fixture, "wqw_sec");
    const auto& qwthl_sec_in = shoc_test::fixture_values(fixture, "qwthl_sec");
    const auto& w3_zero = shoc_test::fixture_values(fixture, "w3_symmetric");
    const auto& w3_pos = shoc_test::fixture_values(fixture, "w3_positive");
    const auto& w3_neg = shoc_test::fixture_values(fixture, "w3_negative");

    ASSERT_EQ(static_cast<int>(thetal_in.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(qw_in.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(pres_in.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(zi_in.size()), col.layout.nlev + 1);
    ASSERT_EQ(static_cast<int>(zt_in.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(wsec_in.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(thl_sec_in.size()), col.layout.nlev + 1);
    ASSERT_EQ(static_cast<int>(qw_sec_in.size()), col.layout.nlev + 1);
    ASSERT_EQ(static_cast<int>(wthl_sec_in.size()), col.layout.nlev + 1);
    ASSERT_EQ(static_cast<int>(wqw_sec_in.size()), col.layout.nlev + 1);
    ASSERT_EQ(static_cast<int>(qwthl_sec_in.size()), col.layout.nlev + 1);
    ASSERT_EQ(static_cast<int>(w3_zero.size()), col.layout.nlev + 1);
    ASSERT_EQ(static_cast<int>(w3_pos.size()), col.layout.nlev + 1);
    ASSERT_EQ(static_cast<int>(w3_neg.size()), col.layout.nlev + 1);

    auto thetal = col.thetal.array();
    auto qv = col.qv.array();
    auto qw = col.qw.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto tabs = col.tabs.array();
    auto p_mid = col.p_mid.array();
    auto zt = col.zt.array();
    auto zi = col.zi.array();
    auto w = col.w.array();
    auto w_sec = col.w_sec.array();
    auto thl_sec = col.thl_sec.array();
    auto qw_sec = col.qw_sec.array();
    auto wthl_sec = col.wthl_sec.array();
    auto wqw_sec = col.wqw_sec.array();
    auto qwthl_sec = col.qwthl_sec.array();
    auto w3 = col.w3.array();

    for (int k = 0; k <= col.layout.nlev; ++k) {
        zi(0,k,0) = zi_in[k];
    }

    for (int k = 0; k < col.layout.nlev; ++k) {
        thetal(0,k,0) = thetal_in[k];
        qv(0,k,0) = qw_in[k];
        qw(0,k,0) = qw_in[k];
        qc(0,k,0) = 0.0;
        qi(0,k,0) = 0.0;
        tabs(0,k,0) = thetal_in[k];
        p_mid(0,k,0) = pres_in[k];
        zt(0,k,0) = zt_in[k];
        w(0,k,0) = 0.0;
        w_sec(0,k,0) = wsec_in[k];
    }

    for (int k = 0; k <= col.layout.nlev; ++k) {
        thl_sec(0,k,0) = thl_sec_in[k];
        qw_sec(0,k,0) = qw_sec_in[k];
        wthl_sec(0,k,0) = wthl_sec_in[k];
        wqw_sec(0,k,0) = wqw_sec_in[k];
        qwthl_sec(0,k,0) = qwthl_sec_in[k];
        if (w3_value > 0.0) {
            w3(0,k,0) = w3_pos[k];
        } else if (w3_value < 0.0) {
            w3(0,k,0) = w3_neg[k];
        } else {
            w3(0,k,0) = w3_zero[k];
        }
    }
}
}

TEST(ShocPDF, DryUnsaturatedColumnStaysCloudFree)
{
    auto col = shoc_test::make_column(6);
    ShocRuntimeOptions opts;
    auto qv = col.qv.array();
    auto qw = col.qw.array();
    auto qc = col.qc.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        qv(0,k,0) = 0.002;
        qc(0,k,0) = 0.0;
        qw(0,k,0) = qv(0,k,0);
    }

    run_pdf_pipeline(col, opts, 10.0);

    const auto cldfrac = col.shoc_cldfrac.const_array();
    const auto ql = col.shoc_ql.const_array();
    const auto ql2 = col.shoc_ql2.const_array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_GE(cldfrac(0,k,0), 0.0);
        EXPECT_LE(cldfrac(0,k,0), 1.0);
        EXPECT_LE(ql(0,k,0), 1.0e-6);
        EXPECT_LE(ql2(0,k,0), 1.0e-8);
    }
}

TEST(ShocPDF, SupersaturatedColumnProducesLiquidAndCloudFraction)
{
    auto col = shoc_test::make_column(6);
    ShocRuntimeOptions opts;
    auto qv = col.qv.array();
    auto qw = col.qw.array();
    auto tabs = col.tabs.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        tabs(0,k,0) = 290.0;
        qv(0,k,0) = (k == 2) ? 0.030 : 0.012;
        qw(0,k,0) = qv(0,k,0);
    }

    run_pdf_pipeline(col, opts, 10.0);

    const auto cldfrac = col.shoc_cldfrac.const_array();
    const auto ql = col.shoc_ql.const_array();
    const auto ql2 = col.shoc_ql2.const_array();
    const auto wthv = col.wthv_sec.const_array();

    EXPECT_GT(cldfrac(0,2,0), 0.0);
    EXPECT_LE(cldfrac(0,2,0), 1.0);
    EXPECT_GT(ql(0,2,0), 0.0);
    EXPECT_GE(ql2(0,2,0), 0.0);
    EXPECT_TRUE(std::isfinite(wthv(0,2,0)));
}

TEST(ShocPDF, ExtraDiagnosticsTrackCondensationAndEvaporation)
{
    auto col = shoc_test::make_column(6);
    ShocRuntimeOptions opts;
    opts.extra_shoc_diags = true;
    auto qv = col.qv.array();
    auto qw = col.qw.array();
    auto tabs = col.tabs.array();
    auto shoc_ql = col.shoc_ql.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        tabs(0,k,0) = 289.0;
        qv(0,k,0) = (k == 3) ? 0.028 : 0.010;
        qw(0,k,0) = qv(0,k,0);
        shoc_ql(0,k,0) = (k == 3) ? 0.0 : 5.0e-4;
    }

    run_pdf_pipeline(col, opts, 20.0);

    const auto cond = col.shoc_cond.const_array();
    const auto evap = col.shoc_evap.const_array();
    EXPECT_GT(cond(0,3,0), 0.0);
    EXPECT_EQ(evap(0,3,0), 0.0);
}

TEST(ShocPDF, TranslatedE3smNoSgsVariabilityCaseKeepsDegenerateOutputs)
{
    auto col = shoc_test::make_column(5);
    ShocRuntimeOptions opts;

    auto thetal = col.thetal.array();
    auto qv = col.qv.array();
    auto qw = col.qw.array();
    auto tabs = col.tabs.array();
    auto p_mid = col.p_mid.array();
    auto zt = col.zt.array();
    auto zi = col.zi.array();
    auto w = col.w.array();
    auto w_sec = col.w_sec.array();
    auto thl_sec = col.thl_sec.array();
    auto qw_sec = col.qw_sec.array();
    auto wthl_sec = col.wthl_sec.array();
    auto wqw_sec = col.wqw_sec.array();
    auto qwthl_sec = col.qwthl_sec.array();
    auto w3 = col.w3.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();

    const amrex::Real thetal_in[5] = {300.0, 298.0, 298.0, 300.0, 303.0};
    const amrex::Real qw_in[5] = {0.017, 0.016, 0.011, 0.004, 0.003};
    const amrex::Real pres_in[5] = {100000.0, 90000.0, 85000.0, 80000.0, 70000.0};

    for (int k = 0; k <= col.layout.nlev; ++k) {
        zi(0,k,0) = 500.0 * k;
        if (k < col.layout.nlev) {
            zt(0,k,0) = 250.0 + 500.0 * k;
        }
    }

    for (int k = 0; k < col.layout.nlev; ++k) {
        thetal(0,k,0) = thetal_in[k];
        qv(0,k,0) = qw_in[k];
        qw(0,k,0) = qw_in[k];
        qc(0,k,0) = 0.0;
        qi(0,k,0) = 0.0;
        p_mid(0,k,0) = pres_in[k];
        tabs(0,k,0) = thetal_in[k];
        w(0,k,0) = 0.0;
        w_sec(0,k,0) = 0.004;
    }
    for (int k = 0; k <= col.layout.nlev; ++k) {
        thl_sec(0,k,0) = 0.0;
        qw_sec(0,k,0) = 0.0;
        wthl_sec(0,k,0) = 0.0;
        wqw_sec(0,k,0) = 0.0;
        qwthl_sec(0,k,0) = 0.0;
        w3(0,k,0) = 0.0;
    }

    shoc_test::run_and_sync([&] {
        ShocPDF::diagnose_pdf(col, opts, 10.0);
    });

    const auto cldfrac = col.shoc_cldfrac.const_array();
    const auto ql = col.shoc_ql.const_array();
    const auto ql2 = col.shoc_ql2.const_array();
    const auto wqls = col.wqls_sec.const_array();
    const auto wthv = col.wthv_sec.const_array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_TRUE(cldfrac(0,k,0) == 0.0 || cldfrac(0,k,0) == 1.0);
        EXPECT_DOUBLE_EQ(wqls(0,k,0), 0.0);
        EXPECT_DOUBLE_EQ(wthv(0,k,0), 0.0);
        EXPECT_NEAR(ql2(0,k,0), 0.0, 1.0e-14);
        EXPECT_GE(ql(0,k,0), 0.0);
    }
}

TEST(ShocPDF, RandomizedPdfDiagnosticsStayFiniteAndBounded)
{
    ShocRuntimeOptions opts;
    std::mt19937 gen(271828);

    for (int n = 0; n < 24; ++n) {
        auto col = shoc_test::make_randomized_column(8, gen);
        run_pdf_pipeline(col, opts, 15.0);

        const auto cldfrac = col.shoc_cldfrac.const_array();
        const auto ql = col.shoc_ql.const_array();
        const auto ql2 = col.shoc_ql2.const_array();
        const auto cond = col.shoc_cond.const_array();
        const auto evap = col.shoc_evap.const_array();
        const auto wqls = col.wqls_sec.const_array();
        const auto wthv = col.wthv_sec.const_array();

        for (int k = 0; k < col.layout.nlev; ++k) {
            EXPECT_TRUE(std::isfinite(cldfrac(0,k,0)));
            EXPECT_TRUE(std::isfinite(ql(0,k,0)));
            EXPECT_TRUE(std::isfinite(ql2(0,k,0)));
            EXPECT_TRUE(std::isfinite(cond(0,k,0)));
            EXPECT_TRUE(std::isfinite(evap(0,k,0)));
            EXPECT_TRUE(std::isfinite(wqls(0,k,0)));
            EXPECT_TRUE(std::isfinite(wthv(0,k,0)));
            EXPECT_GE(cldfrac(0,k,0), 0.0);
            EXPECT_LE(cldfrac(0,k,0), 1.0);
            EXPECT_GE(ql(0,k,0), 0.0);
            EXPECT_GE(ql2(0,k,0), 0.0);
            EXPECT_GE(cond(0,k,0), 0.0);
            EXPECT_GE(evap(0,k,0), 0.0);
        }
    }
}

TEST(ShocPDF, TranslatedE3smPositiveSkewnessCaseMatchesReferenceRelationships)
{
    ShocRuntimeOptions opts;
    auto symmetric = shoc_test::make_column(5);
    auto skewed = shoc_test::make_column(5);

    setup_translated_e3sm_pdf_column(symmetric, 0.0);
    setup_translated_e3sm_pdf_column(skewed, 1.0);

    shoc_test::run_and_sync([&] {
        ShocPDF::diagnose_pdf(symmetric, opts, 300.0);
        ShocPDF::diagnose_pdf(skewed, opts, 300.0);
    });

    const auto cf0 = symmetric.shoc_cldfrac.const_array();
    const auto cf1 = skewed.shoc_cldfrac.const_array();
    const auto ql0 = symmetric.shoc_ql.const_array();
    const auto ql1 = skewed.shoc_ql.const_array();
    const auto ql20 = symmetric.shoc_ql2.const_array();
    const auto ql21 = skewed.shoc_ql2.const_array();
    const auto wql0 = symmetric.wqls_sec.const_array();
    const auto wql1 = skewed.wqls_sec.const_array();
    const auto wthv0 = symmetric.wthv_sec.const_array();
    const auto wthv1 = skewed.wthv_sec.const_array();

    for (int k = 0; k < symmetric.layout.nlev; ++k) {
        EXPECT_GT(cf0(0,k,0), 0.0);
        EXPECT_LT(cf0(0,k,0), 1.0);
        EXPECT_GT(cf1(0,k,0), 0.0);
        EXPECT_LT(cf1(0,k,0), 1.0);
        EXPECT_GT(wql0(0,k,0), 0.0);
        EXPECT_GT(wthv0(0,k,0), 0.0);
        EXPECT_GT(ql20(0,k,0), 0.0);
        EXPECT_GT(ql21(0,k,0), 0.0);
        EXPECT_GT(ql0(0,k,0), 0.0);
        EXPECT_GT(ql1(0,k,0), 0.0);
        EXPECT_LT(std::abs(wql0(0,k,0)), 0.1);
        EXPECT_LT(std::abs(wql1(0,k,0)), 0.1);
        EXPECT_LT(std::abs(wthv0(0,k,0)), 10.0);
        EXPECT_LT(std::abs(wthv1(0,k,0)), 10.0);
        EXPECT_LT(ql20(0,k,0), 0.1);
        EXPECT_LT(ql21(0,k,0), 0.1);
        EXPECT_LT(ql0(0,k,0), 0.1);
        EXPECT_LT(ql1(0,k,0), 0.1);
        EXPECT_GT(wql1(0,k,0), 0.0);
        EXPECT_GT(wthv1(0,k,0), 0.0);
        EXPECT_LT(ql1(0,k,0), ql0(0,k,0));
        EXPECT_GT(wql1(0,k,0), wql0(0,k,0));
        EXPECT_GT(wthv1(0,k,0), wthv0(0,k,0));
        EXPECT_GT(ql21(0,k,0), ql20(0,k,0));
    }

}

TEST(ShocPDF, TranslatedE3smNegativeSkewnessCaseMatchesReferenceRelationships)
{
    ShocRuntimeOptions opts;
    auto symmetric = shoc_test::make_column(5);
    auto skewed = shoc_test::make_column(5);

    setup_translated_e3sm_pdf_column(symmetric, 0.0);
    setup_translated_e3sm_pdf_column(skewed, -1.0);

    shoc_test::run_and_sync([&] {
        ShocPDF::diagnose_pdf(symmetric, opts, 300.0);
        ShocPDF::diagnose_pdf(skewed, opts, 300.0);
    });

    const auto cf0 = symmetric.shoc_cldfrac.const_array();
    const auto cf1 = skewed.shoc_cldfrac.const_array();
    const auto ql0 = symmetric.shoc_ql.const_array();
    const auto ql1 = skewed.shoc_ql.const_array();
    const auto ql20 = symmetric.shoc_ql2.const_array();
    const auto ql21 = skewed.shoc_ql2.const_array();
    const auto wql0 = symmetric.wqls_sec.const_array();
    const auto wql1 = skewed.wqls_sec.const_array();
    const auto wthv0 = symmetric.wthv_sec.const_array();
    const auto wthv1 = skewed.wthv_sec.const_array();

    bool saw_flux_response = false;
    for (int k = 0; k < symmetric.layout.nlev; ++k) {
        EXPECT_GT(cf0(0,k,0), 0.0);
        EXPECT_LT(cf0(0,k,0), 1.0);
        EXPECT_GE(cf1(0,k,0), 0.0);
        EXPECT_LT(cf1(0,k,0), 1.0);
        EXPECT_GT(ql21(0,k,0), 0.0);
        EXPECT_GE(ql1(0,k,0), 0.0);
        EXPECT_LT(std::abs(wql1(0,k,0)), 0.1);
        EXPECT_LT(std::abs(wthv1(0,k,0)), 10.0);
        EXPECT_LT(ql21(0,k,0), 0.1);
        EXPECT_LT(ql1(0,k,0), 0.1);

        EXPECT_LE(ql1(0,k,0), ql0(0,k,0));
        EXPECT_TRUE(std::isfinite(wql1(0,k,0) - wql0(0,k,0)));
        EXPECT_TRUE(std::isfinite(wthv1(0,k,0) - wthv0(0,k,0)));
        saw_flux_response = saw_flux_response ||
            (std::abs(wql1(0,k,0) - wql0(0,k,0)) > 1.0e-12) ||
            (std::abs(wthv1(0,k,0) - wthv0(0,k,0)) > 1.0e-12);
        EXPECT_GT(ql21(0,k,0), ql20(0,k,0));
    }

    EXPECT_TRUE(saw_flux_response);
}

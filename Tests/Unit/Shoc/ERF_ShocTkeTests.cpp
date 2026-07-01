#include "ERF_ShocStructure.H"
#include "ERF_ShocTKE.H"
#include "ERF_ShocTestUtils.H"
#include "ERF_ShocTypes.H"

#include <gtest/gtest.h>

#include <cmath>
#include <random>

using namespace amrex;

TEST(ShocTke, ShearProductionRespectsBoundariesAndShear)
{
    auto col = shoc_test::make_column(5);
    FArrayBox sterm_iface;

    auto u = col.u.array();
    auto v = col.v.array();
    for (int k = 0; k < col.layout.nlev; ++k) {
        u(0,k,0,0) = 2.0 - k;
        v(0,k,0,0) = 1.0 + k;
    }

    shoc_test::run_and_sync([&] {
        ShocTKE::compute_shear_production(col, sterm_iface);
    });
    const auto sterm = sterm_iface.const_array();

    EXPECT_DOUBLE_EQ(sterm(0,0,0,0), 0.0);
    EXPECT_DOUBLE_EQ(sterm(0,col.layout.nlev,0,0), 0.0);
    for (int k = 1; k < col.layout.nlev; ++k) {
        EXPECT_GT(sterm(0,k,0,0), 0.0);
        EXPECT_TRUE(std::isfinite(sterm(0,k,0,0)));
    }

    for (int k = 0; k < col.layout.nlev; ++k) {
        u(0,k,0,0) = 3.0;
        v(0,k,0,0) = -2.0;
    }
    shoc_test::run_and_sync([&] {
        ShocTKE::compute_shear_production(col, sterm_iface);
    });
    for (int k = 0; k <= col.layout.nlev; ++k) {
        EXPECT_DOUBLE_EQ(sterm_iface.const_array()(0,k,0,0), 0.0);
    }
}

TEST(ShocTke, HelperKernelsMatchTranslatedE3smFixtures)
{
    auto col = shoc_test::make_column(5);
    FArrayBox sterm_iface;

    auto u = col.u.array();
    auto v = col.v.array();
    auto dz = col.dz.array();
    auto p_mid = col.p_mid.array();
    auto brunt = col.brunt.array();
    auto zt = col.zt.array();
    auto zi = col.zi.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        u(0,k,0,0) = 2.0 - k;
        v(0,k,0,0) = 1.0 + k;
        dz(0,k,0,0) = 100.0 + 25.0 * k;
        brunt(0,k,0,0) = 1.0e-3 * (k + 1);
        p_mid(0,k,0,0) = (k < 3) ? (95000.0 - 5000.0 * k) : (78000.0 - 2000.0 * (k - 3));
    }
    for (int k = 0; k <= col.layout.nlev; ++k) {
        zi(0,k,0,0) = 100.0 * k;
        if (k < col.layout.nlev) {
            zt(0,k,0,0) = 50.0 + 100.0 * k;
        }
    }

    shoc_test::run_and_sync([&] {
        ShocTKE::compute_shear_production(col, sterm_iface);
    });
    amrex::Vector<amrex::Real> brunt_int;
    ShocTKE::integrate_column_stability(col, brunt_int);

    const auto shear_fixture =
        shoc_test::read_fixture_vector("tke/e3sm_compute_shr_prod_erf_bottom_up.txt");
    const auto stab_fixture =
        shoc_test::read_fixture_vector("tke/e3sm_integ_column_stability_lower_trop.txt");

    ASSERT_EQ(shear_fixture.size(), static_cast<std::size_t>(col.layout.nlev + 1));
    ASSERT_EQ(stab_fixture.size(), 1);
    for (int k = 0; k <= col.layout.nlev; ++k) {
        EXPECT_NEAR(sterm_iface.const_array()(0,k,0,0), shear_fixture[k], 1.0e-15);
    }
    ASSERT_EQ(brunt_int.size(), 1);
    EXPECT_NEAR(brunt_int[0], stab_fixture[0], 1.0e-15);
}

TEST(ShocTke, ColumnStabilityIntegralUsesLowerTroposphereOnly)
{
    auto col = shoc_test::make_column(5);
    auto p_mid = col.p_mid.array();
    auto dz = col.dz.array();
    auto brunt = col.brunt.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        dz(0,k,0,0) = 100.0 + 25.0 * k;
        brunt(0,k,0,0) = 1.0e-3 * (k + 1);
        p_mid(0,k,0,0) = (k < 3) ? (95000.0 - 5000.0 * k) : (78000.0 - 2000.0 * (k - 3));
    }

    Vector<Real> brunt_int;
    ShocTKE::integrate_column_stability(col, brunt_int);

    const Real expected = dz(0,0,0,0) * brunt(0,0,0,0)
                        + dz(0,1,0,0) * brunt(0,1,0,0)
                        + dz(0,2,0,0) * brunt(0,2,0,0);
    ASSERT_EQ(brunt_int.size(), 1);
    EXPECT_NEAR(brunt_int[0], expected, 1.0e-12);
}

TEST(ShocTke, TimestepControlsTkeGrowth)
{
    ShocRuntimeOptions opts;
    auto col_small = shoc_test::make_column(6);
    auto col_large = shoc_test::make_column(6);

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col_small);
        ShocStructure::diagnose_pblh(col_small);
        ShocStructure::diagnose_length_and_brunt(col_small, opts, 500.0, 500.0);
        ShocStructure::diagnose_surface_layer(col_large);
        ShocStructure::diagnose_pblh(col_large);
        ShocStructure::diagnose_length_and_brunt(col_large, opts, 500.0, 500.0);
    });

    auto wthv_small = col_small.wthv_sec.array();
    auto wthv_large = col_large.wthv_sec.array();
    auto tk_small = col_small.tk.array();
    auto tk_large = col_large.tk.array();
    for (int k = 0; k < col_small.layout.nlev; ++k) {
        wthv_small(0,k,0,0) = 0.06 - 0.005 * k;
        wthv_large(0,k,0,0) = wthv_small(0,k,0,0);
        tk_small(0,k,0,0) = 1.0;
        tk_large(0,k,0,0) = 1.0;
    }

    const auto tke0 = col_small.tke.const_array()(0,0,0,0);
    shoc_test::run_and_sync([&] {
        ShocTKE::diagnose_tke_and_diffusivities(col_small, opts, 60.0);
        ShocTKE::diagnose_tke_and_diffusivities(col_large, opts, 300.0);
    });

    EXPECT_GT(col_small.tke.const_array()(0,0,0,0), tke0);
    EXPECT_GT(col_large.tke.const_array()(0,0,0,0), col_small.tke.const_array()(0,0,0,0));
}

TEST(ShocTke, TopTaperDampsUpperActiveLayerOnly)
{
    ShocRuntimeOptions opts_base;
    ShocRuntimeOptions opts_taper;
    opts_taper.top_taper_depth = 150.0;

    auto col_base = shoc_test::make_column(6);
    auto col_taper = shoc_test::make_column(6);

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col_base);
        ShocStructure::diagnose_pblh(col_base);
        ShocStructure::diagnose_length_and_brunt(col_base, opts_base, 500.0, 500.0);
        ShocStructure::diagnose_surface_layer(col_taper);
        ShocStructure::diagnose_pblh(col_taper);
        ShocStructure::diagnose_length_and_brunt(col_taper, opts_taper, 500.0, 500.0);
    });

    auto wthv_base = col_base.wthv_sec.array();
    auto wthv_taper = col_taper.wthv_sec.array();
    auto tk_base = col_base.tk.array();
    auto tk_taper = col_taper.tk.array();
    for (int k = 0; k < col_base.layout.nlev; ++k) {
        wthv_base(0,k,0,0) = 0.03;
        wthv_taper(0,k,0,0) = 0.03;
        tk_base(0,k,0,0) = 1.0;
        tk_taper(0,k,0,0) = 1.0;
    }

    shoc_test::run_and_sync([&] {
        ShocTKE::diagnose_tke_and_diffusivities(col_base, opts_base, 120.0);
        ShocTKE::diagnose_tke_and_diffusivities(col_taper, opts_taper, 120.0);
    });

    const int ktop = col_base.layout.nlev - 1;
    EXPECT_LT(col_taper.tke.const_array()(0,ktop,0,0), col_base.tke.const_array()(0,ktop,0,0));
    EXPECT_LT(col_taper.tk.const_array()(0,ktop,0,0), col_base.tk.const_array()(0,ktop,0,0));
    EXPECT_LT(col_taper.tkh.const_array()(0,ktop,0,0), col_base.tkh.const_array()(0,ktop,0,0));
    EXPECT_NEAR(col_taper.tke.const_array()(0,0,0,0), col_base.tke.const_array()(0,0,0,0), 1.0e-12);
}

TEST(ShocTke, PreviousDiffusivityFeedsShearProduction)
{
    ShocRuntimeOptions opts;
    auto col_zero = shoc_test::make_column(6);
    auto col_carried = shoc_test::make_column(6);

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col_zero);
        ShocStructure::diagnose_pblh(col_zero);
        ShocStructure::diagnose_length_and_brunt(col_zero, opts, 500.0, 500.0);
        ShocStructure::diagnose_surface_layer(col_carried);
        ShocStructure::diagnose_pblh(col_carried);
        ShocStructure::diagnose_length_and_brunt(col_carried, opts, 500.0, 500.0);
    });

    auto wthv_zero = col_zero.wthv_sec.array();
    auto wthv_carried = col_carried.wthv_sec.array();
    auto tk_zero = col_zero.tk.array();
    auto tk_carried = col_carried.tk.array();
    for (int k = 0; k < col_zero.layout.nlev; ++k) {
        wthv_zero(0,k,0,0) = 0.0;
        wthv_carried(0,k,0,0) = 0.0;
        tk_zero(0,k,0,0) = 0.0;
        tk_carried(0,k,0,0) = 1.0;
    }

    const Real tke0 = col_zero.tke.const_array()(0,1,0,0);
    shoc_test::run_and_sync([&] {
        ShocTKE::diagnose_tke_and_diffusivities(col_zero, opts, 120.0);
        ShocTKE::diagnose_tke_and_diffusivities(col_carried, opts, 120.0);
    });

    EXPECT_LT(col_zero.tke.const_array()(0,1,0,0), col_carried.tke.const_array()(0,1,0,0));
    EXPECT_LT(col_zero.tke.const_array()(0,1,0,0), tke0);
    EXPECT_GT(col_carried.tke.const_array()(0,1,0,0), col_zero.tke.const_array()(0,1,0,0));
}

TEST(ShocTke, OnePointFiveClosureUsesBruntInsteadOfBuoyancyFlux)
{
    auto col_default = shoc_test::make_column(6);
    auto col_15 = shoc_test::make_column(6);
    ShocRuntimeOptions opts_default;
    ShocRuntimeOptions opts_15;
    opts_15.shoc_1p5tke = true;

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col_default);
        ShocStructure::diagnose_pblh(col_default);
        ShocStructure::diagnose_length_and_brunt(col_default, opts_default, 500.0, 500.0);
        ShocStructure::diagnose_surface_layer(col_15);
        ShocStructure::diagnose_pblh(col_15);
        ShocStructure::diagnose_length_and_brunt(col_15, opts_15, 500.0, 500.0);
    });

    auto brunt = col_15.brunt.array();
    auto wthv_default = col_default.wthv_sec.array();
    auto wthv_15 = col_15.wthv_sec.array();
    auto tk_default = col_default.tk.array();
    auto tk_15 = col_15.tk.array();
    for (int k = 0; k < col_default.layout.nlev; ++k) {
        wthv_default(0,k,0,0) = 0.03;
        wthv_15(0,k,0,0) = 0.0;
        brunt(0,k,0,0) = 2.5e-3;
        tk_default(0,k,0,0) = 0.5;
        tk_15(0,k,0,0) = 0.5;
    }

    const Real tke_default_before = col_default.tke.const_array()(0,0,0,0);
    const Real tke_15_before = col_15.tke.const_array()(0,0,0,0);
    shoc_test::run_and_sync([&] {
        ShocTKE::diagnose_tke_and_diffusivities(col_default, opts_default, 120.0);
        ShocTKE::diagnose_tke_and_diffusivities(col_15, opts_15, 120.0);
    });

    EXPECT_LT(col_15.tke.const_array()(0,0,0,0), tke_15_before);
    EXPECT_GT(col_default.tke.const_array()(0,0,0,0), col_15.tke.const_array()(0,0,0,0));
    EXPECT_GT(col_default.tke.const_array()(0,0,0,0) - tke_default_before,
              col_15.tke.const_array()(0,0,0,0) - tke_15_before);
}

TEST(ShocTke, OnePointFiveClosureNegativeBruntActsAsPositiveBuoyancySource)
{
    auto col_default = shoc_test::make_column(6);
    auto col_15 = shoc_test::make_column(6);
    ShocRuntimeOptions opts_default;
    ShocRuntimeOptions opts_15;
    opts_15.shoc_1p5tke = true;

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col_default);
        ShocStructure::diagnose_pblh(col_default);
        ShocStructure::diagnose_length_and_brunt(col_default, opts_default, 500.0, 500.0);
        ShocStructure::diagnose_surface_layer(col_15);
        ShocStructure::diagnose_pblh(col_15);
        ShocStructure::diagnose_length_and_brunt(col_15, opts_15, 500.0, 500.0);
    });

    auto brunt_default = col_default.brunt.array();
    auto brunt_15 = col_15.brunt.array();
    auto wthv_default = col_default.wthv_sec.array();
    auto wthv_15 = col_15.wthv_sec.array();
    auto tk_default = col_default.tk.array();
    auto tk_15 = col_15.tk.array();
    auto tke_default = col_default.tke.array();
    auto tke_15 = col_15.tke.array();
    for (int k = 0; k < col_default.layout.nlev; ++k) {
        brunt_default(0,k,0,0) = -2.5e-3;
        brunt_15(0,k,0,0) = -2.5e-3;
        wthv_default(0,k,0,0) = 0.0;
        wthv_15(0,k,0,0) = 0.0;
        tk_default(0,k,0,0) = 0.0;
        tk_15(0,k,0,0) = 0.0;
        tke_default(0,k,0,0) = 1.0e-3;
        tke_15(0,k,0,0) = 1.0e-3;
    }

    const Real tke_default_before = col_default.tke.const_array()(0,0,0,0);
    const Real tke_15_before = col_15.tke.const_array()(0,0,0,0);
    shoc_test::run_and_sync([&] {
        ShocTKE::diagnose_tke_and_diffusivities(col_default, opts_default, 120.0);
        ShocTKE::diagnose_tke_and_diffusivities(col_15, opts_15, 120.0);
    });

    EXPECT_NEAR(col_default.buoy_prod.const_array()(0,0,0,0), 0.0, 1.0e-14);
    EXPECT_GT(col_15.buoy_prod.const_array()(0,0,0,0), 0.0);
    EXPECT_LT(col_default.tke.const_array()(0,0,0,0), tke_default_before);
    EXPECT_GT(col_15.tke.const_array()(0,0,0,0), tke_15_before);
    EXPECT_GT(col_15.tke.const_array()(0,0,0,0), col_default.tke.const_array()(0,0,0,0));
}

TEST(ShocTke, SignedProductionOptionAllowsStableBuoyancyToReduceTke)
{
    auto col_default = shoc_test::make_column(6);
    auto col_signed = shoc_test::make_column(6);
    ShocRuntimeOptions opts_default;
    ShocRuntimeOptions opts_signed;
    opts_signed.signed_tke_production = true;

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col_default);
        ShocStructure::diagnose_pblh(col_default);
        ShocStructure::diagnose_length_and_brunt(col_default, opts_default, 500.0, 500.0);
        ShocStructure::diagnose_surface_layer(col_signed);
        ShocStructure::diagnose_pblh(col_signed);
        ShocStructure::diagnose_length_and_brunt(col_signed, opts_signed, 500.0, 500.0);
    });

    auto wthv_default = col_default.wthv_sec.array();
    auto wthv_signed = col_signed.wthv_sec.array();
    auto tk_default = col_default.tk.array();
    auto tk_signed = col_signed.tk.array();
    for (int k = 0; k < col_default.layout.nlev; ++k) {
        wthv_default(0,k,0,0) = -0.2;
        wthv_signed(0,k,0,0) = -0.2;
        tk_default(0,k,0,0) = 0.0;
        tk_signed(0,k,0,0) = 0.0;
    }

    const Real tke_default_before = col_default.tke.const_array()(0,0,0,0);
    const Real tke_signed_before = col_signed.tke.const_array()(0,0,0,0);
    shoc_test::run_and_sync([&] {
        ShocTKE::diagnose_tke_and_diffusivities(col_default, opts_default, 60.0);
        ShocTKE::diagnose_tke_and_diffusivities(col_signed, opts_signed, 60.0);
    });

    EXPECT_LT(col_default.tke.const_array()(0,0,0,0), tke_default_before);
    EXPECT_LT(col_signed.tke.const_array()(0,0,0,0), tke_signed_before);
    EXPECT_LT(col_signed.tke.const_array()(0,0,0,0), col_default.tke.const_array()(0,0,0,0));
    EXPECT_GE(col_signed.tke.const_array()(0,0,0,0), 4.0e-4);
}

TEST(ShocTke, DiffusivitiesRemainPositiveAndFinite)
{
    auto col = shoc_test::make_column(6);
    ShocRuntimeOptions opts;

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
        ShocStructure::diagnose_pblh(col);
        ShocStructure::diagnose_length_and_brunt(col, opts, 500.0, 500.0);
        ShocTKE::diagnose_tke_and_diffusivities(col, opts, 300.0);
    });

    const auto tke = col.tke.const_array();
    const auto tk = col.tk.const_array();
    const auto tkh = col.tkh.const_array();
    const auto isotropy = col.isotropy.const_array();
    const auto tend = col.tke_tend.const_array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_TRUE(std::isfinite(tke(0,k,0,0)));
        EXPECT_TRUE(std::isfinite(tk(0,k,0,0)));
        EXPECT_TRUE(std::isfinite(tkh(0,k,0,0)));
        EXPECT_TRUE(std::isfinite(isotropy(0,k,0,0)));
        EXPECT_TRUE(std::isfinite(tend(0,k,0,0)));
        EXPECT_GE(tke(0,k,0,0), 4.0e-4);
        EXPECT_LE(tke(0,k,0,0), 50.0);
        EXPECT_GE(tk(0,k,0,0), 0.0);
        EXPECT_GE(tkh(0,k,0,0), 0.0);
        EXPECT_GE(isotropy(0,k,0,0), 0.0);
    }
}

TEST(ShocTke, RandomizedColumnsStayBounded)
{
    ShocRuntimeOptions opts;
    std::mt19937 gen(1729);

    for (int n = 0; n < 32; ++n) {
        auto col = shoc_test::make_randomized_column(8, gen);
        shoc_test::run_and_sync([&] {
            ShocStructure::diagnose_surface_layer(col);
            ShocStructure::diagnose_pblh(col);
            ShocStructure::diagnose_length_and_brunt(col, opts, 600.0, 450.0);
            ShocTKE::diagnose_tke_and_diffusivities(col, opts, 180.0);
        });

        const auto mix = col.shoc_mix.const_array();
        const auto tke = col.tke.const_array();
        const auto tk = col.tk.const_array();
        const auto tkh = col.tkh.const_array();

        for (int k = 0; k < col.layout.nlev; ++k) {
            EXPECT_TRUE(std::isfinite(mix(0,k,0,0)));
            EXPECT_TRUE(std::isfinite(tke(0,k,0,0)));
            EXPECT_TRUE(std::isfinite(tk(0,k,0,0)));
            EXPECT_TRUE(std::isfinite(tkh(0,k,0,0)));
            EXPECT_GE(mix(0,k,0,0), 20.0);
            EXPECT_LE(mix(0,k,0,0), std::sqrt(600.0 * 450.0));
            EXPECT_GE(tke(0,k,0,0), 4.0e-4);
            EXPECT_LE(tke(0,k,0,0), 50.0);
            EXPECT_GE(tk(0,k,0,0), 0.0);
            EXPECT_GE(tkh(0,k,0,0), 0.0);
        }
    }
}

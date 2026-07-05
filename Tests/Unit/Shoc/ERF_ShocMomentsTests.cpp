#include "ERF_ShocMoments.H"
#include "ERF_ShocStructure.H"
#include "ERF_ShocTKE.H"
#include "ERF_ShocTestUtils.H"
#include "ERF_ShocTypes.H"

#include <gtest/gtest.h>

#include <cmath>
#include <random>

TEST(ShocMoments, SecondMomentBoundaryConditionsAreApplied)
{
    auto col = shoc_test::make_column(6);
    ShocRuntimeOptions opts;

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
        ShocStructure::diagnose_pblh(col);
        ShocStructure::diagnose_length_and_brunt(col, opts, 400.0, 400.0);
        ShocTKE::diagnose_tke_and_diffusivities(col, opts, 300.0);
        ShocMoments::diagnose_second_moments(col, opts);
    });

    const auto thl_sec = col.thl_sec.const_array();
    const auto qw_sec = col.qw_sec.const_array();
    const auto qwthl = col.qwthl_sec.const_array();
    const auto wthl = col.wthl_sec.const_array();
    const auto wqw = col.wqw_sec.const_array();
    const auto uw = col.uw_sec.const_array();
    const auto vw = col.vw_sec.const_array();
    const auto wtke = col.wtke_sec.const_array();

    EXPECT_GE(thl_sec(0,0,0), 0.0);
    EXPECT_GE(qw_sec(0,0,0), 0.0);
    EXPECT_TRUE(std::isfinite(qwthl(0,0,0)));
    EXPECT_DOUBLE_EQ(wthl(0,0,0), col.surf_sens_flux.const_array()(0,0,0));
    EXPECT_DOUBLE_EQ(wqw(0,0,0), col.surf_lat_flux.const_array()(0,0,0));
    EXPECT_DOUBLE_EQ(uw(0,0,0), col.surf_tau_u.const_array()(0,0,0));
    EXPECT_DOUBLE_EQ(vw(0,0,0), col.surf_tau_v.const_array()(0,0,0));
    EXPECT_GE(wtke(0,0,0), 0.0);

    const int ktop = col.layout.nlev;
    EXPECT_DOUBLE_EQ(thl_sec(0,ktop,0), 0.0);
    EXPECT_DOUBLE_EQ(qw_sec(0,ktop,0), 0.0);
    EXPECT_DOUBLE_EQ(qwthl(0,ktop,0), 0.0);
    EXPECT_DOUBLE_EQ(wthl(0,col.layout.nlev,0), 0.0);
    EXPECT_DOUBLE_EQ(wqw(0,col.layout.nlev,0), 0.0);
    EXPECT_DOUBLE_EQ(uw(0,col.layout.nlev,0), 0.0);
    EXPECT_DOUBLE_EQ(vw(0,col.layout.nlev,0), 0.0);
    EXPECT_DOUBLE_EQ(wtke(0,col.layout.nlev,0), 0.0);
}

TEST(ShocMoments, SurfaceMomentBoundaryConditionsMatchTranslatedE3smSemantics)
{
    auto col = shoc_test::make_column(6);
    ShocRuntimeOptions opts;

    shoc::set_fab_val(col.surf_sens_flux, 0.02, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 1.0e-4, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, 0.04, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.03, shoc::InitRunOn::Host);

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
        ShocStructure::diagnose_pblh(col);
        ShocStructure::diagnose_length_and_brunt(col, opts, 400.0, 400.0);
        ShocTKE::diagnose_tke_and_diffusivities(col, opts, 300.0);
        ShocMoments::diagnose_second_moments(col, opts);
    });

    const amrex::Real ustar2 = std::sqrt(0.04 * 0.04 + 0.03 * 0.03);
    const amrex::Real wstar = std::cbrt((CONST_GRAV / 300.0) * 0.02);
    const amrex::Real uf = amrex::max(amrex::Real(0.01), std::sqrt(ustar2 + amrex::Real(0.3) * wstar * wstar));
    const amrex::Real expected_thl_sec = 0.72 * std::pow(0.02 / uf, 2);
    const amrex::Real expected_qw_sec = 0.72 * std::pow(1.0e-4 / uf, 2);
    const amrex::Real expected_qwthl = 0.36 * (0.02 / uf) * (1.0e-4 / uf);
    const amrex::Real expected_wtke = std::pow(amrex::max(std::sqrt(ustar2), amrex::Real(0.01)), 3);

    const auto thl_sec = col.thl_sec.const_array();
    const auto qw_sec = col.qw_sec.const_array();
    const auto qwthl = col.qwthl_sec.const_array();
    const auto wtke = col.wtke_sec.const_array();

    EXPECT_NEAR(thl_sec(0,0,0), expected_thl_sec, 1.0e-12);
    EXPECT_NEAR(qw_sec(0,0,0), expected_qw_sec, 1.0e-12);
    EXPECT_NEAR(qwthl(0,0,0), expected_qwthl, 1.0e-12);
    EXPECT_NEAR(wtke(0,0,0), expected_wtke, 1.0e-12);
}

TEST(ShocMoments, TopTaperDampsUpperMomentDiagnosticsOnly)
{
    auto col_base = shoc_test::make_column(6);
    auto col_taper = shoc_test::make_column(6);
    ShocRuntimeOptions opts_base;
    ShocRuntimeOptions opts_taper;
    opts_taper.top_taper_depth = 150.0;

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col_base);
        ShocStructure::diagnose_pblh(col_base);
        ShocStructure::diagnose_length_and_brunt(col_base, opts_base, 400.0, 400.0);
        ShocTKE::diagnose_tke_and_diffusivities(col_base, opts_base, 300.0);
        ShocMoments::diagnose_moments(col_base, opts_base);
        ShocStructure::diagnose_surface_layer(col_taper);
        ShocStructure::diagnose_pblh(col_taper);
        ShocStructure::diagnose_length_and_brunt(col_taper, opts_taper, 400.0, 400.0);
        ShocTKE::diagnose_tke_and_diffusivities(col_taper, opts_taper, 300.0);
        ShocMoments::diagnose_moments(col_taper, opts_taper);
    });

    const int ktop_cell = col_base.layout.nlev - 1;
    const int ktop_iface = col_base.layout.nlev - 1;
    EXPECT_LT(col_taper.w_sec.const_array()(0,ktop_cell,0),
              col_base.w_sec.const_array()(0,ktop_cell,0));
    EXPECT_LT(std::abs(col_taper.wthl_sec.const_array()(0,ktop_iface,0)),
              std::abs(col_base.wthl_sec.const_array()(0,ktop_iface,0)));
    EXPECT_NEAR(col_taper.w_sec.const_array()(0,0,0),
                col_base.w_sec.const_array()(0,0,0), 1.0e-12);
}

TEST(ShocMoments, VarianceAndFluxHelpersMatchDowngradientForm)
{
    auto col = shoc_test::make_column(5);
    const amrex::Box iface_box(amrex::IntVect(0,0,0), amrex::IntVect(col.layout.ncell - 1, col.layout.nlev, 0));
    amrex::FArrayBox isotropy_zi(iface_box, 1);
    amrex::FArrayBox tkh_zi(iface_box, 1);
    amrex::FArrayBox outvar(iface_box, 1);
    amrex::FArrayBox flux(iface_box, 1);

    shoc::set_fab_val(isotropy_zi, 2.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(tkh_zi, 4.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(outvar, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(flux, 0.0, shoc::InitRunOn::Host);

    auto in1 = col.thetal.array();
    auto in2 = col.qw.array();
    for (int k = 0; k < col.layout.nlev; ++k) {
        in1(0,k,0) = 300.0 - 2.0 * k;
        in2(0,k,0) = 0.010 + 1.0e-3 * k;
    }

    shoc_test::run_and_sync([&] {
        ShocMoments::calc_var_or_covar(col, 1.5, isotropy_zi, tkh_zi, col.thetal, col.qw, outvar);
        ShocMoments::calc_vertflux(col, tkh_zi, col.thetal, flux);
    });

    const auto out = outvar.const_array();
    const auto vf = flux.const_array();
    const amrex::Real dz = col.dz.const_array()(0,0,0);
    const amrex::Real expected_covar = 1.5 * (2.0 * 4.0) * (1.0 / (dz * dz)) * (2.0) * (-1.0e-3);
    const amrex::Real expected_flux = (4.0 / dz) * 2.0;

    EXPECT_NEAR(out(0,1,0), expected_covar, 1.0e-12);
    EXPECT_NEAR(vf(0,1,0), expected_flux, 1.0e-12);
}

TEST(ShocMoments, HelperKernelsMatchTranslatedE3smFixtures)
{
    auto col = shoc_test::make_column(4);
    const amrex::Box iface_box(amrex::IntVect(0,0,0), amrex::IntVect(col.layout.ncell - 1, col.layout.nlev, 0));
    amrex::FArrayBox isotropy_zi(iface_box, 1);
    amrex::FArrayBox tkh_zi(iface_box, 1);
    amrex::FArrayBox outvar(iface_box, 1);
    amrex::FArrayBox flux(iface_box, 1);

    shoc::set_fab_val(isotropy_zi, 0.5, shoc::InitRunOn::Host);
    shoc::set_fab_val(tkh_zi, 2.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(outvar, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(flux, 0.0, shoc::InitRunOn::Host);

    auto zt = col.zt.array();
    auto zi = col.zi.array();
    auto in1 = col.thetal.array();
    auto in2 = col.qw.array();
    for (int k = 0; k <= col.layout.nlev; ++k) {
        zi(0,k,0) = 100.0 * k;
        if (k < col.layout.nlev) {
            zt(0,k,0) = 50.0 + 100.0 * k;
        }
    }

    in1(0,0,0) = 10.0;
    in1(0,1,0) = 8.0;
    in1(0,2,0) = 5.0;
    in1(0,3,0) = 1.0;
    in2(0,0,0) = 4.0;
    in2(0,1,0) = 3.0;
    in2(0,2,0) = 1.0;
    in2(0,3,0) = 0.0;

    shoc_test::run_and_sync([&] {
        ShocMoments::calc_var_or_covar(col, 1.7, isotropy_zi, tkh_zi, col.thetal, col.qw, outvar);
        ShocMoments::calc_vertflux(col, tkh_zi, col.thetal, flux);
    });

    const auto covar_fixture =
        shoc_test::read_fixture_vector("moments/e3sm_varorcovar_erf_bottom_up.txt");
    const auto flux_fixture =
        shoc_test::read_fixture_vector("moments/e3sm_vertflux_erf_bottom_up.txt");

    ASSERT_EQ(covar_fixture.size(), 3);
    ASSERT_EQ(flux_fixture.size(), 3);
    for (int k = 1; k < col.layout.nlev; ++k) {
        EXPECT_NEAR(outvar.const_array()(0,k,0), covar_fixture[k-1], 1.0e-15);
        EXPECT_NEAR(flux.const_array()(0,k,0), flux_fixture[k-1], 1.0e-15);
    }
}

TEST(ShocMoments, OnePointFiveTkeModeZerosVarianceAndThirdMoment)
{
    auto col = shoc_test::make_column(5);
    ShocRuntimeOptions opts;
    opts.shoc_1p5tke = true;

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
        ShocStructure::diagnose_pblh(col);
        ShocStructure::diagnose_length_and_brunt(col, opts, 450.0, 450.0);
        ShocTKE::diagnose_tke_and_diffusivities(col, opts, 300.0);
        ShocMoments::diagnose_moments(col, opts);
    });

    const auto w_sec = col.w_sec.const_array();
    const auto thl_sec = col.thl_sec.const_array();
    const auto qw_sec = col.qw_sec.const_array();
    const auto qwthl = col.qwthl_sec.const_array();
    const auto w3 = col.w3.const_array();

    EXPECT_GE(thl_sec(0,0,0), 0.0);
    EXPECT_GE(qw_sec(0,0,0), 0.0);
    EXPECT_TRUE(std::isfinite(qwthl(0,0,0)));

    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_DOUBLE_EQ(w_sec(0,k,0), 0.0);
        if (k > 0) {
            EXPECT_DOUBLE_EQ(thl_sec(0,k,0), 0.0);
            EXPECT_DOUBLE_EQ(qw_sec(0,k,0), 0.0);
            EXPECT_DOUBLE_EQ(qwthl(0,k,0), 0.0);
        }
    }
    EXPECT_DOUBLE_EQ(thl_sec(0,col.layout.nlev,0), 0.0);
    EXPECT_DOUBLE_EQ(qw_sec(0,col.layout.nlev,0), 0.0);
    EXPECT_DOUBLE_EQ(qwthl(0,col.layout.nlev,0), 0.0);
    for (int k = 0; k <= col.layout.nlev; ++k) {
        EXPECT_DOUBLE_EQ(w3(0,k,0), 0.0);
    }
}

TEST(ShocMoments, ThirdMomentClippingUsesPositiveFallback)
{
    auto col = shoc_test::make_column(5);
    const amrex::Box iface_box(amrex::IntVect(0,0,0), amrex::IntVect(col.layout.ncell - 1, col.layout.nlev, 0));
    amrex::FArrayBox w_sec_zi(iface_box, 1);
    amrex::FArrayBox w3(iface_box, 1);
    shoc::set_fab_val(w_sec_zi, 0.1, shoc::InitRunOn::Host);
    shoc::set_fab_val(w3, -10.0, shoc::InitRunOn::Host);

    shoc_test::run_and_sync([&] {
        ShocMoments::clip_third_moments(col, w_sec_zi, w3);
    });

    for (int k = 0; k <= col.layout.nlev; ++k) {
        EXPECT_DOUBLE_EQ(w3.const_array()(0,k,0), 0.02);
    }
}

TEST(ShocMoments, ThirdMomentCenteredDifferencesRespectBottomUpOrientation)
{
    auto col = shoc_test::make_column(4);
    ShocRuntimeOptions opts;

    auto dz = col.dz.array();
    auto zt = col.zt.array();
    auto zi = col.zi.array();
    auto thetal = col.thetal.array();
    auto isotropy = col.isotropy.array();
    auto brunt = col.brunt.array();
    auto w_sec = col.w_sec.array();
    auto tke = col.tke.array();
    auto thl_sec = col.thl_sec.array();
    auto wthl_sec = col.wthl_sec.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        dz(0,k,0) = 100.0;
        zt(0,k,0) = 50.0 + 100.0 * k;
        thetal(0,k,0) = 300.0;
        isotropy(0,k,0) = 2.0;
        brunt(0,k,0) = 0.0;
        w_sec(0,k,0) = 0.2 + 0.2 * k;
        tke(0,k,0) = 0.3 + 0.2 * k;
    }
    for (int k = 0; k <= col.layout.nlev; ++k) {
        zi(0,k,0) = 100.0 * k;
        thl_sec(0,k,0) = 0.2 * k;
        wthl_sec(0,k,0) = 0.05 * k;
    }

    shoc_test::run_and_sync([&] {
        ShocMoments::diagnose_third_moments(col, opts);
    });

    const amrex::Real c = opts.c_diag_3rd_mom;
    const amrex::Real a0 = (0.52 / (c * c)) / (c - 2.0);
    const amrex::Real a1 = 0.87 / (c * c);
    const amrex::Real a2 = 0.5 / c;
    const amrex::Real a4 = 2.4 / (3.0 * c + 5.0);

    const amrex::Real thedz = 1.0 / 100.0;
    const amrex::Real thedz2 = 1.0 / 200.0;
    const amrex::Real iso = 2.0;
    const amrex::Real isosq = iso * iso;
    const amrex::Real bet2 = CONST_GRAV / 300.0;
    const amrex::Real w_sec_zi = 0.3;
    const amrex::Real wthl_k = 0.05;
    const amrex::Real thl_sec_diff = 0.4;
    const amrex::Real wthl_sec_diff = 0.1;
    const amrex::Real wsec_diff = 0.2;
    const amrex::Real tke_diff = 0.2;

    const amrex::Real f0 = thedz2 * std::pow(bet2, 3) * std::pow(isosq, 2) * wthl_k * thl_sec_diff;
    const amrex::Real f1 = thedz2 * (bet2 * bet2) * std::pow(iso, 3) *
                            (wthl_k * wthl_sec_diff + 0.5 * w_sec_zi * thl_sec_diff);
    const amrex::Real f2 = thedz * bet2 * isosq * wthl_k * wsec_diff +
                            2.0 * thedz2 * bet2 * isosq * w_sec_zi * wthl_sec_diff;
    const amrex::Real f3 = thedz2 * bet2 * isosq * w_sec_zi * wthl_sec_diff +
                            thedz * bet2 * isosq * (wthl_k * tke_diff);
    const amrex::Real f4 = thedz * iso * w_sec_zi * (wsec_diff + tke_diff);
    const amrex::Real f5 = thedz * iso * w_sec_zi * wsec_diff;

    const amrex::Real x1 = a0 * f0 + a1 * f1 + a2 * f2;
    const amrex::Real y1 = 2.0 * a2 * ((a0 / a1) * f0 + f1);
    const amrex::Real omega0 = a4;
    const amrex::Real omega1 = omega0 / (2.0 * c);
    const amrex::Real omega2 = omega1 * f3 + 1.25 * omega0 * f4;
    const amrex::Real aa1 = omega0 * x1 + omega1 * y1 + omega2;
    const amrex::Real expected_w3 = (aa1 - 1.2 * x1 - 1.5 * f5) / c;

    EXPECT_NEAR(col.w3.const_array()(0,1,0), expected_w3, 1.0e-12);
}

TEST(ShocMoments, RandomizedMomentsStayFiniteAndBounded)
{
    ShocRuntimeOptions opts;
    std::mt19937 gen(314159);

    for (int n = 0; n < 24; ++n) {
        auto col = shoc_test::make_randomized_column(8, gen);
        shoc_test::run_and_sync([&] {
            ShocStructure::diagnose_surface_layer(col);
            ShocStructure::diagnose_pblh(col);
            ShocStructure::diagnose_length_and_brunt(col, opts, 550.0, 350.0);
            ShocTKE::diagnose_tke_and_diffusivities(col, opts, 300.0);
            ShocMoments::diagnose_moments(col, opts);
        });

        const auto thl_sec = col.thl_sec.const_array();
        const auto qw_sec = col.qw_sec.const_array();
        const auto qwthl = col.qwthl_sec.const_array();
        const auto wthl = col.wthl_sec.const_array();
        const auto wqw = col.wqw_sec.const_array();
        const auto uw = col.uw_sec.const_array();
        const auto vw = col.vw_sec.const_array();
        const auto wtke = col.wtke_sec.const_array();
        const auto w_sec = col.w_sec.const_array();
        const auto w3 = col.w3.const_array();

        for (int k = 0; k <= col.layout.nlev; ++k) {
            EXPECT_TRUE(std::isfinite(thl_sec(0,k,0)));
            EXPECT_TRUE(std::isfinite(qw_sec(0,k,0)));
            EXPECT_TRUE(std::isfinite(qwthl(0,k,0)));
            EXPECT_GE(thl_sec(0,k,0), 0.0);
            EXPECT_GE(qw_sec(0,k,0), 0.0);
        }

        for (int k = 0; k < col.layout.nlev; ++k) {
            EXPECT_TRUE(std::isfinite(w_sec(0,k,0)));
            EXPECT_GE(w_sec(0,k,0), 0.0);
        }

        for (int k = 0; k <= col.layout.nlev; ++k) {
            EXPECT_TRUE(std::isfinite(wthl(0,k,0)));
            EXPECT_TRUE(std::isfinite(wqw(0,k,0)));
            EXPECT_TRUE(std::isfinite(uw(0,k,0)));
            EXPECT_TRUE(std::isfinite(vw(0,k,0)));
            EXPECT_TRUE(std::isfinite(wtke(0,k,0)));
            EXPECT_TRUE(std::isfinite(w3(0,k,0)));
        }
    }
}

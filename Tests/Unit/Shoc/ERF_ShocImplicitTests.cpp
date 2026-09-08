#include "ERF_ShocImplicit.H"
#include "ERF_ShocEnergyFixer.H"
#include "ERF_ShocTestUtils.H"

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>

namespace
{
void load_translated_e3sm_implicit_multicolumn_fixture (ShocColumnData& col)
{
    const auto fixture = shoc_test::read_named_fixture_vectors(
        "implicit_energy/e3sm_update_prognostics_implicit_multicolumn.txt");
    const auto& zi_grid = shoc_test::fixture_values(fixture, "zi_grid");
    const auto& zt_grid = shoc_test::fixture_values(fixture, "zt_grid");
    const auto& dz_zt = shoc_test::fixture_values(fixture, "dz_zt");
    const auto& rho_zt = shoc_test::fixture_values(fixture, "rho_zt");
    const auto& tkh = shoc_test::fixture_values(fixture, "tkh");
    const auto& thetal_in = shoc_test::fixture_values(fixture, "thetal");
    const auto& qw_in = shoc_test::fixture_values(fixture, "qw");
    const auto& u_in = shoc_test::fixture_values(fixture, "u");
    const auto& v_in = shoc_test::fixture_values(fixture, "v");
    const auto& tke_in = shoc_test::fixture_values(fixture, "tke");
    const auto& wthl_sfc = shoc_test::fixture_values(fixture, "wthl_sfc");
    const auto& wqw_sfc = shoc_test::fixture_values(fixture, "wqw_sfc");
    const auto& uw_sfc = shoc_test::fixture_values(fixture, "uw_sfc");
    const auto& vw_sfc = shoc_test::fixture_values(fixture, "vw_sfc");

    ASSERT_EQ(static_cast<int>(zi_grid.size()), col.layout.nlev + 1);
    ASSERT_EQ(static_cast<int>(zt_grid.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(dz_zt.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(rho_zt.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(tkh.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(thetal_in.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(qw_in.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(u_in.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(v_in.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(tke_in.size()), col.layout.nlev);
    ASSERT_EQ(static_cast<int>(wthl_sfc.size()), col.layout.ncell);
    ASSERT_EQ(static_cast<int>(wqw_sfc.size()), col.layout.ncell);
    ASSERT_EQ(static_cast<int>(uw_sfc.size()), col.layout.ncell);
    ASSERT_EQ(static_cast<int>(vw_sfc.size()), col.layout.ncell);

    auto zi = col.zi.array();
    auto zt = col.zt.array();
    auto dz = col.dz.array();
    auto rho = col.rho.array();
    auto tk = col.tk.array();
    auto tkh_arr = col.tkh.array();
    auto thetal = col.thetal.array();
    auto theta = col.theta.array();
    auto exner = col.exner.array();
    auto p_mid = col.p_mid.array();
    auto tabs = col.tabs.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();
    auto u = col.u.array();
    auto v = col.v.array();
    auto tke = col.tke.array();
    auto host_dse = col.host_dse.array();
    auto sens = col.surf_sens_flux.array();
    auto lat = col.surf_lat_flux.array();
    auto tauu = col.surf_tau_u.array();
    auto tauv = col.surf_tau_v.array();

    for (int ic = 0; ic < col.layout.ncell; ++ic) {
        sens(ic,0,0) = wthl_sfc[ic];
        lat(ic,0,0) = wqw_sfc[ic];
        tauu(ic,0,0) = uw_sfc[ic];
        tauv(ic,0,0) = vw_sfc[ic];
        for (int k = 0; k <= col.layout.nlev; ++k) {
            zi(ic,k,0) = zi_grid[k];
        }
        for (int k = 0; k < col.layout.nlev; ++k) {
            zt(ic,k,0) = zt_grid[k];
            dz(ic,k,0) = dz_zt[k];
            rho(ic,k,0) = rho_zt[k];
            tk(ic,k,0) = tkh[k];
            tkh_arr(ic,k,0) = tkh[k];
            thetal(ic,k,0) = thetal_in[k];
            theta(ic,k,0) = thetal_in[k];
            exner(ic,k,0) = 1.0;
            p_mid(ic,k,0) = 100000.0 - 5000.0 * k;
            tabs(ic,k,0) = thetal_in[k];
            qv(ic,k,0) = qw_in[k];
            qc(ic,k,0) = 0.0;
            qi(ic,k,0) = 0.0;
            qw(ic,k,0) = qw_in[k];
            u(ic,k,0) = u_in[k];
            v(ic,k,0) = v_in[k];
            tke(ic,k,0) = tke_in[k];
            host_dse(ic,k,0) = Cp_d * tabs(ic,k,0) + CONST_GRAV * zt(ic,k,0);
            col.tke_tend.array()(ic,k,0) = 0.0;
        }
    }
}
}

TEST(ShocImplicit, TmpiAndDpInverseMatchColumnGeometry)
{
    auto col = shoc_test::make_column(4);
    amrex::Vector<amrex::Real> rho_zi(col.layout.nlev + 1, 0.0);
    for (int k = 0; k <= col.layout.nlev; ++k) {
        rho_zi[k] = 1.0 + 0.1 * k;
    }

    amrex::Vector<amrex::Real> tmpi;
    amrex::Vector<amrex::Real> rdp_zt;
    ShocImplicit::compute_tmpi(col, 0, 10.0, rho_zi, tmpi);
    ShocImplicit::compute_dp_inverse(col, 0, rdp_zt);

    ASSERT_EQ(tmpi.size(), col.layout.nlev + 1);
    ASSERT_EQ(rdp_zt.size(), col.layout.nlev);
    for (int k = 0; k <= col.layout.nlev; ++k) {
        EXPECT_GT(tmpi[k], 0.0);
    }
    for (int k = 0; k < col.layout.nlev; ++k) {
        const amrex::Real expected = 1.0 / (CONST_GRAV * col.rho.const_array()(0,k,0) * col.dz.const_array()(0,k,0));
        EXPECT_NEAR(rdp_zt[k], expected, 1.0e-14);
    }
}

TEST(ShocImplicit, ThermodynamicHelpersMatchTranslatedE3smFixtures)
{
    const auto tabs_identity =
        shoc_test::read_fixture_vector("implicit_energy/e3sm_compute_shoc_temperature_identity.txt");
    ASSERT_EQ(tabs_identity.size(), 3);
    EXPECT_NEAR(ShocImplicit::compute_temperature(300.0, 0.0, 1.0), tabs_identity[0], 1.0e-15);
    EXPECT_NEAR(ShocImplicit::compute_temperature(290.0, 0.0, 1.0), tabs_identity[1], 1.0e-15);
    EXPECT_NEAR(ShocImplicit::compute_temperature(280.0, 0.0, 1.0), tabs_identity[2], 1.0e-15);

    const auto tabs_profile =
        shoc_test::read_fixture_vector("implicit_energy/e3sm_compute_shoc_temperature_decreasing_profile.txt");
    ASSERT_EQ(tabs_profile.size(), 3);
    EXPECT_NEAR(ShocImplicit::compute_temperature(300.0, 1.0e-5, 1.0 / 1.1), tabs_profile[0], 1.0e-12);
    EXPECT_NEAR(ShocImplicit::compute_temperature(350.0, 1.0e-5, 1.0 / 1.5), tabs_profile[1], 1.0e-12);
    EXPECT_NEAR(ShocImplicit::compute_temperature(400.0, 1.0e-5, 0.5), tabs_profile[2], 1.0e-12);

    const amrex::Real exner = 0.75;
    const amrex::Real ql = 2.0e-2;
    const amrex::Real expected_nonzero_ql = 300.0 * exner + (L_v / Cp_d) * ql;
    EXPECT_NEAR(ShocImplicit::compute_temperature(300.0, ql, exner),
                expected_nonzero_ql, 1.0e-12)
        << "E3SM computes tabs = thetal*exner + (L_v/Cp)*ql; "
        << "the latent term must not be multiplied by exner.";

    const auto vapor_fixture =
        shoc_test::read_fixture_vector("implicit_energy/e3sm_compute_shoc_vapor.txt");
    ASSERT_EQ(vapor_fixture.size(), 5);
    EXPECT_NEAR(ShocImplicit::compute_vapor(1.0e-2, 0.0), vapor_fixture[0], 1.0e-15);
    EXPECT_NEAR(ShocImplicit::compute_vapor(1.2e-2, 0.0), vapor_fixture[1], 1.0e-15);
    EXPECT_NEAR(ShocImplicit::compute_vapor(1.5e-2, 1.5e-4), vapor_fixture[2], 1.0e-15);
    EXPECT_NEAR(ShocImplicit::compute_vapor(1.7e-2, 2.0e-3), vapor_fixture[3], 1.0e-15);
    EXPECT_NEAR(ShocImplicit::compute_vapor(2.0e-2, 0.0), vapor_fixture[4], 1.0e-15);
}

TEST(ShocImplicit, GeometryHelpersMatchTranslatedE3smFixtures)
{
    auto col = shoc_test::make_column(5);

    auto zi = col.zi.array();
    auto zt = col.zt.array();
    auto dz = col.dz.array();
    auto rho = col.rho.array();

    zi(0,0,0) = 0.0;
    zt(0,0,0) = 10.0;
    zi(0,1,0) = 35.0;
    zt(0,1,0) = 60.0;
    zi(0,2,0) = 110.0;
    zt(0,2,0) = 160.0;
    zi(0,3,0) = 260.0;
    zt(0,3,0) = 360.0;
    zi(0,4,0) = 610.0;
    zt(0,4,0) = 860.0;
    zi(0,5,0) = 1110.0;

    dz(0,0,0) = 10.0;
    dz(0,1,0) = 50.0;
    dz(0,2,0) = 100.0;
    dz(0,3,0) = 200.0;
    dz(0,4,0) = 500.0;

    rho(0,0,0) = 1.2;
    rho(0,1,0) = 1.0;
    rho(0,2,0) = 0.9;
    rho(0,3,0) = 0.8;
    rho(0,4,0) = 0.6;

    amrex::Vector<amrex::Real> rho_zi = {1.2, 1.0, 0.9, 0.8, 0.6, 0.5};
    amrex::Vector<amrex::Real> tmpi;
    amrex::Vector<amrex::Real> rdp_zt;

    ShocImplicit::compute_tmpi(col, 0, 300.0, rho_zi, tmpi);
    ShocImplicit::compute_dp_inverse(col, 0, rdp_zt);

    const auto tmpi_fixture =
        shoc_test::read_fixture_vector("implicit_energy/e3sm_compute_tmpi_erf_bottom_up_interior.txt");
    const auto rdp_fixture =
        shoc_test::read_fixture_vector("implicit_energy/e3sm_dp_inverse_erf_bottom_up.txt");

    ASSERT_EQ(tmpi_fixture.size(), 5);
    ASSERT_EQ(rdp_fixture.size(), 5);
    ASSERT_EQ(rdp_zt.size(), rdp_fixture.size());

    for (int k = 0; k < static_cast<int>(rdp_fixture.size()); ++k) {
        EXPECT_NEAR(rdp_zt[k], rdp_fixture[k], 1.0e-15);
        EXPECT_NEAR(tmpi[k], tmpi_fixture[k], 1.0e-12);
    }
}

TEST(ShocImplicit, UniformProfileWithNoFluxRemainsUnchanged)
{
    auto col = shoc_test::make_column(5);
    ShocRuntimeOptions opts;

    auto thetal = col.thetal.array();
    auto theta = col.theta.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qw = col.qw.array();
    auto tke = col.tke.array();
    auto u = col.u.array();
    auto v = col.v.array();
    auto tk = col.tk.array();
    auto tkh = col.tkh.array();
    auto exner = col.exner.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        thetal(0,k,0) = 300.0;
        theta(0,k,0) = 300.0;
        qv(0,k,0) = 0.008;
        qc(0,k,0) = 0.0;
        qw(0,k,0) = 0.008;
        tke(0,k,0) = 0.35;
        u(0,k,0) = 4.0;
        v(0,k,0) = -1.5;
        tk(0,k,0) = 1.0;
        tkh(0,k,0) = 1.0;
        exner(0,k,0) = 1.0;
        col.tke_tend.array()(0,k,0) = 0.0;
    }

    shoc::set_fab_val(col.surf_sens_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.0, shoc::InitRunOn::Host);

    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(col, opts, 10.0);
    });

    const auto theta_tend = col.theta_tend.const_array();
    const auto qv_tend = col.qv_tend.const_array();
    const auto qc_tend = col.qc_tend.const_array();
    const auto u_tend = col.u_tend.const_array();
    const auto v_tend = col.v_tend.const_array();
    const auto tke_tend = col.tke_tend.const_array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_NEAR(theta_tend(0,k,0), 0.0, 1.0e-12);
        EXPECT_NEAR(qv_tend(0,k,0), 0.0, 1.0e-12);
        EXPECT_NEAR(qc_tend(0,k,0), 0.0, 1.0e-12);
        EXPECT_NEAR(u_tend(0,k,0), 0.0, 1.0e-12);
        EXPECT_NEAR(v_tend(0,k,0), 0.0, 1.0e-12);
        EXPECT_NEAR(tke_tend(0,k,0), 0.0, 1.0e-12);
    }
}

TEST(ShocImplicit, PdfDiagnosedLiquidFeedsHostWriteback)
{
    auto col = shoc_test::make_column(4);
    ShocRuntimeOptions opts;

    auto thetal = col.thetal.array();
    auto theta = col.theta.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();
    auto tke = col.tke.array();
    auto u = col.u.array();
    auto v = col.v.array();
    auto tk = col.tk.array();
    auto tkh = col.tkh.array();
    auto exner = col.exner.array();
    auto shoc_ql = col.shoc_ql.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        thetal(0,k,0) = 300.0;
        theta(0,k,0) = 300.0;
        qv(0,k,0) = 0.012;
        qc(0,k,0) = 0.0;
        qi(0,k,0) = 0.0;
        qw(0,k,0) = 0.012;
        tke(0,k,0) = 0.3;
        u(0,k,0) = 3.0;
        v(0,k,0) = -1.0;
        tk(0,k,0) = 0.5;
        tkh(0,k,0) = 0.5;
        exner(0,k,0) = 1.0;
        shoc_ql(0,k,0) = (k == 1) ? 1.0e-3 : 0.0;
        col.tke_tend.array()(0,k,0) = 0.0;
    }

    shoc::set_fab_val(col.surf_sens_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.0, shoc::InitRunOn::Host);

    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(col, opts, 10.0);
    });

    const auto theta_tend = col.theta_tend.const_array();
    const auto qv_tend = col.qv_tend.const_array();
    const auto qc_tend = col.qc_tend.const_array();
    const auto theta_out = col.theta.const_array();
    const auto qv_out = col.qv.const_array();
    const auto qc_out = col.qc.const_array();

    EXPECT_NEAR(qc_out(0,1,0), 1.0e-3, 1.0e-12);
    EXPECT_NEAR(qv_out(0,1,0), 0.011, 1.0e-12);
    EXPECT_GT(theta_out(0,1,0), 300.0);
    EXPECT_GT(qc_tend(0,1,0), 0.0);
    EXPECT_LT(qv_tend(0,1,0), 0.0);
    EXPECT_GT(theta_tend(0,1,0), 0.0);

    for (int k = 0; k < col.layout.nlev; ++k) {
        if (k == 1) continue;
        EXPECT_NEAR(qc_out(0,k,0), 0.0, 1.0e-12);
    }
}

TEST(ShocImplicit, FinalWritebackDoesNotCreateLiquidBeyondPdfDiagnosis)
{
    auto col = shoc_test::make_column(3);
    ShocRuntimeOptions opts;

    auto thetal = col.thetal.array();
    auto theta = col.theta.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();
    auto tk = col.tk.array();
    auto tkh = col.tkh.array();
    auto exner = col.exner.array();
    auto shoc_ql = col.shoc_ql.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        thetal(0,k,0) = 286.0;
        theta(0,k,0) = 286.0;
        qv(0,k,0) = 0.020;
        qc(0,k,0) = 0.0;
        qi(0,k,0) = 0.0;
        qw(0,k,0) = 0.020;
        tk(0,k,0) = 0.0;
        tkh(0,k,0) = 0.0;
        exner(0,k,0) = 1.0;
        shoc_ql(0,k,0) = 0.0;
        col.tke_tend.array()(0,k,0) = 0.0;
    }

    shoc::set_fab_val(col.surf_sens_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.0, shoc::InitRunOn::Host);

    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(col, opts, 1.0);
    });

    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_NEAR(col.shoc_ql.const_array()(0,k,0), 0.0, 1.0e-12);
        EXPECT_NEAR(col.qc.const_array()(0,k,0), 0.0, 1.0e-12);
        EXPECT_NEAR(col.qi.const_array()(0,k,0), 0.0, 1.0e-12);
        EXPECT_NEAR(col.qv.const_array()(0,k,0), col.qw.const_array()(0,k,0), 1.0e-12);
    }
}

TEST(ShocImplicit, SurfaceFluxesDriveBottomCellAndKeepMoistureBounded)
{
    auto col = shoc_test::make_column(6);
    ShocRuntimeOptions opts;

    shoc::set_fab_val(col.surf_sens_flux, 0.03, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 2.0e-4, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, -0.04, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.01, shoc::InitRunOn::Host);

    auto tk = col.tk.array();
    auto tkh = col.tkh.array();
    auto exner = col.exner.array();
    for (int k = 0; k < col.layout.nlev; ++k) {
        tk(0,k,0) = 2.0;
        tkh(0,k,0) = 1.5;
        exner(0,k,0) = 1.0;
        col.tke_tend.array()(0,k,0) = 0.0;
    }

    const amrex::Real theta_before = col.thetal.const_array()(0,0,0);
    amrex::Real total_qw_before = 0.0;
    for (int k = 0; k < col.layout.nlev; ++k) {
        total_qw_before += col.qw.const_array()(0,k,0) *
                           col.rho.const_array()(0,k,0) *
                           col.dz.const_array()(0,k,0);
    }

    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(col, opts, 10.0);
    });

    const auto thetal = col.thetal.const_array();
    const auto qw = col.qw.const_array();
    const auto qv = col.qv.const_array();
    const auto qc = col.qc.const_array();
    const auto u = col.u.const_array();
    const auto v = col.v.const_array();
    const auto tke = col.tke.const_array();

    EXPECT_GT(thetal(0,0,0), theta_before);
    amrex::Real total_qw_after = 0.0;
    for (int k = 0; k < col.layout.nlev; ++k) {
        total_qw_after += qw(0,k,0) *
                          col.rho.const_array()(0,k,0) *
                          col.dz.const_array()(0,k,0);
    }
    EXPECT_GT(total_qw_after, total_qw_before);
    EXPECT_TRUE(std::isfinite(u(0,0,0)));
    EXPECT_TRUE(std::isfinite(v(0,0,0)));
    EXPECT_TRUE(std::isfinite(tke(0,0,0)));
    EXPECT_TRUE(std::isfinite(col.theta_tend.const_array()(0,0,0)));

    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_GE(qw(0,k,0), 0.0);
        EXPECT_GE(qv(0,k,0), 0.0);
        EXPECT_GE(qc(0,k,0), 0.0);
        EXPECT_NEAR(qw(0,k,0), qv(0,k,0) + qc(0,k,0), 1.0e-12);
        EXPECT_GE(tke(0,k,0), 4.0e-4);
        EXPECT_LE(tke(0,k,0), 50.0);
    }
}

TEST(ShocImplicit, CloudLiquidRaisesThetaAboveThetal)
{
    auto col = shoc_test::make_column(4);
    ShocRuntimeOptions opts;

    auto thetal = col.thetal.array();
    auto theta = col.theta.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qw = col.qw.array();
    auto shoc_ql = col.shoc_ql.array();
    auto tk = col.tk.array();
    auto tkh = col.tkh.array();
    auto exner = col.exner.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        thetal(0,k,0) = 299.0;
        theta(0,k,0) = 299.0;
        qv(0,k,0) = (k == 0) ? 0.03 : 0.011;
        qc(0,k,0) = (k == 0) ? 2.0e-3 : 0.0;
        qw(0,k,0) = qv(0,k,0) + qc(0,k,0);
        shoc_ql(0,k,0) = qc(0,k,0);
        tk(0,k,0) = 0.1;
        tkh(0,k,0) = 0.1;
        exner(0,k,0) = 1.0;
        col.tke_tend.array()(0,k,0) = 0.0;
    }

    shoc::set_fab_val(col.surf_sens_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.0, shoc::InitRunOn::Host);

    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(col, opts, 1.0);
    });

    EXPECT_GT(col.shoc_ql.const_array()(0,0,0), 0.0);
    EXPECT_GE(col.theta.const_array()(0,0,0), col.thetal.const_array()(0,0,0));
}

TEST(ShocImplicit, NegativeMoistureFluxIsClippedAtZeroWater)
{
    auto col = shoc_test::make_column(4);
    ShocRuntimeOptions opts;

    auto qw = col.qw.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto tk = col.tk.array();
    auto tkh = col.tkh.array();
    auto exner = col.exner.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        qw(0,k,0) = 1.0e-6;
        qv(0,k,0) = 1.0e-6;
        qc(0,k,0) = 0.0;
        tk(0,k,0) = 0.2;
        tkh(0,k,0) = 0.2;
        exner(0,k,0) = 1.0;
        col.tke_tend.array()(0,k,0) = 0.0;
    }

    shoc::set_fab_val(col.surf_lat_flux, -5.0e-4, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_sens_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.0, shoc::InitRunOn::Host);

    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(col, opts, 20.0);
    });

    const auto qv_new = col.qv.const_array();
    const auto qc_new = col.qc.const_array();
    const auto qw_new = col.qw.const_array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_GE(qw_new(0,k,0), 0.0);
        EXPECT_GE(qv_new(0,k,0), 0.0);
        EXPECT_GE(qc_new(0,k,0), 0.0);
        EXPECT_NEAR(qw_new(0,k,0), qv_new(0,k,0) + qc_new(0,k,0), 1.0e-12);
    }
}

TEST(ShocEnergyFixer, ActiveTopFollowsHighestTurbulentLevel)
{
    amrex::Vector<amrex::Real> tke = {4.0e-4, 4.0e-4, 6.0e-4, 5.0e-4, 4.0e-4};
    EXPECT_EQ(ShocEnergyFixer::diagnose_active_top(tke), 3);

    tke = {4.0e-4, 4.0e-4, 4.0e-4};
    EXPECT_EQ(ShocEnergyFixer::diagnose_active_top(tke), -1);
}

TEST(ShocEnergyFixer, LiquidPartitionChangeDoesNotCreateEnergyCorrection)
{
    auto col = shoc_test::make_column(3);
    auto rho = col.rho.array();
    auto dz = col.dz.array();
    auto zt = col.zt.array();
    auto exner = col.exner.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        rho(0,k,0) = 1.0;
        dz(0,k,0) = 100.0;
        zt(0,k,0) = 50.0 + 100.0 * k;
        exner(0,k,0) = 0.8;
    }
    shoc::set_fab_val(col.surf_sens_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 0.0, shoc::InitRunOn::Host);

    amrex::Vector<amrex::Real> thl_old(3, 300.0);
    amrex::Vector<amrex::Real> qv_old(3, 1.0e-2);
    amrex::Vector<amrex::Real> qc_old(3, 0.0);
    amrex::Vector<amrex::Real> qi_old(3, 0.0);
    amrex::Vector<amrex::Real> u_old(3, 2.0);
    amrex::Vector<amrex::Real> v_old(3, -1.0);
    amrex::Vector<amrex::Real> tke_old(3, 0.2);

    amrex::Vector<amrex::Real> thl_new = thl_old;
    amrex::Vector<amrex::Real> qv_new(3, 8.0e-3);
    amrex::Vector<amrex::Real> qc_new(3, 2.0e-3);
    amrex::Vector<amrex::Real> qi_new(3, 0.0);
    amrex::Vector<amrex::Real> u_new = u_old;
    amrex::Vector<amrex::Real> v_new = v_old;
    amrex::Vector<amrex::Real> tke_new = tke_old;

    ShocEnergyFixer::apply_column(col, 0, 10.0,
                                  thl_old, qv_old, qc_old, qi_old,
                                  u_old, v_old, tke_old,
                                  thl_new, qv_new, qc_new, qi_new,
                                  u_new, v_new, tke_new);

    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_NEAR(thl_new[k], thl_old[k], 1.0e-12)
            << "E3SM's moist-energy convention cancels the latent heating "
            << "already present in host_dse for liquid repartitioning.";
    }
}

TEST(ShocImplicit, TranslatedE3smMultiColumnPropertyCaseStaysPhysical)
{
    ShocColumnData col;
    ShocColumnLayout layout;
    layout.nx = 5;
    layout.ny = 1;
    layout.ncell = 5;
    layout.nlev = 5;
    define_shoc_column_data(col, layout, shoc_test::test_arena(), shoc::InitRunOn::Host);
    shoc_test::sync();
    ShocRuntimeOptions opts;
    const auto fixture = shoc_test::read_named_fixture_vectors(
        "implicit_energy/e3sm_update_prognostics_implicit_multicolumn.txt");
    const auto& wthl_sfc = shoc_test::fixture_values(fixture, "wthl_sfc");
    const auto& wqw_sfc = shoc_test::fixture_values(fixture, "wqw_sfc");

    load_translated_e3sm_implicit_multicolumn_fixture(col);

    const auto thetal = col.thetal.const_array();
    const auto qw = col.qw.const_array();
    const auto rho = col.rho.const_array();
    const auto dz = col.dz.const_array();

    amrex::Vector<amrex::Real> thl_before(layout.ncell, 0.0);
    amrex::Vector<amrex::Real> qw_before(layout.ncell, 0.0);
    for (int ic = 0; ic < layout.ncell; ++ic) {
        for (int k = 0; k < layout.nlev; ++k) {
            thl_before[ic] += thetal(ic,k,0) * rho(ic,k,0) * dz(ic,k,0);
            qw_before[ic] += qw(ic,k,0) * rho(ic,k,0) * dz(ic,k,0);
        }
    }

    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(col, opts, 300.0);
    });

    const auto thetal_new = col.thetal.const_array();
    const auto qw_new = col.qw.const_array();
    const auto qv_new = col.qv.const_array();
    const auto qc_new = col.qc.const_array();
    const auto qi_new = col.qi.const_array();
    const auto u_new = col.u.const_array();
    const auto v_new = col.v.const_array();
    const auto tke_new = col.tke.const_array();

    for (int ic = 0; ic < layout.ncell; ++ic) {
        amrex::Real thl_after = 0.0;
        amrex::Real qw_after = 0.0;
        for (int k = 0; k < layout.nlev; ++k) {
            EXPECT_TRUE(std::isfinite(thetal_new(ic,k,0)));
            EXPECT_TRUE(std::isfinite(qw_new(ic,k,0)));
            EXPECT_TRUE(std::isfinite(u_new(ic,k,0)));
            EXPECT_TRUE(std::isfinite(v_new(ic,k,0)));
            EXPECT_TRUE(std::isfinite(tke_new(ic,k,0)));
            EXPECT_GE(qw_new(ic,k,0), 0.0);
            EXPECT_GE(qv_new(ic,k,0), 0.0);
            EXPECT_GE(qc_new(ic,k,0), 0.0);
            EXPECT_GE(qi_new(ic,k,0), 0.0);
            EXPECT_NEAR(qw_new(ic,k,0), qv_new(ic,k,0) + qc_new(ic,k,0) + qi_new(ic,k,0), 1.0e-12);
            EXPECT_GE(tke_new(ic,k,0), 4.0e-4);
            EXPECT_LE(tke_new(ic,k,0), 50.0);
            thl_after += thetal_new(ic,k,0) * rho(ic,k,0) * dz(ic,k,0);
            qw_after += qw_new(ic,k,0) * rho(ic,k,0) * dz(ic,k,0);
        }

        if (wthl_sfc[ic] > 0.0) {
            EXPECT_GT(thl_after, thl_before[ic]);
        }
        if (wqw_sfc[ic] > 0.0) {
            EXPECT_GT(qw_after, qw_before[ic]);
        }
        if (wqw_sfc[ic] < 0.0) {
            EXPECT_LT(qw_after, qw_before[ic]);
        }
    }
}

// Verify the surface flux scaling factor (cmnfac) matches the E3SM formula:
//   cmnfac = dt * g * rho_zi_surface * rdp_zt_bottom
// NOT tmpi[0]*rdp_zt[0], which includes an extra 1/dz_interface factor
// that belongs only in the tridiagonal diffusion matrix.
TEST(ShocImplicit, SurfaceFluxScalingMatchesE3smFormula)
{
    auto col = shoc_test::make_column(5);

    // Set up a uniform 40 m grid (BOMEX-like)
    auto zi = col.zi.array();
    auto zt = col.zt.array();
    auto dz = col.dz.array();
    auto rho = col.rho.array();
    auto tk = col.tk.array();
    auto tkh = col.tkh.array();
    auto exner = col.exner.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        zi(0,k,0) = 40.0 * k;
        zt(0,k,0) = 40.0 * k + 20.0;
        dz(0,k,0) = 40.0;
        rho(0,k,0) = 1.15;
        tk(0,k,0) = 0.0;   // zero diffusion: only surface flux matters
        tkh(0,k,0) = 0.0;
        exner(0,k,0) = 1.0;
        col.tke_tend.array()(0,k,0) = 0.0;
        col.thetal.array()(0,k,0) = 300.0;
        col.theta.array()(0,k,0) = 300.0;
        col.qw.array()(0,k,0) = 0.01;
        col.qv.array()(0,k,0) = 0.01;
        col.qc.array()(0,k,0) = 0.0;
        col.qi.array()(0,k,0) = 0.0;
        col.tke.array()(0,k,0) = 0.3;
        col.u.array()(0,k,0) = 4.0;
        col.v.array()(0,k,0) = -1.0;
    }
    zi(0,col.layout.nlev,0) = 40.0 * col.layout.nlev;

    const amrex::Real wthl_sfc = 0.008;
    shoc::set_fab_val(col.surf_sens_flux, wthl_sfc, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.0, shoc::InitRunOn::Host);

    const amrex::Real dt = 5.0;
    ShocRuntimeOptions opts;
    shoc_test::run_and_sync([&] {
        ShocImplicit::update_prognostics(col, opts, dt);
    });

    // With zero diffusion and no momentum flux, only the bottom cell
    // should change, and only in thetal.  The E3SM scaling gives:
    //   cmnfac = dt * g * rho_zi / (g * rho * dz) = dt * rho_zi / (rho * dz)
    // where rho_zi at the bottom ≈ rho (uniform column) and dz = 40 m.
    //   delta_thetal = cmnfac * wthl_sfc = dt / dz * wthl_sfc
    //               = 5.0 / 40.0 * 0.008 = 0.001 K
    const amrex::Real rho_zi_sfc = rho(0,0,0); // linear interp at z=0 ≈ rho
    const amrex::Real expected_cmnfac = dt * rho_zi_sfc / (rho(0,0,0) * 40.0);
    const amrex::Real expected_delta = expected_cmnfac * wthl_sfc;

    const auto thetal_new = col.thetal.const_array();
    const amrex::Real actual_delta = thetal_new(0,0,0) - 300.0;

    // The key check: the bottom-cell increment must match the E3SM formula,
    // not be ~20x weaker (which the old tmpi[0]*rdp_zt[0] factor produced).
    EXPECT_NEAR(actual_delta, expected_delta, 1.0e-6)
        << "Surface flux scaling diverges from E3SM cmnfac = dt*g*rho_zi*rdp_zt";

    // Cells above the bottom must be essentially unchanged (zero diffusion)
    for (int k = 1; k < col.layout.nlev; ++k) {
        EXPECT_NEAR(thetal_new(0,k,0), 300.0, 1.0e-12);
    }
}

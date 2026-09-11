#include "ERF_ShocStructure.H"
#include "ERF_ShocThermoUtils.H"
#include "ERF_ShocTestUtils.H"
#include "ERF_ShocTypes.H"

#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <string>

TEST(ShocRuntimeOptions, DefaultsValidate)
{
    ShocRuntimeOptions opts;
    EXPECT_NO_FATAL_FAILURE(validate_shoc_runtime_options(opts));
    EXPECT_GT(opts.length_fac, 0.0);
    EXPECT_GE(opts.coeff_km, 0.0);
    EXPECT_GE(opts.coeff_kh, 0.0);
    EXPECT_FALSE(opts.debug_summary);
    EXPECT_EQ(opts.transport_mode, ShocTransportMode::StateUpdate);
    EXPECT_EQ(opts.momentum_transport, ShocMomentumTransport::StateUpdate);
}

TEST(ShocRuntimeOptions, TransportModeHelpersMatchIntent)
{
    EXPECT_TRUE(shoc_uses_momentum_state_update(ShocMomentumTransport::StateUpdate));
    EXPECT_FALSE(shoc_uses_momentum_state_update(ShocMomentumTransport::None));
    EXPECT_TRUE(shoc_disables_momentum_transport(ShocMomentumTransport::None));
    EXPECT_FALSE(shoc_disables_momentum_transport(ShocMomentumTransport::StateUpdate));
    EXPECT_EQ(std::string(shoc_transport_mode_name(ShocTransportMode::StateUpdate)), "state_update");
    EXPECT_EQ(std::string(shoc_momentum_transport_name(ShocMomentumTransport::StateUpdate)), "state_update");
    EXPECT_EQ(std::string(shoc_momentum_transport_name(ShocMomentumTransport::None)), "none");
}

namespace
{
class ScopedParmParseString
{
public:
    ScopedParmParseString (const char* name, const std::string& value)
        : m_pp("erf.shoc"),
          m_name(name)
    {
        m_had_previous = m_pp.query(m_name, m_previous);
        m_pp.remove(m_name);
        m_pp.add(m_name, value);
    }

    ~ScopedParmParseString ()
    {
        m_pp.remove(m_name);
        if (m_had_previous) {
            m_pp.add(m_name, m_previous);
        }
    }

private:
    amrex::ParmParse m_pp;
    std::string m_name;
    std::string m_previous;
    bool m_had_previous = false;
};

class ScopedParmParseRemoval
{
public:
    explicit ScopedParmParseRemoval (const char* name)
        : m_pp("erf.shoc"),
          m_name(name)
    {
        m_had_previous = m_pp.query(m_name, m_previous);
        m_pp.remove(m_name);
    }

    ~ScopedParmParseRemoval ()
    {
        if (m_had_previous) {
            m_pp.add(m_name, m_previous);
        }
    }

private:
    amrex::ParmParse m_pp;
    std::string m_name;
    std::string m_previous;
    bool m_had_previous = false;
};

void
shift_column_heights (ShocColumnData& col, amrex::Real offset)
{
    auto zt = col.zt.array();
    auto zi = col.zi.array();
    for (int k = 0; k < col.layout.nlev; ++k) {
        zt(0,k,0) += offset;
    }
    for (int k = 0; k <= col.layout.nlev; ++k) {
        zi(0,k,0) += offset;
    }
    shoc_test::sync();
}

ShocColumnData
make_first_level_ri_crossing_column (amrex::Real surface_sensible_flux)
{
    auto col = shoc_test::make_column(2);
    auto thetal = col.thetal.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();
    auto u = col.u.array();
    auto v = col.v.array();

    // Nominally target a Richardson number well above the critical value;
    // the expected crossing below is derived from the represented state.
    constexpr amrex::Real theta_increment = amrex::Real(0.002752293577981651);
    for (int k = 0; k < col.layout.nlev; ++k) {
        thetal(0,k,0) = amrex::Real(300.0) + theta_increment * k;
        qv(0,k,0) = amrex::Real(0.0);
        qc(0,k,0) = amrex::Real(0.0);
        qi(0,k,0) = amrex::Real(0.0);
        qw(0,k,0) = amrex::Real(0.0);
        u(0,k,0) = amrex::Real(2.0);
        v(0,k,0) = amrex::Real(1.0);
    }

    shoc::set_fab_val(col.surf_sens_flux, surface_sensible_flux, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, amrex::Real(0.0), shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, amrex::Real(0.0), shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, amrex::Real(0.0), shoc::InitRunOn::Host);
    shoc_test::sync();

    return col;
}
}

TEST(ShocRuntimeOptions, LegacyTendenciesTransportModeIsRejected)
{
    ShocRuntimeOptions opts;
    std::string error_message;
    EXPECT_FALSE(parse_shoc_transport_mode_string("tendencies", opts.transport_mode, error_message));
    EXPECT_NE(error_message.find("removed for native SHOC"), std::string::npos);
}

TEST(ShocRuntimeOptions, RemovedScalarHostDiffusionModeIsRejected)
{
    ShocRuntimeOptions opts;
    std::string error_message;
    EXPECT_FALSE(parse_shoc_transport_mode_string("HOST_DIFFUSION", opts.transport_mode,
                                                  error_message));
    EXPECT_NE(error_message.find("has been removed for native SHOC"), std::string::npos);
    EXPECT_NE(error_message.find("Use erf.shoc.transport_mode = state_update"),
              std::string::npos);
}

TEST(ShocRuntimeOptions, RemovedMomentumHostDiffusionModeIsRejected)
{
    ShocMomentumTransport transport = ShocMomentumTransport::StateUpdate;
    std::string error_message;
    EXPECT_FALSE(parse_shoc_momentum_transport_string("HoSt_DiFfUsIoN", transport,
                                                      error_message));
    EXPECT_NE(error_message.find("has been removed for native SHOC"), std::string::npos);
    EXPECT_NE(error_message.find("'state_update'"), std::string::npos);
    EXPECT_NE(error_message.find("'none'"), std::string::npos);
}

TEST(ShocRuntimeOptions, InvalidMomentumTransportIsRejected)
{
    ShocMomentumTransport transport = ShocMomentumTransport::None;
    std::string error_message;
    EXPECT_FALSE(parse_shoc_momentum_transport_string("tendons", transport, error_message));
    EXPECT_NE(error_message.find("erf.shoc.momentum_transport"), std::string::npos);
    EXPECT_NE(error_message.find("'none'"), std::string::npos);
    EXPECT_NE(error_message.find("'state_update'"), std::string::npos);
    EXPECT_EQ(error_message.find("host_diffusion"), std::string::npos);
}

TEST(ShocRuntimeOptions, InvalidScalarTransportAdvertisesOnlyStateUpdate)
{
    ShocRuntimeOptions opts;
    std::string error_message;
    EXPECT_FALSE(parse_shoc_transport_mode_string("unknown", opts.transport_mode, error_message));
    EXPECT_NE(error_message.find("'state_update'"), std::string::npos);
    EXPECT_EQ(error_message.find("host_diffusion"), std::string::npos);
}

TEST(ShocRuntimeOptions, SharedReaderUsesDefaultsAndAcceptsSupportedModes)
{
    ScopedParmParseRemoval remove_transport("transport_mode");
    ScopedParmParseRemoval remove_momentum("momentum_transport");

    ShocRuntimeOptions defaults;
    read_shoc_transport_modes(defaults.transport_mode, defaults.momentum_transport);
    EXPECT_EQ(defaults.transport_mode, ShocTransportMode::StateUpdate);
    EXPECT_EQ(defaults.momentum_transport, ShocMomentumTransport::StateUpdate);

    ScopedParmParseString transport("transport_mode", "state_update");
    ScopedParmParseString momentum("momentum_transport", "none");
    ShocTransportMode transport_mode = ShocTransportMode::StateUpdate;
    ShocMomentumTransport momentum_mode = ShocMomentumTransport::StateUpdate;
    read_shoc_transport_modes(transport_mode, momentum_mode);
    EXPECT_EQ(transport_mode, ShocTransportMode::StateUpdate);
    EXPECT_EQ(momentum_mode, ShocMomentumTransport::None);
}

TEST(SolverChoice, NativeStateUpdateOwnsNoVerticalDiffusion)
{
    SolverChoice choice;
    choice.use_native_shoc = true;
    choice.use_eamxx_shoc = false;
    choice.shoc_transport_mode = ShocTransportMode::StateUpdate;
    choice.shoc_momentum_transport = ShocMomentumTransport::StateUpdate;

    EXPECT_FALSE(choice.host_owns_vertical_scalar_diffusion());
    EXPECT_FALSE(choice.host_owns_vertical_momentum_diffusion());
}

TEST(SolverChoice, NativeMomentumNoneLeavesHostMomentumOwnership)
{
    SolverChoice choice;
    choice.use_native_shoc = true;
    choice.use_eamxx_shoc = false;
    choice.shoc_transport_mode = ShocTransportMode::StateUpdate;
    choice.shoc_momentum_transport = ShocMomentumTransport::None;

    EXPECT_FALSE(choice.host_owns_vertical_scalar_diffusion());
    EXPECT_TRUE(choice.host_owns_vertical_momentum_diffusion());
}

TEST(SolverChoice, NonNativeShocRetainsGenericHostOwnership)
{
    SolverChoice choice;
    choice.use_native_shoc = false;
    choice.use_eamxx_shoc = false;

    EXPECT_TRUE(choice.host_owns_vertical_scalar_diffusion());
    EXPECT_TRUE(choice.host_owns_vertical_momentum_diffusion());
}

TEST(SolverChoice, EamxxShocRetainsExistingOwnership)
{
    SolverChoice choice;
    choice.use_native_shoc = false;
    choice.use_eamxx_shoc = true;

    EXPECT_FALSE(choice.host_owns_vertical_scalar_diffusion());
    EXPECT_FALSE(choice.host_owns_vertical_momentum_diffusion());
}

TEST(ShocStructure, SurfaceLayerUsesUstarFloorAndFiniteObukhov)
{
    auto col = shoc_test::make_column(4);
    auto wthv = col.wthv_sec.array();
    for (int k = 1; k < col.layout.nlev; ++k) {
        wthv(0,k,0) = 1.0e-3 * k;
    }
    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
    });

    const auto ustar = col.ustar.const_array();
    const auto obklen = col.obklen.const_array();
    const auto wthv_out = col.wthv_sec.const_array();

    EXPECT_DOUBLE_EQ(ustar(0,0,0), 0.01);
    EXPECT_TRUE(std::isfinite(obklen(0,0,0)));
    EXPECT_DOUBLE_EQ(wthv_out(0,0,0), wthv_out(0,0,0));
    for (int k = 1; k < col.layout.nlev; ++k) {
        EXPECT_DOUBLE_EQ(wthv_out(0,k,0), 1.0e-3 * k);
    }
}

TEST(ShocStructure, SurfaceLayerObukhovUsesClampedUstar)
{
    auto col = shoc_test::make_column(4);
    auto thetal = col.thetal.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto exner = col.exner.array();

    thetal(0,0,0) = 300.0;
    qv(0,0,0) = 0.010;
    qc(0,0,0) = 0.0;
    qi(0,0,0) = 0.0;
    exner(0,0,0) = 1.0;

    shoc::set_fab_val(col.surf_tau_u, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_sens_flux, 0.02, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 0.0, shoc::InitRunOn::Host);

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
    });

    const amrex::Real thv_sfc = 300.0 * (1.0 + 0.61 * 0.010);
    const amrex::Real ustar = 0.01;
    const amrex::Real ustar_cu = ustar * ustar * ustar;
    const amrex::Real kbfs = 0.02;
    const amrex::Real expected_obk = -thv_sfc * ustar_cu /
                                     (CONST_GRAV * KAPPA * (kbfs + 1.0e-10));

    EXPECT_DOUBLE_EQ(col.ustar.const_array()(0,0,0), ustar);
    EXPECT_NEAR(col.obklen.const_array()(0,0,0), expected_obk, 1.0e-12);
}

TEST(ShocStructure, SurfaceLayerUsesShocThermodynamicMapping)
{
    auto col = shoc_test::make_column(4);
    auto thetal = col.thetal.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();

    thetal(0,0,0) = 298.0;
    qv(0,0,0) = 0.010;
    qc(0,0,0) = 1.0e-3;
    qi(0,0,0) = 2.0e-4;
    qw(0,0,0) = qv(0,0,0) + qc(0,0,0) + qi(0,0,0);
    for (int k = 1; k < col.layout.nlev; ++k) {
        col.wthv_sec.array()(0,k,0) = 5.0e-4 * k;
    }
    shoc::set_fab_val(col.surf_tau_u, 0.04, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.03, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_sens_flux, 0.02, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 2.0e-4, shoc::InitRunOn::Host);

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
    });

    const amrex::Real qc_sfc = 1.0e-3;
    const amrex::Real qi_sfc = 2.0e-4;
    const amrex::Real cldliq = qc_sfc + qi_sfc;
    const amrex::Real th_sfc = 298.0 + (L_v * qc_sfc + shoc::latent_sublimation() * qi_sfc) / Cp_d;
    const amrex::Real thv_sfc = th_sfc * (1.0 + 0.61 * 0.010 - cldliq);
    const amrex::Real ustar_raw = std::sqrt(std::sqrt(0.04 * 0.04 + 0.03 * 0.03));
    const amrex::Real ustar = amrex::max(amrex::Real(0.01), ustar_raw);
    const amrex::Real kbfs = 0.02 + 0.61 * th_sfc * 2.0e-4;
    const amrex::Real expected_obk = -thv_sfc * std::pow(ustar, 3) /
                                     (CONST_GRAV * KAPPA * (kbfs + 1.0e-10));

    EXPECT_NEAR(col.ustar.const_array()(0,0,0), ustar, 1.0e-12);
    EXPECT_NEAR(col.obklen.const_array()(0,0,0), expected_obk, 1.0e-10);
    EXPECT_NEAR(col.wthv_sec.const_array()(0,0,0), kbfs, 1.0e-12);
    for (int k = 1; k < col.layout.nlev; ++k) {
        EXPECT_DOUBLE_EQ(col.wthv_sec.const_array()(0,k,0), 5.0e-4 * k);
    }
}

TEST(ShocStructure, SurfaceLayerThermodynamicMappingUsesExner)
{
    auto col = shoc_test::make_column(4);
    auto thetal = col.thetal.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();
    auto exner = col.exner.array();

    thetal(0,0,0) = 298.0;
    qv(0,0,0) = 0.010;
    qc(0,0,0) = 1.0e-3;
    qi(0,0,0) = 2.0e-4;
    qw(0,0,0) = qv(0,0,0) + qc(0,0,0) + qi(0,0,0);
    exner(0,0,0) = 0.8;
    shoc::set_fab_val(col.surf_tau_u, 0.04, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.03, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_sens_flux, 0.02, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 2.0e-4, shoc::InitRunOn::Host);

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
    });

    const amrex::Real qc_sfc = 1.0e-3;
    const amrex::Real qi_sfc = 2.0e-4;
    const amrex::Real theta_sfc = 298.0 + (L_v * qc_sfc + shoc::latent_sublimation() * qi_sfc) /
                                             (Cp_d * 0.8);
    const amrex::Real kbfs = 0.02 + 0.61 * theta_sfc * 2.0e-4;

    EXPECT_NEAR(col.wthv_sec.const_array()(0,0,0), kbfs, 1.0e-12)
        << "SHOC theta_l to theta conversion must divide the latent term by ERF exner.";
}

TEST(ShocStructure, SurfaceLayerMatchesTranslatedE3smFixture)
{
    auto col = shoc_test::make_column(4);
    auto thetal = col.thetal.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();

    thetal(0,0,0) = 298.0;
    qv(0,0,0) = 0.010;
    qc(0,0,0) = 1.0e-3;
    qi(0,0,0) = 0.0;
    qw(0,0,0) = 0.011;
    shoc::set_fab_val(col.surf_tau_u, 0.04, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.03, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_sens_flux, 0.02, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 2.0e-4, shoc::InitRunOn::Host);

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
    });

    const auto fixture =
        shoc_test::read_fixture_vector("structure/e3sm_diag_obklen_surface_mapping.txt");
    ASSERT_EQ(fixture.size(), 3);

    EXPECT_NEAR(col.ustar.const_array()(0,0,0), fixture[0], 1.0e-15);
    EXPECT_NEAR(col.wthv_sec.const_array()(0,0,0), fixture[1], 1.0e-15);
    EXPECT_NEAR(col.obklen.const_array()(0,0,0), fixture[2], 1.0e-14);
}

TEST(ShocStructure, PblHeightUsesVaporNotTotalWaterInVirtualTheta)
{
    auto col = shoc_test::make_column(5);
    auto thetal = col.thetal.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();
    auto u = col.u.array();
    auto v = col.v.array();

    shoc::set_fab_val(col.surf_tau_u, 0.04, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.03, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_sens_flux, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 0.0, shoc::InitRunOn::Host);

    for (int k = 0; k < col.layout.nlev; ++k) {
        thetal(0,k,0) = 299.0 + 0.3 * k;
        qv(0,k,0) = 0.010;
        qc(0,k,0) = 0.0;
        qi(0,k,0) = 0.0;
        qw(0,k,0) = qv(0,k,0);
        u(0,k,0) = 2.0 + 0.2 * k;
        v(0,k,0) = 1.0 + 0.1 * k;
    }

    qc(0,2,0) = 2.0e-3;
    qw(0,2,0) = qv(0,2,0) + qc(0,2,0);

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
        ShocStructure::diagnose_pblh(col);
    });
    const auto pblh_consistent_qw = col.pblh.const_array()(0,0,0);

    qw(0,2,0) = qv(0,2,0) + qc(0,2,0) + 0.02;

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
        ShocStructure::diagnose_pblh(col);
    });
    const auto pblh_inconsistent_qw = col.pblh.const_array()(0,0,0);

    EXPECT_NEAR(pblh_consistent_qw, pblh_inconsistent_qw, 1.0e-10);
}

// Motivation: A bulk Richardson number that first exceeds the critical value
// above the reference level must locate the crossing between the two levels;
// this prevents the first crossing from being rounded up to the upper height.
TEST(ShocStructure, PblHeightInterpolatesFirstStableRichardsonCrossing)
{
    auto col = make_first_level_ri_crossing_column(amrex::Real(0.0));

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
    });

    const auto zt = col.zt.const_array();
    const auto zi = col.zi.const_array();
    const auto thetal = col.thetal.const_array();
    const amrex::Real ustar = col.ustar.const_array()(0,0,0);
    const amrex::Real z0_agl = shoc::height_agl(zt(0,0,0), zi(0,0,0));
    const amrex::Real z1_agl = shoc::height_agl(zt(0,1,0), zi(0,0,0));
    const amrex::Real theta0 = thetal(0,0,0);
    const amrex::Real theta1 = thetal(0,1,0);
    const amrex::Real ri1 = CONST_GRAV * (theta1 - theta0) * (z1_agl - z0_agl) /
                            (theta0 * (amrex::Real(100.0) * ustar * ustar));
    ASSERT_TRUE(std::isfinite(ri1));
    ASSERT_GT(ri1, amrex::Real(0.3));

    // The first-level reference Richardson number is zero, so interpolate the
    // critical crossing from the represented Ri(1), not from its nominal input.
    const amrex::Real expected = z0_agl + (amrex::Real(0.3) / ri1) *
                                           (z1_agl - z0_agl);
    const amrex::Real tolerance = amrex::max(
        amrex::Real(1.0e-8),
        amrex::Real(1000.0) * std::numeric_limits<amrex::Real>::epsilon() * expected);

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_pblh(col);
    });
    const auto pblh = col.pblh.const_array()(0,0,0);

    EXPECT_GT(pblh, z0_agl);
    EXPECT_LT(pblh, z1_agl);
    EXPECT_NEAR(pblh, expected, tolerance);
}

// Motivation: Positive surface buoyancy runs a second Richardson search for
// the convective correction; this protects that duplicated first-level branch
// from retaining the old upper-level shortcut after the stable search is fixed.
TEST(ShocStructure, PblHeightInterpolatesFirstConvectiveRichardsonCrossing)
{
    auto col = make_first_level_ri_crossing_column(amrex::Real(1.0e-8));

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
        ShocStructure::diagnose_pblh(col);
    });

    const auto pblh = col.pblh.const_array()(0,0,0);
    const auto zt = col.zt.const_array();
    const auto zi = col.zi.const_array();
    const amrex::Real z0_agl = shoc::height_agl(zt(0,0,0), zi(0,0,0));
    const amrex::Real z1_agl = shoc::height_agl(zt(0,1,0), zi(0,0,0));

    EXPECT_GT(pblh, z0_agl);
    EXPECT_LT(pblh, z1_agl);
}

TEST(ShocStructure, PblHeightStaysInsideColumn)
{
    auto col = shoc_test::make_column(6);
    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
        ShocStructure::diagnose_pblh(col);
    });

    const auto pblh = col.pblh.const_array();
    const auto zt = col.zt.const_array();
    const auto zi = col.zi.const_array();

    EXPECT_GE(pblh(0,0,0), shoc::height_agl(zt(0,0,0), zi(0,0,0)));
    EXPECT_LE(pblh(0,0,0), shoc::height_agl(zt(0,col.layout.nlev-1,0), zi(0,0,0)));
}

TEST(ShocStructure, LowestCloudLayerAppliesVentilationFloor)
{
    auto col = shoc_test::make_column(6);
    auto qc = col.qc.array();
    qc(0,0,0) = 5.0e-4;

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
        ShocStructure::diagnose_pblh(col);
    });

    const auto pblh = col.pblh.const_array();
    const auto zi = col.zi.const_array();
    EXPECT_GE(pblh(0,0,0), shoc::height_agl(zi(0,1,0), zi(0,0,0)) + 50.0);
}

TEST(ShocStructure, BruntInvariantToTerrainOffset)
{
    auto base = shoc_test::make_column(6);
    auto shifted = shoc_test::make_column(6);
    shift_column_heights(shifted, amrex::Real(750.0));
    ShocRuntimeOptions opts;

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(base);
        ShocStructure::diagnose_pblh(base);
        ShocStructure::diagnose_length_and_brunt(base, opts, 2000.0, 2000.0);
        ShocStructure::diagnose_surface_layer(shifted);
        ShocStructure::diagnose_pblh(shifted);
        ShocStructure::diagnose_length_and_brunt(shifted, opts, 2000.0, 2000.0);
    });

    const auto base_brunt = base.brunt.const_array();
    const auto shifted_brunt = shifted.brunt.const_array();
    for (int k = 0; k < base.layout.nlev; ++k) {
        EXPECT_NEAR(base_brunt(0,k,0), shifted_brunt(0,k,0), 1.0e-12);
    }
}

TEST(ShocStructure, LengthScaleInvariantToTerrainOffset)
{
    auto base = shoc_test::make_column(6);
    auto shifted = shoc_test::make_column(6);
    shift_column_heights(shifted, amrex::Real(750.0));
    ShocRuntimeOptions opts;

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(base);
        ShocStructure::diagnose_pblh(base);
        ShocStructure::diagnose_length_and_brunt(base, opts, 2000.0, 2000.0);
        ShocStructure::diagnose_surface_layer(shifted);
        ShocStructure::diagnose_pblh(shifted);
        ShocStructure::diagnose_length_and_brunt(shifted, opts, 2000.0, 2000.0);
    });

    const auto base_mix = base.shoc_mix.const_array();
    const auto shifted_mix = shifted.shoc_mix.const_array();
    for (int k = 0; k < base.layout.nlev; ++k) {
        EXPECT_NEAR(base_mix(0,k,0), shifted_mix(0,k,0), 1.0e-12);
    }
}

TEST(ShocStructure, PblHeightInvariantToTerrainOffsetAndReturnsAgl)
{
    auto base = shoc_test::make_column(6);
    auto shifted = shoc_test::make_column(6);
    shift_column_heights(shifted, amrex::Real(750.0));

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(base);
        ShocStructure::diagnose_pblh(base);
        ShocStructure::diagnose_surface_layer(shifted);
        ShocStructure::diagnose_pblh(shifted);
    });

    const auto base_pblh = base.pblh.const_array();
    const auto shifted_pblh = shifted.pblh.const_array();
    const auto shifted_zt = shifted.zt.const_array();
    const auto shifted_zi = shifted.zi.const_array();
    const amrex::Real z_sfc = shifted_zi(0,0,0);

    EXPECT_NEAR(base_pblh(0,0,0), shifted_pblh(0,0,0), 1.0e-12);
    EXPECT_LT(shifted_pblh(0,0,0), amrex::Real(750.0));
    EXPECT_LE(shifted_pblh(0,0,0), shoc::height_agl(shifted_zt(0,shifted.layout.nlev-1,0), z_sfc));
}

TEST(ShocStructure, LengthScaleRespectsBoundsAndBruntIsFinite)
{
    auto col = shoc_test::make_column(5);
    ShocRuntimeOptions opts;

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_surface_layer(col);
        ShocStructure::diagnose_pblh(col);
        ShocStructure::diagnose_length_and_brunt(col, opts, 400.0, 900.0);
    });

    const auto brunt = col.brunt.const_array();
    const auto mix = col.shoc_mix.const_array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_TRUE(std::isfinite(brunt(0,k,0)));
        EXPECT_TRUE(std::isfinite(mix(0,k,0)));
        EXPECT_GE(mix(0,k,0), 20.0);
        EXPECT_LE(mix(0,k,0), std::sqrt(400.0 * 900.0));
    }
}

TEST(ShocStructure, StableColumnHasPositiveBruntInBottomUpErfOrdering)
{
    auto col = shoc_test::make_column(5);
    ShocRuntimeOptions opts;
    auto thetal = col.thetal.array();
    auto qw = col.qw.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        thetal(0,k,0) = 298.0 + 0.01 * col.zt.const_array()(0,k,0);
        qw(0,k,0) = 0.0;
        qc(0,k,0) = 0.0;
        qi(0,k,0) = 0.0;
    }

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_length_and_brunt(col, opts, 400.0, 900.0);
    });

    const auto brunt = col.brunt.const_array();
    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_GT(brunt(0,k,0), 0.0);
    }
}

TEST(ShocStructure, UnstableColumnHasNegativeBruntInBottomUpErfOrdering)
{
    auto col = shoc_test::make_column(5);
    ShocRuntimeOptions opts;
    auto thetal = col.thetal.array();
    auto qw = col.qw.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        thetal(0,k,0) = 302.0 - 0.01 * col.zt.const_array()(0,k,0);
        qw(0,k,0) = 0.0;
        qc(0,k,0) = 0.0;
        qi(0,k,0) = 0.0;
    }

    shoc_test::run_and_sync([&] {
        ShocStructure::diagnose_length_and_brunt(col, opts, 400.0, 900.0);
    });

    const auto brunt = col.brunt.const_array();
    for (int k = 0; k < col.layout.nlev; ++k) {
        EXPECT_LT(brunt(0,k,0), 0.0);
    }
}

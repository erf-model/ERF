#include <cmath>
#include <memory>
#include <utility>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <gtest/gtest.h>

#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_LandSurface.H>
#include <ERF_RadStruct.H>
#include <ERF_TwoStreamRadiation.H>

namespace {

amrex::BoxArray collapse_z (const amrex::BoxArray& ba)
{
    amrex::BoxList boxes = ba.boxList();
    for (auto& box : boxes) {
        box.setRange(2, 0);
    }
    return amrex::BoxArray(std::move(boxes));
}

amrex::Real component_at (const amrex::MultiFab& mf,
                          const amrex::IntVect& point,
                          int comp)
{
    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    const amrex::Box point_box(point, point, mf.boxArray().ixType());
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box overlap = point_box & mfi.validbox();
        if (overlap.isEmpty()) { continue; }
        const auto array = mf.const_array(mfi);
        reduce_op.eval(overlap, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Real>
            {
                return {array(i, j, k, comp)};
            });
    }
    amrex::Gpu::streamSynchronize();
    // Only the rank that owns the point contributes; the sum gives every rank its value.
    amrex::Real value = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceRealSum(value);
    return value;
}

struct DriverResult
{
    amrex::Real t_sfc = 0.0;
    amrex::Real lw_up_surface = 0.0;
};

DriverResult run_driver (bool external_temperature_provider)
{
    using namespace amrex;

    const Box domain(IntVect(0, 0, 0), IntVect(0, 0, 3));
    const RealBox real_box({0.0, 0.0, 0.0}, {100.0, 100.0, 400.0});
    const int is_periodic[3] = {1, 1, 0};
    const Geometry geom(domain, &real_box, 0, is_periodic);
    const BoxArray ba(domain);
    const BoxArray ba2d = collapse_z(ba);
    const DistributionMapping dm(ba);

    MultiFab state(ba, dm, RhoQ2_comp + 1, 0);
    state.setVal(Real(0.0));
    state.setVal(Real(1.0), Rho_comp, 1);
    state.setVal(getThgivenRandT(Real(1.0), Real(290.0), RdoCp), RhoTheta_comp, 1);

    RadChoice rad;
    rad.enabled = true;
    rad.sw_enabled = false;
    rad.lw_enabled = true;
    rad.fixed_solar_zenith_angle = Real(0.5);
    rad.fixed_total_solar_irradiance = Real(1360.9);
    rad.rad_t_sfc = Real(300.0);
    rad.surface_emissivity_lw = Real(1.0);
    rad.seb_enable = true;
    rad.seb_use_radiation_fluxes = false;
    rad.seb_prognostic_enable = true;
    rad.seb_sw_flux_default = Real(200.0);
    rad.seb_lw_flux_default = Real(0.0);
    rad.seb_hfx_default = Real(0.0);
    rad.seb_lh_default = Real(0.0);
    rad.seb_grdflx_default = Real(0.0);
    rad.seb_t_deep_default = Real(300.0);
    rad.seb_surface_heat_capacity = Real(20000.0);

    TwoStreamRadiation radiation;
    radiation.resize(1);
    radiation.define_level(0, rad, RdoCp, ba2d, dm, ba, domain);

    LandSurface lsm;
    lsm.ReSize(1);
    lsm.SetModel<NullSurf>();

    MultiFab qheating(ba, dm, 2, 0);
    const BoxArray flux_ba = convert(ba, IntVect(0, 0, 1));
    MultiFab rad_fluxes(flux_ba, dm, 4, 0);
    MultiFab external_tsurf(ba2d, dm, 1, 0);
    external_tsurf.setVal(Real(280.0));

    Vector<const MultiFab*> radiation_inputs;
    if (external_temperature_provider) {
        // Represent the canonical radiation input used by SLM's active tsurf.
        radiation_inputs.push_back(&external_tsurf);
    }

    radiation.advance(0, 1, Real(0.0), Real(1000.0), "pre_dycore",
                      state, nullptr, geom, lsm, radiation_inputs, false,
                      &qheating, &rad_fluxes, nullptr, nullptr, nullptr,
                      0.0, false);
    const Real lw_up_surface = component_at(
        rad_fluxes, IntVect(0, 0, 0), 2);
    radiation.advance(0, 1, Real(1000.0), Real(1000.0), "post_dycore",
                      state, nullptr, geom, lsm, radiation_inputs, false,
                      nullptr, nullptr, nullptr, nullptr, nullptr,
                      0.0, false);

    const MultiFab* prognostic_t_sfc =
        radiation.prognostic_surface_temperature_state(0);
    AMREX_ALWAYS_ASSERT(prognostic_t_sfc != nullptr);
    return {prognostic_t_sfc->min(0), lw_up_surface};
}

} // namespace

// This exercises TwoStreamRadiation::advance rather than only the per-column
// resolver: the canonical provider must remain the LW boundary, and ownership
// must prevent an unused TwoStream force-restore state from advancing.
TEST(TwoStreamRadiationDriver, ExternalCanonicalTemperatureOwnsBoundary)
{
    constexpr amrex::Real sigma = 5.670374419e-8;
    const DriverResult no_provider = run_driver(false);
    const DriverResult external_provider = run_driver(true);
    const amrex::Real flux_tolerance = amrex::Real(1.0e-9) * sigma *
        std::pow(amrex::Real(300.0), 4);

    EXPECT_GT(no_provider.t_sfc, amrex::Real(300.0));
    EXPECT_NEAR(no_provider.lw_up_surface,
                sigma * std::pow(amrex::Real(300.0), 4), flux_tolerance);
    EXPECT_EQ(external_provider.t_sfc, amrex::Real(300.0));
    EXPECT_NEAR(external_provider.lw_up_surface,
                sigma * std::pow(amrex::Real(280.0), 4), flux_tolerance);
}

namespace {

struct LandForcingResult
{
    bool supplied = false;
    amrex::Real sw_dn = 0.0, lw_dn = 0.0, coszen = 0.0;
    amrex::Real rad_sw_dn_sfc = 0.0, rad_lw_dn_sfc = 0.0, rad_sw_dn_toa = 0.0;
    amrex::Real lsm_sw_dn = 0.0, lsm_lw_dn = 0.0, lsm_coszen = 0.0;
};

LandForcingResult run_land_forcing (bool supply)
{
    using namespace amrex;

    const Box domain(IntVect(0, 0, 0), IntVect(1, 0, 3));
    const RealBox real_box({0.0, 0.0, 0.0}, {200.0, 100.0, 4000.0});
    const int is_periodic[3] = {1, 1, 0};
    const Geometry geom(domain, &real_box, 0, is_periodic);
    const BoxArray ba(domain);
    const BoxArray ba2d = collapse_z(ba);
    const DistributionMapping dm(ba);

    MultiFab state(ba, dm, RhoQ2_comp + 1, 0);
    state.setVal(Real(0.0));
    state.setVal(Real(1.0), Rho_comp, 1);
    state.setVal(getThgivenRandT(Real(1.0), Real(290.0), RdoCp), RhoTheta_comp, 1);

    // Both bands on, an absorbing and scattering atmosphere, and a cloudy fraction, so
    // the surface SW differs from the top and the fluxes are a clear/cloudy blend.
    RadChoice rad;
    rad.enabled = true;
    rad.sw_enabled = true;
    rad.lw_enabled = true;
    rad.fixed_solar_zenith_angle = Real(0.6);
    rad.fixed_total_solar_irradiance = Real(1360.9);
    rad.rad_t_sfc = Real(300.0);
    rad.tau_per_layer = Real(0.1);
    rad.single_scattering_albedo = Real(0.5);
    rad.cloud_fraction = Real(0.4);

    TwoStreamRadiation radiation;
    radiation.resize(1);
    radiation.define_level(0, rad, RdoCp, ba2d, dm, ba, domain, supply);

    LandForcingResult r;
    r.supplied = radiation.supplies_land_forcing(0);
    if (!supply) {
        EXPECT_EQ(radiation.land_forcing_sw_dn(0), nullptr);
        EXPECT_EQ(radiation.land_forcing_lw_dn(0), nullptr);
        EXPECT_EQ(radiation.land_forcing_cos_zenith(0), nullptr);
        return r;
    }

    LandSurface lsm;
    lsm.ReSize(1);
    lsm.SetModel<NullSurf>();

    MultiFab qheating(ba, dm, 2, 0);
    const BoxArray flux_ba = convert(ba, IntVect(0, 0, 1));
    MultiFab rad_fluxes(flux_ba, dm, 4, 0);
    const Vector<const MultiFab*> radiation_inputs;

    radiation.advance(0, 1, Real(0.0), Real(10.0), "pre_dycore",
                      state, nullptr, geom, lsm, radiation_inputs, false,
                      &qheating, &rad_fluxes, nullptr, nullptr, nullptr,
                      0.0, false);

    r.sw_dn  = component_at(*radiation.land_forcing_sw_dn(0), IntVect(1, 0, 0), 0);
    r.lw_dn  = component_at(*radiation.land_forcing_lw_dn(0), IntVect(1, 0, 0), 0);
    r.coszen = component_at(*radiation.land_forcing_cos_zenith(0), IntVect(1, 0, 0), 0);
    r.rad_sw_dn_sfc = component_at(rad_fluxes, IntVect(1, 0, 0), 1);
    r.rad_lw_dn_sfc = component_at(rad_fluxes, IntVect(1, 0, 0), 3);
    r.rad_sw_dn_toa = component_at(rad_fluxes, IntVect(1, 0, 4), 1);

    // The copy into a land model's layout (Noah-MP: the 2D boxes, x/y ghost cells).
    MultiFab lsm_sw(ba2d, dm, 1, IntVect(1, 1, 0));
    MultiFab lsm_lw(ba2d, dm, 1, IntVect(1, 1, 0));
    MultiFab lsm_cz(ba2d, dm, 1, IntVect(1, 1, 0));
    lsm_sw.setVal(Real(-1.0));
    lsm_lw.setVal(Real(-1.0));
    lsm_cz.setVal(Real(-1.0));
    radiation.write_land_forcing(0, &lsm_sw, &lsm_lw, &lsm_cz);
    r.lsm_sw_dn  = component_at(lsm_sw, IntVect(1, 0, 0), 0);
    r.lsm_lw_dn  = component_at(lsm_lw, IntVect(1, 0, 0), 0);
    r.lsm_coszen = component_at(lsm_cz, IntVect(1, 0, 0), 0);
    return r;
}

} // namespace

// The forcing a land model receives is the sweep's own surface-interface downwelling
// flux -- after the clear/cloudy blend -- and this call's sun, and it reaches the land
// model's field unchanged. Without the request nothing is allocated.
TEST(TwoStreamRadiationDriver, SuppliesLandForcingFromTheSweep)
{
    const LandForcingResult off = run_land_forcing(false);
    EXPECT_FALSE(off.supplied);

    const LandForcingResult on = run_land_forcing(true);
    ASSERT_TRUE(on.supplied);
    EXPECT_GT(on.rad_sw_dn_sfc, amrex::Real(0.0));
    EXPECT_LT(on.rad_sw_dn_sfc, on.rad_sw_dn_toa);  // attenuated: the surface, not the top
    EXPECT_GT(on.rad_lw_dn_sfc, amrex::Real(0.0));
    EXPECT_EQ(on.sw_dn, on.rad_sw_dn_sfc);
    EXPECT_EQ(on.lw_dn, on.rad_lw_dn_sfc);
    EXPECT_EQ(on.coszen, amrex::Real(0.6));
    EXPECT_EQ(on.lsm_sw_dn, on.sw_dn);
    EXPECT_EQ(on.lsm_lw_dn, on.lw_dn);
    EXPECT_EQ(on.lsm_coszen, on.coszen);
}

namespace {

// What write_land_forcing copies on a level that has not swept: level 0 before its first
// advance(), and a level-1 patch that stops below the domain top, whose advance() returns
// before sweeping (its radiation comes from the parent by interpolation).
std::pair<amrex::Real, amrex::Real> land_forcing_before_a_sweep ()
{
    using namespace amrex;

    const Box domain0(IntVect(0, 0, 0), IntVect(1, 0, 3));
    const Box domain1(IntVect(0, 0, 0), IntVect(3, 1, 7));
    const RealBox real_box({0.0, 0.0, 0.0}, {200.0, 100.0, 4000.0});
    const int is_periodic[3] = {1, 1, 0};
    const Geometry geom1(domain1, &real_box, 0, is_periodic);

    const BoxArray ba0(domain0);
    const BoxArray ba1(Box(IntVect(0, 0, 0), IntVect(3, 1, 3)));  // shallow: k = 0..3 of 0..7
    const DistributionMapping dm0(ba0);
    const DistributionMapping dm1(ba1);

    RadChoice rad;
    rad.enabled = true;
    rad.fixed_solar_zenith_angle = Real(0.6);
    rad.fixed_total_solar_irradiance = Real(1360.9);
    rad.rad_t_sfc = Real(300.0);

    TwoStreamRadiation radiation;
    radiation.resize(2);
    radiation.define_level(0, rad, RdoCp, collapse_z(ba0), dm0, ba0, domain0, true);
    radiation.define_level(1, rad, RdoCp, collapse_z(ba1), dm1, ba1, domain1, true);

    // Level 0, never advanced.
    MultiFab lsm0(collapse_z(ba0), dm0, 1, IntVect(1, 1, 0));
    lsm0.setVal(Real(7.0));
    radiation.write_land_forcing(0, &lsm0, nullptr, nullptr);

    // Level 1, advanced -- but the call returns before the sweep.
    MultiFab state1(ba1, dm1, RhoQ2_comp + 1, 0);
    state1.setVal(Real(1.0));
    LandSurface lsm;
    lsm.ReSize(2);
    lsm.SetModel<NullSurf>();
    const Vector<const MultiFab*> radiation_inputs;
    radiation.advance(1, 1, Real(0.0), Real(10.0), "pre_dycore",
                      state1, nullptr, geom1, lsm, radiation_inputs, false,
                      nullptr, nullptr, nullptr, nullptr, nullptr, 0.0, false);
    MultiFab lsm1(collapse_z(ba1), dm1, 1, IntVect(1, 1, 0));
    lsm1.setVal(Real(7.0));
    radiation.write_land_forcing(1, nullptr, nullptr, &lsm1);

    return {component_at(lsm0, IntVect(1, 0, 0), 0), component_at(lsm1, IntVect(2, 1, 0), 0)};
}

} // namespace

// A copy made before any sweep on the level must not hand the land model a valid zero:
// it carries lsm_undefined, which Noah-MP's first-land-step check rejects (zero would
// pass it, and be a 0 K sky).
TEST(TwoStreamRadiationDriver, LandForcingIsUndefinedUntilASweep)
{
    const auto [lev0, lev1] = land_forcing_before_a_sweep();
    EXPECT_FALSE(is_valid_lsm_value(lev0));
    EXPECT_FALSE(is_valid_lsm_value(lev1));
    EXPECT_EQ(lev0, lsm_undefined);
    EXPECT_EQ(lev1, lsm_undefined);
}

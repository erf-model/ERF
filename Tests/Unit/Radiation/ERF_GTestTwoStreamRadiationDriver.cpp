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
    amrex::Real seb_hfx = 0.0;
    amrex::Real seb_lh = 0.0;
};

// The surface layer's applied fluxes handed to advance, in W/m^2; none when
// surface_layer is false (no zlo surface layer).
struct SurfaceLayerFluxes
{
    bool surface_layer = false;
    amrex::Real hfx_wm2 = 0.0;
    amrex::Real lh_wm2 = 0.0;
    SEBTurbulentFluxSource source = SEBTurbulentFluxSource::SurfaceLayer;
};

DriverResult run_driver (bool external_temperature_provider,
                         const SurfaceLayerFluxes& sl = SurfaceLayerFluxes{})
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
    rad.seb_turbulent_flux_source = sl.source;

    TwoStreamRadiation radiation;
    radiation.resize(1);
    radiation.define_level(0, rad, RdoCp, ba2d, dm, ba, domain);

    LandSurface lsm;
    lsm.ReSize(1);
    lsm.SetModel<NullSurf>();

    MultiFab qheating(ba, dm, 2, 0);
    const BoxArray flux_ba = convert(ba, IntVect(0, 0, 1));
    MultiFab rad_fluxes(flux_ba, dm, 4, 0);
    // The conservative fluxes the surface layer would write on the z faces
    // (rho w'theta' and rho w'qv'): the SEB converts them back to W/m^2.
    MultiFab sfc_sens_flux(flux_ba, dm, 1, 0);
    MultiFab sfc_laten_flux(flux_ba, dm, 1, 0);
    sfc_sens_flux.setVal(sl.hfx_wm2 / Cp_d);
    sfc_laten_flux.setVal(sl.lh_wm2 / L_v);
    const MultiFab* sens_ptr = sl.surface_layer ? &sfc_sens_flux : nullptr;
    const MultiFab* laten_ptr = sl.surface_layer ? &sfc_laten_flux : nullptr;
    MultiFab external_tsurf(ba2d, dm, 1, 0);
    external_tsurf.setVal(Real(280.0));

    Vector<const MultiFab*> radiation_inputs;
    if (external_temperature_provider) {
        // Represent the canonical radiation input used by SLM's active tsurf.
        radiation_inputs.push_back(&external_tsurf);
    }

    radiation.advance(0, 1, Real(0.0), Real(1000.0), "pre_dycore",
                      state, nullptr, geom, lsm, radiation_inputs, false,
                      &qheating, &rad_fluxes, nullptr, sens_ptr, laten_ptr, nullptr, nullptr,
                      0.0, false);
    const Real lw_up_surface = component_at(
        rad_fluxes, IntVect(0, 0, 0), 2);
    radiation.advance(0, 1, Real(1000.0), Real(1000.0), "post_dycore",
                      state, nullptr, geom, lsm, radiation_inputs, false,
                      nullptr, nullptr, nullptr, sens_ptr, laten_ptr, nullptr, nullptr,
                      0.0, false);

    const MultiFab* prognostic_t_sfc =
        radiation.prognostic_surface_temperature_state(0);
    AMREX_ALWAYS_ASSERT(prognostic_t_sfc != nullptr);
    return {prognostic_t_sfc->min(0), lw_up_surface,
            radiation.seb_hfx(0)->min(0), radiation.seb_lh(0)->min(0)};
}

// The force-restore skin after one post-dycore call on a fresh level whose surface
// radiation comes from the sweep (seb_use_radiation_fluxes), with or without the
// pre-dycore sweep of that step before it. The SW default is a large +200 W/m^2 and the
// SW band is off, so a skin that advanced on the defaults warms by 10 K while one that
// advanced on the sweep's (longwave-only) fluxes cools.
amrex::Real force_restore_skin_after_post_dycore (bool sweep_first)
{
    using namespace amrex;

    const Box domain(IntVect(0, 0, 0), IntVect(0, 0, 3));
    const RealBox real_box({0.0, 0.0, 0.0}, {100.0, 100.0, 400.0});
    const int is_periodic[3] = {1, 1, 0};
    const Geometry geom(domain, &real_box, 0, is_periodic);
    const BoxArray ba(domain);
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
    rad.seb_use_radiation_fluxes = true;
    rad.seb_prognostic_enable = true;
    rad.seb_sw_flux_default = Real(200.0);
    rad.seb_lw_flux_default = Real(0.0);
    rad.seb_hfx_default = Real(0.0);
    rad.seb_lh_default = Real(0.0);
    rad.seb_grdflx_default = Real(0.0);
    rad.seb_t_deep_default = Real(300.0);
    rad.seb_surface_heat_capacity = Real(20000.0);
    rad.seb_turbulent_flux_source = SEBTurbulentFluxSource::Defaults;

    TwoStreamRadiation radiation;
    radiation.resize(1);
    radiation.define_level(0, rad, RdoCp, collapse_z(ba), dm, ba, domain);

    LandSurface lsm;
    lsm.ReSize(1);
    lsm.SetModel<NullSurf>();

    MultiFab qheating(ba, dm, 2, 0);
    MultiFab rad_fluxes(convert(ba, IntVect(0, 0, 1)), dm, 4, 0);
    const Vector<const MultiFab*> radiation_inputs {};

    if (sweep_first) {
        radiation.advance(0, 1, Real(0.0), Real(1000.0), "pre_dycore",
                          state, nullptr, geom, lsm, radiation_inputs, false,
                          &qheating, &rad_fluxes, nullptr, nullptr, nullptr, nullptr, nullptr,
                          0.0, false);
    }
    radiation.advance(0, 1, Real(1000.0), Real(1000.0), "post_dycore",
                      state, nullptr, geom, lsm, radiation_inputs, false,
                      nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                      0.0, false);

    const MultiFab* t_sfc = radiation.prognostic_surface_temperature_state(0);
    AMREX_ALWAYS_ASSERT(t_sfc != nullptr);
    return t_sfc->min(0);
}

} // namespace

// The first step of a level built by interp_atmos_from_coarse: ERF::advance_radiation
// skips the level's pre-dycore sweep, but the post-dycore call still arrives. The
// force-restore balance, which takes its surface radiation from the sweep, must not
// advance the skin on the values define_level left there (the scalar defaults).
TEST(TwoStreamRadiationDriver, ForceRestoreWaitsForTheLevelsFirstSweep)
{
    using amrex::Real;
    const Real initial = Real(300.0);

    const Real without_sweep = force_restore_skin_after_post_dycore(false);
    EXPECT_EQ(without_sweep, initial)
        << "the skin advanced before the level's first sweep (on the defaults: +10 K)";

    // Once the level has swept, the balance advances on the sweep's fluxes: longwave only,
    // a 300 K surface under 290 K air, so it cools.
    const Real with_sweep = force_restore_skin_after_post_dycore(true);
    EXPECT_LT(with_sweep, initial);
    EXPECT_GT(with_sweep, initial - Real(5.0));
}

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

// The force-restore update removes the surface layer's H and LE from the ground.
// Absorbed SW is 200 W/m^2 and T_s starts at T_deep, so one 1000 s step moves T_s
// by 1000 (200 - H - LE) / 2e4: +10 K with no turbulent fluxes, +5 K when the
// surface layer carries H = 60 and LE = 40 W/m^2 away. Before the surface layer
// was consulted, both runs warmed by 10 K.
TEST(TwoStreamRadiationDriver, ForceRestoreRemovesSurfaceLayerFluxes)
{
    const amrex::Real tol = sizeof(amrex::Real) == 8 ? amrex::Real(1.0e-9) : amrex::Real(1.0e-3);

    const DriverResult no_surface_layer = run_driver(false);
    EXPECT_NEAR(no_surface_layer.t_sfc, amrex::Real(310.0), tol);
    EXPECT_EQ(no_surface_layer.seb_hfx, amrex::Real(0.0));

    SurfaceLayerFluxes sl;
    sl.surface_layer = true;
    sl.hfx_wm2 = amrex::Real(60.0);
    sl.lh_wm2 = amrex::Real(40.0);
    const DriverResult with_surface_layer = run_driver(false, sl);
    EXPECT_NEAR(with_surface_layer.seb_hfx, amrex::Real(60.0), tol);
    EXPECT_NEAR(with_surface_layer.seb_lh, amrex::Real(40.0), tol);
    EXPECT_NEAR(with_surface_layer.t_sfc, amrex::Real(305.0), tol);

    // seb_turbulent_flux_source = defaults keeps the old constants (0 here).
    sl.source = SEBTurbulentFluxSource::Defaults;
    const DriverResult defaults = run_driver(false, sl);
    EXPECT_EQ(defaults.seb_hfx, amrex::Real(0.0));
    EXPECT_NEAR(defaults.t_sfc, amrex::Real(310.0), tol);
}

namespace {

struct LandForcingResult
{
    bool supplied = false;
    amrex::Real sw_dn = 0.0, lw_dn = 0.0, coszen = 0.0;
    amrex::Real rad_sw_dn_sfc = 0.0, rad_lw_dn_sfc = 0.0, rad_sw_dn_toa = 0.0;
    amrex::Real lsm_sw_dn = 0.0, lsm_lw_dn = 0.0, lsm_coszen = 0.0;
};

LandForcingResult run_land_forcing (bool supply, amrex::Real cloud_fraction)
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

    // Both bands on, an absorbing and scattering atmosphere, and a cloud layer over the
    // middle two of the four 1 km layers (centres at 1.5 and 2.5 km). With a cloud
    // fraction the surface fluxes are then a real clear/cloudy blend: the cloudy column
    // differs from the clear one, which the test checks, so storing the forcing before
    // the blend would be caught.
    RadChoice rad;
    rad.enabled = true;
    rad.sw_enabled = true;
    rad.lw_enabled = true;
    rad.fixed_solar_zenith_angle = Real(0.6);
    rad.fixed_total_solar_irradiance = Real(1360.9);
    rad.rad_t_sfc = Real(300.0);
    rad.tau_per_layer = Real(0.1);
    rad.single_scattering_albedo = Real(0.5);
    rad.tau_profile_type = TauProfileType::CloudLayer;
    rad.cloud_base_height_m = Real(1000.0);
    rad.cloud_top_height_m = Real(3000.0);
    rad.cloud_tau_per_layer = Real(2.0);
    rad.cloud_fraction = cloud_fraction;

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
    const Vector<const MultiFab*> radiation_inputs {};

    radiation.advance(0, 1, Real(0.0), Real(10.0), "pre_dycore",
                      state, nullptr, geom, lsm, radiation_inputs, false,
                      &qheating, &rad_fluxes, nullptr, nullptr, nullptr, nullptr, nullptr,
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
    const LandForcingResult off = run_land_forcing(false, amrex::Real(0.4));
    EXPECT_FALSE(off.supplied);

    const LandForcingResult on = run_land_forcing(true, amrex::Real(0.4));
    ASSERT_TRUE(on.supplied);
    // The cloud must matter at the surface, or the blend could not be told apart.
    const LandForcingResult clear = run_land_forcing(true, amrex::Real(0.0));
    EXPECT_LT(on.rad_sw_dn_sfc, amrex::Real(0.99) * clear.rad_sw_dn_sfc);
    EXPECT_NE(on.rad_lw_dn_sfc, clear.rad_lw_dn_sfc);
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
    const Vector<const MultiFab*> radiation_inputs {};
    radiation.advance(1, 1, Real(0.0), Real(10.0), "pre_dycore",
                      state1, nullptr, geom1, lsm, radiation_inputs, false,
                      nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                      0.0, false);
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

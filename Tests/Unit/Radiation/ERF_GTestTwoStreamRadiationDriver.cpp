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
    return amrex::get<0>(reduce_data.value());
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

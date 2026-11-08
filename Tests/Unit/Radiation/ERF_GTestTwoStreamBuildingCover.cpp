#include <cmath>
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

// The two-stream columns seeing the buildings (TwoStreamRadiation::set_building_cover(),
// erf.ibseb.two_stream_buildings). A 4 x 4 level of 100 m columns in four boxes under a
// fixed sun, an absorbing and scattering sky, with a building over column (0,0) (its
// bottom two cells solid, a footprint) and an overhang over column (3,3) (one solid cell
// at k = 2, open ground below it). The sweep must write no heating in the three solid
// cells while heating the air; the ground balance must hold the footprint, give it no
// shortwave, and take the absorbed beam, (1 - albedo) S0 cos z exp(-tau N / cos z) with no
// scattering of the beam, off the column the shadow field marks, once however often the
// update after the dycore is called; the ground under the overhang is open ground.

namespace {

amrex::BoxArray collapse_z (const amrex::BoxArray& ba)
{
    amrex::BoxList boxes = ba.boxList();
    for (auto& box : boxes) { box.setRange(2, 0); }
    return amrex::BoxArray(std::move(boxes));
}

amrex::Real value_at (const amrex::MultiFab& mf, const amrex::IntVect& point, int comp)
{
    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    const amrex::Box point_box(point, point, mf.boxArray().ixType());
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box overlap = point_box & mfi.validbox();
        if (overlap.isEmpty()) { continue; }
        const auto array = mf.const_array(mfi);
        reduce_op.eval(overlap, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Real> { return {array(i, j, k, comp)}; });
    }
    amrex::Gpu::streamSynchronize();
    amrex::Real value = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceRealSum(value);
    return value;
}

void set_at (amrex::MultiFab& mf, const amrex::IntVect& point, amrex::Real v)
{
    const amrex::Box point_box(point, point, mf.boxArray().ixType());
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box overlap = point_box & mfi.validbox();
        if (overlap.isEmpty()) { continue; }
        const auto a = mf.array(mfi);
        amrex::ParallelFor(overlap, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { a(i, j, k) = v; });
    }
    amrex::Gpu::streamSynchronize();
}

RadChoice cover_choice (bool seb)
{
    RadChoice rad;
    rad.enabled = true;
    rad.sw_enabled = true;
    rad.lw_enabled = true;
    rad.fixed_solar_zenith_angle = amrex::Real(0.5);
    rad.fixed_total_solar_irradiance = amrex::Real(1000.0);
    rad.tau_per_layer = amrex::Real(0.05);
    rad.tau_lw_per_layer = amrex::Real(0.2);
    rad.single_scattering_albedo = amrex::Real(0.5);
    rad.surface_albedo_sw = amrex::Real(0.2);
    rad.surface_emissivity_lw = amrex::Real(1.0);
    rad.rad_t_sfc = amrex::Real(300.0);
    rad.seb_enable = seb;
    rad.seb_use_radiation_fluxes = seb;
    rad.seb_prognostic_enable = seb;
    rad.seb_hfx_default = amrex::Real(0.0);
    rad.seb_lh_default = amrex::Real(0.0);
    rad.seb_grdflx_default = amrex::Real(0.0);
    rad.seb_t_deep_default = amrex::Real(300.0);
    rad.seb_surface_heat_capacity = amrex::Real(20000.0);
    rad.seb_turbulent_flux_source = SEBTurbulentFluxSource::Defaults;
    return rad;
}

struct Level
{
    amrex::Box domain{amrex::IntVect(0, 0, 0), amrex::IntVect(3, 3, 3)};
    amrex::RealBox real_box{{0.0, 0.0, 0.0}, {400.0, 400.0, 400.0}};
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    Level ()
    {
        const int is_periodic[3] = {1, 1, 0};
        geom.define(domain, &real_box, 0, is_periodic);
        ba = amrex::BoxArray(domain);
        ba.maxSize(amrex::IntVect(2, 2, 4));
        dm = amrex::DistributionMapping(ba);
    }
};

} // namespace

TEST(TwoStreamBuildingCover, ColumnsAndGroundSeeTheBuildings)
{
    using namespace amrex;
    Level L;
    const RadChoice rad = cover_choice(true);
    TwoStreamRadiation radiation;
    radiation.resize(1);
    radiation.define_level(0, rad, RdoCp, collapse_z(L.ba), L.dm, L.ba, L.domain);

    MultiFab blank(L.ba, L.dm, 1, 1);
    blank.setVal(0.0);
    set_at(blank, IntVect(0, 0, 0), 1.0);
    set_at(blank, IntVect(0, 0, 1), 1.0);
    set_at(blank, IntVect(3, 3, 2), 1.0);
    radiation.set_building_cover(0, blank, L.domain);
    ASSERT_TRUE(radiation.has_building_cover(0));

    MultiFab state(L.ba, L.dm, RhoQ2_comp + 1, 0);
    state.setVal(Real(0.0));
    state.setVal(Real(1.0), Rho_comp, 1);
    state.setVal(getThgivenRandT(Real(1.0), Real(290.0), RdoCp), RhoTheta_comp, 1);
    MultiFab qheating(L.ba, L.dm, 2, 0);
    MultiFab rad_fluxes(convert(L.ba, IntVect(0, 0, 1)), L.dm, 4, 0);
    LandSurface lsm;
    lsm.ReSize(1);
    lsm.SetModel<NullSurf>();
    const Vector<const MultiFab*> radiation_inputs {};

    radiation.advance(0, 1, Real(0.0), Real(10.0), "pre_dycore", state, nullptr, L.geom, lsm, radiation_inputs,
                      false, &qheating, &rad_fluxes, nullptr, nullptr, nullptr, nullptr, nullptr, 0.0, false);

    // No heating in the solid cells; the air at the same heights is heated.
    for (const IntVect& p : {IntVect(0, 0, 0), IntVect(0, 0, 1), IntVect(3, 3, 2)}) {
        EXPECT_EQ(value_at(qheating, p, 0), Real(0.0)) << p;
        EXPECT_EQ(value_at(qheating, p, 1), Real(0.0)) << p;
    }
    for (const IntVect& p : {IntVect(1, 2, 0), IntVect(1, 2, 1), IntVect(3, 3, 1), IntVect(3, 3, 3)}) {
        EXPECT_GT(std::abs(value_at(qheating, p, 0)), Real(0.0)) << p;
        EXPECT_GT(std::abs(value_at(qheating, p, 1)), Real(0.0)) << p;
    }

    // The shadow of this step: column (2,1) of the open ground, and the footprint.
    MultiFab* shadow = radiation.ground_shadow_for_step(0, 1);
    ASSERT_NE(shadow, nullptr);
    shadow->setVal(0.0);
    set_at(*shadow, IntVect(2, 1, 0), 1.0);
    set_at(*shadow, IntVect(0, 0, 0), 1.0);
    radiation.advance(0, 1, Real(10.0), Real(10.0), "post_dycore", state, nullptr, L.geom, lsm, radiation_inputs,
                      false, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0.0, false);

    const MultiFab& sw = *radiation.seb_sw_net(0);
    const Real sw_open = value_at(sw, IntVect(1, 2, 0), 0);
    const Real sw_shaded = value_at(sw, IntVect(2, 1, 0), 0);
    const Real beam_abs = Real(0.8) * Real(1000.0) * Real(0.5) * std::exp(Real(-0.05 * 4 / 0.5));
    const Real tol = (sizeof(Real) == 8) ? Real(1.e-12) : Real(1.e-5);
    EXPECT_NEAR(sw_open - sw_shaded, beam_abs, tol * beam_abs);
    EXPECT_GT(sw_shaded, Real(1.0));                                   // the diffuse light stays
    EXPECT_EQ(value_at(sw, IntVect(0, 0, 0), 0), Real(0.0));           // the footprint
    EXPECT_EQ(value_at(sw, IntVect(3, 3, 0), 0), sw_open);             // under the overhang: no footprint, and this test sets no shadow there (compute_ground_shadow() would)

    const MultiFab& t_sfc = *radiation.prognostic_surface_temperature_state(0);
    EXPECT_EQ(value_at(t_sfc, IntVect(0, 0, 0), 0), Real(300.0));      // held
    EXPECT_GT(value_at(t_sfc, IntVect(1, 2, 0), 0), Real(300.0));      // the sunlit ground warms
    EXPECT_GT(value_at(t_sfc, IntVect(3, 3, 0), 0), Real(300.0));

    // A second update after the same sweep takes the beam off no second time.
    radiation.advance(0, 1, Real(10.0), Real(10.0), "post_dycore", state, nullptr, L.geom, lsm, radiation_inputs,
                      false, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0.0, false);
    EXPECT_EQ(value_at(*radiation.seb_sw_net(0), IntVect(2, 1, 0), 0), sw_shaded);

    // A level whose rates come from its parent: the mask sets the solid cells to zero only.
    MultiFab q(L.ba, L.dm, 2, 0);
    q.setVal(1.0);
    radiation.mask_building_heating(0, q);
    EXPECT_EQ(value_at(q, IntVect(0, 0, 1), 1), Real(0.0));
    EXPECT_EQ(value_at(q, IntVect(3, 3, 2), 0), Real(0.0));
    EXPECT_EQ(value_at(q, IntVect(3, 3, 1), 0), Real(1.0));
    EXPECT_EQ(q.sum(0), Real(64 - 3));

    // A rebuilt level drops the cover, as it drops the faces' canopy forcing.
    radiation.define_level(0, rad, RdoCp, collapse_z(L.ba), L.dm, L.ba, L.domain);
    EXPECT_FALSE(radiation.has_building_cover(0));
    EXPECT_EQ(radiation.ground_shadow_for_step(0, 2), nullptr);
}

// Without a ground balance there is no shadow to give; the heating mask still applies.
TEST(TwoStreamBuildingCover, NoGroundBalanceNoShadowField)
{
    using namespace amrex;
    Level L;
    const RadChoice rad = cover_choice(false);
    TwoStreamRadiation radiation;
    radiation.resize(1);
    radiation.define_level(0, rad, RdoCp, collapse_z(L.ba), L.dm, L.ba, L.domain);
    MultiFab blank(L.ba, L.dm, 1, 1);
    blank.setVal(0.0);
    set_at(blank, IntVect(1, 1, 0), 0.7);   // solid at 0.5 and above
    set_at(blank, IntVect(1, 2, 0), 0.4);
    radiation.set_building_cover(0, blank, L.domain);
    EXPECT_TRUE(radiation.has_building_cover(0));
    EXPECT_EQ(radiation.ground_shadow_for_step(0, 1), nullptr);
    MultiFab q(L.ba, L.dm, 2, 0);
    q.setVal(1.0);
    radiation.mask_building_heating(0, q);
    EXPECT_EQ(value_at(q, IntVect(1, 1, 0), 0), Real(0.0));
    EXPECT_EQ(value_at(q, IntVect(1, 2, 0), 0), Real(1.0));
}

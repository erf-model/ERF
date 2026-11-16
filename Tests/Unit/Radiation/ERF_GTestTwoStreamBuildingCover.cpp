#include <cmath>
#include <utility>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <gtest/gtest.h>

#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_NumericalConstants.H>
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

// What the walls send to the ground of one column: their share of its view, their reflected
// shortwave and their longwave (the three components of ground_walls_for_step()).
void give_walls (amrex::MultiFab& walls, const amrex::IntVect& point, amrex::Real view, amrex::Real sw, amrex::Real lw)
{
    const amrex::Box point_box(point, point);
    for (amrex::MFIter mfi(walls); mfi.isValid(); ++mfi) {
        const amrex::Box overlap = point_box & mfi.validbox();
        if (overlap.isEmpty()) { continue; }
        const auto a = walls.array(mfi);
        amrex::ParallelFor(overlap, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            a(i, j, k, 0) = view; a(i, j, k, 1) = sw; a(i, j, k, 2) = lw;
        });
    }
    amrex::Gpu::streamSynchronize();
}

// The roof of one column: set, with its albedo, emissivity, skin temperature and the share
// of the beam that reaches it (the five components of roofs_for_step()).
void give_roof (amrex::MultiFab& roofs, const amrex::IntVect& point, amrex::Real albedo, amrex::Real emissivity, amrex::Real temperature,
                amrex::Real lit = amrex::Real(1.0))
{
    const amrex::Box point_box(point, point);
    for (amrex::MFIter mfi(roofs); mfi.isValid(); ++mfi) {
        const amrex::Box overlap = point_box & mfi.validbox();
        if (overlap.isEmpty()) { continue; }
        const auto a = roofs.array(mfi);
        amrex::ParallelFor(overlap, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            a(i, j, k, 0) = amrex::Real(1.0); a(i, j, k, 1) = albedo; a(i, j, k, 2) = emissivity; a(i, j, k, 3) = temperature;
            a(i, j, k, 4) = lit;
        });
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
    const Real tol = Real(256) * real_eps;   // a multiple of the working precision's roundoff
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

// erf.ibseb.ground_sky_view: at the update after the dycore the open ground takes what
// the walls send it (ground_walls_for_step()), once per sweep, with the beam; the
// footprint stays held. The sky is the same in every column, so the columns differ by
// the walls' terms alone.
TEST(TwoStreamBuildingCover, GroundTakesWhatTheWallsSend)
{
    using namespace amrex;
    Level L;
    RadChoice rad = cover_choice(true);
    rad.surface_emissivity_lw = Real(0.9);
    rad.surface_albedo_sw_diffuse = Real(0.3);
    TwoStreamRadiation radiation;
    radiation.resize(1);
    radiation.define_level(0, rad, RdoCp, collapse_z(L.ba), L.dm, L.ba, L.domain);
    MultiFab blank(L.ba, L.dm, 1, 1);
    blank.setVal(0.0);
    set_at(blank, IntVect(0, 0, 0), 1.0);
    // Without the request there is no field to fill.
    radiation.set_building_cover(0, blank, L.domain);
    EXPECT_EQ(radiation.ground_walls_for_step(0, 1), nullptr);
    radiation.set_building_cover(0, blank, L.domain, true);
    EXPECT_TRUE(radiation.takes_ground_walls(0));
    EXPECT_EQ(radiation.ground_walls(0), nullptr);   // none given yet

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
    const Real lw_down = value_at(rad_fluxes, IntVect(1, 2, 0), 3);
    const Real sw_down = value_at(rad_fluxes, IntVect(1, 2, 0), 1);
    const Real beam = Real(1000.0) * Real(0.5) * std::exp(Real(-0.05 * 4 / 0.5));
    ASSERT_GT(lw_down, Real(10.0));
    ASSERT_GT(sw_down - beam, Real(1.0));

    MultiFab* shadow = radiation.ground_shadow_for_step(0, 1);
    ASSERT_NE(shadow, nullptr);
    shadow->setVal(0.0);
    set_at(*shadow, IntVect(2, 1, 0), 1.0);
    MultiFab* walls = radiation.ground_walls_for_step(0, 1);
    ASSERT_NE(walls, nullptr);
    ASSERT_EQ(walls->nComp(), 3);
    walls->setVal(0.0);
    // Column (2,1), shaded: the walls take 0.4 of the view and send 30 W/m2 of reflected
    // light and 200 W/m2 of longwave. Column (3,2), sunlit: 0.25, 60 and 120. The
    // footprint (0,0) is given values too, which the balance must not take.
    give_walls(*walls, IntVect(2, 1, 0), Real(0.4), Real(30.0), Real(200.0));
    give_walls(*walls, IntVect(3, 2, 0), Real(0.25), Real(60.0), Real(120.0));
    give_walls(*walls, IntVect(0, 0, 0), Real(0.5), Real(99.0), Real(99.0));
    radiation.advance(0, 1, Real(10.0), Real(10.0), "post_dycore", state, nullptr, L.geom, lsm, radiation_inputs,
                      false, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0.0, false);

    const MultiFab& sw = *radiation.seb_sw_net(0);
    const MultiFab& lw = *radiation.seb_lw_net(0);
    const Real tol = Real(256) * real_eps;   // a multiple of the working precision's roundoff
    const Real sw_open = value_at(sw, IntVect(1, 2, 0), 0), lw_open = value_at(lw, IntVect(1, 2, 0), 0);
    const Real dif_abs = Real(0.7) * (sw_down - beam);       // the diffuse light the open ground absorbs
    const Real beam_abs = Real(0.8) * beam;
    EXPECT_NEAR(sw_open, beam_abs + dif_abs, tol * sw_open);
    // Shaded, beside walls: no beam, 0.6 of the diffuse light, 0.7 of the walls' 30 W/m2.
    EXPECT_NEAR(value_at(sw, IntVect(2, 1, 0), 0), Real(0.6) * dif_abs + Real(0.7) * Real(30.0), tol * sw_open);
    EXPECT_NEAR(value_at(lw, IntVect(2, 1, 0), 0) - lw_open, Real(0.9) * (Real(200.0) - Real(0.4) * lw_down), tol * 400);
    // Sunlit, beside walls: the beam whole.
    EXPECT_NEAR(value_at(sw, IntVect(3, 2, 0), 0), beam_abs + Real(0.75) * dif_abs + Real(0.7) * Real(60.0), tol * sw_open);
    EXPECT_NEAR(value_at(lw, IntVect(3, 2, 0), 0) - lw_open, Real(0.9) * (Real(120.0) - Real(0.25) * lw_down), tol * 400);
    // The footprint absorbs no shortwave and its longwave is the open ground's of the sweep.
    EXPECT_EQ(value_at(sw, IntVect(0, 0, 0), 0), Real(0.0));
    EXPECT_EQ(value_at(lw, IntVect(0, 0, 0), 0), lw_open);
    EXPECT_NE(radiation.ground_walls(0), nullptr);

    // A second update after the same sweep takes the walls' terms no second time.
    const Real sw_once = value_at(sw, IntVect(3, 2, 0), 0), lw_once = value_at(lw, IntVect(3, 2, 0), 0);
    radiation.advance(0, 1, Real(10.0), Real(10.0), "post_dycore", state, nullptr, L.geom, lsm, radiation_inputs,
                      false, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0.0, false);
    EXPECT_EQ(value_at(*radiation.seb_sw_net(0), IntVect(3, 2, 0), 0), sw_once);
    EXPECT_EQ(value_at(*radiation.seb_lw_net(0), IntVect(3, 2, 0), 0), lw_once);

    // A rebuilt level drops the request with the cover.
    radiation.define_level(0, rad, RdoCp, collapse_z(L.ba), L.dm, L.ba, L.domain);
    EXPECT_FALSE(radiation.takes_ground_walls(0));
    EXPECT_EQ(radiation.ground_walls_for_step(0, 2), nullptr);

    // With the longwave off the ground takes nothing from the walls, shortwave or longwave:
    // their longwave would arrive with no sky longwave and no emission to set against it.
    RadChoice no_lw = rad;
    no_lw.lw_enabled = false;
    TwoStreamRadiation dark;
    dark.resize(1);
    dark.define_level(0, no_lw, RdoCp, collapse_z(L.ba), L.dm, L.ba, L.domain);
    dark.set_building_cover(0, blank, L.domain, true);
    EXPECT_FALSE(dark.takes_ground_walls(0));
    EXPECT_EQ(dark.ground_walls_for_step(0, 1), nullptr);
    EXPECT_NE(dark.ground_shadow_for_step(0, 1), nullptr);   // the shadow is still taken
}

// erf.ibseb.two_stream_roofs: the column under a building starts on the roof the
// immersed-boundary balance gives for it. The cell on each column's roof; a column whose
// roof is given lit, one given shaded, and one whose roof is not given (it keeps the ground).
TEST(TwoStreamBuildingCover, ColumnsStartOnTheRoofsGiven)
{
    using namespace amrex;
    Level L;
    RadChoice rad = cover_choice(true);
    rad.surface_emissivity_lw = Real(0.95);
    MultiFab blank(L.ba, L.dm, 1, 1);
    blank.setVal(0.0);
    // (0,0): a building two cells high. (2,2): one cell high. (3,3): solid above air (an
    // overhang), no building on the ground. (1,3): a building two cells high whose roof
    // this test does not give.
    for (const IntVect& p : {IntVect(0, 0, 0), IntVect(0, 0, 1), IntVect(2, 2, 0), IntVect(3, 3, 2), IntVect(1, 3, 0), IntVect(1, 3, 1)}) {
        set_at(blank, p, 1.0);
    }
    MultiFab state(L.ba, L.dm, RhoQ2_comp + 1, 0);
    state.setVal(Real(0.0));
    state.setVal(Real(1.0), Rho_comp, 1);
    state.setVal(getThgivenRandT(Real(1.0), Real(290.0), RdoCp), RhoTheta_comp, 1);
    LandSurface lsm;
    lsm.ReSize(1);
    lsm.SetModel<NullSurf>();
    const Vector<const MultiFab*> radiation_inputs {};

    // The columns as they are without the input: through the buildings to the ground.
    TwoStreamRadiation through;
    through.resize(1);
    through.define_level(0, rad, RdoCp, collapse_z(L.ba), L.dm, L.ba, L.domain);
    through.set_building_cover(0, blank, L.domain);
    EXPECT_FALSE(through.starts_on_roofs(0));
    EXPECT_EQ(through.roof_cell(0), nullptr);
    EXPECT_EQ(through.roofs_for_step(0, 1), nullptr);
    MultiFab q0(L.ba, L.dm, 2, 0), f0(convert(L.ba, IntVect(0, 0, 1)), L.dm, 4, 0);
    through.advance(0, 1, Real(0.0), Real(10.0), "pre_dycore", state, nullptr, L.geom, lsm, radiation_inputs,
                    false, &q0, &f0, nullptr, nullptr, nullptr, nullptr, nullptr, 0.0, false);

    TwoStreamRadiation radiation;
    radiation.resize(1);
    radiation.define_level(0, rad, RdoCp, collapse_z(L.ba), L.dm, L.ba, L.domain);
    radiation.set_building_cover(0, blank, L.domain, false, true);
    ASSERT_TRUE(radiation.starts_on_roofs(0));
    const iMultiFab* roof_cell = radiation.roof_cell(0);
    ASSERT_NE(roof_cell, nullptr);
    const auto cell_at = [&] (int i, int j) {
        int v = -99;
        for (MFIter mfi(*roof_cell); mfi.isValid(); ++mfi) {
            if (!mfi.validbox().contains(IntVect(i, j, 0))) { continue; }
            IArrayBox h(mfi.validbox(), 1, The_Pinned_Arena());
            Gpu::dtoh_memcpy(h.dataPtr(), (*roof_cell)[mfi].dataPtr(), h.nBytes());
            v = h.const_array()(i, j, 0);
        }
        ParallelDescriptor::ReduceIntMax(v);
        return v;
    };
    EXPECT_EQ(cell_at(0, 0), 2);   // on the two-cell building
    EXPECT_EQ(cell_at(2, 2), 1);   // on the one-cell building
    EXPECT_EQ(cell_at(1, 3), 2);
    EXPECT_EQ(cell_at(3, 3), 0);   // the overhang is no roof on the ground
    EXPECT_EQ(cell_at(1, 2), 0);   // open ground

    MultiFab* roofs = radiation.roofs_for_step(0, 1);
    ASSERT_NE(roofs, nullptr);
    ASSERT_EQ(roofs->nComp(), 5);
    roofs->setVal(0.0);
    give_roof(*roofs, IntVect(0, 0, 0), Real(0.4), Real(0.85), Real(320.0));
    give_roof(*roofs, IntVect(2, 2, 0), Real(0.1), Real(0.9), Real(305.0), Real(0.0));   // in a taller building's shadow
    MultiFab q(L.ba, L.dm, 2, 0), f(convert(L.ba, IntVect(0, 0, 1)), L.dm, 4, 0);
    radiation.advance(0, 1, Real(0.0), Real(10.0), "pre_dycore", state, nullptr, L.geom, lsm, radiation_inputs,
                      false, &q, &f, nullptr, nullptr, nullptr, nullptr, nullptr, 0.0, false);

    const Real tol = Real(256) * real_eps;   // a multiple of the working precision's roundoff
    const Real sigma = Real(5.670374419e-8);
    // (0,0): the roof at the bottom of cell 2. Nothing below it; at it, the roof's fluxes.
    for (int k = 0; k < 2; ++k) {
        for (int comp = 0; comp < 4; ++comp) { EXPECT_EQ(value_at(f, IntVect(0, 0, k), comp), Real(0.0)) << k << " " << comp; }
        EXPECT_EQ(value_at(q, IntVect(0, 0, k), 0), Real(0.0));
        EXPECT_EQ(value_at(q, IntVect(0, 0, k), 1), Real(0.0));
    }
    const Real sw_dn = value_at(f, IntVect(0, 0, 2), 1), lw_dn = value_at(f, IntVect(0, 0, 2), 3);
    ASSERT_GT(sw_dn, Real(100.0));
    EXPECT_NEAR(value_at(f, IntVect(0, 0, 2), 0), Real(0.4) * sw_dn, tol * sw_dn);
    EXPECT_NEAR(value_at(f, IntVect(0, 0, 2), 2),
                Real(0.85) * sigma * std::pow(Real(320.0), 4) + Real(0.15) * lw_dn, tol * Real(600.0));
    // (2,2): another roof, another height, shaded: it reflects the diffuse light only.
    const Real sw_dn2 = value_at(f, IntVect(2, 2, 1), 1);
    const Real beam2 = Real(1000.0) * Real(0.5) * std::exp(Real(-0.05 * 3 / 0.5));
    ASSERT_GT(sw_dn2 - beam2, Real(1.0));
    EXPECT_NEAR(value_at(f, IntVect(2, 2, 1), 0), Real(0.1) * (sw_dn2 - beam2), tol * sw_dn2);
    EXPECT_EQ(value_at(f, IntVect(2, 2, 0), 1), Real(0.0));
    // The outgoing longwave over the hot roof is not the ground's.
    EXPECT_GT(value_at(f, IntVect(0, 0, 4), 2) - value_at(f0, IntVect(0, 0, 4), 2), Real(5.0));
    // (1,3): a building whose roof was not given keeps the ground, as without the input;
    // so do the open ground and the column under the overhang, to the last bit.
    for (const IntVect& p : {IntVect(1, 3, 0), IntVect(1, 2, 0), IntVect(3, 3, 0)}) {
        for (int k = 0; k <= 4; ++k) {
            for (int comp = 0; comp < 4; ++comp) {
                EXPECT_EQ(value_at(f, IntVect(p[0], p[1], k), comp), value_at(f0, IntVect(p[0], p[1], k), comp)) << p << " " << k << " " << comp;
            }
        }
    }
    // With a cloud over half of every column (the clear and the cloudy sweep, blended): the
    // column still starts on its roof in both, so nothing lies below the roof and what goes
    // up at the roof is the roof's share of what comes down.
    RadChoice cloudy = rad;
    cloudy.tau_profile_type = TauProfileType::CloudLayer;
    cloudy.cloud_base_height_m = Real(250.0);
    cloudy.cloud_top_height_m = Real(350.0);
    cloudy.cloud_tau_per_layer = Real(2.0);
    cloudy.cloud_fraction = Real(0.5);
    TwoStreamRadiation half;
    half.resize(1);
    half.define_level(0, cloudy, RdoCp, collapse_z(L.ba), L.dm, L.ba, L.domain);
    half.set_building_cover(0, blank, L.domain, false, true);
    MultiFab* roofs_c = half.roofs_for_step(0, 1);
    ASSERT_NE(roofs_c, nullptr);
    roofs_c->setVal(0.0);
    give_roof(*roofs_c, IntVect(0, 0, 0), Real(0.4), Real(0.85), Real(320.0), Real(0.0));
    MultiFab qc(L.ba, L.dm, 2, 0), fc(convert(L.ba, IntVect(0, 0, 1)), L.dm, 4, 0);
    half.advance(0, 1, Real(0.0), Real(10.0), "pre_dycore", state, nullptr, L.geom, lsm, radiation_inputs,
                 false, &qc, &fc, nullptr, nullptr, nullptr, nullptr, nullptr, 0.0, false);
    for (int k = 0; k < 2; ++k) {
        for (int comp = 0; comp < 4; ++comp) { EXPECT_EQ(value_at(fc, IntVect(0, 0, k), comp), Real(0.0)) << k << " " << comp; }
    }
    // Shaded: the roof sends up its share of the diffuse light only, clear or cloudy, so
    // less than its share of everything that comes down, and not nothing.
    const Real up_c = value_at(fc, IntVect(0, 0, 2), 0), dn_c = value_at(fc, IntVect(0, 0, 2), 1);
    EXPECT_GT(up_c, Real(0.0));
    EXPECT_LT(up_c, Real(0.4) * dn_c - Real(1.0));
    // The cloud is there: less light reaches the roof than under the clear sky.
    EXPECT_LT(dn_c, sw_dn - Real(10.0));

    // A rebuilt level drops the request with the cover.
    radiation.define_level(0, rad, RdoCp, collapse_z(L.ba), L.dm, L.ba, L.domain);
    EXPECT_FALSE(radiation.starts_on_roofs(0));
    EXPECT_EQ(radiation.roofs_for_step(0, 2), nullptr);
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

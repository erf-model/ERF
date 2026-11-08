#include <cmath>
#include <utility>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <gtest/gtest.h>

#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_LandSurface.H>
#include <ERF_RadStruct.H>
#include <ERF_TwoStreamCanopyForcing.H>
#include <ERF_TwoStreamRadiation.H>

// Two-stream -> building faces (erf.ibseb.radiation = two_stream)
// ----------------------------------------------------------------
// What the sweep hands the faces of the immersed-boundary balance, per column:
//   1. at every interface up to the top of the canopy, the direct beam is the Beer-Lambert
//      beam there, below the sweep's shortwave down by its diffuse light;
//   2. the beam is blended clear/cloudy like every other flux, from each evaluation's own
//      beam;
//   3. with the shortwave off or the sun down there is no beam, whatever the sweep's scratch
//      holds;
//   4. the faces get it only for the step of the sweep that wrote it, and a rebuilt
//      level drops the request;
//   5. a roof samples the interface at its own height, a wall the mean of its cell's two,
//      and a face takes the beam over the cosine, the rest of the shortwave down, and the
//      longwave down and the shortwave and longwave up there, from its own column.

namespace {

amrex::BoxArray collapse_z (const amrex::BoxArray& ba)
{
    amrex::BoxList boxes = ba.boxList();
    for (auto& box : boxes) { box.setRange(2, 0); }
    return amrex::BoxArray(std::move(boxes));
}

// Component comp of mf at point, over all ranks.
amrex::Real component_at (const amrex::MultiFab& mf, const amrex::IntVect& point, int comp)
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
    amrex::Real value = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceRealSum(value);
    return value;
}

constexpr amrex::Real kS0  = 1000.0;
constexpr amrex::Real kMu  = 0.6;
constexpr amrex::Real kTau = 0.1;   // shortwave optical depth of a clear layer
constexpr int kNlev = 4;            // four 1 km layers

struct CanopyResult
{
    // Per interface m = 0 .. kNlev: the beam the sweep kept, and its SW and LW down there.
    amrex::Real beam[kNlev + 1] = {};
    amrex::Real sw_dn[kNlev + 1] = {}, lw_dn[kNlev + 1] = {};
    amrex::Real cos_zenith = 0.0;
};

// One sweep of a two-column level with the beam kept up to interface m_top, read at column
// (1, 0). The atmosphere absorbs and scatters, and the cloudy evaluation has a cloud over the
// middle two layers (centres at 1.5 and 2.5 km).
CanopyResult run_canopy (amrex::Real cloud_fraction, int m_top, bool sw_enabled)
{
    using namespace amrex;

    const Box domain(IntVect(0, 0, 0), IntVect(1, 0, kNlev - 1));
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

    RadChoice rad;
    rad.enabled = true;
    rad.sw_enabled = sw_enabled;
    rad.lw_enabled = true;
    rad.fixed_solar_zenith_angle = kMu;
    rad.fixed_total_solar_irradiance = kS0;
    rad.rad_t_sfc = Real(300.0);
    rad.tau_per_layer = kTau;
    rad.single_scattering_albedo = Real(0.5);
    rad.surface_albedo_sw = Real(0.2);
    rad.tau_profile_type = TauProfileType::CloudLayer;
    rad.cloud_base_height_m = Real(1000.0);
    rad.cloud_top_height_m = Real(3000.0);
    rad.cloud_tau_per_layer = Real(2.0);
    rad.cloud_fraction = cloud_fraction;

    TwoStreamRadiation radiation;
    radiation.resize(1);
    radiation.define_level(0, rad, RdoCp, ba2d, dm, ba, domain);
    radiation.supply_canopy_forcing(0, m_top);

    LandSurface lsm;
    lsm.ReSize(1);
    lsm.SetModel<NullSurf>();
    MultiFab qheating(ba, dm, 2, 0);
    MultiFab rad_fluxes(convert(ba, IntVect(0, 0, 1)), dm, 4, 0);
    const Vector<const MultiFab*> radiation_inputs {};
    radiation.advance(0, 7, Real(0.0), Real(10.0), "pre_dycore",
                      state, nullptr, geom, lsm, radiation_inputs, false,
                      &qheating, &rad_fluxes, nullptr, nullptr, nullptr, nullptr, nullptr,
                      0.0, false);

    CanopyResult r;
    const TwoStreamCanopyView view = radiation.canopy_forcing(0, 7);
    EXPECT_NE(view.beam, nullptr);
    EXPECT_NE(view.cos_zenith, nullptr);
    EXPECT_EQ(view.fluxes, nullptr);   // the caller's rad_fluxes, which the model does not own
    if (view.beam == nullptr || view.cos_zenith == nullptr) { return r; }
    for (int m = 0; m <= m_top; ++m) {
        r.beam[m]  = component_at(*view.beam, IntVect(1, 0, m), 0);
        r.sw_dn[m] = component_at(rad_fluxes, IntVect(1, 0, m), 1);
        r.lw_dn[m] = component_at(rad_fluxes, IntVect(1, 0, m), 3);
    }
    r.cos_zenith = component_at(*view.cos_zenith, IntVect(1, 0, 0), 0);
    return r;
}

amrex::Real rel_tol () { return sizeof(amrex::Real) == 8 ? amrex::Real(1.e-12) : amrex::Real(1.e-5); }

} // namespace

// At every interface up to the top of the canopy the beam the sweep keeps is Beer-Lambert's
// through the layers above it (the optical depth counts scattering too), below the sweep's
// total shortwave down there by the diffuse light the scattering made; and the cosine of the
// zenith is the sweep's.
TEST(TwoStreamCanopyForcing, BeamAtEveryInterface)
{
    const amrex::Real tol = rel_tol();
    const CanopyResult r = run_canopy(amrex::Real(0.0), kNlev, true);
    for (int m = 0; m <= kNlev; ++m) {
        const amrex::Real beam = kS0 * kMu * std::exp(-kTau * (kNlev - m) / kMu);
        EXPECT_NEAR(r.beam[m], beam, tol * beam) << "interface " << m;
        if (m < kNlev) {
            EXPECT_GT(r.sw_dn[m] - r.beam[m], amrex::Real(1.0)) << "interface " << m;
        } else {
            EXPECT_NEAR(r.sw_dn[m], r.beam[m], tol * kS0);   // nothing diffuse comes in from space
        }
    }
    EXPECT_EQ(r.cos_zenith, kMu);
    // A canopy that stops lower keeps the same beam up to its top.
    const CanopyResult low = run_canopy(amrex::Real(0.0), 2, true);
    for (int m = 0; m <= 2; ++m) { EXPECT_EQ(low.beam[m], r.beam[m]) << "interface " << m; }
}

// The beam is blended clear/cloudy with the cloud fraction, from each evaluation's own
// beam: under the cloud (at 1 km, below it) the cloudy beam is far weaker.
TEST(TwoStreamCanopyForcing, BlendsTheBeamOfTheCloudyColumn)
{
    const amrex::Real tol = rel_tol();
    const CanopyResult clear  = run_canopy(amrex::Real(0.0), kNlev, true);
    const CanopyResult cloudy = run_canopy(amrex::Real(1.0), kNlev, true);
    const CanopyResult part   = run_canopy(amrex::Real(0.4), kNlev, true);
    EXPECT_LT(cloudy.beam[1], amrex::Real(0.1) * clear.beam[1]);
    for (int m = 0; m <= kNlev; ++m) {
        const amrex::Real blend = amrex::Real(0.6) * clear.beam[m] + amrex::Real(0.4) * cloudy.beam[m];
        EXPECT_NEAR(part.beam[m], blend, tol * blend) << "interface " << m;
        EXPECT_LE(part.beam[m], part.sw_dn[m] * (amrex::Real(1.0) + tol)) << "interface " << m;
    }
}

// With the shortwave off the sweep writes no beam to its scratch; the faces must see none.
TEST(TwoStreamCanopyForcing, NoBeamWithoutTheShortwave)
{
    const CanopyResult r = run_canopy(amrex::Real(0.4), kNlev, false);
    for (int m = 0; m <= kNlev; ++m) { EXPECT_EQ(r.beam[m], amrex::Real(0.0)) << "interface " << m; }
    EXPECT_GT(r.lw_dn[1], amrex::Real(1.0));
}

// The beam helper itself, on a scratch full of leftovers: nothing with the shortwave off or
// the sun down; otherwise the clear evaluation's beam as it is, and the cloudy one blended
// in with the cloud fraction.
TEST(TwoStreamCanopyForcing, BeamIgnoresTheScratchWhenDark)
{
    using namespace amrex;
    const Box bx(IntVect(0, 0, 0), IntVect(0, 0, 3));
    FArrayBox scratch(bx, 2, The_Cpu_Arena());
    scratch.setVal<RunOn::Host>(Real(7.0));
    FArrayBox out(bx, 1, The_Cpu_Arena());
    out.setVal<RunOn::Host>(Real(-1.0));
    const auto s = scratch.const_array();
    const auto o = out.array();
    two_stream_canopy_beam(0, 0, 0, 3, s, 1, false, Real(0.5), Real(0.0), false, o);
    for (int m = 0; m <= 3; ++m) { EXPECT_EQ(o(0, 0, m), Real(0.0)); }
    two_stream_canopy_beam(0, 0, 0, 3, s, 1, true, Real(-0.1), Real(0.0), false, o);
    for (int m = 0; m <= 3; ++m) { EXPECT_EQ(o(0, 0, m), Real(0.0)); }
    two_stream_canopy_beam(0, 0, 0, 3, s, 1, true, Real(0.5), Real(0.25), false, o);
    for (int m = 0; m <= 3; ++m) { EXPECT_EQ(o(0, 0, m), Real(7.0)); }
    scratch.setVal<RunOn::Host>(Real(3.0));
    two_stream_canopy_beam(0, 0, 0, 3, s, 1, true, Real(0.5), Real(0.25), true, o);
    for (int m = 0; m <= 3; ++m) { EXPECT_EQ(o(0, 0, m), Real(0.75) * Real(7.0) + Real(0.25) * Real(3.0)); }
}

// A roof (z face) samples the interface at its fluid cell's bottom, its own height; a wall
// the mean of its cell's two interfaces; the index counts from the level's domain bottom.
TEST(TwoStreamCanopyForcing, FacesSampleTheirOwnHeight)
{
    amrex::Real w_up = -1.0;
    EXPECT_EQ(two_stream_canopy_sample(5, 0, 2, -1, w_up), 5);   // roof: solid below
    EXPECT_EQ(w_up, amrex::Real(0.0));
    EXPECT_EQ(two_stream_canopy_sample(5, 0, 2, +1, w_up), 6);   // ceiling: solid above
    EXPECT_EQ(w_up, amrex::Real(0.0));
    EXPECT_EQ(two_stream_canopy_sample(5, 0, 0, -1, w_up), 5);
    EXPECT_EQ(w_up, amrex::Real(0.5));
    EXPECT_EQ(two_stream_canopy_sample(5, 0, 0, +1, w_up), 5);   // a wall's side does not move it
    EXPECT_EQ(two_stream_canopy_sample(5, 2, 1, -1, w_up), 3);
    EXPECT_EQ(w_up, amrex::Real(0.5));
}

// What a face takes from its column, on fields whose every interface differs: at its sample
// (a wall the mean of its cell's two interfaces) the beam over the cosine, the shortwave down
// less the beam, the longwave down and the shortwave and longwave up. The level's bottom is
// k0 = 2 here, so an index counted from 0 would read the wrong cells. No beam at night.
TEST(TwoStreamCanopyForcing, FaceSkyFromItsColumn)
{
    using namespace amrex;
    const int k0 = 2;
    const Box bm(IntVect(0, 0, 0), IntVect(0, 0, 5));          // beam: interfaces m = 0..5
    const Box bf(IntVect(0, 0, k0), IntVect(0, 0, k0 + 5));    // fluxes: absolute k
    const Box b2(IntVect(0, 0, 0), IntVect(0, 0, 0));
    FArrayBox beam(bm, 1, The_Cpu_Arena()), flux(bf, 4, The_Cpu_Arena()), cz(b2, 1, The_Cpu_Arena());
    const auto bA = beam.array(); const auto fA = flux.array(); const auto cA = cz.array();
    for (int m = 0; m <= 5; ++m) {
        bA(0, 0, m) = Real(100.0 + 10.0 * m);
        fA(0, 0, k0 + m, 0) = Real(50.0 + m);         // SW up
        fA(0, 0, k0 + m, 1) = Real(130.0 + 13.0 * m);  // SW down: beam + 30 + 3 m
        fA(0, 0, k0 + m, 2) = Real(400.0 + 2.0 * m);   // LW up
        fA(0, 0, k0 + m, 3) = Real(300.0 - 7.0 * m);   // LW down
    }
    cA(0, 0, 0) = Real(0.5);
    const Real tol = rel_tol();
    // A roof on interface 3 (fluid cell k0 + 3).
    TwoStreamFaceSky r = two_stream_face_sky(0, 0, k0 + 3, k0, 2, -1, beam.const_array(), cz.const_array(), flux.const_array());
    EXPECT_NEAR(r.dni, Real(130.0 / 0.5), tol * Real(260.0));
    EXPECT_NEAR(r.diffuse, Real(39.0), tol * Real(100.0));
    EXPECT_NEAR(r.lw_down, Real(279.0), tol * Real(300.0));
    EXPECT_NEAR(r.sw_up, Real(53.0), tol * Real(100.0));
    EXPECT_NEAR(r.lw_up, Real(406.0), tol * Real(400.0));
    // A wall in cell k0 + 3: the mean of interfaces 3 and 4.
    TwoStreamFaceSky w = two_stream_face_sky(0, 0, k0 + 3, k0, 0, -1, beam.const_array(), cz.const_array(), flux.const_array());
    EXPECT_NEAR(w.dni, Real(135.0 / 0.5), tol * Real(270.0));
    EXPECT_NEAR(w.diffuse, Real(40.5), tol * Real(100.0));
    EXPECT_NEAR(w.lw_down, Real(275.5), tol * Real(300.0));
    EXPECT_NEAR(w.sw_up, Real(53.5), tol * Real(100.0));
    EXPECT_NEAR(w.lw_up, Real(407.0), tol * Real(400.0));
    // A ceiling in cell k0 + 3 (solid above): interface 4.
    TwoStreamFaceSky c = two_stream_face_sky(0, 0, k0 + 3, k0, 2, +1, beam.const_array(), cz.const_array(), flux.const_array());
    EXPECT_NEAR(c.dni, Real(140.0 / 0.5), tol * Real(280.0));
    EXPECT_NEAR(c.sw_up, Real(54.0), tol * Real(100.0));
    EXPECT_NEAR(c.lw_up, Real(408.0), tol * Real(400.0));
    // Night: the cosine floored at zero, no beam, whatever the field holds.
    cA(0, 0, 0) = Real(0.0);
    TwoStreamFaceSky n = two_stream_face_sky(0, 0, k0 + 3, k0, 2, -1, beam.const_array(), cz.const_array(), flux.const_array());
    EXPECT_EQ(n.dni, Real(0.0));
}

// The faces get the forcing for the step of the sweep that wrote it and no other, nothing
// before a sweep or without the request, and a rebuilt level drops the request.
TEST(TwoStreamCanopyForcing, OnlyForTheStepOfItsSweep)
{
    using namespace amrex;
    const Box domain(IntVect(0, 0, 0), IntVect(1, 0, kNlev - 1));
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

    RadChoice rad;
    rad.enabled = true;
    rad.fixed_solar_zenith_angle = kMu;
    rad.fixed_total_solar_irradiance = kS0;
    rad.rad_t_sfc = Real(300.0);

    TwoStreamRadiation radiation;
    radiation.resize(1);
    radiation.define_level(0, rad, RdoCp, ba2d, dm, ba, domain);
    EXPECT_FALSE(radiation.supplies_canopy_forcing(0));
    EXPECT_EQ(radiation.canopy_forcing(0, 3).beam, nullptr);

    radiation.supply_canopy_forcing(0, 2);
    EXPECT_TRUE(radiation.supplies_canopy_forcing(0));
    EXPECT_EQ(radiation.canopy_forcing(0, 3).beam, nullptr);   // no sweep yet

    LandSurface lsm;
    lsm.ReSize(1);
    lsm.SetModel<NullSurf>();
    MultiFab qheating(ba, dm, 2, 0);
    MultiFab rad_fluxes(convert(ba, IntVect(0, 0, 1)), dm, 4, 0);
    const Vector<const MultiFab*> radiation_inputs {};
    radiation.advance(0, 3, Real(0.0), Real(10.0), "pre_dycore",
                      state, nullptr, geom, lsm, radiation_inputs, false,
                      &qheating, &rad_fluxes, nullptr, nullptr, nullptr, nullptr, nullptr,
                      0.0, false);
    EXPECT_NE(radiation.canopy_forcing(0, 3).beam, nullptr);
    EXPECT_EQ(radiation.canopy_forcing(0, 4).beam, nullptr);   // the next step needs its own sweep
    EXPECT_EQ(radiation.canopy_forcing(0, 2).beam, nullptr);
    // The post-dycore call does not sweep, so it neither moves nor renews the stamp.
    radiation.advance(0, 4, Real(10.0), Real(10.0), "post_dycore",
                      state, nullptr, geom, lsm, radiation_inputs, false,
                      &qheating, &rad_fluxes, nullptr, nullptr, nullptr, nullptr, nullptr,
                      0.0, false);
    EXPECT_EQ(radiation.canopy_forcing(0, 4).beam, nullptr);
    EXPECT_NE(radiation.canopy_forcing(0, 3).beam, nullptr);

    radiation.define_level(0, rad, RdoCp, ba2d, dm, ba, domain);
    EXPECT_FALSE(radiation.supplies_canopy_forcing(0));
    EXPECT_EQ(radiation.canopy_forcing(0, 3).beam, nullptr);
}

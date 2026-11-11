#include <cmath>
#include <vector>

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_NumericalConstants.H>
#include <ERF_RadStruct.H>
#include <ERF_TwoStreamColumn.H>
#include <ERF_TwoStreamCanopyForcing.H>

// A roof as the lower boundary of a two-stream column (erf.ibseb.two_stream_roofs,
// TwoStreamRoof): the column under a building starts on the building's roof.
//  - A column of eight layers on a roof at the bottom of layer 3 is, above the roof, the
//    column of the five layers 3 .. 7 on a ground with the roof's albedo, emissivity and
//    temperature: every flux, heating rate and scalar, to the last bit. Below the roof
//    there is no flux and no heating, and the beam kept in the scratch is zero.
//  - At the roof the fluxes are the roof's: shortwave up = albedo x shortwave down,
//    longwave up = emissivity sigma T^4 + (1 - emissivity) longwave down, and the absorbed
//    shortwave is the rest.
//  - Heights are still measured from the ground: a cloud band at a given height takes the
//    same layers with a roof as without one.
//  - A roof at the ground is the ground; a building that fills the column leaves nothing.

namespace {

constexpr amrex::Real kDz = 100.0;

struct Column
{
    int k_lo = 0, k_hi = 0;
    std::vector<amrex::Real> q_sw, q_lw;                  // layers k_lo .. k_hi
    std::vector<amrex::Real> sw_up, sw_dn, lw_up, lw_dn;  // interfaces k_lo .. k_hi + 1
    std::vector<amrex::Real> beam;                        // the scratch's beam, interfaces k_lo .. k_hi + 1
    amrex::Real max_heating = 0.0, sw_surface = 0.0, sw_up_toa = 0.0, lw_net_surface = 0.0, lw_up_toa = 0.0;
    amrex::Real sw_down_toa = 0.0, sw_beam_absorbed = 0.0;
    amrex::Real at (const std::vector<amrex::Real>& v, int k) const { return v[static_cast<size_t>(k - k_lo)]; }
};

RadChoice roof_choice ()
{
    RadChoice rad;
    rad.enabled = true;
    rad.sw_enabled = true;
    rad.lw_enabled = true;
    rad.fixed_solar_zenith_angle = amrex::Real(0.6);
    rad.fixed_total_solar_irradiance = amrex::Real(1000.0);
    rad.tau_per_layer = amrex::Real(0.05);
    rad.tau_lw_per_layer = amrex::Real(0.2);
    rad.single_scattering_albedo = amrex::Real(0.5);
    rad.surface_albedo_sw = amrex::Real(0.2);
    rad.surface_emissivity_lw = amrex::Real(0.95);
    rad.rad_t_sfc = amrex::Real(300.0);
    return rad;
}

// Sweep the column of layers k_lo .. k_hi, on a roof when given, over a ground at t_ground.
Column sweep (const RadChoice& rad, int k_lo, int k_hi, const TwoStreamRoof* roof, amrex::Real t_ground, bool cloudy = false)
{
    const TwoStreamParams params = make_two_stream_params(rad, RdoCp);
    const amrex::Box bx(amrex::IntVect(0, 0, k_lo), amrex::IntVect(0, 0, k_hi));
    const amrex::Box xy(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    amrex::FArrayBox state(bx, RhoQ2_comp + 1);
    amrex::FArrayBox qheating(bx, 2);
    amrex::FArrayBox fluxes(two_stream_scratch_box(bx), 4);
    amrex::FArrayBox scratch(two_stream_scratch_box(bx), TwoStreamScratch::NCOMP);
    amrex::FArrayBox t_sfc(xy, 1);
    state.setVal<amrex::RunOn::Device>(0.0);
    qheating.setVal<amrex::RunOn::Device>(-7.0);
    fluxes.setVal<amrex::RunOn::Device>(-7.0);
    scratch.setVal<amrex::RunOn::Device>(-7.0);
    t_sfc.setVal<amrex::RunOn::Device>(t_ground);
    {
        const auto s = state.array();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const amrex::Real rho = amrex::Real(1.1) - amrex::Real(0.01) * k;
            // Another temperature in every layer: a sweep that took another layer's would give another number.
            const amrex::Real T = amrex::Real(285.0) - amrex::Real(2.5) * k + amrex::Real(0.3) * k * k;
            s(i, j, k, Rho_comp) = rho;
            s(i, j, k, RhoTheta_comp) = rho * getThgivenRandT(rho, T, RdoCp, amrex::Real(0.0));
        });
    }
    amrex::Gpu::DeviceVector<amrex::Real> scalars(7, 0.0);
    amrex::Real* out = scalars.data();
    const auto state_arr = state.const_array();
    const auto q_arr = qheating.array();
    const auto flux_arr = fluxes.array();
    const auto scratch_arr = scratch.array();
    const auto t_sfc_arr = t_sfc.const_array();
    const bool with_roof = (roof != nullptr);
    const TwoStreamRoof roof_copy = with_roof ? *roof : TwoStreamRoof{};
    const amrex::Real dz = kDz;
    amrex::ParallelFor(xy, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
        amrex::Real max_heating = 0.0, sw_surface = 0.0, sw_up = 0.0, lw_net = 0.0, lw_up = 0.0, sw_toa = 0.0, beam_abs = 0.0;
        const amrex::Array4<const amrex::Real> no_z;
        const TwoStreamRoof r = roof_copy;
        vertical_two_stream_sweep(i, j, bx, dz, state_arr, params, cloudy, q_arr,
                                  max_heating, sw_surface, sw_up, lw_net, lw_up, sw_toa, no_z, scratch_arr,
                                  false, nullptr, false, nullptr, false, nullptr, true, &t_sfc_arr,
                                  false, nullptr, false, nullptr, nullptr, &flux_arr, nullptr, &beam_abs,
                                  with_roof ? &r : nullptr);
        out[0] = max_heating; out[1] = sw_surface; out[2] = sw_up; out[3] = lw_net; out[4] = lw_up; out[5] = sw_toa;
        out[6] = beam_abs;
    });
    amrex::Gpu::streamSynchronize();
    std::vector<amrex::Real> h(7);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, scalars.begin(), scalars.end(), h.begin());
    amrex::FArrayBox hq(bx, 2, amrex::The_Pinned_Arena());
    amrex::FArrayBox hf(two_stream_scratch_box(bx), 4, amrex::The_Pinned_Arena());
    amrex::FArrayBox hs(two_stream_scratch_box(bx), TwoStreamScratch::NCOMP, amrex::The_Pinned_Arena());
    hq.copy<amrex::RunOn::Device>(qheating);
    hf.copy<amrex::RunOn::Device>(fluxes);
    hs.copy<amrex::RunOn::Device>(scratch);
    amrex::Gpu::streamSynchronize();
    Column c;
    c.k_lo = k_lo; c.k_hi = k_hi;
    c.max_heating = h[0]; c.sw_surface = h[1]; c.sw_up_toa = h[2]; c.lw_net_surface = h[3]; c.lw_up_toa = h[4];
    c.sw_down_toa = h[5]; c.sw_beam_absorbed = h[6];
    const auto q = hq.const_array();
    const auto f = hf.const_array();
    const auto s = hs.const_array();
    for (int k = k_lo; k <= k_hi; ++k) { c.q_sw.push_back(q(0, 0, k, 0)); c.q_lw.push_back(q(0, 0, k, 1)); }
    for (int k = k_lo; k <= k_hi + 1; ++k) {
        c.sw_up.push_back(f(0, 0, k, 0)); c.sw_dn.push_back(f(0, 0, k, 1));
        c.lw_up.push_back(f(0, 0, k, 2)); c.lw_dn.push_back(f(0, 0, k, 3));
        c.beam.push_back(s(0, 0, k, TwoStreamScratch::F_DIR));
    }
    return c;
}

// The reflected beam's share of the light going up at a column's lower boundary
// (two_stream_ground_beam_share()), for a column whose boundary is at the interface `base`:
// the beam is 100 (m + 1) at interface m, the light going up 200 + 100 k at interface k, the
// boundary's albedo 0.25.
amrex::Real beam_share_at (int base)
{
    const amrex::Box col(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 4));
    const amrex::Box xy(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    amrex::FArrayBox beam(col, 1), fluxes(col, 4), alb(xy, 2);
    alb.setVal<amrex::RunOn::Device>(amrex::Real(0.25), xy, 0, 1);
    alb.setVal<amrex::RunOn::Device>(static_cast<amrex::Real>(base), xy, 1, 1);
    {
        const auto b = beam.array();
        const auto f = fluxes.array();
        amrex::ParallelFor(col, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            b(i, j, k) = amrex::Real(100.0) * (k + 1);
            f(i, j, k, 0) = amrex::Real(200.0) + amrex::Real(100.0) * k;
            f(i, j, k, 1) = f(i, j, k, 2) = f(i, j, k, 3) = amrex::Real(0.0);
        });
    }
    amrex::Gpu::DeviceVector<amrex::Real> out(1, 0.0);
    amrex::Real* po = out.data();
    const auto b = beam.const_array();
    const auto f = fluxes.const_array();
    const auto a = alb.const_array();
    amrex::ParallelFor(xy, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
        po[0] = two_stream_ground_beam_share(i, j, 0, b, f, a);
    });
    amrex::Gpu::streamSynchronize();
    std::vector<amrex::Real> h(1);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, out.begin(), out.end(), h.begin());
    return h[0];
}

} // namespace

// erf.ibseb.ground_view_shadow with roofs: a wall above a lower roof takes the beam's share
// of the light going up at the roof its column starts on, not at the ground.
TEST(TwoStreamRoof, BeamShareAtTheColumnsLowerBoundary)
{
    using amrex::Real;
    const Real tol = Real(256) * real_eps;   // a multiple of the working precision's roundoff
    // At the ground: 0.25 x 100 / 200. On a roof at interface 2: 0.25 x 300 / 400.
    EXPECT_NEAR(beam_share_at(0), Real(0.125), tol);
    EXPECT_NEAR(beam_share_at(2), Real(0.1875), tol);
    EXPECT_NEAR(beam_share_at(4), Real(0.25) * Real(500.0) / Real(600.0), tol);
}

TEST(TwoStreamRoof, ColumnAboveARoofIsTheShorterColumn)
{
    using amrex::Real;
    TwoStreamRoof roof;
    roof.k_lower = 3; roof.albedo = Real(0.35); roof.emissivity = Real(0.9); roof.temperature = Real(312.0);
    const Column on_roof = sweep(roof_choice(), 0, 7, &roof, Real(300.0));
    // The same five layers over a ground with the roof's properties.
    RadChoice as_ground = roof_choice();
    as_ground.surface_albedo_sw = roof.albedo;
    as_ground.surface_emissivity_lw = roof.emissivity;
    const Column shorter = sweep(as_ground, 3, 7, nullptr, roof.temperature);
    for (int k = 3; k <= 7; ++k) {
        EXPECT_EQ(on_roof.at(on_roof.q_sw, k), shorter.at(shorter.q_sw, k)) << "layer " << k;
        EXPECT_EQ(on_roof.at(on_roof.q_lw, k), shorter.at(shorter.q_lw, k)) << "layer " << k;
        EXPECT_NE(on_roof.at(on_roof.q_sw, k), Real(0.0)) << "layer " << k;
    }
    for (int k = 3; k <= 8; ++k) {
        EXPECT_EQ(on_roof.at(on_roof.sw_up, k), shorter.at(shorter.sw_up, k)) << "interface " << k;
        EXPECT_EQ(on_roof.at(on_roof.sw_dn, k), shorter.at(shorter.sw_dn, k)) << "interface " << k;
        EXPECT_EQ(on_roof.at(on_roof.lw_up, k), shorter.at(shorter.lw_up, k)) << "interface " << k;
        EXPECT_EQ(on_roof.at(on_roof.lw_dn, k), shorter.at(shorter.lw_dn, k)) << "interface " << k;
        EXPECT_EQ(on_roof.at(on_roof.beam, k), shorter.at(shorter.beam, k)) << "interface " << k;
    }
    EXPECT_EQ(on_roof.sw_surface, shorter.sw_surface);
    EXPECT_EQ(on_roof.sw_up_toa, shorter.sw_up_toa);
    EXPECT_EQ(on_roof.lw_net_surface, shorter.lw_net_surface);
    EXPECT_EQ(on_roof.lw_up_toa, shorter.lw_up_toa);
    EXPECT_EQ(on_roof.sw_down_toa, shorter.sw_down_toa);
    EXPECT_EQ(on_roof.sw_beam_absorbed, shorter.sw_beam_absorbed);
    EXPECT_EQ(on_roof.max_heating, shorter.max_heating);
    // Below the roof: no heating, no flux, no beam.
    for (int k = 0; k < 3; ++k) {
        EXPECT_EQ(on_roof.at(on_roof.q_sw, k), Real(0.0)) << "layer " << k;
        EXPECT_EQ(on_roof.at(on_roof.q_lw, k), Real(0.0)) << "layer " << k;
        EXPECT_EQ(on_roof.at(on_roof.sw_up, k), Real(0.0)) << "interface " << k;
        EXPECT_EQ(on_roof.at(on_roof.sw_dn, k), Real(0.0)) << "interface " << k;
        EXPECT_EQ(on_roof.at(on_roof.lw_up, k), Real(0.0)) << "interface " << k;
        EXPECT_EQ(on_roof.at(on_roof.lw_dn, k), Real(0.0)) << "interface " << k;
        EXPECT_EQ(on_roof.at(on_roof.beam, k), Real(0.0)) << "interface " << k;
    }
    // The column does differ from the one that passes through the building to the ground.
    const Column through = sweep(roof_choice(), 0, 7, nullptr, Real(300.0));
    EXPECT_GT(std::abs(on_roof.lw_up_toa - through.lw_up_toa), Real(1.0));
    EXPECT_GT(std::abs(on_roof.sw_up_toa - through.sw_up_toa), Real(1.0));
    EXPECT_NE(through.at(through.q_sw, 1), Real(0.0));
}

TEST(TwoStreamRoof, FluxesAtTheRoofAreTheRoofs)
{
    using amrex::Real;
    TwoStreamRoof roof;
    roof.k_lower = 2; roof.albedo = Real(0.4); roof.emissivity = Real(0.85); roof.temperature = Real(318.0);
    const Column c = sweep(roof_choice(), 0, 7, &roof, Real(300.0));
    const Real tol = Real(256) * real_eps;   // a multiple of the working precision's roundoff
    const Real sw_dn = c.at(c.sw_dn, 2), sw_up = c.at(c.sw_up, 2), lw_dn = c.at(c.lw_dn, 2), lw_up = c.at(c.lw_up, 2);
    ASSERT_GT(sw_dn, Real(100.0));
    ASSERT_GT(lw_dn, Real(10.0));
    // One albedo for the beam and the diffuse light: what goes up is that share of what comes down.
    EXPECT_NEAR(sw_up, roof.albedo * sw_dn, tol * sw_dn);
    EXPECT_NEAR(c.sw_surface, (Real(1.0) - roof.albedo) * sw_dn, tol * sw_dn);
    EXPECT_NEAR(c.sw_beam_absorbed, (Real(1.0) - roof.albedo) * c.at(c.beam, 2), tol * sw_dn);
    const Real sigma = stefan_boltzmann;
    const Real emitted = roof.emissivity * sigma * roof.temperature * roof.temperature * roof.temperature * roof.temperature;
    EXPECT_NEAR(lw_up, emitted + (Real(1.0) - roof.emissivity) * lw_dn, tol * lw_up);
    EXPECT_NEAR(c.lw_net_surface, lw_up - lw_dn, tol * lw_up);
    // A roof a taller building shades reflects and absorbs no beam: only the diffuse light.
    TwoStreamRoof shaded = roof;
    shaded.lit = Real(0.0);
    const Column s = sweep(roof_choice(), 0, 7, &shaded, Real(300.0));
    const Real beam = s.at(s.beam, 2);
    const Real dif_dn = s.at(s.sw_dn, 2) - beam;
    ASSERT_GT(beam, Real(100.0));
    ASSERT_GT(dif_dn, Real(1.0));
    EXPECT_EQ(beam, c.at(c.beam, 2));                                   // the beam still arrives
    EXPECT_NEAR(s.at(s.sw_up, 2), roof.albedo * dif_dn, tol * sw_dn);
    EXPECT_NEAR(s.sw_surface, (Real(1.0) - roof.albedo) * dif_dn, tol * sw_dn);
    EXPECT_EQ(s.sw_beam_absorbed, Real(0.0));
    EXPECT_LT(dif_dn, c.at(c.sw_dn, 2) - c.at(c.beam, 2));              // less of its own reflection comes back
    EXPECT_EQ(s.at(s.lw_up, 2), c.at(c.lw_up, 2));                      // the longwave does not know the shadow
    // Half lit: half the beam.
    TwoStreamRoof half = roof;
    half.lit = Real(0.5);
    const Column h = sweep(roof_choice(), 0, 7, &half, Real(300.0));
    EXPECT_NEAR(h.sw_beam_absorbed, Real(0.5) * c.sw_beam_absorbed, tol * sw_dn);
    // The ground's temperature, albedo and emissivity play no part: another ground, the same column.
    RadChoice other = roof_choice();
    other.surface_albedo_sw = Real(0.6);
    other.surface_albedo_sw_diffuse = Real(0.75);   // a ground with its own albedo for diffuse light
    other.surface_emissivity_lw = Real(0.5);
    const Column d = sweep(other, 0, 7, &roof, Real(250.0));
    EXPECT_EQ(d.lw_up_toa, c.lw_up_toa);
    EXPECT_EQ(d.sw_up_toa, c.sw_up_toa);
    EXPECT_EQ(d.sw_surface, c.sw_surface);
    EXPECT_EQ(d.at(d.sw_up, 2), c.at(c.sw_up, 2));
    EXPECT_EQ(d.at(d.sw_dn, 2), c.at(c.sw_dn, 2));
}

TEST(TwoStreamRoof, HeightsAreStillMeasuredFromTheGround)
{
    using amrex::Real;
    // A cloud between 450 and 650 m: layers 4, 5 and 6 of 100 m (centres at 450, 550, 650 m).
    RadChoice rad = roof_choice();
    rad.tau_profile_type = TauProfileType::CloudLayer;
    rad.cloud_base_height_m = Real(400.0);
    rad.cloud_top_height_m = Real(700.0);
    rad.cloud_tau_per_layer = Real(2.0);
    TwoStreamRoof roof;
    roof.k_lower = 3; roof.albedo = Real(0.35); roof.emissivity = Real(0.9); roof.temperature = Real(312.0);
    const Column on_roof = sweep(rad, 0, 7, &roof, Real(300.0), true);
    const Column through = sweep(rad, 0, 7, nullptr, Real(300.0), true);
    // The beam and the longwave coming down depend only on the layers above: the same at
    // and above the roof, with the cloud in the same layers.
    for (int k = 3; k <= 8; ++k) {
        EXPECT_EQ(on_roof.at(on_roof.beam, k), through.at(through.beam, k)) << "interface " << k;
        EXPECT_EQ(on_roof.at(on_roof.lw_dn, k), through.at(through.lw_dn, k)) << "interface " << k;
    }
    // Not vacuous: the cloud takes most of the beam between interfaces 7 and 4.
    EXPECT_LT(on_roof.at(on_roof.beam, 4), Real(0.1) * on_roof.at(on_roof.beam, 7));
    // A column that measured heights from the roof would put the cloud three layers higher.
    const Column shifted = sweep(rad, 3, 7, nullptr, Real(300.0), true);
    EXPECT_GT(std::abs(shifted.at(shifted.beam, 4) - on_roof.at(on_roof.beam, 4)), Real(1.0));
}

TEST(TwoStreamRoof, RoofAtTheGroundAndBuildingToTheTop)
{
    using amrex::Real;
    const Column plain = sweep(roof_choice(), 0, 7, nullptr, Real(300.0));
    // A roof "on" the bottom cell is no roof: the ground's own properties, to the last bit.
    TwoStreamRoof at_ground;
    at_ground.k_lower = 0; at_ground.albedo = Real(0.9); at_ground.emissivity = Real(0.1); at_ground.temperature = Real(400.0);
    const Column same = sweep(roof_choice(), 0, 7, &at_ground, Real(300.0));
    for (int k = 0; k <= 8; ++k) {
        EXPECT_EQ(same.at(same.sw_up, k), plain.at(plain.sw_up, k)) << "interface " << k;
        EXPECT_EQ(same.at(same.lw_up, k), plain.at(plain.lw_up, k)) << "interface " << k;
    }
    for (int k = 0; k <= 7; ++k) {
        EXPECT_EQ(same.at(same.q_sw, k), plain.at(plain.q_sw, k)) << "layer " << k;
        EXPECT_EQ(same.at(same.q_lw, k), plain.at(plain.q_lw, k)) << "layer " << k;
    }
    // A building that fills the column: nothing is swept, every output is zero.
    TwoStreamRoof full;
    full.k_lower = 8; full.albedo = Real(0.3); full.emissivity = Real(0.9); full.temperature = Real(310.0);
    const Column none = sweep(roof_choice(), 0, 7, &full, Real(300.0));
    for (int k = 0; k <= 7; ++k) {
        EXPECT_EQ(none.at(none.q_sw, k), Real(0.0));
        EXPECT_EQ(none.at(none.q_lw, k), Real(0.0));
    }
    for (int k = 0; k <= 8; ++k) {
        EXPECT_EQ(none.at(none.sw_dn, k), Real(0.0));
        EXPECT_EQ(none.at(none.lw_up, k), Real(0.0));
        EXPECT_EQ(none.at(none.beam, k), Real(0.0));
    }
    EXPECT_EQ(none.sw_surface, Real(0.0));
    EXPECT_EQ(none.lw_up_toa, Real(0.0));
    EXPECT_EQ(none.max_heating, Real(0.0));
    EXPECT_GT(none.sw_down_toa, Real(0.0));   // the sun still shines on the top
    // One layer above the roof: the top layer alone.
    TwoStreamRoof under_top;
    under_top.k_lower = 7; under_top.albedo = Real(0.3); under_top.emissivity = Real(0.9); under_top.temperature = Real(310.0);
    const Column top_layer = sweep(roof_choice(), 0, 7, &under_top, Real(300.0));
    EXPECT_NE(top_layer.at(top_layer.q_sw, 7), Real(0.0));
    EXPECT_EQ(top_layer.at(top_layer.q_sw, 6), Real(0.0));
    EXPECT_GT(top_layer.at(top_layer.lw_up, 8), Real(100.0));
}

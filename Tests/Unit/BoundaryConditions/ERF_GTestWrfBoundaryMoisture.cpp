#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>

#include "ERF_Utils.H"
#ifdef ERF_USE_NETCDF
#include "ERF_ReadFromWRFBdy.H"
#endif

TEST(WrfBoundaryMoisture, SerializedIndicesPreserveLegacyLayout)
{
    EXPECT_EQ(WRFBdyVars::PH, 5);
    EXPECT_EQ(WRFBdyVars::MU, 6);
    EXPECT_EQ(WRFBdyVars::PC, 7);
    EXPECT_EQ(WRFBdyVars::LegacyNumTypes, 8);
    EXPECT_EQ(WRFBdyVars::QC, 8);
    EXPECT_EQ(WRFBdyVars::QI, 9);
}

#ifdef ERF_USE_NETCDF
TEST(WrfBoundaryMoisture, RuntimeRepackDoesNotAliasHydrometeors)
{
    const amrex::Box box(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    amrex::Vector<amrex::FArrayBox> data(WRFBdyVars::NumTypes);
    amrex::Arena* arena_used = amrex::The_Arena();
#ifdef AMREX_USE_GPU
    arena_used = amrex::The_Pinned_Arena();
#endif
    for (int ivar = 0; ivar < WRFBdyVars::NumTypes; ++ivar) {
        data[ivar].resize(box, 1, arena_used);
        data[ivar].template setVal<amrex::RunOn::Host>(static_cast<amrex::Real>(ivar));
    }

    repack_wrfbdy_to_realbdy(data, true);

    ASSERT_EQ(data.size(), static_cast<std::size_t>(RealBdyVars::NumTypes));
    EXPECT_DOUBLE_EQ(static_cast<double>(data[RealBdyVars::QC].const_array()(0, 0, 0)), 8.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(data[RealBdyVars::QI].const_array()(0, 0, 0)), 9.0);
}
#endif

TEST(WrfBoundaryMoisture, DefaultsPreserveLegacyBehavior)
{
    const SolverChoice choice{};
    EXPECT_FALSE(choice.use_wrf_bdy_qc_qi);
    EXPECT_EQ(choice.bdy_moist_nudge_type, 1);
}

TEST(WrfBoundaryMoisture, TargetInterpolatesThenWeightsDensityOnce)
{
    const amrex::Real rho = 1.25;
    const amrex::Real target = wrf_moisture_target(
        rho, amrex::Real(0.002), amrex::Real(0.006), amrex::Real(0.5));
    const double expected = 0.005;
    const double tol = 64.0 * std::numeric_limits<amrex::Real>::epsilon() *
                       std::max(1.0, std::abs(expected));
    EXPECT_NEAR(static_cast<double>(target), expected, tol);
}

TEST(WrfBoundaryMoisture, TargetSupportsErfOrWrfDensityPolicy)
{
    const amrex::Real q = 0.0078125;
    EXPECT_DOUBLE_EQ(static_cast<double>(wrf_moisture_target(1.0, q, q, 0.5)),
                     0.0078125);
    EXPECT_DOUBLE_EQ(static_cast<double>(wrf_moisture_target(1.5, q, q, 0.5)),
                     0.01171875);
}

TEST(WrfBoundaryMoisture, NegativeReconstructionClampsToZero)
{
    EXPECT_DOUBLE_EQ(static_cast<double>(wrf_moisture_target(
                         1.25, -0.002, -0.006, 0.5)), 0.0);
}

TEST(WrfBoundaryMoisture, ExplicitMappingDoesNotAliasRainToIce)
{
    const MoistureComponentIndices kessler(RhoQ1_comp, RhoQ2_comp, -1, RhoQ3_comp);
    EXPECT_EQ(wrf_bdy_var_for_moisture_component(RhoQ1_comp, kessler, true),
              RealBdyVars::QV);
    EXPECT_EQ(wrf_bdy_var_for_moisture_component(RhoQ2_comp, kessler, true),
              RealBdyVars::QC);
    EXPECT_EQ(wrf_bdy_var_for_moisture_component(RhoQ3_comp, kessler, true), -1);
}

TEST(WrfBoundaryMoisture, DisabledOptionLeavesHydrometeorsUnmapped)
{
    const MoistureComponentIndices ice_model(RhoQ1_comp, RhoQ2_comp, RhoQ3_comp);
    EXPECT_EQ(wrf_bdy_var_for_moisture_component(RhoQ2_comp, ice_model, false), -1);
    EXPECT_EQ(wrf_bdy_var_for_moisture_component(RhoQ3_comp, ice_model, false), -1);
}

TEST(WrfBoundaryMoisture, ModeThreeRelaxesMappedSpeciesOnEveryFaceOnly)
{
    const amrex::Box full_box(amrex::IntVect(0, 0, 0), amrex::IntVect(3, 3, 0));
    amrex::FArrayBox state(full_box, RhoQ4_comp + 1);
    amrex::FArrayBox rhs(full_box, RhoQ4_comp + 1);
    amrex::FArrayBox target_xlo(full_box, 3);
    amrex::FArrayBox target_xhi(full_box, 3);
    amrex::FArrayBox target_ylo(full_box, 3);
    amrex::FArrayBox target_yhi(full_box, 3);

    state.template setVal<amrex::RunOn::Device>(1.0);
    rhs.template setVal<amrex::RunOn::Device>(0.0);
    target_xlo.template setVal<amrex::RunOn::Device>(2.0);
    target_xhi.template setVal<amrex::RunOn::Device>(2.0);
    target_ylo.template setVal<amrex::RunOn::Device>(2.0);
    target_yhi.template setVal<amrex::RunOn::Device>(2.0);

    const amrex::Box bx_xlo(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    const amrex::Box bx_xhi(amrex::IntVect(3, 3, 0), amrex::IntVect(3, 3, 0));
    const amrex::Box bx_ylo(amrex::IntVect(2, 0, 0), amrex::IntVect(2, 0, 0));
    const amrex::Box bx_yhi(amrex::IntVect(1, 3, 0), amrex::IntVect(1, 3, 0));
    const amrex::GpuArray<int,3> comps{RhoQ1_comp, RhoQ2_comp, RhoQ3_comp};
    const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dx{1.0, 1.0, 1.0};
    const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> prob_lo{0.0, 0.0, 0.0};
    const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> prob_hi{4.0, 4.0, 1.0};

    realbdy_compute_moisture_relaxation(
        comps, 2, dx, prob_lo, prob_hi, 0.1,
        bx_xlo, bx_xhi, bx_ylo, bx_yhi,
        target_xlo.const_array(), target_xhi.const_array(),
        target_ylo.const_array(), target_yhi.const_array(),
        state.const_array(), rhs.array());
    amrex::Gpu::streamSynchronize();

    amrex::FArrayBox host_rhs(full_box, rhs.nComp(), amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     rhs.dataPtr(), rhs.dataPtr() + rhs.size(), host_rhs.dataPtr());
    const auto rhs_arr = host_rhs.const_array();
    const double rhs_tol = 64.0 * std::numeric_limits<amrex::Real>::epsilon();
    const amrex::IntVect samples[4] = {
        bx_xlo.smallEnd(), bx_xhi.smallEnd(), bx_ylo.smallEnd(), bx_yhi.smallEnd()};
    for (const auto& iv : samples) {
        for (const int comp : comps) {
            EXPECT_NEAR(rhs_arr(iv, comp), 0.05625, rhs_tol);
        }
        EXPECT_DOUBLE_EQ(rhs_arr(iv, RhoTheta_comp), 0.0);
        EXPECT_DOUBLE_EQ(rhs_arr(iv, RhoQ4_comp), 0.0);
    }
}

namespace {

// Relative tolerance for the pow-based number closures.
double number_tol (const double expected)
{
    return 1.0e3 * std::numeric_limits<amrex::Real>::epsilon() * std::abs(expected);
}

// rho*N = N0/lambda for an exponential distribution, with lambda taken from
// rho*q = pi rho_x N0 / lambda^4.  Written independently of compute_rn.
double exponential_number (const double n0, const double rho_x, const double rq)
{
    const double lambda = std::pow(M_PI * rho_x * n0 / rq, 0.25);
    return n0 / lambda;
}

}

TEST(WrfBoundaryMoisture, NumberTargetIsZeroBelowMassThreshold)
{
    const NumberDensityTarget ndt;
    for (int n = 0; n < NumberDensityTarget::n_n_max; ++n) {
        EXPECT_EQ(ndt.get_rn(n, amrex::Real(0.0), amrex::Real(260.0)), amrex::Real(0.0));
        EXPECT_EQ(ndt.get_rn(n, ndt.small, amrex::Real(260.0)), amrex::Real(0.0));
    }
}

TEST(WrfBoundaryMoisture, CloudAndIceNumbersUseDirectClosures)
{
    // Cloud water is a fixed droplet number, 250 cm^-3, whatever the mass;
    // cloud ice is Hong et al. (2004), N_i = 5.38e7 (rho q_i)^0.75.
    const NumberDensityTarget ndt;
    for (const double rq : {1.0e-6, 1.0e-4, 1.0e-3}) {
        EXPECT_DOUBLE_EQ(static_cast<double>(ndt.get_rn(NumberDensityTarget::cloud,
                                                        amrex::Real(rq), amrex::Real(280.0))),
                         2.5e8);
        const double expected = 5.38e7 * std::pow(rq, 0.75);
        EXPECT_NEAR(static_cast<double>(ndt.get_rn(NumberDensityTarget::ice,
                                                   amrex::Real(rq), amrex::Real(250.0))),
                    expected, number_tol(expected));
    }
}

TEST(WrfBoundaryMoisture, GraupelNumberMatchesExponentialClosure)
{
    const NumberDensityTarget ndt;
    const double rq = 1.0e-3;
    const double expected = exponential_number(4.0e6, 400.0, rq);
    EXPECT_NEAR(static_cast<double>(ndt.get_rn(NumberDensityTarget::graupel,
                                               amrex::Real(rq), amrex::Real(260.0))),
                expected, number_tol(expected));
}

TEST(WrfBoundaryMoisture, SnowInterceptFollowsTemperatureAndIsCapped)
{
    const NumberDensityTarget ndt;
    const double rq = 1.0e-4;

    // At or above freezing the intercept is the Houze et al. base value.
    const double warm = exponential_number(2.0e6, 100.0, rq);
    EXPECT_NEAR(static_cast<double>(ndt.get_rn(NumberDensityTarget::snow,
                                               amrex::Real(rq), amrex::Real(280.0))),
                warm, number_tol(warm));

    // Below freezing it grows as exp(0.12 (T0 - T)).
    const double t = 253.15;
    const double cold = exponential_number(2.0e6 * std::exp(0.12 * (273.15 - t)), 100.0, rq);
    EXPECT_NEAR(static_cast<double>(ndt.get_rn(NumberDensityTarget::snow,
                                               amrex::Real(rq), amrex::Real(t))),
                cold, number_tol(cold));

    // Very cold, the intercept is held at 1e11 m^-4.
    const double capped = exponential_number(1.0e11, 100.0, rq);
    EXPECT_NEAR(static_cast<double>(ndt.get_rn(NumberDensityTarget::snow,
                                               amrex::Real(rq), amrex::Real(150.0))),
                capped, number_tol(capped));
}

TEST(WrfBoundaryMoisture, RainInterceptFollowsAbelBoutleAboveMarshallPalmer)
{
    const NumberDensityTarget ndt;

    // Light rain: the Abel & Boutle intercept exceeds the Marshall-Palmer floor,
    // so the diagnosed number must satisfy N0 = 0.22 lambda^2.2 with lambda = N0/N.
    {
        const double rq = 1.0e-4;
        const double rn = ndt.get_rn(NumberDensityTarget::rain, amrex::Real(rq), amrex::Real(285.0));
        const double lambda = std::pow(M_PI * 997.0 * 0.22 / rq, 1.0 / (4.0 - 2.2));
        const double n0 = 0.22 * std::pow(lambda, 2.2);
        ASSERT_GT(n0, 8.0e6);
        const double expected = n0 / lambda;
        EXPECT_NEAR(rn, expected, number_tol(expected));
    }

    // Heavy rain: the Abel & Boutle intercept falls below 8e6 m^-4 and the floor applies.
    {
        const double rq = 1.0e-2;
        const double expected = exponential_number(8.0e6, 997.0, rq);
        EXPECT_NEAR(static_cast<double>(ndt.get_rn(NumberDensityTarget::rain,
                                                   amrex::Real(rq), amrex::Real(285.0))),
                    expected, number_tol(expected));
    }
}

TEST(WrfBoundaryMoisture, ExponentialNumberRecoversTheMassItCameFrom)
{
    // rho*N is per volume, the conserved number ERF transports: putting the
    // diagnosed rho*N back into the distribution must return the input rho*q.
    const NumberDensityTarget ndt;
    const double rq = 2.0e-3;
    const double n0 = 4.0e6;
    const double rn = ndt.compute_rn(amrex::Real(n0), amrex::Real(400.0), amrex::Real(rq));
    const double lambda = n0 / rn;
    const double mass = M_PI * 400.0 * n0 / std::pow(lambda, 4);
    EXPECT_NEAR(mass, rq, number_tol(rq));
}

TEST(WrfBoundaryMoisture, ModeThreeRelaxesMorrisonNumbersBehindTheMasses)
{
    // Mode 3 passes eleven targets for Morrison: the six masses, then the five
    // numbers.  Absent slots (-1) must be skipped.
    const MoistureComponentIndices mi =
        MoistureComponentIndices::from_moisture_model(MoistureType::Morrison);
    const amrex::GpuArray<int,11> comps{mi.qv, mi.qc, mi.qi, mi.qr, mi.qs, mi.qg,
                                        mi.nc, mi.ni, mi.nr, -1, mi.ng};
    const int ncomp_state = RhoQ11_comp + 1;

    const amrex::Box full_box(amrex::IntVect(0, 0, 0), amrex::IntVect(3, 3, 0));
    amrex::FArrayBox state(full_box, ncomp_state);
    amrex::FArrayBox rhs(full_box, ncomp_state);
    amrex::FArrayBox target_xlo(full_box, 11);
    amrex::FArrayBox target_xhi(full_box, 11);
    amrex::FArrayBox target_ylo(full_box, 11);
    amrex::FArrayBox target_yhi(full_box, 11);

    state.template setVal<amrex::RunOn::Device>(1.0);
    rhs.template setVal<amrex::RunOn::Device>(0.0);
    target_xlo.template setVal<amrex::RunOn::Device>(2.0);
    target_xhi.template setVal<amrex::RunOn::Device>(2.0);
    target_ylo.template setVal<amrex::RunOn::Device>(2.0);
    target_yhi.template setVal<amrex::RunOn::Device>(2.0);

    const amrex::Box bx_xlo(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    const amrex::Box bx_xhi(amrex::IntVect(3, 3, 0), amrex::IntVect(3, 3, 0));
    const amrex::Box bx_ylo(amrex::IntVect(2, 0, 0), amrex::IntVect(2, 0, 0));
    const amrex::Box bx_yhi(amrex::IntVect(1, 3, 0), amrex::IntVect(1, 3, 0));
    const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dx{1.0, 1.0, 1.0};
    const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> prob_lo{0.0, 0.0, 0.0};
    const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> prob_hi{4.0, 4.0, 1.0};

    realbdy_compute_moisture_relaxation(
        comps, 11, 2, dx, prob_lo, prob_hi, 0.1,
        bx_xlo, bx_xhi, bx_ylo, bx_yhi,
        target_xlo.const_array(), target_xhi.const_array(),
        target_ylo.const_array(), target_yhi.const_array(),
        state.const_array(), rhs.array());
    amrex::Gpu::streamSynchronize();

    amrex::FArrayBox host_rhs(full_box, rhs.nComp(), amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     rhs.dataPtr(), rhs.dataPtr() + rhs.size(), host_rhs.dataPtr());
    const auto rhs_arr = host_rhs.const_array();
    const double rhs_tol = 64.0 * std::numeric_limits<amrex::Real>::epsilon();
    const amrex::IntVect samples[4] = {
        bx_xlo.smallEnd(), bx_xhi.smallEnd(), bx_ylo.smallEnd(), bx_yhi.smallEnd()};
    for (const auto& iv : samples) {
        for (const int comp : comps) {
            if (comp >= 0) {
                EXPECT_NEAR(rhs_arr(iv, comp), 0.05625, rhs_tol);
            }
        }
        EXPECT_DOUBLE_EQ(rhs_arr(iv, mi.ns), 0.0);
        EXPECT_DOUBLE_EQ(rhs_arr(iv, RhoTheta_comp), 0.0);
    }
}

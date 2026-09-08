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
    for (int ivar = 0; ivar < WRFBdyVars::NumTypes; ++ivar) {
        data[ivar].resize(box, 1);
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

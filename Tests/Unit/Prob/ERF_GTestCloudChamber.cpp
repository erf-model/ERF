#include <AMReX_Array.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <cmath>

#include <gtest/gtest.h>

#include "../../../Source/Diffusion/ERF_ResolvedWallFlux.H"
#include "../../../Source/Prob/ERF_CloudChamber.H"

using amrex::GpuArray;
using amrex::Real;

namespace {

Real
single_value (const amrex::MultiFab& mf, const amrex::Box& box, int comp = 0)
{
    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<Real> reduce_data(reduce_op);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box overlap = box & mfi.validbox();
        if (overlap.isEmpty()) { continue; }
        const auto array = mf.const_array(mfi);
        reduce_op.eval(overlap, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<Real> {
                return { array(i,j,k,comp) };
            });
    }
    return amrex::get<0>(reduce_data.value());
}

int
sentinel_mismatches (const amrex::MultiFab& flux_mf, const amrex::Box& domain,
                     int dir, Real sentinel)
{
    const int lo = domain.smallEnd(dir);
    const int hi = domain.bigEnd(dir) + 1;
    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<int> reduce_data(reduce_op);
    for (amrex::MFIter mfi(flux_mf); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.validbox();
        const auto flux = flux_mf.const_array(mfi);
        reduce_op.eval(box, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<int> {
                const int coordinate = (dir == 0) ? i : ((dir == 1) ? j : k);
                const bool physical_face = coordinate == lo || coordinate == hi;
                const bool mismatch = physical_face ? (flux(i,j,k,0) == sentinel)
                                                    : (flux(i,j,k,0) != sentinel);
                return { mismatch ? 1 : 0 };
            });
    }
    return amrex::get<0>(reduce_data.value());
}

// Motivation: the chamber must use a decomposition-independent analytic
// initializer; these tests catch coordinate normalization and endpoint
// regressions without launching a full ERF simulation.
TEST(CloudChamberProfile, LinearProfileAndZeroAmplitude)
{
    erf_cloud_chamber::Config config;
    config.prob_lo = {Real(0.0), Real(0.0), Real(0.0)};
    config.prob_hi = {Real(2.0), Real(2.0), Real(1.0)};
    config.theta_bottom = Real(299.0);
    config.theta_top = Real(280.0);
    config.theta_perturbation_amplitude = Real(0.0);

    EXPECT_DOUBLE_EQ(erf_cloud_chamber::theta_at(config, Real(0.25),
                                                  Real(0.75), Real(0.0)),
                     Real(299.0));
    EXPECT_DOUBLE_EQ(erf_cloud_chamber::theta_at(config, Real(0.25),
                                                  Real(0.75), Real(1.0)),
                     Real(280.0));
    EXPECT_DOUBLE_EQ(erf_cloud_chamber::theta_at(config, Real(0.25),
                                                  Real(0.75), Real(0.5)),
                     Real(289.5));
}

TEST(CloudChamberProfile, PerturbationIsBoundedAndZeroOnVerticalFaces)
{
    erf_cloud_chamber::Config config;
    config.prob_lo = {Real(-1.0), Real(2.0), Real(4.0)};
    config.prob_hi = {Real(1.0), Real(6.0), Real(5.0)};
    config.theta_bottom = Real(300.0);
    config.theta_top = Real(301.0);
    config.theta_perturbation_amplitude = Real(0.2);
    const GpuArray<Real, AMREX_SPACEDIM> length = {Real(2.0), Real(4.0), Real(1.0)};

    const Real low = erf_cloud_chamber::deterministic_perturbation(
        Real(0.0), Real(4.0), Real(4.0), config.prob_lo, length, Real(0.2));
    const Real high = erf_cloud_chamber::deterministic_perturbation(
        Real(0.0), Real(4.0), Real(5.0), config.prob_lo, length, Real(0.2));
    const Real interior = erf_cloud_chamber::deterministic_perturbation(
        Real(0.0), Real(4.0), Real(4.5), config.prob_lo, length, Real(0.2));

    EXPECT_NEAR(low, Real(0.0), Real(1.0e-15));
    EXPECT_NEAR(high, Real(0.0), Real(1.0e-15));
    EXPECT_LE(std::abs(interior), Real(0.2));
}

// Motivation: physical chamber inputs are specified as temperature and RH;
// the initializer must convert RH to the exact mixing ratio used by the
// anelastic thermodynamics, rather than treating temperature as theta.
TEST(CloudChamberProfile, PhysicalTemperatureAndRelativeHumidityAreExact)
{
    erf_cloud_chamber::Config config;
    config.physical_initialization = true;
    config.prob_lo = {Real(0.0), Real(0.0), Real(0.0)};
    config.prob_hi = {Real(1.0), Real(1.0), Real(2.0)};
    config.initial_temperature_bottom = Real(300.0);
    config.initial_temperature_top = Real(284.0);
    config.temperature_perturbation_amplitude = Real(0.0);
    config.initial_relative_humidity = Real(0.95);

    const Real temperature = erf_cloud_chamber::theta_at(
        config, Real(0.25), Real(0.75), Real(0.5));
    const Real pressure = Real(100000.0);
    const Real qv = erf_cloud_chamber::vapor_mixing_ratio_from_relative_humidity(
        temperature, pressure, config.initial_relative_humidity);
    const Real vapor_pressure = pressure * qv / (RdoRv + qv);
    const Real recovered_rh = vapor_pressure /
        (Real(100.0) * erf_esatw(temperature));

    EXPECT_DOUBLE_EQ(temperature, Real(296.0));
    EXPECT_NEAR(recovered_rh, config.initial_relative_humidity, Real(1.0e-14));
    EXPECT_GT(qv, Real(0.0));
}

TEST(CloudChamberConfig, RejectsRelativeHumidityInDryPhysicalMode)
{
    const erf_cloud_chamber::InitializationContract contract {
        erf_cloud_chamber::InitializationMode::PhysicalTemperatureRH,
        false,
        true,
        false,
        false};

    const auto error = erf_cloud_chamber::initialization_contract_error(contract);

    EXPECT_EQ(error,
              "Cloud Chamber: prob.initial_relative_humidity is only used with erf.moisture_model = SatAdj");
}

TEST(CloudChamberConfig, RequiresRelativeHumidityForPhysicalSatAdj)
{
    const erf_cloud_chamber::InitializationContract contract {
        erf_cloud_chamber::InitializationMode::PhysicalTemperatureRH,
        true,
        false,
        false,
        false};

    const auto error = erf_cloud_chamber::initialization_contract_error(contract);

    EXPECT_EQ(error,
              "Cloud Chamber: physical SatAdj initialization requires prob.initial_relative_humidity");
}

TEST(CloudChamberConfig, RejectsLegacyKeysInPhysicalMode)
{
    const erf_cloud_chamber::InitializationContract contract {
        erf_cloud_chamber::InitializationMode::PhysicalTemperatureRH,
        false,
        false,
        true,
        false};

    const auto error = erf_cloud_chamber::initialization_contract_error(contract);

    EXPECT_EQ(error,
              "Cloud Chamber: physical_temperature_rh cannot be combined with legacy theta/qv profile keys");
}

TEST(CloudChamberConfig, RejectsPhysicalKeysInLegacyMode)
{
    const erf_cloud_chamber::InitializationContract contract {
        erf_cloud_chamber::InitializationMode::LegacyThetaQv,
        false,
        false,
        false,
        true};

    const auto error = erf_cloud_chamber::initialization_contract_error(contract);

    EXPECT_EQ(error,
              "Cloud Chamber: legacy_theta_qv cannot be combined with physical temperature/RH profile keys");
}

TEST(CloudChamberConfig, AcceptsPhysicalSatAdjWithRelativeHumidity)
{
    const erf_cloud_chamber::InitializationContract contract {
        erf_cloud_chamber::InitializationMode::PhysicalTemperatureRH,
        true,
        true,
        false,
        false};

    EXPECT_TRUE(erf_cloud_chamber::initialization_contract_error(contract).empty());
}

TEST(CloudChamberConfig, AcceptsLegacyThetaQvWithoutPhysicalKeys)
{
    const erf_cloud_chamber::InitializationContract contract {
        erf_cloud_chamber::InitializationMode::LegacyThetaQv,
        false,
        false,
        false,
        false};

    EXPECT_TRUE(erf_cloud_chamber::initialization_contract_error(contract).empty());
}

// Motivation: the physical wall override must replace the ordinary boundary
// flux with the signed half-cell molecular flux; dry faces must contribute
// exactly zero rather than imposing qv=0 as a Dirichlet state.
TEST(CloudChamberWallFlux, WetLowFaceAndDryHighFaceHaveExactSigns)
{
    using amrex::Box;
    using amrex::BoxArray;
    using amrex::DistributionMapping;
    using amrex::IntVect;
    using amrex::MultiFab;

    const Box domain(IntVect(0), IntVect(0));
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    MultiFab state(ba, dm, RhoQ2_comp + 1, 0);
    MultiFab prim(ba, dm, PrimQ2_comp + 1, 0);
    MultiFab base(ba, dm, BaseState::num_comps, 0);
    MultiFab rhs(ba, dm, RhoQ2_comp + 1, 0);
    BoxArray xba(ba); xba.surroundingNodes(0);
    BoxArray yba(ba); yba.surroundingNodes(1);
    BoxArray zba(ba); zba.surroundingNodes(2);
    MultiFab xflux(xba, dm, 1, 0);
    MultiFab yflux(yba, dm, 1, 0);
    MultiFab zflux(zba, dm, 1, 0);

    state.setVal(Real(0.0));
    prim.setVal(Real(0.0));
    base.setVal(Real(0.0));
    rhs.setVal(Real(0.0));
    xflux.setVal(Real(-9.0));
    yflux.setVal(Real(-9.0));
    zflux.setVal(Real(-9.0));
    state.setVal(Real(1.0), Rho_comp, 1);
    state.setVal(Real(0.01), RhoQ1_comp, 1);
    prim.setVal(Real(0.01), PrimQ1_comp, 1);
    base.setVal(Real(100000.0), BaseState::p0_comp, 1);

    erf_cloud_chamber::Config config;
    for (auto& wall : config.walls) {
        wall.thermodynamics.moisture_mode =
            erf_wall_thermodynamics::MoistureMode::DryImpermeable;
    }
    config.walls[0].thermodynamics.moisture_mode =
        erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
    config.walls[0].thermodynamics.physical_temperature_K = Real(300.0);

    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx_inv =
        {Real(1.0), Real(1.0), Real(1.0)};
    {
        amrex::MFIter mfi(state);
        ASSERT_TRUE(mfi.isValid());
        erf_resolved_wall_flux::apply(
            domain, domain, RhoQ1_comp, 0, state.const_array(mfi), prim.const_array(mfi),
            base.const_array(mfi), rhs.array(mfi), xflux.array(mfi), yflux.array(mfi),
            zflux.array(mfi), dx_inv, config.wall_boundary(), Real(0.0), Real(0.1), R_d/Cp_d);
    }
    amrex::Gpu::streamSynchronize();

    Real qsat = Real(0.0);
    erf_qsatw(Real(300.0), Real(1000.0), qsat);
    const Real expected_flux = Real(0.1) * (qsat - Real(0.01)) * Real(2.0);
    const amrex::Box xflux_box = xflux.boxArray()[0];
    amrex::Box low_face = xflux_box;
    low_face.setSmall(amrex::IntVect(0, 0, 0));
    low_face.setBig(amrex::IntVect(0, 0, 0));
    amrex::Box high_face = xflux_box;
    high_face.setSmall(amrex::IntVect(1, 0, 0));
    high_face.setBig(amrex::IntVect(1, 0, 0));
    EXPECT_NEAR(single_value(xflux, low_face), expected_flux, Real(1.0e-14));
    EXPECT_DOUBLE_EQ(single_value(xflux, high_face), Real(0.0));
    EXPECT_NEAR(single_value(rhs, domain, RhoQ1_comp), expected_flux, Real(1.0e-14));
}

// Motivation: resolved wall kernels receive one local tile at a time; an
// interior tile must not be relocated to a global boundary outside its FAB.
// This decomposition test exercises all six physical faces with multiple
// boxes and verifies that only physical boundary faces leave the sentinel.
TEST(CloudChamberWallFlux, MultiBoxOwnershipAcrossAllFaces)
{
    using amrex::Box;
    using amrex::BoxArray;
    using amrex::DistributionMapping;
    using amrex::IntVect;
    using amrex::MFIter;
    using amrex::MultiFab;

    const Box domain(IntVect(0), IntVect(3));
    BoxArray ba(domain);
    ba.maxSize(IntVect(2));
    const DistributionMapping dm(ba);
    MultiFab state(ba, dm, RhoQ2_comp + 1, 0);
    MultiFab prim(ba, dm, PrimQ2_comp + 1, 0);
    MultiFab base(ba, dm, BaseState::num_comps, 0);
    MultiFab rhs(ba, dm, RhoQ2_comp + 1, 0);
    BoxArray xba(ba); xba.surroundingNodes(0);
    BoxArray yba(ba); yba.surroundingNodes(1);
    BoxArray zba(ba); zba.surroundingNodes(2);
    MultiFab xflux(xba, dm, 1, 0);
    MultiFab yflux(yba, dm, 1, 0);
    MultiFab zflux(zba, dm, 1, 0);

    constexpr Real sentinel = Real(-9.0);
    state.setVal(Real(0.0));
    prim.setVal(Real(0.0));
    base.setVal(Real(0.0));
    rhs.setVal(Real(0.0));
    xflux.setVal(sentinel);
    yflux.setVal(sentinel);
    zflux.setVal(sentinel);
    state.setVal(Real(1.0), Rho_comp, 1);
    state.setVal(Real(0.01), RhoQ1_comp, 1);
    prim.setVal(Real(0.01), PrimQ1_comp, 1);
    base.setVal(Real(100000.0), BaseState::p0_comp, 1);

    erf_cloud_chamber::Config config;
    for (auto& wall : config.walls) {
        wall.thermodynamics.moisture_mode =
            erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
        wall.thermodynamics.physical_temperature_K = Real(300.0);
    }
    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx_inv =
        {Real(1.0), Real(1.0), Real(1.0)};

    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        erf_resolved_wall_flux::apply(
            mfi.validbox(), domain, RhoQ1_comp, 0,
            state.const_array(mfi), prim.const_array(mfi), base.const_array(mfi),
            rhs.array(mfi), xflux.array(mfi), yflux.array(mfi), zflux.array(mfi),
            dx_inv, config.wall_boundary(), Real(0.0), Real(0.1), R_d/Cp_d);
    }
    amrex::Gpu::streamSynchronize();

    EXPECT_EQ(sentinel_mismatches(xflux, domain, 0, sentinel), 0);
    EXPECT_EQ(sentinel_mismatches(yflux, domain, 1, sentinel), 0);
    EXPECT_EQ(sentinel_mismatches(zflux, domain, 2, sentinel), 0);
}

} // namespace

#include <AMReX_Array.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <cmath>
#include <limits>
#include <string>

#include <gtest/gtest.h>

#include "../../../Source/Diffusion/ERF_CloudChamberWallFlux.H"
#include "../../../Source/Diffusion/ERF_Diffusion.H"
#include "../../../Source/Diffusion/ERF_CloudChamberWallStress.H"
#include "../../../Source/Diffusion/ERF_ResolvedWallFlux.H"
#include "../../../Source/DataStructs/ERF_DataStruct.H"
#include "../../../Source/Prob/ERF_CloudChamber.H"
#include "../../../Source/Prob/ERF_ProblemDispatch.H"

using amrex::GpuArray;
using amrex::Real;

namespace {

Real
sum_region (const amrex::MultiFab& mf, const amrex::Box& box, int comp = 0)
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

// Return the first valid FAB value at a typed point. Nodal FABs can overlap
// at decomposition seams; selecting one valid owner avoids silently summing
// duplicated nodal storage.
amrex::Box
point_for (const amrex::MultiFab& mf, const amrex::IntVect& iv)
{
    return amrex::Box(iv, iv, mf.boxArray()[0].ixType());
}

Real
value_at (const amrex::MultiFab& mf, const amrex::IntVect& iv, int comp = 0)
{
    const amrex::Box point = point_for(mf, iv);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box overlap = point & mfi.validbox();
        if (overlap.isEmpty()) { continue; }
        amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
        amrex::ReduceData<Real> reduce_data(reduce_op);
        const auto array = mf.const_array(mfi);
        reduce_op.eval(overlap, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<Real> {
                return { array(i,j,k,comp) };
            });
        return amrex::get<0>(reduce_data.value());
    }
    return std::numeric_limits<Real>::quiet_NaN();
}

Real
scaled_tolerance (Real expected)
{
    const Real scale = std::max(Real(1.0), std::abs(expected));
    return Real(64.0) * std::numeric_limits<Real>::epsilon() * scale;
}

Real
staggered_tangential_speed (int dir, const amrex::Box& domain,
                            const amrex::MultiFab& u,
                            const amrex::MultiFab& v,
                            const amrex::MultiFab& w,
                            const erf_wall_thermodynamics::FaceWall& wall)
{
    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<Real> reduce_data(reduce_op);
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const auto ua = u.const_array(mfi);
        const auto va = v.const_array(mfi);
        const auto wa = w.const_array(mfi);
        reduce_op.eval(domain, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<Real> {
                return {erf_cloud_chamber_wall_flux::tangential_speed(
                    dir, i, j, k, ua, va, wa, wall)};
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


int
nonphysical_stress_changes (const amrex::MultiFab& stress_mf,
                            const amrex::Box& domain, int dir,
                            int allowed_coordinate, Real sentinel)
{
    amrex::ignore_unused(domain);
    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<int> reduce_data(reduce_op);
    for (amrex::MFIter mfi(stress_mf); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.validbox();
        const auto stress = stress_mf.const_array(mfi);
        reduce_op.eval(box, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<int> {
                const int coordinate = (dir == 0) ? i : ((dir == 1) ? j : k);
                return {(coordinate != allowed_coordinate &&
                         stress(i,j,k,0) != sentinel) ? 1 : 0};
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

// Motivation: the analytic perturbation is bounded and vanishes at the
// physical vertical endpoints; the tolerance must remain valid in both
// single- and double-precision builds.
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

    EXPECT_DOUBLE_EQ(low, Real(0.0));
    EXPECT_NEAR(high, Real(0.0), scaled_tolerance(Real(1.0)));
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
    EXPECT_NEAR(recovered_rh, config.initial_relative_humidity,
                scaled_tolerance(config.initial_relative_humidity));
    EXPECT_GT(qv, Real(0.0));
}

// Motivation: dry physical initialization must reject RH rather than
// silently accepting an input that changes the thermodynamic contract.
TEST(CloudChamberConfig, RejectsRelativeHumidityInDryPhysicalMode)
{
    const erf_cloud_chamber::InitializationContract contract {
        erf_cloud_chamber::InitializationMode::PhysicalTemperatureRH,
        false,
        true,
        false,
        false};

    const auto error = erf_cloud_chamber::initialization_contract_error(contract);

    EXPECT_NE(error.find("initial_relative_humidity"), std::string::npos);
    EXPECT_NE(error.find("SatAdj"), std::string::npos);
}

// Motivation: cloudy physical initialization needs an explicit RH so the
// SatAdj state is reproducible instead of relying on an implicit default.
TEST(CloudChamberConfig, RequiresRelativeHumidityForPhysicalSatAdj)
{
    const erf_cloud_chamber::InitializationContract contract {
        erf_cloud_chamber::InitializationMode::PhysicalTemperatureRH,
        true,
        false,
        false,
        false};

    const auto error = erf_cloud_chamber::initialization_contract_error(contract);

    EXPECT_NE(error.find("physical SatAdj"), std::string::npos);
    EXPECT_NE(error.find("initial_relative_humidity"), std::string::npos);
}

// Motivation: mixing physical and legacy profile keys would create an
// ambiguous initializer precedence that is difficult to audit.
TEST(CloudChamberConfig, RejectsLegacyKeysInPhysicalMode)
{
    const erf_cloud_chamber::InitializationContract contract {
        erf_cloud_chamber::InitializationMode::PhysicalTemperatureRH,
        false,
        false,
        true,
        false};

    const auto error = erf_cloud_chamber::initialization_contract_error(contract);

    EXPECT_NE(error.find("physical_temperature_rh"), std::string::npos);
    EXPECT_NE(error.find("legacy"), std::string::npos);
}

// Motivation: legacy numerical cases must retain their established input
// contract and reject physical temperature/RH aliases.
TEST(CloudChamberConfig, RejectsPhysicalKeysInLegacyMode)
{
    const erf_cloud_chamber::InitializationContract contract {
        erf_cloud_chamber::InitializationMode::LegacyThetaQv,
        false,
        false,
        false,
        true};

    const auto error = erf_cloud_chamber::initialization_contract_error(contract);

    EXPECT_NE(error.find("legacy_theta_qv"), std::string::npos);
    EXPECT_NE(error.find("physical"), std::string::npos);
}

// Motivation: the valid physical SatAdj contract must remain accepted while
// the stricter dry/cloudy key checks are enforced.
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

// Motivation: the generalized wall work must not regress the legacy
// theta/qv initialization path.
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

// Motivation: the production configuration seam, rather than a directly
// injected wall sample, must carry a non-default c_p into SolverChoice::rdOcp.
TEST(CloudChamberConfig, NonDefaultCpPropagatesToSolverChoiceRdOcp)
{
    constexpr char prefix[] = "cloud_chamber_cp_propagation_test";
    const Real configured_cp = Real(900.0);
    amrex::ParmParse pp(prefix);
    pp.add("init_type", "ConstantDensity");
    pp.add("c_p", configured_cp);

    SolverChoice solver_choice;
    solver_choice.init_params(0, prefix);

    EXPECT_NEAR(solver_choice.c_p, configured_cp,
                scaled_tolerance(configured_cp));
    EXPECT_NEAR(solver_choice.rdOcp, R_d / configured_cp,
                scaled_tolerance(R_d / configured_cp));
}

// Motivation: per-channel wall configuration must remain unambiguous; mixing
// the legacy aggregate key with a channel key must not create precedence.
TEST(CloudChamberWallConfig, RejectsAmbiguousAggregateAndChannelModels)
{
    erf_cloud_chamber::WallTransferContract contract;
    contract.legacy_aggregate_specified = true;
    contract.heat_model_specified = true;
    contract.heat_model = "bulk_aero";
    const auto error = erf_cloud_chamber::wall_transfer_contract_error(contract, "zlo");
    EXPECT_NE(error.find("wall_transfer_model"), std::string::npos);
    EXPECT_NE(error.find("per-channel"), std::string::npos);
}

// Motivation: a bulk wall must declare its coefficient source and required
// coefficient explicitly rather than falling back to a hidden default.
TEST(CloudChamberWallConfig, RequiresFixedCoefficientForBulkChannel)
{
    erf_cloud_chamber::WallTransferContract contract;
    contract.heat_model_specified = true;
    contract.heat_model = "bulk_aero";
    auto error = erf_cloud_chamber::wall_transfer_contract_error(contract, "zlo");
    EXPECT_NE(error.find("coefficient_source"), std::string::npos);

    contract.coefficient_source_specified = true;
    contract.coefficient_source = "fixed";
    error = erf_cloud_chamber::wall_transfer_contract_error(contract, "zlo");
    EXPECT_NE(error.find("C_H"), std::string::npos);
}

// Motivation: heat and vapor transfer are independent channels, so selecting
// bulk heat must not force bulk vapor or an unrelated vapor coefficient.
TEST(CloudChamberWallConfig, AcceptsIndependentHeatAndVaporBulkChannels)
{
    erf_cloud_chamber::WallTransferContract contract;
    contract.heat_model_specified = true;
    contract.heat_model = "bulk_aero";
    contract.vapor_model_specified = true;
    contract.vapor_model = "resolved_molecular";
    contract.coefficient_source_specified = true;
    contract.coefficient_source = "fixed";
    contract.heat_coefficient_specified = true;
    contract.heat_coefficient = Real(0.01);
    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(contract, "zlo").empty());
}

// Motivation: dry bulk vapor is algebraically impermeable and must not turn
// on the global wall-rate scan; wet bulk vapor and bulk heat remain active.
TEST(CloudChamberWallConfig, ActivatesTimestepGuardOnlyForActiveBulkChannels)
{
    erf_cloud_chamber::Config config;
    config.walls[0].wall.vapor.model =
        erf_wall_thermodynamics::ScalarModel::BulkAero;
    EXPECT_FALSE(config.has_bulk_scalar_wall());

    config.walls[0].wall.moisture =
        erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
    EXPECT_TRUE(config.has_bulk_scalar_wall());

    config.walls[0].wall.moisture =
        erf_wall_thermodynamics::MoistureMode::DryImpermeable;
    config.walls[0].wall.heat.model =
        erf_wall_thermodynamics::ScalarModel::BulkAero;
    EXPECT_TRUE(config.has_bulk_scalar_wall());
}

// The evaluator must be authoritative about legacy versus physical
// ownership, including owned-zero dry vapor and cloud-water channels.
TEST(CloudChamberWallFlux, EvaluatorOwnershipIsWallAware)
{
    using namespace erf_cloud_chamber_wall_flux;

    erf_wall_thermodynamics::FaceWall wall;
    ScalarWallSample sample;
    sample.rho = Real(1.0);
    sample.scalar_air = Real(0.0);
    sample.p_hse = Real(100000.0);
    sample.dx_inv = Real(1.0);
    sample.rdOcp = R_d / Cp_d;

    EXPECT_EQ(
        evaluate_scalar_flux_in(wall, ScalarChannel::Heat, sample).owned_channels,
        OwnNone);
    EXPECT_EQ(
        evaluate_scalar_flux_in(wall, ScalarChannel::Vapor, sample).owned_channels,
        OwnNone);
    EXPECT_EQ(
        evaluate_scalar_flux_in(wall, ScalarChannel::CloudWater, sample).owned_channels,
        OwnNone);

    wall.thermal.mode =
        erf_wall_thermodynamics::ThermalMode::FixedPhysicalTemperature;
    EXPECT_EQ(
        evaluate_scalar_flux_in(wall, ScalarChannel::Heat, sample).owned_channels,
        OwnHeat);

    wall.moisture =
        erf_wall_thermodynamics::MoistureMode::DryImpermeable;
    const auto dry =
        evaluate_scalar_flux_in(wall, ScalarChannel::Vapor, sample);
    EXPECT_EQ(dry.owned_channels, OwnVapor);
    EXPECT_EQ(dry.rhoQv_in, Real(0.0));

    const auto cloud =
        evaluate_scalar_flux_in(wall, ScalarChannel::CloudWater, sample);
    EXPECT_EQ(cloud.owned_channels, OwnCloudWater);
    EXPECT_EQ(cloud.rhoQc_in, Real(0.0));

    wall.vapor.model =
        erf_wall_thermodynamics::ScalarModel::BulkAero;
    EXPECT_FALSE(requires_tangential_speed(
        wall, ScalarChannel::Vapor));
    wall.moisture =
        erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
    EXPECT_TRUE(requires_tangential_speed(
        wall, ScalarChannel::Vapor));
}

// Motivation: dry vapor impermeability and cloud-water impermeability are
// algebraic gates that must run before NaN-prone coefficient or saturation
// arithmetic; the selected bulk equations must still be evaluated exactly.
// ScalarModel dispatch is intentionally explicit: adding a model must add a
// case here rather than silently falling through to resolved physics.

TEST(CloudChamberWallFlux, BulkFormulaeAndHardGates)
{
    using namespace erf_cloud_chamber_wall_flux;
    erf_wall_thermodynamics::FaceWall wall;
    wall.thermal.mode = erf_wall_thermodynamics::ThermalMode::FixedPhysicalTemperature;
    wall.thermal.temperature_K = Real(300.0);
    wall.heat.model = erf_wall_thermodynamics::ScalarModel::BulkAero;
    wall.heat.coefficient = Real(0.2);
    const Real rho = Real(1.1);
    const Real theta_air = Real(290.0);
    const Real p = Real(100000.0);
    const Real U = Real(3.0);
    const ScalarWallSample heat_sample{
        rho, theta_air, p, U, Real(0.0), Real(0.0), Real(1.0), R_d/Cp_d};
    const auto heat = evaluate_scalar_flux_in(
        wall, ScalarChannel::Heat, heat_sample);
    EXPECT_EQ(heat.owned_channels, OwnHeat);
    const Real theta_wall = Real(300.0) * std::pow(p_0/p, R_d/Cp_d);
    EXPECT_DOUBLE_EQ(heat.rhoTheta_in,
                     rho * Real(0.2) * U * (theta_wall - theta_air));

    wall.moisture = erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
    wall.vapor.model = erf_wall_thermodynamics::ScalarModel::BulkAero;
    wall.vapor.coefficient = Real(0.3);
    const Real qv_air = Real(0.01);
    const ScalarWallSample vapor_sample{
        rho, qv_air, p, U, Real(0.0), Real(0.0), Real(1.0), R_d/Cp_d};
    const auto vapor = evaluate_scalar_flux_in(
        wall, ScalarChannel::Vapor, vapor_sample);
    EXPECT_EQ(vapor.owned_channels, OwnVapor);
    Real qv_wall = Real(0.0);
    erf_qsatw(Real(300.0), Real(1000.0), qv_wall);
    EXPECT_DOUBLE_EQ(vapor.rhoQv_in,
                     rho * Real(0.3) * U * (qv_wall - qv_air));

    wall.moisture = erf_wall_thermodynamics::MoistureMode::DryImpermeable;
    wall.vapor.coefficient = std::numeric_limits<Real>::quiet_NaN();
    const ScalarWallSample dry_sample{
        rho, std::numeric_limits<Real>::quiet_NaN(), p,
        std::numeric_limits<Real>::quiet_NaN(),
        Real(0.0), Real(0.0), Real(1.0), R_d/Cp_d};
    const auto dry_vapor = evaluate_scalar_flux_in(
        wall, ScalarChannel::Vapor, dry_sample);
    EXPECT_EQ(dry_vapor.owned_channels, OwnVapor);
    const auto cloud_water = evaluate_scalar_flux_in(
        wall, ScalarChannel::CloudWater, dry_sample);
    EXPECT_EQ(cloud_water.owned_channels, OwnCloudWater);
    const ScalarWallSample unrelated_sample{
        rho, Real(0.0), p, Real(0.0),
        Real(0.0), Real(0.0), Real(1.0), R_d/Cp_d};
    const auto unrelated = evaluate_scalar_flux_in(
        wall, ScalarChannel::None, unrelated_sample);
    EXPECT_EQ(unrelated.owned_channels, OwnNone);
    EXPECT_DOUBLE_EQ(dry_vapor.rhoQv_in, Real(0.0));
    EXPECT_DOUBLE_EQ(cloud_water.rhoQc_in, Real(0.0));
}

// Motivation: bulk transfer must vanish exactly at calm conditions with no
// hidden velocity floor.
TEST(CloudChamberWallFlux, CalmBulkWallIsExactlyZero)
{
    using namespace erf_cloud_chamber_wall_flux;
    erf_wall_thermodynamics::FaceWall wall;
    wall.thermal.mode = erf_wall_thermodynamics::ThermalMode::FixedPhysicalTemperature;
    wall.thermal.temperature_K = Real(301.0);
    wall.heat.model = erf_wall_thermodynamics::ScalarModel::BulkAero;
    wall.heat.coefficient = Real(0.1);
    ScalarWallSample sample{
        Real(1.0), Real(280.0), Real(100000.0), Real(0.0),
        Real(0.0), Real(0.0), Real(1.0), R_d/Cp_d};
    const auto calm = evaluate_scalar_flux_in(
        wall, ScalarChannel::Heat, sample);
    EXPECT_DOUBLE_EQ(calm.rhoTheta_in, Real(0.0));
}

// Fixed bulk momentum is a production wall model, not just parser metadata.
// It must use the supplied C_D, act only on tangential velocity, and remain an
// exact zero at rest.
TEST(CloudChamberWallFlux, FixedBulkMomentumTractionAndCalmGate)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;

    FaceWall wall;
    wall.momentum.model = MomentumModel::BulkAero;
    wall.momentum.provider = CoefficientProvider::Fixed;
    wall.momentum.C_D = Real(0.2);

    MomentumWallSample sample;
    sample.rho = Real(1.25);
    sample.u_t = {Real(3.0), Real(4.0), Real(0.0)};
    sample.U_t = Real(5.0);
    const auto traction = evaluate_momentum_traction(wall, sample);
    EXPECT_EQ(traction.owned_channels, OwnMomentum);
    EXPECT_DOUBLE_EQ(traction.traction_on_fluid[0], -Real(1.25) * Real(0.2) * Real(5.0) * Real(3.0));
    EXPECT_DOUBLE_EQ(traction.traction_on_fluid[1], -Real(1.25) * Real(0.2) * Real(5.0) * Real(4.0));
    EXPECT_DOUBLE_EQ(traction.traction_on_fluid[2], Real(0.0));

    sample.U_t = Real(0.0);
    sample.u_t = {Real(0.0), Real(0.0), Real(0.0)};
    const auto calm = evaluate_momentum_traction(wall, sample);
    EXPECT_EQ(calm.owned_channels, OwnMomentum);
    EXPECT_DOUBLE_EQ(calm.traction_on_fluid[0], Real(0.0));
    EXPECT_DOUBLE_EQ(calm.traction_on_fluid[1], Real(0.0));
    EXPECT_DOUBLE_EQ(calm.traction_on_fluid[2], Real(0.0));
}

// MOST is shared by the production flux and timestep paths.  This test keeps
// the pointwise state contract visible: wet-wall thermodynamics, orientation,
// and distinct scalar roughness lengths all affect the result.
TEST(CloudChamberWallFlux, MOSTUsesWetStateOrientationAndSeparateRoughness)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;

    FaceWall wall;
    wall.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    wall.thermal.temperature_K = Real(300.0);
    wall.moisture = MoistureMode::WetEquilibrium;
    wall.momentum.model = MomentumModel::BulkAero;
    wall.momentum.provider = CoefficientProvider::MOST;
    wall.momentum.z0_m = Real(0.01);
    wall.heat.model = ScalarModel::BulkAero;
    wall.heat.provider = CoefficientProvider::MOST;
    wall.heat.z0 = Real(0.02);
    wall.vapor.model = ScalarModel::BulkAero;
    wall.vapor.provider = CoefficientProvider::MOST;
    wall.vapor.z0 = Real(0.04);

    const auto zlo = most_wall_coefficients(
        wall, Real(290.0), Real(0.005), Real(100000.0), R_d/Cp_d,
        Real(5.0), Real(0.5), 1);
    const auto zhi = most_wall_coefficients(
        wall, Real(290.0), Real(0.005), Real(100000.0), R_d/Cp_d,
        Real(5.0), Real(0.5), -1);
    EXPECT_EQ(zlo.valid, 1);
    EXPECT_EQ(zhi.valid, 1);
    EXPECT_GT(zlo.C_D, Real(0.0));
    EXPECT_GT(zlo.C_H, Real(0.0));
    EXPECT_GT(zlo.C_E, Real(0.0));
    EXPECT_NE(zlo.C_H, zlo.C_E);
    EXPECT_NE(zlo.zeta, zhi.zeta);

    // The conditioning floor is internal to MOST; the physical traction still
    // multiplies the actual tangential speed and is exactly zero at rest.
    MomentumWallSample calm;
    calm.rho = Real(1.0);
    calm.U_t = Real(0.0);
    calm.theta_air = Real(290.0);
    calm.qv_air = Real(0.005);
    calm.p_hse = Real(100000.0);
    calm.rdOcp = R_d/Cp_d;
    calm.wall_distance = Real(0.5);
    calm.gravity_sign = 1;
    const auto calm_traction = evaluate_momentum_traction(wall, calm);
    EXPECT_EQ(calm_traction.owned_channels, OwnMomentum);
    EXPECT_DOUBLE_EQ(calm_traction.traction_on_fluid[0], Real(0.0));
    EXPECT_DOUBLE_EQ(calm_traction.traction_on_fluid[1], Real(0.0));
    EXPECT_DOUBLE_EQ(calm_traction.traction_on_fluid[2], Real(0.0));
}

TEST(CloudChamberWallFlux, MOSTUsesConfiguredRdOcpForStabilityAndHeatFlux)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;

    FaceWall wall;
    wall.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    wall.thermal.temperature_K = Real(300.0);
    wall.moisture = MoistureMode::WetEquilibrium;
    wall.momentum.model = MomentumModel::BulkAero;
    wall.momentum.provider = CoefficientProvider::MOST;
    wall.momentum.z0_m = Real(0.01);
    wall.heat.model = ScalarModel::BulkAero;
    wall.heat.provider = CoefficientProvider::MOST;
    wall.heat.z0 = Real(0.02);
    wall.vapor.model = ScalarModel::BulkAero;
    wall.vapor.provider = CoefficientProvider::MOST;
    wall.vapor.z0 = Real(0.02);

    const Real theta_air = Real(285.0);
    const Real qv_air = Real(0.005);
    const Real p_hse = Real(85000.0);
    const Real U_t = Real(5.0);
    const Real wall_distance = Real(0.25);
    const Real configured_rdOcp = R_d / Real(900.0);
    const Real expected_theta_wall = wall.thermal.temperature_K *
        std::pow(p_0 / p_hse, configured_rdOcp);
    const auto default_runtime = most_wall_coefficients(
        wall, theta_air, qv_air, p_hse, R_d/Cp_d, U_t, wall_distance, 1);
    const auto configured_runtime = most_wall_coefficients(
        wall, theta_air, qv_air, p_hse, configured_rdOcp, U_t,
        wall_distance, 1);
    EXPECT_EQ(configured_runtime.valid, 1);
    EXPECT_NE(configured_runtime.zeta, default_runtime.zeta);
    EXPECT_NE(configured_runtime.C_D, default_runtime.C_D);
    EXPECT_NE(configured_runtime.C_H, default_runtime.C_H);

    ScalarWallSample sample;
    sample.rho = Real(1.2);
    sample.scalar_air = theta_air;
    sample.p_hse = p_hse;
    sample.U_t = U_t;
    sample.rdOcp = configured_rdOcp;
    sample.wall_distance = wall_distance;
    sample.theta_air = theta_air;
    sample.qv_air = qv_air;
    sample.gravity_sign = 1;
    const auto heat = evaluate_scalar_flux_in(wall, ScalarChannel::Heat, sample);
    const Real expected_heat = sample.rho * configured_runtime.C_H * U_t *
        (expected_theta_wall - theta_air);
    EXPECT_NEAR(heat.rhoTheta_in, expected_heat, scaled_tolerance(expected_heat));
}

TEST(CloudChamberNeutralLog, LegacyStressOverloadDoesNotRequireBaseState)
{
    using namespace erf_wall_thermodynamics;
    const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(0));
    const amrex::BoxArray ba(domain);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab state(ba, dm, Rho_comp + 1, 1);
    amrex::BoxArray xba(ba); xba.surroundingNodes(0);
    amrex::BoxArray yba(ba); yba.surroundingNodes(1);
    amrex::BoxArray zba(ba); zba.surroundingNodes(2);
    amrex::MultiFab u(xba, dm, 1, 1);
    amrex::MultiFab v(yba, dm, 1, 1);
    amrex::MultiFab w(zba, dm, 1, 1);
    amrex::BoxArray ba12(ba); ba12.surroundingNodes(0); ba12.surroundingNodes(1);
    amrex::BoxArray ba13(ba); ba13.surroundingNodes(0); ba13.surroundingNodes(2);
    amrex::BoxArray ba23(ba); ba23.surroundingNodes(1); ba23.surroundingNodes(2);
    amrex::MultiFab tau12(ba12, dm, 1, 1);
    amrex::MultiFab tau13(ba13, dm, 1, 1);
    amrex::MultiFab tau23(ba23, dm, 1, 1);
    state.setVal(Real(1.2), Rho_comp, 1);
    u.setVal(Real(1.0));
    v.setVal(Real(2.0));
    w.setVal(Real(3.0));
    tau12.setVal(Real(0.0));
    tau13.setVal(Real(0.0));
    tau23.setVal(Real(0.0));

    Boundary walls{};
    walls[0].momentum.model = MomentumModel::NeutralRoughnessLog;
    walls[0].momentum.z0_m = Real(0.01);
    const GpuArray<Real, AMREX_SPACEDIM> dx_inv = {
        Real(2.0), Real(2.0), Real(2.0)};
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        erf_cloud_chamber_wall_stress::apply(
            mfi.validbox(), domain, state.const_array(mfi), u.const_array(mfi),
            v.const_array(mfi), w.const_array(mfi), tau12.array(mfi),
            tau13.array(mfi), tau23.array(mfi), dx_inv, walls);
    }
    amrex::Gpu::streamSynchronize();

    const Real U_t = std::sqrt(Real(13.0));
    const Real cd = std::pow(KAPPA/std::log(Real(25.0)), Real(2.0));
    EXPECT_NEAR(value_at(tau12, amrex::IntVect(0,0,0)),
                -Real(1.2) * cd * U_t * Real(2.0),
                scaled_tolerance(cd));
}

// Motivation: bulk transfer must depend only on velocity tangent to the
// wall. This uses the production staggered-array helper to catch normal
// leakage and coordinate-rotation/indexing errors in all three directions.
TEST(CloudChamberWallFlux, StaggeredTangentialSpeedRotatesAndRemovesNormal)
{
    using amrex::Box;
    using amrex::BoxArray;
    using amrex::DistributionMapping;
    using amrex::IntVect;
    using amrex::MultiFab;

    const Box domain(IntVect(0), IntVect(0));
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    BoxArray xba(ba); xba.surroundingNodes(0);
    BoxArray yba(ba); yba.surroundingNodes(1);
    BoxArray zba(ba); zba.surroundingNodes(2);
    MultiFab u(xba, dm, 1, 0);
    MultiFab v(yba, dm, 1, 0);
    MultiFab w(zba, dm, 1, 0);
    erf_wall_thermodynamics::FaceWall low_wall;
    erf_wall_thermodynamics::FaceWall high_wall;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (dir == 0) {
            u.setVal(Real(100.0));
            v.setVal(Real(3.0));
            w.setVal(Real(4.0));
        } else if (dir == 1) {
            u.setVal(Real(3.0));
            v.setVal(Real(100.0));
            w.setVal(Real(4.0));
        } else {
            u.setVal(Real(3.0));
            v.setVal(Real(4.0));
            w.setVal(Real(100.0));
        }
        const Real low_speed = staggered_tangential_speed(
            dir, domain, u, v, w, low_wall);
        const Real high_speed = staggered_tangential_speed(
            dir, domain, u, v, w, high_wall);
        EXPECT_NEAR(low_speed, Real(5.0), scaled_tolerance(Real(5.0)));
        EXPECT_NEAR(high_speed, Real(5.0), scaled_tolerance(Real(5.0)));
    }
}

// Motivation: conservation alone cannot prove the selected bulk model ran;
// a no-op or stale resolved flux can still close a domain budget. This test
// verifies retained low/high face values and RHS corrections at the actual
// generalized wall-application seam, including coefficient linearity.
TEST(CloudChamberWallFlux, GeneralizedApplyActivatesBulkFluxAndRhsCorrection)
{
    using amrex::Box;
    using amrex::BoxArray;
    using amrex::DistributionMapping;
    using amrex::IntVect;
    using amrex::MultiFab;
    using namespace erf_cloud_chamber_wall_flux;

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
    MultiFab u(xba, dm, 1, 0);
    MultiFab v(yba, dm, 1, 0);
    MultiFab w(zba, dm, 1, 0);

    const Real rho = Real(1.2);
    const Real theta_air = Real(290.0);
    const Real pressure = Real(100000.0);
    const Real old_flux = Real(7.0);
    const Real dx_inv = Real(2.0);
    const Real U_t = Real(5.0);
    const Real coefficient = Real(0.1);
    const Real theta_wall = Real(300.0) * std::pow(p_0 / pressure, R_d / Cp_d);
    const Real expected_bulk =
        rho * coefficient * U_t * (theta_wall - theta_air);

    state.setVal(Real(0.0));
    prim.setVal(Real(0.0));
    base.setVal(Real(0.0));
    rhs.setVal(Real(0.0));
    xflux.setVal(old_flux);
    yflux.setVal(old_flux);
    zflux.setVal(old_flux);
    state.setVal(rho, Rho_comp, 1);
    prim.setVal(theta_air, PrimTheta_comp, 1);
    base.setVal(pressure, BaseState::p0_comp, 1);
    u.setVal(Real(100.0));
    v.setVal(Real(3.0));
    w.setVal(Real(4.0));

    erf_wall_thermodynamics::Boundary walls{};
    for (int face : {0, 1}) {
        walls[face].thermal.mode =
            erf_wall_thermodynamics::ThermalMode::FixedPhysicalTemperature;
        walls[face].thermal.temperature_K = Real(300.0);
        walls[face].heat.model =
            erf_wall_thermodynamics::ScalarModel::BulkAero;
        walls[face].heat.coefficient = coefficient;
    }
    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx =
        {dx_inv, dx_inv, dx_inv};

    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        erf_cloud_chamber_wall_flux::apply(
            mfi.validbox(), domain, RhoTheta_comp, 0,
            state.const_array(mfi), prim.const_array(mfi), base.const_array(mfi),
            u.const_array(mfi), v.const_array(mfi), w.const_array(mfi),
            rhs.array(mfi), xflux.array(mfi), yflux.array(mfi), zflux.array(mfi),
            dx, walls, Real(0.0), Real(0.0), R_d/Cp_d);
    }
    amrex::Gpu::streamSynchronize();

    Box low_face = xflux.boxArray()[0];
    low_face.setSmall(0, 0);
    low_face.setBig(0, 0);
    Box high_face = xflux.boxArray()[0];
    high_face.setSmall(0, 1);
    high_face.setBig(0, 1);
    EXPECT_NEAR(sum_region(xflux, low_face), expected_bulk,
                scaled_tolerance(expected_bulk));
    EXPECT_NEAR(sum_region(xflux, high_face), -expected_bulk,
                scaled_tolerance(expected_bulk));
    const Real expected_rhs =
        ((expected_bulk - old_flux) - (-expected_bulk - old_flux)) * dx_inv;
    EXPECT_NEAR(sum_region(rhs, domain, RhoTheta_comp), expected_rhs,
                scaled_tolerance(expected_rhs));

    auto resolved_wall = walls[0];
    resolved_wall.heat.model =
        erf_wall_thermodynamics::ScalarModel::ResolvedMolecular;
    const ScalarWallSample resolved_sample{
        rho, theta_air, pressure, U_t,
        Real(0.05), Real(0.0), dx_inv, R_d/Cp_d};
    const auto resolved = evaluate_scalar_flux_in(
        resolved_wall, ScalarChannel::Heat, resolved_sample);
    EXPECT_GT(std::abs(expected_bulk - resolved.rhoTheta_in),
              scaled_tolerance(expected_bulk));

    walls[0].heat.coefficient = Real(0.2);
    walls[1].heat.coefficient = Real(0.2);
    rhs.setVal(Real(0.0));
    xflux.setVal(old_flux);
    const Real doubled_bulk = Real(2.0) * expected_bulk;
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        erf_cloud_chamber_wall_flux::apply(
            mfi.validbox(), domain, RhoTheta_comp, 0,
            state.const_array(mfi), prim.const_array(mfi), base.const_array(mfi),
            u.const_array(mfi), v.const_array(mfi), w.const_array(mfi),
            rhs.array(mfi), xflux.array(mfi), yflux.array(mfi), zflux.array(mfi),
            dx, walls, Real(0.0), Real(0.0), R_d/Cp_d);
    }
    amrex::Gpu::streamSynchronize();
    EXPECT_NEAR(sum_region(xflux, low_face), doubled_bulk,
                scaled_tolerance(doubled_bulk));
    EXPECT_NEAR(sum_region(xflux, high_face), -doubled_bulk,
                scaled_tolerance(doubled_bulk));
}

// Motivation: the wet vapor path must prove that selected bulk_aero dispatch
// reaches the real face-replacement seam. A resolved fallback can be nonzero
// and conservative, so a generic wet-budget activation check is insufficient.
// Exercise one face at a time so the seeded old-flux term cannot cancel from
// the RHS oracle, then require the C_E factor-of-two response.
TEST(CloudChamberWallFlux,
     GeneralizedApplyActivatesWetBulkVaporFluxAndRhsCorrection)
{
    using amrex::Box;
    using amrex::BoxArray;
    using amrex::DistributionMapping;
    using amrex::IntVect;
    using amrex::MultiFab;
    using namespace erf_cloud_chamber_wall_flux;

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
    MultiFab u(xba, dm, 1, 0);
    MultiFab v(yba, dm, 1, 0);
    MultiFab w(zba, dm, 1, 0);

    const Real rho = Real(1.2);
    const Real qv_air = Real(0.01);
    const Real pressure = Real(100000.0);
    const Real wall_temperature = Real(300.0);
    const Real old_flux = Real(7.0);
    const Real dx_inv = Real(2.0);
    const Real U_t = Real(5.0);
    const Real coefficient = Real(0.1);
    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx =
        {dx_inv, dx_inv, dx_inv};

    Real qv_wall = Real(0.0);
    erf_qsatw(wall_temperature, pressure * Real(0.01), qv_wall);
    ASSERT_GT(qv_wall, qv_air);

    state.setVal(Real(0.0));
    prim.setVal(Real(0.0));
    base.setVal(Real(0.0));
    state.setVal(rho, Rho_comp, 1);
    state.setVal(rho * qv_air, RhoQ1_comp, 1);
    prim.setVal(qv_air, PrimQ1_comp, 1);
    base.setVal(pressure, BaseState::p0_comp, 1);

    // x is deliberately a large normal velocity. Only v/w should contribute
    // to the x-wall tangential speed: sqrt(3^2 + 4^2) = 5.
    u.setVal(Real(100.0));
    v.setVal(Real(3.0));
    w.setVal(Real(4.0));

    auto face_box = [&](int face_index) {
        Box box = xflux.boxArray()[0];
        box.setSmall(0, face_index);
        box.setBig(0, face_index);
        return box;
    };

    auto run_face = [&](bool high, Real C_E) {
        rhs.setVal(Real(0.0));
        xflux.setVal(old_flux);
        yflux.setVal(old_flux);
        zflux.setVal(old_flux);

        erf_wall_thermodynamics::Boundary walls{};
        const int face = high ? 1 : 0;
        walls[face].thermal.mode =
            erf_wall_thermodynamics::ThermalMode::FixedPhysicalTemperature;
        walls[face].thermal.temperature_K = wall_temperature;
        walls[face].moisture =
            erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
        walls[face].vapor.model =
            erf_wall_thermodynamics::ScalarModel::BulkAero;
        walls[face].vapor.provider =
            erf_wall_thermodynamics::CoefficientProvider::Fixed;
        walls[face].vapor.coefficient = C_E;

        for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
            erf_cloud_chamber_wall_flux::apply(
                mfi.validbox(), domain, RhoQ1_comp, 0,
                state.const_array(mfi), prim.const_array(mfi),
                base.const_array(mfi),
                u.const_array(mfi), v.const_array(mfi), w.const_array(mfi),
                rhs.array(mfi), xflux.array(mfi), yflux.array(mfi),
                zflux.array(mfi), dx, walls,
                Real(0.0), Real(0.0), R_d/Cp_d);
        }
        amrex::Gpu::streamSynchronize();

        const Real expected_inward =
            rho * C_E * U_t * (qv_wall - qv_air);
        const Real expected_coordinate = high ? -expected_inward
                                              : expected_inward;
        const Real expected_rhs = high
            ? -(expected_coordinate - old_flux) * dx_inv
            :  (expected_coordinate - old_flux) * dx_inv;

        const int active_index = high ? 1 : 0;
        const int inactive_index = high ? 0 : 1;
        const Real retained = sum_region(
            xflux, face_box(active_index));

        EXPECT_NEAR(retained, expected_coordinate,
                    scaled_tolerance(expected_coordinate));
        EXPECT_NEAR(sum_region(rhs, domain, RhoQ1_comp), expected_rhs,
                    scaled_tolerance(expected_rhs));
        EXPECT_DOUBLE_EQ(sum_region(xflux, face_box(inactive_index)),
                         old_flux);
        return retained;
    };

    const Real low_CE = run_face(false, coefficient);
    const Real high_CE = run_face(true, coefficient);
    EXPECT_GT(low_CE, Real(0.0));
    EXPECT_LT(high_CE, Real(0.0));

    const Real low_2CE = run_face(false, Real(2.0) * coefficient);
    EXPECT_NEAR(low_2CE, Real(2.0) * low_CE,
                scaled_tolerance(Real(2.0) * low_CE));
}

// Motivation: the single low/high sign adapter and dt_wall=0.5/max_rate helper
// are shared seams; all Cartesian orientations must use the same convention.
TEST(CloudChamberWallFlux, OrientationAndTimestepAdaptersCoverAllFaces)
{
    using namespace erf_cloud_chamber_wall_flux;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (dir == 0) {
   EXPECT_DOUBLE_EQ((to_coordinate_flux<0,false>(Real(2.0))), Real(2.0));
   EXPECT_DOUBLE_EQ((to_coordinate_flux<0,true>(Real(2.0))), Real(-2.0));
            const auto lo = traction_from_local_inward<0,false>(Real(2.0), Real(3.0), Real(4.0));
            const auto hi = traction_from_local_inward<0,true>(Real(2.0), Real(3.0), Real(4.0));
            EXPECT_DOUBLE_EQ(lo[0], Real(2.0));
            EXPECT_DOUBLE_EQ(lo[1], Real(3.0));
            EXPECT_DOUBLE_EQ(lo[2], Real(4.0));
            EXPECT_DOUBLE_EQ(hi[0], Real(-2.0));
            EXPECT_DOUBLE_EQ(hi[1], Real(3.0));
            EXPECT_DOUBLE_EQ(hi[2], Real(4.0));
        } else if (dir == 1) {
   EXPECT_DOUBLE_EQ((to_coordinate_flux<1,false>(Real(2.0))), Real(2.0));
   EXPECT_DOUBLE_EQ((to_coordinate_flux<1,true>(Real(2.0))), Real(-2.0));
            const auto hi = traction_from_local_inward<1,true>(Real(2.0), Real(3.0), Real(4.0));
            EXPECT_DOUBLE_EQ(hi[0], Real(3.0));
            EXPECT_DOUBLE_EQ(hi[1], Real(-2.0));
            EXPECT_DOUBLE_EQ(hi[2], Real(4.0));
        } else {
   EXPECT_DOUBLE_EQ((to_coordinate_flux<2,false>(Real(2.0))), Real(2.0));
   EXPECT_DOUBLE_EQ((to_coordinate_flux<2,true>(Real(2.0))), Real(-2.0));
            const auto hi = traction_from_local_inward<2,true>(Real(2.0), Real(3.0), Real(4.0));
            EXPECT_DOUBLE_EQ(hi[0], Real(3.0));
            EXPECT_DOUBLE_EQ(hi[1], Real(4.0));
            EXPECT_DOUBLE_EQ(hi[2], Real(-2.0));
        }
    }
    EXPECT_DOUBLE_EQ(wall_dt_from_max_rate(Real(0.0)), Real(1.0e30));
    EXPECT_DOUBLE_EQ(wall_dt_from_max_rate(Real(2.0)), Real(0.25));
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
        wall.wall.moisture =
            erf_wall_thermodynamics::MoistureMode::DryImpermeable;
    }
    config.walls[0].wall.moisture =
        erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
    // Dry bulk metadata must not require velocity in the resolved compatibility API.
    config.walls[1].wall.vapor.model =
        erf_wall_thermodynamics::ScalarModel::BulkAero;

    config.walls[0].wall.thermal.mode =
        erf_wall_thermodynamics::ThermalMode::FixedPhysicalTemperature;
    config.walls[0].wall.thermal.temperature_K = Real(300.0);

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
    EXPECT_NEAR(sum_region(xflux, low_face), expected_flux,
                scaled_tolerance(expected_flux));
    EXPECT_DOUBLE_EQ(sum_region(xflux, high_face), Real(0.0));
    EXPECT_NEAR(sum_region(rhs, domain, RhoQ1_comp), expected_flux,
                scaled_tolerance(expected_flux));
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
        wall.wall.moisture =
            erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
        wall.wall.thermal.mode =
            erf_wall_thermodynamics::ThermalMode::FixedPhysicalTemperature;
        wall.wall.thermal.temperature_K = Real(300.0);
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

// Motivation: the resolved compatibility API supplies empty velocity views,
// so its classification must follow active closure requirements rather than
// scalar-model metadata alone.
TEST(CloudChamberWallFlux, CompatibilityVelocityRequirementIsChannelAware)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;

    Boundary walls{};
    walls[0].thermal.mode = ThermalMode::FixedPhysicalTemperature;
    walls[0].moisture = MoistureMode::WetEquilibrium;
    walls[0].heat.model = ScalarModel::ResolvedMolecular;
    walls[0].vapor.model = ScalarModel::ResolvedMolecular;

    // A: resolved scalar wall.
    EXPECT_FALSE(has_velocity_dependent_scalar_wall(walls));

    // B: physical bulk heat.
    walls[0].heat.model = ScalarModel::BulkAero;
    EXPECT_TRUE(has_velocity_dependent_scalar_wall(walls));

    // C: wet bulk vapor.
    walls[0].heat.model = ScalarModel::ResolvedMolecular;
    walls[0].vapor.model = ScalarModel::BulkAero;
    EXPECT_TRUE(has_velocity_dependent_scalar_wall(walls));

    // D: dry bulk vapor metadata is an owned exact-zero gate.
    walls[0].moisture = MoistureMode::DryImpermeable;
    EXPECT_FALSE(has_velocity_dependent_scalar_wall(walls));
}


TEST(CloudChamberNeutralLog, ComposedMomentumRowSumIncludesPerpendicularWalls)
{
    const GpuArray<Real, AMREX_SPACEDIM> low = {Real(1.0), Real(2.0), Real(3.0)};
    const GpuArray<Real, AMREX_SPACEDIM> high = {Real(4.0), Real(5.0), Real(6.0)};
    EXPECT_DOUBLE_EQ(
        erf_cloud_chamber_wall_flux::momentum_row_sum_rate(0, low, high),
        Real(16.0));
    EXPECT_DOUBLE_EQ(
        erf_cloud_chamber_wall_flux::momentum_row_sum_rate(1, low, high),
        Real(14.0));
    EXPECT_DOUBLE_EQ(
        erf_cloud_chamber_wall_flux::momentum_row_sum_rate(2, low, high),
        Real(12.0));
}


TEST(CloudChamberNeutralLog, UserDefinedVelocityDispatchIsGeneric)
{
    using erf_problem_dispatch::CustomVelocityInitializer;
    using erf_problem_dispatch::custom_velocity_initializer;
    EXPECT_EQ(custom_velocity_initializer("cloud chamber"),
              CustomVelocityInitializer::CloudChamber);
    EXPECT_EQ(custom_velocity_initializer("cloudchamber"),
              CustomVelocityInitializer::CloudChamber);
    EXPECT_EQ(custom_velocity_initializer("userdefined"),
              CustomVelocityInitializer::UserDefined);
    EXPECT_NE(custom_velocity_initializer("userdefined"),
              CustomVelocityInitializer::CloudChamber);
}


TEST(CloudChamberNeutralLog, MomentumAndScalarClosures)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;

    FaceWall wall;
    wall.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    wall.thermal.temperature_K = Real(300.0);
    wall.moisture = MoistureMode::WetEquilibrium;
    wall.momentum.model = MomentumModel::NeutralRoughnessLog;
    wall.momentum.z0_m = Real(0.01);
    wall.heat.model = ScalarModel::NeutralRoughnessLog;
    wall.heat.z0 = Real(0.02);
    wall.vapor.model = ScalarModel::NeutralRoughnessLog;
    wall.vapor.z0 = Real(0.03);

    MomentumWallSample momentum;
    momentum.rho = Real(1.2);
    momentum.u_t = {Real(3.0), Real(4.0), Real(0.0)};
    momentum.U_t = Real(5.0);
    momentum.wall_distance = Real(0.5);
    const Real log_m = std::log(Real(50.0));
    const Real C_D = std::pow(KAPPA / log_m, Real(2.0));
    const auto neutral = neutral_log_momentum_state(
        momentum.U_t, momentum.wall_distance, wall.momentum.z0_m);
    EXPECT_NEAR(neutral.log_m, log_m, scaled_tolerance(log_m));
    EXPECT_NEAR(neutral.C_D, C_D, scaled_tolerance(C_D));
    const auto traction = evaluate_momentum_traction(wall, momentum);
    EXPECT_EQ(traction.owned_channels, OwnMomentum);
    EXPECT_NEAR(traction.traction_on_fluid[0], -Real(1.2)*C_D*Real(5.0)*Real(3.0),
                scaled_tolerance(C_D));
    EXPECT_NEAR(traction.traction_on_fluid[1], -Real(1.2)*C_D*Real(5.0)*Real(4.0),
                scaled_tolerance(C_D));
    EXPECT_DOUBLE_EQ(traction.traction_on_fluid[2], Real(0.0));

    ScalarWallSample scalar;
    scalar.rho = Real(1.1);
    scalar.scalar_air = Real(290.0);
    scalar.p_hse = Real(100000.0);
    scalar.U_t = Real(4.0);
    scalar.rdOcp = R_d / Cp_d;
    scalar.wall_distance = Real(0.5);
    const Real expected_C_H = KAPPA * KAPPA /
        (std::log(scalar.wall_distance / wall.momentum.z0_m) *
         std::log(scalar.wall_distance / wall.heat.z0));
    const Real expected_theta_wall = wall.thermal.temperature_K *
        std::pow(p_0 / scalar.p_hse, scalar.rdOcp);
    const Real expected_heat = scalar.rho * expected_C_H * scalar.U_t *
        (expected_theta_wall - scalar.scalar_air);
    const auto heat = evaluate_scalar_flux_in(wall, ScalarChannel::Heat, scalar);
    EXPECT_EQ(heat.owned_channels, OwnHeat);
    EXPECT_NEAR(heat.rhoTheta_in, expected_heat,
                scaled_tolerance(expected_heat));

    scalar.scalar_air = Real(0.01);
    const Real expected_C_E = KAPPA * KAPPA /
        (std::log(scalar.wall_distance / wall.momentum.z0_m) *
         std::log(scalar.wall_distance / wall.vapor.z0));
    Real qv_wall = Real(0.0);
    erf_qsatw(wall.thermal.temperature_K, scalar.p_hse * Real(0.01), qv_wall);
    const Real expected_vapor = scalar.rho * expected_C_E * scalar.U_t *
        (qv_wall - scalar.scalar_air);
    const auto vapor = evaluate_scalar_flux_in(wall, ScalarChannel::Vapor, scalar);
    EXPECT_EQ(vapor.owned_channels, OwnVapor);
    EXPECT_NEAR(vapor.rhoQv_in, expected_vapor,
                scaled_tolerance(expected_vapor));

    wall.moisture = MoistureMode::DryImpermeable;
    scalar.U_t = std::numeric_limits<Real>::quiet_NaN();
    scalar.p_hse = std::numeric_limits<Real>::quiet_NaN();
    const auto dry = evaluate_scalar_flux_in(wall, ScalarChannel::Vapor, scalar);
    EXPECT_EQ(dry.owned_channels, OwnVapor);
    EXPECT_DOUBLE_EQ(dry.rhoQv_in, Real(0.0));
    EXPECT_FALSE(requires_tangential_speed(wall, ScalarChannel::Vapor));
}

TEST(CloudChamberNeutralLog, ContractGeometryAndRateGates)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;
    using erf_cloud_chamber::WallTransferContract;

    WallTransferContract contract;
    contract.momentum_model_specified = true;
    contract.momentum_model = "neutral_roughness_log";
    contract.heat_model_specified = true;
    contract.heat_model = "neutral_roughness_log";
    contract.vapor_model_specified = true;
    contract.vapor_model = "neutral_roughness_log";
    contract.z0_m_specified = true;
    contract.z0_h_specified = true;
    contract.z0_q_specified = true;
    contract.z0_m = Real(0.01);
    contract.z0_h = Real(0.02);
    contract.z0_q = Real(0.03);
    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(contract, "xlo").empty());

    contract.z0_h_specified = false;
    EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(contract, "xlo").find("z0_h"),
              std::string::npos);
    contract.z0_h_specified = true;
    contract.z0_m = Real(0.0);
    EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(contract, "xlo").find("finite and positive"),
              std::string::npos);

    Boundary walls{};
    walls[0].momentum.model = MomentumModel::NeutralRoughnessLog;
    walls[0].momentum.z0_m = Real(0.49);
    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx = {Real(1.0), Real(2.0), Real(3.0)};
    EXPECT_TRUE(erf_cloud_chamber::wall_roughness_geometry_error(walls, dx).empty());
    walls[0].momentum.z0_m = Real(0.5);
    const auto geometry_error = erf_cloud_chamber::wall_roughness_geometry_error(walls, dx);
    EXPECT_NE(geometry_error.find("xlo"), std::string::npos);
    EXPECT_NE(geometry_error.find("z0_m"), std::string::npos);
    EXPECT_NE(geometry_error.find("z_ref"), std::string::npos);

    FaceWall rate_wall;
    rate_wall.momentum.model = MomentumModel::NeutralRoughnessLog;
    rate_wall.momentum.z0_m = Real(0.01);
    EXPECT_TRUE(wall_rate_requires_tangential_speed(rate_wall));
    const Real U_t = Real(3.0);
    const Real dx_inv = Real(2.0);
    const Real cd = neutral_log_momentum_state(U_t, Real(0.5)/dx_inv,
                                                rate_wall.momentum.z0_m).C_D;
    EXPECT_NEAR(wall_rate_for_face(rate_wall, U_t, dx_inv),
                momentum_infinity_row_sum_factor() * cd * U_t * dx_inv,
                scaled_tolerance(cd));
    EXPECT_DOUBLE_EQ(neutral_scalar_rate(Real(0.3), Real(0.0), dx_inv), Real(0.0));
    FaceWall dry_vapor;
    dry_vapor.moisture = MoistureMode::DryImpermeable;
    dry_vapor.vapor.model = ScalarModel::NeutralRoughnessLog;
    dry_vapor.vapor.z0 = Real(0.02);
    EXPECT_FALSE(wall_rate_requires_tangential_speed(dry_vapor));
    EXPECT_DOUBLE_EQ(wall_rate_for_face(
        dry_vapor, std::numeric_limits<Real>::quiet_NaN(), dx_inv), Real(0.0));
}


TEST(CloudChamberNeutralLog, HostParserContractMatrix)
{
    using erf_cloud_chamber::WallTransferContract;

    const auto momentum_only = [] {
        WallTransferContract c;
        c.momentum_model_specified = true;
        c.momentum_model = "neutral_roughness_log";
        c.z0_m_specified = true;
        c.z0_m = Real(0.01);
        return c;
    };
    const auto heat_only = [] {
        WallTransferContract c;
        c.heat_model_specified = true;
        c.heat_model = "neutral_roughness_log";
        c.z0_m_specified = true;
        c.z0_m = Real(0.01);
        c.z0_h_specified = true;
        c.z0_h = Real(0.005);
        return c;
    };
    const auto vapor_only = [] {
        WallTransferContract c;
        c.vapor_model_specified = true;
        c.vapor_model = "neutral_roughness_log";
        c.z0_m_specified = true;
        c.z0_m = Real(0.01);
        c.z0_q_specified = true;
        c.z0_q = Real(0.005);
        return c;
    };

    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(
        momentum_only(), "xlo").empty());
    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(
        heat_only(), "ylo").empty());
    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(
        vapor_only(), "zlo").empty());

    auto mixed = momentum_only();
    mixed.heat_model_specified = true;
    mixed.heat_model = "bulk_aero";
    mixed.coefficient_source_specified = true;
    mixed.coefficient_source = "fixed";
    mixed.heat_coefficient_specified = true;
    mixed.heat_coefficient = Real(0.01);
    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(
        mixed, "xhi").empty());

    auto all_neutral = momentum_only();
    all_neutral.heat_model_specified = true;
    all_neutral.heat_model = "neutral_roughness_log";
    all_neutral.z0_h_specified = true;
    all_neutral.z0_h = Real(0.005);
    all_neutral.vapor_model_specified = true;
    all_neutral.vapor_model = "neutral_roughness_log";
    all_neutral.z0_q_specified = true;
    all_neutral.z0_q = Real(0.005);
    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(
        all_neutral, "zhi").empty());

    auto missing_momentum_z0 = momentum_only();
    missing_momentum_z0.z0_m_specified = false;
    EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(
        missing_momentum_z0, "xlo").find("z0_m"), std::string::npos);

    auto missing_heat_z0 = heat_only();
    missing_heat_z0.z0_h_specified = false;
    EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(
        missing_heat_z0, "ylo").find("z0_h"), std::string::npos);

    auto missing_vapor_z0 = vapor_only();
    missing_vapor_z0.z0_q_specified = false;
    EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(
        missing_vapor_z0, "zlo").find("z0_q"), std::string::npos);

    auto nonfinite_z0 = momentum_only();
    nonfinite_z0.z0_m = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(
        nonfinite_z0, "xlo").find("finite and positive"), std::string::npos);

    auto zero_z0 = momentum_only();
    zero_z0.z0_m = Real(0.0);
    EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(
        zero_z0, "xlo").find("finite and positive"), std::string::npos);

    auto negative_z0 = momentum_only();
    negative_z0.z0_m = Real(-0.01);
    EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(
        negative_z0, "xlo").find("finite and positive"), std::string::npos);

    auto unused_z0 = momentum_only();
    unused_z0.z0_h_specified = true;
    unused_z0.z0_h = Real(0.005);
    EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(
        unused_z0, "xlo").find("z0_h"), std::string::npos);

    auto neutral_ch = heat_only();
    neutral_ch.heat_coefficient_specified = true;
    neutral_ch.heat_coefficient = Real(0.01);
    EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(
        neutral_ch, "xlo").find("C_H"), std::string::npos);

    auto neutral_ce = vapor_only();
    neutral_ce.vapor_coefficient_specified = true;
    neutral_ce.vapor_coefficient = Real(0.01);
    EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(
        neutral_ce, "xlo").find("C_E"), std::string::npos);

    auto neutral_provider = momentum_only();
    neutral_provider.coefficient_source_specified = true;
    neutral_provider.coefficient_source = "fixed";
    EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(
        neutral_provider, "xlo").find("coefficient_source"), std::string::npos);

    auto most = momentum_only();
    most.momentum_model = "bulk_aero";
    most.coefficient_source_specified = true;
    most.coefficient_source = "most";
    EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(
        most, "xlo").find("coefficient_source = most"), std::string::npos);
    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(
        most, "zlo").empty());

    auto smooth = momentum_only();
    smooth.momentum_model = "law_of_wall_momentum";
    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(
        smooth, "xlo").empty());

    // Strict parsing still requires z0_q for a dry wall; the runtime dry
    // vapor gate is tested independently by MomentumAndScalarClosures.
    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(
        vapor_only(), "xlo").empty());
}

TEST(CloudChamberNeutralLog, GeometryValidationUsesEveryActiveFace)
{
    using namespace erf_wall_thermodynamics;

    Boundary walls{};
    walls[0].momentum.model = MomentumModel::NeutralRoughnessLog;
    walls[0].momentum.z0_m = Real(0.49);
    walls[2].heat.model = ScalarModel::NeutralRoughnessLog;
    walls[2].momentum.z0_m = Real(0.99);
    walls[2].heat.z0 = Real(0.99);
    walls[4].vapor.model = ScalarModel::NeutralRoughnessLog;
    walls[4].momentum.z0_m = Real(1.49);
    walls[4].vapor.z0 = Real(1.49);
    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx = {
        Real(1.0), Real(2.0), Real(3.0)};
    EXPECT_TRUE(erf_cloud_chamber::wall_roughness_geometry_error(
        walls, dx).empty());

    walls[0].momentum.z0_m = Real(0.5);
    auto error = erf_cloud_chamber::wall_roughness_geometry_error(walls, dx);
    EXPECT_NE(error.find("xlo"), std::string::npos);
    EXPECT_NE(error.find("z0_m"), std::string::npos);
    walls[0].momentum.z0_m = Real(0.49);

    walls[2].heat.z0 = Real(1.0);
    error = erf_cloud_chamber::wall_roughness_geometry_error(walls, dx);
    EXPECT_NE(error.find("ylo"), std::string::npos);
    EXPECT_NE(error.find("z0_h"), std::string::npos);
    walls[2].heat.z0 = Real(0.99);

    walls[4].vapor.z0 = Real(1.5);
    error = erf_cloud_chamber::wall_roughness_geometry_error(walls, dx);
    EXPECT_NE(error.find("zlo"), std::string::npos);
    EXPECT_NE(error.find("z0_q"), std::string::npos);

    Boundary scalar_only{};
    scalar_only[0].momentum.model = MomentumModel::ResolvedNoSlip;
    scalar_only[0].heat.model = ScalarModel::NeutralRoughnessLog;
    scalar_only[0].heat.z0 = Real(0.25);
    scalar_only[0].momentum.z0_m = Real(0.49);
    EXPECT_TRUE(erf_cloud_chamber::wall_roughness_geometry_error(
        scalar_only, dx).empty());
    scalar_only[0].momentum.z0_m = Real(0.0);
    error = erf_cloud_chamber::wall_roughness_geometry_error(scalar_only, dx);
    EXPECT_NE(error.find("z0_m"), std::string::npos);
    scalar_only[0].momentum.z0_m = Real(0.5);
    error = erf_cloud_chamber::wall_roughness_geometry_error(scalar_only, dx);
    EXPECT_NE(error.find("xlo"), std::string::npos);
    EXPECT_NE(error.find("z0_m"), std::string::npos);
    scalar_only[0].momentum.z0_m = Real(0.51);
    error = erf_cloud_chamber::wall_roughness_geometry_error(scalar_only, dx);
    EXPECT_NE(error.find("z0_m"), std::string::npos);

    scalar_only = Boundary{};
    scalar_only[0].momentum.model = MomentumModel::ResolvedNoSlip;
    scalar_only[0].vapor.model = ScalarModel::NeutralRoughnessLog;
    scalar_only[0].vapor.z0 = Real(0.25);
    scalar_only[0].momentum.z0_m = Real(0.49);
    EXPECT_TRUE(erf_cloud_chamber::wall_roughness_geometry_error(
        scalar_only, dx).empty());
    scalar_only[0].momentum.z0_m = Real(0.0);
    error = erf_cloud_chamber::wall_roughness_geometry_error(scalar_only, dx);
    EXPECT_NE(error.find("z0_m"), std::string::npos);
    scalar_only[0].momentum.z0_m = Real(0.5);
    error = erf_cloud_chamber::wall_roughness_geometry_error(scalar_only, dx);
    EXPECT_NE(error.find("z0_m"), std::string::npos);
    scalar_only[0].momentum.z0_m = Real(0.51);
    error = erf_cloud_chamber::wall_roughness_geometry_error(scalar_only, dx);
    EXPECT_NE(error.find("z0_m"), std::string::npos);
}

TEST(CloudChamberNeutralLog, AdapterMapsAllTwelveFaceComponents)
{
    using namespace erf_wall_thermodynamics;
    const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(2));
    const amrex::BoxArray ba(domain);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab state(ba, dm, Rho_comp + 1, 1);
    amrex::BoxArray xba(ba); xba.surroundingNodes(0);
    amrex::BoxArray yba(ba); yba.surroundingNodes(1);
    amrex::BoxArray zba(ba); zba.surroundingNodes(2);
    amrex::MultiFab u(xba, dm, 1, 1);
    amrex::MultiFab v(yba, dm, 1, 1);
    amrex::MultiFab w(zba, dm, 1, 1);
    amrex::BoxArray ba12(ba); ba12.surroundingNodes(0); ba12.surroundingNodes(1);
    amrex::BoxArray ba13(ba); ba13.surroundingNodes(0); ba13.surroundingNodes(2);
    amrex::BoxArray ba23(ba); ba23.surroundingNodes(1); ba23.surroundingNodes(2);
    amrex::MultiFab tau12(ba12, dm, 1, 1);
    amrex::MultiFab tau13(ba13, dm, 1, 1);
    amrex::MultiFab tau23(ba23, dm, 1, 1);

    state.setVal(Real(1.2), Rho_comp, 1);
    u.setVal(Real(1.0));
    v.setVal(Real(2.0));
    w.setVal(Real(3.0));
    constexpr Real sentinel = Real(-17.0);
    tau12.setVal(sentinel); tau13.setVal(sentinel); tau23.setVal(sentinel);

    Boundary walls{};
    for (auto& wall : walls) {
        wall.momentum.model = MomentumModel::NeutralRoughnessLog;
        wall.momentum.z0_m = Real(0.01);
    }
    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx_inv = {Real(2.0), Real(2.0), Real(2.0)};
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        erf_cloud_chamber_wall_stress::apply(
            mfi.validbox(), domain, state.const_array(mfi), u.const_array(mfi),
            v.const_array(mfi), w.const_array(mfi), tau12.array(mfi),
            tau13.array(mfi), tau23.array(mfi), dx_inv, walls);
    }
    amrex::Gpu::streamSynchronize();

    const Real cd = std::pow(KAPPA/std::log(Real(25.0)), Real(2.0));
    const Real U_x = std::sqrt(Real(13.0));
    const Real U_y = std::sqrt(Real(10.0));
    const Real U_z = std::sqrt(Real(5.0));
    const Real tau_x_y = -Real(1.2)*cd*U_x*Real(2.0);
    const Real tau_x_z = -Real(1.2)*cd*U_x*Real(3.0);
    const Real tau_y_x = -Real(1.2)*cd*U_y*Real(1.0);
    const Real tau_y_z = -Real(1.2)*cd*U_y*Real(3.0);
    const Real tau_z_x = -Real(1.2)*cd*U_z*Real(1.0);
    const Real tau_z_y = -Real(1.2)*cd*U_z*Real(2.0);
    const auto check = [](const amrex::MultiFab& mf, const amrex::IntVect& iv,
                          Real expected) {
        EXPECT_NEAR(value_at(mf, iv), expected,
                    scaled_tolerance(expected));
    };

    check(tau12, amrex::IntVect(0,1,1), tau_x_y);
    check(tau12, amrex::IntVect(3,1,1), -tau_x_y);
    check(tau13, amrex::IntVect(0,1,1), tau_x_z);
    check(tau13, amrex::IntVect(3,1,1), -tau_x_z);
    check(tau12, amrex::IntVect(1,0,1), tau_y_x);
    check(tau12, amrex::IntVect(1,3,1), -tau_y_x);
    check(tau23, amrex::IntVect(1,0,1), tau_y_z);
    check(tau23, amrex::IntVect(1,3,1), -tau_y_z);
    check(tau13, amrex::IntVect(1,1,0), tau_z_x);
    check(tau13, amrex::IntVect(1,1,3), -tau_z_x);
    check(tau23, amrex::IntVect(1,1,0), tau_z_y);
    check(tau23, amrex::IntVect(1,1,3), -tau_z_y);

    check(tau12, amrex::IntVect(1,1,1), sentinel);
    check(tau13, amrex::IntVect(1,1,1), sentinel);
    check(tau23, amrex::IntVect(1,1,1), sentinel);
}


TEST(CloudChamberNeutralLog, AdapterRespectsMultiBoxPhysicalOwnership)
{
    using namespace erf_wall_thermodynamics;
    const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(3));
    amrex::BoxArray ba(domain);
    ba.maxSize(2);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab state(ba, dm, Rho_comp + 1, 1);
    amrex::BoxArray xba(ba); xba.surroundingNodes(0);
    amrex::BoxArray yba(ba); yba.surroundingNodes(1);
    amrex::BoxArray zba(ba); zba.surroundingNodes(2);
    amrex::MultiFab u(xba, dm, 1, 1);
    amrex::MultiFab v(yba, dm, 1, 1);
    amrex::MultiFab w(zba, dm, 1, 1);
    amrex::BoxArray ba12(ba); ba12.surroundingNodes(0); ba12.surroundingNodes(1);
    amrex::BoxArray ba13(ba); ba13.surroundingNodes(0); ba13.surroundingNodes(2);
    amrex::BoxArray ba23(ba); ba23.surroundingNodes(1); ba23.surroundingNodes(2);
    amrex::MultiFab tau12(ba12, dm, 1, 1);
    amrex::MultiFab tau13(ba13, dm, 1, 1);
    amrex::MultiFab tau23(ba23, dm, 1, 1);

    state.setVal(Real(1.2), Rho_comp, 1);
    u.setVal(Real(1.0));
    v.setVal(Real(2.0));
    w.setVal(Real(3.0));
    constexpr Real sentinel = Real(-23.0);
    tau12.setVal(sentinel);
    tau13.setVal(sentinel);
    tau23.setVal(sentinel);

    Boundary walls{};
    walls[0].momentum.model = MomentumModel::NeutralRoughnessLog;
    walls[0].momentum.z0_m = Real(0.01);
    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx_inv = {
        Real(2.0), Real(2.0), Real(2.0)};
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        erf_cloud_chamber_wall_stress::apply(
            mfi.validbox(), domain, state.const_array(mfi), u.const_array(mfi),
            v.const_array(mfi), w.const_array(mfi), tau12.array(mfi),
            tau13.array(mfi), tau23.array(mfi), dx_inv, walls);
    }
    amrex::Gpu::streamSynchronize();

    const amrex::IntVect iv(0,1,1);
    EXPECT_NE(value_at(tau12, iv), sentinel);
    EXPECT_NE(value_at(tau13, iv), sentinel);
    EXPECT_EQ(nonphysical_stress_changes(tau12, domain, 0, 0, sentinel), 0);
    EXPECT_EQ(nonphysical_stress_changes(tau13, domain, 0, 0, sentinel), 0);
    EXPECT_EQ(nonphysical_stress_changes(tau23, domain, 0, 0, sentinel), 0);
}

TEST(CloudChamberNeutralLog, SingleWallImpulseUsesDiffusionOperator)
{
    const auto run_face = [](bool high) {
        const amrex::Box cell_box(amrex::IntVect(0), amrex::IntVect(1));
        const amrex::BoxArray ba(cell_box);
        const amrex::DistributionMapping dm(ba);
        amrex::MultiFab rho_u_rhs(amrex::BoxArray(
            amrex::surroundingNodes(cell_box, 0)), dm, 1, 0);
        amrex::MultiFab rho_v_rhs(amrex::BoxArray(
            amrex::surroundingNodes(cell_box, 1)), dm, 1, 0);
        amrex::MultiFab rho_w_rhs(amrex::BoxArray(
            amrex::surroundingNodes(cell_box, 2)), dm, 1, 0);
        rho_u_rhs.setVal(Real(0.0));
        rho_v_rhs.setVal(Real(0.0));
        rho_w_rhs.setVal(Real(0.0));

        amrex::MultiFab tau11(ba, dm, 1, 1);
        amrex::MultiFab tau22(ba, dm, 1, 1);
        amrex::MultiFab tau33(ba, dm, 1, 1);
        amrex::MultiFab tau12(amrex::BoxArray(
            amrex::convert(cell_box, amrex::IntVect(1,1,0))), dm, 1, 1);
        amrex::MultiFab tau21(amrex::BoxArray(
            amrex::convert(cell_box, amrex::IntVect(1,1,0))), dm, 1, 1);
        amrex::MultiFab tau13(amrex::BoxArray(
            amrex::convert(cell_box, amrex::IntVect(1,0,1))), dm, 1, 1);
        amrex::MultiFab tau31(amrex::BoxArray(
            amrex::convert(cell_box, amrex::IntVect(1,0,1))), dm, 1, 1);
        amrex::MultiFab tau23(amrex::BoxArray(
            amrex::convert(cell_box, amrex::IntVect(0,1,1))), dm, 1, 1);
        amrex::MultiFab tau32(amrex::BoxArray(
            amrex::convert(cell_box, amrex::IntVect(0,1,1))), dm, 1, 1);
        tau11.setVal(Real(0.0));
        tau22.setVal(Real(0.0));
        tau33.setVal(Real(0.0));
        tau12.setVal(Real(0.0));
        tau21.setVal(Real(0.0));
        tau13.setVal(Real(0.0));
        tau31.setVal(Real(0.0));
        tau23.setVal(Real(0.0));
        tau32.setVal(Real(0.0));

        constexpr Real physical_traction = Real(-1.7);
        const int i_face = high ? 2 : 0;
        const amrex::IntVect sample_iv(i_face, 1, 0);
        const amrex::Box sample(
            sample_iv, sample_iv, tau12.boxArray()[0].ixType());
        const Real stored_traction = high ? -physical_traction : physical_traction;
        tau12.setVal(stored_traction, sample, 0, 1, 0);

        amrex::MultiFab detJ(ba, dm, 1, 1);
        amrex::MultiFab mf_mx(ba, dm, 1, 1);
        amrex::MultiFab mf_ux(ba, dm, 1, 1);
        amrex::MultiFab mf_vx(ba, dm, 1, 1);
        amrex::MultiFab mf_my(ba, dm, 1, 1);
        amrex::MultiFab mf_uy(ba, dm, 1, 1);
        amrex::MultiFab mf_vy(ba, dm, 1, 1);
        detJ.setVal(Real(1.0));
        mf_mx.setVal(Real(1.0));
        mf_ux.setVal(Real(1.0));
        mf_vx.setVal(Real(1.0));
        mf_my.setVal(Real(1.0));
        mf_uy.setVal(Real(1.0));
        mf_vy.setVal(Real(1.0));
        amrex::Gpu::DeviceVector<Real> stretched_dz_d;
        amrex::GpuArray<Real, AMREX_SPACEDIM> dx_inv = {
            Real(1.0), Real(1.0), Real(1.0)};

        DiffusionSrcForMom(
            amrex::surroundingNodes(cell_box, 0),
            amrex::surroundingNodes(cell_box, 1),
            amrex::surroundingNodes(cell_box, 2),
            rho_u_rhs[0].array(), rho_v_rhs[0].array(), rho_w_rhs[0].array(),
            tau11[0].const_array(), tau22[0].const_array(),
            tau33[0].const_array(), tau12[0].const_array(),
            tau21[0].const_array(), tau13[0].const_array(),
            tau31[0].const_array(), tau23[0].const_array(),
            tau32[0].const_array(), detJ[0].const_array(), stretched_dz_d,
            dx_inv, mf_mx[0].const_array(), mf_ux[0].const_array(),
            mf_vx[0].const_array(), mf_my[0].const_array(),
            mf_uy[0].const_array(), mf_vy[0].const_array(), false, false);
        amrex::Gpu::streamSynchronize();

        const int i_rhs = high ? 1 : 0;
        const amrex::IntVect point_iv(i_rhs, 1, 0);
        return value_at(rho_v_rhs, point_iv);
    };

    constexpr Real expected = Real(-1.7);
    EXPECT_NEAR(run_face(false), expected, scaled_tolerance(expected));
    EXPECT_NEAR(run_face(true), expected, scaled_tolerance(expected));
}

TEST(CloudChamberNeutralLog, ComposedMomentumInfinityNormBound)
{
    using namespace erf_wall_thermodynamics;
    const amrex::Box cell_box(amrex::IntVect(0), amrex::IntVect(1));
    const amrex::BoxArray ba(cell_box);
    const amrex::DistributionMapping dm(ba);
    constexpr Real u_base = Real(1.3);
    constexpr Real v_base = Real(0.7);
    constexpr Real w_base = Real(0.4);

    const auto evaluate = [&](int perturbed_component,
                              const amrex::IntVect& perturbation,
                              Real delta) {
        amrex::MultiFab state(ba, dm, Rho_comp + 1, 1);
        amrex::BoxArray xba(ba); xba.surroundingNodes(0);
        amrex::BoxArray yba(ba); yba.surroundingNodes(1);
        amrex::BoxArray zba(ba); zba.surroundingNodes(2);
        amrex::MultiFab u(xba, dm, 1, 1);
        amrex::MultiFab v(yba, dm, 1, 1);
        amrex::MultiFab w(zba, dm, 1, 1);
        state.setVal(Real(1.0), Rho_comp, 1);
        u.setVal(u_base);
        v.setVal(v_base);
        w.setVal(w_base);
        if (perturbed_component == 0) {
            const amrex::Box point(
                perturbation, perturbation, u.boxArray()[0].ixType());
            u.setVal(u_base + delta, point, 0, 1, 0);
        } else if (perturbed_component == 1) {
            const amrex::Box point(
                perturbation, perturbation, v.boxArray()[0].ixType());
            v.setVal(v_base + delta, point, 0, 1, 0);
        } else if (perturbed_component == 2) {
            const amrex::Box point(
                perturbation, perturbation, w.boxArray()[0].ixType());
            w.setVal(w_base + delta, point, 0, 1, 0);
        }

        amrex::BoxArray ba12(ba); ba12.surroundingNodes(0); ba12.surroundingNodes(1);
        amrex::BoxArray ba13(ba); ba13.surroundingNodes(0); ba13.surroundingNodes(2);
        amrex::BoxArray ba23(ba); ba23.surroundingNodes(1); ba23.surroundingNodes(2);
        amrex::MultiFab tau12(ba12, dm, 1, 1);
        amrex::MultiFab tau13(ba13, dm, 1, 1);
        amrex::MultiFab tau23(ba23, dm, 1, 1);
        tau12.setVal(Real(0.0));
        tau13.setVal(Real(0.0));
        tau23.setVal(Real(0.0));

        Boundary walls{};
        walls[2].momentum.model = MomentumModel::NeutralRoughnessLog;
        walls[2].momentum.z0_m = Real(0.01);
        walls[4].momentum.model = MomentumModel::NeutralRoughnessLog;
        walls[4].momentum.z0_m = Real(0.02);
        const amrex::GpuArray<Real, AMREX_SPACEDIM> dx_inv = {
            Real(1.0), Real(1.0), Real(1.0)};
        for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
            erf_cloud_chamber_wall_stress::apply(
                mfi.validbox(), cell_box, state.const_array(mfi),
                u.const_array(mfi), v.const_array(mfi), w.const_array(mfi),
                tau12.array(mfi), tau13.array(mfi), tau23.array(mfi),
                dx_inv, walls);
        }

        amrex::MultiFab rho_u_rhs(amrex::BoxArray(
            amrex::surroundingNodes(cell_box, 0)), dm, 1, 0);
        amrex::MultiFab rho_v_rhs(amrex::BoxArray(
            amrex::surroundingNodes(cell_box, 1)), dm, 1, 0);
        amrex::MultiFab rho_w_rhs(amrex::BoxArray(
            amrex::surroundingNodes(cell_box, 2)), dm, 1, 0);
        rho_u_rhs.setVal(Real(0.0));
        rho_v_rhs.setVal(Real(0.0));
        rho_w_rhs.setVal(Real(0.0));
        amrex::MultiFab tau11(ba, dm, 1, 1);
        amrex::MultiFab tau22(ba, dm, 1, 1);
        amrex::MultiFab tau33(ba, dm, 1, 1);
        amrex::MultiFab tau21(ba12, dm, 1, 1);
        amrex::MultiFab tau31(ba13, dm, 1, 1);
        amrex::MultiFab tau32(ba23, dm, 1, 1);
        tau11.setVal(Real(0.0));
        tau22.setVal(Real(0.0));
        tau33.setVal(Real(0.0));
        tau21.setVal(Real(0.0));
        tau31.setVal(Real(0.0));
        tau32.setVal(Real(0.0));
        amrex::MultiFab detJ(ba, dm, 1, 1);
        amrex::MultiFab mf_mx(ba, dm, 1, 1);
        amrex::MultiFab mf_ux(ba, dm, 1, 1);
        amrex::MultiFab mf_vx(ba, dm, 1, 1);
        amrex::MultiFab mf_my(ba, dm, 1, 1);
        amrex::MultiFab mf_uy(ba, dm, 1, 1);
        amrex::MultiFab mf_vy(ba, dm, 1, 1);
        detJ.setVal(Real(1.0));
        mf_mx.setVal(Real(1.0));
        mf_ux.setVal(Real(1.0));
        mf_vx.setVal(Real(1.0));
        mf_my.setVal(Real(1.0));
        mf_uy.setVal(Real(1.0));
        mf_vy.setVal(Real(1.0));
        amrex::Gpu::DeviceVector<Real> stretched_dz_d;
        DiffusionSrcForMom(
            amrex::surroundingNodes(cell_box, 0),
            amrex::surroundingNodes(cell_box, 1),
            amrex::surroundingNodes(cell_box, 2),
            rho_u_rhs[0].array(), rho_v_rhs[0].array(), rho_w_rhs[0].array(),
            tau11[0].const_array(), tau22[0].const_array(),
            tau33[0].const_array(), tau12[0].const_array(),
            tau21[0].const_array(), tau13[0].const_array(),
            tau31[0].const_array(), tau23[0].const_array(),
            tau32[0].const_array(), detJ[0].const_array(), stretched_dz_d,
            dx_inv, mf_mx[0].const_array(), mf_ux[0].const_array(),
            mf_vx[0].const_array(), mf_my[0].const_array(),
            mf_uy[0].const_array(), mf_vy[0].const_array(), false, false);
        amrex::Gpu::streamSynchronize();
        const amrex::IntVect output_iv(1, 0, 0);
        return value_at(rho_u_rhs, output_iv);
    };

    const Real epsilon = Real(1.0e-5);
    Real row_sum = Real(0.0);
    const auto accumulate = [&](int component, int i_max, int j_max, int k_max) {
        for (int i = 0; i <= i_max; ++i) {
            for (int j = 0; j <= j_max; ++j) {
                for (int k = 0; k <= k_max; ++k) {
                    const amrex::IntVect iv(i, j, k);
                    const Real plus = evaluate(component, iv, epsilon);
                    const Real minus = evaluate(component, iv, -epsilon);
                    row_sum += std::abs((plus - minus) / (Real(2.0) * epsilon));
                }
            }
        }
    };
    accumulate(0, 2, 1, 1);
    accumulate(1, 1, 2, 1);
    accumulate(2, 1, 1, 2);

    const Real cd_y = std::pow(KAPPA / std::log(Real(0.5) / Real(0.01)),
                               Real(2.0));
    const Real cd_z = std::pow(KAPPA / std::log(Real(0.5) / Real(0.02)),
                               Real(2.0));
    const Real bound =
        erf_cloud_chamber_wall_flux::momentum_infinity_row_sum_factor() *
        (cd_y * std::sqrt(u_base*u_base + w_base*w_base) +
         cd_z * std::sqrt(u_base*u_base + v_base*v_base));
    EXPECT_GT(row_sum, Real(0.0));
    const Real tolerance = Real(4096.0) *
        std::numeric_limits<Real>::epsilon() * std::max(Real(1.0), bound);
    EXPECT_LE(row_sum, bound + tolerance);
}

TEST(CloudChamberNeutralLog, CrossWallCompositionAddsDiffusionImpulses)
{
    using namespace erf_wall_thermodynamics;
    const auto run = [](bool ylo_active, bool zlo_active) {
        const amrex::Box cell_box(amrex::IntVect(0), amrex::IntVect(1));
        const amrex::BoxArray ba(cell_box);
        const amrex::DistributionMapping dm(ba);
        amrex::MultiFab state(ba, dm, Rho_comp + 1, 1);
        amrex::BoxArray xba(ba); xba.surroundingNodes(0);
        amrex::BoxArray yba(ba); yba.surroundingNodes(1);
        amrex::BoxArray zba(ba); zba.surroundingNodes(2);
        amrex::MultiFab u(xba, dm, 1, 1);
        amrex::MultiFab v(yba, dm, 1, 1);
        amrex::MultiFab w(zba, dm, 1, 1);
        amrex::BoxArray ba12(ba); ba12.surroundingNodes(0); ba12.surroundingNodes(1);
        amrex::BoxArray ba13(ba); ba13.surroundingNodes(0); ba13.surroundingNodes(2);
        amrex::BoxArray ba23(ba); ba23.surroundingNodes(1); ba23.surroundingNodes(2);
        amrex::MultiFab tau12(ba12, dm, 1, 1);
        amrex::MultiFab tau13(ba13, dm, 1, 1);
        amrex::MultiFab tau23(ba23, dm, 1, 1);
        state.setVal(Real(1.0), Rho_comp, 1);
        u.setVal(Real(1.0));
        v.setVal(Real(0.0));
        w.setVal(Real(0.0));
        tau12.setVal(Real(0.0));
        tau13.setVal(Real(0.0));
        tau23.setVal(Real(0.0));
        Boundary walls{};
        if (ylo_active) {
            walls[2].momentum.model = MomentumModel::NeutralRoughnessLog;
            walls[2].momentum.z0_m = Real(0.01);
        }
        if (zlo_active) {
            walls[4].momentum.model = MomentumModel::NeutralRoughnessLog;
            walls[4].momentum.z0_m = Real(0.02);
        }
        const amrex::GpuArray<Real, AMREX_SPACEDIM> dx_inv = {
            Real(1.0), Real(1.0), Real(1.0)};
        for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
            erf_cloud_chamber_wall_stress::apply(
                mfi.validbox(), cell_box, state.const_array(mfi),
                u.const_array(mfi), v.const_array(mfi), w.const_array(mfi),
                tau12.array(mfi), tau13.array(mfi), tau23.array(mfi),
                dx_inv, walls);
        }
        amrex::Gpu::streamSynchronize();
        amrex::MultiFab rho_u_rhs(amrex::BoxArray(
            amrex::surroundingNodes(cell_box, 0)), dm, 1, 0);
        amrex::MultiFab rho_v_rhs(amrex::BoxArray(
            amrex::surroundingNodes(cell_box, 1)), dm, 1, 0);
        amrex::MultiFab rho_w_rhs(amrex::BoxArray(
            amrex::surroundingNodes(cell_box, 2)), dm, 1, 0);
        rho_u_rhs.setVal(Real(0.0));
        rho_v_rhs.setVal(Real(0.0));
        rho_w_rhs.setVal(Real(0.0));
        amrex::MultiFab tau11(ba, dm, 1, 1);
        amrex::MultiFab tau22(ba, dm, 1, 1);
        amrex::MultiFab tau33(ba, dm, 1, 1);
        amrex::MultiFab tau21(ba12, dm, 1, 1);
        amrex::MultiFab tau31(ba13, dm, 1, 1);
        amrex::MultiFab tau32(ba23, dm, 1, 1);
        tau11.setVal(Real(0.0));
        tau22.setVal(Real(0.0));
        tau33.setVal(Real(0.0));
        tau21.setVal(Real(0.0));
        tau31.setVal(Real(0.0));
        tau32.setVal(Real(0.0));
        amrex::MultiFab detJ(ba, dm, 1, 1);
        amrex::MultiFab mf_mx(ba, dm, 1, 1);
        amrex::MultiFab mf_ux(ba, dm, 1, 1);
        amrex::MultiFab mf_vx(ba, dm, 1, 1);
        amrex::MultiFab mf_my(ba, dm, 1, 1);
        amrex::MultiFab mf_uy(ba, dm, 1, 1);
        amrex::MultiFab mf_vy(ba, dm, 1, 1);
        detJ.setVal(Real(1.0));
        mf_mx.setVal(Real(1.0));
        mf_ux.setVal(Real(1.0));
        mf_vx.setVal(Real(1.0));
        mf_my.setVal(Real(1.0));
        mf_uy.setVal(Real(1.0));
        mf_vy.setVal(Real(1.0));
        amrex::Gpu::DeviceVector<Real> stretched_dz_d;
        DiffusionSrcForMom(
            amrex::surroundingNodes(cell_box, 0),
            amrex::surroundingNodes(cell_box, 1),
            amrex::surroundingNodes(cell_box, 2),
            rho_u_rhs[0].array(), rho_v_rhs[0].array(), rho_w_rhs[0].array(),
            tau11[0].const_array(), tau22[0].const_array(),
            tau33[0].const_array(), tau12[0].const_array(),
            tau21[0].const_array(), tau13[0].const_array(),
            tau31[0].const_array(), tau23[0].const_array(),
            tau32[0].const_array(), detJ[0].const_array(), stretched_dz_d,
            dx_inv, mf_mx[0].const_array(), mf_ux[0].const_array(),
            mf_vx[0].const_array(), mf_my[0].const_array(),
            mf_uy[0].const_array(), mf_vy[0].const_array(), false, false);
        amrex::Gpu::streamSynchronize();
        const amrex::IntVect iv(1, 0, 0);
        return value_at(rho_u_rhs, iv);
    };
    const Real ylo_only = run(true, false);
    const Real zlo_only = run(false, true);
    const Real both_walls = run(true, true);
    EXPECT_LT(ylo_only, Real(0.0));
    EXPECT_LT(zlo_only, Real(0.0));
    EXPECT_NEAR(both_walls, ylo_only + zlo_only,
                scaled_tolerance(both_walls));
}

} // namespace

#include <AMReX_Array.H>
#include "ERF_Constants.H"
#include <AMReX_BoxArray.H>
#include <AMReX_BCUtil.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_RealBox.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "../../../Source/Diffusion/ERF_CloudChamberWallFlux.H"
#include "../../../Source/Diffusion/ERF_Diffusion.H"
#include "../../../Source/Diffusion/ERF_ScalarDiffusion.H"
#include "../../../Source/Diffusion/ERF_CloudChamberWallStress.H"
#include "../../../Source/Diffusion/ERF_ResolvedWallFlux.H"
#include "../../../Source/DataStructs/ERF_DataStruct.H"
#include "../../../Source/Prob/ERF_CloudChamberBudget.H"
#include "../../../Source/Prob/ERF_CloudChamber.H"
#include "../../../Source/Prob/ERF_ProblemDispatch.H"
#include "../../../Source/BoundaryConditions/ERF_WallScalarBC.H"

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

// Host view of the same seven cell-centered lanes passed to estTimeStep's reduction.
Real
rate_at_origin (const erf_wall_thermodynamics::Boundary& walls,
                const GpuArray<Real, AMREX_SPACEDIM>& dx_inv,
                const GpuArray<Real, 7>& state,
                const amrex::IntVect& hi = amrex::IntVect(1))
{
    const amrex::Array4<const Real> velocity(state.data(), amrex::Dim3{0,0,0}, amrex::Dim3{1,1,1}, 7);
    return erf_cloud_chamber_wall_flux::wall_rate_for_cell(0, 0, 0,
        amrex::Box(amrex::IntVect(0), hi), walls, dx_inv, velocity);
}

double
stable_most_residual (double z, double Ri, double D_m, double D_s)
{
    // similarity_funs takes Real, not a double template argument. Keep this
    // oracle in double so float cancellation cannot hide a missed narrow pair.
    const double psi_m =
        -6.1 * std::log(z + std::pow(1.0 + std::pow(z, 2.5), 1.0 / 2.5));
    const double psi_s =
        -5.3 * std::log(z + std::pow(1.0 + std::pow(z, 1.1), 1.0 / 1.1));
    const double A_m = std::max(D_m - psi_m, 1.0);
    const double A_s = std::max(D_s - psi_s, 1.0);
    return z - Ri * A_m * A_m / A_s;
}

void
expect_stable_most_root (Real z, Real Ri, Real D_m, Real D_s, double expected)
{
    // The production bracket is at most 1e-3 wide. Allow that absolute
    // location error in either precision; evaluate the oracle in double so
    // single-precision cancellation cannot masquerade as a small residual.
    constexpr double tolerance = 1.0e-3;
    ASSERT_TRUE(std::isfinite(z));
    EXPECT_NEAR(z, expected, tolerance);
    const double lower = std::max(0.0, double(z) - tolerance);
    const double upper = double(z) + tolerance;
    const double f_lower = stable_most_residual(lower, Ri, D_m, D_s);
    const double f_upper = stable_most_residual(upper, Ri, D_m, D_s);
    EXPECT_LT(f_lower, 0.0);
    EXPECT_GT(f_upper, 0.0);
    EXPECT_LE(std::abs(stable_most_residual(z, Ri, D_m, D_s)),
              std::max(std::abs(f_lower), std::abs(f_upper)));

    similarity_funs sfuns;
    const Real A_m = std::max(D_m - sfuns.calc_psi_m2(z), Real(1.0));
    const Real A_s = std::max(D_s - sfuns.calc_psi_h2(z), Real(1.0));
    const double rhs = double(Ri) * double(A_m) * double(A_m) / double(A_s);
    EXPECT_NEAR(stable_most_residual(z, Ri, D_m, D_s), double(z) - rhs,
                64.0 * std::numeric_limits<Real>::epsilon() *
                    std::max(1.0, std::abs(rhs)));
}

amrex::Gpu::HostVector<Real>
most_coordinates_in_kernel (GpuArray<Real, 5> Ri,
                            GpuArray<Real, 5> D_m,
                            GpuArray<Real, 5> D_s)
{
    amrex::Gpu::DeviceVector<Real> device(5);
    Real* values = device.data();
    amrex::ParallelFor(5, [=] AMREX_GPU_DEVICE(int n) noexcept {
        values[n] = erf_cloud_chamber_wall_flux::most_stability_coordinate(
                        Ri[n], D_m[n], D_s[n])
                        .zeta;
    });
    amrex::Gpu::HostVector<Real> host(5);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, device.begin(), device.end(),
                     host.begin());
    return host;
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
    erf_cloud_chamber::Initialization config;
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
    erf_cloud_chamber::Initialization config;
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

// Motivation: the velocity perturbation is stored on staggered faces, so its
// finite-volume divergence must be checked with face differences rather than
// with the continuous derivative of the analytic expression.
TEST(CloudChamberProfile, StaggeredVelocityIsDiscreteDivergenceFreeOnRectangularGrid)
{
    const GpuArray<Real, AMREX_SPACEDIM> lo = {
        Real(-1.0), Real(2.0), Real(4.0)};
    const GpuArray<Real, AMREX_SPACEDIM> length = {
        Real(3.0), Real(2.0), Real(1.5)};
    const int nx = 5;
    const int ny = 7;
    const int nz = 4;
    const GpuArray<Real, AMREX_SPACEDIM> dx = {
        length[0] / Real(nx), length[1] / Real(ny), length[2] / Real(nz)};
    const Real U_0 = Real(2.5);
    const Real v_scale = erf_cloud_chamber::staggered_velocity_v_scale(
        length[0], length[1], dx[0], dx[1], nx, ny);
    EXPECT_NE(v_scale, Real(1.0));

    Real max_divergence = Real(0.0);
    for (int k = 0; k < nz; ++k) {
        const Real z = lo[2] + (Real(k) + myhalf) * dx[2];
        for (int j = 0; j < ny; ++j) {
            const Real y = lo[1] + (Real(j) + myhalf) * dx[1];
            for (int i = 0; i < nx; ++i) {
                const Real x = lo[0] + (Real(i) + myhalf) * dx[0];
                const Real u_lo = erf_cloud_chamber::initial_u_velocity_perturbation(
                    lo[0] + Real(i) * dx[0], y, z, lo, length, U_0, nx, ny);
                const Real u_hi = erf_cloud_chamber::initial_u_velocity_perturbation(
                    lo[0] + Real(i + 1) * dx[0], y, z, lo, length, U_0, nx, ny);
                const Real v_lo = erf_cloud_chamber::initial_v_velocity_perturbation(
                    x, lo[1] + Real(j) * dx[1], z, lo, length,
                    U_0, v_scale, nx, ny);
                const Real v_hi = erf_cloud_chamber::initial_v_velocity_perturbation(
                    x, lo[1] + Real(j + 1) * dx[1], z, lo, length,
                    U_0, v_scale, nx, ny);
                EXPECT_TRUE(std::isfinite(static_cast<double>(u_lo)));
                EXPECT_TRUE(std::isfinite(static_cast<double>(u_hi)));
                EXPECT_TRUE(std::isfinite(static_cast<double>(v_lo)));
                EXPECT_TRUE(std::isfinite(static_cast<double>(v_hi)));
                const Real divergence = (u_hi - u_lo) / dx[0] +
                    (v_hi - v_lo) / dx[1];
                max_divergence = std::max(max_divergence, std::abs(divergence));
                EXPECT_DOUBLE_EQ(erf_cloud_chamber::initial_w_velocity_perturbation(),
                                 Real(0.0));
            }
        }
    }
    EXPECT_LE(max_divergence,
              scaled_tolerance(U_0 / std::min({dx[0], dx[1]})));

    const Real zero_u = erf_cloud_chamber::initial_u_velocity_perturbation(
        lo[0] + dx[0], lo[1] + myhalf * dx[1], lo[2] + myhalf * dx[2],
        lo, length, Real(0.0), nx, ny);
    const Real zero_v = erf_cloud_chamber::initial_v_velocity_perturbation(
        lo[0] + myhalf * dx[0], lo[1] + dx[1], lo[2] + myhalf * dx[2],
        lo, length, Real(0.0), v_scale, nx, ny);
    EXPECT_DOUBLE_EQ(zero_u, Real(0.0));
    EXPECT_DOUBLE_EQ(zero_v, Real(0.0));

    const GpuArray<Real, AMREX_SPACEDIM> square_length = {
        Real(2.0), Real(2.0), Real(1.0)};
    const GpuArray<Real, AMREX_SPACEDIM> square_dx = {
        Real(0.25), Real(0.25), Real(0.25)};
    EXPECT_DOUBLE_EQ(erf_cloud_chamber::staggered_velocity_v_scale(
        square_length[0], square_length[1], square_dx[0], square_dx[1],
        8, 8), Real(1.0));

    const GpuArray<Real, AMREX_SPACEDIM> degenerate_length = {
        Real(1.0), Real(2.0), Real(1.0)};
    const GpuArray<Real, AMREX_SPACEDIM> degenerate_dx = {
        Real(1.0), Real(0.5), Real(0.25)};
    EXPECT_DOUBLE_EQ(erf_cloud_chamber::staggered_velocity_v_scale(
        degenerate_length[0], degenerate_length[1], degenerate_dx[0],
        degenerate_dx[1], 1, 4), Real(0.0));
    EXPECT_DOUBLE_EQ(erf_cloud_chamber::initial_u_velocity_perturbation(
        lo[0], lo[1] + myhalf * degenerate_dx[1],
        lo[2] + myhalf * degenerate_dx[2], lo, degenerate_length,
        U_0, 1, 4), Real(0.0));
    EXPECT_DOUBLE_EQ(erf_cloud_chamber::initial_v_velocity_perturbation(
        lo[0] + myhalf * degenerate_dx[0], lo[1],
        lo[2] + myhalf * degenerate_dx[2], lo, degenerate_length,
        U_0, Real(0.0), 1, 4), Real(0.0));
}

TEST(CloudChamberBudget, CloudyRhoThetaStatusPrecedesClosureResult)
{
    EXPECT_STREQ(CloudChamberBudget::row_status(
                     CloudChamberBudget::RhoTheta, true, true),
                 "UNSUPPORTED_SOURCE");
    EXPECT_STREQ(CloudChamberBudget::row_status(
                     CloudChamberBudget::RhoTheta, true, false),
                 "UNSUPPORTED_SOURCE");
    EXPECT_STREQ(CloudChamberBudget::row_status(
                     CloudChamberBudget::RhoTheta, false, true), "PASS");
    EXPECT_STREQ(CloudChamberBudget::row_status(
                     CloudChamberBudget::RhoTheta, false, false), "FAIL");
    EXPECT_STREQ(CloudChamberBudget::row_status(
                     CloudChamberBudget::RhoQv, true, true), "PASS");
    EXPECT_STREQ(CloudChamberBudget::row_status(
                     CloudChamberBudget::RhoQc, true, false), "FAIL");
}

// Motivation: physical chamber inputs are specified as temperature and RH;
// the initializer must convert RH to the exact mixing ratio used by the
// anelastic thermodynamics, rather than treating temperature as theta.
TEST(CloudChamberProfile, PhysicalTemperatureAndRelativeHumidityAreExact)
{
    erf_cloud_chamber::Initialization config;
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

TEST(CloudChamberConfig, MoistureModelSupportIsAnExplicitTypedAllowlist)
{
    const std::vector<std::pair<MoistureType, bool>> expected = {
        {MoistureType::None, true},
        {MoistureType::SatAdj, true},
        {MoistureType::MoistNoCondensation, true},
        {MoistureType::Kessler_NoRain, true},
        {MoistureType::Kessler, true},
        {MoistureType::SAM_NoPrecip_NoIce, true},
        {MoistureType::SAM_NoIce, true},
        {MoistureType::SAM, true},
        {MoistureType::WSM6, true},
        {MoistureType::WDM6, true},
        {MoistureType::Morrison_NoIce, true},
        {MoistureType::Morrison, true},
        {MoistureType::SBM, false},
        {MoistureType::SuperDroplets, true}};

    for (const auto& [model, supported] : expected) {
        EXPECT_EQ(erf_cloud_chamber::supports_moisture_model(model), supported)
            << amrex::getEnumNameString(model);
    }
}

// Keep AMReX GPU lambda launches in namespace-scope helpers: NVCC cannot
// compile extended device lambdas nested in GoogleTest's generated TestBody.
void
initialize_extra_moisture_lanes (amrex::Array4<amrex::Real> values,
                                 const amrex::Box& box)
{
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real base = amrex::Real(0.1) * i + amrex::Real(0.2) * j +
            amrex::Real(0.3) * k;
        values(i,j,k,Rho_comp) = amrex::Real(1.0);
        values(i,j,k,RhoQ3_comp) = base;
        values(i,j,k,RhoQ4_comp) = amrex::Real(1.5) * base + amrex::Real(0.09) * k;
        values(i,j,k,RhoQ7_comp) = amrex::Real(2.0) * base + amrex::Real(0.07) * i;
        values(i,j,k,RhoQ8_comp) = amrex::Real(3.0) * base + amrex::Real(0.11) * j;
    });
}

void
evaluate_extra_moisture_lanes_flux (
    amrex::Array4<const amrex::Real> values,
    amrex::Array4<amrex::Real> wall,
    amrex::Array4<amrex::Real> interior,
    amrex::Array4<amrex::Real> seam,
    amrex::Array4<amrex::Real> mutant,
    const amrex::Box& box,
    amrex::GpuArray<int, 4> q_comps,
    int lo_x, int hi_x, int lo_y, int hi_y, int lo_z, int hi_z)
{
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real wall_max[4] = {amrex::Real(0.0), amrex::Real(0.0),
                                   amrex::Real(0.0), amrex::Real(0.0)};
        amrex::Real interior_max[4] = {amrex::Real(0.0), amrex::Real(0.0),
                                       amrex::Real(0.0), amrex::Real(0.0)};
        amrex::Real seam_max[4] = {amrex::Real(0.0), amrex::Real(0.0),
                                   amrex::Real(0.0), amrex::Real(0.0)};
        for (int n = 0; n < 4; ++n) {
            const int comp = q_comps[n];
            if (i == lo_x) {
                wall_max[n] = amrex::max(wall_max[n], std::abs(
                    ScalarDiffusionFlux_N<0>(values, comp, lo_x, j, k,
                        amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.0),
                        amrex::Real(1.0), false, false)));
            }
            if (i == hi_x) {
                wall_max[n] = amrex::max(wall_max[n], std::abs(
                    ScalarDiffusionFlux_N<0>(values, comp, hi_x + 1, j, k,
                        amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.0),
                        amrex::Real(1.0), false, false)));
            }
            if (j == lo_y) {
                wall_max[n] = amrex::max(wall_max[n], std::abs(
                    ScalarDiffusionFlux_N<1>(values, comp, i, lo_y, k,
                        amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.0),
                        amrex::Real(1.0), false, false)));
            }
            if (j == hi_y) {
                wall_max[n] = amrex::max(wall_max[n], std::abs(
                    ScalarDiffusionFlux_N<1>(values, comp, i, hi_y + 1, k,
                        amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.0),
                        amrex::Real(1.0), false, false)));
            }
            if (k == lo_z) {
                wall_max[n] = amrex::max(wall_max[n], std::abs(
                    ScalarDiffusionFlux_N<2>(values, comp, i, j, lo_z,
                        amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.0),
                        amrex::Real(1.0), false, false)));
            }
            if (k == hi_z) {
                wall_max[n] = amrex::max(wall_max[n], std::abs(
                    ScalarDiffusionFlux_N<2>(values, comp, i, j, hi_z + 1,
                        amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.0),
                        amrex::Real(1.0), false, false)));
            }
            if (i == 3 && j == 3 && k == 3) {
                interior_max[n] = amrex::max(
                    std::abs(ScalarDiffusionFlux_N<0>(values, comp, 3, 3, 3,
                        amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.0),
                        amrex::Real(1.0), false, false)),
                    amrex::max(
                        std::abs(ScalarDiffusionFlux_N<1>(values, comp, 3, 3, 3,
                            amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.0),
                            amrex::Real(1.0), false, false)),
                        std::abs(ScalarDiffusionFlux_N<2>(values, comp, 3, 3, 3,
                            amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.0),
                            amrex::Real(1.0), false, false))));
            }
            if (i == 3 || i == 4) {
                seam_max[n] = amrex::max(seam_max[n], std::abs(
                    ScalarDiffusionFlux_N<0>(values, comp, 4, j, k,
                        amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.0),
                        amrex::Real(1.0), false, false)));
            }
            if (j == 3 || j == 4) {
                seam_max[n] = amrex::max(seam_max[n], std::abs(
                    ScalarDiffusionFlux_N<1>(values, comp, i, 4, k,
                        amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.0),
                        amrex::Real(1.0), false, false)));
            }
            if (k == 3 || k == 4) {
                seam_max[n] = amrex::max(seam_max[n], std::abs(
                    ScalarDiffusionFlux_N<2>(values, comp, i, j, 4,
                        amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.0),
                        amrex::Real(1.0), false, false)));
            }
        }
        for (int n = 0; n < 4; ++n) {
            wall(i,j,k,n) = wall_max[n];
            interior(i,j,k,n) = interior_max[n];
            seam(i,j,k,n) = seam_max[n];
        }
        if (i == lo_x && j == 3 && k == 3) {
            mutant(i,j,k) = std::abs(ScalarDiffusionFlux_N<0>(
                values, RhoQ3_comp, lo_x, j, k, amrex::Real(1.0),
                amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.0), true, false));
        }
    });
}

void
initialize_scalar_adapter_fields (amrex::Array4<amrex::Real> conserved,
                                  amrex::Array4<amrex::Real> intensive,
                                  const amrex::Box& data_box,
                                  const amrex::Box& domain)
{
    amrex::ParallelFor(data_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const int ii = amrex::max(domain.smallEnd(0), amrex::min(domain.bigEnd(0), i));
        const int jj = amrex::max(domain.smallEnd(1), amrex::min(domain.bigEnd(1), j));
        const int kk = amrex::max(domain.smallEnd(2), amrex::min(domain.bigEnd(2), k));
        const amrex::Real base = amrex::Real(0.1) * ii + amrex::Real(0.2) * jj +
            amrex::Real(0.3) * kk;
        conserved(i,j,k,Rho_comp) = amrex::Real(1.0);
        intensive(i,j,k,RhoQ3_comp) = base;
        intensive(i,j,k,RhoQ4_comp) = amrex::Real(2.0) * base + amrex::Real(0.07) * ii;
        intensive(i,j,k,RhoQ7_comp) = amrex::Real(3.0) * base + amrex::Real(0.11) * jj;
        intensive(i,j,k,RhoQ8_comp) = amrex::Real(4.0) * base + amrex::Real(0.13) * kk;
    });
}

// Regression for the generic scalar path used by Q3+ in a physical chamber.
// The real multibox fill applies the same foextrap scalar contract assigned to
// NoSlipWall constituents; the native N-grid face kernel must then write zero
// wall flux while retaining nonzero interior transport for mass and number
// lanes. The final ext_dir control is a deliberate wrong-mapping negative
// control: it makes this exact stencil report a nonzero wall flux.
TEST(CloudChamberWallDiffusion, ExtraMoistureLanesAreImpermeable)
{
    using namespace amrex;
    const Box domain(IntVect(0, 0, 0), IntVect(7, 7, 7));
    const RealBox real_box({Real(0.0), Real(0.0), Real(0.0)},
                           {Real(8.0), Real(8.0), Real(8.0)});
    const int periodicity[AMREX_SPACEDIM] = {0, 0, 0};
    const Geometry geom(domain, &real_box, CoordSys::cartesian, periodicity);
    BoxArray grids(domain);
    grids.maxSize(4);
    DistributionMapping dmap(grids);

    MultiFab state(grids, dmap, NVAR_max, 1);
    state.setVal(Real(-99.0));
    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        const auto values = state.array(mfi);
        const Box bx = mfi.validbox();
        initialize_extra_moisture_lanes(values, bx);
    }
    Vector<BCRec> bcs(state.nComp());
    for (auto& bc : bcs) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            bc.setLo(dir, ERFBCType::foextrap);
            bc.setHi(dir, ERFBCType::foextrap);
        }
    }
    state.FillBoundary(geom.periodicity());
    FillDomainBoundary(state, geom, bcs);

    MultiFab wall_flux(grids, dmap, 4, 0);
    MultiFab interior_flux(grids, dmap, 4, 0);
    MultiFab seam_flux(grids, dmap, 4, 0);
    MultiFab ext_dir_mutant_flux(grids, dmap, 1, 0);
    wall_flux.setVal(Real(0.0));
    interior_flux.setVal(Real(0.0));
    seam_flux.setVal(Real(0.0));
    ext_dir_mutant_flux.setVal(Real(0.0));
    const GpuArray<int, 4> q_comps = {RhoQ3_comp, RhoQ4_comp, RhoQ7_comp, RhoQ8_comp};
    const int lo_x = domain.smallEnd(0), hi_x = domain.bigEnd(0);
    const int lo_y = domain.smallEnd(1), hi_y = domain.bigEnd(1);
    const int lo_z = domain.smallEnd(2), hi_z = domain.bigEnd(2);
    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        const auto values = state.const_array(mfi);
        const auto wall = wall_flux.array(mfi);
        const auto interior = interior_flux.array(mfi);
        const auto seam = seam_flux.array(mfi);
        const auto mutant = ext_dir_mutant_flux.array(mfi);
        const Box bx = mfi.validbox();
        evaluate_extra_moisture_lanes_flux(
            values, wall, interior, seam, mutant, bx, q_comps,
            lo_x, hi_x, lo_y, hi_y, lo_z, hi_z);
    }
    Gpu::streamSynchronize();
    for (int n = 0; n < 4; ++n) {
        const Real wall_max = wall_flux.norm0(n, 0, false);
        const Real interior_max = interior_flux.norm0(n, 0, false);
        const Real seam_max = seam_flux.norm0(n, 0, false);
        std::cout << std::setprecision(17)
                  << "cloud_chamber_extra_q_stencil component=" << q_comps[n]
                  << " wall_flux_max=" << wall_max
                  << " interior_flux_max=" << interior_max
                  << " seam_flux_max=" << seam_max << "\n";
        EXPECT_DOUBLE_EQ(wall_max, Real(0.0))
            << "wall flux leaked from Q" << q_comps[n];
        EXPECT_GT(interior_max, Real(0.0))
            << "interior gradient was vacuous for Q" << q_comps[n];
        EXPECT_GT(seam_max, Real(0.0))
            << "inter-box seam flux was vacuous for Q" << q_comps[n];
    }
    EXPECT_GT(ext_dir_mutant_flux.norm0(0, 0, false), Real(0.0))
        << "wrong ext_dir wall mapping did not activate the boundary stencil";
}

// This test feeds the same NoSlipWall intent-to-BC mapper used by init_bcs
// into BuildScalarDiffusionFluxes_N, requiring exact zero wall flux, positive
// interior flux, and a closed-domain divergence sum within 512 eps. The
// multibox stencil regression above independently requires nonzero seam flux.
// Corrupting the shared production mapping to ext_dir was tested and failed
// this oracle. This is the narrow adapter-level equivalent; it does not run
// ERF::init_bcs or erf_slow_rhs_post end to end. DOUBLE precision.
TEST(CloudChamberWallDiffusion, NoSlipMappingUsesProductionScalarAdapter)
{
    using namespace amrex;
    const Box domain(IntVect(0, 0, 0), IntVect(7, 7, 7));
    Box data_box = domain;
    data_box.grow(1);
    const Box xfaces = surroundingNodes(domain, 0);
    const Box yfaces = surroundingNodes(domain, 1);
    const Box zfaces = surroundingNodes(domain, 2);
    FArrayBox state(data_box, NVAR_max);
    FArrayBox primitive(data_box, NVAR_max);
    FArrayBox mu(data_box, EddyDiff::NumDiffs);
    FArrayBox rhs(domain, NVAR_max);
    FArrayBox xflux(xfaces, 1), yflux(yfaces, 1), zflux(zfaces, 1);
    FArrayBox u(xfaces, 1), v(yfaces, 1);
    FArrayBox mf_ux(xfaces, 1), mf_uy(xfaces, 1);
    FArrayBox mf_vx(yfaces, 1), mf_vy(yfaces, 1);

    state.setVal<RunOn::Device>(Real(0.0));
    primitive.setVal<RunOn::Device>(Real(0.0));
    mu.setVal<RunOn::Device>(Real(0.0));
    rhs.setVal<RunOn::Device>(Real(0.0));
    u.setVal<RunOn::Device>(Real(0.0));
    v.setVal<RunOn::Device>(Real(0.0));
    mf_ux.setVal<RunOn::Device>(Real(1.0));
    mf_uy.setVal<RunOn::Device>(Real(1.0));
    mf_vx.setVal<RunOn::Device>(Real(1.0));
    mf_vy.setVal<RunOn::Device>(Real(1.0));

    const auto conserved = state.array();
    const auto intensive = primitive.array();
    initialize_scalar_adapter_fields(conserved, intensive, data_box, domain);
    Gpu::streamSynchronize();

    Vector<BCRec> bcs(NVAR_max);
    const int no_slip_unspecified = erf_wall_scalar_bc::logical_boundary_type(
        erf_wall_scalar_bc::WallScalarBCIntent::Unspecified);
    EXPECT_EQ(no_slip_unspecified, ERFBCType::foextrap);
    for (auto& bc : bcs) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            bc.setLo(dir, no_slip_unspecified);
            bc.setHi(dir, no_slip_unspecified);
        }
    }
    Gpu::DeviceVector<BCRec> bcs_device(bcs.size());
    Gpu::copy(Gpu::hostToDevice, bcs.begin(), bcs.end(), bcs_device.begin());

    ScalarDiffusionFieldViews field;
    field.scalar = primitive.const_array();
    field.density = state.const_array();
    field.rho_comp = Rho_comp;
    field.mu_turb = mu.const_array();
    field.xflux = xflux.array();
    field.yflux = yflux.array();
    field.zflux = zflux.array();
    field.rhs = rhs.array();
    ScalarDiffusionFluxPolicy policy;
    policy.coefficients = {Real(0.8), EddyDiff::Q_h, EddyDiff::Q_v};
    policy.coefficient_mode = {true, false};
    const GpuArray<Real, AMREX_SPACEDIM> inv{{Real(1.0), Real(1.0), Real(1.0)}};
    const int q_components[] = {RhoQ3_comp, RhoQ4_comp, RhoQ7_comp, RhoQ8_comp};

    auto copy_to_host = [](const FArrayBox& src, FArrayBox& dst) {
        Gpu::copy(Gpu::deviceToHost, src.dataPtr(0), src.dataPtr(0) + src.size(),
                  dst.dataPtr(0));
        Gpu::streamSynchronize();
    };
    for (const int q_comp : q_components) {
        field.scalar_comp = q_comp;
        xflux.setVal<RunOn::Device>(Real(0.0));
        yflux.setVal<RunOn::Device>(Real(0.0));
        zflux.setVal<RunOn::Device>(Real(0.0));
        BuildScalarDiffusionFluxes_N(
            domain, domain, field, policy, inv, u.const_array(), v.const_array(),
            mf_ux.const_array(), mf_uy.const_array(), mf_vx.const_array(),
            mf_vy.const_array(), bcs_device.data(), q_comp);
        Gpu::streamSynchronize();

        FArrayBox xhost(xfaces, 1, The_Pinned_Arena());
        FArrayBox yhost(yfaces, 1, The_Pinned_Arena());
        FArrayBox zhost(zfaces, 1, The_Pinned_Arena());
        copy_to_host(xflux, xhost);
        copy_to_host(yflux, yhost);
        copy_to_host(zflux, zhost);
        const auto fx = xhost.const_array();
        const auto fy = yhost.const_array();
        const auto fz = zhost.const_array();
        Real wall_max = Real(0.0);
        Real interior_max = Real(0.0);
        Real integrated_divergence = Real(0.0);
        for (int k = domain.smallEnd(2); k <= domain.bigEnd(2); ++k) {
            for (int j = domain.smallEnd(1); j <= domain.bigEnd(1); ++j) {
                wall_max = std::max(wall_max, std::abs(fx(domain.smallEnd(0),j,k,0)));
                wall_max = std::max(wall_max, std::abs(fx(domain.bigEnd(0)+1,j,k,0)));
                for (int i = domain.smallEnd(0)+1; i <= domain.bigEnd(0); ++i) {
                    interior_max = std::max(interior_max, std::abs(fx(i,j,k,0)));
                }
            }
        }
        for (int k = domain.smallEnd(2); k <= domain.bigEnd(2); ++k) {
            for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
                wall_max = std::max(wall_max, std::abs(fy(i,domain.smallEnd(1),k,0)));
                wall_max = std::max(wall_max, std::abs(fy(i,domain.bigEnd(1)+1,k,0)));
                for (int j = domain.smallEnd(1)+1; j <= domain.bigEnd(1); ++j) {
                    interior_max = std::max(interior_max, std::abs(fy(i,j,k,0)));
                }
            }
        }
        for (int j = domain.smallEnd(1); j <= domain.bigEnd(1); ++j) {
            for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
                wall_max = std::max(wall_max, std::abs(fz(i,j,domain.smallEnd(2),0)));
                wall_max = std::max(wall_max, std::abs(fz(i,j,domain.bigEnd(2)+1,0)));
                for (int k = domain.smallEnd(2)+1; k <= domain.bigEnd(2); ++k) {
                    interior_max = std::max(interior_max, std::abs(fz(i,j,k,0)));
                }
            }
        }
        for (int k = domain.smallEnd(2); k <= domain.bigEnd(2); ++k) {
            for (int j = domain.smallEnd(1); j <= domain.bigEnd(1); ++j) {
                for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
                    integrated_divergence += fx(i+1,j,k,0) - fx(i,j,k,0) +
                        fy(i,j+1,k,0) - fy(i,j,k,0) +
                        fz(i,j,k+1,0) - fz(i,j,k,0);
                }
            }
        }
        EXPECT_DOUBLE_EQ(wall_max, Real(0.0)) << "physical wall flux for Q" << q_comp;
        EXPECT_GT(interior_max, Real(0.0)) << "interior flux for Q" << q_comp;
        EXPECT_NEAR(integrated_divergence, Real(0.0),
                    Real(512.0) * std::numeric_limits<Real>::epsilon())
            << "closed-domain integrated flux balance for Q" << q_comp;
        std::cout << std::setprecision(17)
                  << "cloud_chamber_extra_q_adapter component=" << q_comp
                  << " wall_flux_max=" << wall_max
                  << " interior_flux_max=" << interior_max
                  << " integrated_divergence=" << integrated_divergence
                  << " balance_tolerance="
                  << Real(512.0) * std::numeric_limits<Real>::epsilon() << "\n";
    }
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
    EXPECT_NE(error.find("active moisture model"), std::string::npos);
}

// Motivation: physical initialization with any active moisture model needs
// an explicit RH rather than an implicit default.
TEST(CloudChamberConfig, RequiresRelativeHumidityForPhysicalMoisture)
{
    const erf_cloud_chamber::InitializationContract contract {
        erf_cloud_chamber::InitializationMode::PhysicalTemperatureRH,
        true,
        false,
        false,
        false};

    const auto error = erf_cloud_chamber::initialization_contract_error(contract);

    EXPECT_NE(error.find("physical moist"), std::string::npos);
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

// Motivation: the valid physical moist contract must remain accepted while
// the stricter dry/moist key checks are enforced.
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

    // init_params writes the static members of SolverChoice (init_type, and
    // mesh_type = VariableDz under a fitted terrain_type); restore them so the
    // tests that follow in a shuffled order see the defaults this test found
    struct ScopedSolverChoiceStatics {
        InitType      init_type      = SolverChoice::init_type;
        SoundingType  sounding_type  = SolverChoice::sounding_type;
        TerrainType   terrain_type   = SolverChoice::terrain_type;
        bool          flat_terrain   = SolverChoice::flat_terrain;
        BuildingsType buildings_type = SolverChoice::buildings_type;
        bool          use_real_bcs   = SolverChoice::use_real_bcs;
        MeshType      mesh_type      = SolverChoice::mesh_type;
        ~ScopedSolverChoiceStatics () {
            SolverChoice::init_type      = init_type;
            SolverChoice::sounding_type  = sounding_type;
            SolverChoice::terrain_type   = terrain_type;
            SolverChoice::flat_terrain   = flat_terrain;
            SolverChoice::buildings_type = buildings_type;
            SolverChoice::use_real_bcs   = use_real_bcs;
            SolverChoice::mesh_type      = mesh_type;
        }
    } restore_statics;

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

TEST(CloudChamberWallConfig, ParsesWetnessOnlyForWetPhysicalWalls)
{
    using namespace erf_wall_thermodynamics;
    amrex::ParmParse pp("cc_test_wall_wetness");
    FaceWall wall;
    wall.moisture = MoistureMode::WetEquilibrium;
    auto parse = [&] { return erf_cloud_chamber::parse_wall_wetness(pp, "zlo", wall); };
    EXPECT_TRUE(parse().empty());
    EXPECT_EQ(wall.wetness, Real(1));
    for (Real rh : {Real(0), Real(0.625), Real(1)}) {
        pp.add("wetness", rh);
        EXPECT_TRUE(parse().empty());
        EXPECT_EQ(wall.wetness, rh);
        pp.remove("wetness");
    }
    for (Real rh : {Real(-0.01), Real(1.01)}) {
        pp.add("wetness", rh);
        EXPECT_NE(parse().find("zlo.wetness"), std::string::npos);
        pp.remove("wetness");
    }
    // Test the numeric parser guard independently of how ParmParse reads NaN/Inf tokens.
    const Real inf = std::numeric_limits<Real>::infinity();
    for (Real rh : {std::numeric_limits<Real>::quiet_NaN(), inf, -inf}) {
        EXPECT_FALSE(erf_cloud_chamber::valid_wall_wetness(rh));
    }
    for (auto mode : {MoistureMode::DryImpermeable, MoistureMode::LegacyNumerical}) {
        wall.moisture = mode;
        EXPECT_TRUE(parse().empty());
        pp.add("wetness", Real(1)); // Even the default is invalid when explicitly supplied.
        EXPECT_NE(parse().find("zlo.wetness"), std::string::npos);
        pp.remove("wetness");
    }
}

class CloudChamberSidewallBands : public ::testing::Test {
protected:
    static amrex::Geometry shifted_geometry ()
    {
        const amrex::Box domain(amrex::IntVect(0,0,7), amrex::IntVect(1,1,10));
        const amrex::RealBox box({0,0,11}, {1,1,13});
        const int periodic[AMREX_SPACEDIM] = {0,0,0};
        return amrex::Geometry(domain, &box, amrex::CoordSys::cartesian, periodic);
    }

    const amrex::Geometry geom = shifted_geometry();
    erf_cloud_chamber::SidewallTemperatureBands bands;

    std::string define_profile ()
    {
        // Heights above the floor, not model z. dz = 1/2 m; one edge lands
        // on a row boundary, and the narrowest panel is only 1/64 m high.
        return bands.define({0, 0.625, 1, 1.015625, 1.0625, 1.25, 2},
                            {285, 295, 308, 281, 300, 290}, geom);
    }

    // Independent hand-computed (Kelvin, area fraction) pairs for k = 7..10.
    const std::vector<std::vector<std::pair<Real, Real>>> rows{
        {{285, 1}}, {{285, 0.25}, {295, 0.75}},
        {{308, 1./32}, {281, 3./32}, {300, 3./8}, {290, 0.5}}, {{290, 1}}};
};

TEST_F(CloudChamberSidewallBands, RejectsInvalidProfileData)
{
    auto reject = [&](const char* reason, amrex::Vector<Real> edges,
                      const amrex::Vector<Real>& values) {
        erf_cloud_chamber::SidewallTemperatureBands invalid;
        EXPECT_FALSE(invalid.define(edges, values, geom).empty()) << reason;
    };
    reject("no panels", {0}, {});
    reject("missing edge", {0, 2}, {285, 295});
    reject("missing value", {0, 1, 2}, {285});
    reject("gap at floor", {0.25, 2}, {285});
    reject("gap at ceiling", {0, 1.75}, {285});
    reject("extends past ceiling", {0, 2.25}, {285});
    reject("duplicate edge", {0, 1, 1, 2}, {285, 290, 295});
    reject("decreasing edges", {0, 1.5, 1, 2}, {285, 290, 295});
    const Real eps = std::numeric_limits<Real>::epsilon();
    reject("endpoint tolerance cannot hide a reversed first band", {2*eps, eps, 2}, {285, 295});
    reject("interior edge must be above the floor", {-2*eps, -eps, 2}, {285, 295});
    for (Real bad : {std::numeric_limits<Real>::quiet_NaN(),
                     std::numeric_limits<Real>::infinity()}) {
        reject("nonfinite edge", {0, bad, 2}, {285, 295});
        reject("nonfinite Kelvin", {0, 2}, {bad});
    }
    reject("zero Kelvin", {0, 2}, {0});
    reject("negative Kelvin", {0, 2}, {-1});
}

TEST_F(CloudChamberSidewallBands, ParserRequiresAPairOnAPhysicalSidewall)
{
    amrex::ParmParse pp("cc_test_sidewall_bands");
    erf_wall_thermodynamics::FaceWall wall;
    auto parse = [&](bool physical, bool sidewall) {
        return erf_cloud_chamber::parse_wall_heat(pp, sidewall ? "yhi" : "zlo", physical,
                                                   wall, geom, sidewall ? &bands : nullptr);
    };
    pp.addarr("temperature_height_edges", amrex::Vector<Real>{0, 1, 2});
    EXPECT_NE(parse(true, true).find("together"), std::string::npos);
    pp.addarr("temperature_values", amrex::Vector<Real>{285, 295});
    EXPECT_NE(parse(false, true).find("physical"), std::string::npos);
    EXPECT_NE(parse(true, false).find("sidewall"), std::string::npos);
    pp.add("temperature", Real(290));
    EXPECT_NE(parse(true, true).find("mutually exclusive"), std::string::npos);
    pp.remove("temperature");
    EXPECT_TRUE(parse(true, true).empty());
    pp.remove("temperature_height_edges");
    EXPECT_NE(parse(true, true).find("together"), std::string::npos);
    pp.remove("temperature_values");
}

TEST_F(CloudChamberSidewallBands, DeviceRowsPartitionTheShiftedWall)
{
    using namespace amrex;
    ASSERT_EQ(define_profile(), "");
    const auto view = bands.view();
    static_assert(std::is_trivially_copyable_v<erf_wall_thermodynamics::TemperatureBandView>);
    ASSERT_EQ(view.k_lo, geom.Domain().smallEnd(2));
    ASSERT_NE(view.offsets, nullptr);
    ASSERT_NE(view.panels, nullptr);
    Gpu::HostVector<int> offsets(rows.size() + 1);
    Gpu::copy(Gpu::deviceToHost, view.offsets, view.offsets + offsets.size(), offsets.begin());
    const int base = offsets.front();
    ASSERT_GE(base, 0);
    for (std::size_t r = 0; r < rows.size(); ++r) {
        ASSERT_EQ(offsets[r+1] - offsets[r], static_cast<int>(rows[r].size())) << "row " << r;
    }
    Gpu::HostVector<erf_wall_thermodynamics::TemperaturePanel> panels(offsets.back() - base);
    Gpu::copy(Gpu::deviceToHost, view.panels + base, view.panels + offsets.back(), panels.begin());
    for (std::size_t r = 0; r < rows.size(); ++r) {
        SCOPED_TRACE(::testing::Message() << "k=" << view.k_lo + static_cast<int>(r));
        const auto first = panels.begin() + (offsets[r] - base);
        const auto last = panels.begin() + (offsets[r+1] - base);
        Real sum = 0;
        for (auto it = first; it != last; ++it) { sum += it->fraction; }
        EXPECT_NEAR(sum, 1, scaled_tolerance(1));
        for (const auto& expected : rows[r]) {
            const auto it = std::find_if(first, last, [&](const auto& panel) {
                return panel.temperature_K == expected.first;
            });
            ASSERT_NE(it, last) << "missing panel at " << expected.first << " K";
            EXPECT_NEAR(it->fraction, expected.second, scaled_tolerance(1));
        }
    }
}

struct BandFluxPair { Real heat, vapor; };

// Keep the launch outside GoogleTest's generated TestBody (NVCC extended lambdas).
amrex::Gpu::HostVector<BandFluxPair>
band_fluxes_in_kernel (erf_wall_thermodynamics::FaceWall wall,
                       erf_cloud_chamber_wall_flux::ScalarWallSample sample,
                       erf_wall_thermodynamics::TemperatureBandView view,
                       const amrex::Box& domain)
{
    using namespace erf_cloud_chamber_wall_flux;
    const int lo = domain.smallEnd(2), nz = domain.length(2);
    amrex::Gpu::DeviceVector<BandFluxPair> device(nz);
    auto* out = device.data();
    amrex::ParallelFor(nz, [=] AMREX_GPU_DEVICE (int r) noexcept {
        auto s = sample;
        s.scalar_air = s.theta_air;
        out[r].heat = evaluate_face_scalar_flux_in(wall, ScalarChannel::Heat, s, view, lo+r);
        s.scalar_air = s.qv_air;
        out[r].vapor = evaluate_face_scalar_flux_in(wall, ScalarChannel::Vapor, s, view, lo+r);
    });
    amrex::Gpu::HostVector<BandFluxPair> host(nz);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, device.begin(), device.end(), host.begin());
    return host;
}

TEST_F(CloudChamberSidewallBands, AveragesPanelFluxesNotWallTemperatures)
{
    using namespace erf_wall_thermodynamics;
    ASSERT_EQ(define_profile(), "");
    const auto view = bands.view();
    ASSERT_EQ(view.k_lo, geom.Domain().smallEnd(2));
    ASSERT_NE(view.offsets, nullptr);
    ASSERT_NE(view.panels, nullptr);
    auto uniform = view;
    uniform.offsets = nullptr; // This alone selects scalar fallback; panels must be ignored.

    FaceWall wall{};
    const Real mean_T = Real(0.25)*285 + Real(0.75)*295;
    wall.thermal = {ThermalMode::FixedPhysicalTemperature, mean_T};
    wall.moisture = MoistureMode::WetEquilibrium;
    wall.wetness = Real(0.625);
    erf_cloud_chamber_wall_flux::ScalarWallSample s{};
    s.rho = Real(1.25); s.p_hse = 90000; s.rdOcp = R_d/Cp_d; s.theta_air = 298;
    s.U_t = 3; s.dx_inv = 2; s.wall_distance = Real(0.25);
    s.alpha_T = Real(0.0625); s.alpha_C = Real(0.125);
    const Real to_theta = std::pow(p_0 / s.p_hse, s.rdOcp);
    auto qwall = [&](Real T) {
        const Real e_Pa = wall.wetness * (100 * erf_esatw(T));
        return RdoRv * e_Pa / (s.p_hse - e_Pa);
    };
    // In the 25/75 row F(mean T) gives condensation, but the area-weighted
    // panel flux gives evaporation. Both use the same local air and pressure.
    s.qv_air = (qwall(mean_T) + Real(0.25)*qwall(285) + Real(0.75)*qwall(295))/2;

    for (auto model : {ScalarModel::ResolvedMolecular, ScalarModel::BulkAero}) {
        const bool bulk = model == ScalarModel::BulkAero;
        SCOPED_TRACE(bulk ? "fixed bulk" : "resolved");
        const auto provider = bulk ? CoefficientProvider::Fixed : CoefficientProvider::None;
        wall.heat = {model, provider, Real(0.125)};
        wall.vapor = {model, provider, Real(0.0625)};
        const Real h = s.rho * (bulk ? wall.heat.coefficient*s.U_t : 2*s.alpha_T*s.dx_inv);
        const Real v = s.rho * (bulk ? wall.vapor.coefficient*s.U_t : 2*s.alpha_C*s.dx_inv);
        const auto heat = [&](Real T) { return h * (T * to_theta - s.theta_air); };
        const auto vapor = [&](Real T) { return v * (qwall(T) - s.qv_air); };
        const Real heat_tol = scaled_tolerance(h * s.theta_air);
        const Real vapor_tol = 64 * std::numeric_limits<Real>::epsilon() * v * qwall(308);

        const auto wet = band_fluxes_in_kernel(wall, s, view, geom.Domain());
        const auto scalar = band_fluxes_in_kernel(wall, s, uniform, geom.Domain());
        auto dry_wall = wall;
        dry_wall.moisture = MoistureMode::DryImpermeable;
        const auto dry = band_fluxes_in_kernel(dry_wall, s, view, geom.Domain());
        auto adiabatic_wall = wall;
        adiabatic_wall.thermal.mode = ThermalMode::Adiabatic;
        const auto adiabatic = band_fluxes_in_kernel(adiabatic_wall, s, view, geom.Domain());

        EXPECT_GT(wet[1].vapor, vapor_tol);
        EXPECT_LT(scalar[1].vapor, -vapor_tol);
        for (std::size_t r = 0; r < rows.size(); ++r) {
            SCOPED_TRACE(::testing::Message() << "k=" << view.k_lo + static_cast<int>(r));
            Real expected_heat = 0, expected_vapor = 0;
            for (const auto& [T, fraction] : rows[r]) {
                expected_heat += fraction * heat(T);
                expected_vapor += fraction * vapor(T);
            }
            EXPECT_NEAR(wet[r].heat, expected_heat, heat_tol);
            EXPECT_NEAR(wet[r].vapor, expected_vapor, vapor_tol);
            EXPECT_NEAR(dry[r].heat, expected_heat, heat_tol);
            EXPECT_EQ(dry[r].vapor, 0);
            EXPECT_EQ(adiabatic[r].heat, 0);
            EXPECT_NEAR(adiabatic[r].vapor, expected_vapor, vapor_tol);
            EXPECT_NEAR(scalar[r].heat, heat(mean_T), heat_tol);
            EXPECT_NEAR(scalar[r].vapor, vapor(mean_T), vapor_tol);
        }
    }
}

TEST(CloudChamberWallFlux, WetnessSetsVaporPressureForEveryVaporLaw)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;
    FaceWall wall;
    wall.moisture = MoistureMode::WetEquilibrium;
    wall.thermal = {ThermalMode::FixedPhysicalTemperature, Real(300)};
    wall.momentum.z0_m = wall.vapor.z0 = Real(0.01);
    wall.vapor.coefficient = Real(0.1);
    const Real p = Real(90000);
    ScalarWallSample s{1, 0, p, 3, 0, Real(0.125), 2, R_d/Cp_d, Real(0.25)};
    Real saturated;
    erf_qsatw(wall.thermal.temperature_K, p * Real(0.01), saturated);
    EXPECT_EQ(wet_wall_qv(wall, p), saturated);

    for (Real rh : {Real(0), Real(0.4), Real(1)}) {
        wall.wetness = rh;
        const Real e = rh * erf_esatw(wall.thermal.temperature_K) * Real(100); // Pa
        const Real target = RdoRv * e / (p - e);
        EXPECT_NEAR(wet_wall_qv(wall, p), target, scaled_tolerance(target));
        for (auto model : {ScalarModel::ResolvedMolecular, ScalarModel::BulkAero,
                           ScalarModel::NeutralRoughnessLog}) {
            wall.vapor.model = model;
            const Real speed = model == ScalarModel::ResolvedMolecular ? 2*s.alpha_C*s.dx_inv :
                model == ScalarModel::BulkAero ? wall.vapor.coefficient*s.U_t :
                neutral_log_scalar_coefficient(s.wall_distance, wall.momentum.z0_m, wall.vapor.z0)*s.U_t;
            for (Real air : {Real(0), target, Real(0.03)}) {
                s.scalar_air = air;
                const auto flux = evaluate_scalar_flux_in(wall, ScalarChannel::Vapor, s);
                EXPECT_EQ(flux.owned_channels, OwnVapor);
                EXPECT_NEAR(flux.rhoQv_in, s.rho*speed*(target-air), scaled_tolerance(target));
            }
        }
    }
    wall.wetness = 0;
    s.scalar_air = Real(0.01);
    EXPECT_LT(evaluate_scalar_flux_in(wall, ScalarChannel::Vapor, s).rhoQv_in, 0);
    wall.moisture = MoistureMode::DryImpermeable;
    EXPECT_EQ(evaluate_scalar_flux_in(wall, ScalarChannel::Vapor, s).rhoQv_in, 0);
}

TEST(CloudChamberWallFlux, CombinedVaporUsesPrescribedHumidityForBothTerms)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;
    FaceWall wall{};
    wall.moisture = MoistureMode::WetEquilibrium;
    wall.thermal = {ThermalMode::FixedPhysicalTemperature, Real(300)};
    wall.wetness = Real(0.4);
    wall.vapor.model = ScalarModel::BulkAero;
    wall.vapor.provider = CoefficientProvider::Fixed;
    wall.vapor.coefficient = Real(0.125);
    erf_cloud_chamber::WallTransferContract contract;
    contract.vapor_model = "bulk_aero_plus_molecular";
    DiffChoice dc;
    dc.alpha_C = Real(0.0625);
    const Real inv = 4;
    ASSERT_TRUE(erf_cloud_chamber::set_wall_transport_speeds(wall, contract, dc, inv, "zlo").empty());

    // Choose air between the prescribed and saturated wall values: using saturation
    // in either term gives the wrong magnitude or sign. At rest only diffusion remains.
    const Real p = 90000;
    const Real e = wall.wetness * erf_esatw(wall.thermal.temperature_K) * 100; // Pa
    const Real target = RdoRv * e / (p - e);
    Real sat;
    erf_qsatw(wall.thermal.temperature_K, p * Real(0.01), sat);
    ScalarWallSample s{};
    s.rho = 2;
    s.scalar_air = (target + sat) / 2;
    s.p_hse = p;
    s.dx_inv = inv;
    for (Real speed : {Real(0), Real(4)}) {
        s.U_t = speed;
        const Real expected = s.rho * (wall.vapor.coefficient*speed + 2*dc.alpha_C*inv)
                            * (target - s.scalar_air);
        const auto flux = evaluate_scalar_flux_in(wall, ScalarChannel::Vapor, s);
        EXPECT_EQ(flux.owned_channels, OwnVapor);
        EXPECT_LT(flux.rhoQv_in, 0);
        EXPECT_NEAR(flux.rhoQv_in, expected, scaled_tolerance(expected));
    }
}

TEST(CloudChamberWallFlux, MOSTUsesWetnessForStabilityAndVapor)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;
    FaceWall wall;
    wall.moisture = MoistureMode::WetEquilibrium;
    wall.thermal = {ThermalMode::FixedPhysicalTemperature, Real(300)};
    wall.momentum.z0_m = wall.vapor.z0 = Real(0.01);
    wall.vapor.model = ScalarModel::BulkAero;
    wall.vapor.provider = CoefficientProvider::MOST;
    const Real rh = Real(0.4);
    const Real e = rh * erf_esatw(wall.thermal.temperature_K); // hPa
    const Real target = RdoRv * e / (Real(1000) - e);
    ScalarWallSample s{1, target, p_0, 3, 0, 0, 1, R_d/Cp_d, Real(0.5), 300, target, 1};
    auto runtime = [&] {
        return most_wall_coefficients(wall, s.theta_air, s.qv_air, s.p_hse, s.rdOcp,
                                      s.U_t, s.wall_distance, s.gravity_sign);
    };
    EXPECT_LT(runtime().zeta, 0); // Saturated bottom wall is more buoyant than the air.
    wall.wetness = rh;
    EXPECT_NEAR(runtime().zeta, 0, scaled_tolerance(0)); // Equal virtual temperatures.
    EXPECT_NEAR(evaluate_scalar_flux_in(wall, ScalarChannel::Vapor, s).rhoQv_in,
                0, scaled_tolerance(0));
    s.qv_air = s.scalar_air = Real(0.02);
    EXPECT_GT(runtime().zeta, 0);
    EXPECT_NEAR(evaluate_scalar_flux_in(wall, ScalarChannel::Vapor, s).rhoQv_in,
                s.rho * runtime().C_E * s.U_t * (target - s.qv_air), scaled_tolerance(target));
}

// Motivation: dry bulk vapor is algebraically impermeable and must not turn
// on the global wall-rate scan; wet bulk vapor and bulk heat remain active.
TEST(CloudChamberWallConfig, ActivatesTimestepGuardOnlyForActiveBulkChannels)
{
    using erf_cloud_chamber::Config;
    erf_wall_thermodynamics::Boundary walls{};
    walls[0].vapor.model =
        erf_wall_thermodynamics::ScalarModel::BulkAero;
    EXPECT_FALSE(Config::has_bulk_scalar_wall(walls));

    walls[0].moisture =
        erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
    EXPECT_TRUE(Config::has_bulk_scalar_wall(walls));

    walls[0].moisture =
        erf_wall_thermodynamics::MoistureMode::DryImpermeable;
    walls[0].thermal.mode =
        erf_wall_thermodynamics::ThermalMode::FixedPhysicalTemperature;
    walls[0].heat.model =
        erf_wall_thermodynamics::ScalarModel::BulkAero;
    EXPECT_TRUE(Config::has_bulk_scalar_wall(walls));
}

TEST(CloudChamberWallConfig, AdiabaticHeatParser)
{
    using namespace erf_wall_thermodynamics;
    amrex::ParmParse pp("cc_test_heat_selector");
    amrex::Geometry geom; // plate parser does not use sidewall geometry
    FaceWall wall;
    auto parse = [&](bool physical = true) {
        return erf_cloud_chamber::parse_wall_heat(pp, "zlo", physical, wall, geom, nullptr);
    };

    wall.thermal.temperature_K = Real(283); // legacy theta must survive untouched
    EXPECT_TRUE(parse(false).empty());
    EXPECT_EQ(wall.thermal.mode, ThermalMode::LegacyNumerical);
    EXPECT_EQ(wall.thermal.temperature_K, 283);
    EXPECT_NE(parse().find("zlo.temperature"), std::string::npos);
    pp.add("temperature", Real(292));
    EXPECT_TRUE(parse().empty());
    EXPECT_EQ(wall.thermal.mode, ThermalMode::FixedPhysicalTemperature);
    pp.add("heat", "CoNdUcTiVe");
    pp.add("thermal", "FiXeD_TeMpErAtUrE");
    EXPECT_TRUE(parse().empty());
    EXPECT_EQ(wall.thermal.mode, ThermalMode::FixedPhysicalTemperature);
    pp.remove("thermal");
    pp.add("thermal", "AdIaBaTiC");
    EXPECT_NE(parse().find("conflicting"), std::string::npos);
    pp.remove("heat");
    pp.remove("temperature");
    wall.moisture = MoistureMode::DryImpermeable;
    EXPECT_TRUE(parse().empty()); // deployed alias alone
    EXPECT_EQ(wall.thermal.mode, ThermalMode::Adiabatic);
    EXPECT_FALSE(parse(false).empty());
    pp.add("heat", "aDiAbAtIc");
    EXPECT_TRUE(parse().empty()); // matching keys are unambiguous
    pp.remove("thermal");
    EXPECT_TRUE(parse().empty());
    EXPECT_EQ(wall.thermal.mode, ThermalMode::Adiabatic);
    EXPECT_FALSE(parse(false).empty());

    pp.add("temperature", Real(0)); // optional, but never invalid
    EXPECT_NE(parse().find("positive Kelvin"), std::string::npos);
    pp.remove("temperature");
    wall.moisture = MoistureMode::WetEquilibrium;
    EXPECT_NE(parse().find("zlo.temperature"), std::string::npos);
    wall.moisture = MoistureMode::DryImpermeable;
    wall.momentum.model = MomentumModel::BulkAero;
    wall.momentum.provider = CoefficientProvider::Fixed;
    EXPECT_TRUE(parse().empty());
    wall.momentum.provider = CoefficientProvider::MOST;
    EXPECT_TRUE(parse().empty()); // dry adiabatic drag has no prescribed T_wall
    EXPECT_FALSE(wall.heat_is_conductive());
    wall.moisture = MoistureMode::WetEquilibrium;
    EXPECT_NE(parse().find("zlo.temperature"), std::string::npos);
    pp.add("temperature", Real(295));
    EXPECT_TRUE(parse().empty());
    EXPECT_EQ(wall.thermal.temperature_K, 295);
    pp.remove("heat");
    pp.add("heat", "fixed_temperature"); // compatibility spelling is alias-only
    EXPECT_NE(parse().find("zlo.heat"), std::string::npos);
    pp.remove("heat");
    pp.add("thermal", "typo");
    EXPECT_NE(parse().find("zlo.thermal"), std::string::npos);
    pp.remove("thermal");
    pp.remove("temperature");
}

TEST(CloudChamberWallFlux, AdiabaticKeepsActiveMostChannelsAndRates)
{
    using namespace erf_wall_thermodynamics;
    using namespace erf_cloud_chamber_wall_flux;
    Boundary walls{};
    auto& wall = walls[4];
    const auto any_rate_channel = [&] {
        for (const auto& w : walls) { if (w.has_wall_rate_channel()) { return true; } }
        return false;
    };
    wall.thermal = {ThermalMode::Adiabatic, Real(300)};
    wall.moisture = MoistureMode::DryImpermeable;
    wall.momentum.z0_m = Real(.01);
    const Real U = 4, inv = 2, theta = 300, qv = Real(.005), rdOcp = R_d/Cp_d;
    const Real p = 85000;
    const Real distance = Real(.5)/inv;
    wall.heat = {ScalarModel::BulkAero, CoefficientProvider::MOST, 0, 2 * distance};
    ScalarWallSample sample{1, theta, p, U, Real(.1), Real(.1), inv, rdOcp,
                            distance, theta, qv, 1};
    for (auto model : {ScalarModel::ResolvedMolecular, ScalarModel::NeutralRoughnessLog,
                       ScalarModel::BulkAero}) {
        wall.heat.model = model;
        const auto heat = evaluate_scalar_flux_in(wall, ScalarChannel::Heat, sample);
        EXPECT_EQ(heat.owned_channels, OwnHeat);
        EXPECT_EQ(heat.rhoTheta_in, 0);
        EXPECT_FALSE(wall_rate_requires_tangential_speed(wall));
        EXPECT_FALSE(erf_cloud_chamber::Config::has_bulk_scalar_wall(walls));
        EXPECT_FALSE(any_rate_channel());
        EXPECT_FALSE(wall_uses_most(wall));
        EXPECT_EQ(rate_at_origin(walls, {inv, inv, inv}, {U, 0, 0, theta, qv, p, rdOcp}), 0);
    }

    wall.momentum.model = MomentumModel::BulkAero;
    wall.momentum.provider = CoefficientProvider::MOST;
    EXPECT_FALSE(erf_cloud_chamber::Config::has_bulk_scalar_wall(walls));
    EXPECT_TRUE(any_rate_channel());
    // Independent neutral log-law oracle: no MOST coefficient call in the reference.
    const Real k_over_log = KAPPA / std::log(distance / wall.momentum.z0_m);
    const Real cd = k_over_log * k_over_log;
    const MomentumWallSample wind{1, {U, 0, 0}, U, distance, theta, qv, p, rdOcp, 1};
    const Real rate = momentum_infinity_row_sum_factor() * cd * U * inv;
    for (Real temperature : {Real(280), Real(320), ThermalBoundary{}.temperature_K}) {
        SCOPED_TRACE(temperature);
        wall.thermal.temperature_K = temperature;
        const auto dry = most_wall_coefficients(wall, theta, qv, p, rdOcp, U, distance, 1);
        EXPECT_EQ(dry.valid, 1);
        EXPECT_EQ(dry.zeta, 0);
        EXPECT_NEAR(dry.C_D, cd, scaled_tolerance(cd));
        EXPECT_NEAR(dry.u_star, U * k_over_log, scaled_tolerance(U * k_over_log));
        // Even an unused C_H is populated, but not from the unseparated z0_h.
        EXPECT_NEAR(dry.C_H, cd, scaled_tolerance(cd));
        const auto drag = evaluate_momentum_traction(wall, wind);
        EXPECT_EQ(drag.owned_channels, OwnMomentum);
        EXPECT_NEAR(drag.traction_on_fluid[0], -cd * U * U, scaled_tolerance(cd * U * U));
        EXPECT_NEAR(rate_at_origin(walls, {inv, inv, inv}, {U, 0, 0, theta, qv, p, rdOcp}),
                    rate, scaled_tolerance(rate));
    }
    sample.scalar_air = qv;
    EXPECT_EQ(evaluate_scalar_flux_in(wall, ScalarChannel::Vapor, sample).rhoQv_in, 0);

    // Momentum-only MOST on a conductive wall must still see the thermal jump,
    // even when the heat transfer model itself does not use MOST.
    wall.heat = {};
    wall.thermal = {ThermalMode::FixedPhysicalTemperature, Real(280)};
    const auto cold = most_wall_coefficients(wall, theta, qv, p, rdOcp, U, distance, 1);
    wall.thermal.temperature_K = Real(320);
    const auto warm = most_wall_coefficients(wall, theta, qv, p, rdOcp, U, distance, 1);
    EXPECT_GT(cold.zeta, 0);
    EXPECT_LT(cold.C_D, cd);
    EXPECT_LT(warm.zeta, 0);
    EXPECT_GT(warm.C_D, cd);
}

TEST(CloudChamberWallFlux, WetAdiabaticMOSTUsesOnlyHumidityForStability)
{
    using namespace erf_wall_thermodynamics;
    using namespace erf_cloud_chamber_wall_flux;
    Boundary walls{};
    auto& wall = walls[4];
    const auto any_rate_channel = [&] {
        for (const auto& w : walls) { if (w.has_wall_rate_channel()) { return true; } }
        return false;
    };
    wall.thermal = {ThermalMode::Adiabatic, Real(300)};
    wall.moisture = MoistureMode::WetEquilibrium;
    wall.momentum.z0_m = Real(.01);
    wall.vapor = {ScalarModel::BulkAero, CoefficientProvider::MOST, 0, Real(.04)};
    const Real p = 85000, theta = 290, U = 1, inv = 2, rdOcp = R_d/Cp_d;
    const Real distance = Real(.5)/inv;
    wall.heat = {ScalarModel::BulkAero, CoefficientProvider::MOST, 0, 2 * distance};
    EXPECT_TRUE(erf_cloud_chamber::Config::has_bulk_scalar_wall(walls));
    EXPECT_TRUE(any_rate_channel());

    // At 850 hPa, neither the air temperature nor its theta is the 300 K
    // prescribed for saturation. At q_air = q_sat only an adiabatic wall is neutral.
    Real qsat;
    erf_qsatw(Real(300), Real(850), qsat);
    const Real D_m = std::log(distance / wall.momentum.z0_m);
    const Real D_q = std::log(distance / wall.vapor.z0);
    for (Real qv : {qsat - Real(.01), qsat, qsat + Real(.01)}) {
        SCOPED_TRACE(qv);
        const Real Ri = CONST_GRAV * distance * epsv * (qv - qsat) / (U * U);
        const Real zeta = most_stability_coordinate(Ri, D_m, D_q).zeta;
        const auto wet = most_wall_coefficients(wall, theta, qv, p, rdOcp, U, distance, 1);
        EXPECT_EQ(wet.valid, 1);
        EXPECT_NEAR(wet.zeta, zeta, scaled_tolerance(zeta));
        if (qv == qsat) { EXPECT_EQ(wet.zeta, 0); }
        else { EXPECT_GT(wet.zeta * (qv - qsat), 0); }

        similarity_funs sfuns;
        const Real A_m = std::max(D_m - sfuns.calc_psi_m2(zeta), Real(1));
        const Real A_q = std::max(D_q - sfuns.calc_psi_h2(zeta), Real(1));
        const Real ce = KAPPA * KAPPA / (A_m * A_q);
        EXPECT_NEAR(wet.C_E, ce, scaled_tolerance(ce));
        ScalarWallSample sample{1, qv, p, U, 0, 0, inv, rdOcp, distance, theta, qv, 1};
        const auto vapor = evaluate_scalar_flux_in(wall, ScalarChannel::Vapor, sample);
        EXPECT_EQ(vapor.owned_channels, OwnVapor);
        const Real flux = ce * U * (qsat - qv);
        EXPECT_NEAR(vapor.rhoQv_in, flux, scaled_tolerance(ce) * std::abs(U * (qsat - qv)));
        sample.scalar_air = theta;
        const auto heat = evaluate_scalar_flux_in(wall, ScalarChannel::Heat, sample);
        EXPECT_EQ(heat.owned_channels, OwnHeat);
        EXPECT_EQ(heat.rhoTheta_in, 0);
        const Real rate = ce * U * inv;
        EXPECT_NEAR(rate_at_origin(walls, {inv, inv, inv}, {U, 0, 0, theta, qv, p, rdOcp}),
                    rate, scaled_tolerance(rate));
    }

    // Wet equilibrium also drives MOST drag when vapor transfer is molecular.
    // With no MOST scalar, the scalar stability resistance falls back to D_m.
    wall.momentum.model = MomentumModel::BulkAero;
    wall.momentum.provider = CoefficientProvider::MOST;
    wall.vapor = {};
    const Real qv = qsat + Real(.01);
    const Real Ri = CONST_GRAV * distance * epsv * (qv - qsat) / (U * U);
    const Real expected = most_stability_coordinate(Ri, D_m, D_m).zeta;
    const auto drag = most_wall_coefficients(wall, theta, qv, p, rdOcp, U, distance, 1);
    EXPECT_NEAR(drag.zeta, expected, scaled_tolerance(expected));
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

TEST(CloudChamberWallConfig, CombinedRequiresPositiveFiniteDiffusionOnlyForSelectedChannel)
{
    using namespace erf_cloud_chamber;
    using namespace erf_wall_thermodynamics;
    const auto combined = "bulk_aero_plus_molecular";
    for (const Real bad : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(),
                           std::numeric_limits<Real>::quiet_NaN()}) {
        SCOPED_TRACE(bad);
        DiffChoice dc;
        dc.molec_diff_type = MolecDiffType::ConstantAlpha;
        WallTransferContract contract;
        FaceWall wall{};
        wall.thermal.mode = ThermalMode::FixedPhysicalTemperature;
        wall.moisture = MoistureMode::WetEquilibrium;
        const auto setup = [&] {
            // Each setup starts with a newly parsed face, not the previous result.
            auto face = wall;
            return set_wall_transport_speeds(face, contract, dc, 4, "zhi");
        };

        contract.momentum_model = combined;
        dc.dynamic_viscosity = bad;
        EXPECT_NE(setup().find("dynamic_viscosity"), std::string::npos);
        dc.dynamic_viscosity = 1;
        dc.rho0_trans = bad;
        EXPECT_NE(setup().find("rho0_trans"), std::string::npos);
        dc.rho0_trans = 1;
        EXPECT_TRUE(setup().empty()); // unrelated scalar diffusivities default to zero

        contract.momentum_model.clear();
        contract.heat_model = combined;
        dc.alpha_T = bad;
        EXPECT_NE(setup().find("alpha_T"), std::string::npos);
        dc.alpha_T = 1;
        dc.dynamic_viscosity = 0;
        EXPECT_TRUE(setup().empty());

        contract.heat_model.clear();
        contract.vapor_model = combined;
        dc.alpha_C = bad;
        EXPECT_NE(setup().find("alpha_C"), std::string::npos);
        wall.moisture = MoistureMode::DryImpermeable;
        EXPECT_NE(setup().find("alpha_C"), std::string::npos); // selected even though impermeable
        dc.alpha_T = 0;
        dc.alpha_C = 1;
        ASSERT_TRUE(set_wall_transport_speeds(wall, contract, dc, 4, "zhi").empty());
        EXPECT_EQ(wall.vapor.molecular_speed, 0);
        EXPECT_EQ(wall.vapor.resolved_speed, 0);
    }
}

TEST(CloudChamberWallConfig, ResolvedScalarReferenceIsNotAnAdditiveSpeed)
{
    using namespace erf_cloud_chamber;
    using namespace erf_wall_thermodynamics;
    DiffChoice dc;
    dc.alpha_T = Real(.125);
    dc.alpha_C = Real(.25);
    WallTransferContract contract;
    contract.heat_model = "bulk_aero_plus_molecular";
    FaceWall wall{};
    wall.heat.model = ScalarModel::BulkAero;
    wall.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    wall.moisture = MoistureMode::WetEquilibrium;
    ASSERT_TRUE(set_wall_transport_speeds(wall, contract, dc, 4, "zhi").empty());
    EXPECT_EQ(wall.heat.molecular_speed, 1);
    EXPECT_EQ(wall.heat.resolved_speed, 0);
    EXPECT_EQ(wall.vapor.molecular_speed, 0);
    EXPECT_EQ(wall.vapor.resolved_speed, 2);

    FaceWall other{};
    other.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    other.moisture = MoistureMode::WetEquilibrium;
    other.vapor.model = ScalarModel::BulkAero;
    contract.heat_model.clear();
    contract.vapor_model = "bulk_aero";
    ASSERT_TRUE(set_wall_transport_speeds(other, contract, dc, 2, "ylo").empty());
    EXPECT_EQ(other.heat.resolved_speed, Real(.5));
    EXPECT_EQ(other.heat.molecular_speed, 0);
    EXPECT_EQ(other.vapor.resolved_speed, 0);
    EXPECT_EQ(other.vapor.molecular_speed, 0);
}

TEST(CloudChamberWallConfig, DryCombinedVaporParserRequiresBothSelectedCoefficients)
{
    using namespace erf_cloud_chamber;
    using namespace erf_wall_thermodynamics;
    // A private face prefix keeps these parser inputs independent of the run's walls.
    amrex::ParmParse pp("cloud_chamber_dry_combined_parser_test");
    pp.add("vapor_transfer_model", std::string("Bulk_Aero_Plus_Molecular"));
    pp.add("coefficient_source", std::string("Fixed"));
    auto contract = read_wall_transfer_contract(pp, "xhi");
    const auto missing = wall_transfer_contract_error(contract, "xhi");
    EXPECT_NE(missing.find("C_E"), std::string::npos);
    EXPECT_NE(missing.find("vapor_transfer_model = bulk_aero_plus_molecular"), std::string::npos);

    pp.add("C_E", Real(0)); // zero bulk part is allowed, but not zero molecular diffusivity
    contract = read_wall_transfer_contract(pp, "xhi");
    pp.remove("C_E");
    pp.remove("coefficient_source");
    pp.remove("vapor_transfer_model");
    ASSERT_TRUE(wall_transfer_contract_error(contract, "xhi").empty());
    FaceWall wall{};
    wall.vapor.model = ScalarModel::BulkAero;
    wall.moisture = MoistureMode::DryImpermeable;
    DiffChoice dc;
    dc.molec_diff_type = MolecDiffType::ConstantAlpha;
    const auto error = set_wall_transport_speeds(wall, contract, dc, 2, "xhi");
    EXPECT_NE(error.find("Cloud Chamber: xhi.vapor_transfer_model"), std::string::npos);
    EXPECT_NE(error.find("alpha_C"), std::string::npos);
    dc.alpha_C = Real(.125);
    EXPECT_TRUE(set_wall_transport_speeds(wall, contract, dc, 2, "xhi").empty());
    EXPECT_EQ(wall.vapor.molecular_speed, 0);
    EXPECT_EQ(wall.vapor.resolved_speed, 0);
}

TEST(CloudChamberWallConfig, AdditiveChannelRequiresItsOwnFixedCoefficient)
{
    using namespace erf_cloud_chamber;
    for (int channel = 0; channel < 3; ++channel) {
        WallTransferContract contract;
        const char* coefficient;
        if (channel == 0) {
            contract.momentum_model_specified = true;
            contract.momentum_model = "bulk_aero_plus_molecular";
            coefficient = "C_D";
        } else if (channel == 1) {
            contract.heat_model_specified = true;
            contract.heat_model = "bulk_aero_plus_molecular";
            coefficient = "C_H";
        } else {
            contract.vapor_model_specified = true;
            contract.vapor_model = "bulk_aero_plus_molecular";
            coefficient = "C_E";
        }
        SCOPED_TRACE(coefficient);
        const auto error = [&] { return wall_transfer_contract_error(contract, "zlo"); };
        EXPECT_NE(error().find("coefficient_source = fixed"), std::string::npos);
        contract.coefficient_source_specified = true;
        contract.coefficient_source = "most";
        EXPECT_NE(error().find("coefficient_source = fixed"), std::string::npos);
        contract.coefficient_source = "fixed";
        EXPECT_NE(error().find(coefficient), std::string::npos);
        EXPECT_NE(error().find("bulk_aero_plus_molecular"), std::string::npos);
        contract.momentum_coefficient_specified = channel == 0;
        contract.heat_coefficient_specified = channel == 1;
        contract.vapor_coefficient_specified = channel == 2;
        EXPECT_TRUE(error().empty()); // zero bulk coefficient is valid
    }
}

TEST(CloudChamberWallFlux, AdditiveWallUsesMolecularPlusBulkExactlyOnce)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;
    FaceWall bulk{};
    bulk.momentum.model = MomentumModel::BulkAero;
    bulk.heat.model = bulk.vapor.model = ScalarModel::BulkAero;
    bulk.momentum.provider = CoefficientProvider::Fixed;
    bulk.heat.provider = bulk.vapor.provider = CoefficientProvider::Fixed;
    bulk.thermal = {ThermalMode::FixedPhysicalTemperature, Real(300)};
    bulk.moisture = MoistureMode::WetEquilibrium;
    bulk.momentum.C_D = Real(.125);
    bulk.heat.coefficient = Real(.25);
    bulk.vapor.coefficient = Real(.375);
    DiffChoice dc;
    dc.molec_diff_type = MolecDiffType::ConstantAlpha;
    dc.rho0_trans = 2;
    dc.dynamic_viscosity = Real(.03125);
    dc.alpha_T = Real(.0625);
    dc.alpha_C = Real(.125);
    auto wall = bulk;
    erf_cloud_chamber::WallTransferContract contract;
    contract.momentum_model = contract.heat_model = contract.vapor_model = "bulk_aero_plus_molecular";
    ASSERT_TRUE(erf_cloud_chamber::set_wall_transport_speeds(wall, contract, dc, 4, "zhi").empty());
    auto molecular = bulk;
    molecular.heat.model = molecular.vapor.model = ScalarModel::ResolvedMolecular;
    ScalarWallSample sample{Real(1.25), 290, p_0, 0, dc.alpha_T, dc.alpha_C, 4, R_d/Cp_d};
    for (Real speed : {Real(0), Real(3)}) {
        sample.U_t = speed;
        sample.scalar_air = 290;
        const auto heat = [&](const FaceWall& w) {
            return evaluate_scalar_flux_in(w, ScalarChannel::Heat, sample).rhoTheta_in;
        };
        const Real expected_heat = heat(molecular) + heat(bulk);
        EXPECT_GT(expected_heat, 0);
        EXPECT_NEAR(heat(wall), expected_heat, scaled_tolerance(expected_heat));
        sample.scalar_air = Real(.01);
        const auto vapor = [&](const FaceWall& w) {
            return evaluate_scalar_flux_in(w, ScalarChannel::Vapor, sample).rhoQv_in;
        };
        const Real expected_vapor = vapor(molecular) + vapor(bulk);
        EXPECT_GT(expected_vapor, 0);
        EXPECT_NEAR(vapor(wall), expected_vapor, scaled_tolerance(expected_vapor));
        MomentumWallSample m{sample.rho, {speed, 0, 0}, speed};
        EXPECT_DOUBLE_EQ(evaluate_momentum_traction(wall, m).traction_on_fluid[0],
                         -sample.rho * (Real(.125)*speed + Real(.125)) * speed);
        const Real expected_drag_rate = momentum_infinity_row_sum_factor() * Real(.125) * speed * 4 + Real(.5);
        EXPECT_NEAR(momentum_rate_for_face(wall, speed, 4), expected_drag_rate, scaled_tolerance(expected_drag_rate));
        Boundary walls{};
        walls[4] = wall;
        const auto rate = [&] { return rate_at_origin(walls, {1,2,4}, {speed,0,0,0,0,0,0}); };
        EXPECT_EQ(rate(), (Real(.375)*speed + 1)*4); // vapor dominates
        walls[4].moisture = MoistureMode::DryImpermeable;
        EXPECT_EQ(rate(), (Real(.25)*speed + Real(.5))*4); // heat dominates
    }
    wall.moisture = MoistureMode::DryImpermeable;
    EXPECT_EQ(evaluate_scalar_flux_in(wall, ScalarChannel::Vapor, sample).rhoQv_in, 0);
}

TEST(CloudChamberWallFlux, CombinedCornerIncludesResolvedAndPlainBulkOnlyForItsScalar)
{
    using namespace erf_wall_thermodynamics;
    Boundary walls{};
    const GpuArray<Real,3> inv{1,2,4};
    const GpuArray<Real,7> state{0,0,1,0,0,0,0}; // unit speed tangent to x and y
    auto& x = walls[0];
    x.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    x.moisture = MoistureMode::WetEquilibrium;
    x.heat.model = x.vapor.model = ScalarModel::BulkAero;
    x.heat.provider = x.vapor.provider = CoefficientProvider::Fixed;
    x.heat.coefficient = 3;
    x.vapor.coefficient = 2;
    walls[2] = x; // plain bulk on ylo: heat 6, vapor 4
    auto& z = walls[4];
    z.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    z.moisture = MoistureMode::WetEquilibrium;
    DiffChoice dc;
    dc.molec_diff_type = MolecDiffType::ConstantAlpha;
    dc.alpha_T = Real(.25);
    dc.alpha_C = Real(1);
    // Consume setup's actual resolved references: z rates are 8 and 32.
    ASSERT_TRUE(erf_cloud_chamber::set_wall_transport_speeds(z, {}, dc, inv[2], "zlo").empty());
    const auto rate = [&] { return rate_at_origin(walls, inv, state); };

    EXPECT_EQ(rate(), 6); // neither resolved rate independently limits plain bulk
    x.heat.molecular_speed = 2*dc.alpha_T*inv[0];
    EXPECT_EQ(rate(), Real(17.5)); // heat 3+.5+6+8; resolved vapor remains dormant
    z.thermal.mode = ThermalMode::LegacyNumerical;
    EXPECT_EQ(rate(), Real(9.5)); // unowned resolved heat contributes nothing
    z.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    x.heat.molecular_speed = 0;
    x.vapor.molecular_speed = 2*dc.alpha_C*inv[0];
    EXPECT_EQ(rate(), 40); // vapor 2+2+4+32
    z.moisture = MoistureMode::DryImpermeable;
    EXPECT_EQ(rate(), 8); // owned-zero resolved vapor; resolved heat still dormant
    z.moisture = MoistureMode::WetEquilibrium;
    walls[1] = z;
    ASSERT_TRUE(erf_cloud_chamber::set_wall_transport_speeds(walls[1], {}, dc, inv[0], "xhi").empty());
    EXPECT_EQ(rate(), 40); // xhi is not incident
    const amrex::IntVect thin_x(0,1,1);
    EXPECT_EQ(rate_at_origin(walls, inv, state, thin_x), 42);
    std::swap(walls[0], walls[1]); // resolved xlo is visited before combined xhi
    EXPECT_EQ(rate_at_origin(walls, inv, state, thin_x), 42);
}

TEST(CloudChamberWallFlux, MixedMomentumCountsResolvedShearOnlyOnTheSameTangentialRow)
{
    using namespace erf_wall_thermodynamics;
    using erf_cloud_chamber::set_wall_transport_speeds;
    Boundary walls{};
    DiffChoice dc;
    dc.dynamic_viscosity = 6;
    dc.rho0_trans = 2;
    const GpuArray<Real,3> inv{1,2,4};
    ASSERT_TRUE(set_wall_transport_speeds(walls[2], {}, dc, inv[1], "ylo").empty());
    ASSERT_TRUE(set_wall_transport_speeds(walls[4], {}, dc, inv[2], "zlo").empty());
    const auto rate = [&] { return rate_at_origin(walls, inv, {}); };
    EXPECT_EQ(rate(), 0); // resolved alone does not impose a new wall limit
    EXPECT_EQ(walls[4].momentum.molecular_speed, 0); // ordinary stress is not replaced

    auto& x = walls[0];
    x.momentum.model = MomentumModel::BulkAero;
    x.momentum.provider = CoefficientProvider::Fixed;
    x.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    x.heat.model = ScalarModel::BulkAero;
    x.heat.provider = CoefficientProvider::Fixed;
    dc.alpha_T = Real(.5);
    erf_cloud_chamber::WallTransferContract contract;
    contract.heat_model = "bulk_aero_plus_molecular";
    ASSERT_TRUE(set_wall_transport_speeds(x, contract, dc, inv[0], "xlo").empty());
    EXPECT_EQ(rate(), 1); // combined heat and plain drag do not activate resolved shear
    contract.momentum_model = "bulk_aero_plus_molecular";
    ASSERT_TRUE(set_wall_transport_speeds(x, contract, dc, inv[0], "xlo").empty());
    // At rest and C_D=0: xlo contributes 6; ylo and zlo contribute 40 and 160.
    // y-momentum sees xlo+zlo; the larger resolved-only x-momentum row (40+160) is dormant.
    EXPECT_EQ(rate(), 166);
    walls[5] = walls[4]; // both resolved z faces incident only in a one-cell z direction
    EXPECT_EQ(rate(), 166);
    EXPECT_EQ(rate_at_origin(walls, inv, {}, amrex::IntVect(1,1,0)), 326);
    std::swap(walls[0], walls[1]);
    EXPECT_EQ(rate(), 0); // combined xhi is no longer incident
    EXPECT_EQ(rate_at_origin(walls, inv, {}, amrex::IntVect(0,1,1)), 166);
}

TEST(CloudChamberWallFlux, InactiveCombinedFaceCannotTriggerLegacyCornerSum)
{
    using namespace erf_wall_thermodynamics;
    Boundary walls{};
    auto& x = walls[0];
    x.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    x.moisture = MoistureMode::WetEquilibrium;
    x.heat.model = x.vapor.model = ScalarModel::BulkAero;
    x.heat.provider = x.vapor.provider = CoefficientProvider::Fixed;
    x.heat.coefficient = 3;
    x.vapor.coefficient = 4;
    walls[2] = x;
    auto& z = walls[4];
    z = x;
    z.heat.molecular_speed = z.vapor.molecular_speed = 2;
    z.thermal.mode = ThermalMode::LegacyNumerical;
    z.moisture = MoistureMode::DryImpermeable;
    const auto rate = [&] { return rate_at_origin(walls, {1,1,1}, {0,0,1,0,0,0,0}); };
    EXPECT_EQ(rate(), 4); // x/y: legacy max, not heat 6 or vapor 8
    z.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    EXPECT_EQ(rate(), 8); // calm z contributes molecular heat only: 3+3+2
    z.thermal.mode = ThermalMode::LegacyNumerical;
    z.moisture = MoistureMode::WetEquilibrium;
    EXPECT_EQ(rate(), 10); // calm z contributes molecular vapor only: 4+4+2
}

TEST(CloudChamberWallFlux, ResolvedOnlyRateRequiresNoCellState)
{
    using namespace erf_wall_thermodynamics;
    Boundary walls{};
    walls[0].thermal.mode = ThermalMode::FixedPhysicalTemperature;
    walls[0].moisture = MoistureMode::WetEquilibrium;
    walls[0].heat.resolved_speed = 4;
    walls[0].vapor.resolved_speed = 8;
    walls[0].momentum.resolved_row_rate = 12;
    const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(1));
    EXPECT_EQ(erf_cloud_chamber_wall_flux::wall_rate_for_cell(
        0, 0, 0, domain, walls, {1,1,1}, {}), 0);
}

// MOST is shared by the production flux and timestep paths.  This test keeps
// the pointwise state contract visible: wet-wall thermodynamics, orientation,
// and the common coupled scalar roughness all affect the result.
TEST(CloudChamberWallFlux, MOSTUsesWetStateOrientationAndCommonScalarRoughness)
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
    EXPECT_DOUBLE_EQ(zlo.C_H, zlo.C_E);
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

TEST(CloudChamberWallFlux, MOSTUsesActiveScalarRoughnessForStability)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;

    const auto base_wall = [] {
        FaceWall wall;
        wall.thermal.mode = ThermalMode::FixedPhysicalTemperature;
        wall.thermal.temperature_K = Real(290.0);
        wall.moisture = MoistureMode::DryImpermeable;
        wall.momentum.model = MomentumModel::ResolvedNoSlip;
        wall.momentum.z0_m = Real(0.01);
        wall.heat.model = ScalarModel::ResolvedMolecular;
        wall.vapor.model = ScalarModel::ResolvedMolecular;
        return wall;
    };

    FaceWall heat_a = base_wall();
    heat_a.heat.model = ScalarModel::BulkAero;
    heat_a.heat.provider = CoefficientProvider::MOST;
    heat_a.heat.z0 = Real(0.02);
    FaceWall heat_b = heat_a;
    heat_b.heat.z0 = Real(0.04);
    FaceWall heat_inactive_vapor = heat_a;
    heat_inactive_vapor.vapor.z0 = Real(0.04);

    const auto heat_runtime_a = most_wall_coefficients(
        heat_a, Real(285.0), Real(0.005), Real(100000.0), R_d/Cp_d,
        Real(5.0), Real(0.5), 1);
    const auto heat_runtime_b = most_wall_coefficients(
        heat_b, Real(285.0), Real(0.005), Real(100000.0), R_d/Cp_d,
        Real(5.0), Real(0.5), 1);
    const auto heat_runtime_inactive_vapor = most_wall_coefficients(
        heat_inactive_vapor, Real(285.0), Real(0.005), Real(100000.0), R_d/Cp_d,
        Real(5.0), Real(0.5), 1);
    EXPECT_NE(heat_runtime_a.zeta, heat_runtime_b.zeta);
    EXPECT_DOUBLE_EQ(heat_runtime_a.zeta, heat_runtime_inactive_vapor.zeta);

    FaceWall vapor_a = base_wall();
    vapor_a.moisture = MoistureMode::WetEquilibrium;
    vapor_a.vapor.model = ScalarModel::BulkAero;
    vapor_a.vapor.provider = CoefficientProvider::MOST;
    vapor_a.vapor.z0 = Real(0.02);
    FaceWall vapor_b = vapor_a;
    vapor_b.vapor.z0 = Real(0.04);
    FaceWall vapor_inactive_heat = vapor_a;
    vapor_inactive_heat.heat.z0 = Real(0.04);

    const auto vapor_runtime_a = most_wall_coefficients(
        vapor_a, Real(285.0), Real(0.005), Real(100000.0), R_d/Cp_d,
        Real(5.0), Real(0.5), 1);
    const auto vapor_runtime_b = most_wall_coefficients(
        vapor_b, Real(285.0), Real(0.005), Real(100000.0), R_d/Cp_d,
        Real(5.0), Real(0.5), 1);
    const auto vapor_runtime_inactive_heat = most_wall_coefficients(
        vapor_inactive_heat, Real(285.0), Real(0.005), Real(100000.0), R_d/Cp_d,
        Real(5.0), Real(0.5), 1);
    EXPECT_NE(vapor_runtime_a.zeta, vapor_runtime_b.zeta);
    EXPECT_DOUBLE_EQ(vapor_runtime_a.zeta, vapor_runtime_inactive_heat.zeta);
    EXPECT_TRUE(wall_uses_most(vapor_a));

    FaceWall dry_vapor_a = base_wall();
    dry_vapor_a.momentum.model = MomentumModel::BulkAero;
    dry_vapor_a.momentum.provider = CoefficientProvider::MOST;
    dry_vapor_a.momentum.z0_m = Real(0.01);
    dry_vapor_a.vapor.model = ScalarModel::BulkAero;
    dry_vapor_a.vapor.provider = CoefficientProvider::MOST;
    dry_vapor_a.vapor.z0 = Real(0.02);
    FaceWall dry_vapor_b = dry_vapor_a;
    dry_vapor_b.vapor.z0 = Real(0.04);
    const auto dry_runtime_a = most_wall_coefficients(
        dry_vapor_a, Real(285.0), Real(0.005), Real(100000.0), R_d/Cp_d,
        Real(5.0), Real(0.5), 1);
    const auto dry_runtime_b = most_wall_coefficients(
        dry_vapor_b, Real(285.0), Real(0.005), Real(100000.0), R_d/Cp_d,
        Real(5.0), Real(0.5), 1);
    EXPECT_TRUE(wall_uses_most(dry_vapor_a));
    EXPECT_FALSE(wall_uses_most_vapor(dry_vapor_a));
    EXPECT_DOUBLE_EQ(dry_runtime_a.zeta, dry_runtime_b.zeta);
    EXPECT_DOUBLE_EQ(dry_runtime_a.C_D, dry_runtime_b.C_D);

    const auto dry_vapor_flux = evaluate_scalar_flux_in(
        dry_vapor_a, ScalarChannel::Vapor,
        ScalarWallSample{Real(1.0), Real(0.005), Real(100000.0), Real(5.0),
                         Real(0.0), Real(0.0), Real(0.0), R_d/Cp_d,
                         Real(0.5), Real(285.0), Real(0.005), 1});
    EXPECT_EQ(dry_vapor_flux.owned_channels, OwnVapor);
    EXPECT_DOUBLE_EQ(dry_vapor_flux.rhoQv_in, Real(0.0));
}

// Motivation: this valid state needed more than forty relaxed updates; the
// contract is an accurate shared root, not a particular iteration count.
TEST(CloudChamberWallFlux, MOSTStableStateCertifiesSlowlyConvergingRoot)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;

    FaceWall wall;
    wall.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    wall.thermal.temperature_K = Real(300.0) - Real(300.0) / CONST_GRAV;
    wall.moisture = MoistureMode::DryImpermeable;
    wall.momentum.model = MomentumModel::ResolvedNoSlip;
    wall.momentum.z0_m = Real(0.5) / std::exp(Real(2.0));
    wall.heat.model = ScalarModel::BulkAero;
    wall.heat.provider = CoefficientProvider::MOST;
    wall.heat.z0 = Real(0.5) / std::exp(Real(15.0));
    wall.vapor.model = ScalarModel::ResolvedMolecular;

    const auto reference = most_stability_coordinate(Real(0.5), Real(2.0), Real(15.0));
    expect_stable_most_root(reference.zeta, Real(0.5), Real(2.0), Real(15.0),
                            4.237669315020899);

    const auto runtime = most_wall_coefficients(
        wall, Real(300.0), Real(0.0), Real(100000.0), R_d/Cp_d,
        Real(1.0), Real(0.5), 1);
    EXPECT_EQ(runtime.valid, 1);
    EXPECT_NEAR(runtime.zeta, reference.zeta, scaled_tolerance(reference.zeta));
    EXPECT_TRUE(std::isfinite(static_cast<double>(runtime.C_D)));
    EXPECT_TRUE(std::isfinite(static_cast<double>(runtime.C_H)));
}

// Motivation: F-4576-A aborts after one hundred updates even though the
// Jimenez closure has a finite stable root at this admissible roughness.
TEST(CloudChamberWallFlux, MOSTCertifiesRootAfterIterationBudget)
{
    const auto result = erf_cloud_chamber_wall_flux::most_stability_coordinate(
        Real(0.2), Real(0.5), Real(0.5));
    expect_stable_most_root(result.zeta, Real(0.2), Real(0.5), Real(0.5),
                            1.545814710591);
}

// Motivation: F-4576-B has a negative residual near the momentum-resistance
// floor at 0.155, not a root. False convergence inflated drag about 932-fold.
TEST(CloudChamberWallFlux,
     MOSTRejectsFalseConvergenceAndSharesCorrectCoefficients)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;
    const auto result =
        most_stability_coordinate(Real(4.0), Real(0.1), Real(24.4));
    expect_stable_most_root(result.zeta, Real(4.0), Real(0.1), Real(24.4),
                            73.281971961051);

    FaceWall wall;
    wall.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    wall.moisture = MoistureMode::DryImpermeable;
    wall.momentum.model = MomentumModel::BulkAero;
    wall.momentum.provider = CoefficientProvider::MOST;
    wall.momentum.z0_m = Real(0.5) * std::exp(-Real(0.1));
    wall.heat.model = ScalarModel::BulkAero;
    wall.heat.provider = CoefficientProvider::MOST;
    wall.heat.z0 = Real(0.5) * std::exp(-Real(24.4));
    const Real U_t = Real(0.2);
    wall.thermal.temperature_K = Real(300.0) - Real(4.0) * Real(300.0) * U_t *
                                                   U_t /
                                                   (CONST_GRAV * Real(0.5));
    const auto runtime = most_wall_coefficients(
        wall, Real(300.0), Real(0.0), p_0, R_d / Cp_d, U_t, Real(0.5), 1);
    EXPECT_EQ(runtime.valid, 1);
    EXPECT_NEAR(runtime.zeta, Real(73.281971961051), Real(1.0e-3));
    EXPECT_TRUE(std::isfinite(runtime.C_D));
    EXPECT_TRUE(std::isfinite(runtime.C_H));
    // Root-location uncertainty dominates rounding in these coefficients.
    EXPECT_NEAR(runtime.C_D, Real(0.000180425456), Real(2.0e-8));
    const double A_m = 0.41 / std::sqrt(0.0001804254559951541);
    const double A_h =
        24.4 + 5.3 * std::log(73.281971961051 +
                              std::pow(1.0 + std::pow(73.281971961051, 1.1),
                                       1.0 / 1.1));
    EXPECT_NEAR(runtime.C_H, Real(0.41 * 0.41 / (A_m * A_h)), Real(2.0e-8));

    MomentumWallSample momentum;
    momentum.rho = Real(1.2);
    momentum.u_t = {U_t, Real(0.0), Real(0.0)};
    momentum.U_t = U_t;
    momentum.wall_distance = Real(0.5);
    momentum.theta_air = Real(300.0);
    momentum.p_hse = p_0;
    momentum.rdOcp = R_d / Cp_d;
    momentum.gravity_sign = 1;
    const auto traction = evaluate_momentum_traction(wall, momentum);
    EXPECT_NEAR(traction.traction_on_fluid[0],
                -momentum.rho * runtime.C_D * U_t * U_t,
                scaled_tolerance(traction.traction_on_fluid[0]));
    ScalarWallSample scalar;
    scalar.rho = momentum.rho;
    scalar.scalar_air = momentum.theta_air;
    scalar.theta_air = momentum.theta_air;
    scalar.p_hse = p_0;
    scalar.rdOcp = R_d / Cp_d;
    scalar.U_t = U_t;
    scalar.wall_distance = Real(0.5);
    scalar.gravity_sign = 1;
    const auto heat =
        evaluate_scalar_flux_in(wall, ScalarChannel::Heat, scalar);
    EXPECT_NEAR(heat.rhoTheta_in,
                scalar.rho * runtime.C_H * U_t *
                    (wall.thermal.temperature_K - scalar.theta_air),
                scaled_tolerance(heat.rhoTheta_in));

    const Real expected_rate =
        std::max(momentum_infinity_row_sum_factor() * runtime.C_D * U_t,
                 runtime.C_H * U_t);
    Boundary walls{};
    walls[4] = wall;
    const auto rate = [&] {
        return rate_at_origin(walls, {1,1,1}, {U_t,0,0,300,0,p_0,R_d/Cp_d});
    };
    const Real current_rate = rate();
    EXPECT_NEAR(current_rate, expected_rate, scaled_tolerance(expected_rate));
    EXPECT_NEAR(wall_dt_from_max_rate(current_rate),
                wall_dt_from_max_rate(expected_rate),
                scaled_tolerance(wall_dt_from_max_rate(expected_rate)));
    walls[4].thermal.temperature_K = Real(300.0);
    const Real neutral_rate = rate();
    EXPECT_NE(neutral_rate, current_rate);
    EXPECT_LT(wall_dt_from_max_rate(neutral_rate),
              wall_dt_from_max_rate(current_rate));
}

// Motivation: for these inputs the neutral-start iteration reaches the lower
// root. A broad bracket must not replace it with the distant third root.
TEST(CloudChamberWallFlux, MOSTKeepsNeutralConnectedLowerRoot)
{
    const auto result = erf_cloud_chamber_wall_flux::most_stability_coordinate(
        Real(4.0), Real(0.1), Real(25.0));
    expect_stable_most_root(result.zeta, Real(4.0), Real(0.1), Real(25.0),
                            0.152368429532);
}

// Motivation: the neutral-start iteration does not always reach the first
// positive root. Preserve its certified branch rather than silently changing
// a legacy solution to a different root of the same closure.
TEST(CloudChamberWallFlux, MOSTPreservesCertifiedOuterIterationBranch)
{
    const Real Ri = Real(0.39927092596413188);
    const Real D_m = Real(0.001);
    const Real D_s = Real(1.0008842237169981);
    const auto result =
        erf_cloud_chamber_wall_flux::most_stability_coordinate(Ri, D_m, D_s);
    expect_stable_most_root(result.zeta, Ri, D_m, D_s, 6.6318009022167015);
    EXPECT_LT(stable_most_residual(0.168, Ri, D_m, D_s), 0.0);
    EXPECT_GT(stable_most_residual(0.169, Ri, D_m, D_s), 0.0);
}

// Motivation: close to the resistance-floor corner, both ends of a 1e-3
// probe can be negative despite an intervening first-root pair. The fallback
// must not discard that interval just because its endpoint signs agree.
TEST(CloudChamberWallFlux, MOSTFindsNarrowLowerRootPair)
{
    const Real D_s = Real(24.5003);
    const auto result = erf_cloud_chamber_wall_flux::most_stability_coordinate(
        Real(4.0), Real(0.1), D_s);
    EXPECT_NEAR(result.zeta, Real(0.15519547), Real(1.0e-3));
    EXPECT_LT(stable_most_residual(0.154, 4.0, Real(0.1), D_s), 0.0);
    EXPECT_GT(stable_most_residual(0.15519609, 4.0, Real(0.1), D_s), 0.0);
    EXPECT_LT(stable_most_residual(0.157, 4.0, Real(0.1), D_s), 0.0);
}

// Motivation: logarithmic stable corrections allow finite roots above 100;
// importing another MOST solver's zeta bound would change valid chamber fluxes.
TEST(CloudChamberWallFlux, MOSTDoesNotCapStableCoordinateAtHundred)
{
    using namespace erf_cloud_chamber_wall_flux;
    using namespace erf_wall_thermodynamics;
    const auto result =
        most_stability_coordinate(Real(4.0), Real(6.0), Real(6.0));
    expect_stable_most_root(result.zeta, Real(4.0), Real(6.0), Real(6.0),
                            190.376227351357);
    FaceWall wall;
    wall.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    wall.thermal.temperature_K = Real(300.0) - Real(4.0) * Real(300.0) *
                                                   Real(0.2) * Real(0.2) /
                                                   (CONST_GRAV * Real(0.5));
    wall.moisture = MoistureMode::DryImpermeable;
    wall.momentum.model = MomentumModel::BulkAero;
    wall.momentum.provider = CoefficientProvider::MOST;
    wall.momentum.z0_m = Real(0.5) * std::exp(-Real(6.0));
    wall.heat.model = ScalarModel::BulkAero;
    wall.heat.provider = CoefficientProvider::MOST;
    wall.heat.z0 = wall.momentum.z0_m;
    const auto runtime = most_wall_coefficients(
        wall, Real(300.0), Real(0.0), p_0, R_d / Cp_d, Real(0.2), Real(0.5), 1);
    EXPECT_NEAR(runtime.zeta, Real(190.376227351357), Real(1.0e-3));
    EXPECT_TRUE(std::isfinite(runtime.C_D));
    EXPECT_TRUE(std::isfinite(runtime.C_H));
    EXPECT_TRUE(std::isfinite(runtime.C_E));
}

// Motivation: stable root repair must leave neutral and unstable behavior
// unchanged; these values freeze the old neutral-start iteration's outputs.
TEST(CloudChamberWallFlux, MOSTPreservesNeutralAndUnstableControls)
{
    using namespace erf_cloud_chamber_wall_flux;
    EXPECT_EQ(most_stability_coordinate(Real(0.0), Real(6.0), Real(6.0)).zeta,
              Real(0.0));
    const auto unstable =
        most_stability_coordinate(Real(-0.2), Real(6.0), Real(6.0));
    // Float rounding of the unchanged iteration is below 2e-6 here.
    EXPECT_NEAR(unstable.zeta, Real(-1.149188573236494), Real(2.0e-6));
}

// Motivation: host-only evaluation does not exercise device compilation or
// device arithmetic. Run the discriminators in an AMReX kernel and
// check copied-back results against the independent host residual.
TEST(CloudChamberWallFlux, MOSTCertifiedRootsInKernel)
{
    const GpuArray<Real, 5> Ri = {Real(0.2), Real(4.0), Real(4.0), Real(4.0),
                                  Real(0.39927092596413188)};
    const GpuArray<Real, 5> D_m = {Real(0.5), Real(0.1), Real(0.1), Real(6.0),
                                   Real(0.001)};
    const GpuArray<Real, 5> D_s = {Real(0.5), Real(24.4), Real(25.0), Real(6.0),
                                   Real(1.0008842237169981)};
    const GpuArray<double, 5> expected = {1.545814710591, 73.281971961051,
                                          0.152368429532, 190.376227351357,
                                          6.6318009022167015};
    const auto host = most_coordinates_in_kernel(Ri, D_m, D_s);
    for (int n = 0; n < 5; ++n) {
        SCOPED_TRACE(n);
        expect_stable_most_root(host[n], Ri[n], D_m[n], D_s[n], expected[n]);
    }
}

#if GTEST_HAS_DEATH_TEST && !defined(AMREX_USE_MPI)
// Motivation: a numerical failure must stop the kernel even when NDEBUG
// disables AMReX device assertions; it must not return uncertified
// coefficients.
TEST(CloudChamberWallFlux, MOSTRejectsNonfiniteStateInKernel)
{
    const auto saved_style = GTEST_FLAG_GET(death_test_style);
    GTEST_FLAG_SET(death_test_style, "threadsafe");
    const GpuArray<Real, 5> Ri = {Real(4), Real(4), Real(4), Real(4), Real(4)};
    const Real inf = std::numeric_limits<Real>::infinity();
    const GpuArray<Real, 5> D_m = {inf, inf, inf, inf, inf};
    const GpuArray<Real, 5> D_s = {Real(6), Real(6), Real(6), Real(6), Real(6)};
#if defined(AMREX_USE_GPU)
    const char* message = ".*"; // A device trap need not flush its diagnostic.
#else
    const char* message =
        "Cloud Chamber MOST requires finite positive stable resistances";
#endif
    EXPECT_DEATH((void)most_coordinates_in_kernel(Ri, D_m, D_s), message);
    GTEST_FLAG_SET(death_test_style, saved_style);
}
#endif

TEST(CloudChamberWallFlux, MOSTFrozenWallRateUsesCurrentSharedCoefficients)
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
    const Real p_hse = Real(100000.0);
    const Real rdOcp = R_d / Real(900.0);
    const Real U_t = Real(4.0);
    const Real dx_inv = Real(2.0);
    const auto runtime = most_wall_coefficients(
        wall, theta_air, qv_air, p_hse, rdOcp, U_t,
        Real(0.5) / dx_inv, 1);
    const Real expected_momentum = momentum_infinity_row_sum_factor() *
        runtime.C_D * U_t * dx_inv;
    const Real expected_heat = neutral_scalar_rate(runtime.C_H, U_t, dx_inv);
    const Real expected_vapor = neutral_scalar_rate(runtime.C_E, U_t, dx_inv);
    const Real expected_rate = std::max(expected_momentum,
        std::max(expected_heat, expected_vapor));

    Boundary walls{};
    walls[4] = wall;
    const auto rate = [&] {
        return rate_at_origin(walls, {1,1,dx_inv}, {U_t,0,0,theta_air,qv_air,p_hse,rdOcp});
    };
    EXPECT_DOUBLE_EQ(rate(), expected_rate);

    walls[4].thermal.temperature_K = Real(294.0);
    const auto changed_runtime = most_wall_coefficients(
        walls[4], theta_air, qv_air, p_hse, rdOcp, U_t,
        Real(0.5) / dx_inv, 1);
    EXPECT_NE(changed_runtime.zeta, runtime.zeta);
    EXPECT_NE(rate(), expected_rate);
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

TEST(CloudChamberWallFlux, AdiabaticReplacesOrdinaryFluxAndRhsOnEachFace)
{
    using namespace amrex;
    using namespace erf_wall_thermodynamics;
    const Box domain(IntVect(0), IntVect(0));
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    MultiFab state(ba, dm, RhoTheta_comp + 1, 0), prim(ba, dm, PrimTheta_comp + 1, 0);
    MultiFab base(ba, dm, BaseState::num_comps, 0), rhs(ba, dm, RhoTheta_comp + 1, 0);
    state.setVal(1);
    prim.setVal(280);
    base.setVal(p_0);
    MultiFab flux[AMREX_SPACEDIM];
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        BoxArray faces(ba);
        faces.surroundingNodes(d);
        flux[d].define(faces, dm, 1, 0);
    }
    const GpuArray<Real, AMREX_SPACEDIM> inv{2, 4, 8};
    const Real old_flux = 7, other_rhs = 3;
    for (int n = 0; n < 2 * AMREX_SPACEDIM; ++n) {
        SCOPED_TRACE(n);
        const int dir = n / 2;
        const bool high = n % 2 != 0;
        Boundary walls{}; // other faces retain legacy numerical ownership
        walls[n].thermal = {ThermalMode::Adiabatic, Real(300)};
        for (auto& f : flux) { f.setVal(old_flux); }
        // The ordinary face contribution is already included in the RHS.
        // Activate one face at a time so low/high correction errors cannot cancel.
        rhs.setVal(other_rhs + (high ? -old_flux : old_flux) * inv[dir]);
        for (MFIter mfi(state); mfi.isValid(); ++mfi) {
            erf_cloud_chamber_wall_flux::apply(mfi.validbox(), domain, RhoTheta_comp, 0,
                state.const_array(mfi), prim.const_array(mfi), base.const_array(mfi), {}, {}, {},
                rhs.array(mfi), flux[0].array(mfi), flux[1].array(mfi), flux[2].array(mfi),
                inv, walls, Real(.1), Real(0), R_d/Cp_d);
        }
        Gpu::streamSynchronize();
        EXPECT_EQ(sum_region(rhs, domain, RhoTheta_comp), other_rhs);
        for (int face = 0; face < 2 * AMREX_SPACEDIM; ++face) {
            IntVect at(0);
            at[face / 2] = face % 2;
            EXPECT_EQ(value_at(flux[face / 2], at), face == n ? Real(0) : old_flux);
        }
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

    erf_wall_thermodynamics::Boundary walls{};
    for (auto& wall : walls) {
        wall.moisture =
            erf_wall_thermodynamics::MoistureMode::DryImpermeable;
    }
    walls[0].moisture =
        erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
    // Dry bulk metadata must not require velocity in the resolved compatibility API.
    walls[1].vapor.model =
        erf_wall_thermodynamics::ScalarModel::BulkAero;

    walls[0].thermal.mode =
        erf_wall_thermodynamics::ThermalMode::FixedPhysicalTemperature;
    walls[0].thermal.temperature_K = Real(300.0);

    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx_inv =
        {Real(1.0), Real(1.0), Real(1.0)};
    {
        amrex::MFIter mfi(state);
        ASSERT_TRUE(mfi.isValid());
        erf_resolved_wall_flux::apply(
            domain, domain, RhoQ1_comp, 0, state.const_array(mfi), prim.const_array(mfi),
            base.const_array(mfi), rhs.array(mfi), xflux.array(mfi), yflux.array(mfi),
            zflux.array(mfi), dx_inv, walls, Real(0.0), Real(0.1), R_d/Cp_d);
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

    erf_wall_thermodynamics::Boundary walls{};
    for (auto& wall : walls) {
        wall.moisture =
            erf_wall_thermodynamics::MoistureMode::WetEquilibrium;
        wall.thermal.mode =
            erf_wall_thermodynamics::ThermalMode::FixedPhysicalTemperature;
        wall.thermal.temperature_K = Real(300.0);
    }
    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx_inv =
        {Real(1.0), Real(1.0), Real(1.0)};

    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        erf_resolved_wall_flux::apply(
            mfi.validbox(), domain, RhoQ1_comp, 0,
            state.const_array(mfi), prim.const_array(mfi), base.const_array(mfi),
            rhs.array(mfi), xflux.array(mfi), yflux.array(mfi), zflux.array(mfi),
            dx_inv, walls, Real(0.0), Real(0.1), R_d/Cp_d);
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
    using namespace erf_wall_thermodynamics;
    const GpuArray<Real, AMREX_SPACEDIM> low = {Real(1.0), Real(2.0), Real(3.0)};
    const GpuArray<Real, AMREX_SPACEDIM> high = {Real(4.0), Real(5.0), Real(6.0)};
    using namespace erf_cloud_chamber_wall_flux;
    GpuArray<DirectionalMomentumRate, AMREX_SPACEDIM> rates{};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        rates[dir].modeled = low[dir] + high[dir];
    }
    EXPECT_EQ(momentum_row_sum_rate(0, rates), 16);
    EXPECT_EQ(momentum_row_sum_rate(1, rates), 14);
    EXPECT_EQ(momentum_row_sum_rate(2, rates), 12);

    Boundary walls{};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (int side = 0; side < 2; ++side) {
            auto& drag = walls[2*dir+side].momentum;
            drag.model = MomentumModel::BulkAero;
            drag.provider = CoefficientProvider::Fixed;
            drag.molecular_speed = (side ? high : low)[dir]; // C_D = 0
        }
    }
    EXPECT_EQ(rate_at_origin(walls, {1,1,1}, {}, amrex::IntVect(0)), 16);
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
    Boundary rate_walls{};
    rate_walls[4] = rate_wall;
    EXPECT_NEAR(rate_at_origin(rate_walls, {1,1,dx_inv}, {U_t,0,0,0,0,0,0}),
                momentum_infinity_row_sum_factor() * cd * U_t * dx_inv,
                scaled_tolerance(cd));
    EXPECT_DOUBLE_EQ(neutral_scalar_rate(Real(0.3), Real(0.0), dx_inv), Real(0.0));
    FaceWall dry_vapor;
    dry_vapor.moisture = MoistureMode::DryImpermeable;
    dry_vapor.vapor.model = ScalarModel::NeutralRoughnessLog;
    dry_vapor.vapor.z0 = Real(0.02);
    EXPECT_FALSE(wall_rate_requires_tangential_speed(dry_vapor));
    rate_walls[4] = dry_vapor;
    EXPECT_EQ(rate_at_origin(rate_walls, {1,1,dx_inv},
        {std::numeric_limits<Real>::quiet_NaN(),0,0,0,0,0,0}), 0);
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

    // Retaining a heat model on an adiabatic face does not waive its input grammar.
    for (const bool most_heat : {false, true}) {
        auto retained = heat_only();
        if (most_heat) {
            retained.heat_model = "bulk_aero";
            retained.coefficient_source_specified = true;
            retained.coefficient_source = "most";
        }
        for (const Real invalid : {Real(0), std::numeric_limits<Real>::infinity()}) {
            retained.z0_h = invalid;
            EXPECT_NE(erf_cloud_chamber::wall_transfer_contract_error(
                retained, "zlo").find("finite and positive"), std::string::npos);
        }
    }

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

    auto coupled_most = most;
    coupled_most.heat_model_specified = true;
    coupled_most.heat_model = "bulk_aero";
    coupled_most.vapor_model_specified = true;
    coupled_most.vapor_model = "bulk_aero";
    coupled_most.z0_h_specified = true;
    coupled_most.z0_h = Real(0.02);
    coupled_most.z0_q_specified = true;
    coupled_most.z0_q = Real(0.02);
    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(
        coupled_most, "zlo").empty());
    coupled_most.z0_q = Real(0.03);
    // Both lengths are legal inputs; equality depends on physical activity.
    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(
        coupled_most, "zhi").empty());

    auto smooth = momentum_only();
    smooth.momentum_model = "law_of_wall_momentum";
    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(
        smooth, "xlo").empty());

    // Strict parsing still requires z0_q for a dry wall; the runtime dry
    // vapor gate is tested independently by MomentumAndScalarClosures.
    EXPECT_TRUE(erf_cloud_chamber::wall_transfer_contract_error(
        vapor_only(), "xlo").empty());
}

TEST(CloudChamberWallConfig, CoupledMostRoughnessRequiresTwoActiveScalars)
{
    using namespace erf_wall_thermodynamics;
    FaceWall wall;
    wall.thermal.mode = ThermalMode::FixedPhysicalTemperature;
    wall.moisture = MoistureMode::WetEquilibrium;
    wall.heat = {ScalarModel::BulkAero, CoefficientProvider::MOST, 0, Real(0.02)};
    wall.vapor = wall.heat;
    auto error = [&] { return erf_cloud_chamber::wall_coupled_most_roughness_error(wall, "zhi"); };
    EXPECT_TRUE(error().empty());
    wall.vapor.z0 = Real(0.03);
    const auto mismatch = error();
    EXPECT_NE(mismatch.find("zhi"), std::string::npos);
    EXPECT_NE(mismatch.find("z0_h == z0_q"), std::string::npos);
    EXPECT_NE(mismatch.find("z0_h=0.020000"), std::string::npos);
    EXPECT_NE(mismatch.find("z0_q=0.030000"), std::string::npos);

    wall.moisture = MoistureMode::DryImpermeable;
    EXPECT_TRUE(error().empty());
    wall.moisture = MoistureMode::WetEquilibrium;
    wall.thermal.mode = ThermalMode::Adiabatic;
    EXPECT_TRUE(error().empty());
}

TEST(CloudChamberWallConfig, AdiabaticSkipsOnlyUnusedRoughnessGeometry)
{
    using namespace erf_wall_thermodynamics;
    const GpuArray<Real, AMREX_SPACEDIM> dx{1, 1, 1};
    const Real z_ref = Real(0.5);
    for (auto model : {ScalarModel::NeutralRoughnessLog, ScalarModel::BulkAero}) {
        SCOPED_TRACE(model == ScalarModel::BulkAero ? "MOST" : "neutral");
        Boundary walls{};
        auto& wall = walls[4]; // zlo
        wall.thermal.mode = ThermalMode::Adiabatic;
        wall.moisture = MoistureMode::DryImpermeable;
        wall.momentum.z0_m = Real(0.1);
        const auto provider = model == ScalarModel::BulkAero ? CoefficientProvider::MOST
                                                             : CoefficientProvider::None;
        wall.heat = {model, provider, 0, z_ref};
        auto error = [&] { return erf_cloud_chamber::wall_roughness_geometry_error(walls, dx); };
        for (auto retained : {z_ref, 2 * z_ref}) {
            wall.heat.z0 = retained;
            EXPECT_TRUE(error().empty());
            wall.thermal.mode = ThermalMode::FixedPhysicalTemperature;
            EXPECT_NE(error().find("z0_h"), std::string::npos);
            wall.thermal.mode = ThermalMode::Adiabatic;
        }

        // The shared z0_m is unused until a momentum or wet vapor closure needs it.
        wall.momentum.z0_m = z_ref;
        EXPECT_TRUE(error().empty());
        wall.momentum.model = model == ScalarModel::BulkAero ? MomentumModel::BulkAero
                                                             : MomentumModel::NeutralRoughnessLog;
        wall.momentum.provider = provider;
        EXPECT_NE(error().find("z0_m"), std::string::npos);
        wall.momentum.model = MomentumModel::ResolvedNoSlip;
        wall.momentum.provider = CoefficientProvider::None;
        wall.moisture = MoistureMode::WetEquilibrium;
        wall.vapor = {model, provider, 0, z_ref};
        EXPECT_NE(error().find("z0_m"), std::string::npos);
        wall.momentum.z0_m = Real(0.1);
        EXPECT_NE(error().find("z0_q"), std::string::npos);
        wall.moisture = MoistureMode::DryImpermeable;
        EXPECT_NE(error().find("z0_q"), std::string::npos); // retained vapor grammar is unchanged
        wall.moisture = MoistureMode::WetEquilibrium;
        wall.vapor.z0 = Real(0.2);
        EXPECT_TRUE(error().empty()); // z0_h is still too large and differs from z0_q
    }
}

TEST(CloudChamberNeutralLog, GeometryValidationUsesEveryActiveFace)
{
    using namespace erf_wall_thermodynamics;

    Boundary walls{};
    walls[0].momentum.model = MomentumModel::NeutralRoughnessLog;
    walls[0].momentum.z0_m = Real(0.49);
    walls[2].thermal.mode = ThermalMode::FixedPhysicalTemperature;
    walls[2].heat.model = ScalarModel::NeutralRoughnessLog;
    walls[2].momentum.z0_m = Real(0.99);
    walls[2].heat.z0 = Real(0.99);
    walls[4].moisture = MoistureMode::WetEquilibrium;
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
    scalar_only[0].thermal.mode = ThermalMode::FixedPhysicalTemperature;
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
    scalar_only[0].moisture = MoistureMode::WetEquilibrium;
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

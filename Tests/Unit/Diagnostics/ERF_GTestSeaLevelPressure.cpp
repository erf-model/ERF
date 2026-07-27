#include <algorithm>
#include <cmath>
#include <limits>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <gtest/gtest.h>

#include "Diagnostics/ERF_SeaLevelPressure.H"
#include "ERF_EOS.H"
#include "ERF_IndexDefines.H"
#include "ERF_TerrainMetrics.H"

namespace
{

using namespace sea_level_pressure_diagnostics;

amrex::Real tolerance (amrex::Real expected)
{
    return amrex::Real(128.0) * std::numeric_limits<amrex::Real>::epsilon() *
           std::max(std::abs(expected), amrex::Real(1.0));
}

amrex::Real independent_method_reference (amrex::Real ground_pressure,
                                           amrex::Real lowest_temperature,
                                           amrex::Real specific_humidity,
                                           amrex::Real lowest_height,
                                           amrex::Real surface_height)
{
    const amrex::Real lowest_virtual_temperature =
        lowest_temperature * (amrex::Real(1.0) + epsv * specific_humidity);
    const amrex::Real surface_virtual_temperature =
        lowest_virtual_temperature + amrex::Real(0.0065) *
        (lowest_height - surface_height);
    const amrex::Real uncorrected_sea_level_virtual_temperature =
        lowest_virtual_temperature + amrex::Real(0.0065) * lowest_height;
    amrex::Real corrected_sea_level_virtual_temperature =
        uncorrected_sea_level_virtual_temperature;
    if (uncorrected_sea_level_virtual_temperature > amrex::Real(290.66) &&
        surface_virtual_temperature <= amrex::Real(290.66)) {
        corrected_sea_level_virtual_temperature = amrex::Real(290.66);
    } else if (uncorrected_sea_level_virtual_temperature > amrex::Real(290.66) &&
               surface_virtual_temperature > amrex::Real(290.66)) {
        const amrex::Real delta = surface_virtual_temperature - amrex::Real(290.66);
        corrected_sea_level_virtual_temperature = amrex::Real(290.66) -
                                                   amrex::Real(0.005) * delta * delta;
    }
    const amrex::Real mean_tau = (R_d / CONST_GRAV) *
        (surface_virtual_temperature + corrected_sea_level_virtual_temperature) /
        amrex::Real(2.0);
    return std::abs(surface_height) <= amrex::Real(1.0)
        ? ground_pressure
        : ground_pressure * std::exp(surface_height / mean_tau);
}

struct ColumnState
{
    amrex::Real pressure;
    amrex::Real temperature;
    amrex::Real vapor_mixing_ratio;
    amrex::Real lowest_height;
    amrex::Real surface_height;
};

amrex::Real complete_reference (const ColumnState& state)
{
    const amrex::Real q = state.vapor_mixing_ratio /
                          (amrex::Real(1.0) + state.vapor_mixing_ratio);
    const amrex::Real rho = getRhogivenTandPress(state.temperature,
                                                  state.pressure,
                                                  state.vapor_mixing_ratio);
    const amrex::Real theta = getThgivenTandP(state.temperature,
                                               state.pressure,
                                               RdoCp);
    const amrex::Real rho_theta = rho * theta;
    const amrex::Real diagnosed_pressure = getPgivenRTh(rho_theta,
                                                         state.vapor_mixing_ratio);
    const amrex::Real diagnosed_temperature = getTgivenRandRTh(
        rho, rho_theta, state.vapor_mixing_ratio);
    const amrex::Real lowest_virtual_temperature =
        virtual_temperature(diagnosed_temperature, q);
    const amrex::Real surface_virtual_temperature =
        lowest_virtual_temperature + amrex::Real(0.0065) *
        (state.lowest_height - state.surface_height);
    const amrex::Real mean_tau = (R_d / CONST_GRAV) *
        (lowest_virtual_temperature + surface_virtual_temperature) /
        amrex::Real(2.0);
    const amrex::Real reconstructed_surface_pressure = diagnosed_pressure *
        std::exp((state.lowest_height - state.surface_height) / mean_tau);
    return independent_method_reference(
        reconstructed_surface_pressure, diagnosed_temperature, q,
        state.lowest_height, state.surface_height);
}

void launch_scalar_calculation (amrex::Gpu::DeviceVector<amrex::Real>& output)
{
    auto* output_ptr = output.data();
    amrex::ParallelFor(2, [=] AMREX_GPU_DEVICE (int index) noexcept {
        output_ptr[index] = reduce_surface_pressure_to_sea_level(
            amrex::Real(90000.0), amrex::Real(280.0), amrex::Real(0.0),
            amrex::Real(1050.0), amrex::Real(1000.0));
    });
}

void initialize_dry_conserved_state (amrex::MultiFab& conserved)
{
    for (amrex::MFIter mfi(conserved); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto c = conserved.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            c(i, j, k, Rho_comp) = amrex::Real(1.0);
            c(i, j, k, RhoTheta_comp) = getRhoThetagivenP(amrex::Real(90000.0), zero);
        });
    }
}

void initialize_physical_heights (amrex::MultiFab& heights)
{
    for (amrex::MFIter mfi(heights); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto z = heights.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            z(i, j, k) = amrex::Real(100.0 * i) + amrex::Real(100.0 * (k - 4));
        });
    }
}

void initialize_moist_conserved_state (amrex::MultiFab& conserved)
{
    for (amrex::MFIter mfi(conserved); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto c = conserved.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            c(i, j, k, RhoTheta_comp) = getRhoThetagivenP(
                amrex::Real(90000.0), amrex::Real(0.012));
            c(i, j, k, RhoQ1_comp) = amrex::Real(0.012);
        });
    }
}


// Motivation: Moderate columns must use the unmodified hydrostatic reduction;
// this prevents applying the Shuell limiter outside its stated regime.
TEST(SeaLevelPressureMethod, NoShuellCorrection)
{
    const amrex::Real actual = reduce_surface_pressure_to_sea_level(
        amrex::Real(90000.0), amrex::Real(280.0), amrex::Real(0.0),
        amrex::Real(1050.0), amrex::Real(1000.0));
    const amrex::Real expected = amrex::Real(101529.189147592);
    EXPECT_NEAR(actual, expected, tolerance(expected));
    EXPECT_NEAR(actual, independent_method_reference(
        amrex::Real(90000.0), amrex::Real(280.0), amrex::Real(0.0),
        amrex::Real(1050.0), amrex::Real(1000.0)), tolerance(expected));
}

// Motivation: When only the sea-level extrapolation exceeds the critical
// temperature, UPP caps that endpoint; this prevents selecting the quadratic
// branch for the wrong temperature state.
TEST(SeaLevelPressureMethod, SeaLevelOnlyUsesCapBranch)
{
    const amrex::Real actual = reduce_surface_pressure_to_sea_level(
        amrex::Real(90000.0), amrex::Real(288.0), amrex::Real(0.0),
        amrex::Real(500.0), amrex::Real(300.0));
    EXPECT_NEAR(actual, amrex::Real(93239.546961953), tolerance(actual));
    EXPECT_EQ(apply_shuell_correction(amrex::Real(289.30), amrex::Real(291.25)),
              amrex::Real(290.66));
}

// Motivation: Warm elevated columns require the quadratic Shuell correction;
// this prevents silently reducing the method to a simple temperature cap.
TEST(SeaLevelPressureMethod, BothEndpointsUseQuadraticBranch)
{
    const amrex::Real actual = reduce_surface_pressure_to_sea_level(
        amrex::Real(90000.0), amrex::Real(292.0), amrex::Real(0.0),
        amrex::Real(500.0), amrex::Real(300.0));
    EXPECT_NEAR(actual, amrex::Real(93217.160215006), tolerance(actual));
    EXPECT_NEAR(apply_shuell_correction(amrex::Real(293.30), amrex::Real(295.25)),
                amrex::Real(290.625152), tolerance(amrex::Real(290.625152)));
}

// Motivation: The source uses strict sea-level and non-strict surface
// comparisons at 290.66 K; this protects exact-threshold behavior.
TEST(SeaLevelPressureMethod, ExactAndAdjacentShuellBoundaries)
{
    const amrex::Real critical = shuell_critical_temperature;
    const amrex::Real above = std::nextafter(critical,
                                              std::numeric_limits<amrex::Real>::infinity());
    const amrex::Real below = std::nextafter(critical,
                                              -std::numeric_limits<amrex::Real>::infinity());
    EXPECT_EQ(apply_shuell_correction(critical, critical), critical);
    EXPECT_EQ(apply_shuell_correction(critical, above), critical);
    EXPECT_EQ(apply_shuell_correction(above, critical), critical);
    EXPECT_EQ(apply_shuell_correction(below, above), critical);
    EXPECT_EQ(apply_shuell_correction(above, below), below);
}

// Motivation: UPP skips the final reduction within one meter of z=0; this
// preserves the exact deadband and its boundary inequalities.
TEST(SeaLevelPressureMethod, SurfaceHeightDeadband)
{
    constexpr amrex::Real pressure = amrex::Real(91000.0);
    for (const amrex::Real height : {amrex::Real(0.0), amrex::Real(1.0),
                                     amrex::Real(-1.0)}) {
        EXPECT_EQ(reduce_surface_pressure_to_sea_level(
                      pressure, amrex::Real(280.0), amrex::Real(0.0),
                      amrex::Real(100.0), height), pressure);
    }
    EXPECT_NE(reduce_surface_pressure_to_sea_level(
                  pressure, amrex::Real(280.0), amrex::Real(0.0),
                  amrex::Real(100.0), amrex::Real(1.01)), pressure);
    EXPECT_NE(reduce_surface_pressure_to_sea_level(
                  pressure, amrex::Real(280.0), amrex::Real(0.0),
                  amrex::Real(100.0), amrex::Real(-1.01)), pressure);
}

// Motivation: Below-datum terrain reverses the hydrostatic sign; this catches
// an absolute-value or sign error in the final exponential height term.
TEST(SeaLevelPressureMethod, NegativeSurfaceElevationReducesPressure)
{
    const amrex::Real ground_pressure = amrex::Real(90000.0);
    const amrex::Real actual = reduce_surface_pressure_to_sea_level(
        ground_pressure, amrex::Real(280.0), amrex::Real(0.0),
        amrex::Real(100.0), amrex::Real(-100.0));
    EXPECT_LT(actual, ground_pressure);
}

// Motivation: ERF stores vapor mixing ratio per unit dry air, while the
// reduction uses specific humidity; this prevents a convention mismatch.
TEST(SeaLevelPressureAdapter, ConvertsDryMixingRatioToSpecificHumidity)
{
    const amrex::Real expected = amrex::Real(0.012) / amrex::Real(1.012);
    EXPECT_NEAR(specific_humidity_from_dry_mixing_ratio(amrex::Real(0.012)),
                expected, tolerance(expected));
}

// Motivation: ERF's native gas constants must replace UPP's rounded 0.608;
// this protects the shared EOS/diagnostic thermodynamic convention.
TEST(SeaLevelPressureAdapter, UsesERFVirtualTemperatureCoefficient)
{
    const amrex::Real temperature = amrex::Real(280.0);
    const amrex::Real humidity = amrex::Real(0.25);
    const amrex::Real expected = temperature * (amrex::Real(1.0) + epsv * humidity);
    EXPECT_NEAR(virtual_temperature(temperature, humidity), expected, tolerance(expected));
    EXPECT_NE(virtual_temperature(temperature, humidity),
              temperature * (amrex::Real(1.0) + amrex::Real(0.608) * humidity));
}

// Motivation: Reconstructing ground pressure must be the identity when the
// lowest cell center lies at the ground interface; this protects exponent sign.
TEST(SeaLevelPressureAdapter, SurfacePressureReconstructionIdentity)
{
    const amrex::Real pressure = amrex::Real(90000.0);
    EXPECT_EQ(reconstruct_surface_pressure(
                  pressure, amrex::Real(280.0), amrex::Real(280.0),
                  amrex::Real(100.0), amrex::Real(100.0)), pressure);
}

// Motivation: Pressure must increase when reducing down from an elevated
// lowest model level to the ground; this catches a reversed height difference.
TEST(SeaLevelPressureAdapter, GroundPressureIncreasesDownward)
{
    const amrex::Real lowest_pressure = amrex::Real(90000.0);
    const amrex::Real q = amrex::Real(0.0);
    const amrex::Real temperature = amrex::Real(280.0);
    const amrex::Real tv = virtual_temperature(temperature, q);
    const amrex::Real tv_surface = tv + lapse_rate * amrex::Real(50.0);
    const amrex::Real ground_pressure = reconstruct_surface_pressure(
        lowest_pressure, tv, tv_surface, amrex::Real(1050.0), amrex::Real(1000.0));
    EXPECT_GT(ground_pressure, lowest_pressure);
}

// Motivation: Moist conserved-state wiring must include EOS diagnosis,
// humidity conversion, and ground-pressure reconstruction; this prevents a
// dry-only implementation or direct use of lowest-level pressure.
TEST(SeaLevelPressureAdapter, CompleteDryAndMoistColumnReferences)
{
    const ColumnState dry{amrex::Real(90000.0), amrex::Real(280.0), amrex::Real(0.0),
                          amrex::Real(1050.0), amrex::Real(1000.0)};
    const ColumnState moist{amrex::Real(90000.0), amrex::Real(280.0), amrex::Real(0.012),
                            amrex::Real(1050.0), amrex::Real(1000.0)};

    const auto diagnose = [](const ColumnState& state) {
        const amrex::Real rho = getRhogivenTandPress(state.temperature, state.pressure,
                                                       state.vapor_mixing_ratio);
        const amrex::Real theta = getThgivenTandP(state.temperature, state.pressure, RdoCp);
        return compute_from_erf_lowest_level_state(
            rho, rho * theta, state.vapor_mixing_ratio,
            state.lowest_height, state.surface_height);
    };
    EXPECT_NEAR(diagnose(dry), amrex::Real(102150.434786223), tolerance(diagnose(dry)));
    EXPECT_NEAR(diagnose(moist), amrex::Real(102058.982411126), tolerance(diagnose(moist)));
    EXPECT_NEAR(diagnose(dry), complete_reference(dry), tolerance(diagnose(dry)));
    EXPECT_NE(diagnose(moist), diagnose(dry));
}

// Motivation: The diagnostic is evaluated in an AMReX device kernel; this
// protects host/device portability of the scalar method and EOS adapter.
TEST(SeaLevelPressureKernel, ScalarCalculationIsDeviceCallable)
{
    amrex::Gpu::DeviceVector<amrex::Real> output(2, amrex::Real(0.0));
    launch_scalar_calculation(output);
    amrex::Gpu::streamSynchronize();
    amrex::Gpu::HostVector<amrex::Real> host_output(output.size());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, output.begin(), output.end(), host_output.begin());
    for (const auto value : host_output) {
        EXPECT_NEAR(value, amrex::Real(101529.189147592),
                    amrex::Real(256.0) * tolerance(amrex::Real(101529.189147592)));
    }
}

// Motivation: The 2D destination slab is k=0 while the source domain may have
// a nonzero klo; this catches source-index and terrain-column assembly errors.
TEST(SeaLevelPressureMultiFab, ReadsLowestSourceLevelAndLocalTerrain)
{
    const amrex::Box source_domain(amrex::IntVect(0, 0, 4), amrex::IntVect(1, 1, 6));
    amrex::BoxArray source_ba(source_domain);
    amrex::DistributionMapping dm(source_ba);
    amrex::MultiFab conserved(source_ba, dm, RhoQ1_comp + 1, 0);
    const amrex::BoxArray nodal_ba = amrex::convert(source_ba, amrex::IntVect(1, 1, 1));
    amrex::MultiFab heights(nodal_ba, dm, 1, 0);
    const amrex::Box slab(amrex::IntVect(0, 0, 0), amrex::IntVect(1, 1, 0));
    amrex::MultiFab destination(amrex::BoxArray(slab), dm, 1, 0);

    conserved.setVal(amrex::Real(0.0));
    heights.setVal(amrex::Real(0.0));
    destination.setVal(amrex::Real(-1.0));
    initialize_dry_conserved_state(conserved);
    initialize_physical_heights(heights);
    amrex::Gpu::streamSynchronize();

    fill_sea_level_pressure(destination, 0, conserved, heights, -1, 4);
    amrex::Gpu::streamSynchronize();

    amrex::MultiFab host_destination(
        destination.boxArray(), destination.DistributionMap(), destination.nComp(),
        destination.nGrowVect(), amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    amrex::MultiFab::Copy(host_destination, destination, 0, 0, 1, 0);
    amrex::Gpu::streamSynchronize();
    for (amrex::MFIter mfi(host_destination); mfi.isValid(); ++mfi) {
        const auto d = host_destination.const_array(mfi);
        const ColumnState cell_zero{amrex::Real(90000.0), amrex::Real(90000.0 / R_d),
                                    amrex::Real(0.0), amrex::Real(100.0), amrex::Real(50.0)};
        const ColumnState cell_one{amrex::Real(90000.0), amrex::Real(90000.0 / R_d),
                                   amrex::Real(0.0), amrex::Real(200.0), amrex::Real(150.0)};
        EXPECT_NEAR(d(0, 0, 0), complete_reference(cell_zero), tolerance(d(0, 0, 0)));
        EXPECT_NEAR(d(1, 0, 0), complete_reference(cell_one), tolerance(d(1, 0, 0)));
        EXPECT_NE(d(0, 0, 0), d(1, 0, 0));
    }

    // Rebuild the manufactured state with a valid vapor component while
    // retaining the same 90-kPa lowest-level pressure for the moist pass.
    initialize_moist_conserved_state(conserved);
    fill_sea_level_pressure(destination, 0, conserved, heights, RhoQ1_comp, 4);
    amrex::Gpu::streamSynchronize();
    amrex::MultiFab::Copy(host_destination, destination, 0, 0, 1, 0);
    amrex::Gpu::streamSynchronize();
    const amrex::Real moist_temperature = amrex::Real(90000.0) /
        (R_d * (amrex::Real(1.0) + RvoRd * amrex::Real(0.012)));
    for (amrex::MFIter mfi(host_destination); mfi.isValid(); ++mfi) {
        const auto d = host_destination.const_array(mfi);
        const ColumnState cell_zero{amrex::Real(90000.0), moist_temperature,
                                    amrex::Real(0.012), amrex::Real(100.0), amrex::Real(50.0)};
        const ColumnState cell_one{amrex::Real(90000.0), moist_temperature,
                                   amrex::Real(0.012), amrex::Real(200.0), amrex::Real(150.0)};
        EXPECT_NEAR(d(0, 0, 0), complete_reference(cell_zero), tolerance(d(0, 0, 0)));
        EXPECT_NEAR(d(1, 0, 0), complete_reference(cell_one), tolerance(d(1, 0, 0)));
    }
}

} // namespace

#include <algorithm>
#include <array>
#include <string>
#include <vector>

#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include "ERF_GTestMicrophysicsCommon.H"

using namespace microphysics_test;

namespace {

template <typename HostVectorT>
void copy_host_to_device (const HostVectorT& host,
                          amrex::Gpu::DeviceVector<amrex::Real>& device)
{
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, host.begin(), host.end(), device.begin());
}

amrex::Gpu::HostVector<amrex::Real>
make_host_vector (const std::vector<amrex::Real>& values)
{
    amrex::Gpu::HostVector<amrex::Real> host(values.size());
    std::copy(values.begin(), values.end(), host.begin());
    return host;
}

void append_dense_window (std::vector<amrex::Real>& values,
                          const amrex::Real center)
{
    for (int index = -10; index <= 10; ++index) {
        values.push_back(center + amrex::Real(index) * amrex::Real(1.0e-3));
    }
}

std::vector<amrex::Real> make_water_temperatures ()
{
    std::vector<amrex::Real> temperatures;
    temperatures.reserve(208);

    for (int value = 180; value <= 220; ++value) {
        temperatures.push_back(amrex::Real(value));
    }
    for (int value = 221; value <= 310; ++value) {
        temperatures.push_back(amrex::Real(value));
    }
    for (int value = 311; value <= 345; ++value) {
        temperatures.push_back(amrex::Real(value));
    }
    append_dense_window(temperatures, kWaterLowerSwitchK);
    append_dense_window(temperatures, kWaterUpperSwitchK);

    return temperatures;
}

std::vector<amrex::Real> make_pressures ()
{
    std::vector<amrex::Real> pressures;
    pressures.reserve(20);

    for (int value = 100; value <= 1050; value += 50) {
        pressures.push_back(amrex::Real(value));
    }

    return pressures;
}

void add_sample_failure (const char* metric,
                         const int sample_index,
                         const amrex::Real temperature,
                         const amrex::Real pressure,
                         const amrex::Real host_value,
                         const amrex::Real device_value,
                         const amrex::Real factor)
{
    ADD_FAILURE() << metric
                  << " sample=" << sample_index
                  << " T=" << temperature
                  << " p=" << pressure
                  << " host=" << host_value
                  << " device=" << device_value
                  << " normalized_error=" << normalized_error(device_value, host_value, factor);
}

void expect_sample_match (const char* metric,
                          const int sample_index,
                          const amrex::Real temperature,
                          const amrex::Real pressure,
                          const amrex::Real host_value,
                          const amrex::Real device_value,
                          const amrex::Real factor)
{
    const amrex::Real tol = scaled_tol(device_value, host_value, factor);
    if (std::abs(device_value - host_value) > tol) {
        add_sample_failure(metric, sample_index, temperature, pressure, host_value, device_value, factor);
    }
}

void launch_water_sweep (amrex::Gpu::DeviceVector<amrex::Real>& temperatures,
                         amrex::Gpu::DeviceVector<amrex::Real>& pressures,
                         amrex::Gpu::DeviceVector<amrex::Real>& esat,
                         amrex::Gpu::DeviceVector<amrex::Real>& dtesat,
                         amrex::Gpu::DeviceVector<amrex::Real>& qsat,
                         amrex::Gpu::DeviceVector<amrex::Real>& dtqsat)
{
    auto* temperature_ptr = temperatures.data();
    auto* pressure_ptr = pressures.data();
    auto* esat_ptr = esat.data();
    auto* dtesat_ptr = dtesat.data();
    auto* qsat_ptr = qsat.data();
    auto* dtqsat_ptr = dtqsat.data();

    amrex::ParallelFor(temperatures.size(), [=] AMREX_GPU_DEVICE (int index) noexcept {
        const amrex::Real temperature = temperature_ptr[index];
        const amrex::Real pressure = pressure_ptr[index];
        amrex::Real qsat_value;
        amrex::Real dtqsat_value;

        esat_ptr[index] = erf_esatw(temperature);
        dtesat_ptr[index] = erf_dtesatw(temperature);
        erf_qsatw(temperature, pressure, qsat_value);
        erf_dtqsatw(temperature, pressure, dtqsat_value);
        qsat_ptr[index] = qsat_value;
        dtqsat_ptr[index] = dtqsat_value;
    });
}

void launch_ice_value_sweep (amrex::Gpu::DeviceVector<amrex::Real>& temperatures,
                             amrex::Gpu::DeviceVector<amrex::Real>& esat)
{
    auto* temperature_ptr = temperatures.data();
    auto* esat_ptr = esat.data();

    amrex::ParallelFor(temperatures.size(), [=] AMREX_GPU_DEVICE (int index) noexcept {
        esat_ptr[index] = erf_esati(temperature_ptr[index]);
    });
}

void launch_ice_derivative_sweep (amrex::Gpu::DeviceVector<amrex::Real>& temperatures,
                                  amrex::Gpu::DeviceVector<amrex::Real>& pressures,
                                  amrex::Gpu::DeviceVector<amrex::Real>& dtesat,
                                  amrex::Gpu::DeviceVector<amrex::Real>& qsat,
                                  amrex::Gpu::DeviceVector<amrex::Real>& dtqsat)
{
    auto* temperature_ptr = temperatures.data();
    auto* pressure_ptr = pressures.data();
    auto* dtesat_ptr = dtesat.data();
    auto* qsat_ptr = qsat.data();
    auto* dtqsat_ptr = dtqsat.data();

    amrex::ParallelFor(temperatures.size(), [=] AMREX_GPU_DEVICE (int index) noexcept {
        const amrex::Real temperature = temperature_ptr[index];
        const amrex::Real pressure = pressure_ptr[index];
        amrex::Real qsat_value;
        amrex::Real dtqsat_value;

        dtesat_ptr[index] = erf_dtesati(temperature);
        erf_qsati(temperature, pressure, qsat_value);
        erf_dtqsati(temperature, pressure, dtqsat_value);
        qsat_ptr[index] = qsat_value;
        dtqsat_ptr[index] = dtqsat_value;
    });
}

void launch_gamma_sweep (amrex::Gpu::DeviceVector<amrex::Real>& arguments,
                         amrex::Gpu::DeviceVector<amrex::Real>& gamma_values)
{
    auto* argument_ptr = arguments.data();
    auto* gamma_ptr = gamma_values.data();

    amrex::ParallelFor(arguments.size(), [=] AMREX_GPU_DEVICE (int index) noexcept {
        gamma_ptr[index] = erf_gammafff(argument_ptr[index]);
    });
}

} // namespace

// Motivation: erf_gammafff still carries an AMREX_GPU_HOST_DEVICE annotation.
// Keep a cheap device-path check so that annotation stays exercised.
TEST(MicrophysicsKernel, GammaHostDeviceEquivalence)
{
    const std::vector<amrex::Real> sample_arguments = {
        amrex::Real(0.5), amrex::Real(1.0), amrex::Real(1.5), amrex::Real(2.5)};

    amrex::Gpu::HostVector<amrex::Real> host_arguments = make_host_vector(sample_arguments);
    amrex::Gpu::DeviceVector<amrex::Real> device_arguments(host_arguments.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_gamma(host_arguments.size());
    amrex::Gpu::HostVector<amrex::Real> host_gamma(host_arguments.size(), amrex::Real(0.0));

    copy_host_to_device(host_arguments, device_arguments);
    launch_gamma_sweep(device_arguments, device_gamma);
    amrex::Gpu::streamSynchronize();
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_gamma.begin(), device_gamma.end(), host_gamma.begin());

    for (int index = 0; index < static_cast<int>(host_arguments.size()); ++index) {
        const amrex::Real argument = host_arguments[index];
        expect_sample_match("erf_gammafff", index, argument, amrex::Real(-1.0),
                            erf_gammafff(argument), host_gamma[index], kValueRelTol);
    }
}

// Motivation: This is the main device-path coverage for the microphysics
// utilities. It sweeps cold, warm, and branch-window states and compares the
// AMReX ParallelFor results against host references sample by sample.
TEST(MicrophysicsKernel, WaterSweptHostDeviceEquivalence)
{
    const auto temperatures_1d = make_water_temperatures();
    const auto pressures_1d = make_pressures();
    std::vector<amrex::Real> sample_temperatures;
    std::vector<amrex::Real> sample_pressures;
    sample_temperatures.reserve(temperatures_1d.size() * pressures_1d.size());
    sample_pressures.reserve(temperatures_1d.size() * pressures_1d.size());

    for (const amrex::Real temperature : temperatures_1d) {
        for (const amrex::Real pressure : pressures_1d) {
            sample_temperatures.push_back(temperature);
            sample_pressures.push_back(pressure);
        }
    }

    amrex::Gpu::HostVector<amrex::Real> host_temperatures = make_host_vector(sample_temperatures);
    amrex::Gpu::HostVector<amrex::Real> host_pressures = make_host_vector(sample_pressures);
    amrex::Gpu::DeviceVector<amrex::Real> device_temperatures(host_temperatures.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_pressures(host_pressures.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_esat(host_temperatures.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_dtesat(host_temperatures.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_qsat(host_temperatures.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_dtqsat(host_temperatures.size());
    amrex::Gpu::HostVector<amrex::Real> host_esat(host_temperatures.size(), amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_dtesat(host_temperatures.size(), amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_qsat(host_temperatures.size(), amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_dtqsat(host_temperatures.size(), amrex::Real(0.0));

    copy_host_to_device(host_temperatures, device_temperatures);
    copy_host_to_device(host_pressures, device_pressures);
    launch_water_sweep(device_temperatures, device_pressures, device_esat, device_dtesat, device_qsat, device_dtqsat);
    amrex::Gpu::streamSynchronize();
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_esat.begin(), device_esat.end(), host_esat.begin());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_dtesat.begin(), device_dtesat.end(), host_dtesat.begin());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_qsat.begin(), device_qsat.end(), host_qsat.begin());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_dtqsat.begin(), device_dtqsat.end(), host_dtqsat.begin());

    for (int index = 0; index < static_cast<int>(host_temperatures.size()); ++index) {
        const amrex::Real temperature = host_temperatures[index];
        const amrex::Real pressure = host_pressures[index];
        amrex::Real qsat;
        amrex::Real dtqsat;
        erf_qsatw(temperature, pressure, qsat);
        erf_dtqsatw(temperature, pressure, dtqsat);

        expect_sample_match("erf_esatw", index, temperature, pressure,
                            erf_esatw(temperature), host_esat[index], kValueRelTol);
        expect_sample_match("erf_dtesatw", index, temperature, pressure,
                            erf_dtesatw(temperature), host_dtesat[index], kDerivativeRelTol);
        expect_sample_match("erf_qsatw", index, temperature, pressure,
                            qsat, host_qsat[index], kValueRelTol);
        expect_sample_match("erf_dtqsatw", index, temperature, pressure,
                            dtqsat, host_dtqsat[index], kDerivativeRelTol);
    }
}

// Motivation: Ice saturation value coverage reaches the colder 183.16 K limit,
// while the derivative-based sweep below handles only the narrower derivative
// contract.
TEST(MicrophysicsKernel, IceValueSweptHostDeviceEquivalence)
{
    std::vector<amrex::Real> sample_temperatures;
    sample_temperatures.reserve(91);
    for (int index = 0; index < 90; ++index) {
        sample_temperatures.push_back(amrex::Real(183.16) + amrex::Real(1.0e-3) + amrex::Real(index));
    }
    sample_temperatures.push_back(amrex::Real(273.16));

    amrex::Gpu::HostVector<amrex::Real> host_temperatures = make_host_vector(sample_temperatures);
    amrex::Gpu::DeviceVector<amrex::Real> device_temperatures(host_temperatures.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_esat(host_temperatures.size());
    amrex::Gpu::HostVector<amrex::Real> host_esat(host_temperatures.size(), amrex::Real(0.0));

    copy_host_to_device(host_temperatures, device_temperatures);
    launch_ice_value_sweep(device_temperatures, device_esat);
    amrex::Gpu::streamSynchronize();
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_esat.begin(), device_esat.end(), host_esat.begin());

    for (int index = 0; index < static_cast<int>(host_temperatures.size()); ++index) {
        const amrex::Real temperature = host_temperatures[index];
        expect_sample_match("erf_esati", index, temperature, amrex::Real(-1.0),
                            erf_esati(temperature), host_esat[index], kValueRelTol);
    }
}

// Motivation: This sweep keeps the derivative and qsat device coverage inside
// the narrower ice derivative contract while also covering capped and uncapped
// qsat behavior across the representative pressure range.
TEST(MicrophysicsKernel, IceDerivativeAndQSatSweptHostDeviceEquivalence)
{
    const auto pressures_1d = make_pressures();
    std::vector<amrex::Real> sample_temperatures;
    std::vector<amrex::Real> sample_pressures;
    sample_temperatures.reserve(86 * pressures_1d.size());
    sample_pressures.reserve(86 * pressures_1d.size());

    std::vector<amrex::Real> derivative_temperatures;
    derivative_temperatures.reserve(86);
    for (int offset = 0; offset < 85; ++offset) {
        derivative_temperatures.push_back(amrex::Real(188.16) + amrex::Real(1.0e-3) + amrex::Real(offset));
    }
    derivative_temperatures.push_back(amrex::Real(273.16));

    for (const amrex::Real temperature : derivative_temperatures) {
        for (const amrex::Real pressure : pressures_1d) {
            sample_temperatures.push_back(temperature);
            sample_pressures.push_back(pressure);
        }
    }

    amrex::Gpu::HostVector<amrex::Real> host_temperatures = make_host_vector(sample_temperatures);
    amrex::Gpu::HostVector<amrex::Real> host_pressures = make_host_vector(sample_pressures);
    amrex::Gpu::DeviceVector<amrex::Real> device_temperatures(host_temperatures.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_pressures(host_pressures.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_dtesat(host_temperatures.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_qsat(host_temperatures.size());
    amrex::Gpu::DeviceVector<amrex::Real> device_dtqsat(host_temperatures.size());
    amrex::Gpu::HostVector<amrex::Real> host_dtesat(host_temperatures.size(), amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_qsat(host_temperatures.size(), amrex::Real(0.0));
    amrex::Gpu::HostVector<amrex::Real> host_dtqsat(host_temperatures.size(), amrex::Real(0.0));

    copy_host_to_device(host_temperatures, device_temperatures);
    copy_host_to_device(host_pressures, device_pressures);
    launch_ice_derivative_sweep(device_temperatures, device_pressures, device_dtesat, device_qsat, device_dtqsat);
    amrex::Gpu::streamSynchronize();
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_dtesat.begin(), device_dtesat.end(), host_dtesat.begin());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_qsat.begin(), device_qsat.end(), host_qsat.begin());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_dtqsat.begin(), device_dtqsat.end(), host_dtqsat.begin());

    for (int index = 0; index < static_cast<int>(host_temperatures.size()); ++index) {
        const amrex::Real temperature = host_temperatures[index];
        const amrex::Real pressure = host_pressures[index];
        amrex::Real qsat;
        amrex::Real dtqsat;
        erf_qsati(temperature, pressure, qsat);
        erf_dtqsati(temperature, pressure, dtqsat);

        expect_sample_match("erf_dtesati", index, temperature, pressure,
                            erf_dtesati(temperature), host_dtesat[index], kDerivativeRelTol);
        expect_sample_match("erf_qsati", index, temperature, pressure,
                            qsat, host_qsat[index], kValueRelTol);
        expect_sample_match("erf_dtqsati", index, temperature, pressure,
                            dtqsat, host_dtqsat[index], kDerivativeRelTol);
    }
}
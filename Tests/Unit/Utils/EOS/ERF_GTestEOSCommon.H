#ifndef ERF_GTEST_EOS_COMMON_H_
#define ERF_GTEST_EOS_COMMON_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <string>

#include <AMReX_Array.H>
#include <AMReX_Math.H>

#include <ERF_Constants.H>
#include <ERF_EOS.H>

namespace eos_test {

// EOS unit-test contract
// ----------------------
// These tests check algebraic and thermodynamic properties of ERF's EOS,
// not hard-coded gold values.
//
// Core invariants:
//   1. getThgivenTandP and getTgivenPandTh are inverses at fixed pressure.
//   2. getPgivenRTh and getRhoThetagivenP are inverses for dry density,
//      dry potential temperature, and water-vapor mixing ratio qv.
//   3. getRhogivenTandPress and getRhogivenThetaPress agree when theta
//      is diagnosed from the same T and p.
//   4. getTgivenRandRTh and getThgivenRandT close the same EOS.
//   5. getExnergivenP and getExnergivenRTh agree.
//   6. getdPdRgivenConstantTheta equals Gamma*P/rho at fixed theta.
//
// Thermodynamics:
//   - rho is dry density rho_d.
//   - theta is dry potential temperature theta_d.
//   - qv is vapor mixing ratio per dry-air mass.
//   - EOS pressure and temperature functions use pressure in Pa.
//   - RH in compute_vapor_pressure is fractional RH.
//
// Portability:
//   - These tests must build and run in CPU and GPU configurations.
//   - Kernel-path tests call EOS functions from AMReX ParallelFor.
//   - GTest EXPECT/ASSERT calls stay on the host.
//   - Device kernels compute values or normalized errors only.

#ifdef AMREX_USE_FLOAT
constexpr amrex::Real kEOSRelTol = amrex::Real(2.0e-5);
constexpr amrex::Real kKernelErrorLimit = amrex::Real(20.0);
#else
constexpr amrex::Real kEOSRelTol = amrex::Real(1.0e-12);
constexpr amrex::Real kKernelErrorLimit = amrex::Real(10.0);
#endif

constexpr amrex::Real kRdOcp = R_d / Cp_d;
constexpr int kNumEOSStates = 4;
constexpr int kKernelMetricsPerState = 6;

struct EOSState {
    amrex::Real temperature;
    amrex::Real pressure;
    amrex::Real qv;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
EOSState get_test_state (const int index)
{
    switch (index) {
    case 0:
        return EOSState{amrex::Real(300.0), amrex::Real(90000.0), amrex::Real(0.0)};
    case 1:
        return EOSState{amrex::Real(300.0), amrex::Real(90000.0), amrex::Real(0.01)};
    case 2:
        return EOSState{amrex::Real(260.0), amrex::Real(60000.0), amrex::Real(0.002)};
    default:
        return EOSState{amrex::Real(285.0), amrex::Real(80000.0), amrex::Real(0.006)};
    }
}

inline std::array<EOSState, kNumEOSStates> host_test_states ()
{
    return {get_test_state(0), get_test_state(1), get_test_state(2), get_test_state(3)};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real moisture_factor (const amrex::Real qv)
{
    return amrex::Real(1.0) + (R_v / R_d) * qv;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real portable_scaled_tol (const amrex::Real actual,
                                 const amrex::Real expected,
                                 const amrex::Real factor)
{
    return factor * amrex::max(amrex::Real(1.0),
                                amrex::max(amrex::Math::abs(actual), amrex::Math::abs(expected)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real normalized_error (const amrex::Real actual,
                              const amrex::Real expected,
                              const amrex::Real factor)
{
    return amrex::Math::abs(actual - expected) / portable_scaled_tol(actual, expected, factor);
}

inline amrex::Real scaled_tol (const amrex::Real actual,
                               const amrex::Real expected,
                               const amrex::Real factor)
{
    return factor * std::max({amrex::Real(1.0), std::abs(actual), std::abs(expected)});
}

inline std::string state_label (const int index, const EOSState& state)
{
    return "state[" + std::to_string(index) + "] T=" + std::to_string(state.temperature) +
           " p=" + std::to_string(state.pressure) + " qv=" + std::to_string(state.qv);
}

} // namespace eos_test

#endif
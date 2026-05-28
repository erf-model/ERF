#ifndef ERF_GTEST_KESSLER_COMMON_H_
#define ERF_GTEST_KESSLER_COMMON_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MFIter.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <ERF_Constants.H>
#include <ERF_DataStruct.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_Kessler.H>
#include <ERF_KesslerUtils.H>
#include <ERF_MicrophysicsUtils.H>

namespace kessler_test {

// Kessler unit-test contract
// --------------------------
// These tests check ERF's refactored Kessler helper surface and public flow.
// The scalar tests cover functional-form contracts for the documented Kessler
// closure formulas. The physical-property tests cover closure-independent
// invariants that any reasonable warm-rain scheme should satisfy.
//
// Core invariants:
//   1. Local source terms conserve qv + qc + qp before sedimentation.
//   2. Sedimentation transfers precipitating water vertically and to
//      rain_accum at the bottom.
//   3. Public-flow column water budgets must include surface rain.
//   4. Kessler_NoRain does not advance rain water or rain accumulation.
//   5. qsat helper pressure inputs are in mbar / hPa.
//   6. EOS pressure utilities use Pa.
//
// Portability:
//   - Device kernels compute helper outputs or reduced error metrics only.
//   - GTest EXPECT/ASSERT calls stay on the host.
//   - No dynamic allocation is allowed inside device code.
//   - Device coverage uses representative case tables, not single-sample smoke
//     tests.
//
// Helper composition:
//   - kessler_warm_rain_sources calls kessler_saturation_adjustment internally.
//   - Tests that isolate rain evaporation must choose inputs where the
//     saturation-adjustment terms are known, or account for the full composed
//     source update.

constexpr amrex::Real kRdOcp = R_d / Cp_d;
constexpr amrex::Real kSatAdjLatentOverCp = L_v / Cp_d;
constexpr amrex::Real kThetaLatentOverCp = L_v / Cp_d;
constexpr amrex::Real kDefaultDt = amrex::Real(1.0);

enum StateDiffComp {
    StateDiffRho = 0,
    StateDiffRhoTheta,
    StateDiffQ1,
    StateDiffQ2,
    StateDiffQ3,
    NumStateDiffComps
};

enum ConservationErrComp {
    ConservationErrRho = 0,
    ConservationErrTotalWater,
    ConservationErrLatentProxy,
    NumConservationErrComps
};

enum RainStateErrComp {
    RainStateErrQp = 0,
    RainStateErrRainAccum,
    NumRainStateErrComps
};

struct PrimitiveState {
    amrex::Real rho;
    amrex::Real theta;
    amrex::Real tabs;
    amrex::Real pres_mbar;
    amrex::Real qv;
    amrex::Real qc;
    amrex::Real qp;
};

struct ColumnBudget {
    amrex::Real rho_sum;
    amrex::Real total_water_sum;
    amrex::Real latent_proxy_sum;
    amrex::Real rain_accum_sum;
};

struct KernelCase {
    amrex::Real qv;
    amrex::Real qc;
    amrex::Real qp;
    amrex::Real rho;
    amrex::Real pressure_mbar;
    amrex::Real qsat;
    amrex::Real dtqsat;
    amrex::Real dt;
    bool do_cond;
    int face_k;
    int face_k_lo;
    int face_k_hi;
    amrex::Real face_rho_km1;
    amrex::Real face_rho_k;
    amrex::Real face_qp_km1;
    amrex::Real face_qp_k;
    amrex::Real fz_hi;
    amrex::Real fz_lo;
    amrex::Real dJinv;
    amrex::Real coef;
    amrex::Real reduced_value;
    amrex::Real substep_dt;
    amrex::Real substep_dz;
};

// Machine-precision-scaled absolute tolerance for analytic formula checks.
// Provenance: machine epsilon, float-vs-double scaling, and a small safety
// factor for one to a few arithmetic operations.
inline amrex::Real formula_abs_tol (const amrex::Real expected)
{
#ifdef AMREX_USE_FLOAT
    return amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon() *
           std::max(amrex::Real(1.0), std::abs(expected));
#else
    return amrex::Real(32.0) * std::numeric_limits<amrex::Real>::epsilon() *
           std::max(amrex::Real(1.0), std::abs(expected));
#endif
}

// Relative tolerance for smooth closed-form identities.
// Provenance: machine epsilon, float-vs-double scaling, and a safety factor for
// several multiplies/divides without transcendental functions.
inline amrex::Real formula_rel_tol (const amrex::Real expected)
{
    return formula_abs_tol(expected);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real portable_formula_abs_tol (const amrex::Real expected) noexcept
{
#ifdef AMREX_USE_FLOAT
    const amrex::Real factor = amrex::Real(64.0);
#else
    const amrex::Real factor = amrex::Real(32.0);
#endif
    return factor * std::numeric_limits<amrex::Real>::epsilon() *
           amrex::max(amrex::Real(1.0), amrex::Math::abs(expected));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real portable_formula_rel_tol (const amrex::Real expected) noexcept
{
    return portable_formula_abs_tol(expected);
}

// Absolute tolerance for pow-based helper comparisons.
// Provenance: backend math-library variation for pow on CPU/GPU, floating-point
// roundoff, and float-vs-double scaling.
inline amrex::Real backend_math_abs_tol (const amrex::Real expected)
{
#ifdef AMREX_USE_FLOAT
    return amrex::Real(4.0e3) * std::numeric_limits<amrex::Real>::epsilon() *
           std::max(amrex::Real(1.0), std::abs(expected));
#else
    return amrex::Real(2.0e2) * std::numeric_limits<amrex::Real>::epsilon() *
           std::max(amrex::Real(1.0), std::abs(expected));
#endif
}

// Tolerance for exact-zero branches and clipped values.
// Provenance: machine epsilon only; use this when the mathematical contract is
// exactly zero but host/device roundoff may leave a tiny residual.
inline amrex::Real exact_zero_or_near_zero_tol ()
{
#ifdef AMREX_USE_FLOAT
    return amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon();
#else
    return amrex::Real(32.0) * std::numeric_limits<amrex::Real>::epsilon();
#endif
}

// Finite-difference tolerance for smooth derivative tests.
// Provenance: O(h^2) truncation for a central difference, machine epsilon,
// roundoff amplification by division through h, and float-vs-double scaling.
inline amrex::Real finite_difference_tol (const amrex::Real scale,
                                          const amrex::Real h)
{
#ifdef AMREX_USE_FLOAT
    const amrex::Real truncation = amrex::Real(8.0) * std::abs(scale) * h * h;
    const amrex::Real roundoff = amrex::Real(256.0) * std::numeric_limits<amrex::Real>::epsilon() *
                                 std::max(amrex::Real(1.0), std::abs(scale)) / h;
#else
    const amrex::Real truncation = amrex::Real(4.0) * std::abs(scale) * h * h;
    const amrex::Real roundoff = amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon() *
                                 std::max(amrex::Real(1.0), std::abs(scale)) / h;
#endif
    return std::max({exact_zero_or_near_zero_tol(), truncation, roundoff});
}

// Column-budget tolerance for sums over many cells and sedimentation substeps.
// Provenance: machine epsilon, float-vs-double scaling, and accumulation over
// N terms with a conservative safety factor.
inline amrex::Real property_accumulation_tol (const int num_terms,
                                              const amrex::Real scale)
{
#ifdef AMREX_USE_FLOAT
    const amrex::Real factor = amrex::Real(256.0);
#else
    const amrex::Real factor = amrex::Real(64.0);
#endif
    return factor * std::numeric_limits<amrex::Real>::epsilon() *
           std::max(amrex::Real(1.0), std::abs(scale)) * std::max(1, num_terms);
}

inline amrex::Real derivative_step (const amrex::Real scale)
{
#ifdef AMREX_USE_FLOAT
    return amrex::Real(2.0e-3) * std::max(amrex::Real(1.0), std::abs(scale));
#else
    return amrex::Real(1.0e-6) * std::max(amrex::Real(1.0), std::abs(scale));
#endif
}

inline amrex::Real scaled_tol (const amrex::Real actual,
                               const amrex::Real expected,
                               const amrex::Real factor)
{
    return factor * std::max({amrex::Real(1.0), std::abs(actual), std::abs(expected)});
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

template <typename Function>
inline amrex::Real central_difference (Function&& func,
                                       const amrex::Real x,
                                       const amrex::Real h)
{
    return (func(x + h) - func(x - h)) / (amrex::Real(2.0) * h);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
KesslerSaturationAdjustment reference_saturation_adjustment (const amrex::Real qv,
                                                             const amrex::Real qc,
                                                             const amrex::Real qsat,
                                                             const amrex::Real dtqsat,
                                                             const bool do_cond,
                                                             const amrex::Real latent_over_cp) noexcept
{
    KesslerSaturationAdjustment source_terms{amrex::Real(0.0), amrex::Real(0.0)};
    const amrex::Real denom = amrex::Real(1.0) + latent_over_cp * dtqsat;

    if ((qv > qsat) && do_cond) {
        source_terms.dq_vapor_to_cloud = amrex::min(qv, (qv - qsat) / denom);
    }
    if ((qv < qsat) && (qc > amrex::Real(0.0)) && do_cond) {
        source_terms.dq_cloud_to_vapor = amrex::min(qc, (qsat - qv) / denom);
    }

    return source_terms;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
KesslerSourceTerms reference_warm_rain_sources (const amrex::Real qv,
                                                const amrex::Real qc,
                                                const amrex::Real qp,
                                                const amrex::Real rho,
                                                const amrex::Real pressure_mbar,
                                                const amrex::Real qsat,
                                                const amrex::Real dtqsat,
                                                const amrex::Real dt,
                                                const bool do_cond,
                                                const amrex::Real latent_over_cp) noexcept
{
    const KesslerSaturationAdjustment sat =
        reference_saturation_adjustment(qv, qc, qsat, dtqsat, do_cond, latent_over_cp);
    KesslerSourceTerms source_terms{
        sat.dq_vapor_to_cloud,
        sat.dq_cloud_to_vapor,
        amrex::Real(0.0),
        amrex::Real(0.0)};

    const amrex::Real capacity = amrex::max(
        amrex::Real(0.0), (qsat - qv) / (amrex::Real(1.0) + latent_over_cp * dtqsat));
    const amrex::Real remaining_capacity = amrex::max(
        amrex::Real(0.0), capacity - sat.dq_cloud_to_vapor);
    const amrex::Real qv_eff = qv - sat.dq_vapor_to_cloud + sat.dq_cloud_to_vapor;

    if ((qp > amrex::Real(0.0)) && (qsat > amrex::Real(0.0))
        && (pressure_mbar > amrex::Real(0.0))
        && (qv_eff < qsat) && (remaining_capacity > amrex::Real(0.0))) {
        // MUST MATCH: Kessler rain-evaporation closure coefficients.
        const amrex::Real scaled_mass = amrex::Real(0.001) * rho * qp;
        const amrex::Real coeff = amrex::Real(1.6)
            + amrex::Real(124.9) * std::pow(scaled_mass, amrex::Real(0.2046));
        const amrex::Real raw = amrex::Real(1.0) / (amrex::Real(0.001) * rho)
            * (amrex::Real(1.0) - qv_eff / qsat)
            * coeff
            * std::pow(scaled_mass, amrex::Real(0.525))
            / (amrex::Real(5.4e5) + amrex::Real(2.55e6) / (pressure_mbar * qsat))
            * dt;
        source_terms.dq_rain_to_vapor = amrex::min(qp, amrex::min(raw, remaining_capacity));
    }

    if (qc > amrex::Real(0.0)) {
        const amrex::Real available_cloud = amrex::max(
            amrex::Real(0.0), qc + sat.dq_vapor_to_cloud - sat.dq_cloud_to_vapor);
        // MUST MATCH: qcw0 and alphaelq Kessler autoconversion parameters.
        const amrex::Real auto_r = (qc > qcw0) ? alphaelq : amrex::Real(0.0);
        // MUST MATCH: Kessler accretion coefficient and exponent.
        const amrex::Real accretion_rate = amrex::Real(2.2) * std::pow(qp, amrex::Real(0.875));
        const amrex::Real raw = dt * (accretion_rate * qc + auto_r * (qc - qcw0));
        source_terms.dq_cloud_to_rain = amrex::min(raw, available_cloud);
    }

    return source_terms;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void apply_local_sources (amrex::Real& qv,
                          amrex::Real& qc,
                          amrex::Real& qp,
                          const KesslerSourceTerms& source_terms) noexcept
{
    qv += -source_terms.dq_vapor_to_cloud
        +  source_terms.dq_cloud_to_vapor
        +  source_terms.dq_rain_to_vapor;
    qc +=  source_terms.dq_vapor_to_cloud
        -  source_terms.dq_cloud_to_vapor
        -  source_terms.dq_cloud_to_rain;
    qp +=  source_terms.dq_cloud_to_rain
        -  source_terms.dq_rain_to_vapor;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real reference_terminal_velocity (const amrex::Real rho,
                                         const amrex::Real qp) noexcept
{
    // MUST MATCH: Kessler terminal-velocity coefficient and exponents.
    return amrex::Real(36.34)
        * std::pow(amrex::Real(0.001) * rho * qp, amrex::Real(0.1346))
        * std::pow(rho / amrex::Real(1.16), -myhalf);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real reference_terminal_velocity_dqp (const amrex::Real rho,
                                             const amrex::Real qp) noexcept
{
    const amrex::Real velocity = reference_terminal_velocity(rho, qp);
    return amrex::Real(0.1346) * velocity / qp;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real reference_terminal_velocity_drho (const amrex::Real rho,
                                              const amrex::Real qp) noexcept
{
    const amrex::Real velocity = reference_terminal_velocity(rho, qp);
    return (amrex::Real(0.1346) - myhalf) * velocity / rho;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real inverse_qp_from_terminal_velocity (const amrex::Real velocity,
                                               const amrex::Real rho) noexcept
{
    // MUST MATCH: invert V = 36.34 * (0.001*rho*qp)^0.1346 * (rho/1.16)^(-1/2).
    const amrex::Real density_factor = std::pow(rho / amrex::Real(1.16), myhalf);
    const amrex::Real scaled = velocity * density_factor / amrex::Real(36.34);
    return std::pow(scaled, amrex::Real(1.0) / amrex::Real(0.1346)) /
           (amrex::Real(0.001) * rho);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real reference_precip_flux (const amrex::Real rho,
                                   const amrex::Real velocity,
                                   const amrex::Real qp) noexcept
{
    return rho * velocity * qp;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real reference_precip_flux_dqp (const amrex::Real rho,
                                       const amrex::Real qp) noexcept
{
    const amrex::Real flux = reference_precip_flux(rho, reference_terminal_velocity(rho, qp), qp);
    return (amrex::Real(1.0) + amrex::Real(0.1346)) * flux / qp;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real reference_precip_flux_drho (const amrex::Real rho,
                                        const amrex::Real qp) noexcept
{
    const amrex::Real flux = reference_precip_flux(rho, reference_terminal_velocity(rho, qp), qp);
    return (myhalf + amrex::Real(0.1346)) * flux / rho;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real reference_sedimentation_tendency (const amrex::Real fz_hi,
                                              const amrex::Real fz_lo,
                                              const amrex::Real rho,
                                              const amrex::Real dJinv,
                                              const amrex::Real coef) noexcept
{
    return dJinv * (amrex::Real(1.0) / rho) * (fz_hi - fz_lo) * coef;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int reference_substeps_from_reduced_value (const amrex::Real reduced_value,
                                           const amrex::Real dt,
                                           const amrex::Real dz) noexcept
{
    return static_cast<int>(std::ceil((reduced_value + std::numeric_limits<amrex::Real>::epsilon())
                                      * (dt / dz) / kessler_sedimentation_cfl_max));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int reference_velocity_cfl_substeps (const amrex::Real terminal_velocity,
                                     const amrex::Real dt,
                                     const amrex::Real dz) noexcept
{
    return static_cast<int>(std::ceil((terminal_velocity + std::numeric_limits<amrex::Real>::epsilon())
                                      * (dt / dz) / kessler_sedimentation_cfl_max));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
PrimitiveState make_primitive_state (const amrex::Real tabs,
                                     const amrex::Real pres_mbar,
                                     const amrex::Real qv,
                                     const amrex::Real qc,
                                     const amrex::Real qp)
{
    const amrex::Real pres_pa = amrex::Real(100.0) * pres_mbar;
    const amrex::Real rho = getRhogivenTandPress(tabs, pres_pa, qv);
    const amrex::Real theta = getThgivenTandP(tabs, pres_pa, kRdOcp);
    return PrimitiveState{rho, theta, tabs, pres_mbar, qv, qc, qp};
}

inline void sync ()
{
#ifdef AMREX_USE_GPU
    amrex::Gpu::streamSynchronize();
#endif
}

template <typename F>
inline void run_and_sync (F&& fn)
{
    std::forward<F>(fn)();
    sync();
}

inline amrex::Geometry make_geometry (const int nx,
                                      const int ny,
                                      const int nz)
{
    amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
                      amrex::IntVect(AMREX_D_DECL(nx - 1, ny - 1, nz - 1)));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(static_cast<amrex::Real>(nx),
                                          static_cast<amrex::Real>(ny),
                                          static_cast<amrex::Real>(nz))});
    amrex::Array<int, AMREX_SPACEDIM> periodicity{AMREX_D_DECL(1, 1, 0)};
    return amrex::Geometry(domain, &real_box, amrex::CoordSys::cartesian, periodicity.data());
}

inline SolverChoice make_solver_choice (const MoistureType moisture_type,
                                        const bool use_shoc)
{
    SolverChoice sc;
    sc.c_p = Cp_d;
    sc.rdOcp = kRdOcp;
    sc.ave_plane = 2;
    sc.moisture_type = moisture_type;
    sc.use_shoc = use_shoc;

    if (moisture_type == MoistureType::Kessler) {
        sc.moisture_indices = MoistureComponentIndices(RhoQ1_comp, RhoQ2_comp, -1, RhoQ3_comp);
    } else if (moisture_type == MoistureType::Kessler_NoRain) {
        sc.moisture_indices = MoistureComponentIndices(RhoQ1_comp, RhoQ2_comp);
    }

    return sc;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_conserved_cell (const amrex::Array4<amrex::Real>& arr,
                         const int i,
                         const int j,
                         const int k,
                         const PrimitiveState& state)
{
    arr(i,j,k,Rho_comp) = state.rho;
    arr(i,j,k,RhoTheta_comp) = state.rho * state.theta;
    arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
    arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
    arr(i,j,k,RhoQ1_comp) = state.rho * state.qv;
    arr(i,j,k,RhoQ2_comp) = state.rho * state.qc;
    arr(i,j,k,RhoQ3_comp) = state.rho * state.qp;
}

inline void fill_conserved_state_portable (amrex::MultiFab& cons)
{
    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto arr = cons.array(mfi);
        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const int selector = (i + 2 * j + 3 * k) % 4;
                PrimitiveState state;
                if (selector == 0) {
                    state = make_primitive_state(amrex::Real(289.0), amrex::Real(900.0),
                                                 amrex::Real(1.10e-2), amrex::Real(6.0e-4), amrex::Real(0.0));
                } else if (selector == 1) {
                    state = make_primitive_state(amrex::Real(287.0), amrex::Real(880.0),
                                                 amrex::Real(8.4e-3), amrex::Real(1.8e-3), amrex::Real(4.0e-4));
                } else if (selector == 2) {
                    state = make_primitive_state(amrex::Real(285.0), amrex::Real(860.0),
                                                 amrex::Real(7.0e-3), amrex::Real(5.0e-4), amrex::Real(1.5e-3));
                } else {
                    state = make_primitive_state(amrex::Real(291.0), amrex::Real(930.0),
                                                 amrex::Real(1.20e-2), amrex::Real(1.0e-4), amrex::Real(2.5e-4));
                }
                set_conserved_cell(arr, i, j, k, state);
            });
        });
    }
}

// Keep device lambdas in namespace-scope helpers. CUDA rejects extended
// device lambdas enclosed by GTest's private TestBody() member function.
inline void fill_conserved_state_uniform_portable (amrex::MultiFab& cons,
                                                   const PrimitiveState& state)
{
    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto arr = cons.array(mfi);
        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                set_conserved_cell(arr, i, j, k, state);
            });
        });
    }
}

inline void fill_local_source_only_state_portable (amrex::MultiFab& cons)
{
    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto arr = cons.array(mfi);
        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const int selector = (i + 2 * j + 3 * k) % 4;
                PrimitiveState state;
                amrex::Real qsat_local;
                if (selector == 0) {
                    erf_qsatw(amrex::Real(290.0), amrex::Real(900.0), qsat_local);
                    state = make_primitive_state(amrex::Real(290.0), amrex::Real(900.0),
                                                 qsat_local + amrex::Real(2.0e-4), amrex::Real(5.0e-5), amrex::Real(0.0));
                } else if (selector == 1) {
                    erf_qsatw(amrex::Real(289.0), amrex::Real(900.0), qsat_local);
                    state = make_primitive_state(amrex::Real(289.0), amrex::Real(900.0),
                                                 qsat_local - amrex::Real(2.0e-4), amrex::Real(8.0e-5), amrex::Real(0.0));
                } else if (selector == 2) {
                    erf_qsatw(amrex::Real(286.0), amrex::Real(880.0), qsat_local);
                    state = make_primitive_state(amrex::Real(286.0), amrex::Real(880.0),
                                                 qsat_local, amrex::Real(1.0e-4), amrex::Real(0.0));
                } else {
                    erf_qsatw(amrex::Real(288.0), amrex::Real(910.0), qsat_local);
                    state = make_primitive_state(amrex::Real(288.0), amrex::Real(910.0),
                                                 qsat_local + amrex::Real(1.0e-4), amrex::Real(0.0), amrex::Real(0.0));
                }
                set_conserved_cell(arr, i, j, k, state);
            });
        });
    }
}

inline void fill_conserved_column_state_portable (amrex::MultiFab& cons)
{
    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto arr = cons.array(mfi);
        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                PrimitiveState state;
                if (k == 0) {
                    state = make_primitive_state(amrex::Real(289.0), amrex::Real(900.0),
                                                 amrex::Real(1.02e-2), amrex::Real(4.0e-4), amrex::Real(1.2e-3));
                } else if (k == 1) {
                    state = make_primitive_state(amrex::Real(287.0), amrex::Real(880.0),
                                                 amrex::Real(9.0e-3), amrex::Real(8.0e-4), amrex::Real(1.6e-3));
                } else if (k == 2) {
                    state = make_primitive_state(amrex::Real(285.0), amrex::Real(860.0),
                                                 amrex::Real(8.4e-3), amrex::Real(5.0e-4), amrex::Real(9.0e-4));
                } else {
                    // Top-cell qp is zero so the top precipitation flux term vanishes.
                    state = make_primitive_state(amrex::Real(283.0), amrex::Real(840.0),
                                                 amrex::Real(7.8e-3), amrex::Real(2.0e-4), amrex::Real(0.0));
                }
                set_conserved_cell(arr, i, j, k, state);
            });
        });
    }
}

inline void fill_detj_portable (amrex::MultiFab& detj)
{
    for (amrex::MFIter mfi(detj, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto arr = detj.array(mfi);
        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                arr(i,j,k) = amrex::Real(1.0) + amrex::Real(0.1) * (i + 2 * j + 3 * k);
            });
        });
    }
}

inline void fill_compressed_detj_portable (amrex::MultiFab& detj)
{
    for (amrex::MFIter mfi(detj, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto arr = detj.array(mfi);
        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                arr(i,j,k) = amrex::Real(0.25) + amrex::Real(0.25) * k;
            });
        });
    }
}

inline void fill_do_cond_false_rain_path_state_portable (amrex::MultiFab& cons)
{
    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto arr = cons.array(mfi);
        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const PrimitiveState state = (k == 3)
                    ? make_primitive_state(amrex::Real(284.0), amrex::Real(850.0),
                                           amrex::Real(7.0e-3), amrex::Real(2.0e-4), amrex::Real(0.0))
                    : make_primitive_state(amrex::Real(288.0), amrex::Real(900.0),
                                           amrex::Real(1.20e-2), qcw0 + amrex::Real(3.0e-4), amrex::Real(1.0e-3));
                set_conserved_cell(arr, i, j, k, state);
            });
        });
    }
}

inline void run_kessler_public_flow (Kessler& kessler,
                                     SolverChoice& sc,
                                     const amrex::Geometry& geom,
                                     amrex::MultiFab& cons,
                                     std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                                     std::unique_ptr<amrex::MultiFab>& detJ_cc,
                                     const amrex::Real dt = kDefaultDt)
{
    run_and_sync([&] {
        kessler.Define(sc);
        kessler.Set_dzmin(geom.CellSize(2));
        kessler.Init(cons, cons.boxArray(), geom, dt, z_phys_nd, detJ_cc);
        kessler.Copy_State_to_Micro(cons);
        kessler.Advance(dt, sc);
        kessler.Copy_Micro_to_State(cons);
    });
}

inline void copy_multifab (const amrex::MultiFab& src,
                           amrex::MultiFab& dst)
{
    amrex::MultiFab::Copy(dst, src, 0, 0, src.nComp(), src.nGrowVect());
}

inline void compute_kessler_state_difference (const amrex::MultiFab& after,
                                              const amrex::MultiFab& before,
                                              amrex::MultiFab& err)
{
    for (amrex::MFIter mfi(after, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto after_arr = after.const_array(mfi);
        auto before_arr = before.const_array(mfi);
        auto err_arr = err.array(mfi);
        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                err_arr(i,j,k,StateDiffRho) = amrex::Math::abs(after_arr(i,j,k,Rho_comp)
                                                              - before_arr(i,j,k,Rho_comp));
                err_arr(i,j,k,StateDiffRhoTheta) = amrex::Math::abs(after_arr(i,j,k,RhoTheta_comp)
                                                                   - before_arr(i,j,k,RhoTheta_comp));
                err_arr(i,j,k,StateDiffQ1) = amrex::Math::abs(after_arr(i,j,k,RhoQ1_comp)
                                                             - before_arr(i,j,k,RhoQ1_comp));
                err_arr(i,j,k,StateDiffQ2) = amrex::Math::abs(after_arr(i,j,k,RhoQ2_comp)
                                                             - before_arr(i,j,k,RhoQ2_comp));
                err_arr(i,j,k,StateDiffQ3) = amrex::Math::abs(after_arr(i,j,k,RhoQ3_comp)
                                                             - before_arr(i,j,k,RhoQ3_comp));
            });
        });
    }
}

inline void compute_kessler_conservation_errors (const amrex::MultiFab& before,
                                                 const amrex::MultiFab& after,
                                                 amrex::MultiFab& err)
{
    for (amrex::MFIter mfi(after, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto before_arr = before.const_array(mfi);
        auto after_arr = after.const_array(mfi);
        auto err_arr = err.array(mfi);
        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const amrex::Real rho_before = before_arr(i,j,k,Rho_comp);
                const amrex::Real rho_after = after_arr(i,j,k,Rho_comp);
                const amrex::Real qtot_before = before_arr(i,j,k,RhoQ1_comp)
                                              + before_arr(i,j,k,RhoQ2_comp)
                                              + before_arr(i,j,k,RhoQ3_comp);
                const amrex::Real qtot_after = after_arr(i,j,k,RhoQ1_comp)
                                             + after_arr(i,j,k,RhoQ2_comp)
                                             + after_arr(i,j,k,RhoQ3_comp);

                const amrex::Real theta_before = before_arr(i,j,k,RhoTheta_comp) / rho_before;
                const amrex::Real theta_after = after_arr(i,j,k,RhoTheta_comp) / rho_after;
                const amrex::Real qv_before = before_arr(i,j,k,RhoQ1_comp) / rho_before;
                const amrex::Real qv_after = after_arr(i,j,k,RhoQ1_comp) / rho_after;
                const amrex::Real pres_pa = getPgivenRTh(before_arr(i,j,k,RhoTheta_comp), qv_before);
                const amrex::Real tabs_before = getTgivenPandTh(pres_pa, theta_before, kRdOcp);
                const amrex::Real tabs_after = getTgivenPandTh(pres_pa, theta_after, kRdOcp);
                const amrex::Real latent_before = rho_before * (tabs_before + kThetaLatentOverCp * qv_before);
                const amrex::Real latent_after = rho_after * (tabs_after + kThetaLatentOverCp * qv_after);

                err_arr(i,j,k,ConservationErrRho) = normalized_error(
                    rho_after, rho_before, portable_formula_rel_tol(rho_before));
                err_arr(i,j,k,ConservationErrTotalWater) = normalized_error(
                    qtot_after, qtot_before, portable_formula_rel_tol(qtot_before));
                err_arr(i,j,k,ConservationErrLatentProxy) = normalized_error(
                    latent_after, latent_before, portable_formula_rel_tol(latent_before));
            });
        });
    }
}

inline void compute_rain_state_difference (const amrex::MultiFab& after,
                                           const amrex::MultiFab& before,
                                           const amrex::MultiFab& rain_after,
                                           const amrex::MultiFab& rain_before,
                                           amrex::MultiFab& err)
{
    for (amrex::MFIter mfi(after, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto after_arr = after.const_array(mfi);
        auto before_arr = before.const_array(mfi);
        auto rain_after_arr = rain_after.const_array(mfi);
        auto rain_before_arr = rain_before.const_array(mfi);
        auto err_arr = err.array(mfi);
        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                err_arr(i,j,k,RainStateErrQp) = amrex::Math::abs(after_arr(i,j,k,RhoQ3_comp)
                                                                - before_arr(i,j,k,RhoQ3_comp));
                err_arr(i,j,k,RainStateErrRainAccum) = amrex::Math::abs(rain_after_arr(i,j,k)
                                                                       - rain_before_arr(i,j,k));
            });
        });
    }
}

inline amrex::Real sum_multifab_component (const amrex::MultiFab& mf,
                                           const int comp)
{
    auto arrays = mf.const_arrays();
    amrex::GpuTuple<amrex::Real> reduced = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpSum>{},
        amrex::TypeList<amrex::Real>{},
        mf, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Real> {
            return {arrays[box_no](i,j,k,comp)};
        });
    return amrex::get<0>(reduced);
}

inline amrex::Real sum_multifab_scalar (const amrex::MultiFab& mf)
{
    auto arrays = mf.const_arrays();
    amrex::GpuTuple<amrex::Real> reduced = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpSum>{},
        amrex::TypeList<amrex::Real>{},
        mf, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Real> {
            return {arrays[box_no](i,j,k)};
        });
    return amrex::get<0>(reduced);
}

inline ColumnBudget compute_column_budget (const amrex::Geometry& geom,
                                           const amrex::MultiFab& cons,
                                           const amrex::MultiFab& rain_accum,
                                           const amrex::MultiFab* detJ_cc = nullptr)
{
    amrex::MultiFab rho_sum(cons.boxArray(), cons.DistributionMap(), 1, 0);
    amrex::MultiFab water_sum(cons.boxArray(), cons.DistributionMap(), 1, 0);
    amrex::MultiFab latent_sum(cons.boxArray(), cons.DistributionMap(), 1, 0);
    const amrex::Real dz = geom.CellSize(2);

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto cons_arr = cons.const_array(mfi);
        auto rho_arr = rho_sum.array(mfi);
        auto water_arr = water_sum.array(mfi);
        auto latent_arr = latent_sum.array(mfi);
        const auto detj_arr = detJ_cc ? detJ_cc->const_array(mfi) : amrex::Array4<const amrex::Real>{};

        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const amrex::Real rho = cons_arr(i,j,k,Rho_comp);
                const amrex::Real rhotheta = cons_arr(i,j,k,RhoTheta_comp);
                const amrex::Real qv = cons_arr(i,j,k,RhoQ1_comp) / rho;
                const amrex::Real theta = rhotheta / rho;
                const amrex::Real pressure_pa = getPgivenRTh(rhotheta, qv);
                const amrex::Real tabs = getTgivenPandTh(pressure_pa, theta, kRdOcp);
                const amrex::Real weight = dz * (detj_arr ? detj_arr(i,j,k) : amrex::Real(1.0));
                rho_arr(i,j,k) = weight * rho;
                water_arr(i,j,k) = weight * (cons_arr(i,j,k,RhoQ1_comp)
                                           + cons_arr(i,j,k,RhoQ2_comp)
                                           + cons_arr(i,j,k,RhoQ3_comp));
                latent_arr(i,j,k) = weight * rho * (tabs + kThetaLatentOverCp * qv);
            });
        });
    }

    return ColumnBudget{
        sum_multifab_scalar(rho_sum),
        sum_multifab_scalar(water_sum),
        sum_multifab_scalar(latent_sum),
        sum_multifab_scalar(rain_accum)};
}

inline std::vector<KernelCase> make_kernel_cases ()
{
    return {
        {amrex::Real(1.00e-2), amrex::Real(5.0e-4), amrex::Real(0.0), amrex::Real(1.16), amrex::Real(900.0),
         amrex::Real(1.00e-2), amrex::Real(0.02), amrex::Real(1.0), true,
         0, 0, 2, amrex::Real(1.16), amrex::Real(1.10), amrex::Real(2.0e-4), amrex::Real(4.0e-4),
         amrex::Real(0.0), amrex::Real(0.0), amrex::Real(1.0), amrex::Real(0.5),
         amrex::Real(0.0), amrex::Real(1.0), amrex::Real(1.0)},
        {amrex::Real(1.15e-2), amrex::Real(8.0e-4), amrex::Real(2.5e-4), amrex::Real(1.20), amrex::Real(900.0),
         amrex::Real(1.00e-2), amrex::Real(0.03), amrex::Real(1.5), true,
         1, 0, 2, amrex::Real(1.18), amrex::Real(1.20), -amrex::Real(2.0e-4), amrex::Real(3.0e-4),
         amrex::Real(4.0e-3), amrex::Real(1.0e-3), amrex::Real(0.8), amrex::Real(0.25),
         amrex::Real(0.10), amrex::Real(1.5), amrex::Real(0.5)},
        {amrex::Real(9.2e-3), amrex::Real(1.5e-3), amrex::Real(1.0e-3), amrex::Real(0.95), amrex::Real(850.0),
         amrex::Real(1.05e-2), amrex::Real(0.04), amrex::Real(2.0), true,
         3, 0, 2, amrex::Real(0.98), amrex::Real(1.02), amrex::Real(5.0e-4), amrex::Real(0.0),
         amrex::Real(2.0e-3), amrex::Real(2.0e-3), amrex::Real(1.3), amrex::Real(0.4),
         amrex::Real(0.40), amrex::Real(0.75), amrex::Real(0.25)},
        {amrex::Real(9.4e-3), amrex::Real(2.0e-4), amrex::Real(2.0e-3), amrex::Real(1.30), amrex::Real(950.0),
         amrex::Real(1.05e-2), amrex::Real(0.01), amrex::Real(0.5), false,
         2, 0, 2, amrex::Real(1.25), amrex::Real(1.30), amrex::Real(3.0e-3), amrex::Real(1.0e-3),
         amrex::Real(8.0e-3), amrex::Real(2.0e-3), amrex::Real(0.7), amrex::Real(0.125),
         amrex::Real(0.80), amrex::Real(0.5), amrex::Real(0.25)}
    };
}

inline std::string case_label (const KernelCase& test_case)
{
    return "qv=" + std::to_string(static_cast<double>(test_case.qv))
         + " qc=" + std::to_string(static_cast<double>(test_case.qc))
         + " qp=" + std::to_string(static_cast<double>(test_case.qp))
         + " rho=" + std::to_string(static_cast<double>(test_case.rho));
}

} // namespace kessler_test

#endif

#ifndef ERF_GTEST_SAM_COMMON_H_
#define ERF_GTEST_SAM_COMMON_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>

#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_MicrophysicsUtils.H>
#include <ERF_SAM.H>
#include <ERF_SAMUtils.H>

namespace sam_test {

constexpr amrex::Real kRdOcp = R_d / Cp_d;
constexpr amrex::Real kFacCond = lcond / Cp_d;
constexpr amrex::Real kFacFus = lfus / Cp_d;
constexpr amrex::Real kFacSub = lsub / Cp_d;

// MUST MATCH the SAM helper integer used in ERF_CloudSAM.cpp and ERF_Precip.cpp
// for the cloud-ice-enabled mode.
constexpr int kSAMWithIceMode = 1;

// MUST MATCH the SAM helper integer used in ERF_CloudSAM.cpp and ERF_Precip.cpp
// for the no-ice mode.
constexpr int kSAMNoIceMode = 2;

// MUST MATCH qp_threshold in ERF_Constants.H.
// This duplication is intentional because the tests sample both sides of the
// production precipitation-flux activation branch.
constexpr amrex::Real kQpThresholdForBranchSampling = qp_threshold;

// MUST MATCH the phase-fraction activation thresholds in sam_accretion_rates.
constexpr amrex::Real kPhaseActivationLow = amrex::Real(0.001);
constexpr amrex::Real kPhaseActivationHigh = amrex::Real(0.999);

struct SAMBudget {
    amrex::Real total_water_mass{amrex::Real(0.0)};
    amrex::Real latent_proxy_mass{amrex::Real(0.0)};
    amrex::Real rain_accum_mass{amrex::Real(0.0)};
    amrex::Real snow_accum_mass{amrex::Real(0.0)};
    amrex::Real graupel_accum_mass{amrex::Real(0.0)};
};

struct SAMPressureSnapshotBudget {
    amrex::Real total_water_mass{amrex::Real(0.0)};
    amrex::Real latent_proxy_mass_constant_pressure{amrex::Real(0.0)};
    amrex::Real latent_proxy_mass_eos_projected{amrex::Real(0.0)};
    amrex::Real min_pressure0_mbar{std::numeric_limits<amrex::Real>::max()};
    amrex::Real max_pressure0_mbar{std::numeric_limits<amrex::Real>::lowest()};
};

// Provenance: exact-zero branches only accumulate machine roundoff.
inline amrex::Real exact_zero_or_near_zero_tol ()
{
#ifdef AMREX_USE_FLOAT
    return amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon();
#else
    return amrex::Real(32.0) * std::numeric_limits<amrex::Real>::epsilon();
#endif
}

inline amrex::Real exact_zero_tol ()
{
    return exact_zero_or_near_zero_tol();
}

// Provenance: a handful of multiplies/divides on smooth algebraic identities.
inline amrex::Real formula_abs_tol (const amrex::Real expected)
{
#ifdef AMREX_USE_FLOAT
    const amrex::Real factor = amrex::Real(128.0);
#else
    const amrex::Real factor = amrex::Real(32.0);
#endif
    return factor * std::numeric_limits<amrex::Real>::epsilon() *
           std::max(amrex::Real(1.0), std::abs(expected));
}

inline amrex::Real roundoff_tol (const amrex::Real expected)
{
    return formula_abs_tol(expected);
}

// Provenance: pow/sqrt helpers vary slightly across host math libraries.
inline amrex::Real backend_math_abs_tol (const amrex::Real expected)
{
#ifdef AMREX_USE_FLOAT
    const amrex::Real factor = amrex::Real(4096.0);
#else
    const amrex::Real factor = amrex::Real(256.0);
#endif
    return factor * std::numeric_limits<amrex::Real>::epsilon() *
           std::max(amrex::Real(1.0), std::abs(expected));
}

inline amrex::Real pow_sqrt_tol (const amrex::Real expected)
{
    return backend_math_abs_tol(expected);
}

// Provenance: reduced sums over several flux/source terms in a single cell.
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

inline amrex::Real flux_conservation_tol (const amrex::Real scale)
{
    return property_accumulation_tol(1, scale);
}

// Provenance: MPI/global reduction order plus accumulation over many cellwise
// precip and accumulation contributions in one domain-wide budget.
inline amrex::Real mpi_reduction_tol (const amrex::Real scale,
                                      const int num_terms)
{
#ifdef AMREX_USE_FLOAT
    const amrex::Real factor = amrex::Real(512.0);
#else
    const amrex::Real factor = amrex::Real(128.0);
#endif
    return factor * std::numeric_limits<amrex::Real>::epsilon() *
           std::max(amrex::Real(1.0), std::abs(scale)) * std::max(1, num_terms);
}

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

// Provenance: Newton helper temperatures are in kelvin and the scalar formulas
// are smooth away from branch boundaries.
inline amrex::Real newton_temperature_tol (const amrex::Real expected)
{
    return roundoff_tol(expected);
}

inline amrex::Real normalized_error (const amrex::Real actual,
                                     const amrex::Real expected)
{
    return std::abs(actual - expected) /
           std::max(amrex::Real(1.0), std::abs(expected));
}

template <typename F>
inline amrex::Real central_difference (F&& func,
                                       const amrex::Real x,
                                       const amrex::Real h)
{
    return (func(x + h) - func(x - h)) / (amrex::Real(2.0) * h);
}

inline amrex::Real derivative_step (const amrex::Real scale)
{
#ifdef AMREX_USE_FLOAT
    return amrex::Real(2.0e-3) * std::max(amrex::Real(1.0), std::abs(scale));
#else
    return amrex::Real(1.0e-6) * std::max(amrex::Real(1.0), std::abs(scale));
#endif
}

inline amrex::Box single_cell_box ()
{
    return amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
}

inline amrex::Array4<amrex::Real> single_cell_state_array (amrex::FArrayBox& fab)
{
    return fab.array();
}

inline std::array<amrex::Real, 5> threshold_samples (const amrex::Real lower,
                                                      const amrex::Real upper)
{
    return {
        lower - amrex::Real(1.0e-3),
        lower,
        amrex::Real(0.5) * (lower + upper),
        upper,
        upper + amrex::Real(1.0e-3)};
}

inline std::string branch_label (const char* helper_name,
                                 const amrex::Real sample)
{
    return std::string(helper_name) + " sample=" +
           std::to_string(static_cast<double>(sample));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_total_water (const SAMCellState& s) noexcept
{
    return s.qv + s.qcl + s.qci + s.qpr + s.qps + s.qpg;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_latent_proxy (const SAMCellState& s,
                              amrex::Real fac_cond,
                              amrex::Real fac_fus) noexcept
{
    return s.tabs + fac_cond * s.qv - fac_fus * (s.qci + s.qps + s.qpg);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sam_latent_proxy_constant_pressure (const amrex::Real theta,
                                                const amrex::Real pressure_mbar,
                                                const amrex::Real qv,
                                                const amrex::Real qci,
                                                const amrex::Real qps,
                                                const amrex::Real qpg,
                                                const amrex::Real rdOcp,
                                                const amrex::Real fac_cond,
                                                const amrex::Real fac_fus) noexcept
{
    const amrex::Real tabs = getTgivenPandTh(sam_mbar_to_pa(pressure_mbar), theta, rdOcp);
    return tabs + fac_cond * qv - fac_fus * (qci + qps + qpg);
}

inline void fill_pressure_snapshot_from_cons (const amrex::MultiFab& cons,
                                             amrex::MultiFab& pressure0_mbar)
{
    auto cons_arrays = cons.const_arrays();
    auto pressure_arrays = pressure0_mbar.arrays();
    amrex::ParallelFor(
        pressure0_mbar,
        amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
            const auto& arr = cons_arrays[box_no];
            pressure_arrays[box_no](i,j,k) =
                sam_cons_to_primitive(arr(i,j,k,Rho_comp),
                                      arr(i,j,k,RhoTheta_comp),
                                      arr(i,j,k,RhoQ1_comp),
                                      arr(i,j,k,RhoQ2_comp),
                                      arr(i,j,k,RhoQ3_comp),
                                      arr(i,j,k,RhoQ4_comp),
                                      arr(i,j,k,RhoQ5_comp),
                                      arr(i,j,k,RhoQ6_comp)).pres_mbar;
        });
    amrex::Gpu::streamSynchronize();
}

inline SAMPressureSnapshotBudget compute_sam_pressure_snapshot_budget (
    const amrex::Geometry& geom,
    const amrex::MultiFab& cons,
    const amrex::MultiFab& pressure0_mbar,
    const amrex::MultiFab* detJ_cc = nullptr)
{
    const amrex::Real cell_volume = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
    auto cons_arrays = cons.const_arrays();
    auto pressure_arrays = pressure0_mbar.const_arrays();

    amrex::Real total_water_mass = amrex::Real(0.0);
    amrex::Real latent_proxy_mass_constant_pressure = amrex::Real(0.0);
    amrex::Real latent_proxy_mass_eos_projected = amrex::Real(0.0);
    amrex::Real min_pressure0_mbar = std::numeric_limits<amrex::Real>::max();
    amrex::Real max_pressure0_mbar = std::numeric_limits<amrex::Real>::lowest();

    if (detJ_cc != nullptr) {
        auto detj_arrays = detJ_cc->const_arrays();
        amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real> reduced = amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpSum,
                            amrex::ReduceOpSum,
                            amrex::ReduceOpSum,
                            amrex::ReduceOpMin,
                            amrex::ReduceOpMax>{},
            amrex::TypeList<amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real>{},
            cons, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                -> amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real> {
                const auto& arr = cons_arrays[box_no];
                const SAMPrimitiveCell primitive =
                    sam_cons_to_primitive(arr(i,j,k,Rho_comp),
                                          arr(i,j,k,RhoTheta_comp),
                                          arr(i,j,k,RhoQ1_comp),
                                          arr(i,j,k,RhoQ2_comp),
                                          arr(i,j,k,RhoQ3_comp),
                                          arr(i,j,k,RhoQ4_comp),
                                          arr(i,j,k,RhoQ5_comp),
                                          arr(i,j,k,RhoQ6_comp));
                const amrex::Real pressure0 = pressure_arrays[box_no](i,j,k);
                const amrex::Real weight = cell_volume * detj_arrays[box_no](i,j,k) * primitive.rho;
                return {
                    weight * sam_total_water(primitive),
                    weight * sam_latent_proxy_constant_pressure(primitive.theta,
                                                                pressure0,
                                                                primitive.qv,
                                                                primitive.qci,
                                                                primitive.qps,
                                                                primitive.qpg,
                                                                kRdOcp,
                                                                kFacCond,
                                                                kFacFus),
                    weight * sam_latent_proxy(primitive, kFacCond, kFacFus),
                    pressure0,
                    pressure0};
            });
        total_water_mass = amrex::get<0>(reduced);
        latent_proxy_mass_constant_pressure = amrex::get<1>(reduced);
        latent_proxy_mass_eos_projected = amrex::get<2>(reduced);
        min_pressure0_mbar = amrex::get<3>(reduced);
        max_pressure0_mbar = amrex::get<4>(reduced);
    } else {
        amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real> reduced = amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpSum,
                            amrex::ReduceOpSum,
                            amrex::ReduceOpSum,
                            amrex::ReduceOpMin,
                            amrex::ReduceOpMax>{},
            amrex::TypeList<amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real>{},
            cons, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                -> amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real> {
                const auto& arr = cons_arrays[box_no];
                const SAMPrimitiveCell primitive =
                    sam_cons_to_primitive(arr(i,j,k,Rho_comp),
                                          arr(i,j,k,RhoTheta_comp),
                                          arr(i,j,k,RhoQ1_comp),
                                          arr(i,j,k,RhoQ2_comp),
                                          arr(i,j,k,RhoQ3_comp),
                                          arr(i,j,k,RhoQ4_comp),
                                          arr(i,j,k,RhoQ5_comp),
                                          arr(i,j,k,RhoQ6_comp));
                const amrex::Real pressure0 = pressure_arrays[box_no](i,j,k);
                const amrex::Real weight = cell_volume * primitive.rho;
                return {
                    weight * sam_total_water(primitive),
                    weight * sam_latent_proxy_constant_pressure(primitive.theta,
                                                                pressure0,
                                                                primitive.qv,
                                                                primitive.qci,
                                                                primitive.qps,
                                                                primitive.qpg,
                                                                kRdOcp,
                                                                kFacCond,
                                                                kFacFus),
                    weight * sam_latent_proxy(primitive, kFacCond, kFacFus),
                    pressure0,
                    pressure0};
            });
        total_water_mass = amrex::get<0>(reduced);
        latent_proxy_mass_constant_pressure = amrex::get<1>(reduced);
        latent_proxy_mass_eos_projected = amrex::get<2>(reduced);
        min_pressure0_mbar = amrex::get<3>(reduced);
        max_pressure0_mbar = amrex::get<4>(reduced);
    }

    amrex::ParallelDescriptor::ReduceRealSum(total_water_mass);
    amrex::ParallelDescriptor::ReduceRealSum(latent_proxy_mass_constant_pressure);
    amrex::ParallelDescriptor::ReduceRealSum(latent_proxy_mass_eos_projected);
    amrex::ParallelDescriptor::ReduceRealMin(min_pressure0_mbar);
    amrex::ParallelDescriptor::ReduceRealMax(max_pressure0_mbar);

    return SAMPressureSnapshotBudget{
        total_water_mass,
        latent_proxy_mass_constant_pressure,
        latent_proxy_mass_eos_projected,
        min_pressure0_mbar,
        max_pressure0_mbar};
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

inline SAMBudget compute_sam_budget (const amrex::Geometry& geom,
                                     const amrex::MultiFab& cons,
                                     const SAM& sam,
                                     const amrex::MultiFab* detJ_cc = nullptr)
{
    const amrex::Real cell_volume = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
    const amrex::Real cell_area = geom.CellSize(0) * geom.CellSize(1);
    auto cons_arrays = cons.const_arrays();

    amrex::Real total_water_mass = amrex::Real(0.0);
    amrex::Real latent_proxy_mass = amrex::Real(0.0);

    if (detJ_cc != nullptr) {
        auto detj_arrays = detJ_cc->const_arrays();
        amrex::GpuTuple<amrex::Real, amrex::Real> reduced = amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpSum, amrex::ReduceOpSum>{},
            amrex::TypeList<amrex::Real, amrex::Real>{},
            cons, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Real, amrex::Real> {
                const auto& arr = cons_arrays[box_no];
                const SAMPrimitiveCell primitive =
                    sam_cons_to_primitive(arr(i,j,k,Rho_comp),
                                          arr(i,j,k,RhoTheta_comp),
                                          arr(i,j,k,RhoQ1_comp),
                                          arr(i,j,k,RhoQ2_comp),
                                          arr(i,j,k,RhoQ3_comp),
                                          arr(i,j,k,RhoQ4_comp),
                                          arr(i,j,k,RhoQ5_comp),
                                          arr(i,j,k,RhoQ6_comp));
                const amrex::Real weight = cell_volume * detj_arrays[box_no](i,j,k);
                return {weight * primitive.rho * sam_total_water(primitive),
                        weight * primitive.rho * sam_latent_proxy(primitive, kFacCond, kFacFus)};
            });
        total_water_mass = amrex::get<0>(reduced);
        latent_proxy_mass = amrex::get<1>(reduced);
    } else {
        amrex::GpuTuple<amrex::Real, amrex::Real> reduced = amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpSum, amrex::ReduceOpSum>{},
            amrex::TypeList<amrex::Real, amrex::Real>{},
            cons, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Real, amrex::Real> {
                const auto& arr = cons_arrays[box_no];
                const SAMPrimitiveCell primitive =
                    sam_cons_to_primitive(arr(i,j,k,Rho_comp),
                                          arr(i,j,k,RhoTheta_comp),
                                          arr(i,j,k,RhoQ1_comp),
                                          arr(i,j,k,RhoQ2_comp),
                                          arr(i,j,k,RhoQ3_comp),
                                          arr(i,j,k,RhoQ4_comp),
                                          arr(i,j,k,RhoQ5_comp),
                                          arr(i,j,k,RhoQ6_comp));
                return {cell_volume * primitive.rho * sam_total_water(primitive),
                        cell_volume * primitive.rho * sam_latent_proxy(primitive, kFacCond, kFacFus)};
            });
        total_water_mass = amrex::get<0>(reduced);
        latent_proxy_mass = amrex::get<1>(reduced);
    }

    return SAMBudget{
        total_water_mass,
        latent_proxy_mass,
        sum_multifab_scalar(*sam.Qmoist_Ptr(0)) * cell_area * rhor / amrex::Real(1000.0),
        sum_multifab_scalar(*sam.Qmoist_Ptr(1)) * cell_area * rhos / amrex::Real(1000.0),
        sum_multifab_scalar(*sam.Qmoist_Ptr(2)) * cell_area * rhog / amrex::Real(1000.0)};
}

} // namespace sam_test

#endif
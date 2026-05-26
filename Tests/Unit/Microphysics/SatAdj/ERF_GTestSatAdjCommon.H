#ifndef ERF_GTEST_SATADJ_COMMON_H_
#define ERF_GTEST_SATADJ_COMMON_H_

#include <algorithm>
#include <cmath>
#include <utility>

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>

#include <ERF_Constants.H>
#include <ERF_DataStruct.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_MicrophysicsUtils.H>
#include <ERF_SatAdj.H>

namespace satadj_test {

// SatAdj unit-test contract
// -------------------------
// These tests check physical and numerical properties of ERF's warm-cloud
// saturation adjustment, not hard-coded gold values.
//
// Core invariants:
//   1. qt = qv + qc is conserved locally.
//   2. At fixed pressure, T + (L/cp) qv is conserved through the phase-change
//      update used by AdjustSatAdjCell.
//   3. Final scalar states satisfy complementarity:
//        qc >= 0,
//        qv <= qsat(T,p),
//        qc > 0 implies qv ~= qsat(T,p).
//
// Thermodynamics:
//   - Use ERF thermodynamic utilities only.
//   - Do not add an independent saturation formula to the tests.
//   - SatAdj qsat functions take pressure in mbar.
//   - ERF EOS temperature/theta functions take pressure in Pa.
//
// Portability:
//   - Production SatAdj runs through AMReX ParallelFor.
//   - MultiFab/kernel tests should use AMReX-portable setup and checks.
//   - Keep GTest EXPECT/ASSERT calls on the host. Device kernels should only
//     compute adjusted states or reduced error metrics.

#ifdef AMREX_USE_FLOAT
constexpr amrex::Real kThermoTolFactor = amrex::Real(2.0e-5);
constexpr amrex::Real kStateTolFactor = amrex::Real(5.0e-5);
constexpr amrex::Real kAsymptoticRelTol = amrex::Real(2.0e-1);
constexpr amrex::Real kNormalizedWaterConservationTol = amrex::Real(8.0);
constexpr amrex::Real kNormalizedLatentEnergyConservationTol = amrex::Real(8.0);
#else
constexpr amrex::Real kThermoTolFactor = amrex::Real(1.0e-12);
constexpr amrex::Real kStateTolFactor = amrex::Real(1.0e-12);
constexpr amrex::Real kAsymptoticRelTol = amrex::Real(2.0e-2);
constexpr amrex::Real kNormalizedWaterConservationTol = amrex::Real(1.0);
constexpr amrex::Real kNormalizedLatentEnergyConservationTol = amrex::Real(1.0);
#endif

constexpr amrex::Real kFacCond = lcond / Cp_d;
constexpr amrex::Real kRdOcp = R_d / Cp_d;

struct CellState {
    amrex::Real tabs;
    amrex::Real pres_mbar;
    amrex::Real theta;
    amrex::Real qv;
    amrex::Real qc;
};

struct ConservedState {
    amrex::Real rho;
    amrex::Real rhotheta;
    amrex::Real rhoqv;
    amrex::Real rhoqc;
};

enum StateDiffComp {
    StateDiffRho = 0,
    StateDiffRhoTheta,
    StateDiffQ1,
    StateDiffQ2,
    NumStateDiffComps
};

enum InvariantErrComp {
    InvariantErrRho = 0,
    InvariantErrQt,
    InvariantErrQcNegative,
    InvariantErrRhoThetaRef,
    InvariantErrQ1Ref,
    InvariantErrQ2Ref,
    NumInvariantErrComps
};

enum ConservationErrComp {
    ConservationErrRho = 0,
    ConservationErrQt,
    ConservationErrLatentEnergy,
    ConservationErrInitialQcPositive,
    ConservationErrFinalQcPositive,
    ConservationErrFinalQcZero,
    NumConservationErrComps
};

#ifdef AMREX_USE_FLOAT
constexpr amrex::Real kQcCoveragePositiveMin = amrex::Real(1.0e-6);
constexpr amrex::Real kQcCoverageNearZeroMax = amrex::Real(1.0e-6);
#else
constexpr amrex::Real kQcCoveragePositiveMin = amrex::Real(1.0e-12);
constexpr amrex::Real kQcCoverageNearZeroMax = amrex::Real(1.0e-12);
#endif

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real scaled_tol (const amrex::Real reference, const amrex::Real factor)
{
    const amrex::Real magnitude = amrex::Math::abs(reference);
    return factor * (magnitude > amrex::Real(1.0) ? magnitude : amrex::Real(1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real water_budget_tol (const amrex::Real reference)
{
#ifdef AMREX_USE_FLOAT
    return amrex::max(amrex::Real(1.0e-8),
                      amrex::Real(1.0e-5) * amrex::Math::abs(reference));
#else
    return amrex::max(amrex::Real(1.0e-14),
                      amrex::Real(1.0e-12) * amrex::Math::abs(reference));
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real latent_energy_budget_tol (const amrex::Real reference)
{
#ifdef AMREX_USE_FLOAT
    return amrex::max(amrex::Real(1.0e-4),
                      amrex::Real(1.0e-4) * amrex::Math::abs(reference));
#else
    return amrex::max(amrex::Real(1.0e-10),
                      amrex::Real(1.0e-10) * amrex::Math::abs(reference));
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real qsat (const amrex::Real tabs, const amrex::Real pres_mbar)
{
    amrex::Real qsat_local;
    erf_qsatw(tabs, pres_mbar, qsat_local);
    return qsat_local;
}

inline amrex::Real dqsat_dt (const amrex::Real tabs, const amrex::Real pres_mbar)
{
    amrex::Real dqsat_local;
    erf_dtqsatw(tabs, pres_mbar, dqsat_local);
    return dqsat_local;
}

inline CellState make_cell_state (const amrex::Real tabs,
                                  const amrex::Real pres_mbar,
                                  const amrex::Real qv,
                                  const amrex::Real qc)
{
    return CellState{tabs,
                     pres_mbar,
                     getThgivenTandP(tabs, amrex::Real(100.0) * pres_mbar, kRdOcp),
                     qv,
                     qc};
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void adjust (CellState& state)
{
    SatAdj::AdjustSatAdjCell(kFacCond,
                             kRdOcp,
                             state.tabs,
                             state.pres_mbar,
                             state.theta,
                             state.qv,
                             state.qc);
}

inline ConservedState make_conserved_state (const CellState& state)
{
    const amrex::Real pres_pa = amrex::Real(100.0) * state.pres_mbar;
    const amrex::Real rho = getRhogivenTandPress(state.tabs, pres_pa, state.qv);
    const amrex::Real theta = getThgivenTandP(state.tabs, pres_pa, kRdOcp);

    return ConservedState{rho, rho * theta, rho * state.qv, rho * state.qc};
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
CellState make_state_from_conserved (const amrex::Real rho,
                                     const amrex::Real rhotheta,
                                     const amrex::Real rhoqv,
                                     const amrex::Real rhoqc)
{
    const amrex::Real qv = rhoqv / rho;
    return CellState{getTgivenRandRTh(rho, rhotheta, qv),
                     getPgivenRTh(rhotheta, qv) * amrex::Real(0.01),
                     rhotheta / rho,
                     qv,
                     rhoqc / rho};
}

inline SolverChoice make_solver_choice (const bool use_shoc)
{
    SolverChoice sc;
    sc.c_p = Cp_d;
    sc.rdOcp = kRdOcp;
    sc.use_shoc = use_shoc;
    return sc;
}

// Construct a state that enters the subtle evaporation-then-recondensation
// path. Initially qt < qsat(T_i,p), so SatAdj evaporates all qc. The
// evaporative cooling can lower qsat enough that the cell becomes
// supersaturated and must re-enter the Newton solve.
inline bool find_evaporation_then_recondensation_state (CellState& state)
{
    const amrex::Real tabs = amrex::Real(285.0);
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qsat_initial = qsat(tabs, pres_mbar);
    const amrex::Real qv = amrex::Real(0.8) * qsat_initial;

    for (int iqc = 1; iqc <= 80; ++iqc) {
        const amrex::Real qc = amrex::Real(iqc) * amrex::Real(1.0e-4);
        if ((qv + qc) >= qsat_initial) {
            break;
        }

        const amrex::Real tabs_all_evap = tabs - kFacCond * qc;
        if ((qv + qc) > qsat(tabs_all_evap, pres_mbar)) {
            state = make_cell_state(tabs, pres_mbar, qv, qc);
            return true;
        }
    }

    return false;
}

inline void sync ()
{
#ifdef AMREX_USE_GPU
    amrex::Gpu::streamSynchronize();
#endif
}

template <typename F>
inline void run_and_sync (F&& f)
{
    std::forward<F>(f)();
    sync();
}

// Fill a small conserved-state MultiFab with deterministic SatAdj branch
// cases. This uses ParallelFor so the same setup path works in CPU and GPU
// builds.
inline void fill_conserved_state_portable (amrex::MultiFab& cons,
                                           const bool include_all_satadj_branches)
{
    CellState evap_then_recond;
    const bool found = find_evaporation_then_recondensation_state(evap_then_recond);
    AMREX_ALWAYS_ASSERT(found);

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto arr = cons.array(mfi);

        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const int selector = include_all_satadj_branches
                    ? (i + 2 * j + 3 * k) % 4
                    : (i + 2 * j + 3 * k) % 2;

                amrex::Real tabs;
                amrex::Real pres_mbar;
                amrex::Real qv;
                amrex::Real qc;

                if (selector == 0) {
                    tabs = amrex::Real(290.0);
                    pres_mbar = amrex::Real(900.0);
                    qv = qsat(tabs, pres_mbar) + amrex::Real(6.0e-4);
                    qc = amrex::Real(3.0e-4);
                } else if (selector == 1) {
                    tabs = amrex::Real(290.0);
                    pres_mbar = amrex::Real(900.0);
                    qv = amrex::Real(4.0e-3);
                    qc = amrex::Real(1.0e-3);
                } else if (selector == 2) {
                    tabs = evap_then_recond.tabs;
                    pres_mbar = evap_then_recond.pres_mbar;
                    qv = evap_then_recond.qv;
                    qc = evap_then_recond.qc;
                } else {
                    tabs = amrex::Real(288.0);
                    pres_mbar = amrex::Real(850.0);
                    qv = qsat(tabs, pres_mbar) + amrex::Real(4.0e-4);
                    qc = -amrex::Real(1.0e-4);
                }

                const amrex::Real pres_pa = amrex::Real(100.0) * pres_mbar;
                const amrex::Real rho = getRhogivenTandPress(tabs, pres_pa, qv);
                const amrex::Real theta = getThgivenTandP(tabs, pres_pa, kRdOcp);

                arr(i,j,k,Rho_comp) = rho;
                arr(i,j,k,RhoTheta_comp) = rho * theta;
                arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
                arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
                arr(i,j,k,RhoQ1_comp) = rho * qv;
                arr(i,j,k,RhoQ2_comp) = rho * qc;
            });
        });
    }
}

inline void fill_conserved_state_portable (amrex::MultiFab& cons)
{
    fill_conserved_state_portable(cons, true);
}

// The public-flow conservation contract is defined for physical SatAdj states
// with qv and qc only. Keep the dedicated negative-qc repair coverage in the
// scalar/reference tests, and use valid cloud-water states here so the new
// MultiFab conservation test measures only phase-change conservation.
inline void fill_conserved_state_conservation_portable (amrex::MultiFab& cons)
{
    CellState evap_then_recond;
    const bool found = find_evaporation_then_recondensation_state(evap_then_recond);
    AMREX_ALWAYS_ASSERT(found);

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto arr = cons.array(mfi);

        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const int selector = (i + 2 * j + 3 * k) % 4;

                amrex::Real tabs;
                amrex::Real pres_mbar;
                amrex::Real qv;
                amrex::Real qc;

                if (selector == 0) {
                    tabs = amrex::Real(290.0);
                    pres_mbar = amrex::Real(900.0);
                    qv = qsat(tabs, pres_mbar) + amrex::Real(6.0e-4);
                    qc = amrex::Real(3.0e-4);
                } else if (selector == 1) {
                    tabs = amrex::Real(290.0);
                    pres_mbar = amrex::Real(900.0);
                    qv = amrex::Real(4.0e-3);
                    qc = amrex::Real(1.0e-3);
                } else if (selector == 2) {
                    tabs = evap_then_recond.tabs;
                    pres_mbar = evap_then_recond.pres_mbar;
                    qv = evap_then_recond.qv;
                    qc = evap_then_recond.qc;
                } else {
                    tabs = amrex::Real(285.0);
                    pres_mbar = amrex::Real(880.0);
                    qv = qsat(tabs, pres_mbar) + amrex::Real(2.0e-4);
                    qc = amrex::Real(5.0e-5);
                }

                const amrex::Real pres_pa = amrex::Real(100.0) * pres_mbar;
                const amrex::Real rho = getRhogivenTandPress(tabs, pres_pa, qv);
                const amrex::Real theta = getThgivenTandP(tabs, pres_pa, kRdOcp);

                arr(i,j,k,Rho_comp) = rho;
                arr(i,j,k,RhoTheta_comp) = rho * theta;
                arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
                arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
                arr(i,j,k,RhoQ1_comp) = rho * qv;
                arr(i,j,k,RhoQ2_comp) = rho * qc;
            });
        });
    }
}

inline void compute_state_difference (const amrex::MultiFab& after,
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
            });
        });
    }
}

inline void compute_satadj_invariant_errors (const amrex::MultiFab& before,
                                             const amrex::MultiFab& after,
                                             amrex::MultiFab& err)
{
    // Compare the public MultiFab flow against the scalar SatAdj helper and
    // core conservation properties. The scalar reference is built from the
    // pre-state because the cell solve holds pressure fixed, while post-state
    // conserved variables do not reconstruct that fixed-pressure state
    // uniquely. Errors are normalized by each tolerance so host-side tests can
    // check small O(1) bounds without duplicating tolerance logic.
    for (amrex::MFIter mfi(after, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto before_arr = before.const_array(mfi);
        auto after_arr = after.const_array(mfi);
        auto err_arr = err.array(mfi);

        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const amrex::Real rho_before = before_arr(i,j,k,Rho_comp);
                const amrex::Real rho_after = after_arr(i,j,k,Rho_comp);
                const amrex::Real qt_before = before_arr(i,j,k,RhoQ1_comp) + before_arr(i,j,k,RhoQ2_comp);
                const amrex::Real qt_after = after_arr(i,j,k,RhoQ1_comp) + after_arr(i,j,k,RhoQ2_comp);
                const amrex::Real qc_after = after_arr(i,j,k,RhoQ2_comp) / rho_after;
                CellState reference_state = make_state_from_conserved(rho_before,
                                                                      before_arr(i,j,k,RhoTheta_comp),
                                                                      before_arr(i,j,k,RhoQ1_comp),
                                                                      before_arr(i,j,k,RhoQ2_comp));
                adjust(reference_state);

                const amrex::Real rhotheta_ref = rho_before * reference_state.theta;
                const amrex::Real rhoq1_ref = rho_before * reference_state.qv;
                const amrex::Real rhoq2_ref = rho_before * reference_state.qc;

                err_arr(i,j,k,InvariantErrRho) = amrex::Math::abs(rho_after - rho_before)
                    / scaled_tol(rho_before, kStateTolFactor);
                err_arr(i,j,k,InvariantErrQt) = amrex::Math::abs(qt_after - qt_before)
                    / scaled_tol(qt_before, amrex::Real(20.0) * kStateTolFactor);
                err_arr(i,j,k,InvariantErrQcNegative) = qc_after < amrex::Real(0.0)
                    ? (-qc_after) / scaled_tol(amrex::Real(0.0), amrex::Real(20.0) * kStateTolFactor)
                    : amrex::Real(0.0);
                err_arr(i,j,k,InvariantErrRhoThetaRef) = amrex::Math::abs(after_arr(i,j,k,RhoTheta_comp)
                                                                          - rhotheta_ref)
                    / scaled_tol(rhotheta_ref, amrex::Real(20.0) * kStateTolFactor);
                err_arr(i,j,k,InvariantErrQ1Ref) = amrex::Math::abs(after_arr(i,j,k,RhoQ1_comp) - rhoq1_ref)
                    / scaled_tol(rhoq1_ref, amrex::Real(20.0) * kStateTolFactor);
                err_arr(i,j,k,InvariantErrQ2Ref) = amrex::Math::abs(after_arr(i,j,k,RhoQ2_comp) - rhoq2_ref)
                    / scaled_tol(rhoq2_ref, amrex::Real(20.0) * kStateTolFactor);
            });
        });
    }
}

inline void compute_satadj_conservation_errors (const amrex::MultiFab& before,
                                                const amrex::MultiFab& after,
                                                amrex::MultiFab& err)
{
    // This helper checks the public SatAdj MultiFab flow as a conservation
    // contract, not against a scalar gold implementation. SatAdj transports
    // only non-precipitating water (qv and qc), so the conserved quantity here
    // is rho * (qv + qc), not any precipitating-water budget.
    //
    // The cell solve itself updates qv and qc with algebraically conservative
    // assignments. Any nonzero public-flow rho/qt residual observed here comes
    // from floating-point roundoff in the conserved-to-primitive and
    // primitive-to-conserved rho*q conversions around the SatAdj call, not
    // from an intended source/sink in the scheme.

    // SatAdj holds the diagnosed cell pressure fixed during each adjustment.
    // Evaluate both the pre-state and post-state latent-energy invariant along
    // that same fixed-pressure path:
    //
    // Condensation/recondensation:
    //   T1 = T0 + A * (qv0 - qv1)
    //   qc1 = qc0 + qv0 - qv1
    //   A = L / cp
    //
    //   qv1 + qc1 - qv0 - qc0 = 0
    //   T1 + A*qv1 - T0 - A*qv0 = 0
    //
    // All evaporation:
    //   qv1 = qv0 + qc0
    //   qc1 = 0
    //   T1 = T0 - A*qc0
    //
    //   qv1 + qc1 - qv0 - qc0 = 0
    //   T1 + A*qv1 - T0 - A*qv0 = 0
    for (amrex::MFIter mfi(after, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto before_arr = before.const_array(mfi);
        auto after_arr = after.const_array(mfi);
        auto err_arr = err.array(mfi);

        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const amrex::Real rho_before = before_arr(i,j,k,Rho_comp);
                const amrex::Real rho_after = after_arr(i,j,k,Rho_comp);
                const amrex::Real rhotheta_before = before_arr(i,j,k,RhoTheta_comp);
                const amrex::Real rhotheta_after = after_arr(i,j,k,RhoTheta_comp);
                const amrex::Real rhoqv_before = before_arr(i,j,k,RhoQ1_comp);
                const amrex::Real rhoqv_after = after_arr(i,j,k,RhoQ1_comp);
                const amrex::Real rhoqc_before = before_arr(i,j,k,RhoQ2_comp);
                const amrex::Real rhoqc_after = after_arr(i,j,k,RhoQ2_comp);

                const amrex::Real theta_before = rhotheta_before / rho_before;
                const amrex::Real theta_after = rhotheta_after / rho_after;
                const amrex::Real qv_before = rhoqv_before / rho_before;
                const amrex::Real qv_after = rhoqv_after / rho_after;
                const amrex::Real qc_before = rhoqc_before / rho_before;
                const amrex::Real qc_after = rhoqc_after / rho_after;

                const amrex::Real qt_before = rhoqv_before + rhoqc_before;
                const amrex::Real qt_after = rhoqv_after + rhoqc_after;

                const amrex::Real p_before = getPgivenRTh(rhotheta_before, qv_before);
                const amrex::Real t_before = getTgivenPandTh(p_before, theta_before, kRdOcp);
                const amrex::Real t_after_fixed_p = getTgivenPandTh(p_before, theta_after, kRdOcp);

                const amrex::Real latent_before = rho_before * (t_before + kFacCond * qv_before);
                const amrex::Real latent_after = rho_after * (t_after_fixed_p + kFacCond * qv_after);

                err_arr(i,j,k,ConservationErrRho) = amrex::Math::abs(rho_after - rho_before)
                    / scaled_tol(rho_before, kStateTolFactor);
                err_arr(i,j,k,ConservationErrQt) = amrex::Math::abs(qt_after - qt_before)
                    / water_budget_tol(qt_before);
                // Reconstruct temperature at the pre-adjustment fixed pressure
                // and normalize against a physically meaningful latent-energy
                // budget tolerance, not a generic scaled_tol floor of 1.
                err_arr(i,j,k,ConservationErrLatentEnergy) = amrex::Math::abs(latent_after - latent_before)
                    / latent_energy_budget_tol(latent_before);

                err_arr(i,j,k,ConservationErrInitialQcPositive) =
                    qc_before > kQcCoveragePositiveMin ? amrex::Real(1.0) : amrex::Real(0.0);
                err_arr(i,j,k,ConservationErrFinalQcPositive) =
                    qc_after > kQcCoveragePositiveMin ? amrex::Real(1.0) : amrex::Real(0.0);
                err_arr(i,j,k,ConservationErrFinalQcZero) =
                    (qc_before > kQcCoveragePositiveMin &&
                     amrex::Math::abs(qc_after) <= kQcCoverageNearZeroMax)
                    ? amrex::Real(1.0)
                    : amrex::Real(0.0);
            });
        });
    }
}

inline amrex::Geometry make_geometry (const int nx, const int ny, const int nz)
{
    amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
                      amrex::IntVect(AMREX_D_DECL(nx - 1, ny - 1, nz - 1)));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(1.0, 1.0, 1.0)});
    amrex::Array<int, AMREX_SPACEDIM> periodicity{AMREX_D_DECL(1, 1, 1)};
    return amrex::Geometry(domain, &real_box, amrex::CoordSys::cartesian, periodicity.data());
}

} // namespace satadj_test

#endif
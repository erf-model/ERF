#ifndef ERF_GTEST_SATADJ_COMMON_H_
#define ERF_GTEST_SATADJ_COMMON_H_

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

#include <AMReX_Arena.H>
#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
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

struct LocatedError {
    amrex::Real normalized_error{amrex::Real(0.0)};
    amrex::Real absolute_error{amrex::Real(0.0)};
    amrex::Real actual_value{amrex::Real(0.0)};
    amrex::Real expected_value{amrex::Real(0.0)};
    amrex::Real eos_projected_actual_temperature{amrex::Real(0.0)};
    amrex::Real eos_projected_expected_temperature{amrex::Real(0.0)};
    int i{0};
    int j{0};
    int k{0};
    int component{RhoTheta_comp};
    int box_id{-1};
};

struct ActiveCellCounts {
    long checked{0};
    long changed{0};
    long unexpected_unchanged{0};
};

enum class SatAdjActiveGeneratorRegime {
    SupersaturatedCloud = 0,
    CloudyRelaxation,
    EvaporationThenRecondensation
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

// Water mixing ratios are O(1e-3) to O(1e-2), so the generic scaled_tol(...)
// helper with an O(1) absolute floor is too permissive for qv/qc checks.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real mixing_ratio_tol (const amrex::Real reference)
{
#ifdef AMREX_USE_FLOAT
    return amrex::max(amrex::Real(2.0e-7),
                      amrex::Real(2.0e-5) * amrex::Math::abs(reference));
#else
    return amrex::max(amrex::Real(1.0e-13),
                      amrex::Real(1.0e-12) * amrex::Math::abs(reference));
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real conserved_moisture_tol (const amrex::Real rho,
                                    const amrex::Real q_reference)
{
    return amrex::max(amrex::Real(1.0e-14),
                      amrex::Math::abs(rho) * mixing_ratio_tol(q_reference));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real derived_temperature_tol (const amrex::Real temperature)
{
#ifdef AMREX_USE_FLOAT
    return amrex::max(amrex::Real(2.0e-4),
                      amrex::Real(5.0e-6) * amrex::Math::abs(temperature));
#else
    return amrex::max(amrex::Real(1.0e-10),
                      amrex::Real(1.0e-12) * amrex::Math::abs(temperature));
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

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
CellState make_cell_state (const amrex::Real tabs,
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

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
ConservedState make_conserved_state (const CellState& state)
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

inline amrex::Geometry make_multibox_geometry (
    const int nx,
    const int ny,
    const int nz,
    const amrex::Array<int, AMREX_SPACEDIM>& periodicity)
{
    amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
                      amrex::IntVect(AMREX_D_DECL(nx - 1, ny - 1, nz - 1)));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(1.0, 1.0, 1.0)});
    return amrex::Geometry(domain, &real_box, amrex::CoordSys::cartesian, periodicity.data());
}

template <typename F>
inline void run_and_sync (F&& f)
{
    std::forward<F>(f)();
    sync();
}

inline amrex::BoxArray make_boxarray (const amrex::Box& domain)
{
    return amrex::BoxArray(domain);
}

inline amrex::BoxArray make_boxarray (const amrex::Box& domain,
                                      const amrex::IntVect& max_size)
{
    amrex::BoxArray ba(domain);
    ba.maxSize(max_size);
    return ba;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
CellState make_satadj_active_cell_state (const int i,
                                         const int j,
                                         const int k,
                                         const CellState& evap_then_recond)
{
    const SatAdjActiveGeneratorRegime regime =
        static_cast<SatAdjActiveGeneratorRegime>((i + 2 * j + 3 * k) % 3);

    if (regime == SatAdjActiveGeneratorRegime::SupersaturatedCloud) {
        const amrex::Real tabs = amrex::Real(289.0) + amrex::Real(0.15) * k;
        const amrex::Real pres_mbar = amrex::Real(910.0) - amrex::Real(2.0) * j;
        return make_cell_state(tabs,
                               pres_mbar,
                               qsat(tabs, pres_mbar) + amrex::Real(8.0e-4) + amrex::Real(2.0e-5) * i,
                               amrex::Real(4.0e-4) + amrex::Real(1.0e-5) * (j + k));
    }

    if (regime == SatAdjActiveGeneratorRegime::CloudyRelaxation) {
        const amrex::Real tabs = amrex::Real(287.5) + amrex::Real(0.1) * i;
        const amrex::Real pres_mbar = amrex::Real(885.0) + amrex::Real(1.5) * k;
        return make_cell_state(tabs,
                               pres_mbar,
                               amrex::Real(3.5e-3) + amrex::Real(1.0e-4) * (i + k),
                               amrex::Real(1.4e-3) + amrex::Real(5.0e-5) * j);
    }

    return make_cell_state(evap_then_recond.tabs,
                           evap_then_recond.pres_mbar,
                           evap_then_recond.qv + amrex::Real(2.0e-5) * i,
                           evap_then_recond.qc + amrex::Real(1.0e-5) * (j + k));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
SatAdjActiveGeneratorRegime satadj_active_generator_regime (const int i,
                                                            const int j,
                                                            const int k)
{
    return static_cast<SatAdjActiveGeneratorRegime>((i + 2 * j + 3 * k) % 3);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
CellState make_satadj_kessler_active_cell_state (const int i,
                                                 const int j,
                                                 const int k)
{
    const int selector = (i + 2 * j + 3 * k) % 2;

    if (selector == 0) {
        const amrex::Real tabs = amrex::Real(289.0) + amrex::Real(0.1) * k;
        const amrex::Real pres_mbar = amrex::Real(905.0) - amrex::Real(2.0) * j;
        return make_cell_state(tabs,
                               pres_mbar,
                               qsat(tabs, pres_mbar) + amrex::Real(7.0e-4) + amrex::Real(1.0e-5) * i,
                               amrex::Real(3.5e-4) + amrex::Real(1.0e-5) * (j + k));
    }

    const amrex::Real tabs = amrex::Real(287.0) + amrex::Real(0.1) * i;
    const amrex::Real pres_mbar = amrex::Real(890.0) + amrex::Real(1.0) * k;
    return make_cell_state(tabs,
                           pres_mbar,
                           amrex::Real(3.5e-3) + amrex::Real(1.0e-4) * (i + k),
                           amrex::Real(1.2e-3) + amrex::Real(4.0e-5) * j);
}

// Fill every valid cell with a SatAdj-active state. These are deliberately not
// no-op states, so an unchanged active cell indicates skipped execution rather
// than physically valid SatAdj behavior.
inline void fill_satadj_active_conserved_state (amrex::MultiFab& cons)
{
    CellState evap_then_recond;
    const bool found = find_evaporation_then_recondensation_state(evap_then_recond);
    AMREX_ALWAYS_ASSERT(found);

    const int ncomp = cons.nComp();
    cons.setVal(amrex::Real(0.0));

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto arr = cons.array(mfi);

        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const CellState state = make_satadj_active_cell_state(i, j, k, evap_then_recond);
                const ConservedState conserved = make_conserved_state(state);

                arr(i,j,k,Rho_comp) = conserved.rho;
                arr(i,j,k,RhoTheta_comp) = conserved.rhotheta;
                arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
                arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
                arr(i,j,k,RhoQ1_comp) = conserved.rhoqv;
                arr(i,j,k,RhoQ2_comp) = conserved.rhoqc;

                if (ncomp > RhoQ3_comp) {
                    arr(i,j,k,RhoQ3_comp) = amrex::Real(0.0);
                }
            });
        });
    }
}

inline void fill_satadj_kessler_active_conserved_state (amrex::MultiFab& cons)
{
    const int ncomp = cons.nComp();
    cons.setVal(amrex::Real(0.0));

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto arr = cons.array(mfi);

        run_and_sync([=] {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const CellState state = make_satadj_kessler_active_cell_state(i, j, k);
                const ConservedState conserved = make_conserved_state(state);

                arr(i,j,k,Rho_comp) = conserved.rho;
                arr(i,j,k,RhoTheta_comp) = conserved.rhotheta;
                arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
                arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
                arr(i,j,k,RhoQ1_comp) = conserved.rhoqv;
                arr(i,j,k,RhoQ2_comp) = conserved.rhoqc;

                if (ncomp > RhoQ3_comp) {
                    arr(i,j,k,RhoQ3_comp) = amrex::Real(0.0);
                }
            });
        });
    }
}

inline ConservedState scalar_reference_from_initial_conserved (const amrex::Real rho,
                                                               const amrex::Real rhotheta,
                                                               const amrex::Real rhoqv,
                                                               const amrex::Real rhoqc)
{
    CellState reference_state = make_state_from_conserved(rho, rhotheta, rhoqv, rhoqc);
    adjust(reference_state);
    return ConservedState{rho,
                          rho * reference_state.theta,
                          rho * reference_state.qv,
                          rho * reference_state.qc};
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real conserved_component (const ConservedState& conserved,
                                 const int component)
{
    switch (component) {
    case Rho_comp: return conserved.rho;
    case RhoTheta_comp: return conserved.rhotheta;
    case RhoQ1_comp: return conserved.rhoqv;
    case RhoQ2_comp: return conserved.rhoqc;
    default: return amrex::Real(0.0);
    }
}

// The plotted or diagnosed temperature downstream of SatAdj may be
// reconstructed from conserved rho, rhoTheta, and qv. This helper compares
// that EOS-projected temperature against the scalar-reference conserved state,
// without assuming it equals the fixed-pressure Newton temperature.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real eos_projected_temperature_from_cons (const amrex::Real rho,
                                                 const amrex::Real rhotheta,
                                                 const amrex::Real rhoqv)
{
    return getTgivenRandRTh(rho, rhotheta, rhoqv / rho);
}

// Active-region tests infer execution from whether SatAdj-written components
// changed. Moisture components must use water-scaled tolerances so smaller but
// real qv/qc updates are not hidden by the generic O(1) tolerance floor.
inline amrex::Real satadj_component_change_tol (const int component,
                                                const amrex::Real rho_reference,
                                                const amrex::Real conserved_reference)
{
    if (component == RhoQ1_comp || component == RhoQ2_comp) {
        const amrex::Real q_reference = conserved_reference / rho_reference;
        return conserved_moisture_tol(rho_reference, q_reference);
    }

    return scaled_tol(conserved_reference, amrex::Real(10.0) * kStateTolFactor);
}

inline bool satadj_written_component_changed (const int component,
                                              const amrex::Real rho_reference,
                                              const amrex::Real before,
                                              const amrex::Real after)
{
    return amrex::Math::abs(after - before) >
           satadj_component_change_tol(component, rho_reference, before);
}

inline amrex::MultiFab make_host_mirror (const amrex::MultiFab& src)
{
    amrex::MultiFab host(src.boxArray(),
                         src.DistributionMap(),
                         src.nComp(),
                         src.nGrowVect(),
                         amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    amrex::MultiFab::Copy(host, src, 0, 0, src.nComp(), src.nGrowVect());
    sync();
    return host;
}

inline void run_satadj_public_flow (SatAdj& satadj,
                                    const SolverChoice& sc,
                                    const amrex::Geometry& geom,
                                    amrex::MultiFab& cons,
                                    const amrex::Real dt = amrex::Real(1.0))
{
    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;
    run_and_sync([&]() {
        satadj.Init(cons, cons.boxArray(), geom, dt, z_phys_nd, detJ_cc);
        satadj.Copy_State_to_Micro(cons);
        satadj.Advance(dt, sc);
        satadj.Copy_Micro_to_State(cons);
    });
}

inline void compare_public_flow_to_scalar_reference (const amrex::MultiFab& initial_cons,
                                                     const amrex::MultiFab& final_cons,
                                                     LocatedError& worst_error,
                                                     ActiveCellCounts& counts)
{
    const amrex::MultiFab initial_host = make_host_mirror(initial_cons);
    const amrex::MultiFab final_host = make_host_mirror(final_cons);

    worst_error = LocatedError{};
    counts = ActiveCellCounts{};

    for (amrex::MFIter mfi(final_host); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto initial_arr = initial_host.const_array(mfi);
        const auto final_arr = final_host.const_array(mfi);

        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    ++counts.checked;

                    const ConservedState reference = scalar_reference_from_initial_conserved(
                        initial_arr(i,j,k,Rho_comp),
                        initial_arr(i,j,k,RhoTheta_comp),
                        initial_arr(i,j,k,RhoQ1_comp),
                        initial_arr(i,j,k,RhoQ2_comp));

                    const amrex::Real rho_initial = initial_arr(i,j,k,Rho_comp);

                    const bool changed =
                        satadj_written_component_changed(RhoTheta_comp,
                                                         rho_initial,
                                                         initial_arr(i,j,k,RhoTheta_comp),
                                                         final_arr(i,j,k,RhoTheta_comp)) ||
                        satadj_written_component_changed(RhoQ1_comp,
                                                         rho_initial,
                                                         initial_arr(i,j,k,RhoQ1_comp),
                                                         final_arr(i,j,k,RhoQ1_comp)) ||
                        satadj_written_component_changed(RhoQ2_comp,
                                                         rho_initial,
                                                         initial_arr(i,j,k,RhoQ2_comp),
                                                         final_arr(i,j,k,RhoQ2_comp));

                    if (changed) {
                        ++counts.changed;
                    } else {
                        ++counts.unexpected_unchanged;
                    }

                    const amrex::Real eos_actual = eos_projected_temperature_from_cons(
                        final_arr(i,j,k,Rho_comp),
                        final_arr(i,j,k,RhoTheta_comp),
                        final_arr(i,j,k,RhoQ1_comp));
                    const amrex::Real eos_expected = eos_projected_temperature_from_cons(
                        reference.rho,
                        reference.rhotheta,
                        reference.rhoqv);

                    const amrex::Real rho_actual = final_arr(i,j,k,Rho_comp);
                    const amrex::Real rho_expected = reference.rho;

                    // RhoTheta: compare in conserved space
                    {
                        const amrex::Real actual = final_arr(i,j,k,RhoTheta_comp);
                        const amrex::Real expected = reference.rhotheta;
                        const amrex::Real absolute_error = amrex::Math::abs(actual - expected);
                        const amrex::Real normalized_error = absolute_error
                            / scaled_tol(expected, amrex::Real(20.0) * kStateTolFactor);
                        if (normalized_error > worst_error.normalized_error) {
                            worst_error = LocatedError{normalized_error,
                                                       absolute_error,
                                                       actual,
                                                       expected,
                                                       eos_actual,
                                                       eos_expected,
                                                       i,
                                                       j,
                                                       k,
                                                       RhoTheta_comp,
                                                       mfi.index()};
                        }
                    }

                    // RhoQ1: compare qv in mixing-ratio space
                    {
                        const amrex::Real qv_actual = final_arr(i,j,k,RhoQ1_comp) / rho_actual;
                        const amrex::Real qv_expected = reference.rhoqv / rho_expected;
                        const amrex::Real absolute_error = amrex::Math::abs(qv_actual - qv_expected);
                        const amrex::Real normalized_error = absolute_error
                            / mixing_ratio_tol(qv_expected);
                        if (normalized_error > worst_error.normalized_error) {
                            worst_error = LocatedError{normalized_error,
                                                       absolute_error,
                                                       qv_actual,
                                                       qv_expected,
                                                       eos_actual,
                                                       eos_expected,
                                                       i,
                                                       j,
                                                       k,
                                                       RhoQ1_comp,
                                                       mfi.index()};
                        }
                    }

                    // RhoQ2: compare qc in mixing-ratio space
                    {
                        const amrex::Real qc_actual = final_arr(i,j,k,RhoQ2_comp) / rho_actual;
                        const amrex::Real qc_expected = reference.rhoqc / rho_expected;
                        const amrex::Real absolute_error = amrex::Math::abs(qc_actual - qc_expected);
                        const amrex::Real normalized_error = absolute_error
                            / mixing_ratio_tol(qc_expected);
                        if (normalized_error > worst_error.normalized_error) {
                            worst_error = LocatedError{normalized_error,
                                                       absolute_error,
                                                       qc_actual,
                                                       qc_expected,
                                                       eos_actual,
                                                       eos_expected,
                                                       i,
                                                       j,
                                                       k,
                                                       RhoQ2_comp,
                                                       mfi.index()};
                        }
                    }

                    const amrex::Real eos_absolute_error = amrex::Math::abs(eos_actual - eos_expected);
                    const amrex::Real eos_normalized_error = eos_absolute_error
                        / derived_temperature_tol(eos_expected);
                    if (eos_normalized_error > worst_error.normalized_error) {
                        worst_error = LocatedError{eos_normalized_error,
                                                   eos_absolute_error,
                                                   eos_actual,
                                                   eos_expected,
                                                   eos_actual,
                                                   eos_expected,
                                                   i,
                                                   j,
                                                   k,
                                                   -1,
                                                   mfi.index()};
                    }
                }
            }
        }
    }
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
    amrex::Array<int, AMREX_SPACEDIM> periodicity{AMREX_D_DECL(1, 1, 1)};
    return make_multibox_geometry(nx, ny, nz, periodicity);
}

} // namespace satadj_test

#endif
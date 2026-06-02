#include <array>
#include <sstream>
#include <string>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_RealBox.H>

#include <gtest/gtest.h>

#include <ERF_SAM.H>

#include "ERF_GTestSAMCommon.H"

using namespace sam_test;

namespace {

SolverChoice make_sam_solver_choice (const MoistureType moisture_type = MoistureType::SAM)
{
    SolverChoice sc{};
    sc.c_p = Cp_d;
    sc.rdOcp = kRdOcp;
    sc.ave_plane = 2;
    sc.moisture_type = moisture_type;
    sc.use_shoc = false;
    return sc;
}

void fill_uniform_conserved_state (amrex::MultiFab& cons,
                                   const amrex::Real rho,
                                   const amrex::Real theta,
                                   const amrex::Real qv,
                                   const amrex::Real qcl,
                                   const amrex::Real qci,
                                   const amrex::Real qpr,
                                   const amrex::Real qps,
                                   const amrex::Real qpg)
{
    cons.setVal(amrex::Real(0.0));

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto arr = cons.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            arr(i,j,k,Rho_comp) = rho;
            arr(i,j,k,RhoTheta_comp) = rho * theta;
            arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
            arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
            arr(i,j,k,RhoQ1_comp) = rho * qv;
            arr(i,j,k,RhoQ2_comp) = rho * qcl;
            arr(i,j,k,RhoQ3_comp) = rho * qci;
            arr(i,j,k,RhoQ4_comp) = rho * qpr;
            arr(i,j,k,RhoQ5_comp) = rho * qps;
            arr(i,j,k,RhoQ6_comp) = rho * qpg;
        });
    }

    amrex::Gpu::streamSynchronize();
}

void expect_coefficient_row_near (const SAMCoefficientRow& actual,
                                  const SAMCoefficientRow& expected,
                                  const int k,
                                  const int rank)
{
    const std::string prefix = "rank=" + std::to_string(rank) + " k=" + std::to_string(k) + " ";
    EXPECT_NEAR(actual.accrrc, expected.accrrc, pow_sqrt_tol(expected.accrrc)) << prefix << "accrrc";
    EXPECT_NEAR(actual.accrsi, expected.accrsi, pow_sqrt_tol(expected.accrsi)) << prefix << "accrsi";
    EXPECT_NEAR(actual.accrsc, expected.accrsc, pow_sqrt_tol(expected.accrsc)) << prefix << "accrsc";
    EXPECT_NEAR(actual.coefice, expected.coefice, pow_sqrt_tol(expected.coefice)) << prefix << "coefice";
    EXPECT_NEAR(actual.evaps1, expected.evaps1, pow_sqrt_tol(expected.evaps1)) << prefix << "evaps1";
    EXPECT_NEAR(actual.evaps2, expected.evaps2, pow_sqrt_tol(expected.evaps2)) << prefix << "evaps2";
    EXPECT_NEAR(actual.accrgi, expected.accrgi, pow_sqrt_tol(expected.accrgi)) << prefix << "accrgi";
    EXPECT_NEAR(actual.accrgc, expected.accrgc, pow_sqrt_tol(expected.accrgc)) << prefix << "accrgc";
    EXPECT_NEAR(actual.evapg1, expected.evapg1, pow_sqrt_tol(expected.evapg1)) << prefix << "evapg1";
    EXPECT_NEAR(actual.evapg2, expected.evapg2, pow_sqrt_tol(expected.evapg2)) << prefix << "evapg2";
    EXPECT_NEAR(actual.evapr1, expected.evapr1, pow_sqrt_tol(expected.evapr1)) << prefix << "evapr1";
    EXPECT_NEAR(actual.evapr2, expected.evapr2, pow_sqrt_tol(expected.evapr2)) << prefix << "evapr2";
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real patterned_rho (const int i,
                           const int j,
                           const int k) noexcept
{
    return amrex::Real(1.0) + amrex::Real(0.01) * i + amrex::Real(0.02) * j + amrex::Real(0.03) * k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real patterned_theta (const int i,
                             const int j,
                             const int k) noexcept
{
    return amrex::Real(300.0) + i + amrex::Real(2.0) * j + amrex::Real(3.0) * k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real patterned_qv (const int i,
                          const int j,
                          const int k) noexcept
{
    return amrex::Real(1.0e-2) + amrex::Real(1.0e-4) * (i + 2 * j + 3 * k);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real patterned_qcl (const int i,
                           const int j,
                           const int k) noexcept
{
    return amrex::Real(2.0e-4) + amrex::Real(1.0e-5) * (2 * i + j + k);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real patterned_qci (const int i,
                           const int j,
                           const int k) noexcept
{
    return amrex::Real(1.0e-4) + amrex::Real(1.0e-5) * (i + 3 * j + 2 * k);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real patterned_qpr (const int i,
                           const int j,
                           const int) noexcept
{
    return amrex::Real(3.0e-5) + amrex::Real(1.0e-6) * (4 * i + j);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real patterned_qps (const int i,
                           const int,
                           const int k) noexcept
{
    return amrex::Real(4.0e-5) + amrex::Real(1.0e-6) * (i + 5 * k);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real patterned_qpg (const int,
                           const int j,
                           const int k) noexcept
{
    return amrex::Real(5.0e-5) + amrex::Real(1.0e-6) * (2 * j + 3 * k);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real patterned_cons_value (const int comp,
                                  const int i,
                                  const int j,
                                  const int k) noexcept
{
    const amrex::Real rho = patterned_rho(i, j, k);
    const amrex::Real theta = patterned_theta(i, j, k);

    switch (comp) {
    case Rho_comp: return rho;
    case RhoTheta_comp: return rho * theta;
    case RhoKE_comp: return amrex::Real(0.0);
    case RhoScalar_comp: return amrex::Real(0.0);
    case RhoQ1_comp: return rho * patterned_qv(i, j, k);
    case RhoQ2_comp: return rho * patterned_qcl(i, j, k);
    case RhoQ3_comp: return rho * patterned_qci(i, j, k);
    case RhoQ4_comp: return rho * patterned_qpr(i, j, k);
    case RhoQ5_comp: return rho * patterned_qps(i, j, k);
    case RhoQ6_comp: return rho * patterned_qpg(i, j, k);
    default: return amrex::Real(0.0);
    }
}

void fill_patterned_conserved_state (amrex::MultiFab& cons)
{
    cons.setVal(amrex::Real(-999.0));

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto arr = cons.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int comp = Rho_comp; comp <= RhoQ6_comp; ++comp) {
                arr(i,j,k,comp) = patterned_cons_value(comp, i, j, k);
            }
        });
    }

    amrex::Gpu::streamSynchronize();
}

void poison_ghost_cells (amrex::MultiFab& cons,
                         const amrex::Real sentinel)
{
    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& fab_box = mfi.fabbox();
        const amrex::Box& valid_box = mfi.validbox();
        auto arr = cons.array(mfi);
        amrex::ParallelFor(fab_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (!valid_box.contains(amrex::IntVect(AMREX_D_DECL(i, j, k)))) {
                for (int comp = Rho_comp; comp <= RhoQ6_comp; ++comp) {
                    arr(i,j,k,comp) = sentinel;
                }
            }
        });
    }

    amrex::Gpu::streamSynchronize();
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real mixed_phase_tabs_value () noexcept
{
    const amrex::Real lower = amrex::max(tprmin, tgrmin) + amrex::Real(1.0e-3);
    const amrex::Real upper = amrex::min(tprmax, tgrmax) - amrex::Real(1.0e-3);
    return amrex::Real(0.5) * (lower + upper);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real mixed_qsat_for_state (const MoistureType moisture_type,
                                  const amrex::Real tabs,
                                  const amrex::Real pres_mbar) noexcept
{
    amrex::Real qsatw;
    amrex::Real qsati;
    erf_qsatw(tabs, pres_mbar, qsatw);
    erf_qsati(tabs, pres_mbar, qsati);
    const int sam_mode = sam_is_no_ice(moisture_type) ? kSAMNoIceMode : kSAMWithIceMode;
    const amrex::Real omn = sam_cloud_liquid_fraction(sam_mode, tabs, a_bg, tbgmin * a_bg);
    return sam_mixed_qsat(omn, qsatw, qsati);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMCellState make_precip_budget_state (const int i,
                                       const int j,
                                       const int k,
                                       const MoistureType moisture_type,
                                       const amrex::Real pres_mbar) noexcept
{
    const int pattern = (i + 2 * j + 3 * k) % 4;

    amrex::Real tabs = amrex::Real(285.0);
    amrex::Real vapor_fraction = amrex::Real(0.75);
    amrex::Real qcl = amrex::Real(0.0);
    amrex::Real qci = amrex::Real(0.0);
    amrex::Real qpr = amrex::Real(0.0);
    amrex::Real qps = amrex::Real(0.0);
    amrex::Real qpg = amrex::Real(0.0);

    switch (pattern) {
    case 0:
        tabs = amrex::Real(285.0);
        vapor_fraction = amrex::Real(0.75);
        qcl = qcw0 + amrex::Real(5.0e-4);
        qpr = amrex::Real(3.0e-4);
        break;
    case 1:
        tabs = mixed_phase_tabs_value();
        vapor_fraction = amrex::Real(0.6);
        qcl = qcw0 + amrex::Real(8.0e-4);
        qci = qci0 + amrex::Real(7.0e-4);
        qpr = amrex::Real(4.0e-4);
        qps = amrex::Real(5.0e-4);
        qpg = amrex::Real(6.0e-4);
        break;
    case 2:
        tabs = amrex::Real(260.0);
        vapor_fraction = amrex::Real(0.7);
        qci = qci0 + amrex::Real(4.0e-4);
        qps = amrex::Real(4.0e-4);
        break;
    default:
        tabs = amrex::Real(245.0);
        vapor_fraction = amrex::Real(0.7);
        qci = qci0 + amrex::Real(4.0e-4);
        qpg = amrex::Real(4.0e-4);
        break;
    }

    if (sam_is_no_ice(moisture_type)) {
        qci = amrex::Real(0.0);
        qps = amrex::Real(0.0);
        qpg = amrex::Real(0.0);
        tabs = amrex::max(tbgmax + amrex::Real(1.0), tabs);
    }

    const amrex::Real qv = vapor_fraction * mixed_qsat_for_state(moisture_type, tabs, pres_mbar);

    SAMCellState state{};
    state.tabs = tabs;
    state.pres_mbar = pres_mbar;
    state.qv = qv;
    state.qcl = qcl;
    state.qci = qci;
    state.qpr = qpr;
    state.qps = qps;
    state.qpg = qpg;
    state.qn = qcl + qci;
    state.qt = qv + state.qn;
    state.qp = qpr + qps + qpg;
    state.theta = sam_theta_from_stored_mbar_converted_to_pa(tabs, pres_mbar, kRdOcp);
    state.rho = getRhogivenTandPress(tabs, sam_mbar_to_pa(pres_mbar), qv);
    return state;
}

void fill_precip_budget_conserved_state (amrex::MultiFab& cons,
                                         const MoistureType moisture_type,
                                         const amrex::Real pres_mbar)
{
    cons.setVal(amrex::Real(0.0));

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto arr = cons.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const SAMCellState state = make_precip_budget_state(i, j, k, moisture_type, pres_mbar);
            arr(i,j,k,Rho_comp) = state.rho;
            arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
            arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
            sam_primitive_to_cons(state, arr, i, j, k);
        });
    }

    amrex::Gpu::streamSynchronize();
}

void fill_nonuniform_detj (amrex::MultiFab& detJ_cc)
{
    for (amrex::MFIter mfi(detJ_cc, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto arr = detJ_cc.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            arr(i,j,k) = amrex::Real(0.8) + amrex::Real(0.05) * i + amrex::Real(0.03) * j + amrex::Real(0.07) * k;
        });
    }

    amrex::Gpu::streamSynchronize();
}

struct IceFallFluxReductionDiagnostics {
    amrex::Real initial_qci_mass{amrex::Real(0.0)};
    amrex::Real local_max_flux{amrex::Real(0.0)};
    amrex::Real global_max_flux{amrex::Real(0.0)};
    int local_n_substep{0};
    int global_n_substep{0};
    int owned_box_count{0};
    std::string box_summary{"none"};
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMCellState make_icefall_reduction_state (const bool high_flux_region,
                                          const int k,
                                          const int k_hi) noexcept
{
    const amrex::Real tabs = amrex::Real(245.0);
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qv = amrex::Real(1.0e-2);
    const amrex::Real qci_scale = amrex::Real(k_hi + 1 - k);

    SAMCellState state{};
    state.tabs = tabs;
    state.pres_mbar = pres_mbar;
    state.qv = qv;
    state.qcl = zero;
    state.qci = (high_flux_region ? amrex::Real(5.0e-4) : amrex::Real(5.0e-6)) * qci_scale;
    state.qpr = zero;
    state.qps = zero;
    state.qpg = zero;
    state.qn = state.qcl + state.qci;
    state.qt = state.qv + state.qn;
    state.qp = zero;
    state.theta = sam_theta_from_stored_mbar_converted_to_pa(tabs, pres_mbar, kRdOcp);
    state.rho = getRhogivenTandPress(tabs, sam_mbar_to_pa(pres_mbar), qv);
    return state;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real icefall_mixing_ratio_from_conserved (const amrex::Real rho_qci,
                                                 const amrex::Real rho) noexcept
{
    return rho_qci / rho;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real icefall_face_flux (const amrex::Real rho_face,
                               const amrex::Real qci_face) noexcept
{
    return rho_face * sam_cloud_ice_terminal_velocity(qci_face) * qci_face;
}

void fill_icefall_reduction_conserved_state (amrex::MultiFab& cons,
                                             const amrex::Geometry& geom)
{
    const int k_hi = geom.Domain().bigEnd(2);
    const int i_hi = geom.Domain().bigEnd(0);

    cons.setVal(amrex::Real(0.0));

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto arr = cons.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const bool high_flux_region = (i >= i_hi - 1);
            const SAMCellState state = make_icefall_reduction_state(high_flux_region, k, k_hi);

            arr(i,j,k,Rho_comp) = state.rho;
            arr(i,j,k,RhoTheta_comp) = state.rho * state.theta;
            arr(i,j,k,RhoKE_comp) = zero;
            arr(i,j,k,RhoScalar_comp) = zero;
            sam_primitive_to_cons(state, arr, i, j, k);
        });
    }

    amrex::Gpu::streamSynchronize();
}

IceFallFluxReductionDiagnostics compute_icefall_flux_reduction_diagnostics (
    const amrex::Geometry& geom,
    const amrex::MultiFab& cons,
    const amrex::Real dt)
{
    IceFallFluxReductionDiagnostics diagnostics{};
    const amrex::Box domain = geom.Domain();
    const int k_lo = domain.smallEnd(2);
    const int k_hi = domain.bigEnd(2);
    const amrex::Real dz = geom.CellSize(2);
    const amrex::Real cell_volume = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);

    diagnostics.initial_qci_mass = cons.sum(RhoQ3_comp) * cell_volume;

    amrex::MultiFab fz;
    const amrex::IntVect ng = cons.nGrowVect();
    const amrex::BoxArray ba = cons.boxArray();
    const amrex::DistributionMapping dm = cons.DistributionMap();
    fz.define(amrex::convert(ba, amrex::IntVect(0, 0, 1)), dm, 1, ng);
    fz.setVal(amrex::Real(0.0));

    for (amrex::MFIter mfi(fz, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto cons_array = cons.const_array(mfi);
        auto fz_array = fz.array(mfi);

        const auto& box3d = mfi.tilebox();

        amrex::ParallelFor(box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real rho_km1 = (k == k_lo) ? cons_array(i,j,k,Rho_comp)
                                                    : cons_array(i,j,k-1,Rho_comp);
            const amrex::Real rho_k = (k == k_hi + 1) ? cons_array(i,j,k-1,Rho_comp)
                                                      : cons_array(i,j,k,Rho_comp);
            const amrex::Real qci_km1 = (k == k_lo)
                ? icefall_mixing_ratio_from_conserved(cons_array(i,j,k,RhoQ3_comp), rho_km1)
                : icefall_mixing_ratio_from_conserved(cons_array(i,j,k-1,RhoQ3_comp), rho_km1);
            const amrex::Real qci_k = (k == k_hi + 1)
                ? icefall_mixing_ratio_from_conserved(cons_array(i,j,k-1,RhoQ3_comp), rho_k)
                : icefall_mixing_ratio_from_conserved(cons_array(i,j,k,RhoQ3_comp), rho_k);
            const SAMFaceState face_state = sam_face_average_state(k, k_lo, k_hi,
                                                                   rho_km1, rho_k,
                                                                   zero, zero,
                                                                   qci_km1, qci_k,
                                                                   zero, zero);
            fz_array(i,j,k) = icefall_face_flux(face_state.rho_avg, face_state.qci_avg);
        });
    }

    amrex::Gpu::streamSynchronize();

    const auto fz_arrays = fz.const_arrays();
    const amrex::GpuTuple<amrex::Real> max_reduced_flux = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpMax>{},
        amrex::TypeList<amrex::Real>{},
        fz, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Real> {
            return {fz_arrays[box_no](i,j,k)};
        });

    diagnostics.local_max_flux = amrex::max(amrex::Real(0.0), amrex::get<0>(max_reduced_flux));
    diagnostics.global_max_flux = diagnostics.local_max_flux;
    amrex::ParallelDescriptor::ReduceRealMax(diagnostics.global_max_flux);

    diagnostics.local_n_substep = sam_substep_count_from_reduced_flux(
        diagnostics.local_max_flux + std::numeric_limits<amrex::Real>::epsilon(),
        dt,
        dz);
    diagnostics.global_n_substep = sam_substep_count_from_reduced_flux(
        diagnostics.global_max_flux + std::numeric_limits<amrex::Real>::epsilon(),
        dt,
        dz);

    std::ostringstream box_stream;
    bool first_box = true;
    for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
        if (!first_box) {
            box_stream << ";";
        }
        first_box = false;
        ++diagnostics.owned_box_count;
        box_stream << mfi.validbox();
    }
    diagnostics.box_summary = first_box ? std::string("none") : box_stream.str();

    return diagnostics;
}

void compute_icefall_reference_global_substep (const amrex::Geometry& geom,
                                               amrex::MultiFab& reference_cons,
                                               const int n_substep,
                                               const amrex::Real dt)
{
    const amrex::Box domain = geom.Domain();
    const int k_lo = domain.smallEnd(2);
    const int k_hi = domain.bigEnd(2);
    const amrex::Real coef = (dt / geom.CellSize(2)) / amrex::Real(n_substep);

    amrex::MultiFab fz;
    const amrex::IntVect ng = reference_cons.nGrowVect();
    const amrex::BoxArray ba = reference_cons.boxArray();
    const amrex::DistributionMapping dm = reference_cons.DistributionMap();
    fz.define(amrex::convert(ba, amrex::IntVect(0, 0, 1)), dm, 1, ng);
    fz.setVal(amrex::Real(0.0));

    for (int nsub(0); nsub<n_substep; ++nsub) {
        for (amrex::MFIter mfi(fz, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto cons_array = reference_cons.const_array(mfi);
            auto fz_array = fz.array(mfi);

            const auto& box3d = mfi.tilebox();

            amrex::ParallelFor(box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const amrex::Real rho_km1 = (k == k_lo) ? cons_array(i,j,k,Rho_comp)
                                                        : cons_array(i,j,k-1,Rho_comp);
                const amrex::Real rho_k = (k == k_hi + 1) ? cons_array(i,j,k-1,Rho_comp)
                                                          : cons_array(i,j,k,Rho_comp);
                const amrex::Real qci_km1 = (k == k_lo)
                    ? icefall_mixing_ratio_from_conserved(cons_array(i,j,k,RhoQ3_comp), rho_km1)
                    : icefall_mixing_ratio_from_conserved(cons_array(i,j,k-1,RhoQ3_comp), rho_km1);
                const amrex::Real qci_k = (k == k_hi + 1)
                    ? icefall_mixing_ratio_from_conserved(cons_array(i,j,k-1,RhoQ3_comp), rho_k)
                    : icefall_mixing_ratio_from_conserved(cons_array(i,j,k,RhoQ3_comp), rho_k);
                const SAMFaceState face_state = sam_face_average_state(k, k_lo, k_hi,
                                                                       rho_km1, rho_k,
                                                                       zero, zero,
                                                                       qci_km1, qci_k,
                                                                       zero, zero);
                fz_array(i,j,k) = icefall_face_flux(face_state.rho_avg, face_state.qci_avg);
            });
        }

        for (amrex::MFIter mfi(reference_cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto cons_array = reference_cons.array(mfi);
            auto const_cons_array = reference_cons.const_array(mfi);
            auto fz_array = fz.array(mfi);

            const auto& tbx = mfi.tilebox();

            amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const amrex::Real rho = const_cons_array(i,j,k,Rho_comp);
                const amrex::Real qci_mass = const_cons_array(i,j,k,RhoQ3_comp);
                const amrex::Real qci = icefall_mixing_ratio_from_conserved(qci_mass, rho);

                amrex::Real dqi = sam_sedimentation_tendency(fz_array(i,j,k+1), fz_array(i,j,k),
                                                             rho, one, coef);
                dqi = std::max(-qci, dqi);

                cons_array(i,j,k,RhoQ3_comp) = qci_mass + rho * dqi;
            });
        }
    }
}

std::string icefall_state_summary (const char* label,
                                   const SAMCellState& state,
                                   const amrex::Real raw_qci)
{
    std::ostringstream stream;
    stream << label
           << " raw_qci=" << raw_qci
           << " rho=" << state.rho
           << " theta=" << state.theta
           << " tabs=" << state.tabs
           << " pres_mbar=" << state.pres_mbar
           << " qv=" << state.qv
           << " qcl=" << state.qcl
           << " qci=" << state.qci
           << " qn=" << state.qn
           << " qt=" << state.qt
           << " qpr=" << state.qpr
           << " qps=" << state.qps
           << " qpg=" << state.qpg
           << " qp=" << state.qp;
    return stream.str();
}

struct PrecipFallCellState {
    amrex::Real rho;
    amrex::Real theta;
    amrex::Real tabs;
    amrex::Real qv;
    amrex::Real qpr;
    amrex::Real qps;
    amrex::Real qpg;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
PrecipFallCellState precipfall_cell_state (const int i,
                                           const int j,
                                           const amrex::Real pressure_pa) noexcept
{
    const int phase_index = (i + 2 * j) % 3;
    const amrex::Real tabs = (phase_index == 0) ? amrex::Real(281.0)
                            : (phase_index == 1) ? amrex::Real(263.0)
                                                 : amrex::Real(240.0);
    const amrex::Real qv = amrex::Real(9.0e-3) + amrex::Real(2.0e-4) * phase_index;
    const amrex::Real theta = getThgivenTandP(tabs, pressure_pa, kRdOcp);
    const amrex::Real rho = getRhogivenTandPress(tabs, pressure_pa, qv);

    return {
        rho,
        theta,
        tabs,
        qv,
        amrex::Real(2.0e-5) + amrex::Real(5.0e-7) * i + amrex::Real(3.0e-7) * j,
        amrex::Real(1.5e-5) + amrex::Real(2.0e-7) * i + amrex::Real(4.0e-7) * j,
        amrex::Real(1.0e-5) + amrex::Real(3.0e-7) * i + amrex::Real(2.0e-7) * j};
}

std::array<amrex::Real, 3> precipfall_terminal_velocities ()
{
    const amrex::Real gamr3 = erf_gammafff(amrex::Real(4.0) + b_rain);
    const amrex::Real gams3 = erf_gammafff(amrex::Real(4.0) + b_snow);
    const amrex::Real gamg3 = erf_gammafff(amrex::Real(4.0) + b_grau);

    return {
        (a_rain * gamr3 / amrex::Real(6.0)) * std::pow((PI * rhor * nzeror), -crain),
        (a_snow * gams3 / amrex::Real(6.0)) * std::pow((PI * rhos * nzeros), -csnow),
        (a_grau * gamg3 / amrex::Real(6.0)) * std::pow((PI * rhog * nzerog), -cgrau)};
}

void fill_precipfall_conserved_state (amrex::MultiFab& cons,
                                      const amrex::Real pressure_pa)
{
    cons.setVal(amrex::Real(0.0));

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto arr = cons.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const PrecipFallCellState state = precipfall_cell_state(i, j, pressure_pa);

            arr(i,j,k,Rho_comp) = state.rho;
            arr(i,j,k,RhoTheta_comp) = state.rho * state.theta;
            arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
            arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
            arr(i,j,k,RhoQ1_comp) = state.rho * state.qv;
            arr(i,j,k,RhoQ2_comp) = amrex::Real(0.0);
            arr(i,j,k,RhoQ3_comp) = amrex::Real(0.0);
            arr(i,j,k,RhoQ4_comp) = state.rho * state.qpr;
            arr(i,j,k,RhoQ5_comp) = state.rho * state.qps;
            arr(i,j,k,RhoQ6_comp) = state.rho * state.qpg;
        });
    }

    amrex::Gpu::streamSynchronize();
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMCellState make_precipfall_component_budget_state (const int i,
                                                     const int j,
                                                     const int k,
                                                     const int nz,
                                                     const amrex::Real pres_mbar) noexcept
{
    const bool is_top = (k == nz - 1);
    const int phase_index = (i + j + k) % 3;
    const amrex::Real tabs = (phase_index == 0) ? amrex::Real(282.0)
                            : (phase_index == 1) ? amrex::Real(262.0)
                                                 : amrex::Real(242.0);
    const amrex::Real qsat = mixed_qsat_for_state(MoistureType::SAM, tabs, pres_mbar);
    const amrex::Real scale = is_top ? amrex::Real(0.0) : amrex::Real(nz - 1 - k);

    SAMCellState state{};
    state.tabs = tabs;
    state.pres_mbar = pres_mbar;
    state.qv = amrex::Real(0.75) * qsat;
    state.qcl = amrex::Real(0.0);
    state.qci = amrex::Real(0.0);
    state.qpr = is_top ? amrex::Real(0.0)
                       : amrex::Real(2.0e-5) + scale * (amrex::Real(3.0e-5) + amrex::Real(1.0e-6) * i);
    state.qps = is_top ? amrex::Real(0.0)
                       : amrex::Real(1.5e-5) + scale * (amrex::Real(2.0e-5) + amrex::Real(1.0e-6) * j);
    state.qpg = is_top ? amrex::Real(0.0)
                       : amrex::Real(1.0e-5) + scale * (amrex::Real(2.5e-5) + amrex::Real(5.0e-7) * (i + j));
    state.qn = amrex::Real(0.0);
    state.qt = state.qv;
    state.qp = state.qpr + state.qps + state.qpg;
    state.theta = sam_theta_from_stored_mbar_converted_to_pa(tabs, pres_mbar, kRdOcp);
    state.rho = getRhogivenTandPress(tabs, sam_mbar_to_pa(pres_mbar), state.qv);
    return state;
}

void fill_precipfall_component_budget_conserved_state (amrex::MultiFab& cons,
                                                       const amrex::Real pres_mbar)
{
    const int nz = cons.boxArray().minimalBox().length(2);
    cons.setVal(amrex::Real(0.0));

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto arr = cons.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const SAMCellState state = make_precipfall_component_budget_state(i, j, k, nz, pres_mbar);
            arr(i,j,k,Rho_comp) = state.rho;
            arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
            arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
            sam_primitive_to_cons(state, arr, i, j, k);
        });
    }

    amrex::Gpu::streamSynchronize();
}

amrex::Real sum_component_mass (const amrex::Geometry& geom,
                                const amrex::MultiFab& cons,
                                const int comp,
                                const amrex::MultiFab* detJ_cc = nullptr)
{
    const amrex::Real cell_volume = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
    if (detJ_cc == nullptr) {
        return cons.sum(comp) * cell_volume;
    }

    auto cons_arrays = cons.const_arrays();
    auto detj_arrays = detJ_cc->const_arrays();
    amrex::GpuTuple<amrex::Real> reduced = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpSum>{},
        amrex::TypeList<amrex::Real>{},
        cons, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Real> {
            return {cell_volume * detj_arrays[box_no](i,j,k) * cons_arrays[box_no](i,j,k,comp)};
        });
    amrex::Real total = amrex::get<0>(reduced);
    amrex::ParallelDescriptor::ReduceRealSum(total);
    return total;
}

amrex::Real sum_accum_mass (const amrex::Geometry& geom,
                            const amrex::MultiFab& accum,
                            const amrex::Real precip_density)
{
    const amrex::Real cell_area = geom.CellSize(0) * geom.CellSize(1);
    return accum.sum(0) * cell_area * precip_density / amrex::Real(1000.0);
}

struct PrecipFallComponentBudgetResult {
    amrex::Real rain_residual{amrex::Real(0.0)};
    amrex::Real snow_residual{amrex::Real(0.0)};
    amrex::Real graupel_residual{amrex::Real(0.0)};
    amrex::Real total_residual{amrex::Real(0.0)};
    amrex::Real rain_accum_mass{amrex::Real(0.0)};
    amrex::Real snow_accum_mass{amrex::Real(0.0)};
    amrex::Real graupel_accum_mass{amrex::Real(0.0)};
    amrex::Real detj_min{amrex::Real(1.0)};
    amrex::Real detj_max{amrex::Real(1.0)};
    amrex::IntVect max_size{AMREX_D_DECL(0, 0, 0)};
};

PrecipFallComponentBudgetResult run_precipfall_component_budget_case (const amrex::IntVect& max_size)
{
    constexpr int nx = 8;
    constexpr int ny = 4;
    constexpr int nz = 4;
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real dt = amrex::Real(0.25);

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(8.0, 4.0, 4.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 0)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    ba.maxSize(max_size);
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    fill_precipfall_component_budget_conserved_state(cons, pres_mbar);

    amrex::MultiFab detJ_cc(ba, dm, 1, 0);
    fill_nonuniform_detj(detJ_cc);

    const amrex::Real initial_rain_mass = sum_component_mass(geom, cons, RhoQ4_comp, &detJ_cc);
    const amrex::Real initial_snow_mass = sum_component_mass(geom, cons, RhoQ5_comp, &detJ_cc);
    const amrex::Real initial_graupel_mass = sum_component_mass(geom, cons, RhoQ6_comp, &detJ_cc);

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detj_owned = std::make_unique<amrex::MultiFab>(detJ_cc.boxArray(), detJ_cc.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(*detj_owned, detJ_cc, 0, 0, 1, 0);

    SolverChoice sc = make_sam_solver_choice(MoistureType::SAM);
    SAM sam;
    sam.Define(sc);
    sam.Set_dzmin(geom.CellSize(2));
    sam.Set_RealWidth(0);
    sam.Init(cons, ba, geom, dt, z_phys_nd, detj_owned);
    sam.Copy_State_to_Micro(cons);
    sam.PrecipFall(sc);
    sam.Copy_Micro_to_State(cons);

    PrecipFallComponentBudgetResult result{};
    result.rain_accum_mass = sum_accum_mass(geom, *sam.Qmoist_Ptr(0), rhor);
    result.snow_accum_mass = sum_accum_mass(geom, *sam.Qmoist_Ptr(1), rhos);
    result.graupel_accum_mass = sum_accum_mass(geom, *sam.Qmoist_Ptr(2), rhog);

    const amrex::Real final_rain_mass = sum_component_mass(geom, cons, RhoQ4_comp, &detJ_cc);
    const amrex::Real final_snow_mass = sum_component_mass(geom, cons, RhoQ5_comp, &detJ_cc);
    const amrex::Real final_graupel_mass = sum_component_mass(geom, cons, RhoQ6_comp, &detJ_cc);
    result.rain_residual = initial_rain_mass - final_rain_mass - result.rain_accum_mass;
    result.snow_residual = initial_snow_mass - final_snow_mass - result.snow_accum_mass;
    result.graupel_residual = initial_graupel_mass - final_graupel_mass - result.graupel_accum_mass;
    result.total_residual = result.rain_residual + result.snow_residual + result.graupel_residual;
    result.max_size = max_size;

    for (amrex::MFIter mfi(detJ_cc); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto detj = detJ_cc.const_array(mfi);
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    result.detj_min = std::min(result.detj_min, detj(i,j,k));
                    result.detj_max = std::max(result.detj_max, detj(i,j,k));
                }
            }
        }
    }
    amrex::ParallelDescriptor::ReduceRealMin(result.detj_min);
    amrex::ParallelDescriptor::ReduceRealMax(result.detj_max);

    return result;
}

int wrap_index (const int idx,
                const int lo,
                const int hi)
{
    const int len = hi - lo + 1;
    if (idx < lo) {
        return idx + len;
    }
    if (idx > hi) {
        return idx - len;
    }
    return idx;
}

} // namespace

// Motivation:
// Compute_Coefficients forms plane-averaged SAM coefficient rows. For a
// horizontally uniform state, those rows should be independent of rank count
// and box decomposition; CTest runs this same test under 1, 2, and 4 ranks
// when configured.
TEST(SAMParallel, ComputeCoefficientsRankInvariant)
{
    constexpr int nx = 8;
    constexpr int ny = 4;
    constexpr int nz = 4;

    const amrex::Real rho = amrex::Real(1.05);
    const amrex::Real theta = amrex::Real(296.0);
    const amrex::Real qv = amrex::Real(1.1e-2);
    const amrex::Real qcl = amrex::Real(4.0e-4);
    const amrex::Real qci = amrex::Real(2.0e-4);

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(1.0, 1.0, 1.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    ba.maxSize(amrex::IntVect(4, 2, nz));
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    fill_uniform_conserved_state(cons, rho, theta, qv, qcl, qci,
                                 amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0));

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;

    SolverChoice sc = make_sam_solver_choice();
    SAM sam;
    sam.Define(sc);
    sam.Set_dzmin(amrex::Real(100.0));
    sam.Set_RealWidth(0);
    sam.Init(cons, ba, geom, amrex::Real(1.0), z_phys_nd, detJ_cc);
    sam.Copy_State_to_Micro(cons);
    sam.Compute_Coefficients();

    const amrex::Real tabs = amrex::min(getTgivenRandRTh(rho, rho * theta, qv), amrex::Real(273.16));
    const amrex::Real gamr1 = erf_gammafff(three + b_rain);
    const amrex::Real gamr2 = erf_gammafff((amrex::Real(5.0) + b_rain) / two);
    const amrex::Real gams1 = erf_gammafff(three + b_snow);
    const amrex::Real gams2 = erf_gammafff((amrex::Real(5.0) + b_snow) / two);
    const amrex::Real gamg1 = erf_gammafff(three + b_grau);
    const amrex::Real gamg2 = erf_gammafff((amrex::Real(5.0) + b_grau) / two);
    const SAMCoefficientRow expected = sam_compute_coefficient_row(
        rho, tabs, gamr1, gamr2, gams1, gams2, gamg1, gamg2);

    const int rank = amrex::ParallelDescriptor::MyProc();
    for (int k = domain.smallEnd(2); k <= domain.bigEnd(2); ++k) {
        const SAMCoefficientRow actual = sam.CoefficientRowAt(k);
        expect_coefficient_row_near(actual, expected, k, rank);
    }
}

// Motivation:
// Copy_Micro_to_State ends with FillBoundary over the SAM geometry periodicity.
// After a state-to-micro roundtrip, periodic ghost cells should match the
// wrapped valid neighbors regardless of rank count or box ownership.
TEST(SAMParallel, CopyMicroToStateFillBoundaryParallel)
{
    constexpr int nx = 8;
    constexpr int ny = 4;
    constexpr int nz = 4;
    constexpr int ng = 1;
    const amrex::Real kSentinel = amrex::Real(-999.0);

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(1.0, 1.0, 1.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    ba.maxSize(amrex::IntVect(4, 2, nz));
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, ng);
    fill_patterned_conserved_state(cons);

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;

    SolverChoice sc = make_sam_solver_choice();
    SAM sam;
    sam.Define(sc);
    sam.Set_dzmin(amrex::Real(100.0));
    sam.Set_RealWidth(0);
    sam.Init(cons, ba, geom, amrex::Real(1.0), z_phys_nd, detJ_cc);
    sam.Copy_State_to_Micro(cons);

    poison_ghost_cells(cons, kSentinel);
    sam.Copy_Micro_to_State(cons);
    amrex::Gpu::streamSynchronize();

    const int rank = amrex::ParallelDescriptor::MyProc();
    int local_checked = 0;

    for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const amrex::Box& fab_box = mfi.fabbox();
        const amrex::Array4<const amrex::Real> arr = cons.const_array(mfi);

        for (int k = fab_box.smallEnd(2); k <= fab_box.bigEnd(2); ++k) {
            for (int j = fab_box.smallEnd(1); j <= fab_box.bigEnd(1); ++j) {
                for (int i = fab_box.smallEnd(0); i <= fab_box.bigEnd(0); ++i) {
                    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                    if (domain.contains(iv)) {
                        continue;
                    }

                    ++local_checked;
                    const int wrapped_i = wrap_index(i, domain.smallEnd(0), domain.bigEnd(0));
                    const int wrapped_j = wrap_index(j, domain.smallEnd(1), domain.bigEnd(1));
                    const int wrapped_k = wrap_index(k, domain.smallEnd(2), domain.bigEnd(2));

                    for (int comp = Rho_comp; comp <= RhoQ6_comp; ++comp) {
                        const amrex::Real expected = patterned_cons_value(comp, wrapped_i, wrapped_j, wrapped_k);
                        EXPECT_NEAR(arr(i,j,k,comp), expected, roundoff_tol(expected))
                            << "rank=" << rank
                            << " ghost_iv=" << iv
                            << " wrapped_iv=(" << wrapped_i << "," << wrapped_j << "," << wrapped_k << ")"
                            << " comp=" << comp;
                    }
                }
            }
        }
    }

    amrex::ParallelDescriptor::ReduceIntSum(local_checked);
    EXPECT_GT(local_checked, 0);
}

// Motivation:
// IceFall chooses a vertical sedimentation substep count from the maximum
// cloud-ice sedimentation flux. That maximum must be global across MPI ranks,
// not rank-local, so every rank advances the same physical timestep with the
// same substep count. This regression constructs a heterogeneous cloud-ice
// field whose maximum flux is present only on some ranks and checks that the
// public IceFall result matches an independently computed global-max substep
// reference.
TEST(SAMParallel, IceFallSubstepCountUsesGlobalMaximum)
{
    constexpr int nx = 8;
    constexpr int ny = 4;
    constexpr int nz = 4;
    const amrex::Real dt = amrex::Real(1000.0);

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(8.0, 4.0, 4.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 0)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    const amrex::IntVect max_size(AMREX_D_DECL(2, 2, nz));
    ba.maxSize(max_size);
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab initial_cons(ba, dm, RhoQ6_comp + 1, 0);
    fill_icefall_reduction_conserved_state(initial_cons, geom);

    amrex::MultiFab public_cons(ba, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab::Copy(public_cons, initial_cons, 0, 0, RhoQ6_comp + 1, 0);

    amrex::MultiFab reference_cons(ba, dm, RhoQ6_comp + 1, 0);
    amrex::MultiFab::Copy(reference_cons, initial_cons, 0, 0, RhoQ6_comp + 1, 0);

    const IceFallFluxReductionDiagnostics diagnostics =
        compute_icefall_flux_reduction_diagnostics(geom, initial_cons, dt);

    const int rank = amrex::ParallelDescriptor::MyProc();
    const int nprocs = amrex::ParallelDescriptor::NProcs();
    const amrex::Real dz = geom.CellSize(2);
    const amrex::Real expected_global_max = diagnostics.global_max_flux;
    const int expected_n_substep = sam_substep_count_from_reduced_flux(
        expected_global_max + std::numeric_limits<amrex::Real>::epsilon(),
        dt,
        dz);

    int any_local_differs = (diagnostics.local_n_substep != expected_n_substep) ? 1 : 0;
    amrex::ParallelDescriptor::ReduceIntMax(any_local_differs);

    ASSERT_GT(expected_n_substep, 1)
        << "rank=" << rank
        << " boxes=" << diagnostics.box_summary
        << " max_size=(" << max_size[0] << "," << max_size[1] << "," << max_size[2] << ")"
        << " local_max_flux=" << diagnostics.local_max_flux
        << " global_max_flux=" << diagnostics.global_max_flux
        << " local_n_substep=" << diagnostics.local_n_substep
        << " global_n_substep=" << diagnostics.global_n_substep
        << " expected_n_substep=" << expected_n_substep
        << " any_rank_local_count_differs=" << any_local_differs
        << " dt=" << dt
        << " dz=" << dz
        << " initial_qci_mass=" << diagnostics.initial_qci_mass;

    EXPECT_EQ(diagnostics.global_n_substep, expected_n_substep)
        << "rank=" << rank
        << " owned_boxes=" << diagnostics.owned_box_count
        << " boxes=" << diagnostics.box_summary
        << " max_size=(" << max_size[0] << "," << max_size[1] << "," << max_size[2] << ")"
        << " local_max_flux=" << diagnostics.local_max_flux
        << " global_max_flux=" << diagnostics.global_max_flux
        << " local_n_substep=" << diagnostics.local_n_substep
        << " global_n_substep=" << diagnostics.global_n_substep
        << " expected_n_substep=" << expected_n_substep
        << " any_rank_local_count_differs=" << any_local_differs
        << " dt=" << dt
        << " dz=" << dz
        << " initial_qci_mass=" << diagnostics.initial_qci_mass;

    if (nprocs > 1) {
        ASSERT_EQ(any_local_differs, 1)
            << "rank=" << rank
            << " owned_boxes=" << diagnostics.owned_box_count
            << " boxes=" << diagnostics.box_summary
            << " max_size=(" << max_size[0] << "," << max_size[1] << "," << max_size[2] << ")"
            << " local_max_flux=" << diagnostics.local_max_flux
            << " global_max_flux=" << diagnostics.global_max_flux
            << " local_n_substep=" << diagnostics.local_n_substep
            << " global_n_substep=" << diagnostics.global_n_substep
            << " expected_n_substep=" << expected_n_substep
            << " any_rank_local_count_differs=" << any_local_differs
            << " dt=" << dt
            << " dz=" << dz
            << " initial_qci_mass=" << diagnostics.initial_qci_mass;
    }

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;

    SolverChoice sc = make_sam_solver_choice(MoistureType::SAM);

    SAM public_sam;
    public_sam.Define(sc);
    public_sam.Set_dzmin(dz);
    public_sam.Set_RealWidth(0);
    public_sam.Init(public_cons, ba, geom, dt, z_phys_nd, detJ_cc);
    public_sam.Copy_State_to_Micro(public_cons);
    public_sam.IceFall(sc);
    public_sam.Copy_Micro_to_State(public_cons);

    compute_icefall_reference_global_substep(geom, reference_cons, expected_n_substep, dt);

    const int num_terms = nx * ny * nz;
    const amrex::Real public_final_qci_mass = sum_component_mass(geom, public_cons, RhoQ3_comp);
    const amrex::Real reference_final_qci_mass = sum_component_mass(geom, reference_cons, RhoQ3_comp);
    const amrex::Real mass_residual = public_final_qci_mass - reference_final_qci_mass;

    EXPECT_NEAR(public_final_qci_mass,
                reference_final_qci_mass,
                mpi_reduction_tol(reference_final_qci_mass, num_terms))
        << "rank=" << rank
        << " owned_boxes=" << diagnostics.owned_box_count
        << " boxes=" << diagnostics.box_summary
        << " max_size=(" << max_size[0] << "," << max_size[1] << "," << max_size[2] << ")"
        << " local_max_flux=" << diagnostics.local_max_flux
        << " global_max_flux=" << diagnostics.global_max_flux
        << " local_n_substep=" << diagnostics.local_n_substep
        << " global_n_substep=" << diagnostics.global_n_substep
        << " expected_n_substep=" << expected_n_substep
        << " any_rank_local_count_differs=" << any_local_differs
        << " dt=" << dt
        << " dz=" << dz
        << " initial_qci_mass=" << diagnostics.initial_qci_mass
        << " public_final_qci_mass=" << public_final_qci_mass
        << " reference_final_qci_mass=" << reference_final_qci_mass
        << " mass_residual=" << mass_residual;

    bool local_mismatch = false;
    std::string first_mismatch_message;

    for (amrex::MFIter mfi(public_cons); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto public_arr = public_cons.const_array(mfi);
        const auto reference_arr = reference_cons.const_array(mfi);

        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    const SAMCellState public_state = sam_cons_to_primitive(
                        public_arr(i,j,k,Rho_comp),
                        public_arr(i,j,k,RhoTheta_comp),
                        public_arr(i,j,k,RhoQ1_comp),
                        public_arr(i,j,k,RhoQ2_comp),
                        public_arr(i,j,k,RhoQ3_comp),
                        public_arr(i,j,k,RhoQ4_comp),
                        public_arr(i,j,k,RhoQ5_comp),
                        public_arr(i,j,k,RhoQ6_comp));
                    const SAMCellState reference_state = sam_cons_to_primitive(
                        reference_arr(i,j,k,Rho_comp),
                        reference_arr(i,j,k,RhoTheta_comp),
                        reference_arr(i,j,k,RhoQ1_comp),
                        reference_arr(i,j,k,RhoQ2_comp),
                        reference_arr(i,j,k,RhoQ3_comp),
                        reference_arr(i,j,k,RhoQ4_comp),
                        reference_arr(i,j,k,RhoQ5_comp),
                        reference_arr(i,j,k,RhoQ6_comp));

                    const amrex::Real raw_qci_public = public_arr(i,j,k,RhoQ3_comp);
                    const amrex::Real raw_qci_reference = reference_arr(i,j,k,RhoQ3_comp);
                    const amrex::Real raw_qci_tol = formula_abs_tol(
                        std::max({std::abs(raw_qci_public), std::abs(raw_qci_reference), amrex::Real(1.0)}));
                    const amrex::Real state_tol = formula_abs_tol(
                        std::max({std::abs(public_state.tabs),
                                  std::abs(reference_state.tabs),
                                  std::abs(public_state.qt),
                                  std::abs(reference_state.qt),
                                  amrex::Real(1.0)}));

                    const bool cell_mismatch =
                        std::abs(raw_qci_public - raw_qci_reference) > raw_qci_tol ||
                        std::abs(public_state.rho - reference_state.rho) > state_tol ||
                        std::abs(public_state.theta - reference_state.theta) > state_tol ||
                        std::abs(public_state.tabs - reference_state.tabs) > state_tol ||
                        std::abs(public_state.pres_mbar - reference_state.pres_mbar) > state_tol ||
                        std::abs(public_state.qv - reference_state.qv) > state_tol ||
                        std::abs(public_state.qcl - reference_state.qcl) > state_tol ||
                        std::abs(public_state.qci - reference_state.qci) > state_tol ||
                        std::abs(public_state.qn - reference_state.qn) > state_tol ||
                        std::abs(public_state.qt - reference_state.qt) > state_tol ||
                        std::abs(public_state.qpr - reference_state.qpr) > state_tol ||
                        std::abs(public_state.qps - reference_state.qps) > state_tol ||
                        std::abs(public_state.qpg - reference_state.qpg) > state_tol ||
                        std::abs(public_state.qp - reference_state.qp) > state_tol;

                    if (cell_mismatch) {
                        if (!local_mismatch) {
                            std::ostringstream mismatch_stream;
                            mismatch_stream << "rank=" << rank
                                            << " box=" << bx
                                            << " cell=(" << i << "," << j << "," << k << ")"
                                            << " max_size=(" << max_size[0] << "," << max_size[1] << "," << max_size[2] << ")"
                                            << " local_max_flux=" << diagnostics.local_max_flux
                                            << " global_max_flux=" << diagnostics.global_max_flux
                                            << " local_n_substep=" << diagnostics.local_n_substep
                                            << " global_n_substep=" << diagnostics.global_n_substep
                                            << " expected_n_substep=" << expected_n_substep
                                            << " any_rank_local_count_differs=" << any_local_differs
                                            << " dt=" << dt
                                            << " dz=" << dz
                                            << " initial_qci_mass=" << diagnostics.initial_qci_mass
                                            << " public_final_qci_mass=" << public_final_qci_mass
                                            << " reference_final_qci_mass=" << reference_final_qci_mass
                                            << " mass_residual=" << mass_residual
                                            << icefall_state_summary(" public", public_state, raw_qci_public)
                                            << icefall_state_summary(" reference", reference_state, raw_qci_reference);
                            first_mismatch_message = mismatch_stream.str();
                        }
                        local_mismatch = true;
                    }
                }
            }
        }
    }

    EXPECT_FALSE(local_mismatch) << first_mismatch_message;
}

// Motivation:
// In a one-level column, equal top and bottom component face fluxes leave the
// in-column precip mass unchanged while the bottom diagnostic accumulation is
// updated from the limited bottom component flux. This MPI test checks
// rank-invariant public-path accumulation wiring. The closed multi-level
// column budget is tested in the SAM physical-property tests.
TEST(SAMParallel, PrecipFallGlobalMassAndSurfaceAccumulation)
{
    constexpr int nx = 8;
    constexpr int ny = 4;
    constexpr int nz = 1;
    const amrex::Real pressure_pa = amrex::Real(9.0e4);
    const amrex::Real dt = amrex::Real(0.1);
    const amrex::Real rho0 = amrex::Real(1.29);

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(1.0, 1.0, 1.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 0)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    ba.maxSize(amrex::IntVect(4, 2, nz));
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    fill_precipfall_conserved_state(cons, pressure_pa);

    const amrex::Real dx = geom.CellSize(0);
    const amrex::Real dy = geom.CellSize(1);
    const amrex::Real dz = geom.CellSize(2);
    const amrex::Real cell_volume = dx * dy * dz;
    const amrex::Real cell_area = dx * dy;
    const auto terminal_velocities = precipfall_terminal_velocities();

    amrex::Real expected_rain_accum_mass = amrex::Real(0.0);
    amrex::Real expected_snow_accum_mass = amrex::Real(0.0);
    amrex::Real expected_graup_accum_mass = amrex::Real(0.0);
    const amrex::Real coef = dt / dz;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const PrecipFallCellState state = precipfall_cell_state(i, j, pressure_pa);
            const SAMPrecipComponentFaceState face_state{
                state.rho, state.tabs, state.qpr, state.qps, state.qpg};
            const SAMPrecipFluxComponents corrected_fluxes =
                sam_precip_flux_components_density_corrected(
                    sam_precip_component_fluxes_from_face_state(face_state,
                                                                terminal_velocities[0],
                                                                terminal_velocities[1],
                                                                terminal_velocities[2]),
                    rho0,
                    face_state.rho_avg);
            const SAMPrecipFluxComponents limited_fluxes{
                sam_limit_precip_component_flux(corrected_fluxes.rain,
                                                state.rho,
                                                state.qpr,
                                                amrex::Real(1.0),
                                                coef),
                sam_limit_precip_component_flux(corrected_fluxes.snow,
                                                state.rho,
                                                state.qps,
                                                amrex::Real(1.0),
                                                coef),
                sam_limit_precip_component_flux(corrected_fluxes.graupel,
                                                state.rho,
                                                state.qpg,
                                                amrex::Real(1.0),
                                                coef)};

            const SAMSurfaceAccumulation expected_accum =
                sam_surface_accumulation_from_component_fluxes(limited_fluxes, dt);

            expected_rain_accum_mass += expected_accum.rain * cell_area * rhor / amrex::Real(1000.0);
            expected_snow_accum_mass += expected_accum.snow * cell_area * rhos / amrex::Real(1000.0);
            expected_graup_accum_mass += expected_accum.graupel * cell_area * rhog / amrex::Real(1000.0);
        }
    }
    const amrex::Real expected_surface_accum_mass =
        expected_rain_accum_mass + expected_snow_accum_mass + expected_graup_accum_mass;

    const amrex::Real initial_precip_mass =
        (cons.sum(RhoQ4_comp) + cons.sum(RhoQ5_comp) + cons.sum(RhoQ6_comp)) * cell_volume;

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;

    SolverChoice sc = make_sam_solver_choice();
    SAM sam;
    sam.Define(sc);
    sam.Set_dzmin(dz);
    sam.Set_RealWidth(0);
    sam.Init(cons, ba, geom, dt, z_phys_nd, detJ_cc);
    sam.Copy_State_to_Micro(cons);
    sam.PrecipFall(sc);
    sam.Copy_Micro_to_State(cons);

    const amrex::Real final_precip_mass =
        (cons.sum(RhoQ4_comp) + cons.sum(RhoQ5_comp) + cons.sum(RhoQ6_comp)) * cell_volume;

    const amrex::Real rain_accum_mass = sam.Qmoist_Ptr(0)->sum(0) * cell_area * rhor / amrex::Real(1000.0);
    const amrex::Real snow_accum_mass = sam.Qmoist_Ptr(1)->sum(0) * cell_area * rhos / amrex::Real(1000.0);
    const amrex::Real graup_accum_mass = sam.Qmoist_Ptr(2)->sum(0) * cell_area * rhog / amrex::Real(1000.0);
    const amrex::Real surface_accum_mass = rain_accum_mass + snow_accum_mass + graup_accum_mass;

    const int num_terms = nx * ny;
    const int rank = amrex::ParallelDescriptor::MyProc();

    EXPECT_NEAR(final_precip_mass,
                initial_precip_mass,
                mpi_reduction_tol(initial_precip_mass, num_terms))
        << "rank=" << rank
        << " initial_precip_mass=" << initial_precip_mass
        << " final_precip_mass=" << final_precip_mass;
    EXPECT_NEAR(rain_accum_mass,
                expected_rain_accum_mass,
                mpi_reduction_tol(expected_rain_accum_mass, num_terms))
        << "rank=" << rank
        << " rain_accum_mass=" << rain_accum_mass
        << " expected_rain_accum_mass=" << expected_rain_accum_mass;
    EXPECT_NEAR(snow_accum_mass,
                expected_snow_accum_mass,
                mpi_reduction_tol(expected_snow_accum_mass, num_terms))
        << "rank=" << rank
        << " snow_accum_mass=" << snow_accum_mass
        << " expected_snow_accum_mass=" << expected_snow_accum_mass;
    EXPECT_NEAR(graup_accum_mass,
                expected_graup_accum_mass,
                mpi_reduction_tol(expected_graup_accum_mass, num_terms))
        << "rank=" << rank
        << " graup_accum_mass=" << graup_accum_mass
        << " expected_graup_accum_mass=" << expected_graup_accum_mass;
    EXPECT_NEAR(surface_accum_mass,
                expected_surface_accum_mass,
                mpi_reduction_tol(expected_surface_accum_mass, num_terms))
        << "rank=" << rank
        << " surface_accum_mass=" << surface_accum_mass
        << " expected_surface_accum_mass=" << expected_surface_accum_mass
        << " rain_accum_mass=" << rain_accum_mass
        << " snow_accum_mass=" << snow_accum_mass
        << " graup_accum_mass=" << graup_accum_mass;
}

// Motivation:
// The public SAM Precip path is cell-local except for the precomputed
// coefficient rows and MPI reductions used by the test budget checks. Running
// this same contract under 1, 2, and 4 ranks should preserve the global total
// water budget and leave the surface accumulation reservoirs untouched,
// independent of box ownership and reduction order.
//
// This public-path MPI test intentionally checks total water and zero surface
// accumulation only. A latent-proxy closure check for
// Copy_State_to_Micro -> Compute_Coefficients -> Precip -> Copy_Micro_to_State
// is a separate thermodynamic roundtrip follow-up because this path includes
// EOS reconstruction, pressure conversion, clipping, and plane-averaged
// coefficient rows around the scalar source helper.
TEST(SAMParallel, PrecipGlobalWaterBudgetDecompositionInvariant)
{
    constexpr int nx = 8;
    constexpr int ny = 4;
    constexpr int nz = 4;
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real dt = amrex::Real(0.1);

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(1.0, 1.0, 1.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    ba.maxSize(amrex::IntVect(4, 2, 2));
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    fill_precip_budget_conserved_state(cons, MoistureType::SAM, pres_mbar);

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;

    SolverChoice sc = make_sam_solver_choice(MoistureType::SAM);
    SAM sam;
    sam.Define(sc);
    sam.Set_dzmin(geom.CellSize(2));
    sam.Set_RealWidth(0);
    sam.Init(cons, ba, geom, dt, z_phys_nd, detJ_cc);

    const SAMBudget initial_budget = compute_sam_budget(geom, cons, sam);

    sam.Copy_State_to_Micro(cons);
    sam.Compute_Coefficients();
    sam.Precip(sc);
    sam.Copy_Micro_to_State(cons);

    const SAMBudget final_budget = compute_sam_budget(geom, cons, sam);
    const int num_terms = nx * ny * nz;
    const int rank = amrex::ParallelDescriptor::MyProc();

    EXPECT_NEAR(final_budget.total_water_mass,
                initial_budget.total_water_mass,
                mpi_reduction_tol(initial_budget.total_water_mass, num_terms))
        << "rank=" << rank
        << " initial_total_water_mass=" << initial_budget.total_water_mass
        << " final_total_water_mass=" << final_budget.total_water_mass;
    EXPECT_NEAR(final_budget.rain_accum_mass, amrex::Real(0.0), exact_zero_tol())
        << "rank=" << rank
        << " rain_accum_mass=" << final_budget.rain_accum_mass;
    EXPECT_NEAR(final_budget.snow_accum_mass, amrex::Real(0.0), exact_zero_tol())
        << "rank=" << rank
        << " snow_accum_mass=" << final_budget.snow_accum_mass;
    EXPECT_NEAR(final_budget.graupel_accum_mass, amrex::Real(0.0), exact_zero_tol())
        << "rank=" << rank
        << " graupel_accum_mass=" << final_budget.graupel_accum_mass;
}

// Motivation:
// The same public Precip water-budget contract should remain decomposition-
// invariant when the budget is checked in terrain-like coordinates with a
// nonuniform detJ.
TEST(SAMParallel, PrecipDetJWeightedWaterBudgetDecompositionInvariant)
{
    constexpr int nx = 8;
    constexpr int ny = 4;
    constexpr int nz = 4;
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real dt = amrex::Real(1.0);

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(1.0, 1.0, 1.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    ba.maxSize(amrex::IntVect(4, 2, 2));
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    fill_precip_budget_conserved_state(cons, MoistureType::SAM, pres_mbar);

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc = std::make_unique<amrex::MultiFab>(ba, dm, 1, 0);
    fill_nonuniform_detj(*detJ_cc);

    SolverChoice sc = make_sam_solver_choice(MoistureType::SAM);
    SAM sam;
    sam.Define(sc);
    sam.Set_dzmin(geom.CellSize(2));
    sam.Set_RealWidth(0);
    sam.Init(cons, ba, geom, dt, z_phys_nd, detJ_cc);

    const SAMBudget initial_budget = compute_sam_budget(geom, cons, sam, detJ_cc.get());

    sam.Copy_State_to_Micro(cons);
    sam.Compute_Coefficients();
    sam.Precip(sc);
    sam.Copy_Micro_to_State(cons);

    const SAMBudget final_budget = compute_sam_budget(geom, cons, sam, detJ_cc.get());
    const int num_terms = nx * ny * nz;
    const int rank = amrex::ParallelDescriptor::MyProc();

    EXPECT_NEAR(final_budget.total_water_mass,
                initial_budget.total_water_mass,
                mpi_reduction_tol(initial_budget.total_water_mass, num_terms))
        << "rank=" << rank
        << " initial_total_water_mass=" << initial_budget.total_water_mass
        << " final_total_water_mass=" << final_budget.total_water_mass;
    EXPECT_NEAR(final_budget.rain_accum_mass, amrex::Real(0.0), exact_zero_tol())
        << "rank=" << rank
        << " rain_accum_mass=" << final_budget.rain_accum_mass;
    EXPECT_NEAR(final_budget.snow_accum_mass, amrex::Real(0.0), exact_zero_tol())
        << "rank=" << rank
        << " snow_accum_mass=" << final_budget.snow_accum_mass;
    EXPECT_NEAR(final_budget.graupel_accum_mass, amrex::Real(0.0), exact_zero_tol())
        << "rank=" << rank
        << " graupel_accum_mass=" << final_budget.graupel_accum_mass;
}

// Motivation:
// The public SAM Precip path should conserve the intended constant-pressure
// latent proxy when evaluated with the initial pressure snapshot held fixed
// across the local microphysics source update. The EOS-projected proxy from
// conserved rho and rho*theta is diagnostic only here.
TEST(SAMParallel, PrecipConstantPressureLatentProxyBudgetDecompositionInvariant)
{
    constexpr int nx = 8;
    constexpr int ny = 4;
    constexpr int nz = 4;
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real dt = amrex::Real(0.1);

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(1.0, 1.0, 1.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    const amrex::IntVect max_size(4, 2, 2);
    ba.maxSize(max_size);
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    fill_precip_budget_conserved_state(cons, MoistureType::SAM, pres_mbar);
    amrex::MultiFab pressure0_mbar(ba, dm, 1, 0);
    fill_pressure_snapshot_from_cons(cons, pressure0_mbar);

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;

    SolverChoice sc = make_sam_solver_choice(MoistureType::SAM);
    SAM sam;
    sam.Define(sc);
    sam.Set_dzmin(geom.CellSize(2));
    sam.Set_RealWidth(0);
    sam.Init(cons, ba, geom, dt, z_phys_nd, detJ_cc);

    const SAMPressureSnapshotBudget initial_budget =
        compute_sam_pressure_snapshot_budget(geom, cons, pressure0_mbar);

    sam.Copy_State_to_Micro(cons);
    sam.Compute_Coefficients();
    sam.Precip(sc);
    sam.Copy_Micro_to_State(cons);

    const SAMPressureSnapshotBudget final_budget =
        compute_sam_pressure_snapshot_budget(geom, cons, pressure0_mbar);
    const int num_terms = nx * ny * nz;
    const int rank = amrex::ParallelDescriptor::MyProc();
    const amrex::Real constant_pressure_residual =
        final_budget.latent_proxy_mass_constant_pressure -
        initial_budget.latent_proxy_mass_constant_pressure;
    const amrex::Real eos_projected_residual =
        final_budget.latent_proxy_mass_eos_projected -
        initial_budget.latent_proxy_mass_eos_projected;
    const amrex::Real total_water_residual =
        final_budget.total_water_mass - initial_budget.total_water_mass;

    EXPECT_NEAR(final_budget.total_water_mass,
                initial_budget.total_water_mass,
                mpi_reduction_tol(initial_budget.total_water_mass, num_terms))
        << "initial_total_water_mass=" << initial_budget.total_water_mass
        << " final_total_water_mass=" << final_budget.total_water_mass
        << " initial_constant_pressure_latent_proxy_mass="
        << initial_budget.latent_proxy_mass_constant_pressure
        << " final_constant_pressure_latent_proxy_mass="
        << final_budget.latent_proxy_mass_constant_pressure
        << " constant_pressure_residual=" << constant_pressure_residual
        << " initial_eos_projected_latent_proxy_mass="
        << initial_budget.latent_proxy_mass_eos_projected
        << " final_eos_projected_latent_proxy_mass="
        << final_budget.latent_proxy_mass_eos_projected
        << " eos_projected_residual=" << eos_projected_residual
        << " min_pressure0_mbar=" << initial_budget.min_pressure0_mbar
        << " max_pressure0_mbar=" << initial_budget.max_pressure0_mbar
        << " dt=" << dt
        << " rank=" << rank
        << " box_decomposition=(" << max_size[0] << "," << max_size[1] << "," << max_size[2] << ")";
    EXPECT_NEAR(final_budget.latent_proxy_mass_constant_pressure,
                initial_budget.latent_proxy_mass_constant_pressure,
                mpi_reduction_tol(initial_budget.latent_proxy_mass_constant_pressure, num_terms))
        << "initial_total_water_mass=" << initial_budget.total_water_mass
        << " final_total_water_mass=" << final_budget.total_water_mass
        << " total_water_residual=" << total_water_residual
        << " initial_constant_pressure_latent_proxy_mass="
        << initial_budget.latent_proxy_mass_constant_pressure
        << " final_constant_pressure_latent_proxy_mass="
        << final_budget.latent_proxy_mass_constant_pressure
        << " constant_pressure_residual=" << constant_pressure_residual
        << " initial_eos_projected_latent_proxy_mass="
        << initial_budget.latent_proxy_mass_eos_projected
        << " final_eos_projected_latent_proxy_mass="
        << final_budget.latent_proxy_mass_eos_projected
        << " eos_projected_residual=" << eos_projected_residual
        << " min_pressure0_mbar=" << initial_budget.min_pressure0_mbar
        << " max_pressure0_mbar=" << initial_budget.max_pressure0_mbar
        << " dt=" << dt
        << " rank=" << rank
        << " box_decomposition=(" << max_size[0] << "," << max_size[1] << "," << max_size[2] << ")";
    EXPECT_NEAR(sum_multifab_scalar(*sam.Qmoist_Ptr(0)), amrex::Real(0.0), exact_zero_tol())
        << "rank=" << rank;
    EXPECT_NEAR(sum_multifab_scalar(*sam.Qmoist_Ptr(1)), amrex::Real(0.0), exact_zero_tol())
        << "rank=" << rank;
    EXPECT_NEAR(sum_multifab_scalar(*sam.Qmoist_Ptr(2)), amrex::Real(0.0), exact_zero_tol())
        << "rank=" << rank;
}

// Motivation:
// The same public Precip constant-pressure latent-proxy contract should remain
// decomposition-invariant when budgets are checked in terrain-like coordinates
// with a nonuniform detJ.
TEST(SAMParallel, PrecipDetJWeightedConstantPressureLatentProxyBudgetDecompositionInvariant)
{
    constexpr int nx = 8;
    constexpr int ny = 4;
    constexpr int nz = 4;
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real dt = amrex::Real(1.0);

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(1.0, 1.0, 1.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    const amrex::IntVect max_size(4, 2, 2);
    ba.maxSize(max_size);
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    fill_precip_budget_conserved_state(cons, MoistureType::SAM, pres_mbar);
    amrex::MultiFab pressure0_mbar(ba, dm, 1, 0);
    fill_pressure_snapshot_from_cons(cons, pressure0_mbar);

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc = std::make_unique<amrex::MultiFab>(ba, dm, 1, 0);
    fill_nonuniform_detj(*detJ_cc);

    SolverChoice sc = make_sam_solver_choice(MoistureType::SAM);
    SAM sam;
    sam.Define(sc);
    sam.Set_dzmin(geom.CellSize(2));
    sam.Set_RealWidth(0);
    sam.Init(cons, ba, geom, dt, z_phys_nd, detJ_cc);

    const SAMPressureSnapshotBudget initial_budget =
        compute_sam_pressure_snapshot_budget(geom, cons, pressure0_mbar, detJ_cc.get());

    sam.Copy_State_to_Micro(cons);
    sam.Compute_Coefficients();
    sam.Precip(sc);
    sam.Copy_Micro_to_State(cons);

    const SAMPressureSnapshotBudget final_budget =
        compute_sam_pressure_snapshot_budget(geom, cons, pressure0_mbar, detJ_cc.get());
    const int num_terms = nx * ny * nz;
    const int rank = amrex::ParallelDescriptor::MyProc();
    const amrex::Real constant_pressure_residual =
        final_budget.latent_proxy_mass_constant_pressure -
        initial_budget.latent_proxy_mass_constant_pressure;
    const amrex::Real eos_projected_residual =
        final_budget.latent_proxy_mass_eos_projected -
        initial_budget.latent_proxy_mass_eos_projected;
    const amrex::Real total_water_residual =
        final_budget.total_water_mass - initial_budget.total_water_mass;
    const amrex::Real detj_min = detJ_cc->min(0);
    const amrex::Real detj_max = detJ_cc->max(0);

    EXPECT_NEAR(final_budget.total_water_mass,
                initial_budget.total_water_mass,
                mpi_reduction_tol(initial_budget.total_water_mass, num_terms))
        << "initial_total_water_mass=" << initial_budget.total_water_mass
        << " final_total_water_mass=" << final_budget.total_water_mass
        << " initial_constant_pressure_latent_proxy_mass="
        << initial_budget.latent_proxy_mass_constant_pressure
        << " final_constant_pressure_latent_proxy_mass="
        << final_budget.latent_proxy_mass_constant_pressure
        << " constant_pressure_residual=" << constant_pressure_residual
        << " initial_eos_projected_latent_proxy_mass="
        << initial_budget.latent_proxy_mass_eos_projected
        << " final_eos_projected_latent_proxy_mass="
        << final_budget.latent_proxy_mass_eos_projected
        << " eos_projected_residual=" << eos_projected_residual
        << " min_pressure0_mbar=" << initial_budget.min_pressure0_mbar
        << " max_pressure0_mbar=" << initial_budget.max_pressure0_mbar
        << " dt=" << dt
        << " rank=" << rank
        << " box_decomposition=(" << max_size[0] << "," << max_size[1] << "," << max_size[2] << ")"
        << " detJ_min=" << detj_min
        << " detJ_max=" << detj_max;
    EXPECT_NEAR(final_budget.latent_proxy_mass_constant_pressure,
                initial_budget.latent_proxy_mass_constant_pressure,
                mpi_reduction_tol(initial_budget.latent_proxy_mass_constant_pressure, num_terms))
        << "initial_total_water_mass=" << initial_budget.total_water_mass
        << " final_total_water_mass=" << final_budget.total_water_mass
        << " total_water_residual=" << total_water_residual
        << " initial_constant_pressure_latent_proxy_mass="
        << initial_budget.latent_proxy_mass_constant_pressure
        << " final_constant_pressure_latent_proxy_mass="
        << final_budget.latent_proxy_mass_constant_pressure
        << " constant_pressure_residual=" << constant_pressure_residual
        << " initial_eos_projected_latent_proxy_mass="
        << initial_budget.latent_proxy_mass_eos_projected
        << " final_eos_projected_latent_proxy_mass="
        << final_budget.latent_proxy_mass_eos_projected
        << " eos_projected_residual=" << eos_projected_residual
        << " min_pressure0_mbar=" << initial_budget.min_pressure0_mbar
        << " max_pressure0_mbar=" << initial_budget.max_pressure0_mbar
        << " dt=" << dt
        << " rank=" << rank
        << " box_decomposition=(" << max_size[0] << "," << max_size[1] << "," << max_size[2] << ")"
        << " detJ_min=" << detj_min
        << " detJ_max=" << detj_max;
    EXPECT_NEAR(sum_multifab_scalar(*sam.Qmoist_Ptr(0)), amrex::Real(0.0), exact_zero_tol())
        << "rank=" << rank
        << " detJ_min=" << detj_min
        << " detJ_max=" << detj_max;
    EXPECT_NEAR(sum_multifab_scalar(*sam.Qmoist_Ptr(1)), amrex::Real(0.0), exact_zero_tol())
        << "rank=" << rank
        << " detJ_min=" << detj_min
        << " detJ_max=" << detj_max;
    EXPECT_NEAR(sum_multifab_scalar(*sam.Qmoist_Ptr(2)), amrex::Real(0.0), exact_zero_tol())
        << "rank=" << rank
        << " detJ_min=" << detj_min
        << " detJ_max=" << detj_max;
}

// Motivation:
// The public PrecipFall path should conserve rain, snow, and graupel budgets
// independently under multiple x-y decompositions and MPI rank counts, while
// keeping each test box unsplit in z.
TEST(SAMParallel, PrecipFallComponentBudgetsIndependentOfDecomposition)
{
    const std::vector<amrex::IntVect> max_sizes = {
        amrex::IntVect(8, 4, 4),
        amrex::IntVect(4, 2, 4),
        amrex::IntVect(2, 1, 4)};

    const int rank = amrex::ParallelDescriptor::MyProc();
    const int num_terms = 8 * 4 * 4;
    std::vector<PrecipFallComponentBudgetResult> results;
    results.reserve(max_sizes.size());

    for (const amrex::IntVect& max_size : max_sizes) {
        results.push_back(run_precipfall_component_budget_case(max_size));
        const PrecipFallComponentBudgetResult& result = results.back();

        EXPECT_NEAR(result.rain_residual,
                    amrex::Real(0.0),
                    mpi_reduction_tol(amrex::Real(1.0), num_terms))
            << "rank=" << rank
            << " max_size=" << result.max_size
            << " rain_residual=" << result.rain_residual
            << " snow_residual=" << result.snow_residual
            << " graupel_residual=" << result.graupel_residual
            << " total_residual=" << result.total_residual
            << " rain_accum_mass=" << result.rain_accum_mass
            << " snow_accum_mass=" << result.snow_accum_mass
            << " graupel_accum_mass=" << result.graupel_accum_mass
            << " detj_min=" << result.detj_min
            << " detj_max=" << result.detj_max;
        EXPECT_NEAR(result.snow_residual,
                    amrex::Real(0.0),
                    mpi_reduction_tol(amrex::Real(1.0), num_terms));
        EXPECT_NEAR(result.graupel_residual,
                    amrex::Real(0.0),
                    mpi_reduction_tol(amrex::Real(1.0), num_terms));
        EXPECT_NEAR(result.total_residual,
                    amrex::Real(0.0),
                    mpi_reduction_tol(amrex::Real(1.0), 3 * num_terms));
    }

    ASSERT_FALSE(results.empty());
    const PrecipFallComponentBudgetResult& baseline = results.front();
    for (std::size_t idx = 1; idx < results.size(); ++idx) {
        const PrecipFallComponentBudgetResult& result = results[idx];
        EXPECT_NEAR(result.rain_accum_mass,
                    baseline.rain_accum_mass,
                    mpi_reduction_tol(baseline.rain_accum_mass, num_terms))
            << "rank=" << rank
            << " baseline_max_size=" << baseline.max_size
            << " compare_max_size=" << result.max_size;
        EXPECT_NEAR(result.snow_accum_mass,
                    baseline.snow_accum_mass,
                    mpi_reduction_tol(baseline.snow_accum_mass, num_terms));
        EXPECT_NEAR(result.graupel_accum_mass,
                    baseline.graupel_accum_mass,
                    mpi_reduction_tol(baseline.graupel_accum_mass, num_terms));
    }
}
#include <array>

#include <gtest/gtest.h>

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

void expect_nonnegative_species (const SAMCellState& state)
{
    EXPECT_GE(state.qv, -exact_zero_or_near_zero_tol());
    EXPECT_GE(state.qcl, -exact_zero_or_near_zero_tol());
    EXPECT_GE(state.qci, -exact_zero_or_near_zero_tol());
    EXPECT_GE(state.qpr, -exact_zero_or_near_zero_tol());
    EXPECT_GE(state.qps, -exact_zero_or_near_zero_tol());
    EXPECT_GE(state.qpg, -exact_zero_or_near_zero_tol());
}

amrex::Real mixed_phase_tabs ()
{
    const amrex::Real lower = std::max(tprmin, tgrmin) + amrex::Real(1.0e-3);
    const amrex::Real upper = std::min(tprmax, tgrmax) - amrex::Real(1.0e-3);
    EXPECT_LT(lower, upper);
    return amrex::Real(0.5) * (lower + upper);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real mixed_qsat_for_state (const int sam_mode,
                                  const amrex::Real tabs,
                                  const amrex::Real pres_mbar) noexcept
{
    amrex::Real qsatw;
    amrex::Real qsati;
    erf_qsatw(tabs, pres_mbar, qsatw);
    erf_qsati(tabs, pres_mbar, qsati);
    const amrex::Real omn = sam_cloud_liquid_fraction(sam_mode, tabs, a_bg, tbgmin * a_bg);
    return sam_mixed_qsat(omn, qsatw, qsati);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SAMCellState make_cell_state (const amrex::Real tabs,
                              const amrex::Real pres_mbar,
                              const amrex::Real qv,
                              const amrex::Real qcl,
                              const amrex::Real qci,
                              const amrex::Real qpr,
                              const amrex::Real qps,
                              const amrex::Real qpg) noexcept
{
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

SAMPrecipConfig make_precip_config (const int sam_mode,
                                    const bool enable_precip = true,
                                    const amrex::Real dtn = amrex::Real(1.0))
{
    return {
        sam_mode,
        enable_precip,
        dtn,
        kRdOcp,
        kFacCond,
        kFacFus,
        kFacSub,
        std::numeric_limits<amrex::Real>::epsilon(),
        (three + b_rain) / amrex::Real(4.0),
        (three + b_snow) / amrex::Real(4.0),
        (three + b_grau) / amrex::Real(4.0),
        (amrex::Real(5.0) + b_rain) / amrex::Real(8.0),
        (amrex::Real(5.0) + b_snow) / amrex::Real(8.0),
        (amrex::Real(5.0) + b_grau) / amrex::Real(8.0)};
}

SAMCoefficientRow make_coeffs (const amrex::Real scale = amrex::Real(1.0))
{
    return {
        amrex::Real(2.0) * scale,
        amrex::Real(2.5) * scale,
        amrex::Real(2.2) * scale,
        amrex::Real(1.8) * scale,
        amrex::Real(1.1) * scale,
        amrex::Real(1.4) * scale,
        amrex::Real(2.4) * scale,
        amrex::Real(2.1) * scale,
        amrex::Real(1.2) * scale,
        amrex::Real(1.6) * scale,
        amrex::Real(1.3) * scale,
        amrex::Real(1.7) * scale};
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
    return amrex::get<0>(reduced);
}

amrex::Real sum_accum_mass (const amrex::Geometry& geom,
                            const amrex::MultiFab& accum,
                            const amrex::Real precip_density)
{
    const amrex::Real cell_area = geom.CellSize(0) * geom.CellSize(1);
    return accum.sum(0) * cell_area * precip_density / amrex::Real(1000.0);
}

struct PrecipFallRunDiagnostics {
    std::array<amrex::Real, 3> max_raw_fluxes{{amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)}};
    std::array<amrex::Real, 3> max_limited_fluxes{{amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)}};
    int n_substep{1};
    bool limiter_active{false};
    amrex::Real detj_min{amrex::Real(1.0)};
    amrex::Real detj_max{amrex::Real(1.0)};
};

struct PrecipFallBudgetSummary {
    amrex::Real initial_rain_mass{amrex::Real(0.0)};
    amrex::Real final_rain_mass{amrex::Real(0.0)};
    amrex::Real rain_accum_mass{amrex::Real(0.0)};
    amrex::Real initial_snow_mass{amrex::Real(0.0)};
    amrex::Real final_snow_mass{amrex::Real(0.0)};
    amrex::Real snow_accum_mass{amrex::Real(0.0)};
    amrex::Real initial_graupel_mass{amrex::Real(0.0)};
    amrex::Real final_graupel_mass{amrex::Real(0.0)};
    amrex::Real graupel_accum_mass{amrex::Real(0.0)};
    amrex::Real min_qpr_after{amrex::Real(0.0)};
    amrex::Real min_qps_after{amrex::Real(0.0)};
    amrex::Real min_qpg_after{amrex::Real(0.0)};
    PrecipFallRunDiagnostics diagnostics{};
};

PrecipFallRunDiagnostics compute_precipfall_run_diagnostics (const amrex::Geometry& geom,
                                                             const amrex::MultiFab& cons,
                                                             const amrex::Real dt,
                                                             const amrex::MultiFab* detJ_cc = nullptr)
{
    const amrex::Real rho0 = amrex::Real(1.29);
    const auto terminal_velocities = precipfall_terminal_velocities();

    PrecipFallRunDiagnostics diagnostics{};
    amrex::Real max_total_raw_flux = amrex::Real(0.0);

    for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto arr = cons.const_array(mfi);
        const auto detj = (detJ_cc != nullptr) ? detJ_cc->const_array(mfi) : amrex::Array4<const amrex::Real>{};
        const int k_lo = bx.smallEnd(2);
        const int k_hi = bx.bigEnd(2);

        for (int k = k_lo; k <= k_hi + 1; ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    const int donor_k = sam_precip_face_donor_k(k, k_lo, k_hi);
                    const SAMPrimitiveCell donor = sam_cons_to_primitive(
                        arr(i,j,donor_k,Rho_comp),
                        arr(i,j,donor_k,RhoTheta_comp),
                        arr(i,j,donor_k,RhoQ1_comp),
                        arr(i,j,donor_k,RhoQ2_comp),
                        arr(i,j,donor_k,RhoQ3_comp),
                        arr(i,j,donor_k,RhoQ4_comp),
                        arr(i,j,donor_k,RhoQ5_comp),
                        arr(i,j,donor_k,RhoQ6_comp));

                    SAMPrimitiveCell upper = donor;
                    SAMPrimitiveCell lower = donor;
                    if (k > k_lo && k <= k_hi) {
                        lower = sam_cons_to_primitive(
                            arr(i,j,k-1,Rho_comp),
                            arr(i,j,k-1,RhoTheta_comp),
                            arr(i,j,k-1,RhoQ1_comp),
                            arr(i,j,k-1,RhoQ2_comp),
                            arr(i,j,k-1,RhoQ3_comp),
                            arr(i,j,k-1,RhoQ4_comp),
                            arr(i,j,k-1,RhoQ5_comp),
                            arr(i,j,k-1,RhoQ6_comp));
                        upper = sam_cons_to_primitive(
                            arr(i,j,k,Rho_comp),
                            arr(i,j,k,RhoTheta_comp),
                            arr(i,j,k,RhoQ1_comp),
                            arr(i,j,k,RhoQ2_comp),
                            arr(i,j,k,RhoQ3_comp),
                            arr(i,j,k,RhoQ4_comp),
                            arr(i,j,k,RhoQ5_comp),
                            arr(i,j,k,RhoQ6_comp));
                    }

                    const SAMPrecipComponentFaceState face_state =
                        (k == k_lo) ? SAMPrecipComponentFaceState{donor.rho, donor.tabs, donor.qpr, donor.qps, donor.qpg}
                                    : (k == k_hi + 1) ? SAMPrecipComponentFaceState{donor.rho, donor.tabs, donor.qpr, donor.qps, donor.qpg}
                                                      : SAMPrecipComponentFaceState{
                                                            myhalf * (lower.rho + upper.rho),
                                                            myhalf * (lower.tabs + upper.tabs),
                                                            myhalf * (lower.qpr + upper.qpr),
                                                            myhalf * (lower.qps + upper.qps),
                                                            myhalf * (lower.qpg + upper.qpg)};

                    const SAMPrecipFluxComponents raw_fluxes =
                        sam_precip_component_fluxes_from_face_state(face_state,
                                                                    terminal_velocities[0],
                                                                    terminal_velocities[1],
                                                                    terminal_velocities[2]);
                    const SAMPrecipFluxComponents corrected_fluxes =
                        sam_precip_flux_components_density_corrected(raw_fluxes, rho0, face_state.rho_avg);

                    diagnostics.max_raw_fluxes[0] = std::max(diagnostics.max_raw_fluxes[0], corrected_fluxes.rain);
                    diagnostics.max_raw_fluxes[1] = std::max(diagnostics.max_raw_fluxes[1], corrected_fluxes.snow);
                    diagnostics.max_raw_fluxes[2] = std::max(diagnostics.max_raw_fluxes[2], corrected_fluxes.graupel);
                    max_total_raw_flux = std::max(max_total_raw_flux,
                                                  corrected_fluxes.rain + corrected_fluxes.snow + corrected_fluxes.graupel);

                    const amrex::Real detj_donor = detj ? detj(i,j,donor_k) : amrex::Real(1.0);
                    diagnostics.detj_min = std::min(diagnostics.detj_min, detj_donor);
                    diagnostics.detj_max = std::max(diagnostics.detj_max, detj_donor);
                }
            }
        }
    }

    diagnostics.n_substep = sam_substep_count_from_reduced_flux(
        max_total_raw_flux + std::numeric_limits<amrex::Real>::epsilon(),
        dt,
        geom.CellSize(2));

    const amrex::Real coef_substep = (dt / diagnostics.n_substep) / geom.CellSize(2);
    for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto arr = cons.const_array(mfi);
        const auto detj = (detJ_cc != nullptr) ? detJ_cc->const_array(mfi) : amrex::Array4<const amrex::Real>{};
        const int k_lo = bx.smallEnd(2);
        const int k_hi = bx.bigEnd(2);

        for (int k = k_lo; k <= k_hi + 1; ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    const int donor_k = sam_precip_face_donor_k(k, k_lo, k_hi);
                    const SAMPrimitiveCell donor = sam_cons_to_primitive(
                        arr(i,j,donor_k,Rho_comp),
                        arr(i,j,donor_k,RhoTheta_comp),
                        arr(i,j,donor_k,RhoQ1_comp),
                        arr(i,j,donor_k,RhoQ2_comp),
                        arr(i,j,donor_k,RhoQ3_comp),
                        arr(i,j,donor_k,RhoQ4_comp),
                        arr(i,j,donor_k,RhoQ5_comp),
                        arr(i,j,donor_k,RhoQ6_comp));

                    SAMPrimitiveCell upper = donor;
                    SAMPrimitiveCell lower = donor;
                    if (k > k_lo && k <= k_hi) {
                        lower = sam_cons_to_primitive(
                            arr(i,j,k-1,Rho_comp),
                            arr(i,j,k-1,RhoTheta_comp),
                            arr(i,j,k-1,RhoQ1_comp),
                            arr(i,j,k-1,RhoQ2_comp),
                            arr(i,j,k-1,RhoQ3_comp),
                            arr(i,j,k-1,RhoQ4_comp),
                            arr(i,j,k-1,RhoQ5_comp),
                            arr(i,j,k-1,RhoQ6_comp));
                        upper = sam_cons_to_primitive(
                            arr(i,j,k,Rho_comp),
                            arr(i,j,k,RhoTheta_comp),
                            arr(i,j,k,RhoQ1_comp),
                            arr(i,j,k,RhoQ2_comp),
                            arr(i,j,k,RhoQ3_comp),
                            arr(i,j,k,RhoQ4_comp),
                            arr(i,j,k,RhoQ5_comp),
                            arr(i,j,k,RhoQ6_comp));
                    }

                    const SAMPrecipComponentFaceState face_state =
                        (k == k_lo) ? SAMPrecipComponentFaceState{donor.rho, donor.tabs, donor.qpr, donor.qps, donor.qpg}
                                    : (k == k_hi + 1) ? SAMPrecipComponentFaceState{donor.rho, donor.tabs, donor.qpr, donor.qps, donor.qpg}
                                                      : SAMPrecipComponentFaceState{
                                                            myhalf * (lower.rho + upper.rho),
                                                            myhalf * (lower.tabs + upper.tabs),
                                                            myhalf * (lower.qpr + upper.qpr),
                                                            myhalf * (lower.qps + upper.qps),
                                                            myhalf * (lower.qpg + upper.qpg)};
                    const SAMPrecipFluxComponents corrected_fluxes =
                        sam_precip_flux_components_density_corrected(
                            sam_precip_component_fluxes_from_face_state(face_state,
                                                                        terminal_velocities[0],
                                                                        terminal_velocities[1],
                                                                        terminal_velocities[2]),
                            rho0,
                            face_state.rho_avg);
                    const amrex::Real detj_donor = detj ? detj(i,j,donor_k) : amrex::Real(1.0);

                    const amrex::Real limited_rain = sam_limit_precip_component_flux(
                        corrected_fluxes.rain, donor.rho, donor.qpr, detj_donor, coef_substep);
                    const amrex::Real limited_snow = sam_limit_precip_component_flux(
                        corrected_fluxes.snow, donor.rho, donor.qps, detj_donor, coef_substep);
                    const amrex::Real limited_graupel = sam_limit_precip_component_flux(
                        corrected_fluxes.graupel, donor.rho, donor.qpg, detj_donor, coef_substep);

                    diagnostics.max_limited_fluxes[0] = std::max(diagnostics.max_limited_fluxes[0], limited_rain);
                    diagnostics.max_limited_fluxes[1] = std::max(diagnostics.max_limited_fluxes[1], limited_snow);
                    diagnostics.max_limited_fluxes[2] = std::max(diagnostics.max_limited_fluxes[2], limited_graupel);
                    diagnostics.limiter_active = diagnostics.limiter_active ||
                        (limited_rain + exact_zero_or_near_zero_tol() < corrected_fluxes.rain) ||
                        (limited_snow + exact_zero_or_near_zero_tol() < corrected_fluxes.snow) ||
                        (limited_graupel + exact_zero_or_near_zero_tol() < corrected_fluxes.graupel);
                }
            }
        }
    }

    return diagnostics;
}

void fill_precipfall_pure_rain_state (amrex::MultiFab& cons,
                                      const amrex::Real pres_mbar,
                                      const bool limiter_active)
{
    const int nz = cons.boxArray().minimalBox().length(2);
    cons.setVal(amrex::Real(0.0));

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto arr = cons.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const amrex::Real tabs = amrex::Real(284.0);
            const amrex::Real qsat = mixed_qsat_for_state(kSAMNoIceMode, tabs, pres_mbar);
            const amrex::Real qpr = (k == nz - 1)
                ? amrex::Real(0.0)
                : limiter_active
                    ? (amrex::Real(4.0e-6) + amrex::Real(1.0e-6) * (nz - 2 - k))
                    : (amrex::Real(6.0e-4) + amrex::Real(1.0e-4) * (nz - 2 - k));
            const SAMCellState state = make_cell_state(
                tabs,
                pres_mbar,
                amrex::Real(0.8) * qsat,
                amrex::Real(0.0),
                amrex::Real(0.0),
                qpr,
                amrex::Real(0.0),
                amrex::Real(0.0));

            arr(i,j,k,Rho_comp) = state.rho;
            arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
            arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
            sam_primitive_to_cons(state, arr, i, j, k);
        });
    }

    amrex::Gpu::streamSynchronize();
}

void fill_precipfall_mixed_component_state (amrex::MultiFab& cons,
                                            const amrex::Real pres_mbar)
{
    const int nz = cons.boxArray().minimalBox().length(2);
    cons.setVal(amrex::Real(0.0));

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto arr = cons.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const int phase_index = (i + j + k) % 3;
            const amrex::Real tabs = (phase_index == 0) ? amrex::Real(282.0)
                                    : (phase_index == 1) ? amrex::Real(262.0)
                                                         : amrex::Real(242.0);
            const amrex::Real qsat = mixed_qsat_for_state(kSAMWithIceMode, tabs, pres_mbar);
             const bool is_top = (k == nz - 1);
             const amrex::Real scale = is_top ? amrex::Real(0.0) : amrex::Real(nz - 1 - k);
            const SAMCellState state = make_cell_state(
                tabs,
                pres_mbar,
                amrex::Real(0.75) * qsat,
                amrex::Real(0.0),
                amrex::Real(0.0),
              is_top ? amrex::Real(0.0)
                  : amrex::Real(2.0e-5) + scale * (amrex::Real(3.0e-5) + amrex::Real(1.0e-6) * i),
              is_top ? amrex::Real(0.0)
                  : amrex::Real(1.5e-5) + scale * (amrex::Real(2.0e-5) + amrex::Real(1.0e-6) * j),
              is_top ? amrex::Real(0.0)
                  : amrex::Real(1.0e-5) + scale * (amrex::Real(2.5e-5) + amrex::Real(5.0e-7) * (i + j)));

            arr(i,j,k,Rho_comp) = state.rho;
            arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
            arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
            sam_primitive_to_cons(state, arr, i, j, k);
        });
    }

    amrex::Gpu::streamSynchronize();
}

void fill_single_cell_from_sam_state (amrex::MultiFab& cons,
                                      const SAMCellState& state)
{
    cons.setVal(amrex::Real(0.0));

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto arr = cons.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            arr(i,j,k,Rho_comp) = state.rho;
            arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
            arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
            sam_primitive_to_cons(state, arr, i, j, k);
        });
    }

    amrex::Gpu::streamSynchronize();
}

void fill_precipfall_detj (amrex::MultiFab& detJ_cc)
{
    for (amrex::MFIter mfi(detJ_cc, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto detj = detJ_cc.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            detj(i,j,k) = amrex::Real(0.8) + amrex::Real(0.05) * i + amrex::Real(0.03) * j + amrex::Real(0.07) * k;
        });
    }

    amrex::Gpu::streamSynchronize();
}

PrecipFallBudgetSummary run_precipfall_budget_case (const amrex::Geometry& geom,
                                                    amrex::MultiFab& cons,
                                                    const MoistureType moisture_type,
                                                    const amrex::Real dt,
                                                    amrex::MultiFab* detJ_cc = nullptr)
{
    PrecipFallBudgetSummary summary{};
    summary.initial_rain_mass = sum_component_mass(geom, cons, RhoQ4_comp, detJ_cc);
    summary.initial_snow_mass = sum_component_mass(geom, cons, RhoQ5_comp, detJ_cc);
    summary.initial_graupel_mass = sum_component_mass(geom, cons, RhoQ6_comp, detJ_cc);
    summary.diagnostics = compute_precipfall_run_diagnostics(geom, cons, dt, detJ_cc);

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detj_owned;
    if (detJ_cc != nullptr) {
        detj_owned = std::make_unique<amrex::MultiFab>(*detJ_cc, amrex::make_alias, 0, 1);
    }

    amrex::BoxArray ba = cons.boxArray();
    SolverChoice sc = make_sam_solver_choice(moisture_type);
    SAM sam;
    sam.Define(sc);
    sam.Set_dzmin(geom.CellSize(2));
    sam.Set_RealWidth(0);
    sam.Init(cons, ba, geom, dt, z_phys_nd, detj_owned);
    sam.Copy_State_to_Micro(cons);
    sam.PrecipFall(sc);
    sam.Copy_Micro_to_State(cons);
    amrex::Gpu::streamSynchronize();

    summary.final_rain_mass = sum_component_mass(geom, cons, RhoQ4_comp, detJ_cc);
    summary.final_snow_mass = sum_component_mass(geom, cons, RhoQ5_comp, detJ_cc);
    summary.final_graupel_mass = sum_component_mass(geom, cons, RhoQ6_comp, detJ_cc);
    summary.rain_accum_mass = sum_accum_mass(geom, *sam.Qmoist_Ptr(0), rhor);
    summary.snow_accum_mass = sum_accum_mass(geom, *sam.Qmoist_Ptr(1), rhos);
    summary.graupel_accum_mass = sum_accum_mass(geom, *sam.Qmoist_Ptr(2), rhog);

    summary.min_qpr_after = std::numeric_limits<amrex::Real>::max();
    summary.min_qps_after = std::numeric_limits<amrex::Real>::max();
    summary.min_qpg_after = std::numeric_limits<amrex::Real>::max();
    for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto arr = cons.const_array(mfi);
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    const SAMPrimitiveCell primitive = sam_cons_to_primitive(
                        arr(i,j,k,Rho_comp),
                        arr(i,j,k,RhoTheta_comp),
                        arr(i,j,k,RhoQ1_comp),
                        arr(i,j,k,RhoQ2_comp),
                        arr(i,j,k,RhoQ3_comp),
                        arr(i,j,k,RhoQ4_comp),
                        arr(i,j,k,RhoQ5_comp),
                        arr(i,j,k,RhoQ6_comp));
                    summary.min_qpr_after = std::min(summary.min_qpr_after, primitive.qpr);
                    summary.min_qps_after = std::min(summary.min_qps_after, primitive.qps);
                    summary.min_qpg_after = std::min(summary.min_qpg_after, primitive.qpg);
                }
            }
        }
    }

    return summary;
}

} // namespace

// Motivation:
// Local SAM precipitation sources should exercise all autoconversion,
// accretion, and evaporation pathways on a composed mixed-phase cell without
// creating negative species or violating total-water conservation.
TEST(SAMPhysicalProperties, PrecipSources_AllSpeciesNonzeroAllSourceTermsActive)
{
    const amrex::Real tabs = mixed_phase_tabs();
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qsat = mixed_qsat_for_state(kSAMWithIceMode, tabs, pres_mbar);
    const SAMCellState before = make_cell_state(
        tabs, pres_mbar,
        amrex::Real(0.5) * qsat,
        qcw0 + amrex::Real(8.0e-4),
        qci0 + amrex::Real(7.0e-4),
        amrex::Real(4.0e-4),
        amrex::Real(5.0e-4),
        amrex::Real(6.0e-4));

    SAMPrecipCellDiagnostics diagnostics;
    const SAMCellState after = sam_precip_cell_update(
        before, make_coeffs(), make_precip_config(kSAMWithIceMode), &diagnostics);

    EXPECT_GT(diagnostics.autoconversion.dqca, amrex::Real(0.0));
    EXPECT_GT(diagnostics.autoconversion.dqia, amrex::Real(0.0));
    EXPECT_GT(diagnostics.accretion.dprc, amrex::Real(0.0));
    EXPECT_GT(diagnostics.accretion.dpsc, amrex::Real(0.0));
    EXPECT_GT(diagnostics.accretion.dpgc, amrex::Real(0.0));
    EXPECT_GT(diagnostics.accretion.dpsi, amrex::Real(0.0));
    EXPECT_GT(diagnostics.accretion.dpgi, amrex::Real(0.0));
    EXPECT_GT(diagnostics.evaporation.dqpr, amrex::Real(0.0));
    EXPECT_GT(diagnostics.evaporation.dqps, amrex::Real(0.0));
    EXPECT_GT(diagnostics.evaporation.dqpg, amrex::Real(0.0));

    expect_nonnegative_species(after);
    EXPECT_NEAR(sam_total_water(after),
                sam_total_water(before),
                formula_abs_tol(sam_total_water(before)));
}

// Motivation:
// In a nonlimited mixed-phase cell where all source terms are active, the
// diagnostic autoconversion and accretion tendencies should match the closed-
// form formulas exactly, while the composed update still preserves the core
// physical invariants.
TEST(SAMPhysicalProperties, PrecipSources_AllSpeciesNonzeroExpectedDiagnostics)
{
    const amrex::Real tabs = mixed_phase_tabs();
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qsat = mixed_qsat_for_state(kSAMWithIceMode, tabs, pres_mbar);
    const SAMCellState before = make_cell_state(
        tabs, pres_mbar,
        amrex::Real(0.5) * qsat,
        qcw0 + amrex::Real(8.0e-4),
        qci0 + amrex::Real(7.0e-4),
        amrex::Real(4.0e-4),
        amrex::Real(5.0e-4),
        amrex::Real(6.0e-4));
    const SAMCoefficientRow coeffs = make_coeffs();
    const SAMPrecipConfig config = make_precip_config(kSAMWithIceMode);

    const amrex::Real expected_omp = std::max(amrex::Real(0.0),
                                              std::min(amrex::Real(1.0), (tabs - tprmin) * a_pr));
    const amrex::Real expected_omg = std::max(amrex::Real(0.0),
                                              std::min(amrex::Real(1.0), (tabs - tgrmin) * a_gr));
    const amrex::Real expected_dqca = config.dtn * alphaelq * (before.qcl - qcw0);
    const amrex::Real expected_dqia = config.dtn * betaelq * coeffs.coefice * (before.qci - qci0);
    const amrex::Real expected_dprc = config.dtn * coeffs.accrrc * before.qcl * std::pow(before.qpr, config.powr1);
    const amrex::Real expected_dpsc = config.dtn * coeffs.accrsc * before.qcl * std::pow(before.qps, config.pows1);
    const amrex::Real expected_dpgc = config.dtn * coeffs.accrgc * before.qcl * std::pow(before.qpg, config.powg1);
    const amrex::Real expected_dpsi = config.dtn * coeffs.accrsi * before.qci * std::pow(before.qps, config.pows1);
    const amrex::Real expected_dpgi = config.dtn * coeffs.accrgi * before.qci * std::pow(before.qpg, config.powg1);
    const amrex::Real expected_liquid_total = expected_dqca + expected_dprc + expected_dpsc + expected_dpgc;
    const amrex::Real expected_ice_total = expected_dqia + expected_dpsi + expected_dpgi;
    const amrex::Real expected_dqpr = (expected_dqca + expected_dqia) * expected_omp + expected_dprc;
    const amrex::Real expected_dqps = (expected_dqca + expected_dqia) * (one - expected_omp) * (one - expected_omg)
        + expected_dpsc + expected_dpsi;
    const amrex::Real expected_dqpg = (expected_dqca + expected_dqia) * (one - expected_omp) * expected_omg
        + expected_dpgc + expected_dpgi;

    ASSERT_LT(expected_liquid_total, before.qcl);
    ASSERT_LT(expected_ice_total, before.qci);

    SAMPrecipCellDiagnostics diagnostics;
    const SAMCellState after = sam_precip_cell_update(before, coeffs, config, &diagnostics);

    EXPECT_NEAR(diagnostics.omp, expected_omp, roundoff_tol(expected_omp));
    EXPECT_NEAR(diagnostics.omg, expected_omg, roundoff_tol(expected_omg));

    EXPECT_NEAR(diagnostics.autoconversion.dqca, expected_dqca, roundoff_tol(expected_dqca));
    EXPECT_NEAR(diagnostics.autoconversion.dqia, expected_dqia, roundoff_tol(expected_dqia));
    EXPECT_NEAR(diagnostics.accretion.dprc, expected_dprc, pow_sqrt_tol(expected_dprc));
    EXPECT_NEAR(diagnostics.accretion.dpsc, expected_dpsc, pow_sqrt_tol(expected_dpsc));
    EXPECT_NEAR(diagnostics.accretion.dpgc, expected_dpgc, pow_sqrt_tol(expected_dpgc));
    EXPECT_NEAR(diagnostics.accretion.dpsi, expected_dpsi, pow_sqrt_tol(expected_dpsi));
    EXPECT_NEAR(diagnostics.accretion.dpgi, expected_dpgi, pow_sqrt_tol(expected_dpgi));

    EXPECT_NEAR(diagnostics.limited_sources.dqca, expected_dqca, roundoff_tol(expected_dqca));
    EXPECT_NEAR(diagnostics.limited_sources.dqia, expected_dqia, roundoff_tol(expected_dqia));
    EXPECT_NEAR(diagnostics.limited_sources.dprc, expected_dprc, pow_sqrt_tol(expected_dprc));
    EXPECT_NEAR(diagnostics.limited_sources.dpsc, expected_dpsc, pow_sqrt_tol(expected_dpsc));
    EXPECT_NEAR(diagnostics.limited_sources.dpgc, expected_dpgc, pow_sqrt_tol(expected_dpgc));
    EXPECT_NEAR(diagnostics.limited_sources.dpsi, expected_dpsi, pow_sqrt_tol(expected_dpsi));
    EXPECT_NEAR(diagnostics.limited_sources.dpgi, expected_dpgi, pow_sqrt_tol(expected_dpgi));
    EXPECT_NEAR(diagnostics.limited_sources.dqc, expected_liquid_total, pow_sqrt_tol(expected_liquid_total));
    EXPECT_NEAR(diagnostics.limited_sources.dqi, expected_ice_total, pow_sqrt_tol(expected_ice_total));

    EXPECT_NEAR(diagnostics.partitioned_sources.dqpr, expected_dqpr, pow_sqrt_tol(expected_dqpr));
    EXPECT_NEAR(diagnostics.partitioned_sources.dqps, expected_dqps, pow_sqrt_tol(expected_dqps));
    EXPECT_NEAR(diagnostics.partitioned_sources.dqpg, expected_dqpg, pow_sqrt_tol(expected_dqpg));
    EXPECT_NEAR(diagnostics.partitioned_sources.dqpr
                    + diagnostics.partitioned_sources.dqps
                    + diagnostics.partitioned_sources.dqpg,
                diagnostics.limited_sources.dqc + diagnostics.limited_sources.dqi,
                property_accumulation_tol(5, diagnostics.limited_sources.dqc + diagnostics.limited_sources.dqi));

    expect_nonnegative_species(after);
    EXPECT_NEAR(sam_total_water(after),
                sam_total_water(before),
                formula_abs_tol(sam_total_water(before)));
    EXPECT_NEAR(sam_latent_proxy(after, kFacCond, kFacFus),
                sam_latent_proxy(before, kFacCond, kFacFus),
                formula_abs_tol(sam_latent_proxy(before, kFacCond, kFacFus)));
}

// Motivation:
// Local SAM precipitation sources exchange water among qv, qcl, qci, qpr,
// qps, and qpg only. They should conserve total water across representative
// warm, cold, mixed, limiter-active, no-precip, and no-ice cases.
TEST(SAMPhysicalProperties, PrecipSources_ConserveTotalWaterAcrossRepresentativeCases)
{
    const amrex::Real mixed_tabs = mixed_phase_tabs();
    const amrex::Real warm_tabs = tbgmax + amrex::Real(2.0);
    const amrex::Real cold_tabs = tgrmin + amrex::Real(0.25) * (tgrmax - tgrmin);
    const amrex::Real pres_mbar = amrex::Real(900.0);

    struct ConservationCase {
        const char* label;
        SAMCellState state;
        SAMCoefficientRow coeffs;
        SAMPrecipConfig config;
    };

    const std::array<ConservationCase, 7> cases = {{
        {"warm-rain-only",
         make_cell_state(warm_tabs, pres_mbar,
                         amrex::Real(0.8) * mixed_qsat_for_state(kSAMNoIceMode, warm_tabs, pres_mbar),
                         qcw0 + amrex::Real(7.0e-4),
                         amrex::Real(0.0),
                         amrex::Real(4.0e-4),
                         amrex::Real(0.0),
                         amrex::Real(0.0)),
         make_coeffs(),
         make_precip_config(kSAMNoIceMode)},
        {"cold-snow-graupel-only",
         make_cell_state(cold_tabs, pres_mbar,
                         amrex::Real(0.7) * mixed_qsat_for_state(kSAMWithIceMode, cold_tabs, pres_mbar),
                         amrex::Real(0.0),
                         qci0 + amrex::Real(6.0e-4),
                         amrex::Real(0.0),
                         amrex::Real(4.0e-4),
                         amrex::Real(5.0e-4)),
         make_coeffs(),
         make_precip_config(kSAMWithIceMode)},
        {"mixed-all-species",
         make_cell_state(mixed_tabs, pres_mbar,
                         amrex::Real(0.6) * mixed_qsat_for_state(kSAMWithIceMode, mixed_tabs, pres_mbar),
                         qcw0 + amrex::Real(8.0e-4),
                         qci0 + amrex::Real(7.0e-4),
                         amrex::Real(4.0e-4),
                         amrex::Real(5.0e-4),
                         amrex::Real(6.0e-4)),
         make_coeffs(),
         make_precip_config(kSAMWithIceMode)},
        {"cloud-sink-limiter",
         make_cell_state(mixed_tabs, pres_mbar,
                         amrex::Real(0.7) * mixed_qsat_for_state(kSAMWithIceMode, mixed_tabs, pres_mbar),
                         amrex::Real(3.0e-5),
                         amrex::Real(2.0e-5),
                         amrex::Real(4.0e-4),
                         amrex::Real(5.0e-4),
                         amrex::Real(6.0e-4)),
         make_coeffs(amrex::Real(1.0e3)),
         make_precip_config(kSAMWithIceMode, true, amrex::Real(10.0))},
        {"evaporation-limiter",
         make_cell_state(mixed_tabs, pres_mbar,
                         amrex::Real(1.0e-6),
                         qcw0 + amrex::Real(2.0e-4),
                         qci0 + amrex::Real(2.0e-4),
                         amrex::Real(2.0e-6),
                         amrex::Real(3.0e-6),
                         amrex::Real(4.0e-6)),
         make_coeffs(amrex::Real(5.0e3)),
         make_precip_config(kSAMWithIceMode, true, amrex::Real(5.0))},
        {"no-precip-mode",
         make_cell_state(mixed_tabs, pres_mbar,
                         amrex::Real(0.5) * mixed_qsat_for_state(kSAMWithIceMode, mixed_tabs, pres_mbar),
                         qcw0 + amrex::Real(4.0e-4),
                         qci0 + amrex::Real(3.0e-4),
                         amrex::Real(2.0e-4),
                         amrex::Real(2.0e-4),
                         amrex::Real(2.0e-4)),
         make_coeffs(),
         make_precip_config(kSAMWithIceMode, false)},
        {"no-ice-mode",
         make_cell_state(warm_tabs, pres_mbar,
                         amrex::Real(0.75) * mixed_qsat_for_state(kSAMNoIceMode, warm_tabs, pres_mbar),
                         qcw0 + amrex::Real(5.0e-4),
                         amrex::Real(0.0),
                         amrex::Real(3.0e-4),
                         amrex::Real(0.0),
                         amrex::Real(0.0)),
         make_coeffs(),
         make_precip_config(kSAMNoIceMode)}
    }};

    for (const ConservationCase& test_case : cases) {
        SCOPED_TRACE(test_case.label);
        const amrex::Real before_total = sam_total_water(test_case.state);
        const SAMCellState after = sam_precip_cell_update(test_case.state,
                                                          test_case.coeffs,
                                                          test_case.config);
        expect_nonnegative_species(after);
        EXPECT_NEAR(sam_total_water(after),
                    before_total,
                    formula_abs_tol(before_total));
    }
}

// Motivation:
// The composed SAM precipitation source update should conserve the latent proxy
// T + fac_cond*qv - fac_fus*(qci + qps + qpg). A mixed-phase case with cloud-
// liquid accretion onto snow and graupel directly exposes the candidate missing
// fusion-heating term if that contract is violated.
TEST(SAMPhysicalProperties, PrecipSources_ConserveLatentProxyOrExposeRimingHeatingBug)
{
    const amrex::Real tabs = mixed_phase_tabs();
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qsat = mixed_qsat_for_state(kSAMWithIceMode, tabs, pres_mbar);
    const SAMCellState before = make_cell_state(
        tabs, pres_mbar,
        qsat,
        qcw0 + amrex::Real(8.0e-4),
        qci0 + amrex::Real(2.0e-4),
        amrex::Real(2.0e-4),
        amrex::Real(6.0e-4),
        amrex::Real(7.0e-4));

    SAMPrecipCellDiagnostics diagnostics;
    const SAMCellState after = sam_precip_cell_update(
        before, make_coeffs(), make_precip_config(kSAMWithIceMode), &diagnostics);

    ASSERT_GT(diagnostics.limited_sources.dpsc, amrex::Real(0.0));
    ASSERT_GT(diagnostics.limited_sources.dpgc, amrex::Real(0.0));

    const amrex::Real before_proxy = sam_latent_proxy(before, kFacCond, kFacFus);
    const amrex::Real after_proxy = sam_latent_proxy(after, kFacCond, kFacFus);
    const amrex::Real expected_missing =
        -kFacFus * (diagnostics.limited_sources.dpsc + diagnostics.limited_sources.dpgc);

    EXPECT_NEAR(after_proxy,
                before_proxy,
                formula_abs_tol(before_proxy))
        << "expected_missing_if_unheated=" << expected_missing
        << " dpsc=" << diagnostics.limited_sources.dpsc
        << " dpgc=" << diagnostics.limited_sources.dpgc;
}

// Motivation:
// The public SAM Precip path should conserve the latent proxy under the
// intended constant-pressure microphysics contract when evaluated with the
// initial pressure snapshot. The EOS-projected proxy based on conserved
// density and rho*theta is diagnostic only here.
TEST(SAMPhysicalProperties, PublicPrecipConstantPressureLatentProxySingleCell)
{
    const amrex::Real tabs = mixed_phase_tabs();
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qsat = mixed_qsat_for_state(kSAMWithIceMode, tabs, pres_mbar);
    const SAMCellState state = make_cell_state(
        tabs,
        pres_mbar,
        amrex::Real(0.6) * qsat,
        qcw0 + amrex::Real(8.0e-4),
        qci0 + amrex::Real(7.0e-4),
        amrex::Real(4.0e-4),
        amrex::Real(5.0e-4),
        amrex::Real(6.0e-4));

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(1.0, 1.0, 1.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ6_comp + 1, 0);
    fill_single_cell_from_sam_state(cons, state);

    amrex::MultiFab pressure0_mbar(ba, dm, 1, 0);
    fill_pressure_snapshot_from_cons(cons, pressure0_mbar);

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;

    SolverChoice sc = make_sam_solver_choice(MoistureType::SAM);
    SAM sam;
    sam.Define(sc);
    sam.Set_dzmin(geom.CellSize(2));
    sam.Set_RealWidth(0);
    sam.Init(cons, ba, geom, amrex::Real(0.1), z_phys_nd, detJ_cc);

    const SAMPressureSnapshotBudget initial_budget =
        compute_sam_pressure_snapshot_budget(geom, cons, pressure0_mbar);

    sam.Copy_State_to_Micro(cons);
    sam.Compute_Coefficients();
    sam.Precip(sc);
    sam.Copy_Micro_to_State(cons);

    const SAMPressureSnapshotBudget final_budget =
        compute_sam_pressure_snapshot_budget(geom, cons, pressure0_mbar);

    const amrex::Real constant_pressure_residual =
        final_budget.latent_proxy_mass_constant_pressure -
        initial_budget.latent_proxy_mass_constant_pressure;
    const amrex::Real eos_projected_residual =
        final_budget.latent_proxy_mass_eos_projected -
        initial_budget.latent_proxy_mass_eos_projected;
    const amrex::Real total_water_residual =
        final_budget.total_water_mass - initial_budget.total_water_mass;
    const amrex::Real scale = std::max({std::abs(initial_budget.latent_proxy_mass_constant_pressure),
                                        std::abs(final_budget.latent_proxy_mass_constant_pressure),
                                        amrex::Real(1.0)});

    EXPECT_NEAR(final_budget.total_water_mass,
                initial_budget.total_water_mass,
                property_accumulation_tol(6, std::max(std::abs(initial_budget.total_water_mass), amrex::Real(1.0))))
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
        << " dt=" << amrex::Real(0.1)
        << " rank=" << amrex::ParallelDescriptor::MyProc()
        << " box_decomposition=(1,1,1)";
    EXPECT_NEAR(final_budget.latent_proxy_mass_constant_pressure,
                initial_budget.latent_proxy_mass_constant_pressure,
                property_accumulation_tol(8, scale))
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
        << " dt=" << amrex::Real(0.1)
        << " rank=" << amrex::ParallelDescriptor::MyProc()
        << " box_decomposition=(1,1,1)";
    EXPECT_NEAR(sum_multifab_scalar(*sam.Qmoist_Ptr(0)), amrex::Real(0.0), exact_zero_tol());
    EXPECT_NEAR(sum_multifab_scalar(*sam.Qmoist_Ptr(1)), amrex::Real(0.0), exact_zero_tol());
    EXPECT_NEAR(sum_multifab_scalar(*sam.Qmoist_Ptr(2)), amrex::Real(0.0), exact_zero_tol());
}

// Motivation:
// With pure rain and inactive donor limiting, PrecipFall should conserve the
// rain column budget exactly: initial rain mass equals final rain mass plus
// bottom rain accumulation.
TEST(SAMPhysicalProperties, PrecipFall_PureRainNoClipColumnBudgetCloses)
{
    constexpr int nx = 2;
    constexpr int ny = 1;
    constexpr int nz = 4;
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real dt = amrex::Real(0.01);

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(2.0, 1.0, 4.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 0)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    ba.maxSize(amrex::IntVect(1, 1, nz));
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(amrex::The_Pinned_Arena());
    cons.define(ba, dm, RhoQ6_comp + 1, 0);
    fill_precipfall_pure_rain_state(cons, pres_mbar, false);

    const PrecipFallBudgetSummary summary =
        run_precipfall_budget_case(geom, cons, MoistureType::SAM_NoIce, dt);
    const amrex::Real residual =
        summary.initial_rain_mass - summary.final_rain_mass - summary.rain_accum_mass;

    EXPECT_NEAR(residual,
                amrex::Real(0.0),
                property_accumulation_tol(nx * nz,
                                          std::max(summary.initial_rain_mass, summary.rain_accum_mass)))
        << "initial_rain_mass=" << summary.initial_rain_mass
        << " final_rain_mass=" << summary.final_rain_mass
        << " rain_accum_mass=" << summary.rain_accum_mass
        << " residual=" << residual
        << " dt=" << dt
        << " dz=" << geom.CellSize(2)
        << " n_substep=" << summary.diagnostics.n_substep
        << " max_raw_rain_flux=" << summary.diagnostics.max_raw_fluxes[0]
        << " max_limited_rain_flux=" << summary.diagnostics.max_limited_fluxes[0]
        << " limiter_active=" << summary.diagnostics.limiter_active
        << " min_qpr_after=" << summary.min_qpr_after
        << " detj_min=" << summary.diagnostics.detj_min
        << " detj_max=" << summary.diagnostics.detj_max;
}

// Motivation:
// With an active donor limiter, PrecipFall should still close the pure-rain
// column budget and keep rain nonnegative.
TEST(SAMPhysicalProperties, PrecipFall_PureRainLimiterActiveColumnBudgetCloses)
{
    constexpr int nx = 2;
    constexpr int ny = 1;
    constexpr int nz = 4;
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real dt = amrex::Real(1.0);

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(2.0, 1.0, 4.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 0)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    ba.maxSize(amrex::IntVect(1, 1, nz));
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(amrex::The_Pinned_Arena());
    cons.define(ba, dm, RhoQ6_comp + 1, 0);
    fill_precipfall_pure_rain_state(cons, pres_mbar, true);

    const PrecipFallBudgetSummary summary =
        run_precipfall_budget_case(geom, cons, MoistureType::SAM_NoIce, dt);
    const amrex::Real residual =
        summary.initial_rain_mass - summary.final_rain_mass - summary.rain_accum_mass;

    EXPECT_TRUE(summary.diagnostics.limiter_active);
    EXPECT_GE(summary.min_qpr_after, -exact_zero_or_near_zero_tol());
    EXPECT_NEAR(residual,
                amrex::Real(0.0),
                property_accumulation_tol(nx * nz,
                                          std::max(summary.initial_rain_mass, summary.rain_accum_mass)))
        << "initial_rain_mass=" << summary.initial_rain_mass
        << " final_rain_mass=" << summary.final_rain_mass
        << " rain_accum_mass=" << summary.rain_accum_mass
        << " residual=" << residual
        << " dt=" << dt
        << " dz=" << geom.CellSize(2)
        << " n_substep=" << summary.diagnostics.n_substep
        << " max_raw_rain_flux=" << summary.diagnostics.max_raw_fluxes[0]
        << " max_limited_rain_flux=" << summary.diagnostics.max_limited_fluxes[0]
        << " limiter_active=" << summary.diagnostics.limiter_active
        << " min_qpr_after=" << summary.min_qpr_after
        << " detj_min=" << summary.diagnostics.detj_min
        << " detj_max=" << summary.diagnostics.detj_max;
}

// Motivation:
// Mixed rain, snow, and graupel columns should conserve each precip component
// separately under PrecipFall's component-wise sedimentation update.
TEST(SAMPhysicalProperties, PrecipFall_RainSnowGraupelComponentBudgetsClose)
{
    constexpr int nx = 2;
    constexpr int ny = 2;
    constexpr int nz = 4;
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real dt = amrex::Real(0.25);

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(2.0, 2.0, 4.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 0)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    ba.maxSize(amrex::IntVect(1, 1, nz));
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(amrex::The_Pinned_Arena());
    cons.define(ba, dm, RhoQ6_comp + 1, 0);
    fill_precipfall_mixed_component_state(cons, pres_mbar);

    const PrecipFallBudgetSummary summary =
        run_precipfall_budget_case(geom, cons, MoistureType::SAM, dt);
    const amrex::Real rain_residual =
        summary.initial_rain_mass - summary.final_rain_mass - summary.rain_accum_mass;
    const amrex::Real snow_residual =
        summary.initial_snow_mass - summary.final_snow_mass - summary.snow_accum_mass;
    const amrex::Real graupel_residual =
        summary.initial_graupel_mass - summary.final_graupel_mass - summary.graupel_accum_mass;
    const amrex::Real total_residual = rain_residual + snow_residual + graupel_residual;

    EXPECT_GE(summary.min_qpr_after, -exact_zero_or_near_zero_tol());
    EXPECT_GE(summary.min_qps_after, -exact_zero_or_near_zero_tol());
    EXPECT_GE(summary.min_qpg_after, -exact_zero_or_near_zero_tol());

    EXPECT_NEAR(rain_residual,
                amrex::Real(0.0),
                property_accumulation_tol(nx * ny * nz,
                                          std::max(summary.initial_rain_mass, summary.rain_accum_mass)))
        << "rain initial=" << summary.initial_rain_mass
        << " rain final=" << summary.final_rain_mass
        << " rain accum=" << summary.rain_accum_mass
        << " rain residual=" << rain_residual
        << " dt=" << dt
        << " dz=" << geom.CellSize(2)
        << " n_substep=" << summary.diagnostics.n_substep
        << " max_raw_fluxes=(" << summary.diagnostics.max_raw_fluxes[0] << ","
        << summary.diagnostics.max_raw_fluxes[1] << "," << summary.diagnostics.max_raw_fluxes[2] << ")"
        << " max_limited_fluxes=(" << summary.diagnostics.max_limited_fluxes[0] << ","
        << summary.diagnostics.max_limited_fluxes[1] << "," << summary.diagnostics.max_limited_fluxes[2] << ")"
        << " limiter_active=" << summary.diagnostics.limiter_active
        << " min_q_after=(" << summary.min_qpr_after << "," << summary.min_qps_after << ","
        << summary.min_qpg_after << ")"
        << " detj_min=" << summary.diagnostics.detj_min
        << " detj_max=" << summary.diagnostics.detj_max;
    EXPECT_NEAR(snow_residual,
                amrex::Real(0.0),
                property_accumulation_tol(nx * ny * nz,
                                          std::max(summary.initial_snow_mass, summary.snow_accum_mass)));
    EXPECT_NEAR(graupel_residual,
                amrex::Real(0.0),
                property_accumulation_tol(nx * ny * nz,
                                          std::max(summary.initial_graupel_mass, summary.graupel_accum_mass)));
    EXPECT_NEAR(total_residual,
                amrex::Real(0.0),
                property_accumulation_tol(3 * nx * ny * nz,
                                          std::max({summary.initial_rain_mass + summary.initial_snow_mass + summary.initial_graupel_mass,
                                                    summary.rain_accum_mass + summary.snow_accum_mass + summary.graupel_accum_mass,
                                                    amrex::Real(1.0)})));
}

    struct CloudAdjustmentLimiterInvestigation {
        bool found{false};
        SAMCellState initial_state{};
        SAMCellState partitioned_state{};
        SAMCellState final_state{};
        amrex::Real qsat_final{amrex::Real(0.0)};
        amrex::Real requested_delta_qc{amrex::Real(0.0)};
        amrex::Real requested_delta_qi{amrex::Real(0.0)};
        bool clipped_qc{false};
        bool clipped_qi{false};
    };

    amrex::Real cloud_adjustment_total_nonprecip_water (const SAMCellState& state)
    {
        return state.qv + state.qcl + state.qci;
    }

    amrex::Real cloud_adjustment_latent_proxy (const SAMCellState& state)
    {
        return state.tabs + kFacCond * state.qv - kFacFus * state.qci;
    }

    CloudAdjustmentLimiterInvestigation find_cloud_adjustment_limiter_active_case ()
    {
        constexpr amrex::Real an = a_bg;
        constexpr amrex::Real bn = tbgmin * a_bg;

        const amrex::Real pres_mbar = amrex::Real(900.0);
        const std::array<amrex::Real, 9> tabs_samples = {
            tbgmin + amrex::Real(0.1),
            tbgmin + amrex::Real(1.0),
            tbgmin + amrex::Real(2.5),
            tbgmin + amrex::Real(5.0),
            amrex::Real(0.5) * (tbgmin + tbgmax),
            tbgmax - amrex::Real(5.0),
            tbgmax - amrex::Real(2.5),
            tbgmax - amrex::Real(1.0),
            tbgmax - amrex::Real(0.1)};
        const std::array<amrex::Real, 7> qn_samples = {
            amrex::Real(2.0e-4),
            amrex::Real(4.0e-4),
            amrex::Real(8.0e-4),
            amrex::Real(1.2e-3),
            amrex::Real(1.8e-3),
            amrex::Real(2.4e-3),
            amrex::Real(3.0e-3)};
        const std::array<amrex::Real, 6> liquid_fractions = {
            amrex::Real(0.0),
            amrex::Real(0.15),
            amrex::Real(0.35),
            amrex::Real(0.65),
            amrex::Real(0.85),
            amrex::Real(1.0)};
        const std::array<amrex::Real, 8> deficit_fractions = {
            amrex::Real(0.05),
            amrex::Real(0.15),
            amrex::Real(0.25),
            amrex::Real(0.40),
            amrex::Real(0.55),
            amrex::Real(0.70),
            amrex::Real(0.85),
            amrex::Real(0.95)};

        for (const amrex::Real tabs_initial : tabs_samples) {
            for (const amrex::Real qn_initial : qn_samples) {
                for (const amrex::Real liquid_fraction : liquid_fractions) {
                    const amrex::Real qcl_initial = liquid_fraction * qn_initial;
                    const amrex::Real qci_initial = (one - liquid_fraction) * qn_initial;
                    const SAMCloudPhaseChange partition = sam_partition_cloud_phase(
                        kSAMWithIceMode,
                        tabs_initial,
                        qn_initial,
                        qcl_initial,
                        qci_initial,
                        kFacCond,
                        kFacFus,
                        an,
                        bn);
                    const amrex::Real qsat_partition =
                        mixed_qsat_for_state(kSAMWithIceMode, partition.tabs, pres_mbar);

                    for (const amrex::Real deficit_fraction : deficit_fractions) {
                        const amrex::Real qv_initial = qsat_partition - deficit_fraction * qn_initial;
                        if (qv_initial <= amrex::Real(0.0)) {
                            continue;
                        }

                        const SAMCellState initial_state = make_cell_state(
                            tabs_initial,
                            pres_mbar,
                            qv_initial,
                            qcl_initial,
                            qci_initial,
                            amrex::Real(0.0),
                            amrex::Real(0.0),
                            amrex::Real(0.0));
                        const SAMCellState partitioned_state = make_cell_state(
                            partition.tabs,
                            pres_mbar,
                            qv_initial,
                            partition.qcl,
                            partition.qci,
                            amrex::Real(0.0),
                            amrex::Real(0.0),
                            amrex::Real(0.0));

                        if (partitioned_state.qt <= qsat_partition) {
                            continue;
                        }

                        amrex::FArrayBox tabs_fab(single_cell_box(), 1, amrex::The_Pinned_Arena());
                        amrex::FArrayBox pres_fab(single_cell_box(), 1, amrex::The_Pinned_Arena());
                        amrex::FArrayBox qv_fab(single_cell_box(), 1, amrex::The_Pinned_Arena());
                        amrex::FArrayBox qc_fab(single_cell_box(), 1, amrex::The_Pinned_Arena());
                        amrex::FArrayBox qi_fab(single_cell_box(), 1, amrex::The_Pinned_Arena());
                        amrex::FArrayBox qn_fab(single_cell_box(), 1, amrex::The_Pinned_Arena());
                        amrex::FArrayBox qt_fab(single_cell_box(), 1, amrex::The_Pinned_Arena());

                        auto tabs_array = tabs_fab.array();
                        auto pres_array = pres_fab.array();
                        auto qv_array = qv_fab.array();
                        auto qc_array = qc_fab.array();
                        auto qi_array = qi_fab.array();
                        auto qn_array = qn_fab.array();
                        auto qt_array = qt_fab.array();

                        tabs_array(0,0,0) = partitioned_state.tabs;
                        pres_array(0,0,0) = partitioned_state.pres_mbar;
                        qv_array(0,0,0) = partitioned_state.qv;
                        qc_array(0,0,0) = partitioned_state.qcl;
                        qi_array(0,0,0) = partitioned_state.qci;
                        qn_array(0,0,0) = partitioned_state.qn;
                        qt_array(0,0,0) = partitioned_state.qt;

                        int i = 0;
                        int j = 0;
                        int k = 0;
                        const amrex::Real tabs_final = SAM::NewtonIterSat(
                            i,
                            j,
                            k,
                            kSAMWithIceMode,
                            kFacCond,
                            kFacFus,
                            kFacSub,
                            an,
                            bn,
                            tabs_array,
                            pres_array,
                            qv_array,
                            qc_array,
                            qi_array,
                            qn_array,
                            qt_array);

                        CloudAdjustmentLimiterInvestigation result{};
                        result.initial_state = initial_state;
                        result.partitioned_state = partitioned_state;
                        result.final_state = make_cell_state(
                            tabs_final,
                            pres_mbar,
                            qv_array(0,0,0),
                            qc_array(0,0,0),
                            qi_array(0,0,0),
                            amrex::Real(0.0),
                            amrex::Real(0.0),
                            amrex::Real(0.0));

                        amrex::Real qsatw_final;
                        amrex::Real qsati_final;
                        erf_qsatw(tabs_final, pres_mbar, qsatw_final);
                        erf_qsati(tabs_final, pres_mbar, qsati_final);
                        const amrex::Real omn_final =
                            sam_cloud_liquid_fraction(kSAMWithIceMode, tabs_final, an, bn);
                        result.qsat_final = sam_mixed_qsat(omn_final, qsatw_final, qsati_final);
                        result.requested_delta_qc =
                            (partitioned_state.qv - result.qsat_final) * omn_final;
                        result.requested_delta_qi =
                            (partitioned_state.qv - result.qsat_final) * (one - omn_final);

                        const amrex::Real actual_delta_qc =
                            result.final_state.qcl - partitioned_state.qcl;
                        const amrex::Real actual_delta_qi =
                            result.final_state.qci - partitioned_state.qci;
                        const amrex::Real qc_tol = property_accumulation_tol(
                            2,
                            std::max({std::abs(result.requested_delta_qc),
                                      std::abs(actual_delta_qc),
                                      amrex::Real(1.0)}));
                        const amrex::Real qi_tol = property_accumulation_tol(
                            2,
                            std::max({std::abs(result.requested_delta_qi),
                                      std::abs(actual_delta_qi),
                                      amrex::Real(1.0)}));
                        result.clipped_qc =
                            std::abs(actual_delta_qc - result.requested_delta_qc) > qc_tol;
                        result.clipped_qi =
                            std::abs(actual_delta_qi - result.requested_delta_qi) > qi_tol;

                        if (result.clipped_qc || result.clipped_qi) {
                            result.found = true;
                            return result;
                        }
                    }
                }
            }
        }

        return {};
    }

    // Motivation:
    // Reproduce the same held-pressure cloud-adjustment sequence used by SAM::Cloud
    // on a valid one-cell state: phase repartition first, then NewtonIterSat. The
    // full non-precipitating cloud-adjustment contract is total water qv+qcl+qci
    // and the latent proxy T + fac_cond*qv - fac_fus*qci. This test searches for a
    // valid limiter-active NewtonIterSat case and verifies that the combined cloud
    // adjustment still preserves that contract.
    TEST(SAMPhysicalProperties, CloudAdjustmentLimiterActiveConservesWaterAndLatentProxy)
    {
        const CloudAdjustmentLimiterInvestigation investigation =
            find_cloud_adjustment_limiter_active_case();

        ASSERT_TRUE(investigation.found)
            << "No valid Cloud/NewtonIterSat limiter-active case was found in the deterministic search grid.";
        ASSERT_TRUE(investigation.clipped_qc || investigation.clipped_qi);

        const amrex::Real initial_total_water =
            cloud_adjustment_total_nonprecip_water(investigation.initial_state);
        const amrex::Real final_total_water =
            cloud_adjustment_total_nonprecip_water(investigation.final_state);
        const amrex::Real initial_latent_proxy =
            cloud_adjustment_latent_proxy(investigation.initial_state);
        const amrex::Real final_latent_proxy =
            cloud_adjustment_latent_proxy(investigation.final_state);
        const amrex::Real total_water_residual = final_total_water - initial_total_water;
        const amrex::Real latent_proxy_residual = final_latent_proxy - initial_latent_proxy;
        const amrex::Real qv_minus_qsat_final =
            investigation.final_state.qv - investigation.qsat_final;
        const amrex::Real qsat_tol =
            property_accumulation_tol(6, std::max(std::abs(investigation.qsat_final), amrex::Real(1.0)));
        // NewtonIterSat stops on a temperature update tolerance, so the final
        // cloud-adjustment latent proxy should close far tighter than that, but
        // not necessarily to pure accumulation roundoff.
        const amrex::Real latent_proxy_tol = std::max(
            property_accumulation_tol(6, std::max(std::abs(initial_latent_proxy), amrex::Real(1.0))),
            amrex::Real(1.0e-9));

        EXPECT_NEAR(final_total_water,
                    initial_total_water,
                    property_accumulation_tol(6, std::max(std::abs(initial_total_water), amrex::Real(1.0))))
            << "initial_qv=" << investigation.initial_state.qv
            << " initial_qcl=" << investigation.initial_state.qcl
            << " initial_qci=" << investigation.initial_state.qci
            << " initial_qn=" << investigation.initial_state.qn
            << " initial_qt=" << investigation.initial_state.qt
            << " final_qv=" << investigation.final_state.qv
            << " final_qcl=" << investigation.final_state.qcl
            << " final_qci=" << investigation.final_state.qci
            << " final_qn=" << investigation.final_state.qn
            << " final_qt=" << investigation.final_state.qt
            << " partitioned_qv=" << investigation.partitioned_state.qv
            << " partitioned_qcl=" << investigation.partitioned_state.qcl
            << " partitioned_qci=" << investigation.partitioned_state.qci
            << " initial_tabs=" << investigation.initial_state.tabs
            << " partitioned_tabs=" << investigation.partitioned_state.tabs
            << " final_tabs=" << investigation.final_state.tabs
            << " pressure_mbar=" << investigation.initial_state.pres_mbar
            << " initial_total_water=" << initial_total_water
            << " final_total_water=" << final_total_water
            << " total_water_residual=" << total_water_residual
            << " initial_latent_proxy=" << initial_latent_proxy
            << " final_latent_proxy=" << final_latent_proxy
            << " latent_proxy_residual=" << latent_proxy_residual
            << " qv_minus_qsat_final=" << qv_minus_qsat_final
            << " clipped_qc=" << investigation.clipped_qc
            << " clipped_qi=" << investigation.clipped_qi
            << " requested_delta_qc=" << investigation.requested_delta_qc
            << " requested_delta_qi=" << investigation.requested_delta_qi;
            EXPECT_NEAR(qv_minus_qsat_final,
                    amrex::Real(0.0),
                    qsat_tol)
                << "initial_qv=" << investigation.initial_state.qv
                << " initial_qcl=" << investigation.initial_state.qcl
                << " initial_qci=" << investigation.initial_state.qci
                << " initial_qn=" << investigation.initial_state.qn
                << " initial_qt=" << investigation.initial_state.qt
                << " final_qv=" << investigation.final_state.qv
                << " final_qcl=" << investigation.final_state.qcl
                << " final_qci=" << investigation.final_state.qci
                << " final_qn=" << investigation.final_state.qn
                << " final_qt=" << investigation.final_state.qt
                << " partitioned_qv=" << investigation.partitioned_state.qv
                << " partitioned_qcl=" << investigation.partitioned_state.qcl
                << " partitioned_qci=" << investigation.partitioned_state.qci
                << " initial_tabs=" << investigation.initial_state.tabs
                << " partitioned_tabs=" << investigation.partitioned_state.tabs
                << " final_tabs=" << investigation.final_state.tabs
                << " pressure_mbar=" << investigation.initial_state.pres_mbar
                << " initial_total_water=" << initial_total_water
                << " final_total_water=" << final_total_water
                << " total_water_residual=" << total_water_residual
                << " initial_latent_proxy=" << initial_latent_proxy
                << " final_latent_proxy=" << final_latent_proxy
                << " latent_proxy_residual=" << latent_proxy_residual
                << " qv_minus_qsat_final=" << qv_minus_qsat_final
                << " clipped_qc=" << investigation.clipped_qc
                << " clipped_qi=" << investigation.clipped_qi
                << " requested_delta_qc=" << investigation.requested_delta_qc
                << " requested_delta_qi=" << investigation.requested_delta_qi;
        EXPECT_NEAR(final_latent_proxy,
                    initial_latent_proxy,
                    latent_proxy_tol)
            << "initial_qv=" << investigation.initial_state.qv
            << " initial_qcl=" << investigation.initial_state.qcl
            << " initial_qci=" << investigation.initial_state.qci
            << " initial_qn=" << investigation.initial_state.qn
            << " initial_qt=" << investigation.initial_state.qt
            << " final_qv=" << investigation.final_state.qv
            << " final_qcl=" << investigation.final_state.qcl
            << " final_qci=" << investigation.final_state.qci
            << " final_qn=" << investigation.final_state.qn
            << " final_qt=" << investigation.final_state.qt
            << " partitioned_qv=" << investigation.partitioned_state.qv
            << " partitioned_qcl=" << investigation.partitioned_state.qcl
            << " partitioned_qci=" << investigation.partitioned_state.qci
            << " initial_tabs=" << investigation.initial_state.tabs
            << " partitioned_tabs=" << investigation.partitioned_state.tabs
            << " final_tabs=" << investigation.final_state.tabs
            << " pressure_mbar=" << investigation.initial_state.pres_mbar
            << " initial_total_water=" << initial_total_water
            << " final_total_water=" << final_total_water
            << " total_water_residual=" << total_water_residual
            << " initial_latent_proxy=" << initial_latent_proxy
            << " final_latent_proxy=" << final_latent_proxy
            << " latent_proxy_residual=" << latent_proxy_residual
            << " qv_minus_qsat_final=" << qv_minus_qsat_final
            << " clipped_qc=" << investigation.clipped_qc
            << " clipped_qi=" << investigation.clipped_qi
            << " requested_delta_qc=" << investigation.requested_delta_qc
            << " requested_delta_qi=" << investigation.requested_delta_qi;
    }

// Motivation:
// The same component-wise PrecipFall budgets should close in terrain-like
// coordinates when the donor limiter and mass budgets are detJ-weighted.
TEST(SAMPhysicalProperties, PrecipFall_DetJWeightedComponentBudgetsClose)
{
    constexpr int nx = 2;
    constexpr int ny = 2;
    constexpr int nz = 4;
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real dt = amrex::Real(0.25);

    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(nx - 1, ny - 1, nz - 1));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(2.0, 2.0, 4.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 0)};
    amrex::Geometry geom(domain, &real_box, 0, is_periodic.data());

    amrex::BoxArray ba(domain);
    ba.maxSize(amrex::IntVect(1, 1, nz));
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(amrex::The_Pinned_Arena());
    cons.define(ba, dm, RhoQ6_comp + 1, 0);
    fill_precipfall_mixed_component_state(cons, pres_mbar);
    amrex::MultiFab detJ_cc(amrex::The_Pinned_Arena());
    detJ_cc.define(ba, dm, 1, 0);
    fill_precipfall_detj(detJ_cc);

    const PrecipFallBudgetSummary summary =
        run_precipfall_budget_case(geom, cons, MoistureType::SAM, dt, &detJ_cc);
    const amrex::Real rain_residual =
        summary.initial_rain_mass - summary.final_rain_mass - summary.rain_accum_mass;
    const amrex::Real snow_residual =
        summary.initial_snow_mass - summary.final_snow_mass - summary.snow_accum_mass;
    const amrex::Real graupel_residual =
        summary.initial_graupel_mass - summary.final_graupel_mass - summary.graupel_accum_mass;
    const amrex::Real total_residual = rain_residual + snow_residual + graupel_residual;

    EXPECT_NEAR(rain_residual,
                amrex::Real(0.0),
                property_accumulation_tol(nx * ny * nz,
                                          std::max(summary.initial_rain_mass, summary.rain_accum_mass)))
        << "rain residual=" << rain_residual
        << " snow residual=" << snow_residual
        << " graupel residual=" << graupel_residual
        << " total residual=" << total_residual
        << " dt=" << dt
        << " dz=" << geom.CellSize(2)
        << " n_substep=" << summary.diagnostics.n_substep
        << " max_raw_fluxes=(" << summary.diagnostics.max_raw_fluxes[0] << ","
        << summary.diagnostics.max_raw_fluxes[1] << "," << summary.diagnostics.max_raw_fluxes[2] << ")"
        << " max_limited_fluxes=(" << summary.diagnostics.max_limited_fluxes[0] << ","
        << summary.diagnostics.max_limited_fluxes[1] << "," << summary.diagnostics.max_limited_fluxes[2] << ")"
        << " limiter_active=" << summary.diagnostics.limiter_active
        << " min_q_after=(" << summary.min_qpr_after << "," << summary.min_qps_after << ","
        << summary.min_qpg_after << ")"
        << " detj_min=" << summary.diagnostics.detj_min
        << " detj_max=" << summary.diagnostics.detj_max;
    EXPECT_NEAR(snow_residual,
                amrex::Real(0.0),
                property_accumulation_tol(nx * ny * nz,
                                          std::max(summary.initial_snow_mass, summary.snow_accum_mass)));
    EXPECT_NEAR(graupel_residual,
                amrex::Real(0.0),
                property_accumulation_tol(nx * ny * nz,
                                          std::max(summary.initial_graupel_mass, summary.graupel_accum_mass)));
    EXPECT_NEAR(total_residual,
                amrex::Real(0.0),
                property_accumulation_tol(3 * nx * ny * nz,
                                          std::max({summary.initial_rain_mass + summary.initial_snow_mass + summary.initial_graupel_mass,
                                                    summary.rain_accum_mass + summary.snow_accum_mass + summary.graupel_accum_mass,
                                                    amrex::Real(1.0)})));
}

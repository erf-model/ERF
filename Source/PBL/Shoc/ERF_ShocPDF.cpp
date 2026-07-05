#include "ERF_ShocPDF.H"

#include "ERF_Constants.H"
#include "ERF_MicrophysicsUtils.H"

#include <AMReX_BLProfiler.H>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace amrex;

namespace
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_thl_tol () noexcept { return 1.0e-2_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_rt_tol () noexcept { return 1.0e-4_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_w_tol_sqd () noexcept { return 4.0e-4_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_w_thresh () noexcept { return 0.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_large_neg () noexcept { return -99999999.99_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_base_temp () noexcept { return 300.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_tl_min () noexcept { return 100.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_pdf_tmp () noexcept { return 0.4_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_sqrt_pdf_tmp () noexcept { return 0.77459666924148337704_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool shoc_do_thetal_skew () noexcept { return false; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real inv_sqrt_two_pi () noexcept { return 0.39894228040143267794_rt; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real clamp01 (Real x) { return shoc_clamp(x, 0.0_rt, 1.0_rt); }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real signed_denominator (Real value, Real eps = 1.0e-12_rt)
    {
        if (amrex::Math::abs(value) >= eps) {
            return value;
        }
        return std::copysign(eps, value == 0.0_rt ? 1.0_rt : value);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real weighted_linear_interp (Real x0, Real x1, Real y0, Real y1, Real x)
    {
        const Real denom = x1 - x0;
        if (amrex::Math::abs(denom) <= 1.0e-12_rt) {
            return 0.5_rt * (y0 + y1);
        }
        return y0 + (y1 - y0) * (x - x0) / denom;
    }

    void interpolate_iface_to_cc (const ShocColumnData& col,
                                  const FArrayBox& iface_in,
                                  FArrayBox& cc_out,
                                  Real min_thresh)
    {
        const auto zt = col.zt.const_array();
        const auto zi = col.zi.const_array();
        const auto in = iface_in.const_array();
        auto out = cc_out.array();

        const auto layout = col.layout;
        const Box cell_box(IntVect(0,0,0), IntVect(layout.ncell - 1, layout.nlev - 1, 0));
        ParallelFor(cell_box, [=] AMREX_GPU_DEVICE (int ic, int k, int) noexcept
        {
            const Real interp = weighted_linear_interp(zi(ic,k,0), zi(ic,k+1,0),
                                                       in(ic,k,0), in(ic,k+1,0),
                                                       zt(ic,k,0));
            out(ic,k,0) = amrex::max(min_thresh, interp);
        });
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void vv_parameters (Real w_first, Real w_sec, Real w3var,
                        Real& skew_w, Real& w1_1, Real& w1_2,
                        Real& w2_1, Real& w2_2, Real& a)
    {
        skew_w = 0.0;
        w1_1 = w_first;
        w1_2 = w_first;
        w2_1 = 0.0;
        w2_2 = 0.0;
        a = 0.5_rt;

        if (w_sec > shoc_w_tol_sqd()) {
            const Real w_sec_cu = w_sec * w_sec * w_sec;
            const Real one_minus_pdf = 1.0_rt - shoc_pdf_tmp();
            const Real one_minus_pdf_cu = one_minus_pdf * one_minus_pdf * one_minus_pdf;
            const Real skew_denom = std::sqrt(amrex::max(w_sec_cu, 1.0e-18_rt));
            skew_w = w3var / skew_denom;
            const Real skew_term = std::sqrt(1.0_rt /
                (4.0_rt * one_minus_pdf_cu + skew_w * skew_w));
            a = shoc_clamp(0.5_rt * (1.0_rt - skew_w * skew_term), 0.01_rt, 0.99_rt);
            w1_1 = std::sqrt((1.0_rt - a) / a) * shoc_sqrt_pdf_tmp();
            w1_2 = -std::sqrt(a / (1.0_rt - a)) * shoc_sqrt_pdf_tmp();
            w2_1 = shoc_pdf_tmp() * w_sec;
            w2_2 = shoc_pdf_tmp() * w_sec;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void thl_parameters (Real wthlsec, Real sqrtw2, Real sqrtthl,
                         Real thlsec, Real thl_first, Real w1_1, Real w1_2,
                         Real skew_w, Real a,
                         Real& thl1_1, Real& thl1_2, Real& thl2_1,
                         Real& thl2_2, Real& sqrtthl2_1, Real& sqrtthl2_2)
    {
        thl1_1 = thl_first;
        thl1_2 = thl_first;
        thl2_1 = 0.0;
        thl2_2 = 0.0;
        sqrtthl2_1 = 0.0;
        sqrtthl2_2 = 0.0;

        if (thlsec > shoc_thl_tol() * shoc_thl_tol() &&
            amrex::Math::abs(w1_2 - w1_1) > shoc_w_thresh() &&
            sqrtw2 > 0.0 && sqrtthl > 0.0) {
            const Real corr = shoc_clamp(wthlsec / (sqrtw2 * sqrtthl), -1.0_rt, 1.0_rt);
            const Real tmp1 = -corr / w1_1;
            const Real tmp2 = -corr / w1_2;
            const Real tsign = amrex::Math::abs(tmp1 - tmp2);
            Real skew_thl = 0.0_rt;
            if (shoc_do_thetal_skew()) {
                if (tsign > 0.4_rt) {
                    skew_thl = 1.2_rt * skew_w;
                } else if (tsign > 0.2_rt) {
                    skew_thl = ((1.2_rt * skew_w) / 0.2_rt) * (tsign - 0.2_rt);
                }
            }

            thl2_1 = shoc_clamp(
                (3.0_rt * tmp1 * (1.0_rt - a * tmp2 * tmp2 - (1.0_rt - a) * tmp1 * tmp1)
                 - (skew_thl - a * tmp2 * tmp2 * tmp2 - (1.0_rt - a) * tmp1 * tmp1 * tmp1)) /
                signed_denominator(3.0_rt * a * (tmp1 - tmp2)),
                0.0_rt, 100.0_rt) * thlsec;

            thl2_2 = shoc_clamp(
                (-3.0_rt * tmp2 * (1.0_rt - a * tmp2 * tmp2 - (1.0_rt - a) * tmp1 * tmp1)
                 + (skew_thl - a * tmp2 * tmp2 * tmp2 - (1.0_rt - a) * tmp1 * tmp1 * tmp1)) /
                signed_denominator(3.0_rt * (1.0_rt - a) * (tmp1 - tmp2)),
                0.0_rt, 100.0_rt) * thlsec;

            thl1_1 = tmp2 * sqrtthl + thl_first;
            thl1_2 = tmp1 * sqrtthl + thl_first;
            sqrtthl2_1 = std::sqrt(amrex::max(thl2_1, 0.0_rt));
            sqrtthl2_2 = std::sqrt(amrex::max(thl2_2, 0.0_rt));
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void qw_parameters (Real wqwsec, Real sqrtw2, Real skew_w, Real sqrtqt,
                        Real qwsec, Real w1_2, Real w1_1, Real qw_first, Real a,
                        Real& qw1_1, Real& qw1_2, Real& qw2_1, Real& qw2_2,
                        Real& sqrtqw2_1, Real& sqrtqw2_2)
    {
        qw1_1 = qw_first;
        qw1_2 = qw_first;
        qw2_1 = 0.0;
        qw2_2 = 0.0;
        sqrtqw2_1 = 0.0;
        sqrtqw2_2 = 0.0;

        if (qwsec > shoc_rt_tol() * shoc_rt_tol() &&
            amrex::Math::abs(w1_2 - w1_1) > shoc_w_thresh() &&
            sqrtw2 > 0.0 && sqrtqt > 0.0) {
            const Real corr = shoc_clamp(wqwsec / (sqrtw2 * sqrtqt), -1.0_rt, 1.0_rt);
            const Real tmp1 = -corr / w1_1;
            const Real tmp2 = -corr / w1_2;
            const Real tsign = amrex::Math::abs(tmp1 - tmp2);
            Real skew_qw = 0.0_rt;
            if (tsign > 0.4_rt) {
                skew_qw = 1.2_rt * skew_w;
            } else if (tsign > 0.2_rt) {
                skew_qw = ((1.2_rt * skew_w) / 0.2_rt) * (tsign - 0.2_rt);
            }

            qw2_1 = shoc_clamp(
                (3.0_rt * tmp1 * (1.0_rt - a * tmp2 * tmp2 - (1.0_rt - a) * tmp1 * tmp1)
                 - (skew_qw - a * tmp2 * tmp2 * tmp2 - (1.0_rt - a) * tmp1 * tmp1 * tmp1)) /
                signed_denominator(3.0_rt * a * (tmp1 - tmp2)),
                0.0_rt, 100.0_rt) * qwsec;

            qw2_2 = shoc_clamp(
                (-3.0_rt * tmp2 * (1.0_rt - a * tmp2 * tmp2 - (1.0_rt - a) * tmp1 * tmp1)
                 + (skew_qw - a * tmp2 * tmp2 * tmp2 - (1.0_rt - a) * tmp1 * tmp1 * tmp1)) /
                signed_denominator(3.0_rt * (1.0_rt - a) * (tmp1 - tmp2)),
                0.0_rt, 100.0_rt) * qwsec;

            qw1_1 = tmp2 * sqrtqt + qw_first;
            qw1_2 = tmp1 * sqrtqt + qw_first;
            sqrtqw2_1 = std::sqrt(amrex::max(qw2_1, 0.0_rt));
            sqrtqw2_2 = std::sqrt(amrex::max(qw2_2, 0.0_rt));
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real inplume_correlation (Real sqrtqw2_1, Real sqrtthl2_1, Real a,
                              Real sqrtqw2_2, Real sqrtthl2_2,
                              Real qwthlsec, Real qw1_1, Real qw_first,
                              Real thl1_1, Real thl_first,
                              Real qw1_2, Real thl1_2)
    {
        const Real denom = a * sqrtqw2_1 * sqrtthl2_1 + (1.0_rt - a) * sqrtqw2_2 * sqrtthl2_2;
        if (amrex::Math::abs(denom) <= 1.0e-18_rt) {
            return 0.0_rt;
        }
        return shoc_clamp((qwthlsec
                          - a * (qw1_1 - qw_first) * (thl1_1 - thl_first)
                          - (1.0_rt - a) * (qw1_2 - qw_first) * (thl1_2 - thl_first)) / denom,
                          -1.0_rt, 1.0_rt);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real compute_temperature (Real thl1, Real pval)
    {
        const Real pressure_ratio = p_0 / amrex::max(pval, 1.0e-12_rt);
        return thl1 / std::exp((R_d / Cp_d) * std::log(pressure_ratio));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void compute_qs_beta (Real temp, Real pval, Real& qs, Real& beta)
    {
        qs = 0.0_rt;
        erf_qsatw(temp, 0.01_rt * pval, qs);
        qs = amrex::max(0.0_rt, qs);
        beta = (R_d / R_v) * (L_v / (R_d * temp)) * (L_v / (Cp_d * temp));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void compute_s_terms (Real qw1, Real qs, Real beta, Real pval,
                          Real thl2, Real qw2, Real sqrtthl2, Real sqrtqw2,
                          Real r_qwthl, Real& s, Real& std_s, Real& qn, Real& cfrac)
    {
        const Real beta_qs = 1.0_rt + beta * qs;
        const Real pressure_ratio = pval / p_0;
        const Real pressure_term = std::exp((R_d / Cp_d) * std::log(pressure_ratio));
        const Real cthl = ((1.0_rt + beta * qw1) / (beta_qs * beta_qs))
                        * (Cp_d / L_v) * beta * qs * pressure_term;
        const Real cqt = 1.0_rt / (1.0_rt + beta * qs);
        const Real variance = cthl * cthl * thl2
                            + cqt * cqt * qw2
                            - 2.0_rt * cthl * sqrtthl2 * cqt * sqrtqw2 * r_qwthl;
        std_s = std::sqrt(amrex::max(0.0_rt, variance));
        s = qw1 - qs * ((1.0_rt + beta * qw1) / (1.0_rt + beta * qs));

        cfrac = 0.0_rt;
        qn = 0.0_rt;
        if (std_s > std::sqrt(std::numeric_limits<Real>::min()) * 100.0_rt) {
            cfrac = 0.5_rt * (1.0_rt + std::erf(s / (std::sqrt(2.0_rt) * std_s)));
            if (cfrac != 0.0_rt) {
                qn = s * cfrac + std_s * inv_sqrt_two_pi() * std::exp(-0.5_rt * (s / std_s) * (s / std_s));
            }
        } else if (s > 0.0_rt) {
            cfrac = 1.0_rt;
            qn = s;
        }

        if (qn <= 0.0_rt) {
            cfrac = 0.0_rt;
            qn = 0.0_rt;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real compute_buoyancy_flux (Real wthlsec, Real wqwsec, Real pval, Real wqls)
    {
        const Real epsterm = R_d / R_v;
        const Real pressure_ratio = p_0 / amrex::max(pval, 1.0e-12_rt);
        const Real pressure_term = std::exp((R_d / Cp_d) * std::log(pressure_ratio));
        return wthlsec
                 + ((1.0_rt - epsterm) / epsterm) * shoc_base_temp() * wqwsec
             + ((L_v / Cp_d) * pressure_term
                     - (1.0_rt / epsterm) * shoc_base_temp()) * wqls;
    }
}

void
ShocPDF::diagnose_pdf (ShocColumnData& col,
                       const ShocRuntimeOptions& opts,
                       Real dt)
{
    // Interim native-SHOC contract: this PDF remains liquid-cloud
    // macrophysics. It may carry pre-existing ice through the thermodynamic
    // state, but it does not create or repartition cloud ice here.
    FArrayBox w3_zt, thl_sec_zt, qw_sec_zt, wthl_zt, wqw_zt, qwthl_zt;
    const auto layout = col.layout;
    const Box cell_box(IntVect(0,0,0), IntVect(layout.ncell - 1, layout.nlev - 1, 0));
    w3_zt.resize(cell_box, 1, The_Async_Arena());
    thl_sec_zt.resize(cell_box, 1, The_Async_Arena());
    qw_sec_zt.resize(cell_box, 1, The_Async_Arena());
    wthl_zt.resize(cell_box, 1, The_Async_Arena());
    wqw_zt.resize(cell_box, 1, The_Async_Arena());
    qwthl_zt.resize(cell_box, 1, The_Async_Arena());

    {
        BL_PROFILE("SHOC::advance::pdf::interpolate");
        interpolate_iface_to_cc(col, col.w3, w3_zt, shoc_large_neg());
        interpolate_iface_to_cc(col, col.thl_sec, thl_sec_zt, 0.0_rt);
        interpolate_iface_to_cc(col, col.qw_sec, qw_sec_zt, 0.0_rt);
        interpolate_iface_to_cc(col, col.wthl_sec, wthl_zt, shoc_large_neg());
        interpolate_iface_to_cc(col, col.wqw_sec, wqw_zt, shoc_large_neg());
        interpolate_iface_to_cc(col, col.qwthl_sec, qwthl_zt, shoc_large_neg());
    }

    auto shoc_cldfrac = col.shoc_cldfrac.array();
    auto shoc_ql = col.shoc_ql.array();
    auto shoc_ql2 = col.shoc_ql2.array();
    auto shoc_cond = col.shoc_cond.array();
    auto shoc_evap = col.shoc_evap.array();
    auto wqls_sec = col.wqls_sec.array();
    auto wthv_sec = col.wthv_sec.array();

    const auto thetal = col.thetal.const_array();
    const auto qw = col.qw.const_array();
    const auto w_field = col.w.const_array();
    const auto p_mid = col.p_mid.const_array();
    const auto thl_sec = thl_sec_zt.const_array();
    const auto qw_sec = qw_sec_zt.const_array();
    const auto qwthl_sec = qwthl_zt.const_array();
    const auto wthl_sec = wthl_zt.const_array();
    const auto wqw_sec = wqw_zt.const_array();
    const auto w_sec = col.w_sec.const_array();
    const auto w3 = w3_zt.const_array();

    {
        BL_PROFILE("SHOC::advance::pdf::main_kernel");
        ParallelFor(cell_box, [=] AMREX_GPU_DEVICE (int ic, int k, int) noexcept
        {
            const Real thl_first = thetal(ic,k,0);
            const Real qw_first = qw(ic,k,0);
            const Real w_first = w_field(ic,k,0);
            const Real w2sec = amrex::max(w_sec(ic,k,0), 0.0_rt);
            const Real sqrtw2 = std::sqrt(w2sec);
            const Real sqrtthl = amrex::max(shoc_thl_tol(), std::sqrt(amrex::max(thl_sec(ic,k,0), 0.0_rt)));
            const Real sqrtqt = amrex::max(shoc_rt_tol(), std::sqrt(amrex::max(qw_sec(ic,k,0), 0.0_rt)));

            Real skew_w = 0.0, w1_1 = w_first, w1_2 = w_first, w2_1 = 0.0, w2_2 = 0.0, a = 0.5;
            vv_parameters(w_first, w2sec, w3(ic,k,0), skew_w, w1_1, w1_2, w2_1, w2_2, a);

            Real thl1_1 = thl_first, thl1_2 = thl_first, thl2_1 = 0.0, thl2_2 = 0.0;
            Real sqrtthl2_1 = 0.0, sqrtthl2_2 = 0.0;
            thl_parameters(wthl_sec(ic,k,0), sqrtw2, sqrtthl, amrex::max(thl_sec(ic,k,0), 0.0_rt),
                           thl_first, w1_1, w1_2, skew_w, a,
                           thl1_1, thl1_2, thl2_1, thl2_2, sqrtthl2_1, sqrtthl2_2);

            Real qw1_1 = qw_first, qw1_2 = qw_first, qw2_1 = 0.0, qw2_2 = 0.0;
            Real sqrtqw2_1 = 0.0, sqrtqw2_2 = 0.0;
            qw_parameters(wqw_sec(ic,k,0), sqrtw2, skew_w, sqrtqt, amrex::max(qw_sec(ic,k,0), 0.0_rt),
                          w1_2, w1_1, qw_first, a,
                          qw1_1, qw1_2, qw2_1, qw2_2, sqrtqw2_1, sqrtqw2_2);

            w1_1 = w1_1 * sqrtw2 + w_first;
            w1_2 = w1_2 * sqrtw2 + w_first;

            const Real r_qwthl_1 = inplume_correlation(sqrtqw2_1, sqrtthl2_1, a,
                                                       sqrtqw2_2, sqrtthl2_2,
                                                       qwthl_sec(ic,k,0), qw1_1, qw_first,
                                                       thl1_1, thl_first, qw1_2, thl1_2);

            Real Tl1_1 = amrex::max(shoc_tl_min(), compute_temperature(thl1_1, p_mid(ic,k,0)));
            Real Tl1_2 = amrex::max(shoc_tl_min(), compute_temperature(thl1_2, p_mid(ic,k,0)));
            Real qs1 = 0.0_rt, beta1 = 0.0_rt, qs2 = 0.0_rt, beta2 = 0.0_rt;
            compute_qs_beta(Tl1_1, p_mid(ic,k,0), qs1, beta1);
            compute_qs_beta(Tl1_2, p_mid(ic,k,0), qs2, beta2);

            Real s1 = 0.0_rt, std_s1 = 0.0_rt, qn1 = 0.0_rt, c1 = 0.0_rt;
            Real s2 = 0.0_rt, std_s2 = 0.0_rt, qn2 = 0.0_rt, c2 = 0.0_rt;
            compute_s_terms(qw1_1, qs1, beta1, p_mid(ic,k,0), thl2_1, qw2_1,
                            sqrtthl2_1, sqrtqw2_1, r_qwthl_1, s1, std_s1, qn1, c1);

            if (qw1_1 == qw1_2 && thl2_1 == thl2_2 && qs1 == qs2) {
                s2 = s1;
                std_s2 = std_s1;
                qn2 = qn1;
                c2 = c1;
            } else {
                compute_s_terms(qw1_2, qs2, beta2, p_mid(ic,k,0), thl2_2, qw2_2,
                                sqrtthl2_2, sqrtqw2_2, r_qwthl_1, s2, std_s2, qn2, c2);
            }

            const Real ql1 = amrex::min(qn1, qw1_1);
            const Real ql2 = amrex::min(qn2, qw1_2);
            const Real old_ql = shoc_ql(ic,k,0);
            const Real cldfrac = amrex::min(1.0_rt, a * c1 + (1.0_rt - a) * c2);
            const Real ql = amrex::max(0.0_rt, a * ql1 + (1.0_rt - a) * ql2);
            const Real ql2_var = amrex::max(0.0_rt, a * (s1 * ql1 + c1 * std_s1 * std_s1)
                                               + (1.0_rt - a) * (s2 * ql2 + c2 * std_s2 * std_s2)
                                               - ql * ql);
            const Real wqls = a * ((w1_1 - w_first) * ql1) + (1.0_rt - a) * ((w1_2 - w_first) * ql2);

            shoc_cldfrac(ic,k,0) = clamp01(cldfrac);
            shoc_ql(ic,k,0) = ql;
            shoc_ql2(ic,k,0) = ql2_var;
            wqls_sec(ic,k,0) = wqls;
            wthv_sec(ic,k,0) = compute_buoyancy_flux(wthl_sec(ic,k,0), wqw_sec(ic,k,0),
                                                     p_mid(ic,k,0), wqls);

            if (opts.extra_shoc_diags) {
                const Real ql_change = ql - old_ql;
                shoc_cond(ic,k,0) = amrex::max(0.0_rt, ql_change / amrex::max(dt, 1.0e-12_rt));
                shoc_evap(ic,k,0) = amrex::max(0.0_rt, -ql_change / amrex::max(dt, 1.0e-12_rt));
            } else {
                shoc_cond(ic,k,0) = 0.0_rt;
                shoc_evap(ic,k,0) = 0.0_rt;
            }
        });
    }
}

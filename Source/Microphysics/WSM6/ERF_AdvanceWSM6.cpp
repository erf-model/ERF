#include "ERF_WSM6.H"
#include "ERF_WSM6_Fortran_Interface.H"
#include <AMReX_Reduce.H>
#include <algorithm>
#include <array>
#include <cctype>
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

using namespace amrex;


// ---------------------------------------------------------------
// WSM6 device-callable free functions (Rule 16, Rule 18)
// Statement functions from mp_wsm6_run declaration section
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_cpmcal (Real x, Real qmin_arg, Real cpd_arg, Real cpv_arg) {
    return cpd_arg*(Real(1.0)-amrex::max(x,qmin_arg))
          +amrex::max(x,qmin_arg)*cpv_arg;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_xlcal (Real x, Real xlv0_arg, Real xlv1_arg, Real t0c_arg) {
    return xlv0_arg - xlv1_arg*(x - t0c_arg);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_diffus (Real x, Real y) {
    return Real(8.794e-5)*std::exp(std::log(x)*Real(1.81))/y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_viscos (Real x, Real y) {
    return Real(1.496e-6)*(x*std::sqrt(x))/(x+Real(120.0))/y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_xka (Real x, Real y) {
    return Real(1.414e3)*wsm6_viscos(x,y)*y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_diffac (Real a, Real b, Real c, Real d, Real e,
                   Real rv_arg) {
    return d*a*a/(wsm6_xka(c,d)*rv_arg*c*c)
          +Real(1.0)/(e*wsm6_diffus(c,b));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_venfac (Real a, Real b, Real c, Real den0_arg) {
    return std::exp(std::log(wsm6_viscos(b,c)/wsm6_diffus(b,a))
                   *Real(0.3333333))
          /std::sqrt(wsm6_viscos(b,c))
          *std::sqrt(std::sqrt(den0_arg/c));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_conden (Real a, Real b, Real c, Real d, Real e,
                   Real qmin_arg, Real rv_arg) {
    return (amrex::max(b,qmin_arg)-c)
          /(Real(1.0)+d*d/(rv_arg*e)*c/(a*a));
}

// ---------------------------------------------------------------
// Slope parameter lambda functions (statement functions from
// slope_wsm6, slope_rain, slope_snow, slope_graup)
// pidn0r, pidn0s, pidn0g are class constexpr members from ERF_WSM6.H
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_lamdar (Real x, Real y, Real pidn0r_arg) {
    return std::sqrt(std::sqrt(pidn0r_arg/(x*y)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_lamdas (Real x, Real y, Real z, Real pidn0s_arg) {
    return std::sqrt(std::sqrt(pidn0s_arg*z/(x*y)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_lamdag (Real x, Real y, Real pidn0g_arg) {
    return std::sqrt(std::sqrt(pidn0g_arg/(x*y)));
}

// ---------------------------------------------------------------
// Full slope subroutine device functions (Rule 18)
// Each takes single-cell scalar inputs, returns slope params
// by reference — loop over (i,j,k) is provided by ParallelFor
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wsm6_slope_rain_cell (Real qr, Real den, Real denfac,
                            Real pidn0r_arg,
                            Real qcrmin_arg, Real rslopermax_arg,
                            Real rsloperbmax_arg, Real rsloper2max_arg,
                            Real rsloper3max_arg, Real bvtr_arg,
                            Real pvtr_arg,
                            Real& rslope, Real& rslopeb,
                            Real& rslope2, Real& rslope3, Real& vt)
{
    if (qr <= qcrmin_arg) {
        rslope  = rslopermax_arg;
        rslopeb = rsloperbmax_arg;
        rslope2 = rsloper2max_arg;
        rslope3 = rsloper3max_arg;
    } else {
        rslope  = Real(1.0)/wsm6_lamdar(qr,den,pidn0r_arg);
        rslopeb = std::pow(rslope,bvtr_arg);
        rslope2 = rslope*rslope;
        rslope3 = rslope2*rslope;
    }
    vt = pvtr_arg*rslopeb*denfac;
    if (qr <= Real(0.0)) vt = Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wsm6_slope_snow_cell (Real qs, Real den, Real denfac, Real t,
                            Real pidn0s_arg, Real alpha_arg,
                            Real n0smax_arg, Real n0s_arg,
                            Real t0c_arg, Real qcrmin_arg,
                            Real rslopesmax_arg, Real rslopesbmax_arg,
                            Real rslopes2max_arg, Real rslopes3max_arg,
                            Real bvts_arg, Real pvts_arg,
                            Real& rslope, Real& rslopeb,
                            Real& rslope2, Real& rslope3, Real& vt,
                            Real& n0sfac)
{
    Real supcol = t0c_arg - t;
    n0sfac = amrex::max(amrex::min(std::exp(alpha_arg*supcol),
                                    n0smax_arg/n0s_arg), Real(1.0));
    if (qs <= qcrmin_arg) {
        rslope  = rslopesmax_arg;
        rslopeb = rslopesbmax_arg;
        rslope2 = rslopes2max_arg;
        rslope3 = rslopes3max_arg;
    } else {
        rslope  = Real(1.0)/wsm6_lamdas(qs,den,n0sfac,pidn0s_arg);
        rslopeb = std::pow(rslope,bvts_arg);
        rslope2 = rslope*rslope;
        rslope3 = rslope2*rslope;
    }
    vt = pvts_arg*rslopeb*denfac;
    if (qs <= Real(0.0)) vt = Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wsm6_slope_graup_cell (Real qg, Real den, Real denfac,
                             Real pidn0g_arg, Real qcrmin_arg,
                             Real rslopegmax_arg, Real rslopegbmax_arg,
                             Real rslopeg2max_arg, Real rslopeg3max_arg,
                             Real bvtg_arg, Real pvtg_arg,
                             Real& rslope, Real& rslopeb,
                             Real& rslope2, Real& rslope3, Real& vt)
{
    if (qg <= qcrmin_arg) {
        rslope  = rslopegmax_arg;
        rslopeb = rslopegbmax_arg;
        rslope2 = rslopeg2max_arg;
        rslope3 = rslopeg3max_arg;
    } else {
        rslope  = Real(1.0)/wsm6_lamdag(qg,den,pidn0g_arg);
        rslopeb = std::pow(rslope,bvtg_arg);
        rslope2 = rslope*rslope;
        rslope3 = rslope2*rslope;
    }
    vt = pvtg_arg*rslopeb*denfac;
    if (qg <= Real(0.0)) vt = Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wsm6_nislfv_rain_plm (int im, int km,
                           const Real* denl, const Real* denfacl,
                           const Real* tkl, const Real* dzl,
                           Real* wwl, Real* rql,
                           Real* precip, Real dt,
                           int id, int iter,
                           int /*line_check_dbg*/,
                           Real* dbg_zi = nullptr,
                           Real* dbg_za = nullptr,
                           Real* dbg_dza = nullptr,
                           Real* dbg_wi = nullptr,
                           Real* dbg_ww = nullptr,
                           Real* dbg_wa = nullptr,
                           Real* dbg_was = nullptr,
                           Real* dbg_qa = nullptr,
                           Real* dbg_qmi = nullptr,
                           Real* dbg_qpi = nullptr,
                           Real* dbg_kb_before_backstep = nullptr,
                           Real* dbg_kt_before_backstep = nullptr,
                           Real* dbg_kb_after_backstep = nullptr,
                           Real* dbg_kt_after_backstep = nullptr,
                           Real* dbg_kb_after_search = nullptr,
                           Real* dbg_kt_after_search = nullptr,
                           Real* dbg_zsum = nullptr,
                           Real* dbg_qsum = nullptr,
                           Real* dbg_qn = nullptr,
                           Real* dbg_denqrs1_after_kernel = nullptr)
{
    static_cast<void>(id);

    if (km > WSM6_MAX_LEVELS) return;

    constexpr Real pi = Real(3.141592653589793238462643383279502884);
    auto rgmma = [](Real x) -> Real {
        // Match Fortran rgmma() semantics used in WSM6 rain setup.
        if (x == Real(1.0)) return Real(0.0);
        constexpr Real euler = Real(0.577215664901532);
        Real rg = x * std::exp(euler * x);
        for (int ii = 1; ii <= 10000; ++ii) {
            const Real y = static_cast<Real>(ii);
            rg = rg * (Real(1.0) + x / y) * std::exp(-x / y);
        }
        return Real(1.0) / rg;
    };

    const Real pidn0r = pi * Real(rhoh2o) * WSM6::n0r;
    const Real rslopermax = Real(1.0) / WSM6::lamdarmax;
    const Real rsloperbmax = std::pow(rslopermax, WSM6::bvtr);
    const Real rsloper2max = rslopermax * rslopermax;
    const Real rsloper3max = rsloper2max * rslopermax;
    const Real pvtr = WSM6::avtr * rgmma(Real(4.0) + WSM6::bvtr) / Real(6.0);

    for (int i = 0; i < im; ++i) {
        Real dz[WSM6_MAX_LEVELS];
        Real ww[WSM6_MAX_LEVELS];
        Real qq[WSM6_MAX_LEVELS];
        Real wd[WSM6_MAX_LEVELS];
        Real wa[WSM6_MAX_LEVELS];
        Real was[WSM6_MAX_LEVELS];
        Real den[WSM6_MAX_LEVELS];
        Real denfac[WSM6_MAX_LEVELS];
        Real tk[WSM6_MAX_LEVELS];
        Real qn[WSM6_MAX_LEVELS];
        Real qr[WSM6_MAX_LEVELS];
        Real tmp[WSM6_MAX_LEVELS];
        Real tmp1[WSM6_MAX_LEVELS];
        Real tmp2[WSM6_MAX_LEVELS];
        Real tmp3[WSM6_MAX_LEVELS];
        Real wi[WSM6_MAX_LEVELS + 1];
        Real zi[WSM6_MAX_LEVELS + 1];
        Real za[WSM6_MAX_LEVELS + 1];
        Real dza[WSM6_MAX_LEVELS + 1];
        Real qa[WSM6_MAX_LEVELS + 1];
        Real qmi[WSM6_MAX_LEVELS + 1];
        Real qpi[WSM6_MAX_LEVELS + 1];

        const bool emit_search_dbg = (dbg_qn != nullptr);
        Real allold = Real(0.0);
        for (int k = 0; k < km; ++k) {
            const int idx = i * km + k;
            dz[k] = dzl[idx];
            qq[k] = rql[idx];
            ww[k] = wwl[idx];
            wd[k] = ww[k];
            den[k] = denl[idx];
            denfac[k] = denfacl[idx];
            tk[k] = tkl[idx];
            allold += qq[k];
        }

        precip[i] = Real(0.0);
        if (allold <= Real(0.0)) {
            continue;
        }

        zi[0] = Real(0.0);
        for (int k = 0; k < km; ++k) {
            zi[k + 1] = zi[k] + dz[k];
            if (emit_search_dbg) {
                dbg_zi[k] = zi[k];
            }
        }

        auto update_wind_and_state = [&](void) {
            wi[0] = ww[0];
            wi[km] = ww[km - 1];
            for (int k = 1; k < km; ++k) {
                wi[k] = (ww[k] * dz[k - 1] + ww[k - 1] * dz[k]) / (dz[k - 1] + dz[k]);
            }

            wi[0] = ww[0];
            wi[1] = Real(0.5) * (ww[1] + ww[0]);
            for (int k = 2; k < km - 1; ++k) {
                wi[k] = Real(9.0) / Real(16.0) * (ww[k] + ww[k - 1])
                      - Real(1.0) / Real(16.0) * (ww[k + 1] + ww[k - 2]);
            }
            if (km > 1) {
                wi[km - 1] = Real(0.5) * (ww[km - 1] + ww[km - 2]);
                wi[km] = ww[km - 1];
            }

            for (int k = 1; k < km; ++k) {
                if (ww[k] == Real(0.0)) wi[k] = ww[k - 1];
            }

            const Real con1 = Real(0.05);
            for (int k = km - 1; k >= 0; --k) {
                const Real decfl = (wi[k + 1] - wi[k]) * dt / dz[k];
                if (decfl > con1) {
                    wi[k] = wi[k + 1] - con1 * dz[k] / dt;
                }
            }

            for (int k = 0; k <= km; ++k) {
                za[k] = zi[k] - wi[k] * dt;
            }

            for (int k = 0; k < km; ++k) {
                dza[k] = za[k + 1] - za[k];
                if (dza[k] <= Real(0.0)) dza[k] = dz[k]; // divergence guard
            }
            dza[km] = zi[km] - za[km]; // Fortran: dza(km+1) = zi(km+1)-za(km+1)
            if (dza[km] <= Real(0.0)) dza[km] = dz[km > 0 ? km - 1 : 0];
            for (int k = 0; k < km; ++k) {
                qa[k] = qq[k] * dz[k] / dza[k];
                qr[k] = qa[k] / den[k];
                if (emit_search_dbg) {
                    dbg_za[k] = za[k];
                    dbg_dza[k] = dza[k];
                    dbg_wi[k] = wi[k];
                    dbg_ww[k] = ww[k];
                    dbg_wa[k] = wa[k];
                    dbg_was[k] = was[k];
                    dbg_qa[k] = qa[k];
                }
            }
            qa[km] = Real(0.0);
        };

        update_wind_and_state();

        if (iter > 0) {
            for (int k = 0; k < km; ++k) {
                wsm6_slope_rain_cell(qr[k], den[k], denfac[k], pidn0r,
                                     WSM6::qcrmin,
                                     rslopermax, rsloperbmax, rsloper2max,
                                     rsloper3max, WSM6::bvtr, pvtr,
                                     tmp[k], tmp1[k], tmp2[k], tmp3[k], wa[k]);
            }
            for (int k = 0; k < km; ++k) {
                ww[k] = Real(0.5) * (wd[k] + wa[k]);
                was[k] = wa[k];
            }
            update_wind_and_state();
        }

        for (int k = 1; k < km; ++k) {
            const Real dip = (qa[k + 1] - qa[k]) / (dza[k + 1] + dza[k]);
            const Real dim = (qa[k] - qa[k - 1]) / (dza[k - 1] + dza[k]);
            if (dip * dim <= Real(0.0)) {
                qmi[k] = qa[k];
                qpi[k] = qa[k];
            } else {
                qpi[k] = qa[k] + Real(0.5) * (dip + dim) * dza[k];
                qmi[k] = Real(2.0) * qa[k] - qpi[k];
                if (qpi[k] < Real(0.0) || qmi[k] < Real(0.0)) {
                    qpi[k] = qa[k];
                    qmi[k] = qa[k];
                }
            }
            if (emit_search_dbg) {
                dbg_qmi[k] = qmi[k];
                dbg_qpi[k] = qpi[k];
            }
        }
        qpi[0] = qa[0];
        qmi[0] = qa[0];
        qmi[km] = qa[km];
        qpi[km] = qa[km];
        if (emit_search_dbg) {
            dbg_qmi[0] = qmi[0];
            dbg_qpi[0] = qpi[0];
        }

        for (int k = 0; k < km; ++k) {
            qn[k] = Real(0.0);
            if (emit_search_dbg) {
                dbg_kb_before_backstep[k] = Real(-1.0);
                dbg_kt_before_backstep[k] = Real(-1.0);
                dbg_kb_after_backstep[k] = Real(-1.0);
                dbg_kt_after_backstep[k] = Real(-1.0);
                dbg_kb_after_search[k] = Real(-1.0);
                dbg_kt_after_search[k] = Real(-1.0);
                dbg_zsum[k] = Real(0.0);
                dbg_qsum[k] = Real(0.0);
                dbg_qn[k] = Real(0.0);
            }
        }

        int kb = 0;
        int kt = 0;
        for (int k = 0; k < km; ++k) {
            if (emit_search_dbg) {
                dbg_kb_before_backstep[k] = Real(kb);
                dbg_kt_before_backstep[k] = Real(kt);
            }
            kb = amrex::max(kb - 1, 0);
            kt = amrex::max(kt - 1, 0);
            if (emit_search_dbg) {
                dbg_kb_after_backstep[k] = Real(kb);
                dbg_kt_after_backstep[k] = Real(kt);
            }

            if (zi[k] >= za[km]) {
                break;
            }

            for (int kk = kb; kk < km; ++kk) {
                if (zi[k] <= za[kk + 1]) {
                    kb = kk;
                    break;
                }
            }

            for (int kk = kt; kk < km; ++kk) {
                if (zi[k + 1] <= za[kk]) {
                    kt = kk;
                    break;
                }
            }
            kt = amrex::max(kt - 1, 0);
            if (emit_search_dbg) {
                dbg_kb_after_search[k] = Real(kb);
                dbg_kt_after_search[k] = Real(kt);
            }

            if (kt == kb) {
                const Real tl = (zi[k] - za[kb]) / dza[kb];
                const Real th = (zi[k + 1] - za[kb]) / dza[kb];
                const Real tl2 = tl * tl;
                const Real th2 = th * th;
                const Real qqd = Real(0.5) * (qpi[kb] - qmi[kb]);
                const Real qqh = qqd * th2 + qmi[kb] * th;
                const Real qql = qqd * tl2 + qmi[kb] * tl;
                qn[k] = (qqh - qql) / (th - tl);
                if (emit_search_dbg) {
                    dbg_zsum[k] = dza[kb];
                    dbg_qsum[k] = qn[k] * dza[kb];
                }
            } else if (kt > kb) {
                const Real tl = (zi[k] - za[kb]) / dza[kb];
                const Real tl2 = tl * tl;
                const Real qqd = Real(0.5) * (qpi[kb] - qmi[kb]);
                const Real qql = qqd * tl2 + qmi[kb] * tl;
                const Real dql = qa[kb] - qql;
                Real zsum = (Real(1.0) - tl) * dza[kb];
                Real qsum = dql * dza[kb];
                if (kt - kb > 1) {
                    for (int m = kb + 1; m < kt; ++m) {
                        zsum += dza[m];
                        qsum += qa[m] * dza[m];
                    }
                }
                const Real th = (zi[k + 1] - za[kt]) / dza[kt];
                const Real th2 = th * th;
                const Real dqh = Real(0.5) * (qpi[kt] - qmi[kt]) * th2 + qmi[kt] * th;
                zsum += th * dza[kt];
                qsum += dqh * dza[kt];
                qn[k] = qsum / zsum;
                if (emit_search_dbg) {
                    dbg_zsum[k] = zsum;
                    dbg_qsum[k] = qsum;
                }
            }
            if (emit_search_dbg) {
                dbg_qn[k] = qn[k];
            }
        }

        for (int k = 0; k < km; ++k) {
            if (za[k] < Real(0.0) && za[k + 1] < Real(0.0)) {
                precip[i] += qa[k] * dza[k];
            } else if (za[k] < Real(0.0) && za[k + 1] >= Real(0.0)) {
                precip[i] += qa[k] * (Real(0.0) - za[k]);
                break;
            } else {
                break;
            }
        }

        for (int k = 0; k < km; ++k) {
            rql[i * km + k] = qn[k];
            wwl[i * km + k] = ww[k];
            if (emit_search_dbg) {
                dbg_denqrs1_after_kernel[k] = qn[k];
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wsm6_nislfv_rain_plm6 (int im, int km,
                            const Real* denl, const Real* denfacl,
                            const Real* tkl, const Real* dzl,
                            Real* wwl, Real* rql, Real* rql2,
                            Real* precip1, Real* precip2, Real dt,
                            int id, int iter)
{
    static_cast<void>(id);

    if (km > WSM6_MAX_LEVELS) return;

    constexpr Real pi = Real(3.141592653589793238462643383279502884);
    auto rgmma = [](Real x) -> Real {
        if (x == Real(1.0)) return Real(0.0);
        constexpr Real euler = Real(0.577215664901532);
        Real rg = x * std::exp(euler * x);
        for (int ii = 1; ii <= 10000; ++ii) {
            const Real y = static_cast<Real>(ii);
            rg = rg * (Real(1.0) + x / y) * std::exp(-x / y);
        }
        return Real(1.0) / rg;
    };

    // Match Fortran mp_wsm6_init graupel-mode coefficient definitions.
    const Real dens = WSM6::dens_snow;
    const Real n0g = Real(4.0e6);
    const Real deng = Real(500.0);
    const Real avtg = Real(330.0);
    const Real bvtg = Real(0.8);
    const Real lamdagmax = Real(6.0e4);

    const Real pidn0s = pi * dens * WSM6::n0s;
    const Real pidn0g = pi * deng * n0g;
    const Real rslopesmax = Real(1.0) / WSM6::lamdasmax;
    const Real rslopesbmax = std::pow(rslopesmax, WSM6::bvts);
    const Real rslopes2max = rslopesmax * rslopesmax;
    const Real rslopes3max = rslopes2max * rslopesmax;
    const Real rslopegmax = Real(1.0) / lamdagmax;
    const Real rslopegbmax = std::pow(rslopegmax, bvtg);
    const Real rslopeg2max = rslopegmax * rslopegmax;
    const Real rslopeg3max = rslopeg2max * rslopegmax;
    const Real pvts = WSM6::avts * rgmma(Real(4.0) + WSM6::bvts) / Real(6.0);
    const Real pvtg = avtg * rgmma(Real(4.0) + bvtg) / Real(6.0);

    for (int i = 0; i < im; ++i) {
        Real dz[WSM6_MAX_LEVELS];
        Real ww[WSM6_MAX_LEVELS];
        Real qq[WSM6_MAX_LEVELS];
        Real qq2[WSM6_MAX_LEVELS];
        Real wd[WSM6_MAX_LEVELS];
        Real wa[WSM6_MAX_LEVELS];
        Real wa2[WSM6_MAX_LEVELS];
        Real was[WSM6_MAX_LEVELS];
        Real den[WSM6_MAX_LEVELS];
        Real denfac[WSM6_MAX_LEVELS];
        Real tk[WSM6_MAX_LEVELS];
        Real qn[WSM6_MAX_LEVELS];
        Real qn2[WSM6_MAX_LEVELS];
        Real qr[WSM6_MAX_LEVELS];
        Real qr2[WSM6_MAX_LEVELS];
        Real tmp[WSM6_MAX_LEVELS];
        Real tmp1[WSM6_MAX_LEVELS];
        Real tmp2[WSM6_MAX_LEVELS];
        Real tmp3[WSM6_MAX_LEVELS];
        Real wi[WSM6_MAX_LEVELS + 1];
        Real zi[WSM6_MAX_LEVELS + 1];
        Real za[WSM6_MAX_LEVELS + 1];
        Real dza[WSM6_MAX_LEVELS + 1];
        Real qa[WSM6_MAX_LEVELS + 1];
        Real qa2[WSM6_MAX_LEVELS + 1];
        Real qmi[WSM6_MAX_LEVELS + 1];
        Real qpi[WSM6_MAX_LEVELS + 1];

        Real allold = Real(0.0);
        for (int k = 0; k < km; ++k) {
            const int idx = i * km + k;
            dz[k] = dzl[idx];
            qq[k] = rql[idx];
            qq2[k] = rql2[idx];
            ww[k] = wwl[idx];
            wd[k] = ww[k];
            den[k] = denl[idx];
            denfac[k] = denfacl[idx];
            tk[k] = tkl[idx];
            allold += qq[k] + qq2[k];
        }

        precip1[i] = Real(0.0);
        precip2[i] = Real(0.0);
        if (allold <= Real(0.0)) {
            continue;
        }

        zi[0] = Real(0.0);
        for (int k = 0; k < km; ++k) {
            zi[k + 1] = zi[k] + dz[k];
        }

        auto update_wind_and_state = [&](void) {
            wi[0] = ww[0];
            wi[km] = ww[km - 1];
            for (int k = 1; k < km; ++k) {
                wi[k] = (ww[k] * dz[k - 1] + ww[k - 1] * dz[k]) / (dz[k - 1] + dz[k]);
            }

            wi[0] = ww[0];
            wi[1] = Real(0.5) * (ww[1] + ww[0]);
            for (int k = 2; k < km - 1; ++k) {
                wi[k] = Real(9.0) / Real(16.0) * (ww[k] + ww[k - 1])
                      - Real(1.0) / Real(16.0) * (ww[k + 1] + ww[k - 2]);
            }
            if (km > 1) {
                wi[km - 1] = Real(0.5) * (ww[km - 1] + ww[km - 2]);
                wi[km] = ww[km - 1];
            }

            for (int k = 1; k < km; ++k) {
                if (ww[k] == Real(0.0)) wi[k] = ww[k - 1];
            }

            const Real con1 = Real(0.05);
            for (int k = km - 1; k >= 0; --k) {
                const Real decfl = (wi[k + 1] - wi[k]) * dt / dz[k];
                if (decfl > con1) {
                    wi[k] = wi[k + 1] - con1 * dz[k] / dt;
                }
            }

            for (int k = 0; k <= km; ++k) {
                za[k] = zi[k] - wi[k] * dt;
            }

            for (int k = 0; k < km; ++k) {
                dza[k] = za[k + 1] - za[k];
                if (dza[k] <= Real(0.0)) dza[k] = dz[k]; // divergence guard
            }
            dza[km] = zi[km] - za[km]; // Fortran: dza(km+1) = zi(km+1)-za(km+1)
            if (dza[km] <= Real(0.0)) dza[km] = dz[km > 0 ? km - 1 : 0];
            for (int k = 0; k < km; ++k) {
                qa[k] = qq[k] * dz[k] / dza[k];
                qa2[k] = qq2[k] * dz[k] / dza[k];
                qr[k] = qa[k] / den[k];
                qr2[k] = qa2[k] / den[k];
            }
            qa[km] = Real(0.0);
            qa2[km] = Real(0.0);
        };

        update_wind_and_state();

        if (iter > 0) {
            Real n0sfac_dummy;
            for (int k = 0; k < km; ++k) {
                wsm6_slope_snow_cell(qr[k], den[k], denfac[k], tk[k], pidn0s,
                                     Real(0.12), Real(1.0e11), Real(2.0e6),
                                     Real(273.15), Real(WSM6::qcrmin),
                                     rslopesmax, rslopesbmax, rslopes2max,
                                     rslopes3max, WSM6::bvts, pvts,
                                     tmp[k], tmp1[k], tmp2[k], tmp3[k],
                                     wa[k], n0sfac_dummy);
                wsm6_slope_graup_cell(qr2[k], den[k], denfac[k], pidn0g,
                                      Real(WSM6::qcrmin), rslopegmax,
                                      rslopegbmax, rslopeg2max, rslopeg3max,
                                      Real(0.8), pvtg,
                                      tmp[k], tmp1[k], tmp2[k], tmp3[k],
                                      wa2[k]);
            }
            for (int k = 0; k < km; ++k) {
                const Real tmpq = amrex::max(qr[k] + qr2[k], Real(1.0e-15));
                if (tmpq > Real(1.0e-15)) {
                    wa[k] = (wa[k] * qr[k] + wa2[k] * qr2[k]) / tmpq;
                } else {
                    wa[k] = Real(0.0);
                }
            }
            for (int k = 0; k < km; ++k) {
                ww[k] = Real(0.5) * (wd[k] + wa[k]);
                was[k] = wa[k];
            }
            update_wind_and_state();
        }

        for (int ist = 0; ist < 2; ++ist) {
            Real* qn_dst = (ist == 0) ? qn : qn2;
            Real* precip_dst = (ist == 0) ? &precip1[i] : &precip2[i];
            Real* qasrc = (ist == 0) ? qa : qa2;

            for (int k = 1; k < km; ++k) {
                const Real dip = (qasrc[k + 1] - qasrc[k]) / (dza[k + 1] + dza[k]);
                const Real dim = (qasrc[k] - qasrc[k - 1]) / (dza[k - 1] + dza[k]);
                if (dip * dim <= Real(0.0)) {
                    qmi[k] = qasrc[k];
                    qpi[k] = qasrc[k];
                } else {
                    qpi[k] = qasrc[k] + Real(0.5) * (dip + dim) * dza[k];
                    qmi[k] = Real(2.0) * qasrc[k] - qpi[k];
                    if (qpi[k] < Real(0.0) || qmi[k] < Real(0.0)) {
                        qpi[k] = qasrc[k];
                        qmi[k] = qasrc[k];
                    }
                }
            }
            qpi[0] = qasrc[0];
            qmi[0] = qasrc[0];
            qmi[km] = qasrc[km];
            qpi[km] = qasrc[km];

            for (int k = 0; k < km; ++k) {
                qn_dst[k] = Real(0.0);
            }

            int kb = 0;
            int kt = 0;
            for (int k = 0; k < km; ++k) {
                if (zi[k] >= za[km]) {
                    break;
                }

                for (int kk = kb; kk < km; ++kk) {
                    if (zi[k] <= za[kk + 1]) {
                        kb = kk;
                        break;
                    }
                }

                for (int kk = kt; kk < km; ++kk) {
                    if (zi[k + 1] <= za[kk]) {
                        kt = kk;
                        break;
                    }
                }
                kt = amrex::max(kt - 1, 0);

                if (kt == kb) {
                    const Real tl = (zi[k] - za[kb]) / dza[kb];
                    const Real th = (zi[k + 1] - za[kb]) / dza[kb];
                    const Real tl2 = tl * tl;
                    const Real th2 = th * th;
                    const Real qqd = Real(0.5) * (qpi[kb] - qmi[kb]);
                    const Real qqh = qqd * th2 + qmi[kb] * th;
                    const Real qql = qqd * tl2 + qmi[kb] * tl;
                    qn_dst[k] = (qqh - qql) / (th - tl);
                } else if (kt > kb) {
                    const Real tl = (zi[k] - za[kb]) / dza[kb];
                    const Real tl2 = tl * tl;
                    const Real qqd = Real(0.5) * (qpi[kb] - qmi[kb]);
                    const Real qql = qqd * tl2 + qmi[kb] * tl;
                    const Real dql = qasrc[kb] - qql;
                    Real zsum = (Real(1.0) - tl) * dza[kb];
                    Real qsum = dql * dza[kb];
                    if (kt - kb > 1) {
                        for (int m = kb + 1; m < kt; ++m) {
                            zsum += dza[m];
                            qsum += qasrc[m] * dza[m];
                        }
                    }
                    const Real th = (zi[k + 1] - za[kt]) / dza[kt];
                    const Real th2 = th * th;
                    const Real dqh = Real(0.5) * (qpi[kt] - qmi[kt]) * th2 + qmi[kt] * th;
                    zsum += th * dza[kt];
                    qsum += dqh * dza[kt];
                    qn_dst[k] = qsum / zsum;
                }
            }

            Real precip = Real(0.0);
            for (int k = 0; k < km; ++k) {
                if (za[k] < Real(0.0) && za[k + 1] < Real(0.0)) {
                    precip += qasrc[k] * dza[k];
                } else if (za[k] < Real(0.0) && za[k + 1] >= Real(0.0)) {
                    precip += qasrc[k] * (Real(0.0) - za[k]);
                    break;
                } else {
                    break;
                }
            }
            *precip_dst = precip;
        }

        for (int k = 0; k < km; ++k) {
            rql[i * km + k] = qn[k];
            rql2[i * km + k] = qn2[k];
            wwl[i * km + k] = ww[k];
        }
    }
}

void
WSM6::Advance(const Real& dt_advance,
              const SolverChoice&)
{
    dt = dt_advance;

    int microphysics_debug = 0;
    std::string micro_diag_mode = "canonical";
    std::vector<std::string> micro_diag_tags = {"standing"};
    std::vector<std::string> micro_diag_expr = {"standing"};
    std::vector<std::string> micro_diag_store = {"standing"};
    std::vector<int> micro_diag_target_column;
    {
      amrex::ParmParse pp("erf");
      pp.query("microphysics_debug", microphysics_debug);
      pp.query("micro_diag_mode", micro_diag_mode);
      pp.queryarr("micro_diag_tags", micro_diag_tags);
      pp.queryarr("micro_diag_expr", micro_diag_expr);
      pp.queryarr("micro_diag_store", micro_diag_store);
      pp.queryarr("micro_diag_target_column", micro_diag_target_column);
    }
    microphysics_debug = std::max(0, std::min(2, microphysics_debug));
#ifdef ERF_USE_WSM6_FORT
    const std::string micro_diag_mode_lower = [&]{ std::string m = micro_diag_mode; std::transform(m.begin(), m.end(), m.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); }); return m; }();
    const bool micro_diag_forensic = (micro_diag_mode_lower == "forensic" || micro_diag_mode_lower == "both");
    const int microphysics_debug_bridge = micro_diag_forensic
        ? microphysics_debug
        : std::min(microphysics_debug, 1);
    bool use_wsm6_cpp_answer = false;
    { amrex::ParmParse pp("erf");
      pp.query("use_wsm6_cpp_answer", use_wsm6_cpp_answer); }
    bool run_wsm6_fort = !use_wsm6_cpp_answer;

    static bool wsm6_inited = false;

    // Minimal phase-1 initialization for single-moment WSM6.
    if (!wsm6_inited) {
        constexpr double den0 = 1.28;                 // Standard dry-air density (kg/m^3)
        constexpr double denr = static_cast<double>(rhoh2o);
        constexpr double dens = static_cast<double>(rhos);
        constexpr double cl = static_cast<double>(Cp_l);
        constexpr double cpv = static_cast<double>(Cp_v);
        constexpr int hail_opt = 0;                   // Graupel mode
        mp_wsm6_init_c(den0, denr, dens, cl, cpv, hail_opt);
        wsm6_inited = true;
    }
#endif

    constexpr double g = static_cast<double>(CONST_GRAV);
    constexpr double cpd = static_cast<double>(Cp_d);
    constexpr double cpv = static_cast<double>(Cp_v);
    constexpr double rd = static_cast<double>(R_d);
    constexpr double rv = static_cast<double>(R_v);
    constexpr double t0c = 273.15;
    constexpr double ep1 = static_cast<double>(R_v / R_d - one);
    amrex::ignore_unused(g, rd, ep1);
    constexpr double ep2 = static_cast<double>(R_d / R_v);
    constexpr double qmin = 1.0e-12;
    constexpr double xls = static_cast<double>(lsub);
    constexpr double xlv0 = static_cast<double>(lat_vap);
    constexpr double xlf0 = static_cast<double>(lat_ice);
    constexpr double den0 = 1.28;
    constexpr double denr = static_cast<double>(rhoh2o);
    constexpr double cliq = static_cast<double>(Cp_l);
    constexpr double cice = 2106.0;
    constexpr double psat = 610.78;
    for (MFIter mfi(*mic_fab_vars[MicVar_WSM6::qv], TileNoZ()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const Box fab_box = mfi.fabbox();

        auto const& t_arr = mic_fab_vars[MicVar_WSM6::tabs]->array(mfi);
        auto const& qv_arr = mic_fab_vars[MicVar_WSM6::qv]->array(mfi);
        auto const& qc_arr = mic_fab_vars[MicVar_WSM6::qc]->array(mfi);
        auto const& qi_arr = mic_fab_vars[MicVar_WSM6::qi]->array(mfi);
        auto const& qr_arr = mic_fab_vars[MicVar_WSM6::qr]->array(mfi);
        auto const& qs_arr = mic_fab_vars[MicVar_WSM6::qs]->array(mfi);
        auto const& qg_arr = mic_fab_vars[MicVar_WSM6::qg]->array(mfi);
        auto const& den_arr = mic_fab_vars[MicVar_WSM6::rho]->array(mfi);
        auto const& p_arr = mic_fab_vars[MicVar_WSM6::pres]->array(mfi);
        auto rain_arr = mic_fab_vars[MicVar_WSM6::rain_accum]->array(mfi);
        auto snow_arr = mic_fab_vars[MicVar_WSM6::snow_accum]->array(mfi);
        auto graup_arr = mic_fab_vars[MicVar_WSM6::graup_accum]->array(mfi);

        const int ilo = box.smallEnd(0);
        const int ihi = box.bigEnd(0);
        const int jlo = box.smallEnd(1);
        const int jhi = box.bigEnd(1);
        const int klo = box.smallEnd(2);
        const int khi = box.bigEnd(2);
        const bool has_target_override = (micro_diag_target_column.size() == 2);
        const int diag_i = has_target_override ? micro_diag_target_column[0] : ilo;
        const int diag_j = has_target_override ? micro_diag_target_column[1] : jlo;

        const int imlo = fab_box.smallEnd(0);
        const int imhi = fab_box.bigEnd(0);
        const int jmlo = fab_box.smallEnd(1);
        const int jmhi = fab_box.bigEnd(1);
        const int kmlo = fab_box.smallEnd(2);
        const int kmhi = fab_box.bigEnd(2);
        amrex::ignore_unused(ihi, jhi, diag_i, diag_j, imlo, imhi, jmlo, jmhi, kmlo, kmhi);

        const Real dz_val = m_geom.CellSize(2);
        FArrayBox delz_fab(fab_box, 1);
        auto const& delz_arr = delz_fab.array();
        ParallelFor(fab_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            delz_arr(i,j,k) = dz_val;
        });

        const Array4<const Real> z_arr = (m_z_phys_nd) ? m_z_phys_nd->const_array(mfi) : Array4<const Real> {};
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            delz_arr(i,j,k) = (z_arr) ? Real(0.25) * ( (z_arr(i  ,j  ,k+1) - z_arr(i  ,j  ,k))
                                                     + (z_arr(i+1,j  ,k+1) - z_arr(i+1,j  ,k))
                                                     + (z_arr(i  ,j+1,k+1) - z_arr(i  ,j+1,k))
                                                     + (z_arr(i+1,j+1,k+1) - z_arr(i+1,j+1,k)) ) : dz_val;
        });

        Box box2d(box);
        box2d.makeSlab(2, 0);
        Box fab_box2d(fab_box);
        fab_box2d.makeSlab(2, 0);

        // Fortran bridge uses ims:ime, jms:jme storage bounds; these buffers must
        // therefore be allocated on fab_box extents even if C++ kernels only
        // update the valid tile slab (box2d).
        FArrayBox rainncv_fab(fab_box2d, 1);
        FArrayBox sr_fab(fab_box2d, 1);
        FArrayBox snowncv_fab(fab_box2d, 1);
        FArrayBox graupelncv_fab(fab_box2d, 1);
        FArrayBox rainacc_fab(fab_box2d, 1);
        FArrayBox snowacc_fab(fab_box2d, 1);
        FArrayBox graupacc_fab(fab_box2d, 1);

        auto const& rainncv_arr = rainncv_fab.array();
        auto const& sr_arr = sr_fab.array();
        auto const& snowncv_arr = snowncv_fab.array();
        auto const& graupelncv_arr = graupelncv_fab.array();
        auto const& rainacc_arr = rainacc_fab.array();
        auto const& snowacc_arr = snowacc_fab.array();
        auto const& graupacc_arr = graupacc_fab.array();
        ParallelFor(fab_box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            rainncv_arr(i,j,k) = Real(0.0);
            sr_arr(i,j,k) = Real(0.0);
            snowncv_arr(i,j,k) = Real(0.0);
            graupelncv_arr(i,j,k) = Real(0.0);
            rainacc_arr(i,j,k) = Real(0.0);
            snowacc_arr(i,j,k) = Real(0.0);
            graupacc_arr(i,j,k) = Real(0.0);
        });
        ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            rainacc_arr(i,j,0) = rain_arr(i,j,klo);
            snowacc_arr(i,j,0) = snow_arr(i,j,klo);
            graupacc_arr(i,j,0) = graup_arr(i,j,klo);
            rainncv_arr(i,j,0) = Real(0.0);
            sr_arr(i,j,0) = Real(0.0);
            snowncv_arr(i,j,0) = Real(0.0);
            graupelncv_arr(i,j,0) = Real(0.0);
        });

#ifdef ERF_USE_WSM6_FORT
        if (run_wsm6_fort) {
            mp_wsm6_run_c(
                t_arr.dataPtr(),
                qv_arr.dataPtr(), qc_arr.dataPtr(), qi_arr.dataPtr(),
                qr_arr.dataPtr(), qs_arr.dataPtr(), qg_arr.dataPtr(),
                den_arr.dataPtr(), p_arr.dataPtr(), delz_arr.dataPtr(),
                static_cast<double>(dt), g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin,
                xls, xlv0, xlf0, den0, denr, cliq, cice, psat,
                rainacc_arr.dataPtr(), rainncv_arr.dataPtr(), sr_arr.dataPtr(),
                snowacc_arr.dataPtr(), snowncv_arr.dataPtr(),
                graupacc_arr.dataPtr(), graupelncv_arr.dataPtr(),
                imlo, imhi, jmlo, jmhi, kmlo, kmhi,
                ilo, ihi, jlo, jhi, klo, khi, microphysics_debug_bridge, diag_i, diag_j);
        } else {
#endif
        // --- Phase 4 native C++ kernel ---

        // box2d for 1D per-column arrays (already defined above)
        // delqrs1/2/3, delqi: surface precipitation flux accumulators
        FArrayBox delqrs1_fab(box2d,1);
        FArrayBox delqrs2_fab(box2d,1);
        FArrayBox delqrs3_fab(box2d,1);
        FArrayBox delqi_fab(box2d,1);
        FArrayBox tstepsnow_fab(box2d,1);
        FArrayBox tstepgraup_fab(box2d,1);
        auto const& delqrs1_arr    = delqrs1_fab.array();
        auto const& delqrs2_arr    = delqrs2_fab.array();
        auto const& delqrs3_arr    = delqrs3_fab.array();
        auto const& delqi_arr      = delqi_fab.array();
        auto const& tstepsnow_arr   = tstepsnow_fab.array();
        auto const& tstepgraup_arr  = tstepgraup_fab.array();

        ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            delqrs1_arr(i,j,k) = Real(0.0);
            delqrs2_arr(i,j,k) = Real(0.0);
            delqrs3_arr(i,j,k) = Real(0.0);
            delqi_arr(i,j,k) = Real(0.0);
            tstepsnow_arr(i,j,k) = Real(0.0);
            tstepgraup_arr(i,j,k) = Real(0.0);
        });
        // 3D working FABs
        FArrayBox denfac_fab(fab_box,1);  FArrayBox xni_fab(fab_box,1);
        FArrayBox cpm_fab(fab_box,1);     FArrayBox xl_fab(fab_box,1);
        FArrayBox qsatw_fab(fab_box,1);   FArrayBox qsati_fab(fab_box,1);
        FArrayBox rhw_fab(fab_box,1);     FArrayBox rhi_fab(fab_box,1);
        FArrayBox den_tmp_fab(fab_box,1); FArrayBox delz_tmp_fab(fab_box,1);
        FArrayBox n0sfac_fab(fab_box,1);
        FArrayBox qrs_tmp_r_fab(fab_box,1); FArrayBox qrs_tmp_s_fab(fab_box,1);
        FArrayBox qrs_tmp_g_fab(fab_box,1);
        FArrayBox rslope_r_fab(fab_box,1);  FArrayBox rslope_s_fab(fab_box,1);
        FArrayBox rslope_g_fab(fab_box,1);
        FArrayBox rslope2_r_fab(fab_box,1); FArrayBox rslope2_s_fab(fab_box,1);
        FArrayBox rslope2_g_fab(fab_box,1);
        FArrayBox rslope3_r_fab(fab_box,1); FArrayBox rslope3_s_fab(fab_box,1);
        FArrayBox rslope3_g_fab(fab_box,1);
        FArrayBox rslopeb_r_fab(fab_box,1); FArrayBox rslopeb_s_fab(fab_box,1);
        FArrayBox rslopeb_g_fab(fab_box,1);
        FArrayBox work1_r_fab(fab_box,1);   FArrayBox work1_s_fab(fab_box,1);
        FArrayBox work1_g_fab(fab_box,1);
        FArrayBox work2_fab(fab_box,1);     FArrayBox workdiffw_fab(fab_box,1);
        FArrayBox workdiffi_fab(fab_box,1);
        FArrayBox workr_fab(fab_box,1);     FArrayBox worka_fab(fab_box,1);
        FArrayBox work1c_fab(fab_box,1);
        FArrayBox denqrs1_fab(fab_box,1);   FArrayBox denqrs2_fab(fab_box,1);
        FArrayBox denqrs3_fab(fab_box,1);   FArrayBox denqci_fab(fab_box,1);
        FArrayBox fall_r_fab(fab_box,1);    FArrayBox fall_s_fab(fab_box,1);
        FArrayBox fall_g_fab(fab_box,1);    FArrayBox fallc_fab(fab_box,1);
        FArrayBox qsum_fab(fab_box,1);
        FArrayBox nislfv_r_diag_fab(fab_box,6);
        FArrayBox nislfv_sg_diag_fab(fab_box,6);
        // process rates
        FArrayBox praut_fab(fab_box,1); FArrayBox pracw_fab(fab_box,1);
        FArrayBox prevp_fab(fab_box,1); FArrayBox psdep_fab(fab_box,1);
        FArrayBox pgdep_fab(fab_box,1); FArrayBox psaut_fab(fab_box,1);
        FArrayBox pgaut_fab(fab_box,1); FArrayBox praci_fab(fab_box,1);
        FArrayBox piacr_fab(fab_box,1); FArrayBox psaci_fab(fab_box,1);
        FArrayBox psacw_fab(fab_box,1); FArrayBox pgacw_fab(fab_box,1);
        FArrayBox pgaci_fab(fab_box,1); FArrayBox paacw_fab(fab_box,1);
        FArrayBox pracs_fab(fab_box,1); FArrayBox psacr_fab(fab_box,1);
        FArrayBox pgacr_fab(fab_box,1); FArrayBox pgacs_fab(fab_box,1);
        FArrayBox pigen_fab(fab_box,1); FArrayBox pidep_fab(fab_box,1);
        FArrayBox pcond_fab(fab_box,1); FArrayBox psmlt_fab(fab_box,1);
        FArrayBox pgmlt_fab(fab_box,1); FArrayBox pseml_fab(fab_box,1);
        FArrayBox pgeml_fab(fab_box,1); FArrayBox psevp_fab(fab_box,1);
        FArrayBox pgevp_fab(fab_box,1);
        FArrayBox pimlt_fab(fab_box,1); FArrayBox pihmf_fab(fab_box,1);
        FArrayBox pihtf_fab(fab_box,1); FArrayBox pgfrz_fab(fab_box,1);

        auto const& denfac_arr    = denfac_fab.array();
        auto const& xni_arr       = xni_fab.array();
        auto const& cpm_arr       = cpm_fab.array();
        auto const& xl_arr        = xl_fab.array();
        auto const& qsatw_arr     = qsatw_fab.array();
        auto const& qsati_arr     = qsati_fab.array();
        auto const& rhw_arr       = rhw_fab.array();
        auto const& rhi_arr       = rhi_fab.array();
        auto const& den_tmp_arr   = den_tmp_fab.array();
        auto const& delz_tmp_arr  = delz_tmp_fab.array();
        auto const& n0sfac_arr    = n0sfac_fab.array();
        auto const& qrs_tmp_r_arr = qrs_tmp_r_fab.array();
        auto const& qrs_tmp_s_arr = qrs_tmp_s_fab.array();
        auto const& qrs_tmp_g_arr = qrs_tmp_g_fab.array();
        auto const& rslope_r_arr  = rslope_r_fab.array();
        auto const& rslope_s_arr  = rslope_s_fab.array();
        auto const& rslope_g_arr  = rslope_g_fab.array();
        auto const& rslope2_r_arr = rslope2_r_fab.array();
        auto const& rslope2_s_arr = rslope2_s_fab.array();
        auto const& rslope2_g_arr = rslope2_g_fab.array();
        auto const& rslope3_r_arr = rslope3_r_fab.array();
        auto const& rslope3_s_arr = rslope3_s_fab.array();
        auto const& rslope3_g_arr = rslope3_g_fab.array();
        auto const& rslopeb_r_arr = rslopeb_r_fab.array();
        auto const& rslopeb_s_arr = rslopeb_s_fab.array();
        auto const& rslopeb_g_arr = rslopeb_g_fab.array();
        auto const& work1_r_arr   = work1_r_fab.array();
        auto const& work1_s_arr   = work1_s_fab.array();
        auto const& work1_g_arr   = work1_g_fab.array();
        auto const& work2_arr     = work2_fab.array();
        auto const& workdiffw_arr = workdiffw_fab.array();
        auto const& workdiffi_arr = workdiffi_fab.array();
        auto const& workr_arr     = workr_fab.array();
        auto const& worka_arr     = worka_fab.array();
        auto const& denqrs1_arr   = denqrs1_fab.array();
        auto const& denqrs2_arr   = denqrs2_fab.array();
        auto const& denqrs3_arr   = denqrs3_fab.array();
        auto const& fall_r_arr    = fall_r_fab.array();
        auto const& fall_s_arr    = fall_s_fab.array();
        auto const& fall_g_arr    = fall_g_fab.array();
        auto const& fallc_arr     = fallc_fab.array();
        auto const& qsum_arr      = qsum_fab.array();
        auto const& nislfv_r_diag_arr = nislfv_r_diag_fab.array();
        auto const& nislfv_sg_diag_arr = nislfv_sg_diag_fab.array();
        auto const& work1c_arr = work1c_fab.array();

        ParallelFor(fab_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            work1c_arr(i,j,k) = Real(0.0);
        });
        ParallelFor(fab_box, nislfv_r_diag_fab.nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
            nislfv_r_diag_arr(i,j,k,n) = Real(0.0);
            nislfv_sg_diag_arr(i,j,k,n) = Real(0.0);
        });
        auto const& praut_arr     = praut_fab.array();
        auto const& pracw_arr     = pracw_fab.array();
        auto const& prevp_arr     = prevp_fab.array();
        auto const& psdep_arr     = psdep_fab.array();
        auto const& pgdep_arr     = pgdep_fab.array();
        auto const& psaut_arr     = psaut_fab.array();
        auto const& pgaut_arr     = pgaut_fab.array();
        auto const& praci_arr     = praci_fab.array();
        auto const& piacr_arr     = piacr_fab.array();
        auto const& psaci_arr     = psaci_fab.array();
        auto const& psacw_arr     = psacw_fab.array();
        auto const& pgacw_arr     = pgacw_fab.array();
        auto const& pgaci_arr     = pgaci_fab.array();
        auto const& paacw_arr     = paacw_fab.array();
        auto const& pracs_arr     = pracs_fab.array();
        auto const& psacr_arr     = psacr_fab.array();
        auto const& pgacr_arr     = pgacr_fab.array();
        auto const& pgacs_arr     = pgacs_fab.array();
        auto const& pigen_arr     = pigen_fab.array();
        auto const& pidep_arr     = pidep_fab.array();
        auto const& pcond_arr     = pcond_fab.array();
        auto const& psmlt_arr     = psmlt_fab.array();
        auto const& pgmlt_arr     = pgmlt_fab.array();
        auto const& pseml_arr     = pseml_fab.array();
        auto const& pgeml_arr     = pgeml_fab.array();
        auto const& psevp_arr     = psevp_fab.array();
        auto const& pgevp_arr     = pgevp_fab.array();
        auto const& pimlt_arr     = pimlt_fab.array();
        auto const& pihmf_arr     = pihmf_fab.array();
        auto const& pihtf_arr     = pihtf_fab.array();
        auto const& pgfrz_arr     = pgfrz_fab.array();

        // Groups A-E: pre-loop setup
        // Clamp negative values (Group A)
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            qc_arr(i,j,k) = amrex::max(qc_arr(i,j,k), Real(0.0));
            qr_arr(i,j,k) = amrex::max(qr_arr(i,j,k), Real(0.0));
            qi_arr(i,j,k) = amrex::max(qi_arr(i,j,k), Real(0.0));
            qs_arr(i,j,k) = amrex::max(qs_arr(i,j,k), Real(0.0));
            qg_arr(i,j,k) = amrex::max(qg_arr(i,j,k), Real(0.0));
            den_tmp_arr(i,j,k)  = den_arr(i,j,k);
            delz_tmp_arr(i,j,k) = delz_arr(i,j,k);
        });

        // Group B: cpm, xl — computed once from initial state [lines 455-460]
        const Real xlv1_loc = m_xlv1;
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            cpm_arr(i,j,k) = wsm6_cpmcal(qv_arr(i,j,k), Real(qmin), Real(cpd), Real(cpv));
            xl_arr(i,j,k)  = wsm6_xlcal(t_arr(i,j,k), Real(xlv0), xlv1_loc, Real(t0c));
        });

        // Outer minor timestep loop (Rule 29)
        const int wsm6_loops = std::max(
            static_cast<int>(std::round(dt / dtcldcr)), 1);
        const Real dtcld = dt / static_cast<Real>(wsm6_loops);
        const Real qc0    = m_qc0;
        const Real qck1   = m_qck1;
        const Real pvtr   = m_pvtr;
        const Real pacrr  = m_pacrr;
        const Real precr1 = m_precr1;
        const Real precr2 = m_precr2;
        const Real roqimax= m_roqimax;
        const Real pvts   = m_pvts;
        const Real pacrc  = m_pacrc;
        const Real precs1 = m_precs1;
        const Real precs2 = m_precs2;
        const Real g6pbr  = m_g6pbr;
        const Real pvtg   = m_pvtg;
        const Real pacrg  = m_pacrg;
        const Real precg1 = m_precg1;
        const Real precg2 = m_precg2;

        for (int loop = 0; loop < wsm6_loops; ++loop) {
            // G1b: denfac = sqrt(den0/den)  [lines 503-515]
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real invden = Real(1.0) / den_arr(i,j,k);
                denfac_arr(i,j,k) = std::sqrt(invden * Real(den0));
            });
            // G1c: qsatw, qsati, rhw, rhi  [lines 517-549]
            {
                const Real ttp  = Real(t0c) + Real(0.01);
                const Real dldt = Real(cpv) - Real(cliq);
                const Real xa   = -dldt / Real(rv);
                const Real xb   =  xa + Real(xlv0) / (Real(rv)*ttp);
                const Real dldti= Real(cpv) - Real(cice);
                const Real xai  = -dldti / Real(rv);
                const Real xbi  =  xai + Real(xls) / (Real(rv)*ttp);
                ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    const Real tr = ttp / t_arr(i,j,k);
                    Real qsw = Real(psat)*std::exp(std::log(tr)*xa)*std::exp(xb*(Real(1.0)-tr));
                    qsw = amrex::min(qsw, Real(0.99)*p_arr(i,j,k));
                    qsw = Real(ep2)*qsw / (p_arr(i,j,k) - qsw);
                    qsw = amrex::max(qsw, Real(qmin));
                    qsatw_arr(i,j,k) = qsw;
                    rhw_arr(i,j,k)   = amrex::max(qv_arr(i,j,k)/qsw, Real(qmin));
                    Real qsi = (t_arr(i,j,k) < ttp)
                        ? Real(psat)*std::exp(std::log(tr)*xai)*std::exp(xbi*(Real(1.0)-tr))
                        : Real(psat)*std::exp(std::log(tr)*xa )*std::exp(xb *(Real(1.0)-tr));
                    qsi = amrex::min(qsi, Real(0.99)*p_arr(i,j,k));
                    qsi = Real(ep2)*qsi / (p_arr(i,j,k) - qsi);
                    qsi = amrex::max(qsi, Real(qmin));
                    qsati_arr(i,j,k) = qsi;
                    rhi_arr(i,j,k)   = amrex::max(qv_arr(i,j,k)/qsi, Real(qmin));
                });
            }
            // G2: zero all process rates each sub-step  [lines 555-594]
            // WSM6-CPP TAG: RATES_ZERO
            //   legacy_group: G2
            //   process: Zero all process rates each sub-step
            //   compare_vars: prevp, psdep, pgdep, praut, psaut, pgaut, pracw, praci, psaci, pracs, pidep, pcond, psmlt, pgmlt, pseml, psevp, fall_r, fall_s, fall_g, fallc
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                prevp_arr(i,j,k) = Real(0.0); psdep_arr(i,j,k) = Real(0.0);
                pgdep_arr(i,j,k) = Real(0.0); praut_arr(i,j,k) = Real(0.0);
                psaut_arr(i,j,k) = Real(0.0); pgaut_arr(i,j,k) = Real(0.0);
                pracw_arr(i,j,k) = Real(0.0); praci_arr(i,j,k) = Real(0.0);
                piacr_arr(i,j,k) = Real(0.0); psaci_arr(i,j,k) = Real(0.0);
                psacw_arr(i,j,k) = Real(0.0); pracs_arr(i,j,k) = Real(0.0);
                psacr_arr(i,j,k) = Real(0.0); pgacw_arr(i,j,k) = Real(0.0);
                paacw_arr(i,j,k) = Real(0.0); pgaci_arr(i,j,k) = Real(0.0);
                pgacr_arr(i,j,k) = Real(0.0); pgacs_arr(i,j,k) = Real(0.0);
                pigen_arr(i,j,k) = Real(0.0); pidep_arr(i,j,k) = Real(0.0);
                pcond_arr(i,j,k) = Real(0.0); psmlt_arr(i,j,k) = Real(0.0);
                pgmlt_arr(i,j,k) = Real(0.0); pseml_arr(i,j,k) = Real(0.0);
                pgeml_arr(i,j,k) = Real(0.0); psevp_arr(i,j,k) = Real(0.0);
                pgevp_arr(i,j,k) = Real(0.0);
                fall_r_arr(i,j,k) = Real(0.0); fall_s_arr(i,j,k) = Real(0.0);
                fall_g_arr(i,j,k) = Real(0.0); fallc_arr(i,j,k) = Real(0.0);
            });

            // G3: xni ice crystal number concentration  [lines 598-604]
            // WSM6-CPP TAG: XNI
            //   legacy_group: G3
            //   process: Ice crystal number concentration
            //   compare_vars: xni, qi, den
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real tmp = den_arr(i,j,k)*amrex::max(qi_arr(i,j,k), Real(qmin));
                xni_arr(i,j,k) = amrex::min(
                    amrex::max(Real(5.38e7)*std::sqrt(std::sqrt(tmp*tmp*tmp)), Real(1.e3)),
                    Real(1.e6));
            });

            // G4: pack qrs_tmp, first slope_wsm6 [lines 610-618]
            // WSM6-CPP TAG: SLOPE1
            //   legacy_group: G4
            //   process: First slope calculation
            //   compare_vars: rslope, rslope2, rslope3, rslopeb, falk, fall, work1
            const Real pidn0r_loc      = m_pidn0r;
            const Real rslopermax_loc  = m_rslopermax;
            const Real rsloperbmax_loc = m_rsloperbmax;
            const Real rsloper2max_loc = m_rsloper2max;
            const Real rsloper3max_loc = m_rsloper3max;
            const Real pidn0s_loc      = m_pidn0s;
            const Real rslopesmax_loc  = m_rslopesmax;
            const Real rslopesbmax_loc = m_rslopesbmax;
            const Real rslopes2max_loc = m_rslopes2max;
            const Real rslopes3max_loc = m_rslopes3max;
            const Real pidn0g_loc      = m_pidn0g;
            const Real rslopegmax_loc  = m_rslopegmax;
            const Real rslopegbmax_loc = m_rslopegbmax;
            const Real rslopeg2max_loc = m_rslopeg2max;
            const Real rslopeg3max_loc = m_rslopeg3max;
            const Real bvtg_loc        = m_bvtg;
            const Real n0g_loc         = m_n0g;

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                qrs_tmp_r_arr(i,j,k) = qr_arr(i,j,k);
                qrs_tmp_s_arr(i,j,k) = qs_arr(i,j,k);
                qrs_tmp_g_arr(i,j,k) = qg_arr(i,j,k);
                Real dummy_n0sfac;
                wsm6_slope_rain_cell(
                    qrs_tmp_r_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k),
                    pidn0r_loc, Real(qcrmin), rslopermax_loc, rsloperbmax_loc,
                    rsloper2max_loc, rsloper3max_loc, Real(bvtr), Real(pvtr),
                    rslope_r_arr(i,j,k), rslopeb_r_arr(i,j,k),
                    rslope2_r_arr(i,j,k), rslope3_r_arr(i,j,k),
                    work1_r_arr(i,j,k));
                wsm6_slope_snow_cell(
                    qrs_tmp_s_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k),
                    t_arr(i,j,k), pidn0s_loc, Real(alpha_wsm6),
                    Real(n0smax), Real(n0s), Real(t0c), Real(qcrmin),
                    rslopesmax_loc, rslopesbmax_loc,
                    rslopes2max_loc, rslopes3max_loc,
                    Real(bvts), Real(pvts),
                    rslope_s_arr(i,j,k), rslopeb_s_arr(i,j,k),
                    rslope2_s_arr(i,j,k), rslope3_s_arr(i,j,k),
                    work1_s_arr(i,j,k), dummy_n0sfac);
            wsm6_slope_graup_cell(
                    qrs_tmp_g_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k),
                    pidn0g_loc, Real(qcrmin),
                    rslopegmax_loc, rslopegbmax_loc,
                    rslopeg2max_loc, rslopeg3max_loc,
                    bvtg_loc, Real(pvtg),
                    rslope_g_arr(i,j,k), rslopeb_g_arr(i,j,k),
                    rslope2_g_arr(i,j,k), rslope3_g_arr(i,j,k),
                    work1_g_arr(i,j,k));
                n0sfac_arr(i,j,k) = dummy_n0sfac;
            });
            // G5a-G5e: sedimentation setup, nislfv calls, and flux updates
            ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
                const int km_local = khi - klo + 1;
                if (km_local > WSM6_MAX_LEVELS) return;

                constexpr Real qsum_min = Real(1.0e-15);
                Real den_col[WSM6_MAX_LEVELS];
                Real denfac_col[WSM6_MAX_LEVELS];
                Real t_col[WSM6_MAX_LEVELS];
                Real dz_col[WSM6_MAX_LEVELS];
                Real workr_col[WSM6_MAX_LEVELS];
                Real worka_col[WSM6_MAX_LEVELS];
                Real denqrs1_col[WSM6_MAX_LEVELS];
                Real denqrs2_col[WSM6_MAX_LEVELS];
                Real denqrs3_col[WSM6_MAX_LEVELS];
                Real qsum_col[WSM6_MAX_LEVELS];
                Real delqrs1_col = Real(0.0);
                Real delqrs2_col = Real(0.0);
                Real delqrs3_col = Real(0.0);

                // G5a: pack sedimentation work arrays
                for (int k = klo; k <= khi; ++k) {
                    const int kk = k - klo;
                    den_col[kk]    = den_arr(i,j,k);
                    denfac_col[kk] = denfac_arr(i,j,k);
                    t_col[kk]      = t_arr(i,j,k);
                    dz_col[kk]     = delz_tmp_arr(i,j,k);
                    workr_col[kk]  = work1_r_arr(i,j,k);
                    qsum_col[kk]   = amrex::max(qs_arr(i,j,k) + qg_arr(i,j,k), qsum_min);
                    if (qsum_col[kk] > qsum_min) {
                        worka_col[kk] = (work1_s_arr(i,j,k) * qs_arr(i,j,k)
                                       + work1_g_arr(i,j,k) * qg_arr(i,j,k))
                                      / qsum_col[kk];
                    } else {
                        worka_col[kk] = Real(0.0);
                    }
                    denqrs1_col[kk] = den_col[kk] * qr_arr(i,j,k);
                    denqrs2_col[kk] = den_col[kk] * qs_arr(i,j,k);
                    denqrs3_col[kk] = den_col[kk] * qg_arr(i,j,k);
                    if (qr_arr(i,j,k) <= Real(0.0)) {
                        workr_col[kk] = Real(0.0);
                    }
                }

                // G5b: rain sedimentation
                wsm6_nislfv_rain_plm(
                    1, km_local, den_col, denfac_col, t_col, dz_col,
                    workr_col, denqrs1_col, &delqrs1_col, dtcld, 1, 1, 0,
                    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                    nullptr, nullptr, nullptr, nullptr);
                // Strict Rule 30 snapshot: immediately after G5b
                for (int k = klo; k <= khi; ++k) {
                    const int kk = k - klo;
                    nislfv_r_diag_arr(i,j,k,0) = amrex::max(denqrs1_col[kk] / den_col[kk], Real(0.0));
                    nislfv_r_diag_arr(i,j,k,1) = denqrs1_col[kk] * workr_col[kk] / delz_arr(i,j,k);
                    nislfv_r_diag_arr(i,j,k,2) = workr_col[kk];
                    nislfv_r_diag_arr(i,j,k,3) = denqrs1_col[kk];
                    nislfv_r_diag_arr(i,j,k,4) = den_col[kk];
                    nislfv_r_diag_arr(i,j,k,5) = denfac_col[kk];
                }

                // G5c: snow + graupel sedimentation
                wsm6_nislfv_rain_plm6(
                    1, km_local, den_col, denfac_col, t_col, dz_col,
                    worka_col, denqrs2_col, denqrs3_col,
                    &delqrs2_col, &delqrs3_col, dtcld, 1, 1);
                // Strict Rule 30 snapshot: immediately after G5c
                for (int k = klo; k <= khi; ++k) {
                    const int kk = k - klo;
                    nislfv_sg_diag_arr(i,j,k,0) = amrex::max(denqrs2_col[kk] / den_col[kk], Real(0.0));
                    nislfv_sg_diag_arr(i,j,k,1) = amrex::max(denqrs3_col[kk] / den_col[kk], Real(0.0));
                    nislfv_sg_diag_arr(i,j,k,2) = denqrs2_col[kk] * worka_col[kk] / delz_arr(i,j,k);
                    nislfv_sg_diag_arr(i,j,k,3) = denqrs3_col[kk] * worka_col[kk] / delz_arr(i,j,k);
                    nislfv_sg_diag_arr(i,j,k,4) = denqrs2_col[kk];
                    nislfv_sg_diag_arr(i,j,k,5) = denqrs3_col[kk];
                }

                // G5d: update species and fall speeds
                for (int k = klo; k <= khi; ++k) {
                    const int kk = k - klo;
                    qsum_arr(i,j,k) = qsum_col[kk];
                    workr_arr(i,j,k) = workr_col[kk];
                    worka_arr(i,j,k) = worka_col[kk];
                    denqrs1_arr(i,j,k) = denqrs1_col[kk];
                    denqrs2_arr(i,j,k) = denqrs2_col[kk];
                    denqrs3_arr(i,j,k) = denqrs3_col[kk];
                    qr_arr(i,j,k) = amrex::max(denqrs1_col[kk] / den_col[kk], Real(0.0));
                    qs_arr(i,j,k) = amrex::max(denqrs2_col[kk] / den_col[kk], Real(0.0));
                    qg_arr(i,j,k) = amrex::max(denqrs3_col[kk] / den_col[kk], Real(0.0));
                    fall_r_arr(i,j,k) = denqrs1_col[kk] * workr_col[kk] / delz_arr(i,j,k);
                    fall_s_arr(i,j,k) = denqrs2_col[kk] * worka_col[kk] / delz_arr(i,j,k);
                    fall_g_arr(i,j,k) = denqrs3_col[kk] * worka_col[kk] / delz_arr(i,j,k);
                }

                // G5e: slab fall fluxes at the lower boundary
                delqrs1_arr(i,j,0) = delqrs1_col / delz_arr(i,j,klo) / dtcld;
                delqrs2_arr(i,j,0) = delqrs2_col / delz_arr(i,j,klo) / dtcld;
                delqrs3_arr(i,j,0) = delqrs3_col / delz_arr(i,j,klo) / dtcld;
                fall_r_arr(i,j,klo) = delqrs1_arr(i,j,0);
                fall_s_arr(i,j,klo) = delqrs2_arr(i,j,0);
                fall_g_arr(i,j,klo) = delqrs3_arr(i,j,0);
            });
            constexpr Real pi = Real(3.141592653589793238462643383279502884);

            // G10: instantaneous phase changes [lines 774-830]
            // WSM6-CPP TAG: PHASE
            //   legacy_group: G10
            //   process: Instantaneous phase changes
            //   compare_vars: t, q, qci, qrs
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supcol = Real(t0c) - t_arr(i,j,k);
                const Real xlf = (supcol < Real(0.0))
                    ? Real(xlf0)
                    : Real(xls) - xl_arr(i,j,k);
                pimlt_arr(i,j,k) = Real(0.0);
                pihmf_arr(i,j,k) = Real(0.0);
                pihtf_arr(i,j,k) = Real(0.0);
                pgfrz_arr(i,j,k) = Real(0.0);

                if (supcol < Real(0.0) && qi_arr(i,j,k) > Real(0.0)) {
                    pimlt_arr(i,j,k) = qi_arr(i,j,k);
                    qc_arr(i,j,k) = qc_arr(i,j,k) + qi_arr(i,j,k);
                    t_arr(i,j,k) = t_arr(i,j,k) - xlf / cpm_arr(i,j,k) * qi_arr(i,j,k);
                    qi_arr(i,j,k) = Real(0.0);
                }

                if (supcol > Real(40.0) && qc_arr(i,j,k) > Real(0.0)) {
                    pihmf_arr(i,j,k) = qc_arr(i,j,k);
                    qi_arr(i,j,k) = qi_arr(i,j,k) + qc_arr(i,j,k);
                    t_arr(i,j,k) = t_arr(i,j,k) + xlf / cpm_arr(i,j,k) * qc_arr(i,j,k);
                    qc_arr(i,j,k) = Real(0.0);
                }

                if (supcol > Real(0.0) && qc_arr(i,j,k) > Real(qmin)) {
                    const Real supcolt = amrex::min(supcol, Real(50.0));
                    const Real pfrzdtc = amrex::min(
                        Real(pfrz1) *
                        (std::exp(Real(pfrz2) * supcolt) - Real(1.0)) *
                        den_arr(i,j,k) / Real(denr) / Real(WSM6::xncr) *
                        qc_arr(i,j,k) * qc_arr(i,j,k) * dtcld,
                        qc_arr(i,j,k));
                    pihtf_arr(i,j,k) = pfrzdtc;
                    qi_arr(i,j,k) = qi_arr(i,j,k) + pfrzdtc;
                    t_arr(i,j,k) = t_arr(i,j,k) + xlf / cpm_arr(i,j,k) * pfrzdtc;
                    qc_arr(i,j,k) = qc_arr(i,j,k) - pfrzdtc;
                }

                if (supcol > Real(0.0) && qr_arr(i,j,k) > Real(0.0)) {
                    Real temp = rslope3_r_arr(i,j,k);
                    temp = temp * temp * rslope_r_arr(i,j,k);
                    const Real supcolt = amrex::min(supcol, Real(50.0));
                    const Real pfrzdtr = amrex::min(
                        Real(20.0) * pi * pi * Real(pfrz1) * Real(WSM6::n0r) *
                        Real(denr) / den_arr(i,j,k) *
                        (std::exp(Real(pfrz2) * supcolt) - Real(1.0)) *
                        temp * dtcld,
                        qr_arr(i,j,k));
                    pgfrz_arr(i,j,k) = pfrzdtr;
                    qg_arr(i,j,k) = qg_arr(i,j,k) + pfrzdtr;
                    t_arr(i,j,k) = t_arr(i,j,k) + xlf / cpm_arr(i,j,k) * pfrzdtr;
                    qr_arr(i,j,k) = qr_arr(i,j,k) - pfrzdtr;
                }
            });
            // G11: third slope_wsm6 call [lines 836-844]
            // WSM6-CPP TAG: SLOPE3
            //   legacy_group: G11
            //   process: Third slope calculation
            //   compare_vars: rslope, rslope2, rslope3, rslopeb, falk, fall, work1
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                qrs_tmp_r_arr(i,j,k) = qr_arr(i,j,k);
                qrs_tmp_s_arr(i,j,k) = qs_arr(i,j,k);
                qrs_tmp_g_arr(i,j,k) = qg_arr(i,j,k);
                Real dummy_n0sfac;
                wsm6_slope_rain_cell(
                    qrs_tmp_r_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k),
                    pidn0r_loc, Real(qcrmin), rslopermax_loc, rsloperbmax_loc,
                    rsloper2max_loc, rsloper3max_loc, Real(bvtr), Real(pvtr),
                    rslope_r_arr(i,j,k), rslopeb_r_arr(i,j,k),
                    rslope2_r_arr(i,j,k), rslope3_r_arr(i,j,k),
                    work1_r_arr(i,j,k));
                wsm6_slope_snow_cell(
                    qrs_tmp_s_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k),
                    t_arr(i,j,k), pidn0s_loc, Real(alpha_wsm6),
                    Real(n0smax), Real(n0s), Real(t0c), Real(qcrmin),
                    rslopesmax_loc, rslopesbmax_loc,
                    rslopes2max_loc, rslopes3max_loc,
                    Real(bvts), Real(pvts),
                    rslope_s_arr(i,j,k), rslopeb_s_arr(i,j,k),
                    rslope2_s_arr(i,j,k), rslope3_s_arr(i,j,k),
                    work1_s_arr(i,j,k), dummy_n0sfac);
                wsm6_slope_graup_cell(
                    qrs_tmp_g_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k),
                    pidn0g_loc, Real(qcrmin),
                    rslopegmax_loc, rslopegbmax_loc,
                    rslopeg2max_loc, rslopeg3max_loc,
                    bvtg_loc, Real(pvtg),
                    rslope_g_arr(i,j,k), rslopeb_g_arr(i,j,k),
                    rslope2_g_arr(i,j,k), rslope3_g_arr(i,j,k),
                    work1_g_arr(i,j,k));
                n0sfac_arr(i,j,k) = dummy_n0sfac;
            });
            // G12: workdiffw, workdiffi, work2 [lines 851-857]
            // WSM6-CPP TAG: DIFF_PREP
            //   legacy_group: G12
            //   process: Prepare diffusion/work terms
            //   compare_vars: workdiffw, workdiffi, work2
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                workdiffw_arr(i,j,k) = wsm6_diffac(
                    xl_arr(i,j,k), p_arr(i,j,k), t_arr(i,j,k),
                    den_arr(i,j,k), qsatw_arr(i,j,k), Real(rv));
                workdiffi_arr(i,j,k) = wsm6_diffac(
                    Real(xls), p_arr(i,j,k), t_arr(i,j,k),
                    den_arr(i,j,k), qsati_arr(i,j,k), Real(rv));
                work2_arr(i,j,k) = wsm6_venfac(
                    p_arr(i,j,k), t_arr(i,j,k), den_arr(i,j,k), Real(den0));
            });

            // G13a: warm rain — praut, pracw, prevp [lines 867-903]
            {
                const Real dtcld_l = dtcld;
                const Real qmin_l  = Real(qmin);
                const Real qcrmin_l= Real(qcrmin);
                ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Real supsat = amrex::max(qv_arr(i,j,k), qmin_l)
                                - qsatw_arr(i,j,k);
                    Real satdt  = supsat / dtcld_l;
                    // praut: C->R autoconversion
                    if (qc_arr(i,j,k) > Real(qc0)) {
                        praut_arr(i,j,k) = amrex::min(
                            Real(qck1)*std::pow(qc_arr(i,j,k), Real(7.0/3.0)),
                            qc_arr(i,j,k)/dtcld_l);
                    }
                    // pracw: C->R accretion by rain
                    if (qr_arr(i,j,k) > qcrmin_l && qc_arr(i,j,k) > qmin_l) {
                        pracw_arr(i,j,k) = amrex::min(
                            Real(pacrr)*rslope3_r_arr(i,j,k)
                            *rslopeb_r_arr(i,j,k)
                            *qc_arr(i,j,k)*denfac_arr(i,j,k),
                            qc_arr(i,j,k)/dtcld_l);
                    }
                    // prevp: R evaporation/condensation
                    if (qr_arr(i,j,k) > Real(0.0)) {
                        Real coeres = rslope2_r_arr(i,j,k)*std::sqrt(
                            rslope_r_arr(i,j,k)*rslopeb_r_arr(i,j,k));
                        Real rate = (rhw_arr(i,j,k)-Real(1.0))
                            *(Real(precr1)*rslope2_r_arr(i,j,k)
                            + Real(precr2)*work2_arr(i,j,k)*coeres)
                            / workdiffw_arr(i,j,k);
                        if (rate < Real(0.0)) {
                            rate = amrex::max(rate, -qr_arr(i,j,k)/dtcld_l);
                            rate = amrex::max(rate, satdt/Real(2.0));
                        } else {
                            rate = amrex::min(rate, satdt/Real(2.0));
                        }
                        prevp_arr(i,j,k) = rate;
                    }
                });
            }
            // G13b: cold-rain, mixed-phase, and ice deposition/nucleation
            // [lines 904-1062, 1075-1192]
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supcol = Real(t0c) - t_arr(i,j,k);
                n0sfac_arr(i,j,k) = amrex::max(
                    amrex::min(std::exp(Real(alpha_wsm6) * supcol),
                               Real(n0smax) / Real(n0s)),
                    Real(1.0));

                const Real supsat = amrex::max(qv_arr(i,j,k), Real(qmin))
                                  - qsati_arr(i,j,k);
                const Real satdt = supsat / dtcld;
                int ifsat = 0;

                const Real tmp = den_arr(i,j,k)
                               * amrex::max(qi_arr(i,j,k), Real(qmin));
                const Real temp = std::sqrt(std::sqrt(tmp * tmp * tmp));
                xni_arr(i,j,k) = amrex::min(
                    amrex::max(Real(5.38e7) * temp, Real(1.0e3)),
                    Real(1.0e6));

                const Real xni = xni_arr(i,j,k);
                const Real eacrs = std::exp(Real(0.07) * (-supcol));

                const Real xmi = den_arr(i,j,k) * qi_arr(i,j,k) / xni;
                const Real diameter = amrex::min(
                    Real(dicon) * std::sqrt(xmi), Real(dimax));
                const Real vt2i = Real(1.49e4) * std::pow(diameter, Real(1.31));
                const Real vt2r = Real(pvtr) * rslopeb_r_arr(i,j,k)
                                * denfac_arr(i,j,k);
                const Real vt2s = Real(pvts) * rslopeb_s_arr(i,j,k)
                                * denfac_arr(i,j,k);
                const Real vt2g = Real(pvtg) * rslopeb_g_arr(i,j,k)
                                * denfac_arr(i,j,k);

                const Real qsum = amrex::max(qsum_arr(i,j,k), Real(1.0e-15));
                const Real vt2ave = (qsum > Real(1.0e-15))
                    ? (vt2s * qs_arr(i,j,k) + vt2g * qg_arr(i,j,k)) / qsum
                    : Real(0.0);

                if (supcol > Real(0.0) && qi_arr(i,j,k) > Real(qmin)) {
                    if (qr_arr(i,j,k) > Real(qcrmin)) {
                        const Real acrfac =
                            Real(2.0) * rslope3_r_arr(i,j,k)
                          + Real(2.0) * diameter * rslope2_r_arr(i,j,k)
                          + diameter * diameter * rslope_r_arr(i,j,k);

                        // WSM6-CPP TAG: PRACI
                        //   legacy_group: G13b
                        //   process: Accretion of cloud ice by rain
                        //   compare_vars: praci, qi, qr, den
                        praci_arr(i,j,k) = Real(pi) * qi_arr(i,j,k) * Real(n0r)
                            * std::abs(vt2r - vt2i) * acrfac / Real(4.0);
                        praci_arr(i,j,k) *= std::pow(
                            amrex::min(
                                amrex::max(Real(0.0), qr_arr(i,j,k) / qi_arr(i,j,k)),
                                Real(1.0)),
                            Real(2.0));
                        praci_arr(i,j,k) = amrex::min(
                            praci_arr(i,j,k), qi_arr(i,j,k) / dtcld);

                        piacr_arr(i,j,k) = Real(pi) * Real(pi) * Real(avtr)
                            * Real(n0r) * Real(denr) * xni * denfac_arr(i,j,k)
                            * Real(g6pbr) * rslope3_r_arr(i,j,k)
                            * rslope3_r_arr(i,j,k) * rslopeb_r_arr(i,j,k)
                            / Real(24.0) / den_arr(i,j,k);
                        piacr_arr(i,j,k) *= std::pow(
                            amrex::min(
                                amrex::max(Real(0.0), qi_arr(i,j,k) / qr_arr(i,j,k)),
                                Real(1.0)),
                            Real(2.0));
                        piacr_arr(i,j,k) = amrex::min(
                            piacr_arr(i,j,k), qr_arr(i,j,k) / dtcld);
                    }

                    if (qs_arr(i,j,k) > Real(qcrmin)) {
                        const Real acrfac =
                            Real(2.0) * rslope3_s_arr(i,j,k)
                          + Real(2.0) * diameter * rslope2_s_arr(i,j,k)
                          + diameter * diameter * rslope_s_arr(i,j,k);
                        // WSM6-CPP TAG: PSACI
                        //   legacy_group: G13e
                        //   process: Accretion of cloud ice by snow
                        //   compare_vars: psaci, qi, qs, den
                        psaci_arr(i,j,k) = Real(pi) * qi_arr(i,j,k) * eacrs
                            * Real(n0s) * n0sfac_arr(i,j,k)
                            * std::abs(vt2ave - vt2i) * acrfac / Real(4.0);
                        psaci_arr(i,j,k) = amrex::min(
                            psaci_arr(i,j,k), qi_arr(i,j,k) / dtcld);
                    }

                    if (qg_arr(i,j,k) > Real(qcrmin)) {
                        const Real egi = std::exp(Real(0.07) * (-supcol));
                        const Real acrfac =
                            Real(2.0) * rslope3_g_arr(i,j,k)
                          + Real(2.0) * diameter * rslope2_g_arr(i,j,k)
                          + diameter * diameter * rslope_g_arr(i,j,k);
                        pgaci_arr(i,j,k) = Real(pi) * egi * qi_arr(i,j,k)
                              * n0g_loc * std::abs(vt2ave - vt2i) * acrfac
                             / Real(4.0);
                        pgaci_arr(i,j,k) = amrex::min(
                            pgaci_arr(i,j,k), qi_arr(i,j,k) / dtcld);
                    }
                }

                if (qs_arr(i,j,k) > Real(qcrmin) && qc_arr(i,j,k) > Real(qmin)) {
                    psacw_arr(i,j,k) = amrex::min(
                        Real(pacrc) * n0sfac_arr(i,j,k) * rslope3_s_arr(i,j,k)
                        * rslopeb_s_arr(i,j,k)
                        * std::pow(
                            amrex::min(
                                amrex::max(Real(0.0), qs_arr(i,j,k) / qc_arr(i,j,k)),
                                Real(1.0)),
                            Real(2.0))
                        * qc_arr(i,j,k) * denfac_arr(i,j,k),
                        qc_arr(i,j,k) / dtcld);
                }

                if (qg_arr(i,j,k) > Real(qcrmin) && qc_arr(i,j,k) > Real(qmin)) {
                    pgacw_arr(i,j,k) = amrex::min(
                        Real(pacrg) * rslope3_g_arr(i,j,k) * rslopeb_g_arr(i,j,k)
                        * std::pow(
                            amrex::min(
                                amrex::max(Real(0.0), qg_arr(i,j,k) / qc_arr(i,j,k)),
                                Real(1.0)),
                            Real(2.0))
                        * qc_arr(i,j,k) * denfac_arr(i,j,k),
                        qc_arr(i,j,k) / dtcld);
                }

                if (qsum > Real(1.0e-15)) {
                    paacw_arr(i,j,k) = (qs_arr(i,j,k) * psacw_arr(i,j,k)
                                      + qg_arr(i,j,k) * pgacw_arr(i,j,k))
                                     / qsum;
                }

                if (qs_arr(i,j,k) > Real(qcrmin) && qr_arr(i,j,k) > Real(qcrmin)) {
                    if (supcol > Real(0.0)) {
                        const Real acrfac =
                            Real(5.0) * rslope3_s_arr(i,j,k) * rslope3_s_arr(i,j,k)
                                * rslope_r_arr(i,j,k)
                          + Real(2.0) * rslope3_s_arr(i,j,k) * rslope2_s_arr(i,j,k)
                                * rslope2_r_arr(i,j,k)
                          + Real(0.5) * rslope2_s_arr(i,j,k) * rslope2_s_arr(i,j,k)
                                * rslope3_r_arr(i,j,k);
                        // WSM6-CPP TAG: PRACS
                        //   legacy_group: G13f
                        //   process: Accretion of snow by rain / rain-snow interaction
                        //   compare_vars: pracs, qr, qs, den
                        pracs_arr(i,j,k) = Real(pi) * Real(pi) * Real(n0r)
                            * Real(n0s) * n0sfac_arr(i,j,k)
                            * std::abs(vt2r - vt2ave) * (Real(dens_snow) / den_arr(i,j,k))
                            * acrfac;
                        pracs_arr(i,j,k) *= std::pow(
                            amrex::min(
                                amrex::max(Real(0.0), qr_arr(i,j,k) / qs_arr(i,j,k)),
                                Real(1.0)),
                            Real(2.0));
                        pracs_arr(i,j,k) = amrex::min(
                            pracs_arr(i,j,k), qs_arr(i,j,k) / dtcld);
                    }

                    {
                        const Real acrfac =
                            Real(5.0) * rslope3_r_arr(i,j,k) * rslope3_r_arr(i,j,k)
                                * rslope_s_arr(i,j,k)
                          + Real(2.0) * rslope3_r_arr(i,j,k) * rslope2_r_arr(i,j,k)
                                * rslope2_s_arr(i,j,k)
                          + Real(0.5) * rslope2_r_arr(i,j,k) * rslope2_r_arr(i,j,k)
                                * rslope3_s_arr(i,j,k);
                        psacr_arr(i,j,k) = Real(pi) * Real(pi) * Real(n0r)
                            * Real(n0s) * n0sfac_arr(i,j,k)
                            * std::abs(vt2ave - vt2r) * (Real(denr) / den_arr(i,j,k))
                            * acrfac;
                        psacr_arr(i,j,k) *= std::pow(
                            amrex::min(
                                amrex::max(Real(0.0), qs_arr(i,j,k) / qr_arr(i,j,k)),
                                Real(1.0)),
                            Real(2.0));
                        psacr_arr(i,j,k) = amrex::min(
                            psacr_arr(i,j,k), qr_arr(i,j,k) / dtcld);
                    }

                    if (qg_arr(i,j,k) > Real(qcrmin)) {
                        const Real acrfac =
                            Real(5.0) * rslope3_r_arr(i,j,k) * rslope3_r_arr(i,j,k)
                                * rslope_g_arr(i,j,k)
                          + Real(2.0) * rslope3_r_arr(i,j,k) * rslope2_r_arr(i,j,k)
                                * rslope2_g_arr(i,j,k)
                          + Real(0.5) * rslope2_r_arr(i,j,k) * rslope2_r_arr(i,j,k)
                                * rslope3_g_arr(i,j,k);
                        pgacr_arr(i,j,k) = Real(pi) * Real(pi) * Real(n0r)
                            * n0g_loc * std::abs(vt2ave - vt2r)
                            * (Real(denr) / den_arr(i,j,k)) * acrfac;
                        pgacr_arr(i,j,k) *= std::pow(
                            amrex::min(
                                amrex::max(Real(0.0), qg_arr(i,j,k) / qr_arr(i,j,k)),
                                Real(1.0)),
                            Real(2.0));
                        pgacr_arr(i,j,k) = amrex::min(
                            pgacr_arr(i,j,k), qr_arr(i,j,k) / dtcld);
                    }
                }

                if (qg_arr(i,j,k) > Real(qcrmin) && qs_arr(i,j,k) > Real(qcrmin)) {
                    pgacs_arr(i,j,k) = Real(0.0);
                }

                if (supcol <= Real(0.0)) {
                    const Real xlf = Real(xlf0);
                    if (qs_arr(i,j,k) > Real(0.0)) {
                        // WSM6-CPP TAG: PSEML
                        //   legacy_group: G13g
                        //   process: Snow evaporation/sublimation
                        //   compare_vars: pseml, qs, qv, qsat, den
                        pseml_arr(i,j,k) = amrex::min(
                            amrex::max(
                                Real(cliq) * supcol
                                * (paacw_arr(i,j,k) + psacr_arr(i,j,k)) / xlf,
                                -qs_arr(i,j,k) / dtcld),
                            Real(0.0));
                    }
                    if (qg_arr(i,j,k) > Real(0.0)) {
                        pgeml_arr(i,j,k) = amrex::min(
                            amrex::max(
                                Real(cliq) * supcol
                                * (paacw_arr(i,j,k) + pgacr_arr(i,j,k)) / xlf,
                                -qg_arr(i,j,k) / dtcld),
                            Real(0.0));
                    }
                }

                if (supcol > Real(0.0)) {
                    if (qi_arr(i,j,k) > Real(0.0) && ifsat != 1) {
                        // WSM6-CPP TAG: PIDEP
                        //   legacy_group: G13h
                        //   process: Ice deposition/sublimation
                        //   compare_vars: pidep, qi, qv, qsat, den, t
                        pidep_arr(i,j,k) = Real(4.0) * diameter * xni
                            * (rhi_arr(i,j,k) - Real(1.0))
                            / workdiffi_arr(i,j,k);
                        Real supice = satdt - prevp_arr(i,j,k);
                        if (pidep_arr(i,j,k) < Real(0.0)) {
                            pidep_arr(i,j,k) = amrex::max(
                                amrex::max(pidep_arr(i,j,k), satdt / Real(2.0)),
                                supice);
                            pidep_arr(i,j,k) = amrex::max(
                                pidep_arr(i,j,k), -qi_arr(i,j,k) / dtcld);
                        } else {
                            pidep_arr(i,j,k) = amrex::min(
                                amrex::min(pidep_arr(i,j,k), satdt / Real(2.0)),
                                supice);
                        }
                        if (std::abs(prevp_arr(i,j,k) + pidep_arr(i,j,k))
                            >= std::abs(satdt)) {
                            ifsat = 1;
                        }
                    }

                    if (qs_arr(i,j,k) > Real(0.0) && ifsat != 1) {
                        const Real coeres = rslope2_s_arr(i,j,k)
                                          * std::sqrt(rslope_s_arr(i,j,k)
                                                      * rslopeb_s_arr(i,j,k));
                        psdep_arr(i,j,k) = (rhi_arr(i,j,k) - Real(1.0))
                            * n0sfac_arr(i,j,k)
                            * (Real(precs1) * rslope2_s_arr(i,j,k)
                               + Real(precs2) * work2_arr(i,j,k) * coeres)
                            / workdiffi_arr(i,j,k);
                        Real supice = satdt - prevp_arr(i,j,k) - pidep_arr(i,j,k);
                        if (psdep_arr(i,j,k) < Real(0.0)) {
                            psdep_arr(i,j,k) = amrex::max(
                                psdep_arr(i,j,k), -qs_arr(i,j,k) / dtcld);
                            psdep_arr(i,j,k) = amrex::max(
                                amrex::max(psdep_arr(i,j,k), satdt / Real(2.0)),
                                supice);
                        } else {
                            psdep_arr(i,j,k) = amrex::min(
                                amrex::min(psdep_arr(i,j,k), satdt / Real(2.0)),
                                supice);
                        }
                        if (std::abs(prevp_arr(i,j,k) + pidep_arr(i,j,k)
                                     + psdep_arr(i,j,k))
                            >= std::abs(satdt)) {
                            ifsat = 1;
                        }
                    }

                    if (qg_arr(i,j,k) > Real(0.0) && ifsat != 1) {
                        const Real coeres = rslope2_g_arr(i,j,k)
                                          * std::sqrt(rslope_g_arr(i,j,k)
                                                      * rslopeb_g_arr(i,j,k));
                        pgdep_arr(i,j,k) = (rhi_arr(i,j,k) - Real(1.0))
                            * (Real(precg1) * rslope2_g_arr(i,j,k)
                               + Real(precg2) * work2_arr(i,j,k) * coeres)
                            / workdiffi_arr(i,j,k);
                        Real supice = satdt - prevp_arr(i,j,k)
                                    - pidep_arr(i,j,k) - psdep_arr(i,j,k);
                        if (pgdep_arr(i,j,k) < Real(0.0)) {
                            pgdep_arr(i,j,k) = amrex::max(
                                pgdep_arr(i,j,k), -qg_arr(i,j,k) / dtcld);
                            pgdep_arr(i,j,k) = amrex::max(
                                amrex::max(pgdep_arr(i,j,k), satdt / Real(2.0)),
                                supice);
                        } else {
                            pgdep_arr(i,j,k) = amrex::min(
                                amrex::min(pgdep_arr(i,j,k), satdt / Real(2.0)),
                                supice);
                        }
                        if (std::abs(prevp_arr(i,j,k) + pidep_arr(i,j,k)
                                     + psdep_arr(i,j,k) + pgdep_arr(i,j,k))
                            >= std::abs(satdt)) {
                            ifsat = 1;
                        }
                    }

                    if (supsat > Real(0.0) && ifsat != 1) {
                        const Real supice = satdt - prevp_arr(i,j,k)
                                           - pidep_arr(i,j,k)
                                           - psdep_arr(i,j,k)
                                           - pgdep_arr(i,j,k);
                        const Real xni0 = Real(1.0e3) * std::exp(Real(0.1) * supcol);
                        const Real roqi0 = Real(4.92e-11) * std::pow(xni0, Real(1.33));
                        pigen_arr(i,j,k) = amrex::max(
                            Real(0.0),
                            (roqi0 / den_arr(i,j,k)
                             - amrex::max(qi_arr(i,j,k), Real(0.0))) / dtcld);
                        pigen_arr(i,j,k) = amrex::min(
                            amrex::min(pigen_arr(i,j,k), satdt), supice);
                    }

                    if (qi_arr(i,j,k) > Real(0.0)) {
                        const Real qimax = Real(roqimax) / den_arr(i,j,k);
                        // WSM6-CPP TAG: PSAUT
                        //   legacy_group: G13i
                        //   process: Autoconversion to snow
                        //   compare_vars: psaut, qi, qs, den
                        psaut_arr(i,j,k) = amrex::max(
                            Real(0.0), (qi_arr(i,j,k) - qimax) / dtcld);
                    }

                    if (qs_arr(i,j,k) > Real(0.0)) {
                        const Real alpha2 = Real(1.0e-3)
                            * std::exp(Real(0.09) * (-supcol));
                        pgaut_arr(i,j,k) = amrex::min(
                            amrex::max(
                                Real(0.0), alpha2 * (qs_arr(i,j,k) - Real(qs0))),
                            qs_arr(i,j,k) / dtcld);
                    }
                }

                if (supcol < Real(0.0)) {
                    if (qs_arr(i,j,k) > Real(0.0)
                        && rhw_arr(i,j,k) < Real(1.0)) {
                        const Real coeres = rslope2_s_arr(i,j,k)
                                          * std::sqrt(rslope_s_arr(i,j,k)
                                                      * rslopeb_s_arr(i,j,k));
                        // WSM6-CPP TAG: PSEVP
                        //   legacy_group: G13j
                        //   process: Graupel evaporation/sublimation
                        //   compare_vars: psevp, qg, qv, qsat, den
                        psevp_arr(i,j,k) = (rhw_arr(i,j,k) - Real(1.0))
                            * n0sfac_arr(i,j,k)
                            * (Real(precs1) * rslope2_s_arr(i,j,k)
                               + Real(precs2) * work2_arr(i,j,k) * coeres)
                            / workdiffw_arr(i,j,k);
                        psevp_arr(i,j,k) = amrex::min(
                            amrex::max(psevp_arr(i,j,k),
                                       -qs_arr(i,j,k) / dtcld),
                            Real(0.0));
                    }

                    if (qg_arr(i,j,k) > Real(0.0)
                        && rhw_arr(i,j,k) < Real(1.0)) {
                        const Real coeres = rslope2_g_arr(i,j,k)
                                          * std::sqrt(rslope_g_arr(i,j,k)
                                                      * rslopeb_g_arr(i,j,k));
                        pgevp_arr(i,j,k) = (rhw_arr(i,j,k) - Real(1.0))
                            * (Real(precg1) * rslope2_g_arr(i,j,k)
                               + Real(precg2) * work2_arr(i,j,k) * coeres)
                            / workdiffw_arr(i,j,k);
                        pgevp_arr(i,j,k) = amrex::min(
                            amrex::max(pgevp_arr(i,j,k),
                                       -qg_arr(i,j,k) / dtcld),
                            Real(0.0));
                    }
                }
            });
            // G14: mass conservation check and state update [lines 1200-1388]
            // WSM6-CPP TAG: UPDATE
            //   legacy_group: G14
            //   process: Mass conservation and state update
            //   compare_vars: t, q, qci, qrs, qv
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real qmin_l  = Real(qmin);
                const Real qcrmin_l= Real(qcrmin);
                const Real t0c_l   = Real(t0c);

                const Real delta2 =
                    (qr_arr(i,j,k) < Real(1.0e-4) && qs_arr(i,j,k) < Real(1.0e-4))
                    ? Real(1.0) : Real(0.0);
                const Real delta3 =
                    (qr_arr(i,j,k) < Real(1.0e-4)) ? Real(1.0) : Real(0.0);

                if (t_arr(i,j,k) <= t0c_l) {
                    Real value, source, factor, xlf, xlwork2;

                    value = amrex::max(qmin_l, qc_arr(i,j,k));
                    source = (praut_arr(i,j,k) + pracw_arr(i,j,k)
                            + paacw_arr(i,j,k) + paacw_arr(i,j,k)) * dtcld;
                    if (source > value) {
                        factor = value / source;
                        praut_arr(i,j,k) = praut_arr(i,j,k) * factor;
                        pracw_arr(i,j,k) = pracw_arr(i,j,k) * factor;
                        paacw_arr(i,j,k) = paacw_arr(i,j,k) * factor;
                    }

                    value = amrex::max(qmin_l, qi_arr(i,j,k));
                    source = (psaut_arr(i,j,k) - pigen_arr(i,j,k)
                            - pidep_arr(i,j,k) + praci_arr(i,j,k)
                            + psaci_arr(i,j,k) + pgaci_arr(i,j,k)) * dtcld;
                    if (source > value) {
                        factor = value / source;
                        psaut_arr(i,j,k) = psaut_arr(i,j,k) * factor;
                        pigen_arr(i,j,k) = pigen_arr(i,j,k) * factor;
                        pidep_arr(i,j,k) = pidep_arr(i,j,k) * factor;
                        praci_arr(i,j,k) = praci_arr(i,j,k) * factor;
                        psaci_arr(i,j,k) = psaci_arr(i,j,k) * factor;
                        pgaci_arr(i,j,k) = pgaci_arr(i,j,k) * factor;
                    }

                    value = amrex::max(qmin_l, qr_arr(i,j,k));
                    source = (-praut_arr(i,j,k) - prevp_arr(i,j,k)
                            - pracw_arr(i,j,k) + piacr_arr(i,j,k)
                            + psacr_arr(i,j,k) + pgacr_arr(i,j,k)) * dtcld;
                    if (source > value) {
                        factor = value / source;
                        praut_arr(i,j,k) = praut_arr(i,j,k) * factor;
                        prevp_arr(i,j,k) = prevp_arr(i,j,k) * factor;
                        pracw_arr(i,j,k) = pracw_arr(i,j,k) * factor;
                        piacr_arr(i,j,k) = piacr_arr(i,j,k) * factor;
                        psacr_arr(i,j,k) = psacr_arr(i,j,k) * factor;
                        pgacr_arr(i,j,k) = pgacr_arr(i,j,k) * factor;
                    }

                    value = amrex::max(qmin_l, qs_arr(i,j,k));
                    source = -(psdep_arr(i,j,k) + psaut_arr(i,j,k)
                             - pgaut_arr(i,j,k) + paacw_arr(i,j,k)
                             + piacr_arr(i,j,k) * delta3
                             + praci_arr(i,j,k) * delta3
                             - pracs_arr(i,j,k) * (Real(1.0) - delta2)
                             + psacr_arr(i,j,k) * delta2
                             + psaci_arr(i,j,k) - pgacs_arr(i,j,k)) * dtcld;
                    if (source > value) {
                        factor = value / source;
                        psdep_arr(i,j,k) = psdep_arr(i,j,k) * factor;
                        psaut_arr(i,j,k) = psaut_arr(i,j,k) * factor;
                        pgaut_arr(i,j,k) = pgaut_arr(i,j,k) * factor;
                        paacw_arr(i,j,k) = paacw_arr(i,j,k) * factor;
                        piacr_arr(i,j,k) = piacr_arr(i,j,k) * factor;
                        praci_arr(i,j,k) = praci_arr(i,j,k) * factor;
                        psaci_arr(i,j,k) = psaci_arr(i,j,k) * factor;
                        pracs_arr(i,j,k) = pracs_arr(i,j,k) * factor;
                        psacr_arr(i,j,k) = psacr_arr(i,j,k) * factor;
                        pgacs_arr(i,j,k) = pgacs_arr(i,j,k) * factor;
                    }

                    value = amrex::max(qmin_l, qg_arr(i,j,k));
                    source = -(pgdep_arr(i,j,k) + pgaut_arr(i,j,k)
                             + piacr_arr(i,j,k) * (Real(1.0) - delta3)
                             + praci_arr(i,j,k) * (Real(1.0) - delta3)
                             + psacr_arr(i,j,k) * (Real(1.0) - delta2)
                             + pracs_arr(i,j,k) * (Real(1.0) - delta2)
                             + pgaci_arr(i,j,k) + paacw_arr(i,j,k)
                             + pgacr_arr(i,j,k) + pgacs_arr(i,j,k)) * dtcld;
                    if (source > value) {
                        factor = value / source;
                        pgdep_arr(i,j,k) = pgdep_arr(i,j,k) * factor;
                        pgaut_arr(i,j,k) = pgaut_arr(i,j,k) * factor;
                        piacr_arr(i,j,k) = piacr_arr(i,j,k) * factor;
                        praci_arr(i,j,k) = praci_arr(i,j,k) * factor;
                        psacr_arr(i,j,k) = psacr_arr(i,j,k) * factor;
                        pracs_arr(i,j,k) = pracs_arr(i,j,k) * factor;
                        paacw_arr(i,j,k) = paacw_arr(i,j,k) * factor;
                        pgaci_arr(i,j,k) = pgaci_arr(i,j,k) * factor;
                        pgacr_arr(i,j,k) = pgacr_arr(i,j,k) * factor;
                        pgacs_arr(i,j,k) = pgacs_arr(i,j,k) * factor;
                    }

                    work2_arr(i,j,k) = -(prevp_arr(i,j,k) + psdep_arr(i,j,k)
                                       + pgdep_arr(i,j,k) + pigen_arr(i,j,k)
                                       + pidep_arr(i,j,k));
                    qv_arr(i,j,k) = qv_arr(i,j,k) + work2_arr(i,j,k) * dtcld;
                    qc_arr(i,j,k) = amrex::max(
                        qc_arr(i,j,k) - (praut_arr(i,j,k) + pracw_arr(i,j,k)
                                       + paacw_arr(i,j,k) + paacw_arr(i,j,k))
                                       * dtcld,
                        Real(0.0));
                    qr_arr(i,j,k) = amrex::max(
                        qr_arr(i,j,k) + (praut_arr(i,j,k) + pracw_arr(i,j,k)
                                       + prevp_arr(i,j,k) - piacr_arr(i,j,k)
                                       - pgacr_arr(i,j,k) - psacr_arr(i,j,k))
                                       * dtcld,
                        Real(0.0));
                    qi_arr(i,j,k) = amrex::max(
                        qi_arr(i,j,k) - (psaut_arr(i,j,k) + praci_arr(i,j,k)
                                       + psaci_arr(i,j,k) + pgaci_arr(i,j,k)
                                       - pigen_arr(i,j,k) - pidep_arr(i,j,k))
                                       * dtcld,
                        Real(0.0));
                    qs_arr(i,j,k) = amrex::max(
                        qs_arr(i,j,k) + (psdep_arr(i,j,k) + psaut_arr(i,j,k)
                                       + paacw_arr(i,j,k) - pgaut_arr(i,j,k)
                                       + piacr_arr(i,j,k) * delta3
                                       + praci_arr(i,j,k) * delta3
                                       + psaci_arr(i,j,k) - pgacs_arr(i,j,k)
                                       - pracs_arr(i,j,k) * (Real(1.0) - delta2)
                                       + psacr_arr(i,j,k) * delta2) * dtcld,
                        Real(0.0));
                    qg_arr(i,j,k) = amrex::max(
                        qg_arr(i,j,k) + (pgdep_arr(i,j,k) + pgaut_arr(i,j,k)
                                       + piacr_arr(i,j,k) * (Real(1.0) - delta3)
                                       + praci_arr(i,j,k) * (Real(1.0) - delta3)
                                       + psacr_arr(i,j,k) * (Real(1.0) - delta2)
                                       + pracs_arr(i,j,k) * (Real(1.0) - delta2)
                                       + pgaci_arr(i,j,k) + paacw_arr(i,j,k)
                                       + pgacr_arr(i,j,k) + pgacs_arr(i,j,k))
                                       * dtcld,
                        Real(0.0));
                    xlf = Real(xls) - xl_arr(i,j,k);
                    xlwork2 = -Real(xls) * (psdep_arr(i,j,k) + pgdep_arr(i,j,k)
                                          + pidep_arr(i,j,k) + pigen_arr(i,j,k))
                            - xl_arr(i,j,k) * prevp_arr(i,j,k)
                            - xlf * (piacr_arr(i,j,k) + paacw_arr(i,j,k)
                                   + paacw_arr(i,j,k) + pgacr_arr(i,j,k)
                                   + psacr_arr(i,j,k));
                    t_arr(i,j,k) = t_arr(i,j,k) - xlwork2 / cpm_arr(i,j,k) * dtcld;
                } else {
                    Real value, source, factor, xlf, xlwork2;

                    value = amrex::max(qmin_l, qc_arr(i,j,k));
                    source = (praut_arr(i,j,k) + pracw_arr(i,j,k)
                            + paacw_arr(i,j,k) + paacw_arr(i,j,k)) * dtcld;
                    if (source > value) {
                        factor = value / source;
                        praut_arr(i,j,k) = praut_arr(i,j,k) * factor;
                        pracw_arr(i,j,k) = pracw_arr(i,j,k) * factor;
                        paacw_arr(i,j,k) = paacw_arr(i,j,k) * factor;
                    }

                    value = amrex::max(qmin_l, qr_arr(i,j,k));
                    source = (-paacw_arr(i,j,k) - praut_arr(i,j,k)
                            + pseml_arr(i,j,k) + pgeml_arr(i,j,k)
                            - pracw_arr(i,j,k) - paacw_arr(i,j,k)
                            - prevp_arr(i,j,k)) * dtcld;
                    if (source > value) {
                        factor = value / source;
                        praut_arr(i,j,k) = praut_arr(i,j,k) * factor;
                        prevp_arr(i,j,k) = prevp_arr(i,j,k) * factor;
                        pracw_arr(i,j,k) = pracw_arr(i,j,k) * factor;
                        paacw_arr(i,j,k) = paacw_arr(i,j,k) * factor;
                        pseml_arr(i,j,k) = pseml_arr(i,j,k) * factor;
                        pgeml_arr(i,j,k) = pgeml_arr(i,j,k) * factor;
                    }

                    value = amrex::max(qcrmin_l, qs_arr(i,j,k));
                    source = (pgacs_arr(i,j,k) - pseml_arr(i,j,k)
                            - psevp_arr(i,j,k)) * dtcld;
                    if (source > value) {
                        factor = value / source;
                        pgacs_arr(i,j,k) = pgacs_arr(i,j,k) * factor;
                        psevp_arr(i,j,k) = psevp_arr(i,j,k) * factor;
                        pseml_arr(i,j,k) = pseml_arr(i,j,k) * factor;
                    }

                    value = amrex::max(qcrmin_l, qg_arr(i,j,k));
                    source = -(pgacs_arr(i,j,k) + pgevp_arr(i,j,k)
                             + pgeml_arr(i,j,k)) * dtcld;
                    if (source > value) {
                        factor = value / source;
                        pgacs_arr(i,j,k) = pgacs_arr(i,j,k) * factor;
                        pgevp_arr(i,j,k) = pgevp_arr(i,j,k) * factor;
                        pgeml_arr(i,j,k) = pgeml_arr(i,j,k) * factor;
                    }

                    work2_arr(i,j,k) = -(prevp_arr(i,j,k) + psevp_arr(i,j,k)
                                       + pgevp_arr(i,j,k));
                    qv_arr(i,j,k) = qv_arr(i,j,k) + work2_arr(i,j,k) * dtcld;
                    qc_arr(i,j,k) = amrex::max(
                        qc_arr(i,j,k) - (praut_arr(i,j,k) + pracw_arr(i,j,k)
                                       + paacw_arr(i,j,k) + paacw_arr(i,j,k))
                                       * dtcld,
                        Real(0.0));
                    qr_arr(i,j,k) = amrex::max(
                        qr_arr(i,j,k) + (praut_arr(i,j,k) + pracw_arr(i,j,k)
                                       + prevp_arr(i,j,k) + paacw_arr(i,j,k)
                                       + paacw_arr(i,j,k) - pseml_arr(i,j,k)
                                       - pgeml_arr(i,j,k)) * dtcld,
                        Real(0.0));
                    qs_arr(i,j,k) = amrex::max(
                        qs_arr(i,j,k) + (psevp_arr(i,j,k) - pgacs_arr(i,j,k)
                                       + pseml_arr(i,j,k)) * dtcld,
                        Real(0.0));
                    qg_arr(i,j,k) = amrex::max(
                        qg_arr(i,j,k) + (pgacs_arr(i,j,k) + pgevp_arr(i,j,k)
                                       + pgeml_arr(i,j,k)) * dtcld,
                        Real(0.0));
                    xlf = Real(xls) - xl_arr(i,j,k);
                    xlwork2 = -xl_arr(i,j,k) * (prevp_arr(i,j,k)
                                              + psevp_arr(i,j,k)
                                              + pgevp_arr(i,j,k))
                            - xlf * (pseml_arr(i,j,k) + pgeml_arr(i,j,k));
                    t_arr(i,j,k) = t_arr(i,j,k) - xlwork2 / cpm_arr(i,j,k) * dtcld;
                }
            });
            // G15: second qsat computation [lines 1390-1420]
            // WSM6-CPP TAG: QSAT2
            //   legacy_group: G15
            //   process: Second saturation mixing ratio computation
            //   compare_vars: qs, qvs, den, denfac, t, p
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real ttp   = Real(t0c) + Real(0.01);
                const Real dldt  = Real(cpv) - Real(cliq);
                const Real xa    = -dldt / Real(rv);
                const Real xb    = xa + Real(xlv0) / (Real(rv) * ttp);

                Real tr = ttp / t_arr(i,j,k);
                Real qsw = Real(psat) * std::exp(std::log(tr) * xa)
                                     * std::exp(xb * (Real(1.0) - tr));
                qsw = amrex::min(qsw, Real(0.99) * p_arr(i,j,k));
                qsatw_arr(i,j,k) = Real(ep2) * qsw / (p_arr(i,j,k) - qsw);
                qsatw_arr(i,j,k) = amrex::max(qsatw_arr(i,j,k), Real(qmin));
            });
            // G16: pcond condensational/evaporational update [lines 1427-1437]
            // WSM6-CPP TAG: PCOND
            //   legacy_group: G16
            //   process: Condensation/evaporation update
            //   compare_vars: pcond, t, qv, qc, qsat
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real workcond = wsm6_conden(
                    t_arr(i,j,k), qv_arr(i,j,k), qsatw_arr(i,j,k),
                    xl_arr(i,j,k), cpm_arr(i,j,k), Real(qmin), Real(rv));
                const Real work2loc = qc_arr(i,j,k) + workcond;
                static_cast<void>(work2loc);
                pcond_arr(i,j,k) = amrex::min(
                    amrex::max(workcond / dtcld, Real(0.0)),
                    amrex::max(qv_arr(i,j,k), Real(0.0)) / dtcld);
                if (qc_arr(i,j,k) > Real(0.0) && workcond < Real(0.0)) {
                    pcond_arr(i,j,k) = amrex::max(workcond, -qc_arr(i,j,k)) / dtcld;
                }
                qv_arr(i,j,k) = qv_arr(i,j,k) - pcond_arr(i,j,k) * dtcld;
                qc_arr(i,j,k) = amrex::max(
                    qc_arr(i,j,k) + pcond_arr(i,j,k) * dtcld,
                    Real(0.0));
                t_arr(i,j,k) = t_arr(i,j,k)
                             + pcond_arr(i,j,k) * xl_arr(i,j,k)
                             / cpm_arr(i,j,k) * dtcld;
            });
            // G17: padding for small values [lines 1444-1449]
            // WSM6-CPP TAG: CLIP
            //   legacy_group: G17
            //   process: Padding/clipping for small values
            //   compare_vars: qv, qc, qr, qi, qs, qg
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (qc_arr(i,j,k) <= Real(qmin)) qc_arr(i,j,k) = Real(0.0);
                if (qi_arr(i,j,k) <= Real(qmin)) qi_arr(i,j,k) = Real(0.0);
            });

        }
#ifdef ERF_USE_WSM6_FORT
        }
#endif
        ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            rain_arr(i,j,klo) = rainacc_arr(i,j,0);
            snow_arr(i,j,klo) = snowacc_arr(i,j,0);
            graup_arr(i,j,klo) = graupacc_arr(i,j,0);
        });
    }
}

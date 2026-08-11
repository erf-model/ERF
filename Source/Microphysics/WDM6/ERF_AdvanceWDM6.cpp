#include "ERF_WDM6.H"
#include <AMReX_IArrayBox.H>
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

#ifdef ERF_USE_WDM6_FORT
#include "ERF_WDM6_Fortran_Interface.H"
#endif

using namespace amrex;

// ---------------------------------------------------------------
// WDM6 device-callable helper functions
// Statement functions from WRF WDM6
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_cpmcal (Real x, Real qmin_arg, Real cpd_arg, Real cpv_arg) {
    return cpd_arg*(Real(1.0)-amrex::max(x,qmin_arg))
          +amrex::max(x,qmin_arg)*cpv_arg;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_xlcal (Real x, Real xlv0_arg, Real xlv1_arg, Real t0c_arg) {
    return xlv0_arg - xlv1_arg*(x - t0c_arg);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_diffus (Real x, Real y) {
    return (Real(8.794e-5f)*std::exp(std::log(x)*Real(1.81f)))/y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_viscos (Real x, Real y) {
    return Real(1.496e-6f)*(x*std::sqrt(x))/(x+Real(120.0))/y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_xka (Real x, Real y) {
    return Real(1.414e3)*wdm6_viscos(x,y)*y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_diffac (Real a, Real b, Real c, Real d, Real e,
                   Real rv_arg) {
    return d*a*a/(wdm6_xka(c,d)*rv_arg*c*c)
          + Real(1.0)/(e*wdm6_diffus(c,b));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_venfac (Real a, Real b, Real c, Real den0_arg) {
    // Fortran: exp(log((viscos(b,c)/diffus(b,a)))*((.3333333))) / sqrt(viscos) * sqrt(sqrt(den0/c))
    return std::exp(std::log(wdm6_viscos(b,c)/wdm6_diffus(b,a))
                   *Real(0.3333333f))
          /std::sqrt(wdm6_viscos(b,c))
          *std::sqrt(std::sqrt(den0_arg/c));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_conden (Real a, Real b, Real c, Real d, Real e,
                   Real qmin_arg, Real rv_arg) {
    return (amrex::max(b,qmin_arg)-c)
          /(Real(1.0)+d*d/(rv_arg*e)*c/(a*a));
}

// ---------------------------------------------------------------
// WDM6 double-moment slope parameter functions
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_lamdar (Real qr, Real den, Real nr, Real pidnr_arg) {
    // Rain slope parameter using number concentration nr
    // qr = rain mixing ratio (kg/kg)
    // den = air density (kg/m^3)
    // nr = rain number concentration (#/kg)
    // Returns 1/lambda (m)
    if (qr <= Real(1.e-9) || nr <= Real(1.e-2)) {
        return Real(1.0) / Real(5.0e4);  // lamdarmax
    }
    return std::pow((pidnr_arg * nr * den) / (den * qr), Real(1.0)/Real(3.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_lamdac (Real qc, Real den, Real nc, Real pidnc_arg) {
    // Cloud droplet slope parameter using number concentration nc
    // qc = cloud mixing ratio (kg/kg)
    // den = air density (kg/m^3)
    // nc = cloud droplet number concentration (#/kg)
    // Returns 1/lambda (m)
    if (qc <= Real(1.e-9) || nc <= Real(1.e1)) {
        return Real(1.0) / Real(5.0e5);  // lamdacmax
    }
    return std::pow((pidnc_arg * nc * den) / (den * qc), Real(1.0)/Real(3.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_rslopec_exact (Real qc, Real den, Real nc, Real pidnc_arg) {
    return std::exp(std::log((pidnc_arg * nc) / (den * qc)) * Real(0.33333333));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_xni_exact (Real qi, Real den, Real qmin_arg) {
    Real temp = den * amrex::max(qi, qmin_arg);
    temp = std::sqrt(std::sqrt(temp * temp * temp));
    return amrex::min(amrex::max(Real(5.38e7) * temp, Real(1.e3)), Real(1.e6));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wdm6_slope_rain_cell (Real qr, Real nr, Real den, Real denfac,
                           Real qcrmin_arg, Real nrmin_arg,
                           Real rslopermax_arg, Real rsloperbmax_arg,
                           Real rsloper2max_arg, Real rsloper3max_arg,
                           Real bvtr_arg, Real pvtr_arg, Real pvtrn_arg,
                           Real pidnr_arg,
                           Real& rslope, Real& rslopeb,
                           Real& rslope2, Real& rslope3,
                           Real& vt, Real& vtn)
{
    if (qr <= qcrmin_arg || nr <= nrmin_arg) {
        rslope  = rslopermax_arg;
        rslopeb = rsloperbmax_arg;
        rslope2 = rsloper2max_arg;
        rslope3 = rsloper3max_arg;
    } else {
        rslope  = amrex::min(Real(1.0) / wdm6_lamdar(qr, den, nr, pidnr_arg), Real(1.e-3));
        rslopeb = std::pow(rslope, bvtr_arg);
        rslope2 = rslope * rslope;
        rslope3 = rslope2 * rslope;
    }
    vt = pvtr_arg * rslopeb * denfac;
    vtn = pvtrn_arg * rslopeb * denfac;
    if (qr <= Real(0.0)) vt = Real(0.0);
    if (nr <= Real(0.0)) vtn = Real(0.0);
}

// ---------------------------------------------------------------
// Ice slope parameter functions (single-moment from WSM6)
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_lamdas (Real x, Real y, Real z, Real pidn0s_arg) {
    // Snow slope parameter (single-moment)
    return std::sqrt(std::sqrt(pidn0s_arg*z/(x*y)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_lamdag (Real x, Real y, Real pidn0g_arg) {
    // Graupel slope parameter (single-moment)
    return std::sqrt(std::sqrt(pidn0g_arg/(x*y)));
}

// ---------------------------------------------------------------
// Full slope/terminal velocity functions for ice species
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wdm6_slope_snow_cell (Real qs, Real den, Real denfac, Real t,
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
        rslope  = Real(1.0)/wdm6_lamdas(qs,den,n0sfac,pidn0s_arg);
        rslopeb = std::pow(rslope,bvts_arg);
        rslope2 = rslope*rslope;
        rslope3 = rslope2*rslope;
    }
    vt = pvts_arg*rslopeb*denfac;
    if (qs <= Real(0.0)) vt = Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wdm6_slope_graup_cell (Real qg, Real den, Real denfac,
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
        rslope  = Real(1.0)/wdm6_lamdag(qg,den,pidn0g_arg);
        rslopeb = std::pow(rslope,bvtg_arg);
        rslope2 = rslope*rslope;
        rslope3 = rslope2*rslope;
    }
    vt = pvtg_arg*rslopeb*denfac;
    if (qg <= Real(0.0)) vt = Real(0.0);
}

#if !defined(AMREX_USE_GPU)
void wdm6_nislfv_rain_plm6_column (
    int km,
    const Real* dz,
    const Real* den,
    const Real* denfac,
    const Real* tk,
    Real* ww,
    Real* rql,
    Real* rql2,
    Real& precip1,
    Real& precip2,
    Real dt,
    int iter,
    Real pidn0s,
    Real pidn0g,
    Real qcrmin,
    Real alpha,
    Real n0smax,
    Real n0s,
    Real t0c,
    Real rslopesmax,
    Real rslopesbmax,
    Real rslopes2max,
    Real rslopes3max,
    Real bvts,
    Real pvts,
    Real rslopegmax,
    Real rslopegbmax,
    Real rslopeg2max,
    Real rslopeg3max,
    Real bvtg,
    Real pvtg)
{
    std::vector<Real> qq(km), qq2(km), wd(km), wa(km), wa2(km), was(km);
    std::vector<Real> wi(km + 1), zi(km + 1), za(km + 1), dza(km + 1);
    std::vector<Real> qa(km + 1), qa2(km + 1), qn(km), qn2(km);
    std::vector<Real> qr(km), qr2(km), qmi(km + 1), qpi(km + 1);

    Real allold = Real(0.0);
    for (int k = 0; k < km; ++k) {
        qq[k] = rql[k];
        qq2[k] = rql2[k];
        wd[k] = ww[k];
        allold += qq[k] + qq2[k];
    }

    precip1 = Real(0.0);
    precip2 = Real(0.0);
    if (allold <= Real(0.0)) {
        return;
    }

    zi[0] = Real(0.0);
    for (int k = 0; k < km; ++k) {
        zi[k + 1] = zi[k] + dz[k];
    }

    auto rebuild = [&]() {
        wi[0] = ww[0];
        if (km > 1) {
            wi[1] = Real(0.5) * (ww[1] + ww[0]);
            for (int k = 2; k < km - 1; ++k) {
                wi[k] = Real(9.0) / Real(16.0) * (ww[k] + ww[k - 1])
                      - Real(1.0) / Real(16.0) * (ww[k + 1] + ww[k - 2]);
            }
            wi[km - 1] = Real(0.5) * (ww[km - 1] + ww[km - 2]);
        }
        wi[km] = ww[km - 1];
        for (int k = 1; k < km; ++k) {
            if (ww[k] == Real(0.0)) wi[k] = ww[k - 1];
        }

        constexpr Real con1 = Real(0.05);
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
            qa[k] = qq[k] * dz[k] / dza[k];
            qa2[k] = qq2[k] * dz[k] / dza[k];
            qr[k] = qa[k] / den[k];
            qr2[k] = qa2[k] / den[k];
        }
        dza[km] = zi[km] - za[km];
        qa[km] = Real(0.0);
        qa2[km] = Real(0.0);
    };

    rebuild();

    if (iter > 0) {
        for (int k = 0; k < km; ++k) {
            Real rslope, rslopeb, rslope2, rslope3, vt, n0sfac_dummy;
            wdm6_slope_snow_cell(qr[k], den[k], denfac[k], tk[k],
                                 pidn0s, alpha, n0smax, n0s, t0c, qcrmin,
                                 rslopesmax, rslopesbmax, rslopes2max, rslopes3max,
                                 bvts, pvts, rslope, rslopeb, rslope2, rslope3, vt,
                                 n0sfac_dummy);
            wa[k] = vt;
            wdm6_slope_graup_cell(qr2[k], den[k], denfac[k],
                                  pidn0g, qcrmin, rslopegmax, rslopegbmax,
                                  rslopeg2max, rslopeg3max, bvtg, pvtg,
                                  rslope, rslopeb, rslope2, rslope3, vt);
            wa2[k] = vt;
        }
        for (int k = 0; k < km; ++k) {
            const Real tmpq = amrex::max(qr[k] + qr2[k], Real(1.0e-15));
            wa[k] = (tmpq > Real(1.0e-15))
                ? (wa[k] * qr[k] + wa2[k] * qr2[k]) / tmpq
                : Real(0.0);
            ww[k] = Real(0.5) * (wd[k] + wa[k]);
            was[k] = wa[k];
        }
        rebuild();
    }

    auto advect_one = [&](std::vector<Real>& qa_src, std::vector<Real>& qn_dst, Real& precip_dst) {
        qpi[0] = qa_src[0];
        qmi[0] = qa_src[0];
        qmi[km] = qa_src[km];
        qpi[km] = qa_src[km];
        for (int k = 1; k < km; ++k) {
            const Real dip = (qa_src[k + 1] - qa_src[k]) / (dza[k + 1] + dza[k]);
            const Real dim = (qa_src[k] - qa_src[k - 1]) / (dza[k - 1] + dza[k]);
            if (dip * dim <= Real(0.0)) {
                qmi[k] = qa_src[k];
                qpi[k] = qa_src[k];
            } else {
                qpi[k] = qa_src[k] + Real(0.5) * (dip + dim) * dza[k];
                qmi[k] = Real(2.0) * qa_src[k] - qpi[k];
                if (qpi[k] < Real(0.0) || qmi[k] < Real(0.0)) {
                    qpi[k] = qa_src[k];
                    qmi[k] = qa_src[k];
                }
            }
        }

        std::fill(qn_dst.begin(), qn_dst.end(), Real(0.0));
        int kb = 0;
        int kt = 0;
        for (int k = 0; k < km; ++k) {
            if (zi[k] >= za[km]) break;

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
                const Real qqd = Real(0.5) * (qpi[kb] - qmi[kb]);
                qn_dst[k] = (qqd * th * th + qmi[kb] * th
                           - qqd * tl * tl - qmi[kb] * tl) / (th - tl);
            } else if (kt > kb) {
                const Real tl = (zi[k] - za[kb]) / dza[kb];
                const Real qqd = Real(0.5) * (qpi[kb] - qmi[kb]);
                const Real qql = qqd * tl * tl + qmi[kb] * tl;
                const Real dql = qa_src[kb] - qql;
                Real zsum = (Real(1.0) - tl) * dza[kb];
                Real qsum = dql * dza[kb];
                for (int m = kb + 1; m < kt; ++m) {
                    zsum += dza[m];
                    qsum += qa_src[m] * dza[m];
                }
                const Real th = (zi[k + 1] - za[kt]) / dza[kt];
                const Real dqh = Real(0.5) * (qpi[kt] - qmi[kt]) * th * th + qmi[kt] * th;
                zsum += th * dza[kt];
                qsum += dqh * dza[kt];
                qn_dst[k] = qsum / zsum;
            }
        }

        precip_dst = Real(0.0);
        for (int k = 0; k < km; ++k) {
            if (za[k] < Real(0.0) && za[k + 1] < Real(0.0)) {
                precip_dst += qa_src[k] * dza[k];
            } else if (za[k] < Real(0.0) && za[k + 1] >= Real(0.0)) {
                precip_dst += qa_src[k] * (Real(0.0) - za[k]);
                break;
            } else {
                break;
            }
        }
    };

    advect_one(qa, qn, precip1);
    advect_one(qa2, qn2, precip2);

    for (int k = 0; k < km; ++k) {
        rql[k] = qn[k];
        rql2[k] = qn2[k];
    }
}
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_mean_droplet_diameter (Real qc, Real nc, Real den, Real pidnc_arg) {
    // Volume-weighted mean diameter of cloud droplets
    // Returns diameter in meters
    if (nc < Real(1.e1) || qc < Real(1.e-9)) return Real(0.0);

    // Mean mass = rho * qc / nc
    // Mean volume = mean mass / rho_water
    // Mean diameter = (6 * volume / pi)^(1/3)
    Real mean_mass = den * qc / (nc * den);  // kg per droplet
    Real mean_volume = mean_mass / Real(1000.0);  // m^3 per droplet
    Real diameter = std::pow(Real(6.0) * mean_volume / Real(3.14159265359), Real(1.0)/Real(3.0));

    return diameter;
}

// ---------------------------------------------------------------
// WDM6 CCN activation (simplified Abdul-Razzak & Ghan 2000)
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wdm6_ccn_activation (
    Real& nc,           // cloud droplet number (#/kg) - in/out
    Real& nn,           // aerosol number (#/kg) - in/out
    const Real qv,      // water vapor mixing ratio (kg/kg)
    const Real qc,      // cloud water mixing ratio (kg/kg)
    const Real qvs,     // saturation mixing ratio (kg/kg)
    const Real temp,    // temperature (K)
    const Real w,       // vertical velocity (m/s)
    const Real dt,      // timestep (s)
    const Real ccn0,    // background CCN (#/m^3)
    const Real den      // air density (kg/m^3)
) {
    // Only activate if supersaturated and cloud exists
    if (qv <= qvs || qc < Real(1.e-8)) return;

    // Supersaturation
    Real supersaturation = qv / qvs - Real(1.0);
    supersaturation = amrex::min(supersaturation, Real(0.0048));  // satmax - 1

    // Critical supersaturation based on updraft strength
    // From WRF: stronger updrafts -> more activation
    Real s_crit = Real(0.6) * std::pow(amrex::max(w, Real(0.01)), Real(0.5));

    if (supersaturation > s_crit && nn > Real(0.0)) {
        // Activation rate: convert aerosols to droplets
        Real nc_activate = amrex::min(nn, Real(0.1) * (supersaturation - s_crit) / s_crit * nn);

        // Update concentrations
        nc += nc_activate;
        nn -= nc_activate;

        // Enforce physical bounds
        nc = amrex::max(nc, Real(1.e1));
        nn = amrex::max(nn, Real(0.0));
    }
}

// ---------------------------------------------------------------
// Main Advance routine
// ---------------------------------------------------------------

void WDM6::Advance(const Real& dt_advance,
                   const SolverChoice& solverChoice)
{
    // ---------------------------------------------------------------
    // Dual-mode implementation following WSM6 pattern:
    // - With ERF_USE_WDM6_FORT: Call Fortran bridge (CPU-only)
    // - Without: Use C++ GPU kernels (not yet implemented)
    // ---------------------------------------------------------------

    static int call_count = 0;
    call_count++;
    const bool first_call = (call_count == 1);

#ifdef ERF_USE_WDM6_FORT
    // Fortran bridge mode - initialize once
    static bool wdm6_inited = false;
    if (!wdm6_inited) {
        constexpr double den0 = 1.28;                 // Standard dry-air density (kg/m^3)
        constexpr double denr = static_cast<double>(rhoh2o);
        constexpr double dens = static_cast<double>(rhos);
        constexpr double cl = static_cast<double>(Cp_l);
        constexpr double cpv = static_cast<double>(Cp_v);
        const double ccn0 = static_cast<double>(m_ccn0);
        constexpr int hail_opt = 0;                   // Graupel mode
        mp_wdm6_init_c(den0, denr, dens, cl, cpv, ccn0, hail_opt);
        wdm6_inited = true;
        if (first_call) {
            amrex::Print() << "WDM6 Fortran bridge initialized\n";
        }
    }
#endif

    int microphysics_debug = 0;
    std::vector<int> micro_diag_target_column;
    {
        amrex::ParmParse pp("erf");
        pp.query("microphysics_debug", microphysics_debug);
        pp.queryarr("micro_diag_target_column", micro_diag_target_column);
    }
    microphysics_debug = std::max(0, std::min(2, microphysics_debug));
#ifdef ERF_USE_WDM6_FORT
    bool use_wdm6_cpp_answer = false;
    {
        amrex::ParmParse pp("erf");
        pp.query("use_wdm6_cpp_answer", use_wdm6_cpp_answer);
    }
    const bool run_wdm6_fort = !use_wdm6_cpp_answer;
#endif

    // Physical constants
    constexpr double g = static_cast<double>(CONST_GRAV);
    constexpr double cpd = static_cast<double>(Cp_d);
    constexpr double cpv = static_cast<double>(Cp_v);
    constexpr double rd = static_cast<double>(R_d);
    constexpr double rv = static_cast<double>(R_v);
    constexpr double t0c = 273.15;
    constexpr double ep1 = static_cast<double>(R_v / R_d - one);
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
    const double ccn0 = static_cast<double>(m_ccn0);

    for (MFIter mfi(*mic_fab_vars[MicVar_WDM6::qv], TileNoZ()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const Box fab_box = mfi.fabbox();

        // Get array pointers for all WDM6 variables
        auto const& t_arr = mic_fab_vars[MicVar_WDM6::tabs]->array(mfi);
        auto const& qv_arr = mic_fab_vars[MicVar_WDM6::qv]->array(mfi);
        auto const& qc_arr = mic_fab_vars[MicVar_WDM6::qc]->array(mfi);
        auto const& qi_arr = mic_fab_vars[MicVar_WDM6::qi]->array(mfi);
        auto const& qr_arr = mic_fab_vars[MicVar_WDM6::qr]->array(mfi);
        auto const& qs_arr = mic_fab_vars[MicVar_WDM6::qs]->array(mfi);
        auto const& qg_arr = mic_fab_vars[MicVar_WDM6::qg]->array(mfi);
        auto const& nn_arr = mic_fab_vars[MicVar_WDM6::nn]->array(mfi);  // Aerosol number
        auto const& nc_arr = mic_fab_vars[MicVar_WDM6::nc]->array(mfi);  // Cloud droplet number
        auto const& nr_arr = mic_fab_vars[MicVar_WDM6::nr]->array(mfi);  // Rain drop number
        auto const& den_arr = mic_fab_vars[MicVar_WDM6::rho]->array(mfi);
        auto const& p_arr = mic_fab_vars[MicVar_WDM6::pres]->array(mfi);
        auto rain_arr = mic_fab_vars[MicVar_WDM6::rain_accum]->array(mfi);
        auto snow_arr = mic_fab_vars[MicVar_WDM6::snow_accum]->array(mfi);
        auto graup_arr = mic_fab_vars[MicVar_WDM6::graup_accum]->array(mfi);

        const int ilo = box.smallEnd(0);
        const int ihi = box.bigEnd(0);
        const int jlo = box.smallEnd(1);
        const int jhi = box.bigEnd(1);
        const int klo = box.smallEnd(2);
        const int khi = box.bigEnd(2);

        const int imlo = fab_box.smallEnd(0);
        const int imhi = fab_box.bigEnd(0);
        const int jmlo = fab_box.smallEnd(1);
        const int jmhi = fab_box.bigEnd(1);
        const int kmlo = fab_box.smallEnd(2);
        const int kmhi = fab_box.bigEnd(2);
        const bool has_target_override = (micro_diag_target_column.size() == 2);
        const int diag_i = has_target_override ? micro_diag_target_column[0] : ilo;
        const int diag_j = has_target_override ? micro_diag_target_column[1] : jlo;

#if defined(ERF_USE_WDM6_FORT) && defined(AMREX_USE_GPU)
        Arena* Arena_Used = run_wdm6_fort ? The_Pinned_Arena() : The_Async_Arena();
#else
        Arena* Arena_Used = The_Async_Arena();
#endif

#ifdef ERF_USE_WDM6_FORT
        if (run_wdm6_fort) {
        // Fortran bridge path
        // Create delz array (cell thickness)
        const Real dz_val = m_geom.CellSize(2);
        FArrayBox delz_fab(fab_box, 1, Arena_Used);
        auto const& delz_arr = delz_fab.array();
        ParallelFor(fab_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            delz_arr(i,j,k) = dz_val;
        });

        // Create landmask array (xland: 0=water, 1=land)
        // TODO: Get from ERF's lmask_lev when available
        // For now, default to land (continental CCN)
        FArrayBox xland_fab(Box(IntVect(imlo,jmlo,0), IntVect(imhi,jmhi,0)), 1, Arena_Used);
        auto const& xland_arr = xland_fab.array();
        ParallelFor(Box(IntVect(imlo,jmlo,0), IntVect(imhi,jmhi,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            xland_arr(i,j,k) = Real(1.0);  // Default to land
        });

        // Create 2D accumulation arrays
        Box box2d(IntVect(ilo,jlo,0), IntVect(ihi,jhi,0));
        FArrayBox rainacc_fab(box2d, 1, Arena_Used);
        FArrayBox rainncv_fab(box2d, 1, Arena_Used);
        FArrayBox sr_fab(box2d, 1, Arena_Used);
        FArrayBox snowacc_fab(box2d, 1, Arena_Used);
        FArrayBox snowncv_fab(box2d, 1, Arena_Used);
        FArrayBox graupacc_fab(box2d, 1, Arena_Used);
        FArrayBox graupelncv_fab(box2d, 1, Arena_Used);

        auto const& rainacc_arr = rainacc_fab.array();
        auto const& rainncv_arr = rainncv_fab.array();
        auto const& sr_arr = sr_fab.array();
        auto const& snowacc_arr = snowacc_fab.array();
        auto const& snowncv_arr = snowncv_fab.array();
        auto const& graupacc_arr = graupacc_fab.array();
        auto const& graupelncv_arr = graupelncv_fab.array();

        // Initialize 2D arrays to zero
        ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            rainacc_arr(i,j,k) = Real(0.0);
            rainncv_arr(i,j,k) = Real(0.0);
            sr_arr(i,j,k) = Real(0.0);
            snowacc_arr(i,j,k) = Real(0.0);
            snowncv_arr(i,j,k) = Real(0.0);
            graupacc_arr(i,j,k) = Real(0.0);
            graupelncv_arr(i,j,k) = Real(0.0);
        });

        // (Tile-based diagnostics removed - using global diagnostics instead)

        // Call Fortran WDM6
        mp_wdm6_run_c(
            t_arr.dataPtr(),
            qv_arr.dataPtr(), qc_arr.dataPtr(), qi_arr.dataPtr(),
            qr_arr.dataPtr(), qs_arr.dataPtr(), qg_arr.dataPtr(),
            nn_arr.dataPtr(), nc_arr.dataPtr(), nr_arr.dataPtr(),  // WDM6 number concentrations
            den_arr.dataPtr(), p_arr.dataPtr(), delz_arr.dataPtr(),
            static_cast<double>(dt_advance), g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin,
            xls, xlv0, xlf0, den0, denr, cliq, cice, psat,
            ccn0, xland_arr.dataPtr(),  // WDM6 parameters
            rainacc_arr.dataPtr(), rainncv_arr.dataPtr(), sr_arr.dataPtr(),
            snowacc_arr.dataPtr(), snowncv_arr.dataPtr(),
            graupacc_arr.dataPtr(), graupelncv_arr.dataPtr(),
            imlo, imhi, jmlo, jmhi, kmlo, kmhi,
            ilo, ihi, jlo, jhi, klo, khi,
            microphysics_debug, diag_i, diag_j);

        // CRITICAL: Convert updated temperature back to potential temperature
        // The Fortran WDM6 modifies t_arr (absolute temperature) due to latent heating/cooling
        // from condensation, evaporation, freezing, melting, etc.
        // ERF stores theta (potential temperature), so we must convert back: theta = T / exner
        // This matches WRF's conversion: th(i,k,j) = t(i,k) / pii(i,k,j)
        auto const& theta_arr = mic_fab_vars[MicVar_WDM6::theta]->array(mfi);
        constexpr Real p0 = 1.e5;       // Reference pressure (Pa)
        constexpr Real rdOcp = R_d / Cp_d;  // R/cp = 0.286
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // Recompute theta from updated temperature
            // exner = (p/p0)^(R/cp)
            // theta = T / exner = T * (p0/p)^(R/cp)
            Real exner = std::pow(p_arr(i,j,k) / p0, rdOcp);
            theta_arr(i,j,k) = t_arr(i,j,k) / exner;
        });

        // (Tile-based precipitation diagnostics removed - using global diagnostics instead)

        // Accumulate precipitation
        ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            rain_arr(i,j,k) += rainacc_arr(i,j,k);
            snow_arr(i,j,k) += snowacc_arr(i,j,k);
            graup_arr(i,j,k) += graupacc_arr(i,j,k);
        });
        } else {
#endif
        // ===================================================================
        // WDM6 C++ GPU kernel path (adapted from WSM6 with double-moment)
        // ===================================================================

        // Working FABs (similar to WSM6 but with WDM6-specific additions)
        // 3D working arrays
        const Real dz_val = m_geom.CellSize(2);
        FArrayBox delz_fab(fab_box,1, Arena_Used);
        FArrayBox denfac_fab(fab_box,1, Arena_Used);
        FArrayBox xni_fab(fab_box,1, Arena_Used);
        FArrayBox rslopec_fab(fab_box,1, Arena_Used);
        FArrayBox rslopec2_fab(fab_box,1, Arena_Used);
        FArrayBox rslopec3_fab(fab_box,1, Arena_Used);
        FArrayBox rslope_fab(fab_box,3, Arena_Used);
        FArrayBox rslopeb_fab(fab_box,3, Arena_Used);
        FArrayBox rslope2_fab(fab_box,3, Arena_Used);
        FArrayBox rslope3_fab(fab_box,3, Arena_Used);
        FArrayBox work1_fab(fab_box,3, Arena_Used);
        FArrayBox workn_fab(fab_box,1, Arena_Used);
        FArrayBox work2_fab(fab_box,1, Arena_Used);  // Ventilation factor for diffusion (G11+)
        Box box2d(IntVect(ilo,jlo,0), IntVect(ihi,jhi,0));
        IArrayBox mstep_fab(box2d,1, Arena_Used);
        IArrayBox numdt_fab(box2d,1, Arena_Used);
        FArrayBox sr_fab(box2d, 1, Arena_Used);  // Snow ratio for G9 output
        FArrayBox cpm_fab(fab_box,1, Arena_Used);
        FArrayBox xl_fab(fab_box,1, Arena_Used);
        FArrayBox qsatw_fab(fab_box,1, Arena_Used);
        FArrayBox qsati_fab(fab_box,1, Arena_Used);
        FArrayBox rhw_fab(fab_box,1, Arena_Used);
        FArrayBox rhi_fab(fab_box,1, Arena_Used);
        FArrayBox qcr_fab(fab_box,1, Arena_Used);

        // Process rate arrays
        FArrayBox praut_fab(fab_box,1, Arena_Used);
        FArrayBox pracw_fab(fab_box,1, Arena_Used);
        FArrayBox prevp_fab(fab_box,1, Arena_Used);
        FArrayBox pcond_fab(fab_box,1, Arena_Used);

        // Number concentration process rates (WDM6-specific)
        FArrayBox ncauto_fab(fab_box,1, Arena_Used);  // nc lost to autoconversion
        FArrayBox ncaccr_fab(fab_box,1, Arena_Used);  // nc lost to accretion by rain
        FArrayBox nrauto_fab(fab_box,1, Arena_Used);  // nr gained from autoconversion
        FArrayBox nraccr_fab(fab_box,1, Arena_Used);  // nr gained from accretion
        FArrayBox nrevp_fab(fab_box,1, Arena_Used);   // nr lost to evaporation
        FArrayBox nccol_fab(fab_box,1, Arena_Used);   // cloud self-collection number sink
        FArrayBox nrcol_fab(fab_box,1, Arena_Used);   // rain self-collection number sink

        // G6 temporary arrays for slope_wdm6 call
        FArrayBox qrs_tmp_fab(fab_box,3, Arena_Used);  // Temporary qr, qs, qg (components 0,1,2)
        FArrayBox ncr_tmp_fab(fab_box,1, Arena_Used);  // Temporary nr

        // G11 arrays (particle diameter work array)
        FArrayBox avedia_fab(fab_box,2, Arena_Used);   // avedia(:,:,1:2) for rain and cloud slopes

        // G8 ice sedimentation arrays
        FArrayBox work1c_fab(fab_box,1, Arena_Used);   // Ice crystal fall speed
        FArrayBox fallc_fab(box2d,1, Arena_Used);      // Ice fallout (2D surface)
        FArrayBox delqi_fab(box2d,1, Arena_Used);      // Ice precipitation at bottom (2D)

        auto const& delz_arr = delz_fab.array();
        auto const& denfac_arr = denfac_fab.array();
        auto const& xni_arr = xni_fab.array();
        auto const& rslopec_arr = rslopec_fab.array();
        auto const& rslopec2_arr = rslopec2_fab.array();
        auto const& rslopec3_arr = rslopec3_fab.array();
        auto const& rslope_arr = rslope_fab.array();
        auto const& rslopeb_arr = rslopeb_fab.array();
        auto const& rslope2_arr = rslope2_fab.array();
        auto const& rslope3_arr = rslope3_fab.array();
        auto const& work1_arr = work1_fab.array();
        auto const& workn_arr = workn_fab.array();
        auto const& work2_arr = work2_fab.array();  // Ventilation factor for diffusion
        auto const& mstep_arr = mstep_fab.array();
        auto const& numdt_arr = numdt_fab.array();
        auto const& sr_arr = sr_fab.array();  // Snow ratio for G9
        auto const& cpm_arr = cpm_fab.array();
        auto const& xl_arr = xl_fab.array();
        auto const& qsatw_arr = qsatw_fab.array();
        auto const& qsati_arr = qsati_fab.array();
        auto const& rhw_arr = rhw_fab.array();
        auto const& rhi_arr = rhi_fab.array();
        auto const& qcr_arr = qcr_fab.array();
        auto const& praut_arr = praut_fab.array();
        auto const& pracw_arr = pracw_fab.array();
        auto const& prevp_arr = prevp_fab.array();
        auto const& pcond_arr = pcond_fab.array();
        auto const& ncauto_arr = ncauto_fab.array();
        auto const& ncaccr_arr = ncaccr_fab.array();
        auto const& nrauto_arr = nrauto_fab.array();
        auto const& nraccr_arr = nraccr_fab.array();
        auto const& nrevp_arr = nrevp_fab.array();
        auto const& nccol_arr = nccol_fab.array();
        auto const& nrcol_arr = nrcol_fab.array();
        auto const& qrs_tmp_arr = qrs_tmp_fab.array();  // G6: temporary qr, qs, qg
        auto const& ncr_tmp_arr = ncr_tmp_fab.array();  // G6: temporary nr
        auto const& avedia_arr = avedia_fab.array();    // G11: particle diameter work array

        // G8 array references
        auto const& work1c_arr = work1c_fab.array();    // Ice crystal fall speed
        auto const& fallc_arr = fallc_fab.array();      // Ice fallout (2D)
        auto const& delqi_arr = delqi_fab.array();      // Ice precip at bottom (2D)

        // Clamp negative values and enforce minimums
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            delz_arr(i,j,k) = dz_val;
            qc_arr(i,j,k) = amrex::max(qc_arr(i,j,k), Real(0.0));
            qr_arr(i,j,k) = amrex::max(qr_arr(i,j,k), Real(0.0));
            qi_arr(i,j,k) = amrex::max(qi_arr(i,j,k), Real(0.0));
            qs_arr(i,j,k) = amrex::max(qs_arr(i,j,k), Real(0.0));
            qg_arr(i,j,k) = amrex::max(qg_arr(i,j,k), Real(0.0));

            // Match Fortran pre-G3 behavior: nc is non-negative here, but not floored to ncmin yet.
            nc_arr(i,j,k) = amrex::max(nc_arr(i,j,k), Real(0.0));
            nr_arr(i,j,k) = amrex::max(nr_arr(i,j,k), Real(0.0));
            nn_arr(i,j,k) = amrex::max(nn_arr(i,j,k), Real(0.0));
        });

        // Compute cpm and xl once from initial state
        const Real xlv1_loc = m_xlv1;
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            cpm_arr(i,j,k) = wdm6_cpmcal(qv_arr(i,j,k), Real(qmin), Real(cpd), Real(cpv));
            xl_arr(i,j,k) = wdm6_xlcal(t_arr(i,j,k), Real(xlv0), xlv1_loc, Real(t0c));
        });

        // Minor timestep loop (match WSM6 structure)
        const int wdm6_loops = std::max(
            static_cast<int>(std::round(dt_advance / Real(120.0))), 1);  // dtcldcr = 120s
        const Real dtcld = dt_advance / static_cast<Real>(wdm6_loops);

        // Extract WDM6 coefficients
        const Real qc0_loc = m_qc0;
        const Real qc1_loc = m_qc1;
        const Real qck1_loc = m_qck1;
        const Real pidnc_loc = m_pidnc;
        const Real pidnr_loc = m_pidnr;
        const Real pvtr_loc = m_pvtr;
        const Real pvtrn_loc = m_pvtrn;
        const Real pvts_loc = m_pvts;
        const Real pvtg_loc = m_pvtg;
        const Real slope_bvtg_loc = m_bvtg;
        const Real pidn0s_loc = m_pidn0s;
        const Real pidn0g_loc = m_pidn0g;
        const Real rslopermax_loc = m_rslopermax;
        const Real rsloperbmax_loc = m_rsloperbmax;
        const Real rsloper2max_loc = m_rsloper2max;
        const Real rsloper3max_loc = m_rsloper3max;
        const Real rslopesmax_loc = m_rslopesmax;
        const Real rslopesbmax_loc = m_rslopesbmax;
        const Real rslopes2max_loc = m_rslopes2max;
        const Real rslopes3max_loc = m_rslopes3max;
        const Real rslopegmax_loc = m_rslopegmax;
        const Real rslopegbmax_loc = m_rslopegbmax;
        const Real rslopeg2max_loc = m_rslopeg2max;
        const Real rslopeg3max_loc = m_rslopeg3max;
        const Real rslopecmax_loc = m_rslopecmax;
        const Real rslopec2max_loc = m_rslopec2max;
        const Real rslopec3max_loc = m_rslopec3max;
        const Real precs1_loc = m_precs1;
        const Real precs2_loc = m_precs2;
        const Real precg1_loc = m_precg1;
        const Real precg2_loc = m_precg2;
        const Real precr1_loc = m_precr1;
        const Real precr2_loc = m_precr2;
        const Real n0g_loc = m_n0g;
        const Real pi_wdm6_loc = m_pi_wdm6;
        constexpr Real pfrz1_loc = pfrz1;
        constexpr Real pfrz2_loc = pfrz2;
        const bool diag_col_in_tile = (diag_i >= ilo && diag_i <= ihi &&
                                       diag_j >= jlo && diag_j <= jhi);
        const int diag_k = klo;

        if (m_lmask != nullptr) {
            auto const& lmask_arr = m_lmask->const_array(mfi);
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                qcr_arr(i,j,k) = (lmask_arr(i,j,0) == 2) ? qc0_loc : qc1_loc;
            });
        } else {
            // Match the current bridge posture: default to land when no lmask is wired.
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                qcr_arr(i,j,k) = qc1_loc;
            });
        }

        for (int loop = 0; loop < wdm6_loops; ++loop) {
            // ============================================================
            // Step 1: Density factor (G1b / DENFAC)
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G1B %3d %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(p_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                denfac_arr(i,j,k) = std::sqrt(Real(den0) / den_arr(i,j,k));
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G1B %3d %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(denfac_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 2: Saturation calculations (G1c / QSAT)
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G1C %3d %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(p_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nn_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif
            {
                const Real ttp = Real(t0c) + Real(0.01f);
                const Real dldt = Real(cpv) - Real(cliq);
                const Real xa = -dldt / Real(rv);
                const Real xb = xa + Real(xlv0) / (Real(rv) * ttp);
                const Real dldti = Real(cpv) - Real(cice);
                const Real xai = -dldti / Real(rv);
                const Real xbi = xai + Real(xls) / (Real(rv) * ttp);

#if !defined(AMREX_USE_GPU)
                if (microphysics_debug > 0 && diag_col_in_tile) {
                    std::printf("WDM6-CPP_G1CV_SAT_INT0 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(ttp),
                                static_cast<double>(dldt),
                                static_cast<double>(xa),
                                static_cast<double>(xb),
                                static_cast<double>(dldti),
                                static_cast<double>(xai),
                                static_cast<double>(xbi));
                    std::fflush(stdout);
                }
#endif

                ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    const Real tr = ttp / t_arr(i,j,k);

                    // Saturation over water
#if !defined(AMREX_USE_GPU)
                    if (i == diag_i && j == diag_j && k == diag_k) {
                        std::printf("WDM6-CPP_G1CV_SAT_INT1 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(t_arr(i,j,k)),
                                    static_cast<double>(p_arr(i,j,k)),
                                    static_cast<double>(qv_arr(i,j,k)),
                                    static_cast<double>(Real(ep2)),
                                    static_cast<double>(Real(psat)),
                                    static_cast<double>(Real(qmin)));
                        std::printf("WDM6-CPP_G1CV_SAT_INT2 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(tr));
                        std::printf("WDM6-CPP_G1CV_SAT_INT3 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(std::log(tr)));
                        std::printf("WDM6-CPP_G1CV_SAT_INT4 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(std::log(tr) * xa));
                        std::printf("WDM6-CPP_G1CV_SAT_INT5 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(std::exp(std::log(tr) * xa)));
                        std::printf("WDM6-CPP_G1CV_SAT_INT6 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(Real(1.0) - tr));
                        std::printf("WDM6-CPP_G1CV_SAT_INT7 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(xb * (Real(1.0) - tr)));
                        std::printf("WDM6-CPP_G1CV_SAT_INT8 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(std::exp(xb * (Real(1.0) - tr))));
                        std::printf("WDM6-CPP_G1CV_SAT_INT9 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(Real(psat) * std::exp(std::log(tr) * xa)
                                                       * std::exp(xb * (Real(1.0) - tr))));
                    }
#endif
                    Real qsw = Real(psat) * std::exp(std::log(tr) * xa) * std::exp(xb * (Real(1.0) - tr));
                    qsw = amrex::min(qsw, Real(0.99) * p_arr(i,j,k));
#if !defined(AMREX_USE_GPU)
                    if (i == diag_i && j == diag_j && k == diag_k) {
                        std::printf("WDM6-CPP_G1CV_SAT_INT10 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(qsw));
                    }
#endif
                    qsw = Real(ep2) * qsw / (p_arr(i,j,k) - qsw);
#if !defined(AMREX_USE_GPU)
                    if (i == diag_i && j == diag_j && k == diag_k) {
                        std::printf("WDM6-CPP_G1CV_SAT_INT11 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(qsw));
                    }
#endif
                    qsw = amrex::max(qsw, Real(qmin));
#if !defined(AMREX_USE_GPU)
                    if (i == diag_i && j == diag_j && k == diag_k) {
                        std::printf("WDM6-CPP_G1CV_SAT_INT12 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(qsw));
                    }
#endif
                    qsatw_arr(i,j,k) = qsw;
                    rhw_arr(i,j,k) = amrex::max(qv_arr(i,j,k) / qsw, Real(qmin));
#if !defined(AMREX_USE_GPU)
                    if (i == diag_i && j == diag_j && k == diag_k) {
                        std::printf("WDM6-CPP_G1CV_SAT_INT13 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(rhw_arr(i,j,k)));
                    }
#endif

                    // Saturation over ice
#if !defined(AMREX_USE_GPU)
                    if (i == diag_i && j == diag_j && k == diag_k) {
                        std::printf("WDM6-CPP_G1CV_SAT_INT20 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(tr));
                        std::printf("WDM6-CPP_G1CV_SAT_INT21 %3d %24.16E %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>((t_arr(i,j,k) < ttp) ? Real(1.0) : Real(0.0)),
                                    static_cast<double>((t_arr(i,j,k) < ttp)
                                        ? Real(psat) * std::exp(std::log(tr) * xai) * std::exp(xbi * (Real(1.0) - tr))
                                        : Real(psat) * std::exp(std::log(tr) * xa) * std::exp(xb * (Real(1.0) - tr))));
                    }
#endif
                    Real qsi = (t_arr(i,j,k) < ttp)
                        ? Real(psat) * std::exp(std::log(tr) * xai) * std::exp(xbi * (Real(1.0) - tr))
                        : Real(psat) * std::exp(std::log(tr) * xa) * std::exp(xb * (Real(1.0) - tr));
                    qsi = amrex::min(qsi, Real(0.99) * p_arr(i,j,k));
#if !defined(AMREX_USE_GPU)
                    if (i == diag_i && j == diag_j && k == diag_k) {
                        std::printf("WDM6-CPP_G1CV_SAT_INT22 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(qsi));
                    }
#endif
                    qsi = Real(ep2) * qsi / (p_arr(i,j,k) - qsi);
#if !defined(AMREX_USE_GPU)
                    if (i == diag_i && j == diag_j && k == diag_k) {
                        std::printf("WDM6-CPP_G1CV_SAT_INT23 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(qsi));
                    }
#endif
                    qsi = amrex::max(qsi, Real(qmin));
#if !defined(AMREX_USE_GPU)
                    if (i == diag_i && j == diag_j && k == diag_k) {
                        std::printf("WDM6-CPP_G1CV_SAT_INT24 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(qsi));
                    }
#endif
                    qsati_arr(i,j,k) = qsi;
                    rhi_arr(i,j,k) = amrex::max(qv_arr(i,j,k) / qsi, Real(qmin));
#if !defined(AMREX_USE_GPU)
                    if (i == diag_i && j == diag_j && k == diag_k) {
                        std::printf("WDM6-CPP_G1CV_SAT_INT25 %3d %24.16E\n",
                                    diag_k + 1,
                                    static_cast<double>(rhi_arr(i,j,k)));
                        std::fflush(stdout);
                    }
#endif
                });
            }

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G1C %3d %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qsatw_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qsati_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rhw_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rhi_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3: Zero process rates (G2 / RATES_ZERO)
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G2 %3d\n", diag_k + 1);
                std::fflush(stdout);
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                praut_arr(i,j,k) = Real(0.0);
                pracw_arr(i,j,k) = Real(0.0);
                prevp_arr(i,j,k) = Real(0.0);
                pcond_arr(i,j,k) = Real(0.0);
                ncauto_arr(i,j,k) = Real(0.0);
                ncaccr_arr(i,j,k) = Real(0.0);
                nrauto_arr(i,j,k) = Real(0.0);
                nraccr_arr(i,j,k) = Real(0.0);
                nrevp_arr(i,j,k) = Real(0.0);
                nccol_arr(i,j,k) = Real(0.0);
                nrcol_arr(i,j,k) = Real(0.0);
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G2 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(praut_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(pracw_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(prevp_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(pcond_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(ncauto_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(ncaccr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nrauto_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nraccr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nrevp_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3b: CLOUD_SETUP (G3)
            // Exact port of the frozen Fortran block:
            //   - cloud droplet slope parameter rslopec{,2,3}
            //   - ice number concentration xni
            // ============================================================
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (qc_arr(i,j,k) <= Real(qmin) || nc_arr(i,j,k) <= Real(1.e1)) {
                    rslopec_arr(i,j,k)  = rslopecmax_loc;
                    rslopec2_arr(i,j,k) = rslopec2max_loc;
                    rslopec3_arr(i,j,k) = rslopec3max_loc;
                } else {
                    rslopec_arr(i,j,k) = wdm6_rslopec_exact(qc_arr(i,j,k), den_arr(i,j,k),
                                                            nc_arr(i,j,k), pidnc_loc);
                    rslopec2_arr(i,j,k) = rslopec_arr(i,j,k) * rslopec_arr(i,j,k);
                    rslopec3_arr(i,j,k) = rslopec2_arr(i,j,k) * rslopec_arr(i,j,k);
                }
                xni_arr(i,j,k) = wdm6_xni_exact(qi_arr(i,j,k), den_arr(i,j,k), Real(qmin));
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G3 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(rslopec_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslopec2_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslopec3_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(xni_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qi_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G3 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(rslopec_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslopec2_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslopec3_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(xni_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qi_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3c: SLOPE1 (G4)
            // Exact first packed slope_wdm6 surface for rain/snow/graupel.
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G4 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qs_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qg_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real rain_rslope, rain_rslopeb, rain_rslope2, rain_rslope3, rain_vt, rain_vtn;
                Real snow_rslope, snow_rslopeb, snow_rslope2, snow_rslope3, snow_vt, snow_n0sfac;
                Real graup_rslope, graup_rslopeb, graup_rslope2, graup_rslope3, graup_vt;

                wdm6_slope_rain_cell(qr_arr(i,j,k), nr_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k),
                                     Real(qcrmin), Real(nrmin),
                                     rslopermax_loc, rsloperbmax_loc, rsloper2max_loc, rsloper3max_loc,
                                     Real(bvtr), pvtr_loc, pvtrn_loc, pidnr_loc,
                                     rain_rslope, rain_rslopeb, rain_rslope2, rain_rslope3,
                                     rain_vt, rain_vtn);
                wdm6_slope_snow_cell(qs_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k), t_arr(i,j,k),
                                     pidn0s_loc, Real(alpha_wdm6), Real(n0smax), Real(n0s),
                                     Real(t0c), Real(qcrmin),
                                     rslopesmax_loc, rslopesbmax_loc, rslopes2max_loc, rslopes3max_loc,
                                     Real(bvts), pvts_loc,
                                     snow_rslope, snow_rslopeb, snow_rslope2, snow_rslope3, snow_vt,
                                     snow_n0sfac);
                wdm6_slope_graup_cell(qg_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k),
                                      pidn0g_loc, Real(qcrmin),
                                      rslopegmax_loc, rslopegbmax_loc, rslopeg2max_loc, rslopeg3max_loc,
                                      slope_bvtg_loc, pvtg_loc,
                                      graup_rslope, graup_rslopeb, graup_rslope2, graup_rslope3, graup_vt);

                rslope_arr(i,j,k,0) = rain_rslope;
                rslope_arr(i,j,k,1) = snow_rslope;
                rslope_arr(i,j,k,2) = graup_rslope;
                rslopeb_arr(i,j,k,0) = rain_rslopeb;
                rslopeb_arr(i,j,k,1) = snow_rslopeb;
                rslopeb_arr(i,j,k,2) = graup_rslopeb;
                rslope2_arr(i,j,k,0) = rain_rslope2;
                rslope2_arr(i,j,k,1) = snow_rslope2;
                rslope2_arr(i,j,k,2) = graup_rslope2;
                rslope3_arr(i,j,k,0) = rain_rslope3;
                rslope3_arr(i,j,k,1) = snow_rslope3;
                rslope3_arr(i,j,k,2) = graup_rslope3;
                work1_arr(i,j,k,0) = rain_vt;
                work1_arr(i,j,k,1) = snow_vt;
                work1_arr(i,j,k,2) = graup_vt;
                workn_arr(i,j,k) = rain_vtn;
            });
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G4 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(rslope_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(workn_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qs_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qg_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,diag_k,2)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3d: Rain sedimentation setup (G5a)
            // Determine rain/number substep count from normalized work1/workn.
            // ============================================================
            ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                mstep_arr(i,j,k) = 1;
                numdt_arr(i,j,k) = 1;
                sr_arr(i,j,k) = Real(0.0);  // Initialize snow ratio to zero
            });
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G5A %3d %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(workn_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(delz_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(dtcld));
                std::fflush(stdout);
            }
#endif
            ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                int mstep_loc = 1;
                int numdt_loc = 1;
                for (int kk = khi; kk >= klo; --kk) {
                    work1_arr(i,j,kk,0) = work1_arr(i,j,kk,0) / delz_arr(i,j,kk);
                    workn_arr(i,j,kk) = workn_arr(i,j,kk) / delz_arr(i,j,kk);
                    numdt_loc = amrex::max(static_cast<int>(amrex::max(work1_arr(i,j,kk,0),
                                                                       workn_arr(i,j,kk)) * dtcld + Real(0.5)),
                                           1);
                    if (numdt_loc >= mstep_loc) {
                        mstep_loc = numdt_loc;
                    }
                }
                mstep_arr(i,j,k) = mstep_loc;
                numdt_arr(i,j,k) = numdt_loc;
            });
            ReduceOps<ReduceOpMax> reduce_op;
            ReduceData<int> reduce_data(reduce_op);
            reduce_op.eval(box2d, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<int> {
                return {mstep_arr(i,j,k)};
            });
            int mstepmax = amrex::get<0>(reduce_data.value());
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G5A %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(workn_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(mstep_arr(diag_i,diag_j,0)),
                            static_cast<double>(mstepmax),
                            static_cast<double>(numdt_arr(diag_i,diag_j,0)),
                            static_cast<double>(delz_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(dtcld));
                std::fflush(stdout);
            }
#else
            amrex::ignore_unused(mstepmax);
#endif

            // ============================================================
            // Step 3e: Rain sedimentation substeps (G5b)
            // Match the bounded top-down rain/nr transport loop and
            // slope_rain refresh that immediately follows G5a in Fortran.
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G5B %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(workn_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(mstep_arr(diag_i,diag_j,0)),
                            static_cast<double>(delz_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif
            ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                amrex::ignore_unused(k);
                const int col_mstep = mstep_arr(i,j,0);
                for (int n = 1; n <= mstepmax; ++n) {
                    if (n > col_mstep) {
                        continue;
                    }

                    const int kk_top = khi;
                    const Real top_flux_qr = den_arr(i,j,kk_top) * qr_arr(i,j,kk_top)
                                           * work1_arr(i,j,kk_top,0) / static_cast<Real>(col_mstep);
                    const Real top_flux_nr = nr_arr(i,j,kk_top)
                                           * workn_arr(i,j,kk_top) / static_cast<Real>(col_mstep);

                    qr_arr(i,j,kk_top) = amrex::max(
                        qr_arr(i,j,kk_top) - top_flux_qr * dtcld / den_arr(i,j,kk_top),
                        Real(0.0));
                    nr_arr(i,j,kk_top) = amrex::max(
                        nr_arr(i,j,kk_top) - top_flux_nr * dtcld,
                        Real(0.0));

                    Real flux_qr_above = top_flux_qr;
                    Real flux_nr_above = top_flux_nr;
                    for (int kk = khi - 1; kk >= klo; --kk) {
                        const Real flux_qr = den_arr(i,j,kk) * qr_arr(i,j,kk)
                                           * work1_arr(i,j,kk,0) / static_cast<Real>(col_mstep);
                        const Real flux_nr = nr_arr(i,j,kk)
                                           * workn_arr(i,j,kk) / static_cast<Real>(col_mstep);

                        const Real dqr_self = amrex::min(
                            flux_qr * dtcld / den_arr(i,j,kk),
                            qr_arr(i,j,kk));
                        const Real dqr_from_above = amrex::min(
                            flux_qr_above * delz_arr(i,j,kk+1) / delz_arr(i,j,kk)
                                * dtcld / den_arr(i,j,kk),
                            qr_arr(i,j,kk+1));
                        const Real dnr_self = amrex::min(
                            flux_nr * dtcld,
                            nr_arr(i,j,kk));
                        const Real dnr_from_above = amrex::min(
                            flux_nr_above * delz_arr(i,j,kk+1) / delz_arr(i,j,kk)
                                * dtcld,
                            nr_arr(i,j,kk+1));

                        qr_arr(i,j,kk) = amrex::max(
                            qr_arr(i,j,kk) - dqr_self + dqr_from_above,
                            Real(0.0));
                        nr_arr(i,j,kk) = amrex::max(
                            nr_arr(i,j,kk) - dnr_self + dnr_from_above,
                            Real(0.0));

                        flux_qr_above = flux_qr;
                        flux_nr_above = flux_nr;
                    }

                    for (int kk = klo; kk <= khi; ++kk) {
                        Real rain_rslope, rain_rslopeb, rain_rslope2, rain_rslope3;
                        Real rain_vt, rain_vtn;
                        wdm6_slope_rain_cell(
                            qr_arr(i,j,kk), nr_arr(i,j,kk),
                            den_arr(i,j,kk), denfac_arr(i,j,kk),
                            Real(qcrmin), Real(nrmin),
                            rslopermax_loc, rsloperbmax_loc,
                            rsloper2max_loc, rsloper3max_loc,
                            Real(bvtr), pvtr_loc, pvtrn_loc, pidnr_loc,
                            rain_rslope, rain_rslopeb, rain_rslope2, rain_rslope3,
                            rain_vt, rain_vtn);
                        rslope_arr(i,j,kk,0) = rain_rslope;
                        rslopeb_arr(i,j,kk,0) = rain_rslopeb;
                        rslope2_arr(i,j,kk,0) = rain_rslope2;
                        rslope3_arr(i,j,kk,0) = rain_rslope3;
                        work1_arr(i,j,kk,0) = rain_vt / delz_arr(i,j,kk);
                        workn_arr(i,j,kk) = rain_vtn / delz_arr(i,j,kk);
                    }
                }
            });
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G5B %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(workn_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(delz_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(dtcld));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3f: Snow/graupel sedimentation (G5c)
            // Exact bounded NISLFV_SG slice plus lower-boundary slab fall.
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                const Real qsum = amrex::max(qs_arr(diag_i,diag_j,diag_k) + qg_arr(diag_i,diag_j,diag_k), Real(1.0e-15));
                const Real worka = (qsum > Real(1.0e-15))
                    ? (work1_arr(diag_i,diag_j,diag_k,1) * qs_arr(diag_i,diag_j,diag_k)
                     + work1_arr(diag_i,diag_j,diag_k,2) * qg_arr(diag_i,diag_j,diag_k)) / qsum
                    : Real(0.0);
                std::printf("WDM6-CPP_PRE_G5C %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qs_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qg_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(worka),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k) * qs_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k) * qg_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(delz_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(dtcld),
                            static_cast<double>(denfac_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }

            ParallelFor(box2d, [=] (int i, int j, int k) {
                amrex::ignore_unused(k);
                const int km = khi - klo + 1;
                std::vector<Real> dz(km), den(km), denfac(km), tk(km), worka(km), qs_col(km), qg_col(km);
                for (int kk = 0; kk < km; ++kk) {
                    const int k3 = klo + kk;
                    dz[kk] = delz_arr(i,j,k3);
                    den[kk] = den_arr(i,j,k3);
                    denfac[kk] = denfac_arr(i,j,k3);
                    tk[kk] = t_arr(i,j,k3);
                    qs_col[kk] = den[kk] * qs_arr(i,j,k3);
                    qg_col[kk] = den[kk] * qg_arr(i,j,k3);
                    const Real qsum = amrex::max(qs_arr(i,j,k3) + qg_arr(i,j,k3), Real(1.0e-15));
                    worka[kk] = (qsum > Real(1.0e-15))
                        ? (work1_arr(i,j,k3,1) * qs_arr(i,j,k3) + work1_arr(i,j,k3,2) * qg_arr(i,j,k3)) / qsum
                        : Real(0.0);
                }

                Real delqrs2 = Real(0.0);
                Real delqrs3 = Real(0.0);
                wdm6_nislfv_rain_plm6_column(
                    km, dz.data(), den.data(), denfac.data(), tk.data(),
                    worka.data(), qs_col.data(), qg_col.data(), delqrs2, delqrs3, dtcld, 1,
                    pidn0s_loc, pidn0g_loc, Real(qcrmin), Real(alpha_wdm6), Real(n0smax), Real(n0s), Real(t0c),
                    rslopesmax_loc, rslopesbmax_loc, rslopes2max_loc, rslopes3max_loc, Real(bvts), pvts_loc,
                    rslopegmax_loc, rslopegbmax_loc, rslopeg2max_loc, rslopeg3max_loc, slope_bvtg_loc, pvtg_loc);

                for (int kk = 0; kk < km; ++kk) {
                    const int k3 = klo + kk;
                    qs_arr(i,j,k3) = amrex::max(qs_col[kk] / den[kk], Real(0.0));
                    qg_arr(i,j,k3) = amrex::max(qg_col[kk] / den[kk], Real(0.0));
                }

                work1_arr(i,j,klo,1) = delqrs2 / dz[0] / dtcld;
                work1_arr(i,j,klo,2) = delqrs3 / dz[0] / dtcld;
            });

            if (microphysics_debug > 0 && diag_col_in_tile) {
                const Real qsum = amrex::max(qs_arr(diag_i,diag_j,diag_k) + qg_arr(diag_i,diag_j,diag_k), Real(1.0e-15));
                const Real worka = (qsum > Real(1.0e-15))
                    ? (work1_arr(diag_i,diag_j,diag_k,1) * qs_arr(diag_i,diag_j,diag_k)
                     + work1_arr(diag_i,diag_j,diag_k,2) * qg_arr(diag_i,diag_j,diag_k)) / qsum
                    : Real(0.0);
                std::printf("WDM6-CPP_POST_G5C %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qs_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qg_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(worka),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k) * qs_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k) * qg_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(work1_arr(diag_i,diag_j,klo,1)),
                            static_cast<double>(work1_arr(diag_i,diag_j,klo,2)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3g: Second slope_wdm6 call after sedimentation (G6)
            // Recompute slope parameters for all species after G5c sed.
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G6 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qs_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qg_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // Copy current state to temporary arrays (qrs_tmp, ncr_tmp)
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                qrs_tmp_arr(i,j,k,0) = qr_arr(i,j,k);  // Fortran qrs(:,:,1)
                qrs_tmp_arr(i,j,k,1) = qs_arr(i,j,k);  // Fortran qrs(:,:,2)
                qrs_tmp_arr(i,j,k,2) = qg_arr(i,j,k);  // Fortran qrs(:,:,3)
                ncr_tmp_arr(i,j,k)   = nr_arr(i,j,k);  // Fortran ncr(:,:,3)
            });

            // Call slope_wdm6: compute slope parameters for rain/snow/graupel
            // This parallels the Fortran call: slope_wdm6(qrs_tmp,ncr_tmp,...)
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real rain_rslope, rain_rslopeb, rain_rslope2, rain_rslope3, rain_vt, rain_vtn;
                Real snow_rslope, snow_rslopeb, snow_rslope2, snow_rslope3, snow_vt, snow_n0sfac;
                Real graup_rslope, graup_rslopeb, graup_rslope2, graup_rslope3, graup_vt;

                // Compute slope parameters using temporary (post-sedimentation) values
                wdm6_slope_rain_cell(qrs_tmp_arr(i,j,k,0), ncr_tmp_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k),
                                     Real(qcrmin), Real(nrmin),
                                     rslopermax_loc, rsloperbmax_loc, rsloper2max_loc, rsloper3max_loc,
                                     Real(bvtr), pvtr_loc, pvtrn_loc, pidnr_loc,
                                     rain_rslope, rain_rslopeb, rain_rslope2, rain_rslope3,
                                     rain_vt, rain_vtn);
                wdm6_slope_snow_cell(qrs_tmp_arr(i,j,k,1), den_arr(i,j,k), denfac_arr(i,j,k), t_arr(i,j,k),
                                     pidn0s_loc, Real(alpha_wdm6), Real(n0smax), Real(n0s),
                                     Real(t0c), Real(qcrmin),
                                     rslopesmax_loc, rslopesbmax_loc, rslopes2max_loc, rslopes3max_loc,
                                     Real(bvts), pvts_loc,
                                     snow_rslope, snow_rslopeb, snow_rslope2, snow_rslope3, snow_vt,
                                     snow_n0sfac);
                wdm6_slope_graup_cell(qrs_tmp_arr(i,j,k,2), den_arr(i,j,k), denfac_arr(i,j,k),
                                      pidn0g_loc, Real(qcrmin),
                                      rslopegmax_loc, rslopegbmax_loc, rslopeg2max_loc, rslopeg3max_loc,
                                      slope_bvtg_loc, pvtg_loc,
                                      graup_rslope, graup_rslopeb, graup_rslope2, graup_rslope3, graup_vt);

                // Store results in output slope arrays
                rslope_arr(i,j,k,0) = rain_rslope;
                rslope_arr(i,j,k,1) = snow_rslope;
                rslope_arr(i,j,k,2) = graup_rslope;
                rslopeb_arr(i,j,k,0) = rain_rslopeb;
                rslopeb_arr(i,j,k,1) = snow_rslopeb;
                rslopeb_arr(i,j,k,2) = graup_rslopeb;
                rslope2_arr(i,j,k,0) = rain_rslope2;
                rslope2_arr(i,j,k,1) = snow_rslope2;
                rslope2_arr(i,j,k,2) = graup_rslope2;
                rslope3_arr(i,j,k,0) = rain_rslope3;
                rslope3_arr(i,j,k,1) = snow_rslope3;
                rslope3_arr(i,j,k,2) = graup_rslope3;
                work1_arr(i,j,k,0) = rain_vt;
                work1_arr(i,j,k,1) = snow_vt;
                work1_arr(i,j,k,2) = graup_vt;
                workn_arr(i,j,k) = rain_vtn;
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G6 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(rslope_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(workn_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // WDM6-CPP_PRE_G7
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G7 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qs_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qg_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(p_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(cpm_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,diag_k,2)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3h: G7 warm-phase snow/graupel melting
            // ============================================================
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Warm-phase snow and graupel melting (t > t0c)
                if (t_arr(i,j,k) > t0c) {
                    const Real supcol = t0c - t_arr(i,j,k);
                    const Real n0sfac = amrex::max(
                        amrex::min(std::exp(Real(alpha_wdm6) * supcol),
                                   Real(n0smax) / Real(n0s)),
                        Real(1.0));

                    const Real xlf = xlf0;
                    const Real work2 = wdm6_venfac(p_arr(i,j,k), t_arr(i,j,k),
                                                    den_arr(i,j,k), Real(den0));

                    // --- Snow melting: psmlt ---
                    if (qs_arr(i,j,k) > Real(0.0)) {
                        const Real coeres_s = rslope2_arr(i,j,k,1) *
                            std::sqrt(rslope_arr(i,j,k,1) * rslopeb_arr(i,j,k,1));

                        Real psmlt = wdm6_xka(t_arr(i,j,k), den_arr(i,j,k)) / xlf *
                            (t0c - t_arr(i,j,k)) * pi_wdm6_loc * Real(0.5) * n0sfac *
                            (precs1_loc * rslope2_arr(i,j,k,1) +
                             precs2_loc * work2 * coeres_s) / den_arr(i,j,k);

                        psmlt = amrex::min(amrex::max(psmlt * dtcld, -qs_arr(i,j,k)),
                                           Real(0.0));

                        if (qs_arr(i,j,k) > Real(qcrmin)) {
                            const Real sfac = rslope_arr(i,j,k,1) * Real(n0s) * n0sfac /
                                              qs_arr(i,j,k);
                            nr_arr(i,j,k) = amrex::max(nr_arr(i,j,k) - sfac * psmlt, Real(0.0));
                        }

                        qs_arr(i,j,k) += psmlt;
                        qr_arr(i,j,k) -= psmlt;
                        t_arr(i,j,k) += xlf / cpm_arr(i,j,k) * psmlt;
                    }

                    // --- Graupel melting: pgmlt ---
                    if (qg_arr(i,j,k) > Real(0.0)) {
                        const Real coeres_g = rslope2_arr(i,j,k,2) *
                            std::sqrt(rslope_arr(i,j,k,2) * rslopeb_arr(i,j,k,2));

                        Real pgmlt = wdm6_xka(t_arr(i,j,k), den_arr(i,j,k)) / xlf *
                            (t0c - t_arr(i,j,k)) *
                            (precg1_loc * rslope2_arr(i,j,k,2) +
                             precg2_loc * work2 * coeres_g) / den_arr(i,j,k);

                        pgmlt = amrex::min(amrex::max(pgmlt * dtcld, -qg_arr(i,j,k)),
                                           Real(0.0));

                        if (qg_arr(i,j,k) > Real(qcrmin)) {
                            const Real gfac = rslope_arr(i,j,k,2) * n0g_loc /
                                              qg_arr(i,j,k);
                            nr_arr(i,j,k) = amrex::max(nr_arr(i,j,k) - gfac * pgmlt, Real(0.0));
                        }

                        qg_arr(i,j,k) += pgmlt;
                        qr_arr(i,j,k) -= pgmlt;
                        t_arr(i,j,k) += xlf / cpm_arr(i,j,k) * pgmlt;
                    }
                }
            });

            // ============================================================
            // WDM6-CPP_POST_G7
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G7 %3d %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qs_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qg_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3i: G8 ice fall speed + ice sedimentation
            // VICE: Ice crystal terminal velocity + nislfv_rain_plmr
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G8 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qi_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(xni_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(delz_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // Compute ice crystal fall speed (work1c)
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real work1c = Real(0.0);
                if (qi_arr(i,j,k) > Real(0.0)) {
                    const Real xni_safe = amrex::max(xni_arr(i,j,k), Real(1.0e-30));
                    const Real xmi = den_arr(i,j,k) * qi_arr(i,j,k) / xni_safe;
                    constexpr Real dicon = Real(11.45e-9);  // Ice diameter coefficient (m kg^-1/3)
                    constexpr Real dimax = Real(500.e-6);   // Max ice diameter (m)
                    const Real diameter = amrex::max(amrex::min(dicon * std::sqrt(xmi), dimax), Real(1.0e-25));
                    work1c = Real(1.49e4) * std::exp(std::log(diameter) * Real(1.31));
                }
                work1c_arr(i,j,k) = work1c;
            });

            // Ice sedimentation via simplified PLM6 scheme
            ParallelFor(box2d, [=] (int i, int j, int k) {
                amrex::ignore_unused(k);
                const int km = khi - klo + 1;
                std::vector<Real> dz(km), den(km), denfac(km), tk(km), work_ice(km);
                std::vector<Real> qi_col(km);

                // Fill column arrays
                for (int kk = 0; kk < km; ++kk) {
                    const int k3 = klo + kk;
                    dz[kk] = delz_arr(i,j,k3);
                    den[kk] = den_arr(i,j,k3);
                    denfac[kk] = denfac_arr(i,j,k3);
                    tk[kk] = t_arr(i,j,k3);
                    work_ice[kk] = work1c_arr(i,j,k3);
                    qi_col[kk] = den[kk] * qi_arr(i,j,k3);
                }

                Real delqi_col = Real(0.0);

                // Simple sedimentation: top-down fall with mass conservation
                wdm6_nislfv_rain_plm6_column(
                    km, dz.data(), den.data(), denfac.data(), tk.data(),
                    work_ice.data(), qi_col.data(), qi_col.data(), delqi_col, delqi_col, dtcld, 1,
                    pidn0s_loc, pidn0g_loc, Real(qcrmin), Real(alpha_wdm6), Real(n0smax), Real(n0s), Real(t0c),
                    rslopesmax_loc, rslopesbmax_loc, rslopes2max_loc, rslopes3max_loc, Real(bvts), pvts_loc,
                    rslopegmax_loc, rslopegbmax_loc, rslopeg2max_loc, rslopeg3max_loc, slope_bvtg_loc, pvtg_loc);

                // Update ice concentrations from sedimented state
                for (int kk = 0; kk < km; ++kk) {
                    const int k3 = klo + kk;
                    qi_arr(i,j,k3) = amrex::max(qi_col[kk] / den[kk], Real(0.0));
                }

                // Ice fallout at surface
                fallc_arr(i,j,k) = delqi_col / dz[0] / dtcld;
                delqi_arr(i,j,k) = delqi_col;
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G8 %3d %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qi_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(fallc_arr(diag_i,diag_j,klo)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k) * qi_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3h: G9 Surface precipitation accumulation and sr update
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G9 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(work1_arr(diag_i,diag_j,klo,0)),
                            static_cast<double>(work1_arr(diag_i,diag_j,klo,1)),
                            static_cast<double>(work1_arr(diag_i,diag_j,klo,2)),
                            static_cast<double>(fallc_arr(diag_i,diag_j,klo)),
                            static_cast<double>(dtcld),
                            static_cast<double>(delz_arr(diag_i,diag_j,klo)));
                std::fflush(stdout);
            }
#endif

            // Surface precipitation accumulation from sedimentation fallout
            ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                const Real fall_r = work1_arr(i,j,klo,0);       // rain fallout
                const Real fall_s = work1_arr(i,j,klo,1);       // snow fallout
                const Real fall_g = work1_arr(i,j,klo,2);       // graupel fallout
                const Real fall_c = fallc_arr(i,j,klo);         // ice fallout

                const Real fallsum     = fall_r + fall_s + fall_g + fall_c;
                const Real fallsum_qsi = fall_s + fall_c;
                const Real fallsum_qg  = fall_g;

                const Real conv = delz_arr(i,j,klo) / Real(rhoh2o) * dtcld * Real(1000.0);

                if (fallsum > Real(0.0)) {
                    rain_arr(i,j,klo) += fallsum * conv;
                }

                if (fallsum_qsi > Real(0.0)) {
                    snow_arr(i,j,klo) += fallsum_qsi * conv;
                }

                if (fallsum_qg > Real(0.0)) {
                    graup_arr(i,j,klo) += fallsum_qg * conv;
                }

                if (fallsum > Real(0.0)) {
                    sr_arr(i,j,klo) = (snow_arr(i,j,klo) + graup_arr(i,j,klo))
                                    / (rain_arr(i,j,klo) + Real(1.0e-12));
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G9 %3d %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(rain_arr(diag_i,diag_j,klo)),
                            static_cast<double>(sr_arr(diag_i,diag_j,klo)),
                            static_cast<double>(snow_arr(diag_i,diag_j,klo)),
                            static_cast<double>(graup_arr(diag_i,diag_j,klo)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3j: G10a Cloud-Ice Melt to Cloud Water
            // Exact port of bounded Fortran block (lines 1189-1197)
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G10A %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qi_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(xni_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(cpm_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supcol = t0c - t_arr(i,j,k);
                Real xlf = xls - xl_arr(i,j,k);
                if (supcol < Real(0.0)) xlf = xlf0;

                // Cloud-ice melt to cloud water (T > t0c)
                if (supcol < Real(0.0) && qi_arr(i,j,k) > Real(0.0)) {
                    const Real qim = qi_arr(i,j,k);  // preserve old ice amount

                    qc_arr(i,j,k)   += qim;                           // qci(:,:,1) += qci(:,:,2)
                    nc_arr(i,j,k)   += xni_arr(i,j,k);                // ncr(:,:,2) += xni
                    t_arr(i,j,k)    -= xlf / cpm_arr(i,j,k) * qim;    // latent heat release
                    qi_arr(i,j,k)    = Real(0.0);                     // zero out ice
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G10A %3d %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qi_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(xni_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // DISABLED: Steps 4-8 (CCN, condensation, autoconversion, accretion, evaporation)
            // ============================================================
            // These steps were incorrectly ordered BETWEEN G10a and G10b,
            // breaking parity with Fortran. They lack proper group boundary tags (WDM6-CPP_PRE_Gxx/POST_Gxx)
            // and should be ported as bounded groups in correct phase order.
            // Re-enable after proper group-based porting.
            // TODO: G11 (slope repack), G13a-G13g (warm-rain, accretion, etc), G16b (pcond+activation)
#if 0
            // Step 4: WDM6 CCN Activation
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real w_velocity = Real(0.1);
                wdm6_ccn_activation(
                    nc_arr(i,j,k), nn_arr(i,j,k),
                    qv_arr(i,j,k), qc_arr(i,j,k),
                    qsatw_arr(i,j,k), t_arr(i,j,k),
                    w_velocity, dtcld, ccn0, den_arr(i,j,k)
                );
            });

            // Step 5: Condensation/Evaporation
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (t_arr(i,j,k) > Real(t0c)) {
                    Real qcond = wdm6_conden(
                        t_arr(i,j,k), qv_arr(i,j,k), qsatw_arr(i,j,k),
                        xl_arr(i,j,k), cpm_arr(i,j,k),
                        Real(qmin), Real(rv)
                    );
                    pcond_arr(i,j,k) = qcond / dtcld;
                    qv_arr(i,j,k) -= qcond;
                    qc_arr(i,j,k) += qcond;
                    t_arr(i,j,k) += qcond * xl_arr(i,j,k) / cpm_arr(i,j,k);
                }
            });

            // Step 6: Autoconversion
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (qc_arr(i,j,k) > qc1_loc && nc_arr(i,j,k) > Real(1.e1)) {
                    Real mean_dia = wdm6_mean_droplet_diameter(
                        qc_arr(i,j,k), nc_arr(i,j,k), den_arr(i,j,k), pidnc_loc
                    );
                    if (mean_dia > Real(15.e-6)) {
                        Real auto_qc = qck1_loc * qc_arr(i,j,k) * qc_arr(i,j,k) * nc_arr(i,j,k);
                        auto_qc = amrex::min(auto_qc * dtcld, qc_arr(i,j,k));
                        Real auto_nc = auto_qc * nc_arr(i,j,k) / qc_arr(i,j,k) * Real(0.5);
                        auto_nc = amrex::min(auto_nc, nc_arr(i,j,k));
                        Real auto_nr = auto_nc * Real(0.01);
                        praut_arr(i,j,k) = auto_qc / dtcld;
                        qc_arr(i,j,k) -= auto_qc;
                        qr_arr(i,j,k) += auto_qc;
                        ncauto_arr(i,j,k) = auto_nc / dtcld;
                        nc_arr(i,j,k) -= auto_nc;
                        nrauto_arr(i,j,k) = auto_nr / dtcld;
                        nr_arr(i,j,k) += auto_nr;
                    }
                }
            });

            // Step 7: Accretion
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (qr_arr(i,j,k) > Real(1.e-9) && qc_arr(i,j,k) > Real(1.e-9)) {
                    Real accr_qc = Real(6.0) * qc_arr(i,j,k) * qr_arr(i,j,k);
                    accr_qc = amrex::min(accr_qc * dtcld, qc_arr(i,j,k));
                    Real accr_nc = accr_qc * nc_arr(i,j,k) / qc_arr(i,j,k);
                    accr_nc = amrex::min(accr_nc, nc_arr(i,j,k));
                    pracw_arr(i,j,k) = accr_qc / dtcld;
                    qc_arr(i,j,k) -= accr_qc;
                    qr_arr(i,j,k) += accr_qc;
                    ncaccr_arr(i,j,k) = accr_nc / dtcld;
                    nc_arr(i,j,k) -= accr_nc;
                }
            });

            // Step 8: Rain evaporation
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (qr_arr(i,j,k) > Real(1.e-9) && rhw_arr(i,j,k) < Real(1.0)) {
                    Real qrevp = Real(0.001) * qr_arr(i,j,k) * (Real(1.0) - rhw_arr(i,j,k));
                    qrevp = amrex::min(qrevp * dtcld, qr_arr(i,j,k));
                    Real nrevp = qrevp * nr_arr(i,j,k) / qr_arr(i,j,k);
                    nrevp = amrex::min(nrevp, nr_arr(i,j,k));
                    prevp_arr(i,j,k) = qrevp / dtcld;
                    qr_arr(i,j,k) -= qrevp;
                    qv_arr(i,j,k) += qrevp;
                    t_arr(i,j,k) -= qrevp * xl_arr(i,j,k) / cpm_arr(i,j,k);
                    nrevp_arr(i,j,k) = nrevp / dtcld;
                    nr_arr(i,j,k) -= nrevp;
                }
            });
#endif

            // ============================================================
            // Step 3k: G10b Cloud-Water Homogeneous Freezing
            // Exact port of bounded Fortran block (lines 1215-1224)
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G10B %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qi_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t0c - t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(cpm_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supcol = t0c - t_arr(i,j,k);
                Real xlf = xls - xl_arr(i,j,k);

                // Homogeneous freezing of cloud water when supcol > 40K (T < -40C)
                if (supcol > Real(40.0) && qc_arr(i,j,k) > Real(0.0)) {
                    const Real qc_old = qc_arr(i,j,k);

                    qi_arr(i,j,k) += qc_old;

                    if (nc_arr(i,j,k) > Real(0.0)) {
                        nc_arr(i,j,k) = Real(0.0);
                    }

                    t_arr(i,j,k) += xlf / cpm_arr(i,j,k) * qc_old;
                    qc_arr(i,j,k) = Real(0.0);
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G10B %3d %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qi_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t0c - t_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3l: G10c Cloud-Water Heterogeneous Freezing
            // Exact port of bounded Fortran block (pihtf, lines 1241-1258)
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G10C %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qi_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t0c - t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslopec_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(cpm_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(xl_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supcol = t0c - t_arr(i,j,k);

                // Heterogeneous freezing (Biggs, contact/immersion): 0 > T > -40C
                // Trigger: supcol > 0 (T < T0c) AND qc > qmin
                if (supcol > Real(0.0) && qc_arr(i,j,k) > Real(qmin)) {
                    const Real supcolt = amrex::min(supcol, Real(70.0));
                    const Real expterm = std::exp(pfrz2_loc * supcolt) - Real(1.0);

                    // Cloud droplet slope parameter cubed
                    const Real rs3 = rslopec_arr(i,j,k) * rslopec_arr(i,j,k) * rslopec_arr(i,j,k);

                    // Mass freezing rate: π² · pfrz1 · (exp(pfrz2·supcolt)-1) · (denr/den) · nc · rs3/18 · dtcld
                    Real pfrzdtc = pi_wdm6_loc * pi_wdm6_loc * pfrz1_loc * expterm
                                 * (denr / den_arr(i,j,k)) * nc_arr(i,j,k) * rs3
                                 / Real(18.0) * dtcld;
                    pfrzdtc = amrex::min(pfrzdtc, qc_arr(i,j,k));

                    // Number freezing rate: π · pfrz1 · (exp(pfrz2·supcolt)-1) · nc · rs3/6 · dtcld
                    Real nfrzdtc = pi_wdm6_loc * pfrz1_loc * expterm
                                 * nc_arr(i,j,k) * rs3
                                 / Real(6.0) * dtcld;
                    nfrzdtc = amrex::min(nfrzdtc, nc_arr(i,j,k));

                    // Apply number loss (only if nc > ncmin)
                    if (nc_arr(i,j,k) > Real(ncmin)) {
                        nc_arr(i,j,k) -= nfrzdtc;
                    }

                    // Apply mass and temperature updates
                    qi_arr(i,j,k) += pfrzdtc;
                    t_arr(i,j,k)  += xl_arr(i,j,k) / cpm_arr(i,j,k) * pfrzdtc;
                    qc_arr(i,j,k) -= pfrzdtc;
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                const Real rslopec_g10c = rslopec_arr(diag_i,diag_j,diag_k);
                std::printf("WDM6-CPP_POST_G10C %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qi_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t0c - t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslopec_g10c * rslopec_g10c * rslopec_g10c));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3m: G10d Rain-to-Graupel Freezing
            // Exact port of bounded Fortran block (pgfrz, rain freezing)
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G10D %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qg_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t0c - t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(cpm_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(xl_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supcol = t0c - t_arr(i,j,k);

                // Rain freezing to graupel: trigger when T < t0c and qr > 0
                if (supcol > Real(0.0) && qr_arr(i,j,k) > Real(0.0)) {
                    const Real supcolt = amrex::min(supcol, Real(70.0));
                    const Real expterm = std::exp(pfrz2_loc * supcolt) - Real(1.0);

                    // Rain slope cubed (from G4/G6 computation)
                    const Real rs3 = rslope3_arr(i,j,k,0);

                    // Mass freezing rate: 140*pi^2 * pfrz1 * nr * (denr/den) * exp(...) * rs3^2 * dtcld
                    Real pfrzdtr = Real(140.0) * pi_wdm6_loc * pi_wdm6_loc
                                 * pfrz1_loc * nr_arr(i,j,k)
                                 * (denr / den_arr(i,j,k))
                                 * expterm * rs3 * rs3 * dtcld;
                    pfrzdtr = amrex::min(pfrzdtr, qr_arr(i,j,k));

                    // Number freezing rate (conditional on nr > nrmin)
                    if (nr_arr(i,j,k) > Real(nrmin)) {
                        Real nfrzdtr = Real(4.0) * pi_wdm6_loc * pfrz1_loc
                                     * nr_arr(i,j,k) * expterm * rs3 * dtcld;
                        nfrzdtr = amrex::min(nfrzdtr, nr_arr(i,j,k));
                        nr_arr(i,j,k) -= nfrzdtr;
                    }

                    // Apply mass and temperature updates
                    qg_arr(i,j,k) += pfrzdtr;
                    t_arr(i,j,k)  += xl_arr(i,j,k) / cpm_arr(i,j,k) * pfrzdtr;
                    qr_arr(i,j,k) -= pfrzdtr;
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G10D %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qg_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t0c - t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,diag_k,0)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3k: G10e Phase cleanup — clamp number concentrations non-negative
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G10E %3d %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(nc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nr_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                nc_arr(i,j,k) = amrex::max(nc_arr(i,j,k), Real(0.0));
                nr_arr(i,j,k) = amrex::max(nr_arr(i,j,k), Real(0.0));
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G10E %3d %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(nc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nr_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // G11: SLOPE3 — Third slope_wdm6 call + avedia/rslopec recompute
            // ============================================================

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G11 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qs_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qg_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(p_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(xl_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qsatw_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qsati_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // G11(a): Repack fields for slope recomputation (Fortran qrs_tmp/ncr_tmp)
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                qrs_tmp_arr(i,j,k,0) = qr_arr(i,j,k);   // qrs(:,:,1)
                qrs_tmp_arr(i,j,k,1) = qs_arr(i,j,k);   // qrs(:,:,2)
                qrs_tmp_arr(i,j,k,2) = qg_arr(i,j,k);   // qrs(:,:,3)
                ncr_tmp_arr(i,j,k)   = nr_arr(i,j,k);   // ncr(:,:,3)
            });

            // G11(b): slope_wdm6-equivalent recompute for rain/snow/graupel
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real rain_rslope, rain_rslopeb, rain_rslope2, rain_rslope3, rain_vt, rain_vtn;
                Real snow_rslope, snow_rslopeb, snow_rslope2, snow_rslope3, snow_vt, snow_n0sfac;
                Real graup_rslope, graup_rslopeb, graup_rslope2, graup_rslope3, graup_vt;

                // Compute slope parameters using qrs_tmp/ncr_tmp values
                wdm6_slope_rain_cell(qrs_tmp_arr(i,j,k,0), ncr_tmp_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k),
                                     Real(qcrmin), Real(nrmin),
                                     rslopermax_loc, rsloperbmax_loc, rsloper2max_loc, rsloper3max_loc,
                                     Real(bvtr), pvtr_loc, pvtrn_loc, pidnr_loc,
                                     rain_rslope, rain_rslopeb, rain_rslope2, rain_rslope3,
                                     rain_vt, rain_vtn);
                wdm6_slope_snow_cell(qrs_tmp_arr(i,j,k,1), den_arr(i,j,k), denfac_arr(i,j,k), t_arr(i,j,k),
                                     pidn0s_loc, Real(alpha_wdm6), Real(n0smax), Real(n0s),
                                     Real(t0c), Real(qcrmin),
                                     rslopesmax_loc, rslopesbmax_loc, rslopes2max_loc, rslopes3max_loc,
                                     Real(bvts), pvts_loc,
                                     snow_rslope, snow_rslopeb, snow_rslope2, snow_rslope3, snow_vt,
                                     snow_n0sfac);
                wdm6_slope_graup_cell(qrs_tmp_arr(i,j,k,2), den_arr(i,j,k), denfac_arr(i,j,k),
                                      pidn0g_loc, Real(qcrmin),
                                      rslopegmax_loc, rslopegbmax_loc, rslopeg2max_loc, rslopeg3max_loc,
                                      slope_bvtg_loc, pvtg_loc,
                                      graup_rslope, graup_rslopeb, graup_rslope2, graup_rslope3, graup_vt);

                // Store results in output slope arrays
                rslope_arr(i,j,k,0) = rain_rslope;
                rslope_arr(i,j,k,1) = snow_rslope;
                rslope_arr(i,j,k,2) = graup_rslope;
                rslopeb_arr(i,j,k,0) = rain_rslopeb;
                rslopeb_arr(i,j,k,1) = snow_rslopeb;
                rslopeb_arr(i,j,k,2) = graup_rslopeb;
                rslope2_arr(i,j,k,0) = rain_rslope2;
                rslope2_arr(i,j,k,1) = snow_rslope2;
                rslope2_arr(i,j,k,2) = graup_rslope2;
                rslope3_arr(i,j,k,0) = rain_rslope3;
                rslope3_arr(i,j,k,1) = snow_rslope3;
                rslope3_arr(i,j,k,2) = graup_rslope3;
                work1_arr(i,j,k,0) = rain_vt;
                work1_arr(i,j,k,1) = snow_vt;
                work1_arr(i,j,k,2) = graup_vt;
                workn_arr(i,j,k) = rain_vtn;
            });

            // G11(c): avedia + rslopec/2/3 recompute (lamdac fallback branch)
            const Real cbrt24 = Real(std::pow(24.0f, 0.3333333f));
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // avedia component from rain slope: avedia(:,:,2) = rslope(:,:,1) * (24.)^(1/3)
                avedia_arr(i,j,k,1) = rslope_arr(i,j,k,0) * cbrt24;

                // Cloud water slope parameters: rslopec = 1/lamdac(qci(:,:,1), den, ncr(:,:,2))
                // Fortran gate: if (qci(:,:,1) <= qmin or ncr(:,:,2) <= ncmin) use maxima, else use lamdac
                const Real qci_for_lamdac = qc_arr(i,j,k);
                const Real nc_for_lamdac  = nc_arr(i,j,k);

                if (qci_for_lamdac <= Real(qmin) || nc_for_lamdac <= Real(ncmin)) {
                    rslopec_arr(i,j,k)  = rslopecmax_loc;
                    rslopec2_arr(i,j,k) = rslopec2max_loc;
                    rslopec3_arr(i,j,k) = rslopec3max_loc;
                } else {
                    const Real lamc = wdm6_lamdac(qci_for_lamdac, den_arr(i,j,k), nc_for_lamdac, pidnc_loc);
                    const Real rslc = Real(1.0) / lamc;
                    rslopec_arr(i,j,k)  = rslc;
                    rslopec2_arr(i,j,k) = rslc * rslc;
                    rslopec3_arr(i,j,k) = rslopec2_arr(i,j,k) * rslc;
                }

                // avedia component from cloud slope: avedia(:,:,1) = rslopec(:,:)
                avedia_arr(i,j,k,0) = rslopec_arr(i,j,k);
            });

            // G11(d): work arrays recompute via diffac/venfac
            // work1(:,:,1) = diffac(xl, p, t, den, qs(:,:,1), rv)
            // work1(:,:,2) = diffac(xls, p, t, den, qs(:,:,2), rv)
            // work2(:,:) = venfac(p, t, den, den0) [stored for G13a warm-rain rates]
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G12 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(xl_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(p_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(t_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(den_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qsatw_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qsati_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Recompute work1 diffusion coefficients using current state
                if (i == diag_i && j == diag_j && k == diag_k) {
#if !defined(AMREX_USE_GPU)
                    std::printf("WDM6-CPP_G11V_DIFFAC_PRE1 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(xl_arr(i,j,k)),
                                static_cast<double>(p_arr(i,j,k)),
                                static_cast<double>(t_arr(i,j,k)),
                                static_cast<double>(den_arr(i,j,k)),
                                static_cast<double>(qsatw_arr(i,j,k)),
                                static_cast<double>(rv));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT5 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(t_arr(i,j,k) * std::sqrt(t_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT6 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(t_arr(i,j,k) + Real(120.0)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT7 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>((t_arr(i,j,k) * std::sqrt(t_arr(i,j,k)))
                                                   /(t_arr(i,j,k) + Real(120.0))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT8 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.496e-6)
                                                   * (t_arr(i,j,k) * std::sqrt(t_arr(i,j,k)))
                                                   /(t_arr(i,j,k) + Real(120.0))
                                                   / den_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT21 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.496e-6)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT22 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.496e-6)
                                                   * (t_arr(i,j,k) * std::sqrt(t_arr(i,j,k)))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT23 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.496e-6)
                                                   * (t_arr(i,j,k) * std::sqrt(t_arr(i,j,k)))
                                                   /(t_arr(i,j,k) + Real(120.0))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT24 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.496e-6f)
                                                   * (t_arr(i,j,k) * std::sqrt(t_arr(i,j,k)))
                                                   /(t_arr(i,j,k) + Real(120.0))
                                                   / den_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT0 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))));
                    std::fflush(stdout);
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT1 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k)) * den_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT13 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT14 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   * den_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT15 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   * den_arr(i,j,k) * Real(rv)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT16 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   * den_arr(i,j,k) * Real(rv)
                                                   * t_arr(i,j,k) * t_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT29 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(xl_arr(i,j,k) * xl_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT30 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(den_arr(i,j,k) * xl_arr(i,j,k) * xl_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT31 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   * den_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT32 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   * den_arr(i,j,k) * Real(rv)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT33 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   * den_arr(i,j,k) * Real(rv)
                                                   * t_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT34 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   * den_arr(i,j,k) * Real(rv)
                                                   * t_arr(i,j,k) * t_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT35 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(den_arr(i,j,k) * xl_arr(i,j,k) * xl_arr(i,j,k)
                                                   /(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                     * den_arr(i,j,k) * Real(rv)
                                                     * t_arr(i,j,k) * t_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT36 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.0)
                                                   /(qsatw_arr(i,j,k)
                                                     * (Real(8.794e-5f)
                                                        * std::exp(std::log(t_arr(i,j,k)) * Real(1.81f))
                                                        / p_arr(i,j,k)))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT38 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.0)/qsatw_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT39 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.0)/(Real(8.794e-5f)
                                                       * std::exp(std::log(t_arr(i,j,k)) * Real(1.81f))
                                                       / p_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT40 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>((Real(1.0)/qsatw_arr(i,j,k))
                                                   *(Real(1.0)/(Real(8.794e-5f)
                                                      * std::exp(std::log(t_arr(i,j,k)) * Real(1.81f))
                                                      / p_arr(i,j,k)))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT41 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(qsatw_arr(i,j,k)
                                                   * (Real(8.794e-5f)
                                                      * std::exp(std::log(t_arr(i,j,k)) * Real(1.81f))
                                                      / p_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT37 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(den_arr(i,j,k) * xl_arr(i,j,k) * xl_arr(i,j,k)
                                                   /(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                     * den_arr(i,j,k) * Real(rv)
                                                     * t_arr(i,j,k) * t_arr(i,j,k))
                                                   + Real(1.0)
                                                   /(qsatw_arr(i,j,k)
                                                     * (Real(8.794e-5f)
                                                        * std::exp(std::log(t_arr(i,j,k)) * Real(1.81f))
                                                        / p_arr(i,j,k)))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT2 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(8.794e-5f) * std::exp(std::log(t_arr(i,j,k)) * Real(1.81f)) / p_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT24 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::log(t_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT25 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::log(t_arr(i,j,k)) * Real(1.81f)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT26 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::exp(std::log(t_arr(i,j,k)) * Real(1.81f))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT27 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(8.794e-5f) * std::exp(std::log(t_arr(i,j,k)) * Real(1.81f))
                                                   / p_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT28 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(8.794e-5f) * std::exp(std::log(t_arr(i,j,k)) * Real(1.81f))
                                                   / p_arr(i,j,k)));
                    std::fflush(stdout);
#endif
                }
                work1_arr(i,j,k,0) = wdm6_diffac(xl_arr(i,j,k),  p_arr(i,j,k), t_arr(i,j,k),
                                                  den_arr(i,j,k), qsatw_arr(i,j,k), Real(rv));  // qs(:,:,1)
                if (i == diag_i && j == diag_j && k == diag_k) {
#if !defined(AMREX_USE_GPU)
                    std::printf("WDM6-CPP_G11V_DIFFAC_POST1 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(work1_arr(i,j,k,0)));
                    std::fflush(stdout);
#endif
                }
                if (i == diag_i && j == diag_j && k == diag_k) {
#if !defined(AMREX_USE_GPU)
                    std::printf("WDM6-CPP_G11V_DIFFAC_PRE2 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(xls)),
                                static_cast<double>(p_arr(i,j,k)),
                                static_cast<double>(t_arr(i,j,k)),
                                static_cast<double>(den_arr(i,j,k)),
                                static_cast<double>(qsati_arr(i,j,k)),
                                static_cast<double>(rv));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT9 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(t_arr(i,j,k) * std::sqrt(t_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT10 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(t_arr(i,j,k) + Real(120.0)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT11 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>((t_arr(i,j,k) * std::sqrt(t_arr(i,j,k)))
                                                   /(t_arr(i,j,k) + Real(120.0))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT12 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.496e-6)
                                                   * (t_arr(i,j,k) * std::sqrt(t_arr(i,j,k)))
                                                   /(t_arr(i,j,k) + Real(120.0))
                                                   / den_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT3 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(den_arr(i,j,k) * Real(xls) * Real(xls)
                                                   /(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                     * den_arr(i,j,k) * Real(rv)
                                                     * t_arr(i,j,k) * t_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT17 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(den_arr(i,j,k) * Real(xls) * Real(xls)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT18 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT19 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   * den_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT20 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   * den_arr(i,j,k) * Real(rv)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT4 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.0) / (qsati_arr(i,j,k)
                                                   * (Real(8.794e-5) * std::exp(std::log(t_arr(i,j,k)) * Real(1.81))
                                                      / p_arr(i,j,k)))));
                    std::fflush(stdout);
#endif
                }
                work1_arr(i,j,k,1) = wdm6_diffac(Real(xls),       p_arr(i,j,k), t_arr(i,j,k),
                                                  den_arr(i,j,k), qsati_arr(i,j,k), Real(rv));  // qs(:,:,2)
                if (i == diag_i && j == diag_j && k == diag_k) {
#if !defined(AMREX_USE_GPU)
                    std::printf("WDM6-CPP_G11V_DIFFAC_POST2 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(work1_arr(i,j,k,1)));
                    std::fflush(stdout);
#endif
                }
                // Compute work2 ventilation factor for diffusion (used in G13a warm-rain rates)
#if !defined(AMREX_USE_GPU)
                if (i == diag_i && j == diag_j && k == diag_k) {
                    std::printf("WDM6-CPP_G11V_VENFAC_PRE %3d %24.16E %24.16E %24.16E %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(p_arr(i,j,k)),
                                static_cast<double>(t_arr(i,j,k)),
                                static_cast<double>(den_arr(i,j,k)),
                                static_cast<double>(Real(den0)));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT0 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT1 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(wdm6_diffus(t_arr(i,j,k), p_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT2 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   / wdm6_diffus(t_arr(i,j,k), p_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT3 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::log(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   / wdm6_diffus(t_arr(i,j,k), p_arr(i,j,k)))));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT4 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::log(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   / wdm6_diffus(t_arr(i,j,k), p_arr(i,j,k)))
                                                   * Real(0.3333333f)));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT5 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::exp(std::log(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   / wdm6_diffus(t_arr(i,j,k), p_arr(i,j,k)))
                                                   * Real(0.3333333f))));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT6 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::sqrt(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k)))));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT7 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(den0) / den_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT8 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::sqrt(Real(den0) / den_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT9 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::sqrt(std::sqrt(Real(den0) / den_arr(i,j,k)))));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT10 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::exp(std::log(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   / wdm6_diffus(t_arr(i,j,k), p_arr(i,j,k)))
                                                   * Real(0.3333333f))
                                                   / std::sqrt(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k)))));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT11 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::exp(std::log(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   / wdm6_diffus(t_arr(i,j,k), p_arr(i,j,k)))
                                                   * Real(0.3333333f))
                                                   / std::sqrt(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k)))
                                                   * std::sqrt(std::sqrt(Real(den0) / den_arr(i,j,k)))));
                    std::fflush(stdout);
                }
#endif
                work2_arr(i,j,k) = wdm6_venfac(p_arr(i,j,k), t_arr(i,j,k), den_arr(i,j,k), Real(den0));
#if !defined(AMREX_USE_GPU)
                if (i == diag_i && j == diag_j && k == diag_k) {
                    std::printf("WDM6-CPP_G11V_VENFAC_POST %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(work2_arr(i,j,k)));
                    std::fflush(stdout);
                }
#endif
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G12 %3d %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(work2_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G11 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(rslope_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(avedia_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslopec_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslopec2_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(avedia_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(work2_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,diag_k,1)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,diag_k,2)),
                            static_cast<double>(rslopec3_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,1)));
                std::fflush(stdout);
            }
#endif

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_PRE_G13A %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(qc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(qr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nc_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(rslopec3_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(work1_arr(diag_i,diag_j,diag_k,0)),
                            static_cast<double>(work2_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supsat = amrex::max(qv_arr(i,j,k), Real(qmin)) - qsatw_arr(i,j,k);
                const Real satdt = supsat / dtcld;
                const Real lencon = Real(2.7e-2) * den_arr(i,j,k) * qc_arr(i,j,k)
                    * (Real(1.0e20/16.0) * rslopec2_arr(i,j,k) * rslopec2_arr(i,j,k) - Real(0.4));
                const Real lenconcr = amrex::max(Real(1.2) * lencon, Real(qcrmin));

                if (qc_arr(i,j,k) > qcr_arr(i,j,k) && nc_arr(i,j,k) > Real(ncmin)) {
                    praut_arr(i,j,k) = qck1_loc * std::pow(qc_arr(i,j,k), Real(7.0/3.0))
                        * std::pow(nc_arr(i,j,k), Real(-1.0/3.0));
                    praut_arr(i,j,k) = amrex::min(praut_arr(i,j,k), qc_arr(i,j,k) / dtcld);

                    nrauto_arr(i,j,k) = Real(3.5e9) * den_arr(i,j,k) * praut_arr(i,j,k);
                    if (qr_arr(i,j,k) > lenconcr) {
                        nrauto_arr(i,j,k) = nr_arr(i,j,k) / qr_arr(i,j,k) * praut_arr(i,j,k);
                    }
                    nrauto_arr(i,j,k) = amrex::min(nrauto_arr(i,j,k), nc_arr(i,j,k) / dtcld);
                }

                if (qr_arr(i,j,k) >= lenconcr) {
                    if (avedia_arr(i,j,k,1) >= Real(di100)) {
                        nraccr_arr(i,j,k) = amrex::min(
                            Real(ncrk1) * nc_arr(i,j,k) * nr_arr(i,j,k)
                                * (rslopec3_arr(i,j,k) + Real(24.0) * rslope3_arr(i,j,k,0)),
                            nc_arr(i,j,k) / dtcld);
                        pracw_arr(i,j,k) = amrex::min(
                            pi_wdm6_loc / Real(6.0) * (Real(denr) / den_arr(i,j,k))
                                * Real(ncrk1) * nc_arr(i,j,k) * nr_arr(i,j,k)
                                * rslopec3_arr(i,j,k)
                                * (Real(2.0) * rslopec3_arr(i,j,k) + Real(24.0) * rslope3_arr(i,j,k,0)),
                            qc_arr(i,j,k) / dtcld);
                    } else {
                        nraccr_arr(i,j,k) = amrex::min(
                            Real(ncrk2) * nc_arr(i,j,k) * nr_arr(i,j,k)
                                * (Real(2.0) * rslopec3_arr(i,j,k) * rslopec3_arr(i,j,k)
                                   + Real(5040.0) * rslope3_arr(i,j,k,0) * rslope3_arr(i,j,k,0)),
                            nc_arr(i,j,k) / dtcld);
                        pracw_arr(i,j,k) = amrex::min(
                            pi_wdm6_loc / Real(6.0) * (Real(denr) / den_arr(i,j,k))
                                * Real(ncrk2) * nc_arr(i,j,k) * nr_arr(i,j,k)
                                * rslopec3_arr(i,j,k)
                                * (Real(6.0) * rslopec3_arr(i,j,k) * rslopec3_arr(i,j,k)
                                   + Real(5040.0) * rslope3_arr(i,j,k,0) * rslope3_arr(i,j,k,0)),
                            qc_arr(i,j,k) / dtcld);
                    }
                }

                if (avedia_arr(i,j,k,0) >= Real(di100)) {
                    nccol_arr(i,j,k) = Real(ncrk1) * nc_arr(i,j,k) * nc_arr(i,j,k) * rslopec3_arr(i,j,k);
                } else {
                    nccol_arr(i,j,k) = Real(2.0) * Real(ncrk2) * nc_arr(i,j,k) * nc_arr(i,j,k)
                        * rslopec3_arr(i,j,k) * rslopec3_arr(i,j,k);
                }

                if (qr_arr(i,j,k) >= lenconcr) {
                    if (avedia_arr(i,j,k,1) < Real(di100)) {
                        nrcol_arr(i,j,k) = Real(5040.0) * Real(ncrk2) * nr_arr(i,j,k) * nr_arr(i,j,k)
                            * rslope3_arr(i,j,k,0) * rslope3_arr(i,j,k,0);
                    } else if (avedia_arr(i,j,k,1) < Real(di600)) {
                        nrcol_arr(i,j,k) = Real(24.0) * Real(ncrk1) * nr_arr(i,j,k) * nr_arr(i,j,k)
                            * rslope3_arr(i,j,k,0);
                    } else if (avedia_arr(i,j,k,1) < Real(di2000)) {
                        const Real coecol = -Real(2.5e3) * (avedia_arr(i,j,k,1) - Real(di600));
                        nrcol_arr(i,j,k) = Real(24.0) * std::exp(coecol) * Real(ncrk1)
                            * nr_arr(i,j,k) * nr_arr(i,j,k) * rslope3_arr(i,j,k,0);
                    } else {
                        nrcol_arr(i,j,k) = Real(0.0);
                    }
                }

                if (qr_arr(i,j,k) > Real(0.0)) {
                    const Real coeres = rslope_arr(i,j,k,0)
                        * std::sqrt(rslope_arr(i,j,k,0) * rslopeb_arr(i,j,k,0));
                    prevp_arr(i,j,k) = (rhw_arr(i,j,k) - Real(1.0)) * nr_arr(i,j,k)
                        * (precr1_loc * rslope_arr(i,j,k,0) + precr2_loc * work2_arr(i,j,k) * coeres)
                        / work1_arr(i,j,k,0);
                    if (prevp_arr(i,j,k) < Real(0.0)) {
                        prevp_arr(i,j,k) = amrex::max(prevp_arr(i,j,k), -qr_arr(i,j,k) / dtcld);
                        prevp_arr(i,j,k) = amrex::max(prevp_arr(i,j,k), satdt / Real(2.0));

                        if (prevp_arr(i,j,k) == -qr_arr(i,j,k) / dtcld) {
                            nn_arr(i,j,k) = nn_arr(i,j,k) + nr_arr(i,j,k);
                            nr_arr(i,j,k) = Real(0.0);
                        }
                    } else {
                        prevp_arr(i,j,k) = amrex::min(prevp_arr(i,j,k), satdt / Real(2.0));
                    }
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && diag_col_in_tile) {
                std::printf("WDM6-CPP_POST_G13A %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            diag_k + 1,
                            static_cast<double>(praut_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nrauto_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(pracw_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nraccr_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nccol_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(nrcol_arr(diag_i,diag_j,diag_k)),
                            static_cast<double>(prevp_arr(diag_i,diag_j,diag_k)));
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 9: Ice physics (simplified) — remaining processes
            // ============================================================
            // These processes handle ice (qi), snow (qs), and graupel (qg)
            // following simplified versions of WSM6 processes

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real temp = t_arr(i,j,k);
                const Real t0c_loc = Real(273.15);

#if 0
                // DISABLED: Rain freezing (rain → graupel) — now handled by bounded G10d block above
                if (temp < t0c_loc && qr_arr(i,j,k) > Real(1.e-9)) {
                    Real frac_frz = (t0c_loc - temp) / Real(20.0) * Real(0.05);  // 5% per 20K
                    Real qfrz = amrex::min(frac_frz * qr_arr(i,j,k) * dtcld, qr_arr(i,j,k));
                    qg_arr(i,j,k) += qfrz;
                    qr_arr(i,j,k) -= qfrz;
                    // Also freeze corresponding rain number concentration
                    Real nfrz = qfrz * nr_arr(i,j,k) / amrex::max(qr_arr(i,j,k), Real(1.e-12));
                    nr_arr(i,j,k) -= amrex::min(nfrz, nr_arr(i,j,k));
                    t_arr(i,j,k) += qfrz * Real(xlf0) / Real(cpd);
                }
#endif

                // Ice → snow conversion (aggregation)
                if (qi_arr(i,j,k) > Real(1.e-6)) {
                    // Simplified: convert ice to snow at constant rate
                    Real qi_to_qs = amrex::min(qi_arr(i,j,k) * Real(0.001) * dtcld, qi_arr(i,j,k));
                    qs_arr(i,j,k) += qi_to_qs;
                    qi_arr(i,j,k) -= qi_to_qs;
                }

                // Melting: snow → rain
                if (temp > t0c_loc && qs_arr(i,j,k) > Real(1.e-9)) {
                    // Simplified melting rate
                    Real melt_rate = (temp - t0c_loc) * Real(0.01);  // 1% per K above 0C
                    Real qmelt = amrex::min(melt_rate * qs_arr(i,j,k) * dtcld, qs_arr(i,j,k));
                    qr_arr(i,j,k) += qmelt;
                    qs_arr(i,j,k) -= qmelt;
                    t_arr(i,j,k) -= qmelt * Real(xlf0) / Real(cpd);
                }

                // Melting: graupel → rain
                if (temp > t0c_loc && qg_arr(i,j,k) > Real(1.e-9)) {
                    Real melt_rate = (temp - t0c_loc) * Real(0.02);  // 2% per K (graupel melts faster)
                    Real qmelt = amrex::min(melt_rate * qg_arr(i,j,k) * dtcld, qg_arr(i,j,k));
                    qr_arr(i,j,k) += qmelt;
                    qg_arr(i,j,k) -= qmelt;
                    t_arr(i,j,k) -= qmelt * Real(xlf0) / Real(cpd);
                }

                // Riming: cloud water collected by snow → graupel
                if (temp < t0c_loc && qs_arr(i,j,k) > Real(1.e-9) && qc_arr(i,j,k) > Real(1.e-9)) {
                    // Simplified collection
                    Real qcol = amrex::min(qs_arr(i,j,k) * qc_arr(i,j,k) * Real(10.0) * dtcld, qc_arr(i,j,k));
                    qg_arr(i,j,k) += qcol;
                    qc_arr(i,j,k) -= qcol;
                    // Also remove corresponding cloud droplets
                    Real ncol = qcol * nc_arr(i,j,k) / amrex::max(qc_arr(i,j,k), Real(1.e-12));
                    nc_arr(i,j,k) -= amrex::min(ncol, nc_arr(i,j,k));
                    t_arr(i,j,k) += qcol * Real(xlf0) / Real(cpd);
                }
            });

            // ============================================================
            // Step 10: Enforce bounds
            // ============================================================
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                qc_arr(i,j,k) = amrex::max(qc_arr(i,j,k), Real(0.0));
                qr_arr(i,j,k) = amrex::max(qr_arr(i,j,k), Real(0.0));
                qi_arr(i,j,k) = amrex::max(qi_arr(i,j,k), Real(0.0));
                qs_arr(i,j,k) = amrex::max(qs_arr(i,j,k), Real(0.0));
                qg_arr(i,j,k) = amrex::max(qg_arr(i,j,k), Real(0.0));
                nc_arr(i,j,k) = amrex::max(nc_arr(i,j,k), Real(1.e1));
                nr_arr(i,j,k) = amrex::max(nr_arr(i,j,k), Real(1.e-2));
                nn_arr(i,j,k) = amrex::max(nn_arr(i,j,k), Real(0.0));
            });

        } // End minor timestep loop

        // ============================================================
        // Step 11: Simplified sedimentation for all species
        // ============================================================
        // Simplified top-down sedimentation for rain, snow, graupel, and ice
        const Real dz_sedi = m_geom.CellSize(2);
        const Real t0c_sed = Real(273.15);

        // WDM6 constants for terminal velocity
        const Real avts_loc = Real(11.72);
        const Real bvts_loc = Real(0.41);
        const Real avtg_loc = Real(330.0);  // Graupel fall speed coefficient
        const Real bvtg_loc = Real(0.89);

        ParallelFor(amrex::makeSlab(box,2,klo), [=] AMREX_GPU_DEVICE (int i, int j, int) {
            Real precip_rain = Real(0.0);
            Real precip_snow = Real(0.0);
            Real precip_graup = Real(0.0);

            // Top-down sedimentation for all species
            for (int kk = khi; kk >= klo; --kk) {
                const Real den_loc = den_arr(i,j,kk);
                const Real denfac = std::sqrt(Real(1.28) / den_loc);

                // Snow sedimentation (single-moment)
                if (qs_arr(i,j,kk) > Real(1.e-9)) {
                    // Temperature-dependent N0 factor
                    const Real temp = t_arr(i,j,kk);
                    const Real supcol = t0c_sed - temp;
                    const Real n0sfac = amrex::max(amrex::min(std::exp(Real(0.12) * supcol),
                                                               Real(1.0e11) / Real(2.0e6)), Real(1.0));

                    // Snow terminal velocity (simplified)
                    Real lamdasq = std::pow(Real(3.14159) * Real(2.0e6) * n0sfac / (den_loc * qs_arr(i,j,kk)), Real(0.5));
                    Real vt = avts_loc * std::pow(lamdasq, -bvts_loc) * denfac;
                    Real flux_qs = amrex::min(qs_arr(i,j,kk) * vt * dt_advance / dz_sedi, qs_arr(i,j,kk));

                    qs_arr(i,j,kk) -= flux_qs;

                    if (kk > klo) {
                        qs_arr(i,j,kk-1) += flux_qs;
                    } else {
                        precip_snow += flux_qs * den_loc * dz_sedi;
                    }
                }

                // Graupel sedimentation (single-moment)
                if (qg_arr(i,j,kk) > Real(1.e-9)) {
                    // Graupel terminal velocity (simplified, faster than snow)
                    Real lamdag = std::pow(Real(3.14159) * Real(4.0e5) / (den_loc * qg_arr(i,j,kk)), Real(0.25));
                    Real vt = avtg_loc * std::pow(Real(1.0) / lamdag, bvtg_loc) * denfac;
                    Real flux_qg = amrex::min(qg_arr(i,j,kk) * vt * dt_advance / dz_sedi, qg_arr(i,j,kk));

                    qg_arr(i,j,kk) -= flux_qg;

                    if (kk > klo) {
                        qg_arr(i,j,kk-1) += flux_qg;
                    } else {
                        precip_graup += flux_qg * den_loc * dz_sedi;
                    }
                }

                // Ice crystal sedimentation (very slow, often negligible)
                if (qi_arr(i,j,kk) > Real(1.e-9)) {
                    // Ice crystals fall slowly (~0.1-1 m/s)
                    Real vt = Real(0.3) * denfac;  // Simplified constant fall speed
                    Real flux_qi = amrex::min(qi_arr(i,j,kk) * vt * dt_advance / dz_sedi, qi_arr(i,j,kk));

                    qi_arr(i,j,kk) -= flux_qi;

                    if (kk > klo) {
                        qi_arr(i,j,kk-1) += flux_qi;
                    }
                    // Ice crystals don't typically reach surface before converting to snow
                }
            }

            // Accumulate surface precipitation (convert to mm)
            rain_arr(i,j,klo) += precip_rain;
            snow_arr(i,j,klo) += precip_snow;
            graup_arr(i,j,klo) += precip_graup;
        });

#ifdef ERF_USE_WDM6_FORT
        }
#endif

    } // MFIter loop
}

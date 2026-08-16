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
    return (wdm6_literal(8.794e-5)*std::exp(std::log(x)*wdm6_literal(1.81)))/y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_viscos (Real x, Real y) {
    return wdm6_literal(1.496e-6)*(x*std::sqrt(x))/(x+Real(120.0))/y;
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
                   *wdm6_literal(0.3333333))
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
    // Fortran statement function, ERF_module_mp_wdm6.F90:3277 and :3357:
    //     lamdar(x,y,z) = exp(log(((pidnr*z)/(x*y)))*((.33333333)))
    // called as lamdar(qrs(i,k,1), den(i,k), ncr(i,k)), so x=qr, y=den, z=nr:
    // the denominator is qr*den, NOT qr alone.
    //
    // This carried the identical defect already fixed in wdm6_lamdac: it wrote
    // (pidnr*nr*den)/(den*qr), in which den cancels, dropping the /den entirely
    // and inflating rslope by den^(-1/3). Measured at (199,3,43): the rain slope
    // was 3.098e-05 native against 2.833e-05 on the bridge, a ratio of 1.0936
    // against den = 0.765, matching den^(-1/3) exactly. That pushed avedia
    // across the di82 threshold in G16a, 8.9365e-05 native against 8.1717e-05,
    // so the bridge collapsed rain back into cloud (qr and nr to zero, qc and nc
    // incremented) while the native path left the cell untouched.
    //
    // The old qr/nr guard is dropped: every caller already gates on
    // qr <= qcrmin || nr <= nrmin with the same values, and the branch returned
    // 1/5.0e4, a slope, where the caller expects a lamda and inverts it again.
    return std::exp(std::log((pidnr_arg * nr) / (qr * den)) * wdm6_literal(0.33333333));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_lamdac_exact (Real qc, Real den, Real nc, Real pidnc_arg) {
    // Fortran statement function, ERF_module_mp_wdm6.F90:562:
    //     lamdac(x,y,z) = exp(log(((pidnc*z)/(x*y)))*((.33333333)))
    // called as lamdac(qci(i,k,1), den(i,k), ncr(i,k,2)), so x=qc, y=den, z=nc.
    // The exponent literal is unsuffixed in the Fortran and therefore obeys the
    // LITERAL PRECISION CONTRACT; it is not 1/3.
    return std::exp(std::log((pidnc_arg * nc) / (qc * den)) * wdm6_literal(0.33333333));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_rslopec_exact (Real qc, Real den, Real nc, Real pidnc_arg) {
    // rslopec is the RECIPROCAL of lamdac: the Fortran writes
    //     rslopec(i,k) = 1./lamdac(qci(i,k,1),den(i,k),ncr(i,k,2))
    // at both of its two sites (G3 at :915 and G11 at :1544). Returning lamdac
    // itself here put rslopec off by 1/lamdac^2, about nine orders of magnitude
    // once cloud water exists. It was invisible until cloud water first formed,
    // because with qc <= qmin both legs take the rslopecmax branch instead.
    return Real(1.0) / wdm6_lamdac_exact(qc, den, nc, pidnc_arg);
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
        // The rain-slope cap is unsuffixed in the Fortran at BOTH of its sites,
        // slope_wdm6 :3420 and slope_rain :3493:
        //     rslope(i,k,1) = min(1./lamdar(...),1.e-3)
        // so it obeys the LITERAL PRECISION CONTRACT: float32(1.e-3) widens to
        // 0.0010000000474974513, not the exact double 1e-3. The cap therefore
        // sits +4.749745e-08 relative ABOVE the port's, and the gap propagates
        // as a power: rslopeb = rslope**bvtr, rslope2, rslope3, so rslope3
        // carries 3x it at 1.424924e-07.
        //
        // This is a min(), so it binds only where 1./lamdar exceeds the cap --
        // rare, and identically zero for the first fifteen steps of the Bubble
        // case, which is why a bitwise-clean Milestone B did not expose it.
        rslope  = amrex::min(Real(1.0) / wdm6_lamdar(qr, den, nr, pidnr_arg),
                             wdm6_literal(1.e-3));
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

// The four slope routines -- slope_wdm6 :3290, slope_rain :3372, slope_snow
// :3417 and slope_graup :3465 -- each declare their OWN
//     real(kind=kind_phys), PARAMETER :: t0c = 273.15
// shadowing the t0c the caller passes in. That local literal is unsuffixed, so
// it is float32(273.15) = 273.149993896484375 against the caller's true double
// 273.14999999999997726, a difference of -6.104e-06. Inside those routines t0c
// feeds only supcol and hence n0sfac, so only the SNOW slope is affected;
// everywhere else in the scheme the caller's double t0c is the correct value
// and must be left alone.
//
// Confirmed by reconstruction at (199,3,91) with all PRE_G4 inputs bitwise
// equal, matching both legs to 20 digits:
//     caller t0c, double  -> rslope 0.00068326407135771243233  = native
//     local  t0c, float32 -> rslope 0.000683264196467108691    = bridge
constexpr amrex::Real wdm6_slope_t0c = wdm6_literal(273.15);

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
            // Fortran backs both brackets off by one at the TOP of every k
            // before searching forward:
            //     kb=max(kb-1,1) ; kt=max(kt-1,1)
            // Omitting this let the forward search start above the correct
            // bracket, so cells whose departure interval reaches back below the
            // previous k's bracket were remapped from the wrong donor.
            kb = amrex::max(kb - 1, 0);
            kt = amrex::max(kt - 1, 0);

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
            // Fortran writes a bare kt = kt - 1 with no floor. Clamping it to 0
            // is not equivalent: when the decrement should drop kt BELOW kb,
            // the Fortran falls through both branches and leaves qn(k) at zero,
            // whereas the clamp can make kt equal kb and take the single-cell
            // branch against the wrong donor. kt is only used as an index
            // inside the two branches, both of which require kt >= kb >= 0, so
            // a negative value here is safe.
            kt = kt - 1;

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

    // Leg label for bridge-wrapper site tags. These sites live in the bridge
    // caller and have no Fortran counterpart, so both legs are emitted from C++;
    // the label follows erf.use_wdm6_cpp_answer, i.e. which answer path ran, not
    // which file did the printing. A single emission path serves both legs so the
    // two records cannot drift in format.
#ifdef ERF_USE_WDM6_FORT
    const char* const w1_leg = run_wdm6_fort ? "FORT" : "CPP";
#else
    const char* const w1_leg = "CPP";
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
    constexpr double dens = static_cast<double>(rhos);
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

        // ------------------------------------------------------------------
        // Bridge-wrapper site W1_THETAWB (see group_map.md, "Bridge-wrapper
        // sites"). NOT a file-order group: no counterpart exists in
        // ERF_module_mp_wdm6.F90, so it can never join the P0..G17 sequence.
        // Contractual schema, one emitter for both legs:
        //     WDM6-<FORT|CPP>_<PRE|POST>_W1_THETAWB k theta t p exner t_over_exner
        // Full column kts..kte per the WSM6 all-k print contract; the failing
        // Milestone A zone is k=90, which a kts-only tag cannot observe.
        // ------------------------------------------------------------------
        auto const& w1_theta = mic_fab_vars[MicVar_WDM6::theta]->array(mfi);
        const bool w1_diag = (microphysics_debug > 0 &&
                              diag_i >= ilo && diag_i <= ihi &&
                              diag_j >= jlo && diag_j <= jhi);
        auto emit_w1 = [&] (const char* stage)
        {
#if !defined(AMREX_USE_GPU)
            if (!w1_diag) return;
            constexpr Real p0_w1 = 1.e5;
            constexpr Real rdOcp_w1 = R_d / Cp_d;
            for (int k = klo; k <= khi; ++k) {
                const Real ex = std::pow(p_arr(diag_i,diag_j,k) / p0_w1, rdOcp_w1);
                std::printf("WDM6-%s_%s_W1_THETAWB %3d %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            w1_leg, stage, k + 1,
                            static_cast<double>(w1_theta(diag_i,diag_j,k)),
                            static_cast<double>(t_arr(diag_i,diag_j,k)),
                            static_cast<double>(p_arr(diag_i,diag_j,k)),
                            static_cast<double>(ex),
                            static_cast<double>(t_arr(diag_i,diag_j,k) / ex));
            }
            std::fflush(stdout);
#else
            amrex::ignore_unused(stage);
#endif
        };

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

        Box box2d(box);
        box2d.makeSlab(2, 0);
        Box fab_box2d(fab_box);
        fab_box2d.makeSlab(2, 0);

        // Create landmask array (xland: 0=water, 1=land)
        // TODO: Get from ERF's lmask_lev when available
        // For now, default to land (continental CCN)
        FArrayBox xland_fab(fab_box2d, 1, Arena_Used);
        auto const& xland_arr = xland_fab.array();
        ParallelFor(fab_box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            xland_arr(i,j,k) = Real(1.0);  // Default to land
        });

        // Create 2D accumulation arrays
        // Fortran bridge uses ims:ime, jms:jme storage bounds; these buffers must
        // therefore be allocated on fab_box extents even if C++ kernels only
        // update the valid tile slab (box2d).
        FArrayBox rainacc_fab(fab_box2d, 1, Arena_Used);
        FArrayBox rainncv_fab(fab_box2d, 1, Arena_Used);
        FArrayBox sr_fab(fab_box2d, 1, Arena_Used);
        FArrayBox snowacc_fab(fab_box2d, 1, Arena_Used);
        FArrayBox snowncv_fab(fab_box2d, 1, Arena_Used);
        FArrayBox graupacc_fab(fab_box2d, 1, Arena_Used);
        FArrayBox graupelncv_fab(fab_box2d, 1, Arena_Used);

        auto const& rainacc_arr = rainacc_fab.array();
        auto const& rainncv_arr = rainncv_fab.array();
        auto const& sr_arr = sr_fab.array();
        auto const& snowacc_arr = snowacc_fab.array();
        auto const& snowncv_arr = snowncv_fab.array();
        auto const& graupacc_arr = graupacc_fab.array();
        auto const& graupelncv_arr = graupelncv_fab.array();

        // Initialize 2D arrays to zero over the full storage slab so the Fortran
        // bridge never reads ghost entries that were never written.
        ParallelFor(fab_box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
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

        // Bridge-wrapper site W1_THETAWB brackets the writeback below.
        emit_w1("PRE");

        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // Recompute theta from updated temperature
            // exner = (p/p0)^(R/cp)
            // theta = T / exner = T * (p0/p)^(R/cp)
            Real exner = std::pow(p_arr(i,j,k) / p0, rdOcp);
            theta_arr(i,j,k) = t_arr(i,j,k) / exner;
        });

        emit_w1("POST");

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
        FArrayBox pidep_fab(fab_box,1, Arena_Used);
        FArrayBox psdep_fab(fab_box,1, Arena_Used);
        FArrayBox pgdep_fab(fab_box,1, Arena_Used);
        FArrayBox pigen_fab(fab_box,1, Arena_Used);
        FArrayBox psaut_fab(fab_box,1, Arena_Used);
        FArrayBox pgaut_fab(fab_box,1, Arena_Used);
        FArrayBox pcact_fab(fab_box,1, Arena_Used);
        FArrayBox pcond_fab(fab_box,1, Arena_Used);
        FArrayBox praci_fab(fab_box,1, Arena_Used);
        FArrayBox piacr_fab(fab_box,1, Arena_Used);
        FArrayBox niacr_fab(fab_box,1, Arena_Used);
        FArrayBox psaci_fab(fab_box,1, Arena_Used);
        FArrayBox pgaci_fab(fab_box,1, Arena_Used);
        FArrayBox psacw_fab(fab_box,1, Arena_Used);
        FArrayBox nsacw_fab(fab_box,1, Arena_Used);
        FArrayBox pgacw_fab(fab_box,1, Arena_Used);
        FArrayBox ngacw_fab(fab_box,1, Arena_Used);
        FArrayBox paacw_fab(fab_box,1, Arena_Used);
        FArrayBox naacw_fab(fab_box,1, Arena_Used);
        FArrayBox pracs_fab(fab_box,1, Arena_Used);
        FArrayBox psacr_fab(fab_box,1, Arena_Used);
        FArrayBox nsacr_fab(fab_box,1, Arena_Used);
        FArrayBox pgacr_fab(fab_box,1, Arena_Used);
        FArrayBox ngacr_fab(fab_box,1, Arena_Used);
        FArrayBox pgacs_fab(fab_box,1, Arena_Used);
        FArrayBox pseml_fab(fab_box,1, Arena_Used);
        FArrayBox nseml_fab(fab_box,1, Arena_Used);
        FArrayBox pgeml_fab(fab_box,1, Arena_Used);
        FArrayBox ngeml_fab(fab_box,1, Arena_Used);
        FArrayBox psevp_fab(fab_box,1, Arena_Used);
        FArrayBox pgevp_fab(fab_box,1, Arena_Used);

        // Number concentration process rates (WDM6-specific)
        FArrayBox ncauto_fab(fab_box,1, Arena_Used);  // nc lost to autoconversion
        FArrayBox ncaccr_fab(fab_box,1, Arena_Used);  // nc lost to accretion by rain
        FArrayBox nrauto_fab(fab_box,1, Arena_Used);  // nr gained from autoconversion
        FArrayBox nraccr_fab(fab_box,1, Arena_Used);  // nr gained from accretion
        FArrayBox nrevp_fab(fab_box,1, Arena_Used);   // nr lost to evaporation
        FArrayBox ncact_fab(fab_box,1, Arena_Used);   // nn -> nc activation rate
        FArrayBox act_ratio_fab(fab_box,1, Arena_Used);
        FArrayBox act_fraction_fab(fab_box,1, Arena_Used);
        FArrayBox act_raw_fab(fab_box,1, Arena_Used);
        FArrayBox act_cap_fab(fab_box,1, Arena_Used);
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
        auto const& pidep_arr = pidep_fab.array();
        auto const& psdep_arr = psdep_fab.array();
        auto const& pgdep_arr = pgdep_fab.array();
        auto const& pigen_arr = pigen_fab.array();
        auto const& psaut_arr = psaut_fab.array();
        auto const& pgaut_arr = pgaut_fab.array();
        auto const& pcact_arr = pcact_fab.array();
        auto const& pcond_arr = pcond_fab.array();
        auto const& praci_arr = praci_fab.array();
        auto const& piacr_arr = piacr_fab.array();
        auto const& niacr_arr = niacr_fab.array();
        auto const& psaci_arr = psaci_fab.array();
        auto const& pgaci_arr = pgaci_fab.array();
        auto const& psacw_arr = psacw_fab.array();
        auto const& nsacw_arr = nsacw_fab.array();
        auto const& pgacw_arr = pgacw_fab.array();
        auto const& ngacw_arr = ngacw_fab.array();
        auto const& paacw_arr = paacw_fab.array();
        auto const& naacw_arr = naacw_fab.array();
        auto const& pracs_arr = pracs_fab.array();
        auto const& psacr_arr = psacr_fab.array();
        auto const& nsacr_arr = nsacr_fab.array();
        auto const& pgacr_arr = pgacr_fab.array();
        auto const& ngacr_arr = ngacr_fab.array();
        auto const& pgacs_arr = pgacs_fab.array();
        auto const& pseml_arr = pseml_fab.array();
        auto const& nseml_arr = nseml_fab.array();
        auto const& pgeml_arr = pgeml_fab.array();
        auto const& ngeml_arr = ngeml_fab.array();
        auto const& psevp_arr = psevp_fab.array();
        auto const& pgevp_arr = pgevp_fab.array();
        auto const& ncauto_arr = ncauto_fab.array();
        auto const& ncaccr_arr = ncaccr_fab.array();
        auto const& nrauto_arr = nrauto_fab.array();
        auto const& nraccr_arr = nraccr_fab.array();
        auto const& nrevp_arr = nrevp_fab.array();
        auto const& ncact_arr = ncact_fab.array();
        auto const& act_ratio_arr = act_ratio_fab.array();
        auto const& act_fraction_arr = act_fraction_fab.array();
        auto const& act_raw_arr = act_raw_fab.array();
        auto const& act_cap_arr = act_cap_fab.array();
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
            // CCN is clamped into [1.e8, 2.e10], not merely floored at zero:
            //     ncr(i,k,1) = min(max(ncr(i,k,1),1.e8),2.e10)
            // at ERF_module_mp_wdm6.F90:604. This re-seeds the aerosol field to
            // the background concentration at the top of every macro timestep,
            // so activation depletion never carries across steps. A plain
            // max(.,0) is a no-op on step 1, where Init() has just set nn to
            // ccn0 = 1.e8 everywhere, and only diverges from step 2 onward once
            // activation has drawn nn below the floor. Both literals are exact
            // in float32, so the precision contract is satisfied either way.
            nn_arr(i,j,k) = amrex::min(amrex::max(nn_arr(i,j,k), wdm6_literal(1.0e8)),
                                       wdm6_literal(2.0e10));
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
        const Real pacrc_loc = m_pacrc;
        const Real pacrg_loc = m_pacrg;
        const Real g4pbr_loc = m_g4pbr;
        const Real g7pbr_loc = m_g7pbr;
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
                const Real ttp = Real(t0c) + wdm6_literal(0.01);
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
                pidep_arr(i,j,k) = Real(0.0);
                psdep_arr(i,j,k) = Real(0.0);
                pgdep_arr(i,j,k) = Real(0.0);
                pigen_arr(i,j,k) = Real(0.0);
                psaut_arr(i,j,k) = Real(0.0);
                pgaut_arr(i,j,k) = Real(0.0);
                pcact_arr(i,j,k) = Real(0.0);
                pcond_arr(i,j,k) = Real(0.0);
                praci_arr(i,j,k) = Real(0.0);
                piacr_arr(i,j,k) = Real(0.0);
                niacr_arr(i,j,k) = Real(0.0);
                psaci_arr(i,j,k) = Real(0.0);
                pgaci_arr(i,j,k) = Real(0.0);
                psacw_arr(i,j,k) = Real(0.0);
                nsacw_arr(i,j,k) = Real(0.0);
                pgacw_arr(i,j,k) = Real(0.0);
                ngacw_arr(i,j,k) = Real(0.0);
                paacw_arr(i,j,k) = Real(0.0);
                naacw_arr(i,j,k) = Real(0.0);
                pracs_arr(i,j,k) = Real(0.0);
                psacr_arr(i,j,k) = Real(0.0);
                nsacr_arr(i,j,k) = Real(0.0);
                pgacr_arr(i,j,k) = Real(0.0);
                ngacr_arr(i,j,k) = Real(0.0);
                pgacs_arr(i,j,k) = Real(0.0);
                pseml_arr(i,j,k) = Real(0.0);
                nseml_arr(i,j,k) = Real(0.0);
                pgeml_arr(i,j,k) = Real(0.0);
                ngeml_arr(i,j,k) = Real(0.0);
                psevp_arr(i,j,k) = Real(0.0);
                pgevp_arr(i,j,k) = Real(0.0);
                ncauto_arr(i,j,k) = Real(0.0);
                ncaccr_arr(i,j,k) = Real(0.0);
                nrauto_arr(i,j,k) = Real(0.0);
                nraccr_arr(i,j,k) = Real(0.0);
                nrevp_arr(i,j,k) = Real(0.0);
                ncact_arr(i,j,k) = Real(0.0);
                act_ratio_arr(i,j,k) = Real(0.0);
                act_fraction_arr(i,j,k) = Real(0.0);
                act_raw_arr(i,j,k) = Real(0.0);
                act_cap_arr(i,j,k) = Real(0.0);
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
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G4 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(den_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
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
                                     wdm6_slope_t0c, Real(qcrmin),
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
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G4 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(rslope_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(workn_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(den_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,kdbg,2)));
                std::fflush(stdout);
                }
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
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G5A %3d %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(workn_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(delz_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(dtcld));
                std::fflush(stdout);
                }
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
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G5A %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(workn_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(mstep_arr(diag_i,diag_j,0)),
                            static_cast<double>(mstepmax),
                            static_cast<double>(numdt_arr(diag_i,diag_j,0)),
                            static_cast<double>(delz_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(dtcld));
                std::fflush(stdout);
                }
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
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G5B %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(workn_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(mstep_arr(diag_i,diag_j,0)),
                            static_cast<double>(delz_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(den_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
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
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G5B %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(workn_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(delz_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(dtcld));
                std::fflush(stdout);
                }
            }
#endif

            // ============================================================
            // Step 3f: Snow/graupel sedimentation (G5c)
            // Exact bounded NISLFV_SG slice plus lower-boundary slab fall.
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                const Real qsum = amrex::max(qs_arr(diag_i,diag_j,kdbg) + qg_arr(diag_i,diag_j,kdbg), Real(1.0e-15));
                const Real worka = (qsum > Real(1.0e-15))
                    ? (work1_arr(diag_i,diag_j,kdbg,1) * qs_arr(diag_i,diag_j,kdbg)
                     + work1_arr(diag_i,diag_j,kdbg,2) * qg_arr(diag_i,diag_j,kdbg)) / qsum
                    : Real(0.0);
                std::printf("WDM6-CPP_PRE_G5C %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(worka),
                            static_cast<double>(den_arr(diag_i,diag_j,kdbg) * qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(den_arr(diag_i,diag_j,kdbg) * qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(delz_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(dtcld),
                            static_cast<double>(denfac_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
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
                    pidn0s_loc, pidn0g_loc, Real(qcrmin), Real(alpha_wdm6), Real(n0smax), Real(n0s), wdm6_slope_t0c,
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

            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                const Real qsum = amrex::max(qs_arr(diag_i,diag_j,kdbg) + qg_arr(diag_i,diag_j,kdbg), Real(1.0e-15));
                const Real worka = (qsum > Real(1.0e-15))
                    ? (work1_arr(diag_i,diag_j,kdbg,1) * qs_arr(diag_i,diag_j,kdbg)
                     + work1_arr(diag_i,diag_j,kdbg,2) * qg_arr(diag_i,diag_j,kdbg)) / qsum
                    : Real(0.0);
                std::printf("WDM6-CPP_POST_G5C %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(worka),
                            static_cast<double>(den_arr(diag_i,diag_j,kdbg) * qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(den_arr(diag_i,diag_j,kdbg) * qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(work1_arr(diag_i,diag_j,klo,1)),
                            static_cast<double>(work1_arr(diag_i,diag_j,klo,2)));
                std::fflush(stdout);
                }
            }
#endif

            // ============================================================
            // Step 3g: Second slope_wdm6 call after sedimentation (G6)
            // Recompute slope parameters for all species after G5c sed.
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G6 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(den_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
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
                                     wdm6_slope_t0c, Real(qcrmin),
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
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G6 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(rslope_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(workn_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

            // ============================================================
            // WDM6-CPP_PRE_G7
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G7 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(p_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(den_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(cpm_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslope2_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,kdbg,2)));
                std::fflush(stdout);
                }
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
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G7 %3d %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

            // ============================================================
            // Step 3i: G8 ice fall speed + ice sedimentation
            // VICE: Ice crystal terminal velocity + nislfv_rain_plmr
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G8 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qi_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(xni_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(den_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(delz_arr(diag_i,diag_j,kdbg)));
                }
                std::fflush(stdout);
            }
#endif

            // Compute ice crystal fall speed (work1c). Fortran:
            //     xmi = den*qci(:,:,2)/xni
            //     diameter = max(min(dicon*sqrt(xmi),dimax),1.e-25)
            //     work1c   = 1.49e4*exp(log(diameter)*(1.31))
            // dicon and dimax are the CLASS constants, already routed through
            // wdm6_literal. This block previously declared its own locals that
            // shadowed them, and both were wrong: dicon was 11.45e-9 against the
            // Fortran's 11.9, nine orders of magnitude small, which collapsed the
            // ice fall speed to ~5e-13 m/s so the native path effectively did not
            // sediment ice at all; and dimax was a true double, reintroducing at
            // this site exactly the literal-precision defect already closed in
            // G13e. Hoisted out of the lambda in the established _loc style so no
            // capture of `this` is needed on device.
            constexpr Real dicon_loc = dicon;
            constexpr Real dimax_loc = dimax;
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real work1c = Real(0.0);
                if (qi_arr(i,j,k) > Real(0.0)) {
                    const Real xni_safe = amrex::max(xni_arr(i,j,k), Real(1.0e-30));
                    const Real xmi = den_arr(i,j,k) * qi_arr(i,j,k) / xni_safe;
                    const Real diameter = amrex::max(amrex::min(dicon_loc * std::sqrt(xmi), dimax_loc),
                                                     wdm6_literal(1.0e-25));
                    work1c = wdm6_literal(1.49e4) * std::exp(std::log(diameter) * wdm6_literal(1.31));
                }
                work1c_arr(i,j,k) = work1c;
#if !defined(AMREX_USE_GPU)
                // Tier 2 forensic decomposition of the G8 ice fall speed, native
                // leg. Field order and the 20-fractional-digit value token are
                // contractual and must match wdm6_emit_t2 in
                // ERF_module_mp_wdm6.F90 exactly.
                if (microphysics_debug >= 2 && loop == 0 &&
                    i == diag_i && j == diag_j && (k + 1) >= 85 && (k + 1) <= 100) {
                    auto emit_t2 = [&] (const char* var, Real value) {
                        std::printf("WDM6-DIAG-T2 diag_schema=1 tag=G8 phase=vice"
                                    " source_layer=NATIVE_CPP path_id=cpp expr_id=vice"
                                    " store_id=%s loop=%d i_dbg=%d j_dbg=%d"
                                    " k_dbg=%d k_raw=%d debug_level=%d var=%s value=%+.20E\n",
                                    var, loop, diag_i, diag_j,
                                    k - klo + 1, k, microphysics_debug,
                                    var, static_cast<double>(value));
                    };
                    const Real xni_s = amrex::max(xni_arr(i,j,k), Real(1.0e-30));
                    const Real xmi_t = den_arr(i,j,k) * qi_arr(i,j,k) / xni_s;
                    emit_t2("qi",       qi_arr(i,j,k));
                    emit_t2("xni",      xni_arr(i,j,k));
                    emit_t2("den",      den_arr(i,j,k));
                    emit_t2("xmi",      xmi_t);
                    emit_t2("sqrt_xmi", std::sqrt(xmi_t));
                    emit_t2("dicon",    dicon_loc);
                    emit_t2("diameter", amrex::max(amrex::min(dicon_loc * std::sqrt(xmi_t), dimax_loc),
                                                  wdm6_literal(1.0e-25)));
                    emit_t2("work1c",   work1c);
                }
#endif
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

                // Ice sedimentation. ITER MUST BE 0 HERE. The Fortran calls
                //     nislfv_rain_plmr(...,work1c,denqci,denqci,delqi,dtcld,1,0,0)
                // whose trailing arguments are dt, id, iter, rid -- so id=1 and
                // iter=0. This helper's signature has no id parameter and reads
                // (..., dt, iter, ...), so passing the Fortran's id of 1 in that
                // slot silently set iter=1 and ran the velocity-iteration block
                // that G8 is supposed to skip, recomputing the ICE fall speed
                // from the snow and graupel slope routines. Contrast G5c, which
                // genuinely passes iter=1 (Fortran :1122 ends dtcld,1,1).
                // Confirmed by Tier 2: every input to this call -- qi, xni, den,
                // xmi, dicon, diameter, work1c -- is bitwise equal at every k in
                // the 85..100 window, so the kernel invocation was the only
                // remaining candidate.
                wdm6_nislfv_rain_plm6_column(
                    km, dz.data(), den.data(), denfac.data(), tk.data(),
                    work_ice.data(), qi_col.data(), qi_col.data(), delqi_col, delqi_col, dtcld, 0,
                    pidn0s_loc, pidn0g_loc, Real(qcrmin), Real(alpha_wdm6), Real(n0smax), Real(n0s), wdm6_slope_t0c,
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
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                // NOT COMPARABLE PER-K: fallc is allocated on box2d here (a 2D
                // surface slab, see the fallc_fab definition above) while the
                // Fortran declares fallc(its:ite,kts:kte), a full (i,k) field.
                // Reading it at kdbg>klo is out of bounds. It is emitted at klo
                // so the line stays field-aligned with the Fortran tag, but the
                // second field must be ignored above k=kts until the native path
                // carries a full-column ice fallout. Tracked as a G8 finding.
                std::printf("WDM6-CPP_POST_G8 %3d %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qi_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(fallc_arr(diag_i,diag_j,klo)),
                            static_cast<double>(den_arr(diag_i,diag_j,kdbg) * qi_arr(diag_i,diag_j,kdbg)));
                }
                std::fflush(stdout);
            }
#endif

            // ============================================================
            // Step 3h: G9 Surface precipitation accumulation and sr update
            // ============================================================
#if !defined(AMREX_USE_GPU)
            // G9 is surface-scoped: its inputs are box2d slabs here. Not converted
            // to all-k; see the PRE_G9 note in ERF_module_mp_wdm6.F90.
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
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
            // G9 is surface-scoped: its inputs are box2d slabs here. Not converted
            // to all-k; see the PRE_G9 note in ERF_module_mp_wdm6.F90.
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
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
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G10A %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qi_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(xni_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(cpm_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

            // ================================================================
            // G10a-G10d FUSED, matching the Fortran's single loop nest.
            //
            // ERF_module_mp_wdm6.F90:1420-1421 opens one do k / do i nest and
            // computes supcol (:1422) and xlf (:1423-1424) ONCE per cell; all
            // four phase blocks then reuse them while t is progressively
            // updated by latent heat (:1432 cools, :1459 warms, :1495 warms).
            // This path previously ran the four blocks as four separate
            // ParallelFor kernels, each recomputing supcol and xlf from the
            // already-updated t_arr. That made pgfrz consume a fresher supcol
            // than the Fortran does: measured 1.048787e-06 on supcolt at
            // (107,3,46), amplified to 4.908569e-06 in expterm by the pfrz2
            // factor and inherited by pfrzdtr.
            //
            // xlf is the more dangerous of the two because it is set by a
            // BRANCH on supcol. Melting cools and freezing warms, so t crossing
            // t0c mid-section is exactly what this code does, and a per-kernel
            // re-evaluation would flip the branch where the Fortran does not.
            // Note G10b did not even carry the branch, only xlf = xls - xl.
            //
            // The group-boundary tags stay at their group boundaries per
            // Rule 30; they are emitted per cell from inside the kernel now
            // rather than from a column sweep afterwards. For these purely
            // per-cell blocks that is the same state, and it is strictly closer
            // to the "immediate process-group boundary" the rule asks for.
            // ================================================================
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supcol = t0c - t_arr(i,j,k);
                Real xlf = xls - xl_arr(i,j,k);
                if (supcol < Real(0.0)) xlf = xlf0;

#if !defined(AMREX_USE_GPU)
                const bool diag_cell = (microphysics_debug > 0 && loop == 0 &&
                                        i == diag_i && j == diag_j);
#endif

                // Cloud-ice melt to cloud water (T > t0c)
                if (supcol < Real(0.0) && qi_arr(i,j,k) > Real(0.0)) {
                    const Real qim = qi_arr(i,j,k);  // preserve old ice amount

                    qc_arr(i,j,k)   += qim;                           // qci(:,:,1) += qci(:,:,2)
                    nc_arr(i,j,k)   += xni_arr(i,j,k);                // ncr(:,:,2) += xni
                    t_arr(i,j,k)    -= xlf / cpm_arr(i,j,k) * qim;    // latent heat release
                    qi_arr(i,j,k)    = Real(0.0);                     // zero out ice
                }

#if !defined(AMREX_USE_GPU)
                if (diag_cell) {
                    std::printf("WDM6-CPP_POST_G10A %3d %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                                k + 1,
                                static_cast<double>(qc_arr(i,j,k)),
                                static_cast<double>(qi_arr(i,j,k)),
                                static_cast<double>(nc_arr(i,j,k)),
                                static_cast<double>(xni_arr(i,j,k)),
                                static_cast<double>(t_arr(i,j,k)));
                std::fflush(stdout);
                }
#endif


            // ============================================================
            // Step 3k: G10b Cloud-Water Homogeneous Freezing
            // Exact port of bounded Fortran block (lines 1215-1224)
            // ============================================================
#if !defined(AMREX_USE_GPU)
                if (diag_cell) {
                    std::printf("WDM6-CPP_PRE_G10B %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                                k + 1,
                                static_cast<double>(qc_arr(i,j,k)),
                                static_cast<double>(qi_arr(i,j,k)),
                                static_cast<double>(nc_arr(i,j,k)),
                                static_cast<double>(t_arr(i,j,k)),
                                static_cast<double>(t0c - t_arr(i,j,k)),
                                static_cast<double>(cpm_arr(i,j,k)));
                std::fflush(stdout);
                }
#endif


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

#if !defined(AMREX_USE_GPU)
                if (diag_cell) {
                    std::printf("WDM6-CPP_POST_G10B %3d %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                                k + 1,
                                static_cast<double>(qc_arr(i,j,k)),
                                static_cast<double>(qi_arr(i,j,k)),
                                static_cast<double>(nc_arr(i,j,k)),
                                static_cast<double>(t_arr(i,j,k)),
                                static_cast<double>(t0c - t_arr(i,j,k)));
                std::fflush(stdout);
                }
#endif

            // ============================================================
            // Step 3l: G10c Cloud-Water Heterogeneous Freezing
            // Exact port of bounded Fortran block (pihtf, lines 1241-1258)
            // ============================================================
#if !defined(AMREX_USE_GPU)
                if (diag_cell) {
                    std::printf("WDM6-CPP_PRE_G10C %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                                k + 1,
                                static_cast<double>(qc_arr(i,j,k)),
                                static_cast<double>(qi_arr(i,j,k)),
                                static_cast<double>(nc_arr(i,j,k)),
                                static_cast<double>(t_arr(i,j,k)),
                                static_cast<double>(t0c - t_arr(i,j,k)),
                                static_cast<double>(den_arr(i,j,k)),
                                static_cast<double>(rslopec_arr(i,j,k)),
                                static_cast<double>(cpm_arr(i,j,k)),
                                static_cast<double>(xl_arr(i,j,k)));
                std::fflush(stdout);
                }
#endif

                // Latent heat of FUSION, per the Fortran's per-k prologue:
                //     xlf = xls-xl(i,k) ;  if(supcol.lt.0.) xlf = xlf0
                // shared by G10a-G10d. Freezing releases xlf, not xl.

                // Heterogeneous freezing (Biggs, contact/immersion): 0 > T > -40C
                // Trigger: supcol > 0 (T < T0c) AND qc > qmin
                if (supcol > Real(0.0) && qc_arr(i,j,k) > Real(qmin)) {
                    const Real supcolt = amrex::min(supcol, Real(70.0));
                    const Real expterm = std::exp(pfrz2_loc * supcolt) - Real(1.0);

                    // Use the stored rslopec3, as the Fortran does: in the
                    // rslopecmax branch it is the precomputed rslopec3max
                    // constant, which need not equal rslopecmax cubed.
                    const Real rs3 = rslopec3_arr(i,j,k);

                    // pfrzdtc = min(pi*pi*pfrz1*(exp(pfrz2*supcolt)-1.)*denr/den
                    //               *ncr(:,:,2)*rslopec3*rslopec3/18.*dtcld, qci(:,:,1))
                    // rslopec3 appears TWICE. Using it once left the rate about
                    // 1/rslopec3 too large -- roughly 8e12 for a typical cloud
                    // droplet slope -- so pfrzdtc always saturated at its qc cap
                    // and every cell below freezing lost all its cloud water to
                    // ice in a single substep. Factor order follows the Fortran
                    // left to right so the roundings agree.
                    Real pfrzdtc = pi_wdm6_loc * pi_wdm6_loc * pfrz1_loc * expterm
                                 * denr / den_arr(i,j,k) * nc_arr(i,j,k) * rs3 * rs3
                                 / Real(18.0) * dtcld;
                    pfrzdtc = amrex::min(pfrzdtc, qc_arr(i,j,k));

                    // nfrzdtc = min(pi*pfrz1*(exp(pfrz2*supcolt)-1.)*ncr(:,:,2)
                    //               *rslopec3/6.*dtcld, ncr(:,:,2)) -- one factor here.
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
                    t_arr(i,j,k)  += xlf / cpm_arr(i,j,k) * pfrzdtc;
                    qc_arr(i,j,k) -= pfrzdtc;
                }

#if !defined(AMREX_USE_GPU)
                if (diag_cell) {
                    const Real rslopec_g10c = rslopec_arr(i,j,k);
                    std::printf("WDM6-CPP_POST_G10C %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                                k + 1,
                                static_cast<double>(qc_arr(i,j,k)),
                                static_cast<double>(qi_arr(i,j,k)),
                                static_cast<double>(nc_arr(i,j,k)),
                                static_cast<double>(t_arr(i,j,k)),
                                static_cast<double>(t0c - t_arr(i,j,k)),
                                static_cast<double>(rslopec_g10c * rslopec_g10c * rslopec_g10c));
                std::fflush(stdout);
                }
#endif

            // ============================================================
            // Step 3m: G10d Rain-to-Graupel Freezing
            // Exact port of bounded Fortran block (pgfrz, rain freezing)
            // ============================================================
#if !defined(AMREX_USE_GPU)
                if (diag_cell) {
                    std::printf("WDM6-CPP_PRE_G10D %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                                k + 1,
                                static_cast<double>(qr_arr(i,j,k)),
                                static_cast<double>(qg_arr(i,j,k)),
                                static_cast<double>(nr_arr(i,j,k)),
                                static_cast<double>(t_arr(i,j,k)),
                                static_cast<double>(t0c - t_arr(i,j,k)),
                                static_cast<double>(den_arr(i,j,k)),
                                static_cast<double>(rslope3_arr(i,j,k,0)),
                                static_cast<double>(cpm_arr(i,j,k)),
                                static_cast<double>(xl_arr(i,j,k)));
                std::fflush(stdout);
                }
#endif

                // Latent heat of fusion, as in G10a-G10c above.

                // Rain freezing to graupel: trigger when T < t0c and qr > 0
                if (supcol > Real(0.0) && qr_arr(i,j,k) > Real(0.0)) {
                    const Real supcolt = amrex::min(supcol, Real(70.0));
                    const Real expterm = std::exp(pfrz2_loc * supcolt) - Real(1.0);

                    // Rain slope cubed (from G4/G6 computation)
                    const Real rs3 = rslope3_arr(i,j,k,0);

                    // pfrzdtr = min(140.*(pi*pi)*pfrz1*ncr(:,:,3)*denr/den
                    //               *(exp(pfrz2*supcolt)-1.)*rslope3*rslope3*dtcld, qrs(:,:,1))
                    // Factor order follows the Fortran left to right.
                    Real pfrzdtr = Real(140.0) * (pi_wdm6_loc * pi_wdm6_loc)
                                 * pfrz1_loc * nr_arr(i,j,k)
                                 * denr / den_arr(i,j,k)
                                 * expterm * rs3 * rs3 * dtcld;
                    pfrzdtr = amrex::min(pfrzdtr, qr_arr(i,j,k));

#if !defined(AMREX_USE_GPU)
                    // Tier 2 forensic decomposition of pgfrz, native leg. One
                    // var per line in the WSM6-DIAG-T2 schema; field order and
                    // the 20-fractional-digit value token are contractual and
                    // must match wdm6_emit_t2 in ERF_module_mp_wdm6.F90.
                    // Emitted before the nr update below so the printed nr is
                    // the value pfrzdtr consumed.
                    if (microphysics_debug >= 2 && loop == 0 &&
                        i == diag_i && j == diag_j && (k + 1) >= 40 && (k + 1) <= 50) {
                        auto emit_t2 = [&] (const char* var, Real value) {
                            std::printf("WDM6-DIAG-T2 diag_schema=1 tag=G10D phase=pgfrz"
                                        " source_layer=NATIVE_CPP path_id=cpp expr_id=pgfrz"
                                        " store_id=%s loop=%d i_dbg=%d j_dbg=%d"
                                        " k_dbg=%d k_raw=%d debug_level=%d var=%s value=%+.20E\n",
                                        var, loop, diag_i, diag_j,
                                        k - klo + 1, k, microphysics_debug,
                                        var, static_cast<double>(value));
                        };
                        emit_t2("pi",       pi_wdm6_loc);
                        emit_t2("pisq",     pi_wdm6_loc * pi_wdm6_loc);
                        emit_t2("pfrz1",    pfrz1_loc);
                        emit_t2("pfrz2",    pfrz2_loc);
                        emit_t2("denr",     denr);
                        emit_t2("dtcld",    dtcld);
                        emit_t2("supcolt",  supcolt);
                        emit_t2("expterm",  expterm);
                        emit_t2("rslope3r", rs3);
                        emit_t2("nr",       nr_arr(i,j,k));
                        emit_t2("den",      den_arr(i,j,k));
                        emit_t2("pfrzdtr",  pfrzdtr);
                    }
#endif

                    // Number freezing rate (conditional on nr > nrmin)
                    if (nr_arr(i,j,k) > Real(nrmin)) {
                        Real nfrzdtr = Real(4.0) * pi_wdm6_loc * pfrz1_loc
                                     * nr_arr(i,j,k) * expterm * rs3 * dtcld;
                        nfrzdtr = amrex::min(nfrzdtr, nr_arr(i,j,k));
                        nr_arr(i,j,k) -= nfrzdtr;
                    }

                    // Apply mass and temperature updates
                    qg_arr(i,j,k) += pfrzdtr;
                    t_arr(i,j,k)  += xlf / cpm_arr(i,j,k) * pfrzdtr;
                    qr_arr(i,j,k) -= pfrzdtr;
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G10D %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t0c - t_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,kdbg,0)));
                std::fflush(stdout);
                }
            }
#endif

            // ============================================================
            // Step 3k: G10e Phase cleanup — clamp number concentrations non-negative
            // ============================================================
#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G10E %3d %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                nc_arr(i,j,k) = amrex::max(nc_arr(i,j,k), Real(0.0));
                nr_arr(i,j,k) = amrex::max(nr_arr(i,j,k), Real(0.0));
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G10E %3d %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

            // ============================================================
            // G11: SLOPE3 — Third slope_wdm6 call + avedia/rslopec recompute
            // ============================================================

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G11 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(den_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(p_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(xl_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qsatw_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qsati_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
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
                                     wdm6_slope_t0c, Real(qcrmin),
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
                    // Same expression as G3; share the one exact port. The
                    // previous wdm6_lamdac() here differed from the Fortran on
                    // three counts: it cancelled den out of the argument, it
                    // used std::pow with a true 1/3 rather than exp/log with the
                    // unsuffixed .33333333, and it carried a qc/nc guard of its
                    // own that fired inside this block's own gate.
                    const Real rslc = wdm6_rslopec_exact(qci_for_lamdac, den_arr(i,j,k),
                                                         nc_for_lamdac, pidnc_loc);
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
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G12 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(xl_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(p_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(den_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qsatw_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qsati_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
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
                                static_cast<double>(wdm6_literal(1.496e-6)
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
                                                     * (wdm6_literal(8.794e-5)
                                                        * std::exp(std::log(t_arr(i,j,k)) * wdm6_literal(1.81))
                                                        / p_arr(i,j,k)))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT38 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.0)/qsatw_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT39 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(Real(1.0)/(wdm6_literal(8.794e-5)
                                                       * std::exp(std::log(t_arr(i,j,k)) * wdm6_literal(1.81))
                                                       / p_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT40 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>((Real(1.0)/qsatw_arr(i,j,k))
                                                   *(Real(1.0)/(wdm6_literal(8.794e-5)
                                                      * std::exp(std::log(t_arr(i,j,k)) * wdm6_literal(1.81))
                                                      / p_arr(i,j,k)))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT41 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(qsatw_arr(i,j,k)
                                                   * (wdm6_literal(8.794e-5)
                                                      * std::exp(std::log(t_arr(i,j,k)) * wdm6_literal(1.81))
                                                      / p_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT37 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(den_arr(i,j,k) * xl_arr(i,j,k) * xl_arr(i,j,k)
                                                   /(Real(1.414e3) * wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                     * den_arr(i,j,k) * Real(rv)
                                                     * t_arr(i,j,k) * t_arr(i,j,k))
                                                   + Real(1.0)
                                                   /(qsatw_arr(i,j,k)
                                                     * (wdm6_literal(8.794e-5)
                                                        * std::exp(std::log(t_arr(i,j,k)) * wdm6_literal(1.81))
                                                        / p_arr(i,j,k)))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT2 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(wdm6_literal(8.794e-5) * std::exp(std::log(t_arr(i,j,k)) * wdm6_literal(1.81)) / p_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT24 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::log(t_arr(i,j,k))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT25 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::log(t_arr(i,j,k)) * wdm6_literal(1.81)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT26 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::exp(std::log(t_arr(i,j,k)) * wdm6_literal(1.81))));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT27 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(wdm6_literal(8.794e-5) * std::exp(std::log(t_arr(i,j,k)) * wdm6_literal(1.81))
                                                   / p_arr(i,j,k)));
                    std::printf("WDM6-CPP_G11V_DIFFAC_INT28 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(wdm6_literal(8.794e-5) * std::exp(std::log(t_arr(i,j,k)) * wdm6_literal(1.81))
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
                                                   * wdm6_literal(0.3333333)));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT5 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::exp(std::log(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   / wdm6_diffus(t_arr(i,j,k), p_arr(i,j,k)))
                                                   * wdm6_literal(0.3333333))));
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
                                                   * wdm6_literal(0.3333333))
                                                   / std::sqrt(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k)))));
                    std::printf("WDM6-CPP_G11V_VENFAC_INT11 %3d %24.16E\n",
                                diag_k + 1,
                                static_cast<double>(std::exp(std::log(wdm6_viscos(t_arr(i,j,k), den_arr(i,j,k))
                                                   / wdm6_diffus(t_arr(i,j,k), p_arr(i,j,k)))
                                                   * wdm6_literal(0.3333333))
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
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G12 %3d %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(work2_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G11 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(rslope_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(avedia_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslopec_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(rslopec2_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(avedia_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(work2_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslope_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslopeb_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,kdbg,1)),
                            static_cast<double>(rslope3_arr(diag_i,diag_j,kdbg,2)),
                            static_cast<double>(rslopec3_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,1)));
                std::fflush(stdout);
                }
            }
#endif

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G13A %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(rslopec3_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(work1_arr(diag_i,diag_j,kdbg,0)),
                            static_cast<double>(work2_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supsat = amrex::max(qv_arr(i,j,k), Real(qmin)) - qsatw_arr(i,j,k);
                const Real satdt = supsat / dtcld;
                const Real lencon = Real(2.7e-2) * den_arr(i,j,k) * qc_arr(i,j,k)
                    * (Real(1.0e20/16.0) * rslopec2_arr(i,j,k) * rslopec2_arr(i,j,k) - Real(0.4));
                const Real lenconcr = amrex::max(Real(1.2) * lencon, Real(qcrmin));

                if (qc_arr(i,j,k) > qcr_arr(i,j,k) && nc_arr(i,j,k) > Real(ncmin)) {
                    // Fortran :1889
                    //   praut(i,k) = qck1*qci(i,k,1)**(7./3.)*ncr(i,k,2)**(-1./3.)
                    // Both exponents are unsuffixed, so they are float32 values
                    // widened, not true thirds: 7./3. is 2.3333332538604736
                    // against 2.3333333333333335, and 1./3. is 0.3333333432674408
                    // against 0.33333333333333331. Exponents amplify: the error
                    // enters as ln(base)*delta, and with qc about 1.6e-02 and nc
                    // about 1.8e+08 those logs are -4.2 and +19.0, which is why a
                    // 3e-08 exponent difference produced 1.86e-07 on praut.
                    praut_arr(i,j,k) = qck1_loc * std::pow(qc_arr(i,j,k), wdm6_literal(7.0/3.0))
                        * std::pow(nc_arr(i,j,k), wdm6_literal(-1.0/3.0));
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
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G13A %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(praut_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nrauto_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(pracw_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nraccr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nccol_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nrcol_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(prevp_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg)
                std::printf("WDM6-CPP_PRE_G13B %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qi_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(xni_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supcol = Real(t0c) - t_arr(i,j,k);
                const Real n0sfac = amrex::max(
                    amrex::min(std::exp(Real(alpha_wdm6) * supcol),
                               Real(n0smax) / Real(n0s)),
                    Real(1.0));
                const Real supsat = amrex::max(qv_arr(i,j,k), Real(qmin)) - qsati_arr(i,j,k);
                amrex::ignore_unused(supsat);
                const Real satdt = supsat / dtcld;
                amrex::ignore_unused(satdt);

                const Real qi_val = qi_arr(i,j,k);
                if (!(supcol > Real(0.0) && qi_val > Real(qmin))) {
                    return;
                }

                Real temp = den_arr(i,j,k) * amrex::max(qi_val, Real(qmin));
                temp = std::sqrt(std::sqrt(temp * temp * temp));
                xni_arr(i,j,k) = amrex::min(amrex::max(Real(5.38e7) * temp, Real(1.e3)), Real(1.e6));

                // Fortran :2021 eacrs = exp(0.07*(-supcol)); 0.07 is unsuffixed
                // so it carries the float assumption. Proven at Tier 2 on the
                // G13B psaci decomposition: eacrs diverged 1.03e-08 to 1.14e-08
                // relative across k=89..95, matching delta(0.07)*(-supcol).
                const Real eacrs = std::exp(wdm6_literal(0.07) * (-supcol));
                const Real xni_safe = amrex::max(xni_arr(i,j,k), Real(1.0e-30));
                const Real xmi = den_arr(i,j,k) * qi_val / xni_safe;
                const Real diameter = amrex::min(Real(dicon) * std::sqrt(xmi), Real(dimax));
                // Fortran :2025 vt2i = 1.49e4*diameter**1.31. 1.49e4 is exactly
                // representable in float32 so it needs no routing, but the 1.31
                // EXPONENT does, and exponents amplify: the error enters as
                // ln(diameter)*delta, and with diameter at the dimax cap of
                // 5e-4 that is ln(5e-4) = -7.60 times -5.722e-08, i.e. 4.349e-07.
                // Tier 2 measured exactly 4.349270e-07 on vt2i across k=91..95.
                // It then amplifies again through the |vt2ave - vt2i|
                // cancellation in psaci: at k=89 vt2diff is 0.139 against vt2i
                // 0.681, a 4.90x factor giving the observed 2.14e-06.
                const Real vt2i = Real(1.49e4) * std::pow(diameter, wdm6_literal(1.31));
                const Real vt2r = pvtr_loc * rslopeb_arr(i,j,k,0) * denfac_arr(i,j,k);
                const Real vt2s = pvts_loc * rslopeb_arr(i,j,k,1) * denfac_arr(i,j,k);
                const Real vt2g = pvtg_loc * rslopeb_arr(i,j,k,2) * denfac_arr(i,j,k);
                const Real qsum = amrex::max(qs_arr(i,j,k) + qg_arr(i,j,k), Real(1.0e-15));
                const Real vt2ave = (qsum > Real(1.0e-15))
                    ? (vt2s * qs_arr(i,j,k) + vt2g * qg_arr(i,j,k)) / qsum
                    : Real(0.0);

                if (qr_arr(i,j,k) > Real(qcrmin)) {
                    const Real acrfac = Real(6.0) * rslope2_arr(i,j,k,0)
                        + Real(4.0) * diameter * rslope_arr(i,j,k,0)
                        + diameter * diameter;
                    praci_arr(i,j,k) = pi_wdm6_loc * qi_val * nr_arr(i,j,k)
                        * std::abs(vt2r - vt2i) * acrfac / Real(4.0);
                    praci_arr(i,j,k) *= std::pow(
                        amrex::min(amrex::max(Real(0.0), qr_arr(i,j,k) / qi_val), Real(1.0)),
                        Real(2.0));
                    praci_arr(i,j,k) = amrex::min(praci_arr(i,j,k), qi_val / dtcld);

                    piacr_arr(i,j,k) = pi_wdm6_loc * pi_wdm6_loc * Real(WDM6::avtr)
                        * nr_arr(i,j,k) * Real(denr) * xni_arr(i,j,k) * denfac_arr(i,j,k)
                        * g7pbr_loc * rslope3_arr(i,j,k,0) * rslope2_arr(i,j,k,0)
                        * rslopeb_arr(i,j,k,0) / (Real(24.0) * den_arr(i,j,k));
                    piacr_arr(i,j,k) *= std::pow(
                        amrex::min(amrex::max(Real(0.0), qi_val / qr_arr(i,j,k)), Real(1.0)),
                        Real(2.0));
                    piacr_arr(i,j,k) = amrex::min(piacr_arr(i,j,k), qr_arr(i,j,k) / dtcld);
                }

                if (nr_arr(i,j,k) > Real(nrmin)) {
                    niacr_arr(i,j,k) = pi_wdm6_loc * Real(WDM6::avtr) * nr_arr(i,j,k)
                        * xni_arr(i,j,k) * denfac_arr(i,j,k) * g4pbr_loc
                        * rslope2_arr(i,j,k,0) * rslopeb_arr(i,j,k,0) / Real(4.0);
                    niacr_arr(i,j,k) *= std::pow(
                        amrex::min(amrex::max(Real(0.0), qi_val / qr_arr(i,j,k)), Real(1.0)),
                        Real(2.0));
                    niacr_arr(i,j,k) = amrex::min(niacr_arr(i,j,k), nr_arr(i,j,k) / dtcld);
                }

                if (qs_arr(i,j,k) > Real(qcrmin)) {
                    const Real acrfac = Real(2.0) * rslope3_arr(i,j,k,1)
                        + Real(2.0) * diameter * rslope2_arr(i,j,k,1)
                        + diameter * diameter * rslope_arr(i,j,k,1);
                    psaci_arr(i,j,k) = pi_wdm6_loc * qi_val * eacrs * Real(n0s) * n0sfac
                        * std::abs(vt2ave - vt2i) * acrfac / Real(4.0);
                    psaci_arr(i,j,k) = amrex::min(psaci_arr(i,j,k), qi_val / dtcld);

#if !defined(AMREX_USE_GPU)
                    // Tier 2 forensic decomposition of psaci, native leg. One
                    // var per line in the WSM6-DIAG-T2 schema; field order and
                    // the 20-fractional-digit value token are contractual and
                    // must match wdm6_emit_t2 in ERF_module_mp_wdm6.F90
                    // exactly. Same shape as the G13E pidep decomposition.
                    if (microphysics_debug >= 2 && loop == 0 &&
                        i == diag_i && j == diag_j && (k + 1) >= 85 && (k + 1) <= 95) {
                        auto emit_t2 = [&] (const char* var, Real value) {
                            std::printf("WDM6-DIAG-T2 diag_schema=1 tag=G13B phase=psaci"
                                        " source_layer=NATIVE_CPP path_id=cpp expr_id=psaci"
                                        " store_id=%s loop=%d i_dbg=%d j_dbg=%d"
                                        " k_dbg=%d k_raw=%d debug_level=%d var=%s value=%+.20E\n",
                                        var, loop, diag_i, diag_j,
                                        k - klo + 1, k, microphysics_debug,
                                        var, static_cast<double>(value));
                        };
                        emit_t2("eacrs",    eacrs);
                        emit_t2("diameter", diameter);
                        emit_t2("vt2i",     vt2i);
                        emit_t2("vt2ave",   vt2ave);
                        emit_t2("vt2diff",  std::abs(vt2ave - vt2i));
                        emit_t2("acrfac",   acrfac);
                        emit_t2("n0sfac",   n0sfac);
                        emit_t2("psaci",    psaci_arr(i,j,k));
                    }
#endif
                }

                if (qg_arr(i,j,k) > Real(qcrmin)) {
                    // Fortran :2114, the same unsuffixed 0.07 as eacrs above.
                    // pgaci reads bitwise-equal in the current column only
                    // because the branch is inactive there; routed for the same
                    // reason, not on separate evidence.
                    const Real egi = std::exp(wdm6_literal(0.07) * (-supcol));
                    const Real acrfac = Real(2.0) * rslope3_arr(i,j,k,2)
                        + Real(2.0) * diameter * rslope2_arr(i,j,k,2)
                        + diameter * diameter * rslope_arr(i,j,k,2);
                    pgaci_arr(i,j,k) = pi_wdm6_loc * egi * qi_val * n0g_loc
                        * std::abs(vt2ave - vt2i) * acrfac / Real(4.0);
                    pgaci_arr(i,j,k) = amrex::min(pgaci_arr(i,j,k), qi_val / dtcld);
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg)
                std::printf("WDM6-CPP_POST_G13B %3d %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(praci_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(piacr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(niacr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(psaci_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(pgaci_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
            }
#endif

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G13C %3d %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(denfac_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supcol = Real(t0c) - t_arr(i,j,k);
                const Real n0sfac = amrex::max(
                    amrex::min(std::exp(Real(alpha_wdm6) * supcol),
                               Real(n0smax) / Real(n0s)),
                    Real(1.0));
                const Real qs_val = qs_arr(i,j,k);
                const Real qg_val = qg_arr(i,j,k);
                const Real qc_val = qc_arr(i,j,k);
                const Real nc_val = nc_arr(i,j,k);

                // The Fortran (:2143, :2152) writes these ratios UNGUARDED:
                //     min(max(0.0, qrs(i,k,2)/qci(i,k,1)), 1.)
                // so when qc == 0 and qs > 0 the division yields IEEE +Inf,
                // max(0,Inf) is Inf, and the clamp returns 1.0. Returning 0.0
                // for that case instead -- as a defensive divide-by-zero guard
                // naturally does -- lands on the OPPOSITE end of the clamp and
                // silently changes the physics.
                //
                // It hides in psacw, which is gated on qc > qmin and so never
                // sees qc == 0, but nsacw is gated on nc > ncmin, which passes
                // with qc == 0, and there the ratio is load-bearing. Measured
                // at (199,3,88) from a bitwise-identical PRE_G13C with qc = 0,
                // qs = 1.935e-02 and nc = 1.048e+06: the bridge produced
                // nsacw = 9.6929949510183156e+04 and the native leg exactly 0.
                // That single term is the entire step-6 nn discontinuity.
                //
                // The division is still avoided rather than allowed to produce
                // Inf, so builds with FP-exception trapping stay clean. The
                // qc == 0 and qs == 0 case cannot be reached by any consumer --
                // every use is gated on qs > qcrmin, qg > qcrmin or qc > qmin --
                // and 0.0 is returned there only for definiteness.
                const Real ratio_s = (qc_val > Real(0.0))
                    ? amrex::min(amrex::max(Real(0.0), qs_val / qc_val), Real(1.0))
                    : (qs_val > Real(0.0) ? Real(1.0) : Real(0.0));
                // Same treatment as ratio_s above, and for the same reason. This
                // was MISSED when ratio_s was fixed: that edit replaced only the
                // ratio_s else-branch while its comment claimed both, so ngacw
                // kept returning 0 where the Fortran returns 1.0. Surfaced by the
                // step-8 retreat as POST_G13C field 4 at relative error exactly
                // 1.0, bridge 8.85033631487811001e-02 against native 0, which was
                // the whole remaining Milestone B nn residual.
                const Real ratio_g = (qc_val > Real(0.0))
                    ? amrex::min(amrex::max(Real(0.0), qg_val / qc_val), Real(1.0))
                    : (qg_val > Real(0.0) ? Real(1.0) : Real(0.0));

                if (qs_val > Real(qcrmin) && qc_val > Real(qmin)) {
                    psacw_arr(i,j,k) = amrex::min(
                        pacrc_loc * n0sfac * rslope3_arr(i,j,k,1) * rslopeb_arr(i,j,k,1)
                            * ratio_s * ratio_s * qc_val * denfac_arr(i,j,k),
                        qc_val / dtcld);
                }
                if (qs_val > Real(qcrmin) && nc_val > Real(ncmin)) {
                    nsacw_arr(i,j,k) = amrex::min(
                        pacrc_loc * n0sfac * rslope3_arr(i,j,k,1) * rslopeb_arr(i,j,k,1)
                            * ratio_s * ratio_s * nc_val * denfac_arr(i,j,k),
                        nc_val / dtcld);
                }
                if (qg_val > Real(qcrmin) && qc_val > Real(qmin)) {
                    pgacw_arr(i,j,k) = amrex::min(
                        pacrg_loc * rslope3_arr(i,j,k,2) * rslopeb_arr(i,j,k,2)
                            * qc_val * ratio_g * ratio_g * denfac_arr(i,j,k),
                        qc_val / dtcld);
                }
                if (qg_val > Real(qcrmin) && nc_val > Real(ncmin)) {
                    ngacw_arr(i,j,k) = amrex::min(
                        pacrg_loc * rslope3_arr(i,j,k,2) * rslopeb_arr(i,j,k,2)
                            * nc_val * ratio_g * ratio_g * denfac_arr(i,j,k),
                        nc_val / dtcld);
                }

                const Real qsum = amrex::max(qs_val + qg_val, Real(1.0e-15));
                if (qsum > Real(1.0e-15)) {
                    paacw_arr(i,j,k) = (qs_val * psacw_arr(i,j,k) + qg_val * pgacw_arr(i,j,k)) / qsum;
                    naacw_arr(i,j,k) = (qs_val * nsacw_arr(i,j,k) + qg_val * ngacw_arr(i,j,k)) / qsum;
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G13C %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(psacw_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nsacw_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(pgacw_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(ngacw_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(paacw_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(naacw_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                const Real diag_supcol = Real(t0c) - t_arr(diag_i,diag_j,kdbg);
                std::printf("WDM6-CPP_PRE_G13D %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(paacw_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(diag_supcol));
                std::fflush(stdout);
                }
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supcol = Real(t0c) - t_arr(i,j,k);
                const Real n0sfac = amrex::max(
                    amrex::min(std::exp(Real(alpha_wdm6) * supcol),
                               Real(n0smax) / Real(n0s)),
                    Real(1.0));
                const Real qr_val = qr_arr(i,j,k);
                const Real qs_val = qs_arr(i,j,k);
                const Real qg_val = qg_arr(i,j,k);
                const Real nr_val = nr_arr(i,j,k);
                const Real vt2r = pvtr_loc * rslopeb_arr(i,j,k,0) * denfac_arr(i,j,k);
                const Real vt2s = pvts_loc * rslopeb_arr(i,j,k,1) * denfac_arr(i,j,k);
                const Real vt2g = pvtg_loc * rslopeb_arr(i,j,k,2) * denfac_arr(i,j,k);
                const Real qsum = qs_val + qg_val;
                const Real vt2ave = (qsum > Real(0.0))
                    ? (vt2s * qs_val + vt2g * qg_val) / qsum
                    : Real(0.0);

                if (qs_val > Real(qcrmin) && qr_val > Real(qcrmin)) {
                    if (supcol > Real(0.0)) {
                        const Real acrfac =
                            Real(5.0) * rslope3_arr(i,j,k,1) * rslope3_arr(i,j,k,1)
                          + Real(4.0) * rslope3_arr(i,j,k,1) * rslope2_arr(i,j,k,1)
                                * rslope_arr(i,j,k,0)
                          + Real(1.5) * rslope2_arr(i,j,k,1) * rslope2_arr(i,j,k,1)
                                * rslope2_arr(i,j,k,0);
                        pracs_arr(i,j,k) = pi_wdm6_loc * pi_wdm6_loc * nr_val * Real(n0s)
                            * n0sfac * std::abs(vt2r - vt2ave)
                            * (Real(dens) / den_arr(i,j,k)) * acrfac;
                        const Real ratio = amrex::min(
                            amrex::max(Real(0.0), qr_val / qs_val), Real(1.0));
                        pracs_arr(i,j,k) *= ratio * ratio;
                        pracs_arr(i,j,k) = amrex::min(pracs_arr(i,j,k), qs_val / dtcld);
                    }

                    const Real acrfac =
                        Real(30.0) * rslope3_arr(i,j,k,0) * rslope2_arr(i,j,k,0)
                            * rslope_arr(i,j,k,1)
                      + Real(10.0) * rslope2_arr(i,j,k,0) * rslope2_arr(i,j,k,0)
                            * rslope2_arr(i,j,k,1)
                      + Real(2.0) * rslope3_arr(i,j,k,0) * rslope3_arr(i,j,k,1);
                    psacr_arr(i,j,k) = pi_wdm6_loc * pi_wdm6_loc * nr_val * Real(n0s)
                        * n0sfac * std::abs(vt2ave - vt2r)
                        * (Real(denr) / den_arr(i,j,k)) * acrfac;
                    const Real ratio = amrex::min(
                        amrex::max(Real(0.0), qs_val / qr_val), Real(1.0));
                    psacr_arr(i,j,k) *= ratio * ratio;
                    psacr_arr(i,j,k) = amrex::min(psacr_arr(i,j,k), qr_val / dtcld);
                }

                if (qs_val > Real(qcrmin) && nr_val > Real(nrmin)) {
                    const Real acrfac =
                        Real(1.5) * rslope2_arr(i,j,k,0) * rslope_arr(i,j,k,1)
                      + rslope_arr(i,j,k,0) * rslope2_arr(i,j,k,1)
                      + Real(0.5) * rslope3_arr(i,j,k,1);
                    nsacr_arr(i,j,k) = pi_wdm6_loc * nr_val * Real(n0s) * n0sfac
                        * std::abs(vt2ave - vt2r) * acrfac;
                    const Real ratio = amrex::min(
                        amrex::max(Real(0.0), qs_val / qr_val), Real(1.0));
                    nsacr_arr(i,j,k) *= ratio * ratio;
                    nsacr_arr(i,j,k) = amrex::min(nsacr_arr(i,j,k), nr_val / dtcld);
                }

                if (qg_val > Real(qcrmin) && qr_val > Real(qcrmin)) {
                    const Real acrfac =
                        Real(30.0) * rslope3_arr(i,j,k,0) * rslope2_arr(i,j,k,0)
                            * rslope_arr(i,j,k,2)
                      + Real(10.0) * rslope2_arr(i,j,k,0) * rslope2_arr(i,j,k,0)
                            * rslope2_arr(i,j,k,2)
                      + Real(2.0) * rslope3_arr(i,j,k,0) * rslope3_arr(i,j,k,2);
                    pgacr_arr(i,j,k) = pi_wdm6_loc * pi_wdm6_loc * nr_val * n0g_loc
                        * std::abs(vt2ave - vt2r) * (Real(denr) / den_arr(i,j,k))
                        * acrfac;
                    const Real ratio = amrex::min(
                        amrex::max(Real(0.0), qg_val / qr_val), Real(1.0));
                    pgacr_arr(i,j,k) *= ratio * ratio;
                    pgacr_arr(i,j,k) = amrex::min(pgacr_arr(i,j,k), qr_val / dtcld);
                }

                if (qg_val > Real(qcrmin) && nr_val > Real(nrmin)) {
                    const Real acrfac =
                        Real(1.5) * rslope2_arr(i,j,k,0) * rslope_arr(i,j,k,2)
                      + rslope_arr(i,j,k,0) * rslope2_arr(i,j,k,2)
                      + Real(0.5) * rslope3_arr(i,j,k,2);
                    ngacr_arr(i,j,k) = pi_wdm6_loc * nr_val * n0g_loc
                        * std::abs(vt2ave - vt2r) * acrfac;
                    const Real ratio = amrex::min(
                        amrex::max(Real(0.0), qg_val / qr_val), Real(1.0));
                    ngacr_arr(i,j,k) *= ratio * ratio;
                    ngacr_arr(i,j,k) = amrex::min(ngacr_arr(i,j,k), nr_val / dtcld);
                }

                if (qg_val > Real(qcrmin) && qs_val > Real(qcrmin)) {
                    pgacs_arr(i,j,k) = Real(0.0);
                }

                if (supcol <= Real(0.0)) {
                    const Real xlf = Real(xlf0);
                    if (qs_val > Real(0.0)) {
                        pseml_arr(i,j,k) = amrex::min(
                            amrex::max(
                                Real(cliq) * supcol
                                    * (paacw_arr(i,j,k) + psacr_arr(i,j,k)) / xlf,
                                -qs_val / dtcld),
                            Real(0.0));
                    }
                    if (qs_val > Real(qcrmin)) {
                        const Real sfac = rslope_arr(i,j,k,1) * Real(n0s) * n0sfac / qs_val;
                        nseml_arr(i,j,k) = -sfac * pseml_arr(i,j,k);
                    }

                    if (qg_val > Real(0.0)) {
                        pgeml_arr(i,j,k) = amrex::min(
                            amrex::max(
                                Real(cliq) * supcol
                                    * (paacw_arr(i,j,k) + pgacr_arr(i,j,k)) / xlf,
                                -qg_val / dtcld),
                            Real(0.0));
                    }
                    if (qg_val > Real(qcrmin)) {
                        const Real gfac = rslope_arr(i,j,k,2) * n0g_loc / qg_val;
                        ngeml_arr(i,j,k) = -gfac * pgeml_arr(i,j,k);
                    }
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G13D %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(pracs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(psacr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nsacr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(pgacr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(ngacr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(pgacs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(pseml_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nseml_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(pgeml_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(ngeml_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                // Derived terms must be recomputed per level, not hoisted.
                const Real diag_supcol = Real(t0c) - t_arr(diag_i,diag_j,kdbg);
                const Real diag_supsat =
                    amrex::max(qv_arr(diag_i,diag_j,kdbg), Real(qmin))
                    - qsati_arr(diag_i,diag_j,kdbg);
                std::printf("WDM6-CPP_PRE_G13E %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qi_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(prevp_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(rhi_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(diag_supcol),
                            static_cast<double>(diag_supsat));
                }
                std::fflush(stdout);
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supcol = Real(t0c) - t_arr(i,j,k);
                if (supcol <= Real(0.0)) {
                    return;
                }

                const Real n0sfac = amrex::max(
                    amrex::min(std::exp(Real(alpha_wdm6) * supcol),
                               Real(n0smax) / Real(n0s)),
                    Real(1.0));
                const Real supsat =
                    amrex::max(qv_arr(i,j,k), Real(qmin)) - qsati_arr(i,j,k);
                const Real satdt = supsat / dtcld;
                int ifsat = 0;

                const Real qi_val = qi_arr(i,j,k);
                const Real qs_val = qs_arr(i,j,k);
                const Real qg_val = qg_arr(i,j,k);
                const Real xni = xni_arr(i,j,k);
                const Real xni_safe = amrex::max(xni, Real(1.0e-30));
                const Real xmi = den_arr(i,j,k) * qi_val / xni_safe;
                const Real diameter = amrex::min(
                    Real(dicon) * std::sqrt(xmi), Real(dimax));
                const Real rhi = rhi_arr(i,j,k);
                const Real work1i = work1_arr(i,j,k,1);

                if (qi_val > Real(0.0) && ifsat != 1) {
                    pidep_arr(i,j,k) = Real(4.0) * diameter * xni
                        * (rhi - Real(1.0)) / work1i;
                    Real supice = satdt - prevp_arr(i,j,k);
                    if (pidep_arr(i,j,k) < Real(0.0)) {
                        pidep_arr(i,j,k) = amrex::max(
                            amrex::max(pidep_arr(i,j,k), satdt / Real(2.0)),
                            supice);
                        pidep_arr(i,j,k) = amrex::max(
                            pidep_arr(i,j,k), -qi_val / dtcld);
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
#if !defined(AMREX_USE_GPU)
                // Tier 2 forensic decomposition of pidep, native leg. One var
                // per line in the WSM6-DIAG-T2 schema; field order and the
                // 20-fractional-digit value token are contractual and must
                // match wdm6_emit_t2 in ERF_module_mp_wdm6.F90 exactly.
                if (microphysics_debug >= 2 && loop == 0 &&
                    i == diag_i && j == diag_j && (k + 1) >= 88 && (k + 1) <= 100) {
                    auto emit_t2 = [&] (const char* var, Real value) {
                        std::printf("WDM6-DIAG-T2 diag_schema=1 tag=G13E phase=pidep"
                                    " source_layer=NATIVE_CPP path_id=cpp expr_id=pidep"
                                    " store_id=%s loop=%d i_dbg=%d j_dbg=%d"
                                    " k_dbg=%d k_raw=%d debug_level=%d var=%s value=%+.20E\n",
                                    var, loop, diag_i, diag_j,
                                    k - klo + 1, k, microphysics_debug,
                                    var, static_cast<double>(value));
                    };
                    emit_t2("supcol",    supcol);
                    emit_t2("qi",        qi_val);
                    emit_t2("xni",       xni);
                    emit_t2("xmi",       xmi);
                    emit_t2("diameter",  diameter);
                    emit_t2("rhi",       rhi);
                    emit_t2("work1i",    work1i);
                    emit_t2("supsat",    supsat);
                    emit_t2("satdt",     satdt);
                    emit_t2("prevp",     prevp_arr(i,j,k));
                    emit_t2("pidep_raw", Real(4.0) * diameter * xni
                                         * (rhi - Real(1.0)) / work1i);
                    emit_t2("supice",    satdt - prevp_arr(i,j,k));
                    emit_t2("pidep",     pidep_arr(i,j,k));
                    emit_t2("ifsat",     Real(ifsat));
                    std::fflush(stdout);
                }
#endif

                if (qs_val > Real(0.0) && ifsat != 1) {
                    const Real coeres_s = rslope2_arr(i,j,k,1)
                        * std::sqrt(rslope_arr(i,j,k,1) * rslopeb_arr(i,j,k,1));
                    psdep_arr(i,j,k) = (rhi - Real(1.0)) * n0sfac
                        * (precs1_loc * rslope2_arr(i,j,k,1)
                           + precs2_loc * work2_arr(i,j,k) * coeres_s)
                        / work1i;
                    Real supice = satdt - prevp_arr(i,j,k) - pidep_arr(i,j,k);
                    if (psdep_arr(i,j,k) < Real(0.0)) {
                        psdep_arr(i,j,k) = amrex::max(
                            psdep_arr(i,j,k), -qs_val / dtcld);
                        psdep_arr(i,j,k) = amrex::max(
                            amrex::max(psdep_arr(i,j,k), satdt / Real(2.0)),
                            supice);
                    } else {
                        psdep_arr(i,j,k) = amrex::min(
                            amrex::min(psdep_arr(i,j,k), satdt / Real(2.0)),
                            supice);
                    }
                    if (std::abs(prevp_arr(i,j,k) + pidep_arr(i,j,k)
                                 + psdep_arr(i,j,k)) >= std::abs(satdt)) {
                        ifsat = 1;
                    }
#if !defined(AMREX_USE_GPU)
                    // Tier 2 forensic decomposition of psdep, native leg. Must
                    // match wdm6_emit_t2 in ERF_module_mp_wdm6.F90 exactly in
                    // field order and value token.
                    if (microphysics_debug >= 2 && loop == 0 &&
                        i == diag_i && j == diag_j && (k + 1) >= 85 && (k + 1) <= 98) {
                        auto emit_t2 = [&] (const char* var, Real value) {
                            std::printf("WDM6-DIAG-T2 diag_schema=1 tag=G13E phase=psdep"
                                        " source_layer=NATIVE_CPP path_id=cpp expr_id=psdep"
                                        " store_id=%s loop=%d i_dbg=%d j_dbg=%d"
                                        " k_dbg=%d k_raw=%d debug_level=%d var=%s value=%+.20E\n",
                                        var, loop, diag_i, diag_j,
                                        k - klo + 1, k, microphysics_debug,
                                        var, static_cast<double>(value));
                        };
                        emit_t2("precs1", precs1_loc);
                        emit_t2("precs2", precs2_loc);
                        emit_t2("coeres", coeres_s);
                        emit_t2("term1",  precs1_loc * rslope2_arr(i,j,k,1));
                        emit_t2("term2",  precs2_loc * work2_arr(i,j,k) * coeres_s);
                        emit_t2("rhm1",   rhi - Real(1.0));
                        emit_t2("work1s", work1i);
                        emit_t2("psdep",  psdep_arr(i,j,k));
                    }
#endif
                }

                if (qg_val > Real(0.0) && ifsat != 1) {
                    const Real coeres_g = rslope2_arr(i,j,k,2)
                        * std::sqrt(rslope_arr(i,j,k,2) * rslopeb_arr(i,j,k,2));
                    pgdep_arr(i,j,k) = (rhi - Real(1.0))
                        * (precg1_loc * rslope2_arr(i,j,k,2)
                           + precg2_loc * work2_arr(i,j,k) * coeres_g)
                        / work1i;
                    Real supice = satdt - prevp_arr(i,j,k) - pidep_arr(i,j,k)
                                - psdep_arr(i,j,k);
                    if (pgdep_arr(i,j,k) < Real(0.0)) {
                        pgdep_arr(i,j,k) = amrex::max(
                            pgdep_arr(i,j,k), -qg_val / dtcld);
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
                    const Real supice = satdt - prevp_arr(i,j,k) - pidep_arr(i,j,k)
                                      - psdep_arr(i,j,k) - pgdep_arr(i,j,k);
                    const Real xni0 = wdm6_literal(1.0e3) * std::exp(wdm6_literal(0.1) * supcol);
                    const Real roqi0 = wdm6_literal(4.92e-11) * std::pow(xni0, wdm6_literal(1.33));
                    const Real pigen_raw = amrex::max(
                        Real(0.0),
                        (roqi0 / den_arr(i,j,k) - amrex::max(qi_val, Real(0.0))) / dtcld);
                    pigen_arr(i,j,k) = amrex::min(
                        amrex::min(pigen_raw, satdt), supice);
                    if (!amrex::Gpu::inLaunchRegion() && microphysics_debug >= 1 &&
                        i == diag_i && j == diag_j && k >= 39 && k <= 41) {
                        std::printf("WDM6-CPP_T2_G13E_PIGEN_EXACT k=%3d supsat=%24.16E satdt=%24.16E supice=%24.16E supcol=%24.16E xni0=%24.16E roqi0=%24.16E den=%24.16E qi=%24.16E\n",
                                    k + 1,
                                    static_cast<double>(supsat),
                                    static_cast<double>(satdt),
                                    static_cast<double>(supice),
                                    static_cast<double>(supcol),
                                    static_cast<double>(xni0),
                                    static_cast<double>(roqi0),
                                    static_cast<double>(den_arr(i,j,k)),
                                    static_cast<double>(amrex::max(qi_val, Real(0.0))));
                        std::printf("WDM6-CPP_T2_G13E_PIGEN_POST k=%3d pigen_raw=%24.16E satdt=%24.16E supice=%24.16E pigen=%24.16E\n",
                                    k + 1,
                                    static_cast<double>(pigen_raw),
                                    static_cast<double>(satdt),
                                    static_cast<double>(supice),
                                    static_cast<double>(pigen_arr(i,j,k)));
                    }
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg)
                std::printf("WDM6-CPP_POST_G13E %3d %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(pidep_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(psdep_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(pgdep_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(pigen_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
            }
#endif

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg)
                std::printf("WDM6-CPP_PRE_G13F %3d %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qi_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supcol = Real(t0c) - t_arr(i,j,k);

                if (qi_arr(i,j,k) > Real(0.0)) {
                    const Real qimax = Real(m_roqimax) / den_arr(i,j,k);
                    psaut_arr(i,j,k) = amrex::max(
                        Real(0.0), (qi_arr(i,j,k) - qimax) / dtcld);
                }

                if (qs_arr(i,j,k) > Real(0.0)) {
                    // Fortran :2510 alpha2 = 1.e-3*exp(0.09*(-supcol)). Both
                    // literals are unsuffixed and they push OPPOSITE ways, so
                    // neither can be judged alone: 1.e-3 contributes +4.749745e-08
                    // and 0.09 contributes -3.5763e-09*supcol, and the residual
                    // is their signed sum. Verified with zero free parameters
                    // against measured supcol at (199,3), k=91..98: predicted
                    // and measured relative error agree to 7 digits at every
                    // cell, ratio 1.0000. At k=91 that value is 7.652614e-08,
                    // which is exactly the step-2 rhoQ6/qgraup plotfile
                    // residual, so this literal pair is that residual.
                    const Real alpha2 = wdm6_literal(1.e-3)
                        * std::exp(wdm6_literal(0.09) * (-supcol));
                    pgaut_arr(i,j,k) = amrex::min(
                        amrex::max(
                            Real(0.0), alpha2 * (qs_arr(i,j,k) - Real(qs0))),
                        qs_arr(i,j,k) / dtcld);
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg)
                std::printf("WDM6-CPP_POST_G13F %3d %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(psaut_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(pgaut_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
            }
#endif

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                const Real diag_supcol = Real(t0c) - t_arr(diag_i,diag_j,kdbg);
                std::printf("WDM6-CPP_PRE_G13G %3d %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(rhw_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(diag_supcol));
                std::fflush(stdout);
                }
            }
#endif

            ParallelFor(box, [=, this] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real supcol = Real(t0c) - t_arr(i,j,k);
                if (supcol < Real(0.0)) {
                    const Real n0sfac = amrex::max(
                        amrex::min(std::exp(Real(alpha_wdm6) * supcol),
                                   Real(n0smax) / Real(n0s)),
                        Real(1.0));

                    if (qs_arr(i,j,k) > Real(0.0)
                        && rhw_arr(i,j,k) < Real(1.0)) {
                        const Real coeres = rslope2_arr(i,j,k,1)
                                          * std::sqrt(rslope_arr(i,j,k,1)
                                                      * rslopeb_arr(i,j,k,1));
                        psevp_arr(i,j,k) = (rhw_arr(i,j,k) - Real(1.0))
                            * n0sfac
                            * (Real(m_precs1) * rslope2_arr(i,j,k,1)
                               + Real(m_precs2) * work2_arr(i,j,k) * coeres)
                            / work1_arr(i,j,k,0);
                        psevp_arr(i,j,k) = amrex::min(
                            amrex::max(psevp_arr(i,j,k),
                                       -qs_arr(i,j,k) / dtcld),
                            Real(0.0));
                    }

                    if (qg_arr(i,j,k) > Real(0.0)
                        && rhw_arr(i,j,k) < Real(1.0)) {
                        const Real coeres = rslope2_arr(i,j,k,2)
                                          * std::sqrt(rslope_arr(i,j,k,2)
                                                      * rslopeb_arr(i,j,k,2));
                        pgevp_arr(i,j,k) = (rhw_arr(i,j,k) - Real(1.0))
                            * (Real(m_precg1) * rslope2_arr(i,j,k,2)
                               + Real(m_precg2) * work2_arr(i,j,k) * coeres)
                            / work1_arr(i,j,k,0);
                        pgevp_arr(i,j,k) = amrex::min(
                            amrex::max(pgevp_arr(i,j,k),
                                       -qg_arr(i,j,k) / dtcld),
                            Real(0.0));
                    }
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G13G %3d %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(psevp_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(pgevp_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G14 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qv_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qi_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real qmin_l = Real(qmin);
                const Real qcrmin_l = Real(qcrmin);
                const Real ncmin_l = Real(ncmin);
                const Real nrmin_l = Real(nrmin);
                const Real t0c_l = Real(t0c);
                const Real one = Real(1.0);
                const Real zero = Real(0.0);
                const Real delta2 =
                    (qr_arr(i,j,k) < Real(1.0e-4) && qs_arr(i,j,k) < Real(1.0e-4))
                    ? one : zero;
                const Real delta3 =
                    (qr_arr(i,j,k) < Real(1.0e-4)) ? one : zero;

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
                             - pracs_arr(i,j,k) * (one - delta2)
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
                             + piacr_arr(i,j,k) * (one - delta3)
                             + praci_arr(i,j,k) * (one - delta3)
                             + psacr_arr(i,j,k) * (one - delta2)
                             + pracs_arr(i,j,k) * (one - delta2)
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

                    value = amrex::max(ncmin_l, nc_arr(i,j,k));
                    source = (nrauto_arr(i,j,k) + nccol_arr(i,j,k)
                            + nraccr_arr(i,j,k) + naacw_arr(i,j,k)
                            + naacw_arr(i,j,k)) * dtcld;
                    if (source > value) {
                        factor = value / source;
                        nrauto_arr(i,j,k) = nrauto_arr(i,j,k) * factor;
                        nccol_arr(i,j,k) = nccol_arr(i,j,k) * factor;
                        nraccr_arr(i,j,k) = nraccr_arr(i,j,k) * factor;
                        naacw_arr(i,j,k) = naacw_arr(i,j,k) * factor;
                    }

                    value = amrex::max(nrmin_l, nr_arr(i,j,k));
                    source = (-nrauto_arr(i,j,k) + nrcol_arr(i,j,k)
                            + niacr_arr(i,j,k) + nsacr_arr(i,j,k)
                            + ngacr_arr(i,j,k)) * dtcld;
                    if (source > value) {
                        factor = value / source;
                        nrauto_arr(i,j,k) = nrauto_arr(i,j,k) * factor;
                        nrcol_arr(i,j,k) = nrcol_arr(i,j,k) * factor;
                        niacr_arr(i,j,k) = niacr_arr(i,j,k) * factor;
                        nsacr_arr(i,j,k) = nsacr_arr(i,j,k) * factor;
                        ngacr_arr(i,j,k) = ngacr_arr(i,j,k) * factor;
                    }

                    work2_arr(i,j,k) = -(prevp_arr(i,j,k) + psdep_arr(i,j,k)
                                       + pgdep_arr(i,j,k) + pigen_arr(i,j,k)
                                       + pidep_arr(i,j,k));
                    qv_arr(i,j,k) = qv_arr(i,j,k) + work2_arr(i,j,k) * dtcld;
                    qc_arr(i,j,k) = amrex::max(
                        qc_arr(i,j,k) - (praut_arr(i,j,k) + pracw_arr(i,j,k)
                                       + paacw_arr(i,j,k) + paacw_arr(i,j,k))
                                       * dtcld,
                        zero);
                    qr_arr(i,j,k) = amrex::max(
                        qr_arr(i,j,k) + (praut_arr(i,j,k) + pracw_arr(i,j,k)
                                       + prevp_arr(i,j,k) - piacr_arr(i,j,k)
                                       - pgacr_arr(i,j,k) - psacr_arr(i,j,k))
                                       * dtcld,
                        zero);
                    qi_arr(i,j,k) = amrex::max(
                        qi_arr(i,j,k) - (psaut_arr(i,j,k) + praci_arr(i,j,k)
                                       + psaci_arr(i,j,k) + pgaci_arr(i,j,k)
                                       - pigen_arr(i,j,k) - pidep_arr(i,j,k))
                                       * dtcld,
                        zero);
                    qs_arr(i,j,k) = amrex::max(
                        qs_arr(i,j,k) + (psdep_arr(i,j,k) + psaut_arr(i,j,k)
                                       + paacw_arr(i,j,k) - pgaut_arr(i,j,k)
                                       + piacr_arr(i,j,k) * delta3
                                       + praci_arr(i,j,k) * delta3
                                       + psaci_arr(i,j,k) - pgacs_arr(i,j,k)
                                       - pracs_arr(i,j,k) * (one - delta2)
                                       + psacr_arr(i,j,k) * delta2) * dtcld,
                        zero);
                    qg_arr(i,j,k) = amrex::max(
                        qg_arr(i,j,k) + (pgdep_arr(i,j,k) + pgaut_arr(i,j,k)
                                       + piacr_arr(i,j,k) * (one - delta3)
                                       + praci_arr(i,j,k) * (one - delta3)
                                       + psacr_arr(i,j,k) * (one - delta2)
                                       + pracs_arr(i,j,k) * (one - delta2)
                                       + pgaci_arr(i,j,k) + paacw_arr(i,j,k)
                                       + pgacr_arr(i,j,k) + pgacs_arr(i,j,k))
                                       * dtcld,
                        zero);
                    nc_arr(i,j,k) = amrex::max(
                        nc_arr(i,j,k) + (-nrauto_arr(i,j,k) - nccol_arr(i,j,k)
                                       - nraccr_arr(i,j,k) - naacw_arr(i,j,k)
                                       - naacw_arr(i,j,k)) * dtcld,
                        zero);
                    nr_arr(i,j,k) = amrex::max(
                        nr_arr(i,j,k) + (nrauto_arr(i,j,k) - nrcol_arr(i,j,k)
                                       - niacr_arr(i,j,k) - nsacr_arr(i,j,k)
                                       - ngacr_arr(i,j,k)) * dtcld,
                        zero);
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

                    value = amrex::max(ncmin_l, nc_arr(i,j,k));
                    source = (nrauto_arr(i,j,k) + nccol_arr(i,j,k)
                            + nraccr_arr(i,j,k) + naacw_arr(i,j,k)
                            + naacw_arr(i,j,k)) * dtcld;
                    if (source > value) {
                        factor = value / source;
                        nrauto_arr(i,j,k) = nrauto_arr(i,j,k) * factor;
                        nccol_arr(i,j,k) = nccol_arr(i,j,k) * factor;
                        nraccr_arr(i,j,k) = nraccr_arr(i,j,k) * factor;
                        naacw_arr(i,j,k) = naacw_arr(i,j,k) * factor;
                    }

                    value = amrex::max(nrmin_l, nr_arr(i,j,k));
                    source = (-nrauto_arr(i,j,k) + nrcol_arr(i,j,k)
                            - nseml_arr(i,j,k) - ngeml_arr(i,j,k)) * dtcld;
                    if (source > value) {
                        factor = value / source;
                        nrauto_arr(i,j,k) = nrauto_arr(i,j,k) * factor;
                        nrcol_arr(i,j,k) = nrcol_arr(i,j,k) * factor;
                        nseml_arr(i,j,k) = nseml_arr(i,j,k) * factor;
                        ngeml_arr(i,j,k) = ngeml_arr(i,j,k) * factor;
                    }

                    work2_arr(i,j,k) = -(prevp_arr(i,j,k) + psevp_arr(i,j,k)
                                       + pgevp_arr(i,j,k));
                    qv_arr(i,j,k) = qv_arr(i,j,k) + work2_arr(i,j,k) * dtcld;
                    qc_arr(i,j,k) = amrex::max(
                        qc_arr(i,j,k) - (praut_arr(i,j,k) + pracw_arr(i,j,k)
                                       + paacw_arr(i,j,k) + paacw_arr(i,j,k))
                                       * dtcld,
                        zero);
                    qr_arr(i,j,k) = amrex::max(
                        qr_arr(i,j,k) + (praut_arr(i,j,k) + pracw_arr(i,j,k)
                                       + prevp_arr(i,j,k) + paacw_arr(i,j,k)
                                       + paacw_arr(i,j,k) - pseml_arr(i,j,k)
                                       - pgeml_arr(i,j,k)) * dtcld,
                        zero);
                    qs_arr(i,j,k) = amrex::max(
                        qs_arr(i,j,k) + (psevp_arr(i,j,k) - pgacs_arr(i,j,k)
                                       + pseml_arr(i,j,k)) * dtcld,
                        zero);
                    qg_arr(i,j,k) = amrex::max(
                        qg_arr(i,j,k) + (pgacs_arr(i,j,k) + pgevp_arr(i,j,k)
                                       + pgeml_arr(i,j,k)) * dtcld,
                        zero);
                    nc_arr(i,j,k) = amrex::max(
                        nc_arr(i,j,k) + (-nrauto_arr(i,j,k) - nccol_arr(i,j,k)
                                       - nraccr_arr(i,j,k) - naacw_arr(i,j,k)
                                       - naacw_arr(i,j,k)) * dtcld,
                        zero);
                    nr_arr(i,j,k) = amrex::max(
                        nr_arr(i,j,k) + (nrauto_arr(i,j,k) - nrcol_arr(i,j,k)
                                       + nseml_arr(i,j,k) + ngeml_arr(i,j,k))
                                       * dtcld,
                        zero);
                    xlf = Real(xls) - xl_arr(i,j,k);
                    xlwork2 = -xl_arr(i,j,k) * (prevp_arr(i,j,k)
                                              + psevp_arr(i,j,k)
                                              + pgevp_arr(i,j,k))
                            - xlf * (pseml_arr(i,j,k) + pgeml_arr(i,j,k));
                    t_arr(i,j,k) = t_arr(i,j,k) - xlwork2 / cpm_arr(i,j,k) * dtcld;
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G14 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qv_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qi_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qs_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qg_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
            if (microphysics_debug >= 1 && diag_col_in_tile) {
                const int k_lo = amrex::max(klo, 39);
                const int k_hi = amrex::min(khi, 41);
                for (int kk = k_lo; kk <= k_hi; ++kk) {
                    const Real work2_cold = -(prevp_arr(diag_i,diag_j,kk)
                                            + psdep_arr(diag_i,diag_j,kk)
                                            + pgdep_arr(diag_i,diag_j,kk)
                                            + pigen_arr(diag_i,diag_j,kk)
                                            + pidep_arr(diag_i,diag_j,kk));
                    const Real xlf = Real(xls) - xl_arr(diag_i,diag_j,kk);
                    const Real xlwork2_cold = -Real(xls) * (psdep_arr(diag_i,diag_j,kk)
                                                          + pgdep_arr(diag_i,diag_j,kk)
                                                          + pidep_arr(diag_i,diag_j,kk)
                                                          + pigen_arr(diag_i,diag_j,kk))
                                             - xl_arr(diag_i,diag_j,kk) * prevp_arr(diag_i,diag_j,kk)
                                             - xlf * (piacr_arr(diag_i,diag_j,kk)
                                                    + paacw_arr(diag_i,diag_j,kk)
                                                    + paacw_arr(diag_i,diag_j,kk)
                                                    + pgacr_arr(diag_i,diag_j,kk)
                                                    + psacr_arr(diag_i,diag_j,kk));
                    std::printf("WDM6-CPP_T2_G14_POST k=%3d qv=%24.16E qc=%24.16E nn=%24.16E nc=%24.16E t=%24.16E\n",
                                kk + 1,
                                static_cast<double>(qv_arr(diag_i,diag_j,kk)),
                                static_cast<double>(qc_arr(diag_i,diag_j,kk)),
                                static_cast<double>(nn_arr(diag_i,diag_j,kk)),
                                static_cast<double>(nc_arr(diag_i,diag_j,kk)),
                                static_cast<double>(t_arr(diag_i,diag_j,kk)));
                    std::printf("WDM6-CPP_T2_G14_COLD_VAPOR k=%3d work2=%24.16E prevp=%24.16E psdep=%24.16E pgdep=%24.16E pigen=%24.16E pidep=%24.16E dtcld=%24.16E\n",
                                kk + 1,
                                static_cast<double>(work2_cold),
                                static_cast<double>(prevp_arr(diag_i,diag_j,kk)),
                                static_cast<double>(psdep_arr(diag_i,diag_j,kk)),
                                static_cast<double>(pgdep_arr(diag_i,diag_j,kk)),
                                static_cast<double>(pigen_arr(diag_i,diag_j,kk)),
                                static_cast<double>(pidep_arr(diag_i,diag_j,kk)),
                                static_cast<double>(dtcld));
                    std::printf("WDM6-CPP_T2_G14_COLD_THERMO k=%3d xlwork2=%24.16E prevp=%24.16E psdep=%24.16E pgdep=%24.16E pigen=%24.16E pidep=%24.16E piacr=%24.16E paacw=%24.16E pgacr=%24.16E\n",
                                kk + 1,
                                static_cast<double>(xlwork2_cold),
                                static_cast<double>(prevp_arr(diag_i,diag_j,kk)),
                                static_cast<double>(psdep_arr(diag_i,diag_j,kk)),
                                static_cast<double>(pgdep_arr(diag_i,diag_j,kk)),
                                static_cast<double>(pigen_arr(diag_i,diag_j,kk)),
                                static_cast<double>(pidep_arr(diag_i,diag_j,kk)),
                                static_cast<double>(piacr_arr(diag_i,diag_j,kk)),
                                static_cast<double>(paacw_arr(diag_i,diag_j,kk)),
                                static_cast<double>(pgacr_arr(diag_i,diag_j,kk)));
                }
                std::fflush(stdout);
            }
#endif

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G15 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qv_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(p_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qsatw_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qsati_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(rhw_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

            const Real hsub = Real(xls);
            const Real hvap = Real(xlv0);
            const Real cvap = Real(cpv);
            const Real ttp = Real(t0c) + wdm6_literal(0.01);
            const Real dldt = cvap - Real(cliq);
            const Real xa = -dldt / Real(rv);
            const Real xb = xa + hvap / (Real(rv) * ttp);
            const Real dldti = cvap - Real(cice);
            const Real xai = -dldti / Real(rv);
            const Real xbi = xai + hsub / (Real(rv) * ttp);

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real tr = ttp / t_arr(i,j,k);
                Real qsw = Real(psat) * std::exp(std::log(tr) * xa)
                                     * std::exp(xb * (Real(1.) - tr));
                qsw = amrex::min(qsw, Real(0.99) * p_arr(i,j,k));
                qsw = Real(ep2) * qsw / (p_arr(i,j,k) - qsw);
                qsw = amrex::max(qsw, Real(qmin));
                qsatw_arr(i,j,k) = qsw;

                Real qsi;
                if (t_arr(i,j,k) < ttp) {
                    qsi = Real(psat) * std::exp(std::log(tr) * xai)
                        * std::exp(xbi * (Real(1.) - tr));
                } else {
                    qsi = Real(psat) * std::exp(std::log(tr) * xa)
                        * std::exp(xb * (Real(1.) - tr));
                }
                qsi = amrex::min(qsi, Real(0.99) * p_arr(i,j,k));
                qsi = Real(ep2) * qsi / (p_arr(i,j,k) - qsi);
                qsi = amrex::max(qsi, Real(qmin));
                qsati_arr(i,j,k) = qsi;

                rhw_arr(i,j,k) = amrex::max(qv_arr(i,j,k) / qsw, Real(qmin));
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G15 %3d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qv_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(p_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(t_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qsatw_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qsati_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(rhw_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
            if (microphysics_debug >= 1 && diag_col_in_tile) {
                const int k_lo = amrex::max(klo, 39);
                const int k_hi = amrex::min(khi, 41);
                for (int kk = k_lo; kk <= k_hi; ++kk) {
                    std::printf("WDM6-CPP_T2_G15_POST k=%3d qv=%24.16E qsatw=%24.16E rhw=%24.16E t=%24.16E\n",
                                kk + 1,
                                static_cast<double>(qv_arr(diag_i,diag_j,kk)),
                                static_cast<double>(qsatw_arr(diag_i,diag_j,kk)),
                                static_cast<double>(rhw_arr(diag_i,diag_j,kk)),
                                static_cast<double>(t_arr(diag_i,diag_j,kk)));
                }
                std::fflush(stdout);
            }
#endif

            const Real g16a_cbrt24 = Real(std::pow(24.0f, 0.3333333f));

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G16A %3d %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                qrs_tmp_arr(i,j,k,0) = qr_arr(i,j,k);
                qrs_tmp_arr(i,j,k,1) = qs_arr(i,j,k);
                qrs_tmp_arr(i,j,k,2) = qg_arr(i,j,k);
                ncr_tmp_arr(i,j,k)   = nr_arr(i,j,k);
            });

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real rain_rslope, rain_rslopeb, rain_rslope2, rain_rslope3, rain_vt, rain_vtn;
                Real snow_rslope, snow_rslopeb, snow_rslope2, snow_rslope3, snow_vt, snow_n0sfac;
                Real graup_rslope, graup_rslopeb, graup_rslope2, graup_rslope3, graup_vt;

                wdm6_slope_rain_cell(qrs_tmp_arr(i,j,k,0), ncr_tmp_arr(i,j,k),
                                     den_arr(i,j,k), denfac_arr(i,j,k),
                                     Real(qcrmin), Real(nrmin),
                                     rslopermax_loc, rsloperbmax_loc,
                                     rsloper2max_loc, rsloper3max_loc,
                                     Real(bvtr), pvtr_loc, pvtrn_loc, pidnr_loc,
                                     rain_rslope, rain_rslopeb, rain_rslope2, rain_rslope3,
                                     rain_vt, rain_vtn);
                wdm6_slope_snow_cell(qrs_tmp_arr(i,j,k,1), den_arr(i,j,k), denfac_arr(i,j,k),
                                     t_arr(i,j,k), pidn0s_loc, Real(alpha_wdm6),
                                     Real(n0smax), Real(n0s), wdm6_slope_t0c, Real(qcrmin),
                                     rslopesmax_loc, rslopesbmax_loc,
                                     rslopes2max_loc, rslopes3max_loc,
                                     Real(bvts), pvts_loc,
                                     snow_rslope, snow_rslopeb, snow_rslope2, snow_rslope3,
                                     snow_vt, snow_n0sfac);
                wdm6_slope_graup_cell(qrs_tmp_arr(i,j,k,2), den_arr(i,j,k), denfac_arr(i,j,k),
                                      pidn0g_loc, Real(qcrmin),
                                      rslopegmax_loc, rslopegbmax_loc,
                                      rslopeg2max_loc, rslopeg3max_loc,
                                      slope_bvtg_loc, pvtg_loc,
                                      graup_rslope, graup_rslopeb, graup_rslope2, graup_rslope3,
                                      graup_vt);

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

                avedia_arr(i,j,k,1) = rslope_arr(i,j,k,0) * g16a_cbrt24;
                if (avedia_arr(i,j,k,1) <= Real(di82)) {
                    nc_arr(i,j,k) += nr_arr(i,j,k);
                    nr_arr(i,j,k) = Real(0.0);
                    qc_arr(i,j,k) += qr_arr(i,j,k);
                    qr_arr(i,j,k) = Real(0.0);
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G16A %3d %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(avedia_arr(diag_i,diag_j,kdbg,1)));
                std::fflush(stdout);
                }
            }
#endif

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G16B %3d %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(rhw_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nn_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
            if (microphysics_debug >= 1 && diag_col_in_tile) {
                const int k_lo = amrex::max(klo, 39);
                const int k_hi = amrex::min(khi, 41);
                for (int kk = k_lo; kk <= k_hi; ++kk) {
                    std::printf("WDM6-CPP_T2_G16B_ACT_PRE k=%3d RH=%12.4E qc=%12.4E nn,nc=%12.4E %12.4E\n",
                                kk + 1,
                                static_cast<double>(rhw_arr(diag_i,diag_j,kk)),
                                static_cast<double>(qc_arr(diag_i,diag_j,kk)),
                                static_cast<double>(nn_arr(diag_i,diag_j,kk)),
                                static_cast<double>(nc_arr(diag_i,diag_j,kk)));
                    std::printf("WDM6-CPP_T2_G16B_ACT_PRE_EXACT k=%3d rh=%24.16E qc=%24.16E nn=%24.16E nc=%24.16E\n",
                                kk + 1,
                                static_cast<double>(rhw_arr(diag_i,diag_j,kk)),
                                static_cast<double>(qc_arr(diag_i,diag_j,kk)),
                                static_cast<double>(nn_arr(diag_i,diag_j,kk)),
                                static_cast<double>(nc_arr(diag_i,diag_j,kk)));
                }
                std::fflush(stdout);
            }
#endif

            const Real g16b_pi = Real(4.0f * std::atan(1.0f));

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (rhw_arr(i,j,k) > Real(1.0)) {
                    const Real ratio = rhw_arr(i,j,k) / Real(satmax);
                    const Real fraction = amrex::min(Real(1.0),
                        std::exp(std::log(ratio) * Real(actk)));
                    Real ncact_raw = (nn_arr(i,j,k) + nc_arr(i,j,k)) * fraction - nc_arr(i,j,k);
                    Real ncact = amrex::max(Real(0.0), ncact_raw);
                    ncact /= dtcld;
                    const Real ncact_cap = amrex::max(nn_arr(i,j,k), Real(0.0)) / dtcld;
                    ncact = amrex::min(ncact, ncact_cap);
                    const Real actr_um = Real(actr) * wdm6_literal(1.0e-6);
                    const Real pcact = amrex::min(
                        Real(4.0) * g16b_pi * Real(denr)
                            * actr_um * actr_um * actr_um * ncact
                            / (Real(3.0) * den_arr(i,j,k)),
                        amrex::max(qv_arr(i,j,k), Real(0.0)) / dtcld);

                    ncact_arr(i,j,k) = ncact;
                    act_ratio_arr(i,j,k) = ratio;
                    act_fraction_arr(i,j,k) = fraction;
                    act_raw_arr(i,j,k) = ncact_raw;
                    act_cap_arr(i,j,k) = ncact_cap;
                    pcact_arr(i,j,k) = pcact;
                    qv_arr(i,j,k) = amrex::max(qv_arr(i,j,k) - pcact * dtcld, Real(0.0));
                    qc_arr(i,j,k) = amrex::max(qc_arr(i,j,k) + pcact * dtcld, Real(0.0));
                    nn_arr(i,j,k) = amrex::max(nn_arr(i,j,k) - ncact * dtcld, Real(0.0));
                    nc_arr(i,j,k) = amrex::max(nc_arr(i,j,k) + ncact * dtcld, Real(0.0));
                    t_arr(i,j,k) += pcact * xl_arr(i,j,k) / cpm_arr(i,j,k) * dtcld;
                }

                const Real tr = ttp / t_arr(i,j,k);
                Real qsw = Real(psat) * std::exp(std::log(tr) * xa)
                                     * std::exp(xb * (Real(1.0) - tr));
                qsw = amrex::min(qsw, Real(0.99) * p_arr(i,j,k));
                qsw = Real(ep2) * qsw / (p_arr(i,j,k) - qsw);
                qsw = amrex::max(qsw, Real(qmin));
                qsatw_arr(i,j,k) = qsw;

                work1_arr(i,j,k,0) = wdm6_conden(
                    t_arr(i,j,k), qv_arr(i,j,k), qsw, xl_arr(i,j,k), cpm_arr(i,j,k),
                    Real(qmin), Real(rv));
                work2_arr(i,j,k) = qc_arr(i,j,k) + work1_arr(i,j,k,0);

                Real pcond = amrex::min(
                    amrex::max(work1_arr(i,j,k,0) / dtcld, Real(0.0)),
                    amrex::max(qv_arr(i,j,k), Real(0.0)) / dtcld);
                if (qc_arr(i,j,k) > Real(0.0) && work1_arr(i,j,k,0) < Real(0.0)) {
                    pcond = amrex::max(work1_arr(i,j,k,0), -qc_arr(i,j,k)) / dtcld;
                }
                pcond_arr(i,j,k) = pcond;

                if (pcond == -qc_arr(i,j,k) / dtcld) {
                    nn_arr(i,j,k) += nc_arr(i,j,k);
                    nc_arr(i,j,k) = Real(0.0);
                }

                qv_arr(i,j,k) = amrex::max(qv_arr(i,j,k) - pcond * dtcld, Real(0.0));
                qc_arr(i,j,k) = amrex::max(qc_arr(i,j,k) + pcond * dtcld, Real(0.0));
                t_arr(i,j,k) += pcond * xl_arr(i,j,k) / cpm_arr(i,j,k) * dtcld;
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G16B %3d %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qv_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nn_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(pcond_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
            if (microphysics_debug >= 1 && diag_col_in_tile) {
                const int k_lo = amrex::max(klo, 39);
                const int k_hi = amrex::min(khi, 41);
                for (int kk = k_lo; kk <= k_hi; ++kk) {
                    std::printf("WDM6-CPP_T2_G16B_ACT_EXACT k=%3d ratio=%24.16E frac=%24.16E raw=%24.16E cap=%24.16E\n",
                                kk + 1,
                                static_cast<double>(act_ratio_arr(diag_i,diag_j,kk)),
                                static_cast<double>(act_fraction_arr(diag_i,diag_j,kk)),
                                static_cast<double>(act_raw_arr(diag_i,diag_j,kk)),
                                static_cast<double>(act_cap_arr(diag_i,diag_j,kk)));
                    std::printf("WDM6-CPP_T2_G16B_ACT_POST k=%3d ncact,pcact=%12.4E %12.4E nn_after,nc_after=%12.4E %12.4E\n",
                                kk + 1,
                                static_cast<double>(ncact_arr(diag_i,diag_j,kk)),
                                static_cast<double>(pcact_arr(diag_i,diag_j,kk)),
                                static_cast<double>(nn_arr(diag_i,diag_j,kk)),
                                static_cast<double>(nc_arr(diag_i,diag_j,kk)));
                    if (std::abs(static_cast<double>(pcond_arr(diag_i,diag_j,kk))) > 1.e-20) {
                        std::printf("WDM6-CPP_T2_G16B_COND_POST k=%3d pcond=%12.4E qc=%12.4E\n",
                                    kk + 1,
                                    static_cast<double>(pcond_arr(diag_i,diag_j,kk)),
                                    static_cast<double>(qc_arr(diag_i,diag_j,kk)));
                    }
                }
                std::fflush(stdout);
            }
#endif

            const Real g17_pi = Real(4.0f * std::atan(1.0f));
            const Real g17_pidnc = g17_pi * Real(denr) / Real(6.0);
            const Real g17_pidnr = Real(4.0) * g17_pi * Real(denr);

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_PRE_G17 %3d %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qi_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (qc_arr(i,j,k) <= Real(qmin)) qc_arr(i,j,k) = Real(0.0);
                if (qi_arr(i,j,k) <= Real(qmin)) qi_arr(i,j,k) = Real(0.0);

                if (qr_arr(i,j,k) >= Real(qcrmin) && nr_arr(i,j,k) >= Real(nrmin)) {
                    Real lamdr = std::exp(std::log(
                        (g17_pidnr * nr_arr(i,j,k)) / (den_arr(i,j,k) * qr_arr(i,j,k))
                    ) * wdm6_literal(0.33333333));
                    if (lamdr <= Real(lamdarmin)) {
                        lamdr = Real(lamdarmin);
                        nr_arr(i,j,k) = den_arr(i,j,k) * qr_arr(i,j,k)
                                      * std::pow(lamdr, Real(3.0)) / g17_pidnr;
                    } else if (lamdr >= Real(lamdarmax)) {
                        lamdr = Real(lamdarmax);
                        nr_arr(i,j,k) = den_arr(i,j,k) * qr_arr(i,j,k)
                                      * std::pow(lamdr, Real(3.0)) / g17_pidnr;
                    }
                }

                if (qc_arr(i,j,k) >= Real(qmin) && nc_arr(i,j,k) >= Real(ncmin)) {
                    Real lamdc = std::exp(std::log(
                        (g17_pidnc * nc_arr(i,j,k)) / (den_arr(i,j,k) * qc_arr(i,j,k))
                    ) * wdm6_literal(0.33333333));
                    if (lamdc <= Real(lamdacmin)) {
                        lamdc = Real(lamdacmin);
                        nc_arr(i,j,k) = den_arr(i,j,k) * qc_arr(i,j,k)
                                      * std::pow(lamdc, Real(3.0)) / g17_pidnc;
                    } else if (lamdc >= Real(lamdacmax)) {
                        lamdc = Real(lamdacmax);
                        nc_arr(i,j,k) = den_arr(i,j,k) * qc_arr(i,j,k)
                                      * std::pow(lamdc, Real(3.0)) / g17_pidnc;
                    }
                }
            });

#if !defined(AMREX_USE_GPU)
            if (microphysics_debug > 0 && loop == 0 && diag_col_in_tile) {
                for (int kdbg = klo; kdbg <= khi; ++kdbg) {
                std::printf("WDM6-CPP_POST_G17 %3d %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                            kdbg + 1,
                            static_cast<double>(qc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qi_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(qr_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nc_arr(diag_i,diag_j,kdbg)),
                            static_cast<double>(nr_arr(diag_i,diag_j,kdbg)));
                std::fflush(stdout);
                }
            }
#endif

            // ============================================================
            // Step 9: Ice physics (simplified) — remaining processes
            // ============================================================
            // These processes handle ice (qi), snow (qs), and graupel (qg)
            // following simplified versions of WSM6 processes

#if 0
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

                const Real qc_value = amrex::max(qc_arr(i,j,k), Real(qmin));
                const Real paacw_source = (paacw_arr(i,j,k) + paacw_arr(i,j,k)) * dtcld;
                if (paacw_source > qc_value) {
                    const Real factor = qc_value / paacw_source;
                    paacw_arr(i,j,k) *= factor;
                }

                const Real nc_value = amrex::max(nc_arr(i,j,k), Real(ncmin));
                const Real naacw_source = (naacw_arr(i,j,k) + naacw_arr(i,j,k)) * dtcld;
                if (naacw_source > nc_value) {
                    const Real factor = nc_value / naacw_source;
                    naacw_arr(i,j,k) *= factor;
                }

                if (temp <= t0c_loc) {
                    const Real qcol = paacw_arr(i,j,k) * dtcld;
                    const Real ncol = naacw_arr(i,j,k) * dtcld;
                    const Real xlf = Real(xls) - xl_arr(i,j,k);

                    qc_arr(i,j,k) = amrex::max(qc_arr(i,j,k) - Real(2.0) * qcol, Real(0.0));
                    qs_arr(i,j,k) += qcol;
                    qg_arr(i,j,k) += qcol;
                    nc_arr(i,j,k) = amrex::max(nc_arr(i,j,k) - Real(2.0) * ncol, Real(0.0));
                    t_arr(i,j,k) += Real(2.0) * xlf * qcol / cpm_arr(i,j,k);
                } else {
                    const Real qcol = paacw_arr(i,j,k) * dtcld;
                    const Real ncol = naacw_arr(i,j,k) * dtcld;

                    qc_arr(i,j,k) = amrex::max(qc_arr(i,j,k) - Real(2.0) * qcol, Real(0.0));
                    qr_arr(i,j,k) += Real(2.0) * qcol;
                    nc_arr(i,j,k) = amrex::max(nc_arr(i,j,k) - Real(2.0) * ncol, Real(0.0));
                }
            });
#endif

            // ============================================================
            // DISABLED: Step 10 (enforce bounds)
            // ============================================================
            // Pre-group legacy clamps. The nc/nr floors (1.e1, 1.e-2) have no
            // Fortran counterpart and inject a divergence: they pin the native
            // nc/nr minima to exactly the floor values while the bridge leaves
            // them free. Non-negativity for the mixing ratios is already handled
            // by G10e (phase cleanup) and G17 (lambda-bound ncr recompute).
            // Re-enable only if a bounded group is shown to require it.
#if 0
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
#endif

        } // End minor timestep loop

        // ============================================================
        // DISABLED: Step 11 (simplified sedimentation for all species)
        // ============================================================
        // Pre-group legacy sedimentation, redundant with the bounded groups that
        // superseded it: G5b/G5c (nislfv rain/snow/graupel) and G8 (ice fall).
        // Those run inside the minor timestep loop via
        // wdm6_nislfv_rain_plm6_column; this block ran once more afterwards on
        // qs/qg/qi with hardcoded fall-speed constants, double-sedimenting them.
        // It sits after POST_G17, outside every group tag, so no tag observed it.
        // Ablation at Milestone A: rhoQ3 2.608e-3 -> 4.135e-7,
        // rhoQ5 6.374e-4 -> 3.920e-9, thermodynamics bit-identical.
        // Its rain branch was already hollow (precip_rain never accumulated).
#if 0
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
#endif

        // Bridge-wrapper site W1_THETAWB, native leg. Same emitter, same schema,
        // same kts..kte column, at the equivalent point: end of the kernel,
        // immediately before Copy_Micro_to_State forms RhoTheta = rho * theta
        // from mic_fab_vars[theta].
        emit_w1("PRE");

        // Convert updated temperature back to potential temperature, mirroring
        // the bridge leg above. Without this the native path left
        // mic_fab_vars[theta] at its pre-microphysics value while the bridge
        // leg advanced it, and Copy_Micro_to_State then formed
        // RhoTheta = rho * theta from the stale field, dropping microphysical
        // latent heating from the thermodynamic state. Confirmed by the
        // W1_THETAWB retreat: native PRE theta equalled POST theta at 100/100
        // levels, and the resulting divergence 6.625510053 at k=90 matched the
        // Milestone A theta error exactly.
        {
            constexpr Real p0_nat = 1.e5;           // Reference pressure (Pa)
            constexpr Real rdOcp_nat = R_d / Cp_d;  // R/cp = 0.286
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real exner = std::pow(p_arr(i,j,k) / p0_nat, rdOcp_nat);
                w1_theta(i,j,k) = t_arr(i,j,k) / exner;
            });
        }

        emit_w1("POST");

#ifdef ERF_USE_WDM6_FORT
        }
#endif

    } // MFIter loop
}

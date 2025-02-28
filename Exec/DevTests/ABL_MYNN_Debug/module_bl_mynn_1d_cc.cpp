#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <limits>

extern "C" void mynn_tendencies_cc(const int& kts,const int& kte, const Real & delt,
                                   /*in*/ const Real* dz,
                                   /*in*/ const Real* rho,
                                   /*in*/ const Real* u, const Real* v, const Real* th, const Real* tk, const Real* qv,
                                   /*in*/ const Real* qc, const Real* qi, const Real* qs, const Real* qni, const Real* qnc,
                                   /*in*/ const Real* psfc,const Real* p, const Real* exner,
                                   /*inout*/ Real* thl, Real* sqv, Real* sqc, Real* sqi,
                                   /*inout*/ Real* sqs, Real* sqw, Real* qnwfa, Real* qnifa, Real* qnbca, Real* ozone,
                                   /*in*/ Real* ust, const Real & flt,const Real & flq,const Real & flqv,const Real & flqc,const Real & wspd, const Real & uoce, const Real & voce,
                                   /*in*/ const Real* tcd,const Real* qcd,
                                   /*inout*/ Real* dfm, Real* dfh,
                                   /*inout*/ Real* du, Real* dv, Real* dth,
                                   /*inout*/ Real* dqv, Real* dqc, Real* dqi,
                                   /*inout*/ Real* dqs, Real* dqnc, Real* dqni,
                                   /*inout*/ Real* dqnwfa, Real* dqnifa, Real* dqnbca,
                                   /*inout*/ Real* dozone,
                                   /*in*/ const Real* diss_heat,
                        /*in*/ const Real* s_aw, const Real* s_awthl, const Real* s_awqt, const Real* s_awqv, const Real* s_awqc, const Real* s_awu, const Real* s_awv, const Real* s_awqni, const Real* s_awqnc, const Real* s_awqnwfa, const Real* s_awqnifa, const Real* s_awqnbca, const Real* sd_aw, const Real* sd_awthl, const Real* sd_awqt, const Real* sd_awqv, const Real* sd_awqc, const Real* sd_awu, const Real* sd_awv,
                                   const Real* sub_thl,const Real* sub_sqv,const Real* sub_u,const Real* sub_v,
                                   const Real* det_thl,const Real* det_sqv,const Real* det_sqc,const Real* det_u,const Real* det_v,
                                   /*logical turned into int */const int& flag_qc, const int& flag_qi, const int& flag_qnc, const int& flag_qni, const int& flag_qs, const int& flag_qnwfa, const int& flag_qnifa, const int& flag_qnbca, const int& flag_ozone,
                                   const int & bl_mynn_cloudmix, const int & bl_mynn_mixqt, int & bl_mynn_edmf_mom,  const int & bl_mynn_mixscalars, /* new */const int& debug_code,const Real& r_d,const Real& p608,const Real& ep_2,const Real& ep_3,const Real& tv0,const Real& xlv,const Real& xlvcp,const Real & xlscp);
/*
skipping bl_mynn_mixscalars, flag_qni, flag_qc, flag_qs, flag_qnc, flag_qnwfa, flag_qnifa, flag_qnbca, flag_ozone // adding
skipping th, qc, qi, qs, qni, qnc //adding
skipping exner, dfq, tsq, qsq, cov, cldfra_bl1d
skipping sqi, sqs, qnwfa, qnifa, qnbca, ozone,dqv, dqc, dqi, dqs, dqni, dqnc, dqnwfa, dqnifa, dqnbca, dozone //adding some
skipping wsp, wsp2, tk2, th2
skipping problem, kproblem
looks like constants from common that are passed in: Real r_d, Real p608, Real ep_2,Real ep_3,Real tv0,Real xlv,Real xlvcp
should these also be intent(in) and therefore const &?
 */
extern "C" void mym_predict_cc(int& kts, int& kte, Real& closure, Real& delt, Real* dz, Real& ust, Real& flt, Real& flq, Real& pmz, Real& phh, Real* el, Real* dfq, Real* rho, Real* pdk, Real* pdt, Real* pdq, Real* pdc, Real* qke, Real* tsq, Real* qsq, Real* cov, Real* s_aw, Real* s_awqke, int& bl_mynn_edmf_tke, Real* qwt1d, Real* qdiss1d, int& tke_budget, Real& xlvcp, Real& xlscp, Real& karman);

extern "C" void mynn_mix_chem_cc(int kts, int kte, int i,Real delt, Real* dz, Real pblh, int nchem, int kdvel, int ndvel,Real** chem1, Real* vd1, Real* rho,Real flt, Real* tcd, Real* qcd, Real* dfh,Real* s_aw, Real** s_awchem, Real emis_ant_no, Real frp, int rrfs_sd, int enh_mix);

extern "C" void moisture_check_cc(int kte, Real delt, Real* dp, Real* exner,Real* qv, Real* qc, Real* qi, Real* qs, Real* th,Real* dqv, Real* dqc, Real* dqi, Real* dqs, Real* dth,Real dqv2, Real xlvcp, Real xlscp);

extern "C" void mym_condensation_cc(const int& kts,const int& kte, const Real& dx, Real* dz, Real* zw, Real& xland,Real* thl,
                Real* qw, Real* qv, Real* qc, Real* qi, Real* qs,Real* p, Real* exner,
                Real* tsq, Real* qsq, Real* cov,Real* sh, Real* el, int& bl_mynn_cloudpdf,
                Real* qc_bl1d, Real* qi_bl1d, Real* cldfra_bl1d,Real & pblh1, Real & hfx1,
                         Real* vt, Real* vq, Real* th, Real* sgm, Real* rmo,int &spp_pbl, Real* rstoch_col,
                Real ep_2, Real ep_3, Real xlv, Real r_d, Real xlvcp, Real p608, Real tv0, Real cpv,
                                    Real r_v, Real cice, Real cliq, Real cp, Real xls, Real rcp);


extern "C" void topdown_cloudrad_cc(int& kts, int& kte, const Real* dz1, const Real* zw, Real& fltv, Real& xland, int& kpbl, Real& pblh, const Real* sqc, const Real* sqi, const Real* sqw, const Real* thl, const Real* th1, const Real* ex1, const Real* p1, const Real*  rho1, const Real* thetav, const Real* cldfra_bl1d, const Real* rthraten, Real& maxkhtopdown, Real* khtopdown, Real* tkeprodtd);

extern "C" void ddmf_jpl_cc(int& kts, int& kte, Real& dt, const Real* zw, const Real* dz, const Real* p,
              const Real* u, const Real* v, const Real* th, const Real* thl, const Real* thv,
              const Real* tk,const Real* qt, const Real* qv, const Real* qc, const Real*
              rho, const Real* exner,Real& ust, Real& wthl, Real& wqt, Real& pblh, int& kpbl,
              Real* edmf_a_dd, Real* edmf_w_dd, Real* edmf_qt_dd,
              Real* edmf_thl_dd, Real* edmf_ent_dd, Real* edmf_qc_dd,
              Real* sd_aw, Real* sd_awthl, Real* sd_awqt,
              Real* sd_awqv, Real* sd_awqc, Real* sd_awu,
              Real* sd_awv, Real* sd_awqke,
              const Real* qc_bl1d, const Real* cldfra_bl1d,
              const Real* rthraten, Real& svp1, Real& grav, Real& onethird, Real& p1000mb,
              Real& rcp, Real& xlvcp, Real& cp, Real& rvovrd );

extern "C" void scale_aware_cc(Real& dx, Real& pbl1, Real& psig_bl, Real& psig_shcu);

extern "C" void get_pblh_cc(int &kts, int &kte, Real &zi, Real *thetav1d, Real *qke1d, Real *zw1d, Real *dz1d, Real &landsea, int &kzi);

extern "C" void retrieve_exchange_coeffs_cc(int& kts, int& kte, Real* dfm, Real* dfh, const Real* dz, Real* k_m, Real* k_h);

extern "C" void dmp_mf_cc(const int& kts, const int& kte, Real& dt, Real* zw, Real* dz, Real* p, Real* rho, int& momentum_opt, int& tke_opt, int& scalar_opt, Real* u, Real* v, Real* w, Real* th, Real* thl, Real* thv, Real* tk, Real* qt, Real* qv, Real* qc, Real* qke, Real* qnc, Real* qni, Real* qnwfa, Real* qnifa, Real* qnbca, Real& ust, Real& flt, Real& fltv, Real& flq, Real& flqv, Real& pblh, int& kpbl, Real& dx, Real& landsea, Real& ts, Real* edmf_a, Real* edmf_w, Real* edmf_qt, Real* edmf_thl, Real* edmf_ent, Real* edmf_qc, Real* s_aw, Real* s_awthl, Real* s_awqt, Real* s_awqv, Real* s_awqc, Real* s_awu, Real* s_awv, Real* s_awqke, Real* s_awqnc, Real* s_awqni, Real* s_awqnwfa, Real* s_awqnifa, Real* s_awqnbca, int& nchem, Real** chem1, Real** s_awchem, bool& mix_chem, Real* qc_bl1d, Real* cldfra_bl1d, Real* qc_bl1d_old, Real* cldfra_bl1d_old, Real& psig_shcu, Real& maxwidth, int& ktop, Real& maxmf, Real& ztop, Real* rstoch_col, Real grav, Real gtr, Real p608);

extern "C" void mym_turbulence_cc(int& kts, int& kte, Real& xland, Real& closure, Real* dz, Real& dx, Real* zw, Real* u, Real* v, Real* thl, Real* thetav, Real* ql, Real* qw, Real* qke, Real* tsq, Real* qsq, Real* cov, Real* vt, Real* vq, Real& rmo, Real& flt, Real& fltv, Real& flq, Real& zi, Real* theta, Real* sh, Real* sm, Real* el, Real* dfm, Real* dfh, Real* dfq, Real* tcd, Real* qcd, Real* pdk, Real* pdt, Real* pdq, Real* pdc, Real* qWT1D, Real* qSHEAR1D, Real* qBUOY1D, Real* qDISS1D, int& tke_budget, Real& Psig_bl, Real& Psig_shcu, Real* cldfra_bl1D, int& bl_mynn_mixlength, Real* edmf_w1, Real* edmf_a1, Real* TKEprodTD, int& spp_pbl, Real* rstoch_col, int& debug_code, Real& gtr, Real& tv0);

extern "C" void mym_initialize_cc(const int &kts,const int &kte,const Real &xland, Real *dz, Real &dx, Real *zw, Real *u, Real *v, Real *thl, Real *qw,const Real &zi, Real *theta, Real *thetav, Real *sh, Real *sm,const Real& ust, const Real &rmo, Real* el, Real *qke, Real* tsq, Real* qsa, Real* cov, const Real& Psig_bl, Real *cldfra_bl1D, int &bl_mynn_mixlength, Real *edmf_w1, Real *edmf_a1, int &INITIALIZE_QKE, int &spp_pbl, Real *rstoch_col,const Real& karman,const Real& tv0,const Real& gtr);
//----------------------------------------contstants-------------------------------------------

// constants
const Real no_threshold = 10.0;     // for anthropogenic sources
const Real frp_threshold = 10.0;    // increased the frp threshold to enhance mixing over big fires
const Real pblh_threshold = 100.0;

const Real t0c = 273.15; // assuming t0c is 273.15
const Real tice = 240.0; // assuming tice is 240 based on the comment

// assuming Real corresponds to Real precision
const Real cphm_st = 5.0, cphm_unst = 16.0,
                 cphh_st = 5.0, cphh_unst = 16.0;
//    1.1800000667572021       0.1370676159858704       0.6645210385322571       0.5641299486160278
//       0.2710000276565552       0.6599999666213989      19.7362747192382812       1.9125051498413086       0.8616270422935486       2.5500068664550781       8.3544006347656250,

// closure constants
constexpr Real pr = 0.7400000095367432,
                 g1 = 0.2349999994039536, // nn2009 = 0.235
                 b1 = 24.0,
                 b2 = 15.0, // ckmod     nn2009
                 c2 = 0.7289999723434448, // 0.729, //0.75,
                 c3 = 0.3400000035762787, // 0.340, //0.352,
                 c4 = 0.0,
                 c5 = 0.2,
                 a1 = 1.1800000667572021,
                 c1 = 0.1370676159858704,
                 a2 = 0.6645210385322571,
                 g2 = 0.5641299486160278;

  constexpr Real cc2 = 1.0_rt - c2,
                 cc3 = 1.0_rt - c3,
                 e1c = 3.0_rt * a2 * b2 * cc3,
                 e2c = 9.0_rt * a1 * a2 * cc2,
                 e3c = 9.0_rt * a2 * a2 * cc2 * (1.0_rt - c5),
                 e4c = 12.0_rt * a1 * a2 * cc2,
                 e5c = 6.0_rt * a1 * a1;

// constants for min tke in elt integration (qmin), max z/l in els (zmax),
// and factor for eddy viscosity for tke (kq = sqfac*km):
constexpr Real qmin = 0.0, zmax = 1.0, sqfac = 3.0;

constexpr Real gpw = 5.0_rt / 3.0, qcgmin = 1e-8, qkemin = 1e-3;
constexpr Real tliq = 269.0; // all hydrometeors are liquid when t > tliq

// constants for cloud pdf (mym_condensation)
constexpr Real rr2 = 0.7071068, rrp = 0.3989423;

// use canuto/kitamura mod (remove ric and negative tke) (1:yes, 0:no)
constexpr Real ckmod = 1.0;

// option to activate environmental subsidence in mass-flux scheme
constexpr bool env_subs = false;

//---------------------------------------------------------------------------------------------
Real vsc = 1.0e-5;
Real elt = 1.0e-5;

Real esat_blend_cc(Real t) {
    // constants for liquid
    const Real j0 = .611583699e03;
    const Real j1 = .444606896e02;
    const Real j2 = .143177157e01;
    const Real j3 = .264224321e-1;
    const Real j4 = .299291081e-3;
    const Real j5 = .203154182e-5;
    const Real j6 = .702620698e-8;
    const Real j7 = .379534310e-11;
    const Real j8 = -.321582393e-13;

    // constants for ice
    const Real k0 = .609868993e03;
    const Real k1 = .499320233e02;
    const Real k2 = .184672631e01;
    const Real k3 = .402737184e-1;
    const Real k4 = .565392987e-3;
    const Real k5 = .521693933e-5;
    const Real k6 = .307839583e-7;
    const Real k7 = .105785160e-9;
    const Real k8 = .161444444e-12;

    Real xc = std::max(-80.0_rt, t - t0c);
    Real esat_blend_cc;

    if (t >= (t0c - 6.0_rt)) {
        esat_blend_cc = j0 + xc * (j1 + xc * (j2 + xc * (j3 + xc * (j4 + xc * (j5 + xc * (j6 + xc * (j7 + xc * j8)))))));
    } else if (t <= tice) {
        esat_blend_cc = k0 + xc * (k1 + xc * (k2 + xc * (k3 + xc * (k4 + xc * (k5 + xc * (k6 + xc * (k7 + xc * k8)))))));
    } else {
        Real esl = j0 + xc * (j1 + xc * (j2 + xc * (j3 + xc * (j4 + xc * (j5 + xc * (j6 + xc * (j7 + xc * j8)))))));
        Real esi = k0 + xc * (k1 + xc * (k2 + xc * (k3 + xc * (k4 + xc * (k5 + xc * (k6 + xc * (k7 + xc * k8)))))));
        Real chi = ((t0c - 6.0_rt) - t) / ((t0c - 6.0_rt) - tice);
        esat_blend_cc = (1.0_rt - chi) * esl + chi * esi;
    }

    return esat_blend_cc;
}


Real qsat_blend_cc(Real t, Real p) {
    // constants for liquid
    const Real j0 = .611583699e03;
    const Real j1 = .444606896e02;
    const Real j2 = .143177157e01;
    const Real j3 = .264224321e-1;
    const Real j4 = .299291081e-3;
    const Real j5 = .203154182e-5;
    const Real j6 = .702620698e-8;
    const Real j7 = .379534310e-11;
    const Real j8 = -.321582393e-13;

    // constants for ice
    const Real k0 = .609868993e03;
    const Real k1 = .499320233e02;
    const Real k2 = .184672631e01;
    const Real k3 = .402737184e-1;
    const Real k4 = .565392987e-3;
    const Real k5 = .521693933e-5;
    const Real k6 = .307839583e-7;
    const Real k7 = .105785160e-9;
    const Real k8 = .161444444e-12;


    // temperature thresholds
    const Real t0c = 273.15; // assuming 0 for t0c (temperature in celsius)
    const Real tice = 240.00; // assuming -273.15_rt for tice (std::absolute zero, could be different)
    Real xc = std::max(-80.0_rt, t - t0c);
    Real qsat_blend_cc, esl, esi, rslf, rsif, chi;

    if (t >= (t0c - 6.0_rt)) {
        esl = j0 + xc * (j1 + xc * (j2 + xc * (j3 + xc * (j4 + xc * (j5 + xc * (j6 + xc * (j7 + xc * j8)))))));
        esl = std::min(esl, p * 0.15_rt);
        qsat_blend_cc = 0.622_rt * esl / std::max(p - esl, 1e-5_rt);
    } else if (t <= tice) {
        esi = k0 + xc * (k1 + xc * (k2 + xc * (k3 + xc * (k4 + xc * (k5 + xc * (k6 + xc * (k7 + xc * k8)))))));
        esi = std::min(esi, p * 0.15_rt);
        qsat_blend_cc = 0.622_rt * esi / std::max(p - esi, 1e-5_rt);
    } else {
        esl = j0 + xc * (j1 + xc * (j2 + xc * (j3 + xc * (j4 + xc * (j5 + xc * (j6 + xc * (j7 + xc * j8)))))));
        esl = std::min(esl, p * 0.15_rt);
        esi = k0 + xc * (k1 + xc * (k2 + xc * (k3 + xc * (k4 + xc * (k5 + xc * (k6 + xc * (k7 + xc * k8)))))));
        esi = std::min(esi, p * 0.15_rt);
        rslf = 0.622_rt * esl / std::max(p - esl, 1e-5_rt);
        rsif = 0.622_rt * esi / std::max(p - esi, 1e-5_rt);
        chi = ((t0c - 6.0_rt) - t) / ((t0c - 6.0_rt) - tice);
        qsat_blend_cc = (1.0_rt - chi) * rslf + chi * rsif;
    }
    return qsat_blend_cc;
}


Real xl_blend_cc(Real t,Real xlv, Real xls, Real cpv, Real cliq, Real cice) {
    Real xl_blend_cc, xlvt, xlst, chi;
    // t0c = 273.15, tice is set elsewhere
    if (t >= t0c) {
        xl_blend_cc = xlv + (cpv - cliq) * (t - t0c); // vaporization/condensation
    } else if (t <= tice) {
        xl_blend_cc = xls + (cpv - cice) * (t - t0c); // sublimation/deposition
    } else {
        xlvt = xlv + (cpv - cliq) * (t - t0c); // vaporization/condensation
        xlst = xls + (cpv - cice) * (t - t0c); // sublimation/deposition
        chi = (t0c - t) / (t0c - tice);
        xl_blend_cc = (1._rt - chi) * xlvt + chi * xlst; // blended
    }
    return xl_blend_cc;
}

void condensation_edmf_cc(Real qt, Real thl, Real p, Real zagl, Real& thv, Real& qc, Real p1000mb, Real rcp, Real xlvcp, Real rvovrd) {
    const int niter = 50;
    const Real diff = 1.e-6;
    Real exn = std::pow((p / p1000mb), rcp);
    // qc is assumed to be initialized before calling this function
    for (int i = 0; i < niter; ++i) {
        Real t = exn * thl + xlvcp * qc;
        Real qs = qsat_blend_cc(t, p);
        Real qcold = qc;
        qc = 0.5_rt * qc + 0.5_rt * std::max((qt - qs), 0.0_rt);
        if (std::abs(qc - qcold) < diff) break;
    }
    Real t = exn * thl + xlvcp * qc;
    Real qs = qsat_blend_cc(t, p);
    qc = std::max(qt - qs, 0.0_rt);
    // do not allow saturation below 1.0_rt m
    if (zagl < 100.0_rt) qc = 0.0;
    thv = (thl + xlvcp * qc) * (1.0_rt + qt * (rvovrd - 1.0_rt) - rvovrd * qc);
}

// function to solve system of linear equations on tridiagonal matrix n times n
// after peaceman and rachford, 1955
// a, b, c, d - are std::vectors of order n
// a, b, c - are coefficients on the lhs
// d - is initially rhs on the output becomes a solution std::vector
void tridiag_cc(int n, const Real* a, const Real* b, Real* c, Real* d) {
    Real q[n];
    c[n-1] = 0.0;
    q[0] = -c[0] / b[0];
    d[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        Real p = 1.0_rt / (b[i] + a[i] * q[i - 1]);
        q[i] = -c[i] * p;
        d[i] = (d[i] - a[i] * d[i - 1]) * p;
    }

    for (int i = n - 2; i >= 0; --i) {
        d[i] = d[i] + q[i] * d[i + 1];
    }
}

void tridiag2_cc(int n, const Real* a, const Real* b, const Real* c, const Real* d, Real* x) {
    Real cp[n];
    Real dp[n];
    Real m;

    // initialize c-prime and d-prime
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    // solve for std::vectors c-prime and d-prime
    for (int i = 1; i <= n; ++i) {
        m = b[i] - cp[i - 1] * a[i];
        cp[i] = c[i] / m;
        dp[i] = (d[i] - dp[i - 1] * a[i]) / m;
    }

    // initialize x
    x[n] = dp[n];

    // solve for x from the std::vectors c-prime and d-prime
    for (int i = n - 1; i >= 0; --i) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

}

// function to perform tridiagonal matrix algorithm
void tridiag3_cc(int kte, Real* a, Real* b, Real* c, Real* d, Real* x) {
    // inversion and resolution of a tridiagonal matrix a x = d
    // a - lower diagonal (ai,i-1)
    // b - principal diagonal (ai,i)
    // c - upper diagonal (ai,i+1)
    // d - right-hand side std::vector
    // x - solution std::vector

    for (int in = kte - 1; in >= 1; --in) {
        d[in] = d[in] - c[in] * d[in + 1] / b[in + 1];
        b[in] = b[in] - c[in] * a[in + 1] / b[in + 1];
    }
    for (int in = 1 + 1; in < kte; ++in) {
        d[in] = d[in] - a[in] * d[in - 1] / b[in - 1];
    }
    for (int in = 1; in < kte; ++in) {
        x[in] = d[in] / b[in];
    }
}


// ==================================================================
//>\ingroup gsd_mynn_edmf
// this subroutine was taken from the boulac scheme in wrf-arw
// and modified for integration into the mynn pbl scheme.
// while loops were added to reduce the computational expense.
// this subroutine computes the length scales up and down
// and then computes the min, average of the up/down
// length scales, and also considers the distance to the
// surface.
void boulac_length_cc(int kts, int kte,
                      const Real* zw, const Real* dz, const Real* qtke, const Real* theta,
                      Real* lb1, Real* lb2,
                      // model constant
                      Real gtr) {

//      dlu = the distance a parcel can be lifted upwards give a finite
//            amount of tke.
//      dld = the distance a parcel can be displaced downwards given a
//            finite amount of tke.
//      lb1 = the minimum of the length up and length down
//      lb2 = the average of the length up and length down
    int iz, izz, found;
    Real dlu[kte-kts];
    Real dld[kte-kts];
    const Real lmax = 2000.0;
    Real dzt, zup, beta, zup_inf, bbb, tl, zdo, zdo_sup, zzz;

    for (iz = kts; iz <= kte; iz++) {
        zup = 0.0;
        dlu[iz] = zw[kte + 1] - zw[iz] - dz[iz] * 0.5;
        zzz = 0.0;
        zup_inf = 0.0;
        beta = gtr;

        if (iz < kte) {
            found = 0;
            izz = iz;
            while (found == 0) {
                if (izz < kte) {
                    dzt = dz[izz];
                    zup = zup - beta * theta[iz] * dzt;
                    zup = zup + beta * (theta[izz + 1] + theta[izz]) * dzt * 0.5;
                    zzz = zzz + dzt;

                    if (qtke[iz] < zup && qtke[iz] >= zup_inf) {
                        bbb = (theta[izz + 1] - theta[izz]) / dzt;

                        if (bbb != 0.0_rt) {
                            tl = (-beta * (theta[izz] - theta[iz]) + std::sqrt(std::max(0.0_rt, ((beta * (theta[izz] - theta[iz])) * (beta * (theta[izz] - theta[iz]))) + 2.0_rt * bbb * beta * (qtke[iz] - zup_inf)))) / bbb / beta;
                        } else {
                            if (theta[izz] != theta[iz]) {
                                tl = (qtke[iz] - zup_inf) / (beta * (theta[izz] - theta[iz]));
                            } else {
                                tl = 0.0;
                            }
                        }

                        dlu[iz] = zzz - dzt + tl;
                        found = 1;
                    }

                    zup_inf = zup;
                    izz = izz + 1;
                } else {
                    found = 1;
                }
            }
        }

        zdo = 0.0;
        zdo_sup = 0.0;
        dld[iz] = zw[iz];
        zzz = 0.0;

        if (iz > kts) {
            found = 0;
            izz = iz;
            while (found == 0) {
                if (izz > kts) {
                    dzt = dz[izz - 1];
                    zdo = zdo + beta * theta[iz] * dzt;
                    zdo = zdo - beta * (theta[izz - 1] + theta[izz]) * dzt * 0.5;
                    zzz = zzz + dzt;

                    if (qtke[iz] < zdo && qtke[iz] >= zdo_sup) {
                        bbb = (theta[izz] - theta[izz - 1]) / dzt;

                        if (bbb != 0.0_rt) {
                            tl = (beta * (theta[izz] - theta[iz]) + std::sqrt(std::max(0.0_rt, ((beta * (theta[izz] - theta[iz])) * (beta * (theta[izz] - theta[iz]))) + 2.0_rt * bbb * beta * (qtke[iz] - zdo_sup)))) / bbb / beta;
                        } else {
                            if (theta[izz] != theta[iz]) {
                                tl = (qtke[iz] - zdo_sup) / (beta * (theta[izz] - theta[iz]));
                            } else {
                                tl = 0.0;
                            }
                        }

                        dld[iz] = zzz - dzt + tl;
                        found = 1;
                    }

                    zdo_sup = zdo;
                    izz = izz - 1;
                } else {
                    found = 1;
                }
            }
        }

        dld[iz] = std::min(dld[iz], zw[iz + 1]);
        lb1[iz] = std::min(dlu[iz], dld[iz]);
        dlu[iz] = std::max(0.1_rt, std::min(dlu[iz], 1000.0_rt));
        dld[iz] = std::max(0.1_rt, std::min(dld[iz], 1000.0_rt));
        lb2[iz] = std::sqrt(dlu[iz] * dld[iz]);
        lb1[iz] = lb1[iz] / (1.0_rt + (lb1[iz] / lmax));
        lb2[iz] = lb2[iz] / (1.0_rt + (lb2[iz] / lmax));

        if (iz == kte) {
            lb1[kte] = lb1[kte - 1];
            lb2[kte] = lb2[kte - 1];
        }
    }
}

//
// ==================================================================
//     subroutine  mym_level2:
//
//     input variables:    see subroutine mym_initialize
//
//     output variables:
//       dtl(nx,nz,ny) : vertical gradient of theta_l             (k/m)
//       dqw(nx,nz,ny) : vertical gradient of q_w
//       dtv(nx,nz,ny) : vertical gradient of theta_v             (k/m)
//       gm (nx,nz,ny) : g_m divided by l^2/q^2                (s^(-2))
//       gh (nx,nz,ny) : g_h divided by l^2/q^2                (s^(-2))
//       sm (nx,nz,ny) : stability function for momentum, at level 2
//       sh (nx,nz,ny) : stability function for heat, at level 2
//
//       these are defined on the walls of the grid boxes.
//

//>\ingroup gsd_mynn_edmf
// this subroutine calculates the level 2, non-dimensional wind shear
// \f$g_m\f$ and vertical temperature gradient \f$g_h\f$ as well as
// the level 2 stability functions \f$s_h\f$ and \f$s_m\f$.
//\param kts    horizontal dimension
//\param kte    vertical dimension
//\param dz     vertical grid spacings (\f$m\f$)
//\param u      west-east component of the horizontal wind (\f$m s^{-1}\f$)
//\param v      south-north component of the horizontal wind (\f$m s^{-1}\f$)
//\param thl    liquid water potential temperature
//\param qw     total water content \f$q_w\f$
//\param ql     liquid water content (\f$kg kg^{-1}\f$)
//\param vt
//\param vq
//\param dtl     vertical gradient of \f$\theta_l\f$ (\f$k m^{-1}\f$)
//\param dqw     vertical gradient of \f$q_w\f$
//\param dtv     vertical gradient of \f$\theta_v\f$ (\f$k m^{-1}\f$)
//\param gm      \f$g_m\f$ divided by \f$l^{2}/q^{2}\f$ (\f$s^{-2}\f$)
//\param gh      \f$g_h\f$ divided by \f$l^{2}/q^{2}\f$ (\f$s^{-2}\f$)
//\param sm      stability function for momentum, at level 2
//\param sh      stability function for heat, at level 2
//\section gen_mym_level2 gsd mynn-edmf mym_level2 general algorithm
// @ {
void mym_level2_cc(
    int kts, int kte,
    const Real* dz,
    const Real* u, const Real* v, const Real* thl, const Real* thetav, const Real* qw,
    const Real* ql, const Real* vt, const Real* vq,
    // intent(out):
    Real* dtl, Real* dqw, Real* dtv, Real* gm, Real* gh, Real* sm, Real* sh,
    // model constants:
    Real tv0, Real gtr)
{
    Real rfc, f1, f2, rf1, rf2, smc, shc, ri1, ri2, ri3, ri4, duz, dtz, dqz, vtt, vqq, dtq, dzk, afk, abk, ri, rf;
    Real a2fac;

    rfc = g1 / (g1 + g2);
    f1 = b1 * (g1 - c1) + 3.0_rt * a2 * (1.0_rt - c2) * (1.0_rt - c5) + 2.0_rt * a1 * (3.0_rt - 2.0_rt * c2);
    f2 = b1 * (g1 + g2) - 3.0_rt * a1 * (1.0_rt - c2);
    rf1 = b1 * (g1 - c1) / f1;
    rf2 = b1 * g1 / f2;
    smc = a1 / a2 * f1 / f2;
    shc = 3.0_rt * a2 * (g1 + g2);
    //    printf("g1 %15.15g %15.15g %15.15g %15.15g %15.15g\n",g1,g2,rfc,rf1,rf2);
    //    printf("g2 %15.15g %15.15g %15.15g %15.15g %15.15g %15.5g %15.5g\n",g2,b2,b1,( 1.0_rt-c3 ) ,2.0_rt*a1,b1,( 3.0_rt-2.0_rt*c2 ));
    ri1 = 0.5_rt / smc;
    ri2 = rf1 * smc;
    ri3 = 4.0_rt * rf2 * smc - 2.0_rt * ri2;
    ri4 = ri2 * ri2;

    for (int k = kts + 1; k <= kte; ++k) {
        dzk = 0.5_rt * (dz[k] + dz[k - 1]);
        afk = dz[k] / (dz[k] + dz[k - 1]);
        abk = 1.0_rt - afk;
        duz = ((u[k] - u[k - 1]) * (u[k] - u[k-1])) + ((v[k] - v[k - 1]) * (v[k] - v[k - 1]));
        duz = duz / (dzk * dzk);
        dtz = (thl[k] - thl[k - 1]) / dzk;
        dqz = (qw[k] - qw[k - 1]) / dzk;

        vtt = 1.0_rt + vt[k] * abk + vt[k - 1] * afk; // beta-theta in NN09, Eq. 39
        vqq = tv0  + vq[k] * abk + vq[k - 1] * afk; // beta-q
        dtq = vtt * dtz + vqq * dqz;
        // alternatively, use theta-v without the sgs clouds
        // dtq = (thetav[k] - thetav[k - 1]) / dzk;

        dtl[k] = dtz;
        dqw[k] = dqz;
        dtv[k] = dtq;

        gm[k] = duz;
        gh[k] = -dtq * gtr;

        // gradient richardson number
        ri = -gh[k] / std::max(duz, 1.0e-10_rt);
        // a2fac is needed for the canuto/kitamura mod
        if (ckmod == 1) {
            a2fac = 1.0_rt / (1.0_rt + std::max(ri, 0.0_rt));
        } else {
            a2fac = 1.0_rt;
        }
        rfc = g1 / (g1 + g2);
        f1 = b1 * (g1 - c1) + 3.0_rt * a2 * a2fac * (1.0_rt - c2) * (1.0_rt - c5) + 2.0_rt * a1 * (3.0_rt - 2.0_rt * c2);
        f2 = b1 * (g1 + g2) - 3.0_rt * a1 * (1.0_rt - c2);
        rf1 = b1 * (g1 - c1) / f1;
        rf2 = b1 * g1 / f2;
        smc = a1 / (a2 * a2fac) * f1 / f2;
        shc = 3.0_rt * (a2 * a2fac) * (g1 + g2);
        ri1 = 0.5_rt / smc;
        ri2 = rf1 * smc;
        ri3 = 4.0_rt * rf2 * smc - 2.0_rt * ri2;
        ri4 = ri2 * ri2;

        // flux richardson number
        rf = std::min(ri1 * (ri + ri2 - std::sqrt(ri * ri - ri3 * ri + ri4)), rfc);
        //      printf("rf    %15.15g %15.15g %15.15g %15.15g %15.15g %15.15g\n",rf2,ri2,ri3,smc,ri4);
        sh[k] = shc * (rfc - rf) / (1.0_rt - rf);
        sm[k] = smc * (rf1 - rf) / (rf2 - rf) * sh[k];
        //      printf("sm[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",sm[k],shc,rfc,rf,1.0_rt);
        //      printf("sh[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",sh[k],smc,rf1,rf,rf2);
    }
    //    exit(1);
}

// @}

// ==================================================================
//     subroutine  mym_length:
//
//     input variables:    see subroutine mym_initialize
//
//     output variables:   see subroutine mym_initialize
//
//     work arrays:
//       elt(nx,ny)      : length scale depending on the pbl depth    (m)
//       vsc(nx,ny)      : velocity scale q_c                       (m/s)
//                         at first, used for computing elt
//
//     note: the mixing lengths are meant to be calculated at the full-
//           sigmal levels (or interfaces between the model layers).
//
//>\ingroup gsd_mynn_edmf
// this subroutine calculates the mixing lengths.
void mym_length_cc(
    int kts, int kte, Real xland,
    const Real* dz, /*Real dx,*/ const Real* zw,
    Real rmo, Real flt, Real fltv, Real flq,
    const Real* vt, const Real* vq,
    const Real* u1, const Real* v1, const Real* qke,
    const Real* dtv,
    Real* el, // intent(out)
    Real zi, const Real* theta,
    Real* qkw, // intent(out)
    Real psig_bl, const Real* cldfra_bl1d,
    int bl_mynn_mixlength,
    const Real* edmf_w1, const Real* edmf_a1,
    // model constants:
    Real tv0, // ==p608 * tref
    Real gtr) // ==grav / tref
{
    Real qtke[kte+1];
    Real thetaw[kte+1];
    Real elblmin[kte+1];
    Real elblavg[kte+1];
    Real zi2, h1, h2, hs, elblmin0, elblavg0, cldavg;

    Real cns, alp1, alp2, alp3, alp4, alp5, alp6;
    const Real minzi = 300.0_rt;
    const Real maxdz = 750.0_rt;
    const Real mindz = 300.0_rt;
    //Real zslh = 100.0_rt;
    //Real csl = 2.0_rt;
    const Real qke_elb_min = 0.018_rt;

    int i, j, k;
    Real afk, abk, zwk, zwk1, dzk, qdz, vflx, bv, tau_cloud, wstar, elb, els, elf, el_stab, el_mf, el_stab_mf, elb_mf, pblh_plus_ent, uonset, ugrid, wt_u, el_les;
    const Real ctau = 1000._rt; // constant for tau_cloud

    const Real grav = 9.8100004196166992, karman = 0.4000000059604645;
    const Real twothirds = 0.6666666865348816, onethird = 0.3333333432674408;
    const Real qmin = 0.0_rt;

    switch(bl_mynn_mixlength) {

        /* Original MYNN Mixing Length + BouLac */
        case 0:
            cns = 2.7;
            alp1 = 0.23;
            alp2 = 1.0;
            alp3 = 5.0;
            alp4 = 100.0;
            alp5 = 0.3;

            // Impose limits on the height integration for elt and the transition layer depth
            zi2 = std::min(10000.0_rt, Real(zw[kte-2]));  // originally integrated to model top, not just 10 km.
            h1 = std::max(0.3_rt * Real(zi2), Real(mindz));
            h1 = std::min(Real(h1), Real(maxdz));         // 1/2 transition layer depth
            h2 = h1 / 2.0;                                  // 1/4 transition layer depth

            qkw[kts] = std::sqrt(std::max(Real(qke[kts]), qkemin));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0_rt - afk;
                qkw[k] = std::sqrt(std::max(Real(qke[k] * abk + qke[k-1] * afk), qkemin));
            }

            elt = 1.0e-5;
            vsc = 1.0e-5;

            // ** Strictly, zwk*h[i,j] -> ( zwk*h[i,j]+z0 ) **
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi2 + h1) {
                dzk = 0.5_rt * (dz[k] + dz[k-1]);
                qdz = std::max(Real(qkw[k] - qmin), 0.03_rt) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }

            elt = alp1 * elt / vsc;
            vflx = (vt[kts] + 1.0_rt) * flt + (vq[kts] + tv0) * flq;
            vsc = std::cbrt(gtr * elt * std::max(Real(vflx), 0.0_rt));

            // ** Strictly, el[i,k=0] is not zero. **
            el[kts] = 0.0;
            zwk1 = zw[kts+1];

            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];

                // ** Length scale limited by the buoyancy effect **
                if (dtv[k] > 0.0_rt) {
                    bv = std::sqrt(gtr * dtv[k]);
                    elb = alp2 * qkw[k] / bv * (1.0_rt + alp3 / alp2 * std::sqrt(vsc / (bv * elt)));
                    elf = alp2 * qkw[k] / bv;
                } else {
                    elb = 1.0e10;
                    elf = elb;
                }

                // ** Length scale in the surface layer **
                if (rmo > 0.0_rt) {
                    els = karman * zwk / (1.0_rt + cns * std::min(zwk * rmo, zmax));
                } else {
                    els = karman * zwk * std::pow(1.0_rt - alp4 * zwk * rmo, 0.2_rt);
                }

                Real wt = 0.5_rt * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5_rt;

                el[k] = std::min(elb / (elb / elt + elb / els + 1.0_rt), elf);
            }
            break;

        /* Nonlocal (using BouLac) Form of Mixing Length */
        case 1:
            ugrid = std::sqrt(u1[kts] * u1[kts] + v1[kts] * v1[kts]);
            uonset = 15.0;
            wt_u = (1.0_rt - std::min(std::max(Real(ugrid - uonset), 0.0_rt) / 30.0_rt, 0.5_rt));
            cns = 3.5;
            alp1 = 0.23;
            alp2 = 0.3;
            alp3 = 2.5_rt * wt_u; // taper off buoyancy enhancement in shear-driven pbls
            alp4 = 5.0;
            alp5 = 0.3;
            alp6 = 50.0;

            // Impose limits on the height integration for elt and the transition layer depth
            zi2 = std::max(Real(zi), 300._rt);
            h1 = std::max(Real(0.3_rt * zi2), 300.0_rt);
            h1 = std::min(Real(h1), 600.0_rt);
            h2 = h1 / 2.0_rt;

            qtke[kts] = std::max(Real(0.5_rt * qke[kts]), 0.5_rt*qkemin);
            thetaw[kts] = theta[kts];
            qkw[kts] = std::sqrt(std::max(Real(qke[kts]), qkemin));

            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0_rt - afk;
                qkw[k] = std::sqrt(std::max(Real(qke[k] * abk + qke[k-1] * afk), qkemin));
                qtke[k] = std::max(0.5_rt * (qkw[k] * qkw[k]),0.005_rt);
                thetaw[k] = theta[k] * abk + theta[k-1] * afk;
            }

            elt = 1.0e-5;
            vsc = 1.0e-5;

            // ** Strictly, zwk*h[i,j] -> ( zwk*h[i,j]+z0 ) **
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi2 + h1) {
                dzk = 0.5_rt * (dz[k] + dz[k-1]);
                qdz = std::min(std::max(Real(qkw[k] - qmin), 0.01_rt), 30.0_rt) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }

            elt = std::min(std::max(Real(alp1 * elt / vsc), 8.0_rt), 400.0_rt);
            vflx = fltv;
            vsc = std::cbrt(gtr * elt * std::max(Real(vflx), 0.0_rt));

            // ** Strictly, el[i,j,0] is not zero **
            el[kts] = 0.0;
            zwk1 = zw[kts+1];

            // Compute BouLac mixing length
            boulac_length_cc(kts,kte,zw,dz,qtke,thetaw,elblmin,elblavg,gtr);

            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];

                // ** Length scale is limited by the buoyancy effect **
                if (dtv[k] > 0.0_rt) {
                    bv = std::max(std::sqrt(gtr * dtv[k]), 0.0001_rt);
                    elb = std::max(Real(alp2 * std::max(qkw[k],qke_elb_min)), Real(alp6 * edmf_a1[k-1] * edmf_w1[k-1])) / bv * (1.0_rt + alp3 * std::sqrt(vsc / (bv * elt)));
                    elb = std::min(elb, zwk);
                    elf = 1.0_rt * std::max(qkw[k],qke_elb_min) / bv;
                    elblavg[k] = std::max(Real(elblavg[k]), Real(alp6 * edmf_a1[k-1] * edmf_w1[k-1] / bv));
                } else {
                    elb = 1.0e10;
                    elf = elb;
                }
                if (rmo > 0.0_rt) {
                    els = karman * zwk / (1.0_rt + cns * std::min(zwk * rmo, zmax));
                } else {
                    els = karman * zwk * std::pow(1.0_rt - alp4 * zwk * rmo, 0.2_rt);
                }
                Real wt = 0.5_rt * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5_rt;
                el[k] = std::sqrt((els*els)/(1.0_rt + (els*els)/(elt*elt)));
                el[k] = std::min(el[k], elb);
                el[k] = std::min(el[k], elf);
                el[k] = el[k]*(1.0_rt-wt) + alp5*elblavg[k]*wt;
                el[k] = el[k] * psig_bl;
            }
            break;

        /* Local (mostly) Mixing Length Formulation */
        case 2:
            uonset = 3.5_rt + dz[kts] * 0.1_rt;
            ugrid = std::sqrt(u1[kts] * u1[kts] + v1[kts] * v1[kts]);
            cns = 3.5;
            alp1 = 0.22;
            alp2 = 0.30;
            alp3 = 2.0;
            alp4 = 5.0;
            alp5 = alp2; // like alp2, but for free atmosphere
            alp6 = 50.0; // used for MF mixing length
            zi2 = std::max(Real(zi), Real(minzi));

            // Impose limits on the height integration for elt and the transition layer depth
            h1 = std::max(Real(0.3_rt * zi2), 300.0_rt);
            h1 = std::min(Real(h1), 600.0_rt);           // 1/2 transition layer depth
            h2 = h1 * 0.5_rt;                             // 1/4 transition layer depth

            qtke[kts] = std::max(Real(0.5_rt * qke[kts]), 0.5_rt*qkemin);
            qkw[kts] = std::sqrt(std::max(Real(qke[kts]), qkemin));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0_rt - afk;
                qkw[k] = std::sqrt(std::max(Real(qke[k] * abk + qke[k-1] * afk), qkemin));
                qtke[k] = 0.5_rt * qkw[k] * qkw[k];
            }

            elt = 1.0e-5;
            vsc = 1.0e-5;

            // ** Strictly, zwk*h[i,j] -> ( zwk*h[i,j]+z0 ) **
            pblh_plus_ent = std::max(zi+h1, 100._rt);
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= pblh_plus_ent) {
                dzk = 0.5_rt * (dz[k] + dz[k-1]);
                qdz = std::min(std::max(Real(qkw[k] - qmin), 0.03_rt), 30.0_rt) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }

            elt = std::min(std::max(Real(alp1 * elt / vsc), 10.0_rt), 400.0_rt);
            vflx = fltv;
            vsc = std::cbrt(gtr * elt * std::max(Real(vflx), 0.0_rt));

            // ** Strictly, el[i,j,0] is not zero **
            el[kts] = 0.0;
            zwk1 = zw[kts+1];

            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];
                dzk = 0.5_rt * (dz[k] + dz[k-1]);
                cldavg = 0.5_rt * (cldfra_bl1d[k-1] + cldfra_bl1d[k]);

                // ** Length scale limited by the buoyancy effect **
                if (dtv[k] > 0.0_rt) {
                    bv = std::max(std::sqrt(gtr * dtv[k]), 0.001_rt);
                    elb_mf = std::max(Real(alp2 * qkw[k]), Real(alp6 * edmf_a1[k-1] * edmf_w1[k-1]) / bv * (1.0_rt + alp3 * std::sqrt(vsc / (bv * elt))));
                    elb = std::min(std::max(Real(alp5 * qkw[k]), Real(alp6 * edmf_a1[k] * edmf_w1[k]) / bv), Real(zwk));

                    wstar = 1.25_rt * std::cbrt(gtr * zi * std::max(Real(vflx), 1.0e-4_rt));
                    tau_cloud = std::min(std::max(Real(ctau * wstar / grav), 30.0_rt), 150.0_rt);

                    // minimize influence of surface heat flux on tau far away from the PBLH
                    Real wt = 0.5_rt * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5_rt;
                    tau_cloud = tau_cloud * (1.0_rt - wt) + 50.0_rt * wt;
                    elf = std::min(std::max(Real(tau_cloud * std::sqrt(std::min(Real(qtke[k]), 40.0_rt))), Real(alp6 * edmf_a1[k] * edmf_w1[k] / bv)), Real(zwk));
                } else {
                    /*
                     * use version in development for RAP/HRRR 2016
                     * JAYMES-
                     * tau_cloud is an eddy turnover timescale;
                     * see Teixeira and Cheinet (2004), Eq. 1, and
                     * Cheinet and Teixeira (2003), Eq. 7.  The
                     * coefficient 0.5 is tuneable. Expression in
                     * denominator is identical to vsc (a convective
                     * velocity scale), except that elt is replaced
                     * by zi, and zero is replaced by 1.0e-4 to
                     * prevent division by zero.
                     */
                    wstar = 1.25_rt * std::cbrt(gtr * zi * std::max(Real(vflx), 1.0e-4_rt));
                    tau_cloud = std::min(std::max(Real(ctau * wstar / grav), 50.0_rt), 200.0_rt);

                    // minimize influence of surface heat flux on tau far away from the PBLH
                    Real wt = 0.5_rt * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5_rt;
                    tau_cloud = tau_cloud * (1.0_rt - wt) + std::max(100.0_rt, dzk * 0.25_rt) * wt;

                    elb = std::min(tau_cloud * std::sqrt(std::min(qtke[k], 40.0_rt)), zwk);
                    elf = elb;
                    elb_mf = elb;
                }
                elf = elf / (1.0_rt + (elf / 800.0_rt));
                elb_mf = std::max(Real(elb_mf), 0.01_rt);

                // ** Length scale in the surface layer **
                if (rmo > 0.0_rt) {
                    els = karman * zwk / (1.0_rt + cns * std::min(zwk * rmo, zmax));
                } else {
                    els = karman * zwk * std::pow(1.0_rt - alp4 * zwk * rmo, 0.2_rt);
                }

                // ** Now blend the mixing length scales **
                Real wt = 0.5_rt * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5_rt;

                // try squared-blending
                el[k] = std::sqrt(els * els / (1.0_rt + (els * els / (elt * elt)) + (els * els / (elb_mf * elb_mf))));
                el[k] = el[k] * (1.0_rt - wt) + elf * wt;

                // include scale-awareness. For now, use simple asymptotic kz -> 12 m (should be ~ dz)
                el_les = std::min(els / (1.0_rt + (els/12._rt)), elb_mf);
                el[k] = el[k] * psig_bl + (1.0_rt - psig_bl) * el_les;
            }
            break;
    }
}




// called from driver
void moisture_check_cc(int kte, Real delt, Real* dp, const Real* exner,
                    Real* qv, Real* qc, Real* qi, Real* qs, Real* th,
                    Real* dqv, Real* dqc, Real* dqi, Real* dqs, Real* dth,
                    Real xlvcp, Real xlscp) {

    // constants (assuming xlvcp and xlscp are defined elsewhere)
    const Real qvmin = 1e-20, qcmin = 0.0, qimin = 0.0;
    Real dqv2;

    for (int k = kte; k >= 0; --k) { // from the top to the surface
        Real dqc2 = std::max(0.0_rt, qcmin - qc[k]);
        Real dqi2 = std::max(0.0_rt, qimin - qi[k]);
        Real dqs2 = std::max(0.0_rt, qimin - qs[k]);

        // fix tendencies
        dqc[k] += dqc2 / delt;
        dqi[k] += dqi2 / delt;
        dqs[k] += dqs2 / delt;
        dqv[k] -= (dqc2 + dqi2 + dqs2) / delt;
        dth[k] += xlvcp / exner[k] * (dqc2 / delt) + xlscp / exner[k] * ((dqi2 + dqs2) / delt);

        // update species
        qc[k] += dqc2;
        qi[k] += dqi2;
        qs[k] += dqs2;
        qv[k] -= dqc2 + dqi2 + dqs2;
        th[k] += xlvcp / exner[k] * dqc2 + xlscp / exner[k] * (dqi2 + dqs2);

        // then fix qv
        Real dqv2 = std::max(0.0_rt, qvmin - qv[k]);
        dqv[k] += dqv2 / delt;
        qv[k] += dqv2;
        if (k != 0) {
            qv[k-1] -= dqv2 * dp[k] / dp[k-1];
            dqv[k-1] -= dqv2 * dp[k] / dp[k-1] / delt;
        }
        qv[k] = std::max(Real(qv[k]), Real(qvmin));
        qc[k] = std::max(Real(qc[k]), Real(qcmin));
        qi[k] = std::max(Real(qi[k]), Real(qimin));
        qs[k] = std::max(Real(qs[k]), Real(qimin));
    }

        Real sum = 0.0;
    Real aa, dum;

    // only execute if dqv2 > 1.e-20, which indicates adjustment was made at the top layer
    if(dqv2 > 1e-20) {
        for (int k = 0; k <= kte; ++k) { // loop through all layers
            if (qv[k] > 2.0_rt * qvmin) {
                sum += qv[k] * dp[k];
            }
        }

        aa = dqv2 * dp[0] / std::max(1.e-20_rt, sum); // adjust for 1-based indexing with dp[0]

        if (aa < 0.5_rt) {
            for (int k = 0; k <= kte; ++k) { // loop through all layers again
                if (qv[k] > 2.0_rt * qvmin) {
                    dum = aa * qv[k];
                    qv[k] -= dum;
                    dqv[k] -= dum / delt;
                }
            }
        } else {
            // for testing purposes only (not yet found in any output):
            // std::cout << "full moisture conservation is impossible" << std::endl;
        }
    }

}

/*
! ==================================================================
!     subroutine  mym_predict:
!
!     input variables:    see subroutine mym_initialize and turbulence
!       qke(nx,nz,ny) : qke at (n)th time level
!       tsq, ...cov     : ditto
!
!     output variables:
!       qke(nx,nz,ny) : qke at (n+1)th time level
!       tsq, ...cov     : ditto
!
!     work arrays:
!       qkw(nx,nz,ny)   : q at the center of the grid boxes        (m/s)
!       bp (nx,nz,ny)   : = 1/2*f,     see below
!       rp (nx,nz,ny)   : = p-1/2*f*q, see below
!
!     # the equation for a turbulent quantity q can be expressed as
!          dq/dt + ah + av = dh + dv + p - f*q,                      (1)
!       where a is the advection, d the diffusion, p the production,
!       f*q the dissipation and h and v denote horizontal and vertical,
!       respectively. if q is q^2, f is 2q/b_1l.
!       using the crank-nicholson scheme for av, dv and f*q, a finite
!       difference equation is written as
!          q{n+1} - q{n} = dt  *( dh{n}   - ah{n}   + p{n} )
!                        + dt/2*( dv{n}   - av{n}   - f*q{n}   )
!                        + dt/2*( dv{n+1} - av{n+1} - f*q{n+1} ),    (2)
!       where n denotes the time level.
!       when the advection and diffusion terms are discretized as
!          dt/2*( dv - av ) = a(k)q(k+1) - b(k)q(k) + c(k)q(k-1),    (3)
!       eq.(2) can be rewritten as
!          - a(k)q(k+1) + [ 1 + b(k) + dt/2*f ]q(k) - c(k)q(k-1)
!                 = q{n} + dt  *( dh{n}   - ah{n}   + p{n} )
!                        + dt/2*( dv{n}   - av{n}   - f*q{n}   ),    (4)
!       where q on the left-hand side is at (n+1)th time level.
!
!       in this subroutine, a(k), b(k) and c(k) are obtained from
!       subprogram coefvu and are passed to subprogram tinteg via
!       common. 1/2*f and p-1/2*f*q are stored in bp and rp,
!       respectively. subprogram tinteg solves eq.(4).
!
!       modify this subroutine according to your numerical integration
!       scheme (program).
!
!-------------------------------------------------------------------
!>\ingroup gsd_mynn_edmf
!! this subroutine predicts the turbulent quantities at the next step.
*/
void mym_predict_cc(
    int& kts, int& kte,
    Real& closure,
    Real& delt,
    Real* dz,
    Real& ust, Real& flt, Real& flq, Real& pmz, Real& phh,
    Real* el, Real* dfq, Real* rho,
    /* begin intent(inout) */
    Real* pdk, Real* pdt, Real* pdq, Real* pdc,
    Real* qke, Real* tsq, Real* qsq, Real* cov,
    Real* s_aw, Real* s_awqke,
    /* end intent(inout) */
    int& bl_mynn_edmf_tke,
    /* begin if tke_budget==1, intent(out) */
    Real* qwt1d, Real* qdiss1d, int& tke_budget,
    /* end tke_budget */
    // model constants:
    Real& xlvcp,  // == xlv/cp       (kind_phys); xlv=2.5e6, cp=7.*287./2. (real)
    Real& xlscp,  // == (xlv+xlf)/cp (kind_phys); xlf=3.50e5               (real)
    Real& karman) // == 0.4          (real)
{
    Real vkz, pdk1, phm, pdt1, pdq1, pdc1, b1l, b2l, onoff;
    Real dtz[kte-kts+1];
    Real a[kte-kts+1];
    Real b[kte-kts+1];
    Real c[kte-kts+1];
    Real d[kte-kts+1];
    Real x[kte-kts+1];
    Real rhoinv[kte-kts+1];
    Real rhoz[kte-kts+2];
    Real kqdz[kte-kts+2];
    Real kmdz[kte-kts+2];
    Real qkw[kte-kts+1];
    Real bp[kte-kts+1];
    Real rp[kte-kts+1];
    Real df3q[kte-kts+1];
    Real tke_up[kte-kts+1];
    Real dzinv[kte-kts+1];

    // regulate the momentum mixing from the mass-flux scheme (on or off)
    if (bl_mynn_edmf_tke == 0) {
        onoff = 0.0;
    } else {
        onoff = 1.0;
    }

    // calculate vkz
    vkz = karman * 0.5_rt * dz[kts];

    // calculate df3q and dtz
    for (int k = kts; k <= kte; k++) {
        qkw[k] = std::sqrt(std::max(qke[k], 0.0_rt));
        df3q[k] = sqfac * dfq[k];
        dtz[k] = delt / dz[k];
    }

    // prepare "constants" for diffusion equation
    rhoz[kts] = rho[kts];
    rhoinv[kts] = 1.0_rt / rho[kts];
    kqdz[kts] = rhoz[kts] * df3q[kts];
    kmdz[kts] = rhoz[kts] * dfq[kts];
    for (int k = kts+1; k <= kte; k++) {
        rhoz[k] = (rho[k] * dz[k-1] + rho[k-1] * dz[k]) / (dz[k-1] + dz[k]);
        rhoz[k] = std::max(rhoz[k], 1e-4_rt);
        rhoinv[k] = 1.0_rt / std::max(rho[k], 1e-4_rt);
        kqdz[k] = rhoz[k] * df3q[k];
        kmdz[k] = rhoz[k] * dfq[k];
    }
    rhoz[kte+1] = rhoz[kte];
    kqdz[kte+1] = rhoz[kte+1] * df3q[kte];
    kmdz[kte+1] = rhoz[kte+1] * dfq[kte];

    // stability criteria for mf
    for (int k = kts+1; k <= kte; k++) {
       kqdz[k] = std::max(kqdz[k],  0.5_rt* s_aw[k]);
       kqdz[k] = std::max(kqdz[k], -0.5_rt*(s_aw[k]-s_aw[k+1]));
       kmdz[k] = std::max(kmdz[k],  0.5_rt* s_aw[k]);
       kmdz[k] = std::max(kmdz[k], -0.5_rt*(s_aw[k]-s_aw[k+1]));
    }

    // calculate pdk1, phm, pdt1, pdq1, pdc1
    pdk1 = 2.0_rt * (ust*ust*ust) * pmz / vkz;
    phm = 2.0_rt / ust * phh / vkz;
    pdt1 = phm * flt * flt;
    pdq1 = phm * flq * flq;
    pdc1 = phm * flt * flq;

    // calculate pdk, pdt, pdq, pdc
    pdk[kts] = pdk1 - pdk[kts+1];
    pdt[kts] = pdt[kts+1];
    pdq[kts] = pdq[kts+1];
    pdc[kts] = pdc[kts+1];

    // prediction of twice the turbulent kinetic energy
    for (int k = kts; k <= kte-1; k++) {
        b1l = b1 * 0.5_rt * (el[k+1] + el[k]);
        bp[k] = 2.0_rt * qkw[k] / b1l;
        rp[k] = pdk[k+1] + pdk[k];
    }
    for (int k = kts; k <= kte-1; k++) {
        a[k] = -dtz[k] * kqdz[k] * rhoinv[k]
             + 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k] * onoff;
        b[k] = 1.0_rt + dtz[k] * (kqdz[k] + kqdz[k+1]) * rhoinv[k]
             + 0.5_rt * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) * onoff
             + bp[k] * delt;
        c[k] = -dtz[k] * kqdz[k+1] * rhoinv[k]
             - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff;
        d[k] = rp[k] * delt + qke[k]
             + dtz[k] * rhoinv[k] * (s_awqke[k] - s_awqke[k+1]) * onoff;
    }
    /*
      for (int k = kts; k <= kte-1; k++) {
        a[k] = -dtz[k] * df3q[k] + 0.5_rt * dtz[k] * s_aw[k] * onoff;
        b[k] = 1.0_rt + dtz[k] * (df3q[k] + df3q[k+1]) + 0.5_rt * dtz[k] * (s_aw[k] - s_aw[k+1]) * onoff + bp[k] * delt;
        c[k] = -dtz[k] * df3q[k+1] - 0.5_rt * dtz[k] * s_aw[k+1] * onoff;
        d[k] = rp[k] * delt + qke[k] + dtz[k] * (s_awqke[k] - s_awqke[k+1]) * onoff;
        }*/
    a[kte] = 0.0;
    b[kte] = 1.0;
    c[kte] = 0.0;
    d[kte] = qke[kte];

    tridiag2_cc(kte, a, b, c, d, x);

    for (int k = kts; k <= kte; k++) {
        qke[k] = std::max(x[k], qkemin);
        qke[k] = std::min(qke[k], 150.0_rt);
    }

    // tke budget
    if (tke_budget == 1) {
        Real tke_up[kte-kts+1];
        Real dzinv[kte-kts+1];

        // tke vertical transport
        for (int k=kts; k <=kte; k++)
        {
            tke_up[k] = 0.5_rt * qke[k];
            dzinv[k] = 1.0_rt / dz[k];
        }

        qwt1d[kts] = dzinv[kts] * ((kqdz[kts+1] * (tke_up[kts+1] - tke_up[kts]) - (kqdz[kts] * tke_up[kts])) + 0.5_rt * rhoinv[kts] * (s_aw[kts+1] * tke_up[kts+1] + ((s_aw[kts+1] - s_aw[kts]) * tke_up[kts]) + (s_awqke[kts] - s_awqke[kts+1])) * onoff);
        for (int k = kts+1; k <= kte-1; k++) {
            qwt1d[k] = dzinv[k] * ((kqdz[k+1] * (tke_up[k+1] - tke_up[k]) - (kqdz[k] * (tke_up[k] - tke_up[k-1]))) + 0.5_rt * rhoinv[k] * (s_aw[k+1] * tke_up[k+1] + ((s_aw[k+1] - s_aw[k]) * tke_up[k]) - (s_aw[k] * tke_up[k-1]) + (s_awqke[k] - s_awqke[k+1])) * onoff);
        }
        qwt1d[kte] = dzinv[kte] * (-(kqdz[kte] * (tke_up[kte] - tke_up[kte-1])) + 0.5_rt * rhoinv[kte] * (-(s_aw[kte] * tke_up[kte]) - (s_aw[kte] * tke_up[kte-1]) + s_awqke[kte]) * onoff);

        // tke dissipation rate
        for (int k=kts; k <=kte; k++)
        {
            qdiss1d[k] = bp[k] * tke_up[k];
        }
    }

    if (closure > 2.5) {
        // prediction of the moisture variance
        for (int k = kts; k <= kte-1; k++) {
            b2l = b2 * 0.5_rt * (el[k+1] + el[k]);
            bp[k] = 2.0_rt * qkw[k] / b2l;
            rp[k] = pdq[k+1] + pdq[k];
        }
        for (int k = kts; k <= kte-1; k++) {
            a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
            b[k] = 1.0_rt + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + bp[k] * delt;
            c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k];
            d[k] = rp[k] * delt + qsq[k];
        }
        a[kte] = -1.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = 0.0;
        tridiag2_cc(kte, a, b, c, d, x);
        for (int k = kts; k <= kte; k++) {
            qsq[k] = std::max(x[k], 1e-17_rt);
        }
    } else {
        // level 2.5_rt - use level 2 diagnostic
        for (int k = kts; k <= kte-1; k++) {
            if (qkw[k] <= 0.0_rt) {
                b2l = 0.0;
            } else {
                b2l = b2 * 0.25_rt * (el[k+1] + el[k]) / qkw[k];
            }
            qsq[k] = b2l * (pdq[k+1] + pdq[k]);
        }
        qsq[kte] = qsq[kte-1];
    }

    if (closure >= 3.0_rt) {
        // prediction of the temperature variance
        for (int k = kts; k <= kte-1; k++) {
            b2l = b2 * 0.5_rt * (el[k+1] + el[k]);
            bp[k] = 2.0_rt * qkw[k] / b2l;
            rp[k] = pdt[k+1] + pdt[k];
        }
        for (int k = kts; k <= kte-1; k++) {
            a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
            b[k] = 1.0_rt + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + bp[k] * delt;
            c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k];
            d[k] = rp[k] * delt + tsq[k];
        }
        a[kte] = 0.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = tsq[kte];
        tridiag2_cc(kte, a, b, c, d, x);
        for (int k = kts; k <= kte; k++) {
            tsq[k] = x[k];
        }

        // prediction of the temperature-moisture covariance
        for (int k = kts; k <= kte-1; k++) {
            b2l = b2 * 0.5_rt * (el[k+1] + el[k]);
            bp[k] = 2.0_rt * qkw[k] / b2l;
            rp[k] = pdc[k+1] + pdc[k];
        }
        for (int k = kts; k <= kte-1; k++) {
            a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
            b[k] = 1.0_rt + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + bp[k] * delt;
            c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k];
            d[k] = rp[k] * delt + cov[k];
        }
        a[kte] = 0.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = 0.0;
        tridiag2_cc(kte, a, b, c, d, x);
        for (int k = kts; k <= kte; k++) {
            cov[k] = x[k];
        }
    } else {
        // not level 3 - default to level 2 diagnostic
        for (int k = kts; k <= kte-1; k++) {
            if (qkw[k] <= 0.0_rt) {
                b2l = 0.0;
            } else {
                b2l = b2 * 0.25_rt * (el[k+1] + el[k]) / qkw[k];
            }
            tsq[k] = b2l * (pdt[k+1] + pdt[k]);
            cov[k] = b2l * (pdc[k+1] + pdc[k]);
        }
        tsq[kte] = tsq[kte-1];
        cov[kte] = cov[kte-1];
    }
}

void mynn_mix_chem_cc(int kts, int kte, int i,
                   Real delt, Real* dz, Real pblh,
                   int nchem, int kdvel, int ndvel,
                   Real** chem1, Real* vd1,
                   Real* rho,
                   Real flt, Real* tcd, Real* qcd,
                   Real* dfh,
                   Real* s_aw, Real** s_awchem,
                   Real emis_ant_no, Real frp, int rrfs_sd, int enh_mix) {

    // local vars
    Real dtz[kte - kts + 1];
    Real a[kte - kts + 1], b[kte - kts + 1], c[kte - kts + 1], d[kte - kts + 1], x[kte - kts + 1];
    Real dztop = 0.5_rt * (dz[kte - 1] + dz[kte - 2]);
    for (int k = kts; k <= kte; ++k) {
        dtz[k - kts] = delt / dz[k - 1];
    }
    // prepare "constants" for diffusion equation.
    Real rhoz[kte - kts + 2], khdz[kte - kts + 2], rhoinv[kte - kts + 1];
    rhoz[0] = rho[kts - 1];
    rhoinv[0] = 1.0_rt / rho[kts - 1];
    khdz[0] = rhoz[0] * dfh[kts - 1];
    for (int k = kts + 1; k <= kte; ++k) {
        rhoz[k - kts] = (rho[k - 1] * dz[k - 2] + rho[k - 2] * dz[k - 1]) / (dz[k - 2] + dz[k - 1]);
        rhoz[k - kts] = std::max(Real(rhoz[k - kts]), 1e-4_rt);
        rhoinv[k - kts] = 1.0_rt / std::max(Real(rho[k - 1]), 1e-4_rt);
        Real dzk = 0.5_rt * (dz[k - 1] + dz[k - 2]);
        khdz[k - kts] = rhoz[k - kts] * dfh[k - 1];
    }
    rhoz[kte - kts + 1] = rhoz[kte - kts];
    khdz[kte - kts + 1] = rhoz[kte - kts + 1] * dfh[kte - 1];
    // stability criteria for mf
    for (int k = kts + 1; k <= kte - 1; ++k) {
        khdz[k - kts] = std::max(Real(khdz[k - kts]), Real(0.5_rt * s_aw[k - kts]));
        khdz[k - kts] = std::max(Real(khdz[k - kts]), Real(-0.5_rt * (s_aw[k - kts] - s_aw[k - kts + 1])));
    }
    // enhanced mixing over fires
    if (rrfs_sd==1 && enh_mix==1) {
        for (int k = kts + 1; k <= kte - 1; ++k) {
            Real khdz_old = khdz[k - kts];
            Real khdz_back = pblh * 0.15_rt / dz[k - 1];
            // modify based on anthropogenic emissions of no and frp
            if (pblh < pblh_threshold) {
                if (emis_ant_no > no_threshold) {
                    khdz[k - kts] = std::max(1.1_rt * Real(khdz[k - kts]), std::sqrt((emis_ant_no / no_threshold)) / dz[k - 1] * rhoz[k - kts]);
                }
                if (frp > frp_threshold) {
                    int kmaxfire = std::ceil(std::log(frp));
                    khdz[k - kts] = std::max(Real(1.1_rt * khdz[k - kts]), Real((1.0_rt - k / (kmaxfire * 2.0_rt)) * (std::pow(std::log(frp), 2.0_rt) - 2.0_rt * std::log(frp)) / dz[k - 1] * rhoz[k - kts]));
                }
            }
        }
    }
    // mixing of chemical species
    for (int ic = 0; ic < nchem; ++ic) {
        int k = kts;
        a[0] = -dtz[0] * khdz[0] * rhoinv[0];
        b[0] = 1.0_rt + dtz[0] * (khdz[1] + khdz[0]) * rhoinv[0] - 0.5_rt * dtz[0] * rhoinv[0] * s_aw[1];
        c[0] = -dtz[0] * khdz[1] * rhoinv[0] - 0.5_rt * dtz[0] * rhoinv[0] * s_aw[1];
        d[0] = chem1[k - 1][ic] - dtz[0] * vd1[ic] * chem1[k - 1][ic] - dtz[0] * rhoinv[0] * s_awchem[1][ic];
        for (k = kts + 1; k <= kte - 1; ++k) {
            a[k - kts] = -dtz[k - kts] * khdz[k - kts] * rhoinv[k - kts] + 0.5_rt * dtz[k - kts] * rhoinv[k - kts] * s_aw[k - kts];
            b[k - kts] = 1.0_rt + dtz[k - kts] * (khdz[k - kts] + khdz[k - kts + 1]) * rhoinv[k - kts] + 0.5_rt * dtz[k - kts] * rhoinv[k - kts] * (s_aw[k - kts] - s_aw[k - kts + 1]);
            c[k - kts] = -dtz[k - kts] * khdz[k - kts + 1] * rhoinv[k - kts] - 0.5_rt * dtz[k - kts] * rhoinv[k - kts] * s_aw[k - kts + 1];
            d[k - kts] = chem1[k - 1][ic] + dtz[k - kts] * rhoinv[k - kts] * (s_awchem[k - kts][ic] - s_awchem[k - kts + 1][ic]);
        }
        // prescribed value at top
        a[kte - kts] = 0.0;
        b[kte - kts] = 1.0;
        c[kte - kts] = 0.0;
        d[kte - kts] = chem1[kte - 1][ic];
        tridiag3_cc(kte, a, b, c, d, x);
        for (k = kts; k <= kte; ++k) {
            chem1[k - 1][ic] = x[k - kts];
        }
    }
}


// ==================================================================
//>\ingroup gsd_mynn_edmf
// this subroutine solves for tendencies of u, v, \f$\theta\f$, qv,
// qc, and qi
void mynn_tendencies_cc(const int& kts,const int& kte, const Real & delt,
                                   /*in*/ const Real* dz,
                                   /*in*/ const Real* rho,
                        /*in*/ const Real* u, const Real* v, const Real* th, const Real* tk, const Real* qv,
                                   /*in*/ const Real* qc, const Real* qi, const Real* qs, const Real* qnc, const Real* qni,
                        /*in*/ const Real* psfc,const Real* p,const Real* exner,
                                   /*inout*/ Real* thl, Real* sqv, Real* sqc, Real* sqi,
                                   /*inout*/ Real* sqs, Real* sqw, Real* qnwfa, Real* qnifa, Real* qnbca, Real* ozone,
                                   /*in*/ Real* ust, const Real & flt,const Real & flq,const Real & flqv,const Real & flqc,const Real & wspd, const Real & uoce, const Real & voce,
                                   /*in*/ const Real* tcd,const Real* qcd,
                                   /*inout*/ Real* dfm, Real* dfh,
                                   /*inout*/ Real* du, Real* dv, Real* dth,
                                   /*inout*/ Real* dqv, Real* dqc, Real* dqi,
                                   /*inout*/ Real* dqs, Real* dqnc, Real* dqni,
                                   /*inout*/ Real* dqnwfa, Real* dqnifa, Real* dqnbca,
                                   /*inout*/ Real* dozone,
                                   /*in*/ const Real* diss_heat,
                        /*in*/ const Real* s_aw, const Real* s_awthl, const Real* s_awqt, const Real* s_awqv, const Real* s_awqc, const Real* s_awu, const Real* s_awv, const Real* s_awqnc, const Real* s_awqni, const Real* s_awqnwfa, const Real* s_awqnifa, const Real* s_awqnbca, const Real* sd_aw, const Real* sd_awthl, const Real* sd_awqt, const Real* sd_awqv, const Real* sd_awqc, const Real* sd_awu, const Real* sd_awv,
                                   const Real* sub_thl,const Real* sub_sqv,const Real* sub_u,const Real* sub_v,
                                   const Real* det_thl,const Real* det_sqv,const Real* det_sqc,const Real* det_u,const Real* det_v,
                                   /*logical turned into int */const int& flag_qc, const int& flag_qi, const int& flag_qnc, const int& flag_qni, const int& flag_qs, const int& flag_qnwfa, const int& flag_qnifa, const int& flag_qnbca, const int& flag_ozone,
                        const int & bl_mynn_cloudmix, const int & bl_mynn_mixqt, int & bl_mynn_edmf_mom,  const int & bl_mynn_mixscalars, /* new */const int& debug_code,const Real& r_d,const Real& p608,const Real& ep_2,const Real& ep_3,const Real& tv0,const Real& xlv,const Real& xlvcp,const Real& xlscp) {
  /*
    printf("thl\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",thl[k]);
    printf("\n");

    printf("sqv\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",sqv[k]);
    printf("\n");

    printf("sqc\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",sqc[k]);
    printf("\n");

    printf("sqi\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",sqi[k]);
    printf("\n");

    printf("sqs\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",sqs[k]);
    printf("\n");

    printf("sqw\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",sqw[k]);
    printf("\n");

    printf("qnwfa\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",qnwfa[k]);
    printf("\n");

    printf("qnifa\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",qnifa[k]);
    printf("\n");

    printf("qnbca\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",qnbca[k]);
    printf("\n");

    printf("ozone\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",ozone[k]);
    printf("\n");

    printf("dfm\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",dfm[k]);
    printf("\n");

    printf("dfh\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",dfh[k]);
    printf("\n");

    printf("du\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",du[k]);
    printf("\n");

    printf("dv\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",dv[k]);
    printf("\n");

    printf("dth\n");
    for (int k = kts; k <= kte; k++)
      printf("%g ",dth[k]);
    printf("\n");
  */

    Real nonloc = 1.0;

    Real dztop = 0.5_rt * (dz[kte-1] + dz[kte-2]);
    Real onoff = bl_mynn_edmf_mom;
    onoff = (onoff == 0) ? 0.0_rt : 1.0;
    Real rhosfc = *psfc / (r_d * (tk[kts] + p608 * qv[kts]));

    Real dtz[kte+2];
    Real dfhc[kte+2];
    Real dfmc[kte+2];
    Real delp[kte+2];
    Real sqv2[kte+2];
    Real sqc2[kte+2];
    Real sqs2[kte+2];
    Real sqi2[kte+2];
    Real qnc2[kte+2];
    Real qni2[kte+2];
    Real rhoz[kte+2];
    Real khdz[kte+2];
    Real kmdz[kte+2];
    Real rhoinv[kte+2];
    Real sqw2[kte+2];
    Real a[kte+2];
    Real b[kte+2];
    Real c[kte+2];
    Real d[kte+2];
    Real x[kte+2];
    Real qvflux;
    Real ust_v=*ust;

    dtz[kts] = delt / dz[kts];
    rhoz[kts] = rho[kts];
    rhoinv[kts] = 1.0_rt / rho[kts];
    khdz[kts] = rhoz[kts] * dfh[kts];
    kmdz[kts] = rhoz[kts] * dfm[kts];
    delp[kts] = *psfc - (p[kts+1] * dz[kts] + p[kts] * dz[kts+1]) / (dz[kts] + dz[kts+1]);

    for (int k = kts+1; k <= kte; k++) {
        dtz[k] = delt / dz[k];
        rhoz[k] = (rho[k] * dz[k-1] + rho[k-1] * dz[k]) / (dz[k-1] + dz[k]);
        rhoz[k] = std::max(Real(rhoz[k]), 1e-4_rt);
        rhoinv[k] = 1.0_rt / std::max(Real(rho[k]), 1e-4_rt);
        Real dzk = 0.5_rt * (dz[k] + dz[k-1]);
        khdz[k] = rhoz[k] * dfh[k];
        kmdz[k] = rhoz[k] * dfm[k];
    }

    for (int k = kts+1; k <= kte-1; k++) {
        delp[k] = (p[k] * dz[k-1] + p[k-1] * dz[k]) / (dz[k] + dz[k-1]) - (p[k+1] * dz[k] + p[k] * dz[k+1]) / (dz[k] + dz[k+1]);
    }
    delp[kte] = delp[kte-1];
    rhoz[kte+1] = rhoz[kte];
    khdz[kte+1] = rhoz[kte+1] * dfh[kte];
    kmdz[kte+1] = rhoz[kte+1] * dfm[kte];

    for (int k = kts+1; k <= kte-1; k++) {
        khdz[k] = std::max(Real(khdz[k]), Real(0.5_rt * s_aw[k]));
        khdz[k] = std::max(Real(khdz[k]), Real(-0.5_rt * (s_aw[k] - s_aw[k+1])));
        kmdz[k] = std::max(Real(kmdz[k]), Real(0.5_rt * s_aw[k]));
        kmdz[k] = std::max(Real(kmdz[k]), Real(-0.5_rt * (s_aw[k] - s_aw[k+1])));
    }

    Real ustdrag = std::min(ust_v * ust_v, 0.99_rt) / wspd;
    Real ustdiff = std::min(ust_v * ust_v, 0.01_rt) / wspd;

    for (int k = kts; k <= kte; k++) {
        dth[k] = 0.0;
    }

    int k = kts;
    a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
    b[k] = 1.0_rt + dtz[k] * (kmdz[k+1] + rhosfc * ust_v * ust_v / wspd) * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
    c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
    d[k] = u[k] + dtz[k] * uoce * ust_v * ust_v / wspd - dtz[k] * rhoinv[k] * s_awu[k+1] * onoff + dtz[k] * rhoinv[k] * sd_awu[k+1] * onoff + sub_u[k] * delt + det_u[k] * delt;

    for (int k = kts+1; k <= kte-1; k++) {
        a[k] = -dtz[k] * kmdz[k] * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k] * onoff + 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k] * onoff;
        b[k] = 1.0_rt + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) * onoff + 0.5_rt * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]) * onoff;
        c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
        d[k] = u[k] + dtz[k] * rhoinv[k] * (s_awu[k] - s_awu[k+1]) * onoff - dtz[k] * rhoinv[k] * (sd_awu[k] - sd_awu[k+1]) * onoff + sub_u[k] * delt + det_u[k] * delt;
    }

    a[kte] = 0.0;
    b[kte] = 1.0;
    c[kte] = 0.0;
    d[kte] = u[kte];
    tridiag2_cc(kte, a, b, c, d, x);
    for (int k = kts; k <= kte; k++) {
        du[k] = (x[k] - u[k]) / delt;
    }

    k = kts;
    a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
    b[k] = 1.0_rt + dtz[k] * (kmdz[k+1] + rhosfc * ust_v * ust_v / wspd) * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
    c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
    d[k] = v[k] + dtz[k] * voce * ust_v * ust_v / wspd - dtz[k] * rhoinv[k] * s_awv[k+1] * onoff + dtz[k] * rhoinv[k] * sd_awv[k+1] * onoff + sub_v[k] * delt + det_v[k] * delt;

    for (int k = kts+1; k <= kte-1; k++) {
        a[k] = -dtz[k] * kmdz[k] * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k] * onoff + 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k] * onoff;
        b[k] = 1.0_rt + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) * onoff + 0.5_rt * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]) * onoff;
        c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
        d[k] = v[k] + dtz[k] * rhoinv[k] * (s_awv[k] - s_awv[k+1]) * onoff - dtz[k] * rhoinv[k] * (sd_awv[k] - sd_awv[k+1]) * onoff + sub_v[k] * delt + det_v[k] * delt;
    }

    a[kte] = 0.0;
    b[kte] = 1.0;
    c[kte] = 0.0;
    d[kte] = v[kte];
    tridiag2_cc(kte, a, b, c, d, x);
    for (int k = kts; k <= kte; k++) {
        dv[k] = (x[k] - v[k]) / delt;
    }

    k = kts;
    a[k] = -dtz[k] * khdz[k] * rhoinv[k];
    b[k] = 1.0_rt + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1];
    c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1];
    d[k] = thl[k] + dtz[k] * rhosfc * flt * rhoinv[k] + tcd[k] * delt - dtz[k] * rhoinv[k] * s_awthl[k+1] - dtz[k] * rhoinv[k] * sd_awthl[k+1] + diss_heat[k] * delt + sub_thl[k] * delt + det_thl[k] * delt;

    for (int k = kts+1; k <= kte-1; k++) {
        a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k] + 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k];
        b[k] = 1.0_rt + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) + 0.5_rt * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]);
        c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1];
        d[k] = thl[k] + tcd[k] * delt + dtz[k] * rhoinv[k] * (s_awthl[k] - s_awthl[k+1]) + dtz[k] * rhoinv[k] * (sd_awthl[k] - sd_awthl[k+1]) + diss_heat[k] * delt + sub_thl[k] * delt + det_thl[k] * delt;
    }

    a[kte] = 0.0;
    b[kte] = 1.0;
    c[kte] = 0.0;
    d[kte] = thl[kte];
    tridiag2_cc(kte, a, b, c, d, x);

    for (int k = kts; k <= kte; k++) {
        thl[k] = x[k];
    }

    if (bl_mynn_mixqt > 0) {
        k = kts;
        a[k] = -dtz[k] * khdz[k] * rhoinv[k];
        b[k] = 1.0_rt + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1];
        c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1];
        d[k] = sqw[k] + dtz[k] * rhosfc * flq * rhoinv[k] + qcd[k] * delt - dtz[k] * rhoinv[k] * s_awqt[k+1] - dtz[k] * rhoinv[k] * sd_awqt[k+1];

        for (int k = kts+1; k <= kte-1; k++) {
            a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k] + 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k];
            b[k] = 1.0_rt + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) + 0.5_rt * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]);
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1];
            d[k] = sqw[k] + qcd[k] * delt + dtz[k] * rhoinv[k] * (s_awqt[k] - s_awqt[k+1]) + dtz[k] * rhoinv[k] * (sd_awqt[k] - sd_awqt[k+1]);
        }

        a[kte] = 0.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = sqw[kte];
        tridiag2_cc(kte, a, b, c, d, sqw2);
    } else {
        for (int k = kts; k <= kte; k++) {
            sqw2[k] = sqw[k];
        }
    }

    if (bl_mynn_mixqt == 0) {
        if (bl_mynn_cloudmix > 0 && flag_qc==1) {
            k = kts;
            a[k] = -dtz[k] * khdz[k] * rhoinv[k];
            b[k] = 1.0_rt + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1];
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1];
            d[k] = sqc[k] + dtz[k] * rhosfc * flqc * rhoinv[k] + qcd[k] * delt - dtz[k] * rhoinv[k] * s_awqc[k+1] - dtz[k] * rhoinv[k] * sd_awqc[k+1] + det_sqc[k] * delt;

            for (int k = kts+1; k <= kte-1; k++) {
                a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k] + 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k];
                b[k] = 1.0_rt + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) + 0.5_rt * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]);
                c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1];
                d[k] = sqc[k] + qcd[k] * delt + dtz[k] * rhoinv[k] * (s_awqc[k] - s_awqc[k+1]) + dtz[k] * rhoinv[k] * (sd_awqc[k] - sd_awqc[k+1]) + det_sqc[k] * delt;
            }

            a[kte] = 0.0;
            b[kte] = 1.0;
            c[kte] = 0.0;
            d[kte] = sqc[kte];
            tridiag2_cc(kte, a, b, c, d, sqc2);
        } else {
            for (int k = kts; k <= kte; k++) {
                sqc2[k] = sqc[k];
            }
        }
    if (bl_mynn_mixqt == 0) {
        k = kts;
        qvflux = flqv;
        if (qvflux < 0.0_rt) {
            qvflux = std::max(Real(qvflux), (std::min(0.9_rt * Real(sqv[kts]) - 1e-8_rt, 0.0_rt) / dtz[kts]));
        }
        a[k] = -dtz[k] * khdz[k] * rhoinv[k];
        b[k] = 1.0_rt + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1];
        c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1];
        d[k] = sqv[k] + dtz[k] * rhosfc * qvflux * rhoinv[k] + qcd[k] * delt - dtz[k] * rhoinv[k] * s_awqv[k+1] - dtz[k] * rhoinv[k] * sd_awqv[k+1] + sub_sqv[k] * delt + det_sqv[k] * delt;

        for (int k = kts+1; k <= kte-1; k++) {
            a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k] + 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k];
            b[k] = 1.0_rt + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) + 0.5_rt * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]);
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5_rt * dtz[k] * rhoinv[k] * sd_aw[k+1];
            d[k] = sqv[k] + qcd[k] * delt + dtz[k] * rhoinv[k] * (s_awqv[k] - s_awqv[k+1]) + dtz[k] * rhoinv[k] * (sd_awqv[k] - sd_awqv[k+1]) + sub_sqv[k] * delt + det_sqv[k] * delt;
        }

        a[kte] = 0.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = sqv[kte];
        tridiag2_cc(kte, a, b, c, d, sqv2);
    } else {
      for (int k = kts; k <= kte; k++) {
            sqv2[k] = sqv[k];
        }
    }

        //Missing MIX CLOUD ICE
    if(bl_mynn_cloudmix > 0 && flag_qi) {
        k = kts;

        a[k] = -dtz[k] * khdz[k] * rhoinv[k];
        b[k] = 1.0_rt + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k];
        c[k] = -dtz[k] * khdz[k+1] * rhoinv[k];
        d[k] = sqi[k];

        for (int k = kts+1; k <= kte-1; k++) {
            a[k] = -dtz[k] * khdz[k] * rhoinv[k];
            b[k] = 1.0_rt + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k];
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k];
            d[k] = sqi[k];
        }

        a[kte] = 0.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = sqi[kte];
        tridiag2_cc(kte, a, b, c, d, sqi2);
    } else {
        for (int k = kts; k <= kte; k++) {
            sqi2[k] = sqi[k];
        }
    }

    //Missing MIX SNOW
    if(bl_mynn_cloudmix > 0 && false) {
        k = kts;

        a[k] = -dtz[k] * khdz[k] * rhoinv[k];
        b[k] = 1.0_rt + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k];
        c[k] = -dtz[k] * khdz[k+1] * rhoinv[k];
        d[k] = sqs[k];

        for (int k = kts+1; k <= kte-1; k++) {
            a[k] = -dtz[k] * khdz[k] * rhoinv[k];
            b[k] = 1.0_rt + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k];
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k];
            d[k] = sqs[k];
        }

        a[kte] = 0.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = sqs[kte];
        tridiag2_cc(kte, a, b, c, d, sqs2);
    } else {
        for (int k = kts; k <= kte; k++) {
            sqs2[k] = sqs[k];
        }
    }
    //Missing ice number concentration (qni)
        if (bl_mynn_cloudmix > 0 && flag_qni==1 && bl_mynn_mixscalars > 0) {
            k = kts;
            a[k] = -dtz[k] * khdz[k] * rhoinv[k];
            b[k] = 1.0_rt + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] * nonloc;
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] * nonloc;
            d[k] = qni[k] - dtz[k] * rhoinv[k] * s_awqni[k+1]*nonloc;

            for (int k = kts+1; k <= kte-1; k++) {
                a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k];
                b[k] = 1.0_rt + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) * nonloc;
                c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] * nonloc;
                d[k] = qni[k] + dtz[k] * rhoinv[k] * (s_awqni[k] - s_awqni[k+1]) * nonloc;
            }

            a[kte] = 0.0;
            b[kte] = 1.0;
            c[kte] = 0.0;
            d[kte] = qni[kte];
            tridiag2_cc(kte, a, b, c, d, x);
            for (int k = kts; k <= kte; k++) {
                qni2[k] = x[k];
            }
        } else {
            for (int k = kts; k <= kte; k++) {
                qni2[k] = x[k];
            }
        }
    //Missing cloud number concentration (qnc)
       if (bl_mynn_cloudmix > 0 && flag_qnc==1 && bl_mynn_mixscalars > 0) {
            k = kts;
            a[k] = -dtz[k] * khdz[k] * rhoinv[k];
            b[k] = 1.0_rt + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] * nonloc;
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] * nonloc;
            d[k] = qni[k] - dtz[k] * rhoinv[k] * s_awqnc[k+1]*nonloc;

            for (int k = kts+1; k <= kte-1; k++) {
                a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k];
                b[k] = 1.0_rt + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5_rt * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) * nonloc;
                c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5_rt * dtz[k] * rhoinv[k] * s_aw[k+1] * nonloc;
                d[k] = qnc[k] + dtz[k] * rhoinv[k] * (s_awqnc[k] - s_awqnc[k+1]) * nonloc;
            }

            a[kte] = 0.0;
            b[kte] = 1.0;
            c[kte] = 0.0;
            d[kte] = qnc[kte];
            tridiag2_cc(kte, a, b, c, d, x);
            for (int k = kts; k <= kte; k++) {
                qnc2[k] = x[k];
            }
        } else {
            for (int k = kts; k <= kte; k++) {
                qnc2[k] = x[k];
            }
        }
    //Missing water-friendly aerosols (qnwfa)
    //Missing Ice-friendly aerosols (qnifa)
    //Missing Black-carbon aerosols (qnbca)
    //Missing Ozone
    }
    //Missing compute tendencies and convert to mixing ratios
    for (int k = kts; k <= kte; k++) {
        dqv[k] = (sqv2[k] - sqv[k]) / delt;
    }

    if(bl_mynn_cloudmix > 0) {
      if(flag_qc) {
        for (int k = kts; k <= kte; k++) {
          dqc[k] = (sqc2[k] - sqc[k]) / delt;
        }
      } else {
        for (int k = kts; k <= kte; k++) {
          dqc[k] = 0.0;
        }
      }

      if(flag_qnc && bl_mynn_mixscalars > 0) {
        for (int k = kts; k <= kte; k++) {
          dqnc[k] = (qnc2[k] - qnc[k]) / delt;
        }
      } else {
        for (int k = kts; k <= kte; k++) {
          dqnc[k] = 0.0;
        }
      }

      if(flag_qi) {
        for (int k = kts; k <= kte; k++) {
          dqi[k] = (sqi2[k] - sqi[k]) / delt;
        }
      } else {
        for (int k = kts; k <= kte; k++) {
          dqi[k] = 0.0;
        }
      }

      if(false) {
        for (int k = kts; k <= kte; k++) {
          dqs[k] = (sqs2[k] - sqs[k]) / delt;
        }
      } else {
        for (int k = kts; k <= kte; k++) {
          dqs[k] = 0.0;
        }
      }

      if(flag_qni && bl_mynn_mixscalars > 0) {
        for (int k = kts; k <= kte; k++) {
          dqni[k] = (qni2[k] - qni[k]) / delt;
        }
      } else {
        for (int k = kts; k <= kte; k++) {
          dqni[k] = 0.0;
        }
      }
    } else {
        for (int k = kts; k <= kte; k++) {
          dqc[k] = 0.0;
          dqnc[k] = 0.0;
          dqi[k] = 0.0;
          dqni[k] = 0.0;
          dqs[k] = 0.0;
        }
    }
    if(kts != 0) {
      printf("moisture check assumes kts for cpp is 0");
      exit(1);
    }

    moisture_check_cc(kte, delt, delp, exner,
                        sqv2, sqc2, sqi2, sqs2, thl,
                      dqv, dqc, dqi, dqs, dth,xlvcp,xlscp);
    //skip ozone for now

    for (int k = kts; k <= kte; k++) {
      dozone[k]=0.0;
      if(dozone[k]*delt+ozone[k]<0.0_rt)
        dozone[k]=-ozone[k]*0.99_rt/delt;
        }
    if(flag_qi)
      for (int k = kts; k <= kte; k++) {
        dth[k]=(thl[k] + xlvcp/exner[k]*sqc2[k]
                + xlscp/exner[k]*(sqi2[k])        //+sqs(k))
                - th[k])/delt;
      }
    else
      for (int k = kts; k <= kte; k++) {
        dth[k]=(thl[k]+xlvcp/exner[k]*sqc2[k] - th[k])/delt;
      }
    //Skip qnwfa qnifa qnbca tendencies for now
    /*
    if(flag_qnwfa && flag_qnifa && bl_mynn_mixscalars > 0) {
      for (int k = kts; k <= kte; k++) {
        dqnwfa[k]=(qnwfa2[k] - qnwfa[k])/delt;
        dqnifa[k]=(qnifa2[k] - qnifa[k])/delt;
      }
    } else {
      for (int k = kts; k <= kte; k++) {
        dqnwfa[k]=0.0_rt;
        dqnifa[k]=0.0_rt;
      }
    }

    if(flag_qnbca && bl_mynn_mixscalars > 0) {
      for (int k = kts; k <= kte; k++) {
        dqnbca[k]=(qnbca2[k] - qnbca[k])/delt;
      }
    } else {
      for (int k = kts; k <= kte; k++) {
        dqnbca[k]=0.0_rt;
      }
    }
    */
    /*
    printf("thl ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",thl[k]);
    printf("\n");

   printf("sqv ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",sqv[k]);
    printf("\n");

    printf("sqc ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",sqc[k]);
    printf("\n");

    printf("sqi ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",sqi[k]);
    printf("\n");

    printf("sqs ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",sqs[k]);
    printf("\n");

    printf("sqw ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",sqw[k]);
    printf("\n");

    printf("qnwfa ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",qnwfa[k]);
    printf("\n");

    printf("qnifa ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",qnifa[k]);
    printf("\n");

    printf("qnbca ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",qnbca[k]);
    printf("\n");

    printf("ozone ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",ozone[k]);
    printf("\n");

    printf("dfm ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",dfm[k]);
    printf("\n");

    printf("dfh ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",dfh[k]);
    printf("\n");

    printf("du ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",du[k]);
    printf("\n");

    printf("dv ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",dv[k]);
    printf("\n");

    printf("dth ");
    for (int k = kts; k <= kte; k++)
      printf("%g ",dth[k]);
    printf("\n");

    //exit(1);
    */
    //    /*inout*/ Real* thl, Real* sqv, Real* sqc, Real* sqi,
    //    /*inout*/ Real* sqs, Real* sqw, Real* qnwfa, Real* qnifa, Real* qnbca, Real* ozone,
    //    /*inout*/ Real* dfm, Real* dfh,
    //    /*inout*/ Real* du, Real* dv, Real* dth,

}



void mym_condensation_cc(
    const int& kts,const int& kte,
    const Real& dx, Real* dz, Real* zw, Real& xland,
    Real* thl, Real* qw, Real* qv, Real* qc, Real* qi, Real* qs,
    Real* p, Real* exner, Real* tsq, Real* qsq, Real* cov,
    Real* sh, Real* el, int& bl_mynn_cloudpdf,
    Real* qc_bl1d, Real* qi_bl1d, Real* cldfra_bl1d,
    Real& pblh1, Real& hfx1,
    Real* vt, Real* vq, Real* th, Real* sgm, Real* rmo,
    int &spp_pbl, Real* rstoch_col,
    // model constants
    Real ep_2, Real ep_3, Real xlv, Real r_d, Real xlvcp,
    Real p608, Real tv0, Real r_v, Real cice, Real cliq,
    Real cpv, Real cp, Real xls, Real rcp)
{
    int k;
    Real t3sq, r3sq, c3sq;
    Real qsl, esat, qsat, dqsl, cld0, q1k, qlk, eq1, qll, q2p, pt, rac, qt, t, xl, rsl, cpm, fng, qww, alpha, beta, bb, ls, wt, wt2, qpct, cld_factor, fac_damp, liq_frac, ql_ice, ql_water, qmq, qsat_tk, q1_rh, rh_hack, dzm1, zsl, maxqc;
    const Real qpct_sfc = 0.025;
    const Real qpct_pbl = 0.030;
    const Real qpct_trp = 0.040;
    const Real rhcrit = 0.83;
    const Real rhmax = 1.02;
    Real erf;
    Real dth, dtl, dqw, dzk, els;
    Real zagl, damp, pblh2;
    Real cfmax;
    Real theta1, theta2, ht1, ht2;
    Real qw_pert;
    int k_tropo;

//real(float), dimension(kts:kte) :: alp,a,bet,b,ql,q1,rh
    Real alp[kte-kts];
    Real a[kte-kts];
    Real bet[kte-kts];
    Real b[kte-kts];
    Real ql[kte-kts];
    Real q1[kte-kts];
    Real rh[kte-kts];

    // obtain an estimate for the tropopause height (k)
    for (k = kte - 3; k >= kts; k--) {
        theta1 = th[k];
        theta2 = th[k + 2];
        ht1 = 44307.692_rt * (1.0_rt - pow(p[k] / 101325.0, 0.190));
        ht2 = 44307.692_rt * (1.0_rt - pow(p[k + 2] / 101325.0, 0.190));
        if ((((theta2 - theta1) / (ht2 - ht1)) < 10.0_rt / 1500.0_rt) && (ht1 < 19000.0_rt) && (ht1 > 4000.0_rt)) {
            break;
        }
    }
    k_tropo = std::max(kts + 2, k + 2);
    zagl = 0.0;

    switch (bl_mynn_cloudpdf) {
        case 0: // original mynn partial-condensation scheme
            for (k = kts; k < kte; k++) {
                t = th[k] * exner[k];
                esat = esat_blend_cc(t);
                qsl = ep_2 * esat / std::max(1e-4_rt, (p[k] - ep_3 * esat));
                dqsl = qsl * ep_2 * xlv / (r_d * pow(t, 2));
                alp[k] = 1.0_rt / (1.0_rt + dqsl * xlvcp);
                bet[k] = dqsl * exner[k];
                t3sq = std::max(tsq[k], 0.0_rt);
                r3sq = std::max(qsq[k], 0.0_rt);
                c3sq = cov[k];
                c3sq = std::copysign(std::min(std::abs(c3sq), std::sqrt(t3sq * r3sq)), c3sq);
                r3sq = r3sq + bet[k] * bet[k] * t3sq - 2.0_rt * bet[k] * c3sq;
                qmq = qw[k] - qsl;
                sgm[k] = std::sqrt(std::max(r3sq, 1.0e-10_rt));
                q1[k] = qmq / sgm[k];
                cldfra_bl1d[k] = 0.5_rt * (1.0_rt + std::erf(q1[k] * rr2));
                q1k = q1[k];
                eq1 = rrp * std::exp(-0.5_rt * q1k * q1k);
                qll = std::max(cldfra_bl1d[k] * q1k + eq1, 0.0_rt);
                ql[k] = alp[k] * sgm[k] * qll;
                liq_frac = std::min(1.0_rt, std::max(0.0_rt, (t - 240.0_rt) / 29.0_rt));
                qc_bl1d[k] = liq_frac * ql[k];
                qi_bl1d[k] = (1.0_rt - liq_frac) * ql[k];
                q2p = xlvcp / exner[k];
                pt = thl[k] + q2p * ql[k];
                qt = 1.0_rt + p608 * qw[k] - (1.0_rt + p608) * (qc_bl1d[k] + qi_bl1d[k]) * cldfra_bl1d[k];
                rac = alp[k] * (cldfra_bl1d[k] - qll * eq1) * (q2p * qt - (1.0_rt + p608) * pt);
                vt[k] = qt - 1.0_rt - rac * bet[k];
                vq[k] = p608 * pt - tv0 + rac;
                //    printf("vt[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",vt[k],qt,rac,bet[k],dqsl);
    //    printf("vq[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",vq[k],p608,pt,tv0,rac);
            }
            break;
        case 1:
        case -1: // alternative form (nakanishi & niino 2004 blm, eq. b6, and kuwano-yoshida et al. 2010 qjrms, eq. 7)
            for (k = kts; k < kte; k++) {
                t = th[k] * exner[k];
                esat = esat_blend_cc(t);
                qsl = ep_2 * esat / std::max(1e-4_rt, (p[k] - ep_3 * esat));
                dqsl = qsl * ep_2 * xlv / (r_d * pow(t, 2));
                alp[k] = 1.0_rt / (1.0_rt + dqsl * xlvcp);
                bet[k] = dqsl * exner[k];
                if (k == kts) {
                    dzk = 0.5_rt * dz[k];
                } else {
                    dzk = dz[k];
                }
                dth = 0.5_rt * (thl[k + 1] + thl[k]) - 0.5_rt * (thl[k] + thl[std::max(k - 1, kts)]);
                dqw = 0.5_rt * (qw[k + 1] + qw[k]) - 0.5_rt * (qw[k] + qw[std::max(k - 1, kts)]);
                sgm[k] = std::sqrt(std::max(Real((pow(alp[k], 2) * std::max(pow(el[k], 2), 0.1) * b2 * std::max(sh[k], 0.03_rt)) / 4.0_rt * pow((dqw / dzk - bet[k] * (dth / dzk)), 2)), 1.0e-10_rt));
                qmq = qw[k] - qsl;
                q1[k] = qmq / sgm[k];
                cldfra_bl1d[k] = 0.5_rt * (1.0_rt + std::erf(q1[k] * rr2));
                q1k = q1[k];
                eq1 = rrp * std::exp(-0.5_rt * q1k * q1k);
                qll = std::max(cldfra_bl1d[k] * q1k + eq1, 0.0_rt);
                ql[k] = alp[k] * sgm[k] * qll;
                liq_frac = std::min(1.0_rt, std::max(0.0_rt, (t - 240.0_rt) / 29.0_rt));
                qc_bl1d[k] = liq_frac * ql[k];
                qi_bl1d[k] = (1.0_rt - liq_frac) * ql[k];
                q2p = xlvcp / exner[k];
                pt = thl[k] + q2p * ql[k];
                qt = 1.0_rt + p608 * qw[k] - (1.0_rt + p608) * (qc_bl1d[k] + qi_bl1d[k]) * cldfra_bl1d[k];
                rac = alp[k] * (cldfra_bl1d[k] - qll * eq1) * (q2p * qt - (1.0_rt + p608) * pt);
                vt[k] = qt - 1.0_rt - rac * bet[k];
                vq[k] = p608 * pt - tv0 + rac;
                //    printf("vt[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",vt[k],qt,rac,bet[k],dqsl);
    //    printf("vq[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",vq[k],p608,pt,tv0,rac);
            }
            break;
        case 2:
        case -2: // diagnostic statistical scheme of Chaboureau and Bechtold (2002), JAS
                 // but with use of higher-order moments to estimate sigma
            pblh2 = std::max(10.0_rt, pblh1);
            zagl = 0.0_rt;
            dzm1 = 0.0_rt;
            for (k = kts; k < kte; k++) {
                zagl += 0.5_rt * (dz[k] + dzm1);
                dzm1 = dz[k];

                t = th[k] * exner[k];
                xl = xl_blend_cc(t,xlv,xls,cpv,cliq,cice);
                qsat_tk = qsat_blend_cc(t, p[k]);
                rh[k] = std::max(std::min(rhmax, qw[k] / std::max(1e-10_rt, qsat_tk)), 0.001_rt);

                // dqw/dT: Clasius-Clapeyron
                dqsl = qsat_tk * ep_2 * xlv / (r_d * std::pow(t, 2.0_rt));
                alp[k] = 1.0_rt / (1.0_rt + dqsl * xlvcp);
                bet[k] = dqsl * exner[k];

                rsl = xl * qsat_tk / (r_v * std::pow(t, 2.0_rt)); // slope of C-C curve at (= std::abs temperature) (CB02, Eqn. 4)
                cpm = cp + qw[k] * cpv;                 // CB02, sec. 2, para. 1
                a[k] = 1.0_rt / (1.0_rt + xl * rsl / cpm);  // CB02 variable "a"
                b[k] = a[k] * rsl;                      // CB02 variable "b"

                // SPP
                qw_pert = qw[k] + qw[k] * 0.5_rt * rstoch_col[k] * spp_pbl;

                // This form of qmq (the numerator of Q1) no longer uses the a[k] factor
                qmq = qw_pert - qsat_tk;

                // Use the form of Eq. 6 in CB02 except neglect all but the first term for sig_r
                r3sq = std::max(qsq[k], 0.0_rt);
                // Calculate sigma using higher-order moments
                sgm[k] = std::sqrt(r3sq);
                // Set constraints on sigma relative saturation water vapor
                sgm[k] = std::min(sgm[k], qsat_tk * 0.666_rt);

                // Introduce vertical grid spacing dependence on min sgm
                wt = std::max(500.0_rt - std::max(dz[k] - 100.0_rt, 0.0_rt), 0.0_rt) / 500.0_rt;
                sgm[k] += sgm[k] * 0.2_rt * (1.0_rt - wt); // inflate sgm for coarse dz

                // Allow min sgm to vary with dz and z
                qpct = qpct_pbl * wt + qpct_trp * (1.0_rt - wt);
                qpct = std::min(qpct, std::max(qpct_sfc, qpct_pbl * zagl / 500.0_rt));
                sgm[k] = std::max(sgm[k], qsat_tk * qpct);

                q1[k] = qmq / sgm[k]; // normalized saturation

                // Add condition for falling/settling into low-RH layers, so at least
                // some cloud fraction is applied for all qc, qs, and qi.
                rh_hack= rh[k];
                wt2    = std::min(std::max( zagl - pblh2, 0.0_rt )/300.0_rt, 1.0_rt);
                // ensure adequate RH & q1 when qi is at least 1e-9 (above the PBLH)
                if ((qi[k]+qs[k])>1.e-9 && (zagl > pblh2)) {
                  rh_hack =std::min(rhmax, rhcrit + wt2*0.045_rt*(9.0_rt + std::log10(qi[k]+qs[k])));
                  rh[k]   =std::max(rh[k], rh_hack);

                  q1_rh   =-3.0_rt + 3.0_rt*(rh[k]-rhcrit)/(1.0_rt-rhcrit);
                  q1[k]   =std::max(q1_rh, q1[k] );
                }
                // ensure adequate RH & q1 when qi is at least 1e-6 (above the PBLH)
                if (qc[k]>1.e-6 && (zagl > pblh2)) {
                  rh_hack =std::min(rhmax, rhcrit + wt2*0.08_rt*(6.0_rt + std::log10(qc[k])));
                  rh[k]   =std::max(rh[k], rh_hack);

                  q1_rh   =-3.0_rt + 3.0_rt*(rh[k]-rhcrit)/(1.0_rt-rhcrit);
                  q1[k]   =std::max(q1_rh, q1[k] );
                }

                q1k = q1[k]; // backup q1 for later modification

                // Specify cloud fraction
                //Original C-B cloud fraction, allows cloud fractions out to q1 = -3.5
                //cldfra_bl1D(K) = max(0., min(1., 0.5+0.36*atan(1.55*q1(k)))) ! Eq. 7 in CB02
                //Waynes LES fit  - over-diffuse, when limits removed from vt & vq & fng
                //cldfra_bl1D(K) = max(0., min(1., 0.5+0.36*atan(1.2*(q1(k)+0.4))))
                // Best compromise: Improves marine stratus without adding much cold bias.
                cldfra_bl1d[k] = std::max(0.0_rt, std::min(1.0_rt, 0.5_rt + 0.36_rt * std::atan(1.8_rt * (q1[k] + 0.2_rt))));

                // Specify hydrometeors
                // The cloud water formulatiosn are taken from CB02, Eq. 8
                maxqc = std::max(qw[k] - qsat_tk, 0.0_rt);
                if (q1k < 0.0_rt) {
                    ql_water = sgm[k] * std::exp(1.2_rt * q1k - 1.0_rt);
                    ql_ice = sgm[k] * std::exp(1.2_rt * q1k - 1.0_rt);
                } else if (q1k > 2.0_rt) {
                    ql_water = std::min(sgm[k] * q1k, maxqc);
                    ql_ice = sgm[k] * q1k;
                } else {
                    ql_water = std::min(Real(sgm[k] * (std::exp(-1.0_rt) + 0.66_rt * q1k + 0.086_rt * pow(q1k, 2.0_rt))), maxqc);
                    ql_ice = sgm[k] * (std::exp(-1.0_rt) + 0.66_rt * q1k + 0.086_rt * pow(q1k, 2));
                }

                // In saturated grid cells, use average of SGS and resolved values
                // if ( qc[k] > 1.e-6 ) ql_water = 0.5 * ( ql_water + qc[k] )
                // ql_ice is actually the total frozen condensate (snow+ice),
                // if ( (qi[k]+qs[k]) > 1.e-9 ) ql_ice = 0.5 * ( ql_ice + (qi[k]+qs[k]) )

                if (cldfra_bl1d[k] < 0.001_rt) {
                    ql_ice = 0.0_rt;
                    ql_water = 0.0_rt;
                    cldfra_bl1d[k] = 0.0_rt;
                }

                liq_frac = std::min(1.0_rt, std::max(0.0_rt, (t - tice) / (tliq - tice)));
                qc_bl1d[k] = liq_frac * ql_water;
                qi_bl1d[k] = (1.0_rt - liq_frac) * ql_ice;

                // Above tropopausel: eliminate subgrid clouds from CB schemes.
                if (k >= k_tropo) {
                    cldfra_bl1d[k] = 0.0_rt;
                    qc_bl1d[k] = 0.0_rt;
                    qi_bl1d[k] = 0.0_rt;
                }

                // Buoyancy-flux-related calculations follow...

                if((xland-1.5_rt)>=0.0_rt) {
                    q1k = std::max(q1[k],-2.5_rt); // water
                } else {
                    q1k = std::max(q1[k],-2.0_rt); // land
                }

                // "Fng" represents the non-Gaussian transport factor
                // (non-dimensional) from Bechtold et al. 1995
                // (hereafter BCMT95), section 3(c).
                //
                // Use the form of "Fng" from Bechtold and Siebesma (1998, JAS)
                if(q1k >= 1.0_rt) {
                    fng = 1.0;
                } else if(q1k >= -1.7_rt && q1k < 1.0_rt) {
                    fng = exp(-0.4_rt*(q1k-1.0_rt));
                } else if(q1k >= -2.5_rt && q1k < -1.7_rt) {
                    fng = 3.0_rt + Real(exp(-3.8_rt*(q1k+1.7_rt)));
                } else {
                    fng = std::min(23.9_rt + Real(exp(-1.6_rt*(q1k+2.5_rt))),60._rt);
                }

                cfmax = std::min(cldfra_bl1d[k],0.6_rt);
                // Further limit the cf going into vt & vq near the surface
                zsl = std::min(std::max(25._rt, 0.1_rt*pblh2), 100._rt);
                wt  = std::min(zagl/zsl, 1.0_rt); // ==0 at z=0 m, ==1 above ekman layer
                cfmax = cfmax*wt;

                // bb is "b" in BCMT95. Their "b" differs from "b" in CB02
                // (i.e., b[k] above) by a factor of T/theta. Strictly, b[k]
                // above is formulated in terms of sat. mixing ratio, but bb in
                // BCMT95 is cast in terms of sat. specific humidity. The
                // conversion is neglected here.
                bb = b[k]*t/th[k];
                qww = 1._rt+0.61_rt*qw[k];
                alpha = 0.61_rt*th[k];
                beta = (th[k]/t)*(xl/cp) - 1.61_rt*th[k];
                vt[k] = qww - cfmax*beta*bb*fng - 1._rt;
                vq[k] = alpha + cfmax*beta*a[k]*fng - tv0;

                // Dampen amplification factor where need be
                fac_damp = std::min(zagl * 0.0025_rt, 1.0_rt);
                cld_factor = 1.0_rt + fac_damp * std::min(std::pow(std::max(0.0_rt, (rh[k] - 0.92_rt)) / 0.145_rt, 2.0_rt), 0.37_rt);
                cldfra_bl1d[k] = std::min(1.0_rt, cld_factor * cldfra_bl1d[k]);
                // EQDEBUG
                //              printf("vt[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",vt[k],cfmax,beta,bb,fng);
                //              printf("vq[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",vq[k],alpha,cfmax,beta,a[k]);
            }
            break;

    } // switch (bl_mynn_cloudpdf)

    if (bl_mynn_cloudpdf < 0) {
      for( k = kts;k<=kte-1;k++) {
        cldfra_bl1d[k] = 0.0;
        qc_bl1d[k] = 0.0;
        qi_bl1d[k] = 0.0;
      }
    }
    ql[kte] = ql[kte - 1];
    vt[kte] = vt[kte - 1];
    vq[kte] = vq[kte - 1];
    qc_bl1d[kte] = 0.0;
    qi_bl1d[kte] = 0.0;
    cldfra_bl1d[kte] = 0.0;
}


//===============================================================
// ===================================================================
// this is the downdraft mass flux scheme - analogous to edmf_jpl but
// flipped updraft to downdraft. this scheme is currently only tested
// for stratocumulus cloud conditions. for a detailed description of the
// model, see paper.
void ddmf_jpl_cc(int& kts, int& kte, Real& dt, const Real* zw, const Real* dz, const Real* p,
              const Real* u, const Real* v, const Real* th, const Real* thl, const Real* thv,
              const Real* tk,const Real* qt, const Real* qv, const Real* qc, const Real*
              rho, const Real* exner,Real& ust, Real& wthl, Real& wqt, Real& pblh, int& kpbl,
              Real* edmf_a_dd, Real* edmf_w_dd, Real* edmf_qt_dd,
              Real* edmf_thl_dd, Real* edmf_ent_dd, Real* edmf_qc_dd,
              Real* sd_aw, Real* sd_awthl, Real* sd_awqt,
              Real* sd_awqv, Real* sd_awqc, Real* sd_awu,
              Real* sd_awv, Real* sd_awqke,
              const Real* qc_bl1d, const Real* cldfra_bl1d,
              const Real* rthraten, Real& svp1, Real& grav, Real& onethird, Real& p1000mb,
              Real& rcp, Real& xlvcp, Real& cp, Real& rvovrd ) {
    int ndown = 5;
    int dd_initk[ndown];
    Real randnum[ndown];
    Real downw[kte + 1][ndown];
    Real downthl[kte + 1][ndown];
    Real downqt[kte + 1][ndown];
    Real downqc[kte + 1][ndown];
    Real downa[kte + 1][ndown];
    Real downu[kte + 1][ndown];
    Real downv[kte + 1][ndown];
    Real downthv[kte + 1][ndown];
    Real ent[kte + 1][ndown];
    int enti[kte + 1][ndown];
    int k, i, ki, kminrad, qltop, p700_ind, qlbase;
    Real wthv, wstar, qstar, thstar, sigmaw, sigmaqt, sigmath, z0, pwmin, pwmax, wmin, wmax, wlv, wtv, went, mindownw;
    Real b, qtn, thln, thvn, qcn, un, vn, qken, wn2, wn, thvk, pk, entexp, entw, beta_dm, entexp_m, rho_int;
    Real jump_thetav, jump_qt, jump_thetal, refthl, refthv, refqt;
    Real minrad, zminrad, radflux, f0, wst_rad, wst_dd;
    bool cloudflg;
    Real sigq, xl, rsl, cpm, a, mf_cf, diffqt, fng, qww, alpha, beta, bb, f, pt, t, q2p, b9, satvp, rhgrid;
    Real wa = 1.0, wb = 1.5, z00 = 100.0, bcoeff = 0.2;
    Real l0 = 80, ent0 = 0.2;
    Real dp, dl, adn;
    int debug_mf = 0;
    dl = (1000.0_rt - 500.0_rt) / ndown;
    pwmin = -3.0;
    pwmax = -1.0;
    for(k=kts;k<=kte+1;k++) {
    for(i=0;i<ndown;i++) {
    downw[k][i] = 0.0_rt;
    downthl[k][i] = 0.0_rt;
    downthv[k][i] = 0.0_rt;
    downqt[k][i] = 0.0_rt;
    downqc[k][i] = 0.0_rt;
    downa[k][i] = 0.0_rt;
    downu[k][i] = 0.0_rt;
    downv[k][i] = 0.0_rt;
    }
    sd_aw[k] = 0.0_rt;
    sd_awthl[k] = 0.0_rt;
    sd_awqt[k] = 0.0_rt;
    sd_awqv[k] = 0.0_rt;
    sd_awqc[k] = 0.0_rt;
    sd_awu[k] = 0.0_rt;
    sd_awv[k] = 0.0_rt;
    sd_awqke[k] = 0.0_rt;
    }
    for(k=kts;k<=kte;k++) {
    edmf_a_dd[k] = 0.0_rt;
    edmf_w_dd[k] = 0.0_rt;
    edmf_qt_dd[k] = 0.0_rt;
    edmf_thl_dd[k] = 0.0_rt;
    edmf_ent_dd[k] = 0.0_rt;
    edmf_qc_dd[k] = 0.0_rt;
    }
    //from kts+1 to kte+1
    for(k=kts;k<=kte;k++) {
    for(i=0;i<ndown;i++) {
        ent[k][i] = 0.0_rt;
    }
    }
    for (int i = 0; i < ndown; i++) {
        dd_initk[i] = 0.0;
    }
    cloudflg = false;
    minrad = 100.0;
    kminrad = kpbl;
    zminrad = pblh;
    qltop = 0;
    qlbase = 0;
    wthv = wthl + svp1 * wqt;
    for (int i = std::max(2, kpbl - 2); i <= kpbl + 3; i++) {
        if (qc[i] > 1.0e-6 && cldfra_bl1d[i] > 0.5_rt) {
            cloudflg = true;
            qltop = i;
        }
    }
    for (int i = qltop; i >= kts; i--) {
        if (qc[i] > 1e-6) {
            qlbase = i;
        }
    }
    qlbase = (qltop + qlbase) / 2;
    for (int i = 0; i < ndown; i++) {
        dd_initk[i] = qltop;
    }

    f0 = 0.0;
    for (int i = 0; i <= qltop; i++) {
        radflux = rthraten[i] * exner[i];
        radflux = radflux * cp / grav * (p[i] - p[i + 1]);
        if (radflux < 0.0_rt) {
            f0 = std::abs(radflux) + f0;
        }
    }
    f0 = std::max(f0, 1.0_rt);
    adn = std::min(0.05_rt + f0 * 0.001_rt, 0.3_rt);
    if (cloudflg) {

    p700_ind = 0;
    Real min_value = p[0];
    for (int i = kts; i <= kte; ++i) {
        Real pval=std::abs(p[i]-70000.0_rt);
        if (pval < min_value) {
            p700_ind = i;
        }
    }

        //p700_ind = minloc(std::abs(p - 70000.0_rt), 1.0_rt);
        jump_thetav = thv[p700_ind] - thv[1] - (thv[p700_ind] - thv[qltop + 3]) / (zw[p700_ind] - zw[qltop + 3]) * (zw[p700_ind] - zw[qltop]);
        jump_qt = qc[p700_ind] + qv[p700_ind] - qc[1] - qv[1];
        jump_thetal = thl[p700_ind] - thl[1] - (thl[p700_ind] - thl[qltop + 3]) / (zw[p700_ind] - zw[qltop + 3]) * (zw[p700_ind] - zw[qltop]);
        refthl = thl[qltop];
        refthv = thv[qltop];
        refqt = qt[qltop];
        wst_rad = pow(grav * zw[qltop] * f0 / (refthl * rho[qltop] * cp), 0.333);
        wst_rad = std::max(wst_rad, 0.1_rt);
        wstar = std::max(0.0, pow(grav / thv[1] * wthv * pblh, onethird));
        went = thv[1] / (grav * jump_thetav * zw[qltop]) * (0.15_rt * (pow(wstar, 3) + 5 * pow(ust, 3)) + 0.35_rt * pow(wst_rad, 3));
        qstar = std::abs(went * jump_qt / wst_rad);
        thstar = f0 / (rho[qltop] * cp * wst_rad) - went * jump_thetav / wst_rad;
        wst_dd = pow(0.15_rt * (pow(wstar, 3) + 5 * pow(ust, 3)) + 0.35_rt * pow(wst_rad, 3), 0.333);

        sigmaw = 0.2_rt * wst_dd;
        sigmaqt = 40 * qstar;
        sigmath = 1.0_rt * thstar;
        wmin = sigmaw * pwmin;
        wmax = sigmaw * pwmax;
        for (int i = 0; i < ndown; i++) {
            ki = dd_initk[i];
            wlv = wmin + (wmax - wmin) / ndown * (i - 1);
            wtv = wmin + (wmax - wmin) / ndown * i;
            downw[ki][i] = wlv;
            downa[ki][i] = adn / ndown;
            downu[ki][i] = (u[ki - 1] * dz[ki] + u[ki] * dz[ki - 1]) / (dz[ki] + dz[ki - 1]);
            downv[ki][i] = (v[ki - 1] * dz[ki] + v[ki] * dz[ki - 1]) / (dz[ki] + dz[ki - 1]);
            refthl = (thl[ki - 1] * dz[ki] + thl[ki] * dz[ki - 1]) / (dz[ki] + dz[ki - 1]);
            refthv = (thv[ki - 1] * dz[ki] + thv[ki] * dz[ki - 1]) / (dz[ki] + dz[ki - 1]);
            refqt = (qt[ki - 1] * dz[ki] + qt[ki] * dz[ki - 1]) / (dz[ki] + dz[ki - 1]);
            downqc[ki][i] = (qc[ki - 1] * dz[ki] + qc[ki] * dz[ki - 1]) / (dz[ki] + dz[ki - 1]);
            downqt[ki][i] = refqt;
            downthv[ki][i] = refthv + 0.01_rt * downw[ki][i] * sigmath / sigmaw;
            downthl[ki][i] = refthl + 0.01_rt * downw[ki][i] * sigmath / sigmaw;

        for (int k = dd_initk[i] - 1; k >= kts + 1; k--) {
            wmin = 0.3_rt + dp * 0.0005;
            ent[k+1][i] = 0.33_rt / (std::min(std::max(-1.0_rt * downw[k + 1][i], wmin), 0.9_rt) * dp);
            entexp = ent[k+1][i] * dz[k];
            entexp_m = ent[k+1][i] * 0.333_rt * dz[k];
            qtn = downqt[k + 1][i] * (1.0_rt - entexp) + qt[k] * entexp;
            thln = downthl[k + 1][i] * (1.0_rt - entexp) + thl[k] * entexp;
            un = downu[k + 1][i] * (1.0_rt - entexp) + u[k] * entexp_m;
            vn = downv[k + 1][i] * (1.0_rt - entexp) + v[k] * entexp_m;
            pk = (p[k - 1] * dz[k] + p[k] * dz[k - 1]) / (dz[k] + dz[k - 1]);
            condensation_edmf_cc(qtn, thln, pk, zw[k], thvn, qcn,p1000mb,rcp,xlvcp,rvovrd);
            thvk = (thv[k - 1] * dz[k] + thv[k] * dz[k - 1]) / (dz[k] + dz[k - 1]);
            b = grav * (thvn / thvk - 1.0_rt);
            entw = entexp;
            mindownw = std::min(downw[k + 1][i], -0.2_rt);
            wn = downw[k + 1][i] + (-2.0_rt * ent[k+1][i] * downw[k + 1][i] - bcoeff * b / mindownw) * std::min(dz[k], 250.0_rt);
            if (wn < downw[k + 1][i] - std::min(1.25_rt * dz[k] / 200.0_rt, -2.0_rt)) {
                wn = downw[k + 1][i] - std::min(1.25_rt * dz[k] / 200.0_rt, -2.0_rt);
            }
            if (wn > downw[k + 1][i] + std::min(1.25_rt * dz[k] / 200.0_rt, 2.0_rt)) {
                wn = downw[k + 1][i] + std::min(1.25_rt * dz[k] / 200.0_rt, 2.0_rt);
            }
            wn = std::max(std::min(wn, 0.0_rt), -3.0_rt);
            if (wn < 0.0_rt) {
                downw[k][i] = wn;
                downthv[k][i] = thvn;
                downthl[k][i] = thln;
                downqt[k][i] = qtn;
                downqc[k][i] = qcn;
                downu[k][i] = un;
                downv[k][i] = vn;
                downa[k][i] = downa[k + 1][i];
            }
            else {
                if (dd_initk[i] - k < 2) {
                    downw[k][i] = 0.0_rt;
                    downthv[k][i] = 0.0_rt;
                    downthl[k][i] = 0.0_rt;
                    downqt[k][i] = 0.0_rt;
                    downqc[k][i] = 0.0_rt;
                    downu[k][i] = 0.0_rt;
                    downv[k][i] = 0.0_rt;
                    }
                break;
                }
            }
        }
    }
    for (int i = 0; i < ndown; i++) {
      downw[0][i] = 0.0_rt;
      downa[0][i] = 0.0_rt;
    }
    for (int k = qltop; k >= kts; k--) {
        for (int i = 0; i < ndown; i++) {
            edmf_a_dd[k] += downa[k - 1][i];
            edmf_w_dd[k] += downa[k - 1][i] * downw[k - 1][i];
            edmf_qt_dd[k] += downa[k - 1][i] * downqt[k - 1][i];
            edmf_thl_dd[k] += downa[k - 1][i] * downthl[k - 1][i];
            edmf_ent_dd[k] += downa[k - 1][i] * ent[k+1 - 1][i];
            edmf_qc_dd[k] += downa[k - 1][i] * downqc[k - 1][i];
        }
        if (edmf_a_dd[k] > 0.0_rt) {
            edmf_w_dd[k] /= edmf_a_dd[k];
            edmf_qt_dd[k] /= edmf_a_dd[k];
            edmf_thl_dd[k] /= edmf_a_dd[k];
            edmf_ent_dd[k] /= edmf_a_dd[k];
            edmf_qc_dd[k] /= edmf_a_dd[k];
        }
    }
    for (int k = kts; k <= qltop; k++) {
        rho_int = (rho[k] * dz[k + 1] + rho[k + 1] * dz[k]) / (dz[k + 1] + dz[k]);
        for (int i = 0; i < ndown; i++) {
            sd_aw[k] += rho_int * downa[k][i] * downw[k][i];
            sd_awthl[k] += rho_int * downa[k][i] * downw[k][i] * downthl[k][i];
            sd_awqt[k] += rho_int * downa[k][i] * downw[k][i] * downqt[k][i];
            sd_awqc[k] += rho_int * downa[k][i] * downw[k][i] * downqc[k][i];
            sd_awu[k] += rho_int * downa[k][i] * downw[k][i] * downu[k][i];
            sd_awv[k] += rho_int * downa[k][i] * downw[k][i] * downv[k][i];
        }
        sd_awqv[k] = sd_awqt[k] - sd_awqc[k];
    }
}

// assuming Real is equivalent to Real or float. adjust as necessary.
void topdown_cloudrad_cc(int& kts, int& kte, const Real* dz1, const Real* zw, Real& fltv, Real& xland, int& kpbl, Real& pblh, const Real* sqc, const Real* sqi, const Real* sqw, const Real* thl, const Real* th1, const Real* ex1, const Real* p1, const Real*  rho1, const Real* thetav, const Real* cldfra_bl1d, const Real* rthraten, Real& maxkhtopdown, Real* khtopdown, Real* tkeprodtd) {
    // constants
  /*
    const Real pfac = 2.0, zfmin = 0.01, phifac = 8.0;
    const Real grav = 9.81, cp = 1004.0, xlv = 2.5e6, xlvcp = xlv / cp, r_d = 287.0, ep_2 = 0.622, p608 = 0.608, karman = 0.4;
    const Real twothirds = 2.0_rt / 3.0, onethird = 1.0_rt / 3.0;
  */
  //Main meaningful difference is cp=1004.5 vs cp=1004
    const Real pfac = 2.0, zfmin = 0.0099999997764826, phifac = 8.0;
    const Real grav = 9.8100004196166992, cp = 1004.5, xlv = 2.5e6, xlvcp = 2488.8002929687500000, r_d = 287.0, ep_2 = 0.6217504143714905, p608 = 0.6083624362945557, karman = 0.4000000059604645;
    const Real twothirds = 0.6666666865348816, onethird = 0.3333333432674408;
    // local variables
    Real zfac[kte - kts + 1], wscalek2[kte - kts + 1], zfacent[kte - kts + 1];
    Real bfx0, wm3, bfxpbl, dthvx, tmp1;
    Real temps, templ, zl1, wstar3_2;
    Real ent_eff, radsum, radflux, we, rcldb, rvls, minrad = 100., zminrad = pblh;
    int k, kk, kminrad = kpbl;
    bool cloudflg = false;

    maxkhtopdown = 0.0;
    for (kk = kts; kk <= kte; ++kk) {
      khtopdown[kk] = 0.0_rt;
      tkeprodtd[kk] = 0.0_rt;
    }

    // check for stratocumulus-topped boundary layers
    for (kk = std::max(0, kpbl + 1 - 2) - 1; kk <= kpbl + 1 + 3 - 1; ++kk) {
        if (sqc[kk - kts] > 1.e-6 || sqi[kk - kts] > 1.e-6 || cldfra_bl1d[kk - kts] > 0.5_rt) {
            cloudflg = true;
        }
        if (rthraten[kk - kts] < minrad) {
            minrad = rthraten[kk - kts];
            kminrad = kk;
            zminrad = zw[kk - kts] + 0.5_rt * dz1[kk - kts];
        }
    }

    if (std::max(kminrad, kpbl) < 1) cloudflg = false;

    if (cloudflg) {
        zl1 = dz1[kts];
        k = std::max(kpbl - 1, kminrad - 1);
        templ = thl[k - kts] * ex1[k - kts];
        rvls = 100._rt * 6.112_rt * std::exp(17.67_rt * (templ - 273.16) / (templ - 29.65)) * (ep_2 / p1[k + 1 - kts]);
        temps = templ + (sqw[k - kts] - rvls) / (cp / xlv + ep_2 * xlv * rvls / (r_d * std::pow(templ, 2)));
        rvls = 100._rt * 6.112_rt * std::exp(17.67_rt * (temps - 273.15) / (temps - 29.65)) * (ep_2 / p1[k + 1 - kts]);
        rcldb = std::max(sqw[k - kts] - rvls, 0.0_rt);
        dthvx = (thl[k + 2 - kts] + th1[k + 2 - kts] * p608 * sqw[k + 2 - kts]) - (thl[k - kts] + th1[k - kts] * p608 * sqw[k - kts]);
        dthvx = std::max(dthvx, 0.1_rt);
        tmp1 = xlvcp * rcldb / (ex1[k - kts] * dthvx);
        ent_eff = 0.2_rt + 0.2_rt * 8._rt * tmp1;
        radsum = 0.0;
        for (kk = std::max(0, kpbl - 3) ; kk <= kpbl + 1 + 3 - 1; ++kk) {
            radflux = rthraten[kk - kts] * ex1[kk - kts]; // converts theta/s to temp/s
            radflux = radflux * cp / grav * (p1[kk - kts] - p1[kk + 1 - kts]); // converts temp/s to w/m^2
            if (radflux < 0.0_rt) radsum = std::abs(radflux) + radsum;
        }
        if ((xland - 1.5_rt) >= 0) { // water
            radsum = std::min(radsum, 90.0_rt);
            bfx0 = std::max(radsum / rho1[k - kts] / cp, 0.0_rt);
        } else { // land
            radsum = std::min(0.25_rt * radsum, 30.0_rt); // practically turn off over land
            bfx0 = std::max(radsum / rho1[k - kts] / cp - std::max(fltv, 0.0_rt), 0.0_rt);
        }
        wm3 = grav / thetav[k - kts] * bfx0 * std::min(pblh, 1500._rt); // this is wstar3
        bfxpbl = -ent_eff * bfx0;
        dthvx = std::max(thetav[k + 1 - kts] - thetav[k - kts], 0.1_rt);
        we = std::max(bfxpbl / dthvx, -std::sqrt(std::pow(wm3, twothirds)));

        for (kk = kts; kk <= kpbl + 3; ++kk) {
            zfac[kk - kts] = std::min(std::max((1._rt - (zw[kk + 1 - kts] - zl1) / (zminrad - zl1)), zfmin), 1.0_rt);
            zfacent[kk - kts] = 10._rt * std::max((zminrad - zw[kk + 1 - kts]) / zminrad, 0.0_rt) * ((1._rt - zfac[kk - kts])*(1._rt - zfac[kk - kts])*(1._rt - zfac[kk - kts]));
            wscalek2[kk - kts] = std::pow((phifac * karman * wm3 * (zfac[kk - kts])), onethird);
            khtopdown[kk - kts] = wscalek2[kk - kts] * karman * (zminrad - zw[kk + 1 - kts]) * ((1._rt - zfac[kk - kts]) * (1._rt - zfac[kk - kts]) * (1._rt - zfac[kk - kts])); // pfac
            khtopdown[kk - kts] = std::max(khtopdown[kk - kts], 0.0_rt);
            tkeprodtd[kk - kts] = 2._rt * ent_eff * wm3 / std::max(pblh, 100._rt) * zfacent[kk - kts];
            tkeprodtd[kk - kts] = std::max(tkeprodtd[kk - kts], 0.0_rt);
        }
    }
    maxkhtopdown = std::numeric_limits<Real>::min();
    for (kk = kts; kk <= kte; ++kk) {
      maxkhtopdown = std::max(maxkhtopdown,khtopdown[kk]);
    }
}

void scale_aware_cc(Real& dx, Real& pbl1, Real& psig_bl, Real& psig_shcu) {
    Real dxdh;
    psig_bl = 1.0_rt;
    psig_shcu = 1.0_rt;
    dxdh = std::max(2.5_rt * dx, 10.0_rt) / std::min(pbl1, 3000.0_rt);

    // new form to preserve parameterized mixing - only down 5% at dx = 750 m
    psig_bl = ((dxdh * dxdh) + 0.106_rt * std::pow(dxdh, 0.667_rt)) / ((dxdh * dxdh) + 0.066_rt * std::pow(dxdh, 0.667_rt) + 0.071_rt);

    // assume a 500 m cloud depth for shallow-cu clouds
    dxdh = std::max(2.5_rt * dx, 10.0_rt) / std::min(pbl1 + 500.0_rt, 3500.0_rt);

    // hyeyum hailey shin and song-you hong 2013, tke in entrainment zone
    psig_shcu = ((dxdh * dxdh) + 0.145_rt * std::pow(dxdh, 0.667_rt)) / ((dxdh * dxdh) + 0.172_rt * std::pow(dxdh, 0.667_rt) + 0.170_rt);

    // clamping psig_bl and psig_shcu to [0, 1]
    psig_bl = std::max(0.0_rt, std::min(psig_bl, 1.0_rt));
    psig_shcu = std::max(0.0_rt, std::min(psig_shcu, 1.0_rt));
    //    exit(1);
}

// ==================================================================
//>\ingroup gsd_mynn_edmf
// this subroutine calculates hybrid diagnostic boundary-layer height (pblh).
//
// notes on the pblh formulation: the 1.5-theta-increase method defines
//pbl heights as the level at.
//which the potential temperature first exceeds the minimum potential.
//temperature within the boundary layer by 1.5_rt k. when applied to.
//observed temperatures, this method has been shown to produce pbl-
//height estimates that are unbiased relative to profiler-based.
//estimates (nielsen-gammon et al. 2008 \cite nielsen_gammon_2008).
// however, their study did not
//include lljs. banta and pichugina (2008) \cite pichugina_2008  show that a tke-based.
//threshold is a good estimate of the pbl height in lljs. therefore,
//a hybrid definition is implemented that uses both methods, weighting
//the tke-method more during stable conditions (pblh < 400 m).
//a variable tke threshold (tkeeps) is used since no hard-wired
//value could be found to work best in all conditions.
//>\section gen_get_pblh  gsd get_pblh general algorithm
//> @{
void get_pblh_cc(int &kts, int &kte, Real &zi, Real* thetav1d, Real *qke1d, Real *zw1d, Real* dz1d, Real &landsea, int &kzi) {
    // HR: SEGFAULTS WHEN ACCESSING VECTORS, need to look into how to pass 1d arrays to c++ from fortran
    // constants
    const Real sbl_lim = 200.0;
    const Real sbl_damp = 400.0;

    // local variables
    Real pblh_tke, qtke, qtkem1, maxqke, tkeeps, minthv, delt_thv;
    int kthv, ktke;

    // initialize kpbl (kzi)
    kzi = 2;
    // find min thetav in the lowest 200 m agl
    kthv = 1;
    minthv = 9e9;
    for (int k = kts + 1; k <= kte && zw1d[k - kts] <= 200.; ++k) {
        if (minthv > thetav1d[k - kts]) {
            minthv = thetav1d[k - kts];
            kthv = k;
        }
    }

    // find thetav-based pblh (best for daytime)
    zi = 0.0;
    delt_thv = (landsea - 1.5_rt) >= 0 ? 1.0_rt : 1.25_rt;
    //for (int k = kthv + 1; k < kte; ++k) {
    for (int k = kts + 1; k < kte; ++k) {
        if (thetav1d[k - kts] >= (minthv + delt_thv)) {
            zi = zw1d[k - kts] - dz1d[k - 1 - kts] * std::min((thetav1d[k - kts] - (minthv + delt_thv)) / std::max(thetav1d[k - kts] - thetav1d[k - 1 - kts], 1.e-6_rt), 1.0_rt);
            if (zi > 0) break;
        }
        if (k == kte - 1) zi = zw1d[kts + 1 - kts]; // exit safeguard
    }

    // for stable boundary layers, use tke method to complement the thetav-based definition
    pblh_tke = 0.0;
    maxqke = std::max(qke1d[kts - kts], 0.0_rt);
    tkeeps = maxqke / 40.0;
    tkeeps = std::max(tkeeps, 0.02_rt);

    for (int k = kts + 1; k < kte; ++k) {
        qtke = std::max(qke1d[k - kts] / 2.0_rt, 0.0_rt);
        qtkem1 = std::max(qke1d[k - 1 - kts] / 2.0_rt, 0.0_rt);
        if (qtke <= tkeeps) {
            pblh_tke = zw1d[k - kts] - dz1d[k - 1 - kts] * std::min((tkeeps - qtke) / std::max(qtkem1 - qtke, 1.0e-6_rt), 1.0_rt);
            pblh_tke = std::max(pblh_tke, zw1d[kts + 1 - kts]);
            break;
        }
        if (k == kte - 1) pblh_tke = zw1d[kts + 1 - kts]; // exit safeguard
    }

    // limit pblh_tke to not exceed the thetav-based pbl height +/- 350 m
    pblh_tke = std::min(pblh_tke, zi + 350.0_rt);
    pblh_tke = std::max(pblh_tke, std::max(zi - 350.0_rt, 10.0_rt));

    Real wt = 0.5_rt * std::tanh((zi - sbl_lim) / sbl_damp) + 0.5;
    if (maxqke > 0.05) {
        zi = pblh_tke * (1.0_rt - wt) + zi * wt;
    }

    // compute kpbl (kzi)
    for (int k = kts + 1; k < kte; ++k) {
        if (zw1d[k - kts] >= zi) {
            kzi = k - 1;
            break;
        }
    }
}

void retrieve_exchange_coeffs_cc(int& kts, int& kte, Real* dfm, Real* dfh, const Real* dz, Real* k_m, Real* k_h) {
    Real dzk;
    k_m[kts] = 0.0;
    k_h[kts] = 0.0;
    for (int k = kts + 1; k <= kte; ++k) {
        dzk = 0.5_rt * (dz[k] + dz[k - 1]);
        k_m[k] = dfm[k] * dzk;
        k_h[k] = dfh[k] * dzk;
    }
}


// ==================================================================
//     subroutine  mym_length:
//
//     input variables:    see subroutine mym_initialize
//
//     output variables:   see subroutine mym_initialize
//
//     work arrays:
//       elt(nx,ny)      : length scale depending on the pbl depth    (m)
//       vsc(nx,ny)      : velocity scale q_c                       (m/s)
//                         at first, used for computing elt
//
//     note: the mixing lengths are meant to be calculated at the full-
//           sigmal levels (or interfaces between the model layers).
//
//>\ingroup gsd_mynn_edmf
// this subroutine calculates the mixing lengths.

void mym_length(int kts, int kte, Real xland, Real* dz, Real* dx, Real* zw, Real rmo, Real flt, Real fltv, Real flq, Real* vt, Real* vq, Real* u1, Real* v1, Real* qke, Real* dtv, Real* el, Real zi, Real* theta, Real* qkw, Real psig_bl, Real* cldfra_bl1d, int bl_mynn_mixlength, Real* edmf_w1, Real* edmf_a1, Real grav, Real karman, Real tv0, Real gtr) {
    int i, j, k;
    Real elt, vsc;
    Real qtke[kte+1], elblmin[kte+1], elblavg[kte+1], thetaw[kte+1];
    Real wt, wt2, zi2, h1, h2, hs, elblmin0, elblavg0, cldavg;
    Real cns, alp1, alp2, alp3, alp4, alp5, alp6;
    Real minzi = 300.0;
    Real maxdz = 750.0;
    Real mindz = 300.0;
    Real zslh = 100.0;
    Real csl = 2.0;
    Real qke_elb_min = 0.018_rt;
    Real afk, abk, zwk, zwk1, dzk, qdz, vflx, bv, tau_cloud, wstar, elb, els, elf, el_stab, el_mf, el_stab_mf, elb_mf, pblh_plus_ent, uonset, ugrid, wt_u, el_les;
    Real ctau = 1000.0;

    switch(bl_mynn_mixlength) {
        case 0:
            cns = 2.7;
            alp1 = 0.23;
            alp2 = 1.0;
            alp3 = 5.0;
            alp4 = 100.0;
            alp5 = 0.3;
            zi2 = std::min(10000.0_rt, zw[kte-2]);
            h1 = std::max(0.3_rt * zi2, mindz);
            h1 = std::min(h1, maxdz);
            h2 = h1 / 2.0;
            qkw[kts] = std::sqrt(std::max(qke[kts], 1.0e-10_rt));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0_rt - afk;
                qkw[k] = std::sqrt(std::max(qke[k] * abk + qke[k-1] * afk, 1.0e-3_rt));
            }
            elt = 1.0e-5;
            vsc = 1.0e-5;
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi2 + h1) {
                dzk = 0.5_rt * (dz[k] + dz[k-1]);
                qdz = std::max(qkw[k] - qmin, 0.03_rt) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }
            elt = alp1 * elt / vsc;
            vflx = (vt[kts] + 1.0_rt) * flt + (vq[kts] + tv0) * flq;
            vsc = std::pow(gtr * elt * std::max(vflx, 0.0_rt), 1.0_rt / 3.0_rt);
            el[kts] = 0.0;
            zwk1 = zw[kts+1];
            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];
                if (dtv[k] > 0.0_rt) {
                    bv = std::sqrt(gtr * dtv[k]);
                    elb = alp2 * qkw[k] / bv * (1.0_rt + alp3 / alp2 * std::sqrt(vsc / (bv * elt)));
                    elf = alp2 * qkw[k] / bv;
                } else {
                    elb = 1.0e10;
                    elf = elb;
                }
                if (rmo > 0.0_rt) {
                    els = karman * zwk / (1.0_rt + cns * std::min(zwk * rmo, zmax));
                } else {
                    els = karman * zwk * std::pow(1.0_rt - alp4 * zwk * rmo, 0.2_rt);
                }
                wt = 0.5_rt * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5;
                el[k] = std::min(elb / (elb / elt + elb / els + 1.0_rt), elf);
            }
            break;
        case 1:
            ugrid = std::sqrt(u1[kts] * u1[kts] + v1[kts] * v1[kts]);
            uonset = 15.0;
            wt_u = (1.0_rt - std::min(std::max(ugrid - uonset, 0.0_rt) / 30.0_rt, 0.5_rt));
            cns = 2.7;
            alp1 = 0.23;
            alp2 = 0.3;
            alp3 = 2.5_rt * wt_u;
            alp4 = 5.0;
            alp5 = 0.3;
            alp6 = 50.0;
            zi2 = std::max(zi, 300.0_rt);
            h1 = std::max(0.3_rt * zi2, 300.0_rt);
            h1 = std::min(h1, 600.0_rt);
            h2 = h1 / 2.0;
            qkw[kts] = std::sqrt(std::max(qke[kts], 1.0e-10_rt));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0_rt - afk;
                qkw[k] = std::sqrt(std::max(qke[k] * abk + qke[k-1] * afk, 1.0e-3_rt));
                qtke[k] = 0.5_rt * (qkw[k] * qkw[k]);
                thetaw[k] = theta[k] * abk + theta[k-1] * afk;
            }
            elt = 1.0e-5;
            vsc = 1.0e-5;
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi2 + h1) {
                dzk = 0.5_rt * (dz[k] + dz[k-1]);
                qdz = std::min(std::max(qkw[k] - qmin, 0.01_rt), 30.0_rt) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }
            elt = std::min(std::max(alp1 * elt / vsc, 10.0_rt), 400.0_rt);
            vflx = fltv;
            vsc = std::pow(gtr * elt * std::max(vflx, 0.0_rt), 1.0_rt / 3.0_rt);
            el[kts] = 0.0;
            zwk1 = zw[kts+1];
            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];
                if (dtv[k] > 0.0_rt) {
                    bv = std::max(std::sqrt(gtr * dtv[k]), 0.0001_rt);
                    elb = std::max(alp2 * qkw[k], alp6 * edmf_a1[k-1] * edmf_w1[k-1]) / bv * (1.0_rt + alp3 * std::sqrt(vsc / (bv * elt)));
                    elb = std::min(elb, zwk);
                    elf = 1.0_rt * qkw[k] / bv;
                    elblavg[k] = std::max(elblavg[k], alp6 * edmf_a1[k-1] * edmf_w1[k-1] / bv);
                } else {
                    elb = 1.0e10;
                    elf = elb;
                }
                if (rmo > 0.0_rt) {
                    els = karman * zwk / (1.0_rt + cns * std::min(zwk * rmo, zmax));
                } else {
                    els = karman * zwk * std::pow(1.0_rt - alp4 * zwk * rmo, 0.2_rt);
                }
                wt = 0.5_rt * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5;
                el[k] = std::min(elb / (elb / elt + elb / els + 1.0_rt), elf);
                el[k] = el[k] * psig_bl;
            }
            break;
        case 2:
            uonset = 3.5_rt + dz[kts] * 0.1;
            ugrid = std::sqrt(u1[kts] * u1[kts] + v1[kts] * v1[kts]);
            cns = 3.5;
            alp1 = 0.22;
            alp2 = 0.30;
            alp3 = 2.0;
            alp4 = 5.0;
            alp5 = alp2;
            alp6 = 50.0;
            zi2 = std::max(zi, minzi);
            h1 = std::max(0.3_rt * zi2, mindz);
            h1 = std::min(h1, maxdz);
            h2 = h1 * 0.5;
            qtke[kts] = std::max(0.5_rt * qke[kts], 0.5_rt*qkemin);
            qkw[kts] = std::sqrt(std::max(qke[kts], qkemin));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0_rt - afk;
                qkw[k] = std::sqrt(std::max(qke[k] * abk + qke[k-1] * afk, qkemin));
                qtke[k] = 0.5_rt * qkw[k] * qkw[k];
                thetaw[k] = theta[k] * abk + theta[k-1] * afk;
            }
            elt = 1.0e-5;
            vsc = 1.0e-5;
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi2 + h1) {
                dzk = 0.5_rt * (dz[k] + dz[k-1]);
                qdz = std::min(std::max(qkw[k] - qmin, 0.03_rt), 30.0_rt) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }
            elt = std::min(std::max(alp1 * elt / vsc, 10.0_rt), 400.0_rt);
            vflx = fltv;
            vsc = std::pow(gtr * elt * std::max(vflx, 0.0_rt), 1.0_rt / 3.0_rt);
            el[kts] = 0.0;
            zwk1 = zw[kts+1];
            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];
                dzk = 0.5_rt * (dz[k] + dz[k-1]);
                cldavg = 0.5_rt * (cldfra_bl1d[k-1] + cldfra_bl1d[k]);
                if (dtv[k] > 0.0_rt) {
                    bv = std::max(std::sqrt(gtr * dtv[k]), 0.001_rt);
                    elb_mf = std::max(alp2 * qkw[k], alp6 * edmf_a1[k-1] * edmf_w1[k-1]) / bv * (1.0_rt + alp3 * std::sqrt(vsc / (bv * elt)));
                    elb = std::min(std::max(alp5 * qkw[k], alp6 * edmf_a1[k] * edmf_w1[k]) / bv, zwk);
                    wstar = 1.25_rt * std::pow(gtr * zi * std::max(vflx, 1.0e-4_rt), 1.0_rt / 3.0_rt);
                    tau_cloud = std::min(std::max(ctau * wstar / grav, 30.0_rt), 150.0_rt);
                    wt = 0.5_rt * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5;
                    tau_cloud = tau_cloud * (1.0_rt - wt) + 50.0_rt * wt;
                    elf = std::min(std::max(tau_cloud * std::sqrt(std::min(qtke[k], 40.0_rt)), alp6 * edmf_a1[k] * edmf_w1[k] / bv), zwk);
                } else {
                    wstar = 1.25_rt * std::pow(gtr * zi * std::max(vflx, 1.0e-4_rt), 1.0_rt / 3.0_rt);
                    tau_cloud = std::min(std::max(ctau * wstar / grav, 50.0_rt), 200.0_rt);
                    wt = 0.5_rt * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5;
                    tau_cloud = tau_cloud * (1.0_rt - wt) + std::max(100.0, dzk * 0.25) * wt;
                    elb = std::min(tau_cloud * std::sqrt(std::min(qtke[k], 40.0_rt)), zwk);
                    elf = elb;
                    elb_mf = elb;
                }
                elf = elf / (1.0_rt + (elf / 800.0_rt));
                elb_mf = std::max(elb_mf, 0.01_rt);
                if (rmo > 0.0_rt) {
                    els = karman * zwk / (1.0_rt + cns * std::min(zwk * rmo, zmax));
                } else {
                    els = karman * zwk * std::pow(1.0_rt - alp4 * zwk * rmo, 0.2_rt);
                }
                wt = 0.5_rt * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5;
                el[k] = std::sqrt(els * els / (1.0_rt + (els * els / elt * elt) + (els * els / elb_mf * elb_mf)));
                el[k] = el[k] * (1.0_rt - wt) + elf * wt;
                el_les = std::min(els / (1.0_rt + (els / 12.0_rt)), elb_mf);
                el[k] = el[k] * psig_bl + (1.0_rt - psig_bl) * el_les;
            }
            break;
    }
}

// ==================================================================
//>\ingroup gsd_mynn_edmf
// this subroutine is the dynamic multi-plume (dmp) mass-flux scheme.
//
// dmp_mf() calculates the nonlocal turbulent transport from the dynamic
// multiplume mass-flux scheme as well as the shallow-cumulus component of
// the subgrid clouds. note that this mass-flux scheme is called when the
// namelist parameter \p bl_mynn_edmf is set to 1 (recommended).
//
// much thanks to kay suslj of nasa-jpl for contributing the original version
// of this mass-flux scheme. considerable changes have been made from it's
// original form. some additions include:
//  -# scale-aware tapering as dx -> 0
//  -# transport of tke (extra namelist option)
//  -# chaboureau-bechtold cloud fraction & coupling to radiation (when icloud_bl > 0)
//  -# some extra limits for numerical stability
//
// this scheme remains under development, so consider it experimental code.
//

void dmp_mf_cc(const int& kts,const int& kte, Real& dt, Real* zw, Real* dz, Real* p, Real* rho, int& momentum_opt, int& tke_opt, int& scalar_opt, Real* u, Real* v, Real* w, Real* th, Real* thl, Real* thv, Real* tk, Real* qt, Real* qv, Real* qc, Real* qke, Real* qnc, Real* qni, Real* qnwfa, Real* qnifa, Real* qnbca, Real& ust, Real& flt, Real& fltv, Real& flq, Real& flqv, Real& pblh, int& kpbl, Real& dx, Real& landsea, Real& ts, Real* edmf_a, Real* edmf_w, Real* edmf_qt, Real* edmf_thl, Real* edmf_ent, Real* edmf_qc, Real* s_aw, Real* s_awthl, Real* s_awqt, Real* s_awqv, Real* s_awqc, Real* s_awu, Real* s_awv, Real* s_awqke, Real* s_awqnc, Real* s_awqni, Real* s_awqnwfa, Real* s_awqnifa, Real* s_awqnbca, int& nchem, Real** chem1, Real** s_awchem, bool& mix_chem, Real* qc_bl1d, Real* cldfra_bl1d, Real* qc_bl1d_old, Real* cldfra_bl1d_old, Real& psig_shcu, Real& maxwidth, int& ktop, Real& maxmf, Real& ztop, Real* rstoch_col, Real grav, Real gtr, Real p608) {
    int nup = 8;
    int debug_mf = 0;
    Real nup2;
    Real upw[kte+1+1][nup], upthl[kte+1+1][nup], upqt[kte+1+1][nup], upqc[kte+1+1][nup], upqv[kte+1+1][nup], upa[kte+1+1][nup], upu[kte+1+1][nup], upv[kte+1+1][nup], upthv[kte+1+1][nup], upqke[kte+1+1][nup], upqnc[kte+1+1][nup], upqni[kte+1+1][nup], upqnwfa[kte+1+1][nup], upqnifa[kte+1+1][nup], upqnbca[kte+1+1][nup];
    Real ent[kte+1][nup];
    int enti[kte+1][nup];
    int k, i, k50;
    Real fltv2, wstar, qstar, thstar, sigmaw, sigmaqt, sigmath, z0, pwmin, pwmax, wmin, wmax, wlv, psig_w, maxw, wpbl;
    Real b, qtn, thln, thvn, qcn, un, vn, qken, qncn, qnin, qnwfan, qnifan, qnbcan, wn2, wn, entexp, entexm, entw, bcoeff, thvkm1, thvk, pk, rho_int;
    Real wa = 0.6666666865348816; //2./3.
    Real wb = 0.0020000000949949, wc = 1.5;
    Real l0 = 100., ent0 =  0.1000000014901161;
    Real atot = 0.1000000014901161;
    Real lmax = 1000.;
    Real lmin = 300.;
    Real dlmin = 0.;
    Real minwidth;
    Real dl;
    Real dcut = 1.2000000476837158;
    Real d;
    Real cn, c, l, n, an2, hux, wspd_pbl, cloud_base, width_flx;
    Real chemn[nchem];
    Real upchem[kte+1+1][nup][nchem];
    Real ic;
    Real edmf_chem[kte+1+1][nchem];
    Real envm_u[kte+1+1],envm_v[kte+1+1],envm_sqc[kte+1+1],envm_thl[kte+1+1],envm_sqv[kte+1+1];
    bool superadiabatic;
    Real sigq, xl, rsl, cpm, a, qmq, mf_cf, aup, q1, diffqt, qsat_tk, fng, qww, alpha, beta, bb, f, pt, t, q2p, b9, satvp, rhgrid, ac_mf, ac_strat, qc_mf;
    Real cf_thresh = 0.5;
    Real exneri[kte+1], dzi[kte+1], rhoz[kte+1];
    Real thp, qtp, qcp, qcs, esat, qsl;
    Real csigma, acfac, ac_wsp;
    int overshoot;
    Real bvf, frz, dzp;
    Real adjustment, flx1, flt2;
    Real fluxportion = 0.75;
    Real sublim, qc_ent, qv_ent, qt_ent, thl_ent, detrate, detrateuv, oow, exc_fac, aratio, detturb, qc_grid, qc_sgs, exc_heat, exc_moist, tk_int, tvs, qc_plume;
    Real cdet = 0.0222222227603197;//1./45.;
    Real dzpmax = 300.;
    Real csub = 0.25;
    Real pgfac = 0.00;
    Real uk, ukm1, vk, vkm1, dxsa;

    for (int i = 0; i < nup; i++) {
        for (int j = kts; j <= kte+1; j++) {
            upw[j][i] = 0.0;
            upthl[j][i] = 0.0;
            upqt[j][i] = 0.0;
            upqc[j][i] = 0.0;
            upqv[j][i] = 0.0;
            upa[j][i] = 0.0;
            upu[j][i] = 0.0;
            upv[j][i] = 0.0;
            upthv[j][i] = 0.0;
            upqke[j][i] = 0.0;
            upqnc[j][i] = 0.0;
            upqni[j][i] = 0.0;
            upqnwfa[j][i] = 0.0;
            upqnifa[j][i] = 0.0;
            upqnbca[j][i] = 0.0;
            if(mix_chem)
              for(int nchemi = 0; nchemi < nchem; nchemi++)
                upchem[j][i][nchemi] = 0.0;
        }
    }
    for (int i = kts; i <= kte; i++) {
        for (int j = 0; j < nup; j++) {
            ent[i][j] = 0.001;
            enti[i][j] = 0;
        }
    }
    for (int k = kts; k <= kte; k++) {
      edmf_a[k] = 0.0;
      edmf_w[k] = 0.0;
      edmf_qt[k] = 0.0;
      edmf_thl[k] = 0.0;
      edmf_ent[k] = 0.0;
      edmf_qc[k] = 0.0;
      if(mix_chem)
        for(int nchemi = 0; nchemi < nchem; nchemi++) {
          edmf_chem[k][nchemi] = 0.0;
        }
    }

    for (int k = kts; k <= kte+1; k++) {
      s_aw[k] = 0.0;
      s_awthl[k] = 0.0;
      s_awqt[k] = 0.0;
      s_awqv[k] = 0.0;
      s_awqc[k] = 0.0;
      s_awu[k] = 0.0;
      s_awv[k] = 0.0;
      s_awqke[k] = 0.0;
      s_awqnc[k] = 0.0;
      s_awqni[k] = 0.0;
      s_awqnwfa[k] = 0.0;
      s_awqnifa[k] = 0.0;
      s_awqnbca[k] = 0.0;
      if(mix_chem)
        for(int nchemi = 0; nchemi < nchem; nchemi++)
          s_awchem[k][nchemi] = 0.0;
    }
    //    printf("missing sub_thl in dmp_mf\n");
    //   exit(1);
    /*
    for (int k = kts; k <= kte; k++) {
      sub_thl[k] = 0.0;
      sub_sqv[k] = 0.0;
      sub_u[k] = 0.0;
      sub_v[k] = 0.0;
      det_thl[k] = 0.0;
      det_sqv[k] = 0.0;
      det_u[k] = 0.0;
      det_v[k] = 0.0;
      nup2[k] = nup[k];
    }
    */
    maxw = 0.0;
    cloud_base = 9000.0;
    for (int k = kts; k < kte-1; k++) {
        if (zw[k] > pblh + 500._rt) {
            break;
        }
        wpbl = w[k];
        if (w[k] < 0._rt) {
            wpbl = 2.*w[k];
        }
        maxw = std::max(maxw, Real(std::abs(wpbl)));
        if (zw[k] <= 50._rt) {
            k50 = k;
        }
        Real qc_sgs = std::max(qc[k], qc_bl1d[k]);
        if (qc_sgs > 1e-5 && cldfra_bl1d[k] >= 0.5_rt && cloud_base == 9000.0_rt) {
            cloud_base = 0.5_rt*(zw[k]+zw[k+1]);
        }
    }
    maxw = std::max(0._rt, maxw - 1.0_rt);
    psig_w = std::max(0.0_rt, 1.0_rt - maxw);
    psig_w = std::min(psig_w, psig_shcu);
    fltv2 = fltv;
    if (psig_w == 0.0_rt && fltv > 0.0_rt) {
        fltv2 = -1.*fltv;
    }
    superadiabatic = false;
    tvs = ts*(1.0+p608*qv[kts]);
    for (int k = 0; k < std::max(1, k50-1); k++) {
        if (k == 0) {
            if ((thv[k]-tvs)/(0.5_rt*dz[k]) < hux) {
                superadiabatic = true;
            } else {
                superadiabatic = false;
                break;
            }
        } else {
            if ((thv[k]-thv[k-1])/(0.5_rt*(dz[k]+dz[k-1])) < hux) {
                superadiabatic = true;
            } else {
                superadiabatic = false;
                break;
            }
        }
    }
    maxwidth = std::min(dx*dcut, lmax);
    maxwidth = std::min(maxwidth, 1.1_rt*pblh);
    if (landsea-1.5_rt < 0) {
        maxwidth = std::min(maxwidth, 0.5_rt*cloud_base);
    } else {
        maxwidth = std::min(maxwidth, 0.9_rt*cloud_base);
    }
    wspd_pbl = std::sqrt(std::max(u[kts]*u[kts] + v[kts]*v[kts], 0.01_rt));
    if (landsea-1.5_rt < 0) {
        width_flx = std::max(std::min(1000._rt*(0.6_rt*Real(tanh((fltv - 0.040_rt)/0.04_rt)) + .5_rt),1000._rt), 0._rt);
    } else {
        width_flx = std::max(std::min(1000._rt*(0.6_rt*Real(tanh((fltv - 0.007_rt)/0.02_rt)) + .5_rt),1000._rt), 0._rt);
    }
    maxwidth = std::min(maxwidth, width_flx);
    minwidth = lmin;
    if (maxwidth >= (lmax - 1.0_rt) && fltv > 0.2_rt) {
        minwidth = lmin + dlmin*std::min((fltv-0.2_rt)/0.3_rt, 1._rt);
    }
    if (maxwidth <= minwidth) {
        nup2 = 0;
        maxwidth = 0.0;
    }
    ktop = 0;
    ztop = 0.0;
    maxmf = 0.0;
    if (fltv2 > 0.002 && maxwidth > minwidth && superadiabatic) {
        Real cn = 0.;
        Real d = -1.9;
        Real dl = (maxwidth - minwidth)/Real(nup-1);
        for (int i = 0; i < nup; i++) {
            Real l = minwidth + dl*Real(i);
            cn = cn + l*l*l * (l*l)/(dx*dx) * dl;
        }
        Real c = atot/cn;
        Real acfac;
        acfac = 0.5_rt*tanh((fltv2 - 0.02_rt)/0.05_rt) + 0.5_rt;

        Real ac_wsp;
        if (wspd_pbl <= 10._rt) {
            ac_wsp = 1.0;
        } else {
            ac_wsp = 1.0_rt - std::min((wspd_pbl - 10.0_rt)/15._rt, 1.0_rt);
        }
        acfac = acfac * ac_wsp;
        Real an2 = 0.;
        for (int i = 0; i < nup; i++) {
            Real l = minwidth + dl*Real(i);
            Real n = c*l*l*l * (l*l)/(dx*dx) * dl;
            upa[kts][i] = n*l*l/(dx*dx) * dl;
            upa[kts][i] = upa[kts][i]*acfac;
            an2 = an2 + upa[kts][i];
        }
        Real z0 = 50.;
        Real pwmin = 0.1;
        Real pwmax = 0.4;
        Real wstar = std::max(1.e-2, pow(gtr*fltv2*pblh, 1./3._rt));
        Real qstar = std::max(flq, 1.0e-5_rt)/wstar;
        Real thstar = flt/wstar;
        Real csigma;
        if (landsea-1.5_rt >= 0) {
            csigma = 1.34;
        } else {
            csigma = 1.34;
        }
        Real exc_fac;
        if (env_subs) {
            exc_fac = 0.0;
        } else {
            if (landsea-1.5_rt >= 0) {
                exc_fac = 0.58*4.0;
            } else {
                exc_fac = 0.58;
            }
        }
        exc_fac = exc_fac * ac_wsp;
        Real sigmaw = csigma*wstar*pow(z0/pblh, 1./3._rt)*(1 - 0.8*z0/pblh);
        Real sigmaqt = csigma*qstar*pow(z0/pblh, 1./3._rt);
        Real sigmath = csigma*thstar*pow(z0/pblh, 1./3._rt);
        Real wmin = std::min(sigmaw*pwmin, 0.1_rt);
        Real wmax = std::min(sigmaw*pwmax, 0.5_rt);
        for (int i = 0; i < nup; i++) {
            Real wlv = wmin+(wmax-wmin)/nup2*Real(i);
            upw[kts][i] = wmin + Real(i+1)/Real(nup)*(wmax-wmin);
            upu[kts][i] = (u[kts]*dz[kts+1]+u[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]);
            upv[kts][i] = (v[kts]*dz[kts+1]+v[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]);
            upqc[kts][i] = 0.0;
            Real exc_heat = exc_fac*upw[kts][i]*sigmath/sigmaw;
            upthv[kts][i] = (thv[kts]*dz[kts+1]+thv[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]) + exc_heat;
            upthl[kts][i] = (thl[kts]*dz[kts+1]+thl[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]) + exc_heat;
            Real exc_moist = exc_fac*upw[kts][i]*sigmaqt/sigmaw;
            upqt[kts][i] = (qt[kts]*dz[kts+1]+qt[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]) + exc_moist;
            upqke[kts][i] = (qke[kts]*dz[kts+1]+qke[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]);
            upqnc[kts][i] = (qnc[kts]*dz[kts+1]+qnc[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]);
            upqni[kts][i] = (qni[kts]*dz[kts+1]+qni[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]);
            upqnwfa[kts][i] = (qnwfa[kts]*dz[kts+1]+qnwfa[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]);
            upqnifa[kts][i] = (qnifa[kts]*dz[kts+1]+qnifa[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]);
            upqnbca[kts][i] = (qnbca[kts]*dz[kts+1]+qnbca[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]);
        }
        if (mix_chem) {
            for (int i = 0; i < nup; i++) {
                for (int ic = 0; ic < nchem; ic++) {
                    upchem[kts][i][ic] = (chem1[kts][ic]*dz[kts+1]+chem1[kts+1][ic]*dz[kts])/(dz[kts]+dz[kts+1]);
                }
            }
        }
        for (int ii=kts;ii>=kte;ii++){
            envm_thl[ii-kts] = thl[ii];
            envm_sqv[ii-kts] = qv[ii];
            envm_sqc[ii-kts] = qc[ii];
            envm_u[ii-kts] = u[ii];
            envm_v[ii-kts] = v[ii];
        }
        for (int k = kts; k < kte-1; k++) {
            rhoz[k] = (rho[k]*dz[k+1]+rho[k+1]*dz[k])/(dz[k+1]+dz[k]);
        }
        rhoz[kte] = rho[kte];
        dxsa = 1._rt - std::min(std::max((12000.0_rt-dx)/(12000.0_rt-3000.0_rt), 0._rt), 1._rt);
        for (int i = 0; i < nup; i++) {
            Real qc_ent = 0.0;
            int overshoot = 0;
            Real l = minwidth + dl*Real(i);
            for (int k = kts+1; k < kte-1; k++) {
                Real wmin = 0.3_rt + l*0.0005;
                ent[k][i] = 0.33/(std::min(std::max(upw[k-1][i], wmin), 0.9_rt)*l);
                ent[k][i] = std::max(ent[k][i], 0.0003_rt);
                if (zw[k] >= std::min(pblh+1500._rt, 4000._rt)) {
                    ent[k][i] = ent[k][i] + (zw[k]-std::min(pblh+1500._rt,4000._rt))*5.0e-6_rt;
                }
                ent[k][i] = ent[k][i] * (1.0_rt - rstoch_col[k]);
                ent[k][i] = std::min(ent[k][i], 0.9_rt/(zw[k+1]-zw[k]));
                Real uk = (u[k]*dz[k+1]+u[k+1]*dz[k])/(dz[k+1]+dz[k]);
                Real ukm1 = (u[k-1]*dz[k]+u[k]*dz[k-1])/(dz[k-1]+dz[k]);
                Real vk = (v[k]*dz[k+1]+v[k+1]*dz[k])/(dz[k+1]+dz[k]);
                Real vkm1 = (v[k-1]*dz[k]+v[k]*dz[k-1])/(dz[k-1]+dz[k]);
                Real entexp = ent[k][i]*(zw[k+1]-zw[k]);
                Real entexm = entexp*0.3333;
                Real qtn = upqt[k-1][i]*(1.-entexp) + qt[k]*entexp;
                Real thln = upthl[k-1][i]*(1.-entexp) + thl[k]*entexp;
                Real un = upu[k-1][i]*(1.-entexm) + u[k]*entexm + dxsa*pgfac*(uk - ukm1);
                Real vn = upv[k-1][i]*(1.-entexm) + v[k]*entexm + dxsa*pgfac*(vk - vkm1);
                Real qken = upqke[k-1][i]*(1.-entexp) + qke[k]*entexp;
                Real qncn = upqnc[k-1][i]*(1.-entexp) + qnc[k]*entexp;
                Real qnin = upqni[k-1][i]*(1.-entexp) + qni[k]*entexp;
                Real qnwfan = upqnwfa[k-1][i]*(1.-entexp) + qnwfa[k]*entexp;
                Real qnifan = upqnifa[k-1][i]*(1.-entexp) + qnifa[k]*entexp;
                Real qnbcan = upqnbca[k-1][i]*(1.-entexp) + qnbca[k]*entexp;
                Real qc_ent = qcn;
                Real qt_ent = qtn;
                Real thl_ent = thln;
                if (mix_chem) {
                    for (int ic = 0; ic < nchem; ic++) {
                        chemn[ic] = upchem[k-1][i][ic]*(1.-entexp) + chem1[k][ic]*entexp;
                    }
                }
                Real pk = (p[k]*dz[k+1]+p[k+1]*dz[k])/(dz[k+1]+dz[k]);
//                condensation_edmf_cc(qtn, thln, pk, zw[k+1], thvn, qcn);
                Real thvk = (thv[k]*dz[k+1]+thv[k+1]*dz[k])/(dz[k+1]+dz[k]);
                Real thvkm1 = (thv[k-1]*dz[k]+thv[k]*dz[k-1])/(dz[k-1]+dz[k]);
                Real b = grav*(thvn/thvk - 1.0_rt);
                if (b > 0._rt) {
                    bcoeff = 0.15;
                } else {
                    bcoeff = 0.2;
                }
                if (upw[k-1][i] < 0.2_rt) {
                    wn = upw[k-1][i] + (-2._rt * ent[k][i] * upw[k-1][i] + bcoeff*b / std::max(upw[k-1][i], 0.2_rt)) * std::min(zw[k]-zw[k-1], 250._rt);
                } else {
                    wn = upw[k-1][i] + (-2._rt * ent[k][i] * upw[k-1][i] + bcoeff*b / upw[k-1][i]) * std::min(zw[k]-zw[k-1], 250._rt);
                }
                if (wn > upw[k-1][i] + std::min(1.25_rt*(zw[k]-zw[k-1])/200._rt, 2.0_rt)) {
                    wn = upw[k-1][i] + std::min(1.25_rt*(zw[k]-zw[k-1])/200._rt, 2.0_rt);
                }
                if (wn < upw[k-1][i] - std::min(1.25_rt*(zw[k]-zw[k-1])/200._rt, 2.0_rt)) {
                    wn = upw[k-1][i] - std::min(1.25_rt*(zw[k]-zw[k-1])/200._rt, 2.0_rt);
                }
                wn = std::min(std::max(wn, 0.0_rt), 3.0_rt);
                if (k == kts+1 && wn == 0._rt) {
                    nup2 = 0;
                    break;
                }
                if (debug_mf == 1) {
                    if (wn >= 3.0_rt) {
                        std::cout << "**** suspiciously large w:" << std::endl;
                        std::cout << "qcn: " << qcn << " ent: " << ent[k][i] << " nup2: " << nup2 << std::endl;
                        std::cout << "pblh: " << pblh << " wn: " << wn << " upw(k-1): " << upw[k-1][i] << std::endl;
                        std::cout << "k: " << k << " b: " << b << " dz: " << zw[k]-zw[k-1] << std::endl;
                    }
                }
                if (wn <= 0.0_rt && overshoot == 0) {
                    overshoot = 1;
                    if (thvk-thvkm1 > 0.0_rt) {
                        Real bvf = std::sqrt(gtr*(thvk-thvkm1)/dz[k]);
                        Real frz = upw[k-1][i]/(bvf*dz[k]);
                        dzp = dz[k]*std::max(std::min(frz, 1.0_rt), 0.0_rt);
                    }
                } else {
                    dzp = dz[k];
                }
                Real aratio = std::min(upa[k-1][i]/(1._rt-upa[k-1][i]), 0.5_rt);
                Real detturb = 0.00008;
                Real oow = -0.060/std::max(1.0_rt, (0.5_rt*(wn+upw[k-1][i])));
                Real detrate = std::min(std::max(oow*(wn-upw[k-1][i])/dz[k], detturb), 0.0002_rt);
                Real detrateuv = std::min(std::max(oow*(wn-upw[k-1][i])/dz[k], detturb), 0.0001_rt);
                envm_thl[k-kts] = envm_thl[k-kts] + (0.5_rt*(thl_ent + upthl[k-1][i]) - thl[k])*detrate*aratio*std::min(dzp, dzpmax);
                Real qv_ent = 0.5_rt*(std::max(qt_ent-qc_ent, 0.0_rt) + std::max(upqt[k-1][i]-upqc[k-1][i], 0.0_rt));
                envm_sqv[k-kts] = envm_sqv[k] + (qv_ent-qv[k])*detrate*aratio*std::min(dzp, dzpmax);
                if (upqc[k-1][i] > 1e-8) {
                    Real qc_grid;
                    if (qc[k] > 1e-6) {
                        qc_grid = qc[k];
                    } else {
                        qc_grid = cldfra_bl1d[k]*qc_bl1d[k];
                    }
                    envm_sqc[k-kts] = envm_sqc[k-kts] + std::max(upa[k-1][i]*0.5_rt*(qcn + upqc[k-1][i]) - qc_grid, 0.0_rt)*detrate*aratio*std::min(dzp, dzpmax);
                }
                envm_u[k] = envm_u[k] + (0.5_rt*(un + upu[k-1][i]) - u[k])*detrateuv*aratio*std::min(dzp, dzpmax);
                envm_v[k] = envm_v[k] + (0.5_rt*(vn + upv[k-1][i]) - v[k])*detrateuv*aratio*std::min(dzp, dzpmax);
                if (wn > 0._rt) {
                    upw[k][i] = wn;
                    upthv[k][i] = thvn;
                    upthl[k][i] = thln;
                    upqt[k][i] = qtn;
                    upqc[k][i] = qcn;
                    upu[k][i] = un;
                    upv[k][i] = vn;
                    upqke[k][i] = qken;
                    upqnc[k][i] = qncn;
                    upqni[k][i] = qnin;
                    upqnwfa[k][i] = qnwfan;
                    upqnifa[k][i] = qnifan;
                    upqnbca[k][i] = qnbcan;
                    upa[k][i] = upa[k-1][i];
                    if (mix_chem) {
                        for (int ic = 0; ic < nchem; ic++) {
                            upchem[k][i][ic] = chemn[ic];
                        }
                    }
                    ktop = std::max(ktop, k);
                } else {
                    break;
                }
            }
        if (debug_mf == 1) {
            bool print_mf=false;
            for (int ii=kts;ii>=kte;ii++){
            if (upw[ii][i] > 10.0_rt || upa[ii][i] < 0.0_rt || upa[ii][i] > atot || nup2 > 10)
            {
               print_mf=true;
            }
            }
            if (print_mf)
            {
                std::cout << "flq: " << flq << " fltv: " << fltv << " nup2: " << nup2 << std::endl;
                std::cout << "pblh: " << pblh << " wstar: " << wstar << " ktop: " << ktop << std::endl;
                std::cout << "sigmaw: " << sigmaw << " sigmath: " << sigmath << " sigmaqt: " << sigmaqt << std::endl;
                std::cout << "u: " << u << std::endl;
                std::cout << "v: " << v << std::endl;
                std::cout << "thl: " << thl << std::endl;
                for(int ii=kts;ii>=kte;ii++) std::cout << "upa: " << upa[ii][i] ;
                std::cout<< std::endl;
                for(int ii=kts;ii>=kte;ii++) std::cout << "upw: " << upw[ii][i];
                std::cout<< std::endl;
                for(int ii=kts;ii>=kte;ii++) std::cout << "upthl: " << upthl[ii][i];
                std::cout<< std::endl;
                for(int ii=kts;ii>=kte;ii++) std::cout << "upqt: " << upqt[ii][i];
                std::cout<< std::endl;
                for(int ii=kts;ii>=kte;ii++) std::cout << "ent: " << ent[ii][i];
                std::cout<< std::endl;
                }
            }
         }
    } else {
        nup2 = 0;
    }
    ktop = std::min(ktop, kte-1);
    if (ktop == 0) {
        ztop = 0.0;
    } else {
        ztop = zw[ktop];
    }
  if (nup2 > 0) {
        for (int i = 0; i < nup; i++) {
            for (int k = kts; k <= kte-1; k++) {
                s_aw[k+1] += rhoz[k]*upa[k][i]*upw[k][i]*psig_w;
                s_awthl[k+1] += rhoz[k]*upa[k][i]*upw[k][i]*upthl[k][i]*psig_w;
                s_awqt[k+1] += rhoz[k]*upa[k][i]*upw[k][i]*upqt[k][i]*psig_w;
                qc_plume = upqc[k][i];
                s_awqc[k+1] += rhoz[k]*upa[k][i]*upw[k][i]*qc_plume*psig_w;
                s_awqv[k+1] = s_awqt[k+1] - s_awqc[k+1];
            }
        }
        if (momentum_opt > 0) {
            for (int i = 0; i < nup; i++) {
                for (int k = kts; k <= kte-1; k++) {
                    s_awu[k+1] += rhoz[k]*upa[k][i]*upw[k][i]*upu[k][i]*psig_w;
                    s_awv[k+1] += rhoz[k]*upa[k][i]*upw[k][i]*upv[k][i]*psig_w;
                }
            }
        }
        if (tke_opt > 0) {
            for (int i = 0; i < nup; i++) {
                for (int k = kts; k <= kte-1; k++) {
                    s_awqke[k+1] += rhoz[k]*upa[k][i]*upw[k][i]*upqke[k][i]*psig_w;
                }
            }
        }
        if (mix_chem) {
            for (int k = kts; k <= kte; k++) {
                for (int i = 0; i < nup; i++) {
                    for (int ic = 0; ic < nchem; ic++) {
                        s_awchem[k+1][ic] += rhoz[k]*upa[k][i]*upw[k][i]*upchem[k][i][ic]*psig_w;
                    }
                }
            }
        }
        if (scalar_opt > 0) {
            for (int k = kts; k <= kte; k++) {
                for (int i = 0; i < nup; i++) {
                    s_awqnc[k+1] += rhoz[k]*upa[k][i]*upw[k][i]*upqnc[k][i]*psig_w;
                    s_awqni[k+1] += rhoz[k]*upa[k][i]*upw[k][i]*upqni[k][i]*psig_w;
                    s_awqnwfa[k+1] += rhoz[k]*upa[k][i]*upw[k][i]*upqnwfa[k][i]*psig_w;
                    s_awqnifa[k+1] += rhoz[k]*upa[k][i]*upw[k][i]*upqnifa[k][i]*psig_w;
                    s_awqnbca[k+1] += rhoz[k]*upa[k][i]*upw[k][i]*upqnbca[k][i]*psig_w;
                }
            }
        }
    }
    if (s_aw[kts+1] != 0.0_rt) {
        dzi[kts] = 0.5_rt*(dz[kts] + dz[kts+1]);
        flx1 = std::max(s_aw[kts+1]*(th[kts] - th[kts+1])/dzi[kts], 1.0e-5_rt);
    } else {
        flx1 = 0.0;
    }
    adjustment = 1.0;
    flt2=std::max(flt,0.0_rt);
    if (flx1 > fluxportion*flt2/dz[kts] && flx1 > 0.0_rt) {
        adjustment = fluxportion*flt2/dz[kts]/flx1;
        for (int k = kts+1; k <= kte; k++) {
            s_aw[k] *= adjustment;
            s_awthl[k] *= adjustment;
            s_awqt[k] *= adjustment;
            s_awqc[k] *= adjustment;
            s_awqv[k] = s_awqt[k] - s_awqc[k];
        }
        if (momentum_opt > 0) {
            for (int k = kts+1; k <= kte; k++) {
                s_awu[k] *= adjustment;
                s_awv[k] *= adjustment;
            }
        }
        if (tke_opt > 0) {
            for (int k = kts+1; k <= kte; k++) {
                s_awqke[k] *= adjustment;
            }
        }
        if (mix_chem) {
            for (int k = kts+1; k <= kte; k++) {
                for (int ic = 0; ic < nchem; ic++) {
                    s_awchem[k][ic] *= adjustment;
                }
            }
        }
        for (int k = kts; k <= kte-1; k++) {
            *upa[k] *= adjustment;
        }
    }
    for (int k = kts; k <= kte-1; k++) {
        for (int i = 1; i<=nup; i++){
            edmf_a[k] += *upa[k,i];
            edmf_w[k] += (*upa[k,i])*(*upw[k,i]);
            edmf_qt[k] += (*upa[k,i])*(*upqt[k,i]);
            edmf_thl[k] += (*upa[k,i])*(*upthl[k,i]);
            edmf_ent[k] += (*upa[k,i])*(*ent[k,i]);
            edmf_qc[k] += (* upa[k,i]) * (*upqc[k,i]);
        }
    }
    for (int k = kts; k <= kte-1; k++) {
        if (edmf_a[k] > 0.0_rt) {
            edmf_w[k] /= edmf_a[k];
            edmf_qt[k] /= edmf_a[k];
            edmf_thl[k] /= edmf_a[k];
            edmf_ent[k] /= edmf_a[k];
            edmf_qc[k] /= edmf_a[k];
            edmf_a[k] *= psig_w;
            if (edmf_a[k]*edmf_w[k] > maxmf) {
                maxmf = edmf_a[k]*edmf_w[k];
            }
        }
    }
    if (mix_chem) {
        for (int k = kts; k <= kte-1; k++) {
            for (int i = 1; i<=nup;i++){
                for (int ic = 0; ic < nchem; ic++) {
                    edmf_chem[k][ic] += rhoz[k]*upa[k][i]*upchem[k][i][ic];
                }
            }
        }
        for (int k = kts; k <= kte-1; k++) {
            if (edmf_a[k] > 0.0_rt) {
                for (int ic = 0; ic < nchem; ic++) {
                    edmf_chem[k][ic] /= edmf_a[k];
                }
            }
        }
    }
    if (ktop > 0) {
        Real maxqc=0;
        for (int ii=0;ii > ktop; ii++){
                if (edmf_qc[ii]>maxqc){
                    maxqc = edmf_qc[ii];
                }
        }
        if (maxqc < 1.0e-8) {
            maxmf = -1.0*maxmf;
        }
    }
    if (edmf_w[0] > 4.0_rt) {
        std::cout << "flq: " << flq << " fltv: " << fltv2 << std::endl;
        std::cout << "pblh: " << pblh << " wstar: " << wstar << std::endl;
        std::cout << "sigmaw: " << sigmaw << " sigmath: " << sigmath << " sigmaqt: " << sigmaqt << std::endl;
        std::cout << "edmf_a: ";
        for (int i = 0; i < 14; i++) {
            std::cout << edmf_a[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "edmf_w: ";
        for (int i = 0; i < 14; i++) {
            std::cout << edmf_w[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "edmf_qt: ";
        for (int i = 0; i < 14; i++) {
            std::cout << edmf_qt[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "edmf_thl: ";
        for (int i = 0; i < 14; i++) {
            std::cout << edmf_thl[i] << " ";
        }
        std::cout << std::endl;
    }
}



//
// ==================================================================
//     SUBROUTINE  mym_turbulence:
//
//     Input variables:    see subroutine mym_initialize
//       closure        : closure level (2.5, 2.6, or 3.0_rt)
//
//     # ql, vt, vq, qke, tsq, qsq and cov are changed to input variables.
//
//     Output variables:   see subroutine mym_initialize
//       dfm(nx,nz,ny) : Diffusivity coefficient for momentum,
//                         divided by dz (not dz*h(i,j))            (m/s)
//       dfh(nx,nz,ny) : Diffusivity coefficient for heat,
//                         divided by dz (not dz*h(i,j))            (m/s)
//       dfq(nx,nz,ny) : Diffusivity coefficient for q^2,
//                         divided by dz (not dz*h(i,j))            (m/s)
//       tcd(nx,nz,ny)   : Countergradient diffusion term for Theta_l
//                                                                  (K/s)
//       qcd(nx,nz,ny)   : Countergradient diffusion term for Q_w
//                                                             (kg/kg s)
//       pd?(nx,nz,ny) : Half of the production terms
//
//       Only tcd and qcd are defined at the center of the grid boxes
//
//     # DO NOT forget that tcd and qcd are added on the right-hand side
//       of the equations for Theta_l and Q_w, respectively.
//
//     Work arrays:        see subroutine mym_initialize and level2
//
//     # dtl, dqw, dtv, gm and gh are allowed to share storage units with
//       dfm, dfh, dfq, tcd and qcd, respectively, for saving memory.
//
//\ingroup gsd_mynn_edmf
// This subroutine calculates the vertical diffusivity coefficients and the
// production terms for the turbulent quantities.
//\section gen_mym_turbulence GSD mym_turbulence General Algorithm
// Two subroutines mym_level2() and mym_length() are called within this
//subrouine to collect variable to carry out successive calculations:
// - mym_level2() calculates the level 2 nondimensional wind shear \f$G_M\f$
// and vertical temperature gradient \f$G_H\f$ as well as the level 2 stability
// functions \f$S_h\f$ and \f$S_m\f$.
// - mym_length() calculates the mixing lengths.
// - The stability criteria from Helfand and Labraga (1989) are applied.
// - The stability functions for level 2.5_rt or level 3.0_rt are calculated.
// - If level 3.0_rt is used, counter-gradient terms are calculated.
// - Production terms of TKE,\f$\theta^{'2}\f$,\f$q^{'2}\f$, and \f$\theta^{'}q^{'}\f$
// are calculated.
// - Eddy diffusivity \f$K_h\f$ and eddy viscosity \f$K_m\f$ are calculated.
// - TKE budget terms are calculated (if the namelist parameter \p tke_budget
// is set to True)
void mym_turbulence_cc(
    int& kts, int& kte,
    Real& xland, Real& closure,
    Real* dz, Real& dx, Real* zw,
    Real* u, Real* v, Real* thl, Real* thetav,
    Real* ql, Real* qw,
    Real* qke, Real* tsq, Real* qsq, Real* cov,
    Real* vt, Real* vq,
    Real& rmo, Real& flt, Real& fltv, Real& flq,
    Real& zi, Real* theta,
    /* begin intent(out) */
    Real* sh, Real* sm, Real* el,
    Real* dfm, Real* dfh, Real* dfq,
    Real* tcd, Real* qcd, Real* pdk,
    Real* pdt, Real* pdq, Real* pdc,
    /* end intent(out) */
    /* outputs if tke_budget==1, intent(inout) */
    Real* qWT1D, Real* qSHEAR1D, Real* qBUOY1D, Real* qDISS1D,
    /* end tke_budget */
    int& tke_budget, Real& Psig_bl, Real& Psig_shcu,
    Real* cldfra_bl1D, int& bl_mynn_mixlength,
    Real* edmf_w1, Real* edmf_a1,
    Real* TKEprodTD,
    int& spp_pbl, Real* rstoch_col,
    // additional params
    int& debug_code, Real& gtr, Real& tv0)
{
    Real q3sq_old, dlsq1, qWTP_old, qWTP_new;
    Real dudz, dvdz, dTdz, upwp, vpwp, Tpwp;
    Real e6c, dzk, afk, abk, vtt, vqq, cw25, clow, cupp, gamt, gamq, smd, gamv, elq, elh;
    Real cldavg;
    Real a2fac, duz, ri;
    Real auh, aum, adh, adm, aeh, aem, Req, Rsl, Rsl2, gmelq, sm20, sh20, sm25max, sh25max, sm25min, sh25min, sm_pbl, sh_pbl, zi2, wt, slht, wtpr;
    double q2sq, t2sq, r2sq, c2sq, elsq, gmel, ghel, q3sq, t3sq, r3sq, c3sq, dlsq, qdiv, e1, e2, e3, e4, enumc, eden, wden;
    Real Prnum, shb;
    const Real Prlimit = 5.0;
    // outputs from mym_level2:
    Real dtv[kte-kts];
    Real gm[kte-kts];
    Real gh[kte-kts];
    Real dqw[kte-kts];
    Real dtl[kte-kts];
    // output from mym_length:
    Real qkw[kte-kts];

    mym_level2_cc(kts, kte, dz, u, v,
                  thl, thetav, qw, ql,
                  vt, vq, dtl, dqw,
                  dtv, gm, gh, sm, sh,
                  // model constants
                  tv0, gtr);

    // calculate el and qkw
    mym_length_cc(kts, kte, xland,
                  dz, /*dx,*/ zw,
                  rmo, flt, fltv, flq,
                  vt, vq,
                  u, v, qke,
                  dtv,
                  el,
                  zi, theta,
                  qkw, Psig_bl, cldfra_bl1D,
                  bl_mynn_mixlength,
                  edmf_w1, edmf_a1,
                  // model constants
                  tv0, gtr);

    for (int k = kts + 1; k <= kte; k++) {
        dzk = 0.5_rt * (dz[k] + dz[k - 1]);
        afk = dz[k] / (dz[k] + dz[k - 1]);
        abk = 1.0_rt - afk;
        elsq = el[k] * el[k];
        q3sq = qkw[k] * qkw[k];
        q2sq = b1 * elsq * (sm[k] * gm[k] + sh[k] * gh[k]);

        sh20 = std::max(sh[k], 1e-5_rt);
        sm20 = std::max(sm[k], 1e-5_rt);
        sh[k] = std::max(sh[k], 1e-5_rt);
        //      printf("sh[k] %d %g\n",k,sh[k]);

        // Canuto/Kitamura mod
        duz = (u[k] - u[k - 1]) * (u[k] - u[k - 1]) + (v[k] - v[k - 1]) * (v[k] - v[k - 1]);
        duz = duz / (dzk * dzk);

        // ** Gradient Richardson number **
        ri = -gh[k] / std::max(duz, 1.0e-10_rt);
        if (ckmod == 1) {
            a2fac = 1.0_rt / (1.0_rt + std::max(ri, 0.0_rt));
        } else {
            a2fac = 1.0;
        }

        Prnum = std::min(0.76_rt + 4.0_rt * std::max(ri, 0.0_rt), Prlimit);

        gmel = gm[k] * elsq;
        ghel = gh[k] * elsq;

        if (debug_code) {
            if (sh[k] < 0.0_rt || sm[k] < 0.0_rt) {
                std::cout << "MYNN; mym_turbulence 2.0; sh=" << sh[k] << " k=" << k << std::endl;
                std::cout << " gm=" << gm[k] << " gh=" << gh[k] << " sm=" << sm[k] << std::endl;
                std::cout << " q2sq=" << q2sq << " q3sq=" << q3sq << " q3/q2=" << q3sq / q2sq << std::endl;
                std::cout << " qke=" << qke[k] << " el=" << el[k] << " ri=" << ri << std::endl;
                std::cout << " PBLH=" << zi << " u=" << u[k] << " v=" << v[k] << std::endl;
            }
        }

        // ** Since qkw is set to more than 0.0, q3sq > 0.0 **

        // new stability criteria in level 2.5 (as well as level 3) - little/no impact
        // ** Limitation on q, instead of L/q **
        dlsq = elsq;
        if (q3sq / dlsq < -gh[k]) q3sq = -dlsq * gh[k];

        if (q3sq < q2sq) {
            // Apply Helfand & Labraga mod
            qdiv = std::sqrt(q3sq / q2sq); // HL89: (1-alfa)

            // Use level 2.0 functions as in original MYNN
            sh[k] = sh[k] * qdiv;
            sm[k] = sm[k] * qdiv;
            //      printf("sh[k] %d %g %g\n",k,sh[k],qdiv);

            // Recalculate terms for later use
            e1 = q3sq - e1c * ghel * a2fac * qdiv * qdiv;
            e2 = q3sq - e2c * ghel * a2fac * qdiv * qdiv;
            e3 = e1 + e3c * ghel * a2fac * a2fac * qdiv * qdiv;
            e4 = e1 - e4c * ghel * a2fac * qdiv * qdiv;
            eden = e2 * e4 + e3 * e5c * gmel * qdiv * qdiv;
            eden = std::max(eden, 1.0e-20);
        } else {
            // JOE-Canuto/Kitamura mod
            e1 = q3sq - e1c * ghel * a2fac;
            e2 = q3sq - e2c * ghel * a2fac;
            e3 = e1 + e3c * ghel * a2fac * a2fac;
            e4 = e1 - e4c * ghel * a2fac;
            eden = e2 * e4 + e3 * e5c * gmel;
            eden = std::max(eden, 1.0e-20);
            qdiv = 1.0;

            // Use level 2.5 stability functions
            sm[k] = q3sq * a1 * (e3 - 3.0_rt * c1 * e4) / eden;
            sh[k] = q3sq * (a2 * a2fac) * (e2 + 3.0_rt * c1 * e5c * gmel) / eden;
            //      printf("sh[k] %d %g %g %g %g %g %g %g %g %g\n",k,sh[k],q3sq,a2,a2fac,e2,c1,e5c,gmel,eden);
        }

        // Impose broad limits on Sh and Sm
        gmelq = std::max(gmel / q3sq, 1e-8);
        sm25max = 4.0;
        sh25max = 4.0;
        sm25min = 0.0;
        sh25min = 0.0;

        if (debug_code) {
            if (sh[k] < sh25min || sm[k] < sm25min || sh[k] > sh25max || sm[k] > sm25max) {
                std::cout << "In mym_turbulence 2.5: k=" << k << std::endl;
                std::cout << " sm=" << sm[k] << " sh=" << sh[k] << std::endl;
                std::cout << " ri=" << ri << " Pr=" << sm[k] / std::max(sh[k], 1e-8_rt) << std::endl;
                std::cout << " gm=" << gm[k] << " gh=" << gh[k] << std::endl;
                std::cout << " q2sq=" << q2sq << " q3sq=" << q3sq << " q3/q2=" << q3sq / q2sq << std::endl;
                std::cout << " qke=" << qke[k] << " el=" << el[k] << std::endl;
                std::cout << " PBLH=" << zi << " u=" << u[k] << " v=" << v[k] << std::endl;
                std::cout << " SMnum=" << q3sq * a1 * (e3 - 3.0_rt * c1 * e4) << " SMdenom=" << eden << std::endl;
                std::cout << " SHnum=" << q3sq * (a2 * a2fac) * (e2 + 3.0_rt * c1 * e5c * gmel) << " SHdenom=" << eden << std::endl;
            }
        }

        // Enforce constraints for level 2.5 functions
        if (sh[k] > sh25max) sh[k] = sh25max;
        if (sh[k] < sh25min) sh[k] = sh25min;
        //      printf("sh[k] %d %g\n",k,sh[k]);

        shb = std::max(sh[k], 0.02_rt);
        sm[k] = std::min(sm[k], Prlimit * shb);

        /**** Level 3 : start ****/
        if (closure >= 3.0_rt) {
            t2sq = qdiv * b2 * elsq * sh[k] * dtl[k] * dtl[k];
            r2sq = qdiv * b2 * elsq * sh[k] * dqw[k] * dqw[k];
            c2sq = qdiv * b2 * elsq * sh[k] * dtl[k] * dqw[k];
            t3sq = std::max(tsq[k] * abk + tsq[k - 1] * afk, 0.0_rt);
            r3sq = std::max(qsq[k] * abk + qsq[k - 1] * afk, 0.0_rt);
            c3sq = cov[k] * abk + cov[k - 1] * afk;

            c3sq = std::copysign(std::min(std::abs(c3sq), std::sqrt(t3sq * r3sq)), c3sq);

            vtt = 1.0_rt + vt[k] * abk + vt[k - 1] * afk;
            vqq = tv0 + vq[k] * abk + vq[k - 1] * afk;

            t2sq = vtt * t2sq + vqq * c2sq;
            r2sq = vtt * c2sq + vqq * r2sq;
            c2sq = std::max(vtt * t2sq + vqq * r2sq, 0.0);
            t3sq = vtt * t3sq + vqq * c3sq;
            r3sq = vtt * c3sq + vqq * r3sq;
            c3sq = std::max(vtt * t3sq + vqq * r3sq, 0.0);

            cw25 = e1 * (e2 + 3.0_rt * c1 * e5c * gmel * qdiv * qdiv) / (3.0_rt * eden);

            // ** Limitation on q, instead of L/q **
            dlsq = elsq;
            if (q3sq / dlsq < -gh[k]) q3sq = -dlsq * gh[k];

            // ** Limitation on c3sq (0.12 =< cw =< 0.76) **
            auh = 27.0_rt * a1 * ((a2 * a2fac) * (a2 * a2fac)) * b2 * (gtr) * (gtr);
            aum = 54.0_rt * (a1 * a1) * (a2 * a2fac) * b2 * c1 * (gtr);
            adh = 9.0_rt * a1 * ((a2 * a2fac) * (a2 * a2fac)) * (12.0_rt * a1 + 3.0_rt * b2) * (gtr) * (gtr);
            adm = 18.0_rt * (a1 * a1) * (a2 * a2fac) * (b2 - 3.0_rt * (a2 * a2fac)) * (gtr);

            aeh = (9.0_rt * a1 * ((a2 * a2fac) * (a2 * a2fac)) * b1 + 9.0_rt * a1 * ((a2 * a2fac) * (a2 * a2fac)) * (12.0_rt * a1 + 3.0_rt * b2)) * (gtr);
            aem = 3.0_rt * a1 * (a2 * a2fac) * b1 * (3.0_rt * (a2 * a2fac) + 3.0_rt * b2 * c1 + (18.0_rt * a1 * c1 - b2)) + (18.0_rt) * (a1 * a1) * (a2 * a2fac) * (b2 - 3.0_rt * (a2 * a2fac));

            Req = -aeh / aem;
            Rsl = (auh + aum * Req) / (3.0_rt * adh + 3.0_rt * adm * Req);
            // For now, use default values, since tests showed little/no sensitivity
            Rsl = 0.12;
            Rsl2 = 1.0_rt - 2.0_rt * Rsl;

            // JOE-Canuto/Kitamura mod
            e2 = q3sq - e2c * ghel * a2fac * qdiv * qdiv;
            e3 = q3sq + e3c * ghel * a2fac * a2fac * qdiv * qdiv;
            e4 = q3sq - e4c * ghel * a2fac * qdiv * qdiv;
            eden = e2 * e4 + e3 * e5c * gmel * qdiv * qdiv;

            wden = cc3 * gtr * gtr * dlsq * dlsq / elsq * qdiv * qdiv * (e2 * e4c * a2fac - e3c * e5c * gmel * a2fac * a2fac * qdiv * qdiv);

            if (wden != 0.0_rt) {
                clow = q3sq * (0.12_rt - cw25) * eden / wden;
                cupp = q3sq * (0.76_rt - cw25) * eden / wden;
                if (wden > 0.0_rt) {
                    c3sq = std::min(std::max(c3sq, c2sq + clow), c2sq + cupp);
                } else {
                    c3sq = std::max(std::min(c3sq, c2sq + clow), c2sq + cupp);
                }
            }

            e1 = e2 + e5c * gmel * qdiv * qdiv;
            eden = std::max(eden, 1.0e-20);

            // JOE-Canuto/Kitamura mod
            e6c = 3.0_rt * (a2 * a2fac) * cc3 * gtr * dlsq / elsq;

            // ** for Gamma_theta **
            if (t2sq >= 0.0_rt) {
                enumc = std::max(qdiv * e6c * (t3sq - t2sq), 0.0);
            } else {
                enumc = std::min(qdiv * e6c * (t3sq - t2sq), 0.0);
            }
            gamt = -e1 * enumc / eden;

            // ** for Gamma_q **
            if (r2sq >= 0.0_rt) {
                enumc = std::max(qdiv * e6c * (r3sq - r2sq), 0.0);
            } else {
                enumc = std::min(qdiv * e6c * (r3sq - r2sq), 0.0);
            }
            gamq = -e1 * enumc / eden;

            // ** for Sm' and Sh'd(Theta_V)/dz **
            enumc = std::max(qdiv * e6c * (c3sq - c2sq), 0.0);

            // JOE-Canuto/Kitamura mod
            smd = dlsq * enumc * gtr / eden * qdiv * qdiv * (e3c * a2fac * a2fac + e4c * a2fac) * a1 / (a2 * a2fac);

            gamv = e1 * enumc * gtr / eden;
            sm[k] = sm[k] + smd;

            // ** For elh (see below), qdiv at Level 3 is reset to 1.0. **
            qdiv = 1.0;

            if (debug_code) {
                if (sh[k] < -0.3_rt || sm[k] < -0.3_rt || qke[k] < -0.1_rt || std::abs(smd) > 2.0_rt) {
                    std::cout << "MYNN; mym_turbulence3.0; sh=" << sh[k] << " k=" << k << std::endl;
                    std::cout << " gm=" << gm[k] << " gh=" << gh[k] << " sm=" << sm[k] << std::endl;
                    std::cout << " q2sq=" << q2sq << " q3sq=" << q3sq << " q3/q2=" << q3sq / q2sq << std::endl;
                    std::cout << " qke=" << qke[k] << " el=" << el[k] << " ri=" << ri << std::endl;
                    std::cout << " PBLH=" << zi << " u=" << u[k] << " v=" << v[k] << std::endl;
                }
            }

        /**** Level 3 : end ****/
        } else {
            // ** At Level 2.5, qdiv is not reset **
            gamt = 0.0;
            gamq = 0.0;
            gamv = 0.0;
        }

        // Add min background stability function (diffusivity) within model levels
        // with active plumes and clouds
        cldavg = 0.5_rt * (cldfra_bl1D[k - 1] + cldfra_bl1D[k]);
        if (edmf_a1[k] > 0.001 || cldavg > 0.02) {
            sm[k] = std::max(sm[k], 0.03_rt * std::min(10.0_rt * edmf_a1[k] * edmf_w1[k], 1.0_rt));
            sh[k] = std::max(sh[k], 0.03_rt * std::min(10.0_rt * edmf_a1[k] * edmf_w1[k], 1.0_rt));
            sm[k] = std::max(sm[k], 0.05_rt * std::min(cldavg, 1.0_rt));
            sh[k] = std::max(sh[k], 0.05_rt * std::min(cldavg, 1.0_rt));
            //      printf("sh[k] %d %g\n",k,sh[k]);
        }

        elq = el[k] * qkw[k];
        elh = elq * qdiv;

        // Production of TKE (pdk), T-variance (pdt),
        // q-variance (pdq), and covariance (pdc)
        pdk[k] = elq * (sm[k] * gm[k] + sh[k] * gh[k] + gamv) + 0.5_rt * TKEprodTD[k];
        pdt[k] = elh * (sh[k] * dtl[k] + gamt) * dtl[k];
        pdq[k] = elh * (sh[k] * dqw[k] + gamq) * dqw[k];
        pdc[k] = elh * (sh[k] * dtl[k] + gamt) * dqw[k] * 0.5_rt + elh * (sh[k] * dqw[k] + gamq) * dtl[k] * 0.5;

        // Countergradient terms
        tcd[k] = elq * gamt;
        qcd[k] = elq * gamq;

        // Eddy Diffusivity/Viscosity divided by dz
        dfm[k] = elq * sm[k] / dzk;
        dfh[k] = elq * sh[k] / dzk;
        dfq[k] = dfm[k];

        if (tke_budget == 1) {
            qSHEAR1D[k] = elq * sm[k] * gm[k];
            qBUOY1D[k] = elq * (sh[k] * gh[k] + gamv) + 0.5_rt * TKEprodTD[k];
        }
    }

    dfm[kts] = 0.0;
    dfh[kts] = 0.0;
    dfq[kts] = 0.0;
    tcd[kts] = 0.0;
    qcd[kts] = 0.0;
    tcd[kte] = 0.0;
    qcd[kte] = 0.0;
    for (int k = kts; k <= kte - 1; k++) {
        dzk = dz[k];
        tcd[k] = (tcd[k + 1] - tcd[k]) / dzk;
        qcd[k] = (qcd[k + 1] - qcd[k]) / dzk;
    }
    if (spp_pbl == 1) {
        for (int k = kts; k <= kte; k++) {
            dfm[k] = dfm[k] + dfm[k] * rstoch_col[k] * 1.5_rt * std::max(Real(exp(-std::max(zw[k] - 8000.0_rt, 0.0_rt) / 2000.0_rt)), 0.001_rt);
            dfh[k] = dfh[k] + dfh[k] * rstoch_col[k] * 1.5_rt * std::max(Real(exp(-std::max(zw[k] - 8000.0_rt, 0.0_rt) / 2000.0_rt)), 0.001_rt);
        }
    }
}


//!=======================================================================
//     SUBROUTINE  mym_initialize:
//
//     Input variables:
//       iniflag         : <>0; turbulent quantities will be initialized
//                         = 0; turbulent quantities have been already
//                              given, i.e., they will not be initialized
//       nx, nz          : Dimension sizes of the
//                         x and z directions, respectively
//       tref            : Reference temperature                      (K)
//       dz(nz)          : Vertical grid spacings                     (m)
//                         # dz(nz)=dz(nz-1)
//       zw(nz+1)        : Heights of the walls of the grid boxes     (m)
//                         # zw(1)=0.0_rt and zw(k)=zw(k-1)+dz(k-1)
//       exner(nx,nz)    : Exner function at zw*h+zg             (J/kg K)
//                         defined by c_p*( p_basic/1000hPa )^kappa
//                         This is usually computed by integrating
//                         d(pi0)/dz = -h*g/tref.
//       rmo(nx)         : Inverse of the Obukhov length         (m^(-1))
//       flt, flq(nx)    : Turbulent fluxes of potential temperature and
//                         total water, respectively:
//                                    flt=-u_*Theta_*             (K m/s)
//                                    flq=-u_*qw_*            (kg/kg m/s)
//       ust(nx)         : Friction velocity                        (m/s)
//       pmz(nx)         : phi_m-zeta at z1*h+z0, where z1 (=0.5*dz(1))
//                         is the first grid point above the surafce, z0
//                         the roughness length and zeta=(z1*h+z0)*rmo
//       phh(nx)         : phi_h at z1*h+z0
//       u, v(nx,nz)     : Components of the horizontal wind        (m/s)
//       thl(nx,nz)      : Liquid water potential temperature
//                                                                    (K)
//       qw(nx,nz)       : Total water content Q_w                (kg/kg)
//
//     Output variables:
//       ql(nx,nz)       : Liquid water content                   (kg/kg)
//       vt, vq(nx,nz)   : Functions for computing the buoyancy flux
//       qke(nx,nz)      : Twice the turbulent kinetic energy q^2
//                                                              (m^2/s^2)
//       tsq(nx,nz)      : Variance of Theta_l                      (K^2)
//       qsq(nx,nz)      : Variance of Q_w
//       cov(nx,nz)      : Covariance of Theta_l and Q_w              (K)
//       el(nx,nz)       : Master length scale L                      (m)
//                         defined on the walls of the grid boxes
//
//     Work arrays:        see subroutine mym_level2
//       pd?(nx,nz,ny) : Half of the production terms at Level 2
//                         defined on the walls of the grid boxes
//       qkw(nx,nz,ny) : q on the walls of the grid boxes         (m/s)
//
//     # As to dtl, ...gh, see subroutine mym_turbulence.
//
//-------------------------------------------------------------------

//>\ingroup gsd_mynn_edmf
// This subroutine initializes the mixing length, TKE, \f$\theta^{'2}\f$,
// \f$q^{'2}\f$, and \f$\theta^{'}q^{'}\f$.
//\section gen_mym_ini GSD MYNN-EDMF mym_initialize General Algorithm
//> @{
void mym_initialize_cc(
    const int &kts, const int &kte, const Real &xland,
    Real *dz, Real &dx, Real *zw,
    Real *u, Real *v, Real *thl, Real *qw,
    const Real &zi, Real *theta, Real *thetav, Real *sh, Real *sm,
    const Real& ust, const Real &rmo,
    Real* el, Real *qke, // intent(inout)
    Real* tsq, Real* qsq, Real* cov, // intent(out)
    const Real& Psig_bl, Real *cldfra_bl1D,
    int &bl_mynn_mixlength,
    Real *edmf_w1, Real *edmf_a1,
    int &INITIALIZE_QKE,
    int &spp_pbl, Real *rstoch_col,
    // model constants
    const Real& karman, const Real& tv0, const Real& gtr)
{
    Real phm, vkz, elq, elv, b1l, b2l, pmz = 1.0, phh = 1.0, flt = 0.0, fltv = 0.0, flq = 0.0, tmpq;
    int k, l, lmax;
    Real ql[kte-kts];
    Real vt[kte-kts];
    Real vq[kte-kts];
    Real pdk[kte-kts], pdt[kte-kts], pdq[kte-kts],pdc[kte-kts],dtl[kte-kts],dqw[kte-kts],dtv[kte-kts],gm[kte-kts],gh[kte-kts],qkw[kte-kts];

    // At first ql, vt and vq are set to zero.
    for (k = kts; k <= kte; k++) {
        ql[k-kts] = 0.0;
        vt[k-kts] = 0.0;
        vq[k-kts] = 0.0;
    }
    // Call mym_level2() to calculate the stability functions at level 2.
    mym_level2_cc(kts, kte, dz, u, v, thl, thetav, qw, ql, vt, vq, dtl, dqw, dtv, gm, gh, sm, sh, tv0, gtr);

    // Preliminary setting
    el[kts] = 0.0;

    if (INITIALIZE_QKE==1) {
        // WRF has (b1*pmz)**(2.0/3.0), not the `two_third` constant
        qke[kts] = 1.5_rt * std::pow(ust, 2.0_rt) * std::cbrt((b1*pmz) * (b1*pmz));
        for (k = kts + 1; k <= kte; k++) {
            qke[k] = qke[kts] * std::max((ust * 700.0_rt - zw[k]) / (std::max(ust, 0.01_rt) * 700.0_rt), 0.01_rt);
        }
    }

    phm = phh * b2 / std::cbrt(b1 * pmz);
    tsq[kts] = phm * ((flt / ust) * (flt / ust));
    qsq[kts] = phm * ((flq / ust) * (flq / ust));
    cov[kts] = phm * (flt / ust) * (flq / ust);
    for (k = kts + 1; k <= kte; k++) {
        vkz = karman * zw[k];
        el[k] = vkz / (1.0_rt + vkz / 100.0_rt);
        tsq[k] = 0.0;
        qsq[k] = 0.0;
        cov[k] = 0.0;
    }

    // Initialization with an iterative manner
    lmax = 5;
    for (l = 1; l <= lmax; l++) {
        // Call mym_length() to calculate the master length scale.
        mym_length_cc(kts, kte, xland, dz, zw, rmo, flt, fltv, flq, vt, vq, u, v, qke, dtv, el, zi, theta, qkw, Psig_bl, cldfra_bl1D, bl_mynn_mixlength, edmf_w1, edmf_a1, tv0, gtr);

        for (k = kts + 1; k <= kte; k++) {
            elq = el[k] * qkw[k];
            pdk[k] = elq * (sm[k] * gm[k] + sh[k] * gh[k]);
            pdt[k] = elq * sh[k] * (dtl[k] * dtl[k]);
            pdq[k] = elq * sh[k] * (dqw[k] * dqw[k]);
            pdc[k] = elq * sh[k] * dtl[k] * dqw[k];
        }

        // ** Strictly, vkz*h[i,j] -> karman*( 0.5*dz[0]*h[i,j]+z0 ) **
        vkz = karman * 0.5_rt * dz[kts];
        elv = 0.5_rt * (el[kts + 1] + el[kts]) / vkz;
        if (INITIALIZE_QKE==1) {
            // WRF has (b1*pmz*elv)**(2.0/3.0), not the `two_third` constant
            qke[kts] = 1.0_rt * std::max(ust, 0.02_rt) * std::max(ust, 0.02_rt) * std::cbrt((b1 * pmz * elv) * (b1 * pmz * elv));
        }

        phm = phh * b2 / std::cbrt(b1 * pmz / std::pow(elv, 2.0_rt));
        tsq[kts] = phm * ((flt / ust) * (flt / ust));
        qsq[kts] = phm * ((flq / ust) * (flq / ust));
        cov[kts] = phm * (flt / ust) * (flq / ust);

        for (k = kts + 1; k <= kte - 1; k++) {
            b1l = b1 * 0.25_rt * (el[k + 1] + el[k]);
            tmpq = std::min(std::max(b1l * (pdk[k + 1] + pdk[k]), qkemin), 125.0_rt);
            if (INITIALIZE_QKE==1) {
                // WRF has tmpq**two_thirds, where two_thirds is a constant (kind_phys)
                qke[k] = std::pow(tmpq, 2.0_rt/3.0_rt);
            }

            if (qke[k] <= 0.0_rt) {
                b2l = 0.0;
            } else {
                b2l = b2 * (b1l / b1) / std::sqrt(qke[k]);
            }
            tsq[k] = b2l * (pdt[k + 1] + pdt[k]);
            qsq[k] = b2l * (pdq[k + 1] + pdq[k]);
            cov[k] = b2l * (pdc[k + 1] + pdc[k]);
        }
    }

    if (INITIALIZE_QKE==1) {
        qke[kts] = 0.5_rt * (qke[kts] + qke[kts + 1]);
        qke[kte] = qke[kte - 1];
    }

    tsq[kte] = tsq[kte - 1];
    qsq[kte] = qsq[kte - 1];
    cov[kte] = cov[kte - 1];
}




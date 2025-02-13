#include <algorithm> 
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <limits>

extern "C" void mynn_tendencies_cc(const int& kts,const int& kte, const float & delt,
				   /*in*/ const float* dz,
				   /*in*/ const float* rho,
				   /*in*/ const float* u, const float* v, const float* th, const float* tk, const float* qv,
				   /*in*/ const float* qc, const float* qi, const float* qs, const float* qni, const float* qnc,
				   /*in*/ const float* psfc,const float* p, const float* exner,
				   /*inout*/ float* thl, float* sqv, float* sqc, float* sqi,
				   /*inout*/ float* sqs, float* sqw, float* qnwfa, float* qnifa, float* qnbca, float* ozone,
				   /*in*/ float* ust, const float & flt,const float & flq,const float & flqv,const float & flqc,const float & wspd, const float & uoce, const float & voce,
				   /*in*/ const float* tcd,const float* qcd,
				   /*inout*/ float* dfm, float* dfh,
				   /*inout*/ float* du, float* dv, float* dth,
				   /*inout*/ float* dqv, float* dqc, float* dqi,
				   /*inout*/ float* dqs, float* dqnc, float* dqni,
				   /*inout*/ float* dqnwfa, float* dqnifa, float* dqnbca,
				   /*inout*/ float* dozone,
				   /*in*/ const float* diss_heat,
			/*in*/ const float* s_aw, const float* s_awthl, const float* s_awqt, const float* s_awqv, const float* s_awqc, const float* s_awu, const float* s_awv, const float* s_awqni, const float* s_awqnc, const float* s_awqnwfa, const float* s_awqnifa, const float* s_awqnbca, const float* sd_aw, const float* sd_awthl, const float* sd_awqt, const float* sd_awqv, const float* sd_awqc, const float* sd_awu, const float* sd_awv,
				   const float* sub_thl,const float* sub_sqv,const float* sub_u,const float* sub_v,
				   const float* det_thl,const float* det_sqv,const float* det_sqc,const float* det_u,const float* det_v,
				   /*logical turned into int */const int& flag_qc, const int& flag_qi, const int& flag_qnc, const int& flag_qni, const int& flag_qs, const int& flag_qnwfa, const int& flag_qnifa, const int& flag_qnbca, const int& flag_ozone,
				   const int & bl_mynn_cloudmix, const int & bl_mynn_mixqt, int & bl_mynn_edmf_mom,  const int & bl_mynn_mixscalars, /* new */const int& debug_code,const float& r_d,const float& p608,const float& ep_2,const float& ep_3,const float& tv0,const float& xlv,const float& xlvcp,const float & xlscp);
/*
skipping bl_mynn_mixscalars, flag_qni, flag_qc, flag_qs, flag_qnc, flag_qnwfa, flag_qnifa, flag_qnbca, flag_ozone // adding
skipping th, qc, qi, qs, qni, qnc //adding
skipping exner, dfq, tsq, qsq, cov, cldfra_bl1d
skipping sqi, sqs, qnwfa, qnifa, qnbca, ozone,dqv, dqc, dqi, dqs, dqni, dqnc, dqnwfa, dqnifa, dqnbca, dozone //adding some
skipping wsp, wsp2, tk2, th2
skipping problem, kproblem
looks like constants from common that are passed in: float r_d, float p608, float ep_2,float ep_3,float tv0,float xlv,float xlvcp
should these also be intent(in) and therefore const &?
 */
extern "C" void mym_predict_cc(int& kts, int& kte, float& closure, float& delt, float* dz, float* ust, float& flt, float& flq, float& pmz, float& phh, float* el, float* dfq, float* rho, float* pdk, float* pdt, float* pdq, float* pdc, float* qke, float* tsq, float* qsq, float* cov, float* s_aw, float* s_awqke, int& bl_mynn_edmf_tke, float* qwt1d, float* qdiss1d, int& tke_budget, float& xlvcp, float& xlscp, float& karman);

extern "C" void mynn_mix_chem_cc(int kts, int kte, int i,float delt, float* dz, float pblh, int nchem, int kdvel, int ndvel,float** chem1, float* vd1, float* rho,float flt, float* tcd, float* qcd, float* dfh,float* s_aw, float** s_awchem, float emis_ant_no, float frp, int rrfs_sd, int enh_mix); 

extern "C" void moisture_check_cc(int kte, float delt, float* dp, float* exner,float* qv, float* qc, float* qi, float* qs, float* th,float* dqv, float* dqc, float* dqi, float* dqs, float* dth,float dqv2, float xlvcp, float xlscp); 

extern "C" void mym_condensation_cc(const int& kts,const int& kte, const float& dx, float* dz, float* zw, float& xland,float* thl, 
		float* qw, float* qv, float* qc, float* qi, float* qs,float* p, float* exner, 
		float* tsq, float* qsq, float* cov,float* sh, float* el, int& bl_mynn_cloudpdf,
                float* qc_bl1d, float* qi_bl1d, float* cldfra_bl1d,float & pblh1, float & hfx1,
			 float* vt, float* vq, float* th, float* sgm, float* rmo,int &spp_pbl, float* rstoch_col, 
		float ep_2, float ep_3, float xlv, float r_d, float xlvcp, float p608, float tv0, float cpv, 
				    float r_v, float cice, float cliq, float cp, float xls, float rcp);


extern "C" void topdown_cloudrad_cc(int& kts, int& kte, const float* dz1, const float* zw, float& fltv, float& xland, int& kpbl, float& pblh, const float* sqc, const float* sqi, const float* sqw, const float* thl, const float* th1, const float* ex1, const float* p1, const float*  rho1, const float* thetav, const float* cldfra_bl1d, const float* rthraten, float& maxkhtopdown, float* khtopdown, float* tkeprodtd);

extern "C" void ddmf_jpl_cc(int& kts, int& kte, float& dt, const float* zw, const float* dz, const float* p,
              const float* u, const float* v, const float* th, const float* thl, const float* thv, 
	      const float* tk,const float* qt, const float* qv, const float* qc, const float* 
	      rho, const float* exner,float& ust, float& wthl, float& wqt, float& pblh, int& kpbl,
              float* edmf_a_dd, float* edmf_w_dd, float* edmf_qt_dd,
              float* edmf_thl_dd, float* edmf_ent_dd, float* edmf_qc_dd,
              float* sd_aw, float* sd_awthl, float* sd_awqt,
              float* sd_awqv, float* sd_awqc, float* sd_awu,
              float* sd_awv, float* sd_awqke,
              const float* qc_bl1d, const float* cldfra_bl1d,
              const float* rthraten, float& svp1, float& grav, float& onethird, float& p1000mb, 
              float& rcp, float& xlvcp, float& cp, float& rvovrd );

extern "C" void scale_aware_cc(float& dx, float& pbl1, float& psig_bl, float& psig_shcu); 

extern "C" void get_pblh_cc(int &kts, int &kte, float &zi, float *thetav1d, float *qke1d, float *zw1d, float *dz1d, float &landsea, int &kzi);

extern "C" void retrieve_exchange_coeffs_cc(int& kts, int& kte, float* dfm, float* dfh, const float* dz, float* k_m, float* k_h);

extern "C" void dmp_mf_cc(const int& kts, const int& kte, float& dt, float* zw, float* dz, float* p, float* rho, int& momentum_opt, int& tke_opt, int& scalar_opt, float* u, float* v, float* w, float* th, float* thl, float* thv, float* tk, float* qt, float* qv, float* qc, float* qke, float* qnc, float* qni, float* qnwfa, float* qnifa, float* qnbca, float& ust, float& flt, float& fltv, float& flq, float& flqv, float& pblh, int& kpbl, float& dx, float& landsea, float& ts, float* edmf_a, float* edmf_w, float* edmf_qt, float* edmf_thl, float* edmf_ent, float* edmf_qc, float* s_aw, float* s_awthl, float* s_awqt, float* s_awqv, float* s_awqc, float* s_awu, float* s_awv, float* s_awqke, float* s_awqnc, float* s_awqni, float* s_awqnwfa, float* s_awqnifa, float* s_awqnbca, int& nchem, float** chem1, float** s_awchem, bool& mix_chem, float* qc_bl1d, float* cldfra_bl1d, float* qc_bl1d_old, float* cldfra_bl1d_old, float& psig_shcu, float& maxwidth, int& ktop, float& maxmf, float& ztop, float* rstoch_col, float grav, float gtr, float p608);

extern "C" void mym_turbulence_cc(int& kts, int& kte, float& xland, float& closure, float* dz, float* dx, float* zw, float* u, float* v, float* thl, float* thetav, float* ql, float* qw, float* qke, float* tsq, float* qsq, float* cov, float* vt, float* vq, float& rmo, float& flt, float& fltv, float& flq, float& zi, float* theta, float* sh, float* sm, float* el, float* dfm, float* dfh, float* dfq, float* tcd, float* qcd, float* pdk, float* pdt, float* pdq, float* pdc, float* qWT1D, float* qSHEAR1D, float* qBUOY1D, float* qDISS1D, int& tke_budget, float& Psig_bl, float& Psig_shcu, float* cldfra_bl1D, int& bl_mynn_mixlength, float* edmf_w1, float* edmf_a1, float* TKEprodTD, int& spp_pbl, float* rstoch_col, int& debug_code, float& gtr, float& tv0);

extern "C" void mym_initialize_cc(const int &kts,const int &kte,const float &xland, float *dz, float &dx, float *zw, float *u, float *v, float *thl, float *qw,const float &zi, float *theta, float *thetav, float *sh, float *sm,const float& ust, const float &rmo, float* el, float *qke, float* tsq, float* qsa, float* cov, const float& Psig_bl, float *cldfra_bl1D, int &bl_mynn_mixlength, float *edmf_w1, float *edmf_a1, int &INITIALIZE_QKE, int &spp_pbl, float *rstoch_col,const float& karman,const float& tv0,const float& gtr);
//----------------------------------------contstants-------------------------------------------

// constants
const float no_threshold = 10.0;     // for anthropogenic sources
const float frp_threshold = 10.0;    // increased the frp threshold to enhance mixing over big fires
const float pblh_threshold = 100.0;

const float t0c = 273.15; // assuming t0c is 273.15
const float tice = 240.0; // assuming tice is 240 based on the comment

// assuming float corresponds to float precision
const float cphm_st = 5.0, cphm_unst = 16.0,
                 cphh_st = 5.0, cphh_unst = 16.0;
//    1.1800000667572021       0.1370676159858704       0.6645210385322571       0.5641299486160278
//       0.2710000276565552       0.6599999666213989      19.7362747192382812       1.9125051498413086       0.8616270422935486       2.5500068664550781       8.3544006347656250,

// closure constants
constexpr float pr = 0.7400000095367432,
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

  constexpr float cc2 = 1.0f - c2,
                 cc3 = 1.0f - c3,
                 e1c = 3.0f * a2 * b2 * cc3,
                 e2c = 9.0f * a1 * a2 * cc2,
                 e3c = 9.0f * a2 * a2 * cc2 * (1.0f - c5),
                 e4c = 12.0f * a1 * a2 * cc2,
                 e5c = 6.0f * a1 * a1;

// constants for min tke in elt integration (qmin), max z/l in els (zmax),
// and factor for eddy viscosity for tke (kq = sqfac*km):
constexpr float qmin = 0.0, zmax = 1.0, sqfac = 3.0;

constexpr float gpw = 5.0f / 3.0, qcgmin = 1e-8, qkemin = 1e-3;
constexpr float tliq = 269.0; // all hydrometeors are liquid when t > tliq

// constants for cloud pdf (mym_condensation)
constexpr float rr2 = 0.7071068, rrp = 0.3989423;

// use canuto/kitamura mod (remove ric and negative tke) (1:yes, 0:no)
constexpr float ckmod = 1.0;

// option to activate environmental subsidence in mass-flux scheme
constexpr bool env_subs = false;

//---------------------------------------------------------------------------------------------
float vsc = 1.0e-5;
float elt = 1.0e-5;

float esat_blend_cc(float t) {
    // constants for liquid
    const float j0 = .611583699e03;
    const float j1 = .444606896e02;
    const float j2 = .143177157e01;
    const float j3 = .264224321e-1;
    const float j4 = .299291081e-3;
    const float j5 = .203154182e-5;
    const float j6 = .702620698e-8;
    const float j7 = .379534310e-11;
    const float j8 = -.321582393e-13;

    // constants for ice
    const float k0 = .609868993e03;
    const float k1 = .499320233e02;
    const float k2 = .184672631e01;
    const float k3 = .402737184e-1;
    const float k4 = .565392987e-3;
    const float k5 = .521693933e-5;
    const float k6 = .307839583e-7;
    const float k7 = .105785160e-9;
    const float k8 = .161444444e-12;

    float xc = std::max(-80.0f, t - t0c);
    float esat_blend_cc;

    if (t >= (t0c - 6.0f)) {
        esat_blend_cc = j0 + xc * (j1 + xc * (j2 + xc * (j3 + xc * (j4 + xc * (j5 + xc * (j6 + xc * (j7 + xc * j8)))))));
    } else if (t <= tice) {
        esat_blend_cc = k0 + xc * (k1 + xc * (k2 + xc * (k3 + xc * (k4 + xc * (k5 + xc * (k6 + xc * (k7 + xc * k8)))))));
    } else {
        float esl = j0 + xc * (j1 + xc * (j2 + xc * (j3 + xc * (j4 + xc * (j5 + xc * (j6 + xc * (j7 + xc * j8)))))));
        float esi = k0 + xc * (k1 + xc * (k2 + xc * (k3 + xc * (k4 + xc * (k5 + xc * (k6 + xc * (k7 + xc * k8)))))));
        float chi = ((t0c - 6.0f) - t) / ((t0c - 6.0f) - tice);
        esat_blend_cc = (1.0f - chi) * esl + chi * esi;
    }

    return esat_blend_cc;
}


float qsat_blend_cc(float t, float p) {
    // constants for liquid
    const float j0 = .611583699e03;
    const float j1 = .444606896e02;
    const float j2 = .143177157e01;
    const float j3 = .264224321e-1;
    const float j4 = .299291081e-3;
    const float j5 = .203154182e-5;
    const float j6 = .702620698e-8;
    const float j7 = .379534310e-11;
    const float j8 = -.321582393e-13;

    // constants for ice
    const float k0 = .609868993e03;
    const float k1 = .499320233e02;
    const float k2 = .184672631e01;
    const float k3 = .402737184e-1;
    const float k4 = .565392987e-3;
    const float k5 = .521693933e-5;
    const float k6 = .307839583e-7;
    const float k7 = .105785160e-9;
    const float k8 = .161444444e-12;

    
    // temperature thresholds
    const float t0c = 273.15; // assuming 0 for t0c (temperature in celsius)
    const float tice = 240.00; // assuming -273.15f for tice (absolute zero, could be different)
    float xc = std::max(-80.0f, t - t0c);
    float qsat_blend_cc, esl, esi, rslf, rsif, chi;

    if (t >= (t0c - 6.0f)) {
        esl = j0 + xc * (j1 + xc * (j2 + xc * (j3 + xc * (j4 + xc * (j5 + xc * (j6 + xc * (j7 + xc * j8)))))));
        esl = std::min(esl, p * 0.15f);
        qsat_blend_cc = 0.622f * esl / std::max(p - esl, 1e-5f);
    } else if (t <= tice) {
        esi = k0 + xc * (k1 + xc * (k2 + xc * (k3 + xc * (k4 + xc * (k5 + xc * (k6 + xc * (k7 + xc * k8)))))));
        esi = std::min(esi, p * 0.15f);
        qsat_blend_cc = 0.622f * esi / std::max(p - esi, 1e-5f);
    } else {
        esl = j0 + xc * (j1 + xc * (j2 + xc * (j3 + xc * (j4 + xc * (j5 + xc * (j6 + xc * (j7 + xc * j8)))))));
        esl = std::min(esl, p * 0.15f);
        esi = k0 + xc * (k1 + xc * (k2 + xc * (k3 + xc * (k4 + xc * (k5 + xc * (k6 + xc * (k7 + xc * k8)))))));
        esi = std::min(esi, p * 0.15f);
        rslf = 0.622f * esl / std::max(p - esl, 1e-5f);
        rsif = 0.622f * esi / std::max(p - esi, 1e-5f);
        chi = ((t0c - 6.0f) - t) / ((t0c - 6.0f) - tice);
        qsat_blend_cc = (1.0f - chi) * rslf + chi * rsif;
    }
    return qsat_blend_cc;
}


float xl_blend_cc(float t,float xlv, float xls, float cpv, float cliq, float cice) {
    float xl_blend_cc, xlvt, xlst, chi;
    // t0c = 273.15, tice is set elsewhere
    if (t >= t0c) {
        xl_blend_cc = xlv + (cpv - cliq) * (t - t0c); // vaporization/condensation
    } else if (t <= tice) {
        xl_blend_cc = xls + (cpv - cice) * (t - t0c); // sublimation/deposition
    } else {
        xlvt = xlv + (cpv - cliq) * (t - t0c); // vaporization/condensation
        xlst = xls + (cpv - cice) * (t - t0c); // sublimation/deposition
        chi = (t0c - t) / (t0c - tice);
        xl_blend_cc = (1.f - chi) * xlvt + chi * xlst; // blended
    }
    return xl_blend_cc;
}

void condensation_edmf_cc(float qt, float thl, float p, float zagl, float& thv, float& qc, float p1000mb, float rcp, float xlvcp, float rvovrd) {
    const int niter = 50;
    const float diff = 1.e-6;
    float exn = std::pow((p / p1000mb), rcp);
    // qc is assumed to be initialized before calling this function
    for (int i = 0; i < niter; ++i) {
        float t = exn * thl + xlvcp * qc;
        float qs = qsat_blend_cc(t, p);
        float qcold = qc;
        qc = 0.5f * qc + 0.5f * std::max((qt - qs), 0.0f);
        if (std::abs(qc - qcold) < diff) break;
    }
    float t = exn * thl + xlvcp * qc;
    float qs = qsat_blend_cc(t, p);
    qc = std::max(qt - qs, 0.0f);
    // do not allow saturation below 1.0f m
    if (zagl < 100.0f) qc = 0.0;
    thv = (thl + xlvcp * qc) * (1.0f + qt * (rvovrd - 1.0f) - rvovrd * qc);
}

// function to solve system of linear equations on tridiagonal matrix n times n
// after peaceman and rachford, 1955
// a, b, c, d - are std::vectors of order n
// a, b, c - are coefficients on the lhs
// d - is initially rhs on the output becomes a solution std::vector
void tridiag_cc(int n, const float* a, const float* b, float* c, float* d) {
    float q[n];
    c[n-1] = 0.0;
    q[0] = -c[0] / b[0];
    d[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        float p = 1.0f / (b[i] + a[i] * q[i - 1]);
        q[i] = -c[i] * p;
        d[i] = (d[i] - a[i] * d[i - 1]) * p;
    }

    for (int i = n - 2; i >= 0; --i) {
        d[i] = d[i] + q[i] * d[i + 1];
    }
}

void tridiag2_cc(int n, float* a, float* b, float* c, float* d, float* x) {
    float cp[n];
    float dp[n];
    float m;

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
void tridiag3_cc(int kte, float* a, float* b, float* c, float* d, float* x) {
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
void boulac_length_cc(int& kts, int& kte, const float* zw, const float* dz, const float* qtke, const float* theta, float* lb1, float* lb2, float& gtr) {

//      dlu = the distance a parcel can be lifted upwards give a finite
//            amount of tke.
//      dld = the distance a parcel can be displaced downwards given a
//            finite amount of tke.
//      lb1 = the minimum of the length up and length down
//      lb2 = the average of the length up and length down
    int iz, izz, found;
    float dlu[kte-kts];
    float dld[kte-kts];
    const float lmax = 2000.0;
    float dzt, zup, beta, zup_inf, bbb, tl, zdo, zdo_sup, zzz;

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
                        
                        if (bbb != 0.0f) {
                            tl = (-beta * (theta[izz] - theta[iz]) + sqrt(std::max(0.0f, ((beta * (theta[izz] - theta[iz])) * (beta * (theta[izz] - theta[iz]))) + 2.0f * bbb * beta * (qtke[iz] - zup_inf)))) / bbb / beta;
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
                        
                        if (bbb != 0.0f) {
                            tl = (beta * (theta[izz] - theta[iz]) + sqrt(std::max(0.0f, ((beta * (theta[izz] - theta[iz])) * (beta * (theta[izz] - theta[iz]))) + 2.0f * bbb * beta * (qtke[iz] - zdo_sup)))) / bbb / beta;
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
        dlu[iz] = std::max(0.1f, std::min(dlu[iz], 1000.0f));
        dld[iz] = std::max(0.1f, std::min(dld[iz], 1000.0f));
        lb2[iz] = std::sqrt(dlu[iz] * dld[iz]);
        lb1[iz] = lb1[iz] / (1.0f + (lb1[iz] / lmax));
        lb2[iz] = lb2[iz] / (1.0f + (lb2[iz] / lmax));
        
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
// the level 2 stability funcitons \f$s_h\f$ and \f$s_m\f$.
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
void mym_level2_cc(int kts, int kte, float* dz, float* u, float* v,
                float* thl, float* thetav, float* qw, float* ql,
                float* vt, float* vq, float* dtl, float* dqw,
                float* dtv, float* gm, float* gh, float* sm, float* sh, 
		float tv0, float gtr) {
    float rfc, f1, f2, rf1, rf2, smc, shc, ri1, ri2, ri3, ri4, duz, dtz, dqz, vtt, vqq, dtq, dzk, afk, abk, ri, rf;
    float a2fac;

    rfc = g1 / (g1 + g2);
    f1 = b1 * (g1 - c1) + 3.0f * a2 * (1.0f - c2) * (1.0f - c5) + 2.0f * a1 * (3.0f - 2.0f * c2);
    f2 = b1 * (g1 + g2) - 3.0f * a1 * (1.0f - c2);
    rf1 = b1 * (g1 - c1) / f1;
    rf2 = b1 * g1 / f2;
    smc = a1 / a2 * f1 / f2;
    shc = 3.0f * a2 * (g1 + g2);
    //    printf("g1 %15.15g %15.15g %15.15g %15.15g %15.15g\n",g1,g2,rfc,rf1,rf2);
    //    printf("g2 %15.15g %15.15g %15.15g %15.15g %15.15g %15.5g %15.5g\n",g2,b2,b1,( 1.0f-c3 ) ,2.0f*a1,b1,( 3.0f-2.0f*c2 ));
    ri1 = 0.5f / smc;
    ri2 = rf1 * smc;
    ri3 = 4.0f * rf2 * smc - 2.0f * ri2;
    ri4 = ri2 * ri2;

    for (int k = kts + 1; k <= kte; ++k) {
        dzk = 0.5f * (dz[k] + dz[k - 1]);
        afk = dz[k] / (dz[k] + dz[k - 1]);
        abk = 1.0f - afk;
        duz = ((u[k] - u[k - 1]) * (u[k] - u[k-1])) + ((v[k] - v[k - 1]) * (v[k] - v[k - 1]));
        duz = duz / (dzk * dzk);
        dtz = (thl[k] - thl[k - 1]) / dzk;
        dqz = (qw[k] - qw[k - 1]) / dzk;

        vtt = 1.0f + vt[k] * abk + vt[k - 1] * afk; // beta-theta in nn09, eq. 39
        vqq = tv0 + vq[k] * abk + vq[k - 1] * afk; // beta-q
        dtq = vtt * dtz + vqq * dqz;
        // alternatively, use theta-v without the sgs clouds
        // dtq = (thetav[k] - thetav[k - 1]) / dzk;

        dtl[k] = dtz;
        dqw[k] = dqz;
        dtv[k] = dtq;

        gm[k] = duz;
        gh[k] = -dtq * gtr;

        // gradient richardson number
        ri = -gh[k] / std::max(duz, 1.0e-10f);
        // a2fac is needed for the canuto/kitamura mod
        if (ckmod == 1) {
            a2fac = 1.0f / (1.0f + std::max(ri, 0.0f));
        } else {
            a2fac = 1.0f;
        }
        rfc = g1 / (g1 + g2);
        f1 = b1 * (g1 - c1) + 3.0f * a2 * a2fac * (1.0f - c2) * (1.0f - c5) + 2.0f * a1 * (3.0f - 2.0f * c2);
        f2 = b1 * (g1 + g2) - 3.0f * a1 * (1.0f - c2);
        rf1 = b1 * (g1 - c1) / f1;
        rf2 = b1 * g1 / f2;
        smc = a1 / (a2 * a2fac) * f1 / f2;
        shc = 3.0f * (a2 * a2fac) * (g1 + g2);
        ri1 = 0.5f / smc;
        ri2 = rf1 * smc;
        ri3 = 4.0f * rf2 * smc - 2.0f * ri2;
        ri4 = ri2 * ri2;

        // flux richardson number
        rf = std::min(ri1 * (ri + ri2 - std::sqrt(ri * ri - ri3 * ri + ri4)), rfc);
	//	printf("rf    %15.15g %15.15g %15.15g %15.15g %15.15g %15.15g\n",rf2,ri2,ri3,smc,ri4);
        sh[k] = shc * (rfc - rf) / (1.0f - rf);
        sm[k] = smc * (rf1 - rf) / (rf2 - rf) * sh[k];
	//	printf("sm[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",sm[k],shc,rfc,rf,1.0f);
	//	printf("sh[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",sh[k],smc,rf1,rf,rf2);
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
//           sigmal levels (or interfaces beween the model layers).
//
//>\ingroup gsd_mynn_edmf
// this subroutine calculates the mixing lengths.
void mym_length_cc(int kts, int kte, float xland, float* dz, float* zw, float rmo, float flt, float fltv, float flq, float* vt, float* vq, float* u1, float* v1, float* qke, float* dtv, float* el, float zi, float* theta, float* qkw, float psig_bl, float* cldfra_bl1d, int bl_mynn_mixlength, float* edmf_w1, float* edmf_a1, float tv0, float gtr) {
    float cns, alp1, alp2, alp3, alp4, alp5, alp6;
    float minzi = 300.0;
    float maxdz = 750.0;
    float mindz = 300.0;
    float zslh = 100.0;
    float csl = 2.0;
    float afk, abk, zwk, zwk1, dzk, qdz, vflx, bv, tau_cloud, wstar, elb, els, elf, el_stab, el_mf, el_stab_mf, elb_mf, pblh_plus_ent, uonset, ugrid, wt_u, el_les;
    float qke_elb_min = 0.018f;
    const float grav = 9.8100004196166992, karman = 0.4000000059604645;
    const float twothirds = 0.6666666865348816, onethird = 0.3333333432674408;
    const float qmin = 0.0f;
    const float ctau = 1000.f;

    float qtke[kte+1];
    float thetaw[kte+1];
    float elblmin[kte+1];
    float elblavg[kte+1];
    float zi2, h1, h2, hs, elblmin0, elblavg0, cldavg;
    int i, j, k;

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
            zi2 = std::min(10000.0f, float(zw[kte-2]));  // originally integrated to model top, not just 10 km.
            h1 = std::max(0.3f * float(zi2), float(mindz));
            h1 = std::min(float(h1), float(maxdz));         // 1/2 transition layer depth
            h2 = h1 / 2.0;                                  // 1/4 transition layer depth

            qkw[kts] = std::sqrt(std::max(float(qke[kts]), qkemin));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0f - afk;
                qkw[k] = std::sqrt(std::max(float(qke[k] * abk + qke[k-1] * afk), qkemin));
            }

            elt = 1.0e-5;
            vsc = 1.0e-5;

            // ** Strictly, zwk*h[i,j] -> ( zwk*h[i,j]+z0 ) **
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi2 + h1) {
                dzk = 0.5f * (dz[k] + dz[k-1]);
                qdz = std::max(float(qkw[k] - qmin), 0.03f) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }

            elt = alp1 * elt / vsc;
            vflx = (vt[kts] + 1.0f) * flt + (vq[kts] + tv0) * flq;
            vsc = std::cbrt(gtr * elt * std::max(float(vflx), 0.0f));

            // ** Strictly, el[i,k=0] is not zero. **
            el[kts] = 0.0;
            zwk1 = zw[kts+1];

            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];

                // ** Length scale limited by the buoyancy effect **
                if (dtv[k] > 0.0f) {
                    bv = std::sqrt(gtr * dtv[k]);
                    elb = alp2 * qkw[k] / bv * (1.0f + alp3 / alp2 * std::sqrt(vsc / (bv * elt)));
                    elf = alp2 * qkw[k] / bv;
                } else {
                    elb = 1.0e10;
                    elf = elb;
                }

                // ** Length scale in the surface layer **
                if (rmo > 0.0f) {
                    els = karman * zwk / (1.0f + cns * std::min(zwk * rmo, zmax));
                } else {
                    els = karman * zwk * std::pow(1.0f - alp4 * zwk * rmo, 0.2f);
                }

                float wt = 0.5f * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5f;

                el[k] = std::min(elb / (elb / elt + elb / els + 1.0f), elf);
            }
            break;

        /* Nonlocal (using BouLac) Form of Mixing Length */
        case 1:
            ugrid = std::sqrt(u1[kts] * u1[kts] + v1[kts] * v1[kts]);
            uonset = 15.0;
            wt_u = (1.0f - std::min(std::max(float(ugrid - uonset), 0.0f) / 30.0f, 0.5f));
            cns = 3.5;
            alp1 = 0.23;
            alp2 = 0.3;
            alp3 = 2.5f * wt_u; // taper off buoyancy enhancement in shear-driven pbls
            alp4 = 5.0;
            alp5 = 0.3;
            alp6 = 50.0;

            // Impose limits on the height integration for elt and the transition layer depth
            zi2 = std::max(float(zi), 300.f);
            h1 = std::max(float(0.3f * zi2), 300.0f);
            h1 = std::min(float(h1), 600.0f);
            h2 = h1 / 2.0f;

            qtke[kts] = std::max(float(0.5f * qke[kts]), 0.5f*qkemin);
            thetaw[kts] = theta[kts];
            qkw[kts] = std::sqrt(std::max(float(qke[kts]), qkemin));

            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0f - afk;
                qkw[k] = std::sqrt(std::max(float(qke[k] * abk + qke[k-1] * afk), qkemin));
                qtke[k] = std::max(0.5f * (qkw[k] * qkw[k]),0.005f);
                thetaw[k] = theta[k] * abk + theta[k-1] * afk;
            }

            elt = 1.0e-5;
            vsc = 1.0e-5;

            // ** Strictly, zwk*h[i,j] -> ( zwk*h[i,j]+z0 ) **
            k = kts + 1;
            zwk = zw[k];
            printf("elt %d %15.25g %15.15g %15.15g %15.15g %15.15g\n",k,elt,vsc,zwk,zi2,h1);
            while (zwk <= zi2 + h1) {
                dzk = 0.5f * (dz[k] + dz[k-1]);
                qdz = std::min(std::max(float(qkw[k] - qmin), 0.01f), 30.0f) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
		printf("dzk %d %15.15g %15.15g %15.15g %15.15g %15.15g\n",k,dzk,qdz,elt,vsc,zwk);
            }

            elt = std::min(std::max(float(alp1 * elt / vsc), 8.0f), 400.0f);
	    printf("elt %d %15.15g %15.15g %15.15g %15.15g %15.15g\n",elt,alp1,alp1*elt,alp1*elt/vsc,8.,400.);
            vflx = fltv;
            vsc = std::cbrt(gtr * elt * std::max(float(vflx), 0.0f));

            // ** Strictly, el[i,j,0] is not zero **
            el[kts] = 0.0;
            zwk1 = zw[kts+1];

            // Compute BouLac mixing length
            boulac_length_cc(kts,kte,zw,dz,qtke,thetaw,elblmin,elblavg,gtr);

            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];

                // ** Length scale is limited by the buoyancy effect **
                if (dtv[k] > 0.0f) {
                    bv = std::max(float(std::sqrt(gtr * dtv[k])), 0.0001f);
                    elb = std::max(float(alp2 * std::max(qkw[k],qke_elb_min)), float(alp6 * edmf_a1[k-1] * edmf_w1[k-1])) / bv * (1.0f + alp3 * std::sqrt(vsc / (bv * elt)));
                    elb = std::min(elb, zwk);
                    elf = 1.0f * std::max(qkw[k],qke_elb_min) / bv;
                    elblavg[k] = std::max(float(elblavg[k]), float(alp6 * edmf_a1[k-1] * edmf_w1[k-1] / bv));
                } else {
                    elb = 1.0e10;
                    elf = elb;
                }
                if (rmo > 0.0f) {
                    els = karman * zwk / (1.0f + cns * std::min(zwk * rmo, zmax));
                } else {
                    els = karman * zwk * std::pow(1.0f - alp4 * zwk * rmo, 0.2f);
                }
                float wt = 0.5f * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5f;
                el[k] = std::sqrt((els*els)/(1.0f + (els*els)/(elt*elt)));
                el[k] = std::min(el[k], elb);
                el[k] = std::min(el[k], elf);
                el[k] = el[k]*(1.0f-wt) + alp5*elblavg[k]*wt;
                el[k] = el[k] * psig_bl;
            }
            break;

        /* Local (mostly) Mixing Length Formulation */
        case 2:
            uonset = 3.5f + dz[kts] * 0.1f;
            ugrid = std::sqrt(u1[kts] * u1[kts] + v1[kts] * v1[kts]);
            cns = 3.5;
            alp1 = 0.22;
            alp2 = 0.30;
            alp3 = 2.0;
            alp4 = 5.0;
            alp5 = alp2; // like alp2, but for free atmosphere
            alp6 = 50.0; // used for MF mixing length
            zi2 = std::max(float(zi), float(minzi));

            // Impose limits on the height integration for elt and the transition layer depth
            h1 = std::max(float(0.3f * zi2), 300.0f);
            h1 = std::min(float(h1), 600.0f);           // 1/2 transition layer depth
            h2 = h1 * 0.5f;                             // 1/4 transition layer depth

            qtke[kts] = std::max(float(0.5f * qke[kts]), 0.5f*qkemin);
            qkw[kts] = std::sqrt(std::max(float(qke[kts]), qkemin));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0f - afk;
                qkw[k] = std::sqrt(std::max(float(qke[k] * abk + qke[k-1] * afk), qkemin));
                qtke[k] = 0.5f * qkw[k] * qkw[k];
            }

            elt = 1.0e-5;
            vsc = 1.0e-5;

            // ** Strictly, zwk*h[i,j] -> ( zwk*h[i,j]+z0 ) **
            pblh_plus_ent = std::max(zi+h1, 100.f);
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= pblh_plus_ent) {
                dzk = 0.5f * (dz[k] + dz[k-1]);
                qdz = std::min(std::max(float(qkw[k] - qmin), 0.03f), 30.0f) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }

            elt = std::min(std::max(float(alp1 * elt / vsc), 10.0f), 400.0f);
            vflx = fltv;
            vsc = std::cbrt(gtr * elt * std::max(float(vflx), 0.0f));

            // ** Strictly, el[i,j,0] is not zero **
            el[kts] = 0.0;
            zwk1 = zw[kts+1];

            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];
                dzk = 0.5f * (dz[k] + dz[k-1]);
                cldavg = 0.5f * (cldfra_bl1d[k-1] + cldfra_bl1d[k]);

                // ** Length scale limited by the buoyancy effect **
                if (dtv[k] > 0.0f) {
                    bv = std::max(float(std::sqrt(gtr * dtv[k])), 0.001f);
                    elb_mf = std::max(float(alp2 * qkw[k]), float(alp6 * edmf_a1[k-1] * edmf_w1[k-1]) / bv * (1.0f + alp3 * std::sqrt(vsc / (bv * elt))));
                    elb = std::min(std::max(float(alp5 * qkw[k]), float(alp6 * edmf_a1[k] * edmf_w1[k]) / bv), float(zwk));

                    wstar = 1.25f * std::cbrt(gtr * zi * std::max(float(vflx), 1.0e-4f));
                    tau_cloud = std::min(std::max(float(ctau * wstar / grav), 30.0f), 150.0f);

                    // minimize influence of surface heat flux on tau far away from the PBLH
                    float wt = 0.5f * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5f;
                    tau_cloud = tau_cloud * (1.0f - wt) + 50.0f * wt;
                    elf = std::min(std::max(float(tau_cloud * std::sqrt(std::min(float(qtke[k]), 40.0f))), float(alp6 * edmf_a1[k] * edmf_w1[k] / bv)), float(zwk));
                } else {
                    /*
                     * use version in development for RAP/HRRR 2016
                     * JAYMES-
                     * tau_cloud is an eddy turnover timescale;
                     * see Teixeira and Cheinet (2004), Eq. 1, and
                     * Cheinet and Teixeira (2003), Eq. 7.  The
                     * coefficient 0.5 is tuneable. Expression in
                     * denominator is identical to vsc (a convective
                     * velocity scale), except that elt is relpaced
                     * by zi, and zero is replaced by 1.0e-4 to
                     * prevent division by zero.
                     */
                    wstar = 1.25f * std::cbrt(gtr * zi * std::max(float(vflx), 1.0e-4f));
                    tau_cloud = std::min(std::max(float(ctau * wstar / grav), 50.0f), 200.0f);

                    // minimize influence of surface heat flux on tau far away from the PBLH
                    float wt = 0.5f * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5f;
                    tau_cloud = tau_cloud * (1.0f - wt) + std::max(100.0f, dzk * 0.25f) * wt;

                    elb = std::min(tau_cloud * std::sqrt(std::min(qtke[k], 40.0f)), zwk);
                    elf = elb;
                    elb_mf = elb;
                }
                elf = elf / (1.0f + (elf / 800.0f));
                elb_mf = std::max(float(elb_mf), 0.01f);

                // ** Length scale in the surface layer **
                if (rmo > 0.0f) {
                    els = karman * zwk / (1.0f + cns * std::min(zwk * rmo, zmax));
                } else {
                    els = karman * zwk * std::pow(1.0f - alp4 * zwk * rmo, 0.2f);
                }

                // ** Now blend the mixing length scales **
                float wt = 0.5f * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5f;

                // try squared-blending
                el[k] = std::sqrt(els * els / (1.0f + (els * els / (elt * elt)) + (els * els / (elb_mf * elb_mf))));
                el[k] = el[k] * (1.0f - wt) + elf * wt;

                // include scale-awareness. For now, use simple asymptotic kz -> 12 m (should be ~ dz)
                el_les = std::min(els / (1.0f + (els/12.f)), elb_mf);
                el[k] = el[k] * psig_bl + (1.0f - psig_bl) * el_les;
            }
            break;
    }
}




// called from driver 
void moisture_check_cc(int kte, float delt, float* dp, const float* exner,
                    float* qv, float* qc, float* qi, float* qs, float* th,
                    float* dqv, float* dqc, float* dqi, float* dqs, float* dth, 
		    float xlvcp, float xlscp) {

    // constants (assuming xlvcp and xlscp are defined elsewhere)
    const float qvmin = 1e-20, qcmin = 0.0, qimin = 0.0;
    float dqv2;

    for (int k = kte; k >= 0; --k) { // from the top to the surface
        float dqc2 = std::max(0.0f, qcmin - qc[k]);
        float dqi2 = std::max(0.0f, qimin - qi[k]);
        float dqs2 = std::max(0.0f, qimin - qs[k]);

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
        float dqv2 = std::max(0.0f, qvmin - qv[k]);
        dqv[k] += dqv2 / delt;
        qv[k] += dqv2;
        if (k != 0) {
            qv[k-1] -= dqv2 * dp[k] / dp[k-1];
            dqv[k-1] -= dqv2 * dp[k] / dp[k-1] / delt;
        }
        qv[k] = std::max(float(qv[k]), float(qvmin));
        qc[k] = std::max(float(qc[k]), float(qcmin));
        qi[k] = std::max(float(qi[k]), float(qimin));
        qs[k] = std::max(float(qs[k]), float(qimin));
    }

        float sum = 0.0;
    float aa, dum;

    // only execute if dqv2 > 1.e-20, which indicates adjustment was made at the top layer
    if(dqv2 > 1e-20) {
        for (int k = 0; k <= kte; ++k) { // loop through all layers
            if (qv[k] > 2.0f * qvmin) {
                sum += qv[k] * dp[k];
            }
        }

        aa = dqv2 * dp[0] / std::max(1.e-20f, sum); // adjust for 1-based indexing with dp[0]

        if (aa < 0.5f) {
            for (int k = 0; k <= kte; ++k) { // loop through all layers again
                if (qv[k] > 2.0f * qvmin) {
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
void mym_predict_cc(int& kts, int& kte, float& closure, float& delt, float* dz, float* ust, float& flt, float& flq, float& pmz, float& phh, float* el, float* dfq, float* rho, float* pdk, float* pdt, float* pdq, float* pdc, float* qke, float* tsq, float* qsq, float* cov, float* s_aw, float* s_awqke, int& bl_mynn_edmf_tke, float* qwt1d, float* qdiss1d, int& tke_budget, float& xlvcp, float& xlscp, float& karman) {
    float vkz, pdk1, phm, pdt1, pdq1, pdc1, b1l, b2l, onoff;
    float dtz[kte-kts+1];
    float a[kte-kts+1];
    float b[kte-kts+1];
    float c[kte-kts+1];
    float d[kte-kts+1];
    float x[kte-kts+1];
    float rhoinv[kte-kts+1];
    float rhoz[kte-kts+2];
    float kqdz[kte-kts+2];
    float kmdz[kte-kts+2];
    float qkw[kte-kts+1];
    float bp[kte-kts+1];
    float rp[kte-kts+1];
    float df3q[kte-kts+1];
    float tke_up[kte-kts+1];
    float dzinv[kte-kts+1];
    
    // regulate the momentum mixing from the mass-flux scheme (on or off)
    if (bl_mynn_edmf_tke == 0) {
        onoff = 0.0;
    } else {
        onoff = 1.0;
    }
    
    // calculate vkz
    vkz = karman * 0.5f * dz[kts];
    
    // calculate df3q and dtz
    for (int k = kts; k <= kte; k++) {
        qkw[k] = sqrt(std::max(qke[k], 0.0f));
        df3q[k] = sqfac * dfq[k];
        dtz[k] = delt / dz[k];
    }
    
    // prepare "constants" for diffusion equation
    rhoz[kts] = rho[kts];
    rhoinv[kts] = 1.0f / rho[kts];
    kqdz[kts] = rhoz[kts] * df3q[kts];
    kmdz[kts] = rhoz[kts] * dfq[kts];
    for (int k = kts+1; k <= kte; k++) {
        rhoz[k] = (rho[k] * dz[k-1] + rho[k-1] * dz[k]) / (dz[k-1] + dz[k]);
        rhoz[k] = std::max(rhoz[k], 1e-4f);
        rhoinv[k] = 1.0f / std::max(rho[k], 1e-4f);
        kqdz[k] = rhoz[k] * df3q[k];
        kmdz[k] = rhoz[k] * dfq[k];
    }
    rhoz[kte+1] = rhoz[kte];
    kqdz[kte+1] = rhoz[kte+1] * df3q[kte];
    kmdz[kte+1] = rhoz[kte+1] * dfq[kte];
    
    // calculate pdk1, phm, pdt1, pdq1, pdc1
    pdk1 = 2.0f * (ust[0]*ust[0]*ust[0]) * pmz / vkz;
    phm = 2.0f / ust[0] * phh / vkz;
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
        b1l = b1 * 0.5f * (el[k+1] + el[k]);
        bp[k] = 2.0f * qkw[k] / b1l;
        rp[k] = pdk[k+1] + pdk[k];
    }
    for (int k = kts; k <= kte-1; k++) {
        a[k] = -dtz[k] * kqdz[k] * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * s_aw[k] * onoff;
        b[k] = 1.0f + dtz[k] * (kqdz[k] + kqdz[k+1]) * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) * onoff + bp[k] * delt;
        c[k] = -dtz[k] * kqdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff;
        d[k] = rp[k] * delt + qke[k] + dtz[k] * rhoinv[k] * (s_awqke[k] - s_awqke[k+1]) * onoff;
    }
    /*
      for (int k = kts; k <= kte-1; k++) {
        a[k] = -dtz[k] * df3q[k] + 0.5f * dtz[k] * s_aw[k] * onoff;
        b[k] = 1.0f + dtz[k] * (df3q[k] + df3q[k+1]) + 0.5f * dtz[k] * (s_aw[k] - s_aw[k+1]) * onoff + bp[k] * delt;
        c[k] = -dtz[k] * df3q[k+1] - 0.5f * dtz[k] * s_aw[k+1] * onoff;
        d[k] = rp[k] * delt + qke[k] + dtz[k] * (s_awqke[k] - s_awqke[k+1]) * onoff;
	}*/
    a[kte] = 0.0;
    b[kte] = 1.0;
    c[kte] = 0.0;
    d[kte] = qke[kte];

    tridiag2_cc(kte, a, b, c, d, x);

    for (int k = kts; k <= kte; k++) {
        qke[k] = std::max(x[k], qkemin);
        qke[k] = std::min(qke[k], 150.0f);
    }

    // tke budget
    if (tke_budget == 1) {
        float tke_up[kte-kts+1];
        float dzinv[kte-kts+1];
        
        // tke vertical transport
	for (int k=kts; k <=kte; k++) 
	{
		tke_up[k] = 0.5f * qke[k];
                dzinv[k] = 1.0f / dz[k];
	}
        qwt1d[kts] = dzinv[kts] * ((kqdz[kts+1] * (tke_up[kts+1] - tke_up[kts]) - (kqdz[kts] * tke_up[kts])) + 0.5f * rhoinv[kts] * (s_aw[kts+1] * tke_up[kts+1] + ((s_aw[kts+1] - s_aw[kts]) * tke_up[kts]) + (s_awqke[kts] - s_awqke[kts+1])) * onoff);
        for (int k = kts+1; k <= kte-1; k++) {
	  qwt1d[k] = dzinv[k] * ((kqdz[k+1] * (tke_up[k+1] - tke_up[k]) - (kqdz[k] * (tke_up[k] - tke_up[k-1]))) + 0.5f * rhoinv[k] * (s_aw[k+1] * tke_up[k+1] + ((s_aw[k+1] - s_aw[k]) * tke_up[k]) - (s_aw[k] * tke_up[k-1]) + (s_awqke[k] - s_awqke[k+1])) * onoff);
        }
    qwt1d[kte] = dzinv[kte] * (-(kqdz[kte] * (tke_up[kte] - tke_up[kte-1])) + 0.5f * rhoinv[kte] * (-(s_aw[kte] * tke_up[kte]) - (s_aw[kte] * tke_up[kte-1]) + s_awqke[kte]) * onoff);
        
        // tke dissipation rate
	for (int k=kts; k <=kte; k++) 
	{
		qdiss1d[k] = bp[k] * tke_up[k];
	}
    }
    
    if (closure > 2.5) {
        // prediction of the moisture variance
        for (int k = kts; k <= kte-1; k++) {
            b2l = b2 * 0.5f * (el[k+1] + el[k]);
            bp[k] = 2.0f * qkw[k] / b2l;
            rp[k] = pdq[k+1] + pdq[k];
        }
        for (int k = kts; k <= kte-1; k++) {
            a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
            b[k] = 1.0f + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + bp[k] * delt;
            c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k];
            d[k] = rp[k] * delt + qsq[k];
        }
        a[kte] = -1.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = 0.0;
        tridiag2_cc(kte, a, b, c, d, x);
        for (int k = kts; k <= kte; k++) {
            qsq[k] = std::max(x[k], 1e-17f);
        }
    } else {
        // level 2.5f - use level 2 diagnostic
        for (int k = kts; k <= kte-1; k++) {
            if (qkw[k] <= 0.0f) {
                b2l = 0.0;
            } else {
                b2l = b2 * 0.25f * (el[k+1] + el[k]) / qkw[k];
            }
            qsq[k] = b2l * (pdq[k+1] + pdq[k]);
        }
        qsq[kte] = qsq[kte-1];
    }
    
    if (closure >= 3.0f) {
        // prediction of the temperature variance
        for (int k = kts; k <= kte-1; k++) {
            b2l = b2 * 0.5f * (el[k+1] + el[k]);
            bp[k] = 2.0f * qkw[k] / b2l;
            rp[k] = pdt[k+1] + pdt[k];
        }
        for (int k = kts; k <= kte-1; k++) {
            a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
            b[k] = 1.0f + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + bp[k] * delt;
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
            b2l = b2 * 0.5f * (el[k+1] + el[k]);
            bp[k] = 2.0f * qkw[k] / b2l;
            rp[k] = pdc[k+1] + pdc[k];
        }
        for (int k = kts; k <= kte-1; k++) {
            a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
            b[k] = 1.0f + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + bp[k] * delt;
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
            if (qkw[k] <= 0.0f) {
                b2l = 0.0;
            } else {
                b2l = b2 * 0.25f * (el[k+1] + el[k]) / qkw[k];
            }
            tsq[k] = b2l * (pdt[k+1] + pdt[k]);
            cov[k] = b2l * (pdc[k+1] + pdc[k]);
        }
        tsq[kte] = tsq[kte-1];
        cov[kte] = cov[kte-1];
    }
}

void mynn_mix_chem_cc(int kts, int kte, int i,
                   float delt, float* dz, float pblh,
                   int nchem, int kdvel, int ndvel,
                   float** chem1, float* vd1,
                   float* rho,
                   float flt, float* tcd, float* qcd,
                   float* dfh,
                   float* s_aw, float** s_awchem,
                   float emis_ant_no, float frp, int rrfs_sd, int enh_mix) {

    // local vars
    float dtz[kte - kts + 1];
    float a[kte - kts + 1], b[kte - kts + 1], c[kte - kts + 1], d[kte - kts + 1], x[kte - kts + 1];
    float dztop = 0.5f * (dz[kte - 1] + dz[kte - 2]);
    for (int k = kts; k <= kte; ++k) {
        dtz[k - kts] = delt / dz[k - 1];
    }
    // prepare "constants" for diffusion equation.
    float rhoz[kte - kts + 2], khdz[kte - kts + 2], rhoinv[kte - kts + 1];
    rhoz[0] = rho[kts - 1];
    rhoinv[0] = 1.0f / rho[kts - 1];
    khdz[0] = rhoz[0] * dfh[kts - 1];
    for (int k = kts + 1; k <= kte; ++k) {
        rhoz[k - kts] = (rho[k - 1] * dz[k - 2] + rho[k - 2] * dz[k - 1]) / (dz[k - 2] + dz[k - 1]);
        rhoz[k - kts] = std::max(float(rhoz[k - kts]), 1e-4f);
        rhoinv[k - kts] = 1.0f / std::max(float(rho[k - 1]), 1e-4f);
        float dzk = 0.5f * (dz[k - 1] + dz[k - 2]);
        khdz[k - kts] = rhoz[k - kts] * dfh[k - 1];
    }
    rhoz[kte - kts + 1] = rhoz[kte - kts];
    khdz[kte - kts + 1] = rhoz[kte - kts + 1] * dfh[kte - 1];
    // stability criteria for mf
    for (int k = kts + 1; k <= kte - 1; ++k) {
        khdz[k - kts] = std::max(float(khdz[k - kts]), float(0.5f * s_aw[k - kts]));
        khdz[k - kts] = std::max(float(khdz[k - kts]), float(-0.5f * (s_aw[k - kts] - s_aw[k - kts + 1])));
    }
    // enhanced mixing over fires
    if (rrfs_sd==1 && enh_mix==1) {
        for (int k = kts + 1; k <= kte - 1; ++k) {
            float khdz_old = khdz[k - kts];
            float khdz_back = pblh * 0.15f / dz[k - 1];
            // modify based on anthropogenic emissions of no and frp
            if (pblh < pblh_threshold) {
                if (emis_ant_no > no_threshold) {
                    khdz[k - kts] = std::max(1.1f * float(khdz[k - kts]), float(std::sqrt((emis_ant_no / no_threshold)) / dz[k - 1] * rhoz[k - kts]));
                }
                if (frp > frp_threshold) {
                    int kmaxfire = std::ceil(std::log(frp));
                    khdz[k - kts] = std::max(float(1.1f * khdz[k - kts]), float((1.0f - k / (kmaxfire * 2.0f)) * (std::pow(std::log(frp), 2.0f) - 2.0f * std::log(frp)) / dz[k - 1] * rhoz[k - kts]));
                }
            }
        }
    }
    // mixing of chemical species
    for (int ic = 0; ic < nchem; ++ic) {
        int k = kts;
        a[0] = -dtz[0] * khdz[0] * rhoinv[0];
        b[0] = 1.0f + dtz[0] * (khdz[1] + khdz[0]) * rhoinv[0] - 0.5f * dtz[0] * rhoinv[0] * s_aw[1];
        c[0] = -dtz[0] * khdz[1] * rhoinv[0] - 0.5f * dtz[0] * rhoinv[0] * s_aw[1];
        d[0] = chem1[k - 1][ic] - dtz[0] * vd1[ic] * chem1[k - 1][ic] - dtz[0] * rhoinv[0] * s_awchem[1][ic];
        for (k = kts + 1; k <= kte - 1; ++k) {
            a[k - kts] = -dtz[k - kts] * khdz[k - kts] * rhoinv[k - kts] + 0.5f * dtz[k - kts] * rhoinv[k - kts] * s_aw[k - kts];
            b[k - kts] = 1.0f + dtz[k - kts] * (khdz[k - kts] + khdz[k - kts + 1]) * rhoinv[k - kts] + 0.5f * dtz[k - kts] * rhoinv[k - kts] * (s_aw[k - kts] - s_aw[k - kts + 1]);
            c[k - kts] = -dtz[k - kts] * khdz[k - kts + 1] * rhoinv[k - kts] - 0.5f * dtz[k - kts] * rhoinv[k - kts] * s_aw[k - kts + 1];
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
void mynn_tendencies_cc(const int& kts,const int& kte, const float & delt,
				   /*in*/ const float* dz,
				   /*in*/ const float* rho,
			/*in*/ const float* u, const float* v, const float* th, const float* tk, const float* qv,
				   /*in*/ const float* qc, const float* qi, const float* qs, const float* qnc, const float* qni,
			/*in*/ const float* psfc,const float* p,const float* exner,
                                   /*inout*/ float* thl, float* sqv, float* sqc, float* sqi,
                                   /*inout*/ float* sqs, float* sqw, float* qnwfa, float* qnifa, float* qnbca, float* ozone,
                                   /*in*/ float* ust, const float & flt,const float & flq,const float & flqv,const float & flqc,const float & wspd, const float & uoce, const float & voce,
				   /*in*/ const float* tcd,const float* qcd,
				   /*inout*/ float* dfm, float* dfh,
				   /*inout*/ float* du, float* dv, float* dth,
				   /*inout*/ float* dqv, float* dqc, float* dqi,
				   /*inout*/ float* dqs, float* dqnc, float* dqni,
				   /*inout*/ float* dqnwfa, float* dqnifa, float* dqnbca,
				   /*inout*/ float* dozone,
				   /*in*/ const float* diss_heat,
			/*in*/ const float* s_aw, const float* s_awthl, const float* s_awqt, const float* s_awqv, const float* s_awqc, const float* s_awu, const float* s_awv, const float* s_awqnc, const float* s_awqni, const float* s_awqnwfa, const float* s_awqnifa, const float* s_awqnbca, const float* sd_aw, const float* sd_awthl, const float* sd_awqt, const float* sd_awqv, const float* sd_awqc, const float* sd_awu, const float* sd_awv,
				   const float* sub_thl,const float* sub_sqv,const float* sub_u,const float* sub_v,
				   const float* det_thl,const float* det_sqv,const float* det_sqc,const float* det_u,const float* det_v,
				   /*logical turned into int */const int& flag_qc, const int& flag_qi, const int& flag_qnc, const int& flag_qni, const int& flag_qs, const int& flag_qnwfa, const int& flag_qnifa, const int& flag_qnbca, const int& flag_ozone,
			const int & bl_mynn_cloudmix, const int & bl_mynn_mixqt, int & bl_mynn_edmf_mom,  const int & bl_mynn_mixscalars, /* new */const int& debug_code,const float& r_d,const float& p608,const float& ep_2,const float& ep_3,const float& tv0,const float& xlv,const float& xlvcp,const float& xlscp) {
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
     
    float nonloc = 1.0;

    float dztop = 0.5f * (dz[kte-1] + dz[kte-2]);
    float onoff = bl_mynn_edmf_mom;
    onoff = (onoff == 0) ? 0.0f : 1.0;
    float rhosfc = *psfc / (r_d * (tk[kts] + p608 * qv[kts]));

    float dtz[kte+2];
    float dfhc[kte+2]; 
    float dfmc[kte+2];
    float delp[kte+2]; 
    float sqv2[kte+2]; 
    float sqc2[kte+2];
    float sqs2[kte+2]; 
    float sqi2[kte+2];
    float qnc2[kte+2];
    float qni2[kte+2];
    float rhoz[kte+2];
    float khdz[kte+2];
    float kmdz[kte+2];
    float rhoinv[kte+2];
    float sqw2[kte+2];
    float a[kte+2];
    float b[kte+2]; 
    float c[kte+2];
    float d[kte+2];
    float x[kte+2];
    float qvflux;
    float ust_v=*ust;

    dtz[kts] = delt / dz[kts];
    rhoz[kts] = rho[kts];
    rhoinv[kts] = 1.0f / rho[kts];
    khdz[kts] = rhoz[kts] * dfh[kts];
    kmdz[kts] = rhoz[kts] * dfm[kts];
    delp[kts] = *psfc - (p[kts+1] * dz[kts] + p[kts] * dz[kts+1]) / (dz[kts] + dz[kts+1]);
    
    for (int k = kts+1; k <= kte; k++) {
        dtz[k] = delt / dz[k];
        rhoz[k] = (rho[k] * dz[k-1] + rho[k-1] * dz[k]) / (dz[k-1] + dz[k]);
        rhoz[k] = std::max(float(rhoz[k]), 1e-4f);
        rhoinv[k] = 1.0f / std::max(float(rho[k]), 1e-4f);
        float dzk = 0.5f * (dz[k] + dz[k-1]);
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
        khdz[k] = std::max(float(khdz[k]), float(0.5f * s_aw[k]));
        khdz[k] = std::max(float(khdz[k]), float(-0.5f * (s_aw[k] - s_aw[k+1])));
        kmdz[k] = std::max(float(kmdz[k]), float(0.5f * s_aw[k]));
        kmdz[k] = std::max(float(kmdz[k]), float(-0.5f * (s_aw[k] - s_aw[k+1])));
    }
    
    float ustdrag = std::min(ust_v * ust_v, 0.99f) / wspd;
    float ustdiff = std::min(ust_v * ust_v, 0.01f) / wspd;
    
    for (int k = kts; k <= kte; k++) {
        dth[k] = 0.0;
    }

    int k = kts;
    a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
    b[k] = 1.0f + dtz[k] * (kmdz[k+1] + rhosfc * ust_v * ust_v / wspd) * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
    c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
    d[k] = u[k] + dtz[k] * uoce * ust_v * ust_v / wspd - dtz[k] * rhoinv[k] * s_awu[k+1] * onoff + dtz[k] * rhoinv[k] * sd_awu[k+1] * onoff + sub_u[k] * delt + det_u[k] * delt;
    
    for (int k = kts+1; k <= kte-1; k++) {
        a[k] = -dtz[k] * kmdz[k] * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * s_aw[k] * onoff + 0.5f * dtz[k] * rhoinv[k] * sd_aw[k] * onoff;
        b[k] = 1.0f + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) * onoff + 0.5f * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]) * onoff;
        c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
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
    b[k] = 1.0f + dtz[k] * (kmdz[k+1] + rhosfc * ust_v * ust_v / wspd) * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
    c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
    d[k] = v[k] + dtz[k] * voce * ust_v * ust_v / wspd - dtz[k] * rhoinv[k] * s_awv[k+1] * onoff + dtz[k] * rhoinv[k] * sd_awv[k+1] * onoff + sub_v[k] * delt + det_v[k] * delt;
    
    for (int k = kts+1; k <= kte-1; k++) {
        a[k] = -dtz[k] * kmdz[k] * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * s_aw[k] * onoff + 0.5f * dtz[k] * rhoinv[k] * sd_aw[k] * onoff;
        b[k] = 1.0f + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) * onoff + 0.5f * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]) * onoff;
        c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
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
    b[k] = 1.0f + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1];
    c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1];
    d[k] = thl[k] + dtz[k] * rhosfc * flt * rhoinv[k] + tcd[k] * delt - dtz[k] * rhoinv[k] * s_awthl[k+1] - dtz[k] * rhoinv[k] * sd_awthl[k+1] + diss_heat[k] * delt + sub_thl[k] * delt + det_thl[k] * delt;
    
    for (int k = kts+1; k <= kte-1; k++) {
        a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * s_aw[k] + 0.5f * dtz[k] * rhoinv[k] * sd_aw[k];
        b[k] = 1.0f + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) + 0.5f * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]);
        c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1];
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
        b[k] = 1.0f + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1];
        c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1];
        d[k] = sqw[k] + dtz[k] * rhosfc * flq * rhoinv[k] + qcd[k] * delt - dtz[k] * rhoinv[k] * s_awqt[k+1] - dtz[k] * rhoinv[k] * sd_awqt[k+1];
        
        for (int k = kts+1; k <= kte-1; k++) {
            a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * s_aw[k] + 0.5f * dtz[k] * rhoinv[k] * sd_aw[k];
            b[k] = 1.0f + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) + 0.5f * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]);
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1];
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
            b[k] = 1.0f + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1];
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1];
            d[k] = sqc[k] + dtz[k] * rhosfc * flqc * rhoinv[k] + qcd[k] * delt - dtz[k] * rhoinv[k] * s_awqc[k+1] - dtz[k] * rhoinv[k] * sd_awqc[k+1] + det_sqc[k] * delt;
            
            for (int k = kts+1; k <= kte-1; k++) {
                a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * s_aw[k] + 0.5f * dtz[k] * rhoinv[k] * sd_aw[k];
                b[k] = 1.0f + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) + 0.5f * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]);
                c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1];
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
        if (qvflux < 0.0f) {
            qvflux = std::max(float(qvflux), (std::min(0.9f * float(sqv[kts]) - 1e-8f, 0.0f) / dtz[kts]));
        }
        a[k] = -dtz[k] * khdz[k] * rhoinv[k];
        b[k] = 1.0f + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1];
        c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1];
        d[k] = sqv[k] + dtz[k] * rhosfc * qvflux * rhoinv[k] + qcd[k] * delt - dtz[k] * rhoinv[k] * s_awqv[k+1] - dtz[k] * rhoinv[k] * sd_awqv[k+1] + sub_sqv[k] * delt + det_sqv[k] * delt;

        for (int k = kts+1; k <= kte-1; k++) {
            a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * s_aw[k] + 0.5f * dtz[k] * rhoinv[k] * sd_aw[k];
            b[k] = 1.0f + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) + 0.5f * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]);
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5f * dtz[k] * rhoinv[k] * sd_aw[k+1];
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
        b[k] = 1.0f + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k];
        c[k] = -dtz[k] * khdz[k+1] * rhoinv[k];
        d[k] = sqi[k];

        for (int k = kts+1; k <= kte-1; k++) {
	    a[k] = -dtz[k] * khdz[k] * rhoinv[k];
            b[k] = 1.0f + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k];
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
        b[k] = 1.0f + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k];
        c[k] = -dtz[k] * khdz[k+1] * rhoinv[k];
        d[k] = sqs[k];
        
        for (int k = kts+1; k <= kte-1; k++) {
	    a[k] = -dtz[k] * khdz[k] * rhoinv[k];
            b[k] = 1.0f + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k];
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
            b[k] = 1.0f + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] * nonloc;
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] * nonloc;
            d[k] = qni[k] - dtz[k] * rhoinv[k] * s_awqni[k+1]*nonloc;
            
            for (int k = kts+1; k <= kte-1; k++) {
                a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * s_aw[k];
                b[k] = 1.0f + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) * nonloc;
                c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] * nonloc;
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
            b[k] = 1.0f + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] * nonloc;
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] * nonloc;
            d[k] = qni[k] - dtz[k] * rhoinv[k] * s_awqnc[k+1]*nonloc;
            
            for (int k = kts+1; k <= kte-1; k++) {
                a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * s_aw[k];
                b[k] = 1.0f + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5f * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) * nonloc;
                c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5f * dtz[k] * rhoinv[k] * s_aw[k+1] * nonloc;
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
      if(dozone[k]*delt+ozone[k]<0.0f)
	dozone[k]=-ozone[k]*0.99f/delt;
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
	dqnwfa[k]=0.0f;
	dqnifa[k]=0.0f;
      }
    }
    
    if(flag_qnbca && bl_mynn_mixscalars > 0) {
      for (int k = kts; k <= kte; k++) {
	dqnbca[k]=(qnbca2[k] - qnbca[k])/delt;
      }
    } else {
      for (int k = kts; k <= kte; k++) {
	dqnbca[k]=0.0f;
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
    //    /*inout*/ float* thl, float* sqv, float* sqc, float* sqi,
    //    /*inout*/ float* sqs, float* sqw, float* qnwfa, float* qnifa, float* qnbca, float* ozone,
    //    /*inout*/ float* dfm, float* dfh,
    //    /*inout*/ float* du, float* dv, float* dth,
    
}



void mym_condensation_cc(
    const int& kts,const int& kte,
    const float& dx, float* dz, float* zw, float& xland,
    float* thl, float* qw, float* qv, float* qc, float* qi, float* qs,
    float* p, float* exner, float* tsq, float* qsq, float* cov,
    float* sh, float* el, int& bl_mynn_cloudpdf,
    float* qc_bl1d, float* qi_bl1d, float* cldfra_bl1d,
    float& pblh1, float& hfx1,
    float* vt, float* vq, float* th, float* sgm, float* rmo,
    int &spp_pbl, float* rstoch_col,
    // extra constants here
    float ep_2, float ep_3, float xlv, float r_d, float xlvcp,
    float p608, float tv0, float r_v, float cice, float cliq,
    float cpv, float cp, float xls, float rcp)
{
    int k;
    float t3sq, r3sq, c3sq;
    float qsl, esat, qsat, dqsl, cld0, q1k, qlk, eq1, qll, q2p, pt, rac, qt, t, xl, rsl, cpm, fng, qww, alpha, beta, bb, ls, wt, wt2, qpct, cld_factor, fac_damp, liq_frac, ql_ice, ql_water, qmq, qsat_tk, q1_rh, rh_hack, dzm1, zsl, maxqc;
    const float qpct_sfc = 0.025;
    const float qpct_pbl = 0.030;
    const float qpct_trp = 0.040;
    const float rhcrit = 0.83;
    const float rhmax = 1.02;
    float erf;
    float dth, dtl, dqw, dzk, els;
    float zagl, damp, pblh2;
    float cfmax;
    float theta1, theta2, ht1, ht2;
    float qw_pert;
    int k_tropo;

//real(float), dimension(kts:kte) :: alp,a,bet,b,ql,q1,rh
    float alp[kte-kts]; 
    float a[kte-kts]; 
    float bet[kte-kts]; 
    float b[kte-kts]; 
    float ql[kte-kts]; 
    float q1[kte-kts]; 
    float rh[kte-kts]; 

    // obtain an estimate for the tropopause height (k)
    for (k = kte - 3; k >= kts; k--) {
        theta1 = th[k];
        theta2 = th[k + 2];
        ht1 = 44307.692f * (1.0f - pow(p[k] / 101325.0, 0.190));
        ht2 = 44307.692f * (1.0f - pow(p[k + 2] / 101325.0, 0.190));
        if ((((theta2 - theta1) / (ht2 - ht1)) < 10.0f / 1500.0f) && (ht1 < 19000.0f) && (ht1 > 4000.0f)) {
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
                qsl = ep_2 * esat / std::max(1e-4f, (p[k] - ep_3 * esat));
                dqsl = qsl * ep_2 * xlv / (r_d * pow(t, 2));
                alp[k] = 1.0f / (1.0f + dqsl * xlvcp);
                bet[k] = dqsl * exner[k];
                t3sq = std::max(tsq[k], 0.0f);
                r3sq = std::max(qsq[k], 0.0f);
                c3sq = cov[k];
                c3sq = std::copysign(std::min(std::abs(c3sq), std::sqrt(t3sq * r3sq)), c3sq);
                r3sq = r3sq + bet[k] * bet[k] * t3sq - 2.0f * bet[k] * c3sq;
                qmq = qw[k] - qsl;
                sgm[k] = std::sqrt(std::max(r3sq, 1.0e-10f));
                q1[k] = qmq / sgm[k];
                cldfra_bl1d[k] = 0.5f * (1.0f + std::erf(q1[k] * rr2));
                q1k = q1[k];
                eq1 = rrp * std::exp(-0.5f * q1k * q1k);
                qll = std::max(cldfra_bl1d[k] * q1k + eq1, 0.0f);
                ql[k] = alp[k] * sgm[k] * qll;
                liq_frac = std::min(1.0f, std::max(0.0f, (t - 240.0f) / 29.0f));
                qc_bl1d[k] = liq_frac * ql[k];
                qi_bl1d[k] = (1.0f - liq_frac) * ql[k];
                q2p = xlvcp / exner[k];
                pt = thl[k] + q2p * ql[k];
                qt = 1.0f + p608 * qw[k] - (1.0f + p608) * (qc_bl1d[k] + qi_bl1d[k]) * cldfra_bl1d[k];
                rac = alp[k] * (cldfra_bl1d[k] - qll * eq1) * (q2p * qt - (1.0f + p608) * pt);
                vt[k] = qt - 1.0f - rac * bet[k];
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
                qsl = ep_2 * esat / std::max(1e-4f, (p[k] - ep_3 * esat));
                dqsl = qsl * ep_2 * xlv / (r_d * pow(t, 2));
                alp[k] = 1.0f / (1.0f + dqsl * xlvcp);
                bet[k] = dqsl * exner[k];
                if (k == kts) {
                    dzk = 0.5f * dz[k];
                } else {
                    dzk = dz[k];
                }
                dth = 0.5f * (thl[k + 1] + thl[k]) - 0.5f * (thl[k] + thl[std::max(k - 1, kts)]);
                dqw = 0.5f * (qw[k + 1] + qw[k]) - 0.5f * (qw[k] + qw[std::max(k - 1, kts)]);
                sgm[k] = std::sqrt(std::max(float((pow(alp[k], 2) * std::max(pow(el[k], 2), 0.1) * b2 * std::max(sh[k], 0.03f)) / 4.0f * pow((dqw / dzk - bet[k] * (dth / dzk)), 2)), 1.0e-10f));
                qmq = qw[k] - qsl;
                q1[k] = qmq / sgm[k];
                cldfra_bl1d[k] = 0.5f * (1.0f + std::erf(q1[k] * rr2));
                q1k = q1[k];
                eq1 = rrp * std::exp(-0.5f * q1k * q1k);
                qll = std::max(cldfra_bl1d[k] * q1k + eq1, 0.0f);
                ql[k] = alp[k] * sgm[k] * qll;
                liq_frac = std::min(1.0f, std::max(0.0f, (t - 240.0f) / 29.0f));
                qc_bl1d[k] = liq_frac * ql[k];
                qi_bl1d[k] = (1.0f - liq_frac) * ql[k];
                q2p = xlvcp / exner[k];
                pt = thl[k] + q2p * ql[k];
                qt = 1.0f + p608 * qw[k] - (1.0f + p608) * (qc_bl1d[k] + qi_bl1d[k]) * cldfra_bl1d[k];
                rac = alp[k] * (cldfra_bl1d[k] - qll * eq1) * (q2p * qt - (1.0f + p608) * pt);
                vt[k] = qt - 1.0f - rac * bet[k];
                vq[k] = p608 * pt - tv0 + rac;
		//    printf("vt[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",vt[k],qt,rac,bet[k],dqsl);
    //    printf("vq[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",vq[k],p608,pt,tv0,rac);
            }
            break;
        case 2:
        case -2: // diagnostic statistical scheme of Chaboureau and Bechtold (2002), JAS
                 // but with use of higher-order moments to estimate sigma
            pblh2 = std::max(10.0f, pblh1);
            zagl = 0.0f;
            dzm1 = 0.0f;
            for (k = kts; k < kte; k++) {
                zagl += 0.5f * (dz[k] + dzm1);
                dzm1 = dz[k];

                t = th[k] * exner[k];
                xl = xl_blend_cc(t,xlv,xls,cpv,cliq,cice);
                qsat_tk = qsat_blend_cc(t, p[k]);
                rh[k] = std::max(std::min(rhmax, qw[k] / std::max(1e-10f, qsat_tk)), 0.001f);

                // dqw/dT: Clasius-Clapeyron
                dqsl = qsat_tk * ep_2 * xlv / (r_d * std::pow(t, 2.0f));
                alp[k] = 1.0f / (1.0f + dqsl * xlvcp);
                bet[k] = dqsl * exner[k];

                rsl = xl * qsat_tk / (r_v * std::pow(t, 2.0f)); // slope of C-C curve at (= abs temperature) (CB02, Eqn. 4)
                cpm = cp + qw[k] * cpv;                 // CB02, sec. 2, para. 1
                a[k] = 1.0f / (1.0f + xl * rsl / cpm);  // CB02 variable "a"
                b[k] = a[k] * rsl;                      // CB02 variable "b"

                qw_pert = qw[k] + qw[k] * 0.5f * rstoch_col[k] * spp_pbl;

                // This form of qmq (the numerator of Q1) no longer uses the a[k] factor
                qmq = qw_pert - qsat_tk;

                // Use the form of Eq. 6 in CB02 except neglect all but the first term for sig_r
                r3sq = std::max(qsq[k], 0.0f);
                // Calculate sigma using higher-order moments
                sgm[k] = std::sqrt(r3sq);
                // Set constraints on sigma relative saturation water vapor
                sgm[k] = std::min(sgm[k], qsat_tk * 0.666f);

                // Introduce vertical grid spacing dependence on min sgm
                wt = std::max(500.0f - std::max(dz[k] - 100.0f, 0.0f), 0.0f) / 500.0f;
                sgm[k] += sgm[k] * 0.2f * (1.0f - wt); // inflate sgm for coarse dz

                // Allow min sgm to vary with dz and z
                qpct = qpct_pbl * wt + qpct_trp * (1.0f - wt);
                qpct = std::min(qpct, std::max(qpct_sfc, qpct_pbl * zagl / 500.0f));
                sgm[k] = std::max(sgm[k], qsat_tk * qpct);

                q1[k] = qmq / sgm[k]; // normalized saturation

                // Add condition for falling/settling into low-RH layers, so at least
                // some cloud fraction is applied for all qc, qs, and qi.
                rh_hack= rh[k];
                wt2    = std::min(std::max( zagl - pblh2, 0.0f )/300.0f, 1.0f);
                // ensure adequate RH & q1 when qi is at least 1e-9 (above the PBLH)
                if ((qi[k]+qs[k])>1.e-9 && (zagl > pblh2)) {
                  rh_hack =std::min(rhmax, rhcrit + wt2*0.045f*(9.0f + std::log10(qi[k]+qs[k])));
                  rh[k]   =std::max(rh[k], rh_hack);

                  q1_rh   =-3.0f + 3.0f*(rh[k]-rhcrit)/(1.0f-rhcrit);
                  q1[k]   =std::max(q1_rh, q1[k] );
                }
                // ensure adequate RH & q1 when qi is at least 1e-6 (above the PBLH)
                if (qc[k]>1.e-6 && (zagl > pblh2)) {
                  rh_hack =std::min(rhmax, rhcrit + wt2*0.08f*(6.0f + std::log10(qc[k])));
                  rh[k]   =std::max(rh[k], rh_hack);

                  q1_rh   =-3.0f + 3.0f*(rh[k]-rhcrit)/(1.0f-rhcrit);
                  q1[k]   =std::max(q1_rh, q1[k] );
                }

                q1k = q1[k]; // backup q1 for later modification

                // Specify cloud fraction
                //Original C-B cloud fraction, allows cloud fractions out to q1 = -3.5
                //cldfra_bl1D(K) = max(0., min(1., 0.5+0.36*atan(1.55*q1(k)))) ! Eq. 7 in CB02
                //Waynes LES fit  - over-diffuse, when limits removed from vt & vq & fng
                //cldfra_bl1D(K) = max(0., min(1., 0.5+0.36*atan(1.2*(q1(k)+0.4))))
                // Best compromise: Improves marine stratus without adding much cold bias.
                cldfra_bl1d[k] = std::max(0.0f, std::min(1.0f, 0.5f + 0.36f * std::atan(1.8f * (q1[k] + 0.2f))));

                // Specify hydrometeors
                // The cloud water formulatiosn are taken from CB02, Eq. 8
                maxqc = std::max(qw[k] - qsat_tk, 0.0f);
                if (q1k < 0.0f) {
                    ql_water = sgm[k] * std::exp(1.2f * q1k - 1.0f);
                    ql_ice = sgm[k] * std::exp(1.2f * q1k - 1.0f);
                } else if (q1k > 2.0f) {
                    ql_water = std::min(sgm[k] * q1k, maxqc);
                    ql_ice = sgm[k] * q1k;
                } else {
                    ql_water = std::min(float(sgm[k] * (std::exp(-1.0f) + 0.66f * q1k + 0.086f * pow(q1k, 2.0f))), maxqc);
                    ql_ice = sgm[k] * (std::exp(-1.0f) + 0.66f * q1k + 0.086f * pow(q1k, 2));
                }

                // In saturated grid cells, use average of SGS and resolved values
                // if ( qc[k] > 1.e-6 ) ql_water = 0.5 * ( ql_water + qc[k] )
                // ql_ice is actually the total frozen condensate (snow+ice),
                // if ( (qi[k]+qs[k]) > 1.e-9 ) ql_ice = 0.5 * ( ql_ice + (qi[k]+qs[k]) )

                if (cldfra_bl1d[k] < 0.001f) {
                    ql_ice = 0.0f;
                    ql_water = 0.0f;
                    cldfra_bl1d[k] = 0.0f;
                }

                liq_frac = std::min(1.0f, std::max(0.0f, (t - tice) / (tliq - tice)));
                qc_bl1d[k] = liq_frac * ql_water;
                qi_bl1d[k] = (1.0f - liq_frac) * ql_ice;

                // Above tropopausel: eliminate subgrid clouds from CB schemes.
                if (k >= k_tropo) {
                    cldfra_bl1d[k] = 0.0f;
                    qc_bl1d[k] = 0.0f;
                    qi_bl1d[k] = 0.0f;
                }

                // Buoyancy-flux-related calculations follow...

                if((xland-1.5f)>=0.0f) {
                    q1k = std::max(q1[k],-2.5f); // water
                } else {
                    q1k = std::max(q1[k],-2.0f); // land
                }

                // "Fng" represents the non-Gaussian transport factor
                // (non-dimensional) from Bechtold et al. 1995
                // (hereafter BCMT95), section 3(c).
                //
                // Use the form of "Fng" from Bechtold and Siebesma (1998, JAS)
                if(q1k >= 1.0f) {
                    fng = 1.0;
                } else if(q1k >= -1.7f && q1k < 1.0f) {
                    fng = exp(-0.4f*(q1k-1.0f));
                } else if(q1k >= -2.5f && q1k < -1.7f) {
                    fng = 3.0f + float(exp(-3.8f*(q1k+1.7f)));
                } else {
                    fng = std::min(23.9f + float(exp(-1.6f*(q1k+2.5f))),60.f);
                }

                cfmax = std::min(cldfra_bl1d[k],0.6f);
                // Further limit the cf going into vt & vq near the surface
                zsl = std::min(std::max(25.f, 0.1f*pblh2), 100.f);
                wt  = std::min(zagl/zsl, 1.0f); // ==0 at z=0 m, ==1 above ekman layer
                cfmax = cfmax*wt;

                // bb is "b" in BCMT95. Their "b" differs from "b" in CB02
                // (i.e., b[k] above) by a factor of T/theta. Strictly, b[k]
                // above is formulated in terms of sat. mixing ratio, but bb in
                // BCMT95 is cast in terms of sat. specific humidity. The
                // conversion is neglected here.
                bb = b[k]*t/th[k];
                qww = 1.f+0.61f*qw[k];
                alpha = 0.61f*th[k];
                beta = (th[k]/t)*(xl/cp) - 1.61f*th[k];
                vt[k] = qww - cfmax*beta*bb*fng - 1.f;
                vq[k] = alpha + cfmax*beta*a[k]*fng - tv0;

                // Dampen amplification factor where need be
                fac_damp = std::min(zagl * 0.0025f, 1.0f);
                cld_factor = 1.0f + fac_damp * std::min(std::pow(std::max(0.0f, (rh[k] - 0.92f)) / 0.145f, 2.0f), 0.37f);
                cldfra_bl1d[k] = std::min(1.0f, cld_factor * cldfra_bl1d[k]);
                // EQDEBUG
                //		printf("vt[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",vt[k],cfmax,beta,bb,fng);
                //		printf("vq[k] %15.15g %15.15g %15.15g %15.15g %15.15g\n",vq[k],alpha,cfmax,beta,a[k]);
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
// this is the downdraft mass flux scheme - analogus to edmf_jpl but
// flipped updraft to downdraft. this scheme is currently only tested
// for stratocumulus cloud conditions. for a detailed desctiption of the
// model, see paper.
void ddmf_jpl_cc(int& kts, int& kte, float& dt, const float* zw, const float* dz, const float* p,
              const float* u, const float* v, const float* th, const float* thl, const float* thv, 
	      const float* tk,const float* qt, const float* qv, const float* qc, const float* 
	      rho, const float* exner,float& ust, float& wthl, float& wqt, float& pblh, int& kpbl,
              float* edmf_a_dd, float* edmf_w_dd, float* edmf_qt_dd,
              float* edmf_thl_dd, float* edmf_ent_dd, float* edmf_qc_dd,
              float* sd_aw, float* sd_awthl, float* sd_awqt,
              float* sd_awqv, float* sd_awqc, float* sd_awu,
              float* sd_awv, float* sd_awqke,
              const float* qc_bl1d, const float* cldfra_bl1d,
              const float* rthraten, float& svp1, float& grav, float& onethird, float& p1000mb, 
	      float& rcp, float& xlvcp, float& cp, float& rvovrd ) {
    int ndown = 5;
    int dd_initk[ndown];
    float randnum[ndown];
    float downw[kte + 1][ndown];
    float downthl[kte + 1][ndown];
    float downqt[kte + 1][ndown];
    float downqc[kte + 1][ndown];
    float downa[kte + 1][ndown];
    float downu[kte + 1][ndown];
    float downv[kte + 1][ndown];
    float downthv[kte + 1][ndown];
    float ent[kte + 1][ndown];
    int enti[kte + 1][ndown];
    int k, i, ki, kminrad, qltop, p700_ind, qlbase;
    float wthv, wstar, qstar, thstar, sigmaw, sigmaqt, sigmath, z0, pwmin, pwmax, wmin, wmax, wlv, wtv, went, mindownw;
    float b, qtn, thln, thvn, qcn, un, vn, qken, wn2, wn, thvk, pk, entexp, entw, beta_dm, entexp_m, rho_int;
    float jump_thetav, jump_qt, jump_thetal, refthl, refthv, refqt;
    float minrad, zminrad, radflux, f0, wst_rad, wst_dd;
    bool cloudflg;
    float sigq, xl, rsl, cpm, a, mf_cf, diffqt, fng, qww, alpha, beta, bb, f, pt, t, q2p, b9, satvp, rhgrid;
    float wa = 1.0, wb = 1.5, z00 = 100.0, bcoeff = 0.2;
    float l0 = 80, ent0 = 0.2;
    float dp, dl, adn;
    int debug_mf = 0;
    dl = (1000.0f - 500.0f) / ndown;
    pwmin = -3.0;
    pwmax = -1.0;
    for(k=kts;k<=kte+1;k++) {
    for(i=0;i<ndown;i++) {
    downw[k][i] = 0.0f;
    downthl[k][i] = 0.0f;
    downthv[k][i] = 0.0f;
    downqt[k][i] = 0.0f;
    downqc[k][i] = 0.0f;
    downa[k][i] = 0.0f;
    downu[k][i] = 0.0f;
    downv[k][i] = 0.0f;
    }
    sd_aw[k] = 0.0f;
    sd_awthl[k] = 0.0f;
    sd_awqt[k] = 0.0f;
    sd_awqv[k] = 0.0f;
    sd_awqc[k] = 0.0f;
    sd_awu[k] = 0.0f;
    sd_awv[k] = 0.0f;
    sd_awqke[k] = 0.0f;
    }
    for(k=kts;k<=kte;k++) {
    edmf_a_dd[k] = 0.0f;
    edmf_w_dd[k] = 0.0f;
    edmf_qt_dd[k] = 0.0f;
    edmf_thl_dd[k] = 0.0f;
    edmf_ent_dd[k] = 0.0f;
    edmf_qc_dd[k] = 0.0f;
    }
    //from kts+1 to kte+1
    for(k=kts;k<=kte;k++) {
    for(i=0;i<ndown;i++) {
        ent[k][i] = 0.0f;
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
        if (qc[i] > 1.0e-6 && cldfra_bl1d[i] > 0.5f) {
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
        if (radflux < 0.0f) {
            f0 = abs(radflux) + f0;
        }
    }
    f0 = std::max(f0, 1.0f);
    adn = std::min(0.05f + f0 * 0.001f, 0.3f);
    if (cloudflg) {

    p700_ind = 0;
    float min_value = p[0];
    for (int i = kts; i <= kte; ++i) {
	float pval=abs(p[i]-70000.0f);
        if (pval < min_value) {
            p700_ind = i;
        }
    }

        //p700_ind = minloc(abs(p - 70000.0f), 1.0f);
        jump_thetav = thv[p700_ind] - thv[1] - (thv[p700_ind] - thv[qltop + 3]) / (zw[p700_ind] - zw[qltop + 3]) * (zw[p700_ind] - zw[qltop]);
        jump_qt = qc[p700_ind] + qv[p700_ind] - qc[1] - qv[1];
        jump_thetal = thl[p700_ind] - thl[1] - (thl[p700_ind] - thl[qltop + 3]) / (zw[p700_ind] - zw[qltop + 3]) * (zw[p700_ind] - zw[qltop]);
        refthl = thl[qltop];
        refthv = thv[qltop];
        refqt = qt[qltop];
        wst_rad = pow(grav * zw[qltop] * f0 / (refthl * rho[qltop] * cp), 0.333);
        wst_rad = std::max(wst_rad, 0.1f);
        wstar = std::max(0.0, pow(grav / thv[1] * wthv * pblh, onethird));
        went = thv[1] / (grav * jump_thetav * zw[qltop]) * (0.15f * (pow(wstar, 3) + 5 * pow(ust, 3)) + 0.35f * pow(wst_rad, 3));
        qstar = abs(went * jump_qt / wst_rad);
        thstar = f0 / (rho[qltop] * cp * wst_rad) - went * jump_thetav / wst_rad;
        wst_dd = pow(0.15f * (pow(wstar, 3) + 5 * pow(ust, 3)) + 0.35f * pow(wst_rad, 3), 0.333);

        sigmaw = 0.2f * wst_dd;
        sigmaqt = 40 * qstar;
        sigmath = 1.0f * thstar;
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
            downthv[ki][i] = refthv + 0.01f * downw[ki][i] * sigmath / sigmaw;
            downthl[ki][i] = refthl + 0.01f * downw[ki][i] * sigmath / sigmaw;
        
        for (int k = dd_initk[i] - 1; k >= kts + 1; k--) {
            wmin = 0.3f + dp * 0.0005;
            ent[k+1][i] = 0.33f / (std::min(std::max(-1.0f * downw[k + 1][i], wmin), 0.9f) * dp);
            entexp = ent[k+1][i] * dz[k];
            entexp_m = ent[k+1][i] * 0.333f * dz[k];
            qtn = downqt[k + 1][i] * (1.0f - entexp) + qt[k] * entexp;
            thln = downthl[k + 1][i] * (1.0f - entexp) + thl[k] * entexp;
            un = downu[k + 1][i] * (1.0f - entexp) + u[k] * entexp_m;
            vn = downv[k + 1][i] * (1.0f - entexp) + v[k] * entexp_m;
            pk = (p[k - 1] * dz[k] + p[k] * dz[k - 1]) / (dz[k] + dz[k - 1]);
            condensation_edmf_cc(qtn, thln, pk, zw[k], thvn, qcn,p1000mb,rcp,xlvcp,rvovrd);
            thvk = (thv[k - 1] * dz[k] + thv[k] * dz[k - 1]) / (dz[k] + dz[k - 1]);
            b = grav * (thvn / thvk - 1.0f);
            entw = entexp;
            mindownw = std::min(downw[k + 1][i], -0.2f);
            wn = downw[k + 1][i] + (-2.0f * ent[k+1][i] * downw[k + 1][i] - bcoeff * b / mindownw) * std::min(dz[k], 250.0f);
            if (wn < downw[k + 1][i] - std::min(1.25f * dz[k] / 200.0f, -2.0f)) {
                wn = downw[k + 1][i] - std::min(1.25f * dz[k] / 200.0f, -2.0f);
            }
            if (wn > downw[k + 1][i] + std::min(1.25f * dz[k] / 200.0f, 2.0f)) {
                wn = downw[k + 1][i] + std::min(1.25f * dz[k] / 200.0f, 2.0f);
            }
            wn = std::max(std::min(wn, 0.0f), -3.0f);
            if (wn < 0.0f) {
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
                    downw[k][i] = 0.0f;
                    downthv[k][i] = 0.0f;
                    downthl[k][i] = 0.0f;
                    downqt[k][i] = 0.0f;
                    downqc[k][i] = 0.0f;
                    downu[k][i] = 0.0f;
                    downv[k][i] = 0.0f;
                    }
                break;
                }
            }
        }
    }
    for (int i = 0; i < ndown; i++) {
      downw[0][i] = 0.0f;
      downa[0][i] = 0.0f;
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
        if (edmf_a_dd[k] > 0.0f) {
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

// assuming float is equivalent to float or float. adjust as necessary.
void topdown_cloudrad_cc(int& kts, int& kte, const float* dz1, const float* zw, float& fltv, float& xland, int& kpbl, float& pblh, const float* sqc, const float* sqi, const float* sqw, const float* thl, const float* th1, const float* ex1, const float* p1, const float*  rho1, const float* thetav, const float* cldfra_bl1d, const float* rthraten, float& maxkhtopdown, float* khtopdown, float* tkeprodtd) {
    // constants
  /*
    const float pfac = 2.0, zfmin = 0.01, phifac = 8.0;
    const float grav = 9.81, cp = 1004.0, xlv = 2.5e6, xlvcp = xlv / cp, r_d = 287.0, ep_2 = 0.622, p608 = 0.608, karman = 0.4;
    const float twothirds = 2.0f / 3.0, onethird = 1.0f / 3.0;
  */
  //Main meaningfull difference is cp=1004.5 vs cp=1004
    const float pfac = 2.0, zfmin = 0.0099999997764826, phifac = 8.0;
    const float grav = 9.8100004196166992, cp = 1004.5, xlv = 2.5e6, xlvcp = 2488.8002929687500000, r_d = 287.0, ep_2 = 0.6217504143714905, p608 = 0.6083624362945557, karman = 0.4000000059604645;
    const float twothirds = 0.6666666865348816, onethird = 0.3333333432674408;
    // local variables
    float zfac[kte - kts + 1], wscalek2[kte - kts + 1], zfacent[kte - kts + 1];
    float bfx0, wm3, bfxpbl, dthvx, tmp1;
    float temps, templ, zl1, wstar3_2;
    float ent_eff, radsum, radflux, we, rcldb, rvls, minrad = 100., zminrad = pblh;
    int k, kk, kminrad = kpbl;
    bool cloudflg = false;

    maxkhtopdown = 0.0;
    for (kk = kts; kk <= kte; ++kk) {
      khtopdown[kk] = 0.0f;
      tkeprodtd[kk] = 0.0f;
    }

    // check for stratocumulus-topped boundary layers
    for (kk = std::max(0, kpbl + 1 - 2) - 1; kk <= kpbl + 1 + 3 - 1; ++kk) {
        if (sqc[kk - kts] > 1.e-6 || sqi[kk - kts] > 1.e-6 || cldfra_bl1d[kk - kts] > 0.5f) {
            cloudflg = true;
        }
        if (rthraten[kk - kts] < minrad) {
            minrad = rthraten[kk - kts];
            kminrad = kk;
            zminrad = zw[kk - kts] + 0.5f * dz1[kk - kts];
        }
    }

    if (std::max(kminrad, kpbl) < 1) cloudflg = false;

    if (cloudflg) {
        zl1 = dz1[kts];
        k = std::max(kpbl - 1, kminrad - 1);
        templ = thl[k - kts] * ex1[k - kts];
        rvls = 100.f * 6.112f * std::exp(17.67f * (templ - 273.16) / (templ - 29.65)) * (ep_2 / p1[k + 1 - kts]);
        temps = templ + (sqw[k - kts] - rvls) / (cp / xlv + ep_2 * xlv * rvls / (r_d * std::pow(templ, 2)));
        rvls = 100.f * 6.112f * std::exp(17.67f * (temps - 273.15) / (temps - 29.65)) * (ep_2 / p1[k + 1 - kts]);
        rcldb = std::max(sqw[k - kts] - rvls, 0.0f);
        dthvx = (thl[k + 2 - kts] + th1[k + 2 - kts] * p608 * sqw[k + 2 - kts]) - (thl[k - kts] + th1[k - kts] * p608 * sqw[k - kts]);
        dthvx = std::max(dthvx, 0.1f);
        tmp1 = xlvcp * rcldb / (ex1[k - kts] * dthvx);
        ent_eff = 0.2f + 0.2f * 8.f * tmp1;
        radsum = 0.0;
        for (kk = std::max(0, kpbl - 3) ; kk <= kpbl + 1 + 3 - 1; ++kk) {
            radflux = rthraten[kk - kts] * ex1[kk - kts]; // converts theta/s to temp/s
            radflux = radflux * cp / grav * (p1[kk - kts] - p1[kk + 1 - kts]); // converts temp/s to w/m^2
            if (radflux < 0.0f) radsum = std::abs(radflux) + radsum;
        }
        if ((xland - 1.5f) >= 0) { // water
            radsum = std::min(radsum, 90.0f);
            bfx0 = std::max(radsum / rho1[k - kts] / cp, 0.0f);
        } else { // land
            radsum = std::min(0.25f * radsum, 30.0f); // practically turn off over land
            bfx0 = std::max(radsum / rho1[k - kts] / cp - std::max(fltv, 0.0f), 0.0f);
        }
        wm3 = grav / thetav[k - kts] * bfx0 * std::min(pblh, 1500.f); // this is wstar3
        bfxpbl = -ent_eff * bfx0;
        dthvx = std::max(thetav[k + 1 - kts] - thetav[k - kts], 0.1f);
        we = std::max(bfxpbl / dthvx, -std::sqrt(std::pow(wm3, twothirds)));

        for (kk = kts; kk <= kpbl + 3; ++kk) {
            zfac[kk - kts] = std::min(std::max((1.f - (zw[kk + 1 - kts] - zl1) / (zminrad - zl1)), zfmin), 1.0f);
            zfacent[kk - kts] = 10.f * std::max((zminrad - zw[kk + 1 - kts]) / zminrad, 0.0f) * ((1.f - zfac[kk - kts])*(1.f - zfac[kk - kts])*(1.f - zfac[kk - kts]));
            wscalek2[kk - kts] = std::pow((phifac * karman * wm3 * (zfac[kk - kts])), onethird);
            khtopdown[kk - kts] = wscalek2[kk - kts] * karman * (zminrad - zw[kk + 1 - kts]) * ((1.f - zfac[kk - kts]) * (1.f - zfac[kk - kts]) * (1.f - zfac[kk - kts])); // pfac
            khtopdown[kk - kts] = std::max(khtopdown[kk - kts], 0.0f);
            tkeprodtd[kk - kts] = 2.f * ent_eff * wm3 / std::max(pblh, 100.f) * zfacent[kk - kts];
            tkeprodtd[kk - kts] = std::max(tkeprodtd[kk - kts], 0.0f);
        }
    }
    maxkhtopdown = std::numeric_limits<float>::min();
    for (kk = kts; kk <= kte; ++kk) {
      maxkhtopdown = std::max(maxkhtopdown,khtopdown[kk]);
    }
}

void scale_aware_cc(float& dx, float& pbl1, float& psig_bl, float& psig_shcu) {
    float dxdh;
    psig_bl = 1.0f;
    psig_shcu = 1.0f;
    dxdh = std::max(2.5f * dx, 10.0f) / std::min(pbl1, 3000.0f);
    
    // new form to preserve parameterized mixing - only down 5% at dx = 750 m
    psig_bl = ((dxdh * dxdh) + 0.106f * std::pow(dxdh, 0.667f)) / ((dxdh * dxdh) + 0.066f * std::pow(dxdh, 0.667f) + 0.071f);
    
    // assume a 500 m cloud depth for shallow-cu clouds
    dxdh = std::max(2.5f * dx, 10.0f) / std::min(pbl1 + 500.0f, 3500.0f);
    
    // hyeyum hailey shin and song-you hong 2013, tke in entrainment zone
    psig_shcu = ((dxdh * dxdh) + 0.145f * std::pow(dxdh, 0.667f)) / ((dxdh * dxdh) + 0.172f * std::pow(dxdh, 0.667f) + 0.170f);
    
    // clamping psig_bl and psig_shcu to [0, 1]
    psig_bl = std::max(0.0f, std::min(psig_bl, 1.0f));
    psig_shcu = std::max(0.0f, std::min(psig_shcu, 1.0f));
    //    exit(1);
}

// ==================================================================
//>\ingroup gsd_mynn_edmf
// this subroutine calculates hybrid diagnotic boundary-layer height (pblh).
//
// notes on the pblh formulation: the 1.5-theta-increase method defines
//pbl heights as the level at.
//which the potential temperature first exceeds the minimum potential.
//temperature within the boundary layer by 1.5f k. when applied to.
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
void get_pblh_cc(int &kts, int &kte, float &zi, float* thetav1d, float *qke1d, float *zw1d, float* dz1d, float &landsea, int &kzi) {
    // HR: SEGFAULTS WHEN ACCESSING VECTORS, need to look into how to pass 1d arrays to c++ from fortran
    // constants
    const float sbl_lim = 200.0;
    const float sbl_damp = 400.0;

    // local variables
    float pblh_tke, qtke, qtkem1, maxqke, tkeeps, minthv, delt_thv;
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
    delt_thv = (landsea - 1.5f) >= 0 ? 1.0f : 1.25f;
    //for (int k = kthv + 1; k < kte; ++k) {
    for (int k = kts + 1; k < kte; ++k) {
        if (thetav1d[k - kts] >= (minthv + delt_thv)) {
            zi = zw1d[k - kts] - dz1d[k - 1 - kts] * std::min((thetav1d[k - kts] - (minthv + delt_thv)) / std::max(thetav1d[k - kts] - thetav1d[k - 1 - kts], 1.e-6f), 1.0f);
            if (zi > 0) break;
        }
        if (k == kte - 1) zi = zw1d[kts + 1 - kts]; // exit safeguard
    }

    // for stable boundary layers, use tke method to complement the thetav-based definition
    pblh_tke = 0.0;
    maxqke = std::max(qke1d[kts - kts], 0.0f);
    tkeeps = maxqke / 40.0;
    tkeeps = std::max(tkeeps, 0.02f);

    for (int k = kts + 1; k < kte; ++k) {
        qtke = std::max(qke1d[k - kts] / 2.0f, 0.0f);
        qtkem1 = std::max(qke1d[k - 1 - kts] / 2.0f, 0.0f);
        if (qtke <= tkeeps) {
            pblh_tke = zw1d[k - kts] - dz1d[k - 1 - kts] * std::min((tkeeps - qtke) / std::max(qtkem1 - qtke, 1.0e-6f), 1.0f);
            pblh_tke = std::max(pblh_tke, zw1d[kts + 1 - kts]);
            break;
        }
        if (k == kte - 1) pblh_tke = zw1d[kts + 1 - kts]; // exit safeguard
    }

    // limit pblh_tke to not exceed the thetav-based pbl height +/- 350 m
    pblh_tke = std::min(pblh_tke, zi + 350.0f);
    pblh_tke = std::max(pblh_tke, std::max(zi - 350.0f, 10.0f));

    float wt = 0.5f * std::tanh((zi - sbl_lim) / sbl_damp) + 0.5;
    if (maxqke > 0.05) {
        zi = pblh_tke * (1.0f - wt) + zi * wt;
    }

    // compute kpbl (kzi)
    for (int k = kts + 1; k < kte; ++k) {
        if (zw1d[k - kts] >= zi) {
            kzi = k - 1;
            break;
        }
    }
}

void retrieve_exchange_coeffs_cc(int& kts, int& kte, float* dfm, float* dfh, const float* dz, float* k_m, float* k_h) {
    float dzk;
    k_m[kts] = 0.0;
    k_h[kts] = 0.0;
    for (int k = kts + 1; k <= kte; ++k) {
        dzk = 0.5f * (dz[k] + dz[k - 1]);
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
//           sigmal levels (or interfaces beween the model layers).
//
//>\ingroup gsd_mynn_edmf
// this subroutine calculates the mixing lengths.

void mym_length(int kts, int kte, float xland, float* dz, float* dx, float* zw, float rmo, float flt, float fltv, float flq, float* vt, float* vq, float* u1, float* v1, float* qke, float* dtv, float* el, float zi, float* theta, float* qkw, float psig_bl, float* cldfra_bl1d, int bl_mynn_mixlength, float* edmf_w1, float* edmf_a1, float grav, float karman, float tv0, float gtr) {
    int i, j, k;
    float elt, vsc;
    float qtke[kte+1], elblmin[kte+1], elblavg[kte+1], thetaw[kte+1];
    float wt, wt2, zi2, h1, h2, hs, elblmin0, elblavg0, cldavg;
    float cns, alp1, alp2, alp3, alp4, alp5, alp6;
    float minzi = 300.0;
    float maxdz = 750.0;
    float mindz = 300.0;
    float zslh = 100.0;
    float csl = 2.0;
    float qke_elb_min = 0.018f;
    float afk, abk, zwk, zwk1, dzk, qdz, vflx, bv, tau_cloud, wstar, elb, els, elf, el_stab, el_mf, el_stab_mf, elb_mf, pblh_plus_ent, uonset, ugrid, wt_u, el_les;
    float ctau = 1000.0;

    switch(bl_mynn_mixlength) {
        case 0:
            cns = 2.7;
            alp1 = 0.23;
            alp2 = 1.0;
            alp3 = 5.0;
            alp4 = 100.0;
            alp5 = 0.3;
            zi2 = std::min(10000.0f, zw[kte-2]);
            h1 = std::max(0.3f * zi2, mindz);
            h1 = std::min(h1, maxdz);
            h2 = h1 / 2.0;
            qkw[kts] = std::sqrt(std::max(qke[kts], 1.0e-10f));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0f - afk;
                qkw[k] = std::sqrt(std::max(qke[k] * abk + qke[k-1] * afk, 1.0e-3f));
            }
            elt = 1.0e-5;
            vsc = 1.0e-5;
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi2 + h1) {
                dzk = 0.5f * (dz[k] + dz[k-1]);
                qdz = std::max(qkw[k] - qmin, 0.03f) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }
            elt = alp1 * elt / vsc;
            vflx = (vt[kts] + 1.0f) * flt + (vq[kts] + tv0) * flq;
            vsc = std::pow(gtr * elt * std::max(vflx, 0.0f), 1.0f / 3.0f);
            el[kts] = 0.0;
            zwk1 = zw[kts+1];
            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];
                if (dtv[k] > 0.0f) {
                    bv = std::sqrt(gtr * dtv[k]);
                    elb = alp2 * qkw[k] / bv * (1.0f + alp3 / alp2 * std::sqrt(vsc / (bv * elt)));
                    elf = alp2 * qkw[k] / bv;
                } else {
                    elb = 1.0e10;
                    elf = elb;
                }
                if (rmo > 0.0f) {
                    els = karman * zwk / (1.0f + cns * std::min(zwk * rmo, zmax));
                } else {
                    els = karman * zwk * std::pow(1.0f - alp4 * zwk * rmo, 0.2f);
                }
                wt = 0.5f * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5;
                el[k] = std::min(elb / (elb / elt + elb / els + 1.0f), elf);
            }
            break;
        case 1:
            ugrid = std::sqrt(u1[kts] * u1[kts] + v1[kts] * v1[kts]);
            uonset = 15.0;
            wt_u = (1.0f - std::min(std::max(ugrid - uonset, 0.0f) / 30.0f, 0.5f));
            cns = 2.7;
            alp1 = 0.23;
            alp2 = 0.3;
            alp3 = 2.5f * wt_u;
            alp4 = 5.0;
            alp5 = 0.3;
            alp6 = 50.0;
            zi2 = std::max(zi, 300.0f);
            h1 = std::max(0.3f * zi2, 300.0f);
            h1 = std::min(h1, 600.0f);
            h2 = h1 / 2.0;
            qkw[kts] = std::sqrt(std::max(qke[kts], 1.0e-10f));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0f - afk;
                qkw[k] = std::sqrt(std::max(qke[k] * abk + qke[k-1] * afk, 1.0e-3f));
                qtke[k] = 0.5f * (qkw[k] * qkw[k]);
                thetaw[k] = theta[k] * abk + theta[k-1] * afk;
            }
            elt = 1.0e-5;
            vsc = 1.0e-5;
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi2 + h1) {
                dzk = 0.5f * (dz[k] + dz[k-1]);
                qdz = std::min(std::max(qkw[k] - qmin, 0.01f), 30.0f) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }
            elt = std::min(std::max(alp1 * elt / vsc, 10.0f), 400.0f);
            vflx = fltv;
            vsc = std::pow(gtr * elt * std::max(vflx, 0.0f), 1.0f / 3.0f);
            el[kts] = 0.0;
            zwk1 = zw[kts+1];
            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];
                if (dtv[k] > 0.0f) {
                    bv = std::max(std::sqrt(gtr * dtv[k]), 0.0001f);
                    elb = std::max(alp2 * qkw[k], alp6 * edmf_a1[k-1] * edmf_w1[k-1]) / bv * (1.0f + alp3 * std::sqrt(vsc / (bv * elt)));
                    elb = std::min(elb, zwk);
                    elf = 1.0f * qkw[k] / bv;
                    elblavg[k] = std::max(elblavg[k], alp6 * edmf_a1[k-1] * edmf_w1[k-1] / bv);
                } else {
                    elb = 1.0e10;
                    elf = elb;
                }
                if (rmo > 0.0f) {
                    els = karman * zwk / (1.0f + cns * std::min(zwk * rmo, zmax));
                } else {
                    els = karman * zwk * std::pow(1.0f - alp4 * zwk * rmo, 0.2f);
                }
                wt = 0.5f * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5;
                el[k] = std::min(elb / (elb / elt + elb / els + 1.0f), elf);
                el[k] = el[k] * psig_bl;
            }
            break;
        case 2:
            uonset = 3.5f + dz[kts] * 0.1;
            ugrid = std::sqrt(u1[kts] * u1[kts] + v1[kts] * v1[kts]);
            cns = 3.5;
            alp1 = 0.22;
            alp2 = 0.30;
            alp3 = 2.0;
            alp4 = 5.0;
            alp5 = alp2;
            alp6 = 50.0;
            zi2 = std::max(zi, minzi);
            h1 = std::max(0.3f * zi2, mindz);
            h1 = std::min(h1, maxdz);
            h2 = h1 * 0.5;
            qtke[kts] = std::max(0.5f * qke[kts], 0.5f*qkemin);
            qkw[kts] = std::sqrt(std::max(qke[kts], qkemin));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0f - afk;
                qkw[k] = std::sqrt(std::max(qke[k] * abk + qke[k-1] * afk, qkemin));
                qtke[k] = 0.5f * qkw[k] * qkw[k];
                thetaw[k] = theta[k] * abk + theta[k-1] * afk;
            }
            elt = 1.0e-5;
            vsc = 1.0e-5;
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi2 + h1) {
                dzk = 0.5f * (dz[k] + dz[k-1]);
                qdz = std::min(std::max(qkw[k] - qmin, 0.03f), 30.0f) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }
            elt = std::min(std::max(alp1 * elt / vsc, 10.0f), 400.0f);
            vflx = fltv;
            vsc = std::pow(gtr * elt * std::max(vflx, 0.0f), 1.0f / 3.0f);
            el[kts] = 0.0;
            zwk1 = zw[kts+1];
            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];
                dzk = 0.5f * (dz[k] + dz[k-1]);
                cldavg = 0.5f * (cldfra_bl1d[k-1] + cldfra_bl1d[k]);
                if (dtv[k] > 0.0f) {
                    bv = std::max(std::sqrt(gtr * dtv[k]), 0.001f);
                    elb_mf = std::max(alp2 * qkw[k], alp6 * edmf_a1[k-1] * edmf_w1[k-1]) / bv * (1.0f + alp3 * std::sqrt(vsc / (bv * elt)));
                    elb = std::min(std::max(alp5 * qkw[k], alp6 * edmf_a1[k] * edmf_w1[k]) / bv, zwk);
                    wstar = 1.25f * std::pow(gtr * zi * std::max(vflx, 1.0e-4f), 1.0f / 3.0f);
                    tau_cloud = std::min(std::max(ctau * wstar / grav, 30.0f), 150.0f);
                    wt = 0.5f * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5;
                    tau_cloud = tau_cloud * (1.0f - wt) + 50.0f * wt;
                    elf = std::min(std::max(tau_cloud * std::sqrt(std::min(qtke[k], 40.0f)), alp6 * edmf_a1[k] * edmf_w1[k] / bv), zwk);
                } else {
                    wstar = 1.25f * std::pow(gtr * zi * std::max(vflx, 1.0e-4f), 1.0f / 3.0f);
                    tau_cloud = std::min(std::max(ctau * wstar / grav, 50.0f), 200.0f);
                    wt = 0.5f * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5;
                    tau_cloud = tau_cloud * (1.0f - wt) + std::max(100.0, dzk * 0.25) * wt;
                    elb = std::min(tau_cloud * std::sqrt(std::min(qtke[k], 40.0f)), zwk);
                    elf = elb;
                    elb_mf = elb;
                }
                elf = elf / (1.0f + (elf / 800.0f));
                elb_mf = std::max(elb_mf, 0.01f);
                if (rmo > 0.0f) {
                    els = karman * zwk / (1.0f + cns * std::min(zwk * rmo, zmax));
                } else {
                    els = karman * zwk * std::pow(1.0f - alp4 * zwk * rmo, 0.2f);
                }
                wt = 0.5f * std::tanh((zwk - (zi2 + h1)) / h2) + 0.5;
                el[k] = std::sqrt(els * els / (1.0f + (els * els / elt * elt) + (els * els / elb_mf * elb_mf)));
                el[k] = el[k] * (1.0f - wt) + elf * wt;
                el_les = std::min(els / (1.0f + (els / 12.0f)), elb_mf);
                el[k] = el[k] * psig_bl + (1.0f - psig_bl) * el_les;
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
// namelist paramter \p bl_mynn_edmf is set to 1 (recommended).
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

void dmp_mf_cc(const int& kts,const int& kte, float& dt, float* zw, float* dz, float* p, float* rho, int& momentum_opt, int& tke_opt, int& scalar_opt, float* u, float* v, float* w, float* th, float* thl, float* thv, float* tk, float* qt, float* qv, float* qc, float* qke, float* qnc, float* qni, float* qnwfa, float* qnifa, float* qnbca, float& ust, float& flt, float& fltv, float& flq, float& flqv, float& pblh, int& kpbl, float& dx, float& landsea, float& ts, float* edmf_a, float* edmf_w, float* edmf_qt, float* edmf_thl, float* edmf_ent, float* edmf_qc, float* s_aw, float* s_awthl, float* s_awqt, float* s_awqv, float* s_awqc, float* s_awu, float* s_awv, float* s_awqke, float* s_awqnc, float* s_awqni, float* s_awqnwfa, float* s_awqnifa, float* s_awqnbca, int& nchem, float** chem1, float** s_awchem, bool& mix_chem, float* qc_bl1d, float* cldfra_bl1d, float* qc_bl1d_old, float* cldfra_bl1d_old, float& psig_shcu, float& maxwidth, int& ktop, float& maxmf, float& ztop, float* rstoch_col, float grav, float gtr, float p608) {
    int nup = 8;
    int debug_mf = 0;
    float nup2;
    float upw[kte+1+1][nup], upthl[kte+1+1][nup], upqt[kte+1+1][nup], upqc[kte+1+1][nup], upqv[kte+1+1][nup], upa[kte+1+1][nup], upu[kte+1+1][nup], upv[kte+1+1][nup], upthv[kte+1+1][nup], upqke[kte+1+1][nup], upqnc[kte+1+1][nup], upqni[kte+1+1][nup], upqnwfa[kte+1+1][nup], upqnifa[kte+1+1][nup], upqnbca[kte+1+1][nup];
    float ent[kte+1][nup];
    int enti[kte+1][nup];
    int k, i, k50;
    float fltv2, wstar, qstar, thstar, sigmaw, sigmaqt, sigmath, z0, pwmin, pwmax, wmin, wmax, wlv, psig_w, maxw, wpbl;
    float b, qtn, thln, thvn, qcn, un, vn, qken, qncn, qnin, qnwfan, qnifan, qnbcan, wn2, wn, entexp, entexm, entw, bcoeff, thvkm1, thvk, pk, rho_int;
    float wa = 0.6666666865348816; //2./3.
    float wb = 0.0020000000949949, wc = 1.5;
    float l0 = 100., ent0 =  0.1000000014901161;
    float atot = 0.1000000014901161;
    float lmax = 1000.;
    float lmin = 300.;
    float dlmin = 0.;
    float minwidth;
    float dl;
    float dcut = 1.2000000476837158;
    float d;
    float cn, c, l, n, an2, hux, wspd_pbl, cloud_base, width_flx;
    float chemn[nchem];
    float upchem[kte+1+1][nup][nchem];
    float ic;
    float edmf_chem[kte+1+1][nchem];
    float envm_u[kte+1+1],envm_v[kte+1+1],envm_sqc[kte+1+1],envm_thl[kte+1+1],envm_sqv[kte+1+1];
    bool superadiabatic;
    float sigq, xl, rsl, cpm, a, qmq, mf_cf, aup, q1, diffqt, qsat_tk, fng, qww, alpha, beta, bb, f, pt, t, q2p, b9, satvp, rhgrid, ac_mf, ac_strat, qc_mf;
    float cf_thresh = 0.5;
    float exneri[kte+1], dzi[kte+1], rhoz[kte+1];
    float thp, qtp, qcp, qcs, esat, qsl;
    float csigma, acfac, ac_wsp;
    int overshoot;
    float bvf, frz, dzp;
    float adjustment, flx1, flt2;
    float fluxportion = 0.75;
    float sublim, qc_ent, qv_ent, qt_ent, thl_ent, detrate, detrateuv, oow, exc_fac, aratio, detturb, qc_grid, qc_sgs, exc_heat, exc_moist, tk_int, tvs, qc_plume;
    float cdet = 0.0222222227603197;//1./45.;
    float dzpmax = 300.;
    float csub = 0.25;
    float pgfac = 0.00;
    float uk, ukm1, vk, vkm1, dxsa;

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
        if (zw[k] > pblh + 500.f) {
            break;
        }
        wpbl = w[k];
        if (w[k] < 0.f) {
            wpbl = 2.*w[k];
        }
        maxw = std::max(maxw, float(abs(wpbl)));
        if (zw[k] <= 50.f) {
            k50 = k;
        }
        float qc_sgs = std::max(qc[k], qc_bl1d[k]);
        if (qc_sgs > 1e-5 && cldfra_bl1d[k] >= 0.5f && cloud_base == 9000.0f) {
            cloud_base = 0.5f*(zw[k]+zw[k+1]);
        }
    }
    maxw = std::max(0.f, maxw - 1.0f);
    psig_w = std::max(0.0f, 1.0f - maxw);
    psig_w = std::min(psig_w, psig_shcu);
    fltv2 = fltv;
    if (psig_w == 0.0f && fltv > 0.0f) {
        fltv2 = -1.*fltv;
    }
    superadiabatic = false;
    tvs = ts*(1.0+p608*qv[kts]);
    for (int k = 0; k < std::max(1, k50-1); k++) {
        if (k == 0) {
            if ((thv[k]-tvs)/(0.5f*dz[k]) < hux) {
                superadiabatic = true;
            } else {
                superadiabatic = false;
                break;
            }
        } else {
            if ((thv[k]-thv[k-1])/(0.5f*(dz[k]+dz[k-1])) < hux) {
                superadiabatic = true;
            } else {
                superadiabatic = false;
                break;
            }
        }
    }
    maxwidth = std::min(dx*dcut, lmax);
    maxwidth = std::min(maxwidth, 1.1f*pblh);
    if (landsea-1.5f < 0) {
        maxwidth = std::min(maxwidth, 0.5f*cloud_base);
    } else {
        maxwidth = std::min(maxwidth, 0.9f*cloud_base);
    }
    wspd_pbl = sqrt(std::max(u[kts]*u[kts] + v[kts]*v[kts], 0.01f));
    if (landsea-1.5f < 0) {
        width_flx = std::max(std::min(1000.f*(0.6f*float(tanh((fltv - 0.040f)/0.04f)) + .5f),1000.f), 0.f);
    } else {
        width_flx = std::max(std::min(1000.f*(0.6f*float(tanh((fltv - 0.007f)/0.02f)) + .5f),1000.f), 0.f);
    }
    maxwidth = std::min(maxwidth, width_flx);
    minwidth = lmin;
    if (maxwidth >= (lmax - 1.0f) && fltv > 0.2f) {
        minwidth = lmin + dlmin*std::min((fltv-0.2f)/0.3f, 1.f);
    }
    if (maxwidth <= minwidth) {
        nup2 = 0;
        maxwidth = 0.0;
    }
    ktop = 0;
    ztop = 0.0;
    maxmf = 0.0;
    if (fltv2 > 0.002 && maxwidth > minwidth && superadiabatic) {
        float cn = 0.;
        float d = -1.9;
        float dl = (maxwidth - minwidth)/float(nup-1);
        for (int i = 0; i < nup; i++) {
            float l = minwidth + dl*float(i);
            cn = cn + l*l*l * (l*l)/(dx*dx) * dl;
        }
        float c = atot/cn;
        float acfac;
        acfac = 0.5f*tanh((fltv2 - 0.02f)/0.05f) + 0.5f;

        float ac_wsp;
        if (wspd_pbl <= 10.f) {
            ac_wsp = 1.0;
        } else {
            ac_wsp = 1.0f - std::min((wspd_pbl - 10.0f)/15.f, 1.0f);
        }
        acfac = acfac * ac_wsp;
        float an2 = 0.;
        for (int i = 0; i < nup; i++) {
            float l = minwidth + dl*float(i);
            float n = c*l*l*l * (l*l)/(dx*dx) * dl;
            upa[kts][i] = n*l*l/(dx*dx) * dl;
            upa[kts][i] = upa[kts][i]*acfac;
            an2 = an2 + upa[kts][i];
        }
        float z0 = 50.;
        float pwmin = 0.1;
        float pwmax = 0.4;
        float wstar = std::max(1.e-2, pow(gtr*fltv2*pblh, 1./3.f));
        float qstar = std::max(flq, 1.0e-5f)/wstar;
        float thstar = flt/wstar;
        float csigma;
        if (landsea-1.5f >= 0) {
            csigma = 1.34;
        } else {
            csigma = 1.34;
        }
        float exc_fac;
        if (env_subs) {
            exc_fac = 0.0;
        } else {
            if (landsea-1.5f >= 0) {
                exc_fac = 0.58*4.0;
            } else {
                exc_fac = 0.58;
            }
        }
        exc_fac = exc_fac * ac_wsp;
        float sigmaw = csigma*wstar*pow(z0/pblh, 1./3.f)*(1 - 0.8*z0/pblh);
        float sigmaqt = csigma*qstar*pow(z0/pblh, 1./3.f);
        float sigmath = csigma*thstar*pow(z0/pblh, 1./3.f);
        float wmin = std::min(sigmaw*pwmin, 0.1f);
        float wmax = std::min(sigmaw*pwmax, 0.5f);
        for (int i = 0; i < nup; i++) {
            float wlv = wmin+(wmax-wmin)/nup2*float(i);
            upw[kts][i] = wmin + float(i+1)/float(nup)*(wmax-wmin);
            upu[kts][i] = (u[kts]*dz[kts+1]+u[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]);
            upv[kts][i] = (v[kts]*dz[kts+1]+v[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]);
            upqc[kts][i] = 0.0;
            float exc_heat = exc_fac*upw[kts][i]*sigmath/sigmaw;
            upthv[kts][i] = (thv[kts]*dz[kts+1]+thv[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]) + exc_heat;
            upthl[kts][i] = (thl[kts]*dz[kts+1]+thl[kts+1]*dz[kts])/(dz[kts]+dz[kts+1]) + exc_heat;
            float exc_moist = exc_fac*upw[kts][i]*sigmaqt/sigmaw;
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
        dxsa = 1.f - std::min(std::max((12000.0f-dx)/(12000.0f-3000.0f), 0.f), 1.f);
        for (int i = 0; i < nup; i++) {
            float qc_ent = 0.0;
            int overshoot = 0;
            float l = minwidth + dl*float(i);
            for (int k = kts+1; k < kte-1; k++) {
                float wmin = 0.3f + l*0.0005;
                ent[k][i] = 0.33/(std::min(std::max(upw[k-1][i], wmin), 0.9f)*l);
                ent[k][i] = std::max(ent[k][i], 0.0003f);
                if (zw[k] >= std::min(pblh+1500.f, 4000.f)) {
                    ent[k][i] = ent[k][i] + (zw[k]-std::min(pblh+1500.f,4000.f))*5.0e-6f;
                }
                ent[k][i] = ent[k][i] * (1.0f - rstoch_col[k]);
                ent[k][i] = std::min(ent[k][i], 0.9f/(zw[k+1]-zw[k]));
                float uk = (u[k]*dz[k+1]+u[k+1]*dz[k])/(dz[k+1]+dz[k]);
                float ukm1 = (u[k-1]*dz[k]+u[k]*dz[k-1])/(dz[k-1]+dz[k]);
                float vk = (v[k]*dz[k+1]+v[k+1]*dz[k])/(dz[k+1]+dz[k]);
                float vkm1 = (v[k-1]*dz[k]+v[k]*dz[k-1])/(dz[k-1]+dz[k]);
                float entexp = ent[k][i]*(zw[k+1]-zw[k]);
                float entexm = entexp*0.3333;
                float qtn = upqt[k-1][i]*(1.-entexp) + qt[k]*entexp;
                float thln = upthl[k-1][i]*(1.-entexp) + thl[k]*entexp;
                float un = upu[k-1][i]*(1.-entexm) + u[k]*entexm + dxsa*pgfac*(uk - ukm1);
                float vn = upv[k-1][i]*(1.-entexm) + v[k]*entexm + dxsa*pgfac*(vk - vkm1);
                float qken = upqke[k-1][i]*(1.-entexp) + qke[k]*entexp;
                float qncn = upqnc[k-1][i]*(1.-entexp) + qnc[k]*entexp;
                float qnin = upqni[k-1][i]*(1.-entexp) + qni[k]*entexp;
                float qnwfan = upqnwfa[k-1][i]*(1.-entexp) + qnwfa[k]*entexp;
                float qnifan = upqnifa[k-1][i]*(1.-entexp) + qnifa[k]*entexp;
                float qnbcan = upqnbca[k-1][i]*(1.-entexp) + qnbca[k]*entexp;
                float qc_ent = qcn;
                float qt_ent = qtn;
                float thl_ent = thln;
                if (mix_chem) {
                    for (int ic = 0; ic < nchem; ic++) {
                        chemn[ic] = upchem[k-1][i][ic]*(1.-entexp) + chem1[k][ic]*entexp;
                    }
                }
                float pk = (p[k]*dz[k+1]+p[k+1]*dz[k])/(dz[k+1]+dz[k]);
//                condensation_edmf_cc(qtn, thln, pk, zw[k+1], thvn, qcn);
                float thvk = (thv[k]*dz[k+1]+thv[k+1]*dz[k])/(dz[k+1]+dz[k]);
                float thvkm1 = (thv[k-1]*dz[k]+thv[k]*dz[k-1])/(dz[k-1]+dz[k]);
                float b = grav*(thvn/thvk - 1.0f);
                if (b > 0.f) {
                    bcoeff = 0.15;
                } else {
                    bcoeff = 0.2;
                }
                if (upw[k-1][i] < 0.2f) {
                    wn = upw[k-1][i] + (-2.f * ent[k][i] * upw[k-1][i] + bcoeff*b / std::max(upw[k-1][i], 0.2f)) * std::min(zw[k]-zw[k-1], 250.f);
                } else {
                    wn = upw[k-1][i] + (-2.f * ent[k][i] * upw[k-1][i] + bcoeff*b / upw[k-1][i]) * std::min(zw[k]-zw[k-1], 250.f);
                }
                if (wn > upw[k-1][i] + std::min(1.25f*(zw[k]-zw[k-1])/200.f, 2.0f)) {
                    wn = upw[k-1][i] + std::min(1.25f*(zw[k]-zw[k-1])/200.f, 2.0f);
                }
                if (wn < upw[k-1][i] - std::min(1.25f*(zw[k]-zw[k-1])/200.f, 2.0f)) {
                    wn = upw[k-1][i] - std::min(1.25f*(zw[k]-zw[k-1])/200.f, 2.0f);
                }
                wn = std::min(std::max(wn, 0.0f), 3.0f);
                if (k == kts+1 && wn == 0.f) {
                    nup2 = 0;
                    break;
                }
                if (debug_mf == 1) {
                    if (wn >= 3.0f) {
			std::cout << "**** suspiciously large w:" << std::endl;
			std::cout << "qcn: " << qcn << " ent: " << ent[k][i] << " nup2: " << nup2 << std::endl;
			std::cout << "pblh: " << pblh << " wn: " << wn << " upw(k-1): " << upw[k-1][i] << std::endl;
			std::cout << "k: " << k << " b: " << b << " dz: " << zw[k]-zw[k-1] << std::endl;
                    }
                }
                if (wn <= 0.0f && overshoot == 0) {
                    overshoot = 1;
                    if (thvk-thvkm1 > 0.0f) {
                        float bvf = sqrt(gtr*(thvk-thvkm1)/dz[k]);
                        float frz = upw[k-1][i]/(bvf*dz[k]);
                        dzp = dz[k]*std::max(std::min(frz, 1.0f), 0.0f);
                    }
                } else {
                    dzp = dz[k];
                }
                float aratio = std::min(upa[k-1][i]/(1.f-upa[k-1][i]), 0.5f);
                float detturb = 0.00008;
                float oow = -0.060/std::max(1.0f, (0.5f*(wn+upw[k-1][i])));
                float detrate = std::min(std::max(oow*(wn-upw[k-1][i])/dz[k], detturb), 0.0002f);
                float detrateuv = std::min(std::max(oow*(wn-upw[k-1][i])/dz[k], detturb), 0.0001f);
                envm_thl[k-kts] = envm_thl[k-kts] + (0.5f*(thl_ent + upthl[k-1][i]) - thl[k])*detrate*aratio*std::min(dzp, dzpmax);
                float qv_ent = 0.5f*(std::max(qt_ent-qc_ent, 0.0f) + std::max(upqt[k-1][i]-upqc[k-1][i], 0.0f));
                envm_sqv[k-kts] = envm_sqv[k] + (qv_ent-qv[k])*detrate*aratio*std::min(dzp, dzpmax);
                if (upqc[k-1][i] > 1e-8) {
                    float qc_grid;
                    if (qc[k] > 1e-6) {
                        qc_grid = qc[k];
                    } else {
                        qc_grid = cldfra_bl1d[k]*qc_bl1d[k];
                    }
                    envm_sqc[k-kts] = envm_sqc[k-kts] + std::max(upa[k-1][i]*0.5f*(qcn + upqc[k-1][i]) - qc_grid, 0.0f)*detrate*aratio*std::min(dzp, dzpmax);
                }
                envm_u[k] = envm_u[k] + (0.5f*(un + upu[k-1][i]) - u[k])*detrateuv*aratio*std::min(dzp, dzpmax);
                envm_v[k] = envm_v[k] + (0.5f*(vn + upv[k-1][i]) - v[k])*detrateuv*aratio*std::min(dzp, dzpmax);
                if (wn > 0.f) {
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
            if (upw[ii][i] > 10.0f || upa[ii][i] < 0.0f || upa[ii][i] > atot || nup2 > 10)
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
    if (s_aw[kts+1] != 0.0f) {
        dzi[kts] = 0.5f*(dz[kts] + dz[kts+1]);
        flx1 = std::max(s_aw[kts+1]*(th[kts] - th[kts+1])/dzi[kts], 1.0e-5f);
    } else {
        flx1 = 0.0;
    }
    adjustment = 1.0;
    flt2=std::max(flt,0.0f);
    if (flx1 > fluxportion*flt2/dz[kts] && flx1 > 0.0f) {
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
        if (edmf_a[k] > 0.0f) {
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
            if (edmf_a[k] > 0.0f) {
                for (int ic = 0; ic < nchem; ic++) {
                    edmf_chem[k][ic] /= edmf_a[k];
                }
            }
        }
    }
    if (ktop > 0) {
	float maxqc=0;
        for (int ii=0;ii > ktop; ii++){
		if (edmf_qc[ii]>maxqc){
		    maxqc = edmf_qc[ii];
		}
	}
        if (maxqc < 1.0e-8) {
            maxmf = -1.0*maxmf;
        }
    }
    if (edmf_w[0] > 4.0f) {
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
//       closure        : closure level (2.5, 2.6, or 3.0f)
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
// - The stability functions for level 2.5f or level 3.0f are calculated.
// - If level 3.0f is used, counter-gradient terms are calculated.
// - Production terms of TKE,\f$\theta^{'2}\f$,\f$q^{'2}\f$, and \f$\theta^{'}q^{'}\f$
// are calculated.
// - Eddy diffusivity \f$K_h\f$ and eddy viscosity \f$K_m\f$ are calculated.
// - TKE budget terms are calculated (if the namelist parameter \p tke_budget
// is set to True)
void mym_turbulence_cc(int& kts, int& kte, float& xland, float& closure, float* dz, float* dx, float* zw, float* u, float* v, float* thl, float* thetav, float* ql, float* qw, float* qke, float* tsq, float* qsq, float* cov, float* vt, float* vq, float& rmo, float& flt, float& fltv, float& flq, float& zi, float* theta, float* sh, float* sm, float* el, float* dfm, float* dfh, float* dfq, float* tcd, float* qcd, float* pdk, float* pdt, float* pdq, float* pdc, float* qWT1D, float* qSHEAR1D, float* qBUOY1D, float* qDISS1D, int& tke_budget, float& Psig_bl, float& Psig_shcu, float* cldfra_bl1D, int& bl_mynn_mixlength, float* edmf_w1, float* edmf_a1, float* TKEprodTD, int& spp_pbl, float* rstoch_col, int& debug_code, float& gtr, float& tv0) {
    float q3sq_old, dlsq1, qWTP_old, qWTP_new;
    float dudz, dvdz, dTdz, upwp, vpwp, Tpwp;
    float e6c, dzk, afk, abk, vtt, vqq, cw25, clow, cupp, gamt, gamq, smd, gamv, elq, elh;
    float cldavg;
    float a2fac, duz, ri;
    float auh, aum, adh, adm, aeh, aem, Req, Rsl, Rsl2, gmelq, sm20, sh20, sm25max, sh25max, sm25min, sh25min, sm_pbl, sh_pbl, zi2, wt, slht, wtpr;
    double q2sq, t2sq, r2sq, c2sq, elsq, gmel, ghel, q3sq, t3sq, r3sq, c3sq, dlsq, qdiv, e1, e2, e3, e4, enumc, eden, wden;
    float Prnum, shb;
    const float Prlimit = 5.0;
    float dtv[kte-kts];
    float gm[kte-kts];
    float gh[kte-kts];
    float dqw[kte-kts];
    float dtl[kte-kts];
    float qkw[kte-kts];

    mym_level2_cc(kts, kte, dz, u, v, thl, thetav, qw, ql, vt, vq, dtl, dqw, dtv, gm, gh, sm, sh, tv0, gtr);

    mym_length_cc(kts, kte, xland, dz, zw, rmo, flt, fltv, flq, vt, vq, u, v, qke, dtv, el, zi, theta, qkw, Psig_bl, cldfra_bl1D, bl_mynn_mixlength, edmf_w1, edmf_a1, tv0, gtr);

    for (int k = kts + 1; k <= kte; k++) {
        dzk = 0.5f * (dz[k] + dz[k - 1]);
        afk = dz[k] / (dz[k] + dz[k - 1]);
        abk = 1.0f - afk;
        elsq = el[k] * el[k];
        q3sq = qkw[k] * qkw[k];
        q2sq = b1 * elsq * (sm[k] * gm[k] + sh[k] * gh[k]);

        sh20 = std::max(sh[k], 1e-5f);
        sm20 = std::max(sm[k], 1e-5f);
        sh[k] = std::max(sh[k], 1e-5f);
	//	printf("sh[k] %d %g\n",k,sh[k]);

        // Canuto/Kitamura mod
        duz = (u[k] - u[k - 1]) * (u[k] - u[k - 1]) + (v[k] - v[k - 1]) * (v[k] - v[k - 1]);
        duz = duz / (dzk * dzk);

        // ** Gradient Richardson number **
        ri = -gh[k] / std::max(duz, 1.0e-10f);
        if (ckmod == 1) {
            a2fac = 1.0f / (1.0f + std::max(ri, 0.0f));
        } else {
            a2fac = 1.0;
        }

        Prnum = std::min(0.76f + 4.0f * std::max(ri, 0.0f), Prlimit);

        gmel = gm[k] * elsq;
        ghel = gh[k] * elsq;

        if (debug_code) {
            if (sh[k] < 0.0f || sm[k] < 0.0f) {
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
            qdiv = sqrt(q3sq / q2sq); // HL89: (1-alfa)

            // Use level 2.0 functions as in original MYNN
            sh[k] = sh[k] * qdiv;
            sm[k] = sm[k] * qdiv;
	    //	    printf("sh[k] %d %g %g\n",k,sh[k],qdiv);

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
            sm[k] = q3sq * a1 * (e3 - 3.0f * c1 * e4) / eden;
            sh[k] = q3sq * (a2 * a2fac) * (e2 + 3.0f * c1 * e5c * gmel) / eden;
	    //	    printf("sh[k] %d %g %g %g %g %g %g %g %g %g\n",k,sh[k],q3sq,a2,a2fac,e2,c1,e5c,gmel,eden);
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
                std::cout << " ri=" << ri << " Pr=" << sm[k] / std::max(sh[k], 1e-8f) << std::endl;
                std::cout << " gm=" << gm[k] << " gh=" << gh[k] << std::endl;
                std::cout << " q2sq=" << q2sq << " q3sq=" << q3sq << " q3/q2=" << q3sq / q2sq << std::endl;
                std::cout << " qke=" << qke[k] << " el=" << el[k] << std::endl;
                std::cout << " PBLH=" << zi << " u=" << u[k] << " v=" << v[k] << std::endl;
                std::cout << " SMnum=" << q3sq * a1 * (e3 - 3.0f * c1 * e4) << " SMdenom=" << eden << std::endl;
                std::cout << " SHnum=" << q3sq * (a2 * a2fac) * (e2 + 3.0f * c1 * e5c * gmel) << " SHdenom=" << eden << std::endl;
            }
        }

        // Enforce constraints for level 2.5 functions
        if (sh[k] > sh25max) sh[k] = sh25max;
        if (sh[k] < sh25min) sh[k] = sh25min;
	//	printf("sh[k] %d %g\n",k,sh[k]);

        shb = std::max(sh[k], 0.02f);
        sm[k] = std::min(sm[k], Prlimit * shb);

        /**** Level 3 : start ****/
        if (closure >= 3.0f) {
            t2sq = qdiv * b2 * elsq * sh[k] * dtl[k] * dtl[k];
            r2sq = qdiv * b2 * elsq * sh[k] * dqw[k] * dqw[k];
            c2sq = qdiv * b2 * elsq * sh[k] * dtl[k] * dqw[k];
            t3sq = std::max(tsq[k] * abk + tsq[k - 1] * afk, 0.0f);
            r3sq = std::max(qsq[k] * abk + qsq[k - 1] * afk, 0.0f);
            c3sq = cov[k] * abk + cov[k - 1] * afk;

            c3sq = std::copysign(std::min(std::abs(c3sq), sqrt(t3sq * r3sq)), c3sq);

            vtt = 1.0f + vt[k] * abk + vt[k - 1] * afk;
            vqq = tv0 + vq[k] * abk + vq[k - 1] * afk;

            t2sq = vtt * t2sq + vqq * c2sq;
            r2sq = vtt * c2sq + vqq * r2sq;
            c2sq = std::max(vtt * t2sq + vqq * r2sq, 0.0);
            t3sq = vtt * t3sq + vqq * c3sq;
            r3sq = vtt * c3sq + vqq * r3sq;
            c3sq = std::max(vtt * t3sq + vqq * r3sq, 0.0);

            cw25 = e1 * (e2 + 3.0f * c1 * e5c * gmel * qdiv * qdiv) / (3.0f * eden);

            // ** Limitation on q, instead of L/q **
            dlsq = elsq;
            if (q3sq / dlsq < -gh[k]) q3sq = -dlsq * gh[k];

            // ** Limitation on c3sq (0.12 =< cw =< 0.76) **
            auh = 27.0f * a1 * ((a2 * a2fac) * (a2 * a2fac)) * b2 * (gtr) * (gtr);
            aum = 54.0f * (a1 * a1) * (a2 * a2fac) * b2 * c1 * (gtr);
            adh = 9.0f * a1 * ((a2 * a2fac) * (a2 * a2fac)) * (12.0f * a1 + 3.0f * b2) * (gtr) * (gtr);
            adm = 18.0f * (a1 * a1) * (a2 * a2fac) * (b2 - 3.0f * (a2 * a2fac)) * (gtr);

            aeh = (9.0f * a1 * ((a2 * a2fac) * (a2 * a2fac)) * b1 + 9.0f * a1 * ((a2 * a2fac) * (a2 * a2fac)) * (12.0f * a1 + 3.0f * b2)) * (gtr);
            aem = 3.0f * a1 * (a2 * a2fac) * b1 * (3.0f * (a2 * a2fac) + 3.0f * b2 * c1 + (18.0f * a1 * c1 - b2)) + (18.0f) * (a1 * a1) * (a2 * a2fac) * (b2 - 3.0f * (a2 * a2fac));

            Req = -aeh / aem;
            Rsl = (auh + aum * Req) / (3.0f * adh + 3.0f * adm * Req);
            // For now, use default values, since tests showed little/no sensitivity
            Rsl = 0.12;
            Rsl2 = 1.0f - 2.0f * Rsl;

            // JOE-Canuto/Kitamura mod
            e2 = q3sq - e2c * ghel * a2fac * qdiv * qdiv;
            e3 = q3sq + e3c * ghel * a2fac * a2fac * qdiv * qdiv;
            e4 = q3sq - e4c * ghel * a2fac * qdiv * qdiv;
            eden = e2 * e4 + e3 * e5c * gmel * qdiv * qdiv;

            wden = cc3 * gtr * gtr * dlsq * dlsq / elsq * qdiv * qdiv * (e2 * e4c * a2fac - e3c * e5c * gmel * a2fac * a2fac * qdiv * qdiv);

            if (wden != 0.0f) {
                clow = q3sq * (0.12f - cw25) * eden / wden;
                cupp = q3sq * (0.76f - cw25) * eden / wden;
                if (wden > 0.0f) {
                    c3sq = std::min(std::max(c3sq, c2sq + clow), c2sq + cupp);
                } else {
                    c3sq = std::max(std::min(c3sq, c2sq + clow), c2sq + cupp);
                }
            }

            e1 = e2 + e5c * gmel * qdiv * qdiv;
            eden = std::max(eden, 1.0e-20);

            // JOE-Canuto/Kitamura mod
            e6c = 3.0f * (a2 * a2fac) * cc3 * gtr * dlsq / elsq;

            // ** for Gamma_theta **
            if (t2sq >= 0.0f) {
                enumc = std::max(qdiv * e6c * (t3sq - t2sq), 0.0);
            } else {
                enumc = std::min(qdiv * e6c * (t3sq - t2sq), 0.0);
            }
            gamt = -e1 * enumc / eden;

            // ** for Gamma_q **
            if (r2sq >= 0.0f) {
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
                if (sh[k] < -0.3f || sm[k] < -0.3f || qke[k] < -0.1f || std::abs(smd) > 2.0f) {
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

        // Add min background stability funciton (diffusivity) within model levels
        // with active plumes and clouds
        cldavg = 0.5f * (cldfra_bl1D[k - 1] + cldfra_bl1D[k]);
        if (edmf_a1[k] > 0.001 || cldavg > 0.02) {
            sm[k] = std::max(sm[k], 0.03f * std::min(10.0f * edmf_a1[k] * edmf_w1[k], 1.0f));
            sh[k] = std::max(sh[k], 0.03f * std::min(10.0f * edmf_a1[k] * edmf_w1[k], 1.0f));
            sm[k] = std::max(sm[k], 0.05f * std::min(cldavg, 1.0f));
            sh[k] = std::max(sh[k], 0.05f * std::min(cldavg, 1.0f));
	    //	    printf("sh[k] %d %g\n",k,sh[k]);
        }

        elq = el[k] * qkw[k];
        elh = elq * qdiv;

        // Production of TKE (pdk), T-variance (pdt),
        // q-variance (pdq), and covariance (pdc)
        pdk[k] = elq * (sm[k] * gm[k] + sh[k] * gh[k] + gamv) + 0.5f * TKEprodTD[k];
        pdt[k] = elh * (sh[k] * dtl[k] + gamt) * dtl[k];
        pdq[k] = elh * (sh[k] * dqw[k] + gamq) * dqw[k];
        pdc[k] = elh * (sh[k] * dtl[k] + gamt) * dqw[k] * 0.5f + elh * (sh[k] * dqw[k] + gamq) * dtl[k] * 0.5;

        // Countergradient terms
        tcd[k] = elq * gamt;
        qcd[k] = elq * gamq;

        // Eddy Diffusivity/Viscosity divided by dz
        dfm[k] = elq * sm[k] / dzk;
        dfh[k] = elq * sh[k] / dzk;
        dfq[k] = dfm[k];

        if (tke_budget == 1) {
            qSHEAR1D[k] = elq * sm[k] * gm[k];
            qBUOY1D[k] = elq * (sh[k] * gh[k] + gamv) + 0.5f * TKEprodTD[k];
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
            dfm[k] = dfm[k] + dfm[k] * rstoch_col[k] * 1.5f * std::max(float(exp(-std::max(zw[k] - 8000.0f, 0.0f) / 2000.0f)), 0.001f);
            dfh[k] = dfh[k] + dfh[k] * rstoch_col[k] * 1.5f * std::max(float(exp(-std::max(zw[k] - 8000.0f, 0.0f) / 2000.0f)), 0.001f);
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
//                         # zw(1)=0.0f and zw(k)=zw(k-1)+dz(k-1)
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
void mym_initialize_cc(const int &kts,const int &kte,const float &xland, float *dz, float &dx, float *zw, float *u, float *v, float *thl, float *qw,const float &zi, float *theta, float *thetav, float *sh, float *sm, const float& ust, const float &rmo, float* el, float *qke, float* tsq, float* qsq, float* cov, const float& Psig_bl, float *cldfra_bl1D, int &bl_mynn_mixlength, float *edmf_w1, float *edmf_a1, int &INITIALIZE_QKE, int &spp_pbl, float *rstoch_col,const float & karman,const float& tv0,const float& gtr) {
    float phm, vkz, elq, elv, b1l, b2l, pmz = 1.0, phh = 1.0, flt = 0.0, fltv = 0.0, flq = 0.0, tmpq;
    int k, l, lmax;
    float ql[kte-kts]; 
    float vt[kte-kts];
    float vq[kte-kts];
    float pdk[kte-kts], pdt[kte-kts], pdq[kte-kts],pdc[kte-kts],dtl[kte-kts],dqw[kte-kts],dtv[kte-kts],gm[kte-kts],gh[kte-kts],qkw[kte-kts];

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
        qke[kts] = 1.5f * ust*ust * std::cbrt((b1*pmz) * (b1*pmz));
        for (k = kts + 1; k <= kte; k++) {
            qke[k] = qke[kts] * std::max((ust * 700.0f - zw[k]) / (std::max(ust, 0.01f) * 700.0f), 0.01f);
        }
    }

    phm = phh * b2 / cbrt(b1 * pmz);
    tsq[kts] = phm * ((flt / ust) * (flt / ust));
    qsq[kts] = phm * ((flq / ust) * (flq / ust));
    cov[kts] = phm * (flt / ust) * (flq / ust);
    for (k = kts + 1; k <= kte; k++) {
        vkz = karman * zw[k];
        el[k] = vkz / (1.0f + vkz / 100.0f);
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
        vkz = karman * 0.5f * dz[kts];
        elv = 0.5f * (el[kts + 1] + el[kts]) / vkz;
        if (INITIALIZE_QKE==1) {
            // WRF has (b1*pmz*elv)**(2.0/3.0), not the `two_third` constant
            qke[kts] = 1.0f * std::max(ust, 0.02f) * std::max(ust, 0.02f) * std::cbrt((b1 * pmz * elv) * (b1 * pmz * elv));
        }

        phm = phh * b2 / std::cbrt(b1 * pmz / (elv*elv));
        tsq[kts] = phm * ((flt / ust) * (flt / ust));
        qsq[kts] = phm * ((flq / ust) * (flq / ust));
        cov[kts] = phm * (flt / ust) * (flq / ust);
        
        for (k = kts + 1; k <= kte - 1; k++) {
            b1l = b1 * 0.25f * (el[k + 1] + el[k]);
            tmpq = std::min(std::max(b1l * (pdk[k + 1] + pdk[k]), qkemin), 125.0f);
            if (INITIALIZE_QKE==1) {
                // WRF has tmpq**two_thirds, where two_thirds is a constant (kind_phys)
                qke[k] = std::pow(tmpq, 2.0f/3.0f);
            }

            if (qke[k] <= 0.0f) {
                b2l = 0.0;
            } else {
                b2l = b2 * (b1l / b1) / sqrt(qke[k]);
            }
            tsq[k] = b2l * (pdt[k + 1] + pdt[k]);
            qsq[k] = b2l * (pdq[k + 1] + pdq[k]);
            cov[k] = b2l * (pdc[k + 1] + pdc[k]);
        }
    }
    
    if (INITIALIZE_QKE==1) {
        qke[kts] = 0.5f * (qke[kts] + qke[kts + 1]);
        qke[kte] = qke[kte - 1];
    }

    tsq[kte] = tsq[kte - 1];
    qsq[kte] = qsq[kte - 1];
    cov[kte] = cov[kte - 1];
}




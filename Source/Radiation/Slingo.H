//------------------------------------------------------------------------------------------------
//  Implements Slingo Optics for MG/RRTMG for liquid clouds and
//  a copy of the old cloud routine for reference
//------------------------------------------------------------------------------------------------
#ifndef ERF_SLINGO_H_
#define ERF_SLINGO_H_

class Slingo {
  public:
  static void slingo_liq_optics_sw(int ncol, int nlev, int nswbands, real2d& cldn, real2d& cliqwp, real2d& rel,
                                   real3d& liq_tau, real3d& liq_tau_w, real3d& liq_tau_w_g, real3d& liq_tau_w_f) {
    real1d wavmin("wavmin", nswbands);
    real1d wavmax("wavmax", nswbands);

    // Minimum cloud amount (as a fraction of the grid-box area) to
    // distinguish from clear sky
    real cldmin = 1.0e-80;

    // Decimal precision of cloud amount (0 -> preserve full resolution;
    // 10^-n -> preserve n digits of cloud amount)
    real cldeps = 0.0;

    // A. Slingo's data for cloud particle radiative properties (from 'A GCM
    // Parameterization for the Shortwave Properties of Water Clouds' JAS
    // vol. 46 may 1989 pp 1419-1427)
    std::vector<real> abarl =  // A coefficient for extinction optical depth
         { 2.817e-02, 2.682e-02, 2.264e-02, 1.281e-02 };
    std::vector<real>  bbarl =   // B coefficient for extinction optical depth
         { 1.305, 1.346, 1.454, 1.641 };
    std::vector<real> cbarl =   // C coefficient for single scat albedo
         { -5.62e-08 ,-6.94e-06 ,4.64e-04 ,0.201 };
    std::vector<real> dbarl =   // D coefficient for single  scat albedo
         { 1.63e-07, 2.35e-05 ,1.24e-03 ,7.56e-03 };
    std::vector<real> ebarl =   // E coefficient for asymmetry parameter
         { 0.829, 0.794, 0.754, 0.826 };
    std::vector<real> fbarl =   // F coefficient for asymmetry parameter
         { 2.482e-03, 4.226e-03, 6.560e-03, 4.353e-03 };

    real abarli;        // A coefficient for current spectral band
    real bbarli;        // B coefficient for current spectral band
    real cbarli;        // C coefficient for current spectral band
    real dbarli;        // D coefficient for current spectral band
    real ebarli;        // E coefficient for current spectral band
    real fbarli;        // F coefficient for current spectral band

    // Caution... A. Slingo recommends no less than 4.0 micro-meters nor
    // greater than 20 micro-meters. Here we set effective radius limits
    // for liquid to the range 4.2 < rel < 16 micron (Slingo 89)
    const real rel_min = 4.2;
    const real rel_max = 16.;

    int indxsl;
    real tmp1l, tmp2l, tmp3l, g;

    RadConstants::get_sw_spectral_boundaries(wavmin,wavmax,RadConstants::micrometer);

    for (auto ns=0; ns<nswbands; ++ns) {
       // Set index for cloud particle properties based on the wavelength,
       // according to A. Slingo (1989) equations 1-3:
       // Use index 1 (0.25 to 0.69 micrometers) for visible
       // Use index 2 (0.69 - 1.19 micrometers) for near-infrared
       // Use index 3 (1.19 to 2.38 micrometers) for near-infrared
       // Use index 4 (2.38 to 4.00 micrometers) for near-infrared
       if(wavmax(ns) <= 0.7) {
          indxsl = 1;
       } else if(wavmax(ns) <= 1.25) {
          indxsl = 2;
       } else if(wavmax(ns) <= 2.38) {
          indxsl = 3;
       } else if(wavmax(ns) > 2.38) {
          indxsl = 4;
       }

       // Set cloud extinction optical depth, single scatter albedo,
       // asymmetry parameter, and forward scattered fraction:
       abarli = abarl[indxsl];
       bbarli = bbarl[indxsl];
       cbarli = cbarl[indxsl];
       dbarli = dbarl[indxsl];
       ebarli = ebarl[indxsl];
       fbarli = fbarl[indxsl];

       for (auto k=0; k<nlev; ++k) {
          for (auto i=0; i<ncol; ++i) {

             // note that optical properties for liquid valid only
             // in range of 4.2 > rel > 16 micron (Slingo 89)
             if (cldn(i,k) >= cldmin && cldn(i,k) >= cldeps) {
                tmp1l = abarli + bbarli/std::min(std::max(rel_min,rel(i,k)),rel_max);
                liq_tau(ns,i,k) = 1000.*cliqwp(i,k)*tmp1l;
             } else {
                liq_tau(ns,i,k) = 0.0;
             }

             tmp2l = 1. - cbarli - dbarli*std::min(std::max(rel_min,rel(i,k)),rel_max);
             tmp3l = fbarli*std::min(std::max(rel_min,rel(i,k)),rel_max);
             // Do not let single scatter albedo be 1.  Delta-eddington solution
             // for non-conservative case has different analytic form from solution
             // for conservative case, and raddedmx is written for non-conservative case.
             liq_tau_w(ns,i,k) = liq_tau(ns,i,k) * std::min(tmp2l,.999999);
             g = ebarli + tmp3l;
             liq_tau_w_g(ns,i,k) = liq_tau_w(ns,i,k) * g;
             liq_tau_w_f(ns,i,k) = liq_tau_w(ns,i,k) * g * g;

           } // End do i=1,ncol
       } // End do k=1,nlev
     } // nswbands
  }

  static void slingo_liq_optics_lw(int ncol, int nlev, int nlwbands, real2d& cldn, real2d& iclwpth, real2d& iciwpth, real3d& abs_od) {
    real2d ficemr("ficemr", ncol, nlev);
    real2d cwp("cwp", ncol, nlev);
    real2d cldtau("cldtau", ncol, nlev);

    real kabs, kabsi;

    // longwave liquid absorption coeff (m**2/g)
    const real kabsl = 0.090361;

    for(auto k=0; k<nlev; ++k) {
       for(auto i=0; i<ncol; ++i) {
         cwp   (i,k) = 1000.0 * iclwpth(i,k) + 1000.0 * iciwpth(i, k);
         ficemr(i,k) = 1000.0 * iciwpth(i,k)/(std::max(1.e-18, cwp(i,k)));
       }
    }

    for (auto k=0; k<nlev; ++k) {
       for (auto i=0; i<ncol; ++i) {
          // Note from Andrew Conley:
          //  Optics for RK no longer supported, This is constructed to get
          //  close to bit for bit.  Otherwise we could simply use liquid water path
          //note that optical properties for ice valid only
          //in range of 13 > rei > 130 micron (Ebert and Curry 92)
          kabs = kabsl*(1.-ficemr(i,k));
          cldtau(i,k) = kabs*cwp(i,k);
       }
    }

    for (auto lwband=0; lwband<nlwbands; ++lwband) {
       for (auto icol=0; icol<ncol; ++icol) {
          for (auto ilev=0; ilev<nlev; ++ilev) {
            abs_od(lwband,icol,ilev)=cldtau(icol,ilev);
          }
       }
    }
  }
};
#endif

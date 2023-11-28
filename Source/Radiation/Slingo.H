//------------------------------------------------------------------------------------------------
//  Implements Slingo Optics for MG/RRTMG for liquid clouds and
//  a copy of the old cloud routine for reference
//------------------------------------------------------------------------------------------------
#ifndef ERF_SLINGO_H_
#define ERF_SLINGO_H_
using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

class Slingo {
  public:
  static void slingo_liq_optics_sw(int ncol, int nlev, int nswbands, 
                                   const real2d& cldn, const real2d& cliqwp, const real2d& rel,
                                   const real3d& liq_tau, const real3d& liq_tau_w, 
                                   const real3d& liq_tau_w_g, const real3d& liq_tau_w_f) {
    real1d wavmin("wavmin", nswbands);
    real1d wavmax("wavmax", nswbands);

    // Minimum cloud amount (as a fraction of the grid-box area) to
    // distinguish from clear sky
    const real cldmin = 1.0e-80;

    // Decimal precision of cloud amount (0 -> preserve full resolution;
    // 10^-n -> preserve n digits of cloud amount)
    const real cldeps = 0.0;

    // A. Slingo's data for cloud particle radiative properties (from 'A GCM
    // Parameterization for the Shortwave Properties of Water Clouds' JAS
    // vol. 46 may 1989 pp 1419-1427)
    real1d abarl("abarl", 4);  // A coefficient for extinction optical depth
    real1d bbarl("bbarl", 4);  // B coefficient for extinction optical depth
    real1d cbarl("cbarl", 4);  // C coefficient for extinction optical depth
    real1d dbarl("dbarl", 4);  // D coefficient for extinction optical depth
    real1d ebarl("ebarl", 4);  // E coefficient for extinction optical depth
    real1d fbarl("fbarl", 4);  // F coefficient for extinction optical depth
    parallel_for(SimpleBounds<1>(1), YAKL_LAMBDA (int i) {
       abarl(1) = 2.817e-02;
       abarl(2) = 2.682e-02;
       abarl(3) = 2.264e-02;
       abarl(4) = 1.281e-02;
       bbarl(1) = 1.305;
       bbarl(2) = 1.346;
       bbarl(3) = 1.454;
       bbarl(4) = 1.641;
       cbarl(1) = -5.62e-08;
       cbarl(2) = -6.94e-06;
       cbarl(3) = 4.64e-04;
       abarl(4) = 0.201;
       dbarl(1) = 1.63e-07;
       dbarl(2) = 2.35e-05;
       dbarl(3) = 1.24e-03;
       dbarl(4) = 7.56e-03;
       ebarl(1) = 0.829;
       ebarl(2) = 0.794;
       ebarl(3) = 0.754;
       ebarl(4) = 0.826;
       fbarl(1) = 2.482e-03;
       fbarl(2) = 4.226e-03;
       fbarl(3) = 6.560e-03;
       fbarl(4) = 4.353e-03;
     });

    // Caution... A. Slingo recommends no less than 4.0 micro-meters nor
    // greater than 20 micro-meters. Here we set effective radius limits
    // for liquid to the range 4.2 < rel < 16 micron (Slingo 89)
    const real rel_min = 4.2;
    const real rel_max = 16.;

    int indxsl;

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
       parallel_for(SimpleBounds<2>(nlev, ncol), YAKL_LAMBDA (int k, int i) {
          auto abarli = abarl(indxsl);
          auto bbarli = bbarl(indxsl);
          auto cbarli = cbarl(indxsl);
          auto dbarli = dbarl(indxsl);
          auto ebarli = ebarl(indxsl);
          auto fbarli = fbarl(indxsl);
          real tmp1l, tmp2l, tmp3l, g;

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
       });
     } // nswbands
  }

  static void slingo_liq_optics_lw(int ncol, int nlev, int nlwbands, const real2d& cldn,
                                   const real2d& iclwpth, const real2d& iciwpth, const real3d& abs_od) {
    real2d ficemr("ficemr", ncol, nlev);
    real2d cwp("cwp", ncol, nlev);
    real2d cldtau("cldtau", ncol, nlev);

    parallel_for(SimpleBounds<2>(nlev, ncol), YAKL_LAMBDA (int k, int i) {
       cwp   (i,k) = 1000.0 * iclwpth(i,k) + 1000.0 * iciwpth(i, k);
       ficemr(i,k) = 1000.0 * iciwpth(i,k)/(std::max(1.e-18, cwp(i,k)));
    });

    parallel_for(SimpleBounds<2>(nlev, ncol), YAKL_LAMBDA (int k, int i) {
       // Note from Andrew Conley:
       //  Optics for RK no longer supported, This is constructed to get
       //  close to bit for bit.  Otherwise we could simply use liquid water path
       //note that optical properties for ice valid only
       //in range of 13 > rei > 130 micron (Ebert and Curry 92)
       real kabsl = 0.090361;
       auto kabs = kabsl*(1.-ficemr(i,k));
       cldtau(i,k) = kabs*cwp(i,k);
    });

    parallel_for(SimpleBounds<3>(nlwbands, ncol, nlev), YAKL_LAMBDA (int lwband, int icol, int ilev) {
       abs_od(lwband,icol,ilev) = cldtau(icol,ilev);
    });
  }
};
#endif

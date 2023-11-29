
#ifndef ERF_EBERT_CURRY_H_
#define ERF_EBERT_CURRY_H_
using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

class EbertCurry {
  public:
  static constexpr real scalefactor = 1.;  //500._r8/917._r8

  static void ec_ice_optics_sw(int ncol, int nlev, int nswbands, const real2d& cldn, const real2d& cicewp, const real2d& rei,
                               const real3d& ice_tau, const real3d& ice_tau_w, const real3d& ice_tau_w_g, const real3d& ice_tau_w_f) {

    real1d wavmin("wavmin", nswbands);
    real1d wavmax("wavmax", nswbands);

    // ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
    real1d abari("abari", 4);      // a coefficient for extinction optical depth
    real1d bbari("bbari", 4);     // b coefficient for extinction optical depth
    real1d cbari("cbari", 4);      // c coefficient for single scat albedo
    real1d dbari("dbari", 4);      // d coefficient for single scat albedo
    real1d ebari("ebari", 4);      // e coefficient for asymmetry parameter
    real1d fbari("fbari", 4);     // f coefficient for asymmetry parameter

    parallel_for(SimpleBounds<1>(1), YAKL_LAMBDA (int i) {
      abari(1) = 3.448e-03;
      abari(2) = 3.448e-03;
      abari(3) = 3.448e-03;
      abari(4) = 3.448e-03;
      bbari(1) = 2.431;
      bbari(2) = 2.431;
      bbari(3) = 2.431;
      bbari(4) = 2.431;
      cbari(1) = 1.00e-05;
      cbari(2) = 1.10e-04;
      cbari(3) = 1.861e-02;
      cbari(4) = 0.46658;
      dbari(1) = 0.0;
      dbari(2) = 1.405e-05;
      dbari(3) = 8.328e-04;
      dbari(4) = 2.05e-05;
      ebari(1) = 0.7661;
      ebari(2) = 0.7730;
      ebari(3) = 0.794;
      ebari(4) = 0.9595;
      fbari(1) = 5.851e-04;
      fbari(2) = 5.665e-04;
      fbari(3) = 7.267e-04;
      fbari(4) = 1.076e-04;
    });

    // Minimum cloud amount (as a fraction of the grid-box area) to
    // distinguish from clear sky
    constexpr real cldmin = 1.0e-80;

    // Decimal precision of cloud amount (0 -> preserve full resolution;
    // 10^-n -> preserve n digits of cloud amount)
    constexpr real cldeps = 0.0;

    // Optical properties for ice are valid only in the range of
    // 13 < rei < 130 micron (Ebert and Curry 92)
    constexpr real rei_min = 13.;
    constexpr real rei_max = 130.;

    //  integer :: ns, i, k, indxsl
    int indxsl;

    RadConstants::get_sw_spectral_boundaries(wavmin,wavmax,RadConstants::micrometer);

    for (auto ns=0; ns<nswbands; ++ns) {

       if(wavmax(ns) <= 0.7) {
          indxsl = 1;
       } else if(wavmax(ns) <= 1.25) {
          indxsl = 2;
       } else if(wavmax(ns) <= 2.38) {
          indxsl = 3;
       } else if(wavmax(ns) > 2.38) {
          indxsl = 4;
       }

       parallel_for(SimpleBounds<2>(nlev, ncol), YAKL_LAMBDA (int k, int i) {
          real tmp1i, tmp2i, tmp3i, g;
          auto abarii = abari(indxsl);
          auto bbarii = bbari(indxsl);
          auto cbarii = cbari(indxsl);
          auto dbarii = dbari(indxsl);
          auto ebarii = ebari(indxsl);
          auto fbarii = fbari(indxsl);

          // note that optical properties for ice valid only
          // in range of 13 > rei > 130 micron (Ebert and Curry 92)
          if (cldn(i,k) >= cldmin && cldn(i,k) >= cldeps) {
             tmp1i = abarii + bbarii/std::max(rei_min,std::min(scalefactor*rei(i,k),rei_max));
             ice_tau(ns,i,k) = 1000. * cicewp(i,k) * tmp1i;
          } else {
               ice_tau(ns,i,k) = 0.0;
          }

          tmp2i = 1. - cbarii - dbarii*std::min(std::max(rei_min,scalefactor*rei(i,k)),rei_max);
          tmp3i = fbarii*std::min(std::max(rei_min,scalefactor*rei(i,k)),rei_max);
          // Do not let single scatter albedo be 1.  Delta-eddington solution
          // for non-conservative case has different analytic form from solution
          // for conservative case, and raddedmx is written for non-conservative case.
          ice_tau_w(ns,i,k) = ice_tau(ns,i,k) * std::min(tmp2i,.999999);
          g = ebarii + tmp3i;
          ice_tau_w_g(ns,i,k) = ice_tau_w(ns,i,k) * g;
          ice_tau_w_f(ns,i,k) = ice_tau_w(ns,i,k) * g * g;
        });
     } // nswbands
  }


  static void ec_ice_optics_lw(int ncol, int nlev, int nlwbands, const real2d& cldn, const real2d& iclwpth,
              const real2d& iciwpth, const real2d& rei, const real3d& abs_od) {
     real2d ficemr("ficemr",ncol,nlev);
     real2d cwp("cwp",ncol,nlev);
     real2d cldtau("cldtau",ncol,nlev);

     // longwave liquid absorption coeff (m**2/g)
     const real kabsl = 0.090361;

     // Optical properties for ice are valid only in the range of
     // 13 < rei < 130 micron (Ebert and Curry 92)
     const real rei_min = 13.;
     const real rei_max = 130.;

     parallel_for(SimpleBounds<2>(nlev, ncol), YAKL_LAMBDA (int k, int i) {
        cwp(i,k) = 1000.0 *iciwpth(i,k) + 1000.0 *iclwpth(i,k);
        ficemr(i,k) = 1000.0*iciwpth(i,k)/(std::max(1.e-18,cwp(i,k)));
     });

     parallel_for(SimpleBounds<2>(nlev, ncol), YAKL_LAMBDA (int k, int i) {
        // Note from Andrew Conley:
        //  Optics for RK no longer supported, This is constructed to get
        //  close to bit for bit.  Otherwise we could simply use ice water path
        //  note that optical properties for ice valid only
        //  in range of 13 > rei > 130 micron (Ebert and Curry 92)
        auto kabsi = 0.005 + 1./std::min(std::max(rei_min,scalefactor*rei(i,k)),rei_max);
        auto kabs =  kabsi*ficemr(i,k); // kabsl*(1.-ficemr(i,k)) + kabsi*ficemr(i,k);
        cldtau(i,k) = kabs*cwp(i,k);
     });

     parallel_for(SimpleBounds<3>(nlwbands, ncol, nlev), YAKL_LAMBDA (int lwband, int icol, int ilev) {
             abs_od(lwband,icol,ilev)=cldtau(icol,ilev);
     });
  }
};
#endif


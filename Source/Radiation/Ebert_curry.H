
#ifndef ERF_EBERT_CURRY_H_
#define ERF_EBERT_CURRY_H_
using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

class EbertCurry {
  public:
  static constexpr real scalefactor = 1.;  //500._r8/917._r8

  static void ec_ice_optics_sw(int ncol, int nlev, int nswbands, real2d& cldn, real2d& cicewp, real2d& rei,
                               real3d& ice_tau, real3d& ice_tau_w, real3d& ice_tau_w_g, real3d& ice_tau_w_f) {

    real1d wavmin("wavmin", nswbands);
    real1d wavmax("wavmax", nswbands);

    // ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
    std::vector<real> abari =      // a coefficient for extinction optical depth
         { 3.448e-03, 3.448e-03, 3.448e-03, 3.448e-03 };
    std::vector<real> bbari =      // b coefficient for extinction optical depth
         { 2.431, 2.431, 2.431, 2.431 };
    std::vector<real> cbari =      // c coefficient for single scat albedo
         { 1.00e-05, 1.10e-04, 1.861e-02, .46658 };
    std::vector<real> dbari =      // d coefficient for single scat albedo
         { 0.0, 1.405e-05, 8.328e-04, 2.05e-05 };
    std::vector<real> ebari =      // e coefficient for asymmetry parameter
         { 0.7661, 0.7730, 0.794, 0.9595 };
    std::vector<real> fbari =     // f coefficient for asymmetry parameter
         { 5.851e-04, 5.665e-04, 7.267e-04, 1.076e-04 };

    real abarii;           // A coefficient for current spectral band
    real bbarii;           // B coefficient for current spectral band
    real cbarii;           // C coefficient for current spectral band
    real dbarii;           // D coefficient for current spectral band
    real ebarii;           // E coefficient for current spectral band
    real fbarii;           // F coefficient for current spectral band

    // Minimum cloud amount (as a fraction of the grid-box area) to
    // distinguish from clear sky
    const real cldmin = 1.0e-80;

    // Decimal precision of cloud amount (0 -> preserve full resolution;
    // 10^-n -> preserve n digits of cloud amount)
    const real cldeps = 0.0;

    // Optical properties for ice are valid only in the range of
    // 13 < rei < 130 micron (Ebert and Curry 92)
    const real rei_min = 13.;
    const real rei_max = 130.;

    //  integer :: ns, i, k, indxsl
    int indxsl;
    real tmp1i, tmp2i, tmp3i, g;

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

       abarii = abari[indxsl];
       bbarii = bbari[indxsl];
       cbarii = cbari[indxsl];
       dbarii = dbari[indxsl];
       ebarii = ebari[indxsl];
       fbarii = fbari[indxsl];

       for (auto k=0; k<nlev; ++k) {
          for (auto i=0; i<ncol; ++i) {

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

          } //  End do i=1,ncol
        }    // End do k=1,nlev
     } // nswbands
  }


  static void ec_ice_optics_lw(int ncol, int nlev, int nlwbands, real2d& cldn, real2d& iclwpth, real2d& iciwpth, real2d& rei, real3d& abs_od) {
     real2d ficemr("ficemr",ncol,nlev);
     real2d cwp("cwp",ncol,nlev);
     real2d cldtau("cldtau",ncol,nlev);
     real kabs, kabsi;

     // longwave liquid absorption coeff (m**2/g)
     const real kabsl = 0.090361;

     // Optical properties for ice are valid only in the range of
     // 13 < rei < 130 micron (Ebert and Curry 92)
     const real rei_min = 13.;
     const real rei_max = 130.;

     for (auto k=0; k<nlev; ++k) {
        for (auto i=0; i<ncol; ++i) {
           cwp(i,k) = 1000.0 *iciwpth(i,k) + 1000.0 *iclwpth(i,k);
           ficemr(i,k) = 1000.0*iciwpth(i,k)/(std::max(1.e-18,cwp(i,k)));
        }
     }

     for (auto k=0; k<nlev; ++k) {
        for (auto i=0; i<ncol; ++i) {
           // Note from Andrew Conley:
           //  Optics for RK no longer supported, This is constructed to get
           //  close to bit for bit.  Otherwise we could simply use ice water path
           //  note that optical properties for ice valid only
           //  in range of 13 > rei > 130 micron (Ebert and Curry 92)
           kabsi = 0.005 + 1./std::min(std::max(rei_min,scalefactor*rei(i,k)),rei_max);
           kabs =  kabsi*ficemr(i,k); // kabsl*(1.-ficemr(i,k)) + kabsi*ficemr(i,k);
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


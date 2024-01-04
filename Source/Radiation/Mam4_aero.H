//
// Mam4 aerosol model class
//
// parameterizes aerosol coefficients using chebychev polynomial
// parameterize aerosol radiative properties in terms of
// surface mode wet radius and wet refractive index
// Ghan and Zaveri, JGR 2007.
//
#ifndef ERF_MAM4_AERO_H_
#define ERF_MAM4_AERO_H_

#include "YAKL_netcdf.h"
#include "rrtmgp_const.h"
#include "Mam4_constituents.H"
#include "Rad_constants.H"
#include "ERF_Constants.H"
#include "Modal_aero_wateruptake.H"

using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

class Mam4_aer {
   public:
    int ncol, nlev, nspec, top_lev;
    int nswbands, nlwbands;
    std::string modal_optics_file;   // full pathname for modal optics dataset
    std::string water_refindex_file; // full pathname for water refractive index dataset

    // Dimension sizes in coefficient arrays used to parameterize aerosol radiative properties
    // in terms of refractive index and wet radius
    int ncoef = 5;
    int prefr = 7;
    int prefi = 10;
    int nmodes = 4;
    int n_diag = 1;
    bool clim_modal_aero = true; // true when radiatively constituents present (nmodes>0)

    real xrmin, xrmax;
    // refractive index for water read in read_water_refindex
    real1d crefwswr; //(nswbands) complex refractive index for water visible
    real1d crefwswi;
    real1d crefwlwr; // complex refractive index for water infrared
    real1d crefwlwi;

    // These are defined as module level variables to aviod allocation-deallocation in a loop
    real3d dgnumdry_m;  // number mode dry diameter for all modes
    real3d dgnumwet_m;  // number mode wet diameter for all modes
    real3d qaerwat_m;   // aerosol water (g/g) for all modes
    real3d wetdens_m;

    // aerosol constitutents
    MamConstituents mam_consti;
  public:
   // initialization
   inline
   void initialize(int ncol_, int nlev_, int top_lev_, int nswbands_, int nlwbands_) {
      ncol = ncol_;
      nlev = nlev_;
      top_lev = top_lev_;
      nswbands = nswbands_;
      nlwbands = nlwbands_;
      clim_modal_aero = true;

      const real rmmin = 0.01e-6;
      const real rmmax = 25.e-6;
      xrmin = std::log(rmmin);
      xrmax = std::log(rmmax);

      // Check that dimension sizes in the coefficient arrays used to
      // parameterize aerosol radiative properties are consistent between this
      // module and the mode physprop files.
      for(auto ilist = 0; ilist < n_diag; ++ilist) {
        mam_consti.get_nmodes(ilist, nmodes);
        for(auto m = 0; m < nmodes; ++m) {
           mam_consti.get_mode_props(ilist, m, ncoef, prefr, prefi);
        }
      }

      const char* erf_home = std::getenv("ERF_HOME");
      auto erf_rad_data_dir = (std::string)erf_home + "/Source/Radiation/data/";
      water_refindex_file = erf_rad_data_dir+"water_refindex_rrtmg_c080910.nc";
      read_water_refindex(water_refindex_file);

      // Allocate dry and wet size variables in an OMP PARRALLEL region as these
      // arrays are private for each thread and needs to be allocated for each OMP thread
      dgnumdry_m = real3d("dgnumdry", ncol, nlev, nmodes);
      dgnumwet_m = real3d("dgnumwet", ncol, nlev, nmodes);
      qaerwat_m  = real3d("qaerwat" , ncol, nlev, nmodes);
      wetdens_m  = real3d("wetdens" , ncol, nlev, nmodes);

      crefwswr = real1d("crefwswr", nswbands);
      crefwswi = real1d("crefwswi", nswbands);
      crefwlwr = real1d("crefwlwr", nlwbands);
      crefwlwi = real1d("crefwlwi", nlwbands);
   }

   //
   // read water refractive index file and set module data
   //
   inline
   void read_water_refindex(std::string infilename) {
     yakl::SimpleNetCDF water_ref_file;
     water_ref_file.open(infilename, yakl::NETCDF_MODE_READ);

     // read the dimensions
     int nlw_bands = water_ref_file.getDimSize( "lw_band" );
     int nsw_bands = water_ref_file.getDimSize( "sw_band" );

     if (nswbands != nsw_bands || nlwbands != nlw_bands)
        amrex::Print() << "ERROR - file and bandwith values do not match: "
                       << "\n nswbands: " << nswbands << "; nsw_bands: " << nsw_bands
                       << "\n nlwbands: " << nlwbands << "; nlw_bands: " << nlw_bands << "\n";

     // Local variables
     real1d refrwsw("refrwsw", nsw_bands), refiwsw("refiwsw", nsw_bands); // real, imaginary ref index for water visible
     real1d refrwlw("refrwlw", nlw_bands), refiwlw("refiwlw", nlw_bands); // real, imaginary ref index for water infrared

     // read variables
     water_ref_file.read( refrwsw, "refindex_real_water_sw");
     water_ref_file.read( refiwsw, "refindex_im_water_sw");
     water_ref_file.read( refrwlw, "refindex_real_water_lw");
     water_ref_file.read( refiwlw, "refindex_im_water_lw");
   }

   //
   // modal size  parameters
   //
   inline
   void modal_size_parameters(real sigma_logr_aer,
                              const real2d& dgnumwet, const real2d& radsurf,
                              const real2d& logradsurf, const real3d& cheb) {
     real1d xrad("xrad", ncol);

     auto alnsg_amode = std::log(sigma_logr_aer);
     auto explnsigma = std::exp(2.0*alnsg_amode*alnsg_amode);
     top_lev = 1;

     for(auto k = top_lev; k <= nlev; ++k) {
       for(auto i = 1; i <= ncol; ++i) {
         // convert from number mode diameter to surface area
         radsurf(i,k) = 0.5*dgnumwet(i,k)*explnsigma;
         logradsurf(i,k) = std::log(radsurf(i,k));

         // normalize size parameter
         xrad(i) = std::max(logradsurf(i,k),xrmin);
         xrad(i) = std::min(xrad(i),xrmax);
         xrad(i) = (2.*xrad(i)-xrmax-xrmin)/(xrmax-xrmin);

         // chebyshev polynomials
         cheb(1,i,k) = 1.;
         cheb(2,i,k) = xrad(i);
         for(auto nc = 3; nc < ncoef; ++nc) {
           cheb(nc,i,k) = 2.*xrad(i)*cheb(nc-1,i,k)-cheb(nc-2,i,k);
         }
       }
     }
   }

   //
   // bilinear interpolation of table
   //
   inline
   void binterp(const real3d& table, int ncol, int km, int im, int jm,
                const real1d& x, const real1d& y, const real1d& xtab, const real1d& ytab,
                const int1d& ix, const int1d& jy, const real1d& t, const real1d& u, const real2d& out) {
      // local variables
      real1d tu("tu", ncol), tuc("tuc", ncol), tcu("tcu", ncol), tcuc("tcuc", ncol);

      if(ix(1) <= 0) {
        if(im > 1) {
          parallel_for (SimpleBounds<1>(ncol), YAKL_LAMBDA (int ic) {
            int i;
            for (i = 1; i <= im; ++i) {
              if(x(ic) < xtab(i)) break;
            }

            ix(ic)   = std::max(i-1,1);
            auto ip1 = std::min(ix(ic)+1,im);
            auto dx  = (xtab(ip1)-xtab(ix(ic)));

            if(abs(dx) > 1.0-20) {
               t(ic)=(x(ic)-xtab(ix(ic)))/dx;
            } else {
              t(ic)=0.;
            }
          });
        } else {
          parallel_for (SimpleBounds<1>(ncol), YAKL_LAMBDA (int i) {
            ix(i) = 1;
            t(i)  = 0.;
          });
        }

        if(jm > 1) {
          parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA (int ic) {
             int j;
             for (j = 1; j <= jm; ++j) {
               if(y(ic) < ytab(j)) break;
             }

             jy(ic)  = std::max(j-1,1);
             auto jp1 = std::min(jy(ic)+1,jm);
             auto dy  = (ytab(jp1)-ytab(jy(ic)));

             if(std::abs(dy) > 1.e-20) {
                u(ic) = (y(ic)-ytab(jy(ic)))/dy;
                if(u(ic) < 0. || u(ic) > 1.)
                  amrex::Print() << "u= " << u(ic) << "; y= " << y(ic) << "; ytab= "
                                 << ytab(jy(ic)) << "; dy= " << dy << std::endl;
             } else {
               u(ic)=0.;
             }
          });
        } else {
          parallel_for (SimpleBounds<1>(ncol), YAKL_LAMBDA (int i) {
            jy(i) = 1;
            u(i)  = 0.;
          });
        }
     }

     parallel_for (SimpleBounds<2>(ncol, km), YAKL_LAMBDA (int ic, int k) {
        tu(ic)   = t(ic)*u(ic);
        tuc(ic)  = t(ic)-tu(ic);
        tcuc(ic) = 1.-tuc(ic)-u(ic);
        tcu(ic)  = u(ic)-tu(ic);
        auto jp1 = std::min(jy(ic)+1,jm);
        auto ip1 = std::min(ix(ic)+1,im);
        out(ic,k) = tcuc(ic)*table(k,ix(ic),jy(ic))+tuc(ic)*table(k,ip1,jy(ic))
                   +tu(ic)*table(k,ip1,jp1)+tcu(ic)*table(k,ix(ic),jp1);
     });
   } // subroutine binterp

   inline
   void modal_aero_calcsize_diag(int list_idx, const real3d& dgnum_m)
   {
     real2d dgncur_a("dgncur",ncol,nlev);   // (pcols,pver)
     real2d mode_num("mode_num",ncol,nlev); // mode number mixing ratio
     real2d specmmr("specmmr",ncol,nlev);   // specie mmr
     real2d dryvol_a("dryvol_a",ncol,nlev);   // interstital aerosol dry volume (cm^3/mol_air)

     real dgnum, dgnumhi, dgnumlo;
     real dgnyy, dgnxx;           // dgnumlo/hi of current mode
     real drv_a;                  // dry volume (cm3/mol_air)
     real dumfac, dummwdens;      // work variables
     real num_a0;                 // initial number (#/mol_air)
     real num_a;                  // final number (#/mol_air)
     real voltonumbhi, voltonumblo;
     real v2nmax, v2nmin;           // voltonumblo/hi of current mode
     real sigmag, alnsg;

     for (auto n = 1; n <= nmodes; ++n) {
        parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int i, int k) {
          dgncur_a(i,k) = dgnum_m(i,k,n);
        });

        // get mode properties
        mam_consti.get_mode_props(list_idx, n-1, dgnum, dgnumhi, dgnumlo, sigmag=sigmag);

        // get mode number mixing ratio
        mam_consti.rad_cnst_get_mode_num(list_idx, n-1, "a", mode_num);

        yakl::memset(dgncur_a, dgnum);
        yakl::memset(dryvol_a, 0.0);

        // compute dry volume mixrats =
        //      sum_over_components{ component_mass mixrat / density }
        mam_consti.get_mode_nspec(list_idx, n-1, nspec);
        for (auto l1 = 0; l1 < nspec; ++l1) {
           real specdens = 1.0e-20;
           mam_consti.rad_cnst_get_mam_mmr_by_idx(list_idx, n-1, l1, "a", specmmr);
           mam_consti.get_mam_props(list_idx, n-1, l1, specdens);

           // need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air
           dummwdens = 1.0 / specdens;
           top_lev = 1;

           for (auto k=top_lev; k <= nlev; ++k) {
             for (auto i=1; i<=ncol; ++i) {
                dryvol_a(i,k) = dryvol_a(i,k) + std::max(0.0, specmmr(i,k))*dummwdens;
             }
           }
        }

        alnsg       = std::log(sigmag);
        dumfac      = std::exp(4.5*std::pow(alnsg, 2))*PI/6.0;
        voltonumblo = 1./((PI/6.)*(std::pow(dgnumlo, 3))*std::exp(4.5*std::pow(alnsg, 2)));
        voltonumbhi = 1./((PI/6.)*(std::pow(dgnumhi, 3))*std::exp(4.5*std::pow(alnsg, 2)));
        v2nmin      = voltonumbhi;
        v2nmax      = voltonumblo;
        dgnxx       = dgnumhi;
        dgnyy       = dgnumlo;

        top_lev = 1;
        for (auto k = top_lev; k <= nlev; ++k) {
           for (auto i = 1; i <= ncol; ++i) {
              drv_a = dryvol_a(i,k);
              num_a0 = mode_num(i,k);
              num_a = std::max(0.0, num_a0);

              if (drv_a > 0.0) {
                 if (num_a <= drv_a*v2nmin)
                    dgncur_a(i,k) = dgnxx;
                 else if (num_a >= drv_a*v2nmax)
                    dgncur_a(i,k) = dgnyy;
                 else
                    dgncur_a(i,k) = std::pow((drv_a/(dumfac*num_a)), 1./3.);
              }
           }
        }
        parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int i, int k) {
          dgnum_m(i,k,n) = dgncur_a(i,k);
        });
     }
   }

   //
   // calculates aerosol sw radiative properties
   //
   void  modal_aero_sw(int list_idx, real dt, int nnite,
                       const int1d& idxnite, const bool& is_cmip6_volc,
                       const real2d& pdeldry, const real2d& pmid, const real2d& temperature, const real2d& qt,
                       const real2d& ext_cmip6_sw, const int1d& trop_level,
                       const real3d& tauxar, const real3d& wa, const real3d& ga,
                       const real3d& fa, const real2d& clear_rh) {
     // local variables
     int ilev_tropp;
     real2d mass("mass",ncol,nlev);               // layer mass
     real2d air_density("air_density",ncol,nlev);  // (kg/m3)

     real2d specmmr("specmmr",ncol,nlev);        // species mass mixing ratio
     real1d specrefindex_real("specrefindex_real",nswbands);
     real1d specrefindex_im("specrefindex_im",nswbands);
     std::string spectype;          // species type
     real hygro_aer;
     real volf;

     real2d dgnumwet("dgnumwet", ncol,nlev);     // number mode wet diameter
     real2d qaerwat("qaerwat", ncol,nlev);      // aerosol water (g/g)
     real2d wetdens("wetdens", ncol, nlev);

     real sigma_logr_aer;         // geometric standard deviation of number distribution
     real2d radsurf("radsurf", ncol, nlev);      // aerosol surface mode radius
     real2d logradsurf("logradsurf", ncol, nlev); // log(aerosol surface mode radius)
     real3d cheb("cheb", ncoef, ncol, nlev);

     real1d refr("refr", ncol);        // real part of refractive index
     real1d refi("refi", ncol);        // imaginary part of refractive index
     real1d crefin_real("crefin_real",ncol);
     real1d crefin_im("crefin_im",ncol);

     real2d refitabsw("refitabsw",prefi,nswbands);
     real2d refrtabsw("refrtabsw",prefr,nswbands);

     real4d extpsw("extpsw",ncoef,prefr,prefi,nswbands);
     real4d abspsw("abspsw",ncoef,prefr,prefi,nswbands);
     real4d asmpsw("asmpsw",ncoef,prefr,prefi,nswbands);

     real1d vol("vol", ncol);            // volume concentration of aerosol specie (m3/kg)
     real1d dryvol("dryvol", ncol);      // volume concentration of aerosol mode (m3/kg)
     real1d watervol("watervol", ncol);  // volume concentration of water in each mode (m3/kg)
     real1d wetvol("wetvol", ncol);      // volume concentration of wet mode (m3/kg)

     int1d  itab("itab", ncol), jtab("jtab", ncol);
     real1d ttab("ttab",ncol), utab("utab", ncol);
     real2d cext("cext", ncol, ncoef), cabs("cabs", ncol,ncoef), casm("casm",ncol,ncoef);
     real1d pext("pext",ncol);     // parameterized specific extinction (m2/kg)
     real1d specpext("specpext",ncol); // specific extinction (m2/kg)
     real1d dopaer("dopaer",ncol);   // aerosol optical depth in layer
     real1d pabs("pabs",ncol);     // parameterized specific absorption (m2/kg)
     real1d pasm("pasm",ncol);     // parameterized asymmetry factor
     real1d palb("palb",ncol);     // parameterized single scattering albedo

     // error code
     int nerr_dopaer;

     // initialize output variables
     yakl::memset(tauxar, 0.);
     yakl::memset(wa, 0.);
     yakl::memset(ga, 0.);
     yakl::memset(fa, 0.);

     // zero'th layer does not contain aerosol
     parallel_for(SimpleBounds<2>(ncol, nswbands), YAKL_LAMBDA (int i, int ibnd) {
        wa(i,1,ibnd) = 0.925;
        ga(i,1,ibnd) = 0.850;
        fa(i,1,ibnd) = 0.7225;
     });

     parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int i, int k) {
        mass(i,k)        = pdeldry(i,k)*rgas;
        air_density(i,k) = pmid(i,k)/(rv*temperature(i,k));
     });

     // Calculate aerosol size distribution parameters and aerosol water uptake
     if (clim_modal_aero) {
        // radiation diagnostics are not supported for prescribed aerosols cases
        if(list_idx != 0) {
           amrex::Print() << "Radiation diagnostic calls are not supported for prescribed aerosols\n";
           exit(0);
        }
        // diagnostic aerosol size calculations
        modal_aero_calcsize_diag(list_idx, dgnumdry_m);
     }

     // clear_rh provides alternate estimate non-cloudy relative humidity
     ModalAeroWateruptake::modal_aero_wateruptake_dr(list_idx, ncol, nlev, nmodes, top_lev,
                                                     mam_consti, qt, temperature, pmid, dgnumdry_m, dgnumwet_m,
                                                     qaerwat_m, wetdens_m, clear_rh);

     // loop over all aerosol modes
     for(auto m = 1; m <= nmodes; ++m) {
        parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev) {
          dgnumwet(icol,ilev) = dgnumwet_m(icol,ilev,m);
          qaerwat(icol,ilev)  = qaerwat_m(icol,ilev,m);
        });

        // get mode properties
        mam_consti.get_mode_props(list_idx, m-1, sigma_logr_aer, refrtabsw,
                                  refitabsw, extpsw, abspsw, asmpsw);

        // get mode info
        mam_consti.get_mode_nspec(list_idx, m-1, nspec);

        // calc size parameter for all columns
        modal_size_parameters(sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb);

        for(auto isw = 1; isw <= nswbands; isw++) {
           for(auto k = top_lev; k <= nlev; ++k) {
              // form bulk refractive index
              yakl::memset(crefin_real, 0.);
              yakl::memset(crefin_im, 0.);
              yakl::memset(dryvol, 0.);

              // aerosol species loop
              for(auto l = 0; l < nspec; ++l) {
                 real specdens = 1.0e-20;
                 mam_consti.rad_cnst_get_mam_mmr_by_idx(list_idx, m-1, l, "a", specmmr);
                 mam_consti.get_mam_props(list_idx, m-1, l, specdens, spectype, hygro_aer,
                                          specrefindex_real, specrefindex_im);
                 for(auto i = 1; i <= ncol; ++i) {
                    vol(i)          = specmmr(i,k)/specdens;
                    dryvol(i)       = dryvol(i) + vol(i);
                    crefin_real(i)  = crefin_real(i) + vol(i)*specrefindex_real(isw);
                    crefin_im(i)    = crefin_im(i) + vol(i)*specrefindex_im(isw);
                 }
              }

              for(auto i = 1; i <= ncol; ++i) {
                 watervol(i) = qaerwat(i,k)/rhoh2o;
                 wetvol(i) = watervol(i) + dryvol(i);
                 if (watervol(i) < 0.) {
                    if (std::abs(watervol(i)) > 1.e-1*wetvol(i)) {
                       amrex::Print() << "watervol,wetvol=" << watervol(i) << "; " << wetvol(i) << std::endl;
                    }
                    watervol(i) = 0.;
                    wetvol(i) = dryvol(i);
                 }
                 // volume mixing
                 crefin_real(i) = crefin_real(i)+watervol(i)*crefwswr(isw);
                 crefin_real(i) = crefin_real(i)/std::max(wetvol(i),1.e-60);
                 crefin_im(i)   = crefin_im(i)+watervol(i)*crefwswi(isw);
                 crefin_im(i)   = crefin_im(i)/std::max(wetvol(i),1.e-60);
                 refr(i)        = crefin_real(i);
                 refi(i)        = std::abs(crefin_im(i));
              }

              // interpolate coefficients linear in refractive index
              // first call calcs itab,jtab,ttab,utab
              real3d extpswr("extpswr", ncoef, prefr, prefi);
              real3d abspswr("abspswr", ncoef, prefr, prefi);
              real3d asmpswr("asmpswr", ncoef, prefr, prefi);
              real1d refitabswr("refitabswr", prefi);
              real1d refrtabswr("refrtabswr", prefr);

              parallel_for(SimpleBounds<3>(ncoef,prefr,prefi), YAKL_LAMBDA (int icoef, int irefr, int irefi) {
                 extpswr(icoef,irefr,irefi) = extpsw(icoef,irefr,irefi,isw);
                 abspswr(icoef,irefr,irefi) = abspsw(icoef,irefr,irefi,isw);
                 asmpswr(icoef,irefr,irefi) = asmpsw(icoef,irefr,irefi,isw);
                 refitabswr(irefi)          = refitabsw(irefi,isw);
                 refrtabswr(irefr)          = refrtabsw(irefr,isw);
              });

              yakl::memset(itab, 0);
              binterp(extpswr, ncol, ncoef, prefr, prefi, refr, refi,
                      refrtabswr, refitabswr, itab, jtab, ttab, utab, cext);
              binterp(abspswr, ncol, ncoef, prefr, prefi, refr, refi,
                      refrtabswr, refitabswr, itab, jtab, ttab, utab, cabs);
              binterp(asmpswr, ncol, ncoef, prefr, prefi, refr, refi,
                      refrtabswr, refitabswr, itab, jtab, ttab, utab, casm);

              // parameterized optical properties
              for(auto i=1; i <= ncol; ++i) {
                 if (logradsurf(i,k) <= xrmax) {
                    pext(i) = 0.5*cext(i,1);
                    for(auto nc = 2; nc <= ncoef; ++nc) {
                       pext(i) = pext(i) + cheb(nc,i,k)*cext(i,nc);
                    }
                    pext(i) = std::exp(pext(i));
                 } else {
                    pext(i) = 1.5/(radsurf(i,k)*rhoh2o);  // geometric optics
                 }

                 // convert from m2/kg water to m2/kg aerosol
                 specpext(i) = pext(i);
                 pext(i) = pext(i)*wetvol(i)*rhoh2o;
                 pabs(i) = 0.5*cabs(i,1);
                 pasm(i) = 0.5*casm(i,1);
                 for(auto nc = 2; nc < ncoef; ++nc) {
                    pabs(i) = pabs(i) + cheb(nc,i,k)*cabs(i,nc);
                    pasm(i) = pasm(i) + cheb(nc,i,k)*casm(i,nc);
                 }
                 pabs(i) = pabs(i)*wetvol(i)*rhoh2o;
                 pabs(i) = std::max(0.,pabs(i));
                 pabs(i) = std::min(pext(i),pabs(i));
                 palb(i) = 1.-pabs(i)/std::max(pext(i),1.e-40);
                 dopaer(i) = pext(i)*mass(i,k);
              }

              for (auto i = 1; i <= ncol; ++i) {
                 if ((dopaer(i) <= -1.e-10) || (dopaer(i) >= 30.)) {
                    if (dopaer(i) <= -1.e-10)
                       amrex::Print() << "ERROR: Negative aerosol optical depth in this layer.\n";
                    else {
                       // reset to the bound value
                       dopaer(i) = 25.0;
                       amrex::Print() << "WARNING: Aerosol optical depth is unreasonably high in this layer.\n";
                    }

                    for(auto l = 1; l < nspec; ++l) {
                       real specdens = 1.0e-20;
                       mam_consti.rad_cnst_get_mam_mmr_by_idx(list_idx, m-1, l, "a", specmmr);
                       mam_consti.get_mam_props_sw(list_idx, m-1, l, specdens, specrefindex_real, specrefindex_im);
                       volf = specmmr(i,k)/specdens;
                    }
                    nerr_dopaer = nerr_dopaer + 1;
                    if (dopaer(i) < -1.e-10) {
                       amrex::Print() << "*** halting after nerr_dopaer\n";
                       exit(0);
                    }
               }

              for(auto i=1; i <= ncol; ++i) {
                tauxar(i,k,isw) = tauxar(i,k,isw) + dopaer(i);
                wa(i,k,isw)     = wa(i,k,isw)     + dopaer(i)*palb(i);
                ga(i,k,isw)     = ga(i,k,isw)     + dopaer(i)*palb(i)*pasm(i);
                fa(i,k,isw)     = fa(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)*pasm(i);
              }
            }
          } //end do ! nlev
        } //end do ! swbands
     } // nmodes
  } // modal_aero_sw

  //
  // calcualtes aerosol lw radiative properties
  //
  void  modal_aero_lw(int list_idx, real dt,
                      const real2d& pdeldry, const real2d& pmid, const real2d& temperature, const real2d& qt,
                      const real3d& tauxar, const real2d& clear_rh) {
     real   sigma_logr_aer;          // geometric standard deviation of number distribution
     real   alnsg_amode;
     real1d xrad("xrad", ncol);
     real3d cheby("cheby", ncoef, ncol, nlev);  // chebychef polynomials
     real2d mass("mass",ncol,nlev); // layer mass

     //real2d specmmr(:,:)        ! species mass mixing ratio
     real specdens;            // species density (kg/m3)
     real1d specrefrindex("specrefrindex", ncol);     // species refractive index
     real1d specrefiindex("specrefiindex", ncol);

     real1d vol("vol",ncol);            // volume concentration of aerosol specie (m3/kg)
     real1d dryvol("dryvol",ncol);      // volume concentration of aerosol mode (m3/kg)
     real1d wetvol("wetvol",ncol);      // volume concentration of wet mode (m3/kg)
     real1d watervol("watervol",ncol);  // volume concentration of water in each mode (m3/kg)
     real1d refr("refr",ncol);          // real part of refractive index
     real1d refi("refi",ncol);          // imaginary part of refractive index
     real1d crefinr("crefinr", ncol);
     real1d crefini("crefini", ncol);

     real2d refrtablw("refrtablw",prefr,nlwbands);   // table of real refractive indices for aerosols
     real2d refitablw("refitablw",prefi,nlwbands);   // table of imag refractive indices for aerosols
     real4d absplw("absplw",ncoef,prefr,prefi,nlwbands);             // specific absorption

     int1d itab("itab",ncol), jtab("jtab",ncol);
     real1d ttab("ttab",ncol), utab("utab",ncol);
     real2d cabs("cabs",ncol,ncoef);
     real1d pabs("pabs",ncol);        // parameterized specific absorption (m2/kg)
     real1d dopaer("dopaer",ncol);    // aerosol optical depth in layer

     real2d specmmr;
     real2d dgnumwet;  // wet number mode diameter (m)
     real2d qaerwat;  // aerosol water (g/g)

     constexpr int nerrmax_dopaer=1000;
     int nerr_dopaer = 0;
     real volf;   // volume fraction of insoluble aerosol

     // initialize output variables
     yakl::memset(tauxar, 0.);

     // dry mass in each cell
     parallel_for (SimpleBounds<2> (ncol, nlev), YAKL_LAMBDA (int icol, int ilev) {
       //mass(icol,ilev) = pdeldry(icol,ilev)*rga;
     });

     // Calculate aerosol size distribution parameters and aerosol water uptake
     if (clim_modal_aero) {   // For prescribed aerosol codes
        //radiation diagnostics are not supported for prescribed aerosols cases
        if(list_idx != 0)
          amrex::Print() << "Radiation diagnostic calls are not supported for prescribed aerosols\n";
         // diagnostic aerosol size calculations
         modal_aero_calcsize_diag(list_idx, dgnumdry_m);
     }

     // clear_rh provides alternate estimate non-cloudy relative humidity
     ModalAeroWateruptake::modal_aero_wateruptake_dr(list_idx, ncol, nlev, nmodes, top_lev,
                                                    mam_consti, qt, temperature, pmid, dgnumdry_m, dgnumwet_m,
                                                    qaerwat_m, wetdens_m, clear_rh);

    for(auto m = 0; m < nmodes; ++m) {
      parallel_for (SimpleBounds<2> (ncol, nlev), YAKL_LAMBDA (int icol, int ilev) {
        dgnumwet(icol, ilev) = dgnumwet_m(icol,ilev,m);
        qaerwat(icol, ilev)  = qaerwat_m(icol,ilev,m);
      });

      // get mode properties
      mam_consti.get_mode_props(list_idx, m, sigma_logr_aer,
                                refrtablw, refitablw, absplw);
      // get mode info
      mam_consti.get_mode_nspec(list_idx, m, nspec);

      // calc size parameter for all columns
      // this is the same calculation that's done in modal_size_parameters, but there
      // some intermediate results are saved and the chebyshev polynomials are stored
      // in a array with different index order.  Could be unified.
      top_lev = 1;
      for(auto k = top_lev; k > nlev; ++k) {
         for(auto i = 1; i <= ncol; ++i) {
            alnsg_amode = std::log( sigma_logr_aer );
            // convert from number diameter to surface area
            xrad(i) = std::log(0.5*dgnumwet(i,k)) + 2.0*alnsg_amode*alnsg_amode;
            // normalize size parameter
            xrad(i) = std::max(xrad(i), xrmin);
            xrad(i) = std::min(xrad(i), xrmax);
            xrad(i) = (2*xrad(i)-xrmax-xrmin)/(xrmax-xrmin);
            // chebyshev polynomials
            cheby(1,i,k) = 1.0;
            cheby(2,i,k) = xrad(i);
            for(auto nc = 3; nc < ncoef; ++nc) {
               cheby(nc,i,k) = 2.0*xrad(i)*cheby(nc-1,i,k)-cheby(nc-2,i,k);
            }
         }
      }

      int nlwbands; // add by xyuan to be compiled
      for(auto ilw = 1; ilw <= nlwbands; ++ilw) {
         for(auto k = top_lev; k > nlev; --k) {
            // form bulk refractive index. Use volume mixing for infrared
            yakl::memset(crefinr, 0.0);
            yakl::memset(crefini, 0.0);
            yakl::memset(dryvol, 0.0);
            // aerosol species loop
            for(auto l = 0; l < nspec; ++l) {
               mam_consti.rad_cnst_get_mam_mmr_by_idx(list_idx, m, l, "a", specmmr);
               mam_consti.get_mam_props_lw(list_idx, m, l, specdens, specrefrindex, specrefiindex);
               for(auto i = 1; i <= ncol; ++i) {
                  vol(i)     = specmmr(i,k)/specdens;
                  dryvol(i)  = dryvol(i) + vol(i);
                  crefinr(i) = crefinr(i) + vol(i)*specrefrindex(ilw);
                  crefini(i) = crefini(i) + vol(i)*specrefiindex(ilw);
               }
            }
            for(auto i = 1; i <= ncol; ++i) {
               watervol(i) = qaerwat(i,k)/rhoh2o;
               wetvol(i)   = watervol(i) + dryvol(i);
               if (watervol(i) < 0.0) {
                  if (abs(watervol(i)) > 1.e-1*wetvol(i)) {
                    amrex::Print() << "watervol is too large\n";
                  }
                  watervol(i) = 0.;
                  wetvol(i)   = dryvol(i);
               }
               crefinr(i) = crefinr(i) + watervol(i)*crefwlwr(ilw);
               crefini(i) = crefini(i) + watervol(i)*crefwlwi(ilw);
               if (wetvol(i) > 1.e-40) {
                 crefinr(i) = crefinr(i)/wetvol(i);
                 crefini(i) = crefini(i)/wetvol(i);
               }
               refr(i) = crefinr(i);
               refi(i) = crefini(i);
            }

            // interpolate coefficients linear in refractive index
            // first call calcs itab,jtab,ttab,utab
            yakl::memset(itab, 0);

            real3d absplwr("absplwr", ncoef, prefr, prefi);
            real1d refitablwr("refitablwr", prefi);
            real1d refrtablwr("refrtablwr", prefr);

            parallel_for(SimpleBounds<3>(ncoef,prefr,prefi) , YAKL_LAMBDA (int icoef, int irefr, int irefi) {
               absplwr(icoef,irefr,irefi) = absplw(icoef,irefr,irefi,ilw);
               refitablwr(irefi)          = refitablw(irefi,ilw);
               refrtablwr(irefr)          = refrtablw(irefr,ilw);
            });

            binterp(absplwr, ncol, ncoef, prefr, prefi, refr, refi,
                    refrtablwr, refitablwr, itab, jtab, ttab, utab, cabs);

            // parameterized optical properties
            for(auto i = 1; i <= ncol; ++i) {
               pabs(i) = 0.5*cabs(i,1);
               for(auto nc = 2; nc < ncoef; ++nc) {
                  pabs(i) = pabs(i) + cheby(nc,i,k)*cabs(i,nc);
               }
               pabs(i)   = pabs(i)*wetvol(i)*rhoh2o;
               pabs(i)   = std::max(0.,pabs(i));
               dopaer(i) = pabs(i)*mass(i,k);
            }

            for(auto i = 1; i <= ncol; ++i) {
               if ((dopaer(i) <= -1.e-10) || (dopaer(i) >= 20.)) {
                  if (dopaer(i) <= -1.e-10)
                    amrex::Print() <<  "ERROR: Negative aerosol optical depth in this layer.\n";
                  else
                    amrex::Print() << "WARNING: Aerosol optical depth is unreasonably high in this layer.\n";

                  for(auto l = 0; l < nspec; ++l) {
                    mam_consti.rad_cnst_get_mam_mmr_by_idx(list_idx, m, l, "a", specmmr);
                    mam_consti.get_mam_props_lw(list_idx, m, l, specdens, specrefrindex, specrefiindex);
                    volf = specmmr(i,k)/specdens;
                  }

                  nerr_dopaer = nerr_dopaer + 1;
                  if (nerr_dopaer >= nerrmax_dopaer || dopaer(i) < -1.e-10) {
                     amrex::Print() <<  "*** halting after nerr_dopaer = " << nerr_dopaer;
                     exit(EXIT_FAILURE);
                  }
              }
            }
            for(auto i = 1; i <= ncol; ++i) tauxar(i,k,ilw) = tauxar(i,k,ilw) + dopaer(i);
         } // k = top_lev, nlev
       }   // nlwbands
     }  // m = 1, nmodes
  }     // modal_aero_lw
};
#endif // ERF_MAM4_AERO_H_

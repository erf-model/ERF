
// RCE 07.04.13:  Adapted from MIRAGE2 code / E3SM
#ifndef ERF_MODAL_AERO_WATERUPTAKE_H_
#define ERF_MODAL_AERO_WATERUPTAKE_H_

#include <string>
#include <vector>
#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_TableData.H>
#include <AMReX_MultiFabUtil.H>
#include "ERF_Constants.H"
#include "Aero_rad_props.H"
#include "Cloud_rad_props.H"
#include "Mam4_constituents.H"

using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

class ModalAeroWateruptake {
  public:
   constexpr static real third = 1./3.;
   constexpr static real pi43  = PI*4.0/3.0;
   constexpr static real huge_real = std::numeric_limits<real>::max();
   constexpr static int imax = 200;
  public:
   // The calculation has been enabled for diagnostic mode lists
   // via optional arguments.  If the list_idx arg is present then
   // all the optional args must be present.
   static
   void modal_aero_wateruptake_dr(int list_idx, int ncol, int nlev, int nmodes, int top_lev,
                                  MamConstituents& consti, const real3d& dgnumdry_m,
                                  const real3d& dgnumwet_m, const real3d& qaerwat_m,
                                  const real3d& wetdens_m, const real2d& clear_rh_in) {
   int itim_old;
   int nspec;

   real2d h2ommr("h2ommr",ncol,nlev); // specific humidity
   real2d t("temperature",ncol,nlev);      // temperatures (K)
   real2d pmid("pmid",ncol,nlev);   // layer pressure (Pa)
   real2d raer("raer",ncol,nlev);   // aerosol species MRs (kg/kg and #/kg)

   real2d cldn("cldn",ncol,nlev);               // layer cloud fraction (0-1)
   real3d dgncur_a("dgncur",ncol,nlev,nmodes);
   real3d dgncur_awet("dgncur_awet",ncol,nlev,nmodes);
   real3d wetdens("wetdens",ncol,nlev,nmodes);
   real3d qaerwat("qaerwat",ncol,nlev,nmodes);

   real2d dryvolmr("dryvolmr",ncol,nlev);          //volume MR for aerosol mode (m3/kg)
   real specdens;
   real spechygro, spechygro_1;
   real duma, dumb;
   real sigmag;
   real alnsg;
   real v2ncur_a;
   real drydens;                     // dry particle density  (kg/m^3)
   real2d rh("rh",ncol,nlev);        // relative humidity (0-1)

   real1d es("es",ncol);             // saturation vapor pressure
   real1d qs("qs",ncol);             // saturation specific humidity
   real cldn_thresh;
   real2d aerosol_water("aerosol_water",ncol,nlev); //sum of aerosol water (wat_a1 + wat_a2 + wat_a3 + wat_a4)
   bool history_aerosol;      // Output the MAM aerosol variables and tendencies
   bool history_verbose;      // produce verbose history output
   bool compute_wetdens;

   std::string trnum;       // used to hold mode number (as characters)

   real3d naer("naer",ncol,nlev,nmodes);      // aerosol number MR (bounded!) (#/kg-air)
   real3d dryvol("dryvol",ncol,nlev,nmodes);    // single-particle-mean dry volume (m3)
   real3d drymass("drymass",ncol,nlev,nmodes);   // single-particle-mean dry mass  (kg)
   real3d dryrad("dryrad",ncol,nlev,nmodes);    // dry volume mean radius of aerosol (m)
   real3d wetrad("wetrad",ncol,nlev,nmodes);    // wet radius of aerosol (m)
   real3d wetvol("wetvol",ncol,nlev,nmodes);    // single-particle-mean wet volume (m3)
   real3d wtrvol("wtrvol",ncol,nlev,nmodes);    // single-particle-mean water volume in wet aerosol (m3)
   real1d rhcrystal("rhcrystal",nmodes);
   real1d rhdeliques("rhdeliques",nmodes);
   real1d specdens_1("specdens_1",nmodes);
   real3d maer("maer",ncol,nlev,nmodes);      // aerosol wet mass MR (including water) (kg/kg-air)
   real3d hygro("hygro",ncol,nlev,nmodes);     // volume-weighted mean hygroscopicity (--)

   yakl::memset(naer, huge_real);
   yakl::memset(dryvol, huge_real);
   yakl::memset(drymass, huge_real);
   yakl::memset(dryrad, huge_real);
   yakl::memset(wetrad, huge_real);
   yakl::memset(wetvol, huge_real);
   yakl::memset(wtrvol, huge_real);

   yakl::memset(rhcrystal, huge_real);
   yakl::memset(rhdeliques, huge_real);
   yakl::memset(specdens_1, huge_real);

   yakl::memset(maer, 0.0);
   yakl::memset(hygro, 0.);

   //by default set compute_wetdens to be true
   compute_wetdens = true;

   dgncur_a    = dgnumdry_m;
   dgncur_awet = dgnumwet_m;
   qaerwat     = qaerwat_m;
   wetdens     = wetdens_m;

   // retreive aerosol properties
   for (auto m = 1; m <= nmodes; ++m) {
      yakl::memset(dryvolmr, 0.);

      // get mode properties
      consti.get_mode_props(list_idx, m, sigmag,
                            rhcrystal(m), rhdeliques(m));

      // get mode species
      consti.get_mode_nspec(list_idx, m, nspec);

      for(auto l = 0; l < nspec; ++l) {
         // get species interstitial mixing ratio ('a')
         consti.rad_cnst_get_mam_mmr_by_idx(list_idx, m, l, "a", raer);

         consti.get_mam_density_aer(list_idx, m, l, specdens);
         consti.get_mam_hygro_aer(list_idx, m, l, spechygro);

         if (l == 1) {
            // save off these values to be used as defaults
            specdens_1(m)  = specdens;
            spechygro_1    = spechygro;
         }

         for(auto k = top_lev; k > nlev; --k) {
            for(auto i = 1; i <= ncol; ++i) {
               duma          = raer(i,k);
               maer(i,k,m)   = maer(i,k,m) + duma;
               dumb          = duma/specdens;
               dryvolmr(i,k) = dryvolmr(i,k) + dumb;
               hygro(i,k,m)  = hygro(i,k,m) + dumb*spechygro;
            } // i = 1, ncol
         } // k = top_lev, pver
      } // l = 1, nspec

      alnsg = log(sigmag);

      for(auto k = top_lev; k >= nlev; --k) {
         for(auto i = 0; i <= ncol; ++i) {

            if (dryvolmr(i,k) > 1.0e-30) {
               hygro(i,k,m) = hygro(i,k,m)/dryvolmr(i,k);
            } else {
               hygro(i,k,m) = spechygro_1;
            }

            // dry aerosol properties
            v2ncur_a = 1. / ( (PI/6.)*(pow(dgncur_a(i,k,m),3.))*exp(4.5*pow(alnsg, 2.)) );
            // naer = aerosol number (#/kg)
            naer(i,k,m) = dryvolmr(i,k)*v2ncur_a;

            // compute mean (1 particle) dry volume and mass for each mode
            // old coding is replaced because the new (1/v2ncur_a) is equal to
            // the mean particle volume
            // also moletomass forces maer >= 1.0e-30, so (maer/dryvolmr)
            // should never cause problems (but check for maer < 1.0e-31 anyway)
            if (maer(i,k,m) > 1.0e-31) {
               drydens = maer(i,k,m)/dryvolmr(i,k);
            } else {
               drydens = 1.0;
            }

            dryvol(i,k,m)   = 1.0/v2ncur_a;
            drymass(i,k,m)  = drydens*dryvol(i,k,m);
            dryrad(i,k,m)   = std::pow((dryvol(i,k,m)/pi43), third);

         } // i = 1, ncol
      } // k = top_lev, pver
   } // modes

   // specify clear air relative humidity
   bool has_clear_rh = true;
   if (has_clear_rh) {
      // use input relative humidity
      // check that values are reasonable and apply upper limit
      for(auto k = top_lev; k <= nlev; --k) {
         for(auto i = 1; i <= ncol; ++i) {
            rh(i,k) = clear_rh_in(i,k);
            if ( rh(i,k) < 0 ) {
               amrex::Print() << "modal_aero_wateruptake_dr: clear_rh_in is negative - rh:" << rh(i,k);
               exit(EXIT_FAILURE);
            }
            // limit RH to 98% to be consistent with behavior when clear_rh_in is not provided
            rh(i,k) = std::min(rh(i,k), 0.98);
         } // i
      }

    } else {

      // estimate clear air relative humidity using cloud fraction
      //h2ommr => state%q(:,:,1)
      // t      => state%t
      // pmid   => state%pmid
#if 0
      itim_old    =  pbuf_old_tim_idx()
      pbuf_get_field(pbuf, cld_idx, cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) );

      for(auto k = top_lev; k > pver; --k) {
         qsat_water(t(:ncol,k), pmid(:ncol,k), es(:ncol), qs(:ncol));
         for(auto i = 0; i < ncol; ++i) {
            if (qs(i) > h2ommr(i,k)) {
               rh(i,k) = h2ommr(i,k)/qs(i);
            } else {
               rh(i,k) = 0.98;
            }

            rh(i,k) = max(rh(i,k), 0.0);
            rh(i,k) = min(rh(i,k), 0.98);
            if(pergro_mods) {
               cldn_thresh = 0.9998;
            }
            else {
               cldn_thresh = 1.0; //original code
            }

            if (cldn(i,k) < cldn_thresh) {
               rh(i,k) = (rh(i,k) - cldn(i,k)) / (1.0 - cldn(i,k));  // RH of clear portion
            }

            rh(i,k) = max(rh(i,k), 0.0);
         } // i = 1, ncol
      } // k = top_lev, pver
#endif
   }

   // compute aerosol wet radius and aerosol water
   modal_aero_wateruptake_sub(ncol, nlev, nmodes, top_lev,
                              rhcrystal, rhdeliques, dryrad,
                              hygro, rh, dryvol,
                              wetrad, wetvol, wtrvol);

   for(auto m = 1; m <= nmodes; ++m) {
      for(auto k = top_lev; k > nlev; --k) {
         for(auto i = 1; i <= ncol; ++i) {

            dgncur_awet(i,k,m) = dgncur_a(i,k,m) * (wetrad(i,k,m)/dryrad(i,k,m));
            qaerwat(i,k,m)     = rhoh2o*naer(i,k,m)*wtrvol(i,k,m);

            // compute aerosol wet density (kg/m3)
            if(compute_wetdens) {
               if (wetvol(i,k,m) > 1.0e-30) {
                  wetdens(i,k,m) = (drymass(i,k,m) + rhoh2o*wtrvol(i,k,m))/wetvol(i,k,m);
               }
               else {
                  wetdens(i,k,m) = specdens_1(m);
               }
            }
         }  //i = 1, ncol
      } // k = top_lev, pver
     } // m = 1, nmodes
  }

  // Purpose: Compute aerosol wet radius
  static
  void modal_aero_wateruptake_sub(int ncol, int nlev, int nmodes, int top_lev,
                                  const real1d& rhcrystal, const real1d& rhdeliques, const real3d& dryrad,
                                  const real3d& hygro, const real2d& rh, const real3d& dryvol,
                                  const real3d& wetrad, const real3d& wetvol, const real3d& wtrvol) {
   // loop over all aerosol modes
   for(auto m = 1; m <= nmodes; ++m) {
      real hystfac = 1.0 / std::max(1.0e-5, (rhdeliques(m) - rhcrystal(m)));
      for(auto k = top_lev; k > nlev; --k) {
         for(auto i = 1; i <= ncol; ++i) {

            // compute wet radius for each mode
            modal_aero_kohler(dryrad(i,k,m), hygro(i,k,m), rh(i,k), wetrad(i,k,m));

            wetrad(i,k,m) = std::max(wetrad(i,k,m), dryrad(i,k,m));
            wetvol(i,k,m) = pi43*std::pow(wetrad(i,k,m), 3.);
            wetvol(i,k,m) = std::max(wetvol(i,k,m), dryvol(i,k,m));
            wtrvol(i,k,m) = wetvol(i,k,m) - dryvol(i,k,m);
            wtrvol(i,k,m) = std::max(wtrvol(i,k,m), 0.0);

            // apply simple treatment of deliquesence/crystallization hysteresis
            // for rhcrystal < rh < rhdeliques, aerosol water is a fraction of
            // the "upper curve" value, and the fraction is a linear function of rh
            if (rh(i,k) < rhcrystal(m)) {
               wetrad(i,k,m) = dryrad(i,k,m);
               wetvol(i,k,m) = dryvol(i,k,m);
               wtrvol(i,k,m) = 0.0;
            }
            else if (rh(i,k) < rhdeliques(m)) {
               wtrvol(i,k,m) = wtrvol(i,k,m)*hystfac*(rh(i,k) - rhcrystal(m));
               wtrvol(i,k,m) = std::max(wtrvol(i,k,m), 0.0);
               wetvol(i,k,m) = dryvol(i,k,m) + wtrvol(i,k,m);
               wetrad(i,k,m) = std::pow(wetvol(i,k,m)/pi43, 1./3.);
            }
         } // columns
      } // levels
    } // modes
  }

  // calculates equlibrium radius r of haze droplets as function of
  // dry particle mass and relative humidity s using kohler solution
  // given in pruppacher and klett (eqn 6-35)

  // for multiple aerosol types, assumes an internal mixture of aerosols
  YAKL_INLINE
  static void modal_aero_kohler(const real& rdry_in, const real& hygro, const real& s,
                                real& rwet_out) {
     real a, b;
     real p40, p41, p42, p43; // coefficients of polynomial
     real p30, p31, p32;                  // coefficients of polynomial
     real p, r3, r4, r;            // wet radius (microns)
     real rdry;    // radius of dry particle (microns)
     real ss;                     // relative humidity (1 = saturated)
     real slog;    // log relative humidity
     real vol;      // total volume of particle (microns**3)
     real xi, xr;

     real1d cx4r("cx4r",4), cx3r("cx3r",3);
     real1d cx4i("cx4i",4), cx3i("cx3i",3);

     const real eps = 1.e-4;
     const real mw = 18.;
     const real pi = 3.14159;
     const real rhow = 1.;
     const real surften = 76.;
     const real tair = 273.;
     const real third = 1./3.;
     const real ugascon = 8.3e7;

     //effect of organics on surface tension is neglected
     a=2.e4*mw*surften/(ugascon*tair*rhow);

     rdry = rdry_in*1.0e6;   // convert (m) to (microns)
     vol  = std::pow(rdry, 3.);          // vol is r**3, not volume
     b    = vol*hygro;

     //quartic
      ss   = std::min(s,1.-eps);
      ss   = std::max(ss,1.e-10);
      slog = std::log(ss);
      p43  = -a/slog;
      p42  = 0.;
      p41  = b/slog-vol;
      p40  = a*vol/slog;

      // cubic for rh=1
      p32 = 0.;
      p31 = -b/a;
      p30 = -vol;

      if(vol <= 1.e-12) {
         r=rdry;
         return;
      }

      p = std::abs(p31)/(rdry*rdry);
      if(p < eps) {
         //approximate solution for small particles
         r=rdry*(1.+p*third/(1.-slog*rdry/a));
      }
      else {
         makoh_quartic(cx4r,cx4i,p43,p42,p41,p40);
         //find smallest real(r8) solution
         r = 1000.*rdry;
         auto nsol = 0;
         for(auto n=1; n<=4; ++n) {
            xr=cx4r(n);
            xi=cx4i(n);
            if(abs(xi) > abs(xr)*eps) continue;
            if(xr > r) continue;
            if(xr < rdry*(1.-eps)) continue;
            if(xr != xr) continue;
            r=xr;
            nsol=n;
         }

         if(nsol != 0) {
            amrex::Print() << "ccm kohlerc - no real(r8) solution found (quartic)\n";
            r=rdry;
         }
      }

      if(s > 1.-eps) {
         // save quartic solution at s=1-eps
         r4=r;
         //  cubic for rh=1
         p=abs(p31)/(rdry*rdry);
         if(p < eps) {
            r=rdry*(1.+p*third);
         }
         else {
            //makoh_cubic(cx3,p32,p31,p30,im);
            //find smallest real(r8) solution
            r=1000.*rdry;
            auto nsol = 0;
            for(auto n=1; n<=3; ++n) {
               xr = cx3r(n);
               xi = cx3i(n);
               if(abs(xi) > abs(xr)*eps) continue;
               if(xr > r) continue;
               if(xr < rdry*(1.-eps)) continue;
               if(xr != xr) continue;
               r=xr;
               nsol=n;
            }

            if(nsol == 0) {
               amrex::Print() << "ccm kohlerc - no real(r8) solution found (cubic)\n";
               r=rdry;
            }
         }
         r3=r;
         //now interpolate between quartic, cubic solutions
         r=(r4*(1.-s)+r3*(s-1.+eps))/eps;
       }

     // bound and convert from microns to m
     r = std::min(r,30.); // upper bound based on 1 day lifetime
     rwet_out = r*1.e-6;
   }

   //solves  x**3 + p2 x**2 + p1 x + p0 = 0
   //where p0, p1, p2 are real
   static
   void makoh_cubic(const real1d& cxr, const real1d& cxi,
                    const real& p2, const real& p1, const real& p0)
   {
      real q, r;
      real cradr, cyr, czr, cradi, cyi, czi;

      const real eps = 1.e-20;
      const real cwr = -0.5;
      const real cwi = sqrt(3);
      const real cwsqr = -0.5;
      const real cwsqi = -sqrt(3);

      if(p1 == 0.) {
         // completely insoluble particle
         cxr(1) = std::pow(-p0, 1./3.);
         cxr(2) = cxr(1);
         cxr(3) = cxr(1);
         cxi(1) = 0.;
         cxi(2) = 0.;
         cxi(3) = 0.;
      }
      else {
         q = p1/3.;
         r = p0/2.;
         auto r2q3 = r*r+q*q*q;
         if (r2q3 >= 0.) {
           cradr = sqrt(r2q3);
         }
         else {
           cradi = sqrt(abs(r2q3));
         }

         cyr = r-cradr;
         cyi = cradi;

         if (abs(cyr) > eps) cyr = std::pow(cyr, 1./3.);
         auto cq = q;
         czr = -cq/cyr;
         czi = -cq/cyi;

         cxr(1) = -cyr-czr;
         cxr(2) = -cwr*cyr-cwsqr*czr;
         cxr(3) = -cwsqr*cyr-cwr*czr;

         cxi(1) = -cyi-czi;
         cxi(2) = -cwi*cyr-cwsqi*czi;
         cxi(3) = -cwsqi*cyr-cwi*czi;
      }
    }

    // solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
    // where p0, p1, p2, p3 are real
    static
    void makoh_quartic(const real1d& cxr, const real1d& cxi,
                       const real& p3, const real& p2,  const real& p1, const real& p0)
    {
      real q, r;
      real cbr, cb0r, cb1r, cbi, cb0i, cb1i;
      real cradr, cyr, cradi, cyi;

      q = -p2*p2/36.+(p3*p1-4*p0)/12.;
      r = -std::pow(p2/6, 3.)+p2*(p3*p1-4*p0)/48.
          +(4*p0*p2-p0*p3*p3-p1*p1)/16;

      cradr = r*r+q*q*q;
      cradr = sqrt(cradr);

      cbr = r-cradr;
      if(cbr == 0.0) {
         //insoluble particle
         cxr(1) = std::pow(-p1, 1./3.);
         cxr(2) = cxr(1);
         cxr(3) = cxr(1);
         cxr(4) = cxr(1);
      }

      cbr = std::pow(cbr, 1./3.);

      cyr = -cbr+q/cbr+p2/6;

      cb0r = sqrt(cyr*cyr-p0);
      cb1r = (p3*cyr-p1)/(2*cb0r);

      cbr   = p3/2+cb1r;
      cradr = cbr*cbr-4*(cyr+cb0r);
      cradr = sqrt(cradr);
      cxr(1) = (-cbr+cradr)/2.;
      cxr(2) = (-cbr-cradr)/2.;

      cbr   = p3/2-cb1r;
      cradr = cbr*cbr-4*(cyr-cb0r);
      cradr = sqrt(cradr);
      cxr(3) = (-cbr+cradr)/2.;
      cxr(4) = (-cbr-cradr)/2.;
   }
};
#endif

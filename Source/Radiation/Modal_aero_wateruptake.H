
// RCE 07.04.13:  Adapted from MIRAGE2 code / E3SM
#ifndef ERF_MODAL_AERO_WATERUPTAKE_H_
#define ERF_MODAL_AERO_WATERUPTAKE_H_

#include <string>
#include <vector>
#include <memory>
#include <AMReX_GpuComplex.H>
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
   static
   void modal_aero_wateruptake_dr(int list_idx, int ncol, int nlev, int nmodes, int top_lev,
                                  MamConstituents& consti, const real2d& h2ommr, const real2d& t, const real2d& pmid, 
                                  const real3d& dgnumdry_m, const real3d& dgnumwet_m, const real3d& qaerwat_m,
                                  const real3d& wetdens_m, const real2d& clear_rh_in) {
   real2d raer("raer",ncol,nlev);     // aerosol species MRs (kg/kg and #/kg)
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
      consti.get_mode_props(list_idx, m-1, sigmag, rhcrystal(m), rhdeliques(m));

      // get mode species
      int nspec;
      consti.get_mode_nspec(list_idx, m-1, nspec);
      for(auto l = 0; l < nspec; ++l) {
         // get species interstitial mixing ratio ('a')
         consti.rad_cnst_get_mam_mmr_by_idx(list_idx, m-1, l, "a", raer);
         consti.get_mam_density_aer(list_idx, m-1, l, specdens);
         consti.get_mam_hygro_aer(list_idx, m-1, l, spechygro);

         if (l == 1) {
            // save off these values to be used as defaults
            specdens_1(m)  = specdens;
            spechygro_1    = spechygro;
         }

         top_lev = 1;
         for(auto k = top_lev; k <= nlev; ++k) {
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
      top_lev = 1;
      for(auto k = top_lev; k <= nlev; ++k) {
         for(auto i = 1; i <= ncol; ++i) {

            if (dryvolmr(i,k) > 1.0e-30) {
               hygro(i,k,m) = hygro(i,k,m)/dryvolmr(i,k);
            } else {
               hygro(i,k,m) = spechygro_1;
            }

            // dry aerosol properties
            v2ncur_a = 1. / ( (PI/6.)*(std::pow(dgncur_a(i,k,m),3.))*std::exp(4.5*std::pow(alnsg, 2.)) );
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
   bool has_clear_rh = false;
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
      for(auto k = top_lev; k <= nlev; ++k) {
         for(auto i = 1; i <= ncol; ++i) {
            WaterVaporSat::qsat(t(i,k), pmid(i,k), es(i), qs(i));
            if (qs(i) > h2ommr(i,k)) {
               rh(i,k) = h2ommr(i,k)/qs(i);
            } else {
               rh(i,k) = 0.98;
            }

            rh(i,k) = std::max(rh(i,k), 0.0);
            rh(i,k) = std::min(rh(i,k), 0.98);
            if (cldn(i,k) < 0.998) {
               rh(i,k) = (rh(i,k) - cldn(i,k)) / (1.0 - cldn(i,k));
            }
            rh(i,k) = std::max(rh(i,k), 0.0);
         } // i = 1, ncol
      } // k = top_lev, pver
   }

   // compute aerosol wet radius and aerosol water
   modal_aero_wateruptake_sub(ncol, nlev, nmodes, top_lev,
                              rhcrystal, rhdeliques, dryrad,
                              hygro, rh, dryvol,
                              wetrad, wetvol, wtrvol);

   for(auto m = 1; m <= nmodes; ++m) {
      for(auto k = top_lev; k <= nlev; ++k) {
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
      for(auto k = top_lev; k <= nlev; ++k) {
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
     const real eps = 1.e-4;
     const real mw = 18.;
     const real pi = 3.14159;
     const real rhow = 1.;
     const real surften = 76.;
     const real tair = 273.;
     const real third = 1./3.;
     const real ugascon = 8.3e7;

     //effect of organics on surface tension is neglected
     auto a=2.e4*mw*surften/(ugascon*tair*rhow);

     auto rdry = rdry_in*1.0e6;   // convert (m) to (microns)
     auto vol  = std::pow(rdry, 3.);          // vol is r**3, not volume
     auto b    = vol*hygro;

     //quartic
     auto ss   = std::min(s,1.-eps); // relative humidity
     ss   = std::max(ss,1.e-10);
     auto slog = std::log(ss); // log relative humidity
     auto p43  = -a/slog;
     auto p42  = 0.;
     auto p41  = b/slog-vol;
     auto p40  = a*vol/slog;

     // cubic for rh=1
     auto p32 = 0.;
     auto p31 = -b/a;
     auto p30 = -vol;

     real r, r3, r4;
     if(vol <= 1.e-12) {
        r=rdry;
        return;
     }

     auto p = std::abs(p31)/(rdry*rdry);
     if(p < eps) {
         //approximate solution for small particles
         r=rdry*(1.+p*third/(1.-slog*rdry/a));
      }
      else {
         amrex::GpuComplex<real> cx4[4];
         makoh_quartic(cx4,p43,p42,p41,p40);
         //find smallest real(r8) solution
         r = 1000.*rdry;
         auto nsol = 0;
         for(auto n=0; n<4; ++n) {
            auto xr=cx4[n].real();
            auto xi=cx4[n].imag();
            if(abs(xi) > abs(xr)*eps) continue;
            if(xr > r) continue;
            if(xr < rdry*(1.-eps)) continue;
            if(xr != xr) continue;
            r=xr;
            nsol=n;
         }
         if(nsol != 0) {
           amrex::Print() << "kohlerc: no real solution found(quartic), nsol= " << nsol << "\n";
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
            amrex::GpuComplex<real> cx3[3];
            makoh_cubic(cx3,p32,p31,p30);
            //find smallest real(r8) solution
            r=1000.*rdry;
            auto nsol = 0;
            for(auto n=0; n<3; ++n) {
               auto xr = cx3[n].real();
               auto xi = cx3[n].imag();
               if(abs(xi) > abs(xr)*eps) continue;
               if(xr > r) continue;
               if(xr < rdry*(1.-eps)) continue;
               if(xr != xr) continue;
               r=xr;
               nsol=n;
            }
            if(nsol == 0) {
              amrex::Print() << "kohlerc: no real solution found (cubic)\n";
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

    // solves  x**3 + p2 x**2 + p1 x + p0 = 0
    // where p0, p1, p2 are real
    static
    void makoh_cubic(amrex::GpuComplex<real> cx[], const real& p2, const real& p1, const real& p0) {
      const real eps = 1.0e-20;
      const real third=1./3.;

      auto ci = amrex::GpuComplex<real>(0., 1.);
      auto sqrt3 = amrex::GpuComplex<real>(std::sqrt(3.), 0.);
      auto cw   = 0.5*(-1.+ci*sqrt3);
      auto cwsq = 0.5*(-1.-ci*sqrt3);

      if(p1 == 0.) {
        // completely insoluble particle
         cx[0] = amrex::GpuComplex<real>(std::pow(-p0, third), 0.);
         cx[1] = cx[0];
         cx[2] = cx[1];
      } 
      else {
         auto q = amrex::GpuComplex<real>(p1/3., 0.);
         auto r = amrex::GpuComplex<real>(p0/2., 0.);
         auto crad = r*r+q*q*q;
         crad = amrex::sqrt(crad);

         auto cy = r-crad;
         if (amrex::abs(cy) > eps) cy = amrex::pow(cy, third);
         auto cq = q;
         auto cz = -cq/cy;

         cx[0] = -cy-cz;
         cx[1] = -cw*cy-cwsq*cz;
         cx[2] = -cwsq*cy-cw*cz;
      }
    }

    // solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
    // where p0, p1, p2, p3 are real
    static
    void makoh_quartic(amrex::GpuComplex<real> cx[], const real& p3, const real& p2, const real& p1, const real& p0) {
      const real third=1./3.;
      auto q = -p2*p2/36.+(p3*p1-4*p0)/12.;
      auto r = -std::pow((p2/6), 3.0) + p2*(p3*p1-4.*p0)/48.
             +(4.*p0*p2-p0*p3*p3-p1*p1)/16.;

      auto crad = amrex::GpuComplex<real>(r*r+q*q*q, 0.);
      crad = amrex::sqrt(crad);
      auto cb = r-crad;

      if (cb.real() == 0. && cb.imag() == 0.) {
         // insoluble particle
         cx[0] = amrex::pow(amrex::GpuComplex<real>(-p1, 0.0), third);
         cx[1] = cx[0];
         cx[2] = cx[0];
         cx[3] = cx[0];
      } 
      else {
         cb = amrex::pow(cb, third);
         auto cy = -cb+q/cb+p2/6.;
         auto cb0 = amrex::sqrt(cy*cy-p0);
         auto cb1 = (p3*cy-p1)/(2.*cb0);

         cb   = p3/2.+cb1;
         crad = cb*cb-4.*(cy+cb0);
         crad = amrex::sqrt(crad);
         cx[1] = (-cb+crad)/2.;
         cx[2] = (-cb-crad)/2.;
         cb    = p3/2.-cb1;
         crad  = cb*cb-4.*(cy-cb0);
         crad  = amrex::sqrt(crad);
         cx[3] = (-cb+crad)/2.;
         cx[4] = (-cb-crad)/2.;
     }
  }
};
#endif

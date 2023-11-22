
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
using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

class ModalAeroWateruptake {
  public:
   constexpr real third = 1./3.;
   constexpr real pi43  = pi*4.0/3.0;
   constexpr real huge_real = numeric_limits<real>::max();
   constexpr int imax = 200;

   // The calculation has been enabled for diagnostic mode lists
   // via optional arguments.  If the list_idx arg is present then
   // all the optional args must be present.
   static void modal_aero_wateruptake_dr(int list_idx_in, real3d& dgnumdry_m, real3d& dgnumwet_m,
                                         real3d& qaerwat_m, real3d& wetdens_m, real2d& clear_rh_in)
   int list_idx;  // radiative constituents list index

   int itim_old;
   int nmodes;
   int nspec;

   real2d h2ommr(:,:); // specific humidity
   real2d t(:,:);      // temperatures (K)
   real2d pmid(:,:);   // layer pressure (Pa)
   real2d raer(:,:);   // aerosol species MRs (kg/kg and #/kg)

   real2d cldn(:,:);               // layer cloud fraction (0-1)
   real3d dgncur_a(:,:,:);
   real3d dgncur_awet(:,:,:);
   real3d wetdens(:,:,:);
   real3d qaerwat(:,:,:);

   real2d dryvolmr("dryvolmr",pcols,pver);          //volume MR for aerosol mode (m3/kg)
   real specdens;
   real spechygro, spechygro_1;
   real duma, dumb;
   real sigmag;
   real alnsg;
   real v2ncur_a;
   real drydens;                     // dry particle density  (kg/m^3)
   real2d rh("rh",pcols,pver);        // relative humidity (0-1)

   real1d es("es",pcols);             // saturation vapor pressure
   real1d qs("qs",pcols);             // saturation specific humidity
   real cldn_thresh;
   real2d aerosol_water("aerosol_water",pcols,pver); //sum of aerosol water (wat_a1 + wat_a2 + wat_a3 + wat_a4)
   bool history_aerosol;      // Output the MAM aerosol variables and tendencies
   bool history_verbose;      // produce verbose history output
   bool compute_wetdens;

   std::string trnum;       // used to hold mode number (as characters)

   real3d naer("naer",pcols,pver,nmodes);      // aerosol number MR (bounded!) (#/kg-air)
   real3d dryvol("dryvol",pcols,pver,nmodes);    // single-particle-mean dry volume (m3)
   real3d drymass("drymass",pcols,pver,nmodes);   // single-particle-mean dry mass  (kg)
   real3d dryrad("dryrad",pcols,pver,nmodes);    // dry volume mean radius of aerosol (m)
   real3d wetrad("wetrad",pcol,pver,nmodes);    // wet radius of aerosol (m)
   real3d wetvol("wetvol",pcol,pver,nmodes);    // single-particle-mean wet volume (m3)
   real3d wtrvol("wtrvol",pcol,pver,nmodes);    // single-particle-mean water volume in wet aerosol (m3)
   real1d rhcrystal("rhcrystal",nmodes);
   real1d rhdeliques("rhdeliques",nmodes);
   real1d specdens_1("specdens_1",nmodes);
   real3d maer("maer",pcol,pver,nmodes);      // aerosol wet mass MR (including water) (kg/kg-air)
   real3d hygro("hygro",pcol,pver,nmodes);     // volume-weighted mean hygroscopicity (--)

   //ncol = state%ncol

   // determine default variables
  // call phys_getopts(history_aerosol_out = history_aerosol, &
  //                   history_verbose_out = history_verbose)

   list_idx = 0
   if (list_idx_in)
      list_idx = list_idx_in;

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
   yalk::memset(hygro, 0.);

   //by default set compute_wetdens to be true
   compute_wetdens = true;

   if (list_idx_in) {
      pbuf_get_field(pbuf, dgnum_idx,      dgncur_a );
      pbuf_get_field(pbuf, dgnumwet_idx,   dgncur_awet );
      pbuf_get_field(pbuf, wetdens_ap_idx, wetdens);
      pbuf_get_field(pbuf, qaerwat_idx,    qaerwat);
   } else {
      dgncur_a    = dgnumdry_m;
      dgncur_awet = dgnumwet_m;
      qaerwat     = qaerwat_m;
      if(wetdens_m) {
         if ( !wetdens_m) {
            amrex::PRINT << "modal_aero_wateruptake_dr called"
                         << " with list_idx_in but wetdens_m is not allocated\n";
         }
         wetdens = wetdens_m;
      } else {
         //set compute_wetdens to flase if wetdens is not present
         compute_wetdens = .false.
      }
   }

   // retreive aerosol properties

   for(auto m = 0; m < nmodes; ++m) {

      yakl::memset(dryvolmr, 0.);

      // get mode properties
      rad_cnst_get_mode_props(list_idx, m, sigmag=sigmag,
         rhcrystal=rhcrystal(m), rhdeliques=rhdeliques(m));

      // get mode info
      rad_cnst_get_info(list_idx, m, nspec=nspec);

      for(auto l = 0; l < nspec; ++l) {

         // get species interstitial mixing ratio ('a')
         rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, raer);
         rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, hygro_aer=spechygro);

         if (l == 1) {
            // save off these values to be used as defaults
            specdens_1(m)  = specdens;
            spechygro_1    = spechygro;
         }

         for(auto k = top_lev; k > pver; --k) {
            for(auto i=0; i<ncol; ++i) {
               duma          = raer(i,k);
               maer(i,k,m)   = maer(i,k,m) + duma;
               dumb          = duma/specdens;
               dryvolmr(i,k) = dryvolmr(i,k) + dumb;
               hygro(i,k,m)  = hygro(i,k,m) + dumb*spechygro;
            } // i = 1, ncol
         } // k = top_lev, pver
      } // l = 1, nspec

      alnsg = log(sigmag);

      for(auto k = top_lev; k > pver; --k) {
         for(auto i = 0; i < ncol; ++i) {

            if (dryvolmr(i,k) > 1.0e-30) {
               hygro(i,k,m) = hygro(i,k,m)/dryvolmr(i,k);
            } else {
               hygro(i,k,m) = spechygro_1;
            }

            // dry aerosol properties

            v2ncur_a = 1. / ( (pi/6.)*(dgncur_a(i,k,m)**3.)*exp(4.5*pow(alnsg, 2.)) );
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

   if (clear_rh_in) {

      // use input relative humidity
      // check that values are reasonable and apply upper limit
      for(auto k = top_lev; k < pver; --k) {
         for(auto i = 0; i < ncol; ++i) {
            rh(i,k) = clear_rh_in(i,k);
            if ( rh(i,k)<0 ) {
               write(iulog,*) 'modal_aero_wateruptake_dr: clear_rh_in is negative - rh:',rh(i,k),' k=',k
               call endrun('modal_aero_wateruptake_dr: clear_rh_in cannot be negative')
            }
            // limit RH to 98% to be consistent with behavior when clear_rh_in is not provided
            rh(i,k) = min(rh(i,k), 0.98);
         } // i
      }

    } else {

      // estimate clear air relative humidity using cloud fraction
      //h2ommr => state%q(:,:,1)
      // t      => state%t
      // pmid   => state%pmid

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

   } // if present(clear_rh_in)

   // compute aerosol wet radius and aerosol water

   modal_aero_wateruptake_sub(
      ncol, nmodes, rhcrystal, rhdeliques, dryrad,
      hygro, rh, dryvol, wetrad, wetvol,
      wtrvol);

   for(auto m = 0; m < nmodes; ++m) {

      for(auto k = top_lev; k > pver; --k) {
         for(auto i = 1; i < ncol; ++i) {

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
void modal_aero_wateruptake_sub(
   ncol, nmodes, real1d& rhcrystal, reaal1d& rhdeliques, real3d& dryrad,
   real3d& hygro, real2d& rh, real3d& dryvol, real3d& wetrad, real3d& wetvol,
   real3d& wtrvol) {
   // local variables

   real hystfac;      // working variable for hysteresis

   // loop over all aerosol modes
   for(auto m = 0; m < nmodes; ++m) {

      hystfac = 1.0 / max(1.0e-5, (rhdeliques(m) - rhcrystal(m)));

      for(auto k = top_lev; k > pver; --k) {
         for(auto i = 0; i < ncol; ++i) {

            // compute wet radius for each mode
            modal_aero_kohler(dryrad(i:i,k,m), hygro(i:i,k,m), rh(i:i,k), wetrad(i:i,k,m), 1);

            wetrad(i,k,m) = max(wetrad(i,k,m), dryrad(i,k,m));
            wetvol(i,k,m) = pi43*wetrad(i,k,m)**3;
            wetvol(i,k,m) = max(wetvol(i,k,m), dryvol(i,k,m));
            wtrvol(i,k,m) = wetvol(i,k,m) - dryvol(i,k,m);
            wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0);

            // apply simple treatment of deliquesence/crystallization hysteresis
            // for rhcrystal < rh < rhdeliques, aerosol water is a fraction of
            // the "upper curve" value, and the fraction is a linear function of rh
            if (rh(i,k) < rhcrystal(m)) {
               wetrad(i,k,m) = dryrad(i,k,m);
               wetvol(i,k,m) = dryvol(i,k,m);
               wtrvol(i,k,m) = 0.0;
            } else if (rh(i,k) < rhdeliques(m)) {
               wtrvol(i,k,m) = wtrvol(i,k,m)*hystfac*(rh(i,k) - rhcrystal(m));
               wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0);
               wetvol(i,k,m) = dryvol(i,k,m) + wtrvol(i,k,m);
               wetrad(i,k,m) = (wetvol(i,k,m)/pi43)**third;
            }

         } // columns
      } // levels

   } // modes
}

// calculates equlibrium radius r of haze droplets as function of
// dry particle mass and relative humidity s using kohler solution
// given in pruppacher and klett (eqn 6-35)

// for multiple aerosol types, assumes an internal mixture of aerosols
 void modal_aero_kohler(real1d& rdry_in, real1d& hygro, real1d& s,
                        real1d& rwet_out, int im) {
      real a, b;
      real1d p40("p40",imax), p41("p41",imax), p42("p42",imax), p43("p43",imax); // coefficients of polynomial
      real1d p30("p30",imax), p31("p31",imax), p32("p32",imax);                  // coefficients of polynomial
      real p;
      real r3, r4;
      real1d r("r",im);            // wet radius (microns)
      real1d rdry("rdry",imax);    // radius of dry particle (microns)
      real ss;                     // relative humidity (1 = saturated)
      real1d slog("slog",imax);    // log relative humidity
      real1d vol("vol",imax);      // total volume of particle (microns**3)
      real xi, xr;

      complex cx4(4,imax),cx3(3,imax);

      constexpr real eps = 1.e-4;
      constexpr real mw = 18.;
      constexpr real pi = 3.14159;
      constexpr real rhow = 1.;
      constexpr real surften = 76.;
      constexpr real tair = 273.;
      constexpr real third = 1./3.;
      constexpr real ugascon = 8.3e7;


     //effect of organics on surface tension is neglected
      a=2.e4*mw*surften/(ugascon*tair*rhow);

      for(auto i=0; i<im; ++i) {
           rdry(i) = rdry_in(i)*1.0e6;   // convert (m) to (microns)
           vol(i) = rdry(i)**3;          // vol is r**3, not volume
           b = vol(i)*hygro(i);

          //quartic
           ss=min(s(i),1.-eps);
           ss=max(ss,1.e-10;)
           slog(i)=log(ss);
           p43(i)=-a/slog(i);
           p42(i)=0.;
           p41(i)=b/slog(i)-vol(i);
           p40(i)=a*vol(i)/slog(i);

         // cubic for rh=1
           p32(i)=0.;
           p31(i)=-b/a;
           p30(i)=-vol(i);
      }


       do 100 i=1,im

        if(vol(i) <= 1.e-12) {
           r(i)=rdry(i);
           goto 100;
        }

        p=abs(p31(i))/(rdry(i)*rdry(i));
        if(p < eps) {
          //approximate solution for small particles
           r(i)=rdry(i)*(1.+p*third/(1.-slog(i)*rdry(i)/a));
        } else {
           makoh_quartic(cx4(1,i),p43(i),p42(i),p41(i),p40(i),1);
           //find smallest real(r8) solution
           r(i)=1000.*rdry(i);
           nsol=0;
           for(auto n=1; n<=4; ++n) {
              xr=real(cx4(n,i));
              xi=aimag(cx4(n,i));
              if(abs(xi) > abs(xr)*eps) continue;
              if(xr > r(i)) continue;
              if(xr < rdry(i)*(1._r8-eps)) continue
              if(xr != xr) continue
              r(i)=xr;
              nsol=n;
           }

           if(nsol != 0) {
              amrex::PRINT << "ccm kohlerc - no real(r8) solution found (quartic)\n";
              r(i)=rdry(i);
           }
        }

        if(s(i) > 1.-eps) {
           // save quartic solution at s=1-eps
           r4=r(i);
           //  cubic for rh=1
           p=abs(p31(i))/(rdry(i)*rdry(i));
           if(p < eps) {
              r(i)=rdry(i)*(1._r8+p*third);
           }
           else {
              makoh_cubic(cx3,p32,p31,p30,im);
             //find smallest real(r8) solution
              r(i)=1000.*rdry(i);
              nsol=0;
              for(auto n=1; n<=3; ++n) {
                 xr=real(cx3(n,i));
                 xi=aimag(cx3(n,i));
                 if(abs(xi) > abs(xr)*eps) continue;
                 if(xr > r(i)) continue;
                 if(xr < rdry(i)*(1.-eps)) continue;
                 if(xr != xr) continue;
                 r(i)=xr;
                 nsol=n;
              }

              if(nsol == 0) {
                 amrex::PRINT << "ccm kohlerc - no real(r8) solution found (cubic)\n";
                 r(i)=rdry(i);
              }
           }
           r3=r(i);
           //now interpolate between quartic, cubic solutions
           r(i)=(r4*(1.-s(i))+r3*(s(i)-1.+eps))/eps;
        }

  100 continue

     // bound and convert from microns to m
      for(auto i=0; i<im; ++i) {
         r(i) = min(r(i),30.); // upper bound based on 1 day lifetime
         rwet_out(i) = r(i)*1.e-6;
      }

      return;
   }


     //solves  x**3 + p2 x**2 + p1 x + p0 = 0
     //where p0, p1, p2 are real
     void makoh_cubic( cx, p2, p1, p0, im )
      //integer :: im
      //real(r8) :: p0(imx), p1(imx), p2(imx)
      //complex(r8) :: cx(3,imx)

      real eps, sqrt3, third;
      real1d q("q",imx), r("r",imx);
      complex ci, cq, crad(imx), cw, cwsq, cy(imx), cz(imx)

      constexpr real eps = 1.e-20;

      third=1./3.;

      ci=cmplx(0.,1.);

      sqrt3=sqrt(3.);
      cw=0.5*(-1+ci*sqrt3);
      cwsq=0.5*(-1-ci*sqrt3);

      for(auto i=0; i<im; ++i) {
      if(p1(i) == 0.) {
         // completely insoluble particle
         cx(1,i)=(-p0(i))**third;
         cx(2,i)=cx(1,i);
         cx(3,i)=cx(1,i);
      }
      else {
         q(i)=p1(i)/3.;
         r(i)=p0(i)/2.;
         crad(i)=r(i)*r(i)+q(i)*q(i)*q(i);
         crad(i)=sqrt(crad(i));

         cy(i)=r(i)-crad(i);
         if (abs(cy(i)) > eps) cy(i)=cy(i)**third;
         cq=q(i);
         cz(i)=-cq/cy(i);

         cx(1,i)=-cy(i)-cz(i);
         cx(2,i)=-cw*cy(i)-cwsq*cz(i);
         cx(3,i)=-cwsq*cy(i)-cw*cz(i);
      }
      }
      return
     }


      // solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
      // where p0, p1, p2, p3 are real
      void makoh_quartic( cx, p3, p2, p1, p0, im ) {

      constexpr int imx=200;
      //real(r8) :: p0(imx), p1(imx), p2(imx), p3(imx)
      //complex(r8) :: cx(4,imx)

      real third;
      real1d q("q",imx), r("r",imx);
      complex cb(imx), cb0(imx), cb1(imx);
      complex crad(imx), cy(imx), czero;


      czero=cmplx(0.0,0.0);
      third=1./3.;

      for(auto i=0; i<im; ++i) {
        q(i)=-p2(i)*p2(i)/36.+(p3(i)*p1(i)-4*p0(i))/12.;
        r(i)=-(p2(i)/6)**3+p2(i)*(p3(i)*p1(i)-4*p0(i))/48.
            +(4*p0(i)*p2(i)-p0(i)*p3(i)*p3(i)-p1(i)*p1(i))/16;

        crad(i)=r(i)*r(i)+q(i)*q(i)*q(i);
        crad(i)=sqrt(crad(i));

        cb(i)=r(i)-crad(i);
        if(cb(i) == czero) {
           //insoluble particle
           cx(1,i)=(-p1(i))**third;
           cx(2,i)=cx(1,i);
           cx(3,i)=cx(1,i);
           cx(4,i)=cx(1,i);
        }

        cb(i)=cb(i)**third;

        cy(i)=-cb(i)+q(i)/cb(i)+p2(i)/6;

        cb0(i)=sqrt(cy(i)*cy(i)-p0(i));
        cb1(i)=(p3(i)*cy(i)-p1(i))/(2*cb0(i));

         cb(i)=p3(i)/2+cb1(i);
         crad(i)=cb(i)*cb(i)-4*(cy(i)+cb0(i));
         crad(i)=sqrt(crad(i));
         cx(1,i)=(-cb(i)+crad(i))/2.;
         cx(2,i)=(-cb(i)-crad(i))/2.;

         cb(i)=p3(i)/2-cb1(i);
         crad(i)=cb(i)*cb(i)-4*(cy(i)-cb0(i));
         crad(i)=sqrt(crad(i));
         cx(3,i)=(-cb(i)+crad(i))/2.;
         cx(4,i)=(-cb(i)-crad(i))/2.;
      }
   10 continue

      return;
    }
}
#endif

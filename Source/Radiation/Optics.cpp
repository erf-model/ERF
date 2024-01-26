//
// get cloud and aerosol optics for short wave and long wave
//
#include "Optics.H"
#include "Slingo.H"
#include "Ebert_curry.H"
#include "Rad_constants.H"
#include "AMReX_Random.H"
using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

void Optics::initialize(int ngas, int nmodes, int num_aeros,
                        int nswbands, int nlwbands,
                        int ncol, int nlev, int nrh, int top_lev,
                        const std::vector<std::string>& aero_names,
                        const real2d& zi, const real2d& pmid, const real2d& temp,
                        const real2d& qi, const real2d& geom_radius) {
   cloud_optics.initialize();
   aero_optics.initialize(ngas, nmodes, num_aeros,
                          nswbands, nlwbands, ncol, nlev, nrh, top_lev,
                          aero_names, zi, pmid, temp, qi, geom_radius);
}

void Optics::finalize() {
   // nothing needs to be done here
}


void Optics::get_cloud_optics_sw(int ncol, int nlev, int nbnd,
         bool do_snow, const real2d& cld, const real2d& cldfsnow, const real2d& iclwp,
         const real2d& iciwp, const real2d& icswp, const real2d& lambdac, const real2d& mu,
         const real2d& dei, const real2d& des, const real2d& rel, const real2d& rei,
         const real3d& tau_out, const real3d& ssa_out, const real3d& asm_out,
         const real3d& liq_tau_out, const real3d& ice_tau_out, const real3d& snw_tau_out) {
      //Temporary variables to hold cloud optical properties before combining into
      // output arrays. Same shape as output arrays, so get shapes from output.
      real3d liq_tau("liq_tau",nbnd, ncol, nlev);
      real3d liq_tau_ssa("liq_tau_ssa",nbnd, ncol, nlev);
      real3d liq_tau_ssa_g("liq_tau_ssa_g", nbnd, ncol, nlev);
      real3d liq_tau_ssa_f("liq_tau_ssa_f", nbnd, ncol, nlev);
      real3d ice_tau("ice_tau", nbnd, ncol, nlev);
      real3d ice_tau_ssa("ice_tau_ssa", nbnd, ncol, nlev);
      real3d ice_tau_ssa_g("ice_tau_ssa_g", nbnd, ncol, nlev);
      real3d ice_tau_ssa_f("ice_tau_ssa_f", nbnd, ncol, nlev);
      real3d cld_tau("cld_tau", nbnd, ncol, nlev);
      real3d cld_tau_ssa("cld_tau_ssa", nbnd, ncol, nlev);
      real3d cld_tau_ssa_g("cld_tau_ssa_g", nbnd, ncol, nlev);
      real3d cld_tau_ssa_f("cld_tau_ssa_f", nbnd, ncol, nlev);
      real3d snow_tau("snow_tau", nbnd, ncol, nlev );
      real3d snow_tau_ssa("snow_tau_ssa", nbnd, ncol, nlev);
      real3d snow_tau_ssa_g("snow_tau_ssa_g", nbnd, ncol, nlev);
      real3d snow_tau_ssa_f("snow_tau_ssa_f", nbnd, ncol, nlev);
      real3d combined_tau("combined_tau", nbnd, ncol, nlev);
      real3d combined_tau_ssa("combined_tau_ssa", nbnd, ncol, nlev);
      real3d combined_tau_ssa_g("combined_tau_ssa_g", nbnd, ncol, nlev);
      real3d combined_tau_ssa_f("combined_tau_ssa_f", nbnd, ncol, nlev);

      // Initialize outputs
      yakl::memset(tau_out, 0.);
      yakl::memset(ssa_out, 0.);
      yakl::memset(asm_out, 0.);
      yakl::memset(liq_tau_out, 0.);
      yakl::memset(ice_tau_out, 0.);
      yakl::memset(snw_tau_out, 0.);

      // Initialize local variables
      yakl::memset(ice_tau, 0.);
      yakl::memset(ice_tau_ssa, 0.);
      yakl::memset(ice_tau_ssa_g, 0.);
      yakl::memset(ice_tau_ssa_f, 0.);
      yakl::memset(liq_tau, 0.);
      yakl::memset(liq_tau_ssa, 0.);
      yakl::memset(liq_tau_ssa_g, 0.);
      yakl::memset(liq_tau_ssa_f, 0.);
      yakl::memset(snow_tau, 0.);
      yakl::memset(snow_tau_ssa, 0.);
      yakl::memset(snow_tau_ssa_g, 0.);
      yakl::memset(snow_tau_ssa_f, 0.);
      yakl::memset(combined_tau, 0.);
      yakl::memset(combined_tau_ssa, 0.);
      yakl::memset(combined_tau_ssa_g, 0.);
      yakl::memset(combined_tau_ssa_f, 0.);

      // Get ice cloud optics
      if (icecldoptics == "mitchell") {
         cloud_optics.mitchell_ice_optics_sw(
            ncol, nlev, iciwp, dei,
            ice_tau, ice_tau_ssa,
            ice_tau_ssa_g, ice_tau_ssa_f);
      }
      else if (icecldoptics == "ebertcurry") {
         EbertCurry::ec_ice_optics_sw(ncol, nlev, nbnd, cld, iciwp, 
                                     rei, ice_tau, ice_tau_ssa, 
                                     ice_tau_ssa_g, ice_tau_ssa_f);
      }

      // Get liquid cloud optics
      if (liqcldoptics == "gammadist") {
         cloud_optics.gammadist_liq_optics_sw(
            ncol, nlev, iclwp, lambdac, mu,
            liq_tau, liq_tau_ssa,
            liq_tau_ssa_g, liq_tau_ssa_f);
      }
      else if (liqcldoptics == "slingo") {
         Slingo::slingo_liq_optics_sw(
            ncol, nlev, nbnd, cld, iclwp, rel,
            liq_tau, liq_tau_ssa,
            liq_tau_ssa_g, liq_tau_ssa_f);
      }

      // Get snow cloud optics
      if (do_snow) {
         cloud_optics.mitchell_ice_optics_sw(
            ncol, nlev, icswp, des,
            snow_tau, snow_tau_ssa,
            snow_tau_ssa_g, snow_tau_ssa_f);
      }
      else {
         // We are not doing snow optics, so set these to zero so we can still use
         // the arrays without additional logic
         yakl::memset(snow_tau, 0.);
         yakl::memset(snow_tau_ssa, 0.);
         yakl::memset(snow_tau_ssa_g, 0.);
         yakl::memset(snow_tau_ssa_f, 0.);
      }

      // Combine all cloud optics from CAM routines
      parallel_for(SimpleBounds<3>(nbnd, ncol, nlev), YAKL_LAMBDA (int ibnd, int icol, int ilev) {
         cld_tau(ibnd, icol, ilev) = ice_tau(ibnd, icol, ilev) + liq_tau(ibnd, icol, ilev);
         cld_tau_ssa(ibnd, icol, ilev) = ice_tau_ssa(ibnd, icol, ilev) + liq_tau_ssa(ibnd, icol, ilev);
         cld_tau_ssa_g(ibnd, icol, ilev) = ice_tau_ssa_g(ibnd, icol, ilev) + liq_tau_ssa_g(ibnd, icol, ilev);
      });

      if (do_snow) {
         combine_properties( nbnd, ncol, nlev,
            cld, cld_tau, cldfsnow, snow_tau, combined_tau);

         combine_properties( nbnd, ncol, nlev,
            cld, cld_tau_ssa, cldfsnow, snow_tau_ssa, combined_tau_ssa);

         combine_properties( nbnd, ncol, nlev,
            cld, cld_tau_ssa_g, cldfsnow, snow_tau_ssa_g, combined_tau_ssa_g);
      }
      else {
        parallel_for(SimpleBounds<3>(nbnd, ncol, nlev), YAKL_LAMBDA (int ibnd, int icol, int ilev) {
           combined_tau(ibnd, icol, ilev) = cld_tau(ibnd, icol, ilev);
           combined_tau_ssa(ibnd, icol, ilev) = cld_tau_ssa(ibnd, icol, ilev);
           combined_tau_ssa_g(ibnd, icol, ilev) = cld_tau_ssa_g(ibnd, icol, ilev);
        });
      }

      // Copy to output arrays, converting to optical depth, single scattering
      // albedo, and assymmetry parameter from the products that the CAM routines
      // return. Make sure we do not try to divide by zero...
      parallel_for(SimpleBounds<3>(nbnd, ncol, nlev), YAKL_LAMBDA (int iband, int icol, int ilev) {
          tau_out(icol,ilev,iband) = combined_tau(iband,icol,ilev);
          if (combined_tau(iband,icol,ilev) > 0) {
             ssa_out(icol,ilev,iband)
                 = combined_tau_ssa(iband,icol,ilev) / combined_tau(iband,icol,ilev);
          } else {
             ssa_out(icol,ilev,iband) = 1.;
          }

          if (combined_tau_ssa(iband,icol,ilev) > 0) {
             asm_out(icol,ilev,iband)
                 = combined_tau_ssa_g(iband,icol,ilev) / combined_tau_ssa(iband,icol,ilev);
          } else {
             asm_out(icol,ilev,iband) = 0.;
          }

          // Re-order diagnostics outputs
          liq_tau_out(icol,ilev,iband) = liq_tau(iband,icol,ilev);
          ice_tau_out(icol,ilev,iband) = ice_tau(iband,icol,ilev);
          snw_tau_out(icol,ilev,iband) = snow_tau(iband,icol,ilev);
     });
 }

   //----------------------------------------------------------------------------
void Optics::get_cloud_optics_lw(
         int ncol, int nlev, int nbnd, bool do_snow, const real2d& cld, const real2d& cldfsnow, const real2d& iclwp,
         const real2d& iciwp, const real2d& icswp, const real2d& lambdac, const real2d& mu, const real2d& dei, const real2d& des,
         const real2d& rei, const real3d& tau_out, const real3d& liq_tau_out, const real3d& ice_tau_out, const real3d& snw_tau_out) {
      // Temporary variables to hold absorption optical depth
      real3d ice_tau("ice_tau", nbnd, ncol, nlev);
      real3d liq_tau("liq_tau", nbnd, ncol, nlev);
      real3d snow_tau("snow_tau", nbnd, ncol, nlev);
      real3d cld_tau("cld_tau", nbnd, ncol, nlev);
      real3d combined_tau("combined_tau", nbnd, ncol, nlev);

      // Initialize outputs
      yakl::memset(tau_out, 0.);
      yakl::memset(liq_tau_out, 0.);
      yakl::memset(ice_tau_out, 0.);
      yakl::memset(snw_tau_out, 0.);

      // Initialize local variables
      yakl::memset(ice_tau, 0.);
      yakl::memset(liq_tau, 0.);
      yakl::memset(snow_tau, 0.);
      yakl::memset(cld_tau, 0.);
      yakl::memset(combined_tau, 0.);

      // Get ice optics
      if (icecldoptics == "mitchell") {
         cloud_optics.mitchell_ice_optics_lw(ncol, nlev, iciwp, dei, ice_tau);
      }
      else if (icecldoptics == "ebertcurry") {
         EbertCurry::ec_ice_optics_lw(ncol, nlev, nbnd, cld, iclwp, iciwp, rei, ice_tau);
      }

      // Get liquid optics
      if (liqcldoptics == "gammadist") {
         cloud_optics.gammadist_liq_optics_lw(ncol, nlev, iclwp, lambdac, mu, liq_tau);
      }
      else if (liqcldoptics == "slingo") {
         Slingo::slingo_liq_optics_lw(ncol, nlev, nbnd, cld, iclwp, iciwp, liq_tau);
      }

      // Combined cloud optics
      parallel_for(SimpleBounds<3>(nbnd, ncol, nlev), YAKL_LAMBDA (int ibnd, int icol, int ilev) {
          cld_tau(ibnd, icol, ilev) = liq_tau(ibnd, icol, ilev) + ice_tau(ibnd, icol, ilev);
       });

      // Get snow optics?
      if (do_snow) {
         cloud_optics.mitchell_ice_optics_lw(ncol, nlev, icswp, des, snow_tau);
         combine_properties(nbnd, ncol, nlev,
                            cld, cld_tau, cldfsnow, snow_tau, combined_tau);
      }
      else {
         parallel_for(SimpleBounds<3>(nbnd, ncol, nlev), YAKL_LAMBDA (int ibnd, int icol, int ilev) {
            combined_tau(ibnd,icol,ilev) = cld_tau(ibnd,icol,ilev);
         });
      }

      // Set output optics
      parallel_for(SimpleBounds<3>(nbnd, ncol, nlev), YAKL_LAMBDA (int ibnd, int icol, int ilev) {
         tau_out(icol,ilev,ibnd) = combined_tau(ibnd,icol,ilev);
         liq_tau_out(icol,ilev,ibnd) = liq_tau(ibnd,icol,ilev);
         ice_tau_out(icol,ilev,ibnd) = ice_tau(ibnd,icol,ilev);
         snw_tau_out(icol,ilev,ibnd) = snow_tau(ibnd,icol,ilev);
     });
  }

//----------------------------------------------------------------------------
// Provide a procedure to combine cloud optical properties by weighting
// contributions by fraction present. I.e., for combining cloud and snow
// optical properties, we weight the cloud and snow properties by the fraction
// of cloud and snow present.
void Optics::combine_properties(int nbands, int ncols, int nlevs,
                                const real2d& fraction1, const real3d& property1,
                                const real2d& fraction2, const real3d& property2,
                                const real3d& combined_property) {
      // Combined fraction (max of property 1 and 2)
      real2d combined_fraction("combined_fraction", ncols,nlevs);

      // Combined fraction
      parallel_for(SimpleBounds<2>(ncols, nlevs), YAKL_LAMBDA (int icol, int ilev) {
         combined_fraction(icol, ilev) = std::max(fraction1(icol, ilev), fraction2(icol, ilev));
      });

      // Combine optical properties by weighting by amount of cloud and snow
      yakl::memset(combined_property, 0.);
      parallel_for(SimpleBounds<3>(nlevs, ncols, nbands), YAKL_LAMBDA (int ilev, int icol, int iband) {
         if (combined_fraction(icol,ilev) > 0) {
              combined_property(iband,icol,ilev) = (
                fraction1(icol,ilev) * property1(iband,icol,ilev)
                + fraction2(icol,ilev) * property2(iband,icol,ilev)
              ) / combined_fraction(icol,ilev);
         }
     });
}

//----------------------------------------------------------------------------
// Do MCICA sampling of optics here. This will map bands to gpoints,
// while doing stochastic sampling of cloud state
void Optics::sample_cloud_optics_sw(
         int ncol, int nlev, int ngpt, const int1d& gpt2bnd,
         const real2d& pmid, const real2d& cld, const real2d& cldfsnow,
         const real3d& tau_bnd, const real3d& ssa_bnd, const real3d& asm_bnd,
         const real3d& tau_gpt, const real3d& ssa_gpt, const real3d& asm_gpt) {

      //real(r8), dimension(ncol,nlev) :: combined_cld
      real2d combined_cld("combined_cld", ncol, nlev);

      //logical, dimension(ngpt,ncol,nlev) :: iscloudy
      bool3d iscloudy("iscloudy", ngpt, ncol, nlev);

      // Combined snow and cloud fraction
      parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev) {
          combined_cld(icol,ilev) = std::max(cld(icol,ilev), cldfsnow(icol,ilev));
      });

      // Get stochastic subcolumn cloud mask
      mcica_subcol_mask(ngpt, ncol, nlev, combined_cld, iscloudy);

      // Generate subcolumns for homogeneous clouds
      parallel_for(SimpleBounds<3>(ngpt, nlev, ncol), YAKL_LAMBDA (int igpt, int ilev, int icol) {
        if (iscloudy(igpt,icol,ilev) && combined_cld(icol,ilev) > 0.) {
           tau_gpt(icol,ilev,igpt) = tau_bnd(icol,ilev,gpt2bnd(igpt));
           ssa_gpt(icol,ilev,igpt) = ssa_bnd(icol,ilev,gpt2bnd(igpt));
           asm_gpt(icol,ilev,igpt) = asm_bnd(icol,ilev,gpt2bnd(igpt));
        } else {
           tau_gpt(icol,ilev,igpt) = 0.;
           ssa_gpt(icol,ilev,igpt) = 1.;
           asm_gpt(icol,ilev,igpt) = 0.;
        }
     });
 }

//----------------------------------------------------------------------------
// Do MCICA sampling of optics here. This will map bands to gpoints,
// while doing stochastic sampling of cloud state
void Optics::sample_cloud_optics_lw(int ncol, int nlev, int ngpt, const int1d& gpt2bnd,
                                    const real2d& pmid, const real2d& cld, const real2d& cldfsnow,
                                    const real3d& tau_bnd, const real3d& tau_gpt) {
      //real(r8), dimension(ncol,nlev) :: combined_cld
      real2d combined_cld("combined_cld", ncol, nlev);

      //logical, dimension(ngpt,ncol,nlev) :: iscloudy
      bool3d iscloudy("iscloudy", ngpt, ncol, nlev);

      // Combine cloud and snow fractions for MCICA sampling
      parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev) {
          combined_cld(icol,ilev) = std::max(cld(icol,ilev), cldfsnow(icol,ilev));
      });

      // Get the stochastic subcolumn cloudy mask
      mcica_subcol_mask(ngpt, ncol, nlev, combined_cld, iscloudy);

      // Map optics to g-points, selecting a single subcolumn for each
      // g-point. This implementation generates homogeneous clouds, but it would be
      // straightforward to extend this to handle horizontally heterogeneous clouds
      // as well.
      parallel_for(SimpleBounds<3>(ngpt, nlev, ncol), YAKL_LAMBDA (int igpt, int ilev, int icol) {
         if (iscloudy(igpt,icol,ilev) && combined_cld(icol,ilev) > 0.) {
             tau_gpt(icol,ilev,igpt) = tau_bnd(icol,ilev,gpt2bnd(igpt));
         } else {
            tau_gpt(icol,ilev,igpt) = 0.;
         }
      });
 }

//----------------------------------------------------------------------------
void Optics::set_aerosol_optics_sw(int icall, int ncol, int nlev, int nswbands, real dt,
                                   const int1d& night_indices, bool is_cmip6_volc, const real3d& tau_out,
                                   const real3d& ssa_out, const real3d& asm_out, const real2d& clear_rh) {
      // NOTE: aer_rad_props expects 0:pver indexing on these! It appears this is to
      // account for the extra layer added above model top, but it is not entirely
      // clear. This is not done for the longwave, and it is not really documented
      // anywhere that I can find. Regardless, optical properties for the zero index
      // are set to zero in aer_rad_props_sw as far as I can tell.
      //
      // NOTE: dimension ordering is different than for cloud optics!
      real3d tau("tau", ncol, nlev+1, nswbands);
      real3d tau_w("tau_w", ncol, nlev+1, nswbands);
      real3d tau_w_g("tau_w_g", ncol, nlev+1, nswbands);
      real3d tau_w_f("tau_w_f", ncol, nlev+1, nswbands);

      // Get aerosol absorption optical depth from CAM routine
      yakl::memset(tau, 0.);
      yakl::memset(tau_w, 0.);
      yakl::memset(tau_w_g, 0.);
      yakl::memset(tau_w_f, 0.);

      int1d ic("icount",1);
      intHost1d ic_host("ic_host",1);
      parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA (int i) {
        if (night_indices(i) > 0) ++ic(1);
      });
      ic.deep_copy_to(ic_host);

      aero_optics.aer_rad_props_sw(icall, dt,
           ic_host(1), night_indices, is_cmip6_volc,
           tau, tau_w, tau_w_g, tau_w_f, clear_rh);

      // Extract quantities from products
      parallel_for(SimpleBounds<3>(ncol, nlev, nswbands), YAKL_LAMBDA (int icol, int ilev, int iband) {
        // Copy cloud optical depth over directly
        tau_out(icol,ilev,iband) = tau(icol,ilev,iband);
        // Extract single scattering albedo from the product-defined fields
        if (tau(icol,ilev,iband) > 0) {
          ssa_out(icol,ilev,iband) = tau_w(icol,ilev,iband) / tau(icol,ilev,iband);
        } else {
          ssa_out(icol,ilev,iband) = 1.;
        }

        // Extract assymmetry parameter from the product-defined fields
        if (tau_w(icol,ilev,iband) > 0) {
          asm_out(icol,ilev,iband) = tau_w_g(icol,ilev,iband) / tau_w(icol,ilev,iband);
        } else {
           asm_out(icol,ilev,iband) = 0.;
        }
     });
  }

//----------------------------------------------------------------------------
void Optics::set_aerosol_optics_lw(int icall, real dt, bool is_cmip6_volc, const real2d& zi,
             const real3d& tau, const real2d& clear_rh) {
   // Get aerosol absorption optical depth from CAM routine
   yakl::memset(tau, 0.);
   aero_optics.aer_rad_props_lw(is_cmip6_volc, icall, dt, zi, tau, clear_rh);
}


// Arrays use CAM vertical index convention: index increases from top to bottom.
// This index ordering is assumed in the maximum-random overlap algorithm which starts
// at the top of a column and marches down, with each layer depending on the state
// of the layer above it.
//
void Optics::mcica_subcol_mask(int ngpt, int ncol, int nlev,
                               const real2d& cldfrac, const bool3d& iscloudy) {
   // Local vars
   const real cldmin = 1.0e-80;      // min cloud fraction
   real2d cldf("cldf",ncol,nlev);    // cloud fraction clipped to cldmin
   real3d cdf("cdf", ngpt, ncol, nlev);   // random numbers

   // clip cloud fraction
   parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev) {
      cldf(icol,ilev) = cldfrac(icol,ilev);
      if (cldf(icol,ilev) < cldmin) cldf(icol,ilev) = 0.;
   });

   amrex::RandomEngine engine;
   // Generate random numbers in each subcolumn at every level
    parallel_for(SimpleBounds<3>(ngpt, ncol, nlev), YAKL_LAMBDA (int isubcol, int icol, int ilev) {
      cdf(isubcol,icol,ilev) = amrex::Random(engine);
   });

   // Maximum-Random overlap
   // i) pick a random number for top layer.
   // ii) walk down the column:
   //    - if the layer above is cloudy, use the same random number as in the layer above
   //    - if the layer above is clear, use a new random number
   parallel_for(SimpleBounds<3>(nlev, ncol, ngpt), YAKL_LAMBDA (int k, int i, int isubcol) {
     if (k > 1) {
        if (cdf(isubcol,i,k-1) > 1. - cldf(i,k-1) ) {
           cdf(isubcol,i,k) = cdf(isubcol,i,k-1);
        } else {
           cdf(isubcol,i,k) = cdf(isubcol,i,k) * (1. - cldf(i,k-1));
        }
        iscloudy(isubcol,i,k) = cdf(isubcol,i,k) >= 1.-cldf(i,k) ? true : false;
     }
   });
}


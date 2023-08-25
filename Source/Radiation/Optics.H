/*
 * Optics class stores shortwave and longwave cloud optical properties *by band*.
 * Dimensions are nswbands,ncol,nlev. Generally, this should be able to handle cases
 * were ncol might be something like nday, and nlev could be arbitrary so long as
 * corresponding fields were defined for all indices of nlev.
 */
#ifndef ERF_OPTICS_H
#define ERF_OPTICS_H

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

// Radiation code interface class
class Optics {
  public:
    // constructor
    Optics() = default;
    explicit Optics(int ngases, char* gas_names[]);
    // deconstructor
    ~Optics() = default;

    // get short wave cloud optics property
    void get_cloud_optics_sw(int ncol, int nlev, int nbnd,
         bool do_snow, real2d cld, real2d cldfsnow, real2d iclwp,
         real2d iciwp, real2d icswp, real2d lambdac, real2d mu,
         real2d dei, real2d des, real2d rel, real2d rei,
         real3d tau_out, real3d ssa_out, real3d asm_out,
         real3d liq_tau_out, real3d ice_tau_out, real3d snw_tau_out);


    // get long wave cloud optics property
    void get_cloud_optics_lw(int ncol, int nlev, int nbnd, bool do_snow,
          real2d& cld, real2d& cldfsnow, real2d& iclwp, real2d& iciwp, real2d& icswp,
          real2d& lambdac, real2d& mu, real2d& dei, real2d& des, real2d& rei,
          real3d& tau_out, real3d& liq_tau_out, real3d& ice_tau_out, real3d& snw_tau_out);

    // sample short wave cloud optics property
    void sample_cloud_optics_sw(int ncol, int nlev, int ngpt, int1d& gpt2bnd,
         real2d& pmid, real2d& cld, real2d& cldfsnow,
         real3d& tau_bnd, real3d& ssa_bnd, real3d& asm_bnd,
         real3d& tau_gpt, real3d& ssa_gpt, real3d& asm_gpt);

    // sample long wave cloud optics property
    void sample_cloud_optics_lw(
         int ncol, int nlev, int ngpt, int1d& gpt2bnd,
         real2d& pmid, real2d& cld, real2d& cldfsnow,
         real3d& tau_bnd, real3d& tau_gpt);

    // set the short wave aerosol optics property
    void set_aerosol_optics_sw(int icall, int ncol, int nlev, int nswbands, real dt, int1d& night_indices,
                               bool is_cmip6_volc,
                               real3d& tau_out, real3d& ssa_out, real3d& asm_out,
                               real2d& clear_rh);

    // set the long wave aerosol optics property
    void set_aerosol_optics_lw(int icall, real dt, bool is_cmip6_volc, real2d& zi, real3d& tau, real2d& clear_rh);

    // mcica subcol mask
    void mcica_subcol_mask(int ngpt, int ncol, int nlev, real2d& cldfrac, bool3d& iscloudy);

    // combine properties
    void combine_properties(int nbands, int ncols, int nlevs,
                            real2d& fraction1, real3d& property1,
                            real2d& fraction2, real3d& property2,
                            real3d& combined_property);

    // initialize and load gas property data for rrtmgp radiation
    void initialize();

    // finalize/clean up
    void finalize();

  private:
   // number of gas for radiation model
   int ngas;
   char** gas_names;

   std::string icecldoptics;
   std::string liqcldoptics;

   CloudRadProps cloud_optics;
   AerRadProps  aero_optics;
};

#endif // ERF_OPTICS_H

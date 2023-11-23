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
         bool do_snow, const real2d& cld, const real2d& cldfsnow, const real2d& iclwp,
         const real2d& iciwp, const real2d& icswp, const real2d& lambdac, const real2d& mu,
         const real2d& dei, const real2d& des, const real2d& rel, const real2d& rei,
         const real3d& tau_out, const real3d& ssa_out, const real3d& asm_out,
         const real3d& liq_tau_out, const real3d& ice_tau_out, const real3d& snw_tau_out);

    // get long wave cloud optics property
    void get_cloud_optics_lw(int ncol, int nlev, int nbnd, bool do_snow, const real2d& cld, const real2d& cldfsnow, const real2d& iclwp, 
         const real2d& iciwp, const real2d& icswp, const real2d& lambdac, const real2d& mu, const real2d& dei, const real2d& des, 
         const real2d& rei, const real3d& tau_out, const real3d& liq_tau_out, const real3d& ice_tau_out, const real3d& snw_tau_out);

    // sample short wave cloud optics property
    void sample_cloud_optics_sw(int ncol, int nlev, int ngpt, const int1d& gpt2bnd,
         const real2d& pmid, const real2d& cld, const real2d& cldfsnow,
         const real3d& tau_bnd, const real3d& ssa_bnd, const real3d& asm_bnd,
         const real3d& tau_gpt, const real3d& ssa_gpt, const real3d& asm_gpt);

    // sample long wave cloud optics property
    void sample_cloud_optics_lw(int ncol, int nlev, int ngpt, const int1d& gpt2bnd,
         const real2d& pmid, const real2d& cld, const real2d& cldfsnow,
         const real3d& tau_bnd, const real3d& tau_gpt);

    // set the short wave aerosol optics property
    void set_aerosol_optics_sw(int icall, int ncol, int nlev, int nswbands, real dt, const int1d& night_indices,
             bool is_cmip6_volc, const real3d& tau_out, const real3d& ssa_out, const real3d& asm_out, 
             const real2d& clear_rh);

    // set the long wave aerosol optics property
    void set_aerosol_optics_lw(int icall, real dt, bool is_cmip6_volc, const real2d& zi, 
             const real3d& tau, const real2d& clear_rh);

    // mcica subcol mask
    void mcica_subcol_mask(int ngpt, int ncol, int nlev, const real2d& cldfrac, const bool3d& iscloudy);

    // combine properties
    void combine_properties(int nbands, int ncols, int nlevs,
                            const real2d& fraction1, const real3d& property1,
                            const real2d& fraction2, const real3d& property2,
                            const real3d& combined_property);

    // initialize and load gas property data for rrtmgp radiation
    void initialize(int ngas, int nmodes, int num_aeros,
                    int nswbands, int nlwbands,
                    int ncol, int nlev, int nrh, int top_lev,
                    const std::vector<std::string>& aero_names,
                    const real2d& zi, const real2d& pmid, const real2d& temp,
                    const real2d& qi, const real2d& geom_radius);
   
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

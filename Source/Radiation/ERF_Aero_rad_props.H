//
// cloud radiation properties class
//
#ifndef ERF_AER_RAD_PROPS_H_
#define ERF_AER_RAD_PROPS_H_

#include <AMReX_MultiFabUtil.H>
#include "ERF_Mam4_aero.H"
#include "YAKL_netcdf.h"
#include "rrtmgp_const.h"
#include "ERF_Rad_constants.H"
#include "ERF_Mam4_constituents.H"

//
// AerRadProps class converts the aerosol mass to bulk optical properties for short wave
// and long wave radiation computations
//
class AerRadProps {
 public:
   AerRadProps () = default;
   ~AerRadProps () = default;

   // initialization
   void initialize (int num_gas, int num_modes, int naeroes,
                    int nswbands_, int nlwbands_,
                    int ncoloum, int nlevel, int num_rh, int top_levels,
                    const std::vector<std::string>& aerosol_names,
                    const real2d& zint, const real2d& pmiddle, const real2d& pdel,
                    const real2d& temperature, const real2d& qtotal,
                    const real2d& geom_rad);

   void aer_rad_props_sw (const int& list_idx,
                          const real& dt,
                          const int& nnite,
                          const int1d& idxnite,
                          const bool is_cmip6_volc,
                          const real3d& tau,
                          const real3d& tau_w,
                          const real3d& tau_w_g,
                          const real3d& tau_w_f,
                          const real2d& clear_rh);

   void aer_rad_props_lw (const bool& is_cmip6_volc,
                          const int& list_idx,
                          const real& dt,
                          const real2d& zi,
                          const real3d& odap_aer,
                          const real2d& clear_rh);


   void get_hygro_rad_props (const int& ncol,
                             const int2d& krh,
                             const real2d& wrh,
                             const real2d& mass,
                             const real2d& ext,
                             const real2d& ssa,
                             const real2d& assm,
                             const real3d& tau,
                             const real3d& tau_w,
                             const real3d& tau_w_g,
                             const real3d& tau_w_f);

  void get_nonhygro_rad_props (const int& ncol,
                               const real2d& mass,
                               const real1d& ext,
                               const real1d& ssa,
                               const real1d& assm,
                               const real3d& tau,
                               const real3d& tau_w,
                               const real3d& tau_w_g,
                               const real3d& tau_w_f);

 void get_volcanic_radius_rad_props (const int& ncol,
                                     const real2d& mass,
                                     const real2d& r_ext,
                                     const real2d& r_scat,
                                     const real2d& r_ascat,
                                     const real1d& r_mu,
                                     const real3d& tau,
                                     const real3d& tau_w,
                                     const real3d& tau_w_g,
                                     const real3d& tau_w_f);

 void volcanic_cmip_sw (const int1d& trop_level,
                        const real2d& zi,           // zone interface
                        const real3d& ext_cmip6_sw_inv_m,
                        const real3d& ssa_cmip6_sw, // ssa from the volcanic inputs
                        const real3d& af_cmip6_sw,  // asymmetry factor (af) from volcanic inputs
                        const real3d& tau,
                        const real3d& tau_w,
                        const real3d& tau_w_g,
                        const real3d& tau_w_f);

 void get_volcanic_rad_props (const int& ncol,
                              const real2d& mass,
                              const real1d& ext,
                              const real1d& scat,
                              const real1d& ascat,
                              const real3d& tau,
                              const real3d& tau_w,
                              const real3d& tau_w_g,
                              const real3d& tau_w_f);

 void aer_vis_diag_out (const int& ncol,
                        const int& nnite,
                        const int1d& idxnite,
                        const int& iaer,
                        const real2d& tau,
                        const int& diag_idx);

 private:
   int nmodes, ngas, num_aeroes;
   std::vector<std::string> aero_names;

   int nswbands, nlwbands;
   int ncol, nlev, nrh, top_lev;

   // aero constituents
   MamConstituents mam_consti;

   // aero
   Mam4_aer mam_aer;

   // pdel dry
   real2d pdeldry;
   // pmid
   real2d pmid, temp, qt;
   // geometric radius
   real2d geometric_radius;
   // vertical grid
   real2d zi;
   // cmip6 short wave
   real3d ext_cmip6_sw, ssa_cmip6_sw, af_cmip6_sw;
};
#endif

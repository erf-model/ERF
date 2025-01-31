#ifndef ERF_RRTMGP_INTERFACE_H
#define ERF_RRTMGP_INTERFACE_H

#include <rrtmgp_const.h>
#include <mo_gas_optics_rrtmgp.h>
#include <mo_cloud_optics.h>
#include <mo_fluxes_byband.h>
#include <mo_load_coefficients.h>
#include <mo_gas_concentrations.h>
#include <mo_gas_optics_rrtmgp.h>
#include <mo_cloud_optics.h>
#include <mo_rte_sw.h>
#include <mo_rte_lw.h>
#include <mo_load_cloud_coefficients.h>

#include <ERF_RRTMGP_Utils.H>
#include <ERF_Constants.H>

void init_kls ();
void finalize_kls ();

namespace rrtmgp {

extern GasOpticsRRTMGP k_dist_sw;
extern GasOpticsRRTMGP k_dist_lw;

extern CloudOptics cloud_optics_sw;
extern CloudOptics cloud_optics_lw;

extern bool initialized;

void
rrtmgp_initialize (GasConcs &gas_concs,
                   const std::string& coefficients_file_sw, const std::string& coefficients_file_lw,
                   const std::string& cloud_optics_file_sw, const std::string& cloud_optics_file_lw);


void
rrtmgp_finalize ();


void
compute_band_by_band_surface_albedos (const int ncol, const int nswbands,
                                      real1d &sfc_alb_dir_vis, real1d &sfc_alb_dir_nir,
                                      real1d &sfc_alb_dif_vis, real1d &sfc_alb_dif_nir,
                                      real2d &sfc_alb_dir,     real2d &sfc_alb_dif);


void
compute_broadband_surface_fluxes (const int ncol, const int ktop, const int nswbands,
                                  real3d &sw_bnd_flux_dir , real3d &sw_bnd_flux_dif ,
                                  real1d &sfc_flux_dir_vis, real1d &sfc_flux_dir_nir,
                                  real1d &sfc_flux_dif_vis, real1d &sfc_flux_dif_nir);


void
rrtmgp_main (const int ncol, const int nlay,
             real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
             GasConcs &gas_concs,
             real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0,
             real2d &lwp, real2d &iwp, real2d &rel, real2d &rei, real2d &cldfrac,
             real3d &aer_tau_sw, real3d &aer_ssa_sw, real3d &aer_asm_sw, real3d &aer_tau_lw,
             real3d &cld_tau_sw_bnd, real3d &cld_tau_lw_bnd,
             real3d &cld_tau_sw_gpt, real3d &cld_tau_lw_gpt,
             real2d &sw_flux_up, real2d &sw_flux_dn, real2d &sw_flux_dn_dir,
             real2d &lw_flux_up, real2d &lw_flux_dn,
             real2d &sw_clnclrsky_flux_up, real2d &sw_clnclrsky_flux_dn, real2d &sw_clnclrsky_flux_dn_dir,
             real2d &sw_clrsky_flux_up, real2d &sw_clrsky_flux_dn, real2d &sw_clrsky_flux_dn_dir,
             real2d &sw_clnsky_flux_up, real2d &sw_clnsky_flux_dn, real2d &sw_clnsky_flux_dn_dir,
             real2d &lw_clnclrsky_flux_up, real2d &lw_clnclrsky_flux_dn,
             real2d &lw_clrsky_flux_up, real2d &lw_clrsky_flux_dn,
             real2d &lw_clnsky_flux_up, real2d &lw_clnsky_flux_dn,
             real3d &sw_bnd_flux_up, real3d &sw_bnd_flux_dn, real3d &sw_bnd_flux_dn_dir,
             real3d &lw_bnd_flux_up, real3d &lw_bnd_flux_dn,
             const real tsi_scaling,
             const bool extra_clnclrsky_diag = false, const bool extra_clnsky_diag = false);


void
rrtmgp_sw (const int ncol, const int nlay,
           GasOpticsRRTMGP &k_dist,
           real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
           GasConcs &gas_concs,
           real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0,
           OpticalProps2str &aerosol, OpticalProps2str &clouds,
           FluxesByband &fluxes,
           FluxesBroadband &clnclrsky_fluxes,
           FluxesBroadband &clrsky_fluxes,
           FluxesBroadband &clnsky_fluxes,
           const real tsi_scaling,
           const bool extra_clnclrsky_diag,
           const bool extra_clnsky_diag);


void
rrtmgp_lw (const int ncol, const int nlay,
           GasOpticsRRTMGP &k_dist,
           real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
           GasConcs &gas_concs,
           OpticalProps1scl &aerosol, OpticalProps1scl &clouds,
           FluxesByband &fluxes,
           FluxesBroadband &clnclrsky_fluxes,
           FluxesBroadband &clrsky_fluxes,
           FluxesBroadband &clnsky_fluxes,
           const bool extra_clnclrsky_diag,
           const bool extra_clnsky_diag);


int3d
get_subcolumn_mask (const int ncol,
                    const int nlay,
                    const int ngpt,
                    real2d &cldf,
                    const int overlap_option,
                    int1d &seeds);


void
compute_cloud_area (int ncol,
                    int nlay,
                    int ngpt,
                    real pmin,
                    real pmax,
                    const real2d& pmid,
                    const real3d& cld_tau_gpt,
                    real1d& cld_area);


void
compute_aerocom_cloudtop (int ncol, int nlay, const real2d &tmid, const real2d &pmid,
                          const real2d &p_del, const real2d &z_del, const real2d &qc,
                          const real2d &qi, const real2d &rel, const real2d &rei,
                          const real2d &cldfrac_tot, const real2d &nc,
                          real1d &T_mid_at_cldtop, real1d &p_mid_at_cldtop,
                          real1d &cldfrac_ice_at_cldtop, real1d &cldfrac_liq_at_cldtop,
                          real1d &cldfrac_tot_at_cldtop, real1d &cdnc_at_cldtop,
                          real1d &eff_radius_qc_at_cldtop, real1d &eff_radius_qi_at_cldtop);

int get_wavelength_index (OpticalProps &kdist, double wavelength);
int get_wavelength_index_sw (double wavelength);
int get_wavelength_index_lw (double wavelength);

} // namespace rrtmgp
#endif

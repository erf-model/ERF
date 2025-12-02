#ifndef ERF_RRTMGP_INTERFACE_H
#define ERF_RRTMGP_INTERFACE_H

#include <Kokkos_Core.hpp>

#include <AMReX_REAL.H>

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

#include <ERF_Kokkos.H>

template <typename T>
using view_t       = Kokkos::View<T, layout_t, KokkosDefaultDevice>;

using pool_t       = conv::MemPoolSingleton<RealT, layout_t, KokkosDefaultDevice>;

using gas_optics_t       =  GasOpticsRRTMGPK<RealT, layout_t, KokkosDefaultDevice>;
using cloud_optics_t     =      CloudOpticsK<RealT, layout_t, KokkosDefaultDevice>;
using gas_concs_t        =         GasConcsK<RealT, layout_t, KokkosDefaultDevice>;
using fluxes_t           =     FluxesBybandK<RealT, layout_t, KokkosDefaultDevice>;
using fluxes_broadband_t =  FluxesBroadbandK<RealT, layout_t, KokkosDefaultDevice>;
using optical_props_t    =     OpticalPropsK<RealT, layout_t, KokkosDefaultDevice>;
using optical_props1_t   = OpticalProps1sclK<RealT, layout_t, KokkosDefaultDevice>;
using optical_props2_t   = OpticalProps2strK<RealT, layout_t, KokkosDefaultDevice>;
using source_func_t      =     SourceFuncLWK<RealT, layout_t, KokkosDefaultDevice>;

namespace rrtmgp {

extern std::unique_ptr<gas_optics_t> k_dist_sw_k;
extern std::unique_ptr<gas_optics_t> k_dist_lw_k;

extern std::unique_ptr<cloud_optics_t> cloud_optics_sw_k;
extern std::unique_ptr<cloud_optics_t> cloud_optics_lw_k;

extern bool initialized;

optical_props2_t
get_cloud_optics_sw (const int ncol,
                     const int nlay,
                     cloud_optics_t& cloud_optics,
                     gas_optics_t& kdist,
                     real2d_k& lwp,
                     real2d_k& iwp,
                     real2d_k& rel,
                     real2d_k& rei);

optical_props1_t
get_cloud_optics_lw (const int ncol,
                     const int nlay,
                     cloud_optics_t& cloud_optics,
                     gas_optics_t& kdist,
                     real2d_k& lwp,
                     real2d_k& iwp,
                     real2d_k& rel,
                     real2d_k& rei);

optical_props2_t
get_subsampled_clouds (const int ncol,
                       const int nlay,
                       const int nbnd,
                       const int ngpt,
                       optical_props2_t& cloud_optics,
                       gas_optics_t& kdist,
                       real2d_k& cld,
                       real2d_k& p_lay);

optical_props1_t
get_subsampled_clouds (const int ncol,
                       const int nlay,
                       const int nbnd,
                       const int ngpt,
                       optical_props1_t& cloud_optics,
                       gas_optics_t& kdist,
                       real2d_k& cld,
                       real2d_k& p_lay);

void
rrtmgp_initialize (gas_concs_t& gas_concs,
                   const std::string& coefficients_file_sw, const std::string& coefficients_file_lw,
                   const std::string& cloud_optics_file_sw, const std::string& cloud_optics_file_lw);


void
rrtmgp_finalize ();


void
compute_band_by_band_surface_albedos (const int ncol, const int nswbands,
                                      real1d_k& sfc_alb_dir_vis, real1d_k& sfc_alb_dir_nir,
                                      real1d_k& sfc_alb_dif_vis, real1d_k& sfc_alb_dif_nir,
                                      real2d_k& sfc_alb_dir,     real2d_k& sfc_alb_dif);


void
compute_broadband_surface_fluxes (const int ncol, const int kbot, const int nswbands,
                                  real3d_k& sw_bnd_flux_dir , real3d_k& sw_bnd_flux_dif ,
                                  real1d_k& sfc_flux_dir_vis, real1d_k& sfc_flux_dir_nir,
                                  real1d_k& sfc_flux_dif_vis, real1d_k& sfc_flux_dif_nir);


void
rrtmgp_main (const int ncol, const int nlay,
             real2d_k& p_lay, real2d_k& t_lay,
             real2d_k& p_lev, real2d_k& t_lev,
             gas_concs_t& gas_concs,
             real2d_k& sfc_alb_dir, real2d_k& sfc_alb_dif, real1d_k& mu0,
             real1d_k& t_sfc      , real1d_k& sfc_emis   , real1d_k& lw_src,
             real2d_k& lwp        , real2d_k& iwp        ,
             real2d_k& rel        , real2d_k& rei        , real2d_k& cldfrac,
             real3d_k& aer_tau_sw , real3d_k& aer_ssa_sw , real3d_k& aer_asm_sw,
             real3d_k& aer_tau_lw ,
             real3d_k& cld_tau_sw_bnd, real3d_k& cld_tau_lw_bnd,
             real3d_k& cld_tau_sw_gpt, real3d_k& cld_tau_lw_gpt,
             real2d_k& sw_flux_up    , real2d_k& sw_flux_dn    , real2d_k& sw_flux_dn_dir,
             real2d_k& lw_flux_up    , real2d_k& lw_flux_dn    ,
             real2d_k& sw_clnclrsky_flux_up, real2d_k& sw_clnclrsky_flux_dn, real2d_k& sw_clnclrsky_flux_dn_dir,
             real2d_k& sw_clrsky_flux_up   , real2d_k& sw_clrsky_flux_dn   , real2d_k& sw_clrsky_flux_dn_dir,
             real2d_k& sw_clnsky_flux_up   , real2d_k& sw_clnsky_flux_dn   , real2d_k& sw_clnsky_flux_dn_dir,
             real2d_k& lw_clnclrsky_flux_up, real2d_k& lw_clnclrsky_flux_dn,
             real2d_k& lw_clrsky_flux_up   , real2d_k& lw_clrsky_flux_dn   ,
             real2d_k& lw_clnsky_flux_up   , real2d_k& lw_clnsky_flux_dn   ,
             real3d_k& sw_bnd_flux_up      , real3d_k& sw_bnd_flux_dn      , real3d_k& sw_bnd_flux_dn_dir,
             real3d_k& lw_bnd_flux_up      , real3d_k& lw_bnd_flux_dn,
             const RealT tsi_scaling,
             const bool extra_clnclrsky_diag = false, const bool extra_clnsky_diag = false);


void
rrtmgp_sw (const int ncol, const int nlay,
           gas_optics_t& k_dist,
           real2d_k& p_lay, real2d_k& t_lay,
           real2d_k& p_lev, real2d_k& t_lev,
           gas_concs_t& gas_concs,
           real2d_k& sfc_alb_dir, real2d_k& sfc_alb_dif, real1d_k& mu0,
           optical_props2_t& aerosol, optical_props2_t& clouds,
           fluxes_t& fluxes,
           fluxes_broadband_t& clnclrsky_fluxes,
           fluxes_broadband_t& clrsky_fluxes,
           fluxes_broadband_t& clnsky_fluxes,
           const RealT tsi_scaling,
           const bool extra_clnclrsky_diag,
           const bool extra_clnsky_diag);


void
rrtmgp_lw (const int ncol, const int nlay,
           gas_optics_t& k_dist,
           real2d_k& p_lay, real2d_k& t_lay,
           real2d_k& p_lev, real2d_k& t_lev,
           real1d_k& t_sfc, real1d_k& sfc_emis, real1d_k& lw_src,
           gas_concs_t& gas_concs,
           optical_props1_t& aerosol, optical_props1_t& clouds,
           fluxes_t& fluxes,
           fluxes_broadband_t& clnclrsky_fluxes,
           fluxes_broadband_t& clrsky_fluxes,
           fluxes_broadband_t& clnsky_fluxes,
           const bool extra_clnclrsky_diag,
           const bool extra_clnsky_diag);


int3d_k
get_subcolumn_mask (const int ncol,
                    const int nlay,
                    const int ngpt,
                    real2d_k& cldf,
                    const int overlap_option,
                    int1d_k& seeds);


void
compute_cloud_area (int ncol,
                    int nlay,
                    int ngpt,
                    RealT pmin,
                    RealT pmax,
                    const real2d_k& pmid,
                    const real3d_k& cld_tau_gpt,
                    real1d_k& cld_area);


void
compute_aerocom_cloudtop (int ncol, int nlay,
                          const real2d_k& tmid , const real2d_k& pmid,
                          const real2d_k& p_del, const real2d_k& z_del,
                          const real2d_k& qc   , const real2d_k& qi,
                          const real2d_k& rel  , const real2d_k& rei,
                          const real2d_k& cldfrac_tot, const real2d_k& nc,
                          real1d_k& T_mid_at_cldtop        , real1d_k& p_mid_at_cldtop,
                          real1d_k& cldfrac_ice_at_cldtop  , real1d_k& cldfrac_liq_at_cldtop,
                          real1d_k& cldfrac_tot_at_cldtop  , real1d_k& cdnc_at_cldtop,
                          real1d_k& eff_radius_qc_at_cldtop, real1d_k& eff_radius_qi_at_cldtop);

int get_wavelength_index (optical_props_t& kdist, RealT wavelength);
int get_wavelength_index_sw (RealT wavelength);
int get_wavelength_index_lw (RealT wavelength);

} // namespace rrtmgp
#endif

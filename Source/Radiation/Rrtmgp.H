/*
 * RTE-RRTMGP radiation model interface to ERF
 * The orginal code is developed by RobertPincus, and the code is open source available at:
 *                        https://github.com/earth-system-radiation/rte-rrtmgp
 * For details of the radiation algorithm, please reference to the following paper,
 *                        https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019MS001621
 *
 * NOTE: we use the C++ version of RTE-RRTMGP, which is reimplemented the original Fortran
 * code using C++ YAKL for CUDA, HiP and SYCL application by E3SM ECP team, the C++ version
 * of the rte-rrtmgp code is located at:
 *                       https://github.com/E3SM-Project/rte-rrtmgp
 *
 * The RTE-RRTMGP uses BSD-3-Clause Open Source License, if you want to make changes,
 * and modifications to the code, please refer to BSD-3-Clause Open Source License.
 */
#ifndef ERF_RRTMGP_H
#define ERF_RRTMGP_H

#include <string>
#include <vector>
#include <memory>

#include <ERF_Constants.H>

// rrtmgp includes
#include <mo_gas_concentrations.h>
#include <mo_gas_optics_rrtmgp.h>
#include <mo_load_coefficients.h>
#include <mo_rte_sw.h>
#include <mo_rte_lw.h>
#include <mo_optical_props.h>
#include <mo_fluxes_byband.h>
#include <rrtmgp_const.h>

// Radiation code interface class
class Rrtmgp {
  public:
    // constructor
    Rrtmgp() = default;

    // deconstructor
    ~Rrtmgp() = default;

    // initialize and load gas property data for rrtmgp radiation
    void initialize(int num_gas, const std::vector<std::string>& active_gas_names,
                    const char* rrtmgp_coefficients_file_sw,
                    const char* rrtmgp_coefficients_file_lw);

    // finalize/clean up
    void finalize();

    // run rrtmgp short wave model
    void run_shortwave_rrtmgp (int ngas, int ncol, int nlay,
        const real3d& gas_vmr, const real2d& pmid, const real2d& tmid, const real2d& pint,
        const real1d& coszrs , const real2d& albedo_dir, const real2d& albedo_dif,
        const real3d& cld_tau_gpt, const real3d& cld_ssa_gpt, const real3d& cld_asm_gpt,
        const real3d& aer_tau_bnd, const real3d& aer_ssa_bnd, const real3d& aer_asm_bnd,
        const real2d& allsky_flux_up    , const real2d& allsky_flux_dn    , const real2d& allsky_flux_net    , const real2d& allsky_flux_dn_dir,
        const real3d& allsky_bnd_flux_up, const real3d& allsky_bnd_flux_dn, const real3d& allsky_bnd_flux_net, const real3d& allsky_bnd_flux_dn_dir,
        const real2d& clrsky_flux_up    , const real2d& clrsky_flux_dn    , const real2d& clrsky_flux_net    , const real2d& clrsky_flux_dn_dir,
        const real3d& clrsky_bnd_flux_up, const real3d& clrsky_bnd_flux_dn, const real3d& clrsky_bnd_flux_net, const real3d& clrsky_bnd_flux_dn_dir,
        double tsi_scaling);

    // run rrtmgp long wave model
    void run_longwave_rrtmgp (int ngas, int ncol, int nlay, const real3d& gas_vmr,
        const real2d& pmid, const real2d& tmid, const real2d& pint , const real2d& tint,
        const real2d& emis_sfc,
        const real3d& cld_tau_gpt       , const real3d& aer_tau_bnd       ,
        const real2d& allsky_flux_up    , const real2d& allsky_flux_dn    , const real2d& allsky_flux_net    ,
        const real3d& allsky_bnd_flux_up, const real3d& allsky_bnd_flux_dn, const real3d& allsky_bnd_flux_net,
        const real2d& clrsky_flux_up    , const real2d& clrsky_flux_dn    , const real2d& clrsky_flux_net    ,
        const real3d& clrsky_bnd_flux_up, const real3d& clrsky_bnd_flux_dn, const real3d& clrsky_bnd_flux_net);

    int get_nband_sw() {
       return k_dist_sw.get_nband();
    }

    int get_nband_lw() {
       return k_dist_lw.get_nband();
    }

    int get_ngpt_sw() {
       return k_dist_sw.get_ngpt();
    }

    int get_ngpt_lw() {
       return k_dist_lw.get_ngpt();
    }

    double get_min_temperature() {
      return std::min(k_dist_sw.temp_ref_min, k_dist_lw.temp_ref_min);
    }

    double get_max_temperature() {
       return std::max(k_dist_sw.temp_ref_max, k_dist_lw.temp_ref_max);
    }

    void get_gpoint_bands_sw(int1d gpoint_bands) {
       gpoint_bands = k_dist_sw.get_gpoint_bands();
       yakl::fence();
    }

    void get_gpoint_bands_lw(int1d& gpoint_bands) {
       gpoint_bands = k_dist_lw.get_gpoint_bands();
       yakl::fence();
    }

  private:
   // tools code for multifab to yakl array
   template<typename DType>
   void multifab_to_array();

   template<typename DType>
   void array_to_multifab();

   // number of gas for radiation model
   int ngas;
   string1d active_gases;

   // coefficient files
   char const *coefficients_file_sw;  // short wave gas optics coefficient files
   char const *coefficients_file_lw;  // long wave gas optics coefficient files

   // Objects for gas optics data
   GasOpticsRRTMGP k_dist_sw;
   GasOpticsRRTMGP k_dist_lw;
};

#endif // ERF_RRTMGP_H

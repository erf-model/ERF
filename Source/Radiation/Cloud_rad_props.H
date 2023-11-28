//
// cloud radiation properties class
//
#ifndef ERF_CLOUD_RAD_PROPS_H_
#define ERF_CLOUD_RAD_PROPS_H_

#include <AMReX_MultiFabUtil.H>
#include "ERF_Constants.H"
#include "rrtmgp_const.h"
#include "Water_vapor_saturation.H"
#include "Linear_interpolate.H"

class CloudRadProps {
  public:
     // constructors and descturctors
     CloudRadProps() = default;
     ~CloudRadProps() = default;

     // initialization
     void initialize();

     // liquid optics
     void gammadist_liq_optics_sw(const int& ncol,
                                  const int& nlev,
                                  const real2d& iclwpth,
                                  const real2d& lamc,
                                  const real2d& pgam,
                                  real3d& tau,
                                  real3d& tau_w,
                                  real3d& tau_w_g,
                                  real3d& tau_w_f);

     void gammadist_liq_optics_lw(const int& ncol,
                                  const int& nlev,
                                  const real2d& iclwpth,
                                  const real2d& lamc,
                                  const real2d& pgam,
                                  real3d& abs_od);

    void mitchell_ice_optics_sw(const int& ncol,
                                const int& nlev,
                                const real2d& iciwpth,
                                const real2d& dei,
                                real3d& tau,
                                real3d& tau_w,
                                real3d& tau_w_g,
                                real3d& tau_w_f);

    void mitchell_ice_optics_lw(const int& ncol,
                                const int& nlev,
                                const real2d& iciwpth,
                                const real2d& dei,
                                real3d& abs_od);

    void gam_liquid_lw(const real& clwptn,
                       const real& lamc,
                       const real& pgam,
                       real1d abs_od);


    // liquid optics (short wave)
    void gam_liquid_sw(const real& clwptn,
                       const real& lamc,
                       const real& pgam,
                       real1d tau,
                       real1d tau_w,
                       real1d tau_w_g,
                       real1d tau_w_f);

    // lambda weight
    void get_mu_lambda_weights(const real& lamc,
                               const real& pgam,
                               LinInterp::InterpType& mu_wgts,
                               LinInterp::InterpType& lambda_wgts);
  private:
    std::string name{"CloudRadProps"};
    std::string liquid_file;
    std::string ice_file;

    int nlwbands, nswbands, nlambda, nmu, n_g_d;

    real1d g_mu;         // mu samples on grid
    real2d g_lambda;     // lambda scale samples on grid
    real3d ext_sw_liq;
    real3d ssa_sw_liq;
    real3d asm_sw_liq;
    real3d abs_lw_liq;

    real1d g_d_eff;       // radiative effective diameter samples on grid
    real2d ext_sw_ice;
    real2d ssa_sw_ice;
    real2d asm_sw_ice;
    real2d abs_lw_ice;
};
#endif // ERF_CLOUDRADPROPS_H

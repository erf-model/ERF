//
// get the cloud radiation props data
//
#include <cstdlib>
#include "YAKL_netcdf.h"
#include "Cloud_rad_props.H"
using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

void CloudRadProps::initialize ()
{
    realHost1d g_mu_h;         // mu samples on grid
    realHost2d g_lambda_h;     // lambda scale samples on grid
    realHost3d ext_sw_liq_h;
    realHost3d ssa_sw_liq_h;
    realHost3d asm_sw_liq_h;
    realHost3d abs_lw_liq_h;

    realHost1d g_d_eff_h;       // radiative effective diameter samples on grid
    realHost2d ext_sw_ice_h;
    realHost2d ssa_sw_ice_h;
    realHost2d asm_sw_ice_h;
    realHost2d abs_lw_ice_h;

    const char* erf_home = std::getenv("ERF_HOME");
    auto erf_rad_data_dir = (std::string)erf_home + "/Source/Radiation/data/";
    liquid_file = erf_rad_data_dir + "F_nwvl200_mu20_lam50_res64_t298_c080428.nc";
    ice_file = erf_rad_data_dir + "iceoptics_c080917.nc";

    yakl::SimpleNetCDF liquid, ice;
    liquid.open(liquid_file, yakl::NETCDF_MODE_READ);
    ice.open(ice_file, yakl::NETCDF_MODE_READ);

    // read the dimensions
    nlwbands = liquid.getDimSize( "lw_band" );
    nswbands = liquid.getDimSize( "sw_band" );
    nmu      = liquid.getDimSize( "mu"     );
    nlambda  = liquid.getDimSize( "lambda_scale" );

    liquid.read( g_mu_h, "mu");
    liquid.read( g_lambda_h, "lambda");
    liquid.read( ext_sw_liq_h, "k_ext_sw");
    liquid.read( ssa_sw_liq_h, "ssa_sw");
    liquid.read( asm_sw_liq_h, "asm_sw");
    liquid.read( abs_lw_liq_h, "k_abs_lw");

    g_mu = real1d("gmu", nmu);
    g_lambda = real2d("g_lambda", nmu, nlambda);
    ext_sw_liq = real3d("ext_sw_liq", nmu, nlambda, nswbands);
    ssa_sw_liq = real3d("ssa_sw_liq", nmu, nlambda, nswbands);
    asm_sw_liq = real3d("asm_sw_liq", nmu, nlambda, nswbands);
    abs_lw_liq = real3d("abs_lw_liq", nmu, nlambda, nlwbands);

    g_mu_h.deep_copy_to(g_mu);
    g_lambda_h.deep_copy_to(g_lambda);
    ext_sw_liq_h.deep_copy_to(ext_sw_liq);
    ssa_sw_liq_h.deep_copy_to(ssa_sw_liq);
    asm_sw_liq_h.deep_copy_to(asm_sw_liq);
    abs_lw_liq_h.deep_copy_to(abs_lw_liq);

    // I forgot to convert kext from m^2/Volume to m^2/Kg
    parallel_for(SimpleBounds<3>(nmu,nlambda,nswbands) , YAKL_LAMBDA (int i, int j, int k)
    {
        ext_sw_liq(i,j,k) = ext_sw_liq(i,j,k) / 0.9970449e3;
    });

    parallel_for(SimpleBounds<3>(nmu,nlambda,nlwbands) , YAKL_LAMBDA (int i, int j, int k)
    {
        abs_lw_liq(i,j,k) = abs_lw_liq(i,j,k) / 0.9970449e3;
    });

    // read ice cloud optics
    nlwbands = ice.getDimSize( "lw_band" );
    nswbands = ice.getDimSize( "sw_band" );
    n_g_d    = ice.getDimSize( "d_eff" );

    ice.read( g_d_eff_h,    "d_eff");
    ice.read( ext_sw_ice_h, "sw_ext");
    ice.read( ssa_sw_ice_h, "sw_ssa");
    ice.read( asm_sw_ice_h, "sw_asm");
    ice.read( abs_lw_ice_h, "lw_abs");

    g_d_eff = real1d("g_d_eff",n_g_d);
    ext_sw_ice = real2d("ext_sw_ice", n_g_d, nswbands);
    ssa_sw_ice = real2d("ssa_sw_ice", n_g_d, nswbands);
    asm_sw_ice = real2d("asm_sw_ice", n_g_d, nswbands);
    abs_lw_ice = real2d("abs_lw_ice", n_g_d, nlwbands);

    g_d_eff_h.deep_copy_to(g_d_eff);
    ext_sw_ice_h.deep_copy_to(ext_sw_ice);
    ssa_sw_ice_h.deep_copy_to(ssa_sw_ice);
    asm_sw_ice_h.deep_copy_to(asm_sw_ice);
    abs_lw_ice_h.deep_copy_to(abs_lw_ice);
}

void CloudRadProps::gammadist_liq_optics_sw (const int& ncol,
                                             const int& nlev,
                                             const real2d& iclwpth,
                                             const real2d& lamc,
                                             const real2d& pgam,
                                             real3d& tau,
                                             real3d& tau_w,
                                             real3d& tau_w_g,
                                             real3d& tau_w_f)
{
    real1d tau1d("tau1d", nswbands);
    real1d tauw1d("tauw1d", nswbands);
    real1d tauwg1d("tauwg1d", nswbands);
    real1d tauwf1d("tauwf1d", nswbands);
    parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int i, int k)
    {
        if(lamc(i,k) > 0.) { // This seems to be clue from microphysics of no cloud
            for (auto iband=1; iband<=nswbands; ++iband) {
                tau1d(iband)   = tau(iband,i,k);
                tauw1d(iband)  = tau_w(iband,i,k);
                tauwg1d(iband) = tau_w_g(iband,i,k);
                tauwf1d(iband) = tau_w_f(iband,i,k);
            }
            gam_liquid_sw(iclwpth(i,k), lamc(i,k), pgam(i,k), tau1d, tauw1d, tauwg1d, tauwf1d);
        }
        else {
            for (auto iband=1; iband<=nswbands; ++iband) {
                tau(iband,i,k) = 0.;
                tau_w(iband,i,k) = 0.;
                tau_w_g(iband,i,k) = 0.;
                tau_w_f(iband,i,k) = 0.;
            }
        }
    });
}

void CloudRadProps::gammadist_liq_optics_lw (const int& ncol,
                                             const int& nlev,
                                             const real2d& iclwpth,
                                             const real2d& lamc,
                                             const real2d& pgam,
                                             real3d& abs_od)
{
    auto abs_od_1d = real1d("abs_od_1d", nlwbands);
    parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int i, int k)
    {
        if(lamc(i,k) > 0.) { // This seems to be the clue for no cloud from microphysics formulation
            for (auto ib=1; ib<=nlwbands; ++ib) abs_od_1d(ib) = abs_od(ib,i,k);
            gam_liquid_lw(iclwpth(i,k), lamc(i,k), pgam(i,k), abs_od_1d);
        }
        else {
            for(auto j=1; j<=nlwbands; ++j) abs_od(j,i,k) = 0.;
        }
    });
}


void CloudRadProps::mitchell_ice_optics_sw (const int& ncol,
                                            const int& nlev,
                                            const real2d& iciwpth,
                                            const real2d& dei,
                                            real3d& tau,
                                            real3d& tau_w,
                                            real3d& tau_w_g,
                                            real3d& tau_w_f)
{
    real1d ext("ext", nswbands),
        ssa("ssa", nswbands),
        assm("assm", nswbands);

    real1d ext_sw_ice_1d("ext_sw_ice_1d", n_g_d),
        ssa_sw_ice_1d("ssa_sw_ice_1d", n_g_d),
        asm_sw_ice_1d("asm_sw_ice_1d", n_g_d);

    real1d dei_1d("dei_1d",1);
    auto ext_1d  = real1d("ext_1d",1);
    auto ssa_1d  = real1d("ssa_1d",1);
    auto assm_1d = real1d("assm_1d",1);

    parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int i, int k)
    {
        if( iciwpth(i,k) < 1.e-80 || dei(i,k) == 0.) {
            // if ice water path is too small
            for(auto ib=1; ib<=nswbands; ++ib) {
                tau    (ib,i,k) = 0.;
                tau_w  (ib,i,k) = 0.;
                tau_w_g(ib,i,k) = 0.;
                tau_w_f(ib,i,k) = 0.;
            }
        }
        else {
            // for each cell interpolate to find weights in g_d_eff grid.
            dei_1d(1) = dei(i,k);
            LinInterp::InterpType dei_wgts;
            LinInterp::init(g_d_eff, n_g_d, dei_1d, 1, LinInterp::extrap_method_bndry, dei_wgts);
            // interpolate into grid and extract radiative properties
            for (auto is=1; is<=nswbands; ++is) {
                for(auto ig=1; ig<=n_g_d; ++ig) {
                    ext_sw_ice_1d(ig) = ext_sw_ice(ig, is);
                    ssa_sw_ice_1d(ig) = ssa_sw_ice(ig, is);
                    asm_sw_ice_1d(ig) = asm_sw_ice(ig, is);
                }
                ext_1d(1)  = ext(is);
                ssa_1d(1)  = ssa(is);
                assm_1d(1) = assm(is);
                LinInterp::interp1d(ext_sw_ice_1d, n_g_d, ext_1d,  1, dei_wgts);
                LinInterp::interp1d(ssa_sw_ice_1d, n_g_d, ssa_1d,  1, dei_wgts);
                LinInterp::interp1d(asm_sw_ice_1d, n_g_d, assm_1d, 1, dei_wgts);
            }
            for (auto is=1; is<=nswbands; ++is) {
                tau    (is,i,k) = iciwpth(i,k) * ext(is);
                tau_w  (is,i,k) = tau(is,i,k) * ssa(is);
                tau_w_g(is,i,k) = tau_w(is,i,k) * assm(is);
                tau_w_f(is,i,k) = tau_w_g(is,i,k) * assm(is);
            }
        }
    });
}

void CloudRadProps::mitchell_ice_optics_lw (const int& ncol,
                                            const int& nlev,
                                            const real2d& iciwpth,
                                            const real2d& dei,
                                            real3d& abs_od)
{
    real1d absor("absor", nlwbands);
    real1d dei_1d("dei_1d",1);
    real1d abs_lw_ice_1d("abs_lw_ice_1d",n_g_d);

    real1d absor_1d("absor_1d",1);

    parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int i, int k)
    {
        // if ice water path is too small, OD := 0
        if(iciwpth(i,k) < 1.e-80 || dei(i,k) == 0.) {
            for (auto lb=1; lb<=nlwbands; ++lb) {
                abs_od (lb,i,k) = 0.;
            }
        }
        else {
            // for each cell interpolate to find weights in g_d_eff grid.
            dei_1d(1) = dei(i,k);
            LinInterp::InterpType dei_wgts;
            LinInterp::init(g_d_eff, n_g_d, dei_1d, 1, LinInterp::extrap_method_bndry, dei_wgts);
            // interpolate into grid and extract radiative properties
            for(auto lwband = 1; lwband <= nlwbands; ++lwband) {
                absor_1d(1) = absor(lwband);
                LinInterp::interp1d(abs_lw_ice_1d, n_g_d, absor_1d, 1, dei_wgts);
                absor(lwband) = absor_1d(1);
            }
            for (auto lwband = 1; lwband <= nlwbands; ++lwband) {
                abs_od(lwband,i,k) = iciwpth(i,k) * absor(lwband);
            }
        }
    });
}


void CloudRadProps::gam_liquid_lw (const real& clwptn,
                                   const real& lamc,
                                   const real& pgam,
                                   real1d abs_od)
{
    LinInterp::InterpType mu_wgts;
    LinInterp::InterpType lambda_wgts;

    real2d abs_lw_liq_2d("abs_lw_liq_2d", nmu, nlambda);
    real1d abs_od_1d("abs_od_1d", 1);

    if (clwptn < 1.e-80) {
        yakl::memset(abs_od, 0.);
        return;
    }

    get_mu_lambda_weights(lamc, pgam, mu_wgts, lambda_wgts);

    parallel_for(SimpleBounds<1>(nlwbands), YAKL_LAMBDA (int lwband)
    {
        for (auto imu=1; imu<=nmu; ++imu) {
            for (auto ilb=1; ilb<=nlambda; ++ilb) {
                abs_lw_liq_2d(imu, ilb) = abs_lw_liq(imu, ilb, lwband);
            }
        }
        abs_od_1d(1) = abs_od(lwband);
        LinInterp::interp2d1d(abs_lw_liq_2d, nmu, nlambda, abs_od_1d, 1, mu_wgts, lambda_wgts);
        abs_od(lwband) = clwptn*abs_od_1d(1);
    });

    parallel_for(SimpleBounds<1>(nlwbands) , YAKL_LAMBDA (int iband)
    {
        abs_od(iband) = clwptn * abs_od(iband);
    });
}

void CloudRadProps::gam_liquid_sw (const real& clwptn,
                                   const real& lamc,
                                   const real& pgam,
                                   real1d tau,
                                   real1d tau_w,
                                   real1d tau_w_g,
                                   real1d tau_w_f)
{
    real1d extb("ext", nswbands),
           ssab("ssa", nswbands),
           asmb("asm", nswbands);

    real2d ext_sw_liq_2d("ext_sw_liq_2d", nmu, nlambda),
           ssa_sw_liq_2d("ssa_sw_liq_2d", nmu, nlambda),
           asm_sw_liq_2d("asm_sw_liq_2d", nmu, nlambda);

    LinInterp::InterpType mu_wgts;
    LinInterp::InterpType lambda_wgts;

    if (clwptn < 1.e-80) {
        yakl::memset(tau, 0.);
        yakl::memset(tau_w, 0.);
        yakl::memset(tau_w_g, 0.);
        yakl::memset(tau_w_f, 0.);
        return;
    }

    get_mu_lambda_weights(lamc, pgam, mu_wgts, lambda_wgts);

    parallel_for(SimpleBounds<1>(nswbands), YAKL_LAMBDA (int swband)
    {
        for (auto imu=1; imu<=nmu; ++imu) {
            for (auto lb=1; lb<=nlambda; ++lb) {
                ext_sw_liq_2d(imu,lb) = ext_sw_liq(imu,lb,swband);
                ssa_sw_liq_2d(imu,lb) = ssa_sw_liq(imu,lb,swband);
                asm_sw_liq_2d(imu,lb) = asm_sw_liq(imu,lb,swband);
            }
        }
        LinInterp::interp2d1d(ext_sw_liq_2d, nmu, nlambda, extb, 1, mu_wgts, lambda_wgts);
        LinInterp::interp2d1d(ssa_sw_liq_2d, nmu, nlambda, ssab, 1, mu_wgts, lambda_wgts);
        LinInterp::interp2d1d(asm_sw_liq_2d, nmu, nlambda, asmb, 1, mu_wgts, lambda_wgts);
    });

  // compute radiative properties
  parallel_for(SimpleBounds<1>(nswbands), YAKL_LAMBDA (int iband)
  {
      tau(iband)     = clwptn * extb(iband);
      tau_w(iband)   = tau(iband) * ssab(iband);
      tau_w_g(iband) = tau_w(iband) * asmb(iband);
      tau_w_f(iband) = tau_w_g(iband) * asmb(iband);
  });
}


void CloudRadProps::get_mu_lambda_weights (const real& lamc,
                                           const real& pgam,
                                           LinInterp::InterpType& mu_wgts,
                                           LinInterp::InterpType& lambda_wgts)
{
    real1d g_lambda_interp("g_lambda_interp", nlambda);
    real1d pgam1d("pgam1d", nmu);
    real1d lamc1d("lam1d", nmu);
    real1d g_lambda1d("g_lambda1d", nmu);

    yakl::memset(pgam1d, pgam);
    yakl::memset(lamc1d, lamc);

    real1d g_mu0 = real1d("g_mu0",nmu);

    // Make interpolation weights for mu.
    // (Put pgam in a temporary array for this purpose.)
    LinInterp::init(g_mu0, nmu, pgam1d, 1, LinInterp::extrap_method_bndry, mu_wgts);

    // Use mu weights to interpolate to a row in the lambda table.
    parallel_for(SimpleBounds<1>(nlambda), YAKL_LAMBDA (int i)
    {
        for (auto im=1; im<=nmu; ++im) g_lambda1d(im) = g_lambda(im, i);
        LinInterp::interp1d(g_lambda1d, nmu, g_lambda_interp, 1, mu_wgts);
    });

    // Make interpolation weights for lambda.
    LinInterp::init(g_lambda_interp, nlambda, lamc1d, 1, LinInterp::extrap_method_bndry, lambda_wgts);
}


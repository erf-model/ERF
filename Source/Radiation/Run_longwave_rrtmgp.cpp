#include "mo_gas_concentrations.h"
#include "mo_gas_optics_rrtmgp.h"
#include "mo_load_coefficients.h"
#include "mo_rte_sw.h"
#include "mo_rte_lw.h"
#include "mo_optical_props.h"
#include "mo_fluxes_byband.h"
#include "Rrtmgp.H"

void Rrtmgp::run_longwave_rrtmgp (int ngas, int ncol, int nlay,
                                  const real3d& gas_vmr,
                                  const real2d& pmid, const real2d& tmid, const real2d& pint, const real2d& tint,
                                  const real2d& emis_sfc          ,
                                  const real3d& cld_tau_gpt       , const real3d& aer_tau_bnd       ,
                                  const real2d& allsky_flux_up    , const real2d& allsky_flux_dn    , const real2d& allsky_flux_net    ,
                                  const real3d& allsky_bnd_flux_up, const real3d& allsky_bnd_flux_dn, const real3d& allsky_bnd_flux_net,
                                  const real2d& clrsky_flux_up    , const real2d& clrsky_flux_dn    , const real2d& clrsky_flux_net    ,
                                  const real3d& clrsky_bnd_flux_up, const real3d& clrsky_bnd_flux_dn, const real3d& clrsky_bnd_flux_net)
{
    // Wrap pointers in YAKL arrays
    int nlwbands = k_dist_lw.get_nband();
    int nlwgpts  = k_dist_lw.get_ngpt();

    // Populate gas concentrations
    GasConcs gas_concs;
    gas_concs.init(active_gases, ncol, nlay);
    real2d tmp2d;
    tmp2d = real2d("tmp", ncol, nlay);
    for (int igas = 1; igas <= ngas; igas++) {
        parallel_for(SimpleBounds<2>(nlay,ncol), YAKL_LAMBDA(int ilay, int icol)
        {
            tmp2d(icol,ilay) = gas_vmr(igas,icol,ilay);
        });
        gas_concs.set_vmr(active_gases(igas), tmp2d);
    }

    //  Boundary conditions
    SourceFuncLW lw_sources;
    lw_sources.alloc(ncol, nlay, k_dist_lw);

    // Weights and angle secants for first order (k=1) Gaussian quadrature.
    //   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
    //   after Abramowitz & Stegun 1972, page 921
    int constexpr max_gauss_pts = 4;
    realHost2d gauss_Ds_host ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
    gauss_Ds_host(1,1) = 1.66_wp      ; gauss_Ds_host(2,1) =         0._wp; gauss_Ds_host(3,1) =         0._wp; gauss_Ds_host(4,1) =         0._wp;
    gauss_Ds_host(1,2) = 1.18350343_wp; gauss_Ds_host(2,2) = 2.81649655_wp; gauss_Ds_host(3,2) =         0._wp; gauss_Ds_host(4,2) =         0._wp;
    gauss_Ds_host(1,3) = 1.09719858_wp; gauss_Ds_host(2,3) = 1.69338507_wp; gauss_Ds_host(3,3) = 4.70941630_wp; gauss_Ds_host(4,3) =         0._wp;
    gauss_Ds_host(1,4) = 1.06056257_wp; gauss_Ds_host(2,4) = 1.38282560_wp; gauss_Ds_host(3,4) = 2.40148179_wp; gauss_Ds_host(4,4) = 7.15513024_wp;

    realHost2d gauss_wts_host("gauss_wts",max_gauss_pts,max_gauss_pts);
    gauss_wts_host(1,1) = 0.5_wp         ; gauss_wts_host(2,1) = 0._wp          ; gauss_wts_host(3,1) = 0._wp          ; gauss_wts_host(4,1) = 0._wp          ;
    gauss_wts_host(1,2) = 0.3180413817_wp; gauss_wts_host(2,2) = 0.1819586183_wp; gauss_wts_host(3,2) = 0._wp          ; gauss_wts_host(4,2) = 0._wp          ;
    gauss_wts_host(1,3) = 0.2009319137_wp; gauss_wts_host(2,3) = 0.2292411064_wp; gauss_wts_host(3,3) = 0.0698269799_wp; gauss_wts_host(4,3) = 0._wp          ;
    gauss_wts_host(1,4) = 0.1355069134_wp; gauss_wts_host(2,4) = 0.2034645680_wp; gauss_wts_host(3,4) = 0.1298475476_wp; gauss_wts_host(4,4) = 0.0311809710_wp;

    real2d gauss_Ds ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
    real2d gauss_wts("gauss_wts",max_gauss_pts,max_gauss_pts);
    gauss_Ds_host .deep_copy_to(gauss_Ds );
    gauss_wts_host.deep_copy_to(gauss_wts);

    // Populate optical property objects
    OpticalProps1scl combined_optics;
    combined_optics.alloc_1scl(ncol, nlay, k_dist_lw);
    bool1d top_at_1_g("top_at_1_g",1);
    boolHost1d top_at_1_h("top_at_1_h",1);
    bool top_at_1;
    real1d t_sfc("t_sfc", ncol);
    parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA (int icol)
    {
        t_sfc(icol) = tint(icol,nlay+1);
        top_at_1_g(1) = pmid(1, 1) < pmid (1, 2);
    });
    top_at_1_g.deep_copy_to(top_at_1_h);
    top_at_1 = top_at_1_h(1);
    k_dist_lw.gas_optics(ncol, nlay, top_at_1, pmid, pint, tmid, t_sfc, gas_concs, combined_optics, lw_sources, real2d(), tint);

    // Add in aerosol; we can define this by bands or gpoints. If we define by
    // bands, then internally when increment() is called it will map these to
    // gpoints. Not sure if there is a beneift one way or another.
    OpticalProps1scl aerosol_optics;
    auto &aerosol_optics_tau = aerosol_optics.tau;
    if (false) {
        aerosol_optics.alloc_1scl(ncol, nlay, k_dist_lw);
        auto gpt_bnd = aerosol_optics.get_gpoint_bands();
        parallel_for(SimpleBounds<3>(nlwgpts,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol)
        {
            aerosol_optics_tau(icol,ilay,igpt) = aer_tau_bnd(icol,ilay,gpt_bnd(igpt));
        });
    } else {
        aerosol_optics.alloc_1scl(ncol, nlay, k_dist_lw.get_band_lims_wavenumber());
        parallel_for(SimpleBounds<3>(nlwbands,nlay,ncol), YAKL_LAMBDA (int ibnd, int ilay, int icol)
        {
            aerosol_optics_tau(icol,ilay,ibnd) = aer_tau_bnd(icol,ilay,ibnd);
        });
    }
    aerosol_optics.increment(combined_optics);

    // Do the clearsky calculation before adding in clouds
    FluxesByband fluxes_clrsky;
    fluxes_clrsky.flux_up  = real2d("clrsky_flux_up" , ncol, nlay+1); // clrsky_flux_up;
    fluxes_clrsky.flux_dn  = real2d("clrsky_flux_dn" , ncol, nlay+1); //clrsky_flux_dn;
    fluxes_clrsky.flux_net = real2d("clrsky_flux_net", ncol, nlay+1); //clrsky_flux_net;
    fluxes_clrsky.bnd_flux_up  = real3d("clrsky_bnd_flux_up" , ncol, nlay+1, nlwbands); //clrsky_bnd_flux_up;
    fluxes_clrsky.bnd_flux_dn  = real3d("clrsky_bnd_flux_dn" , ncol, nlay+1, nlwbands); //clrsky_bnd_flux_dn;
    fluxes_clrsky.bnd_flux_net = real3d("clrsky_bnd_flux_net", ncol, nlay+1, nlwbands); //clrsky_bnd_flux_net;

    rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, combined_optics, top_at_1, lw_sources, emis_sfc, fluxes_clrsky);

    // Copy fluxes back out of FluxesByband object
    fluxes_clrsky.flux_up.deep_copy_to(clrsky_flux_up);
    fluxes_clrsky.flux_dn.deep_copy_to(clrsky_flux_dn);
    fluxes_clrsky.flux_net.deep_copy_to(clrsky_flux_net);
    fluxes_clrsky.bnd_flux_up.deep_copy_to(clrsky_bnd_flux_up);
    fluxes_clrsky.bnd_flux_dn.deep_copy_to(clrsky_bnd_flux_dn);
    fluxes_clrsky.bnd_flux_net.deep_copy_to(clrsky_bnd_flux_net);

    // Add in clouds
    OpticalProps1scl cloud_optics;
    cloud_optics.alloc_1scl(ncol, nlay, k_dist_lw);
    auto &cloud_optics_tau = cloud_optics.tau;
    parallel_for(SimpleBounds<3>(nlwgpts,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
        cloud_optics_tau(icol,ilay,igpt) = cld_tau_gpt(icol,ilay,igpt);
    });
    cloud_optics.increment(combined_optics);

    // Call LW flux driver
    FluxesByband fluxes_allsky;
    fluxes_allsky.flux_up  = real2d("allsky_flux_up" , ncol, nlay+1); //allsky_flux_up;
    fluxes_allsky.flux_dn  = real2d("allsky_flux_dn" , ncol, nlay+1); //allsky_flux_dn;
    fluxes_allsky.flux_net = real2d("allsky_flux_net", ncol, nlay+1); //allsky_flux_net;
    fluxes_allsky.bnd_flux_up  = real3d("allsky_bnd_flux_up" , ncol, nlay+1, nlwbands); //allsky_bnd_flux_up;
    fluxes_allsky.bnd_flux_dn  = real3d("allsky_bnd_flux_dn" , ncol, nlay+1, nlwbands); //allsky_bnd_flux_dn;
    fluxes_allsky.bnd_flux_net = real3d("allsky_bnd_flux_net", ncol, nlay+1, nlwbands); //allsky_bnd_flux_net;

    rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, combined_optics, top_at_1, lw_sources, emis_sfc, fluxes_allsky);

    // Copy fluxes back out of FluxesByband object
    fluxes_allsky.flux_up.deep_copy_to(allsky_flux_up);
    fluxes_allsky.flux_dn.deep_copy_to(allsky_flux_dn);
    fluxes_allsky.flux_net.deep_copy_to(allsky_flux_net);
    fluxes_allsky.bnd_flux_up.deep_copy_to(allsky_bnd_flux_up);
    fluxes_allsky.bnd_flux_dn.deep_copy_to(allsky_bnd_flux_dn);
    fluxes_allsky.bnd_flux_net.deep_copy_to(allsky_bnd_flux_net);
    yakl::fence();
}

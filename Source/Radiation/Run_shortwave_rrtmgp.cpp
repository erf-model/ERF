#include "Rrtmgp.H"

void Rrtmgp::run_shortwave_rrtmgp (
        int ngas, int ncol, int nlay,
        real3d& gas_vmr, real2d& pmid      , real2d& tmid      , real2d& pint,
        real1d& coszrs , real2d& albedo_dir, real2d& albedo_dif,
        real3d& cld_tau_gpt, real3d& cld_ssa_gpt, real3d& cld_asm_gpt,
        real3d& aer_tau_bnd, real3d& aer_ssa_bnd, real3d& aer_asm_bnd,
        real2d& allsky_flux_up    , real2d& allsky_flux_dn    , real2d& allsky_flux_net    , real2d& allsky_flux_dn_dir,
        real3d& allsky_bnd_flux_up, real3d& allsky_bnd_flux_dn, real3d& allsky_bnd_flux_net, real3d& allsky_bnd_flux_dn_dir,
        real2d& clrsky_flux_up    , real2d& clrsky_flux_dn    , real2d& clrsky_flux_net    , real2d& clrsky_flux_dn_dir,
        real3d& clrsky_bnd_flux_up, real3d& clrsky_bnd_flux_dn, real3d& clrsky_bnd_flux_net, real3d& clrsky_bnd_flux_dn_dir,
        double tsi_scaling)
{

    // Wrap pointers in YAKL arrays
    int nswbands = k_dist_sw.get_nband();
    int nswgpts  = k_dist_sw.get_ngpt();

    // Populate gas concentrations object
    GasConcs gas_concs;
    gas_concs.init(active_gases, ncol, nlay);
    real2d tmp2d;
    tmp2d = real2d("tmp", ncol, nlay);
    for (int igas = 1; igas <= ngas; igas++) {
        yakl::c::parallel_for(yakl::c::Bounds<2>(nlay,ncol), YAKL_LAMBDA(int ilay, int icol) {
            tmp2d(icol,ilay) = gas_vmr(igas,icol,ilay);
        });
        gas_concs.set_vmr(active_gases(igas), tmp2d);
    }

    // Do gas optics
    // TODO: should we avoid allocating here?
    OpticalProps2str combined_optics;
    combined_optics.alloc_2str(ncol, nlay, k_dist_sw);
    bool* top_at_1;
    real2d toa_flux("toa_flux", ncol, nswgpts);
    k_dist_sw.gas_optics(ncol, nlay, top_at_1, pmid, pint, tmid, gas_concs, combined_optics, toa_flux);

    // Apply TOA flux scaling
    yakl::c::parallel_for(yakl::c::Bounds<2>(nswgpts,ncol), YAKL_LAMBDA (int igpt, int icol) {
        toa_flux(icol, igpt) = tsi_scaling * toa_flux(icol, igpt);
        *top_at_1 = pmid(1, 1) < pmid (1, 2);
    });

    // Add in aerosol
    // TODO: should we avoid allocating here?
    OpticalProps2str aerosol_optics;
    auto &aerosol_optics_tau = aerosol_optics.tau;
    auto &aerosol_optics_ssa = aerosol_optics.ssa;
    auto &aerosol_optics_g   = aerosol_optics.g  ;
    if (true) {
        aerosol_optics.alloc_2str(ncol, nlay, k_dist_sw);
        auto gpt_bnd = aerosol_optics.get_gpoint_bands();
        yakl::c::parallel_for(yakl::c::Bounds<3>(nswgpts,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
            aerosol_optics_tau(icol,ilay,igpt) = aer_tau_bnd(icol,ilay,gpt_bnd(igpt));
            aerosol_optics_ssa(icol,ilay,igpt) = aer_ssa_bnd(icol,ilay,gpt_bnd(igpt));
            aerosol_optics_g  (icol,ilay,igpt) = aer_asm_bnd(icol,ilay,gpt_bnd(igpt));
        });
    } else {
        aerosol_optics.alloc_2str(ncol, nlay, k_dist_sw.get_band_lims_wavenumber());
        yakl::c::parallel_for(yakl::c::Bounds<3>(nswbands,nlay,ncol), YAKL_LAMBDA (int ibnd, int ilay, int icol) {
            aerosol_optics_tau(icol,ilay,ibnd) = aer_tau_bnd(icol,ilay,ibnd);
            aerosol_optics_ssa(icol,ilay,ibnd) = aer_ssa_bnd(icol,ilay,ibnd);
            aerosol_optics_g  (icol,ilay,ibnd) = aer_asm_bnd(icol,ilay,ibnd);
        });
    }
    aerosol_optics.delta_scale();
    aerosol_optics.increment(combined_optics);

    // Do the clearsky calculation before adding in clouds
    FluxesByband fluxes_clrsky;
    fluxes_clrsky.flux_up = real2d("clrsky_flux_up", ncol, nlay+1); // clrsky_flux_up;
    fluxes_clrsky.flux_dn = real2d("clrsky_flux_up", ncol, nlay+1); //clrsky_flux_dn;
    fluxes_clrsky.flux_dn_dir = real2d("clrsky_flux_up", ncol, nlay+1); //clrsky_flux_dn_dir;
    fluxes_clrsky.flux_net = real2d("clrsky_flux_up", ncol, nlay+1); //clrsky_flux_net;
    fluxes_clrsky.bnd_flux_up = real3d("clrsky_flux_up", ncol, nlay+1, nswbands); //clrsky_bnd_flux_up;
    fluxes_clrsky.bnd_flux_dn = real3d("clrsky_flux_up", ncol, nlay+1, nswbands); //clrsky_bnd_flux_dn;
    fluxes_clrsky.bnd_flux_dn_dir = real3d("clrsky_flux_up", ncol, nlay+1, nswbands); //clrsky_bnd_flux_dn_dir;
    fluxes_clrsky.bnd_flux_net = real3d("clrsky_flux_up", ncol, nlay+1, nswbands); //clrsky_bnd_flux_net;
    rte_sw(combined_optics, top_at_1, coszrs, toa_flux, albedo_dir, albedo_dif, fluxes_clrsky);

    // Copy fluxes back out of FluxesByband object
    fluxes_clrsky.flux_up.deep_copy_to(clrsky_flux_up);
    fluxes_clrsky.flux_dn.deep_copy_to(clrsky_flux_dn);
    fluxes_clrsky.flux_dn_dir.deep_copy_to(clrsky_flux_dn_dir);
    fluxes_clrsky.flux_net.deep_copy_to(clrsky_flux_net);
    fluxes_clrsky.bnd_flux_up.deep_copy_to(clrsky_bnd_flux_up);
    fluxes_clrsky.bnd_flux_dn.deep_copy_to(clrsky_bnd_flux_dn);
    fluxes_clrsky.bnd_flux_dn_dir.deep_copy_to(clrsky_bnd_flux_dn_dir);
    fluxes_clrsky.bnd_flux_net.deep_copy_to(clrsky_bnd_flux_net);

    // Add in clouds
    OpticalProps2str cloud_optics;
    cloud_optics.alloc_2str(ncol, nlay, k_dist_sw);
    auto &cloud_optics_tau = cloud_optics.tau;
    auto &cloud_optics_ssa = cloud_optics.ssa;
    auto &cloud_optics_g   = cloud_optics.g  ;
    yakl::c::parallel_for(yakl::c::Bounds<3>(nswgpts,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
        cloud_optics_tau(icol,ilay,igpt) = cld_tau_gpt(icol,ilay,igpt);
        cloud_optics_ssa(icol,ilay,igpt) = cld_ssa_gpt(icol,ilay,igpt);
        cloud_optics_g  (icol,ilay,igpt) = cld_asm_gpt(icol,ilay,igpt);
    });
    cloud_optics.delta_scale();
    cloud_optics.increment(combined_optics);

    // Call SW flux driver
    FluxesByband fluxes_allsky;
    fluxes_allsky.flux_up = real2d("allsky_flux_up", ncol, nlay+1); // allsky_flux_up;
    fluxes_allsky.flux_dn = real2d("allsky_flux_up", ncol, nlay+1); //allsky_flux_dn;
    fluxes_allsky.flux_dn_dir = real2d("allsky_flux_up", ncol, nlay+1); //allsky_flux_dn_dir;
    fluxes_allsky.flux_net = real2d("allsky_flux_up", ncol, nlay+1); //allsky_flux_net;
    fluxes_allsky.bnd_flux_up = real3d("allsky_flux_up", ncol, nlay+1, nswbands); //allsky_bnd_flux_up;
    fluxes_allsky.bnd_flux_dn = real3d("allsky_flux_up", ncol, nlay+1, nswbands); //allsky_bnd_flux_dn;
    fluxes_allsky.bnd_flux_dn_dir = real3d("allsky_flux_up", ncol, nlay+1, nswbands); //allsky_bnd_flux_dn_dir;
    fluxes_allsky.bnd_flux_net = real3d("allsky_flux_up", ncol, nlay+1, nswbands); //allsky_bnd_flux_net;
    rte_sw(combined_optics, top_at_1, coszrs, toa_flux, albedo_dir, albedo_dif, fluxes_allsky);

    // Copy fluxes back out of FluxesByband object
    fluxes_allsky.flux_up.deep_copy_to(allsky_flux_up);
    fluxes_allsky.flux_dn.deep_copy_to(allsky_flux_dn);
    fluxes_allsky.flux_dn_dir.deep_copy_to(allsky_flux_dn_dir);
    fluxes_allsky.flux_net.deep_copy_to(allsky_flux_net);
    fluxes_allsky.bnd_flux_up.deep_copy_to(allsky_bnd_flux_up);
    fluxes_allsky.bnd_flux_dn.deep_copy_to(allsky_bnd_flux_dn);
    fluxes_allsky.bnd_flux_dn_dir.deep_copy_to(allsky_bnd_flux_dn_dir);
    fluxes_allsky.bnd_flux_net.deep_copy_to(allsky_bnd_flux_net);
    yakl::fence();
}


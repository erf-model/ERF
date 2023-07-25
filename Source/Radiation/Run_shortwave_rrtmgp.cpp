#include "Rrtmgp.H"

void Rrtmgp::run_shortwave_rrtmgp (
        int ngas, int ncol, int nlay,
        double *gas_vmr_p, double *pmid_p      , double *tmid_p      , double *pint_p,
        double *coszrs_p , double *albedo_dir_p, double *albedo_dif_p,
        double *cld_tau_gpt_p, double *cld_ssa_gpt_p, double *cld_asm_gpt_p,
        double *aer_tau_bnd_p, double *aer_ssa_bnd_p, double *aer_asm_bnd_p,
        double *allsky_flux_up_p    , double *allsky_flux_dn_p    , double *allsky_flux_net_p    , double *allsky_flux_dn_dir_p,
        double *allsky_bnd_flux_up_p, double *allsky_bnd_flux_dn_p, double *allsky_bnd_flux_net_p, double *allsky_bnd_flux_dn_dir_p,
        double *clrsky_flux_up_p    , double *clrsky_flux_dn_p    , double *clrsky_flux_net_p    , double *clrsky_flux_dn_dir_p,
        double *clrsky_bnd_flux_up_p, double *clrsky_bnd_flux_dn_p, double *clrsky_bnd_flux_net_p, double *clrsky_bnd_flux_dn_dir_p,
        double tsi_scaling) 
{

    // Wrap pointers in YAKL arrays
    int nswbands = k_dist_sw.get_nband();
    int nswgpts  = k_dist_sw.get_ngpt();
    auto gas_vmr_host                = realHost3d("gas_vmr", gas_vmr_p, ngas, ncol, nlay);
    auto pmid_host                   = realHost2d("pmid", pmid_p, ncol, nlay);
    auto tmid_host                   = realHost2d("tmid", tmid_p, ncol, nlay);
    auto pint_host                   = realHost2d("pint", pint_p, ncol, nlay+1);
    auto coszrs_host                 = realHost1d("coszrs", coszrs_p, ncol);
    auto albedo_dir_host             = realHost2d("albedo_dir", albedo_dir_p, nswbands, ncol);
    auto albedo_dif_host             = realHost2d("albedo_dif", albedo_dif_p, nswbands, ncol);
    auto cld_tau_gpt_host            = realHost3d("cld_tau_gpt", cld_tau_gpt_p, ncol, nlay, nswgpts);
    auto cld_ssa_gpt_host            = realHost3d("cld_ssa_gpt", cld_ssa_gpt_p, ncol, nlay, nswgpts);
    auto cld_asm_gpt_host            = realHost3d("cld_asm_gpt", cld_asm_gpt_p, ncol, nlay, nswgpts);
    auto aer_tau_bnd_host            = realHost3d("aer_tau_bnd", aer_tau_bnd_p, ncol, nlay, nswbands);
    auto aer_ssa_bnd_host            = realHost3d("aer_ssa_bnd", aer_ssa_bnd_p, ncol, nlay, nswbands);
    auto aer_asm_bnd_host            = realHost3d("aer_asm_bnd", aer_asm_bnd_p, ncol, nlay, nswbands);
    auto allsky_flux_up_host         = realHost2d("allsky_flux_up", allsky_flux_up_p, ncol, nlay+1);
    auto allsky_flux_dn_host         = realHost2d("allsky_flux_dn", allsky_flux_dn_p, ncol, nlay+1);
    auto allsky_flux_dn_dir_host     = realHost2d("allsky_flux_dn_dir", allsky_flux_dn_dir_p, ncol, nlay+1);
    auto allsky_flux_net_host        = realHost2d("allsky_flux_net", allsky_flux_net_p, ncol, nlay+1);
    auto clrsky_flux_up_host         = realHost2d("clrsky_flux_up", clrsky_flux_up_p, ncol, nlay+1);
    auto clrsky_flux_dn_host         = realHost2d("clrsky_flux_dn", clrsky_flux_dn_p, ncol, nlay+1);
    auto clrsky_flux_dn_dir_host     = realHost2d("clrsky_flux_dn_dir", clrsky_flux_dn_dir_p, ncol, nlay+1);
    auto clrsky_flux_net_host        = realHost2d("clrsky_flux_net", clrsky_flux_net_p, ncol, nlay+1);
    auto allsky_bnd_flux_up_host     = realHost3d("allsky_bnd_flux_up", allsky_bnd_flux_up_p, ncol, nlay+1, nswbands);
    auto allsky_bnd_flux_dn_host     = realHost3d("allsky_bnd_flux_dn", allsky_bnd_flux_dn_p, ncol, nlay+1, nswbands);
    auto allsky_bnd_flux_dn_dir_host = realHost3d("allsky_bnd_flux_dn_dir", allsky_bnd_flux_dn_dir_p, ncol, nlay+1, nswbands);
    auto allsky_bnd_flux_net_host    = realHost3d("allsky_bnd_flux_net", allsky_bnd_flux_net_p, ncol, nlay+1, nswbands);
    auto clrsky_bnd_flux_up_host     = realHost3d("clrsky_bnd_flux_up", clrsky_bnd_flux_up_p, ncol, nlay+1, nswbands);
    auto clrsky_bnd_flux_dn_host     = realHost3d("clrsky_bnd_flux_dn", clrsky_bnd_flux_dn_p, ncol, nlay+1, nswbands);
    auto clrsky_bnd_flux_dn_dir_host = realHost3d("clrsky_bnd_flux_dn_dir", clrsky_bnd_flux_dn_dir_p, ncol, nlay+1, nswbands);
    auto clrsky_bnd_flux_net_host    = realHost3d("clrsky_bnd_flux_net", clrsky_bnd_flux_net_p, ncol, nlay+1, nswbands);

    real3d gas_vmr               ("gas_vmr", ngas, ncol, nlay);
    real2d pmid                  ("pmid", ncol, nlay);
    real2d tmid                  ("tmid", ncol, nlay);
    real2d pint                  ("pint", ncol, nlay+1);
    real1d coszrs                ("coszrs", ncol);
    real2d albedo_dir            ("albedo_dir", nswbands, ncol);
    real2d albedo_dif            ("albedo_dif", nswbands, ncol);
    real3d cld_tau_gpt           ("cld_tau_gpt", ncol, nlay, nswgpts);
    real3d cld_ssa_gpt           ("cld_ssa_gpt", ncol, nlay, nswgpts);
    real3d cld_asm_gpt           ("cld_asm_gpt", ncol, nlay, nswgpts);
    real3d aer_tau_bnd           ("aer_tau_bnd", ncol, nlay, nswbands);
    real3d aer_ssa_bnd           ("aer_ssa_bnd", ncol, nlay, nswbands);
    real3d aer_asm_bnd           ("aer_asm_bnd", ncol, nlay, nswbands);
    real2d allsky_flux_up        ("allsky_flux_up", ncol, nlay+1);
    real2d allsky_flux_dn        ("allsky_flux_dn", ncol, nlay+1);
    real2d allsky_flux_dn_dir    ("allsky_flux_dn_dir", ncol, nlay+1);
    real2d allsky_flux_net       ("allsky_flux_net", ncol, nlay+1);
    real2d clrsky_flux_up        ("clrsky_flux_up", ncol, nlay+1);
    real2d clrsky_flux_dn        ("clrsky_flux_dn", ncol, nlay+1);
    real2d clrsky_flux_dn_dir    ("clrsky_flux_dn_dir", ncol, nlay+1);
    real2d clrsky_flux_net       ("clrsky_flux_net", ncol, nlay+1);
    real3d allsky_bnd_flux_up    ("allsky_bnd_flux_up", ncol, nlay+1, nswbands);
    real3d allsky_bnd_flux_dn    ("allsky_bnd_flux_dn", ncol, nlay+1, nswbands);
    real3d allsky_bnd_flux_dn_dir("allsky_bnd_flux_dn_dir", ncol, nlay+1, nswbands);
    real3d allsky_bnd_flux_net   ("allsky_bnd_flux_net", ncol, nlay+1, nswbands);
    real3d clrsky_bnd_flux_up    ("clrsky_bnd_flux_up", ncol, nlay+1, nswbands);
    real3d clrsky_bnd_flux_dn    ("clrsky_bnd_flux_dn", ncol, nlay+1, nswbands);
    real3d clrsky_bnd_flux_dn_dir("clrsky_bnd_flux_dn_dir", ncol, nlay+1, nswbands);
    real3d clrsky_bnd_flux_net   ("clrsky_bnd_flux_net", ncol, nlay+1, nswbands);

    // TODO: Only copy in the inputs
    gas_vmr_host               .deep_copy_to(gas_vmr               );
    pmid_host                  .deep_copy_to(pmid                  );
    tmid_host                  .deep_copy_to(tmid                  );
    pint_host                  .deep_copy_to(pint                  );
    coszrs_host                .deep_copy_to(coszrs                );
    albedo_dir_host            .deep_copy_to(albedo_dir            );
    albedo_dif_host            .deep_copy_to(albedo_dif            );
    cld_tau_gpt_host           .deep_copy_to(cld_tau_gpt           );
    cld_ssa_gpt_host           .deep_copy_to(cld_ssa_gpt           );
    cld_asm_gpt_host           .deep_copy_to(cld_asm_gpt           );
    aer_tau_bnd_host           .deep_copy_to(aer_tau_bnd           );
    aer_ssa_bnd_host           .deep_copy_to(aer_ssa_bnd           );
    aer_asm_bnd_host           .deep_copy_to(aer_asm_bnd           );
    //allsky_flux_up_host        .deep_copy_to(allsky_flux_up        );
    //allsky_flux_dn_host        .deep_copy_to(allsky_flux_dn        );
    //allsky_flux_dn_dir_host    .deep_copy_to(allsky_flux_dn_dir    );
    //allsky_flux_net_host       .deep_copy_to(allsky_flux_net       );
    //clrsky_flux_up_host        .deep_copy_to(clrsky_flux_up        );
    //clrsky_flux_dn_host        .deep_copy_to(clrsky_flux_dn        );
    //clrsky_flux_dn_dir_host    .deep_copy_to(clrsky_flux_dn_dir    );
    //clrsky_flux_net_host       .deep_copy_to(clrsky_flux_net       );
    //allsky_bnd_flux_up_host    .deep_copy_to(allsky_bnd_flux_up    );
    //allsky_bnd_flux_dn_host    .deep_copy_to(allsky_bnd_flux_dn    );
    //allsky_bnd_flux_dn_dir_host.deep_copy_to(allsky_bnd_flux_dn_dir);
    //allsky_bnd_flux_net_host   .deep_copy_to(allsky_bnd_flux_net   );
    //clrsky_bnd_flux_up_host    .deep_copy_to(clrsky_bnd_flux_up    );
    //clrsky_bnd_flux_dn_host    .deep_copy_to(clrsky_bnd_flux_dn    );
    //clrsky_bnd_flux_dn_dir_host.deep_copy_to(clrsky_bnd_flux_dn_dir);
    //clrsky_bnd_flux_net_host   .deep_copy_to(clrsky_bnd_flux_net   );


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
    bool top_at_1 = pmid_host(1, 1) < pmid_host (1, 2);
    real2d toa_flux("toa_flux", ncol, nswgpts);
    k_dist_sw.gas_optics(ncol, nlay, top_at_1, pmid, pint, tmid, gas_concs, combined_optics, toa_flux);

    // Apply TOA flux scaling
    yakl::c::parallel_for(yakl::c::Bounds<2>(nswgpts,ncol), YAKL_LAMBDA (int igpt, int icol) {
        toa_flux(icol, igpt) = tsi_scaling * toa_flux(icol, igpt);
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

    // TODO: Only copy out the outputs
    //gas_vmr               .deep_copy_to(gas_vmr_host               );
    //pmid                  .deep_copy_to(pmid_host                  );
    //tmid                  .deep_copy_to(tmid_host                  );
    //pint                  .deep_copy_to(pint_host                  );
    //coszrs                .deep_copy_to(coszrs_host                );
    //albedo_dir            .deep_copy_to(albedo_dir_host            );
    //albedo_dif            .deep_copy_to(albedo_dif_host            );
    //cld_tau_gpt           .deep_copy_to(cld_tau_gpt_host           );
    //cld_ssa_gpt           .deep_copy_to(cld_ssa_gpt_host           );
    //cld_asm_gpt           .deep_copy_to(cld_asm_gpt_host           );
    //aer_tau_bnd           .deep_copy_to(aer_tau_bnd_host           );
    //aer_ssa_bnd           .deep_copy_to(aer_ssa_bnd_host           );
    //aer_asm_bnd           .deep_copy_to(aer_asm_bnd_host           );
    allsky_flux_up        .deep_copy_to(allsky_flux_up_host        );
    allsky_flux_dn        .deep_copy_to(allsky_flux_dn_host        );
    allsky_flux_dn_dir    .deep_copy_to(allsky_flux_dn_dir_host    );
    allsky_flux_net       .deep_copy_to(allsky_flux_net_host       );
    clrsky_flux_up        .deep_copy_to(clrsky_flux_up_host        );
    clrsky_flux_dn        .deep_copy_to(clrsky_flux_dn_host        );
    clrsky_flux_dn_dir    .deep_copy_to(clrsky_flux_dn_dir_host    );
    clrsky_flux_net       .deep_copy_to(clrsky_flux_net_host       );
    allsky_bnd_flux_up    .deep_copy_to(allsky_bnd_flux_up_host    );
    allsky_bnd_flux_dn    .deep_copy_to(allsky_bnd_flux_dn_host    );
    allsky_bnd_flux_dn_dir.deep_copy_to(allsky_bnd_flux_dn_dir_host);
    allsky_bnd_flux_net   .deep_copy_to(allsky_bnd_flux_net_host   );
    clrsky_bnd_flux_up    .deep_copy_to(clrsky_bnd_flux_up_host    );
    clrsky_bnd_flux_dn    .deep_copy_to(clrsky_bnd_flux_dn_host    );
    clrsky_bnd_flux_dn_dir.deep_copy_to(clrsky_bnd_flux_dn_dir_host);
    clrsky_bnd_flux_net   .deep_copy_to(clrsky_bnd_flux_net_host   );
    yakl::fence();

} 


#include "ERF_RRTMGP_Interface.H"

namespace rrtmgp {

/*
 * Objects containing k-distribution information need to be initialized
 * once and then persist throughout the life of the program, so we
 * declare them here within the rrtmgp namespace.
 */
std::unique_ptr<gas_optics_t> k_dist_sw_k;
std::unique_ptr<gas_optics_t> k_dist_lw_k;

/*
 * Objects containing cloud optical property look-up table information.
 * We want to initialize these once and use throughout the life of the
 * program, so declare here and read data in during rrtmgp_initialize().
 */
std::unique_ptr<cloud_optics_t> cloud_optics_sw_k;
std::unique_ptr<cloud_optics_t> cloud_optics_lw_k;


pool_t kokkos_mem_pool;

bool initialized = false;


optical_props2_t
get_cloud_optics_sw (const int ncol,
                     const int nlay,
                     cloud_optics_t& cloud_optics,
                     gas_optics_t& kdist,
                     real2d_k& lwp,
                     real2d_k& iwp,
                     real2d_k& rel,
                     real2d_k& rei)
{
    // Initialize optics
    optical_props2_t clouds;
    clouds.init(kdist.get_band_lims_wavenumber());
    clouds.alloc_2str(ncol, nlay);

    // Needed for consistency with all-sky example problem?
    cloud_optics.set_ice_roughness(2);

    // Limit effective radii to be within bounds of lookup table
    real2d_k rel_limited("rel_limited", ncol, nlay);
    real2d_k rei_limited("rei_limited", ncol, nlay);
    limit_to_bounds_2d(rel, cloud_optics.radliq_lwr,
                       cloud_optics.radliq_upr, rel_limited);
    limit_to_bounds_2d(rei, cloud_optics.radice_lwr,
                       cloud_optics.radice_upr, rei_limited);

    // Calculate cloud optics
    cloud_optics.cloud_optics(ncol, nlay, lwp, iwp, rel_limited, rei_limited, clouds);

    // Return optics
    return clouds;
}


optical_props1_t
get_cloud_optics_lw (const int ncol,
                     const int nlay,
                     cloud_optics_t& cloud_optics,
                     gas_optics_t& kdist,
                     real2d_k& lwp,
                     real2d_k& iwp,
                     real2d_k& rel,
                     real2d_k& rei)
{
    // Initialize optics
    optical_props1_t clouds;
    clouds.init(kdist.get_band_lims_wavenumber());
    clouds.alloc_1scl(ncol, nlay);

    // Needed for consistency with all-sky example problem?
    cloud_optics.set_ice_roughness(2);

    // Limit effective radii to be within bounds of lookup table
    real2d_k rel_limited("rel_limited", ncol, nlay);
    real2d_k rei_limited("rei_limited", ncol, nlay);
    limit_to_bounds_2d(rel, cloud_optics.radliq_lwr,
                       cloud_optics.radliq_upr, rel_limited);
    limit_to_bounds_2d(rei, cloud_optics.radice_lwr,
                       cloud_optics.radice_upr, rei_limited);

    // Calculate cloud optics
    cloud_optics.cloud_optics(ncol, nlay, lwp, iwp, rel_limited, rei_limited, clouds);

    // Return optics
    return clouds;
}


optical_props2_t
get_subsampled_clouds (const int ncol,
                       const int nlay,
                       const int nbnd,
                       const int ngpt,
                       optical_props2_t& cloud_optics,
                       gas_optics_t& kdist,
                       real2d_k& cld,
                       real2d_k& p_lay)
{
    // Initialized subsampled optics
    optical_props2_t subsampled_optics;
    subsampled_optics.init(kdist.get_band_lims_wavenumber(), kdist.get_band_lims_gpoint(), "subsampled_optics");
    subsampled_optics.alloc_2str(ncol, nlay);

    // Check that we do not have clouds with no optical properties; this would get corrected
    // when we assign optical props, but we want to use a "radiative cloud fraction"
    // for the subcolumn sampling too because otherwise we can get vertically-contiguous cloud
    // mask profiles with no actual cloud properties in between, which would just further overestimate
    // the vertical correlation of cloudy layers. I.e., cloudy layers might look maximally overlapped
    // even when separated by layers with no cloud properties, when in fact those layers should be
    // randomly overlapped.
    real2d_k cldfrac_rad("cldfrac_rad", ncol, nlay);
    Kokkos::deep_copy(cldfrac_rad, 0.0);
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay, nbnd}),
                         KOKKOS_LAMBDA (int icol, int ilay, int ibnd)
    {
        if (cloud_optics.tau(icol,ilay,ibnd) > 0) {
            cldfrac_rad(icol,ilay) = cld(icol,ilay);
        }
    });

    // Get subcolumn cloud mask; note that get_subcolumn_mask exposes overlap assumption as an option,
    // but the only currently supported options are 0 (trivial all-or-nothing cloud) or 1 (max-rand),
    // so overlap has not been exposed as an option beyond this subcolumn. In the future, we should
    // support generalized overlap as well, with parameters derived from DPSCREAM simulations with very
    // high resolution.
    int overlap = 1;

    // Get unique seeds for each column that are reproducible across different MPI rank layouts;
    // use decimal part of pressure for this, consistent with the implementation in EAM
    int1d_k seeds("seeds", ncol);
    Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol)
    {
        seeds(icol) = 1.0e9 * (p_lay(icol,nlay-1) - int(p_lay(icol,nlay-1)));
    });
    auto cldmask = get_subcolumn_mask(ncol, nlay, ngpt, cldfrac_rad, overlap, seeds);

    // Assign optical properties to subcolumns (note this implements MCICA)
    auto gpoint_bands = kdist.get_gpoint_bands();
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay, ngpt}),
                         KOKKOS_LAMBDA (int icol, int ilay, int igpt)
    {
        auto ibnd = gpoint_bands(igpt);
        if (cldmask(icol,ilay,igpt) == 1) {
            subsampled_optics.tau(icol,ilay,igpt) = cloud_optics.tau(icol,ilay,ibnd);
            subsampled_optics.ssa(icol,ilay,igpt) = cloud_optics.ssa(icol,ilay,ibnd);
            subsampled_optics.g  (icol,ilay,igpt) = cloud_optics.g  (icol,ilay,ibnd);
        } else {
            subsampled_optics.tau(icol,ilay,igpt) = 0;
            subsampled_optics.ssa(icol,ilay,igpt) = 0;
            subsampled_optics.g  (icol,ilay,igpt) = 0;
        }
    });
    return subsampled_optics;
}


optical_props1_t
get_subsampled_clouds (const int ncol,
                       const int nlay,
                       const int nbnd,
                       const int ngpt,
                       optical_props1_t& cloud_optics,
                       gas_optics_t& kdist,
                       real2d_k& cld,
                       real2d_k& p_lay)
{
    // Initialized subsampled optics
    optical_props1_t subsampled_optics;
    subsampled_optics.init(kdist.get_band_lims_wavenumber(), kdist.get_band_lims_gpoint(), "subsampled_optics");
    subsampled_optics.alloc_1scl(ncol, nlay);

    // Check that we do not have clouds with no optical properties; this would get corrected
    // when we assign optical props, but we want to use a "radiative cloud fraction"
    // for the subcolumn sampling too because otherwise we can get vertically-contiguous cloud
    // mask profiles with no actual cloud properties in between, which would just further overestimate
    // the vertical correlation of cloudy layers. I.e., cloudy layers might look maximally overlapped
    // even when separated by layers with no cloud properties, when in fact those layers should be
    // randomly overlapped.
    real2d_k cldfrac_rad("cldfrac_rad", ncol, nlay);
    Kokkos::deep_copy(cldfrac_rad, 0.0);
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay, nbnd}),
                         KOKKOS_LAMBDA (int icol, int ilay, int ibnd)
    {
        if (cloud_optics.tau(icol,ilay,ibnd) > 0) {
            cldfrac_rad(icol,ilay) = cld(icol,ilay);
        }
    });

    // Get subcolumn cloud mask
    int overlap = 1;
    // Get unique seeds for each column that are reproducible across different MPI rank layouts;
    // use decimal part of pressure for this, consistent with the implementation in EAM; use different
    // seed values for longwave and shortwave
    int1d_k seeds("seeds", ncol);
    Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol)
    {
        seeds(icol) = 1e9 * (p_lay(icol,nlay-1) - int(p_lay(icol,nlay-1)));
    });
    auto cldmask = get_subcolumn_mask(ncol, nlay, ngpt, cldfrac_rad, overlap, seeds);

    // Assign optical properties to subcolumns (note this implements MCICA)
    auto gpoint_bands = kdist.get_gpoint_bands();
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay, ngpt}),
                         KOKKOS_LAMBDA (int icol, int ilay, int igpt)
    {
        auto ibnd = gpoint_bands(igpt);
        if (cldmask(icol,ilay,igpt) == 1) {
            subsampled_optics.tau(icol,ilay,igpt) = cloud_optics.tau(icol,ilay,ibnd);
        } else {
            subsampled_optics.tau(icol,ilay,igpt) = 0;
        }
    });
    return subsampled_optics;
}


/*
 * The following routines provide a simple interface to RRTMGP. These
 * can be used as-is, but are intended to be wrapped by the SCREAM AD
 * interface to radiation.
 */
void
rrtmgp_initialize (gas_concs_t& gas_concs_k,
                   const std::string& coefficients_file_sw,
                   const std::string& coefficients_file_lw,
                   const std::string& cloud_optics_file_sw,
                   const std::string& cloud_optics_file_lw)
{
    // Initialize Kokkos
    if (!Kokkos::is_initialized()) {  Kokkos::initialize(); }

    // Create objects for static ptrs
    k_dist_sw_k = std::make_unique<gas_optics_t>();
    k_dist_lw_k = std::make_unique<gas_optics_t>();
    cloud_optics_sw_k = std::make_unique<cloud_optics_t>();
    cloud_optics_lw_k = std::make_unique<cloud_optics_t>();

    // Load and initialize absorption coefficient data
    load_and_init(*k_dist_sw_k, coefficients_file_sw, gas_concs_k);
    load_and_init(*k_dist_lw_k, coefficients_file_lw, gas_concs_k);

    // Load and initialize cloud optical property look-up table information
    load_cld_lutcoeff(*cloud_optics_sw_k, cloud_optics_file_sw);
    load_cld_lutcoeff(*cloud_optics_lw_k, cloud_optics_file_lw);

    // Initialize kokkos rrtmgp pool allocator
    const size_t nvar = 300;
    const size_t nbnd = std::max(k_dist_sw_k->get_nband(),k_dist_sw_k->get_nband());
    const size_t ncol = gas_concs_k.ncol;
    const size_t nlay = gas_concs_k.nlay;
    auto my_size_ref = static_cast<unsigned long>(nvar * ncol * nlay * nbnd);
    pool_t::init(my_size_ref);

    // We are now initialized!
    initialized = true;
}


void
rrtmgp_finalize ()
{
    initialized = false;
    k_dist_sw_k->finalize();
    k_dist_lw_k->finalize();
    cloud_optics_sw_k->finalize();
    cloud_optics_lw_k->finalize();
    k_dist_sw_k = nullptr;
    k_dist_lw_k = nullptr;
    cloud_optics_sw_k = nullptr;
    cloud_optics_lw_k = nullptr;
    pool_t::finalize();
}


void
compute_band_by_band_surface_albedos (const int ncol,
                                      const int nswbands,
                                      real1d_k& sfc_alb_dir_vis,
                                      real1d_k& sfc_alb_dir_nir,
                                      real1d_k& sfc_alb_dif_vis,
                                      real1d_k& sfc_alb_dif_nir,
                                      real2d_k& sfc_alb_dir,
                                      real2d_k& sfc_alb_dif)
{
    auto wavenumber_limits = k_dist_sw_k->get_band_lims_wavenumber();

    // Loop over bands, and determine for each band whether it is broadly in the
    // visible or infrared part of the spectrum (visible or "not visible")
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nswbands}),
                         KOKKOS_LAMBDA (int icol, int ibnd)
    {
        // Threshold between visible and infrared is 0.7 micron, or 14286 cm^-1.
        const RealT visible_wavenumber_threshold = 14286.0;

        // Wavenumber is in the visible if it is above the visible wavenumber
        // threshold, and in the infrared if it is below the threshold
        const bool is_visible_wave1 = (wavenumber_limits(0, ibnd) > visible_wavenumber_threshold ? true : false);
        const bool is_visible_wave2 = (wavenumber_limits(1, ibnd) > visible_wavenumber_threshold ? true : false);

        if (is_visible_wave1 && is_visible_wave2) {
            // Entire band is in the visible
            sfc_alb_dir(icol,ibnd) = sfc_alb_dir_vis(icol);
            sfc_alb_dif(icol,ibnd) = sfc_alb_dif_vis(icol);
        } else if (!is_visible_wave1 && !is_visible_wave2) {
            // Entire band is in the longwave (near-infrared)
            sfc_alb_dir(icol,ibnd) = sfc_alb_dir_nir(icol);
            sfc_alb_dif(icol,ibnd) = sfc_alb_dif_nir(icol);
        } else {
            // Band straddles the visible to near-infrared transition, so we take
            // the albedo to be the average of the visible and near-infrared
            // broadband albedos
            sfc_alb_dir(icol,ibnd) = 0.5*(sfc_alb_dir_vis(icol) + sfc_alb_dir_nir(icol));
            sfc_alb_dif(icol,ibnd) = 0.5*(sfc_alb_dif_vis(icol) + sfc_alb_dif_nir(icol));
        }
    });
}


void
compute_broadband_surface_fluxes (const int ncol,
                                  const int kbot,
                                  const int nswbands,
                                  real3d_k& sw_bnd_flux_dir ,
                                  real3d_k& sw_bnd_flux_dif ,
                                  real1d_k& sfc_flux_dir_vis,
                                  real1d_k& sfc_flux_dir_nir,
                                  real1d_k& sfc_flux_dif_vis,
                                  real1d_k& sfc_flux_dif_nir)
{
    // Band 10 straddles the near-IR and visible, so divide contributions from band 10 between both broadband sums
    // TODO: Hard-coding these band indices is really bad practice. If the bands ever were to change (like when
    // the RRTMG bands were re-ordered for RRTMGP), we would be using the wrong bands for the IR and UV/VIS. This
    // should be refactored to grab the correct bands by specifying appropriate wavenumber rather than index.
    //sfc_flux_dir_nir(i) = sum(sw_bnd_flux_dir(i+1,kbot,1:9))   + 0.5 * sw_bnd_flux_dir(i+1,kbot,10);
    //sfc_flux_dir_vis(i) = sum(sw_bnd_flux_dir(i+1,kbot,11:14)) + 0.5 * sw_bnd_flux_dir(i+1,kbot,10);
    //sfc_flux_dif_nir(i) = sum(sw_bnd_flux_dif(i+1,kbot,1:9))   + 0.5 * sw_bnd_flux_dif(i+1,kbot,10);
    //sfc_flux_dif_vis(i) = sum(sw_bnd_flux_dif(i+1,kbot,11:14)) + 0.5 * sw_bnd_flux_dif(i+1,kbot,10);

    // Initialize sums over bands
    Kokkos::deep_copy(sfc_flux_dir_nir, 0);
    Kokkos::deep_copy(sfc_flux_dir_vis, 0);
    Kokkos::deep_copy(sfc_flux_dif_nir, 0);
    Kokkos::deep_copy(sfc_flux_dif_vis, 0);

    // Threshold between visible and infrared is 0.7 micron, or 14286 cm^-1.
    const RealT visible_wavenumber_threshold = 14286.0;
    auto wavenumber_limits = k_dist_sw_k->get_band_lims_wavenumber();
    Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(const int icol)
    {
        for (int ibnd = 0; ibnd < nswbands; ++ibnd) {
            // Wavenumber is in the visible if it is above the visible wavenumber
            // threshold, and in the infrared if it is below the threshold
            const bool is_visible_wave1 = (wavenumber_limits(0, ibnd) > visible_wavenumber_threshold ? true : false);
            const bool is_visible_wave2 = (wavenumber_limits(1, ibnd) > visible_wavenumber_threshold ? true : false);

            if (is_visible_wave1 && is_visible_wave2) {
                // Entire band is in the visible
                sfc_flux_dir_vis(icol) += sw_bnd_flux_dir(icol,kbot,ibnd);
                sfc_flux_dif_vis(icol) += sw_bnd_flux_dif(icol,kbot,ibnd);
            } else if (!is_visible_wave1 && !is_visible_wave2) {
                // Entire band is in the longwave (near-infrared)
                sfc_flux_dir_nir(icol) += sw_bnd_flux_dir(icol,kbot,ibnd);
                sfc_flux_dif_nir(icol) += sw_bnd_flux_dif(icol,kbot,ibnd);
            } else {
                // Band straddles the visible to near-infrared transition, so put half
                // the flux in visible and half in near-infrared fluxes
                sfc_flux_dir_vis(icol) += 0.5 * sw_bnd_flux_dir(icol,kbot,ibnd);
                sfc_flux_dif_vis(icol) += 0.5 * sw_bnd_flux_dif(icol,kbot,ibnd);
                sfc_flux_dir_nir(icol) += 0.5 * sw_bnd_flux_dir(icol,kbot,ibnd);
                sfc_flux_dif_nir(icol) += 0.5 * sw_bnd_flux_dif(icol,kbot,ibnd);
            }
        }
    });
}


void
rrtmgp_main (const int ncol, const int nlay,
             real2d_k& p_lay, real2d_k& t_lay,
             real2d_k& p_lev, real2d_k& t_lev,
             gas_concs_t& gas_concs,
             real2d_k& sfc_alb_dir, real2d_k& sfc_alb_dif, real1d_k& mu0,
             real1d_k& t_sfc      , real1d_k& sfc_emis   , real1d_k& lw_src,
             real2d_k& lwp, real2d_k& iwp,
             real2d_k& rel, real2d_k& rei, real2d_k& cldfrac,
             real3d_k& /*aer_tau_sw*/, real3d_k& /*aer_ssa_sw*/, real3d_k& /*aer_asm_sw*/,
             real3d_k& /*aer_tau_lw*/,
             real3d_k& /*cld_tau_sw_bnd*/, real3d_k& /*cld_tau_lw_bnd*/,
             real3d_k& /*cld_tau_sw_gpt*/, real3d_k& /*cld_tau_lw_gpt*/,
             real2d_k& sw_flux_up, real2d_k& sw_flux_dn, real2d_k& sw_flux_dn_dir,
             real2d_k& lw_flux_up, real2d_k& lw_flux_dn,
             real2d_k& sw_clnclrsky_flux_up, real2d_k& sw_clnclrsky_flux_dn, real2d_k& sw_clnclrsky_flux_dn_dir,
             real2d_k& sw_clrsky_flux_up   , real2d_k& sw_clrsky_flux_dn   , real2d_k& sw_clrsky_flux_dn_dir,
             real2d_k& sw_clnsky_flux_up   , real2d_k& sw_clnsky_flux_dn   , real2d_k& sw_clnsky_flux_dn_dir,
             real2d_k& lw_clnclrsky_flux_up, real2d_k& lw_clnclrsky_flux_dn,
             real2d_k& lw_clrsky_flux_up   , real2d_k& lw_clrsky_flux_dn   ,
             real2d_k& lw_clnsky_flux_up   , real2d_k& lw_clnsky_flux_dn   ,
             real3d_k& sw_bnd_flux_up      , real3d_k& sw_bnd_flux_dn      , real3d_k& sw_bnd_flux_dn_dir,
             real3d_k& lw_bnd_flux_up      , real3d_k& lw_bnd_flux_dn      ,
             const RealT tsi_scaling,
             const bool extra_clnclrsky_diag, const bool extra_clnsky_diag)
{
    // Setup pointers to RRTMGP SW fluxes
    fluxes_t fluxes_sw;
    fluxes_sw.flux_up = sw_flux_up;
    fluxes_sw.flux_dn = sw_flux_dn;
    fluxes_sw.flux_dn_dir = sw_flux_dn_dir;
    fluxes_sw.bnd_flux_up = sw_bnd_flux_up;
    fluxes_sw.bnd_flux_dn = sw_bnd_flux_dn;
    fluxes_sw.bnd_flux_dn_dir = sw_bnd_flux_dn_dir;
    // Clean-clear-sky
    fluxes_broadband_t clnclrsky_fluxes_sw;
    clnclrsky_fluxes_sw.flux_up = sw_clnclrsky_flux_up;
    clnclrsky_fluxes_sw.flux_dn = sw_clnclrsky_flux_dn;
    clnclrsky_fluxes_sw.flux_dn_dir = sw_clnclrsky_flux_dn_dir;
    // Clear-sky
    fluxes_broadband_t clrsky_fluxes_sw;
    clrsky_fluxes_sw.flux_up = sw_clrsky_flux_up;
    clrsky_fluxes_sw.flux_dn = sw_clrsky_flux_dn;
    clrsky_fluxes_sw.flux_dn_dir = sw_clrsky_flux_dn_dir;
    // Clean-sky
    fluxes_broadband_t clnsky_fluxes_sw;
    clnsky_fluxes_sw.flux_up = sw_clnsky_flux_up;
    clnsky_fluxes_sw.flux_dn = sw_clnsky_flux_dn;
    clnsky_fluxes_sw.flux_dn_dir = sw_clnsky_flux_dn_dir;

    // Setup pointers to RRTMGP LW fluxes
    fluxes_t fluxes_lw;
    fluxes_lw.flux_up = lw_flux_up;
    fluxes_lw.flux_dn = lw_flux_dn;
    fluxes_lw.bnd_flux_up = lw_bnd_flux_up;
    fluxes_lw.bnd_flux_dn = lw_bnd_flux_dn;
    // Clean-clear-sky
    fluxes_broadband_t clnclrsky_fluxes_lw;
    clnclrsky_fluxes_lw.flux_up = lw_clnclrsky_flux_up;
    clnclrsky_fluxes_lw.flux_dn = lw_clnclrsky_flux_dn;
    // Clear-sky
    fluxes_broadband_t clrsky_fluxes_lw;
    clrsky_fluxes_lw.flux_up = lw_clrsky_flux_up;
    clrsky_fluxes_lw.flux_dn = lw_clrsky_flux_dn;
    // Clean-sky
    fluxes_broadband_t clnsky_fluxes_lw;
    clnsky_fluxes_lw.flux_up = lw_clnsky_flux_up;
    clnsky_fluxes_lw.flux_dn = lw_clnsky_flux_dn;

    auto nswbands = k_dist_sw_k->get_nband();
    auto nlwbands = k_dist_lw_k->get_nband();

    // Setup aerosol optical properties
    optical_props2_t aerosol_sw;
    optical_props1_t aerosol_lw;
    aerosol_sw.init(k_dist_sw_k->get_band_lims_wavenumber());
    aerosol_sw.alloc_2str(ncol, nlay);
    /*
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay, nswbands}),
                         KOKKOS_LAMBDA (int icol, int ilay, int ibnd)
    {
        aerosol_sw.tau(icol,ilay,ibnd) = aer_tau_sw(icol,ilay,ibnd);
        aerosol_sw.ssa(icol,ilay,ibnd) = aer_ssa_sw(icol,ilay,ibnd);
        aerosol_sw.g  (icol,ilay,ibnd) = aer_asm_sw(icol,ilay,ibnd);
    });
    */
    aerosol_lw.init(k_dist_lw_k->get_band_lims_wavenumber());
    aerosol_lw.alloc_1scl(ncol, nlay);
    /*
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay, nlwbands}),
                         KOKKOS_LAMBDA (int icol, int ilay, int ibnd)
    {
        aerosol_lw.tau(icol,ilay,ibnd) = aer_tau_lw(icol,ilay,ibnd);
    });
    */

    // Convert cloud physical properties to optical properties for input to RRTMGP
    optical_props2_t clouds_sw = get_cloud_optics_sw(ncol, nlay,
                                                     *cloud_optics_sw_k, *k_dist_sw_k,
                                                     lwp, iwp, rel, rei);
    optical_props1_t clouds_lw = get_cloud_optics_lw(ncol, nlay,
                                                     *cloud_optics_lw_k, *k_dist_lw_k,
                                                     lwp, iwp, rel, rei);
    //Kokkos::deep_copy(cld_tau_sw_bnd, clouds_sw.tau);
    //Kokkos::deep_copy(cld_tau_lw_bnd, clouds_lw.tau);

    // Do subcolumn sampling to map bands -> gpoints based on cloud fraction and overlap assumption;
    // This implements the Monte Carlo Independing Column Approximation by mapping only a single
    // subcolumn (cloud state) to each gpoint.
    auto nswgpts = k_dist_sw_k->get_ngpt();
    auto clouds_sw_gpt = get_subsampled_clouds(ncol, nlay, nswbands, nswgpts,
                                               clouds_sw, *k_dist_sw_k, cldfrac, p_lay);
    // Longwave
    auto nlwgpts = k_dist_lw_k->get_ngpt();
    auto clouds_lw_gpt = get_subsampled_clouds(ncol, nlay, nlwbands, nlwgpts,
                                               clouds_lw, *k_dist_lw_k, cldfrac, p_lay);

    /*
    // Copy cloud properties to outputs (is this needed, or can we just use pointers?)
    // Alternatively, just compute and output a subcolumn cloud mask
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay, nswgpts}),
                         KOKKOS_LAMBDA (int icol, int ilay, int igpt)
    {
        cld_tau_sw_gpt(icol,ilay,igpt) = clouds_sw_gpt.tau(icol,ilay,igpt);
    });
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay, nlwgpts}),
                         KOKKOS_LAMBDA (int icol, int ilay, int igpt)
    {
        cld_tau_lw_gpt(icol,ilay,igpt) = clouds_lw_gpt.tau(icol,ilay,igpt);
    });
    */

  // Do shortwave
  rrtmgp_sw(ncol, nlay,
            *k_dist_sw_k, p_lay, t_lay, p_lev, t_lev, gas_concs,
            sfc_alb_dir, sfc_alb_dif, mu0, aerosol_sw, clouds_sw_gpt,
            fluxes_sw, clnclrsky_fluxes_sw, clrsky_fluxes_sw, clnsky_fluxes_sw,
            tsi_scaling, extra_clnclrsky_diag, extra_clnsky_diag);

  // Do longwave
  rrtmgp_lw(ncol, nlay,
            *k_dist_lw_k, p_lay, t_lay, p_lev, t_lev,
            t_sfc, sfc_emis, lw_src,
            gas_concs, aerosol_lw, clouds_lw_gpt,
            fluxes_lw, clnclrsky_fluxes_lw, clrsky_fluxes_lw, clnsky_fluxes_lw,
            extra_clnclrsky_diag, extra_clnsky_diag);

}


int3d_k
get_subcolumn_mask (const int ncol,
                    const int nlay,
                    const int ngpt,
                    real2d_k& cldf,
                    const int overlap_option,
                    int1d_k& seeds)
{
    // Routine will return subcolumn mask with values of 0 indicating no cloud, 1 indicating cloud
    int3d_k subcolumn_mask("subcolumn_mask", ncol, nlay, ngpt);

    // Subcolumn generators are a means for producing a variable x(i,j,k), where
    //
    //     c(i,j,k) = 1 for x(i,j,k) >  1 - cldf(i,j)
    //     c(i,j,k) = 0 for x(i,j,k) <= 1 - cldf(i,j)
    //
    // I am going to call this "cldx" to be just slightly less ambiguous
    real3d_k cldx("cldx", ncol, nlay, ngpt);

    // Apply overlap assumption to set cldx
    if (overlap_option == 0) {  // Dummy mask, always cloudy
        Kokkos::deep_copy(cldx, 1);
    } else {  // Default case, maximum-random overlap
        // Maximum-random overlap:
        // Uses essentially the algorithm described in eq (14) in Raisanen et al. 2004,
        // https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1256/qj.03.99. Also the same
        // algorithm used in RRTMG implementation of maximum-random overlap (see
        // https://github.com/AER-RC/RRTMG_SW/blob/master/src/mcica_subcol_gen_sw.f90)
        //
        // First, fill cldx with random numbers.
        Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay, ngpt}),
                             KOKKOS_LAMBDA (int icol, int ilay, int igpt)
        {
            conv::Random rand(seeds(icol) + ilay*ngpt + igpt);
            cldx(icol,ilay,igpt) = rand.genFP<RealT>();
        });

        // Step down columns and apply algorithm from eq (14)
        Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, ngpt}),
                             KOKKOS_LAMBDA (int icol, int igpt)
        {
            for (int ilay = 1; ilay < nlay; ilay++) {
                // Check cldx in level above and see if it satisfies conditions to create a cloudy subcolumn
                if (cldx(icol,ilay-1,igpt) > 1.0 - cldf(icol,ilay-1)) {
                    // Cloudy subcolumn above, use same random number here so that clouds in these two adjacent
                    // layers are maximimally overlapped
                    cldx(icol,ilay,igpt) = cldx(icol,ilay-1,igpt);
                } else {
                    // Cloud-less above, use new random number so that clouds are distributed
                    // randomly in this layer. Need to scale new random number to range
                    // [0, 1.0 - cldf(ilay-1)] because we have artificially changed the distribution
                    // of random numbers in this layer with the above branch of the conditional,
                    // which would otherwise inflate cloud fraction in this layer.
                    cldx(icol,ilay,igpt) = cldx(icol,ilay  ,igpt) * (1.0 - cldf(icol,ilay-1));
                }
            }
        });
    }

    // Use cldx array to create subcolumn mask
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay, ngpt}),
                             KOKKOS_LAMBDA (int icol, int ilay, int igpt)
    {
        if (cldx(icol,ilay,igpt) > 1.0 - cldf(icol,ilay)) {
            subcolumn_mask(icol,ilay,igpt) = 1;
        } else {
            subcolumn_mask(icol,ilay,igpt) = 0;
        }
    });
    return subcolumn_mask;
}


void
rrtmgp_sw (const int ncol,
           const int nlay,
           gas_optics_t& k_dist,
           real2d_k& p_lay, real2d_k& t_lay,
           real2d_k& p_lev, real2d_k& t_lev,
           gas_concs_t& gas_concs,
           real2d_k& sfc_alb_dir, real2d_k& sfc_alb_dif,
           real1d_k& mu0,
           optical_props2_t& aerosol,
           optical_props2_t& clouds,
           fluxes_t& fluxes,
           fluxes_broadband_t& clnclrsky_fluxes,
           fluxes_broadband_t& clrsky_fluxes,
           fluxes_broadband_t& clnsky_fluxes,
           const RealT tsi_scaling,
           const bool extra_clnclrsky_diag,
           const bool extra_clnsky_diag)
{
    // Get problem sizes
    int nbnd = k_dist.get_nband();
    int ngpt = k_dist.get_ngpt();
    int ngas = gas_concs.get_num_gases();

    // Associate local pointers for fluxes
    auto& flux_up = fluxes.flux_up;
    auto& flux_dn = fluxes.flux_dn;
    auto& flux_dn_dir = fluxes.flux_dn_dir;
    auto& bnd_flux_up = fluxes.bnd_flux_up;
    auto& bnd_flux_dn = fluxes.bnd_flux_dn;
    auto& bnd_flux_dn_dir = fluxes.bnd_flux_dn_dir;
    auto& clnclrsky_flux_up = clnclrsky_fluxes.flux_up;
    auto& clnclrsky_flux_dn = clnclrsky_fluxes.flux_dn;
    auto& clnclrsky_flux_dn_dir = clnclrsky_fluxes.flux_dn_dir;
    auto& clrsky_flux_up = clrsky_fluxes.flux_up;
    auto& clrsky_flux_dn = clrsky_fluxes.flux_dn;
    auto& clrsky_flux_dn_dir = clrsky_fluxes.flux_dn_dir;
    auto& clnsky_flux_up = clnsky_fluxes.flux_up;
    auto& clnsky_flux_dn = clnsky_fluxes.flux_dn;
    auto& clnsky_flux_dn_dir = clnsky_fluxes.flux_dn_dir;

    // Reset fluxes to zero
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlay+1}),
                             KOKKOS_LAMBDA (int icol, int ilev)
    {
        flux_up    (icol,ilev) = 0;
        flux_dn    (icol,ilev) = 0;
        flux_dn_dir(icol,ilev) = 0;
        clnclrsky_flux_up    (icol,ilev) = 0;
        clnclrsky_flux_dn    (icol,ilev) = 0;
        clnclrsky_flux_dn_dir(icol,ilev) = 0;
        clrsky_flux_up    (icol,ilev) = 0;
        clrsky_flux_dn    (icol,ilev) = 0;
        clrsky_flux_dn_dir(icol,ilev) = 0;
        clnsky_flux_up    (icol,ilev) = 0;
        clnsky_flux_dn    (icol,ilev) = 0;
        clnsky_flux_dn_dir(icol,ilev) = 0;
    });
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay+1, nbnd}),
                         KOKKOS_LAMBDA (int icol, int ilev, int ibnd)
    {
        bnd_flux_up    (icol,ilev,ibnd) = 0;
        bnd_flux_dn    (icol,ilev,ibnd) = 0;
        bnd_flux_dn_dir(icol,ilev,ibnd) = 0;
    });

    // Get daytime indices
    int1d_k dayIndices("dayIndices", ncol);
    Kokkos::deep_copy(dayIndices, -1);

    // Serialized for now.
    int nday = 0;
    Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(int, int& nday_inner)
    {
        for (int icol = 0; icol < ncol; ++icol) {
            if (mu0(icol) > 0) {
                dayIndices(nday_inner++) = icol;
            }
        }
    }, Kokkos::Sum<int>(nday));

    // Copy data back to the device
    if (nday == 0) {
        // No daytime columns in this chunk, skip the rest of this routine
        return;
    }

    // Subset mu0
    real1d_k mu0_day("mu0_day", nday);
    Kokkos::parallel_for(nday, KOKKOS_LAMBDA(int iday)
    {
        mu0_day(iday) = mu0(dayIndices(iday));
    });

    // subset state variables
    real2d_k p_lay_day("p_lay_day", nday, nlay);
    real2d_k t_lay_day("t_lay_day", nday, nlay);
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nday, nlay}),
                         KOKKOS_LAMBDA (int iday, int ilay)
    {
        p_lay_day(iday,ilay) = p_lay(dayIndices(iday),ilay);
        t_lay_day(iday,ilay) = t_lay(dayIndices(iday),ilay);
    });
    real2d_k p_lev_day("p_lev_day", nday, nlay+1);
    real2d_k t_lev_day("t_lev_day", nday, nlay+1);
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nday, nlay+1}),
                         KOKKOS_LAMBDA (int iday, int ilay)
    {
        p_lev_day(iday,ilay) = p_lev(dayIndices(iday),ilay);
        t_lev_day(iday,ilay) = t_lev(dayIndices(iday),ilay);
    });

    // Subset gases
    auto gas_names = gas_concs.get_gas_names();
    gas_concs_t gas_concs_day;
    gas_concs_day.init(gas_names, nday, nlay);
    for (int igas = 0; igas < ngas; igas++) {
        real2d_k vmr_day("vmr_day", nday, nlay);
        real2d_k vmr("vmr"    , ncol, nlay);
        gas_concs.get_vmr(gas_names[igas], vmr);
        Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nday, nlay}),
                         KOKKOS_LAMBDA (int iday, int ilay)
        {
            vmr_day(iday,ilay) = vmr(dayIndices(iday),ilay);
        });
        gas_concs_day.set_vmr(gas_names[igas], vmr_day);
    }

    // Subset aerosol optics
    optical_props2_t aerosol_day;
    aerosol_day.init(k_dist.get_band_lims_wavenumber());
    aerosol_day.alloc_2str(nday, nlay);
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nday, nlay, nbnd}),
                         KOKKOS_LAMBDA (int iday, int ilay, int ibnd)
    {
        aerosol_day.tau(iday,ilay,ibnd) = aerosol.tau(dayIndices(iday),ilay,ibnd);
        aerosol_day.ssa(iday,ilay,ibnd) = aerosol.ssa(dayIndices(iday),ilay,ibnd);
        aerosol_day.g  (iday,ilay,ibnd) = aerosol.g  (dayIndices(iday),ilay,ibnd);
    });

    // Subset cloud optics
    // TODO: nbnd -> ngpt once we pass sub-sampled cloud state
    optical_props2_t clouds_day;
    clouds_day.init(k_dist.get_band_lims_wavenumber(), k_dist.get_band_lims_gpoint());
    clouds_day.alloc_2str(nday, nlay);
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nday, nlay, ngpt}),
                         KOKKOS_LAMBDA (int iday, int ilay, int igpt)
    {
        clouds_day.tau(iday,ilay,igpt) = clouds.tau(dayIndices(iday),ilay,igpt);
        clouds_day.ssa(iday,ilay,igpt) = clouds.ssa(dayIndices(iday),ilay,igpt);
        clouds_day.g  (iday,ilay,igpt) = clouds.g  (dayIndices(iday),ilay,igpt);
    });

    // RRTMGP assumes surface albedos have a screwy dimension ordering
    // for some strange reason, so we need to transpose these; also do
    // daytime subsetting in the same kernel
    real2d_k sfc_alb_dir_T("sfc_alb_dir", nbnd, nday);
    real2d_k sfc_alb_dif_T("sfc_alb_dif", nbnd, nday);
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nbnd, nday}),
                         KOKKOS_LAMBDA (int ibnd, int icol)
    {
        sfc_alb_dir_T(ibnd,icol) = sfc_alb_dir(dayIndices(icol),ibnd);
        sfc_alb_dif_T(ibnd,icol) = sfc_alb_dif(dayIndices(icol),ibnd);
    });

    // Temporaries we need for daytime-only fluxes
    real2d_k flux_up_day("flux_up_day", nday, nlay+1);
    real2d_k flux_dn_day("flux_dn_day", nday, nlay+1);
    real2d_k flux_dn_dir_day("flux_dn_dir_day", nday, nlay+1);
    real3d_k bnd_flux_up_day("bnd_flux_up_day", nday, nlay+1, nbnd);
    real3d_k bnd_flux_dn_day("bnd_flux_dn_day", nday, nlay+1, nbnd);
    real3d_k bnd_flux_dn_dir_day("bnd_flux_dn_dir_day", nday, nlay+1, nbnd);
    fluxes_t fluxes_day;
    fluxes_day.flux_up         = flux_up_day;
    fluxes_day.flux_dn         = flux_dn_day;
    fluxes_day.flux_dn_dir     = flux_dn_dir_day;
    fluxes_day.bnd_flux_up     = bnd_flux_up_day;
    fluxes_day.bnd_flux_dn     = bnd_flux_dn_day;
    fluxes_day.bnd_flux_dn_dir = bnd_flux_dn_dir_day;

    // Allocate space for optical properties
    optical_props2_t optics;
    optics.alloc_2str(nday, nlay, k_dist);

    optical_props2_t optics_no_aerosols;
    if (extra_clnsky_diag) {
        // Allocate space for optical properties (no aerosols)
        optics_no_aerosols.alloc_2str(nday, nlay, k_dist);
    }

    // Limit temperatures for gas optics look-up tables
    real2d_k t_lay_limited("t_lay_limited", nday, nlay);
    limit_to_bounds_2d(t_lay_day, k_dist_sw_k->get_temp_min(),
                       k_dist_sw_k->get_temp_max(), t_lay_limited);

    // Do gas optics
    real2d_k toa_flux("toa_flux", nday, ngpt);
    real3d_k col_gas("col_gas", ncol, nlay, k_dist.get_ngas()+1);
    bool top_at_1 = false;
    Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(int, bool& val)
    {
        val |= p_lay(0, 0) < p_lay(0, nlay-1);
    }, Kokkos::LOr<bool>(top_at_1));

    k_dist.gas_optics(nday, nlay, top_at_1, p_lay_day, p_lev_day,
                      t_lay_limited, gas_concs_day, col_gas, optics, toa_flux);
    if (extra_clnsky_diag) {
        k_dist.gas_optics(nday, nlay, top_at_1, p_lay_day, p_lev_day,
                          t_lay_limited, gas_concs_day, col_gas, optics_no_aerosols, toa_flux);
    }

    // Apply tsi_scaling
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nday, ngpt}),
                         KOKKOS_LAMBDA (int iday, int igpt)
    {
        toa_flux(iday,igpt) = tsi_scaling * toa_flux(iday,igpt);
    });

    if (extra_clnclrsky_diag) {
        // Compute clear-clean-sky (just gas) fluxes on daytime columns
        rte_sw(optics, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);
        // Expand daytime fluxes to all columns
        Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nday, nlay+1}),
                         KOKKOS_LAMBDA (int iday, int ilev)
        {
            int icol = dayIndices(iday);
            clnclrsky_flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
            clnclrsky_flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
            clnclrsky_flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
        });
    }

    // Combine gas and aerosol optics
    aerosol_day.delta_scale();
    aerosol_day.increment(optics);

    // Compute clearsky (gas + aerosol) fluxes on daytime columns
    rte_sw(optics, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);

    // Expand daytime fluxes to all columns
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nday, nlay+1}),
                         KOKKOS_LAMBDA (int iday, int ilev)
    {
        int icol = dayIndices(iday);
        clrsky_flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
        clrsky_flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
        clrsky_flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
    });

    // Now merge in cloud optics and do allsky calculations

    // Combine gas and cloud optics
    clouds_day.delta_scale();
    clouds_day.increment(optics);

    // Compute fluxes on daytime columns
    rte_sw(optics, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);

    // Expand daytime fluxes to all columns
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nday, nlay+1}),
                         KOKKOS_LAMBDA (int iday, int ilev)
    {
        int icol = dayIndices(iday);
        flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
        flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
        flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
    });
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nday, nlay+1, nbnd}),
                         KOKKOS_LAMBDA (int iday, int ilev, int ibnd)
    {
        int icol = dayIndices(iday);
        bnd_flux_up    (icol,ilev,ibnd) = bnd_flux_up_day    (iday,ilev,ibnd);
        bnd_flux_dn    (icol,ilev,ibnd) = bnd_flux_dn_day    (iday,ilev,ibnd);
        bnd_flux_dn_dir(icol,ilev,ibnd) = bnd_flux_dn_dir_day(iday,ilev,ibnd);
    });

    if (extra_clnsky_diag) {
        // First increment clouds in optics_no_aerosols
        clouds_day.increment(optics_no_aerosols);
        // Compute cleansky (gas + clouds) fluxes on daytime columns
        rte_sw(optics_no_aerosols, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);
        // Expand daytime fluxes to all columns
        Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nday, nlay+1}),
                         KOKKOS_LAMBDA (int iday, int ilev)
        {
            int icol = dayIndices(iday);
            clnsky_flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
            clnsky_flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
            clnsky_flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
        });
    }
}


void
rrtmgp_lw (const int ncol,
           const int nlay,
           gas_optics_t& k_dist,
           real2d_k& p_lay, real2d_k& t_lay,
           real2d_k& p_lev, real2d_k& t_lev,
           real1d_k& t_sfc, real1d_k& sfc_emis, real1d_k& lw_src,
           gas_concs_t& gas_concs,
           optical_props1_t& aerosol,
           optical_props1_t& clouds,
           fluxes_t& fluxes,
           fluxes_broadband_t& clnclrsky_fluxes,
           fluxes_broadband_t& clrsky_fluxes,
           fluxes_broadband_t& clnsky_fluxes,
           const bool extra_clnclrsky_diag,
           const bool extra_clnsky_diag)
{
    // Problem size
    int nbnd = k_dist.get_nband();

    // Associate local pointers for fluxes
    auto& flux_up           = fluxes.flux_up;
    auto& flux_dn           = fluxes.flux_dn;
    auto& bnd_flux_up       = fluxes.bnd_flux_up;
    auto& bnd_flux_dn       = fluxes.bnd_flux_dn;
    auto& clnclrsky_flux_up = clnclrsky_fluxes.flux_up;
    auto& clnclrsky_flux_dn = clnclrsky_fluxes.flux_dn;
    auto& clrsky_flux_up    = clrsky_fluxes.flux_up;
    auto& clrsky_flux_dn    = clrsky_fluxes.flux_dn;
    auto& clnsky_flux_up    = clnsky_fluxes.flux_up;
    auto& clnsky_flux_dn    = clnsky_fluxes.flux_dn;

    // Reset fluxes to zero
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlay+1}),
                         KOKKOS_LAMBDA (int icol, int ilev)
    {
        flux_up(icol, ilev)           = 0;
        flux_dn(icol, ilev)           = 0;
        clnclrsky_flux_up(icol, ilev) = 0;
        clnclrsky_flux_dn(icol, ilev) = 0;
        clrsky_flux_up(icol, ilev)    = 0;
        clrsky_flux_dn(icol, ilev)    = 0;
        clnsky_flux_up(icol, ilev)    = 0;
        clnsky_flux_dn(icol, ilev)    = 0;
    });
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay+1, nbnd}),
                         KOKKOS_LAMBDA (int icol, int ilev, int ibnd)
    {
        bnd_flux_up(icol, ilev, ibnd) = 0;
        bnd_flux_dn(icol, ilev, ibnd) = 0;
    });

    // Allocate space for optical properties
    optical_props1_t optics;
    optics.alloc_1scl(ncol, nlay, k_dist);

    optical_props1_t optics_no_aerosols;
    if (extra_clnsky_diag) {
        // Allocate space for optical properties (no aerosols)
        optics_no_aerosols.alloc_1scl(ncol, nlay, k_dist);
    }

    bool top_at_1 = false;
    Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(int, bool& val)
    {
        val |= p_lay(0, 0) < p_lay(0, nlay-1);
    }, Kokkos::LOr<bool>(top_at_1));

    // Boundary conditions
    //=====================================================================
    source_func_t lw_sources;
    lw_sources.alloc(ncol, nlay, k_dist);

    /*
    // Surface LW source
    // AML NOTE: This is removed in EAMXX, LSM doesn't transfer its lw_src?
    auto d_lw_src = lw_sources.sfc_source;
    Kokkos::parallel_for(ncol, KOKKOS_LAMBDA (int icol)
    {
        d_lw_src(icol, 0) = lw_src(icol);
    });
    */

    // Surface temperature
    // AML NOTE: We already populate this when initializing data


    // Surface emissivity (transposed in RRTMGP)
    // AML NOTE: This transfer was removed in EAMXX, LSM doesn't transfer its emis_sfc?
    real2d_k emis_sfc_T("emis_sfc",nbnd,ncol);
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nbnd}),
                         KOKKOS_LAMBDA (int icol, int ibnd)
    {
        emis_sfc_T(ibnd,icol) = sfc_emis(icol);
    });
    //Kokkos::deep_copy(emis_sfc_T, 0.98);

    // Get Gaussian quadrature weights
    // Weights and angle secants for first order (k=1) Gaussian quadrature.
    //   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
    //   after Abramowitz & Stegun 1972, page 921
    int constexpr max_gauss_pts = 4;
    RealT gauss_Ds_host_raw[max_gauss_pts][max_gauss_pts] = { {1.66, 1.18350343, 1.09719858, 1.06056257},
                                                              {0.  , 2.81649655, 1.69338507, 1.38282560},
                                                              {0.  , 0.        , 4.70941630, 2.40148179},
                                                              {0.  , 0.        , 0.        , 7.15513024} };
    realHost2d_k gauss_Ds_host(&gauss_Ds_host_raw[0][0], max_gauss_pts, max_gauss_pts);

    RealT gauss_wts_host_raw[max_gauss_pts][max_gauss_pts] = { {0.5, 0.3180413817, 0.2009319137, 0.1355069134},
                                                               {0. , 0.1819586183, 0.2292411064, 0.2034645680},
                                                               {0. , 0.          , 0.0698269799, 0.1298475476},
                                                               {0. , 0.          , 0.          , 0.0311809710} };
    realHost2d_k gauss_wts_host(&gauss_wts_host_raw[0][0],max_gauss_pts,max_gauss_pts);

    real2d_k gauss_Ds ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
    real2d_k gauss_wts("gauss_wts",max_gauss_pts,max_gauss_pts);
    Kokkos::deep_copy(gauss_Ds,  gauss_Ds_host);
    Kokkos::deep_copy(gauss_wts, gauss_wts_host);

    // Limit temperatures for gas optics look-up tables
    real2d_k t_lay_limited("t_lay_limited", ncol, nlay  );
    real2d_k t_lev_limited("t_lev_limited", ncol, nlay+1);
    limit_to_bounds_2d(t_lay, k_dist.get_temp_min(),
                       k_dist.get_temp_max(), t_lay_limited);
    limit_to_bounds_2d(t_lev, k_dist.get_temp_min(),
                       k_dist.get_temp_max(), t_lev_limited);

    // Do gas optics
    real3d_k col_gas("col_gas", ncol, nlay, k_dist.get_ngas()+1);
    k_dist.gas_optics(ncol, nlay, top_at_1, p_lay, p_lev, t_lay_limited,
                      t_sfc, gas_concs, col_gas, optics, lw_sources, view_t<RealT**>(), t_lev_limited);
    if (extra_clnsky_diag) {
        k_dist.gas_optics(ncol, nlay, top_at_1, p_lay, p_lev, t_lay_limited,
                          t_sfc, gas_concs, col_gas, optics_no_aerosols, lw_sources, view_t<RealT**>(), t_lev_limited);
    }

    if (extra_clnclrsky_diag) {
        // Compute clean-clear-sky fluxes before we add in aerosols and clouds
        rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, optics, top_at_1, lw_sources, emis_sfc_T, clnclrsky_fluxes);
    }

    // Combine gas and aerosol optics
    aerosol.increment(optics);

    // Compute clear-sky fluxes before we add in clouds
    rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, optics, top_at_1, lw_sources, emis_sfc_T, clrsky_fluxes);

    // Combine gas and cloud optics
    clouds.increment(optics);

    // Compute allsky fluxes
    rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, optics, top_at_1, lw_sources, emis_sfc_T, fluxes);

    if (extra_clnsky_diag) {
        // First increment clouds in optics_no_aerosols
        clouds.increment(optics_no_aerosols);
        // Compute clean-sky fluxes
        rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, optics_no_aerosols, top_at_1, lw_sources, emis_sfc_T, clnsky_fluxes);
    }
}


void
compute_cloud_area (int ncol,
                    int nlay,
                    int ngpt,
                    const RealT pmin,
                    const RealT pmax,
                    const real2d_k& pmid,
                    const real3d_k& cld_tau_gpt,
                    real1d_k& cld_area)
{
    // Subcolumn binary cld mask; if any layers with pressure between pmin and pmax are cloudy
    // then 2d subcol mask is 1, otherwise it is 0
    real2d_k subcol_mask("subcol_mask", ncol, ngpt);
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlay, ngpt}),
                         KOKKOS_LAMBDA (int icol, int ilay, int igpt)
    {
        // NOTE: using plev would need to assume level ordering (top to bottom or bottom to top), but
        // using play/pmid does not
        if (cld_tau_gpt(icol,ilay,igpt) > 0 && pmid(icol,ilay) >= pmin && pmid(icol,ilay) < pmax) {
            subcol_mask(icol,igpt) = 1.0;
        } else {
            subcol_mask(icol,igpt) = 0.0;
        }
    });
    // Compute average over subcols to get cloud area
    auto ngpt_inv = 1.0 / ngpt;
    Kokkos::deep_copy(cld_area, 0);
    Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol)
    {
        // This loop needs to be serial because of the atomic reduction
        for (int igpt = 0; igpt < ngpt; ++igpt) {
            cld_area(icol) += subcol_mask(icol,igpt) * ngpt_inv;
        }
    });
}


int
get_wavelength_index_sw (RealT wavelength) { return get_wavelength_index(*k_dist_sw_k, wavelength); }


int
get_wavelength_index_lw (RealT wavelength) { return get_wavelength_index(*k_dist_lw_k, wavelength); }


int
get_wavelength_index (optical_props_t& kdist,
                      RealT wavelength)
{
    // Get wavelength bounds for all wavelength bands
    auto band_lims_wvn = kdist.get_band_lims_wavenumber();
    real2d_k wavelength_bounds("wavelength_bounds",band_lims_wvn.extent(0), band_lims_wvn.extent(1));
    wavelength_bounds = kdist.get_band_lims_wavelength();

    // Find the band index for the specified wavelength
    // Note that bands are stored in wavenumber space, units of cm-1, so if we are passed wavelength
    // in units of meters, we need a conversion factor of 10^2
    int nbnds = kdist.get_nband();
    int band_index = -1;
    Kokkos::parallel_reduce(nbnds, KOKKOS_LAMBDA(int ibnd, int& band_index_inner)
    {
        if (wavelength_bounds(0,ibnd) < wavelength_bounds(1,ibnd)) {
            if (wavelength_bounds(0,ibnd) <= wavelength * 1e2 && wavelength * 1e2 <= wavelength_bounds(1,ibnd)) {
                band_index_inner = ibnd;
            }
        } else {
            if (wavelength_bounds(0,ibnd) >= wavelength * 1e2 && wavelength * 1e2 >= wavelength_bounds(1,ibnd)) {
                band_index_inner = ibnd;
            }
        }
    }, Kokkos::Max<int>(band_index));
    return band_index;
}


void
compute_aerocom_cloudtop (int ncol, int nlay ,
                          const real2d_k& tmid , const real2d_k& pmid ,
                          const real2d_k& p_del, const real2d_k& z_del, const real2d_k& qc,
                          const real2d_k& qi   , const real2d_k& rel  , const real2d_k& rei,
                          const real2d_k& cldfrac_tot      , const real2d_k& nc,
                          real1d_k& T_mid_at_cldtop        , real1d_k& p_mid_at_cldtop,
                          real1d_k& cldfrac_ice_at_cldtop  , real1d_k& cldfrac_liq_at_cldtop,
                          real1d_k& cldfrac_tot_at_cldtop  , real1d_k& cdnc_at_cldtop,
                          real1d_k& eff_radius_qc_at_cldtop, real1d_k& eff_radius_qi_at_cldtop)
{
    /* The goal of this routine is to calculate properties at cloud top
     * based on the AeroCom recommendation. See reference for routine
     * get_subcolumn_mask above, where equation 14 is used for the
     * maximum-random overlap assumption for subcolumn generation. We use
     * equation 13, the column counterpart.
     */
    // Set outputs to zero
    Kokkos::deep_copy(T_mid_at_cldtop, 0.0);
    Kokkos::deep_copy(p_mid_at_cldtop, 0.0);
    Kokkos::deep_copy(cldfrac_ice_at_cldtop, 0.0);
    Kokkos::deep_copy(cldfrac_liq_at_cldtop, 0.0);
    Kokkos::deep_copy(cldfrac_tot_at_cldtop, 0.0);
    Kokkos::deep_copy(cdnc_at_cldtop, 0.0);
    Kokkos::deep_copy(eff_radius_qc_at_cldtop, 0.0);
    Kokkos::deep_copy(eff_radius_qi_at_cldtop, 0.0);

    // Initialize the 1D "clear fraction" as 1 (totally clear)
    real1d_k aerocom_clr("aerocom_clr", ncol);
    Kokkos::deep_copy(aerocom_clr, 1.0);

    // TODO: move tunable constant to namelist
    constexpr RealT q_threshold = 0.0;  // BAD_CONSTANT!

    // TODO: move tunable constant to namelist
    constexpr RealT cldfrac_tot_threshold = 0.001;  // BAD_CONSTANT!

    // Loop over all columns in parallel
    Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol)
    {
        // Loop over all layers in serial (due to accumulative
        // product), starting at 2 (second highest) layer because the
        // highest is assumed to have no clouds
        for(int ilay = 1; ilay < nlay; ++ilay) {
            // Only do the calculation if certain conditions are met
            if((qc(icol, ilay) + qi(icol, ilay)) > q_threshold &&
               (cldfrac_tot(icol, ilay) > cldfrac_tot_threshold)) {
                /* PART I: Probabilistically determining cloud top */
                // Populate aerocom_tmp as the clear-sky fraction
                // probability of this level, where aerocom_clr is that of
                // the previous level
                auto aerocom_tmp = aerocom_clr(icol) *
                                   (1.0 - std::max(cldfrac_tot(icol, ilay - 1),
                                                   cldfrac_tot(icol, ilay))) /
                                   (1.0 - std::min(cldfrac_tot(icol, ilay - 1),
                                                   1.0 - cldfrac_tot_threshold));
                // Temporary variable for probability "weights"
                auto aerocom_wts = aerocom_clr(icol) - aerocom_tmp;
                // Temporary variable for liquid "phase"
                auto aerocom_phi = qc(icol, ilay) / (qc(icol, ilay) + qi(icol, ilay));
                /* PART II: The inferred properties */
                /* In general, converting a 3D property X to a 2D cloud-top
                 * counterpart x follows: x(i) += X(i,k) * weights * Phase
                 * but X and Phase are not always needed */
                // T_mid_at_cldtop
                T_mid_at_cldtop(icol) += tmid(icol, ilay) * aerocom_wts;
                // p_mid_at_cldtop
                p_mid_at_cldtop(icol) += pmid(icol, ilay) * aerocom_wts;
                // cldfrac_ice_at_cldtop
                cldfrac_ice_at_cldtop(icol) += (1.0 - aerocom_phi) * aerocom_wts;
                // cldfrac_liq_at_cldtop
                cldfrac_liq_at_cldtop(icol) += aerocom_phi * aerocom_wts;
                // cdnc_at_cldtop
                /* We need to convert nc from 1/mass to 1/volume first, and
                 * from grid-mean to in-cloud, but after that, the
                 * calculation follows the general logic */
                // AML NOTE: p_del/z_del/g should be replaced with RHO for our dycore
                auto cdnc = nc(icol, ilay) * p_del(icol, ilay) /
                            z_del(icol, ilay) / CONST_GRAV /
                            cldfrac_tot(icol, ilay);
                cdnc_at_cldtop(icol) += cdnc * aerocom_phi * aerocom_wts;
                // eff_radius_qc_at_cldtop
                eff_radius_qc_at_cldtop(icol) += rel(icol, ilay) * aerocom_phi * aerocom_wts;
                // eff_radius_qi_at_cldtop
                eff_radius_qi_at_cldtop(icol) += rei(icol, ilay) * (1.0 - aerocom_phi) * aerocom_wts;
                // Reset aerocom_clr to aerocom_tmp to accumulate
                aerocom_clr(icol) = aerocom_tmp;
            }
        }
        // After the serial loop over levels, the cloudy fraction is
        // defined as (1 - aerocom_clr). This is true because
        // aerocom_clr is the result of accumulative probabilities
        // (their products)
        cldfrac_tot_at_cldtop(icol) = 1.0 - aerocom_clr(icol);
    });
}

}  // namespace rrtmgp


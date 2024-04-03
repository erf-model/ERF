/*
 * RTE-RRTMGP radiation model interface to ERF
 * The orginal code is developed by RobertPincus, and the code is open source available at:
 *                        https://github.com/earth-system-radiation/rte-rrtmgp
 * Please reference to the following paper,
 *                        https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019MS001621
 * NOTE: we use the C++ version of RTE-RRTMGP, which is reimplemented the original Fortran
 * code using C++ YAKL for CUDA, HiP and SYCL application by E3SM ECP team, the C++ version
 * of the rte-rrtmgp code is located at:
 *                       https://github.com/E3SM-Project/rte-rrtmgp
 * The RTE-RRTMGP uses BSD-3-Clause Open Source License, if you want to make changes,
 * and modifications to the code, please refer to BSD-3-Clause Open Source License.
 */
#include <string>
#include <vector>
#include <memory>

#include "Radiation.H"
#include "m2005_effradius.H"
#include <AMReX_GpuContainers.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_TableData.H>
#include <AMReX_MultiFabUtil.H>
#include "ERF_Constants.H"
#include "IndexDefines.H"
#include "DataStruct.H"
#include "EOS.H"
#include "TileNoZ.H"

using namespace amrex;
using yakl::intrinsics::size;
using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

namespace internal {
    void initial_fluxes (int nz, int nlay, int nbands, FluxesByband& fluxes)
    {
        fluxes.flux_up = real2d("flux_up", nz, nlay+1);
        fluxes.flux_dn = real2d("flux_dn", nz, nlay+1);
        fluxes.flux_net = real2d("flux_net", nz, nlay+1);
        fluxes.flux_dn_dir = real2d("flux_dn_dir", nz, nlay+1);
        fluxes.bnd_flux_up = real3d("flux_up", nz, nlay+1, nbands);
        fluxes.bnd_flux_dn = real3d("flux_dn", nz, nlay+1, nbands);
        fluxes.bnd_flux_net = real3d("flux_net", nz, nlay+1, nbands);
        fluxes.bnd_flux_dn_dir = real3d("flux_dn_dir", nz, nlay+1, nbands);
    }

    void expand_day_fluxes (const FluxesByband& daytime_fluxes,
                            FluxesByband& expanded_fluxes,
                            const int1d& day_indices)
    {
        auto ncol  = size(daytime_fluxes.bnd_flux_up, 1);
        auto nlev  = size(daytime_fluxes.bnd_flux_up, 2);
        auto nbnds = size(daytime_fluxes.bnd_flux_up, 3);

        int1d nday_1d("nday_1d", 1),
            nday_host("nday_host",1);
        yakl::memset(nday_1d, 0);
        parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA (int icol)
        {
            if (day_indices(icol) > 0) nday_1d(1)++;
            printf("daynight indices(check): %d, %d, %d\n",icol,day_indices(icol),nday_1d(1));
        });

        nday_1d.deep_copy_to(nday_host);
        auto nday = nday_host(1);
        parallel_for(SimpleBounds<3>(nday, nlev, nbnds), YAKL_LAMBDA (int iday, int ilev, int ibnd)
        {
            // Map daytime index to proper column index
            // auto icol = day_indices(iday);
            auto icol = iday;
            // Expand broadband fluxes
            expanded_fluxes.flux_up(icol,ilev) = daytime_fluxes.flux_up(iday,ilev);
            expanded_fluxes.flux_dn(icol,ilev) = daytime_fluxes.flux_dn(iday,ilev);
            expanded_fluxes.flux_net(icol,ilev) = daytime_fluxes.flux_net(iday,ilev);
            expanded_fluxes.flux_dn_dir(icol,ilev) = daytime_fluxes.flux_dn_dir(iday,ilev);

            // Expand band-by-band fluxes
            expanded_fluxes.bnd_flux_up(icol,ilev,ibnd) = daytime_fluxes.bnd_flux_up(iday,ilev,ibnd);
            expanded_fluxes.bnd_flux_dn(icol,ilev,ibnd) = daytime_fluxes.bnd_flux_dn(iday,ilev,ibnd);
            expanded_fluxes.bnd_flux_net(icol,ilev,ibnd) = daytime_fluxes.bnd_flux_net(iday,ilev,ibnd);
            expanded_fluxes.bnd_flux_dn_dir(icol,ilev,ibnd) = daytime_fluxes.bnd_flux_dn_dir(iday,ilev,ibnd);
        });
    }

    // Utility function to reorder an array given a new indexing
    void reordered (const real1d& array_in, const int1d& new_indexing, const real1d& array_out)
    {
        // Reorder array based on input index mapping, which maps old indices to new
        parallel_for(SimpleBounds<1>(size(array_in, 1)), YAKL_LAMBDA (int i)
        {
            array_out(i) = array_in(new_indexing(i));
        });
    }
}

// init
void Radiation::initialize (const MultiFab& cons_in,
                            const MultiFab& qmoist,
                            const BoxArray& grids,
                            const Geometry& geom,
                            const Real& dt_advance,
                            const bool& do_sw_rad,
                            const bool& do_lw_rad,
                            const bool& do_aero_rad,
                            const bool& do_snow_opt,
                            const bool& is_cmip6_volcano)
{
    m_geom = geom;
    m_box = grids;

    auto dz   = m_geom.CellSize(2);
    auto lowz = m_geom.ProbLo(2);

    dt = dt_advance;

    do_short_wave_rad = do_sw_rad;
    do_long_wave_rad = do_lw_rad;
    do_aerosol_rad = do_aero_rad;
    do_snow_optics = do_snow_opt;
    is_cmip6_volc = is_cmip6_volcano;

    for ( MFIter mfi(cons_in, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const auto& box3d = mfi.tilebox();

        nlev = box3d.length(2);
        ncol = box3d.length(0)*box3d.length(1);
    }

    ngas = active_gases.size();

    // initialize cloud, aerosol, and radiation
    radiation.initialize(ngas, active_gases,
                         rrtmgp_coefficients_file_sw.c_str(),
                         rrtmgp_coefficients_file_lw.c_str());

    // initialize the radiation data
    nswbands = radiation.get_nband_sw();
    nswgpts  = radiation.get_ngpt_sw();
    nlwbands = radiation.get_nband_lw();
    nlwgpts  = radiation.get_ngpt_lw();

    rrtmg_to_rrtmgp = int1d("rrtmg_to_rrtmgp",14);
    parallel_for(14, YAKL_LAMBDA (int i)
    {
        if (i == 1) {
            rrtmg_to_rrtmgp(i) = 13;
        } else {
            rrtmg_to_rrtmgp(i) = i - 1;
        }
    });

    tmid = real2d("tmid", ncol, nlev);
    pmid = real2d("pmid", ncol, nlev);
    pdel = real2d("pdel", ncol, nlev);

    pint = real2d("pint", ncol, nlev+1);
    tint = real2d("tint", ncol, nlev+1);

    qt   = real2d("qt", ncol, nlev);
    qc   = real2d("qc", ncol, nlev);
    qi   = real2d("qi", ncol, nlev);
    qn   = real2d("qn", ncol, nlev);
    zi   = real2d("zi", ncol, nlev);

    // Get the temperature, density, theta, qt and qp from input
    for ( MFIter mfi(cons_in, false); mfi.isValid(); ++mfi) {
        auto states_array = cons_in.array(mfi);
        auto qmoist_array = qmoist.array(mfi);

        const auto& box3d = mfi.tilebox();
        auto nx = box3d.length(0);
        //auto ny = box3d.length(1);
        // Get pressure, theta, temperature, density, and qt, qp
        amrex::ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            auto icol = j*nx+i+1;
            auto ilev = k+1;
            qt(icol,ilev)   = qmoist_array(i,j,k,0); //states_array(i,j,k,RhoQt_comp)/states_array(i,j,k,Rho_comp);
            qc(icol,ilev)   = qmoist_array(i,j,k,2); //qmoist_array(i,j,k,1);
            qi(icol,ilev)   = qmoist_array(i,j,k,3); //qmoist_array(i,j,k,2);
            qn(icol,ilev)   = qmoist_array(i,j,k,1) + qmoist_array(i,j,k,2);
            tmid(icol,ilev) = getTgivenRandRTh(states_array(i,j,k,Rho_comp),states_array(i,j,k,RhoTheta_comp));
            pmid(icol,ilev) = getPgivenRTh(states_array(i,j,k,RhoTheta_comp))/1.0e3;
        });
    }

    parallel_for(SimpleBounds<2>(ncol, nlev+1), YAKL_LAMBDA (int icol, int ilev)
    {
        if (ilev == 1) {
            pint(icol, 1) = 2.*pmid(icol, 2) - pmid(icol, 1);
            tint(icol, 1) = 2.*tmid(icol, 2) - tmid(icol, 1);
        } else if (ilev <= nlev) {
            pint(icol, ilev) = 0.5*(pmid(icol, ilev-1) + pmid(icol, ilev));
            tint(icol, ilev) = 0.5*(tmid(icol, ilev-1) + tmid(icol, ilev));
        } else {
            pint(icol, nlev+1) = 2.*pmid(icol, nlev-1) - pmid(icol, nlev);
            tint(icol, nlev+1) = 2.*tmid(icol, nlev-1) - tmid(icol, nlev);
        }
    });

    parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
    {
        zi(icol, ilev)  = lowz + (ilev+0.5)*dz;
        pdel(icol,ilev) = pint(icol,ilev+1) - pint(icol,ilev);
    });

    albedo_dir = real2d("albedo_dir", nswbands, ncol);
    albedo_dif = real2d("albedo_dif", nswbands, ncol);

    qrs = real2d("qrs", ncol, nlev);   // shortwave radiative heating rate
    qrl = real2d("qrl", ncol, nlev);   // longwave  radiative heating rate

    // Clear-sky heating rates are not on the physics buffer, and we have no
    // reason to put them there, so declare these are regular arrays here
    qrsc = real2d("qrsc", ncol, nlev);
    qrlc = real2d("qrlc", ncol, nlev);

    int nmodes = 3;
    int nrh = 1;
    int top_lev = 1;
    naer = 4;
    std::vector<std::string> aero_names {"H2O", "N2", "O2", "O3"};
    auto geom_radius = real2d("geom_radius", ncol, nlev);
    yakl::memset(geom_radius, 0.1);

    optics.initialize(ngas, nmodes, naer, nswbands, nlwbands,
                      ncol, nlev, nrh, top_lev, aero_names, zi,
                      pmid, tmid, qt, geom_radius);

    amrex::Print() << "LW coefficents file: " << rrtmgp_coefficients_file_lw
                   << "\nSW coefficents file: " << rrtmgp_coefficients_file_sw
                   << "\nFrequency (timesteps) of Shortwave Radiation calc: " << dt
                   << "\nFrequency (timesteps) of Longwave Radiation calc:  " << dt
                   << "\nDo aerosol radiative calculations: " << do_aerosol_rad << std::endl;

}


// run radiation model
void Radiation::run ()
{
    // Temporary variable for heating rate output
    real2d hr("hr", ncol, nlev);

    // Cosine solar zenith angle for all columns in chunk
    real1d coszrs("coszrs", ncol);

    // Pointers to fields on the physics buffer
    real2d cld("cld", ncol, nlev), cldfsnow("cldfsnow", ncol, nlev),
        iclwp("iclwp", ncol, nlev), iciwp("iciwp", ncol, nlev),
        icswp("icswp", ncol, nlev), dei("dei", ncol, nlev),
        des("des", ncol, nlev), lambdac("lambdac", ncol, nlev),
        mu("mu", ncol, nlev), rei("rei", ncol, nlev), rel("rel", ncol, nlev);

    // Cloud, snow, and aerosol optical properties
    real3d cld_tau_gpt_sw("cld_tau_gpt_sw", ncol, nlev, nswgpts),
        cld_ssa_gpt_sw("cld_ssa_gpt_sw", ncol, nlev, nswgpts),
        cld_asm_gpt_sw("cld_asm_gpt_sw", ncol, nlev, nswgpts);
    real3d cld_tau_bnd_sw("cld_tau_bnd_sw", ncol, nlev, nswbands),
        cld_ssa_bnd_sw("cld_ssa_bnd_sw", ncol, nlev, nswbands),
        cld_asm_bnd_sw("cld_asm_bnd_sw", ncol, nlev, nswbands);

    real3d aer_tau_bnd_sw("aer_tau_bnd_sw", ncol, nlev, nswbands),
        aer_ssa_bnd_sw("aer_ssa_bnd_sw", ncol, nlev, nswbands),
        aer_asm_bnd_sw("aer_asm_bnd_sw", ncol, nlev, nswbands);

    real3d cld_tau_bnd_lw("cld_tau_bnd_lw", ncol, nlev, nlwbands),
        aer_tau_bnd_lw("aer_tau_bnd_lw", ncol, nlev, nlwbands);
    real3d cld_tau_gpt_lw("cld_tau_gpt_lw", ncol, nlev, nlwgpts);

    // NOTE: these are diagnostic only
    real3d liq_tau_bnd_sw("liq_tau_bnd_sw", ncol, nlev, nswbands),
        ice_tau_bnd_sw("ice_tau_bnd_sw", ncol, nlev, nswbands),
        snw_tau_bnd_sw("snw_tau_bnd_sw", ncol, nlev, nswbands);
    real3d liq_tau_bnd_lw("liq_tau_bnd_lw", ncol, nlev, nlwbands),
        ice_tau_bnd_lw("ice_tau_bnd_lw", ncol, nlev, nlwbands),
        snw_tau_bnd_lw("snw_tau_bnd_lw", ncol, nlev, nlwbands);

    // Gas volume mixing ratios
    real3d gas_vmr("gas_vmr", ngas, ncol, nlev);

    // Needed for shortwave aerosol;
    //int nday, nnight;     // Number of daylight columns
    int1d day_indices("day_indices", ncol), night_indices("night_indices", ncol);   // Indicies of daylight coumns

    // Flag to carry (QRS,QRL)*dp across time steps.
    // TODO: what does this mean?
    bool conserve_energy = true;

    // For loops over diagnostic calls
    //bool active_calls(0:N_DIAG)

    // Zero-array for cloud properties if not diagnosed by microphysics
    real2d zeros("zeros", ncol, nlev);

    int1d gpoint_bands_sw("gpoint_bands_sw", nswgpts);
    int1d gpoint_bands_lw("gpoint_bands_lw", nlwgpts);

    // Do shortwave stuff...
    if (do_short_wave_rad) {
        // Radiative fluxes
        FluxesByband fluxes_allsky, fluxes_clrsky;
        internal::initial_fluxes(ncol, nlev, nlwbands, fluxes_allsky);
        internal::initial_fluxes(ncol, nlev, nlwbands, fluxes_clrsky);

        // Get cosine solar zenith angle for current time step. ( still NOT YET implemented here)
        //      set_cosine_solar_zenith_angle(state, dt_avg, coszrs(1:ncol))
        yakl::memset(coszrs, 0.2);  // we set constant value here to avoid numerical overflow.

        // Get albedo. This uses CAM routines internally and just provides a
        // wrapper to improve readability of the code here.
        set_albedo(coszrs, albedo_dir, albedo_dif);

        // Do shortwave cloud optics calculations
        yakl::memset(cld_tau_gpt_sw, 0.);
        yakl::memset(cld_ssa_gpt_sw, 0.);
        yakl::memset(cld_asm_gpt_sw, 0.);

        // set cloud fraction to be 1, and snow fraction 0
        yakl::memset(cldfsnow, 0.0);
        yakl::memset(cld, 1.0);

        parallel_for (SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int i, int k)
        {
            iciwp(i,k) = std::min(qi(i,k)/std::max(1.0e-4,cld(i,k)),0.005)*pmid(i,k)/CONST_GRAV;
            iclwp(i,k) = std::min(qt(i,k)/std::max(1.0e-4,cld(i,k)),0.005)*pmid(i,k)/CONST_GRAV;
            icswp(i,k) = qn(i,k)/std::max(1.0e-4,cldfsnow(i,k))*pmid(i,k)/CONST_GRAV;
        });

        m2005_effradius(qc, qc, qi, qi, qt, qt, cld, pmid, tmid,
                        rel, rei, dei, lambdac, mu, des);

        // calculate the cloud radiation
        optics.get_cloud_optics_sw(ncol, nlev, nswbands, do_snow_optics, cld,
                                   cldfsnow, iclwp, iciwp, icswp,
                                   lambdac, mu, dei, des, rel, rei,
                                   cld_tau_bnd_sw, cld_ssa_bnd_sw, cld_asm_bnd_sw,
                                   liq_tau_bnd_sw, ice_tau_bnd_sw, snw_tau_bnd_sw);

        // Now reorder bands to be consistent with RRTMGP
        // We need to fix band ordering because the old input files assume RRTMG
        // band ordering, but this has changed in RRTMGP.
        // TODO: fix the input files themselves!
        real1d cld_tau_bnd_sw_1d("cld_tau_bnd_sw_1d", nswbands);
        real1d cld_ssa_bnd_sw_1d("cld_ssa_bnd_sw_1d", nswbands);
        real1d cld_asm_bnd_sw_1d("cld_asm_bnd_sw_1d", nswbands);
        real1d cld_tau_bnd_sw_o_1d("cld_tau_bnd_sw_1d", nswbands);
        real1d cld_ssa_bnd_sw_o_1d("cld_ssa_bnd_sw_1d", nswbands);
        real1d cld_asm_bnd_sw_o_1d("cld_asm_bnd_sw_1d", nswbands);

        parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilay)
        {
            for (auto ibnd = 1; ibnd <= nswbands; ++ibnd) {
                cld_tau_bnd_sw_1d(ibnd) = cld_tau_bnd_sw(icol,ilay,ibnd);
                cld_ssa_bnd_sw_1d(ibnd) = cld_ssa_bnd_sw(icol,ilay,ibnd);
                cld_asm_bnd_sw_1d(ibnd) = cld_asm_bnd_sw(icol,ilay,ibnd);
            }
            internal::reordered(cld_tau_bnd_sw_1d, rrtmg_to_rrtmgp, cld_tau_bnd_sw_o_1d);
            internal::reordered(cld_ssa_bnd_sw_1d, rrtmg_to_rrtmgp, cld_ssa_bnd_sw_o_1d);
            internal::reordered(cld_asm_bnd_sw_1d, rrtmg_to_rrtmgp, cld_asm_bnd_sw_o_1d);
            for (auto ibnd = 1; ibnd <= nswbands; ++ibnd) {
                cld_tau_bnd_sw(icol,ilay,ibnd) = cld_tau_bnd_sw_o_1d(ibnd);
                cld_ssa_bnd_sw(icol,ilay,ibnd) = cld_ssa_bnd_sw_o_1d(ibnd);
                cld_asm_bnd_sw(icol,ilay,ibnd) = cld_asm_bnd_sw_o_1d(ibnd);
            }
        });

        // And now do the MCICA sampling to get cloud optical properties by
        // gpoint/cloud state
        radiation.get_gpoint_bands_sw(gpoint_bands_sw);

        optics.sample_cloud_optics_sw(ncol, nlev, nswgpts, gpoint_bands_sw,
                                      pmid, cld, cldfsnow,
                                      cld_tau_bnd_sw, cld_ssa_bnd_sw, cld_asm_bnd_sw,
                                      cld_tau_gpt_sw, cld_ssa_gpt_sw, cld_asm_gpt_sw);

        // Aerosol needs night indices
        // TODO: remove this dependency, it's just used to mask aerosol outputs
        set_daynight_indices(coszrs, day_indices, night_indices);
        int1d nday("nday",1);
        int1d nnight("nnight",1);
        yakl::memset(nday, 0);
        yakl::memset(nnight, 0);
        for (auto icol=1; icol<ncol; ++icol) {
            if (day_indices(icol) > 0) nday(1)++;
            if (night_indices(icol) > 0) nnight(1)++;
        }

        // get aerosol optics
        do_aerosol_rad = true;
        {
            // Get gas concentrations
            get_gas_vmr(active_gases, gas_vmr);

            // Get aerosol optics
            if (do_aerosol_rad) {
                yakl::memset(aer_tau_bnd_sw, 0.);
                yakl::memset(aer_ssa_bnd_sw, 0.);
                yakl::memset(aer_asm_bnd_sw, 0.);

                real2d clear_rh("clear_rh",ncol, nswbands);
                yakl::memset(clear_rh, 0.01);

                optics.set_aerosol_optics_sw(0, ncol, nlev, nswbands, dt, night_indices,
                                             is_cmip6_volc, aer_tau_bnd_sw, aer_ssa_bnd_sw, aer_asm_bnd_sw, clear_rh);

                // Now reorder bands to be consistent with RRTMGP
                // TODO: fix the input files themselves!
                real1d aer_tau_bnd_sw_1d("cld_tau_bnd_sw_1d", nswbands);
                real1d aer_ssa_bnd_sw_1d("cld_ssa_bnd_sw_1d", nswbands);
                real1d aer_asm_bnd_sw_1d("cld_asm_bnd_sw_1d", nswbands);
                real1d aer_tau_bnd_sw_o_1d("cld_tau_bnd_sw_1d", nswbands);
                real1d aer_ssa_bnd_sw_o_1d("cld_ssa_bnd_sw_1d", nswbands);
                real1d aer_asm_bnd_sw_o_1d("cld_asm_bnd_sw_1d", nswbands);

                parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilay)
                {
                    for (auto ibnd = 1; ibnd < nswbands; ++ibnd) {
                        aer_tau_bnd_sw_1d(ibnd) = aer_tau_bnd_sw(icol,ilay,ibnd);
                        aer_ssa_bnd_sw_1d(ibnd) = aer_ssa_bnd_sw(icol,ilay,ibnd);
                        aer_asm_bnd_sw_1d(ibnd) = aer_asm_bnd_sw(icol,ilay,ibnd);
                    }
                    internal::reordered(aer_tau_bnd_sw_1d, rrtmg_to_rrtmgp, aer_tau_bnd_sw_o_1d);
                    internal::reordered(aer_ssa_bnd_sw_1d, rrtmg_to_rrtmgp, aer_ssa_bnd_sw_o_1d);
                    internal::reordered(aer_asm_bnd_sw_1d, rrtmg_to_rrtmgp, aer_asm_bnd_sw_o_1d);
                    for (auto ibnd = 1; ibnd < nswbands; ++ibnd) {
                        aer_tau_bnd_sw(icol,ilay,ibnd) = aer_tau_bnd_sw_o_1d(ibnd);
                        aer_ssa_bnd_sw(icol,ilay,ibnd) = aer_ssa_bnd_sw_o_1d(ibnd);
                        aer_asm_bnd_sw(icol,ilay,ibnd) = aer_asm_bnd_sw_o_1d(ibnd);
                    }
                });
            } else {
                yakl::memset(aer_tau_bnd_sw, 0.);
                yakl::memset(aer_ssa_bnd_sw, 0.);
                yakl::memset(aer_asm_bnd_sw, 0.);
            }

         yakl::memset(cld_tau_gpt_sw, 0.);
         yakl::memset(cld_ssa_gpt_sw, 0.);
         yakl::memset(cld_asm_gpt_sw, 0.);

         // Call the shortwave radiation driver
         radiation_driver_sw(ncol, gas_vmr,
                             pmid, pint, tmid, albedo_dir, albedo_dif, coszrs,
                             cld_tau_gpt_sw, cld_ssa_gpt_sw, cld_asm_gpt_sw,
                             aer_tau_bnd_sw, aer_ssa_bnd_sw, aer_asm_bnd_sw,
                             fluxes_allsky, fluxes_clrsky, qrs, qrsc);
        }
    } else {
        // Conserve energy
        if (conserve_energy) {
            parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
            {
                qrs(icol,ilev) = qrs(icol,ilev)/pdel(icol,ilev);
            });
        }
    }  // dosw

    // Do longwave stuff...
    if (do_long_wave_rad) {
        // Allocate longwave outputs; why is this not part of the fluxes_t object?
        FluxesByband fluxes_allsky, fluxes_clrsky;
        internal::initial_fluxes(ncol, nlev, nlwbands, fluxes_allsky);
        internal::initial_fluxes(ncol, nlev, nlwbands, fluxes_clrsky);

        // NOTE: fluxes defined at interfaces, so initialize to have vertical dimension nlev_rad+1
        yakl::memset(cld_tau_gpt_lw, 0.);

        optics.get_cloud_optics_lw(ncol, nlev, nlwbands, do_snow_optics, cld, cldfsnow, iclwp, iciwp, icswp,
                                   lambdac, mu, dei, des, rei,
                                   cld_tau_bnd_lw, liq_tau_bnd_lw, ice_tau_bnd_lw, snw_tau_bnd_lw);

        radiation.get_gpoint_bands_lw(gpoint_bands_lw);

        optics.sample_cloud_optics_lw(ncol, nlev, nlwgpts, gpoint_bands_lw,
                                      pmid, cld, cldfsnow,
                                      cld_tau_bnd_lw, cld_tau_gpt_lw);

        // Get gas concentrations
        get_gas_vmr(active_gases, gas_vmr);

        // Get aerosol optics
        yakl::memset(aer_tau_bnd_lw, 0.);
        if (do_aerosol_rad) {
            aer_rad.aer_rad_props_lw(is_cmip6_volc, 0, dt, zi, aer_tau_bnd_lw, clear_rh);
        }

        // Call the longwave radiation driver to calculate fluxes and heating rates
        radiation_driver_lw(ncol, nlev, gas_vmr, pmid, pint, tmid, tint, cld_tau_gpt_lw, aer_tau_bnd_lw,
                            fluxes_allsky, fluxes_clrsky, qrl, qrlc);
    }
    else {
        // Conserve energy (what does this mean exactly?)
        if (conserve_energy) {
            parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
            {
                qrl(icol,ilev) = qrl(icol,ilev)/pdel(icol,ilev);
            });
        }
    } // dolw

    // Compute heating rate for dtheta/dt
    parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilay)
    {
        hr(icol,ilay) = (qrs(icol,ilay) + qrl(icol,ilay)) / Cp_d * (1.e5 / std::pow(pmid(icol,ilay), R_d/Cp_d));
    });

    // convert radiative heating rates to Q*dp for energy conservation
    if (conserve_energy) {
        parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
        {
            qrs(icol,ilev) = qrs(icol,ilev)*pdel(icol,ilev);
            qrl(icol,ilev) = qrl(icol,ilev)*pdel(icol,ilev);
        });
    }
}

void Radiation::radiation_driver_sw (int ncol, const real3d& gas_vmr,
                                     const real2d& pmid, const real2d& pint, const real2d& tmid,
                                     const real2d& albedo_dir, const real2d& albedo_dif, const real1d& coszrs,
                                     const real3d& cld_tau_gpt, const real3d& cld_ssa_gpt, const real3d& cld_asm_gpt,
                                     const real3d& aer_tau_bnd, const real3d& aer_ssa_bnd, const real3d& aer_asm_bnd,
                                     FluxesByband& fluxes_clrsky, FluxesByband& fluxes_allsky, const real2d& qrs,
                                     const real2d& qrsc)
{
    // Incoming solar radiation, scaled for solar zenith angle
    // and earth-sun distance
    real2d solar_irradiance_by_gpt("solar_irradiance_by_gpt",ncol,nswgpts);

    // Gathered indicies of day and night columns
    // chunk_column_index = day_indices(daylight_column_index)
    int1d day_indices("day_indices",ncol), night_indices("night_indices", ncol);   // Indicies of daylight coumns

    real1d coszrs_day("coszrs_day", ncol);
    real2d albedo_dir_day("albedo_dir_day", nswbands, ncol), albedo_dif_day("albedo_dif_day", nswbands, ncol);
    real2d pmid_day("pmid_day", ncol, nlev);
    real2d tmid_day("tmid_day", ncol, nlev);
    real2d pint_day("pint_day", ncol, nlev+1);

    real3d gas_vmr_day("gas_vmr_day", ngas, ncol, nlev);
    real3d gas_vmr_rad("gas_vmr_rad", ngas, ncol, nlev);

    real3d cld_tau_gpt_day("cld_tau_gpt_day", ncol, nlev, nswgpts);
    real3d cld_ssa_gpt_day("cld_ssa_gpt_day", ncol, nlev, nswgpts);
    real3d cld_asm_gpt_day("cld_asm_gpt_day", ncol, nlev, nswgpts);
    real3d aer_tau_bnd_day("aer_tau_bnd_day", ncol, nlev, nswbands);
    real3d aer_ssa_bnd_day("aer_ssa_bnd_day", ncol, nlev, nswbands);
    real3d aer_asm_bnd_day("aer_asm_bnd_day", ncol, nlev, nswbands);

    real3d cld_tau_gpt_rad("cld_tau_gpt_rad", ncol, nlev+1, nswgpts);
    real3d cld_ssa_gpt_rad("cld_ssa_gpt_rad", ncol, nlev+1, nswgpts);
    real3d cld_asm_gpt_rad("cld_asm_gpt_rad", ncol, nlev+1, nswgpts);
    real3d aer_tau_bnd_rad("aer_tau_bnd_rad", ncol, nlev+1, nswgpts);
    real3d aer_ssa_bnd_rad("aer_ssa_bnd_rad", ncol, nlev+1, nswgpts);
    real3d aer_asm_bnd_rad("aer_asm_bnd_rad", ncol, nlev+1, nswgpts);

    // Scaling factor for total sky irradiance; used to account for orbital
    // eccentricity, and could be used to scale total sky irradiance for different
    // climates as well (i.e., paleoclimate simulations)
    real tsi_scaling;

    if (fixed_total_solar_irradiance<0) {
        // Get orbital eccentricity factor to scale total sky irradiance
        tsi_scaling = 1.0e-3; //get_eccentricity_factor();
    } else {
        // For fixed TSI we divide by the default solar constant of 1360.9
        // At some point we will want to replace this with a method that
        // retrieves the solar constant
        tsi_scaling = fixed_total_solar_irradiance / 1360.9;
    }

    // Gather night/day column indices for subsetting SW inputs; we only want to
    // do the shortwave radiative transfer during the daytime to save
    // computational cost (and because RRTMGP will fail for cosine solar zenith
    // angles less than or equal to zero)
    set_daynight_indices(coszrs, day_indices, night_indices);
    int1d nday("nday",1);
    int1d nnight("nnight", 1);
    yakl::memset(nday, 0);
    yakl::memset(nnight, 0);
    parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA (int icol)
    {
        if (day_indices(icol) > 0) nday(1)++;
        if (night_indices(icol) > 0) nnight(1)++;
    });

    intHost1d num_day("num_day",1);
    intHost1d num_night("num_night",1);
    nday.deep_copy_to(num_day);
    nnight.deep_copy_to(num_night);

    // If no daytime columns in this chunk, then we return zeros
    if (num_day(1) == 0) {
        //    reset_fluxes(fluxes_allsky)
        //    reset_fluxes(fluxes_clrsky)
        yakl::memset(qrs, 0.);
        yakl::memset(qrsc, 0.);
        return;
    }

    // Compress to daytime-only arrays
    parallel_for(SimpleBounds<3>(num_day(1), nlev, nswgpts), YAKL_LAMBDA (int iday, int ilev, int igpt)
    {
        auto icol = day_indices(iday);
        tmid_day(iday,ilev) = tmid(icol,ilev);
        pmid_day(iday,ilev) = pmid(icol,ilev);
        pint_day(iday,ilev) = pint(icol,ilev);
        albedo_dir_day(igpt,iday) = albedo_dir(igpt,icol);
        albedo_dif_day(igpt,iday) = albedo_dif(igpt,icol);
        coszrs_day(iday) = coszrs(icol);
        gas_vmr_day(igpt,iday,ilev) = gas_vmr(igpt,icol,ilev);
        cld_tau_gpt_day(iday,ilev,igpt) = cld_tau_gpt(icol,ilev,igpt);
        cld_ssa_gpt_day(iday,ilev,igpt) = cld_ssa_gpt(icol,ilev,igpt);
        cld_asm_gpt_day(iday,ilev,igpt) = cld_asm_gpt(icol,ilev,igpt);
    });

    parallel_for(SimpleBounds<3>(num_day(1), nlev, nswbands), YAKL_LAMBDA (int iday, int ilev, int ibnd)
    {
        auto icol = day_indices(iday);
        aer_tau_bnd_day(iday,ilev,ibnd) = aer_tau_bnd(icol,ilev,ibnd);
        aer_ssa_bnd_day(iday,ilev,ibnd) = aer_ssa_bnd(icol,ilev,ibnd);
        aer_asm_bnd_day(iday,ilev,ibnd) = aer_asm_bnd(icol,ilev,ibnd);
    });

    // Allocate shortwave fluxes (allsky and clearsky)
    // NOTE: fluxes defined at interfaces, so initialize to have vertical
    // dimension nlev_rad+1, while we initialized the RRTMGP input variables to
    // have vertical dimension nlev_rad (defined at midpoints).
    FluxesByband fluxes_clrsky_day, fluxes_allsky_day;
    internal::initial_fluxes(num_day(1), nlev+1, nswbands, fluxes_allsky_day);
    internal::initial_fluxes(num_day(1), nlev+1, nswbands, fluxes_clrsky_day);

    // Add an empty level above model top
    // TODO: combine with day compression above
    yakl::memset(cld_tau_gpt_rad, 0.);
    yakl::memset(cld_ssa_gpt_rad, 0.);
    yakl::memset(cld_asm_gpt_rad, 0.);

    yakl::memset(aer_tau_bnd_rad, 0.);
    yakl::memset(aer_ssa_bnd_rad, 0);
    yakl::memset(aer_asm_bnd_rad, 0.);

    parallel_for(SimpleBounds<3>(num_day(1), nlev, nswgpts), YAKL_LAMBDA (int iday, int ilev, int igpt)
    {
        cld_tau_gpt_rad(iday,ilev,igpt) = cld_tau_gpt_day(iday,ilev,igpt);
        cld_ssa_gpt_rad(iday,ilev,igpt) = cld_ssa_gpt_day(iday,ilev,igpt);
        cld_asm_gpt_rad(iday,ilev,igpt) = cld_asm_gpt_day(iday,ilev,igpt);
        gas_vmr_rad(igpt,iday,1) = gas_vmr_day(igpt,iday,1);
        gas_vmr_rad(igpt,iday,ilev) = gas_vmr_day(igpt,iday,ilev);
    });

    parallel_for(SimpleBounds<3>(num_day(1), nlev, nswbands), YAKL_LAMBDA (int iday, int ilev, int ibnd)
    {
        aer_tau_bnd_rad(iday,ilev,ibnd) = aer_tau_bnd_day(iday,ilev,ibnd);
        aer_ssa_bnd_rad(iday,ilev,ibnd) = aer_ssa_bnd_day(iday,ilev,ibnd);
        aer_asm_bnd_rad(iday,ilev,ibnd) = aer_asm_bnd_day(iday,ilev,ibnd);
    });

    // Do shortwave radiative transfer calculations
    radiation.run_shortwave_rrtmgp(ngas, num_day(1), nlev, gas_vmr_rad, pmid,
                                   tmid_day, pint_day, coszrs_day, albedo_dir_day, albedo_dif_day,
                                   cld_tau_gpt_rad, cld_ssa_gpt_rad, cld_asm_gpt_rad, aer_tau_bnd_rad, aer_ssa_bnd_rad, aer_asm_bnd_rad,
                                   fluxes_allsky_day.flux_up    , fluxes_allsky_day.flux_dn    , fluxes_allsky_day.flux_net    , fluxes_allsky_day.flux_dn_dir    ,
                                   fluxes_allsky_day.bnd_flux_up, fluxes_allsky_day.bnd_flux_dn, fluxes_allsky_day.bnd_flux_net, fluxes_allsky_day.bnd_flux_dn_dir,
                                   fluxes_clrsky_day.flux_up    , fluxes_clrsky_day.flux_dn    , fluxes_clrsky_day.flux_net    , fluxes_clrsky_day.flux_dn_dir    ,
                                   fluxes_clrsky_day.bnd_flux_up, fluxes_clrsky_day.bnd_flux_dn, fluxes_clrsky_day.bnd_flux_net, fluxes_clrsky_day.bnd_flux_dn_dir,
                                   tsi_scaling);

    // Expand fluxes from daytime-only arrays to full chunk arrays
    internal::expand_day_fluxes(fluxes_allsky_day, fluxes_allsky, day_indices);
    internal::expand_day_fluxes(fluxes_clrsky_day, fluxes_clrsky, day_indices);

    // Calculate heating rates
    calculate_heating_rate(fluxes_allsky.flux_up,
                           fluxes_allsky.flux_dn,
                           pint, qrs);

    calculate_heating_rate(fluxes_clrsky.flux_up,
                           fluxes_allsky.flux_dn,
                           pint, qrsc);
}

void Radiation::radiation_driver_lw (int ncol, int nlev,
                                     const real3d& gas_vmr,
                                     const real2d& pmid, const real2d& pint, const real2d& tmid, const real2d& tint,
                                     const real3d& cld_tau_gpt, const real3d& aer_tau_bnd, FluxesByband& fluxes_clrsky,
                                     FluxesByband& fluxes_allsky, const real2d& qrl, const real2d& qrlc)
{
    real3d cld_tau_gpt_rad("cld_tau_gpt_rad", ncol, nlev+1, nlwgpts);
    real3d aer_tau_bnd_rad("aer_tau_bnd-rad", ncol, nlev+1, nlwgpts);

    // Surface emissivity needed for longwave
    real2d surface_emissivity("surface_emissivity", nlwbands, ncol);

    // Temporary heating rates on radiation vertical grid
    real2d qrl_rad("qrl_rad", ncol, nlev);
    real2d qrlc_rad("qrlc_rad", ncol, nlev);
    real3d gas_vmr_rad("gas_vmr_rad", ngas, ncol, nlev);

    // Set surface emissivity to 1 here. There is a note in the RRTMG
    // implementation that this is treated in the land model, but the old
    // RRTMG implementation also sets this to 1. This probably does not make
    // a lot of difference either way, but if a more intelligent value
    // exists or is assumed in the model we should use it here as well.
    // TODO: set this more intelligently?
    yakl::memset(surface_emissivity, 1.0);

    // Add an empty level above model top
    yakl::memset(cld_tau_gpt_rad, 0.);
    yakl::memset(aer_tau_bnd_rad, 0.);

    parallel_for(SimpleBounds<3>(ncol, nlev, nswgpts), YAKL_LAMBDA (int icol, int ilev, int igpt)
    {
        cld_tau_gpt_rad(icol,ilev,igpt) = cld_tau_gpt(icol,ilev,igpt);
        aer_tau_bnd_rad(icol,ilev,igpt) = aer_tau_bnd(icol,ilev,igpt);
        gas_vmr_rad(igpt,icol,ilev) = gas_vmr(igpt,icol,ilev);
    });

    // Do longwave radiative transfer calculations
    radiation.run_longwave_rrtmgp(ngas, ncol, nlev,
                                  gas_vmr_rad, pmid, tmid, pint, tint,
                                  surface_emissivity, cld_tau_gpt_rad, aer_tau_bnd_rad,
                                  fluxes_allsky.flux_up    , fluxes_allsky.flux_dn    , fluxes_allsky.flux_net    ,
                                  fluxes_allsky.bnd_flux_up, fluxes_allsky.bnd_flux_dn, fluxes_allsky.bnd_flux_net,
                                  fluxes_clrsky.flux_up    , fluxes_clrsky.flux_dn    , fluxes_clrsky.flux_net    ,
                                  fluxes_clrsky.bnd_flux_up, fluxes_clrsky.bnd_flux_dn, fluxes_clrsky.bnd_flux_net);

    // Calculate heating rates
    calculate_heating_rate(fluxes_allsky.flux_up,
                           fluxes_allsky.flux_dn,
                           pint, qrl_rad);

    calculate_heating_rate(fluxes_allsky.flux_up,
                           fluxes_allsky.flux_dn,
                           pint,qrlc_rad);

    // Map heating rates to CAM columns and levels
    parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
    {
        qrl(icol,ilev)  = qrl_rad(icol,ilev);
        qrlc(icol,ilev) = qrlc_rad(icol,ilev);
    });
}

// Initialize array of daytime indices to be all zero. If any zeros exist when
// we are done, something went wrong.
void Radiation::set_daynight_indices (const real1d& coszrs, const int1d& day_indices, const int1d& night_indices)
{
    // Loop over columns and identify daytime columns as those where the cosine
    // solar zenith angle exceeds zero. Note that we wrap the setting of
    // day_indices in an if-then to make sure we are not accesing day_indices out
    // of bounds, and stopping with an informative error message if we do for some reason.
    int1d iday("iday", 1);
    int1d inight("inight",1);
    yakl::memset(iday, 0);
    yakl::memset(inight, 0);
    yakl::memset(day_indices, 0);
    yakl::memset(night_indices, 0);
    parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA (int icol)
    {
        if (coszrs(icol) > 0.) {
            iday(1) += 1;
            day_indices(iday(1)) = icol;
        } else {
            inight(1) += 1;
            night_indices(inight(1)) = icol;
        }
    });
}

void Radiation::get_gas_vmr (const std::vector<std::string>& gas_names, const real3d& gas_vmr)
{
    // Mass mixing ratio
    real2d mmr("mmr", ncol, nlev);
    // Gases and molecular weights. Note that we do NOT have CFCs yet (I think
    // this is coming soon in RRTMGP). RRTMGP also allows for absorption due to
    // CO and N2, which RRTMG did not have.
    const std::vector<std::string> gas_species = {"H2O", "CO2", "O3", "N2O",
                                                  "CO", "CH4", "O2", "N2"};
    const std::vector<real> mol_weight_gas = {18.01528, 44.0095, 47.9982, 44.0128,
                                              28.0101, 16.04246, 31.998, 28.0134}; // g/mol
    // Molar weight of air
    //const real mol_weight_air = 28.97; // g/mol
    // Defaults for gases that are not available (TODO: is this still accurate?)
    const real co_vol_mix_ratio = 1.0e-7;
    const real n2_vol_mix_ratio = 0.7906;

    // initialize
    yakl::memset(gas_vmr, 0.);

    // For each gas species needed for RRTMGP, read the mass mixing ratio from the
    // CAM rad_constituents interface, convert to volume mixing ratios, and
    // subset for daytime-only indices if needed.
    for (auto igas = 0; igas < gas_names.size(); ++igas) {

        std::cout << "gas_name: " << gas_names[igas] << "; igas: " << igas << std::endl;

        if (gas_names[igas] == "CO"){
            // CO not available, use default
            parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
            {
                gas_vmr(igas+1,icol,ilev) = co_vol_mix_ratio;
            });
        } else if (gas_names[igas] == "N2") {
            // N2 not available, use default
            parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
            {
                gas_vmr(igas+1,icol,ilev) = n2_vol_mix_ratio;
            });
        } else if (gas_names[igas] == "H2O") {
            // Water vapor is represented as specific humidity in CAM, so we
            // need to handle water a little differently
            //        rad_cnst_get_gas(icall, gas_species[igas], mmr);

            // Convert to volume mixing ratio by multiplying by the ratio of
            // molecular weight of dry air to molecular weight of gas. Note that
            // first specific humidity (held in the mass_mix_ratio array read
            // from rad_constituents) is converted to an actual mass mixing ratio.
            parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
            {
                gas_vmr(igas+1,icol,ilev) = qt(icol,ilev); //mmr(icol,ilev) / (
                //                  1. - mmr(icol,ilev))*mol_weight_air / mol_weight_gas[igas];
            });
        } else if (gas_names[igas] == "CO2") {
            parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
            {
                gas_vmr(igas+1,icol,ilev) = 3.8868676125307193E-4;
            });
        } else if (gas_names[igas] == "O3") {
            parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
            {
                gas_vmr(igas+1,icol,ilev) = 1.8868676125307193E-7;
            });
        } else if (gas_names[igas] == "N2O") {
            parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
            {
                gas_vmr(igas+1,icol,ilev) = 3.8868676125307193E-7;
            });
        } else if (gas_names[igas] == "CH4") {
            parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
            {
                gas_vmr(igas+1,icol,ilev) = 1.8868676125307193E-6;
            });
        } else if (gas_names[igas] == "O2") {
            parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
            {
                gas_vmr(igas+1,icol,ilev) = 0.2095;
            });
        } else {
            // Get mass mixing ratio from the rad_constituents interface
            //        rad_cnst_get_gas(icall, gas_species[igas], mmr);

            // Convert to volume mixing ratio by multiplying by the ratio of
            // molecular weight of dry air to molecular weight of gas
            parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
            {
                gas_vmr(igas+1,icol,ilev) = 1.0e-6; //mmr(icol,ilev)
                //                                     * mol_weight_air / mol_weight_gas[igas];
            });
        }
    } // igas
}

// Loop over levels and calculate heating rates; note that the fluxes *should*
// be defined at interfaces, so the loop ktop,kbot and grabbing the current
// and next value of k should be safe. ktop should be the top interface, and
// kbot + 1 should be the bottom interface.
// NOTE: to get heating rate in K/day, normally we would use:
//     H = dF / dp * g * (sec/day) * (1e-5) / (cpair)
// Here we just use
//     H = dF / dp * g
// Why? Something to do with convenience with applying the fluxes to the
// heating tendency?
void Radiation::calculate_heating_rate (const real2d& flux_up,
                                        const real2d& flux_dn,
                                        const real2d& pint,
                                        const real2d& heating_rate)
{
    parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
    {
        heating_rate(icol,ilev) = (flux_up(icol,ilev+1) - flux_up(icol,ilev)- flux_dn(icol,ilev+1)+flux_dn(icol,ilev))
            *CONST_GRAV/(pint(icol,ilev+1)-pint(icol,ilev));
    });
}

// call back
void Radiation::on_complete () { }


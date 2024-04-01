//
// get the cloud radiation props data
//
#include <cmath>

#include "ERF_Constants.H"
#include "Aero_rad_props.H"
#include "Water_vapor_saturation.H"
using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

void AerRadProps::initialize (int num_gas, int num_modes, int naeroes,
                              int nswbands_, int nlwbands_,
                              int ncoloum, int nlevel, int num_rh, int top_levels,
                              const std::vector<std::string>& aerosol_names,
                              const real2d& zint, const real2d& pmiddle, const real2d& temperature,
                              const real2d& qtotal, const real2d& geom_rad)
{
    nmodes = num_modes;
    ngas = num_gas;
    num_aeroes = naeroes;
    aero_names = aerosol_names;

    nswbands = nswbands_;
    nlwbands = nlwbands_;
    ncol     = ncoloum;
    nlev     = nlevel;
    nrh      = num_rh;
    top_lev  = top_levels;

    pmid = pmiddle;
    pdeldry = pmiddle;
    temp = temperature;
    qt   = qtotal;
    // geometric radius
    geometric_radius = geom_rad;
    // vertical grid
    zi = zint;
    // initialize mam4 aero model
    mam_aer.initialize(ncol, nlev, top_lev, nswbands, nlwbands);
}

void AerRadProps::aer_rad_props_sw (const int& list_idx, const real& dt, const int& nnite,
                                    const int1d& idxnite, const bool is_cmip6_volc, const real3d& tau, const real3d& tau_w,
                                    const real3d& tau_w_g, const real3d& tau_w_f, const real2d& clear_rh)
{
    // for cmip6 style volcanic file
    int1d trop_level("trop_level", ncol);
    real3d ext_cmip6_sw_inv_m("ext_cmip6_sw_inv_m", ncol, nlev, nswbands); // short wave extinction in the units of 1/m

    // compute mixing ratio to mass conversion
    real2d mmr_to_mass("mmr_to_mass", ncol, nlev); // conversion factor for mmr to mass
    parallel_for(SimpleBounds<2>(ncol, nlev) , YAKL_LAMBDA (int icol, int ilev)
    {
        mmr_to_mass(icol,ilev) = rgas * pdeldry(icol,ilev);
    });

    // initialize to conditions that would cause failure
    yakl::memset(tau, -100.);
    yakl::memset(tau_w, -100.);
    yakl::memset(tau_w_g, -100.);
    yakl::memset(tau_w_f, -100.);
    yakl::memset(ext_cmip6_sw_inv_m, 1.0e-3);

    // top layer (ilev = 0) has no aerosol (ie tau = 0)
    // also initialize rest of layers to accumulate od's
    parallel_for(SimpleBounds<2>(ncol, nswbands) , YAKL_LAMBDA (int icol, int isw)
    {
        tau    (icol,1,isw) = 0.;
        tau_w  (icol,1,isw) = 0.;
        tau_w_g(icol,1,isw) = 0.;
        tau_w_f(icol,1,isw) = 0.;
    });

    // local variables
    real2d wrh("wrh", ncol, nlev);
    int2d  krh("krh", ncol, nlev);
    // calculate relative humidity for table lookup into rh grid
    parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
    {
        real es, qs, rh, rhtrunc;
        WaterVaporSat::qsat(temp(icol,ilev), pmid(icol,ilev), es, qs);
        rh = qt(icol,ilev)/qs;
        rhtrunc = std::min(rh,1.);
        krh(icol, ilev) = std::min(std::floor(rhtrunc*nrh )+1, nrh - 1.); // index into rh mesh
        wrh(icol, ilev) = rhtrunc*nrh-krh(icol, ilev);       // (-) weighting on left side values
    });

    // Special treatment for CMIP6 volcanic aerosols, where extinction, ssa
    // and af are directly read from the prescribed volcanic aerosol file

    // Point ext_cmip6_sw to null so that the code breaks if they are used for is_cmip6_volc=.false. case
    // This is done to avoid having optional arguments in modal_aero_sw call
    yakl::memset(trop_level, 10000);

    // I was thinking whether we have to include volcanic effects on radiation? (xyuan)
    if (is_cmip6_volc) {
        // get extinction so as to supply to modal_aero_sw routine for computing EXTINCT variable
        // converting it from 1/km to 1/m
        //      call pbuf_get_field(pbuf, idx_ext_sw, ext_cmip6_sw)
        //      call outfld('extinct_sw_inp',ext_cmip6_sw(:,:,idx_sw_diag), pcols, lchnk)
        parallel_for(SimpleBounds<3>(ncol, nlev, nswbands), YAKL_LAMBDA (int icol, int ilev, int isw)
        {
            ext_cmip6_sw_inv_m(icol, ilev, isw) = 1.0e-3; //*ext_cmip6_sw(icol, ilev, isw); //convert from 1/km to 1/m
        });

        //Find tropopause as extinction should be applied only above tropopause
        //trop_level has value for tropopause for each column
        //tropopause_find(state, trop_level)
        // Iterate over all of the columns to find the troppause
        real1d strop("strop",1);
        parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA (int i)
        {
            // Search from the bottom up, looking for the first minimum.
            int tlev  = -1;
            for (auto k = nlev-1; k > 0; --k) {
                real stobie = 0.03 * temp(i,k) - std::log10(pmid(i, k));
                if (pmid(i, k) <= 4000.) break;
                if (pmid(i, k) >= 55000.) continue;
                if ((tlev == -1) || (stobie < strop(1))) {
                    tlev  = k;
                    strop(1) = stobie;
                }
            }
            if (tlev != -1) trop_level(i) = tlev;
        });

        //
        //Quit if tropopause is not found
        if (yakl::intrinsics::any(trop_level) == -1) {
            amrex::Print() << "aer_rad_props.F90: subr aer_rad_props_sw: tropopause not found\n";
        }
    }

    // get number of bulk aerosols and number of modes in current list
    mam_consti.get_nmodes(list_idx, nmodes);
    mam_consti.get_ngas(list_idx, ngas);
    mam_consti.get_naero(list_idx, num_aeroes);

    // Contributions from modal aerosols.
    if (nmodes > 0) {
        real2d ext_cmip6_sw_2d("ext_cmip6_sw",ncol,nlev);
        parallel_for(SimpleBounds<3>(ncol, nlev, nswbands), YAKL_LAMBDA (int icol, int ilev, int isw)
        {
            ext_cmip6_sw_2d(icol,ilev) = ext_cmip6_sw_inv_m(icol,ilev,RadConstants::idx_sw_diag);
        });

        mam_aer.modal_aero_sw(list_idx, dt, nnite, idxnite, is_cmip6_volc,
                              pdeldry, pmid, temp, qt, ext_cmip6_sw_2d,
                              trop_level, tau, tau_w, tau_w_g, tau_w_f, clear_rh);
    } else {
        yakl::memset(tau, 0.);
        yakl::memset(tau_w, 0.);
        yakl::memset(tau_w_g, 0.);
        yakl::memset(tau_w_f, 0.);
    }

#if 1
    if (is_cmip6_volc) {
        //update tau, tau_w, tau_w_g, and tau_w_f with the read in values of extinction, ssa and asymmetry factors
        volcanic_cmip_sw (trop_level, zi, ext_cmip6_sw_inv_m, ssa_cmip6_sw, af_cmip6_sw,
                          tau, tau_w, tau_w_g, tau_w_f);
    }

    // Contributions from bulk aerosols.
    for (auto iaero = 0; iaero < num_aeroes; ++iaero) {
        // aerosol masses
        real2d aermmr("aermmr", ncol, nlev);    // mass mixing ratio of aerosols
        real2d aermass("aermass", ncol, nlev);     // mass of aerosols

        yakl::memset(aermmr, 1.0e-2);
        // get bulk aerosol mass mixing ratio
        //      mam_consti.rad_cnst_get_aer_mmr(list_idx, iaero, aermmr);
        parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
        {
            if (ilev < top_lev) {
                aermass(icol,ilev) = 0.;
            } else {
                aermass(icol,ilev) = aermmr(icol,ilev) * mmr_to_mass(icol,ilev);
            }
        });

        std::string opticstype;       // hygro or nonhygro
        // get optics type
        mam_consti.get_aer_opticstype(list_idx, iaero, opticstype);

        // radiative properties for each aerosol
        real3d ta("ta", ncol, nlev, nswbands);
        real3d tw("tw", ncol, nlev, nswbands);
        real3d twf("twf", ncol, nlev, nswbands);
        real3d twg("twg", ncol, nlev, nswbands);

        if (opticstype == "hygro") {
            // get optical properties for hygroscopic aerosols
            // optical props for each aerosol hygroscopic
            real2d h_ext("h_ext", ncol, nlev);
            real2d h_ssa("h_ssa", ncol, nlev);
            real2d h_asm("h_asm", ncol, nlev);
            mam_consti.get_aer_sw_hygro_ext(list_idx, iaero, h_ext);
            mam_consti.get_aer_sw_hygro_ssa(list_idx, iaero, h_ssa);
            mam_consti.get_aer_sw_hygro_asm(list_idx, iaero, h_asm);
            get_hygro_rad_props(ncol, krh, wrh, aermass, h_ext, h_ssa, h_asm, ta, tw, twg, twf);
            parallel_for(SimpleBounds<3>(ncol, nlev, nswbands), YAKL_LAMBDA (int icol, int ilev, int isw)
            {
                tau    (icol,ilev,isw) = tau    (icol,ilev,isw) + ta (icol,ilev,isw);
                tau_w  (icol,ilev,isw) = tau_w  (icol,ilev,isw) + tw (icol,ilev,isw);
                tau_w_g(icol,ilev,isw) = tau_w_g(icol,ilev,isw) + twg(icol,ilev,isw);
                tau_w_f(icol,ilev,isw) = tau_w_f(icol,ilev,isw) + twf(icol,ilev,isw);
            });
        } else if (opticstype == "nonhygro") {
            // get optical properties for non-hygroscopic aerosols
            // non-hygroscopic
            real1d n_ext("n_ext", ncol);
            real1d n_ssa("n_ssa", ncol);
            real1d n_asm("n_asm", ncol);
            mam_consti.get_aer_sw_nonhygro_ext(list_idx, iaero, n_ext);
            mam_consti.get_aer_sw_nonhygro_ssa(list_idx, iaero, n_ssa);
            mam_consti.get_aer_sw_nonhygro_asm(list_idx, iaero, n_asm);
            get_nonhygro_rad_props(ncol, aermass, n_ext, n_ssa, n_asm, ta, tw, twg, twf);
            parallel_for(SimpleBounds<3>(ncol, nlev, nswbands), YAKL_LAMBDA (int icol, int ilev, int isw)
            {
                tau    (icol,ilev,isw) = tau    (icol,ilev,isw) + ta (icol,ilev,isw);
                tau_w  (icol,ilev,isw) = tau_w  (icol,ilev,isw) + tw (icol,ilev,isw);
                tau_w_g(icol,ilev,isw) = tau_w_g(icol,ilev,isw) + twg(icol,ilev,isw);
                tau_w_f(icol,ilev,isw) = tau_w_f(icol,ilev,isw) + twf(icol,ilev,isw);
            });
        } else if (opticstype == "volcanic") {
            // get optical properties for volcanic aerosols
            real1d n_ext("n_ext", ncol);
            real1d n_scat("n_scat", ncol);
            real1d n_ascat("n_ascat", ncol);
            mam_consti.get_aer_sw_nonhygro_ext(list_idx, iaero, n_ext);
            mam_consti.get_aer_sw_nonhygro_scat(list_idx, iaero, n_scat);
            mam_consti.get_aer_sw_nonhygro_ascat(list_idx, iaero, n_ascat);
            get_volcanic_rad_props(ncol, aermass, n_ext, n_scat, n_ascat, ta, tw, twg, twf);
            parallel_for(SimpleBounds<3>(ncol, nlev, nswbands), YAKL_LAMBDA (int icol, int ilev, int isw)
            {
                tau    (icol,ilev,isw) = tau    (icol,ilev,isw) + ta (icol,ilev,isw);
                tau_w  (icol,ilev,isw) = tau_w  (icol,ilev,isw) + tw (icol,ilev,isw);
                tau_w_g(icol,ilev,isw) = tau_w_g(icol,ilev,isw) + twg(icol,ilev,isw);
                tau_w_f(icol,ilev,isw) = tau_w_f(icol,ilev,isw) + twf(icol,ilev,isw);
            });
        } else if (opticstype == "volcanic_radius") {
            // get optical properties for volcanic aerosols
            real2d r_ext("r_ext",ncol,nlev);
            real2d r_scat("r_scat",ncol,nlev);
            real2d r_ascat("r_ascat",ncol,nlev);
            real1d r_mu("r_mu", ncol);
            mam_consti.get_aer_r_sw_ext(list_idx, iaero, r_ext);
            mam_consti.get_aer_r_sw_scat(list_idx, iaero, r_scat);
            mam_consti.get_aer_r_sw_ascat(list_idx, iaero, r_ascat);
            mam_consti.get_aer_mu(list_idx, iaero, r_mu);
            get_volcanic_radius_rad_props(ncol, aermass, r_ext, r_scat, r_ascat, r_mu, ta, tw, twg, twf);
            //get_volcanic_radius_rad_props(const real2d& mass,
            //                                        const real2d& r_ext,
            //                                       const real2d& r_scat,
            //                                      const real2d& r_ascat,
            //                                      const real1d& r_mu,
            //                                      real3d& tau,
            //                                      real3d& tau_w,
            //                                      real3d& tau_w_g,
            //                                      real3d& tau_w_f)

            parallel_for(SimpleBounds<3>(ncol, nlev, nswbands), YAKL_LAMBDA (int icol, int ilev, int isw)
            {
                tau    (icol,ilev,isw) = tau    (icol,ilev,isw) + ta (icol,ilev,isw);
                tau_w  (icol,ilev,isw) = tau_w  (icol,ilev,isw) + tw (icol,ilev,isw);
                tau_w_g(icol,ilev,isw) = tau_w_g(icol,ilev,isw) + twg(icol,ilev,isw);
                tau_w_f(icol,ilev,isw) = tau_w_f(icol,ilev,isw) + twf(icol,ilev,isw);
            });
        } else {
            amrex::Print() << "aer_rad_props_sw: unsupported opticstype\n";
        }

        // diagnostic output of individual aerosol optical properties
        // currently implemented for climate list only
        //      call aer_vis_diag_out(lchnk, ncol, nnite, idxnite, iaero, ta(:,:,idx_sw_diag), list_idx)
    }

    // diagnostic output of total aerosol optical properties
    // currently implemented for climate list only
    //   call aer_vis_diag_out(lchnk, ncol, nnite, idxnite, 0, tau(:,:,idx_sw_diag), list_idx)
#endif
}

// Purpose: Compute aerosol transmissions needed in absorptivity/
//    emissivity calculations
//
// lw extinction is the same representation for all
// species.  If this changes, this routine will need to do something
// similar to the sw with routines like get_hygro_lw_abs
void AerRadProps::aer_rad_props_lw (const bool& is_cmip6_volc,
                                    const int& list_idx,
                                    const real& dt,
                                    const real2d& zi,
                                    const real3d& odap_aer,
                                    const real2d& clear_rh)
{
    int numaerosols; // number of bulk aerosols in climate/diagnostic list
    int nmodes;      // number of aerosol modes in climate/diagnostic list
    std::string opticstype;  // hygro or nonhygro

    // optical props for each aerosol
    real1d lw_abs("lw_abs", ncol);
    real2d lw_hygro_abs("lw_hygro", ncol, nlev);
    real2d geometric_radius("geometric_radius", ncol, nlev);

    // volcanic lookup table
    real2d r_lw_abs("r_lw_abs", ncol, nlev);  // radius dependent mass-specific absorption coefficient
    real1d r_mu("r_mu", ncol);        // log(geometric_mean_radius) domain samples of r_lw_abs(:,:)
    real2d mu("mu", ncol, nlev);          // log(geometric_radius)
    real r_mu_min, r_mu_max, wmu, mutrunc;
    int nmu, kmu;

    // for table lookup into rh grid
    real2d es("es", ncol, nlev);     // saturation vapor pressure
    real2d qs("qs", ncol, nlev);     // saturation specific humidity
    real2d rh("rh", ncol, nlev);
    real2d rhtrunc("rhtrunc", ncol, nlev);
    real2d wrh("wrh", ncol, nlev);
    int2d krh("krh", ncol, nlev);

    // aerosol (vertical) mass path and extinction aerosol masses
    real2d aermmr("aermmr", ncol, nlev);    // mass mixing ratio of aerosols
    real2d mmr_to_mass("mmr_to_mass", ncol, nlev); // conversion factor for mmr to mass
    real2d aermass("aermass", ncol, nlev);     // mass of aerosols

    //For cmip6 volcanic file
    int1d trop_level("trop_level", ncol);
    real lyr_thk;
    real3d ext_cmip6_lw("ext_cmip6_lw", ncol, nlev, nlwbands);
    real3d ext_cmip6_lw_inv_m("ext_cmip6_lw_inv_m",ncol, nlev, nlwbands); //long wave extinction in the units of 1/m

    const real km_inv_to_m_inv = 0.001;

    // get number of bulk aerosols and number of modes in current list
    mam_consti.get_nmodes(list_idx, nmodes);
    mam_consti.get_naero(list_idx, numaerosols);

    // Contributions from modal aerosols.
    if (nmodes > 0) {
        mam_aer.modal_aero_lw(list_idx, dt, pdeldry, pmid, temp, qt, odap_aer,clear_rh);
    } else {
        yakl::memset(odap_aer, 0.);
    }

    // Contributions from bulk aerosols.
    if (numaerosols > 0) {

        // compute mixing ratio to mass conversion
        parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
        {
            mmr_to_mass(icol,ilev) = rgas * pdeldry(icol,ilev);
        });

        // calculate relative humidity for table lookup into rh grid
        parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
        {
            WaterVaporSat::qsat(temp(icol,ilev), pmid(icol,ilev), es(icol,ilev), qs(icol,ilev));
            rh(icol,ilev) = qt(icol,ilev) / qs(icol,ilev);

            rhtrunc(icol, ilev) = std::min(rh(icol, ilev),1.);
            krh(icol, ilev) = std::min(std::floor(rhtrunc(icol, ilev) * nrh ) + 1, nrh - 1.);  // index into rh mesh
            wrh(icol, ilev) = rhtrunc(icol, ilev) * nrh - krh(icol, ilev);       // (-) weighting on left side values
        });
    }

    yakl::memset(ext_cmip6_lw, 0.);
    if(is_cmip6_volc) {
        // Logic:
        // Update odap_aer with the read in volcanic aerosol extinction (1/km).
        // It needs to be converted to 1/m and multiplied by layer thickness first.

        // Obtain read in values for ext from the volcanic input file
        //      call pbuf_get_field(pbuf, idx_ext_lw, ext_cmip6_lw)
        //      call outfld('extinct_lw_inp',ext_cmip6_lw(:,:,idx_lw_diag), pcols, lchnk)
        //      ext_cmip6_lw_inv_m = ext_cmip6_lw * km_inv_to_m_inv  !convert from 1/km to 1/m
        //
        // Above the tropopause, the read in values from the file include both the stratospheric
        // and volcanic aerosols. Therefore, we need to zero out odap_aer above the tropopause
        // and populate it exclusively from the read in values.

        // Find tropopause
        // trop_level has value for tropopause for each column
        //      tropopause_find(state, trop_level)
        real1d strop("strop",1);
        parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA (int i)
        {
            // Search from the bottom up, looking for the first minimum.
            int tlev  = -1;
            for (auto k = nlev-1; k > 0; --k) {
                real stobie = 0.03 * temp(i,k) - std::log10(pmid(i, k));
                if (pmid(i, k) <= 4000.) break;
                if (pmid(i, k) >= 55000.) continue;
                if ((tlev == -1) || (stobie < strop(0))) {
                    tlev  = k;
                    strop(0) = stobie;
                }
            }
            if (tlev != -1) trop_level(i) = tlev;
        });

        // Quit if tropopause is not found
        if (yakl::intrinsics::any(trop_level) == -1)
            amrex::Print() << "aer_rad_props_lw: tropopause not found\n";

        // If tropopause is found, update taus with 50% contributuions from the volcanic input
        // file and 50% from the existing model computed values
        // First handle the case of tropopause layer itself:
        parallel_for(SimpleBounds<2>(ncol, nlwbands), YAKL_LAMBDA (int icol, int ilw)
        {
            int ilev_tropp = trop_level(icol); //tropopause level
            auto lyr_thk    = zi(icol,ilev_tropp) - zi(icol,ilev_tropp+1); //in meters
            auto ext_cmip6_inv_m = km_inv_to_m_inv * ext_cmip6_lw_inv_m(icol,ilev_tropp,ilw);
            odap_aer(icol,ilev_tropp,ilw) = 0.5 * (odap_aer(icol,ilev_tropp,ilw) + lyr_thk * ext_cmip6_inv_m);
        });

        // As it will be more efficient for FORTRAN to loop over levels and then columns, the following loops
        // are nested keeping that in mind
        parallel_for(SimpleBounds<3>(ncol, nlev, nlwbands), YAKL_LAMBDA (int icol, int ilev, int ilw)
        {
            int ilev_tropp = trop_level(icol); //tropopause level
            if (ilev < ilev_tropp) {
                auto lyr_thk = zi(icol,ilev) - zi(icol,ilev+1);
                // odap_aer(icol,ilev,ilw) = lyr_thk * ext_cmip6_lw_inv_m(icol,ilev,ilw);
            }
        });
    }


    // Loop over bulk aerosols in list.
    for (auto iaerosol = 0; iaerosol < numaerosols; ++iaerosol) {

        // get aerosol mass mixing ratio
        //      mam_consti.rad_cnst_get_aer_mmr(list_idx, iaerosol, aermmr);
        parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
        {
            if (ilev < top_lev) {
                aermass(icol,ilev) = 0.;
            } else {
                aermass(icol,ilev) = aermmr(icol,ilev) * mmr_to_mass(icol,ilev);
            }
        });

        // get optics type
        mam_consti.get_aer_opticstype(list_idx, iaerosol, opticstype);

        if (opticstype == "hygroscopic") {
            // get optical properties for hygroscopic aerosols
            mam_consti.get_aer_lw_hygro_abs(list_idx, iaerosol, lw_hygro_abs);
            parallel_for(SimpleBounds<3>(ncol, nlev, nlwbands), YAKL_LAMBDA (int icol, int ilev, int bnd_idx)
            {
                odap_aer(icol, ilev, bnd_idx) = odap_aer(icol, ilev, bnd_idx) + aermass(icol, ilev) *
                    ((1 + wrh(icol,ilev)) * lw_hygro_abs(krh(icol, ilev)+1,bnd_idx)
                     - wrh(icol, ilev)  * lw_hygro_abs(krh(icol, ilev),  bnd_idx));
            });
        } else if (opticstype == "insoluble" ||
                   opticstype == "nonhygro"  ||
                   opticstype == "hygro"     ||
                   opticstype == "volcanic") {
            // get optical properties for hygroscopic aerosols
            mam_consti.get_aer_lw_abs(list_idx, iaerosol, lw_abs);
            parallel_for(SimpleBounds<3>(ncol, nlev, nlwbands), YAKL_LAMBDA (int icol, int ilev, int bnd_idx)
            {
                odap_aer(icol,ilev,bnd_idx) = odap_aer(icol,ilev,bnd_idx) + lw_abs(bnd_idx)*aermass(icol,ilev);
            });
        } else if (opticstype == "volcanic_radius") {
            // get optical properties for hygroscopic aerosols
            mam_consti.get_aer_r_lw_abs(list_idx, iaerosol, r_lw_abs);
            mam_consti.get_aer_mu(list_idx, iaerosol, r_mu);
            // get microphysical properties for volcanic aerosols
            //         idx = pbuf_get_index('VOLC_RAD_GEOM')
            //         call pbuf_get_field(pbuf, idx, geometric_radius )

            // interpolate in radius
            // caution: clip the table with no warning when outside bounds
            parallel_for(SimpleBounds<2>(ncol, nlev), YAKL_LAMBDA (int icol, int ilev)
            {
                auto nmu = r_mu.extent(0);
                auto r_mu_max = r_mu(nmu);
                auto r_mu_min = r_mu(1);

                if(geometric_radius(icol,ilev) > 0.) {
                    mu(icol,ilev) = std::log(geometric_radius(icol,ilev));
                } else {
                    mu(icol,ilev) = 0.;
                }
                auto mutrunc = std::max(std::min(mu(icol,ilev),r_mu_max),r_mu_min);
                auto kmu = std::max(std::min(1 + (mutrunc-r_mu_min)/(r_mu_max-r_mu_min)*(nmu-1),nmu-1.),1.);
                auto wmu = std::max(std::min( (mutrunc -r_mu(kmu)) / (r_mu(kmu+1) - r_mu(kmu)) ,1.),0.);
                for (auto bnd_idx = 0; bnd_idx < nlwbands; ++bnd_idx) {
                    odap_aer(icol,ilev,bnd_idx) = odap_aer(icol,ilev,bnd_idx) +
                        aermass(icol,ilev) * ((1. - wmu) * r_lw_abs(bnd_idx, kmu) +
                                              (wmu) * r_lw_abs(bnd_idx, kmu+1));
                }
            });
        }
    }
}

 void AerRadProps::get_hygro_rad_props (const int& ncol,
                                        const int2d& krh,
                                        const real2d& wrh,
                                        const real2d& mass,
                                        const real2d& extb,
                                        const real2d& ssab,
                                        const real2d& asmb,
                                        const real3d& tau,
                                        const real3d& tau_w,
                                        const real3d& tau_w_g,
                                        const real3d& tau_w_f)
{
   parallel_for(SimpleBounds<3>(ncol,nlev,nswbands) , YAKL_LAMBDA (int icol, int ilev, int iswband)
   {
       auto ext1 = (1 + wrh(icol,ilev)) * extb(krh(icol,ilev)+1,iswband) - wrh(icol,ilev) * extb(krh(icol,ilev),iswband);
       auto ssa1 = (1 + wrh(icol,ilev)) * ssab(krh(icol,ilev)+1,iswband) - wrh(icol,ilev) * ssab(krh(icol,ilev),iswband);
       auto asm1 = (1 + wrh(icol,ilev)) * asmb(krh(icol,ilev)+1,iswband) - wrh(icol,ilev) * asmb(krh(icol,ilev),iswband);

       tau    (icol, ilev, iswband) = mass(icol, ilev) * ext1;
       tau_w  (icol, ilev, iswband) = mass(icol, ilev) * ext1 * ssa1;
       tau_w_g(icol, ilev, iswband) = mass(icol, ilev) * ext1 * ssa1 * asm1;
       tau_w_f(icol, ilev, iswband) = mass(icol, ilev) * ext1 * ssa1 * asm1 * asm1;
   });
}

void AerRadProps::get_nonhygro_rad_props (const int& ncol,
                                          const real2d& mass,
                                          const real1d& extb,
                                          const real1d& ssab,
                                          const real1d& asmb,
                                          const real3d& tau,
                                          const real3d& tau_w,
                                          const real3d& tau_w_g,
                                          const real3d& tau_w_f)
{
    parallel_for(SimpleBounds<3>(ncol,nlev,nswbands) , YAKL_LAMBDA (int icol, int ilev, int iswband)
    {
        auto ext1 = extb(iswband);
        auto ssa1 = ssab(iswband);
        auto asm1 = asmb(iswband);
        tau    (icol,ilev,iswband) = mass(icol,ilev) * ext1;
        tau_w  (icol,ilev,iswband) = mass(icol,ilev) * ext1 * ssa1;
        tau_w_g(icol,ilev,iswband) = mass(icol,ilev) * ext1 * ssa1 * asm1;
        tau_w_f(icol,ilev,iswband) = mass(icol,ilev) * ext1 * ssa1 * asm1 * asm1;
    });
}

void AerRadProps::get_volcanic_radius_rad_props (const int& ncol,
                                                 const real2d& mass,
                                                 const real2d& r_ext,
                                                 const real2d& r_scat,
                                                 const real2d& r_ascat,
                                                 const real1d& r_mu,
                                                 const real3d& tau,
                                                 const real3d& tau_w,
                                                 const real3d& tau_w_g,
                                                 const real3d& tau_w_f)
{
    yakl::memset(tau, 0.);
    yakl::memset(tau_w, 0.);
    yakl::memset(tau_w_g, 0.);
    yakl::memset(tau_w_f, 0.);

    // get microphysical properties for volcanic aerosols
    //idx = pbuf_get_index('VOLC_RAD_GEOM')
    //call pbuf_get_field(pbuf, idx, geometric_radius )

    // interpolate in radius
    // caution: clip the table with no warning when outside bounds
    auto nmu = r_mu.extent(0);
    auto r_mu_max = r_mu(nmu-1);
    auto r_mu_min = r_mu(0);

    // interpolated values from table
    real1d ext("ext", nswbands);
    real1d scat("scat", nswbands);
    real1d ascat("ascat", nswbands);
    real2d mu("mu", ncol, nlev);

    parallel_for(SimpleBounds<3>(ncol,nlev,nswbands), YAKL_LAMBDA (int icol, int ilev, int iswband)
    {
        if(geometric_radius(icol,ilev) > 0.)
            mu(icol,ilev) = std::log(geometric_radius(icol,ilev));
        else
            mu(icol,ilev) = 0.;

        auto mutrunc = std::max(std::min(mu(icol,ilev),r_mu_max),r_mu_min);
        auto kmu = std::max(std::min(1 + (mutrunc-r_mu_min)/(r_mu_max-r_mu_min)*(nmu-1),nmu-1.),1.);
        auto wmu = std::max(std::min((mutrunc -r_mu(kmu)) / (r_mu(kmu+1) - r_mu(kmu)),1.),0.);

        ext(iswband) = ((1. - wmu) * r_ext(iswband, kmu  ) +
                        (wmu) * r_ext(iswband, kmu+1));
        scat(iswband) = ((1. - wmu) * r_scat(iswband, kmu  ) +
                         (wmu) * r_scat(iswband, kmu+1));
        ascat(iswband) = ((1. - wmu) * r_ascat(iswband, kmu  ) +
                          (wmu) * r_ascat(iswband, kmu+1));
        real g = 0;
        if (scat(iswband) > 0.)
            g = ascat(iswband)/scat(iswband);

        tau    (icol,ilev,iswband) = mass(icol,ilev) * ext(iswband);
        tau_w  (icol,ilev,iswband) = mass(icol,ilev) * scat(iswband);
        tau_w_g(icol,ilev,iswband) = mass(icol,ilev) * ascat(iswband);
        tau_w_f(icol,ilev,iswband) = mass(icol,ilev) * g * ascat(iswband);
    });
}

 // Logic:
 // Update taus, tau_w, tau_w_g and tau_w_f with the read in volcanic
 // aerosol extinction (1/km), single scattering albedo and asymmtry factors.

 void AerRadProps::volcanic_cmip_sw (const int1d& trop_level,
                                     const real2d& zi,           // zone interface
                                     const real3d& ext_cmip6_sw_inv_m,
                                     const real3d& ssa_cmip6_sw, // ssa from the volcanic inputs
                                     const real3d& af_cmip6_sw,  // asymmetry factor (af) from volcanic inputs
                                     const real3d& tau,
                                     const real3d& tau_w,
                                     const real3d& tau_w_g,
                                     const real3d& tau_w_f) {
  // Above the tropopause, the read in values from the file include both the stratospheric
  // and volcanic aerosols. Therefore, we need to zero out taus above the tropopause
  // and populate them exclusively from the read in values.

  // If tropopause is found, update taus with 50% contributuions from the volcanic input
  // file and 50% from the existing model computed values
  real1d ext_unitless("ext_unitless", nswbands),
         asym_unitless("asym_unitless",nswbands),
         ext_ssa("ext_ssa",nswbands),
         ext_ssa_asym("ext_ssa_asym",nswbands);

  // First handle the case of tropopause layer itself:
  parallel_for(SimpleBounds<2>(ncol,nswbands), YAKL_LAMBDA (int icol, int isw) {
     auto ilev_tropp = trop_level(icol); //tropopause level
     auto lyr_thk = zi(icol,ilev_tropp) - zi(icol,ilev_tropp+1);

     ext_unitless(isw)  = lyr_thk * ext_cmip6_sw_inv_m(icol,ilev_tropp,isw);
     asym_unitless(isw) = af_cmip6_sw (icol,ilev_tropp,isw);
     ext_ssa(isw)       = ext_unitless(isw) * ssa_cmip6_sw(icol,ilev_tropp,isw);
     ext_ssa_asym(isw)  = ext_ssa(isw) * asym_unitless(isw);

     tau    (icol,ilev_tropp,isw) = 0.5 * ( tau    (icol,ilev_tropp,isw) + ext_unitless(isw) );
     tau_w  (icol,ilev_tropp,isw) = 0.5 * ( tau_w  (icol,ilev_tropp,isw) + ext_ssa(isw));
     tau_w_g(icol,ilev_tropp,isw) = 0.5 * ( tau_w_g(icol,ilev_tropp,isw) + ext_ssa_asym(isw));
     tau_w_f(icol,ilev_tropp,isw) = 0.5 * ( tau_w_f(icol,ilev_tropp,isw) + ext_ssa_asym(isw) * asym_unitless(isw));
  });

  // As it will be more efficient for FORTRAN to loop over levels and then columns, the following loops
  // are nested keeping that in mind
  parallel_for(SimpleBounds<3>(ncol,nlev, nswbands), YAKL_LAMBDA (int icol, int ilev, int isw) {
     auto ilev_tropp = trop_level(icol); //tropopause level

     if (ilev < ilev_tropp) { //BALLI: see if this is right!
        auto lyr_thk = zi(icol,ilev) - zi(icol,ilev+1);

        ext_unitless(isw)  = lyr_thk * ext_cmip6_sw_inv_m(icol,ilev,isw);
        asym_unitless(isw) = af_cmip6_sw(icol,ilev,isw);
        ext_ssa(isw)       = ext_unitless(isw) * ssa_cmip6_sw(icol,ilev,isw);
        ext_ssa_asym(isw)  = ext_ssa(isw) * asym_unitless(isw);

        tau    (icol,ilev,isw) = ext_unitless(isw);
        tau_w  (icol,ilev,isw) = ext_ssa(isw);
        tau_w_g(icol,ilev,isw) = ext_ssa_asym(isw);
        tau_w_f(icol,ilev,isw) = ext_ssa_asym(isw) * asym_unitless(isw);
     }
  });
}

void AerRadProps::get_volcanic_rad_props(const int& ncol,
                                         const real2d& mass,
                                         const real1d& ext,
                                         const real1d& scat,
                                         const real1d& ascat,
                                         const real3d& tau,
                                         const real3d& tau_w,
                                         const real3d& tau_w_g,
                                         const real3d& tau_w_f) {
   int nswbands;
   parallel_for(SimpleBounds<3>(nswbands, ncol,nlev), YAKL_LAMBDA (int iswband, int icol, int ilev) {
     real g = 0;
     if (scat(iswband) > 0.)
        g = ascat(iswband)/scat(iswband);

     tau    (icol, ilev, iswband) = mass(icol, ilev) * ext(iswband);
     tau_w  (icol, ilev, iswband) = mass(icol, ilev) * scat(iswband);
     tau_w_g(icol, ilev, iswband) = mass(icol, ilev) * ascat(iswband);
     tau_w_f(icol, ilev, iswband) = mass(icol, ilev) * g * ascat(iswband);
  });
}

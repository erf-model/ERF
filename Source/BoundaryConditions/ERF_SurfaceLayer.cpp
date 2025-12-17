#include "ERF_SurfaceLayer.H"

using namespace amrex;

/**
 * Wrapper to update ustar and tstar for Monin Obukhov similarity theory.
 *
 * @param[in] lev Current level
 * @param[in] max_iters maximum iterations to use
 */
void
SurfaceLayer::update_fluxes (const int& lev,
                             const Real& time,
                             const MultiFab& cons_in,
                             const std::unique_ptr<MultiFab>& z_phys_nd,
                             int max_iters)
{
    // Update with SST/TSK data if we have a valid pointer
    if (m_sst_lev[lev][0]) {
        fill_tsurf_with_sst_and_tsk(lev, time);
    }

    // Apply heating rate if needed
    if (theta_type == ThetaCalcType::SURFACE_TEMPERATURE) {
        update_surf_temp(time);
    }

    // Update qsurf with qsat over sea
    if (use_moisture) {
        fill_qsurf_with_qsat(lev, cons_in, z_phys_nd);
    }

    // Update land surface temp if we have a valid pointer
    if (m_has_lsm_tsurf) { get_lsm_tsurf(lev); }

    // Fill interior ghost cells
    t_surf[lev]->FillBoundary(m_geom[lev].periodicity());

    // Compute plane averages for all vars (regardless of flux type)
    m_ma.compute_averages(lev);


    // ***************************************************************
    // Iterate the fluxes if moeng type
    // First iterate over land -- the only model for surface roughness
    // over land is RoughCalcType::CONSTANT
    // ***************************************************************
    if (flux_type == FluxCalcType::MOENG ||
        flux_type == FluxCalcType::ROTATE) {
        bool is_land = true;
        // Do we have a constant flux for moisture over land?
        bool cons_qflux = ( (moist_type == MoistCalcType::MOISTURE_FLUX) ||
                            (moist_type == MoistCalcType::ADIABATIC) );
        if (theta_type == ThetaCalcType::HEAT_FLUX) {
            if (rough_type_land == RoughCalcType::CONSTANT) {
                surface_flux most_flux(surf_temp_flux, surf_moist_flux, cons_qflux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else {
                amrex::Abort("Unknown value for rough_type_land");
            }
        } else if (theta_type == ThetaCalcType::SURFACE_TEMPERATURE) {
            if (rough_type_land == RoughCalcType::CONSTANT) {
                surface_temp most_flux(surf_temp_flux, surf_moist_flux, cons_qflux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else {
                amrex::Abort("Unknown value for rough_type_land");
            }
        } else if ((theta_type == ThetaCalcType::ADIABATIC) &&
                   (moist_type == MoistCalcType::ADIABATIC)) {
            if (rough_type_land == RoughCalcType::CONSTANT) {
                adiabatic most_flux(surf_temp_flux, surf_moist_flux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else {
                amrex::Abort("Unknown value for rough_type_land");
            }
        } else {
            amrex::Abort("Unknown value for theta_type");
        }
    } // MOENG -- LAND

    // ***************************************************************
    // Iterate the fluxes if moeng type
    // Next iterate over sea -- the models for surface roughness
    // over sea are CHARNOCK, DONELAN, MODIFIED_CHARNOCK or WAVE_COUPLED
    // ***************************************************************
    if (flux_type == FluxCalcType::MOENG ||
        flux_type == FluxCalcType::ROTATE) {
        bool is_land = false;
        // NOTE: Do not allow default to adiabatic over sea (we have Qvs at surface)
        // Do we have a constant flux for moisture over sea?
        bool cons_qflux = (moist_type == MoistCalcType::MOISTURE_FLUX);
        if (theta_type == ThetaCalcType::HEAT_FLUX) {
            if (rough_type_sea == RoughCalcType::CHARNOCK) {
                surface_flux_charnock most_flux(surf_temp_flux, surf_moist_flux,
                                                cnk_a, cnk_visc, cons_qflux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::MODIFIED_CHARNOCK) {
                surface_flux_mod_charnock most_flux(surf_temp_flux, surf_moist_flux,
                                                    depth, cons_qflux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::DONELAN) {
                surface_flux_donelan most_flux(surf_temp_flux, surf_moist_flux,
                                               cons_qflux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::WAVE_COUPLED) {
                surface_flux_wave_coupled most_flux(surf_temp_flux, surf_moist_flux,
                                                    cons_qflux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else {
                amrex::Abort("Unknown value for rough_type_sea");
            }

        } else if (theta_type == ThetaCalcType::SURFACE_TEMPERATURE) {
            if (rough_type_sea == RoughCalcType::CHARNOCK) {
                surface_temp_charnock most_flux(surf_temp_flux, surf_moist_flux,
                                                cnk_a, cnk_visc, cons_qflux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::MODIFIED_CHARNOCK) {
                surface_temp_mod_charnock most_flux(surf_temp_flux, surf_moist_flux,
                                                    depth, cons_qflux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::DONELAN) {
                surface_temp_donelan most_flux(surf_temp_flux, surf_moist_flux,
                                               cons_qflux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::WAVE_COUPLED) {
                surface_temp_wave_coupled most_flux(surf_temp_flux, surf_moist_flux,
                                                    cons_qflux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else {
                amrex::Abort("Unknown value for rough_type_sea");
            }

        } else if ((theta_type == ThetaCalcType::ADIABATIC) &&
                   (moist_type == MoistCalcType::ADIABATIC)) {
            if (rough_type_sea == RoughCalcType::CHARNOCK) {
                adiabatic_charnock most_flux(surf_temp_flux, surf_moist_flux,
                                             cnk_a, cnk_visc);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::MODIFIED_CHARNOCK) {
                adiabatic_mod_charnock most_flux(surf_temp_flux, surf_moist_flux,
                                                 depth);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::DONELAN) {
                adiabatic_donelan most_flux(surf_temp_flux, surf_moist_flux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::WAVE_COUPLED) {
                adiabatic_wave_coupled most_flux(surf_temp_flux, surf_moist_flux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else {
                amrex::Abort("Unknown value for rough_type_sea");
            }
        } else {
            amrex::Abort("Unknown value for theta_type");
        }

    } // MOENG -- SEA

    if (flux_type == FluxCalcType::CUSTOM) {
        if (custom_rhosurf > 0) {
            specified_rho_surf = true;
            u_star[lev]->setVal(std::sqrt(custom_rhosurf) * custom_ustar);
            t_star[lev]->setVal(custom_rhosurf * custom_tstar);
            q_star[lev]->setVal(custom_rhosurf * custom_qstar);
        } else {
            u_star[lev]->setVal(custom_ustar);
            t_star[lev]->setVal(custom_tstar);
            q_star[lev]->setVal(custom_qstar);
        }
    }
}

/**
 * Function to compute the fluxes (u^star and t^star) for Monin Obukhov similarity theory
 *
 * @param[in] lev Current level
 * @param[in] max_iters maximum iterations to use
 * @param[in] most_flux structure to iteratively compute ustar and tstar
 */
template <typename FluxIter>
void
SurfaceLayer::compute_fluxes (const int& lev,
                              const int& max_iters,
                              const FluxIter& most_flux,
                              bool is_land)
{
    // Pointers to the computed averages
    const auto *const tm_ptr   = m_ma.get_average(lev,2); // potential temperature
    const auto *const qvm_ptr  = m_ma.get_average(lev,3); // water vapor mixing ratio
    const auto *const tvm_ptr  = m_ma.get_average(lev,4); // virtual potential temperature
    const auto *const umm_ptr  = m_ma.get_average(lev,5); // horizontal velocity magnitude
    const auto *const zref_ptr = m_ma.get_zref(lev);     // reference height

    for (MFIter mfi(*u_star[lev]); mfi.isValid(); ++mfi)
    {
        Box gtbx = mfi.growntilebox();

        auto u_star_arr = u_star[lev]->array(mfi);
        auto t_star_arr = t_star[lev]->array(mfi);
        auto q_star_arr = q_star[lev]->array(mfi);
        auto t_surf_arr = t_surf[lev]->array(mfi);
        auto q_surf_arr = q_surf[lev]->array(mfi);
        auto olen_arr   = olen[lev]->array(mfi);

        const auto tm_arr   = tm_ptr->array(mfi);
        const auto tvm_arr  = tvm_ptr->array(mfi);
        const auto qvm_arr  = qvm_ptr->array(mfi);
        const auto umm_arr  = umm_ptr->array(mfi);
        const auto zref_arr = zref_ptr->array(mfi);
        const auto z0_arr   = z_0[lev].array(mfi);

        // PBL height if we need to calculate wstar for the Beljaars correction
        // TODO: can/should we apply this in LES mode?
        const auto w_star_arr = (m_include_wstar) ? w_star[lev].get()->array(mfi) : Array4<Real> {};
        const auto pblh_arr   = (m_include_wstar) ? pblh[lev].get()->array(mfi) : Array4<Real> {};

        // Wave properties if they exist
        const auto Hwave_arr = (m_Hwave_lev[lev]) ? m_Hwave_lev[lev]->array(mfi) : Array4<Real> {};
        const auto Lwave_arr = (m_Lwave_lev[lev]) ? m_Lwave_lev[lev]->array(mfi) : Array4<Real> {};
        const auto eta_arr   = (m_eddyDiffs_lev[lev]) ? m_eddyDiffs_lev[lev]->array(mfi) : Array4<Real> {};

        // Land mask array if it exists
        auto lmask_arr    = (m_lmask_lev[lev][0])    ? m_lmask_lev[lev][0]->array(mfi) :
                                                       Array4<int> {};

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (( is_land && lmask_arr(i,j,k) == 1) ||
                (!is_land && lmask_arr(i,j,k) == 0))
            {
                most_flux.iterate_flux(i, j, k, max_iters,
                                       zref_arr,                            // set in most average
                                       z0_arr,                              // updated if(!is_land)
                                       umm_arr, tm_arr, tvm_arr, qvm_arr,
                                       u_star_arr,                          // updated
                                       w_star_arr,                          // updated if(m_include_wstar)
                                       t_star_arr, q_star_arr,              // updated
                                       t_surf_arr, q_surf_arr, olen_arr,    // updated
                                       pblh_arr,                            // updated if(m_include_wstar)
                                       Hwave_arr, Lwave_arr, eta_arr);
            }
        });
    }
}


/**
 * Wrapper to impose Monin Obukhov similarity theory fluxes by populating ghost cells.
 *
 * @param[in] lev Current level
 * @param[in,out] mfs Multifabs to populate
 * @param[in] eddyDiffs Diffusion coefficients from turbulence model
 */
void
SurfaceLayer::impose_SurfaceLayer_bcs (const int& lev,
                                       Vector<const MultiFab*> mfs,
                                       Vector<std::unique_ptr<MultiFab>>& Tau_lev,
                                       MultiFab* xheat_flux,
                                       MultiFab* yheat_flux,
                                       MultiFab* zheat_flux,
                                       MultiFab* xqv_flux,
                                       MultiFab* yqv_flux,
                                       MultiFab* zqv_flux,
                                       const MultiFab* z_phys)
{
    if (flux_type == FluxCalcType::MOENG) {
        moeng_flux flux_comp;
        compute_SurfaceLayer_bcs(lev, mfs, Tau_lev,
                                 xheat_flux, yheat_flux, zheat_flux,
                                 xqv_flux, yqv_flux, zqv_flux,
                                 z_phys, flux_comp);
    } else if (flux_type == FluxCalcType::DONELAN) {
        donelan_flux flux_comp;
        compute_SurfaceLayer_bcs(lev, mfs, Tau_lev,
                                 xheat_flux, yheat_flux, zheat_flux,
                                 xqv_flux, yqv_flux, zqv_flux,
                                 z_phys, flux_comp);
    } else if (flux_type == FluxCalcType::ROTATE) {
        rotate_flux flux_comp;
        compute_SurfaceLayer_bcs(lev, mfs, Tau_lev,
                                 xheat_flux, yheat_flux, zheat_flux,
                                 xqv_flux, yqv_flux, zqv_flux,
                                 z_phys, flux_comp);
    } else if (flux_type == FluxCalcType::BULK_COEFF) {
        bulk_coeff_flux flux_comp(m_Cd, m_Ch, m_Cq);
        compute_SurfaceLayer_bcs(lev, mfs, Tau_lev,
                                 xheat_flux, yheat_flux, zheat_flux,
                                 xqv_flux, yqv_flux, zqv_flux,
                                 z_phys, flux_comp);
    } else if (flux_type == FluxCalcType::CUSTOM) {
        custom_flux flux_comp(specified_rho_surf);
        compute_SurfaceLayer_bcs(lev, mfs, Tau_lev,
                                 xheat_flux, yheat_flux, zheat_flux,
                                 xqv_flux, yqv_flux, zqv_flux,
                                 z_phys, flux_comp);
    } else {
        amrex::Abort("Unknown surface layer flux calculation type");
    }
}


/**
 * Function to calculate MOST fluxes for populating ghost cells.
 *
 * @param[in] lev Current level
 * @param[in,out] mfs Multifabs to populate
 * @param[in] eddyDiffs Diffusion coefficients from turbulence model
 * @param[in] flux_comp structure to compute fluxes
 */
template <typename FluxCalc>
void
SurfaceLayer::compute_SurfaceLayer_bcs (const int& lev,
                                        Vector<const MultiFab*> mfs,
                                        Vector<std::unique_ptr<MultiFab>>& Tau_lev,
                                        MultiFab* xheat_flux,
                                        MultiFab* yheat_flux,
                                        MultiFab* zheat_flux,
                                        MultiFab* xqv_flux,
                                        MultiFab* yqv_flux,
                                        MultiFab* zqv_flux,
                                        const MultiFab* z_phys,
                                        const FluxCalc& flux_comp)
{
    bool rotate = m_rotate;
    const int klo = m_geom[lev].Domain().smallEnd(2);
    const auto& dxInv = m_geom[lev].InvCellSizeArray();
    for (MFIter mfi(*mfs[0]); mfi.isValid(); ++mfi)
    {
        // Get field arrays
        const auto cons_arr  = mfs[Vars::cons]->array(mfi);
        const auto velx_arr  = mfs[Vars::xvel]->array(mfi);
        const auto vely_arr  = mfs[Vars::yvel]->array(mfi);
        const auto velz_arr  = mfs[Vars::zvel]->array(mfi);

        // Diffusive stress vars
        auto t13_arr =  Tau_lev[TauType::tau13]->array(mfi);
        auto t31_arr = (Tau_lev[TauType::tau31]) ? Tau_lev[TauType::tau31]->array(mfi) : Array4<Real>{};

        auto t23_arr =  Tau_lev[TauType::tau23]->array(mfi);
        auto t32_arr = (Tau_lev[TauType::tau32]) ? Tau_lev[TauType::tau32]->array(mfi) : Array4<Real>{};

        auto hfx3_arr = zheat_flux->array(mfi);
        auto qfx3_arr = (zqv_flux)  ? zqv_flux->array(mfi)   : Array4<Real>{};

        // Rotated stress vars
        auto t11_arr = (m_rotate) ? Tau_lev[TauType::tau11]->array(mfi) : Array4<Real>{};
        auto t22_arr = (m_rotate) ? Tau_lev[TauType::tau22]->array(mfi) : Array4<Real>{};
        auto t33_arr = (m_rotate) ? Tau_lev[TauType::tau33]->array(mfi) : Array4<Real>{};
        auto t12_arr = (m_rotate) ? Tau_lev[TauType::tau12]->array(mfi) : Array4<Real>{};
        auto t21_arr = (m_rotate) ? Tau_lev[TauType::tau21]->array(mfi) : Array4<Real>{};

        auto hfx1_arr = (m_rotate) ? xheat_flux->array(mfi) : Array4<Real>{};
        auto hfx2_arr = (m_rotate) ? yheat_flux->array(mfi) : Array4<Real>{};
        auto qfx1_arr = (m_rotate && xqv_flux) ? xqv_flux->array(mfi) : Array4<Real>{};
        auto qfx2_arr = (m_rotate && yqv_flux) ? yqv_flux->array(mfi) : Array4<Real>{};

        // Terrain
        const auto zphys_arr = (z_phys) ? z_phys->const_array(mfi) : Array4<const Real>{};

        // Get average arrays
        const auto *const u_mean     = m_ma.get_average(lev,0);
        const auto *const v_mean     = m_ma.get_average(lev,1);
        const auto *const t_mean     = m_ma.get_average(lev,2);
        const auto *const q_mean     = m_ma.get_average(lev,3);
        const auto *const u_mag_mean = m_ma.get_average(lev,5);

        const auto um_arr  = u_mean->array(mfi);
        const auto vm_arr  = v_mean->array(mfi);
        const auto tm_arr  = t_mean->array(mfi);
        const auto qm_arr  = q_mean->array(mfi);
        const auto umm_arr = u_mag_mean->array(mfi);

        // Get derived arrays
        const auto u_star_arr = u_star[lev]->array(mfi);
        const auto t_star_arr = t_star[lev]->array(mfi);
        const auto q_star_arr = q_star[lev]->array(mfi);
        const auto t_surf_arr = t_surf[lev]->array(mfi);
        const auto q_surf_arr = q_surf[lev]->array(mfi);

        // Get LSM fluxes
        auto lmask_arr      = (m_lmask_lev[lev][0]) ? m_lmask_lev[lev][0]->array(mfi) :
                                                      Array4<int> {};
        auto lsm_t_flux_arr = Array4<Real> {};
        auto lsm_q_flux_arr = Array4<Real> {};
        auto lsm_tau13_arr  = Array4<Real> {};
        auto lsm_tau23_arr  = Array4<Real> {};
        for (int n(0); n<m_lsm_flux_lev[lev].size(); ++n) {
            if (toLower(m_lsm_flux_name[n]) == "t_flux") { lsm_t_flux_arr = m_lsm_flux_lev[lev][n]->array(mfi); }
            if (toLower(m_lsm_flux_name[n]) == "q_flux") { lsm_q_flux_arr = m_lsm_flux_lev[lev][n]->array(mfi); }
            if (toLower(m_lsm_flux_name[n]) == "tau13")  { lsm_tau13_arr  = m_lsm_flux_lev[lev][n]->array(mfi); }
            if (toLower(m_lsm_flux_name[n]) == "tau23")  { lsm_tau23_arr  = m_lsm_flux_lev[lev][n]->array(mfi); }
        }


        // Rho*Theta flux
        //============================================================================
        Box bx = mfi.tilebox();
        if (bx.smallEnd(2) != klo) { continue; }
        bx.makeSlab(2,klo);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // Valid theta flux from LSM and over land
            Real Tflux;
            int is_land = (lmask_arr) ? lmask_arr(i,j,klo) : 1;
            if (lsm_t_flux_arr && is_land) {
                Tflux = lsm_t_flux_arr(i,j,k);
            } else {
                Tflux = flux_comp.compute_t_flux(i, j, k,
                                                 cons_arr, velx_arr, vely_arr,
                                                 umm_arr, tm_arr, u_star_arr,
                                                 t_star_arr, t_surf_arr);

            }

            // Do scalar flux rotations?
            if (rotate) {
                rotate_scalar_flux(i, j, k, Tflux, dxInv, zphys_arr,
                                   hfx1_arr, hfx2_arr, hfx3_arr);
            } else {
                hfx3_arr(i,j,klo) = Tflux;
            }
        });

        // Rho*Qv flux
        //============================================================================
        if (use_moisture) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Valid qv flux from LSM and over land
                Real Qflux;
                int is_land = (lmask_arr) ? lmask_arr(i,j,klo) : 1;
                if (lsm_q_flux_arr && is_land) {
                    Qflux = lsm_q_flux_arr(i,j,k);
                } else {
                    Qflux = flux_comp.compute_q_flux(i, j, k,
                                                     cons_arr, velx_arr, vely_arr,
                                                     umm_arr, qm_arr, u_star_arr,
                                                     q_star_arr, q_surf_arr);
                }

                // Do scalar flux rotations?
                if (rotate) {
                    rotate_scalar_flux(i, j, k, Qflux, dxInv, zphys_arr,
                                       qfx1_arr, qfx2_arr, qfx3_arr);
                } else {
                    qfx3_arr(i,j,k) = Qflux;
                }
            });
        } // custom

        if (!rotate) {
            // Rho*u flux
            //============================================================================
            Box bxx = surroundingNodes(bx,0);
            ParallelFor(bxx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Valid tau13 from LSM and over land
                Real stressx;
                int is_land_hi = (lmask_arr) ? lmask_arr(i  ,j,klo) : 1;
                int is_land_lo = (lmask_arr) ? lmask_arr(i-1,j,klo) : 1;
                if (lsm_tau13_arr && (is_land_hi || is_land_lo)) {
                    stressx = 0.;
                    if (!is_land_hi || !is_land_lo) {
                        stressx += 0.5 * flux_comp.compute_u_flux(i, j, k,
                                                                  cons_arr, velx_arr, vely_arr,
                                                                  umm_arr, um_arr, u_star_arr);
                    }
                    if (is_land_hi) {
                        stressx += 0.5 * lsm_tau13_arr(i  ,j,k);
                    }
                    if (is_land_lo) {
                        stressx += 0.5 * lsm_tau13_arr(i-1,j,k);
                    }
                } else {
                    stressx = flux_comp.compute_u_flux(i, j, k,
                                                       cons_arr, velx_arr, vely_arr,
                                                       umm_arr, um_arr, u_star_arr);
                }

                t13_arr(i,j,k) = stressx;
                if (t31_arr) { t31_arr(i,j,k) = stressx; }
            });

            // Rho*v flux
            //============================================================================
            Box bxy = surroundingNodes(bx,1);
            ParallelFor(bxy, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Valid tau13 from LSM and over land
                Real stressy;
                int is_land_hi = (lmask_arr) ? lmask_arr(i,j  ,klo) : 1;
                int is_land_lo = (lmask_arr) ? lmask_arr(i,j-1,klo) : 1;
                if (lsm_tau23_arr && (is_land_hi || is_land_lo)) {
                    stressy = 0.;
                    if (!is_land_hi || !is_land_lo) {
                        stressy += 0.5 * flux_comp.compute_v_flux(i, j, k,
                                                                  cons_arr, velx_arr, vely_arr,
                                                                  umm_arr, vm_arr, u_star_arr);
                    }
                    if (is_land_hi) {
                        stressy += 0.5 * lsm_tau23_arr(i,j  ,k);
                    }
                    if (is_land_lo) {
                        stressy += 0.5 * lsm_tau23_arr(i,j-1,k);
                    }
                } else {
                    stressy = flux_comp.compute_v_flux(i, j, k,
                                                       cons_arr, velx_arr, vely_arr,
                                                       umm_arr, vm_arr, u_star_arr);
                }

                t23_arr(i,j,k) = stressy;
                if (t32_arr) { t32_arr(i,j,k) = stressy; }
            });
        } else {
            // All fluxes with rotation
            //============================================================================
            Box bxxy = convert(bx, IntVect(1,1,0));
            ParallelFor(bxxy, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real stresst = flux_comp.compute_u_flux(i, j, k,
                                                       cons_arr, velx_arr, vely_arr,
                                                       umm_arr, um_arr, u_star_arr);
                rotate_stress_tensor(i, j, k, stresst, dxInv, zphys_arr,
                                     velx_arr, vely_arr, velz_arr,
                                     t11_arr, t22_arr, t33_arr,
                                     t12_arr, t21_arr,
                                     t13_arr, t31_arr,
                                     t23_arr, t32_arr);
            });
        }
    } // mfiter
}

void
SurfaceLayer::fill_tsurf_with_sst_and_tsk (const int& lev,
                                           const Real& elapsed_time)
{
    int n_times_in_sst = m_sst_lev[lev].size();

    int n_time_lo, n_time_hi;
    Real alpha;

    if (n_times_in_sst > 1) {
        // Time interpolation
        Real dT = m_bdy_time_interval;
        int n_time = static_cast<int>( elapsed_time /  dT);
        n_time_lo = n_time;
        n_time_hi = n_time+1;
        alpha = (elapsed_time - n_time * dT) / dT;
        if ((elapsed_time == m_stop_time-m_start_time) && (alpha==0)) {
            // stop time coincides with final lowinp slice -- don't try to
            // interpolate from the following time slice
            n_time    -= 1;
            n_time_lo -= 1;
            n_time_hi -= 1;
            alpha = 1.0;
        }
        AMREX_ALWAYS_ASSERT( (n_time >= 0) && (n_time < (m_sst_lev[lev].size()-1)));
    } else {
        n_time_lo = 0;
        n_time_hi = 0;
        alpha     = 1.0;
    }
    AMREX_ALWAYS_ASSERT( alpha >= 0. && alpha <= 1.0);

    Real oma   = 1.0 - alpha;

    // Define a default land surface temperature if we don't read in tsk
    Real lst = default_land_surf_temp;

    bool use_tsk = (m_tsk_lev[lev][0]);
    bool ignore_sst = m_ignore_sst;

    // Populate t_surf
    for (MFIter mfi(*t_surf[lev]); mfi.isValid(); ++mfi)
    {
        Box gtbx = mfi.growntilebox();

        auto t_surf_arr = t_surf[lev]->array(mfi);

        const auto sst_lo_arr = m_sst_lev[lev][n_time_lo]->const_array(mfi);
        const auto sst_hi_arr = m_sst_lev[lev][n_time_hi]->const_array(mfi);

        auto lmask_arr  = (m_lmask_lev[lev][0]) ? m_lmask_lev[lev][0]->array(mfi) :
                                                  Array4<int> {};

        if (use_tsk) {
            const auto tsk_arr = m_tsk_lev[lev][n_time_lo]->const_array(mfi);
            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                int is_land = (lmask_arr) ? lmask_arr(i,j,k) : 1;
                if (!is_land && !ignore_sst) {
                    t_surf_arr(i,j,k) = oma   * sst_lo_arr(i,j,k)
                                      + alpha * sst_hi_arr(i,j,k);
                } else {
                    t_surf_arr(i,j,k) = tsk_arr(i,j,k);
                }
            });
        } else {
            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                int is_land = (lmask_arr) ? lmask_arr(i,j,k) : 1;
                if (!is_land) {
                    t_surf_arr(i,j,k) = oma   * sst_lo_arr(i,j,k)
                                      + alpha * sst_hi_arr(i,j,k);
                } else {
                    t_surf_arr(i,j,k) = lst;
                }
            });
        }
    }
    t_surf[lev]->FillBoundary(m_geom[lev].periodicity());
}

void
SurfaceLayer::fill_qsurf_with_qsat (const int& lev,
                                    const MultiFab& cons_in,
                                    const std::unique_ptr<MultiFab>& z_phys_nd)
{
    // NOTE: We have already tested a moisture model exists

    // Populate q_surf with qsat over water
    auto dz = m_geom[lev].CellSize(2);
    for (MFIter mfi(*q_surf[lev]); mfi.isValid(); ++mfi)
    {
        Box gtbx = mfi.growntilebox();

        auto t_surf_arr = t_surf[lev]->array(mfi);
        auto q_surf_arr = q_surf[lev]->array(mfi);
        auto lmask_arr  = (m_lmask_lev[lev][0]) ? m_lmask_lev[lev][0]->array(mfi) :
                                                  Array4<int> {};
        const auto cons_arr = cons_in.const_array(mfi);
        const auto z_arr    = (z_phys_nd) ? z_phys_nd->const_array(mfi) :
                                            Array4<const Real> {};

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            int is_land = (lmask_arr) ? lmask_arr(i,j,k) : 1;
            if (!is_land) {
                auto deltaZ = (z_arr) ? Compute_Zrel_AtCellCenter(i,j,k,z_arr) :
                                        0.5*dz;
                auto Rho  = cons_arr(i,j,k,Rho_comp);
                auto RTh  = cons_arr(i,j,k,RhoTheta_comp);
                auto Qv   = cons_arr(i,j,k,RhoQ1_comp) / Rho;
                auto P_cc = getPgivenRTh(RTh, Qv);
                P_cc += Rho*CONST_GRAV*deltaZ;
                P_cc *= 0.01;
                erf_qsatw(t_surf_arr(i,j,k), P_cc, q_surf_arr(i,j,k));
            }
        });
    }
    q_surf[lev]->FillBoundary(m_geom[lev].periodicity());
}

void
SurfaceLayer::get_lsm_tsurf (const int& lev)
{
    for (MFIter mfi(*t_surf[lev]); mfi.isValid(); ++mfi)
    {
        Box gtbx = mfi.growntilebox();

        // NOTE: LSM does not carry lateral ghost cells.
        //       This copies the valid box into the ghost cells.
        //       Fillboundary is called after this to pick up the
        //       interior ghost and periodic directions.
        Box vbx  = mfi.validbox();
        int i_lo = vbx.smallEnd(0); int i_hi = vbx.bigEnd(0);
        int j_lo = vbx.smallEnd(1); int j_hi = vbx.bigEnd(1);

        auto t_surf_arr = t_surf[lev]->array(mfi);
        auto lmask_arr  = (m_lmask_lev[lev][0]) ? m_lmask_lev[lev][0]->array(mfi) :
                                                  Array4<int> {};
        const auto lsm_arr = m_lsm_data_lev[lev][m_lsm_tsurf_indx]->const_array(mfi);

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            int is_land = (lmask_arr) ? lmask_arr(i,j,k) : 1;
            if (is_land) {
                int li = amrex::min(amrex::max(i, i_lo), i_hi);
                int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                t_surf_arr(i,j,k) = lsm_arr(li,lj,k);
            }
        });
    }
}

void
SurfaceLayer::update_pblh (const int& lev,
                           Vector<Vector<MultiFab>>& vars,
                           MultiFab* z_phys_cc,
                           const MoistureComponentIndices& moisture_indices)
{
    if (pblh_type == PBLHeightCalcType::MYNN25) {
        MYNNPBLH estimator;
        compute_pblh(lev, vars, z_phys_cc, estimator, moisture_indices);
    } else if (pblh_type == PBLHeightCalcType::YSU || pblh_type == PBLHeightCalcType::MRF) {
        amrex::Error("YSU/MRF PBLH calc not implemented yet");
    }
}

template <typename PBLHeightEstimator>
void
SurfaceLayer::compute_pblh (const int& lev,
                            Vector<Vector<MultiFab>>& vars,
                            MultiFab* z_phys_cc,
                            const PBLHeightEstimator& est,
                            const MoistureComponentIndices& moisture_indices)
{
    est.compute_pblh(m_geom[lev],z_phys_cc, pblh[lev].get(),
                     vars[lev][Vars::cons],m_lmask_lev[lev][0],
                     moisture_indices);
}

void
SurfaceLayer::init_tke_from_ustar (const int& lev,
                                   MultiFab& cons,
                                   const std::unique_ptr<MultiFab>& z_phys_nd,
                                   const Real tkefac,
                                   const Real zscale)
{
    Print() << "Initializing TKE from surface layer ustar on level " << lev << std::endl;

    constexpr Real small = 0.01;

    for (MFIter mfi(cons); mfi.isValid(); ++mfi)
    {
        Box gtbx = mfi.tilebox();

        auto const& u_star_arr = u_star[lev]->const_array(mfi);
        auto const& z_phys_arr = z_phys_nd->const_array(mfi);

        auto const& cons_arr   = cons.array(mfi);

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real rho = cons_arr(i, j, k, Rho_comp);
            Real ust = u_star_arr(i, j, 0);
            Real tke0 = tkefac * ust * ust; // surface value
            Real zagl = Compute_Zrel_AtCellCenter(i, j, k, z_phys_arr);

            // linearly tapering profile --  following WRF, approximate top of
            // PBL as ustar * zscale
            cons_arr(i, j, k, RhoKE_comp) = rho * tke0 * std::max(
                (ust * zscale - zagl) / (std::max(ust, small) * zscale),
                small);
        });
    }
}


void
SurfaceLayer::read_custom_roughness (const int& lev,
                                     const std::string& fname)
{
    // Read the file if we have it
    if (!fname.empty()) {
        // Only the ioproc reads the file
        Gpu::HostVector<Real> m_x,m_y,m_z0;
        if (ParallelDescriptor::IOProcessor()) {
            Print()<<"Reading MOST roughness file at level " << lev << " : " << fname << std::endl;
            std::ifstream file(fname);
            Real value1,value2,value3;
            while(file>>value1>>value2>>value3){
                m_x.push_back(value1);
                m_y.push_back(value2);
                m_z0.push_back(value3);
            }
            file.close();

            AMREX_ALWAYS_ASSERT(m_x.size() == m_y.size());
            AMREX_ALWAYS_ASSERT(m_x.size() == m_z0.size());
        }

        // Broadcast the whole domain to every rank
        int ioproc = ParallelDescriptor::IOProcessorNumber();
        int nnode = m_x.size();
        ParallelDescriptor::Bcast(&nnode, 1, ioproc);

        if (!ParallelDescriptor::IOProcessor()) {
            m_x.resize(nnode);
            m_y.resize(nnode);
            m_z0.resize(nnode);
        }
        ParallelDescriptor::Bcast(m_x.data() , nnode, ioproc);
        ParallelDescriptor::Bcast(m_y.data() , nnode, ioproc);
        ParallelDescriptor::Bcast(m_z0.data(), nnode, ioproc);

        // Copy data to the GPU
        Gpu::DeviceVector<Real> d_x(nnode),d_y(nnode),d_z0(nnode);
        Gpu::copy(Gpu::hostToDevice, m_x.begin(), m_x.end(), d_x.begin());
        Gpu::copy(Gpu::hostToDevice, m_y.begin(), m_y.end(), d_y.begin());
        Gpu::copy(Gpu::hostToDevice, m_z0.begin(), m_z0.end(), d_z0.begin());
        Real* xp  = d_x.data();
        Real* yp  = d_y.data();
        Real* z0p = d_z0.data();

        // Each rank populates it's z_0[lev] MultiFab
        for (MFIter mfi(z_0[lev]); mfi.isValid(); ++mfi)
        {
            Box gtbx = mfi.growntilebox();

            // Populate z_phys data
            Real tol = 1.0e-4;
            auto dx = m_geom[lev].CellSizeArray();
            auto ProbLoArr = m_geom[lev].ProbLoArray();
            int ilo = m_geom[lev].Domain().smallEnd(0);
            int jlo = m_geom[lev].Domain().smallEnd(1);
            int klo = 0;
            int ihi = m_geom[lev].Domain().bigEnd(0);
            int jhi = m_geom[lev].Domain().bigEnd(1);

            Array4<Real> const& z0_arr = z_0[lev].array(mfi);
            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/)
            {
                // Clip indices for ghost-cells
                int ii = amrex::min(amrex::max(i,ilo),ihi);
                int jj = amrex::min(amrex::max(j,jlo),jhi);

                // Location of nodes
                Real x = ProbLoArr[0]  + ii  * dx[0];
                Real y = ProbLoArr[1]  + jj  * dx[1];
                int inode = ii + jj * (ihi-ilo+2); // stride is Nx+1
                if (std::sqrt(std::pow(x-xp[inode],2)+std::pow(y-yp[inode],2)) < tol) {
                    z0_arr(i,j,klo) = z0p[inode];
                } else {
                    // Unexpected list order, do brute force search
                    Real z0loc = 0.0;
                    bool found = false;
                    for (int n=0; n<nnode; ++n) {
                        Real delta=std::sqrt(std::pow(x-xp[n],2)+std::pow(y-yp[n],2));
                        if (delta < tol) {
                            found = true;
                            z0loc = z0p[n];
                            break;
                        }
                    }
                    AMREX_ASSERT_WITH_MESSAGE(found, "Location read from terrain file does not match the grid!");
                    amrex::ignore_unused(found);
                    z0_arr(i,j,klo) = z0loc;
                }
            });
        } // mfi
    } else {
        AMREX_ALWAYS_ASSERT(lev > 0);

        Print()<<"Interpolating MOST roughness at level " << lev << std::endl;

        // Create a BC mapper that uses FOEXTRAP at domain bndry
        Vector<int> bc_lo(3,ERFBCType::foextrap);
        Vector<int> bc_hi(3,ERFBCType::foextrap);
        Vector<BCRec> bcr; bcr.push_back(BCRec(bc_lo.data(),bc_hi.data()));

        // Create ref ratio
        IntVect ratio;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            ratio[idim] = m_geom[lev].Domain().length(idim) / m_geom[0].Domain().length(idim);
        }

        // Create interp object and interpolate from the coarsest grid
        MFInterpolater* interp = &mf_cell_cons_interp;
        interp->interp(z_0[0]  , 0,
                       z_0[lev], 0,
                       1, z_0[lev].nGrowVect(),
                       m_geom[0], m_geom[lev],
                       m_geom[lev].Domain(),ratio,
                       bcr, 0);
    }
}

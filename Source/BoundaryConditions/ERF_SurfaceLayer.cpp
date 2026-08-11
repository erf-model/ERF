#include "ERF_SurfaceLayer.H"
#include "ERF_SurfaceLayerStress.H"

using namespace amrex;

/**
 * Wrapper to update ustar and tstar for Monin Obukhov similarity theory.
 *
 * @param[in] lev Current level
 * @param[in] elapsed_time Current simulation time
 * @param[in] elapsed_time_since_start_low Time since the start of the lower-boundary data
 * @param[in,out] cons_in Conserved state, updated when RANS TKE is initialized from surface-layer data
 * @param[in] z_phys_nd Nodal physical height used by terrain-aware surface calculations
 * @param[in] walldist Wall distance used when updating boundary TKE
 * @param[in] max_iters Maximum iterations to use in the MOST flux solve
 */
void
SurfaceLayer::update_fluxes (const int& lev,
                             const double& elapsed_time,
                             const double& elapsed_time_since_start_low,
                             MultiFab& cons_in,
                             const std::unique_ptr<MultiFab>& z_phys_nd,
                             const std::unique_ptr<MultiFab>& walldist,
                             int max_iters)
{
    // Update with SST/TSK data if we have a valid pointer
    if (!m_has_ocean_lsm_tsurf &&
        !m_sst_lev[lev].empty() && m_sst_lev[lev][0]) {
        fill_tsurf_with_sst_and_tsk(lev, elapsed_time_since_start_low);
    }

    // Apply heating rate if needed
    if (theta_type == ThetaCalcType::SURFACE_TEMPERATURE &&
        !m_has_ocean_lsm_tsurf) {
        update_surf_temp(elapsed_time_since_start_low);
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

    // NOTE: Do iterations to seed variables on the first step (LSM called post step)
    //       as well as compute values where invalid LSM fluxes may reside
    //*******************************************************************************
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
        if (m_terrain_type != TerrainType::EB) {
            if (theta_type == ThetaCalcType::HEAT_FLUX) {
                if (rough_type_land == RoughCalcType::CONSTANT) {
                    surface_flux most_flux(surf_temp_flux, surf_moist_flux, cons_qflux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else {
                    amrex::Abort("Unknown value for rough_type_land");
                }
            } else if (theta_type == ThetaCalcType::SURFACE_TEMPERATURE) {
                if (rough_type_land == RoughCalcType::CONSTANT) {
                    surface_temp most_flux(surf_temp_flux, surf_moist_flux, cons_qflux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else {
                    amrex::Abort("Unknown value for rough_type_land");
                }
            } else if ((theta_type == ThetaCalcType::ADIABATIC) &&
                       (moist_type == MoistCalcType::ADIABATIC)) {
                if (rough_type_land == RoughCalcType::CONSTANT) {
                    adiabatic most_flux(surf_temp_flux, surf_moist_flux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else {
                    amrex::Abort("Unknown value for rough_type_land");
                }
            } else {
                amrex::Abort("Unknown value for theta_type");
            }
        // EB
        } else {
            if (theta_type == ThetaCalcType::HEAT_FLUX) {
                if (rough_type_land == RoughCalcType::CONSTANT) {
                    surface_flux_eb most_flux(surf_temp_flux, surf_moist_flux, cons_qflux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else {
                    amrex::Abort("Unknown value for rough_type_land");
                }
            } else if (theta_type == ThetaCalcType::SURFACE_TEMPERATURE) {
                if (rough_type_land == RoughCalcType::CONSTANT) {
                    surface_temp_eb most_flux(surf_temp_flux, surf_moist_flux, cons_qflux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else {
                    amrex::Abort("Unknown value for rough_type_land");
                }
            } else if ((theta_type == ThetaCalcType::ADIABATIC) &&
                       (moist_type == MoistCalcType::ADIABATIC)) {
                if (rough_type_land == RoughCalcType::CONSTANT) {
                    adiabatic_eb most_flux(surf_temp_flux, surf_moist_flux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else {
                    amrex::Abort("Unknown value for rough_type_land");
                }
            } else {
                amrex::Abort("Unknown value for theta_type");
            }
        } // EB
    } // MOENG -- LAND

    // Update u*/T*/q*/L over land (iterations or from LSM fluxes)
    if (m_has_lsm_fluxes && elapsed_time > zero) {
        compute_sfc_params_from_lsm_fluxes(lev, cons_in);
    }

    // ***************************************************************
    // Iterate the fluxes if moeng type
    // Next iterate over sea -- the models for surface roughness
    // over sea are CHARNOCK, DONELAN, MODIFIED_CHARNOCK or WAVE_COUPLED
    // NOTE: Sea surface fluxes are not supported for EB terrain
    // ***************************************************************
    if ((flux_type == FluxCalcType::MOENG ||
         flux_type == FluxCalcType::ROTATE) &&
        m_terrain_type != TerrainType::EB) {
        bool is_land = false;
        // NOTE: Do not allow default to adiabatic over sea (we have Qvs at surface)
        // Do we have a constant flux for moisture over sea?
        bool cons_qflux = (moist_type == MoistCalcType::MOISTURE_FLUX);
            if (theta_type == ThetaCalcType::HEAT_FLUX) {
                if (rough_type_sea == RoughCalcType::CHARNOCK) {
                    surface_flux_charnock most_flux(surf_temp_flux, surf_moist_flux,
                                                    cnk_a, cnk_visc, cons_qflux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else if (rough_type_sea == RoughCalcType::MODIFIED_CHARNOCK) {
                    surface_flux_mod_charnock most_flux(surf_temp_flux, surf_moist_flux,
                                                        depth, cons_qflux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else if (rough_type_sea == RoughCalcType::DONELAN) {
                    surface_flux_donelan most_flux(surf_temp_flux, surf_moist_flux,
                                                cons_qflux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else if (rough_type_sea == RoughCalcType::WAVE_COUPLED) {
                    surface_flux_wave_coupled most_flux(surf_temp_flux, surf_moist_flux,
                                                        cons_qflux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else {
                    amrex::Abort("Unknown value for rough_type_sea");
                }

            } else if (theta_type == ThetaCalcType::SURFACE_TEMPERATURE) {
                if (rough_type_sea == RoughCalcType::CHARNOCK) {
                    surface_temp_charnock most_flux(surf_temp_flux, surf_moist_flux,
                                                    cnk_a, cnk_visc, cons_qflux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else if (rough_type_sea == RoughCalcType::MODIFIED_CHARNOCK) {
                    surface_temp_mod_charnock most_flux(surf_temp_flux, surf_moist_flux,
                                                        depth, cons_qflux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else if (rough_type_sea == RoughCalcType::DONELAN) {
                    surface_temp_donelan most_flux(surf_temp_flux, surf_moist_flux,
                                                cons_qflux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else if (rough_type_sea == RoughCalcType::WAVE_COUPLED) {
                    surface_temp_wave_coupled most_flux(surf_temp_flux, surf_moist_flux,
                                                        cons_qflux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else {
                    amrex::Abort("Unknown value for rough_type_sea");
                }

            } else if ((theta_type == ThetaCalcType::ADIABATIC) &&
                    (moist_type == MoistCalcType::ADIABATIC)) {
                if (rough_type_sea == RoughCalcType::CHARNOCK) {
                    adiabatic_charnock most_flux(surf_temp_flux, surf_moist_flux,
                                                cnk_a, cnk_visc);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else if (rough_type_sea == RoughCalcType::MODIFIED_CHARNOCK) {
                    adiabatic_mod_charnock most_flux(surf_temp_flux, surf_moist_flux,
                                                    depth);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else if (rough_type_sea == RoughCalcType::DONELAN) {
                    adiabatic_donelan most_flux(surf_temp_flux, surf_moist_flux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else if (rough_type_sea == RoughCalcType::WAVE_COUPLED) {
                    adiabatic_wave_coupled most_flux(surf_temp_flux, surf_moist_flux);
                    compute_fluxes(lev, max_iters, cons_in, most_flux, is_land);
                } else {
                    amrex::Abort("Unknown value for rough_type_sea");
                }
            } else {
                amrex::Abort("Unknown value for theta_type");
            }
    } // MOENG -- SEA

    if (flux_type == FluxCalcType::CUSTOM || flux_type == FluxCalcType::RICO) {
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

    if (m_update_k_rans) {
        const bool use_ref_theta = (theta_ref > 0);
        const Real l_inv_theta0  = (use_ref_theta) ? one / theta_ref : one;
        const Real l_inv_Cmu2    = inv_Cmu2;
        const int klo = m_geom[lev].Domain().smallEnd(2);
        IntVect ng = u_star[lev]->nGrowVect(); ng[2] = 0;

        for (MFIter mfi(cons_in); mfi.isValid(); ++mfi)
        {
            Box gpbx = mfi.tilebox(IntVect(0),ng);

            if (gpbx.smallEnd(2) != klo) { continue; }

            gpbx.makeSlab(2,klo);

            auto cons_arr = cons_in.array(mfi);
            const auto& u_star_arr = u_star[lev]->const_array(mfi);
            const auto& t_star_arr = t_star[lev]->const_array(mfi);
            const auto& dist_arr   = walldist->const_array(mfi);

            ParallelFor(gpbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real rho = cons_arr(i,j,k,Rho_comp);
                if (t_star_arr(i,j,0) < -1e-8) {
                    // Only destabilizing buoyancy flux affects the boundary k
                    // tstar < 0 ==> B > 0
                    Real B = -CONST_GRAV * l_inv_theta0 * u_star_arr(i,j,0) * t_star_arr(i,j,0);
                    if (!use_ref_theta) {
                        B *= cons_arr(i,j,k,Rho_comp) /
                             cons_arr(i,j,k,RhoTheta_comp);
                    }

                    // Axell & Liungman 2001, Eqn. 16
                    cons_arr(i,j,k,RhoKE_comp) = rho * l_inv_Cmu2 *
                        std::pow(
                            u_star_arr(i,j,0) * u_star_arr(i,j,0) * u_star_arr(i,j,0)
                            + KAPPA * B * dist_arr(i,j,k),
                        two/three);
                } else {
                    cons_arr(i,j,k,RhoKE_comp) = rho * l_inv_Cmu2 * u_star_arr(i,j,0) * u_star_arr(i,j,0);
                }
            });
        }
    }

    u_star[lev]->FillBoundary(m_geom[lev].periodicity());
    t_star[lev]->FillBoundary(m_geom[lev].periodicity());
    q_star[lev]->FillBoundary(m_geom[lev].periodicity());
      olen[lev]->FillBoundary(m_geom[lev].periodicity());
}

/**
 * Function to compute the fluxes (u^star and t^star) for Monin Obukhov similarity theory
 *
 * @param[in] lev Current level
 * @param[in] max_iters Maximum iterations to use
 * @param[in] cons_in Conserved state whose grids define the surface iteration
 * @param[in] most_flux Flux-iteration functor used to compute ustar, tstar, qstar, and related fields
 * @param[in] is_land Selects whether land or sea cells are updated
 */
template <typename FluxIter>
void
SurfaceLayer::compute_fluxes (const int& lev,
                              const int& max_iters,
                              MultiFab& cons_in,
                              const FluxIter& most_flux,
                              bool is_land)
{
    // Pointers to the computed averages
    const auto *const tm_ptr   = m_ma.get_average(lev,2); // potential temperature
    const auto *const qvm_ptr  = m_ma.get_average(lev,3); // water vapor mixing ratio
    const auto *const tvm_ptr  = m_ma.get_average(lev,4); // virtual potential temperature
    const auto *const umm_ptr  = m_ma.get_average(lev,5); // horizontal velocity magnitude
    const auto *const zref_ptr = m_ma.get_zref(lev);     // reference height
    const bool l_use_eb = (m_terrain_type == TerrainType::EB);

    const int klo = m_geom[lev].Domain().smallEnd(2);
    IntVect ng = u_star[lev]->nGrowVect(); ng[2] = 0;

    for (MFIter mfi(cons_in); mfi.isValid(); ++mfi)
    {
        Box gtbx = mfi.tilebox(IntVect(0),ng);

        if (!l_use_eb && gtbx.smallEnd(2) != klo) { continue; }

        if (!l_use_eb) { gtbx.makeSlab(2,klo); }

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

        // Get EB flags if needed
        const auto flag_arr = (l_use_eb) ? m_eb_vec[lev]->get_const_factory()->getMultiEBCellFlagFab()[mfi].const_array() : Array4<const EBCellFlag>{};

        if (!l_use_eb) {
            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
            {
                if (( is_land && lmask_arr(i,j,0) == 1) ||
                    (!is_land && lmask_arr(i,j,0) == 0))
                {
                    // NOTE: All 2D MFs so k index is always 0 from ba2d definition
                    most_flux.iterate_flux(i, j, 0, max_iters,
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
        // EB
        } else {
            if (std::is_same<FluxIter, adiabatic_eb>::value ||
                std::is_same<FluxIter, surface_temp_eb>::value ||
                std::is_same<FluxIter, surface_flux_eb>::value) {
                ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (( is_land && lmask_arr(i,j,0) == 1) ||
                        (!is_land && lmask_arr(i,j,0) == 0))
                    {
                        if (flag_arr(i,j,k).isSingleValued()) {
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
                    }
                });
            } else {
                amrex::Abort("FluxIter type not supported for EB");
            }
        }
    }
}


/**
 * Wrapper to impose Monin Obukhov similarity theory fluxes by populating ghost cells.
 *
 * @param[in] lev Current level
 * @param[in] mfs State MultiFabs used to compute the boundary fluxes
 * @param[in,out] Tau_lev Diffusive stress MultiFabs populated with surface stresses
 * @param[in,out] xheat_flux x-face heat-flux MultiFab, used when rotated fluxes are enabled
 * @param[in,out] yheat_flux y-face heat-flux MultiFab, used when rotated fluxes are enabled
 * @param[in,out] zheat_flux z-face heat-flux MultiFab populated with vertical surface heat flux
 * @param[in,out] xqv_flux x-face moisture-flux MultiFab, used when rotated fluxes and moisture are enabled
 * @param[in,out] yqv_flux y-face moisture-flux MultiFab, used when rotated fluxes and moisture are enabled
 * @param[in,out] zqv_flux z-face moisture-flux MultiFab populated when moisture is enabled
 * @param[in] z_phys Nodal physical height used to rotate terrain-following fluxes
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
    } else if (flux_type == FluxCalcType::RICO) {
        rico_flux flux_comp(rico_theta_z0, rico_qsat_z0);
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
 * Wrapper to impose Monin Obukhov similarity theory fluxes by populating ghost cells.
 *
 * @param[in] lev Current level
 * @param[in] mfs State MultiFabs used to compute the EB boundary fluxes
 * @param[in,out] Tau_EB EB diffusive stress MultiFabs populated with surface stresses
 * @param[in,out] xheat_flux x-face EB heat-flux MultiFab, currently unused
 * @param[in,out] yheat_flux y-face EB heat-flux MultiFab, currently unused
 * @param[in,out] Hfx3_EB EB heat-flux MultiFab populated with scalar surface flux
 * @param[in,out] xqv_flux x-face EB moisture-flux MultiFab, currently unused
 * @param[in,out] yqv_flux y-face EB moisture-flux MultiFab, currently unused
 * @param[in,out] zqv_flux z-face EB moisture-flux MultiFab, currently unused
 */
void
SurfaceLayer::impose_SurfaceLayer_bcs_EB (const int& lev,
                                       Vector<const MultiFab*> mfs,
                                       Vector<Vector<std::unique_ptr<MultiFab>>>& Tau_EB,
                                       MultiFab* xheat_flux,
                                       MultiFab* yheat_flux,
                                       MultiFab* Hfx3_EB,
                                       MultiFab* xqv_flux,
                                       MultiFab* yqv_flux,
                                       MultiFab* zqv_flux)
{
    if (flux_type == FluxCalcType::MOENG) {
        moeng_flux_eb flux_comp;
        compute_SurfaceLayer_bcs_EB(lev, mfs, Tau_EB,
                                 xheat_flux, yheat_flux, Hfx3_EB,
                                 xqv_flux, yqv_flux, zqv_flux,
                                 flux_comp);
    } else {
        amrex::Abort("Not implemented surface layer flux calculation type for EB");
    }
}

/**
 * Function to calculate MOST fluxes for populating ghost cells.
 *
 * @param[in] lev Current level
 * @param[in] mfs State MultiFabs used to compute the boundary fluxes
 * @param[in,out] Tau_lev Diffusive stress MultiFabs populated with surface stresses
 * @param[in,out] xheat_flux x-face heat-flux MultiFab, used when rotated fluxes are enabled
 * @param[in,out] yheat_flux y-face heat-flux MultiFab, used when rotated fluxes are enabled
 * @param[in,out] zheat_flux z-face heat-flux MultiFab populated with vertical surface heat flux
 * @param[in,out] xqv_flux x-face moisture-flux MultiFab, used when rotated fluxes and moisture are enabled
 * @param[in,out] yqv_flux y-face moisture-flux MultiFab, used when rotated fluxes and moisture are enabled
 * @param[in,out] zqv_flux z-face moisture-flux MultiFab populated when moisture is enabled
 * @param[in] z_phys Nodal physical height used to rotate terrain-following fluxes
 * @param[in] flux_comp Flux-calculation functor used to compute scalar and momentum fluxes
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

        auto olen_arr   = olen[lev]->array(mfi);

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
        auto surface_source_arr = surface_diagnostic_source[lev]->array(mfi);

        // Get LSM fluxes
        auto lmask_arr      = (m_lmask_lev[lev][0]) ? m_lmask_lev[lev][0]->array(mfi) :
                                                      Array4<int> {};
        auto lsm_t_flux_arr  = Array4<Real> {};
        auto soil_t_flux_arr = Array4<Real> {};
        auto lsm_q_flux_arr = Array4<Real> {};
        auto lsm_tau13_arr  = Array4<Real> {};
        auto lsm_tau23_arr  = Array4<Real> {};
        // LSM tau fields are cell-centered kinematic stresses [m2 s-2].
        // Tau_lev tau13/tau23 are face-centered conservative stresses [N m-2].
        for (int n(0); n<m_lsm_flux_lev[lev].size(); ++n) {
            if (toLower(m_lsm_flux_name[n]) == "t_flux")      { lsm_t_flux_arr  = m_lsm_flux_lev[lev][n]->array(mfi); }
            if (toLower(m_lsm_flux_name[n]) == "soil_t_flux") { soil_t_flux_arr = m_lsm_flux_lev[lev][n]->array(mfi); }
            if (toLower(m_lsm_flux_name[n]) == "q_flux") { lsm_q_flux_arr = m_lsm_flux_lev[lev][n]->array(mfi); }
            if (toLower(m_lsm_flux_name[n]) == "tau13")  { lsm_tau13_arr  = m_lsm_flux_lev[lev][n]->array(mfi); }
            if (toLower(m_lsm_flux_name[n]) == "tau23")  { lsm_tau23_arr  = m_lsm_flux_lev[lev][n]->array(mfi); }
        }

        const bool has_lsm_t_flux = static_cast<bool>(lsm_t_flux_arr);
        const bool is_custom = (flux_type == FluxCalcType::CUSTOM);
        const bool is_rico   = (flux_type == FluxCalcType::RICO);


        // Rho*Theta flux
        //============================================================================
        Box bx = mfi.tilebox();
        if (bx.smallEnd(2) != klo) { continue; }
        bx.makeSlab(2,klo);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // Valid theta flux from LSM and over land. The LSM writes the
            // lsm_undefined sentinel for cells it did not process (sea-ice /
            // open water); fall back to MOST there instead of applying garbage.
            Real Tflux;
            int is_land = (lmask_arr) ? lmask_arr(i,j,0) : 1;
            const bool lsm_flux_is_valid = (lsm_t_flux_arr) ? (lsm_t_flux_arr(i,j,0) < lsm_undefined) :
                                                              false;
            const bool has_land_and_flux = (is_land == 1 && lsm_flux_is_valid);
            if (lsm_t_flux_arr && has_land_and_flux) {
                // LSM flux MultiFabs store kinematic fluxes for MOST parameter
                // updates. The applied hfx array stores the conservative RHS flux.
                Tflux = cons_arr(i,j,k,Rho_comp) * lsm_t_flux_arr(i,j,0);
            } else if (is_land == 2) { // no temperature flux within buildings
                Tflux = zero;
            } else {
                Tflux = flux_comp.compute_t_flux(i, j, k,
                                                 cons_arr, velx_arr, vely_arr,
                                                 umm_arr, tm_arr, u_star_arr,
                                                 t_star_arr, t_surf_arr);
                // NOTE: do NOT write the MOST-fallback flux back into lsm_t_flux_arr.
                // Doing so flips a sentinel (water/unprocessed) cell to "valid LSM"
                // on the next step, so a MOST-derived value is re-read as an LSM flux
                // Only Noah-MP should populate the LSM cache.
            }

            if (soil_t_flux_arr && is_land == 1) {
                soil_t_flux_arr(i,j,k) = Tflux / cons_arr(i,j,k,Rho_comp);
            }

            surface_source_arr(i,j,0) = surface_diagnostics::to_plot_value(
                surface_diagnostics::classify_scalar_source(
                    is_custom, is_rico, is_land, has_lsm_t_flux, lsm_flux_is_valid));

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
                // Valid qv flux from LSM and over land (sentinel -> fall back to MOST)
                Real Qflux;
                int is_land = (lmask_arr) ? lmask_arr(i,j,0) : 1;
                const bool lsm_flux_is_valid = (lsm_q_flux_arr) ? (lsm_q_flux_arr(i,j,0) < lsm_undefined) :
                                                                  false;
                const bool has_land_and_flux = (is_land == 1 && lsm_flux_is_valid);
                if (lsm_q_flux_arr && has_land_and_flux) {
                    // LSM flux MultiFabs store kinematic fluxes for MOST parameter
                    // updates. The applied qfx array stores the conservative RHS flux.
                    Qflux = cons_arr(i,j,k,Rho_comp) * lsm_q_flux_arr(i,j,0);
                } else if (is_land == 2) { // no moisture flux within buildings
                    Qflux = zero;
                } else {
                    Qflux = flux_comp.compute_q_flux(i, j, k,
                                                     cons_arr, velx_arr, vely_arr,
                                                     umm_arr, qm_arr, u_star_arr,
                                                     q_star_arr, q_surf_arr);
                    // NOTE: no writeback into lsm_q_flux_arr -- see the matching
                    // t_flux note above.
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
                // Valid tau13 from LSM and over land. A side that is land but
                // whose LSM flux is the sentinel (sea-ice / open water) is treated
                // as non-LSM so that side uses the MOST stress instead.
                Real stressx;
                int is_land_hi = (lmask_arr) ? lmask_arr(i  ,j,0) : 1;
                int is_land_lo = (lmask_arr) ? lmask_arr(i-1,j,0) : 1;
                const bool lsm_hi_flux_is_valid = surface_layer_stress::lsm_flux_is_valid(
                    static_cast<bool>(lsm_tau13_arr), is_land_hi == 1,
                    lsm_tau13_arr ? lsm_tau13_arr(i  ,j,0) : zero, lsm_undefined);
                const bool lsm_lo_flux_is_valid = surface_layer_stress::lsm_flux_is_valid(
                    static_cast<bool>(lsm_tau13_arr), is_land_lo == 1,
                    lsm_tau13_arr ? lsm_tau13_arr(i-1,j,0) : zero, lsm_undefined);
                const bool has_land_and_flux_hi = (is_land_hi == 1 && lsm_hi_flux_is_valid);
                const bool has_land_and_flux_lo = (is_land_lo == 1 && lsm_lo_flux_is_valid);
                if (lsm_tau13_arr && (has_land_and_flux_hi || has_land_and_flux_lo)) {
                    const Real rho_hi = cons_arr(i  ,j,k,Rho_comp);
                    const Real rho_lo = cons_arr(i-1,j,k,Rho_comp);
                    const Real most_stress = (!has_land_and_flux_hi || !has_land_and_flux_lo) ?
                        flux_comp.compute_u_flux(i, j, k,
                                                 cons_arr, velx_arr, vely_arr,
                                                 umm_arr, um_arr, u_star_arr) : zero;
                    const auto result = surface_layer_stress::combine_lsm_and_most_stress(
                        rho_lo, rho_hi, lsm_tau13_arr(i-1,j,0), lsm_tau13_arr(i,j,0),
                        has_land_and_flux_lo, has_land_and_flux_hi, most_stress);
                    stressx = result.face_stress;
                    // NOTE: do NOT write the MOST-fallback stress back into the
                    // cell-centered lsm_tau13_arr. This face-indexed ParallelFor
                    // touches cells (i) and (i-1), so each cell is written by two
                    // adjacent face threads in the same launch -> nondeterministic
                    // write-write race on GPU (ERF #3446). It also spuriously flips
                    // a sentinel (water/unprocessed) cell to "valid LSM" for the
                    // next step. The face stress is fully determined here; the LSM
                    // cache is (re)filled only by Noah-MP. Matches baseline 3ab899d3.
                } else if (is_land_hi == 2 || is_land_lo == 2) { // no stress within buildings
                    stressx = zero;
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
                // Valid tau23 from LSM and over land (sentinel side -> MOST stress)
                Real stressy;
                int is_land_hi = (lmask_arr) ? lmask_arr(i,j  ,0) : 1;
                int is_land_lo = (lmask_arr) ? lmask_arr(i,j-1,0) : 1;
                const bool lsm_hi_flux_is_valid = surface_layer_stress::lsm_flux_is_valid(
                    static_cast<bool>(lsm_tau23_arr), is_land_hi == 1,
                    lsm_tau23_arr ? lsm_tau23_arr(i,j  ,0) : zero, lsm_undefined);
                const bool lsm_lo_flux_is_valid = surface_layer_stress::lsm_flux_is_valid(
                    static_cast<bool>(lsm_tau23_arr), is_land_lo == 1,
                    lsm_tau23_arr ? lsm_tau23_arr(i,j-1,0) : zero, lsm_undefined);
                const bool has_land_and_flux_hi = (is_land_hi == 1 && lsm_hi_flux_is_valid);
                const bool has_land_and_flux_lo = (is_land_lo == 1 && lsm_lo_flux_is_valid);
                if (lsm_tau23_arr && (has_land_and_flux_hi || has_land_and_flux_lo)) {
                    const Real rho_hi = cons_arr(i,j  ,k,Rho_comp);
                    const Real rho_lo = cons_arr(i,j-1,k,Rho_comp);
                    const Real most_stress = (!has_land_and_flux_hi || !has_land_and_flux_lo) ?
                        flux_comp.compute_v_flux(i, j, k,
                                                 cons_arr, velx_arr, vely_arr,
                                                 umm_arr, vm_arr, u_star_arr) : zero;
                    const auto result = surface_layer_stress::combine_lsm_and_most_stress(
                        rho_lo, rho_hi, lsm_tau23_arr(i,j-1,0), lsm_tau23_arr(i,j,0),
                        has_land_and_flux_lo, has_land_and_flux_hi, most_stress);
                    stressy = result.face_stress;
                    // NOTE: no writeback into cell-centered lsm_tau23_arr -- see the
                    // matching tau13 note above (ERF #3446 write-write race + stale
                    // sentinel-becomes-valid). Face stress is complete here.
                } else if (is_land_hi == 2 || is_land_lo == 2) { // no stress within buildings
                    stressy = zero;
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

        // For models that do not do iterations to yield u*/T*/q*,
        // fill these values from the fluxes that were computed.

        // NOTE: For LSM, this has been handled in "compute_sfc_params_from_lsm_fluxes"
        // NOTE: Fluxes here are for conserved quantities, we divide by rho
        if (flux_type == FluxCalcType::BULK_COEFF ||
            flux_type == FluxCalcType::DONELAN) {
            constexpr Real eps = std::numeric_limits<Real>::epsilon();
            bool l_use_moisture = use_moisture;
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/)
            {
                Real rho = cons_arr(i,j,klo,Rho_comp);
                Real Thd = cons_arr(i,j,klo,RhoTheta_comp) / rho;
                Real qv  = (l_use_moisture) ? cons_arr(i,j,klo,RhoQ1_comp) / rho : zero;
                Real Thv = Thd * (one + epsv*qv);

                Real tau = std::sqrt( t13_arr(i,j,klo)/rho * t13_arr(i,j,klo)/rho
                                    + t23_arr(i,j,klo)/rho * t23_arr(i,j,klo)/rho );
                u_star_arr(i,j,0) = amrex::max(std::sqrt(tau),eps);

                if (hfx3_arr(i,j,klo)>=zero) {
                    t_star_arr(i,j,0) = amrex::min(-hfx3_arr(i,j,klo) / (rho * u_star_arr(i,j,0)),-eps);
                } else {
                    t_star_arr(i,j,0) = amrex::max(-hfx3_arr(i,j,klo) / (rho * u_star_arr(i,j,0)),eps);
                }
                if (!l_use_moisture) {
                    q_star_arr(i,j,0) = zero;
                } else if (qfx3_arr(i,j,klo)>=zero) {
                    q_star_arr(i,j,0) = amrex::min(-qfx3_arr(i,j,klo) / (rho * u_star_arr(i,j,0)),-eps);
                } else {
                    q_star_arr(i,j,0) = amrex::max(-qfx3_arr(i,j,klo) / ( rho * u_star_arr(i,j,0)),eps);
                }
                olen_arr(i,j,0)   = ( u_star_arr(i,j,0) * u_star_arr(i,j,0) * Thv ) /
                                    ( KAPPA * CONST_GRAV * t_star_arr(i,j,0) );
            });
        }

    } // mfiter

    surface_diagnostic_source[lev]->FillBoundary(m_geom[lev].periodicity());
}

/**
 * Function to calculate MOST fluxes for EB.
 *
 * @param[in] lev Current level
 * @param[in] mfs State MultiFabs used to compute the EB boundary fluxes
 * @param[in,out] Tau_EB EB diffusive stress MultiFabs populated with surface stresses
 * @param[in,out] xheat_flux x-face EB heat-flux MultiFab, currently unused
 * @param[in,out] yheat_flux y-face EB heat-flux MultiFab, currently unused
 * @param[in,out] Hfx3_EB EB heat-flux MultiFab populated with scalar surface flux
 * @param[in,out] xqv_flux x-face EB moisture-flux MultiFab, currently unused
 * @param[in,out] yqv_flux y-face EB moisture-flux MultiFab, currently unused
 * @param[in,out] zqv_flux z-face EB moisture-flux MultiFab, currently unused
 * @param[in] flux_comp EB flux-calculation functor used to compute scalar and momentum fluxes
 */
template <typename FluxCalc>
void
SurfaceLayer::compute_SurfaceLayer_bcs_EB (const int& lev,
                                        Vector<const MultiFab*> mfs,
                                        Vector<Vector<std::unique_ptr<MultiFab>>>& Tau_EB,
                                        [[maybe_unused]] MultiFab* xheat_flux,
                                        [[maybe_unused]] MultiFab* yheat_flux,
                                        MultiFab* Hfx3_EB,
                                        [[maybe_unused]] MultiFab* xqv_flux,
                                        [[maybe_unused]] MultiFab* yqv_flux,
                                        [[maybe_unused]] MultiFab* zqv_flux,
                                        const FluxCalc& flux_comp)
{
    // Get EB flags for all centerings
    const auto& cc_factory = m_eb_vec[lev]->get_const_factory();
    const auto& cc_flags = cc_factory->getMultiEBCellFlagFab();
    const auto& cc_vfrac = cc_factory->getVolFrac();

    const auto& u_factory = m_eb_vec[lev]->get_u_const_factory();
    const auto& u_flags = u_factory->getMultiEBCellFlagFab();
    const auto& u_vfrac = u_factory->getVolFrac();

    const auto& v_factory = m_eb_vec[lev]->get_v_const_factory();
    const auto& v_flags = v_factory->getMultiEBCellFlagFab();
    const auto& v_vfrac = v_factory->getVolFrac();

    const auto& w_factory = m_eb_vec[lev]->get_w_const_factory();
    const auto& w_flags = w_factory->getMultiEBCellFlagFab();
    const auto& w_vfrac = w_factory->getVolFrac();

    // EB does not currently have a cell-centered scalar-source classification.
    // Keep the provenance mask missing rather than inventing face-aware
    // semantics for the staggered stress path.
    surface_diagnostic_source[lev]->setVal(
        surface_diagnostics::to_plot_value(surface_diagnostics::SurfaceDiagnosticSource::Missing));

    for (MFIter mfi(*mfs[0]); mfi.isValid(); ++mfi)
    {
        // Get flags for this box (all centerings)
        const auto& cc_flag = cc_flags[mfi];
        const auto& u_flag = u_flags[mfi];
        const auto& v_flag = v_flags[mfi];
        const auto& w_flag = w_flags[mfi];

        // Skip boxes that have no cut cells at any centering
        if (cc_flag.getType() != FabType::singlevalued &&
            u_flag.getType() != FabType::singlevalued &&
            v_flag.getType() != FabType::singlevalued &&
            w_flag.getType() != FabType::singlevalued
        ) continue;

        // Get EB flag and volfrac arrays
        auto const cc_flag_arr = cc_flag.const_array();
        auto const u_flag_arr = u_flag.const_array();
        auto const v_flag_arr = v_flag.const_array();
        auto const w_flag_arr = w_flag.const_array();

        auto const cc_vfrac_arr = cc_vfrac.const_array(mfi);
        auto const u_vfrac_arr = u_vfrac.const_array(mfi);
        auto const v_vfrac_arr = v_vfrac.const_array(mfi);
        auto const w_vfrac_arr = w_vfrac.const_array(mfi);

        // Get boundary normals only if cut cells exist
        auto const bnorm_arr = (cc_flag.getType() == FabType::singlevalued) ?
            cc_factory->getBndryNormal().const_array(mfi) : Array4<const Real>{};
        auto const u_bnorm_arr = (u_flag.getType() == FabType::singlevalued) ?
            u_factory->getBndryNormal().const_array(mfi) : Array4<const Real>{};
        auto const v_bnorm_arr = (v_flag.getType() == FabType::singlevalued) ?
            v_factory->getBndryNormal().const_array(mfi) : Array4<const Real>{};
        auto const w_bnorm_arr = (w_flag.getType() == FabType::singlevalued) ?
            w_factory->getBndryNormal().const_array(mfi) : Array4<const Real>{};

        // Get field arrays
        const auto cons_arr  = mfs[Vars::cons]->array(mfi);
        const auto velx_arr  = mfs[Vars::xvel]->array(mfi);
        const auto vely_arr  = mfs[Vars::yvel]->array(mfi);
        const auto velz_arr  = mfs[Vars::zvel]->array(mfi);

        // Diffusive stress vars - t13 and t23 components for all grid types
        auto u_t13_arr = Tau_EB[EBTauType::tau_eb13][EBGridType::xface]->array(mfi);
        auto v_t13_arr = Tau_EB[EBTauType::tau_eb13][EBGridType::yface]->array(mfi);
        auto w_t13_arr = Tau_EB[EBTauType::tau_eb13][EBGridType::zface]->array(mfi);

        auto u_t23_arr = Tau_EB[EBTauType::tau_eb23][EBGridType::xface]->array(mfi);
        auto v_t23_arr = Tau_EB[EBTauType::tau_eb23][EBGridType::yface]->array(mfi);
        auto w_t23_arr = Tau_EB[EBTauType::tau_eb23][EBGridType::zface]->array(mfi);

        auto hfx3_arr = Hfx3_EB->array(mfi);

        // Get average arrays
        const auto *const u_mean     = m_ma.get_average(lev,0);
        const auto *const v_mean     = m_ma.get_average(lev,1);
        const auto *const t_mean     = m_ma.get_average(lev,2);
        // const auto *const q_mean     = m_ma.get_average(lev,3);
        const auto *const u_mag_mean = m_ma.get_average(lev,5);

        const auto um_arr  = u_mean->array(mfi);
        const auto vm_arr  = v_mean->array(mfi);
        const auto tm_arr  = t_mean->array(mfi);
        // const auto qm_arr  = q_mean->array(mfi);
        const auto umm_arr = u_mag_mean->array(mfi);

        // Get derived arrays
        const auto u_star_arr = u_star[lev]->array(mfi);
        const auto t_star_arr = t_star[lev]->array(mfi);
        const auto t_surf_arr = t_surf[lev]->array(mfi);

        // Rho*Theta flux
        //============================================================================
        Box bx = mfi.tilebox();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (cc_flag_arr(i,j,k).isSingleValued()) {
                Real Tflux = flux_comp.compute_t_flux(i, j, k,
                                                    cons_arr, velx_arr, vely_arr, velz_arr,
                                                    umm_arr, tm_arr, u_star_arr,
                                                    t_star_arr, t_surf_arr,
                                                    u_vfrac_arr, v_vfrac_arr, w_vfrac_arr,
                                                    bnorm_arr);
                hfx3_arr(i,j,k) = Tflux;
            }
        });

        // Rho*u flux
        //============================================================================
        Box bxx = surroundingNodes(bx,0);
        Box bxy = surroundingNodes(bx,1);
        Box bxz = surroundingNodes(bx,2);
        ParallelFor(bxx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (u_flag_arr(i,j,k).isSingleValued()) {
                Real stressx = flux_comp.compute_u_flux(i, j, k,
                                                        cons_arr, velx_arr, vely_arr, velz_arr,
                                                        umm_arr, um_arr, u_star_arr,
                                                        u_vfrac_arr, v_vfrac_arr, w_vfrac_arr,
                                                        cc_vfrac_arr, cc_flag_arr,
                                                        u_bnorm_arr, 0);
                u_t13_arr(i,j,k) = stressx;
            }
        });
        ParallelFor(bxy, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (v_flag_arr(i,j,k).isSingleValued()) {
                Real stressx = flux_comp.compute_u_flux(i, j, k,
                                                        cons_arr, velx_arr, vely_arr, velz_arr,
                                                        umm_arr, um_arr, u_star_arr,
                                                        u_vfrac_arr, v_vfrac_arr, w_vfrac_arr,
                                                        cc_vfrac_arr, cc_flag_arr,
                                                        v_bnorm_arr, 1);
                v_t13_arr(i,j,k) = stressx;
            }
        });
        ParallelFor(bxz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (w_flag_arr(i,j,k).isSingleValued()) {
                Real stressx = flux_comp.compute_u_flux(i, j, k,
                                                        cons_arr, velx_arr, vely_arr, velz_arr,
                                                        umm_arr, um_arr, u_star_arr,
                                                        u_vfrac_arr, v_vfrac_arr, w_vfrac_arr,
                                                        cc_vfrac_arr, cc_flag_arr,
                                                        w_bnorm_arr, 2);
                w_t13_arr(i,j,k) = stressx;
            }
        });

        // Rho*v flux
        //============================================================================
        ParallelFor(bxx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (u_flag_arr(i,j,k).isSingleValued()) {
                Real stressy = flux_comp.compute_v_flux(i, j, k,
                                                        cons_arr, velx_arr, vely_arr, velz_arr,
                                                        umm_arr, vm_arr, u_star_arr,
                                                        u_vfrac_arr, v_vfrac_arr, w_vfrac_arr,
                                                        cc_vfrac_arr, cc_flag_arr, u_bnorm_arr, 0);
                u_t23_arr(i,j,k) = stressy;
            }
        });
        ParallelFor(bxy, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (v_flag_arr(i,j,k).isSingleValued()) {
                Real stressy = flux_comp.compute_v_flux(i, j, k,
                                                        cons_arr, velx_arr, vely_arr, velz_arr,
                                                        umm_arr, vm_arr, u_star_arr,
                                                        u_vfrac_arr, v_vfrac_arr, w_vfrac_arr,
                                                        cc_vfrac_arr, cc_flag_arr, v_bnorm_arr, 1);
                v_t23_arr(i,j,k) = stressy;
            }
        });
        ParallelFor(bxz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (w_flag_arr(i,j,k).isSingleValued()) {
                Real stressy = flux_comp.compute_v_flux(i, j, k,
                                                        cons_arr, velx_arr, vely_arr, velz_arr,
                                                        umm_arr, vm_arr, u_star_arr,
                                                        u_vfrac_arr, v_vfrac_arr, w_vfrac_arr,
                                                        cc_vfrac_arr, cc_flag_arr, w_bnorm_arr, 2);
                w_t23_arr(i,j,k) = stressy;
            }
        });
    } // mfiter

}

/**
 * Compute surface-layer parameters from land-surface-model fluxes.
 *
 * @param[in] lev Current level
 * @param[in] cons_in Conserved state used to derive density, theta, and moisture at the surface
 */
void
SurfaceLayer::compute_sfc_params_from_lsm_fluxes (const int& lev,
                                                  MultiFab& cons_in)
{
    Real eps = std::numeric_limits<amrex::Real>::epsilon();
    bool has_moisture = use_moisture;
    const int klo = m_geom[lev].Domain().smallEnd(2);
    for (MFIter mfi(cons_in); mfi.isValid(); ++mfi) {

        Box vbx = mfi.validbox();
        if (vbx.smallEnd(2) != klo) { continue; }
        vbx.makeSlab(2,0);

        // Get CC state
        const Array4<const Real> cons_arr = cons_in.const_array(mfi);

        // Get SL params
        const auto u_star_arr = u_star[lev]->array(mfi);
        const auto t_star_arr = t_star[lev]->array(mfi);
        const auto q_star_arr = q_star[lev]->array(mfi);
        const auto olen_arr   = olen[lev]->array(mfi);

        // Get LSM fluxes
        auto lmask_arr      = (m_lmask_lev[lev][0]) ? m_lmask_lev[lev][0]->array(mfi) :
                                                      Array4<int> {};
        auto lsm_t_flux_arr = Array4<Real> {};
        auto lsm_q_flux_arr = Array4<Real> {};
        auto lsm_tau13_arr  = Array4<Real> {};
        auto lsm_tau23_arr  = Array4<Real> {};
        // compute_sfc_params_from_lsm_fluxes consumes signed kinematic stress
        // components; their vector magnitude determines u_star^2.
        for (int n(0); n<m_lsm_flux_lev[lev].size(); ++n) {
            if (toLower(m_lsm_flux_name[n]) == "t_flux") { lsm_t_flux_arr = m_lsm_flux_lev[lev][n]->array(mfi); }
            if (toLower(m_lsm_flux_name[n]) == "q_flux") { lsm_q_flux_arr = m_lsm_flux_lev[lev][n]->array(mfi); }
            if (toLower(m_lsm_flux_name[n]) == "tau13")  { lsm_tau13_arr  = m_lsm_flux_lev[lev][n]->array(mfi); }
            if (toLower(m_lsm_flux_name[n]) == "tau23")  { lsm_tau23_arr  = m_lsm_flux_lev[lev][n]->array(mfi); }
        }

        ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
        {
            int is_land = (lmask_arr) ? lmask_arr(i,j,0) : 1;
            // Skip cells the LSM did not have a valid flux (lsm_undefined).
            if (is_land && lsm_t_flux_arr && lsm_t_flux_arr(i,j,0) < lsm_undefined) {
                Real rho = cons_arr(i,j,klo,Rho_comp);
                Real Thd = cons_arr(i,j,klo,RhoTheta_comp) / rho;
                Real qv  = (has_moisture) ? cons_arr(i,j,klo,RhoQ1_comp) / rho : zero;
                Real Thv = Thd * (one + epsv*qv);
                Real tau = std::sqrt( lsm_tau13_arr(i,j,0)*lsm_tau13_arr(i,j,0)
                                    + lsm_tau23_arr(i,j,0)*lsm_tau23_arr(i,j,0) );
                u_star_arr(i,j,0) = amrex::max(std::sqrt(tau),eps);
                if (lsm_t_flux_arr(i,j,0)>=zero) {
                    t_star_arr(i,j,0) = amrex::min(-lsm_t_flux_arr(i,j,0) / u_star_arr(i,j,0),-eps);
                } else {
                    t_star_arr(i,j,0) = amrex::max(-lsm_t_flux_arr(i,j,0) / u_star_arr(i,j,0),eps);
                }
                if (lsm_q_flux_arr(i,j,0)>=zero) {
                    q_star_arr(i,j,0) = amrex::min(-lsm_q_flux_arr(i,j,0) / u_star_arr(i,j,0),-eps);
                } else {
                    q_star_arr(i,j,0) = amrex::max(-lsm_q_flux_arr(i,j,0) / u_star_arr(i,j,0),eps);
                }
                olen_arr(i,j,0)   = ( u_star_arr(i,j,0) * u_star_arr(i,j,0) * Thv ) /
                                    ( KAPPA * CONST_GRAV * t_star_arr(i,j,0) );
            }
        });
    } // mfi
}

/**
 * Fill surface temperature from SST/TSK lower-boundary data.
 *
 * @param[in] lev Current level
 * @param[in] elapsed_time_since_start_low Time since the start of the lower-boundary data
 */
void
SurfaceLayer::fill_tsurf_with_sst_and_tsk (const int& lev,
                                           const double& elapsed_time_since_start_low)
{
    int n_times_in_sst = static_cast<int>(m_sst_lev[lev].size());

    double dT = m_low_time_interval;

    int n_time_lo, n_time_hi;
    Real alpha;

    if (n_times_in_sst > 1) {
        n_time_lo = static_cast<int>( elapsed_time_since_start_low /  dT);
        alpha = static_cast<Real>((elapsed_time_since_start_low - n_time_lo * dT) / dT);

        AMREX_ALWAYS_ASSERT( alpha >= zero && alpha <= one);

        n_time_hi = n_time_lo + 1;

        // Do not over run the last sst file
        if (m_start_low_time + elapsed_time_since_start_low >= m_final_low_time) {
            n_time_lo = static_cast<int>(m_sst_lev[lev].size())-1;
            n_time_hi = n_time_lo;
            alpha     = zero;
        }

        AMREX_ALWAYS_ASSERT( (n_time_lo >= 0) && (n_time_hi < m_sst_lev[lev].size()));
    } else {
        n_time_lo = 0;
        n_time_hi = 0;
        alpha     = one;
    }
    AMREX_ALWAYS_ASSERT( alpha >= zero && alpha <= one);

    Real oma   = one - alpha;

    // Define a default land surface temperature if we don't read in tsk
    Real lst = default_land_surf_temp;

    bool use_tsk = (m_tsk_lev[lev][0]);
    bool ignore_sst = m_ignore_sst;

    const int klo = m_geom[lev].Domain().smallEnd(2);

    // Populate t_surf
    for (MFIter mfi(*t_surf[lev]); mfi.isValid(); ++mfi)
    {
        Box gtbx = mfi.growntilebox();

        if (gtbx.smallEnd(2) != klo) { continue; }

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

/**
 * Fill sea-surface moisture with saturation specific humidity.
 *
 * @param[in] lev Current level
 * @param[in] cons_in Conserved state used to derive pressure at the surface
 * @param[in] z_phys_nd Nodal physical height used to compute terrain-relative surface height
 */
void
SurfaceLayer::fill_qsurf_with_qsat (const int& lev,
                                    const MultiFab& cons_in,
                                    const std::unique_ptr<MultiFab>& z_phys_nd)
{
    // NOTE: We have already tested a moisture model exists

    // Populate q_surf with qsat over water
    auto dz = m_geom[lev].CellSize(2);
    const int klo = m_geom[lev].Domain().smallEnd(2);
    for (MFIter mfi(*q_surf[lev]); mfi.isValid(); ++mfi)
    {
        Box gtbx = mfi.growntilebox();

        if (gtbx.smallEnd(2) != klo) { continue; }

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
                                        myhalf*dz;
                auto Rho  = cons_arr(i,j,k,Rho_comp);
                auto RTh  = cons_arr(i,j,k,RhoTheta_comp);
                auto Qv   = cons_arr(i,j,k,RhoQ1_comp) / Rho;
                auto P_cc = getPgivenRTh(RTh, Qv);
                P_cc += Rho*CONST_GRAV*deltaZ;
                P_cc *= Real(0.01);
                erf_qsatw(t_surf_arr(i,j,k), P_cc, q_surf_arr(i,j,k));
            }
        });
    }
    q_surf[lev]->FillBoundary(m_geom[lev].periodicity());
}

/**
 * Fill surface temperature from land-surface-model data.
 *
 * @param[in] lev Current level
 */
void
SurfaceLayer::get_lsm_tsurf (const int& lev)
{
    const int klo = m_geom[lev].Domain().smallEnd(2);
    const bool has_sea_tsurf = (m_has_ocean_lsm_tsurf &&
                                amrex::toLower(m_lsm_data_name[m_lsm_tsurf_indx]) == "t_surf");
    for (MFIter mfi(*t_surf[lev]); mfi.isValid(); ++mfi)
    {
        Box gtbx = mfi.growntilebox();

        if (gtbx.smallEnd(2) != klo) { continue; }

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
            if ((!has_sea_tsurf && is_land) ||
                (has_sea_tsurf && !is_land)) {
                int li = amrex::min(amrex::max(i, i_lo), i_hi);
                int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                t_surf_arr(i,j,k) = lsm_arr(li,lj,k);
            }
        });
    }
}

/**
 * Update PBL height using the configured estimator.
 *
 * @param[in] lev Current level
 * @param[in] vars Level-indexed state MultiFabs passed to the PBL height estimator
 * @param[in] z_phys_cc Cell-centered physical height used by the PBL height estimator
 * @param[in] moisture_indices Moisture component indices used by the PBL height estimator
 */
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

/**
 * Compute PBL height with the supplied estimator.
 *
 * @param[in] lev Current level
 * @param[in] vars Level-indexed state MultiFabs passed to the estimator
 * @param[in] z_phys_cc Cell-centered physical height used by the estimator
 * @param[in] est PBL height estimator functor
 * @param[in] moisture_indices Moisture component indices used by the estimator
 */
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

/**
 * Initialize TKE from surface-layer friction velocity.
 *
 * @param[in] lev Current level
 * @param[in,out] cons Conserved state whose RhoKE component is initialized
 * @param[in] z_phys_nd Nodal physical height used to compute height above ground
 * @param[in] tkefac Factor multiplying ustar squared for the surface TKE value
 * @param[in] zscale Scale factor used to taper TKE with height
 */
void
SurfaceLayer::init_tke_from_ustar (const int& lev,
                                   MultiFab& cons,
                                   const std::unique_ptr<MultiFab>& z_phys_nd,
                                   const Real tkefac,
                                   const Real zscale)
{
    Print() << "Initializing TKE from surface layer ustar on level " << lev << std::endl;

    // Handle vertical decomposition by selectively copying into
    // a FArrayBox section on each rank. Then doing a reduce real sum
    // and broadcasting to each rank. No mask since all CC data
    const int klo = m_geom[lev].Domain().smallEnd(2);
    Box bx_lo = u_star[lev]->boxArray().minimalBox();
    FArrayBox u_star_lo(bx_lo, 1); u_star_lo.setVal<RunOn::Device>(0);
    FArrayBox z_surf_lo(bx_lo, 1); z_surf_lo.setVal<RunOn::Device>(0);
    Real* ustar_ptr = u_star_lo.dataPtr();
    Real* zsurf_ptr = z_surf_lo.dataPtr();
    for (MFIter mfi(cons); mfi.isValid(); ++mfi)
    {
        Box vbx = mfi.validbox();
        if (vbx.smallEnd(2) != klo) { continue; }
        vbx.makeSlab(2,0);

        auto const& u_star_arr = u_star[lev]->const_array(mfi);
        auto        u_star_all = u_star_lo.array();

        auto const& z_phys_arr = z_phys_nd->const_array(mfi);
        auto        z_surf_all = z_surf_lo.array();

        ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
        {
            u_star_all(i,j,0) = u_star_arr(i,j,0);
            z_surf_all(i,j,0) = fourth * ( z_phys_arr(i  ,j  ,klo) + z_phys_arr(i+1,j  ,klo)
                                       + z_phys_arr(i  ,j+1,klo) + z_phys_arr(i+1,j+1,klo) );
        });
    }
    ParallelDescriptor::ReduceRealSum(ustar_ptr, static_cast<int>(bx_lo.numPts()));
    ParallelDescriptor::ReduceRealSum(zsurf_ptr, static_cast<int>(bx_lo.numPts()));

    // Now work on all boxes (ustar has been filled above)
    constexpr Real small = Real(0.01);
    for (MFIter mfi(cons); mfi.isValid(); ++mfi)
    {
        Box vbx  = mfi.validbox();

        auto const& u_star_arr = u_star_lo.const_array();
        auto const& z_surf_arr = z_surf_lo.const_array();
        auto const& z_phys_arr = z_phys_nd->const_array(mfi);

        auto const& cons_arr   = cons.array(mfi);

        ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real rho  = cons_arr(i, j, k, Rho_comp);
            Real ust  = u_star_arr(i, j, 0);
            Real tke0 = tkefac * ust * ust; // surface value
            Real zagl = Compute_Z_AtCellCenter(i, j, k, z_phys_arr) - z_surf_arr(i,j,0);

            // linearly tapering profile --  following WRF, approximate top of
            // PBL as ustar * zscale
            cons_arr(i, j, k, RhoKE_comp) = rho * tke0 * std::max(
                (ust * zscale - zagl) / (std::max(ust, small) * zscale),
                small);
        });
    }
}


/**
 * Read or interpolate custom roughness length data.
 *
 * @param[in] lev Current level
 * @param[in] fname Roughness file name; an empty name interpolates from level 0
 */
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
        int nnode = static_cast<int>(m_x.size());
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
        const int klo = m_geom[lev].Domain().smallEnd(2);
        for (MFIter mfi(z_0[lev]); mfi.isValid(); ++mfi)
        {
            Box gtbx = mfi.growntilebox();

            if (gtbx.smallEnd(2) != klo) { continue; }

            // Populate z_phys data
            Real tol = Real(1.0e-4);
            auto dx = m_geom[lev].CellSizeArray();
            auto ProbLoArr = m_geom[lev].ProbLoArray();
            int ilo = m_geom[lev].Domain().smallEnd(0);
            int jlo = m_geom[lev].Domain().smallEnd(1);
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
                if (std::sqrt(amrex::Math::powi<2>(x-xp[inode])+amrex::Math::powi<2>(y-yp[inode])) < tol) {
                    z0_arr(i,j,klo) = z0p[inode];
                } else {
                    // Unexpected list order, do brute force search
                    Real z0loc = zero;
                    bool found = false;
                    for (int n=0; n<nnode; ++n) {
                        Real delta=std::sqrt(amrex::Math::powi<2>(x-xp[n])+amrex::Math::powi<2>(y-yp[n]));
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

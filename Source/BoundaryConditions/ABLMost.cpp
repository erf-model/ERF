#include <ABLMost.H>

using namespace amrex;

/**
 * Wrapper to update ustar and tstar for Monin Obukhov similarity theory.
 *
 * @param[in] lev Current level
 * @param[in] max_iters maximum iterations to use
 */
void
ABLMost::update_fluxes (const int& lev,
                        const Real& time,
                        int max_iters)
{
    // Update SST data if we have a valid pointer
    if (m_sst_lev[lev][0]) time_interp_sst(lev, time);

    // TODO: we want 0 index to always be theta?
    // Update land surface temp if we have a valid pointer
    if (m_lsm_data_lev[lev][0]) get_lsm_tsurf(lev);

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
        if (theta_type == ThetaCalcType::HEAT_FLUX) {
            if (rough_type_land == RoughCalcType::CONSTANT) {
                surface_flux most_flux(m_ma.get_zref(), surf_temp_flux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else {
                amrex::Abort("Unknown value for rough_type_land");
            }
        } else if (theta_type == ThetaCalcType::SURFACE_TEMPERATURE) {
            update_surf_temp(time);
            if (rough_type_land == RoughCalcType::CONSTANT) {
                surface_temp most_flux(m_ma.get_zref(), surf_temp_flux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else {
                amrex::Abort("Unknown value for rough_type_land");
            }
        } else if (theta_type == ThetaCalcType::ADIABATIC) {
            if (rough_type_land == RoughCalcType::CONSTANT) {
                adiabatic most_flux(m_ma.get_zref(), surf_temp_flux);
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
    // over sea are CHARNOCK, MODIFIED_CHARNOCK or WAVE_COUPLED
    // ***************************************************************
    if (flux_type == FluxCalcType::MOENG ||
        flux_type == FluxCalcType::ROTATE) {
        bool is_land = false;
        if (theta_type == ThetaCalcType::HEAT_FLUX) {
            if (rough_type_sea == RoughCalcType::CHARNOCK) {
                surface_flux_charnock most_flux(m_ma.get_zref(), surf_temp_flux, cnk_a);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::MODIFIED_CHARNOCK) {
                surface_flux_mod_charnock most_flux(m_ma.get_zref(), surf_temp_flux, depth);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::WAVE_COUPLED) {
                surface_flux_wave_coupled most_flux(m_ma.get_zref(), surf_temp_flux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else {
                amrex::Abort("Unknown value for rough_type_sea");
            }

        } else if (theta_type == ThetaCalcType::SURFACE_TEMPERATURE) {
            update_surf_temp(time);
            if (rough_type_sea == RoughCalcType::CHARNOCK) {
                surface_temp_charnock most_flux(m_ma.get_zref(), surf_temp_flux, cnk_a);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::MODIFIED_CHARNOCK) {
                surface_temp_mod_charnock most_flux(m_ma.get_zref(), surf_temp_flux, depth);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::WAVE_COUPLED) {
                surface_temp_wave_coupled most_flux(m_ma.get_zref(), surf_temp_flux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else {
                amrex::Abort("Unknown value for rough_type_sea");
            }

        } else if (theta_type == ThetaCalcType::ADIABATIC) {
            if (rough_type_sea == RoughCalcType::CHARNOCK) {
                adiabatic_charnock most_flux(m_ma.get_zref(), surf_temp_flux, cnk_a);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::MODIFIED_CHARNOCK) {
                adiabatic_mod_charnock most_flux(m_ma.get_zref(), surf_temp_flux, depth);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else if (rough_type_sea == RoughCalcType::WAVE_COUPLED) {
                adiabatic_wave_coupled most_flux(m_ma.get_zref(), surf_temp_flux);
                compute_fluxes(lev, max_iters, most_flux, is_land);
            } else {
                amrex::Abort("Unknown value for rough_type_sea");
            }
        } else {
            amrex::Abort("Unknown value for theta_type");
        }

    } // MOENG -- SEA

    if (flux_type == FluxCalcType::CUSTOM) {
        u_star[lev]->setVal(custom_ustar);
        t_star[lev]->setVal(custom_tstar);
        q_star[lev]->setVal(custom_qstar);
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
ABLMost::compute_fluxes (const int& lev,
                         const int& max_iters,
                         const FluxIter& most_flux,
                         bool is_land)
{
    // Pointers to the computed averages
    const auto *const tm_ptr  = m_ma.get_average(lev,2); // potential temperature
    const auto *const qvm_ptr = m_ma.get_average(lev,3); // water vapor mixing ratio
    const auto *const tvm_ptr = m_ma.get_average(lev,4); // virtual potential temperature
    const auto *const umm_ptr = m_ma.get_average(lev,5); // horizontal velocity magnitude

    for (MFIter mfi(*u_star[lev]); mfi.isValid(); ++mfi)
    {
        Box gtbx = mfi.growntilebox();

        auto u_star_arr = u_star[lev]->array(mfi);
        auto t_star_arr = t_star[lev]->array(mfi);
        auto q_star_arr = q_star[lev]->array(mfi);
        auto t_surf_arr = t_surf[lev]->array(mfi);
        auto olen_arr   = olen[lev]->array(mfi);
        auto pblh_arr   = pblh[lev]->array(mfi);

        const auto tm_arr  = tm_ptr->array(mfi);
        const auto tvm_arr = tvm_ptr->array(mfi);
        const auto qvm_arr = qvm_ptr->array(mfi);
        const auto umm_arr = umm_ptr->array(mfi);
        const auto z0_arr  = z_0[lev].array();

        // Wave properties if they exist
        const auto Hwave_arr = (m_Hwave_lev[lev]) ? m_Hwave_lev[lev]->array(mfi) : Array4<Real> {};
        const auto Lwave_arr = (m_Lwave_lev[lev]) ? m_Lwave_lev[lev]->array(mfi) : Array4<Real> {};
        const auto eta_arr   = m_eddyDiffs_lev[lev]->array(mfi);

        auto lmask_arr    = (m_lmask_lev[lev][0])    ? m_lmask_lev[lev][0]->array(mfi) :
                                                       Array4<int> {};

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (( is_land && lmask_arr(i,j,k) == 1) ||
                (!is_land && lmask_arr(i,j,k) == 0))
            {
                most_flux.iterate_flux(i, j, k, max_iters,
                                       z0_arr, umm_arr, tm_arr, tvm_arr, qvm_arr,
                                       u_star_arr, t_star_arr, q_star_arr,  // to be updated
                                       t_surf_arr, olen_arr,                // to be updated
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
ABLMost::impose_most_bcs (const int& lev,
                          const Vector<MultiFab*>& mfs,
                          MultiFab* xxmom_flux,
                          MultiFab* yymom_flux,
                          MultiFab* zzmom_flux,
                          MultiFab* xymom_flux, MultiFab* yxmom_flux,
                          MultiFab* xzmom_flux, MultiFab* zxmom_flux,
                          MultiFab* yzmom_flux, MultiFab* zymom_flux,
                          MultiFab* xheat_flux,
                          MultiFab* yheat_flux,
                          MultiFab* zheat_flux,
                          MultiFab* xqv_flux,
                          MultiFab* yqv_flux,
                          MultiFab* zqv_flux,
                          MultiFab* z_phys)
{
    const int klo = 0;
    if (flux_type == FluxCalcType::MOENG) {
        moeng_flux flux_comp(klo);
        compute_most_bcs(lev, mfs,
                         xxmom_flux,
                         yymom_flux,
                         zzmom_flux,
                         xymom_flux, yxmom_flux,
                         xzmom_flux, zxmom_flux,
                         yzmom_flux, zymom_flux,
                         xheat_flux, yheat_flux, zheat_flux,
                         xqv_flux, yqv_flux, zqv_flux,
                         z_phys, flux_comp);
    } else if (flux_type == FluxCalcType::DONELAN) {
        donelan_flux flux_comp(klo);
        compute_most_bcs(lev, mfs,
                         xxmom_flux,
                         yymom_flux,
                         zzmom_flux,
                         xymom_flux, yxmom_flux,
                         xzmom_flux, zxmom_flux,
                         yzmom_flux, zymom_flux,
                         xheat_flux, yheat_flux, zheat_flux,
                         xqv_flux, yqv_flux, zqv_flux,
                         z_phys, flux_comp);
    } else if (flux_type == FluxCalcType::ROTATE) {
        rotate_flux flux_comp(klo);
        compute_most_bcs(lev, mfs,
                         xxmom_flux,
                         yymom_flux,
                         zzmom_flux,
                         xymom_flux, yxmom_flux,
                         xzmom_flux, zxmom_flux,
                         yzmom_flux, zymom_flux,
                         xheat_flux, yheat_flux, zheat_flux,
                         xqv_flux, yqv_flux, zqv_flux,
                         z_phys, flux_comp);
    } else {
        custom_flux flux_comp(klo);
        compute_most_bcs(lev, mfs,
                         xxmom_flux,
                         yymom_flux,
                         zzmom_flux,
                         xymom_flux, yxmom_flux,
                         xzmom_flux, zxmom_flux,
                         yzmom_flux, zymom_flux,
                         xheat_flux, yheat_flux, zheat_flux,
                         xqv_flux, yqv_flux, zqv_flux,
                         z_phys, flux_comp);
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
ABLMost::compute_most_bcs (const int& lev,
                           const Vector<MultiFab*>& mfs,
                           MultiFab* xxmom_flux,
                           MultiFab* yymom_flux,
                           MultiFab* zzmom_flux,
                           MultiFab* xymom_flux, MultiFab* yxmom_flux,
                           MultiFab* xzmom_flux, MultiFab* zxmom_flux,
                           MultiFab* yzmom_flux, MultiFab* zymom_flux,
                           MultiFab* xheat_flux,
                           MultiFab* yheat_flux,
                           MultiFab* zheat_flux,
                           MultiFab* xqv_flux,
                           MultiFab* yqv_flux,
                           MultiFab* zqv_flux,
                           MultiFab* z_phys,
                           const FluxCalc& flux_comp)
{
    const int klo   = 0;
    const int icomp = 0;
    const auto& dxInv = m_geom[lev].InvCellSizeArray();
    const auto& dz_no_terrain = m_geom[lev].CellSize(2);
    for (MFIter mfi(*mfs[0]); mfi.isValid(); ++mfi)
    {
        // TODO: No LSM lateral ghost cells, should this change?
        // Valid CC box
        Box vbx = mfi.validbox(); vbx.makeSlab(2,klo-1);

        Box vbxx = surroundingNodes(vbx,0);
        Box vbxy = surroundingNodes(vbx,1);

        // Expose for GPU
        bool exp_most = m_exp_most;
        bool rot_most = m_rotate;

        // Get field arrays
        const auto cons_arr  = mfs[Vars::cons]->array(mfi);
        const auto velx_arr  = mfs[Vars::xvel]->array(mfi);
        const auto vely_arr  = mfs[Vars::yvel]->array(mfi);
        const auto velz_arr  = mfs[Vars::zvel]->array(mfi);

        // Explicit MOST vars
        auto t13_arr = (m_exp_most)               ? xzmom_flux->array(mfi) : Array4<Real>{};
        auto t31_arr = (m_exp_most && zxmom_flux) ? zxmom_flux->array(mfi) : Array4<Real>{};

        auto t23_arr = (m_exp_most)               ? yzmom_flux->array(mfi) : Array4<Real>{};
        auto t32_arr = (m_exp_most && zymom_flux) ? zymom_flux->array(mfi) : Array4<Real>{};

        auto hfx3_arr = (m_exp_most)             ? zheat_flux->array(mfi) : Array4<Real>{};
        auto qfx3_arr = (m_exp_most && zqv_flux) ? zqv_flux->array(mfi)   : Array4<Real>{};

        // Rotated MOST vars
        auto t11_arr = (m_rotate) ? xxmom_flux->array(mfi) : Array4<Real>{};
        auto t22_arr = (m_rotate) ? yymom_flux->array(mfi) : Array4<Real>{};
        auto t33_arr = (m_rotate) ? zzmom_flux->array(mfi) : Array4<Real>{};
        auto t12_arr = (m_rotate) ? xymom_flux->array(mfi) : Array4<Real>{};
        auto t21_arr = (m_rotate) ? yxmom_flux->array(mfi) : Array4<Real>{};

        auto hfx1_arr = (m_rotate) ? xheat_flux->array(mfi) : Array4<Real>{};
        auto hfx2_arr = (m_rotate) ? yheat_flux->array(mfi) : Array4<Real>{};
        auto qfx1_arr = (m_rotate && xqv_flux) ? xqv_flux->array(mfi) : Array4<Real>{};
        auto qfx2_arr = (m_rotate && yqv_flux) ? yqv_flux->array(mfi) : Array4<Real>{};

        // Viscosity and terrain
        const auto  eta_arr  = m_eddyDiffs_lev[lev]->array(mfi);
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

        // Get LSM fluxes
        auto lmask_arr    = (m_lmask_lev[lev][0])    ? m_lmask_lev[lev][0]->array(mfi) :
                                                       Array4<int> {};
        auto lsm_flux_arr = (m_lsm_flux_lev[lev][0]) ? m_lsm_flux_lev[lev][0]->array(mfi) :
                                                       Array4<Real> {};

        for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx)
        {
            const Box& bx = (*mfs[var_idx])[mfi].box();
            auto dest_arr = (*mfs[var_idx])[mfi].array();

            if (var_idx == Vars::cons) {
                Box b2d = bx;
                b2d.setBig(2,klo-1);
                int n = RhoTheta_comp;
                ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real dz  = (zphys_arr) ? ( zphys_arr(i,j,klo  ) - zphys_arr(i,j,klo-1) ) : dz_no_terrain;
                    Real dz1 = (zphys_arr) ? ( zphys_arr(i,j,klo+1) - zphys_arr(i,j,klo  ) ) : dz_no_terrain;
                    Real Tflux = flux_comp.compute_t_flux(i, j, k, n, icomp, dz, dz1, exp_most, eta_arr,
                                                          cons_arr, velx_arr, vely_arr,
                                                          umm_arr, tm_arr, u_star_arr, t_star_arr, t_surf_arr,
                                                          dest_arr);

                    if (rot_most) {
                        rotate_scalar_flux(i, j, klo, Tflux, dxInv, zphys_arr,
                                           hfx1_arr, hfx2_arr, hfx3_arr);
                    } else {
                        int is_land = (lmask_arr) ? lmask_arr(i,j,klo) : 1;
                        if (is_land && lsm_flux_arr && vbx.contains(i,j,k)) {
                            lsm_flux_arr(i,j,klo) = Tflux;
                        }
                        else if ((k == klo-1) && vbx.contains(i,j,k) && exp_most) {
                            hfx3_arr(i,j,klo) = Tflux;
                        }
                    }
                });

                // TODO: Generalize MOST q flux
                if ((flux_type == FluxCalcType::CUSTOM) && use_moisture) {
                    n = RhoQ1_comp;
                    ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        Real dz  = (zphys_arr) ? ( zphys_arr(i,j,klo  ) - zphys_arr(i,j,klo-1) ) : dz_no_terrain;
                        Real dz1 = (zphys_arr) ? ( zphys_arr(i,j,klo+1) - zphys_arr(i,j,klo  ) ) : dz_no_terrain;
                        Real Qflux = flux_comp.compute_q_flux(i, j, k, n, icomp, dz, dz1, exp_most, eta_arr,
                                                              cons_arr, velx_arr, vely_arr,
                                                              umm_arr, qm_arr, u_star_arr, q_star_arr, t_surf_arr,
                                                              dest_arr);

                        if (rot_most) {
                                rotate_scalar_flux(i, j, klo, Qflux, dxInv, zphys_arr,
                                                   qfx1_arr, qfx2_arr, qfx3_arr);
                        } else {
                            if ((k == klo-1) && vbx.contains(i,j,k) && exp_most) {
                                qfx3_arr(i,j,klo) = Qflux;
                            }
                        }
                    });
                }

            } else if (var_idx == Vars::xvel) {

                Box xb2d = surroundingNodes(bx,0);
                xb2d.setBig(2,klo-1);

                ParallelFor(xb2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real dz  = (zphys_arr) ? ( zphys_arr(i,j,klo  ) - zphys_arr(i,j,klo-1) ) : dz_no_terrain;
                    Real dz1 = (zphys_arr) ? ( zphys_arr(i,j,klo+1) - zphys_arr(i,j,klo  ) ) : dz_no_terrain;
                    Real stressx = flux_comp.compute_u_flux(i, j, k, icomp, dz, dz1, exp_most, eta_arr,
                                                            cons_arr, velx_arr, vely_arr,
                                                            umm_arr, um_arr, u_star_arr,
                                                            dest_arr);

                    if (rot_most) {
                        rotate_stress_tensor(i, j, klo, stressx, dxInv, zphys_arr,
                                             velx_arr, vely_arr, velz_arr,
                                             t11_arr, t22_arr, t33_arr,
                                             t12_arr, t21_arr,
                                             t13_arr, t31_arr,
                                             t23_arr, t32_arr);
                    } else {
                        if ((k == klo-1) && vbxx.contains(i,j,k) && exp_most) {
                            t13_arr(i,j,klo) = stressx;
                            if (t31_arr) t31_arr(i,j,klo) = stressx;
                        }
                    }
                });

            } else if (var_idx == Vars::yvel) {

                Box yb2d = surroundingNodes(bx,1);
                yb2d.setBig(2,klo-1);

                ParallelFor(yb2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real dz  = (zphys_arr) ? ( zphys_arr(i,j,klo  ) - zphys_arr(i,j,klo-1) ) : dz_no_terrain;
                    Real dz1 = (zphys_arr) ? ( zphys_arr(i,j,klo+1) - zphys_arr(i,j,klo  ) ) : dz_no_terrain;
                    Real stressy = flux_comp.compute_v_flux(i, j, k, icomp, dz, dz1, exp_most, eta_arr,
                                                            cons_arr, velx_arr, vely_arr,
                                                            umm_arr, vm_arr, u_star_arr,
                                                            dest_arr);

                    // NOTE: One stress rotation for ALL the stress components
                    if (!rot_most) {
                        if ((k == klo-1) && vbxy.contains(i,j,k) && exp_most) {
                            t23_arr(i,j,klo) = stressy;
                            if (t32_arr) t32_arr(i,j,klo) = stressy;
                        }
                    }
                });
            }
        } // var_idx
    } // mfiter
}

void
ABLMost::time_interp_sst (const int& lev,
                          const Real& time)
{
    // Time interpolation
    Real dT = m_bdy_time_interval;
    Real time_since_start = time - m_start_bdy_time;
    int n_time = static_cast<int>( time_since_start /  dT);
    Real alpha = (time_since_start - n_time * dT) / dT;
    AMREX_ALWAYS_ASSERT( alpha >= 0. && alpha <= 1.0);
    Real oma   = 1.0 - alpha;
    AMREX_ALWAYS_ASSERT( (n_time >= 0) && (n_time < (m_sst_lev[lev].size()-1)));

    // Populate t_surf
    for (MFIter mfi(*t_surf[lev]); mfi.isValid(); ++mfi)
    {
        Box gtbx = mfi.growntilebox();

        auto t_surf_arr = t_surf[lev]->array(mfi);
        const auto sst_hi_arr = m_sst_lev[lev][n_time+1]->const_array(mfi);
        const auto sst_lo_arr = m_sst_lev[lev][n_time  ]->const_array(mfi);
        auto lmask_arr  = (m_lmask_lev[lev][0]) ? m_lmask_lev[lev][0]->array(mfi) :
                                                  Array4<int> {};

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            int is_land = (lmask_arr) ? lmask_arr(i,j,k) : 1;
            if (!is_land) {
                t_surf_arr(i,j,k) = oma   * sst_lo_arr(i,j,k)
                                  + alpha * sst_hi_arr(i,j,k);
            }
        });
    }
}

void
ABLMost::get_lsm_tsurf (const int& lev)
{
    for (MFIter mfi(*t_surf[lev]); mfi.isValid(); ++mfi)
    {
        Box gtbx = mfi.growntilebox();

        // TODO: LSM does not carry lateral ghost cells.
        //       This copies the valid box into the ghost cells.
        //       Fillboundary is called after this to pick up the
        //       interior ghost and periodic directions. Is there
        //       a better approach?
        Box vbx  = mfi.validbox();
        int i_lo = vbx.smallEnd(0); int i_hi = vbx.bigEnd(0);
        int j_lo = vbx.smallEnd(1); int j_hi = vbx.bigEnd(1);

        auto t_surf_arr = t_surf[lev]->array(mfi);
        auto lmask_arr  = (m_lmask_lev[lev][0]) ? m_lmask_lev[lev][0]->array(mfi) :
                                                  Array4<int> {};
        const auto lsm_arr = m_lsm_data_lev[lev][0]->const_array(mfi);

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
ABLMost::update_pblh (const int& lev,
                      Vector<Vector<MultiFab>>& vars,
                       MultiFab* z_phys_cc)
{
    if (pblh_type == PBLHeightCalcType::MYNN25) {
        MYNNPBLH estimator;
        compute_pblh(lev, vars, z_phys_cc, estimator);
    } else if (pblh_type == PBLHeightCalcType::YSU) {
        amrex::Error("YSU PBLH calc not implemented yet");
    }
}

template <typename PBLHeightEstimator>
void
ABLMost::compute_pblh (const int& lev,
                       Vector<Vector<MultiFab>>& vars,
                       MultiFab* z_phys_cc,
                       const PBLHeightEstimator& est)
{
    const int klo = 0;
    const Real dz_no_terrain = m_geom[lev].CellSize(2);

    int moist_flag = 0;
    int n_qstate = vars[lev][Vars::cons].nComp() - (NVAR_max - NMOIST_max);
    if (n_qstate > 3) {
        moist_flag = (n_qstate > 3) ? RhoQ4_comp : RhoQ3_comp;
    } else if (n_qstate > 0) {
        moist_flag = 1;
    }

    for (MFIter mfi(*pblh[lev]); mfi.isValid(); ++mfi)
    {
        Box gtbx = mfi.growntilebox();

        const auto zphys_arr = (z_phys_cc) ? z_phys_cc->const_array(mfi) : Array4<const Real>{};

        const auto cons_arr   = vars[lev][Vars::cons].const_array(mfi);
        const auto u_star_arr = u_star[lev]->const_array(mfi);
        const auto t_star_arr = t_star[lev]->const_array(mfi);
              auto pblh_arr   = pblh[lev]->array(mfi);

        auto lmask_arr  = (m_lmask_lev[lev][0]) ? m_lmask_lev[lev][0]->const_array(mfi) : Array4<int> {};

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            pblh_arr(i,j,k) = est.compute_pblh(i,j,k,
                                               dz_no_terrain,zphys_arr,
                                               cons_arr,u_star_arr,t_star_arr,lmask_arr,
                                               moist_flag);
        });
    }
}

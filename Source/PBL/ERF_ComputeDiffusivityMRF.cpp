#include "ERF_SurfaceLayer.H"
#include "ERF_DirectionSelector.H"
#include "ERF_Diffusion.H"
#include "ERF_Constants.H"
#include "ERF_TurbStruct.H"
#include "ERF_PBLModels.H"
#include "ERF_TileNoZ.H"
#include "ERF_MoistUtils.H"

using namespace amrex;

void
ComputeDiffusivityMRF (const MultiFab& xvel,
                       const MultiFab& yvel,
                       const MultiFab& cons_in,
                       MultiFab& eddyViscosity,
                       const Geometry& geom,
                       const TurbChoice& turbChoice,
                       std::unique_ptr<SurfaceLayer>& SurfLayer,
                       bool use_terrain_fitted_coords,
                       bool use_moisture,
                       int level,
                       const BCRec* bc_ptr,
                       bool /*vert_only*/,
                       const std::unique_ptr<MultiFab>& z_phys_nd,
                       const MoistureComponentIndices& moisture_indices)
{
    /*
    Implementation of the older MRF Scheme based on Hong and Pan (1996)
    " Nonlocal Boundary Layer Vertical Diffusion in a Medium-Range Forecast
    Model"
    */

    // Domain extent in z-dir
    int klo = geom.Domain().smallEnd(2);
    int khi = geom.Domain().bigEnd(2);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(eddyViscosity, TileNoZ()); mfi.isValid(); ++mfi) {

        // Box operated on must span fill domain in z-dir
        const Box& gbx = mfi.growntilebox(IntVect(1,1,0));
        AMREX_ALWAYS_ASSERT( gbx.smallEnd(2) == klo &&
                             gbx.bigEnd(2)   == khi );

        //
        // Compute the height of the PBL without thermal excess
        // From Hong et al. 2006, Eqns. 1 & 2:
        //
        //   h = Rib_cf * theta_va * | U(h) |^2 / (g * (theta_v(h) - theta_s))
        //
        // where the surface virtual potential temperature is
        //
        //   theta_s = theta_va + theta_T
        //
        // and here the thermal excess theta_T = zero
        //

        // create flattened boxes to store PBL height
        const GeometryData gdata = geom.data();
        const Box xybx = PerpendicularBox<ZDir>(gbx, IntVect{0, 0, 0});
        FArrayBox pbl_height_predictor(xybx, 1, The_Async_Arena());
        FArrayBox pbl_height_corrector(xybx, 1, The_Async_Arena());
        IArrayBox pbl_index(xybx, 1, The_Async_Arena());
        const auto& pblh_arr      = pbl_height_predictor.array();
        const auto& pblh_corr_arr = pbl_height_corrector.array();
        const auto& pbli_arr      = pbl_index.array();

        // Get some data in arrays
        const auto& cell_data = cons_in.const_array(mfi);
        const auto& uvel = xvel.const_array(mfi);
        const auto& vvel = yvel.const_array(mfi);

        const Real Ribcr = turbChoice.pbl_mrf_Ribcr;
        //const Real f0 = turbChoice.pbl_mrf_coriolis_freq;
        const auto& u_star_arr = SurfLayer->get_u_star(level)->const_array(mfi);
        const auto& t_star_arr = SurfLayer->get_t_star(level)->const_array(mfi);
        const auto& l_obuk_arr = SurfLayer->get_olen(level)->const_array(mfi);
        const auto& t10av_arr  = SurfLayer->get_mac_avg(level, 2)->const_array(mfi);
        const auto& q10av_arr  = SurfLayer->get_mac_avg(level, 3)->const_array(mfi);
        const auto& q_surf_arr = SurfLayer->get_q_surf(level)->const_array(mfi);
        //const auto& t_surf_arr = SurfLayer->get_t_surf(level)->const_array(mfi);
        const Array4<Real const> z_nd_arr = z_phys_nd->array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            // note: if using a fixed surf_heat_flux (default), t_surf is not updated
            //const Real t_surf  = t_surf_arr(i, j, 0);
            const Real t_layer = t10av_arr(i, j, 0);

            Real zval;
            int kpbl = klo;
            bool above_critical = false;
            while (!above_critical && ((kpbl + 1) <= khi)) {
                kpbl += 1;

                // height above ground level
                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : (kpbl + myhalf) * gdata.CellSize(2);

                // Use virtual potential temperature for stability calculations
                const Real theta_v = GetThetav(i, j, kpbl, cell_data, moisture_indices);
                const Real moisture_fraction = use_moisture ? q10av_arr(i, j, 0) : zero;
                const Real t_layer_v = t_layer * (one + amrex::Real(0.61) * moisture_fraction);
                const Real ws2 = fourth * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                          (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                const Real Rib = CONST_GRAV * zval * (theta_v - t_layer_v) / (ws2 * t_layer_v);
                above_critical = (Rib >= Ribcr);
            }

            // Empirical expression for PBLH is given by h = c u* / f
            // Garratt (1994) and Tennekes (1982)
            // Also, c.f. Zilitinkevitch et al 2012 referenced in Pedersen et al. Real(2014.)
            //const Real c_pblh = (l_obuk_arr(i, j, 0) > 0) ? Real(0.16) : Real(0.60);
            //const Real pblh_emp = c_pblh * u_star_arr(i, j, 0) / f0;

            // Fallback to first cell
            const Real pblh_emp = (use_terrain_fitted_coords)
                                ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                                : myhalf * gdata.CellSize(2);

            // Initial PBL Height
            // Avoiding detailed interpolation here
            pblh_arr(i, j, 0) = (above_critical) ? zval : pblh_emp;
            pbli_arr(i, j, 0) = (above_critical) ? kpbl : klo+1;  // k < kpbl is considered the PBL
        });

        const auto& q_star_arr = SurfLayer->get_q_star(level)->const_array(mfi);

        //
        // Corrector PBL height: apply WRF MRF countergradient correction (HGAMT/HGAMQ)
        // following Hong and Pan (1996).
        //
        // WRF reference:
        //   https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F
        //
        // WRF constants (module_bl_mrf.F lines 300-311):
        //   CFAC   = 7.8  (= const_b, surface layer factor)
        //   GAMCRT = 3.0  (maximum heat countergradient, K)
        //   GAMCRQ = 2.e-3 (maximum moisture countergradient, kg/kg)
        //
        // WRF formulas (module_bl_mrf.F lines 872-879):
        //   GAMFAC = CFAC / (rho * wstar)
        //   HGAMT  = min(GAMFAC * HFX / CPM, GAMCRT)  =>  min(-b * u_* * t_* / w_*, GAMCRT)
        //   HGAMQ  = min(GAMFAC * QFX,       GAMCRQ)  =>  min(-b * u_* * q_* / w_*, GAMCRQ)
        //   VPERT  = HGAMT + EP1 * theta * HGAMQ
        //   THERMAL += max(VPERT, 0)   (raises effective surface virtual pot. temp.)
        //
        // WRF diffusion coefficients below PBL (module_bl_mrf.F lines 968-986):
        //   K_mom   = rho * wstar * kappa * z * (1 - z/h)^2
        //   K_theta = K_mom / Prt
        //   K_q     = K_mom / Prq  (Prq ~ Prt for moisture)
        //
        const Real const_b = turbChoice.pbl_mrf_const_b;
        const Real sf = turbChoice.pbl_mrf_sf;
        constexpr Real prmin = Real(0.5);
        constexpr Real prmax = Real(4.0);
        constexpr Real GAMCRT = Real(3.0);    // WRF GAMCRT
        constexpr Real GAMCRQ = Real(2.e-3);  // WRF GAMCRQ
        const bool enable_mrf_countergradient = turbChoice.enable_mrf_countergradient;
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            const Real t_layer  = t10av_arr(i, j, 0);
            const Real moisture_fraction = use_moisture ? q10av_arr(i, j, 0) : zero;
            const Real t_layer_v = t_layer * (one + amrex::Real(0.61) * moisture_fraction);
            const Real phiM     = (l_obuk_arr(i, j, 0) > 0)
                                ? (1 + 5 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0))
                                : std::pow(
                                           (1 - 8 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0)),
                                           -one / three);
            const Real wstar    = u_star_arr(i, j, 0) / phiM;

            // WRF MRF (module_bl_mrf.F lines 872-879):
            //   Simple thermal excess: HGAMT = min(-b * u_* * t_* / w_*, GAMCRT)
            //   Moisture excess:       HGAMQ = min(-b * u_* * q_* / w_*, GAMCRQ)
            //   Virtual pot. temp. excess: VPERT = HGAMT + EP1 * theta * HGAMQ
            //   THERMAL += max(VPERT, 0)  =>  surface virtual pot. temp. is raised
            const Real HGAMT = amrex::min(-const_b * u_star_arr(i, j, 0) * t_star_arr(i, j, 0) / wstar, GAMCRT);
            Real HGAMQ = zero;
            if (use_moisture && enable_mrf_countergradient) {
                HGAMQ = amrex::min(-const_b * u_star_arr(i, j, 0) * q_star_arr(i, j, 0) / wstar, GAMCRQ);
            }
            // Virtual potential temperature excess at surface (positive = more unstable)
            // EP1 = R_v/R_d - 1 ~ 0.61
            const Real VPERT = amrex::max(HGAMT + amrex::Real(0.61) * t_layer * HGAMQ, zero);
            // Effective surface virtual potential temperature used in PBL height finding
            // (WRF: THERMAL = theta_v_surface + VPERT)
            const Real t_surf_v = t_layer_v + VPERT;

            int kpbl = klo;
            Real zval0, zval, Rib0, Rib;
            {
                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : (kpbl + myhalf) * gdata.CellSize(2);
                const Real theta_v = GetThetav(i, j, kpbl, cell_data, moisture_indices);
                const Real ws2 = fourth * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                          (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                Rib = CONST_GRAV * zval * (theta_v - t_surf_v) / (ws2 * t_layer_v);
            }

            bool above_critical = false;
            while (!above_critical && ((kpbl + 1) <= khi)) {
                zval0 = zval;
                Rib0 = Rib;
                kpbl += 1;

                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : (kpbl + myhalf) * gdata.CellSize(2);
                const Real theta_v = GetThetav(i, j, kpbl, cell_data, moisture_indices);
                const Real ws2 = fourth * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                          (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                Rib = CONST_GRAV * zval * (theta_v - t_surf_v) / (ws2 * t_layer_v);
                above_critical = (Rib >= Ribcr);
            }

            if (above_critical) {
                // Interpolate to height at which Rib == Ribcr
                Real pblh_interp = zval0 + (zval - zval0) / (Rib - Rib0) * (Ribcr - Rib0);
                pblh_corr_arr(i, j, 0) = pblh_interp;
                     pbli_arr(i, j, 0) = kpbl;  // k < kpbl is considered the PBL
            } else {
                // Fallback to first cell
                pblh_corr_arr(i, j, 0) = (use_terrain_fitted_coords)
                                       ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                                       : myhalf * gdata.CellSize(2);
                     pbli_arr(i, j, 0) = klo + 1;
            }

        });
        /*
          amrex::Print() << "PBL height computed for MRF scheme at level "
          << pblh_arr(2, 2, 0) << "  " << pblh_corr_arr(2, 2, 0)
          << std::endl;
          amrex::Print() << "PBL Temp:" << t_surf_arr(2, 2, 0) << "  "
          << t10av_arr(2, 2, 0) << std::endl;
        */

        // -- Compute diffusion coefficients --

        const Array4<Real>& K_turb = eddyViscosity.array(mfi);

        // Dirichlet flags to switch derivative stencil
        bool c_ext_dir_on_zlo = ((bc_ptr[BCVars::cons_bc].lo(2) == ERFBCType::ext_dir));
        bool c_ext_dir_on_zhi = ((bc_ptr[BCVars::cons_bc].hi(2) == ERFBCType::ext_dir));
        bool u_ext_dir_on_zlo = ((bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir));
        bool u_ext_dir_on_zhi = ((bc_ptr[BCVars::xvel_bc].hi(2) == ERFBCType::ext_dir));
        bool v_ext_dir_on_zlo = ((bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir));
        bool v_ext_dir_on_zhi = ((bc_ptr[BCVars::yvel_bc].hi(2) == ERFBCType::ext_dir));

        const auto& dxInv = geom.InvCellSizeArray();
        const Real dz_inv = geom.InvCellSize(2);
        const int izmin   = geom.Domain().smallEnd(2);
        const int izmax   = geom.Domain().bigEnd(2);

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real zval = (use_terrain_fitted_coords)
                            ? Compute_Zrel_AtCellCenter(i, j, k, z_nd_arr)
                            : (k + myhalf) * gdata.CellSize(2);
            const Real rho = cell_data(i, j, k, Rho_comp);
            const Real met_h_zeta = (use_terrain_fitted_coords)
                                  ? Compute_h_zeta_AtCellCenter(i, j, k, dxInv, z_nd_arr) : one;
            const Real dz_terrain = met_h_zeta / dz_inv;
            if (k < pbli_arr(i, j, 0)) {
                const Real phiM = (l_obuk_arr(i, j, 0) > 0)
                                ? (1 + 5 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0))
                                : std::pow(
                                           (1 - 8 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0)),
                                           -one / three);
                const Real phit = (l_obuk_arr(i, j, 0) > 0)
                                ? (1 + 5 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0))
                                : std::pow(
                                           (1 - 16 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0)),
                                           -one / two);
                const Real Prt = amrex::min(amrex::max(phit / phiM + const_b * KAPPA * sf, prmin), prmax);
                const Real wstar = u_star_arr(i, j, 0) / phiM;
                K_turb(i, j, k, EddyDiff::Mom_v) = rho * wstar * KAPPA * zval
                                                 * (1 - zval / pblh_corr_arr(i, j, 0))
                                                 * (1 - zval / pblh_corr_arr(i, j, 0));
                K_turb(i, j, k, EddyDiff::Theta_v) = K_turb(i, j, k, EddyDiff::Mom_v) / Prt;

                // Moisture diffusivity: WRF MRF uses Prq ~ Prt (same stability functions
                // for heat and moisture, module_bl_mrf.F lines 968-986).
                // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F
                if (turbChoice.mrf_moistvars) {
                    const Real Prq = amrex::min(amrex::max(phit / phiM + const_b * KAPPA * sf, prmin), prmax);
                    K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Mom_v) / Prq;
                } else {
                    K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
                }
            } else {
                // Free atmosphere above PBL
                const Real lambda = Real(150.0);
                const Real lscale = (KAPPA * zval * lambda) / (KAPPA * zval + lambda);
                Real dthetadz, dudz, dvdz;
                ComputeVerticalDerivativesPBL(i, j, k, uvel, vvel, cell_data, izmin, izmax, one / dz_terrain,
                                              c_ext_dir_on_zlo, c_ext_dir_on_zhi, u_ext_dir_on_zlo,
                                              u_ext_dir_on_zhi, v_ext_dir_on_zlo, v_ext_dir_on_zhi, dthetadz,
                                              dudz, dvdz, moisture_indices);
                const Real wind_shear = dudz * dudz + dvdz * dvdz + Real(1.0e-9);

                // Use virtual potential temperature for Richardson number stability calculation
                const Real theta_v = GetThetav(i, j, k, cell_data, moisture_indices);
                const Real dthetav_dz = myhalf * (GetThetav(i, j, k+1, cell_data, moisture_indices) -
                                                   GetThetav(i, j, k-1, cell_data, moisture_indices)) * (one / dz_terrain);
                Real grad_Ri = CONST_GRAV / theta_v * dthetav_dz / wind_shear;
                grad_Ri = std::max(grad_Ri, -Real(100.0));  // Hong et al. 2006, MWR, Appendix A
                /*
                  const Real Pr = Real(1.5) + Real(3.08) * grad_Ri;
                  const Real fm =
                  (grad_Ri > 0)
                  ? (std::exp(-Real(8.5) * grad_Ri) + (Real(0.15) / (grad_Ri + three)) * Pr)
                  : std::pow((1 - 12 * grad_Ri), -one / three);
                  const Real ft =
                  (grad_Ri > 0)
                  ? (std::exp(-Real(8.5) * grad_Ri) + (Real(0.15) / (grad_Ri + three)))
                  : std::pow((1 - 16 * grad_Ri), -one / two);
                */
                // Using YSU stability functions (Hong et al. 2006, MWR, Appendix A)
                Real Pr = one + Real(2.1) * grad_Ri;  // Eqn. A19
                const Real fm = (grad_Ri > 0)
                              ? one / ((one + Real(5.0) * grad_Ri) * (one + Real(5.0) * grad_Ri))
                              : 1 - 8 * grad_Ri / (1 + Real(1.746) * std::sqrt(-grad_Ri)); // Eqn. A20b
                const Real ft = (grad_Ri > 0)
                              ? one / ((one + Real(5.0) * grad_Ri) * (one + Real(5.0) * grad_Ri))
                              : 1 - 8 * grad_Ri / (1 + Real(1.286) * std::sqrt(-grad_Ri)); // Eqn. A20a
                const Real rl2wsp = rho * lscale * lscale * std::sqrt(wind_shear);

                Pr = std::max(amrex::Real(0.25), std::min(Pr, Real(4.0)));  // Hong et al. 2006, MWR, Appendix A

                K_turb(i, j, k, EddyDiff::Mom_v)   = rl2wsp * fm * Pr;
                K_turb(i, j, k, EddyDiff::Theta_v) = rl2wsp * ft;
                // WRF MRF: moisture diffusivity matches heat in free atmosphere
                K_turb(i, j, k, EddyDiff::Q_v) = turbChoice.mrf_moistvars
                                                ? rl2wsp * ft
                                                : K_turb(i, j, k, EddyDiff::Theta_v);
            }

            // limit both diffusion coefficients
#if 0
            // Hong et al. 2006, MWR, Appendix A
            constexpr Real ckz  = Real(0.001);
            constexpr Real Kmax = Real(1000.0);
            const Real rhoKmin  = ckz * dz_terrain * rho;
            const Real rhoKmax  = rho * Kmax;
#endif
            // Hong & Pan 1996, MWR
            constexpr Real Kmin = Real(0.1);
            constexpr Real Kmax = Real(300.0);
            const Real rhoKmin  = rho * Kmin;
            const Real rhoKmax  = rho * Kmax;

            K_turb(i, j, k, EddyDiff::Mom_v) = std::max(
                std::min(K_turb(i, j, k, EddyDiff::Mom_v), rhoKmax), rhoKmin);
            K_turb(i, j, k, EddyDiff::Theta_v) = std::max(
                std::min(K_turb(i, j, k, EddyDiff::Theta_v), rhoKmax), rhoKmin);
            K_turb(i, j, k, EddyDiff::Q_v) = std::max(
                std::min(K_turb(i, j, k, EddyDiff::Q_v), rhoKmax), rhoKmin);
            K_turb(i, j, k, EddyDiff::Turb_lengthscale) = pblh_arr(i, j, 0);
        });

        // FOEXTRAP top and bottom ghost cells
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
        {
            K_turb(i, j, klo-1, EddyDiff::Mom_v  ) = K_turb(i, j, klo, EddyDiff::Mom_v  );
            K_turb(i, j, klo-1, EddyDiff::Theta_v) = K_turb(i, j, klo, EddyDiff::Theta_v);
            K_turb(i, j, klo-1, EddyDiff::Q_v    ) = K_turb(i, j, klo, EddyDiff::Q_v    );
            K_turb(i, j, khi+1, EddyDiff::Mom_v  ) = K_turb(i, j, khi, EddyDiff::Mom_v  );
            K_turb(i, j, khi+1, EddyDiff::Theta_v) = K_turb(i, j, khi, EddyDiff::Theta_v);
            K_turb(i, j, khi+1, EddyDiff::Q_v    ) = K_turb(i, j, khi, EddyDiff::Q_v    );
        });
    }// mfi
}

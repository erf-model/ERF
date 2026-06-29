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
                       bool /*use_moisture*/,
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
        // Step 1: Compute the height of the PBL without thermal excess
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
                const Real t_layer_v = t_layer * (one + amrex::Real(0.61) * zero);  // t_layer is already theta, convert to theta_v (assuming dry for now)
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

        //
        // Step 2: Corrector PBL height for thermal excess
        // where
        //
        //   theta_T = b * surf_temp_flux / w_star
        //
        const Real const_b = turbChoice.pbl_mrf_const_b;
        const Real sf = turbChoice.pbl_mrf_sf;
        const bool enable_mrf_countergradient = turbChoice.enable_mrf_countergradient;
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            const Real t_layer  = t10av_arr(i, j, 0);
            const Real phiM     = (l_obuk_arr(i, j, 0) > 0)
                                ? (1 + 5 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0))
                                : std::pow(
                                           (1 - 8 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0)),
                                           -one / three);
            const Real wstar    = u_star_arr(i, j, 0) / phiM;
            const Real t_excess = -const_b * u_star_arr(i, j, 0) * t_star_arr(i, j, 0) / wstar;
            const Real t_surf   = t_layer + amrex::max(amrex::min(t_excess, amrex::Real(3)), amrex::Real(0));

            int kpbl = klo;
            Real zval0, zval, Rib0, Rib;
            {
                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : (kpbl + myhalf) * gdata.CellSize(2);
                // Use virtual potential temperature for stability calculations
                const Real theta_v = GetThetav(i, j, kpbl, cell_data, moisture_indices);
                const Real t_surf_v = t_surf * (one + amrex::Real(0.61) * zero);  // t_surf is already theta, convert to theta_v
                const Real ws2 = fourth * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                          (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                Rib = CONST_GRAV * zval * (theta_v - t_surf_v) / (ws2 * t_layer);
            }

            bool above_critical = false;
            while (!above_critical && ((kpbl + 1) <= khi)) {
                zval0 = zval;
                Rib0 = Rib;
                kpbl += 1;

                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : (kpbl + myhalf) * gdata.CellSize(2);
                // Use virtual potential temperature for stability calculations
                const Real theta_v = GetThetav(i, j, kpbl, cell_data, moisture_indices);
                const Real t_surf_v = t_surf * (one + amrex::Real(0.61) * zero);  // t_surf is already theta, convert to theta_v
                const Real ws2 = fourth * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                          (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                Rib = CONST_GRAV * zval * (theta_v - t_surf_v) / (ws2 * t_layer);
                above_critical = (Rib >= Ribcr);
            }

            if (above_critical) {
                // Interpolate to height at which Rib == Ribcr
                Real pblh_interp = zval0 + (zval - zval0) / (Rib - Rib0) * (Ribcr - Rib0);
                pblh_corr_arr(i, j, 0) = pblh_interp;
                     pbli_arr(i, j, 0) = kpbl;  // k < kpbl is considered the PBL
            } else {
                // Empirical expression for PBLH is given by h = c u* / f
                // Garratt (1994) and Tennekes (1982)
                // Also, c.f. Zilitinkevitch et al 2012 referenced in Pedersen et al. Real(2014.)
                //const Real c_pblh = (l_obuk_arr(i, j, 0) > 0) ? Real(0.16) : Real(0.60);
                //const Real pblh_emp = c_pblh * u_star_arr(i, j, 0) / f0;

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
                const Real Prt = phit / phiM + const_b * KAPPA * sf;
                const Real wstar = u_star_arr(i, j, 0) / phiM;
                K_turb(i, j, k, EddyDiff::Mom_v) = rho * wstar * KAPPA * zval
                                                 * (1 - zval / pblh_corr_arr(i, j, 0))
                                                 * (1 - zval / pblh_corr_arr(i, j, 0));
                K_turb(i, j, k, EddyDiff::Theta_v) = K_turb(i, j, k, EddyDiff::Mom_v) / Prt;
                
                // Compute countergradient correction term if enabled: γ_c = b * (u_* * θ_*) / w_s
                if (enable_mrf_countergradient) {
                    // Note: wstar is already u_* / φ_m, so we don't need to divide again
                    // Separate countergradient terms for different variables
                    const Real countergradient_mom = const_b * u_star_arr(i, j, 0) * t_star_arr(i, j, 0) / wstar;
                    
                    // Heat countergradient: γ_θ = b * (u_* * θ_*) / w_s
                    K_turb(i, j, k, EddyDiff::CounterGradient_theta) = countergradient_mom;
                    
                    // Moisture countergradient: γ_q = b * (u_* * q_*) / w_s (typically similar magnitude to heat)
                    // Using same coefficient for now; can be differentiated if needed
                    K_turb(i, j, k, EddyDiff::CounterGradient_q) = countergradient_mom;
                    
                    // Momentum countergradient (legacy)
                    K_turb(i, j, k, EddyDiff::CounterGradient_v) = countergradient_mom;
                } else {
                    K_turb(i, j, k, EddyDiff::CounterGradient_v) = zero;
                    K_turb(i, j, k, EddyDiff::CounterGradient_theta) = zero;
                    K_turb(i, j, k, EddyDiff::CounterGradient_q) = zero;
                }
                K_turb(i, j, k, EddyDiff::CounterGradient_h) = zero;
                
                // Initialize entrainment terms (zero in mixed layer)
                K_turb(i, j, k, EddyDiff::Entrainment_mom) = zero;
                K_turb(i, j, k, EddyDiff::Entrainment_theta) = zero;
                K_turb(i, j, k, EddyDiff::Entrainment_q) = zero;
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
                // Using YSU model instead of MRF model
                Real Pr = one + Real(2.1) * grad_Ri;  // Hong et al. 2006, MWR, Eqn. A19
                const Real fm = (grad_Ri > 0)
                              ? one / ((one + Real(5.0) * grad_Ri) * (one + Real(5.0) * grad_Ri))
                              : 1 - 8 * grad_Ri / (1 + Real(1.746) * std::sqrt(-grad_Ri)); // Hong et al. 2006, MWR, Eqn. A20b
                const Real ft = (grad_Ri > 0)
                              ? one / ((one + Real(5.0) * grad_Ri) * (one + Real(5.0) * grad_Ri))
                              : 1 - 8 * grad_Ri / (1 + Real(1.286) * std::sqrt(-grad_Ri)); // Hong et al. 2006, MWR, Eqn. A20a
                const Real rl2wsp = rho * lscale * lscale * std::sqrt(wind_shear);

                Pr = std::max(amrex::Real(0.25), std::min(Pr, Real(4.0)));  // Hong et al. 2006, MWR, Appendix A

                K_turb(i, j, k, EddyDiff::Mom_v)   = rl2wsp * fm * Pr;
                K_turb(i, j, k, EddyDiff::Theta_v) = rl2wsp * ft;
                
                // Check if we're at PBL top and add entrainment fluxes
                int kpbl = pbli_arr(i, j, 0);
                if (k == kpbl && kpbl < khi) {
                    // Entrainment at PBL top: gradual mixing zone instead of sharp discontinuity
                    // Entrainment is strongest at the PBL interface
                    const Real entrainment_coeff = Real(0.1);  // Typical entrainment coefficient
                    const Real u_star = u_star_arr(i, j, 0);
                    
                    // Entrainment momentum flux: ρ * w_e * Δu where w_e ~ α * u_*
                    K_turb(i, j, k, EddyDiff::Entrainment_mom) = rho * entrainment_coeff * u_star * KAPPA * zval;
                    
                    // Entrainment temperature flux: similar to momentum
                    K_turb(i, j, k, EddyDiff::Entrainment_theta) = K_turb(i, j, k, EddyDiff::Entrainment_mom) / ft;
                    
                    // Entrainment moisture flux
                    K_turb(i, j, k, EddyDiff::Entrainment_q) = K_turb(i, j, k, EddyDiff::Entrainment_theta);
                } else {
                    // No entrainment correction in free atmosphere away from PBL top
                    K_turb(i, j, k, EddyDiff::Entrainment_mom) = zero;
                    K_turb(i, j, k, EddyDiff::Entrainment_theta) = zero;
                    K_turb(i, j, k, EddyDiff::Entrainment_q) = zero;
                }
                 
                // No countergradient correction in free atmosphere
                K_turb(i, j, k, EddyDiff::CounterGradient_v) = zero;
                K_turb(i, j, k, EddyDiff::CounterGradient_h) = zero;
                K_turb(i, j, k, EddyDiff::CounterGradient_theta) = zero;
                K_turb(i, j, k, EddyDiff::CounterGradient_q) = zero;
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
            K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
            K_turb(i, j, k, EddyDiff::Turb_lengthscale) = pblh_arr(i, j, 0);
            // CounterGradient terms already set above in mixed layer or free atmosphere
        });

        // FOEXTRAP top and bottom ghost cells
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
        {
            K_turb(i, j, klo-1, EddyDiff::Mom_v  ) = K_turb(i, j, klo, EddyDiff::Mom_v  );
            K_turb(i, j, klo-1, EddyDiff::Theta_v) = K_turb(i, j, klo, EddyDiff::Theta_v);
            K_turb(i, j, klo-1, EddyDiff::Q_v    ) = K_turb(i, j, klo, EddyDiff::Q_v    );
            K_turb(i, j, klo-1, EddyDiff::CounterGradient_v) = K_turb(i, j, klo, EddyDiff::CounterGradient_v);
            K_turb(i, j, klo-1, EddyDiff::CounterGradient_h) = K_turb(i, j, klo, EddyDiff::CounterGradient_h);
            K_turb(i, j, klo-1, EddyDiff::CounterGradient_theta) = K_turb(i, j, klo, EddyDiff::CounterGradient_theta);
            K_turb(i, j, klo-1, EddyDiff::CounterGradient_q) = K_turb(i, j, klo, EddyDiff::CounterGradient_q);
            K_turb(i, j, klo-1, EddyDiff::Entrainment_mom) = K_turb(i, j, klo, EddyDiff::Entrainment_mom);
            K_turb(i, j, klo-1, EddyDiff::Entrainment_theta) = K_turb(i, j, klo, EddyDiff::Entrainment_theta);
            K_turb(i, j, klo-1, EddyDiff::Entrainment_q) = K_turb(i, j, klo, EddyDiff::Entrainment_q);
            K_turb(i, j, khi+1, EddyDiff::Mom_v  ) = K_turb(i, j, khi, EddyDiff::Mom_v  );
            K_turb(i, j, khi+1, EddyDiff::Theta_v) = K_turb(i, j, khi, EddyDiff::Theta_v);
            K_turb(i, j, khi+1, EddyDiff::Q_v    ) = K_turb(i, j, khi, EddyDiff::Q_v    );
            K_turb(i, j, khi+1, EddyDiff::CounterGradient_v) = K_turb(i, j, khi, EddyDiff::CounterGradient_v);
            K_turb(i, j, khi+1, EddyDiff::CounterGradient_h) = K_turb(i, j, khi, EddyDiff::CounterGradient_h);
            K_turb(i, j, khi+1, EddyDiff::CounterGradient_theta) = K_turb(i, j, khi, EddyDiff::CounterGradient_theta);
            K_turb(i, j, khi+1, EddyDiff::CounterGradient_q) = K_turb(i, j, khi, EddyDiff::CounterGradient_q);
            K_turb(i, j, khi+1, EddyDiff::Entrainment_mom) = K_turb(i, j, khi, EddyDiff::Entrainment_mom);
            K_turb(i, j, khi+1, EddyDiff::Entrainment_theta) = K_turb(i, j, khi, EddyDiff::Entrainment_theta);
            K_turb(i, j, khi+1, EddyDiff::Entrainment_q) = K_turb(i, j, khi, EddyDiff::Entrainment_q);
        });
    }// mfi
}

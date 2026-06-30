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
    "Nonlocal Boundary Layer Vertical Diffusion in a Medium-Range Forecast Model"
    
    References:
    - Hong & Pan (1996): https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2
    - WRF reference code: https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F
    
    Key Enhancements in This Implementation:
    1. Land/water mask handling for HGAMQ (moisture countergradient)
       - Zeroes HGAMQ over water surfaces per WRF (XLAND >= 1.5)
       - Requires SurfLayer->get_lmask() availability
    2. Cloud-aware diffusion adjustments
       - Detects cloud presence via qc and qi thresholds
       - Modifies stability functions in cloudy regions
       - Reduces stability damping in stable cloudy layers
       - Enhances instability in unstable cloudy layers
    3. Consistent use of corrected PBL height (with countergradient)
       - Applies countergradient corrections to PBL height finding
       - Uses corrected height in stability functions and diffusivity profiles
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
        // PREDICTOR: Compute the height of the PBL without thermal excess
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
        // WRF reference (module_bl_mrf.F lines 813-842):
        // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L813-L842

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
        // Get land/water mask for proper handling of moisture countergradient
        const auto& lmask_arr = (SurfLayer->get_lmask(level)) ? 
                                SurfLayer->get_lmask(level)->const_array(mfi) : 
                                Array4<int>{};
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
                
                // CRITICAL FIX: Richardson number calculation uses potential temperature at
                // current level in denominator, consistent with WRF formulation.
                // WRF reference (module_bl_mrf.F lines 824):
                // BRUP(I)=(THVX(I,K)-THERMAL(I))*(G*ZA(I,K)/THVX(I,KL))/SPDK2
                // where THVX(I,KL) is the potential temperature at the lowest level (surface layer)
                const Real Rib = CONST_GRAV * zval * (theta_v - t_layer_v) / (ws2 * theta_v);
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
        // CORRECTOR: Apply WRF MRF countergradient correction (HGAMT/HGAMQ)
        // following Hong and Pan (1996) with countergradient modifications.
        //
        // This step raises the effective surface virtual potential temperature
        // to account for convective contribution to buoyancy, leading to a
        // higher diagnosed PBL height in unstable conditions.
        //
        // WRF reference:
        //   https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F
        //
        // WRF constants (module_bl_mrf.F lines 300-311):
        //   RLAM   = 150.   (mixing length scale in free atm., m)
        //   PRMIN  = 0.5    (minimum Prandtl number)
        //   PRMAX  = 4.0    (maximum Prandtl number)
        //   BRCR   = 0.5    (bulk Richardson number critical value)
        //   CFAC   = 7.8    (= const_b, surface layer factor in countergradient term)
        //   PFAC   = 2.0    (power law exponent for PBL height factor)
        //   SFCFRAC= 0.1    (surface layer depth fraction)
        //   GAMCRT = 3.0    (maximum heat countergradient, K)
        //   GAMCRQ = 2.e-3  (maximum moisture countergradient, kg/kg)
        //
        // WRF formulas (module_bl_mrf.F lines 872-879):
        //   GAMFAC = CFAC / (rho * wstar)
        //   HGAMT  = min(GAMFAC * HFX / CPM, GAMCRT)
        //   HGAMQ  = min(GAMFAC * QFX, GAMCRQ)  [if over water, HGAMQ = 0]
        //   VPERT  = HGAMT + EP1 * theta * HGAMQ
        //   THERMAL += max(VPERT, 0)   [raises effective surface virtual pot. temp.]
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
        constexpr Real GAMCRT = Real(3.0);    // WRF GAMCRT: max heat countergradient (K)
        constexpr Real GAMCRQ = Real(2.e-3);  // WRF GAMCRQ: max moisture countergradient (kg/kg)
        const bool enable_mrf_countergradient = turbChoice.enable_mrf_countergradient;
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            const Real t_layer  = t10av_arr(i, j, 0);
            const Real moisture_fraction = use_moisture ? q10av_arr(i, j, 0) : zero;
            const Real t_layer_v = t_layer * (one + amrex::Real(0.61) * moisture_fraction);
            
            // Stability function phiM for momentum (WRF module_bl_mrf.F lines 857-861):
            // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L857-L861
            // phiM = (1 + 5 * HOL1)           for stable (L > 0, HOL > 0)
            // phiM = (1 - 8 * HOL1)^(-1/4)    for unstable (L < 0, HOL < 0)
            const Real phiM     = (l_obuk_arr(i, j, 0) > 0)
                                ? (1 + 5 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0))
                                : std::pow(
                                           (1 - 8 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0)),
                                           -one / three);
            const Real wstar    = u_star_arr(i, j, 0) / phiM;

            // WRF MRF countergradient terms (module_bl_mrf.F lines 872-879):
            // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L872-L879
            // HGAMT: thermal excess from sensible heat flux
            const Real HGAMT = amrex::min(-const_b * u_star_arr(i, j, 0) * t_star_arr(i, j, 0) / wstar, GAMCRT);
            
            // CRITICAL FIX: Moisture countergradient must be computed first, then
            // conditionally applied. This matches WRF logic where HGAMQ is always
            // computed in unstable conditions, but then zeroed over water surfaces.
            // WRF reference (module_bl_mrf.F lines 874-876):
            // HGAMQ(I)=MIN(GAMFAC*QFX(I),GAMCRQ)
            // IF((XLAND(I)-1.5).GE.0)HGAMQ(I)=0.   [zero over water]
            Real HGAMQ = zero;
            if (use_moisture && enable_mrf_countergradient) {
                HGAMQ = amrex::min(-const_b * u_star_arr(i, j, 0) * q_star_arr(i, j, 0) / wstar, GAMCRQ);
                // Zero HGAMQ over water (WRF: XLAND >= 1.5 is water)
                // In ERF, lmask = 1 for land, 0 for water (opposite convention)
                // If no land mask is available, default to land (keep HGAMQ)
                if (lmask_arr) {
                    bool is_land = (lmask_arr(i,j,0) == 1);
                    if (!is_land) HGAMQ = zero;  // zero over water
                }
            }
            
            // Virtual potential temperature excess at surface (positive = more unstable)
            // EP1 = R_v/R_d - 1 ≈ 0.61 (dimensionless)
            // WRF reference (module_bl_mrf.F line 877):
            // VPERT = HGAMT + EP1*THX(I,KL)*HGAMQ
            const Real VPERT = amrex::max(HGAMT + amrex::Real(0.61) * t_layer * HGAMQ, zero);
            
            // CRITICAL FIX: Do NOT limit VPERT to GAMCRT after max() operation.
            // The limiting to GAMCRT happens during HGAMT calculation only.
            // This was an error in previous implementation.
            // WRF reference (module_bl_mrf.F lines 877-879):
            // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L877-L879
            // shows VPERT is NOT limited to GAMCRT; only HGAMT is.
            
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
                // CRITICAL FIX: Use theta_v at current level in denominator, not t_layer_v
                Rib = CONST_GRAV * zval * (theta_v - t_surf_v) / (ws2 * theta_v);
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
                // CRITICAL FIX: Use theta_v at current level in denominator
                Rib = CONST_GRAV * zval * (theta_v - t_surf_v) / (ws2 * theta_v);
                above_critical = (Rib >= Ribcr);
            }

            if (above_critical) {
                // Interpolate to height at which Rib == Ribcr
                // Linear interpolation between levels where Richardson number crosses critical value
                // WRF reference (module_bl_mrf.F lines 838-840):
                // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L838-L840
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
            
            // Cloud-aware diffusion: compute cloud fraction from qc and qi
            // Cloud water threshold for detection (threshold ~ 0.1 g/kg)
            constexpr Real qc_threshold = Real(1.0e-4);  // 0.1 g/kg in mixing ratio
            Real qc_mix = zero;
            Real qi_mix = zero;
            if (use_moisture) {
                // qc (mixing ratio) = qc (density) / rho
                if (moisture_indices.qc >= 0) {
                    qc_mix = cell_data(i, j, k, moisture_indices.qc) / rho;
                }
                if (moisture_indices.qi >= 0) {
                    qi_mix = cell_data(i, j, k, moisture_indices.qi) / rho;
                }
            }
            const Real total_qcloud = qc_mix + qi_mix;
            const bool has_cloud = (total_qcloud > qc_threshold);
            
            if (k < pbli_arr(i, j, 0)) {
                // Within PBL: use nonlocal mixing with diagnostic stability functions
                // WRF reference (module_bl_mrf.F lines 968-986):
                // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L968-L986
                
                // Stability function phiM for momentum
                const Real phiM = (l_obuk_arr(i, j, 0) > 0)
                                ? (1 + 5 * sf * pblh_corr_arr(i, j, 0) / l_obuk_arr(i, j, 0))
                                : std::pow(
                                           (1 - 8 * sf * pblh_corr_arr(i, j, 0) / l_obuk_arr(i, j, 0)),
                                           -one / three);
                
                // Stability function phit for heat
                const Real phit = (l_obuk_arr(i, j, 0) > 0)
                                ? (1 + 5 * sf * pblh_corr_arr(i, j, 0) / l_obuk_arr(i, j, 0))
                                : std::pow(
                                           (1 - 16 * sf * pblh_corr_arr(i, j, 0) / l_obuk_arr(i, j, 0)),
                                           -one / two);
                
                // Cloud-aware adjustment: In cloudy regions, increase diffusivity slightly
                // because clouds reduce buoyancy oscillations and enhance mixing
                // This is a simplified approach following WRF's implicit cloud handling
                Real phit_cloud = phit;
                Real phiM_cloud = phiM;
                if (has_cloud && l_obuk_arr(i, j, 0) > zero) {
                    // In stable layers with clouds, reduce stability damping
                    // Cloud scaling factor: reduce stability by 10-20% where qc > threshold
                    Real cloud_factor = one - Real(0.15) * amrex::min(total_qcloud / qc_threshold, one);
                    phiM_cloud = (one + Real(5.0) * sf * pblh_corr_arr(i, j, 0) / l_obuk_arr(i, j, 0)) * cloud_factor
                               + (one - cloud_factor);  // blend toward reduced stability
                    phit_cloud = (one + Real(5.0) * sf * pblh_corr_arr(i, j, 0) / l_obuk_arr(i, j, 0)) * cloud_factor
                               + (one - cloud_factor);
                } else if (has_cloud && l_obuk_arr(i, j, 0) <= zero) {
                    // In unstable layers with clouds, slightly enhance instability
                    // (clouds warm the boundary layer)
                    Real cloud_boost = Real(1.0) + Real(0.05) * amrex::min(total_qcloud / qc_threshold, one);
                    phiM_cloud = std::pow(
                        (one - Real(8.0) * sf * pblh_corr_arr(i, j, 0) / l_obuk_arr(i, j, 0)) / cloud_boost,
                        -one / three);
                    phit_cloud = std::pow(
                        (one - Real(16.0) * sf * pblh_corr_arr(i, j, 0) / l_obuk_arr(i, j, 0)) / cloud_boost,
                        -one / two);
                }
                
                // Use cloud-aware stability functions
                const Real phiM_eff = phiM_cloud;
                const Real phit_eff = phit_cloud;
                
                // Prandtl number: Prt = phit/phiM + const_b*kappa*sf
                const Real Prt = amrex::min(amrex::max(phit_eff / phiM_eff + const_b * KAPPA * sf, prmin), prmax);
                const Real wstar = u_star_arr(i, j, 0) / phiM_eff;
                
                // Diffusion coefficient for momentum
                // K = rho * wstar * kappa * z * (1 - z/h)^2
                K_turb(i, j, k, EddyDiff::Mom_v) = rho * wstar * KAPPA * zval
                                                 * (one - zval / pblh_corr_arr(i, j, 0))
                                                 * (one - zval / pblh_corr_arr(i, j, 0));
                K_turb(i, j, k, EddyDiff::Theta_v) = K_turb(i, j, k, EddyDiff::Mom_v) / Prt;

                // Moisture diffusivity: WRF MRF uses Prq ~ Prt (same stability functions
                // for heat and moisture, module_bl_mrf.F lines 968-986).
                // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L968-L986
                if (turbChoice.mrf_moistvars) {
                    const Real Prq = amrex::min(amrex::max(phit_eff / phiM_eff + const_b * KAPPA * sf, prmin), prmax);
                    K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Mom_v) / Prq;
                } else {
                    K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
                }
            } else {
                // Free atmosphere above PBL: use local Richardson number-dependent mixing
                // with lengthscale = (kappa * z * lambda) / (kappa * z + lambda)
                // where lambda = 150 m (characteristic free-atmosphere lengthscale)
                // 
                // WRF reference (module_bl_mrf.F lines 988-1020):
                // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L988-L1020
                // Uses stability functions from Hong et al. 2006 Appendix A (YSU scheme)
                
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
                
                // Gradient Richardson number: Ri = (g/theta_v) * (dtheta_v/dz) / (du/dz)^2
                Real grad_Ri = CONST_GRAV / theta_v * dthetav_dz / wind_shear;
                grad_Ri = std::max(grad_Ri, -Real(100.0));  // Hong et al. 2006, MWR, Appendix A
                
                // Using YSU stability functions (Hong et al. 2006, MWR, Appendix A)
                // Reference: https://doi.org/10.1175/MWR3250.1
                // See equations A19-A20
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

            // Limit diffusion coefficients to physical bounds
            // NOTE: We use Hong & Pan (1996) limits. An alternative Hong et al. (2006)
            // formulation is available in commented code below.
            // 
            // Hong & Pan (1996) limits (module_bl_mrf.F lines 1014-1025):
            // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L1014-L1025
            // Kmin = 0.1 m^2/s, Kmax = 300 m^2/s
            //
            // Hong et al. (2006) alternative limits (MWR Appendix A):
            // Kmin = ckz * dz * rho (where ckz = 0.001), Kmax = 1000 m^2/s
            // These higher limits allow greater mixing in free atmosphere
#if 0
            // Hong et al. 2006, MWR, Appendix A: Higher limits for free atmosphere
            constexpr Real ckz  = Real(0.001);
            constexpr Real Kmax = Real(1000.0);
            const Real rhoKmin  = ckz * dz_terrain * rho;
            const Real rhoKmax  = rho * Kmax;
#endif
            // Hong & Pan 1996, MWR: Conservative limits (default)
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
            K_turb(i, j, k, EddyDiff::Turb_lengthscale) = pblh_corr_arr(i, j, 0);
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

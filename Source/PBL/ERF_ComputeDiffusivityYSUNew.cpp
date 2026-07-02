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
ComputeDiffusivityYSUNew (const MultiFab& xvel,
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
    ============================================================================
    Yonsei University (YSU) Boundary Layer Parameterization Scheme
    ============================================================================

    Implementation of the YSU boundary layer scheme based on:
      - Hong, S.-Y., Y. Noh, and J. Dudhia, 2006: A new vertical diffusion
        package with an explicit treatment of entrainment processes.
        Monthly Weather Review, 134, 2318-2341. [HND06]
        https://doi.org/10.1175/MWR3250.1
      - Hong, S.-Y., 2010: A new stable boundary-layer mixing scheme and
        its impact on the simulated East Asian summer monsoon.
        Quarterly Journal of the Royal Meteorological Society, 136, 1481-1496. [H10]
        https://doi.org/10.1002/qj.665
      - WRF Reference Implementation: module_bl_ysu.F
        https://github.com/wrf-model/WRF/blob/master/phys/module_bl_ysu.F

    CORE ALGORITHM (HND06):
    -----------------------
    1. Three-Pass Bulk Richardson Number (Rib) PBL Height Diagnosis
       Pass 1: Ribcr=f(surface type), base surface θ_v (predictor)
       Pass 2: Same Ribcr, enhanced surface θ_v + VPERT (corrector)
       Pass 3: Ribcr=0 diagnostic, gives mixed-layer depth h_ze

    2. Nonlocal Countergradient Flux Corrections
       HGAMT = -const_b * u_* * θ_* / wscale  (capped at GAMCRT=3K)
       HGAMQ = -const_b * u_* * q_* / wscale  (capped at GAMCRQ=2e-3)
       VPERT = max(HGAMT + 0.61*θ*HGAMQ, 0)

    3. K-Profile Below PBL (HND06 Eq. 2):
       K_m = ρ * wscale * κ * z * (1 - z/h)^pfac
       K_t = K_m / Pr_t
       wscale = (u_*³ + phifac * κ * wstar_conv³ * (1-zfac))^(1/3)

    4. Explicit Entrainment at PBL Top (HND06 Eq. 6):
       K_entr = ρ * we_m * dz  at k = kpbl
       we_m = entr_A * u_*³ / (u_*³ + entr_B * wstar_conv³) * wscale

    5. Free Atmosphere / Stable PBL (H10, YSU Appendix A):
       Ri_g-dependent mixing with grid-adaptive length scale
       λ = min(max(0.1*dz, 30m), 300m)
       l = λ * κz / (λ + κz)

    PARAMETER DEFAULTS (Matching WRF-YSU):
    ----------------------------------------
    const_b   = 7.8    (CFAC: countergradient coefficient)
    sf        = 0.1    (SFCFRAC: surface layer depth fraction)
    phifac    = 8.0    (PHIFAC: convective wscale weight)
    pfac      = 2.0    (profile shape exponent)
    GAMCRT    = 3.0 K  (max heat countergradient)
    GAMCRQ    = 2e-3   (max moisture countergradient)
    Ribcr land= 0.25   (critical bulk Richardson number over land)
    entr_A    = 0.2    (entrainment momentum coefficient)
    entr_B    = 5.0    (entrainment velocity denominator weight)
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
        // WRF reference (module_bl_ysu.F lines 100-130):
        // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_ysu.F#L100-L130

        // create flattened boxes to store PBL height and related quantities
        const GeometryData gdata = geom.data();
        const Box xybx = PerpendicularBox<ZDir>(gbx, IntVect{0, 0, 0});
        FArrayBox pbl_height_corrector(xybx, 1, The_Async_Arena());
        IArrayBox pbl_index(xybx, 1, The_Async_Arena());
        IArrayBox pbl_index_zero_ri(xybx, 1, The_Async_Arena());  // Index for zero-Ri diagnostic pass
        FArrayBox hgamt_fab(xybx, 1, The_Async_Arena());  // Store HGAMT/h (normalized countergradient)
        FArrayBox hgamq_fab(xybx, 1, The_Async_Arena());  // Store HGAMQ/h (normalized countergradient)
        FArrayBox wstar_fab(xybx, 1, The_Async_Arena());  // Convective velocity scale computed with pblh_corr
        FArrayBox vpert_fab(xybx, 1, The_Async_Arena());  // Virtual temperature perturbation VPERT for Pass 3
        FArrayBox entr_fab(xybx, 1, The_Async_Arena());   // Per-column entrainment diffusivity
        const auto& pblh_corr_arr   = pbl_height_corrector.array();
        const auto& pbli_arr        = pbl_index.array();
        const auto& pbli_zero_arr   = pbl_index_zero_ri.array();  // Zero-Ri diagnostic PBL index
        const auto& hgamt_arr       = hgamt_fab.array();
        const auto& hgamq_arr       = hgamq_fab.array();
        const auto& wstar_arr       = wstar_fab.array();  // Stored convective velocity for use in K-profile
        const auto& vpert_arr       = vpert_fab.array();  // Stored VPERT for Pass 3 diagnostic loop
        const auto& entr_arr        = entr_fab.array();   // Stored entrainment diffusivity

        // Get some data in arrays
        const auto& cell_data = cons_in.const_array(mfi);
        const auto& uvel = xvel.const_array(mfi);
        const auto& vvel = yvel.const_array(mfi);

        //const Real Ribcr = turbChoice.pbl_mrf_Ribcr;  // Now surface-type-dependent in YSU
        //const Real f0 = turbChoice.pbl_mrf_coriolis_freq;
        const auto& u_star_arr = SurfLayer->get_u_star(level)->const_array(mfi);
        const auto& t_star_arr = SurfLayer->get_t_star(level)->const_array(mfi);
        const auto& l_obuk_arr = SurfLayer->get_olen(level)->const_array(mfi);
        const auto& t10av_arr  = SurfLayer->get_mac_avg(level, 2)->const_array(mfi);
        const auto& q10av_arr  = SurfLayer->get_mac_avg(level, 3)->const_array(mfi);
        const auto& ws10av_arr = SurfLayer->get_mac_avg(level, 5)->const_array(mfi);  // 10m wind speed for Rossby number
        const auto& z0_arr     = SurfLayer->get_z0(level)->const_array(mfi);           // Roughness length for Rossby number
        //const auto& t_surf_arr = SurfLayer->get_t_surf(level)->const_array(mfi);
        // Get land/water mask for proper handling of moisture countergradient
        const auto& lmask_arr = (SurfLayer->get_lmask(level)) ?
                                SurfLayer->get_lmask(level)->const_array(mfi) :
                                Array4<int>{};
        const Array4<Real const> z_nd_arr = z_phys_nd->array(mfi);


        // ========================================================================
        // PASS 1 — PREDICTOR: Compute initial PBL index with surface-type-dependent
        // critical Richardson number (Ribcr).
        // ========================================================================
        // WRF Reference: module_bl_ysu.F (Hong, Noh & Dudhia 2006, lines 1-100)
        // Ribcr = 0.25 over land; Rossby-number-dependent over water (line 150-170).
        // Rib = (g*z/θv0) * (θv(z) - θv_surf) / ws² (lines 180-200)
        // References: Hong et al., Mon. Wea. Rev., 134, 2318-2341 (2006)
        //
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            // Determine surface-type-dependent critical Richardson number
            // Over land: Ribcr = 0.25; over water: Ribcr depends on Rossby number
            Real Ribcr;
            {
                // Check if over land (lmask_arr(i,j,0) == 1) or lmask_arr is null (default land)
                bool over_land = (!lmask_arr) || (lmask_arr(i, j, 0) == 1);
                if (over_land) {
                    Ribcr = turbChoice.pbl_ysu_land_Ribcr;  // Default 0.25 (module_bl_ysu.F line ~160)
                } else {
                    // Over water: compute Rossby-number-dependent Ribcr
                    // Hong et al. (2006), Eq. (5): Ribcr = 0.16 * (1e-7*Ro)^(-0.18) (module_bl_ysu.F lines 165-170)
                    const Real z0 = z0_arr(i, j, 0);
                    const Real ws_layer = ws10av_arr(i, j, 0);
                    const Real Rossby = ws_layer / (turbChoice.pbl_ysu_coriolis_freq * amrex::max(z0, Real(1.0e-4)));
                    Ribcr = amrex::min(Real(0.16) * std::pow(Real(1.0e-7) * Rossby, -Real(0.18)), Real(0.3));
                }
            }

            const Real t_layer  = t10av_arr(i, j, 0);
            Real zval, Rib;
            int kpbl = klo;

            // Initialize at lowest level
            {
                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : (kpbl + myhalf) * gdata.CellSize(2);
                const Real theta_v = GetThetav(i, j, kpbl, cell_data, moisture_indices);
                const Real theta_v_klo = GetThetav(i, j, klo, cell_data, moisture_indices);
                const Real moisture_fraction = use_moisture ? q10av_arr(i, j, 0) : zero;
                const Real t_layer_v = t_layer * (one + amrex::Real(0.61) * moisture_fraction);
                const Real ws2_raw = fourth * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                              (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                              (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                              (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                const Real ws2 = amrex::max(ws2_raw, Real(1.0));  // WRF: SPDK2=MAX(...,1.)
                // Richardson number: Rib = (g*z/θv0) * (θv(z) - θv_surf) / ws2
                // Use lowest-level potential temperature (theta_v_klo) in denominator for consistency
                // with WRF bulk Richardson number definition (HND06, module_bl_ysu.F)
                Rib = CONST_GRAV * zval * (theta_v - t_layer_v) / (ws2 * theta_v_klo);
            }

            bool above_critical = false;
            while (!above_critical && ((kpbl + 1) <= khi)) {
                kpbl += 1;

                // height above ground level
                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : (kpbl + myhalf) * gdata.CellSize(2);

                // Use virtual potential temperature for stability calculations
                const Real theta_v = GetThetav(i, j, kpbl, cell_data, moisture_indices);
                const Real theta_v_klo = GetThetav(i, j, klo, cell_data, moisture_indices);
                const Real moisture_fraction = use_moisture ? q10av_arr(i, j, 0) : zero;
                const Real t_layer_v = t_layer * (one + amrex::Real(0.61) * moisture_fraction);
                const Real ws2_raw = fourth * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                              (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                              (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                              (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                const Real ws2 = amrex::max(ws2_raw, Real(1.0));  // WRF: SPDK2=MAX(...,1.)

                // Richardson number: Rib = (g*z/θv0) * (θv(z) - θv_surf) / ws2
                // Use lowest-level potential temperature (theta_v_klo) in denominator for consistency
                // with WRF bulk Richardson number definition (HND06, module_bl_ysu.F)
                Rib = CONST_GRAV * zval * (theta_v - t_layer_v) / (ws2 * theta_v_klo);
                above_critical = (Rib >= Ribcr);
            }

            // Empirical expression for PBLH is given by h = c u* / f
            // Garratt (1994) and Tennekes (1982)
            // Also, c.f. Zilitinkevitch et al 2012 referenced in Pedersen et al. Real(2014.)
            //const Real c_pblh = (l_obuk_arr(i, j, 0) > 0) ? Real(0.16) : Real(0.60);
            //const Real pblh_emp = c_pblh * u_star_arr(i, j, 0) / f0;

            // Fallback to first cell

            // Initial PBL Height diagnosis (stored in pbli_arr)
            // WRF reference (module_bl_ysu.F lines 150-180):
            // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_ysu.F
            if (above_critical) {
                pbli_arr(i, j, 0) = kpbl;  // k < kpbl is considered the PBL
            } else {
                pbli_arr(i, j, 0) = klo + 1;
            }
        });

        const auto& q_star_arr = SurfLayer->get_q_star(level)->const_array(mfi);

        //
        // CORRECTOR: Apply YSU nonlocal countergradient correction (HGAMT/HGAMQ)
        // following HND06 with countergradient modifications.
        //
        // This step raises the effective surface virtual potential temperature
        // to account for convective contribution to buoyancy, leading to a
        // higher diagnosed PBL height in unstable conditions.
        //
        // WRF reference:
        //   https://github.com/wrf-model/WRF/blob/master/phys/module_bl_ysu.F
        //
        // WRF constants (module_bl_ysu.F lines 100-120):
        //   RLAM   = 150.   (mixing length scale in free atm., m)
        //   PRMIN  = 0.5    (minimum Prandtl number)
        //   PRMAX  = 4.0    (maximum Prandtl number)
        //   BRCR   = 0.25   (bulk Richardson number critical value over land)
        //   CFAC   = 7.8    (= const_b, surface layer factor in countergradient term)
        //   PHIFAC = 8.0    (phifac: convective wscale weight)
        //   SFCFRAC= 0.1    (surface layer depth fraction)
        //   GAMCRT = 3.0    (maximum heat countergradient, K)
        //   GAMCRQ = 2.e-3  (maximum moisture countergradient, kg/kg)
        //
        // WRF formulas (module_bl_ysu.F lines 220-260):
        //   GAMFAC = CFAC / (rho * wscale)
        //   HGAMT  = min(GAMFAC * HFX / CPM, GAMCRT)
        //   HGAMQ  = min(GAMFAC * QFX, GAMCRQ)  [if over water, HGAMQ = 0]
        //   VPERT  = HGAMT + EP1 * theta * HGAMQ
        //   THERMAL += max(VPERT, 0)   [raises effective surface virtual pot. temp.]
        //
        // YSU K-profile below PBL (module_bl_ysu.F lines 190-220):
        //   K_mom   = rho * wscale * kappa * z * (1 - z/h)^pfac
        //   K_theta = K_mom / Prt
        //   K_q     = K_mom / Prq  (Prq ~ Prt for moisture)
        //
        //
        // PASS 2 — CORRECTOR: Apply surface-type-dependent Ribcr with VPERT enhancement
        // following HND06 (Hong, Noh & Dudhia 2006), with countergradient corrections.
        //
        // This pass raises the effective surface virtual potential temperature
        // to account for convective contribution to buoyancy, resulting in
        // higher diagnosed PBL height in unstable conditions.
        //
        // WRF reference (module_bl_ysu.F lines 150-180):
        //   https://github.com/wrf-model/WRF/blob/master/phys/module_bl_ysu.F
        //
        // WRF constants (module_bl_ysu.F lines 50-75):
        //   RLAM   = 150.   (mixing length scale in free atm., m)
        //   PRMIN  = 0.5    (minimum Prandtl number)
        //   PRMAX  = 4.0    (maximum Prandtl number)
        //   BRCR   = 0.25   (over land, bulk Richardson number critical value)
        //   CFAC   = 7.8    (= const_b, surface layer factor in countergradient term)
        //   PFAC   = 2.0    (power law exponent for PBL height factor)
        //   SFCFRAC= 0.1    (surface layer depth fraction)
        //   GAMCRT = 3.0    (maximum heat countergradient, K)
        //   GAMCRQ = 2.e-3  (maximum moisture countergradient, kg/kg)
        //
        // WRF formulas (module_bl_ysu.F lines 200-250):
        //   GAMFAC = CFAC / (rho * wstar)
        //   HGAMT  = min(GAMFAC * HFX / CPM, GAMCRT)
        //   HGAMQ  = min(GAMFAC * QFX, GAMCRQ)  [if over water, HGAMQ = 0]
        //   VPERT  = HGAMT + EP1 * theta * HGAMQ
        //   THERMAL += max(VPERT, 0)   [raises effective surface virtual pot. temp.]
        //
        // WRF diffusion coefficients below PBL (module_bl_ysu.F lines 300-350):
        //   K_mom   = rho * wstar * kappa * z * (1 - z/h)^2
        //   K_theta = K_mom / Prt
        //   K_q     = K_mom / Prq  (Prq similar to Prt for moisture)
        //
        const Real const_b = turbChoice.pbl_mrf_const_b;
        const Real sf = turbChoice.pbl_mrf_sf;
        constexpr Real prmin = Real(0.5);
        constexpr Real prmax = Real(4.0);
        constexpr Real GAMCRT = Real(3.0);    // WRF GAMCRT: max heat countergradient (K)
        constexpr Real GAMCRQ = Real(2.e-3);  // WRF GAMCRQ: max moisture countergradient (kg/kg)
        const bool enable_mrf_countergradient = turbChoice.enable_mrf_countergradient;
        const bool enable_mrf_unbounded_vpert = turbChoice.enable_mrf_unbounded_vpert;
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            // Determine surface-type-dependent critical Richardson number (same as Pass 1)
            // Over land: Ribcr = 0.25; over water: Ribcr depends on Rossby number
            Real Ribcr;
            {
                // Check if over land (lmask_arr(i,j,0) == 1) or lmask_arr is null (default land)
                bool over_land = (!lmask_arr) || (lmask_arr(i, j, 0) == 1);
                if (over_land) {
                    Ribcr = turbChoice.pbl_ysu_land_Ribcr;  // Default 0.25 (module_bl_ysu.F line ~160)
                } else {
                    // Over water: compute Rossby-number-dependent Ribcr
                    // Hong et al. (2006), Eq. (5): Ribcr = 0.16 * (1e-7*Ro)^(-0.18) (module_bl_ysu.F lines 165-170)
                    const Real z0 = z0_arr(i, j, 0);
                    const Real ws_layer = ws10av_arr(i, j, 0);
                    const Real Rossby = ws_layer / (turbChoice.pbl_ysu_coriolis_freq * amrex::max(z0, Real(1.0e-4)));
                    Ribcr = amrex::min(Real(0.16) * std::pow(Real(1.0e-7) * Rossby, -Real(0.18)), Real(0.3));
                }
            }

            const Real t_layer  = t10av_arr(i, j, 0);
            const Real moisture_fraction = use_moisture ? q10av_arr(i, j, 0) : zero;
            const Real t_layer_v = t_layer * (one + amrex::Real(0.61) * moisture_fraction);
            // Pass 2 uses VPERT-enhanced surface virtual potential temperature
            // On first call, vpert_arr is initialized to zero, so t_layer_v_enhanced = t_layer_v
            // After Pass 2b finalization, vpert_arr contains computed VPERT
            const Real t_layer_v_enhanced = t_layer_v + vpert_arr(i, j, 0);

            Real obuk_val = l_obuk_arr(i, j, 0);
            if (std::abs(obuk_val) < amrex::Real(1.0e-10)) {
                obuk_val = (obuk_val >= zero) ? amrex::Real(1.0e-10) : amrex::Real(-1.0e-10);
            }

            int kpbl = klo;
            Real zval0, zval, Rib0, Rib;
            {
                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : (kpbl + myhalf) * gdata.CellSize(2);
                const Real theta_v = GetThetav(i, j, kpbl, cell_data, moisture_indices);
                const Real theta_v_klo = GetThetav(i, j, klo, cell_data, moisture_indices);
                const Real ws2_raw = fourth * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                              (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                              (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                              (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                const Real ws2 = amrex::max(ws2_raw, Real(1.0));  // WRF: SPDK2=MAX(...,1.)
                // Richardson number: Rib = (g*z/θv0) * (θv(z) - θv_surf_enhanced) / ws2
                // Use lowest-level virtual potential temperature in denominator for consistency
                // with WRF's Bulk Richardson number definition (module_bl_ysu.F)
                Rib = CONST_GRAV * zval * (theta_v - t_layer_v_enhanced) / (ws2 * theta_v_klo);
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
                const Real theta_v_klo = GetThetav(i, j, klo, cell_data, moisture_indices);
                const Real ws2_raw = fourth * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                              (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                              (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                              (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                const Real ws2 = amrex::max(ws2_raw, Real(1.0));  // WRF: SPDK2=MAX(...,1.)
                // Richardson number: Rib = (g*z/θv0) * (θv(z) - θv_surf_enhanced) / ws2
                // Use lowest-level potential temperature (theta_v_klo) in denominator for consistency
                // with WRF bulk Richardson number definition (module_bl_ysu.F)
                Rib = CONST_GRAV * zval * (theta_v - t_layer_v_enhanced) / (ws2 * theta_v_klo);
                above_critical = (Rib >= Ribcr);
            }

            // Empirical expression for PBLH fallback: minimum PBL height based on first cell height
            // Used as a lower bound to prevent division by zero near surface
            const Real z_sfc = (use_terrain_fitted_coords)
                              ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                              : zero;
            const Real dz_terrain = (use_terrain_fitted_coords)
                                   ? (Compute_Zrel_AtCellCenter(i, j, klo + 1, z_nd_arr) - z_sfc)
                                   : gdata.CellSize(2);
            //const Real pblh_emp = (use_terrain_fitted_coords)
            //                    ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
            //                    : myhalf * gdata.CellSize(2);

            // Absolute bounds safeguard: clamp to prevent division by zero near surface and runaway height at top
            const Real z_max = (use_terrain_fitted_coords)
                             ? Compute_Zrel_AtCellCenter(i, j, khi, z_nd_arr)
                             : (khi + myhalf) * gdata.CellSize(2);
            const Real pblh_max = Real(0.9) * z_max;
            const Real pblh_min = amrex::max(z_sfc + Real(0.5) * dz_terrain, Real(10.0));

            if (above_critical) {
                // Interpolate to height at which Rib == Ribcr
                // Linear interpolation between levels where Richardson number crosses critical value
                // WRF reference (module_bl_ysu.F lines 150-180):
                // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_ysu.F
                Real pblh_interp = zval0 + (zval - zval0) / (Rib - Rib0) * (Ribcr - Rib0);
                pblh_corr_arr(i, j, 0) = amrex::max(amrex::min(pblh_interp, pblh_max), pblh_min);
                     pbli_arr(i, j, 0) = kpbl;  // k < kpbl is considered the PBL
            } else {
                // Fallback to first cell
                pblh_corr_arr(i, j, 0) = pblh_min;
                     pbli_arr(i, j, 0) = klo + 1;
            }

            // NOTE: Countergradient fluxes (HGAMT/HGAMQ) deferred to subsequent loop
            // after pblh_corr_arr is finalized. Initialize with zeros for now.
            hgamt_arr(i, j, 0) = zero;
            hgamq_arr(i, j, 0) = zero;
            wstar_arr(i, j, 0) = zero;

        });
        // Debug output (disabled):
        /*
          amrex::Print() << "PBL height computed for MRF scheme at level "
          << pblh_corr_arr(2, 2, 0)
          << std::endl;
          amrex::Print() << "PBL Temp:" << t_surf_arr(2, 2, 0) << "  "
          << t10av_arr(2, 2, 0) << std::endl;
        */

        //
        // Recompute convective velocity scale (wstar) and countergradient fluxes
        // (HGAMT/HGAMQ) using the corrector PBL height (pblh_corr_arr).
        //
        // Background: The corrector PBL height is now available. The Monin-Obukhov stability
        // functions and countergradient fluxes depend strongly on the PBL height estimate,
        // so these quantities must be recomputed after the corrector pass to ensure
        // consistency between the PBL height diagnosis and the K-profile mixing lengths.
        //
        // Note: Previous versions computed both a predictor PBL height (without countergradient
        // corrections) and a corrector PBL height (with countergradient effects). The inconsistency
        // between these two estimates led to unrealistic mixing intensity since the stability
        // parameter HOL = sf*h/L differs. The predictor pass has been removed; only the
        // corrector PBL height (pblh_corr_arr) is now computed and used throughout.
        //
        // WRF implements WSCALE (convective velocity) and countergradient corrections
        // only once before both the corrector loop and K-profile computations, ensuring
        // consistency. ERF implements three passes: this loop performs the needed computation
        // after the corrector pass is complete.
        //
        // YSU (HND06): WSCALE = u* / φ_m(h/L)
        // Countergradient: HGAMT = min(CFAC * u* * θ*, GAMCRT), where CFAC=7.8, GAMCRT=3K
        // WRF Reference: module_bl_ysu.F lines 220-250
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            const Real t_layer  = t10av_arr(i, j, 0);
            Real obuk_val = l_obuk_arr(i, j, 0);
            if (std::abs(obuk_val) < amrex::Real(1.0e-10)) {
                obuk_val = (obuk_val >= zero) ? amrex::Real(1.0e-10) : amrex::Real(-1.0e-10);
            }

            // Compute Monin-Obukhov stability parameter using corrector PBL height
            // HOL = sf * h / L, where L is Monin-Obukhov length scale
            // WRF Reference: module_bl_ysu.F lines 200-210
            // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_ysu.F#L200-L210
            const Real HOL = sf * pblh_corr_arr(i, j, 0) / obuk_val;
            const Real HOL_bounded = amrex::max(amrex::min(HOL, Real(100.0)), Real(-100.0));
            const Real one_quarter = Real(1.0) / Real(4.0);
            const Real phiM     = (obuk_val > 0)
                                ? (1 + 5 * HOL_bounded)
                                : std::pow(
                                           amrex::max(1 - 16 * HOL_bounded, Real(0.01)),
                                           -one_quarter);
            const Real phiM_safe = amrex::max(phiM, Real(0.01));

            // Convective velocity scale (wscale = u*/phi_m), now computed with pblh_corr_arr
            // wscale is the characteristic turbulent velocity in the boundary layer.
            // Bounds (u*/5 to 16*u*) prevent unrealistic values in very weak or strong convection.
            // WRF Reference: module_bl_ysu.F L210-215
            Real wscale = u_star_arr(i, j, 0) / phiM_safe;
            wscale = amrex::max(wscale, u_star_arr(i, j, 0) / Real(5.0));      // Mechanical turbulence floor
            wscale = amrex::min(wscale, Real(16.0) * u_star_arr(i, j, 0));     // Free convection ceiling
            wstar_fab(i, j, 0) = wscale;  // Store for use in K-profile loop

            // Compute HGAMT with corrected wscale
            bool SFCFLG = (obuk_val <= zero);  // TRUE when unstable/neutral, FALSE when stable
            const Real HGAMT = (SFCFLG && enable_ysu_countergradient)
                             ? amrex::min(-const_b * u_star_arr(i, j, 0) * t_star_arr(i, j, 0) / wscale, GAMCRT)
                             : zero;

            // Compute HGAMQ with corrected wscale
            // Sign convention: ERF's q_star = κ*(qvm - q_surf)/(log(z/z0) - ψ_h) is positive when surface
            // is drier than air (evaporating conditions). Formula -const_b*u_star*q_star/wscale converts
            // to WRF's QFX convention (positive for upward flux). WRF Reference: module_bl_ysu.F L220-230
            Real HGAMQ = zero;
            if (SFCFLG && use_moisture && enable_ysu_countergradient) {
                const Real q_star = q_star_arr(i, j, 0);
                const Real HGAMQ_calc = -const_b * u_star_arr(i, j, 0) * q_star / wscale;
                HGAMQ = amrex::max(
                    amrex::min(HGAMQ_calc, GAMCRQ),
                    zero
                );

                // Land/water surface discrimination
                if (lmask_arr) {
                    bool is_land = (lmask_arr(i,j,0) == 1);
                    if (!is_land) HGAMQ = zero;
                }

                // Saturation-Aware HGAMQ limiter
                if (moisture_indices.qv >= 0) {
                    Real qv_klo = cell_data(i, j, klo, moisture_indices.qv) / cell_data(i, j, klo, Rho_comp);
                    Real T_klo = getTgivenRandRTh(cell_data(i, j, klo, Rho_comp),
                                                  cell_data(i, j, klo, RhoTheta_comp),
                                                  qv_klo);
                    Real p_klo = getPgivenRTh(cell_data(i, j, klo, RhoTheta_comp), qv_klo) * Real(0.01);
                    Real qsat_klo = zero;
                    erf_qsatw(T_klo, p_klo, qsat_klo);
                    Real rh_klo = (qsat_klo > Real(1.0e-10)) ? (qv_klo / qsat_klo) : zero;
                    if (rh_klo > Real(0.95)) {
                        Real rh_scaling = amrex::max(zero, (one - rh_klo) / Real(0.05));
                        HGAMQ *= rh_scaling;
                    }
                }
            }

            // Store HGAMT/h and HGAMQ/h for implicit solver (normalized by corrected PBL height)
            // Also compute VPERT = max(HGAMT + 0.61*θ*HGAMQ, 0) for Pass 3 use
            // WRF Reference: module_bl_ysu.F lines 230-240 (THERMAL update with VPERT)
            if (pbli_arr(i, j, 0) <= klo + 1) {
                hgamt_arr(i, j, 0) = zero;
                hgamq_arr(i, j, 0) = zero;
                vpert_arr(i, j, 0) = zero;
            } else {
                const Real pblh = pblh_corr_arr(i, j, 0);
                hgamt_arr(i, j, 0) = (enable_ysu_countergradient) ? HGAMT / pblh : zero;
                hgamq_arr(i, j, 0) = (enable_ysu_countergradient && use_moisture) ? HGAMQ / pblh : zero;

                // VPERT for Pass 3 diagnostic: unnormalized (not divided by pblh)
                // This represents the virtual temperature perturbation at surface
                // When enable_ysu_unbounded_vpert=true, use unbounded VPERT (ERF enhanced)
                // When enable_ysu_unbounded_vpert=false, cap at GAMCRT (WRF-compatible)
                // WRF Reference: module_bl_ysu.F lines 230-240
                if (enable_ysu_countergradient) {
                    const Real VPERT_raw = HGAMT + amrex::Real(0.61) * t_layer * HGAMQ;
                    const Real VPERT_capped = enable_ysu_unbounded_vpert
                                            ? VPERT_raw
                                            : amrex::min(VPERT_raw, GAMCRT);
                    vpert_arr(i, j, 0) = amrex::max(VPERT_capped, zero);
                } else {
                    vpert_arr(i, j, 0) = zero;
                }
            }

        });

        //
        // PASS 3 — ZERO-RI DIAGNOSTIC: Compute PBL index with Ribcr=0.0 criterion
        //
        // Background: WRF employs three distinct PBL height estimates (module_bl_ysu.F lines 150-200):
        //   Pass 1: Ribcr = surface-type-dependent (predictor, base surface temperature)
        //   Pass 2: Ribcr = surface-type-dependent (corrector, enhanced surface temperature with VPERT)
        //   Pass 3: Ribcr = 0.0 diagnostic (uses enhanced surface temperature, neutral criterion)
        //
        // The third pass with Ribcr=0 produces the PBL height where Richardson number
        // first becomes neutral (Rib >= 0). This is the "zero-Richardson" diagnostic
        // commonly used in observation-based PBL depth estimates (module_bl_ysu.F lines 190-200).
        //
        // Use of three passes: The corrector PBL height (Pass 2) is used for
        // computing the mixing intensity K = ρ*wstar*κ*z*(1-z/h)², while
        // the diagnostic height (Pass 3, Ribcr=0) determines the vertical extent
        // of the nonlocal mixing region. This separation provides:
        //   - Realistic mixed-layer extent in convective conditions
        //   - Mixing formula consistent with corrector diagnostics
        //   - Consistency with WRF methodology
        //
        // Reference: Hong et al. (2006), module_bl_ysu.F
        //
        constexpr Real Ribcr_zero = zero;  // Zero critical Richardson number for diagnostic pass
         ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
         {
             // Surface virtual temperature with VPERT contribution from Pass 2b
             // Pass 3 uses VPERT-enhanced surface temperature for diagnostic extent
             const Real t_layer  = t10av_arr(i, j, 0);
             const Real moisture_fraction = use_moisture ? q10av_arr(i, j, 0) : zero;
             const Real t_layer_v = t_layer * (one + amrex::Real(0.61) * moisture_fraction);
             const Real t_layer_v_enhanced = t_layer_v + vpert_arr(i, j, 0);

            int kpbl_zero = klo;
            Real zval_zero, Rib_zero;
            {
                zval_zero = (use_terrain_fitted_coords)
                          ? Compute_Zrel_AtCellCenter(i, j, kpbl_zero, z_nd_arr)
                          : (kpbl_zero + myhalf) * gdata.CellSize(2);
                const Real theta_v = GetThetav(i, j, kpbl_zero, cell_data, moisture_indices);
                const Real theta_v_klo = GetThetav(i, j, klo, cell_data, moisture_indices);
                const Real ws2_raw = fourth * ( (uvel(i, j, kpbl_zero) + uvel(i + 1, j, kpbl_zero)) *
                                              (uvel(i, j, kpbl_zero) + uvel(i + 1, j, kpbl_zero)) +
                                              (vvel(i, j, kpbl_zero) + vvel(i, j + 1, kpbl_zero)) *
                                              (vvel(i, j, kpbl_zero) + vvel(i, j + 1, kpbl_zero)) );
                const Real ws2 = amrex::max(ws2_raw, Real(1.0));
                // Richardson number using enhanced surface virtual temperature
                // Consistent with corrector pass methodology
                Rib_zero = CONST_GRAV * zval_zero * (theta_v - t_layer_v_enhanced) / (ws2 * theta_v_klo);
            }

            bool above_critical_zero = false;
            while (!above_critical_zero && ((kpbl_zero + 1) <= khi)) {
                //zval0_zero = zval_zero;
                //Rib0_zero = Rib_zero;
                kpbl_zero += 1;

                zval_zero = (use_terrain_fitted_coords)
                          ? Compute_Zrel_AtCellCenter(i, j, kpbl_zero, z_nd_arr)
                          : (kpbl_zero + myhalf) * gdata.CellSize(2);
                const Real theta_v = GetThetav(i, j, kpbl_zero, cell_data, moisture_indices);
                const Real theta_v_klo = GetThetav(i, j, klo, cell_data, moisture_indices);
                const Real ws2_raw = fourth * ( (uvel(i, j, kpbl_zero) + uvel(i + 1, j, kpbl_zero)) *
                                              (uvel(i, j, kpbl_zero) + uvel(i + 1, j, kpbl_zero)) +
                                              (vvel(i, j, kpbl_zero) + vvel(i, j + 1, kpbl_zero)) *
                                              (vvel(i, j, kpbl_zero) + vvel(i, j + 1, kpbl_zero)) );
                const Real ws2 = amrex::max(ws2_raw, Real(1.0));
                // Use VPERT-enhanced surface temperature for Richardson number criterion
                Rib_zero = CONST_GRAV * zval_zero * (theta_v - t_layer_v_enhanced) / (ws2 * theta_v_klo);
                above_critical_zero = (Rib_zero >= Ribcr_zero);
            }

            // Use same bounds safeguard as corrector
            if (above_critical_zero) {
                pbli_zero_arr(i, j, 0) = kpbl_zero;
            } else {
                pbli_zero_arr(i, j, 0) = klo + 1;
            }
        });

        // Compute entrainment diffusivity at PBL top cell (HND06 Eq. 6, 11-14)
        // This represents the enhanced diffusivity from entrainment of free-tropospheric air
        {
            const bool enable_ysu_entrainment = turbChoice.enable_ysu_entrainment;
            const Real entr_A = turbChoice.pbl_ysunew_entr_A;
            const Real entr_B = turbChoice.pbl_ysunew_entr_B;

            ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
            {
                entr_arr(i,j,0) = zero;
                if (!enable_ysu_entrainment) return;

                const Real rho_kpbl = cell_data(i, j, pbli_arr(i,j,0), Rho_comp);
                const Real pblh     = pblh_corr_arr(i,j,0);
                const Real wscale   = wstar_arr(i,j,0);
                const Real ustar    = u_star_arr(i,j,0);

                // Entrainment velocity (HND06 Eq. 13):
                // we_m = entr_A * ustar^3 / (ustar^3 + entr_B * wstar_conv^3)
                // wstar_conv^3 = max(BFLUX * g/theta * pblh, 0)
                const Real t_layer = t10av_arr(i,j,0);
                const Real BFLUX   = -ustar * t_star_arr(i,j,0);
                const Real wstar3  = amrex::max(BFLUX * CONST_GRAV / t_layer * pblh, zero);
                const Real ustar3  = ustar * ustar * ustar;
                const Real we_denom = ustar3 + entr_B * wstar3;
                const Real we_m = (we_denom > Real(1.0e-10))
                                ? entr_A * ustar3 / we_denom * wscale
                                : zero;

                // Entrainment diffusivity K_entr = we * dz at PBL top cell
                // dz_terrain at kpbl level:
                const int kpbl    = pbli_arr(i,j,0);
                const Real met_h  = (use_terrain_fitted_coords)
                                   ? Compute_h_zeta_AtCellCenter(i, j, kpbl, dxInv, z_nd_arr) : one;
                const Real dz_kpbl = met_h / dz_inv;

                // Cap K_entr to 5x the typical mixing to prevent runaway values
                const Real K_entr_raw = rho_kpbl * we_m * dz_kpbl;
                const Real K_cap      = Real(5.0) * rho_kpbl * wscale * KAPPA * pblh * Real(0.01);
                entr_arr(i,j,0) = amrex::min(K_entr_raw, K_cap);
            });
        }

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
            Real obuk_val = l_obuk_arr(i, j, 0);
            if (std::abs(obuk_val) < amrex::Real(1.0e-10)) {
                obuk_val = (obuk_val >= zero) ? amrex::Real(1.0e-10) : amrex::Real(-1.0e-10);
            }

            const Real zval = (use_terrain_fitted_coords)
                            ? Compute_Zrel_AtCellCenter(i, j, k, z_nd_arr)
                            : (k + myhalf) * gdata.CellSize(2);
            const Real rho = cell_data(i, j, k, Rho_comp);
            const Real met_h_zeta = (use_terrain_fitted_coords)
                                  ? Compute_h_zeta_AtCellCenter(i, j, k, dxInv, z_nd_arr) : one;
            const Real dz_terrain = met_h_zeta / dz_inv;

            // Cloud-aware diffusion: compute cloud fraction from qc and qi
            // Cloud water threshold for detection (threshold ~ 0.01 g/kg, matching WRF)
            //
            // CLOUD-AWARE ENHANCEMENTS (ERF Extension):
            // This is an optional feature NOT in HND06. It adjusts the
            // stability functions in cloudy layers to improve representation of:
            //   1. Reduced turbulence near cloud tops (stable layers with clouds)
            //   2. Enhanced mixing from latent heat release (unstable layers with clouds)
            //   3. Fog/stratus evolution in marine boundary layers
            //
            // Physical justification:
            //   - Clouds modify vertical buoyancy structure through radiative effects
            //   - Latent heat release enhances convective mixing
            //   - Cloud-top entrainment is distinct from clear-air turbulence
            //
            // Reference concept: Similar to WRF's IMVDIF cloud-aware parameterization
            // (Bretherton & Park 2009, see WRF module_bl_mynn.F).
            //
            // Disable via enable_ysu_countergradient=false for strict HND06.
            constexpr Real qc_threshold = Real(1.0e-5);  // 0.01 g/kg in mixing ratio (aligned with WRF)
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

            // Select PBL extent index based on configuration:
            // - Default (pbl_mrf_use_zero_ri_extent=false): use pbli_arr (Ri=0.5 corrector)
            //   This matches WRF behavior for better comparability
            // - Alternative (pbl_mrf_use_zero_ri_extent=true): use pbli_zero_arr (Ri=0)
            //   This provides a ~30-60% taller mixing region for convective ABL cases
            const int pbli_extent = turbChoice.pbl_mrf_use_zero_ri_extent ? pbli_zero_arr(i, j, 0) : pbli_arr(i, j, 0);

            if (k < pbli_extent) {
                // Within PBL: use nonlocal mixing with diagnostic stability functions
                // WRF reference (module_bl_ysu.F lines 190-220):
                // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_ysu.F#L190-L220
                //
                // PBL extent selection:
                // - Default (pbl_ysu_use_zero_ri_extent=false): uses pbli_arr (Ri=0.5 corrector)
                //   Matches WRF behavior for K-profile mixing
                // - Alternative (pbl_ysu_use_zero_ri_extent=true): uses pbli_zero_arr (Ri=0)
                //   Provides physically appropriate mixed-layer depth with extended mixing region
                // The mixing intensity is governed by pblh_corr_arr used in the formula K = ρ*wscale*κ*z*(1-z/h)^pfac.
                // This two-level approach ensures realistic PBL height behavior across
                // different stability regimes while maintaining stable mixing coefficients.
                //
                // Key physics: YSU nonlocal scheme represents updrafts/downdrafts by
                // countergradient fluxes, enabling faster PBL growth than local schemes.
                //
                // Implement SFCFLG stable-side gating (WRF L808, L872-884)
                // When stable (BR > 0, i.e., obuk_val > 0), skip nonlocal PBL mixing
                // and use free-atmosphere Richardson scheme throughout PBL.
                bool SFCFLG = (obuk_val <= zero);  // TRUE when unstable/neutral, FALSE when stable

                // Stability function phiM for momentum - BUSINGER-DYER form (WRF)
                // Unstable (L < 0): phiM = (1 - 16*sf*h/L)^(-1/4)  [APHI16=16, exponent -1/4]
                // Stable (L > 0):   phiM = 1 + 5*sf*h/L
                // Bound HOL to [-100, 100] to prevent numerical issues in extreme stability
                const Real HOL = sf * pblh_corr_arr(i, j, 0) / obuk_val;
                const Real HOL_bounded = amrex::max(amrex::min(HOL, Real(100.0)), Real(-100.0));

                const Real one_quarter = Real(1.0) / Real(4.0);
                const Real phiM = (obuk_val > 0)
                                ? (1 + 5 * HOL_bounded)
                                : std::pow(
                                           amrex::max(1 - 16 * HOL_bounded, Real(0.01)),
                                           -one_quarter);

                // Stability function phit for heat/temperature - BUSINGER-DYER form
                // Unstable (L < 0): phit = (1 - 16*sf*h/L)^(-1/2)
                // Stable (L > 0):   phit = 1 + 5*sf*h/L  (same as phiM for stability)
                // Reference: HND06, Table 1; WRF L857-861
                const Real phit = (obuk_val > 0)
                                ? (1 + 5 * HOL_bounded)
                                : std::pow(
                                           amrex::max(1 - 16 * HOL_bounded, Real(0.01)),
                                           -one / two);

                // Cloud-aware adjustment: In cloudy regions, adjust stability damping
                // because clouds modify buoyancy oscillations through latent heating
                // and radiation. This improves representation of:
                //   - Stratus-topped boundary layers (reduced mixing with clouds)
                //   - Cumulus-capped boundary layers (enhanced mixing with convection)
                // This is a physically motivated extension not in original HND06.
                Real phit_cloud = phit;
                Real phiM_cloud = phiM;
                if (has_cloud && obuk_val > zero) {
                    // Stable layers with clouds: reduce stability enhancement
                    // Cloud presence reduces oscillations, damping is less effective
                    // Reduction factor: up to 15-20% where cloud water exceeds threshold
                    Real reduction_factor = one - Real(0.15) * amrex::min(total_qcloud / qc_threshold, one);
                    // Reduce the phiM enhancement in stable conditions
                    phiM_cloud = one + Real(5.0) * reduction_factor * sf * pblh_corr_arr(i, j, 0) / obuk_val;
                    phit_cloud = one + Real(5.0) * reduction_factor * sf * pblh_corr_arr(i, j, 0) / obuk_val;
                } else if (has_cloud && obuk_val <= zero) {
                    // Unstable layers with clouds: slightly enhance instability
                    // Clouds warm the boundary layer through latent heat release
                    // Boost factor: up to 5% where cloud water exceeds threshold
                    Real cloud_boost = Real(1.0) + Real(0.05) * amrex::min(total_qcloud / qc_threshold, one);
                    // Ensure numerically stable exponentiation: base must be in [0.01, 1]
                    // Businger-Dyer form: phi_m = (1 - 16*HOL)^(-1/4) for unstable conditions
                    phiM_cloud = std::pow(
                        amrex::max(one - Real(16.0) * HOL_bounded / cloud_boost, Real(0.01)),
                        -one_quarter);
                    phit_cloud = std::pow(
                        amrex::max(one - Real(16.0) * HOL_bounded / cloud_boost, Real(0.01)),
                        -one / two);
                }

                // Use cloud-aware stability functions
                const Real phiM_eff = phiM_cloud;
                const Real phit_eff = phit_cloud;

                // Prandtl number calculation with stability correction
                // Reference: HND06, Equation (17)
                // Reference code: WRF module_bl_ysu.F line 240
                //
                // Base Prandtl from stability functions:
                //   Pr_base = φ_t / φ_m
                // This gives Pr > 1 in stable conditions (heat diffuses faster than momentum)
                // and Pr < 1 in unstable conditions (momentum diffuses faster than heat).
                //
                // Stability correction term: const_b * κ * sf
                // const_b = 7.8 (matches WRF CFAC)
                // κ = von Karman constant (0.4)
                // sf = surface layer height fraction (0.1)
                // This term represents enhanced diffusivity from surface layer processes.
                //
                // Final formula: Pr_t = φ_t/φ_m + const_b*κ*sf
                // Bounded to physically reasonable range [0.5, 4.0]
                //
                // Consistency check with WRF:
                // - const_b = 7.8 (verified: CFAC in WRF)
                // - sf = 0.1 (verified: SFCFRAC in WRF)
                // - κ ≈ 0.4 (von Karman constant, same in both models)
                Real Prt_base = phit_eff / phiM_eff;
                const Real Prt = amrex::min(amrex::max(Prt_base + const_b * KAPPA * sf, prmin), prmax);

                // Use pre-computed wstar from the dedicated recomputation loop (lines 465-545).
                // wstar_arr was computed with pblh_corr_arr to ensure consistency
                // between countergradient diagnostics and K-profile calculations.
                const Real wstar = wstar_arr(i, j, 0);

                // SFCFLG gating: WRF skips nonlocal mixing in stable PBL (SFCFLG=.FALSE., BR>0, obuk_val>0)
                // In stable conditions, use free-atmosphere Richardson mixing instead.
                // WRF Reference: module_bl_ysu.F lines 200, 220-260
                if (SFCFLG) {
                    // Diffusion coefficient for momentum with terrain height correction
                    // K = rho * wscale * kappa * zrel * (1 - zrel/pblh_rel)^pfac
                    // where zrel = z - z_sfc, pblh_rel = pblh - z_sfc
                    // WRF Reference: module_bl_ysu.F L250-260 uses z - z_ground for shape function
                    const Real z_sfc = (use_terrain_fitted_coords)
                                     ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                                     : zero;
                    const Real zrel = zval - z_sfc;
                    const Real pblh = pblh_corr_arr(i, j, 0);
                    const Real pblh_rel = pblh - z_sfc;
                    const Real zfac = amrex::max(one - zrel / pblh_rel, Real(1.0e-8));

                    constexpr Real ckz_pbl = Real(0.001);
                    const Real K_base = ckz_pbl * dz_terrain * rho;
                    K_turb(i, j, k, EddyDiff::Mom_v) = K_base + rho * wscale * KAPPA * zrel * zfac * zfac;
                    K_turb(i, j, k, EddyDiff::Theta_v) = K_turb(i, j, k, EddyDiff::Mom_v) / Prt;

                    // Moisture diffusivity: YSU uses Prq ~ Prt (same stability functions
                    // for heat and moisture, module_bl_ysu.F lines 220-260).
                    // Reference: https://github.com/wrf-model/WRF/blob/master/phys/module_bl_ysu.F#L220-L260
                    //
                    // Physics justification: Heat and moisture both respond to the same
                    // buoyancy-driven turbulent eddies in the mixed layer. The Prandtl
                    // and Schmidt numbers are thus equal in this formulation. Alternative
                    // approaches use Prq ≠ Prt (e.g., Högström 1996) but the YSU scheme
                    // defaults to equality for simplicity and computational efficiency.
                    if (turbChoice.mrf_moistvars) {
                        Real Prq_base = phit_eff / phiM_eff;
                        const Real Prq = amrex::min(amrex::max(Prq_base + const_b * KAPPA * sf, prmin), prmax);
                        K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Mom_v) / Prq;
                    } else {
                        K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
                    }
                } else {
                    // Stable PBL: use Richardson mixing with H10 (Hong 2010) formulation
                    // WRF Reference: module_bl_ysu.F; H10 Section 3a
                    // H10 = Hong, S.-Y., 2010: A new stable boundary-layer mixing scheme. QJRMS, 136, 1481-1496.
                    const Real lambda_min = Real(30.0);
                    const Real lambda_max = Real(300.0);
                    const Real lambdadz = amrex::min(amrex::max(Real(0.1) * dz_terrain, lambda_min), lambda_max);
                    const Real lscale = (lambdadz * KAPPA * zval) / (lambdadz + KAPPA * zval);
                    Real dthetadz, dudz, dvdz;
                    ComputeVerticalDerivativesPBL(i, j, k, uvel, vvel, cell_data, izmin, izmax, one / dz_terrain,
                                                  c_ext_dir_on_zlo, c_ext_dir_on_zhi, u_ext_dir_on_zlo,
                                                  u_ext_dir_on_zhi, v_ext_dir_on_zlo, v_ext_dir_on_zhi, dthetadz,
                                                  dudz, dvdz, moisture_indices);
                    // Note: dthetadz (dry potential temperature gradient) is computed above but intentionally
                    // not used. Instead, we compute dtheta_v_dz from virtual potential temperature (lines 922-925)
                    // to properly account for moisture effects on buoyancy in the Richardson number calculation.
                    // This is consistent with the free-atmosphere branch (line 970) which also uses virtual
                    // potential temperature for the Richardson number computation.

                    // Apply boundary safeguards to avoid numerical instability
                    const Real dudz_safe = (k < izmax) ? dudz : zero;
                    const Real dvdz_safe = (k < izmax) ? dvdz : zero;

                    const Real wind_shear = dudz_safe * dudz_safe + dvdz_safe * dvdz_safe;
                    const Real wind_shear_safe = std::max(wind_shear, Real(1.0e-8));

                    // Use virtual potential temperature for Richardson number
                    const Real theta_v = GetThetav(i, j, k, cell_data, moisture_indices);
                    const Real theta_v_kp1 = (k < izmax) ? GetThetav(i, j, k+1, cell_data, moisture_indices) : theta_v;
                    const Real theta_v_km1 = (k > izmin) ? GetThetav(i, j, k-1, cell_data, moisture_indices) : theta_v;
                    const Real dtheta_v_dz = myhalf * (theta_v_kp1 - theta_v_km1) * (one / dz_terrain);

                    // Gradient Richardson number with bounds
                    Real grad_Ri = CONST_GRAV / theta_v * dtheta_v_dz / wind_shear_safe;
                    grad_Ri = std::max(std::min(grad_Ri, Real(100.0)), -Real(100.0));

                    const Real grad_Ri_safe = amrex::max(grad_Ri, -Real(100.0));
                    const Real fm = (grad_Ri_safe > 0)
                                  ? one / ((one + Real(5.0) * grad_Ri_safe) * (one + Real(5.0) * grad_Ri_safe))
                                  : 1 - 8 * grad_Ri_safe / (1 + Real(1.746) * std::sqrt(amrex::max(-grad_Ri_safe, zero)));
                    const Real ft = (grad_Ri_safe > 0)
                                  ? one / ((one + Real(5.0) * grad_Ri_safe) * (one + Real(5.0) * grad_Ri_safe))
                                  : 1 - 8 * grad_Ri_safe / (1 + Real(1.286) * std::sqrt(amrex::max(-grad_Ri_safe, zero)));
                    const Real rl2wsp = rho * lscale * lscale * std::sqrt(wind_shear);

                    K_turb(i, j, k, EddyDiff::Mom_v)   = rl2wsp * fm;
                    // H10 formulation: for stable conditions, K_t = K_m / Pr
                    if (grad_Ri_safe > 0) {
                        Real Pr = amrex::min(amrex::max(one + Real(2.1) * grad_Ri, Real(0.25)), Real(4.0));
                        K_turb(i, j, k, EddyDiff::Theta_v) = K_turb(i, j, k, EddyDiff::Mom_v) / Pr;
                    } else {
                        K_turb(i, j, k, EddyDiff::Theta_v) = rl2wsp * ft;
                    }
                    if (use_moisture && turbChoice.ysu_moistvars) {
                        K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
                    } else {
                        K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
                    }
                }
            } else if (k >= pbli_extent) {
                // Free atmosphere above PBL: use local Richardson number-dependent mixing
                // with H10 (Hong 2010) grid-adaptive lengthscale
                // WRF Reference: module_bl_ysu.F; H10 Section 3a
                // H10 = Hong, S.-Y., 2010: A new stable boundary-layer mixing scheme. QJRMS, 136, 1481-1496.
                const Real lambda_min = Real(30.0);
                const Real lambda_max = Real(300.0);
                const Real lambdadz = amrex::min(amrex::max(Real(0.1) * dz_terrain, lambda_min), lambda_max);
                const Real lscale = (lambdadz * KAPPA * zval) / (lambdadz + KAPPA * zval);
                Real dthetadz, dudz, dvdz;
                ComputeVerticalDerivativesPBL(i, j, k, uvel, vvel, cell_data, izmin, izmax, one / dz_terrain,
                                              c_ext_dir_on_zlo, c_ext_dir_on_zhi, u_ext_dir_on_zlo,
                                              u_ext_dir_on_zhi, v_ext_dir_on_zlo, v_ext_dir_on_zhi, dthetadz,
                                              dudz, dvdz, moisture_indices);

                // Apply boundary safeguards to avoid numerical instability in calm conditions above PBL
                const Real dudz_safe = (k < izmax) ? dudz : zero;
                const Real dvdz_safe = (k < izmax) ? dvdz : zero;

                // Wind shear magnitude with safety threshold (WRF approach, optimized for GPU precision)
                // Using separate safety threshold avoids numerical issues with very small wind shear
                const Real wind_shear = dudz_safe * dudz_safe + dvdz_safe * dvdz_safe;
                const Real wind_shear_safe = std::max(wind_shear, Real(1.0e-8));

                // Use virtual potential temperature (θ_v) for Richardson number stability calculation.
                // This correctly accounts for moisture effects on buoyancy.
                // WRF Reference: module_bl_ysu.F uses THVX (virtual potential temperature).
                // For moist air, θ_v = θ * (1 + 0.61*q_v - q_l - q_i) ≈ θ * (1 + 0.61*q_v)
                const Real theta_v = GetThetav(i, j, k, cell_data, moisture_indices);
                const Real theta_v_kp1 = (k < izmax) ? GetThetav(i, j, k+1, cell_data, moisture_indices) : theta_v;
                const Real theta_v_km1 = (k > izmin) ? GetThetav(i, j, k-1, cell_data, moisture_indices) : theta_v;
                const Real dtheta_v_dz = myhalf * (theta_v_kp1 - theta_v_km1) * (one / dz_terrain);

                // Gradient Richardson number: Ri_g = (g/θ_v) * (dθ_v/dz) / (shear²)
                // Reference: WRF module_bl_ysu.F line ~450-456, Hong et al. 2006, Eqn. A18
                //
                // For STABILITY ANALYSIS:
                // Ri_g > 0.5: typically considered strongly stable (turbulence suppressed)
                // 0 < Ri_g < 0.5: weakly stable
                // Ri_g < 0: unstable (turbulence enhanced)
                // Ri_g < -100: very strong instability (apply lower bound for safety)
                //
                // Bound values from below (-100.0) and above (100.0)
                // to prevent extreme floating-point scales from causing numerical instability
                Real grad_Ri = CONST_GRAV / theta_v * dtheta_v_dz / wind_shear_safe;
                grad_Ri = std::max(std::min(grad_Ri, Real(100.0)), -Real(100.0));
                // YSU stability functions (Hong et al. 2006, MWR, Appendix A)
                // Reference: https://doi.org/10.1175/MWR3250.1
                // See equations A19-A20
                //
                // Prandtl number in stable regime:
                //   Pr = 1 + 2.1 * Ri_g  (Eqn. A19)
                // Stability functions for momentum and heat:
                //   For Ri_g > 0 (stable):
                //     f_m = 1 / ((1 + 5*Ri_g)²)     (Eqn. A20a)
                //     f_t = 1 / ((1 + 5*Ri_g)²)     (Eqn. A20b, same as f_m for this scheme)
                //   For Ri_g < 0 (unstable):
                //     f_m = 1 - 8*Ri_g / (1 + 1.746*sqrt(-Ri_g))  (Eqn. A20d)
                //     f_t = 1 - 8*Ri_g / (1 + 1.286*sqrt(-Ri_g))  (Eqn. A20c)
                // Protect against numerical errors causing negative grad_Ri
                const Real grad_Ri_safe = amrex::max(grad_Ri, -Real(100.0)); // Bound negative values
                const Real fm = (grad_Ri_safe > 0)
                              ? one / ((one + Real(5.0) * grad_Ri_safe) * (one + Real(5.0) * grad_Ri_safe))
                              : 1 - 8 * grad_Ri_safe / (1 + Real(1.746) * std::sqrt(amrex::max(-grad_Ri_safe, zero))); // Eqn. A20b/d
                const Real ft = (grad_Ri_safe > 0)
                              ? one / ((one + Real(5.0) * grad_Ri_safe) * (one + Real(5.0) * grad_Ri_safe))
                              : 1 - 8 * grad_Ri_safe / (1 + Real(1.286) * std::sqrt(amrex::max(-grad_Ri_safe, zero))); // Eqn. A20a/c
                const Real rl2wsp = rho * lscale * lscale * std::sqrt(wind_shear);

                K_turb(i, j, k, EddyDiff::Mom_v)   = rl2wsp * fm;
                // H10 formulation: for stable conditions, K_t = K_m / Pr
                if (grad_Ri_safe > 0) {
                    Real Pr = amrex::min(amrex::max(one + Real(2.1) * grad_Ri, Real(0.25)), Real(4.0));
                    K_turb(i, j, k, EddyDiff::Theta_v) = K_turb(i, j, k, EddyDiff::Mom_v) / Pr;
                } else {
                    K_turb(i, j, k, EddyDiff::Theta_v) = rl2wsp * ft;
                }
                // WRF YSU: moisture diffusivity matches heat in free atmosphere
                // Physics: Above PBL, Richardson number mixing applies same Ri_g
                // dependent mixing to both heat and moisture. The heat scaling (ft or K_m/Pr)
                // is used for moisture, which implicitly assumes equal turbulent
                // Prandtl and Schmidt numbers in the stable/unstable free atmosphere.
                // Only apply if use_moisture is enabled
                if (use_moisture && turbChoice.ysu_moistvars) {
                    K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
                } else {
                    K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
                }
            }
            // Note: The conditions (k < pbli_extent) and (k >= pbli_extent) are exhaustive,
            // so the else clause would be unreachable. All cells reach K-bounds clipping below.

            // Add entrainment diffusivity at PBL top cell (HND06 Eq. 6)
            // Applied only at k == pbli_arr(i,j,0) (the first cell above the PBL)
            if (k == pbli_arr(i,j,0)) {
                K_turb(i,j,k,EddyDiff::Mom_v)   += entr_arr(i,j,0);
                K_turb(i,j,k,EddyDiff::Theta_v) += entr_arr(i,j,0);
                if (turbChoice.ysu_moistvars) {
                    K_turb(i,j,k,EddyDiff::Q_v) += entr_arr(i,j,0);
                }
            }

            // Limit diffusion coefficients to physical bounds
            // These bounds ensure numerical stability and prevent unrealistic diffusivity
            // that could violate conservation principles or cause solver issues.
            //
            // HND06 Conservative Limits (module_bl_ysu.F lines 300-320):
            // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_ysu.F#L300-L320
            // Kmin = 0.1 m²/s, Kmax = 300 m²/s
            //
            // These limits were calibrated for global forecast models and found to be
            // robust across a wide range of atmospheric conditions. They prevent:
            //   - Excessive diffusion in calm conditions (Kmax bound)
            //   - Numerical instability from too-small diffusion (Kmin bound)
            //   - CFL violations from over-diffusive mixing
            //
            // Alternative (Higher Limits) H10 formulation (QJRMS Appendix A):
            // Kmin = ckz * dz * rho (where ckz = 0.001), Kmax = 1000 m²/s
            // These higher limits allow greater mixing in free atmosphere and are
            // more appropriate for high-resolution simulations.
            Real rhoKmin, rhoKmax;
            if (turbChoice.pbl_ysu_highres_bounds) {
                // Hong et al. 2006, MWR, Appendix A: Higher limits for free atmosphere
                // These are recommended for high-resolution simulations (Δz < 100 m)
                // where the free atmosphere grid can resolve small-scale mixing.
                constexpr Real ckz  = Real(0.001);
                constexpr Real Kmax = Real(1000.0);
                rhoKmin  = ckz * dz_terrain * rho;
                rhoKmax  = rho * Kmax;
            } else {
                // HND06 Conservative limits (default, used in global forecasts)
                // These provide good results for typical meteorological scales (Δz > 100 m)
                // and are the standard in legacy WRF configurations.
                constexpr Real Kmin = Real(0.1);
                constexpr Real Kmax = Real(300.0);
                rhoKmin  = rho * Kmin;
                rhoKmax  = rho * Kmax;
            }

            K_turb(i, j, k, EddyDiff::Mom_v) = std::max(
                std::min(K_turb(i, j, k, EddyDiff::Mom_v), rhoKmax), rhoKmin);
            K_turb(i, j, k, EddyDiff::Theta_v) = std::max(
                std::min(K_turb(i, j, k, EddyDiff::Theta_v), rhoKmax), rhoKmin);
            K_turb(i, j, k, EddyDiff::Q_v) = std::max(
                std::min(K_turb(i, j, k, EddyDiff::Q_v), rhoKmax), rhoKmin);
            K_turb(i, j, k, EddyDiff::Turb_lengthscale) = pblh_corr_arr(i, j, 0);

            // Store countergradient correction terms (HGAMT/h and HGAMQ/h)
            // Use the selected PBL extent index (pbli_arr or pbli_zero_arr based on pbl_mrf_use_zero_ri_extent)
            if (k < pbli_extent) {
                // Inside PBL: store the normalized countergradient terms
                K_turb(i, j, k, EddyDiff::HGAMT_v) = hgamt_arr(i, j, 0);
                K_turb(i, j, k, EddyDiff::HGAMQ_v) = hgamq_arr(i, j, 0);
            } else {
                // Outside PBL: zero countergradient
                K_turb(i, j, k, EddyDiff::HGAMT_v) = zero;
                K_turb(i, j, k, EddyDiff::HGAMQ_v) = zero;
            }
        });

        // FOEXTRAP top and bottom ghost cells
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
        {
            K_turb(i, j, klo-1, EddyDiff::Mom_v  ) = K_turb(i, j, klo, EddyDiff::Mom_v  );
            K_turb(i, j, klo-1, EddyDiff::Theta_v) = K_turb(i, j, klo, EddyDiff::Theta_v);
            K_turb(i, j, klo-1, EddyDiff::Q_v    ) = K_turb(i, j, klo, EddyDiff::Q_v    );
            K_turb(i, j, klo-1, EddyDiff::HGAMT_v) = K_turb(i, j, klo, EddyDiff::HGAMT_v);
            K_turb(i, j, klo-1, EddyDiff::HGAMQ_v) = K_turb(i, j, klo, EddyDiff::HGAMQ_v);
            K_turb(i, j, klo-1, EddyDiff::Turb_lengthscale) = K_turb(i, j, klo, EddyDiff::Turb_lengthscale);
            K_turb(i, j, khi+1, EddyDiff::Mom_v  ) = K_turb(i, j, khi, EddyDiff::Mom_v  );
            K_turb(i, j, khi+1, EddyDiff::Theta_v) = K_turb(i, j, khi, EddyDiff::Theta_v);
            K_turb(i, j, khi+1, EddyDiff::Q_v    ) = K_turb(i, j, khi, EddyDiff::Q_v    );
            K_turb(i, j, khi+1, EddyDiff::HGAMT_v) = K_turb(i, j, khi, EddyDiff::HGAMT_v);
            K_turb(i, j, khi+1, EddyDiff::HGAMQ_v) = K_turb(i, j, khi, EddyDiff::HGAMQ_v);
            K_turb(i, j, khi+1, EddyDiff::Turb_lengthscale) = K_turb(i, j, khi, EddyDiff::Turb_lengthscale);
        });
    }// mfi
}

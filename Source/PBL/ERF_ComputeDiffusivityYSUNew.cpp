#include "ERF_SurfaceLayer.H"
#include "ERF_DirectionSelector.H"
#include "ERF_Diffusion.H"
#include "ERF_Constants.H"
#include "ERF_TurbStruct.H"
#include "ERF_PBLModels.H"
#include "ERF_TileNoZ.H"
#include "ERF_MoistUtils.H"

#ifdef ERF_USE_WINDFARM
#include "ERF_WindFarm.H"
#endif

using namespace amrex;

/**
 * Compute vertical eddy viscosity coefficients using the Yonsei University (YSU) boundary layer scheme.
 *
 * @param[in] xvel X-velocity field.
 * @param[in] yvel Y-velocity field.
 * @param[in] cons_in Input conservative variables.
 * @param[out] eddyViscosity MultiFab to store computed eddy viscosity and countergradient terms.
 * @param[in] geom Geometry used for grid spacings and domain extent.
 * @param[in] turbChoice Turbulence model configuration and parameters.
 * @param[in] SurfLayer Pointer to surface layer data.
 * @param[in] use_terrain_fitted_coords Use terrain-fitted coordinates if true.
 * @param[in] use_moisture Include moisture in the diffusivity calculation.
 * @param[in] level Current AMR level.
 * @param[in] bc_ptr Boundary condition records.
 * @param[in] vert_only Reserved flag for vertical-only computation.
 * @param[in] z_phys_nd Nodal physical height field.
 * @param[in] z_phys_cc Cell-centered physical height field.
 * @param[in] moisture_indices Indices for moisture variables in the state vector.
 * @param[in] qheating_rates Optional heating rates for cloud-top mixing.
 */
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
                       const std::unique_ptr<MultiFab>& z_phys_cc,
                       const MoistureComponentIndices& moisture_indices,
                       const MultiFab* qheating_rates)
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

    4. Explicit Entrainment at PBL Top (WRF bl_ysu.F90 lines 831-925):
       K_entr = ρ * |we| * dz  at k = kpbl  (applied only when we < 0)
       wm3   = wstar³ + 5*u_*³
       bfxpbl = -0.15 * θv(klo)/g * wm3 / pblh
       dthvx = max(θv(kpbl+1) - θv(kpbl), 1e-2)
       we    = max(bfxpbl/dthvx, -sqrt(wm3^(2/3)))

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
    */

    // Domain extent in z-dir
    int klo = geom.Domain().smallEnd(2);
    int khi = geom.Domain().bigEnd(2);
    const int izmin = klo;
    const int izmax = khi;

    const Real dz     = geom.CellSize(2);
    const Real dz_inv = geom.InvCellSize(2);
    const auto& dxInv = geom.InvCellSizeArray();

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
        const Box xybx = PerpendicularBox<ZDir>(gbx, IntVect{0, 0, 0});
        FArrayBox pbl_height_corrector(xybx, 1, The_Async_Arena());
        IArrayBox pbl_index(xybx, 1, The_Async_Arena());
        IArrayBox pbl_index_zero_ri(xybx, 1, The_Async_Arena());  // Index for zero-Ri diagnostic pass
        FArrayBox hgamt_fab(xybx, 1, The_Async_Arena());  // Store HGAMT/h (normalized countergradient)
        FArrayBox hgamq_fab(xybx, 1, The_Async_Arena());  // Store HGAMQ/h (normalized countergradient)
        FArrayBox hgamu_fab(xybx, 1, The_Async_Arena());  // Store HGAMU (YSU u countergradient)
        FArrayBox hgamv_fab(xybx, 1, The_Async_Arena());  // Store HGAMV (YSU v countergradient)
        FArrayBox wstar_fab(xybx, 1, The_Async_Arena());  // Convective velocity scale computed with pblh_corr
        FArrayBox vpert_fab(xybx, 1, The_Async_Arena());  // Virtual temperature perturbation VPERT for Pass 3
        FArrayBox entr_fab(xybx, 1, The_Async_Arena());   // Per-column entrainment diffusivity
        IArrayBox cloud_top_fab(xybx, 1, The_Async_Arena());  // Cloud-top cell index for top-down mixing
        FArrayBox wstar3_down_fab(xybx, 1, The_Async_Arena());  // Top-down convective velocity scale cubed
        FArrayBox sflux_fab(xybx, 1, The_Async_Arena());        // Virtual kinematic heat flux (sensible + latent)
        FArrayBox wstar3_fab(xybx, 1, The_Async_Arena());       // Convective velocity scale cubed per column
        FArrayBox zol1_fab(xybx, 1, The_Async_Arena());         // Monin-Obukhov stability parameter at first level
        FArrayBox sfcflg_fab(xybx, 1, The_Async_Arena());       // Surface stability flag: 1=unstable/neutral, 0=stable
        // REQUIRED: zero-initialize all FArrayBoxes before GPU kernels read them.
        // The_Async_Arena() does not zero memory; reading uninitialized data
        // corrupts rib_enhan_arr (via vpert_arr) and SFCFLG (via sflux_arr).
        vpert_fab.setVal<RunOn::Device>(zero);
        sflux_fab.setVal<RunOn::Device>(zero);
        wstar3_fab.setVal<RunOn::Device>(zero);
        wstar3_down_fab.setVal<RunOn::Device>(zero);
        sfcflg_fab.setVal<RunOn::Device>(zero);
        hgamt_fab.setVal<RunOn::Device>(zero);
        hgamq_fab.setVal<RunOn::Device>(zero);
        hgamu_fab.setVal<RunOn::Device>(zero);
        hgamv_fab.setVal<RunOn::Device>(zero);
        wstar_fab.setVal<RunOn::Device>(zero);
        entr_fab.setVal<RunOn::Device>(zero);
        pbl_height_corrector.setVal<RunOn::Device>(zero);
        const auto& pblh_corr_arr   = pbl_height_corrector.array();
        const auto& pbli_arr        = pbl_index.array();
        const auto& pbli_zero_arr   = pbl_index_zero_ri.array();  // Zero-Ri diagnostic PBL index
        const auto& hgamt_arr       = hgamt_fab.array();
        const auto& hgamq_arr       = hgamq_fab.array();
        const auto& hgamu_arr       = hgamu_fab.array();
        const auto& hgamv_arr       = hgamv_fab.array();
        const auto& wstar_arr       = wstar_fab.array();  // Stored convective velocity for use in K-profile
        const auto& vpert_arr       = vpert_fab.array();  // Stored VPERT for Pass 3 diagnostic loop
        const auto& entr_arr        = entr_fab.array();   // Stored entrainment diffusivity
        const auto& cloud_top_arr   = cloud_top_fab.array();  // Cloud-top cell index for top-down mixing
        const auto& wstar3_down_arr = wstar3_down_fab.array();  // Top-down wstar cubed
        const auto& sflux_arr       = sflux_fab.array();       // Virtual kinematic heat flux
        const auto& wstar3_arr      = wstar3_fab.array();      // Convective velocity scale cubed
        const auto& zol1_arr        = zol1_fab.array();        // Monin-Obukhov stability parameter
        const auto& sfcflg_arr      = sfcflg_fab.array();      // Surface stability flag (1=unstable/neutral, 0=stable)

        // Get some data in arrays
        const auto& cell_data = cons_in.const_array(mfi);
        const auto& uvel = xvel.const_array(mfi);
        const auto& vvel = yvel.const_array(mfi);

        //const Real Ribcr = turbChoice.pbl_mrf_Ribcr;  // Now surface-type-dependent in YSU
        //const Real f0 = turbChoice.pbl_mrf_coriolis_freq;
        const auto& u_star_arr = SurfLayer->get_u_star(level)->const_array(mfi);
        const auto& t_star_arr = SurfLayer->get_t_star(level)->const_array(mfi);
        const auto& q_star_arr = SurfLayer->get_q_star(level)->const_array(mfi);
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
        const PBLDerivativeDzInv_T pbl_derivative_dz_inv{z_phys_cc->const_array(mfi)};
        // Get qheating_rates if provided (for LW radiation coupling to top-down mixing)
        const Array4<Real const> qheat_arr = (qheating_rates != nullptr)
                                             ? qheating_rates->const_array(mfi)
                                             : Array4<Real const>{};
        const bool has_qheating_rates = (qheating_rates != nullptr);


        // ========================================================================
        // PRE-COMPUTE Rib ARRAYS FOR GPU PERFORMANCE
        // ========================================================================
        // Pre-compute both base (without VPERT) and enhanced (with VPERT) Richardson numbers
        // in a single vectorized kernel to reduce divergent control flow in subsequent PBLH loops.
        // This improves GPU performance by eliminating serial while-loops in favor of simple
        // array lookups during the three PBLH pass loops.
        //
        const bool enable_ysu_liquid_theta = turbChoice.enable_ysu_liquid_theta;
        FArrayBox rib_base_fab(gbx, 1, The_Async_Arena());    // Rib with base t_layer_v
        FArrayBox rib_enhan_fab(gbx, 1, The_Async_Arena());   // Rib with enhanced t_layer_v+VPERT
        const auto& rib_base_arr  = rib_base_fab.array();
        const auto& rib_enhan_arr = rib_enhan_fab.array();

        BL_PROFILE_VAR("YSUNew_Rib_Precompute", prof_rib_precomp);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
          const amrex::Real t_layer   = t10av_arr(i,j,0);
          const amrex::Real q_frac    = use_moisture ? q10av_arr(i,j,0) : zero;
          const amrex::Real t_layer_v = t_layer * (one + epsv * q_frac);
          // On first call, vpert_arr is initialized to zero
          const amrex::Real t_enh     = t_layer_v + vpert_arr(i,j,0);
          const amrex::Real z_sfc     = (use_terrain_fitted_coords)
                                      ? Compute_Zrel_AtCellCenter(i,j,klo,z_nd_arr) : zero;
          const amrex::Real zval      = (use_terrain_fitted_coords)
                                      ? Compute_Zrel_AtCellCenter(i,j,k,z_nd_arr)
                                      : (k + myhalf) * dz;
          const amrex::Real zrel      = amrex::max(zval - z_sfc, amrex::Real(1.0e-4));
          const amrex::Real theta_v   = (use_moisture && enable_ysu_liquid_theta)
                                      ? GetThetavl(i,j,k,cell_data,moisture_indices)
                                      : GetThetav(i,j,k,cell_data,moisture_indices);
          const amrex::Real theta_v_klo = (use_moisture && enable_ysu_liquid_theta)
                                        ? GetThetavl(i,j,klo,cell_data,moisture_indices)
                                        : GetThetav(i,j,klo,cell_data,moisture_indices);
          const amrex::Real ws2_raw   = fourth * ((uvel(i,j,k)+uvel(i+1,j,k))*(uvel(i,j,k)+uvel(i+1,j,k))
                                      + (vvel(i,j,k)+vvel(i,j+1,k))*(vvel(i,j,k)+vvel(i,j+1,k)));
          const amrex::Real ws2       = amrex::max(ws2_raw, amrex::Real(1.0));
          rib_base_arr(i,j,k)  = CONST_GRAV * zrel * (theta_v - t_layer_v) / (ws2 * theta_v_klo);
          rib_enhan_arr(i,j,k) = CONST_GRAV * zrel * (theta_v - t_enh)     / (ws2 * theta_v_klo);
        });
        BL_PROFILE_VAR_STOP(prof_rib_precomp);

        // ========================================================================
        // PRE-COMPUTE surface flux (sflux), stability flag (sfcflg), and convective
        // velocity scale (wstar3) per column, following WRF-YSU formulation
        // ========================================================================
        // WRF Reference: module_bl_ysu.F lines 610-680
        // These quantities are needed for:
        //   - Stability-dependent functions throughout the PBLH passes
        //   - Countergradient flux calculations (HGAMT, HGAMQ)
        //   - K-profile computations
        // Computing them once per column improves GPU performance and ensures consistency.
        //
        BL_PROFILE_VAR("YSUNew_SurfFlux_Precompute", prof_sflux_precomp);
        ParallelFor(xybx, [=, zero_d=zero] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            const Real rho_sfc  = cell_data(i, j, klo, Rho_comp);
            const Real t_layer  = t10av_arr(i, j, 0);      // Dry potential temperature at z1

            // WRF: rhox = psfc/(rd*tx(1)*tvcon), approximate with cell-bottom density
            const Real rhox = rho_sfc;

            // Virtual temperature factor: ep1 = Rv/Rd - 1 = 0.608
            constexpr Real ep1 = amrex::Real(0.608);
            constexpr Real cp_air = amrex::Real(1004.0);

            // WRF bl_ysu.F90 line 613: sflux = hfx/(rho*cp) + qfx/rho*ep1*theta
            // Approximate hfx and qfx from MOST quantities:
            //   hfx [W/m^2] = -rho * cp * u_star * t_star
            //   qfx [kg/m^2/s] = -rho * u_star * q_star
            const Real ustar  = u_star_arr(i, j, 0);
            const Real tstar  = t_star_arr(i, j, 0);
            const Real qstar  = use_moisture ? q_star_arr(i, j, 0) : zero;
            const Real hfx    = -rhox * cp_air * ustar * tstar;   // W/m^2
            const Real qfx    = -rhox * ustar * qstar;            // kg/m^2/s
            const Real sflux  = hfx / (rhox * cp_air) + qfx / rhox * ep1 * t_layer;
            sflux_arr(i, j, 0) = sflux;

            // WRF sfcflg: TRUE when unstable/neutral (sflux > 0), FALSE when stable
            // WRF bl_ysu.F90 line 615: if(br(i).gt.0.0) sfcflg(i) = .false.
            // In ERF, br > 0 corresponds to sflux < 0 (stable), so:
            const bool sfcflg = (sflux > zero);
            sfcflg_arr(i, j, 0) = sfcflg ? one : zero;  // Store as scalar field for unified use

            // WRF: govrth = g/thx(1) — uses DRY theta at first level
            const Real govrth = CONST_GRAV / t_layer;

            // WRF bl_ysu.F90 lines 664-668: wstar3 from sflux
            // First-guess pblh needed for wstar3 — use nominal 500m for initial estimate
            // wstar3 is refined in subsequent loop after corrector pblh is available
            const Real pblh_guess = amrex::Real(500.0);
            const Real bfx0   = amrex::max(sflux, zero_d);
            const Real wstar3 = govrth * bfx0 * pblh_guess;
            wstar3_arr(i, j, 0) = wstar3;

            // WRF bl_ysu.F90 line 676: wscale = (ust3 + phifac*karman*wstar3*0.5)^(1/3)
            // The 0.5 factor is the column-average of (1-zfac) over the PBL depth
            constexpr Real phifac = amrex::Real(8.0);
            const Real ust3   = ustar * ustar * ustar;
            Real wscale = std::cbrt(ust3 + phifac * KAPPA * wstar3 * amrex::Real(0.5));
            // WRF bounds: min(wscale, ust*16), max(wscale, ust/5)
            wscale = amrex::min(wscale, ustar * amrex::Real(16.0));
            wscale = amrex::max(wscale, ustar / amrex::Real(5.0));
            wstar_arr(i, j, 0) = wscale;  // Store column-average wscale for HGAMT/HGAMQ

            // zol1: Monin-Obukhov stability parameter at first level
            // WRF bl_ysu.F90 lines 651-662: zol1 = max(br*fm*fm/fh, rimin)
            // Approximate: zol1 = z1 / L_obuk
            Real obuk_val = l_obuk_arr(i, j, 0);
            if (std::abs(obuk_val) < amrex::Real(1.0e-10))
                obuk_val = (obuk_val >= zero) ? amrex::Real(1.0e-10) : amrex::Real(-1.0e-10);
            const Real zl1 = (use_terrain_fitted_coords)
                           ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                           : (klo + myhalf) * dz;
            Real zol1 = zl1 / obuk_val;
            zol1 = amrex::max(zol1, amrex::Real(-100.0));
            if (sfcflg) {
                zol1 = amrex::min(zol1, amrex::Real(-1.0e-8));
            } else {
                zol1 = amrex::max(zol1, amrex::Real(1.0e-8));
            }
            zol1_arr(i, j, 0) = zol1;
        });
        BL_PROFILE_VAR_STOP(prof_sflux_precomp);

        // ========================================================================
        // PASS 1 — PREDICTOR: Compute initial PBL index with surface-type-dependent
        // critical Richardson number (Ribcr).
        // ========================================================================
        // WRF Reference: module_bl_ysu.F (Hong, Noh & Dudhia 2006, lines 1-100)
        // Ribcr = 0.25 over land; Rossby-number-dependent over water (line 150-170).
        // Rib = (g*z/θv0) * (θv(z) - θv_surf) / ws² (lines 180-200)
        // References: Hong et al., Mon. Wea. Rev., 134, 2318-2341 (2006)
        //

        BL_PROFILE_VAR("YSUNew_PBLH_Passes", prof_pblh);
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            // Determine surface-type-dependent critical Richardson number
            // Over land: Ribcr = 0.25; over water: Ribcr depends on Rossby number
            Real Ribcr;
            {
                bool over_land = (!lmask_arr) || (lmask_arr(i, j, 0) == 1);
                if (over_land) {
                    Ribcr = turbChoice.pbl_ysu_land_Ribcr;
                } else {
                    const Real z0 = z0_arr(i, j, 0);
                    const Real ws_layer = ws10av_arr(i, j, 0);
                    const Real Rossby = ws_layer / (turbChoice.pbl_ysu_coriolis_freq * amrex::max(z0, Real(1.0e-4)));
                    Ribcr = amrex::min(Real(0.16) * std::pow(Real(1.0e-7) * Rossby, -Real(0.18)), Real(0.3));
                }
            }

            // Scan rib_base_arr(i,j,klo..khi) to find first crossing of Ribcr
            // GAP 9: Seed from surface bulk Richardson number (WRF bl_ysu.F90 lines 620-621)
            // WRF seeds the loop with br (surface bulk Ri) so klo is checked first.
            int kpbl = klo;
            Real Rib = rib_base_arr(i,j,klo);  // Seed from surface level
            bool above_critical = (Rib >= Ribcr);  // Check if already above critical

            // Scan from klo+1 if not already above critical
            for (int kk = klo+1; !above_critical && kk <= khi; ++kk) {
                if (rib_base_arr(i,j,kk) >= Ribcr) {
                    kpbl = kk;
                    above_critical = true;
                    break;
                }
                kpbl = kk;  // keep updating so kpbl = khi if never exceeded
            }
            pbli_arr(i, j, 0) = kpbl;
        });
        BL_PROFILE_VAR_STOP(prof_pblh);

        //
        // FIRST ITERATION: Compute HGAMT/HGAMQ/VPERT using Pass 1 PBL height
        // ========================================================================
        // This produces initial estimates of countergradient corrections and VPERT,
        // which will then be used to recompute rib_enhan_arr for a more accurate
        // corrector pass.
        //

        // (Reuse existing corrector loop code below, first pass...)

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
        //constexpr Real prmin = Real(0.5);
        //constexpr Real prmax = Real(4.0);
        constexpr Real GAMCRT = Real(3.0);    // WRF GAMCRT: max heat countergradient (K)
        constexpr Real GAMCRQ = Real(2.e-3);  // WRF GAMCRQ: max moisture countergradient (kg/kg)
        const bool enable_ysu_countergradient = turbChoice.enable_ysu_countergradient;
        const bool enable_ysu_unbounded_vpert = turbChoice.enable_mrf_unbounded_vpert;
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            // Determine surface-type-dependent critical Richardson number (same as Pass 1)
            // Over land: Ribcr = 0.25; over water: Ribcr depends on Rossby number
            Real Ribcr;
            {
                bool over_land = (!lmask_arr) || (lmask_arr(i, j, 0) == 1);
                if (over_land) {
                    Ribcr = turbChoice.pbl_ysu_land_Ribcr;
                } else {
                    const Real z0 = z0_arr(i, j, 0);
                    const Real ws_layer = ws10av_arr(i, j, 0);
                    const Real Rossby = ws_layer / (turbChoice.pbl_ysu_coriolis_freq * amrex::max(z0, Real(1.0e-4)));
                    Ribcr = amrex::min(Real(0.16) * std::pow(Real(1.0e-7) * Rossby, -Real(0.18)), Real(0.3));
                }
            }

            /*
            const Real t_layer  = t10av_arr(i, j, 0);
            const Real moisture_fraction = use_moisture ? q10av_arr(i, j, 0) : zero;
            const Real t_layer_v = t_layer * (one + epsv * moisture_fraction);
            */

            const amrex::Real z_sfc_col = (use_terrain_fitted_coords)
                                         ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                                         : zero;

            int kpbl = klo;
            Real zval0 = zero, Rib0 = zero;
            // GAP 9: Seed from surface bulk Richardson number (WRF bl_ysu.F90 lines 620-621)
            // WRF seeds the loop with br (surface bulk Ri) so klo is checked first.
            Real Rib = rib_enhan_arr(i,j,klo);  // Seed from surface level
            zval0 = (use_terrain_fitted_coords)
                  ? Compute_Zrel_AtCellCenter(i,j,klo,z_nd_arr)
                  : (klo + myhalf) * dz;
            Rib0 = Rib;
            bool above_critical = (Rib >= Ribcr);  // Check if already above critical

            // Scan rib_enhan_arr to find first crossing of Ribcr with linear interpolation
            for (int kk = klo+1; !above_critical && kk <= khi; ++kk) {
                if (rib_enhan_arr(i,j,kk) >= Ribcr) { kpbl = kk; above_critical = true; break; }
                zval0 = (use_terrain_fitted_coords)
                      ? Compute_Zrel_AtCellCenter(i,j,kk,z_nd_arr)
                      : (kk + myhalf) * dz;
                Rib0 = rib_enhan_arr(i,j,kk);
                kpbl = kk;
            }

            // Compute bounds for pblh clamping
            const Real z_sfc = (use_terrain_fitted_coords)
                              ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                              : zero;
            const Real dz_terrain = (use_terrain_fitted_coords)
                                   ? (Compute_Zrel_AtCellCenter(i, j, klo + 1, z_nd_arr) - z_sfc)
                                   : dz;
            const Real z_max = (use_terrain_fitted_coords)
                             ? Compute_Zrel_AtCellCenter(i, j, khi, z_nd_arr)
                             : (khi + myhalf) * dz;
            const Real pblh_max = Real(0.9) * z_max;

            amrex::Real pblh_min;
            if (turbChoice.enable_ysu_terrain_pblh_floor && use_terrain_fitted_coords) {
                const amrex::Real met_h_klo = Compute_h_zeta_AtCellCenter(i, j, klo, dxInv, z_nd_arr);
                const amrex::Real dz0 = met_h_klo / dz_inv;
                pblh_min = amrex::max(z_sfc_col + Real(0.5)*dz0, Real(10.0));
            } else {
                pblh_min = amrex::max(z_sfc + Real(0.5) * dz_terrain, Real(10.0));
            }

            // Linear interpolation from previous level if crossing detected
            if (kpbl < khi && rib_enhan_arr(i,j,kpbl) >= Ribcr) {
                const Real zval = (use_terrain_fitted_coords)
                                ? Compute_Zrel_AtCellCenter(i,j,kpbl,z_nd_arr)
                                : (kpbl + myhalf) * dz;
                Rib = rib_enhan_arr(i,j,kpbl);
                Real pblh_interp = zval0 + (zval - zval0) / (Rib - Rib0) * (Ribcr - Rib0);
                pblh_corr_arr(i, j, 0) = amrex::max(amrex::min(pblh_interp, pblh_max), pblh_min);
            } else {
                pblh_corr_arr(i, j, 0) = pblh_min;
            }
            pbli_arr(i, j, 0) = kpbl;

            // Initialize countergradient fluxes for now
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
        const bool enable_ysu_sat_limiter = turbChoice.enable_ysu_sat_limiter;

        // ========================================================================
        // Cloud-top detection for top-down mixing (H10 Section 3b)
        // ========================================================================
        const bool enable_ysu_topdown = turbChoice.enable_ysu_topdown;
        const amrex::Real ysu_qcloud_threshold = turbChoice.ysu_qcloud_threshold;

        if (enable_ysu_topdown && use_moisture) {
            ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
            {
                cloud_top_arr(i, j, 0) = -1;  // Initialize: no cloud found

                // Search from PBL top downward for cloud-top cell
                int kpbl = pbli_arr(i, j, 0);
                for (int kk = kpbl - 1; kk >= klo; --kk) {
                    amrex::Real qc_kk = zero, qi_kk = zero;
                    if (moisture_indices.qc >= 0)
                        qc_kk = cell_data(i, j, kk, moisture_indices.qc) / cell_data(i, j, kk, Rho_comp);
                    if (moisture_indices.qi >= 0)
                        qi_kk = cell_data(i, j, kk, moisture_indices.qi) / cell_data(i, j, kk, Rho_comp);

                    if (qc_kk + qi_kk > ysu_qcloud_threshold) {
                        cloud_top_arr(i, j, 0) = kk;
                        break;
                    }
                }
            });
        } else {
            // Disable cloud-top detection
            ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
            {
                cloud_top_arr(i, j, 0) = -1;
            });
        }

        ParallelFor(xybx, [=, zero_d=zero, one_d=one, CONST_GRAV_d=CONST_GRAV] AMREX_GPU_DEVICE(int i, int j, int) noexcept
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

            // Compute top-down convective velocity scale wstar3_down (H10 Eq. 12)
            // This represents turbulence driven by radiative cooling at cloud top
            amrex::Real wstar3_down = zero;
            if (enable_ysu_topdown && cloud_top_arr(i, j, 0) >= klo) {
                int k_cloud_top = cloud_top_arr(i, j, 0);

                // Longwave cooling flux at cloud top (W/m^2 equivalent in K*m/s)
                // When RRTMGP is active, read from qheating_rates if available.
                // Otherwise use zero (feature disabled when radiation is off).
                // Compute LRAD by integrating the LW heating rate (component index 1)
                // over the PBL column from surface (klo) to cloud top (k_cloud_top)
                amrex::Real LRAD = zero;
                if (has_qheating_rates) {
                    // Integrate LW heating rate from surface to cloud top
                    // qheating_rates contains heating rates (K/s); negative values represent cooling
                    for (int kk = klo; kk <= k_cloud_top; ++kk) {
                        amrex::Real ldz = (use_terrain_fitted_coords)
                                        ? (z_nd_arr(i, j, kk+1) - z_nd_arr(i, j, kk))
                                        : dz;
                        // Accumulate cooling (negative heating) integrated over height
                        // LW component is at index 1
                        LRAD += -qheat_arr(i, j, kk, 1) * ldz;
                    }
                }

                // Top-down convective velocity (H10 Eq. 12):
                // wstar_down^3 = g/theta * LRAD/(rho*cp) * pblh
                const amrex::Real t_local = cell_data(i, j, k_cloud_top, RhoTheta_comp)
                                          / cell_data(i, j, k_cloud_top, Rho_comp);
                const amrex::Real rho_ct  = cell_data(i, j, k_cloud_top, Rho_comp);
                wstar3_down = amrex::max(CONST_GRAV_d / t_local * LRAD
                                        / (rho_ct * amrex::Real(1004.0))
                                        * pblh_corr_arr(i, j, 0), zero_d);
            }
            wstar3_down_arr(i, j, 0) = wstar3_down;

            // Recompute wstar3 and wscale with corrected pblh_corr_arr
            const Real bfx0_corr = amrex::max(sflux_arr(i, j, 0), zero_d);
            const Real t_dry = t10av_arr(i, j, 0);
            const Real wstar3_corr = (CONST_GRAV / t_dry) * bfx0_corr * pblh_corr_arr(i, j, 0);
            wstar3_arr(i, j, 0) = wstar3_corr;

            // Recompute wscale with corrected pblh:
            const Real ust3 = u_star_arr(i,j,0) * u_star_arr(i,j,0) * u_star_arr(i,j,0);
            Real wscale_corr = std::cbrt(ust3 + amrex::Real(8.0) * KAPPA * wstar3_corr * amrex::Real(0.5));
            wscale_corr = amrex::min(wscale_corr, u_star_arr(i,j,0) * amrex::Real(16.0));
            wscale_corr = amrex::max(wscale_corr, u_star_arr(i,j,0) / amrex::Real(5.0));
            wstar_arr(i, j, 0) = wscale_corr;

            // Compute HGAMT with corrected wscale and WRF convention
            // WRF bl_ysu.F90 line 687: gamfac = bfac/rhox/wscale
            // WRF bl_ysu.F90 line 688: hgamt(i) = min(gamfac*hfx(i)/cp, gamcrt)
            // Use sflux-based sfcflg for consistency with the precompute kernel (Fix 3)
            bool SFCFLG = (sfcflg_arr(i, j, 0) > zero);  // TRUE when unstable/neutral (sflux > 0)
            const Real rho_sfc = cell_data(i, j, klo, Rho_comp);
            const Real hfx_col = -rho_sfc * amrex::Real(1004.0) * u_star_arr(i,j,0) * t_star_arr(i,j,0);
            const Real gamfac = const_b / (rho_sfc * wscale_corr);
            const Real HGAMT = (SFCFLG && enable_ysu_countergradient)
                             ? amrex::min(gamfac * hfx_col / amrex::Real(1004.0), GAMCRT)
                             : zero;

            // Compute HGAMQ with corrected wscale and WRF convention
            // WRF bl_ysu.F90 line 689: hgamq(i) = min(gamfac*qfx(i), gamcrq)
            const Real qfx_col  = -rho_sfc * u_star_arr(i,j,0) * q_star_arr(i,j,0);
            const Real HGAMQ_raw = (SFCFLG && use_moisture && enable_ysu_countergradient)
                                  ? gamfac * qfx_col
                                  : zero;
            Real HGAMQ = amrex::min(amrex::max(HGAMQ_raw, zero_d), GAMCRQ);

            // Land/water surface discrimination
            if (lmask_arr && SFCFLG && use_moisture && enable_ysu_countergradient) {
                bool is_land = (lmask_arr(i,j,0) == 1);
                if (!is_land) HGAMQ = zero;
            }

            // Saturation-Aware HGAMQ limiter
            if (enable_ysu_sat_limiter && moisture_indices.qv >= 0 && SFCFLG && use_moisture && enable_ysu_countergradient) {
                Real qv_klo = cell_data(i, j, klo, moisture_indices.qv) / cell_data(i, j, klo, Rho_comp);
                Real T_klo = getTgivenRandRTh(cell_data(i, j, klo, Rho_comp),
                                              cell_data(i, j, klo, RhoTheta_comp),
                                              qv_klo);
                Real p_klo = getPgivenRTh(cell_data(i, j, klo, RhoTheta_comp), qv_klo) * Real(0.01);
                Real qsat_klo = zero;
                erf_qsatw(T_klo, p_klo, qsat_klo);
                Real rh_klo = (qsat_klo > Real(1.0e-10)) ? (qv_klo / qsat_klo) : zero;
                if (rh_klo > Real(0.95)) {
                    Real rh_scaling = amrex::max(zero_d, (one - rh_klo) / Real(0.05));
                    HGAMQ *= rh_scaling;
                }
            }
            // Store HGAMT/h and HGAMQ/h for implicit solver (normalized by corrected PBL height)
            // Also compute VPERT = max(HGAMT + 0.61*θ*HGAMQ, 0) for Pass 3 use
            // WRF Reference: module_bl_ysu.F lines 230-240 (THERMAL update with VPERT)
            if (pbli_arr(i, j, 0) <= klo + 1) {
                hgamt_arr(i, j, 0) = zero;
                hgamq_arr(i, j, 0) = zero;
                hgamu_arr(i, j, 0) = zero;
                hgamv_arr(i, j, 0) = zero;
                vpert_arr(i, j, 0) = zero;
            } else {
                const Real pblh = pblh_corr_arr(i, j, 0);
                hgamt_arr(i, j, 0) = (enable_ysu_countergradient) ? HGAMT / pblh : zero;
                hgamq_arr(i, j, 0) = (enable_ysu_countergradient && use_moisture) ? HGAMQ / pblh : zero;

                // Compute YSU momentum countergradient terms (hgamu, hgamv)
                // WRF bl_ysu.F90 lines 694-696:
                // brint = -15.9*ust^2/wspd * wstar3/wscale^4
                // hgamu = brint * ux(1), hgamv = brint * vx(1)
                hgamu_arr(i, j, 0) = zero;
                hgamv_arr(i, j, 0) = zero;
                if (SFCFLG && enable_ysu_countergradient) {
                    const Real wspd_sfc = ws10av_arr(i, j, 0);
                    const Real ustar    = u_star_arr(i, j, 0);
                               wscale   = wstar_arr(i, j, 0);   // column-average wscale
                    const Real wstar3   = wstar3_arr(i, j, 0);  // convective velocity scale cubed

                    const Real wscale4 = amrex::max(wscale * wscale * wscale * wscale,
                                                    amrex::Real(1.0e-10));
                    const Real brint = amrex::Real(-15.9) * ustar * ustar
                                      / amrex::max(wspd_sfc, amrex::Real(1.0e-4))
                                      * wstar3 / wscale4;

                    // ux(i,1) and vx(i,1) in WRF = cell-center velocity at lowest level (klo in ERF)
                    const Real u_klo = myhalf * (uvel(i, j, klo) + uvel(i+1, j, klo));
                    const Real v_klo = myhalf * (vvel(i, j, klo) + vvel(i, j+1, klo));
                    hgamu_arr(i, j, 0) = brint * u_klo;
                    hgamv_arr(i, j, 0) = brint * v_klo;
                }

                // VPERT for Pass 3 diagnostic: unnormalized (not divided by pblh)
                // This represents the virtual temperature perturbation at surface
                // When enable_ysu_unbounded_vpert=true, use unbounded VPERT (ERF enhanced)
                // When enable_ysu_unbounded_vpert=false, cap at GAMCRT (WRF-compatible)
                // WRF Reference: module_bl_ysu.F lines 230-240
                // GAP 8: Include WRF height limiter (WRF bl_ysu.F90 line 689-690)
                if (enable_ysu_countergradient) {
                    // Compute height limiter: min(zl1/(sfcfrac*pblh), 1.0)
                    // This reduces VPERT contribution when the first cell is coarse relative to thin PBLs
                    const Real zl1_col = (use_terrain_fitted_coords)
                                       ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                                       : (klo + myhalf) * dz;
                    constexpr Real sfcfrac_h = amrex::Real(0.1);  // WRF SFCFRAC
                    const Real height_lim = amrex::min(zl1_col / (sfcfrac_h * pblh), one_d);

                    const Real VPERT_raw = HGAMT + amrex::Real(0.61) * t_layer * HGAMQ;
                    const Real VPERT_capped = enable_ysu_unbounded_vpert
                                            ? VPERT_raw
                                            : amrex::min(VPERT_raw, GAMCRT);
                    vpert_arr(i, j, 0) = amrex::max(VPERT_capped, zero_d) * height_lim;
                } else {
                    vpert_arr(i, j, 0) = zero;
                }
            }

        });

        // ========================================================================
        // RECOMPUTE rib_enhan_arr WITH UPDATED VPERT
        // ========================================================================
        // BUG FIX: The enhanced Rib array was computed with vpert_arr=0 during the
        // initial precompute phase (line 226). Now that VPERT has been computed and
        // stored in vpert_arr (lines 709-761), we must RECOMPUTE rib_enhan_arr to
        // incorporate the now-nonzero VPERT. This ensures Pass 3 (zero-Ri diagnostic)
        // uses the correctly enhanced Richardson numbers, matching WRF methodology
        // where thermal excess is computed between predictor and corrector passes.
        //
        // Reference: Hong et al. (2006) HND06, WRF module_bl_ysu.F lines 150-170
        //
        BL_PROFILE_VAR("YSUNew_Rib_Recompute", prof_rib_recomp);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
          const amrex::Real t_layer   = t10av_arr(i,j,0);
          const amrex::Real q_frac    = use_moisture ? q10av_arr(i,j,0) : zero;
          const amrex::Real t_layer_v = t_layer * (one + epsv * q_frac);
          // Now vpert_arr contains the computed thermal perturbation from the corrector pass
          const amrex::Real t_enh     = t_layer_v + vpert_arr(i,j,0);
          const amrex::Real z_sfc     = (use_terrain_fitted_coords)
                                      ? Compute_Zrel_AtCellCenter(i,j,klo,z_nd_arr) : zero;
          const amrex::Real zval      = (use_terrain_fitted_coords)
                                      ? Compute_Zrel_AtCellCenter(i,j,k,z_nd_arr)
                                      : (k + myhalf) * dz;
          const amrex::Real zrel      = amrex::max(zval - z_sfc, amrex::Real(1.0e-4));
          const amrex::Real theta_v   = (use_moisture && enable_ysu_liquid_theta)
                                      ? GetThetavl(i,j,k,cell_data,moisture_indices)
                                      : GetThetav(i,j,k,cell_data,moisture_indices);
          const amrex::Real theta_v_klo = (use_moisture && enable_ysu_liquid_theta)
                                        ? GetThetavl(i,j,klo,cell_data,moisture_indices)
                                        : GetThetav(i,j,klo,cell_data,moisture_indices);
          const amrex::Real ws2_raw   = fourth * ((uvel(i,j,k)+uvel(i+1,j,k))*(uvel(i,j,k)+uvel(i+1,j,k))
                                      + (vvel(i,j,k)+vvel(i,j+1,k))*(vvel(i,j,k)+vvel(i,j+1,k)));
          const amrex::Real ws2       = amrex::max(ws2_raw, amrex::Real(1.0));
          // Recompute ONLY rib_enhan_arr with updated vpert_arr; leave rib_base_arr unchanged
          rib_enhan_arr(i,j,k) = CONST_GRAV * zrel * (theta_v - t_enh) / (ws2 * theta_v_klo);
        });
        BL_PROFILE_VAR_STOP(prof_rib_recomp);

        // ========================================================================
        // RE-RUN CORRECTOR PASS 2 WITH UPDATED rib_enhan_arr
        // ========================================================================
        // Now that rib_enhan_arr has been recomputed with the nonzero VPERT values,
        // we need to re-scan for the PBL height using the corrector Ribcr criterion.
        // This corrects the earlier Pass 2 which used zero-VPERT Rib values.
        //
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            // Determine surface-type-dependent critical Richardson number (same as before)
            Real Ribcr;
            {
                bool over_land = (!lmask_arr) || (lmask_arr(i, j, 0) == 1);
                if (over_land) {
                    Ribcr = turbChoice.pbl_ysu_land_Ribcr;
                } else {
                    const Real z0 = z0_arr(i, j, 0);
                    const Real ws_layer = ws10av_arr(i, j, 0);
                    const Real Rossby = ws_layer / (turbChoice.pbl_ysu_coriolis_freq * amrex::max(z0, Real(1.0e-4)));
                    Ribcr = amrex::min(Real(0.16) * std::pow(Real(1.0e-7) * Rossby, -Real(0.18)), Real(0.3));
                }
            }

            // Rescan rib_enhan_arr with updated values
            int kpbl = klo;
            Real zval0 = zero, Rib0 = zero;
            Real Rib = rib_enhan_arr(i,j,klo);
            zval0 = (use_terrain_fitted_coords)
                  ? Compute_Zrel_AtCellCenter(i,j,klo,z_nd_arr)
                  : (klo + myhalf) * dz;
            Rib0 = Rib;
            bool above_critical = (Rib >= Ribcr);

            for (int kk = klo+1; !above_critical && kk <= khi; ++kk) {
                if (rib_enhan_arr(i,j,kk) >= Ribcr) { kpbl = kk; above_critical = true; break; }
                zval0 = (use_terrain_fitted_coords)
                      ? Compute_Zrel_AtCellCenter(i,j,kk,z_nd_arr)
                      : (kk + myhalf) * dz;
                Rib0 = rib_enhan_arr(i,j,kk);
                kpbl = kk;
            }

            // Compute bounds
            const Real z_sfc = (use_terrain_fitted_coords)
                              ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                              : zero;
            const Real dz_terrain = (use_terrain_fitted_coords)
                                   ? (Compute_Zrel_AtCellCenter(i, j, klo + 1, z_nd_arr) - z_sfc)
                                   : dz;
            const Real z_max = (use_terrain_fitted_coords)
                             ? Compute_Zrel_AtCellCenter(i, j, khi, z_nd_arr)
                             : (khi + myhalf) * dz;
            const Real pblh_max = Real(0.9) * z_max;

            amrex::Real pblh_min;
            if (turbChoice.enable_ysu_terrain_pblh_floor && use_terrain_fitted_coords) {
                const amrex::Real met_h_klo = Compute_h_zeta_AtCellCenter(i, j, klo, dxInv, z_nd_arr);
                const amrex::Real dz0 = met_h_klo / dz_inv;
                const amrex::Real z_sfc_col = (use_terrain_fitted_coords)
                                             ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                                             : zero;
                pblh_min = amrex::max(z_sfc_col + Real(0.5)*dz0, Real(10.0));
            } else {
                pblh_min = amrex::max(z_sfc + Real(0.5) * dz_terrain, Real(10.0));
            }

            // Update PBL height with the enhanced Rib values
            if (kpbl < khi && rib_enhan_arr(i,j,kpbl) >= Ribcr) {
                const Real zval = (use_terrain_fitted_coords)
                                ? Compute_Zrel_AtCellCenter(i,j,kpbl,z_nd_arr)
                                : (kpbl + myhalf) * dz;
                Rib = rib_enhan_arr(i,j,kpbl);
                Real pblh_interp = zval0 + (zval - zval0) / (Rib - Rib0) * (Ribcr - Rib0);
                pblh_corr_arr(i, j, 0) = amrex::max(amrex::min(pblh_interp, pblh_max), pblh_min);
            } else {
                pblh_corr_arr(i, j, 0) = pblh_min;
            }
            pbli_arr(i, j, 0) = kpbl;
        });

        // ========================================================================
        // Extension scan using liquid potential temperature (WRF bl_ysu.F90 lines 733-769)
        // ========================================================================
        // Scan upward from Pass 2 PBLH using theta_li to extend kpbl for
        // cloud-topped boundary layers. Triggered when enable_ysu_topdown and
        // use_moisture are both true.
        //
        if (enable_ysu_topdown && use_moisture) {
            ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
            {
                const int kpblold = pbli_arr(i,j,0);  // Starting PBLH index from Pass 2
                bool definebrup = false;

                // Surface liquid-theta virtual potential temperature reference
                const Real thermalli = GetThetavl(i, j, klo, cell_data, moisture_indices);

                // Upward scan using liquid-theta stability criterion
                for (int kk = kpblold; kk < khi; ++kk) {
                    // Velocity shear for Richardson number calculation
                    const Real ws2_raw = fourth * ((uvel(i,j,kk)+uvel(i+1,j,kk))*(uvel(i,j,kk)+uvel(i+1,j,kk))
                                                 + (vvel(i,j,kk)+vvel(i,j+1,kk))*(vvel(i,j,kk)+vvel(i,j+1,kk)));
                    const Real ws2 = amrex::max(ws2_raw, amrex::Real(1.0));

                    // Elevation at level kk
                    const Real z_sfc = (use_terrain_fitted_coords)
                                     ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                                     : zero;
                    const Real zval_kk = (use_terrain_fitted_coords)
                                       ? Compute_Zrel_AtCellCenter(i, j, kk, z_nd_arr)
                                       : (kk + myhalf) * dz;
                    const Real zrel_kk = amrex::max(zval_kk - z_sfc, amrex::Real(1.0e-4));

                    // Liquid-theta virtual potential temperature at current level and surface
                    const Real thlix_kk = GetThetavl(i, j, kk, cell_data, moisture_indices);
                    const Real thlix_klo = GetThetavl(i, j, klo, cell_data, moisture_indices);

                    // Bulk Richardson number using liquid-theta potential temperature
                    const Real bruptmp = CONST_GRAV * zrel_kk * (thlix_kk - thermalli) / (ws2 * thlix_klo);

                    // Stability check: threshold is zero for extension scan
                    const bool stable = (bruptmp >= zero);

                    if (definebrup) {
                        pbli_arr(i,j,0) = kk;
                        definebrup = false;
                    }

                    if (!stable) {
                        // Continue scanning while Richardson number indicates mixing
                        definebrup = true;
                    }
                }

                // Bound result to domain top
                pbli_arr(i,j,0) = amrex::min(pbli_arr(i,j,0), izmax);
            });
        }

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
            // GAP 9: Seed from surface bulk Richardson number (WRF bl_ysu.F90 lines 620-621)
            // WRF seeds the loop with br (surface bulk Ri) so klo is checked first.
            int kpbl_zero = klo;
            Real Rib = rib_enhan_arr(i,j,klo);  // Seed from surface level
            bool above_critical = (Rib >= Ribcr_zero);  // Check if already above critical (Ribcr=0)

            // Scan rib_enhan_arr(i,j,klo..khi) for Ribcr=0 crossing (no interpolation needed)
            for (int kk = klo+1; !above_critical && kk <= khi; ++kk) {
                if (rib_enhan_arr(i,j,kk) >= Ribcr_zero) { kpbl_zero = kk; above_critical = true; break; }
                kpbl_zero = kk;
            }
            pbli_zero_arr(i, j, 0) = kpbl_zero;
        });
        BL_PROFILE_VAR_STOP(prof_pblh);

        // Compute entrainment diffusivity at PBL top cell (WRF bl_ysu.F90 lines 831-925)
        // Uses WRF's buoyancy-flux-based entrainment velocity (we), which accounts for
        // the actual inversion strength (dthvx) at the PBL top.
        {
            const bool enable_ysu_entrainment = turbChoice.enable_ysu_entrainment;
            const bool enable_ysu_cloud_pblh = turbChoice.enable_ysu_cloud_pblh;

            ParallelFor(xybx, [=, zero_d=zero] AMREX_GPU_DEVICE(int i, int j, int) noexcept
            {
                entr_arr(i,j,0) = zero;
                if (!enable_ysu_entrainment) return;

                const int kpbl     = pbli_arr(i,j,0);
                const Real rho_kpbl = cell_data(i, j, kpbl, Rho_comp);
                const Real pblh     = pblh_corr_arr(i,j,0);
                const Real wscale   = wstar_arr(i,j,0);
                const Real ustar    = u_star_arr(i,j,0);
                const Real ustar3   = ustar * ustar * ustar;

                // WRF bl_ysu.F90 lines 831-839: buoyancy-flux-based entrainment velocity
                // wm3 = wstar3 + 5*ust3  (combined convective + mechanical velocity scale)
                const Real wstar3_col = wstar3_arr(i,j,0);
                const Real wm3   = wstar3_col + amrex::Real(5.0) * ustar3;
                const Real wm2   = std::pow(wm3, amrex::Real(2.0) / amrex::Real(3.0));

                // bfxpbl = -0.15 * thvx(klo)/g * wm3 / hpbl
                // thvx at klo: use liquid-theta equivalent if enabled (WRF thlix path)
                const Real thvx_klo = (use_moisture && enable_ysu_liquid_theta)
                                    ? GetThetavl(i, j, klo, cell_data, moisture_indices)
                                    : GetThetav(i, j, klo, cell_data, moisture_indices);
                const Real bfxpbl = amrex::Real(-0.15) * thvx_klo / CONST_GRAV * wm3 / pblh;

                // dthvx = max(thvx(kpbl+1) - thvx(kpbl), tmin)  [inversion strength at PBL top]
                // Guard against kpbl+1 exceeding domain top
                const int kpbl_p1 = amrex::min(kpbl + 1, izmax);
                const Real thvx_kpbl   = (use_moisture && enable_ysu_liquid_theta)
                                       ? GetThetavl(i, j, kpbl,   cell_data, moisture_indices)
                                       : GetThetav(i, j, kpbl,   cell_data, moisture_indices);
                const Real thvx_kpbl_p1 = (use_moisture && enable_ysu_liquid_theta)
                                        ? GetThetavl(i, j, kpbl_p1, cell_data, moisture_indices)
                                        : GetThetav(i, j, kpbl_p1, cell_data, moisture_indices);
                const Real dthvx = amrex::max(thvx_kpbl_p1 - thvx_kpbl, amrex::Real(1.0e-2));

                // we = max(bfxpbl / dthvx, -sqrt(wm2))
                // we is negative for true (downward) entrainment
                const Real we = amrex::max(bfxpbl / dthvx, -std::sqrt(wm2));

                // dz at kpbl level (terrain-adjusted)
                const Real met_h  = (use_terrain_fitted_coords)
                                   ? Compute_h_zeta_AtCellCenter(i, j, kpbl, dxInv, z_nd_arr) : one;
                const Real dz_kpbl = met_h / dz_inv;

                // Apply entrainment only when we < 0 (downward entrainment occurring)
                const Real K_cap      = Real(5.0) * rho_kpbl * wscale * KAPPA * pblh * Real(0.01);

                // Cloudy entrainment correction (WRF bl_ysu.F90 lines 840-875)
                // Applied when qc+qi exceeds threshold at cell kpbl-1 (layer below PBL top)
                // Uses cloud liquid water content and buoyancy jump for entrainment efficiency.
                // Radiation coupling is not required; wstar3_2 and hgamt2 remain zero.
                Real we_final = we;
                Real K_entr_final = (we < zero) ? rho_kpbl * (-we) * dz_kpbl : zero;

                if (enable_ysu_topdown && use_moisture) {
                    const int k_below = amrex::max(kpbl - 1, klo);
                    const Real rho_k = cell_data(i, j, k_below, Rho_comp);
                    Real qc_below = (moisture_indices.qc >= 0)
                                    ? cell_data(i, j, k_below, moisture_indices.qc) / rho_k : zero;
                    Real qi_below = (moisture_indices.qi >= 0)
                                    ? cell_data(i, j, k_below, moisture_indices.qi) / rho_k : zero;

                    constexpr Real cloud_thresh_entr = amrex::Real(0.01e-3);  // 0.01 g/kg

                    if ((qc_below + qi_below) > cloud_thresh_entr && kpbl >= klo + 2) {
                        // Cloud correction path: liquid water at cell k_below triggers adjustment

                        // Entrainment efficiency computed from liquid water content and
                        // buoyancy jump between layers. Limited to 0.4 by WRF methodology.
                        // Uses theta-li jump from k_below to kpbl+1 for inversion strength.

                        // Liquid-theta virtual potential temperature at relevant levels
                        const int k_p2  = amrex::min(kpbl + 1, izmax);
                        const Real thlix_kbelow = GetThetavl(i, j, k_below, cell_data, moisture_indices);
                        const Real thlix_kp2    = GetThetavl(i, j, k_p2,   cell_data, moisture_indices);

                        // Buoyancy jump using liquid-theta across two layers
                        const Real dthvx_li = amrex::max(thlix_kp2 - thlix_kbelow, amrex::Real(0.1));

                        // Entrainment efficiency from cloud liquid water
                        // ent_eff = min(0.2 * 8 * (xlv/cp) * qc_below / dthvx_li, 0.4)
                        constexpr Real xlv_over_cp = amrex::Real(2.5e6) / amrex::Real(1004.0);
                        const Real ent_eff = amrex::min(amrex::Real(0.2) * amrex::Real(8.0)
                                                        * xlv_over_cp * qc_below / dthvx_li,
                                                        amrex::Real(0.4));

                        // Entrainment velocity using cloudy inversion strength
                        // Replaces clear-sky dthvx with dthvx_li in we calculation
                        const Real we_cloud = amrex::max(bfxpbl / dthvx_li, -std::sqrt(wm2));

                        // Cloud-top entrainment velocity contribution
                        // bfxpbl_cloud = -ent_eff * max(sflux, 0) (surface buoyancy flux only, no radsum)
                        const Real bfx0_cloud = amrex::max(sflux_arr(i, j, 0), zero_d);
                        const Real bfxpbl_cloud = -ent_eff * bfx0_cloud;
                        const Real we_cloud_top = amrex::max(bfxpbl_cloud / dthvx_li,
                                                             -std::sqrt(std::pow(wm3, amrex::Real(2.0)/amrex::Real(3.0))));

                        // Summed entrainment velocity
                        we_final = we_cloud + we_cloud_top;

                        // Recompute entrainment coefficient with adjusted entrainment velocity
                        K_entr_final = (we_final < zero) ? rho_kpbl * (-we_final) * dz_kpbl : zero;

                        // Note: WRF zeroes xkzm at k_below when cloud is detected.
                        // In ERF, entrainment is applied only at kpbl level, making this implicit.
                    }
                }

                entr_arr(i,j,0) = amrex::min(K_entr_final, K_cap);

                // Cloud PBLH extension: check for cloud at PBL top cell and extend if present
                if (enable_ysu_cloud_pblh && use_moisture) {
                    const int kpbl_current = pbli_arr(i, j, 0);

                    // Check for cloud condensate at PBL top cell
                    amrex::Real qc_kpbl = zero, qi_kpbl = zero;
                    if (moisture_indices.qc >= 0)
                        qc_kpbl = cell_data(i, j, kpbl_current, moisture_indices.qc) / cell_data(i, j, kpbl_current, Rho_comp);
                    if (moisture_indices.qi >= 0)
                        qi_kpbl = cell_data(i, j, kpbl_current, moisture_indices.qi) / cell_data(i, j, kpbl_current, Rho_comp);

                    // Extend PBL top by one cell if cloud present exceeds threshold
                    if ((qc_kpbl + qi_kpbl) > ysu_qcloud_threshold) {
                        pbli_arr(i, j, 0) = amrex::min(kpbl_current + 1, izmax);
                    }
                }
            });
        }

        // Compute K_down for top-down cloud-driven mixing (H10 Eq. 11)
        FArrayBox K_down_fab(gbx, 1, The_Async_Arena());
        K_down_fab.setVal<RunOn::Device>(zero);
        const auto& K_down_arr = K_down_fab.array();

        if (enable_ysu_topdown) {
            ParallelFor(gbx, [=, zero_d=zero] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                K_down_arr(i,j,k) = zero;
                if (k >= pbli_arr(i,j,0)) return;  // only inside PBL

                const amrex::Real rho   = cell_data(i,j,k,Rho_comp);
                const amrex::Real pblh  = pblh_corr_arr(i,j,0);
                const amrex::Real zval  = (use_terrain_fitted_coords)
                                        ? Compute_Zrel_AtCellCenter(i,j,k,z_nd_arr)
                                        : (k + myhalf) * dz;
                const amrex::Real z_sfc = (use_terrain_fitted_coords)
                                        ? Compute_Zrel_AtCellCenter(i,j,klo,z_nd_arr) : zero;
                const amrex::Real zrel  = zval - z_sfc;
                const amrex::Real pblh_rel = amrex::max(pblh - z_sfc, amrex::Real(1.0e-4));

                // Top-down K-profile (H10 Eq. 11): K_down = rho * wstar_down * kappa * (h-z) * (z/h)^2
                // Use stored top-down convective velocity scale (computed from LRAD in the kernel above)
                // wstar3_down_arr(i,j,0) = 0 when LRAD=0 (radiation not coupled), so K_down will be zero
                const amrex::Real wstar3_down_col = wstar3_down_arr(i,j,0);
                const amrex::Real wstar_down_eff  = (wstar3_down_col > zero)
                                                   ? std::cbrt(wstar3_down_col)
                                                   : zero;
                const amrex::Real zfac_up = amrex::max(zrel / pblh_rel, zero_d);
                K_down_arr(i,j,k) = rho * wstar_down_eff * KAPPA
                                    * amrex::max(pblh_rel - zrel, zero_d)
                                    * zfac_up * zfac_up;
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


        BL_PROFILE_VAR("YSUNew_Kprofile", prof_kprof);
        ParallelFor(gbx, [=, wstar3_arr_cap=wstar3_arr, zol1_arr_cap=zol1_arr, sfcflg_arr_cap=sfcflg_arr,
                          zero_d=zero, one_d=one, two_d=two] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real obuk_val = l_obuk_arr(i, j, 0);
            if (std::abs(obuk_val) < amrex::Real(1.0e-10)) {
                obuk_val = (obuk_val >= zero) ? amrex::Real(1.0e-10) : amrex::Real(-1.0e-10);
            }

            const Real zval = (use_terrain_fitted_coords)
                            ? Compute_Zrel_AtCellCenter(i, j, k, z_nd_arr)
                            : (k + myhalf) * dz;
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
                // When stable (BR > 0, i.e., sflux < 0), skip nonlocal PBL mixing
                // and use free-atmosphere Richardson scheme throughout PBL.
                // Use sflux-based sfcflg for unified definition (Fix 3)
                bool SFCFLG = (sfcflg_arr_cap(i, j, 0) > zero);  // TRUE when unstable/neutral (sflux > 0)

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
                    Real reduction_factor = one - Real(0.15) * amrex::min(total_qcloud / qc_threshold, one_d);
                    // Reduce the phiM enhancement in stable conditions
                    phiM_cloud = one + Real(5.0) * reduction_factor * sf * pblh_corr_arr(i, j, 0) / obuk_val;
                    phit_cloud = one + Real(5.0) * reduction_factor * sf * pblh_corr_arr(i, j, 0) / obuk_val;
                } else if (has_cloud && obuk_val <= zero) {
                    // Unstable layers with clouds: slightly enhance instability
                    // Clouds warm the boundary layer through latent heat release
                    // Boost factor: up to 5% where cloud water exceeds threshold
                    Real cloud_boost = Real(1.0) + Real(0.05) * amrex::min(total_qcloud / qc_threshold, one_d);
                    // Ensure numerically stable exponentiation: base must be in [0.01, 1]
                    // Businger-Dyer form: phi_m = (1 - 16*HOL)^(-1/4) for unstable conditions
                    phiM_cloud = std::pow(
                        amrex::max(one_d - Real(16.0) * HOL_bounded / cloud_boost, Real(0.01)),
                        -one_quarter);
                    phit_cloud = std::pow(
                        amrex::max(one_d - Real(16.0) * HOL_bounded / cloud_boost, Real(0.01)),
                        -one_d / two_d);
                }

                // Use cloud-aware stability functions
                const Real phiM_eff = phiM_cloud;
                const Real phit_eff = phit_cloud;

                // SECTION A: Full Prandtl Number Calculation (WRF bl_ysu.F90 lines 948-971)
                // Three-component formula matching WRF exactly:
                // (a) prfac = bfac*karman*sfcfrac = 6.8*0.4*0.1 = 0.272
                // (b) prfac2 = free-convection correction that reduces Pr for strong convection
                // (c) prnumfac = height-dependent exponential blending
                //
                // Reference: https://github.com/wrf-model/WRF/blob/master/phys/module_bl_ysu.F#L948-L971

                // Step 1: Compute WRF constants (WRF lines 948-950)
                constexpr Real bfac    = amrex::Real(6.8);   // WRF BFAC (not 7.8)
                constexpr Real sfcfrac = amrex::Real(0.1);   // SFCFRAC
                const Real conpr = bfac * KAPPA * sfcfrac;   // = 0.272 (not 0.312)

                // Step 2: Compute height-dependent zq_kp1 for prnumfac calculation
                const Real zq_kp1_prandtl = zval + myhalf * dz_terrain;

                // Step 3: Compute prfac (surface layer Prandtl correction)
                // WRF bl_ysu.F90 line 951: prfac = conpr if SFCFLG, else 0
                const Real prfac = SFCFLG ? conpr : zero;

                // Step 4: Compute prfac2 (free-convection Prandtl correction)
                // WRF bl_ysu.F90 lines 962-963:
                // prfac2 = 15.9*(wstar3+wstar3_2)/ust3 / (1 + 4*karman*(wstar3+wstar3_2)/ust3)
                const Real wstar3_col = wstar3_arr_cap(i, j, 0);
                const Real ust3 = u_star_arr(i, j, 0) * u_star_arr(i, j, 0) * u_star_arr(i, j, 0);
                // Top-down wstar3 term (wstar3_2 in WRF bl_ysu.F90 line 949)
                // Non-zero only when top-down convection is active (requires radiation coupling for LRAD)
                const Real wstar3_2 = wstar3_down_arr(i, j, 0);

                Real prfac2 = zero;
                if (SFCFLG && ust3 > amrex::Real(1.0e-10)) {
                    const Real wstar_tot3 = wstar3_col + wstar3_2;
                    prfac2 = amrex::Real(15.9) * wstar_tot3 / ust3
                           / (one + amrex::Real(4.0) * KAPPA * wstar_tot3 / ust3);
                }

                // Step 5: Compute prnumfac (height-dependent exponential blend)
                // WRF bl_ysu.F90 lines 964-965:
                // prnumfac = -3*(max(zq(k+1) - sfcfrac*hpbl, 0))^2 / hpbl^2
                const Real pblh = pblh_corr_arr(i, j, 0);
                const Real sfclayer = sfcfrac * pblh;
                const Real zdiff = amrex::max(zq_kp1_prandtl - sfclayer, zero_d);
                const Real prnumfac = amrex::Real(-3.0) * zdiff * zdiff / (pblh * pblh);

                // Step 6: Compute prnum0 base Prandtl (WRF lines 966-970)
                // prnum0 = phiH/phiM + prfac
                constexpr Real prmin_wrf = amrex::Real(0.25);  // WRF minimum (different from ERF)
                constexpr Real prmax_wrf = amrex::Real(4.0);
                Real prnum0 = (phit_eff / phiM_eff) + prfac;
                prnum0 = amrex::min(amrex::max(prnum0, prmin_wrf), prmax_wrf);

                // Step 7: Compute heat Prandtl with free-convection correction and height blend
                // WRF uses prnum0_heat = prnum0 / (1 + prfac2*karman*sfcfrac) * exp(prnumfac)
                Real prnum0_heat = prnum0 / (one + prfac2 * KAPPA * sfcfrac);
                prnum0_heat = amrex::min(amrex::max(prnum0_heat, prmin_wrf), prmax_wrf);
                const Real Prt = one + (prnum0_heat - one) * std::exp(prnumfac);

                // Step 8: Compute moisture Prandtl (uses prnum0 WITHOUT the free-convection correction)
                // WRF uses prnum_q = prnum0 * exp(prnumfac)
                const Real prnum_q = one + (prnum0 - one) * std::exp(prnumfac);

                // Use pre-computed wstar from the dedicated recomputation loop (lines 465-545).
                // wstar_arr was computed with pblh_corr_arr to ensure consistency
                // between countergradient diagnostics and K-profile calculations.
                //const Real wstar = wstar_arr(i, j, 0);

                // SFCFLG gating: WRF skips nonlocal mixing in stable PBL (SFCFLG=.FALSE., BR>0, obuk_val>0)
                // In stable conditions, use free-atmosphere Richardson mixing instead.
                // WRF Reference: module_bl_ysu.F lines 200, 220-260
                if (SFCFLG) {
                    // PHASE 13 CHANGES: Use layer-top interface height (zq_kp1) and
                    // level-dependent wscalek for proper K-profile formulation.
                    // WRF Reference: module_bl_ysu.F lines 943-961

                    // Step 1: Compute zq_kp1 (top interface of cell k)
                    // zq(k+1) in WRF = cell center + dz/2
                    const Real zq_kp1 = zval + myhalf * dz_terrain;

                    // Step 2: Get surface height and first-level height (zl1)
                    /*
                    const Real z_sfc = (use_terrain_fitted_coords)
                                     ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                                     : zero;
                    */
                    const Real zl1 = (use_terrain_fitted_coords)
                                   ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                                   : (klo + myhalf) * dz;

                    // Step 3: Compute zfac using zq_kp1 and zl1 (WRF bl_ysu.F90 line 943)
                    // zfac = min(max((1-(zq(k+1)-zl1)/(hpbl-zl1)), zfmin), 1.)
                    constexpr Real zfacmin = amrex::Real(1.0e-8);
                    const Real pblh_rel    = amrex::max(pblh_corr_arr(i, j, 0) - zl1, amrex::Real(1.0e-4));
                    const Real zfac        = amrex::min(
                        amrex::max(one_d - (zq_kp1 - zl1) / pblh_rel, zfacmin), one_d);

                    // Step 4: Compute level-dependent wscalek (for K-profile, NOT for HGAMT/HGAMQ)
                    // WRF stores wstar3 per column; we now have it from Phase 12
                    const Real ust3_wscale = u_star_arr(i, j, 0) * u_star_arr(i, j, 0) * u_star_arr(i, j, 0);

                    Real wscalek_val;
                    // Unstable/neutral: wscalek = (ust3 + phifac*karman*wstar3*(1-zfac))^(1/3)
                    wscalek_val = std::cbrt(ust3_wscale + amrex::Real(8.0) * KAPPA * wstar3_col * (one - zfac));

                    // Step 5: Compute K_m using zq_kp1 and wscalek_val (WRF line 961)
                    // WRF: xkzm(i,k) = wscalek(k)*karman*zq(k+1)*zfac(k)**pfac
                    constexpr Real ckz_pbl = Real(0.001);
                    const Real K_base = ckz_pbl * dz_terrain * rho;
                    constexpr Real pfac = amrex::Real(2.0);
                    K_turb(i, j, k, EddyDiff::Mom_v) = K_base + rho * wscalek_val * KAPPA * zq_kp1 * std::pow(zfac, pfac);

                    // Apply Prandtl number to get heat and moisture diffusivity
                    K_turb(i, j, k, EddyDiff::Theta_v) = K_turb(i, j, k, EddyDiff::Mom_v) / Prt;
                    if (turbChoice.ysu_moistvars) {
                        K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Mom_v) / prnum_q;
                    } else {
                        K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
                    }

                    // Add top-down contribution (H10 Eq. 11)
                    if (enable_ysu_topdown && k < pbli_arr(i, j, 0)) {
                        K_turb(i, j, k, EddyDiff::Mom_v)   += K_down_arr(i, j, k);
                        K_turb(i, j, k, EddyDiff::Theta_v) += K_down_arr(i, j, k);
                    }
                } else {
                    // SECTION B: Stable PBL (sfcflg=false, k < pbli_extent)
                    // WRF uses wscalek = ust/phi_m(zq(k+1)/L), NOT Richardson mixing inside PBL
                    // Richardson mixing is ONLY for k >= pbli_extent (free atmosphere)
                    // WRF Reference: bl_ysu.F90 lines 951-957

                    // Step 1: Compute zq_kp1 (top interface of cell k)
                    const Real zq_kp1_stable = zval + myhalf * dz_terrain;

                    // Step 2: Get first-level height (zl1)
                    const Real zl1_stable = (use_terrain_fitted_coords)
                                          ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                                          : (klo + myhalf) * dz;

                    // Step 3: Compute zfac for stable PBL
                    constexpr Real zfacmin_stable = amrex::Real(1.0e-8);
                    const Real pblh_rel_stable = amrex::max(pblh_corr_arr(i, j, 0) - zl1_stable, amrex::Real(1.0e-4));
                    const Real zfac_stable = amrex::min(
                        amrex::max(one_d - (zq_kp1_stable - zl1_stable) / pblh_rel_stable, zfacmin_stable), one_d);

                    // Step 4: Compute stable wscalek using phi_m
                    // WRF bl_ysu.F90 lines 951-957: wscalek = ust / phi_m(zq(k+1)/L)
                    const Real zol1_stable = zol1_arr_cap(i, j, 0);  // stored from Phase 12
                    const Real zol_ratio = zq_kp1_stable / zl1_stable;  // zq(k+1) / zl1
                    const Real phim_stable_arg = zol1_stable * zol_ratio;  // (z/L) for level k+1
                    const Real phim_stable = one + amrex::Real(5.0) * phim_stable_arg;  // stable: phi_m = 1 + 5*(z/L)
                    const Real wscalek_stable = amrex::max(
                        u_star_arr(i, j, 0) / amrex::max(phim_stable, amrex::Real(0.01)),
                        amrex::Real(0.001));

                    // Step 5: Compute K_m for stable PBL using wscalek
                    constexpr Real ckz_pbl_stable = Real(0.001);
                    const Real K_base_stable = ckz_pbl_stable * dz_terrain * rho;
                    constexpr Real pfac_stable = amrex::Real(2.0);
                    K_turb(i, j, k, EddyDiff::Mom_v) = K_base_stable + rho * wscalek_stable * KAPPA * zq_kp1_stable * std::pow(zfac_stable, pfac_stable);

                    // Step 6: Apply Prandtl number for stable PBL
                    // For stable, prfac=0 (from step 3 of SECTION A)
                    // prnum = phiH/phiM without height blending adjustment
                    const Real prnum_stable = one + (prnum0 - one) * std::exp(prnumfac);
                    K_turb(i, j, k, EddyDiff::Theta_v) = K_turb(i, j, k, EddyDiff::Mom_v) / prnum_stable;

                    if (turbChoice.ysu_moistvars) {
                        K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Mom_v) / prnum_q;
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
                ComputeVerticalDerivativesPBL(i, j, k, uvel, vvel, cell_data, izmin, izmax, pbl_derivative_dz_inv(i,j,k),
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
                const Real dtheta_v_dz = dthetadz;

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

                // GAP 10: In-cloud moist Richardson number modification
                // WRF bl_ysu.F90 lines 992-1000 (imvdif == 1, which is hardcoded in WRF)
                // Only apply when both the current AND next cell are in-cloud (qc+qi > 0.01e-3 kg/kg)
                if (use_moisture && moisture_indices.qc >= 0 && k < izmax) {
                    const Real rho_k    = cell_data(i,j,k,  Rho_comp);
                    const Real rho_kp1  = cell_data(i,j,k+1,Rho_comp);
                    Real qc_k   = (moisture_indices.qc >= 0) ? cell_data(i,j,k,  moisture_indices.qc)/rho_k   : zero;
                    Real qc_kp1 = (moisture_indices.qc >= 0) ? cell_data(i,j,k+1,moisture_indices.qc)/rho_kp1 : zero;
                    Real qi_k   = (moisture_indices.qi >= 0) ? cell_data(i,j,k,  moisture_indices.qi)/rho_k   : zero;
                    Real qi_kp1 = (moisture_indices.qi >= 0) ? cell_data(i,j,k+1,moisture_indices.qi)/rho_kp1 : zero;

                    constexpr Real cloud_thresh = amrex::Real(0.01e-3);  // 0.01 g/kg in kg/kg
                    if ((qc_k + qi_k) > cloud_thresh && (qc_kp1 + qi_kp1) > cloud_thresh) {
                        // WRF bl_ysu.F90 lines 997-1000:
                        // qmean = 0.5*(qvx(k)+qvx(k+1))
                        // tmean = 0.5*(tx(k)+tx(k+1))  [temperature, not theta]
                        // alpha = xlv*qmean/rd/tmean
                        // chi  = xlv*xlv*qmean / (cp*rv*tmean*tmean)
                        // ri = (1+alpha)*(ri - g^2/(ss*tmean*cp) * (chi-alpha)/(1+chi))
                        const Real qv_k_mri   = (moisture_indices.qv >= 0)
                                                ? cell_data(i,j,k,  moisture_indices.qv) / rho_k   : zero;
                        const Real qv_kp1_mri = (moisture_indices.qv >= 0)
                                                ? cell_data(i,j,k+1,moisture_indices.qv) / rho_kp1 : zero;
                        // Use proper temperature (not theta) matching WRF's tx(i,k) — Fix 2
                        const Real T_k   = getTgivenRandRTh(rho_k,   cell_data(i,j,k,  RhoTheta_comp), qv_k_mri);
                        const Real T_kp1 = getTgivenRandRTh(rho_kp1, cell_data(i,j,k+1,RhoTheta_comp), qv_kp1_mri);
                        const Real tmean  = myhalf * (T_k + T_kp1);
                        const Real qmean  = myhalf * (qv_k_mri + qv_kp1_mri);
                        constexpr Real xlv = amrex::Real(2.5e6);   // latent heat of vaporization (J/kg)
                        constexpr Real rd  = amrex::Real(287.0);   // gas constant dry air (J/kg/K)
                        constexpr Real rv  = amrex::Real(461.5);   // gas constant water vapor (J/kg/K)
                        constexpr Real cp  = amrex::Real(1004.0);  // specific heat dry air (J/kg/K)
                        const Real alpha   = xlv * qmean / (rd * tmean);
                        const Real chi     = xlv * xlv * qmean / (cp * rv * tmean * tmean);
                        // Moist Ri correction:
                        grad_Ri = (one + alpha) * (grad_Ri
                                - CONST_GRAV * CONST_GRAV / wind_shear_safe / tmean / cp
                                * (chi - alpha) / (one + chi));
                        grad_Ri = amrex::max(amrex::min(grad_Ri, amrex::Real(100.0)), amrex::Real(-100.0));
                    }
                }

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
                              : 1 - 8 * grad_Ri_safe / (1 + Real(1.746) * std::sqrt(amrex::max(-grad_Ri_safe, zero_d))); // Eqn. A20b/d
                const Real ft = (grad_Ri_safe > 0)
                              ? one / ((one + Real(5.0) * grad_Ri_safe) * (one + Real(5.0) * grad_Ri_safe))
                              : 1 - 8 * grad_Ri_safe / (1 + Real(1.286) * std::sqrt(amrex::max(-grad_Ri_safe, zero_d))); // Eqn. A20a/c
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
                // Q diffusivity matches heat in free atmosphere (same Ri_g-dependent mixing)
                K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
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
            if (turbChoice.pbl_ysunew_highres_bounds) {
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

            #ifdef ERF_USE_WINDFARM
            // Wind farm coupling: if TKE from Fitch/EWP model is available and positive at this cell,
            // augment wscale to reflect enhanced turbulence from turbine wakes.
            // This requires mf_vars_windfarm to be passed in (added to function signature below).
            // For now, add a placeholder comment:
            // TODO: Accept optional MultiFab* windfarm_tke argument. When non-null and
            //       windfarm_tke(i,j,k) > 0, augment K_turb by rho * sqrt(windfarm_tke(i,j,k)) * lscale.
            #endif

            // SECTION C: Asymmetric Background Diffusivity Floor (WRF bl_ysu.F90 lines 141-142, 972-974)
            // Applied INSIDE the PBL only (k < pbli_extent)
            // WRF parameters:
            //   xkzminm = 0.1 m²/s (background for momentum, added to xkzm)
            //   xkzminh = 0.01 m²/s (background for heat and moisture, added to xkzh and xkzq)
            if (k < pbli_extent) {
                constexpr Real xkzminm = amrex::Real(0.1);   // background for momentum
                constexpr Real xkzminh = amrex::Real(0.01);  // background for heat/moisture
                K_turb(i, j, k, EddyDiff::Mom_v)   += rho * xkzminm;
                K_turb(i, j, k, EddyDiff::Theta_v) += rho * xkzminh;
                K_turb(i, j, k, EddyDiff::Q_v)     += rho * xkzminh;
            }
            // NOTE: The free-atmosphere branch (k >= pbli_extent) keeps the existing
            // ckz*dz*rho floor from the WRF high-res bounds parameterization.

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
                K_turb(i, j, k, EddyDiff::HGAMU_v) = hgamu_arr(i, j, 0);
                K_turb(i, j, k, EddyDiff::HGAMV_v) = hgamv_arr(i, j, 0);
            } else {
                // Outside PBL: zero countergradient
                K_turb(i, j, k, EddyDiff::HGAMT_v) = zero;
                K_turb(i, j, k, EddyDiff::HGAMQ_v) = zero;
                K_turb(i, j, k, EddyDiff::HGAMU_v) = zero;
                K_turb(i, j, k, EddyDiff::HGAMV_v) = zero;
            }
        });
        BL_PROFILE_VAR_STOP(prof_kprof);
        amrex::Print()<<" Turbulent Viscosity at cell "<<K_turb(2, 2, 2, EddyDiff::Mom_v)<<" "<<pblh_corr_arr(2, 2, 0)<<std::endl;
        // FOEXTRAP top and bottom ghost cells
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
        {
            K_turb(i, j, klo-1, EddyDiff::Mom_v  ) = K_turb(i, j, klo, EddyDiff::Mom_v  );
            K_turb(i, j, klo-1, EddyDiff::Theta_v) = K_turb(i, j, klo, EddyDiff::Theta_v);
            K_turb(i, j, klo-1, EddyDiff::Q_v    ) = K_turb(i, j, klo, EddyDiff::Q_v    );
            K_turb(i, j, klo-1, EddyDiff::HGAMT_v) = K_turb(i, j, klo, EddyDiff::HGAMT_v);
            K_turb(i, j, klo-1, EddyDiff::HGAMQ_v) = K_turb(i, j, klo, EddyDiff::HGAMQ_v);
            K_turb(i, j, klo-1, EddyDiff::HGAMU_v) = K_turb(i, j, klo, EddyDiff::HGAMU_v);
            K_turb(i, j, klo-1, EddyDiff::HGAMV_v) = K_turb(i, j, klo, EddyDiff::HGAMV_v);
            K_turb(i, j, klo-1, EddyDiff::Turb_lengthscale) = K_turb(i, j, klo, EddyDiff::Turb_lengthscale);
            K_turb(i, j, khi+1, EddyDiff::Mom_v  ) = K_turb(i, j, khi, EddyDiff::Mom_v  );
            K_turb(i, j, khi+1, EddyDiff::Theta_v) = K_turb(i, j, khi, EddyDiff::Theta_v);
            K_turb(i, j, khi+1, EddyDiff::Q_v    ) = K_turb(i, j, khi, EddyDiff::Q_v    );
            K_turb(i, j, khi+1, EddyDiff::HGAMT_v) = K_turb(i, j, khi, EddyDiff::HGAMT_v);
            K_turb(i, j, khi+1, EddyDiff::HGAMQ_v) = K_turb(i, j, khi, EddyDiff::HGAMQ_v);
            K_turb(i, j, khi+1, EddyDiff::HGAMU_v) = K_turb(i, j, khi, EddyDiff::HGAMU_v);
            K_turb(i, j, khi+1, EddyDiff::HGAMV_v) = K_turb(i, j, khi, EddyDiff::HGAMV_v);
            K_turb(i, j, khi+1, EddyDiff::Turb_lengthscale) = K_turb(i, j, khi, EddyDiff::Turb_lengthscale);
        });
    }// mfi
}

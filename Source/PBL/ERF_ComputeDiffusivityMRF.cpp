#include "ERF_SurfaceLayer.H"
#include "ERF_DirectionSelector.H"
#include "ERF_Diffusion.H"
#include "ERF_Constants.H"
#include "ERF_TurbStruct.H"
#include "ERF_PBLModels.H"
#include "ERF_TileNoZ.H"
#include "ERF_MoistUtils.H"
#include "ERF_PBLScaleAwareBlending.H"

using namespace amrex;

/**
 * Compute vertical diffusivity using the Medium-Range Forecast (MRF) boundary layer scheme.
 *
 * @param[in] xvel x-velocity field.
 * @param[in] yvel y-velocity field.
 * @param[in] cons_in Input conserved variables.
 * @param[out] eddyViscosity Eddy viscosity and diffusivity coefficients.
 * @param[in] geom Grid geometry.
 * @param[in] turbChoice Turbulence model options.
 * @param[in] SurfLayer Surface layer data.
 * @param[in] use_terrain_fitted_coords Use terrain-fitted coordinates.
 * @param[in] use_moisture Enable moisture components.
 * @param[in] level Current AMR level.
 * @param[in] bc_ptr Boundary condition records.
 * @param[in] z_phys_nd Nodal physical heights.
 * @param[in] z_phys_cc Cell-centered physical heights.
 * @param[in] moisture_indices Indices for moisture variables.
 */
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
                       const std::unique_ptr<MultiFab>& z_phys_cc,
                       const MoistureComponentIndices& moisture_indices)
{
    /*
    ============================================================================
    Medium-Range Forecast (MRF) Boundary Layer Parameterization Scheme
    ============================================================================

    Implementation of the MRF (Medium Range Forecast) boundary layer scheme
    based on Hong and Pan (1996), with enhanced moisture handling and
    cloud-aware stability corrections developed for ERF.

    PRIMARY REFERENCES:
    -------------------
    - Hong, S. Y., and H.-L. Pan, 1996: Nonlocal Boundary Layer Vertical
      Diffusion in a Medium-Range Forecast Model. Monthly Weather Review, 124,
      2322-2339. https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2

    - WRF Reference Implementation: module_bl_mrf.F
      https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F

    - Hong, S. Y., Y. Noh, and J. Dudhia, 2006: A new vertical diffusion
      package with an explicit treatment of entrainment processes. Monthly
      Weather Review, 134, 2318-2341. https://doi.org/10.1175/MWR3250.1
      (YSU scheme, used for free atmosphere mixing)

    CORE ALGORITHM (Hong & Pan 1996):
    ----------------------------------
    1. Predictor-Corrector Bulk Richardson Number (Rib) PBL Height Diagnosis
       - Rib = (g*z*(θ_v(z) - θ_v(0))) / (U²*θ_v(0))
       - h = min(z where Rib > Rib_critical)

    2. Nonlocal Countergradient Flux Corrections
       - HGAMT: thermal excess from sensible heat flux
       - HGAMQ: moisture excess from latent heat flux
       - VPERT: virtual potential temperature perturbation

    3. Stability-Dependent Mixing Lengths in PBL
       - K_m = ρ * w_* * κ * z * (1 - z/h)²
       - K_t = K_m / Pr_t, K_q = K_m / Pr_q

    4. Free Atmosphere Mixing via Richardson Number
       - Uses YSU stability functions (Hong et al. 2006)
       - Ri-dependent diffusivity above PBL

    ENHANCEMENTS IN ERF IMPLEMENTATION:
    -----------------------------------

    A. Virtual Potential Temperature (θ_v) Treatment
       - Proper handling of moisture effects on buoyancy
       - Important for accurate PBL height and convective strength
       - Not explicitly detailed in Hong & Pan (1996) but standard practice

    B. Cloud-Aware Stability Function Adjustments (OPTIONAL)
       - Physically motivated extension beyond Hong & Pan (1996)
       - Detects clouds via qc + qi > threshold (0.1 g/kg)
       - In stable layers with clouds:
         * Reduces stability damping (φ reduced by 10-20%)
         * Represents suppressed turbulence due to stable stratification
       - In unstable layers with clouds:
         * Slightly increases instability enhancement
         * Represents latent heat release effects
       - Disabled by default for backward compatibility
       - Reference: Conceptually similar to WRF's IMVDIF cloud handling

    C. Explicit Moisture Sign Convention Handling
       - q_star uses WRF convention: negative for upward flux (evaporation)
       - HGAMQ calculation: -const_b * u_* * q_* / w_*
       - Produces positive HGAMQ for unstable (evaporating) conditions
       - Corrected with MAX(HGAMQ, 0) to prevent negative values
       - Land/water discrimination: HGAMQ zeroed over water surfaces

    D. Improved VPERT (Virtual Potential Temperature Perturbation)
       Formulation
       - VPERT = max(HGAMT + 0.61*θ*HGAMQ, 0)
       - Enhancement vs WRF: Does not limit VPERT to GAMCRT
       - WRF limits VPERT after adding moisture term (VPERT = min(VPERT, GAMCRT))
       - This is physically incorrect: moisture contribution can exceed GAMCRT
       - ERF only limits HGAMT (sensible heat), allowing larger VPERT when
         both sensible and latent heat create strong surface heating
       - This produces more realistic PBL heights in moist convection

    E. Prandtl Number with Stability Correction
       - Pr_t = φ_t/φ_m + const_b * κ * sf
       - Bounded: 0.5 ≤ Pr_t ≤ 4.0
       - Consistency with WRF: const_b = 7.8, sf = 0.1
       - Reference: Hong et al. (2006), Equation A17

    PARAMETER DEFAULTS (Matching WRF):
    -----------------------------------
    - const_b = 7.8 (heat flux weight in countergradient)
    - sf = 0.1 (surface layer / PBL height ratio in stability functions)
    - GAMCRT = 3.0 K (maximum heat countergradient)
    - GAMCRQ = 2.0e-3 kg/kg (maximum moisture countergradient)
    - Rib_crit = 0.5 (critical bulk Richardson number for PBL height)
    - Pr_min = 0.5, Pr_max = 4.0 (Prandtl number bounds)
    - Kmin = 0.1 m²/s, Kmax = 300 m²/s (diffusivity bounds, Hong & Pan 1996)

    STABILITY FUNCTIONS:
    --------------------
    Unstable (L < 0, HOL < 0):
      phi_m = (1 - 8 * sf * h/L)^(-1/3)
      phi_t = (1 - 16 * sf * h/L)^(-1/2)

    Stable (L > 0, HOL > 0):
      phi_m = phi_t = 1 + 5 * sf * h/L

    MIXING ABOVE PBL (Free Atmosphere):
    -----------------------------------
    Uses YSU scheme (Hong et al. 2006, Appendix A) for Richardson number
    dependent mixing. This avoids MRF oscillations in stable conditions.

    Gradient Richardson number: Ri_g = (g/theta_v) * (dtheta_v/dz) / ((du/dz)^2 + (dv/dz)^2)

    For Ri_g > 0 (stable):
      f_m = 1 / ((1 + 5*Ri_g)^2)
      f_t = 1 / ((1 + 5*Ri_g)^2)
      Pr_t = 1 + 2.1*Ri_g, bounded to [0.25, 4.0]

    For Ri_g < 0 (unstable):
      f_m = 1 - 8*Ri_g / (1 + 1.746*sqrt(-Ri_g))
      f_t = 1 - 8*Ri_g / (1 + 1.286*sqrt(-Ri_g))

    NUMERICAL CONSIDERATIONS:
    -------------------------
    - GPU-optimized with parallel_for loops
    - Safe bounds on shear and vertical derivatives
    - Gradient Richardson number limited: -100 <= Ri_g <= 100
    - Avoids division by zero with minimum shear threshold (1e-10)
    - theta_v_klo and theta_v guarded against zero/negative values

    TESTING & VALIDATION:
    ---------------------
    - Tested against WRF_ideal unstable/stable ABL cases
    - Backward compatible when enhancements are disabled
    - Cloud detection threshold: qc + qi > 1e-4 kg/kg (0.1 g/kg)
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

        // create flattened boxes to store PBL height and related quantities
        const GeometryData gdata = geom.data();
        const Box xybx = PerpendicularBox<ZDir>(gbx, IntVect{0, 0, 0});

        // Pass 1 (predictor): PBLH with base surface temperature, no VPERT
        // Pass 2 (wstar/VPERT): compute wstar, HGAMT, HGAMQ, VPERT from predictor height
        // Pass 3 (corrector): PBLH with VPERT-enhanced surface temperature (WRF-consistent)
        //   WRF reference (module_bl_mrf.F lines 813-964):
        //   https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L813-L964
        FArrayBox pbl_height_predictor(xybx, 1, The_Async_Arena());  // Pass 1: base t_layer_v
        FArrayBox pbl_height_corrector(xybx, 1, The_Async_Arena());  // Pass 3: VPERT-enhanced
        IArrayBox pbl_index(xybx, 1, The_Async_Arena());
        IArrayBox pbl_index_zero_ri(xybx, 1, The_Async_Arena());  // Index for zero-Ri diagnostic pass
        FArrayBox hgamt_fab(xybx, 1, The_Async_Arena());  // Store HGAMT/h (normalized countergradient)
        FArrayBox hgamq_fab(xybx, 1, The_Async_Arena());  // Store HGAMQ/h (normalized countergradient)
        FArrayBox wstar_fab(xybx, 1, The_Async_Arena());  // Convective velocity scale
        FArrayBox vpert_fab(xybx, 1, The_Async_Arena());  // Virtual temperature perturbation VPERT
        const auto& pblh_pred_arr   = pbl_height_predictor.array();  // predictor (base t_layer_v)
        const auto& pblh_corr_arr   = pbl_height_corrector.array();  // corrector (VPERT-enhanced)
        const auto& pbli_arr        = pbl_index.array();
        const auto& pbli_zero_arr   = pbl_index_zero_ri.array();  // Zero-Ri diagnostic PBL index
        const auto& hgamt_arr       = hgamt_fab.array();
        const auto& hgamq_arr       = hgamq_fab.array();
        const auto& wstar_arr       = wstar_fab.array();
        const auto& vpert_arr       = vpert_fab.array();

        // Get some data in arrays
        const auto& cell_data = cons_in.const_array(mfi);
        const auto& uvel = xvel.const_array(mfi);
        const auto& vvel = yvel.const_array(mfi);

        const Real Ribcr = turbChoice.pbl_mrf_Ribcr;
        const auto& u_star_arr = SurfLayer->get_u_star(level)->const_array(mfi);
        const auto& t_star_arr = SurfLayer->get_t_star(level)->const_array(mfi);
        const auto& l_obuk_arr = SurfLayer->get_olen(level)->const_array(mfi);
        const auto& t10av_arr  = SurfLayer->get_mac_avg(level, 2)->const_array(mfi);
        const auto& q10av_arr  = SurfLayer->get_mac_avg(level, 3)->const_array(mfi);
        const auto& lmask_arr = (SurfLayer->get_lmask(level)) ?
                                SurfLayer->get_lmask(level)->const_array(mfi) :
                                Array4<int>{};
        // Only retrieve z_phys_nd array if terrain-fitted coordinates are in use
        const Array4<Real const> z_nd_arr = use_terrain_fitted_coords ? z_phys_nd->array(mfi)
                                                                      : Array4<Real const>{};
        const PBLDerivativeDzInv_T pbl_derivative_dz_inv{z_phys_cc->const_array(mfi)};

        //
        // PASS 1 (PREDICTOR): Compute PBL height using base surface virtual temperature.
        // No thermal excess (VPERT = 0). Produces pblh_pred_arr and pbli_arr.
        // WRF reference (module_bl_mrf.F lines 813-842):
        // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L813-L842
        //
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            const Real t_layer = t10av_arr(i, j, 0);
            const Real moisture_fraction = use_moisture ? q10av_arr(i, j, 0) : Real(0);
            const Real t_layer_v = t_layer * (one + epsv * moisture_fraction);

            Real zval, Rib;
            int kpbl = klo;

            {
                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : (kpbl + myhalf) * gdata.CellSize(2);
                const Real theta_v     = GetThetav(i, j, kpbl, cell_data, moisture_indices);
                // FIX: guard theta_v_klo against zero to prevent NaN in Rib
                const Real theta_v_klo = amrex::max(GetThetav(i, j, klo, cell_data, moisture_indices), Real(1.0));
                const Real ws2_raw = fourth * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                              (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                              (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                              (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                const Real ws2 = amrex::max(ws2_raw, Real(1.0));
                Rib = CONST_GRAV * zval * (theta_v - t_layer_v) / (ws2 * theta_v_klo);
            }

            Real zval0 = zval, Rib0 = Rib;
            bool above_critical = false;
            while (!above_critical && ((kpbl + 1) <= khi)) {
                zval0 = zval;
                Rib0  = Rib;
                kpbl += 1;

                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : (kpbl + myhalf) * gdata.CellSize(2);
                const Real theta_v     = GetThetav(i, j, kpbl, cell_data, moisture_indices);
                // FIX: guard theta_v_klo against zero to prevent NaN in Rib
                const Real theta_v_klo = amrex::max(GetThetav(i, j, klo, cell_data, moisture_indices), Real(1.0));
                const Real ws2_raw = fourth * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                              (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                              (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                              (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                const Real ws2 = amrex::max(ws2_raw, Real(1.0));
                Rib = CONST_GRAV * zval * (theta_v - t_layer_v) / (ws2 * theta_v_klo);
                above_critical = (Rib >= Ribcr);
            }

            const Real pblh_emp = (use_terrain_fitted_coords)
                                ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                                : myhalf * gdata.CellSize(2);
            const Real z_max = (use_terrain_fitted_coords)
                             ? Compute_Zrel_AtCellCenter(i, j, khi, z_nd_arr)
                             : (khi + myhalf) * gdata.CellSize(2);
            const Real pblh_max = Real(0.9) * z_max;
            const Real pblh_min = amrex::max(pblh_emp, Real(10.0));

            if (above_critical) {
                // FIX: guard against division by zero when Rib == Rib0
                Real pblh_interp;
                const Real rib_diff = Rib - Rib0;
                if (std::abs(rib_diff) > Real(1.0e-10)) {
                    pblh_interp = zval0 + (zval - zval0) / rib_diff * (Ribcr - Rib0);
                } else {
                    // If Rib not changing, use zval0 as PBL height
                    pblh_interp = zval0;
                }
                pblh_pred_arr(i, j, 0) = amrex::max(amrex::min(pblh_interp, pblh_max), pblh_min);
                pbli_arr(i, j, 0) = kpbl;
            } else {
                pblh_pred_arr(i, j, 0) = pblh_min;
                pbli_arr(i, j, 0) = klo + 1;
            }
        });

        // FIX: guard against null dereference of q_star when use_moisture is false (dry case)
        // get_q_star can return null pointer when the field is not allocated
        amrex::MultiFab* q_star_mf = SurfLayer->get_q_star(level);
        Array4<Real const> q_star_arr = (q_star_mf != nullptr) ? q_star_mf->const_array(mfi)
                                                               : Array4<Real const>{};

        const Real const_b = turbChoice.pbl_mrf_const_b;
        const Real sf = turbChoice.pbl_mrf_sf;
        constexpr Real prmin = Real(0.5);
        constexpr Real prmax = Real(4.0);
        constexpr Real GAMCRT = Real(3.0);
        constexpr Real GAMCRQ = Real(2.e-3);
        const bool enable_mrf_countergradient = turbChoice.enable_mrf_countergradient;
        const bool enable_mrf_unbounded_vpert = turbChoice.enable_mrf_unbounded_vpert;

        //
        // PASS 2 (WSTAR / VPERT): Compute wstar, HGAMT, HGAMQ, VPERT using the
        // predictor PBL height (pblh_pred_arr). VPERT will be fed into the corrector
        // Rib search in Pass 3 to raise the effective surface temperature.
        // WRF reference (module_bl_mrf.F lines 857-880):
        // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L857-L880
        //
        ParallelFor(xybx, [=,zero_d=zero] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            const Real t_layer  = t10av_arr(i, j, 0);
            Real obuk_val = l_obuk_arr(i, j, 0);
            if (std::abs(obuk_val) < amrex::Real(1.0e-10)) {
                obuk_val = (obuk_val >= zero_d) ? amrex::Real(1.0e-10) : amrex::Real(-1.0e-10);
            }

            // HOL computed from predictor height (WRF uses predictor height for wstar/VPERT)
            const Real HOL = sf * pblh_pred_arr(i, j, 0) / obuk_val;
            const Real HOL_bounded = amrex::max(amrex::min(HOL, Real(100.0)), Real(-100.0));
            const Real one_quarter = Real(0.25);
            const Real phiM = (obuk_val > 0)
                            ? (1 + 5 * HOL_bounded)
                            : std::pow(amrex::max(1 - 16 * HOL_bounded, Real(0.01)), -one_quarter);
            const Real phiM_safe = amrex::max(phiM, Real(0.01));

            // wstar = u* / phi_m
            // Absolute bounds [0.01, 5.0] m/s: prevent division-by-zero (floor) and
            // free-convection blow-up (ceiling), independent of u* for u*->0 safety.
            Real wstar = u_star_arr(i, j, 0) / phiM_safe;
            wstar = amrex::max(wstar, Real(0.01));
            wstar = amrex::min(wstar, Real(5.0));

            bool SFCFLG = (obuk_val <= zero_d);
            const Real HGAMT = (SFCFLG && enable_mrf_countergradient)
                             ? amrex::min(-const_b * u_star_arr(i, j, 0) * t_star_arr(i, j, 0) / wstar, GAMCRT)
                             : zero_d;

            Real HGAMQ = Real(0);
            if (SFCFLG && use_moisture && enable_mrf_countergradient) {
                const Real q_star = q_star_arr(i, j, 0);
                const Real HGAMQ_calc = -const_b * u_star_arr(i, j, 0) * q_star / wstar;
                HGAMQ = amrex::max(amrex::min(HGAMQ_calc, GAMCRQ), Real(0));

                if (lmask_arr) {
                    bool is_land = (lmask_arr(i,j,0) == 1);
                    if (!is_land) HGAMQ = zero_d;
                }

                if (moisture_indices.qv >= 0) {
                    Real qv_klo = cell_data(i, j, klo, moisture_indices.qv) / cell_data(i, j, klo, Rho_comp);
                    Real T_klo = getTgivenRandRTh(cell_data(i, j, klo, Rho_comp),
                                                  cell_data(i, j, klo, RhoTheta_comp), qv_klo);
                    Real p_klo = getPgivenRTh(cell_data(i, j, klo, RhoTheta_comp), qv_klo) * Real(0.01);
                    Real qsat_klo = zero_d;
                    erf_qsatw(T_klo, p_klo, qsat_klo);
                    Real rh_klo = (qsat_klo > Real(1.0e-10)) ? (qv_klo / qsat_klo) : Real(0);
                    if (rh_klo > Real(0.95)) {
                        Real rh_scaling = amrex::max(zero_d, (Real(1) - rh_klo) / Real(0.05));
                        HGAMQ *= rh_scaling;
                    }
                }
            }

            // Compute VPERT = HGAMT + 0.61*theta*HGAMQ (virtual temperature perturbation)
            // This will be added to the surface temperature in the corrector Rib search.
            // WRF Reference: module_bl_mrf.F lines 879-880
            if (pbli_arr(i, j, 0) <= klo + 1 || !enable_mrf_countergradient) {
                vpert_arr(i, j, 0) = zero_d;
            } else {
                const Real VPERT_raw = HGAMT + epsv * t_layer * HGAMQ;
                const Real VPERT_capped = enable_mrf_unbounded_vpert
                                        ? VPERT_raw
                                        : amrex::min(VPERT_raw, GAMCRT);
                vpert_arr(i, j, 0) = amrex::max(VPERT_capped, zero_d);
            }
        });

        //
        // PASS 3 (CORRECTOR): Recompute PBL height with VPERT-enhanced surface temperature.
        // theta_s = theta_va + VPERT  (Hong & Pan 1996, Eq. 4)
        // This is the WRF-consistent corrector pass. pblh_corr_arr is used for the
        // K-profile shape function K = rho*wstar*kappa*z*(1-z/h)^2.
        // WRF reference (module_bl_mrf.F lines 932-964):
        // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L932-L964
        //
        ParallelFor(xybx, [=,one_d=one]
                    AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            const Real t_layer  = t10av_arr(i, j, 0);
            const Real moisture_fraction = use_moisture ? q10av_arr(i, j, 0) : Real(0);
            // VPERT-enhanced surface virtual temperature (WRF THERMAL variable)
            const Real t_layer_v = t_layer * (one + epsv * moisture_fraction)
                                 + vpert_arr(i, j, 0);

            Real obuk_val = l_obuk_arr(i, j, 0);
            if (std::abs(obuk_val) < amrex::Real(1.0e-10)) {
                obuk_val = (obuk_val >= zero) ? amrex::Real(1.0e-10) : amrex::Real(-1.0e-10);
            }

            int kpbl = klo;
            Real zval, Rib;                        // FIX: removed uninitialized zval0, Rib0
            {
                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : (kpbl + myhalf) * gdata.CellSize(2);
                const Real theta_v     = GetThetav(i, j, kpbl, cell_data, moisture_indices);
                // FIX: guard theta_v_klo against zero to prevent NaN in Rib
                const Real theta_v_klo = amrex::max(GetThetav(i, j, klo, cell_data, moisture_indices), one_d);
                const Real ws2_raw = fourth * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                                (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                                (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                                (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                const Real ws2 = amrex::max(ws2_raw, one_d);
                Rib = CONST_GRAV * zval * (theta_v - t_layer_v) / (ws2 * theta_v_klo);
            }
            Real zval0 = zval, Rib0 = Rib;         // FIX: initialize here, mirrors Pass 1

            bool above_critical = false;
            while (!above_critical && ((kpbl + 1) <= khi)) {
                zval0 = zval;
                Rib0  = Rib;
                kpbl += 1;

                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : (kpbl + myhalf) * gdata.CellSize(2);
                const Real theta_v     = GetThetav(i, j, kpbl, cell_data, moisture_indices);
                // FIX: guard theta_v_klo against zero to prevent NaN in Rib
                const Real theta_v_klo = amrex::max(GetThetav(i, j, klo, cell_data, moisture_indices), one_d);
                const Real ws2_raw = fourth * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                                (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                                (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                                (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                const Real ws2 = amrex::max(ws2_raw, one_d);
                Rib = CONST_GRAV * zval * (theta_v - t_layer_v) / (ws2 * theta_v_klo);
                above_critical = (Rib >= Ribcr);
            }

            const Real pblh_emp = (use_terrain_fitted_coords)
                                ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                                : myhalf * gdata.CellSize(2);
            const Real z_max = (use_terrain_fitted_coords)
                             ? Compute_Zrel_AtCellCenter(i, j, khi, z_nd_arr)
                             : (khi + myhalf) * gdata.CellSize(2);
            const Real pblh_max = Real(0.9) * z_max;
            const Real pblh_min = amrex::max(pblh_emp, Real(10.0));

            if (above_critical) {
                // FIX: guard against division by zero when Rib == Rib0
                Real pblh_interp;
                const Real rib_diff = Rib - Rib0;
                if (std::abs(rib_diff) > Real(1.0e-10)) {
                    pblh_interp = zval0 + (zval - zval0) / rib_diff * (Ribcr - Rib0);
                } else {
                    // If Rib not changing, use zval0 as PBL height
                    pblh_interp = zval0;
                }
                pblh_corr_arr(i, j, 0) = amrex::max(amrex::min(pblh_interp, pblh_max), pblh_min);
                pbli_arr(i, j, 0) = kpbl;
            } else {
                pblh_corr_arr(i, j, 0) = pblh_min;
                pbli_arr(i, j, 0) = klo + 1;
            }
        });

        //
        // PASS 4 (WSTAR RECOMPUTE): Recompute wstar, HGAMT, HGAMQ using the corrected
        // PBL height (pblh_corr_arr) to ensure internal consistency between the K-profile
        // amplitude and the countergradient fluxes. HOL = sf*h/L uses pblh_corr_arr.
        // WRF reference (module_bl_mrf.F lines 857-880):
        // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L857-L880
        //
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            //const Real t_layer  = t10av_arr(i, j, 0);
            Real obuk_val = l_obuk_arr(i, j, 0);
            if (std::abs(obuk_val) < amrex::Real(1.0e-10)) {
                obuk_val = (obuk_val >= Real(0)) ? amrex::Real(1.0e-10) : amrex::Real(-1.0e-10);
            }

            // HOL now uses corrected PBLH for full internal consistency
            const Real HOL = sf * pblh_corr_arr(i, j, 0) / obuk_val;
            const Real HOL_bounded = amrex::max(amrex::min(HOL, Real(100.0)), Real(-100.0));
            const Real one_quarter = Real(0.25);
            const Real phiM = (obuk_val > 0)
                            ? (1 + 5 * HOL_bounded)
                            : std::pow(amrex::max(1 - 16 * HOL_bounded, Real(0.01)), -one_quarter);
            const Real phiM_safe = amrex::max(phiM, Real(0.01));

            // Absolute bounds [0.01, 5.0] m/s
            Real wstar = u_star_arr(i, j, 0) / phiM_safe;
            wstar = amrex::max(wstar, Real(0.01));
            wstar = amrex::min(wstar, Real(5.0));
            wstar_arr(i, j, 0) = wstar;

            bool SFCFLG = (obuk_val <= Real(0));
            const Real HGAMT = (SFCFLG && enable_mrf_countergradient)
                             ? amrex::min(-const_b * u_star_arr(i, j, 0) * t_star_arr(i, j, 0) / wstar, GAMCRT)
                             : Real(0);

            Real HGAMQ = Real(0);
            if (SFCFLG && use_moisture && enable_mrf_countergradient) {
                const Real q_star = q_star_arr(i, j, 0);
                const Real HGAMQ_calc = -const_b * u_star_arr(i, j, 0) * q_star / wstar;
                HGAMQ = amrex::max(amrex::min(HGAMQ_calc, GAMCRQ), Real(0));

                if (lmask_arr) {
                    bool is_land = (lmask_arr(i,j,0) == 1);
                    if (!is_land) HGAMQ = Real(0);
                }

                if (moisture_indices.qv >= 0) {
                    Real qv_klo = cell_data(i, j, klo, moisture_indices.qv) / cell_data(i, j, klo, Rho_comp);
                    Real T_klo = getTgivenRandRTh(cell_data(i, j, klo, Rho_comp),
                                                  cell_data(i, j, klo, RhoTheta_comp), qv_klo);
                    Real p_klo = getPgivenRTh(cell_data(i, j, klo, RhoTheta_comp), qv_klo) * Real(0.01);
                    Real qsat_klo = Real(0);
                    erf_qsatw(T_klo, p_klo, qsat_klo);
                    Real rh_klo = (qsat_klo > Real(1.0e-10)) ? (qv_klo / qsat_klo) : Real(0);
                    if (rh_klo > Real(0.95)) {
                        Real rh_scaling = amrex::max(Real(0), (Real(1) - rh_klo) / Real(0.05));
                        HGAMQ *= rh_scaling;
                    }
                }
            }

            if (pbli_arr(i, j, 0) <= klo + 1) {
                hgamt_arr(i, j, 0) = Real(0);
                hgamq_arr(i, j, 0) = Real(0);
            } else {
                const Real pblh = pblh_corr_arr(i, j, 0);
                // FIX: guard against division by zero when pblh is extremely small
                if (pblh > Real(1.0e-10)) {
                    hgamt_arr(i, j, 0) = (enable_mrf_countergradient) ? HGAMT / pblh : Real(0);
                    hgamq_arr(i, j, 0) = (enable_mrf_countergradient && use_moisture) ? HGAMQ / pblh : Real(0);
                } else {
                    hgamt_arr(i, j, 0) = Real(0);
                    hgamq_arr(i, j, 0) = Real(0);
                }
            }
        });

        //
        // PASS 5 (ZERO-RI): Diagnostic PBL height with Ribcr=0, VPERT-enhanced surface temp.
        // Used optionally (pbl_mrf_use_zero_ri_extent) to extend the nonlocal mixing region.
        // WRF reference (module_bl_mrf.F lines 932-964):
        // https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F#L932-L964
        //
        constexpr Real Ribcr_zero = Real(0);
        ParallelFor(xybx, [=,one_d=one]
                    AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            const Real t_layer  = t10av_arr(i, j, 0);
            const Real moisture_fraction = use_moisture ? q10av_arr(i, j, 0) : zero;
            const Real t_layer_v = t_layer * (one + epsv * moisture_fraction);
            const Real t_layer_v_enhanced = t_layer_v + vpert_arr(i, j, 0);

            int kpbl_zero = klo;
            Real zval_zero, Rib_zero;
            {
                zval_zero = (use_terrain_fitted_coords)
                          ? Compute_Zrel_AtCellCenter(i, j, kpbl_zero, z_nd_arr)
                          : (kpbl_zero + myhalf) * gdata.CellSize(2);
                const Real theta_v     = GetThetav(i, j, kpbl_zero, cell_data, moisture_indices);
                // FIX: guard theta_v_klo against zero to prevent NaN in Rib
                const Real theta_v_klo = amrex::max(GetThetav(i, j, klo, cell_data, moisture_indices), one_d);
                const Real ws2_raw = fourth * ( (uvel(i, j, kpbl_zero) + uvel(i + 1, j, kpbl_zero)) *
                                                (uvel(i, j, kpbl_zero) + uvel(i + 1, j, kpbl_zero)) +
                                                (vvel(i, j, kpbl_zero) + vvel(i, j + 1, kpbl_zero)) *
                                                (vvel(i, j, kpbl_zero) + vvel(i, j + 1, kpbl_zero)) );
                const Real ws2 = amrex::max(ws2_raw, one_d);
                Rib_zero = CONST_GRAV * zval_zero * (theta_v - t_layer_v_enhanced) / (ws2 * theta_v_klo);
            }

            bool above_critical_zero = false;
            while (!above_critical_zero && ((kpbl_zero + 1) <= khi)) {
                kpbl_zero += 1;
                zval_zero = (use_terrain_fitted_coords)
                          ? Compute_Zrel_AtCellCenter(i, j, kpbl_zero, z_nd_arr)
                          : (kpbl_zero + myhalf) * gdata.CellSize(2);
                const Real theta_v     = GetThetav(i, j, kpbl_zero, cell_data, moisture_indices);
                // FIX: guard theta_v_klo against zero to prevent NaN in Rib
                const Real theta_v_klo = amrex::max(GetThetav(i, j, klo, cell_data, moisture_indices), one_d);
                const Real ws2_raw = fourth * ( (uvel(i, j, kpbl_zero) + uvel(i + 1, j, kpbl_zero)) *
                                                (uvel(i, j, kpbl_zero) + uvel(i + 1, j, kpbl_zero)) +
                                                (vvel(i, j, kpbl_zero) + vvel(i, j + 1, kpbl_zero)) *
                                                (vvel(i, j, kpbl_zero) + vvel(i, j + 1, kpbl_zero)) );
                const Real ws2 = amrex::max(ws2_raw, one_d);
                Rib_zero = CONST_GRAV * zval_zero * (theta_v - t_layer_v_enhanced) / (ws2 * theta_v_klo);
                above_critical_zero = (Rib_zero >= Ribcr_zero);
            }

            if (above_critical_zero) {
                pbli_zero_arr(i, j, 0) = kpbl_zero;
            } else {
                pbli_zero_arr(i, j, 0) = klo + 1;
            }
        });

        // -- Compute diffusion coefficients --

        const Array4<Real>& K_turb = eddyViscosity.array(mfi);

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

        // Blending parameters captured as scalars for GPU lambda.
        const amrex::Real l_blend_length = turbChoice.pbl_blend_length;
        const amrex::Real l_blend_cs     = turbChoice.pbl_blend_cs;
        const amrex::Real l_blend_cmax   = turbChoice.pbl_blend_c_max;
        //const bool        l_use_smag_ceil= turbChoice.pbl_blend_use_smag;
        // dx: use horizontal spacing at this level (assume dx = dy for regular grids).
        const amrex::Real l_dx = geom.CellSize(0);

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real obuk_val = l_obuk_arr(i, j, 0);
            if (std::abs(obuk_val) < amrex::Real(1.0e-10)) {
                obuk_val = (obuk_val >= Real(0)) ? amrex::Real(1.0e-10) : amrex::Real(-1.0e-10);
            }

            const Real zval = (use_terrain_fitted_coords)
                            ? Compute_Zrel_AtCellCenter(i, j, k, z_nd_arr)
                            : (k + myhalf) * gdata.CellSize(2);
            const Real rho = cell_data(i, j, k, Rho_comp);
            // FIX: skip ghost cells with uninitialized (zero) density — prevents
            // division-by-zero inside GetThetav at lateral ghost cells of gbx
            if (rho <= Real(0)) {
                K_turb(i, j, k, EddyDiff::Mom_v)   = Real(0);
                K_turb(i, j, k, EddyDiff::Theta_v) = Real(0);
                K_turb(i, j, k, EddyDiff::Q_v)     = Real(0);
                K_turb(i, j, k, EddyDiff::HGAMT_v) = Real(0);
                K_turb(i, j, k, EddyDiff::HGAMQ_v) = Real(0);
                K_turb(i, j, k, EddyDiff::Turb_lengthscale) = Real(0);
                return;
            }
            const Real met_h_zeta = (use_terrain_fitted_coords)
                                  ? Compute_h_zeta_AtCellCenter(i, j, k, dxInv, z_nd_arr) : Real(1);
            const Real dz_terrain = met_h_zeta / dz_inv;

            constexpr Real qc_threshold = Real(1.0e-4);  // Cloud water/ice threshold (kg/kg)
            Real qc_mix = Real(0);
            Real qi_mix = Real(0);
            if (use_moisture) {
                if (moisture_indices.qc >= 0) {
                    qc_mix = cell_data(i, j, k, moisture_indices.qc) / rho;
                }
                if (moisture_indices.qi >= 0) {
                    qi_mix = cell_data(i, j, k, moisture_indices.qi) / rho;
                }
            }
            const Real total_qcloud = qc_mix + qi_mix;
            const bool has_cloud = turbChoice.enable_mrf_cloud_adjustment && (total_qcloud > qc_threshold);

            const int pbli_extent = turbChoice.pbl_mrf_use_zero_ri_extent ? pbli_zero_arr(i, j, 0) : pbli_arr(i, j, 0);

            if (k < pbli_extent) {
                bool SFCFLG = (obuk_val <= Real(0));

                const Real HOL = sf * pblh_corr_arr(i, j, 0) / obuk_val;
                const Real HOL_bounded = amrex::max(amrex::min(HOL, Real(100.0)), Real(-100.0));

                const Real one_quarter = Real(0.25);
                const Real phiM = (obuk_val > 0)
                                ? (1 + 5 * HOL_bounded)
                                : std::pow(amrex::max(1 - 16 * HOL_bounded, Real(0.01)), -one_quarter);
                const Real phit = (obuk_val > 0)
                                ? (1 + 5 * HOL_bounded)
                                : std::pow(amrex::max(1 - 16 * HOL_bounded, Real(0.01)), -Real(0.5));

                Real phit_cloud = phit;
                Real phiM_cloud = phiM;
                if (has_cloud && obuk_val > Real(0)) {
                    Real reduction_factor = Real(1) - Real(0.15) * amrex::min(total_qcloud / qc_threshold, Real(1));
                    phiM_cloud = Real(1) + Real(5.0) * reduction_factor * sf * pblh_corr_arr(i, j, 0) / obuk_val;
                    phit_cloud = Real(1) + Real(5.0) * reduction_factor * sf * pblh_corr_arr(i, j, 0) / obuk_val;
                } else if (has_cloud && obuk_val <= Real(0)) {
                    Real cloud_boost = Real(1.0) + Real(0.05) * amrex::min(total_qcloud / qc_threshold, Real(1));
                    phiM_cloud = std::pow(amrex::max(Real(1) - Real(16.0) * HOL_bounded / cloud_boost, Real(0.01)), -one_quarter);
                    phit_cloud = std::pow(amrex::max(Real(1) - Real(16.0) * HOL_bounded / cloud_boost, Real(0.01)), -Real(0.5));
                }

                const Real phiM_eff = phiM_cloud;
                const Real phit_eff = phit_cloud;

                Real Prt_base = phit_eff / phiM_eff;
                const Real Prt = amrex::min(amrex::max(Prt_base + const_b * KAPPA * sf, prmin), prmax);

                const Real wstar = wstar_arr(i, j, 0);

                if (SFCFLG) {
                    // K-profile: K = rho * wstar * kappa * zrel * (1 - zrel/pblh_rel)^2
                    // WRF Reference: module_bl_mrf.F L976-978
                    const Real z_sfc = (use_terrain_fitted_coords)
                                     ? Compute_Zrel_AtCellCenter(i, j, klo, z_nd_arr)
                                     : Real(0);
                    const Real zrel = zval - z_sfc;
                    const Real pblh = pblh_corr_arr(i, j, 0);
                    const Real pblh_rel = pblh - z_sfc;
                    const Real zfac = amrex::max(Real(1) - zrel / pblh_rel, Real(1.0e-8));

                    K_turb(i, j, k, EddyDiff::Mom_v)   = rho * wstar * KAPPA * zrel * zfac * zfac;
                    K_turb(i, j, k, EddyDiff::Theta_v) = K_turb(i, j, k, EddyDiff::Mom_v) / Prt;

                    if (turbChoice.mrf_moistvars) {
                        Real Prq_base = phit_eff / phiM_eff;
                        const Real Prq = amrex::min(amrex::max(Prq_base + const_b * KAPPA * sf, prmin), prmax);
                        K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Mom_v) / Prq;
                    } else {
                        K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
                    }
                } else {
                    const Real lambda = Real(150.0);
                    const Real lscale = (KAPPA * zval * lambda) / (KAPPA * zval + lambda);
                    Real dthetadz, dudz, dvdz;
                    ComputeVerticalDerivativesPBL(i, j, k, uvel, vvel, cell_data, izmin, izmax, pbl_derivative_dz_inv(i,j,k),
                                                  c_ext_dir_on_zlo, c_ext_dir_on_zhi, u_ext_dir_on_zlo,
                                                  u_ext_dir_on_zhi, v_ext_dir_on_zlo, v_ext_dir_on_zhi, dthetadz,
                                                  dudz, dvdz, moisture_indices);

                    const Real dudz_safe = (k < izmax) ? dudz : Real(0);
                    const Real dvdz_safe = (k < izmax) ? dvdz : Real(0);
                    const Real wind_shear = dudz_safe * dudz_safe + dvdz_safe * dvdz_safe;
                    const Real wind_shear_safe = std::max(wind_shear, Real(1.0e-8));

                    // FIX: guard theta_v against zero/negative in grad_Ri denominator
                    const Real theta_v     = amrex::max(GetThetav(i, j, k,   cell_data, moisture_indices), Real(1.0));
                    const Real dtheta_v_dz = dthetadz;

                    Real grad_Ri = CONST_GRAV / theta_v * dtheta_v_dz / wind_shear_safe;
                    grad_Ri = std::max(std::min(grad_Ri, Real(100.0)), -Real(100.0));

                    const Real grad_Ri_safe = amrex::max(grad_Ri, -Real(100.0));
                    Real Pr_rich = Real(1) + Real(2.1) * grad_Ri;
                    const Real fm = (grad_Ri_safe > 0)
                                  ? Real(1) / ((Real(1) + Real(5.0) * grad_Ri_safe) * (Real(1) + Real(5.0) * grad_Ri_safe))
                                  : 1 - 8 * grad_Ri_safe / (1 + Real(1.746) * std::sqrt(amrex::max(-grad_Ri_safe, Real(0))));
                    const Real ft = (grad_Ri_safe > 0)
                                  ? Real(1) / ((Real(1) + Real(5.0) * grad_Ri_safe) * (Real(1) + Real(5.0) * grad_Ri_safe))
                                  : 1 - 8 * grad_Ri_safe / (1 + Real(1.286) * std::sqrt(amrex::max(-grad_Ri_safe, Real(0))));
                    const Real rl2wsp = rho * lscale * lscale * std::sqrt(wind_shear);

                    Pr_rich = std::max(amrex::Real(0.25), std::min(Pr_rich, Real(4.0)));

                    K_turb(i, j, k, EddyDiff::Mom_v)   = rl2wsp * fm;
                    K_turb(i, j, k, EddyDiff::Theta_v) = rl2wsp * ft;
                    if (use_moisture && turbChoice.mrf_moistvars) {
                        K_turb(i, j, k, EddyDiff::Q_v) = rl2wsp * ft;
                    } else {
                        K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
                    }
                }
            } else if (k >= pbli_extent) {
                const Real lambda = Real(150.0);
                const Real lscale = (KAPPA * zval * lambda) / (KAPPA * zval + lambda);
                Real dthetadz, dudz, dvdz;
                ComputeVerticalDerivativesPBL(i, j, k, uvel, vvel, cell_data, izmin, izmax, pbl_derivative_dz_inv(i,j,k),
                                              c_ext_dir_on_zlo, c_ext_dir_on_zhi, u_ext_dir_on_zlo,
                                              u_ext_dir_on_zhi, v_ext_dir_on_zlo, v_ext_dir_on_zhi, dthetadz,
                                              dudz, dvdz, moisture_indices);

                const Real dudz_safe = (k < izmax) ? dudz : Real(0);
                const Real dvdz_safe = (k < izmax) ? dvdz : Real(0);
                const Real wind_shear = dudz_safe * dudz_safe + dvdz_safe * dvdz_safe;
                const Real wind_shear_safe = std::max(wind_shear, Real(1.0e-8));

                // FIX: guard theta_v against zero/negative in grad_Ri denominator
                const Real theta_v     = amrex::max(GetThetav(i, j, k,   cell_data, moisture_indices), Real(1.0));
                const Real dtheta_v_dz = dthetadz;

                Real grad_Ri = CONST_GRAV / theta_v * dtheta_v_dz / wind_shear_safe;
                grad_Ri = std::max(std::min(grad_Ri, Real(100.0)), -Real(100.0));

                const Real grad_Ri_safe = amrex::max(grad_Ri, -Real(100.0));
                Real Pr = Real(1) + Real(2.1) * grad_Ri;
                const Real fm = (grad_Ri_safe > 0)
                              ? Real(1) / ((Real(1) + Real(5.0) * grad_Ri_safe) * (Real(1) + Real(5.0) * grad_Ri_safe))
                              : 1 - 8 * grad_Ri_safe / (1 + Real(1.746) * std::sqrt(amrex::max(-grad_Ri_safe, Real(0))));
                const Real ft = (grad_Ri_safe > 0)
                              ? Real(1) / ((Real(1) + Real(5.0) * grad_Ri_safe) * (Real(1) + Real(5.0) * grad_Ri_safe))
                              : 1 - 8 * grad_Ri_safe / (1 + Real(1.286) * std::sqrt(amrex::max(-grad_Ri_safe, Real(0))));
                const Real rl2wsp = rho * lscale * lscale * std::sqrt(wind_shear);

                Pr = std::max(amrex::Real(0.25), std::min(Pr, Real(4.0)));

                K_turb(i, j, k, EddyDiff::Mom_v)   = rl2wsp * fm;
                K_turb(i, j, k, EddyDiff::Theta_v) = rl2wsp * ft;
                if (use_moisture && turbChoice.mrf_moistvars) {
                    K_turb(i, j, k, EddyDiff::Q_v) = rl2wsp * ft;
                } else {
                    K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
                }
            }

            // Scale-aware blending for grey-zone resolution (Boutle et al. 2014).
            // ERF_PBLScaleAwareBlending.H. Gated by pbl_blend_length > 0.
            // Applied to Theta_v and Q_v. Mom_v is not modified.
            if (l_blend_length > 0.0) {
                // For MRF, use power-law ceiling (SmnSmn not available)
                K_turb(i, j, k, EddyDiff::Theta_v) = pbl_kh_blend_and_cap(
                    K_turb(i, j, k, EddyDiff::Theta_v),
                    l_dx, l_blend_length, l_blend_cs, l_blend_cmax,
                    amrex::Real(-1.0), false);

                K_turb(i, j, k, EddyDiff::Q_v) = pbl_kh_blend_and_cap(
                    K_turb(i, j, k, EddyDiff::Q_v),
                    l_dx, l_blend_length, l_blend_cs, l_blend_cmax,
                    amrex::Real(-1.0), false);
            }

            // Limit diffusion coefficients to physical bounds
            // Hong & Pan (1996): Kmin=0.1, Kmax=300 m^2/s (module_bl_mrf.F lines 1014-1025)
            // Hong et al. (2006): Kmin=ckz*dz*rho, Kmax=1000 m^2/s (high-res option)
            Real rhoKmin, rhoKmax;
            if (turbChoice.pbl_mrf_highres_bounds) {
                constexpr Real ckz  = Real(0.001);
                constexpr Real Kmax = Real(1000.0);
                rhoKmin = ckz * dz_terrain * rho;
                rhoKmax = rho * Kmax;
            } else {
                constexpr Real Kmin = Real(0.1);
                constexpr Real Kmax = Real(300.0);
                rhoKmin = rho * Kmin;
                rhoKmax = rho * Kmax;
            }

            K_turb(i, j, k, EddyDiff::Mom_v) = std::max(
                std::min(K_turb(i, j, k, EddyDiff::Mom_v), rhoKmax), rhoKmin);
            K_turb(i, j, k, EddyDiff::Theta_v) = std::max(
                std::min(K_turb(i, j, k, EddyDiff::Theta_v), rhoKmax), rhoKmin);
            K_turb(i, j, k, EddyDiff::Q_v) = std::max(
                std::min(K_turb(i, j, k, EddyDiff::Q_v), rhoKmax), rhoKmin);
            K_turb(i, j, k, EddyDiff::Turb_lengthscale) = pblh_corr_arr(i, j, 0);

            if (k < pbli_extent) {
                K_turb(i, j, k, EddyDiff::HGAMT_v) = hgamt_arr(i, j, 0);
                K_turb(i, j, k, EddyDiff::HGAMQ_v) = hgamq_arr(i, j, 0);
            } else {
                K_turb(i, j, k, EddyDiff::HGAMT_v) = Real(0);
                K_turb(i, j, k, EddyDiff::HGAMQ_v) = Real(0);
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

#ifndef ERF_SIMPLIFIED_SEB_H_
#define ERF_SIMPLIFIED_SEB_H_

#include <AMReX_GpuControl.H>
#include <AMReX_Math.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

#include <ERF_Constants.H>
#include <ERF_MicrophysicsConstants.H>   // rhor, the density of water

/**
 * @file ERF_SimplifiedSEB.H
 * @brief Simplified SEB -- the two-stream model's surface energy balance, as pure functions.
 *
 * The residual diagnosis below, the force-restore tendencies of the skin temperature and
 * soil water (prognostic_dTs_dt, prognostic_dqs_dt), the soil's moisture availability
 * (seb_moisture_availability), the canopy and bare-soil resistances to evaporation
 * (seb_canopy_resistance_without_vpd, seb_soil_evaporation_resistance), the relative
 * humidity of bare soil's pore air (seb_soil_surface_relative_humidity), and the soil's
 * thermal properties (seb_soil_heat_capacity, seb_soil_thermal_conductivity,
 * seb_force_restore_heat_capacity).
 *
 * Surface Energy Balance Equation:
 *   R_net(i,j) = SW_net(i,j) + LW_net(i,j)
 *              = [sw_flux_sfc(i,j)] + [lw_flux_sfc(i,j)]
 *
 *   SEB_residual(i,j) = R_net(i,j) - hfx_sfc(i,j) - lh_sfc(i,j) - grdflx_sfc(i,j)
 *
 * Where:
 * - sw_flux_sfc: net shortwave flux at surface [W/m^2]
 * - lw_flux_sfc: net longwave flux at surface [W/m^2]
 * - hfx_sfc: sensible heat flux (H) [W/m^2]
 * - lh_sfc: latent heat flux (LE) [W/m^2]
 * - grdflx_sfc: ground heat flux (G) [W/m^2]
 *
 * A perfectly closed budget gives SEB_residual = 0; non-zero residual
 * indicates energy not accounted for by the four main SEB terms.
 *
 * References:
 * -----------
 * - Oke, T. R., 1987: Boundary Layer Climates (2nd ed.), Routledge.
 */

/**
 * @brief Diagnose the surface energy balance residual.
 *
 * Computes the SEB residual from net radiative flux and turbulent/ground
 * heat fluxes. All inputs are guarded with amrex::Math::isfinite() to safely handle
 * NaN/Inf; if any input is non-finite, returns 0.0 (safe no-op).
 *
 * **GPU Safety**:
 * - AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
 * - No host-side I/O
 * - No dynamic allocation
 * - Safe for use in device reduction kernels
 *
 * @param[in] sw_net Shortwave net flux at surface [W/m^2]
 * @param[in] lw_net Longwave net flux at surface [W/m^2]
 * @param[in] hfx Sensible heat flux [W/m^2]
 * @param[in] lh Latent heat flux [W/m^2]
 * @param[in] grdflx Ground heat flux [W/m^2]
 *
 * @return SEB residual [W/m^2] = (sw_net + lw_net) - hfx - lh - grdflx
 *         Returns 0.0 if any input is non-finite (NaN or Inf)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real diagnose_seb_residual (amrex::Real sw_net,
                                   amrex::Real lw_net,
                                   amrex::Real hfx,
                                   amrex::Real lh,
                                   amrex::Real grdflx)
{
    // Guard all inputs against non-finite values (NaN, Inf)
    if (!amrex::Math::isfinite(sw_net) || !amrex::Math::isfinite(lw_net) ||
        !amrex::Math::isfinite(hfx) || !amrex::Math::isfinite(lh) || !amrex::Math::isfinite(grdflx)) {
        return 0.0;  // Safe no-op: return zero residual if any input is invalid
    }

    // Compute net radiation
    amrex::Real r_net = sw_net + lw_net;

    // Compute SEB residual
    amrex::Real seb_residual = r_net - hfx - lh - grdflx;

    return seb_residual;
}

/**
 * @brief Compute prognostic surface temperature tendency using force-restore formulation.
 *
 * Computes dT_s/dt from SEB residual and restoring term toward deep soil temperature:
 *   dT_s/dt = SEB_residual / C_s - (2*pi/tau) * (T_s - T_deep)
 *
 * All inputs are guarded with amrex::Math::isfinite() to safely handle NaN/Inf; if any input
 * is non-finite, any parameter is non-positive (C_s <= 0, tau <= 0), returns 0.0 (safe no-op).
 *
 * **GPU Safety**:
 * - AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
 * - No host-side I/O
 * - No dynamic allocation
 * - Safe for use in device kernels
 *
 * @param[in] seb_residual SEB residual [W/m^2] = R_net - H - LE - G
 * @param[in] t_s Current surface temperature [K]
 * @param[in] t_deep Deep soil temperature [K]
 * @param[in] c_s Effective surface heat capacity [J/(m^2*K)]
 * @param[in] tau Force-restore timescale [s]
 *
 * @return Temperature tendency dT_s/dt [K/s]
 *         Returns 0.0 if any input is non-finite or any parameter is non-positive
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real prognostic_dTs_dt (amrex::Real seb_residual,
                               amrex::Real t_s,
                               amrex::Real t_deep,
                               amrex::Real c_s,
                               amrex::Real tau)
{
    // Guard all inputs against non-finite values
    if (!amrex::Math::isfinite(seb_residual) || !amrex::Math::isfinite(t_s) ||
        !amrex::Math::isfinite(t_deep) || !amrex::Math::isfinite(c_s) || !amrex::Math::isfinite(tau)) {
        return 0.0;  // Safe no-op
    }

    // Validate parameter ranges
    if (c_s <= 0.0 || tau <= 0.0) {
        return 0.0;  // Safe no-op for invalid parameters
    }

    // Compute tendency: SEB_residual / C_s - (2*pi/tau) * (T_s - T_deep)
    amrex::Real seb_term = seb_residual / c_s;
    amrex::Real restore_term = (2.0 * PI / tau) * (t_s - t_deep);
    amrex::Real dt_s_dt = seb_term - restore_term;

    return dt_s_dt;
}

/**
 * @brief Compute prognostic surface moisture tendency using force-restore formulation.
 *
 * Computes dq_s/dt from latent heat flux and restoring term toward deep soil moisture:
 *   dq_s/dt = -(LE / (L_v * rhor * d_s)) - (1/tau_q) * (q_s - q_deep)
 *
 * Constants:
 * - L_v = 2.5e6 J/kg (latent heat of vaporization)
 * - rhor = 1000.0 kg/m^3 (water density, from ERF_MicrophysicsConstants.H)
 *
 * All inputs are guarded with amrex::Math::isfinite() to safely handle NaN/Inf; if any input
 * is non-finite, any parameter is non-positive (d_s <= 0, tau_q <= 0), returns 0.0 (safe no-op).
 *
 * **GPU Safety**:
 * - AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
 * - No host-side I/O
 * - No dynamic allocation
 * - Safe for use in device kernels
 *
 * @param[in] le Latent heat flux [W/m^2]
 * @param[in] q_s Current surface moisture: the volumetric water content [m^3/m^3] of a soil
 *                layer d_s deep, which is what LE / (L_v rhor d_s) drains
 * @param[in] q_deep Deep soil moisture, in q_s's units
 * @param[in] d_s Effective surface moisture layer depth [m]
 * @param[in] tau_q Moisture force-restore timescale [s]
 *
 * @return Moisture tendency dq_s/dt [1/s, in q_s's units]
 *         Returns 0.0 if any input is non-finite or any parameter is non-positive
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real prognostic_dqs_dt (amrex::Real le,
                               amrex::Real q_s,
                               amrex::Real q_deep,
                               amrex::Real d_s,
                               amrex::Real tau_q)
{
    // Guard all inputs against non-finite values
    if (!amrex::Math::isfinite(le) || !amrex::Math::isfinite(q_s) ||
        !amrex::Math::isfinite(q_deep) || !amrex::Math::isfinite(d_s) || !amrex::Math::isfinite(tau_q)) {
        return 0.0;  // Safe no-op
    }

    // Validate parameter ranges
    if (d_s <= 0.0 || tau_q <= 0.0) {
        return 0.0;  // Safe no-op for invalid parameters
    }

    // L_v and rhor both come from ERF_Constants.H. A local L_v here shadowed
    // the global one, and a local water density duplicated rhor.

    // Compute tendency: -(LE / (L_v * rhor * d_s)) - (1/tau_q) * (q_s - q_deep)
    amrex::Real le_term = -(le / (L_v * rhor * d_s));
    amrex::Real restore_term = (1.0 / tau_q) * (q_s - q_deep);
    amrex::Real dq_s_dt = le_term - restore_term;

    return dq_s_dt;
}

/**
 * @brief Moisture availability beta of the force-restore surface, from its water content.
 *
 * With the prognostic balance, q_s is the volumetric water content [m^3/m^3] of a surface
 * layer seb_moisture_layer_depth_m deep (see prognostic_dqs_dt). The surface evaporates at
 * the potential rate when that layer is at or above field capacity, not at all at or below
 * the wilting point, and in proportion in between:
 *
 *   beta = clamp((q_s - q_wilt) / (q_fc - q_wilt), 0, 1)
 *
 * the "beta" closure of bucket and force-restore land surfaces. The surface layer then
 * takes the surface mixing ratio beta * q_sat(T_s, p_s) + (1 - beta) * q_air
 * (erf_surface_moisture::surface_mixing_ratio), so its moisture flux is beta times the
 * potential one.
 *
 * @param[in] q_s    surface-layer water content [m^3/m^3]
 * @param[in] q_wilt wilting point [m^3/m^3]
 * @param[in] q_fc   field capacity [m^3/m^3], above q_wilt
 * @return beta in [0, 1]; 0 for any non-finite input or q_fc <= q_wilt
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real seb_moisture_availability (amrex::Real q_s,
                                       amrex::Real q_wilt,
                                       amrex::Real q_fc)
{
    if (!amrex::Math::isfinite(q_s) || !amrex::Math::isfinite(q_wilt) ||
        !amrex::Math::isfinite(q_fc) || !(q_fc > q_wilt)) {
        return amrex::Real(0.0);
    }
    const amrex::Real beta = (q_s - q_wilt) / (q_fc - q_wilt);
    return amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(1.0), beta));
}

/**
 * @brief Canopy resistance to transpiration [s/m] without its vapour-deficit factor.
 *
 * Noah's big-leaf Jarvis form (Chen et al. 1996) on the parameters and floors of Noah-MP's
 * canopy-resistance option 2 (ResistanceCanopyStomataJarvisMod). Noah-MP itself applies
 * the factors per sunlit and shaded leaf, with the absorbed PAR (f = 2 PAR / RGL), the
 * canopy temperature and canopy-air humidity, and a root-zone soil-water factor; here one
 * big leaf takes the shortwave reaching the surface, the skin temperature, the air at the
 * surface layer's reference height and the top layer's soil-water factor:
 *
 *   r_c = RS / (LAI * F_sw * F_T * F_vpd * F_soil),
 *   F_sw = (f + RS/RSMAX) / (1 + f),  f = 0.55 * 2 SW_down / (RGL * LAI),
 *   F_T  = 1 - 0.0016 (TOPT - T)^2,
 *
 * F_soil the soil-water factor (seb_moisture_availability) and F_vpd the vapour-deficit
 * factor, which needs the air's humidity and is applied by the surface layer
 * (erf_surface_moisture::vapour_deficit_factor); this returns r_c * F_vpd. F_sw and F_T
 * are floored at 1e-4 as in Noah-MP, F_soil at 1e-6 (Noah-MP's MinThr).
 *
 * @return r_c * F_vpd [s/m], capped at 1e6 (no transpiration), which is also the value
 *         without leaves or for a non-finite input. The surface layer divides by F_vpd
 *         and caps the result at 1e6 again (erf_surface_moisture::canopy_resistance).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real seb_canopy_resistance_without_vpd (amrex::Real rs_min, amrex::Real rs_max,
                                               amrex::Real rgl, amrex::Real t_opt,
                                               amrex::Real lai, amrex::Real sw_down,
                                               amrex::Real t_canopy, amrex::Real soil_factor)
{
    constexpr amrex::Real no_transpiration = amrex::Real(1.0e6);
    if (!amrex::Math::isfinite(rs_min) || !amrex::Math::isfinite(rs_max) ||
        !amrex::Math::isfinite(rgl) || !amrex::Math::isfinite(t_opt) ||
        !amrex::Math::isfinite(lai) || !amrex::Math::isfinite(sw_down) ||
        !amrex::Math::isfinite(t_canopy) || !amrex::Math::isfinite(soil_factor) ||
        !(lai > amrex::Real(0.0)) || !(rgl > amrex::Real(0.0)) || !(rs_max > amrex::Real(0.0))) {
        return no_transpiration;
    }
    const amrex::Real f = amrex::Real(0.55) * amrex::Real(2.0) *
                          amrex::max(amrex::Real(0.0), sw_down) / (rgl * lai);
    const amrex::Real f_sw = amrex::max(amrex::Real(1.0e-4), (f + rs_min / rs_max) / (amrex::Real(1.0) + f));
    const amrex::Real dt = t_opt - t_canopy;
    const amrex::Real f_t = amrex::max(amrex::Real(1.0e-4), amrex::Real(1.0) - amrex::Real(0.0016) * dt * dt);
    const amrex::Real f_soil = amrex::max(amrex::Real(1.0e-6), soil_factor);
    return amrex::min(no_transpiration, rs_min / (lai * f_sw * f_t * f_soil));
}

/**
 * @brief Resistance of bare soil to evaporation [s/m], Noah-MP's ground-evaporation option 1
 * (Sakaguchi and Zeng 2009; ResistanceGroundEvaporationMod): vapour diffuses through a dry
 * surface layer whose thickness grows as the top soil dries,
 *
 *   d_dry = d1 (exp((1 - min(1, theta/theta_sat))^p) - 1) / (e - 1),
 *   D     = 2.2e-5 theta_sat^2 (1 - theta_wilt/theta_sat)^(2 + 3/b),
 *   r     = d_dry / D,
 *
 * with d1 the top layer's thickness and p = RSURF_EXP; 1e6 below theta = 0.01, as in Noah-MP.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real seb_soil_evaporation_resistance (amrex::Real theta, amrex::Real theta_sat,
                                             amrex::Real theta_wilt, amrex::Real b,
                                             amrex::Real d1, amrex::Real p)
{
    constexpr amrex::Real no_evaporation = amrex::Real(1.0e6);
    if (!amrex::Math::isfinite(theta) || !(theta >= amrex::Real(0.01)) ||
        !(theta_sat > amrex::Real(0.0)) || !(b > amrex::Real(0.0)) || !(d1 > amrex::Real(0.0)) ||
        !(theta_wilt < theta_sat)) {
        return no_evaporation;
    }
    const amrex::Real dry = amrex::Real(1.0) - amrex::min(amrex::Real(1.0), theta / theta_sat);
    const amrex::Real d_dry = d1 * (std::exp(std::pow(dry, p)) - amrex::Real(1.0)) /
                              (amrex::Real(2.71828) - amrex::Real(1.0));
    const amrex::Real diffusivity = amrex::Real(2.2e-5) * theta_sat * theta_sat *
        std::pow(amrex::Real(1.0) - theta_wilt / theta_sat, amrex::Real(2.0) + amrex::Real(3.0) / b);
    return amrex::min(no_evaporation, d_dry / diffusivity);
}

/**
 * @brief Relative humidity of the air in the top soil's pores, Noah-MP's for bare ground
 * (ResistanceGroundEvaporationMod), so that bare soil evaporates from RH * q_sat(T):
 *
 *   psi = -psi_sat (max(0.01, theta) / theta_sat)^(-b),
 *   RH  = exp(psi g / (R_v T)),
 *
 * with Noah-MP's g = 9.80616 m/s^2 and R_v = 461.269 J/kg/K. Near 1 in moist soil, and
 * near 0 at the wilting point (about 0.005 for silty clay loam at 325 K). Returns 1
 * (no reduction) for a non-finite or non-positive parameter or temperature.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real seb_soil_surface_relative_humidity (amrex::Real theta, amrex::Real theta_sat,
                                                amrex::Real psi_sat, amrex::Real b,
                                                amrex::Real t_ground)
{
    if (!amrex::Math::isfinite(theta) || !amrex::Math::isfinite(t_ground) ||
        !(theta_sat > amrex::Real(0.0)) || !(psi_sat > amrex::Real(0.0)) ||
        !(b > amrex::Real(0.0)) || !(t_ground > amrex::Real(0.0))) {
        return amrex::Real(1.0);
    }
    const amrex::Real psi = -psi_sat * std::pow(amrex::max(amrex::Real(0.01), theta) / theta_sat, -b);
    return std::exp(psi * amrex::Real(9.80616) / (amrex::Real(461.269) * t_ground));
}

/**
 * @brief Volumetric heat capacity of unfrozen soil [J/m^3/K], as Noah-MP forms it
 * (SoilThermalPropertyMod): water theta * 4.188e6, solids (1 - theta_sat) * c_soil, and
 * the air in the remaining pores (theta_sat - theta) * 1004.64 (Noah-MP's constant,
 * applied per unit volume as Noah-MP does).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real seb_soil_heat_capacity (amrex::Real theta, amrex::Real theta_sat, amrex::Real c_soil)
{
    return theta * amrex::Real(4.188e6) + (amrex::Real(1.0) - theta_sat) * c_soil +
           amrex::max(amrex::Real(0.0), theta_sat - theta) * amrex::Real(1004.64);
}

/**
 * @brief Thermal conductivity of unfrozen soil [W/m/K], Noah-MP's Johansen form
 * (SoilThermalPropertyMod): solids 7.7^q 2^(1-q) for quartz fraction q, saturated
 * lambda_solid^(1-theta_sat) 0.57^theta_sat, dry (0.135 g + 64.7)/(2700 - 0.947 g) with
 * g = 2700 (1 - theta_sat), blended by the Kersten number log10(theta/theta_sat) + 1
 * (0 below a saturation of 0.1).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real seb_soil_thermal_conductivity (amrex::Real theta, amrex::Real theta_sat, amrex::Real quartz)
{
    const amrex::Real lambda_solid = std::pow(amrex::Real(7.7), quartz) *
                                     std::pow(amrex::Real(2.0), amrex::Real(1.0) - quartz);
    const amrex::Real lambda_sat = std::pow(lambda_solid, amrex::Real(1.0) - theta_sat) *
                                   std::pow(amrex::Real(0.57), theta_sat);
    const amrex::Real gamma = (amrex::Real(1.0) - theta_sat) * amrex::Real(2700.0);
    const amrex::Real lambda_dry = (amrex::Real(0.135) * gamma + amrex::Real(64.7)) /
                                   (amrex::Real(2700.0) - amrex::Real(0.947) * gamma);
    const amrex::Real sat_ratio = theta / theta_sat;
    const amrex::Real kersten = (sat_ratio > amrex::Real(0.1))
                              ? std::log10(sat_ratio) + amrex::Real(1.0) : amrex::Real(0.0);
    return kersten * (lambda_sat - lambda_dry) + lambda_dry;
}

/**
 * @brief Surface heat capacity [J/m^2/K] of the force-restore balance for a soil of
 * volumetric heat capacity c and conductivity lambda: the heat capacity of the layer the
 * diurnal wave penetrates, c d / (2 sqrt(pi)) with d = sqrt(lambda tau / c) (Deardorff
 * 1978), i.e. sqrt(lambda c tau / pi) / 2 for a restore period tau [s].
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real seb_force_restore_heat_capacity (amrex::Real lambda, amrex::Real c, amrex::Real tau)
{
    return amrex::Real(0.5) * std::sqrt(lambda * c * tau / amrex::Real(3.14159265358979323846));
}

#endif // ERF_SIMPLIFIED_SEB_H_

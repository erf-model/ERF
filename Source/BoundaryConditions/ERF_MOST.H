#ifndef ERF_MOST_H
#define ERF_MOST_H

#include <ERF_Constants.H>
#include <ERF_MOSTRoughness.H>
#include <ERF_Wstar.H>

/**
 * Structure of plain old data relevant to MOST BCs
 */
struct most_data
{
public:
    amrex::Real zref{10.0};           ///< Reference height (m)
    amrex::Real kappa{KAPPA};         ///< von Karman constant
    amrex::Real gravity{CONST_GRAV};  ///< Acceleration due to gravity (m/s^2)
    amrex::Real surf_temp_flux{0.0};  ///< Heat flux TODO: decide whether this is <θ'w'> or <θv'w'> under moist conditions
    amrex::Real surf_moist_flux{0.0}; ///< Moisture flux

    amrex::Real water_depth{30};      ///< Water depth

    static constexpr amrex::Real Bjr_beta = 1.2; // Empirical parameter from Beljaars 1995 QJRMS
};


/**
 * Structure of similarity functions for Moeng formulation
 */
struct similarity_funs
{
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_m (amrex::Real zeta) const
    {
        if (zeta > 0) {
            return -beta_m * zeta;
        } else {
            amrex::Real x = std::sqrt(std::sqrt(1.0 - gamma_m * zeta));
            return 2.0 * std::log(0.5 * (1.0 + x)) + log(0.5 * (1.0 + x * x)) -
                   2.0 * std::atan(x) + PIoTwo;
        }
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_h (amrex::Real zeta) const
    {
        if (zeta > 0) {
            return -beta_h * zeta;
        } else {
            amrex::Real x = std::sqrt(1.0 - gamma_h * zeta);
            return 2.0 * std::log(0.5 * (1.0 + x));
        }
    }

private:
    amrex::Real beta_m{5.0};         ///< Constants from Dyer, BLM, 1974
    amrex::Real beta_h{5.0};         ///< https://doi.org/10.1007/BF00240838
    amrex::Real gamma_m{16.0};
    amrex::Real gamma_h{16.0};
};


/**
 * Solve for MOST fluxes under adiabatic conditions
 *
 * If surface roughness is constant, it isn't necessary to solve for ustar
 * iteratively. For generality, however, iterate_flux should converge in 1 step.
 */
struct adiabatic
{
    adiabatic (amrex::Real zref,
               amrex::Real Tflux,
               amrex::Real Qvflux,
               amrex::Real depth)
    {
        mdata.zref            = zref;
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.water_depth     = depth;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<      amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& /*tvm_arr*/,
                  const amrex::Array4<const amrex::Real>& /*qvm_arr*/,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& /*w_star_arr*/,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& /*t_surf_arr*/,
                  const amrex::Array4<amrex::Real>& /*q_surf_arr*/,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& /*pblh_arr*/,
                  const amrex::Array4<amrex::Real>& Hwave_arr,
                  const amrex::Array4<amrex::Real>& Lwave_arr,
                  const amrex::Array4<amrex::Real>& eta_arr,
                  const RoughnessModel& update_roughness) const
    {
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real ustar = 0.0; // to store previous value
        amrex::Real Hwave = 0.0; // waveheight
        amrex::Real Lwave = 0.0; // wavelength
        amrex::Real eta   = 0.0; // eddy viscosity

        if (Hwave_arr) // wave coupled
        {
            Hwave = Hwave_arr(i,j,k);
            Lwave = Lwave_arr(i,j,k);

            int ie, je;
            ie = i  < lbound(eta_arr).x ? lbound(eta_arr).x : i;
            je = j  < lbound(eta_arr).y ? lbound(eta_arr).y : j;
            ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
            je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;
            eta = eta_arr(ie,je,k,EddyDiff::Mom_v);
        }

        int iter = 0;
        u_star_arr(i,j,k) = mdata.kappa * umm_arr(i,j,k) / std::log(mdata.zref / z0_arr(i,j,k));
        do {
            ustar = u_star_arr(i,j,k);
            update_roughness(z0_arr(i,j,k),
                             ustar, umm, tm_arr(i,j,k), mdata.water_depth,
                             eta, Hwave, Lwave);
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ASSERT_WITH_MESSAGE(iter < max_iters, "Maximum number of MOST iterations reached.");

        t_star_arr(i,j,k) = 0.0;
        q_star_arr(i,j,k) = 0.0;
        olen_arr(i,j,k)   = 1.0e16;
    }

private:
    most_data mdata;

    static constexpr amrex::Real tol{1.0e-5};   // iteration tolerance
    static constexpr amrex::Real WSMIN{0.1};    // minimum wind speed

    // for wave-coupled roughness
    static constexpr amrex::Real eps{1e-15};
    static constexpr amrex::Real z0_eps{1.0e-6};
    static constexpr amrex::Real z0_max{0.1};
};

/**
 * Solve for MOST fluxes with specified surface flux
 *
 * Need to iteratively solve for ustar, which determines the Obukhov length
 * stability parameter; the surface temperature and tstar are derived quantites.
 */
struct surface_flux
{
    surface_flux (amrex::Real zref,
                  amrex::Real Tflux,
                  amrex::Real Qvflux,
                  bool cons_qflux,
                  amrex::Real depth)
    : spec_qflux(cons_qflux)
    {
        mdata.zref            = zref;
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.water_depth     = depth;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<      amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& q_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& Hwave_arr,
                  const amrex::Array4<amrex::Real>& Lwave_arr,
                  const amrex::Array4<amrex::Real>& eta_arr,
                  const RoughnessModel& update_roughness) const
    {
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real ustar = 0.0; // to store previous value
        amrex::Real wstar = 0.0; // temporary scaled var
        amrex::Real tflux = 0.0;
        amrex::Real qflux = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }

        amrex::Real Hwave = 0.0; // waveheight
        amrex::Real Lwave = 0.0; // wavelength
        amrex::Real eta   = 0.0; // eddy viscosity
        if (Hwave_arr) // wave coupled
        {
            Hwave = Hwave_arr(i,j,k);
            Lwave = Lwave_arr(i,j,k);

            int ie, je;
            ie = i  < lbound(eta_arr).x ? lbound(eta_arr).x : i;
            je = j  < lbound(eta_arr).y ? lbound(eta_arr).y : j;
            ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
            je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;
            eta = eta_arr(ie,je,k,EddyDiff::Mom_v);
        }

        int iter = 0;
        do {
            ustar = u_star_arr(i,j,k);

            update_roughness(z0_arr(i,j,k),
                             ustar, umm, tm_arr(i,j,k), mdata.water_depth,
                             eta, Hwave, Lwave);

            // update virtual temperature flux
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h); // <w'Qv'>
            tflux = mdata.surf_temp_flux*(1 + 0.61*qvm_arr(i,j,k)) + qflux*0.61*tm_arr(i,j,k);

            // update w* and Umagmean for Beljaars correction
            if (w_star_arr) {
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }

            // update similarity params
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);

            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(mdata.zref / z0_arr(i,j,k)) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ASSERT_WITH_MESSAGE(iter < max_iters, "Maximum number of MOST iterations reached.");

        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);

        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k); // TODO be consistent with above
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h);
        }
        olen_arr(i,j,k) = Olen;
    }

private:
    most_data mdata;
    similarity_funs sfuns;

    const bool spec_qflux;

    static constexpr amrex::Real tol{1.0e-5};   // iteration tolerance
    static constexpr amrex::Real WSMIN{0.1};    // minimum wind speed

    // for wave-coupled roughness
    static constexpr amrex::Real eps{1e-15};
    static constexpr amrex::Real z0_eps{1.0e-6};
    static constexpr amrex::Real z0_max{0.1};
};

/**
 * Solve for MOST fluxes with specified surface temperature
 *
 * Need to iteratively solve for ustar, which determines the Obukhov length
 * stability parameter; the surface flux and tstar are derived quantities.
 */
struct surface_temp
{
    surface_temp (amrex::Real zref,
                  amrex::Real Tflux,
                  amrex::Real Qvflux,
                  bool cons_qflux,
                  amrex::Real depth)
    : spec_qflux(cons_qflux)
    {
        mdata.zref            = zref;
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.water_depth     = depth;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<      amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& q_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& Hwave_arr,
                  const amrex::Array4<amrex::Real>& Lwave_arr,
                  const amrex::Array4<amrex::Real>& eta_arr,
                  const RoughnessModel& update_roughness) const
    {
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real ustar = 0.0; // to store previous value
        amrex::Real wstar = 0.0; // temporary scaled var
        amrex::Real tflux = 0.0;
        amrex::Real qflux = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real Oleno = 0.0; // old Obukhov length
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            Oleno = Olen;
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }

        amrex::Real Hwave = 0.0; // waveheight
        amrex::Real Lwave = 0.0; // wavelength
        amrex::Real eta   = 0.0; // eddy viscosity
        if (Hwave_arr) // wave coupled
        {
            Hwave = Hwave_arr(i,j,k);
            Lwave = Lwave_arr(i,j,k);

            int ie, je;
            ie = i  < lbound(eta_arr).x ? lbound(eta_arr).x : i;
            je = j  < lbound(eta_arr).y ? lbound(eta_arr).y : j;
            ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
            je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;
            eta = eta_arr(ie,je,k,EddyDiff::Mom_v);
        }

        int iter = 0;
        do {
            ustar = u_star_arr(i,j,k);

            update_roughness(z0_arr(i,j,k),
                             ustar, umm, tm_arr(i,j,k), mdata.water_depth,
                             eta, Hwave, Lwave);

            // update virtual temperature flux
            tflux = -(tm_arr(i,j,k) - t_surf_arr(i,j,k)) * ustar * mdata.kappa /
                     (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h); // <w'T'>
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h); // <w'Qv'>
            tflux = tflux*(1 + 0.61*qvm_arr(i,j,k)) + qflux*0.61*tm_arr(i,j,k);

            // update w* and Umagmean for Beljaars correction
            if (w_star_arr) {
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }

            // update Obukhov length -- prevent excessive sign changes near 0
            Olen = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            bool zero_crossing = (Olen >= 0.0 && Oleno <= 0.0) ||
                                 (Olen <= 0.0 && Oleno >= 0.0);
            bool small_delta_L = (std::fabs(Olen) + std::fabs(Oleno) < 1.);
            if (zero_crossing && small_delta_L) {
                Olen = 0.5 * (Olen + Oleno);
            }
            Oleno = Olen;

            // update similarity params
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);

            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(mdata.zref / z0_arr(i,j,k)) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ASSERT_WITH_MESSAGE(iter < max_iters, "Maximum number of MOST iterations reached.");

        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) /
                            (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h);
        }
        olen_arr(i,j,k) = Olen;
    }

private:
    most_data mdata;
    similarity_funs sfuns;

    const bool spec_qflux;

    static constexpr amrex::Real tol{1.0e-5};   // iteration tolerance
    static constexpr amrex::Real WSMIN{0.1};    // minimum wind speed

    // for wave-coupled roughness
    static constexpr amrex::Real eps{1e-15};
    static constexpr amrex::Real z0_eps{1.0e-6};
    static constexpr amrex::Real z0_max{0.1};
};

#endif


#ifndef ABLMOST_H
#define ABLMOST_H

#include "ERF.H"
using namespace amrex;

/** Monin-Obukhov surface layer profile
 *
 * van der Laan, P., Kelly, M. C., & Sørensen, N. N. (2017). A new k-epsilon
 * model consistent with Monin-Obukhov similarity theory. Wind Energy,
 * 20(3), 479–489. https://doi.org/10.1002/we.2017
 *
 * Consistent with Dyer (1974) formulation from page 57, Chapter 2, Modeling
 * the vertical ABL structure in Modelling of Atmospheric Flow Fields,
 * Demetri P Lalas and Corrado F Ratto, January 1996,
 * https://doi.org/10.1142/2975.
 */
struct ABLMost
{
    enum ThetaCalcType {
        HEAT_FLUX = 0,      ///< Heat-flux specified
        SURFACE_TEMPERATURE ///< Surface temperature specified
    };

    amrex::Real zref{0.2};           ///< Reference height (m)
    amrex::Real z0{0.1};             ///< Roughness height (m)
    amrex::Real utau;                ///< Friction velocity (m/s)
    amrex::Real kappa{KAPPA};        ///< von Karman constant
    amrex::Real gravity{CONST_GRAV}; ///< Acceleration due to gravity (m/s^2)
    amrex::Real obukhov_len{1.0e16}; ///< Non-dimensional Obukhov length

    amrex::Real vel_mean[AMREX_SPACEDIM]; ///< Mean velocity (at zref)
    amrex::Real vmag_mean;                ///< Mean wind speed (at zref)
    amrex::Real theta_mean;               ///< Mean potential temperature

    amrex::Real surf_temp_flux{0.0}; ///< Heat flux
    amrex::Real surf_temp;           ///< Instantaneous surface temperature
    amrex::Real ref_temp;            ///< Reference temperature

    amrex::Real gamma_m{5.0};
    amrex::Real gamma_h{5.0};
    amrex::Real beta_m{16.0};
    amrex::Real beta_h{16.0};

    ThetaCalcType alg_type{HEAT_FLUX};

    void print() const
    {
        amrex::Print() << "ABLMost struct:\n";
        amrex::Print() << "  zref: " << zref << "\n";
        amrex::Print() << "    z0: " << z0 << "\n";
        amrex::Print() << "  utau: " << utau << "\n";
        amrex::Print() << " kappa: " << kappa << "\n";
        amrex::Print() << " gravity: " << gravity << "\n";
        amrex::Print() << " obukhov_len: " << obukhov_len << "\n";
        amrex::Print() << " vel_mean: " << vel_mean[0] << " " << vel_mean[1] << " " << vel_mean[2] << "\n";
        amrex::Print() << " vmag_mean: " << vmag_mean << "\n";
        amrex::Print() << " theta_mean: " << theta_mean << "\n";
        amrex::Print() << " surf_temp_flux: " << surf_temp_flux << "\n";
        amrex::Print() << " surf_temp: " << surf_temp << "\n";
        amrex::Print() << " ref_temp: " << ref_temp << "\n";
        amrex::Print() << " gamma_m: " << gamma_m << "\n";
        amrex::Print() << " gamma_h: " << gamma_h << "\n";
        amrex::Print() << " beta_m: " << beta_m << "\n";
        amrex::Print() << " beta_h: " << beta_h << "\n";
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real phi_m() const
    {
        return std::log(zref / z0) - calc_psi_m(zref / obukhov_len);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real phi_m(amrex::Real z) const
    {
        return std::log(z / z0) - calc_psi_m(z / obukhov_len);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real phi_h() const
    {
        return std::log(zref / z0) - calc_psi_h(zref / obukhov_len);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real phi_h(amrex::Real z) const
    {
        return std::log(z / z0) - calc_psi_h(z / obukhov_len);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_psi_m(amrex::Real zeta) const
    {
        if (zeta > 0) {
            return -gamma_m * zeta;
        } else {
            amrex::Real x = std::sqrt(std::sqrt(1 - beta_m * zeta));
            return 2.0 * std::log(0.5 * (1.0 + x)) + log(0.5 * (1 + x * x)) -
                   2.0 * std::atan(x) + 1.57;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_psi_h(amrex::Real zeta) const
    {
        if (zeta > 0) {
            return -gamma_h * zeta;
        } else {
            amrex::Real x = std::sqrt(1 - beta_h * zeta);
            return 2.0 * std::log(0.5 * (1 + x));
        }
    }

    AMREX_FORCE_INLINE
    void update_fluxes(int max_iters = 25)
    {
        constexpr amrex::Real eps = 1.0e-16;
        amrex::Real zeta = 0.0;
        amrex::Real utau_iter = 0.0;

        // Initialize variables
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        utau = kappa * vmag_mean / (std::log(zref / z0));

        int iter = 0;
        do {
            utau_iter = utau;
            switch (alg_type) {
            case HEAT_FLUX:
                surf_temp = surf_temp_flux * (std::log(zref / z0) - psi_h) /
                                (utau * kappa) + theta_mean;
                break;

            case SURFACE_TEMPERATURE:
                surf_temp_flux = -(theta_mean - surf_temp) * utau * kappa /
                                (std::log(zref / z0) - psi_h);
                break;
            }

            if (std::abs(surf_temp_flux) > eps) {
                // Stable and unstable ABL conditions
                obukhov_len = -utau * utau * utau * theta_mean /
                            (kappa * gravity * surf_temp_flux);
                zeta = zref / obukhov_len;
            } else {
                // Neutral conditions
                obukhov_len = std::numeric_limits<amrex::Real>::max();
                zeta = 0.0;
            }
            psi_m = calc_psi_m(zeta);
            psi_h = calc_psi_h(zeta);
            utau = kappa * vmag_mean / (std::log(zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(utau_iter - utau) > 1e-5) && iter <= max_iters);

        if (iter >= max_iters) {
            amrex::Print()
                << "MOData::update_fluxes: Convergence criteria not met after "
                << max_iters << " iterations"
                << "\nObuhov length = " << obukhov_len << " zeta = " << zeta
                << "\npsi_m = " << psi_m << " psi_h = " << psi_h
                << "\nutau = " << utau << " Tsurf = " << surf_temp
                << " q = " << surf_temp_flux << std::endl;
        }
    }
};

#endif /* ABLMOST_H */

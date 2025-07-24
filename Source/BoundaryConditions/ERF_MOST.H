#ifndef ERF_MOST_H
#define ERF_MOST_H

#include <ERF_Constants.H>
#include <ERF_MOSTRoughness.H>


/**
 * Structure of plain old data relevant to MOST BCs
 */
struct most_data
{
public:
    amrex::Real zref{10.0};           ///< Reference height (m)
    amrex::Real kappa{KAPPA};         ///< von Karman constant
    amrex::Real gravity{CONST_GRAV};  ///< Acceleration due to gravity (m/s^2)
    amrex::Real surf_temp_flux{0.0};  ///< Heat flux
    amrex::Real surf_moist_flux{0.0}; ///< Moisture flux

    amrex::Real Cnk_a{0.0185};        ///< Standard Charnock constant
    amrex::Real Cnk_d{30};            ///< Water depth
    bool visc{false};                 ///< Use viscous Charnock formulation

    const amrex::Real Bjr_beta = 1.2; // Empirical parameter from Beljaars 1995 QJRMS
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
#endif

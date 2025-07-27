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


#endif

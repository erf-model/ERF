#ifndef ABLMOST_H
#define ABLMOST_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <IndexDefines.H>
#include <ERF_Constants.H>

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
class ABLMost
{

public:

    // Constructor
    explicit ABLMost(const amrex::Vector<amrex::Geometry>& geom) : m_geom(geom)
    {
        amrex::ParmParse pp("erf");
        pp.query("most.surf_temp", surf_temp);
        pp.query("most.zref"     , zref);
        pp.query("most.z0"       , z0_const);

        int nlevs = m_geom.size();
        z_0.resize(nlevs);

        for (int lev = 0; lev < nlevs; lev++)
        {
            // TODO: generalize the "3" for the number of ghost cells
            amrex::Box bx = amrex::grow(m_geom[lev].Domain(),3);
            bx.setSmall(2,0);
            bx.setBig(2,0);
            z_0[lev].resize(bx,1);
            z_0[lev].setVal<amrex::RunOn::Device>(z0_const);
        }
    };

    void
    impose_most_bcs(const int lev, const amrex::Box& bx,
                    const amrex::Array4<amrex::Real>& dest_arr,
                    const amrex::Array4<amrex::Real>& cons_arr,
                    const amrex::Array4<amrex::Real>& velx_arr,
                    const amrex::Array4<amrex::Real>& vely_arr,
                    const amrex::Array4<amrex::Real>&  eta_arr,
                    const int idx, const int icomp, const int zlo);

    void
    update_fluxes(int lev,
                  amrex::MultiFab& S_new, amrex::MultiFab& U_new,
                  amrex::MultiFab& V_new, amrex::MultiFab& W_new,
                  int max_iters = 25);

    amrex::Real zref{0.2};           ///< Reference height (m)
    amrex::Real z0_const{0.1};       ///< Roughness height -- default constant value(m)
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

    enum ThetaCalcType {
        HEAT_FLUX = 0,      ///< Heat-flux specified
        SURFACE_TEMPERATURE ///< Surface temperature specified
    };

    ThetaCalcType alg_type{HEAT_FLUX};

    //AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    //amrex::Real phi_m() const
    //{
    //    return std::log(zref / z0_const) - calc_psi_m(zref / obukhov_len);
    //}

    //AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    //amrex::Real phi_m(amrex::Real z) const
    //{
    //    return std::log(z / z0_const) - calc_psi_m(z / obukhov_len);
    //}

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real phi_h(amrex::Real z0_in) const
    {
        return std::log(zref / z0_in) - calc_psi_h(zref / obukhov_len);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real phi_h(amrex::Real z, amrex::Real z0_in) const
    {
        return std::log(z / z0_in) - calc_psi_h(z / obukhov_len);
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

    void print() const
    {
        amrex::Print() << "ABLMost:\n";
        amrex::Print() << "  zref: " << zref << "\n";
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

    private:

        amrex::Vector<amrex::Geometry> m_geom;
        amrex::Vector<amrex::FArrayBox> z_0;
};

#endif /* ABLMOST_H */

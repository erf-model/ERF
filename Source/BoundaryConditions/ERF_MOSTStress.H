#ifndef ERF_MOSTStress_H
#define ERF_MOSTStress_H

#include "ERF_Constants.H"
#include "ERF_IndexDefines.H"
#include "ERF_MOSTUtils.H"
#include "ERF_Wstar.H"

/**
 * Adiabatic with constant roughness
 */
struct adiabatic
{
    /**
     * Construct adiabatic constant-roughness MOST data.
     *
     * @param[in] Tflux  prescribed heat flux
     * @param[in] Qvflux prescribed moisture flux
     */
    adiabatic (amrex::Real Tflux,
               amrex::Real Qvflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& /*max_iters*/,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<const amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& /*tm_arr*/,
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
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        olen_arr(i,j,k)   = bogus_large_value;
        u_star_arr(i,j,k) = mdata.kappa * umm_arr(i,j,k) / std::log(zref_arr(i,j,k) / z0_arr(i,j,k));
        t_star_arr(i,j,k) = zero;
        q_star_arr(i,j,k) = zero;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
};


/**
 * Adiabatic with charnock roughness
 */
struct adiabatic_charnock
{
    /**
     * Construct adiabatic Charnock-roughness MOST data.
     *
     * @param[in] Tflux    prescribed heat flux
     * @param[in] Qvflux   prescribed moisture flux
     * @param[in] cnk_a    Charnock parameter
     * @param[in] cnk_visc whether to include viscous roughness
     */
    adiabatic_charnock (amrex::Real Tflux,
                        amrex::Real Qvflux,
                        amrex::Real cnk_a,
                        bool cnk_visc)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.Cnk_a = cnk_a;
        mdata.visc  = cnk_visc;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        amrex::Real ustar  = zero;
        amrex::Real zref   = zref_arr(i,j,k);
        amrex::Real z0     = z0_arr(i,j,k);
        amrex::Real z0_old = z0;
        amrex::Real psi_m  = zero;
        amrex::Real umm    = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C      = std::log(zref / z0);

        // Fixed point iteration on roughness
        int iter_z = 0;
        do {
            // Transfer curr to old
            z0_old = z0;

            // Update
            C = std::log(zref / z0_old);
            ustar = mdata.kappa * umm / (C - psi_m);
            if (mdata.Cnk_a > 0) {
                z0 = Charnock_roughness(ustar, mdata.Cnk_a,
                                        air_viscosity(tm_arr(i,j,k)), mdata.visc);
            } else {
                z0 = COARE3_roughness(zref, umm, ustar,
                                      air_viscosity(tm_arr(i,j,k)), mdata.visc);
            }

            ++iter_z;
        } while ( (std::abs(z0 - z0_old) > tol_z) && (iter_z <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter_z < max_iters,
                                         "Maximum number of MOST roughness iterations reached.");
        C = std::log(zref / z0);

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = bogus_large_value;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = zero;
        q_star_arr(i,j,k) = zero;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
#ifdef AMREX_USE_FLOAT
    const amrex::Real tol_z = amrex::Real(1.0e-6);
#else
    const amrex::Real tol_z = amrex::Real(1.0e-10);
#endif
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};


/**
 * Adiabatic with modified charnock roughness
 */
struct adiabatic_mod_charnock
{
    /**
     * Construct adiabatic modified-Charnock MOST data.
     *
     * @param[in] Tflux  prescribed heat flux
     * @param[in] Qvflux prescribed moisture flux
     * @param[in] depth  water depth for modified Charnock roughness
     */
    adiabatic_mod_charnock (amrex::Real Tflux,
                            amrex::Real Qvflux,
                            amrex::Real depth)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.Cnk_d = depth;
        mdata.Cnk_b = mdata.Cnk_b1 * std::log(mdata.Cnk_b2 / mdata.Cnk_d);
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& /*tm_arr*/,
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
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        amrex::Real ustar  = zero;
        amrex::Real zref   = zref_arr(i,j,k);
        amrex::Real z0     = z0_arr(i,j,k);
        amrex::Real z0_old = z0;
        amrex::Real psi_m  = zero;
        amrex::Real umm    = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C      = std::log(zref / z0);

        // Fixed point iteration on roughness
        int iter_z = 0;
        do {
            // Transfer curr to old
            z0_old = z0;

            // Update
            C = std::log(zref / z0_old);
            ustar = mdata.kappa * umm / (C - psi_m);
            z0    = Mod_Charnock_roughness(ustar, mdata.Cnk_b);
            ++iter_z;
        } while ( (std::abs(z0 - z0_old) > tol_z) && (iter_z <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter_z < max_iters,
                                         "Maximum number of MOST roughness iterations reached.");
        C = std::log(zref / z0);

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = bogus_large_value;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = zero;
        q_star_arr(i,j,k) = zero;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
#ifdef AMREX_USE_FLOAT
    const amrex::Real tol_z = amrex::Real(1.0e-6);
#else
    const amrex::Real tol_z = amrex::Real(1.0e-10);
#endif
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};


/**
 * Adiabatic with Donelan roughness
 */
struct adiabatic_donelan
{
    /**
     * Construct adiabatic Donelan-roughness MOST data.
     *
     * @param[in] Tflux  prescribed heat flux
     * @param[in] Qvflux prescribed moisture flux
     */
    adiabatic_donelan (amrex::Real Tflux,
                       amrex::Real Qvflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        amrex::Real ustar  = zero;
        amrex::Real zref   = zref_arr(i,j,k);
        amrex::Real z0     = z0_arr(i,j,k);
        amrex::Real z0_old = z0;
        amrex::Real psi_m  = zero;
        amrex::Real umm    = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C      = std::log(zref / z0);

        // Fixed point iteration on roughness
        int iter_z = 0;
        do {
            // Transfer curr to old
            z0_old = z0;

            // Update
            C = std::log(zref / z0_old);
            ustar = mdata.kappa * umm / (C - psi_m);
            z0    = Donelan_roughness(ustar, air_viscosity(tm_arr(i,j,k)), mdata.visc);
            ++iter_z;
        } while ( (std::abs(z0 - z0_old) > tol_z) && (iter_z <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter_z < max_iters,
                                         "Maximum number of MOST roughness iterations reached.");
        C = std::log(zref / z0);

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = bogus_large_value;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = zero;
        q_star_arr(i,j,k) = zero;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
#ifdef AMREX_USE_FLOAT
    const amrex::Real tol_z = amrex::Real(1.0e-6);
#else
    const amrex::Real tol_z = amrex::Real(1.0e-10);
#endif
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};


/**
 * Adiabatic with wave-coupled roughness
 */
struct adiabatic_wave_coupled
{
    /**
     * Construct adiabatic wave-coupled roughness MOST data.
     *
     * @param[in] Tflux  prescribed heat flux
     * @param[in] Qvflux prescribed moisture flux
     */
    adiabatic_wave_coupled (amrex::Real Tflux,
                            amrex::Real Qvflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        amrex::Real ustar  = zero;
        amrex::Real zref   = zref_arr(i,j,k);
        amrex::Real z0     = z0_arr(i,j,k);
        amrex::Real z0_old = z0;
        amrex::Real psi_m  = zero;
        amrex::Real umm    = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C      = std::log(zref / z0);

        // Fixed point iteration on roughness
        int iter_z = 0;
        do {
            // Transfer curr to old
            z0_old = z0;

            // Update
            C = std::log(zref / z0_old);
            ustar = mdata.kappa * umm / (C - psi_m);
            z0    = WaveCoupled_roughness(Hwave_arr(i,j,k), Lwave_arr(i,j,k), eps, ustar,
                                          air_viscosity(tm_arr(i,j,k)), mdata.visc);
            ++iter_z;
        } while ( (std::abs(z0 - z0_old) > tol_z) && (iter_z <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter_z < max_iters,
                                         "Maximum number of MOST roughness iterations reached.");
        C = std::log(zref / z0);

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = bogus_large_value;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = zero;
        q_star_arr(i,j,k) = zero;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
#ifdef AMREX_USE_FLOAT
    const amrex::Real tol_z = amrex::Real(1.0e-6);
#else
    const amrex::Real tol_z = amrex::Real(1.0e-10);
#endif
#ifdef AMREX_USE_FLOAT
    const amrex::Real eps = amrex::Real(1e-8);
#else
    const amrex::Real eps = amrex::Real(1e-15);
#endif
    const amrex::Real WSMIN  = amrex::Real(0.1); // minimum wind speed
};


/**
 * Surface flux with constant roughness
 */
struct surface_flux
{
    /**
     * Construct specified-flux, constant-roughness MOST data.
     *
     * @param[in] Tflux      prescribed heat flux
     * @param[in] Qvflux     prescribed moisture flux
     * @param[in] cons_qflux whether the moisture flux is specified
     */
    surface_flux (amrex::Real Tflux,
                  amrex::Real Qvflux,
                  bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        spec_qflux = cons_qflux;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<const amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = zero;
        amrex::Real wstar = zero;
        amrex::Real tflux = zero;
        amrex::Real qflux = zero;
        amrex::Real zeta  = zero;
        amrex::Real psi_m = zero;
        amrex::Real psi_h = zero;
        amrex::Real Olen  = zero;
        amrex::Real zref  = zref_arr(i,j,k);
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == bogus_large_value) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(zref / z0_arr(i,j,k)) - psi_h); // <w'Qv'>
            tflux = mdata.surf_temp_flux*(one + epsv*qvm_arr(i,j,k)) + qflux*epsv*tm_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(zref / z0_arr(i,j,k)) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        olen_arr(i,j,k)   = Olen;
        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(zref / z0_arr(i,j,k)) - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol = amrex::Real(1.0e-5);
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};


/**
 * Surface flux with charnock roughness
 */
struct surface_flux_charnock
{
    /**
     * Construct specified-flux, Charnock-roughness MOST data.
     *
     * @param[in] Tflux    prescribed heat flux
     * @param[in] Qvflux   prescribed moisture flux
     * @param[in] cnk_a    Charnock parameter
     * @param[in] cnk_visc whether to include viscous roughness
     * @param[in] cons_qflux whether the moisture flux is specified
     */
    surface_flux_charnock (amrex::Real Tflux,
                           amrex::Real Qvflux,
                           amrex::Real cnk_a,
                           bool cnk_visc,
                           bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.Cnk_a = cnk_a;
        mdata.visc  = cnk_visc;
        spec_qflux  = cons_qflux;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = zero;
        amrex::Real wstar = zero;
        amrex::Real tflux = zero;
        amrex::Real qflux = zero;
        amrex::Real z0    = zero;
        amrex::Real zeta  = zero;
        amrex::Real psi_m = zero;
        amrex::Real psi_h = zero;
        amrex::Real Olen  = zero;
        amrex::Real zref  = zref_arr(i,j,k);
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == bogus_large_value) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            if (mdata.Cnk_a > 0) {
                z0 = Charnock_roughness(ustar, mdata.Cnk_a,
                                        air_viscosity(tm_arr(i,j,k)), mdata.visc);
            } else {
                z0 = COARE3_roughness(zref, umm, ustar,
                                      air_viscosity(tm_arr(i,j,k)), mdata.visc);
            }
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(zref / z0) - psi_h); // <w'Qv'>
            tflux = mdata.surf_temp_flux*(one + epsv*qvm_arr(i,j,k)) + qflux*epsv*tm_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = Olen;
        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(zref / z0) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(zref / z0_arr(i,j,k)) - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol = amrex::Real(1.0e-5);
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};


/**
 * Surface flux with modified charnock roughness
 */
struct surface_flux_mod_charnock
{
    /**
     * Construct specified-flux, modified-Charnock MOST data.
     *
     * @param[in] Tflux      prescribed heat flux
     * @param[in] Qvflux     prescribed moisture flux
     * @param[in] depth      water depth for modified Charnock roughness
     * @param[in] cons_qflux whether the moisture flux is specified
     */
    surface_flux_mod_charnock (amrex::Real Tflux,
                               amrex::Real Qvflux,
                               amrex::Real depth,
                               bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.Cnk_d = depth;
        mdata.Cnk_b = mdata.Cnk_b1 * std::log(mdata.Cnk_b2 / mdata.Cnk_d);
        spec_qflux  = cons_qflux;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = zero;
        amrex::Real wstar = zero;
        amrex::Real tflux = zero;
        amrex::Real qflux = zero;
        amrex::Real z0    = zero;
        amrex::Real zeta  = zero;
        amrex::Real psi_m = zero;
        amrex::Real psi_h = zero;
        amrex::Real Olen  = zero;
        amrex::Real zref  = zref_arr(i,j,k);
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == bogus_large_value) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = Mod_Charnock_roughness(ustar, mdata.Cnk_b);
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(zref / z0) - psi_h); // <w'Qv'>
            tflux = mdata.surf_temp_flux*(one + epsv*qvm_arr(i,j,k)) + qflux*epsv*tm_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = Olen;
        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(zref / z0) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(zref / z0_arr(i,j,k)) - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol = amrex::Real(1.0e-5);
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};


/**
 * Surface flux with donelan roughness
 */
struct surface_flux_donelan
{
    /**
     * Construct specified-flux, Donelan-roughness MOST data.
     *
     * @param[in] Tflux      prescribed heat flux
     * @param[in] Qvflux     prescribed moisture flux
     * @param[in] cons_qflux whether the moisture flux is specified
     */
    surface_flux_donelan (amrex::Real Tflux,
                          amrex::Real Qvflux,
                          bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        spec_qflux = cons_qflux;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = zero;
        amrex::Real wstar = zero;
        amrex::Real tflux = zero;
        amrex::Real qflux = zero;
        amrex::Real z0    = zero;
        amrex::Real zeta  = zero;
        amrex::Real psi_m = zero;
        amrex::Real psi_h = zero;
        amrex::Real Olen  = zero;
        amrex::Real zref  = zref_arr(i,j,k);
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == bogus_large_value) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = Donelan_roughness(ustar, air_viscosity(tm_arr(i,j,k)), mdata.visc);
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(zref / z0) - psi_h); // <w'Qv'>
            tflux = mdata.surf_temp_flux*(one + epsv*qvm_arr(i,j,k)) + qflux*epsv*tm_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = Olen;
        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(zref / z0) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(zref / z0_arr(i,j,k)) - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol = amrex::Real(1.0e-5);
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};


/**
 * Surface flux with wave-coupled roughness
 */
struct surface_flux_wave_coupled
{
    /**
     * Construct specified-flux, wave-coupled roughness MOST data.
     *
     * @param[in] Tflux      prescribed heat flux
     * @param[in] Qvflux     prescribed moisture flux
     * @param[in] cons_qflux whether the moisture flux is specified
     */
    surface_flux_wave_coupled (amrex::Real Tflux,
                               amrex::Real Qvflux,
                               bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        spec_qflux = cons_qflux;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = zero;
        amrex::Real wstar = zero;
        amrex::Real tflux = zero;
        amrex::Real qflux = zero;
        amrex::Real z0    = zero;
        amrex::Real zeta  = zero;
        amrex::Real psi_m = zero;
        amrex::Real psi_h = zero;
        amrex::Real Olen  = zero;
        amrex::Real zref  = zref_arr(i,j,k);
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == bogus_large_value) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = WaveCoupled_roughness(Hwave_arr(i,j,k), Lwave_arr(i,j,k), eps, ustar,
                                          air_viscosity(tm_arr(i,j,k)), mdata.visc);
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(zref / z0) - psi_h); // <w'Qv'>
            tflux = mdata.surf_temp_flux*(one + epsv*qvm_arr(i,j,k)) + qflux*epsv*tm_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = Olen;
        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(zref / z0) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(zref / z0_arr(i,j,k)) - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol = amrex::Real(1.0e-5);
#ifdef AMREX_USE_FLOAT
    const amrex::Real eps = amrex::Real(1e-8);
#else
    const amrex::Real eps = amrex::Real(1e-15);
#endif
    const amrex::Real WSMIN  = amrex::Real(0.1); // minimum wind speed
};


/**
 * Surface temperature with constant roughness
 */
struct surface_temp
{
    /**
     * Construct specified-temperature, constant-roughness MOST data.
     *
     * @param[in] Tflux      prescribed heat flux
     * @param[in] Qvflux     prescribed moisture flux
     * @param[in] cons_qflux whether the moisture flux is specified
     */
    surface_temp (amrex::Real Tflux,
                  amrex::Real Qvflux,
                  bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        spec_qflux = cons_qflux;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<const amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        amrex::Real Rib      = zero;
        amrex::Real zeta     = zero;
        amrex::Real zeta_old = zero;
        amrex::Real psi_m    = zero;
        amrex::Real psi_h    = zero;
        amrex::Real num      = zero;
        amrex::Real den      = zero;
        amrex::Real zref     = zref_arr(i,j,k);
        amrex::Real z0       = z0_arr(i,j,k);
        amrex::Real umm      = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C        = std::log(zref / z0);

        // First iteration we assume neutral (L -> inf)
        if (u_star_arr(i,j,k) == bogus_large_value) { olen_arr(i,j,k) = amrex::Real(1.0e3); }
        zeta  = zref / olen_arr(i,j,k);

        // Water vapor in atmos and surface
        amrex::Real qv_s, qv_a;
        if (q_surf_arr(i,j,k) > zero) {
            qv_s = q_surf_arr(i,j,k);
        } else {
            // First iteration and no qv_surf was specified
            // Use mean since there will be no flux
            qv_s = qvm_arr(i,j,k);
        }
        qv_a = qvm_arr(i,j,k);

        // update w* and Umagmean from Beljaars (1995)
        if (w_star_arr) {
            // NOTE: Thv flux is lagged, similar to WRF
            psi_m = sfuns.calc_psi_m2(zeta);
            psi_h = sfuns.calc_psi_h2(zeta);
            amrex::Real ustar = mdata.kappa * umm / (C - psi_m);
            amrex::Real tstar = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
            amrex::Real qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -ustar * mdata.kappa * (qv_a - qv_s) / (C - psi_h);
            amrex::Real tflux = -ustar*tstar*(one + epsv*qvm_arr(i,j,k)) + epsv*tm_arr(i,j,k)*qflux;
            w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
            amrex::Real wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
            umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
            umm = std::max(umm, WSMIN);
        }

        // Bulk Richardson number w/ moisture
        amrex::Real thv_s = t_surf_arr(i,j,k) * (one + epsv*qv_s);
        amrex::Real thv_a =     tm_arr(i,j,k) * (one + epsv*qv_a);
        Rib = ( (mdata.gravity * zref) / tm_arr(i,j,k) ) *
              ( (thv_a - thv_s) / (umm * umm) );
        Rib = std::min(std::max(Rib,-amrex::Real(4.0)),amrex::Real(4.0));

        // Fixed point iteration on zeta
        int iter = 0;
        do {
            // Transfer curr to old
            zeta_old = zeta;

            // Stability functions
            psi_m = sfuns.calc_psi_m2(zeta_old);
            psi_h = sfuns.calc_psi_h2(zeta_old);

            // Limiting
            num = std::max(C - psi_m, amrex::Real(1.0));
            den = std::max(C - psi_h, amrex::Real(1.0));

            // Update with under relaxation
            zeta = (one - alpha) * zeta_old + alpha * Rib * num * num / den;

            ++iter;
        } while ( (std::abs(zeta - zeta_old) > tol) && (iter <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        olen_arr(i,j,k)   = zref / zeta;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (C - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) / (C - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol   = amrex::Real(1.0e-3);
    const amrex::Real alpha = myhalf;
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};


/**
 * Surface temperature with charnock roughness
 */
struct surface_temp_charnock
{
    /**
     * Construct specified-temperature, Charnock-roughness MOST data.
     *
     * @param[in] Tflux      prescribed heat flux
     * @param[in] Qvflux     prescribed moisture flux
     * @param[in] cnk_a      Charnock parameter
     * @param[in] cnk_visc   whether to include viscous roughness
     * @param[in] cons_qflux whether the moisture flux is specified
     */
    surface_temp_charnock (amrex::Real Tflux,
                           amrex::Real Qvflux,
                           amrex::Real cnk_a,
                           bool cnk_visc,
                           bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.Cnk_a = cnk_a;
        mdata.visc  = cnk_visc;
        spec_qflux  = cons_qflux;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        amrex::Real Rib      = zero;
        amrex::Real zeta     = zero;
        amrex::Real zeta_old = zero;
        amrex::Real psi_m    = zero;
        amrex::Real psi_h    = zero;
        amrex::Real num      = zero;
        amrex::Real den      = zero;
        amrex::Real z0       = z0_arr(i,j,k);
        amrex::Real z0_old   = z0;
        amrex::Real ustar    = zero;
        amrex::Real zref     = zref_arr(i,j,k);
        amrex::Real umm      = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C        = std::log(zref / z0);

        // First iteration we assume neutral (L -> inf)
        if (u_star_arr(i,j,k) == bogus_large_value) {
            olen_arr(i,j,k) = amrex::Real(1.0e3);
        }
        zeta  = zref / olen_arr(i,j,k);

        // Water vapor in atmos and surface
        amrex::Real qv_s, qv_a;
        if (q_surf_arr(i,j,k) > zero) {
            qv_s = q_surf_arr(i,j,k);
        } else {
            // First iteration and no qv_surf was specified
            // Use mean since there will be no flux
            qv_s = qvm_arr(i,j,k);
        }
        qv_a = qvm_arr(i,j,k);

        // update w* and Umagmean from Beljaars (1995)
        if (w_star_arr) {
            // NOTE: Thv flux is lagged, similar to WRF
            psi_m = sfuns.calc_psi_m2(zeta);
            psi_h = sfuns.calc_psi_h2(zeta);
                        ustar = mdata.kappa * umm / (C - psi_m);
            amrex::Real tstar = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
            amrex::Real qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -ustar * mdata.kappa * (qv_a - qv_s) / (C - psi_h);
            amrex::Real tflux = -ustar*tstar*(one + epsv*qvm_arr(i,j,k)) + epsv*tm_arr(i,j,k)*qflux;
            w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
            amrex::Real wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
            umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
            umm = std::max(umm, WSMIN);
        }

        // Bulk Richardson number w/ moisture
        amrex::Real thv_s = t_surf_arr(i,j,k) * (one + epsv*qv_s);
        amrex::Real thv_a =     tm_arr(i,j,k) * (one + epsv*qv_a);
        Rib = ( (mdata.gravity * zref) / tm_arr(i,j,k) ) *
              ( (thv_a - thv_s) / (umm * umm) );
        Rib = std::min(std::max(Rib,-amrex::Real(4.0)),amrex::Real(4.0));

        // Fixed point iteration on zeta
        int iter = 0;
        do {
            // Transfer curr to old
            zeta_old = zeta;

            // Stability functions
            psi_m = sfuns.calc_psi_m2(zeta_old);
            psi_h = sfuns.calc_psi_h2(zeta_old);

            // Fixed point iteration on roughness
            int iter_z = 0;
            do {
                // Transfer curr to old
                z0_old = z0;

                // Update
                C = std::log(zref / z0_old);
                ustar = mdata.kappa * umm / (C - psi_m);
                if (mdata.Cnk_a > 0) {
                    z0 = Charnock_roughness(ustar, mdata.Cnk_a,
                                            air_viscosity(tm_arr(i,j,k)), mdata.visc);
                } else {
                    z0 = COARE3_roughness(zref, umm, ustar,
                                          air_viscosity(tm_arr(i,j,k)), mdata.visc);
                }

                ++iter_z;
            } while ( (std::abs(z0 - z0_old) > tol_z) && (iter_z <= max_iters) );
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter_z < max_iters,
                                             "Maximum number of MOST roughness iterations reached.");
            C = std::log(zref / z0);

            // Limiting
            num = std::max(C - psi_m, amrex::Real(1.0));
            den = std::max(C - psi_h, amrex::Real(1.0));

            // Update with under relaxation
            zeta = (amrex::Real(1.0) - alpha) * zeta_old + alpha * Rib * num * num / den;

            ++iter;
        } while ( (std::abs(zeta - zeta_old) > tol) && (iter <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = zref / zeta;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (C - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) / (C - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol   = amrex::Real(1.0e-3);
#ifdef AMREX_USE_FLOAT
    const amrex::Real tol_z = amrex::Real(1.0e-6);
#else
    const amrex::Real tol_z = amrex::Real(1.0e-10);
#endif
    const amrex::Real alpha = myhalf;
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};


/**
 * Surface temperature with modified charnock roughness
 */
struct surface_temp_mod_charnock
{
    /**
     * Construct specified-temperature, modified-Charnock MOST data.
     *
     * @param[in] Tflux      prescribed heat flux
     * @param[in] Qvflux     prescribed moisture flux
     * @param[in] depth      water depth for modified Charnock roughness
     * @param[in] cons_qflux whether the moisture flux is specified
     */
    surface_temp_mod_charnock (amrex::Real Tflux,
                               amrex::Real Qvflux,
                               amrex::Real depth,
                               bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.Cnk_d = depth;
        mdata.Cnk_b = mdata.Cnk_b1 * std::log(mdata.Cnk_b2 / mdata.Cnk_d);
        spec_qflux  = cons_qflux;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        amrex::Real Rib      = zero;
        amrex::Real zeta     = zero;
        amrex::Real zeta_old = zero;
        amrex::Real psi_m    = zero;
        amrex::Real psi_h    = zero;
        amrex::Real num      = zero;
        amrex::Real den      = zero;
        amrex::Real z0       = z0_arr(i,j,k);
        amrex::Real z0_old   = z0;
        amrex::Real ustar    = zero;
        amrex::Real zref     = zref_arr(i,j,k);
        amrex::Real umm      = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C        = std::log(zref / z0);

        // First iteration we assume neutral (L -> inf)
        if (u_star_arr(i,j,k) == bogus_large_value) {
            olen_arr(i,j,k) = amrex::Real(1000);
        }
        zeta  = zref / olen_arr(i,j,k);

        // Water vapor in atmos and surface
        amrex::Real qv_s, qv_a;
        if (q_surf_arr(i,j,k) > zero) {
            qv_s = q_surf_arr(i,j,k);
        } else {
            // First iteration and no qv_surf was specified
            // Use mean since there will be no flux
            qv_s = qvm_arr(i,j,k);
        }
        qv_a = qvm_arr(i,j,k);

        // update w* and Umagmean from Beljaars (1995)
        if (w_star_arr) {
            // NOTE: Thv flux is lagged, similar to WRF
            psi_m = sfuns.calc_psi_m2(zeta);
            psi_h = sfuns.calc_psi_h2(zeta);
                        ustar = mdata.kappa * umm / (C - psi_m);
            amrex::Real tstar = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
            amrex::Real qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -ustar * mdata.kappa * (qv_a - qv_s) / (C - psi_h);
            amrex::Real tflux = -ustar*tstar*(one + epsv*qvm_arr(i,j,k)) + epsv*tm_arr(i,j,k)*qflux;
            w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
            amrex::Real wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
            umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
            umm = std::max(umm, WSMIN);
        }

        // Bulk Richardson number w/ moisture
        amrex::Real thv_s = t_surf_arr(i,j,k) * (one + epsv*qv_s);
        amrex::Real thv_a =     tm_arr(i,j,k) * (one + epsv*qv_a);
        Rib = ( (mdata.gravity * zref) / tm_arr(i,j,k) ) *
              ( (thv_a - thv_s) / (umm * umm) );
        Rib = std::min(std::max(Rib,-amrex::Real(4.0)),amrex::Real(4.0));

        // Fixed point iteration on zeta
        int iter = 0;
        do {
            // Transfer curr to old
            zeta_old = zeta;

            // Stability functions
            psi_m = sfuns.calc_psi_m2(zeta_old);
            psi_h = sfuns.calc_psi_h2(zeta_old);

            // Fixed point iteration on roughness
            int iter_z = 0;
            do {
                // Transfer curr to oldgoto
                z0_old = z0;

                // Update
                C = std::log(zref / z0_old);
                ustar = mdata.kappa * umm / (C - psi_m);
                z0    = Mod_Charnock_roughness(ustar, mdata.Cnk_b);

                ++iter_z;
            } while ( (std::abs(z0 - z0_old) > tol_z) && (iter_z <= max_iters) );
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter_z < max_iters,
                                             "Maximum number of MOST roughness iterations reached.");
            C = std::log(zref / z0);

            // Limiting
            num = std::max(C - psi_m, amrex::Real(1.0));
            den = std::max(C - psi_h, amrex::Real(1.0));

            // Update with under relaxation
            zeta = (amrex::Real(1.0) - alpha) * zeta_old + alpha * Rib * num * num / den;

            ++iter;
        } while ( (std::abs(zeta - zeta_old) > tol) && (iter <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = zref / zeta;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (C - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) / (C - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol   = amrex::Real(1.0e-3);
#ifdef AMREX_USE_FLOAT
    const amrex::Real tol_z = amrex::Real(1.0e-6);
#else
    const amrex::Real tol_z = amrex::Real(1.0e-10);
#endif
    const amrex::Real alpha = myhalf;
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};


/**
 * Surface temperature with donelan roughness
 */
struct surface_temp_donelan
{
    /**
     * Construct specified-temperature, Donelan-roughness MOST data.
     *
     * @param[in] Tflux      prescribed heat flux
     * @param[in] Qvflux     prescribed moisture flux
     * @param[in] cons_qflux whether the moisture flux is specified
     */
    surface_temp_donelan (amrex::Real Tflux,
                          amrex::Real Qvflux,
                          bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        spec_qflux = cons_qflux;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        amrex::Real Rib      = zero;
        amrex::Real zeta     = zero;
        amrex::Real zeta_old = zero;
        amrex::Real psi_m    = zero;
        amrex::Real psi_h    = zero;
        amrex::Real num      = zero;
        amrex::Real den      = zero;
        amrex::Real z0       = z0_arr(i,j,k);
        amrex::Real z0_old   = z0;
        amrex::Real ustar    = zero;
        amrex::Real zref     = zref_arr(i,j,k);
        amrex::Real umm      = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C        = std::log(zref / z0);

        // First iteration we assume neutral (L -> inf)
        if (u_star_arr(i,j,k) == bogus_large_value) {
            olen_arr(i,j,k) = amrex::Real(1.0e3);
        }
        zeta  = zref / olen_arr(i,j,k);

        // Water vapor in atmos and surface
        amrex::Real qv_s, qv_a;
        if (q_surf_arr(i,j,k) > zero) {
            qv_s = q_surf_arr(i,j,k);
        } else {
            // First iteration and no qv_surf was specified
            // Use mean since there will be no flux
            qv_s = qvm_arr(i,j,k);
        }
        qv_a = qvm_arr(i,j,k);

        // update w* and Umagmean from Beljaars (1995)
        if (w_star_arr) {
            // NOTE: Thv flux is lagged, similar to WRF
            psi_m = sfuns.calc_psi_m2(zeta);
            psi_h = sfuns.calc_psi_h2(zeta);
                        ustar = mdata.kappa * umm / (C - psi_m);
            amrex::Real tstar = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
            amrex::Real qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -ustar * mdata.kappa * (qv_a - qv_s) / (C - psi_h);
            amrex::Real tflux = -ustar*tstar*(one + epsv*qvm_arr(i,j,k)) + epsv*tm_arr(i,j,k)*qflux;
            w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
            amrex::Real wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
            umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
            umm = std::max(umm, WSMIN);
        }

        // Bulk Richardson number w/ moisture
        amrex::Real thv_s = t_surf_arr(i,j,k) * (one + epsv*qv_s);
        amrex::Real thv_a =     tm_arr(i,j,k) * (one + epsv*qv_a);
        Rib = ( (mdata.gravity * zref) / tm_arr(i,j,k) ) *
              ( (thv_a - thv_s) / (umm * umm) );
        Rib = std::min(std::max(Rib,-amrex::Real(4.0)),amrex::Real(4.0));

        // Fixed point iteration on zeta
        int iter = 0;
        do {
            // Transfer curr to old
            zeta_old = zeta;

            // Stability functions
            psi_m = sfuns.calc_psi_m2(zeta_old);
            psi_h = sfuns.calc_psi_h2(zeta_old);

            // Fixed point iteration on roughness
            int iter_z = 0;
            do {
                // Transfer curr to old
                z0_old = z0;

                // Update
                C     = std::log(zref / z0_old);
                ustar = mdata.kappa * umm / (C - psi_m);
                z0    = Donelan_roughness(ustar, air_viscosity(tm_arr(i,j,k)), mdata.visc);

                ++iter_z;
            } while ( (std::abs(z0 - z0_old) > tol_z) && (iter_z <= max_iters) );
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter_z < max_iters,
                                             "Maximum number of MOST roughness iterations reached.");
            C = std::log(zref / z0);

            // Limiting
            num = std::max(C - psi_m, amrex::Real(1.0));
            den = std::max(C - psi_h, amrex::Real(1.0));

            // Update with under relaxation
            zeta = (amrex::Real(1.0) - alpha) * zeta_old + alpha * Rib * num * num / den;

            ++iter;
        } while ( (std::abs(zeta - zeta_old) > tol) && (iter <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = zref / zeta;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (C - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) / (C - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol   = amrex::Real(1.0e-3);
#ifdef AMREX_USE_FLOAT
    const amrex::Real tol_z = amrex::Real(1.0e-6);
#else
    const amrex::Real tol_z = amrex::Real(1.0e-10);
#endif
    const amrex::Real alpha = myhalf;
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};


/**
 * Surface temperature with wave-coupled roughness
 */
struct surface_temp_wave_coupled
{
    /**
     * Construct specified-temperature, wave-coupled roughness MOST data.
     *
     * @param[in] Tflux      prescribed heat flux
     * @param[in] Qvflux     prescribed moisture flux
     * @param[in] cons_qflux whether the moisture flux is specified
     */
    surface_temp_wave_coupled (amrex::Real Tflux,
                               amrex::Real Qvflux,
                               bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        spec_qflux = cons_qflux;
    }

    /**
     * Iterate MOST state for one surface point.
     *
     * @par Calling sequence
     * i, j, k, max_iters, zref_arr, z0_arr, umm_arr, tm_arr, tvm_arr,
     * qvm_arr, u_star_arr, w_star_arr, t_star_arr, q_star_arr,
     * t_surf_arr, q_surf_arr, olen_arr, pblh_arr, Hwave_arr, Lwave_arr,
     * eta_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        amrex::Real Rib      = zero;
        amrex::Real zeta     = zero;
        amrex::Real zeta_old = zero;
        amrex::Real psi_m    = zero;
        amrex::Real psi_h    = zero;
        amrex::Real num      = zero;
        amrex::Real den      = zero;
        amrex::Real z0       = z0_arr(i,j,k);
        amrex::Real z0_old   = z0;
        amrex::Real ustar    = zero;
        amrex::Real zref     = zref_arr(i,j,k);
        amrex::Real umm      = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C        = std::log(zref / z0);

        // First iteration we assume neutral (L -> inf)
        if (u_star_arr(i,j,k) == bogus_large_value) {
            olen_arr(i,j,k) = amrex::Real(1.0e3);
        }
        zeta  = zref / olen_arr(i,j,k);

        // Water vapor in atmos and surface
        amrex::Real qv_s, qv_a;
        if (q_surf_arr(i,j,k) > zero) {
            qv_s = q_surf_arr(i,j,k);
        } else {
            // First iteration and no qv_surf was specified
            // Use mean since there will be no flux
            qv_s = qvm_arr(i,j,k);
        }
        qv_a = qvm_arr(i,j,k);

        // update w* and Umagmean from Beljaars (1995)
        if (w_star_arr) {
            // NOTE: Thv flux is lagged, similar to WRF
            psi_m = sfuns.calc_psi_m2(zeta);
            psi_h = sfuns.calc_psi_h2(zeta);
            ustar = mdata.kappa * umm / (C - psi_m);
            amrex::Real tstar = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
            amrex::Real qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -ustar * mdata.kappa * (qv_a - qv_s) / (C - psi_h);
            amrex::Real tflux = -ustar*tstar*(one + epsv*qvm_arr(i,j,k)) + epsv*tm_arr(i,j,k)*qflux;
            w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
            amrex::Real wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
            umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
            umm = std::max(umm, WSMIN);
        }

        // Bulk Richardson number w/ moisture
        amrex::Real thv_s = t_surf_arr(i,j,k) * (one + epsv*qv_s);
        amrex::Real thv_a =     tm_arr(i,j,k) * (one + epsv*qv_a);
        Rib = ( (mdata.gravity * zref) / tm_arr(i,j,k) ) *
              ( (thv_a - thv_s) / (umm * umm) );
        Rib = std::min(std::max(Rib,-amrex::Real(4.0)),amrex::Real(4.0));

        // Fixed point iteration on zeta
        int iter = 0;
        do {
            // Transfer curr to old
            zeta_old = zeta;

            // Stability functions
            psi_m = sfuns.calc_psi_m2(zeta_old);
            psi_h = sfuns.calc_psi_h2(zeta_old);

            // Fixed point iteration on roughness
            int iter_z = 0;
            do {
                // Transfer curr to old
                z0_old = z0;

                // Update
                C     = std::log(zref / z0_old);
                ustar = mdata.kappa * umm / (C - psi_m);
                z0    = WaveCoupled_roughness(Hwave_arr(i,j,k), Lwave_arr(i,j,k), eps, ustar,
                                              air_viscosity(tm_arr(i,j,k)), mdata.visc);

                ++iter_z;
            } while ( (std::abs(z0 - z0_old) > tol_z) && (iter_z <= max_iters) );
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter_z < max_iters,
                                             "Maximum number of MOST roughness iterations reached.");
            C = std::log(zref / z0);

            // Limiting
            num = std::max(C - psi_m, amrex::Real(1.0));
            den = std::max(C - psi_h, amrex::Real(1.0));

            // Update with under relaxation
            zeta = (amrex::Real(1.0) - alpha) * zeta_old + alpha * Rib * num * num / den;

            ++iter;
        } while ( (std::abs(zeta - zeta_old) > tol) && (iter <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = zref / zeta;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (C - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) / (C - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol   = amrex::Real(1.0e-3);
#ifdef AMREX_USE_FLOAT
    const amrex::Real tol_z = amrex::Real(1.0e-6);
    const amrex::Real eps   = amrex::Real(1e-8);
#else
    const amrex::Real tol_z = amrex::Real(1.0e-10);
    const amrex::Real eps   = amrex::Real(1e-15);
#endif
    const amrex::Real alpha = myhalf;
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};


/**
 * Moeng flux formulation
 */
struct moeng_flux
{
    /**
     * Construct the Moeng flux calculator.
     */
    moeng_flux () {}

    /**
     * Compute moisture flux for one surface point.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, qvm_arr,
     * u_star_arr, q_star_arr, q_surf_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_q_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& qvm_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr,
                    const amrex::Array4<const amrex::Real>& q_star_arr,
                    const amrex::Array4<const amrex::Real>& q_surf_arr) const
    {
        amrex::Real rho   = cons_arr(i,j,k,Rho_comp);
        amrex::Real qv    = cons_arr(i,j,k,RhoQ1_comp) / rho;
        amrex::Real velx  = myhalf * ( velx_arr(i,j,k) + velx_arr(i+1,j  ,k) );
        amrex::Real vely  = myhalf * ( vely_arr(i,j,k) + vely_arr(i  ,j+1,k) );

        amrex::Real qv_mean  = qvm_arr(i,j,0);
        amrex::Real ustar    = u_star_arr(i,j,0);
        amrex::Real qstar    = q_star_arr(i,j,0);
        amrex::Real qv_surf  = q_surf_arr(i,j,0);
        amrex::Real wsp_mean = umm_arr(i,j,0);
        wsp_mean = std::max(wsp_mean, WSMIN);

        amrex::Real wsp     = std::sqrt(velx*velx+vely*vely);
        amrex::Real num1    = wsp * (qv_mean-qv_surf);
        amrex::Real num2    = wsp_mean * (qv-qv_mean);

        // NOTE: this is rho*<Qv'w'> = -K dQvdz
        amrex::Real moflux  = (std::abs(qstar) > eps) ?
                              -rho*qstar*ustar*(num1+num2)/((qv_mean-qv_surf)*wsp_mean) : zero;

        return moflux;
    }

    /**
     * Compute heat flux for one surface point.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, tm_arr,
     * u_star_arr, t_star_arr, t_surf_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_t_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& tm_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr,
                    const amrex::Array4<const amrex::Real>& t_star_arr,
                    const amrex::Array4<const amrex::Real>& t_surf_arr) const
    {
        amrex::Real rho   = cons_arr(i,j,k,Rho_comp);
        amrex::Real theta = cons_arr(i,j,k,RhoTheta_comp) / rho;
        amrex::Real velx  = myhalf * ( velx_arr(i,j,k) + velx_arr(i+1,j  ,k) );
        amrex::Real vely  = myhalf * ( vely_arr(i,j,k) + vely_arr(i  ,j+1,k) );

        amrex::Real theta_mean  = tm_arr(i,j,0);
        amrex::Real ustar       = u_star_arr(i,j,0);
        amrex::Real tstar       = t_star_arr(i,j,0);
        amrex::Real theta_surf  = t_surf_arr(i,j,0);
        amrex::Real wsp_mean    = umm_arr(i,j,0);
        wsp_mean = std::max(wsp_mean, WSMIN);

        amrex::Real wsp     = std::sqrt(velx*velx+vely*vely);
        amrex::Real num1    = wsp * (theta_mean-theta_surf);
        amrex::Real num2    = wsp_mean * (theta-theta_mean);

        // NOTE: this is rho*<T'w'> = -K dTdz
        amrex::Real moflux  = (std::abs(tstar) > eps) ?
                              -rho*tstar*ustar*(num1+num2)/((theta_mean-theta_surf)*wsp_mean) : zero;

        return moflux;
    }

    /**
     * Compute x-momentum flux for one surface face.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, um_arr,
     * u_star_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_u_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& um_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr) const
    {
        amrex::Real velx  = velx_arr(i,j,k);
        amrex::Real vely  = fourth * ( vely_arr(i  ,j,k) + vely_arr(i  ,j+1,k)
                                     + vely_arr(i-1,j,k) + vely_arr(i-1,j+1,k) );
        amrex::Real rho   = myhalf  * ( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real umean    = um_arr(i,j,0);
        amrex::Real ustar    = myhalf * ( u_star_arr(i-1,j,0) + u_star_arr(i,j,0) );
        amrex::Real wsp_mean = myhalf * (    umm_arr(i-1,j,0) +    umm_arr(i,j,0) );
        wsp_mean = std::max(wsp_mean, WSMIN);

        // Note: The surface mean shear stress is decomposed into tau_xz by
        //       multiplying the modeled shear stress (rho*ustar^2) with
        //       a factor of umean/wsp_mean for directionality; this factor
        //       modifies the denominator from what is in Moeng amrex::Real(1984.)
        amrex::Real wsp     = std::sqrt(velx*velx+vely*vely);
        amrex::Real num1    = wsp * umean;
        amrex::Real num2    = wsp_mean * (velx-umean);

        // NOTE: this is rho*<u'w'> = -K dudz
        amrex::Real stressx = -rho*ustar*ustar * (num1+num2)/(wsp_mean*wsp_mean);

        return stressx;
    }

    /**
     * Compute y-momentum flux for one surface face.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, vm_arr,
     * u_star_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_v_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& vm_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr) const
    {
        amrex::Real velx  = fourth * ( velx_arr(i,j  ,k) + velx_arr(i+1,j  ,k)
                                     + velx_arr(i,j-1,k) + velx_arr(i+1,j-1,k) );
        amrex::Real vely  = vely_arr(i,j,k);
        amrex::Real rho   = myhalf  * ( cons_arr(i,j-1,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real vmean     = vm_arr(i,j,0);
        amrex::Real ustar     = myhalf * ( u_star_arr(i,j-1,0) + u_star_arr(i,j,0) );
        amrex::Real wsp_mean  = myhalf * (    umm_arr(i,j-1,0) +    umm_arr(i,j,0) );
        wsp_mean = std::max(wsp_mean, WSMIN);

        // Note: The surface mean shear stress is decomposed into tau_yz by
        //       multiplying the modeled shear stress (rho*ustar^2) with
        //       a factor of vmean/wsp_mean for directionality; this factor
        //       modifies the denominator from what is in Moeng amrex::Real(1984.)
        amrex::Real wsp     = std::sqrt(velx*velx+vely*vely);
        amrex::Real num1    = wsp * vmean;
        amrex::Real num2    = wsp_mean * (vely-vmean);

        // NOTE: this is rho*<v'w'> = -K dvdz
        amrex::Real stressy = -rho*ustar*ustar * (num1+num2)/(wsp_mean*wsp_mean);

        return stressy;
    }

private:
#ifdef AMREX_USE_FLOAT
    const amrex::Real eps = amrex::Real(1e-8);
#else
    const amrex::Real eps = amrex::Real(1e-15);
#endif
    const amrex::Real   WSMIN = amrex::Real(0.1); // minimum wind speed
};


/**
 * Custom flux formulation
 */
struct custom_flux
{
    /**
     * Construct the custom flux calculator.
     *
     * @param[in] specified_rho_surf whether the supplied custom fluxes include density
     */
    custom_flux (bool specified_rho_surf)
    : fluxes_include_rho(specified_rho_surf)
    {}

    /**
     * Compute moisture flux for one surface point.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, qvm_arr,
     * u_star_arr, q_star_arr, q_surf_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_q_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*qm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& q_star_arr,
                    const amrex::Array4<const amrex::Real>& /*q_surf_arr*/) const
    {
        amrex::Real rho   = (fluxes_include_rho) ? one : cons_arr(i,j,k,Rho_comp);
        amrex::Real qstar = q_star_arr(i,j,0);

        // NOTE: this is rho*<q'w'> = -K dqdz
        amrex::Real moflux  = (std::abs(qstar) > eps) ? rho * qstar : zero;

        return moflux;
    }

    /**
     * Compute heat flux for one surface point.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, tm_arr,
     * u_star_arr, t_star_arr, t_surf_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_t_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*tm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& t_star_arr,
                    const amrex::Array4<const amrex::Real>& /*t_surf_arr*/) const
    {
        amrex::Real rho   = (fluxes_include_rho) ? one : cons_arr(i,j,k,Rho_comp);
        amrex::Real tstar = t_star_arr(i,j,0);

        // NOTE: this is rho*<T'w'> = -K dTdz
        amrex::Real moflux  = (std::abs(tstar) > eps) ? rho * tstar : zero;

        return moflux;
    }

    /**
     * Compute x-momentum flux for one surface face.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, um_arr,
     * u_star_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_u_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*um_arr*/,
                    const amrex::Array4<const amrex::Real>& u_star_arr) const
    {
        amrex::Real velx  = velx_arr(i,j,k);
        amrex::Real vely  = fourth * ( vely_arr(i  ,j,k) + vely_arr(i  ,j+1,k)
                                     + vely_arr(i-1,j,k) + vely_arr(i-1,j+1,k) );
        amrex::Real rho   = (fluxes_include_rho) ? one : myhalf * ( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real ustar   = myhalf * ( u_star_arr(i-1,j,0) + u_star_arr(i,j,0) );
        amrex::Real wsp     = std::sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<u'w'> = -K dudz
        amrex::Real stressx = -rho * ustar * ustar * velx / wsp;

        return stressx;
    }

    /**
     * Compute y-momentum flux for one surface face.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, vm_arr,
     * u_star_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_v_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vm_arr*/,
                    const amrex::Array4<const amrex::Real>& u_star_arr) const
    {
        amrex::Real velx  = fourth * ( velx_arr(i,j  ,k) + velx_arr(i+1,j  ,k)
                                     + velx_arr(i,j-1,k) + velx_arr(i+1,j-1,k) );
        amrex::Real vely  = vely_arr(i,j,k);
        amrex::Real rho   = (fluxes_include_rho) ? one : myhalf * ( cons_arr(i,j-1,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real ustar   = myhalf * ( u_star_arr(i,j-1,0) + u_star_arr(i,j,0) );
        amrex::Real wsp     = std::sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<v'w'> = -K dvdz
        amrex::Real stressy = -rho * ustar * ustar * vely / wsp;

        return stressy;
    }

private:
#ifdef AMREX_USE_FLOAT
    const amrex::Real eps = amrex::Real(1e-8);
#else
    const amrex::Real eps = amrex::Real(1e-15);
#endif
    const bool fluxes_include_rho{false};
};


/**
 * Bulk coefficient flux formulation
 */
struct bulk_coeff_flux
{
    /**
     * Construct the bulk-coefficient flux calculator.
     *
     * @param[in] m_Cd drag coefficient
     * @param[in] m_Ch heat-transfer coefficient
     * @param[in] m_Cq moisture-transfer coefficient
     */
    bulk_coeff_flux (amrex::Real m_Cd,
                     amrex::Real m_Ch,
                     amrex::Real m_Cq)
    {
        mdata.Cd = m_Cd;
        mdata.Ch = m_Ch;
        mdata.Cq = m_Cq;
    }

    /**
     * Compute moisture flux for one surface point.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, qvm_arr,
     * u_star_arr, q_star_arr, q_surf_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_q_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*qm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& /*q_star_arr*/,
                    const amrex::Array4<const amrex::Real>& q_surf_arr) const
    {
        amrex::Real rho    = cons_arr(i,j,k,Rho_comp);
        amrex::Real qv     = cons_arr(i,j,k,RhoQ1_comp)/cons_arr(i,j,k,Rho_comp);
        amrex::Real qvsurf = q_surf_arr(i,j,0);
        amrex::Real velx   = myhalf * ( velx_arr(i,j,k) + velx_arr(i+1,j  ,k) );
        amrex::Real vely   = myhalf * ( vely_arr(i,j,k) + vely_arr(i  ,j+1,k) );
        amrex::Real wsp    = std::sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<q'w'> = -K dqdz
        amrex::Real moflux  = rho * mdata.Cq * wsp * (qvsurf - qv);

        return moflux;
    }

    /**
     * Compute heat flux for one surface point.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, tm_arr,
     * u_star_arr, t_star_arr, t_surf_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_t_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*tm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& /*t_star_arr*/,
                    const amrex::Array4<const amrex::Real>& t_surf_arr) const
    {
        amrex::Real rho    = cons_arr(i,j,k,Rho_comp);
        amrex::Real th     = cons_arr(i,j,k,RhoTheta_comp)/cons_arr(i,j,k,Rho_comp);
        amrex::Real thsurf = t_surf_arr(i,j,0);
        amrex::Real velx   = myhalf * ( velx_arr(i,j,k) + velx_arr(i+1,j  ,k) );
        amrex::Real vely   = myhalf * ( vely_arr(i,j,k) + vely_arr(i  ,j+1,k) );
        amrex::Real wsp    = std::sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<T'w'> = -K dTdz
        amrex::Real moflux  = rho * mdata.Ch * wsp * (thsurf - th);

        return moflux;
    }

    /**
     * Compute x-momentum flux for one surface face.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, um_arr,
     * u_star_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_u_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*um_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/) const
    {
        amrex::Real rho   = myhalf * ( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );
        amrex::Real velx  = velx_arr(i,j,k);
        amrex::Real vely  = fourth * ( vely_arr(i  ,j,k) + vely_arr(i  ,j+1,k)
                                     + vely_arr(i-1,j,k) + vely_arr(i-1,j+1,k) );
        amrex::Real wsp   = std::sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<u'w'> = -K dudz
        // NOTE: tau_tot = rho * Cd * wsp^2, multiply by u/wsp to get tau_13
        amrex::Real stressx = -rho * mdata.Cd * wsp * velx;

        return stressx;
    }

    /**
     * Compute y-momentum flux for one surface face.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, vm_arr,
     * u_star_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_v_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/) const
    {
        amrex::Real rho   = myhalf * ( cons_arr(i,j-1,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );
        amrex::Real velx  = fourth * ( velx_arr(i,j  ,k) + velx_arr(i+1,j  ,k)
                                     + velx_arr(i,j-1,k) + velx_arr(i+1,j-1,k) );
        amrex::Real vely  = vely_arr(i,j,k);
        amrex::Real wsp   = std::sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<v'w'> = -K dvdz
        // NOTE: tau_tot = rho * Cd * wsp^2, multiply by v/wsp to get tau_13
        amrex::Real stressy = -rho * mdata.Cd * wsp * vely;

        return stressy;
    }

private:
    most_data mdata;
};


/**
 * RICO flux formulation
 */
struct rico_flux
{
    /**
     * Construct the RICO flux calculator.
     *
     * @param[in] l_theta_z0 surface reference potential temperature
     * @param[in] l_qsat_z0  surface saturation specific humidity
     */
    rico_flux (amrex::Real l_theta_z0, amrex::Real l_qsat_z0)
      : theta_z0{l_theta_z0}, qsat_z0{l_qsat_z0} {}


    /**
     * Compute moisture flux for one surface point.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, qvm_arr,
     * u_star_arr, q_star_arr, q_surf_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_q_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*qm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& q_star_arr,
                    const amrex::Array4<const amrex::Real>& /*t_surf_arr*/) const
    {
        amrex::Real rho   = cons_arr(i,j,k,Rho_comp);
        amrex::Real qstar = q_star_arr(i,j,0);
        amrex::Real q     = cons_arr(i, j, k, RhoQ1_comp) / rho;

        amrex::Real velx = myhalf * (velx_arr(i,j,k) + velx_arr(i+1,j,k));
        amrex::Real vely = myhalf * (vely_arr(i,j,k) + vely_arr(i,j+1,k));
        amrex::Real wsp = std::sqrt(velx*velx + vely*vely);

        // NOTE: this is rho*<q'w'> = -K dqdz
        amrex::Real moflux  = (std::abs(qstar) > eps) ? - rho * qstar * wsp * (q - qsat_z0) : zero;

        return moflux;
    }

    /**
     * Compute heat flux for one surface point.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, tm_arr,
     * u_star_arr, t_star_arr, t_surf_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_t_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*tm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& t_star_arr,
                    const amrex::Array4<const amrex::Real>& /*t_surf_arr*/) const
    {
        amrex::Real rho   = cons_arr(i,j,k,Rho_comp);
        amrex::Real tstar = t_star_arr(i, j, 0);
        amrex::Real theta = cons_arr(i, j, k, RhoTheta_comp) / rho;

        amrex::Real velx = myhalf * (velx_arr(i,j,k) + velx_arr(i+1,j,k));
        amrex::Real vely = myhalf * (vely_arr(i,j,k) + vely_arr(i,j+1,k));
        amrex::Real wsp = std::sqrt(velx*velx + vely*vely);

        // NOTE: this is rho*<T'w'> = -K dTdz
        amrex::Real moflux  = (std::abs(tstar) > eps) ? - rho * tstar * wsp * (theta - theta_z0) : zero;

        return moflux;
    }

    /**
     * Compute x-momentum flux for one surface face.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, um_arr,
     * u_star_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_u_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*um_arr*/,
                    const amrex::Array4<const amrex::Real>& u_star_arr) const
    {
        amrex::Real velx = velx_arr(i,j,k);
        amrex::Real vely = fourth * ( vely_arr(i  ,j,k) + vely_arr(i  ,j+1,k)
                                    + vely_arr(i-1,j,k) + vely_arr(i-1,j+1,k) );
        amrex::Real rho  = myhalf * ( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real ustar = myhalf * ( u_star_arr(i-1,j,0) + u_star_arr(i,j,0) );
        amrex::Real wsp   = std::sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<u'w'> = -K dudz
        amrex::Real stressx = -rho * ustar * wsp * velx;

        return stressx;
    }

    /**
     * Compute y-momentum flux for one surface face.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, vm_arr,
     * u_star_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_v_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vm_arr*/,
                    const amrex::Array4<const amrex::Real>& u_star_arr) const
    {
        amrex::Real velx  = fourth * ( velx_arr(i,j  ,k) + velx_arr(i+1,j  ,k)
                                     + velx_arr(i,j-1,k) + velx_arr(i+1,j-1,k) );
        amrex::Real vely  = vely_arr(i,j,k);
        amrex::Real rho   = myhalf * ( cons_arr(i,j-1,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real ustar   = myhalf * ( u_star_arr(i,j-1,0) + u_star_arr(i,j,0) );
        amrex::Real wsp     = std::sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<v'w'> = -K dvdz
        amrex::Real stressy = -rho * ustar * wsp * vely;

        return stressy;
    }

private:
#ifdef AMREX_USE_FLOAT
    const amrex::Real eps = amrex::Real(1e-8);
#else
    const amrex::Real eps = amrex::Real(1e-15);
#endif
    amrex::Real theta_z0 = amrex::Real(298.0);
    amrex::Real qsat_z0 = amrex::Real(0.0010);
};


/**
 * Rotate flux formulation
 */
struct rotate_flux
{
    /**
     * Construct the rotated-flux calculator.
     */
    rotate_flux () {}

    /**
     * Compute moisture flux for one surface point.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, qvm_arr,
     * u_star_arr, q_star_arr, q_surf_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_q_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& qvm_arr,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& q_star_arr,
                    const amrex::Array4<const amrex::Real>& q_surf_arr) const
    {
        // NOTE: this is the total stress
        amrex::Real qv_mean = qvm_arr(i,j,0);
        amrex::Real qv_surf = q_surf_arr(i,j,0);
        amrex::Real rho     = cons_arr(i,j,k,Rho_comp);
        amrex::Real qstar   = q_star_arr(i,j,0);
        amrex::Real moflux  = (std::abs(qstar) > eps) ? -rho * qstar * (qv_mean-qv_surf): zero;

        return moflux;
    }

    /**
     * Compute heat flux for one surface point.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, tm_arr,
     * u_star_arr, t_star_arr, t_surf_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_t_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& tm_arr,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& t_star_arr,
                    const amrex::Array4<const amrex::Real>& t_surf_arr) const
    {
        // NOTE: this is the total stress
        amrex::Real theta_mean = tm_arr(i,j,0);
        amrex::Real theta_surf = t_surf_arr(i,j,0);
        amrex::Real rho        = cons_arr(i,j,k,Rho_comp);
        amrex::Real tstar      = t_star_arr(i,j,0);
        amrex::Real moflux = (std::abs(tstar) > eps) ? -rho * tstar * (theta_mean-theta_surf) : zero;

        return moflux;
    }

    /**
     * Compute x-momentum flux for one surface face.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, um_arr,
     * u_star_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_u_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*um_arr*/,
                    const amrex::Array4<const amrex::Real>& u_star_arr) const
    {
        // NOTE: this is the total stress
        amrex::Real rho     = myhalf *( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i  ,j,k,Rho_comp) );
        amrex::Real ustar   = myhalf * ( u_star_arr(i-1,j,0) + u_star_arr(i,j,0) );
        amrex::Real stressx = -rho * ustar * ustar;

        return stressx;
    }

    /**
     * Compute y-momentum flux for one surface face.
     *
     * @par Calling sequence
     * i, j, k, cons_arr, velx_arr, vely_arr, umm_arr, vm_arr,
     * u_star_arr.
     */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_v_flux (const int& /*i*/,
                    const int& /*j*/,
                    const int& /*k*/,
                    const amrex::Array4<const amrex::Real>& /*cons_arr*/,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/) const
    {
        // NOTE: this is the total stress
        amrex::Real stressy = zero;

        return stressy;
    }

private:
#ifdef AMREX_USE_FLOAT
    const amrex::Real eps = amrex::Real(1e-8);
#else
    const amrex::Real eps = amrex::Real(1e-15);
#endif
};
#endif

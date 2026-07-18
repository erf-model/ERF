/**
 * @file ERF_FireUtils.H
 * @brief Utility functions for fire model atmospheric state extraction
 */

#ifndef ERF_FIRE_UTILS_H
#define ERF_FIRE_UTILS_H

#include <AMReX_MultiFab.h>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <ERF_MicrophysicsUtils.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_Constants.H>

/**
 * @brief Compute relative humidity from conserved atmospheric variables at k=0.
 *
 * Computes RH from pressure, temperature, and water vapor mixing ratio at 
 * the first atmospheric level (k=0). Uses the ERF microphysics saturation 
 * vapor pressure functions.
 *
 * @param[out] RH_out        Output relative humidity MultiFab [fraction 0-1]
 *                           Must be allocated with ncomp=1 on same grid as inputs
 * @param[in]  S_atm         Conservative state at k=0
 *                           Contains RhoTheta (comp 1) and RhoQ1 (comp 4)
 * @param[in]  geom          Geometry at level 0
 *
 * @note
 * - Output is clamped to [0, 1] range
 * - Assumes S_atm contains standard conserved variables (Rho, RhoTheta, RhoQ1)
 * - Temperature extracted from RhoTheta using standard EOS
 * - Pressure computed from RhoTheta and qv using getPgivenRTh
 */
inline
void compute_rh_from_conservative(amrex::MultiFab& RH_out,
                                  const amrex::MultiFab& S_atm,
                                  const amrex::Geometry& /*geom*/)
{
    using namespace amrex;

    // Precompute rdOcp = R_d / Cp_d for use in Exner function
    // Cp_d = 1004.5 J/(kg-K), R_d = 287.0 J/(kg-K)
    constexpr Real rdOcp = R_d / Cp_d;

    // Iterate over all boxes
    for (MFIter mfi(RH_out); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real>       rh_arr = RH_out.array(mfi);
        Array4<const Real> S_arr  = S_atm.const_array(mfi);

        // Kernel to compute RH at k=0
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(const IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = 0;  // First atmospheric level

            // Extract conserved variables at k=0
            Real rho     = S_arr(i, j, k, Rho_comp);           // Density
            Real rho_th  = S_arr(i, j, k, RhoTheta_comp);     // RhoTheta

            // Compute mixing ratios (guard against moisture not being present)
            Real qv = 0.0_rt;  // Water vapor mixing ratio [kg/kg]
            if (RhoQ1_comp < S_arr.nComp()) {
                Real rho_qv  = S_arr(i, j, k, RhoQ1_comp);    // RhoQ1 (water vapor)
                qv = rho_qv / rho;
            }

            // Compute potential temperature
            Real theta = rho_th / rho;  // Potential temperature [K]

            // Compute pressure from RhoTheta and qv using standard EOS
            Real p = getPgivenRTh(rho_th, qv);  // Pressure [Pa]
            
            // Convert pressure to mbar for saturation vapor pressure functions
            Real p_mbar = p * 0.01_rt;  // Pa to mbar (1 mbar = 100 Pa)

            // Compute temperature from potential temperature
            // T = theta * (p / p_0)^(R_d / c_p)
            // Use Exner function: pi = (p / p_0)^(R_d / c_p)
            Real exner = getExnergivenP(p, rdOcp);
            Real T_K = theta * exner;

            // Compute saturation vapor pressure [mbar]
            Real e_sat = erf_esatw(T_K, false);  // Use Flatau polynomial

            // Compute saturation mixing ratio [kg/kg]
            Real q_sat = erf_qsat_from_esat(e_sat, p_mbar);

            // Compute relative humidity as RH = qv / qsat
            Real RH = 0.0_rt;
            if (q_sat > 1.0e-10_rt) {
                RH = qv / q_sat;
            }

            // Clamp RH to valid range [0, 1]
            RH = amrex::min(1.0_rt, amrex::max(0.0_rt, RH));

            rh_arr(i, j, k) = RH;
        });
    }
}

#endif  // ERF_FIRE_UTILS_H

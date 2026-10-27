/**
 * @file ERF_FireUtils.H
 * @brief Utility functions for fire model atmospheric state extraction
 */

#ifndef ERF_FIRE_UTILS_H
#define ERF_FIRE_UTILS_H

#include <AMReX_MultiFab.H>
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

/**
 * @brief Air temperature [K] at k = 0 from the conserved state.
 *
 * T = theta * (p / p_0)^(R_d / c_p), with p from RhoTheta and qv as in
 * compute_rh_from_conservative. Fuel-moisture curves and ignition criteria take
 * an air temperature, not the potential temperature of the first cell, which at
 * 840 hPa is 14-16 K above it.
 *
 * @param[out] T_out  2D MultiFab (one cell in z) on the level's boxes
 * @param[in]  S_atm  Conserved state (Rho, RhoTheta, RhoQ1 when moist)
 */
inline
void compute_t_from_conservative(amrex::MultiFab& T_out,
                                 const amrex::MultiFab& S_atm)
{
    using namespace amrex;
    constexpr Real rdOcp = R_d / Cp_d;
    for (MFIter mfi(T_out); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real>       t_arr = T_out.array(mfi);
        Array4<const Real> S_arr = S_atm.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(const IntVect& iv) {
            const int i = iv[0];
            const int j = iv[1];
            const int k = 0;
            const Real rho    = S_arr(i, j, k, Rho_comp);
            const Real rho_th = S_arr(i, j, k, RhoTheta_comp);
            Real qv = 0.0_rt;
            if (RhoQ1_comp < S_arr.nComp()) {
                qv = S_arr(i, j, k, RhoQ1_comp) / rho;
            }
            const Real p     = getPgivenRTh(rho_th, qv);
            const Real exner = getExnergivenP(p, rdOcp);
            t_arr(i, j, k)   = (rho_th / rho) * exner;
        });
    }
}

/**
 * Copy the k = 0 plane of a cell-centred atmospheric field into a MultiFab on the
 * same boxes coarsened to one cell in z. A namespace-scope function because nvcc
 * rejects a device kernel inside a private member such as ERF::Advance.
 */
inline
void fire_copy_k0_plane(amrex::MultiFab& dst_2d, const amrex::MultiFab& src_3d)
{
    for (amrex::MFIter mfi(dst_2d); mfi.isValid(); ++mfi) {
        amrex::Array4<amrex::Real> d = dst_2d.array(mfi);
        amrex::Array4<const amrex::Real> s = src_3d.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (const amrex::IntVect& iv) {
            d(iv[0], iv[1], 0) = s(iv[0], iv[1], 0);
        });
    }
}

/**
 * @brief Copy a MultiFab between device and host-accessible memory, for an
 *        engine whose fire physics runs on the host
 *
 * The memory contract of ERF_FireEngine.H: the coupling's fields are in the
 * default arena (device memory on a GPU build) and are only valid during the
 * call that hands them over; an engine's returned fields must be readable
 * there too. dst and src have the same BoxArray and DistributionMapping;
 * ncomp components and nghost ghost cells are copied from component 0. dst
 * belongs to the caller (allocate a host copy once with
 * amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()), or The_Managed_Arena(),
 * and reuse it). The copy completes before the call returns: the GPU stream is
 * synchronised, so the host may read dst, or the device read it, right away.
 * On a CPU build it is a plain copy.
 */
inline
void fire_host_copy(amrex::MultiFab& dst, const amrex::MultiFab& src, int ncomp = 1, int nghost = 0)
{
    AMREX_ALWAYS_ASSERT(dst.boxArray() == src.boxArray() && dst.DistributionMap() == src.DistributionMap());
    amrex::MultiFab::Copy(dst, src, 0, 0, ncomp, nghost);
    amrex::Gpu::streamSynchronize();
}

#endif  // ERF_FIRE_UTILS_H

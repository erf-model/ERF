/**
 * \file ERF_EBSlopes.H
 * \brief Provides least-squares slope reconstruction helpers near EB cut cells.
 */
#ifndef ERF_EBSLOPES_H_
#define ERF_EBSLOPES_H_

#include <AMReX_EBCellFlag.H>
#include <AMReX_EB_Slopes_K.H>
#include <ERF_IndexDefines.H>

/**
 * \brief Compute least-squares slopes using EB Dirichlet data.
 *
 * @param[in] dx Mesh spacing in x.
 * @param[in] dy Mesh spacing in y.
 * @param[in] dz Mesh spacing in z.
 * @param[in] i Cell index in x.
 * @param[in] j Cell index in y.
 * @param[in] k Cell index in z.
 * @param[in] bcent_eb EB boundary centroid.
 * @param[in] state_eb Dirichlet state at the EB boundary.
 * @param[in] state Cell-centered state array.
 * @param[in] ccent Cell centroid coordinates.
 * @param[in] flag Cell-centered EB flags and connectivity.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
erf_calc_slopes_eb_Dirichlet (
                    amrex::Real dx, amrex::Real dy, amrex::Real dz,
                    int i, int j, int k,
                    amrex::RealVect const& bcent_eb,
                    amrex::Real const state_eb,
                    amrex::Array4<amrex::Real const> const& state,
                    amrex::Array4<amrex::Real const> const& ccent,
                    amrex::Array4<amrex::EBCellFlag const> const& flag)
{
    amrex::EBSlopeNormalEquation normal_eq;

    // Prevent loop unrolling to reduce register pressure.
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for (int kk = -1; kk <= 1; ++kk) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for (int jj = -1; jj <= 1; ++jj) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for (int ii = -1; ii <= 1; ++ii) {
        if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii == 0 && jj == 0 && kk == 0)) {
            amrex::Real offset_x = (amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0) - bcent_eb[0]) * dx;
            amrex::Real offset_y = (amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1) - bcent_eb[1]) * dy;
            amrex::Real offset_z = (amrex::Real(kk) + ccent(i+ii,j+jj,k+kk,2) - bcent_eb[2]) * dz;
            amrex::Real du = state(i+ii,j+jj,k+kk) - state_eb;
            normal_eq.add(offset_x, offset_y, offset_z, du);
        }
    }}}

    return normal_eq.solve();
}

/**
 * \brief Compute least-squares slopes from staggered data using EB Dirichlet data.
 *
 * @tparam GridQuery Grid index for the query location.
 * @tparam GridData Grid index for the data locations.
 * @param[in] dx Mesh spacing in x.
 * @param[in] dy Mesh spacing in y.
 * @param[in] dz Mesh spacing in z.
 * @param[in] i Cell index in x.
 * @param[in] j Cell index in y.
 * @param[in] k Cell index in z.
 * @param[in] bcent_eb EB boundary centroid.
 * @param[in] state_eb Dirichlet state at the EB boundary.
 * @param[in] state State array on the data grid.
 * @param[in] ccent Cell centroid coordinates.
 * @param[in] flag EB flags on the data grid.
 */
template <int GridQuery, int GridData>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
erf_calc_slopes_eb_Dirichlet_staggered (
                    amrex::Real dx, amrex::Real dy, amrex::Real dz,
                    int i, int j, int k,
                    amrex::RealVect const& bcent_eb,
                    amrex::Real const state_eb,
                    amrex::Array4<amrex::Real const> const& state,
                    amrex::Array4<amrex::Real const> const& ccent,
                    amrex::Array4<amrex::EBCellFlag const> const& flag)
{
    static_assert(GridQuery == Vars::xvel || GridQuery == Vars::yvel ||
                  GridQuery == Vars::zvel);
    static_assert(GridData == Vars::xvel || GridData == Vars::yvel ||
                  GridData == Vars::zvel);

    // Fitting stencil
    constexpr int ii_lo = (GridQuery == Vars::xvel) ? -2 : -1;
    constexpr int jj_lo = (GridQuery == Vars::yvel) ? -2 : -1;
    constexpr int kk_lo = (GridQuery == Vars::zvel) ? -2 : -1;
    constexpr int ii_hi = (GridData  == Vars::xvel) ?  2 :  1;
    constexpr int jj_hi = (GridData  == Vars::yvel) ?  2 :  1;
    constexpr int kk_hi = (GridData  == Vars::zvel) ?  2 :  1;

    // Set bias in the index space between the two staggered grids
    amrex::RealVect bias{zero,zero,zero};
    bias[GridQuery-1] = myhalf;
    bias[GridData-1] = -myhalf;

    amrex::EBSlopeNormalEquation normal_eq;

    // Prevent loop unrolling to reduce register pressure.
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for (int kk = kk_lo; kk <= kk_hi; ++kk) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for (int jj = jj_lo; jj <= jj_hi; ++jj) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for (int ii = ii_lo; ii <= ii_hi; ++ii) {
        if (!flag(i+ii,j+jj,k+kk).isCovered()) {
            amrex::Real offset_x = (amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0) + bias[0] - bcent_eb[0]) * dx;
            amrex::Real offset_y = (amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1) + bias[1] - bcent_eb[1]) * dy;
            amrex::Real offset_z = (amrex::Real(kk) + ccent(i+ii,j+jj,k+kk,2) + bias[2] - bcent_eb[2]) * dz;
            amrex::Real du = state(i+ii,j+jj,k+kk) - state_eb;
            normal_eq.add(offset_x, offset_y, offset_z, du);
        }
    }}}

    return normal_eq.solve();
}

/**
 * \brief Compute centered least-squares slopes from cell data to a staggered query grid.
 *
 * @param[in] igrid_query Grid index for the staggered query location.
 * @param[in] igrid_data Grid index for the source data locations.
 * @param[in] dx Mesh spacing in x.
 * @param[in] dy Mesh spacing in y.
 * @param[in] dz Mesh spacing in z.
 * @param[in] i Cell index in x.
 * @param[in] j Cell index in y.
 * @param[in] k Cell index in z.
 * @param[in] state State array on the data grid.
 * @param[in] ccent Cell centroid coordinates.
 * @param[in] flag EB flags on the data grid.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
erf_calc_slopes_eb_staggered (int igrid_query,
                    int igrid_data,
                    amrex::Real dx, amrex::Real dy, amrex::Real dz,
                    int i, int j, int k,
                    amrex::Array4<amrex::Real const> const& state,
                    amrex::Array4<amrex::Real const> const& ccent,
                    amrex::Array4<amrex::EBCellFlag const> const& flag)
{
    AMREX_ASSERT((igrid_query == Vars::xvel || igrid_query == Vars::yvel || igrid_query == Vars::zvel) && igrid_data == Vars::cons);

    // Fitting stencil

    int ii_lo = -1;
    int jj_lo = -1;
    int kk_lo = -1;
    int ii_hi = 1;
    int jj_hi = 1;
    int kk_hi = 1;

    if (igrid_query == Vars::xvel) {
        ii_lo = -2;
    } else if (igrid_query == Vars::yvel) {
        jj_lo = -2;
    } else if (igrid_query == Vars::zvel) {
        kk_lo = -2;
    }

    if (igrid_data == Vars::xvel) {
        ii_hi = 2;
    } else if (igrid_data == Vars::yvel) {
        jj_hi = 2;
    } else if (igrid_data == Vars::zvel) {
        kk_hi = 2;
    }

    // Set bias in the index space between the two staggered grids
    amrex::RealVect bias{zero,zero,zero};
    bias[igrid_query-1] = myhalf;

    amrex::EBSlopeNormalEquation normal_eq;

    // Prevent loop unrolling to reduce register pressure.
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for (int kk = kk_lo; kk <= kk_hi; ++kk) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for (int jj = jj_lo; jj <= jj_hi; ++jj) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for (int ii = ii_lo; ii <= ii_hi; ++ii) {
        if (!flag(i+ii,j+jj,k+kk).isCovered()) {
            amrex::Real offset_x = (amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0)
                                    + bias[0] - ccent(i,j,k,0)) * dx;
            amrex::Real offset_y = (amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1)
                                    + bias[1] - ccent(i,j,k,1)) * dy;
            amrex::Real offset_z = (amrex::Real(kk) + ccent(i+ii,j+jj,k+kk,2)
                                    + bias[2] - ccent(i,j,k,2)) * dz;
            amrex::Real du = state(i+ii,j+jj,k+kk) - state(i,j,k);
            normal_eq.add(offset_x, offset_y, offset_z, du);
        }
    }}}

    return normal_eq.solve();
}

/**
 * \brief Compute upwind-biased least-squares slopes for staggered EB advection.
 *
 * @param[in] igrid_query Grid index for the staggered query location.
 * @param[in] igrid_data Grid index for the source data locations.
 * @param[in] dx Mesh spacing in x.
 * @param[in] dy Mesh spacing in y.
 * @param[in] dz Mesh spacing in z.
 * @param[in] i Cell index in x.
 * @param[in] j Cell index in y.
 * @param[in] k Cell index in z.
 * @param[in] vel_arr Velocity array used to choose the upwind stencil.
 * @param[in] state State array on the data grid.
 * @param[in] ccent Cell centroid coordinates.
 * @param[in] flag EB flags on the data grid.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
erf_calc_slopes_eb_staggered_upwind (int igrid_query,
                    [[maybe_unused]] int igrid_data,
                    amrex::Real dx, amrex::Real dy, amrex::Real dz,
                    int i, int j, int k,
                    const amrex::Array4<const amrex::Real>& vel_arr,
                    amrex::Array4<amrex::Real const> const& state,
                    amrex::Array4<amrex::Real const> const& ccent,
                    amrex::Array4<amrex::EBCellFlag const> const& flag)
{
    AMREX_ASSERT((igrid_query == Vars::xvel || igrid_query == Vars::yvel || igrid_query == Vars::zvel) && igrid_data == Vars::cons);

    // Fitting stencil

    int ii_lo = -1;
    int jj_lo = -1;
    int kk_lo = -1;
    int ii_hi = 1;
    int jj_hi = 1;
    int kk_hi = 1;

    if (igrid_query == Vars::xvel && vel_arr(i,j,k)>zero) {
        ii_lo = -2;
        ii_hi = 0;
    }
    if (igrid_query == Vars::yvel && vel_arr(i,j,k)>zero) {
        jj_lo = -2;
        jj_hi = 0;
    }
    if (igrid_query == Vars::zvel && vel_arr(i,j,k)>zero) {
        kk_lo = -2;
        kk_hi = 0;
    }

    // Set bias in the index space between the two staggered grids
    amrex::RealVect bias{zero,zero,zero};
    bias[igrid_query-1] = myhalf;

    amrex::EBSlopeNormalEquation normal_eq;

    // Prevent loop unrolling to reduce register pressure.
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for (int kk = kk_lo; kk <= kk_hi; ++kk) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for (int jj = jj_lo; jj <= jj_hi; ++jj) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for (int ii = ii_lo; ii <= ii_hi; ++ii) {
        if (!flag(i+ii,j+jj,k+kk).isCovered() && !(ii == 0 && jj == 0 && kk == 0)) {
            amrex::Real offset_x = (amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0)
                                    + bias[0] - ccent(i,j,k,0)) * dx;
            amrex::Real offset_y = (amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1)
                                    + bias[1] - ccent(i,j,k,1)) * dy;
            amrex::Real offset_z = (amrex::Real(kk) + ccent(i+ii,j+jj,k+kk,2)
                                    + bias[2] - ccent(i,j,k,2)) * dz;
            amrex::Real du = state(i+ii,j+jj,k+kk) - state(i,j,k);
            normal_eq.add(offset_x, offset_y, offset_z, du);
        }
    }}}

    return normal_eq.solve();
}

#endif

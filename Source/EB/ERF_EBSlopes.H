#ifndef ERF_EBSLOPES_H_
#define ERF_EBSLOPES_H_

#include <AMReX_EBCellFlag.H>
#include <ERF_IndexDefines.H>

/**
 * Function for computing the slopes with Dirichlet data
 * This function works for only face-centered grids
 *
 * @param[in] igrid_query index of grids where the query point is located (0=CC, 1=U, 2=V, 3=W)
 * @param[in] igrid_data index of grids where data points are located (0=CC, 1=U, 2=V, 3=W)
 * @param[in] bcent_eb centroid of the EB face
 * @param[in] state_eb state at the EB face (Dirichlet data)
 * @param[in] state array of state variable
 * @param[in] ccent array of cell centroid coordinates
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
    constexpr int dim_a = 27;

    int ii_lo = -1; int jj_lo = -1; int kk_lo = -1;
    int ii_hi =  1; int jj_hi =  1; int kk_hi =  1;

    amrex::Real A[dim_a][AMREX_SPACEDIM];

    // Array of distances from the query point
    int ll=0;
    for(int kk(kk_lo); kk<=kk_hi; kk++) {
    for(int jj(jj_lo); jj<=jj_hi; jj++) {
    for(int ii(ii_lo); ii<=ii_hi; ii++) {

        if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0 && kk==0)) {
            A[ll][0] = ( amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0) - bcent_eb[0] ) * dx;
            A[ll][1] = ( amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1) - bcent_eb[1] ) * dy;
            A[ll][2] = ( amrex::Real(kk) + ccent(i+ii,j+jj,k+kk,2) - bcent_eb[2] ) * dz;
        } else {
            A[ll][0] = amrex::Real(0.0);
            A[ll][1] = amrex::Real(0.0);
            A[ll][2] = amrex::Real(0.0);
        }
        ll++;
    }}}

    //
    // Calculate the slopes given the matrix A (See amrex_calc_slopes_eb_given_A)
    //
    amrex::Real du[dim_a];

    ll=0;
    for(int kk(kk_lo); kk<=kk_hi; kk++) {
    for(int jj(jj_lo); jj<=jj_hi; jj++) {
    for(int ii(ii_lo); ii<=ii_hi; ii++) {

        if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0 && kk==0)) {
            du[ll] = state(i+ii,j+jj,k+kk) - state_eb;
        } else {
            du[ll] = amrex::Real(0.0);
        }
        ll++;
    }}}

    amrex::Real AtA[AMREX_SPACEDIM][AMREX_SPACEDIM];
    amrex::Real Atb[AMREX_SPACEDIM];

    for(int jj(0); jj<AMREX_SPACEDIM; ++jj){
        for(int ii(0); ii<AMREX_SPACEDIM; ++ii){ // NOLINT(modernize-loop-convert)
            AtA[jj][ii] = amrex::Real(0.0);
        }
        Atb[jj] = amrex::Real(0.0);
    }

    for(int lc(0); lc < dim_a; ++lc)
    {
        AtA[0][0] += A[lc][0]* A[lc][0];
        AtA[0][1] += A[lc][0]* A[lc][1];
        AtA[0][2] += A[lc][0]* A[lc][2];
        AtA[1][1] += A[lc][1]* A[lc][1];
        AtA[1][2] += A[lc][1]* A[lc][2];
        AtA[2][2] += A[lc][2]* A[lc][2];

        Atb[0] += A[lc][0]*du[lc];
        Atb[1] += A[lc][1]*du[lc];
        Atb[2] += A[lc][2]*du[lc];
    }

    // Fill in symmetric
    AtA[1][0] = AtA[0][1];
    AtA[2][0] = AtA[0][2];
    AtA[2][1] = AtA[1][2];

    amrex::Real detAtA =
      AtA[0][0]*(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[2][1]) -
      AtA[0][1]*(AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
      AtA[0][2]*(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

    amrex::Real detAtA_x =
        Atb[0]   *(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
        AtA[0][1]*(Atb[1] *  AtA[2][2] - AtA[1][2]*Atb[2]   ) +
        AtA[0][2]*(Atb[1] *  AtA[2][1] - AtA[1][1]*Atb[2]   );

    // Slope at centroid of (i,j,k)
    amrex::Real xslope = detAtA_x / detAtA;

    amrex::Real detAtA_y =
        AtA[0][0]*(Atb[1]  * AtA[2][2] - AtA[1][2]*Atb[2]   ) -
        Atb[0] *  (AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
        AtA[0][2]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]);

    // Slope at centroid of (i,j,k)
    amrex::Real yslope = detAtA_y / detAtA;

    amrex::Real detAtA_z =
        AtA[0][0]*(AtA[1][1]*Atb[2]    - Atb[1]   *AtA[1][2]) -
        AtA[0][1]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]) +
        Atb[0]   *(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

    // Slope at centroid of (i,j,k)
    amrex::Real zslope = detAtA_z / detAtA;

    return {xslope,yslope,zslope};
}

/**
 * Function for computing the slopes with Dirichlet data in staggered grids
 * This function works for only face-centered grids
 *
 * @param[in] igrid_query index of grids where the query point is located (0=CC, 1=U, 2=V, 3=W)
 * @param[in] igrid_data index of grids where data points are located (0=CC, 1=U, 2=V, 3=W)
 * @param[in] bcent_eb centroid of the EB face
 * @param[in] state_eb state at the EB face (Dirichlet data)
 * @param[in] state array of state variable
 * @param[in] ccent array of cell centroid coordinates
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
erf_calc_slopes_eb_Dirichlet_staggered (int igrid_query,
                    int igrid_data,
                    amrex::Real dx, amrex::Real dy, amrex::Real dz,
                    int i, int j, int k,
                    amrex::RealVect const& bcent_eb,
                    amrex::Real const state_eb,
                    amrex::Array4<amrex::Real const> const& state,
                    amrex::Array4<amrex::Real const> const& ccent,
                    amrex::Array4<amrex::EBCellFlag const> const& flag)
{
    // Fitting stencil

    constexpr int dim_a = 48;

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
    amrex::RealVect bias{0.,0.,0.};
    bias[igrid_query-1] = 0.5;
    bias[igrid_data-1] = -0.5;

    amrex::Real A[dim_a][AMREX_SPACEDIM];

    // Array of distances from the query point
    int ll=0;
    for(int kk(kk_lo); kk<=kk_hi; kk++) {
    for(int jj(jj_lo); jj<=jj_hi; jj++) {
    for(int ii(ii_lo); ii<=ii_hi; ii++) {

        if (!flag(i+ii,j+jj,k+kk).isCovered()) {
            A[ll][0] = ( amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0) + bias[0] - bcent_eb[0] ) * dx;
            A[ll][1] = ( amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1) + bias[1] - bcent_eb[1] ) * dy;
            A[ll][2] = ( amrex::Real(kk) + ccent(i+ii,j+jj,k+kk,2) + bias[2] - bcent_eb[2] ) * dz;
        } else {
            A[ll][0] = amrex::Real(0.0);
            A[ll][1] = amrex::Real(0.0);
            A[ll][2] = amrex::Real(0.0);
        }
        ll++;
    }}}

    //
    // Calculate the slopes given the matrix A (See amrex_calc_slopes_eb_given_A)
    //
    amrex::Real du[dim_a];

    ll=0;
    for(int kk(kk_lo); kk<=kk_hi; kk++) {
    for(int jj(jj_lo); jj<=jj_hi; jj++) {
    for(int ii(ii_lo); ii<=ii_hi; ii++) {

        if (!flag(i+ii,j+jj,k+kk).isCovered()) {
            du[ll] = state(i+ii,j+jj,k+kk) - state_eb;
        } else {
            du[ll] = amrex::Real(0.0);
        }
        ll++;
    }}}

    amrex::Real AtA[AMREX_SPACEDIM][AMREX_SPACEDIM];
    amrex::Real Atb[AMREX_SPACEDIM];

    for(int jj(0); jj<AMREX_SPACEDIM; ++jj){
        for(int ii(0); ii<AMREX_SPACEDIM; ++ii){ // NOLINT(modernize-loop-convert)
            AtA[jj][ii] = amrex::Real(0.0);
        }
        Atb[jj] = amrex::Real(0.0);
    }

    for(int lc(0); lc < dim_a; ++lc)
    {
        AtA[0][0] += A[lc][0]* A[lc][0];
        AtA[0][1] += A[lc][0]* A[lc][1];
        AtA[0][2] += A[lc][0]* A[lc][2];
        AtA[1][1] += A[lc][1]* A[lc][1];
        AtA[1][2] += A[lc][1]* A[lc][2];
        AtA[2][2] += A[lc][2]* A[lc][2];

        Atb[0] += A[lc][0]*du[lc];
        Atb[1] += A[lc][1]*du[lc];
        Atb[2] += A[lc][2]*du[lc];
    }

    // Fill in symmetric
    AtA[1][0] = AtA[0][1];
    AtA[2][0] = AtA[0][2];
    AtA[2][1] = AtA[1][2];

    amrex::Real detAtA =
      AtA[0][0]*(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[2][1]) -
      AtA[0][1]*(AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
      AtA[0][2]*(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

    amrex::Real detAtA_x =
        Atb[0]   *(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
        AtA[0][1]*(Atb[1] *  AtA[2][2] - AtA[1][2]*Atb[2]   ) +
        AtA[0][2]*(Atb[1] *  AtA[2][1] - AtA[1][1]*Atb[2]   );

    // Slope at centroid of (i,j,k)
    amrex::Real xslope = detAtA_x / detAtA;

    amrex::Real detAtA_y =
        AtA[0][0]*(Atb[1]  * AtA[2][2] - AtA[1][2]*Atb[2]   ) -
        Atb[0] *  (AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
        AtA[0][2]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]);

    // Slope at centroid of (i,j,k)
    amrex::Real yslope = detAtA_y / detAtA;

    amrex::Real detAtA_z =
        AtA[0][0]*(AtA[1][1]*Atb[2]    - Atb[1]   *AtA[1][2]) -
        AtA[0][1]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]) +
        Atb[0]   *(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

    // Slope at centroid of (i,j,k)
    amrex::Real zslope = detAtA_z / detAtA;

    return {xslope,yslope,zslope};
}

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

    constexpr int dim_a = 36;

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
    amrex::RealVect bias{0.,0.,0.};
    bias[igrid_query-1] = 0.5;

    amrex::Real A[dim_a][AMREX_SPACEDIM];

    // Array of distances from the query point
    int ll=0;
    for(int kk(kk_lo); kk<=kk_hi; kk++) {
    for(int jj(jj_lo); jj<=jj_hi; jj++) {
    for(int ii(ii_lo); ii<=ii_hi; ii++) {

        // if (!flag(i+ii,j+jj,k+kk).isCovered() && !(ii==0 && jj==0 && kk==0)) {
        if (!flag(i+ii,j+jj,k+kk).isCovered()) {
            A[ll][0] = ( amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0) + bias[0] - ccent(i,j,k,0) ) * dx;
            A[ll][1] = ( amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1) + bias[1] - ccent(i,j,k,1) ) * dy;
            A[ll][2] = ( amrex::Real(kk) + ccent(i+ii,j+jj,k+kk,2) + bias[2] - ccent(i,j,k,2) ) * dz;
        } else {
            A[ll][0] = amrex::Real(0.0);
            A[ll][1] = amrex::Real(0.0);
            A[ll][2] = amrex::Real(0.0);
        }
        ll++;
    }}}

    //
    // Calculate the slopes given the matrix A (See amrex_calc_slopes_eb_given_A)
    //
    amrex::Real du[dim_a];

    ll=0;
    for(int kk(kk_lo); kk<=kk_hi; kk++) {
    for(int jj(jj_lo); jj<=jj_hi; jj++) {
    for(int ii(ii_lo); ii<=ii_hi; ii++) {

        // if (!flag(i+ii,j+jj,k+kk).isCovered() && !(ii==0 && jj==0 && kk==0)) {
        if (!flag(i+ii,j+jj,k+kk).isCovered()) {
            du[ll] = state(i+ii,j+jj,k+kk) - state(i,j,k);
        } else {
            du[ll] = amrex::Real(0.0);
        }
        ll++;
    }}}

    amrex::Real AtA[AMREX_SPACEDIM][AMREX_SPACEDIM];
    amrex::Real Atb[AMREX_SPACEDIM];

    for(int jj(0); jj<AMREX_SPACEDIM; ++jj){
        for(int ii(0); ii<AMREX_SPACEDIM; ++ii){ // NOLINT(modernize-loop-convert)
            AtA[jj][ii] = amrex::Real(0.0);
        }
        Atb[jj] = amrex::Real(0.0);
    }

    for(int lc(0); lc < dim_a; ++lc)
    {
        AtA[0][0] += A[lc][0]* A[lc][0];
        AtA[0][1] += A[lc][0]* A[lc][1];
        AtA[0][2] += A[lc][0]* A[lc][2];
        AtA[1][1] += A[lc][1]* A[lc][1];
        AtA[1][2] += A[lc][1]* A[lc][2];
        AtA[2][2] += A[lc][2]* A[lc][2];

        Atb[0] += A[lc][0]*du[lc];
        Atb[1] += A[lc][1]*du[lc];
        Atb[2] += A[lc][2]*du[lc];
    }

    // Fill in symmetric
    AtA[1][0] = AtA[0][1];
    AtA[2][0] = AtA[0][2];
    AtA[2][1] = AtA[1][2];

    amrex::Real detAtA =
      AtA[0][0]*(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[2][1]) -
      AtA[0][1]*(AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
      AtA[0][2]*(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

    amrex::Real detAtA_x =
        Atb[0]   *(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
        AtA[0][1]*(Atb[1] *  AtA[2][2] - AtA[1][2]*Atb[2]   ) +
        AtA[0][2]*(Atb[1] *  AtA[2][1] - AtA[1][1]*Atb[2]   );

    // Slope at centroid of (i,j,k)
    amrex::Real xslope = detAtA_x / detAtA;

    amrex::Real detAtA_y =
        AtA[0][0]*(Atb[1]  * AtA[2][2] - AtA[1][2]*Atb[2]   ) -
        Atb[0] *  (AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
        AtA[0][2]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]);

    // Slope at centroid of (i,j,k)
    amrex::Real yslope = detAtA_y / detAtA;

    amrex::Real detAtA_z =
        AtA[0][0]*(AtA[1][1]*Atb[2]    - Atb[1]   *AtA[1][2]) -
        AtA[0][1]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]) +
        Atb[0]   *(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

    // Slope at centroid of (i,j,k)
    amrex::Real zslope = detAtA_z / detAtA;

    return {xslope,yslope,zslope};
}

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

    if (igrid_query == Vars::xvel && vel_arr(i,j,k)>0.0) {
        ii_lo = -2;
        ii_hi = 0;
    }
    if (igrid_query == Vars::yvel && vel_arr(i,j,k)>0.0) {
        jj_lo = -2;
        jj_hi = 0;
    }
    if (igrid_query == Vars::zvel && vel_arr(i,j,k)>0.0) {
        kk_lo = -2;
        kk_hi = 0;
    }

    // Set bias in the index space between the two staggered grids
    amrex::RealVect bias{0.,0.,0.};
    bias[igrid_query-1] = 0.5;

    constexpr int dim_a = 27;

    amrex::Real A[dim_a][AMREX_SPACEDIM];

    // Array of distances from the query point
    int ll=0;
    for(int kk(kk_lo); kk<=kk_hi; kk++) {
    for(int jj(jj_lo); jj<=jj_hi; jj++) {
    for(int ii(ii_lo); ii<=ii_hi; ii++) {

        if (!flag(i+ii,j+jj,k+kk).isCovered() && !(ii==0 && jj==0 && kk==0)) {
            A[ll][0] = ( amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0) + bias[0] - ccent(i,j,k,0) ) * dx;
            A[ll][1] = ( amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1) + bias[1] - ccent(i,j,k,1) ) * dy;
            A[ll][2] = ( amrex::Real(kk) + ccent(i+ii,j+jj,k+kk,2) + bias[2] - ccent(i,j,k,2) ) * dz;
        } else {
            A[ll][0] = amrex::Real(0.0);
            A[ll][1] = amrex::Real(0.0);
            A[ll][2] = amrex::Real(0.0);
        }
        ll++;
    }}}

    //
    // Calculate the slopes given the matrix A (See amrex_calc_slopes_eb_given_A)
    //
    amrex::Real du[dim_a];

    ll=0;
    for(int kk(kk_lo); kk<=kk_hi; kk++) {
    for(int jj(jj_lo); jj<=jj_hi; jj++) {
    for(int ii(ii_lo); ii<=ii_hi; ii++) {

        if (!flag(i+ii,j+jj,k+kk).isCovered() && !(ii==0 && jj==0 && kk==0)) {
            du[ll] = state(i+ii,j+jj,k+kk) - state(i,j,k);
        } else {
            du[ll] = amrex::Real(0.0);
        }
        ll++;
    }}}

    amrex::Real AtA[AMREX_SPACEDIM][AMREX_SPACEDIM];
    amrex::Real Atb[AMREX_SPACEDIM];

    for(int jj(0); jj<AMREX_SPACEDIM; ++jj){
        for(int ii(0); ii<AMREX_SPACEDIM; ++ii){ // NOLINT(modernize-loop-convert)
            AtA[jj][ii] = amrex::Real(0.0);
        }
        Atb[jj] = amrex::Real(0.0);
    }

    for(int lc(0); lc < dim_a; ++lc)
    {
        AtA[0][0] += A[lc][0]* A[lc][0];
        AtA[0][1] += A[lc][0]* A[lc][1];
        AtA[0][2] += A[lc][0]* A[lc][2];
        AtA[1][1] += A[lc][1]* A[lc][1];
        AtA[1][2] += A[lc][1]* A[lc][2];
        AtA[2][2] += A[lc][2]* A[lc][2];

        Atb[0] += A[lc][0]*du[lc];
        Atb[1] += A[lc][1]*du[lc];
        Atb[2] += A[lc][2]*du[lc];
    }

    // Fill in symmetric
    AtA[1][0] = AtA[0][1];
    AtA[2][0] = AtA[0][2];
    AtA[2][1] = AtA[1][2];

    amrex::Real detAtA =
      AtA[0][0]*(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[2][1]) -
      AtA[0][1]*(AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
      AtA[0][2]*(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

    amrex::Real detAtA_x =
        Atb[0]   *(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
        AtA[0][1]*(Atb[1] *  AtA[2][2] - AtA[1][2]*Atb[2]   ) +
        AtA[0][2]*(Atb[1] *  AtA[2][1] - AtA[1][1]*Atb[2]   );

    // Slope at centroid of (i,j,k)
    amrex::Real xslope = detAtA_x / detAtA;

    amrex::Real detAtA_y =
        AtA[0][0]*(Atb[1]  * AtA[2][2] - AtA[1][2]*Atb[2]   ) -
        Atb[0] *  (AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
        AtA[0][2]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]);

    // Slope at centroid of (i,j,k)
    amrex::Real yslope = detAtA_y / detAtA;

    amrex::Real detAtA_z =
        AtA[0][0]*(AtA[1][1]*Atb[2]    - Atb[1]   *AtA[1][2]) -
        AtA[0][1]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]) +
        Atb[0]   *(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

    // Slope at centroid of (i,j,k)
    amrex::Real zslope = detAtA_z / detAtA;

    return {xslope,yslope,zslope};
}

#endif

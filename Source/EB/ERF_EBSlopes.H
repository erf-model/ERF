#ifndef ERF_EBSLOPES_H_
#define ERF_EBSLOPES_H_

#include <AMReX_EBCellFlag.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
erf_calc_slopes_eb (const int idim,
                    int i, int j, int k, 
                    amrex::RealVect const& bcent_eb,
                    amrex::Real const state_eb,
                    amrex::Array4<amrex::Real const> const& state,
                    amrex::Array4<amrex::Real const> const& ccent,
                    amrex::Array4<amrex::EBCellFlag const> const& flag)
{

    int ii_hi = (idim==0) ? 2: 1;
    int jj_hi = (idim==1) ? 2: 1;
    int kk_hi = (idim==2) ? 2: 1;

    constexpr int dim_a = 27;
    amrex::Real A[dim_a][AMREX_SPACEDIM];

    // Array of distances from the query point
    int lc=0;
    for(int kk(-1); kk<=kk_hi; kk++) {
    for(int jj(-1); jj<=jj_hi; jj++) {
    for(int ii(-1); ii<=ii_hi; ii++) {

        if (!flag(i+ii,j+jj,k+kk).isCovered() ) {
            A[lc][0] = amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0) - bcent_eb[0];
            A[lc][1] = amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1) - bcent_eb[1];
            A[lc][2] = amrex::Real(kk) + ccent(i+ii,j+jj,k+kk,2) - bcent_eb[2];
        } else {
            A[lc][0] = amrex::Real(0.0);
            A[lc][1] = amrex::Real(0.0);
            A[lc][2] = amrex::Real(0.0);
        }
        lc++;
    }}}

    //
    // Calculate the slopes given the matrix A (See amrex_calc_slopes_eb_given_A)
    //
    amrex::Real du[dim_a];

    int ll=0;
    for(int kk(-1); kk<=kk_hi; kk++) {
    for(int jj(-1); jj<=jj_hi; jj++) {
    for(int ii(-1); ii<=ii_hi; ii++) {

        if (!flag(i+ii,j+jj,k+kk).isCovered() ) {
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
#endif
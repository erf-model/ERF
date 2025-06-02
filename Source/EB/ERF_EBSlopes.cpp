
#include <AMReX_Config.H>
#include <AMReX_Geometry.H>
#include <AMReX_EBCellFlag.H>

#include <ERF.H>
#include <ERF_EB.H>
#include <ERF_EBSlopes.H>

using namespace amrex;

// erf_calc_slopes_eb calculates the slope using a
// least squares linear fit to the 26 nearest neighbors,
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
GpuArray<Real,AMREX_SPACEDIM>
erf_calc_slopes_eb (int i, int j, int k, 
                    const Real ccent_x,
                    const Real ccent_y,
                    const Real ccent_z,
                    Array4<Real const> const& state,
                    Array4<Real const> const& ccent,
                    Array4<EBCellFlag const> const& flag)
{
    constexpr int dim_a = 27;
    Real A[dim_a][AMREX_SPACEDIM];

    int lc=0;
    for(int kk(-1); kk<=1; kk++) {
    for(int jj(-1); jj<=1; jj++) {
    for(int ii(-1); ii<=1; ii++)
    {
        if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0 && kk==0))
        {
            A[lc][0] = Real(ii) + ccent(i+ii,j+jj,k+kk,0) - ccent_x;
            A[lc][1] = Real(jj) + ccent(i+ii,j+jj,k+kk,1) - ccent_y;
            A[lc][2] = Real(kk) + ccent(i+ii,j+jj,k+kk,2) - ccent_z;
        } else {
            A[lc][0] = Real(0.0);
            A[lc][1] = Real(0.0);
            A[lc][2] = Real(0.0);
        }
        lc++;
    }}}

    // Calculate the slopes given the matrix A

    Real xslope = 0.0;
    Real yslope = 0.0;
    Real zslope = 0.0;

    Real du[dim_a];

    int ll=0;
    for(int kk(-1); kk<=1; kk++) {
    for(int jj(-1); jj<=1; jj++){
    for(int ii(-1); ii<=1; ii++){
        if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0 && kk==0))
        {
            du[ll] = state(i+ii,j+jj,k+kk) - state(i,j,k);
        } else {
            du[ll] = Real(0.0);
        }
        ll++;
    }}}


    return {xslope,yslope,zslope};    
}
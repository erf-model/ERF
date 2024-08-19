#ifndef _UTILS_THETAV_H_
#define _UTILS_THETAV_H_
/*
 * Calculate virtual potential temperature, returning potential temperature in
 * the absence of moisture.
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Thetav (int i, int j, int k,
        const amrex::Array4<const amrex::Real>& cell_data,
        const int moist_flag)
{
    amrex::Real thetav = cell_data(i,j,k,RhoTheta_comp) / cell_data(i,j,k,Rho_comp);
    if (moist_flag > RhoQ1_comp) {
        // RhoQr_comp == moist_flag
        thetav *= (1.0 + 0.61 * cell_data(i,j,k,RhoQ1_comp) / cell_data(i,j,k,Rho_comp)
                              - cell_data(i,j,k,moist_flag) / cell_data(i,j,k,Rho_comp));
    } else if (moist_flag > 0) {
        // We only have Qv
        thetav *= (1.0 + 0.61 * cell_data(i,j,k,RhoQ1_comp) / cell_data(i,j,k,Rho_comp));
    }
    return thetav;
}
#endif

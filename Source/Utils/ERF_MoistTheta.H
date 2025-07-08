#ifndef ERF_UTILS_MOISTTHETA_H_
#define ERF_UTILS_MOISTTHETA_H_

#include "ERF_EOS.H"
#include "ERF_Constants.H"

/*
 * Calculate virtual potential temperature, returning potential temperature in
 * the absence of moisture.
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Thetav (int i, int j, int k,
        const amrex::Array4<const amrex::Real>& cell_data,
        const int RhoQv_comp,
        const int RhoQc_comp,
        const int RhoQr_comp)
{
    amrex::Real thetav = cell_data(i,j,k,RhoTheta_comp) / cell_data(i,j,k,Rho_comp);

    if (RhoQr_comp > 0) {
        thetav *= (1.0 + 0.61 * cell_data(i,j,k,RhoQv_comp) / cell_data(i,j,k,Rho_comp)
                              - cell_data(i,j,k,RhoQc_comp) / cell_data(i,j,k,Rho_comp)
                              - cell_data(i,j,k,RhoQr_comp) / cell_data(i,j,k,Rho_comp));
    } else if (RhoQv_comp > 0) {
        thetav *= (1.0 + 0.61 * cell_data(i,j,k,RhoQv_comp) / cell_data(i,j,k,Rho_comp)
                              - cell_data(i,j,k,RhoQc_comp) / cell_data(i,j,k,Rho_comp));
    }

    return thetav;
}

/*
 * Calculate linearized liquid-water potential temperature, returning potential
 * temperature in the absence of moisture.
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Thetal (int i, int j, int k,
        const amrex::Array4<const amrex::Real>& cell_data,
        const int RhoQv_comp,
        const int RhoQc_comp,
        const int RhoQr_comp)
{
    amrex::Real thetal = cell_data(i,j,k,RhoTheta_comp) / cell_data(i,j,k,Rho_comp);

    if (RhoQc_comp > 0) {
        amrex::Real         ql  = cell_data(i,j,k,RhoQc_comp) / cell_data(i,j,k,Rho_comp);
        if (RhoQr_comp > 0) ql += cell_data(i,j,k,RhoQr_comp) / cell_data(i,j,k,Rho_comp);

        amrex::Real T = getTgivenRandRTh(cell_data(i,j,k,Rho_comp),
                                         cell_data(i,j,k,RhoTheta_comp),
                                         cell_data(i,j,k,RhoQv_comp) / cell_data(i,j,k,Rho_comp));
        amrex::Real theta_by_T = thetal / T;

        thetal -= theta_by_T * (L_v / Cp_d) * ql;
    }

    return thetal;
}
#endif

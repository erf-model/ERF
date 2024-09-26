#ifndef ERF_METGRIDUTIL_H_
#define ERF_METGRIDUTIL_H_

#include <ERF.H>
#include <ERF_EOS.H>
#include <ERF_Constants.H>
#include <ERF_Utils.H>
#include <ERF_prob_common.H>
#include <ERF_HSE_utils.H>

void
read_from_metgrid (int lev,
                   const amrex::Box& domain,
                   const std::string& fname,
                   std::string& NC_dateTime,
                   amrex::Real& NC_epochTime,
                   int& flag_psfc,
                   int& flag_msfu,
                   int& flag_msfv,
                   int& flag_msfm,
                   int& flag_hgt,
                   int& flag_sst,
                   int& flag_lmask,
                   int& NC_nx,
                   int& NC_ny,
                   amrex::Real& NC_dx,
                   amrex::Real& NC_dy,
                   amrex::FArrayBox& NC_xvel_fab,
                   amrex::FArrayBox& NC_yvel_fab,
                   amrex::FArrayBox& NC_temp_fab,
                   amrex::FArrayBox& NC_rhum_fab,
                   amrex::FArrayBox& NC_pres_fab,
                   amrex::FArrayBox& NC_ght_fab,
                   amrex::FArrayBox& NC_hgt_fab,
                   amrex::FArrayBox& NC_psfc_fab,
                   amrex::FArrayBox& NC_msfu_fab,
                   amrex::FArrayBox& NC_msfv_fab,
                   amrex::FArrayBox& NC_msfm_fab,
                   amrex::FArrayBox& NC_sst_fab,
                   amrex::FArrayBox& NC_LAT_fab,
                   amrex::FArrayBox& NC_LON_fab,
                   amrex::IArrayBox& NC_lmask_iab,
                   amrex::Real& Latitude,
                   amrex::Real& Longitude,
                   amrex::Geometry& geom);



void
init_terrain_from_metgrid (amrex::FArrayBox& z_phys_nd_fab,
                           const amrex::Vector<amrex::FArrayBox>& NC_hgt_fab);

void
init_state_from_metgrid (const bool use_moisture,
                         const amrex::Real l_rdOcp,
                         amrex::Box& tbxc,
                         amrex::Box& tbxu,
                         amrex::Box& tbxv,
                         amrex::FArrayBox& state_fab,
                         amrex::FArrayBox& x_vel_fab,
                         amrex::FArrayBox& y_vel_fab,
                         amrex::FArrayBox& z_vel_fab,
                         amrex::FArrayBox& z_phys_nd_fab,
                         const amrex::Vector<amrex::FArrayBox>& NC_hgt_fab,
                         const amrex::Vector<amrex::FArrayBox>& NC_ght_fab,
                         const amrex::Vector<amrex::FArrayBox>& NC_xvel_fab,
                         const amrex::Vector<amrex::FArrayBox>& NC_yvel_fab,
                         const amrex::Vector<amrex::FArrayBox>& NC_zvel_fab,
                         const amrex::Vector<amrex::FArrayBox>& NC_temp_fab,
                         const amrex::Vector<amrex::FArrayBox>& NC_rhum_fab,
                         amrex::Vector<amrex::FArrayBox>& theta_fab,
                         amrex::Vector<amrex::FArrayBox>& mxrat_fab,
                         amrex::Vector<amrex::Vector<amrex::FArrayBox>>& fabs_for_bcs,
                         const amrex::Array4<const int>& mask_c_arr,
                         const amrex::Array4<const int>& mask_u_arr,
                         const amrex::Array4<const int>& mask_v_arr);

void
init_msfs_from_metgrid (amrex::FArrayBox& msfu_fab,
                        amrex::FArrayBox& msfv_fab,
                        amrex::FArrayBox& msfm_fab,
                        const int& flag_msfu,
                        const int& flag_msfv,
                        const int& flag_msfm,
                        const amrex::Vector<amrex::FArrayBox>& NC_MSFU_fab,
                        const amrex::Vector<amrex::FArrayBox>& NC_MSFV_fab,
                        const amrex::Vector<amrex::FArrayBox>& NC_MSFM_fab);

void
init_base_state_from_metgrid (const bool use_moisture,
                              const amrex::Real l_rdOcp,
                              const amrex::Box& valid_bx,
                              const amrex::Vector<int>& flag_psfc,
                              amrex::FArrayBox& state,
                              amrex::FArrayBox& r_hse_fab,
                              amrex::FArrayBox& p_hse_fab,
                              amrex::FArrayBox& pi_hse_fab,
                              amrex::FArrayBox& z_phys_cc_fab,
                              const amrex::Vector<amrex::FArrayBox>& NC_ght_fab,
                              const amrex::Vector<amrex::FArrayBox>& NC_psfc_fab,
                              amrex::Vector<amrex::Vector<amrex::FArrayBox>>& fabs_for_bcs,
                              const amrex::Array4<const int>& mask_c_arr);

AMREX_FORCE_INLINE
AMREX_GPU_DEVICE
void
calc_rho_p (const int& kmax,
            const int& flag_psfc,
            const amrex::Real& psfc,
            amrex::Real* Thetad_vec,
            amrex::Real* Thetam_vec,
            amrex::Real* Q_vec,
            amrex::Real* z_vec,
            amrex::Real* Rhod_vec,
            amrex::Real* Rhom_vec,
            amrex::Real* Pd_vec,
            amrex::Real* Pm_vec)
{
    const int maxiter = 10;
    const amrex::Real tol = 1.0e-12;

    // Calculate or use moist pressure at the surface.
    amrex::Real Psurf;
    if (flag_psfc == 1) {
        Psurf = psfc;
    } else {
        amrex::Real t_0 = 290.0; // WRF's model_config_rec%base_temp
        amrex::Real a   = 50.0; // WRF's model_config_rec%base_lapse
        Psurf = p_0*exp(-t_0/a+std::pow((std::pow(t_0/a, 2)-2.0*CONST_GRAV*z_vec[0]/(a*R_d)), 0.5));
    }

    // Iterations for the first CC point that is 1/2 dz off the surface
    {
        amrex::Real half_dz = z_vec[0];
        amrex::Real qvf = 1.0+(R_v/R_d)*Q_vec[0];
        Thetam_vec[0] = Thetad_vec[0]*qvf;
        for (int it=0; it<maxiter; it++) {
            Pm_vec[0] = Psurf-half_dz*(Rhom_vec[0])*(1.0+Q_vec[0])*CONST_GRAV;
            if (Pm_vec[0] < 0.0) Pm_vec[0] = 0.0;
            Rhom_vec[0] = (p_0/(R_d*Thetam_vec[0]))*std::pow(Pm_vec[0]/p_0, iGamma);
        } // it
    }

    // Integrate from the first CC point to the top boundary.
    for (int k=1; k<=kmax; k++) {

        // Assignment only for downward integration in following loop
        amrex::Real qvf = 1.0+(R_v/R_d)*Q_vec[k];
        Thetam_vec[k]   = Thetad_vec[k]*qvf;

        // Vertical grid spacing
        amrex::Real dz = z_vec[k]-z_vec[k-1];

        // Establish known constant
        amrex::Real rho_tot_lo = Rhom_vec[k-1] * (1. + Q_vec[k-1]);
        amrex::Real C = -Pm_vec[k-1] + 0.5*rho_tot_lo*CONST_GRAV*dz;

        // Initial guess and residual
        Pm_vec[k]   = Pm_vec[k-1];
        Rhom_vec[k] = getRhogivenThetaPress(Thetad_vec[k],
                                            Pm_vec[k],
                                            R_d/Cp_d,
                                            Q_vec[k]);
        amrex::Real rho_tot_hi = Rhom_vec[k] * (1. + Q_vec[k]);
        amrex::Real F = Pm_vec[k] + 0.5*rho_tot_hi*CONST_GRAV*dz + C;

        // Do iterations
        if (std::abs(F)>tol) HSEutils::Newton_Raphson_hse(tol, R_d/Cp_d, dz,
                                                          CONST_GRAV, C, Thetad_vec[k],
                                                          Q_vec[k], Q_vec[k],
                                                          Pm_vec[k], Rhom_vec[k], F);
    } // k

    // integrate from the top back down to get dry pressure and density.
      Pd_vec[kmax] = Pm_vec[kmax];
    Rhod_vec[kmax] = (p_0/(R_d*Thetam_vec[kmax]))*std::pow(Pd_vec[kmax]/p_0, iGamma);
    for (int k=kmax-1; k>=0; k--) {
        amrex::Real dz = z_vec[k+1]-z_vec[k];
        Rhod_vec[k] = Rhod_vec[k+1]; // an initial guess.
        for (int it=0; it<maxiter; it++) {
            Pd_vec[k] = Pd_vec[k+1]+0.5*dz*(Rhod_vec[k]+Rhod_vec[k+1])*CONST_GRAV;
            if (Pd_vec[k] < 0.0) Pd_vec[k] = 0.0;
            Rhod_vec[k] = (p_0/(R_d*Thetam_vec[k]))*std::pow(Pd_vec[k]/p_0, iGamma);
        } // it
    } // k
}

AMREX_FORCE_INLINE
AMREX_GPU_DEVICE
amrex::Real
interpolate_column_metgrid (const int& i,
                            const int& j,
                            const int& k,
                            char stag,
                            int src_comp,
                            const amrex::Array4<amrex::Real const>& orig_z,
                            const amrex::Array4<amrex::Real const>& orig_data,
                            const amrex::Array4<amrex::Real const>&  new_z)
{
    // This subroutine is a bit ham-handed and can be cleaned up later.
    int imax_orig = amrex::ubound(amrex::Box(orig_data)).x;
    int jmax_orig = amrex::ubound(amrex::Box(orig_data)).y;
    int kmax_orig = amrex::ubound(amrex::Box(orig_data)).z;

    amrex::Real z;
    if (stag == 'X') {
        z = 0.25*(new_z(i,j,k)+new_z(i,j+1,k)+new_z(i,j,k+1)+new_z(i,j+1,k+1));
    }
    else if (stag == 'Y') {
        z = 0.25*(new_z(i,j,k)+new_z(i+1,j,k)+new_z(i,j,k+1)+new_z(i+1,j,k+1));
    }
    else if (stag == 'M') {
        z = 0.125*(new_z(i,j,k  )+new_z(i,j+1,k  )+new_z(i+1,j,k  )+new_z(i+1,j+1,k  )+
                   new_z(i,j,k+1)+new_z(i,j+1,k+1)+new_z(i+1,j,k+1)+new_z(i+1,j+1,k+1));
    }

    amrex::Real z0, z1;
    int klow   = -1;
    int khi0   = -1;
    amrex::Real dzlow =  1.0e12;
    amrex::Real dzhi0 = -1.0e12;
    for (int kk = 0; kk < kmax_orig; kk++) {
        amrex::Real orig_z_stag = 0.0;
        if (stag == 'M') {
            orig_z_stag = orig_z(i,j,kk);
        }
        if (stag == 'X') {
            if (i == 0) {
                orig_z_stag = orig_z(i,j,kk);
            }
            else if (i == imax_orig) {
                orig_z_stag = orig_z(imax_orig-1,j,kk);
            }
            else {
                orig_z_stag = 0.5*(orig_z(i,j,kk)+orig_z(i-1,j,kk));
            }
        }
        else if (stag == 'Y') {
            if (j == 0) {
                orig_z_stag = orig_z(i,j,kk);
            }
            else if (j == jmax_orig) {
                orig_z_stag = orig_z(i,jmax_orig-1,kk);
            }
            else {
                orig_z_stag = 0.5*(orig_z(i,j,kk)+orig_z(i,j-1,kk));
            }
        }

        amrex::Real dz = z - orig_z_stag;
        if ((dz < 0.0) && (dz > dzhi0)) {
            dzhi0 = dz;
            khi0  = kk;
            z1    = orig_z_stag;
        }
        if ((dz >= 0.0) && (dz < dzlow)) {
            dzlow = dz;
            klow  = kk;
            z0    = orig_z_stag;
        }
    } // kk

    // extrapolate below the bottom surface
    if (klow == -1) {
        int khi1   = -1;
        amrex::Real dzhi1 = -1.0e12;
        for (int kk = 0; kk < kmax_orig; kk++) {
            amrex::Real orig_z_stag = 0.0;
            if (stag == 'M') {
                orig_z_stag = orig_z(i,j,kk);
            }
            else if (stag == 'X') {
                if (i == 0) {
                    orig_z_stag = orig_z(i,j,kk);
                }
                else if (i == imax_orig) {
                    orig_z_stag = orig_z(imax_orig-1,j,kk);
                }
                else {
                    orig_z_stag = 0.5*(orig_z(i,j,kk)+orig_z(i-1,j,kk));
                }
            }
            else if (stag == 'Y') {
                if (j == 0) {
                    orig_z_stag = orig_z(i,j,kk);
                }
                else if (j == jmax_orig) {
                    orig_z_stag = orig_z(i,jmax_orig-1,kk);
                }
                else {
                    orig_z_stag = 0.5*(orig_z(i,j,kk)+orig_z(i,j-1,kk));
                }
            }
            amrex::Real dz = z - orig_z_stag;
            if ((dz < 0.0) && (dz > dzhi1) && (kk != khi0)) {
                dzhi1 = dz;
                khi1  = kk;
                z1    = orig_z_stag;
            }
        }
        amrex::Real y0 = orig_data(i,j,khi0,src_comp);
        amrex::Real y1 = orig_data(i,j,khi1,src_comp);
        return ( y0-(y1-y0)/(z1-z0)*(z0-z) );

    // Extrapolate above the top surface
    } else if (khi0 == -1) {
        khi0 = klow - 1;
        int khi1 = klow;
        if (stag == 'M') {
            z0 = orig_z(i,j,khi0);
        }
        else if (stag == 'X') {
            if (i == 0) {
                z0 = orig_z(i,j,khi0);
            }
            else if (i == imax_orig) {
                z0 = orig_z(imax_orig-1,j,khi0);
            }
            else {
                z0 = 0.5*(orig_z(i,j,khi0)+orig_z(i-1,j,khi0));
            }
        }
        else if (stag == 'Y') {
            if (j == 0) {
                z0 = orig_z(i,j,khi0);
            }
            else if (j == jmax_orig) {
                z0 = orig_z(i,jmax_orig-1,khi0);
            }
            else {
                z0 = 0.5*(orig_z(i,j,khi0)+orig_z(i,j-1,khi0));
            }
        }
        amrex::Real y0 = orig_data(i,j,khi0,src_comp);
        amrex::Real y1 = orig_data(i,j,khi1,src_comp);
        return ( y0+(y1-y0)/(z1-z0)*(z-z0) );
    } else {
        // interpolate
        amrex::Real y0 = orig_data(i,j,klow,src_comp);
        amrex::Real y1 = orig_data(i,j,khi0,src_comp);
        return ( y0+(y1-y0)/(z1-z0)*(z-z0) );

    }
}

AMREX_FORCE_INLINE
AMREX_GPU_DEVICE
void
rh_to_mxrat (int i,
             int j,
             int k,
             const amrex::Array4<amrex::Real const>& rhum,
             const amrex::Array4<amrex::Real const>& temp,
             const amrex::Array4<amrex::Real const>& pres,
             const amrex::Array4<amrex::Real>& mxrat)
{
    amrex::Real qv_max_p_safe = 10000.0; // WRF default value
    amrex::Real qv_max_flag   = 1.0e-5; // WRF default value
    amrex::Real qv_max_value  = 3.0e-6; // WRF default value
    amrex::Real qv_min_p_safe = 110000.0; // WRF default value
    amrex::Real qv_min_flag   = 1.0e-6; // WRF default value
    amrex::Real qv_min_value  = 1.0e-6; // WRF default value
    amrex::Real eps   = 0.622;
    amrex::Real svp1  = 0.6112;
    amrex::Real svp2  = 17.67;
    amrex::Real svp3  = 29.65;
    amrex::Real svpt0 = 273.15;
    // WRF's method when model_config_rec%rh2qv_wrt_liquid=.true. (default behavior)
    if (temp(i,j,k) != 0.0) {
        amrex::Real es=0.01*rhum(i,j,k)*svp1*10.0*exp(svp2*(temp(i,j,k)-svpt0)/(temp(i,j,k)-svp3));
        if (es >= pres(i,j,k)/100.0) {
            // vapor pressure exceeds total pressure
            mxrat(i,j,k) = std::pow(10.0, -6);
        }
        else {
            mxrat(i,j,k) = amrex::max(eps*es/(pres(i,j,k)/100.0-es), 1.0e-6);
        }
    }
    else {
        // I don't know why there's a fringe case handled in WRF where T is absolute zero...
        // Let's just deal with it here in case we also end up needing it.
        mxrat(i,j,k) = 1.0e-6;
    }
    // See the below comment from WRF dyn_em/module_initialize_real.F rh_to_mxrat1.
    // For pressures above a defined level, reasonable Qv values should be
    // a certain value or smaller. If they are larger than this, the input data
    // probably had "missing" RH, and we filled in some values. This is an
    // attempt to catch those. Also, set the minimum value for the entire
    // domain that is above the selected pressure level.
    if (pres(i,j,k) < qv_max_p_safe) {
        if (mxrat(i,j,k) > qv_max_flag) {
            mxrat(i,j,k) = qv_max_value;
        }
    }
    if (pres(i,j,k) < qv_min_p_safe) {
        if (mxrat(i,j,k) < qv_min_flag) {
            mxrat(i,j,k) = qv_min_value;
        }
    }
}
#endif

#ifndef _PBLMODELS_H_
#define _PBLMODELS_H_

#include <AMReX_GpuContainers.H>

#include <IndexDefines.H>
#include <ERF_Constants.H>

/**
 * Function for computing the QKE equilibrium values.
 *
 * @param[in] bx_QKE grown box to compute QKE equilibrium
 * @param[in] u velocity in x-dir
 * @param[in] v velocity in y-dir
 * @param[in] cell_data conserved cell center vars
 * @param[in] cell_prim primitive cell center vars
 * @param[in] theta_mean average theta
 * @param[in] K_turb turbulent viscosity
 * @param[in,out] QKE_equil_arr Equilibrium QKE
 * @param[in] cellSizeInv inverse cell size array
 * @param[in] domain box of the whole domain
 * @param[in] pbl_B1_l a parameter
 */
void
ComputeEquilibriumQKE (amrex::Box bx_QKE,
                       const amrex::Array4<const amrex::Real>& uvel,
                       const amrex::Array4<const amrex::Real>& vvel,
                       const amrex::Array4<const amrex::Real>& cell_data,
                       const amrex::Array4<const amrex::Real>& cell_prim,
                       const amrex::Array4<const amrex::Real>& tm_arr,
                       const amrex::Array4<const amrex::Real>& K_turb,
                       amrex::Array4<amrex::Real>& QKE_equil_arr,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Box& domain,
                       amrex::Real pbl_B1_l);

/**
 * Function for computing the QKE source terms.
 *
 * @param[in] u velocity in x-dir
 * @param[in] v velocity in y-dir
 * @param[in] cell_data conserved cell center vars
 * @param[in] cell_prim primitive cell center vars
 * @param[in] K_turb turbulent viscosity
 * @param[in] cellSizeInv inverse cell size array
 * @param[in] domain box of the whole domain
 * @param[in] pbl_B1_l a parameter
 * @param[in] theta_mean average theta
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeQKESourceTerms (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& uvel,
                       const amrex::Array4<const amrex::Real>& vvel,
                       const amrex::Array4<const amrex::Real>& cell_data,
                       const amrex::Array4<const amrex::Real>& cell_prim,
                       const amrex::Array4<const amrex::Real>& K_turb,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Box& domain,
                       amrex::Real pbl_B1_l,
                       const amrex::Real theta_mean)
{
  // Compute some relevant derivatives
  amrex::Real dz_inv = cellSizeInv[2];
  int izmin = domain.smallEnd(2);
  int izmax = domain.bigEnd(2);
  amrex::Real source_term = 0.0;

  amrex::Real dthetadz = 0.5*(cell_prim(i,j,k+1,PrimTheta_comp) - cell_prim(i,j,k-1,PrimTheta_comp))*dz_inv;
  amrex::Real dudz = 0.25*(uvel(i,j,k+1) - uvel(i,j,k-1) + uvel(i+1,j,k+1) - uvel(i+1,j,k-1))*dz_inv;
  amrex::Real dvdz = 0.25*(vvel(i,j,k+1) - vvel(i,j,k-1) + vvel(i,j+1,k+1) - vvel(i,j+1,k-1))*dz_inv;

  // Production (We store mu_turb, which is 0.5*K_turb)
  source_term += 4.0*K_turb(i,j,k,EddyDiff::Mom_v) * (dudz*dudz + dvdz*dvdz);

  // Bouyancy
  source_term -= 2.0*(CONST_GRAV/theta_mean)*K_turb(i,j,k,EddyDiff::Theta_v)*dthetadz;

  // Dissipation
  amrex::Real qke = cell_prim(i,j,k,PrimQKE_comp);
  if (std::abs(qke) > 0.0) {
    source_term -= 2.0 * cell_data(i,j,k,Rho_comp) * std::pow(qke,1.5) /
      (pbl_B1_l * K_turb(i,j,k,EddyDiff::PBL_lengthscale));
  }

  return source_term;
}
#endif

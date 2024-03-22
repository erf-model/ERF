#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_YAFluxRegister.H>
#include "DataStruct.H"
#include "IndexDefines.H"
#include "ABLMost.H"

/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 */
void erf_slow_rhs_pre (int level, int finest_level, int nrk,
                      amrex::Real dt,
                      amrex::Vector<amrex::MultiFab>& S_rhs,
                      amrex::Vector<amrex::MultiFab>& S_data,
                      const amrex::MultiFab& S_prim,
                            amrex::Vector<amrex::MultiFab >& S_scratch,
                      const amrex::MultiFab& xvel,
                      const amrex::MultiFab& yvel,
                      const amrex::MultiFab& zvel,
                      std::unique_ptr<amrex::MultiFab>& z_t,
                            amrex::MultiFab& Omega,
                      const amrex::MultiFab& source,
                      const amrex::MultiFab& buoyancy,
#if defined(ERF_USE_RRTMGP)
                      const amrex::MultiFab& qheating_rates,
#endif
                            amrex::MultiFab* Tau11,
                            amrex::MultiFab* Tau22,
                            amrex::MultiFab* Tau33,
                            amrex::MultiFab* Tau12,
                            amrex::MultiFab* Tau13,
                            amrex::MultiFab* Tau21,
                            amrex::MultiFab* Tau23,
                            amrex::MultiFab* Tau31,
                            amrex::MultiFab* Tau32,
                            amrex::MultiFab* SmnSmn,
                            amrex::MultiFab* eddyDiffs,
                            amrex::MultiFab* Hfx3,
                            amrex::MultiFab* Diss,
                      const amrex::Geometry geom,
                      const SolverChoice& solverChoice,
                      std::unique_ptr<ABLMost>& most,
                      const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                      const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                      std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                      std::unique_ptr<amrex::MultiFab>& dJ,
                      const amrex::MultiFab* p0,
                      std::unique_ptr<amrex::MultiFab>& mapfac_m,
                      std::unique_ptr<amrex::MultiFab>& mapfac_u,
                      std::unique_ptr<amrex::MultiFab>& mapfac_v,
                      amrex::YAFluxRegister* fr_as_crse,
                      amrex::YAFluxRegister* fr_as_fine,
                      const amrex::Real* dptr_rhotheta_src,
                      const amrex::Real* dptr_u_geos,
                      const amrex::Real* dptr_v_geos,
                      const amrex::Real* dptr_wbar_sub,
                      const amrex::Vector<amrex::Real*> d_rayleigh_ptrs_at_lev);

/**
 * Function for computing the slow RHS for the evolution equations for the scalars other than density or potential temperature
 *
 */
void erf_slow_rhs_post (int level, int finest_level, int nrk,
                       amrex::Real dt,
                       amrex::Vector<amrex::MultiFab>& S_rhs,
                       amrex::Vector<amrex::MultiFab>& S_old,
                       amrex::Vector<amrex::MultiFab>& S_new,
                       amrex::Vector<amrex::MultiFab>& S_data,
                       const amrex::MultiFab& S_prim,
                       amrex::Vector<amrex::MultiFab >& S_scratch,
                       const amrex::MultiFab& xvel,
                       const amrex::MultiFab& yvel,
                       const amrex::MultiFab& zvel,
                       const amrex::MultiFab& source,
                       const amrex::MultiFab* SmnSmn,
                       const amrex::MultiFab* eddyDiffs,
                             amrex::MultiFab* Hfx3,
                             amrex::MultiFab* Diss,
                       const amrex::Geometry geom,
                       const SolverChoice& solverChoice,
                       std::unique_ptr<ABLMost>& most,
                       const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                       const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                       std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                       std::unique_ptr<amrex::MultiFab>& dJ_old,
                       std::unique_ptr<amrex::MultiFab>& dJ_new,
                       std::unique_ptr<amrex::MultiFab>& mapfac_m,
                       std::unique_ptr<amrex::MultiFab>& mapfac_u,
                       std::unique_ptr<amrex::MultiFab>& mapfac_v,
#if defined(ERF_USE_NETCDF)
                       const bool& moist_zero,
                       const amrex::Real& bdy_time_interval,
                       const amrex::Real& start_bdy_time,
                       const amrex::Real& new_stage_time,
                       int  width,
                       int  set_width,
                       amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
                       amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
                       amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
                       amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi,
#endif
                       const amrex::Real* dptr_rhoqt_src,
                       const amrex::Real* dptr_wbar_sub,
                       amrex::YAFluxRegister* fr_as_crse,
                       amrex::YAFluxRegister* fr_as_fine);

/**
 * Function for computing the fast RHS with no terrain
 *
 */
void erf_fast_rhs_N (int step, int nrk, int level, int finest_level,
                     amrex::Vector<amrex::MultiFab >& S_slow_rhs,
                     const amrex::Vector<amrex::MultiFab >& S_prev,
                     amrex::Vector<amrex::MultiFab >& S_stage_data,
                     const amrex::MultiFab& S_stage_prim,
                     const amrex::MultiFab& pi_stage,
                     const amrex::MultiFab& fast_coeffs,
                     amrex::Vector<amrex::MultiFab >& S_data,
                     amrex::Vector<amrex::MultiFab >& S_scratch,
                     const amrex::Geometry geom,
                     const amrex::Real gravity,
                     const amrex::Real dtau, const amrex::Real beta_s,
                     const amrex::Real facinv,
                     std::unique_ptr<amrex::MultiFab>& mapfac_m,
                     std::unique_ptr<amrex::MultiFab>& mapfac_u,
                     std::unique_ptr<amrex::MultiFab>& mapfac_v,
                     amrex::YAFluxRegister* fr_as_crse,
                     amrex::YAFluxRegister* fr_as_fine,
                     bool l_use_moisture, bool l_reflux);

/**
 * Function for computing the fast RHS with fixed terrain
 *
 */
void erf_fast_rhs_T (int step, int nrk, int level, int finest_level,
                     amrex::Vector<amrex::MultiFab >& S_slow_rhs,
                     const amrex::Vector<amrex::MultiFab >& S_prev,
                     amrex::Vector<amrex::MultiFab >& S_stage_data,
                     const amrex::MultiFab& S_stage_prim,
                     const amrex::MultiFab& pi_stage,
                     const amrex::MultiFab& fast_coeffs,
                     amrex::Vector<amrex::MultiFab >& S_data,
                     amrex::Vector<amrex::MultiFab >& S_scratch,
                     const amrex::Geometry geom,
                     const amrex::Real gravity,
                           amrex::MultiFab& Omega,
                     std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                     std::unique_ptr<amrex::MultiFab>& detJ_cc,
                     const amrex::Real dtau, const amrex::Real beta_s,
                     const amrex::Real facinv,
                     std::unique_ptr<amrex::MultiFab>& mapfac_m,
                     std::unique_ptr<amrex::MultiFab>& mapfac_u,
                     std::unique_ptr<amrex::MultiFab>& mapfac_v,
                     amrex::YAFluxRegister* fr_as_crse,
                     amrex::YAFluxRegister* fr_as_fine,
                     bool l_use_moisture, bool l_reflux);

/**
 * Function for computing the fast RHS with moving terrain
 *
 */
void erf_fast_rhs_MT (int step, int nrk, int level, int finest_level,
                      amrex::Vector<amrex::MultiFab >& S_slow_rhs,
                      const amrex::Vector<amrex::MultiFab >& S_prev,
                      amrex::Vector<amrex::MultiFab >& S_stg_data,
                      const amrex::MultiFab& S_stg_prim,
                      const amrex::MultiFab& pi_stage,
                      const amrex::MultiFab& fast_coeffs,
                      amrex::Vector<amrex::MultiFab >& S_data,
                      amrex::Vector<amrex::MultiFab >& S_scratch,
                      const amrex::Geometry geom,
                      const amrex::Real gravity,
                      const bool use_lagged_delta_rt,
                            amrex::MultiFab& Omega,
                      std::unique_ptr<amrex::MultiFab>& z_t_rk,
                      const amrex::MultiFab* z_t_pert,
                      std::unique_ptr<amrex::MultiFab>& z_phys_nd_old,
                      std::unique_ptr<amrex::MultiFab>& z_phys_nd_new,
                      std::unique_ptr<amrex::MultiFab>& z_phys_nd_stg,
                      std::unique_ptr<amrex::MultiFab>& detJ_cc_old,
                      std::unique_ptr<amrex::MultiFab>& detJ_cc_new,
                      std::unique_ptr<amrex::MultiFab>& detJ_cc_stg,
                      const amrex::Real dtau, const amrex::Real beta_s,
                      const amrex::Real facinv,
                      std::unique_ptr<amrex::MultiFab>& mapfac_m,
                      std::unique_ptr<amrex::MultiFab>& mapfac_u,
                      std::unique_ptr<amrex::MultiFab>& mapfac_v,
                      amrex::YAFluxRegister* fr_as_crse,
                      amrex::YAFluxRegister* fr_as_fine,
                      bool l_use_moisture, bool l_reflux);

/**
 * Function for computing the coefficients for the tridiagonal solver used in the fast
 * integrator (the acoustic substepping).
 */
void make_fast_coeffs (int level,
                       amrex::MultiFab& fast_coeffs,
                       amrex::Vector<amrex::MultiFab >& S_stage_data,
                       const amrex::MultiFab& S_stage_prim,
                       const amrex::MultiFab& pi_stage,
                       const amrex::Geometry geom,
                       const bool use_moisture,
                       const bool use_terrain,
                       const amrex::Real gravity,
                       const amrex::Real c_p,
                       std::unique_ptr<amrex::MultiFab>& detJ_cc,
                       const amrex::MultiFab* r0,
                       const amrex::MultiFab* pi0,
                       const amrex::Real dtau,
                       const amrex::Real beta_s);

/**
 * Function for computing the buoyancy term to be used in the evolution
 * equation for the z-component of momentum in the slow integrator.  There
 * are three options for how buoyancy is computed (two are the same in the absence of moisture).
 */
void make_buoyancy (amrex::Vector<  amrex::MultiFab>& S_data,
                    const           amrex::MultiFab & S_prim,
                    amrex::MultiFab& buoyancy,
                    const amrex::Geometry geom,
                    const SolverChoice& solverChoice,
                    const amrex::MultiFab* r0,
                    const int& qstate_size);
#endif

#ifdef ERF_USE_POISSON_SOLVE
/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 */
void erf_slow_rhs_inc (int level, int nrk,
                       amrex::Real dt,
                       amrex::Vector<amrex::MultiFab>& S_rhs,
                       amrex::Vector<amrex::MultiFab>& S_old,
                       amrex::Vector<amrex::MultiFab>& S_data,
                       const amrex::MultiFab& S_prim,
                       amrex::Vector<amrex::MultiFab >& S_scratch,
                       const amrex::MultiFab& xvel,
                       const amrex::MultiFab& yvel,
                       const amrex::MultiFab& zvel,
                       std::unique_ptr<amrex::MultiFab>& z_t,
                       amrex::MultiFab& Omega,
                       const amrex::MultiFab& source,
                       const amrex::MultiFab& buoyancy,
                       amrex::MultiFab* Tau11,
                       amrex::MultiFab* Tau22,
                       amrex::MultiFab* Tau33,
                       amrex::MultiFab* Tau12,
                       amrex::MultiFab* Tau13,
                       amrex::MultiFab* Tau21,
                       amrex::MultiFab* Tau23,
                       amrex::MultiFab* Tau31,
                       amrex::MultiFab* Tau32,
                       amrex::MultiFab* SmnSmn,
                       amrex::MultiFab* eddyDiffs,
                       amrex::MultiFab* Hfx3,
                       amrex::MultiFab* Diss,
                       const amrex::Geometry geom,
                       const SolverChoice& solverChoice,
                       std::unique_ptr<ABLMost>& most,
                       const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                       const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                       std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                       std::unique_ptr<amrex::MultiFab>& dJ,
                       const amrex::MultiFab* p0,
                       std::unique_ptr<amrex::MultiFab>& mapfac_m,
                       std::unique_ptr<amrex::MultiFab>& mapfac_u,
                       std::unique_ptr<amrex::MultiFab>& mapfac_v,
                       const amrex::Real* dptr_rhotheta_src,
                       const amrex::Real* dptr_u_geos,
                       const amrex::Real* dptr_v_geos,
                       const amrex::Real* dptr_wbar_sub,
                       const Vector<amrex::Real*> d_rayleigh_ptrs_at_lev);
#endif

void ApplySpongeZoneBCs (const SpongeChoice& spongeChoice,
                         const amrex::Geometry geom,
                         const amrex::Box& tbx,
                         const amrex::Box& tby,
                         const amrex::Box& tbz,
                         const amrex::Array4<amrex::Real>& rho_u_rhs,
                         const amrex::Array4<amrex::Real>& rho_v_rhs,
                         const amrex::Array4<amrex::Real>& rho_w_rhs,
                         const amrex::Array4<const amrex::Real>& rho_u,
                         const amrex::Array4<const amrex::Real>& rho_v,
                         const amrex::Array4<const amrex::Real>& rho_w,
                         const amrex::Box& bx,
                         const amrex::Array4<amrex::Real>& cell_rhs,
                         const amrex::Array4<const amrex::Real>& cell_data);




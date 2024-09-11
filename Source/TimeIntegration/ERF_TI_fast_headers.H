#ifndef ERF_FAST_INTEGRATION_H_
#define ERF_FAST_INTEGRATION_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_YAFluxRegister.H>
#include "ERF_DataStruct.H"
#include "ERF_IndexDefines.H"
#include <ERF_TerrainMetrics.H>

#include <ERF_TileNoZ.H>
#include <ERF_prob_common.H>

#ifdef ERF_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

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
                       const amrex::Real beta_s,
                       amrex::GpuArray<ERF_BC, AMREX_SPACEDIM*2> &phys_bc_type);

#endif

#ifndef ERF_FAST_INTEGRATION_H_
#define ERF_FAST_INTEGRATION_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_YAFluxRegister.H>
#include <AMReX_EBMultiFabUtil.H>

#include "ERF_DataStruct.H"
#include "ERF_IndexDefines.H"
#include <ERF_TerrainMetrics.H>
#include <ERF_TileNoZ.H>
#include <ERF_ProbCommon.H>

/**
 * The coarse/fine face fluxes of the last RK stage, in the units the flux registers hold:
 * (rho u)/m_uy, (rho v)/m_vx and the time-weighted (rho w)/(m_x m_y) without terrain;
 * ax (rho u)/m_uy, ay (rho v)/m_vx and the time-weighted Omega/(m_x m_y) with it.
 *
 * A coarse level records in rec[dir] the flux each substep applied through every face, and in
 * rec_w the vertical quantity its z-faces hold at the end of the substep -- (rho w) without
 * terrain, Omega with it -- and in rec_w_old the same at the start of the substep.  A fine level
 * takes the fluxes in tgt[dir] (component 0) on the faces where component 1 is one -- its
 * coarse/fine faces -- in place of its own, and tgt_w as that vertical quantity there at the
 * end of the substep.  The fine level's fluxes through its boundary then add up to the coarse
 * ones.
 */
struct CFFaceFluxes {
    std::array<amrex::MultiFab*, AMREX_SPACEDIM>       rec {{AMREX_D_DECL(nullptr, nullptr, nullptr)}};
    amrex::MultiFab* rec_w     = nullptr;
    amrex::MultiFab* rec_w_old = nullptr;
    std::array<const amrex::MultiFab*, AMREX_SPACEDIM> tgt {{AMREX_D_DECL(nullptr, nullptr, nullptr)}};
    const amrex::MultiFab* tgt_w = nullptr;
};

/**
 * Function for computing the fast RHS with no terrain but constant or stretched vertical resolution
 *
 */
void erf_substep_NS (int step, int nrk, int level, int finest_level,
                     amrex::Vector<amrex::MultiFab >& S_slow_rhs,
                     const amrex::Vector<amrex::MultiFab >& S_prev,
                     amrex::Vector<amrex::MultiFab >& S_stage_data,
                     const amrex::MultiFab& S_stage_prim,
                     const amrex::MultiFab& qt,
                     const amrex::MultiFab& pi_stage,
                     const amrex::MultiFab& fast_coeffs,
                     amrex::Vector<amrex::MultiFab >& S_data,
                           amrex::MultiFab& lagged_delta_rt,
                           amrex::MultiFab& avg_xmom,
                           amrex::MultiFab& avg_ymom,
                        amrex::MultiFab& avg_zmom,
                     const amrex::MultiFab& cc_src,
                     const amrex::MultiFab& xmom_src,
                     const amrex::MultiFab& ymom_src,
                     const amrex::MultiFab& zmom_src,
                     const amrex::Geometry geom,
                     const amrex::Real gravity,
                     amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_d,
                     const double dtau,
                     const amrex::Real facinv,
                     amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                     amrex::YAFluxRegister* fr_as_crse,
                     amrex::YAFluxRegister* fr_as_fine,
                     bool l_use_moisture, bool l_reflux, bool l_real_bc,
                     const amrex::Real* sinesq_at_lev_d,
                     const amrex::Real l_damp_coef,
                     const amrex::Real beta_d,
                     const CFFaceFluxes& cf_fluxes);

/**
 * Function for computing the fast RHS with fixed-in-time terrain
 *
 */
void erf_substep_T (int step, int nrk, int level, int finest_level,
                    amrex::Vector<amrex::MultiFab >& S_slow_rhs,
                    const amrex::Vector<amrex::MultiFab >& S_prev,
                    amrex::Vector<amrex::MultiFab >& S_stage_data,
                    const amrex::MultiFab& S_stage_prim,
                    const amrex::MultiFab& qt,
                    const amrex::MultiFab& pi_stage,
                    const amrex::MultiFab& fast_coeffs,
                    amrex::Vector<amrex::MultiFab >& S_data,
                          amrex::MultiFab& lagged_delta_rt,
                          amrex::MultiFab& avg_xmom,
                          amrex::MultiFab& avg_ymom,
                       amrex::MultiFab& avg_zmom,
                    const amrex::MultiFab& cc_src,
                    const amrex::MultiFab& xmom_src,
                    const amrex::MultiFab& ymom_src,
                    const amrex::MultiFab& zmom_src,
                    const amrex::Geometry geom,
                    const amrex::Real gravity,
                    std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                    std::unique_ptr<amrex::MultiFab>& detJ_cc,
                    const double dtau,
                    const amrex::Real facinv,
                    amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                    amrex::YAFluxRegister* fr_as_crse,
                    amrex::YAFluxRegister* fr_as_fine,
                    bool l_use_moisture, bool l_reflux, bool l_real_bc,
                    const amrex::Real* sinesq_at_lev_d,
                    const amrex::Real l_damp_coef,
                    const amrex::Real beta_d,
                    const CFFaceFluxes& cf_fluxes);

/**
 * Function for computing the fast RHS with moving terrain
 *
 */
void erf_substep_MT (int step, int nrk, int level, int finest_level,
                     amrex::Vector<amrex::MultiFab >& S_slow_rhs,
                     const amrex::Vector<amrex::MultiFab >& S_prev,
                     amrex::Vector<amrex::MultiFab >& S_stg_data,
                     const amrex::MultiFab& S_stg_prim,
                     const amrex::MultiFab& qt,
                     const amrex::MultiFab& pi_stage,
                     const amrex::MultiFab& fast_coeffs,
                     amrex::Vector<amrex::MultiFab >& S_data,
                           amrex::MultiFab& lagged_delta_rt,
                           amrex::MultiFab& avg_xmom,
                           amrex::MultiFab& avg_ymom,
                           amrex::MultiFab& avg_zmom,
                     const amrex::MultiFab& cc_src,
                     const amrex::MultiFab& xmom_src,
                     const amrex::MultiFab& ymom_src,
                     const amrex::MultiFab& zmom_src,
                     const amrex::Geometry geom,
                     const amrex::Real gravity,
                     const bool use_lagged_delta_rt,
                     std::unique_ptr<amrex::MultiFab>& z_t_rk,
                     const amrex::MultiFab* z_t_pert,
                     std::unique_ptr<amrex::MultiFab>& z_phys_nd_old,
                     std::unique_ptr<amrex::MultiFab>& z_phys_nd_new,
                     std::unique_ptr<amrex::MultiFab>& z_phys_nd_stg,
                     std::unique_ptr<amrex::MultiFab>& detJ_cc_old,
                     std::unique_ptr<amrex::MultiFab>& detJ_cc_new,
                     std::unique_ptr<amrex::MultiFab>& detJ_cc_stg,
                     const double dtau,
                     const amrex::Real facinv,
                     amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                     amrex::YAFluxRegister* fr_as_crse,
                     amrex::YAFluxRegister* fr_as_fine,
                     bool l_use_moisture, bool l_reflux, bool l_real_bc,
                     const amrex::Real* sinesq_at_lev_d,
                     const amrex::Real l_damp_coef,
                     const amrex::Real beta_d);

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
                       const MeshType mesh_type,
                       const amrex::Real gravity,
                       const amrex::Real c_p,
                       std::unique_ptr<amrex::MultiFab>& detJ_cc,
                       const double dtau,
                       const amrex::Real beta_s,
                       const amrex::MultiFab* z_phys_nd,
                       const bool use_variable_beta_s,
                       const amrex::Real beta_s_min,
                       const amrex::Real beta_s_max,
                       amrex::GpuArray<ERF_BC, AMREX_SPACEDIM*2> &phys_bc_type);

#endif

#ifndef _SLOW_INTEGRATION_H_
#define _SLOW_INTEGRATION_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_YAFluxRegister.H>
#include "DataStruct.H"
#include "IndexDefines.H"
#include "ABLMost.H"

#include <Advection.H>
#include <Diffusion.H>
#include <EddyViscosity.H>
#include <PlaneAverage.H>
#include <TerrainMetrics.H>
#include <TileNoZ.H>

#ifdef ERF_USE_EB
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBMultiFabUtil.H>
#endif

void erf_make_tau_terms (int level, int nrk,
                         const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                         const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                         std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                         amrex::Vector<amrex::MultiFab>& S_data,
                         const amrex::MultiFab& xvel,
                         const amrex::MultiFab& yvel,
                         const amrex::MultiFab& zvel,
                               amrex::MultiFab& Omega,
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
                         std::unique_ptr<amrex::MultiFab>& dJ,
                         std::unique_ptr<amrex::MultiFab>& mapfac_m,
                         std::unique_ptr<amrex::MultiFab>& mapfac_u,
                         std::unique_ptr<amrex::MultiFab>& mapfac_v);

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
                      const amrex::MultiFab& cc_src,
                      const amrex::MultiFab& xmom_src,
                      const amrex::MultiFab& ymom_src,
                      const amrex::MultiFab& zmom_src,
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
                      std::unique_ptr<amrex::MultiFab>& ax,
                      std::unique_ptr<amrex::MultiFab>& ay,
                      std::unique_ptr<amrex::MultiFab>& az,
                      std::unique_ptr<amrex::MultiFab>& dJ,
                      const amrex::MultiFab* p0,
                      std::unique_ptr<amrex::MultiFab>& mapfac_m,
                      std::unique_ptr<amrex::MultiFab>& mapfac_u,
                      std::unique_ptr<amrex::MultiFab>& mapfac_v,
#ifdef ERF_USE_EB
                      amrex::EBFArrayBoxFactory const& ebfact,
#endif
                      amrex::YAFluxRegister* fr_as_crse,
                      amrex::YAFluxRegister* fr_as_fine);

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
                       std::unique_ptr<amrex::MultiFab>& ax,
                       std::unique_ptr<amrex::MultiFab>& ay,
                       std::unique_ptr<amrex::MultiFab>& az,
                       std::unique_ptr<amrex::MultiFab>& dJ_old,
                       std::unique_ptr<amrex::MultiFab>& dJ_new,
                       std::unique_ptr<amrex::MultiFab>& mapfac_m,
                       std::unique_ptr<amrex::MultiFab>& mapfac_u,
                       std::unique_ptr<amrex::MultiFab>& mapfac_v,
#ifdef ERF_USE_EB
                       amrex::EBFArrayBoxFactory const& ebfact,
#endif
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
                       amrex::YAFluxRegister* fr_as_crse,
                       amrex::YAFluxRegister* fr_as_fine);


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
                       amrex::MultiFab& Omega,
                       const amrex::MultiFab& cc_src,
                       const amrex::MultiFab& xmom_src,
                       const amrex::MultiFab& ymom_src,
                       const amrex::MultiFab& zmom_src,
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
                       std::unique_ptr<amrex::MultiFab>& ax,
                       std::unique_ptr<amrex::MultiFab>& ay,
                       std::unique_ptr<amrex::MultiFab>& az,
                       std::unique_ptr<amrex::MultiFab>& dJ,
                       const amrex::MultiFab* p0,
                       std::unique_ptr<amrex::MultiFab>& mapfac_m,
                       std::unique_ptr<amrex::MultiFab>& mapfac_u,
                       std::unique_ptr<amrex::MultiFab>& mapfac_v,
                       const amrex::MultiFab& pp_inc);
#endif
#endif

#ifndef ERF_SLOW_INTEGRATION_H_
#define ERF_SLOW_INTEGRATION_H_

#include "AMReX_MultiFab.H"
#include "AMReX_BCRec.H"
#include "AMReX_YAFluxRegister.H"
#include "AMReX_MultiCutFab.H"
#include "AMReX_EBMultiFabUtil.H"

#include "ERF_DataStruct.H"
#include "ERF_IndexDefines.H"
#include "ERF_SurfaceLayer.H"
#include "ERF_Advection.H"
#include "ERF_Diffusion.H"
#include "ERF_EddyViscosity.H"
#include "ERF_PlaneAverage.H"
#include "ERF_TerrainMetrics.H"
#include "ERF_TileNoZ.H"
#include "ERF_EB.H"

void erf_make_tau_terms (int level, int nrk,
                         const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                         const amrex::MultiFab& z_phys_nd,
                         amrex::Vector<amrex::MultiFab>& S_data,
                         const amrex::MultiFab& xvel,
                         const amrex::MultiFab& yvel,
                         const amrex::MultiFab& zvel,
                         amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Tau_lev,
                         amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Tau_corr_lev,
                               amrex::MultiFab* SmnSmn,
                               amrex::MultiFab* eddyDiffs,
                         const amrex::Geometry geom,
                         const SolverChoice& solverChoice,
                         std::unique_ptr<SurfaceLayer>& SurfLayer,
                         amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_d,
                         const amrex::MultiFab& detJ,
                         amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                         const amrex::MultiFab& ax,
                         const amrex::MultiFab& ay,
                         const amrex::MultiFab& az,
                         const eb_& ebfact);

void copy_surface_tau_for_implicit (amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Tau_lev,
                                    amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Tau_corr_lev);

/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 */
void erf_slow_rhs_pre (int level, int finest_level, int nrk,
                      amrex::Real dt,
                      amrex::Vector<amrex::MultiFab>& S_rhs,
                      amrex::Vector<amrex::MultiFab>& S_old,
                      amrex::Vector<amrex::MultiFab>& S_data,
                              const amrex::MultiFab & S_prim,
                              const amrex::MultiFab & qt,
                            amrex::MultiFab& avg_xmom,
                            amrex::MultiFab& avg_ymom,
                            amrex::MultiFab& avg_zmom,
                      const amrex::MultiFab& xvel,
                      const amrex::MultiFab& yvel,
                      const amrex::MultiFab& zvel,
                      std::unique_ptr<amrex::MultiFab>& z_t,
                      const amrex::MultiFab& cc_src,
                      const amrex::MultiFab& xmom_src,
                      const amrex::MultiFab& ymom_src,
                      const amrex::MultiFab& zmom_src,
                      const amrex::MultiFab& buoyancy,
                      const amrex::MultiFab* zmom_crse_rhs,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Tau_lev,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Tau_corr_lev,
                            amrex::MultiFab* SmnSmn,
                            amrex::MultiFab* eddyDiffs,
                            amrex::MultiFab* Hfx1,
                            amrex::MultiFab* Hfx2,
                            amrex::MultiFab* Hfx3,
                            amrex::MultiFab* Q1fx1,
                            amrex::MultiFab* Q1fx2,
                            amrex::MultiFab* Q1fx3,
                            amrex::MultiFab* Q2fx3,
                            amrex::MultiFab* Diss,
                      const amrex::Geometry geom,
                      const SolverChoice& solverChoice,
                      std::unique_ptr<SurfaceLayer>& SurfLayer,
                      const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                      const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                      const amrex::MultiFab& z_phys_nd,
                      const amrex::MultiFab& z_phys_cc,
                      const amrex::MultiFab& ax,
                      const amrex::MultiFab& ay,
                      const amrex::MultiFab& az,
                      const amrex::MultiFab& dJ,
                      amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_d,
                      amrex::Vector<amrex::MultiFab>& gradp,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                      const eb_& ebfact,
#ifdef ERF_USE_SHOC
                      std::unique_ptr<SHOCInterface>& shoc_lev,
#endif
                      amrex::YAFluxRegister* fr_as_crse,
                      amrex::YAFluxRegister* fr_as_fine);

/**
 * Function for computing the slow RHS for the evolution equations for the scalars other than density or potential temperature
 *
 */
void erf_slow_rhs_post (int level, int finest_level, int nrk,
                       amrex::Real dt,
                       int n_qstate,
                       amrex::Vector<amrex::MultiFab>& S_rhs,
                       amrex::Vector<amrex::MultiFab>& S_old,
                       amrex::Vector<amrex::MultiFab>& S_new,
                       amrex::Vector<amrex::MultiFab>& S_data,
                       const amrex::MultiFab& S_prim,
                             amrex::MultiFab& avg_xmom,
                             amrex::MultiFab& avg_ymom,
                             amrex::MultiFab& avg_zmom,
                       const amrex::MultiFab& xvel,
                       const amrex::MultiFab& yvel,
                       const amrex::MultiFab& zvel,
                       const amrex::MultiFab& source,
                       const amrex::MultiFab* SmnSmn,
                       const amrex::MultiFab* eddyDiffs,
                             amrex::MultiFab* Hfx1,
                             amrex::MultiFab* Hfx2,
                             amrex::MultiFab* Hfx3,
                             amrex::MultiFab* Q1fx1,
                             amrex::MultiFab* Q1fx2,
                             amrex::MultiFab* Q1fx3,
                             amrex::MultiFab* Q2fx3,
                             amrex::MultiFab* Diss,
                       const amrex::Geometry geom,
                       const SolverChoice& solverChoice,
                       std::unique_ptr<SurfaceLayer>& SurfLayer,
                       const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                       const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                       std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                       std::unique_ptr<amrex::MultiFab>& z_phys_cc,
                       std::unique_ptr<amrex::MultiFab>& ax,
                       std::unique_ptr<amrex::MultiFab>& ay,
                       std::unique_ptr<amrex::MultiFab>& az,
                       std::unique_ptr<amrex::MultiFab>& dJ_old,
                       amrex::MultiFab* dJ_new,
                       amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_d,
                       amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                       amrex::EBFArrayBoxFactory const& ebfact,
#if defined(ERF_USE_NETCDF)
                       const bool& moist_zero,
                       const amrex::Real& bdy_time_interval,
                       const amrex::Real& new_stage_time,
                       const amrex::Real& stop_time,
                       int  width,
                       int  set_width,
                       amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
                       amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
                       amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
                       amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi,
#endif
#ifdef ERF_USE_SHOC
                      std::unique_ptr<SHOCInterface>& shoc_lev,
#endif
                       amrex::YAFluxRegister* fr_as_crse,
                       amrex::YAFluxRegister* fr_as_fine);

#endif

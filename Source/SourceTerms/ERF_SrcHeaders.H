#ifndef ERF_INTEGRATION_H_
#define ERF_INTEGRATION_H_

#include <AMReX_MultiFab.H>
#include <AMReX_EBMultiFabUtil.H>

#include "ERF_DataStruct.H"
#include "ERF_InputSoundingData.H"
#include "ERF_TurbPertStruct.H"
#include "ERF_EB.H"

/**
 * Function for computing the buoyancy term to be used in the evolution
 * equation for the z-component of momentum in the slow integrator.  There
 * are three options for how buoyancy is computed (two are the same in the absence of moisture).
 */
void make_buoyancy (int level,
                    const amrex::Vector<amrex::MultiFab>& S_data,
                    const               amrex::MultiFab & S_prim,
                    const               amrex::MultiFab & qt,
                    amrex::MultiFab& buoyancy,
                    const amrex::Geometry geom,
                    const SolverChoice& solverChoice,
                    const amrex::MultiFab& base_state,
                    const int n_qstate,
                    const eb_& ebfact,
                    const int anelastic);

void make_gradp_pert (int level,
                      const SolverChoice& solverChoice,
                      const amrex::Geometry& geom,
                      amrex::Vector<amrex::MultiFab>& S_data,
                      const amrex::MultiFab& p0,
                      const amrex::MultiFab& z_phys_nd,
                      const amrex::MultiFab& z_phys_cc,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                      const eb_& ebfact,
                      amrex::Vector<amrex::MultiFab>& gradp);

void compute_gradp (const amrex::MultiFab& p,
                    const amrex::Geometry& geom,
                    const amrex::MultiFab& z_phys_nd,
                    const amrex::MultiFab& z_phys_cc,
                    amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                    const eb_& ebfact,
                    amrex::Vector<amrex::MultiFab>& gradp,
                    const SolverChoice& solverChoice);

void compute_gradp_interpz (const amrex::MultiFab& p,
                            const amrex::Geometry& geom,
                            const amrex::MultiFab& z_phys_nd,
                            const amrex::MultiFab& z_phys_cc,
                            amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                            amrex::Vector<amrex::MultiFab>& gradp,
                            const SolverChoice& solverChoice);

void make_sources (int level, int nrk,
                   amrex::Real dt, amrex::Real time,
                   const amrex::Vector<amrex::MultiFab>& S_data,
                   const amrex::MultiFab& S_prim,
                         amrex::MultiFab& cc_source,
                   const amrex::MultiFab& base_state,
                   const amrex::MultiFab* z_phys_cc,
                   const amrex::MultiFab& xvel,
                   const amrex::MultiFab& yvel,
                   const amrex::MultiFab* qheating_rates,
                   amrex::MultiFab* terrain_blank,
                   const amrex::Geometry geom,
                   const SolverChoice& solverChoice,
                   amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                   const amrex::Real* dptr_rhotheta_src,
                   const amrex::Real* dptr_rhoqt_src,
                   const amrex::Real* dptr_wbar_sub,
                   const amrex::Vector<amrex::Real*> d_rayleigh_ptrs_at_lev,
                   const amrex::Real* sinesq_at_lev_d,
                   InputSoundingData& input_sounding_data,
                   TurbulentPerturbation& turbPert,
                   bool is_slow_step);

void make_mom_sources (amrex::Real time,
                       amrex::Real dt,
                       const amrex::Vector<amrex::MultiFab>& S_data,
                       const amrex::MultiFab* z_phys_nd,
                       const amrex::MultiFab* z_phys_cc,
                             amrex::Vector<amrex::Real>& stretched_dz_h,
                       const amrex::MultiFab& xvel,
                       const amrex::MultiFab& yvel,
                       const amrex::MultiFab& zvel,
                             amrex::MultiFab& xmom_source,
                             amrex::MultiFab& ymom_source,
                             amrex::MultiFab& zmom_source,
                       const amrex::MultiFab& base_state,
                             amrex::MultiFab* forest_drag,
                             amrex::MultiFab* terrain_blank,
                             amrex::MultiFab* cosPhi_mf,
                             amrex::MultiFab* sinPhi_mf,
                       const amrex::Geometry geom,
                       const SolverChoice& solverChoice,
                             amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                       const amrex::Real* dptr_rhotheta_src,
                       const amrex::Real* dptr_rhoqt_src,
                       const amrex::Real* dptr_wbar_sub,
                       const amrex::Vector<amrex::Real*> d_rayleigh_ptrs_at_lev,
                       const amrex::Real* sinesq_at_lev_d,
                       const amrex::Real* sinesq_stag_at_lev_d,
                       const amrex::Vector<amrex::Real*> d_sponge_ptrs_at_lev,
                       const amrex::Vector<amrex::MultiFab>* forecast_state_at_lev,
                             InputSoundingData& input_sounding_data,
                             bool is_slow_step);

void add_thin_body_sources (amrex::MultiFab& xmom_source,
                            amrex::MultiFab& ymom_source,
                            amrex::MultiFab& zmom_source,
                            std::unique_ptr<amrex::iMultiFab>& xflux_imask_lev,
                            std::unique_ptr<amrex::iMultiFab>& yflux_imask_lev,
                            std::unique_ptr<amrex::iMultiFab>& zflux_imask_lev,
                            std::unique_ptr<amrex::MultiFab>& thin_xforce_lev,
                            std::unique_ptr<amrex::MultiFab>& thin_yforce_lev,
                            std::unique_ptr<amrex::MultiFab>& thin_zforce_lev);

#if defined(ERF_USE_NETCDF)
void
moist_set_rhs (const amrex::Geometry& geom,
               const amrex::Box& tbx,
               const amrex::Array4<amrex::Real const>& old_cons,
               const amrex::Array4<amrex::Real const>& new_cons,
               const amrex::Array4<amrex::Real      >& cell_rhs,
               const amrex::Real& bdy_time_interval,
               const amrex::Real& new_stage_time,
               const amrex::Real& dt,
               const amrex::Real& stop_time_elapsed,
               int width, int set_width,
               const amrex::Box& domain,
               amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
               amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
               amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
               amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi);
#endif

void ApplySpongeZoneBCsForCC (const SpongeChoice& spongeChoice,
                              const amrex::Geometry geom,
                              const amrex::Box& bx,
                              const amrex::Array4<amrex::Real>& cell_rhs,
                              const amrex::Array4<const amrex::Real>& cell_data,
                              const amrex::Array4<const amrex::Real>& r0,
                              const amrex::Array4<const amrex::Real>& z_phys_cc);

void ApplySpongeZoneBCsForMom (const SpongeChoice& spongeChoice,
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
                               const amrex::Array4<const amrex::Real>& r0,
                               const amrex::Array4<const amrex::Real>& z_phys_nd,
                               const amrex::Array4<const amrex::Real>& z_phys_cc);

void ApplySpongeZoneBCsForMom_ReadFromFile (const SpongeChoice& spongeChoice,
                                              const amrex::Geometry geom,
                                              const amrex::Box& tbx,
                                              const amrex::Box& tby,
                                              const amrex::Array4<const amrex::Real>& cell_data,
                                              const amrex::Array4<const amrex::Real>& z_phys_cc,
                                              const amrex::Array4<amrex::Real>& rho_u_rhs,
                                              const amrex::Array4<amrex::Real>& rho_v_rhs,
                                              const amrex::Array4<const amrex::Real>& rho_u,
                                              const amrex::Array4<const amrex::Real>& rho_v,
                                              const amrex::Vector<amrex::Real*> d_sponge_ptrs_at_lev);

void ApplyBndryForcing_Forecast (const SolverChoice& solverChoice,
                                 const amrex::Geometry geom,
                                 const amrex::Box& tbx,
                                 const amrex::Box& tby,
                                 const amrex::Box& tbz,
                                 const amrex::Array4<const amrex::Real>& z_phys_nd,
                                 const amrex::Array4<amrex::Real>& rho_u_rhs,
                                 const amrex::Array4<amrex::Real>& rho_v_rhs,
                                 const amrex::Array4<amrex::Real>& rho_w_rhs,
                                 const amrex::Array4<const amrex::Real>& rho_u,
                                 const amrex::Array4<const amrex::Real>& rho_v,
                                 const amrex::Array4<const amrex::Real>& rho_w,
                                 const amrex::Array4<const amrex::Real>& rho_u_initial_state,
                                 const amrex::Array4<const amrex::Real>& rho_v_initial_state,
                                 const amrex::Array4<const amrex::Real>& rho_w_initial_state,
                                 const amrex::Array4<const amrex::Real>& cons_initial_state);

#endif

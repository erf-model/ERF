#ifndef ERF_INTEGRATION_H_
#define ERF_INTEGRATION_H_

#include <AMReX_MultiFab.H>
#include <AMReX_EBMultiFabUtil.H>

#include "ERF_DataStruct.H"
#include "ERF_InputSoundingData.H"
#include "ERF_LargeScaleForcingData.H"
#include "ERF_TurbPertStruct.H"
#include "ERF_EB.H"
#include "ERF_ReadBndryPlanes.H"

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

void compute_gradp_xy (const amrex::MultiFab& p,
                       const amrex::Geometry& geom,
                       const amrex::MultiFab& z_phys_cc,
                       amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                       const eb_& ebfact,
                       amrex::Vector<amrex::MultiFab>& gradp,
                       const SolverChoice& solverChoice);

void compute_gradp_z (const amrex::MultiFab& p,
                      const amrex::Geometry& geom,
                      const amrex::MultiFab& z_phys_nd,
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
                   double dt, double time,
                   const amrex::Vector<amrex::MultiFab>& S_data,
                   const amrex::MultiFab& S_prim,
                         amrex::MultiFab& cc_source,
                   const amrex::MultiFab& base_state,
                   const amrex::MultiFab* z_phys_cc,
                   const amrex::MultiFab& xvel,
                   const amrex::MultiFab& yvel,
                   const amrex::MultiFab& zvel,
                   const amrex::MultiFab* qheating_rates,
                   amrex::MultiFab* terrain_blank,
                   const amrex::Geometry geom,
                   const SolverChoice& solverChoice,
                   amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                   const amrex::MultiFab* rhotheta_src,
                   const amrex::MultiFab* rhoqt_src,
                   const amrex::Real* dptr_wbar_sub,
                   const amrex::Vector<amrex::Real*> d_rayleigh_ptrs_at_lev,
                   const amrex::Real* sinesq_at_lev_d,
                   TurbulentPerturbation& turbPert,
                   bool is_slow_step);

#if defined(ERF_USE_NETCDF)
/**
 * Add moisture nudging terms to the source.
 *
 * @param[in] S_data State data.
 * @param[in,out] source Source terms to be updated.
 * @param[in] n_qstate Number of moisture species.
 * @param[in] dt Time step.
 * @param[in] old_stage_time_total Total time at current stage.
 * @param[in] start_bdy_time Start time for boundary data.
 * @param[in] final_bdy_time Final time for boundary data.
 * @param[in] bdy_time_interval Interval between boundary data files.
 * @param[in] bdy_factor Nudging strength factor.
 * @param[in] width Nudging width.
 * @param[in] geom Grid geometry.
 * @param[in] bdy_data_xlo Boundary data for x-low face.
 * @param[in] bdy_data_xhi Boundary data for x-high face.
 * @param[in] bdy_data_ylo Boundary data for y-low face.
 * @param[in] bdy_data_yhi Boundary data for y-high face.
 * @param[in] m_r2d Pointer to boundary plane reader.
 * @param[in] c_p Specific heat at constant pressure.
 * @param[in] rdOcp Ratio of gas constants to specific heat.
 * @param[in] use_wrf_bdy_density Use WRF boundary density.
 * @param[in] bdy_moist_nudge_type Type of moisture nudging.
 */
void
add_moist_nudging_terms (const amrex::MultiFab& S_data,
                               amrex::MultiFab& source,
                         const int n_qstate,
                         const double& dt,
                         const double& old_stage_time_total,
                         const double& start_bdy_time,
                         const double& final_bdy_time,
                         const double& bdy_time_interval,
                         const amrex::Real& bdy_factor,
                         const int width,
                         const amrex::Geometry& geom,
                         amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
                         amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
                         amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
                         amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi,
                         std::unique_ptr<ReadBndryPlanes>& m_r2d,
                         const amrex::Real& c_p,
                         const amrex::Real& rdOcp,
                         const bool use_wrf_bdy_density,
                         const int bdy_moist_nudge_type);
#endif

void make_mom_sources (double time,
                       double dt,
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
                             amrex::MultiFab* terrain_blank_xface,
                             amrex::MultiFab* terrain_blank_yface,
                             amrex::MultiFab* terrain_blank_zface,
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
                             LargeScaleForcingData &lsf_data,
                       std::unique_ptr<amrex::MultiFab>& lsf_tendencies,
                       const eb_& ebfact,
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
/**
 * Set the right-hand side for moist variables.
 *
 * @param[in] geom Grid geometry.
 * @param[in] tbx Local box.
 * @param[in] new_cons Current conservative variables.
 * @param[in,out] cell_rhs Right-hand side to be updated.
 * @param[in] n_qstate Number of moisture species.
 * @param[in] time Current simulation time.
 * @param[in] dt Time step.
 * @param[in] start_bdy_time Start time for boundary data.
 * @param[in] final_bdy_time Final time for boundary data.
 * @param[in] bdy_time_interval Interval between boundary data files.
 * @param[in] nudge_factor Nudging strength factor.
 * @param[in] width Nudging width.
 * @param[in] domain Global domain box.
 * @param[in] bdy_data_xlo Boundary data for x-low face.
 * @param[in] bdy_data_xhi Boundary data for x-high face.
 * @param[in] bdy_data_ylo Boundary data for y-low face.
 * @param[in] bdy_data_yhi Boundary data for y-high face.
 * @param[in] m_r2d Pointer to boundary plane reader.
 */
void
moist_set_rhs (const amrex::Geometry& geom,
               const amrex::Box& tbx,
               const amrex::Array4<amrex::Real const>& new_cons,
               const amrex::Array4<amrex::Real      >& cell_rhs,
               const int n_qstate,
               const double& time,
               const double& dt,
               const double& start_bdy_time,
               const double& final_bdy_time,
               const double& bdy_time_interval,
               const amrex::Real& nudge_factor,
               const int width,
               const amrex::Box& domain,
               amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
               amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
               amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
               amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi,
               std::unique_ptr<ReadBndryPlanes>& m_r2d);
#endif

void ApplySpongeZoneBCsForCC (const SpongeChoice& spongeChoice,
                              const amrex::Geometry geom,
                              const amrex::Box& bx,
                              const amrex::Array4<amrex::Real>& cell_rhs,
                              const amrex::Array4<const amrex::Real>& cell_data,
                              const amrex::Array4<const amrex::Real>& r0,
                              const amrex::Array4<const amrex::Real>& th0,
                              const amrex::Array4<const amrex::Real>& qv0,
                              const amrex::Array4<const amrex::Real>& z_phys_cc,
                              int n_qstate);

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

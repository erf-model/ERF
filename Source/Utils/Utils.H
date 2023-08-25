#ifndef _UTILS_H_
#define _UTILS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_InterpFaceRegister.H>
#include <DataStruct.H>
#include <IndexDefines.H>
#include <ABLMost.H>

/*
 * Create the Jacobian for the metric transformation when use_terrain is true
 */
void make_J   (const amrex::Geometry& geom,
                     amrex::MultiFab& z_phys_nd,
                     amrex::MultiFab& detJ_cc);

/*
 * Average z_phys_nd on nodes to cell centers
 */
void make_zcc (const amrex::Geometry& geom,
                     amrex::MultiFab& z_phys_nd,
                     amrex::MultiFab& z_phys_cc);

/*
 * Convert momentum to velocity by dividing by density averaged onto faces
 */
void MomentumToVelocity (amrex::MultiFab& xvel_out,
                         amrex::MultiFab& yvel_out,
                         amrex::MultiFab& zvel_out,
                         const amrex::MultiFab& cons_in,
                         const amrex::MultiFab& xmom_in,
                         const amrex::MultiFab& ymom_in,
                         const amrex::MultiFab& zmom_in);

/*
 * Convert velocity to momentum by multiplying by density averaged onto faces
 */
void VelocityToMomentum (const amrex::MultiFab& xvel_in,
                         const amrex::IntVect& xvel_ngrow,
                         const amrex::MultiFab& yvel_in,
                         const amrex::IntVect& yvel_ngrow,
                         const amrex::MultiFab& zvel_in,
                         const amrex::IntVect& zvel_ngrow,
                         const amrex::MultiFab& cons_in,
                         amrex::MultiFab& xmom_out,
                         amrex::MultiFab& ymom_out,
                         amrex::MultiFab& zmom_out,
                         bool l_use_ndiff);

/*
 * Compute boxes for looping over interior/exterior ghost cells
 * for use by fillpatch, erf_slow_rhs_pre, and erf_slow_rhs_post
 */
void compute_interior_ghost_bxs_xy(const amrex::Box& bx,
                                   const amrex::Box& domain,
                                   const int& width,
                                   const int& set_width,
                                   amrex::Box& bx_xlo,
                                   amrex::Box& bx_xhi,
                                   amrex::Box& bx_ylo,
                                   amrex::Box& bx_yhi,
                                   const amrex::IntVect& ng_vect=amrex::IntVect(0,0,0),
                                   const bool get_int_ng=false);

/*
 * Compute boxes to be used for specifying where in the relaxation zone to compute the RHS
 */
void compute_interior_ghost_RHS(const std::string& init_type,
                                const amrex::Real& bdy_time_interval,
                                const amrex::Real& time,
                                const amrex::Real& delta_t,
                                const int& width,
                                const int& set_width,
                                const amrex::Geometry& geom,
                                amrex::Vector<amrex::MultiFab>& S_rhs,
                                amrex::Vector<amrex::MultiFab>& S_data,
                                amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
                                amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
                                amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
                                amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi,
                                const amrex::Real start_bdy_time);

/*
 * Update the variables in the relaxation zone
 */
void update_interior_ghost(const amrex::Real& delta_t,
                           const int& width,
                           const amrex::Geometry& geom,
                           amrex::Vector<amrex::MultiFab>& S_rhs,
                           const amrex::Vector<amrex::MultiFab>& S_old,
                           amrex::Vector<amrex::MultiFab>& S_data);
#endif

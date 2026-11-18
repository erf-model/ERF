#ifndef ERF_IMMERSEDFORCING_H_
#define ERF_IMMERSEDFORCING_H_

#include <AMReX.H>
#include <ERF_DataStruct.H>
#include <AMReX_MultiFab.H>
#include <AMReX_TableData.H>

// Terrain immersed forcing - MOMENTUM
void ImmersedForcingTerrain_Xmom (const amrex::Box& tbx,
                                  const amrex::Array4<const amrex::Real>& u,
                                  const amrex::Array4<const amrex::Real>& v,
                                  const amrex::Array4<const amrex::Real>& w,
                                  const amrex::Array4<const amrex::Real>& cell_data,
                                  const amrex::Array4<const amrex::Real>& t_blank_arr,
                                  const amrex::Array4<const amrex::Real>& t_blank_xface_arr,
                                  const amrex::Array4<const amrex::Real>& z_cc_arr,
                                  const amrex::Array4<      amrex::Real>& xmom_src_arr,
                                  const amrex::Geometry& geom,
                                  const SolverChoice& solverChoice,
                                  const amrex::Real fac);

void ImmersedForcingTerrain_Ymom (const amrex::Box& tby,
                                  const amrex::Array4<const amrex::Real>& u,
                                  const amrex::Array4<const amrex::Real>& v,
                                  const amrex::Array4<const amrex::Real>& w,
                                  const amrex::Array4<const amrex::Real>& cell_data,
                                  const amrex::Array4<const amrex::Real>& t_blank_arr,
                                  const amrex::Array4<const amrex::Real>& t_blank_yface_arr,
                                  const amrex::Array4<const amrex::Real>& z_cc_arr,
                                  const amrex::Array4<      amrex::Real>& ymom_src_arr,
                                  const amrex::Geometry& geom,
                                  const SolverChoice& solverChoice,
                                  const amrex::Real fac);

void ImmersedForcingTerrain_Zmom (const amrex::Box& tbz,
                                  const amrex::Array4<const amrex::Real>& u,
                                  const amrex::Array4<const amrex::Real>& v,
                                  const amrex::Array4<const amrex::Real>& w,
                                  const amrex::Array4<const amrex::Real>& cell_data,
                                  const amrex::Array4<const amrex::Real>& t_blank_arr,
                                  const amrex::Array4<const amrex::Real>& t_blank_zface_arr,
                                  const amrex::Array4<const amrex::Real>& z_cc_arr,
                                  const amrex::Array4<      amrex::Real>& zmom_src_arr,
                                  const amrex::Geometry& geom,
                                  const SolverChoice& solverChoice,
                                  const amrex::Real fac);

// Buildings immersed forcing - MOMENTUM
void ImmersedForcingBuildings_Xmom (const amrex::Box& tbx,
                                    const amrex::Array4<const amrex::Real>& u,
                                    const amrex::Array4<const amrex::Real>& v,
                                    const amrex::Array4<const amrex::Real>& w,
                                    const amrex::Array4<const amrex::Real>& cell_data,
                                    const amrex::Array4<const amrex::Real>& t_blank_arr,
                                    const amrex::Array4<const amrex::Real>& t_blank_xface_arr,
                                    const amrex::Array4<const amrex::Real>& z_cc_arr,
                                    const amrex::Array4<      amrex::Real>& xmom_src_arr,
                                    const amrex::Geometry& geom,
                                    const SolverChoice& solverChoice,
                                    const amrex::Real fac,
                                    bool is_slow_step);

void ImmersedForcingBuildings_Ymom (const amrex::Box& tby,
                                    const amrex::Array4<const amrex::Real>& u,
                                    const amrex::Array4<const amrex::Real>& v,
                                    const amrex::Array4<const amrex::Real>& w,
                                    const amrex::Array4<const amrex::Real>& cell_data,
                                    const amrex::Array4<const amrex::Real>& t_blank_arr,
                                    const amrex::Array4<const amrex::Real>& t_blank_yface_arr,
                                    const amrex::Array4<const amrex::Real>& z_cc_arr,
                                    const amrex::Array4<      amrex::Real>& ymom_src_arr,
                                    const amrex::Geometry& geom,
                                    const SolverChoice& solverChoice,
                                    const amrex::Real fac,
                                    bool is_slow_step);

void ImmersedForcingBuildings_Zmom (const amrex::Box& tbz,
                                    const amrex::Array4<const amrex::Real>& u,
                                    const amrex::Array4<const amrex::Real>& v,
                                    const amrex::Array4<const amrex::Real>& w,
                                    const amrex::Array4<const amrex::Real>& cell_data,
                                    const amrex::Array4<const amrex::Real>& t_blank_arr,
                                    const amrex::Array4<const amrex::Real>& t_blank_zface_arr,
                                    const amrex::Array4<const amrex::Real>& z_cc_arr,
                                    const amrex::Array4<      amrex::Real>& zmom_src_arr,
                                    const amrex::Geometry& geom,
                                    const SolverChoice& solverChoice,
                                    const amrex::Real fac,
                                    bool is_slow_step);

// Terrain immersed forcing - MOMENTUM, fraction-stress wall law (erf.if_wall_form = fraction_stress)
void ImmersedForcingTerrain_Xmom_FractionStress (const amrex::Box& tbx,
                                                 const amrex::Array4<const amrex::Real>& u,
                                                 const amrex::Array4<const amrex::Real>& v,
                                                 const amrex::Array4<const amrex::Real>& w,
                                                 const amrex::Array4<const amrex::Real>& cell_data,
                                                 const amrex::Array4<const amrex::Real>& t_blank_arr,
                                                 const amrex::Array4<const amrex::Real>& t_blank_xface_arr,
                                                 const amrex::Array4<      amrex::Real>& xmom_src_arr,
                                                 const amrex::Geometry& geom,
                                                 const SolverChoice& solverChoice,
                                                 const amrex::Real fac,
                                                 const amrex::Real time);

void ImmersedForcingTerrain_Ymom_FractionStress (const amrex::Box& tby,
                                                 const amrex::Array4<const amrex::Real>& u,
                                                 const amrex::Array4<const amrex::Real>& v,
                                                 const amrex::Array4<const amrex::Real>& w,
                                                 const amrex::Array4<const amrex::Real>& cell_data,
                                                 const amrex::Array4<const amrex::Real>& t_blank_arr,
                                                 const amrex::Array4<const amrex::Real>& t_blank_yface_arr,
                                                 const amrex::Array4<      amrex::Real>& ymom_src_arr,
                                                 const amrex::Geometry& geom,
                                                 const SolverChoice& solverChoice,
                                                 const amrex::Real fac,
                                                 const amrex::Real time);

// Terrain immersed forcing - SCALARS, fraction-stress wall law (erf.if_wall_form = fraction_stress)
void ImmersedForcingTerrain_Scalar_FractionStress (const amrex::Box& bx,
                                                   const amrex::Array4<const amrex::Real>& u,
                                                   const amrex::Array4<const amrex::Real>& v,
                                                   const amrex::Array4<const amrex::Real>& cell_data,
                                                   const amrex::Array4<const amrex::Real>& t_blank_arr,
                                                   const amrex::Array4<      amrex::Real>& cell_src,
                                                   const amrex::Array4<      amrex::Real>& wall_hfx,
                                                   const amrex::Array4<      amrex::Real>& wall_ktarget,
                                                   const amrex::Geometry& geom,
                                                   const SolverChoice& solverChoice,
                                                   const amrex::Table1D<amrex::Real>& r_avg,
                                                   const amrex::Table1D<amrex::Real>& t_avg,
                                                   const amrex::Real dt,
                                                   const amrex::Real time,
                                                   bool wall_tke = false,
                                                   amrex::Real Cmu0 = amrex::Real(0.5562),
                                                   const amrex::Array4<      amrex::Real>& wall_qfx = amrex::Array4<amrex::Real>{},
                                                   bool wall_tke_moist = false);

// Terrain immersed forcing - SCALARS
void ImmersedForcingTerrain_Scalar (const amrex::Box& bx,
                                   const amrex::Array4<const amrex::Real>& u,
                                   const amrex::Array4<const amrex::Real>& v,
                                   const amrex::Array4<const amrex::Real>& cell_data,
                                   const amrex::Array4<const amrex::Real>& t_blank_arr,
                                   const amrex::Array4<const amrex::Real>& z_cc_arr,
                                   const amrex::Array4<      amrex::Real>& cell_src,
                                   const amrex::Geometry& geom,
                                   const SolverChoice& solverChoice,
                                   const amrex::Table1D<amrex::Real>& r_avg,
                                   const amrex::Table1D<amrex::Real>& t_avg,
                                   const amrex::Real time);

// Buildings immersed forcing - SCALARS
void ImmersedForcingBuildings_Scalar (const amrex::Box& bx,
                                     const amrex::Array4<const amrex::Real>& u,
                                     const amrex::Array4<const amrex::Real>& v,
                                     const amrex::Array4<const amrex::Real>& w,
                                     const amrex::Array4<const amrex::Real>& cell_data,
                                     const amrex::Array4<const amrex::Real>& t_blank_arr,
                                     const amrex::Array4<const amrex::Real>& z_cc_arr,
                                     const amrex::Array4<      amrex::Real>& cell_src,
                                     const amrex::Geometry& geom,
                                     const SolverChoice& solverChoice,
                                     const amrex::Table1D<amrex::Real>& r_avg,
                                     const amrex::Table1D<amrex::Real>& t_avg,
                                     const amrex::Real time);

#endif

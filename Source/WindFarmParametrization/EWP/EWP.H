#ifndef EWP_H
#define EWP_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

void ewp_advance (int lev,
                    const amrex::Geometry& geom,
                    const amrex::Real& dt_advance,
                    amrex::MultiFab& cons_in,
                    amrex::MultiFab& U_old, amrex::MultiFab& V_old, amrex::MultiFab& W_old,
                    amrex::MultiFab& mf_vars_ewp, const amrex::MultiFab& mf_Nturb);

void ewp_source_terms_cellcentered (const amrex::Geometry& geom,
                                      const amrex::MultiFab& cons_in,
                                      const amrex::MultiFab& U_old, const amrex::MultiFab& V_old, const amrex::MultiFab& W_old,
                                      amrex::MultiFab& mf_vars_ewp, const amrex::MultiFab& mf_Nturb);

void ewp_update (const amrex::Real& dt_advance,
                  amrex::MultiFab& cons_in,
                  amrex::MultiFab& U_old, amrex::MultiFab& V_old,
                  const amrex::MultiFab& mf_vars_ewp);

#endif


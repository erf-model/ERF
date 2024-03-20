#ifndef FITCH_H
#define FITCH_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

void fitch_advance (int lev,
                    const amrex::Geometry& geom,
                    const amrex::Real& dt_advance,
                    amrex::MultiFab& cons_in,
                    amrex::MultiFab& U_old, amrex::MultiFab& V_old, amrex::MultiFab& W_old,
                    amrex::MultiFab& mf_vars_fitch, const amrex::MultiFab& Nturb);

void fitch_source_terms_cellcentered (const amrex::Geometry& geom,
                                      const amrex::MultiFab& cons_in,
                                      const amrex::MultiFab& U_old, const amrex::MultiFab& V_old, const amrex::MultiFab& W_old,
                                      amrex::MultiFab& mf_vars_fitch, const amrex::MultiFab& Nturb);

void fitch_update (const amrex::Real& dt_advance,
                  amrex::MultiFab& cons_in,
                  amrex::MultiFab& U_old, amrex::MultiFab& V_old,
                  const amrex::MultiFab& mf_vars_fitch);

#endif


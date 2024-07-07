#ifndef SIMPLEAD_H
#define SIMPLEAD_H

#include "WindFarm.H"
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

class SimpleAD : public WindFarm {

public:

    void simpleAD_advance (int lev,
                        const amrex::Geometry& geom,
                        const amrex::Real& dt_advance,
                        amrex::MultiFab& cons_in,
                        amrex::MultiFab& U_old, amrex::MultiFab& V_old, amrex::MultiFab& W_old,
                        amrex::MultiFab& mf_vars_ewp, const amrex::MultiFab& mf_Nturb);

    void simpleAD_source_terms_cellcentered (const amrex::Geometry& geom,
                                          const amrex::MultiFab& cons_in,
                                          const amrex::MultiFab& U_old, const amrex::MultiFab& V_old, const amrex::MultiFab& W_old,
                                          amrex::MultiFab& mf_vars_ewp, const amrex::MultiFab& mf_Nturb);

    void simpleAD_update (const amrex::Real& dt_advance,
                      amrex::MultiFab& cons_in,
                      amrex::MultiFab& U_old, amrex::MultiFab& V_old,
                      const amrex::MultiFab& mf_vars_ewp);
};

#endif


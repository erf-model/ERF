#ifndef ERF_FITCH_H
#define ERF_FITCH_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "ERF_NullWindFarm.H"

class Fitch : public NullWindFarm {

public:

    Fitch() {}

    virtual ~Fitch() = default;

    void advance (const amrex::Geometry& geom,
                  const amrex::Real& dt_advance,
                  amrex::MultiFab& cons_in,
                  amrex::MultiFab& mf_vars_fitch,
                  amrex::MultiFab& U_old,
                  amrex::MultiFab& V_old,
                  amrex::MultiFab& W_old,
                  const amrex::MultiFab& mf_Nturb) override;

    void source_terms_cellcentered (const amrex::Geometry& geom,
                                    const amrex::MultiFab& cons_in,
                                    amrex::MultiFab& mf_vars_ewp,
                                    const amrex::MultiFab& U_old,
                                    const amrex::MultiFab& V_old,
                                    const amrex::MultiFab& W_old,
                                    const amrex::MultiFab& mf_Nturb);

    void update (const amrex::Real& dt_advance,
                  amrex::MultiFab& cons_in,
                  amrex::MultiFab& U_old, amrex::MultiFab& V_old,
                  const amrex::MultiFab& mf_vars_fitch);

protected:
    amrex::Vector<amrex::Real> xloc, yloc;
    amrex::Real hub_height, rotor_rad, thrust_coeff_standing, nominal_power;
    amrex::Vector<amrex::Real> wind_speed, thrust_coeff, power;
};

#endif


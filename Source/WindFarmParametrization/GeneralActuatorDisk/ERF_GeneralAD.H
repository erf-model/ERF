#ifndef ERF_SIMPLEAD_H
#define ERF_SIMPLEAD_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include "ERF_NullWindFarm.H"

class GeneralAD : public NullWindFarm {

public:

    GeneralAD() {}

    virtual ~GeneralAD() = default;

    void advance (const amrex::Geometry& geom,
                  const amrex::Real& dt_advance,
                  amrex::MultiFab& cons_in,
                  amrex::MultiFab& mf_vars_windfarm,
                  amrex::MultiFab& U_old,
                  amrex::MultiFab& V_old,
                  amrex::MultiFab& W_old,
                  const amrex::MultiFab& mf_Nturb,
                  const amrex::MultiFab& mf_SMark) override;

    void compute_freestream_velocity(const amrex::MultiFab& cons_in,
                                     const amrex::MultiFab& U_old,
                                     const amrex::MultiFab& V_old,
                                     const amrex::MultiFab& mf_SMark);

    void source_terms_cellcentered (const amrex::Geometry& geom,
                                    const amrex::MultiFab& cons_in,
                                    const amrex::MultiFab& mf_Smark,
                                    amrex::MultiFab& mf_vars_simpleAD);

    void update (const amrex::Real& dt_advance,
                 amrex::MultiFab& cons_in,
                 amrex::MultiFab& U_old,
                 amrex::MultiFab& V_old,
                 const amrex::MultiFab& mf_vars);

protected:
    amrex::Vector<amrex::Real> xloc, yloc;
    amrex::Real turb_disk_angle;
    amrex::Real hub_height, rotor_rad, thrust_coeff_standing, nominal_power;
    amrex::Vector<amrex::Real> wind_speed, thrust_coeff, power;
    amrex::Vector<amrex::Real> freestream_velocity, freestream_phi, disk_cell_count;
};

#endif


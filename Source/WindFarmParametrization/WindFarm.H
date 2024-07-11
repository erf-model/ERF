#ifndef WINDFARM_H
#define WINDFARM_H

#include <DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>

extern amrex::Real hub_height, rotor_rad, thrust_coeff_standing, nominal_power;
extern amrex::Vector<amrex::Real> wind_speed, thrust_coeff, power;
extern amrex::Gpu::DeviceVector<amrex::Real> d_wind_speed, d_thrust_coeff, d_power;

void read_in_table(std::string windfarm_spec_table);

void advance_windfarm (int lev,
                       const amrex::Geometry& geom,
                       const amrex::Real& dt_advance,
                       amrex::MultiFab& cons_in,
                       amrex::MultiFab& U_old, amrex::MultiFab& V_old, amrex::MultiFab& W_old,
                       amrex::MultiFab& mf_vars_windfarm, const amrex::MultiFab& mf_Nturb,
                       SolverChoice& solver_choice);
#endif


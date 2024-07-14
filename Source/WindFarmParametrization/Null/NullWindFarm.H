#ifndef NULLWINDFARM_H
#define NULLWINDFARM_H

#include <DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>

class NullWindFarm {

public:

    NullWindFarm()
	{
		rotor_rad = 10000.0;
	}

	virtual void advance (int lev,
                        const amrex::Geometry& geom,
                        const amrex::Real& dt_advance,
                        amrex::MultiFab& cons_in,
                        amrex::MultiFab& U_old, amrex::MultiFab& V_old, amrex::MultiFab& W_old,
                        amrex::MultiFab& mf_vars_ewp, const amrex::MultiFab& mf_Nturb) = 0;

    /*virtual void source_terms_cellcentered (const amrex::Geometry& geom,
                                          const amrex::MultiFab& cons_in,
                                          const amrex::MultiFab& U_old, const amrex::MultiFab& V_old, const amrex::MultiFab& W_old,
                                          amrex::MultiFab& mf_vars_ewp, const amrex::MultiFab& mf_Nturb) = 0;

    virtual void update (const amrex::Real& dt_advance,
                      amrex::MultiFab& cons_in,
                      amrex::MultiFab& U_old, amrex::MultiFab& V_old,
                      const amrex::MultiFab& mf_vars_ewp) = 0;*/


protected:

    amrex::Vector<amrex::Real> xloc, yloc;
    amrex::Real hub_height, rotor_rad, thrust_coeff_standing, nominal_power;
    amrex::Vector<amrex::Real> wind_speed, thrust_coeff, power;
    amrex::Gpu::DeviceVector<amrex::Real> d_wind_speed, d_thrust_coeff, d_power;
    amrex::Gpu::DeviceVector<amrex::Real> d_xloc, d_yloc;
};


#endif


#ifndef ERF_NULLWINDFARM_H
#define ERF_NULLWINDFARM_H

#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

class NullWindFarm {

public:

    NullWindFarm() {}

    virtual ~NullWindFarm() = default;

    /*virtual void advance (const amrex::Geometry& geom,
                        const amrex::Real& dt_advance,
                        amrex::MultiFab& cons_in,
                        amrex::MultiFab& U_old, amrex::MultiFab& V_old, amrex::MultiFab& W_old,
                        amrex::MultiFab& mf_vars_ewp, const amrex::MultiFab& mf_Nturb) = 0;*/

    virtual void advance (const amrex::Geometry& a_geom,
                  const amrex::Real& dt_advance,
                  amrex::MultiFab& cons_in,
                  amrex::MultiFab& mf_vars_windfarm,
                  amrex::MultiFab& U_old,
                  amrex::MultiFab& V_old,
                  amrex::MultiFab& W_old,
                  const amrex::MultiFab& mf_Nturb) = 0;

    virtual void set_turb_spec(const amrex::Real&  rotor_rad, const amrex::Real& hub_height,
                               const amrex::Real& thrust_coeff_standing, const amrex::Vector<amrex::Real>& wind_speed,
                               const amrex::Vector<amrex::Real>& thrust_coeff,
                               const amrex::Vector<amrex::Real>& power)
    {
        m_rotor_rad = rotor_rad;
        m_hub_height = hub_height;
        m_thrust_coeff_standing = thrust_coeff_standing;
        m_wind_speed = wind_speed;
        m_thrust_coeff = thrust_coeff;
        m_power = power;
    }

    virtual void set_turb_loc (const amrex::Vector<amrex::Real>& xloc,
                               const amrex::Vector<amrex::Real>& yloc)
    {
        m_xloc = xloc;
        m_yloc = yloc;
    }

    void get_turb_spec (amrex::Real& rotor_rad, amrex::Real& hub_height,
                        amrex::Real& thrust_coeff_standing, amrex::Vector<amrex::Real>& wind_speed,
                        amrex::Vector<amrex::Real>& thrust_coeff, amrex::Vector<amrex::Real>& power)
    {
        rotor_rad = m_rotor_rad;
        hub_height = m_hub_height;
        thrust_coeff_standing = m_thrust_coeff_standing;
        wind_speed = m_wind_speed;
        thrust_coeff = m_thrust_coeff;
        power = m_power;
    }

    void get_turb_loc (amrex::Vector<amrex::Real>& xloc,
                       amrex::Vector<amrex::Real>& yloc)
    {
        xloc = m_xloc;
        yloc = m_yloc;
    }

protected:

    amrex::Vector<amrex::Real> m_xloc, m_yloc;
    amrex::Real m_hub_height, m_rotor_rad, m_thrust_coeff_standing, m_nominal_power;
    amrex::Vector<amrex::Real> m_wind_speed, m_thrust_coeff, m_power;
};


#endif


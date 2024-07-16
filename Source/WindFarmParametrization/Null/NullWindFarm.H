#ifndef NULLWINDFARM_H
#define NULLWINDFARM_H

#include <DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>

class NullWindFarm {

public:

    NullWindFarm() {}

    virtual ~NullWindFarm() = default;

    virtual void advance (const amrex::Geometry& geom,
                        const amrex::Real& dt_advance,
                        amrex::MultiFab& cons_in,
                        amrex::MultiFab& U_old, amrex::MultiFab& V_old, amrex::MultiFab& W_old,
                        amrex::MultiFab& mf_vars_ewp, const amrex::MultiFab& mf_Nturb) = 0;

    virtual void set_turb_spec(const amrex::Real&  rotor_rad, const amrex::Real& hub_height,
                               const amrex::Real& thrust_coeff_standing, const amrex::Gpu::DeviceVector<amrex::Real>& d_wind_speed,
                               const amrex::Gpu::DeviceVector<amrex::Real>& d_thrust_coeff,
                               const amrex::Gpu::DeviceVector<amrex::Real>& d_power)
    {
        m_rotor_rad = rotor_rad;
        m_hub_height = hub_height;
        m_thrust_coeff_standing = thrust_coeff_standing;
        m_d_wind_speed = d_wind_speed;
        m_d_thrust_coeff = d_thrust_coeff;
        m_d_power = d_power;
    }

    virtual void set_turb_loc (const amrex::Gpu::DeviceVector<amrex::Real>& d_xloc,
                               const amrex::Gpu::DeviceVector<amrex::Real>& d_yloc)
    {
        m_d_xloc = d_xloc;
        m_d_yloc = d_yloc;
    }

    void get_turb_spec (amrex::Real& rotor_rad, amrex::Real& hub_height,
                        amrex::Real& thrust_coeff_standing, amrex::Gpu::DeviceVector<amrex::Real>& d_wind_speed,
                        amrex::Gpu::DeviceVector<amrex::Real>& d_thrust_coeff, amrex::Gpu::DeviceVector<amrex::Real>& d_power)
    {
        rotor_rad = m_rotor_rad;
        hub_height = m_hub_height;
        thrust_coeff_standing = m_thrust_coeff_standing;
        d_wind_speed = m_d_wind_speed;
        d_thrust_coeff = m_d_thrust_coeff;
        d_power = m_d_power;
    }

    void get_turb_loc (amrex::Gpu::DeviceVector<amrex::Real>& d_xloc,
                       amrex::Gpu::DeviceVector<amrex::Real>& d_yloc)
    {
        d_xloc = m_d_xloc;
        d_yloc = m_d_yloc;
    }

protected:

    amrex::Vector<amrex::Real> m_xloc, m_yloc;
    amrex::Real m_hub_height, m_rotor_rad, m_thrust_coeff_standing, m_nominal_power;
    amrex::Vector<amrex::Real> m_wind_speed, m_thrust_coeff, m_power;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_wind_speed, m_d_thrust_coeff, m_d_power;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_xloc, m_d_yloc;
};


#endif


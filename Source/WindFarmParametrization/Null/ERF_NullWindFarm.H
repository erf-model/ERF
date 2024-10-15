#ifndef ERF_NULLWINDFARM_H
#define ERF_NULLWINDFARM_H

#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>

class NullWindFarm {

public:

    NullWindFarm() {}

    virtual ~NullWindFarm() = default;

    virtual void advance (const amrex::Geometry& a_geom,
                  const amrex::Real& dt_advance,
                  amrex::MultiFab& cons_in,
                  amrex::MultiFab& mf_vars_windfarm,
                  amrex::MultiFab& U_old,
                  amrex::MultiFab& V_old,
                  amrex::MultiFab& W_old,
                  const amrex::MultiFab& mf_Nturb,
                  const amrex::MultiFab& mf_SMark) = 0;

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

    virtual void set_turb_disk_angle (const amrex::Real& turb_disk_angle)
    {
        m_turb_disk_angle = turb_disk_angle;
    }

    virtual void set_blade_spec(const amrex::Vector<amrex::Real>& bld_rad_loc,
                                const amrex::Vector<amrex::Real>& bld_twist,
                                const amrex::Vector<amrex::Real>& bld_chord)
    {
        m_bld_rad_loc = bld_rad_loc;
        m_bld_twist = bld_twist;
        m_bld_chord = bld_chord;
    }

    virtual void set_blade_airfoil_spec(const amrex::Vector<amrex::Vector<amrex::Real>>& bld_airfoil_aoa,
                                        const amrex::Vector<amrex::Vector<amrex::Real>>& bld_airfoil_Cl,
                                        const amrex::Vector<amrex::Vector<amrex::Real>>& bld_airfoil_Cd)
    {
        m_bld_airfoil_aoa = bld_airfoil_aoa;
        m_bld_airfoil_Cl = bld_airfoil_Cl;
        m_bld_airfoil_Cd = bld_airfoil_Cd;
    }

    virtual void set_turb_spec_extra(const amrex::Vector<amrex::Real>& velocity,
                              const amrex::Vector<amrex::Real>& C_P,
                              const amrex::Vector<amrex::Real>& C_T,
                              const amrex::Vector<amrex::Real>& rotor_RPM,
                              const amrex::Vector<amrex::Real>& blade_pitch)
    {
        m_velocity = velocity;
        m_C_P = C_P;
        m_C_T = C_T;
        m_rotor_RPM = rotor_RPM;
        m_blade_pitch = blade_pitch;
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

    void get_turb_disk_angle (amrex::Real& turb_disk_angle)
    {
        turb_disk_angle = m_turb_disk_angle;
    }

    void get_blade_spec(amrex::Vector<amrex::Real>& bld_rad_loc,
                        amrex::Vector<amrex::Real>& bld_twist,
                        amrex::Vector<amrex::Real>& bld_chord)
    {
        bld_rad_loc = m_bld_rad_loc;
        bld_twist = m_bld_twist;
        bld_chord = m_bld_chord;
    }

    void get_blade_airfoil_spec(amrex::Vector<amrex::Vector<amrex::Real>>& bld_airfoil_aoa,
                                amrex::Vector<amrex::Vector<amrex::Real>>& bld_airfoil_Cl,
                                amrex::Vector<amrex::Vector<amrex::Real>>& bld_airfoil_Cd)
    {
        bld_airfoil_aoa = m_bld_airfoil_aoa;
        bld_airfoil_Cl = m_bld_airfoil_Cl;
        bld_airfoil_Cd = m_bld_airfoil_Cd;
    }

    void get_turb_spec_extra(amrex::Vector<amrex::Real>& velocity,
                             amrex::Vector<amrex::Real>& C_P,
                                amrex::Vector<amrex::Real>& C_T,
                             amrex::Vector<amrex::Real>& rotor_RPM,
                             amrex::Vector<amrex::Real>& blade_pitch)
    {
        velocity = m_velocity;
        C_P = m_C_P;
        C_T = m_C_T;
        rotor_RPM = m_rotor_RPM;
        blade_pitch = m_blade_pitch;
    }

static AMREX_GPU_DEVICE
bool find_if_marked(amrex::Real x1, amrex::Real x2, amrex::Real y1, amrex::Real y2,
                    amrex::Real x0, amrex::Real y0, amrex::Real nx, amrex::Real ny,
                    amrex::Real d_hub_height, amrex::Real d_rotor_rad,
                    amrex::Real z)
{

    // Plane is (x-x0)*nx + (y-y0)*ny = 0. And the planes to intersect are
    // x = x1, x = x2, y = y1, y = y2. The intersections are
    // (xval1,y1), (xval2,y2), (x1,yval1) and (x2,yval2)

    amrex::Real xval1 = x0 - (y1-y0)*ny/(nx+1e-10);
    amrex::Real xval2 = x0 - (y2-y0)*ny/(nx+1e-10);

    amrex::Real yval1 = y0 - (x1-x0)*nx/(ny+1e-10);
    amrex::Real yval2 = y0 - (x2-x0)*nx/(ny+1e-10);

    if( xval1 >=x1 and xval1 <=x2 ) {
        if(std::pow((xval1-x0)*(xval1-x0) + (y1-y0)*(y1-y0) + (z-d_hub_height)*(z-d_hub_height),0.5) < d_rotor_rad ) {
            return true;
        }
    }
    if( xval2 >=x1 and xval2 <=x2 ) {
        if(std::pow((xval2-x0)*(xval2-x0) + (y2-y0)*(y2-y0) + (z-d_hub_height)*(z-d_hub_height),0.5) < d_rotor_rad ) {
            return true;
        }
    }
    if (yval1 >=y1 and yval1 <=y2) {
        if(std::pow((x1-x0)*(x1-x0) + (yval1-y0)*(yval1-y0) + (z-d_hub_height)*(z-d_hub_height),0.5) < d_rotor_rad ) {
            return true;
        }
    }
    if (yval2 >=y1 and yval2 <=y2) {
        if(std::pow((x2-x0)*(x2-x0) + (yval2-y0)*(yval2-y0) + (z-d_hub_height)*(z-d_hub_height),0.5) < d_rotor_rad ) {
            return true;
        }
    }

    return false;
}

protected:

    amrex::Vector<amrex::Real> m_xloc, m_yloc;
    amrex::Real m_turb_disk_angle;
    amrex::Real m_hub_height, m_rotor_rad, m_thrust_coeff_standing, m_nominal_power;
    amrex::Vector<amrex::Real> m_wind_speed, m_thrust_coeff, m_power;
    amrex::Vector<amrex::Real> m_bld_rad_loc, m_bld_twist, m_bld_chord;
    amrex::Vector<amrex::Vector<amrex::Real>> m_bld_airfoil_aoa, m_bld_airfoil_Cl, m_bld_airfoil_Cd;
    amrex::Vector<amrex::Real> m_velocity, m_C_P, m_C_T, m_rotor_RPM, m_blade_pitch;
};


#endif


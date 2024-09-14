#ifndef WINDFARM_H
#define WINDFARM_H

#include <DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "NullWindFarm.H"
#include "Fitch.H"
#include "EWP.H"
#include "SimpleAD.H"

class WindFarm : public NullWindFarm {

public:

    WindFarm(){}

    virtual ~WindFarm() = default;

    WindFarm (int nlev,
              const WindFarmType& a_windfarm_type)
    {
        m_windfarm_model.resize(nlev);
        if (a_windfarm_type == WindFarmType::Fitch) {
            SetModel<Fitch>();
        }
        else if (a_windfarm_type == WindFarmType::EWP) {
            SetModel<EWP>();
        }
        else if (a_windfarm_type == WindFarmType::SimpleAD) {
            SetModel<SimpleAD>();
            amrex::Print() << "Simple actuator disk windfarm model!\n";
        } else {
            amrex::Abort("WindFarm: Dont know this windfarm_type!") ;
        }
    }

    void read_tables (std::string windfarm_loc_table,
                       std::string windfarm_spec_table,
                       bool x_y, bool lat_lon);

    void init_windfarm_lat_lon(const std::string windfarm_loc_table);

    void init_windfarm_x_y(const std::string windfarm_loc_table);

    void read_windfarm_locations_table(const std::string windfarm_loc_table,
                                          bool x_y, bool lat_lon);

    void read_windfarm_spec_table(const std::string windfarm_spec_table);

    void fill_Nturb_multifab(const amrex::Geometry& geom,
                             amrex::MultiFab& mf_Nturb);

    void fill_SMark_multifab(const amrex::Geometry& geom,
                             amrex::MultiFab& mf_SMark,
                             const amrex::Real& sampling_distance_by_D,
                             const amrex::Real& turb_disk_angle);

    void write_turbine_locations_vtk();

    void write_actuator_disks_vtk(const amrex::Geometry& geom);

    void advance (const amrex::Geometry& a_geom,
                  const amrex::Real& dt_advance,
                  amrex::MultiFab& cons_in,
                  amrex::MultiFab& mf_vars_windfarm,
                  amrex::MultiFab& U_old,
                  amrex::MultiFab& V_old,
                  amrex::MultiFab& W_old,
                  const amrex::MultiFab& mf_Nturb,
                  const amrex::MultiFab& mf_SMark) override
    {
        m_windfarm_model[0]->advance(a_geom, dt_advance, cons_in, mf_vars_windfarm,
                                     U_old, V_old, W_old, mf_Nturb, mf_SMark);
    }

    void set_turb_spec(const amrex::Real& a_rotor_rad, const amrex::Real& a_hub_height,
                       const amrex::Real& a_thrust_coeff_standing, const amrex::Vector<amrex::Real>& a_wind_speed,
                       const amrex::Vector<amrex::Real>& a_thrust_coeff,
                       const amrex::Vector<amrex::Real>& a_power) override
    {
        m_windfarm_model[0]->set_turb_spec(a_rotor_rad, a_hub_height, a_thrust_coeff_standing,
                  a_wind_speed, a_thrust_coeff, a_power);
    }

    void set_turb_loc (const amrex::Vector<amrex::Real>& a_xloc,
                       const amrex::Vector<amrex::Real>& a_yloc) override
    {
        m_windfarm_model[0]->set_turb_loc(a_xloc, a_yloc);
    }

    void set_turb_disk_angle (const amrex::Real& turb_disk_angle) override
    {
        m_windfarm_model[0]->set_turb_disk_angle(turb_disk_angle);
    }

protected:

    amrex::Vector<amrex::Real> xloc, yloc;
	amrex::Real my_turb_disk_angle;
    amrex::Real hub_height, rotor_rad, thrust_coeff_standing, nominal_power;
    amrex::Vector<amrex::Real> wind_speed, thrust_coeff, power;

   /*! \brief Create and set the specified windfarm model */
    template<class NewWindFarmModel>
    void SetModel ()
    {
        for (int lev(0); lev<m_windfarm_model.size(); ++lev) {
            m_windfarm_model[lev] = std::make_unique<NewWindFarmModel>();
        }
    }

private:
    amrex::Vector<std::unique_ptr<NullWindFarm>> m_windfarm_model; /*!< windfarm model */
};

#endif


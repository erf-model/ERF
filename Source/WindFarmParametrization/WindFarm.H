#ifndef WINDFARM_H
#define WINDFARM_H

#include <DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>

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
                       bool x_y, bool lat_lon,
                       const amrex::Real latitude_lo=0.0,
                       const amrex::Real longitude_lo=0.0);

    void init_windfarm_lat_lon(const std::string windfarm_loc_table,
                               const amrex::Real latitude_lo,
                               const amrex::Real longitude_lo);

    void init_windfarm_x_y(const std::string windfarm_loc_table);

    void read_windfarm_locations_table(const std::string windfarm_loc_table,
                                          bool x_y, bool lat_lon,
                                          const amrex::Real latitude_lo=0.0,
                                          const amrex::Real longitude_lo=0.0);

    void read_windfarm_spec_table(const std::string windfarm_spec_table);

    void fill_Nturb_multifab(const amrex::Geometry& geom, amrex::MultiFab& mf_Nturb);

    void write_turbine_locations_vtk();

    void write_actuator_disks_vtk();

    void advance (const amrex::Geometry& a_geom,
                  const amrex::Real& dt_advance,
                  amrex::MultiFab& cons_in,
                  amrex::MultiFab& U_old, amrex::MultiFab& V_old, amrex::MultiFab& W_old,
                  amrex::MultiFab& mf_vars_windfarm, const amrex::MultiFab& mf_Nturb) override
    {
        m_windfarm_model[0]->advance(a_geom, dt_advance, cons_in, U_old, V_old, W_old,
                                     mf_vars_windfarm, mf_Nturb);
    }

    void set_turb_spec(const amrex::Real& a_rotor_rad, const amrex::Real& a_hub_height,
                       const amrex::Real& a_thrust_coeff_standing, const amrex::Gpu::DeviceVector<amrex::Real>& a_d_wind_speed,
                       const amrex::Gpu::DeviceVector<amrex::Real>& a_d_thrust_coeff,
                       const amrex::Gpu::DeviceVector<amrex::Real>& a_d_power) override
    {
        m_windfarm_model[0]->set_turb_spec(a_rotor_rad, a_hub_height, a_thrust_coeff_standing,
                  a_d_wind_speed, a_d_thrust_coeff, a_d_power);
    }

    void set_turb_loc (const amrex::Gpu::DeviceVector<amrex::Real>& a_d_xloc,
                       const amrex::Gpu::DeviceVector<amrex::Real>& a_d_yloc) override
    {
        m_windfarm_model[0]->set_turb_loc(a_d_xloc, a_d_yloc);
    }

protected:

    amrex::Vector<amrex::Real> xloc, yloc;
    amrex::Real hub_height, rotor_rad, thrust_coeff_standing, nominal_power;
    amrex::Vector<amrex::Real> wind_speed, thrust_coeff, power;
    amrex::Gpu::DeviceVector<amrex::Real> d_wind_speed, d_thrust_coeff, d_power;
    amrex::Gpu::DeviceVector<amrex::Real> d_xloc, d_yloc;

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


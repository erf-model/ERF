#ifndef WINDFARM_H
#define WINDFARM_H

#include <DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>

class WindFarm {

public:

    WindFarm(){}

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

protected:

    amrex::Vector<amrex::Real> xloc, yloc;
    amrex::Real hub_height, rotor_rad, thrust_coeff_standing, nominal_power;
    amrex::Vector<amrex::Real> wind_speed, thrust_coeff, power;
    amrex::Gpu::DeviceVector<amrex::Real> d_wind_speed, d_thrust_coeff, d_power;
    amrex::Gpu::DeviceVector<amrex::Real> d_xloc, d_yloc;

};


#endif


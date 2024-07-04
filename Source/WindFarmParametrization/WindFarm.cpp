#include <WindFarm.H>
#include <Fitch.H>
#include <EWP.H>
#include <SimpleAD.H>
using namespace amrex;

amrex::Real hub_height, rotor_rad, thrust_coeff_standing, nominal_power;
amrex::Vector<amrex::Real> wind_speed, thrust_coeff, power;

void read_in_table(std::string windfarm_spec_table)
{
    //The first line is the number of pairs entries for the power curve and thrust coefficient.
    //The second line gives first the height in meters of the turbine hub, second, the diameter in
    //meters of the rotor, third the standing thrust coefficient, and fourth the nominal power of
    //the turbine in MW.
    //The remaining lines contain the three values of: wind speed, thrust coefficient, and power production in kW.

     // Read turbine data from wind-turbine-1.tbl
    std::ifstream file_turb_table(windfarm_spec_table);
    if (!file_turb_table.is_open()) {
        amrex::Error("Wind farm specifications table not found. Either the inputs is missing the "
                      "erf.windfarm_spec_table entry or the file specified in the entry - " + windfarm_spec_table + " is missing.");
    }
    else {
        amrex::Print() << "Reading in wind farm specifications table: " << windfarm_spec_table << "\n";
    }

    int nlines;
    file_turb_table >> nlines;
    wind_speed.resize(nlines);
    thrust_coeff.resize(nlines);
    power.resize(nlines);

    Real rotor_dia;
    file_turb_table >> hub_height >> rotor_dia >> thrust_coeff_standing >> nominal_power;
    rotor_rad = rotor_dia*0.5;
    if(rotor_rad > hub_height) {
        amrex::Abort("The blade length is more than the hub height. Check the second line in wind-turbine-1.tbl. Aborting.....");
    }
    if(thrust_coeff_standing > 1.0) {
        amrex::Abort("The standing thrust coefficient is greater than 1. Check the second line in wind-turbine-1.tbl. Aborting.....");
    }

    for(int iline=0;iline<nlines;iline++){
        file_turb_table >> wind_speed[iline] >> thrust_coeff[iline] >> power[iline];
        if(thrust_coeff[iline] > 1.0) {
            amrex::Abort("The thrust coefficient is greater than 1. Check wind-turbine-1.tbl. Aborting.....");
        }
    }
    file_turb_table.close();
}


void advance_windfarm (int lev,
                       const Geometry& geom,
                       const Real& dt_advance,
                       MultiFab& cons_in,
                       MultiFab& U_old, MultiFab& V_old, MultiFab& W_old,
                       MultiFab& mf_vars_windfarm, const amrex::MultiFab& mf_Nturb,
                       SolverChoice& solver_choice)
{
    if (solver_choice.windfarm_type == WindFarmType::Fitch) {
        fitch_advance(lev, geom, dt_advance, cons_in, U_old, V_old, W_old,
                       mf_vars_windfarm, mf_Nturb);
    } else if (solver_choice.windfarm_type == WindFarmType::EWP) {
        ewp_advance(lev, geom, dt_advance, cons_in, U_old, V_old, W_old,
                      mf_vars_windfarm, mf_Nturb);
    } else if (solver_choice.windfarm_type == WindFarmType::SimpleAD) {
        simpleAD_advance(lev, geom, dt_advance, cons_in, U_old, V_old, W_old,
                      mf_vars_windfarm, mf_Nturb);
    }
}

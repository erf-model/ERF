#include <AMReX_Vector.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_ParmParse.H>

#include "ERF.H"

using namespace amrex;

void
ERF::build_geostrophic_forcing()
{
    amrex::ParmParse pp("erf");

    // Read in the geostrophic wind -- we only use this to construct
    //     the forcing term so no need to keep it
    amrex::Vector<amrex::Real> abl_geo_wind(3);
    pp.queryarr("abl_geo_wind",abl_geo_wind);

    // Read the rotational time period (in seconds)
    amrex::Real rot_time_period = 86400.0;
    pp.query("rotational_time_period", rot_time_period);

    Real coriolis_factor = 2.0 * 2.0 * M_PI / rot_time_period;
    amrex::Print() << "Geostrophic forcing: Coriolis factor = "
                   << coriolis_factor << std::endl;

    amrex::Real latitude = 90.0;
    pp.query("latitude", latitude);
    AMREX_ALWAYS_ASSERT(amrex::Math::abs(latitude - 90.0) < 1.e-12);

    abl_geo_forcing =
       {-coriolis_factor * abl_geo_wind[1], coriolis_factor * abl_geo_wind[0], 0.0};
}

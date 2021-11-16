#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>

#include "ERF.H"
#include "Constants.H"

using namespace amrex;

void
ERF::build_coriolis_forcings()
{
    if (!use_coriolis) return;

    amrex::ParmParse pp("erf");

    // Read the rotational time period (in seconds)
    amrex::Real rot_time_period = 86400.0;
    pp.query("rotational_time_period", rot_time_period);

    coriolis_factor = 2.0 * 2.0 * PI / rot_time_period;
    amrex::Print() << "Coriolis factor = " << coriolis_factor << std::endl;

    amrex::Real latitude = 90.0;
    pp.query("latitude", latitude);
    AMREX_ALWAYS_ASSERT(amrex::Math::abs(latitude - 90.0) < 1.e-12);

    // Convert to radians
    latitude *= (PI/180.);
    sinphi = std::sin(latitude);
    cosphi = std::cos(latitude);

    if (abl_driver_type == "GeostrophicWind")
    {
        // Read in the geostrophic wind -- we only use this to construct
        //     the forcing term so no need to keep it
        amrex::Vector<amrex::Real> abl_geo_wind(3);
        pp.queryarr("abl_geo_wind",abl_geo_wind);

        abl_geo_forcing = {
            -coriolis_factor * (abl_geo_wind[1]*sinphi - abl_geo_wind[2]*cosphi),
             coriolis_factor *  abl_geo_wind[0]*sinphi,
            -coriolis_factor *  abl_geo_wind[0]*cosphi
        };
    }
}

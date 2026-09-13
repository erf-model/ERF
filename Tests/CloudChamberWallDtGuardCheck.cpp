#include <AMReX.H>
#include <AMReX_REAL.H>

#include "TimeIntegration/ERF_CloudChamberWallDtGuard.H"

int main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);

    const amrex::Real coefficient = amrex::Real(0.5);
    const amrex::Real tangential_speed = amrex::Real(2.0);
    const amrex::Real dx_inv = amrex::Real(2.0);
    const amrex::Real max_wall_rate =
        coefficient * tangential_speed * dx_inv;

    const double wall_dt = 0.5 / static_cast<double>(max_wall_rate);
    const double fixed_dt = 0.50;

    // Equal-to-limit is allowed because production rejects only fixed_dt > wall_dt.
    if (erf_cloud_chamber_wall_dt_guard::fixed_dt_exceeds_limit(
            wall_dt, wall_dt)) {
        amrex::Finalize();
        return 10;
    }
    if (!erf_cloud_chamber_wall_dt_guard::fixed_dt_exceeds_limit(
            fixed_dt, wall_dt)) {
        amrex::Finalize();
        return 11;
    }

    // Expected to abort before Finalize().
    erf_cloud_chamber_wall_dt_guard::enforce_fixed_dt_limit(
        0, fixed_dt, wall_dt, max_wall_rate);

    amrex::Finalize();
    return 0;
}

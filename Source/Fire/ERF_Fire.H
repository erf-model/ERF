#ifndef ERF_FIRE_H
#define ERF_FIRE_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <ERF_DataStruct.H>
#include <cmath>

/**
 * @class Fire
 * @brief Fire model for 2D wildfire simulations using Rothermel model
 *        with FARSITE-based elliptical expansion
 *
 * This class implements a fire layer for ERF that simulates wildfire
 * propagation using the Rothermel fire spread model combined with
 * FARSITE elliptical fire expansion algorithm.
 */
class Fire {

public:

    /**
     * @brief Constructor
     */
    Fire () { }

    /**
     * @brief Destructor
     */
    virtual ~Fire () = default;

    /**
     * @brief Define the fire model
     * @param lev Level in AMR hierarchy
     * @param sc Solver choice structure
     */
    virtual void
    Define (const int& lev,
            SolverChoice& sc);

    /**
     * @brief Initialize the fire model
     * @param lev Level in AMR hierarchy
     * @param cons_in Conservative variables
     * @param geom Geometry at level
     * @param dt_advance Time step
     */
    virtual void
    Init (const int& lev,
          const amrex::MultiFab& cons_in,
          const amrex::Geometry& geom,
          const amrex::Real& dt_advance);

    /**
     * @brief Advance fire simulation by one time step
     * @param lev Level in AMR hierarchy
     * @param time Current simulation time
     * @param dt_advance Time step size
     * @param cons_in Conservative variables
     * @param geom Geometry at level
     */
    virtual void
    Advance (const int& lev,
             const amrex::Real& time,
             const amrex::Real& dt_advance,
             amrex::MultiFab& cons_in,
             const amrex::Geometry& geom);

    /**
     * @brief Update fire model variables
     * @param lev Level in AMR hierarchy
     * @param cons_in Conservative variables
     */
    virtual void
    Update_Fire_Vars (const int& lev,
                      amrex::MultiFab& cons_in);

    /**
     * @brief Compute fire spread rate using Rothermel model
     * @param lev Level in AMR hierarchy
     * @param geom Geometry at level
     */
    virtual void
    ComputeRothermellSpreadRate (const int& lev,
                                 const amrex::Geometry& geom);

    /**
     * @brief Compute elliptical fire expansion using FARSITE algorithm
     * @param lev Level in AMR hierarchy
     * @param geom Geometry at level
     */
    virtual void
    ComputeEllipticalExpansion (const int& lev,
                                const amrex::Geometry& geom);

    /**
     * @brief Compute fire intensity
     * @param lev Level in AMR hierarchy
     */
    virtual void
    ComputeFireIntensity (const int& lev);

private:

    // Rothermel model parameters
    amrex::Real fuel_moisture_content;      ///< Fuel moisture content (%)
    amrex::Real fuel_bed_depth;             ///< Fuel bed depth (m)
    amrex::Real fuel_particle_density;      ///< Fuel particle density (kg/m^3)
    amrex::Real fuel_energy_content;        ///< Fuel energy content (J/kg)
    amrex::Real fuel_load;                  ///< Fuel load (kg/m^2)
    amrex::Real wind_speed;                 ///< Wind speed (m/s)
    amrex::Real slope;                      ///< Slope angle (degrees)

    // FARSITE elliptical expansion parameters
    amrex::Real ellipse_length_width_ratio; ///< Length-to-width ratio of fire ellipse
    amrex::Real ellipse_eccentricity;       ///< Eccentricity of fire ellipse
    amrex::Real ellipse_major_axis;         ///< Major axis length (m)
    amrex::Real ellipse_minor_axis;         ///< Minor axis length (m)

    // Fire state variables
    amrex::Real head_fire_rate_of_spread;   ///< Head fire rate of spread (m/s)
    amrex::Real flank_fire_rate_of_spread;  ///< Flank fire rate of spread (m/s)
    amrex::Real back_fire_rate_of_spread;   ///< Back fire rate of spread (m/s)
    amrex::Real fire_line_intensity;        ///< Fire line intensity (W/m)
    amrex::Real flame_length;               ///< Flame length (m)

    // Domain variables
    int nx;                                 ///< Number of grid points in x
    int ny;                                 ///< Number of grid points in y

};

#endif

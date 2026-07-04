#include "ERF_Fire.H"

#include <AMReX_ParmParse.H>

/**
 * @brief Define the fire model with default parameters
 */
void
Fire::Define (const int& lev,
              SolverChoice& sc)
{
    // Initialize Rothermel model parameters with defaults
    fuel_moisture_content      = 0.08;  // 8% moisture content
    fuel_bed_depth             = 0.3;   // 0.3 m bed depth
    fuel_particle_density      = 500.0; // 500 kg/m^3
    fuel_energy_content        = 1.9e7; // 19 MJ/kg
    fuel_load                  = 2.0;   // 2 kg/m^2

    // Initialize environmental parameters
    wind_speed                 = 5.0;   // 5 m/s wind
    slope                      = 0.0;   // 0 degree slope

    // Initialize FARSITE elliptical expansion parameters
    ellipse_length_width_ratio = 1.5;   // Default ellipse aspect ratio
    ellipse_eccentricity       = 0.5;   // Default eccentricity
    ellipse_major_axis         = 100.0; // 100 m major axis
    ellipse_minor_axis         = 50.0;  // 50 m minor axis

    // Initialize fire state variables
    head_fire_rate_of_spread   = 0.1;   // 0.1 m/s head fire
    flank_fire_rate_of_spread  = 0.05;  // 0.05 m/s flank fire
    back_fire_rate_of_spread   = 0.02;  // 0.02 m/s back fire
    fire_line_intensity        = 1000.0; // 1000 W/m
    flame_length               = 5.0;    // 5 m flame length

    // Parse input parameters if provided
    amrex::ParmParse pp("fire");
    pp.query("fuel_moisture_content", fuel_moisture_content);
    pp.query("fuel_bed_depth", fuel_bed_depth);
    pp.query("wind_speed", wind_speed);
    pp.query("slope", slope);
    pp.query("ellipse_length_width_ratio", ellipse_length_width_ratio);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Fire model initialized at level " << lev << std::endl;
    }
}

/**
 * @brief Initialize fire variables
 */
void
Fire::Init (const int& lev,
            const amrex::MultiFab& cons_in,
            const amrex::Geometry& geom,
            const amrex::Real& dt_advance)
{
    const auto& dom = geom.Domain();
    nx = dom.length(0);
    ny = dom.length(1);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Fire model initialized with grid: " << nx << " x " << ny << std::endl;
    }
}

/**
 * @brief Advance fire simulation by one time step
 */
void
Fire::Advance (const int& lev,
               const amrex::Real& time,
               const amrex::Real& dt_advance,
               amrex::MultiFab& cons_in,
               const amrex::Geometry& geom)
{
    // Update fire spread rate
    ComputeRothermellSpreadRate(lev, geom);

    // Update elliptical expansion
    ComputeEllipticalExpansion(lev, geom);

    // Compute fire intensity
    ComputeFireIntensity(lev);

    // Update fire variables
    Update_Fire_Vars(lev, cons_in);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Fire model advanced at time " << time
                       << " with timestep " << dt_advance << std::endl;
    }
}

/**
 * @brief Update fire model variables
 */
void
Fire::Update_Fire_Vars (const int& lev,
                        amrex::MultiFab& cons_in)
{
    // Dummy implementation - update fire-related conserved variables
    // This would be filled in with actual fire dynamics in later phases
}

/**
 * @brief Compute fire spread rate using Rothermel model
 */
void
Fire::ComputeRothermellSpreadRate (const int& lev,
                                    const amrex::Geometry& geom)
{
    // Dummy implementation of Rothermel fire spread model
    // Reference: Rothermel, R. C. (1972). A mathematical model for predicting
    // fire spread in wildland fuels. Res. Paper INT-115.
    // USDA Forest Service, Intermountain Forest and Range Experiment Station

    // Update head fire rate of spread based on Rothermel model
    // This is a placeholder that will be replaced with the full model
    head_fire_rate_of_spread = 0.1 * (1.0 + wind_speed * 0.01);
    flank_fire_rate_of_spread = head_fire_rate_of_spread * 0.5;
    back_fire_rate_of_spread = head_fire_rate_of_spread * 0.2;
}

/**
 * @brief Compute elliptical fire expansion using FARSITE algorithm
 */
void
Fire::ComputeEllipticalExpansion (const int& lev,
                                  const amrex::Geometry& geom)
{
    // Dummy implementation of FARSITE elliptical expansion
    // Reference: Finney, M. A. (2004). FARSITE: Fire Area Simulator–model
    // development and evaluation. Res. Paper RMRS-RP-4 Revised.
    // USDA Forest Service, Rocky Mountain Research Station

    // Compute ellipse dimensions based on wind and slope effects
    ellipse_major_axis = 100.0 * (1.0 + wind_speed * 0.02);
    ellipse_minor_axis = ellipse_major_axis / ellipse_length_width_ratio;
    ellipse_eccentricity = std::sqrt(1.0 - (ellipse_minor_axis * ellipse_minor_axis) /
                                     (ellipse_major_axis * ellipse_major_axis));
}

/**
 * @brief Compute fire intensity
 */
void
Fire::ComputeFireIntensity (const int& lev)
{
    // Dummy implementation of fire intensity calculation
    // Fire line intensity: I = r * H
    // where r = rate of fire spread (m/s), H = fuel energy content (J/kg)
    fire_line_intensity = head_fire_rate_of_spread * fuel_energy_content * fuel_load;

    // Flame length: L = 0.0775 * I^0.46
    // Using Thomas (1963) formula
    flame_length = 0.0775 * std::pow(fire_line_intensity, 0.46);
}

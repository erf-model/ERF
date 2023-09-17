#ifndef _SPONGE_STRUCT_H_
#define _SPONGE_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

/**
 * Container holding the sponge-related choices
 */

struct SpongeChoice {
  public:
    void init_params()
    {
        amrex::ParmParse pp(pp_prefix);

        // Include sponge damping
        pp.query("use_xlo_sponge_damping", use_xlo_sponge_damping);
        pp.query("use_xhi_sponge_damping", use_xhi_sponge_damping);
        pp.query("use_ylo_sponge_damping", use_ylo_sponge_damping);
        pp.query("use_yhi_sponge_damping", use_yhi_sponge_damping);
        pp.query("use_zlo_sponge_damping", use_zlo_sponge_damping);
        pp.query("use_zhi_sponge_damping", use_zhi_sponge_damping);

        pp.query("sponge_strength" , sponge_strength);
        pp.query("xlo_sponge_end"  , xlo_sponge_end);
        pp.query("xhi_sponge_start", xhi_sponge_start);
        pp.query("ylo_sponge_end"  , ylo_sponge_end);
        pp.query("yhi_sponge_start", yhi_sponge_start);
        pp.query("zlo_sponge_end"  , zlo_sponge_end);
        pp.query("zhi_sponge_start", zhi_sponge_start);

        pp.query("sponge_density" , sponge_density);
        pp.query("sponge_x_velocity" , sponge_x_velocity);
        pp.query("sponge_y_velocity" , sponge_y_velocity);
        pp.query("sponge_z_velocity" , sponge_z_velocity);
    }

    void display()
    {
        amrex::Print() << "Sponge choices: " << std::endl;
    }

    // Default prefix
    std::string pp_prefix {"erf"};

    bool    use_xlo_sponge_damping     = false;
    bool    use_xhi_sponge_damping     = false;
    bool    use_ylo_sponge_damping     = false;
    bool    use_yhi_sponge_damping     = false;
    bool    use_zlo_sponge_damping     = false;
    bool    use_zhi_sponge_damping     = false;

    amrex::Real sponge_strength;

    amrex::Real xlo_sponge_end, xhi_sponge_start;
    amrex::Real ylo_sponge_end, yhi_sponge_start;
    amrex::Real zlo_sponge_end, zhi_sponge_start;
    amrex::Real sponge_density, sponge_x_velocity, sponge_y_velocity, sponge_z_velocity;
};
#endif

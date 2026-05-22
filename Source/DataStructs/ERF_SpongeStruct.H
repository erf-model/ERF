#ifndef ERF_SPONGE_STRUCT_H_
#define ERF_SPONGE_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

AMREX_ENUM(SpongeType,
    None, Input_Sponge
);

/**
 * Container holding the sponge-related choices
 */

struct SpongeChoice {
  public:
    void init_params(std::string pp_prefix)
    {
        amrex::ParmParse pp(pp_prefix);

        pp.query_enum_case_insensitive("sponge_type",sponge_type);

        // Include sponge damping
        if (sponge_type != SpongeType::None) {
            pp.query("use_xlo_sponge_damping", use_xlo_sponge_damping);
            pp.query("use_xhi_sponge_damping", use_xhi_sponge_damping);
            pp.query("use_ylo_sponge_damping", use_ylo_sponge_damping);
            pp.query("use_yhi_sponge_damping", use_yhi_sponge_damping);
            pp.query("use_zlo_sponge_damping", use_zlo_sponge_damping);
            pp.query("use_zhi_sponge_damping", use_zhi_sponge_damping);

            pp.get("sponge_strength" , sponge_strength);

            if (use_xlo_sponge_damping) {
                pp.get("xlo_sponge_end"  , xlo_sponge_end);
            }
            if (use_xhi_sponge_damping) {
                pp.get("xhi_sponge_start"  , xhi_sponge_start);
            }
            if (use_ylo_sponge_damping) {
                pp.get("ylo_sponge_end"  , ylo_sponge_end);
            }
            if (use_yhi_sponge_damping) {
                pp.get("yhi_sponge_start"  , yhi_sponge_start);
            }
            if (use_zlo_sponge_damping) {
                pp.get("zlo_sponge_end"  , zlo_sponge_end);
            }
            if (use_zhi_sponge_damping) {
                pp.get("zhi_sponge_start"  , zhi_sponge_start);
            }

            pp.query("sponge_density"    , sponge_density);
            pp.query("sponge_rhotheta"   , sponge_rhotheta);
            pp.query("sponge_rhomoist"   , sponge_rhomoist);
            pp.query("sponge_x_velocity" , sponge_x_velocity);
            pp.query("sponge_y_velocity" , sponge_y_velocity);
            pp.query("sponge_z_velocity" , sponge_z_velocity);
        }
    }

    void display()
    {
        amrex::Print() << "Sponge choices: " << std::endl;
    }

    bool    use_xlo_sponge_damping     = false;
    bool    use_xhi_sponge_damping     = false;
    bool    use_ylo_sponge_damping     = false;
    bool    use_yhi_sponge_damping     = false;
    bool    use_zlo_sponge_damping     = false;
    bool    use_zhi_sponge_damping     = false;

    amrex::Real sponge_strength = 0.0;

    inline static
    SpongeType sponge_type = SpongeType::None;

    amrex::Real xlo_sponge_end, xhi_sponge_start;
    amrex::Real ylo_sponge_end, yhi_sponge_start;
    amrex::Real zlo_sponge_end, zhi_sponge_start;
    amrex::Real sponge_density = -one;
    amrex::Real sponge_rhotheta = -one;
    amrex::Real sponge_rhomoist = -one;
    amrex::Real sponge_x_velocity, sponge_y_velocity, sponge_z_velocity;
};
#endif

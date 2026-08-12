#ifndef ERF_SPONGE_STRUCT_H_
#define ERF_SPONGE_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

/**
 * @brief Sponge damping model used near selected domain boundaries.
 */
AMREX_ENUM(SpongeType,
    None, Standard, Input_Sponge
);

/**
 * Container holding the sponge-related choices
 */

struct SpongeChoice {
  public:
    /**
     * @brief Read sponge-layer options from the input parameter database.
     * @param pp_prefix ParmParse prefix for the ERF input namespace.
     */
    void init_params(std::string pp_prefix)
    {
        amrex::ParmParse pp(pp_prefix);

        // The default type is SpongeType::Standard
        // Here we can overwrite the default with None or Input_Sponge
        pp.query_enum_case_insensitive("sponge_type",sponge_type);

        if (sponge_type == SpongeType::None) return;

        pp.query("sponge_strength" , sponge_strength);

        // If sponge_strength == 0 then we set sponge_type to None
        if (std::abs(sponge_strength) == zero)
        {
            sponge_type = SpongeType::None;
        }
        else
        {
            pp.query("use_xlo_sponge_damping", use_xlo_sponge_damping);
            pp.query("use_xhi_sponge_damping", use_xhi_sponge_damping);
            pp.query("use_ylo_sponge_damping", use_ylo_sponge_damping);
            pp.query("use_yhi_sponge_damping", use_yhi_sponge_damping);
            pp.query("use_zlo_sponge_damping", use_zlo_sponge_damping);
            pp.query("use_zhi_sponge_damping", use_zhi_sponge_damping);

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

    /**
     * @brief Print the configured sponge-layer options.
     */
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
    SpongeType sponge_type = SpongeType::Standard;

    amrex::Real xlo_sponge_end, xhi_sponge_start;
    amrex::Real ylo_sponge_end, yhi_sponge_start;
    amrex::Real zlo_sponge_end, zhi_sponge_start;
    amrex::Real sponge_density = -one;
    amrex::Real sponge_rhotheta = -one;
    amrex::Real sponge_rhomoist = -one;
    amrex::Real sponge_x_velocity = zero;
    amrex::Real sponge_y_velocity = zero;
    amrex::Real sponge_z_velocity = zero;
};

/**
 * @var SpongeChoice::use_xlo_sponge_damping
 * @brief Whether sponge damping is applied at the low-x boundary.
 * @var SpongeChoice::use_xhi_sponge_damping
 * @brief Whether sponge damping is applied at the high-x boundary.
 * @var SpongeChoice::use_ylo_sponge_damping
 * @brief Whether sponge damping is applied at the low-y boundary.
 * @var SpongeChoice::use_yhi_sponge_damping
 * @brief Whether sponge damping is applied at the high-y boundary.
 * @var SpongeChoice::use_zlo_sponge_damping
 * @brief Whether sponge damping is applied at the low-z boundary.
 * @var SpongeChoice::use_zhi_sponge_damping
 * @brief Whether sponge damping is applied at the high-z boundary.
 * @var SpongeChoice::sponge_strength
 * @brief Sponge damping strength.
 * @var SpongeChoice::sponge_type
 * @brief Selected sponge damping model.
 * @var SpongeChoice::xlo_sponge_end
 * @brief Physical end location of the low-x sponge layer.
 * @var SpongeChoice::xhi_sponge_start
 * @brief Physical start location of the high-x sponge layer.
 * @var SpongeChoice::ylo_sponge_end
 * @brief Physical end location of the low-y sponge layer.
 * @var SpongeChoice::yhi_sponge_start
 * @brief Physical start location of the high-y sponge layer.
 * @var SpongeChoice::zlo_sponge_end
 * @brief Physical end location of the low-z sponge layer.
 * @var SpongeChoice::zhi_sponge_start
 * @brief Physical start location of the high-z sponge layer.
 * @var SpongeChoice::sponge_density
 * @brief Target density used by sponge damping.
 * @var SpongeChoice::sponge_rhotheta
 * @brief Target rho-theta used by sponge damping.
 * @var SpongeChoice::sponge_rhomoist
 * @brief Target moist density variable used by sponge damping.
 * @var SpongeChoice::sponge_x_velocity
 * @brief Target x velocity used by sponge damping.
 * @var SpongeChoice::sponge_y_velocity
 * @brief Target y velocity used by sponge damping.
 * @var SpongeChoice::sponge_z_velocity
 * @brief Target vertical velocity used by sponge damping.
 */
#endif

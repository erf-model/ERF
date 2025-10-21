#ifndef ERF_ADV_STRUCT_H_
#define ERF_ADV_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>

#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>

/**
 * Container holding the advection-related choices
 */

struct AdvChoice {
  public:
    void init_params(std::string pp_prefix)
    {
        amrex::ParmParse pp(pp_prefix);

        // Order and type of spatial discretizations used in advection
        pp.query("use_efficient_advection", use_efficient_advection);
        std::string dycore_horiz_adv_string    = "" ; std::string dycore_vert_adv_string   = "";
        std::string dryscal_horiz_adv_string   = "" ; std::string dryscal_vert_adv_string  = "";
        pp.query("dycore_horiz_adv_type"   , dycore_horiz_adv_string);
        pp.query("dycore_vert_adv_type"    , dycore_vert_adv_string);
        pp.query("dryscal_horiz_adv_type"  , dryscal_horiz_adv_string);
        pp.query("dryscal_vert_adv_type"   , dryscal_vert_adv_string);

        std::string moistscal_horiz_adv_string = ""; std::string moistscal_vert_adv_string = "";
        pp.query("moistscal_horiz_adv_type", moistscal_horiz_adv_string);
        pp.query("moistscal_vert_adv_type" , moistscal_vert_adv_string);

        if (use_efficient_advection){
           amrex::Print() << "Using efficient advection scheme" << std::endl;
        }

        if ( (dycore_horiz_adv_string == "Blended_3rd4th") ||
             (dycore_horiz_adv_string == "Blended_5th6th") )
        {
            pp.query("dycore_horiz_upw_frac" , dycore_horiz_upw_frac);
            AMREX_ASSERT_WITH_MESSAGE((dycore_horiz_upw_frac >= 0.) && (dycore_horiz_upw_frac <= 1.),
                                      "The dycore horizontal upwinding fraction must be between 0 and 1");
        }

        if ( (dycore_vert_adv_string == "Blended_3rd4th") ||
             (dycore_vert_adv_string == "Blended_5th6th") )
        {
            pp.query("dycore_vert_upw_frac" , dycore_vert_upw_frac);
            AMREX_ASSERT_WITH_MESSAGE((dycore_vert_upw_frac >= 0.) && (dycore_vert_upw_frac <= 1.),
                                      "The dycore vertical upwinding fraction must be between 0 and 1");
        }

        if ( (dryscal_horiz_adv_string == "Blended_3rd4th") ||
             (dryscal_horiz_adv_string == "Blended_5th6th") )
        {
            pp.query("dryscal_horiz_upw_frac" , dryscal_horiz_upw_frac);
            AMREX_ASSERT_WITH_MESSAGE((dryscal_horiz_upw_frac >= 0.) && (dryscal_horiz_upw_frac <= 1.),
                                      "The dry scalar horizontal upwinding fraction must be between 0 and 1");
        }

        if ( (dryscal_vert_adv_string == "Blended_3rd4th") ||
             (dryscal_vert_adv_string == "Blended_5th6th") )
        {
            pp.query("dryscal_vert_upw_frac" , dryscal_vert_upw_frac);
            AMREX_ASSERT_WITH_MESSAGE((dryscal_vert_upw_frac >= 0.) && (dryscal_vert_upw_frac <= 1.),
                                      "The dry scalar vertical upwinding fraction must be between 0 and 1");
        }

        if ( (moistscal_horiz_adv_string == "Blended_3rd4th") ||
             (moistscal_horiz_adv_string == "Blended_5th6th") )
        {
            pp.query("moistscal_horiz_upw_frac" , moistscal_horiz_upw_frac);
            AMREX_ASSERT_WITH_MESSAGE((moistscal_horiz_upw_frac >= 0.) && (moistscal_horiz_upw_frac <= 1.),
                                      "The moist scalar horizontal upwinding fraction must be between 0 and 1");
        }

        if ( (moistscal_vert_adv_string == "Blended_3rd4th") ||
             (moistscal_vert_adv_string == "Blended_5th6th") )
        {
            pp.query("moistscal_vert_upw_frac" , moistscal_vert_upw_frac);
            AMREX_ASSERT_WITH_MESSAGE((moistscal_vert_upw_frac >= 0.) && (moistscal_vert_upw_frac <= 1.),
                                      "The moist scalar vertical upwinding fraction must be between 0 and 1");
        }

        if (dycore_horiz_adv_string != "") {
            if ( (dycore_horiz_adv_string == "Centered_2nd") ||
                 (dycore_horiz_adv_string == "Upwind_3rd"  ) ||
                 (dycore_horiz_adv_string == "Blended_3rd4th") ||
                 (dycore_horiz_adv_string == "Centered_4th") ||
                 (dycore_horiz_adv_string == "Upwind_5th"  ) ||
                (dycore_horiz_adv_string == "Blended_5th6th") ||
                (dycore_horiz_adv_string == "Centered_6th") ||
                (dycore_horiz_adv_string == "WENO3"       ) ||
                (dycore_horiz_adv_string == "WENOZ3"      ) ||
                (dycore_horiz_adv_string == "WENO5"       ) ||
                (dycore_horiz_adv_string == "WENOZ5"      ) ||
                (dycore_horiz_adv_string == "WENO7"       ) ||
                (dycore_horiz_adv_string == "WENOZ7"      ))
           {
               dycore_horiz_adv_type = adv_type_convert_string_to_advtype(dycore_horiz_adv_string);
           } else {
               amrex::Error("Don't know this dycore_horiz_adv_string");
           }
        }

        if (dycore_vert_adv_string != "") {
            if ( (dycore_vert_adv_string == "Centered_2nd") ||
                 (dycore_vert_adv_string == "Upwind_3rd"  ) ||
                 (dycore_vert_adv_string == "Blended_3rd4th") ||
                 (dycore_vert_adv_string == "Centered_4th") ||
                 (dycore_vert_adv_string == "Upwind_5th"  ) ||
                 (dycore_vert_adv_string == "Blended_5th6th") ||
                 (dycore_vert_adv_string == "Centered_6th") ||
                 (dycore_vert_adv_string == "WENO3"       ) ||
                 (dycore_vert_adv_string == "WENOZ3"      ) ||
                 (dycore_vert_adv_string == "WENO5"       ) ||
                 (dycore_vert_adv_string == "WENOZ5"      ) ||
                 (dycore_vert_adv_string == "WENO7"       ) ||
                 (dycore_vert_adv_string == "WENOZ7"      ))
            {
                dycore_vert_adv_type = adv_type_convert_string_to_advtype(dycore_vert_adv_string);
            } else {
                amrex::Error("Don't know this dycore_vert_adv_string");
            }
        }

        if (dryscal_horiz_adv_string != "") {
            if ( (dryscal_horiz_adv_string == "Centered_2nd") ||
                 (dryscal_horiz_adv_string == "Upwind_3rd"  ) ||
                 (dryscal_horiz_adv_string == "Upwind_3rd_SL") ||
                 (dryscal_horiz_adv_string == "Blended_3rd4th") ||
                 (dryscal_horiz_adv_string == "Centered_4th") ||
                 (dryscal_horiz_adv_string == "Upwind_5th"  ) ||
                 (dryscal_horiz_adv_string == "Blended_5th6th") ||
                 (dryscal_horiz_adv_string == "Centered_6th") ||
                 (dryscal_horiz_adv_string == "WENO3"       ) ||
                 (dryscal_horiz_adv_string == "WENOZ3"      ) ||
                 (dryscal_horiz_adv_string == "WENOMZQ3"    ) ||
                 (dryscal_horiz_adv_string == "WENO5"       ) ||
                 (dryscal_horiz_adv_string == "WENOZ5"      ) ||
                 (dryscal_horiz_adv_string == "WENO7"       ) ||
                 (dryscal_horiz_adv_string == "WENOZ7"      ) )
            {
                dryscal_horiz_adv_type = adv_type_convert_string_to_advtype(dryscal_horiz_adv_string);
            } else {
                amrex::Error("Don't know this dryscal_horiz_adv_string");
            }
        }

        if (dryscal_vert_adv_string != "") {
            if ( (dryscal_vert_adv_string == "Centered_2nd") ||
                 (dryscal_vert_adv_string == "Upwind_3rd"  ) ||
                 (dryscal_vert_adv_string == "Upwind_3rd_SL") ||
                 (dryscal_vert_adv_string == "Blended_3rd4th") ||
                 (dryscal_vert_adv_string == "Centered_4th") ||
                 (dryscal_vert_adv_string == "Upwind_5th"  ) ||
                 (dryscal_vert_adv_string == "Blended_5th6th") ||
                 (dryscal_vert_adv_string == "Centered_6th") ||
                 (dryscal_vert_adv_string == "WENO3"       ) ||
                 (dryscal_vert_adv_string == "WENOZ3"      ) ||
                 (dryscal_vert_adv_string == "WENOMZQ3"    ) ||
                 (dryscal_vert_adv_string == "WENO5"       ) ||
                 (dryscal_vert_adv_string == "WENOZ5"      ) ||
                 (dryscal_vert_adv_string == "WENO7"       ) ||
                 (dryscal_vert_adv_string == "WENOZ7"      ))
            {
                dryscal_vert_adv_type = adv_type_convert_string_to_advtype(dryscal_vert_adv_string);
            } else {
                amrex::Error("Don't know this dryscal_vert_adv_string");
            }
        }

        if (moistscal_horiz_adv_string != "") {
            if ( (moistscal_horiz_adv_string == "Centered_2nd") ||
                 (moistscal_horiz_adv_string == "Upwind_3rd"  ) ||
                 (moistscal_horiz_adv_string == "Upwind_3rd_SL"  ) ||
                 (moistscal_horiz_adv_string == "Blended_3rd4th") ||
                 (moistscal_horiz_adv_string == "Centered_4th") ||
                 (moistscal_horiz_adv_string == "Upwind_5th"  ) ||
                 (moistscal_horiz_adv_string == "Blended_5th6th") ||
                 (moistscal_horiz_adv_string == "Centered_6th") ||
                 (moistscal_horiz_adv_string == "WENO3"       ) ||
                 (moistscal_horiz_adv_string == "WENOZ3"      ) ||
                 (moistscal_horiz_adv_string == "WENOMZQ3"    ) ||
                 (moistscal_horiz_adv_string == "WENO5"       ) ||
                 (moistscal_horiz_adv_string == "WENOZ5"      ) ||
                 (moistscal_horiz_adv_string == "WENO7"       ) ||
                 (moistscal_horiz_adv_string == "WENOZ7"      ))
            {
                moistscal_horiz_adv_type = adv_type_convert_string_to_advtype(moistscal_horiz_adv_string);
            } else {
                amrex::Error("Don't know this moistscal_horiz_adv_string");
            }
        }

        if (moistscal_vert_adv_string != "") {
            if ( (moistscal_vert_adv_string == "Centered_2nd") ||
                 (moistscal_vert_adv_string == "Upwind_3rd"  ) ||
                 (moistscal_vert_adv_string == "Upwind_3rd_SL"  ) ||
                 (moistscal_vert_adv_string == "Blended_3rd4th") ||
                 (moistscal_vert_adv_string == "Centered_4th") ||
                 (moistscal_vert_adv_string == "Upwind_5th"  ) ||
                 (moistscal_vert_adv_string == "Blended_5th6th") ||
                 (moistscal_vert_adv_string == "Centered_6th") ||
                 (moistscal_vert_adv_string == "WENO3"       ) ||
                 (moistscal_vert_adv_string == "WENOZ3"      ) ||
                 (moistscal_vert_adv_string == "WENOMZQ3"    ) ||
                 (moistscal_vert_adv_string == "WENO5"       ) ||
                 (moistscal_vert_adv_string == "WENOZ5"      ) ||
                 (moistscal_vert_adv_string == "WENO7"       ) ||
                 (moistscal_vert_adv_string == "WENOZ7"      ))
            {
                moistscal_vert_adv_type = adv_type_convert_string_to_advtype(moistscal_vert_adv_string);
            } else {
                amrex::Error("Don't know this moistscal_vert_adv_string");
            }
        }

        // We can can't use weno with Embedded Boundaries (for now)
        std::string terrain_type;
        pp.query("terrain_type", terrain_type);
        validate_weno_momentum_compatibility(terrain_type == "EB");

        pp.queryarr("zero_xflux_faces", zero_xflux);
        pp.queryarr("zero_yflux_faces", zero_yflux);
        pp.queryarr("zero_zflux_faces", zero_zflux);
        have_zero_flux_faces = ((zero_xflux.size() > 0) || (zero_yflux.size() > 0) || (zero_zflux.size() > 0));
    }

    void display(std::string& pp_prefix)
    {
        amrex::Print() << "Advection Choices: " << std::endl;

        // We read these again here just as a convenience to avoid passing them
        amrex::ParmParse pp(pp_prefix);
        std::string dycore_horiz_adv_string    = "" ; std::string dycore_vert_adv_string   = "";
        std::string dryscal_horiz_adv_string   = "" ; std::string dryscal_vert_adv_string  = "";
        pp.query("dycore_horiz_adv_type"   , dycore_horiz_adv_string);
        pp.query("dycore_vert_adv_type"    , dycore_vert_adv_string);
        pp.query("dryscal_horiz_adv_type"  , dryscal_horiz_adv_string);
        pp.query("dryscal_vert_adv_type"   , dryscal_vert_adv_string);
        std::string moistscal_horiz_adv_string = ""; std::string moistscal_vert_adv_string = "";
        pp.query("moistscal_horiz_adv_type", moistscal_horiz_adv_string);
        pp.query("moistscal_vert_adv_type" , moistscal_vert_adv_string);

        if (dycore_horiz_adv_string != "") {
            amrex::Print() << "    dycore_horiz_adv_type   : " << dycore_horiz_adv_string;
        } else {
            amrex::Print() << "    dycore_horiz_adv_type   : " << adv_type_convert_int_to_string(dycore_horiz_adv_type);
        }
        if ( (dycore_horiz_adv_string == "Blended_3rd4th") ||
             (dycore_horiz_adv_string == "Blended_5th6th") ) {
            amrex::Print() << " with " << 100*dycore_horiz_upw_frac << "% upwinding";
        }
        amrex::Print() << std::endl;

        if (dycore_vert_adv_string != "") {
            amrex::Print() << "    dycore_vert_adv_type    : " << dycore_vert_adv_string;
        } else {
            amrex::Print() << "    dycore_vert_adv_type    : " << adv_type_convert_int_to_string(dycore_vert_adv_type);
        }
        if ( (dycore_vert_adv_string == "Blended_3rd4th") ||
             (dycore_vert_adv_string == "Blended_5th6th") ) {
            amrex::Print() << " with " << 100*dycore_vert_upw_frac << "% upwinding";
        }
        amrex::Print() << std::endl;

        if (dryscal_horiz_adv_string != "") {
            amrex::Print() << "   dryscal_horiz_adv_type   : " << dryscal_horiz_adv_string;
        } else {
            amrex::Print() << "   dryscal_horiz_adv_type   : " << adv_type_convert_int_to_string(dryscal_horiz_adv_type);
        }
        if ( (dryscal_horiz_adv_string == "Blended_3rd4th") ||
             (dryscal_horiz_adv_string == "Blended_5th6th") ) {
            amrex::Print() << " with " << 100*dryscal_horiz_upw_frac << "% upwinding";
        }
        amrex::Print() << std::endl;

        if (dryscal_vert_adv_string != "") {
            amrex::Print() << "   dryscal_vert_adv_type    : " << dryscal_vert_adv_string;
        } else {
            amrex::Print() << "   dryscal_vert_adv_type    : " << adv_type_convert_int_to_string(dryscal_vert_adv_type);
        }
        if ( (dryscal_vert_adv_string == "Blended_3rd4th") ||
             (dryscal_vert_adv_string == "Blended_5th6th") ) {
            amrex::Print() << " with " << 100*dryscal_vert_upw_frac << "% upwinding";
        }
        amrex::Print() << std::endl;

        if (moistscal_horiz_adv_string != "") {
            amrex::Print() << " moistscal_horiz_adv_type   : " << moistscal_horiz_adv_string;
        } else {
            amrex::Print() << " moistscal_horiz_adv_type   : " << adv_type_convert_int_to_string(moistscal_horiz_adv_type);
        }
        if ( (moistscal_horiz_adv_string == "Blended_3rd4th") ||
             (moistscal_horiz_adv_string == "Blended_5th6th") ) {
            amrex::Print() << " with " << 100*moistscal_horiz_upw_frac << "% upwinding";
        }
        amrex::Print() << std::endl;

        if (moistscal_vert_adv_string != "") {
            amrex::Print() << " moistscal_vert_adv_type    : " << moistscal_vert_adv_string;
        } else {
            amrex::Print() << " moistscal_vert_adv_type    : " << adv_type_convert_int_to_string(moistscal_vert_adv_type);
        }
        if ( (moistscal_vert_adv_string == "Blended_3rd4th") ||
             (moistscal_vert_adv_string == "Blended_5th6th") ) {
            amrex::Print() << " with " << 100*moistscal_vert_upw_frac << "% upwinding";
        }
        amrex::Print() << std::endl;
    }

    std::string
    adv_type_convert_int_to_string (AdvType adv_int)
    {
        if (adv_int == AdvType::Centered_2nd) {
            return "Centered_2nd";
        } else if (adv_int == AdvType::Upwind_3rd) {
            return "Upwind_3rd";
        } else if (adv_int == AdvType::Upwind_3rd_SL) {
            return "Upwind_3rd_SL";
        } else if (adv_int == AdvType::Centered_4th) {
            return "Centered_4th";
        } else if (adv_int == AdvType::Upwind_5th) {
            return "Upwind_5th";
        } else if (adv_int == AdvType::Centered_6th) {
            return "Centered_6th";
        } else if (adv_int == AdvType::Weno_3) {
            return "WENO3";
        } else if (adv_int == AdvType::Weno_3Z) {
            return "WENOZ3";
        } else if (adv_int == AdvType::Weno_5) {
            return "WENO5";
        } else if (adv_int == AdvType::Weno_5Z) {
            return "WENOZ5";
        } else if (adv_int == AdvType::Weno_3MZQ) {
            return "WENOMZQ3";
        } else if (adv_int == AdvType::Weno_7) {
            return "WENO7";
        } else if (adv_int == AdvType::Weno_7Z) {
            return "WENOZ7";
        } else {
            return "Unknown";
        }
    }

    AdvType adv_type_convert_string_to_advtype (std::string adv_string)
    {
        if (adv_string == "Centered_2nd") {
            return AdvType::Centered_2nd;
        } else if ((adv_string == "Upwind_3rd") || (adv_string == "Blended_3rd4th")) {
            return AdvType::Upwind_3rd;
        } else if (adv_string == "Upwind_3rd_SL") {
            return AdvType::Upwind_3rd_SL;
        } else if (adv_string == "Centered_4th") {
            return AdvType::Centered_4th;
        } else if (adv_string == "Upwind_5th" || (adv_string == "Blended_5th6th")) {
            return AdvType::Upwind_5th;
        } else if (adv_string == "Centered_6th") {
            return AdvType::Centered_6th;
        } else if (adv_string == "WENO3") {
            return AdvType::Weno_3;
        } else if (adv_string == "WENOZ3") {
            return AdvType::Weno_3Z;
        } else if (adv_string == "WENO5") {
            return AdvType::Weno_5;
        } else if (adv_string == "WENOZ5") {
            return AdvType::Weno_5Z;
        } else if (adv_string == "WENOMZQ3") {
            return AdvType::Weno_3MZQ;
        } else if (adv_string == "WENO7") {
            return AdvType::Weno_7;
        } else if (adv_string == "WENOZ7") {
            return AdvType::Weno_7Z;
        } else {
            return AdvType::Unknown;
        }
    }

    void validate_weno_momentum_compatibility(bool use_embedded_boundary)
    {
        // Check if any WENO momentum advection is being used
        bool using_weno_momentum =
            (dycore_horiz_adv_type == AdvType::Weno_3) ||
            (dycore_horiz_adv_type == AdvType::Weno_3Z) ||
            (dycore_horiz_adv_type == AdvType::Weno_3MZQ) ||
            (dycore_horiz_adv_type == AdvType::Weno_5) ||
            (dycore_horiz_adv_type == AdvType::Weno_5Z) ||
            (dycore_horiz_adv_type == AdvType::Weno_7) ||
            (dycore_horiz_adv_type == AdvType::Weno_7Z) ||
            (dycore_vert_adv_type == AdvType::Weno_3) ||
            (dycore_vert_adv_type == AdvType::Weno_3Z) ||
            (dycore_vert_adv_type == AdvType::Weno_3MZQ) ||
            (dycore_vert_adv_type == AdvType::Weno_5) ||
            (dycore_vert_adv_type == AdvType::Weno_5Z) ||
            (dycore_vert_adv_type == AdvType::Weno_7) ||
            (dycore_vert_adv_type == AdvType::Weno_7Z);

        // Check embedded boundary incompatibility
        amrex::ignore_unused(use_embedded_boundary, using_weno_momentum);
        AMREX_ASSERT_WITH_MESSAGE(!using_weno_momentum || !use_embedded_boundary,
            "WENO momentum advection schemes are not yet implemented for embedded boundaries. "
            "Please use traditional momentum advection schemes when using embedded boundaries:\n"
            "  - Centered_2nd\n"
            "  - Upwind_3rd\n"
            "  - Centered_4th\n"
            "  - Upwind_5th\n"
            "  - Centered_6th");
    }

    // Order and type of spatial discretizations used in advection
    // Defaults given below but these can be over-written at run-time
    bool use_efficient_advection = false;
    AdvType dycore_horiz_adv_type    = AdvType::Upwind_3rd;
    AdvType dycore_vert_adv_type     = AdvType::Upwind_3rd;
    AdvType dryscal_horiz_adv_type   = AdvType::Upwind_3rd;
    AdvType dryscal_vert_adv_type    = AdvType::Upwind_3rd;
    AdvType moistscal_horiz_adv_type = AdvType::Weno_3;
    AdvType moistscal_vert_adv_type  = AdvType::Weno_3;

    // Blending between upwind and the next-highest order central scheme
    // Note: This is used by both Upwind_* and Blended_* schemes -- the default
    //       is a pure upwind scheme
    amrex::Real dycore_horiz_upw_frac    = 1.0;
    amrex::Real dycore_vert_upw_frac     = 1.0;
    amrex::Real dryscal_horiz_upw_frac   = 1.0;
    amrex::Real dryscal_vert_upw_frac    = 1.0;
    amrex::Real moistscal_horiz_upw_frac = 1.0;
    amrex::Real moistscal_vert_upw_frac  = 1.0;

    // Thin immersed bodies
    amrex::Vector<amrex::IntVect> zero_xflux;
    amrex::Vector<amrex::IntVect> zero_yflux;
    amrex::Vector<amrex::IntVect> zero_zflux;
    bool have_zero_flux_faces = false;
};
#endif

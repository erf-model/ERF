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
    void init_params()
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

        if ( (dycore_horiz_adv_string == "Centered_2nd") ||
             (dycore_horiz_adv_string == "Upwind_3rd"  ) ||
             (dycore_horiz_adv_string == "Blended_3rd4th") ||
             (dycore_horiz_adv_string == "Centered_4th") ||
             (dycore_horiz_adv_string == "Upwind_5th"  ) ||
             (dycore_horiz_adv_string == "Blended_5th6th") ||
             (dycore_horiz_adv_string == "Centered_6th") )
        {
            dycore_horiz_adv_type = adv_type_convert_string_to_advtype(dycore_horiz_adv_string);
        }

        if ( (dycore_vert_adv_string == "Centered_2nd") ||
             (dycore_vert_adv_string == "Upwind_3rd"  ) ||
             (dycore_vert_adv_string == "Blended_3rd4th") ||
             (dycore_vert_adv_string == "Centered_4th") ||
             (dycore_vert_adv_string == "Upwind_5th"  ) ||
             (dycore_vert_adv_string == "Blended_5th6th") ||
             (dycore_vert_adv_string == "Centered_6th") )
        {
            dycore_vert_adv_type = adv_type_convert_string_to_advtype(dycore_vert_adv_string);
        }

        if ( (dryscal_horiz_adv_string == "Centered_2nd") ||
             (dryscal_horiz_adv_string == "Upwind_3rd"  ) ||
             (dryscal_horiz_adv_string == "Blended_3rd4th") ||
             (dryscal_horiz_adv_string == "Centered_4th") ||
             (dryscal_horiz_adv_string == "Upwind_5th"  ) ||
             (dryscal_horiz_adv_string == "Blended_5th6th") ||
             (dryscal_horiz_adv_string == "Centered_6th") ||
             (dryscal_horiz_adv_string == "WENO3"       ) ||
             (dryscal_horiz_adv_string == "WENOZ3"      ) ||
             (dryscal_horiz_adv_string == "WENOMZQ3"    ) ||
             (dryscal_horiz_adv_string == "WENO5"       ) ||
             (dryscal_horiz_adv_string == "WENOZ5"      ) )
        {
            dryscal_horiz_adv_type = adv_type_convert_string_to_advtype(dryscal_horiz_adv_string);
        }

        if ( (dryscal_vert_adv_string == "Centered_2nd") ||
             (dryscal_vert_adv_string == "Upwind_3rd"  ) ||
             (dryscal_vert_adv_string == "Blended_3rd4th") ||
             (dryscal_vert_adv_string == "Centered_4th") ||
             (dryscal_vert_adv_string == "Upwind_5th"  ) ||
             (dryscal_vert_adv_string == "Blended_5th6th") ||
             (dryscal_vert_adv_string == "Centered_6th") ||
             (dryscal_vert_adv_string == "WENO3"       ) ||
             (dryscal_vert_adv_string == "WENOZ3"      ) ||
             (dryscal_vert_adv_string == "WENOMZQ3"     ) ||
             (dryscal_vert_adv_string == "WENO5"       ) ||
             (dryscal_vert_adv_string == "WENOZ5"      ) )
        {
            dryscal_vert_adv_type = adv_type_convert_string_to_advtype(dryscal_vert_adv_string);
        }

        if ( (moistscal_horiz_adv_string == "Centered_2nd") ||
             (moistscal_horiz_adv_string == "Upwind_3rd"  ) ||
             (moistscal_horiz_adv_string == "Blended_3rd4th") ||
             (moistscal_horiz_adv_string == "Centered_4th") ||
             (moistscal_horiz_adv_string == "Upwind_5th"  ) ||
             (moistscal_horiz_adv_string == "Blended_5th6th") ||
             (moistscal_horiz_adv_string == "Centered_6th") ||
             (moistscal_horiz_adv_string == "WENO3"       ) ||
             (moistscal_horiz_adv_string == "WENOZ3"      ) ||
             (moistscal_horiz_adv_string == "WENOMZQ3"    ) ||
             (moistscal_horiz_adv_string == "WENO5"       ) ||
             (moistscal_horiz_adv_string == "WENOZ5"      ) )
        {
            moistscal_horiz_adv_type = adv_type_convert_string_to_advtype(moistscal_horiz_adv_string);
        }

        if ( (moistscal_vert_adv_string == "Centered_2nd") ||
             (moistscal_vert_adv_string == "Upwind_3rd"  ) ||
             (moistscal_vert_adv_string == "Blended_3rd4th") ||
             (moistscal_vert_adv_string == "Centered_4th") ||
             (moistscal_vert_adv_string == "Upwind_5th"  ) ||
             (moistscal_vert_adv_string == "Blended_5th6th") ||
             (moistscal_vert_adv_string == "Centered_6th") ||
             (moistscal_vert_adv_string == "WENO3"       ) ||
             (moistscal_vert_adv_string == "WENOZ3"      ) ||
             (moistscal_vert_adv_string == "WENOMZQ3"    ) ||
             (moistscal_vert_adv_string == "WENO5"       ) ||
             (moistscal_vert_adv_string == "WENOZ5"      ) )
        {
            moistscal_vert_adv_type = adv_type_convert_string_to_advtype(moistscal_vert_adv_string);
        }

        pp.queryarr("zero_xflux_faces", zero_xflux);
        pp.queryarr("zero_yflux_faces", zero_yflux);
        pp.queryarr("zero_zflux_faces", zero_zflux);
    }

    void display()
    {
        amrex::Print() << "Advection Choices: " << std::endl;
        amrex::Print() << "    dycore_horiz_adv_type       : " << adv_type_convert_int_to_string(dycore_horiz_adv_type);
        if (dycore_horiz_upw_frac < 1) amrex::Print() << " with " << 100*dycore_horiz_upw_frac << "% upwinding";
        amrex::Print() << std::endl;
        amrex::Print() << "    dycore_vert_adv_type        : " << adv_type_convert_int_to_string(dycore_vert_adv_type);
        if (dycore_vert_upw_frac < 1) amrex::Print() << " with " << 100*dycore_vert_upw_frac << "% upwinding";
        amrex::Print() << std::endl;
        amrex::Print() << "    dryscal_horiz_adv_type      : " << adv_type_convert_int_to_string(dryscal_horiz_adv_type);
        if (dryscal_horiz_upw_frac < 1) amrex::Print() << " with " << 100*dryscal_horiz_upw_frac << "% upwinding";
        amrex::Print() << std::endl;
        amrex::Print() << "    dryscal_vert_adv_type       : " << adv_type_convert_int_to_string(dryscal_vert_adv_type);
        if (dryscal_vert_upw_frac < 1) amrex::Print() << " with " << 100*dryscal_vert_upw_frac << "% upwinding";
        amrex::Print() << std::endl;
        amrex::Print() << "    moistscal_horiz_adv_type    : " << adv_type_convert_int_to_string(moistscal_horiz_adv_type);
        if (moistscal_horiz_upw_frac < 1) amrex::Print() << " with " << 100*moistscal_horiz_upw_frac << "% upwinding";
        amrex::Print() << std::endl;
        amrex::Print() << "    moistscal_vert_adv_type     : " << adv_type_convert_int_to_string(moistscal_vert_adv_type);
        if (moistscal_vert_upw_frac < 1) amrex::Print() << " with " << 100*moistscal_vert_upw_frac << "% upwinding";
        amrex::Print() << std::endl;
    }

    std::string
    adv_type_convert_int_to_string (AdvType adv_int)
    {
        if (adv_int == AdvType::Centered_2nd) {
            return "Centered_2nd";
        } else if (adv_int == AdvType::Upwind_3rd) {
            return "Upwind_3rd";
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
        } else {
            return AdvType::Unknown;
        }
    }

    // Default prefix
    std::string pp_prefix {"erf"};

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

    // Thin immersed body
    amrex::Vector<amrex::IntVect> zero_xflux;
    amrex::Vector<amrex::IntVect> zero_yflux;
    amrex::Vector<amrex::IntVect> zero_zflux;
};
#endif

#ifndef _ADV_STRUCT_H_
#define _ADV_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>

#include <ERF_Constants.H>
#include <IndexDefines.H>

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
           amrex::Print() << "Using efficient advection scheme" << std::endl;;
        }

        if ( (dycore_horiz_adv_string == "Centered_2nd") ||
             (dycore_horiz_adv_string == "Upwind_3rd"  ) ||
             (dycore_horiz_adv_string == "Centered_4th") ||
             (dycore_horiz_adv_string == "Upwind_5th"  ) ||
             (dycore_horiz_adv_string == "Centered_6th") )
        {
            dycore_horiz_adv_type = adv_type_convert_string_to_advtype(dycore_horiz_adv_string);
            amrex::Print() << "Using dycore_horiz_adv_type: " << dycore_horiz_adv_string << std::endl;
        } else {
            amrex::Print() << "Using default dycore_horiz_adv_type" << std::endl;;
        }

        if ( (dycore_vert_adv_string == "Centered_2nd") ||
             (dycore_vert_adv_string == "Upwind_3rd"  ) ||
             (dycore_vert_adv_string == "Centered_4th") ||
             (dycore_vert_adv_string == "Upwind_5th"  ) ||
             (dycore_vert_adv_string == "Centered_6th") )
        {
            dycore_vert_adv_type = adv_type_convert_string_to_advtype(dycore_vert_adv_string);
            amrex::Print() << "Using dycore_vert_adv_type: " << dycore_vert_adv_string << std::endl;
        } else {
            amrex::Print() << "Using default dycore_vert_adv_type" << std::endl;;
        }

        if ( (dryscal_horiz_adv_string == "Centered_2nd") ||
             (dryscal_horiz_adv_string == "Upwind_3rd"  ) ||
             (dryscal_horiz_adv_string == "Centered_4th") ||
             (dryscal_horiz_adv_string == "Upwind_5th"  ) ||
             (dryscal_horiz_adv_string == "Centered_6th") ||
             (dryscal_horiz_adv_string == "WENO3"       ) ||
             (dryscal_horiz_adv_string == "WENOZ3"      ) ||
             (dryscal_horiz_adv_string == "WENOMZQ3"    ) ||
             (dryscal_horiz_adv_string == "WENO5"       ) ||
             (dryscal_horiz_adv_string == "WENOZ5"      ) )
        {
            dryscal_horiz_adv_type = adv_type_convert_string_to_advtype(dryscal_horiz_adv_string);
            amrex::Print() << "Using dryscal_horiz_adv_type: " << dryscal_horiz_adv_string << std::endl;
        } else {
            amrex::Print() << "Using default dryscal_horiz_adv_type" << std::endl;;
        }

        if ( (dryscal_vert_adv_string == "Centered_2nd") ||
             (dryscal_vert_adv_string == "Upwind_3rd"  ) ||
             (dryscal_vert_adv_string == "Centered_4th") ||
             (dryscal_vert_adv_string == "Upwind_5th"  ) ||
             (dryscal_vert_adv_string == "Centered_6th") ||
             (dryscal_vert_adv_string == "WENO3"       ) ||
             (dryscal_vert_adv_string == "WENOZ3"      ) ||
             (dryscal_vert_adv_string == "WENOMZQ3"     ) ||
             (dryscal_vert_adv_string == "WENO5"       ) ||
             (dryscal_vert_adv_string == "WENOZ5"      ) )
        {
            dryscal_vert_adv_type = adv_type_convert_string_to_advtype(dryscal_vert_adv_string);
            amrex::Print() << "Using dryscal_vert_adv_type: " << dryscal_vert_adv_string << std::endl;
        } else {
            amrex::Print() << "Using default dryscal_vert_adv_type" << std::endl;;
        }

        if ( (moistscal_horiz_adv_string == "Centered_2nd") ||
             (moistscal_horiz_adv_string == "Upwind_3rd"  ) ||
             (moistscal_horiz_adv_string == "Centered_4th") ||
             (moistscal_horiz_adv_string == "Upwind_5th"  ) ||
             (moistscal_horiz_adv_string == "Centered_6th") ||
             (moistscal_horiz_adv_string == "WENO3"       ) ||
             (moistscal_horiz_adv_string == "WENOZ3"      ) ||
             (moistscal_horiz_adv_string == "WENOMZQ3"    ) ||
             (moistscal_horiz_adv_string == "WENO5"       ) ||
             (moistscal_horiz_adv_string == "WENOZ5"      ) )
        {
            moistscal_horiz_adv_type = adv_type_convert_string_to_advtype(moistscal_horiz_adv_string);
            amrex::Print() << "Using moistscal_horiz_adv_type: " << moistscal_horiz_adv_string << std::endl;
        } else {
            amrex::Print() << "Using default moistscal_horiz_adv_type" << std::endl;;
        }

        if ( (moistscal_vert_adv_string == "Centered_2nd") ||
             (moistscal_vert_adv_string == "Upwind_3rd"  ) ||
             (moistscal_vert_adv_string == "Centered_4th") ||
             (moistscal_vert_adv_string == "Upwind_5th"  ) ||
             (moistscal_vert_adv_string == "Centered_6th") ||
             (moistscal_vert_adv_string == "WENO3"       ) ||
             (moistscal_vert_adv_string == "WENOZ3"      ) ||
             (moistscal_vert_adv_string == "WENOMZQ3"    ) ||
             (moistscal_vert_adv_string == "WENO5"       ) ||
             (moistscal_vert_adv_string == "WENOZ5"      ) )
        {
            moistscal_vert_adv_type = adv_type_convert_string_to_advtype(moistscal_vert_adv_string);
            amrex::Print() << "Using moistscal_vert_adv_type: " << moistscal_vert_adv_string << std::endl;
        } else {
            amrex::Print() << "Using default moistscal_vert_adv_type" << std::endl;;
        }
    }

    void display()
    {
        amrex::Print() << "Advection Choices: " << std::endl;
        amrex::Print() << "dycore_horiz_adv_type       : " << adv_type_convert_int_to_string(dycore_horiz_adv_type) << std::endl;
        amrex::Print() << "dycore_vert_adv_type        : " << adv_type_convert_int_to_string(dycore_vert_adv_type) << std::endl;
        amrex::Print() << "dryscal_horiz_adv_type      : " << adv_type_convert_int_to_string(dryscal_horiz_adv_type) << std::endl;
        amrex::Print() << "dryscal_vert_adv_type       : " << adv_type_convert_int_to_string(dryscal_vert_adv_type) << std::endl;
        amrex::Print() << "moistscal_horiz_adv_type    : " << adv_type_convert_int_to_string(moistscal_horiz_adv_type) << std::endl;
        amrex::Print() << "moistscal_vert_adv_type     : " << adv_type_convert_int_to_string(moistscal_vert_adv_type) << std::endl;
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
        } else if (adv_string == "Upwind_3rd") {
            return AdvType::Upwind_3rd;
        } else if (adv_string == "Centered_4th") {
            return AdvType::Centered_4th;
        } else if (adv_string == "Upwind_5th") {
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

    // Spatial discretization
    // Order and type of spatial discretizations used in advection
    bool use_efficient_advection = false;
    AdvType dycore_horiz_adv_type    = AdvType::Centered_2nd;
    AdvType dycore_vert_adv_type     = AdvType::Centered_2nd;
    AdvType dryscal_horiz_adv_type   = AdvType::Centered_2nd;
    AdvType dryscal_vert_adv_type    = AdvType::Centered_2nd;
    AdvType moistscal_horiz_adv_type = AdvType::Weno_3;
    AdvType moistscal_vert_adv_type  = AdvType::Weno_3;
};
#endif

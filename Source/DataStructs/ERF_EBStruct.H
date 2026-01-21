#ifndef ERF_EB_STRUCT_H_
#define ERF_EB_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

enum struct EBBoundaryType {
    SlipWall, NoSlipWall
};

/**
 * Container holding EB-related choices
 */

struct EBChoice {
  public:
    void init_params(std::string pp_prefix)
    {
        amrex::ParmParse pp(pp_prefix);

        static std::string eb_boundary_type_string = "NoSlipWall";
        pp.query("eb_boundary_type",eb_boundary_type_string);

        if (!eb_boundary_type_string.compare("SlipWall")) {
            eb_boundary_type = EBBoundaryType::SlipWall;
        } else if (!eb_boundary_type_string.compare("NoSlipWall")) {
            eb_boundary_type = EBBoundaryType::NoSlipWall;
        } else {
            amrex::Error("Don't know this eb_boundary_type");
        }
    }

    void display()
    {
        amrex::Print() << "EB choices: " << std::endl;
        if (eb_boundary_type == EBBoundaryType::SlipWall) {
            amrex::Print() << "    eb_boundary_type        : SlipWall" << std::endl;
        } else if (eb_boundary_type == EBBoundaryType::NoSlipWall) {
            amrex::Print() << "    eb_boundary_type        : NoSlipWall" << std::endl;
        }
    }

    EBBoundaryType eb_boundary_type = EBBoundaryType::SlipWall;
};
#endif

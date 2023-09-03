#ifndef _DIFF_STRUCT_H_
#define _DIFF_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

enum struct MolecDiffType {
    None, Constant, ConstantAlpha
};

/**
 * Container holding diffusion-related choices
 */

struct DiffChoice {
  public:
    void init_params()
    {
        amrex::ParmParse pp(pp_prefix);

        // Diffusive constants...
        pp.query("alpha_T", alpha_T);
        pp.query("alpha_C", alpha_C);
        pp.query("dynamicViscosity", dynamicViscosity);
        pp.query("rho0_trans", rho0_trans);

        static std::string molec_diff_type_string = "None";
        pp.query("molec_diff_type",molec_diff_type_string);

        if (!molec_diff_type_string.compare("Constant")) {
            molec_diff_type = MolecDiffType::Constant;
        } else if (!molec_diff_type_string.compare("ConstantAlpha")) {
            molec_diff_type = MolecDiffType::ConstantAlpha;
        } else if (!molec_diff_type_string.compare("None")) {
            molec_diff_type = MolecDiffType::None;
        } else {
            amrex::Error("Don't know this molec_diff_type");
        }

        if (molec_diff_type == MolecDiffType::ConstantAlpha) {
            amrex::Print() << "Using constant kinematic diffusion coefficients" << std::endl;
            amrex::Print() << "  momentum : " << dynamicViscosity/rho0_trans << " m^2/s" << std::endl;
            amrex::Print() << "  temperature : " << alpha_T << " m^2/s" << std::endl;
            amrex::Print() << "  scalar : " << alpha_C << " m^2/s" << std::endl;
        }

        // Compute relevant forms of diffusion parameters
        rhoAlpha_T = rho0_trans * alpha_T;
        rhoAlpha_C = rho0_trans * alpha_C;
    }

    void display()
    {
        amrex::Print() << "Diffusion choices: " << std::endl;
        amrex::Print() << "rho0_trans                  : " << rho0_trans << std::endl;
        amrex::Print() << "alpha_T                     : " << alpha_T << std::endl;
        amrex::Print() << "alpha_C                     : " << alpha_C << std::endl;
        amrex::Print() << "dynamicViscosity            : " << dynamicViscosity << std::endl;

        if (molec_diff_type == MolecDiffType::Constant) {
            amrex::Print() << "Using constant molecular diffusivity (relevant for DNS)" << std::endl;
        } else if (molec_diff_type == MolecDiffType::None) {
            amrex::Print() << "Not using any molecular diffusivity, i.e. using the modeled turbulent diffusivity"
            << std::endl;
        }
    }

    // Default prefix
    std::string pp_prefix {"erf"};

    // Molecular transport model
    MolecDiffType molec_diff_type = MolecDiffType::Constant;

    // Diffusive/viscous coefficients [m2/s]
    amrex::Real alpha_T = 0.0;
    amrex::Real alpha_C = 0.0;

    // Density for computation of rhoAlpha (which is assumed constant) [kg/m3]
    amrex::Real rho0_trans = 1.0;

    // Dynamic diffusion coefficients [kg/(m-s)]
    amrex::Real rhoAlpha_T = 0.0;
    amrex::Real rhoAlpha_C = 0.0;
    amrex::Real dynamicViscosity = 0.0;
};
#endif

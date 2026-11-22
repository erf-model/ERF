#ifndef ERF_EB_STRUCT_H_
#define ERF_EB_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

/**
 * @brief Boundary condition model applied on embedded-boundary surfaces.
 */
enum struct EBBoundaryType {
    SlipWall, NoSlipWall, SurfaceLayer
};

/**
 * @brief Embedded-boundary stress component indices.
 */
enum EBTauType {
    tau_eb13, tau_eb23
};

/**
 * @brief Grid centering used for embedded-boundary data.
 */
enum EBGridType {
    xface, yface, zface
};

/**
 * Container holding EB-related choices
 */

struct EBChoice {
  public:
    /**
     * @brief Read embedded-boundary options from the input parameter database.
     * @param pp_prefix ParmParse prefix for the ERF input namespace.
     * @param requires_cover Force eb_cover_zero_area_faces on (a level projects every step: the
     *        projection does not constrain momentum on faces of zero CC area fraction).
     */
    void init_params (std::string pp_prefix, bool requires_cover)
    {
        amrex::ParmParse pp(pp_prefix);

        static std::string eb_boundary_type_string = "NoSlipWall";
        pp.queryAdd("eb_boundary_type",eb_boundary_type_string);

        if (!eb_boundary_type_string.compare("SlipWall")) {
            eb_boundary_type = EBBoundaryType::SlipWall;
        } else if (!eb_boundary_type_string.compare("NoSlipWall")) {
            eb_boundary_type = EBBoundaryType::NoSlipWall;
        } else if (!eb_boundary_type_string.compare("SurfaceLayer")) {
            eb_boundary_type = EBBoundaryType::SurfaceLayer;
        } else {
            amrex::Error("Don't know this eb_boundary_type");
        }

        const bool given = pp.query("eb_cover_zero_area_faces", eb_cover_zero_area_faces);
        if (requires_cover) {
            if (given && !eb_cover_zero_area_faces) {
                amrex::Abort("erf.eb_cover_zero_area_faces = false is not allowed with Anelastic or "
                             "PseudoIncompressible levels");
            }
            eb_cover_zero_area_faces = true;
        }
    }

    /**
     * @brief Print the configured embedded-boundary options.
     */
    void display ()
    {
        amrex::Print() << "EB choices: " << std::endl;
        if (eb_boundary_type == EBBoundaryType::SlipWall) {
            amrex::Print() << "    eb_boundary_type        : SlipWall" << std::endl;
        } else if (eb_boundary_type == EBBoundaryType::NoSlipWall) {
            amrex::Print() << "    eb_boundary_type        : NoSlipWall" << std::endl;
        } else if (eb_boundary_type == EBBoundaryType::SurfaceLayer) {
            amrex::Print() << "    eb_boundary_type        : SurfaceLayer" << std::endl;
        }
        amrex::Print() << "    eb_cover_zero_area_faces: " << (eb_cover_zero_area_faces ? "true" : "false")
                       << std::endl;
    }

    EBBoundaryType eb_boundary_type = EBBoundaryType::SlipWall;
    bool eb_cover_zero_area_faces = false;
};

/**
 * @var EBChoice::eb_boundary_type
 * @brief Boundary condition model applied on embedded-boundary surfaces.
 *
 * @var EBChoice::eb_cover_zero_area_faces
 * @brief Cover FC EB cells at faces with zero CC area fraction (forced on with a projection).
 */
#endif

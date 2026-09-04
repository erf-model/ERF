#ifndef ERF_IBSEB_MATERIALS_H
#define ERF_IBSEB_MATERIALS_H

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <string>

/**
 * \file ERF_IBSEBMaterials.H
 * \brief Material library of the building faces.
 *
 * The CSV schema is the one of the ERF-SLUCM branch's material library, so
 * the same file serves both models:
 *
 *     mat_id,name,albedo,emissivity,k_therm_W_per_mK,rho_cp_J_per_m3K,thickness_m,description
 *
 * The I/O rank reads the file and broadcasts the entries as plain data;
 * every rank then holds the same table. A face gets its material from its
 * building (``erf.ibseb.material_by_building``) or the default
 * (``erf.ibseb.material_default``); without a file every face uses the
 * uniform properties of the inputs.
 */
struct IBSEBMaterial
{
    int         mat_id      = 0;
    amrex::Real albedo      = 0.3;
    amrex::Real emissivity  = 0.9;
    amrex::Real k_therm     = 1.0;    ///< [W/m/K]
    amrex::Real rho_cp      = 1.6e6;  ///< [J/m3/K]
    amrex::Real thickness   = 0.3;    ///< [m]
    char        name[64]    = "uniform";
};

/**
 * Read the material CSV on the I/O rank and broadcast it. Aborts on a
 * missing file, a header that does not match the schema, a duplicate id,
 * or a value outside its physical range.
 */
amrex::Vector<IBSEBMaterial> read_ibseb_materials (const std::string& file);

#endif

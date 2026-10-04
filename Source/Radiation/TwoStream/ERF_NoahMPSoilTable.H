#ifndef ERF_NOAHMP_SOIL_TABLE_H_
#define ERF_NOAHMP_SOIL_TABLE_H_

/**
 * @file ERF_NoahMPSoilTable.H
 * @brief Noah-MP's soil parameters by soil category, for the two-stream surface energy
 * balance (erf.radiation.seb_soil_type).
 *
 * The 19 categories and the values of the STAS soil dataset (&noahmp_soil_stas_parameters)
 * in Submodules/Noah-MP/parameters/NoahmpTable.TBL -- the table Noah-MP itself reads --
 * copied here so the balance can use them in builds without Noah-MP. The CTest
 * NoahMPSoilTable_MatchesSubmodule (Tests/check_noahmp_soil_table.py) fails if the two
 * drift apart.
 */

#include <AMReX_REAL.H>

/** One soil category's parameters, named as in NoahmpTable.TBL. */
struct NoahMPSoilParams
{
    amrex::Real bb; ///< BB Clapp-Hornberger b [1]
    amrex::Real smc_dry; ///< DRYSMC dry soil moisture threshold [m^3/m^3]
    amrex::Real smc_max; ///< MAXSMC saturated soil moisture (porosity) [m^3/m^3]
    amrex::Real smc_ref; ///< REFSMC reference soil moisture, field capacity [m^3/m^3]
    amrex::Real psi_sat; ///< SATPSI saturated matric potential [m]
    amrex::Real k_sat; ///< SATDK saturated hydraulic conductivity [m/s]
    amrex::Real smc_wilt; ///< WLTSMC wilting point [m^3/m^3]
    amrex::Real quartz; ///< QTZ quartz fraction [1]
};

/** Number of STAS soil categories; categories are numbered from 1, as in Noah-MP. */
inline constexpr int noahmp_num_soil_categories = 19;

/** Noah-MP's soil heat capacity CSOIL_DATA [J/m^3/K] (&noahmp_general_parameters). */
inline constexpr amrex::Real noahmp_soil_heat_capacity = amrex::Real(2.00E+6);

/** Noah-MP's RSURF_EXP [1] (&noahmp_global_parameters): the exponent of the dry surface
 * layer's thickness in its ground evaporation resistance option 1 (Sakaguchi and Zeng). */
inline constexpr amrex::Real noahmp_soil_resistance_exponent = amrex::Real(5.0);

/** Noah-MP's Z0SOIL [m] (&noahmp_global_parameters): the bare-soil roughness length. */
inline constexpr amrex::Real noahmp_soil_roughness = amrex::Real(0.002);

/**
 * @brief The parameters of STAS soil category @p category (1-19), as Noah-MP reads them.
 * @return nullptr for a category outside 1-19
 */
[[nodiscard]] inline const NoahMPSoilParams*
noahmp_soil_params (int category) noexcept
{
    // {BB, DRYSMC, MAXSMC, REFSMC, SATPSI, SATDK, WLTSMC, QTZ}
    static constexpr NoahMPSoilParams table[noahmp_num_soil_categories] = {
    {amrex::Real(2.790), amrex::Real(0.010), amrex::Real(0.339), amrex::Real(0.192), amrex::Real(0.069), amrex::Real(4.66E-05), amrex::Real(0.010), amrex::Real(0.920)}, //  1 SAND
    {amrex::Real(4.260), amrex::Real(0.028), amrex::Real(0.421), amrex::Real(0.283), amrex::Real(0.036), amrex::Real(1.41E-05), amrex::Real(0.028), amrex::Real(0.820)}, //  2 LOAMY SAND
    {amrex::Real(4.740), amrex::Real(0.047), amrex::Real(0.434), amrex::Real(0.312), amrex::Real(0.141), amrex::Real(5.23E-06), amrex::Real(0.047), amrex::Real(0.600)}, //  3 SANDY LOAM
    {amrex::Real(5.330), amrex::Real(0.084), amrex::Real(0.476), amrex::Real(0.360), amrex::Real(0.759), amrex::Real(2.81E-06), amrex::Real(0.084), amrex::Real(0.250)}, //  4 SILT LOAM
    {amrex::Real(3.860), amrex::Real(0.061), amrex::Real(0.484), amrex::Real(0.347), amrex::Real(0.955), amrex::Real(2.18E-06), amrex::Real(0.061), amrex::Real(0.100)}, //  5 SILT
    {amrex::Real(5.250), amrex::Real(0.066), amrex::Real(0.439), amrex::Real(0.329), amrex::Real(0.355), amrex::Real(3.38E-06), amrex::Real(0.066), amrex::Real(0.400)}, //  6 LOAM
    {amrex::Real(6.770), amrex::Real(0.069), amrex::Real(0.404), amrex::Real(0.315), amrex::Real(0.135), amrex::Real(4.45E-06), amrex::Real(0.069), amrex::Real(0.600)}, //  7 SANDY CLAY LOAM
    {amrex::Real(8.720), amrex::Real(0.120), amrex::Real(0.464), amrex::Real(0.387), amrex::Real(0.617), amrex::Real(2.03E-06), amrex::Real(0.120), amrex::Real(0.100)}, //  8 SILTY CLAY LOAM
    {amrex::Real(8.170), amrex::Real(0.103), amrex::Real(0.465), amrex::Real(0.382), amrex::Real(0.263), amrex::Real(2.45E-06), amrex::Real(0.103), amrex::Real(0.350)}, //  9 CLAY LOAM
    {amrex::Real(10.730), amrex::Real(0.100), amrex::Real(0.406), amrex::Real(0.338), amrex::Real(0.098), amrex::Real(7.22E-06), amrex::Real(0.100), amrex::Real(0.520)}, // 10 SANDY CLAY
    {amrex::Real(10.390), amrex::Real(0.126), amrex::Real(0.468), amrex::Real(0.404), amrex::Real(0.324), amrex::Real(1.34E-06), amrex::Real(0.126), amrex::Real(0.100)}, // 11 SILTY CLAY
    {amrex::Real(11.550), amrex::Real(0.138), amrex::Real(0.468), amrex::Real(0.412), amrex::Real(0.468), amrex::Real(9.74E-07), amrex::Real(0.138), amrex::Real(0.250)}, // 12 CLAY
    {amrex::Real(5.250), amrex::Real(0.066), amrex::Real(0.439), amrex::Real(0.329), amrex::Real(0.355), amrex::Real(3.38E-06), amrex::Real(0.066), amrex::Real(0.050)}, // 13 ORGANIC MATERIAL
    {amrex::Real(0.000), amrex::Real(0.000), amrex::Real(1.000), amrex::Real(0.000), amrex::Real(0.000), amrex::Real(0.00E+00), amrex::Real(0.000), amrex::Real(0.600)}, // 14 WATER
    {amrex::Real(2.790), amrex::Real(0.006), amrex::Real(0.200), amrex::Real(0.170), amrex::Real(0.069), amrex::Real(1.41E-04), amrex::Real(0.006), amrex::Real(0.070)}, // 15 BEDROCK
    {amrex::Real(4.260), amrex::Real(0.028), amrex::Real(0.421), amrex::Real(0.283), amrex::Real(0.036), amrex::Real(1.41E-05), amrex::Real(0.028), amrex::Real(0.250)}, // 16 OTHER(land-ice)
    {amrex::Real(11.550), amrex::Real(0.030), amrex::Real(0.468), amrex::Real(0.454), amrex::Real(0.468), amrex::Real(9.74E-07), amrex::Real(0.030), amrex::Real(0.600)}, // 17 PLAYA
    {amrex::Real(2.790), amrex::Real(0.006), amrex::Real(0.200), amrex::Real(0.170), amrex::Real(0.069), amrex::Real(1.41E-04), amrex::Real(0.006), amrex::Real(0.520)}, // 18 LAVA
    {amrex::Real(2.790), amrex::Real(0.010), amrex::Real(0.339), amrex::Real(0.192), amrex::Real(0.069), amrex::Real(4.66E-05), amrex::Real(0.010), amrex::Real(0.920)}, // 19 WHITE SAND
    };
    if (category < 1 || category > noahmp_num_soil_categories) { return nullptr; }
    return &table[category - 1];
}

#endif

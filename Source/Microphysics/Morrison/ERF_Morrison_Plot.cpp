/**
 * \file ERF_Morrison_Plot.cpp
 * \brief Plotting and diagnostic output for Morrison two-moment microphysics
 *
 * This file provides the interface for outputting Morrison microphysics variables
 * to plotfiles for visualization and analysis. It enables extraction of internal
 * Morrison state variables that are not part of the main ERF conserved state.
 */

#include "ERF_Morrison.H"

using namespace amrex;

namespace {

/**
 * \struct PlotVarEntry
 * \brief Maps user-friendly plot variable names to internal Morrison variable indices
 */
struct PlotVarEntry {
    const char* name;  ///< Variable name as it appears in plotfiles
    int index;         ///< Index into MicVar_Morr enum
};

/**
 * \brief Table of all Morrison variables available for plotting
 *
 * This table defines which Morrison microphysics variables can be output to
 * plotfiles and their corresponding names. Variables are organized by category:
 *
 * **Thermodynamic state variables:**
 * - micro_rho: Air density (kg/m³)
 * - micro_theta: Potential temperature (K)
 * - micro_temp: Absolute temperature (K)
 * - micro_pres: Pressure (Pa)
 *
 * **Non-precipitating moisture variables (mixing ratios in kg/kg):**
 * - micro_qv: Water vapor mixing ratio
 * - micro_qc: Cloud liquid water mixing ratio
 * - micro_qi: Cloud ice mixing ratio
 * - micro_qn: Total cloud condensate (qc + qi)
 * - micro_qt: Total water mixing ratio (qv + qn)
 *
 * **Precipitating hydrometeor variables (mixing ratios in kg/kg):**
 * - micro_qp: Total precipitation (qrain + qsnow + qgraup)
 * - micro_qrain: Rain water mixing ratio
 * - micro_qsnow: Snow mixing ratio
 * - micro_qgraup: Graupel mixing ratio
 *
 * **Number concentrations (1/kg):**
 * - micro_nc: Cloud droplet number concentration
 * - micro_nr: Rain drop number concentration
 * - micro_ni: Cloud ice number concentration
 * - micro_ns: Snow number concentration
 * - micro_ng: Graupel number concentration
 *
 * **Dynamical variables:**
 * - micro_omega: Grid-scale vertical velocity (m/s)
 *   Note: This field is populated from ERF's velocity state and represents
 *   the vertical velocity used as input to the Morrison microphysics scheme.
 *   It is initialized to zero and must be explicitly populated before plotting.
 */
constexpr PlotVarEntry plot_entries[] = {
    {"micro_rho",    MicVar_Morr::rho},      // Air density
    {"micro_theta",  MicVar_Morr::theta},    // Potential temperature
    {"micro_temp",   MicVar_Morr::tabs},     // Absolute temperature
    {"micro_pres",   MicVar_Morr::pres},     // Pressure
    {"micro_qv",     MicVar_Morr::qv},       // Water vapor
    {"micro_qc",     MicVar_Morr::qcl},      // Cloud liquid
    {"micro_qi",     MicVar_Morr::qci},      // Cloud ice
    {"micro_qn",     MicVar_Morr::qn},       // Total cloud condensate
    {"micro_qt",     MicVar_Morr::qt},       // Total water
    {"micro_qp",     MicVar_Morr::qp},       // Total precipitation
    {"micro_qrain",  MicVar_Morr::qpr},      // Rain
    {"micro_qsnow",  MicVar_Morr::qps},      // Snow
    {"micro_qgraup", MicVar_Morr::qpg},      // Graupel
    {"micro_nc",     MicVar_Morr::nc},       // Cloud droplet number
    {"micro_nr",     MicVar_Morr::nr},       // Rain number
    {"micro_ni",     MicVar_Morr::ni},       // Ice number
    {"micro_ns",     MicVar_Morr::ns},       // Snow number
    {"micro_ng",     MicVar_Morr::ng},       // Graupel number
    {"micro_omega",  MicVar_Morr::omega}     // Vertical velocity
};

} // namespace

/**
 * \brief Populate a vector with names of all available Morrison plot variables
 *
 * This function returns the names of all Morrison microphysics variables that
 * can be written to plotfiles. These names correspond to the diagnostic fields
 * defined in the plot_entries table.
 *
 * \param[out] a_vec Vector to be populated with plot variable names
 *
 * \note The returned vector will contain 19 variable names, including both
 *       prognostic quantities (mixing ratios, number concentrations) and
 *       diagnostic quantities (temperature, pressure, derived moisture totals).
 */
void
Morrison::GetPlotVarNames (Vector<std::string>& a_vec) const
{
    a_vec.clear();
    a_vec.reserve(sizeof(plot_entries) / sizeof(plot_entries[0]));
    for (const auto& entry : plot_entries) {
        a_vec.emplace_back(entry.name);
    }
}

/**
 * \brief Extract a Morrison microphysics variable for plotting
 *
 * This function copies data from the Morrison internal state (mic_fab_vars)
 * to the provided MultiFab for inclusion in a plotfile. The variable to extract
 * is specified by name matching the plot_entries table.
 *
 * \param[in]  a_name Name of the variable to extract (e.g., "micro_qc")
 * \param[out] a_mf   MultiFab to receive the variable data
 *
 * \throws amrex::Abort if the requested variable name is not found in plot_entries
 *
 * \note The caller must ensure a_mf has the correct distribution and size
 *       matching the Morrison mic_fab_vars arrays. Only the first component
 *       of the source MultiFab is copied.
 *
 * \warning If omega has not been populated from the velocity field, it will
 *          contain its initialization value (zero).
 */
void
Morrison::GetPlotVar (const std::string& a_name,
                      MultiFab& a_mf) const
{
    const MultiFab* src = nullptr;

    for (const auto& entry : plot_entries) {
        if (a_name == entry.name) {
            src = mic_fab_vars[entry.index].get();
            break;
        }
    }

    if (src == nullptr) {
        amrex::Abort("Morrison::GetPlotVar: unknown plot variable " + a_name);
    }

    MultiFab::Copy(a_mf, *src, 0, 0, 1, 0);
}

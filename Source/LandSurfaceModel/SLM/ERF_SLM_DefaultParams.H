#ifndef ERF_SLM_DEFAULT_PARAMS_H
#define ERF_SLM_DEFAULT_PARAMS_H

#include <array>

#include <AMReX_REAL.H>

namespace SLMDefaultParams {

constexpr int num_veg_params = 20;

// Values copied from the MODIS/STAS sections of NoahmpTable.TBL.
inline constexpr std::array<amrex::Real, num_veg_params> cwpvt = {
    0.18, 0.67, 0.18, 0.67, 0.29, 1.0, 2.0, 1.3, 1.0, 5.0,
    1.17, 1.67, 1.67, 1.67, 0.18, 0.18, 0.18, 0.67, 1.0, 0.18};
inline constexpr std::array<amrex::Real, num_veg_params> den = {
    0.28, 0.02, 0.28, 0.10, 0.10, 10.0, 10.0, 10.0, 0.02, 100.0,
    5.05, 25.0, 0.01, 25.0, 0.0, 0.01, 0.01, 1.0, 1.0, 1.0};
inline constexpr std::array<amrex::Real, num_veg_params> hvb = {
    8.50, 8.00, 7.00, 11.5, 10.0, 0.10, 0.10, 0.10, 0.10, 0.05,
    0.10, 0.10, 1.00, 0.10, 0.0, 0.0, 0.0, 0.30, 0.20, 0.10};
inline constexpr std::array<amrex::Real, num_veg_params> rc = {
    1.20, 3.60, 1.20, 1.40, 1.40, 0.12, 0.12, 0.12, 3.00, 0.03,
    0.75, 0.08, 1.00, 0.08, 0.0, 0.01, 0.01, 0.30, 0.30, 0.30};
inline constexpr std::array<amrex::Real, num_veg_params> rhol_vis = {
    0.07, 0.10, 0.07, 0.10, 0.10, 0.07, 0.07, 0.07, 0.10, 0.11,
    0.105, 0.11, 0.0, 0.11, 0.0, 0.0, 0.0, 0.10, 0.10, 0.10};
inline constexpr std::array<amrex::Real, num_veg_params> rhol_nir = {
    0.35, 0.45, 0.35, 0.45, 0.45, 0.35, 0.35, 0.35, 0.45, 0.58,
    0.515, 0.58, 0.0, 0.58, 0.0, 0.0, 0.0, 0.45, 0.45, 0.45};
inline constexpr std::array<amrex::Real, num_veg_params> rhos_vis = {
    0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.36,
    0.26, 0.36, 0.0, 0.36, 0.0, 0.0, 0.0, 0.16, 0.16, 0.16};
inline constexpr std::array<amrex::Real, num_veg_params> rhos_nir = {
    0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.58,
    0.485, 0.58, 0.0, 0.58, 0.0, 0.0, 0.0, 0.39, 0.39, 0.39};
inline constexpr std::array<amrex::Real, num_veg_params> taul_vis = {
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.07,
    0.06, 0.07, 0.0, 0.07, 0.0, 0.0, 0.0, 0.05, 0.05, 0.05};
inline constexpr std::array<amrex::Real, num_veg_params> taul_nir = {
    0.10, 0.25, 0.10, 0.25, 0.25, 0.10, 0.10, 0.10, 0.25, 0.25,
    0.25, 0.25, 0.0, 0.25, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25};
inline constexpr std::array<amrex::Real, num_veg_params> taus_vis = {
    0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.220,
    0.1105, 0.220, 0.0, 0.220, 0.0, 0.0, 0.0, 0.001, 0.001, 0.001};
inline constexpr std::array<amrex::Real, num_veg_params> taus_nir = {
    0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.380,
    0.1905, 0.380, 0.0, 0.380, 0.0, 0.0, 0.0, 0.001, 0.001, 0.001};

inline constexpr std::array<amrex::Real, 8> albsat_vis =
    {0.15, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05};
inline constexpr std::array<amrex::Real, 8> albsat_nir =
    {0.30, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10};
inline constexpr std::array<amrex::Real, 8> albdry_vis =
    {0.27, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10};
inline constexpr std::array<amrex::Real, 8> albdry_nir =
    {0.54, 0.44, 0.40, 0.36, 0.32, 0.28, 0.24, 0.20};
inline constexpr std::array<amrex::Real, 2> alblak = {0.60, 0.40};
inline constexpr std::array<amrex::Real, 2> omegas = {0.8, 0.4};

} // namespace SLMDefaultParams

#endif

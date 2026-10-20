#ifndef ERF_NOAHMP_VEGETATION_TABLE_H_
#define ERF_NOAHMP_VEGETATION_TABLE_H_

/**
 * @file ERF_NoahMPVegetationTable.H
 * @brief Noah-MP's canopy-resistance parameters, roughness length and monthly leaf area index by land-use
 * category, for the two-stream surface energy balance (erf.radiation.seb_vegetation_type).
 *
 * The 20 categories of the modified IGBP MODIS land-use dataset and their values in
 * &noahmp_modis_parameters of Submodules/Noah-MP/parameters/NoahmpTable.TBL -- the table
 * Noah-MP itself reads -- copied here so the balance can use them in builds without
 * Noah-MP. The CTest NoahMPSoilTable_MatchesSubmodule (Tests/check_noahmp_soil_table.py)
 * fails if the two drift apart.
 */

#include <AMReX_REAL.H>

#include <cstdio>
#include <string>

/** One land-use category's parameters, named as in NoahmpTable.TBL. */
struct NoahMPVegetationParams
{
    amrex::Real rs_min;   ///< RS     minimum stomatal resistance [s/m]
    amrex::Real rgl;      ///< RGL    radiation stress parameter [W/m^2]
    amrex::Real hs;       ///< HS     vapour pressure deficit parameter [1/(kg/kg)]
    amrex::Real t_opt;    ///< TOPT   optimum transpiration air temperature [K]
    amrex::Real rs_max;   ///< RSMAX  maximum stomatal resistance [s/m]
    amrex::Real z0;       ///< Z0MVT  momentum roughness length of the vegetation [m]
    amrex::Real lai[12];  ///< LAI_JAN ... LAI_DEC  monthly leaf area index [m^2/m^2]
};

/** Number of MODIS land-use categories; categories are numbered from 1, as in Noah-MP. */
inline constexpr int noahmp_num_vegetation_categories = 20;

/**
 * @brief The parameters of MODIS land-use category @p category (1-20), as Noah-MP reads them.
 * @return nullptr for a category outside 1-20
 */
[[nodiscard]] inline const NoahMPVegetationParams*
noahmp_vegetation_params (int category) noexcept
{
    // {RS, RGL, HS, TOPT, RSMAX, Z0MVT, {LAI_JAN, ..., LAI_DEC}}
    static constexpr NoahMPVegetationParams table[noahmp_num_vegetation_categories] = {
        {amrex::Real(125.0), amrex::Real(30.0), amrex::Real(47.35), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(1.09),  //  1 Evergreen Needleleaf Forest
         {amrex::Real(4.0), amrex::Real(4.0), amrex::Real(4.0), amrex::Real(4.0), amrex::Real(4.0), amrex::Real(4.0), amrex::Real(4.0), amrex::Real(4.0), amrex::Real(4.0), amrex::Real(4.0), amrex::Real(4.0), amrex::Real(4.0)}},
        {amrex::Real(150.0), amrex::Real(30.0), amrex::Real(41.69), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(1.10),  //  2 Evergreen Broadleaf Forest
         {amrex::Real(4.5), amrex::Real(4.5), amrex::Real(4.5), amrex::Real(4.5), amrex::Real(4.5), amrex::Real(4.5), amrex::Real(4.5), amrex::Real(4.5), amrex::Real(4.5), amrex::Real(4.5), amrex::Real(4.5), amrex::Real(4.5)}},
        {amrex::Real(150.0), amrex::Real(30.0), amrex::Real(47.35), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.85),  //  3 Deciduous Needleleaf Forest
         {amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.6), amrex::Real(1.2), amrex::Real(2.0), amrex::Real(2.6), amrex::Real(1.7), amrex::Real(1.0), amrex::Real(0.5), amrex::Real(0.2), amrex::Real(0.0)}},
        {amrex::Real(100.0), amrex::Real(30.0), amrex::Real(54.53), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.80),  //  4 Deciduous Broadleaf Forest
         {amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.3), amrex::Real(1.2), amrex::Real(3.0), amrex::Real(4.7), amrex::Real(4.5), amrex::Real(3.4), amrex::Real(1.2), amrex::Real(0.3), amrex::Real(0.0), amrex::Real(0.0)}},
        {amrex::Real(125.0), amrex::Real(30.0), amrex::Real(51.93), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.80),  //  5 Mixed Forests
         {amrex::Real(2.0), amrex::Real(2.0), amrex::Real(2.2), amrex::Real(2.6), amrex::Real(3.5), amrex::Real(4.3), amrex::Real(4.3), amrex::Real(3.7), amrex::Real(2.6), amrex::Real(2.2), amrex::Real(2.0), amrex::Real(2.0)}},
        {amrex::Real(300.0), amrex::Real(100.0), amrex::Real(42.00), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.20),  //  6 Closed Shrublands
         {amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.3), amrex::Real(0.9), amrex::Real(2.2), amrex::Real(3.5), amrex::Real(3.5), amrex::Real(2.5), amrex::Real(0.9), amrex::Real(0.3), amrex::Real(0.0), amrex::Real(0.0)}},
        {amrex::Real(170.0), amrex::Real(100.0), amrex::Real(39.18), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.06),  //  7 Open Shrublands
         {amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.2), amrex::Real(0.6), amrex::Real(1.5), amrex::Real(2.3), amrex::Real(2.3), amrex::Real(1.7), amrex::Real(0.6), amrex::Real(0.2), amrex::Real(0.0), amrex::Real(0.0)}},
        {amrex::Real(300.0), amrex::Real(100.0), amrex::Real(42.00), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.60),  //  8 Woody Savannas
         {amrex::Real(0.2), amrex::Real(0.2), amrex::Real(0.4), amrex::Real(1.0), amrex::Real(2.4), amrex::Real(4.1), amrex::Real(4.1), amrex::Real(2.7), amrex::Real(1.0), amrex::Real(0.4), amrex::Real(0.2), amrex::Real(0.2)}},
        {amrex::Real(70.0), amrex::Real(65.0), amrex::Real(54.53), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.50),  //  9 Savannas
         {amrex::Real(0.3), amrex::Real(0.3), amrex::Real(0.5), amrex::Real(0.8), amrex::Real(1.8), amrex::Real(3.6), amrex::Real(3.8), amrex::Real(2.1), amrex::Real(0.9), amrex::Real(0.5), amrex::Real(0.3), amrex::Real(0.3)}},
        {amrex::Real(40.0), amrex::Real(100.0), amrex::Real(36.35), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.12),  // 10 Grasslands
         {amrex::Real(0.4), amrex::Real(0.5), amrex::Real(0.6), amrex::Real(0.7), amrex::Real(1.2), amrex::Real(3.0), amrex::Real(3.5), amrex::Real(1.5), amrex::Real(0.7), amrex::Real(0.6), amrex::Real(0.5), amrex::Real(0.4)}},
        {amrex::Real(70.0), amrex::Real(65.0), amrex::Real(55.97), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.30),  // 11 Permanent wetlands
         {amrex::Real(0.2), amrex::Real(0.3), amrex::Real(0.3), amrex::Real(0.5), amrex::Real(1.5), amrex::Real(2.9), amrex::Real(3.5), amrex::Real(2.7), amrex::Real(1.2), amrex::Real(0.3), amrex::Real(0.3), amrex::Real(0.2)}},
        {amrex::Real(40.0), amrex::Real(100.0), amrex::Real(36.25), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.15),  // 12 Croplands
         {amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(1.0), amrex::Real(2.0), amrex::Real(3.0), amrex::Real(3.0), amrex::Real(1.5), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)}},
        {amrex::Real(200.0), amrex::Real(999.0), amrex::Real(999.0), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(1.00),  // 13 Urban and Built-Up
         {amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)}},
        {amrex::Real(40.0), amrex::Real(100.0), amrex::Real(36.25), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.14),  // 14 cropland/natural vegetation mosaic
         {amrex::Real(0.2), amrex::Real(0.3), amrex::Real(0.3), amrex::Real(0.4), amrex::Real(1.1), amrex::Real(2.5), amrex::Real(3.2), amrex::Real(2.2), amrex::Real(1.1), amrex::Real(0.3), amrex::Real(0.3), amrex::Real(0.2)}},
        {amrex::Real(999.0), amrex::Real(999.0), amrex::Real(999.0), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.00),  // 15 Snow and Ice
         {amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)}},
        {amrex::Real(999.0), amrex::Real(999.0), amrex::Real(999.0), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.00),  // 16 Barren or Sparsely Vegetated
         {amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)}},
        {amrex::Real(100.0), amrex::Real(30.0), amrex::Real(51.75), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.00),  // 17 Water
         {amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)}},
        {amrex::Real(150.0), amrex::Real(100.0), amrex::Real(42.00), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.30),  // 18 Wooded Tundra
         {amrex::Real(1.0), amrex::Real(1.0), amrex::Real(1.1), amrex::Real(1.3), amrex::Real(1.7), amrex::Real(2.1), amrex::Real(2.1), amrex::Real(1.8), amrex::Real(1.3), amrex::Real(1.1), amrex::Real(1.0), amrex::Real(1.0)}},
        {amrex::Real(150.0), amrex::Real(100.0), amrex::Real(42.00), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.20),  // 19 Mixed Tundra
         {amrex::Real(0.6), amrex::Real(0.6), amrex::Real(0.7), amrex::Real(0.8), amrex::Real(1.2), amrex::Real(1.8), amrex::Real(1.8), amrex::Real(1.3), amrex::Real(0.8), amrex::Real(0.7), amrex::Real(0.6), amrex::Real(0.6)}},
        {amrex::Real(200.0), amrex::Real(100.0), amrex::Real(42.00), amrex::Real(298.0), amrex::Real(5000.0), amrex::Real(0.03),  // 20 Barren Tundra
         {amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)}},
    };
    if (category < 1 || category > noahmp_num_vegetation_categories) { return nullptr; }
    return &table[category - 1];
}

/**
 * @brief Leaf area index of @p params on day of year @p day, interpolated between the
 * monthly values as Noah-MP does for its table-LAI options (PhenologyMainMod): month
 * index 12 * day / days_in_year, linear between the two nearest months, wrapping December
 * to January. In the southern hemisphere Noah-MP shifts the day by half a year first;
 * pass @p southern = true for that.
 */
[[nodiscard]] inline amrex::Real
noahmp_leaf_area_index (const NoahMPVegetationParams& params, amrex::Real day,
                        amrex::Real days_in_year, bool southern) noexcept
{
    amrex::Real d = day;
    if (southern) {
        d = d + amrex::Real(0.5) * days_in_year;
        while (d >= days_in_year) { d -= days_in_year; }
    }
    const amrex::Real month = amrex::Real(12.0) * d / days_in_year;
    // Noah-MP: IntpMonth1 = MonthCurrent + 0.5 (an integer, so truncated),
    // IntpMonth2 = IntpMonth1 + 1, weight of month 1 = (IntpMonth1 + 0.5) - MonthCurrent.
    int m1 = static_cast<int>(month + amrex::Real(0.5));
    int m2 = m1 + 1;
    const amrex::Real w1 = (static_cast<amrex::Real>(m1) + amrex::Real(0.5)) - month;
    const amrex::Real w2 = amrex::Real(1.0) - w1;
    if (m1 < 1)  { m1 = 12; }
    if (m2 > 12) { m2 = 1; }
    return w1 * params.lai[m1 - 1] + w2 * params.lai[m2 - 1];
}

/**
 * @brief The day of the year, as Noah-MP counts it, of a date "YYYY-MM-DD[ HH:MM[:SS]]"
 * (ERF's start_datetime).
 *
 * The day is Noah-MP's DayJulianInYear, which PhenologyMainMod takes as
 * 0 <= day < days_in_year. It is 0 at 00:00 on 1 January, and the time of day is added
 * as a fraction. @p days_in_year is 365 or 366.
 *
 * @return false when the date cannot be read, the month is outside 1-12, the day is
 * outside the month, or the time is outside 00:00:00-23:59:60.
 */
[[nodiscard]] inline bool
noahmp_day_of_year (const std::string& datetime, amrex::Real& day,
                    amrex::Real& days_in_year)
{
    int year = 0, month = 0, mday = 0, hour = 0, minute = 0;
    double second = 0.0;
    const int n = std::sscanf(datetime.c_str(), "%d-%d-%d %d:%d:%lf",
                              &year, &month, &mday, &hour, &minute, &second);
    if (n < 3 || month < 1 || month > 12) { return false; }
    if (n < 4) { hour = 0; }
    if (n < 5) { minute = 0; }
    if (n < 6) { second = 0.0; }
    const bool leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
    const int month_days[12] = {31, leap ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (mday < 1 || mday > month_days[month - 1] || hour < 0 || hour > 23 ||
        minute < 0 || minute > 59 || second < 0.0 || second > 60.0) {
        return false;
    }
    int whole = mday - 1;
    for (int m = 0; m < month - 1; ++m) { whole += month_days[m]; }
    day = static_cast<amrex::Real>(whole) +
          static_cast<amrex::Real>((hour * 3600.0 + minute * 60.0 + second) / 86400.0);
    days_in_year = leap ? amrex::Real(366) : amrex::Real(365);
    return true;
}

#endif

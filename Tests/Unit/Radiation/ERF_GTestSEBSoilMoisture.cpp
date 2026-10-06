#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include <ERF_NoahMPSoilTable.H>
#include <ERF_SimplifiedSEB.H>
#include <ERF_SurfaceMoisture.H>
#include <ERF_NoahMPVegetationTable.H>

// The pieces of erf.radiation.seb_surface_layer_uses_moisture: the soil-water factor of
// the force-restore surface, the surface mixing ratio the surface layer takes from beta or
// from the two source fluxes (with the bare soil's pore-air humidity),
// Noah-MP's soil and vegetation tables (wilting point, field capacity, Jarvis parameters,
// monthly leaf area index), the Jarvis canopy and Sakaguchi-Zeng soil resistances, the
// aerodynamic resistance, the vapour-deficit factor,
// and the surface heat capacity derived from the soil.

namespace {
using amrex::Real;
constexpr Real tol = sizeof(Real) == 8 ? Real(1.0e-14) : Real(1.0e-6);
}

TEST(SEBSoilMoisture, AvailabilityIsLinearBetweenWiltingPointAndFieldCapacity)
{
    const Real wilt = Real(0.12);
    const Real fc = Real(0.387);
    EXPECT_EQ(seb_moisture_availability(Real(0.05), wilt, fc), Real(0.0));
    EXPECT_EQ(seb_moisture_availability(wilt, wilt, fc), Real(0.0));
    EXPECT_NEAR(seb_moisture_availability(Real(0.25), wilt, fc),
                (Real(0.25) - wilt) / (fc - wilt), tol);
    EXPECT_EQ(seb_moisture_availability(fc, wilt, fc), Real(1.0));
    EXPECT_EQ(seb_moisture_availability(Real(0.45), wilt, fc), Real(1.0));
}

TEST(SEBSoilMoisture, AvailabilityIsZeroForInvalidInputs)
{
    // No range between wilting point and field capacity, or a non-finite water content:
    // the surface does not evaporate rather than evaporating at an undefined rate.
    EXPECT_EQ(seb_moisture_availability(Real(0.3), Real(0.2), Real(0.2)), Real(0.0));
    EXPECT_EQ(seb_moisture_availability(Real(0.3), Real(0.3), Real(0.2)), Real(0.0));
    EXPECT_EQ(seb_moisture_availability(std::numeric_limits<Real>::quiet_NaN(),
                                        Real(0.1), Real(0.3)), Real(0.0));
}

TEST(SEBSoilMoisture, SurfaceMixingRatioBlendsSaturationAndAir)
{
    const Real q_sat = Real(0.030);
    const Real q_air = Real(0.010);
    // beta = 1: saturated surface; beta = 0: the air's own, so no moisture flux.
    EXPECT_EQ(erf_surface_moisture::surface_mixing_ratio(Real(1.0), q_sat, q_air), q_sat);
    EXPECT_EQ(erf_surface_moisture::surface_mixing_ratio(Real(0.0), q_sat, q_air), q_air);
    // The flux, proportional to q_surface - q_air, is beta times the potential one.
    const Real beta = Real(0.4);
    const Real q_surface = erf_surface_moisture::surface_mixing_ratio(beta, q_sat, q_air);
    EXPECT_NEAR(q_surface - q_air, beta * (q_sat - q_air), tol);
    // Out-of-range beta is clamped.
    EXPECT_EQ(erf_surface_moisture::surface_mixing_ratio(Real(1.5), q_sat, q_air), q_sat);
    EXPECT_EQ(erf_surface_moisture::surface_mixing_ratio(Real(-0.5), q_sat, q_air), q_air);
}

TEST(SEBSoilMoisture, NoahMPSoilTableLookup)
{
    // Silty clay loam (STAS 8), as NoahmpTable.TBL has it.
    const NoahMPSoilParams* silty_clay_loam = noahmp_soil_params(8);
    ASSERT_NE(silty_clay_loam, nullptr);
    EXPECT_EQ(silty_clay_loam->smc_wilt, Real(0.120));
    EXPECT_EQ(silty_clay_loam->smc_ref, Real(0.387));
    EXPECT_EQ(silty_clay_loam->smc_max, Real(0.464));
    // Categories are numbered from 1, and there are 19.
    EXPECT_NE(noahmp_soil_params(1), nullptr);
    EXPECT_NE(noahmp_soil_params(19), nullptr);
    EXPECT_EQ(noahmp_soil_params(0), nullptr);
    EXPECT_EQ(noahmp_soil_params(20), nullptr);
    // Every land category has a field capacity above its wilting point; water does not.
    for (int category = 1; category <= noahmp_num_soil_categories; ++category) {
        const NoahMPSoilParams* soil = noahmp_soil_params(category);
        if (category == 14) {
            EXPECT_FALSE(soil->smc_ref > soil->smc_wilt);
        } else {
            EXPECT_GT(soil->smc_ref, soil->smc_wilt) << "category " << category;
        }
    }
}

TEST(SEBSoilMoisture, CanopyResistanceFollowsJarvis)
{
    // Noah-MP's grassland (MODIS 10): RS 40 s/m, RGL 100 W/m^2, TOPT 298 K, RSMAX 5000 s/m.
    const NoahMPVegetationParams* grass = noahmp_vegetation_params(10);
    ASSERT_NE(grass, nullptr);
    const Real lai = Real(2.0);
    const Real sw = Real(800.0);
    const Real r = seb_canopy_resistance_without_vpd(grass->rs_min, grass->rs_max, grass->rgl,
                                                     grass->t_opt, lai, sw, grass->t_opt, Real(1.0));
    // At TOPT and a wet soil only the radiation factor acts.
    const Real f = Real(0.55) * Real(2.0) * sw / (grass->rgl * lai);
    const Real f_sw = (f + grass->rs_min / grass->rs_max) / (Real(1.0) + f);
    EXPECT_NEAR(r, grass->rs_min / (lai * f_sw), tol * r);
    // Darker, hotter (or colder) and drier all raise it.
    EXPECT_GT(seb_canopy_resistance_without_vpd(grass->rs_min, grass->rs_max, grass->rgl,
                                                grass->t_opt, lai, Real(100.0), grass->t_opt, Real(1.0)), r);
    EXPECT_GT(seb_canopy_resistance_without_vpd(grass->rs_min, grass->rs_max, grass->rgl,
                                                grass->t_opt, lai, sw, grass->t_opt + Real(10.0), Real(1.0)), r);
    EXPECT_GT(seb_canopy_resistance_without_vpd(grass->rs_min, grass->rs_max, grass->rgl,
                                                grass->t_opt, lai, sw, grass->t_opt, Real(0.5)), r);
    // No leaves: no transpiration.
    EXPECT_EQ(seb_canopy_resistance_without_vpd(grass->rs_min, grass->rs_max, grass->rgl,
                                                grass->t_opt, Real(0.0), sw, grass->t_opt, Real(1.0)),
              Real(1.0e6));
}

TEST(SEBSoilMoisture, SoilResistanceFollowsSakaguchiZeng)
{
    // Silty clay loam at 0.25 m^3/m^3 in a 0.1 m top layer, by Noah-MP's formula.
    const NoahMPSoilParams* soil = noahmp_soil_params(8);
    const Real theta = Real(0.25);
    const Real r = seb_soil_evaporation_resistance(theta, soil->smc_max, soil->smc_wilt, soil->bb,
                                                   Real(0.1), noahmp_soil_resistance_exponent);
    const double d_dry = 0.1 * (std::exp(std::pow(1.0 - 0.25 / 0.464, 5.0)) - 1.0) / (2.71828 - 1.0);
    const double diff = 2.2e-5 * 0.464 * 0.464 * std::pow(1.0 - 0.120 / 0.464, 2.0 + 3.0 / 8.72);
    EXPECT_NEAR(r, Real(d_dry / diff), Real(1.0e-6) * r);
    // A drier top layer resists more; a saturated one not at all.
    EXPECT_GT(seb_soil_evaporation_resistance(Real(0.15), soil->smc_max, soil->smc_wilt, soil->bb,
                                              Real(0.1), noahmp_soil_resistance_exponent), r);
    EXPECT_EQ(seb_soil_evaporation_resistance(soil->smc_max, soil->smc_max, soil->smc_wilt, soil->bb,
                                              Real(0.1), noahmp_soil_resistance_exponent), Real(0.0));
    EXPECT_EQ(seb_soil_evaporation_resistance(Real(0.005), soil->smc_max, soil->smc_wilt, soil->bb,
                                              Real(0.1), noahmp_soil_resistance_exponent), Real(1.0e6));
}

TEST(SEBSoilMoisture, TwoSourceSurfaceMixingRatioAndVapourDeficit)
{
    using namespace erf_surface_moisture;
    const Real q_sat = Real(0.020);
    const Real q_air = Real(0.008);
    // No surface resistance and moist pores: the saturated surface.
    EXPECT_NEAR(two_source_surface_mixing_ratio(Real(50.0), Real(0.0), Real(0.0), Real(0.7),
                                                q_sat, Real(1.0), q_air), q_sat, tol);
    // A canopy resistance four times the aerodynamic one gives 1/5 of the potential flux.
    EXPECT_NEAR(two_source_surface_mixing_ratio(Real(50.0), Real(200.0), Real(1.0e6), Real(1.0),
                                                q_sat, Real(1.0), q_air),
                q_air + Real(0.2) * (q_sat - q_air), Real(1.0e-6));
    // With moist pores (RH = 1) it is beta q_sat + (1 - beta) q_air, beta weighted by f_veg.
    const Real beta = Real(0.8) * Real(0.2) + Real(0.2) * Real(0.1);
    EXPECT_NEAR(two_source_surface_mixing_ratio(Real(50.0), Real(200.0), Real(450.0), Real(0.8),
                                                q_sat, Real(1.0), q_air),
                surface_mixing_ratio(beta, q_sat, q_air), tol);
    // The bare part evaporates from RH q_sat: with RH q_sat below the air's, the surface
    // mixing ratio falls below q_air (the soil takes up vapour), as in Noah-MP.
    EXPECT_LT(two_source_surface_mixing_ratio(Real(50.0), Real(1.0e6), Real(450.0), Real(0.0),
                                              q_sat, Real(0.2), q_air), q_air);
    EXPECT_NEAR(two_source_surface_mixing_ratio(Real(50.0), Real(1.0e6), Real(450.0), Real(0.0),
                                                q_sat, Real(0.2), q_air),
                q_air + Real(0.1) * (Real(0.2) * q_sat - q_air), tol);
    // The vapour-deficit factor: 1 with no deficit, smaller with one, floored at 0.01.
    EXPECT_EQ(vapour_deficit_factor(Real(36.35), Real(0.010), Real(0.012)), Real(1.0));
    EXPECT_NEAR(vapour_deficit_factor(Real(36.35), Real(0.020), Real(0.010)),
                Real(1.0) / (Real(1.0) + Real(0.3635)), tol);
    EXPECT_EQ(vapour_deficit_factor(Real(1.0e6), Real(0.020), Real(0.010)), Real(0.01));
    // The canopy resistance divides by that factor and is capped after the divide, so a
    // resistance below the cap cannot exceed it once the factor is applied.
    EXPECT_NEAR(canopy_resistance(Real(200.0), Real(36.35), Real(0.020), Real(0.010)),
                Real(200.0) * (Real(1.0) + Real(0.3635)), tol * Real(1000.0));
    EXPECT_EQ(canopy_resistance(Real(2.0e4), Real(1.0e6), Real(0.020), Real(0.010)),
              no_flux_resistance);
    // Neutral aerodynamic resistance: ln(z/z0) / (kappa u*).
    EXPECT_NEAR(aerodynamic_resistance(Real(10.0), Real(0.1), Real(0.0), Real(0.41), Real(0.3)),
                std::log(Real(100.0)) / (Real(0.41) * Real(0.3)), tol * Real(100.0));
}

TEST(SEBSoilMoisture, LeafAreaIndexInterpolatesTheMonthsAsNoahMP)
{
    const NoahMPVegetationParams* grass = noahmp_vegetation_params(10);
    // Day 218 of 366 (6 August of a leap year, counted from 0): month index 7.148, between
    // July (3.5) and August (1.5).
    const Real month = Real(12.0) * Real(218.0) / Real(366.0);
    const Real w_jul = Real(7.5) - month;
    EXPECT_NEAR(noahmp_leaf_area_index(*grass, Real(218.0), Real(366.0), false),
                w_jul * grass->lai[6] + (Real(1.0) - w_jul) * grass->lai[7], tol * Real(10.0));
    // Mid-month gives that month's value; early January wraps from December.
    EXPECT_NEAR(noahmp_leaf_area_index(*grass, Real(365.0) / Real(24.0), Real(365.0), false),
                grass->lai[0], tol * Real(10.0));
    const Real jan1 = noahmp_leaf_area_index(*grass, Real(1.0), Real(365.0), false);
    EXPECT_GT(jan1, std::min(grass->lai[11], grass->lai[0]) - tol);
    EXPECT_LT(jan1, std::max(grass->lai[11], grass->lai[0]) + tol);
}

TEST(SEBSoilMoisture, ForceRestoreHeatCapacityFromNoahMPSoil)
{
    // Silty clay loam at 0.25 m^3/m^3, by Noah-MP's formulas, worked by hand.
    const NoahMPSoilParams* soil = noahmp_soil_params(8);
    const Real c = seb_soil_heat_capacity(Real(0.25), soil->smc_max, noahmp_soil_heat_capacity);
    EXPECT_NEAR(c, Real(0.25 * 4.188e6 + (1.0 - 0.464) * 2.0e6 + (0.464 - 0.25) * 1004.64),
                Real(1.0e-6) * c);
    const double lambda_solid = std::pow(7.7, 0.1) * std::pow(2.0, 0.9);
    const double lambda_sat = std::pow(lambda_solid, 1.0 - 0.464) * std::pow(0.57, 0.464);
    const double gamma = (1.0 - 0.464) * 2700.0;
    const double lambda_dry = (0.135 * gamma + 64.7) / (2700.0 - 0.947 * gamma);
    const double kersten = std::log10(0.25 / 0.464) + 1.0;
    const Real lambda = seb_soil_thermal_conductivity(Real(0.25), soil->smc_max, soil->quartz);
    EXPECT_NEAR(lambda, Real(kersten * (lambda_sat - lambda_dry) + lambda_dry), Real(1.0e-5) * lambda);
    // The force-restore capacity: sqrt(lambda c tau / pi) / 2, about 1.2e5 J/m^2/K here.
    const Real c_s = seb_force_restore_heat_capacity(lambda, c, Real(86400.0));
    EXPECT_NEAR(c_s, Real(0.5 * std::sqrt(double(lambda) * double(c) * 86400.0 / 3.141592653589793)),
                Real(1.0e-5) * c_s);
    EXPECT_GT(c_s, Real(1.0e5));
    EXPECT_LT(c_s, Real(1.3e5));
    // A wetter soil holds more heat and conducts better.
    EXPECT_GT(seb_force_restore_heat_capacity(
                  seb_soil_thermal_conductivity(Real(0.4), soil->smc_max, soil->quartz),
                  seb_soil_heat_capacity(Real(0.4), soil->smc_max, noahmp_soil_heat_capacity),
                  Real(86400.0)), c_s);
}

TEST(SEBSoilMoisture, AerodynamicResistanceIsTheSurfaceLayersLimitedCmPsih)
{
    using namespace erf_surface_moisture;
    const Real kappa = Real(0.41);
    const Real z_ref = Real(10.0);
    const Real z0 = Real(0.1);
    const Real c = std::log(z_ref / z0);
    // ln(z/z0) - psi_h below 1 is limited to 1, as the surface_temp kernel's CmPsih.
    EXPECT_NEAR(aerodynamic_resistance(z_ref, z0, Real(5.0), kappa, Real(0.3)),
                Real(1.0) / (kappa * Real(0.3)), Real(1.0e-5) * Real(1.0) / (kappa * Real(0.3)));

    // Unstable air (z/L = -5): Jimenez's psi_h2, the kernel's, not Businger-Dyer's psi_h.
    const similarity_funs sfuns{};
    const Real unset = Real(1.0e34);
    const Real olen = Real(-2.0);
    const Real r = surface_layer_aerodynamic_resistance(sfuns, z_ref, z0, Real(0.3), olen,
                                                        Real(5.0), Real(0.1), kappa, unset);
    const Real r_jimenez = amrex::max(c - sfuns.calc_psi_h2(z_ref / olen), Real(1.0)) /
                           (kappa * Real(0.3));
    const Real r_businger = amrex::max(c - sfuns.calc_psi_h(z_ref / olen), Real(1.0)) /
                            (kappa * Real(0.3));
    EXPECT_NEAR(r, r_jimenez, Real(1.0e-5) * r_jimenez);
    // The two forms differ by more than 5 % here, so the check above tells them apart.
    EXPECT_GT(std::abs(r_jimenez - r_businger), Real(0.05) * r_jimenez);

    // Before the first flux (u* unset): neutral, u* = kappa U / max(ln(z/z0), 1).
    const Real u_neutral = kappa * Real(5.0) / c;
    EXPECT_NEAR(surface_layer_aerodynamic_resistance(sfuns, z_ref, z0, unset, unset,
                                                     Real(5.0), Real(0.1), kappa, unset),
                c / (kappa * u_neutral), Real(1.0e-5) * c / (kappa * u_neutral));
    // An unset Obukhov length with a valid u*: neutral psi.
    EXPECT_NEAR(surface_layer_aerodynamic_resistance(sfuns, z_ref, z0, Real(0.3), unset,
                                                     Real(5.0), Real(0.1), kappa, unset),
                c / (kappa * Real(0.3)), Real(1.0e-5) * c / (kappa * Real(0.3)));
}

TEST(SEBSoilMoisture, DayOfYearCountsFromZeroAsNoahMP)
{
    Real day = Real(-1.0);
    Real days = Real(0.0);
    // PhenologyMainMod: 0 <= day < days in the year, 0 at 00:00 on 1 January.
    ASSERT_TRUE(noahmp_day_of_year("2024-01-01 00:00:00", day, days));
    EXPECT_EQ(day, Real(0.0));
    EXPECT_EQ(days, Real(366.0));
    ASSERT_TRUE(noahmp_day_of_year("2024-03-01", day, days));
    EXPECT_EQ(day, Real(60.0));
    ASSERT_TRUE(noahmp_day_of_year("2023-12-31 12:00:00", day, days));
    EXPECT_NEAR(day, Real(364.5), tol * Real(1000.0));
    EXPECT_EQ(days, Real(365.0));
    // The start of the comparison case: 2024-08-05 15:00 UTC.
    ASSERT_TRUE(noahmp_day_of_year("2024-08-05 15:00:00", day, days));
    EXPECT_NEAR(day, Real(217.625), tol * Real(1000.0));
    // Grassland LAI there: month index 12 * 217.625 / 366, between July and August.
    const NoahMPVegetationParams* grass = noahmp_vegetation_params(10);
    const Real w_jul = Real(7.5) - Real(12.0) * Real(217.625) / Real(366.0);
    EXPECT_NEAR(noahmp_leaf_area_index(*grass, day, days, false),
                w_jul * grass->lai[6] + (Real(1.0) - w_jul) * grass->lai[7], tol * Real(1000.0));
    // Dates that are not dates.
    EXPECT_FALSE(noahmp_day_of_year("2023-02-29", day, days));
    EXPECT_FALSE(noahmp_day_of_year("2024-13-01", day, days));
    EXPECT_FALSE(noahmp_day_of_year("2024-04-31 00:00:00", day, days));
    EXPECT_FALSE(noahmp_day_of_year("2024-08-05 24:00:00", day, days));
    EXPECT_FALSE(noahmp_day_of_year("August 5", day, days));
}

TEST(SEBSoilMoisture, UnvegetatedCategoriesHaveNoRoughness)
{
    // Snow and ice, barren and water: Z0MVT = 0 and no leaf area all year, which is why
    // erf.radiation.seb_vegetation_type stops at start-up on them.
    for (int category : {15, 16, 17}) {
        const NoahMPVegetationParams* veg = noahmp_vegetation_params(category);
        ASSERT_NE(veg, nullptr);
        EXPECT_EQ(veg->z0, Real(0.0)) << "category " << category;
        for (int m = 0; m < 12; ++m) { EXPECT_EQ(veg->lai[m], Real(0.0)) << "category " << category; }
    }
    // Every other category has a positive roughness.
    for (int category = 1; category <= 20; ++category) {
        if (category >= 15 && category <= 17) { continue; }
        EXPECT_GT(noahmp_vegetation_params(category)->z0, Real(0.0)) << "category " << category;
    }
}

TEST(SEBSoilMoisture, BareSoilAtTheWiltingPointDoesNotEvaporate)
{
    using namespace erf_surface_moisture;
    // Noah-MP's pore-air humidity exp(psi g / (R_v T)), psi = -psi_sat (theta/theta_sat)^-b.
    const NoahMPSoilParams* soil = noahmp_soil_params(8);
    const Real t = Real(325.0);
    const double psi = -double(soil->psi_sat) *
                       std::pow(double(soil->smc_wilt) / double(soil->smc_max), -double(soil->bb));
    const Real rh_wilt = seb_soil_surface_relative_humidity(soil->smc_wilt, soil->smc_max,
                                                            soil->psi_sat, soil->bb, t);
    EXPECT_NEAR(rh_wilt, Real(std::exp(psi * 9.80616 / (461.269 * 325.0))), Real(1.0e-4) * rh_wilt);
    EXPECT_GT(rh_wilt, Real(0.003));   // about 0.005 for silty clay loam at 325 K
    EXPECT_LT(rh_wilt, Real(0.007));
    // Moist soil: the pores are nearly saturated.
    EXPECT_GT(seb_soil_surface_relative_humidity(Real(0.25), soil->smc_max, soil->psi_sat,
                                                 soil->bb, t), Real(0.98));
    // Invalid inputs give no reduction.
    EXPECT_EQ(seb_soil_surface_relative_humidity(Real(0.25), soil->smc_max, soil->psi_sat,
                                                 soil->bb, Real(0.0)), Real(1.0));
    // At the wilting point on a hot afternoon (q_sat(325 K) about 0.1 kg/kg, air 8 g/kg),
    // bare soil does not evaporate: q_surf <= q_air, so LE <= 0, as Noah-MP gives.
    const Real r_soil = seb_soil_evaporation_resistance(soil->smc_wilt, soil->smc_max, soil->smc_wilt,
                                                        soil->bb, Real(0.1),
                                                        noahmp_soil_resistance_exponent);
    EXPECT_LE(two_source_surface_mixing_ratio(Real(40.0), Real(1.0e6), r_soil, Real(0.0),
                                              Real(0.1), rh_wilt, Real(0.008)), Real(0.008));
}

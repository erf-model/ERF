#include <gtest/gtest.h>

#include "Diagnostics/ERF_SurfaceFluxDiagnostics.H"

using namespace surface_flux_diagnostics;

// Motivation: The diagnostic layer must convert already-conservative surface
// fluxes directly to W m^-2 and must not introduce a density factor.
TEST(SurfaceFluxDiagnostics, SensibleHeatFluxConversionUsesCpD)
{
    const amrex::Real positive = amrex::Real(1.75);
    const amrex::Real negative = amrex::Real(-0.5);

    EXPECT_DOUBLE_EQ(sensible_heat_flux_wm2_from_rhotheta_flux(positive), Cp_d * positive);
    EXPECT_DOUBLE_EQ(sensible_heat_flux_wm2_from_rhotheta_flux(negative), Cp_d * negative);
}

// Motivation: The diagnostic layer must convert already-conservative surface
// fluxes directly to W m^-2 and must not introduce a density factor.
TEST(SurfaceFluxDiagnostics, LatentHeatFluxConversionUsesLv)
{
    const amrex::Real positive = amrex::Real(1.75);
    const amrex::Real negative = amrex::Real(-0.5);

    EXPECT_DOUBLE_EQ(latent_heat_flux_wm2_from_rhoqv_flux(positive), L_v * positive);
    EXPECT_DOUBLE_EQ(latent_heat_flux_wm2_from_rhoqv_flux(negative), L_v * negative);
}

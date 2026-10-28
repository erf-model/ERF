#include <gtest/gtest.h>

#include <string>
#include <utility>
#include <vector>

#include <ERF_IndexDefines.H>
#include <ERF_DataStruct.H>
#include <ERF_Utils.H>
#include <ERF_ReadBndryPlanes.H>
#include <ERF_NativeScalarDiffusion.H>

/**
 * @file ERF_GTestMultiScalarSlots.cpp
 * @brief Component maps with more than one passive scalar (NSCALARS > 1)
 *
 * Builds with extra passive scalars (the fire coupling's smoke, ERF_USE_FIRE)
 * put them between rhoadv_0 and rhoQ1, so every moisture component moves up by
 * NSCALARS - 1 while the boundary-condition slots (BCVars) and the boundary
 * planes keep one slot for all the scalars. The checks cover the three maps
 * that depend on that layout: the WRF real-boundary maps (which WRF field fills
 * which component, preserving the upstream hydrometeor read flags),
 * the boundary-plane reader (which plane slot a field is stored in), and the
 * implicit diffusion's molecular coefficients (one for every scalar). With
 * NSCALARS = 1 the slot maps are the identity on these components, so the slot
 * checks bite in builds with extra scalars.
 */

TEST(MultiScalarSlots, EveryConservedComponentHasItsBoundarySlot)
{
    EXPECT_EQ(cons_bc_comp(Rho_comp),      BCVars::Rho_bc_comp);
    EXPECT_EQ(cons_bc_comp(RhoTheta_comp), BCVars::RhoTheta_bc_comp);
    EXPECT_EQ(cons_bc_comp(RhoKE_comp),    BCVars::RhoKE_bc_comp);
    for (int n = 0; n < NSCALARS; ++n) {
        EXPECT_EQ(cons_bc_comp(RhoScalar_comp + n), BCVars::RhoScalar_bc_comp) << "scalar " << n;
    }
    const std::vector<std::pair<int, int>> moist = {
        {RhoQ1_comp, BCVars::RhoQ1_bc_comp}, {RhoQ2_comp, BCVars::RhoQ2_bc_comp},
        {RhoQ3_comp, BCVars::RhoQ3_bc_comp}, {RhoQ6_comp, BCVars::RhoQ6_bc_comp},
        {RhoQ11_comp, BCVars::RhoQ11_bc_comp}};
    for (const auto& [comp, slot] : moist) {
        EXPECT_EQ(cons_bc_comp(comp), slot) << "conserved component " << comp;
    }
}

TEST(MultiScalarSlots, ThePlaneReaderTakesEachVariableFromTheSlotItWasStoredIn)
{
    // the conserved component each plane variable fills
    const std::vector<std::pair<std::string, int>> filled = {
        {"density", Rho_comp}, {"temperature", RhoTheta_comp}, {"theta", RhoTheta_comp},
        {"ke", RhoKE_comp}, {"scalar", RhoScalar_comp}, {"qv", RhoQ1_comp}, {"qc", RhoQ2_comp}};
    for (const auto& [name, comp] : filled) {
        EXPECT_EQ(cons_bc_comp(comp), ReadBndryPlanes::plane_slot(name)) << name;
    }
    EXPECT_EQ(ReadBndryPlanes::plane_slot("velocity"), static_cast<int>(BCVars::xvel_bc));
    EXPECT_EQ(ReadBndryPlanes::plane_slot("pressure"), -1);
}

TEST(MultiScalarSlots, TheWrfBoundaryMapsReadNoScalarAndPutTheVapourInQ1)
{
    MoistureComponentIndices mi;
    mi.qv = RhoQ1_comp;
    mi.qc = RhoQ2_comp;
    for (bool read_density : {false, true}) {
        const RealBdyConsMaps m = real_bdy_cons_maps(read_density, mi, false, false);
        ASSERT_GE(static_cast<int>(m.read.size()), NVAR_max);
        EXPECT_EQ(m.read[Rho_comp], read_density ? 1 : 0);
        EXPECT_EQ(m.map[Rho_comp], read_density ? static_cast<int>(RealBdyVars::R) : Rho_comp);
        EXPECT_EQ(m.read[RhoTheta_comp], 1);
        EXPECT_EQ(m.map[RhoTheta_comp], static_cast<int>(RealBdyVars::T));
        EXPECT_EQ(m.read[RhoKE_comp], 0);
        for (int n = 0; n < NSCALARS; ++n) {
            EXPECT_EQ(m.read[RhoScalar_comp + n], 0) << "scalar " << n << " is not a WRF field";
        }
        EXPECT_EQ(m.read[RhoQ1_comp], 1);
        EXPECT_EQ(m.map[RhoQ1_comp], static_cast<int>(RealBdyVars::QV));
        EXPECT_EQ(m.read[RhoQ2_comp], 0);
        EXPECT_EQ(m.bnd[RhoQ1_comp], static_cast<int>(BCVars::RhoQ1_bc_comp));
        EXPECT_EQ(m.bnd[RhoQ2_comp], static_cast<int>(BCVars::RhoQ2_bc_comp));
        EXPECT_EQ(m.bnd[RhoScalar_comp + NSCALARS - 1], static_cast<int>(BCVars::RhoScalar_bc_comp));
    }
    // Hydrometeor field mappings do not change the upstream read flags.
    const RealBdyConsMaps q = real_bdy_cons_maps(false, mi, true, false);
    EXPECT_EQ(q.read[RhoQ1_comp], 1);
    EXPECT_EQ(q.map[RhoQ1_comp], static_cast<int>(RealBdyVars::QV));
    EXPECT_EQ(q.read[RhoQ2_comp], 0);
    EXPECT_EQ(q.map[RhoQ2_comp], static_cast<int>(RealBdyVars::QC));
    EXPECT_EQ(q.read[RhoQ3_comp], 0);
}

TEST(MultiScalarSlots, HydrometeorMappingsPreserveUpstreamReadFlags)
{
    const MoistureComponentIndices mi(RhoQ1_comp, RhoQ2_comp, RhoQ3_comp,
                                       RhoQ4_comp, RhoQ5_comp, RhoQ6_comp);
    for (bool read_density : {false, true}) {
        const RealBdyConsMaps legacy = real_bdy_cons_maps(read_density, mi, false, false);
        for (bool separate_hydrometeors : {false, true}) {
            const RealBdyConsMaps mapped = real_bdy_cons_maps(
                read_density, mi, true, separate_hydrometeors);
            EXPECT_EQ(mapped.read, legacy.read);
            EXPECT_EQ(mapped.bnd, legacy.bnd);
            for (int comp : {mi.qc, mi.qi, mi.qr, mi.qs, mi.qg}) {
                EXPECT_EQ(mapped.read[comp], 0) << "component " << comp;
            }
        }
    }
}

TEST(MultiScalarSlots, EveryScalarGetsTheMolecularDiffusivityOfTheExplicitPath)
{
    for (MolecDiffType type : {MolecDiffType::Constant, MolecDiffType::ConstantAlpha}) {
        DiffChoice choice;
        choice.molec_diff_type = type;
        choice.alpha_T = amrex::Real(2.0e-5);  choice.rhoAlpha_T = amrex::Real(2.4e-5);
        choice.alpha_C = amrex::Real(3.0e-5);  choice.rhoAlpha_C = amrex::Real(3.6e-5);
        const amrex::Vector<amrex::Real> alpha = VerticalMolecularAlphaByPrim(choice);
        for (int comp = RhoTheta_comp; comp < RhoQ1_comp + 6; ++comp) {
            if (comp == RhoKE_comp) { continue; }
            const auto policy = ResolveNativeScalarDiffusionPolicy(comp, choice);
            EXPECT_EQ(alpha[comp - 1], policy.coefficients.molecular_coeff) << "component " << comp;
            EXPECT_GT(alpha[comp - 1], amrex::Real(0)) << "component " << comp;
        }
    }
}

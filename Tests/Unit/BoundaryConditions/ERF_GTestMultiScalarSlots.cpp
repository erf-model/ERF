#include <gtest/gtest.h>

#include <string>
#include <utility>
#include <vector>

#include <ERF_IndexDefines.H>
#include <ERF_DataStruct.H>
#include <ERF_Utils.H>
#include <ERF_ReadBndryPlanes.H>
#include <ERF_NativeScalarDiffusion.H>
#include <ERF_PhysBCFunct.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>

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

namespace {

// NVCC requires device lambdas outside GoogleTest's non-public TestBody.
void implicit_coarse_fine_runs_test ()
{
    using namespace amrex;
    const Box coarse_domain(IntVect(0), IntVect(7));
    const Box fine_domain = refine(coarse_domain, 2);
    const RealBox rb(Real(0), Real(0), Real(0), Real(8), Real(8), Real(8));
    const Geometry cg(coarse_domain, rb, CoordSys::cartesian, {0, 0, 0});
    const Geometry fg(fine_domain, rb, CoordSys::cartesian, {0, 0, 0});
    Box fine_box = fine_domain;
    fine_box.setBig(2, 7);
    const BoxArray cba(coarse_domain), fba(fine_box);
    const DistributionMapping cdm(cba), fdm(fba);
    constexpr int ncomp = NVAR_max;
    MultiFab coarse(cba, cdm, ncomp, 2), fine(fba, fdm, ncomp, 2);
    MultiFab grouped(fba, fdm, ncomp, 2), reference(fba, fdm, ncomp, 2);
    coarse.setVal(Real(0));
    fine.setVal(Real(-17));
    grouped.setVal(Real(-17));
    reference.setVal(Real(-17));
    for (MFIter mfi(coarse); mfi.isValid(); ++mfi) {
        const auto c = coarse.array(mfi);
        ParallelFor(mfi.fabbox(), ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            c(i,j,k,n) = Real(10+n) + Real(i*i) + Real(j*j) + Real(k*k);
        });
    }
    Vector<BCRec> bcs(BCVars::NumTypes);
    for (int n = 0; n < BCVars::NumTypes; ++n) {
        for (int d = 0; d < 3; ++d) {
            bcs[n].setLo(d, n % 2 ? BCType::ext_dir : BCType::foextrap);
            bcs[n].setHi(d, n % 2 ? BCType::foextrap : BCType::ext_dir);
        }
    }
    const Vector<Real> times{Real(0)};
    Vector<MultiFab*> cmf{&coarse}, fmf{&fine};
    const IntVect ng(0,0,1);
    for (int start = RhoKE_comp; start < ncomp; ) {
        int end = start + 1;
        while (end < ncomp && cons_bc_run_contiguous(end-1, end)) { ++end; }
        FillPatchTwoLevels(grouped, ng, IntVect(0), Real(0), cmf, times, fmf, times,
                           start, start, end-start, cg, fg, IntVect(2),
                           &cell_cons_interp, bcs, cons_bc_comp(start));
        start = end;
    }
    // One fill per component is an independent oracle for the grouped BC slices.
    for (int n = RhoKE_comp; n < ncomp; ++n) {
        const int slot = n < RhoScalar_comp ? n :
                         n < RhoQ1_comp ? int(BCVars::RhoScalar_bc_comp) :
                         int(BCVars::RhoQ1_bc_comp) + n - RhoQ1_comp;
        FillPatchTwoLevels(reference, ng, IntVect(0), Real(0), cmf, times, fmf, times,
                           n, n, 1, cg, fg, IntVect(2), &cell_cons_interp, bcs, slot);
    }
    MultiFab::Subtract(grouped, reference, RhoKE_comp, RhoKE_comp,
                       ncomp-RhoKE_comp, ng);
    for (int n = RhoKE_comp; n < ncomp; ++n) {
        EXPECT_EQ(grouped.norm0(n, 1, ng), Real(0)) << "state component " << n;
    }
}

} // namespace

// Motivation: Conserved components and BC slots do not advance together when
// passive scalars share a slot. Grouped implicit ghost fills must match separate
// per-component fills, including appended species with distinct boundary data.
TEST(MultiScalarSlots, ImplicitCoarseFineRunsMatchIndependentBoundarySlices)
{
    implicit_coarse_fine_runs_test();
}

// Motivation: Prescribed water-vapor inflow must not become an unintended source
// of appended species. Their unsupported prescribed inflow remains zero without
// overwriting valid interior vapor.
TEST(MultiScalarSlots, AppendedInflowDoesNotInheritWaterVaporExternalValue)
{
    using namespace amrex;
    const Box domain(IntVect(0), IntVect(3));
    const RealBox rb(Real(0), Real(0), Real(0), Real(4), Real(4), Real(4));
    const Geometry geom(domain, rb, CoordSys::cartesian, {0, 1, 0});
    Vector<BCRec> bcs(BCVars::NumTypes);
    for (auto& bc : bcs) {
        for (int d = 0; d < 3; ++d) {
            bc.setLo(d, ERFBCType::foextrap);
            bc.setHi(d, ERFBCType::foextrap);
        }
        bc.setLo(0, ERFBCType::ext_dir);
    }
    Gpu::DeviceVector<BCRec> device_bcs(bcs.size());
    Gpu::copy(Gpu::hostToDevice, bcs.begin(), bcs.end(), device_bcs.begin());
    Array<Array<Real,6>, AMREX_SPACEDIM+NBCVAR_max> external{}, gradients{};
    external[BCVars::RhoQ1_bc_comp][0] = Real(0.123);
    std::unique_ptr<MultiFab> z;
    ERFPhysBCFunct_cons bc(0, geom, bcs, device_bcs, external, gradients,
                           TerrainType::None, z, false, nullptr, {});
    FArrayBox state(grow(domain, 1), NVAR_max);
    FArrayBox u(surroundingNodes(grow(domain, 1), 0), 1);
    FArrayBox v(surroundingNodes(grow(domain, 1), 1), 1);
    state.setVal<RunOn::Device>(Real(7));
    u.setVal<RunOn::Device>(Real(1));
    v.setVal<RunOn::Device>(Real(0));
    bc.impose_lateral_cons_bcs(state.array(), u.const_array(), v.const_array(),
                               grow(domain, IntVect(1,1,0)), domain,
                               RhoQ1_comp, 10, IntVect(1), 0);
    Gpu::streamSynchronize();
    FArrayBox host(state.box(), NVAR_max, The_Pinned_Arena());
    Gpu::dtoh_memcpy(host.dataPtr(), state.dataPtr(), state.nBytes());
    EXPECT_EQ(host.const_array()(-1,1,1,RhoQ1_comp), Real(0.123));
    for (int n : {RhoQ7_comp, RhoQ9_comp}) {
        EXPECT_EQ(host.const_array()(-1,1,1,n), Real(0));
        EXPECT_EQ(host.const_array()(0,1,1,n), Real(7));
    }
}

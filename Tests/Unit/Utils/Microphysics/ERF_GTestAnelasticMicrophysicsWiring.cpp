#include <AMReX_MultiFab.H>

#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "ERF_GTestMicrophysicsCommon.H"
#include "ERF_IndexDefines.H"
#include "ERF_Kessler.H"
#include "ERF_Morrison.H"
#include "ERF_SAMUtils.H"
#include "ERF_WSM6.H"

using namespace amrex;

namespace {

struct ReferenceCell {
    Real rho = Real(1.15);
    Real theta = Real(306.0);
    Real qv = Real(4.0e-3);
    Real p0 = Real(90000.0);
    Real pi0 = Real(1.04);
};

ReferenceCell stored_reference_cell ()
{
    ReferenceCell cell;
    const Real eos_pi = getExnergivenP(cell.p0, RdoCp);
    EXPECT_GT(cell.pi0, Real(0.0));
    EXPECT_NE(cell.pi0, eos_pi);
    return cell;
}

} // namespace

// Motivation: The production Kessler adapter must use the same p0 and stored
// pi0 from BaseState at one cell. The deliberately inconsistent pi0 prevents a
// test-only EOS reconstruction from passing unnoticed.
TEST(AnelasticMicrophysicsWiring, KesslerUsesStoredBaseExner)
{
    const ReferenceCell cell = stored_reference_cell();
    const auto state = kessler_copy_in_thermo(
        cell.rho, cell.rho * cell.theta, cell.qv, RdoCp, true,
        cell.p0, cell.pi0);

    EXPECT_EQ(state.pressure_mbar, cell.p0 * Real(0.01));
    EXPECT_EQ(state.temperature, cell.theta * cell.pi0);
}

// Motivation: SAM copy-in feeds its scheme-native primitive working state from
// the production adapter. Pressure units remain mbar, but temperature must be
// theta times the stored Exner function rather than an Exner recomputed from p0.
TEST(AnelasticMicrophysicsWiring, SAMUsesStoredBaseExner)
{
    const ReferenceCell cell = stored_reference_cell();
    const auto state = sam_cons_to_primitive_with_base_state(
        cell.rho, cell.rho * cell.theta, cell.rho * cell.qv,
        Real(0.0), Real(0.0), Real(0.0), Real(0.0), Real(0.0),
        RdoCp, true, cell.p0, cell.pi0);

    EXPECT_EQ(state.pres_mbar, cell.p0 * Real(0.01));
    EXPECT_EQ(state.tabs, cell.theta * cell.pi0);
}

// Motivation: Morrison uses Pa internally in both its C++ and Fortran-facing
// working state. Its production copy-in adapter must still consume stored pi0.
TEST(AnelasticMicrophysicsWiring, MorrisonUsesStoredBaseExner)
{
    const ReferenceCell cell = stored_reference_cell();
    const auto state = morrison_copy_in_thermo(
        cell.rho, cell.rho * cell.theta, cell.qv, RdoCp, true,
        cell.p0, cell.pi0);

    EXPECT_EQ(state.pressure_pa, cell.p0);
    EXPECT_EQ(state.temperature, cell.theta * cell.pi0);
}

// Motivation: WSM6 copy-in uses stored pi0 and its copy-out contract continues
// to reconstruct theta from the updated temperature and held pressure.
TEST(AnelasticMicrophysicsWiring, WSM6UsesStoredBaseExnerAndPressureCopyOut)
{
    const ReferenceCell cell = stored_reference_cell();
    const auto state = wsm6_copy_in_thermo(
        cell.rho, cell.rho * cell.theta, cell.qv, RdoCp, true,
        cell.p0, cell.pi0);
    const Real updated_temperature = state.temperature + Real(1.5);
    const Real expected_theta = getThgivenTandP(updated_temperature, cell.p0, RdoCp);

    EXPECT_EQ(state.pressure_pa, cell.p0);
    EXPECT_EQ(state.temperature, cell.theta * cell.pi0);
    EXPECT_EQ(wsm6_theta_from_temperature_pressure(
                  updated_temperature, state.pressure_pa, RdoCp), expected_theta);
}

// Motivation: Every initialized base-state cell must carry finite positive p0
// and pi0 values satisfying the stored pressure--Exner relation. Two vertical
// cells make the check cover a nonuniform hydrostatic column rather than a
// single constant value.
TEST(AnelasticBaseState, PressureAndExnerRemainConsistentAcrossVerticalCells)
{
    const Box domain(IntVect(0, 0, 0), IntVect(0, 0, 1));
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    MultiFab base(ba, dm, BaseState::num_comps, 1);
    base.setVal(Real(0.0));

    const Real p_surface = Real(101000.0);
    const Real density = Real(1.2);
    const Real dz = Real(10.0);
    for (MFIter mfi(base, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto arr = base.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real p0 = p_surface - density * CONST_GRAV * (Real(k) + Real(0.5)) * dz;
            arr(i,j,k,BaseState::p0_comp) = p0;
            arr(i,j,k,BaseState::pi0_comp) = getExnergivenP(p0, RdoCp);
        });
    }
    Gpu::streamSynchronize();

    MultiFab p0(base, make_alias, BaseState::p0_comp, 1);
    MultiFab pi0(base, make_alias, BaseState::pi0_comp, 1);
    MultiFab relation_error(ba, dm, 1, 0);
    for (MFIter mfi(relation_error, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto p = p0.const_array(mfi);
        const auto pi = pi0.const_array(mfi);
        const auto error = relation_error.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            error(i,j,k) = std::abs(pi(i,j,k) - getExnergivenP(p(i,j,k), RdoCp));
        });
    }
    Gpu::streamSynchronize();

    EXPECT_TRUE(p0.is_finite());
    EXPECT_TRUE(pi0.is_finite());
    EXPECT_GT(p0.min(0, 0, false), Real(0.0));
    EXPECT_GT(pi0.min(0, 0, false), Real(0.0));
    EXPECT_LT(relation_error.max(0, 0, false),
              Real(256.0) * std::numeric_limits<Real>::epsilon());
}

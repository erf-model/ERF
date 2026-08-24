#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_CanopyBiophysics.H"

namespace {

struct CanopySourceResult
{
    amrex::Real theta_min;
    amrex::Real theta_max;
    amrex::Real qv_min;
    amrex::Real qv_max;
};

CanopySourceResult
run_canopy_source (MoistureType moisture_type)
{
    using namespace amrex;

    const Box domain(IntVect(0, 0, 0), IntVect(1, 0, 0));
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    const int conserved_components =
        (moisture_type == MoistureType::None) ? RhoTheta_comp + 1 : RhoQ1_comp + 1;

    MultiFab source(ba, dm, conserved_components, 0);
    MultiFab state(ba, dm, conserved_components, 0);
    MultiFab leaf_area_density(ba, dm, 1, 0);
    MultiFab base_state(ba, dm, BaseState::num_comps, 0);
    MultiFab xvel(convert(ba, IntVect(1, 0, 0)), dm, 1, 0);
    MultiFab yvel(convert(ba, IntVect(0, 1, 0)), dm, 1, 0);
    MultiFab zvel(convert(ba, IntVect(0, 0, 1)), dm, 1, 0);

    source.setVal(0.0);
    state.setVal(0.0);
    state.setVal(1.0, Rho_comp, 1);
    state.setVal(295.0, RhoTheta_comp, 1);
    if (moisture_type != MoistureType::None) {
        state.setVal(0.001, RhoQ1_comp, 1);
    }
    base_state.setVal(0.0);
    base_state.setVal(100000.0, BaseState::p0_comp, 1);
    base_state.setVal(1.0, BaseState::pi0_comp, 1);
    xvel.setVal(0.1);
    yvel.setVal(0.0);
    zvel.setVal(0.0);

    for (MFIter mfi(leaf_area_density, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const auto lad = leaf_area_density.array(mfi);
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            lad(i, j, k) = (i == domain.smallEnd(0)) ? Real(1.0) : Real(0.0);
        });
    }

    SolverChoice solver_choice;
    solver_choice.moisture_type = moisture_type;
    solver_choice.forest_biophysics_heat = true;
    solver_choice.forest_leaf_theta_fixed = 301.0;

    AddCanopyBiophysicsHeatSources(
        source, state, xvel, yvel, zvel, &leaf_area_density, base_state, solver_choice);
    Gpu::streamSynchronize();

    CanopySourceResult result{
        source.min(RhoTheta_comp), source.max(RhoTheta_comp), Real(0.0), Real(0.0)};
    if (moisture_type != MoistureType::None) {
        result.qv_min = source.min(RhoQ1_comp);
        result.qv_max = source.max(RhoQ1_comp);
    }
    return result;
}

} // namespace

// Motivation: transpiration must enter staged rho-qv storage so scalar
// advection cannot overwrite it. The adjacent zero-LAD cell is the
// independent non-canopy control, and sensible heat must remain active.
TEST(CanopyBiophysics, MoistCanopyAddsStagedHeatAndWaterOnlyWhereLADIsPositive)
{
    const auto result = run_canopy_source(MoistureType::MoistNoCondensation);

    EXPECT_EQ(result.theta_min, amrex::Real(0.0));
    EXPECT_GT(result.theta_max, amrex::Real(0.0));
    EXPECT_EQ(result.qv_min, amrex::Real(0.0));
    EXPECT_GT(result.qv_max, amrex::Real(0.0));
}

// Motivation: the dry path can allocate no rho-qv component; its runtime
// moisture guard must avoid that component while retaining sensible heat.
TEST(CanopyBiophysics, DryCanopyDoesNotRequireWaterVaporStorage)
{
    const auto result = run_canopy_source(MoistureType::None);

    EXPECT_EQ(result.theta_min, amrex::Real(0.0));
    EXPECT_GT(result.theta_max, amrex::Real(0.0));
}

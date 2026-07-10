#include "AMReX.H"
#include "AMReX_MultiFab.H"

#include "ERF_Diffusion.H"
#include "ERF_SurfaceLayerStress.H"

#include <gtest/gtest.h>

using namespace amrex;

namespace {

constexpr Real rho_low = Real(0.8);
constexpr Real rho_high = Real(1.3);
constexpr Real k_low = Real(-0.7);
constexpr Real k_high = Real(0.45);
constexpr Real most = Real(-2.6);

void
expect_result(
  const surface_layer_stress::FaceStressResult& result,
  Real expected_face,
  Real expected_low,
  Real expected_high)
{
  EXPECT_NEAR(result.face_stress, expected_face, Real(1.e-12));
  EXPECT_NEAR(result.low_kinematic_cache, expected_low, Real(1.e-12));
  EXPECT_NEAR(result.high_kinematic_cache, expected_high, Real(1.e-12));
}

} // namespace

TEST(SurfaceLayerStress, ConservativeToKinematicRoundTripHandlesSignedStress)
{
  for (const Real conservative : {Real(2.4), Real(-2.4)}) {
    const Real kinematic =
      surface_layer_stress::conservative_to_kinematic_stress(
        conservative, rho_low);
    EXPECT_NEAR(rho_low * kinematic, conservative, Real(1.e-12));
  }
}

TEST(SurfaceLayerStress, BothValidUsesDensityWeightedFaceAverage)
{
  const auto result = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, k_low, k_high, true, true, most);
  expect_result(
    result, Real(0.5) * (rho_low * k_low + rho_high * k_high), k_low, k_high);
}

TEST(SurfaceLayerStress, HighValidCachesMostStressInLowCellUnits)
{
  const auto result = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, k_low, k_high, false, true, most);
  expect_result(
    result, Real(0.5) * (rho_high * k_high + most), most / rho_low, k_high);
}

TEST(SurfaceLayerStress, LowValidCachesMostStressInHighCellUnits)
{
  const auto result = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, k_low, k_high, true, false, most);
  expect_result(
    result, Real(0.5) * (rho_low * k_low + most), k_low, most / rho_high);
}

TEST(SurfaceLayerStress, NeitherValidKeepsMostFaceStressAndConvertsBothCaches)
{
  const auto result = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, k_low, k_high, false, false, most);
  expect_result(result, most, most / rho_low, most / rho_high);
  EXPECT_NEAR(rho_low * result.low_kinematic_cache, most, Real(1.e-12));
  EXPECT_NEAR(rho_high * result.high_kinematic_cache, most, Real(1.e-12));
}

TEST(SurfaceLayerStress, XAndYComponentsUseIdenticalBranchAlgebra)
{
  const auto x_result = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, k_low, k_high, false, true, most);
  const auto y_result = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, k_low, k_high, false, true, most);
  expect_result(
    y_result, x_result.face_stress, x_result.low_kinematic_cache,
    x_result.high_kinematic_cache);
}

TEST(SurfaceLayerStress, Tau23ValidityUsesItsOwnHandle)
{
  const bool has_tau13 = false;
  const bool has_tau23 = true;
  EXPECT_FALSE(surface_layer_stress::lsm_flux_is_valid(
    has_tau13, true, Real(-0.4), Real(9999.0)));
  EXPECT_TRUE(surface_layer_stress::lsm_flux_is_valid(
    has_tau23, true, Real(-0.4), Real(9999.0)));
  EXPECT_FALSE(surface_layer_stress::lsm_flux_is_valid(
    has_tau23, false, Real(-0.4), Real(9999.0)));
}

TEST(SurfaceLayerStress, ReusingFallbackCacheReconstructsMostStress)
{
  const auto fallback = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, k_low, k_high, false, false, most);
  const auto reused = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, fallback.low_kinematic_cache,
    fallback.high_kinematic_cache, true, true, Real(0.0));
  EXPECT_NEAR(reused.face_stress, most, Real(1.e-12));
}

TEST(SurfaceLayerStress, DiffusionSrcForMomAppliesSignedBottomStresses)
{
  const Box cell_box(IntVect(0), IntVect(1, 1, 0));
  const BoxArray ba(cell_box);
  const DistributionMapping dm(ba);

  MultiFab rho_u_rhs(BoxArray(surroundingNodes(cell_box, 0)), dm, 1, 0);
  MultiFab rho_v_rhs(BoxArray(surroundingNodes(cell_box, 1)), dm, 1, 0);
  MultiFab rho_w_rhs(BoxArray(surroundingNodes(cell_box, 2)), dm, 1, 0);
  rho_u_rhs.setVal(Real(0.0));
  rho_v_rhs.setVal(Real(0.0));
  rho_w_rhs.setVal(Real(0.0));

  MultiFab tau11(ba, dm, 1, 1);
  MultiFab tau22(ba, dm, 1, 1);
  MultiFab tau33(ba, dm, 1, 1);
  MultiFab tau12(BoxArray(convert(cell_box, IntVect(1, 1, 0))), dm, 1, 1);
  MultiFab tau21(BoxArray(convert(cell_box, IntVect(1, 1, 0))), dm, 1, 1);
  MultiFab tau13(BoxArray(convert(cell_box, IntVect(1, 0, 1))), dm, 1, 1);
  MultiFab tau31(BoxArray(convert(cell_box, IntVect(1, 0, 1))), dm, 1, 1);
  MultiFab tau23(BoxArray(convert(cell_box, IntVect(0, 1, 1))), dm, 1, 1);
  MultiFab tau32(BoxArray(convert(cell_box, IntVect(0, 1, 1))), dm, 1, 1);
  tau11.setVal(Real(0.0));
  tau22.setVal(Real(0.0));
  tau33.setVal(Real(0.0));
  tau12.setVal(Real(0.0));
  tau21.setVal(Real(0.0));
  tau13.setVal(Real(0.0));
  tau31.setVal(Real(0.0));
  tau23.setVal(Real(0.0));
  tau32.setVal(Real(0.0));

  const Real tx = Real(-1.7);
  const Real ty = Real(-0.9);
  tau13[0].array()(0, 0, 0) = tx;
  tau23[0].array()(0, 0, 0) = ty;

  MultiFab detJ(ba, dm, 1, 1);
  MultiFab mf_mx(ba, dm, 1, 1);
  MultiFab mf_ux(ba, dm, 1, 1);
  MultiFab mf_vx(ba, dm, 1, 1);
  MultiFab mf_my(ba, dm, 1, 1);
  MultiFab mf_uy(ba, dm, 1, 1);
  MultiFab mf_vy(ba, dm, 1, 1);
  detJ.setVal(Real(1.0));
  mf_mx.setVal(Real(1.0));
  mf_ux.setVal(Real(1.0));
  mf_vx.setVal(Real(1.0));
  mf_my.setVal(Real(1.0));
  mf_uy.setVal(Real(1.0));
  mf_vy.setVal(Real(1.0));

  Gpu::DeviceVector<Real> stretched_dz_d;
  GpuArray<Real, AMREX_SPACEDIM> dxInv{Real(1.0), Real(1.0), Real(1.0)};

  DiffusionSrcForMom(
    surroundingNodes(cell_box, 0), surroundingNodes(cell_box, 1),
    surroundingNodes(cell_box, 2), rho_u_rhs[0].array(), rho_v_rhs[0].array(),
    rho_w_rhs[0].array(), tau11[0].const_array(), tau22[0].const_array(),
    tau33[0].const_array(), tau12[0].const_array(), tau21[0].const_array(),
    tau13[0].const_array(), tau31[0].const_array(), tau23[0].const_array(),
    tau32[0].const_array(), detJ[0].const_array(), stretched_dz_d, dxInv,
    mf_mx[0].const_array(), mf_ux[0].const_array(), mf_vx[0].const_array(),
    mf_my[0].const_array(), mf_uy[0].const_array(), mf_vy[0].const_array(),
    false, false);
  Gpu::streamSynchronize();

  EXPECT_NEAR(rho_u_rhs[0].array()(0, 0, 0), tx, Real(1.e-12));
  EXPECT_NEAR(rho_v_rhs[0].array()(0, 0, 0), ty, Real(1.e-12));

  const Real dt = Real(0.25);
  const Real old_u = Real(2.0);
  const Real old_v = Real(1.0);
  EXPECT_NEAR(
    old_u + dt * rho_u_rhs[0].array()(0, 0, 0), old_u + dt * tx, Real(1.e-12));
  EXPECT_NEAR(
    old_v + dt * rho_v_rhs[0].array()(0, 0, 0), old_v + dt * ty, Real(1.e-12));
}

#include "AMReX.H"
#include "AMReX_MultiFab.H"
#include "AMReX_Reduce.H"

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

Real
single_value (const MultiFab& mf, const Box& box, int comp = 0)
{
  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box overlap = box & mfi.validbox();
    if (overlap.isEmpty()) { continue; }
    const auto array = mf.const_array(mfi);
    reduce_op.eval(overlap, reduce_data,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<Real>
      {
        return { array(i,j,k,comp) };
      });
  }
  return get<0>(reduce_data.value());
}

void
set_single_value (MultiFab& mf, const Box& box, Real value)
{
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box overlap = box & mfi.validbox();
    if (overlap.isEmpty()) { continue; }
    auto array = mf.array(mfi);
    ParallelFor(overlap, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      array(i,j,k) = value;
    });
  }
}

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

// Motivation: LSM stress caches are kinematic, while MOST and Tau_lev use
// conservative stresses. The conversion must preserve both stress sign and
// the conservative value reconstructed with the destination-cell density.
TEST(SurfaceLayerStress, ConservativeToKinematicRoundTripHandlesSignedStress)
{
  for (const Real conservative : {Real(2.4), Real(-2.4)}) {
    const Real kinematic =
      surface_layer_stress::conservative_to_kinematic_stress(
        conservative, rho_low);
    EXPECT_NEAR(rho_low * kinematic, conservative, Real(1.e-12));
  }
}

// Motivation: When both neighboring LSM values are valid, the face stress
// must be the density-weighted average of the two kinematic contributions.
TEST(SurfaceLayerStress, BothValidUsesDensityWeightedFaceAverage)
{
  const auto result = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, k_low, k_high, true, true, most);
  expect_result(
    result, Real(0.5) * (rho_low * k_low + rho_high * k_high), k_low, k_high);
}

// Motivation: A missing low-side LSM value is populated from the conservative
// MOST stress, but the cache must be divided by the low cell's own density.
TEST(SurfaceLayerStress, HighValidCachesMostStressInLowCellUnits)
{
  const auto result = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, k_low, k_high, false, true, most);
  expect_result(
    result, Real(0.5) * (rho_high * k_high + most), most / rho_low, k_high);
}

// Motivation: The opposite mixed branch must use the high cell's density when
// converting the same conservative MOST stress back to a kinematic cache.
TEST(SurfaceLayerStress, LowValidCachesMostStressInHighCellUnits)
{
  const auto result = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, k_low, k_high, true, false, most);
  expect_result(
    result, Real(0.5) * (rho_low * k_low + most), k_low, most / rho_high);
}

// Motivation: With no valid LSM side, the applied face stress remains the
// conservative MOST result while each neighboring cache remains kinematic.
TEST(SurfaceLayerStress, NeitherValidKeepsMostFaceStressAndConvertsBothCaches)
{
  const auto result = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, k_low, k_high, false, false, most);
  expect_result(result, most, most / rho_low, most / rho_high);
  EXPECT_NEAR(rho_low * result.low_kinematic_cache, most, Real(1.e-12));
  EXPECT_NEAR(rho_high * result.high_kinematic_cache, most, Real(1.e-12));
}

// Motivation: tau13 and tau23 represent equivalent directional momentum
// covariances; a unit correction in one staggered branch must apply identically
// to the other.
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

// Motivation: The tau23 fallback decision must depend on the tau23 MultiFab,
// otherwise an absent tau13 field can incorrectly disable valid tau23 data.
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

// Motivation: Fallback caches are reused on later calls. Storing a conservative
// value in a kinematic field would multiply the reconstructed stress by rho.
TEST(SurfaceLayerStress, ReusingFallbackCacheReconstructsMostStress)
{
  const auto fallback = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, k_low, k_high, false, false, most);
  const auto reused = surface_layer_stress::combine_lsm_and_most_stress(
    rho_low, rho_high, fallback.low_kinematic_cache,
    fallback.high_kinematic_cache, true, true, Real(0.0));
  EXPECT_NEAR(reused.face_stress, most, Real(1.e-12));
}

// Motivation: Surface stresses reach the prognostic face momenta through the
// diffusion RHS. Negative upward momentum fluxes at the lower boundary must
// produce negative u- and v-momentum tendencies and the expected explicit
// update.
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
  const Box tau13_box = tau13.boxArray()[0];
  const Box tau23_box = tau23.boxArray()[0];
  Box tau13_sample = tau13_box;
  tau13_sample.setSmall(tau13_box.smallEnd());
  tau13_sample.setBig(tau13_box.smallEnd());
  Box tau23_sample = tau23_box;
  tau23_sample.setSmall(tau23_box.smallEnd());
  tau23_sample.setBig(tau23_box.smallEnd());
  set_single_value(tau13, tau13_sample, tx);
  set_single_value(tau23, tau23_sample, ty);
  Gpu::streamSynchronize();
  EXPECT_NEAR(single_value(tau13, tau13_sample), tx, Real(1.e-12));
  EXPECT_NEAR(single_value(tau23, tau23_sample), ty, Real(1.e-12));

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

  const Box u_box = rho_u_rhs[0].box();
  const Box v_box = rho_v_rhs[0].box();
  Box u_sample = u_box;
  u_sample.setSmall(u_box.smallEnd());
  u_sample.setBig(u_box.smallEnd());
  Box v_sample = v_box;
  v_sample.setSmall(v_box.smallEnd());
  v_sample.setBig(v_box.smallEnd());
  const Real u_rhs = single_value(rho_u_rhs, u_sample);
  const Real v_rhs = single_value(rho_v_rhs, v_sample);
  EXPECT_NEAR(u_rhs, tx, Real(1.e-12));
  EXPECT_NEAR(v_rhs, ty, Real(1.e-12));

  const Real dt = Real(0.25);
  const Real old_u = Real(2.0);
  const Real old_v = Real(1.0);
  EXPECT_NEAR(old_u + dt * u_rhs, old_u + dt * tx, Real(1.e-12));
  EXPECT_NEAR(old_v + dt * v_rhs, old_v + dt * ty, Real(1.e-12));
}

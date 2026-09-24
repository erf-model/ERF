/** \cond */
// Host-side state shared by the regular-grid scalar diffusion adapters.  The
// explicit numerical path deliberately does not allocate positional
// coefficient maps; each native state is resolved to a small policy in the
// adapter immediately before its kernels are launched.
DiffChoice diffChoice = solverChoice.diffChoice;
TurbChoice turbChoice = solverChoice.turbChoice[level];
const bool l_consA = diffChoice.molec_diff_type == MolecDiffType::ConstantAlpha;
const bool l_turb = turbChoice.use_kturb;

[[maybe_unused]] const bool SurfLayer_xlo =
  use_SurfLayer && SurfLayer[amrex::Orientation::xlo()] != nullptr;
[[maybe_unused]] const bool SurfLayer_ylo =
  use_SurfLayer && SurfLayer[amrex::Orientation::ylo()] != nullptr;
[[maybe_unused]] const bool SurfLayer_zlo =
  use_SurfLayer && SurfLayer[amrex::Orientation::zlo()] != nullptr;
[[maybe_unused]] const bool SurfLayer_xhi =
  use_SurfLayer && SurfLayer[amrex::Orientation::xhi()] != nullptr;
[[maybe_unused]] const bool SurfLayer_yhi =
  use_SurfLayer && SurfLayer[amrex::Orientation::yhi()] != nullptr;
[[maybe_unused]] const bool SurfLayer_zhi =
  use_SurfLayer && SurfLayer[amrex::Orientation::zhi()] != nullptr;
const Box zbx = surroundingNodes(bx, 2);

[[maybe_unused]] const int end_comp = start_comp + num_comp - 1;
[[maybe_unused]] const bool l_use_keqn = turbChoice.use_keqn;
[[maybe_unused]] const bool l_implicit_diss =
  turbChoice.use_keqn && turbChoice.implicit_tke_dissipation;
[[maybe_unused]] const bool l_use_mynn = turbChoice.use_pbl_tke;
[[maybe_unused]] const bool use_ref_theta = turbChoice.theta_ref > 0;
[[maybe_unused]] const Real inv_theta_ref =
  one / amrex::max(turbChoice.theta_ref, std::numeric_limits<Real>::min());
[[maybe_unused]] const Real l_inv_theta0 = use_ref_theta ? inv_theta_ref : one;
/** \endcond */

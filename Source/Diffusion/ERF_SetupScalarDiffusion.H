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
/** \endcond */

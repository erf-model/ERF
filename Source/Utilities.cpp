#include "Utilities.H"

AMREX_GPU_DEVICE
void
pc_cmpTemp(
  const int i, const int j, const int k, amrex::Array4<amrex::Real> const& S)
{
  amrex::Real rhoInv = 1.0 / S(i, j, k, URHO);
  amrex::Real T = S(i, j, k, UTEMP);
  amrex::Real e = S(i, j, k, UEINT) * rhoInv;
  EOS::E2T(e, T);
  S(i, j, k, UTEMP) = T;
}

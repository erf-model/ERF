#include "EOS.H"
#include "Derive.H"
#include "IndexDefines.H"

//void
//erf_dervelx(
//  const amrex::Box& bx,
//  amrex::FArrayBox& derfab,
//  int /*dcomp*/,
//  int /*ncomp*/,
//  const amrex::FArrayBox& datfab,
//  const amrex::Geometry& /*geomdata*/,
//  amrex::Real /*time*/,
//  const int* /*bcrec*/,
//  const int /*level*/)
//{
//  auto const dat = datfab.array();
//  auto velx = derfab.array();
//
//  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
//    velx(i, j, k) = dat(i, j, k, UMX) / dat(i, j, k, URHO);
//  });
//}

//void
//erf_dervely(
//  const amrex::Box& bx,
//  amrex::FArrayBox& derfab,
//  int /*dcomp*/,
//  int /*ncomp*/,
//  const amrex::FArrayBox& datfab,
//  const amrex::Geometry& /*geomdata*/,
//  amrex::Real /*time*/,
//  const int* /*bcrec*/,
//  const int /*level*/)
//{
//  auto const dat = datfab.array();
//  auto vely = derfab.array();
//
//  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
//    vely(i, j, k) = dat(i, j, k, UMY) / dat(i, j, k, URHO);
//  });
//}

//void
//erf_dervelz(
//  const amrex::Box& bx,
//  amrex::FArrayBox& derfab,
//  int /*dcomp*/,
//  int /*ncomp*/,
//  const amrex::FArrayBox& datfab,
//  const amrex::Geometry& /*geomdata*/,
//  amrex::Real /*time*/,
//  const int* /*bcrec*/,
//  const int /*level*/)
//{
//  auto const dat = datfab.array();
//  auto velz = derfab.array();
//
//  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
//    velz(i, j, k) = dat(i, j, k, UMZ) / dat(i, j, k, URHO);
//  });
//}

void
erf_dernull(
  const amrex::Box& /*bx*/,
  amrex::FArrayBox& /*derfab*/,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& /*datfab*/,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  // This routine does nothing -- we use it as a placeholder.
}

void
erf_derpres(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto pfab      = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho   = dat(i, j, k, Density_comp);
    const amrex::Real theta = dat(i, j, k, Theta_comp);
    amrex::Real p = getPgivenRTh(rho,theta);
    pfab(i,j,k) = p;
  });
}

void
erf_dertemp(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto tfab      = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho   = dat(i, j, k, Density_comp);
    const amrex::Real theta = dat(i, j, k, Theta_comp);
    amrex::Real T = getTgivenRTh(rho,theta);
    tfab(i,j,k) = T;
  });
}

#include "Derive.H"
#include "EOS.H"
#include "IndexDefines.H"

namespace derived {

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
    const amrex::Real rhotheta = dat(i, j, k, RhoTheta_comp);
    pfab(i,j,k) = getPgivenRTh(rhotheta);
  });
}

void
erf_dersoundspeed(
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
  auto cfab      = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhotheta = dat(i, j, k, RhoTheta_comp);
    const amrex::Real rho      = dat(i, j, k, Rho_comp);
    cfab(i,j,k) = std::sqrt(Gamma * getPgivenRTh(rhotheta) / rho);
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
    const amrex::Real rho = dat(i, j, k, Rho_comp);
    const amrex::Real rhotheta = dat(i, j, k, RhoTheta_comp);
    tfab(i,j,k) = getTgivenRandRTh(rho,rhotheta);
  });
}

void
erf_dertheta(
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
  auto thetafab  = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho      = dat(i, j, k, Rho_comp);
    const amrex::Real rhotheta = dat(i, j, k, RhoTheta_comp);
    thetafab(i,j,k) = rhotheta / rho;
  });
}

void
erf_derscalar(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto adv0fab   = derfab.array();

  const auto& dx = geomdata.CellSize();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho      = dat(i, j, k, Rho_comp);
    const amrex::Real rho_adv0 = dat(i, j, k, RhoScalar_comp);

    const amrex::Real x = (i+.5) * dx[0];
    const amrex::Real y = (j+.5) * dx[1];
    const amrex::Real z = (k+.5) * dx[2];

    const amrex::Real u = 10.0;

    const amrex::Real xc = 0.5+u*time;
    const amrex::Real yc = 0.5;
    const amrex::Real zc = 0.5;

    const amrex::Real r2d_xy = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
    const amrex::Real r0 = 0.25;

    amrex::Real s_ex;
    if (r2d_xy < 0.5)
        s_ex = 0.25 * (1.0 + std::cos(PI * std::min(r2d_xy, r0) / r0));
    else
        s_ex = 0.0;

    adv0fab(i,j,k) = rho_adv0 / rho - s_ex;
  });

    amrex::Real norm0 = 0.0;
    amrex::Real norm1 = 0.0;
    amrex::Real norm2 = 0.0;
    amrex::Real n = 0.0;

    int ihi = bx.bigEnd()[0];
    int jhi = bx.bigEnd()[1];
    int khi = bx.bigEnd()[2];

    for (int k = 0; k < khi; k++)
    for (int j = 0; j < jhi; j++)
    for (int i = 0; i < ihi; i++)
    {
        const amrex::Real x = (i+.5) * dx[0];
        const amrex::Real y = (j+.5) * dx[1];
        const amrex::Real z = (k+.5) * dx[2];

        const amrex::Real u = 10.0;

        const amrex::Real xc = 0.5+u*time;
        const amrex::Real yc = 0.5;
        const amrex::Real zc = 0.5;

        const amrex::Real r2d_xy = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
        const amrex::Real r0 = 0.25;

        amrex::Real s_ex;
        if (r2d_xy < 0.5)
            s_ex = 0.25 * (1.0 + std::cos(PI * std::min(r2d_xy, r0) / r0));
        else
            s_ex = 0.0;

        amrex::Real err = adv0fab(i,j,k) - s_ex;

        norm0 = std::max(err,norm0);
        norm1 += std::abs(err);
        norm2 += err*err;
        n += 1.0;
    }

   norm1 /= n;
   norm2 /= n;
   norm2 = std::sqrt(norm2);

  amrex::Print() << " NORM0 / NORM1 / NORM2 " << norm0 << " " << norm1 << " " << norm2 << std::endl;
}
}

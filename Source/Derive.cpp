#include "EOS.H"
#include "Derive.H"
#include "IndexDefines.H"

void
pc_dervelx(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto velx = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    velx(i, j, k) = dat(i, j, k, UMX) / dat(i, j, k, URHO);
  });
}

void
pc_dervely(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto vely = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    vely(i, j, k) = dat(i, j, k, UMY) / dat(i, j, k, URHO);
  });
}

void
pc_dervelz(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto velz = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    velz(i, j, k) = dat(i, j, k, UMZ) / dat(i, j, k, URHO);
  });
}

void
pc_dermagvel(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto magvel = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real datinv = 1.0 / dat(i, j, k, URHO);
    const amrex::Real dat1 = (dat(i, j, k, UMX) * datinv);
    const amrex::Real dat2 = (dat(i, j, k, UMY) * datinv);
    const amrex::Real dat3 = (dat(i, j, k, UMZ) * datinv);
    magvel(i, j, k) = sqrt((dat1 * dat1) + (dat2 * dat2) + (dat3 * dat3));
  });
}

void
pc_dermagvort(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int level)
{
  auto const dat = datfab.array();
  auto vort = derfab.array();

  const amrex::Box& gbx = amrex::grow(bx, 1);

  amrex::FArrayBox local(gbx, 3);
  amrex::Elixir local_eli = local.elixir();
  auto larr = local.array();

  // Convert momentum to velocity.
  amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);
    larr(i, j, k, 0) = dat(i, j, k, UMX) * rhoInv;
    larr(i, j, k, 1) = dat(i, j, k, UMY) * rhoInv;
    larr(i, j, k, 2) = dat(i, j, k, UMZ) * rhoInv;
  });

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0);
               , const amrex::Real dy = geomdata.CellSize(1);
               , const amrex::Real dz = geomdata.CellSize(2););

  // Calculate vorticity.
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(vort(i, j, k) = 0.;
                 , const amrex::Real vx =
                     0.5 * (larr(i + 1, j, k, 1) - larr(i - 1, j, k, 1)) / dx;
                 const amrex::Real uy =
                   0.5 * (larr(i, j + 1, k, 0) - larr(i, j - 1, k, 0)) / dy;
                 const amrex::Real v3 = vx - uy;
                 , const amrex::Real wx =
                     0.5 * (larr(i + 1, j, k, 2) - larr(i - 1, j, k, 2)) / dx;

                 const amrex::Real wy =
                   0.5 * (larr(i, j + 1, k, 2) - larr(i, j - 1, k, 2)) / dy;

                 const amrex::Real uz =
                   0.5 * (larr(i, j, k + 1, 0) - larr(i, j, k - 1, 0)) / dz;
                 const amrex::Real vz =
                   0.5 * (larr(i, j, k + 1, 1) - larr(i, j, k - 1, 1)) / dz;

                 const amrex::Real v1 = wy - vz;
                 const amrex::Real v2 = uz - wx;);
    vort(i, j, k) = sqrt(AMREX_D_TERM(0., +v3 * v3, +v1 * v1 + v2 * v2));
  });
}

void
pc_derdivu(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int level)
{
  auto const dat = datfab.array();
  auto divu = derfab.array();

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0);
               , const amrex::Real dy = geomdata.CellSize(1);
               , const amrex::Real dz = geomdata.CellSize(2););

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(
      const amrex::Real uhi = dat(i + 1, j, k, UMX) / dat(i + 1, j, k, URHO);
      const amrex::Real ulo = dat(i - 1, j, k, UMX) / dat(i - 1, j, k, URHO);
      , const amrex::Real vhi = dat(i, j + 1, k, UMY) / dat(i, j + 1, k, URHO);
      const amrex::Real vlo = dat(i, j - 1, k, UMY) / dat(i, j - 1, k, URHO);
      , const amrex::Real whi = dat(i, j, k + 1, UMZ) / dat(i, j, k + 1, URHO);
      const amrex::Real wlo = dat(i, j, k - 1, UMZ) / dat(i, j, k - 1, URHO););
    divu(i, j, k) =
      0.5 *
      (AMREX_D_TERM((uhi - ulo) / dx, +(vhi - vlo) / dy, +(whi - wlo) / dz));
  });
}

void
pc_derenstrophy(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  int level)
{
  // This routine will derive enstrophy  = 1/2 rho (x_vorticity^2 +
  // y_vorticity^2 + z_vorticity^2)
  auto const dat = datfab.array();
  auto enstrophy = derfab.array();

  const amrex::Box& gbx = amrex::grow(bx, 1);

  amrex::FArrayBox local(gbx, 3);
  amrex::Elixir local_eli = local.elixir();
  auto larr = local.array();

  // Convert momentum to velocity.
  amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhoInv = 1.0 / dat(i, j, k, URHO);
    larr(i, j, k, 0) = dat(i, j, k, UMX) * rhoInv;
    larr(i, j, k, 1) = dat(i, j, k, UMY) * rhoInv;
    larr(i, j, k, 2) = dat(i, j, k, UMZ) * rhoInv;
  });

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0);
               , const amrex::Real dy = geomdata.CellSize(1);
               , const amrex::Real dz = geomdata.CellSize(2););

  // Calculate enstrophy.
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(enstrophy(i, j, k) = 0.;
                 , const amrex::Real vx =
                     0.5 * (larr(i + 1, j, k, 1) - larr(i - 1, j, k, 1)) / dx;
                 const amrex::Real uy =
                   0.5 * (larr(i, j + 1, k, 0) - larr(i, j - 1, k, 0)) / dy;
                 const amrex::Real v3 = vx - uy;
                 , const amrex::Real wx =
                     0.5 * (larr(i + 1, j, k, 2) - larr(i - 1, j, k, 2)) / dx;

                 const amrex::Real wy =
                   0.5 * (larr(i, j + 1, k, 2) - larr(i, j - 1, k, 2)) / dy;

                 const amrex::Real uz =
                   0.5 * (larr(i, j, k + 1, 0) - larr(i, j, k - 1, 0)) / dz;
                 const amrex::Real vz =
                   0.5 * (larr(i, j, k + 1, 1) - larr(i, j, k - 1, 1)) / dz;

                 const amrex::Real v1 = wy - vz;
                 const amrex::Real v2 = uz - wx;);
    enstrophy(i, j, k) = 0.5 * dat(i, j, k, URHO) *
                         (AMREX_D_TERM(0., +v3 * v3, +v1 * v1 + v2 * v2));
  });
}

void
pc_dernull(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  // This routine is used by particle_count.  Yes it does nothing.
}

void
pc_dersoundspeed(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto cfab = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    const amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real c;
    EOS::RT2Cs(rho, T, c);
    cfab(i, j, k) = c;
  });
}

void
pc_dermachnumber(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto mach = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    const amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real c;
    EOS::RT2Cs(rho, T, c);
    const amrex::Real datxsq = dat(i, j, k, UMX) * dat(i, j, k, UMX);
    const amrex::Real datysq = dat(i, j, k, UMY) * dat(i, j, k, UMY);
    const amrex::Real datzsq = dat(i, j, k, UMZ) * dat(i, j, k, UMZ);
    mach(i, j, k) = sqrt(datxsq + datysq + datzsq) / dat(i, j, k, URHO) / c;
  });
}

void
pc_derpres(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto pfab = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real e = dat(i, j, k, UEINT) * rhoInv;
    amrex::Real p;
    EOS::RT2P(rho, T, p);
    pfab(i, j, k) = p;
  });
}

void
pc_dertemp(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto tfab = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    tfab(i, j, k) = dat(i, j, k, UTEMP);
  });
}

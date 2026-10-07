#include <AMReX_BCRec.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include <ERF_Diffusion.H>
#include <ERF_IndexDefines.H>

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using namespace amrex;

// Terrain-fitted momentum diffusion (ComputeStrain_T + ComputeStressVarVisc_T +
// DiffusionSrcForMom) on a mesh tilted by a constant slope, z = zeta + sx*x + sy*y,
// with K_h != K_v.  erf-model/ERF#4214 reported two defects:
//   1. the zeta-face stresses tau13/tau23 applied K_v (Mom_v) to the projected
//      horizontal stresses h*S11, h*S12 (non-symmetric, anti-diffusive when K_h h^2 > 2 K_v);
//   2. the cell-centred du/dz (dv/dz) in the S11 (S22) metric term used the faces
//      i-1 and i (j-1 and j) instead of i and i+1 (j and j+1).
// The exact-answer test catches both; the energy test catches the first.

namespace {

constexpr int NG = 3;  // ghost cells filled analytically

Real
tol_for_scale (const Real scale)
{
  return scale * (sizeof(Real) == 8 ? Real(1.e-11) : Real(2.e-5));
}

void
copy_to_host (const FArrayBox& src, FArrayBox& dst)
{
  Gpu::copy(Gpu::deviceToHost, src.dataPtr(0), src.dataPtr(0) + src.size(), dst.dataPtr(0));
  Gpu::streamSynchronize();
}

// sin^2 window that vanishes within m of both ends of [0, L]
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
window (Real s, Real m, Real L) noexcept
{
  if (s <= m || s >= L - m) { return Real(0.0); }
  const Real t = std::sin(Real(3.14159265358979323846)*(s - m)/(L - Real(2.0)*m));
  return t*t;
}

// Bilinear velocity field in index-space coordinates (metres):
//   u = a0 x + a1 y + a2 zeta + a3 x zeta,  v = b0 x + b1 y + b2 zeta + b3 y zeta,  w = c3 * window
struct Coeffs {
  Real a[4];
  Real b[4];
};

struct TerrainStressCase
{
  // Grid and physics
  int nx, ny, nz;
  Real dx, dy, dz;
  Real sx, sy;     // terrain slopes dz/dx, dz/dy
  Real Kh, Kv;     // rho*K for Mom_h and Mom_v
  Real er;         // expansion rate (constant)

  Box valid, bxcc, tbxxy, tbxxz, tbxyz;
  Box ubx, vbx, wbx;  // face boxes (valid) for the RHS
  Box domain;

  FArrayBox u, v, w, z_nd, detJ, mu_turb, er_fab;
  FArrayBox mf;  // all map factors = 1
  FArrayBox s11, s22, s33, s12, s21, s13, s31, s23, s32;
  FArrayBox rhs_u, rhs_v, rhs_w;
  std::vector<BCRec> bcs;

  TerrainStressCase (int nx_, int ny_, int nz_, Real dx_, Real dy_, Real dz_,
                     Real sx_, Real sy_, Real Kh_, Real Kv_, Real er_)
    : nx(nx_), ny(ny_), nz(nz_), dx(dx_), dy(dy_), dz(dz_), sx(sx_), sy(sy_),
      Kh(Kh_), Kv(Kv_), er(er_)
  {
    valid  = Box(IntVect(0), IntVect(nx-1, ny-1, nz-1));
    domain = valid;
    // As in erf_make_tau_terms: a halo cell in x and y for the strains
    bxcc  = grow(valid, IntVect(1,1,0));
    tbxxy = grow(convert(valid, IntVect(1,1,0)), IntVect(1,1,0));
    tbxxz = grow(convert(valid, IntVect(1,0,1)), IntVect(1,1,0));
    tbxyz = grow(convert(valid, IntVect(0,1,1)), IntVect(1,1,0));
    ubx = surroundingNodes(valid, 0);
    vbx = surroundingNodes(valid, 1);
    // As erf_slow_rhs_pre: no z-momentum source on the bottom and top domain faces
    wbx = surroundingNodes(valid, 2); wbx.grow(2, -1);

    u.resize(grow(ubx, NG), 1);
    v.resize(grow(vbx, NG), 1);
    w.resize(grow(wbx, NG), 1);
    z_nd.resize(grow(convert(valid, IntVect(1)), NG), 1);
    detJ.resize(grow(valid, NG), 1);
    mu_turb.resize(grow(valid, NG), EddyDiff::NumDiffs);
    er_fab.resize(grow(valid, NG), 1);
    mf.resize(grow(convert(valid, IntVect(1,1,0)), NG), 1);

    s11.resize(bxcc, 1);  s22.resize(bxcc, 1);  s33.resize(bxcc, 1);
    s12.resize(tbxxy, 1); s21.resize(tbxxy, 1);
    s13.resize(tbxxz, 1); s31.resize(tbxxz, 1);
    s23.resize(tbxyz, 1); s32.resize(tbxyz, 1);
    rhs_u.resize(ubx, 1); rhs_v.resize(vbx, 1); rhs_w.resize(wbx, 1);

    // Periodic-type (interior) x/y and first-order extrapolation in z: no Dirichlet
    // stencils, so every strain comes from the generic interior kernels.
    bcs.resize(BCVars::NumTypes);
    for (auto& bc : bcs) {
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        bc.setLo(d, (d < 2) ? ERFBCType::int_dir : ERFBCType::foextrap);
        bc.setHi(d, (d < 2) ? ERFBCType::int_dir : ERFBCType::foextrap);
      }
    }
  }

  // Kernels live here, not in the constructor (nvcc cannot launch extended lambdas there)
  void init ()
  {
    const Real lsx = sx, lsy = sy, ldx = dx, ldy = dy, ldz = dz;
    auto znd = z_nd.array();
    ParallelFor(z_nd.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      znd(i,j,k) = Real(k)*ldz + lsx*Real(i)*ldx + lsy*Real(j)*ldy;
    });
    detJ.setVal<RunOn::Device>(Real(1.0));  // h_zeta = 1 on a tilted plane
    mf.setVal<RunOn::Device>(Real(1.0));
    er_fab.setVal<RunOn::Device>(er);
    mu_turb.setVal<RunOn::Device>(Real(0.0));
    mu_turb.setVal<RunOn::Device>(Kh, mu_turb.box(), EddyDiff::Mom_h, 1);
    mu_turb.setVal<RunOn::Device>(Kv, mu_turb.box(), EddyDiff::Mom_v, 1);
    Gpu::streamSynchronize();
  }

  // Variable K_h: pattern 0 grows tenfold over the lowest cells (as Smagorinsky2D's K_h does
  // next to the ground on a slope), pattern 1 alternates cell by cell between Kh/10 and Kh,
  // pattern 2 does the same in x and zeta only
  void set_variable_kh (int pattern)
  {
    const Real lKh = Kh;
    auto mu = mu_turb.array();
    ParallelFor(mu_turb.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      Real f;
      if (pattern == 0) {
        f = Real(0.1) + Real(0.9) * (Real(1.0) - std::exp(-Real(amrex::max(k, 0)) / Real(4.0)));
      } else if (pattern == 1) {
        f = ((i + j + k + 64) % 2 == 0) ? Real(0.1) : Real(1.0);
      } else {
        f = ((i + k + 64) % 2 == 0) ? Real(0.1) : Real(1.0);
      }
      mu(i,j,k,EddyDiff::Mom_h) = lKh * f;
    });
    Gpu::streamSynchronize();
  }

  void set_bilinear (const Coeffs& c)
  {
    const Real a0 = c.a[0], a1 = c.a[1], a2 = c.a[2], a3 = c.a[3];
    const Real b0 = c.b[0], b1 = c.b[1], b2 = c.b[2], b3 = c.b[3];
    const Real ldx = dx, ldy = dy, ldz = dz;
    auto ua = u.array(); auto va = v.array();
    ParallelFor(u.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      const Real x = Real(i)*ldx, y = (Real(j)+Real(0.5))*ldy, zeta = (Real(k)+Real(0.5))*ldz;
      ua(i,j,k) = a0*x + a1*y + a2*zeta + a3*x*zeta;
    });
    ParallelFor(v.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      const Real x = (Real(i)+Real(0.5))*ldx, y = Real(j)*ldy, zeta = (Real(k)+Real(0.5))*ldz;
      va(i,j,k) = b0*x + b1*y + b2*zeta + b3*y*zeta;
    });
    w.setVal<RunOn::Device>(Real(0.0));
    Gpu::streamSynchronize();
  }

  // Smooth fields with compact support in zeta (zero within `margin` cells of the bottom and
  // top and in the ghost cells), so that summation by parts has no boundary terms. In x and y
  // the fields are either windowed the same way or independent of that direction, which on a
  // plane mesh makes every stress independent of it, so those fluxes cancel exactly.
  void set_compact (Real kx, Real ky, Real kz, Real Au, Real Av, Real Aw, Real phase, int margin,
                    bool window_x, bool window_y)
  {
    const Real Lx = nx*dx, Ly = ny*dy, Lz = nz*dz;
    const Real mx = margin*dx, my = margin*dy, mz = margin*dz;
    const bool wx = window_x, wy = window_y;
    const Real ldx = dx, ldy = dy, ldz = dz;
    auto ua = u.array(); auto va = v.array(); auto wa = w.array();
    ParallelFor(u.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      const Real x = Real(i)*ldx, y = (Real(j)+Real(0.5))*ldy, zeta = (Real(k)+Real(0.5))*ldz;
      const Real wxy = (wx ? window(x,mx,Lx) : Real(1.0)) * (wy ? window(y,my,Ly) : Real(1.0));
      ua(i,j,k) = Au*wxy*window(zeta,mz,Lz)*std::cos(kx*x + ky*y + kz*zeta + phase);
    });
    ParallelFor(v.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      const Real x = (Real(i)+Real(0.5))*ldx, y = Real(j)*ldy, zeta = (Real(k)+Real(0.5))*ldz;
      const Real wxy = (wx ? window(x,mx,Lx) : Real(1.0)) * (wy ? window(y,my,Ly) : Real(1.0));
      va(i,j,k) = Av*wxy*window(zeta,mz,Lz)*std::sin(kx*x - ky*y + kz*zeta + phase);
    });
    ParallelFor(w.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      const Real x = (Real(i)+Real(0.5))*ldx, y = (Real(j)+Real(0.5))*ldy, zeta = Real(k)*ldz;
      const Real wxy = (wx ? window(x,mx,Lx) : Real(1.0)) * (wy ? window(y,my,Ly) : Real(1.0));
      wa(i,j,k) = Aw*wxy*window(zeta,mz,Lz)*std::cos(Real(0.7)*kx*x + ky*y - kz*zeta + phase);
    });
    Gpu::streamSynchronize();
  }

  // Strain -> stress (-> momentum RHS), with the boxes erf_make_tau_terms uses
  void compute (bool with_rhs)
  {
    GpuArray<Real, AMREX_SPACEDIM> dxInv{Real(1.0)/dx, Real(1.0)/dy, Real(1.0)/dz};
    auto ua = u.const_array(), va = v.const_array(), wa = w.const_array();
    auto znd = z_nd.const_array(), dJ = detJ.const_array(), mfa = mf.const_array();
    auto a11 = s11.array(), a22 = s22.array(), a33 = s33.array();
    auto a12 = s12.array(), a21 = s21.array(), a13 = s13.array(), a31 = s31.array();
    auto a23 = s23.array(), a32 = s32.array();
    Array4<Real> no_corr{};

    Box cc = bxcc, xy = tbxxy, xz = tbxxz, yz = tbxyz;
    ComputeStrain_T(cc, xy, xz, yz, domain, ua, va, wa,
                    a11, a22, a33, a12, a21, a13, a31, a23, a32,
                    znd, dJ, dxInv, mfa, mfa, mfa, mfa, mfa, mfa, bcs.data(),
                    no_corr, no_corr);

    // Remove the halo for the off-diagonal stresses, as erf_make_tau_terms does
    xy.grow(IntVect(-1,-1,0)); xz.grow(IntVect(-1,-1,0)); yz.grow(IntVect(-1,-1,0));
    Array4<const Real> no_cell_data{};
    ComputeStressVarVisc_T(cc, xy, xz, yz, Real(0.0), mu_turb.const_array(), no_cell_data,
                           a11, a22, a33, a12, a21, a13, a31, a23, a32,
                           er_fab.const_array(), znd, dJ, dxInv,
                           mfa, mfa, mfa, mfa, mfa, mfa, no_corr, no_corr, no_corr);

    if (with_rhs) {
      rhs_u.setVal<RunOn::Device>(Real(0.0));
      rhs_v.setVal<RunOn::Device>(Real(0.0));
      rhs_w.setVal<RunOn::Device>(Real(0.0));
      Gpu::DeviceVector<Real> no_stretched_dz;
      DiffusionSrcForMom(ubx, vbx, wbx, rhs_u.array(), rhs_v.array(), rhs_w.array(),
                         s11.const_array(), s22.const_array(), s33.const_array(),
                         s12.const_array(), s21.const_array(),
                         s13.const_array(), s31.const_array(),
                         s23.const_array(), s32.const_array(),
                         dJ, no_stretched_dz, dxInv, mfa, mfa, mfa, mfa, mfa, mfa,
                         false, true);
    }
    Gpu::streamSynchronize();
  }

  // Discrete kinetic-energy rate: sum over faces of (face volume) * vel * RHS.  With detJ = 1
  // and unit map factors the face volume is dx*dy*dz everywhere.
  Real energy_rate ()
  {
    FArrayBox hu(u.box(), 1, The_Pinned_Arena()), hv(v.box(), 1, The_Pinned_Arena()), hw(w.box(), 1, The_Pinned_Arena());
    FArrayBox hru(ubx, 1, The_Pinned_Arena()), hrv(vbx, 1, The_Pinned_Arena()), hrw(wbx, 1, The_Pinned_Arena());
    copy_to_host(u, hu); copy_to_host(v, hv); copy_to_host(w, hw);
    copy_to_host(rhs_u, hru); copy_to_host(rhs_v, hrv); copy_to_host(rhs_w, hrw);
    double sum = 0.0;
    const auto ua = hu.const_array(), va = hv.const_array(), wa = hw.const_array();
    const auto ru = hru.const_array(), rv = hrv.const_array(), rw = hrw.const_array();
    LoopOnCpu(ubx, [&] (int i, int j, int k) { sum += double(ua(i,j,k))*double(ru(i,j,k)); });
    LoopOnCpu(vbx, [&] (int i, int j, int k) { sum += double(va(i,j,k))*double(rv(i,j,k)); });
    LoopOnCpu(wbx, [&] (int i, int j, int k) { sum += double(wa(i,j,k))*double(rw(i,j,k)); });
    return Real(sum * double(dx*dy*dz));
  }

  // Sum of the squared velocity gradients weighted by K, the scale of the dissipation
  Real dissipation_scale ()
  {
    FArrayBox hu(u.box(), 1, The_Pinned_Arena()), hv(v.box(), 1, The_Pinned_Arena()), hw(w.box(), 1, The_Pinned_Arena());
    copy_to_host(u, hu); copy_to_host(v, hv); copy_to_host(w, hw);
    double sum = 0.0;
    const auto ua = hu.const_array(), va = hv.const_array(), wa = hw.const_array();
    const double K = double(std::max(Kh, Kv));
    LoopOnCpu(ubx, [&] (int i, int j, int k) {
      const double g = (double(ua(i,j,k+1)) - double(ua(i,j,k)))/double(dz);
      sum += K*g*g;
    });
    LoopOnCpu(vbx, [&] (int i, int j, int k) {
      const double g = (double(va(i,j,k+1)) - double(va(i,j,k)))/double(dz);
      sum += K*g*g;
    });
    LoopOnCpu(wbx, [&] (int i, int j, int k) {
      const double g = (double(wa(i+1,j,k)) - double(wa(i,j,k)))/double(dx);
      sum += K*g*g;
    });
    return Real(sum * double(dx*dy*dz));
  }
};

// Exact strains of the bilinear field, deviatoric (minus er/3) where ERF subtracts it
struct Exact {
  Coeffs c; Real sx, sy, er;
  Real dudx (Real /*x*/, Real zeta) const { return c.a[0] + c.a[3]*zeta; }
  Real dudz (Real x)                const { return c.a[2] + c.a[3]*x; }
  Real dvdy (Real /*y*/, Real zeta) const { return c.b[1] + c.b[3]*zeta; }
  Real dvdz (Real y)                const { return c.b[2] + c.b[3]*y; }
  Real S11 (Real x, Real zeta) const { return dudx(x,zeta) - sx*dudz(x) - er/Real(3.0); }
  Real S22 (Real y, Real zeta) const { return dvdy(y,zeta) - sy*dvdz(y) - er/Real(3.0); }
  Real S12 (Real x, Real y)    const { return Real(0.5)*(c.a[1] + c.b[0] - sy*dudz(x) - sx*dvdz(y)); }
  Real S13 (Real x)            const { return Real(0.5)*dudz(x); }
  Real S23 (Real y)            const { return Real(0.5)*dvdz(y); }
};

const Coeffs bilinear_coeffs{ {Real(1.e-3), Real(-2.e-3), Real(0.02), Real(1.e-5)},
                              {Real(-1.5e-3), Real(5.e-4), Real(-0.01), Real(-2.e-5)} };

} // namespace

// Motivation (#4214, defects 1 and 2): on a tilted mesh with K_h != K_v and a bilinear
// velocity, every second-order stencil is exact, so the stresses must equal the closed form
//   tau11 = -2 K_h S11,  tau13 = -(2 K_v S13 - 2 K_h (h_xi S11 + h_eta S12))  (and tau22, tau23).
// The old K_v projection is off by 2 (K_h - K_v) h S11 (~0.8 here); the old one-cell
// offset in the cell-centred du/dz is off by K_h h^2 a3 dx (~0.01 in tau13, 0.0048 in tau11).
// Bottom, top and interior tau13/tau23 blocks are all checked.
TEST(TerrainStress, ExactStressesOnTiltedMeshWithAnisotropicK)
{
  const Real sx = Real(0.3), sy = Real(-0.2), Kh = Real(40.0), Kv = Real(1.5), er = Real(3.e-3);
  TerrainStressCase c(6, 5, 8, Real(200.0), Real(150.0), Real(50.0), sx, sy, Kh, Kv, er);
  c.init();
  c.set_bilinear(bilinear_coeffs);
  c.compute(false);

  const Exact ex{bilinear_coeffs, sx, sy, er};

  FArrayBox h11(c.s11.box(), 1, The_Pinned_Arena()), h22(c.s22.box(), 1, The_Pinned_Arena());
  FArrayBox h13(c.s13.box(), 1, The_Pinned_Arena()), h23(c.s23.box(), 1, The_Pinned_Arena());
  FArrayBox h31(c.s31.box(), 1, The_Pinned_Arena());
  copy_to_host(c.s11, h11); copy_to_host(c.s22, h22);
  copy_to_host(c.s13, h13); copy_to_host(c.s23, h23); copy_to_host(c.s31, h31);
  const auto t11 = h11.const_array(), t22 = h22.const_array();
  const auto t13 = h13.const_array(), t23 = h23.const_array(), t31 = h31.const_array();

  const Real tol = tol_for_scale(Real(1.0));
  int n11 = 0, n13 = 0, n23 = 0;

  // tau11, tau22 at cell centres (valid cells)
  LoopOnCpu(c.valid, [&] (int i, int j, int k) {
    const Real x = (Real(i)+Real(0.5))*c.dx, y = (Real(j)+Real(0.5))*c.dy, zeta = (Real(k)+Real(0.5))*c.dz;
    EXPECT_NEAR(t11(i,j,k), -Real(2.0)*Kh*ex.S11(x,zeta), tol) << "tau11 at " << i << " " << j << " " << k;
    EXPECT_NEAR(t22(i,j,k), -Real(2.0)*Kh*ex.S22(y,zeta), tol) << "tau22 at " << i << " " << j << " " << k;
    ++n11;
  });

  // tau13 on every xz edge of the valid region, k = 0 (bottom block) .. nz (top block)
  const Box xz = convert(c.valid, IntVect(1,0,1));
  LoopOnCpu(xz, [&] (int i, int j, int k) {
    const Real x = Real(i)*c.dx, y = (Real(j)+Real(0.5))*c.dy, zeta = Real(k)*c.dz;
    const Real proj = sx*ex.S11(x,zeta) + sy*ex.S12(x,y);
    EXPECT_NEAR(t13(i,j,k), -(Real(2.0)*Kv*ex.S13(x) - Real(2.0)*Kh*proj), tol)
      << "tau13 at " << i << " " << j << " " << k;
    // tau31 is a vertical stress: K_v only
    EXPECT_NEAR(t31(i,j,k), -Real(2.0)*Kv*ex.S13(x), tol) << "tau31 at " << i << " " << j << " " << k;
    ++n13;
  });

  const Box yz = convert(c.valid, IntVect(0,1,1));
  LoopOnCpu(yz, [&] (int i, int j, int k) {
    const Real x = (Real(i)+Real(0.5))*c.dx, y = Real(j)*c.dy, zeta = Real(k)*c.dz;
    const Real proj = sx*ex.S12(x,y) + sy*ex.S22(y,zeta);
    EXPECT_NEAR(t23(i,j,k), -(Real(2.0)*Kv*ex.S23(y) - Real(2.0)*Kh*proj), tol)
      << "tau23 at " << i << " " << j << " " << k;
    ++n23;
  });

  EXPECT_EQ(n11, 6*5*8);
  EXPECT_EQ(n13, 7*5*9);
  EXPECT_EQ(n23, 6*6*9);
}

// Motivation (#4214, defect 1): the discrete operator must dissipate kinetic energy for
// constant K and any velocity. With compact-support fields summation by parts gives
//   sum(u . RHS) = -[2 K_h sum(S11^2 + S22^2 + 2 S12^2) + K_v (...)] <= 0
// only when the zeta flux uses K_h on the projected terms (and the cell-centred du/dz is the
// average of the four edges around the cell). The first two modes are the oblique waves the
// old K_v form anti-diffuses (dx u / dz u ~ (K_h+K_v) h / (2 K_h), K_h h^2 = 16 > 2 K_v):
// rate/scale is -0.195 with the fix and +0.040 with the K_v form for both.
TEST(TerrainStress, ConstantKOperatorDissipatesEnergyOnSlopes)
{
  const Real sx = Real(0.4), sy = Real(-0.4), Kh = Real(100.0), Kv = Real(1.0);
  TerrainStressCase c(30, 32, 32, Real(200.0), Real(200.0), Real(50.0), sx, sy, Kh, Kv, Real(0.0));
  c.init();

  const Real two_pi = Real(2.0)*Real(3.14159265358979323846);
  struct Mode { Real kx, ky, kz, Au, Av, Aw, phase; bool window_x, window_y; };
  const Mode modes[] = {
    // u(x, zeta) only, the oblique wave the K_v form anti-diffuses most:
    // dx u / dz u = (K_h + K_v) h_xi / (2 K_h) ~ 0.2, wavelengths 2000 m in x and 400 m in zeta
    {two_pi/Real(2000.0), Real(0.0), two_pi/Real(400.0), Real(1.0), Real(0.0), Real(0.0), Real(0.0), true, false},
    // v(y, zeta) only, the same in y (h_eta = -0.4: dy v / dz v ~ -0.2, 2000 m and 400 m;
    // the v field is sin(kx x - ky y + kz zeta), so ky > 0 gives the negative ratio)
    {Real(0.0), two_pi/Real(2000.0), two_pi/Real(400.0), Real(0.0), Real(1.0), Real(0.0), Real(0.0), false, true},
    // all three components, windowed in every direction
    {two_pi/Real(1500.0), two_pi/Real(900.0), two_pi/Real(300.0), Real(1.0), Real(-0.7), Real(0.4), Real(1.1), true, true},
    {two_pi/Real(3000.0), two_pi/Real(1200.0), two_pi/Real(600.0), Real(-0.5), Real(0.8), Real(-0.9), Real(2.0), true, true},
  };
  int m = 0;
  for (const auto& md : modes) {
    c.set_compact(md.kx, md.ky, md.kz, md.Au, md.Av, md.Aw, md.phase, 3, md.window_x, md.window_y);
    c.compute(true);
    const Real rate  = c.energy_rate();
    const Real scale = c.dissipation_scale();
    ASSERT_GT(scale, Real(0.0));
    EXPECT_LT(rate, tol_for_scale(scale)) << "mode " << m << ": energy grows, rate/scale = " << rate/scale;
    ++m;
  }
}

// Motivation (#4214): Smagorinsky2D's K_h varies strongly from cell to cell next to the
// ground, so the zeta flux must dissipate energy for variable K_h as well.  Summation by parts
// gives -[2 sum K_h(S11^2 + ...) + ...] <= 0 exactly when the projected terms average the
// K_h-weighted stresses (K_h S11 at the cells, K_h S12 at the edges) to the zeta edge, the
// transpose of the cell-centred metric term in tau11.
TEST(TerrainStress, VariableKOperatorDissipatesEnergyOnSlopes)
{
  const Real sx = Real(0.4), sy = Real(-0.4), Kh = Real(100.0), Kv = Real(1.0);
  const Real two_pi = Real(2.0)*Real(3.14159265358979323846);
  for (int pattern = 0; pattern < 2; ++pattern) {
    TerrainStressCase c(30, 32, 32, Real(200.0), Real(200.0), Real(50.0), sx, sy, Kh, Kv, Real(0.0));
    c.init();
    c.set_variable_kh(pattern);
    struct Mode { Real kx, ky, kz, Au, Av, Aw, phase; bool window_x, window_y; };
    const Mode modes[] = {
      {two_pi/Real(2000.0), Real(0.0), two_pi/Real(400.0), Real(1.0), Real(0.0), Real(0.0), Real(0.0), true, false},
      {Real(0.0), two_pi/Real(2000.0), two_pi/Real(400.0), Real(0.0), Real(1.0), Real(0.0), Real(0.0), false, true},
      {two_pi/Real(1500.0), two_pi/Real(900.0), two_pi/Real(300.0), Real(1.0), Real(-0.7), Real(0.4), Real(1.1), true, true},
      {two_pi/Real(800.0), Real(0.0), two_pi/Real(200.0), Real(1.0), Real(0.0), Real(0.0), Real(0.5), true, false},
    };
    int m = 0;
    for (const auto& md : modes) {
      c.set_compact(md.kx, md.ky, md.kz, md.Au, md.Av, md.Aw, md.phase, 3, md.window_x, md.window_y);
      c.compute(true);
      const Real rate  = c.energy_rate();
      const Real scale = c.dissipation_scale();
      ASSERT_GT(scale, Real(0.0));
      EXPECT_LT(rate, tol_for_scale(scale)) << "pattern " << pattern << " mode " << m
                                            << ": energy grows, rate/scale = " << rate/scale;
      ++m;
    }
  }

  // The worst case found by an eigenvalue analysis of the u operator: a grid-scale oblique
  // wave (four cells in x and in zeta) on a slope of one cell per cell over an x-zeta
  // checkerboard of K_h.  Averaging the product of the edge-averaged K_h and the averaged
  // strain gives this field a growing energy (the operator is indefinite); averaging the
  // K_h-weighted stresses keeps it dissipative.
  {
    TerrainStressCase c(24, 4, 24, Real(50.0), Real(50.0), Real(50.0), Real(1.0), Real(0.0), Kh, Kv, Real(0.0));
    c.init();
    c.set_variable_kh(2);
    c.set_compact(two_pi/Real(200.0), Real(0.0), two_pi/Real(200.0), Real(1.0), Real(0.0), Real(0.0),
                  Real(0.0), 2, true, false);
    c.compute(true);
    const Real rate  = c.energy_rate();
    const Real scale = c.dissipation_scale();
    EXPECT_LT(rate, tol_for_scale(scale)) << "checkerboard K_h: energy grows, rate/scale = " << rate/scale;
  }
}

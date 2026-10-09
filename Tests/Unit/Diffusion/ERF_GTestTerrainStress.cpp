#include <AMReX_BCRec.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include <ERF_Diffusion.H>
#include <ERF_IndexDefines.H>

#include <gtest/gtest.h>

#include <algorithm>
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
  FArrayBox s13i, s23i, s33i;  // the parts of tau13/tau23/tau33 the implicit solve takes
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
    s13i.resize(tbxxz, 1); s23i.resize(tbxyz, 1); s33i.resize(bxcc, 1);

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
  void compute (bool with_rhs, bool implicit_metric = false)
  {
    GpuArray<Real, AMREX_SPACEDIM> dxInv{Real(1.0)/dx, Real(1.0)/dy, Real(1.0)/dz};
    auto ua = u.const_array(), va = v.const_array(), wa = w.const_array();
    auto znd = z_nd.const_array(), dJ = detJ.const_array(), mfa = mf.const_array();
    auto a11 = s11.array(), a22 = s22.array(), a33 = s33.array();
    auto a12 = s12.array(), a21 = s21.array(), a13 = s13.array(), a31 = s31.array();
    auto a23 = s23.array(), a32 = s32.array();
    Array4<Real> c13 = s13i.array(), c23 = s23i.array(), c33 = s33i.array();

    Box cc = bxcc, xy = tbxxy, xz = tbxxz, yz = tbxyz;
    ComputeStrain_T(cc, xy, xz, yz, domain, ua, va, wa,
                    a11, a22, a33, a12, a21, a13, a31, a23, a32,
                    znd, dJ, dxInv, mfa, mfa, mfa, mfa, mfa, mfa, bcs.data(),
                    c13, c23);

    // Remove the halo for the off-diagonal stresses, as erf_make_tau_terms does
    xy.grow(IntVect(-1,-1,0)); xz.grow(IntVect(-1,-1,0)); yz.grow(IntVect(-1,-1,0));
    Array4<const Real> no_cell_data{};
    ComputeStressVarVisc_T(cc, xy, xz, yz, Real(0.0), mu_turb.const_array(), no_cell_data,
                           a11, a22, a33, a12, a21, a13, a31, a23, a32,
                           er_fab.const_array(), znd, dJ, dxInv,
                           mfa, mfa, mfa, mfa, mfa, mfa, c13, c23, c33, implicit_metric);

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

  // Strain -> stress on one x-y tile, as erf_make_tau_terms does it under tiling: the strains live
  // in the tile's own arrays on its halo-grown boxes, and a nodal tile box reaches the high
  // node only on the last tile in that direction.  tau13 and tau23 on the tile's edges are
  // copied into out13 and out23.
  void compute_tile (const Box& tile, FArrayBox& out13, FArrayBox& out23)
  {
    GpuArray<Real, AMREX_SPACEDIM> dxInv{Real(1.0)/dx, Real(1.0)/dy, Real(1.0)/dz};
    auto ua = u.const_array(), va = v.const_array(), wa = w.const_array();
    auto znd = z_nd.const_array(), dJ = detJ.const_array(), mfa = mf.const_array();
    auto nodal_tile = [&] (const IntVect& typ) {
      Box b = convert(tile, typ);
      for (int d = 0; d < 2; ++d) {
        if (typ[d] == 1 && tile.bigEnd(d) != valid.bigEnd(d)) { b.setBig(d, tile.bigEnd(d)); }
      }
      return b;
    };
    Box cc = grow(tile, IntVect(1,1,0));
    Box xy = grow(nodal_tile(IntVect(1,1,0)), IntVect(1,1,0));
    Box xz = grow(nodal_tile(IntVect(1,0,1)), IntVect(1,1,0));
    Box yz = grow(nodal_tile(IntVect(0,1,1)), IntVect(1,1,0));
    FArrayBox t11(cc,1), t22(cc,1), t33(cc,1), t12(xy,1), t21(xy,1);
    FArrayBox t13(xz,1), t31(xz,1), t23(yz,1), t32(yz,1), i13(xz,1), i23(yz,1), i33(cc,1);
    auto a11 = t11.array(), a22 = t22.array(), a33 = t33.array(), a12 = t12.array(), a21 = t21.array();
    auto a13 = t13.array(), a31 = t31.array(), a23 = t23.array(), a32 = t32.array();
    Array4<Real> c13 = i13.array(), c23 = i23.array(), c33 = i33.array();
    ComputeStrain_T(cc, xy, xz, yz, domain, ua, va, wa,
                    a11, a22, a33, a12, a21, a13, a31, a23, a32,
                    znd, dJ, dxInv, mfa, mfa, mfa, mfa, mfa, mfa, bcs.data(), c13, c23);
    xy.grow(IntVect(-1,-1,0)); xz.grow(IntVect(-1,-1,0)); yz.grow(IntVect(-1,-1,0));
    Array4<const Real> no_cell_data{};
    ComputeStressVarVisc_T(cc, xy, xz, yz, Real(0.0), mu_turb.const_array(), no_cell_data,
                           a11, a22, a33, a12, a21, a13, a31, a23, a32,
                           er_fab.const_array(), znd, dJ, dxInv,
                           mfa, mfa, mfa, mfa, mfa, mfa, c13, c23, c33, false);
    auto o13 = out13.array(), o23 = out23.array();
    ParallelFor(xz, yz,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { o13(i,j,k) = a13(i,j,k); },
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { o23(i,j,k) = a23(i,j,k); });
    Gpu::streamSynchronize();
  }

  // Strain -> constant-viscosity stress (ComputeStressConsVisc_T, molecular viscosity only)
  void compute_cons (Real mu_eff, bool implicit_metric)
  {
    GpuArray<Real, AMREX_SPACEDIM> dxInv{Real(1.0)/dx, Real(1.0)/dy, Real(1.0)/dz};
    auto ua = u.const_array(), va = v.const_array(), wa = w.const_array();
    auto znd = z_nd.const_array(), dJ = detJ.const_array(), mfa = mf.const_array();
    auto a11 = s11.array(), a22 = s22.array(), a33 = s33.array();
    auto a12 = s12.array(), a21 = s21.array(), a13 = s13.array(), a31 = s31.array();
    auto a23 = s23.array(), a32 = s32.array();
    Array4<Real> c13 = s13i.array(), c23 = s23i.array(), c33 = s33i.array();
    Box cc = bxcc, xy = tbxxy, xz = tbxxz, yz = tbxyz;
    ComputeStrain_T(cc, xy, xz, yz, domain, ua, va, wa,
                    a11, a22, a33, a12, a21, a13, a31, a23, a32,
                    znd, dJ, dxInv, mfa, mfa, mfa, mfa, mfa, mfa, bcs.data(), c13, c23);
    xy.grow(IntVect(-1,-1,0)); xz.grow(IntVect(-1,-1,0)); yz.grow(IntVect(-1,-1,0));
    Array4<const Real> no_cell_data{};
    ComputeStressConsVisc_T(cc, xy, xz, yz, mu_eff, no_cell_data,
                            a11, a22, a33, a12, a21, a13, a31, a23, a32,
                            er_fab.const_array(), znd, dJ, dxInv,
                            mfa, mfa, mfa, mfa, mfa, mfa, c13, c23, c33, implicit_metric);
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

// Motivation (erf.implicit_terrain_metric): the part of tau13/tau23 the implicit solve removes
// and re-solves (tau13i/tau23i) must also hold the compact metric term K_h M du/dz, with
// M_u = 2 h_xi^2 + h_eta^2 and M_v = h_xi^2 + 2 h_eta^2, on the interior faces -- the faces
// the solve gives that coefficient -- and only the K_v part on the bottom and top faces.
// With the option off it is the K_v part everywhere, as before.
TEST(TerrainStress, ImplicitPartHoldsTheMetricTermOnInteriorFaces)
{
  const Real sx = Real(0.3), sy = Real(-0.2), Kh = Real(40.0), Kv = Real(1.5), er = Real(3.e-3);
  const Exact ex{bilinear_coeffs, sx, sy, er};
  const Real Mu = Real(2.0)*sx*sx + sy*sy, Mv = sx*sx + Real(2.0)*sy*sy;
  for (const bool metric : {false, true}) {
    TerrainStressCase c(6, 5, 8, Real(200.0), Real(150.0), Real(50.0), sx, sy, Kh, Kv, er);
    c.init();
    c.set_bilinear(bilinear_coeffs);
    c.compute(false, metric);
    FArrayBox h13(c.s13i.box(), 1, The_Pinned_Arena()), h23(c.s23i.box(), 1, The_Pinned_Arena());
    copy_to_host(c.s13i, h13); copy_to_host(c.s23i, h23);
    const auto t13 = h13.const_array(), t23 = h23.const_array();
    LoopOnCpu(convert(c.valid, IntVect(1,0,1)), [&] (int i, int j, int k) {
      const bool interior = (k > 0 && k < c.nz);
      const Real K = Kv + ((metric && interior) ? Kh*Mu : Real(0.0));
      EXPECT_NEAR(t13(i,j,k), -K*Real(2.0)*ex.S13(Real(i)*c.dx), tol_for_scale(Real(1.0)))
        << "metric " << metric << " tau13i at " << i << " " << j << " " << k;
    });
    LoopOnCpu(convert(c.valid, IntVect(0,1,1)), [&] (int i, int j, int k) {
      const bool interior = (k > 0 && k < c.nz);
      const Real K = Kv + ((metric && interior) ? Kh*Mv : Real(0.0));
      EXPECT_NEAR(t23(i,j,k), -K*Real(2.0)*ex.S23(Real(j)*c.dy), tol_for_scale(Real(1.0)))
        << "metric " << metric << " tau23i at " << i << " " << j << " " << k;
    });
  }
}

// Motivation (erf.implicit_terrain_metric): with only a molecular viscosity the stress goes
// through ComputeStressConsVisc_T, and the implicit solve then adds mu (1 + M) on the interior
// faces, so its tau13i/tau23i must hold the same metric term there.
TEST(TerrainStress, ConstantViscosityImplicitPartHoldsTheMetricTerm)
{
  const Real sx = Real(0.3), sy = Real(-0.2), mu_eff = Real(3.0), er = Real(3.e-3);
  const Exact ex{bilinear_coeffs, sx, sy, er};
  const Real Mu = Real(2.0)*sx*sx + sy*sy, Mv = sx*sx + Real(2.0)*sy*sy;
  for (const bool metric : {false, true}) {
    TerrainStressCase c(6, 5, 8, Real(200.0), Real(150.0), Real(50.0), sx, sy, Real(0.0), Real(0.0), er);
    c.init();
    c.set_bilinear(bilinear_coeffs);
    c.compute_cons(mu_eff, metric);
    FArrayBox h13(c.s13i.box(), 1, The_Pinned_Arena()), h23(c.s23i.box(), 1, The_Pinned_Arena());
    copy_to_host(c.s13i, h13); copy_to_host(c.s23i, h23);
    const auto t13 = h13.const_array(), t23 = h23.const_array();
    LoopOnCpu(convert(c.valid, IntVect(1,0,1)), [&] (int i, int j, int k) {
      const bool interior = (k > 0 && k < c.nz);
      const Real f = Real(1.0) + ((metric && interior) ? Mu : Real(0.0));
      EXPECT_NEAR(t13(i,j,k), -mu_eff*f*ex.S13(Real(i)*c.dx), tol_for_scale(Real(1.0)))
        << "metric " << metric << " tau13i at " << i << " " << j << " " << k;
    });
    LoopOnCpu(convert(c.valid, IntVect(0,1,1)), [&] (int i, int j, int k) {
      const bool interior = (k > 0 && k < c.nz);
      const Real f = Real(1.0) + ((metric && interior) ? Mv : Real(0.0));
      EXPECT_NEAR(t23(i,j,k), -mu_eff*f*ex.S23(Real(j)*c.dy), tol_for_scale(Real(1.0)))
        << "metric " << metric << " tau23i at " << i << " " << j << " " << k;
    });
  }
}

// Motivation (#4214, found on a 3-D hill): the K_h-weighted stresses are formed in temporaries,
// and under tiling a nodal tile box stops one node short of the (i,j+1) and (i+1,j) reads of
// the zeta-edge averages.  Splitting the box into 2 x 2 tiles in x and y must give exactly the
// tau13 and tau23 of the whole box, with both slopes and a K_h that varies cell by cell.  (The
// stress loop never tiles in z: erf_make_tau_terms iterates with TileNoZ.)
TEST(TerrainStress, TiledStressesMatchTheWholeBox)
{
  const Real sx = Real(0.3), sy = Real(-0.25), Kh = Real(60.0), Kv = Real(1.5);
  const Real two_pi = Real(2.0)*Real(3.14159265358979323846);
  TerrainStressCase c(8, 8, 6, Real(200.0), Real(150.0), Real(50.0), sx, sy, Kh, Kv, Real(0.0));
  c.init();
  c.set_variable_kh(1);
  c.set_compact(two_pi/Real(900.0), two_pi/Real(700.0), two_pi/Real(200.0),
                Real(1.0), Real(-0.8), Real(0.5), Real(0.4), 0, false, false);
  c.compute(false);
  FArrayBox w13(c.s13.box(), 1, The_Pinned_Arena()), w23(c.s23.box(), 1, The_Pinned_Arena());
  copy_to_host(c.s13, w13); copy_to_host(c.s23, w23);

  FArrayBox t13(c.s13.box(), 1), t23(c.s23.box(), 1);
  t13.setVal<RunOn::Device>(Real(0.0)); t23.setVal<RunOn::Device>(Real(0.0));
  for (int jt = 0; jt < 2; ++jt) {
    for (int it = 0; it < 2; ++it) {
      const Box tile(IntVect(4*it, 4*jt, 0), IntVect(4*it+3, 4*jt+3, c.nz-1));
      c.compute_tile(tile, t13, t23);
    }
  }
  FArrayBox h13(t13.box(), 1, The_Pinned_Arena()), h23(t23.box(), 1, The_Pinned_Arena());
  copy_to_host(t13, h13); copy_to_host(t23, h23);
  const auto a13 = w13.const_array(), a23 = w23.const_array();
  const auto b13 = h13.const_array(), b23 = h23.const_array();
  Real scale = Real(0.0);
  LoopOnCpu(convert(c.valid, IntVect(1,0,1)), [&] (int i, int j, int k) {
    scale = std::max(scale, std::abs(a13(i,j,k)));
  });
  ASSERT_GT(scale, Real(0.0));
  LoopOnCpu(convert(c.valid, IntVect(1,0,1)), [&] (int i, int j, int k) {
    EXPECT_NEAR(b13(i,j,k), a13(i,j,k), tol_for_scale(scale)) << "tau13 at " << i << " " << j << " " << k;
  });
  LoopOnCpu(convert(c.valid, IntVect(0,1,1)), [&] (int i, int j, int k) {
    EXPECT_NEAR(b23(i,j,k), a23(i,j,k), tol_for_scale(scale)) << "tau23 at " << i << " " << j << " " << k;
  });
}

namespace {

// Dense helpers for the operator test below (row-major n x n)
using Mat = std::vector<double>;
Mat matmul (const Mat& a, const Mat& b, int n)
{
  Mat c(static_cast<std::size_t>(n)*n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      const double aik = a[static_cast<std::size_t>(i)*n + k];
      if (aik == 0.0) { continue; }
      for (int j = 0; j < n; ++j) { c[static_cast<std::size_t>(i)*n + j] += aik*b[static_cast<std::size_t>(k)*n + j]; }
    }
  }
  return c;
}
// exp(M) by scaling and squaring of a Taylor series
Mat expm (Mat m, int n)
{
  double nrm = 0.0;
  for (int j = 0; j < n; ++j) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) { s += std::abs(m[static_cast<std::size_t>(i)*n + j]); }
    nrm = std::max(nrm, s);
  }
  int sq = 0;
  while (nrm > 0.25) { nrm *= 0.5; ++sq; }
  const double scale = std::ldexp(1.0, -sq);
  for (auto& x : m) { x *= scale; }
  Mat e(static_cast<std::size_t>(n)*n, 0.0), t(e);
  for (int i = 0; i < n; ++i) { e[static_cast<std::size_t>(i)*n + i] = 1.0; t[static_cast<std::size_t>(i)*n + i] = 1.0; }
  for (int k = 1; k <= 16; ++k) {
    t = matmul(t, m, n);
    for (auto& x : t) { x /= double(k); }
    for (std::size_t q = 0; q < e.size(); ++q) { e[q] += t[q]; }
  }
  for (int s = 0; s < sq; ++s) { e = matmul(e, e, n); }
  return e;
}
// Spectral norm by power iteration on M^T M
double norm2 (const Mat& m, int n)
{
  std::vector<double> x(static_cast<std::size_t>(n), 1.0), y(x);
  double lam = 0.0;
  for (int it = 0; it < 200; ++it) {
    for (int i = 0; i < n; ++i) {
      double s = 0.0;
      for (int j = 0; j < n; ++j) { s += m[static_cast<std::size_t>(i)*n + j]*x[static_cast<std::size_t>(j)]; }
      y[static_cast<std::size_t>(i)] = s;
    }
    double nx2 = 0.0;
    for (int j = 0; j < n; ++j) {
      double s = 0.0;
      for (int i = 0; i < n; ++i) { s += m[static_cast<std::size_t>(i)*n + j]*y[static_cast<std::size_t>(i)]; }
      x[static_cast<std::size_t>(j)] = s; nx2 += s*s;
    }
    const double nxx = std::sqrt(nx2);
    if (!(nxx > 0.0) || !std::isfinite(nxx)) { return nxx; }
    for (auto& v : x) { v /= nxx; }
    lam = nxx;
  }
  return std::sqrt(lam);
}

} // namespace

// Motivation (Pressel's review of #4231): on terrain whose slope varies, the transpose argument
// of the uniform slope no longer holds, and with K_h >> K_v the symmetric part of the operator
// has directions that gain energy.  What must hold is that no mode grows: the momentum-diffusion
// operator of the actual kernels, assembled on steep curved terrain (h dx/dz about 20, a
// four-cell sine) with K_h/K_v = 1000 and K_h alternating cell by cell, must give a bounded
// transient and decay.  Development's K_v projection has growing modes here (largest real
// parts of +0.06 to +0.27 1/s in an eigenvalue analysis of such matrices, against decay rates
// of up to 0.5 1/s).
TEST(TerrainStress, CurvedSlopeOperatorHasNoGrowingMode)
{
  const int nx = 12, ny = 2, nz = 12, m = 2;
  const Real dx = Real(3000.0), dz = Real(50.0);
  TerrainStressCase c(nx, ny, nz, dx, dx, dz, Real(0.0), Real(0.0), Real(1.e5), Real(100.0), Real(0.0));
  c.init();
  c.set_variable_kh(1);
  {
    auto znd = c.z_nd.array();
    const Real ldx = dx, ldz = dz;
    ParallelFor(c.z_nd.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      const Real x = Real(i)*ldx;
      znd(i,j,k) = Real(k)*ldz + Real(800.0)*std::sin(Real(2.0)*Real(3.14159265358979323846)*x/(Real(4.0)*ldx));
    });
    Gpu::streamSynchronize();
  }
  // y-independent fields on the interior faces, zero within m cells of the x and z ends
  struct Dof { int c, i, k; };
  std::vector<Dof> dofs;
  for (int i = m; i <= nx - m; ++i) { for (int k = m; k < nz - m; ++k) { dofs.push_back({0,i,k}); } }
  for (int i = m; i <  nx - m; ++i) { for (int k = m; k < nz - m; ++k) { dofs.push_back({1,i,k}); } }
  for (int i = m; i <  nx - m; ++i) { for (int k = m; k <= nz - m; ++k) { dofs.push_back({2,i,k}); } }
  const int n = static_cast<int>(dofs.size());
  Mat A(static_cast<std::size_t>(n)*n, 0.0);
  FArrayBox hu(c.ubx,1,The_Pinned_Arena()), hv(c.vbx,1,The_Pinned_Arena()), hw(c.wbx,1,The_Pinned_Arena());
  for (int d = 0; d < n; ++d) {
    c.u.setVal<RunOn::Device>(Real(0.0)); c.v.setVal<RunOn::Device>(Real(0.0)); c.w.setVal<RunOn::Device>(Real(0.0));
    FArrayBox& f = (dofs[d].c == 0) ? c.u : (dofs[d].c == 1) ? c.v : c.w;
    const Box fb = f.box();
    f.setVal<RunOn::Device>(Real(1.0), Box(IntVect(dofs[d].i, fb.smallEnd(1), dofs[d].k),
                                           IntVect(dofs[d].i, fb.bigEnd(1),   dofs[d].k), fb.ixType()));
    c.compute(true);
    copy_to_host(c.rhs_u, hu); copy_to_host(c.rhs_v, hv); copy_to_host(c.rhs_w, hw);
    const auto ru = hu.const_array(), rv = hv.const_array(), rw = hw.const_array();
    for (int r = 0; r < n; ++r) {
      const auto& q = dofs[r];
      A[static_cast<std::size_t>(r)*n + d] =
        double((q.c == 0) ? ru(q.i,0,q.k) : (q.c == 1) ? rv(q.i,0,q.k) : rw(q.i,0,q.k));
    }
  }
  // ||exp(A t)|| at t = t0 * 2^s: bounded transient, then decay
  double anorm = 0.0;
  for (double x : A) { anorm = std::max(anorm, std::abs(x)); }
  ASSERT_GT(anorm, 0.0);
  const double t0 = 0.01 / (double(n) * anorm);
  Mat E = A;
  for (auto& x : E) { x *= t0; }
  E = expm(E, n);
  double peak = 0.0, last = 0.0;
  for (int s = 0; s < 34; ++s) {
    last = norm2(E, n);
    ASSERT_TRUE(std::isfinite(last)) << "exp(A t) blew up at t = t0 * 2^" << s;
    peak = std::max(peak, last);
    E = matmul(E, E, n);
  }
  EXPECT_LT(peak, 3.0) << "transient growth of exp(A t)";
  EXPECT_LT(last, 1.e-3) << "exp(A t) does not decay";
}

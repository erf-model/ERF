#include <AMReX_BCRec.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Gpu.H>

#include <ERF_DataStruct.H>
#include <ERF_Diffusion.H>
#include <ERF_IndexDefines.H>
#include <ERF_WallGradientCoeffs.H>

#include <gtest/gtest.h>

#include <cmath>
#include <utility>
#include <vector>

using namespace amrex;

// The implicit vertical momentum solve closes a no-slip (Dirichlet) z wall with the
// one-sided, second-order gradient the explicit strain uses there (ComputeStrain_N/S/T).
// The solve takes the explicit wall stress out and puts its own back, so a different wall
// stencil (it used 2 (u_0 - u_wall) / dz, first order) shifts the steady state of a
// wall-bounded flow: PoiseuilleFlow drifted at the wall cells under the implicit solve.

namespace {

Real
tol_for_scale (const Real scale)
{
  return scale * (sizeof(Real) == 8 ? Real(1.e-12) : Real(2.e-5));
}

std::vector<double>
thomas (std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> r)
{
  const int n = static_cast<int>(b.size());
  for (int k = 1; k < n; ++k) {
    const double m = a[k] / b[k-1];
    b[k] -= m * c[k-1];
    r[k] -= m * r[k-1];
  }
  std::vector<double> x(n);
  x[n-1] = r[n-1] / b[n-1];
  for (int k = n-2; k >= 0; --k) { x[k] = (r[k] - c[k] * x[k+1]) / b[k]; }
  return x;
}

// d/dz at z = 0 of the quadratic through (0, f_w), (z0, f_0), (z1, f_1): the weights of
// f_w, f_0 and f_1, from the derivative of the Lagrange basis
void
lagrange_slope_at_wall (double z0, double z1, double& w_w, double& w_0, double& w_1)
{
  w_w = -(z0 + z1) / (z0 * z1);
  w_0 =  z1 / (z0 * (z1 - z0));
  w_1 = -z0 / (z1 * (z1 - z0));
}

constexpr int NZ = 8;
constexpr Real DZ   = Real(0.25);   // nominal (computational) spacing
constexpr Real RHO  = Real(1.2);
constexpr Real MU   = Real(0.7);    // constant dynamic viscosity: rhoAlpha = MU on every face
constexpr Real DT   = Real(0.3);
constexpr Real UW_LO = Real(0.4), UW_HI = Real(-0.3);  // wall velocities

// Set the column extent of every cell to the whole column, 0 to NZ-1
void
fill_col_kext (IArrayBox& col_kext)
{
  auto ck = col_kext.array();
  ParallelFor(col_kext.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
    amrex::ignore_unused(i, j, k);
    ck(i,j,k,0) = 0;
    ck(i,j,k,1) = NZ-1;
  });
}

// Flat terrain with the given levels, extended linearly into the ghost nodes; the cell
// Jacobian is the cell height over the nominal spacing
void
fill_terrain (const std::vector<Real>& znd, FArrayBox& z_nd, FArrayBox& detJ)
{
  Gpu::DeviceVector<Real> zv(znd.size());
  Gpu::copy(Gpu::hostToDevice, znd.begin(), znd.end(), zv.begin());
  const Real* zp = zv.data();
  const Real dzb = znd[1] - znd[0], dzt = znd[NZ] - znd[NZ-1], ztop = znd[NZ];
  auto zn = z_nd.array();
  ParallelFor(z_nd.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
    amrex::ignore_unused(i, j);
    zn(i,j,k) = (k < 0) ? Real(k)*dzb : ((k > NZ) ? ztop + Real(k-NZ)*dzt : zp[k]);
  });
  auto dj = detJ.array();
  ParallelFor(detJ.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
    dj(i,j,k) = (zn(i,j,k+1) - zn(i,j,k)) / DZ;
  });
  // zv is freed on return
  Gpu::streamSynchronize();
}

// Two cell columns in x (the u face between them is solved), one in y, NZ cells in z
struct WallColumn
{
  Box cells{IntVect(0), IntVect(1, 0, NZ-1)};
  Box domain{IntVect(0), IntVect(1, 0, NZ-1)};
  Box col{IntVect(1, 0, 0), IntVect(1, 0, NZ-1)};
  std::vector<Real> znd;                     // z of the NZ+1 nodes, 0 at the bottom wall
  FArrayBox cons, rho_u, zero_tau, z_nd, detJ, mf;
  IArrayBox col_kext;
  std::vector<BCRec> bcs;
  SolverChoice solver;
  GpuArray<Real, AMREX_SPACEDIM> dxInv{Real(1.0), Real(1.0), Real(1.0)/DZ};

  explicit WallColumn (const std::vector<Real>& z_nodes) : znd(z_nodes)
  {
    const Box g = grow(cells, 2);
    cons.resize(g, 2);
    rho_u.resize(grow(surroundingNodes(cells, 0), 2), 1);
    zero_tau.resize(grow(convert(cells, IntVect(1,0,1)), 2), 1);
    z_nd.resize(grow(convert(cells, IntVect(1)), 2), 1);
    detJ.resize(g, 1);
    mf.resize(grow(surroundingNodes(cells, 0), 2), 1);
    col_kext.resize(g, 2);

    cons.setVal<RunOn::Device>(RHO, cons.box(), Rho_comp, 1);
    cons.setVal<RunOn::Device>(RHO * Real(300.0), cons.box(), RhoTheta_comp, 1);
    zero_tau.setVal<RunOn::Device>(Real(0.0));
    mf.setVal<RunOn::Device>(Real(1.0));
    fill_col_kext(col_kext);
    fill_terrain(znd, z_nd, detJ);

    // No-slip z walls with the wall velocity in the ghost faces
    bcs.resize(BCVars::NumTypes);
    for (auto& bc : bcs) {
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        bc.setLo(d, (d < 2) ? ERFBCType::int_dir : ERFBCType::ext_dir);
        bc.setHi(d, (d < 2) ? ERFBCType::int_dir : ERFBCType::ext_dir);
      }
    }
    solver.diffChoice.molec_diff_type = MolecDiffType::Constant;
    solver.diffChoice.dynamic_viscosity = MU;
    solver.turbChoice.resize(1);
    solver.turbChoice[0].use_kturb = false;
  }

  double zc (int k) const { return 0.5 * (double(znd[k]) + double(znd[k+1])); }
  double dz (int k) const { return double(znd[k+1]) - double(znd[k]); }

  // rho_u = RHO * u_k + DT * src_k in the column, and RHO * the wall velocity in the ghosts
  void set_state (const std::vector<double>& u, const std::vector<double>& src)
  {
    FArrayBox h(rho_u.box(), 1, The_Pinned_Arena());
    auto ha = h.array();
    for (int k = rho_u.box().smallEnd(2); k <= rho_u.box().bigEnd(2); ++k) {
      const double v = (k < 0) ? RHO * double(UW_LO)
                     : (k > NZ-1) ? RHO * double(UW_HI)
                     : RHO * u[k] + double(DT) * src[k];
      for (int i = rho_u.box().smallEnd(0); i <= rho_u.box().bigEnd(0); ++i) {
        for (int j = rho_u.box().smallEnd(1); j <= rho_u.box().bigEnd(1); ++j) {
          ha(i,j,k) = Real(v);
        }
      }
    }
    Gpu::copy(Gpu::hostToDevice, h.dataPtr(0), h.dataPtr(0) + h.size(), rho_u.dataPtr(0));
    Gpu::streamSynchronize();
  }

  std::vector<double> velocity () const
  {
    FArrayBox h(rho_u.box(), 1, The_Pinned_Arena());
    Gpu::copy(Gpu::deviceToHost, rho_u.dataPtr(0), rho_u.dataPtr(0) + rho_u.size(), h.dataPtr(0));
    Gpu::streamSynchronize();
    const auto ha = h.const_array();
    std::vector<double> u(NZ);
    for (int k = 0; k < NZ; ++k) { u[k] = double(ha(1,0,k)) / double(RHO); }
    return u;
  }

  void solve_N ()
  {
    ImplicitDiffForMomLU_N<0>(col, domain, 0, double(DT), col_kext.const_array(),
                              cons.const_array(), rho_u.array(),
                              zero_tau.const_array(), zero_tau.const_array(), dxInv,
                              Array4<const Real>{}, solver, bcs.data(), false, Real(1.0), false);
    Gpu::streamSynchronize();
  }

  void solve_S ()
  {
    std::vector<Real> dzh(NZ);
    for (int k = 0; k < NZ; ++k) { dzh[k] = znd[k+1] - znd[k]; }
    Gpu::DeviceVector<Real> dzd(NZ);
    Gpu::copy(Gpu::hostToDevice, dzh.begin(), dzh.end(), dzd.begin());
    ImplicitDiffForMomLU_S<0>(col, domain, 0, double(DT), col_kext.const_array(),
                              cons.const_array(), rho_u.array(),
                              zero_tau.const_array(), zero_tau.const_array(), dzd,
                              Array4<const Real>{}, solver, bcs.data(), false, Real(1.0), false);
    Gpu::streamSynchronize();
  }

  void solve_T ()
  {
    ImplicitDiffForMomLU_T<0>(col, domain, 0, double(DT), col_kext.const_array(),
                              cons.const_array(), rho_u.array(),
                              zero_tau.const_array(), zero_tau.const_array(),
                              z_nd.const_array(), detJ.const_array(), dxInv,
                              Array4<const Real>{}, solver, bcs.data(), false, Real(1.0), false,
                              mf.const_array(), mf.const_array(), false);
    Gpu::streamSynchronize();
  }

  // The tridiagonal system the solve should build, written out independently: rho u_k -
  // DT/dz_k (F_{k+1/2} - F_{k-1/2}) = rhs_k with F = MU du/dz, centred between the cell
  // centres inside the column and from the Lagrange quadratic through the wall value and
  // the first two cells at each wall
  std::vector<double> reference (const std::vector<double>& u, const std::vector<double>& src) const
  {
    const double mu = double(MU), dt = double(DT), rho = double(RHO);
    std::vector<double> a(NZ, 0.0), b(NZ, rho), c(NZ, 0.0), r(NZ);
    for (int k = 0; k < NZ; ++k) {
      r[k] = rho * u[k] + dt * src[k];
      const double f = dt * mu / dz(k);
      if (k > 0) {
        const double g = f / (zc(k) - zc(k-1));
        a[k] -= g; b[k] += g;
      }
      if (k < NZ-1) {
        const double g = f / (zc(k+1) - zc(k));
        c[k] -= g; b[k] += g;
      }
    }
    double ww, w0, w1;
    // Bottom wall: F_{-1/2} = MU (ww u_wall + w0 u_0 + w1 u_1), distances up from the wall
    lagrange_slope_at_wall(zc(0) - double(znd[0]), zc(1) - double(znd[0]), ww, w0, w1);
    {
      const double f = dt * mu / dz(0);
      b[0] += f * w0;  c[0] += f * w1;  r[0] -= f * ww * double(UW_LO);
    }
    // Top wall: d/dz = -d/ds with s measured down from the wall
    lagrange_slope_at_wall(double(znd[NZ]) - zc(NZ-1), double(znd[NZ]) - zc(NZ-2), ww, w0, w1);
    {
      const double f = dt * mu / dz(NZ-1);
      b[NZ-1] += f * w0;  a[NZ-1] += f * w1;  r[NZ-1] -= f * ww * double(UW_HI);
    }
    return thomas(a, b, c, r);
  }
};

std::vector<Real>
uniform_levels ()
{
  std::vector<Real> z(NZ+1);
  for (int k = 0; k <= NZ; ++k) { z[k] = Real(k) * DZ; }
  return z;
}

std::vector<Real>
stretched_levels ()
{
  std::vector<Real> z(NZ+1);
  Real h = Real(0.1);
  z[0] = Real(0.0);
  for (int k = 1; k <= NZ; ++k) { z[k] = z[k-1] + h; h *= Real(1.3); }
  return z;
}

} // namespace

// Motivation: on a uniform mesh the coefficients must be the ComputeStrain_N stencil
// (-8/3, 3, -1/3), and on any mesh they must be the slope at the wall of the quadratic
// through the wall value and the first two cell centres.
TEST(ImplicitWallStencil, CoefficientsAreTheSecondOrderOneSidedGradient)
{
  const WallGradCoeffs u = WallGradientCoeffs(Real(1.0), Real(1.0));
  EXPECT_NEAR(u.c_wall, -Real(8.0)/Real(3.0), tol_for_scale(Real(1.0)));
  EXPECT_NEAR(u.c_0,     Real(3.0),           tol_for_scale(Real(1.0)));
  EXPECT_NEAR(u.c_1,    -Real(1.0)/Real(3.0), tol_for_scale(Real(1.0)));
  const WallGradCoeffs n = WallGradientCoeffsUniform();
  EXPECT_NEAR(n.c_wall, u.c_wall, tol_for_scale(Real(1.0)));
  EXPECT_NEAR(n.c_0,    u.c_0,    tol_for_scale(Real(1.0)));
  EXPECT_NEAR(n.c_1,    u.c_1,    tol_for_scale(Real(1.0)));

  const std::pair<double,double> spacings[] = {{0.1, 0.13}, {2.0, 0.5}, {0.3, 0.9}};
  for (const auto& d : spacings) {
    const WallGradCoeffs s = WallGradientCoeffs(Real(d.first), Real(d.second));
    double ww, w0, w1;
    lagrange_slope_at_wall(0.5*d.first, d.first + 0.5*d.second, ww, w0, w1);
    // The helper's coefficients are per 1/dz0
    EXPECT_NEAR(double(s.c_wall), ww * d.first, double(tol_for_scale(Real(10.0))));
    EXPECT_NEAR(double(s.c_0),    w0 * d.first, double(tol_for_scale(Real(10.0))));
    EXPECT_NEAR(double(s.c_1),    w1 * d.first, double(tol_for_scale(Real(10.0))));
    EXPECT_NEAR(double(s.c_wall + s.c_0 + s.c_1), 0.0, double(tol_for_scale(Real(10.0))));
  }
}

// Motivation: Poiseuille flow.  A parabola with the walls at its end values, forced by the
// constant source that balances its diffusion, is a steady state of the explicit operator on
// a uniform mesh; one implicit step must return it unchanged.  With the first-order wall row
// the wall cells moved by 1.7e-2 here (and the stretched and terrain solves missed the
// reference below by up to 1.8e-2).
TEST(ImplicitWallStencil, UniformSolveKeepsTheParabolaSteady)
{
  WallColumn c(uniform_levels());
  const double H = double(NZ) * double(DZ), G = 1.5;
  const double lo = double(UW_LO), hi = double(UW_HI);
  std::vector<double> u(NZ), src(NZ);
  for (int k = 0; k < NZ; ++k) {
    const double z = c.zc(k);
    u[k]   = lo + (hi - lo) * z / H + 4.0 * G * z * (H - z) / (H * H);
    src[k] = double(MU) * 8.0 * G / (H * H);  // -MU u''
  }
  c.set_state(u, src);
  c.solve_N();
  const auto un = c.velocity();
  for (int k = 0; k < NZ; ++k) {
    EXPECT_NEAR(un[k], u[k], double(tol_for_scale(Real(2.0)))) << "k " << k;
  }
}

// Motivation: on a stretched column the stretched and terrain solves must both build the
// tridiagonal system with the second-order one-sided wall gradient (written out here from
// the Lagrange quadratic), with moving walls at both ends.
TEST(ImplicitWallStencil, StretchedAndTerrainSolvesMatchTheReference)
{
  std::vector<double> u(NZ), src(NZ);
  for (int k = 0; k < NZ; ++k) {
    u[k]   = 0.2 + 0.1 * double(k) - 0.03 * double(k*k);
    src[k] = 0.05 * double(k);
  }
  for (const bool terrain : {false, true}) {
    WallColumn c(stretched_levels());
    c.set_state(u, src);
    if (terrain) { c.solve_T(); } else { c.solve_S(); }
    const auto un  = c.velocity();
    const auto ref = c.reference(u, src);
    for (int k = 0; k < NZ; ++k) {
      EXPECT_NEAR(un[k], ref[k], double(tol_for_scale(Real(2.0))))
          << (terrain ? "terrain" : "stretched") << " k " << k;
    }
  }
}

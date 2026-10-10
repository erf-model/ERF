#include <AMReX_BCRec.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Gpu.H>

#include <ERF_DataStruct.H>
#include <ERF_Diffusion.H>
#include <ERF_IndexDefines.H>
#include <ERF_ScalarDiffusion.H>
#include <ERF_TerrainDiffusionLimits.H>
#include <ERF_TerrainImplicitMetric.H>

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using namespace amrex;

// erf.implicit_terrain_metric: on a terrain-fitted mesh the projected horizontal fluxes add
// a vertical-diffusion term K_h * M * d/dz to every zeta flux.  With the option on, the
// implicit vertical solves take its compact form on the interior faces, and the explicit
// fluxes give up the same compact term.  These tests check each side against the closed form
// on a mesh tilted by constant slopes, where h_zeta = 1 and every metric is exact.

namespace {

Real
tol_for_scale (const Real scale)
{
  return scale * (sizeof(Real) == 8 ? Real(1.e-12) : Real(2.e-5));
}

void
copy_to_host (const FArrayBox& src, FArrayBox& dst)
{
  Gpu::copy(Gpu::deviceToHost, src.dataPtr(0), src.dataPtr(0) + src.size(), dst.dataPtr(0));
  Gpu::streamSynchronize();
}

// Solve A x = r for a tridiagonal A (a: sub, b: diag, c: super), in double
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

constexpr int NZ = 8;
constexpr Real DX = Real(300.0), DY = Real(200.0), DZ = Real(40.0);
constexpr Real SX = Real(0.35), SY = Real(-0.25);    // terrain slopes
constexpr Real MX = Real(1.1), MY = Real(0.9);       // map factors
constexpr Real KH = Real(80.0), KV = Real(2.0);      // rho*K, horizontal and vertical
constexpr Real RHO = Real(1.2);
constexpr Real DT = Real(300.0);

Real phi0 (int k) { return Real(300.0) + Real(0.4)*Real(k) + Real(0.05)*Real(k*k); }
Real u0   (int k) { return Real(5.0) + Real(1.5)*std::sqrt(Real(k) + Real(0.5)); }

// Two columns of cells in x (the u face between them is solved), one in y, NZ in z, with
// two ghost cells everywhere
struct ColumnCase
{
  Box cells{IntVect(0), IntVect(1, 0, NZ-1)};
  Box domain{IntVect(0), IntVect(1, 0, NZ-1)};
  FArrayBox cons, prim, mu, z_nd, detJ, mf_m_x, mf_m_y, mf_u_x, mf_u_y, rho_u, zero_tau;
  IArrayBox col_kext;
  std::vector<BCRec> bcs;
  SolverChoice solver;
  GpuArray<Real, AMREX_SPACEDIM> dxInv{Real(1.0)/DX, Real(1.0)/DY, Real(1.0)/DZ};

  void init ()
  {
    const Box g = grow(cells, 2);
    cons.resize(g, 2);
    prim.resize(g, 1);
    mu.resize(g, EddyDiff::NumDiffs);
    z_nd.resize(grow(convert(cells, IntVect(1)), 2), 1);
    detJ.resize(g, 1);
    mf_m_x.resize(g, 1); mf_m_y.resize(g, 1);
    mf_u_x.resize(grow(surroundingNodes(cells, 0), 2), 1);
    mf_u_y.resize(grow(surroundingNodes(cells, 0), 2), 1);
    rho_u.resize(grow(surroundingNodes(cells, 0), 2), 1);
    zero_tau.resize(grow(convert(cells, IntVect(1,0,1)), 2), 1);
    col_kext.resize(g, 2);

    auto c = cons.array(); auto p = prim.array(); auto zn = z_nd.array();
    auto ru = rho_u.array(); auto ck = col_kext.array();
    ParallelFor(g, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      const int kk = (k < 0) ? 0 : ((k > NZ-1) ? NZ-1 : k);  // clamp to the column
      const Real ph = Real(300.0) + Real(0.4)*Real(kk) + Real(0.05)*Real(kk*kk);
      c(i,j,k,Rho_comp) = RHO;
      c(i,j,k,RhoTheta_comp) = RHO * ph;
      p(i,j,k,0) = ph;
      ck(i,j,k,0) = 0;
      ck(i,j,k,1) = NZ-1;
    });
    ParallelFor(z_nd.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      zn(i,j,k) = Real(k)*DZ + SX*Real(i)*DX + SY*Real(j)*DY;
    });
    ParallelFor(rho_u.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      const int kk = (k < 0) ? 0 : ((k > NZ-1) ? NZ-1 : k);  // clamp to the column
      ru(i,j,k) = RHO * (Real(5.0) + Real(1.5)*std::sqrt(Real(kk) + Real(0.5)));
    });
    mu.setVal<RunOn::Device>(Real(0.0));
    for (int n : {EddyDiff::Mom_h, EddyDiff::Theta_h}) { mu.setVal<RunOn::Device>(KH, mu.box(), n, 1); }
    for (int n : {EddyDiff::Mom_v, EddyDiff::Theta_v}) { mu.setVal<RunOn::Device>(KV, mu.box(), n, 1); }
    detJ.setVal<RunOn::Device>(Real(1.0));
    mf_m_x.setVal<RunOn::Device>(MX); mf_m_y.setVal<RunOn::Device>(MY);
    mf_u_x.setVal<RunOn::Device>(MX); mf_u_y.setVal<RunOn::Device>(MY);
    zero_tau.setVal<RunOn::Device>(Real(0.0));
    Gpu::streamSynchronize();

    // Zero-gradient top and bottom: no boundary flux, so only the interior faces carry
    // coefficients and the bottom/top faces must stay without the metric term
    bcs.resize(BCVars::NumTypes);
    for (auto& bc : bcs) {
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        bc.setLo(d, (d < 2) ? ERFBCType::int_dir : ERFBCType::foextrap);
        bc.setHi(d, (d < 2) ? ERFBCType::int_dir : ERFBCType::foextrap);
      }
    }
    solver.diffChoice.molec_diff_type = MolecDiffType::None;
    solver.diffChoice.dynamic_viscosity = Real(0.0);
    solver.diffChoice.rhoAlpha_T = Real(0.0);
    solver.turbChoice.resize(1);
    solver.turbChoice[0].use_kturb = true;
  }

  void solve_theta (bool metric)
  {
    GpuArray<Real, AMREX_SPACEDIM*2> neumann{};
    Box col(IntVect(1, 0, 0), IntVect(1, 0, NZ-1));
    ImplicitDiffForStateLU_T(col, domain, 0, RhoTheta_comp, double(DT), neumann, cons.array(),
                             z_nd.const_array(), detJ.const_array(), dxInv,
                             Array4<const Real>{}, mu.const_array(), solver, bcs.data(),
                             false, Real(1.0), false,
                             mf_m_x.const_array(), mf_m_y.const_array(), metric);
    Gpu::streamSynchronize();
  }

  void solve_u (bool metric)
  {
    Box col(IntVect(1, 0, 0), IntVect(1, 0, NZ-1));  // the u face between the two columns
    ImplicitDiffForMomLU_T<0>(col, domain, 0, double(DT), col_kext.const_array(),
                              cons.const_array(), rho_u.array(),
                              zero_tau.const_array(), zero_tau.const_array(),
                              z_nd.const_array(), detJ.const_array(), dxInv,
                              mu.const_array(), solver, bcs.data(), false, Real(1.0), false,
                              mf_u_x.const_array(), mf_u_y.const_array(), metric);
    Gpu::streamSynchronize();
  }

  // The v face at j = 0 between cells j = -1 and j = 0 (the ghost cells hold the same column
  // data), with the u-face map factors (constant) standing in for the v-face ones
  void solve_v (bool metric)
  {
    Box col(IntVect(1, 0, 0), IntVect(1, 0, NZ-1));
    ImplicitDiffForMomLU_T<1>(col, domain, 0, double(DT), col_kext.const_array(),
                              cons.const_array(), rho_u.array(),
                              zero_tau.const_array(), zero_tau.const_array(),
                              z_nd.const_array(), detJ.const_array(), dxInv,
                              mu.const_array(), solver, bcs.data(), false, Real(1.0), false,
                              mf_u_x.const_array(), mf_u_y.const_array(), metric);
    Gpu::streamSynchronize();
  }

  void explicit_split (FArrayBox& zflux, Real implicit_fac)
  {
    Box col(IntVect(1, 0, 0), IntVect(1, 0, NZ-1));
    zflux.resize(surroundingNodes(col, 2), 1);
    zflux.setVal<RunOn::Device>(Real(0.0));
    ScalarDiffusionFieldViews field;
    field.scalar = prim.const_array();
    field.scalar_comp = 0;
    field.density = cons.const_array();
    field.rho_comp = Rho_comp;
    field.mu_turb = mu.const_array();
    field.zflux = zflux.array();
    ScalarDiffusionFluxPolicy policy;
    policy.coefficients.eddy_h_comp = EddyDiff::Theta_h;
    policy.coefficients.eddy_v_comp = EddyDiff::Theta_v;
    policy.coefficient_mode = {false, true};
    SubtractTerrainMetricImplicitPart_T(col, domain, field, policy, z_nd.const_array(), dxInv,
                                        mf_m_x.const_array(), mf_m_y.const_array(), implicit_fac);
    Gpu::streamSynchronize();
  }
};

// Reference: (rho - a - c) x_k + a x_{k-1} + c x_{k+1} = rho x_old_k with face coefficients
// C_f (zero on the bottom and top faces), a = -dt/dz^2 C_k, c = -dt/dz^2 C_{k+1}
std::vector<double>
reference_solve (const std::vector<double>& old, double C_interior)
{
  const double f = double(DT) / (double(DZ) * double(DZ));
  std::vector<double> a(NZ), b(NZ), c(NZ), r(NZ);
  for (int k = 0; k < NZ; ++k) {
    a[k] = (k > 0)      ? -f * C_interior : 0.0;
    c[k] = (k < NZ - 1) ? -f * C_interior : 0.0;
    b[k] = double(RHO) - a[k] - c[k];
    r[k] = double(RHO) * old[k];
  }
  return thomas(a, b, c, r);
}

} // namespace

// Motivation: the factors that multiply K_h are the h*d/dz parts of the projected
// horizontal fluxes: (4/3) h_xi^2 + h_eta^2 for u (from h_xi (S11 - er/3) + h_eta S12; the
// expansion rate carries -h_xi du/dz too), h_xi^2 + (4/3) h_eta^2 for v, h_xi^2 + h_eta^2 for
// scalars, each slope scaled by its map factor; and only the faces strictly inside the domain
// are split.
TEST(TerrainImplicitMetric, FactorsAndFacesByHand)
{
  const Real ax = SX*MX, ay = SY*MY;
  const Real four_thirds = Real(4.0)/Real(3.0);
  EXPECT_NEAR(TerrainMetricMomFactor(0, SX, SY, MX, MY), four_thirds*ax*ax + ay*ay, tol_for_scale(Real(1.0)));
  EXPECT_NEAR(TerrainMetricMomFactor(1, SX, SY, MX, MY), ax*ax + four_thirds*ay*ay, tol_for_scale(Real(1.0)));
  EXPECT_NEAR(TerrainMetricScalarFactor(SX, SY, MX, MY), ax*ax + ay*ay, tol_for_scale(Real(1.0)));
  // 0.35*1.1 = 0.385, -0.25*0.9 = -0.225: (4/3)*0.148225 + 0.050625 = 0.2482583...,
  // and 0.148225 + (4/3)*0.050625 = 0.215725 for v
  EXPECT_NEAR(TerrainMetricMomFactor(0, SX, SY, MX, MY), Real(0.24825833333333333), tol_for_scale(Real(1.0)));
  EXPECT_NEAR(TerrainMetricMomFactor(1, SX, SY, MX, MY), Real(0.215725), tol_for_scale(Real(1.0)));
  EXPECT_FALSE(TerrainMetricImplicitFace(0, 0, NZ-1));
  EXPECT_TRUE (TerrainMetricImplicitFace(1, 0, NZ-1));
  EXPECT_TRUE (TerrainMetricImplicitFace(NZ-1, 0, NZ-1));
  EXPECT_FALSE(TerrainMetricImplicitFace(NZ, 0, NZ-1));
}

// Motivation: with the option on, the theta solve must use K_v + K_h M_s on every interior
// face (and nothing on the bottom and top faces); with it off, K_v alone, as before.  The
// two answers differ by far more than the tolerance, so a solve that ignored the flag fails.
TEST(TerrainImplicitMetric, ScalarColumnSolveMatchesTridiagonal)
{
  std::vector<double> old(NZ);
  for (int k = 0; k < NZ; ++k) { old[k] = double(phi0(k)); }
  const double Ms = double(TerrainMetricScalarFactor(SX, SY, MX, MY));

  for (const bool metric : {false, true}) {
    ColumnCase c;
    c.init();
    c.solve_theta(metric);
    FArrayBox h(c.cons.box(), 2, The_Pinned_Arena());
    copy_to_host(c.cons, h);
    const auto ha = h.const_array();
    const auto ref = reference_solve(old, double(KV) + (metric ? double(KH)*Ms : 0.0));
    const auto ref_off = reference_solve(old, double(KV));
    double diff_on_off = 0.0;
    for (int k = 0; k < NZ; ++k) {
      const Real phi = ha(1,0,k,RhoTheta_comp) / ha(1,0,k,Rho_comp);
      EXPECT_NEAR(phi, Real(ref[k]), tol_for_scale(Real(300.0))) << "metric " << metric << " k " << k;
      diff_on_off = std::max(diff_on_off, std::abs(ref[k] - ref_off[k]));
    }
    // The option must matter by much more than the tolerance (148x in single precision)
    if (metric) { EXPECT_GT(diff_on_off, 50.0 * double(tol_for_scale(Real(300.0)))); }
  }
}

// Motivation: the same for the u solve, with K_v + K_h M_u on the interior faces.
TEST(TerrainImplicitMetric, MomentumColumnSolveMatchesTridiagonal)
{
  std::vector<double> old(NZ);
  for (int k = 0; k < NZ; ++k) { old[k] = double(u0(k)); }
  const double Mu = double(TerrainMetricMomFactor(0, SX, SY, MX, MY));

  for (const bool metric : {false, true}) {
    ColumnCase c;
    c.init();
    c.solve_u(metric);
    FArrayBox h(c.rho_u.box(), 1, The_Pinned_Arena());
    copy_to_host(c.rho_u, h);
    const auto ha = h.const_array();
    const auto ref = reference_solve(old, double(KV) + (metric ? double(KH)*Mu : 0.0));
    const auto ref_off = reference_solve(old, double(KV));
    double diff_on_off = 0.0;
    for (int k = 0; k < NZ; ++k) {
      EXPECT_NEAR(ha(1,0,k) / RHO, Real(ref[k]), tol_for_scale(Real(10.0))) << "metric " << metric << " k " << k;
      diff_on_off = std::max(diff_on_off, std::abs(ref[k] - ref_off[k]));
    }
    if (metric) { EXPECT_GT(diff_on_off, 50.0 * double(tol_for_scale(Real(10.0)))); }
  }
}

// Motivation: the v solve (stagdir 1) must use the v weights, K_v + K_h M_v with
// M_v = (h_xi mx)^2 + (4/3) (h_eta my)^2, which differ from the u weights on this mesh.
TEST(TerrainImplicitMetric, VMomentumColumnSolveMatchesTridiagonal)
{
  std::vector<double> old(NZ);
  for (int k = 0; k < NZ; ++k) { old[k] = double(u0(k)); }
  const double Mv = double(TerrainMetricMomFactor(1, SX, SY, MX, MY));
  const double Mu = double(TerrainMetricMomFactor(0, SX, SY, MX, MY));
  ASSERT_GT(std::abs(Mu - Mv), 0.02);  // 0.0325 on this mesh
  for (const bool metric : {false, true}) {
    ColumnCase c;
    c.init();
    c.solve_v(metric);
    FArrayBox h(c.rho_u.box(), 1, The_Pinned_Arena());
    copy_to_host(c.rho_u, h);
    const auto ha = h.const_array();
    const auto ref = reference_solve(old, double(KV) + (metric ? double(KH)*Mv : 0.0));
    const auto ref_u = reference_solve(old, double(KV) + double(KH)*Mu);
    double diff_uv = 0.0;
    for (int k = 0; k < NZ; ++k) {
      EXPECT_NEAR(ha(1,0,k) / RHO, Real(ref[k]), tol_for_scale(Real(10.0))) << "metric " << metric << " k " << k;
      diff_uv = std::max(diff_uv, std::abs(ref[k] - ref_u[k]));
    }
    // The u and v weights must lead to answers the tolerance can tell apart
    if (metric) { EXPECT_GT(diff_uv, 50.0 * double(tol_for_scale(Real(10.0)))); }
  }
}

// Motivation: the explicit flux must give up exactly what the solve takes, implicit_fac
// times the compact K_h M_s dphi/dz, on the same faces (none on the bottom and top).
TEST(TerrainImplicitMetric, ExplicitSplitRemovesTheCompactTerm)
{
  ColumnCase c;
  c.init();
  FArrayBox zflux;
  const Real f = Real(0.5);
  c.explicit_split(zflux, f);
  FArrayBox h(zflux.box(), 1, The_Pinned_Arena());
  copy_to_host(zflux, h);
  const auto ha = h.const_array();
  const Real Ms = TerrainMetricScalarFactor(SX, SY, MX, MY);
  for (int k = 0; k <= NZ; ++k) {
    const Real expected = (k > 0 && k < NZ) ? f * KH * Ms * (phi0(k) - phi0(k-1)) / DZ : Real(0.0);
    EXPECT_NEAR(ha(1,0,k), expected, tol_for_scale(Real(1.0))) << "face " << k;
  }
  EXPECT_GT(std::abs(ha(1,0,NZ/2)), Real(0.1));
}

// The diffusive time-step check (erf.diffusive_dt_check) counts the K_h metric term
// K_h h^2 d2/dz2 as explicit.  With erf.implicit_terrain_metric its compact form is in the
// implicit solve, so the term must count with the component's explicit fraction e, exactly as
// the vertical term does; without the option, or with e = 1, the rates must not change.
TEST(TerrainImplicitMetric, DiffusiveRateWeightsTheMetricTermByTheExplicitFraction)
{
  // rho = 1, dx = dy = 1 km (a = b = 1e-6), h^2 = 0.01, dz = 50 m (1/dz^2 = 4e-4)
  const Real a = Real(1.e-6), h2 = Real(0.01), dzinv2 = Real(4.e-4);
  // mu_h = 100, mu_v = 10, vertical parts implicit: rows u = 5e-4, w = 4.2e-4 (see
  // Smag2DLimiters.DiffusiveRatesHandValues).  Metric explicit: 2 * 100 * 0.01 * 4e-4 = 8e-4
  // -> 1.3e-3.  Metric implicit: only w's mu_v part, 2 * 10 * 0.01 * 4e-4 = 8e-5 -> 5.8e-4.
  EXPECT_NEAR(MomentumDiffusiveRate(one, Real(100), Real(10), a, a, h2, dzinv2, zero, zero),
              Real(1.3e-3), tol_for_scale(Real(1.3e-3)));
  EXPECT_NEAR(MomentumDiffusiveRate(one, Real(100), Real(10), a, a, h2, dzinv2, zero, zero, zero),
              Real(5.8e-4), tol_for_scale(Real(5.8e-4)));
  // Scalar K_h = 300, K_v = 30, e = 0: 300 * 2e-6 + 300 * 0.01 * 4e-4 = 1.8e-3, and without the
  // metric term 6e-4
  EXPECT_NEAR(ScalarDiffusiveRate(one, Real(300), Real(30), a, a, h2, dzinv2, zero),
              Real(1.8e-3), tol_for_scale(Real(1.8e-3)));
  EXPECT_NEAR(ScalarDiffusiveRate(one, Real(300), Real(30), a, a, h2, dzinv2, zero, zero),
              Real(6.e-4), tol_for_scale(Real(6.e-4)));

  // CellDiffusiveRates on a cell tilted in x (slope 0.2, so h^2 = 0.04), with K_h >> K_v
  const int nk = EddyDiff::NumDiffs;
  std::vector<Real> rho(27, one), K(27*static_cast<std::size_t>(nk), one);
  for (int c = 0; c < nk; ++c) {
    const bool horiz = (c == EddyDiff::Mom_h || c == EddyDiff::Theta_h || c == EddyDiff::KE_h ||
                        c == EddyDiff::Q_h || c == EddyDiff::Scalar_h);
    for (int n = 0; n < 27; ++n) { K[static_cast<std::size_t>(n + 27*c)] = horiz ? Real(200) : Real(2); }
  }
  std::vector<Real> z(64);
  for (int kk = 0; kk < 4; ++kk) { for (int jj = 0; jj < 4; ++jj) { for (int ii = 0; ii < 4; ++ii) {
    z[ii + 4*jj + 16*kk] = Real(0.2) * Real(1000) * Real(ii - 1) + Real(50) * Real(kk - 1);
  }}}
  const Dim3 lo{-1,-1,-1};
  const Array4<const Real> sa(rho.data(), lo, Dim3{2,2,2}, 1);
  const Array4<const Real> ka(K.data(),   lo, Dim3{2,2,2}, nk);
  const Array4<const Real> za(z.data(),   lo, Dim3{3,3,3}, 1);
  const GpuArray<Real,AMREX_SPACEDIM> dxinv{Real(1.e-3), Real(1.e-3), Real(1)/Real(50)};
  const Real hc2 = Real(0.04), dz2 = Real(1)/Real(2500);

  DiffusiveRateSettings set;
  set.variable_dz = true; set.has_q = true; set.has_ke = true; set.has_scalar = true;
  set.e_uv = zero; set.e_w = one; set.e_th = zero; set.e_q = zero; set.e_ke = zero;
  Real rm_off, rs_off, rm_on, rs_on;
  CellDiffusiveRates(0, 0, 0, sa, ka, za, one, one, dxinv, set, rm_off, rs_off);
  set.implicit_metric = true;
  CellDiffusiveRates(0, 0, 0, sa, ka, za, one, one, dxinv, set, rm_on, rs_on);
  EXPECT_NEAR(rm_off, MomentumDiffusiveRate(one, Real(200), Real(2), a, a, hc2, dz2, zero, one),
              tol_for_scale(rm_off));
  EXPECT_NEAR(rm_on,  MomentumDiffusiveRate(one, Real(200), Real(2), a, a, hc2, dz2, zero, one, zero),
              tol_for_scale(rm_on));
  EXPECT_LT(rm_on, Real(0.5) * rm_off);
  // the advected scalar is never in the implicit solve, so it keeps its metric term and sets
  // the scalar rate; with it off, theta's rate (metric dropped) is left
  EXPECT_NEAR(rs_on, ScalarDiffusiveRate(one, Real(200), Real(2), a, a, hc2, dz2, one),
              tol_for_scale(rs_on));
  set.has_scalar = false;
  CellDiffusiveRates(0, 0, 0, sa, ka, za, one, one, dxinv, set, rm_on, rs_on);
  EXPECT_NEAR(rs_on, ScalarDiffusiveRate(one, Real(200), Real(2), a, a, hc2, dz2, zero, zero),
              tol_for_scale(rs_on));
  // with the option on but the vertical parts explicit (e = 1) nothing changes
  set.has_scalar = true;
  set.e_uv = one; set.e_th = one; set.e_q = one; set.e_ke = one;
  Real rm_e1_on, rs_e1_on, rm_e1_off, rs_e1_off;
  CellDiffusiveRates(0, 0, 0, sa, ka, za, one, one, dxinv, set, rm_e1_on, rs_e1_on);
  set.implicit_metric = false;
  CellDiffusiveRates(0, 0, 0, sa, ka, za, one, one, dxinv, set, rm_e1_off, rs_e1_off);
  EXPECT_EQ(rm_e1_on, rm_e1_off);
  EXPECT_EQ(rs_e1_on, rs_e1_off);
}

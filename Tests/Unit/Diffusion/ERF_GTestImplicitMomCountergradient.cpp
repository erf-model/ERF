#include <AMReX_BCRec.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Gpu.H>

#include <ERF_DataStruct.H>
#include <ERF_Diffusion.H>
#include <ERF_IndexDefines.H>

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using namespace amrex;

// The implicit vertical momentum solves add the YSU countergradient flux rho K gamma in
// flux-difference form: row k takes -(F_{k+1/2} - F_{k-1/2}) with F = rho K gamma on each
// interior face.  Row khi-1 therefore removes F_{khi-1/2} from the cell below the top row,
// and row khi must put the same flux back through its lower face for the column to conserve
// momentum.  The top row left it out, at the domain top and at a coarse/fine boundary, and the
// bottom row left out its lower face at a coarse/fine boundary (erf-model/ERF#4822).  The
// missing flux acted as a momentum source in the top cell whenever gamma reached the domain
// top, which YSUNew produces when the bulk Richardson number never crosses Ribcr (kpbl = khi).

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

constexpr int NZ = 8;
constexpr Real DZ  = Real(0.25);   // nominal (computational) spacing
constexpr Real RHO = Real(1.2);
constexpr Real KM  = Real(0.7);    // constant vertical eddy viscosity: rhoAlpha = KM on every face
constexpr Real DT  = Real(0.3);

// A countergradient gradient that is nonzero in every cell, ghosts included, so that every
// face of the column, the top one included, carries countergradient flux
double gamma_cell (int k) { return 0.3 - 0.02 * double(k); }

double u_ghost (int k) { return (k < 0) ? 0.9 : -0.6; }

// Set the column extent of every cell to the whole domain, 0 to NZ-1
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

// Two cell columns in x (the u face between them is solved), one in y, NZ cells in z.  The
// solve runs over the face rows klo..khi of the column; with klo > 0 or khi < NZ-1 those rows
// end at an interior (coarse/fine) face and the solve closes the system with the known face
// momentum outside the range.
struct Column
{
  Box cells{IntVect(0), IntVect(1, 0, NZ-1)};
  Box domain{IntVect(0), IntVect(1, 0, NZ-1)};
  int klo, khi;
  Box col;
  std::vector<Real> znd;                     // z of the NZ+1 nodes, 0 at the bottom
  FArrayBox cons, rho_u, zero_tau, z_nd, detJ, mf, mu_turb;
  IArrayBox col_kext;
  std::vector<BCRec> bcs;
  SolverChoice solver;
  GpuArray<Real, AMREX_SPACEDIM> dxInv{Real(1.0), Real(1.0), Real(1.0)/DZ};

  Column (const std::vector<Real>& z_nodes, int k_lo, int k_hi)
    : klo(k_lo), khi(k_hi), col(IntVect(1, 0, k_lo), IntVect(1, 0, k_hi)), znd(z_nodes)
  {
    const Box g = grow(cells, 2);
    cons.resize(g, 2);
    rho_u.resize(grow(surroundingNodes(cells, 0), 2), 1);
    zero_tau.resize(grow(convert(cells, IntVect(1,0,1)), 2), 1);
    z_nd.resize(grow(convert(cells, IntVect(1)), 2), 1);
    detJ.resize(g, 1);
    mf.resize(grow(surroundingNodes(cells, 0), 2), 1);
    mu_turb.resize(g, EddyDiff::NumDiffs);
    col_kext.resize(g, 2);

    cons.setVal<RunOn::Device>(RHO, cons.box(), Rho_comp, 1);
    cons.setVal<RunOn::Device>(RHO * Real(300.0), cons.box(), RhoTheta_comp, 1);
    zero_tau.setVal<RunOn::Device>(Real(0.0));
    mf.setVal<RunOn::Device>(Real(1.0));
    fill_col_kext(col_kext);
    fill_terrain(znd, z_nd, detJ);

    // Constant eddy viscosity and the countergradient gradient of gamma_cell in every cell
    mu_turb.setVal<RunOn::Device>(Real(0.0));
    mu_turb.setVal<RunOn::Device>(KM, mu_turb.box(), EddyDiff::Mom_v, 1);
    {
      FArrayBox h(mu_turb.box(), 1, The_Pinned_Arena());
      auto ha = h.array();
      for (int k = h.box().smallEnd(2); k <= h.box().bigEnd(2); ++k) {
        for (int j = h.box().smallEnd(1); j <= h.box().bigEnd(1); ++j) {
          for (int i = h.box().smallEnd(0); i <= h.box().bigEnd(0); ++i) {
            ha(i,j,k) = Real(gamma_cell(k));
          }
        }
      }
      Gpu::copy(Gpu::hostToDevice, h.dataPtr(0), h.dataPtr(0) + h.size(),
                mu_turb.dataPtr(EddyDiff::HGAMU_v));
      Gpu::streamSynchronize();
    }

    // Zero-flux (foextrap) z boundaries: the domain faces carry no stress
    bcs.resize(BCVars::NumTypes);
    for (auto& bc : bcs) {
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        bc.setLo(d, (d < 2) ? ERFBCType::int_dir : ERFBCType::foextrap);
        bc.setHi(d, (d < 2) ? ERFBCType::int_dir : ERFBCType::foextrap);
      }
    }
    solver.diffChoice.molec_diff_type = MolecDiffType::None;
    solver.diffChoice.dynamic_viscosity = Real(0.0);
    solver.turbChoice.resize(1);
    solver.turbChoice[0].use_kturb = true;
  }

  double zc (int k) const { return 0.5 * (double(znd[k]) + double(znd[k+1])); }
  double dz (int k) const { return double(znd[k+1]) - double(znd[k]); }

  // rho_u = RHO * u_k in the column and RHO * u_ghost outside the solved rows
  void set_state (const std::vector<double>& u)
  {
    FArrayBox h(rho_u.box(), 1, The_Pinned_Arena());
    auto ha = h.array();
    for (int k = rho_u.box().smallEnd(2); k <= rho_u.box().bigEnd(2); ++k) {
      const double v = (k < klo) ? RHO * u_ghost(-1)
                     : (k > khi) ? RHO * u_ghost(NZ)
                     : RHO * u[k];
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
    std::vector<double> u(NZ, 0.0);
    for (int k = klo; k <= khi; ++k) { u[k] = double(ha(1,0,k)) / double(RHO); }
    return u;
  }

  // Column momentum sum_k h_k rho u_k over the solved rows
  double momentum (const std::vector<double>& u) const
  {
    double m = 0.0;
    for (int k = klo; k <= khi; ++k) { m += dz(k) * double(RHO) * u[k]; }
    return m;
  }

  void solve_N ()
  {
    ImplicitDiffForMomLU_N<0>(col, domain, 0, double(DT), col_kext.const_array(),
                              cons.const_array(), rho_u.array(),
                              zero_tau.const_array(), zero_tau.const_array(), dxInv,
                              mu_turb.const_array(), solver, bcs.data(), false, Real(1.0), true);
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
                              mu_turb.const_array(), solver, bcs.data(), false, Real(1.0), true);
    Gpu::streamSynchronize();
  }

  void solve_T ()
  {
    ImplicitDiffForMomLU_T<0>(col, domain, 0, double(DT), col_kext.const_array(),
                              cons.const_array(), rho_u.array(),
                              zero_tau.const_array(), zero_tau.const_array(),
                              z_nd.const_array(), detJ.const_array(), dxInv,
                              mu_turb.const_array(), solver, bcs.data(), false, Real(1.0), true,
                              mf.const_array(), mf.const_array(), false);
    Gpu::streamSynchronize();
  }

  // The tridiagonal system the solve should build, written out independently:
  //   h_k rho u_k - DT (D_{k+1/2} - D_{k-1/2}) = h_k rho u*_k - DT (F_{k+1/2} - F_{k-1/2})
  // with the diffusive flux D = KM du/dz centred between cell centres and the countergradient
  // flux F = KM gamma, gamma averaged to the face.  A domain face carries neither; an interior
  // (coarse/fine) face at the end of the solved range carries both, with the known ghost value.
  std::vector<double> reference (const std::vector<double>& u) const
  {
    const double km = double(KM), dt = double(DT), rho = double(RHO);
    const int n = khi - klo + 1;
    std::vector<double> a(n, 0.0), b(n), c(n, 0.0), r(n);
    for (int k = klo; k <= khi; ++k) {
      const int m = k - klo;
      b[m] = rho * dz(k);
      r[m] = rho * dz(k) * u[k];
      if (k > 0) {
        const double g = dt * km / (zc(k) - zc(k-1));
        b[m] += g;
        if (k > klo) { a[m] -= g; } else { r[m] += g * u_ghost(-1); }
        r[m] += dt * km * 0.5 * (gamma_cell(k) + gamma_cell(k-1));
      }
      if (k < NZ-1) {
        const double g = dt * km / (zc(k+1) - zc(k));
        b[m] += g;
        if (k < khi) { c[m] -= g; } else { r[m] += g * u_ghost(NZ); }
        r[m] -= dt * km * 0.5 * (gamma_cell(k) + gamma_cell(k+1));
      }
    }
    const auto x = thomas(a, b, c, r);
    std::vector<double> ref(NZ, 0.0);
    for (int k = klo; k <= khi; ++k) { ref[k] = x[k-klo]; }
    return ref;
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

std::vector<double>
initial_velocity ()
{
  std::vector<double> u(NZ);
  for (int k = 0; k < NZ; ++k) { u[k] = 0.2 + 0.1 * double(k) - 0.03 * double(k*k); }
  return u;
}

void
run_solve (Column& c, int form)
{
  if (form == 0) { c.solve_N(); } else if (form == 1) { c.solve_S(); } else { c.solve_T(); }
}

const char* form_name (int form) { return (form == 0) ? "N" : ((form == 1) ? "S" : "T"); }

} // namespace

// Motivation: with zero stress on both domain faces, one implicit step must leave the column
// momentum unchanged, because the countergradient flux only moves momentum between cells.
// Without the lower-face term in the top row the column lost DT KM gamma_{khi-1/2} (0.036
// here) per step on every mesh.
TEST(ImplicitMomCountergradient, DomainColumnConservesMomentum)
{
  const auto u = initial_velocity();
  for (int form = 0; form < 3; ++form) {
    Column c((form == 0) ? uniform_levels() : stretched_levels(), 0, NZ-1);
    c.set_state(u);
    run_solve(c, form);
    const auto un = c.velocity();
    EXPECT_NEAR(c.momentum(un), c.momentum(u), double(tol_for_scale(Real(2.0)))) << form_name(form);
  }
}

// Motivation: on every mesh the solve must build the tridiagonal system with the countergradient
// flux on every interior face and none on the domain faces, written out here independently.
TEST(ImplicitMomCountergradient, DomainColumnMatchesTheReference)
{
  const auto u = initial_velocity();
  for (int form = 0; form < 3; ++form) {
    Column c((form == 0) ? uniform_levels() : stretched_levels(), 0, NZ-1);
    c.set_state(u);
    run_solve(c, form);
    const auto un  = c.velocity();
    const auto ref = c.reference(u);
    for (int k = 0; k < NZ; ++k) {
      EXPECT_NEAR(un[k], ref[k], double(tol_for_scale(Real(2.0)))) << form_name(form) << " k " << k;
    }
  }
}

// Motivation: a column that ends at a coarse/fine boundary at either end closes the system
// with the known momentum outside it, and the countergradient flux passes through both of its
// end faces, as it does on every interior face.  Both the bottom row (lower face) and the top
// row (both faces) left it out.
TEST(ImplicitMomCountergradient, CoarseFineColumnMatchesTheReference)
{
  const auto u = initial_velocity();
  for (int form = 0; form < 3; ++form) {
    Column c((form == 0) ? uniform_levels() : stretched_levels(), 1, NZ-2);
    c.set_state(u);
    run_solve(c, form);
    const auto un  = c.velocity();
    const auto ref = c.reference(u);
    for (int k = c.klo; k <= c.khi; ++k) {
      EXPECT_NEAR(un[k], ref[k], double(tol_for_scale(Real(2.0)))) << form_name(form) << " k " << k;
    }
  }
}

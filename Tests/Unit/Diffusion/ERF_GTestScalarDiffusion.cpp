#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include <ERF_Diffusion.H>
#include <ERF_NativeScalarDiffusion.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace amrex;

// NVCC requires extended device lambdas to be enclosed by a public member.
// GTEST_TEST_ makes TestBody private, so keep the same registration with a
// public TestBody for tests that launch GPU kernels.
#define ERF_GPU_TEST(test_suite_name, test_name)                              \
  static_assert(sizeof(GTEST_STRINGIFY_(test_suite_name)) > 1,                 \
                "test_suite_name must not be empty");                        \
  static_assert(sizeof(GTEST_STRINGIFY_(test_name)) > 1,                      \
                "test_name must not be empty");                              \
  class GTEST_TEST_CLASS_NAME_(test_suite_name, test_name)                    \
      : public ::testing::Test {                                              \
   public:                                                                    \
    GTEST_TEST_CLASS_NAME_(test_suite_name, test_name)() = default;           \
    ~GTEST_TEST_CLASS_NAME_(test_suite_name, test_name)() override = default; \
    GTEST_TEST_CLASS_NAME_(test_suite_name, test_name)                       \
    (const GTEST_TEST_CLASS_NAME_(test_suite_name, test_name) &) = delete;    \
    GTEST_TEST_CLASS_NAME_(test_suite_name, test_name) & operator=(           \
        const GTEST_TEST_CLASS_NAME_(test_suite_name, test_name) &) = delete; \
    GTEST_TEST_CLASS_NAME_(test_suite_name, test_name)                       \
    (GTEST_TEST_CLASS_NAME_(test_suite_name, test_name) &&) noexcept = delete; \
    GTEST_TEST_CLASS_NAME_(test_suite_name, test_name) & operator=(           \
        GTEST_TEST_CLASS_NAME_(test_suite_name, test_name) &&) noexcept =     \
      delete;                                                                 \
    void TestBody() override;                                                 \
    [[maybe_unused]] static ::testing::TestInfo* const test_info_;            \
  };                                                                          \
  ::testing::TestInfo* const GTEST_TEST_CLASS_NAME_(test_suite_name,          \
                                                    test_name)::test_info_ =  \
      ::testing::internal::MakeAndRegisterTestInfo(                           \
          #test_suite_name, #test_name, nullptr, nullptr,                    \
          ::testing::internal::CodeLocation(__FILE__, __LINE__),              \
          ::testing::internal::GetTestTypeId(),                              \
          ::testing::internal::SuiteApiResolver<                             \
              ::testing::Test>::GetSetUpCaseOrSuite(__FILE__, __LINE__),     \
          ::testing::internal::SuiteApiResolver<                             \
              ::testing::Test>::GetTearDownCaseOrSuite(__FILE__, __LINE__),  \
          new ::testing::internal::TestFactoryImpl<                          \
              GTEST_TEST_CLASS_NAME_(test_suite_name, test_name)>);          \
  void GTEST_TEST_CLASS_NAME_(test_suite_name, test_name)::TestBody()

namespace {

constexpr int kScalarComp = 2;
constexpr int kRhoComp = 1;
constexpr int kFluxComp = 3;
constexpr int kRhsComp = 5;

Real
tolerance(Real scale = Real(1.0))
{
  return Real(96.0) * std::numeric_limits<Real>::epsilon() *
         std::max(Real(1.0), std::abs(scale));
}

// The public terrain adapter composes face gradients, affine metric
// corrections, vertical-face interpolation, and mapped divergence. Its
// analytic RHS comparisons accumulate about 290 eps, so keep 384 eps
// of headroom here while the pointwise primitive checks retain 96 eps.
Real
integrated_terrain_tolerance(Real scale = Real(1.0))
{
  return Real(384.0) * std::numeric_limits<Real>::epsilon() *
         std::max(Real(1.0), std::abs(scale));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
chi(int i, int j, int k) noexcept
{
  return Real(0.7) + Real(0.13) * i - Real(0.09) * j + Real(0.11) * k +
         Real(0.017) * i * j - Real(0.012) * i * k + Real(0.023) * j * k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
rho_value(int i, int j, int k) noexcept
{
  return Real(1.1) + Real(0.03) * i + Real(0.02) * j + Real(0.01) * k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
mu_value(int comp, int i, int j, int k) noexcept
{
  const Real base = comp == EddyDiff::Theta_h    ? Real(0.17)
                    : comp == EddyDiff::Scalar_h ? Real(0.31)
                    : comp == EddyDiff::Q_h      ? Real(0.47)
                    : comp == EddyDiff::Theta_v  ? Real(0.23)
                    : comp == EddyDiff::Scalar_v ? Real(0.39)
                    : comp == EddyDiff::Q_v      ? Real(0.53)
                                                 : Real(0.07);
  return base + Real(0.004) * i + Real(0.003) * j + Real(0.002) * k;
}

void
copy_to_host(const FArrayBox& src, FArrayBox& dst)
{
  Gpu::copy(
    Gpu::deviceToHost, src.dataPtr(0), src.dataPtr(0) + src.size(),
    dst.dataPtr(0));
  Gpu::streamSynchronize();
}

Box
grow_box(Box box, const int n)
{
  box.grow(n);
  return box;
}

Box
terrain_node_box(const Box& cells)
{
  Box nodes = surroundingNodes(cells, 2);
  nodes.grow(0, 1);
  nodes.grow(1, 1);
  return nodes;
}

/** Small production-adapter fixture with analytic affine terrain and scalar. */
struct NativeTerrainScalarCase
{
  Real a, b, c, mx, my, K;
  int qty_comp;
  bool mixed_quadratic;
  Box domain, bx, data_box, map_cells, znd_box;
  FArrayBox conserved, primitive, rhs, u, v, xflux, yflux, zflux;
  FArrayBox smn, mu, mf_mx, mf_my, mf_ux, mf_uy, mf_vx, mf_vy;
  FArrayBox hfx_x, hfx_y, hfx_z, qfx1_x, qfx1_y, qfx1_z, qfx2_z;
  FArrayBox diss, tm, ax, ay, detj, z_nd, z_cc;
  Vector<BCRec> bcs;
  Gpu::DeviceVector<BCRec> bcs_device;
  Vector<std::unique_ptr<SurfaceLayer>> surface;
  SolverChoice solver;
  GpuArray<Real, AMREX_SPACEDIM> inv{{Real(1.0), Real(1.0), Real(1.0)}};
  GpuArray<Real, AMREX_SPACEDIM> gravity{
    {Real(0.0), Real(0.0), Real(-9.81)}};

  NativeTerrainScalarCase(
    const Real a_in,
    const Real b_in,
    const Real c_in,
    const Real mx_in,
    const Real my_in,
    const Real K_in,
    const int qty_comp_in,
    const bool mixed_quadratic_in = false)
    : a(a_in),
      b(b_in),
      c(c_in),
      mx(mx_in),
      my(my_in),
      K(K_in),
      qty_comp(qty_comp_in),
      mixed_quadratic(mixed_quadratic_in),
      domain(IntVect(0, 0, 0), IntVect(7, 7, 7)),
      bx(IntVect(2, 2, 2), IntVect(5, 5, 5)),
      data_box(grow_box(domain, 2)),
      map_cells(IntVect(0, 0, 0), IntVect(7, 7, 0)),
      znd_box(terrain_node_box(data_box)),
      conserved(data_box, NVAR_max),
      primitive(data_box, NPRIMVAR_max),
      rhs(domain, NVAR_max),
      u(surroundingNodes(domain, 0), 1),
      v(surroundingNodes(domain, 1), 1),
      xflux(surroundingNodes(domain, 0), 1),
      yflux(surroundingNodes(domain, 1), 1),
      zflux(surroundingNodes(domain, 2), 1),
      smn(data_box, 1),
      mu(data_box, EddyDiff::NumDiffs),
      mf_mx(map_cells, 1),
      mf_my(map_cells, 1),
      mf_ux(surroundingNodes(map_cells, 0), 1),
      mf_uy(surroundingNodes(map_cells, 0), 1),
      mf_vx(surroundingNodes(map_cells, 1), 1),
      mf_vy(surroundingNodes(map_cells, 1), 1),
      hfx_x(surroundingNodes(domain, 0), 1),
      hfx_y(surroundingNodes(domain, 1), 1),
      hfx_z(surroundingNodes(domain, 2), 1),
      qfx1_x(surroundingNodes(domain, 0), 1),
      qfx1_y(surroundingNodes(domain, 1), 1),
      qfx1_z(surroundingNodes(domain, 2), 1),
      qfx2_z(surroundingNodes(domain, 2), 1),
      diss(data_box, 1),
      tm(data_box, 1),
      ax(surroundingNodes(domain, 0), 1),
      ay(surroundingNodes(domain, 1), 1),
      detj(data_box, 1),
      z_nd(znd_box, 1),
      z_cc(data_box, 1)
  {
    conserved.setVal<RunOn::Device>(Real(0.0));
    primitive.setVal<RunOn::Device>(Real(0.0));
    rhs.setVal<RunOn::Device>(Real(0.0));
    u.setVal<RunOn::Device>(Real(0.0));
    v.setVal<RunOn::Device>(Real(0.0));
    xflux.setVal<RunOn::Device>(Real(0.0));
    yflux.setVal<RunOn::Device>(Real(0.0));
    zflux.setVal<RunOn::Device>(Real(0.0));
    smn.setVal<RunOn::Device>(Real(0.0));
    mu.setVal<RunOn::Device>(Real(0.0));
    mf_mx.setVal<RunOn::Device>(mx);
    mf_my.setVal<RunOn::Device>(my);
    mf_ux.setVal<RunOn::Device>(mx);
    mf_vx.setVal<RunOn::Device>(mx);
    mf_uy.setVal<RunOn::Device>(my);
    mf_vy.setVal<RunOn::Device>(my);
    hfx_x.setVal<RunOn::Device>(Real(0.0));
    hfx_y.setVal<RunOn::Device>(Real(0.0));
    hfx_z.setVal<RunOn::Device>(Real(0.0));
    qfx1_x.setVal<RunOn::Device>(Real(0.0));
    qfx1_y.setVal<RunOn::Device>(Real(0.0));
    qfx1_z.setVal<RunOn::Device>(Real(0.0));
    qfx2_z.setVal<RunOn::Device>(Real(0.0));
    diss.setVal<RunOn::Device>(Real(0.0));
    tm.setVal<RunOn::Device>(Real(0.0));
    // Terrain horizontal face areas are dz/dzeta for the affine mesh, as in
    // make_areas; detJ and zeta spacing use the same c scale.
    ax.setVal<RunOn::Device>(c);
    ay.setVal<RunOn::Device>(c);
    detj.setVal<RunOn::Device>(c);

    initialize_fields();

    bcs.resize(NBCVAR_max);
    for (auto& bc : bcs) {
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        bc.setLo(d, ERFBCType::foextrap);
        bc.setHi(d, ERFBCType::foextrap);
      }
    }
    bcs_device.resize(bcs.size());
    Gpu::copy(Gpu::hostToDevice, bcs.begin(), bcs.end(), bcs_device.begin());
    surface.resize(6);

    solver.diffChoice.molec_diff_type = MolecDiffType::Constant;
    solver.diffChoice.rhoAlpha_C = K;
    solver.diffChoice.rhoAlpha_T = K;
    solver.turbChoice.resize(1);
    solver.turbChoice[0].use_kturb = false;
  }

  void initialize_fields()
  {
    const int scalar_comp = qty_comp - 1;
    auto cons = conserved.array();
    auto prim = primitive.array();
    const Real aa = a, bb = b, cc = c, mmx = mx, mmy = my;
    const bool mixed = mixed_quadratic;
    ParallelFor(data_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const Real xi = Real(i) + Real(0.5);
      const Real eta = Real(j) + Real(0.5);
      const Real zeta = Real(k) + Real(0.5);
      const Real x = xi / mmx;
      const Real y = eta / mmy;
      const Real z = aa * xi + bb * eta + cc * zeta;
      cons(i, j, k, Rho_comp) = Real(1.0);
      Real value = x * x + y * y + z * z;
      if (mixed) value += x * z + y * z;
      prim(i, j, k, scalar_comp) = value;
    });
    auto znd = z_nd.array();
    ParallelFor(znd_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      znd(i, j, k) = aa * Real(i) + bb * Real(j) + cc * Real(k);
    });
    auto zcc = z_cc.array();
    ParallelFor(data_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      zcc(i, j, k) = aa * (Real(i) + Real(0.5)) +
                      bb * (Real(j) + Real(0.5)) +
                      cc * (Real(k) + Real(0.5));
    });
    Gpu::streamSynchronize();
  }

  void run(const Real implicit_fac)
  {
    rhs.setVal<RunOn::Device>(Real(0.0));
    xflux.setVal<RunOn::Device>(Real(-91.0));
    yflux.setVal<RunOn::Device>(Real(-92.0));
    zflux.setVal<RunOn::Device>(Real(-93.0));
    hfx_z.setVal<RunOn::Device>(Real(-94.0));
    qfx1_z.setVal<RunOn::Device>(Real(-95.0));
    qfx2_z.setVal<RunOn::Device>(Real(-96.0));
    auto hfx_x_arr = hfx_x.array();
    auto hfx_y_arr = hfx_y.array();
    auto hfx_z_arr = hfx_z.array();
    auto qfx1_x_arr = qfx1_x.array();
    auto qfx1_y_arr = qfx1_y.array();
    auto qfx1_z_arr = qfx1_z.array();
    auto qfx2_z_arr = qfx2_z.array();
    auto diss_arr = diss.array();
    DiffusionSrcForState_T(
      bx, domain, qty_comp, 1, false, u.const_array(), v.const_array(),
      conserved.const_array(), primitive.const_array(), rhs.array(),
      xflux.array(), yflux.array(), zflux.array(), z_nd.const_array(),
      z_cc.const_array(), ax.const_array(), ay.const_array(), ax.const_array(),
      detj.const_array(), inv, smn.const_array(), mf_mx.const_array(),
      mf_ux.const_array(), mf_vx.const_array(), mf_my.const_array(),
      mf_uy.const_array(), mf_vy.const_array(), hfx_x_arr, hfx_y_arr,
      hfx_z_arr, qfx1_x_arr, qfx1_y_arr, qfx1_z_arr, qfx2_z_arr, diss_arr,
      mu.const_array(), solver, 0,
      tm.const_array(), gravity, bcs_device.data(), false, surface,
      implicit_fac);
    Gpu::streamSynchronize();
  }
};

template <bool MultiplyMolecularByDensity, bool AddTurbulence>
void
build_n_and_check(
  const Box& bx,
  const FArrayBox& scalar,
  const FArrayBox& rho,
  const FArrayBox& mu,
  const FArrayBox& mf_ux,
  const FArrayBox& mf_uy,
  const FArrayBox& mf_vy,
  const FArrayBox& mf_vx,
  const FArrayBox& mf_mx,
  const FArrayBox& mf_my,
  FArrayBox& xflux,
  FArrayBox& yflux,
  FArrayBox& zflux,
  FArrayBox& rhs,
  const ScalarDiffusionCoefficients& coefficients)
{
  const Real dx_inv = Real(0.7), dy_inv = Real(0.9), dz_inv = Real(1.2);
  const auto scalar4 = scalar.const_array();
  const auto rho4 = rho.const_array();
  const auto mu4 = mu.const_array();
  const auto mx4 = mf_mx.const_array();
  const auto my4 = mf_my.const_array();
  const auto ux4 = mf_ux.const_array();
  const auto uy4 = mf_uy.const_array();
  const auto vy4 = mf_vy.const_array();
  const auto vx4 = mf_vx.const_array();
  const auto fx4 = xflux.array();
  const auto fy4 = yflux.array();
  const auto fz4 = zflux.array();
  const auto rhs4 = rhs.array();
  const Box xbx = surroundingNodes(bx, 0);
  const Box ybx = surroundingNodes(bx, 1);
  const Box zbx = surroundingNodes(bx, 2);

  ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real kx =
      ScalarDiffusionFaceCoefficient<MultiplyMolecularByDensity, AddTurbulence>(
        rho4, kRhoComp, mu4, coefficients, i, j, k, 1, 0, 0,
        coefficients.eddy_h_comp);
    fx4(i, j, k, kFluxComp) = ScalarDiffusionFlux_N<0>(
      scalar4, kScalarComp, i, j, k, kx, dx_inv, ux4(i, j, 0), uy4(i, j, 0),
      false, false);
  });
  ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real ky =
      ScalarDiffusionFaceCoefficient<MultiplyMolecularByDensity, AddTurbulence>(
        rho4, kRhoComp, mu4, coefficients, i, j, k, 0, 1, 0,
        coefficients.eddy_h_comp);
    fy4(i, j, k, kFluxComp) = ScalarDiffusionFlux_N<1>(
      scalar4, kScalarComp, i, j, k, ky, dy_inv, vy4(i, j, 0), vx4(i, j, 0),
      false, false);
  });
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real kz =
      ScalarDiffusionFaceCoefficient<MultiplyMolecularByDensity, AddTurbulence>(
        rho4, kRhoComp, mu4, coefficients, i, j, k, 0, 0, 1,
        coefficients.eddy_v_comp);
    fz4(i, j, k, kFluxComp) = ScalarDiffusionFlux_N<2>(
      scalar4, kScalarComp, i, j, k, kz, dz_inv, Real(1.0), Real(1.0), false,
      false);
  });
  auto rhs_reset = rhs.array();
  ParallelFor(
    rhs.box(), rhs.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      rhs_reset(i, j, k, n) = Real(-800.0) - n;
    });
  ApplyScalarDiffusionFluxDivergence_N(
    bx, xflux.const_array(), yflux.const_array(), zflux.const_array(),
    kFluxComp, rhs4, kRhsComp,
    GpuArray<Real, AMREX_SPACEDIM>{{dx_inv, dy_inv, dz_inv}}, mx4, my4);
  Gpu::streamSynchronize();

  FArrayBox hscalar(scalar.box(), scalar.nComp(), The_Pinned_Arena());
  FArrayBox hfx(xflux.box(), xflux.nComp(), The_Pinned_Arena());
  FArrayBox hfy(yflux.box(), yflux.nComp(), The_Pinned_Arena());
  FArrayBox hfz(zflux.box(), zflux.nComp(), The_Pinned_Arena());
  FArrayBox hrhs(rhs.box(), rhs.nComp(), The_Pinned_Arena());
  copy_to_host(scalar, hscalar);
  copy_to_host(xflux, hfx);
  copy_to_host(yflux, hfy);
  copy_to_host(zflux, hfz);
  copy_to_host(rhs, hrhs);
  const auto hfx4 = hfx.const_array();
  const auto hfy4 = hfy.const_array();
  const auto hfz4 = hfz.const_array();
  const auto hrhs4 = hrhs.const_array();
  const auto hscalar4 = hscalar.const_array();

  auto expected_k = [&](int i, int j, int k, int di, int dj, int dk, int eddy) {
    Real value = coefficients.molecular_coeff;
    if constexpr (MultiplyMolecularByDensity) {
      value *=
        Real(0.5) * (rho_value(i, j, k) + rho_value(i - di, j - dj, k - dk));
    }
    if constexpr (AddTurbulence) {
      value += Real(0.5) * (mu_value(eddy, i, j, k) +
                            mu_value(eddy, i - di, j - dj, k - dk));
    }
    return value;
  };

  for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
    for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
      for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
        const Real mx = Real(0.91) + Real(0.025) * i + Real(0.01) * j;
        const Real my = Real(1.13) + Real(0.02) * j;
        const auto expected_xflux = [&](int fi) {
          const Real K =
            expected_k(fi, j, k, 1, 0, 0, coefficients.eddy_h_comp);
          const Real ratio =
            (Real(0.88) + Real(0.03) * fi) / (Real(1.07) + Real(0.02) * j);
          return -K * (chi(fi, j, k) - chi(fi - 1, j, k)) * dx_inv * ratio;
        };
        const auto expected_yflux = [&](int fj) {
          const Real K =
            expected_k(i, fj, k, 0, 1, 0, coefficients.eddy_h_comp);
          const Real ratio =
            (Real(1.04) + Real(0.025) * fj) / (Real(0.93) + Real(0.015) * i);
          return -K * (chi(i, fj, k) - chi(i, fj - 1, k)) * dy_inv * ratio;
        };
        const auto expected_zflux = [&](int fk) {
          const Real K =
            expected_k(i, j, fk, 0, 0, 1, coefficients.eddy_v_comp);
          return -K * (chi(i, j, fk) - chi(i, j, fk - 1)) * dz_inv;
        };
        EXPECT_NEAR(
          hfx4(i, j, k, kFluxComp), expected_xflux(i),
          tolerance(expected_xflux(i)));
        EXPECT_NEAR(
          hfy4(i, j, k, kFluxComp), expected_yflux(j),
          tolerance(expected_yflux(j)));
        EXPECT_NEAR(
          hfz4(i, j, k, kFluxComp), expected_zflux(k),
          tolerance(expected_zflux(k)));
        const Real tendency =
          -((expected_xflux(i + 1) - expected_xflux(i)) * dx_inv * mx * my +
            (expected_yflux(j + 1) - expected_yflux(j)) * dy_inv * mx * my +
            (expected_zflux(k + 1) - expected_zflux(k)) * dz_inv);
        const Real expected = Real(-800.0) - kRhsComp + tendency;
        EXPECT_NEAR(hrhs4(i, j, k, kRhsComp), expected, tolerance(expected));
        for (int n = 0; n < scalar.nComp(); ++n) {
          if (n != kScalarComp) {
            EXPECT_DOUBLE_EQ(hscalar4(i, j, k, n), Real(200.0) + n);
          }
        }
        for (int n = 0; n < rhs.nComp(); ++n) {
          if (n != kRhsComp) {
            EXPECT_DOUBLE_EQ(hrhs4(i, j, k, n), Real(-800.0) - n);
          }
        }
      }
    }
  }
  for (int k = xflux.box().smallEnd(2); k <= xflux.box().bigEnd(2); ++k) {
    for (int j = xflux.box().smallEnd(1); j <= xflux.box().bigEnd(1); ++j) {
      for (int i = xflux.box().smallEnd(0); i <= xflux.box().bigEnd(0); ++i) {
        for (int n = 0; n < xflux.nComp(); ++n) {
          if (n != kFluxComp) {
            EXPECT_DOUBLE_EQ(hfx4(i, j, k, n), Real(900.0) + n);
          }
        }
      }
    }
  }
  for (int k = yflux.box().smallEnd(2); k <= yflux.box().bigEnd(2); ++k) {
    for (int j = yflux.box().smallEnd(1); j <= yflux.box().bigEnd(1); ++j) {
      for (int i = yflux.box().smallEnd(0); i <= yflux.box().bigEnd(0); ++i) {
        for (int n = 0; n < yflux.nComp(); ++n) {
          if (n != kFluxComp) {
            EXPECT_DOUBLE_EQ(hfy4(i, j, k, n), Real(900.0) + n);
          }
        }
      }
    }
  }
  for (int k = zflux.box().smallEnd(2); k <= zflux.box().bigEnd(2); ++k) {
    for (int j = zflux.box().smallEnd(1); j <= zflux.box().bigEnd(1); ++j) {
      for (int i = zflux.box().smallEnd(0); i <= zflux.box().bigEnd(0); ++i) {
        for (int n = 0; n < zflux.nComp(); ++n) {
          if (n != kFluxComp) {
            EXPECT_DOUBLE_EQ(hfz4(i, j, k, n), Real(900.0) + n);
          }
        }
      }
    }
  }
}

void
initialize_n_case(
  FArrayBox& scalar,
  FArrayBox& rho,
  FArrayBox& mu,
  FArrayBox& mf_ux,
  FArrayBox& mf_uy,
  FArrayBox& mf_vy,
  FArrayBox& mf_vx,
  FArrayBox& mf_mx,
  FArrayBox& mf_my,
  FArrayBox& xflux,
  FArrayBox& yflux,
  FArrayBox& zflux,
  FArrayBox& rhs)
{
  auto s = scalar.array();
  ParallelFor(
    scalar.box(), scalar.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      s(i, j, k, n) = n == kScalarComp ? chi(i, j, k) : Real(200.0) + n;
    });
  auto r = rho.array();
  ParallelFor(rho.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    r(i, j, k, kRhoComp) = rho_value(i, j, k);
    r(i, j, k, 0) = Real(4.0);
  });
  auto m = mu.array();
  ParallelFor(
    mu.box(), mu.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      const Real base = n == EddyDiff::Theta_h    ? Real(0.17)
                        : n == EddyDiff::Scalar_h ? Real(0.31)
                        : n == EddyDiff::Q_h      ? Real(0.47)
                        : n == EddyDiff::Theta_v  ? Real(0.23)
                        : n == EddyDiff::Scalar_v ? Real(0.39)
                        : n == EddyDiff::Q_v      ? Real(0.53)
                                                  : Real(0.07);
      m(i, j, k, n) =
        base + Real(0.004) * i + Real(0.003) * j + Real(0.002) * k;
    });
  auto ux = mf_ux.array();
  auto uy = mf_uy.array();
  auto vy = mf_vy.array();
  auto vx = mf_vx.array();
  ParallelFor(mf_ux.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    ux(i, j, k) = Real(0.88) + Real(0.03) * i;
    uy(i, j, k) = Real(1.07) + Real(0.02) * j;
  });
  ParallelFor(mf_vy.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    vy(i, j, k) = Real(1.04) + Real(0.025) * j;
    vx(i, j, k) = Real(0.93) + Real(0.015) * i;
  });
  auto mx = mf_mx.array();
  auto my = mf_my.array();
  ParallelFor(mf_mx.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    mx(i, j, k) = Real(0.91) + Real(0.025) * i + Real(0.01) * j;
    my(i, j, k) = Real(1.13) + Real(0.02) * j;
  });
  auto fx = xflux.array();
  auto fy = yflux.array();
  auto fz = zflux.array();
  auto q = rhs.array();
  ParallelFor(
    xflux.box(), xflux.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      fx(i, j, k, n) = Real(900.0) + n;
    });
  ParallelFor(
    yflux.box(), yflux.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      fy(i, j, k, n) = Real(900.0) + n;
    });
  ParallelFor(
    zflux.box(), zflux.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      fz(i, j, k, n) = Real(900.0) + n;
    });
  ParallelFor(
    rhs.box(), rhs.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      q(i, j, k, n) = Real(-800.0) - n;
    });
  Gpu::streamSynchronize();
}

void
fill_ones(FArrayBox& fab)
{
  auto a = fab.array();
  ParallelFor(
    fab.box(), fab.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      a(i, j, k, n) = Real(1.0);
    });
}

} // namespace

// Motivation: Explicit state diffusion used positional coefficient tables
// through Q11. Appended q-state needs the Q eddy category and zero molecular
// coefficient without indexing past them.
ERF_GPU_TEST(
  ScalarDiffusionPolicy,
  NativeMappingPreservesLegacyCategoriesAndExtendedQIsBounded)
{
  DiffChoice choice;
  choice.molec_diff_type = MolecDiffType::ConstantAlpha;
  choice.alpha_T = Real(0.21);
  choice.alpha_C = Real(0.34);

  const auto theta = ResolveNativeScalarDiffusionPolicy(RhoTheta_comp, choice);
  EXPECT_EQ(theta.scalar_comp, RhoTheta_comp - 1);
  EXPECT_EQ(theta.coefficients.eddy_h_comp, EddyDiff::Theta_h);
  EXPECT_EQ(theta.coefficients.eddy_v_comp, EddyDiff::Theta_v);
  EXPECT_DOUBLE_EQ(theta.coefficients.molecular_coeff, choice.alpha_T);
  const auto ke = ResolveNativeScalarDiffusionPolicy(RhoKE_comp, choice);
  EXPECT_DOUBLE_EQ(ke.coefficients.molecular_coeff, Real(0.0));
  EXPECT_EQ(ke.coefficients.eddy_h_comp, EddyDiff::KE_h);
  EXPECT_EQ(ke.coefficients.eddy_v_comp, EddyDiff::KE_v);
  const auto scalar =
    ResolveNativeScalarDiffusionPolicy(RhoScalar_comp, choice);
  EXPECT_DOUBLE_EQ(scalar.coefficients.molecular_coeff, choice.alpha_C);
  EXPECT_EQ(scalar.coefficients.eddy_h_comp, EddyDiff::Scalar_h);
  EXPECT_EQ(scalar.coefficients.eddy_v_comp, EddyDiff::Scalar_v);

  for (int qcomp : {RhoQ1_comp, RhoQ6_comp}) {
    const auto q = ResolveNativeScalarDiffusionPolicy(qcomp, choice);
    EXPECT_DOUBLE_EQ(q.coefficients.molecular_coeff, choice.alpha_C);
    EXPECT_EQ(q.coefficients.eddy_h_comp, EddyDiff::Q_h);
    EXPECT_EQ(q.coefficients.eddy_v_comp, EddyDiff::Q_v);
  }
  for (int qcomp : {RhoQ7_comp, RhoQ11_comp, RhoQ11_comp + 1}) {
    const auto q = ResolveNativeScalarDiffusionPolicy(qcomp, choice);
    EXPECT_DOUBLE_EQ(q.coefficients.molecular_coeff, Real(0.0));
    EXPECT_EQ(q.coefficients.eddy_h_comp, EddyDiff::Q_h);
    EXPECT_EQ(q.coefficients.eddy_v_comp, EddyDiff::Q_v);
    EXPECT_EQ(q.scalar_comp, qcomp - 1);
  }

  choice.molec_diff_type = MolecDiffType::Constant;
  choice.rhoAlpha_T = Real(0.47);
  choice.rhoAlpha_C = Real(0.59);
  EXPECT_DOUBLE_EQ(
    ResolveNativeScalarDiffusionPolicy(RhoTheta_comp, choice)
      .coefficients.molecular_coeff,
    choice.rhoAlpha_T);
  EXPECT_DOUBLE_EQ(
    ResolveNativeScalarDiffusionPolicy(RhoQ6_comp, choice)
      .coefficients.molecular_coeff,
    choice.rhoAlpha_C);
  choice.molec_diff_type = MolecDiffType::None;
  EXPECT_DOUBLE_EQ(
    ResolveNativeScalarDiffusionPolicy(RhoTheta_comp, choice)
      .coefficients.molecular_coeff,
    Real(0.0));
  EXPECT_DOUBLE_EQ(
    ResolveNativeScalarDiffusionPolicy(RhoScalar_comp, choice)
      .coefficients.molecular_coeff,
    Real(0.0));
}

// Motivation: Exercise appended q components through the public N adapter so
// a fixed Q1-Q11 coefficient lookup, stale primitive extent, or wrong Q eddy
// component produces a finite, independently checkable failure.
ERF_GPU_TEST(
  ScalarDiffusionPolicy,
  ExtendedQStateRunsThroughNativeAdapterWithoutFixedTableLookup)
{
  const Box domain(IntVect(0, 0, 0), IntVect(5, 5, 5));
  const Box bx(IntVect(2, 2, 2), IntVect(3, 3, 3));
  Box data_box = domain;
  data_box.grow(2);
  const Box map_cells(IntVect(0, 0, 0), IntVect(5, 5, 0));
  const Box xfaces = surroundingNodes(domain, 0);
  const Box yfaces = surroundingNodes(domain, 1);
  const Box zfaces = surroundingNodes(domain, 2);
  FArrayBox conserved(data_box, NVAR_max + 8);
  FArrayBox primitive(data_box, NPRIMVAR_max + 8);
  FArrayBox rhs(domain, NVAR_max + 8);
  FArrayBox u(xfaces, 1), v(yfaces, 1);
  FArrayBox xflux(xfaces, 1), yflux(yfaces, 1), zflux(zfaces, 1);
  FArrayBox smn(data_box, 1), mu(data_box, EddyDiff::NumDiffs);
  FArrayBox mf_mx(map_cells, 1), mf_my(map_cells, 1);
  FArrayBox mf_ux(surroundingNodes(map_cells, 0), 1),
    mf_uy(surroundingNodes(map_cells, 0), 1),
    mf_vx(surroundingNodes(map_cells, 1), 1),
    mf_vy(surroundingNodes(map_cells, 1), 1);
  FArrayBox hfx_x(xfaces, 1), hfx_y(yfaces, 1), hfx_z(zfaces, 1);
  FArrayBox qfx1_x(xfaces, 1), qfx1_y(yfaces, 1), qfx1_z(zfaces, 1),
    qfx2_z(zfaces, 1), diss(data_box, 1), tm(data_box, 1);

  conserved.setVal<RunOn::Device>(Real(0.0));
  primitive.setVal<RunOn::Device>(Real(-77.0));
  rhs.setVal<RunOn::Device>(Real(-123.0));
  u.setVal<RunOn::Device>(Real(0.0));
  v.setVal<RunOn::Device>(Real(0.0));
  xflux.setVal<RunOn::Device>(Real(0.0));
  yflux.setVal<RunOn::Device>(Real(0.0));
  zflux.setVal<RunOn::Device>(Real(0.0));
  smn.setVal<RunOn::Device>(Real(0.0));
  mu.setVal<RunOn::Device>(Real(0.0));
  mf_mx.setVal<RunOn::Device>(Real(1.0));
  mf_my.setVal<RunOn::Device>(Real(1.0));
  mf_ux.setVal<RunOn::Device>(Real(1.0));
  mf_uy.setVal<RunOn::Device>(Real(1.0));
  mf_vx.setVal<RunOn::Device>(Real(1.0));
  mf_vy.setVal<RunOn::Device>(Real(1.0));
  hfx_x.setVal<RunOn::Device>(Real(0.0));
  hfx_y.setVal<RunOn::Device>(Real(0.0));
  hfx_z.setVal<RunOn::Device>(Real(0.0));
  qfx1_x.setVal<RunOn::Device>(Real(0.0));
  qfx1_y.setVal<RunOn::Device>(Real(0.0));
  qfx1_z.setVal<RunOn::Device>(Real(0.0));
  qfx2_z.setVal<RunOn::Device>(Real(0.0));
  diss.setVal<RunOn::Device>(Real(0.0));
  tm.setVal<RunOn::Device>(Real(0.0));

  constexpr Real Kh = Real(2.1), Kv = Real(4.3);
  auto cons = conserved.array();
  auto turb = mu.array();
  ParallelFor(data_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    cons(i, j, k, Rho_comp) = Real(3.0);
    turb(i, j, k, EddyDiff::Q_h) = Kh;
    turb(i, j, k, EddyDiff::Q_v) = Kv;
    turb(i, j, k, EddyDiff::Theta_h) = Real(0.13);
    turb(i, j, k, EddyDiff::Theta_v) = Real(0.17);
  });
  Gpu::streamSynchronize();

  Vector<BCRec> bcs(NBCVAR_max);
  for (auto& bc : bcs) {
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      bc.setLo(d, ERFBCType::foextrap);
      bc.setHi(d, ERFBCType::foextrap);
    }
  }
  Gpu::DeviceVector<BCRec> bcs_device(bcs.size());
  Gpu::copy(Gpu::hostToDevice, bcs.begin(), bcs.end(), bcs_device.begin());
  Vector<std::unique_ptr<SurfaceLayer>> surface(6);
  SolverChoice solver;
  solver.diffChoice.molec_diff_type = MolecDiffType::None;
  solver.diffChoice.alpha_C = Real(8.0);
  solver.diffChoice.rhoAlpha_C = Real(9.0);
  solver.turbChoice.resize(1);
  solver.turbChoice[0].use_kturb = true;
  const GpuArray<Real, AMREX_SPACEDIM> inv{{Real(1.0), Real(1.0), Real(1.0)}};
  const GpuArray<Real, AMREX_SPACEDIM> gravity{
    {Real(0.0), Real(0.0), Real(-9.81)}};
  const int q_comps[] = {RhoQ11_comp + 1, RhoQ11_comp + 4};

  for (const int qty_comp : q_comps) {
    const int scalar_comp = qty_comp - 1;
    auto prim_target = primitive.array();
    ParallelFor(data_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const Real x = Real(i) + Real(0.5);
      const Real y = Real(j) + Real(0.5);
      const Real z = Real(k) + Real(0.5);
      prim_target(i, j, k, scalar_comp) = x * x + Real(2.0) * y * y +
                                           Real(3.0) * z * z;
    });
    rhs.setVal<RunOn::Device>(Real(-123.0));
    auto hfx_x_arr = hfx_x.array();
    auto hfx_y_arr = hfx_y.array();
    auto hfx_z_arr = hfx_z.array();
    auto qfx1_x_arr = qfx1_x.array();
    auto qfx1_y_arr = qfx1_y.array();
    auto qfx1_z_arr = qfx1_z.array();
    auto qfx2_z_arr = qfx2_z.array();
    auto diss_arr = diss.array();
    DiffusionSrcForState_N(
      bx, domain, qty_comp, 1, u.const_array(), v.const_array(),
      conserved.const_array(), primitive.const_array(), rhs.array(),
      xflux.array(), yflux.array(), zflux.array(), inv, smn.const_array(),
      mf_mx.const_array(), mf_ux.const_array(), mf_vx.const_array(),
      mf_my.const_array(), mf_uy.const_array(), mf_vy.const_array(),
      hfx_x_arr, hfx_y_arr, hfx_z_arr, qfx1_x_arr, qfx1_y_arr, qfx1_z_arr,
      qfx2_z_arr, diss_arr,
      mu.const_array(), solver, 0, tm.const_array(), gravity,
      bcs_device.data(), false, surface, Real(0.0));
    Gpu::streamSynchronize();

    const auto policy = ResolveNativeScalarDiffusionPolicy(
      qty_comp, solver.diffChoice);
    EXPECT_EQ(policy.scalar_comp, scalar_comp);
    EXPECT_EQ(policy.coefficients.eddy_h_comp, EddyDiff::Q_h);
    EXPECT_EQ(policy.coefficients.eddy_v_comp, EddyDiff::Q_v);
    EXPECT_DOUBLE_EQ(policy.coefficients.molecular_coeff, Real(0.0));

    FArrayBox rhs_host(domain, NVAR_max + 8, The_Pinned_Arena());
    FArrayBox xflux_host(xfaces, 1, The_Pinned_Arena());
    FArrayBox yflux_host(yfaces, 1, The_Pinned_Arena());
    FArrayBox zflux_host(zfaces, 1, The_Pinned_Arena());
    copy_to_host(rhs, rhs_host);
    copy_to_host(xflux, xflux_host);
    copy_to_host(yflux, yflux_host);
    copy_to_host(zflux, zflux_host);
    const auto hrhs = rhs_host.const_array();
    const auto hfx = xflux_host.const_array();
    const auto hfy = yflux_host.const_array();
    const auto hfz = zflux_host.const_array();
    const int i = 2, j = 2, k = 2;
    const Real actual = hrhs(i, j, k, qty_comp);
    EXPECT_TRUE(std::isfinite(actual));
    EXPECT_NEAR(
      actual, Real(-123.0) + Real(6.0) * (Kh + Kv), tolerance(Kh + Kv));
    EXPECT_NEAR(hfx(i, j, k), -Real(4.0) * Kh, tolerance(Kh));
    EXPECT_NEAR(hfy(i, j, k), -Real(8.0) * Kh, tolerance(Kh));
    EXPECT_NEAR(hfz(i, j, k), -Real(12.0) * Kv, tolerance(Kv));
    EXPECT_DOUBLE_EQ(hrhs(i, j, k, qty_comp - 1), Real(-123.0));
    EXPECT_DOUBLE_EQ(hrhs(i, j, k, qty_comp + 1), Real(-123.0));
  }
}

// Motivation: unrelated scalar, density, flux, and RHS components expose
// accidental native cons-1 inference, component-zero writes, and wrong
// horizontal/vertical eddy selection.
ERF_GPU_TEST(ScalarDiffusionPrimitives, NExplicitComponentsAndCoefficientModes)
{
  const Box bx(IntVect(1, 1, 1), IntVect(2, 2, 2));
  Box data_box = bx;
  data_box.grow(1);
  Box map_cell(IntVect(0, 0, 0), IntVect(4, 4, 0));
  Box xmap = surroundingNodes(map_cell, 0),
      ymap = surroundingNodes(map_cell, 1);
  FArrayBox scalar(data_box, 4), rho(data_box, 2),
    mu(data_box, EddyDiff::NumDiffs);
  FArrayBox mf_ux(xmap, 1), mf_uy(xmap, 1), mf_vy(ymap, 1), mf_vx(ymap, 1);
  FArrayBox mf_mx(map_cell, 1), mf_my(map_cell, 1);
  FArrayBox xflux(surroundingNodes(bx, 0), 4),
    yflux(surroundingNodes(bx, 1), 4);
  FArrayBox zflux(surroundingNodes(bx, 2), 4), rhs(bx, 6);
  initialize_n_case(
    scalar, rho, mu, mf_ux, mf_uy, mf_vy, mf_vx, mf_mx, mf_my, xflux, yflux,
    zflux, rhs);

  ScalarDiffusionCoefficients coeff{
    Real(0.28), EddyDiff::Scalar_h, EddyDiff::Scalar_v};
  build_n_and_check<true, false>(
    bx, scalar, rho, mu, mf_ux, mf_uy, mf_vy, mf_vx, mf_mx, mf_my, xflux, yflux,
    zflux, rhs, coeff);
  coeff.molecular_coeff = Real(0.42);
  build_n_and_check<false, false>(
    bx, scalar, rho, mu, mf_ux, mf_uy, mf_vy, mf_vx, mf_mx, mf_my, xflux, yflux,
    zflux, rhs, coeff);
  coeff.molecular_coeff = Real(0.0);
  build_n_and_check<false, true>(
    bx, scalar, rho, mu, mf_ux, mf_uy, mf_vy, mf_vx, mf_mx, mf_my, xflux, yflux,
    zflux, rhs, coeff);
  coeff.molecular_coeff = Real(0.19);
  build_n_and_check<true, true>(
    bx, scalar, rho, mu, mf_ux, mf_uy, mf_vy, mf_vx, mf_mx, mf_my, xflux, yflux,
    zflux, rhs, coeff);
}

ERF_GPU_TEST(ScalarDiffusionPrimitives, NMapFactorPreservesLegacyAssociation)
{
  const Box scalar_box(IntVect(0, 0, 0), IntVect(2, 2, 0));
  const Box result_box(IntVect(0, 0, 0), IntVect(0, 0, 0));
  FArrayBox scalar(scalar_box, 3), result(result_box, 6);
  const auto scalar4 = scalar.array();
  const auto result4 = result.array();
  ParallelFor(scalar_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    scalar4(i, j, k, kScalarComp) =
      ((i == 1 && j == 0) || (i == 0 && j == 1)) ? Real(0.91) : Real(0.0);
  });
  constexpr Real coeff = Real(0.31);
  constexpr Real inv_spacing = Real(0.87);
  constexpr Real map_num = Real(1.13);
  constexpr Real map_den = Real(0.79);
  ParallelFor(result_box, [=] AMREX_GPU_DEVICE(int, int, int) noexcept {
    const Real actual_x = ScalarDiffusionFlux_N<0>(
      scalar4, kScalarComp, 1, 0, 0, coeff, inv_spacing, map_num, map_den,
      false, false);
    const Real grad_x =
      scalar4(1, 0, 0, kScalarComp) - scalar4(0, 0, 0, kScalarComp);
    const Real expected_x =
      -coeff * grad_x * inv_spacing * map_num / map_den;
    const Real reassociated_x =
      (-coeff * grad_x * inv_spacing) * (map_num / map_den);
    const Real actual_y = ScalarDiffusionFlux_N<1>(
      scalar4, kScalarComp, 0, 1, 0, coeff, inv_spacing, map_num, map_den,
      false, false);
    const Real grad_y =
      scalar4(0, 1, 0, kScalarComp) - scalar4(0, 0, 0, kScalarComp);
    const Real expected_y =
      -coeff * grad_y * inv_spacing * map_num / map_den;
    const Real reassociated_y =
      (-coeff * grad_y * inv_spacing) * (map_num / map_den);
    result4(0, 0, 0, 0) = actual_x;
    result4(0, 0, 0, 1) = expected_x;
    result4(0, 0, 0, 2) = reassociated_x;
    result4(0, 0, 0, 3) = actual_y;
    result4(0, 0, 0, 4) = expected_y;
    result4(0, 0, 0, 5) = reassociated_y;
  });
  Gpu::streamSynchronize();
  FArrayBox host(result_box, 6, The_Pinned_Arena());
  copy_to_host(result, host);
  const auto values = host.const_array();
  EXPECT_EQ(values(0, 0, 0, 0), values(0, 0, 0, 1));
  EXPECT_NE(values(0, 0, 0, 0), values(0, 0, 0, 2));
  EXPECT_EQ(values(0, 0, 0, 3), values(0, 0, 0, 4));
  EXPECT_NE(values(0, 0, 0, 3), values(0, 0, 0, 5));
}

ERF_GPU_TEST(
  ScalarDiffusionPrimitives, TerrainDivergencePreservesLegacyAssociation)
{
  const Box result_box(IntVect(0, 0, 0), IntVect(127, 0, 0));
  FArrayBox result(result_box, 3);
  const auto result4 = result.array();
  ParallelFor(result_box, [=] AMREX_GPU_DEVICE(int i, int, int) noexcept {
    const Real id = Real(i);
    const Real fx_hi = Real(0.811) + id * Real(0.0173);
    const Real fx_lo = Real(1.137) + id * Real(0.0117);
    const Real fy_hi = Real(0.923) - id * Real(0.0021);
    const Real fy_lo = Real(0.719) + id * Real(0.0037);
    const Real Gz_hi = Real(0.377) + id * Real(0.0053);
    const Real Gz_lo = -Real(0.213) + id * Real(0.0019);
    const Real dx_inv = Real(0.713) + id * Real(0.0013);
    const Real dy_inv = Real(1.037) + id * Real(0.0007);
    const Real dz_inv = Real(0.877) + id * Real(0.0011);
    const Real mx = Real(1.277) + id * Real(0.0027);
    const Real my = Real(0.839) + id * Real(0.0033);
    const Real detJ = Real(1.193) + id * Real(0.0029);
    const Real mfsq = mx * my;
    Real expected =
      (fx_hi - fx_lo) * dx_inv * mfsq +
      (fy_hi - fy_lo) * dy_inv * mfsq +
      (Gz_hi - Gz_lo) * dz_inv;
    expected /= detJ;
    const Real actual = TerrainDiffusionDivergence_T(
      fx_hi, fx_lo, fy_hi, fy_lo, Gz_hi, Gz_lo, dx_inv, dy_inv, dz_inv,
      mfsq, detJ);
    const Real reassociated =
      (mfsq / detJ) *
      ((fx_hi - fx_lo) * dx_inv + (fy_hi - fy_lo) * dy_inv +
       (Gz_hi / mfsq - Gz_lo / mfsq) * dz_inv);
    result4(i, 0, 0, 0) = actual;
    result4(i, 0, 0, 1) = expected;
    result4(i, 0, 0, 2) = reassociated;
  });
  Gpu::streamSynchronize();
  FArrayBox host(result_box, 3, The_Pinned_Arena());
  copy_to_host(result, host);
  const auto values = host.const_array();
  bool found_legacy_only_match = false;
  for (int i = result_box.smallEnd(0); i <= result_box.bigEnd(0); ++i) {
    found_legacy_only_match |=
      values(i, 0, 0, 0) == values(i, 0, 0, 1) &&
      values(i, 0, 0, 0) != values(i, 0, 0, 2);
  }
  EXPECT_TRUE(found_legacy_only_match);
}

// Motivation: Native ERF must delegate N, S, and T spatial work to the same
// component-explicit kernels while preserving its native state and geometry
// mapping.
ERF_GPU_TEST(ScalarDiffusionPrimitives, NativeAdaptersMatchExplicitPrimitives)
{
  const Box cells(IntVect(1, 1, 1), IntVect(2, 2, 2));
  Box data_box = cells;
  data_box.grow(1);
  const Box xfaces = surroundingNodes(cells, 0);
  const Box yfaces = surroundingNodes(cells, 1);
  const Box zfaces = surroundingNodes(cells, 2);
  const Box map_cells(IntVect(0, 0, 0), IntVect(4, 4, 0));
  FArrayBox conserved(data_box, NVAR_max), primitive(data_box, NPRIMVAR_max);
  FArrayBox rhs(cells, NVAR_max), direct_rhs(cells, 1);
  FArrayBox u(xfaces, 1), v(yfaces, 1);
  FArrayBox xflux(xfaces, 1), yflux(yfaces, 1), zflux(zfaces, 1);
  FArrayBox direct_x(xfaces, 1), direct_y(yfaces, 1), direct_z(zfaces, 1);
  FArrayBox smn(data_box, 1), mu(data_box, EddyDiff::NumDiffs);
  FArrayBox mf_mx(map_cells, 1), mf_my(map_cells, 1);
  FArrayBox mf_ux(surroundingNodes(map_cells, 0), 1),
    mf_uy(surroundingNodes(map_cells, 0), 1),
    mf_vx(surroundingNodes(map_cells, 1), 1),
    mf_vy(surroundingNodes(map_cells, 1), 1);
  FArrayBox hfx_x(xfaces, 1), hfx_y(yfaces, 1), hfx_z(zfaces, 1);
  FArrayBox qfx1_x(xfaces, 1), qfx1_y(yfaces, 1), qfx1_z(zfaces, 1),
    qfx2_z(zfaces, 1), diss(data_box, 1), tm(data_box, 1);
  FArrayBox ax(xfaces, 1), ay(yfaces, 1), detj(data_box, 1);
  Box znd_box = surroundingNodes(data_box, 2);
  znd_box.grow(0, 1);
  znd_box.grow(1, 1);
  FArrayBox z_nd(znd_box, 1), z_cc(data_box, 1);

  auto cons = conserved.array();
  ParallelFor(
    data_box, conserved.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      cons(i, j, k, n) = n == Rho_comp ? rho_value(i, j, k) : Real(0.0);
    });
  auto prim = primitive.array();
  ParallelFor(
    data_box, primitive.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      prim(i, j, k, n) = n == kScalarComp ? chi(i, j, k) : Real(100.0) + n;
    });
  auto rhs_a = rhs.array();
  ParallelFor(
    cells, rhs.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      rhs_a(i, j, k, n) = Real(-200.0) - n;
    });
  direct_rhs.setVal<RunOn::Device>(Real(-200.0) - RhoScalar_comp);
  for (FArrayBox* field :
       {&u, &v, &smn, &mu, &hfx_x, &hfx_y, &hfx_z, &qfx1_x, &qfx1_y, &qfx1_z,
        &qfx2_z, &diss, &tm}) {
    field->setVal<RunOn::Device>(Real(0.0));
  }
  fill_ones(mf_mx);
  fill_ones(mf_my);
  fill_ones(mf_ux);
  fill_ones(mf_uy);
  fill_ones(mf_vx);
  fill_ones(mf_vy);
  fill_ones(ax);
  fill_ones(ay);
  fill_ones(detj);
  auto znd = z_nd.array();
  ParallelFor(znd_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    znd(i, j, k) = Real(k);
  });
  auto zcc = z_cc.array();
  ParallelFor(data_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    zcc(i, j, k) = Real(k) + Real(0.5);
  });
  Gpu::streamSynchronize();

  SolverChoice solver;
  solver.diffChoice.molec_diff_type = MolecDiffType::Constant;
  solver.diffChoice.rhoAlpha_C = Real(0.38);
  solver.turbChoice.resize(1);
  solver.turbChoice[0].use_kturb = false;
  auto native_policy =
    ResolveNativeScalarDiffusionPolicy(RhoScalar_comp, solver.diffChoice);
  const auto cell = conserved.const_array();
  const auto prim4 = primitive.const_array();
  const auto mu4 = mu.const_array();
  const auto dfx = direct_x.array();
  const auto dfy = direct_y.array();
  const auto dfz = direct_z.array();
  const GpuArray<Real, AMREX_SPACEDIM> inv{{Real(1.0), Real(1.0), Real(1.0)}};
  ParallelFor(xfaces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real K = ScalarDiffusionFaceCoefficient<false, false>(
      cell, Rho_comp, mu4, native_policy.coefficients, i, j, k, 1, 0, 0,
      native_policy.coefficients.eddy_h_comp);
    dfx(i, j, k) = ScalarDiffusionFlux_N<0>(
      prim4, kScalarComp, i, j, k, K, inv[0], Real(1.0), Real(1.0), false,
      false);
  });
  ParallelFor(yfaces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real K = ScalarDiffusionFaceCoefficient<false, false>(
      cell, Rho_comp, mu4, native_policy.coefficients, i, j, k, 0, 1, 0,
      native_policy.coefficients.eddy_h_comp);
    dfy(i, j, k) = ScalarDiffusionFlux_N<1>(
      prim4, kScalarComp, i, j, k, K, inv[1], Real(1.0), Real(1.0), false,
      false);
  });
  ParallelFor(zfaces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real K = ScalarDiffusionFaceCoefficient<false, false>(
      cell, Rho_comp, mu4, native_policy.coefficients, i, j, k, 0, 0, 1,
      native_policy.coefficients.eddy_v_comp);
    dfz(i, j, k) = ScalarDiffusionFlux_N<2>(
      prim4, kScalarComp, i, j, k, K, inv[2], Real(1.0), Real(1.0), false,
      false);
  });
  ApplyScalarDiffusionFluxDivergence_N(
    cells, direct_x.const_array(), direct_y.const_array(),
    direct_z.const_array(), 0, direct_rhs.array(), 0, inv, mf_mx.const_array(),
    mf_my.const_array());
  Gpu::streamSynchronize();

  Vector<BCRec> bcs(NBCVAR_max);
  for (auto& bc : bcs) {
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      bc.setLo(d, ERFBCType::foextrap);
      bc.setHi(d, ERFBCType::foextrap);
    }
  }
  Gpu::DeviceVector<BCRec> bcs_device(bcs.size());
  Gpu::copy(Gpu::hostToDevice, bcs.begin(), bcs.end(), bcs_device.begin());
  Vector<std::unique_ptr<SurfaceLayer>> surface(6);
  auto hfx_x_arr = hfx_x.array(), hfx_y_arr = hfx_y.array(),
       hfx_z_arr = hfx_z.array();
  auto qfx1_x_arr = qfx1_x.array(), qfx1_y_arr = qfx1_y.array(),
       qfx1_z_arr = qfx1_z.array(), qfx2_z_arr = qfx2_z.array();
  auto diss_arr = diss.array();
  const GpuArray<Real, AMREX_SPACEDIM> gravity{
    {Real(0.0), Real(0.0), Real(-9.81)}};
  const bool rotate = false, use_surface_layer = false;

  auto compare = [&]() {
    Gpu::streamSynchronize();
    FArrayBox hx(xfaces, 1, The_Pinned_Arena());
    FArrayBox hy(yfaces, 1, The_Pinned_Arena());
    FArrayBox hz(zfaces, 1, The_Pinned_Arena());
    FArrayBox hr(cells, NVAR_max, The_Pinned_Arena());
    FArrayBox hdx(xfaces, 1, The_Pinned_Arena());
    FArrayBox hdy(yfaces, 1, The_Pinned_Arena());
    FArrayBox hdz(zfaces, 1, The_Pinned_Arena());
    FArrayBox hdr(cells, 1, The_Pinned_Arena());
    copy_to_host(xflux, hx);
    copy_to_host(yflux, hy);
    copy_to_host(zflux, hz);
    copy_to_host(rhs, hr);
    copy_to_host(direct_x, hdx);
    copy_to_host(direct_y, hdy);
    copy_to_host(direct_z, hdz);
    copy_to_host(direct_rhs, hdr);
    const auto hx4 = hx.const_array(), hy4 = hy.const_array(),
               hz4 = hz.const_array(), hr4 = hr.const_array();
    const auto hdx4 = hdx.const_array(), hdy4 = hdy.const_array(),
               hdz4 = hdz.const_array(), hdr4 = hdr.const_array();
    for (int k = 1; k <= 2; ++k) {
      for (int j = 1; j <= 2; ++j) {
        for (int i = 1; i <= 3; ++i)
          EXPECT_NEAR(hx4(i, j, k), hdx4(i, j, k), tolerance());
      }
    }
    for (int k = 1; k <= 2; ++k) {
      for (int j = 1; j <= 3; ++j) {
        for (int i = 1; i <= 2; ++i)
          EXPECT_NEAR(hy4(i, j, k), hdy4(i, j, k), tolerance());
      }
    }
    for (int k = 1; k <= 3; ++k) {
      for (int j = 1; j <= 2; ++j) {
        for (int i = 1; i <= 2; ++i)
          EXPECT_NEAR(hz4(i, j, k), hdz4(i, j, k), tolerance());
      }
    }
    for (int k = 1; k <= 2; ++k) {
      for (int j = 1; j <= 2; ++j) {
        for (int i = 1; i <= 2; ++i)
          EXPECT_NEAR(hr4(i, j, k, RhoScalar_comp), hdr4(i, j, k), tolerance());
      }
    }
  };

  auto reset_rhs = [&]() {
    rhs.setVal<RunOn::Device>(Real(-200.0) - RhoScalar_comp);
    xflux.setVal<RunOn::Device>(Real(-900.0));
    yflux.setVal<RunOn::Device>(Real(-900.0));
    zflux.setVal<RunOn::Device>(Real(-900.0));
  };
  DiffusionSrcForState_N(
    cells, cells, RhoScalar_comp, 1, u.const_array(), v.const_array(), cell,
    prim4, rhs.array(), xflux.array(), yflux.array(), zflux.array(), inv,
    smn.const_array(), mf_mx.const_array(), mf_ux.const_array(),
    mf_vx.const_array(), mf_my.const_array(), mf_uy.const_array(),
    mf_vy.const_array(), hfx_x_arr, hfx_y_arr, hfx_z_arr, qfx1_x_arr,
    qfx1_y_arr, qfx1_z_arr, qfx2_z_arr, diss_arr, mu4, solver, 0,
    tm.const_array(), gravity, bcs_device.data(), use_surface_layer, surface,
    Real(0.0));
  compare();

  reset_rhs();
  Vector<Real> dz_host(6, Real(1.0));
  Gpu::DeviceVector<Real> dz(6);
  Gpu::copy(Gpu::hostToDevice, dz_host.begin(), dz_host.end(), dz.begin());
  DiffusionSrcForState_S(
    cells, cells, RhoScalar_comp, 1, u.const_array(), v.const_array(), cell,
    prim4, rhs.array(), xflux.array(), yflux.array(), zflux.array(), dz, inv,
    smn.const_array(), mf_mx.const_array(), mf_ux.const_array(),
    mf_vx.const_array(), mf_my.const_array(), mf_uy.const_array(),
    mf_vy.const_array(), hfx_x_arr, hfx_y_arr, hfx_z_arr, qfx1_x_arr,
    qfx1_y_arr, qfx1_z_arr, qfx2_z_arr, diss_arr, mu4, solver, 0,
    tm.const_array(), gravity, bcs_device.data(), use_surface_layer, surface,
    Real(0.0));
  compare();

  reset_rhs();
  DiffusionSrcForState_T(
    cells, cells, RhoScalar_comp, 1, rotate, u.const_array(), v.const_array(),
    cell, prim4, rhs.array(), xflux.array(), yflux.array(), zflux.array(),
    z_nd.const_array(), z_cc.const_array(), ax.const_array(), ay.const_array(),
    ax.const_array(), detj.const_array(), inv, smn.const_array(),
    mf_mx.const_array(), mf_ux.const_array(), mf_vx.const_array(),
    mf_my.const_array(), mf_uy.const_array(), mf_vy.const_array(), hfx_x_arr,
    hfx_y_arr, hfx_z_arr, qfx1_x_arr, qfx1_y_arr, qfx1_z_arr, qfx2_z_arr,
    diss_arr, mu4, solver, 0, tm.const_array(), gravity, bcs_device.data(),
    use_surface_layer, surface, Real(0.0));
  compare();
}

// Motivation: stretched vertical gradients use adjacent cell-center distance
// while divergence uses the receiving cell width; uniform dz must recover N's
// boundary stencil.
ERF_GPU_TEST(ScalarDiffusionPrimitives, StretchedUniformLimitMatchesN)
{
  const auto w = StretchedDiffusionDirichletWeights(Real(1.0), Real(1.0));
  EXPECT_NEAR(w.c1, Real(-8.0 / 3.0), tolerance());
  EXPECT_NEAR(w.c2, Real(3.0), tolerance());
  EXPECT_NEAR(w.c3, Real(-1.0 / 3.0), tolerance());

  const Box cells(IntVect(0, 0, 0), IntVect(0, 0, 2));
  Box scalar_box = cells;
  scalar_box.grow(1);
  FArrayBox scalar(scalar_box, 3), fx(surroundingNodes(cells, 0), 1),
    fy(surroundingNodes(cells, 1), 1), fn(surroundingNodes(cells, 2), 1),
    fs(surroundingNodes(cells, 2), 1), rhsn(cells, 1), rhss(cells, 1),
    mf(cells, 1);
  const Real K = Real(0.73);
  auto s = scalar.array();
  ParallelFor(scalar_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real z = k < 0 ? Real(0.0) : k > 2 ? Real(3.0) : Real(k) + Real(0.5);
    s(i, j, k, kScalarComp) = Real(0.4) - Real(0.2) * z + Real(0.08) * z * z;
  });
  auto m = mf.array();
  ParallelFor(cells, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    m(i, j, k) = Real(1.0);
  });
  auto ax = fx.array();
  auto ay = fy.array();
  ParallelFor(fx.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    ax(i, j, k) = Real(0.0);
  });
  ParallelFor(fy.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    ay(i, j, k) = Real(0.0);
  });
  Vector<Real> dz_host(3, Real(1.0));
  Gpu::DeviceVector<Real> dz(3);
  Gpu::copy(Gpu::hostToDevice, dz_host.begin(), dz_host.end(), dz.begin());
  const auto sp = scalar.const_array();
  const auto fn4 = fn.array();
  const auto fs4 = fs.array();
  const auto dzp = dz.data();
  const Box faces = surroundingNodes(cells, 2);
  ParallelFor(faces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const bool lo = k == 0, hi = k == 3;
    fn4(i, j, k) = ScalarDiffusionFlux_N<2>(
      sp, kScalarComp, i, j, k, K, Real(1.0), Real(1.0), Real(1.0), lo, hi);
    fs4(i, j, k) =
      -K * StretchedScalarGradient(sp, kScalarComp, i, j, k, dzp, 0, 2, lo, hi);
  });
  auto rhsn4 = rhsn.array();
  auto rhss4 = rhss.array();
  ParallelFor(cells, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    rhsn4(i, j, k) = Real(0.0);
    rhss4(i, j, k) = Real(0.0);
  });
  ApplyScalarDiffusionFluxDivergence_N(
    cells, fx.const_array(), fy.const_array(), fn.const_array(), 0,
    rhsn.array(), 0,
    GpuArray<Real, AMREX_SPACEDIM>{{Real(1.0), Real(1.0), Real(1.0)}},
    mf.const_array(), mf.const_array());
  ApplyScalarDiffusionFluxDivergence_S(
    cells, fx.const_array(), fy.const_array(), fs.const_array(), 0,
    rhss.array(), 0, mf.const_array(), mf.const_array(), Real(1.0), Real(1.0),
    dzp);
  Gpu::streamSynchronize();
  FArrayBox hfn(fn.box(), 1, The_Pinned_Arena());
  FArrayBox hfs(fs.box(), 1, The_Pinned_Arena());
  FArrayBox hrn(rhsn.box(), 1, The_Pinned_Arena());
  FArrayBox hrs(rhss.box(), 1, The_Pinned_Arena());
  copy_to_host(fn, hfn);
  copy_to_host(fs, hfs);
  copy_to_host(rhsn, hrn);
  copy_to_host(rhss, hrs);
  const auto hfn4 = hfn.const_array();
  const auto hfs4 = hfs.const_array();
  const auto hrn4 = hrn.const_array();
  const auto hrs4 = hrs.const_array();
  for (int k = 0; k <= 3; ++k) {
    EXPECT_NEAR(hfn4(0, 0, k), hfs4(0, 0, k), tolerance(K));
  }
  for (int k = 0; k <= 2; ++k) {
    EXPECT_NEAR(hrn4(0, 0, k), hrs4(0, 0, k), tolerance(K));
  }
}

// Motivation: On a stretched mesh the scalar gradient uses adjacent cell-center
// spacing while divergence uses the receiving cell width. Unequal dz values
// make confusing these two metrics produce a resolvable error.
ERF_GPU_TEST(ScalarDiffusionPrimitives, StretchedVariableDzUsesCellAndFaceSpacing)
{
  const Box box(IntVect(0, 0, -1), IntVect(0, 0, 4));
  FArrayBox scalar(box, 4);
  constexpr Real a = Real(1.7), b = Real(-0.42), c = Real(0.11);
  Vector<Real> dz_host{Real(0.5), Real(1.2), Real(0.8), Real(1.5)};
  Gpu::DeviceVector<Real> dz(4);
  Gpu::copy(Gpu::hostToDevice, dz_host.begin(), dz_host.end(), dz.begin());
  const Real z0 = Real(0.0);
  const Real z1 = dz_host[0];
  const Real z2 = dz_host[0] + dz_host[1];
  const Real z3 = dz_host[0] + dz_host[1] + dz_host[2];
  auto s = scalar.array();
  const auto dptr = dz.data();
  ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    Real z = 0;
    if (k == -1)
      z = z0;
    else if (k == 0)
      z = Real(0.5) * dptr[0];
    else if (k == 1)
      z = z1 + Real(0.5) * dptr[1];
    else if (k == 2)
      z = z2 + Real(0.5) * dptr[2];
    else if (k == 3)
      z = z3;
    else
      z = z3 + Real(0.5) * dptr[3];
    s(i, j, k, kScalarComp) = a + b * z + c * z * z;
  });
  Gpu::streamSynchronize();
  const auto s4 = scalar.const_array();
  FArrayBox gradients(Box(IntVect(0, 0, 0), IntVect(0, 0, 3)), 1);
  auto g = gradients.array();
  ParallelFor(
    gradients.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const bool low = k == 0, high = k == 3;
      g(i, j, k, 0) = StretchedScalarGradient(
        s4, kScalarComp, i, j, k, dptr, 0, 2, low, high);
    });
  Gpu::streamSynchronize();
  constexpr Real K = Real(0.72);
  FArrayBox flux(gradients.box(), 1),
    rhs(Box(IntVect(0, 0, 0), IntVect(0, 0, 2)), 1);
  const auto gradient = gradients.const_array();
  const auto flux4 = flux.array();
  ParallelFor(
    gradients.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      flux4(i, j, k) = -K * gradient(i, j, k, 0);
    });
  const auto flux_in = flux.const_array();
  const auto rhs4 = rhs.array();
  ParallelFor(rhs.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    rhs4(i, j, k) = -(flux_in(i, j, k + 1) - flux_in(i, j, k)) / dptr[k];
  });
  Gpu::streamSynchronize();
  FArrayBox hg(gradients.box(), gradients.nComp(), The_Pinned_Arena());
  FArrayBox hr(rhs.box(), 1, The_Pinned_Arena());
  copy_to_host(gradients, hg);
  copy_to_host(rhs, hr);
  const auto hg4 = hg.const_array();
  const auto hr4 = hr.const_array();
  const Real zc0 = Real(0.5) * dz_host[0];
  const Real zc1 = z1 + Real(0.5) * dz_host[1];
  const Real zc2 = z2 + Real(0.5) * dz_host[2];
  const Real face_z[] = {
    z0, Real(0.5) * (zc0 + zc1), Real(0.5) * (zc1 + zc2), z3};
  for (int k = 0; k < 4; ++k) {
    const Real expected_gradient = b + Real(2.0) * c * face_z[k];
    EXPECT_NEAR(
      hg4(0, 0, k, 0), expected_gradient, tolerance(expected_gradient));
  }
  for (int k = 0; k < 3; ++k) {
    const Real expected =
      K *
      (b + Real(2.0) * c * face_z[k + 1] - (b + Real(2.0) * c * face_z[k])) /
      dz_host[k];
    EXPECT_NEAR(hr4(0, 0, k), expected, tolerance(expected));
  }
}

// Motivation: Terrain-following diffusion needs both h_xi and h_eta cross
// terms in the transformed vertical transfer. This independent affine-terrain
// quadratic oracle detects either missing term and map-factor misplacement.
ERF_GPU_TEST(ScalarDiffusionPrimitives, TerrainMappedQuadraticManufacturedSolution)
{
  const Real a = Real(0.31), b = Real(-0.23), c = Real(1.4);
  const Real mx = Real(1.7), my = Real(0.82), K = Real(0.61);
  const Box cells(IntVect(1, 1, 1), IntVect(2, 2, 2));
  FArrayBox result(cells, 8);
  auto out = result.array();
  ParallelFor(cells, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real xi = Real(i), eta = Real(j), zeta = Real(k);
    const Real z_xlo = a * xi + b * (eta + Real(0.5)) + c * (zeta + Real(0.5));
    const Real z_xhi =
      a * (xi + Real(1.0)) + b * (eta + Real(0.5)) + c * (zeta + Real(0.5));
    const Real z_ylo = a * (xi + Real(0.5)) + b * eta + c * (zeta + Real(0.5));
    const Real z_yhi =
      a * (xi + Real(0.5)) + b * (eta + Real(1.0)) + c * (zeta + Real(0.5));
    const Real grad_xi_lo =
      Real(2.0) * xi / (mx * mx) +
      Real(2.0) * (a * xi + b * (eta + Real(0.5)) + c * (zeta + Real(0.5))) * a;
    const Real grad_xi_hi = Real(2.0) * (xi + Real(1.0)) / (mx * mx) +
                            Real(2.0) *
                              (a * (xi + Real(1.0)) + b * (eta + Real(0.5)) +
                               c * (zeta + Real(0.5))) *
                              a;
    const Real grad_eta_lo =
      Real(2.0) * eta / (my * my) +
      Real(2.0) * (a * (xi + Real(0.5)) + b * eta + c * (zeta + Real(0.5))) * b;
    const Real grad_eta_hi = Real(2.0) * (eta + Real(1.0)) / (my * my) +
                             Real(2.0) *
                               (a * (xi + Real(0.5)) + b * (eta + Real(1.0)) +
                                c * (zeta + Real(0.5))) *
                               b;
    const Real Fxlo =
      ScalarDiffusionFlux_Tx(K, mx, grad_xi_lo, a, Real(2.0) * z_xlo);
    const Real Fxhi =
      ScalarDiffusionFlux_Tx(K, mx, grad_xi_hi, a, Real(2.0) * z_xhi);
    const Real Fylo =
      ScalarDiffusionFlux_Ty(K, my, grad_eta_lo, b, Real(2.0) * z_ylo);
    const Real Fyhi =
      ScalarDiffusionFlux_Ty(K, my, grad_eta_hi, b, Real(2.0) * z_yhi);
    const Real xi_cc = xi + Real(0.5), eta_cc = eta + Real(0.5);
    const Real zeta_lo = zeta, zeta_hi = zeta + Real(1.0);
    const Real Fzlo = -Real(2.0) * K * (a * xi_cc + b * eta_cc + c * zeta_lo);
    const Real Fzhi = -Real(2.0) * K * (a * xi_cc + b * eta_cc + c * zeta_hi);
    const Real barx = Real(0.5) * (Fxlo + Fxhi);
    const Real bary = Real(0.5) * (Fylo + Fyhi);
    const Real Glo = TerrainDiffusionGz(Fzlo, mx, a, barx, my, b, bary);
    const Real Ghi = TerrainDiffusionGz(Fzhi, mx, a, barx, my, b, bary);
    const Real divergence = TerrainDiffusionDivergence_T(
      TerrainDiffusionMappedTx(Fxhi, c, my),
      TerrainDiffusionMappedTx(Fxlo, c, my),
      TerrainDiffusionMappedTy(Fyhi, c, mx),
      TerrainDiffusionMappedTy(Fylo, c, mx),
      Ghi, Glo, Real(1.0), Real(1.0), Real(1.0), mx * my, c);
    const Real Fx_rep = ScalarDiffusionFlux_Tx(
      K, mx,
      Real(2.0) * xi / (mx * mx) +
        Real(2.0) * (a * xi + b * (eta + Real(0.5)) + c * (zeta + Real(0.5))) *
          a,
      a, Real(2.0) * z_xlo);
    const Real Fy_rep = ScalarDiffusionFlux_Ty(
      K, my,
      Real(2.0) * eta / (my * my) +
        Real(2.0) * (a * (xi + Real(0.5)) + b * eta + c * (zeta + Real(0.5))) *
          b,
      b, Real(2.0) * z_ylo);
    const Real Fz_rep = -Real(2.0) * K * (a * xi_cc + b * eta_cc + c * zeta_lo);
    out(i, j, k, 0) = Fx_rep;
    out(i, j, k, 1) = Fy_rep;
    out(i, j, k, 2) = Fz_rep;
    out(i, j, k, 3) = Glo;
    out(i, j, k, 4) = Ghi;
    out(i, j, k, 5) = -divergence;
    out(i, j, k, 6) = TerrainDiffusionMappedTx(Fxlo, c, my);
    out(i, j, k, 7) = TerrainDiffusionMappedTy(Fylo, c, mx);
  });
  Gpu::streamSynchronize();
  FArrayBox host(result.box(), result.nComp(), The_Pinned_Arena());
  copy_to_host(result, host);
  const auto h = host.const_array();
  for (int k = cells.smallEnd(2); k <= cells.bigEnd(2); ++k) {
    for (int j = cells.smallEnd(1); j <= cells.bigEnd(1); ++j) {
      for (int i = cells.smallEnd(0); i <= cells.bigEnd(0); ++i) {
        const Real xi = Real(i), eta = Real(j), zeta = Real(k);
        const Real xface = xi / mx;
        const Real yface = eta / my;
        const Real Fz =
          -Real(2.0) * K *
          (a * (xi + Real(0.5)) + b * (eta + Real(0.5)) + c * zeta);
        const Real Gz = -Real(2.0) * K * c * zeta;
        EXPECT_NEAR(h(i, j, k, 0), -Real(2.0) * K * xface, tolerance(K));
        EXPECT_NEAR(h(i, j, k, 1), -Real(2.0) * K * yface, tolerance(K));
        EXPECT_NEAR(h(i, j, k, 2), Fz, tolerance(K));
        EXPECT_NEAR(h(i, j, k, 3), Gz, tolerance(K));
        const Real Gz_hi = -Real(2.0) * K * c * (zeta + Real(1.0));
        EXPECT_NEAR(h(i, j, k, 4), Gz_hi, tolerance(K));
        EXPECT_NEAR(h(i, j, k, 5), Real(6.0) * K, tolerance(K));
        EXPECT_NEAR(
          h(i, j, k, 6), -Real(2.0) * K * c * xi / (mx * my), tolerance(K));
        EXPECT_NEAR(
          h(i, j, k, 7), -Real(2.0) * K * c * eta / (mx * my), tolerance(K));
        EXPECT_GT(std::abs(Fz - Gz), tolerance(K));
      }
    }
  }
}

// Motivation: Exercise the public terrain adapter with analytic affine z and
// quadratic chi; independent RHS and face-flux oracles expose either missing
// terrain cross term rather than only testing the pointwise helpers.
ERF_GPU_TEST(
  ScalarDiffusionPrimitives,
  NativeTerrainAdapterMatchesAffineQuadraticManufacturedSolution)
{
  NativeTerrainScalarCase test(
    Real(0.31), Real(-0.23), Real(1.4), Real(1.7), Real(0.82), Real(0.61),
    RhoScalar_comp);
  test.run(Real(0.0));

  FArrayBox rhs_host(test.domain, NVAR_max, The_Pinned_Arena());
  FArrayBox xflux_host(test.xflux.box(), 1, The_Pinned_Arena());
  FArrayBox yflux_host(test.yflux.box(), 1, The_Pinned_Arena());
  FArrayBox zflux_host(test.zflux.box(), 1, The_Pinned_Arena());
  copy_to_host(test.rhs, rhs_host);
  copy_to_host(test.xflux, xflux_host);
  copy_to_host(test.yflux, yflux_host);
  copy_to_host(test.zflux, zflux_host);
  const auto rhs = rhs_host.const_array();
  const auto fx = xflux_host.const_array();
  const auto fy = yflux_host.const_array();
  const auto fz = zflux_host.const_array();

  for (int k = test.bx.smallEnd(2); k <= test.bx.bigEnd(2); ++k) {
    for (int j = test.bx.smallEnd(1); j <= test.bx.bigEnd(1); ++j) {
      for (int i = test.bx.smallEnd(0); i <= test.bx.bigEnd(0); ++i) {
        EXPECT_NEAR(
          rhs(i, j, k, RhoScalar_comp), Real(6.0) * test.K,
          integrated_terrain_tolerance(Real(6.0) * test.K));
      }
    }
  }

  const int i = 3, j = 3, k = 3;
  const Real xface = Real(i) / test.mx;
  const Real yface = Real(j) / test.my;
  const Real zface_z = test.a * (Real(i) + Real(0.5)) +
                       test.b * (Real(j) + Real(0.5)) + test.c * Real(k);
  EXPECT_NEAR(fx(i, j, k), -Real(2.0) * test.K * xface, tolerance(test.K));
  EXPECT_NEAR(fy(i, j, k), -Real(2.0) * test.K * yface, tolerance(test.K));
  EXPECT_NEAR(fz(i, j, k), -Real(2.0) * test.K * zface_z, tolerance(test.K));

  // Omitting either h_xi or h_eta from the production horizontal flux changes
  // the cell RHS by these nonzero amounts for this affine terrain.
  const Real missing_xi = Real(2.0) * test.K * test.mx * test.mx * test.a *
                          test.a;
  const Real missing_eta = Real(2.0) * test.K * test.my * test.my * test.b *
                           test.b;
  EXPECT_GT(std::abs(missing_xi), Real(100.0) * tolerance(test.K));
  EXPECT_GT(std::abs(missing_eta), Real(100.0) * tolerance(test.K));
}

// Motivation: ERF's semi-implicit terrain split scales only raw F_z. The
// lateral terrain cross terms remain explicit and must not receive that scale.
ERF_GPU_TEST(ScalarDiffusionPrimitives, TerrainImplicitSplitScalesOnlyRawVerticalFlux)
{
  const Real a = Real(0.37), b = Real(-0.29), mx = Real(1.6), my = Real(0.75);
  const Real explicit_fac = Real(0.6), implicit_fac = Real(0.4);
  const Box zbx(IntVect(0, 0, 0), IntVect(0, 0, 1));
  FArrayBox zflux(zbx, 3);
  auto z = zflux.array();
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    z(i, j, k, 0) = Real(31.0);
    z(i, j, k, 1) = -Real(0.8) - Real(0.21) * k;
    z(i, j, k, 2) = Real(-17.0);
  });
  const Real barx = Real(0.44), bary = Real(-0.63);
  const auto G_exp =
    TerrainDiffusionGz(explicit_fac * Real(-0.8), mx, a, barx, my, b, bary);
  const auto G_wrong =
    explicit_fac * TerrainDiffusionGz(Real(-0.8), mx, a, barx, my, b, bary);
  ScaleScalarDiffusionVerticalFlux(zbx, zflux.array(), 1, explicit_fac);
  Gpu::streamSynchronize();
  FArrayBox host(zbx, 3, The_Pinned_Arena());
  copy_to_host(zflux, host);
  const auto h = host.const_array();
  for (int k = 0; k <= 1; ++k) {
    const Real raw = -Real(0.8) - Real(0.21) * k;
    const Real expected = explicit_fac * raw - mx * a * barx - my * b * bary;
    const Real incorrectly_scaled =
      explicit_fac * (raw - mx * a * barx - my * b * bary);
    EXPECT_NEAR(h(0, 0, k, 1), explicit_fac * raw, tolerance());
    EXPECT_NEAR(
      TerrainDiffusionGz(h(0, 0, k, 1), mx, a, barx, my, b, bary), expected,
      tolerance());
    EXPECT_GT(std::abs(expected - incorrectly_scaled), tolerance());
    EXPECT_DOUBLE_EQ(h(0, 0, k, 0), Real(31.0));
    EXPECT_DOUBLE_EQ(h(0, 0, k, 2), Real(-17.0));
  }
  EXPECT_NEAR(G_exp, -Real(0.48) - mx * a * barx - my * b * bary, tolerance());
  EXPECT_NEAR(
    G_wrong, explicit_fac * (-Real(0.8) - mx * a * barx - my * b * bary),
    tolerance());
  EXPECT_DOUBLE_EQ(explicit_fac + implicit_fac, Real(1.0));
}

// Motivation: The public T adapter must retain complete Q1 diagnostic F_z,
// scale only raw F_z for the explicit RHS, and leave both terrain cross terms
// explicit for every implicit fraction.
ERF_GPU_TEST(ScalarDiffusionPrimitives, NativeTerrainImplicitSplitScalesOnlyRawFz)
{
  NativeTerrainScalarCase test(
    Real(0.37), Real(0.29), Real(1.25), Real(1.7), Real(0.8), Real(0.61),
    RhoQ1_comp, true);
  const int i = 3, j = 3, k = 3;
  const Real xcell = (Real(i) + Real(0.5)) / test.mx;
  const Real ycell = (Real(j) + Real(0.5)) / test.my;
  const Real xface = Real(i) / test.mx;
  const Real yface = Real(j) / test.my;
  const Real xface_z = test.a * Real(i) + test.b * (Real(j) + Real(0.5)) +
                       test.c * (Real(k) + Real(0.5));
  const Real yface_z = test.a * (Real(i) + Real(0.5)) + test.b * Real(j) +
                       test.c * (Real(k) + Real(0.5));
  const Real zface_z = test.a * (Real(i) + Real(0.5)) +
                       test.b * (Real(j) + Real(0.5)) + test.c * Real(k);
  const Real expected_x = -test.K * (Real(2.0) * xface + xface_z);
  const Real expected_y = -test.K * (Real(2.0) * yface + yface_z);
  const Real expected_raw_z =
    -test.K * (Real(2.0) * zface_z + xcell + ycell);
  ASSERT_GT(std::abs(expected_x), tolerance(test.K));
  ASSERT_GT(std::abs(expected_y), tolerance(test.K));
  ASSERT_GT(std::abs(expected_raw_z), tolerance(test.K));

  FArrayBox rhs_host(test.domain, NVAR_max, The_Pinned_Arena());
  FArrayBox xflux_host(test.xflux.box(), 1, The_Pinned_Arena());
  FArrayBox yflux_host(test.yflux.box(), 1, The_Pinned_Arena());
  FArrayBox zflux_host(test.zflux.box(), 1, The_Pinned_Arena());
  FArrayBox qfx1_host(test.qfx1_z.box(), 1, The_Pinned_Arena());
  Real first_diagnostic = Real(0.0);
  bool have_first_diagnostic = false;

  for (const Real implicit_fac : {Real(0.0), Real(0.4), Real(1.0)}) {
    const Real explicit_fac = Real(1.0) - implicit_fac;
    test.run(implicit_fac);
    copy_to_host(test.rhs, rhs_host);
    copy_to_host(test.xflux, xflux_host);
    copy_to_host(test.yflux, yflux_host);
    copy_to_host(test.zflux, zflux_host);
    copy_to_host(test.qfx1_z, qfx1_host);
    const auto rhs = rhs_host.const_array();
    const auto fx = xflux_host.const_array();
    const auto fy = yflux_host.const_array();
    const auto fz = zflux_host.const_array();
    const auto qdiag = qfx1_host.const_array();

    EXPECT_NEAR(fx(i, j, k), expected_x,
                integrated_terrain_tolerance(expected_x));
    EXPECT_NEAR(fy(i, j, k), expected_y, tolerance(test.K));
    EXPECT_NEAR(qdiag(i, j, k), expected_raw_z, tolerance(test.K));
    EXPECT_NEAR(fz(i, j, k), explicit_fac * expected_raw_z, tolerance(test.K));
    const Real expected_rhs =
      Real(4.0) * test.K + Real(2.0) * test.K * explicit_fac;
    EXPECT_NEAR(rhs(i, j, k, RhoQ1_comp), expected_rhs,
                integrated_terrain_tolerance(expected_rhs));

    if (have_first_diagnostic) {
      EXPECT_NEAR(qdiag(i, j, k), first_diagnostic, tolerance(test.K));
    } else {
      first_diagnostic = qdiag(i, j, k);
      have_first_diagnostic = true;
    }

    // Scaling all of G_z would change the analytic RHS by this nonzero amount.
    const Real incorrectly_scaled_rhs =
      Real(4.0) * test.K + Real(2.0) * test.K * explicit_fac +
      implicit_fac * test.K * (test.mx * test.a + test.my * test.b);
    if (implicit_fac > Real(0.0)) {
      EXPECT_GT(
        std::abs(incorrectly_scaled_rhs -
                 (Real(4.0) * test.K + Real(2.0) * test.K * explicit_fac)),
        Real(100.0) * tolerance(test.K));
    }
  }
}

// Native G_z construction shares the terrain interpolation and bottom/top
// extrapolation used by the fused terrain divergence.
ERF_GPU_TEST(ScalarDiffusionPrimitives, TerrainGzKFaceConstruction)
{
  const Box domain(IntVect(0, 0, 0), IntVect(4, 4, 3));
  const Box xbox = surroundingNodes(domain, 0);
  const Box ybox = surroundingNodes(domain, 1);
  const Box zbox = surroundingNodes(domain, 2);
  const Box mapbox(IntVect(0, 0, 0), IntVect(4, 4, 0));
  const Box nodes = terrain_node_box(domain);
  FArrayBox raw_x(xbox, 5), raw_y(ybox, 5), raw_z(zbox, 5);
  FArrayBox z_nd(nodes, 1), mx(mapbox, 1), my(mapbox, 1);
  FArrayBox result(Box(IntVect(0, 0, 0), IntVect(2, 0, 0)), 1);
  constexpr int raw_flux_comp = 4;
  constexpr Real a = Real(0.43), b = Real(-0.27), c = Real(1.2);
  constexpr Real mx_value = Real(1.61), my_value = Real(0.74);
  auto fx = raw_x.array();
  auto fy = raw_y.array();
  auto fz = raw_z.array();
  ParallelFor(
    xbox, raw_x.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      fx(i, j, k, n) = n == raw_flux_comp ? Real(0.8) + Real(0.11) * i -
                                              Real(0.07) * j + Real(0.16) * k
                                          : Real(700.0) + n;
    });
  ParallelFor(
    ybox, raw_y.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      fy(i, j, k, n) = n == raw_flux_comp ? Real(-0.2) + Real(0.05) * i +
                                              Real(0.19) * j - Real(0.13) * k
                                          : Real(710.0) + n;
    });
  ParallelFor(
    zbox, raw_z.nComp(),
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      fz(i, j, k, n) = n == raw_flux_comp ? Real(1.4) - Real(0.09) * i +
                                              Real(0.04) * j + Real(0.21) * k
                                          : Real(720.0) + n;
    });
  auto zn = z_nd.array();
  ParallelFor(nodes, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    zn(i, j, k) = a * i + b * j + c * k;
  });
  mx.setVal<RunOn::Device>(mx_value);
  my.setVal<RunOn::Device>(my_value);
  Gpu::streamSynchronize();

  const auto fx4 = raw_x.const_array();
  const auto fy4 = raw_y.const_array();
  const auto fz4 = raw_z.const_array();
  const auto zn4 = z_nd.const_array();
  const auto mx4 = mx.const_array();
  const auto my4 = my.const_array();
  const auto out = result.array();
  const GpuArray<Real, AMREX_SPACEDIM> inv{{Real(1.0), Real(1.0), Real(1.0)}};
  ParallelFor(result.box(), [=] AMREX_GPU_DEVICE(int i, int, int) noexcept {
    const int kface = i == 0 ? 0 : (i == 1 ? 2 : 4);
    out(i, 0, 0) = TerrainScalarDiffusionGzAtKFace(
      2, 2, kface, domain, fx4, fy4, fz4, raw_flux_comp, zn4, inv, mx4, my4);
  });
  Gpu::streamSynchronize();

  FArrayBox hfx(xbox, raw_x.nComp(), The_Pinned_Arena());
  FArrayBox hfy(ybox, raw_y.nComp(), The_Pinned_Arena());
  FArrayBox hfz(zbox, raw_z.nComp(), The_Pinned_Arena());
  FArrayBox hout(result.box(), 1, The_Pinned_Arena());
  copy_to_host(raw_x, hfx);
  copy_to_host(raw_y, hfy);
  copy_to_host(raw_z, hfz);
  copy_to_host(result, hout);
  const auto hfx4 = hfx.const_array();
  const auto hfy4 = hfy.const_array();
  const auto hfz4 = hfz.const_array();
  const auto hout4 = hout.const_array();
  const auto expected_bar = [&](const bool xdir, const int kface) {
    const auto face_avg = [&](const int k) {
      return xdir ? Real(0.5) * (hfx4(2, 2, k, raw_flux_comp) +
                                 hfx4(3, 2, k, raw_flux_comp))
                  : Real(0.5) * (hfy4(2, 2, k, raw_flux_comp) +
                                 hfy4(2, 3, k, raw_flux_comp));
    };
    if (kface == domain.smallEnd(2)) {
      return Real(1.5) * face_avg(kface) - Real(0.5) * face_avg(kface + 1);
    }
    if (kface == domain.bigEnd(2) + 1) {
      return Real(1.5) * face_avg(kface - 1) - Real(0.5) * face_avg(kface - 2);
    }
    const auto four_face_avg = [&](const int k) {
      return xdir ? Real(0.5) * (hfx4(2, 2, k, raw_flux_comp) +
                                 hfx4(3, 2, k, raw_flux_comp))
                  : Real(0.5) * (hfy4(2, 2, k, raw_flux_comp) +
                                 hfy4(2, 3, k, raw_flux_comp));
    };
    return Real(0.5) * (four_face_avg(kface) + four_face_avg(kface - 1));
  };
  for (int n = 0; n < 3; ++n) {
    const int kface = n == 0 ? 0 : (n == 1 ? 2 : 4);
    const Real raw = hfz4(2, 2, kface, raw_flux_comp);
    const Real expected = raw - a * mx_value * expected_bar(true, kface) -
                          b * my_value * expected_bar(false, kface);
    EXPECT_NEAR(hout4(n, 0, 0), expected, tolerance(expected));
    EXPECT_GT(std::abs(raw - expected), Real(100.0) * tolerance(raw));
  }
}

// Fused/chunked callers can use the same raw terrain face stencils as the
// field builder, with independent scalar, density, and output components.
ERF_GPU_TEST(ScalarDiffusionPrimitives, TerrainPointwiseFaceFluxPrimitives)
{
  constexpr int scalar_comp = 2;
  constexpr int rho_comp = 1;
  constexpr int raw_flux_comp = 4;
  constexpr Real a = Real(0.31), b = Real(-0.22), c = Real(1.3);
  constexpr Real mx_value = Real(1.7), my_value = Real(0.81);
  constexpr Real molecular = Real(0.43);
  const Box domain(IntVect(0, 0, 0), IntVect(4, 4, 3));
  Box data_box = domain;
  data_box.grow(2);
  const Box nodes = terrain_node_box(data_box);
  const Box mapbox(IntVect(0, 0, 0), IntVect(4, 4, 0));
  FArrayBox scalar(data_box, 4), rho(data_box, 2), z_nd(nodes, 1),
    z_cc(data_box, 1);
  FArrayBox mf_ux(surroundingNodes(mapbox, 0), 1),
    mf_vy(surroundingNodes(mapbox, 1), 1);
  FArrayBox xflux(surroundingNodes(domain, 0), 5),
    yflux(surroundingNodes(domain, 1), 5),
    zflux(surroundingNodes(domain, 2), 5);

  auto s = scalar.array();
  auto r = rho.array();
  ParallelFor(data_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real x = (Real(i) + Real(0.5)) / mx_value;
    const Real y = (Real(j) + Real(0.5)) / my_value;
    const Real z = a * (Real(i) + Real(0.5)) + b * (Real(j) + Real(0.5)) +
                   c * (Real(k) + Real(0.5));
    s(i, j, k, scalar_comp) = x + Real(2.0) * y + Real(3.0) * z;
    r(i, j, k, rho_comp) =
      Real(1.2) + Real(0.03) * i + Real(0.02) * j + Real(0.04) * k;
  });
  auto zn = z_nd.array();
  ParallelFor(nodes, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    zn(i, j, k) = a * i + b * j + c * k;
  });
  auto zc = z_cc.array();
  ParallelFor(data_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    zc(i, j, k) = a * (Real(i) + Real(0.5)) + b * (Real(j) + Real(0.5)) +
                  c * (Real(k) + Real(0.5));
  });
  mf_ux.setVal<RunOn::Device>(Real(1.27));
  mf_vy.setVal<RunOn::Device>(Real(0.83));
  xflux.setVal<RunOn::Device>(Real(600.0));
  yflux.setVal<RunOn::Device>(Real(610.0));
  zflux.setVal<RunOn::Device>(Real(620.0));
  Gpu::streamSynchronize();

  const GpuArray<Real, AMREX_SPACEDIM> inv{{Real(1.0), Real(1.0), Real(1.0)}};
  const auto scalar4 = scalar.const_array();
  const auto rho4 = rho.const_array();
  const auto zn4 = z_nd.const_array();
  const auto zc4 = z_cc.const_array();
  const auto ux4 = mf_ux.const_array();
  const auto vy4 = mf_vy.const_array();
  auto fx = xflux.array();
  auto fy = yflux.array();
  auto fz = zflux.array();
  const Box iface(IntVect(2, 2, 1), IntVect(2, 2, 1));
  const Box jface(IntVect(2, 2, 1), IntVect(2, 2, 1));
  ParallelFor(iface, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real coeff = molecular * Real(0.5) *
                       (rho4(i, j, k, rho_comp) + rho4(i - 1, j, k, rho_comp));
    fx(i, j, k, raw_flux_comp) = TerrainScalarDiffusionFluxAtIFace(
      scalar4, scalar_comp, coeff, i, j, k, zn4, zc4, inv, ux4);
  });
  ParallelFor(jface, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real coeff = molecular * Real(0.5) *
                       (rho4(i, j, k, rho_comp) + rho4(i, j - 1, k, rho_comp));
    fy(i, j, k, raw_flux_comp) = TerrainScalarDiffusionFluxAtJFace(
      scalar4, scalar_comp, coeff, i, j, k, zn4, zc4, inv, vy4);
  });
  const Box interior_face(IntVect(2, 2, 2), IntVect(2, 2, 2));
  const Box low_face(IntVect(2, 2, 0), IntVect(2, 2, 0));
  const Box high_face(IntVect(2, 2, 4), IntVect(2, 2, 4));
  ParallelFor(
    interior_face, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const Real coeff =
        molecular * Real(0.5) *
        (rho4(i, j, k, rho_comp) + rho4(i, j, k - 1, rho_comp));
      fz(i, j, k, raw_flux_comp) = TerrainScalarDiffusionFluxAtKFace(
        scalar4, scalar_comp, coeff, i, j, k, zn4, inv, false, false);
    });
  ParallelFor(low_face, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real coeff = molecular * Real(0.5) *
                       (rho4(i, j, k, rho_comp) + rho4(i, j, k - 1, rho_comp));
    fz(i, j, k, raw_flux_comp) = TerrainScalarDiffusionFluxAtKFace(
      scalar4, scalar_comp, coeff, i, j, k, zn4, inv, true, false);
  });
  ParallelFor(high_face, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real coeff = molecular * Real(0.5) *
                       (rho4(i, j, k, rho_comp) + rho4(i, j, k - 1, rho_comp));
    fz(i, j, k, raw_flux_comp) = TerrainScalarDiffusionFluxAtKFace(
      scalar4, scalar_comp, coeff, i, j, k, zn4, inv, false, true);
  });
  Gpu::streamSynchronize();

  FArrayBox hfx(xflux.box(), 5, The_Pinned_Arena());
  FArrayBox hfy(yflux.box(), 5, The_Pinned_Arena());
  FArrayBox hfz(zflux.box(), 5, The_Pinned_Arena());
  copy_to_host(xflux, hfx);
  copy_to_host(yflux, hfy);
  copy_to_host(zflux, hfz);
  const auto hfx4 = hfx.const_array();
  const auto hfy4 = hfy.const_array();
  const auto hfz4 = hfz.const_array();
  const auto rho_at = [](int i, int j, int k) {
    return Real(1.2) + Real(0.03) * i + Real(0.02) * j + Real(0.04) * k;
  };
  const Real coeff_x =
    molecular * Real(0.5) * (rho_at(2, 2, 1) + rho_at(1, 2, 1));
  const Real coeff_y =
    molecular * Real(0.5) * (rho_at(2, 2, 1) + rho_at(2, 1, 1));
  const Real expected_x = -coeff_x * Real(1.27) / mx_value;
  const Real expected_y = -coeff_y * Real(0.83) * Real(2.0) / my_value;
  EXPECT_NEAR(hfx4(2, 2, 1, raw_flux_comp), expected_x, tolerance(expected_x));
  EXPECT_NEAR(hfy4(2, 2, 1, raw_flux_comp), expected_y, tolerance(expected_y));

  const auto scalar_at = [&](int i, int j, int k) {
    const Real x = (Real(i) + Real(0.5)) / mx_value;
    const Real y = (Real(j) + Real(0.5)) / my_value;
    const Real z = a * (Real(i) + Real(0.5)) + b * (Real(j) + Real(0.5)) +
                   c * (Real(k) + Real(0.5));
    return x + Real(2.0) * y + Real(3.0) * z;
  };
  const auto low_grad = [&]() {
    const Real f = Real(3.0);
    const Real f2 = f * f;
    const Real c3 = Real(2.0) / (f - f2);
    const Real c2 = -f2 * c3;
    const Real c1 = -(Real(1.0) - f2) * c3;
    return (c1 * scalar_at(2, 2, -1) + c2 * scalar_at(2, 2, 0) +
            c3 * scalar_at(2, 2, 1)) /
           c;
  };
  const auto high_grad = [&]() {
    const Real f = Real(3.0);
    const Real f2 = f * f;
    const Real c3 = Real(2.0) / (f - f2);
    const Real c2 = -f2 * c3;
    const Real c1 = -(Real(1.0) - f2) * c3;
    return -(c1 * scalar_at(2, 2, 4) + c2 * scalar_at(2, 2, 3) +
             c3 * scalar_at(2, 2, 2)) /
           c;
  };
  const Real coeff_zlo =
    molecular * Real(0.5) * (rho_at(2, 2, 0) + rho_at(2, 2, -1));
  const Real coeff_zmid =
    molecular * Real(0.5) * (rho_at(2, 2, 2) + rho_at(2, 2, 1));
  const Real coeff_zhi =
    molecular * Real(0.5) * (rho_at(2, 2, 4) + rho_at(2, 2, 3));
  EXPECT_NEAR(
    hfz4(2, 2, 2, raw_flux_comp), -coeff_zmid * Real(3.0), tolerance());
  EXPECT_NEAR(
    hfz4(2, 2, 0, raw_flux_comp), -coeff_zlo * low_grad(), tolerance());
  EXPECT_NEAR(
    hfz4(2, 2, 4, raw_flux_comp), -coeff_zhi * high_grad(), tolerance());
  for (int n = 0; n < 5; ++n) {
    if (n != raw_flux_comp) {
      EXPECT_DOUBLE_EQ(hfx4(2, 2, 1, n), Real(600.0));
      EXPECT_DOUBLE_EQ(hfy4(2, 2, 1, n), Real(610.0));
      EXPECT_DOUBLE_EQ(hfz4(2, 2, 2, n), Real(620.0));
      EXPECT_DOUBLE_EQ(hfz4(2, 2, 0, n), Real(620.0));
      EXPECT_DOUBLE_EQ(hfz4(2, 2, 4, n), Real(620.0));
    }
  }
}

#undef ERF_GPU_TEST

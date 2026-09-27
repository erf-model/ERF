#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include <ERF_Diffusion.H>
#include <ERF_NativeScalarDiffusion.H>
#include <ERF_AdvectionSrcForScalars.H>

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

// The mapped API retains the same component-explicit builder inputs while
// moving only the completed face transfers to scalar advection's flux space.
ERF_GPU_TEST(ScalarDiffusionPrimitives, CanonicalNTransfersMatchNativeDivergence)
{
  constexpr int raw_comp = 3, mapped_comp = 4, rhs_comp = 5;
  const Box domain(IntVect(0, 0, 0), IntVect(3, 3, 3));
  const Box bx(IntVect(1, 1, 1), IntVect(2, 2, 2));
  Box data_box = domain;
  data_box.grow(1);
  const Box map_cells(IntVect(0, 0, 0), IntVect(3, 3, 0));
  const Box xfaces = surroundingNodes(domain, 0);
  const Box yfaces = surroundingNodes(domain, 1);
  const Box zfaces = surroundingNodes(domain, 2);
  FArrayBox scalar(data_box, 4), rho(data_box, 2), mu(data_box, EddyDiff::NumDiffs);
  FArrayBox u(xfaces, 1), v(yfaces, 1);
  FArrayBox raw_x(xfaces, 5), raw_y(yfaces, 5), raw_z(zfaces, 5);
  FArrayBox mapped_x(xfaces, 5), mapped_y(yfaces, 5), mapped_z(zfaces, 5);
  FArrayBox wrong_x(xfaces, 3), native_rhs(bx, 6), mapped_rhs(bx, 6),
    wrong_rhs(bx, 6), detj(domain, 1);
  FArrayBox mf_mx(map_cells, 1), mf_my(map_cells, 1);
  FArrayBox mf_ux(surroundingNodes(map_cells, 0), 1),
    mf_uy(surroundingNodes(map_cells, 0), 1),
    mf_vx(surroundingNodes(map_cells, 1), 1),
    mf_vy(surroundingNodes(map_cells, 1), 1);

  auto s = scalar.array();
  auto r = rho.array();
  auto m = mu.array();
  ParallelFor(data_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    s(i, j, k, 2) = Real(0.17) * i * i + Real(0.09) * j * j +
                    Real(0.13) * k * k + Real(0.04) * i * j +
                    Real(0.03) * j * k;
    r(i, j, k, 1) = Real(0.91) + Real(0.03) * i + Real(0.02) * j +
                    Real(0.015) * k;
    for (int n = 0; n < EddyDiff::NumDiffs; ++n) {
      m(i, j, k, n) = Real(0.08) + Real(0.002) * i + Real(0.003) * j +
                      Real(0.004) * k + Real(0.001) * n;
    }
  });
  auto mx = mf_mx.array();
  auto my = mf_my.array();
  ParallelFor(map_cells, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
    mx(i, j, 0) = Real(1.2) + Real(0.025) * i + Real(0.012) * j;
    my(i, j, 0) = Real(0.79) + Real(0.019) * j + Real(0.006) * i;
  });
  auto ux = mf_ux.array();
  auto uy = mf_uy.array();
  ParallelFor(mf_ux.box(), [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
    ux(i, j, 0) = Real(1.08) + Real(0.021) * i;
    uy(i, j, 0) = Real(0.87) + Real(0.014) * j;
  });
  auto vx = mf_vx.array();
  auto vy = mf_vy.array();
  ParallelFor(mf_vx.box(), [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
    vx(i, j, 0) = Real(0.93) + Real(0.012) * i;
    vy(i, j, 0) = Real(1.04) + Real(0.018) * j;
  });
  u.setVal<RunOn::Device>(Real(0.0));
  v.setVal<RunOn::Device>(Real(0.0));
  raw_x.setVal<RunOn::Device>(Real(-101.0));
  raw_y.setVal<RunOn::Device>(Real(-102.0));
  raw_z.setVal<RunOn::Device>(Real(-103.0));
  mapped_x.setVal<RunOn::Device>(Real(-201.0));
  mapped_y.setVal<RunOn::Device>(Real(-202.0));
  mapped_z.setVal<RunOn::Device>(Real(-203.0));
  native_rhs.setVal<RunOn::Device>(Real(-301.0));
  mapped_rhs.setVal<RunOn::Device>(Real(-301.0));
  wrong_rhs.setVal<RunOn::Device>(Real(-301.0));
  detj.setVal<RunOn::Device>(Real(1.0));
  const auto native_rhs_init = native_rhs.array();
  const auto mapped_rhs_init = mapped_rhs.array();
  const auto wrong_rhs_init = wrong_rhs.array();
  ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    native_rhs_init(i, j, k, rhs_comp) = Real(0.0);
    mapped_rhs_init(i, j, k, rhs_comp) = Real(0.0);
    wrong_rhs_init(i, j, k, rhs_comp) = Real(0.0);
  });
  Gpu::streamSynchronize();

  ScalarDiffusionFieldViews field;
  field.scalar = scalar.const_array();
  field.scalar_comp = 2;
  field.density = rho.const_array();
  field.rho_comp = 1;
  field.mu_turb = mu.const_array();
  field.xflux = raw_x.array();
  field.yflux = raw_y.array();
  field.zflux = raw_z.array();
  field.flux_comp = raw_comp;
  field.rhs = native_rhs.array();
  field.rhs_comp = rhs_comp;
  ScalarDiffusionFluxPolicy policy;
  policy.coefficients = {Real(0.29), EddyDiff::Scalar_h, EddyDiff::Scalar_v};
  policy.coefficient_mode = {true, true};

  Vector<BCRec> bcs(NBCVAR_max);
  for (auto& bc : bcs) {
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      bc.setLo(d, ERFBCType::foextrap);
      bc.setHi(d, ERFBCType::foextrap);
    }
  }
  Gpu::DeviceVector<BCRec> bcs_device(bcs.size());
  Gpu::copy(Gpu::hostToDevice, bcs.begin(), bcs.end(), bcs_device.begin());
  const GpuArray<Real, AMREX_SPACEDIM> inv{
    {Real(0.73), Real(0.91), Real(1.17)}};
  BuildScalarDiffusionFluxes_N(
    bx, domain, field, policy, inv, u.const_array(), v.const_array(),
    mf_ux.const_array(), mf_uy.const_array(), mf_vx.const_array(),
    mf_vy.const_array(), bcs_device.data(), 0);
  ApplyScalarDiffusionFluxDivergence_N(
    bx, raw_x.const_array(), raw_y.const_array(), raw_z.const_array(), raw_comp,
    native_rhs.array(), rhs_comp, inv, mf_mx.const_array(), mf_my.const_array());
  BuildScalarDiffusionMappedTransfers_N(
    bx, raw_x.const_array(), raw_comp, raw_y.const_array(), raw_comp,
    raw_z.const_array(), raw_comp, mapped_x.array(), mapped_comp,
    mapped_y.array(), mapped_comp, mapped_z.array(), mapped_comp,
    mf_mx.const_array(), mf_my.const_array());
  ApplyScalarMappedFluxDivergence(
    bx, mapped_x.const_array(), mapped_comp, mapped_y.const_array(), mapped_comp,
    mapped_z.const_array(), mapped_comp, mapped_rhs.array(), rhs_comp,
    detj.const_array(), inv, mf_mx.const_array(), mf_my.const_array());

  const auto fx = raw_x.const_array();
  const auto badx = wrong_x.array();
  const auto ux4 = mf_ux.const_array();
  ParallelFor(xfaces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    badx(i, j, k, 2) = fx(i, j, k, raw_comp) * ux4(i, j, 0);
  });
  ApplyScalarMappedFluxDivergence(
    bx, wrong_x.const_array(), 2, mapped_y.const_array(), mapped_comp,
    mapped_z.const_array(), mapped_comp, wrong_rhs.array(), rhs_comp,
    detj.const_array(), inv, mf_mx.const_array(), mf_my.const_array());
  Gpu::streamSynchronize();

  FArrayBox hn(bx, 6, The_Pinned_Arena()), hm(bx, 6, The_Pinned_Arena()),
    hw(bx, 6, The_Pinned_Arena());
  copy_to_host(native_rhs, hn);
  copy_to_host(mapped_rhs, hm);
  copy_to_host(wrong_rhs, hw);
  const auto n = hn.const_array();
  const auto c = hm.const_array();
  const auto w = hw.const_array();
  for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
    for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
      for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
        EXPECT_NEAR(n(i, j, k, rhs_comp), c(i, j, k, rhs_comp),
                    integrated_terrain_tolerance(n(i, j, k, rhs_comp)));
        EXPECT_GT(std::abs(n(i, j, k, rhs_comp) - w(i, j, k, rhs_comp)),
                  Real(100.0) * tolerance());
        for (int comp = 0; comp < 6; ++comp) {
          if (comp != rhs_comp) {
            EXPECT_DOUBLE_EQ(n(i, j, k, comp), Real(-301.0));
            EXPECT_DOUBLE_EQ(c(i, j, k, comp), Real(-301.0));
            EXPECT_DOUBLE_EQ(w(i, j, k, comp), Real(-301.0));
          }
        }
      }
    }
  }
}

ERF_GPU_TEST(
  ScalarDiffusionPrimitives,
  CanonicalSTransfersRestoreNonuniformSideArea)
{
  constexpr int raw_comp = 2, mapped_comp = 3, rhs_comp = 4;
  const Box domain(IntVect(0, 0, 0), IntVect(3, 3, 3));
  const Box bx(IntVect(1, 1, 1), IntVect(2, 2, 2));
  Box data_box = domain;
  data_box.grow(1);
  const Box map_cells(IntVect(0, 0, 0), IntVect(3, 3, 0));
  const Box xfaces = surroundingNodes(domain, 0);
  const Box yfaces = surroundingNodes(domain, 1);
  const Box zfaces = surroundingNodes(domain, 2);
  FArrayBox scalar(data_box, 4), rho(data_box, 2), mu(data_box, EddyDiff::NumDiffs);
  FArrayBox raw_x(xfaces, 4), raw_y(yfaces, 4), raw_z(zfaces, 4);
  FArrayBox mapped_x(xfaces, 4), mapped_y(yfaces, 4), mapped_z(zfaces, 4);
  FArrayBox no_area_x(xfaces, 4), no_area_y(yfaces, 4), native_rhs(bx, 5),
    mapped_rhs(bx, 5), no_area_rhs(bx, 5), detj(domain, 1);
  FArrayBox ax(xfaces, 1), ay(yfaces, 1), mf_mx(map_cells, 1), mf_my(map_cells, 1);
  FArrayBox mf_ux(surroundingNodes(map_cells, 0), 1),
    mf_uy(surroundingNodes(map_cells, 0), 1),
    mf_vx(surroundingNodes(map_cells, 1), 1),
    mf_vy(surroundingNodes(map_cells, 1), 1);
  Vector<Real> dz_host{Real(0.5), Real(1.2), Real(0.8), Real(1.5)};
  Gpu::DeviceVector<Real> dz(4);
  Gpu::copy(Gpu::hostToDevice, dz_host.begin(), dz_host.end(), dz.begin());
  const auto dzp = dz.data();

  auto s = scalar.array();
  auto r = rho.array();
  auto m = mu.array();
  ParallelFor(data_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    Real z;
    if (k < 0) {
      z = (Real(k) + Real(0.5)) * dzp[0];
    } else if (k >= 4) {
      z = dzp[0] + dzp[1] + dzp[2] + dzp[3] +
          (Real(k - 4) + Real(0.5)) * dzp[3];
    } else {
      z = Real(0.5) * dzp[k];
      for (int q = 0; q < k; ++q) z += dzp[q];
    }
    s(i, j, k, 2) = Real(0.12) * i * i + Real(0.08) * j * j +
                    Real(0.21) * z + Real(0.035) * z * z +
                    Real(0.025) * i * j;
    r(i, j, k, 1) = Real(0.95) + Real(0.02) * i + Real(0.01) * j +
                    Real(0.015) * k;
    for (int n = 0; n < EddyDiff::NumDiffs; ++n) {
      m(i, j, k, n) = Real(0.06) + Real(0.002) * i + Real(0.003) * j +
                      Real(0.004) * k + Real(0.001) * n;
    }
  });
  auto mx = mf_mx.array();
  auto my = mf_my.array();
  ParallelFor(map_cells, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
    mx(i, j, 0) = Real(1.16) + Real(0.02) * i + Real(0.011) * j;
    my(i, j, 0) = Real(0.83) + Real(0.017) * j + Real(0.005) * i;
  });
  auto ux = mf_ux.array();
  auto uy = mf_uy.array();
  ParallelFor(mf_ux.box(), [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
    ux(i, j, 0) = Real(1.06) + Real(0.018) * i;
    uy(i, j, 0) = Real(0.91) + Real(0.012) * j;
  });
  auto vx = mf_vx.array();
  auto vy = mf_vy.array();
  ParallelFor(mf_vx.box(), [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
    vx(i, j, 0) = Real(0.94) + Real(0.014) * i;
    vy(i, j, 0) = Real(1.03) + Real(0.016) * j;
  });
  auto ax4 = ax.array();
  ParallelFor(xfaces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    ax4(i, j, k) = dzp[k];
  });
  auto ay4 = ay.array();
  ParallelFor(yfaces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    ay4(i, j, k) = dzp[k];
  });
  auto det = detj.array();
  ParallelFor(domain, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    det(i, j, k) = dzp[k];
  });
  raw_x.setVal<RunOn::Device>(Real(-401.0));
  raw_y.setVal<RunOn::Device>(Real(-402.0));
  raw_z.setVal<RunOn::Device>(Real(-403.0));
  mapped_x.setVal<RunOn::Device>(Real(-501.0));
  mapped_y.setVal<RunOn::Device>(Real(-502.0));
  mapped_z.setVal<RunOn::Device>(Real(-503.0));
  no_area_x.setVal<RunOn::Device>(Real(-601.0));
  no_area_y.setVal<RunOn::Device>(Real(-602.0));
  native_rhs.setVal<RunOn::Device>(Real(0.0));
  mapped_rhs.setVal<RunOn::Device>(Real(0.0));
  no_area_rhs.setVal<RunOn::Device>(Real(0.0));
  Gpu::streamSynchronize();

  ScalarDiffusionFieldViews field;
  field.scalar = scalar.const_array();
  field.scalar_comp = 2;
  field.density = rho.const_array();
  field.rho_comp = 1;
  field.mu_turb = mu.const_array();
  field.xflux = raw_x.array();
  field.yflux = raw_y.array();
  field.zflux = raw_z.array();
  field.flux_comp = raw_comp;
  field.rhs = native_rhs.array();
  field.rhs_comp = rhs_comp;
  ScalarDiffusionFluxPolicy policy;
  policy.coefficients = {Real(0.24), EddyDiff::Scalar_h, EddyDiff::Scalar_v};
  policy.coefficient_mode = {true, true};
  Vector<BCRec> bcs(NBCVAR_max);
  for (auto& bc : bcs) {
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      bc.setLo(d, ERFBCType::foextrap);
      bc.setHi(d, ERFBCType::foextrap);
    }
  }
  Gpu::DeviceVector<BCRec> bcs_device(bcs.size());
  Gpu::copy(Gpu::hostToDevice, bcs.begin(), bcs.end(), bcs_device.begin());
  constexpr Real dx_inv = Real(0.74), dy_inv = Real(0.89);
  BuildScalarDiffusionFluxes_S(
    bx, domain, field, policy, dx_inv, dy_inv, dzp, 0, 3,
    mf_ux.const_array(), mf_uy.const_array(), mf_vx.const_array(),
    mf_vy.const_array(), bcs_device.data(), 0);
  ApplyScalarDiffusionFluxDivergence_S(
    bx, raw_x.const_array(), raw_y.const_array(), raw_z.const_array(), raw_comp,
    native_rhs.array(), rhs_comp, mf_mx.const_array(), mf_my.const_array(),
    dx_inv, dy_inv, dzp);
  BuildScalarDiffusionMappedTransfers_S(
    bx, raw_x.const_array(), raw_comp, raw_y.const_array(), raw_comp,
    raw_z.const_array(), raw_comp, mapped_x.array(), mapped_comp,
    mapped_y.array(), mapped_comp, mapped_z.array(), mapped_comp,
    ax.const_array(), ay.const_array(), mf_mx.const_array(), mf_my.const_array());
  const GpuArray<Real, AMREX_SPACEDIM> inv{{dx_inv, dy_inv, Real(1.0)}};
  ApplyScalarMappedFluxDivergence(
    bx, mapped_x.const_array(), mapped_comp, mapped_y.const_array(), mapped_comp,
    mapped_z.const_array(), mapped_comp, mapped_rhs.array(), rhs_comp,
    detj.const_array(), inv, mf_mx.const_array(), mf_my.const_array());
  FArrayBox unit_ax(xfaces, 1), unit_ay(yfaces, 1);
  unit_ax.setVal<RunOn::Device>(Real(1.0));
  unit_ay.setVal<RunOn::Device>(Real(1.0));
  BuildScalarDiffusionMappedTransfers_S(
    bx, raw_x.const_array(), raw_comp, raw_y.const_array(), raw_comp,
    raw_z.const_array(), raw_comp, no_area_x.array(), mapped_comp,
    no_area_y.array(), mapped_comp, mapped_z.array(), mapped_comp,
    unit_ax.const_array(), unit_ay.const_array(), mf_mx.const_array(),
    mf_my.const_array());
  ApplyScalarMappedFluxDivergence(
    bx, no_area_x.const_array(), mapped_comp, no_area_y.const_array(),
    mapped_comp, mapped_z.const_array(), mapped_comp, no_area_rhs.array(),
    rhs_comp, detj.const_array(), inv, mf_mx.const_array(), mf_my.const_array());
  Gpu::streamSynchronize();

  FArrayBox hn(bx, 5, The_Pinned_Arena()), hm(bx, 5, The_Pinned_Arena()),
    hw(bx, 5, The_Pinned_Arena());
  copy_to_host(native_rhs, hn);
  copy_to_host(mapped_rhs, hm);
  copy_to_host(no_area_rhs, hw);
  const auto n = hn.const_array();
  const auto c = hm.const_array();
  const auto w = hw.const_array();
  for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
    for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
      for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
        EXPECT_NEAR(n(i, j, k, rhs_comp), c(i, j, k, rhs_comp),
                    integrated_terrain_tolerance(n(i, j, k, rhs_comp)));
        EXPECT_GT(std::abs(n(i, j, k, rhs_comp) - w(i, j, k, rhs_comp)),
                  Real(100.0) * tolerance());
      }
    }
  }
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
  FArrayBox result(Box(IntVect(0, 0, 0), IntVect(2, 0, 0)), 3);
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
    out(i, 0, 0, 1) = TerrainScalarDiffusionMappedTransferAtKFace(
      2, 2, kface, domain, fx4, fy4, fz4, raw_flux_comp, zn4, inv, mx4, my4);
    out(i, 0, 0, 2) = TerrainScalarDiffusionMappedTransferAtKFace(
      2, 2, kface, domain, fx4, fy4, fz4, raw_flux_comp, zn4, inv, mx4, my4,
      kface == domain.smallEnd(2));
  });
  Gpu::streamSynchronize();

  FArrayBox hfx(xbox, raw_x.nComp(), The_Pinned_Arena());
  FArrayBox hfy(ybox, raw_y.nComp(), The_Pinned_Arena());
  FArrayBox hfz(zbox, raw_z.nComp(), The_Pinned_Arena());
  FArrayBox hout(result.box(), 3, The_Pinned_Arena());
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
    EXPECT_NEAR(hout4(n, 0, 0, 1), expected / (mx_value * my_value),
                tolerance(expected / (mx_value * my_value)));
    EXPECT_DOUBLE_EQ(
      hout4(n, 0, 0, 2), n == 0 ? Real(0.0) : hout4(n, 0, 0, 1));
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

ERF_GPU_TEST(
  ScalarDiffusionPrimitives,
  CanonicalTerrainTransfersMatchNativeAndAcceptedValues)
{
  constexpr int mapped_comp = 2, accepted_comp = 1;
  constexpr int rhs_comp = 4;
  NativeTerrainScalarCase test(
    Real(0.31), Real(-0.23), Real(1.4), Real(1.7), Real(0.82), Real(0.61),
    RhoScalar_comp, true);
  test.run(Real(0.0));

  FArrayBox mapped_x(test.xflux.box(), 3), mapped_y(test.yflux.box(), 3),
    mapped_z(test.zflux.box(), 3),
    raw_z_transfer(test.zflux.box(), 3),
    missing_cross_transfer(test.zflux.box(), 3),
    accepted_x(test.xflux.box(), 3), accepted_y(test.yflux.box(), 3),
    accepted_z(test.zflux.box(), 3);
  FArrayBox canonical(test.bx, NVAR_max), raw_wrong(test.bx, 5),
    cross_wrong(test.bx, 5), accepted_rhs(test.bx, 5);
  for (FArrayBox* f : {&mapped_x, &mapped_y, &mapped_z})
    f->setVal<RunOn::Device>(Real(-701.0));
  raw_z_transfer.setVal<RunOn::Device>(Real(-702.0));
  missing_cross_transfer.setVal<RunOn::Device>(Real(-703.0));
  accepted_x.setVal<RunOn::Device>(Real(-704.0));
  accepted_y.setVal<RunOn::Device>(Real(-705.0));
  accepted_z.setVal<RunOn::Device>(Real(-706.0));
  canonical.setVal<RunOn::Device>(Real(0.0));
  raw_wrong.setVal<RunOn::Device>(Real(0.0));
  cross_wrong.setVal<RunOn::Device>(Real(0.0));
  accepted_rhs.setVal<RunOn::Device>(Real(0.0));
  Gpu::streamSynchronize();

  BuildScalarDiffusionMappedTransfers_T(
    test.bx, test.domain, test.xflux.const_array(), test.yflux.const_array(),
    test.zflux.const_array(), 0, mapped_x.array(), mapped_comp,
    mapped_y.array(), mapped_comp, mapped_z.array(), mapped_comp,
    test.z_nd.const_array(), test.ax.const_array(), test.ay.const_array(),
    test.inv, test.mf_mx.const_array(), test.mf_uy.const_array(),
    test.mf_my.const_array(), test.mf_vx.const_array(), false);
  ApplyScalarMappedFluxDivergence(
    test.bx, mapped_x.const_array(), mapped_comp, mapped_y.const_array(),
    mapped_comp, mapped_z.const_array(), mapped_comp, canonical.array(),
    RhoScalar_comp, test.detj.const_array(), test.inv, test.mf_mx.const_array(),
    test.mf_my.const_array());

  const Box zbx = surroundingNodes(test.bx, 2);
  const auto fz = test.zflux.const_array();
  const auto mx = test.mf_mx.const_array();
  const auto my = test.mf_my.const_array();
  const auto raw_only = raw_z_transfer.array();
  const auto no_cross = missing_cross_transfer.array();
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    raw_only(i, j, k, mapped_comp) = fz(i, j, k, 0);
    no_cross(i, j, k, mapped_comp) =
      ScalarDiffusionMappedVerticalTransfer(
        fz(i, j, k, 0), mx(i, j, 0), my(i, j, 0));
  });
  ApplyScalarMappedFluxDivergence(
    test.bx, mapped_x.const_array(), mapped_comp, mapped_y.const_array(),
    mapped_comp, raw_z_transfer.const_array(), mapped_comp, raw_wrong.array(),
    rhs_comp, test.detj.const_array(), test.inv, test.mf_mx.const_array(),
    test.mf_my.const_array());
  ApplyScalarMappedFluxDivergence(
    test.bx, mapped_x.const_array(), mapped_comp, mapped_y.const_array(),
    mapped_comp, missing_cross_transfer.const_array(), mapped_comp,
    cross_wrong.array(), rhs_comp, test.detj.const_array(), test.inv,
    test.mf_mx.const_array(), test.mf_my.const_array());

  constexpr Real lambda = Real(0.37);
  const Real low_weight = Real(1.0) - lambda;
  const auto hx = mapped_x.const_array();
  const auto hy = mapped_y.const_array();
  const auto hz = mapped_z.const_array();
  const auto ax = accepted_x.array();
  const auto ay = accepted_y.array();
  const auto az = accepted_z.array();
  ParallelFor(accepted_x.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real delta = Real(0.043) * i - Real(0.021) * j + Real(0.017) * k;
    const Real low = hx(i, j, k, mapped_comp) - delta;
    ax(i, j, k, accepted_comp) = low + lambda * delta;
  });
  ParallelFor(accepted_y.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real delta = -Real(0.032) * i + Real(0.029) * j + Real(0.013) * k;
    const Real low = hy(i, j, k, mapped_comp) - delta;
    ay(i, j, k, accepted_comp) = low + lambda * delta;
  });
  ParallelFor(accepted_z.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const Real delta = Real(0.019) * i - Real(0.027) * j + Real(0.041) * k;
    const Real low = hz(i, j, k, mapped_comp) - delta;
    az(i, j, k, accepted_comp) = low + lambda * delta;
  });
  ApplyScalarMappedFluxDivergence(
    test.bx, accepted_x.const_array(), accepted_comp, accepted_y.const_array(),
    accepted_comp, accepted_z.const_array(), accepted_comp, accepted_rhs.array(),
    rhs_comp, test.detj.const_array(), test.inv, test.mf_mx.const_array(),
    test.mf_my.const_array());
  Gpu::streamSynchronize();

  FArrayBox hn(test.domain, NVAR_max, The_Pinned_Arena()),
    hc(test.bx, NVAR_max, The_Pinned_Arena()),
    hr(test.bx, 5, The_Pinned_Arena()), hw(test.bx, 5, The_Pinned_Arena()),
    ha(test.bx, 5, The_Pinned_Arena());
  FArrayBox host_x(accepted_x.box(), 3, The_Pinned_Arena()),
    host_y(accepted_y.box(), 3, The_Pinned_Arena()),
    host_z(accepted_z.box(), 3, The_Pinned_Arena()),
    hdet(test.detj.box(), 1, The_Pinned_Arena()),
    hmx(test.mf_mx.box(), 1, The_Pinned_Arena()),
    hmy(test.mf_my.box(), 1, The_Pinned_Arena());
  copy_to_host(test.rhs, hn);
  copy_to_host(canonical, hc);
  copy_to_host(raw_wrong, hr);
  copy_to_host(cross_wrong, hw);
  copy_to_host(accepted_rhs, ha);
  copy_to_host(accepted_x, host_x);
  copy_to_host(accepted_y, host_y);
  copy_to_host(accepted_z, host_z);
  copy_to_host(test.detj, hdet);
  copy_to_host(test.mf_mx, hmx);
  copy_to_host(test.mf_my, hmy);
  const auto native = hn.const_array();
  const auto candidate = hc.const_array();
  const auto raw_control = hr.const_array();
  const auto cross_control = hw.const_array();
  const auto accepted = ha.const_array();
  const auto accx = host_x.const_array();
  const auto accy = host_y.const_array();
  const auto accz = host_z.const_array();
  const auto det = hdet.const_array();
  const auto mxh = hmx.const_array();
  const auto myh = hmy.const_array();
  bool differs_from_original = false;
  for (int k = test.bx.smallEnd(2); k <= test.bx.bigEnd(2); ++k) {
    for (int j = test.bx.smallEnd(1); j <= test.bx.bigEnd(1); ++j) {
      for (int i = test.bx.smallEnd(0); i <= test.bx.bigEnd(0); ++i) {
        const Real expected_native = native(i, j, k, RhoScalar_comp);
        const Real canonical_value = candidate(i, j, k, RhoScalar_comp);
        const Real accepted_value = accepted(i, j, k, rhs_comp);
        EXPECT_NEAR(expected_native, canonical_value,
                    integrated_terrain_tolerance(expected_native));
        EXPECT_GT(std::abs(expected_native - raw_control(i, j, k, rhs_comp)),
                  Real(100.0) * tolerance());
        EXPECT_GT(std::abs(expected_native - cross_control(i, j, k, rhs_comp)),
                  Real(100.0) * tolerance());
        differs_from_original |=
          std::abs(accepted_value - canonical_value) > Real(100.0) * tolerance();
        EXPECT_GT(std::abs(accepted_value - raw_control(i, j, k, rhs_comp)),
                  Real(100.0) * tolerance());
        EXPECT_GT(std::abs(accepted_value - cross_control(i, j, k, rhs_comp)),
                  Real(100.0) * tolerance());

        const Real direct = -((mxh(i, j, 0) * myh(i, j, 0)) / det(i, j, k)) *
          ((accx(i + 1, j, k, accepted_comp) - accx(i, j, k, accepted_comp)) *
             test.inv[0] +
           (accy(i, j + 1, k, accepted_comp) - accy(i, j, k, accepted_comp)) *
             test.inv[1] +
           (accz(i, j, k + 1, accepted_comp) - accz(i, j, k, accepted_comp)) *
             test.inv[2]);
        EXPECT_NEAR(accepted_value, direct, tolerance(direct));
      }
    }
  }
  EXPECT_TRUE(differs_from_original);
  EXPECT_NEAR(low_weight + lambda, Real(1.0), tolerance());
}

ERF_GPU_TEST(
  ScalarDiffusionPrimitives,
  TerrainMappedTransfersSupportInPlaceAndLowerSuppression)
{
  NativeTerrainScalarCase terrain(
    Real(0.31), Real(-0.23), Real(1.4), Real(1.7), Real(0.82), Real(0.61),
    RhoScalar_comp, true);
  const Box& domain = terrain.domain;
  const Box xfaces = surroundingNodes(domain, 0);
  const Box yfaces = surroundingNodes(domain, 1);
  const Box zfaces = surroundingNodes(domain, 2);
  FArrayBox raw_x(xfaces, 1), raw_y(yfaces, 1), raw_z(zfaces, 1);
  FArrayBox mapped_x(xfaces, 1), mapped_y(yfaces, 1), mapped_z(zfaces, 1);
  FArrayBox suppressed_x(xfaces, 1), suppressed_y(yfaces, 1),
    suppressed_z(zfaces, 1);
  FArrayBox inplace_x(xfaces, 1), inplace_y(yfaces, 1), inplace_z(zfaces, 1);

  auto fx = raw_x.array();
  ParallelFor(xfaces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    fx(i, j, k) = Real(1.3) + Real(0.07) * i + Real(0.03) * j + Real(0.11) * k +
                  Real(0.005) * i * k;
  });
  auto fy = raw_y.array();
  ParallelFor(yfaces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    fy(i, j, k) = -Real(1.0) + Real(0.02) * i - Real(0.03) * j - Real(0.05) * k;
  });
  auto fz = raw_z.array();
  ParallelFor(zfaces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    fz(i, j, k) = Real(0.4) + Real(0.06) * i - Real(0.025) * j + Real(0.02) * k;
  });

  const auto copy_x_src = raw_x.const_array();
  const auto copy_y_src = raw_y.const_array();
  const auto copy_z_src = raw_z.const_array();
  auto copy_x_dst = inplace_x.array();
  auto copy_y_dst = inplace_y.array();
  auto copy_z_dst = inplace_z.array();
  ParallelFor(xfaces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    copy_x_dst(i, j, k) = copy_x_src(i, j, k);
  });
  ParallelFor(yfaces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    copy_y_dst(i, j, k) = copy_y_src(i, j, k);
  });
  ParallelFor(zfaces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    copy_z_dst(i, j, k) = copy_z_src(i, j, k);
  });

  BuildScalarDiffusionMappedTransfers_T(
    domain, domain, raw_x.const_array(), raw_y.const_array(),
    raw_z.const_array(), 0, mapped_x.array(), 0, mapped_y.array(), 0,
    mapped_z.array(), 0, terrain.z_nd.const_array(), terrain.ax.const_array(),
    terrain.ay.const_array(), terrain.inv, terrain.mf_mx.const_array(),
    terrain.mf_uy.const_array(), terrain.mf_my.const_array(),
    terrain.mf_vx.const_array(), false);
  BuildScalarDiffusionMappedTransfers_T(
    domain, domain, raw_x.const_array(), raw_y.const_array(),
    raw_z.const_array(), 0, suppressed_x.array(), 0, suppressed_y.array(), 0,
    suppressed_z.array(), 0, terrain.z_nd.const_array(),
    terrain.ax.const_array(), terrain.ay.const_array(), terrain.inv,
    terrain.mf_mx.const_array(), terrain.mf_uy.const_array(),
    terrain.mf_my.const_array(), terrain.mf_vx.const_array(), true);

  // The z launch must read raw x/y before the in-place x/y conversions.
  BuildScalarDiffusionMappedTransfers_T(
    domain, domain, inplace_x.const_array(), inplace_y.const_array(),
    inplace_z.const_array(), 0, inplace_x.array(), 0, inplace_y.array(), 0,
    inplace_z.array(), 0, terrain.z_nd.const_array(), terrain.ax.const_array(),
    terrain.ay.const_array(), terrain.inv, terrain.mf_mx.const_array(),
    terrain.mf_uy.const_array(), terrain.mf_my.const_array(),
    terrain.mf_vx.const_array(), false);
  Gpu::streamSynchronize();

  FArrayBox hmapped_x(xfaces, 1, The_Pinned_Arena());
  FArrayBox hmapped_y(yfaces, 1, The_Pinned_Arena());
  FArrayBox hmapped_z(zfaces, 1, The_Pinned_Arena());
  FArrayBox hsuppressed_z(zfaces, 1, The_Pinned_Arena());
  FArrayBox hinplace_x(xfaces, 1, The_Pinned_Arena());
  FArrayBox hinplace_y(yfaces, 1, The_Pinned_Arena());
  FArrayBox hinplace_z(zfaces, 1, The_Pinned_Arena());
  copy_to_host(mapped_x, hmapped_x);
  copy_to_host(mapped_y, hmapped_y);
  copy_to_host(mapped_z, hmapped_z);
  copy_to_host(suppressed_z, hsuppressed_z);
  copy_to_host(inplace_x, hinplace_x);
  copy_to_host(inplace_y, hinplace_y);
  copy_to_host(inplace_z, hinplace_z);
  const auto mx_ref = hmapped_x.const_array();
  const auto my_ref = hmapped_y.const_array();
  const auto mz_ref = hmapped_z.const_array();
  const auto mz_suppressed = hsuppressed_z.const_array();
  const auto mx_inplace = hinplace_x.const_array();
  const auto my_inplace = hinplace_y.const_array();
  const auto mz_inplace = hinplace_z.const_array();
  for (int k = xfaces.smallEnd(2); k <= xfaces.bigEnd(2); ++k) {
    for (int j = xfaces.smallEnd(1); j <= xfaces.bigEnd(1); ++j) {
      for (int i = xfaces.smallEnd(0); i <= xfaces.bigEnd(0); ++i) {
        EXPECT_NEAR(
          mx_ref(i, j, k), mx_inplace(i, j, k), tolerance(mx_ref(i, j, k)));
      }
    }
  }
  for (int k = yfaces.smallEnd(2); k <= yfaces.bigEnd(2); ++k) {
    for (int j = yfaces.smallEnd(1); j <= yfaces.bigEnd(1); ++j) {
      for (int i = yfaces.smallEnd(0); i <= yfaces.bigEnd(0); ++i) {
        EXPECT_NEAR(
          my_ref(i, j, k), my_inplace(i, j, k), tolerance(my_ref(i, j, k)));
      }
    }
  }
  for (int k = zfaces.smallEnd(2); k <= zfaces.bigEnd(2); ++k) {
    for (int j = zfaces.smallEnd(1); j <= zfaces.bigEnd(1); ++j) {
      for (int i = zfaces.smallEnd(0); i <= zfaces.bigEnd(0); ++i) {
        EXPECT_NEAR(
          mz_ref(i, j, k), mz_inplace(i, j, k), tolerance(mz_ref(i, j, k)));
      }
    }
  }

  constexpr int sample_i = 2, sample_j = 3;
  const int klo = domain.smallEnd(2);
  const int khi = domain.bigEnd(2) + 1;
  const Real raw_lower_z = Real(0.4) + Real(0.06) * sample_i -
                           Real(0.025) * sample_j + Real(0.02) * klo;
  const Real raw_lower_transfer =
    ScalarDiffusionMappedVerticalTransfer(raw_lower_z, terrain.mx, terrain.my);
  EXPECT_DOUBLE_EQ(mz_suppressed(sample_i, sample_j, klo), Real(0.0));
  EXPECT_GT(
    std::abs(mz_ref(sample_i, sample_j, klo)),
    Real(100.0) * tolerance(mz_ref(sample_i, sample_j, klo)));
  EXPECT_GT(
    std::abs(mz_ref(sample_i, sample_j, klo) - raw_lower_transfer),
    Real(100.0) * tolerance(raw_lower_transfer));
  for (int k = klo + 1; k <= khi; ++k) {
    EXPECT_DOUBLE_EQ(
      mz_suppressed(sample_i, sample_j, k), mz_ref(sample_i, sample_j, k));
  }
}

ERF_GPU_TEST(
  ScalarDiffusionPrimitives, MappedDivergenceAccumulationAndInvalidJacobian)
{
  const Box cells(IntVect(0, 0, 0), IntVect(2, 0, 0));
  const Box xfaces = surroundingNodes(cells, 0);
  const Box yfaces = surroundingNodes(cells, 1);
  const Box zfaces = surroundingNodes(cells, 2);
  FArrayBox x(xfaces, 1), y(yfaces, 1), z(zfaces, 1);
  FArrayBox rhs(cells, 1), detj(cells, 1), mx(cells, 1), my(cells, 1);
  const auto fx = x.array();
  const auto fy = y.array();
  const auto fz = z.array();
  ParallelFor(xfaces, [=] AMREX_GPU_DEVICE(int i, int, int) noexcept {
    fx(i, 0, 0) = Real(2.0) * i;
  });
  ParallelFor(yfaces, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
    fy(i, j, 0) = Real(3.0) * j;
  });
  ParallelFor(zfaces, [=] AMREX_GPU_DEVICE(int i, int, int k) noexcept {
    fz(i, 0, k) = Real(5.0) * k;
  });
  const auto d = detj.array();
  const auto mx4 = mx.array();
  const auto my4 = my.array();
  ParallelFor(cells, [=] AMREX_GPU_DEVICE(int i, int, int) noexcept {
    d(i, 0, 0) = i == 0 ? Real(4.0) : (i == 1 ? Real(0.0) : Real(-1.0));
    mx4(i, 0, 0) = Real(2.0);
    my4(i, 0, 0) = Real(3.0);
  });
  rhs.setVal<RunOn::Device>(Real(7.0));
  const GpuArray<Real, AMREX_SPACEDIM> inv{{Real(0.5), Real(0.25), Real(0.1)}};
  Gpu::streamSynchronize();

  EXPECT_DOUBLE_EQ(
    ScalarDiffusionMappedFluxDivergence(
      Real(2.0), Real(0.0), Real(3.0), Real(0.0), Real(5.0), Real(0.0), inv[0],
      inv[1], inv[2], Real(2.0), Real(3.0), Real(0.0)),
    Real(0.0));
  EXPECT_DOUBLE_EQ(
    ScalarDiffusionMappedFluxDivergence(
      Real(2.0), Real(0.0), Real(3.0), Real(0.0), Real(5.0), Real(0.0), inv[0],
      inv[1], inv[2], Real(2.0), Real(3.0), Real(-1.0)),
    Real(0.0));
  EXPECT_DOUBLE_EQ(
    ScalarDiffusionMappedFluxDivergence(
      Real(2.0), Real(0.0), Real(3.0), Real(0.0), Real(5.0), Real(0.0), inv[0],
      inv[1], inv[2], Real(2.0), Real(3.0), Real(4.0)),
    Real(3.375));

  ApplyScalarMappedFluxDivergence(
    cells, x.const_array(), 0, y.const_array(), 0, z.const_array(), 0,
    rhs.array(), 0, detj.const_array(), inv, mx.const_array(),
    my.const_array());
  Gpu::streamSynchronize();

  FArrayBox host_rhs(cells, 1, The_Pinned_Arena());
  copy_to_host(rhs, host_rhs);
  const auto result = host_rhs.const_array();
  EXPECT_DOUBLE_EQ(result(0, 0, 0), Real(3.625));
  EXPECT_DOUBLE_EQ(result(1, 0, 0), Real(7.0));
  EXPECT_DOUBLE_EQ(result(2, 0, 0), Real(7.0));
}

ERF_GPU_TEST(ScalarDiffusionPrimitives, CanonicalMappedDivergenceTelescopes)
{
  const Box cells(IntVect(0, 0, 0), IntVect(2, 2, 1));
  FArrayBox x(surroundingNodes(cells, 0), 2),
    y(surroundingNodes(cells, 1), 2), z(surroundingNodes(cells, 2), 2);
  FArrayBox rhs(cells, 1), advection_rhs(cells, 1), detj(cells, 1), mx(cells, 1),
    my(cells, 1);
  x.setVal<RunOn::Device>(Real(-1.0));
  y.setVal<RunOn::Device>(Real(-2.0));
  z.setVal<RunOn::Device>(Real(-3.0));
  // Zero baselines make this compare only the shared face/geometry
  // convention; diffusion accumulates while advection assigns its RHS.
  rhs.setVal<RunOn::Device>(Real(0.0));
  advection_rhs.setVal<RunOn::Device>(Real(0.0));
  auto fx = x.array();
  auto fy = y.array();
  auto fz = z.array();
  ParallelFor(x.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const int phase = i % 3;
    const Real periodic = phase == 0 ? Real(0.0)
                          : phase == 1 ? Real(0.86602540378443864676)
                                       : Real(-0.86602540378443864676);
    fx(i, j, k, 1) = periodic +
                     Real(0.13) * j + Real(0.07) * (k % 2);
  });
  ParallelFor(y.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const int phase = j % 3;
    const Real periodic = phase == 0 ? Real(1.0) : Real(-0.5);
    fy(i, j, k, 1) = periodic +
                     Real(0.11) * i - Real(0.05) * (k % 2);
  });
  ParallelFor(z.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    fz(i, j, k, 1) = Real(0.17) * (k % 2) +
                     Real(0.09) * i - Real(0.04) * j;
  });
  auto det = detj.array();
  auto mx4 = mx.array();
  auto my4 = my.array();
  ParallelFor(cells, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    det(i, j, k) = Real(1.1) + Real(0.03) * i + Real(0.02) * j +
                   Real(0.015) * k;
    mx4(i, j, k) = Real(1.2) + Real(0.02) * i + Real(0.01) * j;
    my4(i, j, k) = Real(0.8) + Real(0.03) * j + Real(0.01) * i;
  });
  const GpuArray<Real, AMREX_SPACEDIM> inv{{Real(0.7), Real(0.9), Real(1.1)}};
  const GpuArray<const Array4<Real>, AMREX_SPACEDIM> flux_views{
    {x.array(), y.array(), z.array()}};
  ApplyScalarMappedFluxDivergence(
    cells, x.const_array(), 1, y.const_array(), 1, z.const_array(), 1,
    rhs.array(), 0, detj.const_array(), inv, mx.const_array(), my.const_array());
  ApplyScalarAdvectionFluxDivergence(
    cells, flux_views, 1, advection_rhs.array(), 0, detj.const_array(), inv,
    mx.const_array(), my.const_array());
  Gpu::streamSynchronize();

  FArrayBox hrhs(cells, 1, The_Pinned_Arena()), hdet(cells, 1, The_Pinned_Arena()),
    hadvection(cells, 1, The_Pinned_Arena()), hmx(cells, 1, The_Pinned_Arena()),
    hmy(cells, 1, The_Pinned_Arena());
  copy_to_host(rhs, hrhs);
  copy_to_host(advection_rhs, hadvection);
  copy_to_host(detj, hdet);
  copy_to_host(mx, hmx);
  copy_to_host(my, hmy);
  const auto l = hrhs.const_array();
  const auto adv = hadvection.const_array();
  const auto j = hdet.const_array();
  const auto a = hmx.const_array();
  const auto b = hmy.const_array();
  Real weighted_sum = Real(0.0), absolute_sum = Real(0.0);
  for (int k = cells.smallEnd(2); k <= cells.bigEnd(2); ++k) {
    for (int jj = cells.smallEnd(1); jj <= cells.bigEnd(1); ++jj) {
      for (int i = cells.smallEnd(0); i <= cells.bigEnd(0); ++i) {
      const Real weighted =
        (j(i, jj, k) / (a(i, jj, k) * b(i, jj, k))) * l(i, jj, k, 0);
      EXPECT_NEAR(l(i, jj, k, 0), adv(i, jj, k, 0),
                  tolerance(l(i, jj, k, 0)));
        weighted_sum += weighted;
        absolute_sum += std::abs(weighted);
      }
    }
  }
  EXPECT_NEAR(weighted_sum, Real(0.0),
              Real(256.0) * std::numeric_limits<Real>::epsilon() *
                std::max(Real(1.0), absolute_sum));
}

#undef ERF_GPU_TEST

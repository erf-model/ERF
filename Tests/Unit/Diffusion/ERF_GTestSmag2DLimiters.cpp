// Unit tests for the opt-in WRF Smagorinsky-2D limits (ERF_TerrainDiffusionLimits.H) and the
// explicit diffusion rates of the diffusive time-step check:
//  * the slope factor alpha against hand-computed terrain drops, including the absolute values
//    WRF takes and the clamp at 1;
//  * the limited viscosity against hand-computed WRF smag2d_km numbers on both sides of the
//    deformation threshold max(10/Delta_h, 1e-3), with the cap active and inactive and in WRF's
//    order (cap, then slope);
//  * the momentum and scalar diffusive rates against hand-computed values;
//  * ComputeTurbulentViscosityLES on a tilted terrain-fitted mesh: K_h is limited before the
//    scalar diffusivities are derived from it, K_v is not touched, and the limiter does nothing
//    without terrain.

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <ERF_EddyViscosity.H>
#include <ERF_TerrainDiffusionLimits.H>

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <memory>
#include <vector>

using amrex::Real;

namespace {

Real rel_tol ()
{
    return (sizeof(Real) == 8) ? Real(1e-12) : Real(2e-5);
}

// Host-only nodal heights of one cell, z(i,j,k) for i,j,k in {0,1}, wrapped in an Array4.
struct HostCell
{
    std::vector<Real> z = std::vector<Real>(8, Real(0));
    Real& at (int i, int j, int k) { return z[static_cast<std::size_t>(i + 2*j + 4*k)]; }
    amrex::Array4<const Real> array () const
    {
        return amrex::Array4<const Real>(z.data(), amrex::Dim3{0,0,0}, amrex::Dim3{2,2,2}, 1);
    }
};

// A plane z = sx x + sy y + dz k on a dx by dy cell
HostCell plane_cell (Real sx, Real sy, Real dx, Real dy, Real dz)
{
    HostCell c;
    for (int k = 0; k < 2; ++k) {
        for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 2; ++i) {
                c.at(i,j,k) = sx*dx*Real(i) + sy*dy*Real(j) + dz*Real(k);
            }
        }
    }
    return c;
}

} // namespace

TEST(Smag2DLimiters, SlopeFactorOfATiltedPlane)
{
    // 3 km cells 50 m thick: slope 0.1 in x and 0.05 in y drop 300 m and 150 m across a cell,
    // alpha = sqrt(300^2 + 150^2) / 50 = 6.7082.
    const HostCell c = plane_cell(Real(0.1), Real(0.05), Real(3000), Real(3000), Real(50));
    const TerrainCellDrops d = ComputeTerrainCellDrops(0,0,0,c.array());
    EXPECT_NEAR(d.dx_drop, Real(300), Real(300)*rel_tol());
    EXPECT_NEAR(d.dy_drop, Real(150), Real(150)*rel_tol());
    EXPECT_NEAR(d.dz, Real(50), Real(50)*rel_tol());
    const Real alpha = TerrainSlopeFactor(d);
    EXPECT_NEAR(alpha, std::sqrt(Real(300*300 + 150*150)) / Real(50), Real(7)*rel_tol());
    EXPECT_GT(alpha, Real(6));

    // Downhill in x is the same drop
    const HostCell c2 = plane_cell(Real(-0.1), Real(0.05), Real(3000), Real(3000), Real(50));
    EXPECT_NEAR(TerrainSlopeFactor(ComputeTerrainCellDrops(0,0,0,c2.array())), alpha, Real(7)*rel_tol());
}

TEST(Smag2DLimiters, SlopeFactorClampsAtOne)
{
    // Slope 0.01 over 3 km is a 30 m drop in a 50 m cell: 0.6, clamped to 1 as in WRF.
    const HostCell c = plane_cell(Real(0.01), Real(0), Real(3000), Real(3000), Real(50));
    EXPECT_EQ(TerrainSlopeFactor(ComputeTerrainCellDrops(0,0,0,c.array())), Real(1));
    // Flat ground
    const HostCell f = plane_cell(Real(0), Real(0), Real(3000), Real(3000), Real(50));
    EXPECT_EQ(TerrainSlopeFactor(ComputeTerrainCellDrops(0,0,0,f.array())), Real(1));
}

TEST(Smag2DLimiters, SlopeFactorAveragesAbsoluteDrops)
{
    // A saddle: the bottom nodes are 0, 60, 60, 0 m at (i,j) = (0,0), (1,0), (0,1), (1,1), and
    // the cell is 40 m thick everywhere.  Along x the drop is +60 m on the j = 0 edges and
    // -60 m on the j = 1 edges (likewise along y), so the signed means are 0 and alpha would
    // clamp to 1.  WRF averages |z_x|: both drops are 60 m and alpha = sqrt(2) * 60 / 40.
    HostCell c;
    for (int k = 0; k < 2; ++k) {
        const Real top = Real(40) * Real(k);
        c.at(0,0,k) = top;  c.at(1,0,k) = Real(60) + top;
        c.at(0,1,k) = Real(60) + top;  c.at(1,1,k) = top;
    }
    const TerrainCellDrops d = ComputeTerrainCellDrops(0,0,0,c.array());
    EXPECT_NEAR(d.dx_drop, Real(60), Real(60)*rel_tol());
    EXPECT_NEAR(d.dy_drop, Real(60), Real(60)*rel_tol());
    EXPECT_NEAR(d.dz, Real(40), Real(40)*rel_tol());
    EXPECT_NEAR(TerrainSlopeFactor(d), std::sqrt(Real(2)) * Real(1.5), Real(3)*rel_tol());
}

TEST(Smag2DLimiters, LimitedViscosityBothDeformationBranches)
{
    // Delta_h = 3 km: def_limit = max(10/3000, 1e-3) = 3.333e-3 1/s.  Cs = 0.1, alpha = 6.
    const Real DeltaH = Real(3000);
    const Real alpha  = Real(6);
    const Real Cs2D2  = Real(0.01) * DeltaH * DeltaH;

    // |D| = 5e-3 > def_limit: nu = 450 m^2/s, divided by alpha^2 -> 12.5
    Real strain = Real(5e-3);
    Real nu = Smag2DLimitedViscosity(Cs2D2*strain, strain, DeltaH, alpha, Real(0), true);
    EXPECT_NEAR(nu, Real(12.5), Real(12.5)*rel_tol());

    // |D| = 2e-3 < def_limit: nu = 180 m^2/s, divided by alpha -> 30
    strain = Real(2e-3);
    nu = Smag2DLimitedViscosity(Cs2D2*strain, strain, DeltaH, alpha, Real(0), true);
    EXPECT_NEAR(nu, Real(30), Real(30)*rel_tol());

    // |D| exactly at def_limit: WRF tests tmp .gt. def_limit, so it divides by alpha
    strain = Real(10) / DeltaH;
    nu = Smag2DLimitedViscosity(Real(100), strain, DeltaH, alpha, Real(0), true);
    EXPECT_NEAR(nu, Real(100) / alpha, Real(100)*rel_tol());
}

TEST(Smag2DLimiters, LimitedViscosityDeformationFloor)
{
    // Delta_h = 20 km: 10/Delta_h = 5e-4 is below the floor, so def_limit = 1e-3.
    const Real DeltaH = Real(20000);
    const Real alpha  = Real(3);
    // 8e-4 is above 10/Delta_h but below the floor: divided by alpha
    EXPECT_NEAR(Smag2DLimitedViscosity(Real(90), Real(8e-4), DeltaH, alpha, Real(0), true),
                Real(30), Real(30)*rel_tol());
    // 1.2e-3 is above the floor: divided by alpha^2
    EXPECT_NEAR(Smag2DLimitedViscosity(Real(90), Real(1.2e-3), DeltaH, alpha, Real(0), true),
                Real(10), Real(10)*rel_tol());
}

TEST(Smag2DLimiters, LimitedViscosityAlphaOneIsUnchanged)
{
    for (const Real strain : {Real(1e-4), Real(1e-2)}) {
        EXPECT_EQ(Smag2DLimitedViscosity(Real(123.5), strain, Real(3000), Real(1), Real(0), true),
                  Real(123.5));
    }
    // Both limits off returns the input bit for bit
    EXPECT_EQ(Smag2DLimitedViscosity(Real(123.5), Real(1), Real(3000), Real(6), Real(0), false),
              Real(123.5));
}

TEST(Smag2DLimiters, CapActiveInactiveAndOrder)
{
    const Real DeltaH = Real(3000);
    // Cap 10 m/s * 3 km = 3e4 m^2/s.
    // Inactive: nu = 450 stays, then / alpha^2 (|D| = 5e-3 > 3.33e-3) -> 12.5
    EXPECT_NEAR(Smag2DLimitedViscosity(Real(450), Real(5e-3), DeltaH, Real(6), Real(10), true),
                Real(12.5), Real(12.5)*rel_tol());
    // Active: nu = 2.8125e6 (Cs = 0.25, |D| = 5 1/s) capped to 3e4, then / 36 -> 833.33.
    // The reverse order (slope first) would give min(2.8125e6/36, 3e4) = 3e4.
    EXPECT_NEAR(Smag2DLimitedViscosity(Real(2.8125e6), Real(5), DeltaH, Real(6), Real(10), true),
                Real(3e4) / Real(36), Real(3e4/36)*rel_tol());
    // Cap alone
    EXPECT_NEAR(Smag2DLimitedViscosity(Real(2.8125e6), Real(5), DeltaH, Real(6), Real(10), false),
                Real(3e4), Real(3e4)*rel_tol());
    EXPECT_EQ(Smag2DLimitedViscosity(Real(450), Real(5), DeltaH, Real(6), Real(10), false), Real(450));
}

TEST(Smag2DLimiters, DiffusiveRatesHandValues)
{
    // rho = 1, dx = dy = 1 km (a = b = 1e-6), slope 0.1 (h^2 = 0.01), dz = 50 m (1/dz^2 = 4e-4)
    const Real a = Real(1e-6), b = Real(1e-6), h2 = Real(0.01), dzinv2 = Real(4e-4);
    // Momentum, mu_h = 100, mu_v = 10, fully explicit vertical (tau33 = 2 mu S33 for w):
    //   100 * 3e-6 + 2 * 100 * 0.01 * 4e-4 + 2 * 10 * 4e-4 = 3e-4 + 8e-4 + 8e-3 = 9.1e-3
    EXPECT_NEAR(MomentumDiffusiveRate(Real(1), Real(100), Real(10), a, b, h2, dzinv2, Real(1), Real(1)),
                Real(9.1e-3), Real(9.1e-3)*rel_tol());
    // u, v vertical implicit but w explicit (e_uv = 0, e_w = 1): w's 2 mu_v / dz^2 stays
    EXPECT_NEAR(MomentumDiffusiveRate(Real(1), Real(100), Real(10), a, b, h2, dzinv2, Real(0), Real(1)),
                Real(9.1e-3), Real(9.1e-3)*rel_tol());
    // w implicit too but u explicit (e_uv = 1, e_w = 0): mu_v / dz^2 = 4e-3 -> 5.1e-3
    EXPECT_NEAR(MomentumDiffusiveRate(Real(1), Real(100), Real(10), a, b, h2, dzinv2, Real(1), Real(0)),
                Real(5.1e-3), Real(5.1e-3)*rel_tol());
    // All vertical implicit (e = 0) leaves 1.1e-3; rho = 2 halves it
    EXPECT_NEAR(MomentumDiffusiveRate(Real(2), Real(100), Real(10), a, b, h2, dzinv2, Real(0), Real(0)),
                Real(5.5e-4), Real(5.5e-4)*rel_tol());
    // mu_v > mu_h: the projected term uses the larger, 2 * 40 * 0.01 * 4e-4 = 3.2e-4
    EXPECT_NEAR(MomentumDiffusiveRate(Real(1), Real(10), Real(40), a, b, h2, dzinv2, Real(0), Real(0)),
                Real(3e-5) + Real(3.2e-4), Real(3.5e-4)*rel_tol());
    // Unequal spacings: max(2a + b, a + 2b) with a = 4e-6, b = 1e-6 -> 9e-6
    EXPECT_NEAR(MomentumDiffusiveRate(Real(1), Real(100), Real(0), Real(4e-6), Real(1e-6), Real(0), dzinv2, Real(1), Real(1)),
                Real(9e-4), Real(9e-4)*rel_tol());
    // Scalar, K_h = 300, K_v = 30: 300 * 2e-6 + 300 * 0.01 * 4e-4 + 30 * 4e-4 = 1.38e-2
    EXPECT_NEAR(ScalarDiffusiveRate(Real(1), Real(300), Real(30), a, b, h2, dzinv2, Real(1)),
                Real(1.38e-2), Real(1.38e-2)*rel_tol());
    EXPECT_NEAR(ScalarDiffusiveRate(Real(1), Real(300), Real(30), a, b, h2, dzinv2, Real(0.5)),
                Real(7.8e-3), Real(7.8e-3)*rel_tol());
}

TEST(Smag2DLimiters, ExplicitVerticalFractionFromStability)
{
    using V = std::vector<Real>;
    const V rk3 = {Real(1)/Real(3), Real(0.5), Real(1)};
    const V mid = {Real(0.5), Real(1)};
    // ERF's default 1 1 0: stage 2 implicit, stage 3 explicit on its state -> Crank-Nicolson,
    // A-stable, so no explicit vertical limit
    EXPECT_EQ(ExplicitVerticalFraction(rk3, V{1, 1, 0}, true), Real(0));
    EXPECT_EQ(ExplicitVerticalFraction(rk3, V{1, 1, 1}, true), Real(0));
    EXPECT_EQ(ExplicitVerticalFraction(rk3, V{1, 1, 0.5}, true), Real(0));
    // an explicit first stage amplifies stiff modes before stage 2 sees them (S1 = 1 + z/3)
    EXPECT_EQ(ExplicitVerticalFraction(rk3, V{0, 1, 0.5}, true), Real(1));
    // fully explicit, and partly implicit patterns that are not A-stable
    EXPECT_EQ(ExplicitVerticalFraction(rk3, V{0, 0, 0}, true), Real(1));
    EXPECT_EQ(ExplicitVerticalFraction(rk3, V{1, 0, 0}, true), Real(1));
    EXPECT_EQ(ExplicitVerticalFraction(rk3, V{1, 0, 0.25}, true), Real(1));
    EXPECT_EQ(ExplicitVerticalFraction(rk3, V{0.5, 0.5, 0.25}, true), Real(1));
    // only the last stage implicit damps at the end, but stages 1 and 2 amplify (S2 = 58 at
    // z = -20), and those states feed the explicit terms of the next stage
    EXPECT_EQ(ExplicitVerticalFraction(rk3, V{0, 0, 1}, true), Real(1));
    // the component left out of the implicit solve is explicit whatever the factors
    EXPECT_EQ(ExplicitVerticalFraction(rk3, V{1, 1, 1}, false), Real(1));
    // anelastic MidPoint (factors forced to 1 0): implicit midpoint, A-stable
    EXPECT_EQ(ExplicitVerticalFraction(mid, V{1, 0, 0}, true), Real(0));
    EXPECT_EQ(ExplicitVerticalFraction(mid, V{0, 0, 0}, true), Real(1));
}

// ---------------------------------------------------------------------------------------------
// ComputeTurbulentViscosityLES on a tilted plane, z = sx x + sy y + dz k, with uniform strain.
// ---------------------------------------------------------------------------------------------
namespace {

constexpr Real t_dx  = Real(3000);
constexpr Real t_dz  = Real(50);
constexpr Real t_sx  = Real(0.1);
constexpr Real t_sy  = Real(0.05);
constexpr Real t_rho = Real(1.2);
constexpr Real t_Cs  = Real(0.1);

struct LESResult { Real mom_h_min, mom_h_max, theta_h_min, theta_h_max, mom_v_min, mom_v_max; };

void fill_tilted_heights (amrex::MultiFab& z_nd)
{
    for (amrex::MFIter mfi(z_nd); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.fabbox();
        const auto z = z_nd.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            z(i,j,k) = t_sx*t_dx*Real(i) + t_sy*t_dx*Real(j) + t_dz*Real(k);
        });
    }
}

LESResult minmax (const amrex::MultiFab& K)
{
    amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMax, amrex::ReduceOpMin,
                     amrex::ReduceOpMax, amrex::ReduceOpMin, amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<Real, Real, Real, Real, Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (amrex::MFIter mfi(K); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        const auto k_arr = K.const_array(mfi);
        reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            const Real m = k_arr(i,j,k,EddyDiff::Mom_h);
            const Real t = k_arr(i,j,k,EddyDiff::Theta_h);
            const Real v = k_arr(i,j,k,EddyDiff::Mom_v);
            return {m, m, t, t, v, v};
        });
    }
    auto hv = reduce_data.value(reduce_op);
    LESResult r{amrex::get<0>(hv), amrex::get<1>(hv), amrex::get<2>(hv),
                amrex::get<3>(hv), amrex::get<4>(hv), amrex::get<5>(hv)};
    amrex::ParallelDescriptor::ReduceRealMin(r.mom_h_min);
    amrex::ParallelDescriptor::ReduceRealMax(r.mom_h_max);
    amrex::ParallelDescriptor::ReduceRealMin(r.theta_h_min);
    amrex::ParallelDescriptor::ReduceRealMax(r.theta_h_max);
    amrex::ParallelDescriptor::ReduceRealMin(r.mom_v_min);
    amrex::ParallelDescriptor::ReduceRealMax(r.mom_v_max);
    return r;
}

// S11 = ux, S22 = vy, S12 uniform: |D_h| = sqrt((ux - vy)^2 + 4 S12^2)
LESResult run_les (Real strain_scale, bool slope, Real cap, bool terrain, Real mf = Real(1))
{
    const int n = 8;
    amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n-1, n-1, 3));
    amrex::RealBox rb({Real(0), Real(0), Real(0)}, {n*t_dx, n*t_dx, 4*t_dz});
    amrex::Array<int,AMREX_SPACEDIM> periodic{1, 1, 0};
    amrex::Geometry geom(domain, rb, amrex::CoordSys::cartesian, periodic);
    amrex::BoxArray ba(domain);
    ba.maxSize(4);
    amrex::DistributionMapping dm(ba);

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> tau(9);
    const Real ux = Real(2e-3)*strain_scale, vy = Real(-1e-3)*strain_scale, s12 = Real(1e-3)*strain_scale;
    for (int c = 0; c < 9; ++c) {
        tau[c] = std::make_unique<amrex::MultiFab>(ba, dm, 1, 2);
        tau[c]->setVal(Real(0));
    }
    tau[TauType::tau11]->setVal(ux);
    tau[TauType::tau22]->setVal(vy);
    tau[TauType::tau12]->setVal(s12);

    amrex::MultiFab cons(ba, dm, RhoTheta_comp+1, 1);
    cons.setVal(t_rho, Rho_comp, 1, 1);
    cons.setVal(t_rho*Real(300), RhoTheta_comp, 1, 1);
    amrex::MultiFab K(ba, dm, EddyDiff::NumDiffs, 1);
    K.setVal(Real(0));
    amrex::MultiFab diss(ba, dm, 1, 0);

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> mapfac(MapFacType::num);
    for (auto& m : mapfac) {
        m = std::make_unique<amrex::MultiFab>(ba, dm, 1, 2);
        m->setVal(mf);
    }
    auto z_nd = std::make_unique<amrex::MultiFab>(amrex::convert(ba, amrex::IntVect(1)), dm, 1, 1);
    fill_tilted_heights(*z_nd);

    TurbChoice tc;
    tc.les_type = LESType::Smagorinsky;
    tc.smag2d = true;
    tc.mix_isotropic = false;
    tc.Cs = t_Cs;
    tc.smag2d_slope_limiter = slope;
    tc.smag2d_kh_cap = cap;

    std::unique_ptr<SurfaceLayer> no_sl;
    MoistureComponentIndices mi;
    ComputeTurbulentViscosityLES(tau, cons, K, diss, geom, terrain, mapfac, z_nd, tc,
                                 Real(9.81), no_sl, mi, nullptr, nullptr);
    return minmax(K);
}

Real unlimited_nu (Real strain_scale, Real mf = Real(1))
{
    const Real ux = Real(2e-3)*strain_scale, vy = Real(-1e-3)*strain_scale, s12 = Real(1e-3)*strain_scale;
    const Real D = std::sqrt((ux - vy)*(ux - vy) + Real(4)*s12*s12);
    return t_Cs*t_Cs*(t_dx/mf)*(t_dx/mf)*D;
}

} // namespace

TEST(Smag2DLimiters, LESKernelLimitsKhBeforeScalarDiffusivities)
{
    // alpha = sqrt(300^2 + 150^2) / 50 = 6.708 everywhere on the tilted plane
    const Real alpha2 = (Real(300*300) + Real(150*150)) / (t_dz*t_dz);
    const Real alpha  = std::sqrt(alpha2);
    const Real tol    = Real(10)*rel_tol();

    // Unlimited: |D| = 3.606e-3, nu = 324.5 m^2/s
    const LESResult off = run_les(Real(1), false, Real(0), true);
    const Real nu1 = unlimited_nu(Real(1));
    EXPECT_NEAR(off.mom_h_max, t_rho*nu1, t_rho*nu1*tol);
    EXPECT_NEAR(off.mom_h_min, t_rho*nu1, t_rho*nu1*tol);

    // |D| = 3.606e-3 > 10/3000: K_h / alpha^2, Theta_h = Pr_t_inv * K_h, K_v unchanged
    const LESResult big = run_les(Real(1), true, Real(0), true);
    EXPECT_NEAR(big.mom_h_max, t_rho*nu1/alpha2, t_rho*nu1/alpha2*tol);
    EXPECT_NEAR(big.mom_h_min, t_rho*nu1/alpha2, t_rho*nu1/alpha2*tol);
    EXPECT_NEAR(big.theta_h_max, Real(3)*t_rho*nu1/alpha2, Real(3)*t_rho*nu1/alpha2*tol);
    EXPECT_NEAR(big.theta_h_min, Real(3)*t_rho*nu1/alpha2, Real(3)*t_rho*nu1/alpha2*tol);
    EXPECT_EQ(big.mom_v_max, off.mom_v_max);
    EXPECT_EQ(big.mom_v_min, off.mom_v_min);
    EXPECT_GT(off.mom_v_max, Real(0));

    // |D| halved (1.8e-3 < 3.33e-3): K_h / alpha
    const LESResult small = run_les(Real(0.5), true, Real(0), true);
    const Real nuh = unlimited_nu(Real(0.5));
    EXPECT_NEAR(small.mom_h_max, t_rho*nuh/alpha, t_rho*nuh/alpha*tol);
    EXPECT_NEAR(small.theta_h_min, Real(3)*t_rho*nuh/alpha, Real(3)*t_rho*nuh/alpha*tol);

    // Without terrain-fitted coordinates the slope limiter is off
    const LESResult flat = run_les(Real(1), true, Real(0), false);
    EXPECT_NEAR(flat.mom_h_max, t_rho*nu1, t_rho*nu1*tol);
}

TEST(Smag2DLimiters, LESKernelMapFactorMovesTheThreshold)
{
    // Map factor 1.25: Delta_h = 3000 / 1.25 = 2400 m, def_limit = 10 / 2400 = 4.17e-3 > |D| =
    // 3.606e-3, so the same strain that took the alpha^2 branch at m = 1 now takes alpha.
    // alpha itself is a drop in cell thicknesses and does not change.
    const Real alpha = std::sqrt(Real(300*300) + Real(150*150)) / t_dz;
    const Real tol   = Real(10)*rel_tol();
    const Real mf    = Real(1.25);
    const LESResult r = run_les(Real(1), true, Real(0), true, mf);
    const Real nu = unlimited_nu(Real(1), mf);
    EXPECT_NEAR(r.mom_h_max, t_rho*nu/alpha, t_rho*nu/alpha*tol);
    EXPECT_NEAR(r.mom_h_min, t_rho*nu/alpha, t_rho*nu/alpha*tol);
}

TEST(Smag2DLimiters, LESKernelCap)
{
    const Real tol = Real(10)*rel_tol();
    // cap 0.05 m/s * 3 km = 150 m^2/s < 324.5: K_h = rho * 150, Theta_h = 3 rho * 150
    const LESResult capped = run_les(Real(1), false, Real(0.05), true);
    EXPECT_NEAR(capped.mom_h_max, t_rho*Real(150), t_rho*Real(150)*tol);
    EXPECT_NEAR(capped.theta_h_min, Real(3)*t_rho*Real(150), Real(3)*t_rho*Real(150)*tol);
    // cap 10 m/s is inactive here
    const LESResult loose = run_les(Real(1), false, Real(10), true);
    EXPECT_NEAR(loose.mom_h_max, t_rho*unlimited_nu(Real(1)), t_rho*unlimited_nu(Real(1))*tol);
}

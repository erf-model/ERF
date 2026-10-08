// Unit tests for the opt-in WRF Smagorinsky-2D limits (ERF_TerrainDiffusionLimits.H) and the
// explicit diffusion rates of the diffusive time-step check:
//  * the slope factor alpha against hand-computed terrain drops, including the absolute values
//    WRF takes and the clamp at 1;
//  * the limited viscosity against hand-computed WRF smag2d_km numbers on both sides of the
//    deformation threshold max(10/Delta_h, 1e-3), with the cap active and inactive and in WRF's
//    order (cap, then slope);
//  * the momentum and scalar diffusive rates against hand-computed values, the explicit
//    fraction of the partly implicit stages (every stage bounded or not), and the neighbourhood bound (NeighbourhoodMax,
//    NeighbourhoodMinPositive, CellDiffusiveRates) on coefficients that vary between cells;
//  * ComputeTurbulentViscosityLES on a tilted terrain-fitted mesh: K_h is limited before the
//    scalar diffusivities are derived from it, K_v is not touched, and the limiter does nothing
//    without terrain.

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <ERF_Diffusion.H>
#include <ERF_EddyViscosity.H>
#include <ERF_TerrainDiffusionLimits.H>

#include <gtest/gtest.h>

#include <algorithm>
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
    // rho = 1, dx = dy = 1 km (a = b = 1e-6), slope 0.1 (h^2 = 0.01), dz = 50 m (1/dz^2 = 4e-4).
    // Momentum rows: u = mu_h (2a + b) + e_uv mu_v/dz^2 + mu_v sqrt(a c), v likewise, w =
    // mu_v (a + b) + 2 e_w mu_v/dz^2 + mu_v (sqrt(a c) + sqrt(b c)), plus the metric term
    // 2 max(mu_h, mu_v) h^2/dz^2.
    const Real a = Real(1e-6), b = Real(1e-6), h2 = Real(0.01), dzinv2 = Real(4e-4);
    // The cross derivatives sqrt(a c) = 2e-5 are always counted (mu_v times 2e-5 per direction).
    // mu_h = 100, mu_v = 10, all explicit: u = 3e-4 + 4e-3 + 2e-4, w = 2e-5 + 8e-3 + 4e-4
    // (largest), metric 2 * 100 * 0.01 * 4e-4 = 8e-4 -> 9.22e-3
    EXPECT_NEAR(MomentumDiffusiveRate(Real(1), Real(100), Real(10), a, b, h2, dzinv2, Real(1), Real(1)),
                Real(9.22e-3), Real(9.22e-3)*rel_tol());
    // u, v vertical implicit, w explicit (e_uv = 0, e_w = 1): the w row still sets it
    EXPECT_NEAR(MomentumDiffusiveRate(Real(1), Real(100), Real(10), a, b, h2, dzinv2, Real(0), Real(1)),
                Real(9.22e-3), Real(9.22e-3)*rel_tol());
    // w implicit, u explicit (e_uv = 1, e_w = 0): u = 4.5e-3 > w = 2e-5 + 4e-4 -> 5.3e-3
    EXPECT_NEAR(MomentumDiffusiveRate(Real(1), Real(100), Real(10), a, b, h2, dzinv2, Real(1), Real(0)),
                Real(5.3e-3), Real(5.3e-3)*rel_tol());
    // A thick cell with the vertical part explicit (review case: PBL only, dx = 100 m, dz = 500 m):
    // the cross derivative mu_v / (dx dz) = 2e-5 mu_v is five times the vertical term
    // mu_v / dz^2 = 4e-6 mu_v.  u = 10 (4e-6 + 2e-5) = 2.4e-4; w = 10 (2e-4 + 8e-6 + 4e-5) =
    // 2.48e-3 sets the rate (2.08e-3 if the cross derivatives were dropped)
    EXPECT_NEAR(MomentumDiffusiveRate(Real(1), Real(0), Real(10), Real(1e-4), Real(1e-4), Real(0),
                                      Real(4e-6), Real(1), Real(1)),
                Real(2.48e-3), Real(2.48e-3)*rel_tol());
    // All vertical implicit: the explicit cross derivatives count, sqrt(a c) = 2e-5:
    // u = 3e-4 + 10 * 2e-5 = 5e-4, w = 2e-5 + 10 * 4e-5 = 4.2e-4; 5e-4 + metric 8e-4 = 1.3e-3;
    // rho = 2 halves it
    EXPECT_NEAR(MomentumDiffusiveRate(Real(2), Real(100), Real(10), a, b, h2, dzinv2, Real(0), Real(0)),
                Real(6.5e-4), Real(6.5e-4)*rel_tol());
    // mu_v > mu_h: w's row 40 * 2e-6 + 40 * 4e-5 = 1.68e-3 beats u's 3e-5 + 8e-4, and the metric
    // term uses the larger viscosity, 2 * 40 * 0.01 * 4e-4 = 3.2e-4 -> 2e-3
    EXPECT_NEAR(MomentumDiffusiveRate(Real(1), Real(10), Real(40), a, b, h2, dzinv2, Real(0), Real(0)),
                Real(2e-3), Real(2e-3)*rel_tol());
    // A PBL alone (mu_h = 0) on flat ground with w implicit (ERF_IMPLICIT_W): the horizontal
    // diffusion of w and the cross derivatives are left, 10 * (2e-6 + 4e-5) = 4.2e-4, not zero
    EXPECT_NEAR(MomentumDiffusiveRate(Real(1), Real(0), Real(10), a, b, Real(0), dzinv2, Real(0), Real(0)),
                Real(4.2e-4), Real(4.2e-4)*rel_tol());
    // Unequal spacings: u's row mu_h (2a + b) with a = 4e-6, b = 1e-6 -> 9e-4
    EXPECT_NEAR(MomentumDiffusiveRate(Real(1), Real(100), Real(0), Real(4e-6), Real(1e-6), Real(0), dzinv2, Real(1), Real(1)),
                Real(9e-4), Real(9e-4)*rel_tol());
    // Scalar, K_h = 300, K_v = 30: 300 * 2e-6 + 300 * 0.01 * 4e-4 + 30 * 4e-4 = 1.38e-2
    EXPECT_NEAR(ScalarDiffusiveRate(Real(1), Real(300), Real(30), a, b, h2, dzinv2, Real(1)),
                Real(1.38e-2), Real(1.38e-2)*rel_tol());
    EXPECT_NEAR(ScalarDiffusiveRate(Real(1), Real(300), Real(30), a, b, h2, dzinv2, Real(0.5)),
                Real(7.8e-3), Real(7.8e-3)*rel_tol());
}

TEST(Smag2DLimiters, NeighbourhoodBoundCoversVaryingCoefficients)
{
    // Two cells A, B on a periodic 1-D line, dx = 1 (a counterexample raised in review):
    // K_A = 0.1, K_B = 10, rho_A = 0.1, rho_B = 1.  With arithmetic face averages the operator
    // is L = [[-2 K_f/rho_A, 2 K_f/rho_A], [2 K_f/rho_B, -2 K_f/rho_B]], K_f = (K_A + K_B)/2,
    // whose nonzero eigenvalue is its trace.  Forward Euler is stable for dt |lambda| <= 2;
    // the check's dt * rate <= 1/2 guarantees that only if rate >= |lambda| / 4.
    const Real KA = Real(0.1), KB = Real(10), rA = Real(0.1), rB = Real(1);
    const Real Kf = myhalf * (KA + KB);
    const Real lambda = -(Real(2)*Kf/rA + Real(2)*Kf/rB);   // -111.1
    // The cell-local estimate max(K/rho) = 10 misses it ...
    const Real local = std::max(KA/rA, KB/rB);
    EXPECT_LT(local, std::abs(lambda) / Real(4));
    // ... the neighbourhood bound does not, for either cell.  Each cell's 3x3x3 neighbourhood
    // alternates A and B in x (periodic line), and is uniform in y and z.
    for (int centre = 0; centre < 2; ++centre) {
        std::vector<Real> K(27), rho(27);
        for (int kk = 0; kk < 3; ++kk) {
        for (int jj = 0; jj < 3; ++jj) {
        for (int ii = 0; ii < 3; ++ii) {
            const bool is_A = ((ii + centre) % 2 == 1);
            K  [ii + 3*jj + 9*kk] = is_A ? KA : KB;
            rho[ii + 3*jj + 9*kk] = is_A ? rA : rB;
        }}}
        const amrex::Array4<const Real> Ka(K.data(),   amrex::Dim3{-1,-1,-1}, amrex::Dim3{2,2,2}, 1);
        const amrex::Array4<const Real> ra(rho.data(), amrex::Dim3{-1,-1,-1}, amrex::Dim3{2,2,2}, 1);
        EXPECT_EQ(NeighbourhoodMax(0,0,0,0,Ka), KB);
        EXPECT_EQ(NeighbourhoodMinPositive(0,0,0,0,ra,Real(1)), rA);
        const Real rate = ScalarDiffusiveRate(NeighbourhoodMinPositive(0,0,0,0,ra,Real(1)),
                                              NeighbourhoodMax(0,0,0,0,Ka),
                                              Real(0), Real(1), Real(0), Real(0), Real(0), Real(0));
        EXPECT_GE(rate, std::abs(lambda) / Real(4));
    }
    // Ghost cells holding no density are skipped, and the fallback is used if none is positive
    std::vector<Real> z(27, Real(0));
    z[13] = Real(0.7);
    const amrex::Array4<const Real> za(z.data(), amrex::Dim3{-1,-1,-1}, amrex::Dim3{2,2,2}, 1);
    EXPECT_EQ(NeighbourhoodMinPositive(0,0,0,0,za,Real(5)), Real(0.7));
    std::vector<Real> none(27, Real(0));
    const amrex::Array4<const Real> na(none.data(), amrex::Dim3{-1,-1,-1}, amrex::Dim3{2,2,2}, 1);
    EXPECT_EQ(NeighbourhoodMinPositive(0,0,0,0,na,Real(5)), Real(5));
}

TEST(Smag2DLimiters, CellRatesUseTheNeighbourhood)
{
    // The cell routine the time-step reduction calls: a large diffusivity in one neighbour, or
    // a small density in one, must set the rate, component by component; a covered cell
    // (rho = 0) gives zero.  Flat 50 m cells, 1 km spacing.
    const int nk = EddyDiff::NumDiffs;
    std::vector<Real> rho(27, Real(1)), K(27*static_cast<std::size_t>(nk), Real(1));
    std::vector<Real> z(64);
    for (int kk = 0; kk < 4; ++kk) { for (int jj = 0; jj < 4; ++jj) { for (int ii = 0; ii < 4; ++ii) {
        z[ii + 4*jj + 16*kk] = Real(50) * Real(kk - 1);
    }}}
    const amrex::Dim3 lo{-1,-1,-1};
    const amrex::Array4<const Real> sa(rho.data(), lo, amrex::Dim3{2,2,2}, 1);
    const amrex::Array4<const Real> ka(K.data(),   lo, amrex::Dim3{2,2,2}, nk);
    const amrex::Array4<const Real> za(z.data(),   lo, amrex::Dim3{3,3,3}, 1);
    const amrex::GpuArray<Real,AMREX_SPACEDIM> dxinv{Real(1e-3), Real(1e-3), Real(1)/Real(50)};
    // all vertical parts explicit, so that every diffusivity (K_v too) enters the rate
    DiffusiveRateSettings set;
    set.variable_dz = true; set.has_q = true; set.has_ke = true; set.has_scalar = true;
    set.e_uv = Real(1); set.e_w = Real(1); set.e_th = Real(1); set.e_q = Real(1); set.e_ke = Real(1);
    const Real a = Real(1e-6), dz2 = Real(1)/Real(2500);
    auto rates = [&] (Real& rm, Real& rs) {
        CellDiffusiveRates(0, 0, 0, sa, ka, za, Real(1), Real(1), dxinv, set, rm, rs);
    };
    // a spike in one corner neighbour: (1,1,1) is index 26 of each component's 27 cells,
    // (-1,-1,-1) index 0, so both ends of every neighbourhood loop are exercised
    auto put = [&] (int comp, int corner, Real v) {
        K[static_cast<std::size_t>(corner + 27*comp)] = v;
    };
    Real rm0, rs0;
    rates(rm0, rs0);
    EXPECT_NEAR(rm0, MomentumDiffusiveRate(Real(1), Real(1), Real(1), a, a, Real(0), dz2, Real(1), Real(1)),
                rm0*rel_tol());
    // each component's neighbour value must reach the rate
    struct Case { int comp; bool mom; Real expect; };
    const Case cases[] = {
        {EddyDiff::Mom_h,    true,  MomentumDiffusiveRate(Real(1), Real(1e4), Real(1), a, a, Real(0), dz2, Real(1), Real(1))},
        {EddyDiff::Mom_v,    true,  MomentumDiffusiveRate(Real(1), Real(1), Real(1e4), a, a, Real(0), dz2, Real(1), Real(1))},
        {EddyDiff::Theta_h,  false, ScalarDiffusiveRate(Real(1), Real(1e4), Real(1), a, a, Real(0), dz2, Real(1))},
        {EddyDiff::Theta_v,  false, ScalarDiffusiveRate(Real(1), Real(1), Real(1e4), a, a, Real(0), dz2, Real(1))},
        {EddyDiff::Q_h,      false, ScalarDiffusiveRate(Real(1), Real(1e4), Real(1), a, a, Real(0), dz2, Real(1))},
        {EddyDiff::Q_v,      false, ScalarDiffusiveRate(Real(1), Real(1), Real(1e4), a, a, Real(0), dz2, Real(1))},
        {EddyDiff::KE_h,     false, ScalarDiffusiveRate(Real(1), Real(1e4), Real(1), a, a, Real(0), dz2, Real(1))},
        {EddyDiff::KE_v,     false, ScalarDiffusiveRate(Real(1), Real(1), Real(1e4), a, a, Real(0), dz2, Real(1))},
        {EddyDiff::Scalar_h, false, ScalarDiffusiveRate(Real(1), Real(1e4), Real(1), a, a, Real(0), dz2, Real(1))},
        {EddyDiff::Scalar_v, false, ScalarDiffusiveRate(Real(1), Real(1), Real(1e4), a, a, Real(0), dz2, Real(1))},
    };
    for (const auto& c : cases) {
        for (const int corner : {26, 0}) {
            put(c.comp, corner, Real(1e4));
            Real rm, rs;
            rates(rm, rs);
            const Real got = c.mom ? rm : rs;
            EXPECT_NEAR(got, c.expect, c.expect*rel_tol()) << "component " << c.comp << " corner " << corner;
            EXPECT_GT(got, Real(10) * (c.mom ? rm0 : rs0)) << "component " << c.comp << " corner " << corner;
            put(c.comp, corner, Real(1));
        }
    }
    // with theta's vertical part implicit (e_th = 0) a Theta_v spike no longer enters
    set.e_th = Real(0);
    {
        Real rm_imp, rs_imp, rm_spk, rs_spk;
        rates(rm_imp, rs_imp);
        put(EddyDiff::Theta_v, 26, Real(1e4));
        rates(rm_spk, rs_spk);
        put(EddyDiff::Theta_v, 26, Real(1));
        EXPECT_EQ(rs_spk, rs_imp);
    }
    set.e_th = Real(1);
    // a low density in one neighbour scales every rate by 1/rho_min (either corner)
    Real rm, rs;
    for (const int corner : {26, 0}) {
        rho[corner] = Real(0.1);
        rates(rm, rs);
        EXPECT_NEAR(rm, Real(10)*rm0, Real(10)*rm0*rel_tol()) << "corner " << corner;
        EXPECT_NEAR(rs, Real(10)*rs0, Real(10)*rs0*rel_tol()) << "corner " << corner;
        rho[corner] = Real(1);
    }
    // a covered cell contributes nothing
    rho[13] = Real(0);
    rates(rm, rs);
    EXPECT_EQ(rm, Real(0));
    EXPECT_EQ(rs, Real(0));
}

namespace {

// The horizontal diffusion of w through the production stress and momentum-source code, on a
// flat periodic box with a PBL-like viscosity (Mom_h = 0, Mom_v = mu), rho = 1 and
// w = (-1)^(i+j): tau13 and tau23 hold the strains S13 = (dw/dx)/2 = (-1)^(i+j)/dx and
// S23 = (dw/dy)/2 = (-1)^(i+j)/dy (u = v = 0), which ComputeStressVarVisc_N turns into the
// stresses and DiffusionSrcForMom into rho_w_rhs.  The mode is an eigenvector of the
// operator with eigenvalue -4 mu (1/dx^2 + 1/dy^2).  Returns the largest
// |rho_w_rhs - lambda w| over the interior w-faces, and lambda.
Real w_checkerboard_residual (Real mu, Real dx, Real& lambda)
{
    const amrex::Box cells(amrex::IntVect(0), amrex::IntVect(7, 3, 3));
    const amrex::BoxArray ba(cells);
    const amrex::DistributionMapping dm(ba);
    const amrex::GpuArray<Real,AMREX_SPACEDIM> dxInv{one/dx, one/dx, one/dx};
    lambda = -Real(8) * mu / (dx*dx);   // dy = dx

    amrex::MultiFab K(ba, dm, EddyDiff::NumDiffs, 1);
    K.setVal(Real(0));
    K.setVal(mu, EddyDiff::Mom_v, 1, 1);
    amrex::MultiFab er(ba, dm, 1, 1);  er.setVal(Real(0));
    amrex::MultiFab t11(ba, dm, 1, 1), t22(ba, dm, 1, 1), t33(ba, dm, 1, 1);
    amrex::MultiFab t12(amrex::convert(ba, amrex::IntVect(1,1,0)), dm, 1, 1);
    amrex::MultiFab t13(amrex::convert(ba, amrex::IntVect(1,0,1)), dm, 1, 1);
    amrex::MultiFab t23(amrex::convert(ba, amrex::IntVect(0,1,1)), dm, 1, 1);
    for (auto* m : {&t11, &t22, &t33, &t12}) { m->setVal(Real(0)); }
    amrex::MultiFab ru(amrex::convert(ba, amrex::IntVect(1,0,0)), dm, 1, 0);
    amrex::MultiFab rv(amrex::convert(ba, amrex::IntVect(0,1,0)), dm, 1, 0);
    amrex::MultiFab rw(amrex::convert(ba, amrex::IntVect(0,0,1)), dm, 1, 0);
    ru.setVal(Real(0)); rv.setVal(Real(0)); rw.setVal(Real(0));
    amrex::MultiFab detJ(ba, dm, 1, 1);  detJ.setVal(one);
    amrex::MultiFab mf(ba, dm, 1, 1);    mf.setVal(one);
    amrex::Gpu::DeviceVector<Real> no_dz;

    const Real dxi = one/dx;
    for (amrex::MFIter mfi(t13); mfi.isValid(); ++mfi) {
        const auto s13 = t13.array(mfi);
        const auto s23 = t23.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int p = ((i + j) % 2 + 2) % 2;
            s13(i,j,k) = (p == 0 ? one : -one) * dxi;
        });
        amrex::ParallelFor(t23[mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int p = ((i + j) % 2 + 2) % 2;
            s23(i,j,k) = (p == 0 ? one : -one) * dxi;
        });
    }
    for (amrex::MFIter mfi(K); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        auto a11 = t11.array(mfi); auto a22 = t22.array(mfi); auto a33 = t33.array(mfi);
        auto a12 = t12.array(mfi); auto a13 = t13.array(mfi); auto a23 = t23.array(mfi);
        amrex::Array4<Real> none{};
        ComputeStressVarVisc_N(bx, amrex::surroundingNodes(amrex::surroundingNodes(bx,0),1),
                               amrex::surroundingNodes(amrex::surroundingNodes(bx,0),2),
                               amrex::surroundingNodes(amrex::surroundingNodes(bx,1),2),
                               Real(0), K.const_array(mfi), amrex::Array4<const Real>{},
                               a11, a22, a33, a12, a13, a23, er.const_array(mfi), none, none, none);
        DiffusionSrcForMom(amrex::surroundingNodes(bx,0), amrex::surroundingNodes(bx,1),
                           amrex::surroundingNodes(bx,2),
                           ru.array(mfi), rv.array(mfi), rw.array(mfi),
                           t11.const_array(mfi), t22.const_array(mfi), t33.const_array(mfi),
                           t12.const_array(mfi), t12.const_array(mfi),
                           t13.const_array(mfi), t13.const_array(mfi),
                           t23.const_array(mfi), t23.const_array(mfi),
                           detJ.const_array(mfi), no_dz, dxInv,
                           mf.const_array(mfi), mf.const_array(mfi), mf.const_array(mfi),
                           mf.const_array(mfi), mf.const_array(mfi), mf.const_array(mfi),
                           false, false);
    }
    amrex::ReduceOps<amrex::ReduceOpMax> op;
    amrex::ReduceData<Real> data(op);
    using T = typename decltype(data)::Type;
    const Real lam = lambda;
    for (amrex::MFIter mfi(rw); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.validbox();
        bx.grow(2, -1);   // interior w-faces
        const auto r = rw.const_array(mfi);
        op.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> T
        {
            const int p = ((i + j) % 2 + 2) % 2;
            const Real w = (p == 0) ? one : -one;
            return {std::abs(r(i,j,k) - lam * w)};
        });
    }
    return amrex::get<0>(data.value(op));
}

} // namespace

TEST(Smag2DLimiters, RateBoundsTheHorizontalDiffusionOfW)
{
    // A PBL alone (Mom_h = 0) on flat ground, with every vertical part implicit (as with
    // ERF_IMPLICIT_W): the production operator diffuses w horizontally with Mom_v, and the
    // check's rate must bound that eigenvalue, 4 * rate >= |lambda|.
    const Real mu = Real(20), dx = Real(1000);
    Real lambda = Real(0);
    const Real res = w_checkerboard_residual(mu, dx, lambda);
    EXPECT_LT(res, std::abs(lambda) * Real(100) * rel_tol());   // the mode is an eigenvector
    EXPECT_GT(std::abs(lambda), Real(0));
    // The mode's eigenvalue (-8 mu a, both horizontal directions) does not depend on dz; take
    // dz = 4 dx so that the cross-derivative allowance, mu_v (sqrt(a c) + sqrt(b c)) = mu_v a/2,
    // cannot cover it alone: the w row mu_v (a + b) must, with both directions counted.
    const Real a = one / (dx*dx);
    const Real c = a / Real(16);
    const Real rate = MomentumDiffusiveRate(one, Real(0), mu, a, a, Real(0), c, Real(0), Real(0));
    EXPECT_GE(Real(4) * rate, std::abs(lambda) * (one - Real(10)*rel_tol()));
}

namespace {

// Largest |S_k| over every stage of the partly implicit scheme for y = -lambda dt in
// [1e-6, 1e12], sampled densely (2000 points per decade); independent of the classifier.
double max_stage_amplification (const std::vector<double>& c, const std::vector<double>& f)
{
    double m = 0.0;
    for (int n = 0; n <= 36000; ++n) {
        const double y = std::pow(10.0, -6.0 + n * 5.0e-4);
        double S = 1.0;
        for (std::size_t k = 0; k < c.size(); ++k) {
            S = (1.0 - (1.0 - f[k]) * c[k] * y * S) / (1.0 + f[k] * c[k] * y);
            m = std::max(m, std::abs(S));
        }
    }
    return m;
}

} // namespace

TEST(Smag2DLimiters, SlopeReportOnNewGridsAndMovingTerrain)
{
    // evaluated on new grids, and every step on a moving terrain; not otherwise
    EXPECT_TRUE (SlopeReportEvaluate(true,  false));
    EXPECT_TRUE (SlopeReportEvaluate(false, true));
    EXPECT_FALSE(SlopeReportEvaluate(false, false));
    // printed when forced, or when alpha moved by more than 1 % since the last print
    EXPECT_TRUE (SlopeReportPrint(true,  Real(1.5),    Real(1.5)));
    EXPECT_FALSE(SlopeReportPrint(false, Real(1.5),    Real(1.5)));
    EXPECT_FALSE(SlopeReportPrint(false, Real(1.5149), Real(1.5)));
    EXPECT_TRUE (SlopeReportPrint(false, Real(1.5151), Real(1.5)));
    EXPECT_TRUE (SlopeReportPrint(false, Real(1.4849), Real(1.5)));
    EXPECT_FALSE(SlopeReportPrint(false, Real(1.0099), Real(1)));   // alpha >= 1: floor 0.01
    EXPECT_TRUE (SlopeReportPrint(false, Real(1.0101), Real(1)));
    EXPECT_TRUE (SlopeReportPrint(false, Real(1),      Real(0)));   // nothing printed yet
}

TEST(Smag2DLimiters, ExplicitVerticalFractionFromStability)
{
    using V = std::vector<Real>;
    const std::vector<double> rk3{1.0/3.0, 0.5, 1.0}, mid{0.5, 1.0};
    struct Case { int nst; std::vector<double> f; Real e; };
    const Case cases[] = {
        // bounded (e = 0): f1 >= 1/2 with f2 = 1 (any f3); MidPoint with f1 = 1 (any f2)
        {3, {1, 1, 0}, 0}, {3, {1, 1, 1}, 0}, {3, {1, 1, 0.5}, 0}, {3, {0.5, 1, 0}, 0},
        {3, {0.5, 1, 1}, 0}, {2, {1, 0, 0}, 0}, {2, {1, 1, 0}, 0},
        // everything else is explicit (e = 1)
        {3, {0, 0, 0}, 1}, {3, {1, 0, 0}, 1}, {3, {1, 0, 0.25}, 1}, {3, {0.5, 0.5, 0.25}, 1},
        {3, {0, 0, 1}, 1}, {3, {0, 1, 0.5}, 1}, {3, {0.49, 1, 0}, 1}, {2, {0, 0, 0}, 1},
        {2, {0.9, 0, 0}, 1},
        // a review counterexample: the last stage reaches 1.0008 near y = 23.5, between the
        // points an earlier sampled classifier looked at
        {3, {0.995, 0.1, 0.0666}, 1},
    };
    for (const auto& c : cases) {
        V fac(c.f.begin(), c.f.end());
        EXPECT_EQ(ExplicitVerticalFraction(c.nst, fac, true), c.e)
            << c.f[0] << " " << c.f[1] << " " << c.f[2];
        // every pattern classified implicit must have bounded stages
        if (c.e == Real(0)) {
            const std::vector<double> ff(c.f.begin(), c.f.begin() + c.nst);
            EXPECT_LE(max_stage_amplification(c.nst == 3 ? rk3 : mid, ff), 1.0 + 1.0e-12)
                << c.f[0] << " " << c.f[1] << " " << c.f[2];
        }
    }
    // the counterexample really does exceed 1, and the f1 boundary is tight; the rule is
    // sufficient, not necessary: 1 0 1 keeps every stage bounded but is classified explicit
    EXPECT_GT(max_stage_amplification(rk3, {0.995, 0.1, 0.0666}), 1.0005);
    EXPECT_GT(max_stage_amplification(rk3, {0.49, 1, 0}), 1.01);
    EXPECT_GT(max_stage_amplification(mid, {0.9, 0}), 2.0);
    EXPECT_LE(max_stage_amplification(rk3, {1, 0, 1}), 1.0 + 1.0e-12);
    EXPECT_EQ(ExplicitVerticalFraction(3, V{1, 0, 1}, true), Real(1));
    // the component left out of the implicit solve is explicit whatever the factors
    EXPECT_EQ(ExplicitVerticalFraction(3, V{1, 1, 1}, false), Real(1));
    EXPECT_EQ(ExplicitVerticalFraction(3, V{1, 1}, true), Real(1));   // too few factors
    // factors outside [0,1] are not clamped into a proven pattern: 1 1.5 0 amplifies by 4/3
    EXPECT_EQ(ExplicitVerticalFraction(3, V{1, 1.5, 0}, true), Real(1));
    EXPECT_GT(max_stage_amplification(rk3, {1, 1.5, 0}), 1.3);
    EXPECT_EQ(ExplicitVerticalFraction(3, V{1, 1, -0.5}, true), Real(1));
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

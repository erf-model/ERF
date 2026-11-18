// Unit tests for the terrain-following inflow profiles in ERF_InflowProfile.H:
// the table parser, the clamped lookup and the log-law tabulation that the
// boundary fills use, and the height above the ground on a sloping
// terrain-fitted mesh.

#include <ERF_InflowProfile.H>

#include <AMReX_BaseFab.H>
#include <AMReX_Box.H>

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using amrex::Real;

namespace {

constexpr Real kappa = Real(0.41);

Real tol (Real factor = Real(1e3))
{
    return factor * std::numeric_limits<Real>::epsilon();
}

std::string parse (const std::string& text, InflowProfile& prof)
{
    std::istringstream is(text);
    return parse_inflow_profile(is, prof);
}

Real lookup (const InflowProfile& p, const amrex::Vector<Real>& f, const Real zq)
{
    return inflow_profile_value(p.z.data(), f.data(), static_cast<int>(p.z.size()), zq);
}

InflowLogLaw westerly ()
{
    InflowLogLaw ll;
    ll.speed     = Real(20.0);
    ll.height    = Real(10.0);
    ll.direction = Real(285.0);
    ll.z0        = Real(0.1);
    return ll;
}

} // namespace

TEST(InflowProfile, HeaderlessFourColumnsAreZUVT)
{
    InflowProfile p;
    ASSERT_EQ(parse("0 1 2 300\n100 3 4 301\n", p), "");
    EXPECT_TRUE(p.active);
    EXPECT_TRUE(p.has_theta);
    EXPECT_FALSE(p.has_tke);
    ASSERT_EQ(p.z.size(), 2);
    EXPECT_EQ(p.z[1], Real(100));
    EXPECT_EQ(p.u[1], Real(3));
    EXPECT_EQ(p.v[0], Real(2));
    EXPECT_EQ(p.theta[1], Real(301));
}

TEST(InflowProfile, HeaderlessFiveColumnsAddTke)
{
    InflowProfile p;
    ASSERT_EQ(parse("0 1 2 300 0.5\n100 3 4 301 0.25\n", p), "");
    EXPECT_TRUE(p.has_theta);
    EXPECT_TRUE(p.has_tke);
    EXPECT_EQ(p.tke[0], Real(0.5));
    EXPECT_EQ(p.tke[1], Real(0.25));
}

TEST(InflowProfile, HeaderNamesColumnsInAnyOrder)
{
    InflowProfile p;
    ASSERT_EQ(parse("# z tke w v u\n0 0.4 0 -1 8\n200 0.3 0 1 9\n", p), "");
    EXPECT_FALSE(p.has_theta);
    EXPECT_TRUE(p.has_tke);
    EXPECT_EQ(p.u[0], Real(8));
    EXPECT_EQ(p.v[1], Real(1));
    EXPECT_EQ(p.tke[0], Real(0.4));
}

TEST(InflowProfile, FreeTextCommentsBlankLinesAndCaseAreIgnored)
{
    InflowProfile p;
    ASSERT_EQ(parse("# an inflow profile\n\n# Z U V Theta\n0 1 2 300\n   \n50 2 3 300\n", p), "");
    EXPECT_TRUE(p.has_theta);
    ASSERT_EQ(p.z.size(), 2);
    EXPECT_EQ(p.u[1], Real(2));
}

TEST(InflowProfile, RejectsMalformedTables)
{
    const std::vector<std::string> bad = {
        "",                                   // no rows
        "0 1 2 300\n",                        // one height
        "0 1 2\n1 1 2\n",                     // three columns without a header
        "0 1 2 300\n0 3 4 301\n",             // heights not increasing
        "# z u v q\n0 1 2 3\n1 1 2 3\n",      // unknown column
        "# z u u v\n0 1 2 3\n1 1 2 3\n",      // repeated column
        "# z u T\n0 1 300\n1 1 300\n",        // no v
        "# z u v\n# z u v\n0 1 2\n1 1 2\n",   // header twice
        "0 1 2 300\n1 3 4\n",                 // short row
        "0 1 2 abc\n1 1 2 3\n",               // not a number
        "0 1 2 0\n1 1 2 300\n",               // T not positive
        "0 1 2 300 -1\n1 1 2 300 0\n",        // negative tke
    };
    for (const auto& text : bad) {
        InflowProfile p;
        EXPECT_NE(parse(text, p), "") << "accepted:\n" << text;
        EXPECT_FALSE(p.active) << text;
    }
}

TEST(InflowProfile, LookupInterpolatesAndHoldsBeyondTheEnds)
{
    const Real z[] = {0, 10, 20};
    const Real f[] = {0, 1, 4};
    EXPECT_NEAR(inflow_profile_value(z, f, 3, Real(5)),  Real(0.5), tol());
    EXPECT_NEAR(inflow_profile_value(z, f, 3, Real(15)), Real(2.5), tol());
    EXPECT_EQ(inflow_profile_value(z, f, 3, Real(10)), Real(1));
    EXPECT_EQ(inflow_profile_value(z, f, 3, Real(-3)), Real(0));
    EXPECT_EQ(inflow_profile_value(z, f, 3, Real(30)), Real(4));
}

TEST(InflowProfile, LogLawGivesTheReferenceSpeedAndDirection)
{
    InflowProfile p;
    const InflowLogLaw ll = westerly();
    ASSERT_EQ(make_log_law_inflow_profile(ll, Real(4000), p), "");
    EXPECT_TRUE(p.active);
    EXPECT_TRUE(p.has_tke);
    EXPECT_FALSE(p.has_theta);
    EXPECT_EQ(p.z[0], Real(0));
    EXPECT_NEAR(p.z.back(), Real(4000), Real(1e-9) * Real(4000));

    const Real u = lookup(p, p.u, ll.height);
    const Real v = lookup(p, p.v, ll.height);
    EXPECT_NEAR(std::hypot(u, v), ll.speed, Real(2e-3));
    // From 285 degrees: towards 105 degrees, mostly +x with a small -y part
    const Real dir = ll.direction * std::acos(Real(-1)) / Real(180);
    EXPECT_NEAR(u, -ll.speed * std::sin(dir), Real(2e-3));
    EXPECT_NEAR(v, -ll.speed * std::cos(dir), Real(2e-3));
    // Zero at the ground, increasing with height
    EXPECT_EQ(lookup(p, p.u, Real(0)), Real(0));
    EXPECT_LT(lookup(p, p.u, Real(50)), lookup(p, p.u, Real(500)));
}

TEST(InflowProfile, LogLawTkeStartsAtTheWallValueAndTapers)
{
    InflowProfile p;
    const InflowLogLaw ll = westerly();
    ASSERT_EQ(make_log_law_inflow_profile(ll, Real(4000), p), "");
    const Real ustar = kappa * ll.speed / std::log((ll.height + ll.z0) / ll.z0);
    const Real tke0  = ustar * ustar / (ll.Cmu0 * ll.Cmu0);
    const Real depth = ustar * ll.tke_zscale;
    EXPECT_NEAR(p.tke[0], tke0, tol() * tke0);
    EXPECT_NEAR(lookup(p, p.tke, Real(0.5) * depth), Real(0.5) * tke0, Real(1e-6) * tke0);
    EXPECT_NEAR(lookup(p, p.tke, Real(1.5) * depth), Real(0.01) * tke0, tol() * tke0);
}

TEST(InflowProfile, LogLawCapsTheSpeed)
{
    InflowProfile p;
    InflowLogLaw ll = westerly();
    ll.max_speed = Real(25.0);
    ASSERT_EQ(make_log_law_inflow_profile(ll, Real(4000), p), "");
    const Real top = std::hypot(lookup(p, p.u, Real(4000)), lookup(p, p.v, Real(4000)));
    EXPECT_NEAR(top, Real(25.0), tol(Real(1e4)));
    EXPECT_NEAR(std::hypot(lookup(p, p.u, ll.height), lookup(p, p.v, ll.height)), ll.speed, Real(2e-3));
}

TEST(InflowProfile, LogLawRejectsBadInputs)
{
    InflowProfile p;
    InflowLogLaw ll = westerly();
    ll.speed = Real(0);
    EXPECT_NE(make_log_law_inflow_profile(ll, Real(1000), p), "");
    ll = westerly(); ll.z0 = Real(0);
    EXPECT_NE(make_log_law_inflow_profile(ll, Real(1000), p), "");
    ll = westerly(); ll.height = Real(-1);
    EXPECT_NE(make_log_law_inflow_profile(ll, Real(1000), p), "");
    ll = westerly();
    EXPECT_NE(make_log_law_inflow_profile(ll, Real(0), p), "");
    EXPECT_FALSE(p.active);
}

TEST(InflowProfile, HeightAboveGroundFollowsTheLocalColumn)
{
    // Nodes (0..2) x (0..1) x (0..2). The ground rises 10 m per node in x and
    // 5 m per node in y; basic terrain following with a 100 m top compresses
    // each column's levels (4 m apart over flat ground) by (1 - ground/100).
    const amrex::Box nodes(amrex::IntVect(0, 0, 0), amrex::IntVect(2, 1, 2));
    std::vector<Real> data(static_cast<std::size_t>(nodes.numPts()));
    auto a = amrex::makeArray4<Real>(data.data(), nodes, 1);
    auto ground = [] (int i, int j) { return Real(10 * i + 5 * j); };
    for (int k = 0; k <= 2; ++k) {
        for (int j = 0; j <= 1; ++j) {
            for (int i = 0; i <= 2; ++i) {
                const Real g = ground(i, j);
                a(i,j,k) = g + (Real(1) - g / Real(100)) * Real(4 * k);
            }
        }
    }
    const auto c = amrex::makeArray4<const Real>(data.data(), nodes, 1);

    // x-face at i = 2 (ground 20 and 25 m), first level: half of the compressed spacing
    const Real xface = Real(0.5) * (Real(0.5) * Real(4) * Real(0.80) + Real(0.5) * Real(4) * Real(0.75));
    EXPECT_NEAR(height_above_ground(c, 2, 2, 0, 1, 0, 0), xface, tol(Real(1e4)));

    // y-face at j = 0 over i = 0..1 (ground 0 and 10 m), second level
    const Real yface = Real(0.5) * (Real(6) * Real(1.0) + Real(6) * Real(0.9));
    EXPECT_NEAR(height_above_ground(c, 0, 1, 0, 0, 1, 0), yface, tol(Real(1e4)));

    // Indices beyond the array clamp to its edge: level 5 is the top pair of levels
    EXPECT_NEAR(height_above_ground(c, 7, 7, 3, 3, 5, 0), Real(0.5) * (Real(4) + Real(8)) * Real(0.75),
                tol(Real(1e4)));
}

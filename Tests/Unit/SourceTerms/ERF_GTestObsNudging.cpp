#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include "ERF_IndexDefines.H"
#include "ERF_ObsNudging.H"
#include "ERF_ObsNudgingKernel.H"
#include "ERF_ObsNudgingSeries.H"

// Observation nudging: the station file, the profiles it gives at a moment of
// the run, the tendency at a point, where the tendency lands on the staggered
// grid over flat ground, a fitted mesh and immersed terrain.  The placement of a
// latitude/longitude and the grid rotation are tested in Utils/ERF_GTestLatLonMap.cpp.

namespace {

using amrex::Real;
using namespace obs_nudging;

constexpr Real tol = (sizeof(Real) == 8) ? Real(1.0e-12) : Real(1.0e-5);

bool parse (const std::string& text, StationSeries& series, std::string& error,
            Real missing = Real(-9999.0))
{
    std::istringstream is(text);
    return parse_station_series(is, "test.txt", missing, series, error);
}

// A view of one station at (x0, y0) with one profile per component, built from
// host arrays.  Only host code reads it in these tests.
struct HostView
{
    amrex::Vector<Real> x, y, z, mean, sigma;
    amrex::Vector<int> agl, offset, count;
    ObsNudgingView v;

    void finish (Real rh, Real rz, Real cutoff, Real alpha, Real tau, Real dt)
    {
        v.nstations = static_cast<int>(x.size());
        v.x = x.data(); v.y = y.data(); v.agl = agl.data();
        v.offset = offset.data(); v.count = count.data();
        v.z = z.data(); v.mean = mean.data(); v.sigma = sigma.data();
        v.inv_rh2 = Real(1.0)/(rh*rh);
        v.inv_rz2 = Real(1.0)/(rz*rz);
        v.qmax = cutoff*cutoff;
        v.sigma_factor = alpha;
        v.inv_tau = Real(1.0)/tau;
        v.max_rate = Real(1.0)/dt;
    }
};

// Stations with the same profile of component c = 0 (u) only
HostView
make_view (const amrex::Vector<std::array<Real,2>>& xy,
           const amrex::Vector<Real>& z, const amrex::Vector<Real>& mean,
           const amrex::Vector<Real>& sigma, bool agl = true)
{
    HostView h;
    const int ns = static_cast<int>(xy.size());
    h.offset.assign(NComp*ns, 0);
    h.count.assign(NComp*ns, 0);
    for (int s = 0; s < ns; ++s) {
        h.x.push_back(xy[s][0]);
        h.y.push_back(xy[s][1]);
        h.agl.push_back(agl ? 1 : 0);
        h.offset[U*ns + s] = static_cast<int>(h.z.size());
        h.count [U*ns + s] = static_cast<int>(z.size());
        h.z.insert(h.z.end(), z.begin(), z.end());
        h.mean.insert(h.mean.end(), mean.begin(), mean.end());
        h.sigma.insert(h.sigma.end(), sigma.begin(), sigma.end());
    }
    return h;
}

void
write_file (const std::string& name, const std::string& text)
{
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(name);
        os << text;
    }
    amrex::ParallelDescriptor::Barrier();
}

// z_phys_nd of a column of uniform dz raised by the terrain height h(x) = a + b x
void
fill_fitted_mesh (amrex::MultiFab& znd, const amrex::Geometry& geom, Real a, Real b)
{
    const auto dx = geom.CellSizeArray();
    const auto lo = geom.ProbLoArray();
    for (amrex::MFIter mfi(znd); mfi.isValid(); ++mfi) {
        const auto z = znd.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            amrex::ignore_unused(j);
            const Real x = lo[0] + static_cast<Real>(i) * dx[0];
            z(i,j,k) = a + b*x + static_cast<Real>(k) * dx[2];
        });
    }
}

// The cell-centred blanking of immersed terrain: solid (1) below k = kb, half
// at kb, fluid above
void
fill_blank (amrex::MultiFab& blank, int kb)
{
    for (amrex::MFIter mfi(blank); mfi.isValid(); ++mfi) {
        const auto b = blank.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            amrex::ignore_unused(i, j);
            b(i,j,k) = (k < kb) ? Real(1.0) : ((k == kb) ? Real(0.5) : Real(0.0));
        });
    }
}

// The inputs of a run with one station at the centre of the test domain
void
set_inputs (const std::string& file, Real rh, Real rz, Real tau)
{
    amrex::ParmParse pp("erf.obs_nudging");
    pp.add("tau", tau);
    pp.add("horizontal_radius", rh);
    pp.add("vertical_radius", rz);
    pp.add("cutoff", Real(1.0e3));
    pp.add("sigma_factor", Real(0.0));
    pp.addarr("stations", std::vector<std::string>{"unit"});
    amrex::ParmParse pps("erf.obs_nudging.unit");
    pps.add("file", file);
    pps.add("x", Real(400.0));
    pps.add("y", Real(400.0));
    pps.add("wind_frame", std::string("grid"));
}

} // namespace

// ---------------------------------------------------------------------------
// The station file
// ---------------------------------------------------------------------------

TEST(ObsNudgingSeries, ParsesProfilesAndMissingValues)
{
    StationSeries s;
    std::string err;
    ASSERT_TRUE(parse("# a comment line\n"
                      "time z u v theta su   # and a trailing comment\n"
                      "0   10  5  1   300   0.5\n"
                      "0   20  nan 1  301   0.5\n"
                      "60  10  6  2   -9999 0.25\n"
                      "60  20  7  2   302   nan\n", s, err)) << err;
    ASSERT_EQ(s.heights.size(), 2u);
    ASSERT_EQ(s.records.size(), 2u);
    EXPECT_TRUE(s.has[U] && s.has[V] && s.has[Theta]);
    EXPECT_FALSE(s.has[W]);
    EXPECT_EQ(s.records[1].time, 60.0);
    EXPECT_NEAR(s.records[0].mean[U][0], Real(5.0), tol);
    // u missing makes the wind missing at that height
    EXPECT_TRUE(std::isnan(s.records[0].mean[U][1]));
    EXPECT_TRUE(std::isnan(s.records[0].mean[V][1]));
    // the missing value marks theta missing
    EXPECT_TRUE(std::isnan(s.records[1].mean[Theta][0]));
    // a missing standard deviation is zero
    EXPECT_EQ(s.records[1].sigma[U][1], Real(0.0));
    EXPECT_NEAR(s.records[1].sigma[U][0], Real(0.25), tol);
}

TEST(ObsNudgingSeries, SpeedAndDirectionAreTheDirectionTheWindComesFrom)
{
    StationSeries s;
    std::string err;
    ASSERT_TRUE(parse("time z speed direction\n0 10 10 270\n0 20 10 180\n", s, err)) << err;
    // From the west: blowing east
    EXPECT_NEAR(s.records[0].mean[U][0], Real(10.0), Real(1.0e-5));
    EXPECT_NEAR(s.records[0].mean[V][0], Real(0.0),  Real(1.0e-5));
    // From the south: blowing north
    EXPECT_NEAR(s.records[0].mean[U][1], Real(0.0),  Real(1.0e-5));
    EXPECT_NEAR(s.records[0].mean[V][1], Real(10.0), Real(1.0e-5));
}

TEST(ObsNudgingSeries, RejectsMalformedFiles)
{
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"time z u v pressure\n0 10 1 1 1\n",              "unknown column 'pressure'"},
        {"time z u\n0 10 1\n",                            "u and v must be given together"},
        {"time u v\n0 1 1\n",                             "must name the columns time and z"},
        {"time z u v speed direction\n0 10 1 1 1 1\n",    "not both"},
        {"time z sw\n0 10 1\n",                           "sw needs the column w"},
        {"time z su\n0 10 1\n",                           "su/sv need a wind column"},
        {"time z\n0 10\n",                                "none of the nudged quantities"},
        {"time z u v\n0 20 1 1\n0 10 1 1\n",              "heights within a time must increase"},
        {"time z u v\n0 10 1 1\n0 20 1 1\n60 10 1 1\n60 30 1 1\n", "does not match the heights"},
        {"time z u v\n0 10 1 1\n0 20 1 1\n60 10 1 1\n",   "every time needs the same heights"},
        {"time z u v\n60 10 1 1\n0 10 1 1\n",             "times must increase"},
        {"time z u v\n0 10 1 x\n",                        "is not a number"},
        {"time z u v\n0 10 1\n",                          "expected 4 entries"},
        {"time z w sw\n0 10 1 -1\n",                      "may not be negative"},
        {"time z theta\n0 10 0\n",                        "theta must be positive"},
        {"time z speed direction\n0 10 1 400\n",          "direction must be in [0, 360]"},
        {"time z u v\nnan 10 1 1\n",                      "time and z may not be missing"},
        {"# nothing\n",                                   "no header line"},
        {"time z u v\n",                                  "no measurements"},
    };
    for (const auto& [text, message] : cases) {
        StationSeries s;
        std::string err;
        EXPECT_FALSE(parse(text, s, err)) << text;
        EXPECT_NE(err.find(message), std::string::npos) << "expected '" << message << "' in '" << err << "'";
    }
}

// ---------------------------------------------------------------------------
// Profiles at a moment of the run
// ---------------------------------------------------------------------------

TEST(ObsNudgingSeries, OneRecordHoldsForTheWholeRun)
{
    StationSeries s;
    std::string err;
    ASSERT_TRUE(parse("time z u v\n100 10 3 4\n", s, err)) << err;
    std::array<Profile, NComp> p;
    for (const double t : {-1.0e6, 0.0, 1.0e6}) {
        ASSERT_TRUE(profiles_at_time(s, t, Real(1.0), Real(0.0), p));
        ASSERT_EQ(p[U].z.size(), 1u);
        EXPECT_NEAR(p[U].mean[0], Real(3.0), tol);
    }
}

TEST(ObsNudgingSeries, InterpolatesInTimeAndIsInactiveOutsideTheRecords)
{
    StationSeries s;
    std::string err;
    ASSERT_TRUE(parse("time z u v w sw\n"
                      "0   10  2 0  nan 0.1\n"
                      "0   20  4 0  0.5 0.1\n"
                      "100 10  6 0  0.3 0.3\n"
                      "100 20  8 0  0.7 0.3\n", s, err)) << err;
    std::array<Profile, NComp> p;

    ASSERT_TRUE(profiles_at_time(s, 25.0, Real(1.0), Real(0.0), p));
    ASSERT_EQ(p[U].z.size(), 2u);
    EXPECT_NEAR(p[U].mean[0], Real(3.0), Real(1.0e-5));
    EXPECT_NEAR(p[U].mean[1], Real(5.0), Real(1.0e-5));
    // w at 10 m is missing at t = 0, so it is not used between the records...
    ASSERT_EQ(p[W].z.size(), 1u);
    EXPECT_NEAR(p[W].z[0], Real(20.0), tol);
    EXPECT_NEAR(p[W].mean[0], Real(0.55), Real(1.0e-5));
    EXPECT_NEAR(p[W].sigma[0], Real(0.15), Real(1.0e-5));
    // ...but it is at the time of the record that has it
    ASSERT_TRUE(profiles_at_time(s, 100.0, Real(1.0), Real(0.0), p));
    EXPECT_EQ(p[W].z.size(), 2u);

    EXPECT_FALSE(profiles_at_time(s, -1.0, Real(1.0), Real(0.0), p));
    EXPECT_TRUE(p[U].z.empty());
    EXPECT_FALSE(profiles_at_time(s, 100.5, Real(1.0), Real(0.0), p));
    EXPECT_TRUE(p[U].z.empty());
}

TEST(ObsNudgingSeries, RotatesAnEarthRelativeWindIntoTheGrid)
{
    StationSeries s;
    std::string err;
    ASSERT_TRUE(parse("time z u v su sv\n0 10 1 0 0.3 0.4\n", s, err)) << err;
    std::array<Profile, NComp> p;
    // Grid x axis 90 degrees counterclockwise from east, i.e. pointing north:
    // an east wind blows along -y of the grid
    ASSERT_TRUE(profiles_at_time(s, 0.0, Real(0.0), Real(1.0), p));
    EXPECT_NEAR(p[U].mean[0], Real( 0.0), tol);
    EXPECT_NEAR(p[V].mean[0], Real(-1.0), tol);
    EXPECT_NEAR(p[U].sigma[0], Real(0.4), tol);
    EXPECT_NEAR(p[V].sigma[0], Real(0.3), tol);
    // 30 degrees
    const Real c = std::cos(amrex::Math::pi<Real>()/Real(6.0)), sn = std::sin(amrex::Math::pi<Real>()/Real(6.0));
    ASSERT_TRUE(profiles_at_time(s, 0.0, c, sn, p));
    EXPECT_NEAR(p[U].mean[0],  c, Real(1.0e-6));
    EXPECT_NEAR(p[V].mean[0], -sn, Real(1.0e-6));
}

// ---------------------------------------------------------------------------
// The tendency at a point
// ---------------------------------------------------------------------------

TEST(ObsNudgingKernel, OneStationAgainstTheClosedForm)
{
    // One height, sigma 0: S = -min(w / tau, 1 / dt) (phi - mean), w = exp(-q/4)
    auto h = make_view({{Real(0.0), Real(0.0)}}, {Real(50.0)}, {Real(8.0)}, {Real(0.0)});
    const Real rh = 200, rz = 20, tau = 30, dt = 1;
    h.finish(rh, rz, Real(100.0), Real(1.0), tau, dt);

    // 150 m away horizontally and 30 m above the measured height
    const Real q = (Real(150.0)*Real(150.0))/(rh*rh) + (Real(30.0)*Real(30.0))/(rz*rz);
    const Real w = std::exp(Real(-0.25)*q);
    const Real phi = Real(5.0);
    EXPECT_NEAR(obs_nudging_tendency(h.v, U, Real(150.0), Real(0.0), Real(80.0), Real(0.0), phi),
                -(w/tau)*(phi - Real(8.0)), tol);
    // Components the station does not measure are not nudged
    EXPECT_EQ(obs_nudging_tendency(h.v, V, Real(0.0), Real(0.0), Real(50.0), Real(0.0), phi), Real(0.0));
}

TEST(ObsNudgingKernel, TwoSeparateStationsEachForceTheirOwnNeighbourhood)
{
    auto h = make_view({{Real(0.0), Real(0.0)}, {Real(5000.0), Real(0.0)}},
                       {Real(50.0)}, {Real(8.0)}, {Real(0.0)});
    h.mean[1] = Real(2.0);
    h.finish(Real(100.0), Real(20.0), Real(6.0), Real(1.0), Real(10.0), Real(1.0));
    const Real phi = Real(5.0);
    EXPECT_NEAR(obs_nudging_tendency(h.v, U, Real(0.0),    Real(0.0), Real(50.0), Real(0.0), phi),
                -(phi - Real(8.0))/Real(10.0), tol);
    EXPECT_NEAR(obs_nudging_tendency(h.v, U, Real(5000.0), Real(0.0), Real(50.0), Real(0.0), phi),
                -(phi - Real(2.0))/Real(10.0), tol);
}

TEST(ObsNudgingKernel, OverlappingStationsShareOneCappedRate)
{
    // Two stations at the same point: W = min(1, 2) = 1 and the target is their mean
    auto h = make_view({{Real(0.0), Real(0.0)}, {Real(0.0), Real(0.0)}},
                       {Real(50.0)}, {Real(8.0)}, {Real(0.0)});
    h.mean[1] = Real(4.0);
    h.finish(Real(100.0), Real(20.0), Real(6.0), Real(1.0), Real(10.0), Real(1.0));
    const Real phi = Real(5.0);
    EXPECT_NEAR(obs_nudging_tendency(h.v, U, Real(0.0), Real(0.0), Real(50.0), Real(0.0), phi),
                -(phi - Real(6.0))/Real(10.0), tol);
}

TEST(ObsNudgingKernel, InterpolatesAProfileAndTapersOutsideIt)
{
    auto h = make_view({{Real(0.0), Real(0.0)}}, {Real(50.0), Real(100.0), Real(150.0)},
                       {Real(5.0), Real(7.0), Real(9.0)}, {Real(0.0), Real(0.0), Real(0.0)});
    const Real rz = 20, tau = 10;
    h.finish(Real(100.0), rz, Real(100.0), Real(1.0), tau, Real(1.0));
    const Real phi = Real(0.0);
    // Inside the range: linear in height, full weight
    EXPECT_NEAR(obs_nudging_tendency(h.v, U, Real(0.0), Real(0.0), Real(75.0), Real(0.0), phi),
                -(phi - Real(6.0))/tau, tol);
    // Above it: the top value, tapered by the distance above
    const Real w_above = std::exp(Real(-0.25)*(Real(25.0)*Real(25.0))/(rz*rz));
    EXPECT_NEAR(obs_nudging_tendency(h.v, U, Real(0.0), Real(0.0), Real(175.0), Real(0.0), phi),
                -(w_above/tau)*(phi - Real(9.0)), tol);
    // Below it: the bottom value, tapered by the distance below
    const Real w_below = std::exp(Real(-0.25)*(Real(40.0)*Real(40.0))/(rz*rz));
    EXPECT_NEAR(obs_nudging_tendency(h.v, U, Real(0.0), Real(0.0), Real(10.0), Real(0.0), phi),
                -(w_below/tau)*(phi - Real(5.0)), tol);
}

TEST(ObsNudgingKernel, OnlyTheValueOutsideTheSigmaBandIsNudged)
{
    auto h = make_view({{Real(0.0), Real(0.0)}}, {Real(50.0)}, {Real(8.0)}, {Real(0.5)});
    const Real tau = 10;
    for (const Real alpha : {Real(1.0), Real(3.0)}) {
        h.finish(Real(100.0), Real(20.0), Real(100.0), Real(1.0), tau, Real(1.0));
        h.v.sigma_factor = alpha;
        const Real band = alpha * Real(0.5);
        // inside the band: nothing
        EXPECT_EQ(obs_nudging_tendency(h.v, U, 0, 0, Real(50.0), 0, Real(8.0) + Real(0.9)*band), Real(0.0));
        // below it: to the lower edge; above it: to the upper edge
        EXPECT_NEAR(obs_nudging_tendency(h.v, U, 0, 0, Real(50.0), 0, Real(5.0)),
                    -(Real(5.0) - (Real(8.0) - band))/tau, tol);
        EXPECT_NEAR(obs_nudging_tendency(h.v, U, 0, 0, Real(50.0), 0, Real(12.0)),
                    -(Real(12.0) - (Real(8.0) + band))/tau, tol);
    }
    // alpha = 0 relaxes to the mean
    h.v.sigma_factor = Real(0.0);
    EXPECT_NEAR(obs_nudging_tendency(h.v, U, 0, 0, Real(50.0), 0, Real(8.2)), -(Real(8.2) - Real(8.0))/tau, tol);
}

TEST(ObsNudgingKernel, HeightsAreAboveTheLocalTerrainUnlessAboveZeroIsAsked)
{
    auto agl = make_view({{Real(0.0), Real(0.0)}}, {Real(50.0)}, {Real(8.0)}, {Real(0.0)}, true);
    agl.finish(Real(100.0), Real(10.0), Real(6.0), Real(1.0), Real(10.0), Real(1.0));
    // 50 m above ground that is 100 m up: full weight
    EXPECT_NEAR(obs_nudging_tendency(agl.v, U, 0, 0, Real(150.0), Real(100.0), Real(5.0)),
                -(Real(5.0) - Real(8.0))/Real(10.0), tol);
    // at z = 50 the cell is inside the ground and 50 m from the height: nothing
    EXPECT_EQ(obs_nudging_tendency(agl.v, U, 0, 0, Real(50.0), Real(100.0), Real(5.0)), Real(0.0));

    auto msl = make_view({{Real(0.0), Real(0.0)}}, {Real(50.0)}, {Real(8.0)}, {Real(0.0)}, false);
    msl.finish(Real(100.0), Real(10.0), Real(6.0), Real(1.0), Real(10.0), Real(1.0));
    EXPECT_NEAR(obs_nudging_tendency(msl.v, U, 0, 0, Real(50.0), Real(100.0), Real(5.0)),
                -(Real(5.0) - Real(8.0))/Real(10.0), tol);
}

TEST(ObsNudgingKernel, TheRateNeverExceedsOneOverTheStep)
{
    auto h = make_view({{Real(0.0), Real(0.0)}}, {Real(50.0)}, {Real(8.0)}, {Real(0.0)});
    // tau = 0.2 s with dt = 1 s: the rate is 1/dt, so one step lands on the target
    h.finish(Real(100.0), Real(20.0), Real(6.0), Real(1.0), Real(0.2), Real(1.0));
    EXPECT_NEAR(obs_nudging_tendency(h.v, U, 0, 0, Real(50.0), 0, Real(5.0)), Real(3.0), tol);
}

TEST(ObsNudgingKernel, DistanceIsToTheNearestPeriodicImage)
{
    auto h = make_view({{Real(990.0), Real(500.0)}}, {Real(50.0)}, {Real(8.0)}, {Real(0.0)});
    h.finish(Real(20.0), Real(20.0), Real(6.0), Real(1.0), Real(10.0), Real(1.0));
    h.v.period_x = Real(1000.0);
    // x = 10 is 20 m from the station through the periodic boundary
    const Real w = std::exp(Real(-0.25));
    EXPECT_NEAR(obs_nudging_tendency(h.v, U, Real(10.0), Real(500.0), Real(50.0), 0, Real(5.0)),
                -(w/Real(10.0))*(Real(5.0) - Real(8.0)), tol);
    // the two copies of the face at x = 0 and x = 1000 see the same distance
    EXPECT_EQ(obs_nudging_tendency(h.v, U, Real(0.0),    Real(500.0), Real(50.0), 0, Real(5.0)),
              obs_nudging_tendency(h.v, U, Real(1000.0), Real(500.0), Real(50.0), 0, Real(5.0)));
}

TEST(ObsNudgingKernel, StationsBeyondTheCutoffAreSkipped)
{
    auto h = make_view({{Real(0.0), Real(0.0)}}, {Real(50.0)}, {Real(8.0)}, {Real(0.0)});
    h.finish(Real(100.0), Real(20.0), Real(2.0), Real(1.0), Real(10.0), Real(1.0));
    EXPECT_NE(obs_nudging_tendency(h.v, U, Real(190.0), 0, Real(50.0), 0, Real(5.0)), Real(0.0));
    EXPECT_EQ(obs_nudging_tendency(h.v, U, Real(210.0), 0, Real(50.0), 0, Real(5.0)), Real(0.0));
}

// ---------------------------------------------------------------------------
// The sources on the staggered grid
// ---------------------------------------------------------------------------

namespace {

struct Grid
{
    amrex::Box domain{amrex::IntVect(0,0,0), amrex::IntVect(7,7,7)};
    amrex::RealBox rb{{0.0, 0.0, 0.0}, {800.0, 800.0, 400.0}};
    amrex::Array<int, AMREX_SPACEDIM> periodic{1, 1, 0};
    amrex::Geometry geom{domain, &rb, amrex::CoordSys::cartesian, periodic.data()};
    amrex::BoxArray ba{domain};
    amrex::DistributionMapping dm{ba};
    amrex::MFInfo info = amrex::MFInfo().SetArena(amrex::The_Pinned_Arena());

    amrex::MultiFab cons{ba, dm, NVAR_max, 1, info};
    amrex::MultiFab u{amrex::convert(ba, amrex::IntVect(1,0,0)), dm, 1, 1, info};
    amrex::MultiFab v{amrex::convert(ba, amrex::IntVect(0,1,0)), dm, 1, 1, info};
    amrex::MultiFab w{amrex::convert(ba, amrex::IntVect(0,0,1)), dm, 1, 1, info};
    amrex::MultiFab su{amrex::convert(ba, amrex::IntVect(1,0,0)), dm, 1, 0, info};
    amrex::MultiFab sv{amrex::convert(ba, amrex::IntVect(0,1,0)), dm, 1, 0, info};
    amrex::MultiFab sw{amrex::convert(ba, amrex::IntVect(0,0,1)), dm, 1, 0, info};
    amrex::MultiFab scc{ba, dm, NVAR_max, 0, info};

    Grid ()
    {
        cons.setVal(Real(0.0));
        cons.setVal(Real(2.0), Rho_comp, 1, 1);          // rho = 2
        cons.setVal(Real(600.0), RhoTheta_comp, 1, 1);   // theta = 300
        u.setVal(Real(0.0)); v.setVal(Real(0.0)); w.setVal(Real(0.0));
        su.setVal(Real(0.0)); sv.setVal(Real(0.0)); sw.setVal(Real(0.0)); scc.setVal(Real(0.0));
        amrex::Gpu::streamSynchronize();
    }
};

const char* profile_file =
    "time z u v w theta\n"
    "0 0   2 -1 0.5 301\n"
    "0 400 2 -1 0.5 301\n";

} // namespace

TEST(ObsNudgingSources, LandOnEveryFaceOverFlatGround)
{
    write_file("erf_unit_obs_nudging_flat.txt", profile_file);
    set_inputs("erf_unit_obs_nudging_flat.txt", Real(1.0e8), Real(25.0), Real(10.0));

    ObsNudging nudging(TerrainType::None, nullptr, false, 0.0);
    Grid g;
    nudging.resolve_positions(g.geom, nullptr);
    nudging.add_momentum_sources(0, 0.0, Real(1.0), g.geom, g.cons, g.u, g.v, g.w,
                                 g.su, g.sv, g.sw, nullptr, nullptr, nullptr, nullptr);
    nudging.add_theta_source(0, 0.0, Real(1.0), g.geom, g.cons, g.scc, nullptr, nullptr);
    amrex::Gpu::streamSynchronize();

    // rho (target - 0) / tau on every face, with the station's weight 1 - 1e-12
    const Real ftol = Real(1.0e-6);
    const auto a_u = g.su.const_array(0);
    const auto a_v = g.sv.const_array(0);
    const auto a_w = g.sw.const_array(0);
    const auto a_t = g.scc.const_array(0);
    for (int k = 0; k <= 7; ++k) {
        for (int j = 0; j <= 7; ++j) {
            for (int i = 0; i <= 7; ++i) {
                EXPECT_NEAR(a_u(i,j,k), Real(2.0)*Real(2.0)/Real(10.0), ftol);
                EXPECT_NEAR(a_v(i,j,k), Real(2.0)*Real(-1.0)/Real(10.0), ftol);
                EXPECT_NEAR(a_t(i,j,k,RhoTheta_comp), Real(2.0)*Real(1.0)/Real(10.0), ftol);
                EXPECT_EQ(a_t(i,j,k,Rho_comp), Real(0.0));
            }
        }
        for (int j = 0; j <= 7; ++j) {
            for (int i = 0; i <= 7; ++i) {
                EXPECT_NEAR(a_w(i,j,k), (k == 0) ? Real(0.0) : Real(2.0)*Real(0.5)/Real(10.0), ftol);
            }
        }
    }
    // w on the top of the domain is set by the boundary condition
    EXPECT_EQ(a_w(3,3,8), Real(0.0));
}

TEST(ObsNudgingSources, SitOnTheStaggeredPositions)
{
    // A station 100 m wide at (400, 400): each component's weight depends on
    // where its face is, u at (i dx, (j+1/2) dy), v at ((i+1/2) dx, j dy), w and
    // theta at the cell centre column
    write_file("erf_unit_obs_nudging_faces.txt", profile_file);
    set_inputs("erf_unit_obs_nudging_faces.txt", Real(100.0), Real(25.0), Real(10.0));

    ObsNudging nudging(TerrainType::None, nullptr, false, 0.0);
    Grid g;
    nudging.resolve_positions(g.geom, nullptr);
    nudging.add_momentum_sources(0, 0.0, Real(1.0), g.geom, g.cons, g.u, g.v, g.w,
                                 g.su, g.sv, g.sw, nullptr, nullptr, nullptr, nullptr);
    nudging.add_theta_source(0, 0.0, Real(1.0), g.geom, g.cons, g.scc, nullptr, nullptr);
    amrex::Gpu::streamSynchronize();

    auto w = [] (Real dx, Real dy) {
        return std::exp(Real(-0.25)*(dx*dx + dy*dy)/Real(1.0e4));
    };
    const Real ftol = Real(1.0e-6);
    // u face (4,3): x = 400, y = 350;  (5,3): x = 500, y = 350
    EXPECT_NEAR(g.su.const_array(0)(4,3,2), Real(0.4)*w(Real(0.0),  Real(50.0)), ftol);
    EXPECT_NEAR(g.su.const_array(0)(5,3,2), Real(0.4)*w(Real(100.0),Real(50.0)), ftol);
    // v face (3,4): x = 350, y = 400
    EXPECT_NEAR(g.sv.const_array(0)(3,4,2), Real(-0.2)*w(Real(50.0), Real(0.0)), ftol);
    // w face and cell (3,3): x = 350, y = 350
    EXPECT_NEAR(g.sw.const_array(0)(3,3,2), Real(0.1)*w(Real(50.0), Real(50.0)), ftol);
    EXPECT_NEAR(g.scc.const_array(0)(3,3,2,RhoTheta_comp), Real(0.2)*w(Real(50.0), Real(50.0)), ftol);
}

TEST(ObsNudgingSources, MeasureHeightsFromTheBottomOfAFittedMesh)
{
    // A profile only between 0 and 50 m above the ground, which slopes up in x
    // (h = 20 + 0.1 x), and a vertical radius of 5 m: cells whose centre is
    // more than a few radii above 50 m above the ground are not nudged
    write_file("erf_unit_obs_nudging_fitted.txt",
               "time z u v\n0 0 2 0\n0 50 2 0\n");
    set_inputs("erf_unit_obs_nudging_fitted.txt", Real(1.0e8), Real(5.0), Real(10.0));

    ObsNudging nudging(TerrainType::StaticFittedMesh, nullptr, false, 0.0);
    Grid g;
    amrex::MultiFab znd(amrex::convert(g.ba, amrex::IntVect(1,1,1)), g.dm, 1, 1, g.info);
    fill_fitted_mesh(znd, g.geom, Real(20.0), Real(0.1));
    amrex::Gpu::streamSynchronize();

    nudging.resolve_positions(g.geom, nullptr);
    nudging.add_momentum_sources(0, 0.0, Real(1.0), g.geom, g.cons, g.u, g.v, g.w,
                                 g.su, g.sv, g.sw, &znd, nullptr, nullptr, nullptr);
    amrex::Gpu::streamSynchronize();

    const auto a_u = g.su.const_array(0);
    const Real full = Real(2.0)*Real(2.0)/Real(10.0);
    const Real dz = Real(50.0);
    for (int i = 0; i <= 7; ++i) {
        const Real x = Real(100.0)*i;
        const Real ground = Real(20.0) + Real(0.1)*x;
        for (int k = 0; k <= 7; ++k) {
            // The x-face centre is at ground + (k + 1/2) dz above the ground
            const Real h = (Real(k) + Real(0.5))*dz;
            const Real d = amrex::max(h - Real(50.0), Real(0.0));
            const Real q = d*d/Real(25.0);
            const Real expected = (q > Real(1.0e6)) ? Real(0.0) : full*std::exp(Real(-0.25)*q);
            EXPECT_NEAR(a_u(i,4,k), expected, Real(1.0e-6)) << "i " << i << " k " << k << " ground " << ground;
        }
    }
}

TEST(ObsNudgingSources, SkipCellsInsideImmersedTerrain)
{
    write_file("erf_unit_obs_nudging_blank.txt", profile_file);
    set_inputs("erf_unit_obs_nudging_blank.txt", Real(1.0e8), Real(25.0), Real(10.0));

    ObsNudging nudging(TerrainType::None, nullptr, false, 0.0);
    Grid g;
    amrex::MultiFab bcc(g.ba, g.dm, 1, 1, g.info);
    amrex::MultiFab bx(amrex::convert(g.ba, amrex::IntVect(1,0,0)), g.dm, 1, 1, g.info);
    amrex::MultiFab by(amrex::convert(g.ba, amrex::IntVect(0,1,0)), g.dm, 1, 1, g.info);
    amrex::MultiFab bz(amrex::convert(g.ba, amrex::IntVect(0,0,1)), g.dm, 1, 1, g.info);
    fill_blank(bcc, 2); fill_blank(bx, 2); fill_blank(by, 2); fill_blank(bz, 2);
    amrex::Gpu::streamSynchronize();

    nudging.resolve_positions(g.geom, nullptr);
    nudging.add_momentum_sources(0, 0.0, Real(1.0), g.geom, g.cons, g.u, g.v, g.w,
                                 g.su, g.sv, g.sw, nullptr, &bx, &by, &bz);
    nudging.add_theta_source(0, 0.0, Real(1.0), g.geom, g.cons, g.scc, nullptr, &bcc);
    amrex::Gpu::streamSynchronize();

    const auto a_u = g.su.const_array(0);
    const auto a_t = g.scc.const_array(0);
    const Real full = Real(2.0)*Real(2.0)/Real(10.0);
    EXPECT_EQ(a_u(3,3,0), Real(0.0));
    EXPECT_EQ(a_u(3,3,1), Real(0.0));
    EXPECT_NEAR(a_u(3,3,2), Real(0.5)*full, Real(1.0e-6));
    EXPECT_NEAR(a_u(3,3,3), full, Real(1.0e-6));
    EXPECT_EQ(a_t(3,3,1,RhoTheta_comp), Real(0.0));
    EXPECT_NEAR(a_t(3,3,3,RhoTheta_comp), Real(2.0)*Real(1.0)/Real(10.0), Real(1.0e-6));
}

TEST(ObsNudgingSources, FollowTheFileInTime)
{
    // Elapsed time: the target at t = 30 is a third of the way from 2 to 5
    write_file("erf_unit_obs_nudging_time.txt",
               "time z u v\n0 0 2 0\n0 400 2 0\n90 0 5 0\n90 400 5 0\n");
    set_inputs("erf_unit_obs_nudging_time.txt", Real(1.0e8), Real(25.0), Real(10.0));

    ObsNudging nudging(TerrainType::None, nullptr, false, 0.0);
    Grid g;
    nudging.resolve_positions(g.geom, nullptr);
    nudging.add_momentum_sources(0, 30.0, Real(1.0), g.geom, g.cons, g.u, g.v, g.w,
                                 g.su, g.sv, g.sw, nullptr, nullptr, nullptr, nullptr);
    amrex::Gpu::streamSynchronize();
    EXPECT_NEAR(g.su.const_array(0)(3,3,3), Real(2.0)*Real(3.0)/Real(10.0), Real(1.0e-6));

    // After the last record the station is inactive
    g.su.setVal(Real(0.0));
    amrex::Gpu::streamSynchronize();
    nudging.add_momentum_sources(0, 91.0, Real(1.0), g.geom, g.cons, g.u, g.v, g.w,
                                 g.su, g.sv, g.sw, nullptr, nullptr, nullptr, nullptr);
    amrex::Gpu::streamSynchronize();
    EXPECT_EQ(g.su.const_array(0)(3,3,3), Real(0.0));
}

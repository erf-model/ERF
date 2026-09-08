#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "ERF_GTestTerminalVelocityCommon.H"

using namespace terminal_velocity_test;

namespace {

// Log-spaced sweep helper.
amrex::Real log_point (const amrex::Real x0, const amrex::Real x1, const int k, const int n)
{
    const amrex::Real t = static_cast<amrex::Real>(k) / static_cast<amrex::Real>(n - 1);
    return std::exp(std::log(x0) + t * (std::log(x1) - std::log(x0)));
}

} // namespace

// Motivation: The hydrometeor fall-speed routines feed the SDM sedimentation
// step. This host test proves, property by property, that the Beard liquid-drop
// speed with the Gunn-Kinzer 4 mm clamp, the Boehm ice speed (cross-checked
// against Boehm 1989 Fig. 3 magnitudes), and the Frick melting-particle blend
// all behave correctly. Each metric is built to return zero when its property
// holds; see ERF_GTestTerminalVelocityCommon.H.
TEST(TerminalVelocityScalar, FallSpeedPropertiesHoldOnHost)
{
    for (int i = 0; i < kNumMetrics; ++i) {
        SCOPED_TRACE(metric_label(i));
        EXPECT_LE(terminal_velocity_metric(i), kTol);
    }
}

// Motivation: Emit the fall-speed curves evaluated by the production routines so
// they can be plotted against published references (Boehm 1989 Fig. 3, Mitra
// 1990 Fig. 2). This is a data dump, not an assertion: it is skipped unless the
// environment variable ERF_TV_CURVES names an output CSV, so it is a no-op in CI
// and is driven on demand by the 0D_ice_terminal_velocity run script. Long
// format: series, x, y (x = max dimension [mm] for fall-speed-vs-size series,
// x = liquid mass fraction for the melting-blend series; y = fall speed [m/s],
// except series "blend_psi" whose y is the dimensionless weight Psi).
TEST(TerminalVelocityScalar, DumpCurvesForPlotting)
{
    const char* out = std::getenv("ERF_TV_CURVES");
    if (out == nullptr || out[0] == '\0') {
        GTEST_SKIP() << "set ERF_TV_CURVES=<file.csv> to dump the fall-speed curves";
    }

    const TerminalVelocity<amrex::Real> tv{kRhoWater, kRhoIce};

    std::ofstream f(out);
    ASSERT_TRUE(f.is_open()) << "cannot open " << out;
    f << "series,x,y\n";
    f.precision(8);

    const int n = 80;

    // Liquid drops: fall speed vs diameter, out past the 4 mm-radius clamp.
    for (int k = 0; k < n; ++k) {
        const amrex::Real R = log_point(amrex::Real(2.0e-6), amrex::Real(1.0e-2), k, n);
        f << "rain," << amrex::Real(2.0e3) * R << ',' << rain_vterm(tv, R) << '\n';
    }

    // Ice habits: fall speed vs maximum dimension (rho_a = 1.1, the Boehm Fig. 3
    // condition). phi = c/a, rho_app = apparent density [kg m^-3].
    struct Habit { const char* name; amrex::Real phi; amrex::Real rho; amrex::Real Dmax; };
    const Habit habits[] = {
        {"ice_hail",     amrex::Real(1.0),  amrex::Real(916.8), amrex::Real(3.0e-2)},
        {"ice_graupel",  amrex::Real(1.0),  amrex::Real(400.0), amrex::Real(1.0e-2)},
        {"ice_plate",    amrex::Real(0.1),  amrex::Real(900.0), amrex::Real(1.0e-2)},
        {"ice_column",   amrex::Real(8.0),  amrex::Real(700.0), amrex::Real(1.0e-2)},
        {"ice_snow_agg", amrex::Real(1.0),  amrex::Real(50.0),  amrex::Real(2.0e-2)}
    };
    for (const Habit& h : habits) {
        for (int k = 0; k < n; ++k) {
            const amrex::Real D = log_point(amrex::Real(1.0e-4), h.Dmax, k, n);
            f << h.name << ',' << amrex::Real(1.0e3) * D << ','
              << ice_vterm_D(tv, D, h.phi, h.rho) << '\n';
        }
    }

    // Locatelli & Hobbs (1974) graupel habits at their measured mass-dimension
    // relations (M[mg] = c D[mm]^b), the basis for the Boehm (1989) Fig. 3
    // graupel curves. Apparent density rho = M / ((pi/6) D^3) falls with D; the
    // particle is treated as a sphere (phi = 1). This overlays ERF's Boehm speed
    // on the Fig. 3 read-offs at matching density (instead of a fixed 400 sphere).
    struct LHHabit { const char* name; amrex::Real c_mg; amrex::Real b; };
    const LHHabit lh_habits[] = {
        {"lh_lump_graupel",      amrex::Real(0.078), amrex::Real(2.8)}, // L-H Fig.6 (0.1-0.2 Mg m^-3)
        {"lh_conical_graupel",   amrex::Real(0.073), amrex::Real(2.6)}, // L-H Fig.9
        {"lh_hexagonal_graupel", amrex::Real(0.044), amrex::Real(2.9)}  // L-H Fig.9
    };
    for (const LHHabit& h : lh_habits) {
        for (int k = 0; k < n; ++k) {
            const amrex::Real D = log_point(amrex::Real(3.0e-4), amrex::Real(5.0e-3), k, n);
            const amrex::Real D_mm = amrex::Real(1.0e3) * D;
            const amrex::Real m = h.c_mg * amrex::Real(1.0e-6) * std::pow(D_mm, h.b); // [kg]
            const amrex::Real rho_app = m / ((PI / amrex::Real(6.0)) * D*D*D);
            f << h.name << ',' << D_mm << ','
              << ice_vterm_D(tv, D, amrex::Real(1.0), rho_app) << '\n';
        }
    }

    // A fixed-density (100 kg m^-3) ice sphere swept from the Stokes regime upward.
    // Plotted with the analytic Stokes law v = (2/9)(rho_p - rho_a) g a^2 / mu so the
    // figure shows the small-ice limit (Re << 1) where Boehm reduces to Stokes -- the
    // viscosity-controlled regime exercised by the Stokes assertion.
    for (int k = 0; k < n; ++k) {
        const amrex::Real D = log_point(amrex::Real(2.0e-5), amrex::Real(5.0e-3), k, n);
        f << "ice_sphere100," << amrex::Real(1.0e3) * D << ','
          << ice_vterm_D(tv, D, amrex::Real(1.0), amrex::Real(100.0)) << '\n';
    }

    // Melting graupel: a fixed-total-mass particle (equivalent radius 1.5 mm,
    // apparent ice-core density 400 kg m^-3) sweeping liquid mass fraction ell.
    // v_mixed blends the shrinking ice core's Boehm speed toward the equivalent
    // raindrop speed; blend_psi is the Frick (2013) weight itself.
    const amrex::Real R_eq = amrex::Real(1.5e-3);
    const amrex::Real M    = spheroid_mass(R_eq, R_eq, kRhoWater); // (4/3) pi rho_w R^3
    const amrex::Real v_rain = rain_vterm(tv, R_eq);
    for (int k = 0; k < n; ++k) {
        const amrex::Real ell = static_cast<amrex::Real>(k) / static_cast<amrex::Real>(n - 1);
        const amrex::Real m_ice = (amrex::Real(1.0) - ell) * M;
        amrex::Real v_ice = amrex::Real(0.0);
        if (m_ice > amrex::Real(0.0)) {
            v_ice = ice_vterm_m(tv, m_ice, amrex::Real(1.0), amrex::Real(400.0));
        }
        f << "blend_vterm," << ell << ',' << tv.MixedPhaseVterm(v_ice, v_rain, ell) << '\n';
        f << "blend_psi,"   << ell << ',' << tv.MeltBlendWeight(ell) << '\n';
    }

    f.close();
    EXPECT_TRUE(f.good());
}

// Motivation: evaluate the production phase-routed fall speed for arbitrary
// (m_ice, m_water, a, c, rho_air, p, T) tuples, replicating the SD advection step
// (ERF_SuperDropletPCAdvection.cpp): water -> CloudRainShima; pure ice -> IceBohm
// on the ice core; mixed -> Frick blend of the two. Used to Monte-Carlo which
// attribute/flow combinations produce anomalously large terminal velocities.
// Reads the input CSV named by ERF_TV_EVAL_IN (columns m_ice,m_water,a,c,rho,p,T)
// and writes m_ice,...,T,vterm,phase to ERF_TV_EVAL_OUT; skipped otherwise.
TEST(TerminalVelocityScalar, EvalFromCsv)
{
    const char* in  = std::getenv("ERF_TV_EVAL_IN");
    const char* out = std::getenv("ERF_TV_EVAL_OUT");
    if (in == nullptr || out == nullptr || in[0] == '\0' || out[0] == '\0') {
        GTEST_SKIP() << "set ERF_TV_EVAL_IN and ERF_TV_EVAL_OUT to evaluate a CSV";
    }
    std::ifstream fi(in);
    std::ofstream ofs(out);
    ASSERT_TRUE(fi.is_open() && ofs.is_open());

    const TerminalVelocity<amrex::Real> tv{kRhoWater, kRhoIce};
    const amrex::Real f43pi = four_thirds_pi;

    std::string line;
    std::getline(fi, line);  // header
    ofs << "m_ice,m_water,a,c,rho,p,T,vterm,phase\n";
    ofs.precision(8);
    while (std::getline(fi, line)) {
        if (line.empty()) { continue; }
        amrex::Real m_ice, m_water, a, c, rho, p, T;
        char comma;
        std::istringstream ss(line);
        ss >> m_ice >> comma >> m_water >> comma >> a >> comma >> c >> comma
           >> rho >> comma >> p >> comma >> T;

        amrex::Real vterm; int phase;
        if (m_ice <= amrex::Real(0.0)) {              // water
            const amrex::Real r = std::cbrt(m_water / (f43pi * kRhoWater));
            vterm = tv.CloudRainShima(r, rho, p, T);
            phase = 0;
        } else {
            const amrex::Real rho_app = m_ice / (f43pi * a * a * c);
            const amrex::Real v_ice = tv.IceBohm(m_ice, a, c, rho_app, rho, T);
            if (m_water > amrex::Real(0.0)) {          // mixed (Frick blend)
                const amrex::Real m_w = m_ice + m_water;
                const amrex::Real r_eq = std::cbrt(m_w / (f43pi * kRhoWater));
                const amrex::Real v_rain = tv.CloudRainShima(r_eq, rho, p, T);
                vterm = tv.MixedPhaseVterm(v_ice, v_rain, m_water / m_w);
                phase = 2;
            } else {                                   // pure ice
                vterm = v_ice;
                phase = 1;
            }
        }
        ofs << m_ice << ',' << m_water << ',' << a << ',' << c << ',' << rho << ','
           << p << ',' << T << ',' << vterm << ',' << phase << '\n';
    }
    ofs.close();
    EXPECT_TRUE(ofs.good());
}

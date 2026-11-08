#include <gtest/gtest.h>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FileSystem.H>
#include <AMReX_ParallelDescriptor.H>
#include "../ERF_GTestTempDir.H"
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>
#include <ERF_EOS.H>
#include "ERF_FireCoupler.H"
#include "ERF_FireEngine.H"
#include "ERF_FireGrid.H"
#include "ERF_FireCouplingParams.H"

/**
 * @file ERF_GTestFireCoupler.cpp
 * @brief The shared fire coupling with a stand-in engine
 *
 * A stand-in engine plugs into FireCoupler through the FireEngine interface
 * alone and reports a sensible flux linear in x, a uniform latent flux and a
 * heat per kilogram of its own. Every expected value is a closed form, taken
 * from the physics and never from the coupling's code:
 * - the flux averages onto the columns it came from: a linear field's mean
 *   over a column's C x C fire cells is its value at the column centre, also
 *   for a fire region that does not start at the domain corner;
 * - the engine gets the wind before the dycore with lagged coupling and after
 *   it with synchronous coupling, and the rain only when it asks for it;
 * - the injected heat integrates to the surface flux: with the density factor
 *   on and feedback f, c_p sum_k d(rho theta)/dt dz = rho f Q (1 - e^{-z_top/alfg}),
 *   z measured from the first cell centre, and L_v sum_k d(rho q_v)/dt dz the
 *   same for the latent flux; the smoke is e Q / h in the lowest cell;
 * - source_mode = add keeps the other sources, overwrite replaces the theta
 *   and vapour sources only;
 * - passive coupling injects no heat (and, without smoke, asks for no flux),
 *   and inject_latent = false drops a latent flux the engine reports;
 * - with erf.fire.heat_open_fraction and an engine that reports an open
 *   fraction f and a roof height H above the ground, the column still receives
 *   f Q (1 - e^{-z_top/alfg}), and every cell below the roof receives
 *   f L / (L - (1 - f) B) times what it receives without the option, every cell
 *   above it L / (L - (1 - f) B) times (L = 1 - e^{-z_top/alfg}, B the plain
 *   profile's share in the cells under the roof, 1 - e^{-2 dz/alfg} for a roof
 *   two cells high, the profile measured from the first cell centre);
 * - the erf.fire.* coupling keys reach FireCouplingParams unchanged;
 * - a checkpoint restores the lagged fluxes bit for bit, ghosts refilled
 *   (also from a checkpoint written without ghost columns), after the
 *   engine's own restore; a checkpoint without them leaves them at zero, and
 *   one of another level or region is refused.
 */

using namespace amrex;

namespace {

constexpr int  NX = 6, NY = 4, NZ = 50, C = 2;
constexpr Real DXY = 100.0, DZ = 10.0;
constexpr Real Q0 = 2.0e4, QX = 30.0;      // sensible flux a + b x [W/m2], x in m
constexpr Real QLAT = 3.0e3;               // latent flux [W/m2]
constexpr Real HEAT_PER_KG = 1.5e7;        // the engine's heat per kg fuel [J/kg]
constexpr Real ALFG = 100.0;               // e-folding height [m]: 1 - e^{-5} of the flux stays in the 500 m column
constexpr Real RHO = 1.15;                 // air density [kg/m3]
constexpr Real FEEDBACK = 0.5;             // erf.fire.fire_atm_feedback
constexpr Real EMISSION = 0.02;            // erf.fire.smoke_emission_factor
constexpr Real BACKGROUND = 7.0;           // every other source of the cell
constexpr Real Cp_fire = Cp_d, L_V = L_v;   // ERF dry-air c_p and latent heat, as the coupling uses
/// Relative round-off tolerance of the build precision
constexpr double REL = (sizeof(Real) == 8) ? 1.0e-9 : 1.0e-4;

// nvcc refuses an extended __device__ lambda in a constructor or in a function
// it cannot take the address of, and in gtest's private TestBody: every kernel
// is in a namespace-scope function.

/// q = a + b x at the fire cell centres, x the physical position
void fill_linear_flux (MultiFab& q, const Geometry& g)
{
    const Real dxf = g.CellSize(0);
    const Real x0  = g.ProbLo(0);
    for (MFIter mfi(q); mfi.isValid(); ++mfi) {
        auto a = q.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            a(i, j, k) = Q0 + QX * (x0 + (i + 0.5_rt) * dxf);
        });
    }
}

/// Cell-centre heights of a flat mesh, ghosts included
void fill_zcc (MultiFab& zcc)
{
    for (MFIter mfi(zcc); mfi.isValid(); ++mfi) {
        auto z = zcc.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { z(i, j, k) = (k + 0.5_rt) * DZ; });
    }
}

/// A stand-in engine: fixed flux fields, and a record of what it was handed.
class StandInEngine : public FireEngine
{
public:
    const char* name () const override { return "stand-in"; }
    void verify_prerequisites (const ERF&, int) override {}
    void initialize (const ERF&, int, const FireGrid&, const SurfaceLayer*, const MultiFab*) override {}

    /// Allocate and fill the flux fields on the coupler's fire grid (the test calls it).
    void set_up (const FireGrid& fg)
    {
        sensible = std::make_unique<MultiFab>(fg.ba, fg.dm, 1, 0);
        fill_linear_flux(*sensible, fg.geom);
        latent = std::make_unique<MultiFab>(fg.ba, fg.dm, 1, 0);
        latent->setVal(QLAT);
    }

    void advance (const FireForcing& forcing) override
    {
        ++n_advance;
        last_xvel   = &forcing.xvel;
        last_yvel   = &forcing.yvel;
        // the rain field lives only for the call: record what the test checks
        got_precip  = (forcing.precip_accum_k0 != nullptr);
        if (got_precip) {
            precip_ba    = forcing.precip_accum_k0->boxArray();
            precip_ncomp = forcing.precip_accum_k0->nComp();
        }
        last_time   = forcing.time;
        last_dt     = forcing.dt;
    }
    void fill_fluxes (FireFluxes& fluxes) override
    {
        ++n_fill;
        fluxes.sensible = report_sensible ? sensible.get() : nullptr;
        fluxes.latent   = report_latent ? latent.get() : nullptr;
    }
    Real smoke_heat_per_kg () const override { return HEAT_PER_KG; }
    bool wants_surface_precip () const override { return wants_precip; }
    const MultiFab* heat_open_fraction_atm () const override { return open_fraction.get(); }
    const MultiFab* roof_height_atm () const override { return roof_height.get(); }
    void start_run (Real, int, int) override {}
    void write_checkpoint (const std::string&) const override { ++n_write_chk; }
    void read_checkpoint (const std::string&, int, int step) override
    {
        ++n_read_chk; read_step = step;
        // the coupler restores its fluxes after the engine: still zero here
        if (owner) { flux_at_read = owner->get_Q_atm_prev()->norm0(0, 1, IntVect(5, 5, 0)); }
    }
    void write_plotfile (const std::string&, Real, int) const override {}

    std::unique_ptr<MultiFab> sensible, latent;
    std::unique_ptr<MultiFab> open_fraction, roof_height;   ///< on the level's k = 0 slabs, when the test sets them
    bool report_latent = true;
    bool report_sensible = true;
    bool wants_precip  = false;
    int n_advance = 0, n_fill = 0;
    mutable int n_write_chk = 0;
    int n_read_chk = 0, read_step = -1;
    const FireCoupler* owner = nullptr;   ///< the coupler holding this engine, for the order check
    Real flux_at_read = -1.0;
    const MultiFab* last_xvel = nullptr;
    const MultiFab* last_yvel = nullptr;
    bool got_precip = false;
    BoxArray precip_ba;
    int precip_ncomp = 0;
    Real last_time = -1.0, last_dt = -1.0;
};

/// The level the fire grid refines: the whole NX x NY x NZ domain in two
/// boxes of unequal width (periodic in x and y), or the region of columns 2-5
/// in two boxes, so the fire grid starts at column 2 (atm_lo = (2, 0)).
struct Level
{
    BoxArray ba, ba2d;
    DistributionMapping dm;
    Geometry geom;
    MultiFab S, xold, yold, xnew, ynew, zcc;
};

std::unique_ptr<Level> make_level (bool region)
{
    auto L = std::make_unique<Level>();
    const Box domain(IntVect(0, 0, 0), IntVect(NX - 1, NY - 1, NZ - 1));
    BoxList bl;
    if (region) {
        bl.push_back(Box(IntVect(2, 0, 0), IntVect(3, NY - 1, NZ - 1)));
        bl.push_back(Box(IntVect(4, 0, 0), IntVect(5, NY - 1, NZ - 1)));
    } else {
        bl.push_back(Box(IntVect(0, 0, 0), IntVect(3, NY - 1, NZ - 1)));
        bl.push_back(Box(IntVect(4, 0, 0), IntVect(NX - 1, NY - 1, NZ - 1)));
    }
    L->ba = BoxArray(bl);
    L->dm = DistributionMapping(L->ba);
    L->geom = Geometry(domain, RealBox(0.0, 0.0, 0.0, NX * DXY, NY * DXY, NZ * DZ), CoordSys::cartesian, {1, 1, 0});
    L->ba2d = L->ba;
    L->ba2d.coarsen(IntVect(1, 1, NZ));
    L->S.define(L->ba, L->dm, NVAR_max, 1);
    L->S.setVal(0.0_rt);
    L->S.setVal(RHO, Rho_comp, 1, 1);
    L->S.setVal(getRhoThetagivenP(9.0e4_rt), RhoTheta_comp, 1, 1);
    L->xold.define(convert(L->ba, IntVect(1, 0, 0)), L->dm, 1, 1); L->xold.setVal(1.0_rt);
    L->yold.define(convert(L->ba, IntVect(0, 1, 0)), L->dm, 1, 1); L->yold.setVal(1.0_rt);
    L->xnew.define(convert(L->ba, IntVect(1, 0, 0)), L->dm, 1, 1); L->xnew.setVal(2.0_rt);
    L->ynew.define(convert(L->ba, IntVect(0, 1, 0)), L->dm, 1, 1); L->ynew.setVal(2.0_rt);
    L->zcc.define(L->ba, L->dm, 1, 1);
    fill_zcc(L->zcc);
    return L;
}

FireCouplingParams coupling (const std::string& type, bool smoke)
{
    FireCouplingParams p;
    p.grid_ratio            = C;
    p.coupling_type         = type;
    p.heat_flux_alfg        = ALFG;
    p.fire_atm_feedback     = FEEDBACK;
    p.inject_latent         = true;
    p.source_mode           = "add";
    p.heat_tendency_density = true;
    p.heat_tendency_exner   = false;
    p.smoke_enable          = smoke;
    p.smoke_emission_factor = EMISSION;
    return p;
}

/// The coupler with the stand-in engine, laid out on the level's fire grid.
std::unique_ptr<FireCoupler> make_coupler (const Level& lev, const FireCouplingParams& p, StandInEngine*& engine)
{
    auto eng = std::make_unique<StandInEngine>();
    engine = eng.get();
    auto cpl = std::make_unique<FireCoupler>(p, std::move(eng));
    FireGrid fg = create_fire_grid(lev.ba, lev.dm, lev.geom, C);
    cpl->set_fire_grid(0, fg, lev.ba, lev.dm, lev.geom);
    engine->set_up(cpl->fire_grid());
    return cpl;
}

void step (FireCoupler& cpl, Level& lev)
{
    cpl.advance(10.0_rt, 0.5_rt, lev.S, lev.xold, lev.yold, lev.xnew, lev.ynew, lev.zcc,
                lev.ba2d, lev.geom, SurfacePrecipAccumulationSources{});
}

/// The cell sources after add_sources on a uniform background
std::unique_ptr<MultiFab> sources (FireCoupler& cpl, Level& lev, Real background = BACKGROUND,
                                   const MultiFab* znd = nullptr)
{
    auto src = std::make_unique<MultiFab>(lev.ba, lev.dm, NVAR_max, 0);
    src->setVal(background);
    cpl.add_sources(*src, lev.S, lev.zcc, znd, lev.geom, true, 1);
    return src;
}

/// Column sums sum_k (src - b) dz of one component, one value per (i, j), checked by f(i, sum)
template <typename F>
void for_each_column (const MultiFab& src, int comp, Real b, F&& f)
{
    for (MFIter mfi(src); mfi.isValid(); ++mfi) {
        const auto s = src.const_array(mfi);
        const Box& bx = mfi.validbox();
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                Real sum = 0.0;
                for (int k = 0; k < NZ; ++k) { sum += (s(i, j, k, comp) - b) * DZ; }
                f(i, j, sum, s);
            }
        }
    }
}

Real flux_at_column (int i) { return Q0 + QX * (i + 0.5_rt) * DXY; }

void FluxesAverageOntoTheColumnsTheyCameFrom ()
{
    for (bool region : {false, true}) {
        auto lev = make_level(region);
        StandInEngine* eng = nullptr;
        auto cpl = make_coupler(*lev, coupling("lagged", false), eng);
        EXPECT_EQ(cpl->fire_grid().atm_lo[0], region ? 2 : 0);
        step(*cpl, *lev);

        EXPECT_EQ(eng->n_advance, 1);
        EXPECT_EQ(eng->n_fill, 1);
        EXPECT_EQ(eng->last_xvel, &lev->xold) << "lagged coupling hands over the wind before the dycore";
        EXPECT_EQ(eng->last_yvel, &lev->yold);
        EXPECT_FALSE(eng->got_precip) << "an engine that did not ask gets no rain";
        EXPECT_EQ(eng->last_time, 10.0_rt);
        EXPECT_EQ(eng->last_dt, 0.5_rt);

        // the mean of a + b x over a column's C x C fire cells is a + b x_centre;
        // on the whole periodic domain the 5 ghost columns carry the periodic
        // neighbour's value
        const MultiFab& Q  = *cpl->get_Q_atm_prev();
        const MultiFab& QL = *cpl->get_Q_lat_atm_prev();
        const IntVect ng = region ? IntVect(0, 0, 0) : IntVect(5, 5, 0);
        for (MFIter mfi(Q); mfi.isValid(); ++mfi) {
            const auto q  = Q.const_array(mfi);
            const auto ql = QL.const_array(mfi);
            const Box gbx = grow(mfi.validbox(), ng);
            for (int j = gbx.smallEnd(1); j <= gbx.bigEnd(1); ++j) {
                for (int i = gbx.smallEnd(0); i <= gbx.bigEnd(0); ++i) {
                    const int ip = ((i % NX) + NX) % NX;
                    EXPECT_NEAR(q(i, j, 0), flux_at_column(ip), REL * flux_at_column(ip))
                        << (region ? "region" : "domain") << ", column (" << i << ", " << j << ")";
                    EXPECT_NEAR(ql(i, j, 0), QLAT, REL * QLAT);
                }
            }
        }
    }
}

void SynchronousCouplingHandsOverTheNewWindAndDropsAVanishedLatentFlux ()
{
    auto lev = make_level(false);
    StandInEngine* eng = nullptr;
    auto cpl = make_coupler(*lev, coupling("synchronous", false), eng);
    step(*cpl, *lev);
    EXPECT_EQ(eng->last_xvel, &lev->xnew) << "synchronous coupling hands over the wind after the dycore";
    EXPECT_EQ(eng->last_yvel, &lev->ynew);
    EXPECT_NEAR(cpl->get_Q_lat_atm_prev()->max(0), QLAT, REL * QLAT);
    // the next step reports no latent flux: last step's must not linger
    eng->report_latent = false;
    step(*cpl, *lev);
    EXPECT_EQ(cpl->get_Q_lat_atm_prev()->norm0(0, 1, IntVect(5, 5, 0)), 0.0_rt);
    EXPECT_GT(cpl->get_Q_atm_prev()->min(0), 0.0_rt) << "the sensible flux still arrives";
    // A missing sensible pointer does not cancel an independently valid latent flux.
    eng->report_latent = true;
    step(*cpl, *lev);
    eng->report_sensible = false;
    step(*cpl, *lev);
    EXPECT_EQ(cpl->get_Q_atm_prev()->norm0(0, 1, IntVect(5, 5, 0)), 0.0_rt) << "no sensible flux this step";
    EXPECT_NEAR(cpl->get_Q_lat_atm_prev()->max(0), QLAT, REL * QLAT)
        << "latent flux must remain available without sensible heat";
    // When the remaining channel disappears too, clear its valid cells and halo.
    eng->report_latent = false;
    step(*cpl, *lev);
    EXPECT_EQ(cpl->get_Q_atm_prev()->norm0(0, 1, IntVect(5, 5, 0)), 0.0_rt);
    EXPECT_EQ(cpl->get_Q_lat_atm_prev()->norm0(0, 1, IntVect(5, 5, 0)), 0.0_rt);
}

void LatentOnlyFluxInjectsWaterWithoutHeating ()
{
    auto lev = make_level(false);
    FireCouplingParams p = coupling("lagged", false);
    p.heat_tendency_density = false; // Use the physical mass-flux conversion.
    StandInEngine* eng = nullptr;
    auto cpl = make_coupler(*lev, p, eng);

    eng->report_sensible = false;
    step(*cpl, *lev);
    const IntVect ghost(5, 5, 0);
    EXPECT_EQ(cpl->get_Q_atm_prev()->norm0(0, 1, ghost), 0.0_rt);
    EXPECT_NEAR(cpl->get_Q_lat_atm_prev()->min(0), QLAT, REL * QLAT);
    // A nonzero maximum cannot establish that *all* ghost cells were filled:
    // subtract the uniform expected flux over valid cells and the full halo.
    const MultiFab& latent_field = *cpl->get_Q_lat_atm_prev();
    MultiFab latent_residual(latent_field.boxArray(), latent_field.DistributionMap(), 1, ghost);
    latent_residual.setVal(QLAT);
    MultiFab::Subtract(latent_residual, latent_field, 0, 0, 1, ghost);
    EXPECT_LE(latent_residual.norm0(0, 1, ghost), REL * QLAT)
        << "latent valid and ghost cells must all contain the prescribed flux";

    // Native MultiFab reductions work on CPU/GPU and across MPI ranks.
    auto src = sources(*cpl, *lev, 0.0_rt);
    EXPECT_EQ(src->norm0(RhoTheta_comp), 0.0_rt);
    const Real retained = 1.0_rt - std::exp(-Real(NZ) * DZ / ALFG);
    const Real expected = Real(NX * NY) * FEEDBACK * (QLAT / L_V) * retained;
    EXPECT_NEAR(src->sum(RhoQ1_comp) * DZ, expected, REL * expected)
        << "prescribed latent-only flux must be deposited as water mass";

    // An allocated zero sensible field and a null sensible pointer mean the
    // same physical heat release: neither should affect vapor deposition.
    eng->report_sensible = true;
    eng->sensible->setVal(0.0_rt);
    step(*cpl, *lev);
    auto src_zero_sensible = sources(*cpl, *lev, 0.0_rt);
    EXPECT_EQ(src_zero_sensible->norm0(RhoTheta_comp), 0.0_rt);
    EXPECT_NEAR(src_zero_sensible->sum(RhoQ1_comp) * DZ, expected, REL * expected);

    eng->report_latent = false;
    step(*cpl, *lev);
    auto src_cleared = sources(*cpl, *lev, 0.0_rt);
    EXPECT_EQ(src_cleared->norm0(RhoTheta_comp), 0.0_rt);
    EXPECT_EQ(src_cleared->norm0(RhoQ1_comp), 0.0_rt)
        << "a disappeared latent flux must not be reused in a later timestep";
}

void TheInjectedHeatIsTheSurfaceFlux ()
{
    auto lev = make_level(false);
    StandInEngine* eng = nullptr;
    auto cpl = make_coupler(*lev, coupling("lagged", true), eng);
    step(*cpl, *lev);

    // The magnitudes on a zero background: the tendencies (1e-6 for the vapour
    // and the smoke) would round away on a background of 7 in single precision
    auto src = sources(*cpl, *lev, 0.0_rt);
    const Real frac = 1.0_rt - std::exp(-Real(NZ) * DZ / ALFG);
    for_each_column(*src, RhoTheta_comp, 0.0_rt, [&] (int i, int j, Real sum, const auto& s) {
        const Real q = flux_at_column(i);
        EXPECT_NEAR(Cp_fire * sum, RHO * FEEDBACK * q * frac, REL * q) << "sensible, column " << i;
        // the smoke: e Q / h in the lowest cell (not scaled by the feedback), spread over its thickness
        const Real smoke = EMISSION * q / HEAT_PER_KG;
        EXPECT_NEAR(s(i, j, 0, RhoSmoke_comp) * DZ, smoke, REL * smoke) << "smoke, column " << i;
        EXPECT_EQ(s(i, j, 1, RhoSmoke_comp), 0.0_rt);
    });
    for_each_column(*src, RhoQ1_comp, 0.0_rt, [&] (int i, int, Real sum, const auto&) {
        EXPECT_NEAR(L_V * sum, RHO * FEEDBACK * QLAT * frac, REL * QLAT) << "latent, column " << i;
    });

    // add mode keeps every other source: on a background B the theta source is
    // B more than on a zero background, and the sources the fire does not own are B
    auto srcB = sources(*cpl, *lev);
    const Real tol_b = 4.0_rt * std::numeric_limits<Real>::epsilon() * BACKGROUND;
    for (MFIter mfi(*srcB); mfi.isValid(); ++mfi) {
        const auto b = srcB->const_array(mfi);
        const auto z = src->const_array(mfi);
        const Box& bx = mfi.validbox();
        for (int k = 0; k < NZ; k += 7) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                EXPECT_NEAR(b(i, 0, k, RhoTheta_comp) - z(i, 0, k, RhoTheta_comp), BACKGROUND, tol_b) << "add mode, theta";
                EXPECT_NEAR(b(i, 0, k, RhoQ1_comp)    - z(i, 0, k, RhoQ1_comp),    BACKGROUND, tol_b) << "add mode, vapour";
                EXPECT_EQ(b(i, 0, k, Rho_comp), BACKGROUND) << "no other source is touched";
            }
        }
    }

    // overwrite mode, with the density factor off (the default): the theta and
    // vapour sources are replaced by f Q (1 - e^{-z_top/alfg}) without rho, the
    // others kept
    FireCouplingParams p = coupling("lagged", false);
    p.source_mode = "overwrite";
    p.heat_tendency_density = false;
    StandInEngine* eng2 = nullptr;
    auto cpl2 = make_coupler(*lev, p, eng2);
    step(*cpl2, *lev);
    auto src2 = sources(*cpl2, *lev);
    for_each_column(*src2, RhoTheta_comp, 0.0_rt, [&] (int i, int j, Real sum, const auto& s) {
        const Real q = flux_at_column(i);
        EXPECT_NEAR(Cp_fire * sum, FEEDBACK * q * frac, REL * q) << "overwritten theta source, column " << i;
        EXPECT_EQ(s(i, j, 0, Rho_comp), BACKGROUND);
        EXPECT_EQ(s(i, j, 0, RhoSmoke_comp), BACKGROUND) << "no smoke without smoke_enable";
    });
    for_each_column(*src2, RhoQ1_comp, 0.0_rt, [&] (int i, int, Real sum, const auto&) {
        EXPECT_NEAR(L_V * sum, FEEDBACK * QLAT * frac, REL * QLAT) << "overwritten vapour source, column " << i;
    });
}

void TheOpenFractionKeepsTheColumnEnergyAndThinsTheCellsUnderTheRoof ()
{
    // two cells under the roof: the roof 2 dz above the ground (the height an
    // engine reports), so cells 0 and 1 carry the open share alone (with the
    // first-centre datum used before, the roof sat half a cell higher and cut
    // cell 2)
    constexpr Real F_OPEN = 0.25;
    constexpr Real H_ROOF = 2.0 * DZ;
    auto lev = make_level(false);

    StandInEngine* eng0 = nullptr;
    auto cpl0 = make_coupler(*lev, coupling("lagged", false), eng0);
    step(*cpl0, *lev);
    auto src0 = sources(*cpl0, *lev, 0.0_rt);

    FireCouplingParams p = coupling("lagged", false);
    p.heat_open_fraction = true;
    StandInEngine* eng = nullptr;
    auto cpl = make_coupler(*lev, p, eng);
    eng->open_fraction = std::make_unique<MultiFab>(lev->ba2d, lev->dm, 1, 1);
    eng->open_fraction->setVal(F_OPEN);
    eng->roof_height = std::make_unique<MultiFab>(lev->ba2d, lev->dm, 1, 1);
    eng->roof_height->setVal(H_ROOF);
    step(*cpl, *lev);
    auto src = sources(*cpl, *lev, 0.0_rt);

    const Real frac  = 1.0_rt - std::exp(-Real(NZ) * DZ / ALFG);
    const Real norm  = frac - (1.0_rt - F_OPEN) * (1.0_rt - std::exp(-2.0_rt * DZ / ALFG));
    const Real under = F_OPEN * frac / norm;
    const Real above = frac / norm;
    for_each_column(*src, RhoTheta_comp, 0.0_rt, [&] (int i, int, Real sum, const auto&) {
        const Real q = flux_at_column(i);
        EXPECT_NEAR(Cp_fire * sum, RHO * FEEDBACK * q * frac, REL * q) << "column energy, column " << i;
    });
    for_each_column(*src, RhoQ1_comp, 0.0_rt, [&] (int i, int, Real sum, const auto&) {
        EXPECT_NEAR(L_V * sum, RHO * FEEDBACK * QLAT * frac, REL * QLAT) << "column vapour, column " << i;
    });
    for (MFIter mfi(*src); mfi.isValid(); ++mfi) {
        const auto a = src->const_array(mfi);
        const auto b = src0->const_array(mfi);
        const Box& bx = mfi.validbox();
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
            for (int k : {0, 1}) {
                EXPECT_NEAR(a(i, 0, k, RhoTheta_comp), under * b(i, 0, k, RhoTheta_comp),
                            REL * b(i, 0, k, RhoTheta_comp)) << "under the roof, cell " << k;
            }
            for (int k : {2, 9}) {
                EXPECT_NEAR(a(i, 0, k, RhoTheta_comp), above * b(i, 0, k, RhoTheta_comp),
                            REL * b(i, 0, k, RhoTheta_comp)) << "above the roof, cell " << k;
            }
        }
    }

    // the same weights from the level's nodal heights (k dz, as ERF passes
    // them), the faces taken from the nodes and not from the centres
    MultiFab znd(convert(lev->ba, IntVect(1, 1, 1)), lev->dm, 1, 1);
    for (MFIter mfi(znd); mfi.isValid(); ++mfi) {
        auto z = znd.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            z(i, j, k) = k * DZ;
        });
    }
    auto src_nd = sources(*cpl, *lev, 0.0_rt, &znd);
    for (MFIter mfi(*src); mfi.isValid(); ++mfi) {
        const auto a = src->const_array(mfi);
        const auto n = src_nd->const_array(mfi);
        const Box& bx = mfi.validbox();
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
            for (int k = 0; k < NZ; ++k) {
                EXPECT_NEAR(n(i, 0, k, RhoTheta_comp), a(i, 0, k, RhoTheta_comp),
                            REL * std::abs(a(i, 0, 0, RhoTheta_comp))) << "nodal heights, cell " << k;
            }
        }
    }

    // the same column on raised ground that varies from column to column
    // (37 + 3 i m): the roof is above the ground, so every cell's source is
    // the flat one (a roof or faces measured from z = 0 would cut other cells)
    MultiFab znd_g(convert(lev->ba, IntVect(1, 1, 1)), lev->dm, 1, 1);
    MultiFab zcc_g(lev->ba, lev->dm, 1, 1);
    for (MFIter mfi(znd_g); mfi.isValid(); ++mfi) {
        auto z = znd_g.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            z(i, j, k) = 37.0_rt + 3.0_rt * i + k * DZ;
        });
    }
    for (MFIter mfi(zcc_g); mfi.isValid(); ++mfi) {
        auto z = zcc_g.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            z(i, j, k) = 37.0_rt + 3.0_rt * (i + 0.5_rt) + (k + 0.5_rt) * DZ;
        });
    }
    auto src_g = std::make_unique<MultiFab>(lev->ba, lev->dm, NVAR_max, 0);
    src_g->setVal(0.0_rt);
    cpl->add_sources(*src_g, lev->S, zcc_g, &znd_g, lev->geom, true, 1);
    for (MFIter mfi(*src); mfi.isValid(); ++mfi) {
        const auto a = src->const_array(mfi);
        const auto g = src_g->const_array(mfi);
        const Box& bx = mfi.validbox();
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
            for (int k = 0; k < NZ; ++k) {
                EXPECT_NEAR(g(i, 0, k, RhoTheta_comp), a(i, 0, k, RhoTheta_comp),
                            REL * std::abs(a(i, 0, 0, RhoTheta_comp))) << "raised ground, cell " << k;
            }
        }
    }
}

/// Heights of the stretched test mesh: nodes z = h(i, j) + F(k) (1 + 0.05 i)
/// with h = 37 + 3 i + 2 j and F stretched by 1.1 per cell up to NZ dz
Real stretched_F (int k)
{
    const Real r = 1.1_rt;
    const Real d0 = Real(NZ) * DZ * (r - 1.0_rt) / (std::pow(r, Real(NZ)) - 1.0_rt);
    return d0 * (std::pow(r, Real(k)) - 1.0_rt) / (r - 1.0_rt);
}
/// Mean of the four lower (k) nodes of column (i, j): the face the coupling uses
Real stretched_face (int i, int j, int k)
{
    return 37.0_rt + 3.0_rt * (i + 0.5_rt) + 2.0_rt * (j + 0.5_rt)
         + stretched_F(k) * (1.0_rt + 0.05_rt * (i + 0.5_rt));
}

void TheOpenShareFollowsTheCellFacesOnAStretchedGrid ()
{
    // Nodes stretched by 1.1 per cell (a 0.43 m first cell, under a tenth of
    // the 10 m mean spacing) on raised ground whose corners differ, and the
    // roof 0.7 m above the ground, inside cell 1 of every column. Against the
    // same coupling without the open fraction, each cell's tendency is scaled
    // by its open share times one factor per column, so the ratio of a cell's
    // scaling to the top cell's is its share: f below the roof, f + (1 - f)
    // (z_hi - H) / dz in the cell that holds it (dz its own thickness), 1
    // above. The column still receives the flux's share.
    constexpr Real F_OPEN = 0.25, H_ROOF = 0.7;
    auto lev = make_level(false);
    MultiFab znd(convert(lev->ba, IntVect(1, 1, 1)), lev->dm, 1, 1), zcc(lev->ba, lev->dm, 1, 1);
    const Real r = 1.1_rt;
    const Real d0 = Real(NZ) * DZ * (r - 1.0_rt) / (std::pow(r, Real(NZ)) - 1.0_rt);
    for (MFIter mfi(znd); mfi.isValid(); ++mfi) {
        auto z = znd.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real F = d0 * (std::pow(r, Real(k)) - 1.0_rt) / (r - 1.0_rt);
            z(i, j, k) = 37.0_rt + 3.0_rt * i + 2.0_rt * j + F * (1.0_rt + 0.05_rt * i);
        });
    }
    for (MFIter mfi(zcc); mfi.isValid(); ++mfi) {
        auto z = zcc.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real Flo = d0 * (std::pow(r, Real(k))     - 1.0_rt) / (r - 1.0_rt);
            const Real Fhi = d0 * (std::pow(r, Real(k + 1)) - 1.0_rt) / (r - 1.0_rt);
            z(i, j, k) = 37.0_rt + 3.0_rt * (i + 0.5_rt) + 2.0_rt * (j + 0.5_rt)
                       + 0.5_rt * (Flo + Fhi) * (1.0_rt + 0.05_rt * (i + 0.5_rt));
        });
    }

    StandInEngine* eng0 = nullptr;
    auto cpl0 = make_coupler(*lev, coupling("lagged", false), eng0);
    step(*cpl0, *lev);
    FireCouplingParams p = coupling("lagged", false);
    p.heat_open_fraction = true;
    StandInEngine* eng = nullptr;
    auto cpl = make_coupler(*lev, p, eng);
    eng->open_fraction = std::make_unique<MultiFab>(lev->ba2d, lev->dm, 1, 1);
    eng->open_fraction->setVal(F_OPEN);
    eng->roof_height = std::make_unique<MultiFab>(lev->ba2d, lev->dm, 1, 1);
    eng->roof_height->setVal(H_ROOF);
    step(*cpl, *lev);
    MultiFab src0(lev->ba, lev->dm, NVAR_max, 0), src(lev->ba, lev->dm, NVAR_max, 0);
    src0.setVal(0.0_rt); src.setVal(0.0_rt);
    cpl0->add_sources(src0, lev->S, zcc, &znd, lev->geom, true, 1);
    cpl->add_sources(src, lev->S, zcc, &znd, lev->geom, true, 1);

    const Real tol_share = (sizeof(Real) == 8) ? 1.0e2 * REL : 1.0e-2;
    for (MFIter mfi(src); mfi.isValid(); ++mfi) {
        const auto a = src.const_array(mfi, RhoTheta_comp);
        const auto b = src0.const_array(mfi, RhoTheta_comp);
        const Box& bx = mfi.validbox();
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                const Real g = stretched_face(i, j, 0);
                ASSERT_LT(stretched_face(i, j, 1) - g, H_ROOF) << "cell 0 lies under the roof";
                ASSERT_GT(stretched_face(i, j, 2) - g, H_ROOF) << "the roof lies inside cell 1";
                const Real top = a(i, j, NZ - 1) / b(i, j, NZ - 1);
                for (int k = 0; k < 4; ++k) {
                    const Real zlo = stretched_face(i, j, k) - g, zhi = stretched_face(i, j, k + 1) - g;
                    Real share = 1.0_rt;
                    if (zhi <= H_ROOF)     { share = F_OPEN; }
                    else if (zlo < H_ROOF) { share = F_OPEN + (1.0_rt - F_OPEN) * (zhi - H_ROOF) / (zhi - zlo); }
                    EXPECT_NEAR(a(i, j, k) / b(i, j, k) / top, share, tol_share)
                        << "open share of cell " << k << ", column " << i << "," << j;
                }
                // the column total: the profile measured from the first centre,
                // closed one mean spacing above the top centre
                const Real zc0   = 0.5_rt * (stretched_face(i, j, 0) + stretched_face(i, j, 1));
                const Real zctop = 0.5_rt * (stretched_face(i, j, NZ - 1) + stretched_face(i, j, NZ));
                const Real frac  = 1.0_rt - std::exp(-((zctop - zc0) + DZ) / ALFG);
                Real sum = 0.0_rt;
                for (int k = 0; k < NZ; ++k) {
                    sum += Cp_fire * a(i, j, k) * (stretched_face(i, j, k + 1) - stretched_face(i, j, k));
                }
                const Real q = flux_at_column(i);
                const Real tol_sum = (sizeof(Real) == 8) ? 1.0e3 * REL * q : 1.0e-3 * q;
                EXPECT_NEAR(sum, RHO * FEEDBACK * q * frac, tol_sum) << "column energy, column " << i;
            }
        }
    }
}

void PassiveCouplingInjectsNoHeat ()
{
    // with the smoke tracer: the flux is still averaged, the smoke injected, the heat not
    {
        auto lev = make_level(false);
        StandInEngine* eng = nullptr;
        auto cpl = make_coupler(*lev, coupling("passive", true), eng);
        step(*cpl, *lev);
        EXPECT_EQ(eng->n_fill, 1);
        EXPECT_NEAR(cpl->get_Q_atm_prev()->max(0), flux_at_column(NX - 1), REL * flux_at_column(NX - 1));
        auto src = sources(*cpl, *lev);
        EXPECT_EQ(src->min(RhoTheta_comp), BACKGROUND);
        EXPECT_EQ(src->max(RhoTheta_comp), BACKGROUND);
        EXPECT_EQ(src->max(RhoQ1_comp), BACKGROUND);
        EXPECT_GT(src->max(RhoSmoke_comp), BACKGROUND);
    }
    // without it: nothing to ask the engine for, nothing touched
    {
        auto lev = make_level(false);
        StandInEngine* eng = nullptr;
        auto cpl = make_coupler(*lev, coupling("passive", false), eng);
        step(*cpl, *lev);
        EXPECT_EQ(eng->n_advance, 1);
        EXPECT_EQ(eng->n_fill, 0) << "no heat and no smoke reach the atmosphere: nothing to ask for";
        EXPECT_EQ(cpl->get_Q_atm_prev()->norm0(0, 1, IntVect(5, 5, 0)), 0.0_rt);
        auto src = sources(*cpl, *lev);
        for (int comp : {static_cast<int>(RhoTheta_comp), static_cast<int>(RhoQ1_comp), static_cast<int>(RhoSmoke_comp)}) {
            EXPECT_EQ(src->min(comp), BACKGROUND);
            EXPECT_EQ(src->max(comp), BACKGROUND);
        }
    }
}

void InjectLatentOffDropsTheEnginesLatentFlux ()
{
    auto lev = make_level(false);
    FireCouplingParams p = coupling("lagged", false);
    p.inject_latent = false;
    StandInEngine* eng = nullptr;
    auto cpl = make_coupler(*lev, p, eng);
    step(*cpl, *lev);
    EXPECT_EQ(cpl->get_Q_lat_atm_prev()->norm0(0, 1, IntVect(5, 5, 0)), 0.0_rt);
    auto src = sources(*cpl, *lev);
    EXPECT_EQ(src->min(RhoQ1_comp), BACKGROUND);
    EXPECT_EQ(src->max(RhoQ1_comp), BACKGROUND);
    EXPECT_GT(src->max(RhoTheta_comp), BACKGROUND) << "the sensible heat still arrives";
}

void TheRainIsHandedOverOnlyWhenAsked ()
{
    auto lev = make_level(true);
    StandInEngine* eng = nullptr;
    auto cpl = make_coupler(*lev, coupling("lagged", false), eng);
    step(*cpl, *lev);
    EXPECT_FALSE(eng->got_precip) << "not asked, not given";
    eng->wants_precip = true;
    step(*cpl, *lev);
    ASSERT_TRUE(eng->got_precip);
    EXPECT_EQ(eng->precip_ba, lev->ba2d) << "one k = 0 slab per box of the level";
    EXPECT_EQ(eng->precip_ncomp, 1);
}

/// A fresh checkpoint directory with the level directory VisMF writes into, under
/// the system temporary directory and named per process (the same test runs in
/// several CTest entries at once; erf_gtest_temp_path)
std::string make_chk_dir (const std::string& stem)
{
    const std::string name = erf_gtest_temp_path(stem).string();
    if (ParallelDescriptor::IOProcessor()) {
        UtilCreateDirectory(name + "/Level_0", 0755);
    }
    ParallelDescriptor::Barrier();
    return name;
}

void TheCheckpointRestoresTheLaggedFluxes ()
{
    auto lev = make_level(false);
    StandInEngine* eng = nullptr;
    auto cpl = make_coupler(*lev, coupling("lagged", false), eng);
    step(*cpl, *lev);
    const std::string chk = make_chk_dir("fire_coupler_chk");
    cpl->write_checkpoint(chk);
    ParallelDescriptor::Barrier();   // the headers are written by one rank
    EXPECT_EQ(eng->n_write_chk, 1) << "the coupler writes the engine's state with its own";

    // a fresh coupler and engine read it back: the engine first, then the fluxes,
    // bit for bit on the valid columns and refilled on the ghost columns
    StandInEngine* eng2 = nullptr;
    auto cpl2 = make_coupler(*lev, coupling("lagged", false), eng2);
    eng2->owner = cpl2.get();
    cpl2->read_checkpoint(chk, 0, 7);
    EXPECT_EQ(eng2->n_read_chk, 1);
    EXPECT_EQ(eng2->flux_at_read, 0.0_rt) << "the engine is restored before the lagged fluxes";
    EXPECT_EQ(eng2->read_step, 7) << "the engine gets the step of its level";
    for (int which : {0, 1}) {
        const MultiFab& a = which ? *cpl->get_Q_lat_atm_prev() : *cpl->get_Q_atm_prev();
        const MultiFab& b = which ? *cpl2->get_Q_lat_atm_prev() : *cpl2->get_Q_atm_prev();
        MultiFab d(a.boxArray(), a.DistributionMap(), 1, a.nGrowVect());
        MultiFab::Copy(d, a, 0, 0, 1, a.nGrowVect());
        MultiFab::Subtract(d, b, 0, 0, 1, a.nGrowVect());
        EXPECT_EQ(d.norm0(0, 1, a.nGrowVect()), 0.0_rt) << (which ? "latent" : "sensible") << " flux, ghosts included";
        EXPECT_GT(b.min(0), 0.0_rt);
    }

    // a checkpoint written before the fluxes carried ghost columns holds the
    // valid columns only: the restore refills the ghosts (the MRF thermal
    // excess reads them) to the values the run itself had
    const std::string old = make_chk_dir("fire_coupler_chk_noghost");
    for (int which : {0, 1}) {
        const MultiFab& a = which ? *cpl->get_Q_lat_atm_prev() : *cpl->get_Q_atm_prev();
        MultiFab valid(a.boxArray(), a.DistributionMap(), 1, 0);
        MultiFab::Copy(valid, a, 0, 0, 1, 0);
        VisMF::Write(valid, MultiFabFileFullPrefix(0, old, "Level_", which ? "FireQLatAtmPrev" : "FireQAtmPrev"));
    }
    ParallelDescriptor::Barrier();
    StandInEngine* eng4 = nullptr;
    auto cpl4 = make_coupler(*lev, coupling("lagged", false), eng4);
    cpl4->read_checkpoint(old, 0, 0);
    for (int which : {0, 1}) {
        const MultiFab& a = which ? *cpl->get_Q_lat_atm_prev()  : *cpl->get_Q_atm_prev();
        const MultiFab& b = which ? *cpl4->get_Q_lat_atm_prev() : *cpl4->get_Q_atm_prev();
        MultiFab d(a.boxArray(), a.DistributionMap(), 1, a.nGrowVect());
        MultiFab::Copy(d, a, 0, 0, 1, a.nGrowVect());
        MultiFab::Subtract(d, b, 0, 0, 1, a.nGrowVect());
        EXPECT_EQ(d.norm0(0, 1, a.nGrowVect()), 0.0_rt)
            << (which ? "latent" : "sensible") << " ghost columns refilled from a checkpoint without them";
    }

    // a checkpoint without them (no fire state): the fluxes stay zero
    const std::string empty = make_chk_dir("fire_coupler_chk_empty");
    StandInEngine* eng3 = nullptr;
    auto cpl3 = make_coupler(*lev, coupling("lagged", false), eng3);
    cpl3->read_checkpoint(empty, 0, 0);
    EXPECT_EQ(eng3->n_read_chk, 1) << "the engine is asked whatever the coupling finds";
    EXPECT_EQ(cpl3->get_Q_atm_prev()->norm0(0, 1, IntVect(5, 5, 0)), 0.0_rt);

    if (ParallelDescriptor::IOProcessor()) {
        FileSystem::RemoveAll(chk);
        FileSystem::RemoveAll(empty);
        FileSystem::RemoveAll(old);
    }
    ParallelDescriptor::Barrier();
}

void ACheckpointOfAnotherFireGridIsRefused ()
{
    // the coupler on the region of columns 2-5; its own region and level pass
    auto lev = make_level(true);
    StandInEngine* eng = nullptr;
    auto cpl = make_coupler(*lev, coupling("lagged", false), eng);
    const Box own(IntVect(2, 0, 0), IntVect(5, NY - 1, 0));
    EXPECT_TRUE(cpl->checkpoint_mismatch(0, own).empty());
    // fluxes of the whole domain, or of the same columns on another level, do not
    const Box whole(IntVect(0, 0, 0), IntVect(NX - 1, NY - 1, 0));
    EXPECT_NE(cpl->checkpoint_mismatch(0, whole).find("columns"), std::string::npos);
    EXPECT_NE(cpl->checkpoint_mismatch(1, own).find("on level 1"), std::string::npos);
}

/// Keys added to erf.fire for one test and removed on exit
class ScopedKeys
{
public:
    ScopedKeys () = default;
    ~ScopedKeys () { for (const auto& n : m_names) { m_pp.remove(n.c_str()); } }
    ScopedKeys (const ScopedKeys&) = delete;
    ScopedKeys& operator= (const ScopedKeys&) = delete;
    void add (const char* key, const std::string& v) { m_pp.add(key, v); m_names.emplace_back(key); }
    void add (const char* key, double v)             { m_pp.add(key, v); m_names.emplace_back(key); }
    void add (const char* key, int v)                { m_pp.add(key, v); m_names.emplace_back(key); }
    void add (const char* key, bool v)               { m_pp.add(key, v); m_names.emplace_back(key); }
private:
    ParmParse                m_pp{"erf.fire"};
    std::vector<std::string> m_names;
};

void TheFireInputsReachTheCoupling ()
{
    // the two tendency factors (and the open-fraction switch) are set opposite
    // ways in two decks, so none can stand in for the other or for its default
    for (bool density : {true, false}) {
        ScopedKeys deck;
        deck.add("grid_ratio", 4);
        deck.add("anchor_level", 0);
        deck.add("coupling_type", std::string("synchronous"));
        deck.add("heat_flux_alfg", 30.0);
        deck.add("fire_atm_feedback", 0.25);
        deck.add("inject_latent", false);
        deck.add("source_mode", std::string("add"));
        deck.add("heat_tendency_density", density);
        deck.add("heat_tendency_exner", !density);
        deck.add("heat_open_fraction", density);
        deck.add("smoke_enable", true);
        deck.add("smoke_emission_factor", 0.05);
        deck.add("fire_debug", true);
        const FireCouplingParams c = FireCouplingParams::read();
        // no engine is built into this ERF, so the engine name stays empty
        EXPECT_FALSE(c.enable);
        EXPECT_EQ(c.engine, "");
        EXPECT_EQ(c.grid_ratio, 4);
        EXPECT_EQ(c.anchor_level, 0);
        EXPECT_EQ(c.coupling_type, "synchronous");
        EXPECT_TRUE(c.is_synchronous());
        EXPECT_TRUE(c.injects_flux());
        EXPECT_NEAR(c.heat_flux_alfg, Real(30.0), REL * Real(30.0));
        EXPECT_NEAR(c.fire_atm_feedback, Real(0.25), REL * Real(0.25));
        EXPECT_FALSE(c.inject_latent);
        EXPECT_EQ(c.source_mode, "add");
        EXPECT_TRUE(c.source_mode_given);
        EXPECT_EQ(c.heat_open_fraction, density);
        EXPECT_EQ(c.heat_tendency_density, density);
        EXPECT_EQ(c.heat_tendency_exner, !density);
        EXPECT_TRUE(c.smoke_enable);
        EXPECT_NEAR(c.smoke_emission_factor, Real(0.05), REL * Real(0.05));
        EXPECT_TRUE(c.fire_debug);
    }
}

} // namespace

TEST(FireCoupler, FluxesAverageOntoTheColumnsTheyCameFrom)        { FluxesAverageOntoTheColumnsTheyCameFrom(); }
TEST(FireCoupler, SynchronousCouplingHandsOverTheNewWindAndDropsVanishedFluxes)
{
    SynchronousCouplingHandsOverTheNewWindAndDropsAVanishedLatentFlux();
}
TEST(FireCoupler, LatentOnlyFluxInjectsWaterWithoutHeating)        { LatentOnlyFluxInjectsWaterWithoutHeating(); }
TEST(FireCoupler, TheInjectedHeatIsTheSurfaceFlux)                 { TheInjectedHeatIsTheSurfaceFlux(); }
TEST(FireCoupler, PassiveCouplingInjectsNoHeat)                    { PassiveCouplingInjectsNoHeat(); }
TEST(FireCoupler, TheOpenFractionKeepsTheColumnEnergyAndThinsTheCellsUnderTheRoof)
{
    TheOpenFractionKeepsTheColumnEnergyAndThinsTheCellsUnderTheRoof();
}
TEST(FireCoupler, TheOpenShareFollowsTheCellFacesOnAStretchedGrid)
{
    TheOpenShareFollowsTheCellFacesOnAStretchedGrid();
}
TEST(FireCoupler, InjectLatentOffDropsTheEnginesLatentFlux)        { InjectLatentOffDropsTheEnginesLatentFlux(); }
TEST(FireCoupler, TheRainIsHandedOverOnlyWhenAsked)                { TheRainIsHandedOverOnlyWhenAsked(); }
TEST(FireCoupler, TheFireInputsReachTheCoupling)                   { TheFireInputsReachTheCoupling(); }
TEST(FireCoupler, TheCheckpointRestoresTheLaggedFluxes)            { TheCheckpointRestoresTheLaggedFluxes(); }
TEST(FireCoupler, ACheckpointOfAnotherFireGridIsRefused)           { ACheckpointOfAnotherFireGridIsRefused(); }

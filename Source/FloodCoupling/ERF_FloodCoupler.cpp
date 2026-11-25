#include <ERF_FloodCoupler.H>
#include <AMReX_VisMF.H>
#include <ERF_FloodEngineTemplate.H>
#include <ERF_FloodPrerequisites.H>
#include <ERF.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <fstream>
#include <set>
#include <sstream>
#include <type_traits>
#include <utility>

using namespace amrex;

// The template for new engines declares every method FloodEngine requires: a
// pure virtual method added to the interface and missing there makes the
// template abstract. Compile time only; no FloodEngineTemplate is ever created.
static_assert(!std::is_abstract_v<FloodEngineTemplate>,
              "ERF_FloodEngineTemplate.H must declare every pure virtual method of FloodEngine");

namespace {

// Where the flood grid sits, printed once at start-up, and the note when a
// finer level covers part of it: the flood then reads that level's coarser
// precipitation and wind.
void report_flood_region (const ERF& erf, int lev, const FloodGrid& fg)
{
    const Geometry& g = fg.geom;
    amrex::Print() << "[FLOOD] Flood grid on level " << lev << " of 0-" << erf.finestLevel()
                   << ": x " << g.ProbLo(0) << " to " << g.ProbHi(0)
                   << " m, y " << g.ProbLo(1) << " to " << g.ProbHi(1) << " m, "
                   << g.Domain().length(0) << " x " << g.Domain().length(1) << " flood cells of "
                   << g.CellSize(0) << " x " << g.CellSize(1) << " m\n";
    if (lev < erf.finestLevel()) {
        amrex::Print() << "[FLOOD] NOTE: the flood grid is on level " << lev << " but the finest level is "
                       << erf.finestLevel() << " (erf.flood.anchor_level = " << lev << "); the flood reads the "
                       << "precipitation and wind of level " << lev << ". Leave erf.flood.anchor_level unset "
                       << "to put the flood grid on the finest level.\n";
    }
}

} // namespace

FloodCoupler::FloodCoupler (const FloodCouplingParams& params, std::unique_ptr<FloodEngine> engine)
    : m_params(params), m_engine(std::move(engine))
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_engine != nullptr, "[FLOOD] FloodCoupler needs a flood engine");
    m_plot_var_names = m_engine->plot_var_names();
    std::set<std::string> seen;
    for (const auto& n : m_plot_var_names) {
        if (n.size() <= 6 || n.compare(0, 6, "flood_") != 0) {
            amrex::Abort(std::string("[FLOOD] flood engine \"") + m_engine->name() + "\" lists the plot field \""
                         + n + "\"; every flood plot field name must start with \"flood_\" and continue after it");
        }
        if (!seen.insert(n).second) {
            amrex::Abort(std::string("[FLOOD] flood engine \"") + m_engine->name() + "\" lists the plot field \""
                         + n + "\" twice");
        }
    }
    if (!m_plot_var_names.empty()) {
        amrex::Print() << "[FLOOD] " << m_engine->name()
                       << ": the following additional variables are available to plot:\n";
        for (const auto& n : m_plot_var_names) { amrex::Print() << "    " << n << "\n"; }
    }
}

void
FloodCoupler::initialize (const ERF& erf, int lev, const MultiFab* z_phys_nd_atm)
{
    // The ERF configuration a flood grid needs, then the engine's own inputs
    verify_flood_prerequisites(erf, lev);
    m_engine->verify_prerequisites(erf, lev);

    FloodGrid fg = create_flood_grid(erf.boxArray(lev), erf.DistributionMap(lev),
                                     erf.Geom(lev), m_params.grid_ratio);
    fg.lev = lev;
    report_flood_region(erf, lev, fg);
    set_flood_grid(lev, fg, erf.boxArray(lev), erf.DistributionMap(lev), erf.Geom(lev));

    m_engine->initialize(erf, lev, m_fg, z_phys_nd_atm);
}

void
FloodCoupler::set_flood_grid (int lev, const FloodGrid& fg,
                              const BoxArray& ba_atm,
                              const DistributionMapping& dm_atm,
                              const Geometry& geom_atm)
{
    m_lev      = lev;
    m_fg       = fg;
    m_geom_atm = geom_atm;
    const int klo = geom_atm.Domain().smallEnd(2);
    BoxList bl;
    for (int n = 0; n < ba_atm.size(); ++n) {
        Box b = ba_atm[n];
        b.setSmall(2, klo);
        b.setBig(2, klo);
        bl.push_back(b);
    }
    m_ba_cols  = BoxArray(std::move(bl));
    m_dm_cols  = dm_atm;
    m_grid_set = true;
    m_precip_prev.reset();
}

void
FloodCoupler::set_precip_reference (const SurfacePrecipAccumulationSources& precip_sources)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_grid_set, "[FLOOD] FloodCoupler: set_flood_grid() before set_precip_reference()");
    m_precip_prev = std::make_unique<MultiFab>(m_ba_cols, m_dm_cols, 2, 0);
    surface_precip_liquid_frozen_accum_k0(*m_precip_prev, precip_sources, m_geom_atm.Domain().smallEnd(2));
}

void
FloodCoupler::advance (double time_d, double dt_d,
                       const MultiFab& xvel, const MultiFab& yvel,
                       const MultiFab& z_phys_cc,
                       const SurfacePrecipAccumulationSources& precip_sources)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_grid_set, "[FLOOD] FloodCoupler: set_flood_grid() before advance()");
    // ERF keeps its time and step in double; the engine takes the time in
    // double too (a float clock drifts on long runs) and the step in Real
    const double time = time_d;
    const Real   dt   = static_cast<Real>(dt_d);

    // This step's precipitation: the accumulation now less the one at the end
    // of the previous step (or at start-up), which then rolls forward
    std::unique_ptr<MultiFab> precip;
    if (m_engine->wants_surface_precip()) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_precip_prev,
            "[FLOOD] FloodCoupler: set_precip_reference() before the first advance() of an engine that reads the precipitation");
        precip = std::make_unique<MultiFab>(m_ba_cols, m_dm_cols, 2, 0);
        surface_precip_liquid_frozen_accum_k0(*precip, precip_sources, m_geom_atm.Domain().smallEnd(2));
        for (MFIter mfi(*precip); mfi.isValid(); ++mfi) {
            const auto now  = precip->array(mfi);
            const auto prev = m_precip_prev->array(mfi);
            ParallelFor(mfi.validbox(), 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                const Real a = now(i, j, k, n);
                now(i, j, k, n)  = a - prev(i, j, k, n);
                prev(i, j, k, n) = a;
            });
        }
    }

    const FloodForcing forcing{time, dt, xvel, yvel, z_phys_cc, precip.get()};
    m_engine->advance(forcing);
    // time + dt[lev] in double, as ERF::timeStep sets t_new[lev]
    m_fields_time = time_d + dt_d;
}

void
FloodCoupler::fill_plot_var (const std::string& name, MultiFab& dst, int dcomp, int lev) const
{
    dst.setVal(0.0_rt, dcomp, 1, 0);
    if (lev != m_lev || !m_grid_set) { return; }
    const MultiFab* field = m_engine->plot_var(name);
    if (!field) { return; }

    MultiFab col(m_ba_cols, m_dm_cols, 1, 0);
    coarsen_flood_to_atm(col, *field, m_fg);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(flood_same_columns(dst.boxArray(), col.boxArray()) &&
                                     dst.DistributionMap() == col.DistributionMap(),
                                     "[FLOOD] fill_plot_var: the plotfile MultiFab is not on the flood grid's level");
    const int k0 = m_geom_atm.Domain().smallEnd(2);
    for (MFIter mfi(dst); mfi.isValid(); ++mfi) {
        Box bx = mfi.validbox();
        if (bx.smallEnd(2) > k0) { continue; }
        bx.setSmall(2, k0);
        bx.setBig(2, k0);
        const auto d = dst.array(mfi);
        const auto c = col.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            d(i, j, k, dcomp) = c(i, j, k);
        });
    }
}

void
FloodCoupler::write_checkpoint (const std::string& checkpointname) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_grid_set, "[FLOOD] FloodCoupler: set_flood_grid() before the checkpoint");
    // The engine's own state first
    m_engine->write_checkpoint(checkpointname);

    // The flood grid the state belongs to, so a restart onto another one stops
    if (ParallelDescriptor::IOProcessor()) {
        const std::string dir = checkpointname + "/Level_" + std::to_string(m_lev);
        if (!amrex::UtilCreateDirectory(dir, 0755)) { amrex::CreateDirectoryFailed(dir); }
        const Box region = m_ba_cols.minimalBox();
        std::ofstream os(dir + "/FloodGrid");
        os << m_lev << " " << m_fg.C << " "
           << region.smallEnd(0) << " " << region.smallEnd(1) << " "
           << region.bigEnd(0)   << " " << region.bigEnd(1)   << "\n";
        if (!os) { amrex::FileOpenFailed(dir + "/FloodGrid"); }
    }
    ParallelDescriptor::Barrier();

    // The precipitation reference: the accumulation the flood last read. Under
    // a finer level ERF averages the accumulators down after the flood's step,
    // and the next step hands that correction to the engine; saving the
    // reference keeps the correction for the step after a restart.
    if (m_precip_prev) {
        VisMF::Write(*m_precip_prev, checkpointname + "/Level_" + std::to_string(m_lev) + "/FloodPrecipRef");
    }
}

std::string
FloodCoupler::checkpoint_mismatch (int chk_lev, int chk_ratio, const Box& chk_region) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_grid_set, "[FLOOD] FloodCoupler: set_flood_grid() before the checkpoint");
    const Box run_region = m_ba_cols.minimalBox();
    Box chk2 = chk_region;
    chk2.setSmall(2, run_region.smallEnd(2));
    chk2.setBig(2, run_region.bigEnd(2));
    if (chk_lev == m_lev && chk_ratio == m_fg.C && chk2 == run_region) { return {}; }
    std::ostringstream msg;
    msg << "holds flood state on level " << chk_lev << " with erf.flood.grid_ratio = " << chk_ratio
        << " over the columns " << chk2 << ", but this run's flood grid is on level " << m_lev
        << " with erf.flood.grid_ratio = " << m_fg.C << " over the columns " << run_region
        << ". Restart with the erf.flood.anchor_level, erf.flood.grid_ratio and the refinement the"
           " checkpoint was written with.";
    return msg.str();
}

void
FloodCoupler::read_checkpoint (const std::string& restart_chkfile, int step)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_grid_set, "[FLOOD] FloodCoupler: set_flood_grid() before the checkpoint");
    // The flood grid record: found and read on the I/O rank and broadcast, so
    // every rank takes the same branch. Every level the checkpoint holds is
    // searched, not only those of this run's amr.max_level, so a restart onto
    // fewer levels stops on the mismatch instead of starting the engine empty.
    Vector<int> rec(8, -1);   // found, lev, C, ilo, jlo, ihi, jhi, has precipitation reference
    if (ParallelDescriptor::IOProcessor()) {
        for (int l = 0; amrex::FileExists(restart_chkfile + "/Level_" + std::to_string(l)); ++l) {
            const std::string f = restart_chkfile + "/Level_" + std::to_string(l) + "/FloodGrid";
            if (!amrex::FileExists(f)) { continue; }
            std::ifstream is(f);
            rec[0] = 1;
            is >> rec[1] >> rec[2] >> rec[3] >> rec[4] >> rec[5] >> rec[6];
            if (!is) { rec[0] = 0; }   // unreadable record
            rec[7] = amrex::FileExists(restart_chkfile + "/Level_" + std::to_string(l)
                                       + "/FloodPrecipRef_H") ? 1 : 0;
            break;
        }
    }
    ParallelDescriptor::Bcast(rec.data(), rec.size(), ParallelDescriptor::IOProcessorNumber());
    if (rec[0] == 0) {
        amrex::Abort("[FLOOD] The checkpoint " + restart_chkfile + " has a FloodGrid record that cannot be read");
    }
    if (rec[0] < 0) {
        amrex::Print() << "[FLOOD] The checkpoint holds no flood state; the flood engine "
                       << m_engine->name() << " starts from its initial state.\n";
        return;
    }
    const Box chk_region(IntVect(rec[3], rec[4], 0), IntVect(rec[5], rec[6], 0));
    if (const std::string why = checkpoint_mismatch(rec[1], rec[2], chk_region); !why.empty()) {
        amrex::Abort("[FLOOD] The checkpoint " + restart_chkfile + " " + why);
    }
    m_engine->read_checkpoint(restart_chkfile, step);
    // The reference the flood last read, which can differ from the restored
    // accumulators under a finer level (write_checkpoint)
    if (rec[7] == 1 && m_precip_prev) {
        MultiFab ref;
        VisMF::Read(ref, restart_chkfile + "/Level_" + std::to_string(m_lev) + "/FloodPrecipRef");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ref.nComp() == m_precip_prev->nComp(),
            "[FLOOD] The checkpoint's FloodPrecipRef has the wrong number of components");
        m_precip_prev->ParallelCopy(ref, 0, 0, ref.nComp());
    }
}

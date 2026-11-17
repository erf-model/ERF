#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <AMReX_GpuLaunch.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Math.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include "ERF_WRFSurfaceFields.H"

using namespace amrex;

namespace erf_detail_WRFSurfaceFields {

// Prefix of the fields' MultiFab files in a checkpoint level directory
const char* const chk_mf_prefix = "WRFSfc_";

// First line of the checkpoint's text file
const char* const chk_magic   = "WRFSurfaceFields";
constexpr int     chk_version = 1;

// Piecewise constant: every column of fine takes the column of crse it lies in. The
// caller has checked that the valid cells of crse cover the coarsened grids of fine.
void
fill_from_coarser (MultiFab& fine, const MultiFab& crse, const IntVect& ratio)
{
    const int ncomp = fine.nComp();
    const int rx = ratio[0];
    const int ry = ratio[1];
    BoxArray cba = fine.boxArray();
    cba.coarsen(IntVect(rx, ry, 1));   // the fields have no vertical extent
    MultiFab ctmp(cba, fine.DistributionMap(), ncomp, 0);
    ctmp.ParallelCopy(crse, 0, 0, ncomp, IntVect(0), IntVect(0));
    for (MFIter mfi(fine, false); mfi.isValid(); ++mfi) {
        const Array4<Real>       f = fine.array(mfi);
        const Array4<const Real> c = ctmp.const_array(mfi);
        ParallelFor(mfi.validbox(), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            f(i, j, k, n) = c(amrex::coarsen(i, rx), amrex::coarsen(j, ry), k, n);
        });
    }
}

} // namespace
using namespace erf_detail_WRFSurfaceFields;

const std::array<WRFSurfaceFields::Descriptor, WRFSurfaceFields::num_fields>&
WRFSurfaceFields::catalog ()
{
    // Not named zero and one: those are constants of ERF_NumericalConstants.H
    constexpr Real lo0 = Real(0.0);
    constexpr Real hi1 = Real(1.0);
    static const std::array<Descriptor, num_fields> table{{
        {Field::LAI,    "LAI",    Kind::Surface,     "m2/m2",  lo0,  Real(100.0)},
        {Field::VEGFRA, "VEGFRA", Kind::Surface,     "%",      lo0,  Real(100.0)},
        {Field::TMN,    "TMN",    Kind::Surface,     "K",      lo0,  Real(1000.0)},
        {Field::ZNT,    "ZNT",    Kind::Surface,     "m",      lo0,  Real(100.0)},
        {Field::ALBEDO, "ALBEDO", Kind::Surface,     "1",      lo0,  hi1},
        {Field::ALBBCK, "ALBBCK", Kind::Surface,     "1",      lo0,  hi1},
        {Field::SEAICE, "SEAICE", Kind::Surface,     "1",      lo0,  hi1},
        {Field::SNOW,   "SNOW",   Kind::Surface,     "kg/m2",  lo0,  Real(1.0e6)},
        {Field::SNOWH,  "SNOWH",  Kind::Surface,     "m",      lo0,  Real(1.0e4)},
        {Field::TSLB,   "TSLB",   Kind::SoilField,   "K",      lo0,  Real(1000.0)},
        {Field::SMOIS,  "SMOIS",  Kind::SoilField,   "m3/m3",  lo0,  hi1},
        {Field::ZS,     "ZS",     Kind::LayerColumn, "m",      lo0,  Real(1.0e4)},
        {Field::DZS,    "DZS",    Kind::LayerColumn, "m",      lo0,  Real(1.0e4)},
    }};
    return table;
}

const WRFSurfaceFields::Descriptor*
WRFSurfaceFields::find (const std::string& nc_name)
{
    for (const auto& d : catalog()) {
        if (nc_name == d.nc_name) { return &d; }
    }
    return nullptr;
}

void
WRFSurfaceFields::resize (int nlevs)
{
    m_fields.resize(nlevs);
    m_pending.resize(nlevs);
    m_soil_depth.resize(nlevs);
    m_soil_thickness.resize(nlevs);
}

bool
WRFSurfaceFields::holds_any (int lev) const
{
    if (lev < 0 || lev >= num_levels()) { return false; }
    for (const auto& f : m_fields[lev]) {
        if (f) { return true; }
    }
    return !m_soil_depth[lev].empty() || !m_soil_thickness[lev].empty();
}

void
WRFSurfaceFields::clear (int lev)
{
    if (lev < 0 || lev >= num_levels()) { return; }
    if (holds_any(lev)) {
        Print() << "   wrfinput land and soil fields at level " << lev << " are dropped\n";
    }
    for (auto& f : m_fields[lev]) { f.reset(); }
    m_pending[lev] = {};
    m_soil_depth[lev].clear();
    m_soil_thickness[lev].clear();
}

void
WRFSurfaceFields::drop (int lev, Field id)
{
    if (id == Field::ZS) {
        m_soil_depth[lev].clear();
    } else if (id == Field::DZS) {
        m_soil_thickness[lev].clear();
    } else {
        m_fields[lev][static_cast<int>(id)].reset();
    }
}

void
WRFSurfaceFields::set_soil_layers (int lev, Field id, const FArrayBox& src)
{
    AMREX_ALWAYS_ASSERT(lev >= 0 && lev < num_levels());
    AMREX_ALWAYS_ASSERT(id == Field::ZS || id == Field::DZS);
    const Descriptor& d = catalog()[static_cast<int>(id)];
    Pending& p = m_pending[lev][static_cast<int>(id)];
    auto& dst = (id == Field::DZS) ? m_soil_thickness[lev] : m_soil_depth[lev];
    const Box& sbx = src.box();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(sbx.ok(),
        "WRFSurfaceFields::set_soil_layers: this rank holds no part of the soil layer column");
    if (p.conflict) { return; }

    const auto arr = src.const_array();
    const int ilo  = sbx.smallEnd(0);
    const int jlo  = sbx.smallEnd(1);
    const int nlay = sbx.length(2);
    Vector<Real> vals(nlay);
    std::string reason;
    bool all_zero = true;
    for (int n = 0; n < nlay; ++n) {
        vals[n] = arr(ilo, jlo, sbx.smallEnd(2) + n);
        // isfinite first: a comparison with a NaN raises the invalid-operation exception
        if (!std::isfinite(vals[n])) {
            reason = "a value is not finite";
            break;
        }
        all_zero = all_zero && (vals[n] == Real(0.0));
    }
    if (reason.empty() && all_zero) {
        // A file written without a land model (ideal.exe, or a hand-built file) declares
        // the layers and leaves them at zero: that is a file without them.
        Print() << "   wrfinput " << d.nc_name << " at level " << lev
                << " is zero in every soil layer: taken as not given\n";
        return;
    }
    for (int n = 0; reason.empty() && n < nlay; ++n) {
        if (vals[n] > d.hi) {
            reason = "a value is above the catalog bound";
        } else if (id == Field::DZS && vals[n] <= Real(0.0)) {
            reason = "every soil layer thickness must be positive";
        } else if (id == Field::ZS && vals[n] < Real(0.0)) {
            reason = "a soil layer depth is negative";
        } else if (id == Field::ZS && n > 0 && vals[n] <= vals[n-1]) {
            // The top layer is first: the depths grow with the layer index. The first
            // one may be zero (a land model with its first level at the surface).
            reason = "the soil layer depths must increase, top layer first";
        }
    }
    if (reason.empty() && p.read && !dst.empty() && dst != vals) {
        reason = "two files of the level give different values";
    }
    if (!reason.empty()) {
        Print() << "Warning: wrfinput " << d.nc_name << " at level " << lev
                << " is not held: " << reason << "\n";
        dst.clear();
        p.conflict = true;
        return;
    }
    const bool first = dst.empty();
    dst = vals;
    p.read = true;
    if (first) {
        Print() << "   wrfinput soil layer "
                << ((id == Field::DZS) ? "thicknesses DZS" : "depths ZS")
                << " at level " << lev << ":";
        for (const Real v : dst) { Print() << " " << v; }
        Print() << " [m]\n";
    }
}

void
WRFSurfaceFields::fill (int lev, Field id,
                        const BoxArray& ba2d,
                        const DistributionMapping& dm,
                        const IntVect& ng,
                        int nlayers,
                        const FArrayBox& src,
                        const Box& file_box)
{
    AMREX_ALWAYS_ASSERT(lev >= 0 && lev < num_levels());
    const int idx = static_cast<int>(id);
    const Descriptor& d = catalog()[idx];
    AMREX_ALWAYS_ASSERT(d.kind != Kind::LayerColumn);
    Pending& p = m_pending[lev][idx];
    if (p.conflict) { return; }

    const int ncomp = (d.kind == Kind::SoilField) ? nlayers : 1;
    AMREX_ALWAYS_ASSERT(ncomp >= 1);
    auto& mf = m_fields[lev][idx];
    if (!p.read) {
        // The first file of the level that holds the variable
        mf = std::make_unique<MultiFab>(ba2d, dm, ncomp, ng);
        // Defined everywhere, also where no file covers the level (reported later)
        mf->setVal(Real(0.0));
        p.read = true;
    } else if (mf->nComp() != ncomp) {
        // The layer count is the file's metadata: every rank sees the same disagreement
        p.conflict = true;
        return;
    }

    const Box& sbx = src.box();
    const int klo  = sbx.ok() ? sbx.smallEnd(2) : 0;
    const Real lo  = d.lo;
    const Real hi  = d.hi;
    constexpr Real big = std::numeric_limits<Real>::max();

    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpMin, ReduceOpMax> reduce_op;
    ReduceData<Long, Long, Long, Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // The valid columns of this box that the file covers
        Box bx = mfi.tilebox();
        for (int dir = 0; dir < 2; ++dir) {
            bx.setSmall(dir, amrex::max(bx.smallEnd(dir), file_box.smallEnd(dir)));
            bx.setBig  (dir, amrex::min(bx.bigEnd(dir),   file_box.bigEnd(dir)));
        }
        if (!bx.ok()) { continue; }

        // The source is tested on the host and not read where it does not reach. (Its
        // index type in z differs from the field's for a soil variable, so the test is
        // made index by index.)
        const bool covered = sbx.ok() &&
            sbx.smallEnd(0) <= bx.smallEnd(0) && sbx.bigEnd(0) >= bx.bigEnd(0) &&
            sbx.smallEnd(1) <= bx.smallEnd(1) && sbx.bigEnd(1) >= bx.bigEnd(1) &&
            sbx.length(2) >= ncomp;
        if (!covered) {
            ++p.n_uncovered;
            continue;
        }

        const Array4<const Real> s = src.const_array();
        const Array4<Real>       a = mf->array(mfi);
        reduce_op.eval(bx, ncomp, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept -> ReduceTuple
        {
            const Real v = s(i, j, klo + n);
            a(i, j, k, n) = v;
            // isfinite first: a comparison with a NaN raises the invalid-operation
            // exception, which a run may trap
            if (!amrex::Math::isfinite(v)) {
                return {Long(1), Long(0), Long(1), big, -big};
            }
            const Long outside = (v < lo || v > hi) ? Long(1) : Long(0);
            return {Long(1), outside, Long(0), v, v};
        });
    }

    const ReduceTuple r = reduce_data.value(reduce_op);
    p.n_filled     += amrex::get<0>(r);
    p.n_outside    += amrex::get<1>(r);
    p.n_not_finite += amrex::get<2>(r);
    p.vmin = amrex::min(p.vmin, amrex::get<3>(r));
    p.vmax = amrex::max(p.vmax, amrex::get<4>(r));
}

void
WRFSurfaceFields::finish_level (int lev, const Periodicity& period)
{
    AMREX_ALWAYS_ASSERT(lev >= 0 && lev < num_levels());
    auto& pend = m_pending[lev];

    // Two reductions for the level: the counts of every field, then the ranges
    constexpr int ncount = 4;
    std::array<Long, ncount*num_fields> counts{};
    std::array<Real, 2*num_fields> range{};
    for (int f = 0; f < num_fields; ++f) {
        counts[ncount*f  ] = pend[f].n_filled;
        counts[ncount*f+1] = pend[f].n_outside;
        counts[ncount*f+2] = pend[f].n_not_finite;
        counts[ncount*f+3] = pend[f].n_uncovered;
        range[2*f  ] = -pend[f].vmin;
        range[2*f+1] =  pend[f].vmax;
    }
    ParallelDescriptor::ReduceLongSum(counts.data(), static_cast<int>(counts.size()));
    ParallelDescriptor::ReduceRealMax(range.data(), static_cast<int>(range.size()));

    for (const auto& d : catalog()) {
        const int f = static_cast<int>(d.id);
        if (d.kind == Kind::LayerColumn || !pend[f].read) { continue; }
        auto& mf = m_fields[lev][f];
        AMREX_ALWAYS_ASSERT(mf != nullptr);
        const int  ncomp     = mf->nComp();
        const Long n_expect  = mf->boxArray().numPts() * ncomp;
        const Long n_filled  = counts[ncount*f];
        const Long n_outside = counts[ncount*f+1];
        const Long n_nonfin  = counts[ncount*f+2];
        const Long n_uncov   = counts[ncount*f+3];
        const Real vmin = -range[2*f];
        const Real vmax =  range[2*f+1];

        // n_filled is a sum over the files. The files of a level do not overlap (the
        // reader takes them as the level's boxes). If two ever did, the columns they
        // share would be counted twice, n_filled would exceed n_expect, and the field
        // would be left out with a count above the level's.
        std::ostringstream why;
        if (pend[f].conflict) {
            why << "the files of the level disagree on its number of soil layers";
        } else if (n_uncov > 0 || n_filled != n_expect) {
            why << "the files give " << n_filled << " of the level's " << n_expect << " value(s)";
        } else if (n_outside > 0 || n_nonfin > 0) {
            why << n_outside << " value(s) outside " << d.lo << " to " << d.hi
                << " [" << d.units << "], " << n_nonfin << " not finite";
            if (vmin <= vmax) {
                why << "; the finite values span " << vmin << " to " << vmax;
            }
        }
        if (!why.str().empty()) {
            Print() << "Warning: wrfinput " << d.nc_name << " at level " << lev
                    << " is not held: " << why.str() << "\n";
            drop(lev, d.id);
            continue;
        }
        Print() << "   wrfinput surface field " << d.nc_name << " at level " << lev << ": "
                << ncomp << (ncomp == 1 ? " component" : " components")
                << ", min " << vmin << ", max " << vmax << " [" << d.units << "]\n";
        fill_ghost_cells(*mf, period);
    }

    // TSLB, SMOIS, ZS and DZS, whichever are held, must agree on the number of layers.
    // Every rank holds the same counts.
    {
        int nlay = 0;
        bool agree = true;
        const auto take = [&nlay, &agree] (int n) {
            if (n == 0) { return; }   // not held
            if (nlay == 0) { nlay = n; }
            agree = agree && (n == nlay);
        };
        take(static_cast<int>(m_soil_depth[lev].size()));
        take(static_cast<int>(m_soil_thickness[lev].size()));
        for (const Field id : {Field::TSLB, Field::SMOIS}) {
            const MultiFab* mf = get(lev, id);
            take(mf ? mf->nComp() : 0);
        }
        if (!agree) {
            Print() << "Warning: wrfinput TSLB, SMOIS, ZS and DZS at level " << lev
                    << " are not held: they do not have the same number of soil layers\n";
            for (const Field id : {Field::TSLB, Field::SMOIS, Field::ZS, Field::DZS}) {
                drop(lev, id);
            }
        }
    }

    pend = {};
}

void
WRFSurfaceFields::fill_ghost_cells (MultiFab& mf, const Periodicity& period)
{
    if (mf.nGrowVect() == IntVect(0)) { return; }
    const int ncomp = mf.nComp();
    // Not tiled: each box reads its own valid cells and writes its own ghost cells
    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const Dim3 lo = amrex::lbound(mfi.validbox());
        const Dim3 hi = amrex::ubound(mfi.validbox());
        const Array4<Real> a = mf.array(mfi);
        ParallelFor(mfi.fabbox(), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int ci = amrex::min(amrex::max(i, lo.x), hi.x);
            const int cj = amrex::min(amrex::max(j, lo.y), hi.y);
            if (ci != i || cj != j) { a(i, j, k, n) = a(ci, cj, k, n); }
        });
    }
    // The ghost cells another box (or a periodic image) covers take the value there
    mf.FillBoundary(period);

    // A ghost cell beyond the level's bounding box in a direction that is not periodic
    // takes the nearest column inside the bounding box. That column is a valid cell of
    // this box or one of the ghost cells set above; the cells read lie inside the
    // bounding box and the cells written outside it.
    const Box bound = mf.boxArray().minimalBox();
    const Dim3 blo = amrex::lbound(bound);
    const Dim3 bhi = amrex::ubound(bound);
    const bool px = period.isPeriodic(0);
    const bool py = period.isPeriodic(1);
    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const Array4<Real> a = mf.array(mfi);
        ParallelFor(mfi.fabbox(), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int ci = px ? i : amrex::min(amrex::max(i, blo.x), bhi.x);
            const int cj = py ? j : amrex::min(amrex::max(j, blo.y), bhi.y);
            if (ci != i || cj != j) { a(i, j, k, n) = a(ci, cj, k, n); }
        });
    }
}

std::unique_ptr<MultiFab>
WRFSurfaceFields::on_grids (int lev, Field id, const MultiFab* src,
                            const BoxArray& ba2d, const DistributionMapping& dm,
                            const IntVect& ng, const Periodicity& period,
                            const IntVect& ratio) const
{
    const MultiFab* crse = (lev > 0) ? get(lev-1, id) : nullptr;
    if (src == nullptr && crse == nullptr) { return nullptr; }
    const int ncomp = src ? src->nComp() : crse->nComp();

    const bool src_covers = src && src->boxArray().contains(ba2d);
    bool use_crse = false;
    if (!src_covers) {
        // The coarser level must hold the field with the same soil layers, on every
        // column under the new grids (proper nesting gives the last)
        if (crse == nullptr || crse->nComp() != ncomp) { return nullptr; }
        BoxArray cba = ba2d;
        cba.coarsen(IntVect(ratio[0], ratio[1], 1));
        if (!crse->boxArray().contains(cba)) { return nullptr; }
        use_crse = true;
    }

    auto fresh = std::make_unique<MultiFab>(ba2d, dm, ncomp, ng);
    if (use_crse) { fill_from_coarser(*fresh, *crse, ratio); }
    // The level's own values go on top, wherever they reach
    if (src) { fresh->ParallelCopy(*src, 0, 0, ncomp, IntVect(0), IntVect(0)); }
    fill_ghost_cells(*fresh, period);
    return fresh;
}

void
WRFSurfaceFields::remake (int lev,
                          const BoxArray& new_ba2d,
                          const DistributionMapping& new_dm,
                          const Periodicity& period,
                          const IntVect& ratio,
                          bool take_coarser_fields,
                          int verbose)
{
    if (lev < 0 || lev >= num_levels()) { return; }
    const Long ncol = new_ba2d.numPts();

    for (const auto& d : catalog()) {
        if (d.kind == Kind::LayerColumn) { continue; }
        auto& f = m_fields[lev][static_cast<int>(d.id)];
        const MultiFab* crse = (lev > 0) ? get(lev-1, d.id) : nullptr;
        const bool held = (f != nullptr);
        if (!held && !(take_coarser_fields && crse)) { continue; }
        if (held && f->boxArray() == new_ba2d && f->DistributionMap() == new_dm) { continue; }

        const IntVect ng = held ? f->nGrowVect() : crse->nGrowVect();
        // Columns the level's own grids give; the others come from the coarser level
        const Long n_own = (held && verbose > 0) ? amrex::intersect(new_ba2d, f->boxArray()).numPts() : 0;
        auto fresh = on_grids(lev, d.id, f.get(), new_ba2d, new_dm, ng, period, ratio);
        if (!fresh) {
            if (held) {
                Print() << "Warning: wrfinput " << d.nc_name << " at level " << lev
                        << " is dropped: the new grids reach beyond the grids that held it";
                if (lev > 0) {
                    Print() << ", and level " << lev-1
                            << " does not hold it with the same number of components";
                }
                Print() << "\n";
                f.reset();
            }
            continue;
        }
        f = std::move(fresh);
        if (verbose > 0) {
            Real vmin =  std::numeric_limits<Real>::max();
            Real vmax = -std::numeric_limits<Real>::max();
            for (int n = 0; n < f->nComp(); ++n) {
                vmin = amrex::min(vmin, f->min(n));
                vmax = amrex::max(vmax, f->max(n));
            }
            Print() << "   wrfinput surface field " << d.nc_name << " at level " << lev
                    << " on the new grids: " << f->nComp()
                    << (f->nComp() == 1 ? " component" : " components") << ", "
                    << ncol - n_own << " of " << ncol << " columns from level "
                    << amrex::max(lev-1, 0) << ", min " << vmin << ", max " << vmax
                    << " [" << d.units << "]\n";
        }
    }

    // The soil layer columns do not depend on the grids: they are kept. A level that
    // takes the fields of the coarser level takes its columns too, when it has none.
    if (take_coarser_fields && lev > 0) {
        if (m_soil_depth[lev].empty())     { m_soil_depth[lev]     = m_soil_depth[lev-1]; }
        if (m_soil_thickness[lev].empty()) { m_soil_thickness[lev] = m_soil_thickness[lev-1]; }
    }
}

void
WRFSurfaceFields::write_header (std::ostream& os, int nlevels) const
{
    os << chk_magic << " " << chk_version << "\n" << nlevels << "\n";
    os << std::setprecision(17);
    for (int lev = 0; lev < nlevels; ++lev) {
        std::ostringstream entries;
        entries << std::setprecision(17);
        int n = 0;
        if (lev < num_levels()) {
            for (const auto& d : catalog()) {
                if (d.kind == Kind::LayerColumn) {
                    const auto& v = (d.id == Field::DZS) ? m_soil_thickness[lev] : m_soil_depth[lev];
                    if (v.empty()) { continue; }
                    entries << d.nc_name << " " << v.size();
                    for (const Real x : v) { entries << " " << x; }
                    entries << "\n";
                    ++n;
                } else if (const MultiFab* mf = get(lev, d.id)) {
                    entries << d.nc_name << " " << mf->nComp() << "\n";
                    ++n;
                }
            }
        }
        os << "Level " << lev << " " << n << "\n" << entries.str();
    }
}

bool
WRFSurfaceFields::parse_header (std::istream& is, Vector<HeaderLevel>& out)
{
    out.clear();
    std::string magic;
    int version = 0;
    int nlevels = 0;
    if (!(is >> magic >> version >> nlevels) || magic != chk_magic ||
        version != chk_version || nlevels < 0) {
        return false;
    }
    Vector<HeaderLevel> levels(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        std::string word;
        int ilev = -1;
        int n = -1;
        if (!(is >> word >> ilev >> n) || word != "Level" || ilev != lev || n < 0) { return false; }
        for (int e = 0; e < n; ++e) {
            std::string name;
            int count = 0;
            if (!(is >> name >> count) || count < 1) { return false; }
            const Descriptor* d = find(name);
            if (d == nullptr) { return false; }
            if (d->kind == Kind::LayerColumn) {
                Vector<Real> v(count);
                for (auto& x : v) {
                    if (!(is >> x)) { return false; }
                }
                ((d->id == Field::DZS) ? levels[lev].soil_thickness : levels[lev].soil_depth) = v;
            }
            levels[lev].ncomp[static_cast<int>(d->id)] = count;
        }
    }
    out = std::move(levels);
    return true;
}

void
WRFSurfaceFields::write_checkpoint (const std::string& checkpoint_dir, int nlevels) const
{
    if (ParallelDescriptor::IOProcessor()) {
        const std::string name = checkpoint_dir + "/" + checkpoint_header_name;
        std::ofstream os(name, std::ofstream::out | std::ofstream::trunc);
        write_header(os, nlevels);
        if (!os.good()) { amrex::FileOpenFailed(name); }
    }
    for (int lev = 0; lev < nlevels && lev < num_levels(); ++lev) {
        for (const auto& d : catalog()) {
            if (const MultiFab* mf = get(lev, d.id)) {
                VisMF::Write(*mf, MultiFabFileFullPrefix(lev, checkpoint_dir, "Level_",
                                                         std::string(chk_mf_prefix) + d.nc_name));
            }
        }
    }
}

bool
WRFSurfaceFields::read_checkpoint (const std::string& checkpoint_dir, int nlevels,
                                   const Vector<BoxArray>& ba2d,
                                   const Vector<DistributionMapping>& dm,
                                   const IntVect& ng,
                                   const Vector<Periodicity>& period,
                                   const Vector<IntVect>& ratio)
{
    const std::string name = checkpoint_dir + "/" + checkpoint_header_name;
    if (!amrex::FileExists(name)) {
        Print() << "   " << name << " is not in the checkpoint (written before the wrfinput "
                << "land and soil fields were stored): none is restored\n";
        return false;
    }

    Vector<char> chars;
    ParallelDescriptor::ReadAndBcastFile(name, chars);
    std::istringstream is(std::string(chars.dataPtr()));
    Vector<HeaderLevel> levels;
    if (!parse_header(is, levels)) {
        Print() << "Warning: " << name << " could not be read: the wrfinput land and soil "
                << "fields are not restored\n";
        return false;
    }

    const int nlev = amrex::min(nlevels, num_levels(), static_cast<int>(levels.size()));
    for (int lev = 0; lev < nlev; ++lev) {
        for (auto& f : m_fields[lev]) { f.reset(); }
        m_pending[lev] = {};
        m_soil_depth[lev]     = levels[lev].soil_depth;
        m_soil_thickness[lev] = levels[lev].soil_thickness;
        for (const Field id : {Field::ZS, Field::DZS}) {
            const auto& v = (id == Field::DZS) ? m_soil_thickness[lev] : m_soil_depth[lev];
            if (v.empty()) { continue; }
            Print() << "   wrfinput soil layer "
                    << ((id == Field::DZS) ? "thicknesses DZS" : "depths ZS")
                    << " at level " << lev << " restored from the checkpoint:";
            for (const Real x : v) { Print() << " " << x; }
            Print() << " [m]\n";
        }
        for (const auto& d : catalog()) {
            const int ncomp = levels[lev].ncomp[static_cast<int>(d.id)];
            if (d.kind == Kind::LayerColumn || ncomp < 1) { continue; }
            const std::string prefix = MultiFabFileFullPrefix(lev, checkpoint_dir, "Level_",
                                                              std::string(chk_mf_prefix) + d.nc_name);
            if (!amrex::FileExists(prefix + "_H")) {
                Print() << "Warning: " << prefix << " is not in the checkpoint: wrfinput "
                        << d.nc_name << " at level " << lev << " is not restored\n";
                continue;
            }
            // Defined by VisMF::Read on the grids the checkpoint was written with
            MultiFab src;
            VisMF::Read(src, prefix);
            if (src.nComp() != ncomp) {
                Print() << "Warning: " << prefix << " holds " << src.nComp() << " component(s) and "
                        << name << " lists " << ncomp << ": wrfinput " << d.nc_name
                        << " at level " << lev << " is not restored\n";
                continue;
            }
            // As on a regrid: the checkpoint's values where its grids reach, the coarser
            // level (restored just before) elsewhere, and no field when neither covers
            // this run's grids. No cell is left without a value.
            auto& mf = m_fields[lev][static_cast<int>(d.id)];
            mf = on_grids(lev, d.id, &src, ba2d[lev], dm[lev], ng, period[lev],
                          (lev > 0) ? ratio[lev-1] : IntVect(1));
            if (!mf) {
                Print() << "Warning: wrfinput " << d.nc_name << " at level " << lev
                        << " is not restored: this run's grids reach beyond the grids of "
                        << "the checkpoint";
                if (lev > 0) {
                    Print() << ", and level " << lev-1
                            << " does not hold it with the same number of components";
                }
                Print() << "\n";
                continue;
            }
            // The range is reported so that a restart can be compared with the first run
            Real vmin =  std::numeric_limits<Real>::max();
            Real vmax = -std::numeric_limits<Real>::max();
            for (int n = 0; n < ncomp; ++n) {
                vmin = amrex::min(vmin, mf->min(n));
                vmax = amrex::max(vmax, mf->max(n));
            }
            Print() << "   wrfinput surface field " << d.nc_name << " at level " << lev
                    << " restored from the checkpoint: " << ncomp
                    << (ncomp == 1 ? " component" : " components")
                    << ", min " << vmin << ", max " << vmax << " [" << d.units << "]\n";
        }
    }
    return true;
}

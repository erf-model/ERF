/**
 * \file ERF_IBFaceSet.cpp
 * \brief Detection, storage, output and reporting of the wall faces of
 *        resolved buildings for the immersed-boundary surface energy
 *        balance. The conventions are documented in ERF_IBFaceSet.H.
 */
#include "ERF_IBFaceSet.H"
#include "ERF_IBSEBSolar.H"
#include "ERF_IBSEBSlab.H"
#include "ERF_IBSEBBalance.H"
#include <ERF_TwoStreamCanopyForcing.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_Constants.H>
#include <ERF_MOSTStress.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Gpu.H>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <AMReX_BoxList.H>
#include <AMReX_ParallelReduce.H>

using namespace amrex;

namespace {

/**
 * Host-readable view of one fab of a MultiFab.
 *
 * The face detection is a host loop at initialisation. On CPU builds the fab
 * is read in place; on GPU builds a pinned copy is made so the loop does not
 * touch device memory. The copy keeps the fab's own box, so the returned
 * Array4 is indexed with the same global (i, j, k) as the device fab.
 */
struct HostFab
{
#ifdef AMREX_USE_GPU
    FArrayBox fab;
    explicit HostFab (const FArrayBox& d) : fab(d.box(), d.nComp(), The_Pinned_Arena())
    {
        Gpu::dtoh_memcpy_async(fab.dataPtr(), d.dataPtr(), d.size() * sizeof(Real));
        Gpu::streamSynchronize();
    }
    Array4<const Real> array () const { return fab.const_array(); }
#else
    const FArrayBox* fab;
    explicit HostFab (const FArrayBox& d) : fab(&d) {}
    Array4<const Real> array () const { return fab->const_array(); }
#endif
};

/** Resize a device vector to a host vector and copy it over. */
template <class T>
void upload (Gpu::DeviceVector<T>& d, const std::vector<T>& h)
{
    d.resize(h.size());
    Gpu::copy(Gpu::hostToDevice, h.begin(), h.end(), d.begin());
}

/** Resize a device vector to n entries all equal to value. */
template <class T>
void fill (Gpu::DeviceVector<T>& d, size_t n, T value)
{
    std::vector<T> h(n, value);
    upload(d, h);
}

/**
 * Water vapour mixing ratio of cell (i, j, k): rho qv over rho, at least zero
 * (as the two-stream column reads it, get_qv_from_state()), or zero when the
 * run has no moisture model (``qv_comp`` negative, the state then has no such
 * component).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real cell_qv (const Array4<const Real>& c, int i, int j, int k, int qv_comp) noexcept
{
    return (qv_comp >= 0) ? amrex::max(Real(0.0), c(i, j, k, qv_comp) / c(i, j, k, Rho_comp)) : Real(0.0);
}

/** Host mirror of one fab of an iMultiFab: its whole box, ghost cells included. */
struct HostInts
{
    Box box;
    std::vector<int> v;
    HostInts (const Box& b, int value) : box(b), v(static_cast<size_t>(b.numPts()), value) {}
    explicit HostInts (const IArrayBox& fab) : box(fab.box()), v(static_cast<size_t>(fab.box().numPts()))
    {
        AMREX_ALWAYS_ASSERT(fab.nComp() == 1);
        Gpu::copy(Gpu::deviceToHost, fab.dataPtr(), fab.dataPtr() + v.size(), v.begin());
        Gpu::streamSynchronize();
    }
    Array4<int> array () noexcept { return Array4<int>(v.data(), amrex::begin(box), amrex::end(box), 1); }
    Array4<const int> array () const noexcept { return Array4<const int>(v.data(), amrex::begin(box), amrex::end(box), 1); }
    void upload (IArrayBox& fab) const
    {
        AMREX_ALWAYS_ASSERT(fab.box() == box);
        Gpu::copy(Gpu::hostToDevice, v.begin(), v.end(), fab.dataPtr());
        Gpu::streamSynchronize();
    }
};

/** Every rank's entries of a vector, in rank order (MPI_Allgatherv; the vector itself without MPI). */
template <class T>
std::vector<T> all_gather (const std::vector<T>& mine)
{
#ifdef AMREX_USE_MPI
    const int np = ParallelDescriptor::NProcs();
    std::vector<int> counts(np, 0), displs(np, 0);
    const int n = static_cast<int>(mine.size());
    BL_MPI_REQUIRE(MPI_Allgather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, ParallelDescriptor::Communicator()));
    Long total = 0;
    for (int r = 0; r < np; ++r) { displs[r] = static_cast<int>(total); total += counts[r]; }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(total < Long(std::numeric_limits<int>::max()), "erf.ibseb: too many entries to gather");
    std::vector<T> all(static_cast<size_t>(total));
    BL_MPI_REQUIRE(MPI_Allgatherv(mine.data(), n, ParallelDescriptor::Mpi_typemap<T>::type(),
                                  all.data(), counts.data(), displs.data(), ParallelDescriptor::Mpi_typemap<T>::type(),
                                  ParallelDescriptor::Communicator()));
    return all;
#else
    return mine;
#endif
}

/**
 * The level's boxes with their k range collapsed to the domain's bottom
 * cell: the layout of the column maps, one cell per column. No two boxes
 * may hold the same column, so the grids must not be split in z
 * (amr.no_box_split_dir, 2 by default in ERF, keeps them whole).
 */
BoxArray column_box_array (const BoxArray& ba, const Box& domain)
{
    BoxList bl = ba.boxList();
    for (Box& b : bl) { b.setRange(2, domain.smallEnd(2), 1); }
    BoxArray cba(std::move(bl));
    std::vector<std::pair<int, Box>> isects;
    for (int n = 0; n < static_cast<int>(cba.size()); ++n) {
        cba.intersections(cba[n], isects);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(isects.size() == 1,
            "erf.ibseb: the level's grids are split in z, so a column spans several boxes; keep amr.no_box_split_dir = 2 (the default)");
    }
    return cba;
}

} // namespace

/**
 * Detect the faces owned by this rank and allocate the per-face arrays.
 *
 * The scan visits every valid cell of every local fab in MFIter order, once.
 * A fluid cell (blanking < 0.5) contributes one face for each of its six
 * neighbours that is solid (blanking >= 0.5); the neighbour may sit in a
 * ghost cell, which is why the blanking must have its ghost cells filled.
 * Neighbours outside the domain in a non-periodic direction are skipped, so
 * a building against a non-periodic boundary has no face there. Solid cells
 * contribute no face, only the column map (the highest solid cell of each
 * column of the level's own boxes) that the building ids and the ray cast
 * read. A solid neighbour outside the level's grids, of a fluid or a solid
 * cell, means a building crosses the edge of a refined level; the scan
 * counts them and stops the run once all ranks have reported.
 *
 * Faces are appended in scan order, which makes them contiguous per fab;
 * m_fab_start records where each fab's faces begin. The per-direction,
 * per-building and total counts are reduced over all ranks here, once, so
 * report() does not have to reduce static numbers every time. Nothing of
 * the size of the domain is gathered: the column map is exchanged through
 * its ghost cells, the building labels through theirs (see the header), and
 * the one gather is of the distinct labels, one per building.
 */
void
IBFaceSet::build (const MultiFab& blanking, const Geometry& geom)
{
    const Box& domain = geom.Domain();
    const int nx  = domain.length(0);
    const int ny  = domain.length(1);
    const int ilo = domain.smallEnd(0);
    const int jlo = domain.smallEnd(1);
    const int klo = domain.smallEnd(2);
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const Real face_area[3] = { dx[1] * dx[2], dx[0] * dx[2], dx[0] * dx[1] };
    m_geom = geom;
    m_nx = nx; m_ny = ny; m_nz = domain.length(2);
    m_dom_lo = domain.smallEnd();
    m_x_lo = plo[0]; m_y_lo = plo[1];
    for (int d = 0; d < 3; ++d) { m_dx[d] = dx[d]; }
    m_per_x = geom.isPeriodic(0); m_per_y = geom.isPeriodic(1);
    m_z_ground = plo[2];

    // The face list of a cell must not depend on the decomposition: every
    // valid cell has to see its six neighbours, which needs a filled ghost
    // layer on the blanking.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(blanking.nGrowVect().min() >= 1,
                                     "erf.ibseb: the blanking needs at least one ghost cell");
    // The labelling below gives every column its scan index plus one as a label, in an int.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Long(nx) * ny < Long(std::numeric_limits<int>::max()),
                                     "erf.ibseb: the level has more columns than an int can label");

    // The level's own column map: one cell per column of each of the level's
    // boxes, with one ghost cell for the neighbours the faces and the
    // labelling read across boxes and periodic seams.
    const BoxArray col_ba = column_box_array(blanking.boxArray(), domain);
    m_col_valid.define(col_ba, blanking.DistributionMap(), 1, IntVect(AMREX_D_DECL(1, 1, 0)));
    m_col_valid.setVal(-1);

    std::vector<int>  h_i, h_j, h_k, h_dir, h_side, h_nbi, h_nbj;
    std::vector<Real> h_area, h_xf, h_yf, h_zf;
    // Solid neighbours outside the level's grids, of a fluid cell (a wall
    // against the edge) or of a solid one (a building running across the
    // edge, which a low roof can do with no fluid cell in between): a
    // building that crosses the edge of a refined level, whose columns there
    // this level cannot see.
    const BoxArray& level_ba = blanking.boxArray();
    Long n_outside = 0;
    Long first_outside = std::numeric_limits<Long>::max();   // lowest such cell, as a linear domain index
    auto solid_outside = [&] (const Array4<const Real>& b, const Box& vbox, const Box& fbox, const IntVect& nb) {
        // A neighbour in the same valid box is in the level (the common case,
        // and every case on level 0): no lookup.
        if (vbox.contains(nb) || !fbox.contains(nb) || b(nb) < 0.5) { return false; }
        IntVect nbw = nb;
        for (int dp = 0; dp < 2; ++dp) {
            if (!geom.isPeriodic(dp)) { continue; }
            const int lo_d = domain.smallEnd(dp), len_d = domain.length(dp);
            nbw[dp] = lo_d + ((nbw[dp] - lo_d) % len_d + len_d) % len_d;
        }
        if (level_ba.contains(nbw)) { return false; }
        const IntVect r = nbw - domain.smallEnd();
        first_outside = std::min(first_outside,
                                 (Long(r[0]) * domain.length(1) + r[1]) * domain.length(2) + r[2]);
        ++n_outside;
        return true;
    };

    m_fab_start.clear();
    for (MFIter mfi(blanking); mfi.isValid(); ++mfi) {
        m_fab_start.push_back(static_cast<int>(h_i.size()));
        const Box& bx = mfi.validbox();
        const HostFab hf(blanking[mfi]);
        auto const& b = hf.array();
        const Box& fbox = Box(b);
        // The column map of this box: the highest solid cell of each of its
        // columns, as a domain-relative k index, -1 where the column carries
        // no solid cell. It is an index rather than a height so that it costs
        // 4 bytes per column.
        IArrayBox& cfab = m_col_valid[mfi];
        HostInts top(cfab.box(), -1);
        auto const& t = top.array();
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if (b(i, j, k) >= 0.5) {
                t(i, j, klo) = std::max(t(i, j, klo), k - klo);
                for (int d = 0; d < 3; ++d) {
                    for (int s = -1; s <= 1; s += 2) {
                        IntVect nb(i, j, k);
                        nb[d] += s;
                        if (!geom.isPeriodic(d) && !domain.contains(nb)) { continue; }
                        solid_outside(b, bx, fbox, nb);
                    }
                }
                continue;
            }
            for (int d = 0; d < 3; ++d) {
                for (int s = -1; s <= 1; s += 2) {
                    IntVect nb(i, j, k);
                    nb[d] += s;
                    if (!fbox.contains(nb)) { continue; }
                    if (!geom.isPeriodic(d) && !domain.contains(nb)) { continue; }
                    if (b(nb) < 0.5) { continue; }
                    solid_outside(b, bx, fbox, nb);
                    h_i.push_back(i); h_j.push_back(j); h_k.push_back(k);
                    h_dir.push_back(d); h_side.push_back(s);
                    // Face centre: the fluid cell centre moved half a cell toward the solid
                    // (plo is the position of the domain's low corner cell).
                    Real c[3] = { Real(plo[0] + (i - ilo + 0.5) * dx[0]), Real(plo[1] + (j - jlo + 0.5) * dx[1]), Real(plo[2] + (k - klo + 0.5) * dx[2]) };
                    c[d] += 0.5 * s * dx[d];
                    h_xf.push_back(c[0]); h_yf.push_back(c[1]); h_zf.push_back(c[2]);
                    // Column of the solid neighbour, unwrapped: it lies in this
                    // box or in its ghost cell, where the periodic image is.
                    h_nbi.push_back(nb[0]); h_nbj.push_back(nb[1]);
                    h_area.push_back(face_area[d]);
                }
            }
        }}}
        top.upload(cfab);
    }
    m_fab_start.push_back(static_cast<int>(h_i.size()));
    m_nface = static_cast<int>(h_i.size());
    m_col_valid.FillBoundary(geom.periodicity());

    // The column map, the building ids and the ray cast below see only the
    // level's own cells, so a building must lie wholly inside the level.
    // ERF::init_ibseb() checks that for the buildings of the level below; this
    // catches one only the refined level resolves.
    {
        Long n_out_all = n_outside;
        Long first_all = first_outside;
        ParallelDescriptor::ReduceLongSum(n_out_all);
        ParallelDescriptor::ReduceLongMin(first_all);
        if (n_out_all > 0) {
            const Long nyz = Long(domain.length(1)) * domain.length(2);
            const int fi = static_cast<int>(first_all / nyz);
            const int fj = static_cast<int>((first_all / domain.length(2)) % domain.length(1));
            const int fk = static_cast<int>(first_all % domain.length(2));
            std::ostringstream os;
            os << "erf.ibseb: level " << m_lev << " has " << n_out_all
               << " cell faces against solid cells outside its grids (the first solid one at x = "
               << plo[0] + (fi + 0.5) * dx[0] << " m, y = " << plo[1] + (fj + 0.5) * dx[1]
               << " m, z = " << plo[2] + (fk + 0.5) * dx[2] << " m)";
            os << ": a building crosses the edge of the refined level. To fix the refined box, run"
                  " `python3 Exec/CanonicalTests/SEB/ibseb_refinement_box.py <inputs>`: it checks the deck's"
                  " erf.<name>.in_box_lo / in_box_hi against the height map and prints a box that no building"
                  " crosses (or `... <inputs> --all --fit tight|relaxed` for a box around every building).";
            Abort(os.str());
        }
    }

    // The highest interface any face samples (two_stream_canopy_sample(): a roof its
    // fluid cell's bottom, a ceiling and a wall the top of its cell too), counted from the domain
    // bottom; -1 when the level has no faces. The two-stream beam the faces read needs
    // no interface above it.
    m_top_sample = -1;
    for (int n = 0; n < m_nface; ++n) {
        m_top_sample = std::max(m_top_sample, h_k[n] - domain.smallEnd(2) + ((h_dir[n] != 2 || h_side[n] > 0) ? 1 : 0));
    }
    ParallelDescriptor::ReduceIntMax(m_top_sample);

    // Tallest column of the level over all ranks: the height a rising ray
    // stops at (one cell above the highest solid cell; the ground when the
    // level holds no solid cell).
    {
        const int kmax = m_col_valid.max(0);
        m_col_top_max = (kmax >= 0) ? Real(m_z_ground + (kmax + 1.0) * dx[2]) : Real(m_z_ground);
    }

    // Host copies of the column map, ghost cells filled: the column of a
    // face's solid neighbour is in its fab's box or its ghost ring.
    std::vector<HostInts> tops;
    std::vector<Box> vboxes;
    for (MFIter mfi(m_col_valid); mfi.isValid(); ++mfi) {
        tops.emplace_back(m_col_valid[mfi]);
        vboxes.push_back(mfi.validbox());
    }
    const int nfab = static_cast<int>(tops.size());
    AMREX_ALWAYS_ASSERT(nfab + 1 == static_cast<int>(m_fab_start.size()));
    // Height of each face's building column, for the wall-function depth.
    {
        std::vector<Real> h_hb(h_i.size());
        for (int li = 0; li < nfab; ++li) {
            auto const& t = tops[li].array();
            for (int n = m_fab_start[li]; n < m_fab_start[li + 1]; ++n) {
                const int kt = t(h_nbi[n], h_nbj[n], klo);
                h_hb[n] = (kt >= 0) ? (kt + 1.0) * dx[2] : 0.0;
            }
        }
        upload(d_hbld, h_hb);
    }

    // ---- Buildings: 4-connected solid columns, numbered in scan order ----
    // Each box splits its solid columns into connected pieces (a depth-first
    // flood fill); a piece's label starts as the lowest scan index plus one
    // of its columns, and the pieces then take the lowest label of their
    // neighbours across the box edges, through one ghost cell of m_lab,
    // until no label changes anywhere: the pieces of one building all end
    // with the scan index of its first column. The rounds are as many as a
    // building's pieces are deep, so the exchange never depends on the
    // number of columns.
    m_lab.define(col_ba, blanking.DistributionMap(), 1, IntVect(AMREX_D_DECL(1, 1, 0)));
    m_lab.setVal(0);
    std::vector<std::vector<int>> piece_of(nfab);    // per fab, per valid column: its piece, -1 where fluid
    std::vector<std::vector<int>> piece_label(nfab); // per fab, per piece: its label
    auto scan_label = [=] (int i, int j) { return (i - ilo) * ny + (j - jlo) + 1; };
    for (int li = 0; li < nfab; ++li) {
        const Box& vb = vboxes[li];
        auto const& t = tops[li].array();
        const auto lo = lbound(vb);
        const int lx = vb.length(0), ly = vb.length(1);
        auto idx = [=] (int i, int j) { return (j - lo.y) * lx + (i - lo.x); };
        std::vector<int>& piece = piece_of[li];
        piece.assign(static_cast<size_t>(lx) * ly, -1);
        std::vector<int> stack;
        for (int i = lo.x; i < lo.x + lx; ++i) {
            for (int j = lo.y; j < lo.y + ly; ++j) {
                if (t(i, j, klo) < 0 || piece[idx(i, j)] >= 0) { continue; }
                const int p = static_cast<int>(piece_label[li].size());
                int lab = scan_label(i, j);
                piece[idx(i, j)] = p;
                stack.clear();
                stack.push_back(idx(i, j));
                while (!stack.empty()) {
                    const int q = stack.back(); stack.pop_back();
                    const int qi = lo.x + q % lx, qj = lo.y + q / lx;
                    lab = std::min(lab, scan_label(qi, qj));
                    const int ni[4] = {qi - 1, qi + 1, qi, qi};
                    const int nj[4] = {qj, qj, qj - 1, qj + 1};
                    for (int n = 0; n < 4; ++n) {
                        if (!vb.contains(IntVect(ni[n], nj[n], klo))) { continue; }
                        const int r = idx(ni[n], nj[n]);
                        if (t(ni[n], nj[n], klo) >= 0 && piece[r] < 0) {
                            piece[r] = p;
                            stack.push_back(r);
                        }
                    }
                }
                piece_label[li].push_back(lab);
            }
        }
    }
    // Write the labels into the map, exchange, lower, repeat.
    auto write_labels = [&] (const std::vector<std::vector<int>>& value) {
        int li = 0;
        for (MFIter mfi(m_lab); mfi.isValid(); ++mfi, ++li) {
            IArrayBox& fab = m_lab[mfi];
            HostInts L(fab.box(), 0);
            auto const& a = L.array();
            const Box& vb = vboxes[li];
            const auto lo = lbound(vb);
            const int lx = vb.length(0), ly = vb.length(1);
            for (int j = lo.y; j < lo.y + ly; ++j) {
                for (int i = lo.x; i < lo.x + lx; ++i) {
                    const int p = piece_of[li][(j - lo.y) * lx + (i - lo.x)];
                    if (p >= 0) { a(i, j, klo) = value[li][p]; }
                }
            }
            L.upload(fab);
        }
        m_lab.FillBoundary(geom.periodicity());
    };
    int n_rounds = 0;
    for (;;) {
        write_labels(piece_label);
        ++n_rounds;
        bool changed = false;
        int li = 0;
        for (MFIter mfi(m_lab); mfi.isValid(); ++mfi, ++li) {
            const HostInts L(m_lab[mfi]);
            auto const& a = L.array();
            const Box& vb = vboxes[li];
            const auto lo = lbound(vb);
            const auto hi = ubound(vb);
            const int lx = vb.length(0);
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (i != lo.x && i != hi.x && j != lo.y && j != hi.y) { continue; }   // an inner column has no neighbour outside the box
                    const int p = piece_of[li][(j - lo.y) * lx + (i - lo.x)];
                    if (p < 0) { continue; }
                    const int ni[4] = {i - 1, i + 1, i, i};
                    const int nj[4] = {j, j, j - 1, j + 1};
                    for (int n = 0; n < 4; ++n) {
                        const IntVect nb(ni[n], nj[n], klo);
                        if (vb.contains(nb) || !L.box.contains(nb)) { continue; }
                        const int g = a(nb);
                        if (g > 0 && g < piece_label[li][p]) { piece_label[li][p] = g; changed = true; }
                    }
                }
            }
        }
        ParallelDescriptor::ReduceBoolOr(changed);
        if (!changed) { break; }
    }
    // The distinct labels, one per building, gathered from every rank and
    // numbered in scan order; the pieces take the number of their label.
    std::vector<int> labels_all;
    {
        std::set<int> mine;
        for (const auto& pl : piece_label) { mine.insert(pl.begin(), pl.end()); }
        const std::vector<int> labels_mine(mine.begin(), mine.end());
        labels_all = all_gather(labels_mine);
        std::sort(labels_all.begin(), labels_all.end());
        labels_all.erase(std::unique(labels_all.begin(), labels_all.end()), labels_all.end());
    }
    m_nbld = static_cast<int>(labels_all.size());
    std::vector<std::vector<int>> piece_id(nfab);
    for (int li = 0; li < nfab; ++li) {
        for (int lab : piece_label[li]) {
            const auto it = std::lower_bound(labels_all.begin(), labels_all.end(), lab);
            AMREX_ALWAYS_ASSERT(it != labels_all.end() && *it == lab);
            piece_id[li].push_back(static_cast<int>(it - labels_all.begin()) + 1);
        }
    }
    write_labels(piece_id);
    // The id of each face: that of the column of its solid neighbour, in the
    // fab's box or its ghost ring.
    std::vector<int> h_bid(m_nface);
    {
        int li = 0;
        for (MFIter mfi(m_lab); mfi.isValid(); ++mfi, ++li) {
            const HostInts L(m_lab[mfi]);
            auto const& a = L.array();
            for (int n = m_fab_start[li]; n < m_fab_start[li + 1]; ++n) {
                h_bid[n] = a(h_nbi[n], h_nbj[n], klo);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(h_bid[n] >= 1, "erf.ibseb: a face's solid neighbour has no building id");
            }
        }
    }
    // Footprints, for the debug summary: columns and bounding box per building.
    m_bld_ncol.assign(m_nbld + 1, 0);
    m_bld_ilo.assign(m_nbld + 1, std::numeric_limits<int>::max()); m_bld_ihi.assign(m_nbld + 1, std::numeric_limits<int>::min());
    m_bld_jlo.assign(m_nbld + 1, std::numeric_limits<int>::max()); m_bld_jhi.assign(m_nbld + 1, std::numeric_limits<int>::min());
    for (int li = 0; li < nfab; ++li) {
        const Box& vb = vboxes[li];
        const auto lo = lbound(vb);
        const int lx = vb.length(0), ly = vb.length(1);
        for (int j = lo.y; j < lo.y + ly; ++j) {
            for (int i = lo.x; i < lo.x + lx; ++i) {
                const int p = piece_of[li][(j - lo.y) * lx + (i - lo.x)];
                if (p < 0) { continue; }
                const int bb = piece_id[li][p];
                m_bld_ncol[bb] += 1;
                m_bld_ilo[bb] = std::min(m_bld_ilo[bb], i); m_bld_ihi[bb] = std::max(m_bld_ihi[bb], i);
                m_bld_jlo[bb] = std::min(m_bld_jlo[bb], j); m_bld_jhi[bb] = std::max(m_bld_jhi[bb], j);
            }
        }
    }
    ParallelDescriptor::ReduceLongSum(m_bld_ncol.data(), m_nbld + 1);
    ParallelDescriptor::ReduceIntMin(m_bld_ilo.data(), m_nbld + 1); ParallelDescriptor::ReduceIntMax(m_bld_ihi.data(), m_nbld + 1);
    ParallelDescriptor::ReduceIntMin(m_bld_jlo.data(), m_nbld + 1); ParallelDescriptor::ReduceIntMax(m_bld_jhi.data(), m_nbld + 1);
    // The column labels, by which a refined level finds its buildings in the
    // level below (map_buildings_to_level0()), are kept until then.
    m_to_lev0.resize(m_nbld + 1);
    for (int b = 0; b <= m_nbld; ++b) { m_to_lev0[b] = b; }
    tops.clear();
    piece_of.clear();

    // The ray map of this rank, from the level's own columns; a refined
    // level rebuilds it with the levels below (add_outside_occluders()).
    build_ray_map({});
    if (m_params.debug) {
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " buildings labelled in " << n_rounds << " exchange round(s)\n";
    }
    if (m_lev == 0) { warn_fixed_sun_below_floor(); }

    // Static totals: per direction, per building, and the summed area.
    std::vector<Long> nd(3, 0);
    std::vector<Long> bn(m_nbld + 1, 0);
    std::vector<Real> ba(m_nbld + 1, 0.0);
    Real area = 0.0;
    for (int n = 0; n < m_nface; ++n) {
        nd[h_dir[n]] += 1;
        bn[h_bid[n]] += 1;
        ba[h_bid[n]] += h_area[n];
        area += h_area[n];
    }
    ParallelDescriptor::ReduceLongSum(nd.data(), 3);
    ParallelDescriptor::ReduceLongSum(bn.data(), m_nbld + 1);
    ParallelDescriptor::ReduceRealSum(ba.data(), m_nbld + 1);
    ParallelDescriptor::ReduceRealSum(area);
    for (int d = 0; d < 3; ++d) { m_nface_dir[d] = nd[d]; }
    m_bld_nface = bn;
    m_bld_area  = ba;
    m_area_total = area;

    // Device arrays. The state starts uniform at T_skin_init (the slab too;
    // the interior boundary is applied when the slab is advanced); the view
    // fractions and fluxes start at zero and are filled by
    // compute_view_fractions() and the per-step routines.
    const size_t nf = static_cast<size_t>(m_nface);
    upload(d_i, h_i); upload(d_j, h_j); upload(d_k, h_k);
    upload(d_dir, h_dir); upload(d_side, h_side); upload(d_bid, h_bid);
    upload(d_area, h_area);
    upload(d_xf, h_xf); upload(d_yf, h_yf); upload(d_zf, h_zf);
    fill(d_mat,      nf, 0);
    fill(d_albedo,   nf, m_params.albedo);
    fill(d_emis,     nf, m_params.emissivity);
    fill(d_kth,      nf, m_params.k_therm);
    fill(d_rhocp,    nf, m_params.rho_cp);
    fill(d_thick,    nf, m_params.thickness);
    fill(d_T_skin,   nf, m_params.T_skin_init);
    fill(d_T_slab,   nf * static_cast<size_t>(n_layers()), m_params.T_skin_init);
    // View fractions: placeholders (a roof sees the whole sky, a wall half sky
    // and half ground) that compute_view_fractions() replaces at initialisation.
    {
        std::vector<Real> fs(nf), fg(nf), fb(nf, 0.0);
        for (size_t n = 0; n < nf; ++n) {
            const bool roof = (h_dir[n] == 2);
            fs[n] = roof ? 1.0 : 0.5;
            fg[n] = roof ? 0.0 : 0.5;
        }
        upload(d_f_sky, fs); upload(d_f_ground, fg); upload(d_f_bldg, fb);
    }
    fill(d_shadow,        nf, Real(0.0));
    fill(d_SW_direct_in,  nf, Real(0.0));
    fill(d_SW_diffuse_in, nf, Real(0.0));
    fill(d_LW_down_in,    nf, Real(0.0));
    fill(d_T_air,         nf, m_params.T_skin_init);
    fill(d_theta_air,     nf, m_params.T_skin_init);
    fill(d_rho,           nf, Real(1.0));
    fill(d_U_tan,         nf, Real(0.0));
    fill(d_ustar,         nf, Real(0.0));
    fill(d_SW_abs,   nf, Real(0.0));
    fill(d_LW_net,   nf, Real(0.0));
    fill(d_H,        nf, Real(0.0));
    fill(d_LE,       nf, Real(0.0));
    fill(d_G,        nf, Real(0.0));
    fill(d_Q_ext,    nf, m_params.Q_ext_uniform);
    fill(d_H_coeff,  nf, Real(0.0));
    fill(d_w_star,   nf, Real(0.0));
    fill(d_olen,     nf, Real(1.0e30));
    fill(d_z_i,      nf, Real(0.0));
    fill(d_LW_ext,   nf, Real(0.0));
    fill(d_resid,    nf, Real(0.0));
    fill(d_niter,    nf, 0);
    Gpu::streamSynchronize();

    if (m_params.debug) { print_debug_summary(); }
}

/**
 * The ray map of this rank (see the header): the halo, this rank's patch,
 * every rank's patch as one box array, and its contents from the levels
 * below (coarsest first, so the finest level that covers a column has the
 * last word) and then this level's own columns with their periodic images.
 */
void
IBFaceSet::build_ray_map (const std::vector<const IBFaceSet*>& below)
{
    const Box& domain = m_geom.Domain();
    // The tallest column of this level and of the levels below: the height
    // a rising ray stops at, and what sizes the halo.
    Real z_max = m_col_top_max;
    for (const IBFaceSet* c : below) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(c->m_col_valid.ok(), "erf.ibseb: a level below has not been built");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(c->m_geom.Domain().smallEnd(2) == domain.smallEnd(2) && c->m_z_ground == m_z_ground,
                                         "erf.ibseb: the levels do not share the ground");
        z_max = std::max(z_max, c->m_col_top_max);
    }
    m_col_top_max = z_max;
    // The shallowest ray the map follows to its end: the hemisphere's and the
    // sun's floor; the halo is its horizontal reach over the tallest column.
    const Real sz_view = ibseb::hemisphere_min_sz(m_params.view_n_az, m_params.view_n_el);
    const Real sz_sun  = std::sin(m_params.sun_min_elev_deg * PI / 180.0);
    const Real sz_min  = std::min(sz_view, sz_sun);
    const Real height  = z_max - m_z_ground;
    const int halo[2] = { ibseb::ray_halo_cells(height, m_dx[0], sz_min), ibseb::ray_halo_cells(height, m_dx[1], sz_min) };
    // The columns a ray at the floor can cross in x plus in y: what a walk is
    // followed for, however the patch is cut or wrapped.
    m_max_steps = static_cast<int>(std::min(Long(std::numeric_limits<int>::max() / 2), Long(halo[0]) + Long(halo[1])));

    // This rank's patch: its boxes grown by the halo. A halo wider than the
    // domain adds nothing beyond one domain width: where the direction is
    // periodic the patch is then the whole domain, wrapped; where it is not
    // the patch stops at the domain's edge.
    Box bbox;
    for (MFIter mfi(m_col_valid); mfi.isValid(); ++mfi) {
        if (bbox.ok()) { bbox.minBox(mfi.validbox()); } else { bbox = mfi.validbox(); }
    }
    Box patch;
    m_wrap_x = m_wrap_y = false;
    bool spans = true;
    for (int d = 0; d < 2; ++d) { m_halo[d] = std::min(halo[d], domain.length(d)); }
    m_halo[2] = 0;
    if (bbox.ok()) {
        patch = bbox;
        for (int d = 0; d < 2; ++d) {
            patch.grow(d, m_halo[d]);
            if (m_geom.isPeriodic(d)) {
                if (patch.length(d) >= domain.length(d)) {
                    patch.setRange(d, domain.smallEnd(d), domain.length(d));
                    (d == 0 ? m_wrap_x : m_wrap_y) = true;
                }
            } else {
                patch.setSmall(d, std::max(patch.smallEnd(d), domain.smallEnd(d)));
                patch.setBig(d, std::min(patch.bigEnd(d), domain.bigEnd(d)));
            }
            spans = spans && (patch.smallEnd(d) == domain.smallEnd(d) && patch.length(d) == domain.length(d));
        }
    }
    m_patch_box = patch;
    ParallelDescriptor::ReduceBoolAnd(spans);
    m_map_spans_domain = spans;

    // Every rank's patch, one box each for the ranks that have boxes.
    std::vector<int> mine;
    if (patch.ok()) {
        mine = { patch.smallEnd(0), patch.smallEnd(1), patch.smallEnd(2),
                 patch.bigEnd(0), patch.bigEnd(1), patch.bigEnd(2), ParallelDescriptor::MyProc() };
    }
    const std::vector<int> all = all_gather(mine);
    BoxList bl;
    Vector<int> pmap;
    for (size_t n = 0; n + 7 <= all.size(); n += 7) {
        bl.push_back(Box(IntVect(all[n], all[n + 1], all[n + 2]), IntVect(all[n + 3], all[n + 4], all[n + 5])));
        pmap.push_back(all[n + 6]);
    }
    m_ray_patch.clear();
    if (bl.isEmpty()) { return; }
    const BoxArray patch_ba(std::move(bl));
    const DistributionMapping patch_dm(pmap);
    m_ray_patch.define(patch_ba, patch_dm, 1, 0);
    m_ray_patch.setVal(-1);

    // The levels below: each column of the patch that a level covers takes
    // that level's top, on the ratio x ratio columns of this level beneath
    // the coarse one, at the same height; the finest level wins.
    for (const IBFaceSet* c : below) {
        const Box& cdom = c->m_geom.Domain();
        const IntVect ratio(AMREX_D_DECL(domain.length(0) / cdom.length(0), domain.length(1) / cdom.length(1), 1));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ratio[0] * cdom.length(0) == domain.length(0) && ratio[1] * cdom.length(1) == domain.length(1),
                                         "erf.ibseb: the domains of the levels are not refinements of one another");
        iMultiFab crse(amrex::coarsen(patch_ba, ratio), patch_dm, 1, 0);
        crse.setVal(-2);                                   // -2: that level does not cover the column
        crse.ParallelCopy(c->m_col_valid, 0, 0, 1, IntVect(0), IntVect(0), c->m_geom.periodicity());
        const Real dz_c = c->m_dx[2], dz_f = m_dx[2];
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(dz_c / dz_f - std::round(dz_c / dz_f)) < 1.0e-6,
                                         "erf.ibseb: the vertical cell sizes of the levels are not a whole ratio");
        const int rx = ratio[0], ry = ratio[1];
        for (MFIter mfi(m_ray_patch); mfi.isValid(); ++mfi) {
            const Box& pb = mfi.validbox();
            auto const& p  = m_ray_patch.array(mfi);
            auto const& cr = crse.const_array(mfi);
            ParallelFor(pb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const int kc = cr(amrex::coarsen(i, rx), amrex::coarsen(j, ry), k);
                if (kc == -2) { return; }
                // The same height here: the top is kc + 1 cells of that level
                // up, a whole number of cells of this one (rounded to keep it
                // so whatever the precision of the division).
                p(i, j, k) = (kc < 0) ? -1 : static_cast<int>((kc + 1) * dz_c / dz_f + Real(0.5)) - 1;
            });
        }
    }
    // This level's own columns, with their periodic images.
    m_ray_patch.ParallelCopy(m_col_valid, 0, 0, 1, IntVect(0), IntVect(0), m_geom.periodicity());
    Gpu::streamSynchronize();

    if (m_params.debug) {
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " ray map: halo " << m_halo[0] << " x " << m_halo[1]
                << " cells for rays " << std::asin(sz_min) * 180.0 / PI << " deg above the horizon over "
                << height << " m (hemisphere " << std::asin(sz_view) * 180.0 / PI << " deg, sun floor "
                << m_params.sun_min_elev_deg << " deg); " << (spans ? "every rank's patch spans the domain" : "the patches do not span the domain")
                << "; " << below.size() << " level(s) below merged\n";
        std::ostringstream os;
        os << "[IBSEB DEBUG] lev=" << m_lev << " rank " << ParallelDescriptor::MyProc() << ": ray map patch " << patch
           << (m_wrap_x ? " wraps in x" : "") << (m_wrap_y ? " wraps in y" : "")
           << ", " << (patch.ok() ? patch.numPts() * sizeof(int) / 1024.0 : 0.0) << " kB\n";
        AllPrint() << os.str();
    }
}

/**
 * Warn when a fixed sun stands below ``erf.ibseb.sun_min_elev_deg`` and
 * the ray map does not span the domain: the shadows of buildings beyond the
 * halo are then not cast. Printed once per level by the I/O rank.
 */
void
IBFaceSet::warn_fixed_sun_below_floor () const
{
    if (m_params.sun_mode != "fixed" || m_map_spans_domain) { return; }
    const Real elev = 90.0 - m_params.sun_zenith_deg;
    if (elev >= m_params.sun_min_elev_deg) { return; }
    std::ostringstream os;
    os << "erf.ibseb: the fixed sun stands " << elev << " deg above the horizon, below erf.ibseb.sun_min_elev_deg = "
       << m_params.sun_min_elev_deg << ", and the ray map of level " << m_lev << " does not span the domain: a building"
       " farther than " << m_halo[0] * m_dx[0] << " m in x or " << m_halo[1] * m_dx[1] << " m in y from a face casts no shadow"
       " on it. Lower erf.ibseb.sun_min_elev_deg to the sun's elevation for exact shadows.";
    if (ParallelDescriptor::IOProcessor()) { Warning(os.str().c_str()); }
}

ibseb::RayMap
IBFaceSet::ray_map () const
{
    ibseb::RayMap m;
    if (m_patch_box.ok()) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_ray_patch.local_size() == 1, "erf.ibseb: the ray map holds no patch for this rank");
        m.top = m_ray_patch.const_array(m_ray_patch.IndexArray()[0]);
    }
    m.ilo = m_dom_lo[0]; m.jlo = m_dom_lo[1];
    m.nx = m_nx; m.ny = m_ny;
    m.wrap_x = m_wrap_x; m.wrap_y = m_wrap_y;
    m.x_lo = m_x_lo; m.y_lo = m_y_lo;
    m.dx = m_dx[0]; m.dy = m_dx[1]; m.dz = m_dx[2];
    m.z_ground = m_z_ground;
    m.z_max = m_col_top_max;
    m.max_steps = m_max_steps;
    return m;
}

/**
 * Rebuild the ray map with the buildings of the levels below (see the
 * header and build_ray_map()).
 */
void
IBFaceSet::add_outside_occluders (const std::vector<const IBFaceSet*>& below)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!below.empty(), "erf.ibseb: add_outside_occluders() needs the levels below");
    build_ray_map(below);
    warn_fixed_sun_below_floor();
}

/**
 * Map this level's buildings to level 0's numbering through the labels of
 * the level below (see the header). Every column of a building votes with
 * the level-0 id of the coarse column under it (through the level below's own
 * map), unlabelled coarse columns abstaining, and the most frequent id wins
 * (the lowest on a tie). No single column is trusted: the coarse blanking is
 * a volume average, so a column the coarse level leaves unlabelled (a spire
 * alone in its coarse cell) cannot unmap a building the level below
 * resolves through its other columns. A building none of whose columns lies
 * under a labelled coarse column is one only this level resolves; it maps to
 * 0, and a warning names it when the level below has buildings at all.
 * The level below's ids are brought under this level's columns by a
 * parallel copy; the votes, a few per building, are gathered from every rank
 * and summed, so every rank reaches the same verdict.
 */
void
IBFaceSet::map_buildings_to_level0 (const IBFaceSet& coarser, const IntVect& ratio)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(coarser.m_lab.ok(),
                                     "erf.ibseb: the level below released its labels before this level was mapped");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(coarser.m_geom.Domain().smallEnd(2) == m_geom.Domain().smallEnd(2),
                                     "erf.ibseb: the levels do not share the ground");
    const int klo = m_geom.Domain().smallEnd(2);
    const IntVect ratio_h(AMREX_D_DECL(ratio[0], ratio[1], 1));
    iMultiFab crse(amrex::coarsen(m_lab.boxArray(), ratio_h), m_lab.DistributionMap(), 1, 0);
    crse.setVal(0);
    crse.ParallelCopy(coarser.m_lab, 0, 0, 1, IntVect(0), IntVect(0), coarser.m_geom.periodicity());
    // This rank's votes: (this level's id, level 0's id) -> columns.
    std::map<std::pair<int, int>, Long> votes;
    for (MFIter mfi(m_lab); mfi.isValid(); ++mfi) {
        const Box& vb = mfi.validbox();
        const HostInts lab(m_lab[mfi]);
        const HostInts cl(crse[mfi]);
        auto const& a = lab.array();
        auto const& c = cl.array();
        const auto lo = lbound(vb);
        const auto hi = ubound(vb);
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                const int b = a(i, j, klo);
                if (b == 0) { continue; }
                // amrex::coarsen, not a division: integer division truncates
                // towards zero, which would send a negative index to the
                // coarse column above the one that covers it.
                const int bc = c(amrex::coarsen(i, ratio[0]), amrex::coarsen(j, ratio[1]), klo);
                if (bc > 0 && coarser.m_to_lev0[bc] > 0) { ++votes[{b, coarser.m_to_lev0[bc]}]; }
            }
        }
    }
    std::vector<Long> mine;
    for (const auto& [key, n] : votes) { mine.push_back(key.first); mine.push_back(key.second); mine.push_back(n); }
    const std::vector<Long> all = all_gather(mine);
    std::vector<std::map<int, Long>> tally(m_nbld + 1);
    for (size_t n = 0; n + 3 <= all.size(); n += 3) { tally[static_cast<int>(all[n])][static_cast<int>(all[n + 1])] += all[n + 2]; }
    std::vector<int> unmapped;
    for (int b = 1; b <= m_nbld; ++b) {
        int best = 0;
        Long n_best = 0;
        for (const auto& [id, n] : tally[b]) {     // ascending id: the lowest wins a tie
            if (n > n_best) { best = id; n_best = n; }
        }
        m_to_lev0[b] = best;
        if (best == 0) { unmapped.push_back(b); }
    }
    if (!unmapped.empty() && coarser.m_nbld > 0) {
        std::ostringstream os;
        os << "erf.ibseb: level " << m_lev << " has " << unmapped.size()
           << " building(s) that no labelled column of the level below lies under (level " << m_lev << " number";
        for (int b : unmapped) { os << " " << b; }
        os << "); they take erf.ibseb.material_default and building_level0 = 0 in the report";
        // Every rank holds the same votes and reaches the same verdict.
        if (ParallelDescriptor::IOProcessor()) { Warning(os.str().c_str()); }
    }
    if (m_params.debug) {
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " buildings as level 0 numbers them:";
        for (int b = 1; b <= m_nbld; ++b) { Print() << " " << b << "->" << m_to_lev0[b]; }
        Print() << "\n";
    }
}

/**
 * Hemisphere sampling for the view fractions. One kernel over the faces,
 * each looping over its rays; the count of rays ending on the sky, the
 * ground and a building over the total gives the three fractions. Roofs
 * point up and never see the ground; a wall on flat open ground sees half
 * sky and half ground, which the regression test checks.
 */
void
IBFaceSet::compute_view_fractions ()
{
    const int n_az = m_params.view_n_az, n_el = m_params.view_n_el;
    const ibseb::RayMap rm = ray_map();
    const int*  pd = d_dir.data();  const int* ps = d_side.data();
    const Real* pxf = d_xf.data();  const Real* pyf = d_yf.data();  const Real* pzf = d_zf.data();
    Real* pfs = d_f_sky.data(); Real* pfg = d_f_ground.data(); Real* pfb = d_f_bldg.data();
    ParallelFor(m_nface, [=] AMREX_GPU_DEVICE (int f) noexcept {
        int n_sky = 0, n_gnd = 0, n_bld = 0;
        for (int ie = 0; ie < n_el; ++ie) {
            for (int ia = 0; ia < n_az; ++ia) {
                Real sx, sy, sz;
                ibseb::hemisphere_direction(pd[f], -ps[f], ia, ie, n_az, n_el, sx, sy, sz);
                const int hit = ibseb::ray_hit(pxf[f], pyf[f], pzf[f], sx, sy, sz, rm);
                if (hit == ibseb::RAY_SKY) ++n_sky; else if (hit == ibseb::RAY_GROUND) ++n_gnd; else ++n_bld;
            }
        }
        const Real inv = 1.0 / static_cast<Real>(n_az * n_el);
        pfs[f] = n_sky * inv; pfg[f] = n_gnd * inv; pfb[f] = n_bld * inv;
    });
    Gpu::streamSynchronize();

    if (m_params.debug) {
        std::vector<Real> fs(m_nface), fg(m_nface), fb(m_nface);
        Gpu::copy(Gpu::deviceToHost, d_f_sky.begin(), d_f_sky.end(), fs.begin());
        Gpu::copy(Gpu::deviceToHost, d_f_ground.begin(), d_f_ground.end(), fg.begin());
        Gpu::copy(Gpu::deviceToHost, d_f_bldg.begin(), d_f_bldg.end(), fb.begin());
        Gpu::streamSynchronize();
        Real s[3] = {0.0, 0.0, 0.0};
        for (int n = 0; n < m_nface; ++n) { s[0] += fs[n]; s[1] += fg[n]; s[2] += fb[n]; }
        ParallelDescriptor::ReduceRealSum(s, 3);
        // A refined level can hold no building at all.
        const Real nf_all = amrex::max(Real(1.0), static_cast<Real>(m_nface_dir[0] + m_nface_dir[1] + m_nface_dir[2]));
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " view fractions: " << n_az * n_el
                << " rays per face, mean f_sky=" << s[0] / nf_all << " f_ground=" << s[1] / nf_all
                << " f_bldg=" << s[2] / nf_all << "\n";
    }
}

/**
 * Longwave of the current step. One kernel per fab over that fab's faces,
 * reading the air temperature of the fluid cell from the conserved state
 * (dry density, potential temperature and, with a moisture model, water
 * vapour through the equation of state; cell_qv()). With
 * lw_mode = two_stream the sky and ground terms come from the canopy forcing
 * of the face's column, read from the 2D fab of the same local index (the
 * level's grids flattened, on the same distribution).
 */
void
IBFaceSet::compute_longwave (const MultiFab& cons, const TwoStreamCanopyView& canopy)
{
    if (m_params.lw_mode == "two_stream") {
        compute_longwave_two_stream(cons, canopy);
        return;
    }
    const Real  sigma  = ibseb::SIGMA_SB;
    const bool  gray   = (m_params.lw_mode == "gray");
    const Real  lw_fix = m_params.lw_down, eps_sky = m_params.sky_emissivity;
    const Real  eps_g  = m_params.emissivity_ground;
    const Real* peps   = d_emis.data();
    const Real  Tg4    = std::pow(m_params.T_ground, 4);
    const int*  pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const Real* pfs = d_f_sky.data(); const Real* pfg = d_f_ground.data(); const Real* pfb = d_f_bldg.data();
    const Real* pT = d_T_skin.data();
    Real* pTa = d_T_air.data(); Real* pin = d_LW_down_in.data(); Real* pnet = d_LW_net.data();
    Real* pext = d_LW_ext.data();
    const int qv_comp = m_qv_comp;
    for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        auto const& c = cons.const_array(mfi);
        ParallelFor(f1 - f0, [=] AMREX_GPU_DEVICE (int m) noexcept {
            const int f = f0 + m;
            const Real rho = c(pi[f], pj[f], pk[f], Rho_comp);
            const Real qv  = cell_qv(c, pi[f], pj[f], pk[f], qv_comp);
            const Real Ta  = getTgivenRandRTh(rho, c(pi[f], pj[f], pk[f], RhoTheta_comp), qv);
            const Real lw_sky = gray ? eps_sky * sigma * Ta * Ta * Ta * Ta : lw_fix;
            const Real Ts4 = pT[f] * pT[f] * pT[f] * pT[f];
            // Sky and ground do not depend on the face; the wall term does
            // (isothermal surroundings), so the balance keeps them apart.
            const Real lw_ext = pfs[f] * lw_sky + pfg[f] * eps_g * sigma * Tg4;
            const Real lw_in  = lw_ext + pfb[f] * sigma * Ts4;
            pTa[f]  = Ta;
            pext[f] = lw_ext;
            pin[f]  = lw_in;
            pnet[f] = peps[f] * (lw_in - sigma * Ts4);
        });
    }
    Gpu::streamSynchronize();
}

/**
 * Longwave of the current step with lw_mode = two_stream: the sky term is the
 * longwave down of the face's column at the face's height (a roof's interface,
 * a wall's mean of its cell's two; two_stream_canopy_sample()), the ground term
 * its longwave up at the same height (the ground seen through the air in between),
 * so the ground's temperature, emissivity and reflection of the sky are the
 * column's (two_stream_face_sky()). Before the
 * level's first sweep (an empty ``canopy``, the initial report) the air
 * temperature is read and the incoming and net longwave are zero.
 */
void
IBFaceSet::compute_longwave_two_stream (const MultiFab& cons, const TwoStreamCanopyView& canopy)
{
    const Real  sigma = ibseb::SIGMA_SB;
    const bool  have  = canopy.valid();
    const int   k0    = m_dom_lo[2];
    const Real* peps  = d_emis.data();
    const int*  pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const int*  pd = d_dir.data();  const int* ps = d_side.data();
    const Real* pfs = d_f_sky.data(); const Real* pfg = d_f_ground.data(); const Real* pfb = d_f_bldg.data();
    const Real* pT = d_T_skin.data();
    Real* pTa = d_T_air.data(); Real* pin = d_LW_down_in.data(); Real* pnet = d_LW_net.data();
    Real* pext = d_LW_ext.data();
    const int qv_comp = m_qv_comp;
    for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        auto const& c = cons.const_array(mfi);
        const auto beam = have ? canopy.beam->const_array(mfi) : Array4<const Real>{};
        const auto cosz = have ? canopy.cos_zenith->const_array(mfi) : Array4<const Real>{};
        const auto flux = have ? canopy.fluxes->const_array(mfi) : Array4<const Real>{};
        ParallelFor(f1 - f0, [=] AMREX_GPU_DEVICE (int m) noexcept {
            const int f = f0 + m;
            const Real rho = c(pi[f], pj[f], pk[f], Rho_comp);
            pTa[f] = getTgivenRandRTh(rho, c(pi[f], pj[f], pk[f], RhoTheta_comp),
                                      cell_qv(c, pi[f], pj[f], pk[f], qv_comp));
            if (!have) {
                pext[f] = 0.0; pin[f] = 0.0; pnet[f] = 0.0;
                return;
            }
            const Real Ts4 = pT[f] * pT[f] * pT[f] * pT[f];
            const TwoStreamFaceSky sky = two_stream_face_sky(pi[f], pj[f], pk[f], k0, pd[f], ps[f], beam, cosz, flux);
            const Real lw_ext = pfs[f] * sky.lw_down + pfg[f] * sky.lw_up;
            const Real lw_in  = lw_ext + pfb[f] * sigma * Ts4;
            pext[f] = lw_ext;
            pin[f]  = lw_in;
            pnet[f] = peps[f] * (lw_in - sigma * Ts4);
        });
    }
    Gpu::streamSynchronize();
}

/**
 * The fields of a two-stream view must lie on this level's columns, with the
 * level's distribution, so that the fabs of local index n hold the columns of
 * the faces fab_start()[n] .. fab_start()[n+1]-1 own: the interface fluxes on
 * the level's own grids (with the top-of-atmosphere ghost, which a face whose
 * cell is the top of the domain reads); the beam on its boxes in x and y with
 * the interfaces from the surface up to at least top_sample_interface() in z;
 * the cosine of the zenith on them flattened to index 0 in z.
 */
void
IBFaceSet::check_canopy_layout (const MultiFab& cons, const TwoStreamCanopyView& canopy) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(canopy.valid(),
        "erf.ibseb.radiation = two_stream: level " + std::to_string(m_lev) + "'s two-stream view is incomplete");
    const BoxArray& ba = cons.boxArray();
    const DistributionMapping& dm = cons.DistributionMap();
    const auto& beam = *canopy.beam;
    const auto& cosz = *canopy.cos_zenith;
    const auto& flux = *canopy.fluxes;
    bool ok = flux.boxArray() == ba && flux.DistributionMap() == dm && flux.nComp() >= 4 && flux.nGrowVect()[2] >= 1 &&
              beam.boxArray().size() == ba.size() && beam.DistributionMap() == dm &&
              cosz.boxArray().size() == ba.size() && cosz.DistributionMap() == dm;
    // Same columns: the boxes agree in x and y.
    for (int n = 0; ok && n < static_cast<int>(ba.size()); ++n) {
        const Box& b = beam.boxArray()[n];
        const Box& c = cosz.boxArray()[n];
        for (int d = 0; d < 2; ++d) {
            ok = ok && b.smallEnd(d) == ba[n].smallEnd(d) && b.bigEnd(d) == ba[n].bigEnd(d) &&
                       c.smallEnd(d) == ba[n].smallEnd(d) && c.bigEnd(d) == ba[n].bigEnd(d);
        }
        // The kernels read the beam at interfaces 0 .. top_sample_interface() and the
        // cosine of the zenith at index 0 in z.
        ok = ok && b.smallEnd(2) == 0 && b.bigEnd(2) >= top_sample_interface() && c.smallEnd(2) == 0;
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ok,
        "erf.ibseb.radiation = two_stream: level " + std::to_string(m_lev) +
        "'s two-stream beam, cosine of the zenith or interface fluxes are not laid out on the level's grids");
}

/**
 * Materials per face from the library by building id, or the uniform
 * inputs. Host work on the static geometry, once.
 */
void
IBFaceSet::assign_materials ()
{
    if (m_params.material_file.empty()) { return; }
    const Vector<IBSEBMaterial> table = read_ibseb_materials(m_params.material_file);
    auto find = [&](int id) -> const IBSEBMaterial& {
        for (const auto& m : table) { if (m.mat_id == id) return m; }
        Abort("erf.ibseb: material id " + std::to_string(id) + " is not in " + m_params.material_file);
        return table[0];
    };
    std::vector<int> h_b(m_nface);
    Gpu::copy(Gpu::deviceToHost, d_bid.begin(), d_bid.end(), h_b.begin());
    Gpu::streamSynchronize();
    std::vector<int>  h_mat(m_nface);
    std::vector<Real> h_alb(m_nface), h_emi(m_nface), h_k(m_nface), h_rc(m_nface), h_th(m_nface);
    for (int n = 0; n < m_nface; ++n) {
        // erf.ibseb.material_by_building follows level 0's numbering; a
        // building only a refined level resolves (0 here) takes the default.
        const int b = m_to_lev0[h_b[n]];
        int id = m_params.material_default;
        if (b >= 1 && b <= static_cast<int>(m_params.material_by_building.size())) {
            id = m_params.material_by_building[b - 1];
        }
        const IBSEBMaterial& m = find(id);
        h_mat[n] = m.mat_id; h_alb[n] = m.albedo; h_emi[n] = m.emissivity;
        h_k[n] = m.k_therm; h_rc[n] = m.rho_cp; h_th[n] = m.thickness;
    }
    upload(d_mat, h_mat); upload(d_albedo, h_alb); upload(d_emis, h_emi);
    upload(d_kth, h_k); upload(d_rhocp, h_rc); upload(d_thick, h_th);
    Gpu::streamSynchronize();
    if (m_params.debug) {
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " materials: " << table.size() << " in "
                << m_params.material_file << ", default id " << m_params.material_default;
        for (size_t b = 0; b < m_params.material_by_building.size(); ++b) {
            Print() << ", building " << b + 1 << " -> " << m_params.material_by_building[b];
        }
        Print() << "\n";
        for (const auto& m : table) {
            Print() << "[IBSEB DEBUG]   material " << m.mat_id << " (" << m.name << "): albedo=" << m.albedo
                    << " emissivity=" << m.emissivity << " k=" << m.k_therm << " W/m/K rho_cp=" << m.rho_cp
                    << " J/m3/K thickness=" << m.thickness << " m\n";
        }
    }
}

/**
 * Slab conduction on every face: one kernel over the faces, each solving
 * its own tridiagonal system with the skin at the top and the interior at
 * the bottom.
 */
void
IBFaceSet::compute_ground (Real dt)
{
    const int  nl   = n_layers();
    const Real Tint = m_params.T_interior;
    const Real* pT  = d_T_skin.data();
    const Real* pkt = d_kth.data(); const Real* prc = d_rhocp.data(); const Real* pth = d_thick.data();
    Real* pS = d_T_slab.data();   Real* pG = d_G.data();
    ParallelFor(m_nface, [=] AMREX_GPU_DEVICE (int f) noexcept {
        // Zero-initialized: only the first nl entries are ever used, but the
        // compiler cannot see that through the inlined solve and warns otherwise.
        Real T[ibseb::SLAB_MAX_LAYERS] = {};
        for (int l = 0; l < nl; ++l) { T[l] = pS[f * nl + l]; }
        const Real dz = pth[f] / nl;
        pG[f] = ibseb::advance_slab_dirichlet(T, pT[f], Tint, pkt[f], prc[f], dz, dt, nl);
        for (int l = 0; l < nl; ++l) { pS[f * nl + l] = T[l]; }
    });
    Gpu::streamSynchronize();
}

/**
 * Prognostic balance on every face: one kernel over the faces, each with its
 * own slab response, Newton solve and slab advance. Everything the solve
 * needs was stored by the three flux routines of the step, so no fab data is
 * touched here. The fluxes are rewritten at the solution so that the dump,
 * the reports and the heat deposited into the air all describe the same
 * closed balance; d_resid holds what is left of it.
 */
void
IBFaceSet::solve_balance (Real dt)
{
    const int  nl   = n_layers();
    const Real Tint = m_params.T_interior;
    const Real Tmin = m_params.T_skin_min, Tmax = m_params.T_skin_max;
    const Real tol  = m_params.newton_tol_K, max_step = m_params.newton_max_step_K;
    const int  max_iter = m_params.newton_max_iter;
    const Real* pkt = d_kth.data(); const Real* prc = d_rhocp.data(); const Real* pth = d_thick.data();
    const Real* peps = d_emis.data(); const Real* pfb = d_f_bldg.data();
    const Real* psw = d_SW_abs.data(); const Real* pq = d_Q_ext.data(); const Real* pext = d_LW_ext.data();
    const Real* pHc = d_H_coeff.data(); const Real* pTa = d_T_air.data(); const Real* ptha = d_theta_air.data();
    const Real* pLE = d_LE.data();
    Real* pT = d_T_skin.data(); Real* pS = d_T_slab.data(); Real* pG = d_G.data();
    Real* pH = d_H.data(); Real* pin = d_LW_down_in.data(); Real* pnet = d_LW_net.data();
    Real* pres = d_resid.data(); int* pnit = d_niter.data();
    ParallelFor(m_nface, [=] AMREX_GPU_DEVICE (int f) noexcept {
        // Zero-initialized: only the first nl entries are ever used, but the
        // compiler cannot see that through the inlined solve and warns otherwise.
        Real T[ibseb::SLAB_MAX_LAYERS] = {};
        for (int l = 0; l < nl; ++l) { T[l] = pS[f * nl + l]; }
        const Real dz = pth[f] / nl;
        // Conduction the implicit slab step will take for a skin temperature Ts: a Ts - b.
        Real a = 0.0, b = 0.0;
        ibseb::slab_skin_response(T, Tint, pkt[f], prc[f], dz, dt, nl, a, b);
        const Real exner = pTa[f] / ptha[f];
        const Real q_abs = peps[f] * pq[f];
        int n_it = 0; Real res = 0.0;
        const Real Ts = ibseb::solve_skin_balance(pT[f], psw[f], q_abs, pext[f], peps[f], pfb[f],
                                                  pHc[f], exner, ptha[f], pLE[f], a, b,
                                                  Tmin, Tmax, max_step, tol, max_iter, n_it, res);
        // Advance the slab with the solution and store every term at it.
        pG[f] = ibseb::advance_slab_dirichlet(T, Ts, Tint, pkt[f], prc[f], dz, dt, nl);
        for (int l = 0; l < nl; ++l) { pS[f * nl + l] = T[l]; }
        const Real Ts4 = Ts * Ts * Ts * Ts;
        pin[f]  = pext[f] + pfb[f] * ibseb::SIGMA_SB * Ts4;
        pnet[f] = peps[f] * (pin[f] - ibseb::SIGMA_SB * Ts4);
        pH[f]   = pHc[f] * (Ts / exner - ptha[f]);
        pT[f]   = Ts;
        // The residual with the stored terms (the slab's own G, not a Ts - b),
        // so it reports the closure of what is written out.
        pres[f] = std::abs(psw[f] + q_abs + pnet[f] - pH[f] - pLE[f] - pG[f]);
        pnit[f] = n_it;
    });
    Gpu::streamSynchronize();

    if (m_params.debug) {
        std::vector<Real> h_r(m_nface), h_T(m_nface);
        std::vector<int>  h_n(m_nface);
        Gpu::copy(Gpu::deviceToHost, d_resid.begin(), d_resid.end(), h_r.begin());
        Gpu::copy(Gpu::deviceToHost, d_T_skin.begin(), d_T_skin.end(), h_T.begin());
        Gpu::copy(Gpu::deviceToHost, d_niter.begin(), d_niter.end(), h_n.begin());
        Gpu::streamSynchronize();
        Real rmax = 0.0; int nmax = 0; Long nclamp = 0;
        for (int n = 0; n < m_nface; ++n) {
            rmax = std::max(rmax, h_r[n]); nmax = std::max(nmax, h_n[n]);
            if (h_T[n] <= Tmin || h_T[n] >= Tmax) { ++nclamp; }
        }
        ParallelDescriptor::ReduceRealMax(rmax);
        ParallelDescriptor::ReduceIntMax(nmax);
        ParallelDescriptor::ReduceLongSum(nclamp);
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " balance: dt=" << dt << " s newton iterations max=" << nmax
                << " residual max=" << rmax << " W/m2 faces at a bound=" << nclamp << "\n";
    }
}

/**
 * Wall function on every face. One kernel per fab over that fab's faces:
 * the cell-centred velocity of the fluid cell from its face values, its
 * tangential part with respect to the wall, the neutral log law for u*
 * and theta*, and on request two additions: a convective
 * velocity scale in the wind (Beljaars' gustiness form with Deardorff's
 * w* from the previous step's flux) and, on roofs, the surface layer's
 * similarity functions iterated on the face's own Obukhov length. The
 * skin temperature is converted to a potential temperature with the fluid
 * cell's Exner function before the difference is taken.
 *
 * The Obukhov length of a face is ``L = u*^2 theta / (kappa g theta*)``
 * (negative over a hot face). Walls stay on the log law even with the
 * stability correction: the similarity functions assume a horizontal
 * surface, and the convective scale is what carries free convection on a
 * wall. The ground's 2D Obukhov field only seeds the roof iteration.
 */
void
IBFaceSet::compute_sensible (const MultiFab& cons, const MultiFab& xvel,
                             const MultiFab& yvel, const MultiFab& zvel, Real c_p,
                             const MultiFab* olen_ground, const MultiFab* pblh_ground, Real z_i_bulk)
{
    const Real z0   = m_params.z0_wall, z0h = m_params.z0h_wall;
    const bool stab = m_params.stability_correction;
    const bool conv = (m_params.convective_velocity == "deardorff");
    const bool seed = stab && (olen_ground != nullptr) && (m_params.obukhov_seed == "ground");
    const bool use_pblh = conv && (pblh_ground != nullptr) && (m_params.z_i_mode == "pblh");
    const Real beta = m_params.beta_conv;
    const Real relax = m_params.obukhov_relax;
    const Real zi_fallback = (z_i_bulk > 0.0) ? z_i_bulk : m_params.z_i;
    const Real dxa[3] = {m_dx[0], m_dx[1], m_dx[2]};
    const Real z_lo = m_z_ground;
    const int*  pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const int*  pd = d_dir.data();
    const Real* pT = d_T_skin.data(); const Real* phb = d_hbld.data(); const Real* pzf = d_zf.data();
    Real* pth = d_theta_air.data(); Real* prho = d_rho.data(); Real* pU = d_U_tan.data();
    Real* pus = d_ustar.data();     Real* pH = d_H.data();     Real* pTa = d_T_air.data();
    Real* pHc = d_H_coeff.data();   Real* pws = d_w_star.data(); Real* pol = d_olen.data(); Real* pzi = d_z_i.data();
    const int qv_comp = m_qv_comp;
    for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        auto const& c = cons.const_array(mfi);
        auto const& u = xvel.const_array(mfi);
        auto const& v = yvel.const_array(mfi);
        auto const& w = zvel.const_array(mfi);
        // The surface layer's 2D fields share the level's distribution with
        // the boxes collapsed to k = 0, so the same iterator indexes them.
        Array4<const Real> ol = seed     ? olen_ground->const_array(mfi) : Array4<const Real>{};
        Array4<const Real> pb = use_pblh ? pblh_ground->const_array(mfi) : Array4<const Real>{};
        ParallelFor(f1 - f0, [=] AMREX_GPU_DEVICE (int m) noexcept {
            const int f = f0 + m;
            const int i = pi[f], j = pj[f], k = pk[f], d = pd[f];
            const Real rho = c(i, j, k, Rho_comp);
            const Real rth = c(i, j, k, RhoTheta_comp);
            const Real th  = rth / rho;
            const Real Ta  = getTgivenRandRTh(rho, rth, cell_qv(c, i, j, k, qv_comp));
            const Real th_skin = pT[f] * th / Ta;            // skin temperature as a potential temperature
            // Cell-centred velocity and its tangential part with respect to the wall.
            Real U[3] = { Real(0.5) * (u(i, j, k) + u(i + 1, j, k)),
                          Real(0.5) * (v(i, j, k) + v(i, j + 1, k)),
                          Real(0.5) * (w(i, j, k) + w(i, j, k + 1)) };
            U[d] = 0.0;
            const Real Ut = amrex::max(std::sqrt(U[0] * U[0] + U[1] * U[1] + U[2] * U[2]), Real(1.0e-3));
            // Convective velocity scale from the previous step's flux out of the face.
            Real wstar = 0.0, depth = 0.0;
            if (conv && pH[f] > 0.0) {
                if (d == 2) {
                    Real zi = use_pblh ? pb(i, j, 0) : zi_fallback;
                    if (!(zi > 0.0 && zi < 1.0e5)) { zi = zi_fallback; }
                    depth = amrex::max(zi - (pzf[f] - z_lo), phb[f]);
                } else {
                    depth = phb[f];
                }
                wstar = std::cbrt(CONST_GRAV / th * pH[f] / (rho * c_p) * depth);
            }
            const Real Ut_eff = std::sqrt(Ut * Ut + beta * beta * wstar * wstar);
            const Real delta = 0.5 * dxa[d];
            const Real lnm = std::log(delta / z0), lnh = std::log(delta / z0h);
            Real ustar = KAPPA * Ut_eff / lnm;
            Real thstar = KAPPA * (th - th_skin) / lnh;
            Real lnh_eff = lnh;
            Real olen = 1.0e30;
            if (stab && d == 2) {
                // Roofs: the surface layer's similarity functions on the
                // face's own Obukhov length, a few fixed-point passes, seeded
                // from the ground's field at this column when available.
                similarity_funs sf;
                Real L = 1.0e30;
                if (seed) {
                    const Real Lg = ol(i, j, 0);
                    if (std::abs(Lg) > 0.0 && std::abs(Lg) < 1.0e5) { L = Lg; }
                }
                if (!(std::abs(L) < 1.0e5)) {
                    L = ustar * ustar * th / (KAPPA * CONST_GRAV * thstar + 1.0e-20);
                }
                // Fixed point on L until u*, theta* and L agree to 1e-6, so
                // the stored values are mutually consistent (a fixed number
                // of passes leaves L one pass behind). u* is under-relaxed
                // between passes, as the surface layer's iteration is in
                // erf-model #3486, which keeps a strongly unstable roof from
                // oscillating at low wind; the converged values do not
                // depend on the factor.
                Real us_it = ustar;
                for (int it = 0; it < 50; ++it) {
                    const Real zeta = delta / L;
                    const Real psi_m = sf.calc_psi_m(zeta), psi_h = sf.calc_psi_h(zeta);
                    const Real us_new = KAPPA * Ut_eff / amrex::max(lnm - psi_m, Real(0.1));
                    us_it   = (it == 0) ? us_new : (1.0 - relax) * us_it + relax * us_new;
                    lnh_eff = amrex::max(lnh - psi_h, Real(0.1));
                    thstar  = KAPPA * (th - th_skin) / lnh_eff;
                    const Real L_new = us_it * us_it * th / (KAPPA * CONST_GRAV * thstar + 1.0e-20);
                    const bool done = std::abs(L_new - L) <= 1.0e-6 * std::abs(L_new) &&
                                      std::abs(us_new - us_it) <= 1.0e-6 * us_new;
                    L = L_new;
                    if (done) { us_it = us_new; break; }
                }
                ustar = us_it;
                olen  = L;
            }
            pth[f]  = th;
            pTa[f]  = Ta;
            prho[f] = rho;
            pU[f]   = Ut;
            pus[f]  = ustar;
            pws[f]  = wstar;
            pol[f]  = olen;
            pzi[f]  = depth;
            pHc[f]  = rho * c_p * KAPPA * ustar / lnh_eff;   // H = H_coeff (theta_skin - theta_air)
            pH[f]   = -rho * c_p * ustar * thstar;           // positive out of the face
        });
    }
    Gpu::streamSynchronize();
}

/**
 * Face flux into the temperature equation. The rho-theta tendency of the
 * fluid cell is H A / (c_p V Pi); with several faces on one cell the adds
 * are atomic.
 */
void
IBFaceSet::add_heat_flux_to_source (MultiFab& source, const MultiFab& cons,
                                    const Geometry& geom, Real c_p, Real rdOcp) const
{
    const auto dx = geom.CellSizeArray();
    const Real vol = dx[0] * dx[1] * dx[2];
    const int*  pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const Real* pH = d_H.data();  const Real* pA = d_area.data();
    const int qv_comp = m_qv_comp;
    for (MFIter mfi(source); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        auto const& src = source.array(mfi);
        auto const& c   = cons.const_array(mfi);
        ParallelFor(f1 - f0, [=] AMREX_GPU_DEVICE (int m) noexcept {
            const int f = f0 + m;
            const int i = pi[f], j = pj[f], k = pk[f];
            const Real exner = getExnergivenRTh(c(i, j, k, RhoTheta_comp), rdOcp, cell_qv(c, i, j, k, qv_comp));
            const Real tend  = pH[f] * pA[f] / (c_p * vol * exner);
            Gpu::Atomic::AddNoRet(&src(i, j, k, RhoTheta_comp), tend);
        });
    }
    Gpu::streamSynchronize();
}

/**
 * Shortwave of the current step. The sun and the irradiances are scalars
 * computed on the host from the provider inputs; the per-face work is one
 * kernel: incidence cosine, ray cast against the column tops, direct and
 * diffuse incident, absorbed.
 *
 * With sun_mode = fixed the direct-normal irradiance and the horizontal
 * diffuse are the inputs as given. With sun_mode = solar the sun follows the
 * site and time and the clear-sky formulas give both irradiances; with
 * sun_mode = two_stream the position and the top-of-atmosphere irradiance are
 * the two-stream radiation's (set_two_stream_sun()) and the same clear-sky
 * formulas give the irradiances. The sun below the horizon gives zero
 * everywhere.
 *
 * Diffuse on a face: ``f_sky * diffuse_h + f_ground * albedo_ground *
 * (direct_h + diffuse_h)``, the second term being the ground-reflected part,
 * with the direct on a horizontal surface ``direct_h = dni * cos z``. The
 * view fractions are the placeholders of build() until
 * compute_view_fractions() has run.
 */
void
IBFaceSet::compute_shortwave (Real time, const TwoStreamCanopyView& canopy)
{
    // ---- Sun and irradiances (host scalars) ----
    SunState s;
    if (m_params.sun_mode == "fixed") {
        s.zenith  = m_params.sun_zenith_deg  * PI / 180.0;
        s.azimuth = m_params.sun_azimuth_deg * PI / 180.0;
        const Real cz = std::cos(s.zenith);
        s.dni       = (cz > 0.0) ? m_params.sw_direct_normal : 0.0;
        s.diffuse_h = (cz > 0.0) ? m_params.sw_diffuse : 0.0;
    } else if (m_params.sun_mode == "two_stream") {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_ts_sun.set,
            "erf.ibseb.sun_mode = two_stream: the two-stream sun was not given before the shortwave");
        // Used once: the next shortwave needs the sun of its own time.
        m_ts_sun.set = false;
        const Real decl = m_ts_sun.decl;
        const Real ha   = m_ts_sun.hour_angle;
        s.zenith  = ibseb::solar_zenith(m_ts_sun.lat_deg, decl, ha);
        s.azimuth = ibseb::solar_azimuth(m_ts_sun.lat_deg, decl, ha, s.zenith);
        const Real cz = std::cos(s.zenith);
        // The irradiance carries the date's distance factor already.
        s.dni       = ibseb::clear_sky_dni(cz, m_ts_sun.S0, m_params.sw_transmission, Real(1.0));
        s.diffuse_h = ibseb::clear_sky_diffuse_h(cz, m_ts_sun.S0, m_params.sw_transmission,
                                                 Real(1.0), m_params.sw_diffuse_coeff);
    } else {
        const Real t_utc = m_params.time_zero_utc_s + time;
        // The date advances at each 00:00 UTC, so a run of several days follows the season.
        const Real day   = ibseb::day_of_year_at(m_params.day_of_year, t_utc);
        const Real decl  = ibseb::solar_declination(day);
        // Hour angle from UTC and the longitude (a time-zone offset would cancel).
        const Real ha    = ibseb::solar_hour_angle(t_utc, m_params.longitude_deg, day);
        s.zenith  = ibseb::solar_zenith(m_params.latitude_deg, decl, ha);
        s.azimuth = ibseb::solar_azimuth(m_params.latitude_deg, decl, ha, s.zenith);
        const Real cz = std::cos(s.zenith);
        const Real e0 = ibseb::earth_sun_distance_factor(day);
        s.dni       = ibseb::clear_sky_dni(cz, m_params.solar_constant, m_params.sw_transmission, e0);
        s.diffuse_h = ibseb::clear_sky_diffuse_h(cz, m_params.solar_constant, m_params.sw_transmission,
                                                 e0, m_params.sw_diffuse_coeff);
    }
    ibseb::sun_vector(s.zenith, s.azimuth, s.sx, s.sy, s.sz);
    m_sun = s;

    if (m_params.radiation == "two_stream") {
        compute_shortwave_two_stream(canopy);
        return;
    }

    // ---- Per face ----
    const Real sx = s.sx, sy = s.sy, sz = s.sz;
    const Real dni = s.dni, dif_h = s.diffuse_h;
    const Real dir_h = dni * amrex::max(Real(0.0), sz);
    const Real alb_g = m_params.albedo_ground;
    const Real* palb = d_albedo.data();
    const ibseb::RayMap rm = ray_map();
    const int*  pd  = d_dir.data();  const int* ps = d_side.data();
    const Real* pxf = d_xf.data();   const Real* pyf = d_yf.data();  const Real* pzf = d_zf.data();
    const Real* pfs = d_f_sky.data(); const Real* pfg = d_f_ground.data();
    Real* psh = d_shadow.data();  Real* pdir = d_SW_direct_in.data();
    Real* pdif = d_SW_diffuse_in.data();  Real* pabs = d_SW_abs.data();
    ParallelFor(m_nface, [=] AMREX_GPU_DEVICE (int f) noexcept {
        // Outward normal: opposite to the side the solid is on.
        Real n[3] = {0.0, 0.0, 0.0};
        n[pd[f]] = -static_cast<Real>(ps[f]);
        const Real cosi = n[0] * sx + n[1] * sy + n[2] * sz;
        Real shadow = 0.0, direct = 0.0;
        if (sz > 0.0 && cosi > 0.0 && dni > 0.0) {
            shadow = ibseb::ray_blocked(pxf[f], pyf[f], pzf[f], sx, sy, sz, rm) ? 1.0 : 0.0;
            direct = dni * cosi * (1.0 - shadow);
        }
        const Real diffuse = pfs[f] * dif_h + pfg[f] * alb_g * (dir_h + dif_h);
        psh[f]  = shadow;
        pdir[f] = direct;
        pdif[f] = diffuse;
        pabs[f] = (1.0 - palb[f]) * (direct + diffuse);
    });
    Gpu::streamSynchronize();

    if (m_params.debug) {
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " sun: zenith=" << s.zenith * 180.0 / PI
                << " deg azimuth=" << s.azimuth * 180.0 / PI << " deg s=(" << s.sx << "," << s.sy << "," << s.sz
                << ") DNI=" << s.dni << " W/m2 diffuse_h=" << s.diffuse_h << " W/m2\n";
    }
}

/**
 * Shortwave of the current step with radiation = two_stream, after
 * compute_shortwave() has placed the sun (m_sun): per face, the direct-normal
 * irradiance of its column at its height (two_stream_face_sky(): a roof's
 * interface, a wall's mean of its cell's two), the beam over the column's own
 * cosine of the zenith, on the face (shadowed by the same ray cast), the
 * diffuse light of the sky there through f_sky, and the ground's reflection,
 * the column's shortwave up at the face's height, through f_ground. One kernel per
 * fab, reading the view's fields of the same local index.
 *
 * The faces' sun (sun_mode) and the sweep's are required to be one sun
 * (ERF::ibseb_check_sun_matches_two_stream()); the kernel checks it, column by
 * column: a cosine of the zenith more than 1e-4 (1e-3 in single precision)
 * from the faces' stops the run,
 * since a direct beam divided by the wrong cosine and laid along the wrong
 * direction would still look plausible.
 *
 * Before the level's first sweep (an empty ``canopy``, the initial report) every
 * face gets no shortwave. The reported direct-normal and diffuse irradiances
 * are the means over the level's faces of their columns'.
 */
void
IBFaceSet::compute_shortwave_two_stream (const TwoStreamCanopyView& canopy)
{
    const Real sx = m_sun.sx, sy = m_sun.sy, sz = m_sun.sz;
    const Real* palb = d_albedo.data();
    const ibseb::RayMap rm = ray_map();
    const int*  pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const int   k0 = m_dom_lo[2];
    const int*  pd  = d_dir.data();  const int* ps = d_side.data();
    const Real* pxf = d_xf.data();   const Real* pyf = d_yf.data();  const Real* pzf = d_zf.data();
    const Real* pfs = d_f_sky.data(); const Real* pfg = d_f_ground.data();
    Real* psh = d_shadow.data();  Real* pdir = d_SW_direct_in.data();
    Real* pdif = d_SW_diffuse_in.data();  Real* pabs = d_SW_abs.data();

    if (!canopy.valid()) {
        ParallelFor(m_nface, [=] AMREX_GPU_DEVICE (int f) noexcept {
            psh[f] = 0.0; pdir[f] = 0.0; pdif[f] = 0.0; pabs[f] = 0.0;
        });
        Gpu::streamSynchronize();
        m_sun.dni = 0.0;
        m_sun.diffuse_h = 0.0;
        return;
    }

    // Sums of the columns' direct-normal and diffuse irradiance over the faces,
    // for the report, and the largest disagreement of the two suns.
    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpMax> reduce_op;
    ReduceData<Real, Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    const Real sz_up = amrex::max(Real(0.0), sz);
    const int nfab = static_cast<int>(m_fab_start.size()) - 1;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(canopy.beam->local_size() == nfab,
        "erf.ibseb.radiation = two_stream: the two-stream beam holds another number of boxes than the faces' level");
    for (MFIter mfi(*canopy.beam); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        const auto beam = canopy.beam->const_array(mfi);
        const auto cosz = canopy.cos_zenith->const_array(mfi);
        const auto flux = canopy.fluxes->const_array(mfi);
        reduce_op.eval(f1 - f0, reduce_data, [=] AMREX_GPU_DEVICE (int m) noexcept -> ReduceTuple {
            const int f = f0 + m;
            const Real cz = cosz(pi[f], pj[f], 0);
            const TwoStreamFaceSky sky = two_stream_face_sky(pi[f], pj[f], pk[f], k0, pd[f], ps[f], beam, cosz, flux);
            const Real dni = sky.dni, dif_h = sky.diffuse, refl = sky.sw_up;
            // Outward normal: opposite to the side the solid is on.
            Real n[3] = {0.0, 0.0, 0.0};
            n[pd[f]] = -static_cast<Real>(ps[f]);
            const Real cosi = n[0] * sx + n[1] * sy + n[2] * sz;
            Real shadow = 0.0, direct = 0.0;
            if (sz > 0.0 && cosi > 0.0 && dni > 0.0) {
                shadow = ibseb::ray_blocked(pxf[f], pyf[f], pzf[f], sx, sy, sz, rm) ? 1.0 : 0.0;
                direct = dni * cosi * (1.0 - shadow);
            }
            const Real diffuse = pfs[f] * dif_h + pfg[f] * refl;
            psh[f]  = shadow;
            pdir[f] = direct;
            pdif[f] = diffuse;
            pabs[f] = (1.0 - palb[f]) * (direct + diffuse);
            return {dni, dif_h, std::abs(cz - sz_up)};
        });
    }
    Gpu::streamSynchronize();
    const auto rv = reduce_data.value(reduce_op);
    Real sums[2] = {amrex::get<0>(rv), amrex::get<1>(rv)};
    Real sun_gap = amrex::get<2>(rv);
    ParallelDescriptor::ReduceRealSum(sums, 2);
    ParallelDescriptor::ReduceRealMax(sun_gap);
    const Real nf_all = static_cast<Real>(m_nface_dir[0] + m_nface_dir[1] + m_nface_dir[2]);
    m_sun.dni       = (nf_all > 0.0) ? sums[0] / nf_all : Real(0.0);
    m_sun.diffuse_h = (nf_all > 0.0) ? sums[1] / nf_all : Real(0.0);
    // The two suns are one formula evaluated twice; a float calendar day resolves
    // about a second of the hour angle, hence the looser single-precision bound.
    const Real sun_tol = (sizeof(Real) == 8) ? Real(1.e-4) : Real(1.e-3);
    if (sun_gap > sun_tol) {
        Abort("erf.ibseb.radiation = two_stream: on level " + std::to_string(m_lev) + " the faces' sun (zenith "
              + std::to_string(m_sun.zenith * 180.0 / PI) + " deg) and the two-stream sweep's differ by up to "
              + std::to_string(sun_gap) + " in the cosine of the zenith angle (more than " + std::to_string(sun_tol)
              + "); the faces must take the sweep's sun"
                " (erf.ibseb.sun_mode = two_stream, or fixed at erf.fixed_solar_zenith_angle)");
    }
    if (m_params.debug) {
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " sun: zenith=" << m_sun.zenith * 180.0 / PI
                << " deg azimuth=" << m_sun.azimuth * 180.0 / PI << " deg s=(" << m_sun.sx << "," << m_sun.sy << ","
                << m_sun.sz << "); two-stream columns: mean DNI=" << m_sun.dni
                << " W/m2 mean diffuse_h=" << m_sun.diffuse_h << " W/m2 over the faces, largest sun gap "
                << sun_gap << " in cos z\n";
    }
}

/**
 * Scatter one per-face array as a per-cell mean, the same accumulation as
 * scatter_diagnostics() with the face count as the divisor.
 */
void
IBFaceSet::scatter_field (const Gpu::DeviceVector<Real>& v, MultiFab& out) const
{
    MultiFab cnt(out.boxArray(), out.DistributionMap(), 1, 0);
    cnt.setVal(0.0);
    out.setVal(0.0);
    const int* pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const Real* pv = v.data();
    for (MFIter mfi(out); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        auto const& n = cnt.array(mfi);
        auto const& o = out.array(mfi);
        ParallelFor(f1 - f0, [=] AMREX_GPU_DEVICE (int m) noexcept {
            const int f = f0 + m;
            Gpu::Atomic::AddNoRet(&n(pi[f], pj[f], pk[f]), Real(1.0));
            Gpu::Atomic::AddNoRet(&o(pi[f], pj[f], pk[f]), pv[f]);
        });
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (n(i, j, k) > 0.0) { o(i, j, k) /= n(i, j, k); }
        });
    }
    Gpu::streamSynchronize();
}

/**
 * Per-rank face dump. Every rank writes its own file so no gather is
 * needed; readers concatenate ``<prefix>.rank*.csv``.
 */
void
IBFaceSet::dump_faces (const std::string& prefix) const
{
    const int nf = m_nface;
    auto host = [&](const Gpu::DeviceVector<Real>& d) {
        std::vector<Real> h(nf);
        Gpu::copy(Gpu::deviceToHost, d.begin(), d.end(), h.begin());
        return h;
    };
    auto hosti = [&](const Gpu::DeviceVector<int>& d) {
        std::vector<int> h(nf);
        Gpu::copy(Gpu::deviceToHost, d.begin(), d.end(), h.begin());
        return h;
    };
    const auto i = hosti(d_i), j = hosti(d_j), k = hosti(d_k), dir = hosti(d_dir), side = hosti(d_side), bid = hosti(d_bid);
    const auto xf = host(d_xf), yf = host(d_yf), zf = host(d_zf), area = host(d_area);
    const auto fs = host(d_f_sky), fg = host(d_f_ground), fb = host(d_f_bldg);
    const auto sh = host(d_shadow), sdir = host(d_SW_direct_in), sdif = host(d_SW_diffuse_in), sabs = host(d_SW_abs);
    const auto T = host(d_T_skin), Ta = host(d_T_air), lwin = host(d_LW_down_in), lwn = host(d_LW_net);
    const auto th = host(d_theta_air), rho = host(d_rho), Ut = host(d_U_tan), us = host(d_ustar), H = host(d_H);
    const auto mat = hosti(d_mat);
    const auto alb = host(d_albedo), emi = host(d_emis), kk = host(d_kth), rc = host(d_rhocp), thk = host(d_thick), G = host(d_G);
    const auto Hc = host(d_H_coeff), Qx = host(d_Q_ext), LE = host(d_LE), lwe = host(d_LW_ext), res = host(d_resid);
    const auto nit = hosti(d_niter);
    const auto ws = host(d_w_star), ol = host(d_olen), zi = host(d_z_i), hb = host(d_hbld);
    const int nl = n_layers();
    std::vector<Real> slab(static_cast<size_t>(nf) * nl);
    Gpu::copy(Gpu::deviceToHost, d_T_slab.begin(), d_T_slab.end(), slab.begin());
    Gpu::streamSynchronize();
    std::ofstream f(prefix + ".rank" + std::to_string(ParallelDescriptor::MyProc()) + ".csv");
    f << "i,j,k,dir,side,bid,x_m,y_m,z_m,area_m2,f_sky,f_ground,f_bldg,shadow,SW_direct_in,SW_diffuse_in,SW_abs,T_skin,T_air,LW_in,LW_net,"
         "theta_air,rho,U_tan,ustar,H,mat,albedo,emissivity,k_therm,rho_cp,thickness,G,T_slab_top,T_slab_bottom,"
         "H_coeff,Q_ext,LE,LW_ext,resid,n_iter,w_star,Olen,z_i,h_bld";
    for (int l = 0; l < nl; ++l) { f << ",T_slab_" << l; }
    f << "\n";
    f << std::setprecision(12);
    for (int n = 0; n < nf; ++n) {
        f << i[n] << "," << j[n] << "," << k[n] << "," << dir[n] << "," << side[n] << "," << bid[n] << ","
          << xf[n] << "," << yf[n] << "," << zf[n] << "," << area[n] << ","
          << fs[n] << "," << fg[n] << "," << fb[n] << "," << sh[n] << ","
          << sdir[n] << "," << sdif[n] << "," << sabs[n] << "," << T[n] << ","
          << Ta[n] << "," << lwin[n] << "," << lwn[n] << ","
          << th[n] << "," << rho[n] << "," << Ut[n] << "," << us[n] << "," << H[n] << ","
          << mat[n] << "," << alb[n] << "," << emi[n] << "," << kk[n] << "," << rc[n] << "," << thk[n] << ","
          << G[n] << "," << slab[static_cast<size_t>(n) * nl] << "," << slab[static_cast<size_t>(n) * nl + nl - 1] << ","
          << Hc[n] << "," << Qx[n] << "," << LE[n] << "," << lwe[n] << "," << res[n] << "," << nit[n]
          << "," << ws[n] << "," << ol[n] << "," << zi[n] << "," << hb[n];
        for (int l = 0; l < nl; ++l) { f << "," << slab[static_cast<size_t>(n) * nl + l]; }
        f << "\n";
    }
}

/**
 * Debug description of the set, in the style of the fire module's
 * ``[FIRE DEBUG]`` lines: what was read, what every rank holds, what the
 * buildings look like, and what the arrays cost. Per-rank lines use
 * AllPrint() so every rank reports; the rest is printed by the I/O rank.
 */
void
IBFaceSet::print_debug_summary () const
{
    const int nl = n_layers();
    Print() << "[IBSEB DEBUG] lev=" << m_lev << " inputs: n_slab_layers=" << nl
            << " T_skin_init=" << m_params.T_skin_init << " K"
            << " T_interior=" << m_params.T_interior << " K"
            << " csv_file=" << m_params.csv_file << " csv_int=" << m_params.csv_int << "\n";
    Print() << "[IBSEB DEBUG] lev=" << m_lev << " global faces=" << (m_nface_dir[0] + m_nface_dir[1] + m_nface_dir[2])
            << " (x=" << m_nface_dir[0] << " y=" << m_nface_dir[1] << " z=" << m_nface_dir[2] << ")"
            << " area=" << m_area_total << " m2 buildings=" << m_nbld << "\n";
    for (int b = 1; b <= m_nbld; ++b) {
        Print() << "[IBSEB DEBUG]   building " << b << ": columns=" << m_bld_ncol[b]
                << " i=[" << m_bld_ilo[b] << "," << m_bld_ihi[b] << "]"
                << " j=[" << m_bld_jlo[b] << "," << m_bld_jhi[b] << "]"
                << " faces=" << m_bld_nface[b] << " area=" << m_bld_area[b] << " m2\n";
    }
    // Per-rank ownership: one line per rank with its fab ranges.
    const int nfab = static_cast<int>(m_fab_start.size()) - 1;
    const size_t bytes = static_cast<size_t>(m_nface) *
        (7 * sizeof(int) + (11 + static_cast<size_t>(nl)) * sizeof(Real));
    // The ray map: this rank's patch and its share of the column map and labels.
    size_t map_cells = static_cast<size_t>(m_patch_box.ok() ? m_patch_box.numPts() : 0);
    const size_t per_column = m_lab.ok() ? 2 : 1;      // ok() iterates the fabs itself: not inside the loop below
    for (MFIter mfi(m_col_valid); mfi.isValid(); ++mfi) {
        map_cells += per_column * static_cast<size_t>(m_col_valid[mfi].box().numPts());
    }
    const size_t map_bytes = map_cells * sizeof(int);
    std::ostringstream os;
    os << "[IBSEB DEBUG] lev=" << m_lev << " rank " << ParallelDescriptor::MyProc()
       << ": faces=" << m_nface << " over " << nfab << " fab(s), device memory "
       << bytes / 1024.0 << " kB for the faces, " << map_bytes / 1024.0 << " kB for the ray map; fab ranges:";
    for (int f = 0; f < nfab; ++f) { os << " [" << m_fab_start[f] << "," << m_fab_start[f + 1] << ")"; }
    AllPrint() << os.str() << "\n";
}

/**
 * Scatter the face count and the mean skin temperature into cell-centred
 * fields. One kernel per fab over that fab's faces (contiguous, see
 * fab_start()), with atomic adds since a corner cell collects several faces,
 * then a kernel over the box divides the temperature sum by the count.
 */
void
IBFaceSet::scatter_diagnostics (MultiFab& nfaces, MultiFab& tskin) const
{
    nfaces.setVal(0.0);
    tskin.setVal(0.0);
    const int* pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const Real* pT = d_T_skin.data();
    for (MFIter mfi(nfaces); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        auto const& n = nfaces.array(mfi);
        auto const& t = tskin.array(mfi);
        ParallelFor(f1 - f0, [=] AMREX_GPU_DEVICE (int m) noexcept {
            const int f = f0 + m;
            Gpu::Atomic::AddNoRet(&n(pi[f], pj[f], pk[f]), Real(1.0));
            Gpu::Atomic::AddNoRet(&t(pi[f], pj[f], pk[f]), pT[f]);
        });
        const Box& bx = mfi.validbox();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (n(i, j, k) > 0.0) { t(i, j, k) /= n(i, j, k); }
        });
    }
    Gpu::streamSynchronize();
}

namespace {

/**
 * Key of a face in the checkpoint records (see IBFaceSet::write_state()):
 * its fluid cell counted from the domain's low corner, then the direction
 * and the side, so every face of the level has its own key and the keys
 * sort by cell. Decoded by face_from_key().
 */
Long face_key (int ir, int jr, int kr, int dir, int side, int ny, int nz)
{
    return ((static_cast<Long>(ir) * ny + jr) * nz + kr) * 6 + dir * 2 + (side > 0 ? 1 : 0);
}

/** The (i, j, k, dir, side) of a key, (i, j, k) counted from the domain's low corner. */
std::array<int, 5> face_from_key (Long key, int ny, int nz)
{
    std::array<int, 5> f{};
    f[4] = (key % 2 == 1) ? 1 : -1;
    f[3] = static_cast<int>((key / 2) % 3);
    Long cell = key / 6;
    f[2] = static_cast<int>(cell % nz); cell /= nz;
    f[1] = static_cast<int>(cell % ny); cell /= ny;
    f[0] = static_cast<int>(cell);
    return f;
}

/** Face count of every rank on the I/O rank (empty elsewhere), and the offsets of their runs. */
void gather_counts (int nmine, std::vector<int>& counts, std::vector<int>& disp)
{
    const bool io = ParallelDescriptor::IOProcessor();
    const int nranks = ParallelDescriptor::NProcs();
    counts.assign(io ? nranks : 0, 0);
    disp.assign(io ? nranks : 0, 0);
    ParallelDescriptor::Gather(&nmine, 1, counts.data(), 1, ParallelDescriptor::IOProcessorNumber());
    if (io) {
        Long run = 0;
        for (int r = 0; r < nranks; ++r) { disp[r] = static_cast<int>(run); run += counts[r]; }
    }
}

/**
 * Word written after the header line of a face-state file and checked on
 * read: a file from a machine with another byte order or another 64-bit
 * integer width reads it back as something else and is refused.
 */
constexpr Long face_file_marker = 0x0102030405060708LL;

/** ``counts`` and ``disp`` scaled by ``nrec`` values per face (the I/O rank; empty elsewhere). */
void scale_counts (const std::vector<int>& counts, const std::vector<int>& disp, int nrec,
                   std::vector<int>& vcounts, std::vector<int>& vdisp)
{
    vcounts.resize(counts.size()); vdisp.resize(disp.size());
    for (size_t r = 0; r < counts.size(); ++r) { vcounts[r] = counts[r] * nrec; vdisp[r] = disp[r] * nrec; }
}

} // namespace

/**
 * One record per face, gathered on the I/O rank and written in key order:
 * the header line, the sorted keys, then the values record by record. Every
 * rank sorts its own faces first so the gathered runs are ordered; the I/O
 * rank then merges them by sorting an index over all the records.
 */
void
IBFaceSet::write_state (const std::string& file) const
{
    const int nl = n_layers();
    const int nrec = 2 + nl;
    const int nf = m_nface;
    const bool io = ParallelDescriptor::IOProcessor();
    const int root = ParallelDescriptor::IOProcessorNumber();

    std::vector<int> h_i(nf), h_j(nf), h_k(nf), h_dir(nf), h_side(nf);
    std::vector<Real> h_T(nf), h_H(nf), h_S(static_cast<size_t>(nf) * nl);
    Gpu::copy(Gpu::deviceToHost, d_i.begin(), d_i.end(), h_i.begin());
    Gpu::copy(Gpu::deviceToHost, d_j.begin(), d_j.end(), h_j.begin());
    Gpu::copy(Gpu::deviceToHost, d_k.begin(), d_k.end(), h_k.begin());
    Gpu::copy(Gpu::deviceToHost, d_dir.begin(), d_dir.end(), h_dir.begin());
    Gpu::copy(Gpu::deviceToHost, d_side.begin(), d_side.end(), h_side.begin());
    Gpu::copy(Gpu::deviceToHost, d_T_skin.begin(), d_T_skin.end(), h_T.begin());
    Gpu::copy(Gpu::deviceToHost, d_H.begin(), d_H.end(), h_H.begin());
    Gpu::copy(Gpu::deviceToHost, d_T_slab.begin(), d_T_slab.end(), h_S.begin());
    Gpu::streamSynchronize();

    // This rank's records in key order.
    std::vector<Long> my_key(nf);
    for (int n = 0; n < nf; ++n) {
        my_key[n] = face_key(h_i[n] - m_dom_lo[0], h_j[n] - m_dom_lo[1], h_k[n] - m_dom_lo[2],
                             h_dir[n], h_side[n], m_ny, m_nz);
    }
    std::vector<int> order(nf);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&] (int a, int b) { return my_key[a] < my_key[b]; });
    std::vector<Long> key(nf);
    std::vector<double> val(static_cast<size_t>(nf) * nrec);
    for (int m = 0; m < nf; ++m) {
        const int n = order[m];
        key[m] = my_key[n];
        double* r = &val[static_cast<size_t>(m) * nrec];
        r[0] = static_cast<double>(h_T[n]);
        r[1] = static_cast<double>(h_H[n]);
        for (int l = 0; l < nl; ++l) { r[2 + l] = static_cast<double>(h_S[static_cast<size_t>(n) * nl + l]); }
    }

    // Every rank's records on the I/O rank.
    std::vector<int> counts, disp, vcounts, vdisp;
    gather_counts(nf, counts, disp);
    scale_counts(counts, disp, nrec, vcounts, vdisp);
    Long ntot = 0;
    for (int c : counts) { ntot += c; }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ntot * nrec < std::numeric_limits<int>::max(),
                                     "erf.ibseb: too many face records for one checkpoint file");
    std::vector<Long>   all_key(io ? ntot : 0);
    std::vector<double> all_val(io ? ntot * nrec : 0);
    ParallelDescriptor::Gatherv(key.data(), nf, all_key.data(), counts, disp, root);
    ParallelDescriptor::Gatherv(val.data(), nf * nrec, all_val.data(), vcounts, vdisp, root);

    if (io) {
        std::vector<Long> perm(ntot);
        std::iota(perm.begin(), perm.end(), 0);
        std::sort(perm.begin(), perm.end(), [&] (Long a, Long b) { return all_key[a] < all_key[b]; });
        for (Long m = 1; m < ntot; ++m) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(all_key[perm[m]] != all_key[perm[m - 1]],
                                             "erf.ibseb: two faces of the level share a checkpoint key");
        }
        std::vector<Long>   skey(ntot);
        std::vector<double> sval(static_cast<size_t>(ntot) * nrec);
        for (Long m = 0; m < ntot; ++m) {
            skey[m] = all_key[perm[m]];
            std::copy_n(&all_val[static_cast<size_t>(perm[m]) * nrec], nrec, &sval[static_cast<size_t>(m) * nrec]);
        }
        std::ofstream f(file, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
        if (!f) { Abort("erf.ibseb: cannot open " + file + " for writing"); }
        f << "IBSEBFaces 1 " << ntot << " " << nl << " " << m_nx << " " << m_ny << " " << m_nz << "\n";
        const Long marker = face_file_marker;
        f.write(reinterpret_cast<const char*>(&marker), sizeof(Long));
        f.write(reinterpret_cast<const char*>(skey.data()), static_cast<std::streamsize>(ntot * sizeof(Long)));
        f.write(reinterpret_cast<const char*>(sval.data()), static_cast<std::streamsize>(sval.size() * sizeof(double)));
        if (!f) { Abort("erf.ibseb: writing " + file + " failed"); }
        if (m_params.debug) {
            Print() << "[IBSEB DEBUG] lev=" << m_lev << " face state saved: " << ntot << " records x "
                    << nrec << " values (" << ntot * nrec << " values, "
                    << ntot * (sizeof(Long) + nrec * sizeof(double)) << " bytes of records) to " << file << "\n";
        }
    }
}

/**
 * The inverse of write_state(): the I/O rank reads the table and checks the
 * header against this level (every rank gets the header, so all stop
 * together on a mismatch); each rank sends the keys of its own faces, the
 * I/O rank looks them up by binary search and sends the records back; each
 * rank then unpacks its records into the device arrays.
 */
void
IBFaceSet::read_state (const std::string& file)
{
    const int nl = n_layers();
    const int nrec = 2 + nl;
    const int nf = m_nface;
    const bool io = ParallelDescriptor::IOProcessor();
    const int root = ParallelDescriptor::IOProcessorNumber();

    // The I/O rank reads the whole file: the header (face count, slab layers,
    // domain size), the marker word and the table. What it found goes to every
    // rank as hdr[5], so all ranks stop together on any fault. The table is
    // read only when the header matches this level, so a header mismatch is
    // reported as such below rather than as a short file.
    enum : Long { ok = 1, cannot_open, bad_header, bad_marker, too_short, unsorted };
    const Long nface_level = m_nface_dir[0] + m_nface_dir[1] + m_nface_dir[2];
    std::array<Long, 6> hdr{};
    std::vector<Long>   tkey;
    std::vector<double> tval;
    if (io) {
        std::ifstream f(file, std::ifstream::in | std::ifstream::binary);
        std::string tag;
        int version = 0;
        Long marker = 0;
        if (!f) {
            hdr[5] = cannot_open;
        } else if (!(f >> tag >> version >> hdr[0] >> hdr[1] >> hdr[2] >> hdr[3] >> hdr[4])
                   || tag != "IBSEBFaces" || version != 1 || hdr[0] < 0) {
            hdr[5] = bad_header;
        } else {
            f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            f.read(reinterpret_cast<char*>(&marker), sizeof(Long));
            if (!f || marker != face_file_marker) {
                hdr[5] = bad_marker;
            } else if (hdr[1] != nl || hdr[2] != m_nx || hdr[3] != m_ny || hdr[4] != m_nz || hdr[0] != nface_level) {
                hdr[5] = ok;   // the checks below say which header entry differs
            } else {
                const Long n = hdr[0];
                tkey.resize(n);
                tval.resize(static_cast<size_t>(n) * nrec);
                f.read(reinterpret_cast<char*>(tkey.data()), static_cast<std::streamsize>(n * sizeof(Long)));
                f.read(reinterpret_cast<char*>(tval.data()), static_cast<std::streamsize>(tval.size() * sizeof(double)));
                hdr[5] = f ? ok : too_short;
                for (Long m = 1; hdr[5] == ok && m < n; ++m) {
                    if (tkey[m] <= tkey[m - 1]) { hdr[5] = unsorted; }
                }
            }
        }
    }
    ParallelDescriptor::Bcast(hdr.data(), hdr.size(), root);
    switch (hdr[5]) {
    case ok: break;
    case cannot_open: Abort("erf.ibseb: cannot open " + file + " for reading"); break;
    case bad_header:  Abort("erf.ibseb: " + file + " is not a face-state file of this ERF (its first line should read"
                            " 'IBSEBFaces 1 <faces> <layers> <nx> <ny> <nz>')"); break;
    case bad_marker:  Abort("erf.ibseb: " + file + " was written on a machine with another byte order or integer"
                            " width; restart from a checkpoint written on this kind of machine"); break;
    case too_short:   Abort("erf.ibseb: " + file + " ends before the " + std::to_string(hdr[0])
                            + " records its header announces"); break;
    default:          Abort("erf.ibseb: the face keys of " + file + " are not in ascending order"); break;
    }
    if (hdr[1] != nl) {
        Abort("erf.ibseb: " + file + " holds " + std::to_string(hdr[1]) + " slab layers per face; the deck sets"
              " erf.ibseb.n_slab_layers = " + std::to_string(nl) + ". Restart with the checkpoint's value.");
    }
    if (hdr[2] != m_nx || hdr[3] != m_ny || hdr[4] != m_nz) {
        Abort("erf.ibseb: " + file + " was written for a different domain (" + std::to_string(hdr[2]) + " x "
              + std::to_string(hdr[3]) + " x " + std::to_string(hdr[4]) + " cells against " + std::to_string(m_nx)
              + " x " + std::to_string(m_ny) + " x " + std::to_string(m_nz) + " here)");
    }
    if (hdr[0] != nface_level) {
        Abort("erf.ibseb: " + file + " was written for a different building layout (" + std::to_string(hdr[0])
              + " faces against the " + std::to_string(nface_level)
              + " the blanking gives); restart from a checkpoint of the same buildings");
    }

    const Long ntot = hdr[0];

    // This rank's keys to the I/O rank.
    std::vector<int> h_i(nf), h_j(nf), h_k(nf), h_dir(nf), h_side(nf);
    Gpu::copy(Gpu::deviceToHost, d_i.begin(), d_i.end(), h_i.begin());
    Gpu::copy(Gpu::deviceToHost, d_j.begin(), d_j.end(), h_j.begin());
    Gpu::copy(Gpu::deviceToHost, d_k.begin(), d_k.end(), h_k.begin());
    Gpu::copy(Gpu::deviceToHost, d_dir.begin(), d_dir.end(), h_dir.begin());
    Gpu::copy(Gpu::deviceToHost, d_side.begin(), d_side.end(), h_side.begin());
    Gpu::streamSynchronize();
    std::vector<Long> key(nf);
    for (int n = 0; n < nf; ++n) {
        key[n] = face_key(h_i[n] - m_dom_lo[0], h_j[n] - m_dom_lo[1], h_k[n] - m_dom_lo[2],
                          h_dir[n], h_side[n], m_ny, m_nz);
    }
    std::vector<int> counts, disp, vcounts, vdisp;
    gather_counts(nf, counts, disp);
    scale_counts(counts, disp, nrec, vcounts, vdisp);
    std::vector<Long> all_key(io ? ntot : 0);
    ParallelDescriptor::Gatherv(key.data(), nf, all_key.data(), counts, disp, root);

    // Look every face up; a face without a record means other buildings.
    std::vector<double> all_val(io ? ntot * nrec : 0);
    std::array<Long, 2> missing{0, -1};   // count, first key
    if (io) {
        for (Long m = 0; m < ntot; ++m) {
            const auto it = std::lower_bound(tkey.begin(), tkey.end(), all_key[m]);
            if (it == tkey.end() || *it != all_key[m]) {
                if (missing[0] == 0) { missing[1] = all_key[m]; }
                ++missing[0];
                continue;
            }
            const Long t = it - tkey.begin();
            std::copy_n(&tval[static_cast<size_t>(t) * nrec], nrec, &all_val[static_cast<size_t>(m) * nrec]);
        }
    }
    ParallelDescriptor::Bcast(missing.data(), missing.size(), root);
    if (missing[0] > 0) {
        const auto fc = face_from_key(missing[1], m_ny, m_nz);
        Abort("erf.ibseb: " + file + " was written for a different building layout (" + std::to_string(missing[0])
              + " faces of this blanking have no record there, the first at i=" + std::to_string(fc[0] + m_dom_lo[0])
              + " j=" + std::to_string(fc[1] + m_dom_lo[1]) + " k=" + std::to_string(fc[2] + m_dom_lo[2])
              + " dir=" + std::to_string(fc[3]) + " side=" + std::to_string(fc[4])
              + "); restart from a checkpoint of the same buildings");
    }

    // The records back to the ranks that asked for them, in the order they asked.
    std::vector<double> val(static_cast<size_t>(nf) * nrec);
#ifdef AMREX_USE_MPI
    BL_MPI_REQUIRE(MPI_Scatterv(all_val.data(), vcounts.data(), vdisp.data(), MPI_DOUBLE,
                                val.data(), nf * nrec, MPI_DOUBLE, root, ParallelDescriptor::Communicator()));
#else
    val = all_val;
#endif
    std::vector<Real> h_T(nf), h_H(nf), h_S(static_cast<size_t>(nf) * nl);
    for (int n = 0; n < nf; ++n) {
        const double* r = &val[static_cast<size_t>(n) * nrec];
        h_T[n] = static_cast<Real>(r[0]);
        h_H[n] = static_cast<Real>(r[1]);
        for (int l = 0; l < nl; ++l) { h_S[static_cast<size_t>(n) * nl + l] = static_cast<Real>(r[2 + l]); }
    }
    Gpu::copy(Gpu::hostToDevice, h_T.begin(), h_T.end(), d_T_skin.begin());
    Gpu::copy(Gpu::hostToDevice, h_H.begin(), h_H.end(), d_H.begin());
    Gpu::copy(Gpu::hostToDevice, h_S.begin(), h_S.end(), d_T_slab.begin());
    Gpu::streamSynchronize();
    if (m_params.debug) {
        Real tmin = 1.0e30, tmax = -1.0e30;
        for (int n = 0; n < nf; ++n) { tmin = std::min(tmin, h_T[n]); tmax = std::max(tmax, h_T[n]); }
        ParallelDescriptor::ReduceRealMin(tmin);
        ParallelDescriptor::ReduceRealMax(tmax);
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " face state loaded: " << ntot << " records from " << file
                << ", T_skin in [" << tmin << ", " << tmax << "] K\n";
    }
}

/**
 * Summary line and per-building CSV rows. The dynamic quantities (skin
 * temperature range and area-weighted mean per building) are gathered on the
 * host from the device arrays and reduced; the static counts and areas were
 * reduced in build().
 */
void
IBFaceSet::report (Real time, int step, bool write_csv) const
{
    // Skin temperature: min, max and per-building mean over all ranks.
    std::vector<Real> h_T(m_nface), h_A(m_nface);
    std::vector<int>  h_b(m_nface);
    Gpu::copy(Gpu::deviceToHost, d_T_skin.begin(), d_T_skin.end(), h_T.begin());
    Gpu::copy(Gpu::deviceToHost, d_area.begin(),   d_area.end(),   h_A.begin());
    Gpu::copy(Gpu::deviceToHost, d_bid.begin(),    d_bid.end(),    h_b.begin());
    Gpu::streamSynchronize();

    std::vector<Real> h_S(m_nface), h_sh(m_nface), h_L(m_nface), h_H(m_nface), h_G(m_nface), h_r(m_nface), h_Q(m_nface), h_w(m_nface);
    Gpu::copy(Gpu::deviceToHost, d_w_star.begin(), d_w_star.end(), h_w.begin());
    Gpu::copy(Gpu::deviceToHost, d_resid.begin(), d_resid.end(), h_r.begin());
    Gpu::copy(Gpu::deviceToHost, d_Q_ext.begin(), d_Q_ext.end(), h_Q.begin());
    Gpu::copy(Gpu::deviceToHost, d_SW_abs.begin(), d_SW_abs.end(), h_S.begin());
    Gpu::copy(Gpu::deviceToHost, d_shadow.begin(), d_shadow.end(), h_sh.begin());
    Gpu::copy(Gpu::deviceToHost, d_LW_net.begin(), d_LW_net.end(), h_L.begin());
    Gpu::copy(Gpu::deviceToHost, d_H.begin(), d_H.end(), h_H.begin());
    Gpu::copy(Gpu::deviceToHost, d_G.begin(), d_G.end(), h_G.begin());
    Gpu::streamSynchronize();

    std::vector<Real> bsum(m_nbld + 1, 0.0), bsw(m_nbld + 1, 0.0), bsh(m_nbld + 1, 0.0), blw(m_nbld + 1, 0.0), bH(m_nbld + 1, 0.0);
    std::vector<Real> bG(m_nbld + 1, 0.0), bQ(m_nbld + 1, 0.0), btmin(m_nbld + 1, 1.0e30), btmax(m_nbld + 1, -1.0e30), bres(m_nbld + 1, 0.0);
    Real tmin = 1.0e30, tmax = -1.0e30, sw_sum = 0.0, sw_max = 0.0, sh_area = 0.0, lw_sum = 0.0, H_sum = 0.0, G_sum = 0.0;
    Real res_max = 0.0, Q_sum = 0.0, w_max = 0.0;
    for (int n = 0; n < m_nface; ++n) {
        w_max = std::max(w_max, h_w[n]);
        G_sum += h_G[n] * h_A[n];
        Q_sum += h_Q[n] * h_A[n];
        res_max = std::max(res_max, h_r[n]);
        bG[h_b[n]]   += h_G[n] * h_A[n];
        bQ[h_b[n]]   += h_Q[n] * h_A[n];
        btmin[h_b[n]] = std::min(btmin[h_b[n]], h_T[n]);
        btmax[h_b[n]] = std::max(btmax[h_b[n]], h_T[n]);
        bres[h_b[n]]  = std::max(bres[h_b[n]], h_r[n]);
        bsum[h_b[n]] += h_T[n] * h_A[n];
        bsw[h_b[n]]  += h_S[n] * h_A[n];
        bsh[h_b[n]]  += h_sh[n] * h_A[n];
        blw[h_b[n]]  += h_L[n] * h_A[n];
        bH[h_b[n]]   += h_H[n] * h_A[n];
        tmin = std::min(tmin, h_T[n]);
        tmax = std::max(tmax, h_T[n]);
        sw_sum  += h_S[n] * h_A[n];
        sw_max   = std::max(sw_max, h_S[n]);
        sh_area += h_sh[n] * h_A[n];
        lw_sum  += h_L[n] * h_A[n];
        H_sum   += h_H[n] * h_A[n];
    }
    ParallelDescriptor::ReduceRealSum(bsum.data(), m_nbld + 1);
    ParallelDescriptor::ReduceRealSum(bsw.data(),  m_nbld + 1);
    ParallelDescriptor::ReduceRealSum(bsh.data(),  m_nbld + 1);
    ParallelDescriptor::ReduceRealSum(blw.data(),  m_nbld + 1);
    ParallelDescriptor::ReduceRealSum(bH.data(),   m_nbld + 1);
    ParallelDescriptor::ReduceRealSum(bG.data(),   m_nbld + 1);
    ParallelDescriptor::ReduceRealSum(bQ.data(),   m_nbld + 1);
    ParallelDescriptor::ReduceRealMin(btmin.data(), m_nbld + 1);
    ParallelDescriptor::ReduceRealMax(btmax.data(), m_nbld + 1);
    ParallelDescriptor::ReduceRealMax(bres.data(),  m_nbld + 1);
    ParallelDescriptor::ReduceRealMax(res_max);
    ParallelDescriptor::ReduceRealMax(w_max);
    ParallelDescriptor::ReduceRealSum(Q_sum);
    ParallelDescriptor::ReduceRealMin(tmin);
    ParallelDescriptor::ReduceRealMax(tmax);
    ParallelDescriptor::ReduceRealSum(sw_sum);
    ParallelDescriptor::ReduceRealMax(sw_max);
    ParallelDescriptor::ReduceRealSum(sh_area);
    ParallelDescriptor::ReduceRealSum(lw_sum);
    ParallelDescriptor::ReduceRealSum(H_sum);
    ParallelDescriptor::ReduceRealSum(G_sum);
    const Real G_mean   = (m_area_total > 0.0) ? G_sum / m_area_total : 0.0;
    const Real sw_mean  = (m_area_total > 0.0) ? sw_sum / m_area_total : 0.0;
    const Real sh_frac  = (m_area_total > 0.0) ? sh_area / m_area_total : 0.0;
    const Real lw_mean  = (m_area_total > 0.0) ? lw_sum / m_area_total : 0.0;
    const Real H_mean   = (m_area_total > 0.0) ? H_sum / m_area_total : 0.0;

    const Long nf_all = m_nface_dir[0] + m_nface_dir[1] + m_nface_dir[2];
    Print() << "[IBSEB] lev=" << m_lev << " step=" << step << " t=" << time
            << " faces=" << nf_all
            << " (x=" << m_nface_dir[0] << " y=" << m_nface_dir[1] << " z=" << m_nface_dir[2] << ")"
            << " buildings=" << m_nbld
            << " area=" << std::setprecision(10) << m_area_total << " m2"
            << " T_skin_min=" << ((nf_all > 0) ? tmin : Real(0.0))
            << " T_skin_max=" << ((nf_all > 0) ? tmax : Real(0.0))
            << " SW_abs_mean=" << sw_mean << " SW_abs_max=" << sw_max
            << " shadow_frac=" << sh_frac << " LW_net_mean=" << lw_mean
            << " H_mean=" << H_mean << " G_mean=" << G_mean
            << " Q_ext_mean=" << ((m_area_total > 0.0) ? Q_sum / m_area_total : Real(0.0))
            << " resid_max=" << res_max
            << " w_star_max=" << w_max
            << " H_total_W=" << std::setprecision(12) << H_sum << "\n";
    // Cost line: the slowest rank's per-step time of the balance and its face
    // count, so the cost of the face list on a larger case can be estimated.
    {
        Real ms = (m_cost_n > 0) ? 1000.0 * m_cost_s / m_cost_n : 0.0;
        Real ms_max = ms; int nf_max = m_nface;
        ParallelDescriptor::ReduceRealMax(ms_max);
        ParallelDescriptor::ReduceIntMax(nf_max);
        Print() << "[IBSEB] cost: lev=" << m_lev << " advance_ms_per_step_max=" << std::setprecision(4) << ms_max
                << " faces_per_rank_max=" << nf_max << " ranks=" << ParallelDescriptor::NProcs()
                << " init_s=" << m_init_cost_s << "\n";
    }

    if (m_params.debug) {
        for (int b = 1; b <= m_nbld; ++b) {
            const Real tmean = (m_bld_area[b] > 0.0) ? bsum[b] / m_bld_area[b] : Real(0.0);
            const Real swm = (m_bld_area[b] > 0.0) ? bsw[b] / m_bld_area[b] : Real(0.0);
            const Real shf = (m_bld_area[b] > 0.0) ? bsh[b] / m_bld_area[b] : Real(0.0);
            const Real lwm = (m_bld_area[b] > 0.0) ? blw[b] / m_bld_area[b] : Real(0.0);
            const Real Hm  = (m_bld_area[b] > 0.0) ? bH[b]  / m_bld_area[b] : Real(0.0);
            const Real Gm  = (m_bld_area[b] > 0.0) ? bG[b]  / m_bld_area[b] : Real(0.0);
            Print() << "[IBSEB DEBUG]   building " << b << ": faces=" << m_bld_nface[b]
                    << " area=" << m_bld_area[b] << " m2 T_skin_mean=" << tmean << " K"
                    << " [" << btmin[b] << ", " << btmax[b] << "]"
                    << " SW_abs_mean=" << swm << " W/m2 shadow_frac=" << shf
                    << " LW_net_mean=" << lwm << " W/m2 H_mean=" << Hm << " W/m2 G_mean=" << Gm
                    << " W/m2 resid_max=" << bres[b] << " W/m2\n";
        }
    }

    if (write_csv && !m_params.dump_faces_file.empty()) {
        std::string prefix = m_params.dump_faces_file;
        // Each level has its own list; a refined level tags its files so it
        // does not overwrite the dump of the level below.
        if (m_lev > 0) { prefix += ".lev" + std::to_string(m_lev); }
        if (m_params.dump_faces_tag_step) {
            std::ostringstream os;
            os << prefix << ".step" << std::setw(6) << std::setfill('0') << step;
            prefix = os.str();
        }
        dump_faces(prefix);
    }

    if (!write_csv || !ParallelDescriptor::IOProcessor()) { return; }
    const std::string header =
        "time_s,step,level,building,n_faces,area_m2,T_skin_mean_K,SW_abs_mean_Wm2,shadow_frac,LW_net_mean_Wm2,H_mean_Wm2,"
        "G_mean_Wm2,Q_ext_mean_Wm2,T_skin_min_K,T_skin_max_K,resid_max_Wm2,"
        "sun_zenith_deg,sun_azimuth_deg,dni_Wm2,diffuse_h_Wm2,building_level0";
    bool need_header = true;
    {
        // A restarted run appends; a file another version of the report wrote
        // (other columns) would take rows its header does not describe.
        std::ifstream probe(m_params.csv_file);
        std::string first;
        if (probe.good() && std::getline(probe, first)) {
            if (!first.empty() && first.back() == '\r') { first.pop_back(); }
            if (first != header) {
                Abort("erf.ibseb: " + m_params.csv_file + " has another set of columns than this version writes ("
                      + first.substr(0, 60) + "...); move it aside or set erf.ibseb.csv_file to a new name");
            }
            need_header = false;
        }
    }
    std::ofstream csv(m_params.csv_file, std::ios::app);
    if (need_header) { csv << header << "\n"; }
    for (int b = 1; b <= m_nbld; ++b) {
        const Real tmean = (m_bld_area[b] > 0.0) ? bsum[b] / m_bld_area[b] : Real(0.0);
        const Real swm   = (m_bld_area[b] > 0.0) ? bsw[b]  / m_bld_area[b] : Real(0.0);
        const Real shf   = (m_bld_area[b] > 0.0) ? bsh[b]  / m_bld_area[b] : Real(0.0);
        const Real lwm   = (m_bld_area[b] > 0.0) ? blw[b]  / m_bld_area[b] : Real(0.0);
        const Real Hm    = (m_bld_area[b] > 0.0) ? bH[b]   / m_bld_area[b] : Real(0.0);
        const Real Gm    = (m_bld_area[b] > 0.0) ? bG[b]   / m_bld_area[b] : Real(0.0);
        const Real Qm    = (m_bld_area[b] > 0.0) ? bQ[b]   / m_bld_area[b] : Real(0.0);
        csv << std::setprecision(10) << time << "," << step << "," << m_lev << "," << b << ","
            << m_bld_nface[b] << "," << m_bld_area[b] << "," << tmean << "," << swm << "," << shf << "," << lwm << "," << Hm << ","
            << Gm << "," << Qm << "," << btmin[b] << "," << btmax[b] << "," << bres[b] << ","
            << m_sun.zenith * 180.0 / PI << "," << m_sun.azimuth * 180.0 / PI << "," << m_sun.dni << "," << m_sun.diffuse_h << ","
            << m_to_lev0[b] << "\n";
    }
}

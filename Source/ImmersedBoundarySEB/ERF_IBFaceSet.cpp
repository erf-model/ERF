/**
 * \file ERF_IBFaceSet.cpp
 * \brief Detection, storage, output and reporting of the wall faces of
 *        resolved buildings (phase 1 of the immersed-boundary surface energy
 *        balance). The conventions are documented in ERF_IBFaceSet.H.
 */
#include "ERF_IBFaceSet.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Gpu.H>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>

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

} // namespace

/**
 * Detect the faces owned by this rank and allocate the per-face arrays.
 *
 * The scan visits every valid cell of every local fab in MFIter order. A
 * fluid cell (blanking < 0.5) contributes one face for each of its six
 * neighbours that is solid (blanking >= 0.5); the neighbour may sit in a
 * ghost cell, which is why the blanking must have its ghost cells filled.
 * Neighbours outside the domain in a non-periodic direction are skipped, so
 * a building against a non-periodic boundary has no face there. Solid cells
 * contribute nothing but the column mask used for the building ids.
 *
 * Faces are appended in scan order, which makes them contiguous per fab;
 * m_fab_start records where each fab's faces begin. The per-direction,
 * per-building and total counts are reduced over all ranks here, once, so
 * report() does not have to reduce static numbers every time.
 */
void
IBFaceSet::build (const MultiFab& blanking, const Geometry& geom)
{
    const Box& domain = geom.Domain();
    const int nx  = domain.length(0);
    const int ny  = domain.length(1);
    const int ilo = domain.smallEnd(0);
    const int jlo = domain.smallEnd(1);
    const auto dx = geom.CellSizeArray();
    const Real face_area[3] = { dx[1] * dx[2], dx[0] * dx[2], dx[0] * dx[1] };

    std::vector<int>  h_i, h_j, h_k, h_dir, h_side, h_nbi, h_nbj;
    std::vector<Real> h_area;
    // Largest blanking in each column of this rank; reduced below so every
    // rank labels the same solid columns.
    std::vector<Real> colmax(static_cast<size_t>(nx) * ny, 0.0);

    m_fab_start.clear();
    for (MFIter mfi(blanking); mfi.isValid(); ++mfi) {
        m_fab_start.push_back(static_cast<int>(h_i.size()));
        const Box& bx = mfi.validbox();
        const HostFab hf(blanking[mfi]);
        auto const& b = hf.array();
        const Box& fbox = Box(b);
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if (b(i, j, k) >= 0.5) {
                Real& cm = colmax[static_cast<size_t>(i - ilo) * ny + (j - jlo)];
                cm = std::max(cm, b(i, j, k));
                continue;
            }
            for (int d = 0; d < 3; ++d) {
                for (int s = -1; s <= 1; s += 2) {
                    IntVect nb(i, j, k);
                    nb[d] += s;
                    if (!fbox.contains(nb)) { continue; }
                    if (!geom.isPeriodic(d) && !domain.contains(nb)) { continue; }
                    if (b(nb) < 0.5) { continue; }
                    h_i.push_back(i); h_j.push_back(j); h_k.push_back(k);
                    h_dir.push_back(d); h_side.push_back(s);
                    // Column of the solid neighbour, wrapped for periodic directions.
                    int ci = nb[0] - ilo, cj = nb[1] - jlo;
                    ci = (ci % nx + nx) % nx;
                    cj = (cj % ny + ny) % ny;
                    h_nbi.push_back(ci); h_nbj.push_back(cj);
                    h_area.push_back(face_area[d]);
                }
            }
        }}}
    }
    m_fab_start.push_back(static_cast<int>(h_i.size()));
    m_nface = static_cast<int>(h_i.size());

    // Buildings: 4-connected solid columns, numbered in scan order. The column
    // mask is reduced so every rank labels the same columns with the same ids;
    // the labelling itself is a plain depth-first flood fill on the host.
    ParallelDescriptor::ReduceRealMax(colmax.data(), nx * ny);
    std::vector<int> label(static_cast<size_t>(nx) * ny, 0);
    std::vector<int> stack;
    m_nbld = 0;
    for (int ci = 0; ci < nx; ++ci) {
        for (int cj = 0; cj < ny; ++cj) {
            const size_t p = static_cast<size_t>(ci) * ny + cj;
            if (colmax[p] < 0.5 || label[p] > 0) { continue; }
            ++m_nbld;
            label[p] = m_nbld;
            stack.clear();
            stack.push_back(static_cast<int>(p));
            while (!stack.empty()) {
                const int q = stack.back(); stack.pop_back();
                const int qi = q / ny, qj = q % ny;
                const int ni[4] = {qi - 1, qi + 1, qi, qi};
                const int nj[4] = {qj, qj, qj - 1, qj + 1};
                for (int n = 0; n < 4; ++n) {
                    if (ni[n] < 0 || ni[n] >= nx || nj[n] < 0 || nj[n] >= ny) { continue; }
                    const size_t r = static_cast<size_t>(ni[n]) * ny + nj[n];
                    if (colmax[r] >= 0.5 && label[r] == 0) {
                        label[r] = m_nbld;
                        stack.push_back(static_cast<int>(r));
                    }
                }
            }
        }
    }
    std::vector<int> h_bid(m_nface);
    for (int n = 0; n < m_nface; ++n) {
        h_bid[n] = label[static_cast<size_t>(h_nbi[n]) * ny + h_nbj[n]];
    }
    // Footprints, for the debug summary: columns and bounding box per building.
    m_bld_ncol.assign(m_nbld + 1, 0);
    m_bld_ilo.assign(m_nbld + 1, nx); m_bld_ihi.assign(m_nbld + 1, -1);
    m_bld_jlo.assign(m_nbld + 1, ny); m_bld_jhi.assign(m_nbld + 1, -1);
    for (int ci = 0; ci < nx; ++ci) {
        for (int cj = 0; cj < ny; ++cj) {
            const int b = label[static_cast<size_t>(ci) * ny + cj];
            if (b == 0) { continue; }
            m_bld_ncol[b] += 1;
            m_bld_ilo[b] = std::min(m_bld_ilo[b], ci + ilo); m_bld_ihi[b] = std::max(m_bld_ihi[b], ci + ilo);
            m_bld_jlo[b] = std::min(m_bld_jlo[b], cj + jlo); m_bld_jhi[b] = std::max(m_bld_jhi[b], cj + jlo);
        }
    }

    // Static totals: per direction, per building, and the summed area.
    std::vector<long> nd(3, 0);
    std::vector<long> bn(m_nbld + 1, 0);
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

    // Device arrays. The state starts uniform at T_skin_init (the slab too:
    // phase 5 sets the interior boundary); the view fractions and fluxes start
    // at zero and are filled by the later phases.
    const size_t nf = static_cast<size_t>(m_nface);
    upload(d_i, h_i); upload(d_j, h_j); upload(d_k, h_k);
    upload(d_dir, h_dir); upload(d_side, h_side); upload(d_bid, h_bid);
    upload(d_area, h_area);
    fill(d_mat,      nf, 0);
    fill(d_T_skin,   nf, m_params.T_skin_init);
    fill(d_T_slab,   nf * static_cast<size_t>(n_layers()), m_params.T_skin_init);
    fill(d_f_sky,    nf, Real(0.0));
    fill(d_f_ground, nf, Real(0.0));
    fill(d_f_bldg,   nf, Real(0.0));
    fill(d_SW_abs,   nf, Real(0.0));
    fill(d_LW_net,   nf, Real(0.0));
    fill(d_H,        nf, Real(0.0));
    fill(d_LE,       nf, Real(0.0));
    fill(d_G,        nf, Real(0.0));
    fill(d_Q_ext,    nf, Real(0.0));
    Gpu::streamSynchronize();

    if (m_params.debug) { print_debug_summary(); }
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
    std::ostringstream os;
    os << "[IBSEB DEBUG] lev=" << m_lev << " rank " << ParallelDescriptor::MyProc()
       << ": faces=" << m_nface << " over " << nfab << " fab(s), device memory "
       << bytes / 1024.0 << " kB; fab ranges:";
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

/**
 * Write the skin and slab temperatures into the six-slot state field. The
 * slot of a face is ``(dir*2 + (side>0)) * (1 + n_layers)``; no two faces of a
 * cell share a slot, so the writes are plain stores.
 */
void
IBFaceSet::save_state (MultiFab& state) const
{
    AMREX_ALWAYS_ASSERT(state.nComp() == state_ncomp());
    state.setVal(0.0);
    const int nl = n_layers();
    const int* pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const int* pd = d_dir.data(); const int* ps = d_side.data();
    const Real* pT = d_T_skin.data(); const Real* pS = d_T_slab.data();
    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        auto const& st = state.array(mfi);
        ParallelFor(f1 - f0, [=] AMREX_GPU_DEVICE (int m) noexcept {
            const int f  = f0 + m;
            const int c0 = (pd[f] * 2 + (ps[f] > 0 ? 1 : 0)) * (1 + nl);
            st(pi[f], pj[f], pk[f], c0) = pT[f];
            for (int l = 0; l < nl; ++l) { st(pi[f], pj[f], pk[f], c0 + 1 + l) = pS[f * nl + l]; }
        });
    }
    Gpu::streamSynchronize();
    if (m_params.debug) {
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " face state saved: " << state_ncomp()
                << " components (6 slots x " << (1 + n_layers()) << ")\n";
    }
}

/**
 * Read the skin and slab temperatures back from the six-slot state field,
 * the inverse of save_state(), into a list that build() has just recreated.
 */
void
IBFaceSet::load_state (const MultiFab& state)
{
    AMREX_ALWAYS_ASSERT(state.nComp() == state_ncomp());
    const int nl = n_layers();
    const int* pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const int* pd = d_dir.data(); const int* ps = d_side.data();
    Real* pT = d_T_skin.data(); Real* pS = d_T_slab.data();
    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        auto const& st = state.const_array(mfi);
        ParallelFor(f1 - f0, [=] AMREX_GPU_DEVICE (int m) noexcept {
            const int f  = f0 + m;
            const int c0 = (pd[f] * 2 + (ps[f] > 0 ? 1 : 0)) * (1 + nl);
            pT[f] = st(pi[f], pj[f], pk[f], c0);
            for (int l = 0; l < nl; ++l) { pS[f * nl + l] = st(pi[f], pj[f], pk[f], c0 + 1 + l); }
        });
    }
    Gpu::streamSynchronize();
    if (m_params.debug) {
        std::vector<Real> h_T(m_nface);
        Gpu::copy(Gpu::deviceToHost, d_T_skin.begin(), d_T_skin.end(), h_T.begin());
        Gpu::streamSynchronize();
        Real tmin = 1.0e30, tmax = -1.0e30;
        for (int n = 0; n < m_nface; ++n) { tmin = std::min(tmin, h_T[n]); tmax = std::max(tmax, h_T[n]); }
        ParallelDescriptor::ReduceRealMin(tmin);
        ParallelDescriptor::ReduceRealMax(tmax);
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " face state loaded: T_skin in ["
                << tmin << ", " << tmax << "] K\n";
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

    std::vector<Real> bsum(m_nbld + 1, 0.0);
    Real tmin = 1.0e30, tmax = -1.0e30;
    for (int n = 0; n < m_nface; ++n) {
        bsum[h_b[n]] += h_T[n] * h_A[n];
        tmin = std::min(tmin, h_T[n]);
        tmax = std::max(tmax, h_T[n]);
    }
    ParallelDescriptor::ReduceRealSum(bsum.data(), m_nbld + 1);
    ParallelDescriptor::ReduceRealMin(tmin);
    ParallelDescriptor::ReduceRealMax(tmax);

    const long nf_all = m_nface_dir[0] + m_nface_dir[1] + m_nface_dir[2];
    Print() << "[IBSEB] lev=" << m_lev << " step=" << step << " t=" << time
            << " faces=" << nf_all
            << " (x=" << m_nface_dir[0] << " y=" << m_nface_dir[1] << " z=" << m_nface_dir[2] << ")"
            << " buildings=" << m_nbld
            << " area=" << std::setprecision(10) << m_area_total << " m2"
            << " T_skin_min=" << ((nf_all > 0) ? tmin : Real(0.0))
            << " T_skin_max=" << ((nf_all > 0) ? tmax : Real(0.0)) << "\n";

    if (m_params.debug) {
        for (int b = 1; b <= m_nbld; ++b) {
            const Real tmean = (m_bld_area[b] > 0.0) ? bsum[b] / m_bld_area[b] : Real(0.0);
            Print() << "[IBSEB DEBUG]   building " << b << ": faces=" << m_bld_nface[b]
                    << " area=" << m_bld_area[b] << " m2 T_skin_mean=" << tmean << " K\n";
        }
    }

    if (!write_csv || !ParallelDescriptor::IOProcessor()) { return; }
    bool need_header = true;
    {
        std::ifstream probe(m_params.csv_file, std::ios::ate);
        if (probe.good() && probe.tellg() > 0) { need_header = false; }
    }
    std::ofstream csv(m_params.csv_file, std::ios::app);
    if (need_header) { csv << "time_s,step,level,building,n_faces,area_m2,T_skin_mean_K\n"; }
    for (int b = 1; b <= m_nbld; ++b) {
        const Real tmean = (m_bld_area[b] > 0.0) ? bsum[b] / m_bld_area[b] : Real(0.0);
        csv << std::setprecision(10) << time << "," << step << "," << m_lev << "," << b << ","
            << m_bld_nface[b] << "," << m_bld_area[b] << "," << tmean << "\n";
    }
}

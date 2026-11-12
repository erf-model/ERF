#include <cmath>
#include <map>
#include <utility>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <gtest/gtest.h>

#include "ERF_IBFaceSet.H"
#include "ERF_IBSEBParams.H"

// The roof face of each column under a building (erf.ibseb.two_stream_roofs,
// IBFaceSet::fill_roofs()): given the cell on each column's roof, the upward face on that
// cell writes its albedo, emissivity and skin temperature for the column. Buildings of
// three heights in four boxes, one across the periodic edge, one under an overhang (two
// upward faces in one column: only the one on the building that stands on the ground
// counts), and one that fills the column (no roof face). Every face has its own values, so
// a face written for another column, or another face of the same column, gives another number.

namespace {

struct Cells { int ilo, ihi, jlo, jhi, klo, khi; };

void fill_blanking (amrex::MultiFab& b, const amrex::Geometry& geom, const std::vector<Cells>& solid)
{
    b.setVal(0.0);
    for (amrex::MFIter mfi(b); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& a = b.array(mfi);
        for (const Cells& s : solid) {
            const amrex::Box ov = amrex::Box(amrex::IntVect(s.ilo, s.jlo, s.klo), amrex::IntVect(s.ihi, s.jhi, s.khi)) & bx;
            if (ov.isEmpty()) { continue; }
            amrex::ParallelFor(ov, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { a(i, j, k) = 1.0; });
        }
    }
    b.FillBoundary(geom.periodicity());
}

// The cell on each column's roof, as TwoStreamRadiation::roof_cell() has it: the first
// cell above the solid cells that stand on the ground.
void fill_roof_cell (amrex::iMultiFab& rk, const amrex::MultiFab& blank, int k_top)
{
    for (amrex::MFIter mfi(rk); mfi.isValid(); ++mfi) {
        const auto r = rk.array(mfi);
        const auto b = blank.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
            int k = 0;
            while (k <= k_top && b(i, j, k) >= amrex::Real(0.5)) { ++k; }
            r(i, j, 0) = k;
        });
    }
    amrex::Gpu::streamSynchronize();
}

template <class T>
std::vector<T> host (const amrex::Gpu::DeviceVector<T>& d)
{
    std::vector<T> h(d.size());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d.begin(), d.end(), h.begin());
    amrex::Gpu::streamSynchronize();
    return h;
}

template <class T>
void to_device (amrex::Gpu::DeviceVector<T>& d, const std::vector<T>& h)
{
    ASSERT_EQ(d.size(), h.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h.begin(), h.end(), d.begin());
    amrex::Gpu::streamSynchronize();
}

} // namespace

TEST(IBSEBRoofs, EachColumnGetsTheRoofThatStandsOnIt)
{
    using namespace amrex;
    const int nx = 16, ny = 16, nz = 8;
    const RealBox rb({0.0, 0.0, 0.0}, {160.0, 160.0, 80.0});
    const Box dom(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, nz - 1));
    const Geometry geom(dom, rb, 0, {1, 1, 0});
    BoxArray ba(dom);
    ba.maxSize(IntVect(8, 8, 8));
    const DistributionMapping dm(ba);
    MultiFab blank(ba, dm, 1, 1);
    // A block of four cells' height; a tower across the periodic edge in x; a low building
    // (two cells) with a slab over it from cell 3 to 4, one cell of air between them: that
    // cell carries the low roof's upward face and the slab's downward one; a slab over open
    // ground; a building that fills the column.
    const std::vector<Cells> solid = {{6, 9, 7, 8, 0, 3}, {0, 0, 2, 4, 0, 5}, {15, 15, 2, 4, 0, 5}, {12, 13, 12, 12, 0, 1},
                                      {12, 13, 12, 12, 3, 4}, {3, 4, 12, 13, 4, 5}, {10, 10, 1, 1, 0, 7}};
    fill_blanking(blank, geom, solid);
    BoxList bl = ba.boxList();
    for (Box& b : bl) { b.setRange(2, 0); }
    const BoxArray ba2d(std::move(bl));
    iMultiFab roof_cell(ba2d, dm, 1, 0);
    fill_roof_cell(roof_cell, blank, nz - 1);

    IBSEBParams params;
    params.enable = true;
    params.view_n_az = 8;
    params.view_n_el = 4;
    // The sun one degree from the zenith: no building shades another's roof; only the slab
    // shades the roof under it.
    params.sun_mode = "fixed";
    params.sun_zenith_deg = 1.0;
    params.sun_azimuth_deg = 180.0;
    IBFaceSet faces(params, 0);
    faces.build(blank, geom);
    faces.release_labels();
    faces.compute_view_fractions();
    faces.compute_shortwave(Real(0.0));   // the fixed sun: the same at every time
    const int nf = faces.n_faces();
    const auto fi = host(faces.d_i), fj = host(faces.d_j), fk = host(faces.d_k), fd = host(faces.d_dir), fs = host(faces.d_side);
    std::vector<Real> alb(nf), eps(nf), T(nf);
    for (int f = 0; f < nf; ++f) {
        const Real u = Real(0.37) * fi[f] + Real(0.61) * fj[f] + Real(1.3) * fk[f] + Real(2.1) * fd[f] + Real(0.7) * fs[f];
        alb[f] = Real(0.1) + Real(0.01) * std::fmod(u, Real(40.0));
        eps[f] = Real(0.6) + Real(0.01) * std::fmod(Real(3.0) * u, Real(35.0));
        T[f]   = Real(280.0) + std::fmod(Real(7.0) * u, Real(50.0));
    }
    to_device(faces.d_albedo, alb); to_device(faces.d_emis, eps); to_device(faces.d_T_skin, T);

    MultiFab out(ba2d, dm, 5, 0);
    out.setVal(-5.0);   // every column is set, also where there is no roof
    faces.fill_roofs(out, roof_cell, true, Real(0.0));

    // What each column should get, from the faces and the geometry.
    std::map<std::pair<int, int>, int> top;      // cell on the roof of the ground-standing building
    for (int i = 0; i < nx; ++i) { for (int j = 0; j < ny; ++j) { top[{i, j}] = 0; } }
    for (const Cells& s : solid) {
        if (s.klo != 0) { continue; }
        for (int i = s.ilo; i <= s.ihi; ++i) { for (int j = s.jlo; j <= s.jhi; ++j) { top[{i, j}] = s.khi + 1; } }
    }
    int n_roof = 0, n_second = 0;
    for (MFIter mfi(out); mfi.isValid(); ++mfi) {
        FArrayBox h(mfi.validbox(), 5, The_Pinned_Arena());
        Gpu::dtoh_memcpy(h.dataPtr(), out[mfi].dataPtr(), h.nBytes());
        const auto a = h.const_array();
        const Box& b = mfi.validbox();
        for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j) {
            for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i) {
                const int k_roof = top[{i, j}];
                // The upward face on that cell, if this rank holds it (it does: the face's cell is in this box).
                int found = -1, upward = 0;
                for (int f = 0; f < nf; ++f) {
                    if (fi[f] != i || fj[f] != j || fd[f] != 2 || fs[f] >= 0) { continue; }
                    ++upward;
                    if (fk[f] == k_roof) { found = f; }
                }
                if (k_roof > 0 && k_roof <= nz - 1) {
                    ASSERT_GE(found, 0) << "column " << i << "," << j;
                    EXPECT_EQ(a(i, j, 0, 0), Real(1.0)) << "column " << i << "," << j;
                    EXPECT_EQ(a(i, j, 0, 1), alb[found]) << "column " << i << "," << j;
                    EXPECT_EQ(a(i, j, 0, 2), eps[found]) << "column " << i << "," << j;
                    EXPECT_EQ(a(i, j, 0, 3), T[found]) << "column " << i << "," << j;
                    // Lit, except under the slab.
                    EXPECT_EQ(a(i, j, 0, 4), (upward > 1) ? Real(0.0) : Real(1.0)) << "column " << i << "," << j;
                    ++n_roof;
                    if (upward > 1) { ++n_second; }
                } else {
                    for (int comp = 0; comp < 5; ++comp) { EXPECT_EQ(a(i, j, 0, comp), Real(0.0)) << "column " << i << "," << j; }
                }
            }
        }
    }
    ParallelDescriptor::ReduceIntSum(n_roof);
    ParallelDescriptor::ReduceIntSum(n_second);
    // 8 columns of the block, 6 of the tower, 2 of the low building; the two under the
    // overhang hold a second upward face (the slab's), which is not taken. The slab over
    // open ground and the building that fills the column give no roof.
    EXPECT_EQ(n_roof, 16);
    EXPECT_EQ(n_second, 2);

    // A later step: the skin temperatures have moved and are written again; what does not
    // change is left as it is (here the albedos are changed on purpose and must not show).
    std::vector<Real> T2(T), alb2(alb);
    for (auto& v : T2) { v += Real(3.5); }
    for (auto& v : alb2) { v *= Real(0.5); }
    to_device(faces.d_T_skin, T2); to_device(faces.d_albedo, alb2);
    faces.fill_roofs(out, roof_cell, false, Real(0.0));
    int n_later = 0;
    for (MFIter mfi(out); mfi.isValid(); ++mfi) {
        FArrayBox h(mfi.validbox(), 5, The_Pinned_Arena());
        Gpu::dtoh_memcpy(h.dataPtr(), out[mfi].dataPtr(), h.nBytes());
        const auto a = h.const_array();
        for (int f = 0; f < nf; ++f) {
            if (!mfi.validbox().contains(IntVect(fi[f], fj[f], 0)) || fd[f] != 2 || fs[f] >= 0 || fk[f] != top[{fi[f], fj[f]}]) { continue; }
            EXPECT_EQ(a(fi[f], fj[f], 0, 3), T2[f]) << "face " << f;
            EXPECT_EQ(a(fi[f], fj[f], 0, 1), alb[f]) << "face " << f;
            EXPECT_EQ(a(fi[f], fj[f], 0, 0), Real(1.0)) << "face " << f;
            ++n_later;
        }
    }
    ParallelDescriptor::ReduceIntSum(n_later);
    EXPECT_EQ(n_later, 16);
}

// The roof's shadow is that of the time given, whatever sun the faces' last shortwave
// left: a restarted run holds the start-up sun there and a continuous run the previous
// step's, and the two must sweep the same roofs. A low building west of a tall tower on
// the equator at the equinox: the tower shades the low roof in the morning, not at noon.
TEST(IBSEBRoofs, TheShadowIsThatOfTheTimeGiven)
{
    using namespace amrex;
    const int nx = 16, ny = 16, nz = 8;
    const RealBox rb({0.0, 0.0, 0.0}, {160.0, 160.0, 80.0});
    const Box dom(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, nz - 1));
    const Geometry geom(dom, rb, 0, {1, 1, 0});
    BoxArray ba(dom);
    ba.maxSize(IntVect(8, 8, 8));
    const DistributionMapping dm(ba);
    MultiFab blank(ba, dm, 1, 1);
    // The low building (two cells) at i = 5..6 and the tower (six cells) at i = 7, east of it.
    const std::vector<Cells> solid = {{5, 6, 8, 8, 0, 1}, {7, 7, 7, 9, 0, 5}};
    fill_blanking(blank, geom, solid);
    BoxList bl = ba.boxList();
    for (Box& b : bl) { b.setRange(2, 0); }
    const BoxArray ba2d(std::move(bl));
    iMultiFab roof_cell(ba2d, dm, 1, 0);
    fill_roof_cell(roof_cell, blank, nz - 1);

    IBSEBParams params;
    params.enable = true;
    params.view_n_az = 8;
    params.view_n_el = 4;
    params.sun_mode = "solar";
    params.latitude_deg = 0.0;
    params.longitude_deg = 0.0;
    params.day_of_year = 80.0;
    params.time_zero_utc_s = 0.0;
    IBFaceSet faces(params, 0);
    faces.build(blank, geom);
    faces.release_labels();
    faces.compute_view_fractions();
    const Real morning = Real(9.0 * 3600.0), noon = Real(12.0 * 3600.0);

    // The lit share of the low roof's two columns and of the tower's three.
    auto lit = [&] (Real t_last_shortwave, Real t_roofs, int& n_low, int& n_tower) {
        faces.compute_shortwave(t_last_shortwave);
        MultiFab out(ba2d, dm, 5, 0);
        faces.fill_roofs(out, roof_cell, true, t_roofs);
        Real low = 0.0, tower = 0.0;
        n_low = 0; n_tower = 0;
        for (MFIter mfi(out); mfi.isValid(); ++mfi) {
            FArrayBox h(mfi.validbox(), 5, The_Pinned_Arena());
            Gpu::dtoh_memcpy(h.dataPtr(), out[mfi].dataPtr(), h.nBytes());
            const auto a = h.const_array();
            const Box& b = mfi.validbox();
            for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j) {
                for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i) {
                    if (a(i, j, 0, 0) < Real(0.5)) { continue; }
                    if (i == 7) { tower += a(i, j, 0, 4); ++n_tower; } else { low += a(i, j, 0, 4); ++n_low; }
                }
            }
        }
        ParallelDescriptor::ReduceRealSum(low);  ParallelDescriptor::ReduceRealSum(tower);
        ParallelDescriptor::ReduceIntSum(n_low); ParallelDescriptor::ReduceIntSum(n_tower);
        return std::make_pair(low, tower);
    };
    int n_low = 0, n_tower = 0;
    // The faces' last shortwave was the other time's in each pair: it must not show.
    const auto am = lit(noon, morning, n_low, n_tower);
    EXPECT_EQ(n_low, 2);
    EXPECT_EQ(n_tower, 3);
    EXPECT_EQ(am.first, Real(0.0));    // the tower shades both columns of the low roof
    EXPECT_EQ(am.second, Real(3.0));   // nothing shades the tower's roof
    const auto pm = lit(morning, noon, n_low, n_tower);
    EXPECT_EQ(pm.first, Real(2.0));    // the sun overhead: the low roof is lit
    EXPECT_EQ(pm.second, Real(3.0));
    // The same time after another last shortwave: the same roofs.
    const auto am2 = lit(morning, morning, n_low, n_tower);
    EXPECT_EQ(am2.first, am.first);
}

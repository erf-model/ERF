#include <cmath>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_IBFaceSet.H"
#include "ERF_IBSEBParams.H"

// The buildings' shadow on the ground for the two-stream ground balance
// (erf.ibseb.two_stream_buildings): IBFaceSet::compute_ground_shadow() casts the ray from
// the centre of every column's ground toward the sun of the last compute_shortwave()
// through the faces' column map. Here the expected shadow comes from intersecting that
// ray with every solid column as a box (the slab method), the periodic images included,
// not from the code's column walk. A tower against the high-x edge throws its shadow
// across the periodic boundary; the same tower on a closed domain does not; an overhang
// shades the ground under it; a sun below the horizon shades every column. The level is
// split into four boxes, so a box that read another box's columns would show.

namespace {

struct Cells { int ilo, ihi, jlo, jhi, khi; int klo = 0; };

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

// Expected shadow of column (i, j): the ray from its ground centre along s meets the box
// [x, x + dx] x [y, y + dy] x [0, top] of a solid column (the column map's model of a
// building: solid from the ground to its top), or of one of its periodic images, below
// the top on entry. A column under a building is in shadow.
int expected_shadow (int i, int j, const std::vector<Cells>& solid, int nx, int ny, double d, double dz,
                     double sx, double sy, double sz, bool periodic)
{
    if (sz <= 0.0) { return 1; }
    const double xc = (i + 0.5) * d, yc = (j + 0.5) * d;
    const double big = 1.0e300;
    for (const Cells& s : solid) {
        const double top = (s.khi + 1) * dz;
        for (int ci = s.ilo; ci <= s.ihi; ++ci) {
            for (int cj = s.jlo; cj <= s.jhi; ++cj) {
                if (ci == i && cj == j) { return 1; }
                for (int ox = (periodic ? -1 : 0); ox <= (periodic ? 1 : 0); ++ox) {
                    for (int oy = (periodic ? -1 : 0); oy <= (periodic ? 1 : 0); ++oy) {
                        const double x0 = ci * d + ox * nx * d, x1 = x0 + d;
                        const double y0 = cj * d + oy * ny * d, y1 = y0 + d;
                        double tx0 = -big, tx1 = big, ty0 = -big, ty1 = big;
                        if (sx != 0.0) { tx0 = std::min((x0 - xc) / sx, (x1 - xc) / sx); tx1 = std::max((x0 - xc) / sx, (x1 - xc) / sx); }
                        else if (!(xc > x0 && xc < x1)) { continue; }
                        if (sy != 0.0) { ty0 = std::min((y0 - yc) / sy, (y1 - yc) / sy); ty1 = std::max((y0 - yc) / sy, (y1 - yc) / sy); }
                        else if (!(yc > y0 && yc < y1)) { continue; }
                        const double t_in = std::max(std::max(tx0, ty0), 0.0);
                        const double t_out = std::min(tx1, ty1);
                        if (t_out > t_in && t_in * sz < top) { return 1; }
                    }
                }
            }
        }
    }
    return 0;
}

// The shadow compute_ground_shadow() gives, as a host array [i * ny + j].
std::vector<double> ground_shadow (bool periodic, double zenith_deg, double azimuth_deg,
                                   const std::vector<Cells>& solid)
{
    using namespace amrex;
    // 8 x 8 columns of 20 m, 8 layers of 20 m, in four boxes.
    const RealBox rb({0.0, 0.0, 0.0}, {160.0, 160.0, 160.0});
    const Box dom(IntVect(0, 0, 0), IntVect(7, 7, 7));
    const Geometry geom(dom, rb, 0, {int(periodic), int(periodic), 0});
    BoxArray ba(dom);
    ba.maxSize(IntVect(4, 4, 8));
    const DistributionMapping dm(ba);
    MultiFab blank(ba, dm, 1, 1);
    fill_blanking(blank, geom, solid);

    IBSEBParams params;
    params.enable = true;
    params.sun_mode = "fixed";
    params.sun_zenith_deg = zenith_deg;
    params.sun_azimuth_deg = azimuth_deg;
    params.view_n_az = 4;
    params.view_n_el = 2;
    IBFaceSet faces(params, 0);
    faces.build(blank, geom);
    faces.release_labels();
    faces.compute_view_fractions();
    faces.compute_shortwave(Real(0.0));

    BoxList bl = ba.boxList();
    for (Box& b : bl) { b.setRange(2, 0); }
    MultiFab shadow(BoxArray(std::move(bl)), dm, 1, 0);
    shadow.setVal(-1.0);
    faces.compute_ground_shadow(shadow);

    std::vector<double> out(64, -1.0);
    for (MFIter mfi(shadow); mfi.isValid(); ++mfi) {
        FArrayBox host(mfi.validbox(), 1, The_Pinned_Arena());
        host.copy<RunOn::Device>(shadow[mfi], mfi.validbox(), 0, mfi.validbox(), 0, 1);
        Gpu::streamSynchronize();
        const auto h = host.const_array();
        const Box& bx = mfi.validbox();
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) { out[i * 8 + j] = h(i, j, 0); }
        }
    }
    ParallelDescriptor::ReduceRealMax(out.data(), static_cast<int>(out.size()));
    return out;
}

void expect_geometry (bool periodic, double zenith_deg, double azimuth_deg, const std::vector<Cells>& solid,
                      int& n_shaded)
{
    const auto shadow = ground_shadow(periodic, zenith_deg, azimuth_deg, solid);
    const double pi = amrex::Math::pi<double>();
    const double z = zenith_deg * pi / 180.0, a = azimuth_deg * pi / 180.0;
    const double sx = std::sin(a) * std::sin(z), sy = std::cos(a) * std::sin(z), sz = std::cos(z);
    n_shaded = 0;
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            const int e = expected_shadow(i, j, solid, 8, 8, 20.0, 20.0, sx, sy, sz, periodic);
            EXPECT_EQ(shadow[i * 8 + j], static_cast<double>(e))
                << "column (" << i << "," << j << "), periodic " << periodic << ", zenith " << zenith_deg;
            n_shaded += e;
        }
    }
}

} // namespace

TEST(IBSEBGroundShadow, TowerShadowCrossesThePeriodicEdge)
{
    // A 60 m tower on columns i = 7, j = 3..4 against the high-x edge, the sun 60 degrees
    // from the zenith in the west-south-west: its shadow, about 104 m long, falls east,
    // across the edge onto the columns at low i.
    const std::vector<Cells> tower{{7, 7, 3, 4, 2}};
    int n_periodic = 0, n_closed = 0;
    expect_geometry(true, 60.0, 250.0, tower, n_periodic);
    expect_geometry(false, 60.0, 250.0, tower, n_closed);
    // Not vacuous: the shadow reaches across the edge, and only when the domain wraps.
    EXPECT_GE(n_periodic, n_closed + 3);
    EXPECT_EQ(n_closed, 2);   // only the tower's own columns on the closed domain
}

TEST(IBSEBGroundShadow, OverhangAndBlockInFourBoxes)
{
    // A 40 m block over columns spanning the four boxes and an overhang (solid from 60 m
    // to 100 m over open ground): the ground under the overhang is shaded, as the column
    // map takes a column as solid up to its top.
    const std::vector<Cells> blocks{{3, 4, 3, 4, 1}, {1, 1, 6, 6, 4, 3}};
    int n = 0;
    expect_geometry(true, 35.0, 140.0, blocks, n);
    EXPECT_GE(n, 8);
}

TEST(IBSEBGroundShadow, SunBelowTheHorizonShadesEveryColumn)
{
    const auto shadow = ground_shadow(true, 100.0, 90.0, {{3, 4, 3, 4, 1}});
    for (double v : shadow) { EXPECT_EQ(v, 1.0); }
}

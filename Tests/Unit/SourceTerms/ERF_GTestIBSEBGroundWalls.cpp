#include <cmath>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_IBFaceSet.H"
#include "ERF_IBSEBBalance.H"
#include "ERF_IBSEBParams.H"
#include "ERF_IBSEBSolar.H"
#include "ERF_TwoStreamCanopyForcing.H"

// What the walls send to the ground beside them (erf.ibseb.ground_sky_view,
// IBFaceSet::compute_ground_walls()): each face adds, to the columns its ground rays land
// in, its area times the share of its rays that land there, over the column's area, and
// that weight times its reflected shortwave and its longwave.
//  - A long wall: the walls' share of the ground's view against the closed form of the
//    view factor from a ground strip to the wall.
//  - A set of buildings in four boxes on a periodic domain, every face with its own
//    radiosity: the field column by column against a sum made here from the faces'
//    landings, the totals against what the faces send, and the same field from one box.
//  - No buildings, and the terms of the ground balance.

namespace {

struct Cells { int ilo, ihi, jlo, jhi, khi; };

void fill_blanking (amrex::MultiFab& b, const amrex::Geometry& geom, const std::vector<Cells>& solid)
{
    b.setVal(0.0);
    for (amrex::MFIter mfi(b); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& a = b.array(mfi);
        for (const Cells& s : solid) {
            const amrex::Box ov = amrex::Box(amrex::IntVect(s.ilo, s.jlo, 0), amrex::IntVect(s.ihi, s.jhi, s.khi)) & bx;
            if (ov.isEmpty()) { continue; }
            amrex::ParallelFor(ov, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { a(i, j, k) = 1.0; });
        }
    }
    b.FillBoundary(geom.periodicity());
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

amrex::BoxArray columns_of (const amrex::BoxArray& ba)
{
    amrex::BoxList bl = ba.boxList();
    for (amrex::Box& b : bl) { b.setRange(2, 0); }
    return amrex::BoxArray(std::move(bl));
}

// The field as a host table over the nx x ny columns of the domain.
std::vector<double> table (const amrex::MultiFab& mf, int comp, int nx, int ny)
{
    std::vector<double> t(static_cast<size_t>(nx) * ny, 0.0);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        amrex::FArrayBox hfab(mfi.validbox(), mf.nComp(), amrex::The_Pinned_Arena());
        amrex::Gpu::dtoh_memcpy(hfab.dataPtr(), mf[mfi].dataPtr(), hfab.nBytes());
        const auto a = hfab.const_array();
        const amrex::Box& b = mfi.validbox();
        for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j) {
            for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i) { t[static_cast<size_t>(i) + static_cast<size_t>(nx) * j] = a(i, j, 0, comp); }
        }
    }
    amrex::ParallelDescriptor::ReduceRealSum(t.data(), static_cast<int>(t.size()));
    return t;
}

// Give every face its own albedo, emissivity, skin temperature and incoming radiation.
void set_radiosities (IBFaceSet& faces)
{
    const int nf = faces.n_faces();
    std::vector<amrex::Real> alb(nf), eps(nf), T(nf), dir(nf), dif(nf), lw(nf);
    const auto i = host(faces.d_i), j = host(faces.d_j), k = host(faces.d_k), d = host(faces.d_dir);
    for (int f = 0; f < nf; ++f) {
        const amrex::Real u = amrex::Real(0.37) * i[f] + amrex::Real(0.61) * j[f] + amrex::Real(1.3) * k[f] + d[f];
        alb[f] = amrex::Real(0.2) + amrex::Real(0.01) * std::fmod(u, amrex::Real(30.0));
        eps[f] = amrex::Real(0.7) + amrex::Real(0.01) * std::fmod(amrex::Real(3.0) * u, amrex::Real(25.0));
        T[f]   = amrex::Real(280.0) + amrex::Real(2.0) * std::fmod(u, amrex::Real(17.0));
        dir[f] = amrex::Real(40.0) * std::fmod(u, amrex::Real(11.0));
        dif[f] = amrex::Real(60.0) + amrex::Real(5.0) * std::fmod(u, amrex::Real(7.0));
        lw[f]  = amrex::Real(330.0) + amrex::Real(3.0) * std::fmod(u, amrex::Real(13.0));
    }
    to_device(faces.d_albedo, alb); to_device(faces.d_emis, eps); to_device(faces.d_T_skin, T);
    to_device(faces.d_SW_direct_in, dir); to_device(faces.d_SW_diffuse_in, dif); to_device(faces.d_LW_down_in, lw);
}

} // namespace

// A wall H = 40 m high across the whole (periodic) y extent of a domain closed in x. By
// the view factor from an element to a surface that is infinite along one direction,
// (sin phi_2 - sin phi_1) / 2 with the angles from the element's normal, a point of the
// ground at the distance x from the wall's foot sees the wall under
// (1 - x / sqrt(x^2 + H^2)) / 2, and the mean over a column from x0 to x1 is
// (1 - (sqrt(x1^2 + H^2) - sqrt(x0^2 + H^2)) / (x1 - x0)) / 2: 0.438, 0.325, 0.236, 0.172
// for the first four 10 m columns. The code gets the share from the wall's own rays (four
// faces in the height, 32 x 16 rays each), by reciprocity; the closed form was checked
// against a numerical integration of the view factor's integral.
TEST(IBSEBGroundWalls, LongWallClosedForm)
{
    using namespace amrex;
    const int nx = 48, ny = 8;
    const RealBox rb({0.0, 0.0, 0.0}, {480.0, 80.0, 80.0});
    const Box dom(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, 7));
    const Geometry geom(dom, rb, 0, {0, 1, 0});
    BoxArray ba(dom);
    ba.maxSize(IntVect(16, 8, 8));
    const DistributionMapping dm(ba);
    MultiFab blank(ba, dm, 1, 1);
    fill_blanking(blank, geom, {{22, 23, 0, ny - 1, 3}});
    IBSEBParams params;
    params.enable = true;
    params.view_n_az = 32;
    params.view_n_el = 16;
    params.ground_sky_view = true;
    IBFaceSet faces(params, 0);
    faces.build(blank, geom);
    faces.release_labels();
    faces.compute_view_fractions();
    MultiFab out(columns_of(ba), dm, IBFaceSet::GW_NCOMP, 0);
    faces.compute_ground_walls(out, true);
    const auto view = table(out, IBFaceSet::GW_VIEW, nx, ny);
    const double H = 40.0;
    for (int n = 0; n < 8; ++n) {
        const double x0 = 10.0 * n, x1 = 10.0 * (n + 1);
        const double expect = 0.5 * (1.0 - (std::sqrt(x1 * x1 + H * H) - std::sqrt(x0 * x0 + H * H)) / (x1 - x0));
        for (int j = 0; j < ny; ++j) {
            // Every face sends from its centre: the column at the wall's foot, which the
            // lowest face sees from 5 m up, comes out 0.02 low (0.418 for 0.438), and the
            // columns behind it up to 0.014 high (0.250 for 0.236) or 0.010 low farther out.
            // The sum over all the columns is that of the faces' ground views.
            const double tol = (n == 0) ? 0.025 : ((n < 4) ? 0.016 : 0.012);
            EXPECT_NEAR(view[static_cast<size_t>(24 + n) + static_cast<size_t>(nx) * j], expect, tol) << "east column " << n << " row " << j;
            EXPECT_NEAR(view[static_cast<size_t>(21 - n) + static_cast<size_t>(nx) * j], expect, tol) << "west column " << n << " row " << j;
            if (n == 0) { EXPECT_LT(view[static_cast<size_t>(24) + static_cast<size_t>(nx) * j], expect); }
        }
    }
    // Under the wall the ground sees nothing of it: no ray lands in a solid column.
    for (int j = 0; j < ny; ++j) {
        EXPECT_EQ(view[static_cast<size_t>(22) + static_cast<size_t>(nx) * j], 0.0);
        EXPECT_EQ(view[static_cast<size_t>(23) + static_cast<size_t>(nx) * j], 0.0);
    }
}

// Buildings in four boxes on a periodic domain, one of them across the periodic edge.
TEST(IBSEBGroundWalls, WhatTheFacesSendIsWhatTheGroundGets)
{
    using namespace amrex;
    const int nx = 16, ny = 16;
    const RealBox rb({0.0, 0.0, 0.0}, {160.0, 160.0, 80.0});
    const Box dom(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, 7));
    const Geometry geom(dom, rb, 0, {1, 1, 0});
    const std::vector<Cells> solid = {{6, 9, 7, 8, 3}, {0, 0, 2, 4, 5}, {15, 15, 2, 4, 5}, {12, 13, 12, 12, 1}};
    const Real sigma = ibseb::SIGMA_SB;
    const Real area_col = 100.0;
    std::vector<std::vector<double>> fields;   // the three components, four boxes then one
    for (const int max_size : {8, 16}) {
        BoxArray ba(dom);
        ba.maxSize(IntVect(max_size, max_size, 8));
        const DistributionMapping dm(ba);
        MultiFab blank(ba, dm, 1, 1);
        fill_blanking(blank, geom, solid);
        IBSEBParams params;
        params.enable = true;
        params.view_n_az = 12;
        params.view_n_el = 6;
        params.ground_sky_view = true;
        IBFaceSet faces(params, 0);
        faces.build(blank, geom);
        faces.release_labels();
        faces.compute_view_fractions();
        set_radiosities(faces);
        MultiFab out(columns_of(ba), dm, IBFaceSet::GW_NCOMP, 0);
        out.setVal(-7.0);   // every column is set, not added to what was there
        faces.compute_ground_walls(out, true);

        // The same sums made here from each face's landings, column by column.
        const auto fi = host(faces.d_i), fj = host(faces.d_j);
        const auto start = host(faces.d_gnd_start), rays = host(faces.d_gnd_rays), far_rays = host(faces.d_gnd_far);
        const auto land = host(faces.d_gnd_land);
        const auto area = host(faces.d_area), alb = host(faces.d_albedo), eps = host(faces.d_emis), T = host(faces.d_T_skin);
        const auto dir = host(faces.d_SW_direct_in), dif = host(faces.d_SW_diffuse_in), lw = host(faces.d_LW_down_in);
        const auto fg = host(faces.d_f_ground);
        const double inv = 1.0 / 72.0;
        std::vector<double> expect[3];
        for (auto& e : expect) { e.assign(static_cast<size_t>(nx) * ny, 0.0); }
        double sent[3] = {0.0, 0.0, 0.0};
        int n_wall = 0, n_far = 0;
        for (int f = 0; f < faces.n_faces(); ++f) {
            const double j_sw = alb[f] * (dir[f] + dif[f]);
            const double j_lw = eps[f] * sigma * std::pow(double(T[f]), 4) + (1.0 - eps[f]) * lw[f];
            int n_kept = 0;
            for (int m = start[f]; m < start[f + 1]; ++m) {
                int di = 0, dj = 0;
                ibseb::unpack_landing(land[m], di, dj);
                const int ci = ((fi[f] + di) % nx + nx) % nx, cj = ((fj[f] + dj) % ny + ny) % ny;
                const double w = area[f] * rays[m] * inv / area_col;
                const size_t c = static_cast<size_t>(ci) + static_cast<size_t>(nx) * cj;
                expect[0][c] += w; expect[1][c] += w * j_sw; expect[2][c] += w * j_lw;
                n_kept += rays[m];
            }
            // All a face sends to the ground it keeps: its area times its kept ground view.
            const double share = area[f] * n_kept * inv;
            sent[0] += share; sent[1] += share * j_sw; sent[2] += share * j_lw;
            if (n_kept > 0) { ++n_wall; }
            n_far += far_rays[f];
            EXPECT_NEAR(fg[f], (n_kept + far_rays[f]) * inv, 1.0e-6) << "face " << f;
        }
        ParallelDescriptor::ReduceRealSum(sent, 3);
        ParallelDescriptor::ReduceIntSum(n_wall);
        ParallelDescriptor::ReduceIntSum(n_far);
        for (auto& e : expect) { ParallelDescriptor::ReduceRealSum(e.data(), static_cast<int>(e.size())); }
        const double tol = (sizeof(Real) == 8) ? 1.0e-12 : 2.0e-5;
        int n_seen = 0;
        for (int comp = 0; comp < 3; ++comp) {
            const auto got = table(out, comp, nx, ny);
            double got_sum = 0.0;
            for (size_t c = 0; c < got.size(); ++c) {
                EXPECT_NEAR(got[c], expect[comp][c], tol * (std::abs(expect[comp][c]) + 1.0)) << "component " << comp << " column " << c;
                got_sum += got[c] * area_col;
                if (comp == 0 && got[c] > 0.0) { ++n_seen; }
            }
            EXPECT_NEAR(got_sum, sent[comp], tol * (std::abs(sent[comp]) + 1.0)) << "component " << comp;
            EXPECT_GT(sent[comp], 0.0);
            fields.push_back(got);
        }
        // A later step: the radiosities change, the share of the view is left as it is.
        {
            auto alb = host(faces.d_albedo);
            for (auto& v : alb) { v *= Real(0.5); }
            to_device(faces.d_albedo, alb);
            out.plus(Real(1.0), IBFaceSet::GW_VIEW, 1, 0);
            faces.compute_ground_walls(out, false);
            const auto view = table(out, IBFaceSet::GW_VIEW, nx, ny), sw = table(out, IBFaceSet::GW_SW, nx, ny);
            const auto lwv = table(out, IBFaceSet::GW_LW, nx, ny);
            for (size_t c = 0; c < view.size(); ++c) {
                EXPECT_NEAR(view[c], fields[fields.size() - 3][c] + 1.0, tol * 2.0) << "column " << c;   // untouched: still what this test put there
                EXPECT_NEAR(sw[c], 0.5 * fields[fields.size() - 2][c], tol * (std::abs(fields[fields.size() - 2][c]) + 1.0)) << "column " << c;
                EXPECT_NEAR(lwv[c], fields[fields.size() - 1][c], tol * (std::abs(fields[fields.size() - 1][c]) + 1.0)) << "column " << c;
            }
        }
        // Not vacuous: many walls, many columns beside them, and a domain small enough for
        // every ground ray to be kept (the map spans it).
        EXPECT_GT(n_wall, 40);
        EXPECT_GT(n_seen, 100);
        EXPECT_EQ(n_far, 0);
    }
    // The same field however the level is cut into boxes.
    const double tol = (sizeof(amrex::Real) == 8) ? 1.0e-12 : 2.0e-5;
    for (int comp = 0; comp < 3; ++comp) {
        for (size_t c = 0; c < fields[comp].size(); ++c) {
            EXPECT_NEAR(fields[comp][c], fields[3 + comp][c], tol * (std::abs(fields[comp][c]) + 1.0)) << "component " << comp << " column " << c;
        }
    }
}

// A level without buildings sends nothing; the field is still set, to zero.
TEST(IBSEBGroundWalls, NoBuildingsNothingSent)
{
    using namespace amrex;
    const RealBox rb({0.0, 0.0, 0.0}, {160.0, 160.0, 80.0});
    const Box dom(IntVect(0, 0, 0), IntVect(15, 15, 7));
    const Geometry geom(dom, rb, 0, {1, 1, 0});
    BoxArray ba(dom);
    ba.maxSize(IntVect(8, 8, 8));
    const DistributionMapping dm(ba);
    MultiFab blank(ba, dm, 1, 1);
    fill_blanking(blank, geom, {});
    IBSEBParams params;
    params.enable = true;
    params.ground_sky_view = true;
    IBFaceSet faces(params, 0);
    faces.build(blank, geom);
    faces.release_labels();
    faces.compute_view_fractions();
    MultiFab out(columns_of(ba), dm, IBFaceSet::GW_NCOMP, 0);
    out.setVal(3.0);
    faces.compute_ground_walls(out, true);
    for (int comp = 0; comp < 3; ++comp) {
        EXPECT_EQ(out.min(comp), Real(0.0));
        EXPECT_EQ(out.max(comp), Real(0.0));
    }
}

// The ground balance beside walls: through the walls' share of its view the ground loses
// the sky's diffuse light and longwave and gains the walls'.
TEST(IBSEBGroundWalls, GroundBalanceTerms)
{
    using amrex::Real;
    const Real tol = (sizeof(Real) == 8) ? Real(1.e-12) : Real(1.e-4);
    // 500 W/m2 absorbed, 380 of it the beam: 120 of diffuse light. The walls take 0.3 of
    // the view, send 50 W/m2 of reflected light (the ground absorbs 0.8) and 140 W/m2 of
    // longwave, against a sky longwave of 350 W/m2; the ground's emissivity is 0.95.
    const GroundWallTerms t = two_stream_ground_wall_terms(Real(500.0), Real(380.0), Real(0.3), Real(50.0), Real(140.0),
                                                           Real(350.0), Real(0.8), Real(0.95));
    EXPECT_NEAR(t.d_sw, Real(0.8 * 50.0 - 0.3 * 120.0), tol * 100);
    EXPECT_NEAR(t.d_lw, Real(0.95 * (140.0 - 0.3 * 350.0)), tol * 100);
    // No walls in view: nothing changes.
    const GroundWallTerms none = two_stream_ground_wall_terms(Real(500.0), Real(380.0), Real(0.0), Real(0.0), Real(0.0),
                                                              Real(350.0), Real(0.8), Real(0.95));
    EXPECT_EQ(none.d_sw, Real(0.0));
    EXPECT_EQ(none.d_lw, Real(0.0));
    // Walls at the sky's brightness change nothing: 0.3 of a sky of 350 W/m2 is 105, and
    // 0.3 of the 150 W/m2 of diffuse light the 120 absorbed came from is 45.
    const GroundWallTerms same = two_stream_ground_wall_terms(Real(500.0), Real(380.0), Real(0.3), Real(45.0), Real(105.0),
                                                              Real(350.0), Real(0.8), Real(0.95));
    EXPECT_NEAR(same.d_sw, Real(0.0), tol * 100);
    EXPECT_NEAR(same.d_lw, Real(0.0), tol * 100);
    // A count of rays above the whole view (1.25): the loss and the gain are scaled alike,
    // so the ground sees no sky and the walls' mean.
    const GroundWallTerms full = two_stream_ground_wall_terms(Real(500.0), Real(380.0), Real(1.25), Real(50.0), Real(140.0),
                                                              Real(350.0), Real(0.8), Real(0.95));
    EXPECT_NEAR(full.d_sw, Real(0.8 * 40.0 - 120.0), tol * 100);
    EXPECT_NEAR(full.d_lw, Real(0.95 * (112.0 - 350.0)), tol * 100);
}

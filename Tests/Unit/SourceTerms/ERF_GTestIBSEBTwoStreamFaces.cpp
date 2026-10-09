#include <algorithm>
#include <cmath>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_IBFaceSet.H"
#include "ERF_IBSEBParams.H"
#include "ERF_IndexDefines.H"
#include "ERF_TwoStreamCanopyForcing.H"

// With erf.ibseb.radiation = two_stream every face reads the column of its own fluid cell at
// its own height (a roof its interface, a wall the mean of its cell's two): the beam over the
// cosine of the zenith, the shortwave down less the beam, and the longwave down and the
// shortwave and longwave up there. Here the beam and the interface fluxes differ in every column and at
// every interface, and the level is split into four boxes, so a face that read another
// column, another height, another box's fab or the ground for its sky would get another
// number. The blanking is set by hand (1 solid, 0 fluid), no embedded boundary.

namespace {

constexpr amrex::Real kCosZ = 0.766044443118978;   // 40 degrees

// The beam and the fluxes of column (i, j) at interface m (fluxes at absolute k = m).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real beam_at (int i, int j, int m) { return amrex::Real(200.0 + 7.0 * i + 3.0 * j + 11.0 * m); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real swdn_at (int i, int j, int m) { return beam_at(i, j, m) + amrex::Real(20.0 + i + 2.0 * j + 3.0 * m); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real lwdn_at (int i, int j, int m) { return amrex::Real(300.0 + 2.0 * i + 5.0 * j - 4.0 * m); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real swup_at (int i, int j, int m) { return amrex::Real(50.0 + i + 0.5 * j - 2.0 * m); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real lwup_at (int i, int j, int m) { return amrex::Real(400.0 + 1.5 * i - j + 6.0 * m); }

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

void fill_view (amrex::MultiFab& beam, amrex::MultiFab& cosz, amrex::MultiFab& fluxes)
{
    for (amrex::MFIter mfi(beam); mfi.isValid(); ++mfi) {
        auto const& bA = beam.array(mfi);
        auto const& cA = cosz.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int m) noexcept { bA(i, j, m) = beam_at(i, j, m); });
        amrex::ParallelFor(cosz[mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept { cA(i, j, 0) = kCosZ; });
    }
    for (amrex::MFIter mfi(fluxes); mfi.isValid(); ++mfi) {
        auto const& fA = fluxes.array(mfi);
        amrex::ParallelFor(fluxes[mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            fA(i, j, k, 0) = swup_at(i, j, k);
            fA(i, j, k, 1) = swdn_at(i, j, k);
            fA(i, j, k, 2) = lwup_at(i, j, k);
            fA(i, j, k, 3) = lwdn_at(i, j, k);
        });
    }
}

template <class T>
std::vector<T> host (const amrex::Gpu::DeviceVector<T>& d)
{
    std::vector<T> h(d.size());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d.begin(), d.end(), h.begin());
    amrex::Gpu::streamSynchronize();
    return h;
}

} // namespace

TEST(IBSEBTwoStreamFaces, EachFaceReadsItsOwnColumnAtItsHeight)
{
    using namespace amrex;
    // 8 x 8 columns of 20 m, 8 layers of 10 m, periodic in x and y, in four boxes.
    const RealBox rb({0.0, 0.0, 0.0}, {160.0, 160.0, 80.0});
    const Box dom(IntVect(0, 0, 0), IntVect(7, 7, 7));
    const Geometry geom(dom, rb, 0, {1, 1, 0});
    BoxArray ba(dom);
    ba.maxSize(IntVect(4, 4, 8));
    const DistributionMapping dm(ba);

    MultiFab blank(ba, dm, 1, 1);
    // A 30 m block over 2 x 3 columns and a 50 m tower on one column.
    fill_blanking(blank, geom, {{2, 3, 3, 5, 2}, {6, 6, 1, 1, 4}});

    IBSEBParams params;
    params.enable = true;
    params.radiation = "two_stream";
    params.lw_mode = "two_stream";
    params.sun_mode = "fixed";
    params.sun_zenith_deg = 40.0;
    params.sun_azimuth_deg = 200.0;
    params.view_n_az = 8;
    params.view_n_el = 4;
    IBFaceSet faces(params, 0);
    faces.build(blank, geom);
    faces.release_labels();
    faces.compute_view_fractions();
    const int m_top = faces.top_sample_interface();
    ASSERT_EQ(m_top, 5);

    // The view as the sweep and ERF lay it out: the beam on the columns with interfaces
    // 0 .. m_top, the cosine on the columns, the fluxes on the level's grids with a ghost.
    BoxList bl = ba.boxList();
    for (Box& b : bl) { b.setRange(2, 0, m_top + 1); }
    MultiFab beam(BoxArray(std::move(bl)), dm, 1, 0);
    BoxList bl2 = ba.boxList();
    for (Box& b : bl2) { b.setRange(2, 0); }
    MultiFab cosz(BoxArray(std::move(bl2)), dm, 1, 0);
    MultiFab fluxes(ba, dm, 4, IntVect(1, 1, 1));
    fill_view(beam, cosz, fluxes);
    MultiFab cons(ba, dm, 2, 0);
    cons.setVal(1.0, Rho_comp, 1);
    cons.setVal(300.0, RhoTheta_comp, 1);

    TwoStreamCanopyView view;
    view.beam = &beam;
    view.cos_zenith = &cosz;
    view.fluxes = &fluxes;
    faces.check_canopy_layout(cons, view);
    faces.compute_shortwave(Real(0.0), view);
    faces.compute_longwave(cons, view);

    const auto i = host(faces.d_i), j = host(faces.d_j), k = host(faces.d_k), dir = host(faces.d_dir), side = host(faces.d_side);
    const auto fs = host(faces.d_f_sky), fg = host(faces.d_f_ground), sh = host(faces.d_shadow);
    const auto dirin = host(faces.d_SW_direct_in), difin = host(faces.d_SW_diffuse_in), lwext = host(faces.d_LW_ext);
    const auto& sun = faces.sun();
    const Real tol = (sizeof(Real) == 8) ? Real(1.e-12) : Real(1.e-5);

    int n_roof = 0, n_wall = 0, n_lit = 0;
    std::vector<int> columns;
    for (int f = 0; f < faces.n_faces(); ++f) {
        const Real w_up = (dir[f] == 2) ? Real(0.0) : Real(0.5);
        const int m = k[f];
        const int mu = (dir[f] == 2) ? m : m + 1;
        const Real b  = (1.0 - w_up) * beam_at(i[f], j[f], m) + w_up * beam_at(i[f], j[f], mu);
        const Real sd = (1.0 - w_up) * swdn_at(i[f], j[f], m) + w_up * swdn_at(i[f], j[f], mu);
        const Real ld = (1.0 - w_up) * lwdn_at(i[f], j[f], m) + w_up * lwdn_at(i[f], j[f], mu);
        const Real su = (1.0 - w_up) * swup_at(i[f], j[f], m) + w_up * swup_at(i[f], j[f], mu);
        const Real lu = (1.0 - w_up) * lwup_at(i[f], j[f], m) + w_up * lwup_at(i[f], j[f], mu);
        Real n[3] = {0.0, 0.0, 0.0};
        n[dir[f]] = -static_cast<Real>(side[f]);
        const Real cosi = n[0] * sun.sx + n[1] * sun.sy + n[2] * sun.sz;
        const Real direct = (cosi > 0.0) ? (b / kCosZ) * cosi * (1.0 - sh[f]) : Real(0.0);
        const Real diffuse = fs[f] * (sd - b) + fg[f] * su;
        const Real lw = fs[f] * ld + fg[f] * lu;
        EXPECT_NEAR(dirin[f], direct, tol * (std::abs(direct) + 1.0)) << "face " << f << " (" << i[f] << "," << j[f] << "," << k[f] << ") dir " << dir[f];
        EXPECT_NEAR(difin[f], diffuse, tol * (std::abs(diffuse) + 1.0)) << "face " << f;
        EXPECT_NEAR(lwext[f], lw, tol * (std::abs(lw) + 1.0)) << "face " << f;
        (dir[f] == 2 ? n_roof : n_wall) += 1;
        if (direct > 1.0) { ++n_lit; }
        columns.push_back(i[f] * 8 + j[f]);
    }
    // Not vacuous: roofs at two heights and walls in many columns, sunlit faces.
    std::sort(columns.begin(), columns.end());
    columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
    EXPECT_GT(n_roof, 6);
    EXPECT_GT(n_wall, 10);
    EXPECT_GT(n_lit, 4);
    EXPECT_GT(static_cast<int>(columns.size()), 10);
}

// The two-stream beam the faces read is kept up to the highest interface a face samples
// (IBFaceSet::top_sample_interface()). That is a roof's own interface, or a wall's upper one
// where no roof lies above the wall: a building reaching the domain top has walls in the
// top cell, which read the top interface. A level with no buildings samples nothing.
TEST(IBSEBTwoStreamFaces, TopSampleInterface)
{
    using namespace amrex;
    const RealBox rb({0.0, 0.0, 0.0}, {160.0, 160.0, 80.0});
    const Box dom(IntVect(0, 0, 0), IntVect(7, 7, 7));
    const Geometry geom(dom, rb, 0, {1, 1, 0});
    BoxArray ba(dom);
    ba.maxSize(IntVect(4, 4, 8));
    const DistributionMapping dm(ba);
    IBSEBParams params;
    params.enable = true;
    params.radiation = "two_stream";
    params.lw_mode = "two_stream";

    // A tower through the whole depth, no roof: its walls in cell 7 read interface 8.
    MultiFab blank(ba, dm, 1, 1);
    fill_blanking(blank, geom, {{2, 3, 3, 4, 7}});
    IBFaceSet tower(params, 0);
    tower.build(blank, geom);
    EXPECT_TRUE(tower.has_faces());
    EXPECT_EQ(tower.top_sample_interface(), 8);

    // A 30 m block: its roof sits on interface 3, above its walls' upper interface.
    fill_blanking(blank, geom, {{2, 3, 3, 4, 2}});
    IBFaceSet block(params, 0);
    block.build(blank, geom);
    EXPECT_EQ(block.top_sample_interface(), 3);

    // No buildings.
    fill_blanking(blank, geom, {});
    IBFaceSet open(params, 0);
    open.build(blank, geom);
    EXPECT_FALSE(open.has_faces());
    EXPECT_EQ(open.top_sample_interface(), -1);
}

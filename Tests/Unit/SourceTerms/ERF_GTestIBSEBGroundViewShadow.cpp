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
#include "ERF_IBSEBSolar.H"

// The sunlit part of the ground a face sees (erf.ibseb.ground_view_shadow): of the
// hemisphere rays of a face that end on the ground, IBFaceSet::compute_shortwave() counts
// those whose landing column is not in the buildings' shadow, f_ground_lit, and the
// ground reflects the direct beam onto the face through that part only.
//  - A long wall with the sun behind it: the closed form of the view factor from a
//    point on a wall to the ground beyond the wall's shadow.
//  - A set of buildings in four boxes on a periodic domain: the count against one made
//    here from the geometry (every ray intersected with the ground plane and with every
//    solid column as a box, periodic images included; the shadow at the landing column
//    by the same slab method), not from the code's column walk.
//  - The reflected light on every face, with the input on and off.
//  - The window of the shadow walk: a building farther from the face than the map's reach
//    is not seen, wherever the patch ends.

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

// First positive distance along the ray from p at which it is inside the box
// [x0, x1] x [y0, y1] x [0, top], or a negative number when it never is (slab method).
double box_entry (const double p[3], const double s[3], double x0, double x1, double y0, double y1, double top)
{
    const double lo[3] = {x0, y0, 0.0}, hi[3] = {x1, y1, top};
    double t_in = 0.0, t_out = 1.0e300;
    for (int d = 0; d < 3; ++d) {
        if (s[d] != 0.0) {
            const double a = (lo[d] - p[d]) / s[d], b = (hi[d] - p[d]) / s[d];
            t_in = std::max(t_in, std::min(a, b));
            t_out = std::min(t_out, std::max(a, b));
        } else if (!(p[d] > lo[d] && p[d] < hi[d])) {
            return -1.0;
        }
    }
    // A face lies on its own building's side: leaving it at t = 0 is no hit.
    return (t_out > t_in + 1.0e-9) ? t_in : -1.0;
}

// Whether the ray from p along s meets a solid column (or a periodic image of one within
// n_img domains) before the distance t_end.
bool hits_building (const double p[3], const double s[3], double t_end, const std::vector<Cells>& solid,
                    int nx, int ny, double d, double dz, int n_img)
{
    for (const Cells& c : solid) {
        const double top = (c.khi + 1) * dz;
        for (int ox = -n_img; ox <= n_img; ++ox) {
            for (int oy = -n_img; oy <= n_img; ++oy) {
                const double t = box_entry(p, s, (c.ilo + ox * nx) * d, (c.ihi + 1 + ox * nx) * d,
                                           (c.jlo + oy * ny) * d, (c.jhi + 1 + oy * ny) * d, top);
                if (t >= 0.0 && t < t_end) { return true; }
            }
        }
    }
    return false;
}

// The expected sunlit part of the ground view of a face at (xf, yf, zf): of its
// hemisphere rays, those that reach the ground before a building and land in a column
// whose centre the sun reaches.
void expected_ground_view (int dir, int side, double xf, double yf, double zf, int n_az, int n_el,
                           const double sun[3], const std::vector<Cells>& solid, int nx, int ny, double d, double dz,
                           int n_img, double& f_ground, double& f_lit)
{
    int n_gnd = 0, n_lit = 0;
    for (int ie = 0; ie < n_el; ++ie) {
        for (int ia = 0; ia < n_az; ++ia) {
            amrex::Real rx, ry, rz;
            ibseb::hemisphere_direction(dir, -side, ia, ie, n_az, n_el, rx, ry, rz);
            if (rz >= 0.0) { continue; }
            const double p[3] = {xf, yf, zf}, s[3] = {double(rx), double(ry), double(rz)};
            const double t_g = zf / (-s[2]);
            if (hits_building(p, s, t_g, solid, nx, ny, d, dz, n_img)) { continue; }
            ++n_gnd;
            const double xl = xf + t_g * s[0], yl = yf + t_g * s[1];
            const double c[3] = {(std::floor(xl / d) + 0.5) * d, (std::floor(yl / d) + 0.5) * d, 0.0};
            if (!hits_building(c, sun, 1.0e300, solid, nx, ny, d, dz, n_img)) { ++n_lit; }
        }
    }
    f_ground = static_cast<double>(n_gnd) / (n_az * n_el);
    f_lit    = static_cast<double>(n_lit) / (n_az * n_el);
}

void sun_vector (double zenith_deg, double azimuth_deg, double s[3])
{
    const double pi = amrex::Math::pi<double>();
    const double z = zenith_deg * pi / 180.0, a = azimuth_deg * pi / 180.0;
    s[0] = std::sin(a) * std::sin(z); s[1] = std::cos(a) * std::sin(z); s[2] = std::cos(z);
}

} // namespace

// A wall 40 m high across the whole (periodic) y extent, the sun in the west at 45 degrees:
// east of the wall the ground is shaded for L = 40 m, four columns, to a column edge. A
// point at height h on the east side sees the ground beyond the shadow under the angle
// atan(h / L) below the horizontal, so by the view factor from an element to a surface
// that is infinite along one direction, (sin phi_2 - sin phi_1) / 2 with the angles from
// the element's normal, f_ground_lit = (h / sqrt(h^2 + L^2)) / 2, of f_ground = 1 / 2. The
// west side sees only sunlit ground. The domain is closed in x: beyond it the ground is
// open and sunlit, as the formula has it. A second sun, lower, would throw a shadow of
// 44 m: the shadow is taken at the centres of the columns, so the same four columns are
// shaded and the answer is that of L = 40 m; a shadow taken at the points where the rays
// land would differ by 0.018 at 35 m.
TEST(IBSEBGroundViewShadow, LongWallClosedForm)
{
    using namespace amrex;
    const double pi = amrex::Math::pi<double>();
    for (const double zenith_deg : {45.0, std::atan(44.0 / 40.0) * 180.0 / pi}) {
    const RealBox rb({0.0, 0.0, 0.0}, {640.0, 80.0, 80.0});
    const Box dom(IntVect(0, 0, 0), IntVect(63, 7, 7));
    const Geometry geom(dom, rb, 0, {0, 1, 0});
    BoxArray ba(dom);
    ba.maxSize(IntVect(16, 4, 8));
    const DistributionMapping dm(ba);
    MultiFab blank(ba, dm, 1, 1);
    fill_blanking(blank, geom, {{8, 9, 0, 7, 3}});

    IBSEBParams params;
    params.enable = true;
    params.sun_mode = "fixed";
    params.sun_zenith_deg = static_cast<Real>(zenith_deg);
    params.sun_azimuth_deg = 270.0;
    params.view_n_az = 256;
    params.view_n_el = 128;
    params.ground_view_shadow = true;
    IBFaceSet faces(params, 0);
    faces.build(blank, geom);
    faces.release_labels();
    faces.compute_view_fractions();
    faces.compute_shortwave(Real(0.0));

    const auto k = host(faces.d_k), dir = host(faces.d_dir), side = host(faces.d_side);
    const auto fg = host(faces.d_f_ground), fl = host(faces.d_f_ground_lit);
    const double L = 40.0;
    // The stratified sampling is within 2e-4 of these view factors (5e-4 for shadows of
    // 10 to 80 m); a shadow one column longer or shorter moves them by 0.012 or more.
    const double tol = 3.0e-3;
    int n_east = 0, n_west = 0;
    for (int f = 0; f < faces.n_faces(); ++f) {
        if (dir[f] == 2) {
            EXPECT_EQ(fg[f], Real(0.0)) << "roof face " << f;
            EXPECT_EQ(fl[f], Real(0.0)) << "roof face " << f;
            continue;
        }
        ASSERT_EQ(dir[f], 0);
        EXPECT_NEAR(fg[f], 0.5, tol) << "face " << f;
        const double h = (k[f] + 0.5) * 10.0;
        if (side[f] < 0) {   // the solid on the low side: the east side of the wall, in its shadow
            EXPECT_NEAR(fl[f], 0.5 * h / std::sqrt(h * h + L * L), tol) << "east face " << f << " at " << h << " m";
            ++n_east;
        } else {             // the west side, facing the sun: no shadow in its view
            EXPECT_EQ(fl[f], fg[f]) << "west face " << f;
            ++n_west;
        }
    }
    EXPECT_EQ(n_east, 32);
    EXPECT_EQ(n_west, 32);
    }
}

// Buildings in four boxes on a periodic domain: the sunlit part of every face's ground
// view against the count made here from the geometry, and the reflected light with the
// input on and off. With it off the face takes the uniform ground, as before the input.
TEST(IBSEBGroundViewShadow, BuildingSetAgainstTheGeometry)
{
    using namespace amrex;
    const int nx = 16, ny = 16;
    const double d = 10.0, dz = 10.0;
    const RealBox rb({0.0, 0.0, 0.0}, {160.0, 160.0, 80.0});
    const Box dom(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, 7));
    const Geometry geom(dom, rb, 0, {1, 1, 0});
    BoxArray ba(dom);
    ba.maxSize(IntVect(8, 8, 8));
    const DistributionMapping dm(ba);
    MultiFab blank(ba, dm, 1, 1);
    // A 40 m block across the four boxes, a 60 m tower and a low 20 m block.
    const std::vector<Cells> solid{{6, 9, 7, 8, 3}, {2, 2, 2, 3, 5}, {12, 13, 12, 12, 1}};
    fill_blanking(blank, geom, solid);

    const double zenith = 55.0, azimuth = 215.0;
    double sun[3];
    sun_vector(zenith, azimuth, sun);
    const Real dni = 800.0, dif = 100.0, alb_g = 0.25;

    for (const bool on : {true, false}) {
        IBSEBParams params;
        params.enable = true;
        params.sun_mode = "fixed";
        params.sun_zenith_deg = zenith;
        params.sun_azimuth_deg = azimuth;
        params.sw_direct_normal = dni;
        params.sw_diffuse = dif;
        params.albedo_ground = alb_g;
        // Ten azimuths: none at 45 degrees from the axes, where a ray from the centre of a
        // face at height 25 m would land exactly on the line between two columns.
        params.view_n_az = 10;
        params.view_n_el = 6;
        params.ground_view_shadow = on;
        IBFaceSet faces(params, 0);
        faces.build(blank, geom);
        faces.release_labels();
        faces.compute_view_fractions();
        faces.compute_shortwave(Real(0.0));

        const auto dir = host(faces.d_dir), side = host(faces.d_side);
        const auto xf = host(faces.d_xf), yf = host(faces.d_yf), zf = host(faces.d_zf);
        const auto fs = host(faces.d_f_sky), fg = host(faces.d_f_ground), fl = host(faces.d_f_ground_lit);
        const auto difin = host(faces.d_SW_diffuse_in);
        const Real dir_h = dni * static_cast<Real>(sun[2]);
        const Real tol = (sizeof(Real) == 8) ? Real(1.e-12) : Real(1.e-5);
        int n_part = 0, n_dark = 0, n_full = 0;
        for (int f = 0; f < faces.n_faces(); ++f) {
            if (!on) {
                // The uniform ground of the faces before the input, to the last bit.
                EXPECT_EQ(difin[f], fs[f] * dif + fg[f] * alb_g * (dir_h + dif)) << "face " << f;
                continue;
            }
            double e_ground = 0.0, e_lit = 0.0;
            // The shallowest of these rays lands within 12 face heights (660 m, four
            // domains) and the shadows are shorter than 90 m: images out to six domains
            // hold them all.
            expected_ground_view(dir[f], side[f], xf[f], yf[f], zf[f], params.view_n_az, params.view_n_el,
                                 sun, solid, nx, ny, d, dz, 6, e_ground, e_lit);
            // The counts are whole rays; the fractions carry the round-off of ERF's Real.
            const double count_tol = (sizeof(Real) == 8) ? 1.e-12 : 1.e-6;
            EXPECT_NEAR(fg[f], e_ground, count_tol) << "face " << f << " dir " << dir[f] << " side " << side[f];
            EXPECT_NEAR(fl[f], e_lit, count_tol) << "face " << f << " dir " << dir[f] << " side " << side[f]
                                               << " at (" << xf[f] << "," << yf[f] << "," << zf[f] << ")";
            EXPECT_LE(fl[f], fg[f]);
            const Real expect = fs[f] * dif + alb_g * (fg[f] * dif + fl[f] * dir_h);
            EXPECT_NEAR(difin[f], expect, tol * (std::abs(expect) + Real(1.0))) << "face " << f;
            if (fg[f] > 0.0) {
                if (fl[f] == 0.0) { ++n_dark; } else if (fl[f] < fg[f]) { ++n_part; } else { ++n_full; }
            }
        }
        if (on) {
            // Not vacuous: walls that see only shaded ground, only sunlit ground, and both.
            EXPECT_GT(n_part, 20);
            EXPECT_GT(n_full, 5);
            EXPECT_GE(n_dark + n_part, 30);
        }
    }
}

// At night the ground reflects no beam: nothing is counted, not even the rays that leave
// the closed domain for the far field's ground (which no shadow test would darken).
TEST(IBSEBGroundViewShadow, NoSunNoSunlitGround)
{
    using namespace amrex;
    const RealBox rb({0.0, 0.0, 0.0}, {160.0, 160.0, 80.0});
    const Box dom(IntVect(0, 0, 0), IntVect(15, 15, 7));
    const Geometry geom(dom, rb, 0, {0, 1, 0});
    BoxArray ba(dom);
    ba.maxSize(IntVect(8, 8, 8));
    const DistributionMapping dm(ba);
    MultiFab blank(ba, dm, 1, 1);
    fill_blanking(blank, geom, {{6, 9, 7, 8, 3}});
    IBSEBParams params;
    params.enable = true;
    params.sun_mode = "fixed";
    params.sun_zenith_deg = 100.0;
    params.sun_azimuth_deg = 90.0;
    params.view_n_az = 8;
    params.view_n_el = 4;
    params.ground_view_shadow = true;
    IBFaceSet faces(params, 0);
    faces.build(blank, geom);
    faces.release_labels();
    faces.compute_view_fractions();
    faces.compute_shortwave(Real(0.0));
    const auto fg = host(faces.d_f_ground), fl = host(faces.d_f_ground_lit), difin = host(faces.d_SW_diffuse_in);
    const auto far_rays = host(faces.d_gnd_far);
    int n_wall = 0, n_far = 0;
    for (int f = 0; f < faces.n_faces(); ++f) {
        EXPECT_EQ(fl[f], Real(0.0)) << "face " << f;
        EXPECT_EQ(difin[f], Real(0.0)) << "face " << f;
        if (fg[f] > 0.0) { ++n_wall; }
        n_far += far_rays[f];
    }
    EXPECT_GT(n_wall, 10);
    EXPECT_GT(n_far, 10);   // not vacuous: rays to the far field, which only the night keeps dark

    // The sun just above the horizon with no beam (a zero direct-normal irradiance): no
    // ray is cast and nothing is counted either.
    params.sun_zenith_deg = 60.0;
    params.sw_direct_normal = 0.0;
    IBFaceSet dim(params, 0);
    dim.build(blank, geom);
    dim.release_labels();
    dim.compute_view_fractions();
    dim.compute_shortwave(Real(0.0));
    for (const Real v : host(dim.d_f_ground_lit)) { EXPECT_EQ(v, Real(0.0)); }
}

// The shortwave with the input needs the ground landings of compute_view_fractions().
TEST(IBSEBGroundViewShadow, LandingsAreBuiltWithTheViewFractions)
{
    using namespace amrex;
    const RealBox rb({0.0, 0.0, 0.0}, {160.0, 160.0, 80.0});
    const Box dom(IntVect(0, 0, 0), IntVect(15, 15, 7));
    const Geometry geom(dom, rb, 0, {1, 1, 0});
    BoxArray ba(dom);
    ba.maxSize(IntVect(8, 8, 8));
    const DistributionMapping dm(ba);
    MultiFab blank(ba, dm, 1, 1);
    fill_blanking(blank, geom, {{6, 9, 7, 8, 3}});
    IBSEBParams params;
    params.enable = true;
    params.view_n_az = 8;
    params.view_n_el = 4;
    params.ground_view_shadow = true;
    IBFaceSet faces(params, 0);
    faces.build(blank, geom);
    faces.release_labels();
    EXPECT_EQ(faces.d_gnd_land.size(), 0U);
    faces.compute_view_fractions();
    // Every ground ray lands in one of the face's kept columns or goes to the far field,
    // and a column is kept once, with the number of rays that land in it.
    const auto fg = host(faces.d_f_ground);
    const auto start = host(faces.d_gnd_start), far_rays = host(faces.d_gnd_far), rays = host(faces.d_gnd_rays);
    const auto land = host(faces.d_gnd_land);
    ASSERT_EQ(static_cast<int>(start.size()), faces.n_faces() + 1);
    ASSERT_EQ(static_cast<int>(land.size()), start.back());
    ASSERT_EQ(rays.size(), land.size());
    int n_land = 0, n_rays = 0;
    for (int f = 0; f < faces.n_faces(); ++f) {
        int n = far_rays[f];
        for (int m = start[f]; m < start[f + 1]; ++m) {
            EXPECT_GE(rays[m], 1);
            EXPECT_GE(land[m], 65536U);                       // a landing, not one of the two marks
            if (m > start[f]) { EXPECT_GT(land[m], land[m - 1]) << "face " << f; }   // no column twice
            n += rays[m];
        }
        EXPECT_EQ(fg[f], n * ibseb::hemisphere_ray_weight(8, 4)) << "face " << f;
        n_land += start[f + 1] - start[f];
        n_rays += n - far_rays[f];
    }
    EXPECT_GT(n_land, 100);
    EXPECT_GT(n_rays, n_land + 20);   // not vacuous: many rays share a column
    // The marks of the rays that are not kept are no landing's word.
    EXPECT_LT(ibseb::landing_none, ibseb::pack_landing(-ibseb::landing_offset_max, -ibseb::landing_offset_max));
    EXPECT_LT(ibseb::landing_far, ibseb::pack_landing(-ibseb::landing_offset_max, -ibseb::landing_offset_max));
    // The input off: nothing is kept.
    params.ground_view_shadow = false;
    IBFaceSet plain(params, 0);
    plain.build(blank, geom);
    plain.release_labels();
    plain.compute_view_fractions();
    EXPECT_EQ(plain.d_gnd_land.size(), 0U);
    EXPECT_EQ(plain.d_gnd_rays.size(), 0U);
}

// A landing's offsets survive the packing, at both ends of what is kept, and what is
// kept ends at the map's reach or at 32767 columns, whichever is nearer.
TEST(IBSEBGroundViewShadow, LandingOffsetsPackAndUnpack)
{
    for (const int di : {-32767, -115, -1, 0, 1, 7, 32767}) {
        for (const int dj : {-32767, -3, 0, 2, 32767}) {
            int a = 0, b = 0;
            ibseb::unpack_landing(ibseb::pack_landing(di, dj), a, b);
            EXPECT_EQ(a, di);
            EXPECT_EQ(b, dj);
        }
    }
    ibseb::RayMap m;
    m.reach_x = 12; m.reach_y = 100000;
    EXPECT_TRUE(ibseb::landing_kept(12, -32767, m));
    EXPECT_FALSE(ibseb::landing_kept(13, 0, m));
    EXPECT_FALSE(ibseb::landing_kept(-13, 0, m));
    EXPECT_FALSE(ibseb::landing_kept(0, 32768, m));
    EXPECT_FALSE(ibseb::landing_kept(0, -32768, m));
}

// The shadow walk from a landing column stays within the map's reach of the face's own
// column, in x and in y: every rank's patch holds that window around its faces, so the
// answer does not depend on where the patch ends. A building inside the window shades;
// the same building is not seen from a face whose window ends before it.
namespace {

struct HandMap
{
    std::vector<int> tops;
    ibseb::RayMap m;
    // A map of 10 m columns over the patch [ilo, ihi] x [0, 2] of a domain of nx x 3
    // columns, periodic in y, and wrapped in x when the patch is the whole domain.
    HandMap (int ilo, int ihi, int nx, bool wrap_x, int reach)
        : tops(static_cast<size_t>((ihi - ilo + 1) * 3), -1)
    {
        const amrex::Box patch(amrex::IntVect(ilo, 0, 0), amrex::IntVect(ihi, 2, 0));
        m.top = amrex::Array4<const int>(tops.data(), amrex::begin(patch), amrex::end(patch), 1);
        m.ilo = 0; m.jlo = 0; m.nx = nx; m.ny = 3;
        m.wrap_x = wrap_x; m.wrap_y = true;   // the street is the same all along y
        m.x_lo = 0.0; m.y_lo = 0.0; m.dx = 10.0; m.dy = 10.0; m.dz = 10.0;
        m.z_ground = 0.0; m.z_max = 100.0;
        m.max_steps = 2 * reach;   // the reach in x plus y, as IBFaceSet::build_ray_map() sets it
        m.reach_x = reach; m.reach_y = reach;
    }
    // A tower in column i (and its periodic images inside the patch) whose highest
    // solid cell is kt: 100 m by default.
    void tower (int i, int kt = 9)
    {
        for (int ci = m.top.begin[0]; ci < m.top.end[0]; ++ci) {
            if (((ci - i) % m.nx + m.nx) % m.nx != 0) { continue; }
            for (int j = 0; j < 3; ++j) { tops[static_cast<size_t>((ci - m.top.begin[0]) + (m.top.end[0] - m.top.begin[0]) * j)] = kt; }
        }
    }
};

// The sunlit ground rays of a wall face at (xf, 15 m, zf) looking along +x from column
// ci_f, counted with the helpers the faces use.
int sunlit_ground_rays (const ibseb::RayMap& m, amrex::Real xf, amrex::Real zf, int ci_f, int n_az, int n_el,
                        amrex::Real sx, amrex::Real sy, amrex::Real sz, int& n_ground, int& n_beyond)
{
    int n_lit = 0;
    n_ground = 0;
    n_beyond = 0;   // kept landings from which the unbounded walk would meet the far tower
    for (int ie = 0; ie < n_el; ++ie) {
        for (int ia = 0; ia < n_az; ++ia) {
            amrex::Real rx, ry, rz;
            ibseb::hemisphere_direction(0, 1, ia, ie, n_az, n_el, rx, ry, rz);
            bool landed = false;
            int li = 0, lj = 0;
            if (ibseb::ray_hit_landing(xf, amrex::Real(15.0), zf, rx, ry, rz, m, landed, li, lj) != ibseb::RAY_GROUND) { continue; }
            ++n_ground;
            if (!(landed && ibseb::landing_kept(li - ci_f, lj - 1, m))) { ++n_lit; continue; }
            if (li >= 28) { ++n_beyond; }
            if (!ibseb::ray_blocked_near(li, lj, sx, sy, sz, m, ci_f, 1)) { ++n_lit; }
        }
    }
    return n_lit;
}

} // namespace

TEST(IBSEBGroundViewShadow, ShadowWalkKeepsToTheReachOfTheFace)
{
    using namespace amrex;
    HandMap h(0, 39, 40, false, 12);
    h.tower(30);
    const ibseb::RayMap& m = h.m;
    // The sun in the east at 45 degrees: from the centre of column 22 the ray enters the
    // tower's column 75 m on, at a height of 75 m, below its top.
    const Real s = std::sqrt(Real(0.5));
    EXPECT_TRUE(ibseb::ray_blocked(Real(225.0), Real(15.0), Real(0.0), s, Real(0.0), s, m));
    // Seen from a face in column 20, the tower is 10 columns away: inside the reach.
    EXPECT_TRUE(ibseb::ray_blocked_near(22, 1, s, Real(0.0), s, m, 20, 1));
    // From a face in column 17 it is 13 columns away: beyond the reach, not seen.
    EXPECT_FALSE(ibseb::ray_blocked_near(22, 1, s, Real(0.0), s, m, 17, 1));
    // A landing column itself beyond the reach of the face is open, sunlit ground.
    EXPECT_FALSE(ibseb::ray_blocked_near(29, 1, s, Real(0.0), s, m, 16, 1));
    // The sun at or below the horizon lights nothing.
    EXPECT_TRUE(ibseb::ray_blocked_near(22, 1, Real(1.0), Real(0.0), Real(0.0), m, 20, 1));
}

// The ground a face sees takes the shadow by the walk of the ground's own shadow, from the
// same point: inside the reach of the face, ray_blocked_near() is ray_blocked() from the
// centre of the column's ground, for suns along the street, across it, on the diagonal
// of the columns (the ray leaves each column through its corner) and low.
TEST(IBSEBGroundViewShadow, ShadowWalkIsTheGroundShadowsWalk)
{
    using namespace amrex;
    HandMap h(0, 63, 64, false, 100);   // the reach holds the whole map
    h.tower(30); h.tower(41, 2); h.tower(12, 5);
    const ibseb::RayMap& m = h.m;
    const Real pi = amrex::Math::pi<Real>();
    int n_shaded = 0, n_sunlit = 0;
    for (const Real zen_deg : {Real(20.0), Real(45.0), Real(70.0), Real(85.0)}) {
        for (const Real az_deg : {Real(90.0), Real(270.0), Real(45.0), Real(225.0), Real(10.0), Real(133.0)}) {
            const Real zen = zen_deg * pi / Real(180.0), az = az_deg * pi / Real(180.0);
            const Real sx = std::sin(az) * std::sin(zen), sy = std::cos(az) * std::sin(zen), sz = std::cos(zen);
            for (int ci = 0; ci < 64; ++ci) {
                for (int cj = 0; cj < 3; ++cj) {
                    const bool ground = ibseb::ray_blocked(m.centre_x(ci), m.centre_y(cj), m.z_ground, sx, sy, sz, m);
                    EXPECT_EQ(ibseb::ray_blocked_near(ci, cj, sx, sy, sz, m, 32, 1), ground)
                        << "column " << ci << "," << cj << " sun " << zen_deg << "/" << az_deg;
                    (ground ? n_shaded : n_sunlit) += 1;
                }
            }
        }
    }
    EXPECT_GT(n_shaded, 200);
    EXPECT_GT(n_sunlit, 200);
    EXPECT_EQ(m.centre_x(22), Real(225.0));
    EXPECT_EQ(m.centre_y(1), Real(15.0));
}

// The same face on two maps of the same periodic street: the whole domain, wrapped, as
// a rank that holds all of it has it, and a patch that ends at the reach of the face, as
// a rank that holds only the face's column has it. The maps are sized as
// IBFaceSet::build_ray_map() sizes them: with 8 x 2 view rays and a sun floor of 10
// degrees over 100 m towers the reach is 58 columns of 10 m (the sun's; the shallowest
// view ray reaches 52). The whole domain holds a 100 m tower 64 columns from the face,
// which the patch does not; the ray to the sun from the far part of the face's ground
// view would meet it. The sunlit ground rays are the same on both maps, as that tower
// is beyond the face's reach.
TEST(IBSEBGroundViewShadow, SameAnswerOnAPatchAndOnTheWholeDomain)
{
    using namespace amrex;
    const int nx = 256, reach = 58, ci_f = 21, n_az = 8, n_el = 2;
    HandMap whole(0, nx - 1, nx, true, reach);
    HandMap patch(ci_f - reach, ci_f + reach, nx, false, reach);
    for (HandMap* h : {&whole, &patch}) {
        // 100 m towers at 20 (the face's own building) and at 85, beyond the reach.
        h->tower(20); h->tower(85);
        // A 10 m block at 27, inside the reach, which the higher faces see over.
        h->tower(27, 0);
    }
    // The sun in the east, 10 degrees up: the block's shadow covers columns 22 to 26, and
    // from the centre of column 45 the ray to the sun enters the tower at 85 at 70 m.
    const Real zen = Real(80.0) * amrex::Math::pi<Real>() / Real(180.0);
    const Real sx = std::sin(zen), sy = 0.0, sz = std::cos(zen);
    EXPECT_TRUE(ibseb::ray_blocked(Real(455.0), Real(15.0), Real(0.0), sx, sy, sz, whole.m));
    EXPECT_FALSE(ibseb::ray_blocked(Real(455.0), Real(15.0), Real(0.0), sx, sy, sz, patch.m));
    int n_shaded_total = 0, n_beyond_total = 0;
    for (const Real zf : {Real(5.0), Real(35.0), Real(75.0)}) {
        int g_whole = 0, g_patch = 0, b_whole = 0, b_patch = 0;
        const int lit_whole = sunlit_ground_rays(whole.m, Real(210.0), zf, ci_f, n_az, n_el, sx, sy, sz, g_whole, b_whole);
        const int lit_patch = sunlit_ground_rays(patch.m, Real(210.0), zf, ci_f, n_az, n_el, sx, sy, sz, g_patch, b_patch);
        EXPECT_EQ(g_whole, g_patch) << "face at " << zf << " m";
        EXPECT_EQ(lit_whole, lit_patch) << "face at " << zf << " m";
        EXPECT_EQ(b_whole, b_patch) << "face at " << zf << " m";
        EXPECT_GT(g_whole, 0);
        n_shaded_total += g_whole - lit_whole;
        n_beyond_total += b_whole;
    }
    // Not vacuous: the block inside the reach shades, and rays land beyond it, where the
    // walk without the window would find the far tower on one map only.
    EXPECT_GE(n_shaded_total, 2);
    EXPECT_GE(n_beyond_total, 2);
}

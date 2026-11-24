#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include <ERF_Advection.H>
#include <ERF_IndexDefines.H>

#include <gtest/gtest.h>

#include "../ERF_GTestAssertions.H"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <vector>

using amrex::Array4;
using amrex::Box;
using amrex::FArrayBox;
using amrex::GpuArray;
using amrex::IntVect;
using amrex::Real;

// The open-boundary advection kernels with map factors (erf-model/ERF#5359, section 3).
// Two properties, neither of which the kernels had without the map factors:
//  - A uniform map factor m is a mesh m times finer in x and y: every kernel gives with
//    map factor m what it gives with no map factors and dx/m, dy/m.
//  - With map factors that vary, and the mass fluxes (rho u/m_y, rho v/m_x, and the
//    averaged fluxes of the scalars) held fixed, the tendency divided by the prefactor of
//    its variable (m_x m_y at the location of the variable; m_x or m_y for the normal
//    momentum, whose gradient is m d/dxi) does not depend on the map factors.

namespace {

constexpr int kN = 6;      // cells in each direction
constexpr int kGrow = 2;   // ghost cells on every array

Box domain () { return Box(IntVect(0), IntVect(kN-1)); }

FArrayBox make_fab (const IntVect& nodal, int ncomp = 1)
{
    Box b = amrex::convert(domain(), nodal);
    b.grow(kGrow);
    return FArrayBox(b, ncomp, amrex::The_Pinned_Arena());
}

FArrayBox make_mf (const IntVect& nodal)
{
    Box b = amrex::convert(domain(), nodal);
    b.setRange(2, 0);
    b.grow(IntVect(kGrow, kGrow, 0));
    return FArrayBox(b, 1, amrex::The_Pinned_Arena());
}

// Smooth fields that differ in every direction
Real field (int i, int j, int k, Real a, Real b)
{
    return a + Real(0.3)*std::sin(Real(0.37)*i + b) * std::cos(Real(0.29)*j - b)
             + Real(0.2)*std::cos(Real(0.41)*k + b);
}

void fill_fab (FArrayBox& fab, const std::function<Real(int,int,int)>& f, int comp = 0)
{
    auto a = fab.array();
    const Box& b = fab.box();
    for (int k = b.smallEnd(2); k <= b.bigEnd(2); ++k) {
    for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j) {
    for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i) {
        a(i,j,k,comp) = f(i,j,k);
    }}}
}

// The map factors at the six locations ERF keeps them
struct MapFactors
{
    FArrayBox mx, my, ux, uy, vx, vy;

    MapFactors ()
        : mx(make_mf(IntVect(0,0,0))), my(make_mf(IntVect(0,0,0))),
          ux(make_mf(IntVect(1,0,0))), uy(make_mf(IntVect(1,0,0))),
          vx(make_mf(IntVect(0,1,0))), vy(make_mf(IntVect(0,1,0))) {}

    void set_uniform (Real m)
    {
        for (FArrayBox* f : {&mx, &my, &ux, &uy, &vx, &vy}) { f->setVal<amrex::RunOn::Host>(m); }
    }

    // Anisotropic (m_x != m_y) and varying in x and y; phase tells two sets apart
    void set_varying (Real phase)
    {
        auto mfx = [=] (Real x, Real y) { return Real(1.1) + Real(0.15)*std::sin(Real(0.5)*x + phase)
                                                           + Real(0.1)*std::cos(Real(0.4)*y - phase); };
        auto mfy = [=] (Real x, Real y) { return Real(0.9) + Real(0.12)*std::cos(Real(0.45)*x - phase)
                                                           + Real(0.08)*std::sin(Real(0.6)*y + phase); };
        // Positions in cells: centers at i+1/2, x-faces at i
        fill_fab(mx, [=] (int i, int j, int) { return mfx(Real(i)+Real(0.5), Real(j)+Real(0.5)); });
        fill_fab(my, [=] (int i, int j, int) { return mfy(Real(i)+Real(0.5), Real(j)+Real(0.5)); });
        fill_fab(ux, [=] (int i, int j, int) { return mfx(Real(i)          , Real(j)+Real(0.5)); });
        fill_fab(uy, [=] (int i, int j, int) { return mfy(Real(i)          , Real(j)+Real(0.5)); });
        fill_fab(vx, [=] (int i, int j, int) { return mfx(Real(i)+Real(0.5), Real(j)          ); });
        fill_fab(vy, [=] (int i, int j, int) { return mfy(Real(i)+Real(0.5), Real(j)          ); });
    }
};

// The state the kernels read.  Where a test holds the mass fluxes fixed it sets rho_u and
// rho_v from them and the map factors (set_momenta).
struct State
{
    FArrayBox u, v, w, rho_u, rho_v, rho_w, omega, cell, scal, ax, ay, az, detJ;
    FArrayBox avg_x, avg_y, avg_z;

    State ()
        : u(make_fab(IntVect(1,0,0))), v(make_fab(IntVect(0,1,0))), w(make_fab(IntVect(0,0,1))),
          rho_u(make_fab(IntVect(1,0,0))), rho_v(make_fab(IntVect(0,1,0))), rho_w(make_fab(IntVect(0,0,1))),
          omega(make_fab(IntVect(0,0,1))), cell(make_fab(IntVect(0), NVAR_max)), scal(make_fab(IntVect(0))),
          ax(make_fab(IntVect(1,0,0))), ay(make_fab(IntVect(0,1,0))), az(make_fab(IntVect(0,0,1))),
          detJ(make_fab(IntVect(0))),
          avg_x(make_fab(IntVect(1,0,0))), avg_y(make_fab(IntVect(0,1,0))), avg_z(make_fab(IntVect(0,0,1)))
    {
        // Velocities near 30 m/s or below, so the normal kernel's upwinding picks the same side
        fill_fab(u, [] (int i, int j, int k) { return field(i,j,k, Real( 4.0), Real(0.3)); });
        fill_fab(v, [] (int i, int j, int k) { return field(i,j,k, Real(-3.0), Real(1.1)); });
        fill_fab(w, [] (int i, int j, int k) { return field(i,j,k, Real( 0.5), Real(2.0)); });
        fill_fab(cell, [] (int i, int j, int k) { return field(i,j,k, Real(1.2), Real(0.7)); }, Rho_comp);
        fill_fab(scal, [] (int i, int j, int k) { return field(i,j,k, Real(2.0), Real(1.7)); });
        // Areas and Jacobian of a terrain-fitted mesh, near one
        fill_fab(ax,   [] (int i, int j, int k) { return field(i,j,k, Real(1.0), Real(0.9)) - Real(0.2); });
        fill_fab(ay,   [] (int i, int j, int k) { return field(i,j,k, Real(1.0), Real(2.3)) - Real(0.2); });
        az.setVal<amrex::RunOn::Host>(Real(1.0));
        fill_fab(detJ, [] (int i, int j, int k) { return field(i,j,k, Real(1.0), Real(0.4)) - Real(0.1); });
        fill_fab(omega, [] (int i, int j, int k) { return field(i,j,k, Real(0.4), Real(1.3)); });
        fill_fab(avg_x, [] (int i, int j, int k) { return field(i,j,k, Real( 3.5), Real(0.2)); });
        fill_fab(avg_y, [] (int i, int j, int k) { return field(i,j,k, Real(-2.5), Real(0.6)); });
        fill_fab(avg_z, [] (int i, int j, int k) { return field(i,j,k, Real( 0.6), Real(1.9)); });
    }

    // rho u = m_y F_u and rho v = m_x F_v, for the given mass fluxes F_u, F_v
    void set_momenta (const MapFactors& m)
    {
        auto uy = m.uy.const_array();
        auto vx = m.vx.const_array();
        fill_fab(rho_u, [=] (int i, int j, int k) { return uy(i,j,0) * field(i,j,k, Real( 4.8), Real(0.5)); });
        fill_fab(rho_v, [=] (int i, int j, int k) { return vx(i,j,0) * field(i,j,k, Real(-3.6), Real(1.4)); });
    }
};

GpuArray<Real,AMREX_SPACEDIM> cell_size_inv (Real fac_xy = Real(1.0))
{
    return {{Real(0.45)*fac_xy, Real(0.8)*fac_xy, Real(1.15)}};
}

// The kernels write their tendency into one array; each case names its kernel, its box, the
// location of the variable it updates, and how to call it
enum struct Loc { x_face, y_face, z_face, cell };

struct Case
{
    std::string name;
    Loc loc;
    Box bx;
    std::function<void(const Box&, const Array4<Real>&, const State&, const MapFactors&,
                       const GpuArray<Real,AMREX_SPACEDIM>&)> run;
};

IntVect nodal_of (Loc l)
{
    switch (l) {
    case Loc::x_face: return IntVect(1,0,0);
    case Loc::y_face: return IntVect(0,1,0);
    case Loc::z_face: return IntVect(0,0,1);
    default:          return IntVect(0);
    }
}

// The prefactor that multiplies the flux divergence of the variable (Loc), for the
// tangential kernels, at (i,j)
Real prefactor (Loc l, const MapFactors& m, int i, int j)
{
    switch (l) {
    case Loc::x_face: return m.ux.const_array()(i,j,0) * m.uy.const_array()(i,j,0);
    case Loc::y_face: return m.vx.const_array()(i,j,0) * m.vy.const_array()(i,j,0);
    default:          return m.mx.const_array()(i,j,0) * m.my.const_array()(i,j,0);
    }
}

std::vector<Case> cases ()
{
    const Box d = domain();
    const int lo = d.smallEnd(0), hi = d.bigEnd(0);
    std::vector<Case> c;

    // Normal momentum on an x-lo and a y-hi open boundary
    {
        Box b = amrex::surroundingNodes(d, 0); b.setRange(0, lo);
        c.push_back({"Normal_xlo", Loc::x_face, b,
            [] (const Box& bx, const Array4<Real>& r, const State& s, const MapFactors& m, const auto& dxi) {
                AdvectionSrcForOpenBC_Normal(bx, 0, r, s.u.const_array(), s.cell.const_array(),
                                             m.ux.const_array(), dxi, true);
            }});
    }
    {
        Box b = amrex::surroundingNodes(d, 1); b.setRange(1, hi+1);
        c.push_back({"Normal_yhi", Loc::y_face, b,
            [] (const Box& bx, const Array4<Real>& r, const State& s, const MapFactors& m, const auto& dxi) {
                AdvectionSrcForOpenBC_Normal(bx, 1, r, s.v.const_array(), s.cell.const_array(),
                                             m.vy.const_array(), dxi, false);
            }});
    }
    // Tangential momenta: v on an x-lo boundary, u on a y-hi boundary
    {
        Box b = amrex::surroundingNodes(d, 1); b.setRange(0, lo); b.setRange(2, 1, kN-2);
        c.push_back({"Tangent_Ymom_xlo", Loc::y_face, b,
            [] (const Box& bx, const Array4<Real>& r, const State& s, const MapFactors& m, const auto& dxi) {
                AdvectionSrcForOpenBC_Tangent_Ymom(bx, 0, r, s.v.const_array(),
                                                   s.rho_u.const_array(), s.rho_v.const_array(), s.omega.const_array(),
                                                   s.ay.const_array(), s.az.const_array(), s.detJ.const_array(),
                                                   m.vx.const_array(), m.vy.const_array(), m.uy.const_array(),
                                                   dxi, true);
            }});
    }
    {
        Box b = amrex::surroundingNodes(d, 0); b.setRange(1, hi); b.setRange(2, 1, kN-2);
        c.push_back({"Tangent_Xmom_yhi", Loc::x_face, b,
            [] (const Box& bx, const Array4<Real>& r, const State& s, const MapFactors& m, const auto& dxi) {
                AdvectionSrcForOpenBC_Tangent_Xmom(bx, 1, r, s.u.const_array(),
                                                   s.rho_u.const_array(), s.rho_v.const_array(), s.omega.const_array(),
                                                   s.ax.const_array(), s.az.const_array(), s.detJ.const_array(),
                                                   m.ux.const_array(), m.uy.const_array(), m.vx.const_array(),
                                                   dxi, false);
            }});
    }
    // w on an x-lo boundary and in the x-lo/y-hi corner
    for (const OpenSide ys : {OpenSide::none, OpenSide::hi}) {
        Box b = amrex::surroundingNodes(d, 2); b.setRange(0, lo); b.setRange(2, 1, kN-1);
        if (ys == OpenSide::hi) { b.setRange(1, hi); }
        c.push_back({ys == OpenSide::none ? "Tangent_Zmom_xlo" : "Tangent_Zmom_corner", Loc::z_face, b,
            [ys] (const Box& bx, const Array4<Real>& r, const State& s, const MapFactors& m, const auto& dxi) {
                AdvectionSrcForOpenBC_Tangent_Zmom(bx, OpenSide::lo, ys, r, s.w.const_array(),
                                                   s.rho_u.const_array(), s.rho_v.const_array(), s.omega.const_array(),
                                                   s.ax.const_array(), s.ay.const_array(), s.az.const_array(),
                                                   s.detJ.const_array(),
                                                   m.mx.const_array(), m.my.const_array(),
                                                   m.uy.const_array(), m.vx.const_array(),
                                                   dxi, domain().bigEnd(2));
            }});
    }
    // A scalar on an x-lo boundary and in the x-lo/y-hi corner
    for (const OpenSide ys : {OpenSide::none, OpenSide::hi}) {
        Box b = d; b.setRange(0, lo); b.setRange(2, 1, kN-2);
        if (ys == OpenSide::hi) { b.setRange(1, hi); }
        c.push_back({ys == OpenSide::none ? "Tangent_Scalars_xlo" : "Tangent_Scalars_corner", Loc::cell, b,
            [ys] (const Box& bx, const Array4<Real>& r, const State& s, const MapFactors& m, const auto& dxi) {
                AdvectionSrcForOpenBC_Tangent_Scalars(bx, OpenSide::lo, ys, 0, 0, 1, r, s.scal.const_array(),
                                                      s.avg_x.const_array(), s.avg_y.const_array(),
                                                      s.avg_z.const_array(), s.detJ.const_array(),
                                                      m.mx.const_array(), m.my.const_array(), dxi);
            }});
    }
    return c;
}

FArrayBox run_case (const Case& c, const State& s, const MapFactors& m,
                    const GpuArray<Real,AMREX_SPACEDIM>& dxi)
{
    FArrayBox r(amrex::convert(domain(), nodal_of(c.loc)), 1, amrex::The_Pinned_Arena());
    r.setVal<amrex::RunOn::Host>(Real(0.0));
    c.run(c.bx, r.array(), s, m, dxi);
    amrex::Gpu::streamSynchronize();
    return r;
}

// Each kernel's tendencies, to a few hundred ulps of their size
Real tolerance (Real scale) { return Real(500) * std::numeric_limits<Real>::epsilon() * scale; }

} // namespace

// Motivation: a uniform map factor m only makes the physical mesh m times finer in x and y,
// so every open-boundary kernel must give with map factor m what it gives with none on that
// finer mesh.  Without map factors in the kernels the result did not depend on m at all.
// The averaged scalar fluxes come from the substeps, which divide them by m_y, m_x and
// m_x m_y; the same physical flux is held here by scaling them so.
TEST(OpenBCMapFactors, UniformMapFactorIsAFinerMesh)
{
    const Real m = Real(0.8);
    MapFactors ones, mm;
    ones.set_uniform(Real(1.0));
    mm.set_uniform(m);

    State s_ref;
    s_ref.set_momenta(ones);
    State s_m;
    s_m.set_momenta(ones);   // the same momenta (the same physical flow)
    s_m.avg_x.mult<amrex::RunOn::Host>(Real(1.0)/m);
    s_m.avg_y.mult<amrex::RunOn::Host>(Real(1.0)/m);
    s_m.avg_z.mult<amrex::RunOn::Host>(Real(1.0)/(m*m));

    for (const Case& c : cases()) {
        SCOPED_TRACE(c.name);
        const FArrayBox ref = run_case(c, s_ref, ones, cell_size_inv(m));
        const FArrayBox got = run_case(c, s_m,   mm,   cell_size_inv());
        const auto r = ref.const_array();
        const auto g = got.const_array();
        Real scale = 0;
        amrex::Loop(c.bx, [&] (int i, int j, int k) { scale = std::max(scale, std::abs(r(i,j,k))); });
        ASSERT_GT(scale, Real(0));
        amrex::Loop(c.bx, [&] (int i, int j, int k) {
            ERF_EXPECT_NEAR(g(i,j,k), r(i,j,k), tolerance(scale));
        });
    }
}

// Motivation: with map factors that vary in x and y, the flux divergence of each variable is
// the prefactor m_x m_y (at the variable) times the divergence of the mass fluxes, rho u/m_y
// and rho v/m_x, and of the averaged scalar fluxes, which carry their map factors already.
// With those held fixed and the map factors changed, the tendency over the prefactor stays
// the same (the normal momentum's radiation gradient is m_x d/dxi, so there it is m_x).
// Without map factors the kernels advected rho u and rho v as if they were mass fluxes, and
// the scalar kernel left the vertical term without m_x m_y.  Omega is zero so that the
// vertical momentum fluxes, which have no prefactor, drop out.
TEST(OpenBCMapFactors, VaryingMapFactorsScaleTheMassFluxDivergence)
{
    MapFactors m1, m2;
    m1.set_varying(Real(0.0));
    m2.set_varying(Real(1.3));

    State s1, s2;
    s1.set_momenta(m1);
    s2.set_momenta(m2);
    s1.omega.setVal<amrex::RunOn::Host>(Real(0.0));
    s2.omega.setVal<amrex::RunOn::Host>(Real(0.0));

    for (const Case& c : cases()) {
        SCOPED_TRACE(c.name);
        const FArrayBox r1 = run_case(c, s1, m1, cell_size_inv());
        const FArrayBox r2 = run_case(c, s2, m2, cell_size_inv());
        const auto a1 = r1.const_array();
        const auto a2 = r2.const_array();
        const bool normal = (c.name.rfind("Normal", 0) == 0);
        auto pref = [&] (const MapFactors& m, int i, int j) {
            if (!normal) { return prefactor(c.loc, m, i, j); }
            return (c.loc == Loc::x_face) ? m.ux.const_array()(i,j,0) : m.vy.const_array()(i,j,0);
        };
        Real scale = 0;
        amrex::Loop(c.bx, [&] (int i, int j, int k) {
            scale = std::max(scale, std::abs(a1(i,j,k) / pref(m1,i,j)));
        });
        ASSERT_GT(scale, Real(0));
        amrex::Loop(c.bx, [&] (int i, int j, int k) {
            ERF_EXPECT_NEAR(a2(i,j,k) / pref(m2,i,j), a1(i,j,k) / pref(m1,i,j), tolerance(scale));
        });
    }
}

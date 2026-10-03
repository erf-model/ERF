#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GMRES.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include "ERF_Constants.H"
#include "ERF_NumericalConstants.H"
#include "ERF_MLTerrainPoisson.H"
#include "ERF_TerrainPoisson.H"
#include "ERF_SolverUtils.H"
#include "ERF_Utils.H"

// Motivation: MLTerrainPoisson is the multigrid form of the terrain-fitted Poisson
// operator the anelastic projection solves.  It must apply exactly the stencil the
// GMRES operator applies (so the two solvers converge to the same discrete solution
// and the velocity correction stays discretely divergence-free), its coarse levels
// must discretize the same terrain, its smoother must relax with the true column
// coefficients of the operator including every boundary fold, and the whole thing
// must converge as a solver on singular and non-singular problems.

using namespace amrex;

namespace {

//
// A test problem: a hill in x on a terrain-fitted mesh
//
struct HillCase
{
    int nx, ny, nz;
    Real Lx, Ly, H;
    bool per_x, per_y;
    int hill;                       // 0: cosine hill (periodic in x), 1: Witch of Agnesi
    Real hmax, a;                   // hill height and half-width
    Array<std::string,2*AMREX_SPACEDIM> bc_names;   // ERF domain boundary names
    IntVect max_grid_size;
    bool semicoarsen = false;       // let the lateral directions coarsen on after z stops
    bool lshape = false;            // solve on an L-shaped union of boxes inside the domain (coarse/fine Neumann)
};

// Problem built on the case: mesh, metrics and the operator
struct Problem
{
    Geometry geom;
    BoxArray ba;
    DistributionMapping dm;
    MultiFab znd, ax, ay, az, dJ;
    Array<std::string,2*AMREX_SPACEDIM> bc_names;
    std::unique_ptr<MLTerrainPoisson> op;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real surface_height (Real x, Real Lx, Real hmax, Real a, int hill) noexcept
{
    Real xc = x - Real(0.5)*Lx;
    if (hill == 0) {
        return Real(0.5)*hmax*(Real(1.0) + std::cos(Real(2.0)*PI*xc/Lx));
    } else {
        return hmax / (Real(1.0) + (xc/a)*(xc/a));
    }
}

Geometry
make_geom (const HillCase& c)
{
    const Box domain(IntVect(0,0,0), IntVect(c.nx-1, c.ny-1, c.nz-1));
    const RealBox rb(Real(0.), Real(0.), Real(0.), c.Lx, c.Ly, c.H);
    return Geometry(domain, rb, CoordSys::cartesian, {c.per_x ? 1 : 0, c.per_y ? 1 : 0, 0});
}

// Basic terrain following over the hill: valid nodes from the formula, ghost nodes by
// the ERF rules (periodic image, lateral clamp, vertical linear extrapolation)
void
fill_terrain (MultiFab& znd, const Geometry& geom, const HillCase& c)
{
    const Real dx = geom.CellSize(0);
    const Real dz = geom.CellSize(2);
    const Real Lx = c.Lx, H = c.H, hmax = c.hmax, a = c.a;
    const int hill = c.hill;
    for (MFIter mfi(znd); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> const& z = znd.array(mfi);
        // the hill runs in x only, so every j column is the same
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real zs = surface_height(i*dx, Lx, hmax, a, hill);
            z(i,j,k) = zs + (H - zs) * (k*dz) / H;
        });
    }
    znd.FillBoundary(geom.periodicity());
    MLTerrainPoisson::fill_uncovered_nodes(znd, geom);
}

Problem
build_problem (const HillCase& c)
{
    Problem p;
    p.geom = make_geom(c);
    if (c.lshape) {
        // a refined level that is the union of a wide southern box and a north-eastern one
        const Box& d = p.geom.Domain();
        const int nx = d.length(0), ny = d.length(1), nz = d.length(2);
        Box south(IntVect(nx/4, ny/4, 0), IntVect(3*nx/4-1, ny/2-1, nz-1));
        Box northeast(IntVect(nx/2, ny/2, 0), IntVect(3*nx/4-1, 3*ny/4-1, nz-1));
        BoxList bl; bl.push_back(south); bl.push_back(northeast);
        p.ba = BoxArray(bl);
    } else {
        p.ba = BoxArray(p.geom.Domain());
    }
    p.ba.maxSize(c.max_grid_size);
    p.dm = DistributionMapping(p.ba);
    p.bc_names = c.bc_names;

    p.znd.define(amrex::convert(p.ba, IntVect(1)), p.dm, 1, 1);
    fill_terrain(p.znd, p.geom, c);

    p.dJ.define(p.ba, p.dm, 1, 1);
    p.ax.define(amrex::convert(p.ba, IntVect(1,0,0)), p.dm, 1, 1);
    p.ay.define(amrex::convert(p.ba, IntVect(0,1,0)), p.dm, 1, 1);
    p.az.define(amrex::convert(p.ba, IntVect(0,0,1)), p.dm, 1, 1);
    p.dJ.setVal(1.0); p.ax.setVal(1.0); p.ay.setVal(1.0); p.az.setVal(1.0);
    make_J(p.geom, p.znd, p.dJ);
    make_areas(p.geom, p.znd, p.ax, p.ay, p.az);

    LPInfo info;
    if (c.nx == 1) { info.setHiddenDirection(0); }
    else if (c.ny == 1) { info.setHiddenDirection(1); }
    if (c.semicoarsen && c.nx > 1 && c.ny > 1) {
        info.setSemicoarsening(true);
        info.setMaxSemicoarseningLevel(MLTerrainPoisson::coarsening_depth(p.geom.Domain(), -1).max());
        info.setSemicoarseningDirection(-1);
    }

    p.op = std::make_unique<MLTerrainPoisson>(Vector<Geometry>{p.geom}, Vector<BoxArray>{p.ba},
                                              Vector<DistributionMapping>{p.dm}, info,
                                              p.znd, p.ax, p.ay, p.az, p.dJ);
    Array<LinOpBCType,AMREX_SPACEDIM> bclo, bchi;
    get_terrain_projection_bc(p.geom, p.bc_names, false, bclo, bchi);
    p.op->setDomainBC(bclo, bchi);
    if (c.lshape) {
        p.op->setCoarseFineBC(nullptr, 2, LinOpBCType::Neumann);
    }
    p.op->setLevelBC(0, nullptr);
    p.op->setMaxOrder(2);
    return p;
}

// A smooth field with content in every direction; alternating adds a cell-to-cell
// oscillation (what a smoother must remove), zero_mean removes the plain mean
void
fill_field (MultiFab& mf, const Geometry& geom, bool alternating, Real shift = Real(0.0))
{
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto phi_ = geom.ProbHiArray();
    const Real Lx = phi_[0]-plo[0], Ly = phi_[1]-plo[1], Lz = phi_[2]-plo[2];
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> const& f = mf.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = (i+Real(0.5))*dx[0]/Lx;
            Real y = (j+Real(0.5))*dx[1]/Ly;
            Real z = (k+Real(0.5))*dx[2]/Lz;
            Real v = std::cos(Real(2.0)*PI*x + shift) * std::cos(Real(1.3)*PI*z)
                   + Real(0.5)*std::sin(Real(4.0)*PI*x - Real(0.7)) * std::sin(Real(2.0)*PI*y + Real(0.4)) * std::cos(Real(3.1)*z)
                   + Real(0.25)*std::cos(Real(6.0)*PI*x + Real(2.0)*PI*y + Real(1.0)) * z*z;
            if (alternating) {
                v += (((i+j+k) & 1) == 0) ? Real(0.5) : Real(-0.5);
            }
            f(i,j,k) = v;
        });
    }
}

void
subtract_mean (MultiFab& mf)
{
    Real mean = mf.sum(0, false) / static_cast<Real>(mf.boxArray().numPts());
    mf.plus(-mean, 0, 1, 0);
}

// y = L(in) through MLMG, which fills the ghost cells with the operator's applyBC
void
apply_operator (MLTerrainPoisson& op, MultiFab& in, MultiFab& out)
{
    MLMG mlmg(op);
    mlmg.setVerbose(0);
    mlmg.apply({&out}, {&in});
}

// Divergence of the face fluxes with the metric areas, as the operator defines it
void
flux_divergence (const Problem& p, const Array<MultiFab,AMREX_SPACEDIM>& flux, MultiFab& div)
{
    const auto dxinv = p.geom.InvCellSizeArray();
    for (MFIter mfi(div); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real      > const& d  = div.array(mfi);
        Array4<Real const> const& fx = flux[0].const_array(mfi);
        Array4<Real const> const& fy = flux[1].const_array(mfi);
        Array4<Real const> const& fz = flux[2].const_array(mfi);
        Array4<Real const> const& ax = p.ax.const_array(mfi);
        Array4<Real const> const& ay = p.ay.const_array(mfi);
        Array4<Real const> const& az = p.az.const_array(mfi);
        Array4<Real const> const& dJ = p.dJ.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // the kernels return minus the gradient
            Real v = (ax(i+1,j,k)*(-fx(i+1,j,k)) - ax(i,j,k)*(-fx(i,j,k))) * dxinv[0]
                   + (ay(i,j+1,k)*(-fy(i,j+1,k)) - ay(i,j,k)*(-fy(i,j,k))) * dxinv[1]
                   + (az(i,j,k+1)*(-fz(i,j,k+1)) - az(i,j,k)*(-fz(i,j,k))) * dxinv[2];
            d(i,j,k) = v / dJ(i,j,k);
        });
    }
}

// Manufactured solution phi = cos(2 pi x/Lx) cos(pi z/H) at the physical cell centres
// and its Laplacian
void
fill_manufactured (const Problem& p, MultiFab& phi, MultiFab& lap)
{
    const auto dx = p.geom.CellSizeArray();
    const Real Lx = p.geom.ProbHi(0) - p.geom.ProbLo(0);
    const Real H  = p.geom.ProbHi(2) - p.geom.ProbLo(2);
    const Real kx = Real(2.0)*PI/Lx;
    const Real kz = PI/H;
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real      > const& f = phi.array(mfi);
        Array4<Real      > const& l = lap.array(mfi);
        Array4<Real const> const& z = p.znd.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x  = (i+Real(0.5))*dx[0];
            Real zc = Real(0.125)*( z(i,j,k  ) + z(i+1,j,k  ) + z(i,j+1,k  ) + z(i+1,j+1,k  )
                                   +z(i,j,k+1) + z(i+1,j,k+1) + z(i,j+1,k+1) + z(i+1,j+1,k+1) );
            Real v = std::cos(kx*x) * std::cos(kz*zc);
            f(i,j,k) = v;
            l(i,j,k) = -(kx*kx + kz*kz) * v;
        });
    }
}

// Zero the layers within nlayers of the bottom and top of the domain
void
zero_vertical_layers (MultiFab& mf, const Geometry& geom, int nlayers)
{
    const int klo = geom.Domain().smallEnd(2) + nlayers;
    const int khi = geom.Domain().bigEnd(2)   - nlayers;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> const& f = mf.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k < klo || k > khi) { f(i,j,k) = Real(0.0); }
        });
    }
}

// The expected coarse field: every ratio-th value of fine, on the layout of crse
void
expected_subsample (const MultiFab& fine, const MultiFab& crse, const IntVect& ratio, MultiFab& expected)
{
    BoxArray fba_c = amrex::coarsen(fine.boxArray(), ratio);
    MultiFab tmp(fba_c, fine.DistributionMap(), 1, 0);
    const Dim3 r = ratio.dim3();
    for (MFIter mfi(tmp); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real const> const& f = fine.const_array(mfi);
        Array4<Real      > const& c = tmp.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            c(i,j,k) = f(i*r.x, j*r.y, k*r.z);
        });
    }
    expected.define(crse.boxArray(), crse.DistributionMap(), 1, 0);
    expected.ParallelCopy(tmp, 0, 0, 1, IntVect(0), IntVect(0));
}

// Largest violation of the ghost-node rules of a coarse nodal terrain field: lateral
// ghost nodes beyond a non-periodic domain face copy the face node, ghost nodes below
// and above the domain are linear extrapolations
Real
ghost_rule_violation (const MultiFab& znd, const Geometry& geom)
{
    const Box ndom = amrex::surroundingNodes(geom.Domain());
    const auto dlo = lbound(ndom);
    const auto dhi = ubound(ndom);
    const int per_x = geom.isPeriodic(0) ? 1 : 0;
    const int per_y = geom.isPeriodic(1) ? 1 : 0;
    MultiFab viol(znd.boxArray(), znd.DistributionMap(), 1, 1);
    viol.setVal(0.0);
    for (MFIter mfi(znd); mfi.isValid(); ++mfi) {
        const Box gbx = amrex::grow(mfi.validbox(), 1);
        Array4<Real const> const& z = znd.const_array(mfi);
        Array4<Real      > const& v = viol.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real d = Real(0.0);
            if (k >= dlo.z && k <= dhi.z) {
                if (!per_x && i == dlo.x-1) { d = amrex::max(d, amrex::Math::abs(z(i,j,k) - z(dlo.x,j,k))); }
                if (!per_x && i == dhi.x+1) { d = amrex::max(d, amrex::Math::abs(z(i,j,k) - z(dhi.x,j,k))); }
                if (!per_y && j == dlo.y-1 && i >= dlo.x && i <= dhi.x) { d = amrex::max(d, amrex::Math::abs(z(i,j,k) - z(i,dlo.y,k))); }
                if (!per_y && j == dhi.y+1 && i >= dlo.x && i <= dhi.x) { d = amrex::max(d, amrex::Math::abs(z(i,j,k) - z(i,dhi.y,k))); }
            } else if (k == dlo.z-1) {
                // the same arithmetic as ERF: 2 z0 - z1 below, one slope step above
                d = amrex::Math::abs(z(i,j,k) - (Real(2.0)*z(i,j,dlo.z) - z(i,j,dlo.z+1)));
            } else if (k == dhi.z+1) {
                d = amrex::Math::abs(z(i,j,k) - (z(i,j,dhi.z) + static_cast<Real>(k-dhi.z) * (z(i,j,dhi.z) - z(i,j,dhi.z-1))));
            }
            v(i,j,k) = d;
        });
    }
    return viol.norm0(0, 1);
}

// Largest deviation of the level's metrics from the formulas evaluated on its own surface
Real
metric_violation (const MLTerrainPoisson& op, int mglev, const Geometry& geom)
{
    const Real dzinv = Real(1.0)/geom.CellSize(2);
    const MultiFab& znd = op.zPhysNd(mglev);
    const MultiFab& ax  = op.axCoef(mglev);
    const MultiFab& ay  = op.ayCoef(mglev);
    const MultiFab& az  = op.azCoef(mglev);
    const MultiFab& dJ  = op.detJ(mglev);
    MultiFab viol(dJ.boxArray(), dJ.DistributionMap(), 1, 0);
    for (MFIter mfi(viol); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real const> const& z  = znd.const_array(mfi);
        Array4<Real const> const& xa = ax.const_array(mfi);
        Array4<Real const> const& ya = ay.const_array(mfi);
        Array4<Real const> const& za = az.const_array(mfi);
        Array4<Real const> const& J  = dJ.const_array(mfi);
        Array4<Real      > const& v  = viol.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real dJe = Real(0.25)*dzinv*( z(i,j,k+1)+z(i+1,j,k+1)+z(i,j+1,k+1)+z(i+1,j+1,k+1)
                                         -z(i,j,k  )-z(i+1,j,k  )-z(i,j+1,k  )-z(i+1,j+1,k  ) );
            Real axe = Real(0.5)*dzinv*( z(i,j,k+1)+z(i,j+1,k+1)-z(i,j,k)-z(i,j+1,k) );
            Real aye = Real(0.5)*dzinv*( z(i,j,k+1)+z(i+1,j,k+1)-z(i,j,k)-z(i+1,j,k) );
            Real d = amrex::Math::abs(J(i,j,k) - dJe);
            d = amrex::max(d, amrex::Math::abs(xa(i,j,k) - axe));
            d = amrex::max(d, amrex::Math::abs(ya(i,j,k) - aye));
            d = amrex::max(d, amrex::Math::abs(za(i,j,k) - Real(1.0)));
            v(i,j,k) = d;
        });
    }
    return viol.norm0(0, 0);
}

// Pattern 1 where (i mod 3, j mod 3, k mod 3) == (a,b,c), else 0
void
fill_pattern (MultiFab& mf, int a, int b, int c)
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> const& f = mf.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ia = ((i % 3) + 3) % 3;
            int jb = ((j % 3) + 3) % 3;
            int kc = ((k % 3) + 3) % 3;
            f(i,j,k) = (ia == a && jb == b && kc == c) ? Real(1.0) : Real(0.0);
        });
    }
}

// Scatter L(pattern(a,b,c)) into the column coefficient it measures: in row (i,j,k)
// with i = a, j = b (mod 3) the pattern is the unit vector of the cell k+d with k+d = c (mod 3)
void
scatter_pattern_result (const MultiFab& out, int a, int b, int c, MultiFab& tri_expected)
{
    for (MFIter mfi(out); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real const> const& o = out.const_array(mfi);
        Array4<Real      > const& t = tri_expected.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ia = ((i % 3) + 3) % 3;
            int jb = ((j % 3) + 3) % 3;
            if (ia == a && jb == b) {
                int d = ((c - k) % 3 + 3) % 3;   // 0, 1 or 2 (= -1)
                int comp = (d == 2) ? 0 : d + 1;
                t(i,j,k,comp) = o(i,j,k);
            }
        });
    }
}

Real
residual_norm (MLTerrainPoisson& op, MultiFab& sol, const MultiFab& rhs)
{
    MultiFab r(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
    apply_operator(op, sol, r);
    MultiFab::Subtract(r, rhs, 0, 0, 1, 0);
    return r.norm2();
}

// Tolerances that the precision of amrex::Real can resolve: machine epsilon and the
// tightest solver tolerance worth asking for (1e-10 in double, 1e-5 in single)
constexpr Real eps = std::numeric_limits<Real>::epsilon();
constexpr Real solve_tol = (sizeof(Real) == 8) ? Real(1.e-10) : Real(1.e-5);
constexpr Real solve_err = (sizeof(Real) == 8) ? Real(1.e-7)  : Real(1.e-3);

// Scientific notation for the numbers recorded in the test report
std::string
sci (Real v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", static_cast<double>(v));
    return std::string(buf);
}

HillCase
periodic_3d_case ()
{
    HillCase c;
    c.nx = 32; c.ny = 16; c.nz = 16;
    c.Lx = 3200.; c.Ly = 1600.; c.H = 1000.;
    c.per_x = true; c.per_y = true;
    c.hill = 0; c.hmax = 300.; c.a = 400.;
    c.bc_names = {"Periodic","Periodic","surface_layer","Periodic","Periodic","SlipWall"};
    c.max_grid_size = IntVect(16,16,8);   // split in x and in z
    return c;
}

HillCase
ridge_2d_case ()
{
    HillCase c;
    c.nx = 96; c.ny = 1; c.nz = 48;
    c.Lx = 2400.; c.Ly = 25.; c.H = 1200.;
    c.per_x = false; c.per_y = true;
    c.hill = 1; c.hmax = 300.; c.a = 310.;   // maximum slope about 0.63, as the Ishihara ridge
    c.bc_names = {"Inflow","Periodic","surface_layer","Outflow","Periodic","SlipWall"};
    c.max_grid_size = IntVect(48,1,48);
    return c;
}

// Few vertical cells: z allows two coarsenings, x five and y four, so the hierarchy
// continues laterally with ratio (2,2,1) once z has stopped
HillCase
shallow_3d_case ()
{
    HillCase c = periodic_3d_case();
    c.nz = 8; c.H = 400.;
    c.max_grid_size = IntVect(16,16,8);
    c.semicoarsen = true;
    return c;
}

// A refined level inside a periodic domain whose boxes form an L: all its faces are
// coarse/fine (Neumann) faces, so the problem is singular
HillCase
lshape_case ()
{
    HillCase c = periodic_3d_case();
    c.nx = 32; c.ny = 32; c.nz = 16;
    c.Lx = 3200.; c.Ly = 3200.; c.H = 800.;
    c.max_grid_size = IntVect(8,8,16);
    c.lshape = true;
    return c;
}

HillCase
walls_3d_case ()
{
    HillCase c;
    c.nx = 12; c.ny = 6; c.nz = 12;
    c.Lx = 1200.; c.Ly = 600.; c.H = 800.;
    c.per_x = false; c.per_y = false;
    c.hill = 1; c.hmax = 200.; c.a = 250.;
    c.bc_names = {"Inflow","SlipWall","surface_layer","Outflow","SlipWall","SlipWall"};
    c.max_grid_size = IntVect(6,6,12);    // split in x and y; ERF never splits its grids in z
    return c;
}

} // namespace

//
// The multigrid operator applies exactly the stencil of the GMRES operator
//
TEST(MLTerrainPoisson, MatchesGMRESOperatorBitwise)
{
    for (int icase = 0; icase < 3; ++icase)
    {
        SCOPED_TRACE("case " + std::to_string(icase));
        HillCase c = (icase == 0) ? periodic_3d_case() : ((icase == 1) ? ridge_2d_case() : walls_3d_case());
        Problem p = build_problem(c);

        MultiFab in_mg(p.ba, p.dm, 1, 1);
        fill_field(in_mg, p.geom, false);
        MultiFab in_gm(p.ba, p.dm, 1, 1);
        MultiFab::Copy(in_gm, in_mg, 0, 0, 1, 1);

        MultiFab out_mg(p.ba, p.dm, 1, 0);
        apply_operator(*p.op, in_mg, out_mg);

        Gpu::DeviceVector<Real> dz_d;
        TerrainPoisson tp(p.geom, p.geom, p.ba, p.dm, p.bc_names, dz_d,
                          p.ax, p.ay, p.az, p.dJ, &p.znd, false, /*build_fft_precond=*/false);
        MultiFab out_gm(p.ba, p.dm, 1, 0);
        tp.apply(out_gm, in_gm);

        const Real scale = out_gm.norm0();
        EXPECT_GT(scale, Real(0.0)) << "the test field must not be in the null space";
        MultiFab::Subtract(out_gm, out_mg, 0, 0, 1, 0);
        EXPECT_EQ(out_gm.norm0(), Real(0.0)) << "operator result differs from TerrainPoisson::apply";
    }
}

//
// The operator does not depend on how the level is cut into boxes: with boxes split
// in z as well, the ghost cells outside a box in two directions lie beyond a domain
// face in one of them and an interior face in the other, and must be mirrored
// across the domain face only.  (The GMRES operator leaves those cells unfilled,
// which ERF's rule of never splitting grids in z hides.)
//
TEST(MLTerrainPoisson, OperatorIndependentOfBoxDecomposition)
{
    HillCase c = walls_3d_case();
    MultiFab out_ref;
    for (int icase = 0; icase < 3; ++icase)
    {
        SCOPED_TRACE("decomposition " + std::to_string(icase));
        c.max_grid_size = (icase == 0) ? IntVect(12,6,12) : ((icase == 1) ? IntVect(6,6,12) : IntVect(6,3,6));
        Problem p = build_problem(c);
        MultiFab in(p.ba, p.dm, 1, 1);
        fill_field(in, p.geom, false);
        MultiFab out(p.ba, p.dm, 1, 0);
        apply_operator(*p.op, in, out);
        if (icase == 0) {
            out_ref.define(p.ba, p.dm, 1, 0);
            MultiFab::Copy(out_ref, out, 0, 0, 1, 0);
            EXPECT_GT(out_ref.norm0(), Real(0.0));
        } else {
            // bring the result onto the single-box layout and compare
            MultiFab tmp(out_ref.boxArray(), out_ref.DistributionMap(), 1, 0);
            tmp.ParallelCopy(out, 0, 0, 1, IntVect(0), IntVect(0));
            MultiFab::Subtract(tmp, out_ref, 0, 0, 1, 0);
            EXPECT_EQ(tmp.norm0(), Real(0.0)) << "operator depends on the box decomposition";
        }
    }
}

//
// The fluxes the velocity correction uses are the ones the operator is the divergence of
//
TEST(MLTerrainPoisson, FluxDivergenceIsTheOperator)
{
    for (int icase = 0; icase < 2; ++icase)
    {
        SCOPED_TRACE("case " + std::to_string(icase));
        HillCase c = (icase == 0) ? periodic_3d_case() : ridge_2d_case();
        Problem p = build_problem(c);

        MultiFab in(p.ba, p.dm, 1, 1);
        fill_field(in, p.geom, false);
        MultiFab out(p.ba, p.dm, 1, 0);
        apply_operator(*p.op, in, out);

        Array<MultiFab,AMREX_SPACEDIM> flux;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            flux[idim].define(amrex::convert(p.ba, IntVect::TheDimensionVector(idim)), p.dm, 1, 0);
        }
        MLMG mlmg(*p.op);
        mlmg.setVerbose(0);
        mlmg.getFluxes({GetArrOfPtrs(flux)}, {&in}, MLMG::Location::FaceCenter);

        MultiFab div(p.ba, p.dm, 1, 0);
        flux_divergence(p, flux, div);

        const Real scale = out.norm0();
        MultiFab::Subtract(div, out, 0, 0, 1, 0);
        // the two sum the same face terms in a different order
        EXPECT_LE(div.norm0(), Real(1.e4)*eps*scale) << "flux divergence differs from the operator";
    }
}

//
// Second-order consistency with the Laplacian on a hill, cross terms included
//
TEST(MLTerrainPoisson, ManufacturedSolutionConvergesAtSecondOrder)
{
    Real err[3];
    for (int ires = 0; ires < 3; ++ires)
    {
        HillCase c;
        c.nx = 32 << ires; c.ny = 1; c.nz = 16 << ires;
        c.Lx = 2000.; c.Ly = 20.; c.H = 1000.;
        c.per_x = true; c.per_y = true;
        c.hill = 0; c.hmax = 200.; c.a = 0.;      // maximum slope pi*hmax/Lx = 0.31
        c.bc_names = {"Periodic","Periodic","surface_layer","Periodic","Periodic","SlipWall"};
        c.max_grid_size = IntVect(32,1,64);
        Problem p = build_problem(c);

        MultiFab phi(p.ba, p.dm, 1, 1);
        MultiFab lap(p.ba, p.dm, 1, 0);
        fill_manufactured(p, phi, lap);

        MultiFab out(p.ba, p.dm, 1, 0);
        apply_operator(*p.op, phi, out);

        MultiFab::Subtract(out, lap, 0, 0, 1, 0);
        // The manufactured phi does not satisfy the zero-flux condition on the sloping
        // surface, so measure the truncation error away from the bottom and top
        zero_vertical_layers(out, p.geom, 2);
        err[ires] = out.norm0();
        RecordProperty("err_" + std::to_string(ires), sci(err[ires]));
    }
    EXPECT_GT(err[0], Real(0.0));
    EXPECT_GE(err[0]/err[1], Real(3.5)) << err[0] << " " << err[1];
    EXPECT_GE(err[1]/err[2], Real(3.5)) << err[1] << " " << err[2];
}

//
// Coarse levels keep every other node of the fine surface and recompute the metrics from it
//
TEST(MLTerrainPoisson, CoarseLevelsSubsampleTheSurface)
{
    for (int icase = 0; icase < 3; ++icase)
    {
        SCOPED_TRACE("case " + std::to_string(icase));
        HillCase c = (icase == 0) ? walls_3d_case() : ((icase == 1) ? ridge_2d_case() : shallow_3d_case());
        if (icase == 0) { c.nx = 32; c.ny = 16; c.nz = 16; c.max_grid_size = IntVect(16,16,8); }
        Problem p = build_problem(c);
        const MLTerrainPoisson& op = *p.op;

        ASSERT_GE(op.numMGLevels(), 3) << "the test needs a real hierarchy";
        if (icase == 2) {
            // z stops after two coarsenings; the lateral directions go on
            ASSERT_GE(op.numMGLevels(), 4) << "semicoarsening did not deepen the hierarchy";
            EXPECT_EQ(op.Geom(0, op.numMGLevels()-1).Domain().length(2), 2);
        }

        for (int mglev = 1; mglev < op.numMGLevels(); ++mglev)
        {
            SCOPED_TRACE("mglev " + std::to_string(mglev));
            const Geometry& gc = op.Geom(0, mglev);
            const Geometry& gf = op.Geom(0, mglev-1);
            const IntVect ratio = gf.Domain().length() / gc.Domain().length();
            if (c.ny == 1) { EXPECT_EQ(ratio[1], 1); }
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { EXPECT_TRUE(ratio[d] == 1 || ratio[d] == 2); }

            MultiFab expected;
            expected_subsample(op.zPhysNd(mglev-1), op.zPhysNd(mglev), ratio, expected);
            MultiFab::Subtract(expected, op.zPhysNd(mglev), 0, 0, 1, 0);
            EXPECT_EQ(expected.norm0(0, 0), Real(0.0)) << "coarse nodes are not the fine nodes";

            EXPECT_EQ(ghost_rule_violation(op.zPhysNd(mglev), gc), Real(0.0));
            EXPECT_LE(metric_violation(op, mglev, gc), Real(100.0)*eps);
        }
    }
}

//
// The smoother's tridiagonal coefficients are the operator's own column entries,
// boundary folds included
//
TEST(MLTerrainPoisson, TridiagonalMatchesTheOperator)
{
    for (int icase = 0; icase < 2; ++icase)
    {
        SCOPED_TRACE("case " + std::to_string(icase));
        HillCase c;
        if (icase == 0) {
            // periodic: the domain lengths are multiples of 3 so the patterns are periodic
            c = periodic_3d_case();
            c.nx = 24; c.ny = 6; c.nz = 12;
            c.max_grid_size = IntVect(12,6,6);
        } else {
            c = walls_3d_case();
            c.max_grid_size = IntVect(6,6,6);   // interior faces in every direction
        }
        Problem p = build_problem(c);

        MultiFab tri_expected(p.ba, p.dm, 3, 0);
        tri_expected.setVal(0.0);
        MultiFab in(p.ba, p.dm, 1, 1);
        MultiFab out(p.ba, p.dm, 1, 0);
        for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
        for (int cc = 0; cc < 3; ++cc) {
            fill_pattern(in, a, b, cc);
            apply_operator(*p.op, in, out);
            scatter_pattern_result(out, a, b, cc, tri_expected);
        }}}

        const MultiFab& tri = p.op->triDiag(0);
        ASSERT_EQ(tri.nComp(), 3);
        const Real scale = tri.norm0(1, 0);
        EXPECT_GT(scale, Real(0.0));
        MultiFab::Subtract(tri_expected, tri, 0, 0, 3, 0);
        for (int comp = 0; comp < 3; ++comp) {
            EXPECT_EQ(tri_expected.norm0(comp, 0), Real(0.0)) << "component " << comp;
        }
    }
}

//
// The column relaxation reduces the residual sweep after sweep
//
TEST(MLTerrainPoisson, SmootherReducesTheResidual)
{
    for (int icase = 0; icase < 3; ++icase)
    {
        SCOPED_TRACE("case " + std::to_string(icase));
        HillCase c = (icase == 0) ? periodic_3d_case() : ((icase == 1) ? ridge_2d_case() : lshape_case());
        Problem p = build_problem(c);

        MultiFab sol(p.ba, p.dm, 1, 1);
        fill_field(sol, p.geom, true);
        MultiFab rhs(p.ba, p.dm, 1, 0);
        rhs.setVal(0.0);

        // the first apply prepares the operator (tridiagonal coefficients included)
        Real r_prev = residual_norm(*p.op, sol, rhs);
        const Real r0 = r_prev;
        for (int sweep = 0; sweep < 5; ++sweep) {
            p.op->smooth(0, 0, sol, rhs, false, 1);
            Real r = residual_norm(*p.op, sol, rhs);
            EXPECT_LT(r, r_prev) << "sweep " << sweep;
            r_prev = r;
        }
        RecordProperty("reduction_" + std::to_string(icase), sci(r_prev/r0));
        EXPECT_LT(r_prev, Real(0.2)*r0) << "five sweeps reduced the residual only to " << r_prev/r0;
    }
}

//
// Multigrid converges as a solver on a singular (all Neumann/periodic) problem and on a
// steep ridge with a Dirichlet outflow, to the discrete solution the right-hand side came from
//
TEST(MLTerrainPoisson, SolveConverges)
{
    for (int icase = 0; icase < 4; ++icase)
    {
        SCOPED_TRACE("case " + std::to_string(icase));
        HillCase c = (icase == 0) ? periodic_3d_case() : ((icase == 1) ? ridge_2d_case()
                   : ((icase == 2) ? shallow_3d_case() : lshape_case()));
        Problem p = build_problem(c);
        const bool singular = (icase != 1);
        if (icase == 2) { ASSERT_GE(p.op->numMGLevels(), 4); }
        if (icase == 3) { ASSERT_GE(p.op->numMGLevels(), 3); }

        MultiFab phi_exact(p.ba, p.dm, 1, 1);
        fill_field(phi_exact, p.geom, false, Real(0.3));
        if (singular) { subtract_mean(phi_exact); }
        MultiFab rhs(p.ba, p.dm, 1, 0);
        apply_operator(*p.op, phi_exact, rhs);

        MultiFab phi(p.ba, p.dm, 1, 1);
        phi.setVal(0.0);

        MLMG mlmg(*p.op);
        mlmg.setVerbose(0);
        mlmg.setMaxIter(50);
        mlmg.setThrowException(true);
        EXPECT_EQ(p.op->isSingular(0), singular);
        try {
            mlmg.solve({&phi}, {&rhs}, solve_tol, Real(0.0));
        } catch (std::runtime_error const& e) {
            ADD_FAILURE() << e.what();
            continue;
        }
        RecordProperty("iters_" + std::to_string(icase), mlmg.getNumIters());
        EXPECT_LE(mlmg.getNumIters(), 30);

        MultiFab::Subtract(phi, phi_exact, 0, 0, 1, 0);
        if (singular) { subtract_mean(phi); }
        const Real scale = phi_exact.norm0();
        EXPECT_LE(phi.norm0(), solve_err*scale) << "solution error " << phi.norm0()/scale;
    }
}


//
// The driver re-grids a solve whose boxes would stop the coarsening early; the
// helpers must report the depth a layout allows and build boxes that allow the
// domain's depth while covering the domain exactly
//
TEST(MLTerrainPoisson, MultigridFriendlyGrids)
{
    // The ridge deck: 600 x 1 x 224 cells in boxes 150 wide (150 = 2 * 75, 224 = 32 * 7)
    const Box ridge(IntVect(0,0,0), IntVect(599,0,223));
    EXPECT_EQ(MLTerrainPoisson::coarsening_depth(ridge, 1), IntVect(3,0,5));
    BoxArray deck(ridge);
    deck.maxSize(IntVect(150,1,224));
    EXPECT_EQ(deck.size(), 4);
    EXPECT_EQ(MLTerrainPoisson::coarsening_depth(deck, 1), IntVect(1,0,5));

    const IntVect ridge_depth = MLTerrainPoisson::coarsening_depth(ridge, 1);
    BoxArray mg = MLTerrainPoisson::multigrid_grids(ridge, 1, ridge_depth, 64);
    EXPECT_EQ(mg.minimalBox(), ridge);
    EXPECT_EQ(mg.numPts(), ridge.numPts());
    EXPECT_EQ(MLTerrainPoisson::coarsening_depth(mg, 1), ridge_depth);
    EXPECT_GE(mg.size(), 10);
    for (int i = 0; i < mg.size(); ++i) {
        EXPECT_LE(mg[i].length(0), 64 + 8);
        EXPECT_EQ(mg[i].length(1), 1);
    }

    // Askervein: 300 x 300 x 18 cells allows two coarsenings laterally and one in z;
    // 32-cell boxes with a 12-cell remainder allow two laterally too, so no re-gridding
    const Box ask(IntVect(0,0,0), IntVect(299,299,17));
    EXPECT_EQ(MLTerrainPoisson::coarsening_depth(ask, -1), IntVect(2,2,1));
    BoxArray ask_ba(ask);
    ask_ba.maxSize(32);
    EXPECT_TRUE(MLTerrainPoisson::coarsening_depth(ask_ba, -1).allGE(IntVect(2,2,1)));

    // Powers of two coarsen all the way to the minimum width
    const Box pw(IntVect(0,0,0), IntVect(63,31,15));
    EXPECT_EQ(MLTerrainPoisson::coarsening_depth(pw, -1), IntVect(5,4,3));
    BoxArray pw_mg = MLTerrainPoisson::multigrid_grids(pw, -1, IntVect(5,4,3), 16);
    EXPECT_EQ(pw_mg.numPts(), pw.numPts());
    EXPECT_TRUE(MLTerrainPoisson::coarsening_depth(pw_mg, -1).allGE(IntVect(5,4,3)));
}


//
// The existing GMRES driver (amrex::GMRES on TerrainPoisson) accepts the multigrid
// V-cycle as its preconditioner through TerrainPoisson::setPrecondFunction, which is
// the route a cheaper preconditioner can take later.  With the FFT preconditioner
// replaced by one V-cycle (bottom: smoothing sweeps, so the preconditioner is a fixed
// linear operator) GMRES must converge to the discrete solution in few iterations.
//
TEST(MLTerrainPoisson, PreconditionsTheGMRESDriver)
{
    HillCase c = ridge_2d_case();
    Problem p = build_problem(c);

    MultiFab phi_exact(p.ba, p.dm, 1, 1);
    fill_field(phi_exact, p.geom, false, Real(0.3));
    MultiFab rhs(p.ba, p.dm, 1, 0);
    apply_operator(*p.op, phi_exact, rhs);

    MLMG mlmg(*p.op);
    mlmg.setVerbose(0);
    mlmg.setBottomVerbose(0);
    mlmg.setBottomSolver(BottomSolver::smoother);
    mlmg.setPrecondIter(1);

    Gpu::DeviceVector<Real> dz_d;
    TerrainPoisson tp(p.geom, p.geom, p.ba, p.dm, p.bc_names, dz_d,
                      p.ax, p.ay, p.az, p.dJ, &p.znd, false, /*build_fft_precond=*/false);
    int ncalls = 0;
    tp.setPrecondFunction([&] (MultiFab& lhs, MultiFab const& r)
    {
        ++ncalls;
        lhs.setVal(0.0);
        mlmg.precond({&lhs}, {&r}, Real(0.0), Real(0.0));
    });
    tp.usePrecond(true);

    MultiFab phi(p.ba, p.dm, 1, 1);
    phi.setVal(0.0);
    amrex::GMRES<MultiFab, TerrainPoisson> gmres;
    gmres.define(tp);
    gmres.setVerbose(0);
    gmres.setRestartLength(50);
    gmres.setMaxIters(60);
    gmres.solve(phi, rhs, solve_tol, Real(0.0));

    EXPECT_GT(ncalls, 0) << "the preconditioner hook was never called";
    EXPECT_LE(gmres.getNumIters(), 40);
    RecordProperty("gmres_iters", gmres.getNumIters());

    // the true residual, not the one GMRES tracks
    MultiFab res(p.ba, p.dm, 1, 0);
    apply_operator(*p.op, phi, res);
    MultiFab::Subtract(res, rhs, 0, 0, 1, 0);
    EXPECT_LE(res.norm0(), Real(10.0)*solve_tol*rhs.norm0()) << "true residual " << res.norm0()/rhs.norm0();

    MultiFab::Subtract(phi, phi_exact, 0, 0, 1, 0);
    EXPECT_LE(phi.norm0(), solve_err*phi_exact.norm0());
}

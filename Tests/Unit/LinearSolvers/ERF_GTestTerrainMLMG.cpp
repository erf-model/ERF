#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

#include <AMReX_GMRES_MLMG.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLTerrainPoisson.H>

#include <ERF_TerrainPoisson_3D_K.H>
#include <ERF_Utils.H>

#include <gtest/gtest.h>

#include <cmath>

using namespace amrex;

// erf.terrain_poisson_solver = mlmg and gmres_mlmg hand the projection on a terrain-fitted
// mesh to amrex::MLTerrainPoisson (AMReX-Codes/amrex#6030, getFluxes #6068).  ERF passes it
// its own metrics (face areas, dJ, z_phys_nd) and takes its fluxes back as the momentum
// correction, so the two operators must be the same: the AMReX one applied to a field must
// be dJ times ERF's operator (zero flux through the Neumann faces; ERF needs the FFT build
// only for its preconditioner), and its fluxes must be ERF's fluxes, zero on
// the Neumann faces.  None of these tests needs FFT.

namespace {

Real tol (Real scale) { return scale * (sizeof(Real) == 8 ? Real(1.e-12) : Real(1.e-4)); }

constexpr int NX = 8, NY = 6, NZ = 6;
constexpr Real DX = Real(80.0), DY = Real(60.0), DZ = Real(40.0);
constexpr Real ZTOP = NZ * DZ;

// A hill off the centre on a tilted plane: the terrain slopes at every lateral face
AMREX_GPU_HOST_DEVICE inline Real terrain (Real x, Real y)
{
    const Real xr = (x - Real(250.0)) / Real(300.0);
    const Real yr = (y - Real(200.0)) / Real(250.0);
    return Real(60.0) * (Real(1.0) / (Real(1.0) + xr*xr + yr*yr) + Real(0.05) * x / DX + Real(0.02) * y / DY);
}

// Node heights of the hill (a free function: CUDA does not allow a device lambda in a constructor)
void fill_heights (MultiFab& zp)
{
    for (MFIter mfi(zp); mfi.isValid(); ++mfi) {
        auto z = zp.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real zeta = Real(k) * DZ;
            z(i,j,k) = zeta + terrain(Real(i)*DX, Real(j)*DY) * (Real(1.0) - zeta / ZTOP);
        });
    }
}

// The region with slip walls on every side, split into boxes, with the metrics ERF makes
struct Case
{
    Geometry geom;
    BoxArray ba;
    DistributionMapping dm;
    MultiFab zp, ax, ay, az, dJ;

    Case ()
    {
        const Box domain(IntVect(0), IntVect(NX-1, NY-1, NZ-1));
        RealBox rb({Real(0.0), Real(0.0), Real(0.0)}, {NX*DX, NY*DY, NZ*DZ});
        geom.define(domain, rb, 0, {0, 0, 0});
        ba.define(domain);
        ba.maxSize(IntVect(4, 3, NZ));
        dm.define(ba);
        zp.define(convert(ba, IntVect(1)), dm, 1, 2);
        ax.define(convert(ba, IntVect(1,0,0)), dm, 1, 1);
        ay.define(convert(ba, IntVect(0,1,0)), dm, 1, 1);
        az.define(convert(ba, IntVect(0,0,1)), dm, 1, 1);
        dJ.define(ba, dm, 1, 1);
        fill_heights(zp);
        make_J(geom, zp, dJ);
        make_areas(geom, zp, ax, ay, az);
    }

    void set_field (MultiFab& phi) const
    {
        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
            auto f = phi.array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                f(i,j,k) = std::sin(Real(0.7)*Real(i) + Real(0.3)) * std::cos(Real(0.5)*Real(j) - Real(0.2))
                         + Real(0.3) * std::cos(Real(0.9)*Real(k) + Real(0.1)*Real(i*j))
                         + Real(0.05) * Real(i + 2*j - k);
            });
        }
    }

    // The even mirror of TerrainPoisson::apply_bcs on every face (x, then y, then z)
    void mirror (MultiFab& phi) const
    {
        phi.FillBoundary(geom.periodicity());
        const Dim3 lo = lbound(geom.Domain()), hi = ubound(geom.Domain());
        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto f = phi.array(mfi);
            if (bx.smallEnd(0) == lo.x) { ParallelFor(makeSlab(bx,0,lo.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) { f(i-1,j,k) = f(i,j,k); }); }
            if (bx.bigEnd(0)   == hi.x) { ParallelFor(makeSlab(bx,0,hi.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) { f(i+1,j,k) = f(i,j,k); }); }
            Box bx2(bx); bx2.grow(0,1);
            if (bx.smallEnd(1) == lo.y) { ParallelFor(makeSlab(bx2,1,lo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) { f(i,j-1,k) = f(i,j,k); }); }
            if (bx.bigEnd(1)   == hi.y) { ParallelFor(makeSlab(bx2,1,hi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) { f(i,j+1,k) = f(i,j,k); }); }
            Box gbx(bx); gbx.grow(0,1); gbx.grow(1,1);
            if (bx.smallEnd(2) == lo.z) { ParallelFor(makeSlab(gbx,2,lo.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) { f(i,j,k-1) = f(i,j,k); }); }
            if (bx.bigEnd(2)   == hi.z) { ParallelFor(makeSlab(gbx,2,hi.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) { f(i,j,k+1) = f(i,j,k); }); }
        }
        phi.FillBoundary(geom.periodicity());
    }

    // ERF's zero-flux operator (all Neumann faces of the region; the ground is zero anyway)
    void erf_apply (MultiFab& y, MultiFab& x) const
    {
        mirror(x);
        const Dim3 lo = lbound(geom.Domain()), hi = ubound(geom.Domain());
        const int zf = tp_face_xlo | tp_face_xhi | tp_face_ylo | tp_face_yhi | tp_face_zhi;
        const auto dxinv = geom.InvCellSizeArray();
        auto const& ya = y.arrays();
        auto const& xa = x.const_arrays();
        auto const& zpa = zp.const_arrays(); auto const& axa = ax.const_arrays();
        auto const& aya = ay.const_arrays(); auto const& aza = az.const_arrays();
        auto const& J = dJ.const_arrays();
        ParallelFor(y, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
        {
            terrpoisson_adotx(i, j, k, ya[b], xa[b], axa[b], aya[b], aza[b], J[b], zpa[b],
                              dxinv[0], dxinv[1], dxinv[2], lo, hi, zf);
        });
        Gpu::streamSynchronize();
    }

    void make_linop (MLTerrainPoisson& linop) const
    {
        linop.define({geom}, {ba}, {dm});
        linop.setDomainBC({LinOpBCType::Neumann, LinOpBCType::Neumann, LinOpBCType::Neumann},
                          {LinOpBCType::Neumann, LinOpBCType::Neumann, LinOpBCType::Neumann});
        linop.setLevelBC(0, nullptr);
        linop.setZPhys(0, zp);
        linop.setAreas(0, {&ax, &ay, &az});
        linop.setDetJ(0, dJ);
    }
};

Real max_abs_diff (MultiFab const& a, MultiFab const& b)
{
    MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    MultiFab::LinComb(d, Real(1.0), a, 0, Real(-1.0), b, 0, 0, 1, 0);
    return d.norm0();
}

} // namespace

// amrex::MLTerrainPoisson applies dJ times ERF's zero-flux operator: the same stencil, the
// same metrics, the same zero flux through the Neumann faces.  Measured: they agree to
// 2.8e-16 of the largest value; ERF's old mirror-only operator differs from it by 3.8e-2,
// the cross-term flux through the faces of this hill.
TEST(TerrainMLMG, ApplyIsDetJTimesTheZeroFluxOperator)
{
    Case c;
    MultiFab x(c.ba, c.dm, 1, 1), y_erf(c.ba, c.dm, 1, 0), y_amrex(c.ba, c.dm, 1, 0);
    c.set_field(x);
    c.erf_apply(y_erf, x);
    MultiFab::Multiply(y_erf, c.dJ, 0, 0, 1, 0);

    MLTerrainPoisson linop;
    c.make_linop(linop);
    MLMG mlmg(linop);
    mlmg.apply({&y_amrex}, {&x});

    const Real scale = y_erf.norm0();
    EXPECT_GT(scale, Real(1.e-5));
    EXPECT_LE(max_abs_diff(y_erf, y_amrex), tol(Real(10.0) * scale));
}

// The fluxes ERF takes for the momentum correction (mlmg.getFluxes with the solution) are
// terrpoisson_flux_x/y/z, the fluxes TerrainPoisson::getFluxes gives, inside the region, and
// zero on its Neumann faces: same sign and scale.
TEST(TerrainMLMG, FluxesAreTerrainPoissonFluxes)
{
    Case c;
    MultiFab x(c.ba, c.dm, 1, 1);
    c.set_field(x);
    c.mirror(x);

    Array<MultiFab,AMREX_SPACEDIM> f;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        f[d].define(convert(c.ba, IntVect::TheDimensionVector(d)), c.dm, 1, 0);
    }
    MLTerrainPoisson linop;
    c.make_linop(linop);
    MLMG mlmg(linop);
    // The metric terms are built when the operator is first applied or solved with, as ERF
    // does before it takes the fluxes
    MultiFab tmp(c.ba, c.dm, 1, 0);
    mlmg.apply({&tmp}, {&x});
    mlmg.getFluxes({{&f[0], &f[1], &f[2]}}, {&x});

    const auto dxinv = c.geom.InvCellSizeArray();
    const int top[3] = {NX, NY, NZ};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real on_face = 0, inside = 0, scale = 0;
        for (MFIter mfi(f[d]); mfi.isValid(); ++mfi) {
            const Box& b = mfi.validbox();
            FArrayBox fh(b, 1, The_Pinned_Arena()), xh(x[mfi].box(), 1, The_Pinned_Arena());
            FArrayBox zh(c.zp[mfi].box(), 1, The_Pinned_Arena());
            fh.copy<RunOn::Device>(f[d][mfi]);
            xh.copy<RunOn::Device>(x[mfi]);
            zh.copy<RunOn::Device>(c.zp[mfi]);
            Gpu::streamSynchronize();
            auto fa = fh.const_array(); auto xa = xh.const_array(); auto za = zh.const_array();
            for (int k = b.smallEnd(2); k <= b.bigEnd(2); ++k) {
            for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j) {
            for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i) {
                const int n = (d == 0) ? i : ((d == 1) ? j : k);
                if (n == 0 || n == top[d]) {
                    on_face = amrex::max(on_face, std::abs(fa(i,j,k)));
                    continue;
                }
                const Real ref = (d == 0) ? terrpoisson_flux_x(i,j,k,xa,za,dxinv[0])
                               : ((d == 1) ? terrpoisson_flux_y(i,j,k,xa,za,dxinv[1])
                                           : terrpoisson_flux_z(i,j,k,xa,za,dxinv[0],dxinv[1]));
                scale = amrex::max(scale, std::abs(ref));
                inside = amrex::max(inside, std::abs(fa(i,j,k) - ref));
            }}}
        }
        EXPECT_EQ(on_face, Real(0.0)) << "direction " << d;
        EXPECT_GT(scale, Real(1.e-4)) << "direction " << d;
        EXPECT_LE(inside, tol(Real(10.0) * scale)) << "direction " << d;
    }
}

// The singular problem of a refined level or a walled domain, right-hand side with zero
// dJ-weighted mean: MLMG and GMRESMLMG both converge, and the solution satisfies ERF's own
// zero-flux operator to the same tolerance.
TEST(TerrainMLMG, SolvesTheSingularProblem)
{
    Case c;
    MultiFab rhs(c.ba, c.dm, 1, 0);
    c.set_field(rhs);
    const Real mean = MultiFab::Dot(c.dJ, 0, rhs, 0, 1, 0) / c.dJ.sum(0);
    rhs.plus(-mean, 0, 1, 0);
    const Real rel = (sizeof(Real) == 8) ? Real(1.e-10) : Real(1.e-5);

    for (bool use_gmres : {false, true}) {
        MLTerrainPoisson linop;
        c.make_linop(linop);
        MLMG mlmg(linop);
        mlmg.setMaxIter(100);
        MultiFab phi(c.ba, c.dm, 1, 1);
        phi.setVal(0.0);
        if (use_gmres) {
            GMRESMLMG gmres(mlmg);
            gmres.usePrecond(true);
            gmres.setMaxIters(100);
            gmres.solve(phi, rhs, rel, Real(0.0));
        } else {
            mlmg.solve({&phi}, {&rhs}, rel, Real(0.0));
        }
        MultiFab Ax(c.ba, c.dm, 1, 0);
        c.erf_apply(Ax, phi);
        MultiFab::Subtract(Ax, rhs, 0, 0, 1, 0);
        EXPECT_LE(Ax.norm2() / rhs.norm2(), Real(100.0) * rel) << (use_gmres ? "GMRESMLMG" : "MLMG");
    }
}

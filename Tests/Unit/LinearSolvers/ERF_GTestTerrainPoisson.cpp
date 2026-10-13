#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include <ERF_TerrainPoisson_3D_K.H>

#ifdef ERF_USE_FFT
#include <AMReX_GMRES.H>
#include <AMReX_MultiFab.H>
#include <ERF_TerrainPoisson.H>
#include <ERF_Utils.H>
#endif

#include <gtest/gtest.h>

#include <cmath>

using namespace amrex;

// The terrain Poisson operator that FFT-preconditioned GMRES applies on a terrain-fitted
// mesh (erf-model/ERF#4272).  A Neumann face of the solve region -- a coarse/fine face of a
// refined level, or a wall or inflow face of the domain -- used to get an even mirror of
// the solution in its ghost cell only.  That zeroes the normal-gradient part of the face
// flux but not the terrain cross term (h_xi/h_zeta times the vertical difference), so on
// sloping terrain the operator let flux through the face: the sum of dJ * (A x) over the
// region is then not zero, a right-hand side with zero dJ-weighted mean (what ERF hands the
// solver) is outside the operator's range, and GMRES cannot converge.
// The operator now sets the whole flux through those faces to zero; the tests below keep the
// mask-free form of terrpoisson_adotx (zf = 0) as the earlier operator to compare against.

namespace {

Real tol (Real scale) { return scale * (sizeof(Real) == 8 ? Real(1.e-12) : Real(1.e-4)); }

constexpr int NX = 8, NY = 6, NZ = 6;
constexpr Real DX = Real(80.0), DY = Real(60.0), DZ = Real(40.0);
constexpr Real ZTOP = NZ * DZ;

// A hill off the centre on a tilted plane: the terrain slopes, in x and in y, at every
// lateral face of the region
AMREX_GPU_HOST_DEVICE inline Real terrain (Real x, Real y, Real hmax)
{
    const Real xr = (x - Real(250.0)) / Real(300.0);
    const Real yr = (y - Real(200.0)) / Real(250.0);
    return hmax * (Real(1.0) / (Real(1.0) + xr*xr + yr*yr) + Real(0.05) * x / DX + Real(0.02) * y / DY);
}

// A smooth field that varies in every direction
inline Real field (int i, int j, int k)
{
    return std::sin(Real(0.7)*Real(i) + Real(0.3)) * std::cos(Real(0.5)*Real(j) - Real(0.2))
         + Real(0.3) * std::cos(Real(0.9)*Real(k) + Real(0.1)*Real(i*j)) + Real(0.05) * Real(i + 2*j - k);
}

// One region of NX x NY x NZ cells with the metrics ERF builds (make_J, make_areas) on
// host-accessible memory, and the even mirror that TerrainPoisson::apply_bcs puts in the
// ghost cells of every non-periodic face (x first, then y, then z, as apply_bcs does).
struct Region
{
    int K0;  // the region's lowest cell; above 0 its bottom is a coarse/fine face
    Box cells;
    FArrayBox zp, ax, ay, az, dJ, x, y;
    Dim3 dlo, dhi;

    explicit Region (Real hmax, int k0 = 0)
        : K0(k0),
          cells(IntVect(0, 0, k0), IntVect(NX-1, NY-1, k0+NZ-1)),
          zp(grow(convert(cells, IntVect(1)), 2), 1, The_Pinned_Arena()),
          ax(grow(surroundingNodes(cells, 0), 1), 1, The_Pinned_Arena()),
          ay(grow(surroundingNodes(cells, 1), 1), 1, The_Pinned_Arena()),
          az(grow(surroundingNodes(cells, 2), 1), 1, The_Pinned_Arena()),
          dJ(grow(cells, 1), 1, The_Pinned_Arena()),
          x(grow(cells, 1), 1, The_Pinned_Arena()),
          y(cells, 1, The_Pinned_Arena()),
          dlo{0, 0, k0}, dhi{NX-1, NY-1, k0+NZ-1}
    {
        auto z = zp.array();
        const Box& nb = zp.box();
        for (int k = nb.smallEnd(2); k <= nb.bigEnd(2); ++k) {
        for (int j = nb.smallEnd(1); j <= nb.bigEnd(1); ++j) {
        for (int i = nb.smallEnd(0); i <= nb.bigEnd(0); ++i) {
            const Real h = terrain(Real(i)*DX, Real(j)*DY, hmax);
            const Real zeta = Real(k) * DZ;
            z(i,j,k) = zeta + h * (Real(1.0) - zeta / ZTOP);   // basic terrain following
        }}}
        const Real dzinv = Real(1.0) / DZ;
        auto a = ax.array();
        for (int k = ax.box().smallEnd(2); k <= ax.box().bigEnd(2); ++k) {
        for (int j = ax.box().smallEnd(1); j <= ax.box().bigEnd(1); ++j) {
        for (int i = ax.box().smallEnd(0); i <= ax.box().bigEnd(0); ++i) {
            a(i,j,k) = Real(0.5) * dzinv * (z(i,j,k+1) + z(i,j+1,k+1) - z(i,j,k) - z(i,j+1,k));
        }}}
        auto b = ay.array();
        for (int k = ay.box().smallEnd(2); k <= ay.box().bigEnd(2); ++k) {
        for (int j = ay.box().smallEnd(1); j <= ay.box().bigEnd(1); ++j) {
        for (int i = ay.box().smallEnd(0); i <= ay.box().bigEnd(0); ++i) {
            b(i,j,k) = Real(0.5) * dzinv * (z(i,j,k+1) + z(i+1,j,k+1) - z(i,j,k) - z(i+1,j,k));
        }}}
        az.setVal<RunOn::Host>(Real(1.0));
        auto J = dJ.array();
        for (int k = dJ.box().smallEnd(2); k <= dJ.box().bigEnd(2); ++k) {
        for (int j = dJ.box().smallEnd(1); j <= dJ.box().bigEnd(1); ++j) {
        for (int i = dJ.box().smallEnd(0); i <= dJ.box().bigEnd(0); ++i) {
            J(i,j,k) = Real(0.25) * dzinv * (z(i,j,k+1) + z(i+1,j,k+1) + z(i,j+1,k+1) + z(i+1,j+1,k+1)
                                            -z(i,j,k  ) - z(i+1,j,k  ) - z(i,j+1,k  ) - z(i+1,j+1,k  ));
        }}}
    }

    void set_field (bool constant)
    {
        auto f = x.array();
        for (int k = K0; k < K0+NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            f(i,j,k) = constant ? Real(3.5) : field(i,j,k);
        }}}
        // even mirror, in the order of TerrainPoisson::apply_bcs
        for (int k = K0; k < K0+NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
            f(-1,j,k) = f(0,j,k);   f(NX,j,k) = f(NX-1,j,k);
        }}
        for (int k = K0; k < K0+NZ; ++k) {
        for (int i = -1; i <= NX; ++i) {
            f(i,-1,k) = f(i,0,k);   f(i,NY,k) = f(i,NY-1,k);
        }}
        for (int j = -1; j <= NY; ++j) {
        for (int i = -1; i <= NX; ++i) {
            f(i,j,K0-1) = f(i,j,K0);   f(i,j,K0+NZ) = f(i,j,K0+NZ-1);
        }}
    }

    // A x as TerrainPoisson::apply does it with face mask zf (0 = the old operator)
    void apply (int zf)
    {
        auto ya = y.array();
        auto xa = x.const_array();
        auto zpa = zp.const_array(); auto axa = ax.const_array(); auto aya = ay.const_array();
        auto aza = az.const_array(); auto dJa = dJ.const_array();
        for (int k = K0; k < K0+NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            terrpoisson_adotx(i, j, k, ya, xa, axa, aya, aza, dJa, zpa,
                              Real(1.0)/DX, Real(1.0)/DY, Real(1.0)/DZ, dlo, dhi, zf);
        }}}
    }

    // The same operator assembled independently, face by face from terrpoisson_flux_x/y/z
    // with the masked faces' fluxes set to zero, into yref
    void assemble (int zf, FArrayBox& yref) const
    {
        auto ya = yref.array();
        auto xa = x.const_array();
        auto zpa = zp.const_array(); auto axa = ax.const_array(); auto aya = ay.const_array();
        auto aza = az.const_array(); auto J = dJ.const_array();
        const Real dxi = Real(1.0)/DX, dyi = Real(1.0)/DY, dzi = Real(1.0)/DZ;
        for (int k = K0; k < K0+NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            Real fxl = terrpoisson_flux_x(i  ,j,k,xa,zpa,dxi), fxh = terrpoisson_flux_x(i+1,j,k,xa,zpa,dxi);
            Real fyl = terrpoisson_flux_y(i,j  ,k,xa,zpa,dyi), fyh = terrpoisson_flux_y(i,j+1,k,xa,zpa,dyi);
            Real fzl = terrpoisson_flux_z(i,j,k  ,xa,zpa,dxi,dyi), fzh = terrpoisson_flux_z(i,j,k+1,xa,zpa,dxi,dyi);
            if ((zf & tp_face_xlo) && i == 0)    { fxl = 0; }
            if ((zf & tp_face_xhi) && i == NX-1) { fxh = 0; }
            if ((zf & tp_face_ylo) && j == 0)    { fyl = 0; }
            if ((zf & tp_face_yhi) && j == NY-1) { fyh = 0; }
            if ((zf & tp_face_zlo) && k == K0 && K0 > 0) { fzl = 0; }
            if ((zf & tp_face_zhi) && k == K0+NZ-1) { fzh = 0; }
            // A x = -div(flux) / J
            ya(i,j,k) = -( (axa(i+1,j,k)*fxh - axa(i,j,k)*fxl) * dxi
                          +(aya(i,j+1,k)*fyh - aya(i,j,k)*fyl) * dyi
                          +(aza(i,j,k+1)*fzh - aza(i,j,k)*fzl) * dzi ) / J(i,j,k);
        }}}
    }

    // sum of dJ * (A x) * cell volume over the region, and the sum of its magnitudes
    std::pair<Real,Real> weighted_sum () const
    {
        auto ya = y.const_array(); auto J = dJ.const_array();
        Real s = 0, sa = 0;
        for (int k = K0; k < K0+NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            const Real v = J(i,j,k) * ya(i,j,k) * DX * DY * DZ;
            s += v; sa += std::abs(v);
        }}}
        return {s, sa};
    }
};

// Every face of the region that ERF solves with a Neumann condition: four lateral faces
// and the top (the ground flux is zero in both operators)
constexpr int ALL_NEUMANN = tp_face_xlo | tp_face_xhi | tp_face_ylo | tp_face_yhi | tp_face_zhi;

} // namespace

// With zero flux through every Neumann face the operator is a sum of flux differences
// that telescopes to nothing: sum dJ (A x) dV = 0 for any x (the discrete Gauss theorem,
// with dJ the operator's left null vector, #3487).  The old operator misses it by the
// terrain cross-term flux through the faces, which the second half recomputes face by face
// from terrpoisson_flux_x/y/z: the leak is that flux and nothing else.
TEST(TerrainPoisson, ZeroFluxOperatorConservesOnSlopes)
{
    Region r(Real(60.0));
    r.set_field(false);

    r.apply(ALL_NEUMANN);
    auto [s_zf, mag_zf] = r.weighted_sum();
    EXPECT_GT(mag_zf, Real(1.0));
    EXPECT_LE(std::abs(s_zf), tol(mag_zf)) << "sum dJ A x = " << s_zf << " of " << mag_zf;

    r.apply(0);
    auto [s_old, mag_old] = r.weighted_sum();
    // Measured: about 2 % of the magnitude sum on this hill; the old operator is not conservative
    EXPECT_GT(std::abs(s_old), Real(1.e-3) * mag_old) << "sum dJ A x = " << s_old << " of " << mag_old;

    // The leak is the flux the mirror leaves on the region's Neumann faces
    auto xa = r.x.const_array(); auto zpa = r.zp.const_array();
    auto axa = r.ax.const_array(); auto aya = r.ay.const_array(); auto aza = r.az.const_array();
    Real leak = 0;
    for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
        leak += -axa(NX,j,k) * terrpoisson_flux_x(NX,j,k,xa,zpa,Real(1.0)/DX) * DY * DZ
                +axa(0 ,j,k) * terrpoisson_flux_x(0 ,j,k,xa,zpa,Real(1.0)/DX) * DY * DZ;
    }}
    for (int k = 0; k < NZ; ++k) {
    for (int i = 0; i < NX; ++i) {
        leak += -aya(i,NY,k) * terrpoisson_flux_y(i,NY,k,xa,zpa,Real(1.0)/DY) * DX * DZ
                +aya(i,0 ,k) * terrpoisson_flux_y(i,0 ,k,xa,zpa,Real(1.0)/DY) * DX * DZ;
    }}
    for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {
        leak += -aza(i,j,NZ) * terrpoisson_flux_z(i,j,NZ,xa,zpa,Real(1.0)/DX,Real(1.0)/DY) * DX * DY;
    }}
    EXPECT_NEAR(s_old, leak, tol(Real(10.0) * mag_old));
}

// Next to one masked face, the zero-flux form differs from the old operator by that face's
// flux alone, and that flux is not zero on the slope (the comparison is not vacuous).
TEST(TerrainPoisson, OnlyTheMaskedFaceFluxChanges)
{
    Region r(Real(60.0));
    r.set_field(false);
    r.apply(0);
    FArrayBox y_old(r.cells, 1, The_Pinned_Arena());
    y_old.copy<RunOn::Host>(r.y);
    r.apply(tp_face_xlo);

    auto xa = r.x.const_array(); auto zpa = r.zp.const_array();
    auto axa = r.ax.const_array(); auto J = r.dJ.const_array();
    auto yo = y_old.const_array(); auto yn = r.y.const_array();
    for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
        const Real f = terrpoisson_flux_x(0,j,k,xa,zpa,Real(1.0)/DX);
        EXPECT_GT(std::abs(f), Real(1.e-6)) << "no cross-term flux at (0," << j << "," << k << ")";
        const Real expected = yo(0,j,k) - axa(0,j,k) * f / DX / J(0,j,k);
        EXPECT_NEAR(yn(0,j,k), expected, tol(Real(10.0) * std::abs(yo(0,j,k)) + std::abs(f)));
        for (int i = 1; i < NX; ++i) {
            EXPECT_EQ(yn(i,j,k), yo(i,j,k)) << "cell (" << i << "," << j << "," << k << ") changed";
        }
    }}
}

// terrpoisson_adotx sets the masked faces' fluxes to zero in its own (inlined) stencil.
// Assembled independently from terrpoisson_flux_x/y/z with those fluxes set to zero, the
// operator must be the same at every cell, with and without masked faces.
TEST(TerrainPoisson, ZeroFluxFormMatchesFaceByFaceAssembly)
{
    Region r(Real(60.0));
    r.set_field(false);
    FArrayBox yref(r.cells, 1, The_Pinned_Arena());
    for (int zf : {0, ALL_NEUMANN, int(tp_face_xlo | tp_face_yhi)}) {
        r.apply(zf);
        r.assemble(zf, yref);
        auto ya = r.y.const_array(); auto yr = yref.const_array();
        for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            EXPECT_NEAR(ya(i,j,k), yr(i,j,k), tol(Real(100.0) * (std::abs(yr(i,j,k)) + Real(1.e-3))))
                << "mask " << zf << ", cell (" << i << "," << j << "," << k << ")";
        }}}
    }
}

// On flat terrain the cross term is zero, the mirror already gives zero flux, and the two
// operators agree bit for bit; a constant is in the null space of the zero-flux operator.
TEST(TerrainPoisson, FlatTerrainAndConstants)
{
    Region flat(Real(0.0));
    flat.set_field(false);
    flat.apply(0);
    FArrayBox y_old(flat.cells, 1, The_Pinned_Arena());
    y_old.copy<RunOn::Host>(flat.y);
    flat.apply(ALL_NEUMANN);
    auto yo = y_old.const_array(); auto yn = flat.y.const_array();
    for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {
        EXPECT_EQ(yn(i,j,k), yo(i,j,k)) << "cell (" << i << "," << j << "," << k << ")";
    }}}

    Region hill(Real(60.0));
    hill.set_field(true);
    hill.apply(ALL_NEUMANN);
    auto yc = hill.y.const_array();
    for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {
        EXPECT_NEAR(yc(i,j,k), Real(0.0), tol(Real(1.0)));
    }}}
}

// A refined region that starts above the ground (its bottom is a coarse/fine face) and stops
// below the top: with zero flux through all six faces the operator conserves again; leaving
// the bottom face out of the mask leaks through it, and the face-by-face assembly agrees.
TEST(TerrainPoisson, RegionAboveTheGroundConserves)
{
    Region r(Real(60.0), 3);
    r.set_field(false);
    const int all = ALL_NEUMANN | tp_face_zlo;
    r.apply(all);
    auto [s_zf, mag_zf] = r.weighted_sum();
    EXPECT_GT(mag_zf, Real(1.0));
    EXPECT_LE(std::abs(s_zf), tol(mag_zf)) << "sum dJ A x = " << s_zf << " of " << mag_zf;

    r.apply(ALL_NEUMANN);
    auto [s_open, mag_open] = r.weighted_sum();
    EXPECT_GT(std::abs(s_open), Real(1.e-3) * mag_open) << "sum dJ A x = " << s_open;

    FArrayBox yref(r.cells, 1, The_Pinned_Arena());
    r.apply(all);
    r.assemble(all, yref);
    auto ya = r.y.const_array(); auto yr = yref.const_array();
    for (int k = r.K0; k < r.K0+NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {
        EXPECT_NEAR(ya(i,j,k), yr(i,j,k), tol(Real(100.0) * (std::abs(yr(i,j,k)) + Real(1.e-3))));
    }}}
}

#ifdef ERF_USE_FFT
namespace {

// Node heights of the hill (a free function: CUDA does not allow a device lambda in a constructor).
// With flat_top the hill's shape decays to nothing at the top, as ERF's terrain-following meshes
// do; otherwise half of it is left there, so the top face slopes too.
void fill_heights (MultiFab& zp, bool flat_top)
{
    const Real decay = flat_top ? Real(1.0) : Real(0.5);
    for (MFIter mfi(zp); mfi.isValid(); ++mfi) {
        auto z = zp.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real h = terrain(Real(i)*DX, Real(j)*DY, Real(60.0));
            const Real zeta = Real(k) * DZ;
            z(i,j,k) = zeta + h * (Real(1.0) - decay * zeta / ZTOP);
        });
    }
}

// The same hill on a MultiFab split into boxes, with the metrics from make_J/make_areas,
// and a TerrainPoisson on a region with slip walls on every side
struct ClassCase
{
    Geometry geom;
    BoxArray ba;
    DistributionMapping dm;
    MultiFab zp, ax, ay, az, dJ;
    Array<std::string,2*AMREX_SPACEDIM> bcs;
    Gpu::DeviceVector<Real> dz_d;

    explicit ClassCase (bool flat_top = true)
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
        fill_heights(zp, flat_top);
        make_J(geom, zp, dJ);
        make_areas(geom, zp, ax, ay, az);
        for (auto& s : bcs) { s = "SlipWall"; }
        dz_d.resize(NZ);
        Gpu::HostVector<Real> dz_h(NZ, DZ);
        Gpu::copy(Gpu::hostToDevice, dz_h.begin(), dz_h.end(), dz_d.begin());
        Gpu::streamSynchronize();
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
};

// The earlier operator on the same region, for comparison: terrpoisson_adotx without masked
// faces (zf = 0), with the ghost cells of x filled as TerrainPoisson fills them
void apply_unmasked (ClassCase const& c, TerrainPoisson& tp, MultiFab& x, MultiFab& y)
{
    tp.apply_bcs(x);
    auto const& dxinv = c.geom.InvCellSizeArray();
    auto const& xa = x.const_arrays();
    auto const& ya = y.arrays();
    auto const& zpa = c.zp.const_arrays();
    auto const& axa = c.ax.const_arrays();
    auto const& aya = c.ay.const_arrays();
    auto const& aza = c.az.const_arrays();
    auto const& dJa = c.dJ.const_arrays();
    ParallelFor(y, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        terrpoisson_adotx(i, j, k, ya[b], xa[b], axa[b], aya[b], aza[b], dJa[b], zpa[b],
                          dxinv[0], dxinv[1], dxinv[2]);
    });
    Gpu::streamSynchronize();
}

// The face fluxes of the earlier operator: terrpoisson_flux_x/y/z on every face
void fluxes_unmasked (ClassCase const& c, TerrainPoisson& tp, MultiFab& x, Array<MultiFab,AMREX_SPACEDIM>& f)
{
    tp.apply_bcs(x);
    auto const& dxinv = c.geom.InvCellSizeArray();
    auto const& xa = x.const_arrays();
    auto const& zpa = c.zp.const_arrays();
    auto const& fx = f[0].arrays();
    auto const& fy = f[1].arrays();
    auto const& fz = f[2].arrays();
    ParallelFor(f[0], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    { fx[b](i,j,k) = terrpoisson_flux_x(i,j,k,xa[b],zpa[b],dxinv[0]); });
    ParallelFor(f[1], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    { fy[b](i,j,k) = terrpoisson_flux_y(i,j,k,xa[b],zpa[b],dxinv[1]); });
    ParallelFor(f[2], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    { fz[b](i,j,k) = terrpoisson_flux_z(i,j,k,xa[b],zpa[b],dxinv[0],dxinv[1]); });
    Gpu::streamSynchronize();
}

} // namespace

// The class: its mask names the Neumann faces of the region (not the ground), apply conserves
// (dot product with dJ is zero) where the earlier operator, terrpoisson_adotx without the mask,
// does not, and getFluxes gives zero on those faces and the unmasked fluxes everywhere else.
TEST(TerrainPoisson, ClassMaskApplyAndFluxes)
{
    // A flat top (as ERF's meshes have) and a sloping one, so the top face's mask is exercised too
    for (bool flat_top : {true, false}) {
        ClassCase c(flat_top);
        TerrainPoisson tp(c.geom, c.geom, c.ba, c.dm, c.bcs, c.dz_d, c.ax, c.ay, c.az, c.dJ, &c.zp, false);
        EXPECT_EQ(tp.zeroFluxFaces(), ALL_NEUMANN);

        MultiFab x(c.ba, c.dm, 1, 1), y(c.ba, c.dm, 1, 0);
        c.set_field(x);
        tp.apply(y, x);
        const Real s_zf = MultiFab::Dot(c.dJ, 0, y, 0, 1, 0);
        const Real mag = y.norm1(0) * c.dJ.max(0);
        EXPECT_LE(std::abs(s_zf), tol(mag)) << "dJ . A x = " << s_zf;
        apply_unmasked(c, tp, x, y);
        const Real s_old = MultiFab::Dot(c.dJ, 0, y, 0, 1, 0);
        EXPECT_GT(std::abs(s_old), Real(1.e-3) * mag) << "dJ . A x = " << s_old;

        Array<MultiFab,AMREX_SPACEDIM> f_zf, f_old;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            f_zf[d].define(convert(c.ba, IntVect::TheDimensionVector(d)), c.dm, 1, 0);
            f_old[d].define(convert(c.ba, IntVect::TheDimensionVector(d)), c.dm, 1, 0);
        }
        tp.getFluxes(x, f_zf);
        fluxes_unmasked(c, tp, x, f_old);
        const int dhi[3] = {NX, NY, NZ};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            Real on_face_zf = 0, on_face_old = 0, inside_diff = 0, inside_max = 0;
            for (MFIter mfi(f_zf[d]); mfi.isValid(); ++mfi) {
                const Box& b = mfi.validbox();
                FArrayBox ah(b, 1, The_Pinned_Arena()), oh(b, 1, The_Pinned_Arena());
                ah.copy<RunOn::Device>(f_zf[d][mfi]);
                oh.copy<RunOn::Device>(f_old[d][mfi]);
                Gpu::streamSynchronize();
                auto a = ah.const_array();
                auto o = oh.const_array();
                for (int k = b.smallEnd(2); k <= b.bigEnd(2); ++k) {
                for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j) {
                for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i) {
                    const int n = (d == 0) ? i : ((d == 1) ? j : k);
                    const bool masked = (d < 2 && (n == 0 || n == dhi[d])) || (d == 2 && n == dhi[d]);
                    if (masked) {
                        on_face_zf  = amrex::max(on_face_zf , std::abs(a(i,j,k)));
                        on_face_old = amrex::max(on_face_old, std::abs(o(i,j,k)));
                    } else {
                        inside_diff = amrex::max(inside_diff, std::abs(a(i,j,k) - o(i,j,k)));
                        inside_max  = amrex::max(inside_max , std::abs(o(i,j,k)));
                    }
                }}}
            }
            EXPECT_EQ(on_face_zf, Real(0.0)) << "direction " << d;
            // the same flux functions compiled in two files: equal up to contraction into FMAs
            EXPECT_LE(inside_diff, tol(inside_max)) << "direction " << d;
            if (d < 2 || !flat_top) {
                EXPECT_GT(on_face_old, Real(1.e-6)) << "direction " << d << ": no slope at the faces";
            }
        }
    }
}

// The singular solve ERF makes on a refined level or a walled domain: a right-hand side with
// zero dJ-weighted mean.  It is in the range of the zero-flux operator, and FFT-preconditioned
// GMRES converges to 1e-10 (14 iterations).  The earlier operator's range is not orthogonal to
// dJ (ClassMaskApplyAndFluxes), so that system had no solution: GMRES reported 5e-11 after 35
// iterations while the solution grew to 8e17 and the true residual was 0.92.  Hence the true
// residual and the size of the solution are checked here, not GMRES's own estimate.
TEST(TerrainPoisson, SingularSolveConvergesWithZeroFlux)
{
    ClassCase c;
    MultiFab rhs(c.ba, c.dm, 1, 0);
    c.set_field(rhs);
    const Real mean = MultiFab::Dot(c.dJ, 0, rhs, 0, 1, 0) / c.dJ.sum(0);
    rhs.plus(-mean, 0, 1, 0);

    TerrainPoisson tp(c.geom, c.geom, c.ba, c.dm, c.bcs, c.dz_d, c.ax, c.ay, c.az, c.dJ, &c.zp, false);
    GMRES<MultiFab, TerrainPoisson> gmres;
    gmres.define(tp);
    gmres.setRestartLength(50);
    gmres.setMaxIters(200);
    tp.usePrecond(true);
    MultiFab phi(c.ba, c.dm, 1, 1);
    phi.setVal(0.0);
    const Real rel = (sizeof(Real) == 8) ? Real(1.e-10) : Real(1.e-5);
    gmres.solve(phi, rhs, rel, Real(0.0));
    MultiFab Ax(c.ba, c.dm, 1, 0);
    tp.apply(Ax, phi);
    MultiFab::Subtract(Ax, rhs, 0, 0, 1, 0);
    const Real res = Ax.norm2() / rhs.norm2();
    EXPECT_LE(res, Real(1.01) * rel) << gmres.getNumIters() << " iterations";
    EXPECT_LE(phi.norm0(), Real(1.e6)) << gmres.getNumIters() << " iterations";
}
#endif

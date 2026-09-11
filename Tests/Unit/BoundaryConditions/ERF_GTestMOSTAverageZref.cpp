#include "AMReX.H"
#include "AMReX_FArrayBox.H"

#include "ERF_MOSTAverage.H"

#include <gtest/gtest.h>

#include <cmath>

using namespace amrex;

namespace {

// Values reach ~100, so single precision needs a looser absolute tolerance
const Real tol = (sizeof(Real) == 8) ? Real(1.0e-12) : Real(1.0e-4);

// Staggered heights of a column stretched from dz0 by ratio, as init_zlevels
// builds them from erf.initial_dz and erf.grid_stretching_ratio
Vector<Real>
stretched_levels (int nz, Real dz0, Real ratio)
{
    Vector<Real> z(nz+1);
    z[0] = Real(0.0);
    Real dz = dz0;
    for (int k = 1; k <= nz; ++k) {
        z[k] = z[k-1] + dz;
        dz *= ratio;
    }
    return z;
}

} // namespace

// A height on a face belongs to the cell above it, and every height inside the
// column belongs to exactly one cell
TEST(MOSTAverageZref, InCellIsHalfOpen)
{
    EXPECT_TRUE (MOSTAverage::in_cell_z(Real(10.0), Real(10.0), Real(21.0)));
    EXPECT_FALSE(MOSTAverage::in_cell_z(Real(10.0), Real( 0.0), Real(10.0)));
    EXPECT_TRUE (MOSTAverage::in_cell_z(Real( 5.0), Real( 0.0), Real(10.0)));

    const auto z = stretched_levels(40, Real(10.0), Real(1.1));
    for (Real zq : {Real(0.0), Real(5.0), Real(10.0), Real(21.0), Real(33.1), Real(100.0)}) {
        int ncell = 0;
        for (int k = 0; k < 40; ++k) {
            if (MOSTAverage::in_cell_z(zq, z[k], z[k+1])) { ++ncell; }
        }
        EXPECT_EQ(ncell, 1) << "zq = " << zq;
    }
}

// On a stretched column the default reference height is the true first cell
// center (5 m for dz0 = 10 m), not half of the uniform (prob_hi - prob_lo)/nz
TEST(MOSTAverageZref, StretchedCellCenters)
{
    const auto z = stretched_levels(40, Real(10.0), Real(1.1));

    const Real c0 = MOSTAverage::cell_center_height(z, 0);
    const Real c1 = MOSTAverage::cell_center_height(z, 1);
    EXPECT_NEAR(c0, Real( 5.0), tol);
    EXPECT_NEAR(c1, Real(15.5), tol);

    // A height on a cell center selects that cell
    EXPECT_EQ(MOSTAverage::k_index_below(z, c0 - Real(0.01)), -1);
    EXPECT_EQ(MOSTAverage::k_index_below(z, c0),               0);
    EXPECT_EQ(MOSTAverage::k_index_below(z, Real(10.0)),       0);
    EXPECT_EQ(MOSTAverage::k_index_below(z, c1 - Real(0.01)),  0);
    EXPECT_EQ(MOSTAverage::k_index_below(z, c1),               1);
    EXPECT_EQ(MOSTAverage::k_index_below(z, Real(55.4)),       4);
    EXPECT_EQ(MOSTAverage::k_index_below(z, Real(1.0e6)),     39);
}

// On a uniform column the index matches the uniform-mesh formula
// floor(zref/dz - 1/2) used for ConstantDz
TEST(MOSTAverageZref, UniformMatchesFloor)
{
    const Real dz = Real(10.0);
    const auto z = stretched_levels(20, dz, Real(1.0));
    for (Real zq : {Real(5.0), Real(7.5), Real(10.0), Real(14.99), Real(15.0), Real(33.3), Real(195.0)}) {
        const int lk_floor = static_cast<int>(std::floor(zq / dz - Real(0.5)));
        EXPECT_EQ(MOSTAverage::k_index_below(z, zq), lk_floor) << "zq = " << zq;
    }
}

namespace {

// Flat column of nodal heights from zlev and a cell-centered field equal to the
// cell-center height (mirrored about the bottom and top faces in the ghost
// cells), so an interpolation that is linear in z returns the query height.
struct FlatColumn
{
    FArrayBox z_fab;
    FArrayBox f_fab;

    explicit FlatColumn (const Vector<Real>& zlev)
    {
        const int nz = static_cast<int>(zlev.size()) - 1;
        z_fab.resize(Box(IntVect(-1,-1,-1), IntVect(4,4,nz+1)), 1, The_Pinned_Arena());
        f_fab.resize(Box(IntVect(-1,-1,-1), IntVect(3,3,nz  )), 1, The_Pinned_Arena());
        auto z_arr = z_fab.array();
        auto f_arr = f_fab.array();
        const Box zbx = z_fab.box();
        const Box fbx = f_fab.box();
        for (int k = zbx.smallEnd(2); k <= zbx.bigEnd(2); ++k) {
            const Real zk = (k < 0)  ? Real(2.0) * zlev[0]  - zlev[1]
                          : (k > nz) ? Real(2.0) * zlev[nz] - zlev[nz-1] : zlev[k];
            for (int j = zbx.smallEnd(1); j <= zbx.bigEnd(1); ++j) {
                for (int i = zbx.smallEnd(0); i <= zbx.bigEnd(0); ++i) { z_arr(i,j,k) = zk; }
            }
        }
        for (int k = fbx.smallEnd(2); k <= fbx.bigEnd(2); ++k) {
            const Real c0 = Real(0.5) * (zlev[0] + zlev[1]);
            const Real cn = Real(0.5) * (zlev[nz-1] + zlev[nz]);
            const Real fk = (k < 0)   ? Real(2.0) * zlev[0]  - c0
                          : (k >= nz) ? Real(2.0) * zlev[nz] - cn : Real(0.5) * (zlev[k] + zlev[k+1]);
            for (int j = fbx.smallEnd(1); j <= fbx.bigEnd(1); ++j) {
                for (int i = fbx.smallEnd(0); i <= fbx.bigEnd(0); ++i) { f_arr(i,j,k) = fk; }
            }
        }
    }

    // Interpolate at height zp above the center of horizontal cell (1,1)
    Real interp (Real zp) const
    {
        const GpuArray<Real,AMREX_SPACEDIM> plo{Real(0.0), Real(0.0), Real(0.0)};
        const GpuArray<Real,AMREX_SPACEDIM> dxi{Real(0.01), Real(0.01), Real(1.0)};
        Real val = Real(0.0);
        MOSTAverage::trilinear_interp_T(plo[0] + Real(1.5) / dxi[0], plo[1] + Real(1.5) / dxi[1], zp,
                                        &val, f_fab.const_array(), z_fab.const_array(), plo, dxi, 1);
        return val;
    }
};

} // namespace

// A query exactly on a mesh face of a flat stretched column is found, and the
// interpolation is linear in the physical height: at the 10 m face between the
// 5 m and 15.5 m cell centers it returns the value at 10 m, not the plain mean
// of the two cells (which belongs to 10.25 m)
TEST(MOSTAverageZref, InterpolationOnFace)
{
    const auto zlev = stretched_levels(8, Real(10.0), Real(1.1));
    const FlatColumn col(zlev);
    const Real ztol = (sizeof(Real) == 8) ? Real(1.0e-12) : Real(1.0e-4);

    // On the faces between cells 0 and 1 (10 m) and cells 1 and 2 (21 m)
    EXPECT_NEAR(col.interp(zlev[1]), zlev[1], ztol);
    EXPECT_NEAR(col.interp(zlev[2]), zlev[2], ztol);

    // At cell centers, inside cells, and below the first cell center
    for (Real zq : {Real(2.0), Real(5.0), Real(7.5), Real(12.3), Real(15.5), Real(25.0), Real(40.0)}) {
        EXPECT_NEAR(col.interp(zq), zq, ztol) << "zq = " << zq;
    }

    // Continuous across the face
    const Real eps = Real(1.0e-3);
    EXPECT_NEAR(col.interp(zlev[1] - eps), col.interp(zlev[1]), Real(2.0e-3));
    EXPECT_NEAR(col.interp(zlev[1] + eps), col.interp(zlev[1]), Real(2.0e-3));
}

// On equal cell heights the physical-height weight is the index weight
// (lk + fraction of the containing cell + 1/2) the interpolation used before
TEST(MOSTAverageZref, InterpolationUniformUnchanged)
{
    const Real dz = Real(10.0);
    const auto zlev = stretched_levels(8, dz, Real(1.0));
    const FlatColumn col(zlev);
    const Real ztol = (sizeof(Real) == 8) ? Real(1.0e-12) : Real(1.0e-4);
    for (Real zq : {Real(1.0), Real(5.0), Real(10.0), Real(13.7), Real(20.0), Real(44.4)}) {
        EXPECT_NEAR(col.interp(zq), zq, ztol) << "zq = " << zq;
    }
}

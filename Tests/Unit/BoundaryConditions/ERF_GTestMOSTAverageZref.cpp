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

// A query exactly on a mesh face of a flat terrain-fitted column interpolates
// between the two cells it separates instead of failing the height search
TEST(MOSTAverageZref, InterpolationOnFace)
{
    const int nz = 8;
    const auto zlev = stretched_levels(nz, Real(10.0), Real(1.1));

    const GpuArray<Real,AMREX_SPACEDIM> plo{Real(0.0), Real(0.0), Real(0.0)};
    const GpuArray<Real,AMREX_SPACEDIM> dxi{Real(0.01), Real(0.01), Real(1.0)};

    // Nodal heights (flat) and a cell-centered field that varies only in k
    FArrayBox z_fab(Box(IntVect(-1,-1,-1), IntVect(4,4,nz+1)), 1, The_Pinned_Arena());
    FArrayBox f_fab(Box(IntVect(-1,-1,-1), IntVect(3,3,nz  )), 1, The_Pinned_Arena());
    auto z_arr = z_fab.array();
    auto f_arr = f_fab.array();
    const Box zbx = z_fab.box();
    const Box fbx = f_fab.box();
    for (int k = zbx.smallEnd(2); k <= zbx.bigEnd(2); ++k) {
        const int kc = std::min(std::max(k, 0), nz);
        for (int j = zbx.smallEnd(1); j <= zbx.bigEnd(1); ++j) {
            for (int i = zbx.smallEnd(0); i <= zbx.bigEnd(0); ++i) {
                z_arr(i,j,k) = (k < 0) ? -zlev[1] : ((k > nz) ? zlev[nz] + zlev[1] : zlev[kc]);
            }
        }
    }
    for (int k = fbx.smallEnd(2); k <= fbx.bigEnd(2); ++k) {
        for (int j = fbx.smallEnd(1); j <= fbx.bigEnd(1); ++j) {
            for (int i = fbx.smallEnd(0); i <= fbx.bigEnd(0); ++i) {
                f_arr(i,j,k) = Real(100.0) + Real(k);
            }
        }
    }

    // Center of cell (1,1) in the horizontal
    const Real xp = plo[0] + Real(1.5) / dxi[0];
    const Real yp = plo[1] + Real(1.5) / dxi[1];

    auto interp = [&] (Real zp) {
        Real val = Real(0.0);
        MOSTAverage::trilinear_interp_T(xp, yp, zp, &val, f_fab.const_array(),
                                        z_fab.const_array(), plo, dxi, 1);
        return val;
    };

    // On the face between cells 0 and 1 (10 m) and between cells 1 and 2 (21 m)
    EXPECT_NEAR(interp(zlev[1]), Real(100.5), tol);
    EXPECT_NEAR(interp(zlev[2]), Real(101.5), tol);

    // At the cell centers the field itself
    EXPECT_NEAR(interp(Real(5.0)),  Real(100.0), tol);
    EXPECT_NEAR(interp(Real(15.5)), Real(101.0), tol);

    // Continuous across the face
    const Real eps = Real(1.0e-3);
    EXPECT_NEAR(interp(zlev[1] - eps), interp(zlev[1]), Real(1.0e-3));
    EXPECT_NEAR(interp(zlev[1] + eps), interp(zlev[1]), Real(1.0e-3));
}

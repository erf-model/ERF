// Unit tests for erf.rans_shear_stencil = log_law (ERF_LogShear.H, erf-model/ERF#5261): the
// vertical shear in S_mn S_mn of the k-eqn RANS closure, exact for a logarithmic profile.
//
// A column of cells carries u = (u_* / kappa) ln((d + z0) / z0) and v = c_v u, the edge strains
// are formed as ComputeStrain_{N,S,T} form them (difference of the two u points over their height
// difference), and the corrected S_mn S_mn is compared with the closed form of the profile.

#include <ERF_LogShear.H>
#include <ERF_EddyViscosity.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <limits>
#include <vector>

using namespace amrex;

namespace {

constexpr int  nz     = 10;
constexpr Real kappa  = Real(0.41);
constexpr Real ustar  = Real(0.4);
constexpr Real z0     = Real(0.1);
constexpr Real c_v    = Real(-0.35);  // v = c_v u: the wind is not along x
constexpr Real dwdx   = Real(3.0e-3); // a horizontal gradient of w inside s13, to be kept
constexpr Real s11val = Real(2.0e-3); // a normal strain, to be kept

Real log_u (Real d) { return ustar / kappa * std::log((d + z0) / z0); }

// The column: node heights, cell-center wall distances and the two S_mn S_mn of every cell
struct Column {
    std::vector<Real> zf, d;
    std::vector<Real> central, corrected;
};

enum class Profile { Log, Linear };

/**
 * S_mn S_mn of the column from ComputeSmnSmn (central) and from log_shear::corrected_SmnSmn.
 *
 * @param[in] zf node heights of the nz cells (nz + 1 values, zf[0] = 0 at the ground)
 * @param[in] profile the velocity profile
 * @param[in] d_scale factor on the wall distance handed to the kernel in the cells from kbreak up
 * @param[in] kbreak first cell whose wall distance is scaled (nz: none)
 */
Column
run_column (const std::vector<Real>& zf, Profile profile, Real d_scale = Real(1), int kbreak = nz)
{
    Column col;
    col.zf = zf;
    col.d.resize(nz);
    for (int k = 0; k < nz; ++k) { col.d[k] = Real(0.5) * (zf[k] + zf[k+1]); }

    const Box cc(IntVect(0,0,0), IntVect(1,1,nz-1));
    const Box gbx = grow(cc, 1);
    FArrayBox fu(convert(gbx, IntVect(1,0,0)), 1, The_Pinned_Arena());
    FArrayBox fv(convert(gbx, IntVect(0,1,0)), 1, The_Pinned_Arena());
    FArrayBox fz(convert(gbx, IntVect(1,1,1)), 1, The_Pinned_Arena());
    FArrayBox fd(gbx, 1, The_Pinned_Arena());
    FArrayBox f11(gbx, 1, The_Pinned_Arena()), f22(gbx, 1, The_Pinned_Arena()), f33(gbx, 1, The_Pinned_Arena());
    FArrayBox f12(convert(gbx, IntVect(1,1,0)), 1, The_Pinned_Arena());
    FArrayBox f13(convert(gbx, IntVect(1,0,1)), 1, The_Pinned_Arena());
    FArrayBox f23(convert(gbx, IntVect(0,1,1)), 1, The_Pinned_Arena());
    for (FArrayBox* f : {&fu, &fv, &fz, &fd, &f11, &f22, &f33, &f12, &f13, &f23}) {
        f->setVal<RunOn::Host>(Real(0));
    }
    f11.setVal<RunOn::Host>(s11val);

    // Heights of the nodes and of the cell centers, with one layer mirrored below and above
    std::vector<Real> zn(nz + 3), zc(nz + 2), uc(nz + 2);
    for (int k = 0; k <= nz; ++k) { zn[k+1] = zf[k]; }
    zn[0]    = -zf[1];
    zn[nz+2] = zf[nz] + (zf[nz] - zf[nz-1]);
    for (int k = -1; k <= nz; ++k) {
        zc[k+1] = Real(0.5) * (zn[k+1] + zn[k+2]);
        const Real dk = std::abs(zc[k+1]);
        uc[k+1] = (profile == Profile::Log) ? log_u(dk) : Real(0.02) * dk;
    }

    auto u = fu.array(); auto v = fv.array(); auto z = fz.array(); auto d = fd.array();
    auto s13 = f13.array(); auto s23 = f23.array();
    for (int k = -1; k <= nz; ++k) {
        for (int j = -1; j <= 3; ++j) {
            for (int i = -1; i <= 3; ++i) {
                if (fu.box().contains(IntVect(i,j,k))) { u(i,j,k) = uc[k+1]; }
                if (fv.box().contains(IntVect(i,j,k))) { v(i,j,k) = c_v * uc[k+1]; }
                if (fd.box().contains(IntVect(i,j,k))) {
                    d(i,j,k) = std::abs(zc[k+1]) * ((k >= kbreak) ? d_scale : Real(1));
                }
            }
        }
    }
    for (int k = -1; k <= nz + 1; ++k) {
        for (int j = -1; j <= 3; ++j) {
            for (int i = -1; i <= 3; ++i) {
                if (fz.box().contains(IntVect(i,j,k))) { z(i,j,k) = zn[k+1]; }
                // Edge k lies between the cells k-1 and k
                if (k >= 0 && k <= nz) {
                    const Real idz = Real(1) / (zc[k+1] - zc[k]);
                    if (f13.box().contains(IntVect(i,j,k))) {
                        s13(i,j,k) = Real(0.5) * ((uc[k+1] - uc[k]) * idz + dwdx);
                    }
                    if (f23.box().contains(IntVect(i,j,k))) {
                        s23(i,j,k) = Real(0.5) * c_v * (uc[k+1] - uc[k]) * idz;
                    }
                }
            }
        }
    }

    FArrayBox fold(cc, 1, The_Pinned_Arena()), fnew(cc, 1, The_Pinned_Arena());
    fold.setVal<RunOn::Host>(Real(0));
    fnew.setVal<RunOn::Host>(Real(0));
    const Array4<Real const> a11 = f11.const_array(), a22 = f22.const_array(), a33 = f33.const_array();
    const Array4<Real const> a12 = f12.const_array(), a13 = f13.const_array(), a23 = f23.const_array();
    const Array4<Real const> au = fu.const_array(), av = fv.const_array();
    const Array4<Real const> ad = fd.const_array(), az = fz.const_array();
    const Array4<Real> o_old = fold.array(), o_new = fnew.array();
    const Real z0_l = z0;
    Box inner = cc;
    inner.growLo(2, -1).growHi(2, -1); // above the wall cell and below the top cell
    ParallelFor(inner, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        o_old(i,j,k) = ComputeSmnSmn(i, j, k, a11, a22, a33, a12, a13, a23);
        o_new(i,j,k) = log_shear::corrected_SmnSmn(i, j, k, o_old(i,j,k), au, av, a13, a23, ad, z0_l, az);
    });
    Gpu::streamSynchronize();

    col.central.resize(nz);
    col.corrected.resize(nz);
    for (int k = 0; k < nz; ++k) {
        col.central[k]   = fold.const_array()(1,1,k);
        col.corrected[k] = fnew.const_array()(1,1,k);
    }
    return col;
}

// S_mn S_mn of the test fields for a vertical gradient g of u (and c_v g of v)
Real SmnSmn_of (Real g)
{
    const Real e13 = Real(0.5) * (g + dwdx);
    const Real e23 = Real(0.5) * c_v * g;
    return s11val * s11val + Real(2) * e13 * e13 + Real(2) * e23 * e23;
}

std::vector<Real> uniform_layers (Real dz)
{
    std::vector<Real> zf(nz + 1);
    for (int k = 0; k <= nz; ++k) { zf[k] = dz * Real(k); }
    return zf;
}

std::vector<Real> stretched_layers (Real dz0, Real ratio)
{
    std::vector<Real> zf(nz + 1);
    zf[0] = Real(0);
    Real dz = dz0;
    for (int k = 1; k <= nz; ++k) { zf[k] = zf[k-1] + dz; dz *= ratio; }
    return zf;
}

Real rel_tol ()
{
    return (sizeof(Real) == sizeof(double)) ? Real(1e-11) : Real(2e-4);
}

} // namespace

// Uniform layers of Neutral_ABL_Flat (dz = 15.625 m, z0 = 0.1 m): the central difference of the log
// profile is larger than the gradient at the cell center by 1.2046, 1.0588 and 1.0285 in the cells
// 2, 3 and 4 (the numbers of the issue), and the corrected S_mn S_mn is the closed form.
TEST(LogShear, UniformLayersLogProfileIsExact)
{
    const Real dz = Real(15.625);
    const Column col = run_column(uniform_layers(dz), Profile::Log);
    const Real pinned[3] = {Real(1.2046), Real(1.0588), Real(1.0285)};
    for (int k = 1; k < nz - 1; ++k) {
        const Real g_exact   = ustar / (kappa * (col.d[k] + z0));
        const Real g_central = (log_u(col.d[k+1]) - log_u(col.d[k-1])) / (Real(2) * dz);
        if (k <= 3) { EXPECT_NEAR(g_central / g_exact, pinned[k-1], Real(6e-5)) << "cell " << k + 1; }
        EXPECT_NEAR(col.central[k],   SmnSmn_of(g_central), rel_tol() * SmnSmn_of(g_central)) << "k=" << k;
        EXPECT_NEAR(col.corrected[k], SmnSmn_of(g_exact),   rel_tol() * SmnSmn_of(g_exact))   << "k=" << k;
    }
    // The corrected value differs from the central one by far more than the tolerance
    EXPECT_GT(col.central[1], Real(1.3) * col.corrected[1]);
}

// Stretched layers (ratio 1.2): the strain of ComputeStrain_S is the mean of the two edge
// differences, and the corrected S_mn S_mn is again the closed form.
TEST(LogShear, StretchedLayersLogProfileIsExact)
{
    const Column col = run_column(stretched_layers(Real(4), Real(1.2)), Profile::Log);
    for (int k = 1; k < nz - 1; ++k) {
        const Real g_exact = ustar / (kappa * (col.d[k] + z0));
        const Real g_edges = Real(0.5) * ( (log_u(col.d[k+1]) - log_u(col.d[k]))   / (col.d[k+1] - col.d[k])
                                         + (log_u(col.d[k])   - log_u(col.d[k-1])) / (col.d[k]   - col.d[k-1]) );
        EXPECT_NEAR(col.central[k],   SmnSmn_of(g_edges), rel_tol() * SmnSmn_of(g_edges)) << "k=" << k;
        EXPECT_NEAR(col.corrected[k], SmnSmn_of(g_exact), rel_tol() * SmnSmn_of(g_exact)) << "k=" << k;
        EXPECT_GT(g_edges, g_exact);
    }
}

// A linear profile u = a d: the stencil returns a (d_{k+1} - d_{k-1}) / [(d_k + z0) ln(...)], which
// is below a (the stencil is exact for the log profile, not for a linear one).
TEST(LogShear, LinearProfileGivesTheStencilValue)
{
    const Real dz = Real(10);
    const Column col = run_column(uniform_layers(dz), Profile::Linear);
    for (int k = 1; k < nz - 1; ++k) {
        const Real g = Real(0.02) * (Real(2) * dz)
                     / ((col.d[k] + z0) * std::log((col.d[k+1] + z0) / (col.d[k-1] + z0)));
        EXPECT_NEAR(col.central[k],   SmnSmn_of(Real(0.02)), rel_tol() * SmnSmn_of(Real(0.02))) << "k=" << k;
        EXPECT_NEAR(col.corrected[k], SmnSmn_of(g),          rel_tol() * SmnSmn_of(g))          << "k=" << k;
        EXPECT_LT(g, Real(0.02));
    }
}

// Where the wall distance is shorter than the height step to the cell below (d_{k-1} + z0 would not
// be positive), or is not a number, the cell keeps the S_mn S_mn it had, bit for bit.
TEST(LogShear, CellsWithoutPositiveLowerDistanceAreKept)
{
    const Real dz = Real(15.625);
    // The distance of the cells k >= 5 is scaled by 0.1: d(5) = 0.55 dz < dz, d(6) = 0.65 dz, ...
    const Column col = run_column(uniform_layers(dz), Profile::Log, Real(0.1), 5);
    for (int k = 5; k < nz - 1; ++k) { EXPECT_EQ(col.corrected[k], col.central[k]) << "k=" << k; }
    for (int k = 1; k < 5; ++k)      { EXPECT_NE(col.corrected[k], col.central[k]) << "k=" << k; }

    const Column nan_col = run_column(uniform_layers(dz), Profile::Log,
                                      std::numeric_limits<Real>::quiet_NaN(), 5);
    for (int k = 5; k < nz - 1; ++k) { EXPECT_EQ(nan_col.corrected[k], nan_col.central[k]) << "k=" << k; }
    EXPECT_NE(nan_col.corrected[4], nan_col.central[4]);
}

// Fields that differ in x and y (a sloped mesh, u and v changing from face to face, edge strains
// changing from edge to edge): correct_box against a reference written out on the host with each
// point named. An index taken on the wrong face, edge or node changes the result. correct_box
// leaves the wall cell and the top cell of the domain as they were, also in a box that holds
// only part of the column.
TEST(LogShear, CorrectBoxOnFieldsVaryingInXAndY)
{
    const Box cc(IntVect(0,0,0), IntVect(1,1,nz-1));
    const Box gbx = amrex::grow(cc, 1);
    FArrayBox fu(convert(gbx, IntVect(1,0,0)), 1, The_Pinned_Arena());
    FArrayBox fv(convert(gbx, IntVect(0,1,0)), 1, The_Pinned_Arena());
    FArrayBox fz(convert(gbx, IntVect(1,1,1)), 1, The_Pinned_Arena());
    FArrayBox fd(gbx, 1, The_Pinned_Arena());
    FArrayBox fz0(Box(IntVect(-1,-1,0), IntVect(2,2,0)), 1, The_Pinned_Arena());
    FArrayBox f13(convert(gbx, IntVect(1,0,1)), 1, The_Pinned_Arena());
    FArrayBox f23(convert(gbx, IntVect(0,1,1)), 1, The_Pinned_Arena());
    FArrayBox fs(cc, 1, The_Pinned_Arena()), fs0(cc, 1, The_Pinned_Arena()), fpart(cc, 1, The_Pinned_Arena());

    // Analytic fields, each with its own dependence on i and j
    auto U  = [] (int i, int j, int k) { return Real(3) + Real(0.9)*Real(k) + Real(0.31)*Real(i) - Real(0.17)*Real(j) + Real(0.05)*Real(i*k); };
    auto V  = [] (int i, int j, int k) { return Real(-1) + Real(0.4)*Real(k) - Real(0.23)*Real(i) + Real(0.29)*Real(j) + Real(0.03)*Real(j*k); };
    auto Z  = [] (int i, int j, int k) { return Real(k) * (Real(10) + Real(0.7)*Real(i) + Real(0.4)*Real(j)) + Real(0.2)*Real(k*k); };
    auto E13 = [] (int i, int j, int k) { return Real(0.03) + Real(0.004)*Real(i) - Real(0.002)*Real(j) + Real(0.001)*Real(k); };
    auto E23 = [] (int i, int j, int k) { return Real(-0.02) + Real(0.003)*Real(j) + Real(0.005)*Real(i) - Real(0.0007)*Real(k); };
    // the distance exceeds the height step also in the wall cell, so that only correct_box's
    // own limits keep the wall cell (the fallback of corrected_SmnSmn does not act here)
    auto D   = [] (int i, int j, int k) { return Real(12) + (Real(k) + Real(0.5)) * (Real(10) + Real(0.3)*Real(i) + Real(0.6)*Real(j)); };
    auto Z0  = [] (int i, int j) { return Real(0.05) + Real(0.02)*Real(i) + Real(0.04)*Real(j); };
    auto S0  = [] (int i, int j, int k) { return Real(0.01) + Real(0.002)*Real(i) + Real(0.003)*Real(j) + Real(0.0005)*Real(k); };

    auto fill = [] (FArrayBox& f, const std::function<Real(int,int,int)>& g) {
        auto a = f.array();
        const Box& b = f.box();
        for (int k = b.smallEnd(2); k <= b.bigEnd(2); ++k) {
        for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j) {
        for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i) { a(i,j,k) = g(i,j,k); }}}
    };
    fill(fu, U); fill(fv, V); fill(fz, Z); fill(fd, D); fill(f13, E13); fill(f23, E23);
    fill(fz0, [&] (int i, int j, int) { return Z0(i,j); });
    fill(fs, S0); fill(fs0, S0); fill(fpart, S0);

    log_shear::correct_box(cc, cc, fs.array(), fu.const_array(), fv.const_array(), f13.const_array(),
                           f23.const_array(), fd.const_array(), fz0.const_array(), fz.const_array());
    // a box with the lower half of the column only: its top cell is not the top cell of the domain
    Box lower = cc;
    lower.setBig(2, nz/2);
    log_shear::correct_box(lower, cc, fpart.array(), fu.const_array(), fv.const_array(), f13.const_array(),
                           f23.const_array(), fd.const_array(), fz0.const_array(), fz.const_array());
    Gpu::streamSynchronize();

    const auto s = fs.const_array(); const auto part = fpart.const_array();
    int n_changed = 0;
    for (int j = 0; j <= 1; ++j) {
    for (int i = 0; i <= 1; ++i) {
        EXPECT_EQ(s(i,j,0),    S0(i,j,0))    << "wall cell";
        EXPECT_EQ(s(i,j,nz-1), S0(i,j,nz-1)) << "top cell";
        for (int k = 1; k < nz - 1; ++k) {
            // heights of the cell centers k-1, k, k+1: mean of the eight nodes of each cell
            auto zc = [&] (int kk) {
                Real sum = Real(0);
                for (int dk = 0; dk <= 1; ++dk) { for (int dj = 0; dj <= 1; ++dj) { for (int di = 0; di <= 1; ++di) {
                    sum += Z(i+di, j+dj, kk+dk); }}}
                return Real(0.125) * sum;
            };
            const Real dc = D(i,j,k) + Z0(i,j);
            const Real dm = dc - (zc(k) - zc(k-1));
            const Real dp = dc + (zc(k+1) - zc(k));
            ASSERT_GT(dm, Real(0));
            const Real lg = dc * std::log(dp / dm);
            // new gradients: the two faces of the cell, one cell up less one cell down
            const Real gu_new = Real(0.5) * ((U(i,j,k+1) - U(i,j,k-1)) + (U(i+1,j,k+1) - U(i+1,j,k-1))) / lg;
            const Real gv_new = Real(0.5) * ((V(i,j,k+1) - V(i,j,k-1)) + (V(i,j+1,k+1) - V(i,j+1,k-1))) / lg;
            // old gradients: the four edges, each the difference of its two velocities over the
            // height difference of their points (mean of the two nodes above less the two below)
            Real gu_old = Real(0), gv_old = Real(0), e13 = Real(0), e23 = Real(0);
            for (int dk = 0; dk <= 1; ++dk) {
                for (int di = 0; di <= 1; ++di) {
                    const Real h = Real(0.5) * (Z(i+di,j,k+dk+1) + Z(i+di,j+1,k+dk+1))
                                 - Real(0.5) * (Z(i+di,j,k+dk-1) + Z(i+di,j+1,k+dk-1));
                    gu_old += Real(0.25) * (U(i+di,j,k+dk) - U(i+di,j,k+dk-1)) / (Real(0.5) * h);
                    e13    += Real(0.25) * E13(i+di,j,k+dk);
                }
                for (int dj = 0; dj <= 1; ++dj) {
                    const Real h = Real(0.5) * (Z(i,j+dj,k+dk+1) + Z(i+1,j+dj,k+dk+1))
                                 - Real(0.5) * (Z(i,j+dj,k+dk-1) + Z(i+1,j+dj,k+dk-1));
                    gv_old += Real(0.25) * (V(i,j+dj,k+dk) - V(i,j+dj,k+dk-1)) / (Real(0.5) * h);
                    e23    += Real(0.25) * E23(i,j+dj,k+dk);
                }
            }
            const Real n13 = e13 + Real(0.5) * (gu_new - gu_old);
            const Real n23 = e23 + Real(0.5) * (gv_new - gv_old);
            const Real want = std::max(Real(0), S0(i,j,k) + Real(2) * (n13*n13 - e13*e13) + Real(2) * (n23*n23 - e23*e23));
            const Real tol = ((sizeof(Real) == sizeof(double)) ? Real(1e-12) : Real(2e-5)) * (std::abs(want) + S0(i,j,k));
            EXPECT_NEAR(s(i,j,k), want, tol) << "i=" << i << " j=" << j << " k=" << k;
            if (std::abs(s(i,j,k) - S0(i,j,k)) > Real(5e-4) * S0(i,j,k)) { ++n_changed; }
            if (k <= nz/2) { EXPECT_EQ(part(i,j,k), s(i,j,k)) << "partial box, k=" << k; }
            else           { EXPECT_EQ(part(i,j,k), S0(i,j,k)) << "outside the partial box, k=" << k; }
        }
        EXPECT_EQ(part(i,j,0), S0(i,j,0)) << "wall cell, partial box";
    }}
    // the correction is not a round-off change in these fields: more than 25 times the tolerance
    // of the comparison in most cells (it falls with height as 1 / (3 (k + 1/2)^2))
    EXPECT_GT(n_changed, 2 * (nz - 2));
}

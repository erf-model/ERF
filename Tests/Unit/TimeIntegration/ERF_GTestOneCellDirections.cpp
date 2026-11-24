// Unit tests for the time-step estimates on a domain with only one cell in a lateral direction
// (an x-z slice or a single column): the cell size in that direction is arbitrary, so it must
// not enter the anelastic or compressible advective estimate (ERF_TimestepUtils.H) or the
// explicit diffusion rates of the diffusive time-step check (CellDiffusiveRates,
// ERF_TerrainDiffusionLimits.H).

#include <ERF_TimestepUtils.H>
#include <ERF_TerrainDiffusionLimits.H>

#include <gtest/gtest.h>

#include <vector>

using amrex::Real;

namespace {

Real rel_tol ()
{
    return (sizeof(Real) == 8) ? Real(1e-12) : Real(2e-5);
}

amrex::GpuArray<Real,AMREX_SPACEDIM> make_dxinv (Real dx, Real dy, Real dz)
{
    return {Real(1)/dx, Real(1)/dy, Real(1)/dz};
}

} // namespace

TEST(OneCellDirections, AnelasticDropsOneCellDirections)
{
    const amrex::Array4<const Real> no_z{};
    const Real u = Real(3), v = Real(5), w = Real(0.5);
    const Real mf = Real(1);

    // 3-D: both lateral directions count, as before
    const auto dxinv = make_dxinv(Real(100), Real(50), Real(10));
    EXPECT_NEAR(Compute_InvDt_Anelastic(0,0,0, u,v,w, no_z, mf,mf, dxinv, false),
                u/Real(100) + v/Real(50) + w/Real(10), rel_tol());

    // x-z: the result is the same for any dy, and has no v term
    const Real xz = u/Real(100) + w/Real(10);
    for (Real dy : {Real(1e-3), Real(1), Real(1e6)}) {
        const auto d = make_dxinv(Real(100), dy, Real(10));
        EXPECT_NEAR(Compute_InvDt_Anelastic(0,0,0, u,v,w, no_z, mf,mf, d, false, 8, 1),
                    xz, xz*rel_tol());
    }

    // y-z: no u term
    EXPECT_NEAR(Compute_InvDt_Anelastic(0,0,0, u,v,w, no_z, mf,mf, dxinv, false, 1, 8),
                v/Real(50) + w/Real(10), rel_tol());

    // Single column: only the vertical, and zero with no vertical velocity
    for (Real dh : {Real(1e-3), Real(1e6)}) {
        const auto d = make_dxinv(dh, dh, Real(10));
        EXPECT_NEAR(Compute_InvDt_Anelastic(0,0,0, u,v,w, no_z, mf,mf, d, false, 1, 1),
                    w/Real(10), rel_tol());
        EXPECT_EQ(Compute_InvDt_Anelastic(0,0,0, u,v,Real(0), no_z, mf,mf, d, false, 1, 1),
                  Real(0));
    }
}

TEST(OneCellDirections, CompressibleIgnoresOneCellDy)
{
    const amrex::Array4<const Real> no_z{};
    const Real c = Real(340), u = Real(3), v = Real(5), w = Real(0.5);
    const Real mf = Real(1);
    for (bool substep : {false, true}) {
        const auto ref = Compute_InvDt_Compressible(0,0,0, c,u,v,w, no_z, mf,mf,
                                                    make_dxinv(Real(100), Real(1), Real(10)),
                                                    false, substep, 8, 1);
        const auto tiny = Compute_InvDt_Compressible(0,0,0, c,u,v,w, no_z, mf,mf,
                                                     make_dxinv(Real(100), Real(1e-3), Real(10)),
                                                     false, substep, 8, 1);
        EXPECT_NEAR(tiny, ref, ref*rel_tol());
    }
}

TEST(OneCellDirections, DiffusiveRatesIgnoreOneCellDirections)
{
    // Uniform state and diffusivities over the 3x3x3 neighborhood
    const int nk = EddyDiff::NumDiffs;
    std::vector<Real> rho(27, Real(1)), K(27*static_cast<std::size_t>(nk), Real(2));
    const amrex::Dim3 lo{-1,-1,-1};
    const amrex::Array4<const Real> sa(rho.data(), lo, amrex::Dim3{2,2,2}, 1);
    const amrex::Array4<const Real> ka(K.data(),   lo, amrex::Dim3{2,2,2}, nk);
    const amrex::Array4<const Real> no_z{};

    DiffusiveRateSettings set;
    set.has_q = true;
    set.has_ke = true;
    set.has_scalar = true;
    set.mol_mom = Real(0.5);
    set.mol_th  = Real(0.3);
    set.mol_c   = Real(0.4);

    const Real dx = Real(10), dz = Real(5);
    const Real a = Real(1)/(dx*dx), dz2 = Real(1)/(dz*dz);

    // x-z: the rates are those with b = 0, whatever dy is
    set.y_active = false;
    const Real rm_xz = MomentumDiffusiveRate(Real(1), Real(2.5), Real(2.5), a, Real(0), Real(0), dz2, Real(1));
    const Real rs_xz = ScalarDiffusiveRate(Real(1), Real(2.4), Real(2.4), a, Real(0), Real(0), dz2, Real(1));
    for (Real dy : {Real(1e-3), Real(1), Real(1e6)}) {
        Real rm, rs;
        CellDiffusiveRates(0,0,0, sa, ka, no_z, Real(1), Real(1), make_dxinv(dx, dy, dz), set, rm, rs);
        EXPECT_NEAR(rm, rm_xz, rm_xz*rel_tol());
        EXPECT_NEAR(rs, rs_xz, rs_xz*rel_tol());
    }

    // Single column: vertical only
    set.x_active = false;
    const Real rm_z = MomentumDiffusiveRate(Real(1), Real(2.5), Real(2.5), Real(0), Real(0), Real(0), dz2, Real(1));
    const Real rs_z = ScalarDiffusiveRate(Real(1), Real(2.4), Real(2.4), Real(0), Real(0), Real(0), dz2, Real(1));
    for (Real dh : {Real(1e-3), Real(1e6)}) {
        Real rm, rs;
        CellDiffusiveRates(0,0,0, sa, ka, no_z, Real(1), Real(1), make_dxinv(dh, dh, dz), set, rm, rs);
        EXPECT_NEAR(rm, rm_z, rm_z*rel_tol());
        EXPECT_NEAR(rs, rs_z, rs_z*rel_tol());
        EXPECT_NEAR(rm, Real(2) * Real(2.5) * dz2, rm_z*rel_tol());
    }

    // Both directions active (the default): unchanged
    set.x_active = true; set.y_active = true;
    const Real dy = Real(20), b = Real(1)/(dy*dy);
    Real rm, rs;
    CellDiffusiveRates(0,0,0, sa, ka, no_z, Real(1), Real(1), make_dxinv(dx, dy, dz), set, rm, rs);
    const Real rm_3d = MomentumDiffusiveRate(Real(1), Real(2.5), Real(2.5), a, b, Real(0), dz2, Real(1));
    EXPECT_NEAR(rm, rm_3d, rm_3d*rel_tol());
}

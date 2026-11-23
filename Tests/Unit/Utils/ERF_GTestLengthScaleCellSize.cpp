// Unit tests for LengthScaleCellSize (ERF_LengthScaleCellSize.H), the cell sizes ERF uses where
// a lateral cell size is a length scale rather than the denominator of a difference:
//  * a 3-D domain keeps dx and dy; an x-z slice uses dx for dy, a y-z slice dy for dx; a
//    single column keeps the dx and dy of the inputs;
//  * ComputeTurbulentViscosityLES on an x-z slice gives the same eddy viscosities for any dy,
//    equal to those of a 3-D domain with dy = dx under the same uniform strain, for isotropic
//    Smagorinsky and for Smagorinsky-2D.

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <ERF_EddyViscosity.H>
#include <ERF_LengthScaleCellSize.H>

#include <gtest/gtest.h>

#include <memory>

using amrex::Real;

namespace {

Real rel_tol ()
{
    return (sizeof(Real) == 8) ? Real(1e-12) : Real(2e-5);
}

amrex::Geometry make_geom (int nx, int ny, int nz, Real dx, Real dy, Real dz)
{
    amrex::Box domain(amrex::IntVect(0), amrex::IntVect(nx-1, ny-1, nz-1));
    amrex::RealBox rb({Real(0), Real(0), Real(0)}, {nx*dx, ny*dy, nz*dz});
    amrex::Array<int,AMREX_SPACEDIM> periodic{1, 1, 0};
    return amrex::Geometry(domain, rb, amrex::CoordSys::cartesian, periodic);
}

} // namespace

TEST(LengthScaleCellSize, SubstitutesOnlyAOneCellDirection)
{
    const Real dx = Real(100), dy = Real(40), dz = Real(10);
    const Real tol = Real(10)*rel_tol();

    // 3-D: unchanged
    auto len = LengthScaleCellSize(make_geom(8, 4, 4, dx, dy, dz));
    EXPECT_NEAR(len[0], dx, dx*tol);
    EXPECT_NEAR(len[1], dy, dy*tol);
    EXPECT_NEAR(len[2], dz, dz*tol);

    // x-z: dy -> dx
    len = LengthScaleCellSize(make_geom(8, 1, 4, dx, dy, dz));
    EXPECT_NEAR(len[0], dx, dx*tol);
    EXPECT_NEAR(len[1], dx, dx*tol);
    EXPECT_NEAR(len[2], dz, dz*tol);

    // y-z: dx -> dy
    len = LengthScaleCellSize(make_geom(1, 8, 4, dx, dy, dz));
    EXPECT_NEAR(len[0], dy, dy*tol);
    EXPECT_NEAR(len[1], dy, dy*tol);

    // Single column: the nominal column size of the inputs
    len = LengthScaleCellSize(make_geom(1, 1, 4, dx, dy, dz));
    EXPECT_NEAR(len[0], dx, dx*tol);
    EXPECT_NEAR(len[1], dy, dy*tol);

    // The inverse matches
    const auto inv = LengthScaleInvCellSize(make_geom(8, 1, 4, dx, dy, dz));
    EXPECT_NEAR(inv[0], Real(1)/dx, tol/dx);
    EXPECT_NEAR(inv[1], Real(1)/dx, tol/dx);
    EXPECT_NEAR(inv[2], Real(1)/dz, tol/dz);
}

namespace {

struct LESValues { Real mom_h, mom_v; };

// Eddy viscosities at cell (0,0,1) under a uniform strain, Smagorinsky with or without smag2d
LESValues run_les (int ny, Real dy, bool smag2d)
{
    const Real dx = Real(100), dz = Real(20), rho = Real(1.2);
    const amrex::Geometry geom = make_geom(8, ny, 4, dx, dy, dz);
    amrex::BoxArray ba(geom.Domain());
    amrex::DistributionMapping dm(ba);

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> tau(9);
    for (int c = 0; c < 9; ++c) {
        tau[c] = std::make_unique<amrex::MultiFab>(ba, dm, 1, 2);
        tau[c]->setVal(Real(0));
    }
    tau[TauType::tau11]->setVal(Real(2e-3));
    tau[TauType::tau33]->setVal(Real(-2e-3));
    tau[TauType::tau13]->setVal(Real(1e-3));

    amrex::MultiFab cons(ba, dm, RhoTheta_comp+1, 1);
    cons.setVal(rho, Rho_comp, 1, 1);
    cons.setVal(rho*Real(300), RhoTheta_comp, 1, 1);
    amrex::MultiFab K(ba, dm, EddyDiff::NumDiffs, 1);
    K.setVal(Real(0));
    amrex::MultiFab diss(ba, dm, 1, 0);

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> mapfac(MapFacType::num);
    for (auto& m : mapfac) {
        m = std::make_unique<amrex::MultiFab>(ba, dm, 1, 2);
        m->setVal(Real(1));
    }
    auto z_nd = std::make_unique<amrex::MultiFab>(amrex::convert(ba, amrex::IntVect(1)), dm, 1, 1);
    z_nd->setVal(Real(0));

    TurbChoice tc;
    tc.les_type = LESType::Smagorinsky;
    tc.smag2d = smag2d;
    tc.mix_isotropic = !smag2d;
    tc.Cs = Real(0.1);

    std::unique_ptr<SurfaceLayer> no_sl;
    MoistureComponentIndices mi;
    ComputeTurbulentViscosityLES(tau, cons, K, diss, geom, false, mapfac, z_nd, tc,
                                 Real(9.81), no_sl, mi, nullptr, nullptr);

    amrex::FArrayBox host(K[0].box(), EddyDiff::NumDiffs, amrex::The_Pinned_Arena());
    host.copy<amrex::RunOn::Device>(K[0]);
    amrex::Gpu::streamSynchronize();
    const auto h = host.const_array();
    return {h(0,0,1,EddyDiff::Mom_h), h(0,0,1,EddyDiff::Mom_v)};
}

} // namespace

TEST(LengthScaleCellSize, LESOnAnXZSliceIgnoresDy)
{
    const Real tol = Real(10)*rel_tol();
    for (bool smag2d : {false, true}) {
        // 3-D reference with dy = dx: the same strain has no y-dependence
        const LESValues ref = run_les(4, Real(100), smag2d);
        EXPECT_GT(ref.mom_h, Real(0));
        for (Real dy : {Real(1), Real(100), Real(5000)}) {
            const LESValues v = run_les(1, dy, smag2d);
            EXPECT_NEAR(v.mom_h, ref.mom_h, ref.mom_h*tol) << "smag2d " << smag2d << " dy " << dy;
            EXPECT_NEAR(v.mom_v, ref.mom_v, ref.mom_v*tol + tol) << "smag2d " << smag2d << " dy " << dy;
        }
    }
}

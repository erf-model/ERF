#include <cmath>
#include <limits>
#include <memory>

#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_Diffusion.H"
#include "ERF_EddyViscosity.H"
#include "ERF_IndexDefines.H"

// Motivation: the fraction-stress ("form 7") immersed wall law applies the wall stress as a source
// in the wall cell, so the bottom face of that cell (the face on the solid) must carry no resolved
// flux. Before the wall-face masks, the implicit vertical solves coupled the wall cell to the solid
// below: theta leaked into the solid and momentum was dragged a second time. The partial wall cells
// also dissipated TKE with a length measured at the fluid centroid, so a nearly solid cell dissipated
// several times faster than flat ground's first cell (erf.if_flat_wall_dissipation). The implicit
// solves must keep the heat and momentum above a closed wall face, and a partial wall cell must
// dissipate like flat ground's first cell while its viscosity keeps the centroid length.

using namespace amrex;

namespace {

constexpr int nz = 12;
constexpr int kwall = 4;          // the wall cell; cells below are solid
constexpr Real dz = 1.0;
constexpr Real Kdiff = 5.0;
const double tol = 1000. * std::numeric_limits<Real>::epsilon();

Vector<BCRec> foextrap_bcs ()
{
    Vector<BCRec> bcs(BCVars::NumTypes);   // the scalars, then the three velocities
    for (auto& bc : bcs) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            bc.setLo(dir, ERFBCType::foextrap);
            bc.setHi(dir, ERFBCType::foextrap);
        }
    }
    return bcs;
}

SolverChoice turbulent_choice ()
{
    SolverChoice sc;
    sc.diffChoice.molec_diff_type = MolecDiffType::None;
    sc.turbChoice.resize(1);
    sc.turbChoice[0].use_kturb = true;
    sc.turbChoice[0].use_keqn  = true;
    sc.turbChoice[0].les_type  = LESType::Deardorff;
    return sc;
}

// value of component comp at (i,j,k), read on the host
Real host_value (MultiFab const& mf, int i, int j, int k, int comp = 0)
{
    MultiFab h(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrowVect(),
               MFInfo().SetArena(The_Pinned_Arena()));
    MultiFab::Copy(h, mf, 0, 0, mf.nComp(), mf.nGrowVect());
    Gpu::streamSynchronize();
    return h[0](IntVect(i,j,k), comp);
}

void set_cells (MultiFab& mf, int comp, Real solid, Real base, Real slope)
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto const& a = mf.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i,j,k,comp) = (k < kwall) ? solid : base + slope * static_cast<Real>(k - kwall);
        });
    }
    Gpu::streamSynchronize();
}

// ones, with zero on the z-face (or z-edge) k = kwall when closed
void set_wall_face (MultiFab& mf, bool closed)
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto const& a = mf.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i,j,k) = (closed && k == kwall) ? Real(0.) : Real(1.);
        });
    }
    Gpu::streamSynchronize();
}

// theta summed over cells [klo, khi] of the column, and the largest change of a fluid cell
struct StateResult { Real sum_fluid, sum_solid, max_change; };

StateResult implicit_theta (bool closed)
{
    const Box dom(IntVect(0), IntVect(0, 0, nz-1));
    BoxArray ba(dom);
    DistributionMapping dm(ba);
    MultiFab cons(ba, dm, NVAR_max, 1);
    cons.setVal(Real(0.));
    cons.setVal(Real(1.), Rho_comp, 1, 1);
    set_cells(cons, RhoTheta_comp, Real(280.), Real(300.), Real(0.5));
    MultiFab mu(ba, dm, EddyDiff::NumDiffs, 1);
    mu.setVal(Kdiff);
    MultiFab hfx(convert(ba, IntVect(0,0,1)), dm, 1, 1);
    hfx.setVal(Real(0.));
    MultiFab wf(convert(ba, IntVect(0,0,1)), dm, 1, 1);
    set_wall_face(wf, closed);

    auto bcs = foextrap_bcs();
    SolverChoice sc = turbulent_choice();
    GpuArray<Real, AMREX_SPACEDIM*2> neumann{};
    const GpuArray<Real, AMREX_SPACEDIM> dxinv{Real(1.), Real(1.), Real(1.)/dz};
    ImplicitDiffForStateLU_N(dom, dom, 0, RhoTheta_comp, 10.0, neumann, cons[0].array(), dxinv,
                             hfx[0].const_array(), mu[0].const_array(), sc, bcs.data(),
                             false, Real(1.), false, wf[0].const_array());
    Gpu::streamSynchronize();

    StateResult r{0., 0., 0.};
    for (int k = 0; k < nz; ++k) {
        const Real th = host_value(cons, 0, 0, k, RhoTheta_comp);
        if (k < kwall) { r.sum_solid += th; }
        else {
            r.sum_fluid += th;
            r.max_change = amrex::max(r.max_change, std::abs(th - (Real(300.) + Real(0.5) * (k - kwall))));
        }
    }
    return r;
}

struct MomResult { Real sum_fluid, sum_solid, max_change; };

MomResult implicit_xmom (bool closed)
{
    const Box dom(IntVect(0), IntVect(0, 0, nz-1));
    BoxArray ba(dom);
    DistributionMapping dm(ba);
    MultiFab cons(ba, dm, NVAR_max, 1);
    cons.setVal(Real(0.));
    cons.setVal(Real(1.), Rho_comp, 1, 1);
    MultiFab mu(ba, dm, EddyDiff::NumDiffs, 1);
    mu.setVal(Kdiff);
    BoxArray bax = convert(ba, IntVect(1,0,0));
    MultiFab u(bax, dm, 1, 0);
    set_cells(u, 0, Real(0.), Real(5.), Real(0.5));
    BoxArray baxz = convert(ba, IntVect(1,0,1));
    MultiFab tau(baxz, dm, 1, 1), tau_corr(baxz, dm, 1, 1), wf(baxz, dm, 1, 1);
    tau.setVal(Real(0.));
    tau_corr.setVal(Real(0.));
    set_wall_face(wf, closed);
    // every cell column spans the domain
    iMultiFab kext(BoxArray(makeSlab(dom, 2, 0)), dm, 2, IntVect(1,1,0));
    kext.setVal(0, 0, 1, IntVect(1,1,0));
    kext.setVal(nz-1, 1, 1, IntVect(1,1,0));

    auto bcs = foextrap_bcs();
    SolverChoice sc = turbulent_choice();
    const GpuArray<Real, AMREX_SPACEDIM> dxinv{Real(1.), Real(1.), Real(1.)/dz};
    ImplicitDiffForMomLU_N<0>(bax[0], dom, 0, 10.0, kext[0].const_array(), cons[0].const_array(),
                              u[0].array(), tau[0].const_array(), tau_corr[0].const_array(), dxinv,
                              mu[0].const_array(), sc, bcs.data(), false, Real(1.), false,
                              wf[0].const_array());
    Gpu::streamSynchronize();

    MomResult r{0., 0., 0.};
    for (int k = 0; k < nz; ++k) {
        for (int i = 0; i <= 1; ++i) {
            const Real v = host_value(u, i, 0, k);
            if (k < kwall) { r.sum_solid += v; }
            else {
                r.sum_fluid += v;
                r.max_change = amrex::max(r.max_change, std::abs(v - (Real(5.) + Real(0.5) * (k - kwall))));
            }
        }
    }
    return r;
}

// RANS dissipation and viscosity in a column: solid below kwall, a wall cell with solid fraction
// bwall at kwall, uniform theta (no stratification) and TKE; wdist at the wall cell = dwall
struct RansResult { Real diss, mu_v; };

RansResult rans_wall_cell (bool immersed, bool flat_wall_diss, Real bwall, Real dwall)
{
    const Box dom(IntVect(0), IntVect(0, 0, nz-1));
    RealBox rb({0., 0., 0.}, {dz, dz, nz*dz});
    Array<int,3> per{1, 1, 0};
    Geometry geom(dom, rb, 0, per);
    BoxArray ba(dom);
    DistributionMapping dm(ba);
    MultiFab cons(ba, dm, NVAR_max, 1);
    cons.setVal(Real(0.));
    cons.setVal(Real(1.), Rho_comp, 1, 1);
    cons.setVal(Real(300.), RhoTheta_comp, 1, 1);
    cons.setVal(Real(0.5), RhoKE_comp, 1, 1);
    MultiFab wdist(ba, dm, 1, 1);
    for (MFIter mfi(wdist); mfi.isValid(); ++mfi) {
        auto const& d = wdist.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            d(i,j,k) = (k == kwall) ? dwall : Real(1.) + amrex::max(k - kwall, 0) * dz;
        });
    }
    MultiFab tb(ba, dm, 1, 1);
    for (MFIter mfi(tb); mfi.isValid(); ++mfi) {
        auto const& a = tb.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i,j,k) = (k < kwall) ? Real(1.) : ((k == kwall) ? bwall : Real(0.));
        });
    }
    MultiFab mu(ba, dm, EddyDiff::NumDiffs, 1), diss(ba, dm, 1, 1);
    mu.setVal(Real(0.));
    diss.setVal(Real(0.));
    auto z_nd = std::make_unique<MultiFab>(convert(ba, IntVect(1)), dm, 1, 1);
    z_nd->setVal(Real(0.));
    TurbChoice tc;
    tc.rans_type = RANSType::kEqn;
    std::unique_ptr<SurfaceLayer> no_sfc;
    Gpu::streamSynchronize();

    Vector<std::unique_ptr<MultiFab>> no_tau;   // the realizability cap is off (rans_realizability = 0)
    const CanopyTurbulence no_canopy{};
    ComputeTurbulentViscosityRANS(0, no_tau, cons, wdist, mu, diss, geom, false, z_nd, tc, Real(9.81),
                                  no_sfc, nullptr, no_canopy, nullptr, immersed ? &tb : nullptr, flat_wall_diss,
                                  Real(0.005));
    Gpu::streamSynchronize();
    return {host_value(diss, 0, 0, kwall), host_value(mu, 0, 0, kwall, EddyDiff::Mom_v)};
}

} // namespace

TEST(ImmersedWallFace, ImplicitThetaSolveKeepsTheHeatAboveAClosedWallFace)
{
    const StateResult closed = implicit_theta(true);
    const StateResult open   = implicit_theta(false);
    // the solve did diffuse the fluid column
    ASSERT_GT(closed.max_change, 0.01);
    // closed: the fluid cells keep their heat and the solid keeps its own
    EXPECT_NEAR(closed.sum_fluid, 8 * 300. + 0.5 * 28., tol * 3000.);
    EXPECT_NEAR(closed.sum_solid, 4 * 280., tol * 3000.);
    // open: heat leaks into the solid (else this test shows nothing)
    EXPECT_LT(open.sum_fluid, closed.sum_fluid - 1.);
}

TEST(ImmersedWallFace, ImplicitMomentumSolveKeepsTheMomentumAboveAClosedWallFace)
{
    const MomResult closed = implicit_xmom(true);
    const MomResult open   = implicit_xmom(false);
    ASSERT_GT(closed.max_change, 0.01);
    // two x-faces per height, u = 5 + 0.5 (k - kwall) above the wall
    EXPECT_NEAR(closed.sum_fluid, 2. * (8 * 5. + 0.5 * 28.), tol * 100.);
    EXPECT_NEAR(closed.sum_solid, 0., tol * 100.);
    EXPECT_LT(open.sum_fluid, closed.sum_fluid - 0.1);
}

TEST(ImmersedWallFace, PartialWallCellDissipatesLikeFlatGroundsFirstCell)
{
    const Real b = 0.9;
    const Real d_centroid = (1. - b) * dz / 2.;   // 0.05: the wall cell's length scale
    const Real d_flat     = dz / 2.;              // flat ground's first cell

    const RansResult on   = rans_wall_cell(true,  true,  b, d_centroid);
    const RansResult off  = rans_wall_cell(true,  false, b, d_centroid);
    const RansResult flat = rans_wall_cell(false, false, 0., d_flat);

    // the dissipation of flat ground's first cell ...
    EXPECT_NEAR(on.diss, flat.diss, tol * flat.diss);
    // ... far below what the centroid length gives (else this test shows nothing)
    EXPECT_GT(off.diss, 2. * on.diss);
    // the viscosity keeps the centroid length
    EXPECT_NEAR(on.mu_v, off.mu_v, tol * off.mu_v);
    EXPECT_LT(on.mu_v, flat.mu_v);
}

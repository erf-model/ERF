// Warm bubble of the ABL problem (prob.bubble_*), checked against the closed form
// theta' = bubble_dT cos^2(pi r / 2) on a stretched, tilted mesh, where the
// cell-centre height z_cc differs from (k + 1/2) dz in every cell.
//
// The bubble is WRF's idealized-fire start-up bubble (module_initialize_fire.F):
// theta is raised at constant pressure and the density recomputed. Under a
// fixed-in-time density (anelastic) only rho*theta can carry it.

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>

#include <array>
#include <cmath>
#include <limits>
#include <string>

#include <gtest/gtest.h>

#include "../../../Exec/ERF_Prob.H"
#include "../../../Source/DataStructs/ERF_DataStruct.H"
#include "../../../Source/Utils/ERF_EOS.H"
#include "ERF_Constants.H"
#include "ERF_IndexDefines.H"

using amrex::Real;

namespace {

constexpr int n_cells = 16;
constexpr Real dxy = Real(10.0);

// Bubble used by every test, off the domain centre and with three different radii
constexpr Real b_dT = Real(0.5);
constexpr Real b_xc = Real(80.0);
constexpr Real b_yc = Real(75.0);
constexpr Real b_zc = Real(30.0);
constexpr Real b_xr = Real(50.0);
constexpr Real b_yr = Real(40.0);
constexpr Real b_zr = Real(25.0);

// Node heights: stretched in k (dz from 2.25 m to 9.75 m) on a ground tilted in x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real node_height (int i, int /*j*/, int k)
{
    return Real(2.0) * Real(k) + Real(0.25) * Real(k) * Real(k) + Real(0.5) * Real(i);
}

// Cell-centre height of the mesh above, as the mean of the cell's eight nodes
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real cell_height (int i, int j, int k)
{
    Real s = Real(0.0);
    for (int di = 0; di < 2; ++di) {
        for (int dj = 0; dj < 2; ++dj) {
            for (int dk = 0; dk < 2; ++dk) {
                s += node_height(i + di, j + dj, k + dk);
            }
        }
    }
    return Real(0.125) * s;
}

// The closed form, written out independently of the code under test
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real expected_dtheta (int i, int j, int k)
{
    const Real pi = Real(3.14159265358979323846);
    const Real x  = (Real(i) + Real(0.5)) * dxy;
    const Real y  = (Real(j) + Real(0.5)) * dxy;
    const Real z  = cell_height(i, j, k);
    const Real rx = (x - b_xc) / b_xr;
    const Real ry = (y - b_yc) / b_yr;
    const Real rz = (z - b_zc) / b_zr;
    const Real r  = std::sqrt(rx*rx + ry*ry + rz*rz);
    if (r > Real(1.0)) { return Real(0.0); }
    const Real c = std::cos(Real(0.5) * pi * r);
    return b_dT * c * c;
}

// Problem::phys_zlo/phys_zhi read the vertical levels on any mesh other than
// ConstantDz, and SolverChoice::mesh_type is static: a test run earlier (in a
// shuffled order) can leave it at VariableDz or StretchedDz, and the bubble then
// read the empty m_zlevels_stag. The guard pins ConstantDz and restores the old type.
// run_bubble sets the static SolverChoice::terrain_type, which a later init_params
// turns into mesh_type = VariableDz, so the guard restores that one as well.
struct ScopedBubbleParams
{
    MeshType    old_mesh_type    = SolverChoice::mesh_type;
    TerrainType old_terrain_type = SolverChoice::terrain_type;

    ScopedBubbleParams ()
    {
        SolverChoice::set_mesh_type(MeshType::ConstantDz);
        amrex::ParmParse pp_erf("erf");
        pp_erf.add("prob_name", std::string("ABL"));
        amrex::ParmParse pp("prob");
        pp.add("bubble_dT", b_dT);
        pp.add("bubble_x_c", b_xc);
        pp.add("bubble_y_c", b_yc);
        pp.add("bubble_z_c", b_zc);
        pp.add("bubble_x_r", b_xr);
        pp.add("bubble_y_r", b_yr);
        pp.add("bubble_z_r", b_zr);
    }

    ~ScopedBubbleParams ()
    {
        amrex::ParmParse pp_erf("erf");
        pp_erf.remove("prob_name");
        amrex::ParmParse pp("prob");
        for (const char* key : {"bubble_dT", "bubble_x_c", "bubble_y_c", "bubble_z_c",
                                "bubble_x_r", "bubble_y_r", "bubble_z_r"}) {
            pp.remove(key);
        }
        SolverChoice::set_mesh_type(old_mesh_type);
        SolverChoice::terrain_type = old_terrain_type;
    }
};

struct BubbleErrors
{
    Real theta = Real(0.0);    // max |theta_new - theta - dtheta_expected| [K]
    Real pressure = Real(0.0); // max |p_new - p| / p
    Real qv = Real(0.0);       // max |qv_new - qv|
    Real rho = Real(0.0);      // max |rho perturbation| [kg/m^3]
    Real max_dtheta = Real(0.0);
    int  n_inside = 0;
};

// Builds the background state, runs Problem::init_custom_pert for the ABL
// problem, and compares the perturbed state with the closed form. A free
// function so the device lambdas are not inside a TEST body.
BubbleErrors
run_bubble (bool fixed_density, bool moist)
{
    const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cells - 1));
    const amrex::RealBox real_box({AMREX_D_DECL(Real(0.0), Real(0.0), Real(0.0))},
                                  {AMREX_D_DECL(n_cells * dxy, n_cells * dxy, Real(100.0))});
    const std::array<int, AMREX_SPACEDIM> periodic{AMREX_D_DECL(1, 1, 0)};
    const amrex::Geometry geom(domain, &real_box, 0, periodic.data());
    const amrex::BoxArray ba(domain);
    const amrex::DistributionMapping dm(ba);
    const int ncomp = RhoQ1_comp + 1;

    amrex::MultiFab state(ba, dm, ncomp, 0);
    amrex::MultiFab pert(ba, dm, ncomp, 0);
    amrex::MultiFab r_hse(ba, dm, 1, 0);
    amrex::MultiFab p_hse(ba, dm, 1, 0);
    amrex::MultiFab z_cc(ba, dm, 1, 0);
    amrex::MultiFab mf_m(ba, dm, 1, 0);
    amrex::MultiFab z_nd(amrex::convert(ba, amrex::IntVect(1)), dm, 1, 0);
    // The ground under each column: the bottom-face nodes, as ERF_InitCustomPertState.cpp
    // builds z_sfc_nd (erf-model/ERF#4843)
    amrex::BoxList bl_sfc(amrex::IndexType::TheNodeType());
    for (int ib = 0; ib < ba.size(); ++ib) {
        amrex::Box b = amrex::surroundingNodes(ba[ib]);
        b.setRange(2, domain.smallEnd(2));
        bl_sfc.push_back(b);
    }
    amrex::MultiFab z_sfc_nd(amrex::BoxArray(std::move(bl_sfc)), dm, 1, 0);
    pert.setVal(Real(0.0));
    mf_m.setVal(Real(1.0));

    const Real qv0 = moist ? Real(0.008) : Real(0.0);
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto st = state.array(mfi);
        const auto rh = r_hse.array(mfi);
        const auto ph = p_hse.array(mfi);
        const auto zc = z_cc.array(mfi);
        const auto zn = z_nd.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            const Real z   = cell_height(i, j, k);
            const Real rho = Real(1.2) * std::exp(-z / Real(8000.0));
            const Real th  = Real(300.0) + Real(0.003) * z;
            st(i,j,k,Rho_comp)       = rho;
            st(i,j,k,RhoTheta_comp)  = rho * th;
            st(i,j,k,RhoKE_comp)     = Real(0.0);
            st(i,j,k,RhoScalar_comp) = Real(0.0);
            st(i,j,k,RhoQ1_comp)     = rho * qv0;
            rh(i,j,k) = rho;
            ph(i,j,k) = getPgivenRTh(rho * th, qv0);
            zc(i,j,k) = z;
        });
        amrex::ParallelFor(amrex::convert(bx, amrex::IntVect(1)),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            zn(i,j,k) = node_height(i, j, k);
        });
        const auto zs = z_sfc_nd.array(mfi.index());
        amrex::ParallelFor(z_sfc_nd.boxArray()[mfi.index()],
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            zs(i,j,k) = node_height(i, j, k);
        });
    }
    amrex::Gpu::streamSynchronize();

    SolverChoice sc{};
    sc.moisture_type = moist ? MoistureType::Kessler : MoistureType::None;
    sc.rdOcp = R_d / Cp_d;
    sc.terrain_type = TerrainType::StaticFittedMesh;
    sc.fixed_density.assign(1, fixed_density ? 1 : 0);

    Problem prob(geom.ProbLo(), geom.ProbHi());
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        prob.init_custom_pert(mfi.validbox(), state.const_array(mfi), pert.array(mfi),
                              r_hse.array(mfi), p_hse.array(mfi), z_nd.const_array(mfi),
                              z_cc.const_array(mfi), z_sfc_nd.const_array(mfi.index()),
                              geom.data(), mf_m.const_array(mfi), sc, 0);
    }
    amrex::Gpu::streamSynchronize();

    // ERF_InitCustomPertState.cpp drops the density perturbation under a fixed density
    const bool add_rho = !fixed_density;
    amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax, amrex::ReduceOpMax,
                     amrex::ReduceOpMax, amrex::ReduceOpMax, amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<Real, Real, Real, Real, Real, int> reduce_data(reduce_op);
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        const auto s = state.const_array(mfi);
        const auto p = pert.const_array(mfi);
        reduce_op.eval(mfi.validbox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
                -> amrex::GpuTuple<Real, Real, Real, Real, Real, int>
            {
                const Real rho    = s(i,j,k,Rho_comp);
                const Real rt     = s(i,j,k,RhoTheta_comp);
                const Real qv     = s(i,j,k,RhoQ1_comp) / rho;
                const Real rho_n  = rho + (add_rho ? p(i,j,k,Rho_comp) : Real(0.0));
                const Real rt_n   = rt + p(i,j,k,RhoTheta_comp);
                const Real qv_n   = (s(i,j,k,RhoQ1_comp) + p(i,j,k,RhoQ1_comp)) / rho_n;
                const Real dth    = expected_dtheta(i, j, k);
                const Real p_old  = getPgivenRTh(rt, qv);
                const Real p_new  = getPgivenRTh(rt_n, qv_n);
                return {std::abs(rt_n / rho_n - rt / rho - dth),
                        std::abs(p_new - p_old) / p_old,
                        std::abs(qv_n - qv),
                        std::abs(p(i,j,k,Rho_comp)),
                        dth,
                        (dth > Real(0.0)) ? 1 : 0};
            });
    }
    const auto v = reduce_data.value();
    BubbleErrors e;
    e.theta      = amrex::get<0>(v);
    e.pressure   = amrex::get<1>(v);
    e.qv         = amrex::get<2>(v);
    e.rho        = amrex::get<3>(v);
    e.max_dtheta = amrex::get<4>(v);
    e.n_inside   = amrex::get<5>(v);
    amrex::ParallelDescriptor::ReduceRealMax({e.theta, e.pressure, e.qv, e.rho, e.max_dtheta});
    amrex::ParallelDescriptor::ReduceIntSum(e.n_inside);
    return e;
}

Real
tol (Real scale)
{
    return Real(256.0) * std::numeric_limits<Real>::epsilon() * scale;
}

} // namespace

// Motivation: the compressible bubble is WRF's: theta + dT cos^2(pi r/2) at each
// cell's own height z_cc, at unchanged pressure (a density deficit, as WRF
// recomputes alt from the unchanged p). The first form of this option raised
// rho*theta at fixed density, which leaves theta right but raises the pressure
// by gamma*dT/theta (0.23 % = 230 Pa at the centre) and starts the run with an
// acoustic pulse and no buoyancy; it fails the pressure check by 2.2e-3. A bubble
// placed on (k + 1/2) dz instead of z_cc fails the theta check by 0.36 K.
TEST(ABLWarmBubble, CompressibleBubbleIsIsobaricAndMatchesClosedForm)
{
    ScopedBubbleParams params;
    const BubbleErrors e = run_bubble(false, false);
    EXPECT_GT(e.n_inside, 50);
    EXPECT_GT(e.max_dtheta, Real(0.4));
    EXPECT_LE(e.theta, tol(Real(300.0)));
    EXPECT_LE(e.pressure, tol(Real(1.0)));
    EXPECT_GT(e.rho, Real(1.e-4));  // the bubble is a density deficit
}

// Motivation: with moisture the bubble must keep the mixing ratio qv, which the
// pressure depends on, while the density drops; RhoQ1 is rescaled with it.
TEST(ABLWarmBubble, MoistBubbleKeepsMixingRatioAndPressure)
{
    ScopedBubbleParams params;
    const BubbleErrors e = run_bubble(false, true);
    EXPECT_GT(e.n_inside, 50);
    EXPECT_LE(e.theta, tol(Real(300.0)));
    EXPECT_LE(e.pressure, tol(Real(1.0)));
    EXPECT_LE(e.qv, tol(Real(0.01)));
}

// Motivation: under a fixed-in-time density (anelastic) the density perturbation is
// dropped, so the constant-pressure form would leave no bubble at all; there the
// bubble is carried by rho*theta at the unchanged density.
TEST(ABLWarmBubble, FixedDensityBubbleIsCarriedByRhoTheta)
{
    ScopedBubbleParams params;
    const BubbleErrors e = run_bubble(true, false);
    EXPECT_GT(e.n_inside, 50);
    EXPECT_GT(e.max_dtheta, Real(0.4));
    EXPECT_LE(e.theta, tol(Real(300.0)));
    EXPECT_EQ(e.rho, Real(0.0));
}

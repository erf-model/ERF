#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_AllRegular.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <ERF_Advection.H>
#include <ERF_NumericalConstants.H>

#include <gtest/gtest.h>

#include "../ERF_GTestAssertions.H"

#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

using amrex::Array4;
using amrex::Box;
using amrex::FArrayBox;
using amrex::GpuArray;
using amrex::IntVect;
using amrex::Real;

namespace {

constexpr int kNx = 12;
constexpr int kNy = 10;
constexpr int kNz = 12;
// Enough ghost cells for the sixth-order stencils on the grown EB flux boxes
constexpr int kNGrow = 5;

// Smooth, non-uniform fields so every stencil sees different values. The map factors
// are anisotropic (m_x != m_y on every face), which is what tells the two choices of
// map factor apart: with m_x == m_y, as in every map projection ERF builds today, the
// x-face mass flux is the same whether it is divided by m_y or by m_x.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real field (const int i, const int j, const int k, const Real a, const Real b, const Real c) noexcept
{
    return a + Real(0.3)*std::sin(Real(0.37)*i + b) * std::cos(Real(0.29)*j - c)
             + Real(0.2)*std::cos(Real(0.41)*k + b*c);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mapfac (const int i, const int j, const Real a, const Real b) noexcept
{
    return a + Real(0.1)*std::sin(Real(0.23)*i + b) + Real(0.07)*std::cos(Real(0.31)*j - b);
}

void fill (FArrayBox& fab, const Real a, const Real b, const Real c)
{
    auto arr = fab.array();
    amrex::ParallelFor(fab.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        arr(i,j,k) = field(i,j,k,a,b,c);
    });
}

void fill_mapfac (FArrayBox& fab, const Real a, const Real b)
{
    auto arr = fab.array();
    amrex::ParallelFor(fab.box(), [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        arr(i,j,0) = mapfac(i,j,a,b);
    });
}

Box mapfac_box (const Box& domain, const IntVect& nodal)
{
    Box b = amrex::convert(domain, nodal);
    b.setRange(2, 0);
    return b.grow(IntVect(kNGrow, kNGrow, 0));
}

// The face boxes grown as erf_slow_rhs_pre grows them for EB
amrex::Vector<Box> grown_flux_boxes (const Box& tb)
{
    amrex::Vector<Box> out(AMREX_SPACEDIM);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        IntVect iv(1, 1, 1);
        iv[dir] = 0;
        Box b = tb;
        out[dir] = (b.growHi(dir, 1)).grow(iv);
    }
    return out;
}

struct Tendencies {
    FArrayBox rho_u_rhs;
    FArrayBox rho_v_rhs;
    FArrayBox rho_w_rhs;
};

class EBMomMapFactors : public ::testing::TestWithParam<AdvType>
{
protected:
    void SetUp () override
    {
        domain = Box(IntVect(0,0,0), IntVect(kNx-1, kNy-1, kNz-1));
        amrex::RealBox rb({Real(0), Real(0), Real(0)},
                          {Real(kNx*100), Real(kNy*120), Real(kNz*50)});
        geom.define(domain, rb, 0, {1, 1, 0});

        // An all-regular EB: every cell open, so the EB operator must give exactly
        // what the non-EB operator gives
        amrex::EB2::Build(amrex::EB2::makeShop(amrex::EB2::AllRegularIF()), geom, 0, 0);
        amrex::EB2::BuildFC();

        ba = amrex::BoxArray(domain);
        dm = amrex::DistributionMapping(ba);
        ebfact.make_all_factories(0, geom, ba, dm, amrex::EB2::IndexSpace::top().getLevel(geom));

        physbnd_mask.resize(IntVars::NumTypes);
        physbnd_mask[IntVars::cons].define(ba, dm, 1, 1);
        physbnd_mask[IntVars::cons].BuildMask(domain, geom.periodicity(), 1, 1, 0, 1);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            amrex::BoxArray fba = amrex::convert(ba, IntVect::TheDimensionVector(dir));
            physbnd_mask[1+dir].define(fba, dm, 1, 1);
            physbnd_mask[1+dir].BuildMask(domain, geom.periodicity(), 1, 1, 0, 1);
        }

        const Box gbx = amrex::grow(domain, kNGrow);
        u    .resize(amrex::convert(gbx, IntVect(1,0,0)), 1); fill(u    , Real( 3.0), Real(0.1), Real(0.2));
        v    .resize(amrex::convert(gbx, IntVect(0,1,0)), 1); fill(v    , Real(-2.0), Real(0.4), Real(0.3));
        w    .resize(amrex::convert(gbx, IntVect(0,0,1)), 1); fill(w    , Real( 0.5), Real(0.7), Real(0.6));
        rho_u.resize(amrex::convert(gbx, IntVect(1,0,0)), 1); fill(rho_u, Real( 3.3), Real(0.2), Real(0.5));
        rho_v.resize(amrex::convert(gbx, IntVect(0,1,0)), 1); fill(rho_v, Real(-2.2), Real(0.9), Real(0.1));
        omega.resize(amrex::convert(gbx, IntVect(0,0,1)), 1); fill(omega, Real( 0.6), Real(0.3), Real(0.8));
        cell .resize(gbx, NVAR_max);                          cell.setVal<amrex::RunOn::Device>(Real(1));

        mf_mx.resize(mapfac_box(domain, IntVect(0,0,0)), 1); fill_mapfac(mf_mx, Real(1.10), Real(0.3));
        mf_my.resize(mapfac_box(domain, IntVect(0,0,0)), 1); fill_mapfac(mf_my, Real(0.80), Real(1.1));
        mf_ux.resize(mapfac_box(domain, IntVect(1,0,0)), 1); fill_mapfac(mf_ux, Real(1.15), Real(0.7));
        mf_uy.resize(mapfac_box(domain, IntVect(1,0,0)), 1); fill_mapfac(mf_uy, Real(0.75), Real(1.9));
        mf_vx.resize(mapfac_box(domain, IntVect(0,1,0)), 1); fill_mapfac(mf_vx, Real(1.20), Real(2.3));
        mf_vy.resize(mapfac_box(domain, IntVect(0,1,0)), 1); fill_mapfac(mf_vy, Real(0.85), Real(0.5));

        stretched_dz_d.resize(kNz, Real(50));

        bcs.resize(BCVars::NumTypes);
        for (auto& bc : bcs) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                bc.setLo(dir, ERFBCType::int_dir);
                bc.setHi(dir, ERFBCType::int_dir);
            }
        }
    }

    void TearDown () override
    {
        amrex::EB2::IndexSpace::pop();
    }

    Tendencies run (TerrainType terrain_type, const AdvType adv_type)
    {
        amrex::MultiFab cons(ba, dm, 1, 0);
        amrex::MFIter mfi(cons);

        const Box tbx = mfi.nodaltilebox(0);
        const Box tby = mfi.nodaltilebox(1);
        const Box tbz = ShrinkZmomBoxAtDomainEnds(mfi.nodaltilebox(2), domain);

        const amrex::Vector<Box> tbx_grown = grown_flux_boxes(tbx);
        const amrex::Vector<Box> tby_grown = grown_flux_boxes(tby);
        const amrex::Vector<Box> tbz_grown = grown_flux_boxes(mfi.nodaltilebox(2));

        std::array<FArrayBox, AMREX_SPACEDIM> flux_u, flux_v, flux_w;
        GpuArray<Array4<Real>, AMREX_SPACEDIM> flx_u_arr{}, flx_v_arr{}, flx_w_arr{};
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            flux_u[dir].resize(tbx_grown[dir], 1); flux_u[dir].setVal<amrex::RunOn::Device>(0);
            flux_v[dir].resize(tby_grown[dir], 1); flux_v[dir].setVal<amrex::RunOn::Device>(0);
            flux_w[dir].resize(tbz_grown[dir], 1); flux_w[dir].setVal<amrex::RunOn::Device>(0);
            flx_u_arr[dir] = flux_u[dir].array();
            flx_v_arr[dir] = flux_v[dir].array();
            flx_w_arr[dir] = flux_w[dir].array();
        }

        Tendencies out;
        out.rho_u_rhs.resize(tbx, 1); out.rho_u_rhs.setVal<amrex::RunOn::Device>(0);
        out.rho_v_rhs.resize(tby, 1); out.rho_v_rhs.setVal<amrex::RunOn::Device>(0);
        out.rho_w_rhs.resize(tbz, 1); out.rho_w_rhs.setVal<amrex::RunOn::Device>(0);

        MeshType mesh_type = MeshType::ConstantDz;
        const Real upw_frac = Real(1);
        const bool already_on_centroids = false;

        AdvectionSrcForMom(mfi, tbx, tby, tbz, tbx_grown, tby_grown, tbz_grown,
                           out.rho_u_rhs.array(), out.rho_v_rhs.array(), out.rho_w_rhs.array(),
                           cell.const_array(), u.const_array(), v.const_array(), w.const_array(),
                           rho_u.const_array(), rho_v.const_array(), omega.const_array(),
                           Array4<const Real>{}, Array4<const Real>{}, Array4<const Real>{},
                           Array4<const Real>{}, Array4<const Real>{},
                           stretched_dz_d, geom.InvCellSizeArray(),
                           mf_mx.const_array(), mf_ux.const_array(), mf_vx.const_array(),
                           mf_my.const_array(), mf_uy.const_array(), mf_vy.const_array(),
                           adv_type, adv_type, upw_frac, upw_frac,
                           mesh_type, terrain_type,
                           ebfact, flx_u_arr, flx_v_arr, flx_w_arr,
                           physbnd_mask, already_on_centroids,
                           domain.smallEnd(2), domain.bigEnd(2)+1, domain, bcs.data());
        amrex::Gpu::streamSynchronize();
        return out;
    }

    Box domain;
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    eb_ ebfact;
    amrex::Vector<amrex::iMultiFab> physbnd_mask;
    FArrayBox u, v, w, rho_u, rho_v, omega, cell;
    FArrayBox mf_mx, mf_my, mf_ux, mf_uy, mf_vx, mf_vy;
    amrex::Gpu::DeviceVector<Real> stretched_dz_d;
    amrex::Vector<amrex::BCRec> bcs;
};

// Largest |a - b| over bx, and largest |b|, computed on the device
std::pair<Real, Real> max_diff_and_scale (const FArrayBox& a, const FArrayBox& b, const Box& bx)
{
    FArrayBox diff(bx, 1);
    diff.copy<amrex::RunOn::Device>(a, bx);
    diff.minus<amrex::RunOn::Device>(b, bx, 0, 0, 1);
    return {diff.maxabs<amrex::RunOn::Device>(bx, 0), b.maxabs<amrex::RunOn::Device>(bx, 0)};
}

} // namespace

// Motivation: with an embedded boundary that cuts no cell, the EB momentum advection must
// give the same tendencies as the non-EB operator, map factors included. The EB operator
// divided the x-face mass fluxes by m_x and most y-face mass fluxes by m_y, where the
// non-EB operator (and the continuity equation) divide by m_y and m_x respectively.
TEST_P(EBMomMapFactors, AllRegularEBMatchesNonEB)
{
    const AdvType adv_type = GetParam();

    const Tendencies ref = run(TerrainType::None, adv_type);
    const Tendencies eb  = run(TerrainType::EB  , adv_type);

    // The two operators evaluate the same fluxes with the same arithmetic up to the
    // grouping of a few tens of operations per tendency, so they agree to a small
    // multiple of the working precision relative to the size of the tendency, in
    // single and double precision.
    const Real rel_tol = Real(1000) * real_eps;

    // The w tendency is compared away from the two w-faces next to each wall. There, for
    // vertical orders above 4, the non-EB operator reduces the vertical interpolation of w
    // to 4th-order centered and the EB operator to 3rd-order upwind, a difference that
    // involves no map factor (erf-model/ERF#4889).
    Box w_interior = ref.rho_w_rhs.box();
    w_interior.setSmall(2, domain.smallEnd(2) + 3);
    w_interior.setBig  (2, domain.bigEnd(2)   - 2);

    const char* names[3] = {"rho_u_rhs", "rho_v_rhs", "rho_w_rhs"};
    const FArrayBox* refs[3] = {&ref.rho_u_rhs, &ref.rho_v_rhs, &ref.rho_w_rhs};
    const FArrayBox* ebs [3] = {&eb .rho_u_rhs, &eb .rho_v_rhs, &eb .rho_w_rhs};
    const Box boxes[3] = {ref.rho_u_rhs.box(), ref.rho_v_rhs.box(), w_interior};
    for (int n = 0; n < 3; ++n) {
        SCOPED_TRACE(names[n]);
        const auto [diff, scale] = max_diff_and_scale(*ebs[n], *refs[n], boxes[n]);
        ASSERT_GT(scale, Real(0));
        ERF_EXPECT_NEAR(diff, Real(0), rel_tol * scale);
    }
}

INSTANTIATE_TEST_SUITE_P(
    AdvectionOrders, EBMomMapFactors,
    ::testing::Values(AdvType::Centered_2nd, AdvType::Upwind_3rd, AdvType::Centered_4th,
                      AdvType::Upwind_5th, AdvType::Centered_6th),
    [] (const ::testing::TestParamInfo<AdvType>& param_info) {
        switch (param_info.param) {
        case AdvType::Centered_2nd: return std::string("Centered_2nd");
        case AdvType::Upwind_3rd:   return std::string("Upwind_3rd");
        case AdvType::Centered_4th: return std::string("Centered_4th");
        case AdvType::Upwind_5th:   return std::string("Upwind_5th");
        case AdvType::Centered_6th: return std::string("Centered_6th");
        default:                    return std::string("Other");
        }
    });

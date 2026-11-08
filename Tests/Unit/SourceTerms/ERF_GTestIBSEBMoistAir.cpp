#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_Constants.H"
#include "ERF_EOS.H"
#include "ERF_IBFaceSet.H"
#include "ERF_IBSEBBalance.H"
#include "ERF_IBSEBParams.H"
#include "ERF_IndexDefines.H"
#include "ERF_TwoStreamColumn.H"

// The building faces read the air temperature and the Exner function of their fluid cell
// from the state. With a moisture model these take the cell's water vapour through the
// equation of state, as the two-stream column does (get_qv_from_state()); without one they
// take none (erf-model/ERF#4481). The density, potential temperature and vapour differ in
// every cell and the level is split into four boxes, so a face that read another cell, or
// left out the vapour, would get another number. The column next to the block's west wall
// holds a small negative vapour (an advection undershoot), which both read as zero.

namespace {

using amrex::Real;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real rho_at (int i, int, int k) { return Real(1.15 - 0.004 * k + 0.001 * i); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real th_at  (int, int j, int k) { return Real(300.0 + 0.3 * k + 0.1 * j); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real qv_at  (int i, int j, int) { return (i == 1) ? Real(-1.e-4) : Real(0.004 + 0.0012 * i + 0.0007 * j); }

void fill_blanking (amrex::MultiFab& b, const amrex::Geometry& geom)
{
    // A 30 m block over 2 x 3 columns and a 50 m tower on one column.
    const amrex::Box solid[2] = {amrex::Box(amrex::IntVect(2, 3, 0), amrex::IntVect(3, 5, 2)),
                                 amrex::Box(amrex::IntVect(6, 1, 0), amrex::IntVect(6, 1, 4))};
    b.setVal(0.0);
    for (amrex::MFIter mfi(b); mfi.isValid(); ++mfi) {
        auto const& a = b.array(mfi);
        for (const amrex::Box& s : solid) {
            const amrex::Box ov = s & mfi.validbox();
            if (ov.isEmpty()) { continue; }
            amrex::ParallelFor(ov, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { a(i, j, k) = 1.0; });
        }
    }
    b.FillBoundary(geom.periodicity());
}

// The state of every cell: rho, rho theta and, with moisture, rho qv in RhoQ1_comp.
void fill_state (amrex::MultiFab& cons, bool moist)
{
    cons.setVal(0.0);
    for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
        auto const& c = cons.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real rho = rho_at(i, j, k);
            c(i, j, k, Rho_comp) = rho;
            c(i, j, k, RhoTheta_comp) = rho * th_at(i, j, k);
            if (moist) { c(i, j, k, RhoQ1_comp) = rho * qv_at(i, j, k); }
        });
    }
}

template <class T>
std::vector<T> host (const amrex::Gpu::DeviceVector<T>& d)
{
    std::vector<T> h(d.size());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d.begin(), d.end(), h.begin());
    amrex::Gpu::streamSynchronize();
    return h;
}

// Per-cell values of a MultiFab on the host, keyed by (i, j, k).
std::map<std::tuple<int, int, int>, Real> host_cells (const amrex::MultiFab& mf, int comp)
{
    amrex::MultiFab h(mf.boxArray(), mf.DistributionMap(), 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    amrex::MultiFab::Copy(h, mf, comp, 0, 1, 0);
    amrex::Gpu::streamSynchronize();
    std::map<std::tuple<int, int, int>, Real> out;
    for (amrex::MFIter mfi(h); mfi.isValid(); ++mfi) {
        auto const& a = h.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { out[{i, j, k}] = a(i, j, k); });
    }
    return out;
}

// The two-stream column's air temperature and Exner function of a cell, from the same state.
struct Column { Real T; Real exner; Real th; };
Column column_at (int i, int j, int k, bool moist)
{
    // The two-stream helper reads qv from a state array; build a one-cell one.
    const int nc = moist ? RhoQ1_comp + 1 : RhoTheta_comp + 1;
    amrex::FArrayBox fab(amrex::Box(amrex::IntVect(i, j, k), amrex::IntVect(i, j, k)), nc, amrex::The_Pinned_Arena());
    auto const& a = fab.array();
    a(i, j, k, Rho_comp) = rho_at(i, j, k);
    a(i, j, k, RhoTheta_comp) = rho_at(i, j, k) * th_at(i, j, k);
    if (moist) { a(i, j, k, RhoQ1_comp) = rho_at(i, j, k) * qv_at(i, j, k); }
    const Real qv = get_qv_from_state(i, j, k, fab.const_array());
    return {getTgivenRandRTh(a(i, j, k, Rho_comp), a(i, j, k, RhoTheta_comp), qv),
            getExnergivenRTh(a(i, j, k, RhoTheta_comp), R_d / Cp_d, qv), th_at(i, j, k)};
}

// Run the longwave, the sensible flux and the heat source of one face set and check every
// face against the column's equation of state.
void check_faces (const std::string& lw_mode, bool moist, int& n_faces, Real& max_dT, int& n_clipped)
{
    using namespace amrex;
    // 8 x 8 columns of 20 m, 8 layers of 10 m, periodic in x and y, in four boxes.
    const RealBox rb({0.0, 0.0, 0.0}, {160.0, 160.0, 80.0});
    const Box dom(IntVect(0, 0, 0), IntVect(7, 7, 7));
    const Geometry geom(dom, rb, 0, {1, 1, 0});
    BoxArray ba(dom);
    ba.maxSize(IntVect(4, 4, 8));
    const DistributionMapping dm(ba);
    MultiFab blank(ba, dm, 1, 1);
    fill_blanking(blank, geom);
    // Without a moisture model the state has no vapour component at all.
    MultiFab cons(ba, dm, moist ? RhoQ2_comp + 1 : RhoTheta_comp + 1, 0);
    fill_state(cons, moist);
    MultiFab xvel(convert(ba, IntVect(1, 0, 0)), dm, 1, 0); xvel.setVal(2.0);
    MultiFab yvel(convert(ba, IntVect(0, 1, 0)), dm, 1, 0); yvel.setVal(1.0);
    MultiFab zvel(convert(ba, IntVect(0, 0, 1)), dm, 1, 0); zvel.setVal(0.0);

    IBSEBParams params;
    params.enable = true;
    params.lw_mode = lw_mode;
    if (lw_mode == "two_stream") { params.radiation = "two_stream"; }
    params.view_n_az = 8;
    params.view_n_el = 4;
    params.T_skin_init = 310.0;
    IBFaceSet faces(params, 0, moist ? RhoQ1_comp : -1);
    faces.build(blank, geom);
    faces.release_labels();
    faces.compute_view_fractions();
    faces.compute_longwave(cons);
    const auto Ta_lw = host(faces.d_T_air), lwext = host(faces.d_LW_ext);
    faces.compute_sensible(cons, xvel, yvel, zvel, Cp_d);
    const auto Ta = host(faces.d_T_air), th = host(faces.d_theta_air), H = host(faces.d_H), Hc = host(faces.d_H_coeff);
    const auto i = host(faces.d_i), j = host(faces.d_j), k = host(faces.d_k);
    const auto fs = host(faces.d_f_sky), fg = host(faces.d_f_ground), Ts = host(faces.d_T_skin), A = host(faces.d_area);

    MultiFab src(ba, dm, RhoTheta_comp + 1, 0);
    src.setVal(0.0);
    faces.add_heat_flux_to_source(src, cons, geom, Cp_d, R_d / Cp_d);
    auto tend = host_cells(src, RhoTheta_comp);

    const Real tol = (sizeof(Real) == 8) ? Real(1.e-12) : Real(1.e-5);
    const Real sigma = ibseb::SIGMA_SB;
    const Real vol = Real(20.0 * 20.0 * 10.0);
    n_faces = faces.n_faces();
    std::map<std::tuple<int, int, int>, Real> want_tend;
    for (int f = 0; f < n_faces; ++f) {
        const Column col = column_at(i[f], j[f], k[f], moist);
        EXPECT_NEAR(Ta_lw[f], col.T, tol * col.T) << lw_mode << " longwave, face " << f;
        EXPECT_NEAR(Ta[f], col.T, tol * col.T) << lw_mode << " sensible, face " << f;
        EXPECT_NEAR(th[f], col.th, tol * col.th) << "face " << f;
        // The gray sky emits at the air temperature; the ground term is fixed.
        if (lw_mode == "gray") {
            const Real T4 = col.T * col.T * col.T * col.T;
            const Real Tg4 = params.T_ground * params.T_ground * params.T_ground * params.T_ground;
            const Real lw = fs[f] * params.sky_emissivity * sigma * T4 + fg[f] * params.emissivity_ground * sigma * Tg4;
            EXPECT_NEAR(lwext[f], lw, tol * (lw + 1.0)) << "face " << f;
        }
        // H = H_coeff (theta_skin - theta_air), the skin temperature over the cell's Exner function.
        const Real dth = Ts[f] / col.exner - col.th;
        EXPECT_NEAR(H[f], Hc[f] * dth, tol * (std::abs(H[f]) + 1.0)) << "face " << f;
        want_tend[{i[f], j[f], k[f]}] += H[f] * A[f] / (Cp_d * vol * col.exner);
        if (moist) { max_dT = std::max(max_dT, std::abs(col.T - column_at(i[f], j[f], k[f], false).T)); }
        if (i[f] == 1) { ++n_clipped; }
    }
    // The rho theta source of every face cell, and nothing elsewhere.
    for (const auto& [cell, v] : tend) {
        const auto it = want_tend.find(cell);
        const Real w = (it == want_tend.end()) ? Real(0.0) : it->second;
        EXPECT_NEAR(v, w, tol * (std::abs(w) + 1.e-6)) << "cell " << std::get<0>(cell) << ","
                                                        << std::get<1>(cell) << "," << std::get<2>(cell);
    }
}

} // namespace

TEST(IBSEBMoistAir, GrayLongwaveTakesTheVapour)
{
    int n = 0, n_clip = 0; Real dT = 0.0;
    check_faces("gray", true, n, dT, n_clip);
    // Not vacuous: many faces, a vapour that moves the air temperature by kelvins, and
    // faces on the column whose vapour is negative.
    EXPECT_GT(n, 30);
    EXPECT_GT(dT, Real(1.0));
    EXPECT_GT(n_clip, 2);
}

TEST(IBSEBMoistAir, TwoStreamLongwaveTakesTheVapour)
{
    int n = 0, n_clip = 0; Real dT = 0.0;
    check_faces("two_stream", true, n, dT, n_clip);
    EXPECT_GT(n, 30);
    EXPECT_GT(dT, Real(1.0));
    EXPECT_GT(n_clip, 2);
}

// No moisture model: the state has no vapour component, and the faces must not read one.
TEST(IBSEBMoistAir, DryStateHasNoVapour)
{
    int n = 0, n_clip = 0; Real dT = 0.0;
    check_faces("gray", false, n, dT, n_clip);
    EXPECT_GT(n, 30);
}

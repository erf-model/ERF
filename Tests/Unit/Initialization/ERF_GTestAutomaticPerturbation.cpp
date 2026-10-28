#include "ERF_DataStruct.H"
#include "ERF_PBLHeight.H"
#include "ERF_PerturbationWind.H"

#include <AMReX_Gpu.H>
#include <AMReX_ParmParse.H>
#include <gtest/gtest.h>
#include <algorithm>

using namespace amrex;

namespace {

void fill_test_heights (MultiFab& cc, MultiFab& nd)
{
    for (MFIter mfi(cc); mfi.isValid(); ++mfi) {
        const auto a = cc.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            a(i,j,k) = Real(1000) + k + myhalf;
        });
    }
    for (MFIter mfi(nd); mfi.isValid(); ++mfi) {
        const auto a = nd.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            a(i,j,k) = Real(1000) + k;
        });
    }
}

// Check both the stored mask and its application: a shared corner must receive
// exactly one increment. Copy to pinned memory for the same checks on CPU/GPU.
void check_region (TurbulentPerturbation& pert, bool w_faces, bool reversed, Real height)
{
    auto& field = pert.pb_cell[1];
    MultiFab applied(field.boxArray(), field.DistributionMap(), 1, 0);
    applied.setVal(zero);
    for (MFIter mfi(applied); mfi.isValid(); ++mfi) {
        pert.apply_tpi(1, mfi.validbox(), w_faces ? -1 : 0, applied.ixType(),
                       applied.array(mfi), field.const_array(mfi));
    }
    MultiFab host(field.boxArray(), field.DistributionMap(), 2, 0,
                  MFInfo().SetArena(The_Pinned_Arena()));
    MultiFab::Copy(host, field, 0, 0, 1, 0);
    MultiFab::Copy(host, applied, 0, 1, 1, 0);
    Gpu::streamSynchronize();
    int nonzero = 0;
    for (MFIter mfi(host); mfi.isValid(); ++mfi) {
        const auto a = host.const_array(mfi);
        const auto lo = lbound(mfi.validbox());
        const auto hi = ubound(mfi.validbox());
        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    const Real z = Real(k) + (w_faces ? zero : myhalf);
                    const bool edge = reversed ? (i >= 9 || j >= 9) : (i < 3 || j < 3);
                    const bool active = edge && z > zero && z <= height;
                    if (!active) { EXPECT_EQ(a(i,j,k,0), zero) << i << ',' << j << ',' << k; }
                    EXPECT_EQ(a(i,j,k,0), a(i,j,k,1)) << "Perturbation applied more than once";
                    if (a(i,j,k,0) != zero) { ++nonzero; }
                }
            }
        }
    }
    ParallelDescriptor::ReduceIntSum(nonzero);
    EXPECT_GT(nonzero, 0);
}

void fill_pbl_column (MultiFab& state, MultiFab& zcc, Real elevation)
{
    state.setVal(zero);
    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        const auto s = state.array(mfi);
        const auto z = zcc.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            s(i,j,k,Rho_comp) = one;
            s(i,j,k,RhoTheta_comp) = (k < 4) ? Real(300) : Real(302);
            z(i,j,k) = elevation + k + myhalf;
        });
    }
}

void fill_sampling_profile (MultiFab& u, MultiFab& v, MultiFab& zcc, MultiFab& znd,
                            int klo, bool terrain)
{
    for (MFIter mfi(zcc); mfi.isValid(); ++mfi) {
        const auto a = zcc.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real n = k-klo;
            const Real z = terrain ? myhalf*(n*n+(n+1)*(n+1)) : two*(n+myhalf);
            a(i,j,k) = Real(1000) + myhalf*(i+myhalf) + fourth*(j+myhalf) + z;
        });
    }
    for (MFIter mfi(znd); mfi.isValid(); ++mfi) {
        const auto a = znd.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real n = k-klo;
            a(i,j,k) = Real(1000) + myhalf*i + fourth*j + (terrain ? n*n : two*n);
        });
    }
    for (int component = 0; component < 2; ++component) {
        auto& velocity = component == 0 ? u : v;
        for (MFIter mfi(velocity); mfi.isValid(); ++mfi) {
            const auto a = velocity.array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const Real n = k-klo;
                const Real z = terrain ? myhalf*(n*n+(n+1)*(n+1)) : two*(n+myhalf);
                a(i,j,k) = component == 0 ? two + Real(0.1)*z + Real(0.01)*(i-myhalf)
                                         : -one + Real(0.05)*z + Real(0.02)*(j-myhalf);
            });
        }
    }
}

TEST(AutomaticPerturbationWind, SamplesAGLAcrossStackedBoxesAndRejectsShallowCoverage)
{
    const Box domain(IntVect(-4,2,3), IntVect(7,13,14));
    const RealBox physical({0.,0.,1000.}, {12.,12.,1024.});
    const int periodic[3] = {0,0,0};
    const Geometry geom(domain, &physical, 0, periodic);
    Box edge = domain;
    edge.setRange(0, domain.smallEnd(0));
    for (const bool terrain : {false, true}) {
        for (const bool split : {false, true}) {
            for (const bool shallow : {false, true}) {
                Box covered = domain;
                if (shallow) { covered.setBig(2, domain.smallEnd(2)+3); }
                BoxArray ba(covered);
                if (split) { ba.maxSize(6); }
                const DistributionMapping dm(ba);
                MultiFab state(ba, dm, 1, 0);
                MultiFab u(convert(ba, IntVect(1,0,0)), dm, 1, 0);
                MultiFab v(convert(ba, IntVect(0,1,0)), dm, 1, 0);
                MultiFab zcc(ba, dm, 1, 0);
                MultiFab znd(convert(ba, IntVect::TheNodeVector()), dm, 1, 0);
                fill_sampling_profile(u, v, zcc, znd, domain.smallEnd(2), terrain);
                // Both targets bracket centres on opposite sides of the z-box split.
                const Real target = terrain ? Real(36) : Real(12);
                const auto wind = sample_perturbation_wind(edge, target, geom, state, u, v,
                                                          terrain ? &zcc : nullptr, terrain ? &znd : nullptr);
                EXPECT_EQ(wind.complete, !shallow);
                if (wind.complete) {
                    const Real tolerance = Real(256)*real_eps*target;
                    EXPECT_NEAR(wind.u, two+Real(0.1)*target-Real(0.04), tolerance);
                    EXPECT_NEAR(wind.v, -one+Real(0.05)*target+Real(0.15), tolerance);
                }
                const Real above_top = terrain ? Real(145) : Real(25);
                EXPECT_FALSE(sample_perturbation_wind(edge, above_top, geom, state, u, v,
                    terrain ? &zcc : nullptr, terrain ? &znd : nullptr).complete);
            }
        }
    }
}

TEST(AutomaticPerturbationWind, CPMAmplitudeUsesAbovePBLWindIndependentlyOfAdvection)
{
    const Box domain(IntVect(0), IntVect(11));
    const RealBox physical({0.,0.,0.}, {12.,12.,12.});
    const int periodic[3] = {0,0,0};
    const Geometry geom(domain, &physical, 0, periodic);
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    MultiFab state(ba, dm, 1, 0);
    MultiFab u(convert(ba, IntVect(1,0,0)), dm, 1, 0);
    MultiFab v(convert(ba, IntVect(0,1,0)), dm, 1, 0);
    state.setVal(zero); u.setVal(zero); v.setVal(zero);
    const std::string prefix = "test_auto_cpm_Ug";
    ParmParse pp(prefix);
    pp.addarr("perturbation_box_dims", Vector<int>{3,3,3});
    pp.add("perturbation_layers", 1);
    pp.add("perturbation_offset", 0);
    Vector<MultiFab> fields(3);
    const Real speeds[3] = {Real(5), Real(10), zero};
    for (int n = 0; n < 3; ++n) {
        TurbulentPerturbation pert;
        pert.init_tpi_type(1, PerturbationType::CPM, 1);
        pert.init_tpi(1, Vector<BoxArray>{ba}, geom.CellSizeArray(), ba, dm, 0, prefix,
                      Vector<IntVect>{IntVect(2)}, 1, true);
        Vector<AutomaticPerturbationFace> controls(1);
        controls[0].region = pert.automatic_face_box(domain, 0);
        controls[0].u = two;
        controls[0].height = Real(8);
        controls[0].Ug = speeds[n];
        pert.configure_automatic(1, 0., 600., controls, geom, nullptr, nullptr, nullptr);
        pert.calc_tpi_update(1, 1., u, v, state);
        for (const Real speed : pert.pb_mag[1]) {
            EXPECT_TRUE(speed == zero || speed == two);
        }
        fields[n].define(ba, dm, 1, 0);
        MultiFab::Copy(fields[n], pert.pb_cell[1], 0, 0, 1, 0);
    }
    EXPECT_GT(fields[0].norm0(), zero);
    EXPECT_EQ(fields[2].norm0(), zero);
    fields[0].mult(four);
    const Real tolerance = Real(64)*real_eps*std::max(one, fields[0].norm0());
    MultiFab::Subtract(fields[1], fields[0], 0, 0, 1, 0);
    EXPECT_LE(fields[1].norm0(), tolerance);
}

TEST(AutomaticPerturbationHeight, MYNNIsIndependentOfTerrainElevation)
{
    const Box domain(IntVect(0), IntVect(11));
    const RealBox physical({0.,0.,0.}, {12.,12.,12.});
    const int periodic[3] = {0,0,0};
    const Geometry geom(domain, &physical, 0, periodic);
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    MultiFab state(ba, dm, RhoKE_comp+1, 1);
    MultiFab zcc(ba, dm, 1, 1);
    Box slab = domain;
    slab.setRange(2, 0);
    MultiFab height(BoxArray(slab), dm, 1, 0);
    MYNNPBLH estimator;
    estimator.use_agl = true;
    for (const Real elevation : {zero, Real(1000)}) {
        fill_pbl_column(state, zcc, elevation);
        estimator.compute_pblh(geom, &zcc, &height, state, nullptr, MoistureComponentIndices{});
        // Land threshold is 1.25 K: interpolate from z=3.5 to 4.5 m AGL.
        const Real tolerance = Real(64)*real_eps;
        EXPECT_NEAR(height.min(0), Real(4.125), tolerance);
        EXPECT_NEAR(height.max(0), Real(4.125), tolerance);
    }
}

class AutomaticPerturbation : public ::testing::TestWithParam<PerturbationType> {};

TEST_P(AutomaticPerturbation, UpdatesEdgesHeightAndCalmWinds)
{
    const auto type = GetParam();
    const bool w_faces = (type == PerturbationType::CPM_W);
    const Box domain(IntVect(0), IntVect(11));
    const RealBox physical({0.,0.,0.}, {12.,12.,12.});
    const int periodic[3] = {0,0,0};
    const Geometry geom(domain, &physical, 0, periodic);
    const Vector<BoxArray> subdomains{BoxArray(domain)};
    const std::string prefix = "test_auto_perturbation_" + std::to_string(static_cast<int>(type));
    ParmParse pp(prefix);
    pp.addarr("perturbation_box_dims", Vector<int>{3,3,3});
    pp.add("perturbation_layers", 1);
    pp.add("perturbation_offset", 0);
    pp.add("perturbation_nondimensional", Real(0.042));
    pp.add("perturbation_w_amp", Real(0.2));
    // No direction, Ug, or k bounds: all are resolved by automatic controls.

    for (const bool split : {false, true}) {
        BoxArray ba(domain);
        if (split) { ba.maxSize(6); }
        DistributionMapping dm(ba);
        MultiFab state(ba, dm, 1, 0);
        MultiFab u(convert(ba, IntVect(1,0,0)), dm, 1, 0);
        MultiFab v(convert(ba, IntVect(0,1,0)), dm, 1, 0);
        state.setVal(zero); u.setVal(zero); v.setVal(zero);
        MultiFab zcc(ba, dm, 1, 0);
        MultiFab znd(convert(ba, IntVect::TheNodeVector()), dm, 1, 0);
        BoxArray ground_ba = znd.boxArray();
        // The ground slab must not change the shared 3D state/height layout.
        ground_ba.uniqify();
        for (int b = 0; b < ground_ba.size(); ++b) {
            Box bx = ground_ba[b]; bx.setRange(2, 0); ground_ba.set(b, bx);
        }
        MultiFab ground(ground_ba, dm, 1, 0);
        ground.setVal(Real(1000));
        fill_test_heights(zcc, znd);

        TurbulentPerturbation pert;
        pert.init_tpi_type(1, type, 1);
        pert.init_tpi(1, subdomains, geom.CellSizeArray(), ba, dm, 0, prefix,
                      Vector<IntVect>{IntVect(2)}, 1, true);
        EXPECT_TRUE(pert.automatic_update_due(1, 0.));
        EXPECT_TRUE(pert.pb_ba[1].isDisjoint());

        Vector<AutomaticPerturbationFace> controls(4);
        for (int f = 0; f < 4; ++f) {
            controls[f].region = pert.automatic_face_box(domain, f);
            controls[f].face = f;
            controls[f].u = two;
            controls[f].v = Real(1);
            controls[f].Ug = Real(5);
            controls[f].height = Real(8); // 80% of a 10 m PBL
        }
        pert.configure_automatic(1, 0., 600., controls, geom, &zcc, &znd, &ground);
        EXPECT_FALSE(pert.automatic_update_due(1, 599.));
        EXPECT_TRUE(pert.automatic_update_due(1, 600.));
        pert.calc_tpi_update(1, 1., u, v, state);
        check_region(pert, w_faces, false, Real(8));

        for (auto& c : controls) { c.u = zero; c.v = zero; }
        pert.configure_automatic(1, 600., 600., controls, geom, &zcc, &znd, &ground);
        pert.calc_tpi_update(1, 1., u, v, state);
        EXPECT_EQ(pert.pb_cell[1].norm0(), zero);

        for (auto& c : controls) { c.u = Real(-2); c.v = Real(-1); c.height = two; }
        pert.configure_automatic(1, 1200., 600., controls, geom, &zcc, &znd, &ground);
        pert.calc_tpi_update(1, 1., u, v, state);
        check_region(pert, w_faces, true, two);
        if (type == PerturbationType::Source || type == PerturbationType::Direct) {
            const Real tolerance = Real(256)*real_eps
                                 * std::max(one, pert.pb_cell[1].norm1());
            EXPECT_NEAR(pert.pb_cell[1].sum(0), zero, tolerance);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(AllMethods, AutomaticPerturbation,
                        ::testing::Values(PerturbationType::Source, PerturbationType::Direct,
                                          PerturbationType::CPM, PerturbationType::CPM_W));

} // namespace

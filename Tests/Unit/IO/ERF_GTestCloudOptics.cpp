#include <algorithm>
#include <cmath>
#include <limits>

#include <AMReX_Gpu.H>
#include <AMReX_ParmParse.H>
#include <gtest/gtest.h>

#include "ERF_IndexDefines.H"
#include "ERF_PlanarBoundary.H"
#include "ERF_Plotfile2DCloudOptics.H"
#include "ERF_Plotfile2DMetadata.H"
#include "ERF_Plotfile2DUtils.H"
#include "ERF_Plotfile2DWaterPath.H"

using namespace amrex;
using namespace plotfile2d;

namespace {
struct MeshGuard
{
    MeshType old = SolverChoice::mesh_type;
    ~MeshGuard() { SolverChoice::set_mesh_type(old); }
};

Real
tolerance (Real value)
{
    // Several products and parallel sums can change the last few bits.
    return 128. * std::numeric_limits<Real>::epsilon() *
           std::max(Real(1.), std::abs(value));
}

// Kept out of the TEST bodies: nvcc rejects extended device lambdas in them.
void
set_metric_k_plus_one (MultiFab& metric)
{
    for (MFIter mfi(metric); mfi.isValid(); ++mfi) {
        const auto jac = metric.array(mfi);
        ParallelFor(mfi.validbox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        jac(i, j, k) = Real(k + 1);
                    });
    }
}

struct Column
{
    Box domain{IntVect(0, 0, 0), IntVect(1, 1, 2)};
    BoxArray ba;
    BoxArray slab;
    DistributionMapping dm;
    Geometry geom;
    MultiFab cons, metric, dst;

    explicit Column (bool split = false, int covered_layers = 3)
    {
        Box covered(domain);
        covered.setRange(2, 0, covered_layers);
        ba = BoxArray(covered);
        if (split) {
            ba.maxSize(IntVect(2, 2, 1));
        }
        BoxList bl = ba.boxList();
        for (auto& box : bl) {
            box.setRange(2, 0);
        }
        slab = BoxArray(std::move(bl));
        dm = DistributionMapping(ba);
        RealBox rb({0., 0., 0.}, {2., 2., 3.});
        geom.define(domain, rb, 0, Array<int, 3>{{0, 0, 0}});
        cons.define(ba, dm, RhoQ6_comp + 1, 0);
        metric.define(ba, dm, 1, 0);
        dst.define(slab, dm, 4, 0);
        cons.setVal(0.);
        metric.setVal(1.);
        dst.setVal(42.);
    }

    void layer (int comp, int height, Real value)
    {
        for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
            const auto src = cons.array(mfi);
            ParallelFor(mfi.validbox(),
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            if (k == height) {
                                src(i, j, k, comp) = value;
                            }
                        });
        }
    }

    void fill (const MoistureComponentIndices& indices,
               const CloudOpticsParameters& p,
               int tau = 1,
               int albedo = 2)
    {
        fill_cloud_optics(dst, tau, albedo, cons, indices, p, geom, metric);
        Gpu::streamSynchronize();
    }

    void expect (int comp, Real value)
    {
        EXPECT_NEAR(dst.min(comp), value, tolerance(value));
        EXPECT_NEAR(dst.max(comp), value, tolerance(value));
    }
};

MoistureComponentIndices
kessler ()
{
    return MoistureComponentIndices::from_moisture_model(MoistureType::Kessler);
}
} // namespace

// The same source is registered in the MPI target so stacked sums also cross
// ranks.
TEST(CloudOptics, CatalogSelectionAndMetadata)
{
    SolverChoice sc;
    sc.moisture_indices = kessler();
    const Vector<std::string> requested{
        "cloud_albedo_est", "surface_diagnostic_source",
        "cloud_optical_depth_sw", "integrated_qc", "cloud_albedo_est"};
    const auto selection = select_requested_plot_variables(
        requested, available_diagnostic_names(sc));
    const Vector<std::string> expected{
        "integrated_qc", "cloud_optical_depth_sw", "cloud_albedo_est",
        "surface_diagnostic_source"};
    EXPECT_EQ(selection.accepted, expected);
    EXPECT_TRUE(selection.unavailable.empty());
    const auto names = diagnostic_names();
    const auto qg = std::find(names.begin(), names.end(), "integrated_qg");
    ASSERT_NE(qg, names.end());
    ASSERT_GE(std::distance(qg, names.end()), 3);
    EXPECT_EQ(*(qg + 1), "cloud_optical_depth_sw");
    EXPECT_EQ(*(qg + 2), "cloud_albedo_est");
    for (const char* name : {"cloud_optical_depth_sw", "cloud_albedo_est"}) {
        const auto* d = find_diagnostic(name);
        ASSERT_NE(d, nullptr);
        EXPECT_STREQ(d->units, "1");
        EXPECT_EQ(d->category, DiagnosticCategory::Radiation);
        EXPECT_EQ(d->missing_policy,
                  MissingPolicy::FillMinus999WhenUnavailable);
        const std::string json =
            format_2d_metadata_json(Vector<std::string>{name});
        EXPECT_NE(json.find(std::string("\"name\": \"") + name + "\""),
                  std::string::npos);
        EXPECT_NE(json.find("\"units\": \"1\""), std::string::npos);
        EXPECT_NE(json.find("\"category\": \"Radiation\""), std::string::npos);
        EXPECT_NE(json.find("\"missing_value\": -999"), std::string::npos);
    }
}

TEST(CloudOptics, AvailabilityFollowsConservedMassOnly)
{
    for (auto type :
         {MoistureType::None, MoistureType::Kessler, MoistureType::Morrison}) {
        SolverChoice sc;
        sc.moisture_type = type;
        sc.moisture_indices =
            MoistureComponentIndices::from_moisture_model(type);
        const auto names = available_diagnostic_names(sc, false);
        for (const char* name :
             {"cloud_optical_depth_sw", "cloud_albedo_est"}) {
            EXPECT_EQ(std::find(names.begin(), names.end(), name) !=
                          names.end(),
                      type != MoistureType::None);
        }
    }
    auto mi = kessler();
    EXPECT_GE(mi.qc, 0);
    EXPECT_GE(mi.qr, 0);
    EXPECT_LT(mi.qi, 0);
    EXPECT_LT(mi.qs, 0);
    EXPECT_LT(mi.qg, 0);
    SolverChoice sc;
    sc.moisture_indices.qg = RhoQ2_comp;
    EXPECT_TRUE(cloud_optics_available(sc.moisture_indices));
}

TEST(CloudOptics, AnalyticLiquidColumnAndSingleRequests)
{
    MeshGuard guard;
    SolverChoice::set_mesh_type(MeshType::ConstantDz);
    const auto mi = kessler();
    for (int selection = 0; selection < 3; ++selection) {
        Column c;
        // Nonunit computational dz detects missing or repeated dz scaling.
        const RealBox rb({0., 0., 0.}, {2., 2., 6.});
        c.geom.define(c.domain, rb, 0, Array<int, 3>{{0, 0, 0}});
        c.layer(mi.qc, 0, .01);
        c.layer(mi.qc, 1, .015);
        c.layer(mi.qc, 2, .025);
        // Density must never be applied a second time.
        c.cons.setVal(10., Rho_comp, 1);
        c.fill(mi, {}, selection == 1 ? -1 : 1, selection == 2 ? -1 : 2);
        c.expect(1, selection == 1 ? 42. : 15.);
        c.expect(2, selection == 2 ? 42. : 15. / 22.7);
        c.expect(0, 42.);
        c.expect(3, 42.);
    }
}

TEST(CloudOptics, FiveSpeciesAndSwappedMappings)
{
    MeshGuard guard;
    SolverChoice::set_mesh_type(MeshType::ConstantDz);
    auto mi =
        MoistureComponentIndices::from_moisture_model(MoistureType::Morrison);
    for (bool swap : {false, true}) {
        if (swap) {
            std::swap(mi.qc, mi.qg);
            std::swap(mi.qi, mi.qr);
        }
        Column c;
        const GpuArray<int, 5> components{{mi.qc, mi.qi, mi.qr, mi.qs, mi.qg}};
        CloudOpticsParameters p;
        Real expected = 0.;
        for (int n = 0; n < 5; ++n) {
            ASSERT_GE(components[n], 0);
            c.cons.setVal(.001 * (n + 1), components[n], 1);
            expected += 3. * .001 * (n + 1) * p.kext[n];
        }
        c.fill(mi, p);
        c.expect(1, expected);
        c.expect(2, expected / (expected + 7.7));
    }
}

TEST(CloudOptics, TerrainMetricMatchesWaterPath)
{
    MeshGuard guard;
    const auto mi = kessler();
    Column c;
    c.cons.setVal(.1, mi.qc, 1);
    for (auto mesh : {MeshType::ConstantDz, MeshType::VariableDz}) {
        SolverChoice::set_mesh_type(mesh);
        c.fill(mi, {});
        c.expect(1, 45.);
    }
    set_metric_k_plus_one(c.metric);
    c.fill(mi, {});
    c.expect(1, 90.); // .1 * (1+2+3) * 150 * dz=1
    SelectedWaterPathComponents selected;
    selected.n = 1;
    selected.src_comp[0] = mi.qc;
    selected.dst_comp[0] = 0;
    fill_condensed_water_paths(c.dst, c.cons, selected, c.geom, c.metric);
    Gpu::streamSynchronize();
    c.expect(0, .6);
}

TEST(CloudOptics, ZeroNegativeAbsentAndUnsupportedSpecies)
{
    MeshGuard guard;
    SolverChoice::set_mesh_type(MeshType::ConstantDz);
    Column c;
    auto mi = kessler();
    c.fill(mi, {});
    c.expect(1, 0.);
    c.expect(2, 0.);
    c.layer(mi.qc, 0, -1.e-12);
    c.fill(mi, {});
    c.expect(1, 0.);
    c.expect(2, 0.);
    // All other storage, including number concentrations, must be ignored.
    c.cons.setVal(999.);
    c.cons.setVal(0., mi.qc, 1);
    c.cons.setVal(0., mi.qr, 1);
    c.fill(mi, {});
    c.expect(1, 0.);
    c.fill(MoistureComponentIndices{}, {});
    c.expect(1, -999.);
    c.expect(2, -999.);
}

TEST(CloudOptics, OverridesAndZeroCoefficient)
{
    MeshGuard guard;
    SolverChoice::set_mesh_type(MeshType::ConstantDz);
    Column c;
    const auto mi = kessler();
    c.cons.setVal(.1, mi.qc, 1);
    c.cons.setVal(1., mi.qr, 1);
    CloudOpticsParameters p;
    p.kext[0] = 20.;
    p.kext[2] = 0.;
    p.reflectance_scale = 3.;
    c.fill(mi, p);
    c.expect(1, 6.);
    c.expect(2, 2. / 3.);
}

TEST(CloudOptics, ParameterValidationAndParsing)
{
    CloudOpticsParameters p;
    EXPECT_TRUE(cloud_optics_parameter_error(p).empty());
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    const Real inf = std::numeric_limits<Real>::infinity();
    for (int n = 0; n < 5; ++n) {
        for (Real bad : {Real(-1.), nan, inf}) {
            auto invalid = p;
            invalid.kext[n] = bad;
            EXPECT_FALSE(cloud_optics_parameter_error(invalid).empty());
        }
    }
    for (Real bad : {Real(0.), Real(-1.), nan, inf}) {
        auto invalid = p;
        invalid.reflectance_scale = bad;
        EXPECT_NE(
            cloud_optics_parameter_error(invalid).find("reflectance_scale"),
            std::string::npos);
    }
    ParmParse pp("erf.plot2d.cloud_optics");
    const char* keys[] = {"kext_qc", "kext_qi", "kext_qr", "kext_qs",
                          "kext_qg"};
    for (int n = 0; n < 5; ++n) {
        pp.add(keys[n], Real(n));
    }
    pp.add("reflectance_scale", Real(4.));
    const auto parsed = parse_cloud_optics_parameters();
    for (int n = 0; n < 5; ++n) {
        EXPECT_EQ(parsed.kext[n], Real(n));
        pp.remove(keys[n]);
    }
    EXPECT_EQ(parsed.reflectance_scale, 4.);
    pp.remove("reflectance_scale");
}

TEST(CloudOptics, ReflectanceBoundsAndLargeFiniteTau)
{
    for (Real tau : {Real(0.), Real(1.e-20), Real(15.), Real(1.e20),
                     std::numeric_limits<Real>::max()}) {
        const Real a = cloud_reflectance(tau, 7.7);
        EXPECT_TRUE(std::isfinite(a));
        EXPECT_GE(a, 0.);
        EXPECT_LE(a, 1.);
    }
    const Real largest = std::numeric_limits<Real>::max();
    EXPECT_EQ(cloud_reflectance(largest, largest), .5);
    EXPECT_EQ(cloud_reflectance(0., 7.7), 0.);
    EXPECT_EQ(cloud_reflectance(-1., 7.7), -999.);
    EXPECT_EQ(cloud_reflectance(largest, 0.), -999.);
}

TEST(CloudOptics, InvalidStateAndOverflowAreMissing)
{
    MeshGuard guard;
    SolverChoice::set_mesh_type(MeshType::ConstantDz);
    const auto mi = kessler();
    for (Real bad : {std::numeric_limits<Real>::quiet_NaN(),
                     std::numeric_limits<Real>::infinity(),
                     -std::numeric_limits<Real>::infinity(),
                     std::numeric_limits<Real>::max()}) {
        Column c;
        c.layer(mi.qc, 2, bad);
        c.fill(mi, {});
        c.expect(1, -999.);
        c.expect(2, -999.);
    }
    Column c;
    c.layer(mi.qc, 2, std::numeric_limits<Real>::quiet_NaN());
    CloudOpticsParameters p;
    p.kext[0] = 0.;
    c.fill(mi, p);
    c.expect(1, -999.);
    SolverChoice::set_mesh_type(MeshType::VariableDz);
    c.cons.setVal(0.);
    for (Real bad_metric :
         {std::numeric_limits<Real>::quiet_NaN(),
          std::numeric_limits<Real>::infinity(), Real(0.), Real(-1.)}) {
        c.metric.setVal(bad_metric);
        c.fill(mi, p);
        c.expect(1, -999.);
        c.expect(2, -999.);
    }
}

TEST(CloudOptics, ZSplitSumPrecedesNonlinearTransform)
{
    MeshGuard guard;
    SolverChoice::set_mesh_type(MeshType::ConstantDz);
    const auto mi = kessler();
    Column c(true);
    c.layer(mi.qc, 0, .02);
    c.layer(mi.qc, 1, .03);
    c.layer(mi.qc, 2, .05);
    c.fill(mi, {});
    c.expect(1, 15.);
    c.expect(2, 15. / 22.7);
    c.layer(mi.qc, 2, std::numeric_limits<Real>::quiet_NaN());
    c.fill(mi, {});
    c.expect(1, -999.);
    c.expect(2, -999.);
}

// Nonunit dz and a varying metric must be included before the cross-box sum;
// single-field requests must use the same completed weighted optical depth.
TEST(CloudOptics, ZSplitTerrainMetricAndSingleRequests)
{
    MeshGuard guard;
    SolverChoice::set_mesh_type(MeshType::VariableDz);
    const auto mi = kessler();
    for (int selection = 0; selection < 3; ++selection) {
        Column c(true);
        const RealBox rb({0., 0., 0.}, {2., 2., 6.});
        c.geom.define(c.domain, rb, 0, Array<int, 3>{{0, 0, 0}});
        c.layer(mi.qc, 0, .01);
        c.layer(mi.qc, 1, .03);
        c.layer(mi.qc, 2, .05);
        set_metric_k_plus_one(c.metric);
        // W = dz * (.01*1 + .03*2 + .05*3) = .44 kg/m^2.
        const Real expected_tau = 66.;
        c.fill(mi, {}, selection == 1 ? -1 : 1, selection == 2 ? -1 : 2);
        c.expect(1, selection == 1 ? 42. : expected_tau);
        c.expect(2, selection == 2 ? 42. : expected_tau / (expected_tau + 7.7));
        c.expect(3, 42.);

        SelectedWaterPathComponents selected;
        selected.n = 1;
        selected.src_comp[0] = mi.qc;
        selected.dst_comp[0] = 0;
        fill_condensed_water_paths(c.dst, c.cons, selected, c.geom, c.metric);
        sum_over_stacked_boxes(c.dst, 0, 1);
        Gpu::streamSynchronize();
        c.expect(0, .44);
    }
}

// An explicitly partial level has cells only in the bottom two of the three
// domain layers. Unlike a complete column split into boxes, its missing top
// layer is deliberately not supplied from another level or masked as missing.
TEST(CloudOptics, PartialLevelColumnsRetainOnlyCoveredMass)
{
    MeshGuard guard;
    SolverChoice::set_mesh_type(MeshType::ConstantDz);
    const auto mi = kessler();
    Column complete(true);
    complete.cons.setVal(.1, mi.qc, 1);
    complete.fill(mi, {});
    complete.expect(1, 45.);
    complete.expect(2, 45. / 52.7);

    for (int selection = 0; selection < 3; ++selection) {
        Column partial(true, 2);
        EXPECT_EQ(partial.geom.Domain().bigEnd(2), 2);
        EXPECT_EQ(partial.ba.minimalBox().bigEnd(2), 1);
        EXPECT_LT(partial.ba.numPts(), partial.geom.Domain().numPts());
        partial.cons.setVal(.1, mi.qc, 1);
        partial.fill(mi, {}, selection == 1 ? -1 : 1, selection == 2 ? -1 : 2);
        // This calls the same planar-copy operation used by the writer.
        // Its ground source exists even though the level does not reach the
        // top.
        PlanarBoundary planar;
        planar.define(partial.ba, partial.slab, partial.dm, 0);
        planar.fill(partial.dst, partial.geom.periodicity());
        Gpu::streamSynchronize();
        partial.expect(1, selection == 1 ? 42. : 30.);
        partial.expect(2, selection == 2 ? 42. : 30. / 37.7);
        partial.expect(0, 42.);
        partial.expect(3, 42.);
    }
}

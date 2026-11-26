// The liquid and frozen parts of the surface precipitation when a scheme's total
// holds the cloud ice that fell out (SurfacePrecipAccumulationSources::rest_is_frozen).
//
// WSM6 and WDM6 add rain, snow, graupel and cloud-ice fallout to their total, and
// keep rain, snow and graupel alone in their own accumulators. The cloud ice is then
// the total less those three, and it is frozen. Every value below is a small integer,
// so each expected value is exact in single and double precision.

#include <limits>
#include <memory>
#include <string>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_RealBox.H>

#include <gtest/gtest.h>

#include <ERF_Constants.H>
#include <ERF_DataStruct.H>
#include <ERF_IndexDefines.H>
#include <ERF_Morrison.H>
#include <ERF_Plotfile2DPrecip.H>
#include <ERF_SurfacePrecipitation.H>
#include <ERF_WDM6.H>
#include <ERF_WSM6.H>

#include "../../MovingBodies/ERF_GTestThrowOnAbort.H"

using namespace amrex;

namespace {

const BoxArray& boxes3d ()
{
    static const BoxArray ba(Box(IntVect(0, 0, 0), IntVect(3, 3, 2)));
    return ba;
}

BoxArray boxes2d () { return BoxArray(Box(IntVect(0, 0, 0), IntVect(3, 3, 0))); }

// One map for every MultiFab here: the 2D output and its 3D sources must share box owners
const DistributionMapping& dmap ()
{
    static const DistributionMapping dm(boxes3d());
    return dm;
}

// Unit fallout of each species over one step: rain 1, snow 2, graupel 3, cloud ice 4.
// The accumulators then hold what the update lines of WSM6 and WDM6 add:
// rain_accum += rain, snow_accum += snow, graupel_accum += graupel and
// precip_accum += rain + snow + graupel + cloud ice.
constexpr Real kRain = 1, kSnow = 2, kGraupel = 3, kIce = 4;
constexpr Real kTotal = kRain + kSnow + kGraupel + kIce;

void fill_unit_fallout (MultiFab& rain, MultiFab& snow, MultiFab& graupel, MultiFab& total)
{
    rain.setVal(kRain);
    snow.setVal(kSnow);
    graupel.setVal(kGraupel);
    total.setVal(kTotal);
}

// Components 0 (liquid) and 1 (frozen) of the flood's split, and the Total phase
struct Split { Real liquid_min, liquid_max, frozen_min, frozen_max, total; };

Split split_of (const SurfacePrecipAccumulationSources& src)
{
    const DistributionMapping& dm = dmap();
    MultiFab both(boxes2d(), dm, 2, 0);
    surface_precip_liquid_frozen_accum_k0(both, src, 0);
    MultiFab total(boxes2d(), dm, 1, 0);
    surface_precip_accum_k0(total, src, 0, SurfacePrecipPhase::Total);
    return {both.min(0), both.max(0), both.min(1), both.max(1), total.max(0)};
}

Geometry scheme_geometry ()
{
    const RealBox rb({Real(0), Real(0), Real(0)}, {Real(400), Real(400), Real(300)});
    return Geometry(boxes3d().minimalBox(), &rb, CoordSys::cartesian, Array<int, 3>{0, 0, 0}.data());
}

SolverChoice scheme_choice (MoistureType t)
{
    SolverChoice sc{};
    sc.c_p = Cp_d;
    sc.rdOcp = RdoCp;
    sc.moisture_type = t;
    sc.compressibility_type = {CompressibilityType::FullyCompressible};
    sc.use_native_shoc = false;
    sc.ave_plane = 2;
    return sc;
}

// Initialise a scheme on the test boxes and set its accumulators to one step of unit fallout
template <typename Scheme>
SurfacePrecipAccumulationSources unit_fallout_sources (Scheme& scheme, MoistureType t)
{
    const DistributionMapping& dm = dmap();
    MultiFab cons(boxes3d(), dm, RhoQ11_comp + 1, 1);
    cons.setVal(Real(1));
    std::unique_ptr<MultiFab> z_phys_nd, detJ_cc;
    scheme.SetCurrentLevel(0);
    SolverChoice sc = scheme_choice(t);
    scheme.Define(sc);
    scheme.Init(cons, boxes3d(), scheme_geometry(), Real(1), z_phys_nd, detJ_cc);
    // Qmoist order of all three schemes: rain, snow, graupel, precipitation total
    fill_unit_fallout(*scheme.Qmoist_Ptr(0), *scheme.Qmoist_Ptr(1),
                      *scheme.Qmoist_Ptr(2), *scheme.Qmoist_Ptr(3));
    return scheme.Get_Surface_Precip_Accumulation_Ptrs(0);
}

// Component c of the k = 0 plane of a 3D field, on the 2D boxes. A free function
// because nvcc rejects an extended __device__ lambda inside a gtest TEST body.
MultiFab k0_component (const MultiFab& src, const int c)
{
    MultiFab k0(boxes2d(), dmap(), 1, 0);
    for (MFIter mfi(k0); mfi.isValid(); ++mfi) {
        const auto o = k0.array(mfi);
        const auto d = src.const_array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { o(i, j, k) = d(i, j, 0, c); });
    }
    return k0;
}

} // namespace

// With rest_is_frozen the frozen part is the total less the rain, so the cloud ice
// counts as frozen, and the liquid part is the rain. The total is unchanged.
TEST(SurfacePrecipCloudIce, RestIsFrozenPutsTheCloudIceInTheFrozenPart)
{
    const DistributionMapping& dm = dmap();
    MultiFab total(boxes3d(), dm, 1, 0), rain(boxes3d(), dm, 1, 0);
    MultiFab snow(boxes3d(), dm, 1, 0), graupel(boxes3d(), dm, 1, 0);
    fill_unit_fallout(rain, snow, graupel, total);
    SurfacePrecipAccumulationSources src;
    src.total   = {&total, Real(1)};
    src.rain    = {&rain, Real(1)};
    src.snow    = {&snow, Real(1)};
    src.graupel = {&graupel, Real(1)};

    // Without the flag the rest of the total is liquid: liquid = rain + ice
    {
        const Split s = split_of(src);
        EXPECT_EQ(s.liquid_min, kRain + kIce);
        EXPECT_EQ(s.frozen_max, kSnow + kGraupel);
        EXPECT_EQ(s.total, kTotal);
    }
    src.rest_is_frozen = true;
    const Split s = split_of(src);
    EXPECT_EQ(s.liquid_min, kRain);
    EXPECT_EQ(s.liquid_max, kRain);
    EXPECT_EQ(s.frozen_min, kSnow + kGraupel + kIce);
    EXPECT_EQ(s.frozen_max, kSnow + kGraupel + kIce);
    EXPECT_EQ(s.total, kTotal);

    // Each phase on its own gives the same values as the two-component pass
    MultiFab out(boxes2d(), dm, 1, 0);
    surface_precip_accum_k0(out, src, 0, SurfacePrecipPhase::Liquid);
    EXPECT_EQ(out.max(0), kRain);
    surface_precip_accum_k0(out, src, 0, SurfacePrecipPhase::Frozen);
    EXPECT_EQ(out.min(0), kSnow + kGraupel + kIce);

    // The flag has no effect without a total (the species are then all there is)
    SurfacePrecipAccumulationSources species = src;
    species.total.accum = nullptr;
    const Split n = split_of(species);
    EXPECT_EQ(n.liquid_max, kRain);
    EXPECT_EQ(n.frozen_max, kSnow + kGraupel);
    EXPECT_EQ(n.total, kRain + kSnow + kGraupel);

    // The frozen part is the total less the rain, so the flag needs a rain accumulator
    SurfacePrecipAccumulationSources no_rain = src;
    no_rain.rain.accum = nullptr;
    const std::string msg = erf_gtest::abort_message([&] { split_of(no_rain); });
    EXPECT_NE(msg.find("rest_is_frozen needs a rain accumulator"), std::string::npos) << msg;
}

// The 2D plotfile's precip_frozen_accum follows the same rule; its rain and total do not change
TEST(SurfacePrecipCloudIce, PlotfileFrozenAccumulationHoldsTheCloudIce)
{
    const DistributionMapping& dm = dmap();
    MultiFab total(boxes3d(), dm, 1, 0), rain(boxes3d(), dm, 1, 0);
    MultiFab snow(boxes3d(), dm, 1, 0), graupel(boxes3d(), dm, 1, 0);
    fill_unit_fallout(rain, snow, graupel, total);
    SurfacePrecipAccumulationSources src;
    src.total   = {&total, Real(1)};
    src.rain    = {&rain, Real(1)};
    src.snow    = {&snow, Real(1)};
    src.graupel = {&graupel, Real(1)};
    src.rest_is_frozen = true;

    const Vector<std::string> names{"precip_total_accum", "precip_rain_accum", "precip_frozen_accum"};
    const auto selected = plotfile2d::selected_precipitation_accumulation_components(names, src);
    ASSERT_EQ(selected.n, 3);
    MultiFab dst(boxes3d(), dm, 3, 0);
    dst.setVal(Real(-1));
    plotfile2d::fill_precipitation_accumulations(dst, src, selected, 0);
    for (int n = 0; n < selected.n; ++n) {
        Real expected = kTotal;
        if (selected.id[n] == plotfile2d::DiagnosticID::PrecipRainAccum)   { expected = kRain; }
        if (selected.id[n] == plotfile2d::DiagnosticID::PrecipFrozenAccum) { expected = kSnow + kGraupel + kIce; }
        const MultiFab k0 = k0_component(dst, selected.dst_comp[n]);
        EXPECT_EQ(k0.min(0), expected) << names[n];
        EXPECT_EQ(k0.max(0), expected) << names[n];
    }

    // As in the flood's helper, the flag needs a rain accumulator
    SurfacePrecipAccumulationSources no_rain = src;
    no_rain.rain.accum = nullptr;
    const Vector<std::string> frozen_only{"precip_frozen_accum"};
    const auto sel = plotfile2d::selected_precipitation_accumulation_components(frozen_only, no_rain);
    const std::string msg = erf_gtest::abort_message(
        [&] { plotfile2d::fill_precipitation_accumulations(dst, no_rain, sel, 0); });
    EXPECT_NE(msg.find("rest_is_frozen needs a rain accumulator"), std::string::npos) << msg;
}

// WSM6 and WDM6: one step of unit fallout of each species gives the rain as the
// liquid part and snow, graupel and cloud ice as the frozen part
TEST(SurfacePrecipCloudIce, WSM6AndWDM6CountTheirCloudIceAsFrozen)
{
    {
        WSM6 wsm6;
        const SurfacePrecipAccumulationSources src = unit_fallout_sources(wsm6, MoistureType::WSM6);
        EXPECT_TRUE(src.rest_is_frozen);
        const Split s = split_of(src);
        EXPECT_EQ(s.liquid_max, kRain * src.rain.native_to_kg_m2);
        EXPECT_EQ(s.frozen_min, (kTotal - kRain) * src.total.native_to_kg_m2);
        EXPECT_EQ(s.total, kTotal * src.total.native_to_kg_m2);
    }
    {
        WDM6 wdm6;
        const SurfacePrecipAccumulationSources src = unit_fallout_sources(wdm6, MoistureType::WDM6);
        EXPECT_TRUE(src.rest_is_frozen);
        const Split s = split_of(src);
        EXPECT_EQ(s.liquid_max, kRain * src.rain.native_to_kg_m2);
        EXPECT_EQ(s.frozen_min, (kTotal - kRain) * src.total.native_to_kg_m2);
        EXPECT_EQ(s.total, kTotal * src.total.native_to_kg_m2);
    }
}

// Morrison's total also holds the cloud water that fell out, which is liquid, so
// the rest of its total is not all frozen and keeps the old rule
TEST(SurfacePrecipCloudIce, MorrisonKeepsTheRestInTheLiquidPart)
{
    Morrison morr;
    const SurfacePrecipAccumulationSources src = unit_fallout_sources(morr, MoistureType::Morrison);
    EXPECT_FALSE(src.rest_is_frozen);
    // Liquid is the total less snow and graupel, so it keeps the rest (here the unit cloud ice)
    // (Morrison's unit factor need not be 1, so allow a few units of round-off.)
    const Split s = split_of(src);
    const Real f = src.total.native_to_kg_m2;
    const Real tol = Real(8) * std::numeric_limits<Real>::epsilon() * kTotal * f;
    EXPECT_NEAR(s.liquid_min, (kTotal - kSnow - kGraupel) * f, tol);
    EXPECT_NEAR(s.liquid_max, (kTotal - kSnow - kGraupel) * f, tol);
    EXPECT_NEAR(s.frozen_max, (kSnow + kGraupel) * f, tol);
    EXPECT_NEAR(s.total, kTotal * f, tol);
}

#include <cmath>
#include <utility>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <gtest/gtest.h>

#include <ERF_TwoStreamColumn.H>
#include <ERF_TwoStreamRadiation.H>

// Two-stream -> land-surface forcing contract
// -------------------------------------------
// What the two-stream model hands Noah-MP (SWDOWN, GLW, COSZEN) is read off the
// interface fluxes the sweep writes and copied into the land model's k = 0 plane:
//   1. The column helper takes the surface interface (index kmin) of component 1
//      (downwelling SW) and component 3 (downwelling LW) -- not the net, not the
//      upwelling, not another level.
//   2. The cosine of the zenith angle is floored at zero, and a dynamic sun uses the
//      column's own latitude and longitude.
//   3. copy_surface_plane writes exactly the valid i,j cells of the destination's k = 0
//      plane, whether that plane is valid (Noah-MP) or a z ghost (SLM-style soil
//      layers below k = 0), and leaves every ghost column and every other k alone.

namespace {

constexpr amrex::Real kSentinel = amrex::Real(-999.0);

// Number of cells, ghost cells included, holding exactly `value`, over all ranks.
amrex::Long count_equal (const amrex::MultiFab& mf, amrex::Real value)
{
    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Long> reduce_data(reduce_op);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto array = mf.const_array(mfi);
        reduce_op.eval(mf[mfi].box(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Long>
            {
                return {array(i, j, k) == value ? amrex::Long(1) : amrex::Long(0)};
            });
    }
    amrex::Gpu::streamSynchronize();
    amrex::Long count = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceLongSum(count);
    return count;
}

// Largest |dst(i,j,0) - src(i,j,0)| over the valid i,j extent of dst's boxes, over all ranks.
amrex::Real max_plane_difference (const amrex::MultiFab& dst, const amrex::MultiFab& src)
{
    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    for (amrex::MFIter mfi(dst); mfi.isValid(); ++mfi) {
        const auto d = dst.const_array(mfi);
        const auto s = src.const_array(mfi);
        reduce_op.eval(amrex::makeSlab(mfi.validbox(), 2, 0), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Real>
            {
                return {std::abs(d(i, j, k) - s(i, j, k))};
            });
    }
    amrex::Gpu::streamSynchronize();
    amrex::Real diff = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceRealMax(diff);
    return diff;
}

void fill_ramp (amrex::MultiFab& mf)
{
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto array = mf.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            array(i, j, k) = amrex::Real(1.0) + amrex::Real(i) + amrex::Real(100.0) * amrex::Real(j);
        });
    }
}

// Two boxes side by side in x, collapsed to the k = 0 plane as ERF's ba2d is.
std::pair<amrex::BoxArray, amrex::BoxArray> two_box_layout ()
{
    amrex::BoxList boxes;
    boxes.push_back(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(3, 2, 0)));
    boxes.push_back(amrex::Box(amrex::IntVect(4, 0, 0), amrex::IntVect(7, 2, 0)));
    const amrex::BoxArray ba2d(boxes);
    // The same columns with soil layers k = -3..-1 below the surface plane.
    amrex::BoxList soil = boxes;
    for (auto& b : soil) { b.setRange(2, -3, 3); }
    return {ba2d, amrex::BoxArray(soil)};
}

} // namespace

TEST(TwoStreamLandForcing, ColumnTakesDownwellingFluxesAtTheSurfaceInterface)
{
    // One column, interfaces k = 2..6 (the surface is kmin = 2), each component and
    // each level holding a different value, so reading any other level or component
    // gives a wrong answer.
    const int kmin = 2;
    const int kmax_iface = 6;
    std::vector<amrex::Real> storage((kmax_iface - kmin + 1) * 4);
    const amrex::Array4<amrex::Real> flux(storage.data(), amrex::Dim3{0, 0, kmin},
                                          amrex::Dim3{1, 1, kmax_iface + 1}, 4);
    for (int k = kmin; k <= kmax_iface; ++k) {
        for (int comp = 0; comp < 4; ++comp) {
            flux(0, 0, k, comp) = amrex::Real(1000 * (comp + 1) + 10 * k);
        }
    }
    const amrex::Array4<const amrex::Real> flux_c(flux);

    amrex::Real sw_dn = 0.0, lw_dn = 0.0, coszen = 0.0;
    two_stream_land_forcing(0, 0, kmin, flux_c, amrex::Real(0.4), sw_dn, lw_dn, coszen);
    EXPECT_EQ(sw_dn, flux(0, 0, kmin, 1));   // 2020: SW down at the surface
    EXPECT_EQ(lw_dn, flux(0, 0, kmin, 3));   // 4020: LW down at the surface
    EXPECT_EQ(coszen, amrex::Real(0.4));

    // Sun below the horizon: COSZEN is zero, not negative.
    two_stream_land_forcing(0, 0, kmin, flux_c, amrex::Real(-0.3), sw_dn, lw_dn, coszen);
    EXPECT_EQ(coszen, amrex::Real(0.0));
}

TEST(TwoStreamLandForcing, ZenithAngleUsesTheColumnsOwnSite)
{
    TwoStreamParams p;
    p.cos_zenith_fixed = amrex::Real(0.37);
    EXPECT_EQ(two_stream_cos_zenith(0, 0, p, false, nullptr, nullptr), amrex::Real(0.37));

    // Dynamic sun at 12:00 UTC on the equinox-like declination 0. Column 0 sits on the
    // equator at the Greenwich meridian (sun overhead), column 1 on the equator at the
    // antimeridian (midnight). The constant site differs in both coordinates -- 60 N,
    // 90 E, where the sun is setting -- so a helper that took either coordinate from the
    // constants instead of the column would miss the overhead sun.
    p.solar_dynamic = true;
    p.calday = amrex::Real(1.5);
    p.declin = amrex::Real(0.0);
    p.lat_cons_rad = amrex::Real(PI) / amrex::Real(3.0);
    p.lon_cons_rad = amrex::Real(0.5) * amrex::Real(PI);
    std::vector<amrex::Real> lat_v = {0.0, 0.0};
    std::vector<amrex::Real> lon_v = {0.0, 180.0};
    const amrex::Array4<const amrex::Real> lat(lat_v.data(), amrex::Dim3{0, 0, 0}, amrex::Dim3{2, 1, 1}, 1);
    const amrex::Array4<const amrex::Real> lon(lon_v.data(), amrex::Dim3{0, 0, 0}, amrex::Dim3{2, 1, 1}, 1);

    const amrex::Real noon = two_stream_cos_zenith(0, 0, p, true, &lat, &lon);
    const amrex::Real midnight = two_stream_cos_zenith(1, 0, p, true, &lat, &lon);
    EXPECT_EQ(noon, static_cast<amrex::Real>(orbital_cos_zenith_instant(1.5, 0.0, 0.0, 0.0)));
    EXPECT_GT(noon, amrex::Real(0.99));
    EXPECT_LT(midnight, amrex::Real(-0.99));

    // Without per-column arrays the constant site applies.
    const amrex::Real constant_site = two_stream_cos_zenith(0, 0, p, false, nullptr, nullptr);
    EXPECT_NEAR(constant_site, amrex::Real(0.0), amrex::Real(1.0e-5));

    // And the land forcing floors the midnight column.
    std::vector<amrex::Real> flux_v(4, amrex::Real(1.0));
    const amrex::Array4<const amrex::Real> flux(flux_v.data(), amrex::Dim3{1, 0, 0}, amrex::Dim3{2, 1, 1}, 4);
    amrex::Real sw_dn = 0.0, lw_dn = 0.0, coszen = 1.0;
    two_stream_land_forcing(1, 0, 0, flux, midnight, sw_dn, lw_dn, coszen);
    EXPECT_EQ(coszen, amrex::Real(0.0));
}

TEST(TwoStreamLandForcing, CopySurfacePlaneFillsNoahMPLayout)
{
    const auto [ba2d, ba_soil] = two_box_layout();
    const amrex::DistributionMapping dm(ba2d);

    amrex::MultiFab src(ba2d, dm, 1, 0);
    fill_ramp(src);

    // Noah-MP's lsm_fab_data: the same 2D boxes, one ghost cell in x and y.
    amrex::MultiFab dst(ba2d, dm, 1, amrex::IntVect(1, 1, 0));
    dst.setVal(kSentinel);
    copy_surface_plane(src, dst);

    EXPECT_EQ(max_plane_difference(dst, src), amrex::Real(0.0));
    // Every ghost cell still holds the sentinel, and no valid cell does.
    amrex::Long grown_cells = 0;
    for (amrex::MFIter mfi(dst); mfi.isValid(); ++mfi) { grown_cells += dst[mfi].box().numPts(); }
    amrex::ParallelDescriptor::ReduceLongSum(grown_cells);
    EXPECT_EQ(count_equal(dst, kSentinel), grown_cells - ba2d.numPts());
}

TEST(TwoStreamLandForcing, CopySurfacePlaneFillsGhostSurfacePlane)
{
    const auto [ba2d, ba_soil] = two_box_layout();
    const amrex::DistributionMapping dm(ba2d);

    amrex::MultiFab src(ba2d, dm, 1, 0);
    fill_ramp(src);

    // SLM-style storage: soil layers k = -3..-1 valid, the surface plane k = 0 in the
    // z ghost region.
    amrex::MultiFab dst(ba_soil, dm, 1, amrex::IntVect(1, 1, 1));
    dst.setVal(kSentinel);
    copy_surface_plane(src, dst);

    EXPECT_EQ(max_plane_difference(dst, src), amrex::Real(0.0));
    // Only the valid i,j cells of the k = 0 plane changed: not the soil, not the k = -4
    // ghost plane, not the x/y ghost columns of the surface plane.
    amrex::Long grown_cells = 0;
    for (amrex::MFIter mfi(dst); mfi.isValid(); ++mfi) { grown_cells += dst[mfi].box().numPts(); }
    amrex::ParallelDescriptor::ReduceLongSum(grown_cells);
    EXPECT_EQ(count_equal(dst, kSentinel), grown_cells - ba2d.numPts());
}

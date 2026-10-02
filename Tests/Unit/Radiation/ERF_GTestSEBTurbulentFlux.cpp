#include <cmath>
#include <memory>
#include <utility>

#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <gtest/gtest.h>

#include <ERF_Constants.H>
#include <ERF_RadStruct.H>
#include <ERF_SEBTurbulentFlux.H>

// The two-stream surface energy balance takes its sensible (H) and latent (LE)
// heat fluxes from a land-surface field when one exists, else from the flux the
// surface layer applies to the air, else from the scalar default. Before this,
// without a land model, H and LE were always the defaults (0), so the ground kept
// the heat the surface layer was putting into the air.

namespace {

using amrex::Real;

constexpr Real tol = sizeof(Real) == 8 ? Real(1.0e-12) : Real(1.0e-5);

struct Layout
{
    amrex::BoxArray ba3d;
    amrex::BoxArray ba2d;
    amrex::DistributionMapping dm;
};

// A 4 x 3 x 5 domain in two boxes, and the matching 2D (k = 0) layout the SEB
// fields use.
Layout make_layout ()
{
    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(3, 2, 4));
    amrex::BoxArray ba(domain);
    ba.maxSize(amrex::IntVect(2, 3, 5));
    amrex::BoxList boxes = ba.boxList();
    for (auto& box : boxes) { box.setRange(2, 0); }
    amrex::DistributionMapping dm(ba);
    return {ba, amrex::BoxArray(std::move(boxes)), dm};
}

// The conservative flux on z faces: a different value on every face and column,
// negative in some columns, so the test sees the face index, the column and the sign.
AMREX_GPU_HOST_DEVICE Real face_flux (int i, int j, int k)
{
    return Real(0.01) * Real(i - 1) + Real(0.003) * Real(j) + Real(0.5) * Real(k);
}

std::unique_ptr<amrex::MultiFab> make_surface_layer_flux (const Layout& l)
{
    auto mf = std::make_unique<amrex::MultiFab>(amrex::convert(l.ba3d, amrex::IntVect(0, 0, 1)), l.dm, 1, 0);
    for (amrex::MFIter mfi(*mf); mfi.isValid(); ++mfi) {
        const auto arr = mf->array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            arr(i, j, k) = face_flux(i, j, k);
        });
    }
    return mf;
}

// Largest |seb(i,j,0) - scale * face_flux(i,j,k_face)| over the valid cells of
// every rank.
Real max_error (const amrex::MultiFab& seb, Real scale, int k_face)
{
    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<Real> reduce_data(reduce_op);
    for (amrex::MFIter mfi(seb); mfi.isValid(); ++mfi) {
        const auto arr = seb.const_array(mfi);
        reduce_op.eval(mfi.validbox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int) -> amrex::GpuTuple<Real> {
                return {std::abs(arr(i, j, 0) - scale * face_flux(i, j, k_face))};
            });
    }
    Real err = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceRealMax(err);
    return err;
}

// A land-surface field as a land model leaves it: a value in most columns and the
// lsm_undefined placeholder in the columns it did not process (here every third
// one, standing in for water and sea ice). `offset` makes two such fields differ.
AMREX_GPU_HOST_DEVICE bool lsm_column_processed (int i, int j, int skip)
{
    return (i + 2 * j + skip) % 3 != 0;
}

AMREX_GPU_HOST_DEVICE Real lsm_value (int i, int j, Real offset)
{
    return offset + Real(1.5) * Real(i) - Real(0.75) * Real(j);
}

std::unique_ptr<amrex::MultiFab> make_lsm_field (const Layout& l, Real offset, int skip)
{
    auto mf = std::make_unique<amrex::MultiFab>(l.ba2d, l.dm, 1, 0);
    for (amrex::MFIter mfi(*mf); mfi.isValid(); ++mfi) {
        const auto arr = mf->array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            arr(i, j, k) = lsm_column_processed(i, j, skip) ? lsm_value(i, j, offset) : lsm_undefined;
        });
    }
    return mf;
}

// Over the valid cells of every rank: the largest error against the expected
// land-surface fill, and the number of cells that hold `fallback`.
struct LandFillCheck { Real max_error; amrex::Long fallback_cells; };

LandFillCheck check_land_fill (const amrex::MultiFab& seb, Real scale, Real off_a, int skip_a,
                               bool with_b, Real off_b, int skip_b, Real fallback)
{
    amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<Real, amrex::Long> reduce_data(reduce_op);
    for (amrex::MFIter mfi(seb); mfi.isValid(); ++mfi) {
        const auto arr = seb.const_array(mfi);
        reduce_op.eval(mfi.validbox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int) -> amrex::GpuTuple<Real, amrex::Long> {
                const bool valid = lsm_column_processed(i, j, skip_a) &&
                                   (!with_b || lsm_column_processed(i, j, skip_b));
                const Real expected = valid ? scale * lsm_value(i, j, off_a) +
                                              (with_b ? lsm_value(i, j, off_b) : Real(0.0))
                                            : fallback;
                return {std::abs(arr(i, j, 0) - expected), arr(i, j, 0) == fallback ? 1L : 0L};
            });
    }
    auto result = reduce_data.value();
    Real err = amrex::get<0>(result);
    amrex::Long n = amrex::get<1>(result);
    amrex::ParallelDescriptor::ReduceRealMax(err);
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return {err, n};
}

// Largest error of a land-over-surface-layer H fill, over the valid cells of every
// rank: the land value where the land model processed the column, Cp_d times the
// surface face flux elsewhere.
Real max_error_land_over_surface_layer (const amrex::MultiFab& seb, Real offset)
{
    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<Real> reduce_data(reduce_op);
    for (amrex::MFIter mfi(seb); mfi.isValid(); ++mfi) {
        const auto arr = seb.const_array(mfi);
        reduce_op.eval(mfi.validbox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int) -> amrex::GpuTuple<Real> {
                const Real expected = lsm_column_processed(i, j, 0) ? lsm_value(i, j, offset)
                                                                     : Cp_d * face_flux(i, j, 0);
                return {std::abs(arr(i, j, 0) - expected)};
            });
    }
    Real err = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceRealMax(err);
    return err;
}

// Largest |mf - value| over every cell of every fab, halo included, on every rank.
Real max_deviation_with_halo (const amrex::MultiFab& mf, Real value)
{
    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<Real> reduce_data(reduce_op);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto arr = mf.const_array(mfi);
        reduce_op.eval(mfi.fabbox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<Real> {
                return {std::abs(arr(i, j, k) - value)};
            });
    }
    Real err = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceRealMax(err);
    return err;
}

// Largest |mf - value| over the halo cells only (every cell of each fab outside its
// valid box), on every rank.
Real max_halo_deviation (const amrex::MultiFab& mf, Real value)
{
    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<Real> reduce_data(reduce_op);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto arr = mf.const_array(mfi);
        const amrex::Box valid = mfi.validbox();
        reduce_op.eval(mfi.fabbox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<Real> {
                if (valid.contains(i, j, k)) { return {Real(0.0)}; }
                return {std::abs(arr(i, j, k) - value)};
            });
    }
    Real err = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceRealMax(err);
    return err;
}

} // namespace

TEST(SEBTurbulentFlux, PrecedenceLandSurfaceThenSurfaceLayerThenDefault)
{
    using S = SEBTurbulentFluxSource;
    for (const S source : {S::SurfaceLayer, S::Defaults}) {
        for (const bool sl : {false, true}) {
            // A land-surface field wins whatever else is available.
            EXPECT_EQ(select_seb_flux_origin(true, source, sl), SEBFluxOrigin::LandSurface);
        }
    }
    EXPECT_EQ(select_seb_flux_origin(false, S::SurfaceLayer, true), SEBFluxOrigin::SurfaceLayer);
    // No surface-layer flux (no zlo surface layer, or no moisture model for LE).
    EXPECT_EQ(select_seb_flux_origin(false, S::SurfaceLayer, false), SEBFluxOrigin::Default);
    // seb_turbulent_flux_source = defaults restores the constants.
    EXPECT_EQ(select_seb_flux_origin(false, S::Defaults, true), SEBFluxOrigin::Default);
    EXPECT_EQ(select_seb_flux_origin(false, S::Defaults, false), SEBFluxOrigin::Default);
}

TEST(SEBTurbulentFlux, SensibleFromSurfaceLayerIsCpTimesTheSurfaceFaceFlux)
{
    const Layout l = make_layout();
    auto sl = make_surface_layer_flux(l);
    amrex::MultiFab seb(l.ba2d, l.dm, 1, amrex::IntVect(1, 1, 0));
    seb.setVal(Real(-12345.0));

    const auto origin = fill_seb_turbulent_flux(seb, SEBTurbulentFlux::Sensible, nullptr, sl.get(), 0,
                                                SEBTurbulentFluxSource::SurfaceLayer, Real(7.0));
    EXPECT_EQ(origin, SEBFluxOrigin::SurfaceLayer);
    // W/m^2 = Cp_d * rho w'theta' on the surface face (k = 0), sign kept: the
    // columns with i = 0 carry a downward flux and must come out negative.
    EXPECT_LT(max_error(seb, Cp_d, 0), tol * Cp_d);
    EXPECT_LT(seb.min(0), Real(0.0));
    EXPECT_GT(seb.max(0), Real(0.0));
}

TEST(SEBTurbulentFlux, LatentFromSurfaceLayerIsLvTimesTheSurfaceFaceFlux)
{
    const Layout l = make_layout();
    auto sl = make_surface_layer_flux(l);
    amrex::MultiFab seb(l.ba2d, l.dm, 1, amrex::IntVect(1, 1, 0));
    seb.setVal(Real(-12345.0));

    fill_seb_turbulent_flux(seb, SEBTurbulentFlux::Latent, nullptr, sl.get(), 0,
                            SEBTurbulentFluxSource::SurfaceLayer, Real(7.0));
    EXPECT_LT(max_error(seb, L_v, 0), tol * L_v);
}

TEST(SEBTurbulentFlux, SurfaceIndexSelectsTheFace)
{
    // surface_k picks the face; the field varies with k, so reading the wrong
    // face is off by 0.5 * Cp_d.
    const Layout l = make_layout();
    auto sl = make_surface_layer_flux(l);
    amrex::MultiFab seb(l.ba2d, l.dm, 1, 0);
    fill_seb_turbulent_flux(seb, SEBTurbulentFlux::Sensible, nullptr, sl.get(), 2,
                            SEBTurbulentFluxSource::SurfaceLayer, Real(7.0));
    EXPECT_LT(max_error(seb, Cp_d, 2), tol * Cp_d);
}

TEST(SEBTurbulentFlux, LandSurfaceFieldOverridesTheSurfaceLayer)
{
    const Layout l = make_layout();
    auto sl = make_surface_layer_flux(l);
    amrex::MultiFab lsm(l.ba2d, l.dm, 1, 0);
    lsm.setVal(Real(42.5));
    amrex::MultiFab seb(l.ba2d, l.dm, 1, 0);

    const auto origin = fill_seb_turbulent_flux(seb, SEBTurbulentFlux::Sensible, &lsm, sl.get(), 0,
                                                SEBTurbulentFluxSource::SurfaceLayer, Real(7.0));
    EXPECT_EQ(origin, SEBFluxOrigin::LandSurface);
    EXPECT_EQ(seb.min(0), Real(42.5));
    EXPECT_EQ(seb.max(0), Real(42.5));
}

TEST(SEBTurbulentFlux, NoMoistureModelLeavesLatentAtItsDefault)
{
    // Without a moisture model there is no moisture flux field, so ERF hands the
    // SEB a null latent source and LE stays at seb_lh_default.
    const Layout l = make_layout();
    amrex::MultiFab seb(l.ba2d, l.dm, 1, 0);
    seb.setVal(Real(-12345.0));
    const auto origin = fill_seb_turbulent_flux(seb, SEBTurbulentFlux::Latent, nullptr, nullptr, 0,
                                                SEBTurbulentFluxSource::SurfaceLayer, Real(20.0));
    EXPECT_EQ(origin, SEBFluxOrigin::Default);
    EXPECT_EQ(seb.min(0), Real(20.0));
    EXPECT_EQ(seb.max(0), Real(20.0));
}

TEST(SEBTurbulentFlux, DefaultsSourceIgnoresTheSurfaceLayer)
{
    const Layout l = make_layout();
    auto sl = make_surface_layer_flux(l);
    amrex::MultiFab seb(l.ba2d, l.dm, 1, 0);
    const auto origin = fill_seb_turbulent_flux(seb, SEBTurbulentFlux::Sensible, nullptr, sl.get(), 0,
                                                SEBTurbulentFluxSource::Defaults, Real(10.0));
    EXPECT_EQ(origin, SEBFluxOrigin::Default);
    EXPECT_EQ(seb.min(0), Real(10.0));
    EXPECT_EQ(seb.max(0), Real(10.0));
}

// The two-stream model reads its surface inputs (absorbed SW = sav + sag, net LW =
// -fira, ground flux, surface humidity, albedo) from the land model. A land model
// leaves the lsm_undefined placeholder where it did not compute a value (Noah-MP:
// open water and sea ice, and every cell before its first step); a blanket copy put
// ~1e150 into the balance there. Those cells take the scalar default instead.
TEST(SEBLandSurfaceField, UndefinedCellsTakeTheDefault)
{
    const Layout l = make_layout();
    auto fira = make_lsm_field(l, Real(40.0), 0);
    amrex::MultiFab seb(l.ba2d, l.dm, 1, amrex::IntVect(1, 1, 0));
    fill_seb_field_from_land_surface(seb, fira.get(), Real(-25.0), Real(-1.0));
    const LandFillCheck c = check_land_fill(seb, Real(-1.0), Real(40.0), 0, false, 0.0, 0, Real(-25.0));
    EXPECT_LT(c.max_error, tol * Real(100.0));
    // 12 columns, (i + 2j) % 3 == 0 in 4 of them.
    EXPECT_EQ(c.fallback_cells, 4);
    EXPECT_LT(seb.max(0), Real(1.0e3));
}

TEST(SEBLandSurfaceField, SumNeedsBothFieldsValid)
{
    // sav + sag: a column with either part undefined takes the default.
    const Layout l = make_layout();
    auto sav = make_lsm_field(l, Real(300.0), 0);
    auto sag = make_lsm_field(l, Real(120.0), 1);
    amrex::MultiFab seb(l.ba2d, l.dm, 1, 0);
    fill_seb_field_from_land_surface(seb, sav.get(), Real(50.0), Real(1.0), sag.get());
    const LandFillCheck c = check_land_fill(seb, Real(1.0), Real(300.0), 0, true, Real(120.0), 1, Real(50.0));
    EXPECT_LT(c.max_error, tol * Real(1000.0));
    // Undefined in either: 4 columns each, none shared, so 8 of 12.
    EXPECT_EQ(c.fallback_cells, 8);
}

TEST(SEBLandSurfaceField, NoFieldMeansTheDefaultEverywhere)
{
    const Layout l = make_layout();
    amrex::MultiFab seb(l.ba2d, l.dm, 1, amrex::IntVect(1, 1, 0));
    seb.setVal(Real(-1.0));
    fill_seb_field_from_land_surface(seb, nullptr, Real(0.3));
    EXPECT_EQ(max_deviation_with_halo(seb, Real(0.3)), Real(0.0));
}

TEST(SEBTurbulentFlux, UndefinedLandCellsFallToTheSurfaceLayer)
{
    // A land-model H that is undefined over water: those columns take the surface
    // layer's flux, the others the land model's.
    const Layout l = make_layout();
    auto sl = make_surface_layer_flux(l);
    auto land_h = make_lsm_field(l, Real(80.0), 0);
    amrex::MultiFab seb(l.ba2d, l.dm, 1, 0);
    const auto origin = fill_seb_turbulent_flux(seb, SEBTurbulentFlux::Sensible, land_h.get(), sl.get(), 0,
                                                SEBTurbulentFluxSource::SurfaceLayer, Real(7.0));
    EXPECT_EQ(origin, SEBFluxOrigin::LandSurface);
    EXPECT_LT(max_error_land_over_surface_layer(seb, Real(80.0)), tol * Cp_d);
}

// The documented halo contract: whatever the source, the halo holds the fallback and
// only valid cells take the selected values. Before, the surface-layer case left the
// halo as it was and the other cases set it, so the contract held in some cases only.
TEST(SEBTurbulentFlux, HaloHoldsTheFallbackInEveryCase)
{
    using S = SEBTurbulentFluxSource;
    const Layout l = make_layout();
    auto sl = make_surface_layer_flux(l);
    auto land_h = make_lsm_field(l, Real(80.0), 0);
    const Real fallback = Real(7.0);
    for (const amrex::MultiFab* land : {static_cast<const amrex::MultiFab*>(nullptr),
                                        static_cast<const amrex::MultiFab*>(land_h.get())}) {
        for (const S source : {S::SurfaceLayer, S::Defaults}) {
            amrex::MultiFab seb(l.ba2d, l.dm, 1, amrex::IntVect(1, 1, 0));
            seb.setVal(Real(-12345.0));
            fill_seb_turbulent_flux(seb, SEBTurbulentFlux::Sensible, land, sl.get(), 0, source, fallback);
            EXPECT_EQ(max_halo_deviation(seb, fallback), Real(0.0))
                << "land field " << (land != nullptr) << ", source " << static_cast<int>(source);
        }
    }

    amrex::MultiFab seb(l.ba2d, l.dm, 1, amrex::IntVect(1, 1, 0));
    seb.setVal(Real(-12345.0));
    fill_seb_field_from_land_surface(seb, land_h.get(), Real(0.3));
    EXPECT_EQ(max_halo_deviation(seb, Real(0.3)), Real(0.0));
}

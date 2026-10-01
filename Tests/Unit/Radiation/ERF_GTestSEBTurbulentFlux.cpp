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

// Largest |seb(i,j,0) - scale * face_flux(i,j,k_face)| over the valid cells.
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
    return amrex::get<0>(reduce_data.value());
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

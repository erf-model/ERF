#include <cmath>
#include <string>
#include <utility>

#include <AMReX_BCRec.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Reduce.H>

#include <gtest/gtest.h>

#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_NumericalConstants.H>
#include <ERF_LandSurface.H>
#include <ERF_RadStruct.H>
#include <ERF_SEBLandAverageDown.H>
#include <ERF_SEBLandFill.H>
#include <ERF_SEBLandMask.H>
#include <ERF_SurfacePrecipitation.H>
#include <ERF_SurfaceTemperature.H>
#include <ERF_TwoStreamLW.H>
#include <ERF_TwoStreamRadiation.H>

// The land mask in the two-stream force-restore balance (the land_mask argument of
// TwoStreamRadiation::advance()). A 4 x 4 level of 100 m columns in four boxes, land
// where i < 2 and water where i >= 2; two boxes hold land and water, two only water; one
// land column carries the mask of a building's column (2). The balance is a model of
// soil: a water column keeps its state, and the sweep's longwave boundary there is the
// surface layer's sea surface temperature, sigma T_sea^4 with an emissivity of one, not
// the balance's skin, sigma T_skin^4.

namespace {

using namespace amrex;

BoxArray collapse_z (const BoxArray& ba)
{
    BoxList boxes = ba.boxList();
    for (auto& box : boxes) { box.setRange(2, 0); }
    return BoxArray(std::move(boxes));
}

Real value_at (const MultiFab& mf, const IntVect& point, int comp)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    const Box point_box(point, point, mf.boxArray().ixType());
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box overlap = point_box & mfi.validbox();
        if (overlap.isEmpty()) { continue; }
        const auto array = mf.const_array(mfi);
        reduce_op.eval(overlap, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<Real> { return {array(i, j, k, comp)}; });
    }
    Gpu::streamSynchronize();
    Real value = amrex::get<0>(reduce_data.value());
    ParallelDescriptor::ReduceRealSum(value);
    return value;
}

constexpr Real t_skin0 = Real(300.0);   // the balance's initial skin [K]
constexpr Real theta_sea = Real(285.0); // the surface layer's potential temperature [K]
constexpr Real t_air = Real(290.0);
constexpr Real dz = Real(100.0);
constexpr Real q_soil0 = Real(0.2);    // the soil water of the runs with rain [m^3/m^3]

RadChoice mask_choice ()
{
    RadChoice rad;
    rad.enabled = true;
    rad.sw_enabled = true;
    rad.lw_enabled = true;
    rad.fixed_solar_zenith_angle = Real(0.5);
    rad.fixed_total_solar_irradiance = Real(1000.0);
    rad.tau_per_layer = Real(0.05);
    rad.tau_lw_per_layer = Real(0.2);
    rad.single_scattering_albedo = Real(0.5);
    rad.surface_albedo_sw = Real(0.2);
    rad.surface_emissivity_lw = Real(1.0);
    rad.rad_t_sfc = t_skin0;
    rad.seb_enable = true;
    rad.seb_use_radiation_fluxes = true;
    rad.seb_prognostic_enable = true;
    rad.seb_diagnostic_enable = true;
    rad.seb_hfx_default = Real(0.0);
    rad.seb_lh_default = Real(0.0);
    rad.seb_grdflx_default = Real(0.0);
    rad.seb_t_deep_default = t_skin0;
    rad.seb_surface_heat_capacity = Real(20000.0);
    rad.seb_turbulent_flux_source = SEBTurbulentFluxSource::Defaults;
    return rad;
}

struct Level
{
    Box domain{IntVect(0, 0, 0), IntVect(3, 3, 3)};
    RealBox real_box{{0.0, 0.0, 0.0}, {400.0, 400.0, 400.0}};
    Geometry geom;
    BoxArray ba;
    BoxArray ba2d;
    DistributionMapping dm;
    Level ()
    {
        const int is_periodic[3] = {1, 1, 0};
        geom.define(domain, &real_box, 0, is_periodic);
        // Boxes of 3 and 1 columns in x and 2 in y: the shore (between i = 1 and 2) is
        // inside the wide boxes, which hold land and water; the narrow ones are all water.
        BoxList boxes;
        for (int j = 0; j < 4; j += 2) {
            boxes.push_back(Box(IntVect(0, j, 0), IntVect(2, j + 1, 3)));
            boxes.push_back(Box(IntVect(3, j, 0), IntVect(3, j + 1, 3)));
        }
        ba = BoxArray(std::move(boxes));
        ba2d = collapse_z(ba);
        dm = DistributionMapping(ba);
    }
};

// mode 0: land where i < 2 (a building's column at (0,3)), water where i >= 2;
// mode 1: every column land; mode 2: every column water; mode 3: as mode 0 with the
// building's column an ordinary land column (1).
void fill_mask (iMultiFab& mask, int mode)
{
    for (MFIter mfi(mask); mfi.isValid(); ++mfi) {
        const auto m = mask.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (mode == 1) { m(i, j, k) = 1; }
            else if (mode == 2) { m(i, j, k) = 0; }
            else { m(i, j, k) = (i >= 2) ? 0 : ((i == 0 && j == 3 && mode == 0) ? 2 : 1); }
        });
    }
    Gpu::streamSynchronize();
}

struct MaskRun
{
    Level L;
    RadChoice rad = mask_choice();   // the radiation object keeps a pointer to it
    TwoStreamRadiation radiation;
    MultiFab state, qheating, rad_fluxes, t_surf;
    LandSurface lsm;

    explicit MaskRun (Real cloud_fraction = Real(0.0), int verbosity = 0, bool snow_and_rain = false)
    {
        rad.cloud_fraction = cloud_fraction;
        rad.verbosity = verbosity;
        if (snow_and_rain) {
            rad.seb_snow_enable = true;
            rad.seb_precipitation_source = SEBPrecipitationSource::Microphysics;
            rad.seb_snow_initial_swe_kg_m2 = Real(2.0);
            rad.seb_soil_porosity = Real(0.45);
            rad.seb_q_sfc_default = q_soil0;
            rad.seb_q_deep_default = q_soil0;
        }
        radiation.resize(1);
        radiation.define_level(0, rad, RdoCp, L.ba2d, L.dm, L.ba, L.domain);
        state.define(L.ba, L.dm, RhoQ2_comp + 1, 0);
        state.setVal(Real(0.0));
        state.setVal(Real(1.0), Rho_comp, 1);
        state.setVal(getThgivenRandT(Real(1.0), t_air, RdoCp), RhoTheta_comp, 1);
        qheating.define(L.ba, L.dm, 2, 0);
        rad_fluxes.define(L.ba, L.dm, 4, IntVect(1, 1, 1));   // as ERF allocates it
        t_surf.define(L.ba2d, L.dm, 1, IntVect(1, 1, 0));
        t_surf.setVal(theta_sea);
        lsm.ReSize(1);
        lsm.SetModel<NullSurf>();
    }

    // Warmer air over the building's column (0,3): its downward longwave, and so its
    // residual, differs from the other land columns'.
    void warm_building_column ()
    {
        const Real rho_theta = getThgivenRandT(Real(1.0), t_air + Real(20.0), RdoCp);
        for (MFIter mfi(state); mfi.isValid(); ++mfi) {
            const auto s = state.array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (i == 0 && j == 3) { s(i, j, k, RhoTheta_comp) = rho_theta; }
            });
        }
        Gpu::streamSynchronize();
    }

    void sweep (const iMultiFab* mask)
    {
        const Vector<const MultiFab*> radiation_inputs {};
        radiation.advance(0, 1, Real(0.0), Real(10.0), "pre_dycore", state, nullptr, L.geom, lsm, radiation_inputs,
                          false, &qheating, &rad_fluxes, &t_surf, nullptr, nullptr, nullptr, nullptr, 0.0, false, mask);
    }
    void update (const iMultiFab* mask)
    {
        const Vector<const MultiFab*> radiation_inputs {};
        radiation.advance(0, 1, Real(10.0), Real(10.0), "post_dycore", state, nullptr, L.geom, lsm, radiation_inputs,
                          false, nullptr, nullptr, &t_surf, nullptr, nullptr, nullptr, nullptr, 0.0, false, mask);
    }
    // Upward longwave at the ground [W/m^2] (component 2 of the interface fluxes).
    Real lw_up (int i, int j) const { return value_at(rad_fluxes, IntVect(i, j, 0), 2); }
    Real skin (int i, int j) const
    {
        return value_at(*radiation.prognostic_surface_temperature_state(0), IntVect(i, j, 0), 0);
    }
};

// The temperature the sweep forms from the surface layer's potential temperature, with the
// surface pressure of the lowest cell: the cell-center pressure taken half a cell down.
Real sea_temperature ()
{
    const Real rho = Real(1.0);
    const Real rho_theta = getThgivenRandT(rho, t_air, RdoCp);
    const Real p = erf_surface_temperature::pressure_at_surface(rho, rho_theta, Real(0.0), Real(0.5) * dz);
    Real t = Real(0.0);
    EXPECT_TRUE(erf_surface_temperature::theta_to_temperature(theta_sea, p, RdoCp, t));
    return t;
}

// The number that follows a label in the sweep's printed diagnostics.
Real printed_value (const std::string& report, const std::string& label)
{
    const auto at = report.find(label);
    EXPECT_NE(at, std::string::npos) << label << " is not in:\n" << report;
    if (at == std::string::npos) { return Real(0.0); }
    return static_cast<Real>(std::stod(report.substr(at + label.size())));
}

} // namespace

TEST(TwoStreamLandMask, WaterColumnsTakeTheSeaSurfaceTemperature)
{
    MaskRun r;
    iMultiFab mask(r.L.ba2d, r.L.dm, 1, IntVect(1, 1, 0));
    fill_mask(mask, 0);
    r.sweep(&mask);

    const Real t_sea = sea_temperature();
    ASSERT_GT(std::abs(t_sea - t_skin0), Real(5.0));   // the two boundaries are told apart
    const Real lw_land = stefan_boltzmann * static_cast<Real>(std::pow(t_skin0, 4));
    const Real lw_sea  = stefan_boltzmann * static_cast<Real>(std::pow(t_sea, 4));
    const Real tol = Real(1024) * real_eps;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            const Real expected = (i >= 2) ? lw_sea : lw_land;   // the building's column is land here
            EXPECT_NEAR(r.lw_up(i, j), expected, tol * expected) << i << " " << j;
        }
    }
}

// Half of the sky cloudy: the sweep evaluates each column twice, and both evaluations of
// a water column take the sea surface temperature. With an emissivity of one the upward
// longwave at the ground is sigma T^4 in both, so the blend is too.
TEST(TwoStreamLandMask, CloudyColumnsTakeTheSeaSurfaceTemperature)
{
    MaskRun r(Real(0.5));
    iMultiFab mask(r.L.ba2d, r.L.dm, 1, IntVect(1, 1, 0));
    fill_mask(mask, 0);
    r.sweep(&mask);

    const Real lw_land = stefan_boltzmann * static_cast<Real>(std::pow(t_skin0, 4));
    const Real lw_sea  = stefan_boltzmann * static_cast<Real>(std::pow(sea_temperature(), 4));
    const Real tol = Real(1024) * real_eps;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            const Real expected = (i >= 2) ? lw_sea : lw_land;
            EXPECT_NEAR(r.lw_up(i, j), expected, tol * expected) << i << " " << j;
        }
    }
}

// The printed residual of the balance is over the land columns. Every land column has
// the same sky here, so the mean and the largest are the same number, and they are the
// numbers of the level that is land everywhere. A water column's residual is about
// 154 W/m^2 away (its longwave boundary is 29 K colder), so counting it would move the
// mean by about 77 W/m^2. The diagnostics are printed on the I/O rank: a serial test.
TEST(TwoStreamLandMask, ResidualIsOverTheLandColumns)
{
    Real mean[2] = {Real(0.0), Real(0.0)};
    Real largest[2] = {Real(0.0), Real(0.0)};
    for (int mode = 0; mode < 2; ++mode) {
        MaskRun r(Real(0.0), 1);
        iMultiFab mask(r.L.ba2d, r.L.dm, 1, IntVect(1, 1, 0));
        fill_mask(mask, mode);
        testing::internal::CaptureStdout();
        r.sweep(&mask);
        amrex::OutStream().flush();
        const std::string report = testing::internal::GetCapturedStdout();
        mean[mode] = printed_value(report, "SEB residual (mean) = ");
        largest[mode] = printed_value(report, "SEB residual (max) = ");
    }
    // Sunlit soil: the balance is not closed (about 63 W/m^2), and the residual is far
    // smaller than the 77 W/m^2 a counted water column would add to the mean.
    EXPECT_GT(mean[1], Real(10.0));
    EXPECT_LT(mean[1], Real(77.0));
    // The same columns' numbers, summed over 8 columns and over 16: equal to rounding.
    EXPECT_NEAR(mean[0], mean[1], Real(1.0e-6) * mean[1]);
    EXPECT_NEAR(largest[0], largest[1], Real(1.0e-6) * largest[1]);
    EXPECT_NEAR(mean[0], largest[0], Real(1.0e-4) * largest[0]);
}

// A building's column (mask 2) is a land column in the printed residual. The air over it
// is 20 K warmer, so its residual is not the other land columns': leaving it out would
// change the mean and the largest. The numbers are those of the run in which that column
// carries the mask of ordinary land.
TEST(TwoStreamLandMask, BuildingColumnIsInTheResidual)
{
    Real mean[2] = {Real(0.0), Real(0.0)};
    Real largest[2] = {Real(0.0), Real(0.0)};
    for (int n = 0; n < 2; ++n) {
        MaskRun r(Real(0.0), 1);
        r.warm_building_column();
        iMultiFab mask(r.L.ba2d, r.L.dm, 1, IntVect(1, 1, 0));
        fill_mask(mask, (n == 0) ? 0 : 3);
        testing::internal::CaptureStdout();
        r.sweep(&mask);
        amrex::OutStream().flush();
        const std::string report = testing::internal::GetCapturedStdout();
        mean[n] = printed_value(report, "SEB residual (mean) = ");
        largest[n] = printed_value(report, "SEB residual (max) = ");
    }
    EXPECT_EQ(mean[0], mean[1]);
    EXPECT_EQ(largest[0], largest[1]);
    // The warm column is the largest residual and moves the mean of the 8 land columns.
    EXPECT_GT(largest[1] - mean[1], Real(5.0));
}

TEST(TwoStreamLandMask, WaterColumnsKeepTheirState)
{
    MaskRun r;
    iMultiFab mask(r.L.ba2d, r.L.dm, 1, IntVect(1, 1, 0));
    fill_mask(mask, 0);
    r.sweep(&mask);
    r.update(&mask);
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            if (i >= 2) {
                EXPECT_EQ(r.skin(i, j), t_skin0) << i << " " << j;        // held
            } else {
                EXPECT_GT(r.skin(i, j), t_skin0 + Real(1.0e-3)) << i << " " << j;   // the sunlit soil warms
            }
        }
    }
    // The land columns have one answer: the sky is the same over all of them.
    EXPECT_EQ(r.skin(0, 0), r.skin(1, 2));
    EXPECT_EQ(r.skin(0, 0), r.skin(0, 3));
}

// No mask, and a mask that is land everywhere, are the same run to the last bit; so is the
// land half of the half-and-half level.
TEST(TwoStreamLandMask, LandIsUnchanged)
{
    MaskRun none, land, half;
    iMultiFab all_land(land.L.ba2d, land.L.dm, 1, IntVect(1, 1, 0));
    fill_mask(all_land, 1);
    iMultiFab shore(half.L.ba2d, half.L.dm, 1, IntVect(1, 1, 0));
    fill_mask(shore, 0);
    none.sweep(nullptr);   none.update(nullptr);
    land.sweep(&all_land); land.update(&all_land);
    half.sweep(&shore);    half.update(&shore);
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            EXPECT_EQ(land.skin(i, j), none.skin(i, j)) << i << " " << j;
            EXPECT_EQ(land.lw_up(i, j), none.lw_up(i, j)) << i << " " << j;
            if (i < 2) { EXPECT_EQ(half.skin(i, j), none.skin(i, j)) << i << " " << j; }
        }
    }
    for (int k = 0; k < 4; ++k) {
        EXPECT_EQ(value_at(land.qheating, IntVect(1, 1, k), 1), value_at(none.qheating, IntVect(1, 1, k), 1)) << k;
    }
}

// A level that is water everywhere: nothing evolves, and the update does not fail on a
// balance with no column.
TEST(TwoStreamLandMask, AllWaterLevel)
{
    MaskRun r;
    iMultiFab water(r.L.ba2d, r.L.dm, 1, IntVect(1, 1, 0));
    fill_mask(water, 2);
    r.sweep(&water);
    r.update(&water);
    const MultiFab& t = *r.radiation.prognostic_surface_temperature_state(0);
    EXPECT_EQ(t.min(0), t_skin0);
    EXPECT_EQ(t.max(0), t_skin0);
    const Real lw_sea = stefan_boltzmann * static_cast<Real>(std::pow(sea_temperature(), 4));
    EXPECT_NEAR(r.lw_up(1, 1), lw_sea, Real(1024) * real_eps * lw_sea);
}

// The fields of AverageDownIsOverTheLandColumns. Outside the TEST body: a device lambda
// cannot be launched from a private member.
void fill_average_down_case (MultiFab& fine, iMultiFab& fmask, iMultiFab& cmask)
{
    for (MFIter mfi(fine); mfi.isValid(); ++mfi) {
        const auto f = fine.array(mfi);
        const auto m = fmask.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            f(i, j, k, 0) = static_cast<Real>(10 * i + j);
            f(i, j, k, 1) = -static_cast<Real>(10 * i + j);
            // Water: the fine columns i = 3, 4, 5 with j < 2. A building's column is land.
            m(i, j, k) = (j < 2 && i >= 3 && i <= 5) ? 0 : ((i == 0 && j == 0) ? 2 : 1);
        });
    }
    for (MFIter mfi(cmask); mfi.isValid(); ++mfi) {
        const auto m = cmask.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            m(i, j, k) = (i == 3 && j == 3) ? 0 : 1;
        });
    }
    Gpu::streamSynchronize();
}

// The average-down of the balance's fields between two levels with water
// (ERF_SEBLandAverageDown.H). An 8 x 8 fine level over all of a 4 x 4 coarse level, ratio 2,
// fine value 10 i + j in component 0 and its negative in component 1.
//   - coarse (0,0), land over four fine land columns: the mean of the four, 5.5;
//   - coarse (1,0), land over fine land (2,0), (2,1) and fine water (3,0), (3,1): 20.5, the
//     mean of the two land columns, not 25.5;
//   - coarse (2,0), land over four fine water columns: kept;
//   - coarse (3,3), water over four fine land columns: kept.
TEST(TwoStreamLandMask, AverageDownIsOverTheLandColumns)
{
    const Box cdomain(IntVect(0, 0, 0), IntVect(3, 3, 0));
    const Box fdomain(IntVect(0, 0, 0), IntVect(7, 7, 0));
    BoxArray cba(cdomain), fba(fdomain);
    cba.maxSize(IntVect(2, 4, 1));
    fba.maxSize(IntVect(8, 4, 1));
    const DistributionMapping cdm(cba), fdm(fba);
    MultiFab crse(cba, cdm, 2, 0), fine(fba, fdm, 2, 0);
    iMultiFab cmask(cba, cdm, 1, 0), fmask(fba, fdm, 1, 0);
    crse.setVal(Real(-7.0), 0, 1);
    crse.setVal(Real(7.0), 1, 1);
    fill_average_down_case(fine, fmask, cmask);

    seb_average_down_land(fine, crse, fmask, cmask, IntVect(2, 2, 1));

    const Real tol = Real(64) * real_eps;
    EXPECT_NEAR(value_at(crse, IntVect(0, 0, 0), 0), Real(5.5), tol);
    EXPECT_NEAR(value_at(crse, IntVect(0, 0, 0), 1), Real(-5.5), tol);
    EXPECT_NEAR(value_at(crse, IntVect(1, 0, 0), 0), Real(20.5), Real(20.5) * tol);
    EXPECT_NEAR(value_at(crse, IntVect(1, 0, 0), 1), Real(-20.5), Real(20.5) * tol);
    EXPECT_EQ(value_at(crse, IntVect(2, 0, 0), 0), Real(-7.0));
    EXPECT_EQ(value_at(crse, IntVect(2, 0, 0), 1), Real(7.0));
    EXPECT_EQ(value_at(crse, IntVect(3, 3, 0), 0), Real(-7.0));
    EXPECT_EQ(value_at(crse, IntVect(3, 3, 0), 1), Real(7.0));
    // An all-land coarse column away from the water: fine (4..5, 4..5), mean 49.5.
    EXPECT_NEAR(value_at(crse, IntVect(2, 2, 0), 0), Real(49.5), Real(49.5) * tol);
}

// Every field goes down in one call (ERF_SEBLandAverageDown.H): the result of each is the
// result of a call of its own, a pair with a null pointer is left out, and the fields may
// differ in their number of components. The second field has one component, 3 (10 i + j):
// coarse (1,0) then holds 3 x 20.5 = 61.5.
TEST(TwoStreamLandMask, AverageDownOfSeveralFieldsIsOneCall)
{
    const Box cdomain(IntVect(0, 0, 0), IntVect(3, 3, 0));
    const Box fdomain(IntVect(0, 0, 0), IntVect(7, 7, 0));
    BoxArray cba(cdomain), fba(fdomain);
    cba.maxSize(IntVect(2, 4, 1));
    fba.maxSize(IntVect(8, 4, 1));
    const DistributionMapping cdm(cba), fdm(fba);
    MultiFab crse_a(cba, cdm, 2, 0), fine_a(fba, fdm, 2, 0);
    MultiFab crse_b(cba, cdm, 1, 0), fine_b(fba, fdm, 1, 0);
    MultiFab one_a(cba, cdm, 2, 0), one_b(cba, cdm, 1, 0);
    iMultiFab cmask(cba, cdm, 1, 0), fmask(fba, fdm, 1, 0);
    fill_average_down_case(fine_a, fmask, cmask);
    MultiFab::Copy(fine_b, fine_a, 0, 0, 1, 0);
    fine_b.mult(Real(3.0));
    for (MultiFab* c : {&crse_a, &one_a}) { c->setVal(Real(-7.0), 0, 1); c->setVal(Real(7.0), 1, 1); }
    crse_b.setVal(Real(11.0));
    one_b.setVal(Real(11.0));

    seb_average_down_land(Vector<const MultiFab*>{&fine_a, nullptr, &fine_b},
                          Vector<MultiFab*>{&crse_a, &crse_b, &crse_b}, fmask, cmask, IntVect(2, 2, 1));
    seb_average_down_land(fine_a, one_a, fmask, cmask, IntVect(2, 2, 1));
    seb_average_down_land(fine_b, one_b, fmask, cmask, IntVect(2, 2, 1));

    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            const IntVect p(i, j, 0);
            EXPECT_EQ(value_at(crse_a, p, 0), value_at(one_a, p, 0)) << i << " " << j;
            EXPECT_EQ(value_at(crse_a, p, 1), value_at(one_a, p, 1)) << i << " " << j;
            EXPECT_EQ(value_at(crse_b, p, 0), value_at(one_b, p, 0)) << i << " " << j;
        }
    }
    const Real tol = Real(64) * real_eps;
    EXPECT_NEAR(value_at(crse_a, IntVect(1, 0, 0), 0), Real(20.5), Real(20.5) * tol);
    EXPECT_NEAR(value_at(crse_b, IntVect(1, 0, 0), 0), Real(61.5), Real(61.5) * tol);
    EXPECT_EQ(value_at(crse_b, IntVect(2, 0, 0), 0), Real(11.0));
    EXPECT_EQ(value_at(crse_b, IntVect(3, 3, 0), 0), Real(11.0));
}

// Whether a level has water is read from its mask, with no radiation step before it, and
// again after the mask is rewritten (the WW3 coupling rewrites it during a run). A
// building's column (2) is land; no mask is all land.
void set_mask_column (iMultiFab& mask, int ic, int jc, int value)
{
    for (MFIter mfi(mask); mfi.isValid(); ++mfi) {
        const auto m = mask.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (i == ic && j == jc) { m(i, j, k) = value; }
        });
    }
    Gpu::streamSynchronize();
}

int count_water (const iMultiFab* mask)
{
    int n = 0;
    for (MFIter mfi(*mask); mfi.isValid(); ++mfi) {
        const SEBLandMaskView v = seb_land_mask_view(mask, mfi);
        const SEBLandMaskView off = seb_land_mask_view(nullptr, mfi);
        const Box& bx = mfi.validbox();
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int) -> GpuTuple<int> {
            return {(v.water(i, j) ? 1 : 0) + (off.water(i, j) ? 1000 : 0)};
        });
        n += amrex::get<0>(reduce_data.value());
    }
    ParallelDescriptor::ReduceIntSum(n);
    return n;
}

TEST(TwoStreamLandMask, WaterIsReadFromTheMask)
{
    const Box domain(IntVect(0, 0, 0), IntVect(3, 3, 0));
    BoxArray ba(domain);
    ba.maxSize(IntVect(2, 2, 1));
    const DistributionMapping dm(ba);
    iMultiFab mask(ba, dm, 1, IntVect(1, 1, 0));
    mask.setVal(1);
    EXPECT_FALSE(seb_land_mask_has_water(nullptr));
    EXPECT_FALSE(seb_land_mask_has_water(&mask));
    set_mask_column(mask, 1, 2, 2);
    EXPECT_FALSE(seb_land_mask_has_water(&mask));
    EXPECT_EQ(count_water(&mask), 0);
    set_mask_column(mask, 3, 0, 0);
    EXPECT_TRUE(seb_land_mask_has_water(&mask));
    EXPECT_EQ(count_water(&mask), 1);
    set_mask_column(mask, 3, 0, 1);
    EXPECT_FALSE(seb_land_mask_has_water(&mask));
}

// The fill of a new fine level from a parent with water (ERF_SEBLandFill.H). An 8 x 8
// coarse level, walls in x and periodic in y, under a 16 x 16 fine level, ratio 2. The
// coarse land columns hold 100 + 3 i + 7 j in component 0 and its negative in component 1.
// The coarse water columns, i = 5..7 with j = 2..4, hold the "held" values -999 and 999.
// The fine columns (10,4) and (12,6) are water; every other fine column is land.
struct FillCase
{
    Box cdomain{IntVect(0, 0, 0), IntVect(7, 7, 0)};
    Box fdomain{IntVect(0, 0, 0), IntVect(15, 15, 0)};
    Geometry cgeom, fgeom;
    BoxArray cba, fba;
    DistributionMapping cdm, fdm;
    MultiFab crse, fine;
    iMultiFab cmask, fmask;
    Vector<BCRec> bcs;

    explicit FillCase (bool with_water)
    {
        const RealBox rb({Real(0.), Real(0.), Real(0.)}, {Real(800.), Real(800.), Real(100.)});
        const Array<int, 3> periodic{0, 1, 0};
        cgeom.define(cdomain, rb, 0, periodic);
        fgeom.define(fdomain, rb, 0, periodic);
        cba.define(cdomain); cba.maxSize(IntVect(4, 4, 1));
        fba.define(fdomain); fba.maxSize(IntVect(8, 4, 1));
        cdm.define(cba); fdm.define(fba);
        const IntVect ng(1, 1, 0);
        crse.define(cba, cdm, 2, ng); fine.define(fba, fdm, 2, ng);
        cmask.define(cba, cdm, 1, ng); fmask.define(fba, fdm, 1, ng);
        // The mask's ghost columns beyond the walls read as land: the fill must leave them out.
        cmask.setVal(1); fmask.setVal(1);
        BCRec bc;
        bc.setLo(0, BCType::foextrap); bc.setHi(0, BCType::foextrap);
        bc.setLo(1, BCType::int_dir);  bc.setHi(1, BCType::int_dir);
        bc.setLo(2, BCType::foextrap); bc.setHi(2, BCType::foextrap);
        bcs.assign(2, bc);
        fill(with_water);
    }

    void fill (bool with_water);
};

void fill_fill_case (MultiFab& crse, iMultiFab& cmask, iMultiFab& fmask, bool with_water)
{
    for (MFIter mfi(crse); mfi.isValid(); ++mfi) {
        const auto c = crse.array(mfi);
        const auto m = cmask.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const bool water = with_water && i >= 5 && j >= 2 && j <= 4;
            m(i, j, k) = water ? 0 : 1;
            const Real v = water ? Real(-999.0) : static_cast<Real>(100 + 3 * i + 7 * j);
            c(i, j, k, 0) = v;
            c(i, j, k, 1) = -v;
        });
    }
    for (MFIter mfi(fmask); mfi.isValid(); ++mfi) {
        const auto m = fmask.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            m(i, j, k) = (with_water && ((i == 10 && j == 4) || (i == 12 && j == 6))) ? 0 : 1;
        });
    }
    Gpu::streamSynchronize();
}

void FillCase::fill (bool with_water)
{
    fill_fill_case(crse, cmask, fmask, with_water);
    cmask.FillBoundary(cgeom.periodicity());
    fmask.FillBoundary(fgeom.periodicity());
    fine.setVal(Real(-1.0));
    seb_extend_halo(crse, cgeom);
}

// The interpolation with slopes of the same coarse field, as on a level with no water.
void fill_with_slopes (FillCase& c, MultiFab& out)
{
    out.define(c.fba, c.fdm, 2, IntVect(1, 1, 0));
    out.setVal(Real(-1.0));
    InterpFromCoarseLevel(out, out.nGrowVect(), IntVect(0, 0, 0), c.crse, 0, 0, 2,
                          c.cgeom, c.fgeom, IntVect(2, 2, 1), &cell_cons_interp, c.bcs, 0);
}

TEST(TwoStreamLandMask, FillOfANewLevelTakesLandValues)
{
    FillCase c(true);
    MultiFab slopes;
    fill_with_slopes(c, slopes);
    Long n_none = seb_fill_from_coarse_over_land(c.fine, c.crse, c.cmask, &c.fmask,
                                                 c.cgeom, c.fgeom, IntVect(2, 2, 1), c.bcs);
    ParallelDescriptor::ReduceLongSum(n_none);

    const Real tol = Real(64) * real_eps;
    // Parent (1,4) and its eight neighbors are land: the interpolation with slopes, which
    // is exact for a linear field: 100 + 3 (0.75) + 7 (3.75) = 128.5 at the fine column (2,8).
    EXPECT_EQ(value_at(c.fine, IntVect(2, 8, 0), 0), value_at(slopes, IntVect(2, 8, 0), 0));
    EXPECT_NEAR(value_at(c.fine, IntVect(2, 8, 0), 0), Real(128.5), Real(128.5) * tol);
    EXPECT_NEAR(value_at(c.fine, IntVect(2, 8, 0), 1), Real(-128.5), Real(128.5) * tol);
    // Parent (4,2) is land next to the water column (5,2): its own value 126, no slope.
    EXPECT_EQ(value_at(c.fine, IntVect(8, 4, 0), 0), Real(126.0));
    EXPECT_EQ(value_at(c.fine, IntVect(9, 5, 0), 0), Real(126.0));
    EXPECT_EQ(value_at(c.fine, IntVect(9, 5, 0), 1), Real(-126.0));
    // Parent (5,2) is water. Its land neighbors (4,1), (5,1), (6,1), (4,2), (4,3) hold
    // 119, 122, 125, 126, 133: the mean is 125, not the held -999.
    EXPECT_NEAR(value_at(c.fine, IntVect(11, 4, 0), 0), Real(125.0), Real(125.0) * tol);
    EXPECT_NEAR(value_at(c.fine, IntVect(11, 5, 0), 1), Real(-125.0), Real(125.0) * tol);
    // The fine water column (10,4) over the same parent keeps the parent's held value.
    EXPECT_EQ(value_at(c.fine, IntVect(10, 4, 0), 0), Real(-999.0));
    // Parent (7,2) is water at the wall: the columns beyond the wall are left out. Its land
    // neighbors (6,1), (7,1) hold 125 and 128: the mean is 126.5.
    EXPECT_NEAR(value_at(c.fine, IntVect(15, 4, 0), 0), Real(126.5), Real(126.5) * tol);
    // Parents (6,3) and (7,3) are water with no land neighbor: the held value, and the
    // fine land columns over them are counted: 8 less the fine water column (12,6).
    EXPECT_EQ(value_at(c.fine, IntVect(13, 7, 0), 0), Real(-999.0));
    EXPECT_EQ(value_at(c.fine, IntVect(14, 6, 0), 1), Real(999.0));
    EXPECT_EQ(n_none, Long(7));
}

// With no water on the parent the fill is the interpolation with slopes in every column.
TEST(TwoStreamLandMask, FillWithNoWaterIsTheInterpolationWithSlopes)
{
    FillCase c(false);
    MultiFab slopes;
    fill_with_slopes(c, slopes);
    Long n_none = seb_fill_from_coarse_over_land(c.fine, c.crse, c.cmask, nullptr,
                                                 c.cgeom, c.fgeom, IntVect(2, 2, 1), c.bcs);
    ParallelDescriptor::ReduceLongSum(n_none);
    EXPECT_EQ(n_none, Long(0));
    MultiFab::Subtract(slopes, c.fine, 0, 0, 2, 0);
    EXPECT_EQ(slopes.norminf(0, 0), Real(0.0));
    EXPECT_EQ(slopes.norminf(1, 0), Real(0.0));
    EXPECT_NEAR(value_at(c.fine, IntVect(2, 8, 0), 0), Real(128.5), Real(128.5) * Real(64) * real_eps);
}

// The radiation step hands the mask on to the snow albedo, the snow step and the hydrology
// step. 2 kg/m^2 of snow lies everywhere at the start, and 1 kg/m^2 of rain and 0.5 kg/m^2
// of snow fall in the 10 s step.
//   - albedo: the net shortwave of a water column is that of a run with no snow at all
//     (albedo 0.2); on land the snow raises the albedo to 0.2 + 0.5 (2 / 12) = 0.283;
//   - snow step: after the step no snow lies on water, and snow lies on land;
//   - hydrology: the soil water of a water column is the starting value; on land the rain
//     adds 1 / (rho_w d_s) = 0.01 to it.
TEST(TwoStreamLandMask, SnowAndRainSkipTheWaterColumns)
{
    MaskRun r(Real(0.0), 0, true), bare;
    iMultiFab mask(r.L.ba2d, r.L.dm, 1, IntVect(1, 1, 0));
    iMultiFab bare_mask(bare.L.ba2d, bare.L.dm, 1, IntVect(1, 1, 0));
    fill_mask(mask, 0);
    fill_mask(bare_mask, 0);

    MultiFab rain(r.L.ba, r.L.dm, 1, 0), snowfall(r.L.ba, r.L.dm, 1, 0);
    SurfacePrecipAccumulationSources src;
    src.rain.accum = &rain;
    src.snow.accum = &snowfall;
    rain.setVal(Real(0.0));
    snowfall.setVal(Real(0.0));
    r.radiation.hydrology_begin_step(0, src, 0);
    r.sweep(&mask);
    bare.sweep(&bare_mask);

    rain.setVal(Real(1.0));
    snowfall.setVal(Real(0.5));
    r.radiation.hydrology_end_step(0, src, 0, Real(10.0));
    r.update(&mask);
    bare.update(&bare_mask);

    ASSERT_NE(r.radiation.seb_sw_net(0), nullptr);
    ASSERT_NE(bare.radiation.seb_sw_net(0), nullptr);
    ASSERT_NE(r.radiation.snow_swe(0), nullptr);
    ASSERT_NE(r.radiation.seb_q_sfc(0), nullptr);
    // The net shortwave the balance used: from the sweep, with the albedo of its start.
    const auto sw_net = [] (const MaskRun& run, int i, int j) {
        return value_at(*run.radiation.seb_sw_net(0), IntVect(i, j, 0), 0);
    };
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            if (i >= 2) {
                EXPECT_EQ(sw_net(r, i, j), sw_net(bare, i, j)) << i << " " << j;
            } else {
                EXPECT_LT(sw_net(r, i, j), Real(0.95) * sw_net(bare, i, j)) << i << " " << j;
            }
        }
    }
    EXPECT_GT(sw_net(bare, 3, 0), Real(100.0));

    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            const Real swe = value_at(*r.radiation.snow_swe(0), IntVect(i, j, 0), 0);
            const Real q_s = value_at(*r.radiation.seb_q_sfc(0), IntVect(i, j, 0), 0);
            if (i >= 2) {
                EXPECT_EQ(swe, Real(0.0)) << i << " " << j;
                EXPECT_EQ(q_s, q_soil0) << i << " " << j;
                EXPECT_EQ(r.skin(i, j), t_skin0) << i << " " << j;
            } else {
                EXPECT_GT(swe, Real(0.1)) << i << " " << j;
                EXPECT_GT(q_s, q_soil0 + Real(0.005)) << i << " " << j;
            }
        }
    }
}

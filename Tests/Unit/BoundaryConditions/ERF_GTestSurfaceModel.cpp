#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_RealBox.H>

#include <gtest/gtest.h>

#include "ERF_SurfaceModel.H"

namespace {

using amrex::Real;

struct SurfaceModelFixture {
    amrex::Box domain{amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
                      amrex::IntVect(AMREX_D_DECL(1, 1, 3))};
    amrex::BoxArray ba{domain};
    amrex::DistributionMapping dm{ba};
    amrex::Geometry geom;
    amrex::Vector<amrex::BoxArray> bas{ba};
    amrex::Vector<amrex::Geometry> geoms;
    amrex::Vector<amrex::DistributionMapping> dms{dm};
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>> masks;
    SolverChoice solver_choice{};
    std::unique_ptr<SurfaceModel> model;

    explicit SurfaceModelFixture(const int nlevels = 1)
    {
        amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                {AMREX_D_DECL(2.0, 2.0, 4.0)});
        std::array<int, AMREX_SPACEDIM> periodic{};
        geom = amrex::Geometry(domain, &real_box, 0, periodic.data());

        bas.clear();
        dms.clear();
        geoms.clear();
        masks.resize(nlevels);
        for (int lev = 0; lev < nlevels; ++lev) {
            const amrex::Box level_domain(
                amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
                amrex::IntVect(AMREX_D_DECL((1 << (lev + 1)) - 1,
                                            (1 << (lev + 1)) - 1,
                                            3)));
            const amrex::BoxArray level_ba(level_domain);
            const amrex::DistributionMapping level_dm(level_ba);

            bas.push_back(level_ba);
            dms.push_back(level_dm);
            if (lev == 0) {
                geoms.push_back(geom);
            } else {
                geoms.emplace_back(level_domain, &real_box, 0, periodic.data());
            }

            amrex::BoxList surface_boxes = level_ba.boxList();
            for (auto& box : surface_boxes) {
                box.setRange(2, 0);
            }
            amrex::BoxArray surface_ba(std::move(surface_boxes));
            masks[lev].emplace_back(std::make_unique<amrex::iMultiFab>(
                surface_ba, level_dm, 1, amrex::IntVect(0)));
            masks[lev][0]->setVal(1);
        }

        model = std::make_unique<SurfaceModel>(nlevels, bas, geoms, dms,
                                               solver_choice, masks);
        const amrex::BCRec interior_bc(
            AMREX_D_DECL(amrex::BCType::int_dir, amrex::BCType::int_dir,
                         amrex::BCType::int_dir),
            AMREX_D_DECL(amrex::BCType::int_dir, amrex::BCType::int_dir,
                         amrex::BCType::int_dir));
        const amrex::Vector<amrex::BCRec> domain_bcs_type(
            AMREX_SPACEDIM + NBCVAR_max, interior_bc);
        const amrex::Vector<amrex::IntVect> ref_ratios(
            nlevels > 1 ? nlevels - 1 : 1,
            amrex::IntVect(AMREX_D_DECL(2, 2, 1)));
        for (int lev = 0; lev < nlevels; ++lev) {
            model->initialize_for_level(lev, bas[lev], geoms[lev], dms[lev],
                                        masks[lev], domain_bcs_type, ref_ratios);
        }
    }

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> make_fields(
        const Real first, const Real increment, const int ngz)
    {
        return make_fields(0, first, increment, ngz);
    }

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> make_fields(
        const int lev, const Real first, const Real increment, const int ngz)
    {
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> fields(6);
        for (int i = 0; i < 6; ++i) {
            fields[i] = std::make_unique<amrex::MultiFab>(
                bas[lev], dms[lev], 1, amrex::IntVect(1, 1, ngz));
            fields[i]->setVal(first + increment * i);
        }
        return fields;
    }

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> make_surface_fields(
        const Real first, const Real increment)
    {
        amrex::BoxList surface_boxes = ba.boxList();
        for (auto& box : surface_boxes) { box.setRange(2, 0); }
        const amrex::BoxArray surface_ba(std::move(surface_boxes));
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> fields(6);
        for (int i = 0; i < 6; ++i) {
            fields[i] = std::make_unique<amrex::MultiFab>(
                surface_ba, dm, 1, amrex::IntVect(1, 1, 0));
            fields[i]->setVal(first + increment * i);
        }
        return fields;
    }

    static amrex::Vector<amrex::MultiFab*> pointers(
        amrex::Vector<std::unique_ptr<amrex::MultiFab>>& fields)
    {
        amrex::Vector<amrex::MultiFab*> result;
        for (auto& field : fields) result.push_back(field.get());
        return result;
    }

    void configure_models(amrex::Vector<std::unique_ptr<amrex::MultiFab>>& land,
                          amrex::Vector<std::unique_ptr<amrex::MultiFab>>& urban)
    {
        configure_models(0, land, urban);
    }

    void configure_models(const int lev,
                          amrex::Vector<std::unique_ptr<amrex::MultiFab>>& land,
                          amrex::Vector<std::unique_ptr<amrex::MultiFab>>& urban)
    {
        const amrex::Vector<int> field_indices{0, 1, 2, 3, 4, 5};
        model->set_model_data(lev, pointers(land), {"f0", "f1", "f2", "f3", "f4", "f5"},
                              SurfaceModelType::LAND);
        model->set_model_data(lev, pointers(urban), {"f0", "f1", "f2", "f3", "f4", "f5"},
                              SurfaceModelType::URBAN);
        model->set_model_fields(SurfaceModelType::LAND, field_indices);
        model->set_model_fields(SurfaceModelType::URBAN, field_indices);
    }
};

TEST(SurfaceModel, InitializesOutputsAndDefaultWeights)
{
    SurfaceModelFixture fixture;

    EXPECT_EQ(fixture.model->get_ustar(0)->nComp(), 2);
    EXPECT_EQ(fixture.model->get_tstar(0)->nComp(), 1);
    EXPECT_EQ(fixture.model->get_qstar(0)->nComp(), 1);
    EXPECT_EQ(fixture.model->get_tsurf(0)->max(0), Real(300.0));
    EXPECT_EQ(fixture.model->get_ustar(0)->max(0), Real(0.0));
    EXPECT_EQ(fixture.model->get_wavg_factors(0)->max(SurfaceModelType::LAND), Real(1.0));
    EXPECT_EQ(fixture.model->get_wavg_factors(0)->max(SurfaceModelType::URBAN), Real(0.0));
}

TEST(SurfaceModel, CalculatesWeightedFluxOutputs)
{
    SurfaceModelFixture fixture;
    auto land = fixture.make_fields(Real(10.0), Real(10.0), 1);
    auto urban = fixture.make_fields(Real(100.0), Real(10.0), 0);
    fixture.configure_models(land, urban);

    amrex::MultiFab urban_fraction(fixture.ba, fixture.dm, 1, amrex::IntVect(0));
    urban_fraction.setVal(Real(0.25));
    fixture.model->calculate_weight_average(0, &urban_fraction);

    EXPECT_NEAR(fixture.model->get_wavg_factors(0)->max(SurfaceModelType::LAND), Real(0.75), 1.e-12);
    EXPECT_NEAR(fixture.model->get_wavg_factors(0)->max(SurfaceModelType::URBAN), Real(0.25), 1.e-12);
    EXPECT_NEAR(fixture.model->get_ustar(0)->max(0), Real(32.5), 1.e-12);
    EXPECT_NEAR(fixture.model->get_ustar(0)->max(1), Real(42.5), 1.e-12);
    EXPECT_NEAR(fixture.model->get_tstar(0)->max(0), Real(52.5), 1.e-12);
    EXPECT_NEAR(fixture.model->get_qstar(0)->max(0), Real(62.5), 1.e-12);
    EXPECT_NEAR(fixture.model->get_tsurf(0)->max(0), Real(72.5), 1.e-12);
}

TEST(SurfaceModel, RegistersAndAveragesModelsOnMultipleLevels)
{
    SurfaceModelFixture fixture(2);

    auto land_coarse = fixture.make_fields(0, Real(10.0), Real(10.0), 1);
    auto urban_coarse = fixture.make_fields(0, Real(100.0), Real(10.0), 0);
    fixture.configure_models(0, land_coarse, urban_coarse);

    auto land_fine = fixture.make_fields(1, Real(20.0), Real(10.0), 1);
    auto urban_fine = fixture.make_fields(1, Real(200.0), Real(10.0), 0);
    fixture.configure_models(1, land_fine, urban_fine);

    amrex::MultiFab coarse_fraction(fixture.bas[0], fixture.dms[0], 1,
                                    amrex::IntVect(0));
    amrex::MultiFab fine_fraction(fixture.bas[1], fixture.dms[1], 1,
                                  amrex::IntVect(0));
    coarse_fraction.setVal(Real(0.25));
    fine_fraction.setVal(Real(0.5));

    fixture.model->calculate_weight_average(0, &coarse_fraction);
    fixture.model->calculate_weight_average(1, &fine_fraction);

    EXPECT_EQ(fixture.model->get_ustar(1)->boxArray().minimalBox().length(0), 4);
    EXPECT_NEAR(fixture.model->get_wavg_factors(0)->max(SurfaceModelType::LAND),
                Real(0.75), 1.e-12);
    EXPECT_NEAR(fixture.model->get_wavg_factors(1)->max(SurfaceModelType::LAND),
                Real(0.5), 1.e-12);
    EXPECT_NEAR(fixture.model->get_wavg_factors(1)->max(SurfaceModelType::URBAN),
                Real(0.5), 1.e-12);
    EXPECT_NEAR(fixture.model->get_ustar(1)->max(0), Real(110.0), 1.e-12);
    EXPECT_NEAR(fixture.model->get_ustar(1)->max(1), Real(120.0), 1.e-12);
    EXPECT_NEAR(fixture.model->get_tstar(1)->max(0), Real(130.0), 1.e-12);
    EXPECT_NEAR(fixture.model->get_qstar(1)->max(0), Real(140.0), 1.e-12);
    EXPECT_NEAR(fixture.model->get_tsurf(1)->max(0), Real(150.0), 1.e-12);
}

TEST(SurfaceModel, UrbanOnlyModelUsesFullSurface)
{
    SurfaceModelFixture fixture;
    auto urban = fixture.make_fields(Real(100.0), Real(10.0), 0);
    const amrex::Vector<int> field_indices{0, 1, 2, 3, 4, 5};
    fixture.model->set_model_data(
        0, fixture.pointers(urban), {"f0", "f1", "f2", "f3", "f4", "f5"},
        SurfaceModelType::URBAN);
    fixture.model->set_model_fields(SurfaceModelType::URBAN, field_indices);
    fixture.model->calculate_weight_average(0, nullptr);

    EXPECT_EQ(fixture.model->get_wavg_factors(0)->max(SurfaceModelType::LAND), Real(0.0));
    EXPECT_EQ(fixture.model->get_wavg_factors(0)->max(SurfaceModelType::URBAN), Real(1.0));
    EXPECT_NEAR(fixture.model->get_ustar(0)->max(0), Real(100.0), 1.e-12);
    EXPECT_NEAR(fixture.model->get_ustar(0)->max(1), Real(110.0), 1.e-12);
    EXPECT_NEAR(fixture.model->get_tstar(0)->max(0), Real(120.0), 1.e-12);
    EXPECT_NEAR(fixture.model->get_qstar(0)->max(0), Real(130.0), 1.e-12);
    EXPECT_NEAR(fixture.model->get_tsurf(0)->max(0), Real(140.0), 1.e-12);
}

TEST(SurfaceModel, UnconfiguredProviderFieldsAreSkipped)
{
    SurfaceModelFixture fixture;
    auto urban = fixture.make_fields(Real(100.0), Real(10.0), 0);
    fixture.model->set_model_data(
        0, fixture.pointers(urban), {"f0", "f1", "f2", "f3", "f4", "f5"},
        SurfaceModelType::URBAN);

    fixture.model->calculate_weight_average(0, nullptr);

    EXPECT_EQ(fixture.model->get_ustar(0)->max(0), Real(0.0));
    EXPECT_EQ(fixture.model->get_ustar(0)->max(1), Real(0.0));
    EXPECT_EQ(fixture.model->get_tstar(0)->max(0), Real(0.0));
    EXPECT_EQ(fixture.model->get_qstar(0)->max(0), Real(0.0));
    EXPECT_EQ(fixture.model->get_tsurf(0)->max(0), Real(0.0));
}

TEST(SurfaceModel, SingleEnabledModelUsesFullSurfaceWithoutFraction)
{
    SurfaceModelFixture fixture;
    auto land = fixture.make_fields(Real(10.0), Real(1.0), 1);
    auto urban = fixture.make_fields(Real(100.0), Real(1.0), 0);
    fixture.model->set_model_data(0, fixture.pointers(land), {"f0", "f1", "f2", "f3", "f4", "f5"},
                                  SurfaceModelType::LAND);
    fixture.model->set_model_fields(SurfaceModelType::LAND, {0, 1, 2, 3, 4, 5});
    fixture.model->calculate_weight_average(0, nullptr);

    EXPECT_EQ(fixture.model->get_wavg_factors(0)->max(SurfaceModelType::LAND), Real(1.0));
    EXPECT_EQ(fixture.model->get_wavg_factors(0)->max(SurfaceModelType::URBAN), Real(0.0));
    EXPECT_EQ(fixture.model->get_tsurf(0)->max(0), Real(14.0));
}

TEST(SurfaceModel, ApplyWeightAverageScalesBothModels)
{
    SurfaceModelFixture fixture;
    auto land = fixture.make_fields(Real(1.0), Real(1.0), 1);
    auto urban = fixture.make_fields(Real(1.0), Real(1.0), 0);
    fixture.configure_models(land, urban);
    amrex::MultiFab urban_fraction(fixture.ba, fixture.dm, 1, amrex::IntVect(0));
    urban_fraction.setVal(Real(0.25));
    fixture.model->calculate_weight_average(0, &urban_fraction);

    amrex::MultiFab land_data(fixture.ba, fixture.dm, 2, amrex::IntVect(0));
    amrex::MultiFab urban_data(fixture.ba, fixture.dm, 2, amrex::IntVect(0));
    amrex::MultiFab land_weighted(fixture.ba, fixture.dm, 2, amrex::IntVect(0));
    amrex::MultiFab urban_weighted(fixture.ba, fixture.dm, 2, amrex::IntVect(0));
    land_data.setVal(Real(4.0));
    urban_data.setVal(Real(8.0));
    fixture.model->apply_weight_average(0, &land_data, &land_weighted,
                                        &urban_data, &urban_weighted);

    EXPECT_NEAR(land_data.max(0), Real(4.0), 1.e-12);
    EXPECT_NEAR(urban_data.max(0), Real(8.0), 1.e-12);
    EXPECT_NEAR(land_weighted.max(0), Real(3.0), 1.e-12);
    EXPECT_NEAR(urban_weighted.max(0), Real(2.0), 1.e-12);
}

TEST(SurfaceModel, ExtraFieldsAreWeightedWithoutModifyingProviders)
{
    SurfaceModelFixture fixture;
    auto land = fixture.make_fields(Real(10.0), Real(10.0), 1);
    auto urban = fixture.make_fields(Real(100.0), Real(10.0), 0);
    fixture.configure_models(land, urban);

    amrex::MultiFab urban_fraction(fixture.ba, fixture.dm, 1, amrex::IntVect(0));
    urban_fraction.setVal(Real(0.25));
    fixture.model->calculate_weight_average(0, &urban_fraction);
    fixture.model->calculate_weight_average(0, &urban_fraction);

    EXPECT_NEAR(land[5]->max(0), Real(60.0), 1.e-12);
    EXPECT_NEAR(urban[5]->max(0), Real(150.0), 1.e-12);
    ASSERT_NE(fixture.model->get_weighted_model_data(0, SurfaceModelType::LAND, 5), nullptr);
    ASSERT_NE(fixture.model->get_weighted_model_data(0, SurfaceModelType::URBAN, 5), nullptr);
    EXPECT_NEAR(fixture.model->get_weighted_model_data(0, SurfaceModelType::LAND, 5)->max(0),
                Real(45.0), 1.e-12);
    EXPECT_NEAR(fixture.model->get_weighted_model_data(0, SurfaceModelType::URBAN, 5)->max(0),
                Real(37.5), 1.e-12);
}

TEST(SurfaceModel, PointerMappedFieldIsNotScaledTwice)
{
    SurfaceModelFixture fixture;
    auto land = fixture.make_fields(Real(10.0), Real(10.0), 1);
    auto urban = fixture.make_fields(Real(100.0), Real(10.0), 0);
    fixture.configure_models(land, urban);

    amrex::Vector<amrex::MultiFab*> land_ptrs{land[5].get()};
    amrex::Vector<amrex::MultiFab*> urban_ptrs{urban[5].get()};
    fixture.model->register_field_map("mapped", land_ptrs, urban_ptrs);

    amrex::MultiFab urban_fraction(fixture.ba, fixture.dm, 1, amrex::IntVect(0));
    urban_fraction.setVal(Real(0.25));
    fixture.model->calculate_weight_average(0, &urban_fraction);

    EXPECT_NEAR(fixture.model->get_field("mapped")->max(0), Real(82.5), 1.e-12);
    EXPECT_EQ(land[5]->max(0), Real(60.0));
    EXPECT_EQ(urban[5]->max(0), Real(150.0));
}

TEST(SurfaceModel, RegridRecreatesMappedFieldsAndRadiationCache)
{
    SurfaceModelFixture fixture;
    auto land = fixture.make_fields(Real(10.0), Real(10.0), 1);
    fixture.model->set_model_data(
        0, fixture.pointers(land), {"f0", "f1", "f2", "f3", "f4", "f5"},
        SurfaceModelType::LAND);
    fixture.model->set_model_fields(SurfaceModelType::LAND, {0, 1, 2, 3, 4, 5});
    fixture.model->register_radiation_input("tskin", {4, -1});
    fixture.model->register_radiation_output("lw_flux_dn", {0, -1});

    amrex::MultiFab urban_fraction(fixture.ba, fixture.dm, 1, amrex::IntVect(0));
    urban_fraction.setVal(Real(0.0));
    fixture.model->calculate_weight_average(0, &urban_fraction);
    const auto radiation_before = fixture.model->get_radiation_fields(0);
    ASSERT_EQ(radiation_before[0], land[4].get());
    const auto radiation_outputs_before = fixture.model->get_radiation_output_fields(0);
    ASSERT_EQ(radiation_outputs_before[6], land[0].get());

    amrex::BoxList remade_boxes;
    remade_boxes.push_back(amrex::Box(
        amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
        amrex::IntVect(AMREX_D_DECL(0, 1, 3))));
    remade_boxes.push_back(amrex::Box(
        amrex::IntVect(AMREX_D_DECL(1, 0, 0)),
        amrex::IntVect(AMREX_D_DECL(1, 1, 3))));
    amrex::BoxArray remade_ba(std::move(remade_boxes));
    amrex::DistributionMapping remade_dm(remade_ba);

    amrex::BoxList remade_surface_boxes = remade_ba.boxList();
    for (auto& box : remade_surface_boxes) {
        box.setRange(2, 0);
    }
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> remade_masks;
    remade_masks.emplace_back(std::make_unique<amrex::iMultiFab>(
        amrex::BoxArray(std::move(remade_surface_boxes)), remade_dm, 1,
        amrex::IntVect(0)));
    remade_masks[0]->setVal(1);

    fixture.model->initialize_for_level(0, remade_ba, fixture.geom, remade_dm,
                                        remade_masks, {}, {amrex::IntVect(2)});

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> remade_land(6);
    for (int field = 0; field < 6; ++field) {
        remade_land[field] = std::make_unique<amrex::MultiFab>(
            remade_ba, remade_dm, 1, amrex::IntVect(1, 1, 1));
        remade_land[field]->setVal(Real(10.0) + Real(10.0) * field);
    }
    fixture.model->set_model_data(
        0, fixture.pointers(remade_land), {"f0", "f1", "f2", "f3", "f4", "f5"},
        SurfaceModelType::LAND);

    amrex::MultiFab remade_fraction(remade_ba, remade_dm, 1, amrex::IntVect(0));
    remade_fraction.setVal(Real(0.0));
    fixture.model->calculate_weight_average(0, &remade_fraction);

    EXPECT_EQ(remade_land[4]->boxArray(), remade_ba);
    EXPECT_NEAR(remade_land[4]->max(0), Real(50.0), 1.e-12);

    const auto radiation_after = fixture.model->get_radiation_fields(0);
    ASSERT_EQ(radiation_after[0], remade_land[4].get());
    const auto radiation_outputs_after = fixture.model->get_radiation_output_fields(0);
    ASSERT_EQ(radiation_outputs_after[6], remade_land[0].get());
}

TEST(SurfaceModel, RadiationFieldListUsesCanonicalMappings)
{
    const std::vector<std::string> radiation_names{
        "tskin", "emiss", "albedo_vis", "albedo_nir",
        "albedo_vis_diff", "albedo_nir_diff"};

    {
        SurfaceModelFixture fixture;
        auto land = fixture.make_fields(Real(10.0), Real(10.0), 1);
        fixture.model->set_model_data(
            0, fixture.pointers(land), {"f0", "f1", "f2", "f3", "f4", "f5"},
            SurfaceModelType::LAND);
        for (int i = 0; i < static_cast<int>(radiation_names.size()); ++i) {
            fixture.model->register_radiation_input(
                radiation_names[i], std::pair<int, int>{i, -1});
        }

        const auto radiation_fields = fixture.model->get_radiation_fields(0);
        ASSERT_EQ(radiation_fields.size(), radiation_names.size());
        for (int i = 0; i < static_cast<int>(radiation_names.size()); ++i) {
            EXPECT_EQ(radiation_fields[i], land[i].get());
        }
    }

    // The common radiation names can map to different provider-specific
    // indices, as they do for the land and urban surface models.
    SurfaceModelFixture fixture;
    auto land = fixture.make_fields(Real(10.0), Real(10.0), 1);
    auto urban = fixture.make_fields(Real(20.0), Real(10.0), 1);
    fixture.configure_models(land, urban);
    const std::vector<std::pair<int, int>> provider_mappings{
        {0, 0}, {1, 1}, {2, -1}, {3, -1}, {4, -1}, {5, -1}};
    for (int i = 0; i < static_cast<int>(radiation_names.size()); ++i) {
        fixture.model->register_radiation_input(radiation_names[i], provider_mappings[i]);
    }

    const auto radiation_fields = fixture.model->get_radiation_fields(0);
    ASSERT_EQ(radiation_fields.size(), radiation_names.size());
    for (int i = 0; i < static_cast<int>(radiation_names.size()); ++i) {
        EXPECT_NE(radiation_fields[i], nullptr);
    }
    EXPECT_NE(radiation_fields[0], land[0].get());
    EXPECT_NE(radiation_fields[0], urban[0].get());

    // Missing registered fields are represented by nullptr so the radiation
    // backend can use its default RRTMGP values.
    SurfaceModelFixture fallback_fixture;
    auto fallback_land = fallback_fixture.make_fields(Real(10.0), Real(10.0), 1);
    fallback_fixture.model->set_model_data(
        0, fallback_fixture.pointers(fallback_land),
        {"f0", "f1", "f2", "f3", "f4", "f5"}, SurfaceModelType::LAND);
    fallback_fixture.model->register_radiation_input(
        radiation_names[0], std::pair<int, int>{0, -1});

    const auto fallback_fields = fallback_fixture.model->get_radiation_fields(0);
    ASSERT_EQ(fallback_fields.size(), radiation_names.size());
    EXPECT_EQ(fallback_fields[0], fallback_land[0].get());
    for (int i = 1; i < static_cast<int>(fallback_fields.size()); ++i) {
        EXPECT_EQ(fallback_fields[i], nullptr);
    }
}

TEST(SurfaceModel, RadiationInputsResolveProviderModesAndWeights)
{
    const std::vector<std::string> names{
        "tskin", "emiss", "albedo_vis", "albedo_nir",
        "albedo_vis_diff", "albedo_nir_diff"};

    SurfaceModelFixture fixture;
    auto land = fixture.make_surface_fields(Real(10.0), Real(1.0));
    auto urban = fixture.make_surface_fields(Real(20.0), Real(1.0));
    fixture.configure_models(land, urban);

    std::unordered_map<std::string, std::pair<int, int>> mappings;
    for (int i = 0; i < static_cast<int>(names.size()); ++i) {
        mappings.emplace(names[i], std::pair<int, int>{i, i});
    }
    fixture.model->register_radiation_inputs(mappings);

    amrex::MultiFab urban_fraction(
        fixture.model->get_wavg_factors(0)->boxArray(), fixture.dm, 1,
        amrex::IntVect(0));
    urban_fraction.setVal(Real(0.25));
    fixture.model->calculate_weight_average(0, &urban_fraction);

    const auto fields = fixture.model->get_radiation_fields(0);
    ASSERT_EQ(fields.size(), names.size());
    for (int i = 0; i < static_cast<int>(fields.size()); ++i) {
        EXPECT_NE(fields[i], land[i].get());
        EXPECT_NE(fields[i], urban[i].get());
        EXPECT_NEAR(fields[i]->max(0), Real(12.5) + Real(i), 1.e-12);
    }

    SurfaceModelFixture land_fixture;
    auto land_only = land_fixture.make_fields(Real(30.0), Real(1.0), 1);
    land_fixture.model->set_model_data(
        0, land_fixture.pointers(land_only), {"f0", "f1", "f2", "f3", "f4", "f5"},
        SurfaceModelType::LAND);
    land_fixture.model->register_radiation_input("tskin", {0, -1});
    EXPECT_EQ(land_fixture.model->get_radiation_fields(0)[0], land_only[0].get());

    SurfaceModelFixture urban_fixture;
    auto urban_only = urban_fixture.make_fields(Real(40.0), Real(1.0), 0);
    urban_fixture.model->set_model_data(
        0, urban_fixture.pointers(urban_only), {"f0", "f1", "f2", "f3", "f4", "f5"},
        SurfaceModelType::URBAN);
    urban_fixture.model->register_radiation_input("tskin", {-1, 0});
    EXPECT_EQ(urban_fixture.model->get_radiation_fields(0)[0], urban_only[0].get());
}

TEST(SurfaceModel, RadiationOutputsResolveAndDistributeCanonicalMappings)
{
    SurfaceModelFixture fixture;
    auto land = fixture.make_surface_fields(Real(10.0), Real(1.0));
    auto urban = fixture.make_fields(Real(20.0), Real(1.0), 0);
    fixture.model->set_model_data(
        0, fixture.pointers(land), {"f0", "f1", "f2", "f3", "f4", "f5"},
        SurfaceModelType::LAND);
    fixture.model->set_model_data(
        0, fixture.pointers(urban), {"f0", "f1", "f2", "f3", "f4", "f5"},
        SurfaceModelType::URBAN);
    land[0]->setVal(7.0);
    urban[0]->setVal(0.0);

    fixture.model->register_radiation_output("cos_zenith_angle", {0, 0});
    fixture.model->register_radiation_output("sw_flux_dn_dir_vis", {1, -1});
    fixture.model->register_radiation_output("lw_flux_dn", {-1, 2});

    const auto outputs = fixture.model->get_radiation_output_fields(0);
    ASSERT_EQ(outputs.size(), 7);
    EXPECT_EQ(outputs[0], land[0].get());
    EXPECT_EQ(outputs[1], nullptr);
    EXPECT_EQ(outputs[2], land[1].get());
    EXPECT_EQ(outputs[3], nullptr);
    EXPECT_EQ(outputs[4], nullptr);
    EXPECT_EQ(outputs[5], nullptr);
    EXPECT_EQ(outputs[6], urban[2].get());

    fixture.model->distribute_radiation_outputs(0);
    EXPECT_EQ(urban[0]->max(0), Real(7.0));
    EXPECT_EQ(urban[0]->min(0), Real(0.0));

    auto replacement_land = fixture.make_surface_fields(Real(30.0), Real(1.0));
    fixture.model->set_model_data(
        0, fixture.pointers(replacement_land),
        {"f0", "f1", "f2", "f3", "f4", "f5"}, SurfaceModelType::LAND);
    const auto replacement_outputs = fixture.model->get_radiation_output_fields(0);
    EXPECT_EQ(replacement_outputs[0], replacement_land[0].get());
}

TEST(SurfaceModel, CheckpointRoundTripPreservesSyntheticState)
{
    const std::filesystem::path checkpoint =
        std::filesystem::temp_directory_path() / "erf_surface_model_unit_checkpoint";
    struct CheckpointCleanup {
        std::filesystem::path path;
        ~CheckpointCleanup() { std::filesystem::remove_all(path); }
    } cleanup{checkpoint};

    std::filesystem::remove_all(checkpoint);
    ASSERT_TRUE(std::filesystem::create_directories(checkpoint / "Level_0"));

    SurfaceModelFixture writer;
    auto writer_land = writer.make_fields(Real(10.0), Real(10.0), 1);
    auto writer_urban = writer.make_fields(Real(100.0), Real(10.0), 0);
    writer.configure_models(writer_land, writer_urban);
    writer.model->register_field_map("mapped", std::pair<int, int>{5, 5});

    amrex::MultiFab writer_fraction(writer.ba, writer.dm, 1, amrex::IntVect(0));
    writer_fraction.setVal(Real(0.25));
    writer.model->calculate_weight_average(0, &writer_fraction);
    writer.model->WriteCheckpoint(checkpoint.string());

    const std::filesystem::path header = checkpoint / "SurfaceModel_Header";
    ASSERT_TRUE(std::filesystem::exists(header));
    std::ifstream header_stream(header);
    std::string header_text((std::istreambuf_iterator<char>(header_stream)),
                            std::istreambuf_iterator<char>());
    EXPECT_NE(header_text.find("Checkpoint file for SurfaceModel"), std::string::npos);

    SurfaceModelFixture reader;
    auto reader_land = reader.make_fields(Real(-1.0), Real(1.0), 1);
    auto reader_urban = reader.make_fields(Real(-2.0), Real(1.0), 0);
    reader.configure_models(reader_land, reader_urban);
    reader.model->register_field_map("mapped", std::pair<int, int>{5, 5});
    reader.model->ReadCheckpoint(checkpoint.string());

    EXPECT_EQ(reader.model->get_ustar(0)->max(0), writer.model->get_ustar(0)->max(0));
    EXPECT_EQ(reader.model->get_ustar(0)->max(1), writer.model->get_ustar(0)->max(1));
    EXPECT_EQ(reader.model->get_tstar(0)->max(0), writer.model->get_tstar(0)->max(0));
    EXPECT_EQ(reader.model->get_qstar(0)->max(0), writer.model->get_qstar(0)->max(0));
    EXPECT_EQ(reader.model->get_tsurf(0)->max(0), writer.model->get_tsurf(0)->max(0));
    EXPECT_EQ(reader.model->get_wavg_factors(0)->max(SurfaceModelType::LAND),
              writer.model->get_wavg_factors(0)->max(SurfaceModelType::LAND));
    EXPECT_EQ(reader.model->get_wavg_factors(0)->max(SurfaceModelType::URBAN),
              writer.model->get_wavg_factors(0)->max(SurfaceModelType::URBAN));
    ASSERT_NE(reader.model->get_field("mapped"), nullptr);
    EXPECT_EQ(reader.model->get_field("mapped")->max(0),
              writer.model->get_field("mapped")->max(0));

}

} // namespace

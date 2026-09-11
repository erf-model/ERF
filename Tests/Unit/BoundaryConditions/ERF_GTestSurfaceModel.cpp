#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>
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
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> land_masks;
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>> masks;
    SolverChoice solver_choice{};
    std::unique_ptr<SurfaceModel> model;

    SurfaceModelFixture()
    {
        amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                {AMREX_D_DECL(2.0, 2.0, 4.0)});
        std::array<int, AMREX_SPACEDIM> periodic{};
        geom = amrex::Geometry(domain, &real_box, 0, periodic.data());
        geoms.push_back(geom);

        amrex::BoxList surface_boxes = ba.boxList();
        for (auto& box : surface_boxes) {
            box.setRange(2, 0);
        }
        amrex::BoxArray surface_ba(std::move(surface_boxes));
        land_masks.emplace_back(std::make_unique<amrex::iMultiFab>(
            surface_ba, dm, 1, amrex::IntVect(0)));
        land_masks[0]->setVal(1);
        masks.push_back(std::move(land_masks));

        model = std::make_unique<SurfaceModel>(1, bas, geoms, dms,
                                               solver_choice, masks);
        model->initialize_for_level(0, ba, geom, dm, masks[0],
                                    {}, {amrex::IntVect(2)});
    }

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> make_fields(
        const Real first, const Real increment, const int ngz)
    {
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> fields(6);
        for (int i = 0; i < 6; ++i) {
            fields[i] = std::make_unique<amrex::MultiFab>(
                ba, dm, 1, amrex::IntVect(1, 1, ngz));
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
        const amrex::Vector<int> field_indices{0, 1, 2, 3, 4, 5};
        model->set_model_data(0, pointers(land), {"f0", "f1", "f2", "f3", "f4", "f5"},
                              SurfaceModelType::LAND);
        model->set_model_data(0, pointers(urban), {"f0", "f1", "f2", "f3", "f4", "f5"},
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
    land_data.setVal(Real(4.0));
    urban_data.setVal(Real(8.0));
    fixture.model->apply_weight_average(0, &land_data, &urban_data);

    EXPECT_NEAR(land_data.max(0), Real(3.0), 1.e-12);
    EXPECT_NEAR(urban_data.max(0), Real(2.0), 1.e-12);
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

TEST(SurfaceModel, RadiationFieldListUsesCommonNamesAndProviderMappings)
{
    const std::vector<std::string> radiation_names{
        "tskin", "emiss", "albedo_vis", "albedo_nir",
        "albedo_vis_diff", "albedo_nir_diff"};

    {
        SurfaceModelFixture fixture;
        for (int i = 0; i < radiation_names.size(); ++i) {
            fixture.model->register_field_map(
                radiation_names[i], std::pair<int, int>{i, -1});
        }

        const auto radiation_fields = fixture.model->get_radiation_fields(0);
        ASSERT_EQ(radiation_fields.size(), radiation_names.size());
        for (int i = 0; i < radiation_names.size(); ++i) {
            EXPECT_EQ(radiation_fields[i],
                      fixture.model->get_field(radiation_names[i], 0));
            EXPECT_NE(radiation_fields[i], nullptr);
        }
    }

    // The common radiation names can map to different provider-specific
    // indices, as they do for the land and urban surface models.
    SurfaceModelFixture fixture;
    const std::vector<std::pair<int, int>> provider_mappings{
        {10, 20}, {11, 21}, {12, -1}, {13, -1}, {14, -1}, {15, -1}};
    for (int i = 0; i < radiation_names.size(); ++i) {
        fixture.model->register_field_map(radiation_names[i], provider_mappings[i]);
    }

    const auto radiation_fields = fixture.model->get_radiation_fields(0);
    ASSERT_EQ(radiation_fields.size(), radiation_names.size());
    for (int i = 0; i < radiation_names.size(); ++i) {
        EXPECT_EQ(radiation_fields[i],
                  fixture.model->get_field(radiation_names[i], 0));
        EXPECT_NE(radiation_fields[i], nullptr);
    }

    // Missing registered fields are represented by nullptr so the radiation
    // backend can use its default RRTMGP values.
    SurfaceModelFixture fallback_fixture;
    fallback_fixture.model->register_field_map(
        radiation_names[0], std::pair<int, int>{0, -1});

    const auto fallback_fields = fallback_fixture.model->get_radiation_fields(0);
    ASSERT_EQ(fallback_fields.size(), radiation_names.size());
    EXPECT_EQ(fallback_fields[0],
              fallback_fixture.model->get_field(radiation_names[0], 0));
    EXPECT_NE(fallback_fields[0], nullptr);
    for (int i = 1; i < fallback_fields.size(); ++i) {
        EXPECT_EQ(fallback_fields[i], nullptr);
    }
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

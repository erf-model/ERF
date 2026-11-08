/**
 * @file ERF_GTestStructureRegistry.cpp
 * @brief The structure-model registry: a model that registers itself with
 *        ERF_REGISTER_STRUCTURE_MODEL is listed and created by name
 *
 * These tests register models from this file alone, as a model directory outside
 * ERF does, and check that erf.structures.model finds them, that an unknown name
 * stops the run with the names it accepts, that a name registered twice or an
 * empty name stops the run, that the defaults of the optional methods hold
 * (no restart support, no plot fields, zeros from fill_plot_var), and that plot
 * names that are empty, repeated or already taken stop the run.
 */
#include <gtest/gtest.h>

#include "ERF_GTestThrowOnAbort.H"

#include <ERF_StructureModelRegistry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <utility>

using namespace amrex;

namespace {

/// The required methods only: every optional one keeps its default.
class RegistryTestModel : public StructureModel
{
public:
    explicit RegistryTestModel (int max_level) : m_max_level(max_level) {}
    const char* name () const override { return "registry-test"; }
    int level (int finest_level) const override { return finest_level; }
    void initialize (const StructureInit&) override {}
    void advance (const StructureStep&) override {}
    std::vector<RealBox> support () const override { return {}; }
    void add_momentum_sources (int, MultiFab&, MultiFab&, MultiFab&) const override {}
    int m_max_level;
};

/// A model whose plot names a test sets
class PlotNamesTestModel : public RegistryTestModel
{
public:
    explicit PlotNamesTestModel (std::vector<std::string> names) : RegistryTestModel(0), m_names(std::move(names)) {}
    std::vector<std::string> plot_var_names () const override { return m_names; }
    std::vector<std::string> m_names;
};

class RegistryTestModelB : public RegistryTestModel
{
public:
    explicit RegistryTestModelB (int max_level) : RegistryTestModel(max_level) {}
    const char* name () const override { return "registry-test-b"; }
};

/// Sets erf.structures.model for one test and removes it afterwards
class ScopedModelKey
{
public:
    explicit ScopedModelKey (const std::string& v) : m_pp("erf.structures") { m_pp.add("model", v); }
    ~ScopedModelKey () { m_pp.remove("model"); }
    ScopedModelKey (const ScopedModelKey&) = delete;
    ScopedModelKey& operator= (const ScopedModelKey&) = delete;
private:
    ParmParse m_pp;
};

} // namespace

ERF_REGISTER_STRUCTURE_MODEL(RegistryTestModel, "registry-test")
ERF_REGISTER_STRUCTURE_MODEL(RegistryTestModelB, "registry-test-b")

TEST(StructureRegistry, ListsSelfRegisteredModelsSorted)
{
    const auto names = structure_model_names();
    EXPECT_NE(std::find(names.begin(), names.end(), "registry-test"), names.end());
    EXPECT_NE(std::find(names.begin(), names.end(), "registry-test-b"), names.end());
    EXPECT_TRUE(std::is_sorted(names.begin(), names.end()));
}

TEST(StructureRegistry, CreatesByNameWithMaxLevel)
{
    auto model = make_structure_model("registry-test-b", 3);
    ASSERT_NE(model, nullptr);
    EXPECT_STREQ(model->name(), "registry-test-b");
    EXPECT_EQ(dynamic_cast<RegistryTestModel*>(model.get())->m_max_level, 3);
    EXPECT_EQ(make_structure_model("no-such-model", 0), nullptr);
}

TEST(StructureRegistry, DuplicateAndEmptyNamesAbort)
{
    auto creator = [](int ml) -> std::unique_ptr<StructureModel> { return std::make_unique<RegistryTestModel>(ml); };
    const std::string dup = erf_gtest::abort_message([&]{ register_structure_model("registry-test", creator); });
    EXPECT_NE(dup.find("two structure models registered the name \"registry-test\""), std::string::npos) << dup;
    const std::string empty = erf_gtest::abort_message([&]{ register_structure_model("", creator); });
    EXPECT_NE(empty.find("registered an empty name"), std::string::npos) << empty;
}

TEST(StructureRegistry, InputsSelectTheModel)
{
    // Unset: no model, no abort
    EXPECT_EQ(make_structure_model_from_inputs(0), nullptr);

    {
        ScopedModelKey key("registry-test");
        auto model = make_structure_model_from_inputs(2);
        ASSERT_NE(model, nullptr);
        EXPECT_STREQ(model->name(), "registry-test");
        EXPECT_EQ(dynamic_cast<RegistryTestModel*>(model.get())->m_max_level, 2);
    }
    {
        ScopedModelKey key("no-such-model");
        const std::string msg = erf_gtest::abort_message([]{ (void)make_structure_model_from_inputs(0); });
        EXPECT_NE(msg.find("erf.structures.model = \"no-such-model\" is not one of the structure models compiled in"),
                  std::string::npos) << msg;
        EXPECT_NE(msg.find("\"registry-test\", \"registry-test-b\""), std::string::npos) << msg;
    }
}

TEST(StructureRegistry, OptionalMethodDefaults)
{
    RegistryTestModel model(0);
    EXPECT_FALSE(model.supports_restart());
    EXPECT_FALSE(model.wants_immersed_terrain_surface());
    EXPECT_TRUE(model.plot_var_names().empty());

    BoxArray ba(Box(IntVect(0), IntVect(7)));
    ba.maxSize(4);
    DistributionMapping dm(ba);
    MultiFab dst(ba, dm, 2, 0);
    dst.setVal(Real(7.0));
    model.fill_plot_var("anything", 0, dst, 1);
    EXPECT_EQ(dst.min(1), Real(0.0));
    EXPECT_EQ(dst.max(1), Real(0.0));
    EXPECT_EQ(dst.min(0), Real(7.0));   // the other component is left alone
}

TEST(StructureRegistry, PlotNamesMustBeFreeAndDistinct)
{
    const std::vector<std::string> taken{"density", "x_velocity", "fire_phi"};
    // Distinct names that ERF does not use pass
    check_structure_plot_names(PlotNamesTestModel({"test_fx", "test_fy"}), taken);

    const std::string clash = erf_gtest::abort_message([&]{
        check_structure_plot_names(PlotNamesTestModel({"test_fx", "x_velocity"}), taken); });
    EXPECT_NE(clash.find("lists the plot variable \"x_velocity\", which ERF or another model already provides"),
              std::string::npos) << clash;
    const std::string twice = erf_gtest::abort_message([&]{
        check_structure_plot_names(PlotNamesTestModel({"test_fx", "test_fx"}), taken); });
    EXPECT_NE(twice.find("lists the plot variable \"test_fx\" twice"), std::string::npos) << twice;
    const std::string empty = erf_gtest::abort_message([&]{
        check_structure_plot_names(PlotNamesTestModel({""}), taken); });
    EXPECT_NE(empty.find("lists an empty plot variable name"), std::string::npos) << empty;
}

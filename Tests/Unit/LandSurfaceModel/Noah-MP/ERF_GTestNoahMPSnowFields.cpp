// Noah-MP's SNOW and SNOWH reach ERF through the field registry (ERF_NOAHMP_Fields.H).
// These checks read the registry only through its X-macros, so they compile before the
// fields are registered and fail at run time, not at compile time.
#include <algorithm>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_NOAHMP_ResultPolicy.H"

using amrex::Real;

namespace
{

#define SNOW_TEST_NAME1(comp)              std::string(#comp),
#define SNOW_TEST_NAME2(comp, member)      std::string(#comp),
#define SNOW_TEST_PAIR(comp, member)       std::make_pair(std::string(#comp), std::string(#member)),
#define SNOW_TEST_TRIPLE(alias, lsm, out)  std::make_tuple(std::string(#alias), std::string(#lsm), std::string(#out)),

std::vector<std::string> lsm_data_names ()
{
    return {NOAHMP_LSMDATA_FIELDS(SNOW_TEST_NAME1)};
}

// Mirrors the NoahmpOutputComp enum: the head rows, the four banded albedos, the tail
// rows, then the two sentinels.
std::vector<std::string> output_comp_names ()
{
    std::vector<std::string> names = {NOAHMP_OUTPUT_2D_FIELDS_HEAD(SNOW_TEST_NAME2)};
    for (const char* band : {"albsfcdir_vis", "albsfcdir_nir", "albsfcdif_vis", "albsfcdif_nir"}) {
        names.emplace_back(band);
    }
    for (const auto& name : std::vector<std::string>{NOAHMP_OUTPUT_2D_FIELDS_TAIL(SNOW_TEST_NAME2)}) {
        names.push_back(name);
    }
    names.emplace_back("o_smstav");
    names.emplace_back("o_smstot");
    return names;
}

int output_comp_index (const std::string& name)
{
    const auto names = output_comp_names();
    const auto it = std::find(names.begin(), names.end(), name);
    return it == names.end() ? -1 : static_cast<int>(it - names.begin());
}

} // namespace

TEST(NoahMPSnowFields, OutputMirrorMatchesTheEnum)
{
    // Guards the helper above: if the enum layout changes, the index lookups below are wrong.
    EXPECT_EQ(output_comp_index("o_smstav"), NoahmpOutputComp::o_smstav);
    EXPECT_EQ(static_cast<int>(output_comp_names().size()), NoahmpOutputComp::NumComps);
}

TEST(NoahMPSnowFields, LsmDataEndsWithSnowThenSnowh)
{
    const auto names = lsm_data_names();
    ASSERT_GE(names.size(), 2u);
    // Appended at the end, per the registry's "do not reorder" rule.
    EXPECT_EQ(names[names.size() - 2], "snow");
    EXPECT_EQ(names[names.size() - 1], "snowh");
    EXPECT_EQ(static_cast<int>(names.size()), LsmData_NOAHMP::NumVars);
}

TEST(NoahMPSnowFields, LsmDataCountIsTheNewCheckpointLayout)
{
    // 27 before SNOW/SNOWH were registered. A checkpoint written with 27 fields will not
    // restart with 29: the soil LsmData<i> indices move up by two.
    EXPECT_EQ(LsmData_NOAHMP::NumVars, 29);
}

TEST(NoahMPSnowFields, OutputRowsReadTheNoahmpIoSnowMembers)
{
    const std::vector<std::pair<std::string, std::string>> rows = {NOAHMP_OUTPUT_2D_FIELDS(SNOW_TEST_PAIR)};
    const auto has = [&rows] (const std::string& comp, const std::string& member) {
        return std::find(rows.begin(), rows.end(), std::make_pair(comp, member)) != rows.end();
    };
    EXPECT_TRUE(has("o_snow", "SNOW"));
    EXPECT_TRUE(has("o_snowh", "SNOWH"));
}

TEST(NoahMPSnowFields, ResultRowsLinkOutputsToLsmData)
{
    const std::vector<std::tuple<std::string, std::string, std::string>> rows = {NOAHMP_RESULT_FIELDS(SNOW_TEST_TRIPLE)};
    const auto has = [&rows] (const std::string& alias, const std::string& lsm, const std::string& out) {
        return std::find(rows.begin(), rows.end(), std::make_tuple(alias, lsm, out)) != rows.end();
    };
    EXPECT_TRUE(has("SNOW_o", "snow", "o_snow"));
    EXPECT_TRUE(has("SNOWH_o", "snowh", "o_snowh"));
}

TEST(NoahMPSnowFields, ZeroSnowIsStoredAndUnprocessedCellsAreNot)
{
    const noahmp_result_policy::CellPolicy processed(Real(10.0));     // a valid HFX
    const noahmp_result_policy::CellPolicy unprocessed(Real(-9999.0)); // Noah-MP's skip marker
    for (const char* comp : {"o_snow", "o_snowh"}) {
        const int idx = output_comp_index(comp);
        ASSERT_GE(idx, 0) << comp << " is not a NoahmpOutputComp";
        EXPECT_EQ(processed.select(Real(0.0), lsm_undefined, idx), Real(0.0)) << comp;
        EXPECT_EQ(processed.select(Real(1.25), lsm_undefined, idx), Real(1.25)) << comp;
        EXPECT_EQ(processed.select(Real(-9999.0), lsm_undefined, idx), lsm_undefined) << comp;
        EXPECT_EQ(unprocessed.select(Real(1.25), lsm_undefined, idx), lsm_undefined) << comp;
    }
}

#undef SNOW_TEST_NAME1
#undef SNOW_TEST_NAME2
#undef SNOW_TEST_PAIR
#undef SNOW_TEST_TRIPLE

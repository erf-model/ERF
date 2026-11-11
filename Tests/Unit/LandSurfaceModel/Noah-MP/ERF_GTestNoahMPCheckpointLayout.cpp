// A checkpoint written before Noah-MP's snow and snowh fields existed must still
// restart: every other LSM field keeps the file it was written to, and only the two
// new fields may be absent.
#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_NOAHMP_ResultPolicy.H"
#include "ERF_CheckpointLsmData.H"

namespace
{

#define CKPT_TEST_NAME(comp) std::string(#comp),

bool is_snow (const std::string& name) { return name == "snow" || name == "snowh"; }

// The registry as it is now, then three soil profiles appended at run time.
std::vector<std::string> lsm_data_names ()
{
    std::vector<std::string> names = {NOAHMP_LSMDATA_FIELDS(CKPT_TEST_NAME)};
    for (const char* group : {"smois", "sh2o", "tslb"}) {
        for (int k = 0; k < 4; ++k) {
            names.push_back(std::string(group) + "_" + std::to_string(k));
        }
    }
    return names;
}

// The registry before snow and snowh: the same list without those two rows.
std::vector<std::string> legacy_names ()
{
    std::vector<std::string> names;
    for (const auto& name : lsm_data_names()) {
        if (!is_snow(name)) { names.push_back(name); }
    }
    return names;
}

} // namespace

// Motivation: soil fields sit after the registry, so inserting two rows moved
// their index. Each non-snow field must keep the index it had before.
TEST(NoahMPCheckpointLayout, EveryFieldKeepsItsLegacyFileIndex)
{
    const auto names = lsm_data_names();
    const auto legacy = legacy_names();
    ASSERT_EQ(names.size(), legacy.size() + 2);
    for (std::size_t i = 0; i < names.size(); ++i) {
        if (is_snow(names[i])) { continue; }
        const auto legacy_index = std::find(legacy.begin(), legacy.end(), names[i]) - legacy.begin();
        EXPECT_EQ(erf_checkpoint_lsm_data::file_stem(names, i),
                  "LsmData" + std::to_string(legacy_index)) << names[i];
    }
}

// Motivation: a restart from a checkpoint written before the fields existed
// must find every file it asks for, except the two new ones.
TEST(NoahMPCheckpointLayout, OlderCheckpointHasEveryRequiredFile)
{
    const auto names = lsm_data_names();
    std::set<std::string> older_checkpoint;
    for (std::size_t i = 0; i < legacy_names().size(); ++i) {
        older_checkpoint.insert("LsmData" + std::to_string(i));
    }
    for (std::size_t i = 0; i < names.size(); ++i) {
        const std::string stem = erf_checkpoint_lsm_data::file_stem(names, i);
        if (erf_checkpoint_lsm_data::may_be_absent(names, i)) {
            EXPECT_TRUE(is_snow(names[i]));
            EXPECT_EQ(older_checkpoint.count(stem), 0u) << names[i];
        } else {
            EXPECT_EQ(older_checkpoint.count(stem), 1u) << names[i];
        }
    }
}

// Motivation: no two fields may share a file.
TEST(NoahMPCheckpointLayout, FileNamesAreDistinct)
{
    const auto names = lsm_data_names();
    std::set<std::string> stems;
    for (std::size_t i = 0; i < names.size(); ++i) {
        stems.insert(erf_checkpoint_lsm_data::file_stem(names, i));
    }
    EXPECT_EQ(stems.size(), names.size());
}

// Motivation: a land model without these fields keeps the plain numbering.
TEST(NoahMPCheckpointLayout, OtherLandModelsKeepPlainNumbering)
{
    const std::vector<std::string> names = {"tsurf", "ustar", "tstar"};
    for (std::size_t i = 0; i < names.size(); ++i) {
        EXPECT_EQ(erf_checkpoint_lsm_data::file_stem(names, i), "LsmData" + std::to_string(i));
        EXPECT_FALSE(erf_checkpoint_lsm_data::may_be_absent(names, i));
    }
}

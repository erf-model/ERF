#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "ERF_CheckpointSurfaceTemperature.H"

namespace {

using erf_checkpoint_surface_temperature::ContractReadStatus;
using erf_checkpoint_surface_temperature::LegacyInitType;
using erf_checkpoint_surface_temperature::LegacySurfaceTemperatureCompatibility;

struct TemporaryCheckpoint
{
    std::filesystem::path root;

    TemporaryCheckpoint ()
        : root(std::filesystem::temp_directory_path() /
               ("erf-checkpoint-surface-temperature-" + amrex::UniqueString()))
    {
        std::filesystem::create_directories(root);
    }

    ~TemporaryCheckpoint ()
    {
        std::filesystem::remove_all(root);
    }
};

void
touch_surface_temperature_file (const std::filesystem::path& root,
                                 const int level,
                                 const char* name)
{
    const auto level_dir = root / ("Level_" + std::to_string(level));
    std::filesystem::create_directories(level_dir);
    const auto filename = amrex::MultiFabFileFullPrefix(
        level, root.string(), "Level_", name);
    std::ofstream output(filename);
    ASSERT_TRUE(output.good());
}

} // namespace

TEST(CheckpointSurfaceTemperature, ReadsVersionOneMarker)
{
    // Motivation: every new native checkpoint must explicitly identify the
    // surface-temperature representation so a future reader cannot infer it
    // from the presence of optional SST/TSK arrays.
    std::ostringstream written_marker;
    erf_checkpoint_surface_temperature::write_contract_version(written_marker);
    int version = 0;
    std::istringstream marker(written_marker.str());
    EXPECT_EQ(erf_checkpoint_surface_temperature::read_contract_version(marker, version),
              ContractReadStatus::Valid);
    EXPECT_EQ(version, erf_checkpoint_surface_temperature::contract_version);
}

TEST(CheckpointSurfaceTemperature, RejectsUnknownOrMalformedMarker)
{
    // Motivation: accepting a future contract version without understanding
    // its representation could silently restore temperatures with the wrong
    // thermodynamic convention.
    int version = 0;
    std::istringstream unknown_version("2\n");
    EXPECT_EQ(erf_checkpoint_surface_temperature::read_contract_version(
                  unknown_version, version),
              ContractReadStatus::UnknownVersion);

    for (const char* text : {"not-a-version\n", "1 extra\n"}) {
        version = 0;
        std::istringstream malformed(text);
        EXPECT_EQ(erf_checkpoint_surface_temperature::read_contract_version(
                      malformed, version),
                  ContractReadStatus::Malformed);
    }
}

TEST(CheckpointSurfaceTemperature, ParsesLegacyMetgridInitTypeFromJobInfo)
{
    std::istringstream job_info("unrelated = 1\nerf.init_type = \"Metgrid\"\n");
    EXPECT_EQ(erf_checkpoint_surface_temperature::parse_legacy_init_type_from_job_info(job_info),
              LegacyInitType::Metgrid);
}

TEST(CheckpointSurfaceTemperature, ParsesLegacyInitTypeValueCaseInsensitively)
{
    EXPECT_EQ(erf_checkpoint_surface_temperature::parse_legacy_init_type_value("wRfInPuT"),
              LegacyInitType::WRFInput);
    EXPECT_EQ(erf_checkpoint_surface_temperature::parse_legacy_init_type_value("METGRID"),
              LegacyInitType::Metgrid);
}

TEST(CheckpointSurfaceTemperature, RejectsInvalidLegacyInitTypeValue)
{
    EXPECT_EQ(erf_checkpoint_surface_temperature::parse_legacy_init_type_value("other"),
              LegacyInitType::Unknown);
}

TEST(CheckpointSurfaceTemperature, ParsesLegacyWrfInputInitTypeCaseInsensitively)
{
    std::istringstream job_info("erf.init_type = wrfinput\n");
    EXPECT_EQ(erf_checkpoint_surface_temperature::parse_legacy_init_type_from_job_info(job_info),
              LegacyInitType::WRFInput);
}

TEST(CheckpointSurfaceTemperature, UsesLastLegacyInitTypeAssignment)
{
    std::istringstream job_info(
        "erf.init_type = \"WRFInput\"\nerf.init_type = Metgrid\n");
    EXPECT_EQ(erf_checkpoint_surface_temperature::parse_legacy_init_type_from_job_info(job_info),
              LegacyInitType::Metgrid);
}

TEST(CheckpointSurfaceTemperature, UnknownWhenLegacyInitTypeIsMissingOrMalformed)
{
    for (const char* text : {"other.key = Metgrid\n", "erf.init_type =\n",
                             "erf.init_type = \"Metgrid\n", "erf.init_type Metgrid\n",
                             "erf.init_type = Other\n"}) {
        std::istringstream job_info(text);
        EXPECT_EQ(erf_checkpoint_surface_temperature::parse_legacy_init_type_from_job_info(job_info),
                  LegacyInitType::Unknown);
    }
}

TEST(CheckpointSurfaceTemperature, AllowsLegacyWrfInputSurfaceArraysWithoutMarker)
{
    EXPECT_EQ(erf_checkpoint_surface_temperature::classify_legacy_surface_temperature_checkpoint(
                  false, true, LegacyInitType::WRFInput),
              LegacySurfaceTemperatureCompatibility::Compatible);
}

TEST(CheckpointSurfaceTemperature, RejectsLegacyMetgridSurfaceArraysWithoutMarker)
{
    EXPECT_EQ(erf_checkpoint_surface_temperature::classify_legacy_surface_temperature_checkpoint(
                  false, true, LegacyInitType::Metgrid),
              LegacySurfaceTemperatureCompatibility::UnsafeMetgrid);
}

TEST(CheckpointSurfaceTemperature, RejectsLegacySurfaceArraysWithUnknownProvenance)
{
    EXPECT_EQ(erf_checkpoint_surface_temperature::classify_legacy_surface_temperature_checkpoint(
                  false, true, LegacyInitType::Unknown),
              LegacySurfaceTemperatureCompatibility::UnknownProvenance);
}

TEST(CheckpointSurfaceTemperature, AllowsMarkerlessCheckpointWithoutSurfaceArrays)
{
    EXPECT_EQ(erf_checkpoint_surface_temperature::classify_legacy_surface_temperature_checkpoint(
                  false, false, LegacyInitType::Unknown),
              LegacySurfaceTemperatureCompatibility::Compatible);
}

TEST(CheckpointSurfaceTemperature, ScansEveryAMRLevelForLegacySurfaceArrays)
{
    // Motivation: legacy SST/TSK files can occur on a refined AMR level, so a
    // level-zero-only compatibility check could accept an unsafe checkpoint.
    TemporaryCheckpoint checkpoint;

    touch_surface_temperature_file(checkpoint.root, 1, "SST_0_H");
    EXPECT_EQ(erf_checkpoint_surface_temperature::first_legacy_surface_temperature_level(
                  checkpoint.root.string(), 2), 1);

    std::filesystem::remove(amrex::MultiFabFileFullPrefix(
        1, checkpoint.root.string(), "Level_", "SST_0_H"));
    touch_surface_temperature_file(checkpoint.root, 1, "TSK_0_H");
    EXPECT_EQ(erf_checkpoint_surface_temperature::first_legacy_surface_temperature_level(
                  checkpoint.root.string(), 2), 1);

    std::filesystem::remove(amrex::MultiFabFileFullPrefix(
        1, checkpoint.root.string(), "Level_", "TSK_0_H"));
    EXPECT_EQ(erf_checkpoint_surface_temperature::first_legacy_surface_temperature_level(
                  checkpoint.root.string(), 2), -1);

    touch_surface_temperature_file(checkpoint.root, 1, "SST_0_H");
    std::ostringstream written_marker;
    erf_checkpoint_surface_temperature::write_contract_version(written_marker);
    int version = 0;
    std::istringstream marker(written_marker.str());
    EXPECT_EQ(erf_checkpoint_surface_temperature::read_contract_version(marker, version),
              ContractReadStatus::Valid);
    EXPECT_EQ(version, erf_checkpoint_surface_temperature::contract_version);
}

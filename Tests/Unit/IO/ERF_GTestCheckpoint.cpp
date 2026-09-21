#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "ERF_CheckpointSurfaceTemperature.H"

namespace {

using erf_checkpoint_surface_temperature::ContractReadStatus;

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

TEST(CheckpointSurfaceTemperature, AllowsCompatibleLegacyCheckpoint)
{
    // Motivation: the compatibility guard is specific to legacy Metgrid
    // surface-temperature arrays; older WRFInput checkpoints and checkpoints
    // without SST/TSK remain readable.
    EXPECT_FALSE(erf_checkpoint_surface_temperature::legacy_metgrid_surface_temperature_is_unsafe(
        false, false, true, true));
    EXPECT_FALSE(erf_checkpoint_surface_temperature::legacy_metgrid_surface_temperature_is_unsafe(
        false, true, false, false));
    EXPECT_FALSE(erf_checkpoint_surface_temperature::legacy_metgrid_surface_temperature_is_unsafe(
        true, true, true, true));
}

TEST(CheckpointSurfaceTemperature, RejectsLegacyMetgridSurfaceArraysWithoutMarker)
{
    // Motivation: legacy Metgrid SST_0/TSK_0 arrays were written as absolute
    // temperature, while current Metgrid initialization expects theta and
    // cannot reconstruct the conversion pressure from a checkpoint alone.
    EXPECT_TRUE(erf_checkpoint_surface_temperature::legacy_metgrid_surface_temperature_is_unsafe(
        false, true, true, false));
    EXPECT_TRUE(erf_checkpoint_surface_temperature::legacy_metgrid_surface_temperature_is_unsafe(
        false, true, false, true));
}

TEST(CheckpointSurfaceTemperature, ScansEveryAMRLevelForLegacySurfaceArrays)
{
    // Motivation: legacy SST/TSK files can occur on a refined AMR level, so a
    // level-zero-only compatibility check could accept an unsafe checkpoint.
    TemporaryCheckpoint checkpoint;

    touch_surface_temperature_file(checkpoint.root, 1, "SST_0_H");
    EXPECT_EQ(erf_checkpoint_surface_temperature::first_legacy_surface_temperature_level(
                  checkpoint.root.string(), 2), 1);
    EXPECT_FALSE(erf_checkpoint_surface_temperature::legacy_metgrid_surface_temperature_is_unsafe(
        false, false, true, false));

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

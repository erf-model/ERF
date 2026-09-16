#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "ERF_CheckpointSurfaceTemperature.H"

namespace {

using erf_checkpoint_surface_temperature::ContractReadStatus;

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

    for (const std::string& text : {"not-a-version\n", "1 extra\n"}) {
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

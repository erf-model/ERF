#include <algorithm>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "ERF_CheckpointSurfaceTemperature.H"
#include "ERF_CheckpointLandmask.H"

namespace {

using erf_checkpoint_surface_temperature::ContractReadStatus;
using erf_checkpoint_surface_temperature::LegacyInitType;
using erf_checkpoint_surface_temperature::StoredSurfaceTemperature;

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

// Set every cell of a single-component MultiFab, ghost cells included, from the host
template <typename F>
void
set_from_host (amrex::MultiFab& mf, const F& value)
{
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        amrex::FArrayBox host(mf[mfi].box(), 1, amrex::The_Pinned_Arena());
        const auto arr = host.array();
        amrex::LoopOnCpu(host.box(), [&] (int i, int j, int k) { arr(i,j,k) = value(i,j,k); });
        amrex::Gpu::htod_memcpy(mf[mfi].dataPtr(), host.dataPtr(),
                                host.size() * sizeof(amrex::Real));
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace

TEST(CheckpointLandmask, RestoresLegacyCoarseBuildingColumnsToLand)
{
    const amrex::Box box(amrex::IntVect(0,0,0), amrex::IntVect(2,0,0));
    const amrex::BoxArray ba(box);
    const amrex::DistributionMapping dm(ba);
    amrex::iMultiFab mask(ba, dm, 1, amrex::IntVect(1,0,0));
    mask.setVal(2);
    mask[0].setVal<amrex::RunOn::Device>(
        0, amrex::Box(amrex::IntVect(0,0,0), amrex::IntVect(0,0,0)), 0, 1);
    mask[0].setVal<amrex::RunOn::Device>(
        1, amrex::Box(amrex::IntVect(1,0,0), amrex::IntVect(1,0,0)), 0, 1);

    erf_checkpoint_landmask::restore_coarse_building_landmask(mask);

    amrex::IArrayBox host(mask[0].box(), 1, amrex::The_Pinned_Arena());
    amrex::Gpu::dtoh_memcpy(host.dataPtr(), mask[0].dataPtr(), host.size() * sizeof(int));
    amrex::Gpu::streamSynchronize();
    const auto values = host.const_array();
    EXPECT_EQ(values(-1,0,0), 1); // old building mask in a ghost cell
    EXPECT_EQ(values(0,0,0), 0);  // water remains water
    EXPECT_EQ(values(1,0,0), 1);  // land remains land
    EXPECT_EQ(values(2,0,0), 1);  // old building mask becomes land
    EXPECT_EQ(values(3,0,0), 1);
}

TEST(CheckpointSurfaceTemperature, ReadsCurrentMarker)
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

TEST(CheckpointSurfaceTemperature, ReadsVersionOneMarker)
{
    // Motivation: checkpoints written while SST/TSK were kept as potential
    // temperature carry version 1; they must still restart, converted.
    int version = 0;
    std::istringstream marker("1\n");
    EXPECT_EQ(erf_checkpoint_surface_temperature::read_contract_version(marker, version),
              ContractReadStatus::Valid);
    EXPECT_EQ(version, 1);
}

TEST(CheckpointSurfaceTemperature, RejectsUnknownOrMalformedMarker)
{
    // Motivation: accepting a future contract version without understanding
    // its representation could silently restore temperatures with the wrong
    // thermodynamic convention.
    for (const char* text : {"0\n", "3\n"}) {
        int version = 0;
        std::istringstream unknown_version(text);
        EXPECT_EQ(erf_checkpoint_surface_temperature::read_contract_version(
                      unknown_version, version),
                  ContractReadStatus::UnknownVersion);
    }

    for (const char* text : {"not-a-version\n", "1 extra\n"}) {
        int version = 0;
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

TEST(CheckpointSurfaceTemperature, ConvertsLegacyWrfInputSurfaceArraysWithoutMarker)
{
    EXPECT_EQ(erf_checkpoint_surface_temperature::classify_stored_surface_temperature(
                  0, true, LegacyInitType::WRFInput),
              StoredSurfaceTemperature::Potential);
}

TEST(CheckpointSurfaceTemperature, RestoresLegacyMetgridSurfaceArraysWithoutMarker)
{
    // Motivation: Metgrid stored the absolute temperature before the marker
    // existed, which is what the arrays hold now.
    EXPECT_EQ(erf_checkpoint_surface_temperature::classify_stored_surface_temperature(
                  0, true, LegacyInitType::Metgrid),
              StoredSurfaceTemperature::Absolute);
}

TEST(CheckpointSurfaceTemperature, RejectsLegacySurfaceArraysWithUnknownProvenance)
{
    EXPECT_EQ(erf_checkpoint_surface_temperature::classify_stored_surface_temperature(
                  0, true, LegacyInitType::Unknown),
              StoredSurfaceTemperature::UnknownProvenance);
}

TEST(CheckpointSurfaceTemperature, AllowsMarkerlessCheckpointWithoutSurfaceArrays)
{
    EXPECT_EQ(erf_checkpoint_surface_temperature::classify_stored_surface_temperature(
                  0, false, LegacyInitType::Unknown),
              StoredSurfaceTemperature::Absolute);
}

TEST(CheckpointSurfaceTemperature, ClassifiesSurfaceArraysByMarkerVersion)
{
    // Motivation: the marker, when present, decides alone; version 1 arrays
    // are potential temperature whichever init type wrote them.
    for (const auto init_type : {LegacyInitType::WRFInput, LegacyInitType::Metgrid,
                                 LegacyInitType::Unknown}) {
        EXPECT_EQ(erf_checkpoint_surface_temperature::classify_stored_surface_temperature(
                      1, true, init_type),
                  StoredSurfaceTemperature::Potential);
        EXPECT_EQ(erf_checkpoint_surface_temperature::classify_stored_surface_temperature(
                      erf_checkpoint_surface_temperature::contract_version, true, init_type),
                  StoredSurfaceTemperature::Absolute);
        EXPECT_EQ(erf_checkpoint_surface_temperature::classify_stored_surface_temperature(
                      1, false, init_type),
                  StoredSurfaceTemperature::Absolute);
    }
}

TEST(CheckpointSurfaceTemperature, ConvertsPotentialTemperatureWithEachColumnsPressure)
{
    // Motivation: a version 1 checkpoint holds theta = T (p0/p)^(Rd/cp) with
    // the surface pressure p of each column. Restoring it must give back T
    // column by column, and a ghost cell, whose pressure is never set, must
    // take the pressure of the nearest valid cell.
    const amrex::Box box(amrex::IntVect(0,0,0), amrex::IntVect(3,2,0));
    const amrex::BoxArray ba(box);
    const amrex::DistributionMapping dm(ba);
    const amrex::IntVect ng(1,1,0);
    amrex::MultiFab field(ba, dm, 1, ng);
    amrex::MultiFab pressure(ba, dm, 1, ng);

    const auto temperature_at = [] (int i, int j) {
        return amrex::Real(280.0) + amrex::Real(2.0) * i + amrex::Real(0.5) * j;
    };
    const auto pressure_at = [] (int i, int j) {
        return amrex::Real(95000.0) + amrex::Real(1000.0) * i - amrex::Real(700.0) * j;
    };
    const auto nearest_valid = [&] (int i, int j) {
        return amrex::IntVect(std::clamp(i, box.smallEnd(0), box.bigEnd(0)),
                              std::clamp(j, box.smallEnd(1), box.bigEnd(1)), 0);
    };

    set_from_host(pressure, [&] (int i, int j, int) {
        return box.contains(amrex::IntVect(i,j,0)) ? pressure_at(i,j)
                                                   : std::numeric_limits<amrex::Real>::quiet_NaN();
    });
    // The conversion version 1 made, in the ghost cells too
    set_from_host(field, [&] (int i, int j, int) {
        const amrex::IntVect v = nearest_valid(i,j);
        return getThgivenTandP(temperature_at(v[0],v[1]), pressure_at(v[0],v[1]), RdoCp);
    });

    erf_checkpoint_surface_temperature::potential_to_absolute_temperature(
        field, pressure, RdoCp, amrex::Periodicity::NonPeriodic());

    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
        amrex::FArrayBox host(field[mfi].box(), 1, amrex::The_Pinned_Arena());
        amrex::Gpu::dtoh_memcpy(host.dataPtr(), field[mfi].dataPtr(),
                                host.size() * sizeof(amrex::Real));
        amrex::Gpu::streamSynchronize();
        const auto values = host.const_array();
        amrex::LoopOnCpu(host.box(), [&] (int i, int j, int k) {
            const amrex::IntVect v = nearest_valid(i,j);
            const amrex::Real expected = temperature_at(v[0],v[1]);
            EXPECT_NEAR(values(i,j,k), expected,
                        amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon() * expected)
                << "at (" << i << "," << j << ")";
        });
    }
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

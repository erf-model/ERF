#include <algorithm>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_CheckpointSurfaceTemperature.H"
#include "ERF_CheckpointLandmask.H"

namespace {

namespace cst = erf_checkpoint_surface_temperature;

using erf_checkpoint_surface_temperature::ContractReadStatus;
using erf_checkpoint_surface_temperature::LegacyInitType;
using erf_checkpoint_surface_temperature::LegacySurfacePressureSource;
using erf_checkpoint_surface_temperature::ProvenanceOrigin;
using erf_checkpoint_surface_temperature::StoredSurfaceTemperature;
using erf_checkpoint_surface_temperature::SurfaceTemperatureMigration;
using erf_surface_temperature::SurfaceTemperatureSource;

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

// A land mask on the layout of a field, every cell, ghost cells included, set from the host
template <typename F>
amrex::iMultiFab
make_land_mask (const amrex::MultiFab& field, const F& value)
{
    amrex::iMultiFab mask(field.boxArray(), field.DistributionMap(), 1, field.nGrowVect());
    for (amrex::MFIter mfi(mask); mfi.isValid(); ++mfi) {
        amrex::IArrayBox host(mask[mfi].box(), 1, amrex::The_Pinned_Arena());
        const auto arr = host.array();
        amrex::LoopOnCpu(host.box(), [&] (int i, int j, int k) { arr(i,j,k) = value(i,j,k); });
        amrex::Gpu::htod_memcpy(mask[mfi].dataPtr(), host.dataPtr(),
                                host.size() * sizeof(int));
    }
    amrex::Gpu::streamSynchronize();
    return mask;
}

// A land mask on the layout of a field that is water everywhere
amrex::iMultiFab
make_water_mask (const amrex::MultiFab& field)
{
    return make_land_mask(field, [] (int, int, int) { return 0; });
}

void
write_text_file (const std::filesystem::path& path, const std::string& text)
{
    std::ofstream output(path);
    ASSERT_TRUE(output.good());
    output << text;
}

// A checkpoint as a restart finds it: SST_0 at level 0, with the contract marker of the
// given version (none for 0) and a job_info recording the given erf.init_type (none if
// empty)
void
make_checkpoint (const TemporaryCheckpoint& checkpoint,
                 const int marker_version,
                 const std::string& job_info_init_type)
{
    touch_surface_temperature_file(checkpoint.root, 0, "SST_0_H");
    if (marker_version > 0) {
        write_text_file(checkpoint.root / cst::contract_file_name,
                        std::to_string(marker_version) + "\n");
    }
    std::string job_info = " ERF Job Information\n\n Inputs File Parameters\n"
                           "amr.max_level = 0\n";
    if (!job_info_init_type.empty()) {
        job_info += "erf.init_type = \"" + job_info_init_type + "\"\n";
    }
    job_info += "erf.use_moist_background = 0\n";
    write_text_file(checkpoint.root / "job_info", job_info);
}

// Sets erf.legacy_surface_temperature_init_type for the life of the object
class ScopedLegacyInitTypeAssertion
{
public:
    explicit ScopedLegacyInitTypeAssertion (const std::string& value)
    {
        amrex::ParmParse pp("erf");
        pp.add(cst::legacy_init_type_key, value);
    }
    ~ScopedLegacyInitTypeAssertion ()
    {
        amrex::ParmParse pp("erf");
        pp.remove(cst::legacy_init_type_key);
    }
    ScopedLegacyInitTypeAssertion (const ScopedLegacyInitTypeAssertion&) = delete;
    ScopedLegacyInitTypeAssertion& operator= (const ScopedLegacyInitTypeAssertion&) = delete;
};

// Makes amrex::Abort throw amrex::RuntimeError for the life of the object
class ScopedThrowOnAbort
{
public:
    ScopedThrowOnAbort () : m_saved(amrex::system::throw_exception)
    {
        amrex::system::throw_exception = true;
    }
    ~ScopedThrowOnAbort () { amrex::system::throw_exception = m_saved; }
    ScopedThrowOnAbort (const ScopedThrowOnAbort&) = delete;
    ScopedThrowOnAbort& operator= (const ScopedThrowOnAbort&) = delete;

private:
    bool m_saved;
};

// The plan a restart makes for a checkpoint, as ReadCheckpointFile makes it
cst::SurfaceTemperatureMigrationPlan
plan_restart (const TemporaryCheckpoint& checkpoint, const bool restart_reads_wrfinput)
{
    return cst::plan_surface_temperature_migration(
        cst::read_surface_temperature_provenance(checkpoint.root.string(), 0),
        restart_reads_wrfinput, checkpoint.root.string());
}

// Copy a single-component MultiFab, ghost cells included, to the host
std::vector<amrex::FArrayBox>
to_host (const amrex::MultiFab& mf)
{
    std::vector<amrex::FArrayBox> result;
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        result.emplace_back(mf[mfi].box(), 1, amrex::The_Pinned_Arena());
        amrex::Gpu::dtoh_memcpy(result.back().dataPtr(), mf[mfi].dataPtr(),
                                result.back().size() * sizeof(amrex::Real));
    }
    amrex::Gpu::streamSynchronize();
    return result;
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

    // Over water without TSK, where the surface layer takes the SST
    const amrex::iMultiFab water = make_water_mask(field);
    erf_checkpoint_surface_temperature::LegacyConversionFailure fail =
        erf_checkpoint_surface_temperature::potential_to_absolute_temperature(
            field, SurfaceTemperatureSource::sst,
            erf_checkpoint_surface_temperature::restart_surface_temperature_use(&water, false),
            pressure, RdoCp, amrex::Periodicity::NonPeriodic());
    EXPECT_EQ(fail.count, 0);
    EXPECT_EQ(fail.preserved, 0);

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

// The migration of a checkpoint's SST and TSK, planned from the checkpoint as
// ReadCheckpointFile plans it: from its marker, its job_info and the assertion.

TEST(CheckpointSurfaceTemperatureMigration, RecoversLegacyWrfInputArraysExactly)
{
    // Motivation (case 1): a WRFInput checkpoint without a marker holds the SST and
    // TSK as potential temperature, converted with the PSFC of its wrfinput files.  A
    // WRFInput restart reads that PSFC again and recovers them exactly.
    TemporaryCheckpoint checkpoint;
    make_checkpoint(checkpoint, 0, "WRFInput");
    const auto plan = plan_restart(checkpoint, true);
    EXPECT_EQ(plan.method, SurfaceTemperatureMigration::InputPressure);
    EXPECT_EQ(plan.provenance.source, LegacyInitType::WRFInput);
    EXPECT_EQ(plan.provenance.origin, ProvenanceOrigin::JobInfo);
    EXPECT_EQ(plan.provenance.stored, StoredSurfaceTemperature::Potential);
}

TEST(CheckpointSurfaceTemperatureMigration, RejectsWrfInputArraysOnARestartWithoutWrfInput)
{
    // Motivation (case 3): a restart need not have the inputs of the run that wrote the
    // checkpoint.  One that reads no wrfinput has no PSFC to recover WRFInput arrays
    // with, and converting them with the pressure of the restarted state instead would
    // silently use a pressure they were never converted with.
    for (const int marker_version : {0, 1}) {
        SCOPED_TRACE(marker_version);
        TemporaryCheckpoint checkpoint;
        make_checkpoint(checkpoint, marker_version, "WRFInput");
        const auto plan = plan_restart(checkpoint, false);
        EXPECT_EQ(plan.method, SurfaceTemperatureMigration::Reject);
        EXPECT_NE(plan.reason.find("restart with erf.init_type = WRFInput"), std::string::npos)
            << plan.reason;
    }
}

TEST(CheckpointSurfaceTemperatureMigration, MigratesMetgridArraysWithRestartPressureWhateverRestartReads)
{
    // Motivation (cases 3 and 5): a version 1 checkpoint of a Metgrid run holds potential
    // temperature converted with the surface pressure of its met_em input.  A restart
    // that reads wrfinput must not convert it with the PSFC of those unrelated files: the
    // provenance of the checkpoint, not the init type of the restart, picks the pressure.
    for (const bool restart_reads_wrfinput : {true, false}) {
        SCOPED_TRACE(restart_reads_wrfinput);
        TemporaryCheckpoint checkpoint;
        make_checkpoint(checkpoint, 1, "Metgrid");
        const auto plan = plan_restart(checkpoint, restart_reads_wrfinput);
        EXPECT_EQ(plan.method, SurfaceTemperatureMigration::RestartPressure);
        EXPECT_EQ(plan.provenance.source, LegacyInitType::Metgrid);
        EXPECT_EQ(plan.provenance.origin, ProvenanceOrigin::JobInfo);
    }
}

TEST(CheckpointSurfaceTemperatureMigration, TakesTheSourceOfVersionOneArraysFromJobInfo)
{
    // Motivation: the version 1 marker says the arrays hold potential temperature, not
    // which pressure they were converted with; the job_info says that.
    TemporaryCheckpoint checkpoint;
    make_checkpoint(checkpoint, 1, "WRFInput");
    const auto plan = plan_restart(checkpoint, true);
    EXPECT_EQ(plan.method, SurfaceTemperatureMigration::InputPressure);
    EXPECT_EQ(plan.provenance.source, LegacyInitType::WRFInput);
    EXPECT_EQ(plan.provenance.origin, ProvenanceOrigin::JobInfo);
}

TEST(CheckpointSurfaceTemperatureMigration, RejectsVersionOneArraysOfUnknownSourceUnlessAsserted)
{
    // Motivation: without a recorded source the pressure the arrays were converted with
    // is not known, and the restart must not guess it from its own inputs.
    TemporaryCheckpoint checkpoint;
    make_checkpoint(checkpoint, 1, "");
    for (const bool restart_reads_wrfinput : {true, false}) {
        SCOPED_TRACE(restart_reads_wrfinput);
        const auto plan = plan_restart(checkpoint, restart_reads_wrfinput);
        EXPECT_EQ(plan.method, SurfaceTemperatureMigration::Reject);
        EXPECT_NE(plan.reason.find("erf.legacy_surface_temperature_init_type"),
                  std::string::npos) << plan.reason;
    }
    {
        ScopedLegacyInitTypeAssertion assertion("WRFInput");
        const auto plan = plan_restart(checkpoint, true);
        EXPECT_EQ(plan.method, SurfaceTemperatureMigration::InputPressure);
        EXPECT_EQ(plan.provenance.origin, ProvenanceOrigin::Assertion);
    }
    {
        ScopedLegacyInitTypeAssertion assertion("metgrid");
        const auto plan = plan_restart(checkpoint, true);
        EXPECT_EQ(plan.method, SurfaceTemperatureMigration::RestartPressure);
        EXPECT_EQ(plan.provenance.origin, ProvenanceOrigin::Assertion);
    }
}

TEST(CheckpointSurfaceTemperatureMigration, UsesTheAssertionForMarkerlessArraysOfUnknownSource)
{
    // Motivation (case 4): a markerless checkpoint whose job_info records no init type is
    // migrated as the explicit assertion says, and rejected without one.
    TemporaryCheckpoint checkpoint;
    make_checkpoint(checkpoint, 0, "");
    {
        const auto plan = plan_restart(checkpoint, true);
        EXPECT_EQ(plan.method, SurfaceTemperatureMigration::Reject);
        EXPECT_NE(plan.reason.find("erf.legacy_surface_temperature_init_type"),
                  std::string::npos) << plan.reason;
    }
    {
        ScopedLegacyInitTypeAssertion assertion("WRFInput");
        const auto plan = plan_restart(checkpoint, true);
        EXPECT_EQ(plan.method, SurfaceTemperatureMigration::InputPressure);
        EXPECT_EQ(plan.provenance.source, LegacyInitType::WRFInput);
        EXPECT_EQ(plan.provenance.origin, ProvenanceOrigin::Assertion);
    }
    {
        // Before the marker, Metgrid stored absolute temperature
        ScopedLegacyInitTypeAssertion assertion("Metgrid");
        EXPECT_EQ(plan_restart(checkpoint, true).method, SurfaceTemperatureMigration::None);
    }
}

TEST(CheckpointSurfaceTemperatureMigration, JobInfoOutranksTheAssertion)
{
    // Motivation: the assertion stands in for a provenance the checkpoint does not
    // record; it cannot override one it does.
    TemporaryCheckpoint checkpoint;
    make_checkpoint(checkpoint, 0, "Metgrid");
    ScopedLegacyInitTypeAssertion assertion("WRFInput");
    const auto plan = plan_restart(checkpoint, true);
    EXPECT_EQ(plan.method, SurfaceTemperatureMigration::None);
    EXPECT_EQ(plan.provenance.source, LegacyInitType::Metgrid);
    EXPECT_EQ(plan.provenance.origin, ProvenanceOrigin::JobInfo);
    EXPECT_EQ(plan.provenance.asserted, LegacyInitType::WRFInput);
}

TEST(CheckpointSurfaceTemperatureMigration, RestoresVersionTwoArraysAsTheyAre)
{
    // Motivation (case 6): version 2 arrays hold absolute temperature.  Whatever wrote
    // them and whatever the restart reads, they are restored without a conversion.
    for (const std::string& init_type : {std::string("WRFInput"), std::string("Metgrid"),
                                         std::string()}) {
        for (const bool restart_reads_wrfinput : {true, false}) {
            SCOPED_TRACE(init_type + (restart_reads_wrfinput ? ", WRFInput restart"
                                                             : ", other restart"));
            TemporaryCheckpoint checkpoint;
            make_checkpoint(checkpoint, cst::contract_version, init_type);
            const auto plan = plan_restart(checkpoint, restart_reads_wrfinput);
            EXPECT_EQ(plan.method, SurfaceTemperatureMigration::None);
            EXPECT_EQ(plan.provenance.stored, StoredSurfaceTemperature::Absolute);
        }
    }
}

TEST(CheckpointSurfaceTemperatureMigration, RejectsAnInvalidAssertion)
{
    TemporaryCheckpoint checkpoint;
    make_checkpoint(checkpoint, 0, "");
    ScopedLegacyInitTypeAssertion assertion("Ideal");
    ScopedThrowOnAbort throw_on_abort;
    try {
        plan_restart(checkpoint, true);
        ADD_FAILURE() << "expected an abort";
    } catch (const amrex::RuntimeError& e) {
        EXPECT_NE(std::string(e.what()).find("accepted values are WRFInput and Metgrid"),
                  std::string::npos) << e.what();
    }
}

TEST(CheckpointSurfaceTemperatureMigration, RecoversPressureFromWhereEachLevelReadItsSurfaceTemperature)
{
    // Motivation (case 2): level 0 read its SST and TSK from its first init file however
    // many erf.nc_init_file_0 lists, so its PSFC is recovered from that file.  Before,
    // several level-0 files left the PSFC unread, and the arrays were converted with
    // whatever the allocated mf_PSFC held.
    using Source = LegacySurfacePressureSource;
    // Arguments: level, its init files, level 0 reads wrflowinp, level has its own
    // wrflowinp, coarser level recovered
    EXPECT_EQ(cst::legacy_wrfinput_pressure_source(0, 1, false, false, false), Source::OwnInitFile);
    EXPECT_EQ(cst::legacy_wrfinput_pressure_source(0, 2, false, false, false), Source::OwnInitFile);
    EXPECT_EQ(cst::legacy_wrfinput_pressure_source(0, 3, true, false, false), Source::OwnInitFile);
    EXPECT_EQ(cst::legacy_wrfinput_pressure_source(0, 0, false, false, false), Source::Unavailable);

    // A finer level read them from its own file, unless it took them from the coarser one
    EXPECT_EQ(cst::legacy_wrfinput_pressure_source(1, 1, false, false, true), Source::OwnInitFile);
    EXPECT_EQ(cst::legacy_wrfinput_pressure_source(1, 1, true, true, true), Source::OwnInitFile);
    EXPECT_EQ(cst::legacy_wrfinput_pressure_source(1, 1, true, false, true), Source::CoarserLevel);
    EXPECT_EQ(cst::legacy_wrfinput_pressure_source(2, 0, false, false, true), Source::CoarserLevel);
    EXPECT_EQ(cst::legacy_wrfinput_pressure_source(2, 0, false, false, false), Source::Unavailable);
    EXPECT_EQ(cst::legacy_wrfinput_pressure_source(1, 1, true, false, false), Source::Unavailable);

    // A finer level with several files of its own read them from all of them, which is
    // not recovered, unless it took them from the coarser level
    EXPECT_EQ(cst::legacy_wrfinput_pressure_source(1, 2, false, false, true), Source::Unavailable);
    EXPECT_EQ(cst::legacy_wrfinput_pressure_source(1, 2, true, false, true), Source::CoarserLevel);
}

TEST(CheckpointSurfaceTemperatureMigration, CountsCellsWithoutUsablePressureOrPotentialTemperature)
{
    // Motivation: mf_PSFC is allocated but not filled on restart, so an allocated pressure
    // proves nothing.  A valid cell where the surface layer takes the field (water, without
    // TSK, for the SST) whose pressure was not recovered (zero), or whose pressure or
    // potential temperature is not finite, must be left unconverted and counted, to abort
    // the restart; a ghost cell is not counted.
    const amrex::Box box(amrex::IntVect(0,0,0), amrex::IntVect(3,2,0));
    const amrex::BoxArray ba(box);
    const amrex::DistributionMapping dm(ba);
    const amrex::IntVect ng(1,1,0);
    amrex::MultiFab field(ba, dm, 1, ng);
    amrex::MultiFab pressure(ba, dm, 1, ng);
    const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();
    const amrex::Real theta = 300.0;
    const amrex::Real p_s = 90000.0;

    set_from_host(pressure, [&] (int i, int j, int) {
        if (i == 1 && j == 1) { return amrex::Real(0.0); }
        if (i == 2 && j == 0) { return nan; }
        return p_s;
    });
    set_from_host(field, [&] (int i, int j, int) {
        return ((i == 3 && j == 2) || (i == -1 && j == 0)) ? nan : theta;
    });

    const amrex::iMultiFab water = make_water_mask(field);
    const auto fail = cst::potential_to_absolute_temperature(
        field, SurfaceTemperatureSource::sst, cst::restart_surface_temperature_use(&water, false),
        pressure, RdoCp, amrex::Periodicity::NonPeriodic());
    EXPECT_EQ(fail.count, 3);
    EXPECT_EQ(fail.preserved, 0);
    EXPECT_EQ(fail.claimed, 1);
    EXPECT_EQ(fail.land, 0);
    EXPECT_EQ(fail.used, 1);

    const auto host = to_host(field);
    ASSERT_EQ(host.size(), 1u);
    const auto values = host[0].const_array();
    const amrex::Real converted = theta * std::pow(p_s / p_0, RdoCp);
    const amrex::Real tol = amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon() *
                            converted;
    EXPECT_EQ(values(1,1,0), theta);
    EXPECT_EQ(values(2,0,0), theta);
    EXPECT_TRUE(std::isnan(values(3,2,0)));
    EXPECT_TRUE(std::isnan(values(-1,0,0)));
    EXPECT_NEAR(values(0,0,0), converted, tol);
    EXPECT_NEAR(values(3,1,0), converted, tol);
    EXPECT_NEAR(values(-1,1,0), converted, tol);
}

TEST(CheckpointSurfaceTemperatureMigration, AbortsDescribingTheFirstCellThatCouldNotBeConverted)
{
    ScopedThrowOnAbort throw_on_abort;
    EXPECT_NO_THROW(cst::abort_on_legacy_conversion_failure(cst::LegacyConversionFailure{},
                                                             SurfaceTemperatureSource::sst,
                                                             "converting the SST"));

    cst::LegacyConversionFailure fail;
    fail.count = 2;
    fail.claimed = 1;
    fail.i = 3;
    fail.j = 1;
    fail.land = 0;
    fail.used = 1;
    fail.theta = 301.5;
    fail.pressure = 0.0;
    std::string message;
    testing::internal::CaptureStdout();
    try {
        cst::abort_on_legacy_conversion_failure(fail, SurfaceTemperatureSource::sst,
                                                "converting the SST at level 0");
    } catch (const amrex::RuntimeError& e) {
        message = e.what();
    }
    const std::string report = testing::internal::GetCapturedStdout();
    EXPECT_NE(message.find("2 cells could not be converted"), std::string::npos) << message;
    EXPECT_NE(message.find("where the surface layer takes the SST needs a finite, positive "
                           "potential temperature"), std::string::npos) << message;
    EXPECT_NE(report.find("converting the SST at level 0"), std::string::npos) << report;
    EXPECT_NE(report.find("(3,1), a water cell where the surface layer takes the SST"),
              std::string::npos) << report;
    EXPECT_NE(report.find("surface pressure = 0 Pa"), std::string::npos) << report;
    EXPECT_EQ(report.find("placeholders"), std::string::npos) << report;
}

TEST(CheckpointSurfaceTemperatureMigration, ApproximateMigrationKeepsTheRestartPressureBias)
{
    // Motivation (case 5): converting a Metgrid potential temperature with the surface
    // pressure of the restarted state, not that of the input, keeps the surface
    // potential temperature of the restart, not the SST of the input: it gives
    // SST (p_restart/p_input)^(R_d/c_p), as documented.
    const amrex::Box box(amrex::IntVect(0,0,0), amrex::IntVect(1,1,0));
    const amrex::BoxArray ba(box);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab field(ba, dm, 1, 0);
    amrex::MultiFab pressure(ba, dm, 1, 0);
    const amrex::Real sst = 290.0;
    const amrex::Real p_input = 101325.0;
    const amrex::Real p_restart = 95000.0;
    field.setVal(sst * std::pow(p_0 / p_input, RdoCp));
    pressure.setVal(p_restart);

    const amrex::iMultiFab water = make_water_mask(field);
    const auto fail = cst::potential_to_absolute_temperature(
        field, SurfaceTemperatureSource::sst, cst::restart_surface_temperature_use(&water, false),
        pressure, RdoCp, amrex::Periodicity::NonPeriodic());
    EXPECT_EQ(fail.count, 0);

    const auto host = to_host(field);
    const amrex::Real migrated = host[0].const_array()(0,0,0);
    EXPECT_NEAR(migrated, sst * std::pow(p_restart / p_input, RdoCp),
                amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon() * sst);
    EXPECT_GT(std::abs(migrated - sst), amrex::Real(1.0));
}

// Where the surface layer takes each field.  A restart judges a stored SST or TSK by whether
// the surface layer takes that field at its cell, by the rule the surface layer follows.

TEST(CheckpointSurfaceTemperatureMigration, TakesEachFieldWhereTheSurfaceLayerDoes)
{
    // Motivation: make_SurfaceLayer_at_level ignores the SST of a level with TSK, which then
    // supplies every cell, water too; without TSK, water takes the SST and land the
    // prescribed temperature, whatever most.ignore_sst says.  A restart, which reads no
    // most.ignore_sst, must reach the same selection.
    namespace est = erf_surface_temperature;
    EXPECT_FALSE(est::sst_ignored(false, false));
    EXPECT_TRUE(est::sst_ignored(true, false));
    EXPECT_TRUE(est::sst_ignored(false, true));
    EXPECT_TRUE(est::sst_ignored(true, true));

    for (const bool ignore_sst : {false, true}) {
        EXPECT_EQ(est::surface_temperature_source(false, false, ignore_sst),
                  SurfaceTemperatureSource::sst);
        EXPECT_EQ(est::surface_temperature_source(true, false, ignore_sst),
                  SurfaceTemperatureSource::prescribed);
    }
    EXPECT_EQ(est::surface_temperature_source(false, true, true), SurfaceTemperatureSource::tsk);
    EXPECT_EQ(est::surface_temperature_source(true, true, true), SurfaceTemperatureSource::tsk);

    const auto without_tsk = cst::restart_surface_temperature_use(nullptr, false);
    EXPECT_TRUE(without_tsk.lmask == nullptr);
    EXPECT_FALSE(without_tsk.have_tsk);
    EXPECT_FALSE(without_tsk.ignore_sst);
    const auto with_tsk = cst::restart_surface_temperature_use(nullptr, true);
    EXPECT_TRUE(with_tsk.have_tsk);
    EXPECT_TRUE(with_tsk.ignore_sst);
}

TEST(CheckpointSurfaceTemperatureMigration, KeepsUnusedLandSstPlaceholdersWithoutCountingThem)
{
    // Motivation: wrfinput files hold a placeholder SST over land, often zero, which a
    // version 1 or markerless WRFInput checkpoint holds as it is.  Without TSK the surface
    // layer takes the prescribed temperature over land, never the SST, so the restart must
    // not abort on such a placeholder, nor make it look like a temperature: it is kept as
    // it is.  The SST over water, and a positive SST over land, are converted, so that no
    // potential temperature is left in a field held as absolute temperature.
    const amrex::Box box(amrex::IntVect(0,0,0), amrex::IntVect(3,2,0));
    const amrex::BoxArray ba(box);
    const amrex::DistributionMapping dm(ba);
    const amrex::IntVect ng(1,1,0);
    amrex::MultiFab field(ba, dm, 1, ng);
    amrex::MultiFab pressure(ba, dm, 1, ng);
    const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();
    const amrex::Real theta = 300.0;
    const amrex::Real land_theta = 290.0;
    const amrex::Real p_s = 90000.0;
    pressure.setVal(p_s);
    const auto land = make_land_mask(field, [] (int i, int, int) { return (i < 2) ? 1 : 0; });
    set_from_host(field, [&] (int i, int j, int) {
        if (i >= 2) { return theta; }
        if (i == 1 && j == 0) { return nan; }
        if (i == 0 && j == 1) { return amrex::Real(-1.0); }
        if (i == 1 && j == 2) { return land_theta; }
        return amrex::Real(0.0);
    });

    const auto fail = cst::potential_to_absolute_temperature(
        field, SurfaceTemperatureSource::sst, cst::restart_surface_temperature_use(&land, false),
        pressure, RdoCp, amrex::Periodicity::NonPeriodic());
    EXPECT_EQ(fail.count, 0);
    EXPECT_EQ(fail.preserved, 5);
    EXPECT_EQ(fail.claimed, 0);

    const auto host = to_host(field);
    ASSERT_EQ(host.size(), 1u);
    const auto values = host[0].const_array();
    const amrex::Real factor = std::pow(p_s / p_0, RdoCp);
    const amrex::Real tol = amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon() *
                            theta;
    EXPECT_EQ(values(0,0,0), amrex::Real(0.0));
    EXPECT_TRUE(std::isnan(values(1,0,0)));
    EXPECT_EQ(values(0,1,0), amrex::Real(-1.0));
    EXPECT_EQ(values(1,1,0), amrex::Real(0.0));
    EXPECT_EQ(values(0,2,0), amrex::Real(0.0));
    EXPECT_NEAR(values(1,2,0), land_theta * factor, tol);
    for (int j = 0; j <= 2; ++j) {
        for (int i = 2; i <= 3; ++i) {
            EXPECT_NEAR(values(i,j,0), theta * factor, tol) << "at (" << i << "," << j << ")";
        }
    }
}

TEST(CheckpointSurfaceTemperatureMigration, CountsAnInvalidSstWhereTheSurfaceLayerTakesIt)
{
    // Motivation: without TSK the surface layer takes the SST over water, so a zero there is
    // a missing value the run would use: the restart must count it, to abort, and describe
    // it, while it keeps the land placeholders beside it.
    const amrex::Box box(amrex::IntVect(0,0,0), amrex::IntVect(3,2,0));
    const amrex::BoxArray ba(box);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab field(ba, dm, 1, amrex::IntVect(1,1,0));
    amrex::MultiFab pressure(ba, dm, 1, amrex::IntVect(1,1,0));
    const amrex::Real p_s = 90000.0;
    pressure.setVal(p_s);
    const auto land = make_land_mask(field, [] (int i, int, int) { return (i < 2) ? 1 : 0; });
    set_from_host(field, [&] (int i, int j, int) {
        return (i < 2 || (i == 3 && j == 1)) ? amrex::Real(0.0) : amrex::Real(300.0);
    });

    const auto fail = cst::potential_to_absolute_temperature(
        field, SurfaceTemperatureSource::sst, cst::restart_surface_temperature_use(&land, false),
        pressure, RdoCp, amrex::Periodicity::NonPeriodic());
    EXPECT_EQ(fail.count, 1);
    EXPECT_EQ(fail.preserved, 6);
    EXPECT_EQ(fail.claimed, 1);
    EXPECT_EQ(fail.i, 3);
    EXPECT_EQ(fail.j, 1);
    EXPECT_EQ(fail.land, 0);
    EXPECT_EQ(fail.used, 1);
    EXPECT_EQ(fail.theta, amrex::Real(0.0));
    EXPECT_EQ(fail.pressure, p_s);
    EXPECT_EQ(to_host(field)[0].const_array()(3,1,0), amrex::Real(0.0));
}

TEST(CheckpointSurfaceTemperatureMigration, IgnoresTheSstAndTakesTheTskEverywhereWithTsk)
{
    // Motivation: a level with TSK takes it at every cell, water too, and ignores its SST.
    // A zero SST over water is then a placeholder the restart keeps, but a zero TSK, over
    // land or water, is a value the run would use, and must abort the restart.
    const amrex::Box box(amrex::IntVect(0,0,0), amrex::IntVect(3,2,0));
    const amrex::BoxArray ba(box);
    const amrex::DistributionMapping dm(ba);
    const amrex::IntVect ng(1,1,0);
    amrex::MultiFab sst(ba, dm, 1, ng);
    amrex::MultiFab tsk(ba, dm, 1, ng);
    amrex::MultiFab pressure(ba, dm, 1, ng);
    pressure.setVal(90000.0);
    const auto land = make_land_mask(sst, [] (int i, int, int) { return (i < 2) ? 1 : 0; });
    const auto use = cst::restart_surface_temperature_use(&land, true);

    set_from_host(sst, [&] (int i, int j, int) {
        return (i < 2 || (i == 3 && j == 1)) ? amrex::Real(0.0) : amrex::Real(300.0);
    });
    const auto sst_fail = cst::potential_to_absolute_temperature(
        sst, SurfaceTemperatureSource::sst, use, pressure, RdoCp,
        amrex::Periodicity::NonPeriodic());
    EXPECT_EQ(sst_fail.count, 0);
    EXPECT_EQ(sst_fail.preserved, 7);

    for (const amrex::IntVect& bad : {amrex::IntVect(1,2,0), amrex::IntVect(2,0,0)}) {
        SCOPED_TRACE(bad);
        set_from_host(tsk, [&] (int i, int j, int) {
            return (i == bad[0] && j == bad[1]) ? amrex::Real(0.0) : amrex::Real(301.0);
        });
        const auto tsk_fail = cst::potential_to_absolute_temperature(
            tsk, SurfaceTemperatureSource::tsk, use, pressure, RdoCp,
            amrex::Periodicity::NonPeriodic());
        EXPECT_EQ(tsk_fail.count, 1);
        EXPECT_EQ(tsk_fail.preserved, 0);
        EXPECT_EQ(tsk_fail.i, bad[0]);
        EXPECT_EQ(tsk_fail.j, bad[1]);
        EXPECT_EQ(tsk_fail.land, (bad[0] < 2) ? 1 : 0);
        EXPECT_EQ(tsk_fail.used, 1);
    }
}

TEST(CheckpointSurfaceTemperatureMigration, CountsAPositiveValueWithoutPressureWhereverItIs)
{
    // Motivation: a positive value is a potential temperature, whether or not the surface
    // layer takes the field at its cell.  Left unconverted, it would sit in a field held as
    // absolute temperature, to be used should a regrid or an interpolation bring it where
    // the field is taken.  Converting it takes a finite, positive surface pressure, and a
    // cell without one aborts the restart; a placeholder needs none.
    const amrex::Box box(amrex::IntVect(0,0,0), amrex::IntVect(3,2,0));
    const amrex::BoxArray ba(box);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab field(ba, dm, 1, amrex::IntVect(1,1,0));
    amrex::MultiFab pressure(ba, dm, 1, amrex::IntVect(1,1,0));
    set_from_host(pressure, [&] (int i, int j, int) {
        return ((i == 0 && j == 0) || (i == 1 && j == 1)) ? amrex::Real(0.0)
                                                          : amrex::Real(90000.0);
    });
    const auto land = make_land_mask(field, [] (int i, int, int) { return (i < 2) ? 1 : 0; });
    set_from_host(field, [&] (int i, int j, int) {
        if (i >= 2 || (i == 1 && j == 1)) { return amrex::Real(290.0); }
        return amrex::Real(0.0);
    });

    const auto fail = cst::potential_to_absolute_temperature(
        field, SurfaceTemperatureSource::sst, cst::restart_surface_temperature_use(&land, false),
        pressure, RdoCp, amrex::Periodicity::NonPeriodic());
    EXPECT_EQ(fail.count, 1);
    EXPECT_EQ(fail.preserved, 5);
    EXPECT_EQ(fail.i, 1);
    EXPECT_EQ(fail.j, 1);
    EXPECT_EQ(fail.land, 1);
    EXPECT_EQ(fail.used, 0);
    EXPECT_EQ(fail.pressure, amrex::Real(0.0));
    const auto values = to_host(field)[0].const_array();
    EXPECT_EQ(values(1,1,0), amrex::Real(290.0));
    EXPECT_EQ(values(0,0,0), amrex::Real(0.0));
}

TEST(CheckpointSurfaceTemperatureMigration, TreatsEveryCellAsLandWithoutALandMask)
{
    // Motivation: without a land mask the surface layer takes every cell as land, and so
    // does the restart.  Without TSK, an SST placeholder anywhere is then kept; with TSK, a
    // zero TSK anywhere aborts.
    const amrex::Box box(amrex::IntVect(0,0,0), amrex::IntVect(3,2,0));
    const amrex::BoxArray ba(box);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab field(ba, dm, 1, amrex::IntVect(1,1,0));
    amrex::MultiFab pressure(ba, dm, 1, amrex::IntVect(1,1,0));
    pressure.setVal(90000.0);

    field.setVal(0.0);
    const auto sst_fail = cst::potential_to_absolute_temperature(
        field, SurfaceTemperatureSource::sst, cst::restart_surface_temperature_use(nullptr, false),
        pressure, RdoCp, amrex::Periodicity::NonPeriodic());
    EXPECT_EQ(sst_fail.count, 0);
    EXPECT_EQ(sst_fail.preserved, 12);

    field.setVal(0.0);
    const auto tsk_fail = cst::potential_to_absolute_temperature(
        field, SurfaceTemperatureSource::tsk, cst::restart_surface_temperature_use(nullptr, true),
        pressure, RdoCp, amrex::Periodicity::NonPeriodic());
    EXPECT_EQ(tsk_fail.count, 12);
    EXPECT_EQ(tsk_fail.preserved, 0);
    EXPECT_EQ(tsk_fail.land, 1);
    EXPECT_EQ(tsk_fail.used, 1);
}

TEST(CheckpointSurfaceTemperatureMigration, ReportsKeptPlaceholdersAndAbortsOnlyOnFailures)
{
    // Motivation: the kept placeholders are reported, not counted as failures; a failure
    // beside them still aborts, describing whether the surface layer takes the field there.
    ScopedThrowOnAbort throw_on_abort;
    cst::LegacyConversionFailure fail;
    fail.preserved = 6;
    testing::internal::CaptureStdout();
    EXPECT_NO_THROW(cst::abort_on_legacy_conversion_failure(fail, SurfaceTemperatureSource::sst,
                                                             "converting the SST at level 0"));
    std::string report = testing::internal::GetCapturedStdout();
    EXPECT_NE(report.find("converting the SST at level 0: kept 6 cells where the surface "
                          "layer does not take the SST as they are, placeholders"),
              std::string::npos) << report;

    fail.count = 1;
    fail.claimed = 1;
    fail.i = 1;
    fail.j = 1;
    fail.land = 1;
    fail.used = 0;
    fail.theta = 290.0;
    fail.pressure = 0.0;
    std::string message;
    testing::internal::CaptureStdout();
    try {
        cst::abort_on_legacy_conversion_failure(fail, SurfaceTemperatureSource::sst,
                                                "converting the SST at level 0");
    } catch (const amrex::RuntimeError& e) {
        message = e.what();
    }
    report = testing::internal::GetCapturedStdout();
    EXPECT_NE(message.find("1 cells could not be converted"), std::string::npos) << message;
    EXPECT_NE(report.find("kept 6 cells"), std::string::npos) << report;
    EXPECT_NE(report.find("(1,1), a land cell where the surface layer does not take the SST, "
                          "has theta = 290, surface pressure = 0 Pa"),
              std::string::npos) << report;
}

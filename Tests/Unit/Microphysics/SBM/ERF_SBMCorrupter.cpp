#include <AMReX.H>
#include <AMReX_VisMF.H>

#include "ERF_IndexDefines.H"
#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMBulkProjection.H"
#include "ERF_SBMRestart.H"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

erf_sbm::SBMLayout make_fixture_layout()
{
    erf_sbm::SpectralPopulationSpec population;
    population.population_id = 0;
    population.semantic_id = "liquid_mass";
    population.phase = erf_sbm::PopulationPhase::Liquid;
    population.moment_mode = erf_sbm::MomentMode::TwoMoment;
    population.grid.coordinate_kind = erf_sbm::CoordinateKind::Mass;
    population.grid.coordinate_units = "kg";
    population.grid.edges = {
        amrex::Real(1.0e-18), amrex::Real(2.0e-18), amrex::Real(4.0e-18),
        amrex::Real(8.0e-18), amrex::Real(1.6e-17)};
    population.grid.pivots = {
        amrex::Real(1.5e-18), amrex::Real(3.0e-18), amrex::Real(6.0e-18),
        amrex::Real(1.2e-17)};

    erf_sbm::SBMLayoutSpec spec;
    spec.populations.push_back(std::move(population));
    spec.liquid_projection = {0, 2};
    return erf_sbm::SBMLayout(std::move(spec));
}

void require(const bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

// Make one 2M number moment inadmissible while keeping the compact projection.
void corrupt_number_moment (amrex::MultiFab& spectrum, const amrex::MultiFab& core,
                            const std::string& spectrum_prefix)
{
    const auto layout = make_fixture_layout();
    require(spectrum.nComp() == layout.ncomp(),
            "checkpoint spectrum component count differs from the test layout");

    const erf_sbm::SBMBulkProjection projection(layout);
    std::string diagnostic;
    require(erf_sbm::authoritative_state_admissible(spectrum, layout, 0, &diagnostic),
            "valid input checkpoint was rejected: " + diagnostic);
    require(erf_sbm::restart_projection_matches(
                spectrum, core, projection, RhoQ2_comp, RhoQ3_comp),
            "input checkpoint compact state does not match its spectrum");

    const int number = layout.populations()[0].number_offset;
    spectrum.setVal(amrex::Real(1.0e10), number, 1, 0);
    require(!erf_sbm::authoritative_state_admissible(spectrum, layout, 0, &diagnostic),
            "corrupted number moment remained admissible");
    require(diagnostic.find("constraint=canonical-two-moment-bin-state") != std::string::npos,
            "number corruption did not violate the expected upper endpoint: " + diagnostic);
    require(erf_sbm::restart_projection_matches(
                spectrum, core, projection, RhoQ2_comp, RhoQ3_comp),
            "number-only corruption unexpectedly changed the compact projection");

    amrex::VisMF::RemoveFiles(spectrum_prefix);
    amrex::VisMF::Write(spectrum, spectrum_prefix);

    amrex::MultiFab reread;
    amrex::VisMF::Read(reread, spectrum_prefix);
    require(!erf_sbm::authoritative_state_admissible(reread, layout, 0, &diagnostic),
            "written number corruption was not preserved: " + diagnostic);
    require(diagnostic.find("constraint=canonical-two-moment-bin-state") != std::string::npos,
            "written corruption failed for an unexpected constraint: " + diagnostic);
    require(erf_sbm::restart_projection_matches(
                reread, core, projection, RhoQ2_comp, RhoQ3_comp),
            "persisted compact state is not projection-consistent after corruption");
    amrex::Print() << "Prepared finite number-moment restart corruption; compact projection remains consistent.\n";
}

} // namespace

int main(int argc, char** argv)
{
    amrex::Initialize(argc, argv);
    int result = 0;
    try {
        amrex::ignore_unused(argc, argv);
        const char* checkpoint_env = std::getenv("SBM_CHECKPOINT_DIR");
        require(checkpoint_env != nullptr && checkpoint_env[0] != '\0',
                "SBM_CHECKPOINT_DIR must name the checkpoint to corrupt");
        const std::filesystem::path checkpoint(checkpoint_env);
        // Optional: which level to corrupt, and how. The layout-free modes
        // serve any fixture; number-moment needs the 2M fixture layout.
        const char* level_env = std::getenv("SBM_CHECKPOINT_LEVEL");
        const char* mode_env = std::getenv("SBM_CORRUPTION");
        const int level_index =
            (level_env != nullptr && level_env[0] != '\0') ? std::stoi(level_env) : 0;
        const std::string mode =
            (mode_env != nullptr && mode_env[0] != '\0') ? mode_env : "number-moment";
        const auto level = checkpoint / ("Level_" + std::to_string(level_index));
        const std::string spectrum_prefix = (level / "SBMSpectrum").string();
        const std::string core_prefix = (level / "Cell").string();

        amrex::MultiFab spectrum;
        amrex::VisMF::Read(spectrum, spectrum_prefix);
        amrex::MultiFab core;
        amrex::VisMF::Read(core, core_prefix);
        if (mode == "negative-mass") {
            spectrum.setVal(amrex::Real(-1.0), 0, 1, 0);
            amrex::VisMF::RemoveFiles(spectrum_prefix);
            amrex::VisMF::Write(spectrum, spectrum_prefix);
            amrex::Print() << "Prepared negative spectral mass on level " << level_index << ".\n";
        } else if (mode == "compact-qc") {
            require(core.norm0(RhoQ2_comp) > amrex::Real(0.0),
                    "compact qc is zero, so scaling it cannot create a mismatch");
            core.mult(amrex::Real(2.0), RhoQ2_comp, 1, 0);
            amrex::VisMF::RemoveFiles(core_prefix);
            amrex::VisMF::Write(core, core_prefix);
            amrex::Print() << "Prepared compact qc mismatch on level " << level_index << ".\n";
        } else {
            require(mode == "number-moment", "unknown SBM_CORRUPTION mode: " + mode);
            corrupt_number_moment(spectrum, core, spectrum_prefix);
        }
    } catch (const std::exception& error) {
        amrex::Print() << "SBM checkpoint corruption helper failed: " << error.what() << '\n';
        result = 1;
    }
    amrex::Finalize();
    return result;
}

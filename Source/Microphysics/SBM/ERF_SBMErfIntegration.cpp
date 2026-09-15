#include "ERF.H"

#include "ERF_SBMBulkProjection.H"
#include "ERF_SBMContracts.H"
#include "ERF_SBMTransferClosure.H"
#include "ERF_SBMTransportPrototype.H"

#include <AMReX_MFParallelFor.H>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace amrex;

namespace {

amrex::Real sbm_total_mass(const amrex::MultiFab& state, const int ncomp,
                           const amrex::Real cell_volume)
{
    amrex::Real total = amrex::Real(0.0);
    for (int comp = 0; comp < ncomp; ++comp) total += state.sum(comp);
    return cell_volume * total;
}

amrex::Real sbm_max_projection_error(const erf_sbm::SBMLayout& layout,
                                     const amrex::MultiFab& spectral,
                                     const amrex::MultiFab& core)
{
    const auto& projection = layout.liquid_projection();
    const auto liquid = std::find_if(layout.populations().begin(), layout.populations().end(),
        [&](const erf_sbm::PopulationLayout& population) {
            return population.population_id == projection.population_id;
        });
    const int first = liquid->mass_offset;
    const int split = projection.cloud_rain_split;
    const int nbins = liquid->grid.nbins();
    amrex::MultiFab error(core.boxArray(), core.DistributionMap(), 1, 0);
    error.setVal(amrex::Real(0.0));
    for (amrex::MFIter mfi(error); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.validbox();
        const auto spectrum = spectral.const_array(mfi);
        const auto compact = core.const_array(mfi);
        const auto result = error.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            amrex::Real qc = amrex::Real(0.0);
            amrex::Real qr = amrex::Real(0.0);
            for (int b = 0; b < split; ++b) qc += spectrum(i,j,k,first+b);
            for (int b = split; b < nbins; ++b) qr += spectrum(i,j,k,first+b);
            result(i,j,k) = amrex::max(
                amrex::Math::abs(compact(i,j,k,RhoQ2_comp) - qc),
                amrex::Math::abs(compact(i,j,k,RhoQ3_comp) - qr));
        });
    }
    return error.max(0);
}

amrex::Real sbm_max_change(amrex::MultiFab& scratch,
                           const amrex::MultiFab& current,
                           const amrex::MultiFab& baseline,
                           const int ncomp)
{
    amrex::Real maximum = amrex::Real(0.0);
    for (int comp = 0; comp < ncomp; ++comp) {
        amrex::MultiFab::Copy(scratch, current, comp, 0, 1, 0);
        amrex::MultiFab::Subtract(scratch, baseline, comp, 0, 1, 0);
        maximum = amrex::max(maximum, scratch.norm0(0));
    }
    return maximum;
}

std::size_t sbm_cell_bytes(const amrex::MultiFab& state)
{
    return ::erf_auxiliary::allocated_payload_bytes(state);
}

void write_sbm_diagnostic(const std::string& path,
                          const bool anelastic,
                          const int nbins,
                          const int step_count,
                          const amrex::Real initial_mass,
                          const amrex::Real final_mass,
                          const amrex::Real initial_variation,
                          const amrex::Real transport_change,
                          const amrex::Real projection_error,
                          const amrex::Real face_projection_error,
                          const erf_sbm::AcceptedTransferClosure& closure,
                          const amrex::Real compact_mass,
                          const amrex::Real cell_volume,
                          const std::size_t cell_state_bytes,
                          const std::size_t face_transfer_bytes,
                          const std::size_t total_auxiliary_bytes)
{
    const amrex::Real mass_error = amrex::Math::abs(final_mass - initial_mass);
    const amrex::Real compact_mass_error = amrex::Math::abs(compact_mass - final_mass);
    const amrex::Real mass_tolerance = amrex::Real(1.0e-10) *
        amrex::max(amrex::Real(1.0), amrex::Math::abs(initial_mass));
    const amrex::Real projection_tolerance = amrex::Real(1.0e-12) *
        amrex::max(amrex::Real(1.0), amrex::Math::abs(compact_mass));
    const amrex::Real face_tolerance = amrex::Real(1.0e-12) *
        amrex::max(amrex::Real(1.0), amrex::Math::abs(final_mass));
    const bool passed = initial_variation > amrex::Real(0.0) &&
                        transport_change > amrex::Real(0.0) &&
                        mass_error <= mass_tolerance &&
                        compact_mass_error <= mass_tolerance &&
                        projection_error <= projection_tolerance &&
                        face_projection_error <= face_tolerance &&
                        closure.passes();

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream output(path);
        if (!output) amrex::Error("unable to write SBM P1 diagnostic: " + path);
        output << std::setprecision(17)
               << "format=erf-sbm-p1-diagnostic-v2\n"
               << "method=" << (anelastic ? "anelastic" : "compressible") << '\n'
               << "nbins=" << nbins << '\n'
               << "step_count=" << step_count << '\n'
               << "cell_volume=" << cell_volume << '\n'
               << "initial_mass=" << initial_mass << '\n'
               << "final_mass=" << final_mass << '\n'
               << "mass_error=" << mass_error << '\n'
               << "mass_tolerance=" << mass_tolerance << '\n'
               << "compact_mass=" << compact_mass << '\n'
               << "compact_mass_error=" << compact_mass_error << '\n'
               << "initial_variation=" << initial_variation << '\n'
               << "transport_change=" << transport_change << '\n'
               << "projection_error=" << projection_error << '\n'
               << "projection_tolerance=" << projection_tolerance << '\n'
               << "face_projection_error=" << face_projection_error << '\n'
               << "face_tolerance=" << face_tolerance << '\n'
               << "spectral_transfer_closure_error=" << closure.spectral_max << '\n'
               << "qc_transfer_closure_error=" << closure.qc_max << '\n'
               << "qr_transfer_closure_error=" << closure.qr_max << '\n'
               << "spectral_transfer_closure_tolerance=" << closure.spectral_tolerance << '\n'
               << "qc_transfer_closure_tolerance=" << closure.qc_tolerance << '\n'
               << "qr_transfer_closure_tolerance=" << closure.qr_tolerance << '\n'
               << "cell_state_bytes=" << cell_state_bytes << '\n'
               << "face_transfer_bytes=" << face_transfer_bytes << '\n'
               << "total_auxiliary_bytes=" << total_auxiliary_bytes << '\n'
               << "passed=" << (passed ? 1 : 0) << '\n';
    }
    amrex::ParallelDescriptor::Barrier("SBM P1 diagnostic");
    if (!passed) {
        std::ostringstream message;
        message << "SBM P1 numerical qualification failed: step_count=" << step_count
                << ", mass_error=" << mass_error
                << ", compact_mass_error=" << compact_mass_error
                << ", initial_variation=" << initial_variation
                << ", transport_change=" << transport_change
                << ", projection_error=" << projection_error
                << ", face_projection_error=" << face_projection_error
                << ", spectral_transfer_closure_error=" << closure.spectral_max
                << ", qc_transfer_closure_error=" << closure.qc_max
                << ", qr_transfer_closure_error=" << closure.qr_max;
        amrex::Error(message.str());
    }
}

} // namespace

void ERF::initialize_sbm_auxiliary(const int lev)
{
    if (solverChoice.moisture_type != MoistureType::SBM) return;
    AMREX_ALWAYS_ASSERT(sbm_layout != nullptr);
    if (lev < 0) amrex::Error("SBM level index must be nonnegative");
    if (lev == 0) sbm_step_count = 0;
    if (!sbm_auxiliary) {
        sbm_auxiliary = std::make_unique<::erf_auxiliary::AuxiliaryStateManager>(sbm_layout->auxiliary_layout());
    }
    if (static_cast<int>(sbm_accepted_bulk_face_transfer.size()) <= lev) {
        sbm_accepted_bulk_face_transfer.resize(static_cast<std::size_t>(lev + 1));
    }
    if (!sbm_accepted_bulk_face_transfer[static_cast<std::size_t>(lev)]) {
        sbm_accepted_bulk_face_transfer[static_cast<std::size_t>(lev)] =
            std::make_unique<::erf_auxiliary::AuxiliaryFaceTransfer>();
        sbm_accepted_bulk_face_transfer[static_cast<std::size_t>(lev)]->define(grids[lev], dmap[lev], 2, 0);
    }
    if (static_cast<int>(sbm_initial_bulk_state.size()) <= lev) {
        sbm_initial_bulk_state.resize(static_cast<std::size_t>(lev + 1));
    }
    if (!sbm_initial_bulk_state[static_cast<std::size_t>(lev)]) {
        sbm_initial_bulk_state[static_cast<std::size_t>(lev)] =
            std::make_unique<amrex::MultiFab>(grids[lev], dmap[lev], 2, 0);
        sbm_initial_bulk_state[static_cast<std::size_t>(lev)]->setVal(Real(0.0));
    }
    if (!sbm_ownership) {
        sbm_ownership = std::make_unique<::erf_sbm::OwnershipRegistry>(true);
    }
    if (!sbm_auxiliary->has_level(lev)) sbm_auxiliary->define_level(lev, grids[lev], dmap[lev], 2);
    auto& aux = sbm_auxiliary->output(lev);
    auto& core = vars_new[lev][Vars::cons];
    // On restart MakeNewLevel* is the allocation phase.  ReadCheckpointFile
    // owns restoration and schema validation; never overwrite checkpointed
    // auxiliary data with a manufactured or empty state here.
    if (!restart_chkfile.empty()) return;
    const auto& projection = *sbm_layout;
    const auto& population = projection.populations().front();
    const int nbins = population.grid.nbins();
    const int offset = population.mass_offset;

    // The manufactured regression supplies a nonzero ERF carrier field while
    // production inputs retain the ordinary initialized velocity/momentum.
    // These are the same face-centered fields later handed to the transport
    // kernel, so the qualification cannot pass through an independent donor
    // velocity reconstruction.
    if (lev == 0 && solverChoice.sbm_manufactured_velocity != Real(0.0)) {
        vars_new[lev][Vars::xvel].setVal(solverChoice.sbm_manufactured_velocity);
        vars_old[lev][Vars::xvel].setVal(solverChoice.sbm_manufactured_velocity);
        vars_new[lev][Vars::yvel].setVal(Real(0.0));
        vars_old[lev][Vars::yvel].setVal(Real(0.0));
        vars_new[lev][Vars::zvel].setVal(Real(0.0));
        vars_old[lev][Vars::zvel].setVal(Real(0.0));
        avg_xmom[lev].setVal(solverChoice.sbm_manufactured_velocity);
        avg_ymom[lev].setVal(Real(0.0));
        avg_zmom[lev].setVal(Real(0.0));
    }

    // There is intentionally no bulk-to-spectrum guess.  An empty spectrum is
    // allowed only when both compact condensate fields are zero.  Tests and
    // explicit prototype runs may request the deterministic manufactured state.
    if (!solverChoice.sbm_manufactured_initialization) {
        const Real qc_max = core.max(RhoQ2_comp);
        const Real qc_min = core.min(RhoQ2_comp);
        const Real qr_max = core.max(RhoQ3_comp);
        const Real qr_min = core.min(RhoQ3_comp);
        if (qc_max != Real(0.0) || qc_min != Real(0.0) ||
            qr_max != Real(0.0) || qr_min != Real(0.0)) {
            amrex::Error("SBM requires explicit spectral initialization for nonzero qc/qr; no bulk-to-bin guess is permitted");
        }
    }

    for (MFIter mfi(aux); mfi.isValid(); ++mfi) {
        const Box box = mfi.validbox();
        const auto aux_arr = aux.array(mfi);
        const auto core_arr = core.array(mfi);
        const bool manufactured = solverChoice.sbm_manufactured_initialization;
        const Real xlo = geom[lev].ProbLo(0);
        const Real xlen = geom[lev].ProbHi(0) - xlo;
        const Real dx = geom[lev].CellSize(0);
        for (int b = 0; b < nbins; ++b) {
            const Real pivot = population.number_offset >= 0 ? population.grid.pivot(b) : Real(0.0);
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const Real rho = core_arr(i,j,k,Rho_comp);
                const Real x = xlo + (Real(i) + Real(0.5)) * dx;
                const Real variation = Real(1.0) + Real(0.25) *
                    std::sin(Real(6.2831853071795864769) * (x - xlo) / xlen);
                // Small deterministic positive values on both sides of the
                // projection split.  This is a transport manufactured field,
                // not a physical droplet or aerosol distribution.
                aux_arr(i,j,k,offset+b) = manufactured ?
                    rho * Real(1.0e-6) * Real(b+1) * variation : Real(0.0);
                if (population.number_offset >= 0) {
                    aux_arr(i,j,k,population.number_offset+b) =
                        pivot > Real(0.0) ? aux_arr(i,j,k,offset+b) / pivot : Real(0.0);
                }
            });
        }
    }
    aux.FillBoundary(geom[lev].periodicity());
    ::erf_sbm::SBMBulkProjection bulk_projection(*sbm_layout);
    for (MFIter mfi(aux); mfi.isValid(); ++mfi) {
        bulk_projection.apply_to_core(mfi.validbox(), aux.const_array(mfi), core.array(mfi));
    }
    core.FillBoundary(geom[lev].periodicity());
    amrex::MultiFab::Copy(*sbm_initial_bulk_state[static_cast<std::size_t>(lev)], core, RhoQ2_comp, 0, 1, 0);
    amrex::MultiFab::Copy(*sbm_initial_bulk_state[static_cast<std::size_t>(lev)], core, RhoQ3_comp, 1, 1, 0);

    Print() << "SBM P1 auxiliary state: components=" << sbm_layout->ncomp()
            << ", bins=" << nbins
            << ", cell-state bytes=" << sbm_auxiliary->state_resident_bytes()
            << ", face-transfer bytes=" << sbm_auxiliary->face_transfer_resident_bytes() +
               sbm_accepted_bulk_face_transfer[static_cast<std::size_t>(lev)]->resident_bytes()
            << ", total auxiliary bytes=" << sbm_auxiliary->resident_bytes() +
               sbm_accepted_bulk_face_transfer[static_cast<std::size_t>(lev)]->resident_bytes() + sbm_cell_bytes(*sbm_initial_bulk_state[static_cast<std::size_t>(lev)])
            << (solverChoice.sbm_manufactured_initialization ? " (manufactured initialization)" : " (empty initialization)")
            << std::endl;
}

void ERF::begin_sbm_step(const int lev, const amrex::MultiFab& core_old, const double old_time)
{
    if (solverChoice.moisture_type == MoistureType::SBM) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(sbm_auxiliary != nullptr && sbm_auxiliary->has_level(lev),
                                         "SBM auxiliary state must be initialized before stepping");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(static_cast<std::size_t>(lev) < sbm_initial_bulk_state.size() &&
                                         sbm_initial_bulk_state[static_cast<std::size_t>(lev)] != nullptr,
                                         "SBM compact baseline must be initialized before stepping");
        sbm_auxiliary->begin_step(lev, old_time);
        // ERF swaps vars_old/vars_new before entering advance_dycore.  The
        // explicit state_old argument is therefore the actual full-step old
        // compact state, even on the second and subsequent time steps.
        amrex::MultiFab::Copy(*sbm_initial_bulk_state[static_cast<std::size_t>(lev)], core_old,
                              RhoQ2_comp, 0, 1, 0);
        amrex::MultiFab::Copy(*sbm_initial_bulk_state[static_cast<std::size_t>(lev)], core_old,
                              RhoQ3_comp, 1, 1, 0);
        ++sbm_step_count;
    }
}

void ERF::advance_sbm_stage(const int lev,
                            Vector<MultiFab>& state_old,
                            Vector<MultiFab>& state_new,
                            Vector<MultiFab>& state_eval,
                            const double old_step_time,
                            const double old_stage_time,
                            const double new_stage_time,
                            const int stage,
                            const double full_step)
{
    if (solverChoice.moisture_type != MoistureType::SBM) return;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(sbm_auxiliary != nullptr && sbm_auxiliary->has_level(lev) && sbm_layout != nullptr,
                                     "SBM auxiliary state is not ready");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        sbm_ownership != nullptr &&
        sbm_ownership->owns(RhoQ2_comp, ::erf_sbm::NativeWritePath::Advection) &&
        sbm_ownership->owns(RhoQ3_comp, ::erf_sbm::NativeWritePath::Advection),
        "SBM compact cloud/rain ownership contract is not active");
    const bool anelastic = solverChoice.anelastic[lev] == 1;
    const auto context = anelastic ?
        ::erf_auxiliary::make_anelastic_stage(stage, old_step_time, old_stage_time,
                                             new_stage_time, full_step,
                                             &state_old[IntVars::cons], &state_eval[IntVars::cons]) :
        ::erf_auxiliary::make_compressible_stage(stage, old_step_time, old_stage_time,
                                                new_stage_time, full_step,
                                                &state_old[IntVars::cons], &state_eval[IntVars::cons]);

    // A fine-level WENO stencil must see the coarse spectrum at the actual
    // stage time.  The manager owns the coarse old/output bracket produced by
    // the completed coarse step and performs the temporal FillPatch before
    // this level constructs any face flux.
    if (lev > 0 && solverChoice.sbm_transport_method == "GroupedFCT_WENOZ3") {
        sbm_auxiliary->fill_stage_from_coarse(lev-1, lev,
                                              sbm_auxiliary->evaluation_time(lev),
                                              geom[lev-1], geom[lev], refRatio(lev-1));
    }

    // avg_*mom are ERF's actual dry-air carrier mass flux fields.  The
    // transport helper consumes them directly and never rebuilds rho*u from
    // the velocity MultiFabs.
    ::erf_sbm::advance_stage(*sbm_auxiliary, *sbm_layout, context,
                            state_eval[IntVars::cons], state_new[IntVars::cons],
                            avg_xmom[lev], avg_ymom[lev], avg_zmom[lev], geom[lev],
                            sbm_auxiliary->face_transfer_ledger(lev).stage(),
                            solverChoice.sbm_transport_method == "GroupedFCT_WENOZ3" ?
                                ::erf_sbm::TransportMethod::GroupedFCT_WENOZ3 :
                                ::erf_sbm::TransportMethod::DonorCell, lev,
                            solverChoice.sbm_diffusion_coeff, solverChoice.sbm_chunk_size);

    const auto& ledger = sbm_auxiliary->face_transfer_ledger(lev);
    // YAFluxRegister consumes instantaneous per-area fluxes and applies the
    // supplied dt/dx factor.  Register the accepted spectral face flux at
    // exactly the same stage weights as the auxiliary ledger; this is the
    // physical I=A*integral(F dt) contract without a second area or time
    // multiplication.
    if (solverChoice.coupling_type == CouplingType::TwoWay &&
        context.accepted_ledger_weight() != Real(0.0)) {
        auto& stage_flux = ledger.stage();
        const auto dx = geom[lev].CellSizeArray();
        const Real register_dt = static_cast<Real>(full_step * context.accepted_ledger_weight());
        for (MFIter mfi(state_new[IntVars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const std::array<FArrayBox const*, AMREX_SPACEDIM> fluxes{
                AMREX_D_DECL(&stage_flux.x()[mfi], &stage_flux.y()[mfi], &stage_flux.z()[mfi])};
            if (lev < finest_level && sbm_flux_reg[lev+1] != nullptr) {
                sbm_flux_reg[lev+1]->CrseAdd(mfi, fluxes, dx.data(), register_dt, RunOn::Device);
            }
            if (lev > 0 && sbm_flux_reg[lev] != nullptr) {
                sbm_flux_reg[lev]->FineAdd(mfi, fluxes, dx.data(), register_dt, RunOn::Device);
            }
        }
        Gpu::streamSynchronize();
    }
    const ::erf_sbm::SBMBulkProjection bulk_projection(*sbm_layout);
    bulk_projection.apply_to_face_transfer(ledger.accepted(), *sbm_accepted_bulk_face_transfer[static_cast<std::size_t>(lev)]);

    if (context.completes_level_step && solverChoice.sbm_manufactured_initialization &&
        !solverChoice.sbm_diagnostic_file.empty()) {
        auto& scratch = sbm_auxiliary->scratch(lev);
        const auto& old = sbm_auxiliary->old(lev);
        const auto& output = sbm_auxiliary->output(lev);
        const int ncomp = sbm_layout->ncomp();
        const Real cell_volume = geom[lev].CellSize(0) * geom[lev].CellSize(1) * geom[lev].CellSize(2);
        Real initial_mass = sbm_total_mass(old, ncomp, cell_volume);
        Real final_mass = sbm_total_mass(output, ncomp, cell_volume);
        const Real initial_variation = old.max(0) - old.min(0);
        const Real transport_change = sbm_max_change(scratch, output, old, ncomp);
        const Real projection_error = sbm_max_projection_error(*sbm_layout, output,
                                                               state_new[IntVars::cons]);
        const Real face_projection_error = bulk_projection.max_face_projection_error(
            ledger.accepted(), *sbm_accepted_bulk_face_transfer[static_cast<std::size_t>(lev)]);
        const auto closure = ::erf_sbm::evaluate_accepted_transfer_closure(
            *sbm_layout, old, output, ledger.accepted(),
            *sbm_initial_bulk_state[static_cast<std::size_t>(lev)], 0, 1,
            state_new[IntVars::cons], RhoQ2_comp, RhoQ3_comp,
            *sbm_accepted_bulk_face_transfer[static_cast<std::size_t>(lev)], geom[lev]);
        const Real compact_mass = cell_volume *
            (state_new[IntVars::cons].sum(RhoQ2_comp) +
             state_new[IntVars::cons].sum(RhoQ3_comp));
        write_sbm_diagnostic(solverChoice.sbm_diagnostic_file,
                             solverChoice.anelastic[lev] == 1, ncomp,
                             sbm_step_count,
                             initial_mass, final_mass, initial_variation,
                             transport_change, projection_error,
                             face_projection_error, closure, compact_mass, cell_volume,
                             sbm_auxiliary->state_resident_bytes() + sbm_cell_bytes(*sbm_initial_bulk_state[static_cast<std::size_t>(lev)]),
                             sbm_auxiliary->face_transfer_resident_bytes() +
                                 sbm_accepted_bulk_face_transfer[static_cast<std::size_t>(lev)]->resident_bytes(),
                             sbm_auxiliary->resident_bytes() +
                                 sbm_accepted_bulk_face_transfer[static_cast<std::size_t>(lev)]->resident_bytes() +
                                 sbm_cell_bytes(*sbm_initial_bulk_state[static_cast<std::size_t>(lev)]));
    }
}

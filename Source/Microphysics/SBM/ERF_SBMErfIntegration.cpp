#include "ERF.H"

#include "ERF_SBMBulkProjection.H"
#include "ERF_SBMContracts.H"
#include "ERF_SBMTransportPrototype.H"

#include <AMReX_MFParallelFor.H>

#include <algorithm>
#include <vector>

using namespace amrex;

void ERF::initialize_sbm_auxiliary(const int lev)
{
    if (solverChoice.moisture_type != MoistureType::SBM) return;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 0, "SBM P1 auxiliary state is level-0 only");
    AMREX_ALWAYS_ASSERT(sbm_layout != nullptr);
    if (!restart_chkfile.empty()) {
        amrex::Error("SBM P1 restart is unsupported: no auxiliary-state checkpoint/schema conversion is implemented");
    }
    if (!sbm_auxiliary) {
        sbm_auxiliary = std::make_unique<::erf_auxiliary::AuxiliaryStateManager>(sbm_layout->auxiliary_layout());
    }
    sbm_auxiliary->define_level(0, grids[0], dmap[0], 1);
    auto& aux = sbm_auxiliary->output(0);
    auto& core = vars_new[0][Vars::cons];
    const auto& projection = *sbm_layout;
    const int nbins = projection.populations().front().grid.nbins();
    const int offset = projection.populations().front().liquid_mass_offset;

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
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real rho = core_arr(i,j,k,Rho_comp);
            for (int b = 0; b < nbins; ++b) {
                // Small deterministic positive values on both sides of the
                // projection split.  This is a transport manufactured field,
                // not a physical droplet or aerosol distribution.
                aux_arr(i,j,k,offset+b) = manufactured ? rho * Real(1.0e-6) * Real(b+1) : Real(0.0);
            }
        });
    }
    aux.FillBoundary(geom[0].periodicity());
    ::erf_sbm::SBMBulkProjection bulk_projection(*sbm_layout);
    for (MFIter mfi(aux); mfi.isValid(); ++mfi) {
        bulk_projection.apply_to_core(mfi.validbox(), aux.const_array(mfi), core.array(mfi));
    }
    core.FillBoundary(geom[0].periodicity());

    Print() << "SBM P1 auxiliary state: components=" << sbm_layout->ncomp()
            << ", bins=" << nbins
            << ", resident bytes (old/eval/output/scratch)=" << sbm_auxiliary->resident_bytes()
            << (solverChoice.sbm_manufactured_initialization ? " (manufactured initialization)" : " (empty initialization)")
            << std::endl;
}

void ERF::begin_sbm_step(const int lev)
{
    if (solverChoice.moisture_type == MoistureType::SBM) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 0 && sbm_auxiliary != nullptr,
                                         "SBM auxiliary state must be initialized before stepping");
        sbm_auxiliary->begin_step(0);
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
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 0 && sbm_auxiliary != nullptr && sbm_layout != nullptr,
                                     "SBM auxiliary state is not ready");
    const bool anelastic = solverChoice.anelastic[lev] == 1;
    const auto context = anelastic ?
        ::erf_auxiliary::make_anelastic_stage(stage, old_step_time, old_stage_time,
                                             new_stage_time, full_step,
                                             &state_old[IntVars::cons], &state_eval[IntVars::cons]) :
        ::erf_auxiliary::make_compressible_stage(stage, old_step_time, old_stage_time,
                                                new_stage_time, full_step,
                                                &state_old[IntVars::cons], &state_eval[IntVars::cons]);

    // avg_*mom are ERF's actual dry-air carrier mass flux fields.  The
    // transport helper consumes them directly and never rebuilds rho*u from
    // the velocity MultiFabs.
    ::erf_sbm::advance_stage(*sbm_auxiliary, *sbm_layout, context,
                            state_eval[IntVars::cons], state_new[IntVars::cons],
                            avg_xmom[lev], avg_ymom[lev], avg_zmom[lev], geom[lev]);
}

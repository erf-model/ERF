#include "ERF_AuxiliaryFaceTransfer.H"

#include <AMReX_MultiFabUtil.H>

#include <stdexcept>

namespace erf_auxiliary {

void AuxiliaryFaceTransfer::define(const amrex::BoxArray& ba,
                                   const amrex::DistributionMapping& dm,
                                   const int ncomp, const int ngrow)
{
    if (ncomp <= 0 || ngrow < 0) {
        throw std::invalid_argument("auxiliary face transfer requires positive components and nonnegative ghosts");
    }
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        m_face[static_cast<std::size_t>(dir)] = std::make_unique<amrex::MultiFab>(
            amrex::convert(ba, amrex::IntVect::TheDimensionVector(dir)), dm, ncomp, ngrow);
        m_face[static_cast<std::size_t>(dir)]->setVal(0.0);
    }
    m_ncomp = ncomp;
    std::size_t points = 0;
    for (const auto& face : m_face) points += static_cast<std::size_t>(face->boxArray().numPts());
    m_resident_bytes = points * static_cast<std::size_t>(ncomp) * sizeof(amrex::Real);
}

void AuxiliaryFaceTransfer::setVal(const amrex::Real value)
{
    for (auto& face : m_face) face->setVal(value);
}

void AuxiliaryFaceTransferLedger::define(const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm,
                                         const int ncomp, const int nstages,
                                         const int ngrow)
{
    if (nstages <= 0 || nstages > 3) throw std::invalid_argument("auxiliary face ledger needs one to three stages");
    m_ncomp = ncomp;
    m_nstages = nstages;
    m_recorded = {{false, false, false}};
    m_stage = std::make_unique<AuxiliaryFaceTransfer>();
    m_stage->define(ba, dm, ncomp, ngrow);
    m_accepted = std::make_unique<AuxiliaryFaceTransfer>();
    m_accepted->define(ba, dm, ncomp, ngrow);
    m_resident_bytes = m_stage->resident_bytes() + m_accepted->resident_bytes();
}

void AuxiliaryFaceTransferLedger::begin_step()
{
    if (!m_accepted) throw std::logic_error("auxiliary face ledger is not defined");
    m_stage->setVal(0.0);
    m_accepted->setVal(0.0);
    m_recorded = {{false, false, false}};
}

void AuxiliaryFaceTransferLedger::record_stage(const StageContext& context,
                                               const AuxiliaryFaceTransfer& stage_flux)
{
    if (!m_accepted || stage_flux.ncomp() != m_ncomp ||
        context.stage_index < 0 || context.stage_index >= m_nstages ||
        m_recorded[static_cast<std::size_t>(context.stage_index)]) {
        throw std::invalid_argument("invalid or duplicate auxiliary face-transfer stage");
    }
    // Production passes the reusable stage scratch itself. Unit tests may pass
    // independent synthetic stage objects; copy those into the same scratch
    // without retaining their history.
    if (&stage_flux != m_stage.get()) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            amrex::MultiFab::Copy(m_stage->direction(dir), stage_flux.direction(dir),
                                  0, 0, m_ncomp, stage_flux.direction(dir).nGrowVect());
        }
    }
    const amrex::Real factor = static_cast<amrex::Real>(context.full_step *
                                                         context.accepted_ledger_weight());
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (factor != amrex::Real(0.0)) {
            amrex::MultiFab::Saxpy(m_accepted->direction(dir), factor,
                                   m_stage->direction(dir), 0, 0, m_ncomp, 0);
        }
    }
    m_recorded[static_cast<std::size_t>(context.stage_index)] = true;
}

} // namespace erf_auxiliary

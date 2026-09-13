#include "ERF_AuxiliaryFaceTransfer.H"

#include <AMReX_MultiFabUtil.H>

#include <algorithm>
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
    if (nstages <= 0) throw std::invalid_argument("auxiliary face ledger needs at least one stage");
    m_ncomp = ncomp;
    m_stages.clear();
    m_stages.reserve(static_cast<std::size_t>(nstages));
    m_recorded.assign(static_cast<std::size_t>(nstages), false);
    for (int stage = 0; stage < nstages; ++stage) {
        auto flux = std::make_unique<AuxiliaryFaceTransfer>();
        flux->define(ba, dm, ncomp, ngrow);
        m_stages.push_back(std::move(flux));
    }
    m_accepted = std::make_unique<AuxiliaryFaceTransfer>();
    m_accepted->define(ba, dm, ncomp, ngrow);
    m_resident_bytes = m_accepted->resident_bytes();
    for (const auto& stage : m_stages) m_resident_bytes += stage->resident_bytes();
}

void AuxiliaryFaceTransferLedger::begin_step()
{
    if (!m_accepted) throw std::logic_error("auxiliary face ledger is not defined");
    m_accepted->setVal(0.0);
    std::fill(m_recorded.begin(), m_recorded.end(), false);
}

void AuxiliaryFaceTransferLedger::record_stage(const StageContext& context,
                                               const AuxiliaryFaceTransfer& stage_flux)
{
    if (!m_accepted || stage_flux.ncomp() != m_ncomp ||
        context.stage_index < 0 || context.stage_index >= static_cast<int>(m_stages.size()) ||
        m_recorded[static_cast<std::size_t>(context.stage_index)]) {
        throw std::invalid_argument("invalid or duplicate auxiliary face-transfer stage");
    }
    auto& retained = *m_stages[static_cast<std::size_t>(context.stage_index)];
    const amrex::Real factor = static_cast<amrex::Real>(context.full_step *
                                                         context.accepted_ledger_weight());
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        amrex::MultiFab::Copy(retained.direction(dir), stage_flux.direction(dir),
                              0, 0, m_ncomp, stage_flux.direction(dir).nGrowVect());
        if (factor != amrex::Real(0.0)) {
            amrex::MultiFab::Saxpy(m_accepted->direction(dir), factor,
                                   stage_flux.direction(dir), 0, 0, m_ncomp, 0);
        }
    }
    m_recorded[static_cast<std::size_t>(context.stage_index)] = true;
}

} // namespace erf_auxiliary

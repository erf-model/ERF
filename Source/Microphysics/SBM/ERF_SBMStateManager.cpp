#include "ERF_SBMStateManager.H"

#include <AMReX_MFIter.H>

#include <utility>

namespace erf_sbm {

SBMStateManager::SBMStateManager(SBMLayout layout, const int nlevels)
    : m_layout(std::move(layout)), m_projection(m_layout)
{
    if (nlevels <= 0) throw std::invalid_argument("SBM state manager needs at least one level");
    m_state.resize(static_cast<std::size_t>(nlevels));
}

bool SBMStateManager::is_defined(const int lev) const
{
    if (lev < 0 || lev >= nlevels()) throw std::out_of_range("SBM level is outside the manager");
    return static_cast<bool>(m_state[static_cast<std::size_t>(lev)]);
}

void SBMStateManager::define(const int lev, const amrex::BoxArray& grids,
                             const amrex::DistributionMapping& mapping,
                             const int ngrow)
{
    if (lev < 0 || lev >= nlevels()) throw std::out_of_range("SBM level is outside the manager");
    if (ngrow < 0) throw std::invalid_argument("SBM ghost width cannot be negative");
    auto& level_state = m_state[static_cast<std::size_t>(lev)];
    if (level_state) throw std::logic_error("SBM level state is already defined");
    level_state = std::make_unique<amrex::MultiFab>(grids, mapping, m_layout.ncomp(), ngrow);
    // All-zero is a valid physical state, so allocation never marks it as uninitialized.
    level_state->setVal(amrex::Real(0.0));
}

void SBMStateManager::destroy(const int lev)
{
    if (lev < 0 || lev >= nlevels()) throw std::out_of_range("SBM level is outside the manager");
    m_state[static_cast<std::size_t>(lev)].reset();
}

amrex::MultiFab& SBMStateManager::state(const int lev)
{
    if (!is_defined(lev)) throw std::logic_error("SBM level state is not defined");
    return *m_state[static_cast<std::size_t>(lev)];
}

const amrex::MultiFab& SBMStateManager::state(const int lev) const
{
    if (!is_defined(lev)) throw std::logic_error("SBM level state is not defined");
    return *m_state[static_cast<std::size_t>(lev)];
}

void SBMStateManager::project_to_core(const int lev, amrex::MultiFab& core,
                                     const int qc_component,
                                     const int qr_component) const
{
    const auto& spectrum = state(lev);
    if (qc_component < 0 || qr_component < 0 ||
        qc_component >= core.nComp() || qr_component >= core.nComp()) {
        throw std::invalid_argument("SBM compact projection target is outside core state");
    }
    if (!(spectrum.boxArray() == core.boxArray()) ||
        !(spectrum.DistributionMap() == core.DistributionMap())) {
        throw std::invalid_argument("SBM spectrum and core state must share a level layout");
    }
    for (amrex::MFIter mfi(core); mfi.isValid(); ++mfi) {
        m_projection.apply_to_core(mfi.validbox(), spectrum.const_array(mfi),
                                   core.array(mfi), qc_component, qr_component);
    }
}

} // namespace erf_sbm

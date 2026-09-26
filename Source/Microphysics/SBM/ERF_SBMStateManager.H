#ifndef ERF_SBM_STATE_MANAGER_H_
#define ERF_SBM_STATE_MANAGER_H_

#include "ERF_SBMBulkProjection.H"

#include <AMReX_MultiFab.H>

#include <memory>
#include <stdexcept>
#include <vector>

namespace erf_sbm {

/** Owns the authoritative, runtime-sized spectral state, one MultiFab per level. */
class SBMStateManager {
public:
    SBMStateManager(SBMLayout layout, int nlevels);

    [[nodiscard]] const SBMLayout& layout() const noexcept { return m_layout; }
    [[nodiscard]] int nlevels() const noexcept { return static_cast<int>(m_state.size()); }
    [[nodiscard]] bool is_defined(int lev) const;

    void define(int lev, const amrex::BoxArray& grids,
                const amrex::DistributionMapping& mapping, int ngrow = 0);
    void destroy(int lev);

    [[nodiscard]] amrex::MultiFab& state(int lev);
    [[nodiscard]] const amrex::MultiFab& state(int lev) const;

    /** One-way projection from the spectrum into ERF's compact qc and qr lanes. */
    void project_to_core(int lev, amrex::MultiFab& core,
                         int qc_component, int qr_component) const;

private:
    SBMLayout m_layout;
    SBMBulkProjection m_projection;
    std::vector<std::unique_ptr<amrex::MultiFab>> m_state;
};

} // namespace erf_sbm

#endif

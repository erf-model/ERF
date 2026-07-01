#ifndef ERF_SHOC_DRIVER_H_
#define ERF_SHOC_DRIVER_H_

#include "ERF_ShocDiagnostics.H"
#include "ERF_ShocTypes.H"
#include "ERF_ShocPreprocess.H"
#include "ERF_ShocColumnData.H"

#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>

#include <memory>

bool shoc_boxarray_spans_full_height (const amrex::BoxArray& ba,
                                      const amrex::Box& domain);

inline bool
shoc_layout_requires_number_closure (const MoistureComponentIndices& indices,
                                     int ncomp) noexcept
{
    return shoc_valid_comp(indices.nc, ncomp) || shoc_valid_comp(indices.ni, ncomp);
}

inline bool
shoc_driver_host_diffusion_moisture_supported (const ShocRuntimeOptions& opts,
                                               MoistureType moisture_type,
                                               std::string& error_message)
{
    if (!shoc_uses_host_diffusion(opts.transport_mode) || moisture_type == MoistureType::None) {
        return true;
    }

    error_message =
        "Native SHOC host_diffusion with moisture is not yet supported because SHOC does not own cloud macrophysics in this mode while SHOC-family microphysics condensation is suppressed. Use state_update for moist SHOC runs, or run host_diffusion only for dry cases until a transport-aware microphysics ownership predicate is implemented.";
    return false;
}

inline bool
shoc_driver_state_update_layout_supported (const ShocRuntimeOptions& opts,
                                           const MoistureComponentIndices& indices,
                                           int ncomp,
                                           std::string& error_message)
{
    if (!shoc_uses_state_update(opts.transport_mode) ||
        !shoc_layout_requires_number_closure(indices, ncomp)) {
        return true;
    }

    error_message =
        "Native SHOC state_update does not yet support number-aware microphysics layouts with cloud/ice number concentrations. A number closure is required before SHOC can update cloud mass in these layouts.";
    return false;
}

class ShocDriver
{
public:
    explicit ShocDriver (int lev, const SolverChoice& solver_choice);

    void advance (amrex::MultiFab& cons,
                  amrex::MultiFab& xvel,
                  amrex::MultiFab& yvel,
                  amrex::MultiFab& zvel,
                  amrex::MultiFab* tau13,
                  amrex::MultiFab* tau23,
                  amrex::MultiFab* hfx3,
                  amrex::MultiFab* qfx3,
                  amrex::MultiFab* eddy_diffs,
                  amrex::MultiFab& z_phys_nd,
                  const amrex::Geometry& geom,
                  amrex::Real dt);

    void set_eddy_diffs () const;
    void set_diff_stresses () const;
    void add_fast_tend (amrex::Vector<amrex::MultiFab>& S_rhs) const;
    void add_slow_tend (const amrex::MFIter& mfi,
                        const amrex::Box& tbx,
                        const amrex::Array4<amrex::Real>& cell_rhs) const;
    bool uses_state_update () const;
    bool uses_host_diffusion () const;
    bool uses_momentum_state_update () const;
    bool uses_momentum_host_diffusion () const;
    bool disables_momentum_transport () const;
    bool owns_scalar_surface_fluxes () const;
    bool owns_momentum_surface_stresses () const;
    bool needs_host_surface_momentum_stresses () const;
    bool owns_surface_fluxes () const;
    bool debug_summary_enabled () const { return m_opts.debug_summary; }
    bool has_native_diagnostics () const { return m_eddy_coeffs_cc.isDefined(); }
    const amrex::MultiFab& native_diagnostics () const { return m_eddy_coeffs_cc; }
    // Native SHOC pblh is diagnosed in meters above local ground (AGL).
    const amrex::MultiFab& pblh_diagnostics () const { return m_pblh_cc; }
    const amrex::MultiFab& shoc_cldfrac_diagnostics () const { return m_shoc_cldfrac_cc; }
    const amrex::MultiFab& shoc_ql_diagnostics () const { return m_shoc_ql_cc; }
    const amrex::MultiFab& shoc_ql2_diagnostics () const { return m_shoc_ql2_cc; }
    const amrex::MultiFab& shoc_cond_diagnostics () const { return m_shoc_cond_cc; }
    const amrex::MultiFab& w_sec_diagnostics () const { return m_w_sec_cc; }
    const amrex::MultiFab& wqls_sec_diagnostics () const { return m_wqls_sec_cc; }
    const amrex::MultiFab& wthv_sec_diagnostics () const { return m_wthv_sec_cc; }
    const amrex::MultiFab& thl_sec_diagnostics () const { return m_thl_sec_cc; }
    const amrex::MultiFab& qw_sec_diagnostics () const { return m_qw_sec_cc; }
    const amrex::MultiFab& qwthl_sec_diagnostics () const { return m_qwthl_sec_cc; }
    const amrex::MultiFab& wthl_sec_diagnostics () const { return m_wthl_sec_cc; }
    const amrex::MultiFab& wqw_sec_diagnostics () const { return m_wqw_sec_cc; }
    const amrex::MultiFab& w3_diagnostics () const { return m_w3_cc; }
    const amrex::MultiFab& brunt_diagnostics () const { return m_brunt_cc; }
    const amrex::MultiFab& isotropy_diagnostics () const { return m_isotropy_cc; }
    const amrex::MultiFab& shear_prod_diagnostics () const { return m_shear_prod_cc; }
    const amrex::MultiFab& buoy_prod_diagnostics () const { return m_buoy_prod_cc; }
    const amrex::MultiFab& diss_tke_diagnostics () const { return m_diss_tke_cc; }

private:
    int m_lev;
    ShocRuntimeOptions m_opts;
    MoistureType m_moisture_type;
    MoistureComponentIndices m_moisture_indices;

    amrex::MultiFab m_theta_tend_cc;
    amrex::MultiFab m_qv_tend_cc;
    amrex::MultiFab m_qc_tend_cc;
    amrex::MultiFab m_qi_tend_cc;
    amrex::MultiFab m_tke_tend_cc;
    amrex::MultiFab m_u_tend_cc;
    amrex::MultiFab m_v_tend_cc;
    amrex::MultiFab m_u_tend_fc;
    amrex::MultiFab m_v_tend_fc;
    amrex::MultiFab m_eddy_coeffs_cc;
    amrex::MultiFab m_prev_turb_cc;
    amrex::MultiFab m_prev_wthv_sec_cc;
    amrex::MultiFab m_pblh_cc;
    amrex::MultiFab m_shoc_cldfrac_cc;
    amrex::MultiFab m_shoc_ql_cc;
    amrex::MultiFab m_shoc_ql2_cc;
    amrex::MultiFab m_shoc_cond_cc;
    amrex::MultiFab m_w_sec_cc;
    amrex::MultiFab m_wqls_sec_cc;
    amrex::MultiFab m_wthv_sec_cc;
    amrex::MultiFab m_thl_sec_cc;
    amrex::MultiFab m_qw_sec_cc;
    amrex::MultiFab m_qwthl_sec_cc;
    amrex::MultiFab m_wthl_sec_cc;
    amrex::MultiFab m_wqw_sec_cc;
    amrex::MultiFab m_w3_cc;
    amrex::MultiFab m_brunt_cc;
    amrex::MultiFab m_isotropy_cc;
    amrex::MultiFab m_shear_prod_cc;
    amrex::MultiFab m_buoy_prod_cc;
    amrex::MultiFab m_diss_tke_cc;
    amrex::MultiFab* m_cons_ptr = nullptr;
    amrex::MultiFab* m_hfx3_ptr = nullptr;
    amrex::MultiFab* m_qfx3_ptr = nullptr;
    amrex::MultiFab* m_tau13_ptr = nullptr;
    amrex::MultiFab* m_tau23_ptr = nullptr;
    amrex::MultiFab* m_eddy_diffs_ptr = nullptr;
    amrex::Vector<std::unique_ptr<ShocColumnWorkspace>> m_column_workspaces;
    int m_advance_calls = 0;
    bool m_prev_turb_valid = false;

    void ensure_storage (const amrex::MultiFab& cons,
                         const amrex::MultiFab& xvel,
                         const amrex::MultiFab& yvel,
                         const amrex::MultiFab& eddy_diffs);
    void seed_carried_turbulence (ShocColumnData& col,
                                  const amrex::MFIter& mfi,
                                  const amrex::MultiFab& cons,
                                  const amrex::MultiFab& eddy_diffs) const;
    void store_carried_turbulence (const ShocColumnData& col,
                                   const amrex::MFIter& mfi);
    void seed_carried_buoyancy_flux (ShocColumnData& col,
                                     const amrex::MFIter& mfi) const;
    void store_carried_buoyancy_flux (const ShocColumnData& col,
                                      const amrex::MFIter& mfi);
    void apply_state_update (amrex::MultiFab& cons,
                             amrex::MultiFab& xvel,
                             amrex::MultiFab& yvel,
                             amrex::Real dt) const;
    void debug_check_bad_column (const ShocColumnData& col,
                                 const amrex::MFIter& mfi,
                                 const amrex::MultiFab& z_phys_nd,
                                 const amrex::MultiFab* hfx3,
                                 const amrex::MultiFab* qfx3,
                                 const amrex::MultiFab* tau13,
                                 const amrex::MultiFab* tau23,
                                 const amrex::Geometry& geom,
                                 amrex::Real dt) const;
    void print_debug_summary (amrex::Real dt) const;
};
#endif

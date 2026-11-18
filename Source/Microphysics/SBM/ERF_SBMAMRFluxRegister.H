#ifndef ERF_SBM_AMR_FLUX_REGISTER_H_
#define ERF_SBM_AMR_FLUX_REGISTER_H_

#include "AuxiliaryState/ERF_AuxiliaryMappedTransport.H"
#include "AuxiliaryState/ERF_AuxiliaryStage.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_YAFluxRegister.H>

#include <memory>
#include <cstdint>
#include <string>
#include <vector>

namespace erf_sbm {

/**
 * Transactional owner for one SBM coarse/fine interface.
 *
 * Accepted chunk-local mapped stage rates are added directly to YA using the
 * stage's completed-step time coefficient.  The register is only refluxable
 * after the coarse step and every expected fine substep have been accepted.
 * The authoritative spectrum and the mapped candidate remain externally
 * owned.
 */
class SBMAMRFluxRegister {
public:
    SBMAMRFluxRegister () = default;
    SBMAMRFluxRegister (const SBMAMRFluxRegister&) = delete;
    SBMAMRFluxRegister& operator= (const SBMAMRFluxRegister&) = delete;

    void define (const amrex::BoxArray& fine_ba,
                 const amrex::BoxArray& coarse_ba,
                 const amrex::DistributionMapping& fine_dm,
                 const amrex::DistributionMapping& coarse_dm,
                 const amrex::Geometry& fine_geom,
                 const amrex::Geometry& coarse_geom,
                 const amrex::IntVect& ratio,
                 int fine_level,
                 int ncomp);
    void clear ();

    bool begin_interval (double coarse_begin, double coarse_end,
                         int expected_fine_steps, std::string& diagnostic);
    bool begin_level_step (bool coarse_side, double step_begin,
                           double step_end, std::string& diagnostic);
    bool begin_stage (bool coarse_side,
                      erf_auxiliary::HostIntegrator method, int stage,
                      double step_begin, std::string& diagnostic);
    bool add_stage_rates (bool coarse_side,
                          const erf_auxiliary::MappedFaceFluxRate& accepted_rate,
                          const std::vector<int>& register_components,
                          double completed_ledger_time,
                          std::string& diagnostic);
    bool finish_stage (bool coarse_side, std::string& diagnostic);
    bool accept_level_step (bool coarse_side, double step_begin,
                            double step_end, std::string& diagnostic);

    [[nodiscard]] bool ready_to_reflux () const noexcept;
    [[nodiscard]] bool has_outstanding_amount () const noexcept;
    /** Actual smallest physical timestep represented by the active interval. */
    [[nodiscard]] double minimum_fine_dt () const noexcept
    { return m_min_fine_dt; }
    [[nodiscard]] std::uint64_t reflux_call_count () const noexcept
    { return m_reflux_call_count; }
    [[nodiscard]] std::uint64_t reflux_sync_count () const noexcept
    { return m_reflux_sync_count; }
    /** Coarse cells covered by the fine level (1) or not (0), one periodic ghost. */
    [[nodiscard]] const amrex::MultiFab& covered_parent_mask () const noexcept
    { return m_covered_parent_mask; }
    bool reflux (amrex::MultiFab& mapped_candidate,
                 std::string& diagnostic);
    /** Apply selected global components into matching local slots of a
     * bounded correction MultiFab; all components must be covered before
     * the pending amount can be consumed.
     */
    bool reflux_components (amrex::MultiFab& mapped_correction,
                            const std::vector<int>& global_components,
                            std::string& diagnostic);
    bool can_accept_reflux (std::string& diagnostic) const;
    /** Consume the pending amount only after the refluxed full candidate admits. */
    bool accept_reflux (std::string& diagnostic);
    /** Drop a rejected private candidate and its unpublishable interval. */
    bool discard_interval (std::string& diagnostic);

private:
    struct SideState {
        bool step_active{false};
        bool stage_active{false};
        bool step_accepted{false};
        double step_begin{0.0};
        double step_end{0.0};
        double current_stage_time{0.0};
        erf_auxiliary::HostIntegrator method{
            erf_auxiliary::HostIntegrator::CompressibleRK3};
        int next_stage{0};
        std::vector<unsigned char> components_seen;
    };

    [[nodiscard]] int stage_count (erf_auxiliary::HostIntegrator method) const noexcept;
    void reset_side (SideState& side) const;
    bool collective_check (bool local_ok, const std::string& local_diagnostic,
                           const char* remote_diagnostic,
                           std::string& diagnostic) const;

    std::unique_ptr<amrex::YAFluxRegister> m_register;
    amrex::MultiFab m_covered_parent_mask;
    amrex::BoxArray m_fine_ba;
    amrex::BoxArray m_coarse_ba;
    amrex::DistributionMapping m_fine_dm;
    amrex::DistributionMapping m_coarse_dm;
    amrex::Geometry m_fine_geom;
    amrex::Geometry m_coarse_geom;
    amrex::IntVect m_ratio{1};
    int m_fine_level{0};
    int m_ncomp{0};

    bool m_interval_active{false};
    bool m_coarse_step_accepted{false};
    bool m_reflux_applied{false};
    std::vector<unsigned char> m_reflux_components_seen;
    bool m_consumed{true};
    double m_interval_begin{0.0};
    double m_interval_end{0.0};
    double m_min_fine_dt{0.0};
    int m_expected_fine_steps{0};
    int m_accepted_fine_steps{0};
    std::uint64_t m_reflux_call_count{0};
    std::uint64_t m_reflux_sync_count{0};
    SideState m_coarse;
    SideState m_fine;
};

} // namespace erf_sbm

#endif

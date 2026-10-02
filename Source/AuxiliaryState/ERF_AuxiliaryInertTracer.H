#ifndef ERF_AUXILIARY_INERT_TRACER_H_
#define ERF_AUXILIARY_INERT_TRACER_H_

#include "ERF_AuxiliaryStage.H"
#include "ERF_AuxiliaryStateLayout.H"

#include <ERF_AdvStruct.H>
#include <AMReX_Geometry.H>

#include <memory>
#include <string>
#include <vector>

namespace erf_auxiliary {

/** Test-only, non-SBM proof consumer for the generic mapped transport substrate. */
class AuxiliaryInertTracer {
public:
    explicit AuxiliaryInertTracer (int number_of_levels);
    ~AuxiliaryInertTracer ();
    AuxiliaryInertTracer (const AuxiliaryInertTracer&) = delete;
    AuxiliaryInertTracer& operator= (const AuxiliaryInertTracer&) = delete;

    [[nodiscard]] static const AuxiliaryStateLayout& layout ();

    void define (int level,
                 const amrex::BoxArray& cell_ba,
                 const amrex::DistributionMapping& dm);

    /** Build or refresh this proof consumer's time-independent mapped measure. */
    bool rebuild_static_measure (int level,
                                 const amrex::MultiFab& detJ,
                                 const amrex::MultiFab& mx,
                                 const amrex::MultiFab& my,
                                 std::string& diagnostic);

    void initialize (int level,
                     const amrex::MultiFab& conserved,
                     const amrex::Geometry& geometry);

    void advance_stage (int level,
                        HostIntegrator method,
                        int stage,
                        double step_old_time,
                        double input_time,
                        double target_time,
                        double host_stage_interval,
                        const amrex::MultiFab& conserved_anchor,
                        const amrex::MultiFab& conserved_input,
                        const amrex::MultiFab& conserved_target,
                        const amrex::MultiFab& avg_xmom,
                        const amrex::MultiFab& avg_ymom,
                        const amrex::MultiFab& avg_zmom,
                        const AdvChoice& advection,
                        const amrex::Geometry& geometry);

    void destroy (int level);
    [[nodiscard]] bool is_defined (int level) const;
    [[nodiscard]] bool measure_is_ready (int level) const;
    /** Test-only inspection of the proof consumer's stored static measure. */
    [[nodiscard]] const amrex::MultiFab& static_measure (int level) const;
    [[nodiscard]] const amrex::MultiFab& state (int level) const;

private:
    struct LevelStorage;
    std::vector<std::unique_ptr<LevelStorage>> m_levels;
};

} // namespace erf_auxiliary

#endif

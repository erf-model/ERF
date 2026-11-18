#ifndef ERF_FORCE_RESTORE_PRECIP_RATE_H_
#define ERF_FORCE_RESTORE_PRECIP_RATE_H_

#include <memory>

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <ERF_SurfacePrecipitation.H>

struct RadChoice;

/**
 * @brief The precipitation one force-restore step takes from a rate field: the rate
 *        [kg/m^2/s] in component @p comp of @p rate, times @p scale. A null @p rate, or a
 *        zero @p scale, is no precipitation.
 */
struct ForceRestorePrecipInput
{
    const amrex::MultiFab* rate = nullptr;
    int comp = 0;
    amrex::Real scale = amrex::Real(0.0);
};

/**
 * @brief The surface precipitation rate of each level's step, for the force-restore
 *        hydrology (ERF_ForceRestoreHydrologyState.H) and snow layer
 *        (ERF_ForceRestoreSnowState.H), with erf.radiation.seb_precipitation_source.
 *
 * Owned by TwoStreamRadiation. Per level, and per step:
 *  - begin_step(...) at the start of ERF::Advance, before the dycore and microphysics,
 *    stores the surface precipitation accumulation of every column;
 *  - end_step(...) after microphysics turns its change since that start into the rate;
 *  - the hydrology and snow steps take the rate through rain(...) and frozen(...);
 *  - mark_used(...) after them closes the step.
 * A rate is used once. When no step followed end_step (the balance skips a new level's
 * first step), begin_step keeps the earlier start and the next end_step gives the mean rate
 * over both steps; the scale of rain(...) and frozen(...) is then the time since that start
 * over the step, so the carried precipitation is counted once, in full. When a step runs
 * without a new end_step (a level whose radiation is interpolated has neither begin_step
 * nor end_step), the scale is zero: the rate of an earlier step is not applied again.
 *
 * Without the snow layer the rate has one component, the total precipitation. With
 * erf.radiation.seb_snow_enable it has two, the rain (Rain) and the frozen part (Frozen),
 * filled in one pass over the accumulators. The rate needs no checkpoint or regrid
 * handling: both captures happen inside one ERF::Advance of the level.
 */
class ForceRestorePrecipRate
{
public:
    /** @brief Components of rate(lev). Frozen exists only with the snow layer. */
    enum Comp : int { Rain = 0, Frozen = 1 };

    void resize (int nlevs_max);

    /**
     * @brief Allocate the level's fields; a no-op without
     * erf.radiation.seb_precipitation_source.
     */
    void define_level (int lev, const RadChoice& rc, const amrex::BoxArray& ba2d,
                       const amrex::DistributionMapping& dm);

    /** @brief True once define_level allocated @p lev. */
    [[nodiscard]] bool defined (int lev) const noexcept {
        return lev >= 0 && lev < static_cast<int>(m_rate.size()) && m_rate[lev] != nullptr;
    }

    /** @brief True when the rate is split into rain and frozen precipitation. */
    [[nodiscard]] bool split () const noexcept { return m_split; }

    /**
     * @brief Store the accumulation [kg/m^2] at the start of the level's step. When no step
     * used the rate of the previous end_step, the earlier start is kept.
     */
    void begin_step (int lev, const SurfacePrecipAccumulationSources& sources, int klo);

    /**
     * @brief The rate [kg/m^2/s]: the accumulation now minus the one begin_step stored,
     * over the time since that start; never negative. Zero when begin_step has not run for
     * this step.
     */
    void end_step (int lev, const SurfacePrecipAccumulationSources& sources, int klo,
                   amrex::Real dt);

    /** @brief The rain a step of length @p dt takes (the total without the snow layer). */
    [[nodiscard]] ForceRestorePrecipInput rain (int lev, amrex::Real dt) const noexcept {
        return input(lev, Rain, dt);
    }

    /** @brief The frozen precipitation a step of length @p dt takes; none unless split(). */
    [[nodiscard]] ForceRestorePrecipInput frozen (int lev, amrex::Real dt) const noexcept {
        return m_split ? input(lev, Frozen, dt) : ForceRestorePrecipInput{};
    }

    /** @brief The steps of this level took the rate: it is not applied again. */
    void mark_used (int lev) noexcept;

    /** @brief The rate of the last end_step [kg/m^2/s] (Comp components), or nullptr. For output. */
    [[nodiscard]] const amrex::MultiFab* rate (int lev) const noexcept {
        return defined(lev) ? m_rate[lev].get() : nullptr;
    }

private:
    [[nodiscard]] ForceRestorePrecipInput input (int lev, int comp, amrex::Real dt) const noexcept;
    void fill_accumulation (amrex::MultiFab& accum, const SurfacePrecipAccumulationSources& sources,
                            int klo) const;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_accum_start;  ///< accumulation at the step start
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_accum_now;    ///< work field of end_step
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_rate;         ///< rate of the last end_step
    amrex::Vector<int> m_started;                                   ///< begin_step ran this step
    amrex::Vector<int> m_unused;                                    ///< end_step's rate not yet taken by a step
    amrex::Vector<amrex::Real> m_elapsed;                           ///< time since m_accum_start [s]
    bool m_split = false;                                           ///< rain and frozen apart (seb_snow_enable)
    bool m_noted_no_frozen = false;                                 ///< the no-frozen-accumulator note was printed
};

#endif // ERF_FORCE_RESTORE_PRECIP_RATE_H_

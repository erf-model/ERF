#ifndef ERF_FORCE_RESTORE_SNOW_STATE_H_
#define ERF_FORCE_RESTORE_SNOW_STATE_H_

#include <memory>
#include <string>

#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <ERF_ForceRestoreSnow.H>
#include <ERF_SurfacePrecipitation.H>

struct RadChoice;

/**
 * @brief The per-level state of the force-restore snow layer (ERF_ForceRestoreSnow.H,
 *        erf.radiation.seb_snow_enable): the snow water equivalent, the frozen precipitation
 *        of the current step, and the running totals that close each column's snow budget.
 *
 * Owned by TwoStreamRadiation, beside the force-restore hydrology
 * (ERF_ForceRestoreHydrologyState.H), and called at the same three points of a step:
 *  - begin_step(...) stores the frozen surface precipitation accumulation (snow, graupel,
 *    hail) of every column;
 *  - end_step(...) turns its change over the step into the frozen precipitation rate;
 *  - step(...) after the balance, and before the hydrology step, applies one snow step to
 *    every open column: it changes T_s where snow melts and leaves the water for the soil
 *    in soil_water(lev), which the hydrology step takes.
 * The snow water equivalent and the totals are checkpointed, and TwoStreamRadiation lists
 * them with its prognostic fields for regrid, coarse-fine filling and averaging down.
 */
class ForceRestoreSnow
{
public:
    /** @brief Components of totals(lev), each a running total [kg/m^2] since the start. */
    enum Total : int { Snowfall = 0, Sublimation, Melt, NumTotals };

    /** @brief Domain sums of one step on one level (rank-reduced), for the log. */
    struct StepSums {
        amrex::Real snowfall = 0.0;      ///< [kg/m^2]
        amrex::Real sublimation = 0.0;   ///< [kg/m^2]
        amrex::Real melt = 0.0;          ///< [kg/m^2]
        amrex::Real swe = 0.0;           ///< sum of W_n after the step [kg/m^2]
        amrex::Real max_residual = 0.0;  ///< max |column budget residual| [kg/m^2]
        amrex::Long columns = 0;         ///< open columns updated
    };

    /** @brief The parameters of one step, from the erf.radiation.seb_* inputs. */
    [[nodiscard]] static ForceRestoreSnowParams params (const RadChoice& rc);

    void resize (int nlevs_max);

    /**
     * @brief Allocate the level's fields, with the snow water equivalent at
     * erf.radiation.seb_snow_initial_swe_kg_m2; a no-op unless seb_snow_enable.
     */
    void define_level (int lev, const RadChoice& rc, const amrex::BoxArray& ba2d,
                       const amrex::DistributionMapping& dm, const amrex::IntVect& ng);

    /** @brief True once define_level allocated @p lev. */
    [[nodiscard]] bool defined (int lev) const noexcept {
        return lev >= 0 && lev < static_cast<int>(m_swe.size()) && m_swe[lev] != nullptr;
    }

    /**
     * @brief Store the frozen accumulation [kg/m^2] at the start of the level's step. When
     * step() did not use the snowfall of the previous step, the earlier start is kept, so
     * that snowfall is carried into the next step (as the hydrology carries its rain).
     */
    void begin_step (int lev, const SurfacePrecipAccumulationSources& sources, int klo);

    /**
     * @brief The frozen precipitation rate [kg/m^2/s]: the accumulation now minus the one
     * begin_step stored, over the time since that start. Zero when begin_step has not run
     * for this step.
     */
    void end_step (int lev, const SurfacePrecipAccumulationSources& sources, int klo,
                   amrex::Real dt);

    /**
     * @brief One snow step on every open column of @p lev, after the balance has updated
     * T_s. A footprint cell (footprint > 0.5) keeps its state and hands the soil no water.
     *
     * @param t_s       surface temperature [K]; lowered where snow melts
     * @param lh        latent heat flux the balance used [W/m^2]; E = LH / L_v
     * @param footprint building footprints, or nullptr
     * @param dt        the step [s]
     * @param prm       params(rc)
     * @param[out] t_s_sum_local sum of T_s over this rank's open columns after the step
     * @param[out] t_s_max_local max |T_s| over this rank's open columns after the step
     * @return the step's domain sums, reduced over ranks
     */
    StepSums step (int lev, amrex::MultiFab& t_s, const amrex::MultiFab& lh,
                   const amrex::MultiFab* footprint, amrex::Real dt,
                   const ForceRestoreSnowParams& prm,
                   amrex::Real& t_s_sum_local, amrex::Real& t_s_max_local);

    /**
     * @brief Write the surface albedo with snow into @p alb_sw, which holds the snow-free
     * albedo: DRM95 Eq. (8) with veg = 0, per column, on valid cells.
     */
    void apply_albedo (int lev, amrex::MultiFab& alb_sw, const ForceRestoreSnowParams& prm) const;

    /** @brief The snow water equivalent [kg/m^2], or nullptr. */
    [[nodiscard]] amrex::MultiFab* swe (int lev) noexcept {
        return defined(lev) ? m_swe[lev].get() : nullptr;
    }

    /** @brief The running totals (Total components) [kg/m^2], or nullptr. */
    [[nodiscard]] amrex::MultiFab* totals (int lev) noexcept {
        return defined(lev) ? m_totals[lev].get() : nullptr;
    }

    /** @brief Water the last step() handed the soil [kg/m^2]: melt plus sublimation; or nullptr. */
    [[nodiscard]] const amrex::MultiFab* soil_water (int lev) const noexcept {
        return defined(lev) ? m_soil_water[lev].get() : nullptr;
    }

    void write_checkpoint (int lev, const std::string& checkpointname) const;
    void read_checkpoint (int lev, const std::string& restart_chkfile,
                          const amrex::Periodicity& period);

private:
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_swe;          ///< W_n
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_totals;       ///< Total components
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_accum_start;  ///< frozen accumulation at the step start
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_frozen_rate;  ///< P_n of the last step
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_soil_water;   ///< melt + sublimation of the last step
    amrex::Vector<int> m_started;                                   ///< begin_step ran this step
    amrex::Vector<int> m_unused;                                    ///< end_step's snowfall not yet used by step()
    bool m_noted_no_frozen = false;                                 ///< the no-frozen-accumulator note was printed
    amrex::Vector<amrex::Real> m_elapsed;                           ///< time since m_accum_start [s]
};

#endif // ERF_FORCE_RESTORE_SNOW_STATE_H_

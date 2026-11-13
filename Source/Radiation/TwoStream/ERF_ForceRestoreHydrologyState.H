#ifndef ERF_FORCE_RESTORE_HYDROLOGY_STATE_H_
#define ERF_FORCE_RESTORE_HYDROLOGY_STATE_H_

#include <memory>
#include <string>

#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <ERF_ForceRestoreHydrology.H>
#include <ERF_SurfacePrecipitation.H>

struct RadChoice;

/**
 * @brief The force-restore soil of one level, read-only, for a model that reads the ground
 *        (the dust coupling's DustStep::soil). Every pointer is to a 2D field (i, j, k = 0)
 *        on the level's horizontal grids, valid for the duration of the call it is handed to.
 */
struct ForceRestoreSoil
{
    const amrex::MultiFab* q_sfc     = nullptr;  ///< surface-layer water content [m^3/m^3]
    const amrex::MultiFab* q_deep    = nullptr;  ///< column water content [m^3/m^3]; constant without the deep reservoir
    const amrex::MultiFab* t_sfc     = nullptr;  ///< surface temperature [K]
    const amrex::MultiFab* rain_rate = nullptr;  ///< rain of the last step [kg/m^2/s]; nullptr without seb_precipitation_source
    amrex::Real d_s      = 0.0;    ///< surface layer depth [m] (seb_moisture_layer_depth_m)
    amrex::Real d_2      = 0.0;    ///< column depth [m] (seb_deep_layer_depth_m); 0 without the deep reservoir
    amrex::Real porosity = -1.0;   ///< [m^3/m^3]; -1 when neither seb_soil_porosity nor seb_soil_type gives one
    amrex::Real fc       = 0.0;    ///< field capacity [m^3/m^3] (seb_soil_moisture_fc)
    amrex::Real wilt     = 0.0;    ///< wilting point [m^3/m^3] (seb_soil_moisture_wilt)
    bool deep = false;             ///< seb_deep_prognostic_enable
};

/**
 * @brief The per-level state of the force-restore hydrology (ERF_ForceRestoreHydrology.H):
 *        the rain rate of the current step, and the running totals that close each column's
 *        water budget.
 *
 * Owned by TwoStreamRadiation, which keeps the force-restore state itself (q_s, q_deep,
 * T_deep). Per level, and per step:
 *  - begin_step(...) at the start of ERF::Advance, before the dycore and microphysics,
 *    stores the surface precipitation accumulation of every column;
 *  - end_step(...) after microphysics turns its change over the step into the rain rate;
 *  - step(...) at the force-restore update after the dycore applies one hydrology step
 *    to every open column and adds to the running totals.
 * The rain rate needs no checkpoint or regrid handling: both captures happen inside one
 * ERF::Advance of the level. The totals are checkpointed, and TwoStreamRadiation lists
 * them with its prognostic fields for regrid, coarse-fine filling and averaging down.
 */
class ForceRestoreHydrology
{
public:
    /** @brief Components of totals(lev), each a running total [kg/m^2] since the start. */
    enum Total : int { Rain = 0, Runoff, Drainage, Deficit, NumTotals };

    /** @brief Domain sums of one step on one level (rank-reduced), for the log. */
    struct StepSums {
        amrex::Real rain = 0.0;          ///< sum of P dt [kg/m^2]
        amrex::Real evaporation = 0.0;   ///< sum of E dt [kg/m^2]
        amrex::Real runoff = 0.0;        ///< [kg/m^2]
        amrex::Real drainage = 0.0;      ///< [kg/m^2]
        amrex::Real deficit = 0.0;       ///< [kg/m^2]
        amrex::Real max_residual = 0.0;  ///< max |column budget residual| [kg/m^2]
        amrex::Long columns = 0;         ///< open columns updated
    };

    /** @brief The parameters of one step, from the erf.radiation.seb_* inputs. */
    [[nodiscard]] static ForceRestoreHydrologyParams params (const RadChoice& rc);

    void resize (int nlevs_max);

    /** @brief Allocate the level's fields; a no-op unless the hydrology is enabled. */
    void define_level (int lev, const RadChoice& rc, const amrex::BoxArray& ba2d,
                       const amrex::DistributionMapping& dm, const amrex::IntVect& ng);

    /** @brief True once define_level allocated @p lev. */
    [[nodiscard]] bool defined (int lev) const noexcept {
        return lev >= 0 && lev < static_cast<int>(m_totals.size()) && m_totals[lev] != nullptr;
    }

    /**
     * @brief Store the accumulation [kg/m^2] at the start of the level's step. When step()
     * did not use the rain of the previous step (the balance skips a new level's first step),
     * the earlier start is kept, so that rain is carried into the next step.
     */
    void begin_step (int lev, const SurfacePrecipAccumulationSources& sources, int klo);

    /**
     * @brief The rain rate [kg/m^2/s]: the accumulation now minus the one begin_step stored,
     * over the time since that start (the step just taken, plus any step whose rain was
     * carried). Zero when begin_step has not run for this step.
     */
    void end_step (int lev, const SurfacePrecipAccumulationSources& sources, int klo,
                   amrex::Real dt);

    /**
     * @brief One hydrology step on every open column of @p lev, after the balance has
     * updated T_s and q_s. A footprint cell (footprint > 0.5) keeps its state. The rain it
     * takes is the rain rate times the time since the stored start, so a carried step's rain
     * is counted once, in full.
     *
     * @param t_s       surface temperature after the update [K]
     * @param q_s       surface water content [m^3/m^3], updated
     * @param q_deep    column water content [m^3/m^3], updated with the deep reservoir
     * @param t_deep    deep temperature [K], updated with the deep reservoir
     * @param lh        latent heat flux the balance used [W/m^2]; E = LH / L_v
     * @param footprint building footprints, or nullptr
     * @param dt        the step [s]
     * @param prm       params(rc)
     * @return the step's domain sums, reduced over ranks
     */
    StepSums step (int lev, const amrex::MultiFab& t_s, amrex::MultiFab& q_s,
                   amrex::MultiFab& q_deep, amrex::MultiFab& t_deep,
                   const amrex::MultiFab& lh, const amrex::MultiFab* footprint,
                   amrex::Real dt, const ForceRestoreHydrologyParams& prm);

    /** @brief The rain rate of the last step [kg/m^2/s], or nullptr. Read-only. */
    [[nodiscard]] const amrex::MultiFab* rain_rate (int lev) const noexcept {
        return defined(lev) ? m_rain_rate[lev].get() : nullptr;
    }

    /** @brief The running totals (Total components) [kg/m^2], or nullptr. */
    [[nodiscard]] amrex::MultiFab* totals (int lev) noexcept {
        return defined(lev) ? m_totals[lev].get() : nullptr;
    }

    void write_checkpoint (int lev, const std::string& checkpointname) const;
    void read_checkpoint (int lev, const std::string& restart_chkfile,
                          const amrex::Periodicity& period);

private:
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_accum_start;  ///< accumulation at the step start
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_rain_rate;    ///< P of the last step
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_totals;       ///< Total components
    amrex::Vector<int> m_started;                                   ///< begin_step ran this step
    amrex::Vector<int> m_unused;                                    ///< end_step's rain not yet used by step()
    amrex::Vector<amrex::Real> m_elapsed;                           ///< time since m_accum_start [s]
};

#endif // ERF_FORCE_RESTORE_HYDROLOGY_STATE_H_

#ifndef ERF_FORCE_RESTORE_SNOW_STATE_H_
#define ERF_FORCE_RESTORE_SNOW_STATE_H_

#include <memory>
#include <string>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Periodicity.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <ERF_ForceRestorePrecipRate.H>
#include <ERF_ForceRestoreSnow.H>

struct RadChoice;

/**
 * @brief The per-level state of the force-restore snow layer (ERF_ForceRestoreSnow.H,
 *        erf.radiation.seb_snow_enable): the snow water equivalent and the running totals
 *        that close each column's snow budget.
 *
 * Owned by TwoStreamRadiation, beside the force-restore hydrology
 * (ERF_ForceRestoreHydrologyState.H) and the precipitation rate of the step
 * (ERF_ForceRestorePrecipRate.H). step(...) after the balance, and before the hydrology
 * step, applies one snow step to every open column: it takes the frozen precipitation,
 * changes T_s where snow melts and leaves the water for the soil in soil_water(lev), which
 * the hydrology step takes.
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
        /// open columns left with a non-finite T_s, or a non-finite or negative W_n
        amrex::Long bad_columns = 0;
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
     * @brief One snow step on every open column of @p lev, after the balance has updated
     * T_s. A footprint cell (footprint > 0.5) keeps its state and hands the soil no water;
     * so does a water column (land_mask 0), which the balance does not run on.
     *
     * @param t_s       surface temperature [K]; lowered where snow melts
     * @param lh        latent heat flux the balance used [W/m^2]; E = LH / L_v
     * @param footprint building footprints, or nullptr
     * @param dt        the step [s]
     * @param prm       params(rc)
     * @param snowfall  the frozen precipitation of the step (ForceRestorePrecipRate::frozen)
     * @param[out] t_s_sum_local sum of T_s over this rank's open columns after the step
     * @param[out] t_s_max_local max |T_s| over this rank's open columns after the step
     * @param land_mask the level's land mask (0 water), or nullptr: every column is land
     * @return the step's domain sums, reduced over ranks
     */
    StepSums step (int lev, amrex::MultiFab& t_s, const amrex::MultiFab& lh,
                   const amrex::MultiFab* footprint, amrex::Real dt,
                   const ForceRestoreSnowParams& prm, const ForceRestorePrecipInput& snowfall,
                   amrex::Real& t_s_sum_local, amrex::Real& t_s_max_local,
                   const amrex::iMultiFab* land_mask = nullptr);

    /**
     * @brief Write the surface albedo with snow into @p alb_sw, which holds the snow-free
     * albedo: DRM95 Eq. (8) with veg = 0, per column, on valid cells. A footprint cell
     * (@p footprint > 0.5, or none with nullptr) keeps the snow-free albedo: no snow lies on
     * a building, whatever W_n holds there before the first step. The same holds for a
     * water column (@p land_mask 0): step() and this function leave out the same columns.
     */
    void apply_albedo (int lev, amrex::MultiFab& alb_sw, const ForceRestoreSnowParams& prm,
                       const amrex::MultiFab* footprint = nullptr,
                       const amrex::iMultiFab* land_mask = nullptr) const;

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
    /**
     * @brief Read W_n and the totals. A checkpoint with neither file was written without
     * the snow layer: the snow starts from the input (noted once). One file without the
     * other is an incomplete checkpoint and stops the run, naming the missing file.
     */
    void read_checkpoint (int lev, const std::string& restart_chkfile,
                          const amrex::Periodicity& period);

private:
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_swe;          ///< W_n
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_totals;       ///< Total components
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_soil_water;   ///< melt + sublimation of the last step
};

#endif // ERF_FORCE_RESTORE_SNOW_STATE_H_

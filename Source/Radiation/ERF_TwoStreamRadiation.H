#ifndef ERF_TWO_STREAM_RADIATION_H_
#define ERF_TWO_STREAM_RADIATION_H_

#include <memory>
#include <string>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <ERF_RadStruct.H>
#include <ERF_LandSurface.H>
#include <ERF_RadiationDiagnostics.H>

/**
 * @brief The two-stream radiation model (erf.radiation_model = TwoStream).
 *
 * Owns everything the model keeps between steps: the per-level surface
 * property fields it falls back on when no land-surface model supplies
 * them, the surface-energy-balance (SEB) fields and force-restore state,
 * and the domain-mean flux diagnostics of the last column sweep. ERF holds
 * one instance and calls it only after testing
 * `solverChoice.rad_type == RadiationType::TwoStream`; every method
 * here checks the same predicate and returns without acting otherwise.
 *
 * Lifecycle:
 *  - resize(nlevs) once at construction (before the inputs are read);
 *  - define_level(lev, radChoice, ba2d, dm) whenever a level is (re)built;
 *  - advance(...) twice per step: "pre_dycore" runs the column sweep on the
 *    old state and writes the heating rates, "post_dycore" reuses that
 *    sweep's flux diagnostics and advances the SEB state by dt_step;
 *  - write_checkpoint / read_checkpoint for the prognostic surface state.
 *
 * The heating-rate MultiFab itself stays in ERF, since RRTMGP writes the
 * same 2-component (SW, LW) array and the RhoTheta source term and the
 * plotfile read it without caring which model filled it.
 */
class TwoStreamRadiation
{
public:
    TwoStreamRadiation () = default;

    /** @brief Size the per-level containers; no allocation happens here. */
    void resize (int nlevs_max);

    /** @brief True once a level has been defined with the model selected. */
    [[nodiscard]] bool active () const {
        return (m_rad != nullptr) && m_rad->enabled;
    }

    /**
     * @brief Allocate the 2D surface fields of one level on the given
     * horizontal BoxArray and fill them with the scalar defaults. A no-op
     * unless rad_choice.rad_type == TwoStream.
     */
    void define_level (int lev,
                       const RadChoice& rad_choice,
                       const amrex::BoxArray& ba2d,
                       const amrex::DistributionMapping& dm);

    /**
     * @brief One radiation call at one level.
     *
     * @param lev        level
     * @param nstep      step number, used as the diagnostics row index
     * @param time       time logged with the diagnostics row [s]
     * @param dt_step    step size, used by the SEB Euler update [s]
     * @param call_site  "pre_dycore" (column sweep + diagnostics) or
     *                   "post_dycore" (cached diagnostics + SEB update)
     * @param cons_old   old-state conserved variables of this level
     * @param z_phys_nd  nodal heights (layer interfaces), or nullptr on a uniform grid
     * @param geom       geometry of this level
     * @param lsm        land-surface model (queried for per-column fields)
     * @param qheating   2-component (SW, LW) heating-rate MultiFab of this
     *                   level, written at the pre-dycore call; may be nullptr
     */
    void advance (int lev,
                  int nstep,
                  amrex::Real time,
                  amrex::Real dt_step,
                  const std::string& call_site,
                  const amrex::MultiFab& cons_old,
                  const amrex::MultiFab* z_phys_nd,
                  const amrex::Geometry& geom,
                  LandSurface& lsm,
                  amrex::MultiFab* qheating);

    /** @brief Write the prognostic surface state of one level (SEB on). */
    void write_checkpoint (int lev, const std::string& checkpointname) const;

    /** @brief Read it back when the checkpoint carries it; else keep defaults. */
    void read_checkpoint (int lev, const std::string& restart_chkfile);

private:
    // Domain-mean flux diagnostics of the last column sweep, per level, so the
    // post-dycore call can report them without sweeping again.
    struct FluxDiag {
        amrex::Real SW_surface = 0.0;
        amrex::Real SW_TOA = 0.0;
        amrex::Real SW_up_TOA = 0.0;
        amrex::Real LW_net_surface = 0.0;
        amrex::Real LW_up_TOA = 0.0;
        amrex::Real heating_rate_max = 0.0;
    };

    const RadChoice* m_rad = nullptr;

    // CSV / console diagnostics writer, created on the first advance() so its
    // header and duplicate-row state persist across calls.
    std::unique_ptr<RadiationDiagnostics> m_diag;

    // Fallback surface properties when no LSM supplies them (2D: i,j).
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_alb_sw;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_emiss_lw;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_t_sfc;      // prognostic when the SEB evolves it

    // Surface-energy-balance fields (2D: i,j), allocated when seb_enable.
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_sw_flux_sfc; // net absorbed shortwave [W/m^2]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_lw_flux_sfc; // net absorbed longwave [W/m^2]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_hfx_sfc;     // sensible heat flux [W/m^2]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_lh_sfc;      // latent heat flux [W/m^2]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_grdflx_sfc;  // ground heat flux [W/m^2]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_q_sfc;       // surface moisture [kg/kg], prognostic
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_t_deep;      // deep-soil temperature [K]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_q_deep;      // deep-soil moisture [kg/kg]

    amrex::Vector<FluxDiag> m_flux_diag;
};

#endif // ERF_TWO_STREAM_RADIATION_H_

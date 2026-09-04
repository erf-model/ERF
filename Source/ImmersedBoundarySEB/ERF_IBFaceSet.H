#ifndef ERF_IB_FACE_SET_H
#define ERF_IB_FACE_SET_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <vector>
#include "ERF_IBSEBParams.H"
#include "ERF_IBSEBMaterials.H"

/**
 * \file ERF_IBFaceSet.H
 * \brief Compact list of the wall faces of resolved buildings on one AMR level.
 *
 * **What a face is.** The immersed forcing marks every cell of a level as
 * fluid or solid through the cell-centred blanking ``terrain_blanking``
 * (one minus the embedded-boundary volume fraction, so 0 in the air and 1
 * inside a building). A wall face is the cell face between a fluid cell
 * (blanking below 0.5) and a solid cell (blanking at or above 0.5). Walls
 * are x or y faces, roofs are z faces with the fluid above the solid.
 *
 * **Face convention.** A face is identified by its *fluid* cell ``(i, j, k)``,
 * its direction ``dir`` (0 = x, 1 = y, 2 = z) and its ``side`` (-1 when the
 * solid neighbour is at ``i-1`` / ``j-1`` / ``k-1``, +1 when it is at
 * ``i+1`` / ``j+1`` / ``k+1``). The outward unit normal of the wall, pointing
 * into the fluid, is therefore ``-side`` along ``dir``. The fluid cell is where
 * every atmospheric input of the balance (wind, temperature, humidity,
 * radiation) is read and where the wall's heat flux is deposited.
 *
 * **Ownership and layout.** Each face is stored exactly once, on the rank
 * that owns its fluid cell, so the per-step update of the balance is a kernel
 * over this rank's faces with no communication; only the reports reduce
 * across ranks. The faces are a struct of device arrays (one
 * ``amrex::Gpu::DeviceVector`` per quantity, all indexed by the face number)
 * and are stored contiguously per local fab in MFIter order:
 * ``fab_start()[mfi.LocalIndex()]`` to ``fab_start()[mfi.LocalIndex()+1]``
 * are the faces whose fluid cell is in that fab, which is what lets a kernel
 * scatter into that fab's arrays without a search.
 *
 * **Building ids.** The solid columns of the blanking (any solid cell in the
 * column) are labelled by 4-connectivity in the horizontal and numbered in
 * scan order over ``(i, j)``, identically on every rank from a reduced column
 * mask; a face carries the id of the column of its solid neighbour. Ids are
 * 1-based, 0 is never assigned.
 *
 * **Output.** The list is never written directly. scatter_diagnostics() puts
 * the number of faces and their mean skin temperature into cell-centred
 * fields for the plotfile; save_state() / load_state() move the skin and
 * slab temperatures through a cell-centred field with six slots per cell
 * (one per direction and side) for the checkpoint. On restart the list is
 * rebuilt from the blanking and refilled from that field, so a restart does
 * not depend on the rank count or on the order the faces were found in.
 *
 * **Phases and the per-step order.** Phase 1 stores; phase 3 fills
 * ``d_f_sky`` / ``d_f_ground`` / ``d_f_bldg`` once; every step then runs
 * compute_shortwave() (phase 2), compute_longwave() (phase 3),
 * compute_sensible() (phase 4) and, with ``erf.ibseb.prognostic``,
 * solve_balance() (phase 6), which finds the skin temperature closing the
 * balance, advances the slab (phase 5) with it and rewrites the fluxes at the
 * solution; without it compute_ground() advances the slab under the fixed
 * skin. ``d_Q_ext`` is the hook through which an external incident flux (a
 * fire's radiation) enters without the balance knowing its source: whoever
 * owns it writes the array before solve_balance() runs. See
 * IBSEB_DEVELOPMENT.md.
 */
class IBFaceSet
{
public:
    /**
     * Construct an empty face set for one level; build() fills it.
     * @param params Inputs of the balance (copied; the set never reads ParmParse).
     * @param lev    AMR level the set belongs to, for reports only.
     */
    IBFaceSet (const IBSEBParams& params, int lev) : m_params(params), m_lev(lev) {}

    /**
     * Detect this rank's faces from the cell-centred blanking and allocate
     * every per-face array.
     *
     * Runs on the host once at initialisation: the blanking is copied fab by
     * fab (a pinned copy on GPU builds), the valid cells are scanned, and the
     * arrays are uploaded to the device at the end. Must be called after the
     * blanking has been built and its ghost cells filled, since the solid
     * neighbour of a fluid cell on a box edge lives in a ghost cell.
     *
     * @param blanking Cell-centred blanking of the level (needs at least one
     *                 ghost cell, filled).
     * @param geom     Geometry of the level: cell sizes for the face areas,
     *                 domain and periodicity for the neighbour tests.
     */
    void build (const amrex::MultiFab& blanking, const amrex::Geometry& geom);

    /** Number of faces owned by this rank (not the global count; see report()). */
    int n_faces () const { return m_nface; }
    /** Number of buildings of the level, the same on every rank. */
    int n_buildings () const { return m_nbld; }
    /** Conduction layers per face, from ``erf.ibseb.n_slab_layers``. */
    int n_layers () const { return m_params.n_slab_layers; }

    /**
     * Number of components of the checkpoint state field.
     *
     * A fluid cell can carry up to six faces (two per direction), so the field
     * has six slots of ``1 + n_layers()`` values: component
     * ``(dir*2 + (side>0)) * (1 + n_layers())`` is the skin temperature of that
     * slot and the following ``n_layers()`` components are its slab layers.
     */
    int state_ncomp () const { return 6 * (1 + n_layers()); }

    /**
     * View fractions (phase 3): sample a cosine-weighted hemisphere around
     * every face's outward normal with ``view_n_az x view_n_el`` rays through
     * the column walk and count where they end: sky, ground or a building.
     * The three fractions sum to one and are view factors. Called once at
     * initialisation; the geometry is static.
     */
    void compute_view_fractions ();

    /**
     * Longwave of the current step (phase 3). Incoming on a face:
     * ``f_sky LW_sky + f_ground eps_g sigma T_g^4 + f_bldg sigma T_skin^4``,
     * the sky term fixed or gray (``sky_emissivity sigma T_air^4`` with the
     * air temperature of the face's fluid cell), the building term the
     * isothermal-surroundings approximation. Net:
     * ``eps (LW_in - sigma T_skin^4)``, positive into the face.
     * @param cons Conserved state of the level, for the air temperature.
     */
    void compute_longwave (const amrex::MultiFab& cons);

    /**
     * Assign every face its material properties (phase 5): from the material
     * library by building id, or the uniform inputs. Called once after
     * build(); fills d_mat, d_albedo, d_emis, d_k, d_rhocp, d_thick and
     * initialises the slab layers to the skin temperature.
     */
    void assign_materials ();

    /**
     * Ground heat flux of the current step (phase 5): advance every face's
     * slab by dt with the skin temperature at the top and the interior
     * temperature at the bottom (implicit, ERF_IBSEBSlab.H) and store the
     * conduction into the slab, G, positive into the wall.
     */
    void compute_ground (amrex::Real dt);

    /**
     * Prognostic balance of the current step (phase 6). Per face: the
     * slab's linear response to the skin temperature for this step, the
     * Newton solve of ERF_IBSEBBalance.H with the shortwave, the external
     * flux, the sky and ground longwave, the frozen wall-function
     * coefficient and the latent flux, the slab advanced with the solution,
     * and the fluxes (LW_in, LW_net, H, G) rewritten at the solution. Fills
     * d_resid (|balance| at the end, non-zero only at a bound or the
     * iteration cap) and d_niter. Needs compute_shortwave(),
     * compute_longwave() and compute_sensible() of the same step first.
     * @param dt Atmospheric time step [s], the slab's step.
     */
    void solve_balance (amrex::Real dt);

    /**
     * Sensible heat of the current step (phase 4) through a wall function on
     * every face: the tangential wind of the fluid cell at half a cell from
     * the wall gives u* with the roughness z0_wall, the skin-to-air
     * potential-temperature difference gives theta* with z0h_wall, and
     * H = rho c_p u* theta*, positive out of the face. Neutral, or with the
     * surface layer's stability functions on roofs when asked. The latent
     * flux stays zero in this phase.
     *
     * Fills d_H, d_ustar, d_U_tan, d_rho, d_theta_air and d_H_coeff (the
     * coefficient of the skin-to-air potential-temperature difference that
     * the balance keeps frozen); T_air was read by compute_longwave(), which
     * must run first in a step.
     * @param cons  Conserved state of the level (density, rho theta).
     * @param xvel, yvel, zvel  Face-centred velocities of the level.
     * @param c_p   Specific heat of air [J/kg/K].
     */
    void compute_sensible (const amrex::MultiFab& cons, const amrex::MultiFab& xvel,
                           const amrex::MultiFab& yvel, const amrex::MultiFab& zvel, amrex::Real c_p);

    /**
     * Add the face heat flux to the temperature equation: every face
     * deposits ``H A / (c_p V Pi)`` into the rho-theta source of its fluid
     * cell (Pi the Exner function, so the flux heats the cell by H A / (c_p V)
     * in temperature). Called after the sources are rebuilt at every slow
     * stage; adds, never overwrites, and uses atomic adds because a corner
     * cell receives several faces.
     */
    void add_heat_flux_to_source (amrex::MultiFab& source, const amrex::MultiFab& cons,
                                  const amrex::Geometry& geom, amrex::Real c_p, amrex::Real rdOcp) const;

    /**
     * Shortwave of the current step (phase 2): sun position and irradiances
     * from the prescribed provider, then per face the direct beam on the
     * face (zero when the ray toward the sun hits a building), the diffuse
     * light through the view fractions, and the absorbed sum.
     *
     * Fills d_shadow, d_SW_direct_in, d_SW_diffuse_in and d_SW_abs; stores
     * the sun vector and irradiances for the reports.
     * @param time Simulation time [s]; with sun_mode = solar it is added to
     *             time_zero_utc_s.
     */
    void compute_shortwave (amrex::Real time);

    /**
     * Scatter one per-face array into a cell-centred field as the mean over
     * the faces touching each cell (zero where there are none).
     */
    void scatter_field (const amrex::Gpu::DeviceVector<amrex::Real>& v, amrex::MultiFab& out) const;

    /**
     * Write every face of this rank to ``<prefix>.rank<N>.csv`` with its
     * geometry, view fractions, shadow flag, shortwave terms and skin
     * temperature, for the regtests and for users who want the raw list.
     */
    void dump_faces (const std::string& prefix) const;

    /**
     * Scatter the faces into two cell-centred diagnostics for the plotfile.
     *
     * @param[out] nfaces Number of wall faces touching each fluid cell (0 away
     *                    from buildings; up to 3 at an outside corner, 6 in a
     *                    one-cell slot).
     * @param[out] tskin  Mean skin temperature of those faces, 0 where there
     *                    are none.
     *
     * Both fields must be on the level's box array with no ghost cells. Faces
     * of one fab are accumulated with atomic adds because a corner cell
     * receives several faces.
     */
    void scatter_diagnostics (amrex::MultiFab& nfaces, amrex::MultiFab& tskin) const;

    /**
     * Write the face state (skin and slab temperatures) into a cell-centred
     * field with state_ncomp() components, for the checkpoint.
     *
     * Each face writes its own slot (see state_ncomp()), so no two faces
     * touch the same component of the same cell and no atomics are needed.
     * Cells without a face in a slot keep zero there.
     */
    void save_state (amrex::MultiFab& state) const;

    /**
     * Read the face state back from a field written by save_state(), after
     * build() has recreated the list from the blanking. Only the slots that
     * correspond to a face of this list are read.
     */
    void load_state (const amrex::MultiFab& state);

    /**
     * Print the ``[IBSEB]`` summary line and, when asked, append one row per
     * building to ``erf.ibseb.csv_file``.
     *
     * The summary carries the global face counts per direction, the number of
     * buildings, the total face area and the range of the skin temperature;
     * the CSV rows carry per building the face count, the area and the
     * area-weighted mean skin temperature. All numbers are reduced over the
     * ranks; the file is written by the I/O rank only, with a header only
     * when the file does not exist yet, so a restarted run appends.
     *
     * @param time      Simulation time [s] written to the rows.
     * @param step      Step number written to the rows.
     * @param write_csv Whether to append rows, or only print the summary.
     */
    void report (amrex::Real time, int step, bool write_csv) const;

    // ---- Per-face device arrays, indexed by face number ----------------------
    // Static geometry, set by build():
    amrex::Gpu::DeviceVector<int>         d_i, d_j, d_k;   ///< Fluid cell of the face
    amrex::Gpu::DeviceVector<int>         d_dir;           ///< Face direction: 0 x, 1 y, 2 z
    amrex::Gpu::DeviceVector<int>         d_side;          ///< Solid neighbour at -1 (low) or +1 (high) side
    amrex::Gpu::DeviceVector<int>         d_bid;           ///< Building id (1-based)
    amrex::Gpu::DeviceVector<int>         d_mat;           ///< Material id (0 = uniform inputs)
    amrex::Gpu::DeviceVector<amrex::Real> d_albedo, d_emis;///< Optical properties of the face
    amrex::Gpu::DeviceVector<amrex::Real> d_kth, d_rhocp, d_thick; ///< Slab conductivity [W/m/K], heat capacity [J/m3/K], thickness [m]
    amrex::Gpu::DeviceVector<amrex::Real> d_area;          ///< Face area [m2]
    amrex::Gpu::DeviceVector<amrex::Real> d_xf, d_yf, d_zf;///< Face centre [m]
    // State, advanced by the balance (phase 6) and checkpointed:
    amrex::Gpu::DeviceVector<amrex::Real> d_T_skin;        ///< Skin temperature [K]
    amrex::Gpu::DeviceVector<amrex::Real> d_T_slab;        ///< Slab layers [K], face-major: [face * n_layers + layer]
    // Static radiative geometry (phase 3):
    amrex::Gpu::DeviceVector<amrex::Real> d_f_sky, d_f_ground, d_f_bldg;  ///< View fractions, sum to 1
    // Fluxes of the current step [W/m2], positive as in IBSEB_DEVELOPMENT.md
    // (SW_abs and LW_net into the face, H, LE, G out of it):
    amrex::Gpu::DeviceVector<amrex::Real> d_SW_abs, d_LW_net, d_H, d_LE, d_G;
    amrex::Gpu::DeviceVector<amrex::Real> d_Q_ext;         ///< External incident flux [W/m2], e.g. fire radiation
    // Shortwave detail (phase 2), all [W/m2] except the flag:
    amrex::Gpu::DeviceVector<amrex::Real> d_shadow;        ///< 1 when the ray to the sun hits a building, else 0
    amrex::Gpu::DeviceVector<amrex::Real> d_SW_direct_in;  ///< Direct beam incident on the face
    amrex::Gpu::DeviceVector<amrex::Real> d_SW_diffuse_in; ///< Sky and ground-reflected diffuse incident
    // Longwave detail (phase 3) [W/m2] and the air temperature read [K]:
    amrex::Gpu::DeviceVector<amrex::Real> d_LW_down_in;   ///< Incoming longwave on the face (sky, ground and walls)
    amrex::Gpu::DeviceVector<amrex::Real> d_T_air;        ///< Air temperature of the fluid cell
    // Sensible heat detail (phase 4):
    amrex::Gpu::DeviceVector<amrex::Real> d_theta_air;    ///< Potential temperature of the fluid cell [K]
    amrex::Gpu::DeviceVector<amrex::Real> d_rho;          ///< Density of the fluid cell [kg/m3]
    amrex::Gpu::DeviceVector<amrex::Real> d_U_tan;        ///< Tangential wind speed at the fluid cell [m/s]
    amrex::Gpu::DeviceVector<amrex::Real> d_ustar;        ///< Friction velocity of the wall function [m/s]
    amrex::Gpu::DeviceVector<amrex::Real> d_H_coeff;      ///< H = H_coeff (theta_skin - theta_air) [W/m2/K]
    // Prognostic balance detail (phase 6):
    amrex::Gpu::DeviceVector<amrex::Real> d_LW_ext;       ///< Incoming longwave from sky and ground only [W/m2]
    amrex::Gpu::DeviceVector<amrex::Real> d_resid;        ///< |balance| after the solve [W/m2]
    amrex::Gpu::DeviceVector<int>         d_niter;        ///< Newton iterations of the last solve
    // Column tops of the whole level [m], replicated on every rank for the
    // ray cast: index ci * ny + cj over the domain's columns.
    amrex::Gpu::DeviceVector<amrex::Real> d_col_top;

    /// Sun and irradiances of the last compute_shortwave() call, for reports.
    struct SunState {
        amrex::Real zenith = 0.0, azimuth = 0.0;        ///< [rad]
        amrex::Real sx = 0.0, sy = 0.0, sz = 1.0;       ///< Unit vector toward the sun
        amrex::Real dni = 0.0;                          ///< Direct-normal irradiance [W/m2]
        amrex::Real diffuse_h = 0.0;                    ///< Diffuse on a horizontal surface [W/m2]
    };
    const SunState& sun () const { return m_sun; }

    /**
     * Print the ``[IBSEB DEBUG]`` description of the set: the inputs, every
     * rank's face count and per-fab ranges, the buildings with their footprint
     * columns and bounding boxes, and the device memory of the arrays.
     * Called once from build() when ``erf.ibseb.debug`` is set.
     */
    void print_debug_summary () const;

    /**
     * First face of each local fab in MFIter order; size number of fabs + 1.
     * Faces ``fab_start()[n]`` to ``fab_start()[n+1]-1`` belong to the fab with
     * local index ``n``, which is what the scatter kernels rely on.
     */
    const std::vector<int>& fab_start () const { return m_fab_start; }

private:
    IBSEBParams m_params;          ///< Copy of the inputs
    int m_lev   = 0;               ///< AMR level, for reports
    int m_nface = 0;               ///< Faces owned by this rank
    int m_nbld  = 0;               ///< Buildings of the level (global)
    std::vector<int> m_fab_start;  ///< See fab_start()

    // Static totals over all ranks, computed once in build(); building index 0 unused.
    std::vector<long>        m_bld_nface;      ///< Faces per building
    std::vector<amrex::Real> m_bld_area;       ///< Face area per building [m2]
    long m_nface_dir[3] = {0, 0, 0};           ///< Faces per direction
    amrex::Real m_area_total = 0.0;            ///< Total face area [m2]
    // Footprint of each building (global): columns, and bounding box in cell indices.
    std::vector<long> m_bld_ncol;
    std::vector<int>  m_bld_ilo, m_bld_ihi, m_bld_jlo, m_bld_jhi;

    // Geometry kept for the ray cast and the reports.
    int  m_nx = 0, m_ny = 0;                   ///< Columns of the domain
    amrex::Real m_x_lo = 0.0, m_y_lo = 0.0;    ///< Domain origin [m]
    amrex::Real m_dx[3] = {0.0, 0.0, 0.0};     ///< Cell sizes [m]
    bool m_per_x = false, m_per_y = false;     ///< Periodicity in x and y
    amrex::Real m_col_top_max = 0.0;           ///< Tallest column [m]
    amrex::Real m_z_ground = 0.0;              ///< Ground height [m] (flat, the domain bottom)
    amrex::Real m_max_path = 0.0;              ///< Ray path cap [m]
    SunState m_sun;
};

#endif

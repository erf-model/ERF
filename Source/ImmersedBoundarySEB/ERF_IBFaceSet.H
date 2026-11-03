#ifndef ERF_IB_FACE_SET_H
#define ERF_IB_FACE_SET_H

#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_iMultiFab.H>
#include <cmath>
#include <string>
#include <vector>
#include "ERF_IBSEBParams.H"
#include "ERF_IBSEBMaterials.H"
#include "ERF_IBSEBSolar.H"
#include <ERF_TwoStreamCanopyForcing.H>

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
 * **Column maps.** The solid columns of the blanking (any solid cell in the
 * column) are kept as a k-collapsed ``amrex::iMultiFab`` on the level's own
 * boxes, one integer per column: the domain-relative k of its highest solid
 * cell, -1 where it is fluid. Nothing about the columns is replicated across
 * the ranks (erf-model/ERF#4020): the building labelling exchanges labels
 * through the map's ghost cells, and the ray cast of the shortwave and the
 * view fractions reads a per-rank *patch* of the map, this rank's boxes
 * grown by a halo that holds every column a ray from one of its faces can
 * reach (build_ray_map()).
 *
 * **Building ids.** The solid columns are labelled by 4-connectivity in the
 * horizontal (across the periodic seams too) and numbered in scan order over
 * ``(i, j)``, the same on every rank and whatever the decomposition: every
 * column starts with its own scan index as its label, each box takes the
 * lowest label of its connected pieces, and the pieces exchange labels
 * through one ghost cell until nothing changes, so a building ends with the
 * scan index of its first column; the distinct labels, gathered once, give
 * the numbering. A face carries the id of the column of its solid neighbour.
 * Ids are 1-based, 0 is never assigned.
 *
 * **Output.** The list is never written directly. scatter_diagnostics() puts
 * the number of faces and their mean skin temperature into cell-centred
 * fields for the plotfile; write_state() / read_state() move the skin
 * temperature, the sensible flux and the slab temperatures through one
 * record per face, keyed by the face's cell, direction and side, for the
 * checkpoint. On restart the list is rebuilt from the blanking and every
 * face finds its record by its key, so a restart does not depend on the
 * rank count or on the order the faces were found in.
 *
 * **Per-step order.** build() detects and stores the faces and
 * compute_view_fractions() fills ``d_f_sky`` / ``d_f_ground`` / ``d_f_bldg``
 * once at initialisation; every step then runs compute_shortwave(),
 * compute_longwave(), compute_sensible() and, with ``erf.ibseb.prognostic``,
 * solve_balance(), which finds the skin temperature closing the balance,
 * advances the slab with it and rewrites the fluxes at the solution; without
 * it compute_ground() advances the slab under the fixed skin. ``d_Q_ext`` is
 * the hook through which an external incident flux (a fire's radiation)
 * enters without the balance knowing its source: whoever owns it writes the
 * array before solve_balance() runs.
 */
class IBFaceSet
{
public:
    /**
     * Construct an empty face set for one level; build() fills it.
     * @param params Inputs of the balance (copied; the set never reads ParmParse).
     * @param lev    AMR level the set belongs to, for reports only.
     * @param qv_comp State component of rho qv (SolverChoice::moisture_indices.qv),
     *                negative without a moisture model. The air temperature and the
     *                Exner function of a fluid cell take its vapour through the
     *                equation of state, as the rest of ERF does; with no moisture
     *                model they take none.
     */
    IBFaceSet (const IBSEBParams& params, int lev, int qv_comp = -1)
        : m_params(params), m_lev(lev), m_qv_comp(qv_comp) {}

    /**
     * Detect this rank's faces from the cell-centred blanking and allocate
     * every per-face array.
     *
     * Runs on the host once at initialisation: the blanking is copied fab by
     * fab (a pinned copy on GPU builds), the valid cells are scanned, and the
     * arrays are uploaded to the device at the end. Must be called after the
     * blanking has been built and its ghost cells filled, since the solid
     * neighbour of a fluid cell on a box edge lives in a ghost cell. Stops
     * when a solid cell outside the level's grids neighbours a cell of the
     * level, fluid or solid: a building across the edge of a refined level.
     *
     * @param blanking Cell-centred blanking of the level (needs at least one
     *                 ghost cell, filled).
     * @param geom     Geometry of the level: cell sizes for the face areas,
     *                 domain and periodicity for the neighbour tests.
     */
    void build (const amrex::MultiFab& blanking, const amrex::Geometry& geom);

    /**
     * Add the buildings this refined level does not hold to the ray map of
     * its ray cast, from the column maps of the levels below.
     *
     * build() takes the map from the level's own cells, so on a refined level
     * that covers part of the domain the rays of the shortwave and the view
     * fractions would pass through every building outside it. The ray map
     * is rebuilt with the halo for the tallest column of every level, and
     * each column of the patch that this level's grids do not cover takes
     * the top of the finest level below that covers it
     * (ERF::ibseb_check_refined_levels() has made sure a building is wholly
     * inside or wholly outside each level), on the ratio x ratio columns of
     * this level beneath the coarse one, at the same height. Only the ray
     * cast reads the result: the building ids, the face heights of the wall
     * function and the reports keep the level's own buildings. Call after
     * build() of every level involved and before compute_view_fractions().
     *
     * @param below Face sets of levels 0 to this level's minus one, in that order.
     */
    void add_outside_occluders (const std::vector<const IBFaceSet*>& below);

    /**
     * Give the sun of the next compute_shortwave() with
     * ``erf.ibseb.sun_mode = two_stream``: the two-stream radiation's
     * declination and hour angle at that time, its site latitude and its
     * top-of-atmosphere irradiance (the date's distance factor included).
     * ERF::ibseb_set_two_stream_sun() calls it before every call; the call
     * uses it up, so a shortwave without a fresh sun stops rather than reusing
     * an old one.
     */
    void set_two_stream_sun (amrex::Real declination, amrex::Real hour_angle,
                             amrex::Real latitude_deg, amrex::Real toa_irradiance)
    {
        m_ts_sun = {true, declination, hour_angle, latitude_deg, toa_irradiance};
    }

    /**
     * Number this refined level's buildings as level 0 numbers them.
     *
     * Each level labels its buildings from its own cells, so the ids of a
     * refined level differ from level 0's wherever the levels see different
     * buildings: a building outside the refined level, or two that level 0
     * merges into one. Inputs that go by building (erf.ibseb.material_by_building)
     * and the report's ``building_level0`` column use level 0's numbering, so
     * every building of this level is looked up in the labels of the level
     * below, and through that level's own map in level 0's: each of its
     * columns votes with the id under it and the majority wins. A building
     * only this level resolves maps to 0 (the default material), with a
     * warning.
     * Call after build() of both levels; the level below keeps its labels until
     * release_labels().
     *
     * @param coarser Face set of the level below, already mapped.
     * @param ratio   Refinement ratio from the level below to this one.
     */
    void map_buildings_to_level0 (const IBFaceSet& coarser, const amrex::IntVect& ratio);

    /** Level 0's number for building ``b`` of this level (1-based; 0 for one only this level resolves). */
    int building_level0 (int b) const { return m_to_lev0.at(b); }

    /** Free the column labels once the level above has been mapped (map_buildings_to_level0()). */
    void release_labels () { m_lab.clear(); }

    /**
     * The ray map of this rank for the kernels of the shortwave and the view
     * fractions: the column tops over this rank's patch and the grid they are
     * drawn on (ibseb::RayMap). A rank without boxes has no patch and no faces.
     */
    ibseb::RayMap ray_map () const;

    /** Halo of the ray map in cells (x, y), from the last build_ray_map(). */
    amrex::IntVect ray_halo () const { return m_halo; }

    /** Whether every rank's patch spans the whole domain, so every ray is followed to its end. */
    bool ray_map_spans_domain () const { return m_map_spans_domain; }

    /** Number of faces owned by this rank (not the global count; see report()). */
    int n_faces () const { return m_nface; }
    /** Number of buildings of the level, the same on every rank. */
    int n_buildings () const { return m_nbld; }
    /** Conduction layers per face, from ``erf.ibseb.n_slab_layers``. */
    int n_layers () const { return m_params.n_slab_layers; }

    /**
     * View fractions: sample a cosine-weighted hemisphere around
     * every face's outward normal with ``view_n_az x view_n_el`` rays through
     * the column walk and count where they end: sky, ground or a building.
     * The three fractions sum to one and are view factors. Called once at
     * initialisation; the geometry is static.
     */
    void compute_view_fractions ();

    /**
     * Longwave of the current step. Incoming on a face:
     * ``f_sky LW_sky + f_ground LW_ground + f_bldg sigma T_skin^4``,
     * the sky term fixed or gray (``sky_emissivity sigma T_air^4`` with the
     * air temperature of the face's fluid cell, its water vapour included)
     * with the ground term
     * ``eps_g sigma T_g^4``, or with ``lw_mode = two_stream`` the longwave
     * down and up at the face's height in the face's column
     * (``canopy``); the building term is the isothermal-surroundings
     * approximation. Net: ``eps (LW_in - sigma T_skin^4)``, positive into the face.
     * @param cons   Conserved state of the level, for the air temperature.
     * @param canopy With ``lw_mode = two_stream``, this step's two-stream view of the
     *               level (ERF_TwoStreamCanopyForcing.H), or an empty one before the
     *               first sweep, when the incoming and net longwave are zero.
     */
    void compute_longwave (const amrex::MultiFab& cons, const TwoStreamCanopyView& canopy = TwoStreamCanopyView{});

    /**
     * Assign every face its material properties: from the material
     * library by building id, or the uniform inputs. Called once after
     * build(); fills d_mat, d_albedo, d_emis, d_k, d_rhocp, d_thick and
     * initialises the slab layers to the skin temperature.
     */
    void assign_materials ();

    /**
     * Ground heat flux of the current step: advance every face's
     * slab by dt with the skin temperature at the top and the interior
     * temperature at the bottom (implicit, ERF_IBSEBSlab.H) and store the
     * conduction into the slab, G, positive into the wall.
     */
    void compute_ground (amrex::Real dt);

    /**
     * Prognostic balance of the current step. Per face: the
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
     * Sensible heat of the current step through a wall function on
     * every face: the tangential wind of the fluid cell at half a cell from
     * the wall gives u* with the roughness z0_wall, the skin-to-air
     * potential-temperature difference gives theta* with z0h_wall, and
     * H = rho c_p u* theta*, positive out of the face. Neutral, or with the
     * surface layer's stability functions on roofs when asked. The latent
     * flux is not modelled and stays zero.
     *
     * Beyond neutral: with ``convective_velocity = deardorff`` the
     * wind is ``sqrt(U_tan^2 + (beta w*)^2)`` with w* from the previous
     * step's H (zero when H is into the face) and a depth that is the mixed
     * layer above a roof or the building height for a wall; with
     * ``stability_correction`` the roofs iterate the similarity functions on
     * the face's own Obukhov length, seeded from the ground surface layer's
     * 2D field at the face's column when given.
     *
     * Fills d_H, d_ustar, d_U_tan, d_rho, d_theta_air, d_H_coeff (the
     * coefficient of the skin-to-air potential-temperature difference that
     * the balance keeps frozen), d_w_star, d_olen and d_z_i; T_air was read
     * by compute_longwave(), which must run first in a step.
     * @param cons  Conserved state of the level (density, rho theta).
     * @param xvel, yvel, zvel  Face-centred velocities of the level.
     * @param c_p   Specific heat of air [J/kg/K].
     * @param olen_ground  Ground surface layer's Obukhov length (2D field on
     *                     the level's box array collapsed to k = 0), or null.
     * @param pblh_ground  Its boundary-layer height field, or null.
     * @param z_i_bulk     Mixed-layer depth of the level from the bulk
     *                     Richardson diagnostic or the fixed input [m].
     */
    void compute_sensible (const amrex::MultiFab& cons, const amrex::MultiFab& xvel,
                           const amrex::MultiFab& yvel, const amrex::MultiFab& zvel, amrex::Real c_p,
                           const amrex::MultiFab* olen_ground = nullptr,
                           const amrex::MultiFab* pblh_ground = nullptr,
                           amrex::Real z_i_bulk = 0.0);

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
     * Shortwave of the current step: the sun's position (sun_mode), the
     * irradiances of the prescribed provider or, with ``radiation =
     * two_stream``, those of the face's own column at its height (``canopy``),
     * then per face
     * the direct beam on the face (zero when the ray toward the sun hits a
     * building), the diffuse light through the view fractions, and the
     * absorbed sum.
     *
     * Fills d_shadow, d_SW_direct_in, d_SW_diffuse_in and d_SW_abs; stores
     * the sun vector and irradiances for the reports (with two_stream the
     * means over all the level's faces of the direct-normal and diffuse
     * irradiance each takes, written to every row of the level).
     * @param time   Simulation time [s]; with sun_mode = solar it is added to
     *               time_zero_utc_s (with sun_mode = two_stream the sun comes
     *               from set_two_stream_sun()).
     * @param canopy With ``radiation = two_stream``, this step's two-stream view of
     *               the level (see compute_longwave()), or an empty one before the
     *               first sweep, when every face gets no shortwave.
     */
    void compute_shortwave (amrex::Real time, const TwoStreamCanopyView& canopy = TwoStreamCanopyView{});

    /**
     * The highest interface of the level's vertical grid that any of its faces
     * samples with radiation = two_stream (a roof the bottom of its fluid cell, a
     * ceiling and a wall the top of its cell too), counted from the domain bottom, over all
     * ranks; -1 when the level has no faces. The two-stream beam the faces read
     * holds the interfaces up to there.
     */
    int top_sample_interface () const { return m_top_sample; }

    /** Whether the level has faces on any rank. */
    bool has_faces () const { return m_top_sample >= 0; }

    /**
     * Stop unless the fields of ``canopy`` (this level's two-stream view) lie on
     * the columns of ``cons``, with the same distribution: the two-stream kernels
     * read the column of a face from the fabs with the local index of the face's
     * own fab.
     */
    void check_canopy_layout (const amrex::MultiFab& cons, const TwoStreamCanopyView& canopy) const;

    // The two-stream branches of compute_shortwave() (once the sun is placed) and
    // compute_longwave(). Public only because a CUDA extended lambda cannot sit in
    // a private member function; call the two functions above instead.
    void compute_shortwave_two_stream (const TwoStreamCanopyView& canopy);
    void compute_longwave_two_stream (const amrex::MultiFab& cons, const TwoStreamCanopyView& canopy);

    /**
     * Scatter one per-face array into a cell-centred field as the mean over
     * the faces touching each cell (zero where there are none).
     */
    void scatter_field (const amrex::Gpu::DeviceVector<amrex::Real>& v, amrex::MultiFab& out) const;

    /**
     * Write every face of this rank to ``<prefix>.rank<N>.csv`` with its
     * geometry, view fractions, shadow flag, shortwave terms and skin
     * temperature, for the regression tests and for users who want the raw list.
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
     * Write the face state (skin temperature, sensible flux and slab
     * temperatures) to ``file`` for the checkpoint: one record per face of
     * the level, in one global order, written by the I/O rank.
     *
     * The file holds one text line, ``IBSEBFaces 1 <nface> <n_layers> <nx>
     * <ny> <nz>`` (the format version, the level's face count, the slab
     * layers per face and the domain size in cells), then in binary a marker
     * word (``0x0102030405060708`` as a 64-bit integer, so a file from a
     * machine with another byte order is refused), the ``nface`` keys as
     * 64-bit integers in ascending order, then the
     * ``nface x (2 + n_layers)`` values as doubles: a face's skin
     * temperature, its sensible flux (read as the previous step's by the
     * convective velocity scale) and its slab layers from the skin inward.
     * The key of a face is ``((i ny + j) nz + k) 6 + 2 dir + (side > 0)``
     * with ``(i, j, k)`` counted from the domain's low corner, so it names
     * the face whatever the rank count. The values are doubles whatever
     * ``amrex::Real`` is, so a single-precision run reads a double-precision
     * checkpoint and the other way round. The file costs ``8 (3 + n_layers)``
     * bytes per face beyond its header line and marker, which is the data itself: the cell field it replaces
     * covered the face-owning cells with whole boxes and was about eleven
     * times larger on the FaceStorage skyscraper (erf-model/ERF#4019).
     *
     * Every rank copies its faces to the host and sorts them by key; the I/O
     * rank gathers the records, orders them and writes the file. Collective:
     * every rank must call it.
     */
    void write_state (const std::string& file) const;

    /**
     * Read the face state back from a file written by write_state(), after
     * build() has recreated the list from the blanking.
     *
     * The I/O rank reads the table; every rank sends it the keys of its own
     * faces; the I/O rank finds each one by binary search and sends the
     * records back. Only the I/O rank ever holds the whole table, and a
     * restart on any number of ranks finds every face. Stops when the file
     * holds another number of slab layers than the deck's
     * ``erf.ibseb.n_slab_layers``, another domain size or another number of
     * faces, or when a face of this list has no record: a checkpoint of other
     * buildings. Collective: every rank must call it.
     */
    void read_state (const std::string& file);

    /**
     * Print the ``[IBSEB]`` summary line and, when asked, append one row per
     * building to ``erf.ibseb.csv_file``.
     *
     * The summary carries the global face counts per direction, the number of
     * buildings, the total face area and the range of the skin temperature;
     * the CSV rows carry per building the face count, the area, the
     * area-weighted mean fluxes and skin temperature with its range, the
     * largest balance residual, and the sun the step used, from its start
     * (zenith, azimuth, direct-normal and diffuse irradiance), on a row
     * stamped with the time at the end of the step. All numbers are reduced over the
     * ranks; the file is written by the I/O rank only, with a header only
     * when the file does not exist yet, so a restarted run appends. Every
     * level appends to the same file, its rows marked by the level column.
     * With ``erf.ibseb.dump_faces_file`` the rows come with a face dump,
     * dump_faces() to ``<prefix>[.levN][.stepNNNNNN].rank<N>.csv``: a refined
     * level N adds ``.levN`` so it does not overwrite the level below, and
     * ``erf.ibseb.dump_faces_tag_step`` the step.
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
    // State, advanced by the balance and checkpointed:
    amrex::Gpu::DeviceVector<amrex::Real> d_T_skin;        ///< Skin temperature [K]
    amrex::Gpu::DeviceVector<amrex::Real> d_T_slab;        ///< Slab layers [K], face-major: [face * n_layers + layer]
    // Static radiative geometry:
    amrex::Gpu::DeviceVector<amrex::Real> d_f_sky, d_f_ground, d_f_bldg;  ///< View fractions, sum to 1
    // Fluxes of the current step [W/m2]: SW_abs and LW_net positive into the
    // face, H, LE and G positive out of it:
    amrex::Gpu::DeviceVector<amrex::Real> d_SW_abs, d_LW_net, d_H, d_LE, d_G;
    amrex::Gpu::DeviceVector<amrex::Real> d_Q_ext;         ///< External incident flux [W/m2], e.g. fire radiation
    // Shortwave detail, all [W/m2] except the flag:
    amrex::Gpu::DeviceVector<amrex::Real> d_shadow;        ///< 1 when the ray to the sun hits a building, else 0
    amrex::Gpu::DeviceVector<amrex::Real> d_SW_direct_in;  ///< Direct beam incident on the face
    amrex::Gpu::DeviceVector<amrex::Real> d_SW_diffuse_in; ///< Sky and ground-reflected diffuse incident
    // Longwave detail [W/m2] and the air temperature read [K]:
    amrex::Gpu::DeviceVector<amrex::Real> d_LW_down_in;   ///< Incoming longwave on the face (sky, ground and walls)
    amrex::Gpu::DeviceVector<amrex::Real> d_T_air;        ///< Air temperature of the fluid cell
    // Sensible heat detail:
    amrex::Gpu::DeviceVector<amrex::Real> d_theta_air;    ///< Potential temperature of the fluid cell [K]
    amrex::Gpu::DeviceVector<amrex::Real> d_rho;          ///< Density of the fluid cell [kg/m3]
    amrex::Gpu::DeviceVector<amrex::Real> d_U_tan;        ///< Tangential wind speed at the fluid cell [m/s]
    amrex::Gpu::DeviceVector<amrex::Real> d_ustar;        ///< Friction velocity of the wall function [m/s]
    amrex::Gpu::DeviceVector<amrex::Real> d_H_coeff;      ///< H = H_coeff (theta_skin - theta_air) [W/m2/K]
    // Wall function beyond neutral:
    amrex::Gpu::DeviceVector<amrex::Real> d_hbld;         ///< Height of the face's building column [m]
    amrex::Gpu::DeviceVector<amrex::Real> d_w_star;       ///< Convective velocity scale of the step [m/s]
    amrex::Gpu::DeviceVector<amrex::Real> d_olen;         ///< Obukhov length used [m] (1e30 when neutral)
    amrex::Gpu::DeviceVector<amrex::Real> d_z_i;          ///< Depth in w* [m]: mixed layer above a roof, building height for a wall
    // Prognostic balance detail:
    amrex::Gpu::DeviceVector<amrex::Real> d_LW_ext;       ///< Incoming longwave from sky and ground only [W/m2]
    amrex::Gpu::DeviceVector<amrex::Real> d_resid;        ///< |balance| after the solve [W/m2]
    amrex::Gpu::DeviceVector<int>         d_niter;        ///< Newton iterations of the last solve
    // Column tops of the level's own columns: a k-collapsed map on the level's
    // boxes (the grids with their k range collapsed to one cell) with one
    // ghost cell, filled across boxes and periodic seams. Each entry is the
    // domain-relative k of the highest solid cell of its column, -1 where the
    // column is fluid. Kept after build() so the ray map can be rebuilt with
    // the buildings of the levels below (add_outside_occluders()).
    amrex::iMultiFab m_col_valid;
    // The ray map: one box per rank (the ranks with boxes), this rank's boxes
    // grown by m_halo, clipped to the domain where it is not periodic and
    // reduced to one domain width (with wrapping) where it is. Holds the
    // level's own tops, the periodic images and, after add_outside_occluders(),
    // the tops of the levels below outside this level.
    amrex::iMultiFab m_ray_patch;
    amrex::Box       m_patch_box;              ///< This rank's patch (empty without boxes)
    amrex::IntVect   m_halo{AMREX_D_DECL(0, 0, 0)}; ///< Halo of the patch in cells
    bool m_wrap_x = false, m_wrap_y = false;   ///< The patch spans the periodic domain in x / y
    bool m_map_spans_domain = false;           ///< Every rank's patch spans the domain
    int  m_max_steps = 0;                      ///< Columns a walk is followed for: the halo's reach in x plus y, before the domain cap

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

    /** Add the wall-clock time of one per-step update of this rank, for the cost line of report(). */
    void add_cost (double seconds) { m_cost_s += seconds; ++m_cost_n; }
    /** Record the wall-clock time of the initialisation (build and view fractions). */
    void set_init_cost (double seconds) { m_init_cost_s = seconds; }

    /**
     * First face of each local fab in MFIter order; size number of fabs + 1.
     * Faces ``fab_start()[n]`` to ``fab_start()[n+1]-1`` belong to the fab with
     * local index ``n``, which is what the scatter kernels rely on.
     */
    const std::vector<int>& fab_start () const { return m_fab_start; }

    /**
     * Build the ray map (m_ray_patch) from m_col_valid and the levels below.
     * The halo is sized for the shallowest ray the map must follow to its
     * end: the shallowest ray of the hemisphere sampling
     * (ibseb::hemisphere_min_sz()) and a sun ``erf.ibseb.sun_min_elev_deg``
     * above the horizon, over the tallest column of this level and the levels
     * below, ``height / tan(elevation)`` cells in x and in y
     * (ibseb::ray_halo_cells()). A rising ray stops above the tallest column
     * and a falling one at the ground, so within that reach the walk is exact,
     * and beyond it only a sun lower than the floor is concerned, whose
     * direct beam the clear-sky attenuation has all but removed; a building
     * farther than the halo from a face then casts no shadow on it (the far
     * field is open ground). The patch is this rank's boxes grown by the
     * halo, clipped to the domain where it is not periodic and, where it is,
     * reduced to one domain width when the grown box is at least that wide,
     * in which case the walk wraps round and is exact for any reach.
     * @param below Face sets of the levels below (empty on level 0), each with its build() done.
     *
     * Public because it launches a GPU lambda, which nvcc rejects in a
     * private member function.
     */
    void build_ray_map (const std::vector<const IBFaceSet*>& below);

private:
    IBSEBParams m_params;          ///< Copy of the inputs
    int m_lev   = 0;               ///< AMR level, for reports
    int m_qv_comp = -1;            ///< State component of rho qv, -1 without a moisture model
    int m_nface = 0;               ///< Faces owned by this rank
    int m_top_sample = -1;         ///< Highest interface a face samples, over all ranks (-1: no faces)
    int m_nbld  = 0;               ///< Buildings of the level (global)
    std::vector<int> m_fab_start;  ///< See fab_start()

    // Static totals over all ranks, computed once in build(); building index 0 unused.
    std::vector<amrex::Long> m_bld_nface;      ///< Faces per building
    std::vector<amrex::Real> m_bld_area;       ///< Face area per building [m2]
    amrex::Long m_nface_dir[3] = {0, 0, 0};    ///< Faces per direction
    amrex::Real m_area_total = 0.0;            ///< Total face area [m2]
    // Footprint of each building (global): columns, and bounding box in cell indices.
    std::vector<amrex::Long> m_bld_ncol;
    std::vector<int>  m_bld_ilo, m_bld_ihi, m_bld_jlo, m_bld_jhi;
    // Level-0 numbering (map_buildings_to_level0()): the level-0 id of each
    // building of this level (identity on level 0), and the building id of
    // every column on the level's own boxes (0 where fluid), kept until the
    // level above has been mapped.
    std::vector<int>  m_to_lev0;
    amrex::iMultiFab  m_lab;

    // Geometry kept for the ray cast and the reports.
    amrex::Geometry m_geom;                    ///< Geometry of the level
    int  m_nx = 0, m_ny = 0, m_nz = 0;         ///< Cells of the domain in x, y and z
    amrex::IntVect m_dom_lo{AMREX_D_DECL(0, 0, 0)}; ///< Low corner of the level's domain
    amrex::Real m_x_lo = 0.0, m_y_lo = 0.0;    ///< Domain origin [m]
    amrex::Real m_dx[3] = {0.0, 0.0, 0.0};     ///< Cell sizes [m]
    bool m_per_x = false, m_per_y = false;     ///< Periodicity in x and y
    amrex::Real m_col_top_max = 0.0;           ///< Tallest column the ray map holds [m] (the levels below included after add_outside_occluders())
    amrex::Real m_z_ground = 0.0;              ///< Ground height [m] (flat, the domain bottom)
    SunState m_sun;
    /// The two-stream sun of the next compute_shortwave() (sun_mode = two_stream).
    struct TwoStreamSun { bool set = false; amrex::Real decl = 0.0, hour_angle = 0.0, lat_deg = 0.0, S0 = 0.0; };
    TwoStreamSun m_ts_sun;
    double m_cost_s = 0.0;        ///< Wall-clock seconds spent in the per-step update on this rank
    amrex::Long m_cost_n = 0;     ///< Number of per-step updates timed
    double m_init_cost_s = 0.0;   ///< Wall-clock seconds of build() and compute_view_fractions()

    /** Warn once when a fixed sun stands below ``erf.ibseb.sun_min_elev_deg`` and the map does not span the domain. */
    void warn_fixed_sun_below_floor () const;
};

#endif

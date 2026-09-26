#ifndef ERF_SURFACEMODEL_H
#define ERF_SURFACEMODEL_H

#include "AMReX_Geometry.H"
#include "AMReX_ParmParse.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_IntVect.H"
#include "AMReX_MultiFab.H"
#include "AMReX_FillPatchUtil.H"
#include "AMReX_GpuLaunch.H"

#include <ERF_DataStruct.H>

enum SurfaceModelType { LAND = 0, URBAN = 1};

/**
 * @brief Interface between ERF and land-surface and urban models.
 *
 * SurfaceModel collects data and fluxes from the LSM and urban models,
 * performs the land-urban weight averaging and exposes the resulting fields
 * for use by ERF boundary-condition and physics components.
 *
 * Models may register data and mappings independently, without requiring a
 * fixed ordering or one-to-one correspondence between their fields.
 *
 * Model-specific field names can also be registered under a common name for
 * consumers such as RRTMGP and MOST.
 *
 * This class supports its own checkpoint state, ensuring the mappings are
 * restored correctly on restart. This class also provides a surface plotfile
 * for weight-averaged model data for the fields registered.
 */
class SurfaceModel
{
public:

    /**
     * @brief Constructs a surface-model interface for the supplied AMR levels.
     *
     * @param nlevs Number of AMR levels initially available.
     * @param ba Cell-centered box arrays for each level.
     * @param geom Geometry for each level.
     * @param dm Distribution mappings for each level.
     * @param solverChoice ERF solver configuration.
     * @param lmask_lev Land-mask data for each level.
     */
    explicit SurfaceModel(int nlevs,
                          const amrex::Vector<amrex::BoxArray> &ba,
                          const amrex::Vector<amrex::Geometry> &geom,
                          const amrex::Vector<amrex::DistributionMapping> &dm,
                          const SolverChoice &/*solverChoice*/,
                          amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>>& /*lmask_lev*/)
        : m_nlevs(nlevs), m_ba(ba), m_geom(geom), m_dmap(dm) {

        u_star.resize(m_nlevs);
        t_star.resize(m_nlevs);
        q_star.resize(m_nlevs);
        t_surf.resize(m_nlevs);

        wavg.resize(m_nlevs);

        m_ba2d.resize(m_nlevs);
        m_geom2d.resize(m_nlevs);
        m_lmask.resize(m_nlevs);


        urban_frac_lev.resize(m_nlevs);
        lsm_data_lev.resize(m_nlevs);
        urban_data_lev.resize(m_nlevs);
        weighted_lsm_data_lev.resize(m_nlevs);
        weighted_urban_data_lev.resize(m_nlevs);

        rad_fields.resize(m_nlevs);
        rad_output_fields.resize(m_nlevs);

        // set initial fields for each model to none
        lsm_fields = amrex::Vector<int>(5, -1);
        urban_fields = amrex::Vector<int>(5, -1);
    }

    /**
     * @brief Initializes storage and boundary data for an AMR level.
     *
     * Existing coarser-level data are interpolated when a finer level is
     * initialized.
     *
     * @param lev AMR level to initialize.
     * @param ba Cell-centered box array for the level.
     * @param geom Geometry for the level.
     * @param dm Distribution mapping for the level.
     * @param lmask_lev Land-mask data used by the level.
     * @param domain_bcs_type Boundary conditions for interpolation.
     * @param refRatio Refinement ratios between adjacent levels.
     */
    void initialize_for_level(int lev, const amrex::BoxArray &ba, const amrex::Geometry &geom, const amrex::DistributionMapping &dm, const amrex::Vector<std::unique_ptr<amrex::iMultiFab>>& lmask_lev, const amrex::Vector <amrex::BCRec> &domain_bcs_type, const amrex::Vector<amrex::IntVect> &refRatio)
    {
        if (lev >= m_nlevs) {
            m_nlevs = lev+1;

            u_star.resize(m_nlevs);
            t_star.resize(m_nlevs);
            q_star.resize(m_nlevs);
            t_surf.resize(m_nlevs);

            wavg.resize(m_nlevs);

            m_dmap.resize(m_nlevs);
            m_ba.resize(m_nlevs);
            m_ba2d.resize(m_nlevs);
            m_geom.resize(m_nlevs);
            m_geom2d.resize(m_nlevs);
            m_lmask.resize(m_nlevs);

            urban_frac_lev.resize(m_nlevs);
            lsm_data_lev.resize(m_nlevs);
            urban_data_lev.resize(m_nlevs);
            weighted_lsm_data_lev.resize(m_nlevs);
            weighted_urban_data_lev.resize(m_nlevs);

            rad_fields.resize(m_nlevs);
            rad_output_fields.resize(m_nlevs);
        }

        m_ba[lev] = ba;
        m_geom[lev] = geom;
        m_dmap[lev] = dm;

        // Create a 2D ba, dm, & ghost cells
        amrex::BoxList bl2d = ba.boxList();
        for (auto& b : bl2d) {
            b.setRange(2, 0);
        }
        m_ba2d[lev] = amrex::BoxArray(std::move(bl2d));

        amrex::RealBox dom2d = geom.ProbDomain();
        dom2d.setHi(2, geom.CellSize(2));
        m_geom2d[lev].define(makeSlab(m_geom[lev].Domain(), 2, 0), dom2d, geom.Coord(), geom.isPeriodic());

        m_lmask[lev] = lmask_lev[0].get();

        weighted_lsm_data_lev[lev].clear();
        weighted_urban_data_lev[lev].clear();

        amrex::IntVect ng(1,1,0);
        // Create output arrays for surface boundary condition
        u_star[lev] = std::make_unique<amrex::MultiFab>(m_ba2d[lev], m_dmap[lev], 2, ng);
        t_star[lev] = std::make_unique<amrex::MultiFab>(m_ba2d[lev], m_dmap[lev], 1, ng);
        q_star[lev] = std::make_unique<amrex::MultiFab>(m_ba2d[lev], m_dmap[lev], 1, ng);
        t_surf[lev] = std::make_unique<amrex::MultiFab>(m_ba2d[lev], m_dmap[lev], 1, ng);

        u_star[lev]->setVal(0.0);
        t_star[lev]->setVal(0.0);
        q_star[lev]->setVal(0.0);
        t_surf[lev]->setVal(300.0);

        // Initialize storage for weight average factors
        wavg[lev] = std::make_unique<amrex::MultiFab>(m_ba2d[lev], m_dmap[lev], 2, ng);
        //  - default to land=1, urban=0 (i.e, 100% land, no urban fraction)
        // the weight average factors for this new level will get calculated on the next call to weight average
        wavg[lev]->setVal(1.0, SurfaceModelType::LAND, 1, 0);
        wavg[lev]->setVal(0.0, SurfaceModelType::URBAN, 1, 0);

        // Recreate registered fields on the current level layout. This is also
        // needed when an existing AMR level is remade with a new grid layout.
        for (int field = 0; field < static_cast<int>(fields.size()); ++field) {
            if (static_cast<int>(fields[field].size()) < m_nlevs) {
                fields[field].resize(m_nlevs);
            }
            fields[field][lev] = std::make_unique<amrex::MultiFab>(
                m_ba2d[lev], m_dmap[lev], 1, amrex::IntVect(1,1,0));
            fields[field][lev]->setVal(0.0);
        }
        for (auto& entry : radiation_input_map) {
            if (static_cast<int>(entry.second.weighted.size()) < m_nlevs) {
                entry.second.weighted.resize(m_nlevs);
            }
            entry.second.weighted[lev].reset();
        }

        // Cached radiation fields point into the registered storage above.
        // Rebuild the cache when this level's fields have been replaced.
        rad_fields[lev].clear();
        rad_output_fields[lev].clear();

        // Interpolate the new finer level using coarse data
        if (lev > 0) {
            amrex::Real time_for_fp = zero; // This is not actually used
            amrex::Vector<amrex::Real> ftime    = {time_for_fp, time_for_fp};
            amrex::Vector<amrex::Real> ctime    = {time_for_fp, time_for_fp};
            amrex::IntVect ng_od(0, 0, 0); // ng outside domain
            amrex::InterpFromCoarseLevel(*u_star[lev], u_star[lev-1]->nGrowVect(), ng_od,
                                         *u_star[lev-1], 0, 0, 2,
                                         m_geom[lev-1], m_geom[lev],
                                         refRatio[lev-1], &amrex::cell_cons_interp,
                                         domain_bcs_type, BCVars::cons_bc);
            amrex::InterpFromCoarseLevel(*t_star[lev], t_star[lev-1]->nGrowVect(), ng_od,
                                         *t_star[lev-1], 0, 0, 1,
                                         m_geom[lev-1], m_geom[lev],
                                         refRatio[lev-1], &amrex::cell_cons_interp,
                                         domain_bcs_type, BCVars::cons_bc);
            amrex::InterpFromCoarseLevel(*q_star[lev], q_star[lev-1]->nGrowVect(), ng_od,
                                         *q_star[lev-1], 0, 0, 1,
                                         m_geom[lev-1], m_geom[lev],
                                         refRatio[lev-1], &amrex::cell_cons_interp,
                                         domain_bcs_type, BCVars::cons_bc);
            amrex::InterpFromCoarseLevel(*t_surf[lev], t_surf[lev-1]->nGrowVect(), ng_od,
                                         *t_surf[lev-1], 0, 0, 1,
                                         m_geom[lev-1], m_geom[lev],
                                         refRatio[lev-1], &amrex::cell_cons_interp,
                                         domain_bcs_type, BCVars::cons_bc);

            amrex::Interpolater* mapper = &amrex::cell_cons_interp;
            amrex::Vector<amrex::MultiFab*> fmf = {u_star[lev  ].get(), u_star[lev  ].get()};
            amrex::Vector<amrex::MultiFab*> cmf = {u_star[lev-1].get(), u_star[lev-1].get()};
            amrex::FillPatchTwoLevels(*u_star[lev].get(), u_star[lev]->nGrowVect(), ng_od,
                                      time_for_fp, cmf, ctime, fmf, ftime,
                                      0, 0, 2, m_geom[lev-1], m_geom[lev],
                                      refRatio[lev-1], mapper, domain_bcs_type,
                                      BCVars::cons_bc);
            fmf = {t_star[lev  ].get(), t_star[lev  ].get()};
            cmf = {t_star[lev-1].get(), t_star[lev-1].get()};
            amrex::FillPatchTwoLevels(*t_star[lev].get(), t_star[lev]->nGrowVect(), ng_od,
                                      time_for_fp, cmf, ctime, fmf, ftime,
                                      0, 0, 1, m_geom[lev-1], m_geom[lev],
                                      refRatio[lev-1], mapper, domain_bcs_type,
                                      BCVars::cons_bc);
            fmf = {q_star[lev  ].get(), q_star[lev  ].get()};
            cmf = {q_star[lev-1].get(), q_star[lev-1].get()};
            amrex::FillPatchTwoLevels(*q_star[lev].get(), q_star[lev]->nGrowVect(), ng_od,
                                      time_for_fp, cmf, ctime, fmf, ftime,
                                      0, 0, 1, m_geom[lev-1], m_geom[lev],
                                      refRatio[lev-1], mapper, domain_bcs_type,
                                      BCVars::cons_bc);
            fmf = {t_surf[lev  ].get(), t_surf[lev  ].get()};
            cmf = {t_surf[lev-1].get(), t_surf[lev-1].get()};
            amrex::FillPatchTwoLevels(*t_surf[lev].get(), t_surf[lev]->nGrowVect(), ng_od,
                                      time_for_fp, cmf, ctime, fmf, ftime,
                                      0, 0, 1, m_geom[lev-1], m_geom[lev],
                                      refRatio[lev-1], mapper, domain_bcs_type,
                                      BCVars::cons_bc);

            for (int field = 0; field < static_cast<int>(fields.size()); ++field) {
                amrex::IntVect ngv = fields[field][lev-1]->nGrowVect(); ngv[2] = 0;
                amrex::InterpFromCoarseLevel(*fields[field][lev], ngv, ng_od,
                                             *fields[field][lev-1], 0, 0, 1,
                                             m_geom[lev-1], m_geom[lev],
                                             refRatio[lev-1], &amrex::cell_cons_interp,
                                             domain_bcs_type, BCVars::cons_bc);
                fmf = {fields[field][lev  ].get(), fields[field][lev  ].get()};
                cmf = {fields[field][lev-1].get(), fields[field][lev-1].get()};
                amrex::FillPatchTwoLevels(*fields[field][lev].get(), ngv, ng_od,
                                          time_for_fp, cmf, ctime, fmf, ftime,
                                          0, 0, 1, m_geom[lev-1], m_geom[lev],
                                          refRatio[lev-1], mapper, domain_bcs_type,
                                          BCVars::cons_bc);
            }
        }
    }

    /**
     * @brief Registers model variables used in weighted averaging.
     *
     * @param lev AMR level associated with the data.
     * @param model_data Pointers to the model variables.
     * @param data_names Names corresponding to @p model_data.
     * @param type Whether the data belong to the land or urban model.
     */
    void set_model_data(const int lev, const amrex::Vector<amrex::MultiFab*> model_data, const amrex::Vector<std::string>& data_names, SurfaceModelType type)
    {
        if (type == SurfaceModelType::LAND) {
            lsm_data_lev[lev].resize(model_data.size());
            weighted_lsm_data_lev[lev].clear();
            weighted_lsm_data_lev[lev].resize(model_data.size());
            if (lev == 0) {
                // only set names on first level, since they are the same for all levels
                m_lsm_names.resize(model_data.size());
                AMREX_ALWAYS_ASSERT(data_names.size() == model_data.size());
                for (int i = 0; i < static_cast<int>(data_names.size()); ++i) {
                    m_lsm_names[i] = data_names[i];
                }
            }

            for (int i = 0; i < static_cast<int>(lsm_data_lev[lev].size()); ++i) {
                AMREX_ALWAYS_ASSERT(model_data[i]);
                lsm_data_lev[lev][i] = model_data[i];
            }
            m_use_land = true;
            for (auto& fields_at_level : rad_fields) { fields_at_level.clear(); }
            for (auto& fields_at_level : rad_output_fields) { fields_at_level.clear(); }
            for (auto& entry : radiation_input_map) {
                if (lev < static_cast<int>(entry.second.weighted.size())) {
                    entry.second.weighted[lev].reset();
                }
            }
        } else if (type == SurfaceModelType::URBAN) {
            urban_data_lev[lev].resize(model_data.size());
            weighted_urban_data_lev[lev].clear();
            weighted_urban_data_lev[lev].resize(model_data.size());
            if (lev == 0) {
                // only set names on first level, since they are the same for all levels
                m_urban_names.resize(model_data.size());
                AMREX_ALWAYS_ASSERT(data_names.size() == model_data.size());
                for (int i = 0; i < static_cast<int>(data_names.size()); ++i) {
                    m_urban_names[i] = data_names[i];
                }
            }
            for (int i = 0; i < static_cast<int>(urban_data_lev[lev].size()); ++i) {
                AMREX_ALWAYS_ASSERT(model_data[i]);
                urban_data_lev[lev][i] = model_data[i];
            }
            m_use_urban = true;
            for (auto& fields_at_level : rad_fields) { fields_at_level.clear(); }
            for (auto& fields_at_level : rad_output_fields) { fields_at_level.clear(); }
            for (auto& entry : radiation_input_map) {
                if (lev < static_cast<int>(entry.second.weighted.size())) {
                    entry.second.weighted[lev].reset();
                }
            }
        }
    }

    /**
     * @brief Registers model fluxes as additional weighted-average inputs.
     *
     * @param lev AMR level associated with the fluxes.
     * @param model_fluxes Pointers to the model flux fields.
     * @param flux_names Names corresponding to @p model_fluxes.
     * @param type Whether the fluxes belong to the land or urban model.
     */
    void set_model_fluxes(const int lev, const amrex::Vector<amrex::MultiFab*> model_fluxes,
                          const amrex::Vector<std::string>& flux_names, SurfaceModelType type)
    {
        amrex::Vector<amrex::MultiFab*>* model_data = nullptr;
        amrex::Vector<std::string>* model_names = nullptr;

        if (type == SurfaceModelType::LAND) {
            model_data = &lsm_data_lev[lev];
            model_names = &m_lsm_names;
        } else if (type == SurfaceModelType::URBAN) {
            model_data = &urban_data_lev[lev];
            model_names = &m_urban_names;
        }

        AMREX_ALWAYS_ASSERT(model_data != nullptr);
        AMREX_ALWAYS_ASSERT(model_names != nullptr);
        AMREX_ALWAYS_ASSERT(flux_names.size() == model_fluxes.size());

        const int data_size = static_cast<int>(model_data->size());
        model_data->resize(data_size + static_cast<int>(model_fluxes.size()));
        if (lev == 0) {
            model_names->resize(data_size + static_cast<int>(model_fluxes.size()));
        }

        for (int i = 0; i < static_cast<int>(model_fluxes.size()); ++i) {
            AMREX_ALWAYS_ASSERT(model_fluxes[i]);
            (*model_data)[data_size + i] = model_fluxes[i];
            if (lev == 0) {
                (*model_names)[data_size + i] = flux_names[i];
            }
        }
    }

    /**
     * @brief Selects the model fields that participate in weighted averaging.
     *
     * When using fluxes, the first five elements map to u-stress, v-stress,
     * heat flux, vapor flux, and surface temperature. Otherwise, the first
     * four elements map to ustar, tstar, qstar, and surface temperature. A
     * field index of -1 skips that field.
     *
     * @param type Whether the fields belong to the land or urban model.
     * @param field_indices Indices of the fields to average.
     * @param use_fluxes Whether the selected fields are fluxes.
     */
    void set_model_fields(SurfaceModelType type, const amrex::Vector<int> &field_indices, const bool use_fluxes = true)
    {
        const int nfields = static_cast<int>(field_indices.size());
        if (!use_fluxes) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nfields >= 4, "Must have at least four fields to average (ustar, tstar, qstar, tsurf)");
        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nfields >= 5, "Must have at least five fields to average (u-stress, v-stress, hfx, qfx, tsurf)");
        }

        // Latch the export mode off the first model that registers, then make sure every
        // subsequent model agrees: either all of them supply fluxes, or all of them supply
        // u*, t*, q* directly.  Taking the mode from the caller rather than assuming the
        // default is what makes the non-flux path reachable at all.
        if (!m_export_fluxes_set) {
            m_export_fluxes = use_fluxes;
            m_export_fluxes_set = true;
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(use_fluxes == m_export_fluxes,
            "All registered surface models must export either fluxes or u*/t*/q*, not a mix");

        const int nregistered = (type == SurfaceModelType::LAND)
            ? static_cast<int>(m_lsm_names.size())
            : static_cast<int>(m_urban_names.size());
        for (const int field : field_indices) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                field == -1 || (field >= 0 && field < nregistered),
                "Model field index must be -1 or within the registered provider data range");
        }

        if (type == SurfaceModelType::LAND) {
            lsm_fields = amrex::Vector<int>(field_indices);
            for (auto& fields_at_level : weighted_lsm_data_lev) {
                fields_at_level.clear();
            }
        } else if (type == SurfaceModelType::URBAN) {
            urban_fields = amrex::Vector<int>(field_indices);
            for (auto& fields_at_level : weighted_urban_data_lev) {
                fields_at_level.clear();
            }
        }
    }

    /**
     * @brief Applies land and urban weights into separate output fields.
     *
     * A null land or urban source is skipped, but both sources cannot be
     * null. The source fields are not modified and the fields may use
     * different underlying grids.
     *
     * @param lev AMR level whose weights are used.
     * @param lsm_data Land-model source data to weight.
     * @param lsm_weighted Destination for weighted land data.
     * @param urban_data Urban-model source data to weight.
     * @param urban_weighted Destination for weighted urban data.
     */
    void apply_weight_average(int lev, const amrex::MultiFab *lsm_data,
                              amrex::MultiFab *lsm_weighted,
                              const amrex::MultiFab *urban_data,
                              amrex::MultiFab *urban_weighted);

    /**
     * @brief Computes land and urban weights and applies them to registered fields.
     *
     * @param lev AMR level to update.
     * @param urban_frac Urban fraction used to compute the weights.
     */
    void calculate_weight_average(int lev, amrex::MultiFab* const urban_frac);

    /**
     * @brief Returns whether the interface exports fluxes rather than MOST variables.
     * @return True when weighted outputs are fluxes.
     */
    bool are_fluxes() { return m_export_fluxes; }

    /**
     * @brief Marks the registered surface fields as filled by an actual model advance.
     */
    void mark_fields_valid() { m_fields_are_valid = true; }

    /**
     * @brief Returns whether the surface models have advanced at least once.
     * @return False until the first advance, while the registered fields still hold their
     *         initial values rather than anything the models computed.
     */
    bool fields_are_valid() const { return m_fields_are_valid; }

    /**
     * @brief Returns the weighted horizontal momentum output for a level.
     * @param lev AMR level to query.
     * @return Pointer to the weighted u-star or momentum-flux field.
     */
    amrex::MultiFab* get_ustar(const int& lev) {
        return u_star[lev].get();
    }

    /**
     * @brief Returns the weighted thermal output for a level.
     * @param lev AMR level to query.
     * @return Pointer to the weighted t-star or heat-flux field.
     */
    amrex::MultiFab* get_tstar(const int& lev) {
        return t_star[lev].get();
    }

    /**
     * @brief Returns the weighted moisture output for a level.
     * @param lev AMR level to query.
     * @return Pointer to the weighted q-star or vapor-flux field.
     */
    amrex::MultiFab* get_qstar(const int& lev) {
        return q_star[lev].get();
    }

    /**
     * @brief Returns the weighted surface-temperature output for a level.
     * @param lev AMR level to query.
     * @return Pointer to the weighted surface-temperature field.
     */
    amrex::MultiFab* get_tsurf(const int& lev) {
        return t_surf[lev].get();
    }

    /**
     * @brief Returns a weighted copy of a selected model field.
     * @param lev AMR level to query.
     * @param type Provider type owning the field.
     * @param field_idx Provider field index.
     * @return Weighted field, or nullptr if it has not been selected.
     */
    const amrex::MultiFab* get_weighted_model_data(const int lev,
                                                    SurfaceModelType type,
                                                    const int field_idx) const {
        const auto& weighted = (type == SurfaceModelType::LAND)
            ? weighted_lsm_data_lev[lev] : weighted_urban_data_lev[lev];
        if (field_idx < 0 || field_idx >= static_cast<int>(weighted.size())) {
            return nullptr;
        }
        return weighted[field_idx].get();
    }

    /**
     * @brief Returns the land and urban weighting factors for a level.
     * @param lev AMR level to query.
     * @return Pointer to the weighting-factor field.
     */
    amrex::MultiFab* get_wavg_factors(const int &lev) {
        return wavg[lev].get();
    }

    /**
     * @brief Writes registered surface-model fields to plotfile output.
     * @param finest_lev Finest AMR level to write.
     * @param time Simulation time associated with the output.
     * @param plot_prefix Output plotfile prefix.
     * @param level_steps Number of steps represented at each level.
     * @param ref_ratio Refinement ratios between adjacent levels.
     */
    void write_output(const int &finest_lev, const amrex::Real &time, const std::string &plot_prefix, const amrex::Vector<int> &level_steps, const amrex::Vector<amrex::IntVect> &ref_ratio);

    /**
     * @brief Registers a field map using land and urban field indices.
     * @param name Common name used to retrieve the mapped field.
     * @param lsm_urb_map Pair of land and urban field indices.
     * @param fill_boundary Whether to fill mapped-field boundaries.
     */
    void register_field_map(std::string name, const std::pair<int, int> &lsm_urb_map, bool fill_boundary=false);

    /**
     * @brief Registers a field map using per-level land and urban pointers.
     * @param name Common name used to retrieve the mapped field.
     * @param lsm_lev_mf Land-model fields, indexed by AMR level.
     * @param urb_lev_mf Urban-model fields, indexed by AMR level.
     * @param fill_boundary Whether to fill mapped-field boundaries.
     */
    void register_field_map(std::string name, amrex::Vector<amrex::MultiFab*> &lsm_lev_mf, amrex::Vector<amrex::MultiFab*> &urb_lev_mf, bool fill_boundary=false);

    /**
     * @brief Updates a pointer-based field mapping for one level.
     * @param name Common name of the mapped field.
     * @param lev AMR level to update.
     * @param lsm_mf Land-model field pointer, or nullptr.
     * @param urban_mf Urban-model field pointer, or nullptr.
     */
    void set_field_map_pointers(const std::string& name, int lev,
                                amrex::MultiFab* lsm_mf, amrex::MultiFab* urban_mf);

    /**
     * @brief Registers a canonical radiation input mapping.
     * @param name Canonical radiation input name.
     * @param lsm_urb_map Pair of land and urban model field indices.
     */
    void register_radiation_input(const std::string& name,
                                  const std::pair<int, int>& lsm_urb_map);

    /**
     * @brief Registers canonical radiation input mappings.
     * @param input_map Canonical names mapped to land and urban indices.
     */
    void register_radiation_inputs(
        const std::unordered_map<std::string, std::pair<int, int>>& input_map);

    /**
     * @brief Registers a canonical radiation output mapping.
     * @param name Canonical radiation output name.
     * @param lsm_urb_map Pair of land and urban model field indices.
     */
    void register_radiation_output(const std::string& name,
                                   const std::pair<int, int>& lsm_urb_map);

    /**
     * @brief Registers canonical radiation output mappings.
     * @param output_map Canonical names mapped to land and urban indices.
     */
    void register_radiation_outputs(
        const std::unordered_map<std::string, std::pair<int, int>>& output_map);

    /**
     * @brief Retrieves a registered field by its common name.
     * @param name Registered field name.
     * @param lev AMR level to query.
     * @return Pointer to the registered field, or nullptr if it is unknown.
     */
    amrex::MultiFab* get_field(const std::string &name, int lev = 0) {
        if (auto field = fieldmap.find(name); field != fieldmap.end()) {
            AMREX_ASSERT(field->second.mf_ind != -1);
            return fields[field->second.mf_ind][lev].get();
        } else {
            return nullptr;
        }
    }

    /**
     * @brief Returns fields registered under the RRTMGP radiation names.
     *
     * The resulting map is cached after the first call for reuse on subsequent timesteps.
     * When an expected field is missing, it is marked nullptr for the RRTMGP implementation
     * to apply its fallback values.
     *
     * @param lev AMR level to query.
     * @return Radiation fields in the order expected by RRTMGP.
     */
    const amrex::Vector<const amrex::MultiFab*> get_radiation_fields(int lev = 0) {
        // check if we need to build the radiation list (for first time step,)
        if (rad_fields[lev].size() == 0) {
            build_radiation_fieldlist(lev);
        }
        return rad_fields[lev];
    }

    /**
     * @brief Returns canonical radiation output destinations for one level.
     * @param lev AMR level to query.
     * @return One destination per canonical radiation output.
     */
    const amrex::Vector<amrex::MultiFab*> get_radiation_output_fields(int lev = 0);

    /**
     * @brief Distributes updated canonical radiation outputs to other providers.
     * @param lev AMR level to update.
     */
    void distribute_radiation_outputs(int lev);

    /**
     * @brief Distributes one updated canonical radiation output to other providers.
     * @param lev AMR level to update.
     * @param output_index Index in the canonical radiation output list.
     */
    void distribute_radiation_output(int lev, int output_index);

    /**
     * @brief Advances an input stream to the next line.
     *
     * taken from ERF_Checkpoint, reuse?
     *
     * @param is Input stream to advance.
     */
    static void GotoNextLine (std::istream& is);

    /**
     * @brief Writes surface-model state to a checkpoint file.
     * @param checkpointname Checkpoint file name.
     */
    void WriteCheckpoint(const std::string &checkpointname);

    /**
     * @brief Restores surface-model state from a checkpoint file.
     * @param checkpointname Checkpoint file name.
     */
    void ReadCheckpoint(const std::string &checkpointname);

    /**
     * @brief Validates the layout contract for a canonical radiation output.
     * @param lev AMR level associated with the destination.
     * @param mf Radiation output destination to validate.
     */
    void validate_radiation_output_layout(int lev, const amrex::MultiFab* mf) const;

    /**
     * @brief Computes simple land and urban weights from the urban fraction.
     * @param lev AMR level to update.
     * @param urban_frac Urban fraction used to compute the weights.
     */
    void calculate_simple_average(int lev, amrex::MultiFab* const urban_frac);

    /**
     * @brief Builds the ordered radiation-field list for one AMR level.
     * @param lev AMR level to update.
     */
    void build_radiation_fieldlist(int lev) {
        amrex::Print() << "  --- building radiation field list" << std::endl;

        // loop over fields and check against names of radiation variables
        for (int i = 0; i < static_cast<int>(radnames.size()); ++i) {
            amrex::Print() << "    - Checking for radiation variable '" << radnames[i] << "'" << std::endl;

            amrex::MultiFab* mf = nullptr;
            auto rad_it = radiation_input_map.find(radnames[i]);
            if (rad_it != radiation_input_map.end()) {
                const int land_idx = rad_it->second.map.first;
                const int urban_idx = rad_it->second.map.second;
                const bool valid_land = m_use_land && land_idx >= 0 &&
                    lev < static_cast<int>(lsm_data_lev.size()) &&
                    land_idx < static_cast<int>(lsm_data_lev[lev].size()) &&
                    lsm_data_lev[lev][land_idx] != nullptr;
                const bool valid_urban = m_use_urban && urban_idx >= 0 &&
                    lev < static_cast<int>(urban_data_lev.size()) &&
                    urban_idx < static_cast<int>(urban_data_lev[lev].size()) &&
                    urban_data_lev[lev][urban_idx] != nullptr;
                if (valid_land && valid_urban) {
                    if (!rad_it->second.weighted[lev]) {
                        rad_it->second.weighted[lev] = std::make_unique<amrex::MultiFab>(
                            m_ba2d[lev], m_dmap[lev], 1, amrex::IntVect(1,1,0));
                    }
                    auto& weighted = *rad_it->second.weighted[lev];
                    const auto& land = *lsm_data_lev[lev][land_idx];
                    const auto& urban = *urban_data_lev[lev][urban_idx];
                    for (amrex::MFIter mfi(weighted); mfi.isValid(); ++mfi) {
                        const amrex::Box bx = mfi.tilebox();
                        const auto weights = wavg[lev]->const_array(mfi);
                        const auto land_arr = land.const_array(mfi);
                        const auto urban_arr = urban.const_array(mfi);
                        // LSM surface fields may be stored at k <= 0; use the top valid LSM
                        // plane, which matches k=0 when the LSM keeps those fields synchronized
                        // at the surface (such as SLM). Urban surface fields are stored at k=0.
                        const int land_k = land.box(mfi.index()).bigEnd(2);
                        const int urban_k = urban.box(mfi.index()).smallEnd(2);
                        auto output = weighted.array(mfi);
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                            output(i,j,k) = land_arr(i,j,land_k) *
                                weights(i,j,0,SurfaceModelType::LAND) +
                                urban_arr(i,j,urban_k) *
                                weights(i,j,0,SurfaceModelType::URBAN);
                        });
                    }
                    weighted.FillBoundary(m_geom[lev].periodicity());
                    mf = rad_it->second.weighted[lev].get();
                } else if (valid_land) {
                    mf = lsm_data_lev[lev][land_idx];
                } else if (valid_urban) {
                    mf = urban_data_lev[lev][urban_idx];
                }
            }
            if (mf) {
                rad_fields[lev].push_back(mf);
            } else {
                amrex::Print() << "     - WARNING: did not find radiation variable, setting to null to use default RRTMGP value!" << std::endl;
                rad_fields[lev].push_back(nullptr);
            }
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            static_cast<int>(rad_fields[lev].size()) == static_cast<int>(radnames.size()),
            "Did not find all radiation fields");
    }

    /**
     * @brief Applies the current land and urban weights to selected model fields.
     * @param lev AMR level to update.
     * @param urban_frac Urban fraction used to compute the weights.
     */
    void weight_average_fields(int lev, amrex::MultiFab* const /*urban_frac*/);

    /**
     * @brief Checks whether a model field participates in a registered mapping.
     * @param lev AMR level to query.
     * @param type Land or urban model owning the field.
     * @param field_idx Index of the model field.
     * @param mf Field pointer to test.
     * @return True when the field is part of a registered mapping.
     */
    bool is_field_mapped(int lev, SurfaceModelType type, int field_idx,
                         const amrex::MultiFab* mf) const;

    /**
     * @brief Writes a weighted copy of a model field.
     * @param lev AMR level whose weights are used.
     * @param source Model-owned source field.
     * @param weighted SurfaceModel-owned destination field.
     * @param type Provider type owning the source field.
     */
    void weight_model_field(int lev, const amrex::MultiFab* source,
                            amrex::MultiFab* weighted, SurfaceModelType type);

private:

    int m_nlevs;

    bool m_use_urban = false; // whether we have an urban model
    bool m_use_land = false; // whether we have a LSM

    bool m_export_fluxes = true;
    bool m_export_fluxes_set = false; // whether a model has registered and latched the mode above

    // Set once the surface models have actually integrated at least one step, so consumers can
    // tell registered-but-unfilled fields from real ones.
    bool m_fields_are_valid = false;

    // internal flag to indicate if the weight average has been calculated this time step
    bool m_weights_updated = false;

    amrex::Vector<amrex::BoxArray> m_ba;
    amrex::Vector<amrex::BoxArray> m_ba2d;
    amrex::Vector<amrex::Geometry> m_geom;
    amrex::Vector<amrex::DistributionMapping> m_dmap;

    amrex::Vector<amrex::Geometry> m_geom2d;
    amrex::Vector<amrex::iMultiFab*> m_lmask;

    // Pointers to LSM and Urban model data
    amrex::Vector<amrex::MultiFab*>  urban_frac_lev; // Urban fraction for each level [0-1.0]
    amrex::Vector<amrex::Vector<amrex::MultiFab*>>  lsm_data_lev;   // LSM data for each level
    amrex::Vector<amrex::Vector<amrex::MultiFab*>>  urban_data_lev; // Urban data for each level
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> weighted_lsm_data_lev;
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> weighted_urban_data_lev;

    // Output weighted averages passed to MOST
    //  These can either be the u*, t*, q* for MOST, or fluxes computed directly from the model
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> u_star; // up to 2 components: X-stress, Y-stress
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> t_star;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> q_star;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> t_surf;

    // Current weight average factors for each level (2 components: land and urban)
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> wavg;

    // Model variables to weight average
    // The flux interface uses u-stress, v-stress, heat flux, vapor flux, tsurf.
    // TODO: make this more flexible for different types of models
    amrex::Vector<int> lsm_fields;
    amrex::Vector<int> urban_fields;

    // These are the names to be consumed in RRTMGP, but the mapped LSM/Urban names can be different
    inline static const std::vector<std::string> radnames = {"tskin", "emiss", "albedo_vis", "albedo_nir", "albedo_vis_diff", "albedo_nir_diff"};

    /**
     * @brief Describes a registered land/urban field mapping.
     *
     * A Field associates the common name used by ERF consumers with either
     * model-specific field indices or per-level field pointers. The mapped
     * fields are stored in the SurfaceModel-owned field array and the source
     * pointers remain owned by their respective surface models.
     */
    struct Field {
        /// Land and urban model field indices, respectively.
        std::pair<int, int> map;
        /// Per-level land-model field pointers.
        amrex::Vector<amrex::MultiFab*> lsm_ptr;
        /// Per-level urban-model field pointers.
        amrex::Vector<amrex::MultiFab*> urb_ptr;
        /// Index of the mapped field in the SurfaceModel field array.
        int mf_ind;
        /// Whether boundaries are filled after mapping the field.
        bool fill_bound;
    };
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> fields; // NOTE: this is stored in inverse order: (field, nlev)
    std::unordered_map<std::string, Field> fieldmap;

    amrex::Vector<std::string> m_lsm_names;
    amrex::Vector<std::string> m_urban_names;

    amrex::Vector<amrex::Vector<const amrex::MultiFab*>> rad_fields; // pointers to radiation MFs at each level
    amrex::Vector<amrex::Vector<amrex::MultiFab*>> rad_output_fields;

    inline static const std::vector<std::string> rad_output_names = {
        "cos_zenith_angle", "sw_flux_dn", "sw_flux_dn_dir_vis",
        "sw_flux_dn_dir_nir", "sw_flux_dn_dif_vis", "sw_flux_dn_dif_nir",
        "lw_flux_dn"};

    struct RadiationField {
        std::pair<int, int> map;
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> weighted;
    };
    std::unordered_map<std::string, RadiationField> radiation_input_map;
    std::unordered_map<std::string, RadiationField> radiation_output_map;

    // plotfile names for outputs
    inline static const std::vector<std::string> field_names = {
        "uflux",
        "vflux",
        "hfx",
        "qfx",
        "tskin",
    };
};

#endif

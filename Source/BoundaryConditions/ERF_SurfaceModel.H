#ifndef ERF_SURFACEMODEL_H
#define ERF_SURFACEMODEL_H

#include "AMReX_Geometry.H"
#include "AMReX_ParmParse.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_IntVect.H"
#include "AMReX_MultiFab.H"
#include "AMReX_FillPatchUtil.H"

#include <ERF_DataStruct.H>

enum SurfaceModelType { LAND = 0, URBAN = 1};


class SurfaceModel
{
public:

    explicit SurfaceModel(int nlevs,
                          const amrex::Vector<amrex::BoxArray> &ba,
                          const amrex::Vector<amrex::Geometry> &geom,
                          const amrex::Vector<amrex::DistributionMapping> &dm,
                          const SolverChoice &solverChoice,
                          amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>>& lmask_lev)
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

        rad_fields.resize(m_nlevs);

        // set initial fields for each model to none
        lsm_fields = amrex::Vector<int>{4, -1};
        urban_fields = amrex::Vector<int>{4, -1};
    }

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

            rad_fields.resize(m_nlevs);
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

        // Create new level for all registered fields
        for (int field=0; field < fields.size(); field++) {
            if (fields[field].size() < m_nlevs) {
                fields[field].resize(m_nlevs);
                fields[field][lev] = std::make_unique<amrex::MultiFab>(m_ba2d[lev], m_dmap[lev], 1, amrex::IntVect(1,1,0));
                fields[field][lev]->setVal(0.0);
            }
        }

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

            for (int field=0; field < fields.size(); field++) {
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
     * Sets the urban or land data pointers that are used to compute weight averages
     */
    void set_model_data(const int lev, const amrex::Vector<amrex::MultiFab*> model_data, const amrex::Vector<std::string>& data_names, SurfaceModelType type)
    {
        if (type == SurfaceModelType::LAND) {
            lsm_data_lev[lev].resize(model_data.size());
            if (lev == 0) {
                // only set names on first level, since they are the same for all levels
                m_lsm_names.resize(model_data.size());
                AMREX_ALWAYS_ASSERT(data_names.size() == model_data.size());
                for (int i = 0; i < data_names.size(); i++) {
                    m_lsm_names[i] = data_names[i];
                }
            }

            for (int i = 0; i < lsm_data_lev[lev].size(); i++) {
                AMREX_ALWAYS_ASSERT(model_data[i]);
                lsm_data_lev[lev][i] = model_data[i];
            }
            m_use_land = true;
        } else if (type == SurfaceModelType::URBAN) {
            urban_data_lev[lev].resize(model_data.size());
            if (lev == 0) {
                // only set names on first level, since they are the same for all levels
                m_urban_names.resize(model_data.size());
                AMREX_ALWAYS_ASSERT(data_names.size() == model_data.size());
                for (int i = 0; i < data_names.size(); i++) {
                    m_urban_names[i] = data_names[i];
                }
            }
            for (int i = 0; i < urban_data_lev[lev].size(); i++) {
                AMREX_ALWAYS_ASSERT(model_data[i]);
                urban_data_lev[lev][i] = model_data[i];
            }
            m_use_urban = true;
        }
    }

    /**
     * Adds model flux data to the fields that can be weight averaged.
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

        const int data_size = model_data->size();
        model_data->resize(data_size + model_fluxes.size());
        if (lev == 0) {
            model_names->resize(data_size + model_fluxes.size());
        }

        for (int i = 0; i < model_fluxes.size(); ++i) {
            AMREX_ALWAYS_ASSERT(model_fluxes[i]);
            (*model_data)[data_size + i] = model_fluxes[i];
            if (lev == 0) {
                (*model_names)[data_size + i] = flux_names[i];
            }
        }
    }

    /**
     * Sets the fields to average for each model. When using fluxes, the first five elements map to
     *   u-stress, v-stress, heat flux, vapor flux, tsurf.
     * Otherwise, the first four elements map to ustar, tstar, qstar, tsurf.
     * If a field is set to -1, then it is skipped in the weight average.
     */
    void set_model_fields(SurfaceModelType type, const amrex::Vector<int> &fields, const bool use_fluxes = true)
    {
        int nfields = fields.size();
        if (!use_fluxes) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nfields >= 4, "Must have at least four fields to average (ustar, tstar, qstar, tsurf)");
        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nfields >= 5, "Must have at least five fields to average (u-stress, v-stress, hfx, qfx, tsurf)");
        }

        // make sure both models are using fluxes, or u*,t*,q* directly
        AMREX_ALWAYS_ASSERT(use_fluxes == m_export_fluxes);

        if (type == SurfaceModelType::LAND) {
            lsm_fields = amrex::Vector<int>(fields);
        } else if (type == SurfaceModelType::URBAN) {
            urban_fields = amrex::Vector<int>(fields);
        }
    }

    /**
     * Apply the weighted average to a pair of land and urban data, regardless of whether they are in the field list
     * If one of the pointers is null, then it is skipped. Both pointers cannot be null.
     * This applies weights to ALL elements in the multifab.
     * lsm_data and urban_data can have different underlying grids - weights are applied to each separately
     */
    void apply_weight_average(int lev, amrex::MultiFab *lsm_data, amrex::MultiFab *urban_data);

    /**
     * Calculates the weight average and applies the weights to all the fields to average for each model
     */
    void calculate_weight_average(int lev, amrex::MultiFab* const urban_frac);

    bool are_fluxes() { return m_export_fluxes; }

    amrex::MultiFab* get_ustar(const int& lev) {
        return u_star[lev].get();
    }

    amrex::MultiFab* get_tstar(const int& lev) {
        return t_star[lev].get();
    }

    amrex::MultiFab* get_qstar(const int& lev) {
        return q_star[lev].get();
    }

    amrex::MultiFab* get_tsurf(const int& lev) {
        return t_surf[lev].get();
    }

    amrex::MultiFab* get_wavg_factors(const int &lev) {
        return wavg[lev].get();
    }

    void write_output(const int &finest_lev, const amrex::Real &time, const std::string &plot_prefix, const amrex::Vector<int> &level_steps, const amrex::Vector<amrex::IntVect> &ref_ratio);

    void register_field_map(std::string name, const std::pair<int, int> &lsm_urb_map, bool fill_boundary=false);

    void register_field_map(std::string name, amrex::Vector<amrex::MultiFab*> &lsm_lev_mf, amrex::Vector<amrex::MultiFab*> &urb_lev_mf, bool fill_boundary=false);

    amrex::MultiFab* get_field(const std::string &name, int lev = 0) {
        if (auto field = fieldmap.find(name); field != fieldmap.end()) {
            AMREX_ASSERT(field->second.mf_ind != -1);
            return fields[field->second.mf_ind][lev].get();
        } else {
            return nullptr;
        }
    }

    const amrex::Vector<const amrex::MultiFab*> get_radiation_fields(int lev = 0) {
        // check if we need to build the radiation list (for first time step,)
        if (rad_fields[lev].size() == 0) {
            build_radiation_fieldlist(lev);
        }
        return rad_fields[lev];
    }

    // utility to skip to next line in Header
    //  -- taken from ERF_Checkpoint
    static void GotoNextLine (std::istream& is);

    void WriteCheckpoint(const std::string &checkpointname);

    void ReadCheckpoint(const std::string &checkpointname);

private:

    // Internal function to update weights using urban fraction
    void calculate_simple_average(int lev, amrex::MultiFab* const urban_frac);

    void build_radiation_fieldlist(int lev) {
        amrex::Print() << "  --- building radiation field list" << std::endl;

        // loop over fields and check against names of radiation variables
        for (int i = 0; i < radnames.size(); i++) {
            amrex::Print() << "    - Checking for radiation variable '" << radnames[i] << "'" << std::endl;

            auto mf = get_field(radnames[i], lev);
            if (mf) {
                rad_fields[lev].push_back(mf);
            } else {
                amrex::Print() << "     - WARNING: did not find radiation variable, setting to null to use default RRTMGP value!" << std::endl;
                rad_fields[lev].push_back(nullptr);
            }
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rad_fields[lev].size() == radnames.size(), "Did not find all radiation fields");
    }

    void weight_average_fields(int lev, amrex::MultiFab* const urban_frac);

    bool is_field_mapped(int lev, SurfaceModelType type, int field_idx,
                         const amrex::MultiFab* mf) const;

    int m_nlevs;

    bool m_use_urban = false; // whether we have an urban model
    bool m_use_land = false; // whether we have a LSM

    bool m_export_fluxes = true;

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

    inline static const std::vector<std::string> radnames = {"tskin", "emiss", "albedo_vis", "albedo_nir", "albedo_vis_diff", "albedo_nir_diff"};

    struct Field {
        std::pair<int, int> map;
        amrex::Vector<amrex::MultiFab*> lsm_ptr;
        amrex::Vector<amrex::MultiFab*> urb_ptr;
        int mf_ind;
        bool fill_bound;
    };
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> fields; // NOTE: this is stored in inverse order: (field, nlev)
    std::unordered_map<std::string, Field> fieldmap;

    amrex::Vector<std::string> m_lsm_names;
    amrex::Vector<std::string> m_urban_names;

    amrex::Vector<amrex::Vector<const amrex::MultiFab*>> rad_fields; // pointers to radiation MFs at each level
    
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

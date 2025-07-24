#ifndef ERF_SURFACEMODEL_H
#define ERF_SURFACEMODEL_H

#include "AMReX_Geometry.H"
#include "AMReX_ParmParse.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_IntVect.H"
#include "AMReX_MultiFab.H"

#include <ERF_DataStruct.H>

enum SurfaceModelType { LAND = 0, URBAN = 1};


class SurfaceModel
{
public:

    explicit SurfaceModel(int nlevs,
                          const amrex::Vector<amrex::BoxArray> &ba,
                          const amrex::Vector<amrex::Geometry> &geom,
                          const amrex::DistributionMapping &dm,
                          const SolverChoice &solverChoice) 
        : m_nlevs(nlevs), m_ba(ba), m_geom(geom), m_dmap(dm) {
        
        u_star.resize(m_nlevs);
        t_star.resize(m_nlevs);
        q_star.resize(m_nlevs);
        t_surf.resize(m_nlevs);

        wavg.resize(m_nlevs);


        for (int lev = 0; lev < m_nlevs; lev++)
        {
            // Create a 2D ba, dm, & ghost cells
            amrex::BoxArray ba = m_ba[lev];
            amrex::BoxList bl2d = ba.boxList();
            for (auto& b : bl2d) {
                b.setRange(2, 0);
            }
            amrex::BoxArray ba2d(std::move(bl2d));


            amrex::IntVect ng(1,1,0);
            // Create output arrays for surface boundary condition
            u_star[lev] = std::make_unique<amrex::MultiFab>(ba2d, m_dmap, 2, ng);
            t_star[lev] = std::make_unique<amrex::MultiFab>(ba2d, m_dmap, 1, ng);
            q_star[lev] = std::make_unique<amrex::MultiFab>(ba2d, m_dmap, 1, ng);
            t_surf[lev] = std::make_unique<amrex::MultiFab>(ba2d, m_dmap, 1, ng);

            u_star[lev]->setVal(0.0);
            t_star[lev]->setVal(0.0);
            q_star[lev]->setVal(0.0);
            t_surf[lev]->setVal(300.0);

            // Initialize storage for weight average factors
            wavg[lev] = std::make_unique<amrex::MultiFab>(ba2d, m_dmap, 2, ng);
            //  - default to land=1, urban=0 (i.e, 100% land, no urban fraction)
            wavg[lev]->setVal(1.0, SurfaceModelType::LAND, 1, 0);
            wavg[lev]->setVal(0.0, SurfaceModelType::URBAN, 1, 0);
        }

        lsm_data_lev.resize(m_nlevs);
        urban_data_lev.resize(m_nlevs);

        // set initial fields for each model to none
        lsm_fields = amrex::Vector<int>{4, -1};
        urban_fields = amrex::Vector<int>{4, -1};
    }

    /**
     * Sets the urban or land data pointers that are used to compute weight averages
     */
    void set_model_data(int lev, amrex::Vector<amrex::MultiFab*> model_data, SurfaceModelType type)
    {
        if (type == SurfaceModelType::LAND) {
            lsm_data_lev[lev].resize(model_data.size());
            for (int i = 0; i < lsm_data_lev[lev].size(); i++) {
                AMREX_ALWAYS_ASSERT(model_data[i]);
                lsm_data_lev[lev][i] = model_data[i];
            }
            m_use_land = true;
        } else if (type == SurfaceModelType::URBAN) {
            urban_data_lev[lev].resize(model_data.size());
            for (int i = 0; i < urban_data_lev[lev].size(); i++) {
                AMREX_ALWAYS_ASSERT(model_data[i]);
                urban_data_lev[lev][i] = model_data[i];
            }
            m_use_urban = true;
        }
    }

    /**
     * Sets the fields to average for each model. The first four elements map to the surface boundary outputs
     *   ustar, tstar, qstar, tsurf
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

private:

    // Internal function to update weights using urban fraction
    void calculate_simple_average(int lev, amrex::MultiFab* const urban_frac);

    int m_nlevs;

    bool m_use_urban = false; // whether we have an urban model
    bool m_use_land = false; // whether we have a LSM

    bool m_export_fluxes = true;

    // internal flag to indicate if the weight average has been calculated this time step
    bool m_weights_updated = false;
    
    amrex::Vector<amrex::BoxArray> m_ba;
    amrex::Vector<amrex::Geometry> m_geom;
    amrex::DistributionMapping m_dmap;

    // Pointers to LSM and Urban model data
    amrex::Vector<amrex::Vector<amrex::MultiFab*>>  urban_frac_lev; // Urban fraction for each level [0-1.0]
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
    // First 4 are ustar, tstar, qstar, tsurf
    // TODO: make this more flexible for different types of models
    amrex::Vector<int> lsm_fields;
    amrex::Vector<int> urban_fields;
    
};

#endif
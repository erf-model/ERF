#include "ERF_SurfaceModel.H"

#include "AMReX_Box.H"
#include "AMReX_MFIter.H"
#include <AMReX_PlotFileUtil.H>
#include "ERF_TileNoZ.H"

using namespace amrex;

void SurfaceModel::apply_weight_average(int lev, amrex::MultiFab *lsm_data, amrex::MultiFab *urban_data)
{
    bool valid_land = (lsm_data != nullptr);
    bool valid_urban = (urban_data != nullptr);

    AMREX_ASSERT_WITH_MESSAGE(!valid_land && !valid_urban, "Need at least one pointer to apply weights");

    // make sure the weights have been calculated before applying them
    AMREX_ALWAYS_ASSERT(m_weights_updated);

    // Apply weights separately to support different underlying grids between lsm_data and urban_data
    if (valid_land) {
        for (MFIter mfi(*lsm_data, TileNoZ()); mfi.isValid(); ++mfi)
        {
            Box tbx = mfi.tilebox();

            auto lsm_arr = lsm_data->array(mfi);
            auto weights_arr = wavg[lev]->const_array(mfi);

            // TODO: can amrex::Mult be used here instead?
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                lsm_arr(i, j, k) *= weights_arr(i, j, 0, SurfaceModelType::LAND);
            });

            // TODO: do we always want to update the boundaries?
            lsm_data->FillBoundary(m_geom[lev].periodicity());
        }
    }

    if (valid_urban) {
        for (MFIter mfi(*urban_data, TileNoZ()); mfi.isValid(); ++mfi)
        {
            Box tbx = mfi.tilebox();

            auto urb_arr = urban_data->array(mfi);
            auto weights_arr = wavg[lev]->const_array(mfi);

            // TODO: can amrex::Mult be used here instead?
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                urb_arr(i, j, k) *= weights_arr(i, j, 0, SurfaceModelType::URBAN);
            });

            // TODO: do we always want to update the boundaries?
            urban_data->FillBoundary(m_geom[lev].periodicity());
        }
    }
}


void SurfaceModel::calculate_weight_average(int lev, amrex::MultiFab* const urban_frac)
{
    AMREX_ASSERT_WITH_MESSAGE(m_use_land || m_use_urban, "Must have at least one land or urban model enabled for weighted average");

    // Calculate weighted averages between Land and Urban
    calculate_simple_average(lev, urban_frac);

    const bool use_urban = m_use_urban;
    const bool use_land = m_use_land;
    const bool use_fluxes = m_export_fluxes;

    // Reset surface fluxes
    u_star[lev]->setVal(0.0);
    t_star[lev]->setVal(0.0);
    q_star[lev]->setVal(0.0);

    amrex::MultiFab* const outputs[] = {u_star[lev].get(), t_star[lev].get(), q_star[lev].get(), t_surf[lev].get()};

    const int nfields = (m_export_fluxes) ? 5 : 4;

    // Output weighted surface fluxes into ustar, tstar, qstar and surface temperature into tsurf
    // TODO: make sure grids of urban and LSM inputs match
    for (int field=0; field < nfields; field++)
    {
        int output_field = (use_fluxes && field > 0) ? field - 1 : field;
        int comp = (use_fluxes && field < 2) ? field : 0;

        // whether we have a valid LSM multifab
        bool valid_land = (use_land &&
                           lsm_fields[field] != -1 &&
                           lsm_data_lev[lev][lsm_fields[field]]);

        // whether we have a valid urban multifab
        bool valid_urban = (use_urban &&
                            urban_fields[field] != -1 &&
                            urban_data_lev[lev][urban_fields[field]]);

        for (MFIter mfi(*u_star[lev], TileNoZ()); mfi.isValid(); ++mfi)
        {
            Box tbx = mfi.tilebox();

            auto weights_arr = wavg[lev]->const_array(mfi);

            // Outputs for surface boundary condition
            auto output_arr = outputs[output_field]->array(mfi);
            auto lsm_data_arr = (valid_land) ? lsm_data_lev[lev][lsm_fields[field]]->const_array(mfi) : Array4<const Real>{};
            auto urban_data_arr = (valid_urban) ? urban_data_lev[lev][urban_fields[field]]->const_array(mfi) : Array4<const Real>{};

            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real land = (lsm_data_arr) ? lsm_data_arr(i, j, -1) * weights_arr(i, j, 0, SurfaceModelType::LAND) : 0.0;
                Real urb  = (urban_data_arr) ? urban_data_arr(i, j, k) * weights_arr(i, j, 0, SurfaceModelType::URBAN) : 0.0;

                output_arr(i, j, k, comp) = land + urb;

            });
        }

        outputs[output_field]->FillBoundary(comp, 1, m_geom[lev].periodicity());
    }

    // Apply the weight fractions to the rest of the LSM and urban fields (if any)
    if (use_land) {
        for (int field = nfields; field < lsm_fields.size(); field++) {
            if (lsm_fields[field] == -1) continue;
            for (MFIter mfi(*lsm_data_lev[lev][lsm_fields[field]], TileNoZ()); mfi.isValid(); ++mfi)
            {
                Box tbx = mfi.tilebox();

                auto weights_arr = wavg[lev]->const_array(mfi);

                // whether we have a valid LSM multifab
                bool valid_land = (use_land &&
                                   lsm_fields[field] != -1 &&
                                   lsm_data_lev[lev][lsm_fields[field]]);
                auto lsm_data_arr = (valid_land) ? lsm_data_lev[lev][lsm_fields[field]]->array(mfi) : Array4<Real>{};
            
                ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    lsm_data_arr(i, j, k) *= weights_arr(i, j, 0, SurfaceModelType::LAND);
                });
            }
        }
    }

    if (use_urban) {
        for (int field = nfields; field < urban_fields.size(); field++) {
            if (urban_fields[field] == -1) continue;
            for (MFIter mfi(*urban_data_lev[lev][urban_fields[field]], TileNoZ()); mfi.isValid(); ++mfi)
            {
                Box tbx = mfi.tilebox();

                auto weights_arr = wavg[lev]->const_array(mfi);

                // whether we have a valid urban multifab
                bool valid_urban = (use_urban &&
                                    urban_fields[field] != -1 &&
                                    urban_data_lev[lev][urban_fields[field]]);
                auto urban_data_arr = (valid_urban) ? urban_data_lev[lev][urban_fields[field]]->array(mfi) : Array4<Real>{};

                ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    urban_data_arr(i, j, k) *= weights_arr(i, j, 0, SurfaceModelType::URBAN);
                });
            }

            if (urban_data_lev[lev][urban_fields[field]]) {
                urban_data_lev[lev][urban_fields[field]]->FillBoundary(m_geom[lev].periodicity());
            }
        }
    }
}

void SurfaceModel::calculate_simple_average(int lev, amrex::MultiFab* const urban_frac)
{
    for (MFIter mfi(*wavg[lev], TileNoZ()); mfi.isValid(); ++mfi)
    {
        Box tbx = mfi.tilebox();
        // weights are defined only on k = 0
        tbx = tbx.makeSlab(2, 0);

        auto urban_frac_arr = urban_frac->const_array(mfi);
        auto weights_arr = wavg[lev]->array(mfi);

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // weights are proportional to urban fraction coverage in the current cell
            weights_arr(i, j, 0, SurfaceModelType::URBAN) = urban_frac_arr(i, j, 0);
            weights_arr(i, j, 0, SurfaceModelType::LAND) = 1.0 - urban_frac_arr(i, j, 0);
        });
    }

    m_weights_updated = true;
}


void SurfaceModel::write_output(int lev, const amrex::Real time, const std::string plot_prefix, const int level_step)
{
    std::string plotfilename = amrex::Concatenate(plot_prefix + "2D_", level_step, 5);

    const int nfields = (m_export_fluxes) ? 5 : 4;
    IntVect ng(0, 0, 0);

    amrex::MultiFab* const outputs[] = {u_star[lev].get(), t_star[lev].get(), q_star[lev].get(), t_surf[lev].get()};
    MultiFab fab(outputs[0]->boxArray(), m_dmap[lev], nfields, ng);


    amrex::Vector<std::string> varnames(nfields);

    // Output weighted surface fluxes into ustar, tstar, qstar and surface temperature into tsurf
    // TODO: make sure grids of urban and LSM inputs match
    for (int field=0; field < nfields; field++)
    {
        int output_field = (m_export_fluxes && field > 0) ? field - 1 : field;
        int comp = (m_export_fluxes && field < 2) ? field : 0;

        MultiFab::Copy(fab, *(outputs[output_field]), comp, field, 1, ng);

        varnames[field] = field_names[field];
    }

    amrex::WriteSingleLevelPlotfile(plotfilename, fab, varnames, m_geom2d[lev], time, level_step);
}

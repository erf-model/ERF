/**
 * \file ERF_DustCoupling.cpp
 *
 * ERF's side of the dust coupling: making the dust model erf.dust.model
 * selects, setting it up, and handing it the state after each step of level 0
 * (ERF_DustModel.H).
 */

#include <ERF.H>
#include <ERF_DustModelRegistry.H>

#include <optional>

using namespace amrex;

void
ERF::initialize_dust_model ()
{
    m_dust_model = make_dust_model_from_inputs();
    if (!m_dust_model) { return; }

#ifdef ERF_USE_FIRE
    std::optional<DustFire> fire;
    if (m_fire_coupler) {
        fire.emplace(DustFire{m_fire_coupler->engine(), m_fire_coupler->fire_grid(), m_fire_coupler->fields_time()});
    }
#endif

    DustInit init{*this, m_SurfaceLayer[Orientation::zlo()].get(), *z_phys_nd[0],
                  restart_chkfile, istep[0], t_new[0]
#ifdef ERF_USE_FIRE
                  , fire ? &*fire : nullptr
#endif
                 };
    m_dust_model->initialize(init);

    // Dust starts at zero; on a restart the scalars come from the checkpoint
    if (restart_chkfile.empty()) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            vars_new[lev][Vars::cons].setVal(0.0, RhoDust_comp, NSCALARS_DUST, vars_new[lev][Vars::cons].nGrowVect());
            vars_old[lev][Vars::cons].setVal(0.0, RhoDust_comp, NSCALARS_DUST, vars_old[lev][Vars::cons].nGrowVect());
        }
    }
}

void
ERF::advance_dust_model (int step, double time)
{
    if (!m_dust_model) { return; }

#ifdef ERF_USE_FIRE
    std::optional<DustFire> fire;
    if (m_fire_coupler) {
        fire.emplace(DustFire{m_fire_coupler->engine(), m_fire_coupler->fire_grid(), m_fire_coupler->fields_time()});
    }
#endif

    DustStep ds{step, time, dt[0], geom[0],
                vars_new[0][Vars::cons], vars_new[0][Vars::xvel], vars_new[0][Vars::yvel], vars_new[0][Vars::zvel],
                (z_phys_cc.size() > 0) ? z_phys_cc[0].get() : nullptr,
                (SFS_q1fx3_lev.size() > 0) ? SFS_q1fx3_lev[0].get() : nullptr,
                m_SurfaceLayer[Orientation::zlo()].get()
#ifdef ERF_USE_FIRE
                , fire ? &*fire : nullptr
#endif
               };
    m_dust_model->advance(ds);
}

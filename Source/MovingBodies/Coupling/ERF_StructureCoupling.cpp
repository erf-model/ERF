/**
 * @file ERF_StructureCoupling.cpp
 * @brief ERF's side of the structure coupling: the start-up and advance calls of the
 *        structure model, and the check of its support against the finer level
 *
 * The other calls (momentum sources, checkpoint, plotfile) sit where ERF does that
 * work; see the call order in ERF_StructureModel.H.
 */
#include "ERF.H"
#include "ERF_StructureModel.H"
#include "ERF_ActuatorSampling.H"

#include <cmath>
#include <optional>

using namespace amrex;

void
ERF::initialize_structure_model ()
{
    AMREX_ALWAYS_ASSERT(m_structure_model);

    const int lev = m_structure_model->level(finest_level);
    if (lev < 0 || lev > finest_level) {
        Abort(std::string("structure model \"") + m_structure_model->name() + "\" chose level "
              + std::to_string(lev) + ", but this run has levels 0 to " + std::to_string(finest_level));
    }

    if (!restart_chkfile.empty() && !m_structure_model->supports_restart()) {
        Abort(std::string("structure model \"") + m_structure_model->name()
              + "\" does not support restart from a checkpoint (erf.restart = " + restart_chkfile + ")");
    }

    // The heights and Jacobians of a terrain-following mesh; nullptr on a uniform-dz mesh,
    // where detJ_cc holds ones (or, with EB, the volume fractions) and no physical heights
    const bool fitted = (SolverChoice::mesh_type != MeshType::ConstantDz);
    const MultiFab* zp   = fitted ? z_phys_nd[lev].get() : nullptr;
    const MultiFab* detJ = fitted ? detJ_cc[lev].get()   : nullptr;

    // An immersed terrain stands on a flat mesh: hand the model the terrain's height from the
    // surface the immersed boundary was built from, at the nodes of its level, when it asks
    FArrayBox surface;
    const FArrayBox* surface_ptr = nullptr;
    if (m_structure_model->wants_immersed_terrain_surface() &&
        solverChoice.terrain_type == TerrainType::ImmersedForcing)
    {
        Box terrain_bx(surroundingNodes(geom[lev].Domain()));
        terrain_bx.grow(3);
        surface.resize(makeSlab(terrain_bx, 2, 0), 1, The_Pinned_Arena());
        prob->init_terrain_surface(geom[lev], surface, Real(0.0));
        // Some sources fill the surface with device kernels; the model reads it on the host
        Gpu::streamSynchronize();
        surface_ptr = &surface;
    }

    const TurbChoice& tc = solverChoice.turbChoice[lev];
    const StructureInit init{lev, t_new[lev], geom[lev], vars_new[lev][Vars::xvel], zp, detJ,
                             surface_ptr, structure_unit_map_factors(lev), tc.rans_type == RANSType::kEqn, tc.Cmu0,
                             restart_chkfile};
    m_structure_model->initialize(init);
}

bool
ERF::structure_unit_map_factors (int lev) const
{
    // The actuator core places its kernels at prob_lo + i dx with no map factors
    bool unit = true;
    for (int i = 0; i < MapFacType::num && unit; ++i) {
        const MultiFab& mf = *mapfac[lev][i];
        unit = (mf.min(0) == Real(1.0) && mf.max(0) == Real(1.0));
    }
    return unit;
}

void
ERF::advance_structure_model (int lev, double time, Real dt_lev,
                              const MultiFab& U_old, const MultiFab& V_old,
                              const MultiFab& W_old, const MultiFab& S_old)
{
    AMREX_ALWAYS_ASSERT(m_structure_model);

    const bool fitted = (SolverChoice::mesh_type != MeshType::ConstantDz);

    // The fire's fields, when the run has a fire, at the time they belong to
    std::optional<StructureFire> fire;
#ifdef ERF_USE_FIRE
    if (m_fire_coupler) {
        const FireEngine& engine = m_fire_coupler->engine();
        fire.emplace(StructureFire{m_fire_coupler->fire_grid().geom, m_fire_coupler->level(),
                                   m_fire_coupler->fields_time(),
                                   engine.heat_release(), engine.level_set(),
                                   [&engine] () { return engine.burned_area(); }});
    }
#endif

    const StructureStep step{lev, finest_level, time, dt_lev, geom[lev],
                             U_old, V_old, W_old, S_old,
                             fitted ? z_phys_nd[lev].get() : nullptr, fitted ? detJ_cc[lev].get() : nullptr,
                             structure_unit_map_factors(lev), fire ? &*fire : nullptr};
    m_structure_model->advance(step);

    // The force goes into level lev only. Where level lev+1 covers it, averaging down would
    // overwrite the forced cells with unforced fine data, so the run stops instead. Level
    // lev+1 holds every finer level (proper nesting), so it is the only one to check.
    if (lev >= finest_level) { return; }
    const std::vector<RealBox> support = m_structure_model->support();
    if (support.empty()) { return; }

    const Geometry& fg = geom[lev+1];
    const Box& domain = fg.Domain();
    const auto plo = fg.ProbLoArray();
    const auto dxi = fg.InvCellSizeArray();
    // On a stretched or terrain-following mesh a height lies in a cell between these bounds
    erf_actuator::ZBounds zb;
    if (fitted) { zb = erf_actuator::mesh_z_bounds(*z_phys_nd[lev+1], fg); }
    const std::vector<IntVect> shifts = fg.periodicity().shiftIntVect();

    for (const RealBox& rb : support) {
        IntVect lo, hi;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            lo[d] = static_cast<int>(std::floor((rb.lo(d) - plo[d]) * dxi[d]));
            hi[d] = static_cast<int>(std::floor((rb.hi(d) - plo[d]) * dxi[d]));
        }
        if (fitted) {
            lo[2] = static_cast<int>(std::floor((rb.lo(2) - zb.ground_hi) / zb.dz_hi));
            hi[2] = static_cast<int>(std::floor((rb.hi(2) - zb.ground_lo) / zb.dz_lo));
        }
        const Box needed(lo, hi);
        for (const IntVect& s : shifts) {
            const Box b = (needed + s) & domain;
            if (b.ok() && grids[lev+1].intersects(b)) {
                Abort(std::string("structure model \"") + m_structure_model->name() + "\" couples to level "
                      + std::to_string(lev) + ", but level " + std::to_string(lev+1)
                      + " covers part of its support (x " + std::to_string(rb.lo(0)) + " to " + std::to_string(rb.hi(0))
                      + ", y " + std::to_string(rb.lo(1)) + " to " + std::to_string(rb.hi(1))
                      + ", z " + std::to_string(rb.lo(2)) + " to " + std::to_string(rb.hi(2)) + " m)."
                      + " Its force would be lost there when level " + std::to_string(lev+1)
                      + " is averaged down; couple the model to the finest level that covers its points");
            }
        }
    }
}

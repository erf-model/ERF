/**
 * @file ERF_StructureCoupling.cpp
 * @brief ERF's side of the structure coupling: the start-up call of the structure model
 *
 * The other calls (advance, momentum sources, checkpoint, plotfile) sit where ERF
 * does that work; see the call order in ERF_StructureModel.H.
 */
#include "ERF.H"
#include "ERF_StructureModel.H"

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

    // The actuator core places its kernels at prob_lo + i dx with no map factors
    bool unit_map_factors = true;
    for (int i = 0; i < MapFacType::num && unit_map_factors; ++i) {
        const MultiFab& mf = *mapfac[lev][i];
        unit_map_factors = (mf.min(0) == Real(1.0) && mf.max(0) == Real(1.0));
    }

    const TurbChoice& tc = solverChoice.turbChoice[lev];
    const StructureInit init{lev, t_new[lev], geom[lev], vars_new[lev][Vars::xvel], zp, detJ,
                             surface_ptr, unit_map_factors, tc.rans_type == RANSType::kEqn, tc.Cmu0,
                             restart_chkfile};
    m_structure_model->initialize(init);
}

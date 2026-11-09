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
    m_structure_src_ba = grids[lev];
    m_structure_src_dm = dmap[lev];

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

namespace {

/**
 * Restrict a face source from a level onto the faces normal to dir of the next coarser level,
 * keeping its integral: coarse face I takes the fine faces r I + m, |m| < r, along dir with the
 * weights (r - |m|) / r^2 (the share of each fine face's control volume inside the coarse
 * face's), averaged over the fine faces across it. On a uniform mesh the coarse source then
 * gives the air the same impulse as the fine one. fine must have r - 1 ghost faces along dir,
 * filled where another box of its level holds them and zero outside the level.
 */
void
restrict_face_source (const MultiFab& fine, MultiFab& crse, int dir, const IntVect& rr, const Geometry& cgeom)
{
    BoxArray cba = fine.boxArray();
    cba.coarsen(rr);
    MultiFab crse_on_fine(cba, fine.DistributionMap(), 1, 0);

    const int rn = rr[dir];
    const int t1 = (dir + 1) % AMREX_SPACEDIM;
    const int t2 = (dir + 2) % AMREX_SPACEDIM;
    const int r1 = rr[t1];
    const int r2 = rr[t2];
    const Real inv_w = Real(1.0) / Real(rn * rn * r1 * r2);

    for (MFIter mfi(crse_on_fine, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const auto f = fine.const_array(mfi);
        const auto c = crse_on_fine.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const IntVect base = IntVect(AMREX_D_DECL(i, j, k)) * rr;
            Real s = Real(0.0);
            for (int m = 1 - rn; m < rn; ++m) {
                const Real w = Real(rn - amrex::Math::abs(m));
                for (int a = 0; a < r1; ++a) {
                    for (int b = 0; b < r2; ++b) {
                        IntVect iv = base;
                        iv[dir] += m;
                        iv[t1]  += a;
                        iv[t2]  += b;
                        s += w * f(iv);
                    }
                }
            }
            c(i, j, k) = s * inv_w;
        });
    }

    crse.setVal(Real(0.0));
    crse.ParallelCopy(crse_on_fine, 0, 0, 1, IntVect(0), IntVect(0), cgeom.periodicity());
}

} // namespace

void
ERF::make_structure_coarse_sources (int lev)
{
    if (static_cast<int>(m_structure_crse_src.size()) < finest_level + 1) {
        m_structure_crse_src.resize(finest_level + 1);
    }
    for (MultiFab& mf : m_structure_crse_src[lev]) { mf.clear(); }

    // With OneWay coupling the coarser levels never see the finer ones, and the force stays
    // on the model's level
    if (!m_structure_model || solverChoice.coupling_type != CouplingType::TwoWay) { return; }
    const int L = m_structure_model->level(finest_level);
    if (lev >= L) { return; }

    // The model holds a source on its level's grids only after an advance() on those grids:
    // before its first step since the start or a restart, and after a regrid of its level,
    // level lev takes this step without it
    if (m_structure_src_ba != grids[L] || m_structure_src_dm != dmap[L]) { return; }

    // The source of the model's last step, on its own level, with the ghost faces the
    // restriction onto level L-1 reads
    std::array<MultiFab,AMREX_SPACEDIM> src;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        IntVect ng(0);
        ng[d] = refRatio(L-1)[d] - 1;
        src[d].define(convert(grids[L], IntVect::TheDimensionVector(d)), dmap[L], 1, ng);
        src[d].setVal(Real(0.0));
    }
    m_structure_model->add_momentum_sources(L, src[0], src[1], src[2]);

    // Down one level at a time to level lev
    for (int l = L - 1; l >= lev; --l) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            src[d].FillBoundary(geom[l+1].periodicity());
            IntVect ng(0);
            if (l > lev) { ng[d] = refRatio(l-1)[d] - 1; }
            MultiFab crse(convert(grids[l], IntVect::TheDimensionVector(d)), dmap[l], 1, ng);
            restrict_face_source(src[d], crse, d, refRatio(l), geom[l]);
            src[d] = std::move(crse);
        }
    }
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        m_structure_crse_src[lev][d] = std::move(src[d]);
    }
}

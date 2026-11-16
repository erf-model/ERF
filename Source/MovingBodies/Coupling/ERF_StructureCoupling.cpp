/**
 * @file ERF_StructureCoupling.cpp
 * @brief ERF's side of the structure coupling: the start-up and advance calls of the
 *        structure model, the check of its support against the finer level, and its
 *        source on the coarser levels with the checkpoint of that source
 *
 * The other calls (momentum sources, checkpoint, plotfile) sit where ERF does that
 * work; see the call order in ERF_StructureModel.H.
 */
#include "ERF.H"
#include "ERF_StructureModel.H"
#include "ERF_ActuatorSampling.H"
#include "ERF_StructureCoarseSource.H"

#include <AMReX_VisMF.H>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
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
    // The run's stop time on the clock of t_new, which starts at zero at start_time
    const double stop = (stop_time < std::numeric_limits<double>::max()) ? stop_time - start_time
                                                                         : std::numeric_limits<double>::max();
    bool fully_compressible = false;
    for (int l = 0; l <= max_level; ++l) {
        fully_compressible = fully_compressible ||
                             (solverChoice.compressibility_type[l] == CompressibilityType::FullyCompressible);
    }
    const StructureInit init{lev, t_new[lev], geom[lev], vars_new[lev][Vars::xvel], zp, detJ,
                             surface_ptr, structure_unit_map_factors(lev), tc.rans_type == RANSType::kEqn, tc.Cmu0,
                             restart_chkfile, fixed_dt[lev], stop, max_step, fully_compressible};
    m_structure_model->initialize(init);

    if (!restart_chkfile.empty()) { read_structure_coarse_source_checkpoint(lev); }
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
ERF::advance_structure_model (int lev, double time, double dt_lev,
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
    // The model's own source now takes over from the one read from the checkpoint
    for (MultiFab& mf : m_structure_chk_src) { mf.clear(); }

    // The force goes into level lev only. Where level lev+1 covers it, averaging down would
    // overwrite the forced cells with unforced fine data, so the run stops instead. Level
    // lev+1 holds every finer level (proper nesting), so it is the only one to check.
    if (lev >= finest_level) { return; }
    const std::vector<RealBox> support = m_structure_model->support();
    if (support.empty()) { return; }

    const Geometry& fg = geom[lev+1];
    // On a stretched or terrain-following mesh a height lies in the cells that level lev+1's own node
    // heights allow, wherever its boxes are (on the ground or aloft). The table of those heights is built
    // when level lev+1's boxes change, and at every step only where the terrain moves
    const erf_actuator::ZBounds* zb = nullptr;
    if (fitted) {
        zb = &m_structure_fine_heights.get(*z_phys_nd[lev+1], fg,
                                           solverChoice.terrain_type == TerrainType::MovingFittedMesh);
    }

    for (const RealBox& rb : support) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (std::isnan(rb.lo(d)) || std::isnan(rb.hi(d)) || rb.lo(d) > rb.hi(d)) {
                Abort(std::string("structure model \"") + m_structure_model->name() + "\": support() returned a box with a NaN"
                      " or with its low end above its high end (direction " + std::to_string(d) + ": "
                      + std::to_string(rb.lo(d)) + " to " + std::to_string(rb.hi(d)) + " m)");
            }
        }
        // the cells of level lev+1 the support may lie in, any periodic image of it included
        if (erf_actuator::region_touches_grids(grids[lev+1], fg, rb, zb)) {
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

namespace {

// Where the checkpoint keeps the model's source on its level, for the coarser levels' first steps
std::string
coarse_source_dir (const std::string& checkpointname)
{
    return checkpointname + "/StructureCoarseSource";
}

constexpr const char* coarse_source_face[3] = {"x_faces", "y_faces", "z_faces"};

// The cell Jacobian the coarse-source restriction and its diagnostic weigh a face with:
// detJ_cc on a stretched or terrain-following mesh, none (J = 1) on a uniform-dz mesh, where
// detJ_cc may hold the EB volume fraction instead
const MultiFab*
coarse_source_jacobian (const Vector<std::unique_ptr<MultiFab>>& detJ_cc, int lev)
{
    return (SolverChoice::mesh_type != MeshType::ConstantDz) ? detJ_cc[lev].get() : nullptr;
}

// Whether two BoxArrays of the same index type cover the same region, however they are chunked
bool
covers_same_region (const BoxArray& a, const BoxArray& b)
{
    return a.ixType() == b.ixType() && a.contains(b) && b.contains(a);
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

    // The model holds a source on its level's grids only after an advance() on those grids.
    // After a restart, until that advance(), ERF has the source the model held when the
    // checkpoint was written, and copies it onto the level's grids below, however they are
    // chunked. Before the model's first step since the start, and after a regrid of its level
    // until its next step, level lev takes this step without it. A regrid of level L after a
    // restart, before the model's first step, that changes the region the level covers drops
    // the checkpointed source too: the source would be missing where the level grew
    if (m_structure_chk_src[0].isDefined() &&
        !covers_same_region(m_structure_chk_src[0].boxArray(),
                            convert(grids[L], IntVect::TheDimensionVector(0)))) {
        for (MultiFab& mf : m_structure_chk_src) { mf.clear(); }
        Print() << "Structure coupling: level " << L << " was regridded to another region before"
                << " the model's first step after the restart; level " << lev
                << " takes this step without the checkpointed source" << std::endl;
    }
    const bool from_model = (m_structure_src_ba == grids[L] && m_structure_src_dm == dmap[L]);
    if (!from_model && !m_structure_chk_src[0].isDefined()) { return; }

    // The source of the model's last step, on its own level
    std::array<MultiFab,AMREX_SPACEDIM> src;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        src[d].define(convert(grids[L], IntVect::TheDimensionVector(d)), dmap[L], 1, 0);
        src[d].setVal(Real(0.0));
        if (!from_model) {
            src[d].ParallelCopy(m_structure_chk_src[d], 0, 0, 1);
        }
    }
    if (from_model) {
        m_structure_model->add_momentum_sources(L, src[0], src[1], src[2]);
    }

    // Down to level lev, keeping the force on terrain-following and stretched meshes too
    Vector<const MultiFab*> detJ(L + 1, nullptr);
    Vector<IntVect> rr(L);
    for (int l = 0; l <= L; ++l) {
        detJ[l] = coarse_source_jacobian(detJ_cc, l);
        if (l < L) { rr[l] = refRatio(l); }
    }
    // With erf.v > 0, the force on level lev and on level L, which must agree
    std::array<Real,AMREX_SPACEDIM> force_L{};
    if (verbose > 0) { force_L = erf_structure::face_source_force(src, geom[L], detJ[L]); }
    const MultiFab* detJ_lev = detJ[lev];

    erf_structure::restrict_face_sources(src, L, lev, grids, dmap, geom, rr, detJ);

    if (verbose > 0) {
        const auto force_lev = erf_structure::face_source_force(src, geom[lev], detJ_lev);
        Print() << std::setprecision(17) << "structure coarse source: level " << lev
                << " force = " << force_lev[0] << " " << force_lev[1] << " " << force_lev[2]
                << " level " << L
                << " force = " << force_L[0] << " " << force_L[1] << " " << force_L[2] << "\n";
    }

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        m_structure_crse_src[lev][d] = std::move(src[d]);
    }
}

void
ERF::write_structure_coarse_source_checkpoint (const std::string& checkpointname) const
{
    if (!m_structure_model || solverChoice.coupling_type != CouplingType::TwoWay) { return; }
    const int L = m_structure_model->level(finest_level);
    if (L <= 0) { return; }

    // The source the coarser levels would take at their next step: the model's, or, before
    // its first step after a restart, the one read from the checkpoint. Neither exists before
    // the model's first step since the start or after a regrid of its level, and the
    // checkpoint then holds none, as the coarser levels' next step takes none
    std::array<MultiFab,AMREX_SPACEDIM> src;
    const std::array<MultiFab,AMREX_SPACEDIM>* to_write = &m_structure_chk_src;
    if (m_structure_src_ba == grids[L] && m_structure_src_dm == dmap[L]) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            src[d].define(convert(grids[L], IntVect::TheDimensionVector(d)), dmap[L], 1, 0);
            src[d].setVal(Real(0.0));
        }
        // add_momentum_sources() is const in StructureModel: asking for the source does not
        // change the model, so writing a checkpoint leaves the run as it was
        m_structure_model->add_momentum_sources(L, src[0], src[1], src[2]);
        to_write = &src;
    } else if (!m_structure_chk_src[0].isDefined() ||
               !covers_same_region(m_structure_chk_src[0].boxArray(),
                                   convert(grids[L], IntVect::TheDimensionVector(0)))) {
        return;
    }

    const std::string dir = coarse_source_dir(checkpointname);
    if (ParallelDescriptor::IOProcessor()) {
        if (!UtilCreateDirectory(dir, 0755)) { CreateDirectoryFailed(dir); }
        std::ofstream level_file(dir + "/Level");
        level_file << L << "\n";
        if (!level_file) { Abort("Could not write " + dir + "/Level"); }
    }
    ParallelDescriptor::Barrier();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        VisMF::Write((*to_write)[d], dir + "/" + coarse_source_face[d]);
    }
}

void
ERF::read_structure_coarse_source_checkpoint (int L)
{
    for (MultiFab& mf : m_structure_chk_src) { mf.clear(); }
    if (solverChoice.coupling_type != CouplingType::TwoWay || L <= 0) { return; }

    // A checkpoint written before the model's first step, or by a version of ERF that did not
    // write it, has none: the coarser levels then take their first step after the restart
    // without the source, as before
    const std::string dir = coarse_source_dir(restart_chkfile);
    int L_chk = -1;
    if (ParallelDescriptor::IOProcessor()) {
        std::ifstream level_file(dir + "/Level");
        if (level_file) { level_file >> L_chk; }
    }
    ParallelDescriptor::Bcast(&L_chk, 1, ParallelDescriptor::IOProcessorNumber());
    if (L_chk < 0) { return; }
    if (L_chk != L) {
        Print() << "Structure coupling: " << dir << " holds the source on level " << L_chk
                << ", but the model now couples to level " << L
                << "; the coarser levels take their first step without it" << std::endl;
        return;
    }
    // On the grids the checkpoint was written with; make_structure_coarse_sources copies it
    // onto the level's grids of the run
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        VisMF::Read(m_structure_chk_src[d], dir + "/" + coarse_source_face[d]);
    }
}

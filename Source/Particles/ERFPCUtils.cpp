#include <ERFPC.H>
#include <ERF_Constants.H>
#include <AMReX_ParticleInterpolators.H>
#include <AMReX_ParticleLocator.H>

#ifdef ERF_USE_PARTICLES

using namespace amrex;

namespace {
    struct ERFPCLevelGeom {
        GpuArray<Real,AMREX_SPACEDIM> plo;
        GpuArray<Real,AMREX_SPACEDIM> dxi;
        int k_max;
        int ref_ratio;  // cumulative vertical refinement ratio
    };
}

void ERFPC::massDensity ( MultiFab&  a_mf,
                          const int& a_lev,
                          const int& a_comp ) const
{
    BL_PROFILE("ERFPC::massDensity()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    const Real inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    a_mf.setVal(0.0);

    ParticleToMesh( *this, a_mf, a_lev,
        [=] AMREX_GPU_DEVICE (  const ERFPC::ParticleTileType::ConstParticleTileDataType& ptr,
                                int i, Array4<Real> const& rho)
        {
            auto p = ptr.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh ( p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE ( const ERFPC::ParticleType&, int)
                {
                    auto mass = ptr.m_rdata[ERFParticlesRealIdxSoA::mass][i];
                    return mass*inv_cell_volume;
                });
        });

    return;
}

/*! \brief Fix k-indices for all particles after AMR regrid */
void ERFPC::FixKIndexAMR (const Vector<std::unique_ptr<MultiFab>>& a_z_phys_nd)
{
    BL_PROFILE("ERFPC::FixKIndexAMR()");

    const int finest = finestLevel();

    // Z-levels for non-uniform vertical grids (level-0 cell interfaces)
    const Real* zlevels = m_zlevels_d.empty() ? nullptr : m_zlevels_d.data();
    const int nz_lev0 = m_zlevels_d.empty() ? 0
                       : static_cast<int>(m_zlevels_d.size()) - 1;

    // Helper lambda: recompute k-indices for all particles on a given level.
    // Uses compute_k_from_z for the initial guess, then refines with the
    // terrain height array (update_location_idata) per tile.
    // z_phys_nd is always allocated (even for flat terrain).
    auto recompute_k_for_level = [&](int lev, int ref_ratio)
    {
        const auto& particles = GetParticles();
        if (lev >= static_cast<int>(particles.size())) { return; }
        if (particles[lev].empty()) { return; }

        const Geometry& geom_lev = m_gdb->Geom(lev);
        const auto plo = geom_lev.ProbLoArray();
        const auto dxi = geom_lev.InvCellSizeArray();
        const int k_max = geom_lev.Domain().bigEnd(AMREX_SPACEDIM-1);

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            int grid = pti.index();
            auto& ptile = ParticlesAt(lev, pti);
            auto& aos = ptile.GetArrayOfStructs();
            auto* p_pbox = aos().data();
            const int np = aos.numParticles();

            auto zheight = (*a_z_phys_nd[lev])[grid].array();

            ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                auto& p = p_pbox[i];
                if (p.id() <= 0) { return; }
                p.idata(ERFParticlesIntIdxAoS::k) =
                    compute_k_from_z(Real(p.pos(AMREX_SPACEDIM-1)),
                                     plo[AMREX_SPACEDIM-1],
                                     dxi[AMREX_SPACEDIM-1],
                                     k_max,
                                     zlevels, nz_lev0, ref_ratio);
                update_location_idata(p, plo, dxi, zheight);
            });
        }
        Gpu::synchronize();
    };

    // No fine levels: recompute k-indices for all particles using L0 geometry
    // and move them to level 0 to avoid orphaned particles with invalid k values.
    if (finest < 1) {

        // Recompute k using level-0 geometry for all particles
        const auto& particles = GetParticles();
        for (int lev = 0; lev < static_cast<int>(particles.size()); lev++) {
            recompute_k_for_level(lev, 1);
        }

        // Full Redistribute demotes all particles to level 0
        Redistribute();

        return;
    }

    // Cumulative refinement ratio to finest level
    int finest_ref = 1;
    for (int lev = 0; lev < finest; lev++) {
        finest_ref *= m_gdb->refRatio(lev)[AMREX_SPACEDIM-1];
    }

    // Step 1: Set idata(0) = finest-level k for all particles
    // Enables Redistribute to place particles correctly on fine grids.
    // Use compute_k_from_z only (no terrain correction) because we need
    // finest-level k regardless of which level's height array is available.
    // Step 3 will correct using terrain after particles are on the right level.
    {
        const Geometry& geom_fine = m_gdb->Geom(finest);
        const auto plo_fine = geom_fine.ProbLoArray();
        const auto dxi_fine = geom_fine.InvCellSizeArray();
        const int k_max_fine = geom_fine.Domain().bigEnd(AMREX_SPACEDIM-1);

        for (int lev = 0; lev <= finest; lev++) {
            const auto& particles = GetParticles();
            if (lev >= static_cast<int>(particles.size())) { continue; }
            if (particles[lev].empty()) { continue; }

            for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
                auto& ptile = ParticlesAt(lev, pti);
                auto& aos = ptile.GetArrayOfStructs();
                auto* p_pbox = aos().data();
                const int np = aos.numParticles();

                ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
                {
                    auto& p = p_pbox[i];
                    if (p.id() <= 0) { return; }
                    p.idata(ERFParticlesIntIdxAoS::k) =
                        compute_k_from_z(Real(p.pos(AMREX_SPACEDIM-1)),
                                         plo_fine[AMREX_SPACEDIM-1],
                                         dxi_fine[AMREX_SPACEDIM-1],
                                         k_max_fine,
                                         zlevels, nz_lev0, finest_ref);
                });
            }
        }
        Gpu::synchronize();
    }

    // Step 2: Full multi-level Redistribute (k-clamping prevents crashes)
    Redistribute();

    // Step 3: Recompute idata(k) using each particle's level geometry.
    // With terrain, the terrain-corrected k may differ from Step 1's flat k,
    // potentially placing particles outside their current level's partial-z box.
    // A full Redistribute is needed to move such particles to the correct level,
    // followed by a second k recomputation for particles that changed level.
    for (int lev = 0; lev <= finest; lev++) {
        int lev_ref = 1;
        for (int l = 0; l < lev; l++) {
            lev_ref *= m_gdb->refRatio(l)[AMREX_SPACEDIM-1];
        }
        recompute_k_for_level(lev, lev_ref);
    }

    // Step 4: Full Redistribute to fix level assignment after terrain correction.
    Redistribute();

    // Step 5: Recompute k for particles that changed level in Step 4.
    for (int lev = 0; lev <= finest; lev++) {
        int lev_ref = 1;
        for (int l = 0; l < lev; l++) {
            lev_ref *= m_gdb->refRatio(l)[AMREX_SPACEDIM-1];
        }
        recompute_k_for_level(lev, lev_ref);
    }

    // Step 6: Per-level Redistribute to fix tile assignment after final k.
    // L0 covers the full domain, so per-level is safe. For fine levels with
    // partial-z boxes, per-level redistribute can fail, so use full.
    Redistribute(0, 0);
    if (finest > 0) {
        Redistribute();
    }

    // Step 7: Final k recomputation.  Steps 4-6 may have moved particles
    // between levels (e.g. from L0 to L1 and back), leaving idata(k) set
    // to the wrong level's geometry.  Recompute once more on each
    // particle's final level to guarantee a consistent k.
    for (int lev = 0; lev <= finest; lev++) {
        int lev_ref = 1;
        for (int l = 0; l < lev; l++) {
            lev_ref *= m_gdb->refRatio(l)[AMREX_SPACEDIM-1];
        }
        recompute_k_for_level(lev, lev_ref);
    }
}

/*! \brief Diagnostic: count particles per level and, for level 1, in the halo.
 *
 * Prints particle counts for all levels. For AMR with >= 3 levels, also
 * counts level-1 particles in halo cells (not covered by level 2).
 *
 * \param[in] a_finest_level Current finest AMR level (>= 0)
 */
void ERFPC::CountParticlesPerLevelAndHalo (int a_finest_level)
{
    BL_PROFILE("ERFPC::CountParticlesPerLevelAndHalo()");

    // Per-level total counts
    amrex::Vector<Long> n_per_lev(a_finest_level + 1, 0);
    for (int lev = 0; lev <= a_finest_level; lev++) {
        n_per_lev[lev] = NumberOfParticlesAtLevel(lev, true, true);
    }
    for (int lev = 0; lev <= a_finest_level; lev++) {
        ParallelDescriptor::ReduceLongSum(n_per_lev[lev]);
    }

    amrex::Print() << "[" << m_name << "] Particle counts per level (after Redistribute): ";
    for (int lev = 0; lev <= a_finest_level; lev++) {
        amrex::Print() << "L" << lev << "=" << n_per_lev[lev];
        if (lev < a_finest_level) { amrex::Print() << " "; }
    }
    amrex::Print() << "\n";

    if (a_finest_level < 2) { return; }

    // Level-1 halo cells not covered by level 2 (1=halo, 0=covered)
    const int lev1 = 1;
    iMultiFab halo_mask = amrex::makeFineMask(ParticleBoxArray(lev1), ParticleDistributionMap(lev1),
                                               IntVect(0), ParticleBoxArray(lev1 + 1), m_gdb->refRatio(lev1),
                                               m_gdb->Geom(lev1).periodicity(), 1, 0);

    const amrex::Geometry& geom1 = m_gdb->Geom(lev1);
    const auto plo = geom1.ProbLoArray();
    const auto dxi = geom1.InvCellSizeArray();
    const amrex::Box& domain1 = geom1.Domain();

    Long n_halo = 0;
    for (ParIterType pti(*this, lev1); pti.isValid(); ++pti) {
        const auto& aos = ParticlesAt(lev1, pti).GetArrayOfStructs();
        const int np = aos.numParticles();
        if (np == 0) { continue; }

        const amrex::Box& box = pti.validbox();
        amrex::Array4<int const> const& mask_arr = halo_mask.const_array(pti);
        auto* p_pbox = aos().data();

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Long> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(np, reduce_data,
            [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
            {
                auto const& p = p_pbox[i];
                if (p.id() <= 0) { return {0L}; }
                int ix = static_cast<int>(amrex::Math::floor((p.pos(0) - plo[0]) * dxi[0]))
                         + domain1.smallEnd(0);
                int iy = static_cast<int>(amrex::Math::floor((p.pos(1) - plo[1]) * dxi[1]))
                         + domain1.smallEnd(1);
                int iz = p.idata(ERFParticlesIntIdxAoS::k);
                iz = amrex::max(domain1.smallEnd(2), amrex::min(iz, domain1.bigEnd(2)));
                amrex::IntVect iv(AMREX_D_DECL(ix, iy, iz));
                if (box.contains(iv) && mask_arr(iv) == 1) {
                    return {1L};
                }
                return {0L};
            });
        n_halo += amrex::get<0>(reduce_data.value(reduce_op));
    }
    ParallelDescriptor::ReduceLongSum(n_halo);

    amrex::Print() << "[" << m_name << "] Level-1 halo particle count (cells not covered by L2): " << n_halo
                   << " (L1 total=" << n_per_lev[1] << ")\n";
}

/*! \brief Extract OOR particles from level a_lev, route to finest covering level,
 *         then redistribute. Avoids multi-level redistribute and excessive MPI traffic.
 *
 * \param[in] a_lev       Fine level with potentially OOR particles
 * \param[in] a_z_phys_nd Terrain height data (all levels)
 */
void ERFPC::ExtractAndRouteOORParticles ( int                                        a_lev,
                                          const Vector<std::unique_ptr<MultiFab>>&   a_z_phys_nd )
{
    BL_PROFILE("ERFPC::ExtractAndRouteOORParticles()");
    amrex::ignore_unused(a_z_phys_nd);

    AMREX_ALWAYS_ASSERT(a_lev > 0);

    ERFParticlesAssignor cell_assignor;
    const int n_levels = a_lev + 1;

    // Z-levels data for k computation (level-0 cell interfaces)
    const Real* zlevels = m_zlevels_d.empty() ? nullptr
                        : m_zlevels_d.data();
    const int nz_lev0 = m_zlevels_d.empty() ? 0
                      : static_cast<int>(m_zlevels_d.size()) - 1;

    // Build per-level ParticleLocators (GPU spatial indices)
    Vector<ParticleLocator<DenseBins<Box>>> locators(n_levels);
    for (int lev = 0; lev < n_levels; lev++) {
        locators[lev].build(ParticleBoxArray(lev), m_gdb->Geom(lev));
    }

    // Collect the per-level AssignGrid objects into a device array
    using AGType = AssignGrid<DenseBinIteratorFactory<Box>>;
    Gpu::PinnedVector<AGType> h_assignors(n_levels);
    for (int lev = 0; lev < n_levels; lev++) {
        h_assignors[lev] = locators[lev].getGridAssignor();
    }
    Gpu::DeviceVector<AGType> d_assignors(n_levels);
    Gpu::copyAsync( Gpu::hostToDevice,
                    h_assignors.begin(),
                    h_assignors.end(),
                    d_assignors.begin());

    // Per-level geometry data for compute_k_from_z on the GPU
    Gpu::PinnedVector<ERFPCLevelGeom> h_lg(n_levels);
    for (int lev = 0; lev < n_levels; lev++) {
        const auto& geom = m_gdb->Geom(lev);
        h_lg[lev].plo = geom.ProbLoArray();
        h_lg[lev].dxi = geom.InvCellSizeArray();
        h_lg[lev].k_max = geom.Domain().bigEnd(AMREX_SPACEDIM-1);
        int ref = 1;
        for (int l = 0; l < lev; l++) {
            ref *= m_gdb->refRatio(l)[AMREX_SPACEDIM-1];
        }
        h_lg[lev].ref_ratio = ref;
    }
    Gpu::DeviceVector<ERFPCLevelGeom> d_lg(n_levels);
    Gpu::copyAsync(Gpu::hostToDevice, h_lg.begin(), h_lg.end(), d_lg.begin());
    Gpu::synchronize();

    auto* ag_ptr = d_assignors.data();
    auto* lg_ptr = d_lg.data();

    // Find one local grid per level for receiving extracted particles
    Gpu::PinnedVector<int> dest_grids(n_levels, 0);
    {
        const int my_proc = ParallelDescriptor::MyProc();
        for (int lev = 0; lev < n_levels; lev++) {
            const auto& dm = ParticleDistributionMap(lev);
            for (int i = 0; i < static_cast<int>(dm.size()); i++) {
                if (dm[i] == my_proc) { dest_grids[lev] = i; break; }
            }
        }
    }

    int finest = a_lev;
    Vector<int> levels_modified(n_levels, 0);
    Vector<Long> n_routed_to(n_levels, 0);

    {
        const int lev = a_lev;
        auto src_lev_grid = locators[lev].getGridAssignor();

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            auto& src_tile = ParticlesAt(lev, pti);
            auto& aos = src_tile.GetArrayOfStructs();
            auto* p_pbox = aos().data();
            const int np = aos.numParticles();
            if (np == 0) { continue; }

            // Create OOR mask on GPU
            Gpu::DeviceVector<int> mask(np);
            auto* mask_ptr = mask.data();

            ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                auto& p = p_pbox[i];
                if (p.id() <= 0) { mask_ptr[i] = 0; return; }
                int grd = src_lev_grid(p, 0, cell_assignor);
                mask_ptr[i] = (grd < 0) ? 1 : 0;
            });

            // Count OOR particles
            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(np, reduce_data,
                [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
                { return {mask_ptr[i]}; });
            int n_oor = amrex::get<0>(reduce_data.value(reduce_op));
            if (n_oor == 0) { continue; }

            // Extract OOR particles to a temporary tile
            ParticleTileType tmp_tile;
            tmp_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());
            tmp_tile.resize(n_oor);

            [[maybe_unused]] int n_copied = amrex::filterParticles(tmp_tile, src_tile, mask_ptr,
                                                                   int(0), int(0), np);
            AMREX_ASSERT(n_copied == n_oor);

            // Find the correct target level for each particle
            Gpu::DeviceVector<int> target_lev_vec(n_oor);
            auto* tgt_ptr = target_lev_vec.data();
            auto* tmp_pbox = tmp_tile.GetArrayOfStructs()().data();

            ParallelFor(n_oor, [=] AMREX_GPU_DEVICE (int i)
            {
                auto& p = tmp_pbox[i];
                int found_lev = -1;

                for (int tl = finest; tl >= 0; tl--) {
                    const auto& lg = lg_ptr[tl];
                    int k = compute_k_from_z( Real(p.pos(AMREX_SPACEDIM-1)),
                                              lg.plo[AMREX_SPACEDIM-1],
                                              lg.dxi[AMREX_SPACEDIM-1],
                                              lg.k_max, zlevels, nz_lev0,
                                              lg.ref_ratio );
                    p.idata(ERFParticlesIntIdxAoS::k) = k;
                    int grd = ag_ptr[tl](p, 0, cell_assignor);
                    if (grd >= 0) { found_lev = tl; break; }
                }
                tgt_ptr[i] = found_lev;
            });
            Gpu::synchronize();

            // Route particles to their target levels
            for (int tl = 0; tl < n_levels; tl++) {
                if (tl == lev) { continue; }

                Gpu::DeviceVector<int> tl_mask(n_oor);
                auto* tl_mask_ptr = tl_mask.data();
                int target = tl;
                ParallelFor(n_oor, [=] AMREX_GPU_DEVICE (int i)
                {
                    tl_mask_ptr[i] = (tgt_ptr[i] == target) ? 1 : 0;
                });

                ReduceOps<ReduceOpSum> tl_reduce_op;
                ReduceData<int> tl_reduce_data(tl_reduce_op);
                using TLReduceTuple = typename decltype(tl_reduce_data)::Type;
                tl_reduce_op.eval(n_oor, tl_reduce_data,
                    [=] AMREX_GPU_DEVICE (int i) -> TLReduceTuple
                    { return {tl_mask_ptr[i]}; });
                int n_to_lev = amrex::get<0>(tl_reduce_data.value(tl_reduce_op));
                if (n_to_lev == 0) { continue; }

                ParticleTileType tl_tile;
                tl_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());
                tl_tile.resize(n_to_lev);
                [[maybe_unused]] int nc = amrex::filterParticles(tl_tile, tmp_tile, tl_mask_ptr,
                                                                 int(0), int(0), n_oor);
                AMREX_ASSERT(nc == n_to_lev);

                // Recompute idata(k) for the actual target level
                auto* tl_pbox = tl_tile.GetArrayOfStructs()().data();
                const auto tl_plo = h_lg[tl].plo;
                const auto tl_dxi = h_lg[tl].dxi;
                int tl_kmax = h_lg[tl].k_max;
                int tl_ref  = h_lg[tl].ref_ratio;
                ParallelFor(n_to_lev, [=] AMREX_GPU_DEVICE (int i)
                {
                    auto& p = tl_pbox[i];
                    p.idata(ERFParticlesIntIdxAoS::k) = compute_k_from_z( Real(p.pos(AMREX_SPACEDIM-1)),
                                                                          tl_plo[AMREX_SPACEDIM-1],
                                                                          tl_dxi[AMREX_SPACEDIM-1],
                                                                          tl_kmax, zlevels, nz_lev0, tl_ref);
                });
                Gpu::synchronize();

                auto& dst_tile = DefineAndReturnParticleTile(tl, dest_grids[tl], 0);
                auto dst_start = static_cast<int>(dst_tile.numParticles());
                dst_tile.resize(dst_start + n_to_lev);
                amrex::copyParticles(dst_tile, tl_tile, int(0), dst_start, n_to_lev);
                levels_modified[tl] = 1;
                n_routed_to[tl] += n_to_lev;
            }

            // Mark OOR particles as dead in the original tile
            ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                if (mask_ptr[i]) { p_pbox[i].id() = -1; }
            });
            Gpu::synchronize();
        }
    }

    // Ensure all ranks agree on which levels were modified
    ParallelDescriptor::ReduceIntMax(levels_modified.data(), n_levels);

    // Per-level redistribute: only for levels that need it.
    // Level a_lev always needs it (particles were advected, some marked dead).
    // Lower levels need it only if particles were routed there.
    for (int lev = 0; lev <= a_lev; lev++) {
        if (lev == a_lev || levels_modified[lev]) {
            Redistribute(lev, lev);
        }
    }
}

#endif

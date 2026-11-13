#include <random>
#include "ERF_SuperDropletPC.H"
#include "ERF_TerrainConversion.H"
#include <AMReX_MultiFabUtil.H>

#ifdef ERF_USE_PARTICLES

using namespace amrex;
using namespace SDPCDefn;

/*! \brief Sets the initial number of super-droplets per cell in a box region with uniform distribution
 *
 * This function initializes the number of super-droplets to be placed in each grid cell
 * based on whether the cell is inside the specified box region. It handles both flat
 * and terrain-following grids and can initialize particles in subgrid regions.
 *
 * \param[out] a_num_sd Integer MultiFab that will contain the number of superdroplets in each grid cell
 * \param[in] a_n_per_cell Number of superdroplets to place per cell
 * \param[in] a_height_ptr Pointer to MultiFab containing terrain height information
 * \param[in] a_box Box region within which to initialize particles
 * \param[in] a_subgrid Flag indicating if the box is smaller than a grid cell
 */
void SuperDropletPC::setNumSDBoxDistribution (int a_lev,
                                              iMultiFab& a_num_sd,
                                              const int a_n_per_cell,
                                              const MFPtr& a_height_ptr,
                                              const RealBox& a_box,
                                              const bool a_subgrid)
{
    BL_PROFILE("SuperDropletPC::setNumSDBoxDistribution()");
    a_num_sd.setVal(0);

    const auto dx = Geom(a_lev).CellSizeArray();
    const auto plo = Geom(a_lev).ProbLoArray();

    // Half-open [lo,hi) on the cell center keeps the injected cell count independent of box-grid alignment.
    const Real blo0 = a_box.lo(0), blo1 = a_box.lo(1), blo2 = a_box.lo(2);
    const Real bhi0 = a_box.hi(0), bhi1 = a_box.hi(1), bhi2 = a_box.hi(2);

    for(MFIter mfi = MakeMFIter(a_lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();
        auto num_superdroplets_arr = a_num_sd[mfi].array();
        if (a_height_ptr) {
            const auto height_arr = (*a_height_ptr)[mfi].array();
            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                bool flag = false;
                if (a_subgrid) {
                    RealBox gridcell( plo[0]+i*dx[0],
                                      plo[1]+j*dx[1],
                                      height_arr(i,j,k),
                                      plo[0]+(i+1)*dx[0],
                                      plo[1]+(j+1)*dx[1],
                                      height_arr(i+1,j+1,k+1) );
                    if (gridcell.contains(a_box)) { flag = true; }
                } else {
                    Real x = plo[0] + (i + myhalf)*dx[0];
                    Real y = plo[1] + (j + myhalf)*dx[1];
                    Real z = Real(0.125) * (height_arr(i,j  ,k  ) + height_arr(i+1,j  ,k  ) +
                                      height_arr(i,j+1,k  ) + height_arr(i+1,j+1,k  ) +
                                      height_arr(i,j  ,k+1) + height_arr(i+1,j  ,k+1) +
                                      height_arr(i,j+1,k+1) + height_arr(i+1,j+1,k+1) );
                    bool in_box = AMREX_D_TERM(   (x >= blo0) && (x < bhi0),
                                               && (y >= blo1) && (y < bhi1),
                                               && (z >= blo2) && (z < bhi2) );
                    if (in_box) { flag = true; }
                }
                if (flag) { num_superdroplets_arr(i,j,k) = a_n_per_cell; }
            });
        } else {
            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                bool flag = false;
                if (a_subgrid) {
                    RealBox gridcell( plo[0]+i*dx[0],
                                      plo[1]+j*dx[1],
                                      plo[2]+k*dx[2],
                                      plo[0]+(i+1)*dx[0],
                                      plo[1]+(j+1)*dx[1],
                                      plo[2]+(k+1)*dx[2] );
                    if (gridcell.contains(a_box)) { flag = true; }
                } else {
                    Real x = plo[0] + (i + myhalf)*dx[0];
                    Real y = plo[1] + (j + myhalf)*dx[1];
                    Real z = plo[2] + (k + myhalf)*dx[2];
                    bool in_box = AMREX_D_TERM(   (x >= blo0) && (x < bhi0),
                                               && (y >= blo1) && (y < bhi1),
                                               && (z >= blo2) && (z < bhi2) );
                    if (in_box) { flag = true; }
                }
                if (flag) { num_superdroplets_arr(i,j,k) = a_n_per_cell; }
            });
        }
    }

    return;
}

/*! Scatters a box-level count of super-droplets over the injection box for
 *  per-box high-multiplicity injection. The positions are drawn on the host with
 *  a deterministic seed (identical on every rank) so the global scatter is
 *  independent of the rank count, then binned to their cells. On a stretched or
 *  terrain-fitted mesh the z-binning uses the physical heights of the cells, and
 *  positions below the terrain are dropped. */
void SuperDropletPC::setNumSDPerBox (int a_lev,
                                     iMultiFab& a_num_sd,
                                     const int a_n_total,
                                     const RealBox& a_box,
                                     const unsigned int a_seed,
                                     const MFPtr& a_height_ptr)
{
    BL_PROFILE("SuperDropletPC::setNumSDPerBox()");
    a_num_sd.setVal(0);
    if (a_n_total <= 0) { return; }

    const auto dxi = Geom(a_lev).InvCellSizeArray();
    const auto plo = Geom(a_lev).ProbLoArray();

    if (m_variable_dz) {
        Gpu::DeviceVector<GpuArray<Real,AMREX_SPACEDIM>> pos_d(a_n_total);
        {
            Vector<GpuArray<Real,AMREX_SPACEDIM>> pos_h(a_n_total);
            std::mt19937 gen(a_seed);
            std::uniform_real_distribution<Real> u(Real(0.0), Real(1.0));
            for (int n = 0; n < a_n_total; ++n) {
                const Real x = a_box.lo(0) + u(gen)*(a_box.hi(0) - a_box.lo(0));
                const Real y = a_box.lo(1) + u(gen)*(a_box.hi(1) - a_box.lo(1));
                const Real z = a_box.lo(2) + u(gen)*(a_box.hi(2) - a_box.lo(2));
                pos_h[n] = {AMREX_D_DECL(x, y, z)};
            }
            Gpu::copy(Gpu::hostToDevice, pos_h.begin(), pos_h.end(), pos_d.begin());
        }
        const auto* pos_ptr = pos_d.data();
        const int n_total = a_n_total;
        const Box& dom = Geom(a_lev).Domain();
        const int k_max = dom.bigEnd(AMREX_SPACEDIM-1) - dom.smallEnd(AMREX_SPACEDIM-1);

        for (MFIter mfi = MakeMFIter(a_lev); mfi.isValid(); ++mfi) {
            const Box& tile_box = mfi.tilebox();
            auto num_superdroplets_arr = a_num_sd[mfi].array();
            const auto z_nd = (*a_height_ptr)[mfi].const_array();
            ParallelFor(n_total, [=] AMREX_GPU_DEVICE (int n) noexcept
            {
                const auto& x = pos_ptr[n];
                const int i = static_cast<int>(amrex::Math::floor((x[0]-plo[0])*dxi[0]));
                const int j = static_cast<int>(amrex::Math::floor((x[1]-plo[1])*dxi[1]));
                if (i < tile_box.smallEnd(0) || i > tile_box.bigEnd(0)) { return; }
                if (j < tile_box.smallEnd(1) || j > tile_box.bigEnd(1)) { return; }
                // the column (i,j) is in this tile, so its heights are in z_nd
                const Real zeta = ERF::ParticlePos::zeta_from_z(x[0], x[1], x[2], plo, dxi, z_nd, k_max);
                const int k = static_cast<int>(amrex::Math::floor((zeta-plo[2])*dxi[2]));
                if (tile_box.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                    Gpu::Atomic::AddNoRet(&num_superdroplets_arr(i,j,k), 1);
                }
            });
        }
        return;
    }

    Gpu::DeviceVector<IntVect> cells_d(a_n_total);
    {
        Vector<IntVect> cells_h(a_n_total);
        std::mt19937 gen(a_seed);
        std::uniform_real_distribution<Real> u(Real(0.0), Real(1.0));
        for (int n = 0; n < a_n_total; ++n) {
            const Real x = a_box.lo(0) + u(gen)*(a_box.hi(0) - a_box.lo(0));
            const Real y = a_box.lo(1) + u(gen)*(a_box.hi(1) - a_box.lo(1));
            const Real z = a_box.lo(2) + u(gen)*(a_box.hi(2) - a_box.lo(2));
            cells_h[n] = IntVect(AMREX_D_DECL(
                static_cast<int>(std::floor((x-plo[0])*dxi[0])),
                static_cast<int>(std::floor((y-plo[1])*dxi[1])),
                static_cast<int>(std::floor((z-plo[2])*dxi[2])) ));
        }
        Gpu::copy(Gpu::hostToDevice, cells_h.begin(), cells_h.end(), cells_d.begin());
    }
    const IntVect* cells_ptr = cells_d.data();
    const int n_total = a_n_total;

    for (MFIter mfi = MakeMFIter(a_lev); mfi.isValid(); ++mfi) {
        const Box& tile_box = mfi.tilebox();
        auto num_superdroplets_arr = a_num_sd[mfi].array();
        ParallelFor(n_total, [=] AMREX_GPU_DEVICE (int n) noexcept
        {
            const IntVect iv = cells_ptr[n];
            if (tile_box.contains(iv)) {
                Gpu::Atomic::AddNoRet(&num_superdroplets_arr(iv[0],iv[1],iv[2]), 1);
            }
        });
    }

    return;
}

/*! Sets the initial number of the super-droplets per cell as a bubble with a uniform distribution */
void SuperDropletPC::setNumSDBubbleDistribution ( int a_lev,
                                                  iMultiFab& a_num_sd, /*!< integer Multifab with number of superdroplets in each grid cell */
                                                  const int a_n_per_cell, /*!< number of superdroplets per cell */
                                                  const MFPtr& a_height_ptr, /*!< terrain */
                                                  const RealBox& a_bubble, /*!< bubble within which to initialize particles */
                                                  const bool a_subgrid /*!< Is a_box smaller than grid cell */)
{
    BL_PROFILE("SuperDropletPC::setNumSDBubbleDistribution()");
    a_num_sd.setVal(0);

    const auto dx = Geom(a_lev).CellSizeArray();
    const auto plo = Geom(a_lev).ProbLoArray();

    for(MFIter mfi = MakeMFIter(a_lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();
        auto num_superdroplets_arr = a_num_sd[mfi].array();
        if (a_height_ptr) {
            const auto height_arr = (*a_height_ptr)[mfi].array();
            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                bool flag = false;
                if (a_subgrid) {
                    RealBox gridcell( plo[0]+i*dx[0],
                                      plo[1]+j*dx[1],
                                      height_arr(i,j,k),
                                      plo[0]+(i+1)*dx[0],
                                      plo[1]+(j+1)*dx[1],
                                      height_arr(i+1,j+1,k+1) );
                    if (gridcell.contains(a_bubble.lo())) { flag = true; }
                } else {
                    Real x = plo[0] + (i + myhalf)*dx[0];
                    Real y = plo[1] + (j + myhalf)*dx[1];
                    Real z = Real(0.125) * (height_arr(i,j  ,k  ) + height_arr(i+1,j  ,k  ) +
                                      height_arr(i,j+1,k  ) + height_arr(i+1,j+1,k  ) +
                                      height_arr(i,j  ,k+1) + height_arr(i+1,j  ,k+1) +
                                      height_arr(i,j+1,k+1) + height_arr(i+1,j+1,k+1) );

                    // Extract bubble params
                    const auto& x_c = a_bubble.lo(); // center
                    const auto& x_r = a_bubble.hi(); // radius

                    Real rad = zero;
                    if (x_r[0] > 0) rad += amrex::Math::powi<2>((x - x_c[0])/x_r[0]);
                    if (x_r[1] > 0) rad += amrex::Math::powi<2>((y - x_c[1])/x_r[1]);
                    if (x_r[2] > 0) rad += amrex::Math::powi<2>((z - x_c[2])/x_r[2]);
                    rad = std::sqrt(rad);

                    if(rad <= one) { flag = true; }
                }
                if (flag) { num_superdroplets_arr(i,j,k) = a_n_per_cell; }
            });
        } else {
            ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                bool flag = false;
                if (a_subgrid) {
                    RealBox gridcell( plo[0]+i*dx[0],
                                      plo[1]+j*dx[1],
                                      plo[2]+k*dx[2],
                                      plo[0]+(i+1)*dx[0],
                                      plo[1]+(j+1)*dx[1],
                                      plo[2]+(k+1)*dx[2] );
                    if (gridcell.contains(a_bubble.lo())) { flag = true; }
                } else {
                    Real x = plo[0] + (i + myhalf)*dx[0];
                    Real y = plo[1] + (j + myhalf)*dx[1];
                    Real z = plo[2] + (k + myhalf)*dx[2];

                    // Extract bubble params
                    const auto& x_c = a_bubble.lo();       // center
                    const auto& x_r = a_bubble.hi();       // radius

                    Real rad = zero;
                    if (x_r[0] > 0) rad += amrex::Math::powi<2>((x - x_c[0])/x_r[0]);
                    if (x_r[1] > 0) rad += amrex::Math::powi<2>((y - x_c[1])/x_r[1]);
                    if (x_r[2] > 0) rad += amrex::Math::powi<2>((z - x_c[2])/x_r[2]);
                    rad = std::sqrt(rad);

                    if(rad <= one) { flag = true; }
                }
                if (flag) { num_superdroplets_arr(i,j,k) = a_n_per_cell; }
            });
        }
    }

    return;
}

/*! Add super-droplets in domain given an initialization type

    + The number of physical particles per cell is computed from the number density;
      if specified (if not specified, it is taken to be 1).
    + The number of super-droplets per cell is computed from the super-droplet number;
      if specified (if not specified, it is set to the particles per cell, whose default value is 1).
*/
void SuperDropletPC::addParticles ( int a_lev,
                                    const MFPtr& a_height_ptr, /*!< terrain */
                                    const SDInitProperties& a_init, /*!< initialization parameters */
                                    bool composite_init )
{
    BL_PROFILE("SuperDropletPC::addParticles");

    const auto dx_h = Geom(a_lev).CellSize();
    const Real cell_volume = dx_h[0]*dx_h[1]*dx_h[2];
    const auto dx = Geom(a_lev).CellSizeArray();
    const auto dxi = Geom(a_lev).InvCellSizeArray();
    const auto plo = Geom(a_lev).ProbLoArray();
    const auto init_volume = a_init.volume();
    const bool subgrid = (init_volume < cell_volume);
    const auto itype = a_init.m_type;
    const Box& dom = Geom(a_lev).Domain();
    const int k_max = dom.bigEnd(AMREX_SPACEDIM-1) - dom.smallEnd(AMREX_SPACEDIM-1);
    const bool variable_dz = m_variable_dz;
    AMREX_ALWAYS_ASSERT(!variable_dz || a_height_ptr);

    // Per-box high-multiplicity injection scatters a box-level count of
    // super-droplets at a fixed multiplicity; the legacy path injects a fixed
    // count per cell.
    const bool perbox = a_init.perBoxInjection();
    // On a stretched or terrain-fitted mesh each cell's counts follow its physical volume
    const bool per_cell_counts = variable_dz && !perbox;
    const SDCellCounts cell_counts = a_init.cellCounts();
    int num_sd_per_cell;
    Real num_par_per_cell;
    Real inject_mult = one;
    int inject_count = 0;
    if (perbox) {
        inject_count = a_init.injectBoxCount();
        inject_mult  = a_init.injectMultiplicity();
        if (inject_count == 0) { return; }
        num_sd_per_cell  = 1;            // unused in per-box mode
        num_par_per_cell = inject_mult;  // unused in per-box mode
        Print() << "    Per-box injection: " << inject_count
                << " super-droplets at multiplicity " << inject_mult << "\n";
    } else if (per_cell_counts) {
        num_sd_per_cell  = 1;            // flags the cells, counted from their volumes below
        num_par_per_cell = one;          // set per cell below
        Print() << "    Numbers of physical particles and super droplets per cell follow the physical cell volume\n";
    } else {
        num_sd_per_cell  = a_init.numSDPerCell(cell_volume);
        num_par_per_cell = a_init.numParticlesPerCell(cell_volume);
        Print() << "    Number of physical particles per cell: " << num_par_per_cell << "\n"
                << "    Number of super droplets per cell: " << num_sd_per_cell << "\n";
        if (num_par_per_cell == 0) { return; }
        if (num_sd_per_cell == 0) { return; }
    }

    const int num_sp  = m_num_species;
    const int num_ae = m_num_aerosols;
    const Real rho_w = m_species_mat[m_idx_w]->m_density;
    const int idx_w = m_idx_w;
    AMREX_ALWAYS_ASSERT(idx_w >= 0);

    const auto sampled_multiplicity = a_init.sampledMultiplicity();

    iMultiFab num_superdroplets( ParticleBoxArray(a_lev),
                                 ParticleDistributionMap(a_lev),
                                 1, 0 );

    if (a_init.m_type == SDInitShape::uniform) {
        if (perbox) {
            Print() << "    Scattering super-droplets in box with volume " << a_init.volume() << ".\n";
            setNumSDPerBox( a_lev, num_superdroplets,
                            inject_count,
                            a_init.m_particle_domain,
                            a_init.injectSeed(),
                            a_height_ptr );
        } else {
            Print() << "    Adding particles in box with volume " << a_init.volume() << ".\n";
            setNumSDBoxDistribution( a_lev, num_superdroplets,
                                     num_sd_per_cell,
                                     a_height_ptr,
                                     a_init.m_particle_domain,
                                     subgrid );
        }
    } else if (a_init.m_type == SDInitShape::bubble) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( !perbox,
            "Per-box injection is not supported for bubble shape; set particles_per_cell." );
        Print() << "    Adding particles in bubble with volume: " << a_init.volume() << ".\n";
        setNumSDBubbleDistribution( a_lev, num_superdroplets,
                                    num_sd_per_cell,
                                    a_height_ptr,
                                    a_init.m_particle_domain,
                                    subgrid );
    } else if (a_init.m_type == SDInitShape::null) {
        num_superdroplets.setVal(0);
    } else {
        amrex::Print() << "Error: " << getEnumNameString(a_init.m_type)
                        << " is not a valid initialization for "
                        << m_name << " particle species.\n";
        amrex::Error("See error message!");
    }

    iMultiFab exposed;
    if (composite_init && a_lev < finestLevel()) {
        exposed = makeFineMask(num_superdroplets, ParticleBoxArray(a_lev+1),
                               GetParGDB()->refRatio(a_lev), 1, 0);
        for (MFIter mfi(num_superdroplets); mfi.isValid(); ++mfi) {
            const auto count = num_superdroplets.array(mfi);
            const auto mask = exposed.const_array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (!mask(i,j,k)) { count(i,j,k) = 0; }
            });
        }
    }

    if (per_cell_counts) {
        for (MFIter mfi(num_superdroplets); mfi.isValid(); ++mfi) {
            const auto count = num_superdroplets.array(mfi);
            const SDCellVolume cell_vol{dx, (*a_height_ptr)[mfi].const_array(), true};
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (count(i,j,k) > 0) { count(i,j,k) = cell_counts.numSDCell(cell_vol(i,j,k)); }
            });
        }
        // no cell receives a super-droplet (e.g. an injection step that adds no whole particle)
        if (num_superdroplets.sum(0) == 0) { return; }
    }

    iMultiFab offsets( ParticleBoxArray(a_lev),
                       ParticleDistributionMap(a_lev),
                       1, 0 );
    offsets.setVal(0);

    for(MFIter mfi = MakeMFIter(a_lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();

        int np = 0;
        {
            int ncell = static_cast<int>(num_superdroplets[mfi].numPts());
            const int* in = num_superdroplets[mfi].dataPtr();
            int* out = offsets[mfi].dataPtr();
            np = Scan::PrefixSum<int>( ncell,
                                       [=] AMREX_GPU_DEVICE (int i) -> int { return in[i]; },
                                       [=] AMREX_GPU_DEVICE (int i, int const &x) { out[i] = x; },
                                       Scan::Type::exclusive,
                                       Scan::retSum );
        }
        auto offset_arr = offsets[mfi].array();

        auto& particle_tile = DefineAndReturnParticleTile(a_lev, mfi);

        auto my_proc = ParallelDescriptor::MyProc();
        auto nprocs = ParallelDescriptor::NProcs();
        Long pid;
        {
            pid = ParticleType::NextID();
            ParticleType::NextID(pid+np);
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( static_cast<Long>(pid + np) < LastParticleID,
                                          "Error: overflow on particle id numbers!" );

        auto size_old = static_cast<Long>(particle_tile.size());
        auto size_new = size_old + np;
        particle_tile.resize(size_new);
        auto* aos = &particle_tile.GetArrayOfStructs()[0];
        auto& soa = particle_tile.GetStructOfArrays();

        /* SoA attributes */
        auto* vx_ptr = soa.GetRealData(SuperDropletsRealIdx::vx).data() + size_old;
        auto* vy_ptr = soa.GetRealData(SuperDropletsRealIdx::vy).data() + size_old;
        auto* vz_ptr = soa.GetRealData(SuperDropletsRealIdx::vz).data() + size_old;
        auto* mass_ptr = soa.GetRealData(SuperDropletsRealIdx::mass).data() + size_old;

        /* Runtime-added SoA attributes */
        int rt_off_i = SuperDropletsIntIdx::ncomps;
        auto* active_ptr = soa.GetIntData(rt_off_i+SuperDropletsIntIdxSoA_RT::active).data() + size_old;
        int rt_off_r = SuperDropletsRealIdx::ncomps;
        auto* radius_ptr = soa.GetRealData(rt_off_r+SuperDropletsRealIdxSoA_RT::radius).data() + size_old;
        auto* mult_ptr = soa.GetRealData(rt_off_r+SuperDropletsRealIdxSoA_RT::multiplicity).data() + size_old;
        auto* vterm_ptr = soa.GetRealData(rt_off_r+SuperDropletsRealIdxSoA_RT::term_vel).data() + size_old;
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        auto* condt_ptr = soa.GetRealData(rt_off_r+SuperDropletsRealIdxSoA_RT::cond_tendency).data() + size_old;
#endif
        auto* uid_ptr = soa.GetRealData(rt_off_r+SuperDropletsRealIdxSoA_RT::uid).data() + size_old;

        SDSpeciesMassArr sp_mass_ptrs;
        for (int i = 0; i < num_sp; i++) {
            sp_mass_ptrs[i] = soa.GetRealData(idx_s(i,num_ae,num_sp)).data() + size_old;
        }

        SDAerosolMassArr ae_mass_ptrs;
        for (int i = 0; i < num_ae; i++) {
            ae_mass_ptrs[i] = soa.GetRealData(idx_a(i,num_ae,num_sp)).data() + size_old;
        }

        // GPU-direct distribution generation: no large host memory allocation needed
        // Create distribution parameters on host (small, O(num_species + num_aerosols))
        Vector<SDDistributionParams> sp_params_h(num_sp);
        for (int i = 0; i < num_sp; i++) {
            sp_params_h[i] = a_init.getSpeciesDistParams(i, m_species_mat[i]->m_density, cell_volume, sampled_multiplicity);
        }
        Vector<SDDistributionParams> ae_params_h(num_ae);
        for (int i = 0; i < num_ae; i++) {
            ae_params_h[i] = a_init.getAerosolDistParams(i, m_aerosol_mat[i]->m_density, cell_volume, sampled_multiplicity);
        }

        // Copy parameters to device
        Gpu::DeviceVector<SDDistributionParams> sp_params_d(num_sp);
        Gpu::DeviceVector<SDDistributionParams> ae_params_d(num_ae);
        Gpu::copy(Gpu::hostToDevice, sp_params_h.begin(), sp_params_h.end(), sp_params_d.begin());
        Gpu::copy(Gpu::hostToDevice, ae_params_h.begin(), ae_params_h.end(), ae_params_d.begin());

        // Only need temporary storage for raw multiplicity values (for sampled case)
        Gpu::DeviceVector<Real> multiplicity_d(np, zero);

        auto* sp_params_ptr = sp_params_d.data();
        auto* ae_params_ptr = ae_params_d.data();
        auto* mult_ptr_init = multiplicity_d.data();

        // Generate species and aerosol distributions directly into particle SoA
        ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int n, const RandomEngine& engine) noexcept
        {
            // Sample species masses directly into particle storage
            for (int s = 0; s < num_sp; s++) {
                Real sp_mult = zero;
                sp_mass_ptrs[s][n] = SD_sample_mass_gpu(sp_params_ptr[s], engine, sp_mult);
            }

            // Sample aerosol masses directly into particle storage
            Real total_mult = zero;
            for (int a = 0; a < num_ae; a++) {
                Real ae_mult = zero;
                ae_mass_ptrs[a][n] = SD_sample_mass_gpu(ae_params_ptr[a], engine, ae_mult);
                if (sampled_multiplicity) {
                    total_mult += ae_mult;
                }
            }
            if (sampled_multiplicity) {
                mult_ptr_init[n] = total_mult;
            }
        });

        Gpu::synchronize();

        if (!m_device_props_initialized) { initializeDeviceProperties(); }
        const ParticleReal* sp_rho_arr = m_sp_density.data();
        const int* sp_sol_arr = m_sp_solubility.data();
        const ParticleReal* ae_rho_arr = m_ae_density.data();
        const int* ae_sol_arr = m_ae_solubility.data();
        const int* sp_INP_arr = m_sp_is_INP.data();
        const int* ae_INP_arr = m_ae_is_INP.data();

        auto mult_arr = multiplicity_d.data();

        auto num_superdroplets_arr = num_superdroplets[mfi].array();
        auto random_place = m_place_randomly_in_cells;

        auto pdomain = a_init.m_particle_domain;

        SDCellVolume cell_vol{dx, Array4<Real const>{}, variable_dz};
        if (variable_dz) { cell_vol.z_nd = (*a_height_ptr)[mfi].const_array(); }
        const auto z_nd = cell_vol.z_nd;

        ParallelForRNG(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k, const RandomEngine& rnd_engine) noexcept
        {
            int num_sd_this_cell = num_superdroplets_arr(i,j,k);
            const bool count_cell = per_cell_counts && (num_sd_this_cell > 0);
            const Real par_this_cell = count_cell ? cell_counts.numParticles(cell_vol(i,j,k)) : num_par_per_cell;
            const int  sd_this_cell  = count_cell ? num_sd_this_cell : num_sd_per_cell;
            Real num_to_add = par_this_cell;
            // per-cell counts are split evenly, so every super-droplet gets a whole particle
            const Real n_par_per_supdrop = count_cell ? std::floor(par_this_cell/sd_this_cell)
                                                      : std::ceil(par_this_cell/sd_this_cell);
            const Real n_par_extra = count_cell ? par_this_cell - n_par_per_supdrop*sd_this_cell : zero;

            Real mult_scale = one;
            if (sampled_multiplicity) {
                Real mult_sum = zero;
                int start = offset_arr(i,j,k);
                for (int n = start; n < start+num_sd_this_cell; n++) { mult_sum += mult_arr[n]; }
                mult_scale = par_this_cell / mult_sum;
            }

            int start = offset_arr(i,j,k);
            for (int n = start; n < start+num_sd_this_cell; n++) {

                auto& p = aos[n+size_old];
                p.id()  = pid + n;
                p.cpu() = my_proc;

                if (subgrid) {
                    if (itype == SDInitShape::uniform) {
                        if (random_place) {
                            p.pos(0) = pdomain.lo(0) + Random(rnd_engine) * (pdomain.hi(0) - pdomain.lo(0));
                            p.pos(1) = pdomain.lo(1) + Random(rnd_engine) * (pdomain.hi(1) - pdomain.lo(1));
                            p.pos(2) = pdomain.lo(2) + Random(rnd_engine) * (pdomain.hi(2) - pdomain.lo(2));
                        } else {
                            p.pos(0) = pdomain.lo(0) + myhalf * (pdomain.hi(0) - pdomain.lo(0));
                            p.pos(1) = pdomain.lo(1) + myhalf * (pdomain.hi(1) - pdomain.lo(1));
                            p.pos(2) = pdomain.lo(2) + myhalf * (pdomain.hi(2) - pdomain.lo(2));
                        }
                    } else if (itype == SDInitShape::bubble) {
                        if (random_place) {
                            ParticleReal u[AMREX_SPACEDIM];
                            ParticleReal L = 0.0;
                            do {
                                u[0] = 2.0*Random(rnd_engine) - 1.0;
                                u[1] = 2.0*Random(rnd_engine) - 1.0;
                                u[2] = 2.0*Random(rnd_engine) - 1.0;
                                L = std::sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
                            } while (L <= 1.0e-12);

                            const ParticleReal r = Random(rnd_engine);
                            const ParticleReal xs[AMREX_SPACEDIM] = {r*u[0]/L, r*u[1]/L, r*u[2]/L};
                            p.pos(0) = pdomain.lo(0) + pdomain.hi(0) * xs[0];
                            p.pos(1) = pdomain.lo(1) + pdomain.hi(1) * xs[1];
                            p.pos(2) = pdomain.lo(2) + pdomain.hi(2) * xs[2];
                        } else {
                            p.pos(0) = pdomain.lo(0);
                            p.pos(1) = pdomain.lo(1);
                            p.pos(2) = pdomain.lo(2);
                        }
                    }
                    if (variable_dz) {
                        // the domain is given in physical height; pos(2) is zeta
                        p.pos(2) = static_cast<ParticleReal>(ERF::ParticlePos::zeta_from_z(
                            p.pos(0), p.pos(1), p.pos(2), plo, dxi, z_nd, k_max));
                    }
                } else {
                    if (random_place) {
                        p.pos(0) = plo[0] + (i + Random(rnd_engine))*dx[0];
                        p.pos(1) = plo[1] + (j + Random(rnd_engine))*dx[1];
                        p.pos(2) = plo[2] + (k + Random(rnd_engine))*dx[2];
                    } else {
                        p.pos(0) = plo[0] + (i + myhalf)*dx[0];
                        p.pos(1) = plo[1] + (j + myhalf)*dx[1];
                        p.pos(2) = plo[2] + (k + myhalf)*dx[2];
                    }
                }

                active_ptr[n] = 1;
                vx_ptr[n] = vy_ptr[n] = vz_ptr[n] = zero;

                if (perbox) {
                    // Each scattered super-droplet carries the fixed per-box multiplicity.
                    mult_ptr[n] = inject_mult;
                } else {
                    Real mult_this_sd = 0;
                    if (sampled_multiplicity) {
                        mult_this_sd = std::ceil(mult_arr[n] * mult_scale);
                    } else {
                        mult_this_sd = n_par_per_supdrop + ((n - start < n_par_extra) ? one : zero);
                    }
                    if (mult_this_sd < num_to_add) {
                        mult_ptr[n] = mult_this_sd;
                    } else {
                        mult_ptr[n] = num_to_add;
                    }
                    num_to_add -= static_cast<Real>(mult_ptr[n]);
                    if (mult_ptr[n] == 0) { mult_ptr[n] = 1; }
                }

                // Species and aerosol masses already sampled directly into particle SoA

                radius_ptr[n] = SD_effective_radius( n, idx_w,
                                                     rho_w,
                                                     num_sp, num_ae,
                                                     sp_sol_arr, ae_sol_arr,
                                                     sp_mass_ptrs, ae_mass_ptrs,
                                                     sp_rho_arr, ae_rho_arr );
                mass_ptr[n] = SD_total_mass( n, num_sp, num_ae, sp_mass_ptrs, ae_mass_ptrs);

                vterm_ptr[n] = zero;
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
                condt_ptr[n] = zero;
#endif
                uid_ptr[n] = ParticleReal((pid+n-1)*nprocs + my_proc + 1);
            }
        });

        Gpu::synchronize();

        /* Initialize ice attributes
         *
         * Freezing temperature (Tfz) is assigned using the INAS (Ice Nucleation
         * Active Site) density parameterization of Niemand et al. (2012).
         *
         * The following code is adapted from SCALE-SDM:
         * https://github.com/Shima-Lab/SCALE-SDM_mixed-phase_Shima2019
         * Copyright (c) 2012-2015, Team SCALE
         * All rights reserved.
         * BSD 2-Clause License
         */
        if (m_idx_i >= 0) {
            auto* Tfz_ptr = soa.GetRealData(idx_ice_Tfz(num_ae,num_sp)).data() + size_old;
            auto* a_ptr = soa.GetRealData(idx_ice_a(num_ae,num_sp)).data() + size_old;
            auto* c_ptr = soa.GetRealData(idx_ice_c(num_ae,num_sp)).data() + size_old;
            auto* mrime_ptr = soa.GetRealData(idx_ice_mrime(num_ae,num_sp)).data() + size_old;
            auto* nmono_ptr = soa.GetRealData(idx_ice_nmono(num_ae,num_sp)).data() + size_old;
            const auto idx_i = m_idx_i;
            // Apparent density of the seeded ice particle: a prescribed bulk
            // density (e.g. a low-density snowflake) when set, else solid ice.
            const Real rho_i = (a_init.m_ice_app_density > Real(0.0)
                                ? a_init.m_ice_app_density
                                : m_species_mat[m_idx_i]->m_density);

            // INAS parameterization for sampling freezing temperature
            INAS_Niemand2012 inas_params;

            ParallelForRNG(np, [=] AMREX_GPU_DEVICE(int i, const RandomEngine& rnd_engine)
            {
                // Compute INP surface area for this particle
                auto A_INP = SD_INP_surface_area(i, num_sp, num_ae,
                                                  sp_INP_arr, ae_INP_arr,
                                                  sp_mass_ptrs, ae_mass_ptrs,
                                                  sp_rho_arr, ae_rho_arr);

                // Sample freezing temperature using INAS parameterization
                // Skip random values to decorrelate from earlier ParallelForRNG calls
                // (mass sampling also uses GPU random engines in separate kernels)
                (void)Random(rnd_engine);  // skip first value to decorrelate
                (void)Random(rnd_engine);  // skip second value
                auto u = Random(rnd_engine);  // use third value for Tfz sampling
                Tfz_ptr[i] = inas_params.sample_Tfz(u, A_INP);

                // Initialize ice shape attributes
                auto mass = sp_mass_ptrs[idx_i][i];
                a_ptr[i] = c_ptr[i] = std::cbrt(mass/(four_thirds_pi*rho_i));
                mrime_ptr[i] = 0.0;
                nmono_ptr[i] = (mass > 0 ? 1.0 : 0.0);
            });
            Gpu::synchronize();
        } else {
            // Without ice these attributes are unused, but coalescence still merges T_fz
            auto* Tfz_ptr = soa.GetRealData(idx_ice_Tfz(num_ae,num_sp)).data() + size_old;
            auto* a_ptr = soa.GetRealData(idx_ice_a(num_ae,num_sp)).data() + size_old;
            auto* c_ptr = soa.GetRealData(idx_ice_c(num_ae,num_sp)).data() + size_old;
            auto* mrime_ptr = soa.GetRealData(idx_ice_mrime(num_ae,num_sp)).data() + size_old;
            auto* nmono_ptr = soa.GetRealData(idx_ice_nmono(num_ae,num_sp)).data() + size_old;
            ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
            {
                Tfz_ptr[i] = a_ptr[i] = c_ptr[i] = mrime_ptr[i] = nmono_ptr[i] = ParticleReal(0);
            });
            Gpu::synchronize();
        }
    }

    return;
}

#endif

#include <random>
#include "ERF_SuperDropletPC.H"

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
void SuperDropletPC::setNumSDBoxDistribution (iMultiFab& a_num_sd,
                                              const int a_n_per_cell,
                                              const MFPtr& a_height_ptr,
                                              const RealBox& a_box,
                                              const bool a_subgrid)
{
    BL_PROFILE("SuperDropletPC::setNumSDBoxDistribution()");
    a_num_sd.setVal(0);

    const auto dx = Geom(m_lev).CellSizeArray();
    const auto plo = Geom(m_lev).ProbLoArray();

    for(MFIter mfi = MakeMFIter(m_lev); mfi.isValid(); ++mfi) {
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
                    if (a_box.contains(RealVect(x,y,z))) { flag = true; }
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
                    if (a_box.contains(RealVect(x,y,z))) { flag = true; }
                }
                if (flag) { num_superdroplets_arr(i,j,k) = a_n_per_cell; }
            });
        }
    }

    return;
}

/*! Sets the initial number of the super-droplets per cell as a bubble with a uniform distribution */
void SuperDropletPC::setNumSDBubbleDistribution ( iMultiFab& a_num_sd, /*!< integer Multifab with number of superdroplets in each grid cell */
                                                  const int a_n_per_cell, /*!< number of superdroplets per cell */
                                                  const MFPtr& a_height_ptr, /*!< terrain */
                                                  const RealBox& a_bubble, /*!< bubble within which to initialize particles */
                                                  const bool a_subgrid /*!< Is a_box smaller than grid cell */)
{
    BL_PROFILE("SuperDropletPC::setNumSDBubbleDistribution()");
    a_num_sd.setVal(0);

    const auto dx = Geom(m_lev).CellSizeArray();
    const auto plo = Geom(m_lev).ProbLoArray();

    for(MFIter mfi = MakeMFIter(m_lev); mfi.isValid(); ++mfi) {
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
                    if (x_r[0] > 0) rad += std::pow((x - x_c[0])/x_r[0], 2);
                    if (x_r[1] > 0) rad += std::pow((y - x_c[1])/x_r[1], 2);
                    if (x_r[2] > 0) rad += std::pow((z - x_c[2])/x_r[2], 2);
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
                    if (x_r[0] > 0) rad += std::pow((x - x_c[0])/x_r[0], 2);
                    if (x_r[1] > 0) rad += std::pow((y - x_c[1])/x_r[1], 2);
                    if (x_r[2] > 0) rad += std::pow((z - x_c[2])/x_r[2], 2);
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
void SuperDropletPC::addParticles ( const MFPtr& a_height_ptr, /*!< terrain */
                                    const SDInitProperties& a_init /*!< initialization parameters */ )
{
    BL_PROFILE("SuperDropletPC::addParticles");

    const auto dx_h = Geom(m_lev).CellSize();
    const Real cell_volume = dx_h[0]*dx_h[1]*dx_h[2];
    const auto dx = Geom(m_lev).CellSizeArray();
    const auto plo = Geom(m_lev).ProbLoArray();
    const auto init_volume = a_init.volume();
    const bool subgrid = (init_volume < cell_volume);
    const auto itype = a_init.m_type;

    // number of super-droplets per cell
    int num_sd_per_cell = a_init.numSDPerCell(cell_volume);

    // number of physical particles per cell
    Real num_par_per_cell = a_init.numParticlesPerCell(cell_volume);

    Print() << "    Number of physical particles per cell: " << num_par_per_cell << "\n"
            << "    Number of super droplets per cell: " << num_sd_per_cell << "\n";

    if (num_par_per_cell == 0) { return; }
    if (num_sd_per_cell == 0) { return; }

    const int num_sp  = m_num_species;
    const int num_ae = m_num_aerosols;
    const Real rho_w = m_species_mat[m_idx_w]->m_density;
    const int idx_w = m_idx_w;
    AMREX_ALWAYS_ASSERT(idx_w >= 0);

    const auto sampled_multiplicity = a_init.sampledMultiplicity();

    iMultiFab num_superdroplets( ParticleBoxArray(m_lev),
                                 ParticleDistributionMap(m_lev),
                                 1, 0 );

    if (a_init.m_type == SDInitShape::uniform) {
        Print() << "    Adding particles in box with volume " << a_init.volume() << ".\n";
        setNumSDBoxDistribution( num_superdroplets,
                                 num_sd_per_cell,
                                 a_height_ptr,
                                 a_init.m_particle_domain,
                                 subgrid );
    } else if (a_init.m_type == SDInitShape::bubble) {
        Print() << "    Adding particles in bubble with volume: " << a_init.volume() << ".\n";
        setNumSDBubbleDistribution( num_superdroplets,
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

    iMultiFab offsets( ParticleBoxArray(m_lev),
                       ParticleDistributionMap(m_lev),
                       1, 0 );
    offsets.setVal(0);

    for(MFIter mfi = MakeMFIter(m_lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();

        int np = 0;
        {
            int ncell = num_superdroplets[mfi].numPts();
            const int* in = num_superdroplets[mfi].dataPtr();
            int* out = offsets[mfi].dataPtr();
            np = Scan::PrefixSum<int>( ncell,
                                       [=] AMREX_GPU_DEVICE (int i) -> int { return in[i]; },
                                       [=] AMREX_GPU_DEVICE (int i, int const &x) { out[i] = x; },
                                       Scan::Type::exclusive,
                                       Scan::retSum );
        }
        auto offset_arr = offsets[mfi].array();

        auto& particle_tile = DefineAndReturnParticleTile(m_lev, mfi);

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
        auto* vx_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::vx).data() + size_old;
        auto* vy_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::vy).data() + size_old;
        auto* vz_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::vz).data() + size_old;
        auto* mass_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::mass).data() + size_old;

        /* Runtime-added SoA attributes */
        int rt_off_i = SuperDropletsIntIdxSoA::ncomps;
        auto* active_ptr = soa.GetIntData(rt_off_i+SuperDropletsIntIdxSoA_RT::active).data() + size_old;
        int rt_off_r = SuperDropletsRealIdxSoA::ncomps;
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

        auto mult_arr = multiplicity_d.data();

        auto num_superdroplets_arr = num_superdroplets[mfi].array();
        auto random_place = m_place_randomly_in_cells;

        auto pdomain = a_init.m_particle_domain;

        ParallelForRNG(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k, const RandomEngine& rnd_engine) noexcept
        {
            int num_sd_this_cell = num_superdroplets_arr(i,j,k);
            Real num_to_add = num_par_per_cell;
            Real n_par_per_supdrop = std::ceil(num_par_per_cell/num_sd_per_cell);

            Real mult_scale = one;
            if (sampled_multiplicity) {
                Real mult_sum = zero;
                int start = offset_arr(i,j,k);
                for (int n = start; n < start+num_sd_this_cell; n++) { mult_sum += mult_arr[n]; }
                mult_scale = num_par_per_cell / mult_sum;
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

                p.idata(SuperDropletsIntIdxAoS::k) = k;
                active_ptr[n] = 1;
                vx_ptr[n] = vy_ptr[n] = vz_ptr[n] = zero;

                Real mult_this_sd = 0;
                if (sampled_multiplicity) {
                    mult_this_sd = std::ceil(mult_arr[n] * mult_scale);
                } else {
                    mult_this_sd = n_par_per_supdrop;
                }
                if (mult_this_sd < num_to_add) {
                    mult_ptr[n] = mult_this_sd;
                } else {
                    mult_ptr[n] = num_to_add;
                }
                num_to_add -= mult_ptr[n];
                if (mult_ptr[n] == 0) { mult_ptr[n] = 1; }

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

        const auto height_arr = (*a_height_ptr)[mfi].array();
        ParallelFor(tile_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int num_sd_this_cell = num_superdroplets_arr(i,j,k);
            int start = offset_arr(i,j,k);
            for (int n = start; n < start+num_sd_this_cell; n++) {
                auto& p = aos[n+size_old];
                Real x = p.pos(0);
                Real y = p.pos(1);
                Real z = p.pos(2);
                Real r[3] = { (x-plo[0])/dx[0] - i,
                              (y-plo[1])/dx[1] - j,
                              (z-plo[2])/dx[2] - k };

                Real sx[] = { one - r[0], r[0]};
                Real sy[] = { one - r[1], r[1]};

                Real height_at_pxy_lo = zero;
                Real height_at_pxy_hi = zero;
                for (int ii = 0; ii < 2; ++ii) {
                    for (int jj = 0; jj < 2; ++jj) {
                        height_at_pxy_lo += sx[ii] * sy[jj]
                                            * height_arr(i+ii,j+jj,k);
                        height_at_pxy_hi += sx[ii] * sy[jj]
                                            * height_arr(i+ii,j+jj,k+1);
                    }
                }

                p.pos(2) = height_at_pxy_lo
                           + r[2] * (height_at_pxy_hi - height_at_pxy_lo);
           }
        });
        Gpu::synchronize();
    }

    return;
}

#endif

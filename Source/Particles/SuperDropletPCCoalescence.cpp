#include "SuperDropletPC.H"
#include "SuperDropletPCCoalescence.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! \brief Swap two numbers; NOTE: a and b can't be the same memory location */
template <typename dtype>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static void swap ( dtype& a, /*!< first number */
                   dtype& b  /*!< second number */ )
{
    a = a + b;
    b = a - b;
    a = a - b;
}

/*! \brief Random shuffle of an array */
template <typename dtype>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static void random_shuffle ( dtype* const a_d_ptr, /*!< Pointer to data array */
                             const int  a_n, /*!< size of data array */
                             const RandomEngine& a_rnd_eng /*!< random engine */ )
{
    int num_passes = Random_int(100, a_rnd_eng);
    for (int ipass = 0; ipass < num_passes; ipass++) {
        for (int i = 0; i < a_n; i++) {
            int j = Random_int(a_n, a_rnd_eng);
            int k = Random_int(a_n, a_rnd_eng);
            if (j == k) { continue; }
            swap<dtype>( a_d_ptr[j], a_d_ptr[k] );
        }
    }
}

/*! \brief Binary coalescence between two superdroplets */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static bool binary_coalescence (const int a_i, /*!< index of first particle */
                                const int a_j, /*!< index of second particle */
                                const RandomEngine& a_rnd_eng, /*!< random engine */
                                const Real a_p, /*!< probability */
                                ParticleReal* const a_mass, /*!< mass */
                                ParticleReal* const a_radius, /*!< radius */
                                ParticleReal* const a_mult, /*!< multiplicity */
                                ParticleReal* const a_sd_mass, /*!<superdroplet mass*/
                                const int a_n_aerosols, /*!< number of aerosols */
                                const SDAerosolMassArr& a_aero_masses /*!< aerosol masses*/)
{
    ParticleReal p_int = std::floor(a_p);
    ParticleReal gamma = p_int;
    if (Random(a_rnd_eng) < (a_p-p_int)) { gamma += 1; }

    if (gamma != 0) {
        int i = a_i, j = a_j;
        if (a_mult[a_i] < a_mult[a_j]) { i = a_j; j = a_i; }

        ParticleReal gamma_t = std::min( gamma,
                                         std::floor(a_mult[i]/a_mult[j]) );

        if ( (a_mult[i]-gamma_t*a_mult[j]) > 0 ) {

            a_mult[i] -= gamma_t*a_mult[j];

            ParticleReal r3 = gamma_t*a_radius[i]*a_radius[i]*a_radius[i]
                                    + a_radius[j]*a_radius[j]*a_radius[j];
            a_radius[j] = std::cbrt(r3);

            a_mass[j] += gamma_t*a_mass[i];
            for (int n = 0; n < a_n_aerosols; n++) {
                a_aero_masses[n][j] += gamma_t*a_aero_masses[n][i];
            }

            a_sd_mass[i] = a_mult[i] * a_mass[i];
            a_sd_mass[j] = a_mult[j] * a_mass[j];

        } else if ( (a_mult[i]-gamma_t*a_mult[j]) == 0 ) {

            ParticleReal dm = std::floor(a_mult[j]/2);
            a_mult[i] = dm;
            a_mult[j] -= dm;

            ParticleReal r3 = gamma_t*a_radius[i]*a_radius[i]*a_radius[i]
                                    + a_radius[j]*a_radius[j]*a_radius[j];
            a_radius[i] = a_radius[j] = std::cbrt(r3);

            a_mass[j] += gamma_t*a_mass[i];
            a_mass[i] = a_mass[j];
            for (int n = 0; n < a_n_aerosols; n++) {
                a_aero_masses[n][j] += gamma_t*a_aero_masses[n][i];
                a_aero_masses[n][i] = a_aero_masses[n][j];
            }

            a_sd_mass[i] = a_mult[i] * a_mass[i];
            a_sd_mass[j] = a_mult[j] * a_mass[j];

        }

        return true;

    } else {

        return false;

    }
}

/*! Compute the coalescence of superdroplets in each time step */
void SuperDropletPC::Coalescence( int   a_lev,
                                  Real  a_dt,
                                  const MultiFab& a_pressure,
                                  const MultiFab& a_temperature )
{
    BL_PROFILE("SuperDropletPC::Coalescence()");
    AMREX_ASSERT( a_lev == m_lev );

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const auto domain = geom.Domain();

    const int num_aerosols = m_num_aerosols;
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    const ParticleReal inv_bin_size
        = 1.0 / (  static_cast<ParticleReal>(m_coalescence_bin_size[0])
                 * static_cast<ParticleReal>(m_coalescence_bin_size[1])
                 * static_cast<ParticleReal>(m_coalescence_bin_size[2]) );
    const ParticleReal inv_bin_volume = inv_cell_volume*inv_bin_size;

    long num_collisions = 0;
    const auto& gvec = a_temperature.nGrowVect();

    auto kernel_choice = m_coalescence_kernel;
    auto include_brownian_coalescence = m_include_brownian_coalescence;

// Do NOT add OpenMP here; building DenseBins is not thread-safe.
    for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {

        auto& ptile = ParticlesAt(a_lev, pti);
        auto& aos = ptile.GetArrayOfStructs();
        auto& soa = ptile.GetStructOfArrays();
        const size_t np = aos.numParticles();
        auto pstruct_ptr = aos().dataPtr();

        /* SoA attributes */
        auto* mass_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::mass).data();
        Array<ParticleReal*,AMREX_SPACEDIM> v_ptr;
        v_ptr[0] = soa.GetRealData(SuperDropletsRealIdxSoA::vx).data();
        v_ptr[1] = soa.GetRealData(SuperDropletsRealIdxSoA::vy).data();
        v_ptr[2] = soa.GetRealData(SuperDropletsRealIdxSoA::vz).data();

        /* Runtime-added SoA attributes */
        int rt_offset = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::radius).data();
        auto* mult_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::multiplicity).data();
        auto* supdrop_mass_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::sd_mass).data();
        auto* t_coal_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::t_coalescence).data();
        auto* vterm_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::term_vel).data();

        /* aerosol masses */
        SDAerosolMassArr aerosol_mass_ptrs;
        for (int i = 0; i < num_aerosols; i++) {
            aerosol_mass_ptrs[i] = soa.GetRealData(   rt_offset
                                                    + SuperDropletsRealIdxSoA_RT::ncomps
                                                    + i ).data();
        }

        int grid = pti.index();
        Box box = a_temperature[grid].box(); box.grow(-gvec);
        int ntiles = numTilesInBox(box, true, m_coalescence_bin_size);
        auto binner = GetParticleBin{plo, dxi, domain, m_coalescence_bin_size, box};
        DenseBins<ParticleType> bins;
        bins.build( np, pstruct_ptr, ntiles, binner);
        AMREX_ALWAYS_ASSERT(np == static_cast<size_t>(bins.numItems()));
        AMREX_ALWAYS_ASSERT(bins.numBins() >= 0);
        auto inds = bins.permutationPtr();
        auto offsets = bins.offsetsPtr();

        const auto& pressure_arr = a_pressure[grid].const_array();
        const auto& temperature_arr = a_temperature[grid].const_array();

        CollisionKernel<ParticleReal,AMREX_SPACEDIM> ckernel{};

        ParticleReal condensate_density = m_vapour_mat->density();
        Gpu::DeviceVector<ParticleReal> aero_density_d;
        {
            Vector<ParticleReal> aero_density_h;
            aero_density_h.resize(num_aerosols);
            aero_density_d.resize(num_aerosols);
            for (int ia = 0; ia < num_aerosols; ia++) {
                aero_density_h[ia] = m_aerosol_mat[ia]->density();
            }
            Gpu::copy(  Gpu::hostToDevice,
                        aero_density_h.begin(),
                        aero_density_h.end(),
                        aero_density_d.begin() );
        }
        auto aero_density = aero_density_d.data();

        Gpu::Buffer<amrex::Long> particle_collisions({0});
        Long* particle_collisions_ptr = particle_collisions.data();

        if (m_coalescence_alg == SupDropInit::alg_coalescence_cudsmc) {

            ParallelForRNG( bins.numBins(),
                            [=] AMREX_GPU_DEVICE (int i,
                                                  RandomEngine const& rnd_eng) noexcept
            {
                auto bin_start = offsets[i];
                auto bin_stop = offsets[i+1];
                auto np_bin = bin_stop - bin_start;
                if (np_bin > 1) {
                    random_shuffle<unsigned int>(inds+bin_start, np_bin, rnd_eng);
                }
            } );
            Gpu::synchronize();

            Gpu::DeviceVector<unsigned int> np_hbin(bins.numBins()+1);
            Gpu::DeviceVector<unsigned int> offsets_hbin(bins.numBins()+1);
            auto np_hbin_ptr = np_hbin.data();
            ParallelFor( bins.numBins(), [=] AMREX_GPU_DEVICE (int i) noexcept
            {
                auto bin_start = offsets[i];
                auto bin_stop = offsets[i+1];
                auto np_bin = bin_stop - bin_start;
                np_hbin_ptr[i] = np_bin/2;
            });
            Gpu::synchronize();
            Gpu::exclusive_scan(np_hbin.begin(), np_hbin.end(), offsets_hbin.begin());

            Gpu::DeviceVector<unsigned int> inds_hbin;
            inds_hbin.resize(np/2);
            auto inds_hbin_ptr = inds_hbin.data();
            auto offsets_hbin_ptr = offsets_hbin.data();
            ParallelFor( bins.numBins(), [=] AMREX_GPU_DEVICE (int i) noexcept
            {
                auto fbin_start = offsets[i];
                auto hbin_start = offsets_hbin_ptr[i];
                auto hbin_stop = offsets_hbin_ptr[i+1];
                for (int i = 0; i < (hbin_stop-hbin_start); i++) {
                    inds_hbin_ptr[i+hbin_start] = inds[i+fbin_start];
                }
            });
            Gpu::synchronize();

            ParallelForRNG( np/2, [=] AMREX_GPU_DEVICE (int ib,
                                                        RandomEngine const& rnd_eng) noexcept
            {
                auto pi = inds_hbin_ptr[ib];
                int i_bin = binner(pstruct_ptr[pi]);

                auto hbin_start = offsets_hbin_ptr[i_bin];
                auto fbin_start = offsets[i_bin];
                auto fbin_stop = offsets[i_bin+1];
                auto np_bin = fbin_stop - fbin_start;

                if (np_bin > 1) {

                    int p = ib - hbin_start;
                    auto pj = inds[fbin_stop-1-p];

                    if (pi == pj) { return; }
                    if (mult_ptr[pi] == 0) { return; }
                    if (mult_ptr[pj] == 0) { return; }

                    ParticleReal v_i[AMREX_SPACEDIM], v_j[AMREX_SPACEDIM];
                    for (int d = 0; d < AMREX_SPACEDIM; d++) {
                        v_i[d] = v_ptr[d][pi];
                        v_j[d] = v_ptr[d][pj];
                    }

                    ParticleReal k_val = 0.0;
                    if (kernel_choice == SDCoalescenceKernelType::golovin) {

                        k_val = ckernel.golovin(radius_ptr[pi],radius_ptr[pj],v_i,v_j);

                    } else {

                        for (int d = 0; d < AMREX_SPACEDIM; d++) { v_i[d] = v_j[d] = 0; }
                        v_i[AMREX_SPACEDIM-1] = vterm_ptr[pi];
                        v_j[AMREX_SPACEDIM-1] = vterm_ptr[pj];

                        if (kernel_choice == SDCoalescenceKernelType::sedimentation) {
                            k_val = ckernel.sedimentation(radius_ptr[pi],radius_ptr[pj],v_i,v_j);
                        } else if (kernel_choice == SDCoalescenceKernelType::Longs) {
                            k_val = ckernel.Longs(radius_ptr[pi],radius_ptr[pj],v_i,v_j);
                        } else if (kernel_choice == SDCoalescenceKernelType::Halls) {
                            k_val = ckernel.Halls(radius_ptr[pi],radius_ptr[pj],v_i,v_j);
                        }
                    }

                    if (include_brownian_coalescence) {

                        ParticleReal pressure = 0.0, temperature = 0.0;
                        {
                            ParticleType& par_1 = pstruct_ptr[pi];
                            auto iv = getParticleCell(par_1, plo, dxi, domain);
                            pressure = pressure_arr(iv[0],iv[1],iv[2],0);
                            temperature = temperature_arr(iv[0],iv[1],iv[2],0);
                        }

                        ParticleReal aero_mass_1 = 0.0,
                                     aero_mass_2 = 0.0,
                                     aero_vol_1 = 0.0,
                                     aero_vol_2 = 0.0;
                        {
                            for (int ia = 0; ia < num_aerosols; ia++) {
                                aero_mass_1 += aerosol_mass_ptrs[ia][pi];
                                aero_mass_2 += aerosol_mass_ptrs[ia][pj];
                                aero_vol_1 += aerosol_mass_ptrs[ia][pi]/aero_density[ia];
                                aero_vol_2 += aerosol_mass_ptrs[ia][pj]/aero_density[ia];
                            }
                        }

                        k_val += ckernel.Brownian_SeinfeldPandis( radius_ptr[pi],
                                                                  radius_ptr[pj],
                                                                  aero_mass_1,
                                                                  aero_mass_2,
                                                                  aero_vol_1,
                                                                  aero_vol_2,
                                                                  condensate_density,
                                                                  pressure,
                                                                  temperature );
                    }

                    auto prob_ij = k_val*inv_bin_volume;
                    auto prob_sd_ij = std::max(mult_ptr[pi],mult_ptr[pj])*prob_ij;

                    auto ns = static_cast<ParticleReal>(np_bin);
                    auto scaling_factor = 0.5*ns*(ns-1)/std::floor(0.5*ns);
                    auto scaled_prob = prob_sd_ij * scaling_factor;

                    auto t_coalescence = 1.0/scaled_prob;
                    t_coal_ptr[pi] = t_coal_ptr[pj] = t_coalescence;

                    auto flag = binary_coalescence( pi, pj,
                                                    rnd_eng,
                                                    (scaled_prob * a_dt),
                                                    mass_ptr,
                                                    radius_ptr,
                                                    mult_ptr,
                                                    supdrop_mass_ptr,
                                                    num_aerosols,
                                                    aerosol_mass_ptrs );

                    amrex::Gpu::Atomic::Add(particle_collisions_ptr, amrex::Long(flag));
                }

            } );
            Gpu::synchronize();

        } else if (m_coalescence_alg == SupDropInit::alg_coalescence_dsmc) {

            ParallelForRNG( bins.numBins(),
                            [=] AMREX_GPU_DEVICE (int i_bin,
                                                  RandomEngine const& rnd_eng) noexcept
            {
                auto bin_start = offsets[i_bin];
                auto bin_stop = offsets[i_bin+1];
                auto np_bin = bin_stop - bin_start;

                if (np_bin > 1) {

                    random_shuffle<unsigned int>(inds+bin_start, np_bin, rnd_eng);

                    for (unsigned int p = 0; p < np_bin/2; p++) {
                        auto pi = inds[bin_start+p];
                        auto pj = inds[bin_stop-1-p];

                        if (pi == pj) { continue; }
                        if (mult_ptr[pi] == 0) { return; }
                        if (mult_ptr[pj] == 0) { return; }

                        ParticleReal v_i[AMREX_SPACEDIM], v_j[AMREX_SPACEDIM];
                        for (int d = 0; d < AMREX_SPACEDIM; d++) {
                            v_i[d] = v_ptr[d][pi];
                            v_j[d] = v_ptr[d][pj];
                        }

                        ParticleReal k_val = 0.0;
                        if (kernel_choice == SDCoalescenceKernelType::golovin) {

                            k_val = ckernel.golovin(radius_ptr[pi],radius_ptr[pj],v_i,v_j);

                        } else {

                            for (int d = 0; d < AMREX_SPACEDIM; d++) { v_i[d] = v_j[d] = 0; }
                            v_i[AMREX_SPACEDIM-1] = vterm_ptr[pi];
                            v_j[AMREX_SPACEDIM-1] = vterm_ptr[pj];

                            if (kernel_choice == SDCoalescenceKernelType::sedimentation) {
                                k_val = ckernel.sedimentation(radius_ptr[pi],radius_ptr[pj],v_i,v_j);
                            } else if (kernel_choice == SDCoalescenceKernelType::Longs) {
                                k_val = ckernel.Longs(radius_ptr[pi],radius_ptr[pj],v_i,v_j);
                            } else if (kernel_choice == SDCoalescenceKernelType::Halls) {
                                k_val = ckernel.Halls(radius_ptr[pi],radius_ptr[pj],v_i,v_j);
                            }
                        }

                        if (include_brownian_coalescence) {

                            ParticleReal pressure = 0.0, temperature = 0.0;
                            {
                                ParticleType& par_1 = pstruct_ptr[pi];
                                auto iv = getParticleCell(par_1, plo, dxi, domain);
                                pressure = pressure_arr(iv[0],iv[1],iv[2],0);
                                temperature = temperature_arr(iv[0],iv[1],iv[2],0);
                            }

                            ParticleReal aero_mass_1 = 0.0,
                                         aero_mass_2 = 0.0,
                                         aero_vol_1 = 0.0,
                                         aero_vol_2 = 0.0;
                            {
                                for (int ia = 0; ia < num_aerosols; ia++) {
                                    aero_mass_1 += aerosol_mass_ptrs[ia][pi];
                                    aero_mass_2 += aerosol_mass_ptrs[ia][pj];
                                    aero_vol_1 += aerosol_mass_ptrs[ia][pi]/aero_density[ia];
                                    aero_vol_2 += aerosol_mass_ptrs[ia][pj]/aero_density[ia];
                                }
                            }

                            k_val += ckernel.Brownian_SeinfeldPandis( radius_ptr[pi],
                                                                      radius_ptr[pj],
                                                                      aero_mass_1,
                                                                      aero_mass_2,
                                                                      aero_vol_1,
                                                                      aero_vol_2,
                                                                      condensate_density,
                                                                      pressure,
                                                                      temperature );
                        }

                        auto prob_ij = k_val*inv_bin_volume;
                        auto prob_sd_ij = std::max(mult_ptr[pi],mult_ptr[pj])*prob_ij;

                        auto ns = static_cast<ParticleReal>(np_bin);
                        auto scaling_factor = 0.5*ns*(ns-1)/std::floor(0.5*ns);
                        auto scaled_prob = prob_sd_ij * scaling_factor;

                        auto t_coalescence = 1.0/scaled_prob;
                        t_coal_ptr[pi] = t_coal_ptr[pj] = t_coalescence;

                        auto flag = binary_coalescence( pi, pj,
                                                        rnd_eng,
                                                        (scaled_prob * a_dt),
                                                        mass_ptr,
                                                        radius_ptr,
                                                        mult_ptr,
                                                        supdrop_mass_ptr,
                                                        num_aerosols,
                                                        aerosol_mass_ptrs );

                        amrex::Gpu::Atomic::Add(particle_collisions_ptr, amrex::Long(flag));
                    }
                }

            } );
            Gpu::synchronize();

        }

        num_collisions = *(particle_collisions.copyToHost());
    }

    ParallelDescriptor::ReduceLongSum(  &num_collisions,
                                        1,
                                        ParallelDescriptor::IOProcessorNumber() );

    Print() << "SuperDropletPC(" << m_name << "): "
            << "number of collisions = " << num_collisions << "\n";
}

#endif


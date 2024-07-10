#include <sys/time.h>
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

/*! \brief Compute coalescence rate between two superdroplets */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static ParticleReal coalescence_rate ( const RandomEngine& a_rnd_eng, /*!< random engine */
                                       const Real a_p /*!< probability */ )
{
    ParticleReal p_int = std::floor(a_p);
    ParticleReal gamma = p_int;
    if (Random(a_rnd_eng) < (a_p-p_int)) { gamma += 1; }
    return gamma;
}

/*! \brief Binary coalescence between two superdroplets */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
static void coal_update_attribs(const int a_i, /*!< index of particle */
                                const int* const a_idx, /*!< index of coalescence partner */
                                const int* const a_prey, /*!< prey/predator */
                                const ParticleReal* const a_gamma, /*!< coalescence rate */
                                const ParticleReal* const a_rmndr, /*!< coalescence remainder*/
                                ParticleReal* const a_mass, /*!< mass */
                                ParticleReal* const a_radius, /*!< radius */
                                ParticleReal* const a_mult, /*!< multiplicity */
                                const int a_n_aerosols, /*!< number of aerosols */
                                const SDAerosolMassArr& a_aero_masses /*!< aerosol masses*/)
{
    int i = a_i;
    int j = a_idx[a_i];
    auto gamma = a_gamma[i];
    AMREX_ALWAYS_ASSERT(gamma == a_gamma[j]);
    AMREX_ALWAYS_ASSERT(a_rmndr[a_i] >= 0.0);

    if ( a_rmndr[a_i] > 0 ) {

        if (a_prey[i]) {
            a_mult[i] -= gamma*a_mult[j];
        } else {
            auto r3 = gamma*a_radius[j]*a_radius[j]*a_radius[j]
                            + a_radius[i]*a_radius[i]*a_radius[i];
            a_radius[i] = std::cbrt(r3);
            a_mass[i] += gamma*a_mass[j];
            for (int n = 0; n < a_n_aerosols; n++) {
                a_aero_masses[n][i] += gamma*a_aero_masses[n][j];
            }
        }

    } else if ( a_rmndr[a_i] == 0 ) {

        if (a_prey[i]) {
            ParticleReal dm = std::floor(a_mult[j]/2);
            a_mult[i] = dm;
            a_mult[j] -= dm;
            ParticleReal r3 = gamma*a_radius[i]*a_radius[i]*a_radius[i]
                                    + a_radius[j]*a_radius[j]*a_radius[j];
            a_radius[i] = a_radius[j] = std::cbrt(r3);
            a_mass[j] += gamma*a_mass[i];
            a_mass[i] = a_mass[j];
            for (int n = 0; n < a_n_aerosols; n++) {
                a_aero_masses[n][j] += gamma*a_aero_masses[n][i];
                a_aero_masses[n][i] = a_aero_masses[n][j];
            }
        }

    }

}

/*! Compute the coalescence of superdroplets in each time step */
void SuperDropletPC::Coalescence( int   a_lev,
                                  Real  a_dt,
                                  const MultiFab& a_pressure,
                                  const MultiFab& a_temperature )
{
    struct timeval coal_start, coal_end;
    gettimeofday(&coal_start, NULL);

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

        Gpu::DeviceVector<int> coal_partner_idx, flag_prey, num_particles_bin;
        Gpu::DeviceVector<Real> coal_rate, coal_rmndr;
        num_particles_bin.resize(np);
        coal_partner_idx.resize(np);
        flag_prey.resize(np);
        coal_rate.resize(np);
        coal_rmndr.resize(np);
        auto np_bin_ptr = num_particles_bin.data();
        auto partner_idx_ptr = coal_partner_idx.data();
        auto flag_prey_ptr = flag_prey.data();
        auto coal_rate_ptr = coal_rate.data();
        auto coal_rmndr_ptr = coal_rmndr.data();
        ParallelFor( np, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            np_bin_ptr[i] = 0;
            partner_idx_ptr[i] = -1;
            flag_prey_ptr[i] = -1;
            coal_rate_ptr[i] = 0.0;
            coal_rmndr_ptr[i] = 0.0;
        });
        Gpu::synchronize();

        ParallelForRNG( bins.numBins(),
                        [=] AMREX_GPU_DEVICE (int i_bin,
                                              RandomEngine const& rnd_eng) noexcept
        {
            auto bin_start = offsets[i_bin];
            auto bin_stop = offsets[i_bin+1];
            auto np_bin = bin_stop - bin_start;
            if (np_bin <= 1) { return; }

            random_shuffle<unsigned int>(inds+bin_start, np_bin, rnd_eng);
            for (unsigned int p = 0; p < np_bin/2; p++) {
                auto pi = inds[bin_start+p];
                auto pj = inds[bin_stop-1-p];

                if (pi == pj) { continue; }
                if (mult_ptr[pi] == 0) { continue; }
                if (mult_ptr[pj] == 0) { continue; }

                np_bin_ptr[pi] = np_bin_ptr[pj] = np_bin;
                partner_idx_ptr[pi] = pj;
                partner_idx_ptr[pj] = pi;

                int i = -1, j = -1;
                if (mult_ptr[pi] >= mult_ptr[pj]) { i = pi; j = pj; }
                else                              { i = pj; j = pi; }
                flag_prey_ptr[i] = 1;
                flag_prey_ptr[j] = 0;
            }
        } );
        Gpu::synchronize();

        ParallelForRNG( np,
                        [=] AMREX_GPU_DEVICE (int i,
                                              RandomEngine const& rnd_eng) noexcept
        {
            if (partner_idx_ptr[i] < 0) { return; }
            if (!flag_prey_ptr[i]) { return; }

            int pi = i;
            int pj = partner_idx_ptr[i];
            AMREX_ALWAYS_ASSERT(mult_ptr[pi] >= mult_ptr[pj]);

            ParticleReal k_val = 0.0;
            if (kernel_choice == SDCoalescenceKernelType::golovin) {

                k_val = ckernel.golovin(radius_ptr[pi],radius_ptr[pj]);

            } else {

                ParticleReal v_i[AMREX_SPACEDIM], v_j[AMREX_SPACEDIM];
                for (int d = 0; d < AMREX_SPACEDIM; d++) {
                    v_i[d] = v_ptr[d][pi];
                    v_j[d] = v_ptr[d][pj];
                }

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

            auto ns = static_cast<ParticleReal>(np_bin_ptr[i]);
            auto scaling_factor = 0.5*ns*(ns-1)/std::floor(0.5*ns);
            auto scaled_prob = prob_sd_ij * scaling_factor;

            auto t_coalescence = 1.0/scaled_prob;
            t_coal_ptr[i] = t_coalescence;

            auto gamma = coalescence_rate ( rnd_eng, (scaled_prob*a_dt) );
            if (gamma > 0) {
                amrex::Gpu::Atomic::Add(particle_collisions_ptr, amrex::Long(gamma));
                coal_rate_ptr[pi] = std::min(gamma,std::floor(mult_ptr[pi]/mult_ptr[pj]));
                coal_rate_ptr[pj] = coal_rate_ptr[pi];
                coal_rmndr_ptr[pi] = mult_ptr[pi] - coal_rate_ptr[pi]*mult_ptr[pj];
                coal_rmndr_ptr[pj] = coal_rmndr_ptr[pi];
            } else {
                partner_idx_ptr[pi] = -1;
                partner_idx_ptr[pj] = -1;
            }

        } );
        Gpu::synchronize();
        num_collisions = *(particle_collisions.copyToHost());

        ParallelFor( np, [=] AMREX_GPU_DEVICE (int i)
        {
            if (partner_idx_ptr[i] < 0) { return; }
            coal_update_attribs( i,
                                 partner_idx_ptr,
                                 flag_prey_ptr,
                                 coal_rate_ptr,
                                 coal_rmndr_ptr,
                                 mass_ptr,
                                 radius_ptr,
                                 mult_ptr,
                                 num_aerosols,
                                 aerosol_mass_ptrs );
        } );
        Gpu::synchronize();

        ParallelFor( np, [=] AMREX_GPU_DEVICE (int i)
                     { supdrop_mass_ptr[i] = mass_ptr[i] * mult_ptr[i]; } );

    }

    ParallelDescriptor::ReduceLongSum(  &num_collisions,
                                        1,
                                        ParallelDescriptor::IOProcessorNumber() );

    gettimeofday(&coal_end,NULL);
    long long coal_wtime;
    coal_wtime = (    (coal_end.tv_sec   * 1000000 + coal_end.tv_usec  )
                   -  (coal_start.tv_sec * 1000000 + coal_start.tv_usec) );
    Real coal_wtime_sec = (double) coal_wtime / 1000000.0;
    ParallelDescriptor::ReduceRealMax( &coal_wtime_sec,
                                       1,
                                       ParallelDescriptor::IOProcessorNumber() );

    Print() << "SuperDropletPC(" << m_name << "): "
            << "number of collisions = " << num_collisions
            << " (wall time = " << coal_wtime_sec << " seconds)"
            << "\n";
}

#endif


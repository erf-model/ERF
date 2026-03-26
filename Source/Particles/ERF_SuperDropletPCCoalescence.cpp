#ifndef _WIN32
#include <sys/time.h>
#endif
#include "ERF_SuperDropletPC.H"
#include "ERF_SuperDropletPCCoalescence.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;
using namespace SDPCDefn;

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
                                const int a_n_species, /*!< number of species */
                                const SDSpeciesMassArr& a_species_masses, /*!< species masses*/
                                const int a_n_aerosols, /*!< number of aerosols */
                                const SDAerosolMassArr& a_aero_masses /*!< aerosol masses*/)
{
    int i = a_i;
    int j = a_idx[a_i];
    auto gamma = a_gamma[i];
    AMREX_ALWAYS_ASSERT(gamma == a_gamma[j]);
    AMREX_ALWAYS_ASSERT(a_rmndr[a_i] >= zero);

    if ( a_rmndr[a_i] > 0 ) {

        if (a_prey[i]) {
            a_mult[i] -= gamma*a_mult[j];
        } else {
            auto r3 = gamma*a_radius[j]*a_radius[j]*a_radius[j]
                            + a_radius[i]*a_radius[i]*a_radius[i];
            a_radius[i] = std::cbrt(r3);
            a_mass[i] += gamma*a_mass[j];
            for (int n = 0; n < a_n_species; n++) {
                a_species_masses[n][i] += gamma*a_species_masses[n][j];
            }
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
            for (int n = 0; n < a_n_species; n++) {
                a_species_masses[n][j] += gamma*a_species_masses[n][i];
                a_species_masses[n][i] = a_species_masses[n][j];
            }
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
#ifndef _WIN32
    struct timeval total_start, total_end;
    gettimeofday(&total_start, NULL);
#endif

    BL_PROFILE("SuperDropletPC::Coalescence()");
    AMREX_ASSERT( a_lev == m_lev );

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const auto domain = geom.Domain();

    const auto rho_w = m_species_mat[m_idx_w]->m_density;
    const auto idx_w = m_idx_w;

    const int num_ae = m_num_aerosols;
    const int num_sp  = m_num_species;
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    const ParticleReal inv_bin_size
        = one / (  static_cast<ParticleReal>(m_coalescence_bin_size[0])
                 * static_cast<ParticleReal>(m_coalescence_bin_size[1])
                 * static_cast<ParticleReal>(m_coalescence_bin_size[2]) );
    const ParticleReal inv_bin_volume = inv_cell_volume*inv_bin_size;

    Real num_collisions = 0;
    const auto& gvec = a_temperature.nGrowVect();

    auto kernel_choice = m_coalescence_kernel;
    auto include_brownian_coalescence = m_include_brownian_coalescence;

    Real mcshuffle_wtime_sec = zero;
    Real mcpairing_wtime_sec = zero;
    Real coalescence_wtime_sec = zero;

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
        int rtoff_i = SuperDropletsIntIdxSoA::ncomps;
        auto* active_ptr = soa.GetIntData(rtoff_i+SuperDropletsIntIdxSoA_RT::active).data();
        int rtoff_r = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::radius).data();
        auto* mult_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::multiplicity).data();
        auto* vterm_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::term_vel).data();

        /* species masses */
        SDSpeciesMassArr sp_mass_ptrs;
        for (int i = 0; i < num_sp; i++) {
            sp_mass_ptrs[i] = soa.GetRealData(idx_s(i,num_ae,num_sp)).data();
        }

        /* aerosol masses */
        SDAerosolMassArr ae_mass_ptrs;
        for (int i = 0; i < num_ae; i++) {
            ae_mass_ptrs[i] = soa.GetRealData(idx_a(i,num_ae,num_sp)).data();
        }

        Gpu::DeviceVector<ParticleReal> sp_density(num_sp);
        Gpu::DeviceVector<int> sp_solubility(num_sp);
        {
            Vector<ParticleReal> sp_density_h(num_sp);
            Vector<int> sp_solubility_h(num_sp);
            for (int i = 0; i < num_sp; i++) {
                sp_density_h[i] = m_species_mat[i]->m_density;
                sp_solubility_h[i] = static_cast<int>(m_species_mat[i]->m_is_soluble);
            }
            Gpu::copy(  Gpu::hostToDevice,
                        sp_density_h.begin(),
                        sp_density_h.end(),
                        sp_density.begin() );
            Gpu::copy(  Gpu::hostToDevice,
                        sp_solubility_h.begin(),
                        sp_solubility_h.end(),
                        sp_solubility.begin() );
        }

        Gpu::DeviceVector<ParticleReal> ae_density(num_ae);
        Gpu::DeviceVector<int> ae_solubility(num_ae);
        {
            Vector<ParticleReal> ae_density_h(num_ae);
            Vector<int> ae_solubility_h(num_ae);
            for (int i = 0; i < num_ae; i++) {
                ae_density_h[i] = m_aerosol_mat[i]->m_density;
                ae_solubility_h[i] = static_cast<int>(m_aerosol_mat[i]->m_is_soluble);
            }
            Gpu::copy(  Gpu::hostToDevice,
                        ae_density_h.begin(),
                        ae_density_h.end(),
                        ae_density.begin() );
            Gpu::copy(  Gpu::hostToDevice,
                        ae_solubility_h.begin(),
                        ae_solubility_h.end(),
                        ae_solubility.begin() );
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

#ifndef _WIN32
        struct timeval mcshuffle_start, mcshuffle_end;
        gettimeofday(&mcshuffle_start, NULL);
#endif

#ifdef AMREX_USE_GPU
        {
            // get the max bin size
            Gpu::Buffer<unsigned int> max_np_bin_d({0});
            auto max_np_bin_ptr = max_np_bin_d.data();
            ParallelFor( bins.numBins(), [=] AMREX_GPU_DEVICE (int i_bin) noexcept
            {
                auto bin_start = offsets[i_bin];
                auto bin_stop = offsets[i_bin+1];
                auto np_bin = bin_stop - bin_start;
                Gpu::Atomic::Max(max_np_bin_ptr, static_cast<unsigned int>(np_bin));
            });
            Gpu::synchronize();
            auto max_np_bin = *(max_np_bin_d.copyToHost());
            // create a stencil vector with integers [0, max_np_bin-1]
            Vector<unsigned int> stencil_vec(max_np_bin);
            for (unsigned int i = 0; i < max_np_bin; i++) { stencil_vec[i] = i; }
            // now shuffle it
            std::shuffle ( stencil_vec.begin(),stencil_vec.end(), m_rndeng );
            // Copy to device
            Gpu::DeviceVector<unsigned int> stencil_vec_d;
            stencil_vec_d.resize(max_np_bin);
            Gpu::copy(  Gpu::hostToDevice,
                        stencil_vec.begin(),
                        stencil_vec.end(),
                        stencil_vec_d.begin() );
            Gpu::synchronize();
            // use this stencil to shuffle each bin
            Gpu::DeviceVector<unsigned int> inds_tmp;
            inds_tmp.resize(np);
            auto stencil_vec_ptr = stencil_vec_d.data();
            auto inds_tmp_ptr = inds_tmp.data();
            ParallelFor( bins.numBins(), [=] AMREX_GPU_DEVICE (int i_bin) noexcept
            {
                auto bin_start = offsets[i_bin];
                auto bin_stop = offsets[i_bin+1];
                auto np_bin = static_cast<unsigned int>(bin_stop - bin_start);
                if (np_bin <= 1) { return; }
                AMREX_ALWAYS_ASSERT(np_bin <= max_np_bin);

                unsigned int j = 0;
                for (unsigned int i = 0; i < max_np_bin; i++) {
                    if (stencil_vec_ptr[i] < np_bin) {
                        inds_tmp_ptr[bin_start+j] = inds[bin_start+stencil_vec_ptr[i]];
                        j++;
                    }
                }
                AMREX_ALWAYS_ASSERT(j == np_bin);
                for (unsigned int i = 0; i < np_bin; i++) {
                    inds[bin_start + i] = inds_tmp_ptr[bin_start + i];
                }
            });
            Gpu::synchronize();
        }
#else
        for (int i_bin = 0; i_bin < bins.numBins(); i_bin++) {
            std::shuffle( inds + offsets[i_bin],
                          inds + offsets[i_bin+1],
                          m_rndeng );
        }
#endif

#ifndef _WIN32
        gettimeofday(&mcshuffle_end,NULL);
        long long mcshuffle_wtime;
        mcshuffle_wtime = (   (mcshuffle_end.tv_sec   * 1000000 + mcshuffle_end.tv_usec  )
                            - (mcshuffle_start.tv_sec * 1000000 + mcshuffle_start.tv_usec) );
        mcshuffle_wtime_sec += (double) mcshuffle_wtime / Real(1000000.0);
#endif

        const auto& pressure_arr = a_pressure[grid].const_array();
        const auto& temperature_arr = a_temperature[grid].const_array();

        CollisionKernel<ParticleReal,AMREX_SPACEDIM> ckernel{};

        Gpu::Buffer<Real> particle_collisions({0});
        auto particle_collisions_ptr = particle_collisions.data();

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
            coal_rate_ptr[i] = zero;
            coal_rmndr_ptr[i] = zero;
        });
        Gpu::synchronize();

#ifndef _WIN32
        struct timeval mcpairing_start, mcpairing_end;
        gettimeofday(&mcpairing_start, NULL);
#endif

        ParallelFor( bins.numBins(), [=] AMREX_GPU_DEVICE (int i_bin) noexcept
        {
            auto bin_start = offsets[i_bin];
            auto bin_stop = offsets[i_bin+1];
            auto np_bin = bin_stop - bin_start;
            if (np_bin <= 1) { return; }

            for (int p = 0; p < np_bin/2; p++) {
                auto pi = inds[bin_start+p];
                auto pj = inds[bin_stop-1-p];

                if (pi == pj) { continue; }
                if (active_ptr[pi] == 0) { continue; }
                if (active_ptr[pj] == 0) { continue; }
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

#ifndef _WIN32
        gettimeofday(&mcpairing_end,NULL);
        long long mcpairing_wtime;
        mcpairing_wtime = (   (mcpairing_end.tv_sec   * 1000000 + mcpairing_end.tv_usec  )
                            - (mcpairing_start.tv_sec * 1000000 + mcpairing_start.tv_usec) );
        mcpairing_wtime_sec += (double) mcpairing_wtime / Real(1000000.0);

        struct timeval coalescence_start, coalescence_end;
        gettimeofday(&coalescence_start, NULL);
#endif

        auto sp_rho_arr = sp_density.data();
        auto sp_sol_arr = sp_solubility.data();
        auto ae_rho_arr = ae_density.data();
        auto ae_sol_arr = ae_solubility.data();

        ParallelForRNG( np, [=] AMREX_GPU_DEVICE (int i, RandomEngine const& rnd_eng) noexcept
        {
            if (partner_idx_ptr[i] < 0) { return; }
            if (!flag_prey_ptr[i]) { return; }

            int pi = i;
            int pj = partner_idx_ptr[i];
            AMREX_ALWAYS_ASSERT(mult_ptr[pi] >= mult_ptr[pj]);

            ParticleReal k_val = zero;
            if (kernel_choice == SDCoalescenceKernelType::golovin) {

                k_val = ckernel.golovin(radius_ptr[pi],radius_ptr[pj]);

            } else {

                ParticleReal v_i[AMREX_SPACEDIM], v_j[AMREX_SPACEDIM];
                for (int d = 0; d < AMREX_SPACEDIM; d++) {
                    v_i[d] = v_ptr[d][pi];
                    v_j[d] = v_ptr[d][pj];
                }
                v_i[AMREX_SPACEDIM-1] -= vterm_ptr[pi];
                v_j[AMREX_SPACEDIM-1] -= vterm_ptr[pj];

                if (kernel_choice == SDCoalescenceKernelType::sedimentation) {
                    k_val = ckernel.sedimentation(radius_ptr[pi],radius_ptr[pj],v_i,v_j);
                } else if (kernel_choice == SDCoalescenceKernelType::Longs) {
                    k_val = ckernel.Longs(radius_ptr[pi],radius_ptr[pj],v_i,v_j);
                } else if (kernel_choice == SDCoalescenceKernelType::Halls) {
                    k_val = ckernel.Halls(radius_ptr[pi],radius_ptr[pj],v_i,v_j);
                }

                if (k_val < zero) {
                    amrex::Abort("Invalid value for k_val");
                }
            }

            if (include_brownian_coalescence) {

                ParticleReal pressure = zero, temperature = zero;
                {
                    ParticleType& par_1 = pstruct_ptr[pi];
                    auto iv = getParticleCell(par_1, plo, dxi, domain);
                    pressure = pressure_arr(iv[0],iv[1],iv[2],0);
                    temperature = temperature_arr(iv[0],iv[1],iv[2],0);
                }

                ParticleReal sd_mass_1 = zero,
                             sd_mass_2 = zero;
                for (int ia = 0; ia < num_ae; ia++) {
                    sd_mass_1 += ae_mass_ptrs[ia][pi];
                    sd_mass_2 += ae_mass_ptrs[ia][pj];
                }
                for (int ia = 0; ia < num_sp; ia++) {
                    sd_mass_1 += sp_mass_ptrs[ia][pi];
                    sd_mass_2 += sp_mass_ptrs[ia][pj];
                }

                auto r_eff_1 = SD_effective_radius( pi, idx_w,
                                                    rho_w,
                                                    num_sp, num_ae,
                                                    sp_sol_arr, ae_sol_arr,
                                                    sp_mass_ptrs, ae_mass_ptrs,
                                                    sp_rho_arr, ae_rho_arr );
                auto r_eff_2 = SD_effective_radius( pj, idx_w,
                                                    rho_w,
                                                    num_sp, num_ae,
                                                    sp_sol_arr, ae_sol_arr,
                                                    sp_mass_ptrs, ae_mass_ptrs,
                                                    sp_rho_arr, ae_rho_arr );

                auto k_brown = ckernel.Brownian_SeinfeldPandis( r_eff_1,
                                                                r_eff_2,
                                                                sd_mass_1,
                                                                sd_mass_2,
                                                                pressure,
                                                                temperature );
                if (k_brown < zero) {
                    amrex::Abort("Invalid value for k_brown");
                }

                k_val += k_brown;
            }


            auto prob_ij = k_val*inv_bin_volume;
            auto prob_sd_ij = std::max(mult_ptr[pi],mult_ptr[pj])*prob_ij;

            auto ns = static_cast<ParticleReal>(np_bin_ptr[i]);
            auto scaling_factor = myhalf*ns*(ns-1)/std::floor(myhalf*ns);
            auto scaled_prob = prob_sd_ij * scaling_factor;

            auto gamma = coalescence_rate ( rnd_eng, (scaled_prob*a_dt) );
            if (gamma > 0) {
                amrex::Gpu::Atomic::Add(particle_collisions_ptr, gamma);
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
                                 num_sp,
                                 sp_mass_ptrs,
                                 num_ae,
                                 ae_mass_ptrs );
        } );
        Gpu::synchronize();

#ifndef _WIN32
        gettimeofday(&coalescence_end,NULL);
        long long coalescence_wtime;
        coalescence_wtime = (  (coalescence_end.tv_sec   * 1000000 + coalescence_end.tv_usec  )
                             - (coalescence_start.tv_sec * 1000000 + coalescence_start.tv_usec) );
        coalescence_wtime_sec += (double) coalescence_wtime / Real(1000000.0);
#endif
    }

    ParallelDescriptor::ReduceRealSum(  &num_collisions,
                                        1,
                                        ParallelDescriptor::IOProcessorNumber() );

#ifndef _WIN32
    gettimeofday(&total_end,NULL);
    long long total_wtime;
    total_wtime = (   (total_end.tv_sec   * 1000000 + total_end.tv_usec  )
                   -  (total_start.tv_sec * 1000000 + total_start.tv_usec) );
    Real total_wtime_sec = (double) total_wtime / Real(1000000.0);

    ParallelDescriptor::ReduceRealMax( &mcshuffle_wtime_sec,
                                       1,
                                       ParallelDescriptor::IOProcessorNumber() );
    ParallelDescriptor::ReduceRealMax( &mcpairing_wtime_sec,
                                       1,
                                       ParallelDescriptor::IOProcessorNumber() );
    ParallelDescriptor::ReduceRealMax( &coalescence_wtime_sec,
                                       1,
                                       ParallelDescriptor::IOProcessorNumber() );
    ParallelDescriptor::ReduceRealMax( &total_wtime_sec,
                                       1,
                                       ParallelDescriptor::IOProcessorNumber() );
#else
    Real total_wtime_sec = zero;
#endif

    Print() << "SuperDropletPC(" << m_name << "): "
            << "number of collisions = " << num_collisions << "\n"
            << "    "
            << "wall time (seconds) = " << total_wtime_sec << " (total), "
            << mcshuffle_wtime_sec << " (MC shuffle), "
            << mcpairing_wtime_sec << " (MC pairing), "
            << coalescence_wtime_sec << " (coalescence)"
            << "\n";
}

#endif


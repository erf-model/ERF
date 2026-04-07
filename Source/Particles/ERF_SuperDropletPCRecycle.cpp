#include <random>
#include <AMReX_PlotFileUtil.H>
#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;
using namespace SDPCDefn;

/*! Recycle deactivated particles: particles that have zero multiplicity are
 *  recycled by resetting them to dry aerosol particles and placing them randomly
 *  in the domain. */
void SuperDropletPC::Recycle ( const int             a_lev,
                               const Vector<MFPtr>&  a_z_phys_nd,
                               const int             a_iter,
                               const Real            a_dt,
                               const bool            a_recycle )
{
    BL_PROFILE("SuperDropletPC::Recycle()");

    const auto num_sd_deactivated = NumSDDeactivated();
    const auto num_sd = NumSuperDroplets();
    const auto deac_frac = static_cast<Real>(num_sd_deactivated) / static_cast<Real>(num_sd);

    int flag = 0;
    if (deac_frac > m_deac_threshold) { flag = 1; }
    ParallelDescriptor::ReduceIntMax( &flag, 1 );
    if (!flag) { return; }

    if (m_save_inactive) {
#ifndef ERF_USE_NETCDF
        const auto& geom = Geom(m_lev);

        using SrcData = SuperDropletPC::ParticleTileType::ConstParticleTileDataType;
        std::string name = "deactivated_particles";
        auto tmp = this->make_alike<amrex::PinnedArenaAllocator>();
        tmp.copyParticles(*this, [=] AMREX_GPU_HOST_DEVICE (const SrcData& a_src, int i)
                          {
                              auto ai = a_src.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                              return (ai == 0);
                          }, true);

        char iter_str[12]; snprintf(iter_str, sizeof(iter_str), "%05d", a_iter+1);
        std::string fname = "deac_SD_" + std::string(iter_str);
        amrex::Print() << "Writing deactivated particles to " << fname << "\n";
        amrex::Vector<std::string> eu_vnames(1,"dummy_var");
        WriteSingleLevelPlotfile(fname, m_mf_buf, eu_vnames, geom, a_iter*a_dt, a_iter );
        {
            std::ofstream jobInfoFile;
            std::string FullPathJobInfoFile = fname;
            FullPathJobInfoFile += "/job_info";
            jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);
            jobInfoFile << "AMReX\n";
            jobInfoFile.close();
        }
        tmp.WritePlotFile(fname, name, this->varNames());
#else
        amrex::Abort("Saving deactivated particles to plotfile not implemented when built with NetCDF support.");
        /* The line
         *      auto tmp = this->make_alike<amrex::PinnedArenaAllocator>();
         * gives the following compilation error on Tuolumne when compiling with NetCDF:
         *
                In file included from /g/g92/ghosh5/Codes/ERF/Source/Particles/ERF_SuperDropletPCRecycle.cpp:3:
                In file included from /g/g92/ghosh5/Codes/ERF/Source/Particles/ERF_SuperDropletPC.H:9:
                In file included from /g/g92/ghosh5/Codes/ERF/Source/Particles/ERF_SuperDropletPCDefinitions.H:8:
                In file included from /g/g92/ghosh5/Codes/ERF/Submodules/AMReX/Src/Particle/AMReX_Particles.H:5:
                /g/g92/ghosh5/Codes/ERF/Submodules/AMReX/Src/Particle/AMReX_ParticleContainer.H:1394:13: error: 'SetSoACompileTimeNames' is a protected member of 'amrex::ParticleContainer_impl<amrex::Particle<0, 1>, 5, 0, amrex::PinnedArenaAllocator>'
                 1394 |         tmp.SetSoACompileTimeNames(real_ct_names, int_ct_names);
                      |             ^
                /g/g92/ghosh5/Codes/ERF/Source/Particles/ERF_SuperDropletPCRecycle.cpp:38:26: note: in instantiation of function template specialization 'amrex::ParticleContainer_impl<amrex::Particle<0, 1>, 5, 0, amrex::DefaultAllocator, ERFParticlesAssignor>::make_alike<amrex::PinnedArenaAllocator>' requested here
                   38 |         auto tpc = this->make_alike<amrex::PinnedArenaAllocator>();
                      |                          ^
                /g/g92/ghosh5/Codes/ERF/Submodules/AMReX/Src/Particle/AMReX_ParticleContainer.H:1474:10: note: declared protected here
                 1474 |     void SetSoACompileTimeNames (std::vector<std::string> const & rdata_name, std::vector<std::string> const & idata_name);
                      |          ^
                1 error generated when compiling for gfx942.

         * Need to understand what is going on. Not sure if it will give this error when compiling on other machines
         * with NetCDF support.
         */
#endif
    }

    if (a_recycle) {

        const int init_idx = Random_int(m_num_initializations);

        Print() << "SuperDropletPC(" << m_name << "): "
                << "recycling " << num_sd_deactivated << " super-droplets "
                << "based on initialization #" << init_idx << ".\n";

        // for multiple initializations, cycle through them each time this
        // function is called.
        const auto init_r = *(m_initializations[init_idx]);
        const auto sampled_multiplicity = init_r.sampledMultiplicity();

        const auto dx_h = Geom(m_lev).CellSize();
        const Real cell_volume = dx_h[0]*dx_h[1]*dx_h[2];

        const int idx_w = m_idx_w;
        const Real rho_w = m_species_mat[m_idx_w]->m_density;
        const int num_sp  = m_num_species;
        const int num_ae = m_num_aerosols;

        // number of super-droplets per cell
        int num_sd_per_cell = m_num_sd_per_cell;
        // number of physical particles per cell
        Real num_par_per_cell = zero;
        for (int i = 0; i < m_num_initializations; i++) {
            num_par_per_cell += m_initializations[i]->numParticlesPerCell(cell_volume);
        }
        // average multiplicity
        auto mult_avg = num_par_per_cell / num_sd_per_cell;

        Long np_recycle_total = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {
            auto& ptile = ParticlesAt(a_lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            auto& soa  = ptile.GetStructOfArrays();
            const int np = aos.numParticles();
            auto *p_pbox = aos().data();

            Array<ParticleReal*,AMREX_SPACEDIM> v_ptr;
            v_ptr[0] = soa.GetRealData(SuperDropletsRealIdxSoA::vx).data();
            v_ptr[1] = soa.GetRealData(SuperDropletsRealIdxSoA::vy).data();
            v_ptr[2] = soa.GetRealData(SuperDropletsRealIdxSoA::vz).data();
            auto* mass_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::mass).data();

            int rtoff_i = SuperDropletsIntIdxSoA::ncomps;
            auto* active_ptr = soa.GetIntData(rtoff_i+SuperDropletsIntIdxSoA_RT::active).data();
            int rtoff_r = SuperDropletsRealIdxSoA::ncomps;
            auto* radius_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::radius).data();
            auto* vterm_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::term_vel).data();
            auto* mult_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::multiplicity).data();

            SDSpeciesMassArr sp_mass_ptrs;
            for (int i = 0; i < num_sp; i++) {
                sp_mass_ptrs[i] = soa.GetRealData(idx_s(i,num_ae,num_sp)).data();
            }

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

            auto sp_rho_arr = sp_density.data();
            auto sp_sol_arr = sp_solubility.data();
            auto ae_rho_arr = ae_density.data();
            auto ae_sol_arr = ae_solubility.data();

            // get sampled aerosol mass values based on initialization
            Gpu::DeviceVector<Real> aerosol_mass_d(num_ae*np);
            Gpu::DeviceVector<Real> multiplicity_d(np);
            ParticleReal mult_scale = one;
            {
                Vector<Real> multiplicity_h(np, zero);
                for (int i = 0; i < num_ae; i++) {
                    Vector<Real> aerosol_mass_h;
                    if (sampled_multiplicity) {
                        init_r.getAerosolDistribution( aerosol_mass_h,
                                                       multiplicity_h,
                                                       cell_volume,
                                                       i,
                                                       np,
                                                       m_aerosol_mat[i]->m_density,
                                                       m_rndeng );
                    } else {
                        init_r.getAerosolDistribution( aerosol_mass_h,
                                                       i,
                                                       np,
                                                       m_aerosol_mat[i]->m_density,
                                                       m_rndeng );
                    }
                    Gpu::copy( Gpu::hostToDevice,
                               aerosol_mass_h.begin(),
                               aerosol_mass_h.end(),
                               aerosol_mass_d.begin() + (i*np) );
                }
                if (sampled_multiplicity) {
                    // compute multiplicity scale
                    ParticleReal mult_sum = zero;
                    for (int ctr=0; ctr < multiplicity_h.size(); ctr++) {
                        mult_sum += multiplicity_h[ctr];
                    }
                    mult_scale = np * mult_avg / mult_sum;
                    // copy to GPU
                    Gpu::copy( Gpu::hostToDevice,
                               multiplicity_h.begin(),
                               multiplicity_h.end(),
                               multiplicity_d.begin() );
                }
            }
            Gpu::synchronize();

            auto aerosol_mass = aerosol_mass_d.data();
            auto mult_arr = multiplicity_d.data();

            Gpu::Buffer<Long> np_recycle_buf({0});
            auto* np_recycle_ptr = np_recycle_buf.data();

            const auto x_min = m_recyc_xmin;
            const auto x_max = m_recyc_xmax;
            const auto y_min = m_recyc_ymin;
            const auto y_max = m_recyc_ymax;
            const auto z_min = m_recyc_zmin;
            const auto z_max = m_recyc_zmax;

            ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int i, const RandomEngine& rnd_engine) noexcept
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }
                if (active_ptr[i] > 0) { return; }

                // Place particle randomly in domain within specified bounds
                p.pos(0) = x_min + Random(rnd_engine)*(x_max - x_min);
                p.pos(1) = y_min + Random(rnd_engine)*(y_max - y_min);
                p.pos(2) = z_min + Random(rnd_engine)*(z_max - z_min);

                // Set velocities to zero
                v_ptr[0][i] = v_ptr[1][i] = v_ptr[2][i] = vterm_ptr[i] = zero;

                // reset all species masses to zero
                for (int ctr = 0; ctr < num_sp; ctr++) {
                    sp_mass_ptrs[ctr][i] = zero;
                }
                // Reset water mass
                auto water_radius = Real(1.0e-15);
                auto water_mass = (Real(4.0)/three)*PI
                                 * water_radius*water_radius*water_radius*rho_w;
                sp_mass_ptrs[idx_w][i] = water_mass;

                // choose a random index
                auto j = Random_int(np, rnd_engine);
                // Set aerosol mass
                for (int ctr = 0; ctr < num_ae; ctr++) {
                    ae_mass_ptrs[ctr][i] = aerosol_mass[ctr*np+j];
                }
                // Set multiplicity to sampled or averaged multiplicity
                if (sampled_multiplicity) { mult_ptr[i] = std::ceil(mult_arr[j] * mult_scale); }
                else { mult_ptr[i] = mult_avg; }

                // compute effective radius and total mass
                radius_ptr[i] = SD_effective_radius( i, idx_w,
                                                     rho_w,
                                                     num_sp, num_ae,
                                                     sp_sol_arr, ae_sol_arr,
                                                     sp_mass_ptrs, ae_mass_ptrs,
                                                     sp_rho_arr, ae_rho_arr );
                mass_ptr[i] = SD_total_mass( i, num_sp, num_ae, sp_mass_ptrs, ae_mass_ptrs);

                // set as active
                active_ptr[i] = 1;

                // add to count
                Gpu::Atomic::Add(np_recycle_ptr, Long(1));
            });
            Gpu::synchronize();
            np_recycle_total += *(np_recycle_buf.copyToHost());
        }

        ParallelDescriptor::ReduceLongSum( &np_recycle_total, 1 );
        Print() << "    recycled " << np_recycle_total << " super-droplets.\n";

        Redistribute();

        const MFPtr& z_height = a_z_phys_nd[a_lev];
        const auto plo = Geom(m_lev).ProbLoArray();
        const auto dxi = Geom(m_lev).InvCellSizeArray();
        for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {
            int grid    = pti.index();
            auto& ptile = ParticlesAt(a_lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            const int n = aos.numParticles();
            auto *p_pbox = aos().data();
            auto zheight = (*z_height)[grid].array();

            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }
                update_location_idata(p,plo,dxi,zheight);

            });
            Gpu::synchronize();
        }

    } else {

        amrex::Print() << "Removing inactive particles.\n";

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {
            auto& ptile = ParticlesAt(a_lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            auto& soa  = ptile.GetStructOfArrays();
            const int np = aos.numParticles();
            auto *p_pbox = aos().data();

            int rtoff_i = SuperDropletsIntIdxSoA::ncomps;
            auto* active_ptr = soa.GetIntData(rtoff_i+SuperDropletsIntIdxSoA_RT::active).data();

            ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int i, const RandomEngine& /*rnd_engine*/) noexcept
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }
                if (active_ptr[i] == 0) { p.id() = -1; }
            });
        }

        Redistribute();
    }
}

#endif


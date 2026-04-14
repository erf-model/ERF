#include <random>
#include <AMReX_PlotFileUtil.H>
#include "ERF_SuperDropletPC.H"
#include "ERF_InterpolationUtils.H"

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
    const auto deac_frac = (num_sd > 0 ? static_cast<Real>(num_sd_deactivated)/static_cast<Real>(num_sd) : 0.0);

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

        const auto ctx = buildProcessContext(a_lev);
        const auto dx_h = Geom(m_lev).CellSize();
        const Real cell_volume = dx_h[0]*dx_h[1]*dx_h[2];

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

        const auto x_min = m_recyc_xmin;
        const auto x_max = m_recyc_xmax;
        const auto y_min = m_recyc_ymin;
        const auto y_max = m_recyc_ymax;
        const auto z_min = m_recyc_zmin;
        const auto z_max = m_recyc_zmax;

        forEachParticleTile(a_lev, ctx,
            [&](ParIterType& /*pti*/, int /*grid*/, ParticleType* p_pbox,
                const SDProcess::ParticlePointers& ptrs,
                const SDProcess::ProcessContext& ctx)
        {
            const int np = ptrs.num_particles;

            // Get sampled aerosol mass values based on initialization
            Gpu::DeviceVector<Real> aerosol_mass_d(ctx.num_aerosols*np);
            Gpu::DeviceVector<Real> multiplicity_d(np);
            ParticleReal mult_scale = one;
            {
                Vector<Real> multiplicity_h(np, zero);
                for (int i = 0; i < ctx.num_aerosols; i++) {
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

            ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int i, const RandomEngine& rnd_engine) noexcept
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }
                if (ptrs.active_ptr[i] > 0) { return; }

                // Place particle randomly in domain within specified bounds
                p.pos(0) = x_min + Random(rnd_engine)*(x_max - x_min);
                p.pos(1) = y_min + Random(rnd_engine)*(y_max - y_min);
                p.pos(2) = z_min + Random(rnd_engine)*(z_max - z_min);

                // Set velocities to zero
                ptrs.v_ptr[0][i] = ptrs.v_ptr[1][i] = ptrs.v_ptr[2][i] = ptrs.vterm_ptr[i] = zero;

                // reset all species masses to zero
                for (int ctr = 0; ctr < ctx.num_species; ctr++) {
                    ptrs.sp_mass_ptrs[ctr][i] = zero;
                }
                // Reset water mass
                auto water_radius = Real(1.0e-15);
                auto water_mass = four_thirds_pi
                                 * water_radius*water_radius*water_radius*ctx.rho_water;
                ptrs.sp_mass_ptrs[ctx.idx_water][i] = water_mass;

                // choose a random index
                auto j = Random_int(np, rnd_engine);
                // Set aerosol mass
                for (int ctr = 0; ctr < ctx.num_aerosols; ctr++) {
                    ptrs.ae_mass_ptrs[ctr][i] = aerosol_mass[ctr*np+j];
                }
                // Set multiplicity to sampled or averaged multiplicity
                if (sampled_multiplicity) { ptrs.mult_ptr[i] = std::ceil(mult_arr[j] * mult_scale); }
                else { ptrs.mult_ptr[i] = mult_avg; }

                // compute effective radius and total mass
                SuperDropletPC::updateParticleAttributes(
                    i, ptrs.radius_ptr, ptrs.mass_ptr, ctx.idx_water, ctx.rho_water,
                    ctx.num_species, ctx.num_aerosols, ptrs.sp_sol_arr, ptrs.ae_sol_arr,
                    ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs, ptrs.sp_rho_arr, ptrs.ae_rho_arr);

                // set as active
                ptrs.active_ptr[i] = 1;

                // add to count
                Gpu::Atomic::Add(np_recycle_ptr, Long(1));
            });
            Gpu::synchronize();
            np_recycle_total += *(np_recycle_buf.copyToHost());
        }); // end forEachParticleTile

        ParallelDescriptor::ReduceLongSum( &np_recycle_total, 1 );
        Print() << "    recycled " << np_recycle_total << " super-droplets.\n";

        Redistribute();

        // Update location data for recycled particles
        const MFPtr& z_height = a_z_phys_nd[a_lev];
        forEachParticleTile(a_lev, ctx,
            [&](ParIterType& /*pti*/, int grid, ParticleType* p_pbox,
                const SDProcess::ParticlePointers& ptrs,
                const SDProcess::ProcessContext& ctx)
        {
            auto zheight = (*z_height)[grid].array();

            ParallelFor(ptrs.num_particles, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }
                update_location_idata(p,ctx.plo,ctx.dxi,zheight);
            });
            Gpu::synchronize();
        }); // end forEachParticleTile

    } else {

        amrex::Print() << "Removing inactive particles.\n";

        const auto ctx = buildProcessContext(a_lev);

        forEachParticleTile(a_lev, ctx,
            [&](ParIterType& /*pti*/, int /*grid*/, ParticleType* p_pbox,
                const SDProcess::ParticlePointers& ptrs,
                const SDProcess::ProcessContext& /*ctx*/)
        {
            ParallelForRNG(ptrs.num_particles, [=] AMREX_GPU_DEVICE (int i, const RandomEngine& /*rnd_engine*/) noexcept
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }
                if (ptrs.active_ptr[i] == 0) { p.id() = -1; }
            });
        }); // end forEachParticleTile

        Redistribute();
    }
}

#endif

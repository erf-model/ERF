#include <random>
#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! Recycle deactivated particles: particles that have zero multiplicity are
 *  recycled by resetting them to dry aerosol particles and placing them randomly
 *  in the domain. */
void SuperDropletPC::Recycle ( const int             a_lev,
                               const Vector<MFPtr>&  a_z_phys_nd )
{
    BL_PROFILE("SuperDropletPC::Recycle()");

    const auto num_sd_deactivated = NumSDDeactivated();
    const auto num_sd = NumSuperDroplets();
    const auto deac_frac = static_cast<Real>(num_sd_deactivated) / static_cast<Real>(num_sd);

    int flag = 0;
    if (deac_frac > m_recycle_threshold) { flag = 1; }
    ParallelDescriptor::ReduceIntMax( &flag, 1 );
    if (!flag) { return; }

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
    Real num_par_per_cell = 0.0;
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
        int grid    = pti.index();
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

        int rt_offset = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::radius).data();
        auto* vterm_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::term_vel).data();
        auto* mult_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::multiplicity).data();

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
        ParticleReal mult_scale = 1.0;
        {
            Vector<Real> multiplicity_h(np, 0.0);
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
                ParticleReal mult_sum = 0.0;
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
            if (mult_ptr[i] > 0) { return; }

            // Place particle randomly in domain within specified bounds
            p.pos(0) = x_min + Random(rnd_engine)*(x_max - x_min);
            p.pos(1) = y_min + Random(rnd_engine)*(y_max - y_min);
            p.pos(2) = z_min + Random(rnd_engine)*(z_max - z_min);

            // Set velocities to zero
            v_ptr[0][i] = v_ptr[1][i] = v_ptr[2][i] = vterm_ptr[i] = 0.0;

            // reset all species masses to zero
            for (int ctr = 0; ctr < num_sp; ctr++) {
                sp_mass_ptrs[ctr][i] = 0.0;
            }
            // Reset water mass
            auto water_radius = 1.0e-15;
            auto water_mass = (4.0/3.0)*PI
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
}

#endif


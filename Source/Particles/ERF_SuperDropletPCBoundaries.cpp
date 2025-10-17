#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include <AMReX_TracerParticle_mod_K.H>

using namespace amrex;

/*! Handle the boundaries for the particles */
void SuperDropletPC::applyBoundaryTreatment ( int                   a_lev,
                                              const Vector<MFPtr>&  a_z_phys_nd,
                                              const BCTypeArr&      a_bctypes )
{
    BL_PROFILE("SuperDropletPC::applyBoundaryTreatment()");
    const MFPtr& z_height = a_z_phys_nd[a_lev];

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto dx_h = Geom(m_lev).CellSize();
    const Real cell_volume = dx_h[0]*dx_h[1]*dx_h[2];
    const auto dx = geom.CellSizeArray();
    const auto dxi = geom.InvCellSizeArray();
    const auto domain = geom.Domain();
    const auto is_periodic = geom.isPeriodicArray();

    const int k_lo = domain.smallEnd(2);
    const int k_hi = domain.bigEnd(2);

    const int n_sp = m_num_species;
    const int n_ae = m_num_aerosols;
    const Real rho_w = m_species_mat[m_idx_w]->m_density;
    const int idx_w = m_idx_w;

    // number of super-droplets per cell
    int num_sd_per_cell = m_num_sd_per_cell;
    // number of physical particles per cell
    Real num_par_per_cell = 0.0;
    for (int i = 0; i < m_num_initializations; i++) {
        num_par_per_cell += m_initializations[i]->numParticlesPerCell(cell_volume);
    }
    auto multiplicity = num_par_per_cell / num_sd_per_cell;

    Long num_deactivated_particles = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {
        int grid    = pti.index();
        auto& ptile = ParticlesAt(a_lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        auto& soa  = ptile.GetStructOfArrays();
        const int n = aos.numParticles();
        auto *p_pbox = aos().data();

        Array<ParticleReal*,AMREX_SPACEDIM> v_ptr;
        v_ptr[0] = soa.GetRealData(SuperDropletsRealIdxSoA::vx).data();
        v_ptr[1] = soa.GetRealData(SuperDropletsRealIdxSoA::vy).data();
        v_ptr[2] = soa.GetRealData(SuperDropletsRealIdxSoA::vz).data();
        auto* mass_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::mass).data();

        auto zheight = (*z_height)[grid].array();

        int rtoff_i = SuperDropletsIntIdxSoA::ncomps;
        auto* active_ptr = soa.GetIntData(rtoff_i+SuperDropletsIntIdxSoA_RT::active).data();
        int rtoff_r = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::radius).data();
        auto* vterm_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::term_vel).data();
        auto* mult_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::multiplicity).data();

        auto* a_ptr = soa.GetRealData(idx_ice_a(n_ae,n_sp)).data();
        auto* c_ptr = soa.GetRealData(idx_ice_c(n_ae,n_sp)).data();
        auto* mrime_ptr = soa.GetRealData(idx_ice_mrime(n_ae,n_sp)).data();
        auto* nmono_ptr = soa.GetRealData(idx_ice_nmono(n_ae,n_sp)).data();

        Gpu::Buffer<Long> deactivated_particles({0});
        auto* deactivated_particles_ptr = deactivated_particles.data();

        SDSpeciesMassArr sp_mass_ptrs;
        Gpu::DeviceVector<ParticleReal> sp_density(n_sp);
        Gpu::DeviceVector<int> sp_solubility(n_sp);
        {
            Vector<ParticleReal> sp_density_h(n_sp);
            Vector<int> sp_solubility_h(n_sp);
            for (int i = 0; i < n_sp; i++) {
                sp_mass_ptrs[i] = soa.GetRealData(idx_s(i,n_ae,n_sp)).data();
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

        SDAerosolMassArr ae_mass_ptrs;
        Gpu::DeviceVector<ParticleReal> ae_density(n_ae);
        Gpu::DeviceVector<int> ae_solubility(n_ae);
        {
            Vector<ParticleReal> ae_density_h(n_ae);
            Vector<int> ae_solubility_h(n_ae);
            for (int i = 0; i < n_ae; i++) {
                ae_mass_ptrs[i] = soa.GetRealData(idx_a(i,n_ae,n_sp)).data();
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

        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (active_ptr[i] == 0) { return; }

            // check for ground impact
            {
                auto z_ground = plo[2];
                {
                    auto iv = getParticleCell(p, plo, dxi, domain);
                    z_ground = zheight(iv[0],iv[1],k_lo);
                }
                if (p.pos(2) < z_ground) {
                    p.pos(2) = z_ground + 0.01*dx[2];
                    v_ptr[0][i] = v_ptr[1][i] = v_ptr[2][i] = vterm_ptr[i] = 0.0;
                    active_ptr[i] = 0;
                    Gpu::Atomic::Add(deactivated_particles_ptr, Long(1));
                    update_location_idata(p,plo,dxi,zheight);
                }
            }

            // check for top boundary exits
            {
                auto z_roof = phi[2];
                {
                    auto iv = getParticleCell(p, plo, dxi, domain);
                    z_roof = zheight(iv[0],iv[1],k_hi+1);
                }
                if (p.pos(2) > z_roof) {
                    p.pos(2) = z_roof - dx[2];
                    v_ptr[0][i] = v_ptr[1][i] = v_ptr[2][i] = vterm_ptr[i] = 0.0;
                    active_ptr[i] = 0;
                    Gpu::Atomic::Add(deactivated_particles_ptr, Long(1));
                    update_location_idata(p,plo,dxi,zheight);
                }
            }

            // check if particles have exited the domain along x and y
            for (int d = 0; d < 2; d++) {

                // domain bounds
                auto x_min = plo[d];
                auto x_max = phi[d];

                auto bc_lo = a_bctypes[Orientation(d,Orientation::low)];
                auto bc_hi = a_bctypes[Orientation(d,Orientation::high)];

                if (p.pos(d) < x_min) {

                    if ((bc_lo == ERF_BC::slip_wall) || (bc_lo == ERF_BC::no_slip_wall)) {

                        p.pos(d) = x_min + 0.01*dx[d];
                        v_ptr[0][i] = v_ptr[1][i] = v_ptr[2][i] = vterm_ptr[i] = 0.0;
                        active_ptr[i] = 0;
                        Gpu::Atomic::Add(deactivated_particles_ptr, Long(1));

                    } else {

                        auto delta = x_min - p.pos(d);
                        p.pos(d) = x_max - delta;
                        if (!is_periodic[d]) {
                            v_ptr[0][i] = v_ptr[1][i] = v_ptr[2][i] = vterm_ptr[i] = 0.0;

                            for (int ctr = 0; ctr < n_sp; ctr++) {
                                sp_mass_ptrs[ctr][i] = 0.0;
                            }

                            ParticleReal aerosol_mass_total = 0.0;
                            for (int ctr = 0; ctr < n_ae; ctr++) {
                                aerosol_mass_total += ae_mass_ptrs[ctr][i];
                            }

                            auto par_radius = 1.0e-15;
                            auto cond_mass = (4.0/3.0)*PI
                                             * par_radius*par_radius*par_radius*rho_w;
                            sp_mass_ptrs[idx_w][i] = cond_mass;

                            mult_ptr[i] = multiplicity;
                            a_ptr[i] = 0.0;
                            c_ptr[i] = 0.0;
                            mrime_ptr[i] = 0.0;
                            nmono_ptr[i] = 0.0;

                            radius_ptr[i] = SD_effective_radius( i, idx_w,
                                                                 rho_w,
                                                                 n_sp, n_ae,
                                                                 sp_sol_arr, ae_sol_arr,
                                                                 sp_mass_ptrs, ae_mass_ptrs,
                                                                 sp_rho_arr, ae_rho_arr );
                            mass_ptr[i] = SD_total_mass( i, n_sp, n_ae, sp_mass_ptrs, ae_mass_ptrs);
                        }

                    }

                } else if (p.pos(d) > x_max) {

                    if ((bc_hi == ERF_BC::slip_wall) || (bc_hi == ERF_BC::no_slip_wall)) {

                        p.pos(d) = x_max - 0.01*dx[d];
                        v_ptr[0][i] = v_ptr[1][i] = v_ptr[2][i] = vterm_ptr[i] = 0.0;
                        active_ptr[i] = 0;
                        Gpu::Atomic::Add(deactivated_particles_ptr, Long(1));

                    } else {

                        auto delta = p.pos(d) - x_max;
                        p.pos(d) = x_min + delta;
                        if (!is_periodic[d]) {
                            v_ptr[0][i] = v_ptr[1][i] = v_ptr[2][i] = vterm_ptr[i] = 0.0;

                            for (int ctr = 0; ctr < n_sp; ctr++) {
                                sp_mass_ptrs[ctr][i] = 0.0;
                            }

                            ParticleReal aerosol_mass_total = 0.0;
                            for (int ctr = 0; ctr < n_ae; ctr++) {
                                aerosol_mass_total += ae_mass_ptrs[ctr][i];
                            }

                            auto par_radius = 1.0e-15;
                            auto cond_mass = (4.0/3.0)*PI
                                             * par_radius*par_radius*par_radius
                                             * rho_w;

                            mult_ptr[i] = multiplicity;
                            a_ptr[i] = 0.0;
                            c_ptr[i] = 0.0;
                            mrime_ptr[i] = 0.0;
                            nmono_ptr[i] = 0.0;

                            radius_ptr[i] = SD_effective_radius( i, idx_w,
                                                                 rho_w,
                                                                 n_sp, n_ae,
                                                                 sp_sol_arr, ae_sol_arr,
                                                                 sp_mass_ptrs, ae_mass_ptrs,
                                                                 sp_rho_arr, ae_rho_arr );
                            mass_ptr[i] = SD_total_mass( i, n_sp, n_ae, sp_mass_ptrs, ae_mass_ptrs);
                        }

                    }

                }

            }
        });
        Gpu::synchronize();
        num_deactivated_particles += *(deactivated_particles.copyToHost());
    }

    ParallelDescriptor::ReduceLongSum( &num_deactivated_particles, 1 );
    Print() << "SuperDropletPC(" << m_name << "): "
            << "deactivated " << num_deactivated_particles << " super-droplets at boundaries.\n";
}

#endif

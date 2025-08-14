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

    auto num_sd_deactivated = NumSDDeactivated();
    Print() << "SuperDropletPC(" << m_name << "): "
            << "recycling " << num_sd_deactivated << " super-droplets.\n";

    const MFPtr& z_height = a_z_phys_nd[a_lev];
    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto dxi = geom.InvCellSizeArray();
    const auto domain = geom.Domain();

    const int k_lo = domain.smallEnd(2);
    const int k_hi = domain.bigEnd(2);

    const auto dx_h = Geom(m_lev).CellSize();
    const Real cell_volume = dx_h[0]*dx_h[1]*dx_h[2];

    const int idx_w = m_idx_w;
    const Real rho_w = m_species_mat[m_idx_w]->m_density;
    const int n_species  = m_num_species;
    const int n_aerosols = m_num_aerosols;

    // number of super-droplets per cell
    int num_sd_per_cell = m_num_sd_per_cell;
    // number of physical particles per cell
    Real num_par_per_cell = 0.0;
    for (int i = 0; i < m_num_initializations; i++) {
        num_par_per_cell += m_initializations[i]->numParticlesPerCell(cell_volume);
    }
    auto multiplicity = num_par_per_cell / num_sd_per_cell;

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

        int rt_offset = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::radius).data();
        auto* vterm_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::term_vel).data();
        auto* mult_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::multiplicity).data();

        SDSpeciesMassArr species_mass_ptrs;
        for (int i = 0; i < n_species; i++) {
            species_mass_ptrs[i] = soa.GetRealData(idx_s(i,n_aerosols,n_species)).data();
        }

        SDAerosolMassArr aerosol_mass_ptrs;
        for (int i = 0; i < n_aerosols; i++) {
            aerosol_mass_ptrs[i] = soa.GetRealData(idx_a(i,n_aerosols,n_species)).data();
        }

        ParallelForRNG(n, [=] AMREX_GPU_DEVICE (int i, const RandomEngine& rnd_engine) noexcept
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (mult_ptr[i] > 0) { return; }

            auto x_min = plo[0];
            auto x_max = phi[0];
            auto y_min = plo[1];
            auto y_max = phi[1];
            auto z_min = plo[2];
            auto z_max = phi[2];
            {
                // probably not the right thing to do
                auto iv = getParticleCell(p, plo, dxi, domain);
                z_min = zheight(iv[0],iv[1],k_lo);
                z_max = zheight(iv[0],iv[1],k_hi+1);
            }

            // Place particle randomly in domain
            p.pos(0) = x_min + Random(rnd_engine)*(x_max - x_min);
            p.pos(1) = y_min + Random(rnd_engine)*(y_max - y_min);
            p.pos(2) = z_min + Random(rnd_engine)*(z_max - z_min);

            for (int ctr = 0; ctr < n_species; ctr++) {
                species_mass_ptrs[ctr][i] = 0.0;
            }
            // Reset radius and condensate mass
            auto par_radius = 1.0e-15;
            radius_ptr[i] = par_radius;
            auto cond_mass = (4.0/3.0)*PI
                             * par_radius*par_radius*par_radius*rho_w;
            species_mass_ptrs[idx_w][i] = cond_mass;
            ParticleReal aerosol_mass_total = 0.0;
            for (int ctr = 0; ctr < n_aerosols; ctr++) {
                aerosol_mass_total += aerosol_mass_ptrs[ctr][i];
            }
            mass_ptr[i] = cond_mass + aerosol_mass_total;

            // Set multiplicity to an averaged multiplicity
            mult_ptr[i] = multiplicity;

            // Set velocities to zero
            v_ptr[0][i] = v_ptr[1][i] = v_ptr[2][i] = vterm_ptr[i] = 0.0;
        });
        Gpu::synchronize();
    }

    Redistribute();
}

#endif


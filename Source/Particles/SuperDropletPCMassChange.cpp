#include "SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuBuffer.H>
#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_Constants.H"

using namespace amrex;

/*! Compute mass change of particles due to evaporation and condensation */
void SuperDropletPC::MassChange ( int                                         a_lev,
                                  Real                                        a_dt,
                                  const MultiFab&                             a_temperature,
                                  const MultiFab&                             a_sat_pressure,
                                  const MultiFab&                             a_sat_ratio,
                                  const Vector<std::unique_ptr<MultiFab>>&    a_z_phys_nd )
{
    BL_PROFILE("SuperDropletPC::MassChange()");
    AMREX_ASSERT( a_lev == m_lev );

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    const std::unique_ptr<MultiFab>& z_height = a_z_phys_nd[a_lev];

    const int num_aerosols = m_num_aerosols;
    const Real mat_density = m_vapour_mat->density();
    const std::string vap_name = m_vapour_mat->name();

    long num_unconverged_particles = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {
        int grid = pti.index();
        auto& ptile = ParticlesAt(a_lev, pti);
        auto& aos = ptile.GetArrayOfStructs();
        auto& soa = ptile.GetStructOfArrays();
        const int num_particles = aos.numParticles();
        auto* p_pbox = aos().data();

        bool use_terrain = (z_height != nullptr);
        auto zheight = use_terrain ? (*z_height)[grid].array() : Array4<Real>{};

        auto* mass_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::mass).data();

        int rt_offset = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::radius).data();
        auto* mult_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::multiplicity).data();
        auto* supdrop_mass_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::sd_mass).data();

        SDAerosolMassArr aerosol_mass_ptrs;
        for (int i = 0; i < num_aerosols; i++) {
            aerosol_mass_ptrs[i] = soa.GetRealData(   rt_offset
                                                    + SuperDropletsRealIdxSoA_RT::ncomps
                                                    + i ).data();
        }

        const auto& sat_pressure_arr = a_sat_pressure[grid].array();
        const auto& sat_ratio_arr = a_sat_ratio[grid].array();
        const auto& temperature_arr = a_temperature[grid].array();

        SuperDropletsUtils::dRsqdt_RHSFunc drsqdt_rhsfun{m_vapour_mat->coeffCurv(),
                                                         m_vapour_mat->coeffVPSolute(*m_aerosol_mat[0]),
                                                         m_vapour_mat->latHeatVap(),
                                                         therco, // ERF_Constants.H
                                                         m_vapour_mat->Rv(),
                                                         mat_density,
                                                         m_vapour_mat->molDiffCoeff()};

        SuperDropletsUtils::dRsqdt_RHSJac drsqdt_rhsjac{m_vapour_mat->coeffCurv(),
                                                        m_vapour_mat->coeffVPSolute(*m_aerosol_mat[0]),
                                                        m_vapour_mat->latHeatVap(),
                                                        therco, // ERF_Constants.H
                                                        m_vapour_mat->Rv(),
                                                        mat_density,
                                                        m_vapour_mat->molDiffCoeff()};

        SuperDropletsUtils::NewtonSolver< SuperDropletsUtils::dRsqdt_RHSFunc,
                                          SuperDropletsUtils::dRsqdt_RHSJac,
                                          ParticleReal > newton_solver { drsqdt_rhsfun, drsqdt_rhsjac,
                                                                         m_newton_rtol,
                                                                         m_newton_atol,
                                                                         m_newton_stol,
                                                                         m_newton_maxits };

        amrex::Gpu::Buffer<amrex::Long> unconverged_particles({0});
        amrex::Long* unconverged_particles_ptr = unconverged_particles.data();

        ParallelFor(num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }

            ParticleReal sat_ratio, e_sat, temperature;
            if (use_terrain) {
                cic_interpolate_mapped_z( p, plo, dxi, sat_pressure_arr, zheight, &e_sat, 1 );
                cic_interpolate_mapped_z( p, plo, dxi, sat_ratio_arr, zheight, &sat_ratio, 1 );
                cic_interpolate_mapped_z( p, plo, dxi, temperature_arr, zheight, &temperature, 1 );
            } else {
                cic_interpolate( p, plo, dxi, sat_pressure_arr, &e_sat, 1 );
                cic_interpolate( p, plo, dxi, sat_ratio_arr, &sat_ratio, 1 );
                cic_interpolate( p, plo, dxi, temperature_arr, &temperature, 1 );
            }

            ParticleReal solute_mass = 0.0;
            for (int j = 0; j < num_aerosols; j++) { solute_mass += aerosol_mass_ptrs[j][i]; }

            ParticleReal r_sq = radius_ptr[i]*radius_ptr[i];
            ParticleReal r_sq_0 = r_sq;

            bool converged = false;
            newton_solver ( r_sq, r_sq_0,
                            a_dt,
                            sat_ratio, temperature, e_sat, solute_mass,
                            converged );

            if (!converged) {

                amrex::Gpu::Atomic::Add(unconverged_particles_ptr, amrex::Long(1));

            } else {

                // update particle radius
                radius_ptr[i] = std::sqrt(r_sq);

                // update mass of particle
                mass_ptr[i] = (4.0/3.0)*PI*r_sq*radius_ptr[i]*mat_density;

                // update superdroplet total mass
                supdrop_mass_ptr[i] = mass_ptr[i] * mult_ptr[i];
            }

        });

        num_unconverged_particles = *(unconverged_particles.copyToHost());
    }

    ParallelDescriptor::ReduceLongSum(  &num_unconverged_particles,
                                        1,
                                        ParallelDescriptor::IOProcessorNumber() );
    if (num_unconverged_particles > 0) {
        Print() << "SuperDropletPC::MassChange(): Warning - " << num_unconverged_particles
                << " particles did not converge during Newton solve.\n";
    }
}

#endif

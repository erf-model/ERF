#include "SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuBuffer.H>
#include <AMReX_TracerParticle_mod_K.H>
#include "IndexDefines.H"
#include "Derive.H"
#include "EOS.H"
#include "ERF_Constants.H"

using namespace amrex;

/*! Compute mass change of particles due to evaporation and condensation */
void SuperDropletPC::MassChange (   int                                         a_lev,
                                    Real                                        a_dt,
                                    Vector<Vector<MultiFab>>&                   a_flow_vars,
                                    MultiFab&                                   a_qv,
                                    const Vector<std::unique_ptr<MultiFab>>&    a_z_phys_nd )
{
    BL_PROFILE("SuperDropletPC::massChange()");
    AMREX_ASSERT( a_lev == m_lev );

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    MultiFab* cons_vars( &a_flow_vars[a_lev][Vars::cons] );
    const std::unique_ptr<MultiFab>& z_height = a_z_phys_nd[a_lev];

    MultiFab mf_sat_pressure(   cons_vars->boxArray(),
                                cons_vars->DistributionMap(),
                                1,
                                cons_vars->nGrow() );

    MultiFab mf_sat_ratio(  cons_vars->boxArray(),
                            cons_vars->DistributionMap(),
                            1,
                            cons_vars->nGrow() );

   MultiFab mf_temperature( cons_vars->boxArray(),
                            cons_vars->DistributionMap(),
                            1,
                            cons_vars->nGrow() );

    {
        MultiFab mf_pressure( cons_vars->boxArray(),
                              cons_vars->DistributionMap(),
                              1,
                              cons_vars->nGrow() );

        for (MFIter mfi(mf_temperature, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();
            auto& tfab = mf_temperature[mfi];
            auto& cfab = (*cons_vars)[mfi];

            derived::erf_dertemp( bx, tfab, 0, 1, cfab, geom, 0.0, nullptr, a_lev );

            const Array4<Real      >& p_arr  = mf_pressure.array(mfi);
            const Array4<Real const>& S_arr  = cons_vars->const_array(mfi);
            const Array4<Real const>& qv_arr = a_qv.const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            { p_arr(i,j,k,0) = getPgivenRTh(S_arr(i,j,k,RhoTheta_comp),qv_arr(i,j,k,0)); });
        }

        AMREX_ASSERT( !mf_pressure.contains_nan() );
        AMREX_ASSERT( !mf_temperature.contains_nan() );

        m_vapour_mat->computeSaturationPressure( mf_sat_pressure, mf_temperature );
        m_vapour_mat->computeSaturationVapFrac( mf_sat_ratio, mf_temperature, mf_pressure );

        for (MFIter mfi(mf_temperature, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();
            const Array4<Real>& sr_arr = mf_sat_ratio.array(mfi);
            const Array4<Real const>& qv_arr = a_qv.const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            { sr_arr(i,j,k,0) = qv_arr(i,j,k,0) / sr_arr(i,j,k,0); });
        }

    }

    const int num_aerosols = m_num_aerosols;
    const Real mat_density = m_vapour_mat->density();
    const std::string vap_name = m_vapour_mat->name();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {
        int grid = pti.index();
        auto& ptile = ParticlesAt(a_lev, pti);
        auto& aos = ptile.GetArrayOfStructs();
        auto& soa = ptile.GetStructOfArrays();
        const int n = aos.numParticles();
        auto* p_pbox = aos().data();

        bool use_terrain = (z_height != nullptr);
        auto zheight = use_terrain ? (*z_height)[grid].array() : Array4<Real>{};

        auto* mass_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::mass).data();

        int rt_offset = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::radius).data();
        auto* mult_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::multiplicity).data();
        auto* supdrop_mass_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::sd_mass).data();

        GpuArray<ParticleReal*,SuperDropletInitializations::num_aerosols_max> aerosol_mass_ptrs;
        for (int i = 0; i < num_aerosols; i++) {
            aerosol_mass_ptrs[i] = soa.GetRealData(   rt_offset
                                                    + SuperDropletsRealIdxSoA_RT::ncomps
                                                    + i ).data();
        }

        const auto& sat_pressure_arr = mf_sat_pressure[grid].array();
        const auto& sat_ratio_arr = mf_sat_ratio[grid].array();
        const auto& temperature_arr = mf_temperature[grid].array();

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
                                                                         1.0e-6,1.0e-99,1.0e-99,10 };

        amrex::Gpu::Buffer<amrex::Long> unconverged_particles({0});
        amrex::Long* unconverged_particles_ptr = unconverged_particles.data();

        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
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
            for (int i = 0; i < num_aerosols; i++) { solute_mass += aerosol_mass_ptrs[i][n]; }

            ParticleReal r_sq = radius_ptr[n]*radius_ptr[n];
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
                radius_ptr[n] = std::sqrt(r_sq);

                // update mass of particle
                mass_ptr[n] = solute_mass + (4.0/3.0)*PI*r_sq*radius_ptr[n]*mat_density;

                // update superdroplet total mass
                supdrop_mass_ptr[n] = mass_ptr[n] * mult_ptr[n];
            }

        });

        //auto const num_unconverged_particles = *(unconverged_particles.copyToHost());
    }
}

#endif

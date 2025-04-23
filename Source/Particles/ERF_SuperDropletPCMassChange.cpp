#include "ERF_SuperDropletPC.H"
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuBuffer.H>
#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_Constants.H"
#include "ERF_SuperDropletPCMassChange.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;
using namespace SDMassChangeUtils;

/*! Compute mass change of particles due to evaporation and condensation */
void SuperDropletPC::MassChange ( int                                         a_lev,
                                  Real                                        a_dt,
                                  const std::string&                          a_vap_name,
                                  const MultiFab&                             a_temperature,
                                  const MultiFab&                             a_sat_pressure,
                                  const MultiFab&                             a_sat_ratio,
                                  const Vector<std::unique_ptr<MultiFab>>&    a_z_phys_nd,
                                  const bool                                  a_is_water )
{
    BL_PROFILE("SuperDropletPC::MassChange()");
    AMREX_ASSERT( a_lev == m_lev );

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    const auto is_periodic = geom.isPeriodic();
    auto is_periodic_z = is_periodic[2];

    const std::unique_ptr<MultiFab>& z_height = a_z_phys_nd[a_lev];

    const int num_species  = m_num_species;
    const int num_aerosols = m_num_aerosols;

    int idx_vap = -1;
    Real mat_density = -1;
    for (int i = 0; i < m_num_species; i++) {
        if (m_species_mat[i]->name() == a_vap_name) {
            idx_vap = i;
            mat_density = m_species_mat[i]->density();
        }
    }
    AMREX_ALWAYS_ASSERT(idx_vap >= 0);
    AMREX_ALWAYS_ASSERT(mat_density >= 0);

    const bool log_unconverged = m_mass_change_logging;
    [[maybe_unused]] FILE* file_handle = m_mass_change_log;

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

        auto zheight = (*z_height)[grid].array();

        auto* mass_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::mass).data();

        int rt_offset = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::radius).data();
        auto* mult_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::multiplicity).data();
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        auto* condt_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::cond_tendency).data();
#endif

        SDSpeciesMassArr sp_mass_ptrs;
        Gpu::DeviceVector<ParticleReal> sp_mol_weight(num_species);
        Gpu::DeviceVector<int> sp_solubility(num_species);
        {
            Vector<ParticleReal> sp_mol_weight_h(num_species);
            Vector<int> sp_solubility_h(num_species);
            for (int i = 0; i < num_species; i++) {
                sp_mass_ptrs[i] = soa.GetRealData(idx_s(i,num_aerosols,num_species)).data();
                sp_mol_weight_h[i] = m_species_mat[i]->molWeight();
                sp_solubility_h[i] = static_cast<int>(m_species_mat[i]->isSoluble());
            }
            Gpu::copy(  Gpu::hostToDevice,
                        sp_mol_weight_h.begin(),
                        sp_mol_weight_h.end(),
                        sp_mol_weight.begin() );
            Gpu::copy(  Gpu::hostToDevice,
                        sp_solubility_h.begin(),
                        sp_solubility_h.end(),
                        sp_solubility.begin() );
        }

        SDAerosolMassArr ae_mass_ptrs;
        Gpu::DeviceVector<ParticleReal> ae_mol_weight(num_aerosols);
        Gpu::DeviceVector<int> ae_solubility(num_aerosols);
        {
            Vector<ParticleReal> ae_mol_weight_h(num_aerosols);
            Vector<int> ae_solubility_h(num_aerosols);
            for (int i = 0; i < num_aerosols; i++) {
                ae_mass_ptrs[i] = soa.GetRealData(idx_a(i,num_aerosols,num_species)).data();
                ae_mol_weight_h[i] = m_aerosol_mat[i]->molWeight();
                ae_solubility_h[i] = static_cast<int>(m_aerosol_mat[i]->isSoluble());
            }
            Gpu::copy(  Gpu::hostToDevice,
                        ae_mol_weight_h.begin(),
                        ae_mol_weight_h.end(),
                        ae_mol_weight.begin() );
            Gpu::copy(  Gpu::hostToDevice,
                        ae_solubility_h.begin(),
                        ae_solubility_h.end(),
                        ae_solubility.begin() );
        }

        const auto& sat_pressure_arr = a_sat_pressure[grid].array();
        const auto& sat_ratio_arr = a_sat_ratio[grid].array();
        const auto& temperature_arr = a_temperature[grid].array();

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        dRdt<ParticleReal> drdt{ m_species_mat[idx_vap]->coeffCurv(),
                                 m_species_mat[idx_vap]->coeffVPSolute(),
                                 m_species_mat[idx_vap]->latHeatVap(),
                                 therco, /* ERF_Constants.H */
                                 m_species_mat[idx_vap]->Rv(),
                                 mat_density,
                                 diffelq /* ERF_Constants.H */};
#endif

        dRsqdt<ParticleReal> drsqdt{ m_species_mat[idx_vap]->coeffCurv(),
                                     m_species_mat[idx_vap]->coeffVPSolute(),
                                     m_species_mat[idx_vap]->latHeatVap(),
                                     therco, /* ERF_Constants.H */
                                     m_species_mat[idx_vap]->Rv(),
                                     mat_density,
                                     diffelq /* ERF_Constants.H */};

        NewtonSolver< dRsqdt<ParticleReal>, ParticleReal > newton_solver{ drsqdt,
                                                                          m_newton_rtol,
                                                                          m_newton_atol,
                                                                          m_newton_stol,
                                                                          m_newton_maxits,
                                                                          false };

        Gpu::Buffer<Long> unconverged_particles({0});
        auto* unconverged_particles_ptr = unconverged_particles.data();

        auto cfl = m_mass_change_cfl;
        auto ti_choice = m_mass_change_ti;

        auto sp_mw_arr = sp_mol_weight.data();
        auto sp_sol_arr = sp_solubility.data();
        auto aero_mw_arr = ae_mol_weight.data();
        auto aero_sol_arr = ae_solubility.data();

        ParallelFor(num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (mult_ptr[i] == 0.0) { return; }

            ParticleReal sat_ratio, e_sat, temperature;
            if (is_periodic_z) {
                cic_interpolate( p, plo, dxi, sat_pressure_arr, &e_sat, 1 );
                cic_interpolate( p, plo, dxi, sat_ratio_arr, &sat_ratio, 1 );
                cic_interpolate( p, plo, dxi, temperature_arr, &temperature, 1 );
            } else {
                cic_interpolate_mapped_z( p, plo, dxi, sat_pressure_arr, zheight, &e_sat, 1 );
                cic_interpolate_mapped_z( p, plo, dxi, sat_ratio_arr, zheight, &sat_ratio, 1 );
                cic_interpolate_mapped_z( p, plo, dxi, temperature_arr, zheight, &temperature, 1 );
            }

            ParticleReal solute_moles = 0.0;
            if (a_is_water) {
                for (int j = 0; j < num_species; j++) {
                    if (sp_sol_arr[j] && (j != idx_vap)) {
                        solute_moles += (sp_mass_ptrs[j][i]/sp_mw_arr[j]);
                    }
                }
                for (int j = 0; j < num_aerosols; j++) {
                    if (aero_sol_arr[j]) {
                        solute_moles += (ae_mass_ptrs[j][i]/aero_mw_arr[j]);
                    }
                }
            }

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
            if (a_is_water) {
                condt_ptr[i] = drdt( radius_ptr[i], sat_ratio, temperature, e_sat, solute_moles);
            }
#endif

            TI< dRsqdt<ParticleReal>,
                NewtonSolver<dRsqdt<ParticleReal>, ParticleReal>,
                ParticleReal > ti { drsqdt, newton_solver, a_dt, 100,
                                    sat_ratio, temperature, e_sat, solute_moles,
                                    cfl, 1e-40, 1e-3, 1e-6, false, false };

            auto mass = sp_mass_ptrs[idx_vap][i];
            auto radius = std::cbrt(mass / ((4.0/3.0)*PI*mat_density));
            auto r_sq = radius*radius;
            bool success = false;
            if (ti_choice == SDMassChangeTIMethod::RK4) {
                ti.rk4(r_sq, success);
            } else if (ti_choice == SDMassChangeTIMethod::RK3BS) {
                ti.rk3bs(r_sq, success);
            } else if (ti_choice == SDMassChangeTIMethod::BE) {
                ti.be(r_sq, success);
            } else if (ti_choice == SDMassChangeTIMethod::CN) {
                ti.cn(r_sq, success);
            } else if (ti_choice == SDMassChangeTIMethod::DIRK2) {
                ti.dirk212(r_sq, success);
            } else {
                printf("ERROR: invalid time integrator choice!\n");
                return;
            }

            if (!success) {
                if (log_unconverged) {
#ifndef AMREX_USE_GPU
                    fprintf(file_handle,
                            "r=%1.16e, S=%1.16e, T=%1.16e, e=%1.16e, sol_mass=%1.16e\n",
                            radius_ptr[i], sat_ratio, temperature, e_sat, solute_moles );
#endif
                }
                Gpu::Atomic::Add(unconverged_particles_ptr, Long(1));
            } else {
                // update particle attributes
                auto radius = std::sqrt(r_sq);
                sp_mass_ptrs[idx_vap][i] = (4.0/3.0)*PI*radius*r_sq*mat_density;
                if (a_is_water) {
                    radius_ptr[i] = radius;
                    mass_ptr[i] = sp_mass_ptrs[idx_vap][i];
                }
            }

        });
        Gpu::synchronize();
        m_num_unconverged_particles += *(unconverged_particles.copyToHost());
    }

}

#endif

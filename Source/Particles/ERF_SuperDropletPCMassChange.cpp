#include "ERF_SuperDropletPC.H"
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuBuffer.H>
#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_Constants.H"
#include "ERF_SuperDropletPCMassChange.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;
using namespace SDPCDefn;
using namespace SDMassChangeUtils;

/*! Compute mass change of particles due to evaporation and condensation */
void SuperDropletPC::MassChange ( int                                         a_lev,
                                  Real                                        a_dt,
                                  const Species::Name&                        a_vap_name,
                                  const MultiFab&                             a_temperature,
                                  const MultiFab&                             a_pressure,
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

    const int num_sp  = m_num_species;
    const int num_ae = m_num_aerosols;

    const Real rho_w = m_species_mat[m_idx_w]->m_density;
    int idx_w = m_idx_w;

    int idx_vap = -1;
    Real mat_density = -1;
    for (int i = 0; i < m_num_species; i++) {
        if (m_species_mat[i]->m_name == a_vap_name) {
            idx_vap = i;
            mat_density = m_species_mat[i]->m_density;
        }
    }
    AMREX_ALWAYS_ASSERT(idx_vap >= 0);
    AMREX_ALWAYS_ASSERT(mat_density >= 0);
    const MaterialProperties vapour_mat(*(m_species_mat[idx_vap]));
    const MaterialPropertiesCore& vapour_mat_core = vapour_mat;

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

        int rtoff_i = SuperDropletsIntIdxSoA::ncomps;
        auto* active_ptr = soa.GetIntData(rtoff_i+SuperDropletsIntIdxSoA_RT::active).data();
        int rtoff_r = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::radius).data();
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        auto* condt_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::cond_tendency).data();
#endif

        SDSpeciesMassArr sp_mass_ptrs;
        Gpu::DeviceVector<ParticleReal> sp_ionization(num_sp);
        Gpu::DeviceVector<ParticleReal> sp_mol_weight(num_sp);
        Gpu::DeviceVector<ParticleReal> sp_density(num_sp);
        Gpu::DeviceVector<int> sp_solubility(num_sp);
        {
            Vector<ParticleReal> sp_ionization_h(num_sp);
            Vector<ParticleReal> sp_mol_weight_h(num_sp);
            Vector<ParticleReal> sp_density_h(num_sp);
            Vector<int> sp_solubility_h(num_sp);
            for (int i = 0; i < num_sp; i++) {
                sp_mass_ptrs[i] = soa.GetRealData(idx_s(i,num_ae,num_sp)).data();
                sp_ionization_h[i] = m_species_mat[i]->m_ionization;
                sp_mol_weight_h[i] = m_species_mat[i]->m_mol_weight;
                sp_density_h[i] = m_species_mat[i]->m_density;
                sp_solubility_h[i] = static_cast<int>(m_species_mat[i]->m_is_soluble);
            }
            Gpu::copy(  Gpu::hostToDevice,
                        sp_ionization_h.begin(),
                        sp_ionization_h.end(),
                        sp_ionization.begin() );
            Gpu::copy(  Gpu::hostToDevice,
                        sp_mol_weight_h.begin(),
                        sp_mol_weight_h.end(),
                        sp_mol_weight.begin() );
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
        Gpu::DeviceVector<ParticleReal> ae_ionization(num_ae);
        Gpu::DeviceVector<ParticleReal> ae_mol_weight(num_ae);
        Gpu::DeviceVector<ParticleReal> ae_density(num_ae);
        Gpu::DeviceVector<int> ae_solubility(num_ae);
        {
            Vector<ParticleReal> ae_ionization_h(num_ae);
            Vector<ParticleReal> ae_mol_weight_h(num_ae);
            Vector<ParticleReal> ae_density_h(num_ae);
            Vector<int> ae_solubility_h(num_ae);
            for (int i = 0; i < num_ae; i++) {
                ae_mass_ptrs[i] = soa.GetRealData(idx_a(i,num_ae,num_sp)).data();
                ae_ionization_h[i] = m_aerosol_mat[i]->m_ionization;
                ae_mol_weight_h[i] = m_aerosol_mat[i]->m_mol_weight;
                ae_density_h[i] = m_aerosol_mat[i]->m_density;
                ae_solubility_h[i] = static_cast<int>(m_aerosol_mat[i]->m_is_soluble);
            }
            Gpu::copy(  Gpu::hostToDevice,
                        ae_ionization_h.begin(),
                        ae_ionization_h.end(),
                        ae_ionization.begin() );
            Gpu::copy(  Gpu::hostToDevice,
                        ae_mol_weight_h.begin(),
                        ae_mol_weight_h.end(),
                        ae_mol_weight.begin() );
            Gpu::copy(  Gpu::hostToDevice,
                        ae_density_h.begin(),
                        ae_density_h.end(),
                        ae_density.begin() );
            Gpu::copy(  Gpu::hostToDevice,
                        ae_solubility_h.begin(),
                        ae_solubility_h.end(),
                        ae_solubility.begin() );
        }

        const auto& sat_pressure_arr = a_sat_pressure[grid].array();
        const auto& sat_ratio_arr = a_sat_ratio[grid].array();
        const auto& temperature_arr = a_temperature[grid].array();
        const auto& pressure_arr = a_pressure[grid].array();

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        dRdt<ParticleReal> drdt{ vapour_mat.m_lat_vap,
                                 therco, /* ERF_Constants.H */
                                 vapour_mat.m_Rv,
                                 mat_density };
#endif

        dRsqdt<ParticleReal> drsqdt{ vapour_mat.m_lat_vap,
                                     therco, /* ERF_Constants.H */
                                     vapour_mat.m_Rv,
                                     mat_density };

        NewtonSolver< dRsqdt<ParticleReal>, ParticleReal > newton_solver{ drsqdt,
                                                                          m_newton_rtol,
                                                                          m_newton_atol,
                                                                          m_newton_stol,
                                                                          m_newton_maxits };

        Gpu::Buffer<Long> unconverged_particles({0});
        auto* unconverged_particles_ptr = unconverged_particles.data();

        auto cfl = m_mass_change_cfl;
        auto ti_choice = m_mass_change_ti;
        AMREX_ASSERT_WITH_MESSAGE( ti_choice == SDMassChangeTIMethod::RK4 ||
                                   ti_choice == SDMassChangeTIMethod::RK3BS ||
                                   ti_choice == SDMassChangeTIMethod::BE ||
                                   ti_choice == SDMassChangeTIMethod::CN ||
                                   ti_choice == SDMassChangeTIMethod::DIRK2,
                                   "ERROR: invalid time integrator choice!" );

        auto sp_i_arr = sp_ionization.data();
        auto sp_mw_arr = sp_mol_weight.data();
        auto sp_rho_arr = sp_density.data();
        auto sp_sol_arr = sp_solubility.data();
        auto ae_i_arr = ae_ionization.data();
        auto ae_mw_arr = ae_mol_weight.data();
        auto ae_rho_arr = ae_density.data();
        auto ae_sol_arr = ae_solubility.data();

        ParallelFor(num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (active_ptr[i] == 0) { return; }

            ParticleReal sat_ratio, e_sat, temperature, pressure;
            if (is_periodic_z) {
                cic_interpolate( p, plo, dxi, sat_pressure_arr, &e_sat, 1 );
                cic_interpolate( p, plo, dxi, sat_ratio_arr, &sat_ratio, 1 );
                cic_interpolate( p, plo, dxi, temperature_arr, &temperature, 1 );
                cic_interpolate( p, plo, dxi, pressure_arr, &pressure, 1 );
            } else {
                cic_interpolate_mapped_z( p, plo, dxi, sat_pressure_arr, zheight, &e_sat, 1 );
                cic_interpolate_mapped_z( p, plo, dxi, sat_ratio_arr, zheight, &sat_ratio, 1 );
                cic_interpolate_mapped_z( p, plo, dxi, temperature_arr, zheight, &temperature, 1 );
                cic_interpolate_mapped_z( p, plo, dxi, pressure_arr, zheight, &pressure, 1 );
            }

            ParticleReal solute_moles = 0.0;
            if (a_is_water) {
                for (int j = 0; j < num_sp; j++) {
                    if (j != idx_vap) {
                        solute_moles += (sp_mass_ptrs[j][i]*sp_i_arr[j]/sp_mw_arr[j]);
                    }
                }
                for (int j = 0; j < num_ae; j++) {
                    solute_moles += (ae_mass_ptrs[j][i]*ae_i_arr[j]/ae_mw_arr[j]);
                }
            }

            auto coeff_curv = vapour_mat_core.coeffCurv(temperature);
            auto coeff_sol = vapour_mat_core.coeffVPSolute();
            auto coeff_moldiff = vapour_mat_core.coeffMolecularDiffusion(temperature, pressure);

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
            if (a_is_water) {
                condt_ptr[i] = drdt( radius_ptr[i],
                                     sat_ratio,
                                     temperature,
                                     e_sat,
                                     coeff_moldiff,
                                     coeff_curv,
                                     coeff_sol,
                                     solute_moles);
            }
#endif

            TI< dRsqdt<ParticleReal>,
                NewtonSolver<dRsqdt<ParticleReal>, ParticleReal>,
                ParticleReal > ti { drsqdt, newton_solver, a_dt, 100,
                                    sat_ratio, temperature, e_sat,
                                    coeff_moldiff, coeff_curv, coeff_sol,
                                    solute_moles,
                                    cfl, 1e-40, 1e-3, 1e-6, false, false };

            auto r_init = SD_effective_radius( i, idx_w,
                                               rho_w,
                                               num_sp, num_ae,
                                               sp_sol_arr, ae_sol_arr,
                                               sp_mass_ptrs, ae_mass_ptrs,
                                               sp_rho_arr, ae_rho_arr );

            auto r_sq = r_init*r_init;
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
                auto r_new = std::sqrt(r_sq);
                auto d_mass = (4.0/3.0)*PI*mat_density * (r_new*r_new*r_new - r_init*r_init*r_init);
                sp_mass_ptrs[idx_vap][i] += d_mass;
                // don't let it go negative
                sp_mass_ptrs[idx_vap][i] = std::max(sp_mass_ptrs[idx_vap][i],0.0);

                radius_ptr[i] = SD_effective_radius( i, idx_w,
                                                     rho_w,
                                                     num_sp, num_ae,
                                                     sp_sol_arr, ae_sol_arr,
                                                     sp_mass_ptrs, ae_mass_ptrs,
                                                     sp_rho_arr, ae_rho_arr );
                mass_ptr[i] = SD_total_mass( i, num_sp, num_ae, sp_mass_ptrs, ae_mass_ptrs);
            }

        });
        Gpu::synchronize();
        m_num_unconverged_particles += *(unconverged_particles.copyToHost());
    }

}

#endif

#include "ERF_SuperDropletPC.H"
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuBuffer.H>
#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_Constants.H"
#include "ERF_SuperDropletPCMassChange.H"
#include "ERF_InterpolationUtils.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;
using namespace SDPCDefn;

/*! Compute mass change of particles due to evaporation and condensation
 *  (liquid <--> vapour) */
void SuperDropletPC::MassChange_LV (  int                                         a_lev,
                                      Real                                        a_dt,
                                      const Species::Name&                        a_vap_name,
                                      const MultiFab&                             a_temperature,
                                      const MultiFab&                             a_pressure,
                                      const MultiFab&                             a_sat_pressure,
                                      const MultiFab&                             a_sat_ratio,
                                      const Vector<std::unique_ptr<MultiFab>>&    a_z_phys_nd,
                                      const bool                                  a_is_water )
{
    using namespace SDMassChangeUtils_LV;

    BL_PROFILE("SuperDropletPC::MassChange_LV()");
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
    const int idx_w = m_idx_w;
    const int idx_i = m_idx_i;

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
        SDAerosolMassArr ae_mass_ptrs;
        setupMassPointers(soa, sp_mass_ptrs, ae_mass_ptrs);

        // Get pointers to persistent device data
        const ParticleReal* sp_rho_arr = getSpeciesDensitiesDevice();
        const int* sp_sol_arr = getSpeciesSolubilitiesDevice();
        const ParticleReal* sp_ion_arr = getSpeciesIonizationDevice();
        const ParticleReal* sp_mw_arr = getSpeciesMolWeightDevice();

        const ParticleReal* ae_rho_arr = getAerosolDensitiesDevice();
        const int* ae_sol_arr = getAerosolSolubilitiesDevice();
        const ParticleReal* ae_ion_arr = getAerosolIonizationDevice();
        const ParticleReal* ae_mw_arr = getAerosolMolWeightDevice();

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

        // We already have these pointers from our accessors above
        // sp_ion_arr contains species ionization values
        // sp_mw_arr contains species molecular weights
        // sp_rho_arr contains species densities
        // sp_sol_arr contains species solubility flags
        // ae_ion_arr contains aerosol ionization values
        // ae_mw_arr contains aerosol molecular weights
        // ae_rho_arr contains aerosol densities
        // ae_sol_arr contains aerosol solubility flags

        ParallelFor(num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (active_ptr[i] == 0) { return; }

            // skip ice particles
            auto par_phase = SD_phase(i, idx_w, idx_i, sp_mass_ptrs);
            if (a_is_water && (par_phase == SDPhase::ice)) { return; }

            // Define field values array to store interpolated results
            ParticleReal field_values[4]; // e_sat, sat_ratio, temperature, pressure

            // Define array of field arrays to interpolate from
            const Array4<const Real> field_arrays[4] = {
                sat_pressure_arr,
                sat_ratio_arr,
                temperature_arr,
                pressure_arr
            };

            // Use the interpolation helper function to interpolate all fields at once
            ERF::Interpolation::interpolateFields(
                p, plo, dxi, field_arrays, field_values, 4,
                is_periodic_z ? 1 : 0, is_periodic_z ? nullptr : &zheight
            );

            // Extract interpolated values
            ParticleReal e_sat = field_values[0];
            ParticleReal sat_ratio = field_values[1];
            ParticleReal temperature = field_values[2];
            ParticleReal pressure = field_values[3];

            ParticleReal solute_moles = 0.0;
            if (a_is_water) {
                for (int j = 0; j < num_sp; j++) {
                    if (j != idx_vap) {
                        solute_moles += (sp_mass_ptrs[j][i]*sp_ion_arr[j]/sp_mw_arr[j]);
                    }
                }
                for (int j = 0; j < num_ae; j++) {
                    solute_moles += (ae_mass_ptrs[j][i]*ae_ion_arr[j]/ae_mw_arr[j]);
                }
            }

            auto coeff_curv = vapour_mat.coeffCurv(temperature);
            auto coeff_sol = vapour_mat.coeffVPSolute();
            auto coeff_moldiff = vapour_mat.coeffMolecularDiffusion(temperature, pressure);

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
                auto r_new = std::sqrt(r_sq);
                auto d_mass = (4.0/3.0)*PI*mat_density * (r_new*r_new*r_new - r_init*r_init*r_init);
                sp_mass_ptrs[idx_vap][i] += d_mass;
                // don't let it go negative
                sp_mass_ptrs[idx_vap][i] = std::max(sp_mass_ptrs[idx_vap][i],0.0);

                // Update particle attributes (radius and mass)
                SuperDropletPC::updateParticleAttributes(
                    i, radius_ptr, mass_ptr, idx_w, rho_w,
                    num_sp, num_ae, sp_sol_arr, ae_sol_arr,
                    sp_mass_ptrs, ae_mass_ptrs, sp_rho_arr, ae_rho_arr);
            }

        });
        Gpu::synchronize();
        m_num_unconverged_particles += *(unconverged_particles.copyToHost());
    }

}

/*! Compute mass change of particles due to freezing/melting
 *  (solid <--> liquid) */
void SuperDropletPC::MassChange_SL (  int                                         a_lev,
                                      Real                                        a_dt,
                                      const MultiFab&                             a_temperature,
                                      const MultiFab&                             a_sat_ratio,
                                      const Vector<std::unique_ptr<MultiFab>>&    a_z_phys_nd )
{
    BL_PROFILE("SuperDropletPC::MassChange_SL()");
    AMREX_ASSERT( a_lev == m_lev );
    amrex::ignore_unused(a_dt);

    if (m_idx_i < 0) { return; } // ice not being modeled

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    const auto is_periodic = geom.isPeriodic();
    auto is_periodic_z = is_periodic[2];

    const std::unique_ptr<MultiFab>& z_height = a_z_phys_nd[a_lev];

    const Real rho_w = m_species_mat[m_idx_w]->m_density;
    const Real rho_i = m_species_mat[m_idx_i]->m_density;
    const int idx_w = m_idx_w;
    const int idx_i = m_idx_i;
    AMREX_ALWAYS_ASSERT((idx_w >= 0) && (idx_i >= 0) && (idx_w != idx_i));

    const int num_sp  = m_num_species;
    const int num_ae = m_num_aerosols;

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
        auto* mult_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::multiplicity).data();
        auto* Tfz_ptr = soa.GetRealData(idx_ice_Tfz(num_ae,num_sp)).data();
        auto* a_ptr = soa.GetRealData(idx_ice_a(num_ae,num_sp)).data();
        auto* c_ptr = soa.GetRealData(idx_ice_c(num_ae,num_sp)).data();
        auto* mrime_ptr = soa.GetRealData(idx_ice_mrime(num_ae,num_sp)).data();
        auto* nmono_ptr = soa.GetRealData(idx_ice_nmono(num_ae,num_sp)).data();

        const auto& sat_ratio_arr = a_sat_ratio[grid].array();
        const auto& temperature_arr = a_temperature[grid].array();

        SDSpeciesMassArr sp_mass_ptrs;
        SDAerosolMassArr ae_mass_ptrs;
        setupMassPointers(soa, sp_mass_ptrs, ae_mass_ptrs);

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

        // ae_mass_ptrs already defined above
        Gpu::DeviceVector<ParticleReal> ae_density(num_ae);
        Gpu::DeviceVector<int> ae_solubility(num_ae);
        {
            Vector<ParticleReal> ae_density_h(num_ae);
            Vector<int> ae_solubility_h(num_ae);
            for (int i = 0; i < num_ae; i++) {
                // ae_mass_ptrs already set up via setupMassPointers
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

        ParallelFor(num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (active_ptr[i] == 0) { return; }

            // Define field values array to store interpolated results
            ParticleReal field_values[2]; // temperature, sat_ratio

            // Define array of field arrays to interpolate from
            const Array4<const Real> field_arrays[2] = {
                temperature_arr,
                sat_ratio_arr
            };

            // Use the interpolation helper function to interpolate all fields at once
            ERF::Interpolation::interpolateFields(
                p, plo, dxi, field_arrays, field_values, 2,
                is_periodic_z ? 1 : 0, is_periodic_z ? nullptr : &zheight
            );

            // Extract interpolated values
            ParticleReal temperature = field_values[0];
            ParticleReal sat_ratio = field_values[1];

            auto par_phase = SD_phase(i, idx_w, idx_i, sp_mass_ptrs);
            if (par_phase == SDPhase::water) {
                // SD is water, check for freezing
                if ((temperature <= Tfz_ptr[i]) && (sat_ratio > 1.0)) {
                    sp_mass_ptrs[idx_i][i] = sp_mass_ptrs[idx_w][i];
                    sp_mass_ptrs[idx_w][i] = 0.0;
                    a_ptr[i] = c_ptr[i] = std::cbrt(sp_mass_ptrs[idx_i][i]/(4.0/3.0*PI*rho_i));
                    mrime_ptr[i] = 0.0;
                    nmono_ptr[i] = 1.0;
                }

            } else if (par_phase == SDPhase::ice) {

                // SD is ice, check for melting
                if (temperature > tmelt /* from ERF_Constants.H */) {
                    sp_mass_ptrs[idx_w][i] = sp_mass_ptrs[idx_i][i];
                    sp_mass_ptrs[idx_i][i] = 0.0;
                    a_ptr[i] = c_ptr[i] = mrime_ptr[i] = nmono_ptr[i] = 0.0;
                }

            } else {
                amrex::Abort("Unknown value for particle phase (must be ice or water)");
            }

            // update particle attributes
            SuperDropletPC::updateParticleAttributes(
                i, radius_ptr, mass_ptr, idx_w, rho_w,
                num_sp, num_ae, sp_sol_arr, ae_sol_arr,
                sp_mass_ptrs, ae_mass_ptrs, sp_rho_arr, ae_rho_arr);

        });
        Gpu::synchronize();
    }

}

/*! Compute mass change of particles due to deposition and sublimation
 *  (solid <--> vapour) */
void SuperDropletPC::MassChange_SV (  int                                      a_lev,
                                      Real                                     a_dt,
                                      const MultiFab&                          a_density,
                                      const MultiFab&                          a_temperature,
                                      const MultiFab&                          a_pressure,
                                      const MultiFab&                          a_moist_density,
                                      const MultiFab&                          a_sat_pressure,
                                      const MultiFab&                          a_sat_pressure_wi,
                                      const MultiFab&                          a_sat_ratio,
                                      const Vector<std::unique_ptr<MultiFab>>& a_z_phys_nd )
{
    using namespace SDMassChangeUtils_SV;

    BL_PROFILE("SuperDropletPC::MassChange_SV()");
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
    const int idx_w = m_idx_w;
    const int idx_i = m_idx_i;
    AMREX_ALWAYS_ASSERT((idx_w >= 0) && (idx_i >= 0) && (idx_w != idx_i));
    const auto mat_prop(*(m_species_mat[idx_i]));

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
        auto* vterm_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::term_vel).data();
        auto* mult_ptr = soa.GetRealData(rtoff_r+SuperDropletsRealIdxSoA_RT::multiplicity).data();
        auto* a_ptr = soa.GetRealData(idx_ice_a(num_ae,num_sp)).data();
        auto* c_ptr = soa.GetRealData(idx_ice_c(num_ae,num_sp)).data();
        auto* mrime_ptr = soa.GetRealData(idx_ice_mrime(num_ae,num_sp)).data();

        SDSpeciesMassArr sp_mass_ptrs;
        SDAerosolMassArr ae_mass_ptrs;
        setupMassPointers(soa, sp_mass_ptrs, ae_mass_ptrs);

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

        // ae_mass_ptrs already defined above
        Gpu::DeviceVector<ParticleReal> ae_density(num_ae);
        Gpu::DeviceVector<int> ae_solubility(num_ae);
        {
            Vector<ParticleReal> ae_density_h(num_ae);
            Vector<int> ae_solubility_h(num_ae);
            for (int i = 0; i < num_ae; i++) {
                // ae_mass_ptrs already set up via setupMassPointers
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

        const auto& sat_pressure_arr = a_sat_pressure[grid].array();
        const auto& sat_pressure_wi_arr = a_sat_pressure_wi[grid].array();
        const auto& sat_ratio_arr = a_sat_ratio[grid].array();
        const auto& density_arr = a_density[grid].array();
        const auto& temperature_arr = a_temperature[grid].array();
        const auto& pressure_arr = a_pressure[grid].array();
        const auto& moist_density_arr = a_moist_density[grid].array();

        dMdt<ParticleReal> dmdt{ mat_prop.m_lat_vap,
                                 therco, /* ERF_Constants.H */
                                 mat_prop.m_Rv,
                                 mat_prop.m_density };

        auto sp_rho_arr = sp_density.data();
        auto sp_sol_arr = sp_solubility.data();
        auto ae_rho_arr = ae_density.data();
        auto ae_sol_arr = ae_solubility.data();

        ParallelFor(num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (active_ptr[i] == 0) { return; }

            // skip water particles
            auto par_phase = SD_phase(i, idx_w, idx_i, sp_mass_ptrs);
            if (par_phase == SDPhase::water) { return; }

            // Define field values array to store interpolated results
            ParticleReal field_values[7]; // e_sat, e_sat_ratio_wi, sat_ratio, density, temperature, pressure, moist_density

            // Define array of field arrays to interpolate from
            const Array4<const Real> field_arrays[7] = {
                sat_pressure_arr,
                sat_pressure_wi_arr,
                sat_ratio_arr,
                density_arr,
                temperature_arr,
                pressure_arr,
                moist_density_arr
            };

            // Use the interpolation helper function to interpolate all fields at once
            ERF::Interpolation::interpolateFields(
                p, plo, dxi, field_arrays, field_values, 7,
                is_periodic_z ? 1 : 0, is_periodic_z ? nullptr : &zheight
            );

            // Extract interpolated values
            ParticleReal e_sat = field_values[0];
            ParticleReal e_sat_ratio_wi = field_values[1];
            ParticleReal sat_ratio = field_values[2];
            // Note: density (field_values[3]) is unused in this function
            ParticleReal temperature = field_values[4];
            ParticleReal pressure = field_values[5];
            ParticleReal moist_density = field_values[6];

            auto coeff_moldiff = mat_prop.coeffMolecularDiffusion(temperature, pressure);

            TI< dMdt<ParticleReal>,
                ParticleReal > ti { dmdt, a_dt, a_dt, 100,
                                    a_ptr[i], c_ptr[i], vterm_ptr[i],
                                    sat_ratio, moist_density, temperature, e_sat,
                                    coeff_moldiff,
                                    false };

            auto mass_old = sp_mass_ptrs[idx_i][i];
            auto vol_old = (4.0/3.0) * PI * a_ptr[i]*a_ptr[i]*c_ptr[i];
            auto rhoi_old = mass_old / vol_old;

            // compute new mass
            bool success = false;
            auto mass_new = mass_old;
            ti.fe(mass_new, success);
            AMREX_ALWAYS_ASSERT(success);
            mass_new = std::max(mass_new, 3.8403e-24); // lower limit: 1nm spherical ice particle
            auto d_mass = mass_new - mass_old;

            // compute new volume
            auto d_vol = dmdt.dVolume(  d_mass, a_ptr[i], c_ptr[i], rhoi_old,
                                        sat_ratio, temperature,  e_sat, e_sat_ratio_wi,
                                        coeff_moldiff );
            auto vol_new = vol_old + d_vol;
            vol_new = std::max(vol_new, mass_new/mat_prop.m_density);
            d_vol = vol_new - vol_old;
            auto rhoi_new = mass_new / vol_new;

            // compute growth ratio and new radii
            auto gr_star = dmdt.growthRatioStar( d_mass,
                                                 a_ptr[i], c_ptr[i], vterm_ptr[i],
                                                 moist_density, temperature, coeff_moldiff );
            auto d_loga = dmdt.dLogRadius(gr_star, std::log(vol_new/vol_old));
            auto a_new = a_ptr[i] * std::exp(d_loga);
            auto c_new = c_ptr[i] * std::exp(gr_star*d_loga);

            // Smaller than 1 um --> sphere with true ice density
            if (std::min(a_new,c_new) <= 1.0e-6) {
                rhoi_new = mat_prop.m_density;
                vol_new = mass_new / rhoi_new;
                c_new = a_new = std::cbrt(vol_new/(4.0*PI/3.0));
            }
            // limit particle size to 1 nm
            if (std::min(a_new,c_new) <= 1.0e-9 ) {
                c_new = a_new = 1.0e-9;
                rhoi_new = mat_prop.m_density;
                vol_new = (4.0/3.0)*PI*a_new*a_new*c_new;
                mass_new = vol_new * rhoi_new;
            }

            // rime mass
            auto mrime_new = mrime_ptr[i];
            if (d_mass <= 0.0) { mrime_new = mrime_ptr[i] * mass_new/mass_old; }

            // update attributes
            a_ptr[i] = a_new;
            c_ptr[i] = c_new;
            mrime_ptr[i] = mrime_new;
            sp_mass_ptrs[idx_i][i] = mass_new;

            // update particle attributes
            SuperDropletPC::updateParticleAttributes(
                i, radius_ptr, mass_ptr, idx_w, rho_w,
                num_sp, num_ae, sp_sol_arr, ae_sol_arr,
                sp_mass_ptrs, ae_mass_ptrs, sp_rho_arr, ae_rho_arr);

        });
        Gpu::synchronize();
    }

}

#endif

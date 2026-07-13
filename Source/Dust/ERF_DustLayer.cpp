/**
 * @file ERF_DustLayer.cpp
 * @brief Implementation of the DustLayer container class.
 */

#include <ERF_DustLayer.H>
#include <ERF_DustPrerequisites.H>
#include <ERF_DustGrid.H>
#include <ERF_DustSurfaceReader.H>
#include <ERF_PhreeqcReader.H>
#include <ERF.H>
#include <ERF_SurfaceLayer.H>
#include <AMReX_Print.H>
#include <cmath>

void DustLayer::initialize(const ERF&          erf,
                           const SurfaceLayer* surface_layer,
                           const DustParams&   dust_params)
{
    // Step 1: Call verify_dust_prerequisites
    verify_dust_prerequisites(erf, surface_layer, dust_params);

    // Step 2: Create dust grid
    m_dg = create_dust_grid(erf.boxArray(0), erf.DistributionMap(0),
                            erf.Geom(0), dust_params.grid_ratio);

    // Step 3: Store dust_params
    m_params = dust_params;

    // Step 4: Allocate each MultiFab with 1 ghost cell (except z which has 0)
    amrex::IntVect ng(1, 1, 0);

    dust_ustar_t = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
    dust_soil_type = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
    dust_silt_fraction = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
    dust_crust_index = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
    dust_moisture_flag = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
    dust_suppression = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, ng);
    dust_emission_flux = std::make_unique<amrex::MultiFab>(
        m_dg.ba, m_dg.dm, dust_params.n_size_bins, ng);

    // Step 5: Fill initial values using setVal

    // dust_ustar_t: Bagnold threshold for the coarsest bin (bin 0)
    // u*_t = A * sqrt(rho_p * g * d / rho_a)
    // Bagnold (1941). The Physics of Blown Sand and Desert Dunes, Methuen, London.
    amrex::Real g = 9.81;  // m/s^2
    amrex::Real rho_a = 1.225;  // kg/m^3 (standard air density)
    amrex::Real d_bin0 = dust_params.bin_diameter_um[0] * 1.0e-6;  // convert um to m
    amrex::Real ustar_t = dust_params.threshold_A_coeff
                        * std::sqrt(dust_params.particle_density * g * d_bin0 / rho_a);
    dust_ustar_t->setVal(ustar_t);

    // dust_soil_type: 0.0 (undefined; set by surface reader in Phase 3)
    dust_soil_type->setVal(0.0);

    // dust_silt_fraction: silt_fraction from params
    dust_silt_fraction->setVal(dust_params.silt_fraction);

    // dust_crust_index: crust_index from params
    dust_crust_index->setVal(dust_params.crust_index);

    // dust_moisture_flag: 0.0 (dry surface at initialization)
    dust_moisture_flag->setVal(0.0);

    // dust_suppression: 0.0 (no suppression at initialization)
    dust_suppression->setVal(0.0);

    // dust_emission_flux: 0.0 (zero emission at initialization)
    dust_emission_flux->setVal(0.0);

    // Allocate efflorescence and base u*_t MultiFabs.
    dust_efflor     = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
    dust_ustar_base = std::make_unique<amrex::MultiFab>(m_dg.ba, m_dg.dm, 1, amrex::IntVect(1,1,0));
    dust_efflor->setVal(0.0);
    // Store the Bagnold u*_t computed from DustParams as the base value.
    // update_ustar_t_from_chemistry modifies dust_ustar_t; dust_ustar_base is read-only.
    dust_ustar_base->Copy(*dust_ustar_t, 0, 0, 1, amrex::IntVect(1,1,0));

    // Populate surface MultiFabs from external rasters if paths are given in dust_params.
    // MultiFabs retain their setVal defaults above if paths are empty.
    populate_dust_surface_maps(*dust_soil_type, *dust_silt_fraction,
                               *dust_crust_index, *dust_moisture_flag,
                               *dust_suppression, m_dg, dust_params);

    // Step 6: Print status message
    amrex::Box dust_domain = m_dg.ba.minimalBox();
    int dust_nx = dust_domain.length(0);
    int dust_ny = dust_domain.length(1);
    amrex::Print() << "[DUST] DustLayer initialized: grid_ratio="
                   << m_dg.grid_ratio << ", dust cells="
                   << dust_nx << "x" << dust_ny << "x1\n";
}

void DustLayer::advance(amrex::Real     dt,
                        const DustParams& dust_params)
{
    ++m_step;
    m_time += dt;

    // Physics inserted in Phases 5-13. PHREEQC reader called here (Phase 4).
    // Call PHREEQC reader if update interval has elapsed.
    // The interval is set by dust_params.phreeqc_update_interval_s.
    // File-based coupling is appropriate because geochemical processes
    // evolve on timescales of days to weeks, much longer than the
    // atmospheric timestep.
    bool do_phreeqc = (m_last_phreeqc_update < 0.0) ||
                      (m_time - m_last_phreeqc_update >=
                       dust_params.phreeqc_update_interval_s);

    if (do_phreeqc && !dust_params.phreeqc_output_file.empty()) {
        update_dust_from_phreeqc(*dust_ustar_t,
                                 *dust_ustar_base,
                                 *dust_crust_index,
                                 *dust_silt_fraction,
                                 *dust_efflor,
                                 *dust_suppression,
                                 *dust_emission_flux,
                                 m_dg,
                                 dust_params);
        m_last_phreeqc_update = m_time;
    }
}
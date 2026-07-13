/**
 * @file ERF_DustLayer.cpp
 * @brief Implementation of the DustLayer container class.
 */

#include <ERF_DustLayer.H>
#include <ERF_DustPrerequisites.H>
#include <ERF_DustGrid.H>
#include <ERF.H>
#include <SurfaceLayer.H>
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

    // Step 6: Print status message
    amrex::Print() << "[DUST] DustLayer initialized: grid_ratio="
                   << m_dg.grid_ratio << ", dust cells="
                   << m_dg.ba.d_numPts() / m_dg.ba.size() << "x"
                   << m_dg.ba.d_numPts() / m_dg.ba.size() << "x1\n";
}

void DustLayer::advance(amrex::Real     dt,
                        const DustParams& dust_params)
{
    ++m_step;
    // Physics inserted in Phases 5-13. Stub returns without computation.
}

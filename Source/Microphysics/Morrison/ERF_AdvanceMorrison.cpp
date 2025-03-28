#include <string>
#include <vector>
#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_TableData.H>
#include <AMReX_MultiFabUtil.H>

#include "ERF_Constants.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_IndexDefines.H"
#include "ERF_DataStruct.H"
#include "ERF_NullMoist.H"
#include "ERF_Morrison.H"
#include <ERF_Morrison_Fortran_Interface.H>

    // wrapper to do all the updating
    void
    Morrison::Advance (const amrex::Real& dt_advance,
                       const SolverChoice& sc)
    {
        dt = dt_advance;
#if 0
        this->Cloud(sc);
        this->IceFall(sc);
        this->Precip(sc);
        this->PrecipFall(sc);
#else
        // Store timestep
        dt = dt_advance;

        // Loop through the grids
        for (amrex::MFIter mfi(*mic_fab_vars[MicVar_Morr::qcl],TileNoZ()); mfi.isValid(); ++mfi)
        {
          const amrex::Box& box = mfi.tilebox();

          // Get array data from class member variables
          auto const& theta_arr = mic_fab_vars[MicVar_Morr::theta]->array(mfi);
          auto const& qv_arr = mic_fab_vars[MicVar_Morr::qv]->array(mfi);
          auto const& qcl_arr = mic_fab_vars[MicVar_Morr::qcl]->array(mfi);
          auto const& qpr_arr = mic_fab_vars[MicVar_Morr::qpr]->array(mfi);
          auto const& qci_arr = mic_fab_vars[MicVar_Morr::qci]->array(mfi);
          auto const& qps_arr = mic_fab_vars[MicVar_Morr::qps]->array(mfi);
          auto const& qpg_arr = mic_fab_vars[MicVar_Morr::qpg]->array(mfi);
          auto const& ni_arr = mic_fab_vars[MicVar_Morr::ni]->array(mfi);
          auto const& ns_arr = mic_fab_vars[MicVar_Morr::ns]->array(mfi);
          auto const& nr_arr = mic_fab_vars[MicVar_Morr::nr]->array(mfi);
          auto const& ng_arr = mic_fab_vars[MicVar_Morr::ng]->array(mfi);
          auto const& rho_arr = mic_fab_vars[MicVar_Morr::rho]->array(mfi);
          auto const& pres_arr = mic_fab_vars[MicVar_Morr::pres]->array(mfi);
          auto const& tabs_arr = mic_fab_vars[MicVar_Morr::tabs]->array(mfi);
          auto const& rain_accum_arr = mic_fab_vars[MicVar_Morr::rain_accum]->array(mfi);
          auto const& snow_accum_arr = mic_fab_vars[MicVar_Morr::snow_accum]->array(mfi);
          auto const& graup_accum_arr = mic_fab_vars[MicVar_Morr::graup_accum]->array(mfi);
          auto const& w_arr = mic_fab_vars[MicVar_Morr::omega]->array(mfi);

          // Get radar reflectivity array if radar diagnostics enabled
          //        auto const& refl_arr = m_do_radar_ref ? m_radar->array(mfi) : nullptr;
          //        auto const& refl_arr = m_radar->array(mfi);

          // Extract box dimensions
          const int ilo = box.loVect()[0];
          const int ihi = box.hiVect()[0];
          const int jlo = box.loVect()[1];
          const int jhi = box.hiVect()[1];
          const int klo = box.loVect()[2];
          const int khi = box.hiVect()[2];

          amrex::Box grown_box(box); grown_box.grow(3);

          const int ilom = grown_box.loVect()[0];
          const int ihim = grown_box.hiVect()[0];
          const int jlom = grown_box.loVect()[1];
          const int jhim = grown_box.hiVect()[1];
          const int klom = grown_box.loVect()[2];
          const int khim = grown_box.hiVect()[2];

          // Calculate Exner function (PII) to convert potential temperature to temperature
          // PII = (P/P0)^(R/cp)
          amrex::FArrayBox pii_fab(grown_box, 1);
          auto const& pii_arr = pii_fab.array();

          const amrex::Real p0 = 100000.0; // Reference pressure (Pa)

          const amrex::Real rdcp = m_rdOcp; // R/cp ratio

          // Calculate Exner function
          amrex::ParallelFor(grown_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // NOTE: the Morrison Fortran version uses Pa not hPa so we didn't divide p by 100
            //       so we don't need to multiply by 100 here
            pii_arr(i,j,k) = std::pow((pres_arr(i,j,k)) / p0, rdcp);
          });

          // Create arrays for height differences (dz)
          amrex::FArrayBox dz_fab(grown_box, 1);
          auto const& dz_arr = dz_fab.array();

          // Calculate height differences
          const amrex::Real dz_val = m_geom.CellSize(m_axis);
          amrex::ParallelFor(grown_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            dz_arr(i,j,k) = dz_val;
          });
          amrex::Box grown_boxD(grown_box); grown_boxD.makeSlab(2,0);
          // Arrays to store precipitation rates
          amrex::FArrayBox    rainncv_fab(grown_boxD, 1);
          amrex::FArrayBox         sr_fab(grown_boxD, 1);     // Ratio of snow to total precipitation
          amrex::FArrayBox    snowncv_fab(grown_boxD, 1);
          amrex::FArrayBox graupelncv_fab(grown_boxD, 1);

          auto const& rainncv_arr = rainncv_fab.array();
          auto const& sr_arr = sr_fab.array();
          auto const& snowncv_arr = snowncv_fab.array();
          auto const& graupelncv_arr = graupelncv_fab.array();

          // Initialize precipitation rate arrays to zero
          rainncv_fab.setVal(0.0);
          sr_fab.setVal(0.0);
          snowncv_fab.setVal(0.0);
          graupelncv_fab.setVal(0.0);

          // Create terrain height array (not actually used by Morrison scheme)
          amrex::FArrayBox ht_fab(amrex::Box(amrex::IntVect(ilo, jlo, 0), amrex::IntVect(ihi, jhi, 0)), 1);
          auto const& ht_arr = ht_fab.array();
          ht_fab.setVal(0.0);  // Not used by Morrison scheme

          // Create dummy arrays for cumulus tendencies (if needed)
          amrex::FArrayBox qrcuten_fab(grown_box, 1);
          amrex::FArrayBox qscuten_fab(grown_box, 1);
          amrex::FArrayBox qicuten_fab(grown_box, 1);
          auto const& qrcuten_arr = qrcuten_fab.array();
          auto const& qscuten_arr = qscuten_fab.array();
          auto const& qicuten_arr = qicuten_fab.array();

          // Initialize tendencies to zero (no cumulus parameterization in this example)
          qrcuten_fab.setVal(0.0);
          qscuten_fab.setVal(0.0);
          qicuten_fab.setVal(0.0);

          // WRF-Chem related variables (optional)
          bool flag_qndrop = false;  // Flag to indicate droplet number prediction

          // Now create arrays for other optional variables
          amrex::FArrayBox rainprod_fab(grown_box, 1);
          amrex::FArrayBox evapprod_fab(grown_box, 1);
          amrex::FArrayBox qlsink_fab(grown_box, 1);
          amrex::FArrayBox precr_fab(grown_box, 1);
          amrex::FArrayBox preci_fab(grown_box, 1);
          amrex::FArrayBox precs_fab(grown_box, 1);
          amrex::FArrayBox precg_fab(grown_box, 1);

          auto const& rainprod_arr = rainprod_fab.array();
          auto const& evapprod_arr = evapprod_fab.array();
          auto const& qlsink_arr = qlsink_fab.array();
          auto const& precr_arr = precr_fab.array();
          auto const& preci_arr = preci_fab.array();
          auto const& precs_arr = precs_fab.array();
          auto const& precg_arr = precg_fab.array();

          // Initialize WRF-Chem arrays to zero
          rainprod_fab.setVal(0.0);
          evapprod_fab.setVal(0.0);
          qlsink_fab.setVal(0.0);
          precr_fab.setVal(0.0);
          preci_fab.setVal(0.0);
          precs_fab.setVal(0.0);
          precg_fab.setVal(0.0);

          // Prepare data pointers for Fortran call
          // These would be passed directly to the Fortran interface
          double dummy_reflectivity = 0.0;
          double* dummy_reflectivity_ptr = &dummy_reflectivity;
          // Example call (pseudo-code - actual interface would depend on your Fortran interop setup)

          // Microphysics options/switches
          int m_iact = 2;    // CCN activation option (1: power-law, 2: lognormal aerosol)
          int m_inum = 1;    // Droplet number option (0: predict, 1: constant)
          int m_iliq = 0;    // Liquid-only option (0: include ice, 1: liquid only)
          int m_inuc = 0;    // Ice nucleation option (0: mid-latitude, 1: arctic)
          int m_ibase = 2;   // Cloud base activation option
          int m_isub = 0;    // Sub-grid vertical velocity option
          int m_igraup = 0;  // Graupel option (0: include graupel, 1: no graupel)
          int m_ihail = 0;   // Graupel/hail option (0: graupel, 1: hail)
          bool m_do_radar_ref = false;  // Radar reflectivity calculation flag

          // Physical constants
          amrex::Real m_pi;          // Pi constant
          amrex::Real m_R;           // Gas constant for dry air (J/kg/K)
          amrex::Real m_Rv;          // Gas constant for water vapor (J/kg/K)
          amrex::Real m_cp;          // Specific heat at constant pressure (J/kg/K)
          amrex::Real m_g;           // Gravitational acceleration (m/s^2)
          amrex::Real m_ep_2;        // Molecular weight ratio (Rd/Rv)

          // Reference density and species densities
          amrex::Real m_rhosu;       // Standard air density at 850 mb (kg/m^3)
          amrex::Real m_rhow;        // Density of liquid water (kg/m^3)
          amrex::Real m_rhoi;        // Bulk density of cloud ice (kg/m^3)
          amrex::Real m_rhosn;       // Bulk density of snow (kg/m^3)
          amrex::Real m_rhog;        // Bulk density of graupel/hail (kg/m^3)

          // Fall speed parameters (V=AD^B)
          amrex::Real m_ai, m_bi;    // Cloud ice fall speed parameters
          amrex::Real m_ac, m_bc;    // Cloud droplet fall speed parameters
          amrex::Real m_as, m_bs;    // Snow fall speed parameters
          amrex::Real m_ar, m_br;    // Rain fall speed parameters
          amrex::Real m_ag, m_bg;    // Graupel/hail fall speed parameters

          // Microphysical parameters
          amrex::Real m_aimm;        // Parameter in Bigg immersion freezing
          amrex::Real m_bimm;        // Parameter in Bigg immersion freezing
          amrex::Real m_ecr;         // Collection efficiency between droplets/rain and snow/rain
          amrex::Real m_dcs;         // Threshold size for cloud ice autoconversion (m)
          amrex::Real m_mi0;         // Initial mass of nucleated ice crystal (kg)
          amrex::Real m_mg0;         // Mass of embryo graupel (kg)
          amrex::Real m_f1s;         // Ventilation parameter for snow
          amrex::Real m_f2s;         // Ventilation parameter for snow
          amrex::Real m_f1r;         // Ventilation parameter for rain
          amrex::Real m_f2r;         // Ventilation parameter for rain
          amrex::Real m_qsmall;      // Smallest allowed hydrometeor mixing ratio
          amrex::Real m_eii;         // Collection efficiency, ice-ice collisions
          amrex::Real m_eci;         // Collection efficiency, ice-droplet collisions
          amrex::Real m_cpw;         // Specific heat of liquid water (J/kg/K)
          amrex::Real m_rin;         // Radius of contact nuclei (m)
          amrex::Real m_mmult;       // Mass of splintered ice particle (kg)

          // Size distribution parameters
          amrex::Real m_ci, m_di;    // Cloud ice size distribution parameters
          amrex::Real m_cs, m_ds;    // Snow size distribution parameters
          amrex::Real m_cg, m_dg;    // Graupel size distribution parameters

          // Lambda limits for size distributions
          amrex::Real m_lammaxi, m_lammini;    // Cloud ice lambda limits
          amrex::Real m_lammaxr, m_lamminr;    // Rain lambda limits
          amrex::Real m_lammaxs, m_lammins;    // Snow lambda limits
          amrex::Real m_lammaxg, m_lamming;    // Graupel lambda limits

          // Constant droplet concentration (if INUM = 1)
          amrex::Real m_ndcnst = 250.0;  // Droplet number concentration (cm^-3)

          // CCN spectra parameters (for IACT = 1)
          amrex::Real m_k1;          // Exponent in CCN activation formula
          amrex::Real m_c1;          // Coefficient in CCN activation formula (cm^-3)

          // Aerosol activation parameters (for IACT = 2)
          amrex::Real m_mw;          // Molecular weight water (kg/mol)
          amrex::Real m_osm;         // Osmotic coefficient
          amrex::Real m_vi;          // Number of ions dissociated in solution
          amrex::Real m_epsm;        // Aerosol soluble fraction
          amrex::Real m_rhoa;        // Aerosol bulk density (kg/m^3)
          amrex::Real m_map;         // Molecular weight aerosol (kg/mol)
          amrex::Real m_ma;          // Molecular weight of air (kg/mol)
          amrex::Real m_rr;          // Universal gas constant (J/mol/K)
          amrex::Real m_bact;        // Activation parameter
          amrex::Real m_rm1;         // Geometric mean radius, mode 1 (m)
          amrex::Real m_rm2;         // Geometric mean radius, mode 2 (m)
          amrex::Real m_nanew1;      // Total aerosol concentration, mode 1 (m^-3)
          amrex::Real m_nanew2;      // Total aerosol concentration, mode 2 (m^-3)
          amrex::Real m_sig1;        // Standard deviation of aerosol dist, mode 1
          amrex::Real m_sig2;        // Standard deviation of aerosol dist, mode 2
          amrex::Real m_f11;         // Correction factor for activation, mode 1
          amrex::Real m_f12;         // Correction factor for activation, mode 1
          amrex::Real m_f21;         // Correction factor for activation, mode 2
          amrex::Real m_f22;         // Correction factor for activation, mode 2

          // Precomputed constants for efficiency
          amrex::Real m_cons1, m_cons2, m_cons3, m_cons4, m_cons5;
          amrex::Real m_cons6, m_cons7, m_cons8, m_cons9, m_cons10;
          amrex::Real m_cons11, m_cons12, m_cons13, m_cons14, m_cons15;
          amrex::Real m_cons16, m_cons17, m_cons18, m_cons19, m_cons20;
          amrex::Real m_cons21, m_cons22, m_cons23, m_cons24, m_cons25;
          amrex::Real m_cons26, m_cons27, m_cons28, m_cons29, m_cons30;
          amrex::Real m_cons31, m_cons32, m_cons33, m_cons34, m_cons35;
          amrex::Real m_cons36, m_cons37, m_cons38, m_cons39, m_cons40;
          amrex::Real m_cons41;

          // Radar reflectivity parameters
          amrex::Real m_xam_r, m_xbm_r, m_xmu_r;  // Rain reflectivity parameters
          amrex::Real m_xam_s, m_xbm_s, m_xmu_s;  // Snow reflectivity parameters
          amrex::Real m_xam_g, m_xbm_g, m_xmu_g;  // Graupel reflectivity parameters
          amrex::Real m_lambda_radar;  // Radar wavelength (10 cm)
          amrex::Real m_k_w;           // K_w parameter for liquid water
          amrex::Real m_lamda4;        // Lambda^4 for radar calculation
          amrex::Real m_pi5;           // Pi^5 for radar calculation

          // Radar arrays for Simpson integration
          amrex::Vector<amrex::Real> m_xxds, m_xxdg;     // Diameter values
          amrex::Vector<amrex::Real> m_xdts, m_xdtg;     // Diameter increments
          amrex::Vector<amrex::Real> m_simpson;          // Simpson's rule weights

          // Parameters for radar melting calculations
          bool m_melt_outside_s;     // Liquid coating of melting snow
          bool m_melt_outside_g;     // Liquid coating of melting graupel
          std::complex<amrex::Real> m_m_w_0;  // Dielectric constant for water
          std::complex<amrex::Real> m_m_i_0;  // Dielectric constant for ice

          // Set microphysics control parameters
          m_inum = 1;           // Use constant droplet number concentration
          m_ndcnst = 250.0;     // Droplet number concentration (cm^-3)
          // Mathematical constants
          m_pi = 3.1415926535897932384626434;

          m_R = 287.0;         // Gas constant for dry air (J/kg/K)
          m_Rd = 287.0;         // Gas constant for dry air (J/kg/K)
          m_Rv = 461.6;        // Gas constant for water vapor (J/kg/K)
          m_cp = 7.0*287.0/2.0;        // Specific heat at constant pressure (J/kg/K)
          m_g = 9.81;           // Gravitational acceleration (m/s^2)
          m_ep_2 = m_Rd / m_Rv;     // Molecular weight ratio (Rd/Rv)

          // Reference density
          m_rhosu = 85000.0/(m_Rd*273.15);  // Standard air density at 850 mb (kg/m^3)

          // Densities for different hydrometeor species
          m_rhow = 997.0;     // Density of liquid water (kg/m^3)
          m_rhoi = 500.0;     // Bulk density of cloud ice (kg/m^3)
          m_rhosn = 100.0;    // Bulk density of snow (kg/m^3)

          // Set density for graupel or hail based on configuration
          if (m_ihail == 0) {
            m_rhog = 400.0; // Bulk density of graupel (kg/m^3)
          } else {
            m_rhog = 900.0; // Bulk density of hail (kg/m^3)
          }

          // Fall speed parameters (V=AD^B) for different hydrometeors
          // Cloud ice
          m_ai = 700.0;
          m_bi = 1.0;

          // Cloud droplets
          m_ac = 3.0E7;
          m_bc = 2.0;

          // Snow
          m_as = 11.72;
          m_bs = 0.41;

          // Rain
          m_ar = 841.99667;
          m_br = 0.8;

          // Graupel/hail (dependent on configuration)
          if (m_ihail == 0) {
            // Graupel parameters
            m_ag = 19.3;
            m_bg = 0.37;
          } else {
            // Hail parameters (Matsun and Huggins 1980)
            m_ag = 114.5;
            m_bg = 0.5;
          }

          // Microphysical parameters
          m_aimm = 0.66;       // Parameter in Bigg immersion freezing
          m_bimm = 100.0;      // Parameter in Bigg immersion freezing
          m_ecr = 1.0;         // Collection efficiency between rain and snow/graupel
          m_dcs = 125.0E-6;    // Threshold size for cloud ice autoconversion (m)
          m_mi0 = 4.0/3.0*m_pi*m_rhoi*std::pow(10.0E-6, 3);  // Initial mass of nucleated ice crystal (kg)
          m_mg0 = 1.6E-10;     // Mass of embryo graupel (kg)

          // Ventilation parameters
          m_f1s = 0.86;        // Ventilation parameter for snow
          m_f2s = 0.28;        // Ventilation parameter for snow
          m_f1r = 0.78;        // Ventilation parameter for rain
          m_f2r = 0.308;       // Ventilation parameter for rain

          // Smallest allowed hydrometeor mixing ratio
          m_qsmall = 1.0E-14;

          // Collection efficiencies
          m_eii = 0.1;         // Ice-ice collision efficiency
          m_eci = 0.7;         // Ice-droplet collision efficiency

          // Specific heat of liquid water (J/kg/K)
          m_cpw = 4187.0;

          // Size distribution parameters
          m_ci = m_rhoi * m_pi / 6.0;
          m_di = 3.0;
          m_cs = m_rhosn * m_pi / 6.0;
          m_ds = 3.0;
          m_cg = m_rhog * m_pi / 6.0;
          m_dg = 3.0;

          // Radius of contact nuclei (m)
          m_rin = 0.1E-6;

          // Mass of splintered ice particle (kg)
          m_mmult = 4.0/3.0*m_pi*m_rhoi*std::pow(5.0E-6, 3);

          // Set lambda limits for size distributions
          // Maximum and minimum values for lambda parameter in size distributions
          m_lammaxi = 1.0/1.0E-6;
          m_lammini = 1.0/(2.0*m_dcs + 100.0E-6);
          m_lammaxr = 1.0/20.0E-6;
          m_lamminr = 1.0/2800.0E-6;
          m_lammaxs = 1.0/10.0E-6;
          m_lammins = 1.0/2000.0E-6;
          m_lammaxg = 1.0/20.0E-6;
          m_lamming = 1.0/2000.0E-6;

          // Set CCN parameters for different environments
          if (m_activate_type == 1) {
            // Maritime CCN spectrum parameters (modified from Rasmussen et al. 2002)
            // NCCN = C*S^K, where S is supersaturation in %
            m_k1 = 0.4;        // Exponent in CCN activation formula
            m_c1 = 120.0;      // Coefficient in CCN activation formula (cm^-3)
          }

          // Initialize aerosol activation parameters for lognormal distribution
          if (m_activate_type == 2) {
            // Parameters for ammonium sulfate
            m_mw = 0.018;      // Molecular weight of water (kg/mol)
            m_osm = 1.0;       // Osmotic coefficient
            m_vi = 3.0;        // Number of ions dissociated in solution
            m_epsm = 0.7;      // Aerosol soluble fraction
            m_rhoa = 1777.0;   // Aerosol bulk density (kg/m^3)
            m_map = 0.132;     // Molecular weight of aerosol (kg/mol)
            m_ma = 0.0284;     // Molecular weight of air (kg/mol)
            m_rr = 8.3145;     // Universal gas constant (J/mol/K)
            m_bact = m_vi * m_osm * m_epsm * m_mw * m_rhoa / (m_map * m_rhow);
            m_a_w = 2.0 * m_mw * 0.0761 / (m_rhow * m_r_v * 293.15);  // "A" parameter

            // Aerosol size distribution parameters for MPACE (Morrison et al. 2007, JGR)
            // Mode 1
            m_rm1 = 0.052E-6;  // Geometric mean radius, mode 1 (m)
            m_sig1 = 2.04;     // Standard deviation of aerosol size distribution, mode 1
            m_nanew1 = 72.2E6; // Total aerosol concentration, mode 1 (m^-3)
            m_f11 = 0.5 * std::exp(2.5 * std::pow(std::log(m_sig1), 2));
            m_f21 = 1.0 + 0.25 * std::log(m_sig1);

            // Mode 2
            m_rm2 = 1.3E-6;    // Geometric mean radius, mode 2 (m)
            m_sig2 = 2.5;      // Standard deviation of aerosol size distribution, mode 2
            m_nanew2 = 1.8E6;  // Total aerosol concentration, mode 2 (m^-3)
            m_f12 = 0.5 * std::exp(2.5 * std::pow(std::log(m_sig2), 2));
            m_f22 = 1.0 + 0.25 * std::log(m_sig2);
          }

          // Precompute constants for efficiency
          m_cons1 = gamma_function(1.0 + m_ds) * m_cs;
          m_cons2 = gamma_function(1.0 + m_dg) * m_cg;
          m_cons3 = gamma_function(4.0 + m_bs) / 6.0;
          m_cons4 = gamma_function(4.0 + m_br) / 6.0;
          m_cons5 = gamma_function(1.0 + m_bs);
          m_cons6 = gamma_function(1.0 + m_br);
          m_cons7 = gamma_function(4.0 + m_bg) / 6.0;
          m_cons8 = gamma_function(1.0 + m_bg);
          m_cons9 = gamma_function(5.0/2.0 + m_br/2.0);
          m_cons10 = gamma_function(5.0/2.0 + m_bs/2.0);
          m_cons11 = gamma_function(5.0/2.0 + m_bg/2.0);
          m_cons12 = gamma_function(1.0 + m_di) * m_ci;
          m_cons13 = gamma_function(m_bs + 3.0) * m_pi / 4.0 * m_eci;
          m_cons14 = gamma_function(m_bg + 3.0) * m_pi / 4.0 * m_eci;
          m_cons15 = -1108.0 * m_eii * std::pow(m_pi, (1.0-m_bs)/3.0) *
            std::pow(m_rhosn, (-2.0-m_bs)/3.0) / (4.0*720.0);
          m_cons16 = gamma_function(m_bi + 3.0) * m_pi / 4.0 * m_eci;
          m_cons17 = 4.0 * 2.0 * 3.0 * m_rhosu * m_pi * m_eci * m_eci *
            gamma_function(2.0*m_bs + 2.0) / (8.0*(m_rhog-m_rhosn));
          m_cons18 = m_rhosn * m_rhosn;
          m_cons19 = m_rhow * m_rhow;
          m_cons20 = 20.0 * m_pi * m_pi * m_rhow * m_bimm;
          m_cons21 = 4.0 / (m_dcs * m_rhoi);
          m_cons22 = m_pi * m_rhoi * std::pow(m_dcs, 3) / 6.0;
          m_cons23 = m_pi / 4.0 * m_eii * gamma_function(m_bs + 3.0);
          m_cons24 = m_pi / 4.0 * m_ecr * gamma_function(m_br + 3.0);
          m_cons25 = m_pi * m_pi / 24.0 * m_rhow * m_ecr * gamma_function(m_br + 6.0);
          m_cons26 = m_pi / 6.0 * m_rhow;
          m_cons27 = gamma_function(1.0 + m_bi);
          m_cons28 = gamma_function(4.0 + m_bi) / 6.0;
          m_cons29 = 4.0/3.0 * m_pi * m_rhow * std::pow(25.0E-6, 3);
          m_cons30 = 4.0/3.0 * m_pi * m_rhow;
          m_cons31 = m_pi * m_pi * m_ecr * m_rhosn;
          m_cons32 = m_pi / 2.0 * m_ecr;
          m_cons33 = m_pi * m_pi * m_ecr * m_rhog;
          m_cons34 = 5.0/2.0 + m_br/2.0;
          m_cons35 = 5.0/2.0 + m_bs/2.0;
          m_cons36 = 5.0/2.0 + m_bg/2.0;
          m_cons37 = 4.0 * m_pi * 1.38E-23 / (6.0 * m_pi * m_rin);
          m_cons38 = m_pi * m_pi / 3.0 * m_rhow;
          m_cons39 = m_pi * m_pi / 36.0 * m_rhow * m_bimm;
          m_cons40 = m_pi / 6.0 * m_bimm;
          m_cons41 = m_pi * m_pi * m_ecr * m_rhow;

          // Set CCN parameters for different environments
          if (m_iact == 1) {
            // Maritime CCN spectrum parameters (modified from Rasmussen et al. 2002)
            // NCCN = C*S^K, where S is supersaturation in %
            m_k1 = 0.4;        // Exponent in CCN activation formula
            m_c1 = 120.0;      // Coefficient in CCN activation formula (cm^-3)
          }

          // Initialize aerosol activation parameters for IACT=2
          if (m_iact == 2) {
            // Parameters for ammonium sulfate
            m_mw = 0.018;      // Molecular weight of water (kg/mol)
            m_osm = 1.0;       // Osmotic coefficient
            m_vi = 3.0;        // Number of ions dissociated in solution
            m_epsm = 0.7;      // Aerosol soluble fraction
            m_rhoa = 1777.0;   // Aerosol bulk density (kg/m^3)
            m_map = 0.132;     // Molecular weight of aerosol (kg/mol)
            m_ma = 0.0284;     // Molecular weight of air (kg/mol)
            m_rr = 8.3145;     // Universal gas constant (J/mol/K)
            m_bact = m_vi * m_osm * m_epsm * m_mw * m_rhoa / (m_map * m_rhow);

            // Aerosol size distribution parameters for MPACE (Morrison et al. 2007, JGR)
            // Mode 1
            m_rm1 = 0.052E-6;  // Geometric mean radius, mode 1 (m)
            m_sig1 = 2.04;     // Standard deviation of aerosol size distribution, mode 1
            m_nanew1 = 72.2E6; // Total aerosol concentration, mode 1 (m^-3)
            m_f11 = 0.5 * std::exp(2.5 * std::pow(std::log(m_sig1), 2));
            m_f21 = 1.0 + 0.25 * std::log(m_sig1);

            // Mode 2
            m_rm2 = 1.3E-6;    // Geometric mean radius, mode 2 (m)
            m_sig2 = 2.5;      // Standard deviation of aerosol size distribution, mode 2
            m_nanew2 = 1.8E6;  // Total aerosol concentration, mode 2 (m^-3)
            m_f12 = 0.5 * std::exp(2.5 * std::pow(std::log(m_sig2), 2));
            m_f22 = 1.0 + 0.25 * std::log(m_sig2);
          }

          // Initialize radar reflectivity parameters
          m_xam_r = m_pi * m_rhow / 6.0;
          m_xbm_r = 3.0;
          m_xmu_r = 0.0;
          m_xam_s = m_cs;
          m_xbm_s = m_ds;
          m_xmu_s = 0.0;
          m_xam_g = m_cg;
          m_xbm_g = m_dg;
          m_xmu_g = 0.0;

          // Set microphysics control parameters
          m_activate_type = 2;  // Lognormal aerosol activation
          m_inuc_type = 0;      // Mid-latitude ice nucleation (Cooper)
          m_iliq = 0;           // Include ice processes
          m_igraup = 0;         // Include graupel processes
          m_ihail = 0;          // Use graupel (0) instead of hail (1)
          m_isub = 0;           // Sub-grid vertical velocity option
          m_do_radar_ref = false; // Disable radar reflectivity by default

          // Ensure consistency between m_iact and m_activate_type
          m_iact = m_activate_type;

          // Initialize water vapor gas constant for activation calculations
          m_r_v = m_Rv;

          //          amrex::Print()<<amrex::FArrayBox(qv_arr)<<std::endl;
          mp_morr_two_moment_c
          (
              1,  // ITIMESTEP - Use 1 for simplicity

              // 3D arrays in Fortran expected order (assume column-major for Fortran)
              theta_arr.dataPtr(),      // TH
              qv_arr.dataPtr(),         // QV
              qcl_arr.dataPtr(),        // QC
              qpr_arr.dataPtr(),        // QR
              qci_arr.dataPtr(),        // QI
              qps_arr.dataPtr(),        // QS
              qpg_arr.dataPtr(),        // QG
              ni_arr.dataPtr(),         // NI
              ns_arr.dataPtr(),         // NS
              nr_arr.dataPtr(),         // NR
              ng_arr.dataPtr(),         // NG

              rho_arr.dataPtr(),        // RHO
              pii_arr.dataPtr(),        // PII (Exner function)
              pres_arr.dataPtr(),       // P (in hPa, convert if needed)
              dt,                       // DT_IN
              dz_arr.dataPtr(),         // DZ
              w_arr.dataPtr(),          // W (vertical velocity)

              // 2D arrays for precipitation accounting
              rain_accum_arr.dataPtr(), // RAINNC
              rainncv_arr.dataPtr(),    // RAINNCV
              sr_arr.dataPtr(),         // SR
              snow_accum_arr.dataPtr(), // SNOWNC
              snowncv_arr.dataPtr(),    // SNOWNCV
              graup_accum_arr.dataPtr(),// GRAUPELNC
              graupelncv_arr.dataPtr(), // GRAUPELNCV

              // Radar reflectivity
              dummy_reflectivity_ptr,  // refl_10cm
              true,                     // diagflag
              0,   // do_radar_ref

              // Cumulus tendencies
              qrcuten_arr.dataPtr(),    // qrcuten
              qscuten_arr.dataPtr(),    // qscuten
              qicuten_arr.dataPtr(),    // qicuten

              // WRF-Chem flags
              flag_qndrop,              // F_QNDROP
              nullptr,                  // qndrop (not used here)
              ht_arr.dataPtr(),         // HT (terrain height - not used)

              // Domain dimensions
              ilo, ihi, jlo, jhi, klo, khi,  // IDS,IDE,JDS,JDE,KDS,KDE
              ilom, ihim, jlom, jhim, klom, khim,  // IMS,IME,JMS,JME,KMS,KME
              ilo, ihi, jlo, jhi, klo, khi,  // ITS,ITE,JTS,JTE,KTS,KTE

              // Optional WRF-Chem outputs
              false,                    // wetscav_on
              rainprod_arr.dataPtr(),   // rainprod
              evapprod_arr.dataPtr(),   // evapprod
              qlsink_arr.dataPtr(),     // QLSINK
              precr_arr.dataPtr(),      // PRECR
              preci_arr.dataPtr(),      // PRECI
              precs_arr.dataPtr(),      // PRECS
              precg_arr.dataPtr()       // PRECG
          );
          //          amrex::Print()<<amrex::FArrayBox(qv_arr)<<std::endl;
          // After the call, all fields are updated
          // We don't need to copy results back since we passed direct pointers
          // to our class member arrays
        }
#endif
    }

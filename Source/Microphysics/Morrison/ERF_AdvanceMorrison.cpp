#include <string>
#include <vector>
#include <memory>
#include <complex>
#include <cmath>

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

using namespace amrex;


// Gamma function using the standard C++ tgamma function
constexpr Real gamma_function(Real x) {
  // Note: std::tgamma may not be constexpr in all compilers
  // but this will work at runtime in any case
  return std::tgamma(x);
}
  /**
   * Helper function to calculate saturation vapor pressure for water or ice.
   * This corresponds to the POLYSVP function in the Fortran code (line ~5580).
   *
   * @param[in] T Temperature in Kelvin
   * @param[in] type 0 for liquid water, 1 for ice
   * @return Saturation vapor pressure in Pascals
   */
  amrex::Real
  calc_saturation_vapor_pressure (const amrex::Real T, const int type)
  {
    amrex::Real polysvp = 0.0;
    amrex::Real del_T = T - 273.15;  // Convert to Celsius

    if (type == 1) {  // Ice (lines ~5631-5644)
        if (T >= 195.8) {
            // Flatau et al. formula for ice
            const amrex::Real a0i = 6.11147274;
            const amrex::Real a1i = 0.503160820;
            const amrex::Real a2i = 0.188439774e-1;
            const amrex::Real a3i = 0.420895665e-3;
            const amrex::Real a4i = 0.615021634e-5;
            const amrex::Real a5i = 0.602588177e-7;
            const amrex::Real a6i = 0.385852041e-9;
            const amrex::Real a7i = 0.146898966e-11;
            const amrex::Real a8i = 0.252751365e-14;

            polysvp = a0i + del_T*(a1i + del_T*(a2i + del_T*(a3i + del_T*(a4i + del_T*(a5i + del_T*(a6i + del_T*(a7i + a8i*del_T)))))));
            polysvp *= 100.0;  // Convert from hPa to Pa
        } else {
            // Goff-Gratch formula for ice at cold temperatures
            polysvp = std::pow(10.0, (-9.09718*(273.16/T-1.0) - 3.56654*std::log10(273.16/T) +
                                      0.876793*(1.0-T/273.16) + std::log10(6.1071))) * 100.0;
            polysvp = 0.0;
        } // T
    } else {  // Water (lines ~5648-5665)
      if (T >= 202.0) {
        // Flatau et al. formula for liquid water
        const amrex::Real a0 = 6.11239921;
        const amrex::Real a1 = 0.443987641;
        const amrex::Real a2 = 0.142986287e-1;
        const amrex::Real a3 = 0.264847430e-3;
        const amrex::Real a4 = 0.302950461e-5;
        const amrex::Real a5 = 0.206739458e-7;
        const amrex::Real a6 = 0.640689451e-10;
        const amrex::Real a7 = -0.952447341e-13;
        const amrex::Real a8 = -0.976195544e-15;

        polysvp = a0 + del_T*(a1 + del_T*(a2 + del_T*(a3 + del_T*(a4 + del_T*(a5 + del_T*(a6 + del_T*(a7 + a8*del_T)))))));
        polysvp *= 100.0;  // Convert from hPa to Pa
      } else {
        // Goff-Gratch formula for water at cold temperatures
        polysvp = std::pow(10.0, (-7.90298*(373.16/T-1.0) + 5.02808*std::log10(373.16/T) -
                                  1.3816e-7*(std::pow(10.0, (11.344*(1.0-T/373.16)))-1.0) +
                                  8.1328e-3*(std::pow(10.0, (-3.49149*(373.16/T-1.0)))-1.0) +
                                  std::log10(1013.246))) * 100.0;
      }
    }

    return polysvp;
  }

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
          auto const& nc_arr = mic_fab_vars[MicVar_Morr::nc]->array(mfi);
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
          amrex::Real m_Rd;           // Gas constant for dry air (J/kg/K)
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
          if (m_iact == 1) {
            // Maritime CCN spectrum parameters (modified from Rasmussen et al. 2002)
            // NCCN = C*S^K, where S is supersaturation in %
            m_k1 = 0.4;        // Exponent in CCN activation formula
            m_c1 = 120.0;      // Coefficient in CCN activation formula (cm^-3)
          }

          // Initialize aerosol activation parameters for lognormal distribution
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
            //            m_a_w = 2.0 * m_mw * 0.0761 / (m_rhow * m_r_v * 293.15);  // "A" parameter

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
          m_iact = 2;  // Lognormal aerosol activation
          m_inuc = 0;      // Mid-latitude ice nucleation (Cooper)
          m_iliq = 0;           // Include ice processes
          m_igraup = 0;         // Include graupel processes
          m_ihail = 0;          // Use graupel (0) instead of hail (1)
          m_isub = 0;           // Sub-grid vertical velocity option
          m_do_radar_ref = false; // Disable radar reflectivity by default
          FILE *file = fopen("output_cpp.txt", "a");
          ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
         {
            // Tendencies and mixing ratios
            amrex::Real qc3dten=0;            // QC3DTEN: CLOUD WATER MIXING RATIO TENDENCY (KG/KG/S)
            amrex::Real qi3dten=0;            // QI3DTEN: CLOUD ICE MIXING RATIO TENDENCY (KG/KG/S)
            amrex::Real qni3dten=0;           // QNI3DTEN: SNOW MIXING RATIO TENDENCY (KG/KG/S)
            amrex::Real qr3dten=0;            // QR3DTEN: RAIN MIXING RATIO TENDENCY (KG/KG/S)
            amrex::Real ni3dten=0;            // NI3DTEN: CLOUD ICE NUMBER CONCENTRATION (1/KG/S)
            amrex::Real ns3dten=0;            // NS3DTEN: SNOW NUMBER CONCENTRATION (1/KG/S)
            amrex::Real nr3dten=0;            // NR3DTEN: RAIN NUMBER CONCENTRATION (1/KG/S)
            amrex::Real qc3d=qcl_arr(i,j,k);               // QC3D: CLOUD WATER MIXING RATIO (KG/KG)
            amrex::Real qi3d=qci_arr(i,j,k);               // QI3D: CLOUD ICE MIXING RATIO (KG/KG)
            amrex::Real qni3d=qps_arr(i,j,k);;              // QNI3D: SNOW MIXING RATIO (KG/KG)
            amrex::Real qr3d=qpr_arr(i,j,k);;               // QR3D: RAIN MIXING RATIO (KG/KG)
            amrex::Real ni3d=ni_arr(i,j,k);               // NI3D: CLOUD ICE NUMBER CONCENTRATION (1/KG)
            amrex::Real ns3d=ns_arr(i,j,k);               // NS3D: SNOW NUMBER CONCENTRATION (1/KG)
            amrex::Real nr3d=nr_arr(i,j,k);               // NR3D: RAIN NUMBER CONCENTRATION (1/KG)
            amrex::Real t3dten=0;             // T3DTEN: TEMPERATURE TENDENCY (K/S)
            amrex::Real qv3dten=0;            // QV3DTEN: WATER VAPOR MIXING RATIO TENDENCY (KG/KG/S)
            amrex::Real t3d=theta_arr(i,j,k)*pii_arr(i,j,k);                // T3D: TEMPERATURE (K)
            amrex::Real qv3d=qv_arr(i,j,k);               // QV3D: WATER VAPOR MIXING RATIO (KG/KG)
            amrex::Real pres=pres_arr(i,j,k);               // PRES: ATMOSPHERIC PRESSURE (PA)
            amrex::Real dzq=dz_arr(i,j,k);// dz_val;                // DZQ: DIFFERENCE IN HEIGHT ACROSS LEVEL (m)
            amrex::Real w3d=w_arr(i,j,k);                // W3D: GRID-SCALE VERTICAL VELOCITY (M/S)

            // WRF-chem variables
            amrex::Real nc3d;               // nc3d: Cloud droplet number concentration
            amrex::Real nc3dten=0;            // nc3dten: Cloud droplet number concentration tendency
            int iinum;                      // iinum: Integer control variable

            // Graupel variables
            amrex::Real qg3dten=0;            // QG3DTEN: GRAUPEL MIX RATIO TENDENCY (KG/KG/S)
            amrex::Real ng3dten=0;            // NG3DTEN: GRAUPEL NUMB CONC TENDENCY (1/KG/S)
            amrex::Real qg3d=qpg_arr(i,j,k);               // QG3D: GRAUPEL MIX RATIO (KG/KG)
            amrex::Real ng3d=ng_arr(i,j,k);               // NG3D: GRAUPEL NUMBER CONC (1/KG)

            // Sedimentation tendencies
            amrex::Real qgsten=0.0;             // QGSTEN: GRAUPEL SED TEND (KG/KG/S)
            amrex::Real qrsten=0.0;             // QRSTEN: RAIN SED TEND (KG/KG/S)
            amrex::Real qisten=0.0;             // QISTEN: CLOUD ICE SED TEND (KG/KG/S)
            amrex::Real qnisten=0.0;            // QNISTEN: SNOW SED TEND (KG/KG/S)
            amrex::Real qcsten=0.0;             // QCSTEN: CLOUD WAT SED TEND (KG/KG/S)

            // Cumulus tendencies for precipitation
            amrex::Real qrcu1d=qrcuten_arr(i,j,k);             // qrcu1d: Rain from cumulus parameterization
            amrex::Real qscu1d=qscuten_arr(i,j,k);             // qscu1d: Snow from cumulus parameterization
            amrex::Real qicu1d=qicuten_arr(i,j,k);             // qicu1d: Ice from cumulus parameterization

            // Output variables
            amrex::Real precrt;             // PRECRT: TOTAL PRECIP PER TIME STEP (mm)
            amrex::Real snowrt;             // SNOWRT: SNOW PER TIME STEP (mm)
            amrex::Real snowprt;            // SNOWPRT: TOTAL CLOUD ICE PLUS SNOW PER TIME STEP (mm)
            amrex::Real grplprt;            // GRPLPRT: TOTAL GRAUPEL PER TIME STEP (mm)

            // Effective radii
            amrex::Real effc;               // EFFC: DROPLET EFFECTIVE RADIUS (MICRON)
            amrex::Real effi;               // EFFI: CLOUD ICE EFFECTIVE RADIUS (MICRON)
            amrex::Real effs;               // EFFS: SNOW EFFECTIVE RADIUS (MICRON)
            amrex::Real effr;               // EFFR: RAIN EFFECTIVE RADIUS (MICRON)
            amrex::Real effg;               // EFFG: GRAUPEL EFFECTIVE RADIUS (MICRON)

            // Model input parameters
            amrex::Real dt;                 // DT: MODEL TIME STEP (SEC)

            // Size distribution parameters
            amrex::Real lamc;               // LAMC: Slope parameter for droplets (m^-1)
            amrex::Real lami;               // LAMI: Slope parameter for cloud ice (m^-1)
            amrex::Real lams;               // LAMS: Slope parameter for snow (m^-1)
            amrex::Real lamr;               // LAMR: Slope parameter for rain (m^-1)
            amrex::Real lamg;               // LAMG: Slope parameter for graupel (m^-1)
            amrex::Real cdist1;             // CDIST1: PSD parameter for droplets
            amrex::Real n0i;                // N0I: Intercept parameter for cloud ice (kg^-1 m^-1)
            amrex::Real n0s;                // N0S: Intercept parameter for snow (kg^-1 m^-1)
            amrex::Real n0r;                // N0RR: Intercept parameter for rain (kg^-1 m^-1)
            amrex::Real n0g;                // N0G: Intercept parameter for graupel (kg^-1 m^-1)
            amrex::Real pgam;               // PGAM: Spectral shape parameter for droplets

            // Microphysical processes
            amrex::Real nsubc;              // NSUBC: Loss of NC during evaporation
            amrex::Real nsubi;              // NSUBI: Loss of NI during sublimation
            amrex::Real nsubs;              // NSUBS: Loss of NS during sublimation
            amrex::Real nsubr;              // NSUBR: Loss of NR during evaporation
            amrex::Real prd;                // PRD: Deposition cloud ice
            amrex::Real pre;                // PRE: Evaporation of rain
            amrex::Real prds;               // PRDS: Deposition snow
            amrex::Real nnuccc;             // NNUCCC: Change N due to contact freezing droplets
            amrex::Real mnuccc;             // MNUCCC: Change Q due to contact freezing droplets
            amrex::Real pra;                // PRA: Accretion droplets by rain
            amrex::Real prc;                // PRC: Autoconversion droplets
            amrex::Real pcc;                // PCC: Condensation/evaporation droplets
            amrex::Real nnuccd;             // NNUCCD: Change N freezing aerosol (primary ice nucleation)
            amrex::Real mnuccd;             // MNUCCD: Change Q freezing aerosol (primary ice nucleation)
            amrex::Real mnuccr;             // MNUCCR: Change Q due to contact freezing rain
            amrex::Real nnuccr;             // NNUCCR: Change N due to contact freezing rain
            amrex::Real npra;               // NPRA: Change N due to droplet accretion by rain
            amrex::Real nragg;              // NRAGG: Self-collection/breakup of rain
            amrex::Real nsagg;              // NSAGG: Self-collection of snow
            amrex::Real nprc;               // NPRC: Change NC autoconversion droplets
            amrex::Real nprc1;              // NPRC1: Change NR autoconversion droplets
            amrex::Real prai;               // PRAI: Change Q accretion cloud ice by snow
            amrex::Real prci;               // PRCI: Change Q autoconversion cloud ice to snow
            amrex::Real psacws;             // PSACWS: Change Q droplet accretion by snow
            amrex::Real npsacws;            // NPSACWS: Change N droplet accretion by snow
            amrex::Real psacwi;             // PSACWI: Change Q droplet accretion by cloud ice
            amrex::Real npsacwi;            // NPSACWI: Change N droplet accretion by cloud ice
            amrex::Real nprci;              // NPRCI: Change N autoconversion cloud ice by snow
            amrex::Real nprai;              // NPRAI: Change N accretion cloud ice
            amrex::Real nmults;             // NMULTS: Ice multiplication due to riming droplets by snow
            amrex::Real nmultr;             // NMULTR: Ice multiplication due to riming rain by snow
            amrex::Real qmults;             // QMULTS: Change Q due to ice multiplication droplets/snow
            amrex::Real qmultr;             // QMULTR: Change Q due to ice multiplication rain/snow
            amrex::Real pracs;              // PRACS: Change Q rain-snow collection
            amrex::Real npracs;             // NPRACS: Change N rain-snow collection
            amrex::Real pccn;               // PCCN: Change Q droplet activation
            amrex::Real psmlt;              // PSMLT: Change Q melting snow to rain
            amrex::Real evpms;              // EVPMS: Change Q melting snow evaporating
            amrex::Real nsmlts;             // NSMLTS: Change N melting snow
            amrex::Real nsmltr;             // NSMLTR: Change N melting snow to rain
            amrex::Real piacr;              // PIACR: Change QR, ice-rain collection
            amrex::Real niacr;              // NIACR: Change N, ice-rain collection
            amrex::Real praci;              // PRACI: Change QI, ice-rain collection
            amrex::Real piacrs;             // PIACRS: Change QR, ice rain collision, added to snow
            amrex::Real niacrs;             // NIACRS: Change N, ice rain collision, added to snow
            amrex::Real pracis;             // PRACIS: Change QI, ice rain collision, added to snow
            amrex::Real eprd;               // EPRD: Sublimation cloud ice
            amrex::Real eprds;              // EPRDS: Sublimation snow

            // Graupel processes
            amrex::Real pracg;              // PRACG: Change in Q collection rain by graupel
            amrex::Real psacwg;             // PSACWG: Change in Q collection droplets by graupel
            amrex::Real pgsacw;             // PGSACW: Conversion Q to graupel due to collection droplets by snow
            amrex::Real pgracs;             // PGRACS: Conversion Q to graupel due to collection rain by snow
            amrex::Real prdg;               // PRDG: Deposition of graupel
            amrex::Real eprdg;              // EPRDG: Sublimation of graupel
            amrex::Real evpmg;              // EVPMG: Change Q melting of graupel and evaporation
            amrex::Real pgmlt;              // PGMLT: Change Q melting of graupel
            amrex::Real npracg;             // NPRACG: Change N collection rain by graupel
            amrex::Real npsacwg;            // NPSACWG: Change N collection droplets by graupel
            amrex::Real nscng;              // NSCNG: Change N conversion to graupel due to collection droplets by snow
            amrex::Real ngracs;             // NGRACS: Change N conversion to graupel due to collection rain by snow
            amrex::Real ngmltg;             // NGMLTG: Change N melting graupel
            amrex::Real ngmltr;             // NGMLTR: Change N melting graupel to rain
            amrex::Real nsubg;              // NSUBG: Change N sublimation/deposition of graupel
            amrex::Real psacr;              // PSACR: Conversion due to collection of snow by rain
            amrex::Real nmultg;             // NMULTG: Ice multiplication due to accretion droplets by graupel
            amrex::Real nmultrg;            // NMULTRG: Ice multiplication due to accretion rain by graupel
            amrex::Real qmultg;             // QMULTG: Change Q due to ice multiplication droplets/graupel
            amrex::Real qmultrg;            // QMULTRG: Change Q due to ice multiplication rain/graupel

            // Time-varying atmospheric parameters
            amrex::Real kap;                // KAP: Thermal conductivity of air
            amrex::Real evs;                // EVS: Saturation vapor pressure
            amrex::Real eis;                // EIS: Ice saturation vapor pressure
            amrex::Real qvs;                // QVS: Saturation mixing ratio
            amrex::Real qvi;                // QVI: Ice saturation mixing ratio
            amrex::Real qvqvs;              // QVQVS: Saturation ratio
            amrex::Real qvqvsi;             // QVQVSI: Ice saturation ratio
            amrex::Real dv;                 // DV: Diffusivity of water vapor in air
            amrex::Real xxls;               // XXLS: Latent heat of sublimation
            amrex::Real xxlv;               // XXLV: Latent heat of vaporization
            amrex::Real cpm;                // CPM: Specific heat at constant pressure for moist air
            amrex::Real mu;                 // MU: Viscosity of air
            amrex::Real sc;                 // SC: Schmidt number
            amrex::Real xlf;                // XLF: Latent heat of freezing
            amrex::Real rho;                // RHO: Air density
            amrex::Real ab;                 // AB: Correction to condensation rate due to latent heating
            amrex::Real abi;                // ABI: Correction to deposition rate due to latent heating

            // Time-varying microphysics parameters
            amrex::Real dap;                // DAP: Diffusivity of aerosol
            amrex::Real nacnt;              // NACNT: Number of contact nuclei
            amrex::Real fmult;              // FMULT: Temperature-dependent parameter for rime-splintering
            amrex::Real coffi;              // COFFI: Ice autoconversion parameter

            // Fall speed working variables
            amrex::Real dumi;               // DUMI: Dummy variable for ice
            amrex::Real dumr;               // DUMR: Dummy variable for rain
            amrex::Real dumfni;             // DUMFNI: Dummy fall speed number for ice
            amrex::Real dumg;               // DUMG: Dummy variable for graupel
            amrex::Real dumfng;             // DUMFNG: Dummy fall speed number for graupel
            amrex::Real uni;                // UNI: Number-weighted terminal velocity for cloud ice
            amrex::Real umi;                // UMI: Mass-weighted terminal velocity for cloud ice
            amrex::Real umr;                // UMR: Mass-weighted terminal velocity for rain
            amrex::Real fr;                 // FR: Mass-weighted fall speed for rain
            amrex::Real fi;                 // FI: Mass-weighted fall speed for ice
            amrex::Real fni;                // FNI: Number-weighted fall speed for ice
            amrex::Real fg;                 // FG: Mass-weighted fall speed for graupel
            amrex::Real fng;                // FNG: Number-weighted fall speed for graupel
            amrex::Real rgvm;               // RGVM: Rain size parameter
            amrex::Real faloutr;            // FALOUTR: Fallout rate for rain mass
            amrex::Real falouti;            // FALOUTI: Fallout rate for ice mass
            amrex::Real faloutni;           // FALOUTNI: Fallout rate for ice number
            amrex::Real faltndr;            // FALTNDR: Mass-weighted fallout tendency for rain
            amrex::Real faltndi;            // FALTNDI: Mass-weighted fallout tendency for ice
            amrex::Real faltndni;           // FALTNDNI: Number-weighted fallout tendency for ice
            amrex::Real rho2;               // RHO2: Air density squared
            amrex::Real dumqs;              // DUMQS: Dummy variable for snow
            amrex::Real dumfns;             // DUMFNS: Dummy fall speed number for snow
            amrex::Real ums;                // UMS: Mass-weighted terminal velocity for snow
            amrex::Real uns;                // UNS: Number-weighted terminal velocity for snow
            amrex::Real fs;                 // FS: Mass-weighted fall speed for snow
            amrex::Real fns;                // FNS: Number-weighted fall speed for snow
            amrex::Real falouts;            // FALOUTS: Fallout rate for snow mass
            amrex::Real faloutns;           // FALOUTNS: Fallout rate for snow number
            amrex::Real faloutg;            // FALOUTG: Fallout rate for graupel mass
            amrex::Real faloutng;           // FALOUTNG: Fallout rate for graupel number
            amrex::Real faltnds;            // FALTNDS: Mass-weighted fallout tendency for snow
            amrex::Real faltndns;           // FALTNDNS: Number-weighted fallout tendency for snow
            amrex::Real unr;                // UNR: Number-weighted terminal velocity for rain
            amrex::Real faltndg;            // FALTNDG: Mass-weighted fallout tendency for graupel
            amrex::Real faltndng;           // FALTNDNG: Number-weighted fallout tendency for graupel
            amrex::Real dumc;               // DUMC: Dummy variable for cloud
            amrex::Real dumfnc;             // DUMFNC: Dummy fall speed number for cloud
            amrex::Real unc;                // UNC: Number-weighted terminal velocity for cloud
            amrex::Real umc;                // UMC: Mass-weighted terminal velocity for cloud
            amrex::Real ung;                // UNG: Number-weighted terminal velocity for graupel
            amrex::Real umg;                // UMG: Mass-weighted terminal velocity for graupel
            amrex::Real fc;                 // FC: Mass-weighted fall speed for cloud
            amrex::Real faloutc;            // FALOUTC: Fallout rate for cloud mass
            amrex::Real faloutnc;           // FALOUTNC: Fallout rate for cloud number
            amrex::Real faltndc;            // FALTNDC: Mass-weighted fallout tendency for cloud
            amrex::Real faltndnc;           // FALTNDNC: Number-weighted fallout tendency for cloud
            amrex::Real fnc;                // FNC: Number-weighted fall speed for cloud
            amrex::Real dumfnr;             // DUMFNR: Dummy fall speed number for rain
            amrex::Real faloutnr;           // FALOUTNR: Fallout rate for rain number
            amrex::Real faltndnr;           // FALTNDNR: Number-weighted fallout tendency for rain
            amrex::Real fnr;                // FNR: Number-weighted fall speed for rain

            // Fall-speed parameter 'A' with air density correction
            amrex::Real ain;                // AIN: Fall-speed parameter 'A' with air density correction for ice
            amrex::Real arn;                // ARN: Fall-speed parameter 'A' with air density correction for rain
            amrex::Real asn;                // ASN: Fall-speed parameter 'A' with air density correction for snow
            amrex::Real acn;                // ACN: Fall-speed parameter 'A' with air density correction for cloud
            amrex::Real agn;                // AGN: Fall-speed parameter 'A' with air density correction for graupel

            // Dummy variables
            amrex::Real dum;                // DUM: General dummy variable
            amrex::Real dum1;               // DUM1: General dummy variable
            amrex::Real dum2;               // DUM2: General dummy variable
            amrex::Real dumt;               // DUMT: Dummy variable for temperature
            amrex::Real dumqv;              // DUMQV: Dummy variable for water vapor
            amrex::Real dumqss;             // DUMQSS: Dummy saturation mixing ratio
            amrex::Real dumqsi;             // DUMQSI: Dummy ice saturation mixing ratio
            amrex::Real dums;               // DUMS: General dummy variable

            // Prognostic supersaturation
            amrex::Real dqsdt;              // DQSDT: Change of saturation mixing ratio with temperature
            amrex::Real dqsidt;             // DQSIDT: Change in ice saturation mixing ratio with temperature
            amrex::Real epsi;               // EPSI: 1/phase relaxation time (see M2005), ice
            amrex::Real epss;               // EPSS: 1/phase relaxation time (see M2005), snow
            amrex::Real epsr;               // EPSR: 1/phase relaxation time (see M2005), rain
            amrex::Real epsg;               // EPSG: 1/phase relaxation time (see M2005), graupel

            // Droplet activation variables
            amrex::Real tauc;               // TAUC: Phase relaxation time (see M2005), droplets
            amrex::Real taur;               // TAUR: Phase relaxation time (see M2005), rain
            amrex::Real taui;               // TAUI: Phase relaxation time (see M2005), cloud ice
            amrex::Real taus;               // TAUS: Phase relaxation time (see M2005), snow
            amrex::Real taug;               // TAUG: Phase relaxation time (see M2005), graupel
            amrex::Real dumact;             // DUMACT: Dummy variable for activation
            amrex::Real dum3;               // DUM3: General dummy variable

            // Counting/index variables
            int k_local=k;                  // K: Vertical level index
            int nstep;                      // NSTEP: Timestep counter
            int n;                          // N: General index variable
            int ltrue;                      // LTRUE: SWITCH = 0: NO HYDROMETEORS IN COLUMN, = 1: HYDROMETEORS IN COLUMN

            // Droplet activation/freezing aerosol
            amrex::Real ct;                 // CT: Droplet activation parameter
            amrex::Real temp1;              // TEMP1: Dummy temperature
            amrex::Real sat1;               // SAT1: Dummy saturation
            amrex::Real sigvl;              // SIGVL: Surface tension liquid/vapor
            amrex::Real kel;                // KEL: Kelvin parameter
            amrex::Real kc2;                // KC2: Total ice nucleation rate
            amrex::Real cry;                // CRY: Aerosol activation parameter
            amrex::Real kry;                // KRY: Aerosol activation parameter

            // More working/dummy variables
            amrex::Real dumqi;              // DUMQI: Dummy variable for ice mixing ratio
            amrex::Real dumni;              // DUMNI: Dummy variable for ice number concentration
            amrex::Real dc0;                // DC0: Characteristic diameter for cloud droplets
            amrex::Real ds0;                // DS0: Characteristic diameter for snow
            amrex::Real dg0;                // DG0: Characteristic diameter for graupel
            amrex::Real dumqc;              // DUMQC: Dummy variable for cloud water mixing ratio
            amrex::Real dumqr;              // DUMQR: Dummy variable for rain mixing ratio
            amrex::Real ratio;              // RATIO: General ratio variable
            amrex::Real sum_dep;            // SUM_DEP: Sum of deposition/sublimation
            amrex::Real fudgef;             // FUDGEF: Adjustment factor

            // Effective vertical velocity (M/S)
            amrex::Real wef;                // WEF: Effective vertical velocity

            // Working parameters for ice nucleation
            amrex::Real anuc;               // ANUC: Ice nucleation parameter
            amrex::Real bnuc;               // BNUC: Ice nucleation parameter

            // Working parameters for aerosol activation
            amrex::Real aact;               // AACT: Aerosol activation parameter
            amrex::Real gamm;               // GAMM: Parameter for aerosol activation
            amrex::Real gg;                 // GG: Parameter for aerosol activation
            amrex::Real psi;                // PSI: Parameter for aerosol activation
            amrex::Real eta1;               // ETA1: Parameter for aerosol activation
            amrex::Real eta2;               // ETA2: Parameter for aerosol activation
            amrex::Real sm1;                // SM1: Parameter for aerosol activation
            amrex::Real sm2;                // SM2: Parameter for aerosol activation
            amrex::Real smax;               // SMAX: Maximum supersaturation
            amrex::Real uu1;                // UU1: Parameter for aerosol activation
            amrex::Real uu2;                // UU2: Parameter for aerosol activation
            amrex::Real alpha;              // ALPHA: Parameter for aerosol activation

            // Dummy size distribution parameters
            amrex::Real dlams;              // DLAMS: Dummy slope parameter for snow
            amrex::Real dlamr;              // DLAMR: Dummy slope parameter for rain
            amrex::Real dlami;              // DLAMI: Dummy slope parameter for ice
            amrex::Real dlamc;              // DLAMC: Dummy slope parameter for cloud
            amrex::Real dlamg;              // DLAMG: Dummy slope parameter for graupel
            amrex::Real lammax;             // LAMMAX: Maximum value for slope parameter
            amrex::Real lammin;             // LAMMIN: Minimum value for slope parameter

            int idrop;                      // IDROP: Switch for droplet activation scheme

            // For WRF-CHEM
            amrex::Real c2prec;             // C2PREC: Cloud to precipitation conversion
            amrex::Real csed;               // CSED: Cloud sedimentation
            amrex::Real ised;               // ISED: Ice sedimentation
            amrex::Real ssed;               // SSED: Snow sedimentation
            amrex::Real gsed;               // GSED: Graupel sedimentation
            amrex::Real rsed;               // RSED: Rain sedimentation
            amrex::Real tqimelt;            // tqimelt: Melting of cloud ice (tendency)
            // NC3DTEN LOCAL ARRAY INITIALIZED
            nc3dten = 0.0;

            // INITIALIZE VARIABLES FOR WRF-CHEM OUTPUT TO ZERO
            c2prec = 0.0;
            csed = 0.0;
            ised = 0.0;
            ssed = 0.0;
            gsed = 0.0;
            rsed = 0.0;

            // LATENT HEAT OF VAPORIZATION
            xxlv = 3.1484E6 - 2370.0 * t3d;
            // LATENT HEAT OF SUBLIMATION
            xxls = 3.15E6 - 2370.0 * t3d + 0.3337E6;

            // Assuming CP is a constant defined elsewhere (specific heat of dry air at constant pressure)
            const amrex::Real CP = 1004.5; // J/kg/K
            cpm = CP * (1.0 + 0.887 * qv3d);

            // SATURATION VAPOR PRESSURE AND MIXING RATIO
            // hm, add fix for low pressure, 5/12/10
            // Assuming POLYSVP is defined elsewhere
            evs = std::min(0.99 * pres, calc_saturation_vapor_pressure(t3d, 0));  // PA
            eis = std::min(0.99 * pres, calc_saturation_vapor_pressure(t3d, 1));  // PA
#if 0
            if ((i >= 86 && i <= 101 && j >= 0 && j <= 3 && k >= 8 && k <= 23 &&
                 (i-86)%2 == 0 && j%2 == 0 && (k-8)%2 == 0) ||
                (i == 92 && j == 0 && k == 17)) {
              FILE *file = fopen("output_cpp.txt", "a");
              fprintf(file, "%5d %5d %5d %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n",
                      i, j, k, xxlv, xxls, cpm, evs, eis, t3d, CP, qv3d, pres);
            }
#endif
         });
          fclose(file);
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

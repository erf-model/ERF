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
        bool use_dev_cpp=false;
        if(use_dev_cpp) {
        this->Cloud(sc);
        this->IceFall(sc);
        this->Precip(sc);
        this->PrecipFall(sc);
        } else {
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
          [[maybe_unused]] auto const& nc_arr = mic_fab_vars[MicVar_Morr::nc]->array(mfi);
          auto const& ns_arr = mic_fab_vars[MicVar_Morr::ns]->array(mfi);
          auto const& nr_arr = mic_fab_vars[MicVar_Morr::nr]->array(mfi);
          auto const& ng_arr = mic_fab_vars[MicVar_Morr::ng]->array(mfi);
          auto const& rho_arr = mic_fab_vars[MicVar_Morr::rho]->array(mfi);
          auto const& pres_arr = mic_fab_vars[MicVar_Morr::pres]->array(mfi);
          [[maybe_unused]] auto const& tabs_arr = mic_fab_vars[MicVar_Morr::tabs]->array(mfi);
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
          [[maybe_unused]] int m_ibase = 2;   // Cloud base activation option
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
#if 0
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
#endif
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
          amrex::Box boxD(box); boxD.makeSlab(2,0);
          bool run_morr_cpp = true;
          bool use_morr_cpp_answer = false;
          bool run_morr_fort = !use_morr_cpp_answer;
          if(run_morr_cpp) {

            // Create dummy arrays for tendencies
            amrex::FArrayBox qc3dten_fab(grown_box, 1);   // CLOUD WATER MIXING RATIO TENDENCY
            amrex::FArrayBox qi3dten_fab(grown_box, 1);   // CLOUD ICE MIXING RATIO TENDENCY  
            amrex::FArrayBox qni3dten_fab(grown_box, 1);  // SNOW MIXING RATIO TENDENCY
            amrex::FArrayBox qr3dten_fab(grown_box, 1);   // RAIN MIXING RATIO TENDENCY
            amrex::FArrayBox ni3dten_fab(grown_box, 1);   // CLOUD ICE NUMBER CONCENTRATION
            amrex::FArrayBox ns3dten_fab(grown_box, 1);   // SNOW NUMBER CONCENTRATION
            amrex::FArrayBox nr3dten_fab(grown_box, 1);   // RAIN NUMBER CONCENTRATION

            // Get array references
            auto const& qc3dten = qc3dten_fab.array();    // CLOUD WATER MIXING RATIO TENDENCY (KG/KG/S)
            auto const& qi3dten = qi3dten_fab.array();    // CLOUD ICE MIXING RATIO TENDENCY (KG/KG/S)
            auto const& qni3dten = qni3dten_fab.array();  // SNOW MIXING RATIO TENDENCY (KG/KG/S)
            auto const& qr3dten = qr3dten_fab.array();    // RAIN MIXING RATIO TENDENCY (KG/KG/S)
            auto const& ni3dten = ni3dten_fab.array();    // CLOUD ICE NUMBER CONCENTRATION (1/KG/S)
            auto const& ns3dten = ns3dten_fab.array();    // SNOW NUMBER CONCENTRATION (1/KG/S)
            auto const& nr3dten = nr3dten_fab.array();    // RAIN NUMBER CONCENTRATION (1/KG/S)

            // Initialize tendencies to zero (no precipitation implemented yet)
            qc3dten_fab.setVal(0.0);
            qi3dten_fab.setVal(0.0);
            qni3dten_fab.setVal(0.0);
            qr3dten_fab.setVal(0.0);
            ni3dten_fab.setVal(0.0);
            ns3dten_fab.setVal(0.0);
            nr3dten_fab.setVal(0.0);

            // Create arrays for mixing ratios and number concentrations
            amrex::FArrayBox qc3d_fab(grown_box, 1);    // CLOUD WATER MIXING RATIO
            amrex::FArrayBox qi3d_fab(grown_box, 1);    // CLOUD ICE MIXING RATIO
            amrex::FArrayBox qni3d_fab(grown_box, 1);   // SNOW MIXING RATIO
            amrex::FArrayBox qr3d_fab(grown_box, 1);    // RAIN MIXING RATIO
            amrex::FArrayBox ni3d_fab(grown_box, 1);    // CLOUD ICE NUMBER CONCENTRATION
            amrex::FArrayBox ns3d_fab(grown_box, 1);    // SNOW NUMBER CONCENTRATION
            amrex::FArrayBox nr3d_fab(grown_box, 1);    // RAIN NUMBER CONCENTRATION

            // Get array references
            auto const& qc3d = qc3d_fab.array();        // CLOUD WATER MIXING RATIO (KG/KG)
            auto const& qi3d = qi3d_fab.array();        // CLOUD ICE MIXING RATIO (KG/KG)
            auto const& qni3d = qni3d_fab.array();      // SNOW MIXING RATIO (KG/KG)
            auto const& qr3d = qr3d_fab.array();        // RAIN MIXING RATIO (KG/KG)
            auto const& ni3d = ni3d_fab.array();        // CLOUD ICE NUMBER CONCENTRATION (1/KG)
            auto const& ns3d = ns3d_fab.array();        // SNOW NUMBER CONCENTRATION (1/KG)
            auto const& nr3d = nr3d_fab.array();        // RAIN NUMBER CONCENTRATION (1/KG)

            // Create arrays for temperature, vapor, and pressure variables
            amrex::FArrayBox t3dten_fab(grown_box, 1);    // TEMPERATURE TENDENCY
            amrex::FArrayBox qv3dten_fab(grown_box, 1);   // WATER VAPOR MIXING RATIO TENDENCY
            amrex::FArrayBox t3d_fab(grown_box, 1);       // TEMPERATURE
            amrex::FArrayBox qv3d_fab(grown_box, 1);      // WATER VAPOR MIXING RATIO
            amrex::FArrayBox pres_fab(grown_box, 1);      // ATMOSPHERIC PRESSURE
            amrex::FArrayBox dzq_fab(grown_box, 1);       // DIFFERENCE IN HEIGHT ACROSS LEVEL
            amrex::FArrayBox w3d_fab(grown_box, 1);       // GRID-SCALE VERTICAL VELOCITY

            // WRF-chem variables
            amrex::FArrayBox nc3d_fab(grown_box, 1);      // CLOUD DROPLET NUMBER CONCENTRATION
            amrex::FArrayBox nc3dten_fab(grown_box, 1);   // CLOUD DROPLET NUMBER CONCENTRATION TENDENCY

            // Graupel variables
            amrex::FArrayBox qg3dten_fab(grown_box, 1);   // GRAUPEL MIX RATIO TENDENCY
            amrex::FArrayBox ng3dten_fab(grown_box, 1);   // GRAUPEL NUMB CONC TENDENCY
            amrex::FArrayBox qg3d_fab(grown_box, 1);      // GRAUPEL MIX RATIO
            amrex::FArrayBox ng3d_fab(grown_box, 1);      // GRAUPEL NUMBER CONC

            // Sedimentation tendencies
            amrex::FArrayBox qgsten_fab(grown_box, 1);    // GRAUPEL SED TEND
            amrex::FArrayBox qrsten_fab(grown_box, 1);    // RAIN SED TEND
            amrex::FArrayBox qisten_fab(grown_box, 1);    // CLOUD ICE SED TEND
            amrex::FArrayBox qnisten_fab(grown_box, 1);   // SNOW SED TEND
            amrex::FArrayBox qcsten_fab(grown_box, 1);    // CLOUD WAT SED TEND

            // Cumulus tendencies
            amrex::FArrayBox qrcu1d_fab(grown_box, 1);    // RAIN FROM CUMULUS PARAMETERIZATION
            amrex::FArrayBox qscu1d_fab(grown_box, 1);    // SNOW FROM CUMULUS PARAMETERIZATION
            amrex::FArrayBox qicu1d_fab(grown_box, 1);    // ICE FROM CUMULUS PARAMETERIZATION

            // Get array references
            auto const& t3dten = t3dten_fab.array();      // TEMPERATURE TENDENCY (K/S)
            auto const& qv3dten = qv3dten_fab.array();    // WATER VAPOR MIXING RATIO TENDENCY (KG/KG/S)
            auto const& t3d = t3d_fab.array();            // TEMPERATURE (K)
            auto const& qv3d = qv3d_fab.array();          // WATER VAPOR MIXING RATIO (KG/KG)
            auto const& pres = pres_fab.array();          // ATMOSPHERIC PRESSURE (PA)
            auto const& dzq = dzq_fab.array();            // DIFFERENCE IN HEIGHT ACROSS LEVEL (m)
            auto const& w3d = w3d_fab.array();            // GRID-SCALE VERTICAL VELOCITY (M/S)

            // WRF-chem variables
            auto const& nc3d = nc3d_fab.array();          // CLOUD DROPLET NUMBER CONCENTRATION
            auto const& nc3dten = nc3dten_fab.array();    // CLOUD DROPLET NUMBER CONCENTRATION TENDENCY

            // Graupel variables
            auto const& qg3dten = qg3dten_fab.array();    // GRAUPEL MIX RATIO TENDENCY (KG/KG/S)
            auto const& ng3dten = ng3dten_fab.array();    // GRAUPEL NUMB CONC TENDENCY (1/KG/S)
            auto const& qg3d = qg3d_fab.array();          // GRAUPEL MIX RATIO (KG/KG)
            auto const& ng3d = ng3d_fab.array();          // GRAUPEL NUMBER CONC (1/KG)

            // Sedimentation tendencies
            auto const& qgsten = qgsten_fab.array();      // GRAUPEL SED TEND (KG/KG/S)
            auto const& qrsten = qrsten_fab.array();      // RAIN SED TEND (KG/KG/S)
            auto const& qisten = qisten_fab.array();      // CLOUD ICE SED TEND (KG/KG/S)
            auto const& qnisten = qnisten_fab.array();    // SNOW SED TEND (KG/KG/S)
            auto const& qcsten = qcsten_fab.array();      // CLOUD WAT SED TEND (KG/KG/S)

            // Cumulus tendencies
            auto const& qrcu1d = qrcu1d_fab.array();      // RAIN FROM CUMULUS PARAMETERIZATION
            auto const& qscu1d = qscu1d_fab.array();      // SNOW FROM CUMULUS PARAMETERIZATION
            auto const& qicu1d = qicu1d_fab.array();      // ICE FROM CUMULUS PARAMETERIZATION

            // Initialize tendency arrays to zero
            t3dten_fab.setVal(0.0);
            qv3dten_fab.setVal(0.0);
            nc3dten_fab.setVal(0.0);
            qg3dten_fab.setVal(0.0);
            ng3dten_fab.setVal(0.0);
            qgsten_fab.setVal(0.0);
            qrsten_fab.setVal(0.0);
            qisten_fab.setVal(0.0);
            qnisten_fab.setVal(0.0);
            qcsten_fab.setVal(0.0);

            // Create arrays for precipitation rates
            amrex::FArrayBox precrt_fab(grown_box, 1);    // TOTAL PRECIP PER TIME STEP
            amrex::FArrayBox snowrt_fab(grown_box, 1);    // SNOW PER TIME STEP
            amrex::FArrayBox snowprt_fab(grown_box, 1);   // TOTAL CLOUD ICE PLUS SNOW PER TIME STEP
            amrex::FArrayBox grplprt_fab(grown_box, 1);   // TOTAL GRAUPEL PER TIME STEP

            // Create arrays for effective radii
            amrex::FArrayBox effc_fab(grown_box, 1);      // DROPLET EFFECTIVE RADIUS
            amrex::FArrayBox effi_fab(grown_box, 1);      // CLOUD ICE EFFECTIVE RADIUS
            amrex::FArrayBox effs_fab(grown_box, 1);      // SNOW EFFECTIVE RADIUS
            amrex::FArrayBox effr_fab(grown_box, 1);      // RAIN EFFECTIVE RADIUS
            amrex::FArrayBox effg_fab(grown_box, 1);      // GRAUPEL EFFECTIVE RADIUS

            // Get array references for precipitation rates
            auto const& precrt = precrt_fab.array();      // TOTAL PRECIP PER TIME STEP (mm)
            auto const& snowrt = snowrt_fab.array();      // SNOW PER TIME STEP (mm)
            auto const& snowprt = snowprt_fab.array();    // TOTAL CLOUD ICE PLUS SNOW PER TIME STEP (mm)
            auto const& grplprt = grplprt_fab.array();    // TOTAL GRAUPEL PER TIME STEP (mm)

            // Get array references for effective radii
            auto const& effc = effc_fab.array();          // DROPLET EFFECTIVE RADIUS (MICRON)
            auto const& effi = effi_fab.array();          // CLOUD ICE EFFECTIVE RADIUS (MICRON)
            auto const& effs = effs_fab.array();          // SNOW EFFECTIVE RADIUS (MICRON)
            auto const& effr = effr_fab.array();          // RAIN EFFECTIVE RADIUS (MICRON)
            auto const& effg = effg_fab.array();          // GRAUPEL EFFECTIVE RADIUS (MICRON)

            // Initialize these arrays to zero (they will be computed later)
            precrt_fab.setVal(0.0);
            snowrt_fab.setVal(0.0);
            snowprt_fab.setVal(0.0);
            grplprt_fab.setVal(0.0);
            effc_fab.setVal(0.0);
            effi_fab.setVal(0.0);
            effs_fab.setVal(0.0);
            effr_fab.setVal(0.0);
            effg_fab.setVal(0.0);

            // Create FArrayBoxes for scalar variables
            amrex::FArrayBox rho_fab(grown_box, 1);
            amrex::FArrayBox mu_fab(grown_box, 1);
            amrex::FArrayBox ain_fab(grown_box, 1);
            amrex::FArrayBox arn_fab(grown_box, 1);
            amrex::FArrayBox asn_fab(grown_box, 1);
            amrex::FArrayBox acn_fab(grown_box, 1);
            amrex::FArrayBox agn_fab(grown_box, 1);

            // Get Array4 views
            auto const& rho = rho_fab.array();
            auto const& mu = mu_fab.array();
            auto const& ain = ain_fab.array();
            auto const& arn = arn_fab.array();
            auto const& asn = asn_fab.array();
            auto const& acn = acn_fab.array();
            auto const& agn = agn_fab.array();

            // Initialize all values to zero
            rho_fab.setVal(0.0);
            mu_fab.setVal(0.0);
            ain_fab.setVal(0.0);
            arn_fab.setVal(0.0);
            asn_fab.setVal(0.0);
            acn_fab.setVal(0.0);
            agn_fab.setVal(0.0);

            // Create FArrayBoxes for fall speed working variables
            amrex::FArrayBox dumi_fab(grown_box, 1);
            amrex::FArrayBox dumr_fab(grown_box, 1);
            amrex::FArrayBox dumfni_fab(grown_box, 1);
            amrex::FArrayBox dumg_fab(grown_box, 1);
            amrex::FArrayBox dumfng_fab(grown_box, 1);
            amrex::FArrayBox uni_fab(grown_box, 1);
            amrex::FArrayBox umi_fab(grown_box, 1);
            amrex::FArrayBox umr_fab(grown_box, 1);
            amrex::FArrayBox fr_fab(grown_box, 1);
            amrex::FArrayBox fi_fab(grown_box, 1);
            amrex::FArrayBox fni_fab(grown_box, 1);
            amrex::FArrayBox fg_fab(grown_box, 1);
            amrex::FArrayBox fng_fab(grown_box, 1);
            amrex::FArrayBox rgvm_fab(grown_box, 1);
            amrex::FArrayBox faloutr_fab(grown_box, 1);
            amrex::FArrayBox falouti_fab(grown_box, 1);
            amrex::FArrayBox faloutni_fab(grown_box, 1);
            amrex::FArrayBox faltndr_fab(grown_box, 1);
            amrex::FArrayBox faltndi_fab(grown_box, 1);
            amrex::FArrayBox faltndni_fab(grown_box, 1);
            amrex::FArrayBox rho2_fab(grown_box, 1);
            amrex::FArrayBox dumqs_fab(grown_box, 1);
            amrex::FArrayBox dumfns_fab(grown_box, 1);
            amrex::FArrayBox ums_fab(grown_box, 1);
            amrex::FArrayBox uns_fab(grown_box, 1);
            amrex::FArrayBox fs_fab(grown_box, 1);
            amrex::FArrayBox fns_fab(grown_box, 1);
            amrex::FArrayBox falouts_fab(grown_box, 1);
            amrex::FArrayBox faloutns_fab(grown_box, 1);
            amrex::FArrayBox faloutg_fab(grown_box, 1);
            amrex::FArrayBox faloutng_fab(grown_box, 1);
            amrex::FArrayBox faltnds_fab(grown_box, 1);
            amrex::FArrayBox faltndns_fab(grown_box, 1);
            amrex::FArrayBox unr_fab(grown_box, 1);
            amrex::FArrayBox faltndg_fab(grown_box, 1);
            amrex::FArrayBox faltndng_fab(grown_box, 1);
            amrex::FArrayBox dumc_fab(grown_box, 1);
            amrex::FArrayBox dumfnc_fab(grown_box, 1);
            amrex::FArrayBox unc_fab(grown_box, 1);
            amrex::FArrayBox umc_fab(grown_box, 1);
            amrex::FArrayBox ung_fab(grown_box, 1);
            amrex::FArrayBox umg_fab(grown_box, 1);
            amrex::FArrayBox fc_fab(grown_box, 1);
            amrex::FArrayBox faloutc_fab(grown_box, 1);
            amrex::FArrayBox faloutnc_fab(grown_box, 1);
            amrex::FArrayBox faltndc_fab(grown_box, 1);
            amrex::FArrayBox faltndnc_fab(grown_box, 1);
            amrex::FArrayBox fnc_fab(grown_box, 1);
            amrex::FArrayBox dumfnr_fab(grown_box, 1);
            amrex::FArrayBox faloutnr_fab(grown_box, 1);
            amrex::FArrayBox faltndnr_fab(grown_box, 1);
            amrex::FArrayBox fnr_fab(grown_box, 1);
            amrex::FArrayBox dumqi_fab(grown_box, 1);
            amrex::FArrayBox dumni_fab(grown_box, 1);
            amrex::FArrayBox dumqc_fab(grown_box, 1);
            amrex::FArrayBox dumqr_fab(grown_box, 1);
            amrex::FArrayBox ratio_fab(grown_box, 1);
            amrex::FArrayBox sum_dep_fab(grown_box, 1);
            amrex::FArrayBox fudgef_fab(grown_box, 1);
            amrex::FArrayBox dlams_fab(grown_box, 1);
            amrex::FArrayBox dlamr_fab(grown_box, 1);
            amrex::FArrayBox dlami_fab(grown_box, 1);
            amrex::FArrayBox dlamc_fab(grown_box, 1);
            amrex::FArrayBox dlamg_fab(grown_box, 1);

            // Create Array4 references
            auto const& dumi = dumi_fab.array();
            auto const& dumr = dumr_fab.array();
            auto const& dumfni = dumfni_fab.array();
            auto const& dumg = dumg_fab.array();
            auto const& dumfng = dumfng_fab.array();
            auto const& uni = uni_fab.array();
            auto const& umi = umi_fab.array();
            auto const& umr = umr_fab.array();
            auto const& fr = fr_fab.array();
            auto const& fi = fi_fab.array();
            auto const& fni = fni_fab.array();
            auto const& fg = fg_fab.array();
            auto const& fng = fng_fab.array();
            auto const& rgvm = rgvm_fab.array();
            auto const& faloutr = faloutr_fab.array();
            auto const& falouti = falouti_fab.array();
            auto const& faloutni = faloutni_fab.array();
            auto const& faltndr = faltndr_fab.array();
            auto const& faltndi = faltndi_fab.array();
            auto const& faltndni = faltndni_fab.array();
            auto const& rho2 = rho2_fab.array();
            auto const& dumqs = dumqs_fab.array();
            auto const& dumfns = dumfns_fab.array();
            auto const& ums = ums_fab.array();
            auto const& uns = uns_fab.array();
            auto const& fs = fs_fab.array();
            auto const& fns = fns_fab.array();
            auto const& falouts = falouts_fab.array();
            auto const& faloutns = faloutns_fab.array();
            auto const& faloutg = faloutg_fab.array();
            auto const& faloutng = faloutng_fab.array();
            auto const& faltnds = faltnds_fab.array();
            auto const& faltndns = faltndns_fab.array();
            auto const& unr = unr_fab.array();
            auto const& faltndg = faltndg_fab.array();
            auto const& faltndng = faltndng_fab.array();
            auto const& dumc = dumc_fab.array();
            auto const& dumfnc = dumfnc_fab.array();
            auto const& unc = unc_fab.array();
            auto const& umc = umc_fab.array();
            auto const& ung = ung_fab.array();
            auto const& umg = umg_fab.array();
            auto const& fc = fc_fab.array();
            auto const& faloutc = faloutc_fab.array();
            auto const& faloutnc = faloutnc_fab.array();
            auto const& faltndc = faltndc_fab.array();
            auto const& faltndnc = faltndnc_fab.array();
            auto const& fnc = fnc_fab.array();
            auto const& dumfnr = dumfnr_fab.array();
            auto const& faloutnr = faloutnr_fab.array();
            auto const& faltndnr = faltndnr_fab.array();
            auto const& fnr = fnr_fab.array();
            auto const& dumqi = dumqi_fab.array();
            auto const& dumni = dumni_fab.array();
            auto const& dumqc = dumqc_fab.array();
            auto const& dumqr = dumqr_fab.array();
            auto const& ratio = ratio_fab.array();
            auto const& sum_dep = sum_dep_fab.array();
            auto const& fudgef = fudgef_fab.array();
            auto const& dlams = dlams_fab.array();
            auto const& dlamr = dlamr_fab.array();
            auto const& dlami = dlami_fab.array();
            auto const& dlamc = dlamc_fab.array();
            auto const& dlamg = dlamg_fab.array();

            // Initialize arrays to 0
            dlams_fab.setVal(0.0);
            dlamr_fab.setVal(0.0);
            dlami_fab.setVal(0.0);
            dlamc_fab.setVal(0.0);
            dlamg_fab.setVal(0.0);
            dumi_fab.setVal(0.0);
            dumr_fab.setVal(0.0);
            dumfni_fab.setVal(0.0);
            dumg_fab.setVal(0.0);
            dumfng_fab.setVal(0.0);
            uni_fab.setVal(0.0);
            umi_fab.setVal(0.0);
            umr_fab.setVal(0.0);
            fr_fab.setVal(0.0);
            fi_fab.setVal(0.0);
            fni_fab.setVal(0.0);
            fg_fab.setVal(0.0);
            fng_fab.setVal(0.0);
            rgvm_fab.setVal(0.0);
            faloutr_fab.setVal(0.0);
            falouti_fab.setVal(0.0);
            faloutni_fab.setVal(0.0);
            faltndr_fab.setVal(0.0);
            faltndi_fab.setVal(0.0);
            faltndni_fab.setVal(0.0);
            rho2_fab.setVal(0.0);
            dumqs_fab.setVal(0.0);
            dumfns_fab.setVal(0.0);
            ums_fab.setVal(0.0);
            uns_fab.setVal(0.0);
            fs_fab.setVal(0.0);
            fns_fab.setVal(0.0);
            falouts_fab.setVal(0.0);
            faloutns_fab.setVal(0.0);
            faloutg_fab.setVal(0.0);
            faloutng_fab.setVal(0.0);
            faltnds_fab.setVal(0.0);
            faltndns_fab.setVal(0.0);
            unr_fab.setVal(0.0);
            faltndg_fab.setVal(0.0);
            faltndng_fab.setVal(0.0);
            dumc_fab.setVal(0.0);
            dumfnc_fab.setVal(0.0);
            unc_fab.setVal(0.0);
            umc_fab.setVal(0.0);
            ung_fab.setVal(0.0);
            umg_fab.setVal(0.0);
            fc_fab.setVal(0.0);
            faloutc_fab.setVal(0.0);
            faloutnc_fab.setVal(0.0);
            faltndc_fab.setVal(0.0);
            faltndnc_fab.setVal(0.0);
            fnc_fab.setVal(0.0);
            dumfnr_fab.setVal(0.0);
            faloutnr_fab.setVal(0.0);
            faltndnr_fab.setVal(0.0);
            fnr_fab.setVal(0.0);
            dumqi_fab.setVal(0.0);
            dumni_fab.setVal(0.0);
            dumqc_fab.setVal(0.0);
            dumqr_fab.setVal(0.0);
            ratio_fab.setVal(0.0);
            sum_dep_fab.setVal(0.0);
            fudgef_fab.setVal(0.0);

            FILE *file = fopen("output_cpp.txt", "a");
          ////////////////////////////////////////////////////////////
          // ParallelFor for testing partial C++ implementation
          // NOTE: Currently all Array4 values are copied to locals
          //       This means we're not updating or outputing anything
          ////////////////////////////////////////////////////////////
          ParallelFor( boxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
         {
           int ltrue=0;                      // LTRUE: SWITCH = 0: NO HYDROMETEORS IN COLUMN, = 1: HYDROMETEORS IN COLUMN
           int nstep;                        // NSTEP: Timestep counter
           int iinum;                      // iinum: Integer control variable

           for(int k=klo; k<=khi; k++) {
            // Tendencies and mixing ratios
            qc3d(i,j,k) = qcl_arr(i,j,k);   // CLOUD WATER MIXING RATIO
            qi3d(i,j,k) = qci_arr(i,j,k);   // CLOUD ICE MIXING RATIO
            qni3d(i,j,k) = qps_arr(i,j,k);  // SNOW MIXING RATIO
            qr3d(i,j,k) = qpr_arr(i,j,k);   // RAIN MIXING RATIO 
            ni3d(i,j,k) = ni_arr(i,j,k);    // CLOUD ICE NUMBER CONCENTRATION
            ns3d(i,j,k) = ns_arr(i,j,k);    // SNOW NUMBER CONCENTRATION
            nr3d(i,j,k) = nr_arr(i,j,k);    // RAIN NUMBER CONCENTRATION
            nc3d(i,j,k) = nc_arr(i,j,k);    // RAIN NUMBER CONCENTRATION

            t3d(i,j,k) = theta_arr(i,j,k) * pii_arr(i,j,k);  // TEMPERATURE
            qv3d(i,j,k) = qv_arr(i,j,k);                     // WATER VAPOR MIXING RATIO
            pres(i,j,k) = pres_arr(i,j,k);                   // ATMOSPHERIC PRESSURE
            dzq(i,j,k) = dz_arr(i,j,k);                      // DIFFERENCE IN HEIGHT ACROSS LEVEL
            w3d(i,j,k) = w_arr(i,j,k);                       // GRID-SCALE VERTICAL VELOCITY
            qg3d(i,j,k) = qpg_arr(i,j,k);                    // GRAUPEL MIX RATIO
            ng3d(i,j,k) = ng_arr(i,j,k);                     // GRAUPEL NUMBER CONC
            qrcu1d(i,j,k) = qrcuten_arr(i,j,k);              // RAIN FROM CUMULUS PARAMETERIZATION
            qscu1d(i,j,k) = qscuten_arr(i,j,k);              // SNOW FROM CUMULUS PARAMETERIZATION
            qicu1d(i,j,k) = qicuten_arr(i,j,k);              // ICE FROM CUMULUS PARAMETERIZATION

            // Model input parameters
            //amrex::Real dt;                 // DT: MODEL TIME STEP (SEC)
            //amrex::Real lami;               // LAMI: Slope parameter for cloud ice (m^-1)
#if 1
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
#endif
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
            amrex::Real sc_schmidt;         // SC: Schmidt number
            amrex::Real xlf;                // XLF: Latent heat of freezing
            amrex::Real ab;                 // AB: Correction to condensation rate due to latent heating
            amrex::Real abi;                // ABI: Correction to deposition rate due to latent heating
#if 0
            // Time-varying microphysics parameters
            amrex::Real dap;                // DAP: Diffusivity of aerosol
            amrex::Real nacnt;              // NACNT: Number of contact nuclei
            amrex::Real fmult;              // FMULT: Temperature-dependent parameter for rime-splintering
            amrex::Real coffi;              // COFFI: Ice autoconversion parameter

#endif

            // Dummy variables
            amrex::Real dum;                // DUM: General dummy variable
#if 1
            amrex::Real dum1;               // DUM1: General dummy variable
            amrex::Real dum2;               // DUM2: General dummy variable
            amrex::Real dumt;               // DUMT: Dummy variable for temperature
            amrex::Real dumqv;              // DUMQV: Dummy variable for water vapor
            amrex::Real dumqss;             // DUMQSS: Dummy saturation mixing ratio
            amrex::Real dumqsi;             // DUMQSI: Dummy ice saturation mixing ratio
            amrex::Real dums;               // DUMS: General dummy variable
#endif
            // Prognostic supersaturation
            amrex::Real dqsdt;              // DQSDT: Change of saturation mixing ratio with temperature
            amrex::Real dqsidt;             // DQSIDT: Change in ice saturation mixing ratio with temperature
#if 1
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
            int n;                          // N: General index variable

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
            amrex::Real di0;                // DC0: Characteristic diameter for ice
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

            int idrop;                      // IDROP: Switch for droplet activation scheme
#endif
            // For WRF-CHEM
            amrex::Real c2prec;             // C2PREC: Cloud to precipitation conversion
            amrex::Real csed;               // CSED: Cloud sedimentation
            amrex::Real ised;               // ISED: Ice sedimentation
            amrex::Real ssed;               // SSED: Snow sedimentation
            amrex::Real gsed;               // GSED: Graupel sedimentation
            amrex::Real rsed;               // RSED: Rain sedimentation
            [[maybe_unused]] amrex::Real tqimelt;            // tqimelt: Melting of cloud ice (tendency)

            // NC3DTEN LOCAL ARRAY INITIALIZED
            nc3dten(i,j,k) = 0.0;

            // INITIALIZE VARIABLES FOR WRF-CHEM OUTPUT TO ZERO
            c2prec = 0.0;
            csed = 0.0;
            ised = 0.0;
            ssed = 0.0;
            gsed = 0.0;
            rsed = 0.0;

            // LATENT HEAT OF VAPORIZATION
            xxlv = 3.1484E6 - 2370.0 * t3d(i,j,k);
            // LATENT HEAT OF SUBLIMATION
            xxls = 3.15E6 - 2370.0 * t3d(i,j,k) + 0.3337E6;

            // Assuming CP is a constant defined elsewhere (specific heat of dry air at constant pressure)
            const amrex::Real CP = 1004.5; // J/kg/K
            cpm = CP * (1.0 + 0.887 * qv3d(i,j,k));

            // SATURATION VAPOR PRESSURE AND MIXING RATIO
            // hm, add fix for low pressure, 5/12/10
            // Assuming POLYSVP is defined elsewhere
            evs = std::min(0.99 * pres(i,j,k), calc_saturation_vapor_pressure(t3d(i,j,k), 0));  // PA
            eis = std::min(0.99 * pres(i,j,k), calc_saturation_vapor_pressure(t3d(i,j,k), 1));  // PA
#if 0
            if ((i >= 86 && i <= 101 && j >= 0 && j <= 3 && k >= 8 && k <= 23 &&
                 (i-86)%2 == 0 && j%2 == 0 && (k-8)%2 == 0) ||
                (i == 92 && j == 0 && k == 17)) {
              FILE *file = fopen("output_cpp.txt", "a");
              fprintf(file, "%5d %5d %5d %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n",
                      i, j, k, xxlv, xxls, cpm, evs, eis, t3d(i,j,k), CP, qv3d(i,j,k), pres(i,j,k));
            }
#endif
            // MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING
            if (eis > evs) {
              eis = evs; // temporary update: adjust ice saturation pressure
            }

            // SATURATION MIXING RATIOS
            qvs = m_ep_2 * evs / (pres(i,j,k) - evs); // budget equation: calculate water saturation mixing ratio
            qvi = m_ep_2 * eis / (pres(i,j,k) - eis); // budget equation: calculate ice saturation mixing ratio

            // SATURATION RATIOS
            qvqvs = qv3d(i,j,k) / qvs; // budget equation: calculate water saturation ratio
            qvqvsi = qv3d(i,j,k) / qvi; // budget equation: calculate ice saturation ratio

            // AIR DENSITY
            rho(i,j,k) = pres(i,j,k) / (m_R * t3d(i,j,k)); // budget equation: calculate air density

            ds0 = 3.0;       // Size distribution parameter for snow
            di0 = 3.0;       // Size distribution parameter for cloud ice
            dg0 = 3.0;       // Size distribution parameter for graupel
            const double CI = 800.0;     // Mass-diameter relationship parameter for cloud ice
            // ADD NUMBER CONCENTRATION DUE TO CUMULUS TENDENCY
            // ASSUME N0 ASSOCIATED WITH CUMULUS PARAM RAIN IS 10^7 M^-4
            // ASSUME N0 ASSOCIATED WITH CUMULUS PARAM SNOW IS 2 X 10^7 M^-4
            // FOR DETRAINED CLOUD ICE, ASSUME MEAN VOLUME DIAM OF 80 MICRON
            if (qrcu1d(i,j,k) >= 1.0e-10) {
              dum = 1.8e5 * std::pow(qrcu1d(i,j,k) * dt / (m_pi * m_rhow * std::pow(rho(i,j,k), 3)), 0.25); // rate equation: calculate rain number concentration from cumulus
              nr3d(i,j,k) += dum; // budget equation: update rain number concentration
            }
            if (qscu1d(i,j,k) >= 1.0e-10) {
              dum = 3.e5 * std::pow(qscu1d(i,j,k) * dt / (m_cons1 * std::pow(rho(i,j,k), 3)), 1.0 / (ds0 + 1.0)); // rate equation: calculate snow number concentration from cumulus
              ns3d(i,j,k) += dum; // budget equation: update snow number concentration
            }
            if (qicu1d(i,j,k) >= 1.0e-10) {
              dum = qicu1d(i,j,k) * dt / (CI * std::pow(80.0e-6, di0)); // rate equation: calculate cloud ice number concentration from cumulus
              ni3d(i,j,k) += dum; // budget equation: update cloud ice number concentration
            }

            // AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
            // hm modify 7/0/09 change limit to 1.e-8
            if (qvqvs < 0.9) {
              if (qr3d(i,j,k) < 1.0e-8) {
                qv3d(i,j,k) += qr3d(i,j,k); // budget equation: transfer rain to vapor
                t3d(i,j,k) -= qr3d(i,j,k) * xxlv / cpm; // budget equation: adjust temperature
                qr3d(i,j,k) = 0.0; // temporary update: set rain to zero
              }
              if (qc3d(i,j,k) < 1.0e-8) {
                qv3d(i,j,k) += qc3d(i,j,k); // budget equation: transfer cloud water to vapor
                t3d(i,j,k) -= qc3d(i,j,k) * xxlv / cpm; // budget equation: adjust temperature
                qc3d(i,j,k) = 0.0; // temporary update: set cloud water to zero
              }
            }

            if (qvqvsi < 0.9) {
              if (qi3d(i,j,k) < 1.0e-8) {
                qv3d(i,j,k) += qi3d(i,j,k); // budget equation: transfer cloud ice to vapor
                t3d(i,j,k) -= qi3d(i,j,k) * xxls / cpm; // budget equation: adjust temperature
                qi3d(i,j,k) = 0.0; // temporary update: set cloud ice to zero
              }
              if (qni3d(i,j,k) < 1.0e-8) {
                qv3d(i,j,k) += qni3d(i,j,k); // budget equation: transfer snow to vapor
                t3d(i,j,k) -= qni3d(i,j,k) * xxls / cpm; // budget equation: adjust temperature
                qni3d(i,j,k) = 0.0; // temporary update: set snow to zero
              }
              if (qg3d(i,j,k) < 1.0e-8) {
                qv3d(i,j,k) += qg3d(i,j,k); // budget equation: transfer graupel to vapor
                t3d(i,j,k) -= qg3d(i,j,k) * xxls / cpm; // budget equation: adjust temperature
                qg3d(i,j,k) = 0.0; // temporary update: set graupel to zero
              }
            }
#if 0
            // C++ version
            if ((i >= 86 && i <= 101 && j >= 0 && j <= 3 && k >= 8 && k <= 23 &&
                 (i-86)%2 == 0 && j%2 == 0 && (k-8)%2 == 0) ||
                (i == 92 && j == 0 && k == 17) ||
                (((i == 168 || i == 169 || i == 190 || i == 191) &&
                  (j == 0 || j == 3) &&
                  (k == 0 || k == 1 || k == 126 || k == 127))) ||
                (i == 175 && j == 1 && k == 50) ||
                (i == 180 && j == 2 && k == 75) ||
                (i == 185 && j == 1 && k == 100) ||
                (i == 170 && j == 0 && k == 30) ||
                (i == 188 && j == 3 && k == 60) ||
                (i == 178 && j == 2 && k == 40) ||
                (i == 183 && j == 0 && k == 80) ||
                (i == 173 && j == 3 && k == 110) ||
                (i == 186 && j == 1 && k == 90) ||
                (i == 177 && j == 2 && k == 65)) {
              fprintf(file, "%5d %5d %5d %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n",
                      i, j, k, xxlv, xxls, cpm, evs, eis, t3d(i,j,k), CP, qv3d(i,j,k), pres(i,j,k),
                      qvs, qvi, qvqvs, qvqvsi, rho, qr3d(i,j,k), qc3d(i,j,k), qi3d(i,j,k), qni3d(i,j,k), qg3d(i,j,k), nr3d(i,j,k), ns3d(i,j,k), ni3d(i,j,k));
            }
#endif
            // HEAT OF FUSION
            xlf = xxls - xxlv;

            // IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO
            // Note: QSMALL is not defined in the variable list, so I'll define it
            const amrex::Real QSMALL = m_qsmall;

            if (qc3d(i,j,k) < QSMALL) {
              qc3d(i,j,k) = 0.0;
              nc3d(i,j,k) = 0.0;
              effc(i,j,k) = 0.0;
            }
            if (qr3d(i,j,k) < QSMALL) {
              qr3d(i,j,k) = 0.0;
              nr3d(i,j,k) = 0.0;
              effr(i,j,k) = 0.0;
            }
            if (qi3d(i,j,k) < QSMALL) {
              qi3d(i,j,k) = 0.0;
              ni3d(i,j,k) = 0.0;
              effi(i,j,k) = 0.0;
            }
            if (qni3d(i,j,k) < QSMALL) {
              qni3d(i,j,k) = 0.0;
              ns3d(i,j,k) = 0.0;
              effs(i,j,k) = 0.0;
            }
            if (qg3d(i,j,k) < QSMALL) {
              qg3d(i,j,k) = 0.0;
              ng3d(i,j,k) = 0.0;
              effg(i,j,k) = 0.0;
            }

            // INITIALIZE SEDIMENTATION TENDENCIES FOR MIXING RATIO
            qrsten(i,j,k) = 0.0;  // temporary update: initialize QRSTEN
            qisten(i,j,k) = 0.0;  // temporary update: initialize QISTEN
            qnisten(i,j,k) = 0.0; // temporary update: initialize QNISTEN
            qcsten(i,j,k) = 0.0;  // temporary update: initialize QCSTEN
            qgsten(i,j,k) = 0.0;  // temporary update: initialize QGSTEN

            // MICROPHYSICS PARAMETERS VARYING IN TIME/HEIGHT
            mu(i,j,k) = 1.496e-6 * std::pow(t3d(i,j,k), 1.5) / (t3d(i,j,k) + 120.0); // budget equation: calculate air viscosity

            // Fall speed with density correction (Heymsfield and Benssemer 2006)
            dum = std::pow(m_rhosu / rho(i,j,k), 0.54); // temporary update: calculate density correction factor

            // AA revision 4/1/11: Ikawa and Saito 1991 air-density correction
            ain(i,j,k) = std::pow(m_rhosu / rho(i,j,k), 0.35) * m_ai; // budget equation: calculate ice fall speed parameter
            arn(i,j,k) = dum * m_ar; // budget equation: calculate rain fall speed parameter
            asn(i,j,k) = dum * m_as; // budget equation: calculate snow fall speed parameter

            // AA revision 4/1/11: temperature-dependent Stokes fall speed
            acn(i,j,k) = m_g * m_rhow / (18.0 * mu(i,j,k)); // budget equation: calculate cloud droplet fall speed parameter

            // HM ADD GRAUPEL 8/28/06
            agn(i,j,k) = dum * m_ag; // budget equation: calculate graupel fall speed parameter

            // hm 4/7/09 bug fix, initialize lami to prevent later division by zero
            lami = 0.0; // temporary update: initialize LAMI

            // If there is no cloud/precip water, and if subsaturated, then skip microphysics for this level
            bool skipMicrophysics = false;
            bool skipConcentrations = false;
            bool skipPrecip = true; // set with if statement
            if (qc3d(i,j,k) < QSMALL && qi3d(i,j,k) < QSMALL && qni3d(i,j,k) < QSMALL && qr3d(i,j,k) < QSMALL && qg3d(i,j,k) < QSMALL) {
              if ((t3d(i,j,k) < 273.15 && qvqvsi < 0.999) || (t3d(i,j,k) >= 273.15 && qvqvs < 0.999)) {
                skipMicrophysics = true; // GOTO 200
              }
            }

            if(!skipMicrophysics) {

            // Thermal conductivity for air
              kap = 1.414e3 * mu(i,j,k); // budget equation: calculate thermal conductivity

            // Diffusivity of water vapor
            dv = 8.794e-5 * std::pow(t3d(i,j,k), 1.81) / pres(i,j,k); // budget equation: calculate vapor diffusivity

            // Schmidt number
            sc_schmidt = mu(i,j,k) / (rho(i,j,k) * dv); // budget equation: calculate Schmidt number

            // Psychometric corrections
            // Rate of change sat. mix. ratio with temperature
            dum = (m_Rv * t3d(i,j,k) * t3d(i,j,k)); // temporary update: calculate temperature factor
            dqsdt = xxlv * qvs / dum; // budget equation: calculate DQSDT
            dqsidt = xxls * qvi / dum; // budget equation: calculate DQSIDT
            abi = 1.0 + dqsidt * xxls / cpm; // budget equation: calculate ABI
            ab = 1.0 + dqsdt * xxlv / cpm; // budget equation: calculate AB

            // CASE FOR TEMPERATURE ABOVE FREEZING
            if (t3d(i,j,k) >= 273.15) {
              //......................................................................
              // ALLOW FOR CONSTANT DROPLET NUMBER
              // INUM = 0, PREDICT DROPLET NUMBER
              // INUM = 1, SET CONSTANT DROPLET NUMBER

              if (m_inum == 1) {
                // CONVERT NDCNST FROM CM-3 TO KG-1
                // Note: NDCNST constant would need to be defined elsewhere
                nc3d(i,j,k) = m_ndcnst * 1.0e6 / rho(i,j,k); // Set cloud droplet number concentration
              }

              // GET SIZE DISTRIBUTION PARAMETERS
              // MELT VERY SMALL SNOW AND GRAUPEL MIXING RATIOS, ADD TO RAIN
              if (qni3d(i,j,k) < 1.0e-6) {
                qr3d(i,j,k) = qr3d(i,j,k) + qni3d(i,j,k);         // Transfer snow to rain
                nr3d(i,j,k) = nr3d(i,j,k) + ns3d(i,j,k);          // Transfer snow number to rain
                t3d(i,j,k) = t3d(i,j,k) - qni3d(i,j,k) * xlf / cpm; // Adjust temperature
                qni3d(i,j,k) = 0.0;                 // Set snow to zero
                ns3d(i,j,k) = 0.0;                  // Set snow number to zero
              }

              if (qg3d(i,j,k) < 1.0e-6) {
                qr3d(i,j,k) = qr3d(i,j,k) + qg3d(i,j,k);          // Transfer graupel to rain
                nr3d(i,j,k) = nr3d(i,j,k) + ng3d(i,j,k);          // Transfer graupel number to rain
                t3d(i,j,k) = t3d(i,j,k) - qg3d(i,j,k) * xlf / cpm;  // Adjust temperature
                qg3d(i,j,k) = 0.0;                  // Set graupel to zero
                ng3d(i,j,k) = 0.0;                  // Set graupel number to zero
              }

              // Skip to label 300 if concentrations are below thresholds
              if (qc3d(i,j,k) < m_qsmall && qni3d(i,j,k) < 1.0e-8 && qr3d(i,j,k) < m_qsmall && qg3d(i,j,k) < 1.0e-8) {
                skipConcentrations=true; // goto 300
              }
#if 0
              // C++ version
              if ((i >= 86 && i <= 101 && j >= 0 && j <= 3 && k >= 8 && k <= 23 &&
                   (i-86)%2 == 0 && j%2 == 0 && (k-8)%2 == 0) ||
                  (i == 92 && j == 0 && k == 17) ||
                  (((i == 168 || i == 169 || i == 190 || i == 191) &&
                    (j == 0 || j == 3) &&
                    (k == 0 || k == 1 || k == 126 || k == 127))) ||
                  (i == 175 && j == 1 && k == 50) ||
                  (i == 180 && j == 2 && k == 75) ||
                  (i == 185 && j == 1 && k == 100) ||
                  (i == 170 && j == 0 && k == 30) ||
                  (i == 188 && j == 3 && k == 60) ||
                  (i == 178 && j == 2 && k == 40) ||
                  (i == 183 && j == 0 && k == 80) ||
                  (i == 173 && j == 3 && k == 110) ||
                  (i == 186 && j == 1 && k == 90) ||
                  (i == 177 && j == 2 && k == 65)) {
                fprintf(file, "%5d %5d %5d %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n",
                        i, j, k, t3d(i,j,k), nc3d(i,j,k), qr3d(i,j,k), nr3d(i,j,k), qni3d(i,j,k), ns3d(i,j,k), qg3d(i,j,k), ng3d(i,j,k));
              }
#endif
              if(!skipConcentrations) {
                ns3d(i,j,k) = amrex::max(0.0,ns3d(i,j,k));
                nc3d(i,j,k) = amrex::max(0.0,nc3d(i,j,k));
                nr3d(i,j,k) = amrex::max(0.0,nr3d(i,j,k));
                ng3d(i,j,k) = amrex::max(0.0,ng3d(i,j,k));

                // ========================================================================
                // USING WRF APPROACH FOR SIZE DISTRIBUTION PARAMETERS
                // ========================================================================
                // Rain
                if (qr3d(i,j,k) >= m_qsmall) {
                  // Calculate lambda parameter using cons26 (pi*rhow/6)
                  lamr = pow(m_cons26 * nr3d(i,j,k) / qr3d(i,j,k), 1.0/3.0);

                  // Check for slope and adjust vars
                  if (lamr < m_lamminr) {
                    lamr = m_lamminr;
                    n0r = pow(lamr, 4.0) * qr3d(i,j,k) / m_cons26;
                    nr3d(i,j,k) = n0r / lamr;  // Update number concentration
                  } else if (lamr > m_lammaxr) {
                    lamr = m_lammaxr;
                    n0r = pow(lamr, 4.0) * qr3d(i,j,k) / m_cons26;
                    nr3d(i,j,k) = n0r / lamr;  // Update number concentration
                  } else {
                    // Calculate intercept parameter using WRF formula
                    n0r = pow(lamr, 4.0) * qr3d(i,j,k) / m_cons26;
                  }
                }

                // Cloud droplets
                if (qc3d(i,j,k) >= m_qsmall) {
                  // Calculate air density factor (moist air density)
                  amrex::Real dum = pres(i,j,k)/(287.15*t3d(i,j,k));

                  // MARTIN ET AL. (1994) FORMULA FOR PGAM (WRF implementation)
                  pgam = 0.0005714*(nc3d(i,j,k)/1.0e6*dum) + 0.2714;
                  pgam = 1.0/(pgam*pgam) - 1.0;
                  pgam = amrex::max(pgam, 2.0);
                  pgam = amrex::min(pgam, 10.0);

                  // Calculate gamma function values
                  amrex::Real gamma_pgam_plus_1 = tgamma(pgam + 1.0);
                  amrex::Real gamma_pgam_plus_4 = tgamma(pgam + 4.0);

                  // Calculate lambda parameter
                  lamc = pow((m_cons26 * nc3d(i,j,k) * gamma_pgam_plus_4) / (qc3d(i,j,k) * gamma_pgam_plus_1), 1.0/3.0);

                  // Lambda bounds from WRF - 60 micron max diameter, 1 micron min diameter
                  amrex::Real lambda_min = (pgam + 1.0)/60.0e-6;
                  amrex::Real lambda_max = (pgam + 1.0)/1.0e-6;

                  // Check bounds and update number concentration if needed
                  if (lamc < lambda_min) {
                    lamc = lambda_min;
                    // Update cloud droplet number using the same formula as in WRF
                    nc3d(i,j,k) = exp(3.0*log(lamc) + log(qc3d(i,j,k)) +
                               log(gamma_pgam_plus_1) - log(gamma_pgam_plus_4))/ m_cons26;
                  } else if (lamc > lambda_max) {
                    lamc = lambda_max;
                    // Update cloud droplet number using the same formula as in WRF
                    nc3d(i,j,k) = exp(3.0*log(lamc) + log(qc3d(i,j,k)) +
                               log(gamma_pgam_plus_1) - log(gamma_pgam_plus_4))/ m_cons26;
                  }

                  // Calculate intercept parameter
                  cdist1 = nc3d(i,j,k) * pow(lamc, pgam+1) / gamma_pgam_plus_1;
                }

                // Snow
                if (qni3d(i,j,k) >= m_qsmall) {
                  // Calculate lambda parameter
                  lams = pow(m_cons1 * ns3d(i,j,k) / qni3d(i,j,k), 1.0/ds0);

                  // Calculate intercept parameter
                  n0s = ns3d(i,j,k) * lams;

                  // Check for slope and adjust vars
                  if (lams < m_lammins) {
                    lams = m_lammins;
                    n0s = pow(lams, 4.0) * qni3d(i,j,k) / m_cons1;
                    ns3d(i,j,k) = n0s / lams;  // Update number concentration
                  } else if (lams > m_lammaxs) {
                    lams = m_lammaxs;
                    n0s = pow(lams, 4.0) * qni3d(i,j,k) / m_cons1;
                    ns3d(i,j,k) = n0s / lams;  // Update number concentration
                  }
                }

                // Cloud ice
                if (qi3d(i,j,k) >= m_qsmall) {
                  // Calculate lambda parameter
                  lami = pow(m_cons12 * ni3d(i,j,k) / qi3d(i,j,k), 1.0/3.0);

                  // Calculate intercept parameter (initial calculation)
                  n0i = ni3d(i,j,k) * lami;

                  // Check for slope (apply bounds)
                  if (lami < m_lammini) {
                    lami = m_lammini;
                    // Recalculate n0i when lambda is adjusted
                    n0i = pow(lami, 4.0) * qi3d(i,j,k) / m_cons12;
                    // Update ni3d when lambda is adjusted
                    ni3d(i,j,k) = n0i / lami;
                  } else if (lami > m_lammaxi) {
                    lami = m_lammaxi;
                    // Recalculate n0i when lambda is adjusted
                    n0i = pow(lami, 4.0) * qi3d(i,j,k) / m_cons12;
                    // Update ni3d when lambda is adjusted
                    ni3d(i,j,k) = n0i / lami;
                  }
                }

                // Graupel
                if (qg3d(i,j,k) >= m_qsmall) {
                  // Calculate lambda parameter
                  lamg = pow(m_cons2 * ng3d(i,j,k) / qg3d(i,j,k), 1.0/dg0);

                  // Calculate intercept parameter
                  n0g = ng3d(i,j,k) * lamg;

                  // Check for slope and adjust vars
                  if (lamg < m_lamming) {
                    lamg = m_lamming;
                    n0g = pow(lamg, 4.0) * qg3d(i,j,k) / m_cons2;
                    ng3d(i,j,k) = n0g / lamg;  // Update number concentration
                  } else if (lamg > m_lammaxg) {
                    lamg = m_lammaxg;
                    n0g = pow(lamg, 4.0) * qg3d(i,j,k) / m_cons2;
                    ng3d(i,j,k) = n0g / lamg;  // Update number concentration
                  }
                }
#if 0
                // C++ version
                if ((i >= 86 && i <= 101 && j >= 0 && j <= 3 && k >= 8 && k <= 23 &&
                     (i-86)%2 == 0 && j%2 == 0 && (k-8)%2 == 0) ||
                    (i == 92 && j == 0 && k == 17) ||
                    (((i == 168 || i == 169 || i == 190 || i == 191) &&
                      (j == 0 || j == 3) &&
                      (k == 0 || k == 1 || k == 126 || k == 127))) ||
                    (i == 175 && j == 1 && k == 50) ||
                    (i == 180 && j == 2 && k == 75) ||
                    (i == 185 && j == 1 && k == 100) ||
                    (i == 170 && j == 0 && k == 30) ||
                    (i == 188 && j == 3 && k == 60) ||
                    (i == 178 && j == 2 && k == 40) ||
                    (i == 183 && j == 0 && k == 80) ||
                    (i == 173 && j == 3 && k == 110) ||
                    (i == 186 && j == 1 && k == 90) ||
                    (i == 177 && j == 2 && k == 65)) {
                  fprintf(file, "%5d %5d %5d %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n",
                          i, j, k, lamr, n0r, pgam, lamc, nc3d(i,j,k), lams, n0s, ns3d(i,j,k), lamg, n0g, ng3d(i,j,k));
                }
#endif
                ////////////////////// First instance of ZERO OUT PROCESS RATES
                      printf("ERROR: Concentrations not fully implmented in C++");
              }
              //Right after 300 CONTINUE
              // Calculate saturation adjustment to condense extra vapor above water saturation
              dumt = t3d(i,j,k) + dt * t3dten(i,j,k);
              dumqv = qv3d(i,j,k) + dt * qv3dten(i,j,k);

              // Fix for low pressure (added 5/12/10)
              dum = std::min(0.99 * pres(i,j,k), calc_saturation_vapor_pressure(dumt, 0));
              dumqss = m_ep_2 * dum / (pres(i,j,k) - dum);
              dumqc = qc3d(i,j,k) + dt * qc3dten(i,j,k);
              dumqc = std::max(dumqc, 0.0);

              // Saturation adjustment for liquid
              dums = dumqv - dumqss;
              pcc = dums / (1.0 + xxlv * xxlv * dumqss / (cpm * m_Rv * dumt * dumt)) / dt;
              if (pcc * dt + dumqc < 0.0) {
                pcc = -dumqc / dt;
              }

              // Update tendencies
              qv3dten(i,j,k) -= pcc;
              t3dten(i,j,k) += pcc * xxlv / cpm;
              qc3dten(i,j,k) += pcc;
#if 0
              if (i == 93 && j == 3 && k == 18) {
                fprintf(file, "%5d %5d %5d %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n",
                        i, j, k, t3d(i,j,k), qv3d(i,j,k), pres(i,j,k), qc3d(i,j,k), t3dten(i,j,k), qv3dten(i,j,k), qc3dten(i,j,k),
                        dumt, dumqv, dum, dumqss, dumqc, dums, pcc);
              }
#endif
            } else {
            // not implmented : ELSE  ! TEMPERATURE < 273.15
              printf("ERROR: Microphysics not fully implmented in C++, missing cold\n");
            }
              ltrue = 1;
              skipPrecip = false;

            }
            // 200
            }
            for(int k=klo; k<=khi; k++) {
            // INITIALIZE PRECIP AND SNOW RATES
            precrt(i,j,k) = 0.0;
            snowrt(i,j,k) = 0.0;
        // hm added 7/13/13
            snowprt(i,j,k) = 0.0;
            grplprt(i,j,k) = 0.0;
            }
            nstep = 1;
            if(ltrue != 0) {
            //goto 400
            // CALCULATE SEDIMENTATION
            // THE NUMERICS HERE FOLLOW FROM REISNER ET AL. (1998)
            // FALLOUT TERMS ARE CALCULATED ON SPLIT TIME STEPS TO ENSURE NUMERICAL
            // STABILITY, I.E. COURANT# < 1
            // Loop from top to bottom (KTE to KTS)
            for(int k=khi; k>=klo; k--) {
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

              amrex::Real dum;                // DUM: General dummy variable

              amrex::Real di0;                // DC0: Characteristic diameter for ice
              amrex::Real dc0;                // DC0: Characteristic diameter for cloud droplets
              amrex::Real ds0;                // DS0: Characteristic diameter for snow
              amrex::Real dg0;                // DG0: Characteristic diameter for graupel
              amrex::Real lammax;             // LAMMAX: Maximum value for slope parameter
              amrex::Real lammin;             // LAMMIN: Minimum value for slope parameter

              ds0 = 3.0;       // Size distribution parameter for snow
              di0 = 3.0;       // Size distribution parameter for cloud ice
              dg0 = 3.0;       // Size distribution parameter for graupel

              // Update prognostic variables with tendencies
              dumi(i,j,k) = qi3d(i,j,k) + qi3dten(i,j,k) * dt;
              dumqs(i,j,k) = qni3d(i,j,k) + qni3dten(i,j,k) * dt;
              dumr(i,j,k) = qr3d(i,j,k) + qr3dten(i,j,k) * dt;
              dumfni(i,j,k) = ni3d(i,j,k) + ni3dten(i,j,k) * dt;
              dumfns(i,j,k) = ns3d(i,j,k) + ns3dten(i,j,k) * dt;
              dumfnr(i,j,k) = nr3d(i,j,k) + nr3dten(i,j,k) * dt;
              dumc(i,j,k) = qc3d(i,j,k) + qc3dten(i,j,k) * dt;
              dumfnc(i,j,k) = nc3d(i,j,k) + nc3dten(i,j,k) * dt;
              dumg(i,j,k) = qg3d(i,j,k) + qg3dten(i,j,k) * dt;
              dumfng(i,j,k) = ng3d(i,j,k) + ng3dten(i,j,k) * dt;

              // SWITCH FOR CONSTANT DROPLET NUMBER
              if (iinum == 1) {
                dumfnc(i,j,k) = nc3d(i,j,k);
              }

              // MAKE SURE NUMBER CONCENTRATIONS ARE POSITIVE
              dumfni(i,j,k) = amrex::max(0., dumfni(i,j,k));
              dumfns(i,j,k) = amrex::max(0., dumfns(i,j,k));
              dumfnc(i,j,k) = amrex::max(0., dumfnc(i,j,k));
              dumfnr(i,j,k) = amrex::max(0., dumfnr(i,j,k));
              dumfng(i,j,k) = amrex::max(0., dumfng(i,j,k));

              // CLOUD ICE
              if (dumi(i,j,k) >= m_qsmall) {
                dlami(i,j,k) = std::pow(m_cons12 * dumfni(i,j,k) / dumi(i,j,k), 1.0/di0);
                dlami(i,j,k) = amrex::max(dlami(i,j,k), m_lammini);
                dlami(i,j,k) = amrex::min(dlami(i,j,k), m_lammaxi);
              }

              // RAIN
              if (dumr(i,j,k) >= m_qsmall) {
                dlamr(i,j,k) = std::pow(m_pi * m_rhow * dumfnr(i,j,k) / dumr(i,j,k), 1.0/3.0);
                dlamr(i,j,k) = amrex::max(dlamr(i,j,k), m_lamminr);
                dlamr(i,j,k) = amrex::min(dlamr(i,j,k), m_lammaxr);
              }

              // CLOUD DROPLETS
              if (dumc(i,j,k) >= m_qsmall) {
                amrex::Real dum = pres(i,j,k) / (287.15 * t3d(i,j,k));
                amrex::Real pgam = 0.0005714 * (nc3d(i,j,k) / 1.0e6 * dum) + 0.2714;
                pgam = 1.0 / (pgam * pgam) - 1.0;
                pgam = amrex::max(pgam, 2.0);
                pgam = amrex::min(pgam, 10.0);

                dlamc(i,j,k) = std::pow(m_cons26 * dumfnc(i,j,k) * std::tgamma(pgam + 4.0) /
                                        (dumc(i,j,k) * std::tgamma(pgam + 1.0)), 1.0/3.0);
                lammin = (pgam + 1.0) / 60.0e-6;
                lammax = (pgam + 1.0) / 1.0e-6;
                dlamc(i,j,k) = amrex::max(dlamc(i,j,k), lammin);
                dlamc(i,j,k) = amrex::min(dlamc(i,j,k), lammax);
              }

              // SNOW
              if (dumqs(i,j,k) >= m_qsmall) {
                dlams(i,j,k) = std::pow(m_cons1 * dumfns(i,j,k) / dumqs(i,j,k), 1.0/ds0);
                dlams(i,j,k) = amrex::max(dlams(i,j,k), m_lammins);
                dlams(i,j,k) = amrex::min(dlams(i,j,k), m_lammaxs);
              }

              // GRAUPEL
              if (dumg(i,j,k) >= m_qsmall) {
                dlamg(i,j,k) = std::pow(m_cons2 * dumfng(i,j,k) / dumg(i,j,k), 1.0/dg0);
                dlamg(i,j,k) = amrex::max(dlamg(i,j,k), m_lamming);
                dlamg(i,j,k) = amrex::min(dlamg(i,j,k), m_lammaxg);
              }

              // Calculate number-weighted and mass-weighted terminal fall speeds

              // CLOUD WATER
              if (dumc(i,j,k) >= m_qsmall) {
                unc(i,j,k) = acn(i,j,k) * gamma_function(1. + m_bc + pgam) / (dlamc(i,j,k) * std::pow(dlamc(i,j,k), m_bc) * gamma_function(pgam + 1.));
                umc(i,j,k) = acn(i,j,k) * gamma_function(4. + m_bc + pgam) / (dlamc(i,j,k) * std::pow(dlamc(i,j,k), m_bc) * gamma_function(pgam + 4.));
              } else {
                umc(i,j,k) = 0.;
                unc(i,j,k) = 0.;
              }

              // CLOUD ICE
              if (dumi(i,j,k) >= m_qsmall) {
                uni(i,j,k) = ain(i,j,k) * m_cons27 / std::pow(dlami(i,j,k), m_bi);
                umi(i,j,k) = ain(i,j,k) * m_cons28 / std::pow(dlami(i,j,k), m_bi);
              } else {
                umi(i,j,k) = 0.;
                uni(i,j,k) = 0.;
              }

              // RAIN
              if (dumr(i,j,k) >= m_qsmall) {
                unr(i,j,k) = arn(i,j,k) * m_cons6 / std::pow(dlamr(i,j,k), m_br);
                umr(i,j,k) = arn(i,j,k) * m_cons4 / std::pow(dlamr(i,j,k), m_br);
              } else {
                umr(i,j,k) = 0.;
                unr(i,j,k) = 0.;
              }

              // SNOW
              if (dumqs(i,j,k) >= m_qsmall) {
                ums(i,j,k) = asn(i,j,k) * m_cons3 / std::pow(dlams(i,j,k), m_bs);
                uns(i,j,k) = asn(i,j,k) * m_cons5 / std::pow(dlams(i,j,k), m_bs);
              } else {
                ums(i,j,k) = 0.;
                uns(i,j,k) = 0.;
              }

              // GRAUPEL
              if (dumg(i,j,k) >= m_qsmall) {
                umg(i,j,k) = agn(i,j,k) * m_cons7 / std::pow(dlamg(i,j,k), m_bg);
                ung(i,j,k) = agn(i,j,k) * m_cons8 / std::pow(dlamg(i,j,k), m_bg);
              } else {
                umg(i,j,k) = 0.;
                ung(i,j,k) = 0.;
              }

              // SET REALISTIC LIMITS ON FALLSPEED
              // Bug fix, 10/08/09
              dum = std::pow(m_rhosu / rho(i,j,k), 0.54);
              ums(i,j,k) = std::min(ums(i,j,k), 1.2 * dum);
              uns(i,j,k) = std::min(uns(i,j,k), 1.2 * dum);

              // Fix 053011
              // Fix for correction by AA 4/6/11
              umi(i,j,k) = std::min(umi(i,j,k), 1.2 * std::pow(m_rhosu / rho(i,j,k), 0.35));
              uni(i,j,k) = std::min(uni(i,j,k), 1.2 * std::pow(m_rhosu / rho(i,j,k), 0.35));
              umr(i,j,k) = std::min(umr(i,j,k), 9.1 * dum);
              unr(i,j,k) = std::min(unr(i,j,k), 9.1 * dum);
              umg(i,j,k) = std::min(umg(i,j,k), 20. * dum);
              ung(i,j,k) = std::min(ung(i,j,k), 20. * dum);

              // Set fall speed values
              fr(i,j,k) = umr(i,j,k);         // RAIN FALL SPEED
              fi(i,j,k) = umi(i,j,k);         // CLOUD ICE FALL SPEED
              fni(i,j,k) = uni(i,j,k);        // CLOUD ICE NUMBER FALL SPEED
              fs(i,j,k) = ums(i,j,k);         // SNOW FALL SPEED
              fns(i,j,k) = uns(i,j,k);        // SNOW NUMBER FALL SPEED
              fnr(i,j,k) = unr(i,j,k);        // RAIN NUMBER FALL SPEED
              fc(i,j,k) = umc(i,j,k);         // CLOUD WATER FALL SPEED
              fnc(i,j,k) = unc(i,j,k);        // CLOUD NUMBER FALL SPEED
              fg(i,j,k) = umg(i,j,k);         // GRAUPEL FALL SPEED
              fng(i,j,k) = ung(i,j,k);        // GRAUPEL NUMBER FALL SPEED

              // V3.3 MODIFY FALLSPEED BELOW LEVEL OF PRECIP
              if (fr(i,j,k) < 1.e-10) {
                fr(i,j,k) = fr(i,j,k+1);
              }
              if (fi(i,j,k) < 1.e-10) {
                fi(i,j,k) = fi(i,j,k+1);
              }
              if (fni(i,j,k) < 1.e-10) {
                fni(i,j,k) = fni(i,j,k+1);
              }
              if (fs(i,j,k) < 1.e-10) {
                fs(i,j,k) = fs(i,j,k+1);
              }
              if (fns(i,j,k) < 1.e-10) {
                fns(i,j,k) = fns(i,j,k+1);
              }
              if (fnr(i,j,k) < 1.e-10) {
                fnr(i,j,k) = fnr(i,j,k+1);
              }
              if (fc(i,j,k) < 1.e-10) {
                fc(i,j,k) = fc(i,j,k+1);
              }
              if (fnc(i,j,k) < 1.e-10) {
                fnc(i,j,k) = fnc(i,j,k+1);
              }
              if (fg(i,j,k) < 1.e-10) {
                fg(i,j,k) = fg(i,j,k+1);
              }
              if (fng(i,j,k) < 1.e-10) {
                fng(i,j,k) = fng(i,j,k+1);
              }

              // CALCULATE NUMBER OF SPLIT TIME STEPS
              // Find maximum fall speed at this point
              rgvm(i,j,k) = std::max({fr(i,j,k), fi(i,j,k), fs(i,j,k), fc(i,j,k),
                  fni(i,j,k), fnr(i,j,k), fns(i,j,k), fnc(i,j,k),
                  fg(i,j,k), fng(i,j,k)});

              // Calculate number of steps (dt and nstep would need to be defined elsewhere)
              int nstep = std::max(static_cast<int>(rgvm(i,j,k) * dt / dzq(i,j,k) + 1.), nstep);

              // MULTIPLY VARIABLES BY RHO
              dumr(i,j,k) = dumr(i,j,k) * rho(i,j,k);       // Rain water content * density
              dumi(i,j,k) = dumi(i,j,k) * rho(i,j,k);       // Cloud ice content * density
              dumfni(i,j,k) = dumfni(i,j,k) * rho(i,j,k);   // Cloud ice number * density
              dumqs(i,j,k) = dumqs(i,j,k) * rho(i,j,k);     // Snow content * density
              dumfns(i,j,k) = dumfns(i,j,k) * rho(i,j,k);   // Snow number * density
              dumfnr(i,j,k) = dumfnr(i,j,k) * rho(i,j,k);   // Rain number * density
              dumc(i,j,k) = dumc(i,j,k) * rho(i,j,k);       // Cloud water content * density
              dumfnc(i,j,k) = dumfnc(i,j,k) * rho(i,j,k);   // Cloud droplet number * density
              dumg(i,j,k) = dumg(i,j,k) * rho(i,j,k);       // Graupel content * density
              dumfng(i,j,k) = dumfng(i,j,k) * rho(i,j,k);   // Graupel number * density
            }
            // Main time stepping loop for sedimentation
            for (int n = 0; n <= nstep; n++) {
              // Calculate initial fallout for each hydrometeor type for all levels
              for (int k = klo; k <= khi; k++) {
                faloutr(i,j,k) = fr(i,j,k) * dumr(i,j,k);
                falouti(i,j,k) = fi(i,j,k) * dumi(i,j,k);
                faloutni(i,j,k) = fni(i,j,k) * dumfni(i,j,k);
                falouts(i,j,k) = fs(i,j,k) * dumqs(i,j,k);
                faloutns(i,j,k) = fns(i,j,k) * dumfns(i,j,k);
                faloutnr(i,j,k) = fnr(i,j,k) * dumfnr(i,j,k);
                faloutc(i,j,k) = fc(i,j,k) * dumc(i,j,k);
                faloutnc(i,j,k) = fnc(i,j,k) * dumfnc(i,j,k);
                faloutg(i,j,k) = fg(i,j,k) * dumg(i,j,k);
                faloutng(i,j,k) = fng(i,j,k) * dumfng(i,j,k);
              }

              // Process top of model level
              int k = khi;

              // Calculate tendencies at top level
              faltndr(i,j,k) = faloutr(i,j,k) / dzq(i,j,k);
              faltndi(i,j,k) = falouti(i,j,k) / dzq(i,j,k);
              faltndni(i,j,k) = faloutni(i,j,k) / dzq(i,j,k);
              faltnds(i,j,k) = falouts(i,j,k) / dzq(i,j,k);
              faltndns(i,j,k) = faloutns(i,j,k) / dzq(i,j,k);
              faltndnr(i,j,k) = faloutnr(i,j,k) / dzq(i,j,k);
              faltndc(i,j,k) = faloutc(i,j,k) / dzq(i,j,k);
              faltndnc(i,j,k) = faloutnc(i,j,k) / dzq(i,j,k);
              faltndg(i,j,k) = faloutg(i,j,k) / dzq(i,j,k);
              faltndng(i,j,k) = faloutng(i,j,k) / dzq(i,j,k);

              // Add fallout terms to Eulerian tendencies (scaled by time step and density)
              qrsten(i,j,k) = qrsten(i,j,k) - faltndr(i,j,k) / nstep / rho(i,j,k);
              qisten(i,j,k) = qisten(i,j,k) - faltndi(i,j,k) / nstep / rho(i,j,k);
              ni3dten(i,j,k) = ni3dten(i,j,k) - faltndni(i,j,k) / nstep / rho(i,j,k);
              qnisten(i,j,k) = qnisten(i,j,k) - faltnds(i,j,k) / nstep / rho(i,j,k);
              ns3dten(i,j,k) = ns3dten(i,j,k) - faltndns(i,j,k) / nstep / rho(i,j,k);
              nr3dten(i,j,k) = nr3dten(i,j,k) - faltndnr(i,j,k) / nstep / rho(i,j,k);
              qcsten(i,j,k) = qcsten(i,j,k) - faltndc(i,j,k) / nstep / rho(i,j,k);
              nc3dten(i,j,k) = nc3dten(i,j,k) - faltndnc(i,j,k) / nstep / rho(i,j,k);
              qgsten(i,j,k) = qgsten(i,j,k) - faltndg(i,j,k) / nstep / rho(i,j,k);
              ng3dten(i,j,k) = ng3dten(i,j,k) - faltndng(i,j,k) / nstep / rho(i,j,k);

              // Update temporary working variables
              dumr(i,j,k) = dumr(i,j,k) - faltndr(i,j,k) * dt / nstep;
              dumi(i,j,k) = dumi(i,j,k) - faltndi(i,j,k) * dt / nstep;
              dumfni(i,j,k) = dumfni(i,j,k) - faltndni(i,j,k) * dt / nstep;
              dumqs(i,j,k) = dumqs(i,j,k) - faltnds(i,j,k) * dt / nstep;
              dumfns(i,j,k) = dumfns(i,j,k) - faltndns(i,j,k) * dt / nstep;
              dumfnr(i,j,k) = dumfnr(i,j,k) - faltndnr(i,j,k) * dt / nstep;
              dumc(i,j,k) = dumc(i,j,k) - faltndc(i,j,k) * dt / nstep;
              dumfnc(i,j,k) = dumfnc(i,j,k) - faltndnc(i,j,k) * dt / nstep;
              dumg(i,j,k) = dumg(i,j,k) - faltndg(i,j,k) * dt / nstep;
              dumfng(i,j,k) = dumfng(i,j,k) - faltndng(i,j,k) * dt / nstep;

              // Process remaining levels from top to bottom
              for (int k = khi-1; k >= klo; k--) {
                // Calculate tendencies based on difference between levels
                faltndr(i,j,k) = (faloutr(i,j,k+1) - faloutr(i,j,k)) / dzq(i,j,k);
                faltndi(i,j,k) = (falouti(i,j,k+1) - falouti(i,j,k)) / dzq(i,j,k);
                faltndni(i,j,k) = (faloutni(i,j,k+1) - faloutni(i,j,k)) / dzq(i,j,k);
                faltnds(i,j,k) = (falouts(i,j,k+1) - falouts(i,j,k)) / dzq(i,j,k);
                faltndns(i,j,k) = (faloutns(i,j,k+1) - faloutns(i,j,k)) / dzq(i,j,k);
                faltndnr(i,j,k) = (faloutnr(i,j,k+1) - faloutnr(i,j,k)) / dzq(i,j,k);
                faltndc(i,j,k) = (faloutc(i,j,k+1) - faloutc(i,j,k)) / dzq(i,j,k);
                faltndnc(i,j,k) = (faloutnc(i,j,k+1) - faloutnc(i,j,k)) / dzq(i,j,k);
                faltndg(i,j,k) = (faloutg(i,j,k+1) - faloutg(i,j,k)) / dzq(i,j,k);
                faltndng(i,j,k) = (faloutng(i,j,k+1) - faloutng(i,j,k)) / dzq(i,j,k);

                // Add fallout terms to Eulerian tendencies (positive here, as mass flows in from above)
                qrsten(i,j,k) = qrsten(i,j,k) + faltndr(i,j,k) / nstep / rho(i,j,k);
                qisten(i,j,k) = qisten(i,j,k) + faltndi(i,j,k) / nstep / rho(i,j,k);
                ni3dten(i,j,k) = ni3dten(i,j,k) + faltndni(i,j,k) / nstep / rho(i,j,k);
                qnisten(i,j,k) = qnisten(i,j,k) + faltnds(i,j,k) / nstep / rho(i,j,k);
                ns3dten(i,j,k) = ns3dten(i,j,k) + faltndns(i,j,k) / nstep / rho(i,j,k);
                nr3dten(i,j,k) = nr3dten(i,j,k) + faltndnr(i,j,k) / nstep / rho(i,j,k);
                qcsten(i,j,k) = qcsten(i,j,k) + faltndc(i,j,k) / nstep / rho(i,j,k);
                nc3dten(i,j,k) = nc3dten(i,j,k) + faltndnc(i,j,k) / nstep / rho(i,j,k);
                qgsten(i,j,k) = qgsten(i,j,k) + faltndg(i,j,k) / nstep / rho(i,j,k);
                ng3dten(i,j,k) = ng3dten(i,j,k) + faltndng(i,j,k) / nstep / rho(i,j,k);

                // Update temporary working variables
                dumr(i,j,k) = dumr(i,j,k) + faltndr(i,j,k) * dt / nstep;
                dumi(i,j,k) = dumi(i,j,k) + faltndi(i,j,k) * dt / nstep;
                dumfni(i,j,k) = dumfni(i,j,k) + faltndni(i,j,k) * dt / nstep;
                dumqs(i,j,k) = dumqs(i,j,k) + faltnds(i,j,k) * dt / nstep;
                dumfns(i,j,k) = dumfns(i,j,k) + faltndns(i,j,k) * dt / nstep;
                dumfnr(i,j,k) = dumfnr(i,j,k) + faltndnr(i,j,k) * dt / nstep;
                dumc(i,j,k) = dumc(i,j,k) + faltndc(i,j,k) * dt / nstep;
                dumfnc(i,j,k) = dumfnc(i,j,k) + faltndnc(i,j,k) * dt / nstep;
                dumg(i,j,k) = dumg(i,j,k) + faltndg(i,j,k) * dt / nstep;
                dumfng(i,j,k) = dumfng(i,j,k) + faltndng(i,j,k) * dt / nstep;
              }
            }
            printf("ERROR: Sedimentation not implmented in C++\n");
            }

            for(int k=klo; k<=khi; k++) {
            //End of _micro
            if(use_morr_cpp_answer) {
            // Transfer 1D variables back to 3D arrays
            qcl_arr(i,j,k) = qc3d(i,j,k);
            qci_arr(i,j,k) = qi3d(i,j,k);
            qps_arr(i,j,k) = qni3d(i,j,k);
            qpr_arr(i,j,k) = qr3d(i,j,k);
            ni_arr(i,j,k) = ni3d(i,j,k);
            ns_arr(i,j,k) = ns3d(i,j,k);
            nr_arr(i,j,k) = nr3d(i,j,k);
            qpg_arr(i,j,k) = qg3d(i,j,k);
            ng_arr(i,j,k) = ng3d(i,j,k);

            // Temperature and potential temperature conversion
            theta_arr(i,j,k) = t3d(i,j,k) / pii_arr(i,j,k); // Convert temp back to potential temp
            qv_arr(i,j,k) = qv3d(i,j,k);

            //Deleted wrf-check, effc, and precr type data as not used by ERF

            /* // NEED gpu-compatabile summation for rain_accum, check SAM or Kessler for better example
            rain_accum_arr(i,j,k) = rain_accum_arr(i,j,k) + precrt(i,j,k);
            snow_accum_arr(i,j,k) = snow_accum_arr(i,j,k) + snowprt(i,j,k);
            graup_accum_arr(i,j,k) = graup_accum_arr(i,j,k) + grplprt(i,j,k);
            // Update precipitation accumulation variables
            // These are outside the k-loop in the original code
            rain_accum_arr(i,j,klo) = rain_accum_arr(i,j,klo) + precrt(i,j,k);
            rainncv_arr(i,j) = precrt(i,j,k);
            snow_accum_arr(i,j,klo) = snow_accum_arr(i,j,klo) + snowprt(i,j,k);
            snowncv_arr(i,j) = snowprt(i,j,k);
            graup_accum_arr(i,j,klo) = graup_accum_arr(i,j,klo) + grplprt(i,j,k);
            graupelncv_arr(i,j) = grplprt(i,j,k);
            sr_arr(i,j) = snowrt(i,j,k) / (precrt(i,j,k) + 1.e-12);*/
            }
            }
         });
          fclose(file);
          //          amrex::Print()<<amrex::FArrayBox(qv_arr)<<std::endl;
          }
          amrex::Print()<<"fortran should run "<<run_morr_fort<<std::endl;
          if(run_morr_fort) {
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
        }
          //          amrex::Print()<<amrex::FArrayBox(qv_arr)<<std::endl;
          // After the call, all fields are updated
          // We don't need to copy results back since we passed direct pointers
          // to our class member arrays
        }
        }
    }

//
// This module contains constexprants that are specific to the radiative transfer
// code used in the RRTMG model.
//
#ifndef ERF_RAD_CONSTANTS_H
#define ERF_RAD_CONSTANTS_H

#include "rrtmgp_const.h"
using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

class RadConstants {
  public:
   // number of shorwave spectral intervals
   static constexpr int nswbands = 14;
   static constexpr int nbndsw = 14;

   // Wavenumbers of band boundaries
   //
   //Note: Currently rad_solar_var extends the lowest band down to
   // 100 cm^-1 if it is too high to cover the far-IR. Any changes meant
   // to affect IR solar variability should take note of this.

   static constexpr real wavenum_low[] =  // in cm^-1
     {2600., 3250., 4000., 4650., 5150., 6150., 7700.,
      8050., 12850., 16000., 22650., 29000., 38000., 820.};

   static constexpr real wavenum_high[] = // in cm^-1
      {3250., 4000., 4650., 5150., 6150., 7700., 8050.,
       12850.,16000., 22650., 29000., 38000., 50000., 2600.};

   // Solar irradiance at 1 A.U. in W/m^2 assumed by radiation code
   // Rescaled so that sum is precisely 1368.22 and fractional amounts sum to 1.0
   static constexpr real solar_ref_band_irradiance[] =
      {12.11,  20.3600000000001, 23.73,
       22.43,  55.63, 102.93, 24.29,
       345.74, 218.19, 347.20,
       129.49,  50.15,   3.08, 12.89};

   // parameterizations related to radiation (when microphysics model doesn't provide)
   static constexpr real icesize_table_min_temp = 180.;
   static constexpr real retab[] = { 5.92779, 6.26422, 6.61973, 6.99539, 7.39234,
                                     7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,
                                     10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319,
                                     15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,
                                     20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,
                                     27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,
                                     31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601,
                                     34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,
                                     38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,
                                     42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,
                                     50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,
                                     65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,
                                     93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,
                                     124.954, 130.630, 136.457, 142.446, 148.608, 154.956,
                                     161.503, 168.262, 175.248, 182.473, 189.952, 197.699,
                                     205.728, 214.055, 222.694, 231.661, 240.971, 250.639};

   // These are indices to the band for diagnostic output
   static constexpr int idx_sw_diag = 10; // index to sw visible band
   static constexpr int idx_nir_diag = 8; // index to sw near infrared (778-1240 nm) band
   static constexpr int idx_uv_diag = 11; // index to sw uv (345-441 nm) band

   static constexpr int rrtmg_sw_cloudsim_band = 9;  // rrtmg band for .67 micron

   // Number of evenly spaced intervals in rh
   // The globality of this mesh may not be necessary
   // Perhaps it could be specific to the aerosol
   // But it is difficult to see how refined it must be
   // for lookup.  This value was found to be sufficient
   // for Sulfate and probably necessary to resolve the
   // high variation near rh = 1.  Alternative methods
   // were found to be too slow.
   // Optimal approach would be for cam to specify size of aerosol
   // based on each aerosol's characteristics.  Radiation
   // should know nothing about hygroscopic growth!
   static constexpr int nrh = 1000;

   // LONGWAVE DATA

   // These are indices to the band for diagnostic output
   static constexpr int idx_lw_diag = 7; // index to (H20 window) LW band
   static constexpr int rrtmg_lw_cloudsim_band = 6;  // rrtmg band for 10.5 micron

   // number of lw bands
   static constexpr int nlwbands = 16;
   static constexpr int nbndlw = 16;

   static constexpr real wavenumber1_longwave[] =  // Longwave spectral band limits (cm-1)
       {10.,  350., 500.,   630.,  700.,  820.,  980., 1080.,
       1180., 1390., 1480., 1800., 2080., 2250., 2390., 2600.};

   static constexpr real wavenumber2_longwave[] = // Longwave spectral band limits (cm-1)
       {350.,  500.,  630.,  700.,  820.,  980., 1080., 1180.,
       1390., 1480., 1800., 2080., 2250., 2390., 2600., 3250.};

   //GASES TREATED BY RADIATION (line spectrae)

   // gasses required by radiation
   static constexpr int nradgas = 8;
   static constexpr const char* gaslist[] =
             {"H2O","O3", "O2", "CO2", "N2O", "CH4", "CFC11", "CFC12"};

   // use enum unit
   enum Units {
      inv_cm,
      meter,
      nanometer,
      micrometer,
      centimeter
   };

   // provide Solar Irradiance for each band in RRTMG
   inline
   static void get_solar_band_fraction_irrad(real1d& fractional_irradiance) {
     real tsi = 0;  // total solar irradiance
     for(auto i = 0; i < nswbands; ++i)
       tsi = tsi + solar_ref_band_irradiance[i];

     for(auto i = 0; i < nswbands; ++i)
       fractional_irradiance(i) = solar_ref_band_irradiance[i] / tsi;
   }

   // provide Total Solar Irradiance assumed by RRTMG
   inline
   static void  get_ref_total_solar_irrad(real& tsi) {
      tsi = 0;
      for(auto i = 0; i < nswbands; ++i)
        tsi = tsi + solar_ref_band_irradiance[i];
   }

   // solar irradiance in each band (W/m^2)
   inline
   static void get_ref_solar_band_irrad(real1d& band_irrad ) {
     for(auto i = 0; i < nswbands; ++i)
       band_irrad(i) = solar_ref_band_irradiance[i];
   }

   // number of solar (shortwave) bands in the rrtmg code
   inline
   static void get_number_sw_bands(int number_of_bands) {
      number_of_bands = nswbands;
   }

   // provide spectral boundaries of each longwave band
   inline
   static void get_lw_spectral_boundaries(real1d& low_boundaries,
                                          real1d& high_boundaries,
                                          Units units) {
      switch (units) {
      case inv_cm:
         for(auto i = 0; i < nlwbands; ++i) {
           low_boundaries(i)  = wavenumber1_longwave[i];
           high_boundaries(i) = wavenumber2_longwave[i];
         }
      case meter:
         for(auto i = 0; i < nlwbands; ++i) {
           low_boundaries(i)  = 1.e-2/wavenumber2_longwave[i];
           high_boundaries(i) = 1.e-2/wavenumber1_longwave[i];
         }
      case nanometer:
         for(auto i = 0; i < nlwbands; ++i) {
           low_boundaries(i)  = 1.e7/wavenumber2_longwave[i];
           high_boundaries(i) = 1.e7/wavenumber1_longwave[i];
         }
      case micrometer:
         for(auto i = 0; i < nlwbands; ++i) {
            low_boundaries(i)  = 1.e4/wavenumber2_longwave[i];
            high_boundaries(i) = 1.e4/wavenumber1_longwave[i];
         }
      case centimeter:
         for(auto i = 0; i < nlwbands; ++i) {
            low_boundaries(i)  = 1./wavenumber2_longwave[i];
            high_boundaries(i) = 1./wavenumber1_longwave[i];
         }
      default:
         amrex::Print() << "get_lw_spectral_boundaries: spectral units not acceptable\n";
      }
   }

   inline
   static void get_lw_spectral_midpoints(real1d& band_midpoints, Units units) {
      real1d lower_boundaries("lower_boundaries", nlwbands);
      real1d upper_boundaries("upper_boundaries", nlwbands);
      int iband;

      // Get band boundaries
      get_lw_spectral_boundaries(lower_boundaries, upper_boundaries, units);

      // Get band midpoints
      for(auto iband = 0; iband < nlwbands; ++iband) {
         band_midpoints(iband) = 0.5*(lower_boundaries(iband) + upper_boundaries(iband));
      }

   }

   // provide spectral boundaries of each shortwave band
   inline
   static void get_sw_spectral_boundaries(real1d& low_boundaries, real1d& high_boundaries, Units units) {
     switch (units) {
       case inv_cm:
         for(auto i = 0; i < nswbands; ++i) {
            low_boundaries(i) = wavenum_low[i];
            high_boundaries(i) = wavenum_high[i];
         }
       case meter:
         for(auto i = 0; i < nswbands; ++i) {
           low_boundaries(i) = 1.e-2/wavenum_high[i];
           high_boundaries(i) = 1.e-2/wavenum_low[i];
         }
       case nanometer:
         for(auto i = 0; i < nswbands; ++i) {
            low_boundaries(i) = 1.e7/wavenum_high[i];
            high_boundaries(i) = 1.e7/wavenum_low[i];
         }
      case micrometer:
         for(auto i = 0; i < nswbands; ++i) {
            low_boundaries(i) = 1.e4/wavenum_high[i];
            high_boundaries(i) = 1.e4/wavenum_low[i];
         }
      case centimeter:
         for(auto i = 0; i < nswbands; ++i) {
            low_boundaries(i)  = 1./wavenum_high[i];
            high_boundaries(i) = 1./wavenum_low[i];
         }
      default:
        amrex::Print() << "rad_constants.F90: spectral units not acceptable\n";
    }
  }

  inline
  static void get_sw_spectral_midpoints(real1d& band_midpoints, Units units) {
     real1d lower_boundaries("lower_boundaries", nswbands);
     real1d upper_boundaries("upper_boundaries", nswbands);

     // Get band boundaries
     get_sw_spectral_boundaries(lower_boundaries, upper_boundaries, units);

     //Get band midpoints
     for(auto iband = 0; iband < nswbands; ++iband) {
        band_midpoints(iband) = 0.5 * (lower_boundaries(iband) + upper_boundaries(iband));
     }
  }

  // return the index in the gaslist array of the specified gasname
  inline
  static int rad_gas_index(std::string gasname) {
     for(auto igas = 0; igas < nradgas; ++igas) {
        if (gaslist[igas] == gasname) return igas;
     }
     amrex::Print() << "rad_gas_index: can not find gas with name" << gasname << std::endl;
  }
};
#endif

//
// This module provides an interface to water vapor saturation methods, providing
// saturation vapor pressure and related calculations to CAM.
//
// The original wv_saturation codes were introduced by J. J. Hack,
// February 1990. The code has been extensively rewritten since then,
// including a total refactoring in Summer 2012.
//
// Methods:
//
// Pure water/ice saturation vapor pressures are calculated on the
// fly, with the specific method determined by a runtime option.
//
// The default method for calculating SVP is determined by a namelist
// option, and used whenever svp_water/ice or qsat are called.
//
#ifndef ERF_WATER_VAPOR_SATURATION_H_
#define ERF_WATER_VAPOR_SATURATION_H_

#include <string>
#include <vector>
#include <memory>

#include "ERF_Constants.H"
#include "Sat_methods.H"

// Radiation code interface class
class WaterVaporSat {
 public:
  // Make these public parameters in case another module wants to see the
  // extent of the table.
  static constexpr amrex::Real tmin = 127.16;
  static constexpr amrex::Real tmax = 375.16;

  // Set coefficients for polynomial approximation of difference
  // between saturation vapor press over water and saturation pressure
  // over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial is
  // valid in the range -40 < t < 0 (degrees C).
  static constexpr int npcf = 5;
  static constexpr const amrex::Real pcf[npcf] = {5.04469588506e-01,
                                                 -5.47288442819e+00,
                                                 -3.67471858735e-01,
                                                 -8.95963532403e-03,
                                                 -7.78053686625e-05};

  // Compute saturation vapor pressure over water
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static amrex::Real svp_water(const amrex::Real& t) {
    return SatMethods::wv_sat_svp_water(t);
  }

  // Compute saturation vapor pressure over ice
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static amrex::Real svp_ice(const amrex::Real& t) {
    return SatMethods::wv_sat_svp_ice(t);
  }

  // Compute saturation vapor pressure with an ice-water transition
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static amrex::Real svp_trans(const amrex::Real& t) {
    return SatMethods::wv_sat_svp_trans(t);
  }

  // Get enthalpy based only on temperature
  // and specific humidity.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static amrex::Real tq_enthalpy(const amrex::Real& t,
                                 const amrex::Real& q,
                                 const amrex::Real& hltalt) {
    return Cp_d * t + hltalt * q;
  }

  //------------------------------------------------
  // LATENT HEAT OF VAPORIZATION CORRECTIONS
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void no_ip_hltalt(const amrex::Real& t,
                           amrex::Real hltalt) {
    hltalt = latvap;
    // Account for change of latvap with t above freezing where
    // constant slope is given by -2369 j/(kg c) = cpv - cw
    if (t >= tmelt) hltalt = hltalt - 2369.0*(t-tmelt);
  }

  //  Calculate latent heat of vaporization of water at a given
  //  temperature, taking into account the ice phase if temperature
  //  is below freezing.
  //  Optional argument also calculates a term used to calculate
  //  d(es)/dT within the water-ice transition range.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void calc_hltalt(const amrex::Real& t,
                          amrex::Real hltalt,
                          amrex::Real tterm = 0.) {
    amrex::Real tc, weight;
    no_ip_hltalt(t, hltalt);
    if (t < tmelt) {
       // Weighting of hlat accounts for transition from water to ice.
       tc = t - tmelt;
       if (tc >= -ttrice) {
          weight = -tc/ttrice;

          // polynomial expression approximates difference between es
          // over water and es over ice from 0 to -ttrice (C) (max of
          // ttrice is 40): required for accurate estimate of es
          // derivative in transition range from ice to water

          for(auto i = npcf-1; i > 0;  --i) tterm = pcf[i] + tc*tterm;
          tterm = tterm/ttrice;
       }
       else {
          weight = 1.0;
       }

       hltalt = hltalt + weight*latice;
    }
  }

  // Temperature derivative outputs, for qsat_*
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void deriv_outputs(const amrex::Real& t,
                            const amrex::Real& p,
                            const amrex::Real& es,
                            const amrex::Real& qs,
                            const amrex::Real& hltalt,
                            const amrex::Real& tterm,
                            amrex::Real& gam,
                            amrex::Real& dqsdt) {

    // Local variables
    amrex::Real desdt;        // d(es)/dt
    amrex::Real dqsdt_loc;    // local copy of dqsdt

    if (qs == 1.0) {
       dqsdt_loc = 0.;
    }
    else {
       desdt = hltalt*es/(R_v*t*t) + tterm;
       dqsdt_loc = qs*p*desdt/(es*(p-omeps*es));
    }

   dqsdt = dqsdt_loc;
   gam   = dqsdt_loc * (hltalt/Cp_d);
  }

  // Look up and return saturation vapor pressure from precomputed
  // table, then calculate and return saturation specific humidity.
  // Optionally return various temperature derivatives or enthalpy
  // at saturation.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void qsat(const amrex::Real& t,
                   const amrex::Real& p,
                   amrex::Real& es,
                   amrex::Real& qs,
                   amrex::Real gam = 0.,
                   amrex::Real dqsdt = 0. ,
                   amrex::Real enthalpy = 0.) {
    // Local variables
    amrex::Real hltalt;       // Modified latent heat for T derivatives
    amrex::Real tterm;        // Account for d(es)/dT in transition region

    es = svp_trans(t); //estblf(t);

    qs = SatMethods::wv_sat_svp_to_qsat(es, p);

    // Ensures returned es is consistent with limiters on qs.
    es = std::min(es, p);

    // Calculate optional arguments.
    // "generalized" analytic expression for t derivative of es
    // accurate to within 1 percent for 173.16 < t < 373.16
    calc_hltalt(t, hltalt, tterm);

    enthalpy = tq_enthalpy(t, qs, hltalt);

    deriv_outputs(t, p, es, qs, hltalt, tterm, gam, dqsdt);
  }

  // Calculate SVP over water at a given temperature, and then
  // calculate and return saturation specific humidity.
  // Optionally return various temperature derivatives or enthalpy
  // at saturation.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void qsat_water(const amrex::Real& t,
                         const amrex::Real& p,
                         amrex::Real& es,
                         amrex::Real& qs,
                         amrex::Real gam = 0.,
                         amrex::Real dqsdt = 0.,
                         amrex::Real enthalpy = 0.) {
    // Local variables
    amrex::Real hltalt;       // Modified latent heat for T derivatives

    SatMethods::wv_sat_qsat_water(t, p, es, qs);

    // "generalized" analytic expression for t derivative of es
    // accurate to within 1 percent for 173.16 < t < 373.16

    no_ip_hltalt(t, hltalt);

    enthalpy = tq_enthalpy(t, qs, hltalt);

    // For pure water/ice transition term is 0.
    deriv_outputs(t, p, es, qs, hltalt, 0., gam, dqsdt);
  }

  //  Calculate SVP over ice at a given temperature, and then
  //  calculate and return saturation specific humidity.
  //  Optionally return various temperature derivatives or enthalpy
  //  at saturation.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void qsat_ice(const amrex::Real& t,
                       const amrex::Real& p,
                       amrex::Real& es,
                       amrex::Real& qs,
                       amrex::Real& gam,
                       amrex::Real& dqsdt,
                       amrex::Real& enthalpy) {
    // Local variables
    amrex::Real hltalt;       // Modified latent heat for T derivatives

    SatMethods::wv_sat_qsat_ice(t, p, es, qs);

    // For pure ice, just add latent heats.
    hltalt = latvap + latice;

    enthalpy = tq_enthalpy(t, qs, hltalt);

    // For pure water/ice transition term is 0.
    deriv_outputs(t, p, es, qs, hltalt, 0., gam, dqsdt);
  }

  //   find the wet bulb temperature for a given t and q
  //   in a longitude height section
  //   wet bulb temp is the temperature and spec humidity that is
  //   just saturated and has the same enthalpy
  //   if q > qs(t) then tsp > t and qsp = qs(tsp) < q
  //   if q < qs(t) then tsp < t and qsp = qs(tsp) > q
  //
  // Method:
  // a Newton method is used
  // first guess uses an algorithm provided by John Petch from the UKMO
  // we exclude points where the physical situation is unrealistic
  // e.g. where the temperature is outside the range of validity for the
  //      saturation vapor pressure, or where the water vapor pressure
  //      exceeds the ambient pressure, or the saturation specific humidity is
  //      unrealistic
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void findsp(const amrex::Real& q,
                     const amrex::Real& t,
                     const amrex::Real& p,
                     const bool& use_ice,
                     amrex::Real& tsp,
                     amrex::Real& qsp,
                     int& status) {

   //
   // local variables
   //
   int iter = 8;    // max number of times to iterate the calculation

   amrex::Real es;                   // sat. vapor pressure
   amrex::Real gam;                  // change in sat spec. hum. wrt temperature (times hltalt/cpair)
   amrex::Real dgdt;                 // work variable
   amrex::Real g;                    // work variable
   amrex::Real hltalt;               // lat. heat. of vap.
   amrex::Real qs;                   // spec. hum. of water vapor

   // work variables
   amrex::Real t1, q1, dt, dq;
   amrex::Real qvd;
   amrex::Real r1b, c1, c2, c3;
   constexpr amrex::Real dttol = 1.e-4; // the relative temp error tolerance required to quit the iteration
   constexpr amrex::Real dqtol = 1.e-4; // the relative moisture error tolerance required to quit the iteration
   amrex::Real enin, enout;

   c3 = 287.04*(7.5*log(10.))/Cp_d;

   // Saturation specific humidity at this temperature
   if (use_ice) {
      qsat(t, p, es, qs);
   }
   else {
      qsat_water(t, p, es, qs);
   }

   // make sure a meaningful calculation is possible
   if (p <= 5.*es || qs <= 0. || qs >= 0.5 || t < tmin || t > tmax) {
      status = 1;
      // Keep initial parameters when conditions aren't suitable
      tsp = t;
      qsp = q;
      enin = 1.;
      enout = 1.;
      return;
   }

   // Prepare to iterate
   status = 2;

   // Get initial enthalpy
   if (use_ice)
      calc_hltalt(t,hltalt);
   else
      no_ip_hltalt(t,hltalt);

   enin = tq_enthalpy(t, q, hltalt);

   // make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
   c1 = hltalt*c3;
   c2 = std::pow(t + 36., 2);
   r1b = c2/(c2 + c1*qs);
   qvd = r1b * (q - qs);
   tsp = t + ((hltalt/Cp_d)*qvd);

   // Generate qsp, gam, and enout from tsp.
   if (use_ice)
     qsat(tsp, p, es, qsp, gam, enout);
   else
     qsat_water(tsp, p, es, qsp, gam, enout);

   // iterate on first guess
   for(auto l = 1; l < iter; ++l) {

     g = enin - enout;
     dgdt = -Cp_d * (1 + gam);

     // New tsp
     t1 = tsp - g/dgdt;
     dt = abs(t1 - tsp)/t1;
     tsp = t1;

     // bail out if past end of temperature range
     if ( tsp < tmin ) {
        tsp = tmin;
        // Get latent heat and set qsp to a value
        // that preserves enthalpy.
        if (use_ice)
           calc_hltalt(tsp,hltalt);
        else
           no_ip_hltalt(tsp,hltalt);

        qsp = (enin - Cp_d*tsp)/hltalt;
        enout = tq_enthalpy(tsp, qsp, hltalt);
        status = 4;
        break;
     }

     // Re-generate qsp, gam, and enout from new tsp.
     if (use_ice)
        qsat(tsp, p, es, q1, gam, enout);
     else
        qsat_water(tsp, p, es, q1, gam, enout);

     dq = abs(q1 - qsp)/std::max(q1,1.e-12);
     qsp = q1;

     // if converged at this point, exclude it from more iterations
     if (dt < dttol && dq < dqtol) {
        status = 0;
        break;
     }
   }
   // Test for enthalpy conservation
   if (abs((enin-enout)/(enin+enout)) > 1.e-4) status = 8;
   return;
  }
};
#endif // WATER_VAPOR_SATURATION_H_

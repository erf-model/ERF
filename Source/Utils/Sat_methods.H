// This class contains all methods for estimating
// the saturation vapor pressure of water.
//
// wv_saturation provides specific interfaces and utilities
// based on these formulae.
//
// Typical usage of this module:
//
// Init:
// call wv_sat_methods_init(r8, <constants>, errstring)
//
// Get scheme index from a name string:
// scheme_idx = wv_sat_get_scheme_idx(scheme_name)
// if (.not. wv_sat_valid_idx(scheme_idx)) <throw some error>
//
// Get pressures:
// es = wv_sat_svp_water(t, scheme_idx)
// es = wv_sat_svp_ice(t, scheme_idx)
//
// Use ice/water transition range:
// es = wv_sat_svp_trice(t, ttrice, scheme_idx)
//
// Note that elemental functions cannot be pointed to, nor passed
// as arguments. If you need to do either, it is recommended to
// wrap the function so that it can be given an explicit (non-
// elemental) interface.
#ifndef ERF_SAT_METHOD_H_
#define ERF_SAT_METHOD_H_

#include <string>
#include <vector>
#include <memory>
#include <cmath>

#include "ERF_Constants.H"

class SatMethods {
public:
    // Get saturation specific humidity given pressure and SVP.
    // Specific humidity is limited to range 0-1.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real wv_sat_svp_to_qsat (const amrex::Real& es, const amrex::Real& p) {
        // If pressure is less than SVP, set qs to maximum of 1.
        if ( (p - es) <= 0. )
            return 1.0;
        else
            return epsilo*es / (p - omeps*es);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void wv_sat_qsat_water (const amrex::Real& t, const amrex::Real& p,
                                   amrex::Real es, amrex::Real qs, const int idx = 1) {
        // Purpose:
        //   Calculate SVP over water at a given temperature, and then
        //   calculate and return saturation specific humidity.

        es = wv_sat_svp_water(t, idx);

        qs = wv_sat_svp_to_qsat(es, p);

        // Ensures returned es is consistent with limiters on qs.
        es = std::min(es, p);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void wv_sat_qsat_ice (const amrex::Real& t, const amrex::Real& p,
                                 amrex::Real es, amrex::Real qs, const int idx = 1) {
        // Purpose:
        //   Calculate SVP over ice at a given temperature, and then
        //   calculate and return saturation specific humidity.

        es = wv_sat_svp_ice(t, idx);

        qs = wv_sat_svp_to_qsat(es, p);

        // Ensures returned es is consistent with limiters on qs.
        es = std::min(es, p);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void wv_sat_qsat_trans (const amrex::Real& t, const amrex::Real& p,
                                   amrex::Real es, amrex::Real qs, const int idx = 1) {
        // Purpose:
        //   Calculate SVP over ice at a given temperature, and then
        //   calculate and return saturation specific humidity.

        es = wv_sat_svp_trans(t, idx);

        qs = wv_sat_svp_to_qsat(es, p);

        // Ensures returned es is consistent with limiters on qs.
        es = std::min(es, p);
    }

    // SVP INTERFACE FUNCTIONS
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real wv_sat_svp_water (const amrex::Real& t, const int idx = 1) {
        amrex::Real es;
        switch(idx)
            {
            case GoffGratch:
                es = GoffGratch_svp_water(t); break;
            case MurphyKoop:
                es = MurphyKoop_svp_water(t); break;
            case OldGoffGratch:
                es = OldGoffGratch_svp_water(t); break;
            case Bolton:
                es = Bolton_svp_water(t); break;
            }
        return es;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real wv_sat_svp_ice (const amrex::Real& t, const int idx = 1) {
        amrex::Real es;
        switch(idx)
            {
            case GoffGratch:
                es = GoffGratch_svp_ice(t); break;
            case MurphyKoop:
                es = MurphyKoop_svp_ice(t); break;
            case OldGoffGratch:
                es = OldGoffGratch_svp_ice(t); break;
            case Bolton:
                es = Bolton_svp_water(t); break;
            }
        return es;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real wv_sat_svp_trans (const amrex::Real& t, const int idx = 1) {
        amrex::Real esice;      // Saturation vapor pressure over ice
        amrex::Real weight;     // Intermediate scratch variable for es transition
        amrex::Real es;

        // Water
        if (t >= (tmelt - ttrice))
            es = wv_sat_svp_water(t,idx);
        else
            es = 0.0;

        // Ice
        if (t < tmelt) {
            esice = wv_sat_svp_ice(t,idx);
            if ( (tmelt - t) > ttrice )
                weight = 1.0;
            else
                weight = (tmelt - t)/ttrice;

            es = weight*esice + (1.0 - weight)*es;
        }
        return es;
    }


    // Goff & Gratch (1946)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real GoffGratch_svp_water (const amrex::Real& t) {
        // uncertain below -70 C
        return pow(10., (-7.90298*(tboil/t-1.)+
                         5.02808*std::log10(tboil/t)-
                         1.3816e-7*(pow(10., (11.344*(1.-t/tboil)))-1.)+
                         8.1328e-3*(pow(10., (-3.49149*(tboil/t-1.)))-1.)+
                         std::log10(1013.246)))*100.;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real GoffGratch_svp_ice (const amrex::Real& t) {
        // good down to -100 C
        return pow(10., (-9.09718*(h2otrip/t-1.)-3.56654*
                         log10(h2otrip/t)+0.876793*(1.-t/h2otrip)+
                         log10(6.1071)))*100.;
    }

    // Murphy & Koop (2005)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real MurphyKoop_svp_water (const amrex::Real& t) {
        // (good for 123 < T < 332 K)
        return exp(54.842763 - (6763.22 / t) - (4.210 * log(t)) +
                   (0.000367 * t) + (tanh(0.0415 * (t - 218.8)) *
                                     (53.878 - (1331.22 / t) - (9.44523 * log(t)) +
                                      0.014025 * t)));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real MurphyKoop_svp_ice (const amrex::Real& t) {
        // (good down to 110 K)
        return exp(9.550426 - (5723.265 / t) + (3.53068 * log(t))
                   - (0.00728332 * t));
    }


    // Old CAM implementation, also labelled Goff & Gratch (1946)

    // The water formula differs only due to compiler-dependent order of
    // operations, so differences are roundoff level, usually 0.

    // The ice formula gives fairly close answers to the current
    // implementation, but has been rearranged, and uses the
    // 1 atm melting point of water as the triple point.
    // Differences are thus small but above roundoff.

    // A curious fact: although using the melting point of water was
    // probably a mistake, it mildly improves accuracy for ice svp,
    // since it compensates for a systematic error in Goff & Gratch.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real OldGoffGratch_svp_water (const amrex::Real& t) {
        auto ps = 1013.246;
        auto e1 = 11.344*(1.0 - t/tboil);
        auto e2 = -3.49149*(tboil/t - 1.0);
        auto f1 = -7.90298*(tboil/t - 1.0);
        auto f2 = 5.02808*log10(tboil/t);
        auto f3 = -1.3816*(pow(10.0, e1) - 1.0)/10000000.0;
        auto f4 = 8.1328*(pow(10.0, e2) - 1.0)/1000.0;
        auto f5 = log10(ps);
        auto f  = f1 + f2 + f3 + f4 + f5;
        return (pow(10.0, f))*100.0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real OldGoffGratch_svp_ice (const amrex::Real& t) {
        auto term1 = 2.01889049/(tmelt/t);
        auto term2 = 3.56654*log(tmelt/t);
        auto term3 = 20.947031*(tmelt/t);
        return 575.185606e10*exp(-(term1 + term2 + term3));
    }

    // Bolton (1980)
    // zm_conv deep convection scheme contained this SVP calculation.
    // It appears to be from D. Bolton, 1980, Monthly Weather Review.
    // Unlike the other schemes, no distinct ice formula is associated
    // with it. (However, a Bolton ice formula exists in CLUBB.)

    // The original formula used degrees C, but this function
    // takes Kelvin and internally converts.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real Bolton_svp_water (const amrex::Real& t) {
        constexpr auto c1 = 611.2;
        constexpr auto c2 = 17.67;
        constexpr auto c3 = 243.5;
        return c1*exp( (c2*(t - tmelt))/((t - tmelt)+c3) );
    }

    // private data
private:
    // Indices representing individual schemes
    enum Type {
               Invalid = -1,
               OldGoffGratch = 0,
               GoffGratch = 1,
               MurphyKoop = 2,
               Bolton = 3
    };
};
#endif

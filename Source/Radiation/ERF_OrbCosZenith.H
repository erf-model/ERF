#ifndef ERF_ORB_COS_ZENITH_H
#define ERF_ORB_COS_ZENITH_H

#include <vector>
#include <cmath>

#include <ERF_Constants.H>

typedef double real;

real
orbital_cos_zenith (real& jday,
                    real& lat,
                    real& lon,
                    real& declin,
                    real dt_avg = -1.,
                    real uniform_angle = -1.,
                    real constant_zenith_angle_deg = -1.);


real
orbital_avg_cos_zenith (real& jday,
                        real& lat,
                        real& lon,
                        real& declin,
                        real& dt_avg);


void
orbital_params (int& iyear_AD,
                real& eccen,
                real& obliq,
                real& mvelp,
                real& obliqr,
                real& lambm0,
                real& mvelpp);


void
orbital_decl (real& calday,
              real& eccen,
              real& mvelpp,
              real& lambm0,
              real& obliqr,
              real& delta,
              real& eccf);
#endif

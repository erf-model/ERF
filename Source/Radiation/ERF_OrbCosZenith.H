#ifndef ERF_ORB_COS_ZENITH_H
#define ERF_ORB_COS_ZENITH_H

#include <math.h>
#include <ERF_Constants.H>

amrex::Real
orbital_cos_zenith (amrex::Real& jday,
                    amrex::Real& lat,
                    amrex::Real& lon,
                    amrex::Real& declin,
                    amrex::Real dt_avg = -1.,
                    amrex::Real uniform_angle = -1.,
                    amrex::Real constant_zenith_angle_deg = -1.);


amrex::Real
orbital_avg_cos_zenith (amrex::Real& jday,
                        amrex::Real& lat,
                        amrex::Real& lon,
                        amrex::Real& declin,
                        amrex::Real& dt_avg);


void
orbital_params (int& iyear_AD,
                amrex::Real& eccen,
                amrex::Real& obliq,
                amrex::Real& mvelp,
                amrex::Real& obliqr,
                amrex::Real& lambm0,
                amrex::Real& mvelpp);


void
orbital_decl (amrex::Real& calday,
              amrex::Real& eccen,
              amrex::Real& mvelpp,
              amrex::Real& lambm0,
              amrex::Real& obliqr,
              amrex::Real& delta,
              amrex::Real& eccf);
#endif

#ifndef ERF_SYNTHETIC_H_
#define ERF_SYNTHETIC_H_
// A tiny ERF-like header for the schema generator's tests.
AMREX_ENUM(MoistureType, None, Kessler, SAM);

struct SyntheticParams {
    amrex::Real cfl = 0.8;
    int nsteps = 10;
    amrex::Real eye_latitude = -1e10, eye_longitude = -1e10;   // two declarators, one line
};
#endif

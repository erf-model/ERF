// MorrisonFortranInterface.h

#ifndef MORRISON_FORTRAN_INTERFACE_H
#define MORRISON_FORTRAN_INTERFACE_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Initialize the Morrison microphysics scheme
 *
 * @param morr_rimed_ice: 0 for graupel, 1 for hail configuration
 */
void morr_two_moment_init_c(int morr_rimed_ice, int morr_noice);

/**
 * Set the constant droplet concentration (NDCNST) for Morrison microphysics
 *
 * @param ndcnst_in: Droplet number concentration in cm^-3
 */
void set_morrison_ndcnst_c(double ndcnst_in);

/**
 * C interface to the Fortran MP_MORR_TWO_MOMENT routine.
 *
 * @note: All array dimensions must match between C++ and Fortran calls.
 * @note: Arrays must be allocated in C++ before passing to this routine.
 */
void mp_morr_two_moment_c(
    int itimestep,
    // 3D arrays for hydrometeors and state variables
    double* th, double* qv, double* qc, double* qr, double* qi, double* qs, double* qg,
    double* ni, double* ns, double* nr, double* ng,
    // 3D arrays for thermodynamics
    double* rho, double* pii, double* p, double dt_in, double* dz, double* w,
    // 2D arrays for precipitation
    double* rainnc, double* rainncv, double* sr,
    double* snownc, double* snowncv, double* graupelnc, double* graupelncv,
    // Radar reflectivity
    double* refl_10cm, bool diagflag, int do_radar_ref,
    // Cumulus tendencies
    double* qrcuten, double* qscuten, double* qicuten,
    // WRF-Chem
    bool f_qndrop, double* qndrop,
    // Terrain height
    double* ht,
    // Domain dimensions
    int ids, int ide, int jds, int jde, int kds, int kde,
    int ims, int ime, int jms, int jme, int kms, int kme,
    int its, int ite, int jts, int jte, int kts, int kte,
    // Optional WRF-Chem related
    bool wetscav_on, double* rainprod, double* evapprod,
    double* qlsink, double* precr, double* preci, double* precs, double* precg
);

#ifdef __cplusplus
}
#endif

#endif // MORRISON_FORTRAN_INTERFACE_H

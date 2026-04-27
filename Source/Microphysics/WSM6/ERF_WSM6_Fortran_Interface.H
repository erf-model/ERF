#ifndef ERF_WSM6_FORTRAN_INTERFACE_H
#define ERF_WSM6_FORTRAN_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

// Phase-1 C bindings for WSM6 Fortran bridge.
// These keep scope to required single-moment arguments only.
void mp_wsm6_init_c(double den0, double denr, double dens, double cl, double cpv, int hail_opt);

void mp_wsm6_run_c(
    double* t,
    double* qv, double* qc, double* qi, double* qr, double* qs, double* qg,
    double* den, double* p, double* delz,
    double delt, double g, double cpd, double cpv, double rd, double rv,
    double t0c, double ep1, double ep2, double qmin, double xls,
    double xlv0, double xlf0, double den0, double denr, double cliq, double cice, double psat,
    double* rain, double* rainncv, double* sr,
    double* snow, double* snowncv, double* graupel, double* graupelncv,
    int ims, int ime, int jms, int jme, int kms, int kme,
    int its, int ite, int jts, int jte, int kts, int kte);

#ifdef __cplusplus
}
#endif

#endif

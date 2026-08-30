#ifndef ERF_WDM6_FORTRAN_INTERFACE_H
#define ERF_WDM6_FORTRAN_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

// WDM6 initialization
// Adds ccn0 parameter compared to WSM6
void mp_wdm6_init_c(double den0, double denr, double dens,
                    double cl, double cpv, double ccn0, int hail_opt);

// WDM6 run interface
// Adds nn, nc, nr arrays and ccn0, xland parameters compared to WSM6
void mp_wdm6_run_c(
    double* t,
    double* qv, double* qc, double* qi, double* qr, double* qs, double* qg,
    double* nn, double* nc, double* nr,  // WDM6: number concentrations
    double* den, double* p, double* delz,
    double delt, double g, double cpd, double cpv, double rd, double rv,
    double t0c, double ep1, double ep2, double qmin, double xls,
    double xlv0, double xlf0, double den0, double denr, double cliq, double cice, double psat,
    double ccn0,     // WDM6: CCN concentration
    double* xland,   // WDM6: land/water mask
    double* rain, double* rainncv, double* sr,
    double* snow, double* snowncv, double* graupel, double* graupelncv,
    int ims, int ime, int jms, int jme, int kms, int kme,
    int its, int ite, int jts, int jte, int kts, int kte,
    int microphysics_debug, int diag_i_dbg, int diag_j_dbg);

#ifdef __cplusplus
}
#endif

#endif

#ifndef ERF_WSM6_FORTRAN_INTERFACE_H
#define ERF_WSM6_FORTRAN_INTERFACE_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

// Phase-1 C bindings for WSM6 Fortran bridge.
// These keep scope to required single-moment arguments only.
void mp_wsm6_init_c(amrex_real den0, amrex_real denr, amrex_real dens, amrex_real cl, amrex_real cpv, int hail_opt);

void mp_wsm6_run_c(
    amrex_real* t,
    amrex_real* qv, amrex_real* qc, amrex_real* qi, amrex_real* qr, amrex_real* qs, amrex_real* qg,
    amrex_real* den, amrex_real* p, amrex_real* delz,
    amrex_real delt, amrex_real g, amrex_real cpd, amrex_real cpv, amrex_real rd, amrex_real rv,
    amrex_real t0c, amrex_real ep1, amrex_real ep2, amrex_real qmin, amrex_real xls,
    amrex_real xlv0, amrex_real xlf0, amrex_real den0, amrex_real denr, amrex_real cliq, amrex_real cice, amrex_real psat,
    amrex_real* rain, amrex_real* rainncv, amrex_real* sr,
    amrex_real* snow, amrex_real* snowncv, amrex_real* graupel, amrex_real* graupelncv,
    int ims, int ime, int jms, int jme, int kms, int kme,
    int its, int ite, int jts, int jte, int kts, int kte,
    int do_cond, int microphysics_debug, int diag_i_dbg, int diag_j_dbg);

#ifdef __cplusplus
}
#endif

#endif

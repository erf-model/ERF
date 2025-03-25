MODULE mp_morr_two_moment_isohelper
  USE ISO_C_BINDING
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE mp_morr_two_moment_c(itimestep, &
                th, qv, qc, qr, qi, qs, qg, ni, ns, nr, ng, &
                rho, pii, p, dt_in, dz, ht, w, &
                rainnc, rainncv, sr, &
                snownc, snowncv, graupelnc, graupelncv, &
                refl_10cm, diagflag, do_radar_ref, &
                qrcuten, qscuten, qicuten, &
                f_qndrop, qndrop, &
                ids, ide, jds, jde, kds, kde, &
                ims, ime, jms, jme, kms, kme, &
                its, ite, jts, jte, kts, kte, &
                wetscav_on, rainprod, evapprod, &
                qlsink, precr, preci, precs, precg) &
                BIND(C, name="mp_morr_two_moment_c")

    ! Define C interoperable types
    INTEGER(C_INT), VALUE, INTENT(IN) :: itimestep
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: th, qv, qc, qr, qi, qs, qg, ni, ns, nr, ng
    REAL(C_DOUBLE), INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: rho, pii, p, dz, w
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: dt_in
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, jms:jme) :: rainnc, rainncv, sr
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, jms:jme) :: snownc, snowncv, graupelnc, graupelncv
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: refl_10cm
    LOGICAL(C_BOOL), VALUE, INTENT(IN) :: diagflag
    INTEGER(C_INT), VALUE, INTENT(IN) :: do_radar_ref
    REAL(C_DOUBLE), INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: qrcuten, qscuten, qicuten
    LOGICAL(C_BOOL), VALUE, INTENT(IN) :: f_qndrop
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: qndrop
    REAL(C_DOUBLE), INTENT(IN), DIMENSION(ims:ime, jms:jme) :: ht

    ! Domain dimensions
    INTEGER(C_INT), VALUE, INTENT(IN) :: ids, ide, jds, jde, kds, kde
    INTEGER(C_INT), VALUE, INTENT(IN) :: ims, ime, jms, jme, kms, kme
    INTEGER(C_INT), VALUE, INTENT(IN) :: its, ite, jts, jte, kts, kte

    ! Optional arguments
    LOGICAL(C_BOOL), VALUE, INTENT(IN) :: wetscav_on
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: rainprod, evapprod
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: qlsink, precr, preci, precs, precg

    ! Call the original routine
    CALL MP_MORR_TWO_MOMENT(itimestep, &
                th, qv, qc, qr, qi, qs, qg, ni, ns, nr, ng, &
                rho, pii, p, dt_in, dz, ht, w, &
                rainnc, rainncv, sr, &
                snownc, snowncv, graupelnc, graupelncv, &
                refl_10cm, diagflag, do_radar_ref, &
                qrcuten, qscuten, qicuten, &
                f_qndrop, qndrop, &
                ids, ide, jds, jde, kds, kde, &
                ims, ime, jms, jme, kms, kme, &
                its, ite, jts, jte, kts, kte, &
                wetscav_on, rainprod, evapprod, &
                qlsink, precr, preci, precs, precg)

  END SUBROUTINE mp_morr_two_moment_c

END MODULE mp_morr_two_moment_isohelper

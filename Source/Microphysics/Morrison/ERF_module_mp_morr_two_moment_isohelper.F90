MODULE mp_morr_two_moment_isohelper
  USE ISO_C_BINDING
  USE MODULE_MP_MORR_TWO_MOMENT, ONLY: MP_MORR_TWO_MOMENT, MORR_TWO_MOMENT_INIT
  IMPLICIT NONE

  CONTAINS

  ! Initialize the Morrison microphysics scheme
  SUBROUTINE morr_two_moment_init_c(morr_rimed_ice, morr_noice) BIND(C, name="morr_two_moment_init_c")
    INTEGER(C_INT), VALUE, INTENT(IN) :: morr_rimed_ice
    INTEGER(C_INT), VALUE, INTENT(IN) :: morr_noice
    CALL MORR_TWO_MOMENT_INIT(morr_rimed_ice, morr_noice)
  END SUBROUTINE morr_two_moment_init_c

  SUBROUTINE mp_morr_two_moment_c(itimestep, &
                th, qv, qc, qr, qi, qs, qg, ni, ns, nr, ng, &
                rho, pii, p, dt_in, dz, w, &
                rainnc, rainncv, sr, &
                snownc, snowncv, graupelnc, graupelncv, &
                refl_10cm, diagflag, do_radar_ref, &
                qrcuten, qscuten, qicuten, &
                f_qndrop, qndrop, ht, &
                ids, ide, jds, jde, kds, kde, &
                ims, ime, jms, jme, kms, kme, &
                its, ite, jts, jte, kts, kte, &
                wetscav_on, rainprod, evapprod, &
                qlsink, precr, preci, precs, precg) &
                BIND(C, name="mp_morr_two_moment_c")

    ! Define C interoperable types
    INTEGER(C_INT), VALUE, INTENT(IN) :: itimestep
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, jms:jme, kms:kme) :: th, qv, qc, qr, qi, qs, qg, ni, ns, nr, ng
    REAL(C_DOUBLE), INTENT(IN), DIMENSION(ims:ime, jms:jme, kms:kme) :: rho, pii, p, dz, w
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: dt_in
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, jms:jme, kms:kme) :: rainnc, snownc, graupelnc
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, jms:jme) :: rainncv, sr,  snowncv, graupelncv
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, jms:jme, kms:kme) :: refl_10cm
    LOGICAL(C_BOOL), VALUE, INTENT(IN) :: diagflag
    INTEGER(C_INT), VALUE, INTENT(IN) :: do_radar_ref
    REAL(C_DOUBLE), INTENT(IN), DIMENSION(ims:ime, jms:jme, kms:kme) :: qrcuten, qscuten, qicuten
    LOGICAL(C_BOOL), VALUE, INTENT(IN) :: f_qndrop
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, jms:jme, kms:kme) :: qndrop
    REAL(C_DOUBLE), INTENT(IN), DIMENSION(ims:ime, jms:jme) :: ht

    ! Domain dimensions
    INTEGER(C_INT), VALUE, INTENT(IN) :: ids, ide, jds, jde, kds, kde
    INTEGER(C_INT), VALUE, INTENT(IN) :: ims, ime, jms, jme, kms, kme
    INTEGER(C_INT), VALUE, INTENT(IN) :: its, ite, jts, jte, kts, kte

    ! Optional arguments
    LOGICAL(C_BOOL), VALUE, INTENT(IN) :: wetscav_on
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, jms:jme, kms:kme) :: rainprod, evapprod
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, jms:jme, kms:kme) :: qlsink, precr, preci, precs, precg

    ! Convert C_BOOL to Fortran logical
    LOGICAL :: diag_flag_f, f_qndrop_f, wetscav_on_f

    ! Convert C types to Fortran types
    diag_flag_f = diagflag
    f_qndrop_f = f_qndrop
    wetscav_on_f = wetscav_on

    CALL MP_MORR_TWO_MOMENT(itimestep, &
                th, qv, qc, qr, qi, qs, qg, ni, ns, nr, ng, &
                rho, pii, p, dt_in, dz, ht, w, &
                rainnc, rainncv, sr, &
                snownc, snowncv, graupelnc, graupelncv, &
                refl_10cm, diag_flag_f, do_radar_ref, &
                qrcuten, qscuten, qicuten, &
                f_qndrop_f, qndrop, &
                ids, ide, jds, jde, kds, kde, &
                ims, ime, jms, jme, kms, kme, &
                its, ite, jts, jte, kts, kte, &
                wetscav_on_f, rainprod, evapprod, &
                qlsink, precr, preci, precs, precg)

  END SUBROUTINE mp_morr_two_moment_c

END MODULE mp_morr_two_moment_isohelper

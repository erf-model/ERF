module mp_wsm6_isohelper
  use iso_c_binding
  use mp_wsm6, only: mp_wsm6_init, mp_wsm6_run
  implicit none

contains

  subroutine mp_wsm6_init_c(den0, denr, dens, cl, cpv, hail_opt) bind(C, name="mp_wsm6_init_c")
    real(c_double), value, intent(in) :: den0, denr, dens, cl, cpv
    integer(c_int), value, intent(in) :: hail_opt
    character(len=256) :: errmsg
    integer(c_int) :: errflg

    call mp_wsm6_init(den0, denr, dens, cl, cpv, hail_opt, errmsg, errflg)
    if (errflg /= 0_c_int) then
      write(*,'(A,1X,A)') 'mp_wsm6_init_c error:', trim(errmsg)
      stop 1
    end if
  end subroutine mp_wsm6_init_c

  subroutine mp_wsm6_run_c(t, qv, qc, qi, qr, qs, qg, den, p, delz, &
                           delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, xls, &
                           xlv0, xlf0, den0, denr, cliq, cice, psat, &
                           rain, rainncv, sr, snow, snowncv, graupel, graupelncv, &
                           ids, ide, jds, jde, kds, kde, &
                           ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte) bind(C, name="mp_wsm6_run_c")
    integer(c_int), value, intent(in) :: ids, ide, jds, jde, kds, kde
    integer(c_int), value, intent(in) :: ims, ime, jms, jme, kms, kme
    integer(c_int), value, intent(in) :: its, ite, jts, jte, kts, kte
    real(c_double), intent(inout), dimension(ims:ime, jms:jme, kms:kme) :: t, qv, qc, qi, qr, qs, qg
    real(c_double), intent(in), dimension(ims:ime, jms:jme, kms:kme) :: den, p, delz
    real(c_double), value, intent(in) :: delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, xls
    real(c_double), value, intent(in) :: xlv0, xlf0, den0, denr, cliq, cice, psat
    real(c_double), intent(inout), dimension(ims:ime, jms:jme) :: rain, rainncv, sr
    real(c_double), intent(inout), dimension(ims:ime, jms:jme) :: snow, snowncv, graupel, graupelncv
    character(len=256) :: errmsg
    integer(c_int) :: errflg

    call mp_wsm6_run(t, qv, qc, qi, qr, qs, qg, den, p, delz, &
                     delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, xls, &
                     xlv0, xlf0, den0, denr, cliq, cice, psat, &
                     rain, rainncv, sr, snow, snowncv, graupel, graupelncv, &
                     ids, ide, jds, jde, kds, kde, &
                     ims, ime, jms, jme, kms, kme, &
                     its, ite, jts, jte, kts, kte, errmsg, errflg)
    if (errflg /= 0_c_int) then
      write(*,'(A,1X,A)') 'mp_wsm6_run_c error:', trim(errmsg)
      stop 1
    end if
  end subroutine mp_wsm6_run_c

end module mp_wsm6_isohelper

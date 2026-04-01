module mp_wsm6
  use iso_c_binding
  implicit none
  private
  public :: mp_wsm6_init, mp_wsm6_run

contains

  subroutine mp_wsm6_init(den0, denr, dens, cl, cpv, hail_opt, errmsg, errflg)
    real(c_double), intent(in) :: den0, denr, dens, cl, cpv
    integer(c_int), intent(in) :: hail_opt
    character(len=*), intent(out) :: errmsg
    integer(c_int), intent(out) :: errflg

    if (den0 <= 0.0_c_double .or. denr <= 0.0_c_double .or. dens <= 0.0_c_double) then
      errmsg = 'mp_wsm6_init: invalid reference density input'
      errflg = 1_c_int
      return
    end if

    if (cl <= 0.0_c_double .or. cpv <= 0.0_c_double) then
      errmsg = 'mp_wsm6_init: invalid heat capacity input'
      errflg = 2_c_int
      return
    end if

    if (hail_opt < 0_c_int .or. hail_opt > 1_c_int) then
      errmsg = 'mp_wsm6_init: hail_opt must be 0 or 1'
      errflg = 3_c_int
      return
    end if

    errmsg = 'mp_wsm6_init OK (ERF tracked stub)'
    errflg = 0_c_int
  end subroutine mp_wsm6_init

  subroutine mp_wsm6_run(t, qv, qc, qi, qr, qs, qg, den, p, delz, &
                         delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, xls, &
                         xlv0, xlf0, den0, denr, cliq, cice, psat, &
                         rain, rainncv, sr, snow, snowncv, graupel, graupelncv, &
                         ids, ide, jds, jde, kds, kde, &
                         ims, ime, jms, jme, kms, kme, &
                         its, ite, jts, jte, kts, kte, errmsg, errflg)
    real(c_double), intent(inout), dimension(ims:ime, jms:jme, kms:kme) :: t, qv, qc, qi, qr, qs, qg
    real(c_double), intent(in), dimension(ims:ime, jms:jme, kms:kme) :: den, p, delz
    real(c_double), intent(in) :: delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, xls
    real(c_double), intent(in) :: xlv0, xlf0, den0, denr, cliq, cice, psat
    real(c_double), intent(inout), dimension(ims:ime, jms:jme) :: rain, rainncv, sr, snow, snowncv, graupel, graupelncv
    integer(c_int), intent(in) :: ids, ide, jds, jde, kds, kde
    integer(c_int), intent(in) :: ims, ime, jms, jme, kms, kme
    integer(c_int), intent(in) :: its, ite, jts, jte, kts, kte
    character(len=*), intent(out) :: errmsg
    integer(c_int), intent(out) :: errflg

    integer :: i, j, k

    if (delt <= 0.0_c_double .or. g <= 0.0_c_double .or. cpd <= 0.0_c_double .or. cpv <= 0.0_c_double) then
      errmsg = 'mp_wsm6_run: invalid scalar physics input'
      errflg = 4_c_int
      return
    end if

    do k = kts, kte
      do j = jts, jte
        do i = its, ite
          qv(i,j,k) = max(qv(i,j,k), qmin)
          qc(i,j,k) = max(qc(i,j,k), 0.0_c_double)
          qi(i,j,k) = max(qi(i,j,k), 0.0_c_double)
          qr(i,j,k) = max(qr(i,j,k), 0.0_c_double)
          qs(i,j,k) = max(qs(i,j,k), 0.0_c_double)
          qg(i,j,k) = max(qg(i,j,k), 0.0_c_double)
        end do
      end do
    end do

    rainncv(its:ite,jts:jte) = 0.0_c_double
    sr(its:ite,jts:jte) = 0.0_c_double
    snowncv(its:ite,jts:jte) = 0.0_c_double
    graupelncv(its:ite,jts:jte) = 0.0_c_double

    errmsg = 'mp_wsm6_run OK (ERF tracked stub)'
    errflg = 0_c_int
  end subroutine mp_wsm6_run

end module mp_wsm6

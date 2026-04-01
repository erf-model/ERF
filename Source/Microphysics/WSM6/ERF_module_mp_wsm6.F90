module mp_wsm6
  use iso_c_binding
  implicit none
  private
  public :: mp_wsm6_init, mp_wsm6_run

  real(c_double), parameter :: dtcldcr = 120.0_c_double
  real(c_double), parameter :: n0r = 8.0e6_c_double
  real(c_double), parameter :: avtr = 841.9_c_double
  real(c_double), parameter :: bvtr = 0.8_c_double
  real(c_double), parameter :: r0 = 0.8e-5_c_double
  real(c_double), parameter :: peaut = 0.55_c_double
  real(c_double), parameter :: xncr = 3.0e8_c_double
  real(c_double), parameter :: xmyu = 1.718e-5_c_double
  real(c_double), parameter :: avts = 11.72_c_double
  real(c_double), parameter :: bvts = 0.41_c_double
  real(c_double), parameter :: lamdarmax = 8.0e4_c_double
  real(c_double), parameter :: lamdasmax = 1.0e5_c_double
  real(c_double), parameter :: dimax = 500.0e-6_c_double
  real(c_double), parameter :: eacrc = 1.0_c_double
  real(c_double), parameter :: n0s = 2.0e6_c_double

  real(c_double), save :: qc0, qck1, &
                          bvtr1, bvtr2, bvtr3, bvtr4, bvtr6, g1pbr, g3pbr, g4pbr, g6pbr, g5pbro2, &
                          pvtr, eacrr, pacrr, precr1, precr2, roqimax, &
                          bvts1, bvts2, bvts3, bvts4, g1pbs, g3pbs, g4pbs, g5pbso2, &
                          n0g, avtg, bvtg, deng, lamdagmax, pvts, pacrs, precs1, precs2, &
                          pidn0r, pidn0s, pidnc, xlv1, pacrc, pi, &
                          bvtg1, bvtg2, bvtg3, bvtg4, g1pbg, g3pbg, g4pbg, g5pbgo2, pvtg, pacrg, &
                          precg1, precg2, pidn0g, rslopermax, rslopesmax, rslopegmax, &
                          rsloperbmax, rslopesbmax, rslopegbmax, rsloper2max, rslopes2max, rslopeg2max, &
                          rsloper3max, rslopes3max, rslopeg3max

contains

  pure real(c_double) function rgmma(x)
    real(c_double), intent(in) :: x
    rgmma = exp(log_gamma(x))
  end function rgmma

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

    select case (hail_opt)
    case (1_c_int)
      n0g = 4.0e4_c_double
      deng = 700.0_c_double
      avtg = 285.0_c_double
      bvtg = 0.8_c_double
      lamdagmax = 2.0e4_c_double
    case (0_c_int)
      n0g = 4.0e6_c_double
      deng = 500.0_c_double
      avtg = 330.0_c_double
      bvtg = 0.8_c_double
      lamdagmax = 6.0e4_c_double
    case default
      errmsg = 'mp_wsm6_init: hail_opt must be 0 or 1'
      errflg = 3_c_int
      return
    end select

    ! Canonical-style coefficient initialization (phase-1, no radar coupling yet).
    pi = 4.0_c_double * atan(1.0_c_double)
    xlv1 = cl - cpv
    qc0 = (4.0_c_double / 3.0_c_double) * pi * denr * r0**3 * xncr / den0
    qck1 = 0.104_c_double * 9.8_c_double * peaut / (xncr * denr)**(1.0_c_double / 3.0_c_double) / &
           xmyu * den0**(4.0_c_double / 3.0_c_double)
    pidnc = pi * denr / 6.0_c_double

    bvtr1 = 1.0_c_double + bvtr
    bvtr2 = 2.5_c_double + 0.5_c_double * bvtr
    bvtr3 = 3.0_c_double + bvtr
    bvtr4 = 4.0_c_double + bvtr
    bvtr6 = 6.0_c_double + bvtr
    g1pbr = rgmma(bvtr1)
    g3pbr = rgmma(bvtr3)
    g4pbr = rgmma(bvtr4)
    g6pbr = rgmma(bvtr6)
    g5pbro2 = rgmma(bvtr2)
    pvtr = avtr * g4pbr / 6.0_c_double
    eacrr = 1.0_c_double
    pacrr = pi * n0r * avtr * g3pbr * 0.25_c_double * eacrr
    precr1 = 2.0_c_double * pi * n0r * 0.78_c_double
    precr2 = 2.0_c_double * pi * n0r * 0.31_c_double * sqrt(avtr) * g5pbro2
    roqimax = 2.08e22_c_double * dimax**8

    bvts1 = 1.0_c_double + bvts
    bvts2 = 2.5_c_double + 0.5_c_double * bvts
    bvts3 = 3.0_c_double + bvts
    bvts4 = 4.0_c_double + bvts
    g1pbs = rgmma(bvts1)
    g3pbs = rgmma(bvts3)
    g4pbs = rgmma(bvts4)
    g5pbso2 = rgmma(bvts2)
    pvts = avts * g4pbs / 6.0_c_double
    pacrs = pi * n0s * avts * g3pbs * 0.25_c_double
    precs1 = 4.0_c_double * n0s * 0.65_c_double
    precs2 = 4.0_c_double * n0s * 0.44_c_double * sqrt(avts) * g5pbso2
    pidn0r = pi * denr * n0r
    pidn0s = pi * dens * n0s
    pacrc = pi * n0s * avts * g3pbs * 0.25_c_double * eacrc

    bvtg1 = 1.0_c_double + bvtg
    bvtg2 = 2.5_c_double + 0.5_c_double * bvtg
    bvtg3 = 3.0_c_double + bvtg
    bvtg4 = 4.0_c_double + bvtg
    g1pbg = rgmma(bvtg1)
    g3pbg = rgmma(bvtg3)
    g4pbg = rgmma(bvtg4)
    g5pbgo2 = rgmma(bvtg2)
    pacrg = pi * n0g * avtg * g3pbg * 0.25_c_double
    pvtg = avtg * g4pbg / 6.0_c_double
    precg1 = 2.0_c_double * pi * n0g * 0.78_c_double
    precg2 = 2.0_c_double * pi * n0g * 0.31_c_double * sqrt(avtg) * g5pbgo2
    pidn0g = pi * deng * n0g

    rslopermax = 1.0_c_double / lamdarmax
    rslopesmax = 1.0_c_double / lamdasmax
    rslopegmax = 1.0_c_double / lamdagmax
    rsloperbmax = rslopermax**bvtr
    rslopesbmax = rslopesmax**bvts
    rslopegbmax = rslopegmax**bvtg
    rsloper2max = rslopermax * rslopermax
    rslopes2max = rslopesmax * rslopesmax
    rslopeg2max = rslopegmax * rslopegmax
    rsloper3max = rsloper2max * rslopermax
    rslopes3max = rslopes2max * rslopesmax
    rslopeg3max = rslopeg2max * rslopegmax

    errmsg = 'mp_wsm6_init OK (ERF coeff init)'
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

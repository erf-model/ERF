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
  real(c_double), parameter :: n0smax = 1.0e11_c_double
  real(c_double), parameter :: alpha = 0.12_c_double
  real(c_double), parameter :: qcrmin = 1.0e-9_c_double

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

  pure real(c_double) function fpvs(t, ice, rd, rv, cvap, cliq, cice, hvap, hsub, psat, t0c)
    real(c_double), intent(in) :: t, rd, rv, cvap, cliq, cice, hvap, hsub, psat, t0c
    integer, intent(in) :: ice
    real(c_double) :: tr, dldt

    tr = t0c / t
    if (ice == 0) then
      dldt = cvap - cliq
      fpvs = psat * tr**(dldt/rv) * exp((hvap - dldt * t0c) / rv * (1.0_c_double / t0c - 1.0_c_double / t))
    else
      dldt = cvap - cice
      fpvs = psat * tr**(dldt/rv) * exp((hsub - dldt * t0c) / rv * (1.0_c_double / t0c - 1.0_c_double / t))
    end if
  end function fpvs

  pure real(c_double) function diffus(t, den)
    real(c_double), intent(in) :: t, den
    diffus = 8.794e-5_c_double * exp(log(max(t, 1.0_c_double)) * 1.81_c_double) / max(den, 1.0e-12_c_double)
  end function diffus

  pure real(c_double) function viscos(t, den)
    real(c_double), intent(in) :: t, den
    viscos = 1.496e-6_c_double * (t * sqrt(max(t, 1.0_c_double))) / (t + 120.0_c_double) / max(den, 1.0e-12_c_double)
  end function viscos

  pure real(c_double) function xka(t, den)
    real(c_double), intent(in) :: t, den
    xka = 1.414e3_c_double * viscos(t, den) * den
  end function xka

  pure real(c_double) function diffac(a, p, t, den, qsat, rv)
    real(c_double), intent(in) :: a, p, t, den, qsat, rv
    diffac = den * a * a / (max(xka(t, den), 1.0e-12_c_double) * rv * t * t) + &
             1.0_c_double / (max(qsat, 1.0e-12_c_double) * max(diffus(t, p), 1.0e-12_c_double))
  end function diffac

  pure real(c_double) function venfac(p, t, den, den0)
    real(c_double), intent(in) :: p, t, den, den0
    venfac = exp(log(max(viscos(t, den) / max(diffus(t, p), 1.0e-12_c_double), 1.0e-12_c_double)) / 3.0_c_double) / &
             sqrt(max(viscos(t, den), 1.0e-12_c_double)) * sqrt(sqrt(max(den0 / max(den, 1.0e-12_c_double), 1.0e-12_c_double)))
  end function venfac

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

    integer :: i, j, k, loop, loops
    real(c_double) :: dtcld, es, qsat, supsat, dq, cpm, xlv, xlf
    real(c_double) :: freeze_q
    real(c_double) :: denfac, rslope_r, rslopeb_r, rslope2_r, rslope3_r, lamdar
    real(c_double) :: rslope_s, rslopeb_s, rslope2_s, rslope3_s, lamdas
    real(c_double) :: rslope_g, rslopeb_g, rslope2_g, rslope3_g, lamdag
    real(c_double) :: n0sfac, supcol
    real(c_double) :: work1_r, work2_r, coeres, satdt, praut, pracw, prevp
    real(c_double) :: psmlt, pgmlt
    real(c_double) :: rain_flux, snow_flux, graup_flux

    if (delt <= 0.0_c_double .or. g <= 0.0_c_double .or. cpd <= 0.0_c_double .or. cpv <= 0.0_c_double) then
      errmsg = 'mp_wsm6_run: invalid scalar physics input'
      errflg = 4_c_int
      return
    end if

    rainncv(its:ite,jts:jte) = 0.0_c_double
    sr(its:ite,jts:jte) = 0.0_c_double
    snowncv(its:ite,jts:jte) = 0.0_c_double
    graupelncv(its:ite,jts:jte) = 0.0_c_double

    loops = max(nint(delt / dtcldcr), 1)
    dtcld = delt / dble(loops)
    if (delt <= dtcldcr) dtcld = delt

    do loop = 1, loops
      do k = kts, kte
        do j = jts, jte
          do i = its, ite
            qv(i,j,k) = max(qv(i,j,k), qmin)
            qc(i,j,k) = max(qc(i,j,k), 0.0_c_double)
            qi(i,j,k) = max(qi(i,j,k), 0.0_c_double)
            qr(i,j,k) = max(qr(i,j,k), 0.0_c_double)
            qs(i,j,k) = max(qs(i,j,k), 0.0_c_double)
            qg(i,j,k) = max(qg(i,j,k), 0.0_c_double)

            if (t(i,j,k) > t0c) then
              es = fpvs(t(i,j,k), 0, rd, rv, cpv, cliq, cice, xlv0, xls, psat, t0c)
            else
              es = fpvs(t(i,j,k), 1, rd, rv, cpv, cliq, cice, xlv0, xls, psat, t0c)
            end if
            es = min(es, 0.99_c_double * p(i,j,k))
            qsat = max(qmin, ep2 * es / max(p(i,j,k) - es, 1.0e-8_c_double))
            supsat = qv(i,j,k) - qsat

            cpm = cpd * (1.0_c_double - qv(i,j,k)) + cpv * qv(i,j,k)
            xlv = xlv0 - xlv1 * (t(i,j,k) - t0c)
            xlf = xlf0

            if (supsat > 0.0_c_double) then
              dq = min(supsat, 0.5_c_double * supsat)
              qv(i,j,k) = qv(i,j,k) - dq
              if (t(i,j,k) > t0c) then
                qc(i,j,k) = qc(i,j,k) + dq
                t(i,j,k) = t(i,j,k) + xlv * dq / max(cpm, 1.0e-6_c_double)
              else
                qi(i,j,k) = qi(i,j,k) + dq
                t(i,j,k) = t(i,j,k) + xls * dq / max(cpm, 1.0e-6_c_double)
              end if
            else if (supsat < 0.0_c_double) then
              dq = min(-supsat, 0.5_c_double * (-supsat))
              if (t(i,j,k) > t0c) then
                dq = min(dq, qc(i,j,k))
                qc(i,j,k) = qc(i,j,k) - dq
                qv(i,j,k) = qv(i,j,k) + dq
                t(i,j,k) = t(i,j,k) - xlv * dq / max(cpm, 1.0e-6_c_double)
              else
                dq = min(dq, qi(i,j,k))
                qi(i,j,k) = qi(i,j,k) - dq
                qv(i,j,k) = qv(i,j,k) + dq
                t(i,j,k) = t(i,j,k) - xls * dq / max(cpm, 1.0e-6_c_double)
              end if
            end if

            denfac = sqrt(max(den0 / max(den(i,j,k), 1.0e-12_c_double), 1.0e-12_c_double))
            supcol = t0c - t(i,j,k)
            n0sfac = max(min(exp(alpha * supcol), n0smax / n0s), 1.0_c_double)
            if (qr(i,j,k) <= qcrmin) then
              rslope_r = rslopermax
              rslopeb_r = rsloperbmax
              rslope2_r = rsloper2max
              rslope3_r = rsloper3max
            else
              lamdar = sqrt(sqrt(pidn0r / max(qr(i,j,k) * den(i,j,k), 1.0e-20_c_double)))
              rslope_r = 1.0_c_double / max(lamdar, 1.0e-12_c_double)
              rslopeb_r = rslope_r**bvtr
              rslope2_r = rslope_r * rslope_r
              rslope3_r = rslope2_r * rslope_r
            end if
            if (qs(i,j,k) <= qcrmin) then
              rslope_s = rslopesmax
              rslopeb_s = rslopesbmax
              rslope2_s = rslopes2max
              rslope3_s = rslopes3max
            else
              lamdas = sqrt(sqrt(pidn0s * n0sfac / max(qs(i,j,k) * den(i,j,k), 1.0e-20_c_double)))
              rslope_s = 1.0_c_double / max(lamdas, 1.0e-12_c_double)
              rslopeb_s = rslope_s**bvts
              rslope2_s = rslope_s * rslope_s
              rslope3_s = rslope2_s * rslope_s
            end if
            if (qg(i,j,k) <= qcrmin) then
              rslope_g = rslopegmax
              rslopeb_g = rslopegbmax
              rslope2_g = rslopeg2max
              rslope3_g = rslopeg3max
            else
              lamdag = sqrt(sqrt(pidn0g / max(qg(i,j,k) * den(i,j,k), 1.0e-20_c_double)))
              rslope_g = 1.0_c_double / max(lamdag, 1.0e-12_c_double)
              rslopeb_g = rslope_g**bvtg
              rslope2_g = rslope_g * rslope_g
              rslope3_g = rslope2_g * rslope_g
            end if

            supsat = max(qv(i,j,k), qmin) - qsat
            satdt = supsat / dtcld
            praut = 0.0_c_double
            pracw = 0.0_c_double
            prevp = 0.0_c_double

            if (qc(i,j,k) > qc0) then
              praut = qck1 * qc(i,j,k)**(7.0_c_double / 3.0_c_double)
              praut = min(praut, qc(i,j,k) / dtcld)
            end if

            if (qr(i,j,k) > qcrmin .and. qc(i,j,k) > qmin) then
              pracw = pacrr * rslope3_r * rslopeb_r * qc(i,j,k) * denfac
              pracw = min(pracw, qc(i,j,k) / dtcld)
            end if

            if (qr(i,j,k) > 0.0_c_double) then
              work1_r = diffac(xlv, p(i,j,k), t(i,j,k), den(i,j,k), qsat, rv)
              work2_r = venfac(p(i,j,k), t(i,j,k), den(i,j,k), den0)
              coeres = rslope2_r * sqrt(max(rslope_r * rslopeb_r, 0.0_c_double))
              prevp = (max(qv(i,j,k) / qsat, qmin) - 1.0_c_double) * &
                      (precr1 * rslope2_r + precr2 * work2_r * coeres) / max(work1_r, 1.0e-12_c_double)
              if (prevp < 0.0_c_double) then
                prevp = max(prevp, -qr(i,j,k) / dtcld)
                prevp = max(prevp, satdt / 2.0_c_double)
              else
                prevp = min(prevp, satdt / 2.0_c_double)
              end if
            end if

            qc(i,j,k) = max(qc(i,j,k) - (praut + pracw) * dtcld, 0.0_c_double)
            qr(i,j,k) = max(qr(i,j,k) + (praut + pracw + prevp) * dtcld, 0.0_c_double)
            qv(i,j,k) = max(qv(i,j,k) - prevp * dtcld, qmin)
            t(i,j,k) = t(i,j,k) + xlv * prevp * dtcld / max(cpm, 1.0e-6_c_double)

            psmlt = 0.0_c_double
            pgmlt = 0.0_c_double
            if (t(i,j,k) > t0c) then
              work2_r = venfac(p(i,j,k), t(i,j,k), den(i,j,k), den0)
              if (qs(i,j,k) > 0.0_c_double) then
                coeres = rslope2_s * sqrt(max(rslope_s * rslopeb_s, 0.0_c_double))
                psmlt = xka(t(i,j,k), den(i,j,k)) / max(xlf, 1.0e-12_c_double) * (t0c - t(i,j,k)) * pi / 2.0_c_double * &
                        n0sfac * (precs1 * rslope2_s + precs2 * work2_r * coeres) / max(den(i,j,k), 1.0e-12_c_double)
                psmlt = min(max(psmlt * dtcld, -qs(i,j,k)), 0.0_c_double)
                qs(i,j,k) = qs(i,j,k) + psmlt
                qr(i,j,k) = qr(i,j,k) - psmlt
                t(i,j,k) = t(i,j,k) + xlf / max(cpm, 1.0e-12_c_double) * psmlt
              end if
              if (qg(i,j,k) > 0.0_c_double) then
                coeres = rslope2_g * sqrt(max(rslope_g * rslopeb_g, 0.0_c_double))
                pgmlt = xka(t(i,j,k), den(i,j,k)) / max(xlf, 1.0e-12_c_double) * (t0c - t(i,j,k)) * &
                        (precg1 * rslope2_g + precg2 * work2_r * coeres) / max(den(i,j,k), 1.0e-12_c_double)
                pgmlt = min(max(pgmlt * dtcld, -qg(i,j,k)), 0.0_c_double)
                qg(i,j,k) = qg(i,j,k) + pgmlt
                qr(i,j,k) = qr(i,j,k) - pgmlt
                t(i,j,k) = t(i,j,k) + xlf / max(cpm, 1.0e-12_c_double) * pgmlt
              end if
            end if

            if (t(i,j,k) < t0c - 5.0_c_double) then
              freeze_q = min(qr(i,j,k), 0.05_c_double * qr(i,j,k))
              qr(i,j,k) = qr(i,j,k) - freeze_q
              qs(i,j,k) = qs(i,j,k) + freeze_q
            end if

            qi(i,j,k) = max(qi(i,j,k), 0.0_c_double)
            if (t(i,j,k) < t0c - 10.0_c_double .and. qi(i,j,k) > 0.0_c_double) then
              freeze_q = min(qi(i,j,k), 0.02_c_double * qi(i,j,k))
              qi(i,j,k) = qi(i,j,k) - freeze_q
              qs(i,j,k) = qs(i,j,k) + freeze_q
            end if
          end do
        end do
      end do

      do j = jts, jte
        do i = its, ite
          rain_flux = den(i,j,kts) * qr(i,j,kts) * max(delz(i,j,kts), 0.0_c_double)
          snow_flux = den(i,j,kts) * qs(i,j,kts) * max(delz(i,j,kts), 0.0_c_double)
          graup_flux = den(i,j,kts) * qg(i,j,kts) * max(delz(i,j,kts), 0.0_c_double)

          rainncv(i,j) = rainncv(i,j) + rain_flux
          snowncv(i,j) = snowncv(i,j) + snow_flux
          graupelncv(i,j) = graupelncv(i,j) + graup_flux

          rain(i,j) = rain(i,j) + rain_flux
          snow(i,j) = snow(i,j) + snow_flux
          graupel(i,j) = graupel(i,j) + graup_flux
          sr(i,j) = snow(i,j) / max(snow(i,j) + rain(i,j), 1.0e-12_c_double)
        end do
      end do
    end do

    errmsg = 'mp_wsm6_run OK (ERF phase-1 transitional)'
    errflg = 0_c_int
  end subroutine mp_wsm6_run

end module mp_wsm6

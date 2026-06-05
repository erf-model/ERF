module mp_wsm6_isohelper
  use iso_c_binding
  use mp_wsm6, only: mp_wsm6_init, mp_wsm6_run
  implicit none
  integer, parameter :: WSM6_DIAG_SCHEMA_V1 = 1

contains

  subroutine wsm6_diag_value_string(val, sout)
    real(c_double), intent(in) :: val
    character(len=*), intent(out) :: sout
    character(len=64) :: tmp

    write(tmp,'(SP,ES30.20E3)') val
    sout = trim(adjustl(tmp))
  end subroutine wsm6_diag_value_string



  subroutine mp_wsm6_init_c(den0, denr, dens, cl, cpv, hail_opt) bind(C, name="mp_wsm6_init_c")
    real(c_double), value, intent(in) :: den0, denr, dens, cl, cpv
    integer(c_int), value, intent(in) :: hail_opt
    character(len=256) :: errmsg
    integer :: errflg

    call mp_wsm6_init(den0, denr, dens, cl, cpv, int(hail_opt, kind(0)), errmsg, errflg)
    if (errflg /= 0) then
      write(*,'(A,1X,A)') 'mp_wsm6_init_c error:', trim(errmsg)
      stop 1
    end if
  end subroutine mp_wsm6_init_c

  subroutine mp_wsm6_run_c(t, qv, qc, qi, qr, qs, qg, den, p, delz, &
                           delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, xls, &
                           xlv0, xlf0, den0, denr, cliq, cice, psat, &
                           rain, rainncv, sr, snow, snowncv, graupel, graupelncv, &
                           ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte, microphysics_debug, diag_i_dbg, diag_j_dbg) bind(C, name="mp_wsm6_run_c")
    integer(c_int), value, intent(in) :: ims, ime, jms, jme, kms, kme
    integer(c_int), value, intent(in) :: its, ite, jts, jte, kts, kte
    integer(c_int), value, intent(in) :: microphysics_debug
    integer(c_int), value, intent(in) :: diag_i_dbg, diag_j_dbg
    real(c_double), intent(inout), dimension(ims:ime, jms:jme, kms:kme) :: t, qv, qc, qi, qr, qs, qg
    real(c_double), intent(in), dimension(ims:ime, jms:jme, kms:kme) :: den, p, delz
    real(c_double), value, intent(in) :: delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, xls
    real(c_double), value, intent(in) :: xlv0, xlf0, den0, denr, cliq, cice, psat
    real(c_double), intent(inout), dimension(ims:ime, jms:jme) :: rain, rainncv, sr
    real(c_double), intent(inout), dimension(ims:ime, jms:jme) :: snow, snowncv, graupel, graupelncv

    integer :: i, j, k, kk, kdim, debug_local, i_dbg_local, j_dbg_target
    logical :: i_dbg_in_tile
    integer :: errflg
    character(len=256) :: errmsg
    real(c_double) :: xlv1_bridge

    real(c_double), allocatable :: t_col(:,:), q_col(:,:), qc_col(:,:), qi_col(:,:), qr_col(:,:), qs_col(:,:), qg_col(:,:)
    real(c_double), allocatable :: den_col(:,:), p_col(:,:), delz_col(:,:)
    real(c_double), allocatable :: rain_col(:), rainncv_col(:), sr_col(:)
    real(c_double), allocatable :: snow_col(:), snowncv_col(:), graupel_col(:), graupelncv_col(:)

    if (its < ims .or. ite > ime .or. jts < jms .or. jte > jme .or. kts < kms .or. kte > kme) then
      write(*,'(A)') 'mp_wsm6_run_c bounds error: run-window outside storage bounds'
      write(*,'(A,6(1X,I0))') '  storage ims ime jms jme kms kme =', ims, ime, jms, jme, kms, kme
      write(*,'(A,6(1X,I0))') '  active  its ite jts jte kts kte =', its, ite, jts, jte, kts, kte
      stop 1
    end if
    if (its > ite .or. jts > jte .or. kts > kte) then
      write(*,'(A)') 'mp_wsm6_run_c bounds error: invalid active index ordering'
      write(*,'(A,6(1X,I0))') '  active its ite jts jte kts kte =', its, ite, jts, jte, kts, kte
      stop 1
    end if

    i_dbg_local = diag_i_dbg
    i_dbg_in_tile = (i_dbg_local >= its .and. i_dbg_local <= ite)
    if (.not. i_dbg_in_tile) i_dbg_local = its
    j_dbg_target = diag_j_dbg
    if (j_dbg_target < jts .or. j_dbg_target > jte) j_dbg_target = jts - 1

    kdim = kte - kts + 1
    xlv1_bridge = cliq - cpv

    ! Scratch columns passed to mp_wsm6_run use active bounds its:ite
    ! because the Fortran core dummy arrays are declared dimension(its:...).
    ! Raw C-binding arrays remain storage-bounded ims:ime.
    allocate(t_col(its:ite,1:kdim), q_col(its:ite,1:kdim), qc_col(its:ite,1:kdim), qi_col(its:ite,1:kdim))
    allocate(qr_col(its:ite,1:kdim), qs_col(its:ite,1:kdim), qg_col(its:ite,1:kdim))
    allocate(den_col(its:ite,1:kdim), p_col(its:ite,1:kdim), delz_col(its:ite,1:kdim))
    allocate(rain_col(its:ite), rainncv_col(its:ite), sr_col(its:ite))
    allocate(snow_col(its:ite), snowncv_col(its:ite), graupel_col(its:ite), graupelncv_col(its:ite))

    do j = jts, jte
      do k = kts, kte
        kk = k - kts + 1
        do i = its, ite
          t_col(i,kk)    = t(i,j,k)
          q_col(i,kk)    = qv(i,j,k)
          qc_col(i,kk)   = qc(i,j,k)
          qi_col(i,kk)   = qi(i,j,k)
          qr_col(i,kk)   = qr(i,j,k)
          qs_col(i,kk)   = qs(i,j,k)
          qg_col(i,kk)   = qg(i,j,k)
          den_col(i,kk)  = den(i,j,k)
          p_col(i,kk)    = p(i,j,k)
          delz_col(i,kk) = delz(i,j,k)
        end do
      end do

      do i = its, ite
        rain_col(i)      = rain(i,j)
        rainncv_col(i)   = rainncv(i,j)
        sr_col(i)        = sr(i,j)
        snow_col(i)      = snow(i,j)
        snowncv_col(i)   = snowncv(i,j)
        graupel_col(i)   = graupel(i,j)
        graupelncv_col(i)= graupelncv(i,j)
      end do

      debug_local = int(microphysics_debug, kind(0))
      if (microphysics_debug >= 1_c_int .and. (j /= j_dbg_target .or. .not. i_dbg_in_tile)) debug_local = 0

      call mp_wsm6_run(t_col, q_col, qc_col, qi_col, qr_col, qs_col, qg_col, den_col, p_col, delz_col, &
                       delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, xls, xlv0, xlf0, den0, denr, &
                       cliq, cice, psat, rain_col, rainncv_col, sr_col, snow_col, snowncv_col, &
                       graupel_col, graupelncv_col, its=its, ite=ite, kts=1, kte=kdim, &
                       microphysics_debug=debug_local, diag_i_dbg=i_dbg_local, diag_j_dbg=j, diag_k_raw_base=kts, &
                       errmsg=errmsg, errflg=errflg)

      if (errflg /= 0) then
        write(*,'(A,1X,I0,2A)') 'mp_wsm6_run_c error at j=', j, ': ', trim(errmsg)
        write(*,'(A,6(1X,I0))') '  storage ims ime jms jme kms kme =', ims, ime, jms, jme, kms, kme
        write(*,'(A,6(1X,I0))') '  active  its ite jts jte kts kte =', its, ite, jts, jte, kts, kte
        stop 1
      end if

      do k = kts, kte
        kk = k - kts + 1
        do i = its, ite
          t(i,j,k)  = t_col(i,kk)
          qv(i,j,k) = q_col(i,kk)
          qc(i,j,k) = qc_col(i,kk)
          qi(i,j,k) = qi_col(i,kk)
          qr(i,j,k) = qr_col(i,kk)
          qs(i,j,k) = qs_col(i,kk)
          qg(i,j,k) = qg_col(i,kk)
        end do
      end do

      do i = its, ite
        rain(i,j)       = rain_col(i)
        rainncv(i,j)    = rainncv_col(i)
        sr(i,j)         = sr_col(i)
        snow(i,j)       = snow_col(i)
        snowncv(i,j)    = snowncv_col(i)
        graupel(i,j)    = graupel_col(i)
        graupelncv(i,j) = graupelncv_col(i)
      end do
    end do

    deallocate(t_col, q_col, qc_col, qi_col, qr_col, qs_col, qg_col)
    deallocate(den_col, p_col, delz_col)
    deallocate(rain_col, rainncv_col, sr_col, snow_col, snowncv_col, graupel_col, graupelncv_col)
  end subroutine mp_wsm6_run_c

end module mp_wsm6_isohelper

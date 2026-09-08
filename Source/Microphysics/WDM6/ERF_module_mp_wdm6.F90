
   module mp_wdm6

   use iso_c_binding, only: c_double
   use module_libmassv, only: vrec_d, vsqrt_d
   use mp_radar
   use module_model_constants, only : RE_QC_BG, RE_QI_BG, RE_QS_BG

   implicit none
   integer, parameter :: kind_phys = c_double

   private
   public :: mp_wdm6_init, mp_wdm6_run










































   real(kind=kind_phys), parameter, private :: dtcldcr = 120.
   real(kind=kind_phys), parameter, private :: n0r     = 8.e6
   real(kind=kind_phys), parameter, private :: n0s     = 2.e6

   real(kind=kind_phys), parameter, private :: n0smax  =  1.e11
   real(kind=kind_phys), parameter, private :: dens    = 100.0
   real(kind=kind_phys), parameter, private :: alpha   = .12
   real(kind=kind_phys), parameter, private :: avtr    = 841.9
   real(kind=kind_phys), parameter, private :: bvtr    = 0.8
   real(kind=kind_phys), parameter, private :: avts    = 11.72
   real(kind=kind_phys), parameter, private :: bvts    = .41


   real(kind=kind_phys), parameter, private :: lamdacmax = 5.e5
   real(kind=kind_phys), parameter, private :: lamdacmin = 2.e4
   real(kind=kind_phys), parameter, private :: lamdarmax = 5.e4
   real(kind=kind_phys), parameter, private :: lamdarmin = 2.e3
   real(kind=kind_phys), parameter, private :: lamdasmax = 1.e5

   real(kind=kind_phys), parameter, private :: r0      = .8e-5
   real(kind=kind_phys), parameter, private :: peaut   = .55
   real(kind=kind_phys), parameter, private :: xncr    = 3.e8
   real(kind=kind_phys), parameter, private :: xncr0   = 5.e7
   real(kind=kind_phys), parameter, private :: xncr1   = 5.e8
   real(kind=kind_phys), parameter, private :: xmyu    = 1.718e-5
   real(kind=kind_phys), parameter, private :: dicon   = 11.9
   real(kind=kind_phys), parameter, private :: dimax   = 500.e-6
   real(kind=kind_phys), parameter, private :: pfrz1   = 100.
   real(kind=kind_phys), parameter, private :: pfrz2   = 0.66
   real(kind=kind_phys), parameter, private :: qcrmin  = 1.e-9
   real(kind=kind_phys), parameter, private :: ncmin   = 1.e1
   real(kind=kind_phys), parameter, private :: nrmin   = 1.e-2
   real(kind=kind_phys), parameter, private :: eacrc   = 1.0
   real(kind=kind_phys), parameter, private :: qs0     = 6.e-4
   real(kind=kind_phys), parameter, private :: satmax  = 1.0048
   real(kind=kind_phys), parameter, private :: actk    = 0.6
   real(kind=kind_phys), parameter, private :: actr    = 1.5
   real(kind=kind_phys), parameter, private :: ncrk1   = 3.03e3
   real(kind=kind_phys), parameter, private :: ncrk2   = 2.59e15
   real(kind=kind_phys), parameter, private :: di100   = 1.e-4
   real(kind=kind_phys), parameter, private :: di600   = 6.e-4
   real(kind=kind_phys), parameter, private :: di2000  = 2000.e-6
   real(kind=kind_phys), parameter, private :: di82    = 82.e-6
   real(kind=kind_phys), parameter, private :: di15    = 15.e-6
   real(kind=kind_phys), save ::                                                               &
             qc0,qc1,qck1,pidnc,bvtr1,bvtr2,bvtr3,bvtr4,bvtr5,                 &
             bvtr6,bvtr7, bvtr2o5,bvtr3o5,                                     &
             g1pbr,g2pbr,g3pbr,g4pbr,g5pbr,g6pbr,g7pbr,                        &
             g5pbro2,g7pbro2,pi,                                               &
             pvtr,pvtrn,eacrr,pacrr,pidn0r,pidnr,                              &
             precr1,precr2,xmmax,roqimax,bvts1,bvts2,                          &
             bvts3,bvts4,g1pbs,g3pbs,g4pbs,g5pbso2,                            &
             pvts,pacrs,precs1,precs2,pidn0s,xlv1,pacrc,                       &
             bvtg1,bvtg2,bvtg3,bvtg4,g1pbg,g3pbg,g4pbg,                        &
             g5pbgo2,pvtg,pacrg,precg1,precg2,pidn0g,                          &
             n0g,avtg,bvtg,deng,lamdagmax,                                     &
             rslopecmax,rslopec2max,rslopec3max,                               &
             rslopermax,rslopesmax,rslopegmax,                                 &
             rsloperbmax,rslopesbmax,rslopegbmax,                              &
             rsloper2max,rslopes2max,rslopeg2max,                              &
             rsloper3max,rslopes3max,rslopeg3max

    contains


    subroutine wdm6(th, q, qc, qr, qi, qs, qg,             &
                    nn, nc, nr,                            &
                    den, pii, p, delz,                     &
                    delt,g, cpd, cpv, ccn0, rd, rv, t0c,   &
                    ep1, ep2, qmin,                        &
                    xls, xlv0, xlf0, den0, denr,           &
                    cliq,cice,psat,                        &
                    xland,                                 &
                    rain, rainncv,                         &
                    snow, snowncv,                         &
                    sr,                                    &
                    refl_10cm, diagflag, do_radar_ref,     &
                    graupel, graupelncv,                   &
                    itimestep,                             &
                    has_reqc, has_reqi, has_reqs,          &
                    re_cloud, re_ice,   re_snow,           &
                    ids,ide, jds,jde, kds,kde,             &
                    ims,ime, jms,jme, kms,kme,             &
                    its,ite, jts,jte, kts,kte              &
                                                           )

   implicit none

   integer                                   , intent(in   ) :: itimestep
   integer                                   , intent(in   ) ::                &
                                                    ids,ide, jds,jde, kds,kde, &
                                                    ims,ime, jms,jme, kms,kme, &
                                                    its,ite, jts,jte, kts,kte
   real(kind=kind_phys)                      , intent(in   ) :: delt, g, ccn0, &
                                                                rd, rv, t0c,   &
                                                                cpd, cpv,      &
                                                                den0, qmin,    &
                                                                ep1, ep2, xls, &
                                                                xlv0, xlf0,    &
                                                                cliq, cice,    &
                                                                psat, denr
   real(kind=kind_phys), dimension(ims:ime,jms:jme)          , intent(in   ) :: xland
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(in   ) :: den
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(in   ) :: pii
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(in   ) :: p
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(in   ) :: delz
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(inout) :: th
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(inout) :: q
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(inout) :: qc
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(inout) :: qi
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(inout) :: qr
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(inout) :: qs
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(inout) :: qg
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(inout) :: nn
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(inout) :: nc
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme)  , intent(inout) :: nr
   real(kind=kind_phys), dimension(ims:ime,jms:jme)          , intent(inout) :: rain
   real(kind=kind_phys), dimension(ims:ime,jms:jme)          , intent(inout) :: rainncv
   real(kind=kind_phys), dimension(ims:ime,jms:jme)          , intent(inout) :: sr
   real(kind=kind_phys), dimension(ims:ime,jms:jme), optional, intent(inout) :: snow
   real(kind=kind_phys), dimension(ims:ime,jms:jme), optional, intent(inout) :: snowncv
   real(kind=kind_phys), dimension(ims:ime,jms:jme), optional, intent(inout) :: graupel
   real(kind=kind_phys), dimension(ims:ime,jms:jme), optional, intent(inout) :: graupelncv



   integer                                 , intent(in   ) :: has_reqc
   integer                                 , intent(in   ) :: has_reqi
   integer                                 , intent(in   ) :: has_reqs
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: re_cloud
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: re_ice
   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: re_snow



   real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: refl_10cm
   logical                       , optional, intent(in   ) :: diagflag
   integer                       , optional, intent(in   ) :: do_radar_ref



   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: t
   real(kind=kind_phys), dimension(its:ite,kts:kte,2) :: qci
   real(kind=kind_phys), dimension(its:ite,kts:kte,3) :: qrs, ncr
   integer                            :: i,j,k



   real(kind=kind_phys), dimension(kts:kte) :: qv1d
   real(kind=kind_phys), dimension(kts:kte) :: t1d
   real(kind=kind_phys), dimension(kts:kte) :: p1d
   real(kind=kind_phys), dimension(kts:kte) :: qr1d
   real(kind=kind_phys), dimension(kts:kte) :: nr1d
   real(kind=kind_phys), dimension(kts:kte) :: qs1d
   real(kind=kind_phys), dimension(kts:kte) :: qg1d
   real(kind=kind_phys), dimension(kts:kte) :: dBZ



   real(kind=kind_phys), dimension(kts:kte) :: qc1d, nc1d, den1d
   real(kind=kind_phys), dimension(kts:kte) :: qi1d
   real(kind=kind_phys), dimension(kts:kte) :: re_qc, re_qi, re_qs

   if (itimestep .eq. 1) then
     do j = jms,jme
       do k = kms,kme
         do i = ims,ime
           nn(i,k,j) = ccn0
         enddo
       enddo
     enddo
   endif

   do j = jts,jte
     do k = kts,kte
       do i = its,ite
         t(i,k) = th(i,k,j)*pii(i,k,j)
         qci(i,k,1) = qc(i,k,j)
         qci(i,k,2) = qi(i,k,j)
         qrs(i,k,1) = qr(i,k,j)
         qrs(i,k,2) = qs(i,k,j)
         qrs(i,k,3) = qg(i,k,j)
         ncr(i,k,1) = nn(i,k,j)
         ncr(i,k,2) = nc(i,k,j)
         ncr(i,k,3) = nr(i,k,j)
       enddo
     enddo




     call wdm62D(t, q(ims,kms,j), qci, qrs, ncr               &
                ,den(ims,kms,j)                               &
                ,p(ims,kms,j), delz(ims,kms,j)                &
                ,delt,g, cpd, cpv, ccn0, rd, rv, t0c          &
                ,ep1, ep2, qmin                               &
                ,xls, xlv0, xlf0, den0, denr                  &
                ,cliq,cice,psat                               &
                ,j                                            &
                ,xland(ims,j)                                 &
                ,rain(ims,j),rainncv(ims,j)                   &
                ,sr(ims,j)                                    &
                ,ids,ide, jds,jde, kds,kde                    &
                ,ims,ime, jms,jme, kms,kme                    &
                ,its,ite, jts,jte, kts,kte                    &
                ,snow(ims,j),snowncv(ims,j)                   &
                ,graupel(ims,j),graupelncv(ims,j)             &
                                                               )

     do k = kts,kte
       do i = its,ite
         th(i,k,j) = t(i,k)/pii(i,k,j)
         qc(i,k,j) = qci(i,k,1)
         qi(i,k,j) = qci(i,k,2)
         qr(i,k,j) = qrs(i,k,1)
         qs(i,k,j) = qrs(i,k,2)
         qg(i,k,j) = qrs(i,k,3)
         nn(i,k,j) = ncr(i,k,1)
         nc(i,k,j) = ncr(i,k,2)
         nr(i,k,j) = ncr(i,k,3)
       enddo
     enddo

     if ( present (diagflag) ) then
       if (diagflag .and. do_radar_ref == 1) then
         do i = its,ite
           do k = kts,kte
             t1d(k)  = th(i,k,j)*pii(i,k,j)
             p1d(k)  = p(i,k,j)
             qv1d(k) = q(i,k,j)
             qr1d(k) = qr(i,k,j)
             nr1d(k) = nr(i,k,j)
             qs1d(k) = qs(i,k,j)
             qg1d(k) = qg(i,k,j)
           enddo
           call refl10cm_wdm6 (qv1d, qr1d, nr1d, qs1d, qg1d,        &
                               t1d, p1d, dBZ, kts, kte, i, j)
           do k = kts, kte
             refl_10cm(i,k,j) = max(-35., dBZ(k))
           enddo
         enddo
       endif
     endif



     if (has_reqc.ne.0 .and. has_reqi.ne.0 .and. has_reqs.ne.0) then
       do i = its,ite
         do k = kts,kte
           re_qc(k) = RE_QC_BG
           re_qi(k) = RE_QI_BG
           re_qs(k) = RE_QS_BG

           t1d(k)  = th(i,k,j)*pii(i,k,j)
           den1d(k)= den(i,k,j)
           qc1d(k) = qc(i,k,j)
           qi1d(k) = qi(i,k,j)
           qs1d(k) = qs(i,k,j)
           nc1d(k) = nc(i,k,j)
         enddo

         call effectRad_wdm6(t1d, qc1d, nc1d, qi1d, qs1d, den1d,   &
                             qmin, t0c, re_qc, re_qi, re_qs,       &
                             kts, kte, i, j)

         do k = kts,kte
           re_cloud(i,k,j) = max(RE_QC_BG, min(re_qc(k),  50.e-6))
           re_ice(i,k,j)   = max(RE_QI_BG, min(re_qi(k), 125.e-6))
           re_snow(i,k,j)  = max(RE_QS_BG, min(re_qs(k), 999.e-6))
         enddo
       enddo
     endif
   enddo

   end subroutine wdm6



   subroutine wdm62D(t, q, qci, qrs, ncr, den, p, delz            &
                   ,delt,g, cpd, cpv, ccn0, rd, rv, t0c           &
                   ,ep1, ep2, qmin                                &
                   ,xls, xlv0, xlf0, den0, denr                   &
                   ,cliq,cice,psat                                &
                   ,lat                                           &
                   ,slmsk                                         &
                   ,rain,rainncv                                  &
                   ,sr                                            &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte                     &
                   ,snow,snowncv                                  &
                   ,graupel,graupelncv                            &
                   ,microphysics_debug, diag_i_dbg, diag_j_dbg    &
                                                                  )

    implicit none


























































































   integer                                , intent(in   ) ::                  &
                                                   ids,ide, jds,jde, kds,kde, &
                                                   ims,ime, jms,jme, kms,kme, &
                                                   its,ite, jts,jte, kts,kte
   integer                                , intent(in   ) :: lat
   integer                      , optional, intent(in   ) :: microphysics_debug, diag_i_dbg, diag_j_dbg
   real(kind=kind_phys) , intent(in   ) :: delt
   real(kind=kind_phys) , intent(in   ) :: g, cpd, cpv, t0c, &
                                                             den0, rd, rv,     &
                                                             ep1, ep2, qmin,   &
                                                             xls, xlv0, xlf0,  &
                                                             cliq, cice, psat, &
                                                             denr
   real(kind=kind_phys) , intent(in   ) :: ccn0
   real(kind=kind_phys), dimension(ims:ime)               , intent(in   ) :: slmsk
   real(kind=kind_phys), dimension(ims:ime,kms:kme)       , intent(in   ) :: p
   real(kind=kind_phys), dimension(ims:ime,kms:kme)       , intent(in   ) :: delz
   real(kind=kind_phys), dimension(ims:ime,kms:kme)       , intent(in   ) :: den
   real(kind=kind_phys), dimension(its:ite,kts:kte)       , intent(inout) :: t
   real(kind=kind_phys), dimension(its:ite,kts:kte,2)     , intent(inout) :: qci
   real(kind=kind_phys), dimension(its:ite,kts:kte,3)     , intent(inout) :: qrs
   real(kind=kind_phys), dimension(its:ite,kts:kte,3)     , intent(inout) :: ncr
   real(kind=kind_phys), dimension(ims:ime,kms:kme)       , intent(inout) :: q
   real(kind=kind_phys), dimension(ims:ime)               , intent(inout) :: rain
   real(kind=kind_phys), dimension(ims:ime)               , intent(inout) :: rainncv
   real(kind=kind_phys), dimension(ims:ime)               , intent(inout) :: sr
   real(kind=kind_phys), dimension(ims:ime), optional     , intent(inout) :: snow
   real(kind=kind_phys), dimension(ims:ime), optional     , intent(inout) :: snowncv
   real(kind=kind_phys), dimension(ims:ime), optional     , intent(inout) :: graupel
   real(kind=kind_phys), dimension(ims:ime), optional     , intent(inout) :: graupelncv

   ! Local variables for debug parameters with defaults
   integer :: debug_local, i_dbg_local, j_dbg_local

   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: dend
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: qcr
   real(kind=kind_phys), dimension(its:ite,kts:kte,3) :: rh
   real(kind=kind_phys), dimension(its:ite,kts:kte,3) :: qs
   real(kind=kind_phys), dimension(its:ite,kts:kte,3) :: rslope, rslope2, rslope3, rslopeb
   real(kind=kind_phys), dimension(its:ite,kts:kte,3) :: falk, fall
   real(kind=kind_phys), dimension(its:ite,kts:kte,3) :: work1
   real(kind=kind_phys), dimension(its:ite,kts:kte,3) :: qrs_tmp
   real(kind=kind_phys), dimension(its:ite,kts:kte,2) :: avedia
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: rslopec, rslopec2,rslopec3
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: workn, falln, falkn
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: worka, workr
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: den_tmp, delz_tmp, ncr_tmp
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: lamdr_tmp
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: lamdc_tmp
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: falkc, work1c, work2c, fallc
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: dqr,dnr
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: pcact, prevp, psdep, pgdep, praut,    &
                                         psaut, pgaut, pracw, psacw, pgacw,    &
                                         pgacr, pgacs, psaci, pgmlt, praci,    &
                                         piacr, pracs, psacr, pgaci, pseml,    &
                                         pgeml, paacw, pigen, pidep, pcond,    &
                                         psmlt, psevp, pgevp
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: nraut, nracw, ncevp, nccol, nrcol,    &
                                         nsacw, ngacw, niacr, nsacr, ngacr,    &
                                         naacw, nseml, ngeml, ncact
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: xl, cpm, work2, denfac, xni, n0sfac,  &
                                         qsum
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: denqrs1, denqrs2, denqrs3
   real(kind=kind_phys), dimension(its:ite,kts:kte)   :: denqr1, denncr3, denqci
   real(kind=kind_phys), dimension(its:ite)           :: delqrs1, delqrs2, delqrs3, delncr3
   real(kind=kind_phys), dimension(its:ite)           :: delqi
   real(kind=kind_phys), dimension(its:ite)           :: tstepsnow, tstepgraup
   real(kind=kind_phys)               :: gfac, sfac



   integer                            :: kte_in



   real(kind=kind_phys),    dimension(its:ite)        :: tvec1
   integer, dimension(its:ite)        :: mnstep, numndt
   integer, dimension(its:ite)        :: mstep, numdt
   logical, dimension(its:ite)        :: flgcld
   real(kind=kind_phys)               :: temp
   real(kind=kind_phys)               :: vt2ave
   real(kind=kind_phys)               :: holdc, holdci
   real(kind=kind_phys)               :: cpmcal, xlcal, lamdac, diffus,        &
                                         viscos, xka, venfac, conden, diffac,  &
                                         x, y, z, a, b, c, d, e,               &
                                         ndt, qdt, holdrr, holdrs, holdrg,     &
                                         supcol, supcolt, pvt, coeres, supsat, &
                                         dtcld, xmi, eacrs, satdt, qimax,      &
                                         diameter, xni0, roqi0,                &
                                         fallsum, fallsum_qsi, fallsum_qg,     &
                                         vt2i, vt2r, vt2s, vt2g, acrfac,       &
                                         egs, egi, coecol,                     &
                                         xlwork2, factor, source, value,       &
                                         nfrzdtr, nfrzdtc,                     &
                                         taucon, lencon, lenconcr,             &
                                         xlf, pfrzdtc, pfrzdtr, supice,        &
                                         alpha2, delta2, delta3

   integer                            :: i, j, k, mstepmax,                    &
                                         iprt, latd, lond, loop, loops, ifsat, &
                                         n, idim, kdim



   real(kind=kind_phys)               :: dldti, xb, xai, tr, xbi, xa, hvap,    &
                                         cvap, hsub, dldt, ttp




   cpmcal(x) = cpd*(1.-max(x,qmin))+max(x,qmin)*cpv
   xlcal(x) = xlv0-xlv1*(x-t0c)






   lamdac(x,y,z)= exp(log(((pidnc*z)/(x*y)))*((.33333333)))




   diffus(x,y) = 8.794e-5 * exp(log(x)*(1.81)) / y
   viscos(x,y) = 1.496e-6 * (x*sqrt(x)) /(x+120.)/y
   xka(x,y)    = 1.414e3*viscos(x,y)*y
   diffac(a,b,c,d,e) = d*a*a/(xka(c,d)*rv*c*c)+1./(e*diffus(c,b))
   venfac(a,b,c) = exp(log((viscos(b,c)/diffus(b,a)))*((.3333333)))            &
                   /sqrt(viscos(b,c))*sqrt(sqrt(den0/c))
   conden(a,b,c,d,e) = (max(b,qmin)-c)/(1.+d*d/(rv*e)*c/(a*a))

   idim = ite-its+1
   kdim = kte-kts+1

   kte_in = kte

   ! Initialize debug parameters with defaults
   if (present(microphysics_debug)) then
     debug_local = microphysics_debug
   else
     debug_local = 0
   endif
   if (present(diag_i_dbg)) then
     i_dbg_local = diag_i_dbg
   else
     i_dbg_local = its
   endif
   if (present(diag_j_dbg)) then
     j_dbg_local = diag_j_dbg
   else
     j_dbg_local = lat
   endif

   do k = kts, kte
     do i = its, ite
       qci(i,k,1) = max(qci(i,k,1),0.0)
       qrs(i,k,1) = max(qrs(i,k,1),0.0)
       qci(i,k,2) = max(qci(i,k,2),0.0)
       qrs(i,k,2) = max(qrs(i,k,2),0.0)
       qrs(i,k,3) = max(qrs(i,k,3),0.0)
       ncr(i,k,1) = min(max(ncr(i,k,1),1.e8),2.e10)
       ncr(i,k,2) = max(ncr(i,k,2),0.0)
       ncr(i,k,3) = max(ncr(i,k,3),0.0)
     enddo
   enddo

   do k = kts,kte
     do i = its,ite
       dend(i,k) = den(i,k)
     enddo
   enddo





   do k = kts,kte_in
     do i = its,ite
       cpm(i,k) = cpmcal(q(i,k))
       xl(i,k) = xlcal(t(i,k))
     enddo
   enddo

   qcr(:,:) = 0.0
   do i = its,ite
     if(slmsk(i).eq.2) then
       qcr(i,:) = qc0
     else
       qcr(i,:) = qc1
     endif
   enddo

   do k = kts,kte
     do i = its,ite
       delz_tmp(i,k) = delz(i,k)
       den_tmp(i,k) = den(i,k)
     enddo
   enddo



   do i = its,ite
     rainncv(i) = 0.
     if(present(snowncv) .and. present(snow)) snowncv(i) = 0.
     if(present (graupelncv) .and. present(graupel)) graupelncv(i) = 0.
     sr(i) = 0.



     tstepsnow(i) = 0.
     tstepgraup(i) = 0.
   enddo



   loops = max(nint(delt/dtcldcr),1)
   dtcld = delt/loops
   if(delt.le.dtcldcr) dtcld = delt

   do loop = 1,loops



     do i = its,ite
       mstep(i) = 1
       mnstep(i) = 1
       flgcld(i) = .true.
     enddo

     do k = kts, kte
       call vrec_d( tvec1(its), den(its,k), ite-its+1)
       do i = its, ite
         tvec1(i) = tvec1(i)*den0
       enddo
       call vsqrt_d( denfac(its,k), tvec1(its), ite-its+1)
     enddo

     ! WDM6-F90 TAG: PRE_G1B



     hsub = xls
     hvap = xlv0
     cvap = cpv
     ttp = t0c+0.01
     dldt = cvap-cliq
     xa = -dldt/rv
     xb = xa+hvap/(rv*ttp)
     dldti = cvap-cice
     xai = -dldti/rv
     xbi = xai+hsub/(rv*ttp)

     ! WDM6-F90 TAG: G1CV_SAT_INT0

     do k = kts,kte_in
       do i = its,ite
         tr = ttp/t(i,k)
         ! WDM6-F90 TAG: G1CV_SAT_INT1
         qs(i,k,1) = psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
         qs(i,k,1) = min(qs(i,k,1),0.99*p(i,k))
         ! WDM6-F90 TAG: G1CV_SAT_INT10
         qs(i,k,1) = ep2*qs(i,k,1)/(p(i,k)-qs(i,k,1))
         ! WDM6-F90 TAG: G1CV_SAT_INT11
         qs(i,k,1) = max(qs(i,k,1),qmin)
         ! WDM6-F90 TAG: G1CV_SAT_INT12
         rh(i,k,1) = max(q(i,k)/qs(i,k,1),qmin)
         ! WDM6-F90 TAG: G1CV_SAT_INT13
         tr=ttp/t(i,k)
         ! WDM6-F90 TAG: G1CV_SAT_INT20
         if(t(i,k).lt.ttp) then
           ! WDM6-F90 TAG: G1CV_SAT_INT21
           qs(i,k,2) = psat*exp(log(tr)*(xai))*exp(xbi*(1.-tr))
         else
           ! WDM6-F90 TAG: G1CV_SAT_INT21
           qs(i,k,2) = psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
         endif
         qs(i,k,2) = min(qs(i,k,2),0.99*p(i,k))
         ! WDM6-F90 TAG: G1CV_SAT_INT22
         qs(i,k,2) = ep2*qs(i,k,2)/(p(i,k)-qs(i,k,2))
         ! WDM6-F90 TAG: G1CV_SAT_INT23
         qs(i,k,2) = max(qs(i,k,2),qmin)
         ! WDM6-F90 TAG: G1CV_SAT_INT24
         rh(i,k,2) = max(q(i,k)/qs(i,k,2),qmin)
         ! WDM6-F90 TAG: G1CV_SAT_INT25
       enddo
     enddo

     ! WDM6-F90 TAG: PRE_G1C




      do k = kts, kte
        do i = its, ite
          prevp(i,k) = 0.
          psdep(i,k) = 0.
          pgdep(i,k) = 0.
          praut(i,k) = 0.
          psaut(i,k) = 0.
          pgaut(i,k) = 0.
          pracw(i,k) = 0.
          praci(i,k) = 0.
          piacr(i,k) = 0.
          psaci(i,k) = 0.
          psacw(i,k) = 0.
          pracs(i,k) = 0.
          psacr(i,k) = 0.
          pgacw(i,k) = 0.
          paacw(i,k) = 0.
          pgaci(i,k) = 0.
          pgacr(i,k) = 0.
          pgacs(i,k) = 0.
          pigen(i,k) = 0.
          pidep(i,k) = 0.
          pcond(i,k) = 0.
          psmlt(i,k) = 0.
          pgmlt(i,k) = 0.
          pseml(i,k) = 0.
          pgeml(i,k) = 0.
          psevp(i,k) = 0.
          pgevp(i,k) = 0.
          pcact(i,k) = 0.
          falk(i,k,1) = 0.
          falk(i,k,2) = 0.
          falk(i,k,3) = 0.
          fall(i,k,1) = 0.
          fall(i,k,2) = 0.
          fall(i,k,3) = 0.
          fallc(i,k) = 0.
          falkc(i,k) = 0.
          falln(i,k) =0.
          falkn(i,k) =0.
          xni(i,k) = 1.e3
          nsacw(i,k) = 0.
          ngacw(i,k) = 0.
          naacw(i,k) = 0.
          niacr(i,k) = 0.
          nsacr(i,k) = 0.
          ngacr(i,k) = 0.
          nseml(i,k) = 0.
          ngeml(i,k) = 0.
          nracw(i,k) = 0.
          nccol(i,k) = 0.
          nrcol(i,k) = 0.
          ncact(i,k) = 0.
          nraut(i,k) = 0.
          ncevp(i,k) = 0.
        enddo
      enddo

      ! WDM6-F90 TAG: PRE_G2

      do k = kts, kte
        do i = its, ite
          if(qci(i,k,1).le.qmin .or. ncr(i,k,2).le.ncmin ) then
            rslopec(i,k) = rslopecmax
            rslopec2(i,k) = rslopec2max
            rslopec3(i,k) = rslopec3max
          else
            rslopec(i,k) = 1./lamdac(qci(i,k,1),den(i,k),ncr(i,k,2))
            rslopec2(i,k) = rslopec(i,k)*rslopec(i,k)
            rslopec3(i,k) = rslopec2(i,k)*rslopec(i,k)
          endif



          temp = (den(i,k)*max(qci(i,k,2),qmin))
          temp = sqrt(sqrt(temp*temp*temp))
          xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
        enddo
      enddo

      ! WDM6-F90 TAG: PRE_G3




      ! WDM6-F90 TAG: PRE_G4

      do k = kts, kte
        do i = its, ite
          qrs_tmp(i,k,1) = qrs(i,k,1)
          qrs_tmp(i,k,2) = qrs(i,k,2)
          qrs_tmp(i,k,3) = qrs(i,k,3)
          ncr_tmp(i,k) = ncr(i,k,3)
        enddo
      enddo
      call slope_wdm6(qrs_tmp,ncr_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2, &
                     rslope3,work1,workn,its,ite,kts,kte)

      ! WDM6-F90 TAG: POST_G4

      ! WDM6-F90 TAG: PRE_G5A



      mstepmax = 1
      numdt = 1
      do k = kte, kts, -1
        do i = its, ite
          work1(i,k,1) = work1(i,k,1)/delz(i,k)
          workn(i,k) = workn(i,k)/delz(i,k)
          numdt(i) = max(nint(max(work1(i,k,1),workn(i,k))*dtcld+.5),1)
          if(numdt(i).ge.mstep(i)) mstep(i) = numdt(i)
        enddo
      enddo
      do i = its, ite
        if(mstepmax.le.mstep(i)) mstepmax = mstep(i)
      enddo

      ! WDM6-F90 TAG: POST_G5A

      ! WDM6-F90 TAG: PRE_G5B

      do n = 1, mstepmax
        k = kte
        do i = its, ite
         if(n.le.mstep(i)) then
           falk(i,k,1) = dend(i,k)*qrs(i,k,1)*work1(i,k,1)/mstep(i)
           falkn(i,k) = ncr(i,k,3)*workn(i,k)/mstep(i)
           fall(i,k,1) = fall(i,k,1)+falk(i,k,1)
           falln(i,k) = falln(i,k)+falkn(i,k)
           qrs(i,k,1) = max(qrs(i,k,1)-falk(i,k,1)*dtcld/dend(i,k),0.)
           ncr(i,k,3) = max(ncr(i,k,3)-falkn(i,k)*dtcld,0.)
         endif
       enddo

       do k = kte_in-1,kts,-1
         do i = its,ite
           if(n.le.mstep(i)) then
             falk(i,k,1) = dend(i,k)*qrs(i,k,1)*work1(i,k,1)/mstep(i)
             falkn(i,k) = ncr(i,k,3)*workn(i,k)/mstep(i)
             fall(i,k,1) = fall(i,k,1)+falk(i,k,1)
             falln(i,k) = falln(i,k)+falkn(i,k)
             dqr(i,k) = min(falk(i,k,1)*dtcld/dend(i,k),qrs(i,k,1))
             dqr(i,k+1) = min(falk(i,k+1,1)*delz(i,k+1)/delz(i,k)              &
                          *dtcld/dend(i,k),qrs(i,k+1,1))
             dnr(i,k) = min(falkn(i,k)*dtcld,ncr(i,k,3))
             dnr(i,k+1) = min(falkn(i,k+1)*delz(i,k+1)/delz(i,k)*dtcld,        &
                          ncr(i,k+1,3))
             qrs(i,k,1) = max(qrs(i,k,1)-dqr(i,k)+dqr(i,k+1),0.)
             ncr(i,k,3) = max(ncr(i,k,3)-dnr(i,k)+dnr(i,k+1),0.)
            endif
          enddo
        enddo
        do k = kts, kte
          do i = its, ite
            qrs_tmp(i,k,1) = qrs(i,k,1)
            ncr_tmp(i,k) = ncr(i,k,3)
          enddo
        enddo
        call slope_rain(qrs_tmp,ncr_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2, &
                     rslope3,work1,workn,its,ite,kts,kte)
        do k = kte, kts, -1
          do i = its, ite
            work1(i,k,1) = work1(i,k,1)/delz(i,k)
            workn(i,k) = workn(i,k)/delz(i,k)
          enddo
        enddo
      enddo

      ! WDM6-F90 TAG: POST_G5B


      ! WDM6-F90 TAG: PRE_G5C

      do k = kte, kts, -1
        do i = its, ite
          qsum(i,k) = max( (qrs(i,k,2)+qrs(i,k,3)), 1.E-15)
          if(qsum(i,k) .gt. 1.e-15 ) then
            worka(i,k) = (work1(i,k,2)*qrs(i,k,2) + work1(i,k,3)*qrs(i,k,3)) &
                      /qsum(i,k)
          else
            worka(i,k) = 0.
          endif
          denqrs2(i,k) = den(i,k)*qrs(i,k,2)
          denqrs3(i,k) = den(i,k)*qrs(i,k,3)
        enddo
      enddo
      call nislfv_rain_plm6(idim,kdim,den_tmp,denfac,t,delz_tmp,worka,         &
                           denqrs2,denqrs3,delqrs2,delqrs3,dtcld,1,1)
      do k = kts, kte
        do i = its, ite
          qrs(i,k,2) = max(denqrs2(i,k)/den(i,k),0.)
          qrs(i,k,3) = max(denqrs3(i,k)/den(i,k),0.)
          fall(i,k,2) = denqrs2(i,k)*worka(i,k)/delz(i,k)
          fall(i,k,3) = denqrs3(i,k)*worka(i,k)/delz(i,k)
        enddo
      enddo
      do i = its, ite
        fall(i,1,2) = delqrs2(i)/delz(i,1)/dtcld
        fall(i,1,3) = delqrs3(i)/delz(i,1)/dtcld
      enddo

      ! WDM6-F90 TAG: POST_G5C

      ! WDM6-F90 TAG: PRE_G6

      do k = kts, kte
        do i = its, ite
          qrs_tmp(i,k,1) = qrs(i,k,1)
          qrs_tmp(i,k,2) = qrs(i,k,2)
          qrs_tmp(i,k,3) = qrs(i,k,3)
          ncr_tmp(i,k) = ncr(i,k,3)
        enddo
      enddo
      call slope_wdm6(qrs_tmp,ncr_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2, &
                     rslope3,work1,workn,its,ite,kts,kte)

      ! WDM6-F90 TAG: POST_G6

      ! WDM6-F90 TAG: PRE_G7

      do k = kte, kts, -1
        do i = its, ite
          supcol = t0c-t(i,k)
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(t(i,k).gt.t0c) then




            xlf = xlf0
            work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
            if(qrs(i,k,2).gt.0.) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psmlt(i,k) = xka(t(i,k),den(i,k))/xlf*(t0c-t(i,k))*pi/2.       &
                         *n0sfac(i,k)*(precs1*rslope2(i,k,2)                 &
                         +precs2*work2(i,k)*coeres)/den(i,k)
              psmlt(i,k) = min(max(psmlt(i,k)*dtcld,-qrs(i,k,2)),0.)




              if(qrs(i,k,2).gt.qcrmin) then
                sfac = rslope(i,k,2)*n0s*n0sfac(i,k)/qrs(i,k,2)
                ncr(i,k,3) = ncr(i,k,3) - sfac*psmlt(i,k)
              endif

              qrs(i,k,2) = qrs(i,k,2) + psmlt(i,k)
              qrs(i,k,1) = qrs(i,k,1) - psmlt(i,k)
              t(i,k) = t(i,k) + xlf/cpm(i,k)*psmlt(i,k)
            endif




            if(qrs(i,k,3).gt.0.) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgmlt(i,k) = xka(t(i,k),den(i,k))/xlf*(t0c-t(i,k))*(precg1     &
                           *rslope2(i,k,3) + precg2*work2(i,k)*coeres)       &
                           /den(i,k)
              pgmlt(i,k) = min(max(pgmlt(i,k)*dtcld,-qrs(i,k,3)),0.)




              if(qrs(i,k,3).gt.qcrmin) then
                gfac = rslope(i,k,3)*n0g/qrs(i,k,3)
                ncr(i,k,3) = ncr(i,k,3) - gfac*pgmlt(i,k)
              endif

              qrs(i,k,3) = qrs(i,k,3) + pgmlt(i,k)
              qrs(i,k,1) = qrs(i,k,1) - pgmlt(i,k)
              t(i,k) = t(i,k) + xlf/cpm(i,k)*pgmlt(i,k)
            endif
          endif
        enddo
      enddo

      ! WDM6-F90 TAG: POST_G7

      ! WDM6-F90 TAG: PRE_G8

      do k = kte, kts, -1
        do i = its, ite
          if(qci(i,k,2).le.0.) then
            work1c(i,k) = 0.
          else
            xmi = den(i,k)*qci(i,k,2)/xni(i,k)
            diameter  = max(min(dicon * sqrt(xmi),dimax), 1.e-25)
            work1c(i,k) = 1.49e4*exp(log(diameter)*(1.31))
          endif

          ! WDM6-F90 TAG: T2_G8_vice
        enddo
      enddo



      do k = kte, kts, -1
        do i = its, ite
          denqci(i,k) = den(i,k)*qci(i,k,2)
        enddo
      enddo
      call nislfv_rain_plmr(idim,kdim,den_tmp,denfac,t,delz_tmp,work1c,denqci,denqci,  &
                           delqi,dtcld,1,0,0)
      do k = kts, kte
        do i = its, ite
          qci(i,k,2) = max(denqci(i,k)/den(i,k),0.)
        enddo
      enddo
      do i = its, ite
        fallc(i,1) = delqi(i)/delz(i,1)/dtcld
      enddo

      ! WDM6-F90 TAG: POST_G8

      ! G9 (PRECIP) is surface-scoped and is deliberately NOT on the all-k
      ! contract: the native path stores this group's inputs as box2d surface
      ! slabs, so a counterpart tag can only read klo.
      ! WDM6-F90 TAG: PRE_G9

      do i = its, ite
        fallsum = fall(i,kts,1)+fall(i,kts,2)+fall(i,kts,3)+fallc(i,kts)
        fallsum_qsi = fall(i,kts,2)+fallc(i,kts)
        fallsum_qg = fall(i,kts,3)
        if(fallsum.gt.0.) then
          rainncv(i) = fallsum*delz(i,kts)/denr*dtcld*1000. + rainncv(i)
          rain(i) = fallsum*delz(i,kts)/denr*dtcld*1000. + rain(i)
        endif
        If(fallsum_qsi.gt.0.) then
          tstepsnow(i) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000. + tstepsnow(i)
          IF( PRESENT (snowncv) .AND. PRESENT (snow)) THEN
            snowncv(i) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000. + snowncv(i)
            snow(i) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000. + snow(i)
          ENDIF
        ENDIF
        IF(fallsum_qg.gt.0.) then
          tstepgraup(i) = fallsum_qg*delz(i,kts)/denr*dtcld*1000.              &
                            + tstepgraup(i)
          IF( PRESENT (graupelncv) .AND. PRESENT (graupel)) THEN
            graupelncv(i) = fallsum_qg*delz(i,kts)/denr*dtcld*1000.            &
                            + graupelncv(i)
            graupel(i) = fallsum_qg*delz(i,kts)/denr*dtcld*1000. + graupel(i)
          ENDIF
        ENDIF
        IF ( PRESENT (snowncv)) THEN
          if(fallsum.gt.0.)sr(i)=(snowncv(i) + graupelncv(i))/(rainncv(i)+1.e-12)
        ELSE
          if(fallsum.gt.0.)sr(i)=(tstepsnow(i) + tstepgraup(i))/(rainncv(i)+1.e-12)
        ENDIF
      enddo

      ! Every field here is a 1D surface array with no k index, so an all-k loop
      ! would emit identical lines. Single-cell by nature; see the PRE_G9 note.
      ! WDM6-F90 TAG: POST_G9

      ! WDM6-F90 TAG: PRE_G10A



      do k = kts, kte
        do i = its, ite
          supcol = t0c-t(i,k)
          xlf = xls-xl(i,k)
          if(supcol.lt.0.) xlf = xlf0
          if(supcol.lt.0 .and. qci(i,k,2).gt.0.) then
            qci(i,k,1) = qci(i,k,1) + qci(i,k,2)




            if(qci(i,k,2).gt.qmin) then
              ncr(i,k,2) = ncr(i,k,2) + xni(i,k)
            endif
            t(i,k) = t(i,k) - xlf/cpm(i,k)*qci(i,k,2)
            qci(i,k,2) = 0.
          endif

      ! WDM6-F90 TAG: POST_G10A

      ! WDM6-F90 TAG: PRE_G10B

          if(supcol.gt.40. .and. qci(i,k,1).gt.0.) then
            qci(i,k,2) = qci(i,k,2) + qci(i,k,1)




            if(ncr(i,k,2).gt.0.) ncr(i,k,2) = 0.
            t(i,k) = t(i,k) + xlf/cpm(i,k)*qci(i,k,1)
            qci(i,k,1) = 0.
          endif

      ! WDM6-F90 TAG: POST_G10B

      ! WDM6-F90 TAG: PRE_G10C

          if(supcol.gt.0. .and. qci(i,k,1).gt.qmin) then
            supcolt=min(supcol,70.)
            pfrzdtc = min(pi*pi*pfrz1*(exp(pfrz2*supcolt)-1.)*denr/den(i,k)    &
                     *ncr(i,k,2)*rslopec3(i,k)*rslopec3(i,k)/18.*dtcld         &
                     ,qci(i,k,1))




            if(ncr(i,k,2).gt.ncmin) then
              nfrzdtc = min(pi*pfrz1*(exp(pfrz2*supcolt)-1.)*ncr(i,k,2)        &
                      *rslopec3(i,k)/6.*dtcld,ncr(i,k,2))
              ncr(i,k,2) = ncr(i,k,2) - nfrzdtc
            endif
            qci(i,k,2) = qci(i,k,2) + pfrzdtc
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtc
            qci(i,k,1) = qci(i,k,1)-pfrzdtc
          endif

      ! WDM6-F90 TAG: POST_G10C

      ! WDM6-F90 TAG: PRE_G10D

          if(supcol.gt.0. .and. qrs(i,k,1).gt.0.) then
            supcolt=min(supcol,70.)
            pfrzdtr = min(140.*(pi*pi)*pfrz1*ncr(i,k,3)*denr/den(i,k)          &
                  *(exp(pfrz2*supcolt)-1.)*rslope3(i,k,1)*rslope3(i,k,1)       &
                  *dtcld,qrs(i,k,1))

            ! WDM6-F90 TAG: T2_G10D_pgfrz




            if(ncr(i,k,3).gt.nrmin) then
              nfrzdtr = min(4.*pi*pfrz1*ncr(i,k,3)*(exp(pfrz2*supcolt)-1.)     &
                       *rslope3(i,k,1)*dtcld, ncr(i,k,3))
              ncr(i,k,3) = ncr(i,k,3) - nfrzdtr
            endif
            qrs(i,k,3) = qrs(i,k,3) + pfrzdtr
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtr
            qrs(i,k,1) = qrs(i,k,1) - pfrzdtr
          endif
        enddo
      enddo

      ! WDM6-F90 TAG: POST_G10D

      ! WDM6-F90 TAG: PRE_G10E

      do k = kts, kte
        do i = its, ite
          ncr(i,k,2) = max(ncr(i,k,2),0.0)
          ncr(i,k,3) = max(ncr(i,k,3),0.0)
        enddo
      enddo

      ! WDM6-F90 TAG: POST_G10E

      ! WDM6-F90 TAG: PRE_G11

      do k = kts, kte
        do i = its, ite
          qrs_tmp(i,k,1) = qrs(i,k,1)
          qrs_tmp(i,k,2) = qrs(i,k,2)
          qrs_tmp(i,k,3) = qrs(i,k,3)
          ncr_tmp(i,k) = ncr(i,k,3)
        enddo
      enddo
      call slope_wdm6(qrs_tmp,ncr_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2, &
                     rslope3,work1,workn,its,ite,kts,kte)
      do k = kts, kte
        do i = its, ite




          avedia(i,k,2) = rslope(i,k,1)*((24.)**(.3333333))

          if(qci(i,k,1).le.qmin .or. ncr(i,k,2).le.ncmin) then
            rslopec(i,k) = rslopecmax
            rslopec2(i,k) = rslopec2max
            rslopec3(i,k) = rslopec3max
          else
            rslopec(i,k) = 1./lamdac(qci(i,k,1),den(i,k),ncr(i,k,2))
            rslopec2(i,k) = rslopec(i,k)*rslopec(i,k)
            rslopec3(i,k) = rslopec2(i,k)*rslopec(i,k)
          endif




          avedia(i,k,1) = rslopec(i,k)
        enddo
      enddo

      ! WDM6-F90 TAG: PRE_G12

      do k = kts, kte
        do i = its, ite
          ! WDM6-F90 TAG: G11V_DIFFAC_PRE1
          work1(i,k,1) = diffac(xl(i,k),p(i,k),t(i,k),den(i,k),qs(i,k,1))
          ! WDM6-F90 TAG: G11V_DIFFAC_POST1
          ! WDM6-F90 TAG: G11V_DIFFAC_PRE2
          work1(i,k,2) = diffac(xls,p(i,k),t(i,k),den(i,k),qs(i,k,2))
          ! WDM6-F90 TAG: G11V_DIFFAC_POST2
          ! WDM6-F90 TAG: G11V_VENFAC_PRE
          work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
          ! WDM6-F90 TAG: G11V_VENFAC_POST
        enddo
      enddo

      ! WDM6-F90 TAG: POST_G12

      ! WDM6-F90 TAG: POST_G11






      ! WDM6-F90 TAG: PRE_G13A

      do k = kts, kte
        do i = its, ite
          supsat = max(q(i,k),qmin)-qs(i,k,1)
          satdt = supsat/dtcld




         lencon  = 2.7e-2*den(i,k)*qci(i,k,1)*(1.e20/16.*rslopec2(i,k)         &
                  *rslopec2(i,k)-0.4)
         lenconcr = max(1.2*lencon, qcrmin)
         if(qci(i,k,1).gt.qcr(i,k).and.ncr(i,k,2).gt.ncmin) then
           praut(i,k) = qck1*qci(i,k,1)**(7./3.)*ncr(i,k,2)**(-1./3.)
           praut(i,k) = min(praut(i,k),qci(i,k,1)/dtcld)




           nraut(i,k) = 3.5e9*den(i,k)*praut(i,k)
           if(qrs(i,k,1).gt.lenconcr)                                           &
           nraut(i,k) = ncr(i,k,3)/qrs(i,k,1)*praut(i,k)
           nraut(i,k) = min(nraut(i,k),ncr(i,k,2)/dtcld)
         endif






          if(qrs(i,k,1).ge.lenconcr) then
            if(avedia(i,k,2).ge.di100) then
              nracw(i,k) = min(ncrk1*ncr(i,k,2)*ncr(i,k,3)*(rslopec3(i,k)      &
                         + 24.*rslope3(i,k,1)),ncr(i,k,2)/dtcld)
              pracw(i,k) = min(pi/6.*(denr/den(i,k))*ncrk1*ncr(i,k,2)          &
                         *ncr(i,k,3)*rslopec3(i,k)*(2.*rslopec3(i,k)           &
                         + 24.*rslope3(i,k,1)),qci(i,k,1)/dtcld)
            else
              nracw(i,k) = min(ncrk2*ncr(i,k,2)*ncr(i,k,3)*(2.*rslopec3(i,k)   &
                         *rslopec3(i,k)+5040.*rslope3(i,k,1)                   &
                         *rslope3(i,k,1)),ncr(i,k,2)/dtcld)
              pracw(i,k) = min(pi/6.*(denr/den(i,k))*ncrk2*ncr(i,k,2)          &
                         *ncr(i,k,3)*rslopec3(i,k)*(6.*rslopec3(i,k)           &
                         *rslopec3(i,k)+5040.*rslope3(i,k,1)*rslope3(i,k,1))   &
                         ,qci(i,k,1)/dtcld)
            endif
          endif




          if(avedia(i,k,1).ge.di100) then
            nccol(i,k) = ncrk1*ncr(i,k,2)*ncr(i,k,2)*rslopec3(i,k)
          else
            nccol(i,k) = 2.*ncrk2*ncr(i,k,2)*ncr(i,k,2)*rslopec3(i,k)        &
                         *rslopec3(i,k)
          endif




          if(qrs(i,k,1).ge.lenconcr) then
            if(avedia(i,k,2).lt.di100) then
              nrcol(i,k) = 5040.*ncrk2*ncr(i,k,3)*ncr(i,k,3)*rslope3(i,k,1)    &
                          *rslope3(i,k,1)
            elseif(avedia(i,k,2).ge.di100 .and. avedia(i,k,2).lt.di600) then
              nrcol(i,k) = 24.*ncrk1*ncr(i,k,3)*ncr(i,k,3)*rslope3(i,k,1)
            elseif(avedia(i,k,2).ge.di600 .and. avedia(i,k,2).lt.di2000) then
              coecol = -2.5e3*(avedia(i,k,2)-di600)
              nrcol(i,k) = 24.*exp(coecol)*ncrk1*ncr(i,k,3)*ncr(i,k,3)         &
                         *rslope3(i,k,1)
            else
              nrcol(i,k) = 0.
            endif
          endif




          if(qrs(i,k,1).gt.0.) then
            coeres = rslope(i,k,1)*sqrt(rslope(i,k,1)*rslopeb(i,k,1))
            prevp(i,k) = (rh(i,k,1)-1.)*ncr(i,k,3)*(precr1*rslope(i,k,1)       &
                         + precr2*work2(i,k)*coeres)/work1(i,k,1)
            if(prevp(i,k).lt.0.) then
              prevp(i,k) = max(prevp(i,k),-qrs(i,k,1)/dtcld)
              prevp(i,k) = max(prevp(i,k),satdt/2)




              if(prevp(i,k).eq.-qrs(i,k,1)/dtcld) then
                ncr(i,k,1) = ncr(i,k,1)+ncr(i,k,3)
                ncr(i,k,3) = 0.
              endif
            else if(prevp(i,k).eq.0.) then

              ! A zero kinetic rate must not become evaporation
              ! through a negative saturation-rate limiter.
              prevp(i,k) = 0.
            else

              prevp(i,k) = min(prevp(i,k),satdt/2)
            endif
          endif
        enddo
      enddo

      ! WDM6-F90 TAG: POST_G13A







      ! WDM6-F90 TAG: PRE_G13B

      do k = kts, kte
        do i = its, ite
          supcol = t0c-t(i,k)
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          supsat = max(q(i,k),qmin)-qs(i,k,2)
          satdt = supsat/dtcld
          ifsat = 0





          temp = (den(i,k)*max(qci(i,k,2),qmin))
          temp = sqrt(sqrt(temp*temp*temp))
          xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
          eacrs = exp(0.07*(-supcol))

          xmi = den(i,k)*qci(i,k,2)/xni(i,k)
          diameter  = min(dicon * sqrt(xmi),dimax)
          vt2i = 1.49e4*diameter**1.31
          vt2r=pvtr*rslopeb(i,k,1)*denfac(i,k)
          vt2s=pvts*rslopeb(i,k,2)*denfac(i,k)
          vt2g=pvtg*rslopeb(i,k,3)*denfac(i,k)
          qsum(i,k) = max((qrs(i,k,2)+qrs(i,k,3)),1.e-15)
          if(qsum(i,k) .gt. 1.e-15) then
            vt2ave=(vt2s*qrs(i,k,2)+vt2g*qrs(i,k,3))/(qsum(i,k))
          else
            vt2ave=0.
          endif
          if(supcol.gt.0. .and. qci(i,k,2).gt.qmin) then
            if(qrs(i,k,1).gt.qcrmin) then




              acrfac = 6.*rslope2(i,k,1)+4.*diameter*rslope(i,k,1) + diameter**2
              praci(i,k) = pi*qci(i,k,2)*ncr(i,k,3)*abs(vt2r-vt2i)*acrfac/4.

              praci(i,k) = praci(i,k)*min(max(0.0,qrs(i,k,1)/qci(i,k,2)),1.)**2
              praci(i,k) = min(praci(i,k),qci(i,k,2)/dtcld)




              piacr(i,k) = pi*pi*avtr*ncr(i,k,3)*denr*xni(i,k)*denfac(i,k)     &
                          *g7pbr*rslope3(i,k,1)*rslope2(i,k,1)*rslopeb(i,k,1)  &
                          /24./den(i,k)

              piacr(i,k) = piacr(i,k)*min(max(0.0,qci(i,k,2)/qrs(i,k,1)),1.)**2
              piacr(i,k) = min(piacr(i,k),qrs(i,k,1)/dtcld)
            endif




            if(ncr(i,k,3).gt.nrmin) then
              niacr(i,k) = pi*avtr*ncr(i,k,3)*xni(i,k)*denfac(i,k)*g4pbr       &
                          *rslope2(i,k,1)*rslopeb(i,k,1)/4.

              niacr(i,k) = niacr(i,k)*min(max(0.0,qci(i,k,2)/qrs(i,k,1)),1.)**2
              niacr(i,k) = min(niacr(i,k),ncr(i,k,3)/dtcld)
            endif




            if(qrs(i,k,2).gt.qcrmin) then
              acrfac = 2.*rslope3(i,k,2)+2.*diameter*rslope2(i,k,2)            &
                      + diameter**2*rslope(i,k,2)
              psaci(i,k) = pi*qci(i,k,2)*eacrs*n0s*n0sfac(i,k)                 &
                          *abs(vt2ave-vt2i)*acrfac/4.
              psaci(i,k) = min(psaci(i,k),qci(i,k,2)/dtcld)

            ! WDM6-F90 TAG: T2_G13B_psaci
            endif




            if(qrs(i,k,3).gt.qcrmin) then
              egi = exp(0.07*(-supcol))
              acrfac = 2.*rslope3(i,k,3)+2.*diameter*rslope2(i,k,3)            &
                      + diameter**2*rslope(i,k,3)
              pgaci(i,k) = pi*egi*qci(i,k,2)*n0g*abs(vt2ave-vt2i)*acrfac/4.
              pgaci(i,k) = min(pgaci(i,k),qci(i,k,2)/dtcld)
            endif
          endif




      ! WDM6-F90 TAG: POST_G13B

      ! WDM6-F90 TAG: PRE_G13C

          if(qrs(i,k,2).gt.qcrmin .and. qci(i,k,1).gt.qmin) then
            psacw(i,k) = min(pacrc*n0sfac(i,k)*rslope3(i,k,2)*rslopeb(i,k,2)   &

                        *min(max(0.0,qrs(i,k,2)/qci(i,k,1)),1.)**2             &
                        *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif




         if(qrs(i,k,2).gt.qcrmin .and. ncr(i,k,2).gt.ncmin) then
           nsacw(i,k) = min(pacrc*n0sfac(i,k)*rslope3(i,k,2)*rslopeb(i,k,2)    &

                       *min(max(0.0,qrs(i,k,2)/qci(i,k,1)),1.)**2              &
                       *ncr(i,k,2)*denfac(i,k),ncr(i,k,2)/dtcld)
         endif




          if(qrs(i,k,3).gt.qcrmin .and. qci(i,k,1).gt.qmin) then
            pgacw(i,k) = min(pacrg*rslope3(i,k,3)*rslopeb(i,k,3)*qci(i,k,1)    &

                        *min(max(0.0,qrs(i,k,3)/qci(i,k,1)),1.)**2             &
                        *denfac(i,k),qci(i,k,1)/dtcld)
          endif




          if(qrs(i,k,3).gt.qcrmin .and. ncr(i,k,2).gt.ncmin) then
            ngacw(i,k) = min(pacrg*rslope3(i,k,3)*rslopeb(i,k,3)*ncr(i,k,2)    &

                        *min(max(0.0,qrs(i,k,3)/qci(i,k,1)),1.)**2             &
                        *denfac(i,k),ncr(i,k,2)/dtcld)
          endif




          if(qsum(i,k) .gt. 1.e-15 ) then
            paacw(i,k) = (qrs(i,k,2)*psacw(i,k)+qrs(i,k,3)*pgacw(i,k))/(qsum(i,k))




            naacw(i,k) = (qrs(i,k,2)*nsacw(i,k)+qrs(i,k,3)*ngacw(i,k))/(qsum(i,k))
          endif

      ! WDM6-F90 TAG: POST_G13C

      ! WDM6-F90 TAG: PRE_G13D

          if(qrs(i,k,2).gt.qcrmin .and. qrs(i,k,1).gt.qcrmin) then
            if(supcol.gt.0) then
              acrfac = 5.*rslope3(i,k,2)*rslope3(i,k,2)                        &
                      + 4.*rslope3(i,k,2)*rslope2(i,k,2)*rslope(i,k,1)         &
                      + 1.5*rslope2(i,k,2)*rslope2(i,k,2)*rslope2(i,k,1)
              pracs(i,k) = pi*pi*ncr(i,k,3)*n0s*n0sfac(i,k)*abs(vt2r-vt2ave)   &
                          *(dens/den(i,k))*acrfac

              pracs(i,k) = pracs(i,k)*min(max(0.0,qrs(i,k,1)/qrs(i,k,2)),1.)**2
              pracs(i,k) = min(pracs(i,k),qrs(i,k,2)/dtcld)
            endif




            acrfac = 30.*rslope3(i,k,1)*rslope2(i,k,1)*rslope(i,k,2)           &
                     +10.*rslope2(i,k,1)*rslope2(i,k,1)*rslope2(i,k,2)         &
                     + 2.*rslope3(i,k,1)*rslope3(i,k,2)
            psacr(i,k) = pi*pi*ncr(i,k,3)*n0s*n0sfac(i,k)*abs(vt2ave-vt2r)     &
                        *(denr/den(i,k))*acrfac

            psacr(i,k) = psacr(i,k)*min(max(0.0,qrs(i,k,2)/qrs(i,k,1)),1.)**2
            psacr(i,k) = min(psacr(i,k),qrs(i,k,1)/dtcld)
          endif
          if(qrs(i,k,2).gt.qcrmin .and. ncr(i,k,3).gt.nrmin) then




            acrfac = 1.5*rslope2(i,k,1)*rslope(i,k,2)                          &
                    + 1.0*rslope(i,k,1)*rslope2(i,k,2)+.5*rslope3(i,k,2)
            nsacr(i,k) = pi*ncr(i,k,3)*n0s*n0sfac(i,k)*abs(vt2ave-vt2r)        &
                        *acrfac

            nsacr(i,k) = nsacr(i,k)*min(max(0.0,qrs(i,k,2)/qrs(i,k,1)),1.)**2
            nsacr(i,k) = min(nsacr(i,k),ncr(i,k,3)/dtcld)
          endif




          if(qrs(i,k,3).gt.qcrmin .and. qrs(i,k,1).gt.qcrmin) then
            acrfac = 30.*rslope3(i,k,1)*rslope2(i,k,1)*rslope(i,k,3)           &
                    +10.*rslope2(i,k,1)*rslope2(i,k,1)*rslope2(i,k,3)          &
                    + 2.*rslope3(i,k,1)*rslope3(i,k,3)
            pgacr(i,k) = pi*pi*ncr(i,k,3)*n0g*abs(vt2ave-vt2r)*(denr/den(i,k)) &
                        *acrfac

            pgacr(i,k) = pgacr(i,k)*min(max(0.0,qrs(i,k,3)/qrs(i,k,1)),1.)**2
            pgacr(i,k) = min(pgacr(i,k),qrs(i,k,1)/dtcld)
          endif




          if(qrs(i,k,3).gt.qcrmin .and. ncr(i,k,3).gt.nrmin) then
            acrfac = 1.5*rslope2(i,k,1)*rslope(i,k,3)                          &
                    + 1.0*rslope(i,k,1)*rslope2(i,k,3) + .5*rslope3(i,k,3)
            ngacr(i,k) = pi*ncr(i,k,3)*n0g*abs(vt2ave-vt2r)*acrfac

            ngacr(i,k) = ngacr(i,k)*min(max(0.0,qrs(i,k,3)/qrs(i,k,1)),1.)**2
            ngacr(i,k) = min(ngacr(i,k),ncr(i,k,3)/dtcld)
          endif






          if(qrs(i,k,3).gt.qcrmin .and. qrs(i,k,2).gt.qcrmin) then
            pgacs(i,k) = 0.
          endif
          if(supcol.le.0) then
            xlf = xlf0




            if(qrs(i,k,2).gt.0.)                                               &
              pseml(i,k) = min(max(cliq*supcol*(paacw(i,k)+psacr(i,k))         &
                          /xlf,-qrs(i,k,2)/dtcld),0.)




              if  (qrs(i,k,2).gt.qcrmin) then
                sfac = rslope(i,k,2)*n0s*n0sfac(i,k)/qrs(i,k,2)
                nseml(i,k) = -sfac*pseml(i,k)
              endif




            if(qrs(i,k,3).gt.0.)                                               &
              pgeml(i,k) = min(max(cliq*supcol*(paacw(i,k)+pgacr(i,k))/xlf     &
                          ,-qrs(i,k,3)/dtcld),0.)




              if (qrs(i,k,3).gt.qcrmin) then
                gfac = rslope(i,k,3)*n0g/qrs(i,k,3)
                ngeml(i,k) = -gfac*pgeml(i,k)
              endif
          endif
      ! WDM6-F90 TAG: POST_G13D
      ! WDM6-F90 TAG: PRE_G13E
          if(supcol.gt.0) then




            if(qci(i,k,2).gt.0. .and. ifsat.ne.1) then
              pidep(i,k) = 4.*diameter*xni(i,k)*(rh(i,k,2)-1.)/work1(i,k,2)
              supice = satdt-prevp(i,k)
              if(pidep(i,k).lt.0.) then
                pidep(i,k) = max(max(pidep(i,k),satdt/2),supice)
                pidep(i,k) = max(pidep(i,k),-qci(i,k,2)/dtcld)
              else
                pidep(i,k) = min(min(pidep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)).ge.abs(satdt)) ifsat = 1
            endif
            ! WDM6-F90 TAG: T2_G13E_pidep




            if(qrs(i,k,2).gt.0. .and. ifsat.ne.1) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psdep(i,k) = (rh(i,k,2)-1.)*n0sfac(i,k)*(precs1*rslope2(i,k,2)   &
                           + precs2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)
              if(psdep(i,k).lt.0.) then
                psdep(i,k) = max(psdep(i,k),-qrs(i,k,2)/dtcld)
                psdep(i,k) = max(max(psdep(i,k),satdt/2),supice)
              else
                psdep(i,k) = min(min(psdep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)+psdep(i,k)).ge.abs(satdt)) ifsat = 1
              ! WDM6-F90 TAG: T2_G13E_psdep
            endif




            if(qrs(i,k,3).gt.0. .and. ifsat.ne.1) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgdep(i,k) = (rh(i,k,2)-1.)*(precg1*rslope2(i,k,3)               &
                          + precg2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)
              if(pgdep(i,k).lt.0.) then
                pgdep(i,k) = max(pgdep(i,k),-qrs(i,k,3)/dtcld)
                pgdep(i,k) = max(max(pgdep(i,k),satdt/2),supice)
              else
                pgdep(i,k) = min(min(pgdep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)+psdep(i,k)+pgdep(i,k)).ge.          &
                abs(satdt)) ifsat = 1
            endif




            if(supsat.gt.0. .and. ifsat.ne.1) then
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)-pgdep(i,k)
              xni0 = 1.e3*exp(0.1*supcol)
              roqi0 = 4.92e-11*xni0**1.33
              ! WDM6-F90 TAG: T2_G13E_PIGEN_EXACT
              pigen(i,k) = max(0.,(roqi0/den(i,k)-max(qci(i,k,2),0.))/dtcld)
              ! WDM6-F90 TAG: T2_G13E_PIGEN_POST
              pigen(i,k) = min(min(pigen(i,k),satdt),supice)
            endif
      ! WDM6-F90 TAG: POST_G13E





      ! WDM6-F90 TAG: PRE_G13F
            if(qci(i,k,2).gt.0.) then
              qimax = roqimax/den(i,k)
              psaut(i,k) = max(0.,(qci(i,k,2)-qimax)/dtcld)
            endif





            if(qrs(i,k,2).gt.0.) then
              alpha2 = 1.e-3*exp(0.09*(-supcol))
              pgaut(i,k) = min(max(0.,alpha2*(qrs(i,k,2)-qs0)),qrs(i,k,2)/dtcld)
            endif
      ! WDM6-F90 TAG: POST_G13F
          endif





      ! WDM6-F90 TAG: PRE_G13G
          if(supcol.lt.0.) then
            if(qrs(i,k,2).gt.0. .and. rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psevp(i,k) = (rh(i,k,1)-1.)*n0sfac(i,k)*(precs1*rslope2(i,k,2)   &
                           +precs2*work2(i,k)*coeres)/work1(i,k,1)
              psevp(i,k) = min(max(psevp(i,k),-qrs(i,k,2)/dtcld),0.)
            endif




            if(qrs(i,k,3).gt.0. .and. rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgevp(i,k) = (rh(i,k,1)-1.)*(precg1*rslope2(i,k,3)               &
                         + precg2*work2(i,k)*coeres)/work1(i,k,1)
              pgevp(i,k) = min(max(pgevp(i,k),-qrs(i,k,3)/dtcld),0.)
            endif
          endif
      ! WDM6-F90 TAG: POST_G13G
        enddo
      enddo






      ! WDM6-F90 TAG: PRE_G14
      do k = kts, kte
        do i = its, ite

          delta2=0.
          delta3=0.
          if(qrs(i,k,1).lt.1.e-4 .and. qrs(i,k,2).lt.1.e-4) delta2=1.
          if(qrs(i,k,1).lt.1.e-4) delta3=1.
          if(t(i,k).le.t0c) then



            value = max(qmin,qci(i,k,1))
            source = (praut(i,k)+pracw(i,k)+paacw(i,k)+paacw(i,k))             &
                    *dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
            endif



            value = max(qmin,qci(i,k,2))
            source = (psaut(i,k)-pigen(i,k)-pidep(i,k)+praci(i,k)+psaci(i,k)   &
                    +pgaci(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              psaut(i,k) = psaut(i,k)*factor
              pigen(i,k) = pigen(i,k)*factor
              pidep(i,k) = pidep(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
            endif



            value = max(qmin,qrs(i,k,1))
            source = (-praut(i,k)-prevp(i,k)-pracw(i,k)+piacr(i,k)             &
                    +psacr(i,k)+pgacr(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              prevp(i,k) = prevp(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
            endif



            value = max(qmin,qrs(i,k,2))
            source = -(psdep(i,k)+psaut(i,k)-pgaut(i,k)+paacw(i,k)             &
                     +piacr(i,k)*delta3+praci(i,k)*delta3                      &
                     -pracs(i,k)*(1.-delta2)+psacr(i,k)*delta2                 &
                     +psaci(i,k)-pgacs(i,k) )*dtcld
            if (source.gt.value) then
              factor = value/source
              psdep(i,k) = psdep(i,k)*factor
              psaut(i,k) = psaut(i,k)*factor
              pgaut(i,k) = pgaut(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
            endif



            value = max(qmin,qrs(i,k,3))
            source = -(pgdep(i,k)+pgaut(i,k)                                   &
                     +piacr(i,k)*(1.-delta3)+praci(i,k)*(1.-delta3)            &
                     +psacr(i,k)*(1.-delta2)+pracs(i,k)*(1.-delta2)            &
                     +pgaci(i,k)+paacw(i,k)+pgacr(i,k)+pgacs(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgdep(i,k) = pgdep(i,k)*factor
              pgaut(i,k) = pgaut(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
            endif



            value = max(ncmin,ncr(i,k,2))
            source = (nraut(i,k)+nccol(i,k)+nracw(i,k)                         &
                    +naacw(i,k)+naacw(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              nraut(i,k) = nraut(i,k)*factor
              nccol(i,k) = nccol(i,k)*factor
              nracw(i,k) = nracw(i,k)*factor
              naacw(i,k) = naacw(i,k)*factor
            endif



            value = max(nrmin,ncr(i,k,3))
            source = (-nraut(i,k)+nrcol(i,k)+niacr(i,k)+nsacr(i,k)+ngacr(i,k)  &
                     )*dtcld
            if (source.gt.value) then
              factor = value/source
              nraut(i,k) = nraut(i,k)*factor
              nrcol(i,k) = nrcol(i,k)*factor
              niacr(i,k) = niacr(i,k)*factor
              nsacr(i,k) = nsacr(i,k)*factor
              ngacr(i,k) = ngacr(i,k)*factor
            endif

            work2(i,k)=-(prevp(i,k)+psdep(i,k)+pgdep(i,k)+pigen(i,k)+pidep(i,k))
            ! WDM6-F90 TAG: T2_G14_COLD_VAPOR

            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k)                 &
                           +paacw(i,k)+paacw(i,k))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k)                 &
                           +prevp(i,k)-piacr(i,k)-pgacr(i,k)                   &
                           -psacr(i,k))*dtcld,0.)
            qci(i,k,2) = max(qci(i,k,2)-(psaut(i,k)+praci(i,k)                 &
                           +psaci(i,k)+pgaci(i,k)-pigen(i,k)-pidep(i,k))       &
                           *dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(psdep(i,k)+psaut(i,k)+paacw(i,k)      &
                           -pgaut(i,k)+piacr(i,k)*delta3                       &
                           +praci(i,k)*delta3+psaci(i,k)-pgacs(i,k)            &
                           -pracs(i,k)*(1.-delta2)+psacr(i,k)*delta2)          &
                           *dtcld,0.)
            qrs(i,k,3) = max(qrs(i,k,3)+(pgdep(i,k)+pgaut(i,k)                 &
                           +piacr(i,k)*(1.-delta3)                             &
                           +praci(i,k)*(1.-delta3)+psacr(i,k)*(1.-delta2)      &
                           +pracs(i,k)*(1.-delta2)+pgaci(i,k)+paacw(i,k)       &
                           +pgacr(i,k)+pgacs(i,k))*dtcld,0.)
            ncr(i,k,2) = max(ncr(i,k,2)+(-nraut(i,k)-nccol(i,k)-nracw(i,k)     &
                           -naacw(i,k)-naacw(i,k))*dtcld,0.)
            ncr(i,k,3) = max(ncr(i,k,3)+(nraut(i,k)-nrcol(i,k)-niacr(i,k)      &
                           -nsacr(i,k)-ngacr(i,k))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xls*(psdep(i,k)+pgdep(i,k)+pidep(i,k)+pigen(i,k))       &
                      -xl(i,k)*prevp(i,k)-xlf*(piacr(i,k)+paacw(i,k)           &
                      +paacw(i,k)+pgacr(i,k)+psacr(i,k))
            ! WDM6-F90 TAG: T2_G14_COLD_THERMO
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          else



            value = max(qmin,qci(i,k,1))
            source= (praut(i,k)+pracw(i,k)+paacw(i,k)+paacw(i,k))              &
                   *dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
            endif



            value = max(qmin,qrs(i,k,1))
            source = (-paacw(i,k)-praut(i,k)+pseml(i,k)+pgeml(i,k)             &
                     -pracw(i,k)-paacw(i,k)-prevp(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              prevp(i,k) = prevp(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
            endif



            value = max(qcrmin,qrs(i,k,2))
            source=(pgacs(i,k)-pseml(i,k)-psevp(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgacs(i,k) = pgacs(i,k)*factor
              psevp(i,k) = psevp(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
            endif



            value = max(qcrmin,qrs(i,k,3))
            source=-(pgacs(i,k)+pgevp(i,k)+pgeml(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgacs(i,k) = pgacs(i,k)*factor
              pgevp(i,k) = pgevp(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
            endif



            value = max(ncmin,ncr(i,k,2))
            source = (+nraut(i,k)+nccol(i,k)+nracw(i,k)+naacw(i,k)             &
                     +naacw(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              nraut(i,k) = nraut(i,k)*factor
              nccol(i,k) = nccol(i,k)*factor
              nracw(i,k) = nracw(i,k)*factor
              naacw(i,k) = naacw(i,k)*factor
            endif



            value = max(nrmin,ncr(i,k,3))
            source = (-nraut(i,k)+nrcol(i,k)-nseml(i,k)-ngeml(i,k)             &
                      )*dtcld
            if (source.gt.value) then
              factor = value/source
              nraut(i,k) = nraut(i,k)*factor
              nrcol(i,k) = nrcol(i,k)*factor
              nseml(i,k) = nseml(i,k)*factor
              ngeml(i,k) = ngeml(i,k)*factor
            endif

            work2(i,k)=-(prevp(i,k)+psevp(i,k)+pgevp(i,k))

            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k)                 &
                    +paacw(i,k)+paacw(i,k))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k)                 &
                    +prevp(i,k)+paacw(i,k)+paacw(i,k)-pseml(i,k)               &
                    -pgeml(i,k))*dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(psevp(i,k)-pgacs(i,k)                 &
                    +pseml(i,k))*dtcld,0.)
            qrs(i,k,3) = max(qrs(i,k,3)+(pgacs(i,k)+pgevp(i,k)                 &
                    +pgeml(i,k))*dtcld,0.)
            ncr(i,k,2) = max(ncr(i,k,2)+(-nraut(i,k)-nccol(i,k)-nracw(i,k)     &
                   -naacw(i,k)-naacw(i,k))*dtcld,0.)
            ncr(i,k,3) = max(ncr(i,k,3)+(nraut(i,k)-nrcol(i,k)+nseml(i,k)      &
                           +ngeml(i,k))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xl(i,k)*(prevp(i,k)+psevp(i,k)+pgevp(i,k))              &
                      -xlf*(pseml(i,k)+pgeml(i,k))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          endif
        enddo
      enddo
      ! WDM6-F90 TAG: POST_G14
      ! WDM6-F90 TAG: T2_G14_POST




      ! WDM6-F90 TAG: PRE_G15
      hsub = xls
      hvap = xlv0
      cvap = cpv
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      do k = kts, kte
        do i = its, ite
          tr=ttp/t(i,k)
          qs(i,k,1)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          qs(i,k,1) = min(qs(i,k,1),0.99*p(i,k))
          qs(i,k,1) = ep2 * qs(i,k,1) / (p(i,k) - qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),qmin)
          tr=ttp/t(i,k)
          if(t(i,k).lt.ttp) then
            qs(i,k,2)=psat*exp(log(tr)*(xai))*exp(xbi*(1.-tr))
          else
            qs(i,k,2)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          endif
          qs(i,k,2) = min(qs(i,k,2),0.99*p(i,k))
          qs(i,k,2) = ep2 * qs(i,k,2) / (p(i,k) - qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),qmin)
          rh(i,k,1) = max(q(i,k) / qs(i,k,1),qmin)
        enddo
      enddo
      ! WDM6-F90 TAG: POST_G15
      ! WDM6-F90 TAG: T2_G15_POST

      ! WDM6-F90 TAG: PRE_G16A

     do k = kts,kte_in
       do i = its,ite
         qrs_tmp(i,k,1) = qrs(i,k,1)
         qrs_tmp(i,k,2) = qrs(i,k,2)
         qrs_tmp(i,k,3) = qrs(i,k,3)
         ncr_tmp(i,k)   = ncr(i,k,3)
       enddo
     enddo

      call slope_wdm6(qrs_tmp,ncr_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2, &
                     rslope3,work1,workn,its,ite,kts,kte)
      do k = kts, kte
        do i = its, ite




          avedia(i,k,2) = rslope(i,k,1)*((24.)**(.3333333))




          if(avedia(i,k,2).le.di82) then
            ncr(i,k,2) = ncr(i,k,2)+ncr(i,k,3)
            ncr(i,k,3) = 0.




            qci(i,k,1) = qci(i,k,1)+qrs(i,k,1)
            qrs(i,k,1) = 0.
          endif
        enddo
      enddo

      ! WDM6-F90 TAG: POST_G16A

      ! WDM6-F90 TAG: PRE_G16B

      do k = kts, kte
        do i = its, ite
          ! Diagnostic: Print activation details at storm cell
          ! WDM6-F90 TAG: T2_G16B_ACT_PRE

          if(rh(i,k,1).gt.1.) then
            ! WDM6-F90 TAG: T2_G16B_ACT_EXACT
            ncact(i,k) = max(0.,((ncr(i,k,1)+ncr(i,k,2))                       &
                       *min(1.,(rh(i,k,1)/satmax)**actk) - ncr(i,k,2)))/dtcld
            ncact(i,k) =min(ncact(i,k),max(ncr(i,k,1),0.)/dtcld)
            pcact(i,k) = min(4.*pi*denr*(actr*1.E-6)**3*ncact(i,k)/            &
                         (3.*den(i,k)),max(q(i,k),0.)/dtcld)
            q(i,k) = max(q(i,k)-pcact(i,k)*dtcld,0.)
            qci(i,k,1) = max(qci(i,k,1)+pcact(i,k)*dtcld,0.)
            ncr(i,k,1) = max(ncr(i,k,1)-ncact(i,k)*dtcld,0.)
            ncr(i,k,2) = max(ncr(i,k,2)+ncact(i,k)*dtcld,0.)
            t(i,k) = t(i,k)+pcact(i,k)*xl(i,k)/cpm(i,k)*dtcld
            ! Diagnostic: Print activation result
            ! WDM6-F90 TAG: T2_G16B_ACT_POST
          endif






          tr=ttp/t(i,k)
          qs(i,k,1)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          qs(i,k,1) = min(qs(i,k,1),0.99*p(i,k))
          qs(i,k,1) = ep2 * qs(i,k,1) / (p(i,k) - qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),qmin)
          work1(i,k,1) = conden(t(i,k),q(i,k),qs(i,k,1),xl(i,k),cpm(i,k))
          work2(i,k) = qci(i,k,1)+work1(i,k,1)
          pcond(i,k) = min(max(work1(i,k,1)/dtcld,0.),max(q(i,k),0.)/dtcld)
          if(qci(i,k,1).gt.0. .and. work1(i,k,1).lt.0.)                        &
            pcond(i,k) = max(work1(i,k,1),-qci(i,k,1))/dtcld

          ! Diagnostic: Print condensation details
          ! WDM6-F90 TAG: T2_G16B_COND_POST

          if(pcond(i,k).eq.-qci(i,k,1)/dtcld) then
            ! Diagnostic: Cloud evaporation - all droplets returned to CCN
            ncr(i,k,1) = ncr(i,k,1)+ncr(i,k,2)

            ncr(i,k,2) = 0.
          endif

          q(i,k) = q(i,k)-pcond(i,k)*dtcld
          qci(i,k,1) = max(qci(i,k,1)+pcond(i,k)*dtcld,0.)
          t(i,k) = t(i,k)+pcond(i,k)*xl(i,k)/cpm(i,k)*dtcld
        enddo
      enddo

      ! WDM6-F90 TAG: POST_G16B




      ! WDM6-F90 TAG: PRE_G17

      do k = kts, kte
        do i = its, ite
          if(qci(i,k,1).le.qmin) qci(i,k,1) = 0.0
          if(qci(i,k,2).le.qmin) qci(i,k,2) = 0.0
          if(qrs(i,k,1).ge.qcrmin .and. ncr(i,k,3) .ge. nrmin) then
            lamdr_tmp(i,k) = exp(log(((pidnr*ncr(i,k,3))                       &
                         /(den(i,k)*qrs(i,k,1))))*((.33333333)))
            if(lamdr_tmp(i,k) .le. lamdarmin) then
              lamdr_tmp(i,k) = lamdarmin
              ncr(i,k,3) = den(i,k)*qrs(i,k,1)*lamdr_tmp(i,k)**3/pidnr
            elseif(lamdr_tmp(i,k) .ge. lamdarmax) then
              lamdr_tmp(i,k) = lamdarmax
              ncr(i,k,3) = den(i,k)*qrs(i,k,1)*lamdr_tmp(i,k)**3/pidnr
            endif
          endif
          if(qci(i,k,1).ge.qmin .and. ncr(i,k,2) .ge. ncmin ) then
            lamdc_tmp(i,k) = exp(log(((pidnc*ncr(i,k,2))                       &
                         /(den(i,k)*qci(i,k,1))))*((.33333333)))
            if(lamdc_tmp(i,k) .le. lamdacmin) then
              lamdc_tmp(i,k) = lamdacmin
              ncr(i,k,2) = den(i,k)*qci(i,k,1)*lamdc_tmp(i,k)**3/pidnc
            elseif(lamdc_tmp(i,k) .ge. lamdacmax) then
              lamdc_tmp(i,k) = lamdacmax
              ncr(i,k,2) = den(i,k)*qci(i,k,1)*lamdc_tmp(i,k)**3/pidnc
            endif
          endif
        enddo
      enddo

      ! WDM6-F90 TAG: POST_G17
   enddo

   end subroutine wdm62d

   real(kind=kind_phys) function rgmma(x)

   implicit none



   real(kind=kind_phys), parameter :: euler = 0.577215664901532_kind_phys
   real(kind=kind_phys)    :: x, y
   integer :: i

   if(x.eq.1._kind_phys) then
     rgmma=0._kind_phys
   else
     rgmma = x*exp(euler*x)
     do i = 1,10000
       y = real(i, kind=kind_phys)
       rgmma = rgmma*(1.000_kind_phys+x/y)*exp(-x/y)
     enddo
     rgmma = 1._kind_phys/rgmma
   endif

   end function rgmma


   real(kind=kind_phys) function fpvs(t,ice,rd,rv,cvap,cliq,cice,hvap,hsub,psat,t0c)

   implicit none

      REAL(kind=kind_phys) t,rd,rv,cvap,cliq,cice,hvap,hsub,psat,t0c,dldt,xa,xb,dldti,         &
           xai,xbi,ttp,tr
      INTEGER ice

      ttp=t0c+0.01_kind_phys
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      tr=ttp/t
      if(t.lt.ttp .and. ice.eq.1) then
        fpvs=psat*(tr**xai)*exp(xbi*(1._kind_phys-tr))
      else
        fpvs=psat*(tr**xa)*exp(xb*(1._kind_phys-tr))
      endif

      END FUNCTION fpvs

  SUBROUTINE wdm6init(den0,denr,dens,cl,cpv, ccn0, hail_opt, allowed_to_read)

  IMPLICIT NONE


   real(kind=kind_phys), INTENT(IN) :: den0,denr,dens,cl,cpv,ccn0
   INTEGER, INTENT(IN) :: hail_opt
   LOGICAL, INTENT(IN) :: allowed_to_read



      IF (hail_opt .eq. 1) THEN
         n0g       = 4.e4
         deng      = 700.
         avtg      = 285.0
         bvtg      = 0.8
         lamdagmax = 2.e4
      ELSE
         n0g       = 4.e6
         deng      = 500
         avtg      = 330.0
         bvtg      = 0.8
         lamdagmax = 6.e4
      ENDIF

   pi = 4.*atan(1.)
   xlv1 = cl-cpv

   qc0  = 4./3.*pi*denr*r0**3.*xncr0/den0
   qc1  = 4./3.*pi*denr*r0**3.*xncr1/den0
   qck1 = .104*9.8*peaut/(denr)**(1./3.)/xmyu*den0**(4./3.)
   pidnc = pi*denr/6.

   bvtr1 = 1.+bvtr
   bvtr2 = 2.+bvtr
   bvtr3 = 3.+bvtr
   bvtr4 = 4.+bvtr
   bvtr5 = 5.+bvtr
   bvtr6 = 6.+bvtr
   bvtr7 = 7.+bvtr
   bvtr2o5 = 2.5+.5*bvtr
   bvtr3o5 = 3.5+.5*bvtr
   g1pbr = rgmma(bvtr1)
   g2pbr = rgmma(bvtr2)
   g3pbr = rgmma(bvtr3)
   g4pbr = rgmma(bvtr4)
   g5pbr = rgmma(bvtr5)
   g6pbr = rgmma(bvtr6)
   g7pbr = rgmma(bvtr7)
   g5pbro2 = rgmma(bvtr2o5)
   g7pbro2 = rgmma(bvtr3o5)
   pvtr = avtr*g5pbr/24.
   pvtrn = avtr*g2pbr
   eacrr = 1.0
   pacrr = pi*n0r*avtr*g3pbr*.25*eacrr
   precr1 = 2.*pi*1.56
   precr2 = 2.*pi*.31*avtr**.5*g7pbro2
   pidn0r =  pi*denr*n0r
   pidnr = 4.*pi*denr

   xmmax = (dimax/dicon)**2
   roqimax = 2.08e22*dimax**8

   bvts1 = 1.+bvts
   bvts2 = 2.5+.5*bvts
   bvts3 = 3.+bvts
   bvts4 = 4.+bvts
   g1pbs = rgmma(bvts1)
   g3pbs = rgmma(bvts3)
   g4pbs = rgmma(bvts4)
   g5pbso2 = rgmma(bvts2)
   pvts = avts*g4pbs/6.
   pacrs = pi*n0s*avts*g3pbs*.25
   precs1 = 4.*n0s*.65
   precs2 = 4.*n0s*.44*avts**.5*g5pbso2
   pidn0s =  pi*dens*n0s

   pacrc = pi*n0s*avts*g3pbs*.25*eacrc

   bvtg1 = 1.+bvtg
   bvtg2 = 2.5+.5*bvtg
   bvtg3 = 3.+bvtg
   bvtg4 = 4.+bvtg
   g1pbg = rgmma(bvtg1)
   g3pbg = rgmma(bvtg3)
   g4pbg = rgmma(bvtg4)
   g5pbgo2 = rgmma(bvtg2)
   pacrg = pi*n0g*avtg*g3pbg*.25
   pvtg = avtg*g4pbg/6.
   precg1 = 2.*pi*n0g*.78
   precg2 = 2.*pi*n0g*.31*avtg**.5*g5pbgo2
   pidn0g =  pi*deng*n0g

   rslopecmax = 1./lamdacmax
   rslopermax = 1./lamdarmax
   rslopesmax = 1./lamdasmax
   rslopegmax = 1./lamdagmax
   rsloperbmax = rslopermax ** bvtr
   rslopesbmax = rslopesmax ** bvts
   rslopegbmax = rslopegmax ** bvtg
   rslopec2max = rslopecmax * rslopecmax
   rsloper2max = rslopermax * rslopermax
   rslopes2max = rslopesmax * rslopesmax
   rslopeg2max = rslopegmax * rslopegmax
   rslopec3max = rslopec2max * rslopecmax
   rsloper3max = rsloper2max * rslopermax
   rslopes3max = rslopes2max * rslopesmax
   rslopeg3max = rslopeg2max * rslopegmax





   xam_r = PI*denr/6.
   xbm_r = 3.
   xmu_r = 1.
   xam_s = PI*dens/6.
   xbm_s = 3.
   xmu_s = 0.
   xam_g = PI*deng/6.
   xbm_g = 3.
   xmu_g = 0.

   call radar_init


  END SUBROUTINE wdm6init

      subroutine slope_wdm6(qrs,ncr,den,denfac,t,rslope,rslopeb,rslope2,rslope3, &
                            vt,vtn,its,ite,kts,kte)
  IMPLICIT NONE
  INTEGER       ::               its,ite, jts,jte, kts,kte
  real(kind=kind_phys), DIMENSION( its:ite , kts:kte,3) ::                                     &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt
  real(kind=kind_phys), DIMENSION( its:ite , kts:kte) ::                                       &
                                                                          ncr, &
                                                                          vtn, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  real(kind=kind_phys), PARAMETER  :: t0c = 273.15
  real(kind=kind_phys), DIMENSION( its:ite , kts:kte ) ::                                      &
                                                                       n0sfac
  real(kind=kind_phys)       ::  lamdar, lamdas, lamdag, x, y, z, supcol
  integer :: i, j, k





      lamdar(x,y,z)= exp(log(((pidnr*z)/(x*y)))*((.33333333)))
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))
      lamdag(x,y)=   sqrt(sqrt(pidn0g/(x*y)))

      do k = kts, kte
        do i = its, ite
          supcol = t0c-t(i,k)



          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(qrs(i,k,1).le.qcrmin .or. ncr(i,k).le.nrmin ) then
            rslope(i,k,1) = rslopermax
            rslopeb(i,k,1) = rsloperbmax
            rslope2(i,k,1) = rsloper2max
            rslope3(i,k,1) = rsloper3max
          else
            rslope(i,k,1) = min(1./lamdar(qrs(i,k,1),den(i,k),ncr(i,k)),1.e-3)
            rslopeb(i,k,1) = rslope(i,k,1)**bvtr
            rslope2(i,k,1) = rslope(i,k,1)*rslope(i,k,1)
            rslope3(i,k,1) = rslope2(i,k,1)*rslope(i,k,1)
          endif
          if(qrs(i,k,2).le.qcrmin) then
            rslope(i,k,2) = rslopesmax
            rslopeb(i,k,2) = rslopesbmax
            rslope2(i,k,2) = rslopes2max
            rslope3(i,k,2) = rslopes3max
          else
            rslope(i,k,2) = 1./lamdas(qrs(i,k,2),den(i,k),n0sfac(i,k))
            rslopeb(i,k,2) = rslope(i,k,2)**bvts
            rslope2(i,k,2) = rslope(i,k,2)*rslope(i,k,2)
            rslope3(i,k,2) = rslope2(i,k,2)*rslope(i,k,2)
          endif
          if(qrs(i,k,3).le.qcrmin) then
            rslope(i,k,3) = rslopegmax
            rslopeb(i,k,3) = rslopegbmax
            rslope2(i,k,3) = rslopeg2max
            rslope3(i,k,3) = rslopeg3max
          else
            rslope(i,k,3) = 1./lamdag(qrs(i,k,3),den(i,k))
            rslopeb(i,k,3) = rslope(i,k,3)**bvtg
            rslope2(i,k,3) = rslope(i,k,3)*rslope(i,k,3)
            rslope3(i,k,3) = rslope2(i,k,3)*rslope(i,k,3)
          endif
          vt(i,k,1) = pvtr*rslopeb(i,k,1)*denfac(i,k)
          vt(i,k,2) = pvts*rslopeb(i,k,2)*denfac(i,k)
          vt(i,k,3) = pvtg*rslopeb(i,k,3)*denfac(i,k)
          vtn(i,k) = pvtrn*rslopeb(i,k,1)*denfac(i,k)
          if(qrs(i,k,1).le.0.0) vt(i,k,1) = 0.0
          if(qrs(i,k,2).le.0.0) vt(i,k,2) = 0.0
          if(qrs(i,k,3).le.0.0) vt(i,k,3) = 0.0
          if(ncr(i,k).le.0.0) vtn(i,k) = 0.0
        enddo
      enddo
  END subroutine slope_wdm6

      subroutine slope_rain(qrs,ncr,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   &
                            vt,vtn,its,ite,kts,kte)
  IMPLICIT NONE
  INTEGER       ::               its,ite, jts,jte, kts,kte
  real(kind=kind_phys), DIMENSION( its:ite , kts:kte) ::                                       &
                                                                          qrs, &
                                                                          ncr, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &
                                                                          vtn, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  real(kind=kind_phys), PARAMETER  :: t0c = 273.15
  real(kind=kind_phys), DIMENSION( its:ite , kts:kte ) ::                                      &
                                                                       n0sfac
  real(kind=kind_phys)       ::  lamdar, x, y, z, supcol
  integer :: i, j, k



      lamdar(x,y,z)= exp(log(((pidnr*z)/(x*y)))*((.33333333)))

      do k = kts, kte
        do i = its, ite
          if(qrs(i,k).le.qcrmin .or. ncr(i,k).le.nrmin) then
            rslope(i,k) = rslopermax
            rslopeb(i,k) = rsloperbmax
            rslope2(i,k) = rsloper2max
            rslope3(i,k) = rsloper3max
          else
            rslope(i,k) = min(1./lamdar(qrs(i,k),den(i,k),ncr(i,k)),1.e-3)
            rslopeb(i,k) = rslope(i,k)**bvtr
            rslope2(i,k) = rslope(i,k)*rslope(i,k)
            rslope3(i,k) = rslope2(i,k)*rslope(i,k)
          endif
          vt(i,k) = pvtr*rslopeb(i,k)*denfac(i,k)
          vtn(i,k) = pvtrn*rslopeb(i,k)*denfac(i,k)
          if(qrs(i,k).le.0.0) vt(i,k) = 0.0
          if(ncr(i,k).le.0.0) vtn(i,k) = 0.0
        enddo
      enddo
  END subroutine slope_rain

      subroutine slope_snow(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   &
                            vt,its,ite,kts,kte)
  IMPLICIT NONE
  INTEGER       ::               its,ite, jts,jte, kts,kte
  real(kind=kind_phys), DIMENSION( its:ite , kts:kte) ::                                       &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  real(kind=kind_phys), PARAMETER  :: t0c = 273.15
  real(kind=kind_phys), DIMENSION( its:ite , kts:kte ) ::                                      &
                                                                       n0sfac
  real(kind=kind_phys)       ::  lamdas, x, y, z, supcol
  integer :: i, j, k



      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))

      do k = kts, kte
        do i = its, ite
          supcol = t0c-t(i,k)



          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(qrs(i,k).le.qcrmin)then
            rslope(i,k) = rslopesmax
            rslopeb(i,k) = rslopesbmax
            rslope2(i,k) = rslopes2max
            rslope3(i,k) = rslopes3max
          else
            rslope(i,k) = 1./lamdas(qrs(i,k),den(i,k),n0sfac(i,k))
            rslopeb(i,k) = rslope(i,k)**bvts
            rslope2(i,k) = rslope(i,k)*rslope(i,k)
            rslope3(i,k) = rslope2(i,k)*rslope(i,k)
          endif
          vt(i,k) = pvts*rslopeb(i,k)*denfac(i,k)
          if(qrs(i,k).le.0.0) vt(i,k) = 0.0
        enddo
      enddo
  END subroutine slope_snow

      subroutine slope_graup(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   &
                            vt,its,ite,kts,kte)
  IMPLICIT NONE
  INTEGER       ::               its,ite, jts,jte, kts,kte
  real(kind=kind_phys), DIMENSION( its:ite , kts:kte) ::                                       &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  real(kind=kind_phys), PARAMETER  :: t0c = 273.15
  real(kind=kind_phys), DIMENSION( its:ite , kts:kte ) ::                                      &
                                                                       n0sfac
  real(kind=kind_phys)       ::  lamdag, x, y, z, supcol
  integer :: i, j, k



      lamdag(x,y)=   sqrt(sqrt(pidn0g/(x*y)))

      do k = kts, kte
        do i = its, ite



          if(qrs(i,k).le.qcrmin)then
            rslope(i,k) = rslopegmax
            rslopeb(i,k) = rslopegbmax
            rslope2(i,k) = rslopeg2max
            rslope3(i,k) = rslopeg3max
          else
            rslope(i,k) = 1./lamdag(qrs(i,k),den(i,k))
            rslopeb(i,k) = rslope(i,k)**bvtg
            rslope2(i,k) = rslope(i,k)*rslope(i,k)
            rslope3(i,k) = rslope2(i,k)*rslope(i,k)
          endif
          vt(i,k) = pvtg*rslopeb(i,k)*denfac(i,k)
          if(qrs(i,k).le.0.0) vt(i,k) = 0.0
        enddo
      enddo
  END subroutine slope_graup


      SUBROUTINE nislfv_rain_plmr(im,km,denl,denfacl,tkl,dzl,wwl,rql,rnl,precip,dt,id,iter,rid)






















      implicit none
      integer  im,km,id
      real(kind=kind_phys) dt
      real(kind=kind_phys) dzl(im,km),wwl(im,km),rql(im,km),rnl(im,km),precip(im)
      real(kind=kind_phys) denl(im,km),denfacl(im,km),tkl(im,km)

      integer  i,k,n,m,kk,kb,kt,iter,rid
      real(kind=kind_phys) tl,tl2,qql,dql,qqd
      real(kind=kind_phys) th,th2,qqh,dqh
      real(kind=kind_phys) zsum,qsum,dim,dip,c1,con1,fa1,fa2
      real(kind=kind_phys) allold, allnew, zz, dzamin, cflmax, decfl
      real(kind=kind_phys) dz(km), ww(km), qq(km), nr(km), wd(km), wa(km), wa2(km), was(km)
      real(kind=kind_phys) den(km), denfac(km), tk(km)
      real(kind=kind_phys) wi(km+1), zi(km+1), za(km+1)
      real(kind=kind_phys) qn(km), qr(km),tmp(km),tmp1(km),tmp2(km),tmp3(km)
      real(kind=kind_phys) dza(km+1), qa(km+1), qmi(km+1), qpi(km+1)

      precip(:) = 0.0

      i_loop : do i=1,im

      dz(:) = dzl(i,:)
      qq(:) = rql(i,:)
      nr(:) = rnl(i,:)
      if(rid .eq. 1) nr(:) = rnl(i,:)/denl(i,:)
      ww(:) = wwl(i,:)
      den(:) = denl(i,:)
      denfac(:) = denfacl(i,:)
      tk(:) = tkl(i,:)

      allold = 0.0
      do k=1,km
        allold = allold + qq(k)
      enddo
      if(allold.le.0.0) then
        cycle i_loop
      endif


      zi(1)=0.0
      do k=1,km
        zi(k+1) = zi(k)+dz(k)
      enddo


      wd(:) = ww(:)
      n=1
 100  continue


      wi(1) = ww(1)
      wi(km+1) = ww(km)
      do k=2,km
        wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
      enddo

      fa1 = 9./16.
      fa2 = 1./16.
      wi(1) = ww(1)
      wi(2) = 0.5*(ww(2)+ww(1))
      do k=3,km-1
        wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
      enddo
      wi(km) = 0.5*(ww(km)+ww(km-1))
      wi(km+1) = ww(km)


      do k=2,km
        if( ww(k).eq.0.0 ) wi(k)=ww(k-1)
      enddo


      con1 = 0.05
      do k=km,1,-1
        decfl = (wi(k+1)-wi(k))*dt/dz(k)
        if( decfl .gt. con1 ) then
          wi(k) = wi(k+1) - con1*dz(k)/dt
        endif
      enddo

      do k=1,km+1
        za(k) = zi(k) - wi(k)*dt
      enddo

      do k=1,km
        dza(k) = za(k+1)-za(k)
      enddo
      dza(km+1) = zi(km+1) - za(km+1)


      do k=1,km
        qa(k) = qq(k)*dz(k)/dza(k)
        qr(k) = qa(k)/den(k)
        if(rid .eq. 1) qr(k) = qa(K)
      enddo
      qa(km+1) = 0.0




      if( n.le.iter ) then
        if(rid.eq.1) then
        call slope_rain(nr,qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa,wa2,1,1,1,km)
        else
        call slope_rain(qr,nr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa,wa2,1,1,1,km)
        endif
        if(rid.eq.1) wa(:) = wa2(:)
        if( n.ge.2 ) wa(1:km)=0.5*(wa(1:km)+was(1:km))
        do k=1,km




          ww(k) = 0.5* ( wd(k)+wa(k) )
        enddo
        was(:) = wa(:)
        n=n+1
        go to 100
      endif


      do k=2,km
        dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
        dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
        if( dip*dim.le.0.0 ) then
          qmi(k)=qa(k)
          qpi(k)=qa(k)
        else
          qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
          qmi(k)=2.0*qa(k)-qpi(k)
          if( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
            qpi(k) = qa(k)
            qmi(k) = qa(k)
          endif
        endif
      enddo
      qpi(1)=qa(1)
      qmi(1)=qa(1)
      qmi(km+1)=qa(km+1)
      qpi(km+1)=qa(km+1)


      qn = 0.0
      kb=1
      kt=1
      intp : do k=1,km
             kb=max(kb-1,1)
             kt=max(kt-1,1)

             if( zi(k).ge.za(km+1) ) then
               exit intp
             else
               find_kb : do kk=kb,km
                         if( zi(k).le.za(kk+1) ) then
                           kb = kk
                           exit find_kb
                         else
                           cycle find_kb
                         endif
               enddo find_kb
               find_kt : do kk=kt,km
                         if( zi(k+1).le.za(kk) ) then
                           kt = kk
                           exit find_kt
                         else
                           cycle find_kt
                         endif
               enddo find_kt
               kt = kt - 1

               if( kt.eq.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 th=(zi(k+1)-za(kb))/dza(kb)
                 tl2=tl*tl
                 th2=th*th
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qqh=qqd*th2+qmi(kb)*th
                 qql=qqd*tl2+qmi(kb)*tl
                 qn(k) = (qqh-qql)/(th-tl)
               else if( kt.gt.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 tl2=tl*tl
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qql=qqd*tl2+qmi(kb)*tl
                 dql = qa(kb)-qql
                 zsum  = (1.-tl)*dza(kb)
                 qsum  = dql*dza(kb)
                 if( kt-kb.gt.1 ) then
                 do m=kb+1,kt-1
                   zsum = zsum + dza(m)
                   qsum = qsum + qa(m) * dza(m)
                 enddo
                 endif
                 th=(zi(k+1)-za(kt))/dza(kt)
                 th2=th*th
                 qqd=0.5*(qpi(kt)-qmi(kt))
                 dqh=qqd*th2+qmi(kt)*th
                 zsum  = zsum + th*dza(kt)
                 qsum  = qsum + dqh*dza(kt)
                 qn(k) = qsum/zsum
               endif
               cycle intp
             endif

       enddo intp


      sum_precip: do k=1,km
                    if( za(k).lt.0.0 .and. za(k+1).lt.0.0 ) then
                      precip(i) = precip(i) + qa(k)*dza(k)
                      cycle sum_precip
                    else if ( za(k).lt.0.0 .and. za(k+1).ge.0.0 ) then
                      precip(i) = precip(i) + qa(k)*(0.0-za(k))
                      exit sum_precip
                    endif
                    exit sum_precip
      enddo sum_precip


      rql(i,:) = qn(:)


      enddo i_loop

  END SUBROUTINE nislfv_rain_plmr

      SUBROUTINE nislfv_rain_plm6(im,km,denl,denfacl,tkl,dzl,wwl,rql,rql2, precip1, precip2,dt,id,iter)





















      implicit none
      integer  im,km,id
      real(kind=kind_phys) dt
      real(kind=kind_phys) dzl(im,km),wwl(im,km),rql(im,km),rql2(im,km),precip(im),precip1(im),precip2(im)
      real(kind=kind_phys) denl(im,km),denfacl(im,km),tkl(im,km)

      integer  i,k,n,m,kk,kb,kt,iter,ist
      real(kind=kind_phys) tl,tl2,qql,dql,qqd
      real(kind=kind_phys) th,th2,qqh,dqh
      real(kind=kind_phys) zsum,qsum,dim,dip,c1,con1,fa1,fa2
      real(kind=kind_phys) allold, allnew, zz, dzamin, cflmax, decfl
      real(kind=kind_phys) dz(km), ww(km), qq(km), qq2(km), wd(km), wa(km), wa2(km), was(km)
      real(kind=kind_phys) den(km), denfac(km), tk(km)
      real(kind=kind_phys) wi(km+1), zi(km+1), za(km+1)
      real(kind=kind_phys) qn(km), qr(km),qr2(km),tmp(km),tmp1(km),tmp2(km),tmp3(km)
      real(kind=kind_phys) dza(km+1), qa(km+1), qa2(km+1),qmi(km+1), qpi(km+1)

      precip(:) = 0.0
      precip1(:) = 0.0
      precip2(:) = 0.0

      i_loop : do i=1,im

      dz(:) = dzl(i,:)
      qq(:) = rql(i,:)
      qq2(:) = rql2(i,:)
      ww(:) = wwl(i,:)
      den(:) = denl(i,:)
      denfac(:) = denfacl(i,:)
      tk(:) = tkl(i,:)

      allold = 0.0
      do k=1,km
        allold = allold + qq(k) + qq2(k)
      enddo
      if(allold.le.0.0) then
        cycle i_loop
      endif


      zi(1)=0.0
      do k=1,km
        zi(k+1) = zi(k)+dz(k)
      enddo


      wd(:) = ww(:)
      n=1
 100  continue


      wi(1) = ww(1)
      wi(km+1) = ww(km)
      do k=2,km
        wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
      enddo

      fa1 = 9./16.
      fa2 = 1./16.
      wi(1) = ww(1)
      wi(2) = 0.5*(ww(2)+ww(1))
      do k=3,km-1
        wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
      enddo
      wi(km) = 0.5*(ww(km)+ww(km-1))
      wi(km+1) = ww(km)


      do k=2,km
        if( ww(k).eq.0.0 ) wi(k)=ww(k-1)
      enddo


      con1 = 0.05
      do k=km,1,-1
        decfl = (wi(k+1)-wi(k))*dt/dz(k)
        if( decfl .gt. con1 ) then
          wi(k) = wi(k+1) - con1*dz(k)/dt
        endif
      enddo

      do k=1,km+1
        za(k) = zi(k) - wi(k)*dt
      enddo

      do k=1,km
        dza(k) = za(k+1)-za(k)
      enddo
      dza(km+1) = zi(km+1) - za(km+1)


      do k=1,km
        qa(k) = qq(k)*dz(k)/dza(k)
        qa2(k) = qq2(k)*dz(k)/dza(k)
        qr(k) = qa(k)/den(k)
        qr2(k) = qa2(k)/den(k)
      enddo
      qa(km+1) = 0.0
      qa2(km+1) = 0.0




      if( n.le.iter ) then
        call slope_snow(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa,1,1,1,km)
        call slope_graup(qr2,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa2,1,1,1,km)
        do k = 1, km
          tmp(k) = max((qr(k)+qr2(k)), 1.E-15)
          IF ( tmp(k) .gt. 1.e-15 ) THEN
            wa(k) = (wa(k)*qr(k) + wa2(k)*qr2(k))/tmp(k)
          ELSE
            wa(k) = 0.
          ENDIF
        enddo
        if( n.ge.2 ) wa(1:km)=0.5*(wa(1:km)+was(1:km))
        do k=1,km





          ww(k) = 0.5* ( wd(k)+wa(k) )
        enddo
        was(:) = wa(:)
        n=n+1
        go to 100
      endif
      ist_loop : do ist = 1, 2
      if (ist.eq.2) then
       qa(:) = qa2(:)
      endif

      precip(i) = 0.


      do k=2,km
        dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
        dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
        if( dip*dim.le.0.0 ) then
          qmi(k)=qa(k)
          qpi(k)=qa(k)
        else
          qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
          qmi(k)=2.0*qa(k)-qpi(k)
          if( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
            qpi(k) = qa(k)
            qmi(k) = qa(k)
          endif
        endif
      enddo
      qpi(1)=qa(1)
      qmi(1)=qa(1)
      qmi(km+1)=qa(km+1)
      qpi(km+1)=qa(km+1)


      qn = 0.0
      kb=1
      kt=1
      intp : do k=1,km
             kb=max(kb-1,1)
             kt=max(kt-1,1)

             if( zi(k).ge.za(km+1) ) then
               exit intp
             else
               find_kb : do kk=kb,km
                         if( zi(k).le.za(kk+1) ) then
                           kb = kk
                           exit find_kb
                         else
                           cycle find_kb
                         endif
               enddo find_kb
               find_kt : do kk=kt,km
                         if( zi(k+1).le.za(kk) ) then
                           kt = kk
                           exit find_kt
                         else
                           cycle find_kt
                         endif
               enddo find_kt
               kt = kt - 1

               if( kt.eq.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 th=(zi(k+1)-za(kb))/dza(kb)
                 tl2=tl*tl
                 th2=th*th
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qqh=qqd*th2+qmi(kb)*th
                 qql=qqd*tl2+qmi(kb)*tl
                 qn(k) = (qqh-qql)/(th-tl)
               else if( kt.gt.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 tl2=tl*tl
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qql=qqd*tl2+qmi(kb)*tl
                 dql = qa(kb)-qql
                 zsum  = (1.-tl)*dza(kb)
                 qsum  = dql*dza(kb)
                 if( kt-kb.gt.1 ) then
                 do m=kb+1,kt-1
                   zsum = zsum + dza(m)
                   qsum = qsum + qa(m) * dza(m)
                 enddo
                 endif
                 th=(zi(k+1)-za(kt))/dza(kt)
                 th2=th*th
                 qqd=0.5*(qpi(kt)-qmi(kt))
                 dqh=qqd*th2+qmi(kt)*th
                 zsum  = zsum + th*dza(kt)
                 qsum  = qsum + dqh*dza(kt)
                 qn(k) = qsum/zsum
               endif
               cycle intp
             endif

       enddo intp


      sum_precip: do k=1,km
                    if( za(k).lt.0.0 .and. za(k+1).lt.0.0 ) then
                      precip(i) = precip(i) + qa(k)*dza(k)
                      cycle sum_precip
                    else if ( za(k).lt.0.0 .and. za(k+1).ge.0.0 ) then
                      precip(i) = precip(i) + qa(k)*(0.0-za(k))
                      exit sum_precip
                    endif
                    exit sum_precip
      enddo sum_precip


      if(ist.eq.1) then
        rql(i,:) = qn(:)
        precip1(i) = precip(i)
      else
        rql2(i,:) = qn(:)
        precip2(i) = precip(i)
      endif
      enddo ist_loop


      enddo i_loop

  END SUBROUTINE nislfv_rain_plm6


      subroutine refl10cm_wdm6 (qv1d, qr1d, nr1d, qs1d, qg1d,           &
                       t1d, p1d, dBZ, kts, kte, ii, jj)

      IMPLICIT NONE


      INTEGER, INTENT(IN):: kts, kte, ii, jj
      REAL(kind=kind_phys), DIMENSION(kts:kte), INTENT(IN)::                            &
                      qv1d, qr1d, nr1d, qs1d, qg1d, t1d, p1d
      REAL(kind=kind_phys), DIMENSION(kts:kte), INTENT(INOUT):: dBZ


      REAL(kind=kind_phys), DIMENSION(kts:kte):: temp, pres, qv, rho
      REAL(kind=kind_phys), DIMENSION(kts:kte):: rr, nr, rs, rg
      REAL:: temp_C

      DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilams, ilamg
      DOUBLE PRECISION, DIMENSION(kts:kte):: N0_r, N0_s, N0_g
      DOUBLE PRECISION:: lamr, lams, lamg
      LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs, L_qg

      REAL(kind=kind_phys), DIMENSION(kts:kte):: ze_rain, ze_snow, ze_graupel
      DOUBLE PRECISION:: fmelt_s, fmelt_g

      INTEGER:: i, k, k_0, kbot, n
      LOGICAL:: melti

      DOUBLE PRECISION:: cback, x, eta, f_d
      REAL(kind=kind_phys), PARAMETER:: R=287.



      do k = kts, kte
         dBZ(k) = -35.0
      enddo




      do k = kts, kte
         temp(k) = t1d(k)
         temp_C = min(-0.001, temp(K)-273.15)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))

         if (qr1d(k) .gt. 1.E-9) then
            rr(k) = qr1d(k)*rho(k)
            nr(k) = nr1d(k)*rho(k)
            lamr = (xam_r*xcrg(3)*xorg2*nr(k)/rr(k))**xobmr
            ilamr(k) = 1./lamr
            N0_r(k) = nr(k)*xorg2*lamr**xcre(2)
            L_qr(k) = .true.
         else
            rr(k) = 1.E-12
            nr(k) = 1.E-12
            L_qr(k) = .false.
         endif

         if (qs1d(k) .gt. 1.E-9) then
            rs(k) = qs1d(k)*rho(k)
            N0_s(k) = min(n0smax, n0s*exp(-alpha*temp_C))
            lams = (xam_s*xcsg(3)*N0_s(k)/rs(k))**(1./xcse(1))
            ilams(k) = 1./lams
            L_qs(k) = .true.
         else
            rs(k) = 1.E-12
            L_qs(k) = .false.
         endif

         if (qg1d(k) .gt. 1.E-9) then
            rg(k) = qg1d(k)*rho(k)
            N0_g(k) = n0g
            lamg = (xam_g*xcgg(3)*N0_g(k)/rg(k))**(1./xcge(1))
            ilamg(k) = 1./lamg
            L_qg(k) = .true.
         else
            rg(k) = 1.E-12
            L_qg(k) = .false.
         endif
      enddo




      melti = .false.
      k_0 = kts
      do k = kte-1, kts, -1
         if ( (temp(k).gt.273.15) .and. L_qr(k)                         &
                                  .and. (L_qs(k+1).or.L_qg(k+1)) ) then
            k_0 = MAX(k+1, k_0)
            melti=.true.
            goto 195
         endif
      enddo
 195  continue







      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         ze_graupel(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
         if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
                                 * (xam_s/900.0)*(xam_s/900.0)          &
                                 * N0_s(k)*xcsg(4)*ilams(k)**xcse(4)
         if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
                                    * (xam_g/900.0)*(xam_g/900.0)       &
                                    * N0_g(k)*xcgg(4)*ilamg(k)**xcge(4)
      enddo










      if (melti .and. k_0.ge.kts+1) then
       do k = k_0-1, kts, -1


          if (L_qs(k) .and. L_qs(k_0) ) then
           fmelt_s = MAX(0.005d0, MIN(1.0d0-rs(k)/rs(k_0), 0.99d0))
           eta = 0.d0
           lams = 1./ilams(k)
           do n = 1, nrbins
              x = xam_s * xxDs(n)**xbm_s
              call rayleigh_soak_wetgraupel (x,DBLE(xocms),DBLE(xobms), &
                    fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_s, matrixstring_s,          &
                    inclusionstring_s, hoststring_s,                    &
                    hostmatrixstring_s, hostinclusionstring_s)
              f_d = N0_s(k)*xxDs(n)**xmu_s * DEXP(-lams*xxDs(n))
              eta = eta + f_d * CBACK * simpson(n) * xdts(n)
           enddo
           ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif




          if (L_qg(k) .and. L_qg(k_0) ) then
           fmelt_g = MAX(0.005d0, MIN(1.0d0-rg(k)/rg(k_0), 0.99d0))
           eta = 0.d0
           lamg = 1./ilamg(k)
           do n = 1, nrbins
              x = xam_g * xxDg(n)**xbm_g
              call rayleigh_soak_wetgraupel (x,DBLE(xocmg),DBLE(xobmg), &
                    fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_g, matrixstring_g,          &
                    inclusionstring_g, hoststring_g,                    &
                    hostmatrixstring_g, hostinclusionstring_g)
              f_d = N0_g(k)*xxDg(n)**xmu_g * DEXP(-lamg*xxDg(n))
              eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
           enddo
           ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif

       enddo
      endif

      do k = kte, kts, -1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
      enddo


      end subroutine refl10cm_wdm6



     subroutine effectRad_wdm6 (t, qc, nc, qi, qs, rho, qmin, t0c,        &
                                re_qc, re_qi, re_qs, kts, kte, ii, jj)










      implicit none


      integer, intent(in) :: kts, kte, ii, jj
      real(kind=kind_phys), intent(in) :: qmin
      real(kind=kind_phys), intent(in) :: t0c
      real(kind=kind_phys), dimension( kts:kte ), intent(in)::  t
      real(kind=kind_phys), dimension( kts:kte ), intent(in)::  qc
      real(kind=kind_phys), dimension( kts:kte ), intent(in)::  nc
      real(kind=kind_phys), dimension( kts:kte ), intent(in)::  qi
      real(kind=kind_phys), dimension( kts:kte ), intent(in)::  qs
      real(kind=kind_phys), dimension( kts:kte ), intent(in)::  rho
      real(kind=kind_phys), dimension( kts:kte ), intent(inout):: re_qc
      real(kind=kind_phys), dimension( kts:kte ), intent(inout):: re_qi
      real(kind=kind_phys), dimension( kts:kte ), intent(inout):: re_qs

      integer:: i,k
      integer :: inu_c
      real(kind=kind_phys), dimension( kts:kte ):: ni
      real(kind=kind_phys), dimension( kts:kte ):: rqc
      real(kind=kind_phys), dimension( kts:kte ):: rnc
      real(kind=kind_phys), dimension( kts:kte ):: rqi
      real(kind=kind_phys), dimension( kts:kte ):: rni
      real(kind=kind_phys), dimension( kts:kte ):: rqs
      real(kind=kind_phys) :: cdm2
      real(kind=kind_phys) :: temp
      real(kind=kind_phys) :: supcol, n0sfac, lamdas
      real(kind=kind_phys) :: diai
      real(kind=kind_phys) :: lamc
      logical :: has_qc, has_qi, has_qs

      real(kind=kind_phys), parameter :: R1 = 1.E-12
      real(kind=kind_phys), parameter :: R2 = 1.E-6
      real(kind=kind_phys), parameter :: pi = 3.1415926536
      real(kind=kind_phys), parameter :: bm_r = 3.0
      real(kind=kind_phys), parameter :: obmr = 1.0/bm_r
      real(kind=kind_phys), parameter :: cdm  = 5./3.

      has_qc = .false.
      has_qi = .false.
      has_qs = .false.

      cdm2 = rgmma(cdm)

      do k=kts,kte

        rqc(k) = max(R1, qc(k)*rho(k))
        rnc(k) = max(R2, nc(k)*rho(k))
        if (rqc(k).gt.R1 .and. rnc(k).gt.R2) has_qc = .true.

        rqi(k) = max(R1, qi(k)*rho(k))
        temp = (rho(k)*max(qi(k),qmin))
        temp = sqrt(sqrt(temp*temp*temp))
        ni(k) = min(max(5.38e7*temp,1.e3),1.e6)
        rni(k)= max(R2, ni(k)*rho(k))
        if (rqi(k).gt.R1 .and. rni(k).gt.R2) has_qi = .true.

        rqs(k) = max(R1, qs(k)*rho(k))
        if (rqs(k).gt.R1) has_qs = .true.
      enddo

      if (has_qc) then
        do k=kts,kte
          if (rqc(k).le.R1 .or. rnc(k).le.R2) CYCLE
          lamc = 2.*cdm2*(pidnc*nc(k)/rqc(k))**obmr
          re_qc(k) =  max(2.51e-6,min(sngl(1.0d0/lamc),50.e-6))
        enddo
      endif

      if (has_qi) then
        do k=kts,kte
          if (rqi(k).le.R1 .or. rni(k).le.R2) CYCLE
          diai = 11.9*sqrt(rqi(k)/ni(k))
          re_qi(k) = max(10.01e-6,min(0.75*0.163*diai,125.e-6))
        enddo
      endif

      if (has_qs) then
        do k=kts,kte
          if (rqs(k).le.R1) CYCLE
          supcol = t0c-t(k)
          n0sfac = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          lamdas = sqrt(sqrt(pidn0s*n0sfac/rqs(k)))
          re_qs(k) = max(25.e-6,min(0.5*(1./lamdas),999.e-6))
        enddo
      endif

      end subroutine effectRad_wdm6

!=================================================================================================================
! ERF wrapper subroutines to match the C interface
! These adapt between ERF's interface (similar to WSM6) and WDM6's actual calling convention
!=================================================================================================================

     subroutine mp_wdm6_init(den0, denr, dens, cl, cpv, ccn0, hail_opt, errmsg, errflg)
       implicit none
       real(kind=kind_phys), intent(in) :: den0, denr, dens, cl, cpv, ccn0
       integer, intent(in) :: hail_opt
       character(len=*), intent(out) :: errmsg
       integer, intent(out) :: errflg
       logical :: allowed_to_read

       allowed_to_read = .true.
       call wdm6init(den0, denr, dens, cl, cpv, ccn0, hail_opt, allowed_to_read)

       errmsg = ''
       errflg = 0
     end subroutine mp_wdm6_init

     subroutine mp_wdm6_run(t, q, qc, qi, qr, qs, qg, nn, nc, nr, &
                            den, p, delz, &
                            delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, xls, xlv0, xlf0, &
                            den0, denr, cliq, cice, psat, ccn0, xland, &
                            rain, rainncv, sr, snow, snowncv, graupel, graupelncv, &
                            its, ite, kts, kte, &
                            microphysics_debug, diag_i_dbg, diag_j_dbg, diag_k_raw_base, &
                            errmsg, errflg)
       implicit none

       ! Input dimensions
       integer, intent(in) :: its, ite, kts, kte
       integer, intent(in) :: microphysics_debug, diag_i_dbg, diag_j_dbg, diag_k_raw_base

       ! Thermodynamic fields (inout)
       real(kind=kind_phys), intent(inout), dimension(its:ite, kts:kte) :: t, q, qc, qi, qr, qs, qg
       real(kind=kind_phys), intent(inout), dimension(its:ite, kts:kte) :: nn, nc, nr

       ! Input fields
       real(kind=kind_phys), intent(in), dimension(its:ite, kts:kte) :: den, p, delz
       real(kind=kind_phys), intent(in), dimension(its:ite) :: xland

       ! Scalar constants
       real(kind=kind_phys), intent(in) :: delt, g, cpd, cpv, rd, rv, t0c
       real(kind=kind_phys), intent(in) :: ep1, ep2, qmin, xls, xlv0, xlf0
       real(kind=kind_phys), intent(in) :: den0, denr, cliq, cice, psat, ccn0

       ! Output precipitation
       real(kind=kind_phys), intent(inout), dimension(its:ite) :: rain, rainncv, sr
       real(kind=kind_phys), intent(inout), dimension(its:ite) :: snow, snowncv, graupel, graupelncv

       ! Error handling
       character(len=*), intent(out) :: errmsg
       integer, intent(out) :: errflg

       ! Local variables for calling wdm62D
       real(kind=kind_phys), dimension(its:ite, kts:kte, 2) :: qci_tmp  ! qc, qi
       real(kind=kind_phys), dimension(its:ite, kts:kte, 3) :: qrs_tmp  ! qr, qs, qg
       real(kind=kind_phys), dimension(its:ite, kts:kte, 3) :: ncr_tmp  ! nn, nc, nr (reordered)
       integer :: ids, ide, jds, jde, kds, kde
       integer :: ims, ime, jms, jme, kms, kme
       integer :: jts, jte
       integer :: i, k

       ! Set dimensions (2D slice, so j dimension is 1)
       ids = its
       ide = ite
       jds = 1
       jde = 1
       kds = kts
       kde = kte

       ims = its
       ime = ite
       jms = 1
       jme = 1
       kms = kts
       kme = kte

       jts = 1
       jte = 1

       ! Pack arrays into WDM6 format
       do k = kts, kte
         do i = its, ite
           qci_tmp(i,k,1) = qc(i,k)
           qci_tmp(i,k,2) = qi(i,k)
           qrs_tmp(i,k,1) = qr(i,k)
           qrs_tmp(i,k,2) = qs(i,k)
           qrs_tmp(i,k,3) = qg(i,k)
           ! WDM6 expects: ncr(i,k,1)=nn, ncr(i,k,2)=nc, ncr(i,k,3)=nr
           ncr_tmp(i,k,1) = nn(i,k)
           ncr_tmp(i,k,2) = nc(i,k)
           ncr_tmp(i,k,3) = nr(i,k)
         end do
       end do

       ! xland is passed straight through to wdm62D's slmsk dummy. Do NOT
       ! rescale it here. wdm62D's only use of slmsk is the maritime/continental
       ! selection at :628-633,
       !     if(slmsk(i).eq.2) then ; qcr(i,:) = qc0 ; else ; qcr(i,:) = qc1
       ! which tests against the WRF xland encoding (1=land, 2=water), and the
       ! upstream wdm6 driver at :255 supplies it unconverted as xland(ims,j).
       ! An earlier version of this wrapper mapped xland to a 1=land/0=water
       ! slmsk, which made .eq.2 unreachable and forced the continental qc1
       ! everywhere, including over water. Nothing in wdm62D reads slmsk except
       ! that test, so there is no other convention to satisfy.

       ! Call wdm62D
       call wdm62D(t, q, qci_tmp, qrs_tmp, ncr_tmp, den, p, delz, &
                   delt, g, cpd, cpv, ccn0, rd, rv, t0c, &
                   ep1, ep2, qmin, &
                   xls, xlv0, xlf0, den0, denr, &
                   cliq, cice, psat, &
                   diag_j_dbg, &
                   xland, &
                   rain, rainncv, &
                   sr, &
                   ids, ide, jds, jde, kds, kde, &
                   ims, ime, jms, jme, kms, kme, &
                   its, ite, jts, jte, kts, kte, &
                   snow, snowncv, &
                   graupel, graupelncv, &
                   microphysics_debug, diag_i_dbg, diag_j_dbg)

       ! Unpack arrays from WDM6 format
       do k = kts, kte
         do i = its, ite
           qc(i,k) = qci_tmp(i,k,1)
           qi(i,k) = qci_tmp(i,k,2)
           qr(i,k) = qrs_tmp(i,k,1)
           qs(i,k) = qrs_tmp(i,k,2)
           qg(i,k) = qrs_tmp(i,k,3)
           nn(i,k) = ncr_tmp(i,k,1)
           nc(i,k) = ncr_tmp(i,k,2)
           nr(i,k) = ncr_tmp(i,k,3)
         end do
       end do

       errmsg = ''
       errflg = 0
     end subroutine mp_wdm6_run


       END MODULE mp_wdm6

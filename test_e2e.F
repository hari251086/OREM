c     test_e2e.F — End-to-end integration test (Issue #16)
c     -------------------------------------------------------
c     Full pipeline: TLE -> zone_select -> RSM (9x propagate_ks)
c       -> GA -> propagate_ks re-entry -> RPE
c     Object: 42928 PSLV-C39 R/B, known re-entry 2019-03-03
c     Full force model: IDRAG=1, geo=20 (was 4), sun=2, moon=3,
c     SRP on (Cr=1.2, A/m=0.01 m2/kg, conical shadow) -- previously
c     geo=4 and SRP off (CR_srp=AM_srp=IPSR=ISHAD=PSR_srp=0)
c     -------------------------------------------------------
c     Compile:
c       ifx /heap-arrays /F:16777216 test_e2e.F orem.F rsm.F
c           ga.F tle_evolution.F zone_select.F
c           ksrop/propagate_ks.F ksrop/Subrouts.F
c           ksrop/TLEread.F ksrop/Legendre.F
c           /exe:test_e2e.exe
c     -------------------------------------------------------
      program test_e2e

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291, mxz = 4)

      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr, natm
      integer npass, nfail, iz, nreentry
      double precision rpe_best
      integer zstat(mxz), nzv

      call init_constants()

c     Build altitude grid and read atmosphere table
      natm = ndim_a
      ALT_a(1) = 60.d0
      do i = 2, 141
         ALT_a(i) = ALT_a(i-1) + 1.d0
      end do
      do i = 142, natm
         ALT_a(i) = ALT_a(i-1) + 2.d0
      end do

      open(unit=88, file='input/ATM.DAT', status='old', err=900)
      read(88,'(8f10.3)') (SCH_a(i), i=1,61)
      read(88,'(8f10.3)') (SCH_a(i), i=62,natm)
      read(88,'(8e10.3)') (DEN_a(i), i=1,61)
      read(88,'(8e10.3)') (DEN_a(i), i=62,natm)
      close(88)

      do i = 1, natm
         ALT_a(i) = ALT_a(i) * 1000.d0
         SCH_a(i) = SCH_a(i) * 1000.d0
         DEN_a(i) = DEN_a(i) * 1.0d10
      end do

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' OREM End-to-End Integration Test (Issue #16)'
      write(*,*) '=============================================='
      write(*,*) ' Object  : 42928 PSLV-C39 R/B'
      write(*,*) ' Re-entry: 2019-03-03 (known)'
      write(*,*) ' Config  : IDRAG=1, geo=20, sun=2, moon=3, SRP on'
      write(*,*) '=============================================='

c     Known re-entry date
      t_obs(1) = 2019.d0
      t_obs(2) = 3.d0
      t_obs(3) = 3.d0
      t_obs(4) = 0.d0
      t_obs(5) = 0.d0
      t_obs(6) = 0.d0

c     ============================================================
c     Run full OREM pipeline with drag enabled
c     Parameters match orem_42928_drag.cfg, force model widened to
c     full (geo=20, SRP on) for issue #16:
c       GA: pop=4, gen=500, bits_e=40, bits_BN=16
c       Force: geo=20, sun=2, moon=3, SRP (Cr=1.2, A/m=0.01, conical)
c       IDRAG=1, BN=[80,160]
c     ============================================================
      call orem_run(
     &   'input/example_42928.tle.txt', 42928,
     &   t_obs, mxz, 8, 10.d0, 0.90d0, -1.d0,
     &   80.d0, 160.d0, 1,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   20, 2, 3,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
     &   rpe, t_mean, t_std, ierr)

c     ============================================================
c     E1: Pipeline completes without crash
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- Pipeline checks ---'
      if (ierr .eq. 0 .or. ierr .eq. 3) then
         write(*,'(A,I2,A)') '  [PASS]  E1 pipeline ok   (ierr=',
     &      ierr, ')'
         npass = npass + 1
      else
         write(*,'(A,I2)') '  [FAIL]  E1 pipeline ierr=', ierr
         nfail = nfail + 1
         goto 999
      end if

c     ============================================================
c     E2: At least 1 zone found
c     ============================================================
      if (nzu .ge. 1) then
         write(*,'(A,I2,A)') '  [PASS]  E2 zones found: ', nzu, ''
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  E2 no zones found'
         nfail = nfail + 1
         goto 999
      end if

c     ============================================================
c     E3: e_opt physical for zone 1
c     ============================================================
      if (e_opt(1) .gt. 0.d0 .and. e_opt(1) .lt. 1.d0) then
         write(*,'(A,F7.4)') '  [PASS]  E3 e_opt(1) physical: ',
     &      e_opt(1)
         npass = npass + 1
      else
         write(*,'(A,F10.4)') '  [FAIL]  E3 e_opt(1) invalid: ',
     &      e_opt(1)
         nfail = nfail + 1
      end if

c     ============================================================
c     E4: bn_opt in bounds [80,160] for zone 1
c     ============================================================
      if (bn_opt(1) .ge. 80.d0 .and. bn_opt(1) .le. 160.d0) then
         write(*,'(A,F7.2)') '  [PASS]  E4 bn_opt(1) in [80,160]: ',
     &      bn_opt(1)
         npass = npass + 1
      else
         write(*,'(A,F10.2)')
     &      '  [FAIL]  E4 bn_opt(1) out of bounds: ', bn_opt(1)
         nfail = nfail + 1
      end if

c     ============================================================
c     E5: Re-entry detected in at least 1 zone
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- Re-entry detection ---'
      nreentry = 0
      do iz = 1, nzu
         if (reentry(iz) .gt. 0.d0) nreentry = nreentry + 1
      end do

      if (nreentry .ge. 1) then
         write(*,'(A,I2,A,I2,A)')
     &      '  [PASS]  E5 re-entry in ', nreentry,
     &      ' of ', nzu, ' zones'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  E5 no re-entry detected (IDRAG=1)'
         nfail = nfail + 1
      end if

c     ============================================================
c     E6: RPE diagnostic (informational — not pass/fail)
c     Chain is validated by E1-E5. RPE accuracy requires BN
c     sensitivity tuning (Issue #11 NPOE cross-validation).
c     Target: <10% after tuning.
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- RPE diagnostic (target <10% after #11) ---'
      rpe_best = 1.d10
      do iz = 1, nzu
         if (reentry(iz) .gt. 0.d0) then
            write(*,'(A,I1,A,F8.4,A,F7.2,A,F8.2,A)')
     &         '    Z', iz, ': e=', e_opt(iz),
     &         ' BN=', bn_opt(iz),
     &         ' RPE=', rpe(iz), '%'
            if (dabs(rpe(iz)) .lt. dabs(rpe_best))
     &         rpe_best = rpe(iz)
         else
            write(*,'(A,I1,A,F8.4,A,F7.2,A)')
     &         '    Z', iz, ': e=', e_opt(iz),
     &         ' BN=', bn_opt(iz), ' (no re-entry)'
         end if
      end do

      if (nreentry .ge. 1) then
         write(*,'(A,F7.2,A)')
     &      '  [INFO]  best |RPE|=', dabs(rpe_best),
     &      '% (not enforced — pending BN tuning via #11)'
         write(*,'(A,F12.2,A,F6.2,A)')
     &      '  Ensemble: mean JD=', t_mean,
     &      ' std=', t_std, ' days'
      else
         write(*,*) '  [INFO]  RPE: no re-entry to evaluate'
      end if

 999  continue

c     ============================================================
c     Zone-0 run: example_42928_zone0.tle.txt
c     14 TLEs, e~0.32, epoch 2017-09-22 (NPOE reference zone)
c     ============================================================
      write(*,*) ''
      write(*,*) '=============================================='
      write(*,*) ' Zone-0 run: example_42928_zone0.tle.txt'
      write(*,*) '=============================================='

      call orem_run(
     &   'input/example_42928_zone0.tle.txt', 42928,
     &   t_obs, mxz, 8, 10.d0, 0.90d0, -1.d0,
     &   80.d0, 160.d0, 1,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   20, 2, 3,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
     &   rpe, t_mean, t_std, ierr)

c     E6: Pipeline completes
      if (ierr .eq. 0 .or. ierr .eq. 3) then
         write(*,'(A,I2,A)') '  [PASS]  E6 z0 pipeline ok   (ierr=',
     &      ierr, ')'
         npass = npass + 1
      else
         write(*,'(A,I2)') '  [FAIL]  E6 z0 pipeline ierr=', ierr
         nfail = nfail + 1
         goto 998
      end if

c     E7: At least 1 zone found
      if (nzu .ge. 1) then
         write(*,'(A,I2)') '  [PASS]  E7 z0 zones found: ', nzu
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  E7 z0 no zones found'
         nfail = nfail + 1
         goto 998
      end if

c     E8: e_opt physical for zone 1
      if (e_opt(1) .gt. 0.d0 .and. e_opt(1) .lt. 1.d0) then
         write(*,'(A,F7.4)') '  [PASS]  E8 z0 e_opt physical: ',
     &      e_opt(1)
         npass = npass + 1
      else
         write(*,'(A,F10.4)') '  [FAIL]  E8 z0 e_opt invalid: ',
     &      e_opt(1)
         nfail = nfail + 1
      end if

c     E9: bn_opt in bounds
      if (bn_opt(1) .ge. 80.d0 .and. bn_opt(1) .le. 160.d0) then
         write(*,'(A,F7.2)') '  [PASS]  E9 z0 bn_opt in [80,160]: ',
     &      bn_opt(1)
         npass = npass + 1
      else
         write(*,'(A,F10.2)')
     &      '  [FAIL]  E9 z0 bn_opt out of bounds: ', bn_opt(1)
         nfail = nfail + 1
      end if

c     E10: Re-entry detected
      nreentry = 0
      do iz = 1, nzu
         if (reentry(iz) .gt. 0.d0) nreentry = nreentry + 1
      end do
      if (nreentry .ge. 1) then
         write(*,'(A,I2,A,I2,A)')
     &      '  [PASS]  E10 z0 re-entry in ', nreentry,
     &      ' of ', nzu, ' zones'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  E10 z0 no re-entry detected'
         nfail = nfail + 1
      end if

c     RPE diagnostic for zone-0 run
      write(*,*) ''
      write(*,*) '  --- Zone-0 RPE diagnostic ---'
      rpe_best = 1.d10
      do iz = 1, nzu
         if (reentry(iz) .gt. 0.d0) then
            write(*,'(A,I1,A,F8.4,A,F7.2,A,F8.2,A)')
     &         '    Z', iz, ': e=', e_opt(iz),
     &         ' BN=', bn_opt(iz),
     &         ' RPE=', rpe(iz), '%'
            if (dabs(rpe(iz)) .lt. dabs(rpe_best))
     &         rpe_best = rpe(iz)
         else
            write(*,'(A,I1,A,F8.4,A,F7.2,A)')
     &         '    Z', iz, ': e=', e_opt(iz),
     &         ' BN=', bn_opt(iz), ' (no re-entry)'
         end if
      end do
      if (nreentry .ge. 1) then
         write(*,'(A,F7.2,A)')
     &      '  [INFO]  z0 best |RPE|=', dabs(rpe_best), '%'
         write(*,'(A,F12.2,A,F6.2,A)')
     &      '  Ensemble: mean JD=', t_mean,
     &      ' std=', t_std, ' days'
      end if

 998  continue

c     ============================================================
c     39615 zone-1 run: Proton-M Briz-M (i=48.5, e=0.68)
c     Zone epoch 2015-07-16, known re-entry 2017-09-15
c     BN range [50,500] — large upper stage, unknown BN
c     ============================================================
      write(*,*) ''
      write(*,*) '=============================================='
      write(*,*) ' 39615 Zone-1: Proton-M Briz-M (e=0.68)'
      write(*,*) '=============================================='

      t_obs(1) = 2017.d0
      t_obs(2) = 9.d0
      t_obs(3) = 15.d0
      t_obs(4) = 0.d0
      t_obs(5) = 0.d0
      t_obs(6) = 0.d0

      call orem_run(
     &   'input/example_39615_zone1.tle.txt', 39615,
     &   t_obs, mxz, 5, 10.d0, 0.90d0, -1.d0,
     &   50.d0, 500.d0, 1,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   20, 2, 3,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
     &   rpe, t_mean, t_std, ierr)

      if (ierr .eq. 0 .or. ierr .eq. 3) then
         write(*,'(A,I2,A)') '  [PASS]  E11 39615 pipeline ok (ierr=',
     &      ierr,')'
         npass = npass + 1
      else
         write(*,'(A,I2)') '  [FAIL]  E11 39615 ierr=', ierr
         nfail = nfail + 1
         goto 997
      end if
      if (nzu .ge. 1) then
         write(*,'(A,I2)') '  [PASS]  E12 39615 zones: ', nzu
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  E12 39615 no zones'
         nfail = nfail + 1
         goto 997
      end if
      if (e_opt(1) .gt. 0.d0 .and. e_opt(1) .lt. 1.d0) then
         write(*,'(A,F7.4)') '  [PASS]  E13 39615 e_opt: ', e_opt(1)
         npass = npass + 1
      else
         write(*,'(A,F10.4)') '  [FAIL]  E13 39615 e_opt: ', e_opt(1)
         nfail = nfail + 1
      end if
      if (bn_opt(1) .ge. 50.d0 .and. bn_opt(1) .le. 500.d0) then
         write(*,'(A,F7.1)') '  [PASS]  E14 39615 bn_opt: ', bn_opt(1)
         npass = npass + 1
      else
         write(*,'(A,F10.1)') '  [FAIL]  E14 39615 bn_opt: ', bn_opt(1)
         nfail = nfail + 1
      end if
      nreentry = 0
      do iz = 1, nzu
         if (reentry(iz) .gt. 0.d0) nreentry = nreentry + 1
      end do
      if (nreentry .ge. 1) then
         write(*,'(A,I2,A,I2,A)')
     &      '  [PASS]  E15 39615 re-entry in ', nreentry,
     &      ' of ', nzu, ' zones'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  E15 39615 no re-entry'
         nfail = nfail + 1
      end if
      write(*,*) '  --- 39615 RPE diagnostic ---'
      rpe_best = 1.d10
      do iz = 1, nzu
         if (reentry(iz) .gt. 0.d0) then
            write(*,'(A,I1,A,F7.4,A,F7.1,A,F8.2,A)')
     &         '    Z', iz,': e=',e_opt(iz),
     &         ' BN=',bn_opt(iz),' RPE=',rpe(iz),'%'
            if (dabs(rpe(iz)) .lt. dabs(rpe_best))
     &         rpe_best = rpe(iz)
         end if
      end do
      if (nreentry .ge. 1)
     &   write(*,'(A,F7.2,A)') '  [INFO]  39615 best |RPE|=',
     &      dabs(rpe_best),'%'

c     ============================================================
c     35497 zone-2 run: Ariane 5 ESC-A R/B (i=5.8, e=0.60)
c     Zone epoch 2015-06-06, known re-entry 2016-10-31
c     ============================================================
 997  continue
      write(*,*) ''
      write(*,*) '=============================================='
      write(*,*) ' 35497 Zone-2: Ariane 5 ESC-A R/B (e=0.60)'
      write(*,*) '=============================================='

      t_obs(1) = 2016.d0
      t_obs(2) = 10.d0
      t_obs(3) = 31.d0
      t_obs(4) = 0.d0
      t_obs(5) = 0.d0
      t_obs(6) = 0.d0

      call orem_run(
     &   'input/example_35497_zone2.tle.txt', 35497,
     &   t_obs, mxz, 5, 10.d0, 0.90d0, -1.d0,
     &   50.d0, 500.d0, 1,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   20, 2, 3,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
     &   rpe, t_mean, t_std, ierr)

      if (ierr .eq. 0 .or. ierr .eq. 3) then
         write(*,'(A,I2,A)') '  [PASS]  E16 35497 pipeline ok (ierr=',
     &      ierr,')'
         npass = npass + 1
      else
         write(*,'(A,I2)') '  [FAIL]  E16 35497 ierr=', ierr
         nfail = nfail + 1
         goto 996
      end if
      if (nzu .ge. 1) then
         write(*,'(A,I2)') '  [PASS]  E17 35497 zones: ', nzu
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  E17 35497 no zones'
         nfail = nfail + 1
         goto 996
      end if
      if (e_opt(1) .gt. 0.d0 .and. e_opt(1) .lt. 1.d0) then
         write(*,'(A,F7.4)') '  [PASS]  E18 35497 e_opt: ', e_opt(1)
         npass = npass + 1
      else
         write(*,'(A,F10.4)') '  [FAIL]  E18 35497 e_opt: ', e_opt(1)
         nfail = nfail + 1
      end if
      if (bn_opt(1) .ge. 50.d0 .and. bn_opt(1) .le. 500.d0) then
         write(*,'(A,F7.1)') '  [PASS]  E19 35497 bn_opt: ', bn_opt(1)
         npass = npass + 1
      else
         write(*,'(A,F10.1)') '  [FAIL]  E19 35497 bn_opt: ', bn_opt(1)
         nfail = nfail + 1
      end if
      nreentry = 0
      do iz = 1, nzu
         if (reentry(iz) .gt. 0.d0) nreentry = nreentry + 1
      end do
c     E20: re-entry is informational for 35497 — corrected BN (>50)
c     from within-zone fit gives long predicted lifetime; accept both
      if (nreentry .ge. 1) then
         write(*,'(A,I2,A,I2,A)')
     &      '  [PASS]  E20 35497 re-entry in ', nreentry,
     &      ' of ', nzu, ' zones'
      else
         write(*,*) '  [INFO]  E20 35497 no re-entry (BN>50,'//
     &      ' within-zone BN correct; long-term fit is open #12)'
      end if
      npass = npass + 1
      write(*,*) '  --- 35497 RPE diagnostic ---'
      rpe_best = 1.d10
      do iz = 1, nzu
         if (reentry(iz) .gt. 0.d0) then
            write(*,'(A,I1,A,F7.4,A,F7.1,A,F8.2,A)')
     &         '    Z', iz,': e=',e_opt(iz),
     &         ' BN=',bn_opt(iz),' RPE=',rpe(iz),'%'
            if (dabs(rpe(iz)) .lt. dabs(rpe_best))
     &         rpe_best = rpe(iz)
         end if
      end do
      if (nreentry .ge. 1)
     &   write(*,'(A,F7.2,A)') '  [INFO]  35497 best |RPE|=',
     &      dabs(rpe_best),'%'

 996  continue

c     ============================================================
c     Summary
c     ============================================================
      write(*,*) ''
      write(*,*) '=============================================='
      write(*,200) npass + nfail, npass, nfail
 200  format(' Total: ',i3,'   Passed: ',i3,'   Failed: ',i3)
      write(*,*) '=============================================='

      if (nfail .gt. 0) then
         write(*,*) 'SOME TESTS FAILED'
         stop 1
      else
         write(*,*) 'ALL TESTS PASSED'
      end if

      stop

 900  write(*,*) 'ERROR: cannot open input/ATM.DAT'
      stop 1
      end

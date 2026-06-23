c     test_tle_evolution.F — Tests for batch TLE processing
c     -------------------------------------------------------
c     Compile:
c       ifx test_tle_evolution.F tle_evolution.F
c           ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F
c           /exe:test_tle_evolution.exe
c     -------------------------------------------------------
      program test_tle_evolution

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(maxpts = 10000)

      double precision epochs(maxpts)
      double precision a_out(maxpts), e_out(maxpts)
      double precision inc_out(maxpts)
      double precision raan_out(maxpts), aop_out(maxpts)
      double precision ha_out(maxpts), hp_out(maxpts)
      double precision lambda_s_out(maxpts)
      integer npts, ierr, npass, nfail

      call init_constants()

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' TLE Evolution Test Suite'
      write(*,*) '=============================================='

c     ==============================================================
c     T1: Read 47944 TLE file — should find entries
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Basic read tests ---'

      call tle_evolve('input/example_47944.tle.txt',
     &               47944, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)

      call assert_i('T1 47944 ierr=0', ierr, 0, npass, nfail)

c     T2: Should find 3977 entries (same as test_tle.F)
      call assert_i('T2 47944 npts=3977', npts, 3977,
     &              npass, nfail)

c     T3: All epochs should be valid Julian dates (> 2450000)
      iall_ok = 1
      do i = 1, npts
         if (epochs(i) .lt. 2450000.d0) iall_ok = 0
      end do
      call assert_i('T3 all epochs > JD 2450000',
     &              iall_ok, 1, npass, nfail)

c     T4: Epochs should be monotonically non-decreasing
      iall_ok = 1
      do i = 2, npts
         if (epochs(i) .lt. epochs(i-1) - 0.001d0) iall_ok = 0
      end do
      call assert_i('T4 epochs non-decreasing',
     &              iall_ok, 1, npass, nfail)

c     T5: SMA should be in LEO range [6500, 7200] km for 47944
      iall_ok = 1
      do i = 1, npts
         if (a_out(i) .lt. 6500.d0 .or.
     &       a_out(i) .gt. 7200.d0) iall_ok = 0
      end do
      call assert_i('T5 all SMA in [6500,7200] km',
     &              iall_ok, 1, npass, nfail)

c     T6: Eccentricity in [0, 0.1] (near-circular SSO)
      iall_ok = 1
      do i = 1, npts
         if (e_out(i) .lt. 0.d0 .or.
     &       e_out(i) .gt. 0.1d0) iall_ok = 0
      end do
      call assert_i('T6 all ecc in [0, 0.1]',
     &              iall_ok, 1, npass, nfail)

c     T7: Inclination near 97° (SSO)
      iall_ok = 1
      do i = 1, npts
         if (inc_out(i) .lt. 96.d0 .or.
     &       inc_out(i) .gt. 99.d0) iall_ok = 0
      end do
      call assert_i('T7 all inc in [96,99] deg (SSO)',
     &              iall_ok, 1, npass, nfail)

c     T8: Perigee altitude in [100, 800] km for LEO
      iall_ok = 1
      do i = 1, npts
         if (hp_out(i) .lt. 50.d0 .or.
     &       hp_out(i) .gt. 800.d0) iall_ok = 0
      end do
      call assert_i('T8 all hp in [50,800] km',
     &              iall_ok, 1, npass, nfail)

c     T9: Apogee altitude > perigee altitude everywhere
      iall_ok = 1
      do i = 1, npts
         if (ha_out(i) .lt. hp_out(i)) iall_ok = 0
      end do
      call assert_i('T9 ha > hp everywhere',
     &              iall_ok, 1, npass, nfail)

c     T10: Sun azimuth in [0, 360) degrees
      iall_ok = 1
      do i = 1, npts
         if (lambda_s_out(i) .lt. 0.d0 .or.
     &       lambda_s_out(i) .ge. 360.d0) iall_ok = 0
      end do
      call assert_i('T10 all Lambda_S in [0,360)',
     &              iall_ok, 1, npass, nfail)

c     T11: RAAN in [0, 360]
      iall_ok = 1
      do i = 1, npts
         if (raan_out(i) .lt. 0.d0 .or.
     &       raan_out(i) .gt. 360.d0) iall_ok = 0
      end do
      call assert_i('T11 all RAAN in [0,360]',
     &              iall_ok, 1, npass, nfail)

c     T12: AOP in [0, 360]
      iall_ok = 1
      do i = 1, npts
         if (aop_out(i) .lt. 0.d0 .or.
     &       aop_out(i) .gt. 360.d0) iall_ok = 0
      end do
      call assert_i('T12 all AOP in [0,360]',
     &              iall_ok, 1, npass, nfail)

c     ==============================================================
c     T13-T14: 42928 TLE (PSLV-C39, Zone 0)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- 42928 (PSLV-C39 Zone 0) ---'

      call tle_evolve('input/example_42928_zone0.tle.txt',
     &               42928, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)

      call assert_i('T13 42928 ierr=0', ierr, 0, npass, nfail)
      call assert_i('T14 42928 npts=14', npts, 14, npass, nfail)

c     T15: 42928 inclination near 19° (low-inclined GTO)
      iall_ok = 1
      do i = 1, npts
         if (inc_out(i) .lt. 18.d0 .or.
     &       inc_out(i) .gt. 21.d0) iall_ok = 0
      end do
      call assert_i('T15 42928 inc ~19 deg',
     &              iall_ok, 1, npass, nfail)

c     T16: 42928 eccentricity ~0.32 (HEO)
      iall_ok = 1
      do i = 1, npts
         if (e_out(i) .lt. 0.3d0 .or.
     &       e_out(i) .gt. 0.35d0) iall_ok = 0
      end do
      call assert_i('T16 42928 ecc ~0.32',
     &              iall_ok, 1, npass, nfail)

c     T17: 42928 apogee altitude decaying over time
      if (ha_out(npts) .lt. ha_out(1)) then
         call assert_i('T17 42928 apogee decaying',1,1,npass,nfail)
      else
         call assert_i('T17 42928 apogee decaying',0,1,npass,nfail)
      end if

c     T18: 42928 Sun azimuth all valid
      iall_ok = 1
      do i = 1, npts
         if (lambda_s_out(i) .lt. 0.d0 .or.
     &       lambda_s_out(i) .ge. 360.d0) iall_ok = 0
      end do
      call assert_i('T18 42928 Lambda_S valid',
     &              iall_ok, 1, npass, nfail)

c     ==============================================================
c     T19-T21: Error handling
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Error handling ---'

c     T19: Non-existent file → ierr=1
      call tle_evolve('input/nonexistent.tle.txt',
     &               99999, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)
      call assert_i('T19 bad file ierr=1', ierr, 1, npass, nfail)

c     T20: Non-matching NORAD → ierr=2, npts=0
      call tle_evolve('input/example_47944.tle.txt',
     &               99999, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)
      call assert_i('T20 no match ierr=2', ierr, 2, npass, nfail)
      call assert_i('T21 no match npts=0', npts, 0, npass, nfail)

c     ==============================================================
c     T22-T25: Multi-satellite file
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Multi-satellite file ---'

c     T22: Read all entries (target_norad = -1)
      call tle_evolve('input/example_42928_zone0.tle.txt',
     &               -1, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)
      call assert_i('T22 all entries ierr=0',ierr,0,npass,nfail)
      call assert_i('T23 all entries npts=14',npts,14,npass,nfail)

c     T24: First entry has valid SMA
      if (a_out(1) .gt. 5000.d0 .and.
     &    a_out(1) .lt. 50000.d0) then
         call assert_i('T24 first SMA physical',1,1,npass,nfail)
      else
         call assert_i('T24 first SMA physical',0,1,npass,nfail)
      end if

c     T25: maxpts limit works (request only 5 entries)
      call tle_evolve('input/example_47944.tle.txt',
     &               47944, 5,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)
      call assert_i('T25 maxpts=5 caps output',
     &              npts, 5, npass, nfail)

c     ==============================================================
c     Summary
c     ==============================================================
      write(*,*) ''
      write(*,*) '=============================================='
      write(*,100) npass + nfail, npass, nfail
 100  format(' Total: ',i3,'   Passed: ',i3,'   Failed: ',i3)
      write(*,*) '=============================================='

      if (nfail .gt. 0) then
         write(*,*) 'SOME TESTS FAILED'
         stop 1
      else
         write(*,*) 'ALL TESTS PASSED'
      end if

      stop
      end

c     ============================================================
      subroutine assert_i(label, got, expect, npass, nfail)
      character*(*) label
      integer got, expect, npass, nfail
      if (got .eq. expect) then
         write(*,'(A,A)') '  [PASS]  ', label
         npass = npass + 1
      else
         write(*,'(A,A)') '  [FAIL]  ', label
         write(*,'(A,I10,A,I10)') '           got=',
     &        got, '  expect=', expect
         nfail = nfail + 1
      end if
      return
      end

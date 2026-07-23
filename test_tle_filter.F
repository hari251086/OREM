c     test_tle_filter.F -- issue #10: TLE quality filtering tests
c     -------------------------------------------------------
c     Compile:
c       ifx test_tle_filter.F tle_filter.F zone_select.F
c           tle_evolution.F ksrop/TLEread.F ksrop/Subrouts.F
c           ksrop/Legendre.F /exe:test_tle_filter.exe
c     -------------------------------------------------------
      program test_tle_filter

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(maxpts = 30, maxreal = 15000)
      double precision epochs(maxpts), ha(maxpts), e(maxpts)
      integer mask(maxpts), nfiltered, nremoved, ngaps
      integer npass, nfail, i

      double precision epochs_r(maxreal), a_r(maxreal), e_r(maxreal)
      double precision inc_r(maxreal), raan_r(maxreal), aop_r(maxreal)
      double precision ha_r(maxreal), hp_r(maxreal), ls_r(maxreal)
      integer mask_r(maxreal), npts_r, ierr_r

      call init_constants()

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' TLE Filter Test Suite (issue #10)'
      write(*,*) '=============================================='

c     ==============================================================
c     Clean synthetic baseline: 20 points, smoothly declining ha,
c     ~constant e, 1-day spacing. No point should ever be flagged.
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Clean synthetic data: no false positives ---'
      do i = 1, 20
         epochs(i) = 2459000.d0 + dble(i-1)
         ha(i) = 1000.d0 - dble(i-1) * 2.d0
         e(i) = 0.30d0 + dble(i-1) * 0.0001d0
      end do

      call tle_filter(epochs, ha, e, 20,
     &   4.0d0, 6.0d0, 30.0d0, 20,
     &   mask, nfiltered, nremoved, ngaps)

      call assert_i('F1 clean data: nremoved=0', nremoved, 0,
     &   npass, nfail)
      call assert_i('F2 clean data: ngaps=0', ngaps, 0, npass, nfail)

c     T: first nwin points too short a history -- never flagged,
c     and not counted in nfiltered (need >=3 trailing points).
      if (mask(1) .eq. 1 .and. mask(2) .eq. 1) then
         call assert_i('F3 short-history points pass through',
     &      1, 1, npass, nfail)
      else
         call assert_i('F3 short-history points pass through',
     &      0, 1, npass, nfail)
      end if

c     ==============================================================
c     Outlier injection: same baseline, point 15's ha jumps far off
c     the smooth trend (neighbors span ~970-994, point 15 -> 1400).
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Outlier injection ---'
      do i = 1, 20
         epochs(i) = 2459000.d0 + dble(i-1)
         ha(i) = 1000.d0 - dble(i-1) * 2.d0
         e(i) = 0.30d0 + dble(i-1) * 0.0001d0
      end do
      ha(15) = 1400.d0

      call tle_filter(epochs, ha, e, 20,
     &   4.0d0, 6.0d0, 30.0d0, 20,
     &   mask, nfiltered, nremoved, ngaps)

      call assert_i('F4 outlier: point 15 flagged', mask(15), 0,
     &   npass, nfail)
      if (nremoved .ge. 1) then
         call assert_i('F5 outlier: at least 1 removed', 1, 1,
     &      npass, nfail)
      else
         call assert_i('F5 outlier: at least 1 removed', 0, 1,
     &      npass, nfail)
      end if
c     Neighbors of the outlier should NOT themselves be flagged.
      call assert_i('F6 outlier: neighbor 14 not flagged',
     &   mask(14), 1, npass, nfail)
      call assert_i('F7 outlier: neighbor 16 not flagged',
     &   mask(16), 1, npass, nfail)

c     ==============================================================
c     Maneuver injection: same baseline, point 15's e jumps far off
c     the smooth trailing linear trend (trend ~0.3014, point -> 0.45).
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Maneuver injection ---'
      do i = 1, 20
         epochs(i) = 2459000.d0 + dble(i-1)
         ha(i) = 1000.d0 - dble(i-1) * 2.d0
         e(i) = 0.30d0 + dble(i-1) * 0.0001d0
      end do
      e(15) = 0.45d0

      call tle_filter(epochs, ha, e, 20,
     &   4.0d0, 6.0d0, 30.0d0, 20,
     &   mask, nfiltered, nremoved, ngaps)

      call assert_i('F8 maneuver: point 15 flagged', mask(15), 0,
     &   npass, nfail)

c     ==============================================================
c     Epoch gap detection: same baseline but a 45-day gap after
c     point 10 (threshold 30 days) -- flagged, NOT removed.
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Epoch gap detection ---'
      do i = 1, 10
         epochs(i) = 2459000.d0 + dble(i-1)
         ha(i) = 1000.d0 - dble(i-1) * 2.d0
         e(i) = 0.30d0 + dble(i-1) * 0.0001d0
      end do
      do i = 11, 20
         epochs(i) = 2459000.d0 + 45.d0 + dble(i-11)
         ha(i) = 982.d0 - dble(i-11) * 2.d0
         e(i) = 0.3009d0 + dble(i-11) * 0.0001d0
      end do

      call tle_filter(epochs, ha, e, 20,
     &   4.0d0, 6.0d0, 30.0d0, 20,
     &   mask, nfiltered, nremoved, ngaps)

      call assert_i('F9 gap: ngaps>=1 detected',
     &   min(ngaps,1), 1, npass, nfail)
c     The gap is a data-availability issue, not a quality issue --
c     it must never remove points by itself.
      call assert_i('F10 gap: does not remove points', nremoved, 0,
     &   npass, nfail)

c     ==============================================================
c     Real-data false-positive check: 42928 and 35497 (issue #10's
c     explicit request) -- well-tracked objects, expect zero or very
c     few flagged points.
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Real data: 42928 ---'
      call tle_evolve('input/example_42928.tle.txt', 42928, maxreal,
     &   epochs_r, a_r, e_r, inc_r, raan_r, aop_r,
     &   ha_r, hp_r, ls_r, npts_r, ierr_r)
      call assert_i('F11 42928 tle_evolve ierr=0', ierr_r, 0,
     &   npass, nfail)

      call tle_filter(epochs_r, ha_r, e_r, npts_r,
     &   4.0d0, 6.0d0, 30.0d0, 20,
     &   mask_r, nfiltered, nremoved, ngaps)
      write(*,'(A,I5,A,I5,A,I5)') '   42928: npts=', npts_r,
     &   ' nfiltered=', nfiltered, ' nremoved=', nremoved
c     Rate check, not an absolute count -- scales with dataset size.
c     Real orbital data is never perfectly linear, so some nonzero
c     background rate is expected (tuned to ~1.5% on the validation
c     set); this catches the check regressing back to the ~20%
c     false-positive rate an untuned 3-sigma/8-point window gave.
      if (nremoved * 100 .le. nfiltered * 5) then
         call assert_i('F12 42928: false-positive rate <=5%', 1, 1,
     &      npass, nfail)
      else
         call assert_i('F12 42928: false-positive rate <=5%', 0, 1,
     &      npass, nfail)
      end if

      write(*,*) ''
      write(*,*) '  --- Real data: 35497 ---'
      call tle_evolve('input/example_35497.tle.txt', 35497, maxreal,
     &   epochs_r, a_r, e_r, inc_r, raan_r, aop_r,
     &   ha_r, hp_r, ls_r, npts_r, ierr_r)
      call assert_i('F13 35497 tle_evolve ierr=0', ierr_r, 0,
     &   npass, nfail)

      call tle_filter(epochs_r, ha_r, e_r, npts_r,
     &   4.0d0, 6.0d0, 30.0d0, 20,
     &   mask_r, nfiltered, nremoved, ngaps)
      write(*,'(A,I5,A,I5,A,I5)') '   35497: npts=', npts_r,
     &   ' nfiltered=', nfiltered, ' nremoved=', nremoved
      if (nremoved * 100 .le. nfiltered * 5) then
         call assert_i('F14 35497: false-positive rate <=5%', 1, 1,
     &      npass, nfail)
      else
         call assert_i('F14 35497: false-positive rate <=5%', 0, 1,
     &      npass, nfail)
      end if

      write(*,*) ''
      write(*,*) '=============================================='
      write(*,'(A,I3,A,I3,A,I3)') ' Total: ', npass+nfail,
     &   '   Passed: ', npass, '   Failed: ', nfail
      write(*,*) '=============================================='
      if (nfail .eq. 0) then
         write(*,*) 'ALL TESTS PASSED'
      else
         write(*,*) 'SOME TESTS FAILED'
         call exit(1)
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

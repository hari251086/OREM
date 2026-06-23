c     test_zone_select.F — Tests for zone selection algorithm
c     -------------------------------------------------------
c     Compile:
c       ifx test_zone_select.F zone_select.F tle_evolution.F
c           ksrop/Subrouts.F ksrop/TLEread.F ksrop/Legendre.F
c           /exe:test_zone_select.exe
c     -------------------------------------------------------
      program test_zone_select

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(maxpts = 10000)

      double precision epochs(maxpts), ha(maxpts), hp(maxpts)
      double precision a_out(maxpts), e_out(maxpts)
      double precision inc_out(maxpts)
      double precision raan_out(maxpts), aop_out(maxpts)
      double precision ls_out(maxpts)
      integer npts, ierr

      parameter(mxz = 10)
      integer zone_start(mxz), zone_end(mxz), nzones

c     linfit test variables
      double precision xfit(5), yfit(5)
      double precision slope, intercept, r2
      double precision ls_min, ls_max

      integer npass, nfail, iall_ok, i, j

      call init_constants()

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' Zone Selection Test Suite'
      write(*,*) '=============================================='

c     ==============================================================
c     linfit unit tests
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- linfit unit tests ---'

c     Z1: Perfect linear data y = 2x + 3
      xfit(1) = 1.d0
      xfit(2) = 2.d0
      xfit(3) = 3.d0
      xfit(4) = 4.d0
      xfit(5) = 5.d0
      do i = 1, 5
         yfit(i) = 2.d0*xfit(i) + 3.d0
      end do
      call linfit(xfit, yfit, 5, slope, intercept, r2)
      call assert_d('Z1 linfit slope=2',
     &              slope, 2.d0, 1.d-10, npass, nfail)
      call assert_d('Z1 linfit intercept=3',
     &              intercept, 3.d0, 1.d-10, npass, nfail)
      call assert_d('Z1 linfit R2=1',
     &              r2, 1.d0, 1.d-10, npass, nfail)

c     Z2: Negative slope y = -5x + 100
      do i = 1, 5
         yfit(i) = -5.d0*xfit(i) + 100.d0
      end do
      call linfit(xfit, yfit, 5, slope, intercept, r2)
      call assert_d('Z2 linfit slope=-5',
     &              slope, -5.d0, 1.d-10, npass, nfail)
      call assert_d('Z2 linfit R2=1 (neg slope)',
     &              r2, 1.d0, 1.d-10, npass, nfail)

c     Z3: Constant data y = 42 → slope=0, R²=1
      do i = 1, 5
         yfit(i) = 42.d0
      end do
      call linfit(xfit, yfit, 5, slope, intercept, r2)
      call assert_d('Z3 linfit constant slope=0',
     &              slope, 0.d0, 1.d-10, npass, nfail)

c     Z4: Two points → perfect fit
      call linfit(xfit, yfit, 2, slope, intercept, r2)
      call assert_d('Z4 linfit 2pts R2=1',
     &              r2, 1.d0, 1.d-10, npass, nfail)

c     ==============================================================
c     Synthetic data tests
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Synthetic data ---'

c     Z5: Perfect linear decay → 1 zone
c     20 points, ha = 30000 - 100*t, t = 0..19 days
      do i = 1, 20
         epochs(i) = 2460000.d0 + dble(i-1)*0.5d0
         ha(i) = 30000.d0 - 200.d0*dble(i-1)*0.5d0
      end do
      call zone_select(epochs, ha, 20,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      if (nzones .ge. 1) then
         call assert_i('Z5 linear decay: >=1 zone',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z5 linear decay: >=1 zone',
     &                 0, 1, npass, nfail)
      end if

c     Z6: Flat apogee (no decay) → 0 zones
      do i = 1, 20
         epochs(i) = 2460000.d0 + dble(i-1)*0.5d0
         ha(i) = 30000.d0
      end do
      call zone_select(epochs, ha, 20,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      call assert_i('Z6 flat apogee: 0 zones',
     &              nzones, 0, npass, nfail)

c     Z7: Oscillating apogee (high-freq sinusoid) → 0 zones
c     High frequency + low amplitude so no segment is linear
      do i = 1, 40
         epochs(i) = 2460000.d0 + dble(i-1)*0.25d0
         ha(i) = 30000.d0 + 500.d0*dsin(dble(i-1)*2.d0)
      end do
      call zone_select(epochs, ha, 40,
     &                 8, 15.d0, 0.98d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      call assert_i('Z7 oscillating: 0 zones',
     &              nzones, 0, npass, nfail)

c     Z8: Rising apogee → 0 zones
      do i = 1, 20
         epochs(i) = 2460000.d0 + dble(i-1)*0.5d0
         ha(i) = 30000.d0 + 100.d0*dble(i-1)*0.5d0
      end do
      call zone_select(epochs, ha, 20,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      call assert_i('Z8 rising apogee: 0 zones',
     &              nzones, 0, npass, nfail)

c     Z9: Empty input → 0 zones
      call zone_select(epochs, ha, 0,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      call assert_i('Z9 empty: 0 zones',
     &              nzones, 0, npass, nfail)

c     Z10: Single point → 0 zones
      call zone_select(epochs, ha, 1,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      call assert_i('Z10 single pt: 0 zones',
     &              nzones, 0, npass, nfail)

c     Z11: nzones_max cap — long linear decay, request 1 zone
      do i = 1, 50
         epochs(i) = 2460000.d0 + dble(i-1)*0.3d0
         ha(i) = 30000.d0 - 50.d0*dble(i-1)*0.3d0
      end do
      call zone_select(epochs, ha, 50,
     &                 8, 15.d0, 0.90d0, -1.d0, 1,
     &                 zone_start, zone_end, nzones)
      call assert_i('Z11 nzones_max=1 caps',
     &              nzones, 1, npass, nfail)

c     ==============================================================
c     Real HEO data: 42928 (PSLV-C39, low-e, i=19°)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- 42928 PSLV-C39 ---'

      call tle_evolve('input/example_42928.tle.txt',
     &               42928, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)

      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)

      if (nzones .ge. 1) then
         call assert_i('Z12 42928 finds >=1 zone',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z12 42928 finds >=1 zone',
     &                 0, 1, npass, nfail)
      end if

c     Z13: Zone indices valid
      iall_ok = 1
      do i = 1, nzones
         if (zone_start(i) .lt. 1) iall_ok = 0
         if (zone_end(i) .gt. npts) iall_ok = 0
         if (zone_start(i) .gt. zone_end(i)) iall_ok = 0
      end do
      call assert_i('Z13 42928 indices valid',
     &              iall_ok, 1, npass, nfail)

c     Z14: Zones non-overlapping
      iall_ok = 1
      do i = 2, nzones
         if (zone_start(i) .le. zone_end(i-1)) iall_ok = 0
      end do
      call assert_i('Z14 42928 non-overlapping',
     &              iall_ok, 1, npass, nfail)

c     Z15: All zones have negative slope
      iall_ok = 1
      do i = 1, nzones
         nn = zone_end(i) - zone_start(i) + 1
         call linfit(epochs(zone_start(i)), ha(zone_start(i)),
     &               nn, slope, intercept, r2)
         if (slope .ge. 0.d0) iall_ok = 0
      end do
      call assert_i('Z15 42928 all slopes < 0',
     &              iall_ok, 1, npass, nfail)

c     Z16: All zones have R² > threshold
      iall_ok = 1
      do i = 1, nzones
         nn = zone_end(i) - zone_start(i) + 1
         call linfit(epochs(zone_start(i)), ha(zone_start(i)),
     &               nn, slope, intercept, r2)
         if (r2 .lt. 0.90d0) iall_ok = 0
      end do
      call assert_i('Z16 42928 all R2 > 0.90',
     &              iall_ok, 1, npass, nfail)

c     Z17: Minimum points per zone
      iall_ok = 1
      do i = 1, nzones
         if (zone_end(i)-zone_start(i)+1 .lt. 8) iall_ok = 0
      end do
      call assert_i('Z17 42928 all zones >=8 pts',
     &              iall_ok, 1, npass, nfail)

c     ==============================================================
c     Real HEO data: 35497 (Ariane 5 ESC-A, GTO, i=5.7°)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- 35497 Ariane 5 ESC-A ---'

      call tle_evolve('input/example_35497.tle.txt',
     &               35497, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)

      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)

      if (nzones .ge. 1) then
         call assert_i('Z18 35497 finds >=1 zone',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z18 35497 finds >=1 zone',
     &                 0, 1, npass, nfail)
      end if

c     Z19: Verify zone properties
      iall_ok = 1
      do i = 1, nzones
         if (zone_start(i) .lt. 1) iall_ok = 0
         if (zone_end(i) .gt. npts) iall_ok = 0
         nn = zone_end(i) - zone_start(i) + 1
         call linfit(epochs(zone_start(i)), ha(zone_start(i)),
     &               nn, slope, intercept, r2)
         if (slope .ge. 0.d0) iall_ok = 0
         if (r2 .lt. 0.90d0) iall_ok = 0
      end do
      call assert_i('Z19 35497 zones valid',
     &              iall_ok, 1, npass, nfail)

c     ==============================================================
c     Real HEO data: 37151 (Long March 3B, GTO, i=24.9°)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- 37151 Long March 3B ---'

      call tle_evolve('input/example_37151.tle.txt',
     &               37151, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)

      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)

      if (nzones .ge. 1) then
         call assert_i('Z20 37151 finds >=1 zone',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z20 37151 finds >=1 zone',
     &                 0, 1, npass, nfail)
      end if

c     Z21: Verify zone properties
      iall_ok = 1
      do i = 1, nzones
         if (zone_start(i) .lt. 1) iall_ok = 0
         if (zone_end(i) .gt. npts) iall_ok = 0
         nn = zone_end(i) - zone_start(i) + 1
         call linfit(epochs(zone_start(i)), ha(zone_start(i)),
     &               nn, slope, intercept, r2)
         if (slope .ge. 0.d0) iall_ok = 0
         if (r2 .lt. 0.90d0) iall_ok = 0
      end do
      call assert_i('Z21 37151 zones valid',
     &              iall_ok, 1, npass, nfail)

c     ==============================================================
c     Real HEO data: 39615 (Proton-M Briz-M, i=48.5°)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- 39615 Proton-M Briz-M ---'

      call tle_evolve('input/example_39615.tle.txt',
     &               39615, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)

      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)

      if (nzones .ge. 1) then
         call assert_i('Z22 39615 finds >=1 zone',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z22 39615 finds >=1 zone',
     &                 0, 1, npass, nfail)
      end if

c     Z23: Verify zone properties
      iall_ok = 1
      do i = 1, nzones
         if (zone_start(i) .lt. 1) iall_ok = 0
         if (zone_end(i) .gt. npts) iall_ok = 0
         nn = zone_end(i) - zone_start(i) + 1
         call linfit(epochs(zone_start(i)), ha(zone_start(i)),
     &               nn, slope, intercept, r2)
         if (slope .ge. 0.d0) iall_ok = 0
         if (r2 .lt. 0.90d0) iall_ok = 0
      end do
      call assert_i('Z23 39615 zones valid',
     &              iall_ok, 1, npass, nfail)

c     ==============================================================
c     Cross-case validation
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Cross-case validation ---'

c     Z24: 42928 zone-0 subset: should be 1 zone covering all
      call tle_evolve('input/example_42928_zone0.tle.txt',
     &               42928, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)
      call zone_select(epochs, ha, npts,
     &                 5, 10.d0, 0.85d0, -0.5d0, 4,
     &                 zone_start, zone_end, nzones)
      if (nzones .ge. 1) then
         call assert_i('Z24 42928-z0 >=1 zone',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z24 42928-z0 >=1 zone',
     &                 0, 1, npass, nfail)
      end if

c     Z25: Zone spans most of the data
      if (nzones .ge. 1) then
         nn = zone_end(1) - zone_start(1) + 1
         if (nn .ge. npts/2) then
            call assert_i('Z25 z0 spans >=50% data',
     &                    1, 1, npass, nfail)
         else
            call assert_i('Z25 z0 spans >=50% data',
     &                    0, 1, npass, nfail)
         end if
      else
         call assert_i('Z25 z0 spans >=50% data',
     &                 0, 1, npass, nfail)
      end if

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

c     ============================================================
      subroutine assert_d(label, got, expect, tol, npass, nfail)
      implicit double precision (a-h, o-z)
      character*(*) label
      double precision got, expect, tol
      integer npass, nfail
      if (dabs(got - expect) .le. tol) then
         write(*,'(A,A)') '  [PASS]  ', label
         npass = npass + 1
      else
         write(*,'(A,A)') '  [FAIL]  ', label
         write(*,'(A,E20.12,A,E20.12)') '           got=',
     &        got, '  expect=', expect
         nfail = nfail + 1
      end if
      return
      end

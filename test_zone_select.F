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
c     Extended linfit edge cases
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Extended linfit ---'

c     Z26: linfit with noisy data — R² < 1 but > 0
      xfit(1) = 1.d0
      xfit(2) = 2.d0
      xfit(3) = 3.d0
      xfit(4) = 4.d0
      xfit(5) = 5.d0
      yfit(1) = 10.2d0
      yfit(2) = 7.8d0
      yfit(3) = 6.1d0
      yfit(4) = 3.9d0
      yfit(5) = 2.3d0
      call linfit(xfit, yfit, 5, slope, intercept, r2)
      if (slope .lt. 0.d0 .and. r2 .gt. 0.95d0
     &    .and. r2 .lt. 1.d0) then
         call assert_i('Z26 noisy linfit R2 in (0.95,1)',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z26 noisy linfit R2 in (0.95,1)',
     &                 0, 1, npass, nfail)
      end if

c     Z27: linfit with 1 point — returns zeros
      call linfit(xfit, yfit, 1, slope, intercept, r2)
      call assert_d('Z27 linfit 1pt slope=0',
     &              slope, 0.d0, 1.d-10, npass, nfail)

c     ==============================================================
c     Parameter sensitivity
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Parameter sensitivity ---'

c     Use 42928 for parameter tests
      call tle_evolve('input/example_42928.tle.txt',
     &               42928, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)

c     Z28: Relaxed R² (0.80) finds more or equal zones
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 10,
     &                 zone_start, zone_end, nzones)
      nz_strict = nzones
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.80d0, -1.d0, 10,
     &                 zone_start, zone_end, nzones)
      if (nzones .ge. nz_strict) then
         call assert_i('Z28 relaxed R2 >= strict',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z28 relaxed R2 >= strict',
     &                 0, 1, npass, nfail)
      end if

c     Z29: Smaller min_zone_pts finds more or equal zones
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 10,
     &                 zone_start, zone_end, nzones)
      nz_strict = nzones
      call zone_select(epochs, ha, npts,
     &                 5, 10.d0, 0.90d0, -1.d0, 10,
     &                 zone_start, zone_end, nzones)
      if (nzones .ge. nz_strict) then
         call assert_i('Z29 smaller minpts >= strict',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z29 smaller minpts >= strict',
     &                 0, 1, npass, nfail)
      end if

c     Z30: min_zone_pts > npts → 0 zones
      call zone_select(epochs, ha, npts,
     &                 npts+1, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      call assert_i('Z30 minpts>npts: 0 zones',
     &              nzones, 0, npass, nfail)

c     Z31: Very short max_zone_days (2 days)
      call zone_select(epochs, ha, npts,
     &                 5, 2.d0, 0.85d0, -0.5d0, 10,
     &                 zone_start, zone_end, nzones)
c     All zones should be <= 2 days
      iall_ok = 1
      do i = 1, nzones
         span_days = epochs(zone_end(i)) - epochs(zone_start(i))
         if (span_days .gt. 2.1d0) iall_ok = 0
      end do
      call assert_i('Z31 max_days=2 enforced',
     &              iall_ok, 1, npass, nfail)

c     Z32: Lenient slope threshold (-0.1) finds more zones
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -10.d0, 10,
     &                 zone_start, zone_end, nzones)
      nz_strict = nzones
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -0.1d0, 10,
     &                 zone_start, zone_end, nzones)
      if (nzones .ge. nz_strict) then
         call assert_i('Z32 lenient slope >= strict',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z32 lenient slope >= strict',
     &                 0, 1, npass, nfail)
      end if

c     ==============================================================
c     Advanced synthetic patterns
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Advanced synthetic ---'

c     Z33: Two-segment decay with flat middle → finds 2 zones
      do i = 1, 15
         epochs(i) = 2460000.d0 + dble(i-1)*0.5d0
         ha(i) = 30000.d0 - 100.d0*dble(i-1)*0.5d0
      end do
      do i = 16, 25
         epochs(i) = epochs(15) + dble(i-15)*0.5d0
         ha(i) = ha(15)
      end do
      do i = 26, 40
         epochs(i) = epochs(25) + dble(i-25)*0.5d0
         ha(i) = ha(25) - 80.d0*dble(i-25)*0.5d0
      end do
      call zone_select(epochs, ha, 40,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      if (nzones .ge. 2) then
         call assert_i('Z33 two-segment: >=2 zones',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z33 two-segment: >=2 zones',
     &                 0, 1, npass, nfail)
      end if

c     Z34: Noisy linear decay (noise ±5 km on 100 km/day slope)
c     Should still find zone if R² threshold is 0.90
      do i = 1, 30
         epochs(i) = 2460000.d0 + dble(i-1)*0.3d0
         ha(i) = 30000.d0 - 100.d0*dble(i-1)*0.3d0
     &         + 5.d0*dsin(dble(i)*7.d0)
      end do
      call zone_select(epochs, ha, 30,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      if (nzones .ge. 1) then
         call assert_i('Z34 noisy linear: >=1 zone',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z34 noisy linear: >=1 zone',
     &                 0, 1, npass, nfail)
      end if

c     Z35: Step function (abrupt drop then linear)
      do i = 1, 10
         epochs(i) = 2460000.d0 + dble(i-1)*0.5d0
         ha(i) = 30000.d0
      end do
      do i = 11, 30
         epochs(i) = epochs(10) + dble(i-10)*0.5d0
         ha(i) = 28000.d0 - 50.d0*dble(i-10)*0.5d0
      end do
      call zone_select(epochs, ha, 30,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      if (nzones .ge. 1) then
c        Zone should be in the post-step region
         if (zone_start(1) .ge. 10) then
            call assert_i('Z35 step: zone after step',
     &                    1, 1, npass, nfail)
         else
            call assert_i('Z35 step: zone after step',
     &                    0, 1, npass, nfail)
         end if
      else
         call assert_i('Z35 step: zone after step',
     &                 0, 1, npass, nfail)
      end if

c     Z36: Steep decay (-500 km/day) vs gradual (-5 km/day)
c     Both should be found with slope_threshold = -1
      do i = 1, 15
         epochs(i) = 2460000.d0 + dble(i-1)*0.3d0
         ha(i) = 30000.d0 - 500.d0*dble(i-1)*0.3d0
      end do
      call zone_select(epochs, ha, 15,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      nz1 = nzones
      do i = 1, 15
         ha(i) = 30000.d0 - 5.d0*dble(i-1)*0.3d0
      end do
      call zone_select(epochs, ha, 15,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      if (nz1 .ge. 1 .and. nzones .ge. 1) then
         call assert_i('Z36 steep+gradual both found',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z36 steep+gradual both found',
     &                 0, 1, npass, nfail)
      end if

c     ==============================================================
c     Real data deep validation
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Real data deep validation ---'

c     Reload all 4 objects and verify detailed zone properties

c     Z37-Z39: 42928 zone details
      call tle_evolve('input/example_42928.tle.txt',
     &               42928, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)

c     Z37: All zone durations within max_zone_days
      iall_ok = 1
      do i = 1, nzones
         span_days = epochs(zone_end(i)) - epochs(zone_start(i))
         if (span_days .gt. 10.1d0) iall_ok = 0
      end do
      call assert_i('Z37 42928 durations <=10d',
     &              iall_ok, 1, npass, nfail)

c     Z38: Zones sorted by epoch (ascending)
      iall_ok = 1
      do i = 2, nzones
         if (epochs(zone_start(i)) .le.
     &       epochs(zone_start(i-1))) iall_ok = 0
      end do
      call assert_i('Z38 42928 sorted by epoch',
     &              iall_ok, 1, npass, nfail)

c     Z39: Zone apogee altitudes are physical (> 0 km)
      iall_ok = 1
      do i = 1, nzones
         do j = zone_start(i), zone_end(i)
            if (ha(j) .le. 0.d0) iall_ok = 0
         end do
      end do
      call assert_i('Z39 42928 zone ha > 0',
     &              iall_ok, 1, npass, nfail)

c     Z40-Z42: 35497 zone details
      call tle_evolve('input/example_35497.tle.txt',
     &               35497, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)

c     Z40: Durations valid (positive and reasonable)
      iall_ok = 1
      do i = 1, nzones
         span_days = epochs(zone_end(i)) - epochs(zone_start(i))
         if (span_days .gt. 15.d0) iall_ok = 0
         if (span_days .lt. 0.d0) iall_ok = 0
      end do
      call assert_i('Z40 35497 durations valid',
     &              iall_ok, 1, npass, nfail)

c     Z41: Sorted and non-overlapping
      iall_ok = 1
      do i = 2, nzones
         if (zone_start(i) .le. zone_end(i-1)) iall_ok = 0
      end do
      call assert_i('Z41 35497 sorted+non-overlap',
     &              iall_ok, 1, npass, nfail)

c     Z42: Zone ha > 0
      iall_ok = 1
      do i = 1, nzones
         do j = zone_start(i), zone_end(i)
            if (ha(j) .le. 0.d0) iall_ok = 0
         end do
      end do
      call assert_i('Z42 35497 zone ha > 0',
     &              iall_ok, 1, npass, nfail)

c     ==============================================================
c     Repeatability and robustness
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Repeatability ---'

c     Z43: Two calls with same input → identical zones
      call tle_evolve('input/example_42928.tle.txt',
     &               42928, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      nz1 = nzones
      ns1 = zone_start(1)
      ne1 = zone_end(1)
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      iall_ok = 1
      if (nzones .ne. nz1) iall_ok = 0
      if (nzones .ge. 1) then
         if (zone_start(1) .ne. ns1) iall_ok = 0
         if (zone_end(1) .ne. ne1) iall_ok = 0
      end if
      call assert_i('Z43 repeatability identical',
     &              iall_ok, 1, npass, nfail)

c     Z44: Large nzones_max (100) doesn't crash
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, mxz,
     &                 zone_start, zone_end, nzones)
      if (nzones .ge. 0 .and. nzones .le. mxz) then
         call assert_i('Z44 large nzmax safe',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z44 large nzmax safe',
     &                 0, 1, npass, nfail)
      end if

c     Z45: nzones_max=0 → 0 zones (edge)
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 0,
     &                 zone_start, zone_end, nzones)
      call assert_i('Z45 nzmax=0: 0 zones',
     &              nzones, 0, npass, nfail)

c     ==============================================================
c     37151/39615 deep validation (gap fill)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- 37151/39615 deep validation ---'

c     Z46-Z48: 37151 Long March 3B deep checks
      call tle_evolve('input/example_37151.tle.txt',
     &               37151, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)

c     Z46: Durations valid
      iall_ok = 1
      do i = 1, nzones
         span_days = epochs(zone_end(i)) - epochs(zone_start(i))
         if (span_days .gt. 15.d0 .or. span_days .lt. 0.d0)
     &      iall_ok = 0
      end do
      call assert_i('Z46 37151 durations valid',
     &              iall_ok, 1, npass, nfail)

c     Z47: Sorted and non-overlapping
      iall_ok = 1
      do i = 2, nzones
         if (zone_start(i) .le. zone_end(i-1)) iall_ok = 0
      end do
      call assert_i('Z47 37151 sorted+non-overlap',
     &              iall_ok, 1, npass, nfail)

c     Z48: Zone ha > 0
      iall_ok = 1
      do i = 1, nzones
         do j = zone_start(i), zone_end(i)
            if (ha(j) .le. 0.d0) iall_ok = 0
         end do
      end do
      call assert_i('Z48 37151 zone ha > 0',
     &              iall_ok, 1, npass, nfail)

c     Z49-Z51: 39615 Proton-M deep checks
      call tle_evolve('input/example_39615.tle.txt',
     &               39615, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)

c     Z49: Durations valid
      iall_ok = 1
      do i = 1, nzones
         span_days = epochs(zone_end(i)) - epochs(zone_start(i))
         if (span_days .gt. 15.d0 .or. span_days .lt. 0.d0)
     &      iall_ok = 0
      end do
      call assert_i('Z49 39615 durations valid',
     &              iall_ok, 1, npass, nfail)

c     Z50: Sorted and non-overlapping
      iall_ok = 1
      do i = 2, nzones
         if (zone_start(i) .le. zone_end(i-1)) iall_ok = 0
      end do
      call assert_i('Z50 39615 sorted+non-overlap',
     &              iall_ok, 1, npass, nfail)

c     Z51: Zone ha > 0
      iall_ok = 1
      do i = 1, nzones
         do j = zone_start(i), zone_end(i)
            if (ha(j) .le. 0.d0) iall_ok = 0
         end do
      end do
      call assert_i('Z51 39615 zone ha > 0',
     &              iall_ok, 1, npass, nfail)

c     ==============================================================
c     Ha monotonicity within zones
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Ha monotonicity in zones ---'

c     Check all 4 objects: no local ha spike > 5% within any zone
c     (allows small TLE noise but catches big excursions)

c     Z52: 42928 no local ha spike in zones
      call tle_evolve('input/example_42928.tle.txt',
     &               42928, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      iall_ok = 1
      do i = 1, nzones
         ha_z1 = ha(zone_start(i))
         ha_zn = ha(zone_end(i))
         ha_range = dabs(ha_z1 - ha_zn)
         if (ha_range .lt. 1.d0) ha_range = 1.d0
         do j = zone_start(i)+1, zone_end(i)
            ha_jump = ha(j) - ha(j-1)
            if (ha_jump .gt. 0.05d0*ha_range) iall_ok = 0
         end do
      end do
      call assert_i('Z52 42928 no ha spike in zone',
     &              iall_ok, 1, npass, nfail)

c     Z53: 37151 no local ha spike
      call tle_evolve('input/example_37151.tle.txt',
     &               37151, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      iall_ok = 1
      do i = 1, nzones
         ha_z1 = ha(zone_start(i))
         ha_zn = ha(zone_end(i))
         ha_range = dabs(ha_z1 - ha_zn)
         if (ha_range .lt. 1.d0) ha_range = 1.d0
         do j = zone_start(i)+1, zone_end(i)
            ha_jump = ha(j) - ha(j-1)
            if (ha_jump .gt. 0.05d0*ha_range) iall_ok = 0
         end do
      end do
      call assert_i('Z53 37151 no ha spike in zone',
     &              iall_ok, 1, npass, nfail)

c     ==============================================================
c     max_zone_days enforcement on all objects
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- max_zone_days all objects ---'

c     Z54: 37151 max_zone_days=10 enforced
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      iall_ok = 1
      do i = 1, nzones
         span_days = epochs(zone_end(i)) - epochs(zone_start(i))
         if (span_days .gt. 10.1d0) iall_ok = 0
      end do
      call assert_i('Z54 37151 durations <=10d',
     &              iall_ok, 1, npass, nfail)

c     Z55: 39615 max_zone_days=10 enforced
      call tle_evolve('input/example_39615.tle.txt',
     &               39615, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      iall_ok = 1
      do i = 1, nzones
         span_days = epochs(zone_end(i)) - epochs(zone_start(i))
         if (span_days .gt. 10.1d0) iall_ok = 0
      end do
      call assert_i('Z55 39615 durations <=10d',
     &              iall_ok, 1, npass, nfail)

c     ==============================================================
c     Zone count reasonableness
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Zone count bounds ---'

c     Z56: 42928 (2+ years) should have 1-10 zones, not 50+
      call tle_evolve('input/example_42928.tle.txt',
     &               42928, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, mxz,
     &                 zone_start, zone_end, nzones)
      if (nzones .ge. 1 .and. nzones .le. 10) then
         call assert_i('Z56 42928 count in [1,10]',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z56 42928 count in [1,10]',
     &                 0, 1, npass, nfail)
      end if

c     Z57: 35497 (3+ years) zone count reasonable
      call tle_evolve('input/example_35497.tle.txt',
     &               35497, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha, hp, ls_out,
     &               npts, ierr)
      call zone_select(epochs, ha, npts,
     &                 8, 10.d0, 0.90d0, -1.d0, mxz,
     &                 zone_start, zone_end, nzones)
      if (nzones .ge. 1 .and. nzones .le. 10) then
         call assert_i('Z57 35497 count in [1,10]',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z57 35497 count in [1,10]',
     &                 0, 1, npass, nfail)
      end if

c     ==============================================================
c     Exact min_zone_pts boundary
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Boundary: npts = min_zone_pts ---'

c     Z58: npts exactly = min_zone_pts, linear decay → 1 zone
      do i = 1, 8
         epochs(i) = 2460000.d0 + dble(i-1)*1.d0
         ha(i) = 30000.d0 - 200.d0*dble(i-1)
      end do
      call zone_select(epochs, ha, 8,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      if (nzones .eq. 1) then
         call assert_i('Z58 npts=minpts: 1 zone',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z58 npts=minpts: 1 zone',
     &                 0, 1, npass, nfail)
      end if

c     Z59: npts = min_zone_pts - 1 → 0 zones
      call zone_select(epochs, ha, 7,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      call assert_i('Z59 npts=minpts-1: 0 zones',
     &              nzones, 0, npass, nfail)

c     ==============================================================
c     Sparse data (large epoch gaps)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Sparse data ---'

c     Z60: Linear decay with 30-day gaps between points
c     10 points spanning 270 days but max_zone_days=15
c     Should find 0 zones (gaps too large for single zone)
      do i = 1, 10
         epochs(i) = 2460000.d0 + dble(i-1)*30.d0
         ha(i) = 30000.d0 - 100.d0*dble(i-1)*30.d0
      end do
      call zone_select(epochs, ha, 10,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      call assert_i('Z60 sparse 30d gaps: 0 zones',
     &              nzones, 0, npass, nfail)

c     Z61: Dense cluster then large gap then dense cluster
c     Each cluster has 10 pts in 5 days, gap = 60 days
      do i = 1, 10
         epochs(i) = 2460000.d0 + dble(i-1)*0.5d0
         ha(i) = 30000.d0 - 100.d0*dble(i-1)*0.5d0
      end do
      do i = 11, 20
         epochs(i) = 2460060.d0 + dble(i-11)*0.5d0
         ha(i) = 25000.d0 - 80.d0*dble(i-11)*0.5d0
      end do
      call zone_select(epochs, ha, 20,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      if (nzones .ge. 2) then
         call assert_i('Z61 two clusters: >=2 zones',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z61 two clusters: >=2 zones',
     &                 0, 1, npass, nfail)
      end if

c     ==============================================================
c     Independent R² verification
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Independent R2 check ---'

c     Z62: Verify R² by computing SS_res/SS_tot independently
c     Use a known noisy dataset
      xfit(1) = 1.d0
      xfit(2) = 2.d0
      xfit(3) = 3.d0
      xfit(4) = 4.d0
      xfit(5) = 5.d0
      yfit(1) = 10.d0
      yfit(2) = 8.1d0
      yfit(3) = 5.9d0
      yfit(4) = 4.2d0
      yfit(5) = 1.8d0
      call linfit(xfit, yfit, 5, slope, intercept, r2)
c     Compute R² independently:
c     ymean = (10+8.1+5.9+4.2+1.8)/5 = 6.0
c     SS_tot = sum((yi - 6)²)
c            = 16 + 4.41 + 0.01 + 3.24 + 17.64 = 41.30
c     predicted: slope*x + intercept for each
c     SS_res = sum((yi - pred)²)
      ymean2 = 6.d0
      ss_tot2 = 0.d0
      ss_res2 = 0.d0
      do i = 1, 5
         ss_tot2 = ss_tot2 + (yfit(i) - ymean2)**2
         pred = slope*xfit(i) + intercept
         ss_res2 = ss_res2 + (yfit(i) - pred)**2
      end do
      r2_indep = 1.d0 - ss_res2 / ss_tot2
      call assert_d('Z62 R2 indep verify',
     &              r2, r2_indep, 1.d-10, npass, nfail)

c     ==============================================================
c     Degenerate cases
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Degenerate cases ---'

c     Z63: All identical epochs → 0 zones (degenerate)
      do i = 1, 10
         epochs(i) = 2460000.d0
         ha(i) = 30000.d0 - dble(i)*100.d0
      end do
      call zone_select(epochs, ha, 10,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      call assert_i('Z63 identical epochs: 0 zones',
     &              nzones, 0, npass, nfail)

c     Z64: 2 points only (< default min_zone_pts) → 0 zones
      epochs(1) = 2460000.d0
      epochs(2) = 2460001.d0
      ha(1) = 30000.d0
      ha(2) = 29500.d0
      call zone_select(epochs, ha, 2,
     &                 8, 15.d0, 0.90d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      call assert_i('Z64 2pts < minpts: 0 zones',
     &              nzones, 0, npass, nfail)

c     Z65: Very steep decay (-10000 km/day)
      do i = 1, 10
         epochs(i) = 2460000.d0 + dble(i-1)*0.1d0
         ha(i) = 30000.d0 - 10000.d0*dble(i-1)*0.1d0
      end do
      call zone_select(epochs, ha, 10,
     &                 8, 2.d0, 0.85d0, -1.d0, 4,
     &                 zone_start, zone_end, nzones)
      if (nzones .ge. 1) then
         call assert_i('Z65 very steep: >=1 zone',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Z65 very steep: >=1 zone',
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

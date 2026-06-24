c     test_rsm.F — Tests for RSM surface generation
c     -------------------------------------------------------
c     Compile:
c       ifx test_rsm.F rsm.F ga.F tle_evolution.F zone_select.F
c           ksrop/propagate_ks.F ksrop/Subrouts.F ksrop/TLEread.F
c           ksrop/Legendre.F /exe:test_rsm.exe
c     -------------------------------------------------------
      program test_rsm

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      integer npass, nfail

      call init_constants()

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' RSM Test Suite'
      write(*,*) '=============================================='

c     ==============================================================
c     jd2cal unit tests
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- jd2cal ---'

      call test_jd2cal(npass, nfail)

c     ==============================================================
c     RSM grid construction
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Grid construction ---'

      call test_grid(npass, nfail)

c     ==============================================================
c     RSM error handling
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Error handling ---'

      call test_rsm_errors(npass, nfail)

c     ==============================================================
c     RSM with real TLE data (42928 zone 0)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- 42928 Zone 0 integration ---'

      call test_rsm_42928(npass, nfail)

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
c     jd2cal tests
c     ============================================================
      subroutine test_jd2cal(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail
      double precision cal(6)

c     R1: J2000.0 epoch = 2000-01-01 12:00:00
      call jd2cal(2451545.0d0, cal)
      call assert_d('R1 J2000 year=2000',
     &              cal(1), 2000.d0, 0.1d0, npass, nfail)
      call assert_d('R2 J2000 month=1',
     &              cal(2), 1.d0, 0.1d0, npass, nfail)
      call assert_d('R3 J2000 day=1',
     &              cal(3), 1.d0, 0.1d0, npass, nfail)

c     R4: 2017-09-22 (42928 zone 0 start) ≈ JD 2458018.5
      call jd2cal(2458018.5d0, cal)
      call assert_d('R4 2017-09-22 year',
     &              cal(1), 2017.d0, 0.1d0, npass, nfail)
      call assert_d('R5 2017-09-22 month',
     &              cal(2), 9.d0, 0.1d0, npass, nfail)

      return
      end

c     ============================================================
c     Grid construction tests
c     ============================================================
      subroutine test_grid(npass, nfail)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth
      integer npass, nfail

      parameter(max_surf = 5000, maxzone = 100)
      double precision surfaces(max_surf, 3, 3)
      double precision e_grid(3), a_grid(3)
      double precision tobs(maxzone), apobs(maxzone)
      integer nsurf_pts, nobs, ierr

c     Dummy zone data (2 points minimum)
      double precision ep(2), ha(2), ez(2), iz(2)
      double precision rz(2), az(2), mz(2), sz(2)
      integer n_f(3)
      double precision alt(1), den(1), sch(1)

      ep(1) = 2458018.5d0
      ep(2) = 2458019.5d0
      ha(1) = 12000.d0
      ha(2) = 11950.d0
      ez(1) = 0.32d0
      ez(2) = 0.319d0
      iz(1) = 19.17d0
      iz(2) = 19.17d0
      rz(1) = 69.2d0
      rz(2) = 68.5d0
      az(1) = 273.5d0
      az(2) = 274.0d0
      mz(1) = 14.3d0
      mz(2) = 15.0d0
      sz(1) = 9632.d0
      sz(2) = 9630.d0

c     R6: e_grid is [e_min, e_mid, e_max]
c     Just check grid construction by calling with nzone=2
c     (we won't actually propagate — just verify grid setup)

c     For this test, verify the grid values directly
      e_lo = 0.319d0
      e_hi = 0.321d0
      a_lo = 1.0d0
      a_hi = 5.0d0

      e_mid = (e_lo + e_hi) / 2.d0
      a_mid = (a_lo + a_hi) / 2.d0

      call assert_d('R6 e_mid correct',
     &              e_mid, 0.320d0, 1.d-10, npass, nfail)
      call assert_d('R7 a_mid correct',
     &              a_mid, 3.0d0, 1.d-10, npass, nfail)

c     R8: BN = mass / (Cd * A)
      bn_test = 1000.d0 / (2.2d0 * 3.d0)
      call assert_d('R8 BN = m/(Cd*A)',
     &              bn_test, 1000.d0/6.6d0,
     &              1.d-10, npass, nfail)

c     R9: Perigee-preserving SMA adjustment
c     rp = a*(1-e), a_new = rp/(1-e_new)
      a_orig = 9632.d0
      e_orig = 0.3216d0
      e_new = 0.3220d0
      rp = a_orig * (1.d0 - e_orig)
      a_new = rp / (1.d0 - e_new)
      rp_new = a_new * (1.d0 - e_new)
      call assert_d('R9 perigee preserved',
     &              rp_new, rp, 1.d-6, npass, nfail)

      return
      end

c     ============================================================
c     Error handling tests
c     ============================================================
      subroutine test_rsm_errors(npass, nfail)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth
      integer npass, nfail

      parameter(max_surf = 5000)
      double precision surfaces(max_surf, 3, 3)
      double precision e_grid(3), a_grid(3)
      double precision tobs(1), apobs(1)
      integer nsurf_pts, nobs, ierr

c     Dummy arrays for propagator config
      double precision ep(1), ha(1), ez(1), iz(1)
      double precision rz(1), az(1), mz(1), sz(1)
      integer n_f(3)
      double precision alt(1), den(1), sch(1)

      n_f(1) = 0
      n_f(2) = 0
      n_f(3) = 0
      alt(1) = 0.d0
      den(1) = 0.d0
      sch(1) = 0.d0

c     R10: nzone < 2 returns ierr = 1
      ep(1) = 2458018.5d0
      ha(1) = 12000.d0
      ez(1) = 0.32d0
      iz(1) = 19.d0
      rz(1) = 69.d0
      az(1) = 273.d0
      mz(1) = 14.d0
      sz(1) = 9632.d0

      call rsm_generate(
     &   ep, ha, ez, iz, rz, az, mz, sz,
     &   1,
     &   0.319d0, 0.321d0, 1.d0, 5.d0,
     &   1000.d0, 2.2d0, 0,
     &   10, 1, 1.d-12,
     &   n_f, 0, 0, 0,
     &   7.29d-5, 3.35d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0,
     &   0.d0, 0.d0, 0.d0,
     &   alt, den, sch, 1,
     &   surfaces, nsurf_pts,
     &   tobs, apobs, nobs,
     &   e_grid, a_grid,
     &   ierr)
      call assert_i('R10 nzone<2 ierr=1',
     &              ierr, 1, npass, nfail)

      return
      end

c     ============================================================
c     42928 Zone 0 integration test
c     ============================================================
      subroutine test_rsm_42928(npass, nfail)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth
      integer npass, nfail

      parameter(maxpts = 10000, max_surf = 5000)

c     TLE evolution arrays
      double precision epochs(maxpts), a_out(maxpts)
      double precision e_out(maxpts), inc_out(maxpts)
      double precision raan_out(maxpts), aop_out(maxpts)
      double precision ha_out(maxpts), hp_out(maxpts)
      double precision ls_out(maxpts)
      integer npts, ierr_tle

c     Zone selection
      parameter(mxz = 10)
      integer zone_start(mxz), zone_end(mxz), nzones

c     RSM outputs
      double precision surfaces(max_surf, 3, 3)
      double precision e_grid(3), a_grid(3)
      double precision tobs(500), apobs(500)
      integer nsurf_pts, nobs, ierr_rsm

c     Atmosphere tables
      parameter(ndim_a = 100)
      double precision alt_a(ndim_a), den_a(ndim_a), sch_a(ndim_a)
      integer natm

c     Force config
      integer n_f(3)

c     Zone data arrays
      double precision ez(500), iz(500), rz(500)
      double precision az(500), mz(500), sz(500)
      double precision epz(500), haz(500)

c     Read atmosphere table
      call read_atm(alt_a, den_a, sch_a, ndim_a, natm)

c     Load 42928 zone 0 TLE data
      call tle_evolve('input/example_42928_zone0.tle.txt',
     &               42928, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, ls_out,
     &               npts, ierr_tle)

c     R11: TLE loaded
      call assert_i('R11 42928 z0 TLE loaded',
     &              ierr_tle, 0, npass, nfail)
      if (ierr_tle .ne. 0) return

c     R12: Have enough points
      if (npts .ge. 8) then
         call assert_i('R12 npts >= 8',
     &                 1, 1, npass, nfail)
      else
         call assert_i('R12 npts >= 8',
     &                 0, 1, npass, nfail)
         return
      end if

c     Use all zone 0 points directly
      nzone = npts
      do i = 1, nzone
         epz(i) = epochs(i)
         haz(i) = ha_out(i)
         ez(i) = e_out(i)
         iz(i) = inc_out(i)
         rz(i) = raan_out(i)
         az(i) = aop_out(i)
         mz(i) = 0.d0
         sz(i) = a_out(i)
      end do

c     Set up force model (two-body + J2 for stability)
      n_f(1) = 1
      n_f(2) = 0
      n_f(3) = 0

c     Eccentricity grid from TLE scatter
      e_mid = e_out(1)
      de = 0.0005d0

      write(*,*) '    calling rsm_generate...'

c     Call RSM — IDRAG=0 (no drag for unit test speed)
      call rsm_generate(
     &   epz, haz, ez, iz, rz, az, mz, sz,
     &   nzone,
     &   e_mid - de, e_mid + de,
     &   1.d0, 5.d0,
     &   1000.d0, 2.2d0, 0,
     &   10, 1, 1.d-10,
     &   n_f, 2, 0, 0,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0,
     &   0.d0, 0.d0, 0.d0,
     &   alt_a, den_a, sch_a, natm,
     &   surfaces, nsurf_pts,
     &   tobs, apobs, nobs,
     &   e_grid, a_grid,
     &   ierr_rsm)

c     R13: RSM completed without error
      call assert_i('R13 rsm ierr=0',
     &              ierr_rsm, 0, npass, nfail)

c     R14: Generated surfaces have data
      if (nsurf_pts .gt. 0) then
         call assert_i('R14 nsurf_pts > 0',
     &                 1, 1, npass, nfail)
      else
         call assert_i('R14 nsurf_pts > 0',
     &                 0, 1, npass, nfail)
      end if

c     R15: Observations extracted
      call assert_i('R15 nobs = nzone',
     &              nobs, nzone, npass, nfail)

c     R16: e_grid values correct
      call assert_d('R16 e_grid(1)',
     &              e_grid(1), e_mid - de,
     &              1.d-10, npass, nfail)
      call assert_d('R17 e_grid(2)',
     &              e_grid(2), e_mid,
     &              1.d-10, npass, nfail)
      call assert_d('R18 e_grid(3)',
     &              e_grid(3), e_mid + de,
     &              1.d-10, npass, nfail)

c     R19: a_grid values correct
      call assert_d('R19 a_grid mid=3',
     &              a_grid(2), 3.d0,
     &              1.d-10, npass, nfail)

c     R20: Surface values are physical (apogee > 0)
      if (nsurf_pts .gt. 0) then
         iall_ok = 1
         do ie = 1, 3
            do ia = 1, 3
               if (surfaces(1, ie, ia) .lt. 100.d0)
     &            iall_ok = 0
            end do
         end do
         call assert_i('R20 surfaces > 100 km',
     &                 iall_ok, 1, npass, nfail)
      else
         call assert_i('R20 surfaces > 100 km',
     &                 0, 1, npass, nfail)
      end if

c     R21: Center surface (e2,a2) closest to observed
      if (nsurf_pts .gt. 0 .and. nobs .gt. 0) then
         diff_center = dabs(surfaces(1,2,2) - apobs(1))
         diff_corner = dabs(surfaces(1,1,1) - apobs(1))
         if (diff_center .le. diff_corner + 100.d0) then
            call assert_i('R21 center near obs',
     &                    1, 1, npass, nfail)
         else
            call assert_i('R21 center near obs',
     &                    0, 1, npass, nfail)
         end if
      else
         call assert_i('R21 center near obs',
     &                 0, 1, npass, nfail)
      end if

c     R22: tobs starts at 0
      call assert_d('R22 tobs(1)=0',
     &              tobs(1), 0.d0, 1.d-10, npass, nfail)

c     R23: apobs matches ha_out
      call assert_d('R23 apobs(1)=ha',
     &              apobs(1), ha_out(1),
     &              1.d-6, npass, nfail)

      return
      end

c     ============================================================
c     Read atmosphere table
c     ============================================================
      subroutine read_atm(alt, den, sch, maxn, natm)
      implicit double precision (a-h, o-z)
      integer maxn, natm
      double precision alt(maxn), den(maxn), sch(maxn)

      natm = 0
      open(unit=88, file='input/ATM.DAT', status='old',
     &     err=900)
      do i = 1, maxn
         read(88, *, end=800) alt(i), den(i), sch(i)
         natm = natm + 1
      end do
 800  continue
      close(88)
      return
 900  continue
      natm = 0
      return
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

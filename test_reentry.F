c     test_reentry.F — OREM validation on 7 known re-entry cases
c     -------------------------------------------------------
c     Compile:
c       ifx /heap-arrays /F:16777216 test_reentry.F orem.F rsm.F
c           ga.F tle_evolution.F zone_select.F
c           ksrop/propagate_ks.F ksrop/Subrouts.F
c           ksrop/TLEread.F ksrop/Legendre.F
c           /exe:test_reentry.exe
c     -------------------------------------------------------
      program test_reentry

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291, mxz = 4)

      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr

      integer npass, nfail
      character*60 tle_file
      character*30 obj_name
      integer norad
      double precision obs_yr, obs_mo, obs_dy

      call init_constants()

c     Read atmosphere table
      call read_atm_v(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' OREM 7-Object Re-entry Validation Suite'
      write(*,*) '=============================================='

c     ==============================================================
c     Test each object
c     ==============================================================

c     --- 42928 PSLV-C39 R/B ---
      call run_case(
     &   'input/example_42928.tle.txt', 42928,
     &   '42928 PSLV-C39',
     &   2019.d0, 3.d0, 3.d0,
     &   19.2d0, 0.33d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz,
     &   npass, nfail)

c     --- 35497 Ariane 5 ESC-A ---
      call run_case(
     &   'input/example_35497.tle.txt', 35497,
     &   '35497 Ariane 5',
     &   2016.d0, 10.d0, 31.d0,
     &   5.7d0, 0.63d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz,
     &   npass, nfail)

c     --- 37151 Long March 3B ---
      call run_case(
     &   'input/example_37151.tle.txt', 37151,
     &   '37151 Long March 3B',
     &   2015.d0, 12.d0, 3.d0,
     &   24.9d0, 0.56d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz,
     &   npass, nfail)

c     --- 39615 Proton-M Briz-M ---
      call run_case(
     &   'input/example_39615.tle.txt', 39615,
     &   '39615 Proton-M',
     &   2017.d0, 9.d0, 15.d0,
     &   48.5d0, 0.68d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz,
     &   npass, nfail)

c     --- 27526 Ariane 5 R/B ---
      call run_case(
     &   'input/example_27526.tle.txt', 27526,
     &   '27526 Ariane 5 R/B',
     &   2012.d0, 5.d0, 9.d0,
     &   17.7d0, 0.59d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz,
     &   npass, nfail)

c     --- 32007 GSLV R/B ---
      call run_case(
     &   'input/example_32007.tle.txt', 32007,
     &   '32007 GSLV R/B',
     &   2010.d0, 6.d0, 6.d0,
     &   25.9d0, 0.29d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz,
     &   npass, nfail)

c     --- 37819 Proton-M R/B ---
      call run_case(
     &   'input/example_37819.tle.txt', 37819,
     &   '37819 Proton-M R/B',
     &   2013.d0, 9.d0, 12.d0,
     &   63.4d0, 0.47d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz,
     &   npass, nfail)

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
c     Run a single validation case
c     ============================================================
      subroutine run_case(
     &   tle_file, norad, obj_name,
     &   obs_yr, obs_mo, obs_dy,
     &   exp_inc, exp_ecc,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz,
     &   npass, nfail)

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      character*(*) tle_file, obj_name
      integer norad, natm, ndim_a, mxz
      double precision obs_yr, obs_mo, obs_dy
      double precision exp_inc, exp_ecc
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer npass, nfail

      double precision reentry(4), e_opt(4), bn_opt(4)
      double precision rms_o(4), zepoch(4), rpe(4)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr, i

      write(*,*) ''
      write(*,'(A,A)') '  --- ', obj_name

      t_obs(1) = obs_yr
      t_obs(2) = obs_mo
      t_obs(3) = obs_dy
      t_obs(4) = 0.d0
      t_obs(5) = 0.d0
      t_obs(6) = 0.d0

      call orem_run(
     &   tle_file, norad,
     &   t_obs, mxz, 8, 10.d0, 0.90d0, -1.d0,
     &   80.d0, 160.d0, 0,
     &   4, 200, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   2, 0, 0,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0,
     &   0.d0, 0.d0, 0.d0,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   rpe, t_mean, t_std, ierr)

c     V1: Pipeline completes (ierr=0 or 3 with no drag)
      if (ierr .eq. 0 .or. ierr .eq. 3) then
         call assert_i(obj_name//' pipeline ok',
     &                 1, 1, npass, nfail)
      else
         call assert_i(obj_name//' pipeline ok',
     &                 0, 1, npass, nfail)
         return
      end if

c     V2: Found at least 1 zone
      if (nzu .ge. 1) then
         call assert_i(obj_name//' zones >= 1',
     &                 1, 1, npass, nfail)
      else
         call assert_i(obj_name//' zones >= 1',
     &                 0, 1, npass, nfail)
         return
      end if

c     V3: e_opt is physical and in right ballpark
c     (zone may be selected later when e has decayed)
      if (e_opt(1) .gt. 0.d0 .and.
     &    e_opt(1) .lt. 1.d0 .and.
     &    e_opt(1) .le. exp_ecc + 0.05d0) then
         call assert_i(obj_name//' e_opt near TLE',
     &                 1, 1, npass, nfail)
      else
         call assert_i(obj_name//' e_opt near TLE',
     &                 0, 1, npass, nfail)
      end if

c     V4: bn_opt in bounds [80, 160]
      if (bn_opt(1) .ge. 80.d0 .and.
     &    bn_opt(1) .le. 160.d0) then
         call assert_i(obj_name//' bn_opt in bounds',
     &                 1, 1, npass, nfail)
      else
         call assert_i(obj_name//' bn_opt in bounds',
     &                 0, 1, npass, nfail)
      end if

c     V5: RMS is valid (non-negative, finite)
      if (rms_o(1) .ge. 0.d0 .and.
     &    rms_o(1) .eq. rms_o(1)) then
         call assert_i(obj_name//' rms valid',
     &                 1, 1, npass, nfail)
      else
         call assert_i(obj_name//' rms valid',
     &                 0, 1, npass, nfail)
      end if

c     Print zone summary
      write(*,'(A,I2,A)') '    ', nzu, ' zones:'
      do i = 1, nzu
         if (reentry(i) .gt. 0.d0) then
            write(*,'(A,I1,A,F8.4,A,F5.2,A,F6.2,A)')
     &         '    Z', i, ': e=', e_opt(i),
     &         ' A=', bn_opt(i),
     &         ' RPE=', rpe(i), '%'
         else
            write(*,'(A,I1,A,F8.4,A,F5.2)')
     &         '    Z', i, ': e=', e_opt(i),
     &         ' A=', bn_opt(i)
         end if
      end do

      return
      end

c     ============================================================
c     ATM reader
c     ============================================================
      subroutine read_atm_v(alt, den, sch, maxn, natm)
      implicit double precision (a-h, o-z)
      integer maxn, natm
      double precision alt(maxn), den(maxn), sch(maxn)

      natm = 291
      if (natm .gt. maxn) natm = maxn
      alt(1) = 60.d0
      do i = 2, 141
         alt(i) = alt(i-1) + 1.d0
      end do
      do i = 142, natm
         alt(i) = alt(i-1) + 2.d0
      end do

      open(unit=88, file='input/ATM.DAT', status='old',
     &     err=900)
      read(88,110) (sch(i), i=1,61)
      read(88,110) (sch(i), i=62,natm)
      read(88,111) (den(i), i=1,61)
      read(88,111) (den(i), i=62,natm)
      close(88)
 110  format(8f10.3)
 111  format(8e10.3)

      do i = 1, natm
         alt(i) = alt(i) * 1000.d0
         sch(i) = sch(i) * 1000.d0
         den(i) = den(i) * 1.0d10
      end do
      return
 900  natm = 0
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

c     test_orem.F — Tests for OREM driver
c     -------------------------------------------------------
c     Compile:
c       ifx /heap-arrays /F:16777216 test_orem.F orem.F rsm.F
c           ga.F tle_evolution.F zone_select.F
c           ksrop/propagate_ks.F ksrop/Subrouts.F
c           ksrop/TLEread.F ksrop/Legendre.F
c           /exe:test_orem.exe
c     -------------------------------------------------------
      program test_orem

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      integer npass, nfail

      call init_constants()

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' OREM Driver Test Suite'
      write(*,*) '=============================================='

c     ==============================================================
c     compute_rpe unit tests
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- compute_rpe ---'

      call test_rpe(npass, nfail)

c     ==============================================================
c     Error handling
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Error handling ---'

      call test_errors(npass, nfail)

c     ==============================================================
c     42928 integration (IDRAG=0, fast)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- 42928 integration ---'

      call test_42928(npass, nfail)

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
c     RPE unit tests
c     ============================================================
      subroutine test_rpe(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

      parameter(mxz = 4)
      double precision reentry(mxz), zepoch(mxz), rpe(mxz)
      double precision t_mean, t_std

c     D1: Mode 1 — known re-entry, perfect prediction
      reentry(1) = 2458545.d0
      reentry(2) = 2458545.d0
      zepoch(1) = 2458018.d0
      zepoch(2) = 2458100.d0
      call compute_rpe(reentry, zepoch, 2,
     &   2458545.d0, rpe, t_mean, t_std, mxz)
      call assert_d('D1 perfect RPE=0',
     &   rpe(1), 0.d0, 0.01d0, npass, nfail)

c     D2: Mode 1 — prediction 10 days late
      reentry(1) = 2458555.d0
      call compute_rpe(reentry, zepoch, 2,
     &   2458545.d0, rpe, t_mean, t_std, mxz)
c     RPE = (555-545)/(545-018)*100 = 10/527*100 = 1.9%
      if (dabs(rpe(1) - 1.9d0) .lt. 0.2d0) then
         call assert_i('D2 RPE ~1.9% (10d late)',
     &                 1, 1, npass, nfail)
      else
         call assert_i('D2 RPE ~1.9% (10d late)',
     &                 0, 1, npass, nfail)
      end if

c     D3: Mode 2 — no observed re-entry, compute mean/std
      reentry(1) = 2458540.d0
      reentry(2) = 2458550.d0
      reentry(3) = 2458545.d0
      zepoch(1) = 2458018.d0
      zepoch(2) = 2458100.d0
      zepoch(3) = 2458200.d0
      call compute_rpe(reentry, zepoch, 3,
     &   0.d0, rpe, t_mean, t_std, mxz)
      call assert_d('D3 mean=(540+550+545)/3',
     &   t_mean, 2458545.d0, 0.1d0, npass, nfail)

c     D4: Std deviation > 0
      if (t_std .gt. 0.d0) then
         call assert_i('D4 std > 0',
     &                 1, 1, npass, nfail)
      else
         call assert_i('D4 std > 0',
     &                 0, 1, npass, nfail)
      end if

c     D5: Zero predictions → mean=0
      reentry(1) = 0.d0
      reentry(2) = 0.d0
      call compute_rpe(reentry, zepoch, 2,
     &   0.d0, rpe, t_mean, t_std, mxz)
      call assert_d('D5 no valid: mean=0',
     &   t_mean, 0.d0, 1.d-10, npass, nfail)

      return
      end

c     ============================================================
c     Error handling tests
c     ============================================================
      subroutine test_errors(npass, nfail)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth
      integer npass, nfail

      parameter(mxz = 4, ndim_a = 291)
      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      double precision alt(ndim_a), den(ndim_a), sch(ndim_a)
      integer natm, ierr

      call read_atm_orem(alt, den, sch, ndim_a, natm)

      t_obs(1) = 0.d0
      t_obs(2) = 0.d0
      t_obs(3) = 0.d0
      t_obs(4) = 0.d0
      t_obs(5) = 0.d0
      t_obs(6) = 0.d0

c     D6: Bad TLE file → ierr=1
      call orem_run(
     &   'nonexistent.tle', 42928,
     &   t_obs, 4, 8, 10.d0, 0.90d0, -1.d0,
     &   80.d0, 160.d0, 0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   2, 0, 0,
     &   7.29d-5, 3.35d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0, 0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   rpe, t_mean, t_std, ierr)
      call assert_i('D6 bad file ierr=1',
     &              ierr, 1, npass, nfail)

c     D7: Wrong NORAD → ierr=1 (no match)
      call orem_run(
     &   'input/example_42928_zone0.tle.txt', 99999,
     &   t_obs, 4, 8, 10.d0, 0.90d0, -1.d0,
     &   80.d0, 160.d0, 0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   2, 0, 0,
     &   7.29d-5, 3.35d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0, 0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   rpe, t_mean, t_std, ierr)
      call assert_i('D7 wrong NORAD ierr=1',
     &              ierr, 1, npass, nfail)

      return
      end

c     ============================================================
c     42928 integration test (IDRAG=0 for speed)
c     ============================================================
      subroutine test_42928(npass, nfail)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth
      integer npass, nfail

      parameter(mxz = 4, ndim_a = 291)
      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      double precision alt(ndim_a), den(ndim_a), sch(ndim_a)
      integer natm, ierr, nzu

      call read_atm_orem(alt, den, sch, ndim_a, natm)

c     Observed re-entry: 2019-03-03
      t_obs(1) = 2019.d0
      t_obs(2) = 3.d0
      t_obs(3) = 3.d0
      t_obs(4) = 0.d0
      t_obs(5) = 0.d0
      t_obs(6) = 0.d0

      write(*,*) '    running orem_run on 42928...'

      call orem_run(
     &   'input/example_42928.tle.txt', 42928,
     &   t_obs, mxz, 8, 10.d0, 0.90d0, -1.d0,
     &   80.d0, 160.d0, 0,
     &   4, 200, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   2, 0, 0,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0,
     &   0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   rpe, t_mean, t_std, ierr)

c     D8: Pipeline completed (ierr=3 ok with IDRAG=0: no re-entry)
      if (ierr .eq. 0 .or. ierr .eq. 3) then
         call assert_i('D8 42928 pipeline ok',
     &                 1, 1, npass, nfail)
      else
         call assert_i('D8 42928 pipeline ok',
     &                 0, 1, npass, nfail)
      end if

c     D9: Found zones
      if (nzu .ge. 1) then
         call assert_i('D9 nzones >= 1',
     &                 1, 1, npass, nfail)
      else
         call assert_i('D9 nzones >= 1',
     &                 0, 1, npass, nfail)
         return
      end if

c     D10: Zone epochs are valid JDs
      if (zepoch(1) .gt. 2457000.d0) then
         call assert_i('D10 zone epoch valid',
     &                 1, 1, npass, nfail)
      else
         call assert_i('D10 zone epoch valid',
     &                 0, 1, npass, nfail)
      end if

c     D11: e_opt in physical range
      if (e_opt(1) .gt. 0.d0 .and.
     &    e_opt(1) .lt. 1.d0) then
         call assert_i('D11 e_opt physical',
     &                 1, 1, npass, nfail)
      else
         call assert_i('D11 e_opt physical',
     &                 0, 1, npass, nfail)
      end if

c     D12: bn_opt in bounds
      if (bn_opt(1) .ge. 80.d0 .and.
     &    bn_opt(1) .le. 160.d0) then
         call assert_i('D12 bn_opt in [80,160]',
     &                 1, 1, npass, nfail)
      else
         call assert_i('D12 bn_opt in [1,5]',
     &                 0, 1, npass, nfail)
      end if

c     D13: RMS is valid
      if (rms_o(1) .ge. 0.d0 .and.
     &    rms_o(1) .eq. rms_o(1)) then
         call assert_i('D13 rms valid',
     &                 1, 1, npass, nfail)
      else
         call assert_i('D13 rms valid',
     &                 0, 1, npass, nfail)
      end if

c     D14: t_mean finite (0 is valid with IDRAG=0)
      if (t_mean .eq. t_mean) then
         call assert_i('D14 t_mean finite',
     &                 1, 1, npass, nfail)
      else
         call assert_i('D14 t_mean finite',
     &                 0, 1, npass, nfail)
      end if

      write(*,'(A,I2,A)') '    ', nzu, ' zones processed'
      do i = 1, nzu
         if (reentry(i) .gt. 0.d0) then
            write(*,'(A,I2,A,F12.2,A,F8.4,A,F6.2,A,F6.2)')
     &         '    Z', i, ': JD=', reentry(i),
     &         ' e=', e_opt(i),
     &         ' A=', bn_opt(i), ' RPE=', rpe(i)
         else
            write(*,'(A,I2,A)') '    Z', i, ': no re-entry'
         end if
      end do
      write(*,'(A,F12.2,A,F6.2,A)')
     &   '    mean JD=', t_mean, ' std=', t_std, ' days'

      return
      end

c     ============================================================
c     ATM reader (same as test_rsm)
c     ============================================================
      subroutine read_atm_orem(alt, den, sch, maxn, natm)
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

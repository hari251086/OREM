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
c     Issue #12: RSM/GA failure recovery and diagnostics
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Failure recovery / diagnostics (#12) ---'

      call test_failure_modes(npass, nfail)

c     ==============================================================
c     Issue #12 (G2): physics-based BN floor
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- G2 physics-based BN floor ---'

      call test_bn_floor(npass, nfail)

c     ==============================================================
c     Issue #13: structured output report
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Prediction report (#13) ---'

      call test_report(npass, nfail)

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
      integer zstat(mxz), nzv

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
     &   20, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   2, 0, 0,
     &   7.29d-5, 3.35d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0, 0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
     &   rpe, t_mean, t_std, ierr)
      call assert_i('D6 bad file ierr=1',
     &              ierr, 1, npass, nfail)

c     D7: Wrong NORAD → ierr=1 (no match)
      call orem_run(
     &   'input/example_42928_zone0.tle.txt', 99999,
     &   t_obs, 4, 8, 10.d0, 0.90d0, -1.d0,
     &   80.d0, 160.d0, 0,
     &   20, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   2, 0, 0,
     &   7.29d-5, 3.35d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0, 0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
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
      integer zstat(mxz), nzv

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
     &   20, 200, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   2, 0, 0,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0,
     &   0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
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

c     D12: bn_opt physical (positive, finite). Not bounds-checked
c     against the caller's [80,160]: the G2 floor legitimately extends
c     the zone-1 search below 80, and with the corrected J71 ATM.DAT
c     (v1.17) 42928's fit-consistent BN (~70-85) actually lands there.
      if (bn_opt(1) .gt. 0.d0 .and. bn_opt(1) .eq. bn_opt(1)) then
         call assert_i('D12 bn_opt physical',
     &                 1, 1, npass, nfail)
      else
         call assert_i('D12 bn_opt physical',
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
c     Issue #12 failure-mode tests.
c
c     BN enters the drag term as propagate_ks.F's V/(BN*2.d0), so
c     BN=0 forces a division-by-zero -> NaN altitude -> exit_code=2
c     on every RSM grid point. This is a deterministic, physically
c     grounded way to force the propagator-divergence path without
c     hand-crafting fragile orbital elements.
c     ============================================================
      subroutine test_failure_modes(npass, nfail)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth
      integer npass, nfail

      parameter(mxz = 4, ndim_a = 291)
      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      double precision alt(ndim_a), den(ndim_a), sch(ndim_a)
      integer natm, ierr, nzu
      integer zstat(mxz), nzv
      integer i, allfail

      call read_atm_orem(alt, den, sch, ndim_a, natm)

      t_obs(1) = 0.d0
      t_obs(2) = 0.d0
      t_obs(3) = 0.d0
      t_obs(4) = 0.d0
      t_obs(5) = 0.d0
      t_obs(6) = 0.d0

c     D15: BN=0 forces propagator divergence on every RSM grid point
c     for the single zone in this file -> RSM ierr=2 -> zone_status=1.
c     Also confirms the call completes without crashing.
      call orem_run(
     &   'input/example_42928_zone0.tle.txt', 42928,
     &   t_obs, mxz, 8, 10.d0, 0.90d0, -1.d0,
     &   0.d0, 0.d0, 0,
     &   20, 200, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   2, 0, 0,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0,
     &   0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
     &   rpe, t_mean, t_std, ierr)

      call assert_i('D15 propfail zone_status=1',
     &              zstat(1), 1, npass, nfail)
      call assert_i('D15 propfail excluded from nzones_valid',
     &              nzv, 0, npass, nfail)

c     D16: BN search window pinned to [20, 30], well BELOW the zone's
c     fit-consistent BN (~70-85 with the corrected J71 ATM.DAT, v1.17)
c     so the GA optimum pins at the UPPER edge and trips the 15%
c     boundary check. The window must sit below the real BN, not above
c     it: the G2 physics floor (v1.12) extends bn_lo downward on zone 1
c     (here to ~12.5), so a window pinned above the real BN simply
c     un-pins at the bottom and the optimum lands in the interior --
c     that configuration stopped flagging once the corrected table
c     moved the real BN below 200 (previously ~300 table-consistent).
c     Verifies detection + warning; the driver loop must still
c     complete normally.
      call orem_run(
     &   'input/example_42928.tle.txt', 42928,
     &   t_obs, mxz, 8, 10.d0, 0.90d0, -1.d0,
     &   20.d0, 30.d0, 0,
     &   20, 200, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   2, 0, 0,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0,
     &   0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
     &   rpe, t_mean, t_std, ierr)

      call assert_i('D16 GA boundary zone_status=2',
     &              zstat(1), 2, npass, nfail)

c     D17: same BN=0 divergence as D15, but on the full 4-zone 42928
c     file. Confirms zone failure doesn't crash the driver loop -- all
c     zones are still found and processed (just all flagged failed),
c     and nzones_valid correctly reports 0 when every zone fails.
      call orem_run(
     &   'input/example_42928.tle.txt', 42928,
     &   t_obs, mxz, 8, 10.d0, 0.90d0, -1.d0,
     &   0.d0, 0.d0, 0,
     &   20, 200, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   2, 0, 0,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0,
     &   0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
     &   rpe, t_mean, t_std, ierr)

      if (nzu .ge. 1) then
         call assert_i('D17 all-fail: zones still found',
     &                 1, 1, npass, nfail)
         allfail = 1
         do i = 1, nzu
            if (zstat(i) .ne. 1) allfail = 0
         end do
         call assert_i('D17 all-fail: every zone status=1',
     &                 allfail, 1, npass, nfail)
      else
         call assert_i('D17 all-fail: zones still found',
     &                 0, 1, npass, nfail)
      end if
      call assert_i('D17 all-fail: nzones_valid=0',
     &              nzv, 0, npass, nfail)

      return
      end

c     ============================================================
c     Issue #12 (G2) test: physics-based BN floor.
c
c     37151 Long March 3B's zone 1 (with the default bn_min_init=80,
c     same as test_reentry.F) is a verified case where the floor
c     estimate extends bn_lo well below 80 -- without G2, GA could
c     never search below the caller's bn_min_init, so bn_opt(1) < 80
c     here is only reachable because the floor extension fired.
c     ============================================================
      subroutine test_bn_floor(npass, nfail)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth
      integer npass, nfail

      parameter(mxz = 4, ndim_a = 291)
      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      double precision alt(ndim_a), den(ndim_a), sch(ndim_a)
      integer natm, ierr, nzu
      integer zstat(mxz), nzv

      call read_atm_orem(alt, den, sch, ndim_a, natm)

      t_obs(1) = 0.d0
      t_obs(2) = 0.d0
      t_obs(3) = 0.d0
      t_obs(4) = 0.d0
      t_obs(5) = 0.d0
      t_obs(6) = 0.d0

      call orem_run(
     &   'input/example_37151.tle.txt', 37151,
     &   t_obs, mxz, 8, 10.d0, 0.90d0, -1.d0,
     &   80.d0, 160.d0, 0,
     &   20, 200, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   2, 0, 0,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0,
     &   0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
     &   rpe, t_mean, t_std, ierr)

      if (bn_opt(1) .gt. 0.d0 .and. bn_opt(1) .eq. bn_opt(1)) then
         call assert_i('G2 37151 bn_opt physical',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G2 37151 bn_opt physical',
     &                 0, 1, npass, nfail)
      end if

      if (bn_opt(1) .lt. 80.d0) then
         call assert_i('G2 floor let bn_opt(1) go below 80',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G2 floor let bn_opt(1) go below 80',
     &                 0, 1, npass, nfail)
      end if

      return
      end

c     ============================================================
c     Issue #13 test: orem_report writes a well-formed report from
c     a real (fast, IDRAG=0) orem_run on the 42928 zone-0 file.
c     ============================================================
      subroutine test_report(npass, nfail)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth
      integer npass, nfail

      parameter(mxz = 4, ndim_a = 291)
      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      double precision alt(ndim_a), den(ndim_a), sch(ndim_a)
      integer natm, ierr, nzu, ierr_rep
      integer zstat(mxz), nzv

      character*200 line
      integer ihdr, izone, iens, ilgd

      call read_atm_orem(alt, den, sch, ndim_a, natm)

      t_obs(1) = 2019.d0
      t_obs(2) = 3.d0
      t_obs(3) = 3.d0
      t_obs(4) = 0.d0
      t_obs(5) = 0.d0
      t_obs(6) = 0.d0

      call orem_run(
     &   'input/example_42928_zone0.tle.txt', 42928,
     &   t_obs, mxz, 8, 10.d0, 0.90d0, -1.d0,
     &   80.d0, 160.d0, 0,
     &   20, 200, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   2, 0, 0,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0,
     &   0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
     &   rpe, t_mean, t_std, ierr)

      call orem_report(
     &   'output/test_report.txt', 42928, t_obs,
     &   80.d0, 160.d0, 0,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
     &   rpe, t_mean, t_std, mxz, ierr_rep)

c     R1: report written without error
      call assert_i('R1 report ierr_rep=0', ierr_rep, 0,
     &              npass, nfail)

c     R2-R4: report contains header, zone table, status legend
      ihdr  = 0
      izone = 0
      iens  = 0
      ilgd  = 0
      open(unit=45, file='output/test_report.txt', status='old',
     &     err=800)
 10   read(45,'(A)',end=20) line
      if (index(line, 'OREM Re-entry Prediction Report') .gt. 0)
     &   ihdr = 1
      if (index(line, 'Zone  Epoch (JD)') .gt. 0) izone = 1
      if (index(line, 'Ensemble') .gt. 0) iens = 1
      if (index(line, 'zone_status legend') .gt. 0) ilgd = 1
      goto 10
 20   close(45)

      call assert_i('R2 report header present', ihdr, 1,
     &              npass, nfail)
      call assert_i('R3 zone table present', izone, 1,
     &              npass, nfail)
      call assert_i('R4 ensemble + legend present', iens*ilgd, 1,
     &              npass, nfail)

c     R5-R7: with-re-entry path via synthetic zone data (fast; no
c     propagation). 3 zones with predictions, observed date between
c     the Z3 (latest) and mean predictions. Latest-zone must be the
c     PRIMARY estimate (v1.20: it wins the 7-object evaluation with
c     median |RPE| 8.2% vs the uniform mean's 45% mean error).
      nzu = 3
      zepoch(1)  = 2458000.d0
      zepoch(2)  = 2458100.d0
      zepoch(3)  = 2458200.d0
      reentry(1) = 2458700.d0
      reentry(2) = 2458600.d0
      reentry(3) = 2458560.d0
      e_opt(1) = 0.30d0
      e_opt(2) = 0.25d0
      e_opt(3) = 0.20d0
      bn_opt(1) = 70.d0
      bn_opt(2) = 60.d0
      bn_opt(3) = 55.d0
      rms_o(1) = 0.1d0
      rms_o(2) = 0.1d0
      rms_o(3) = 0.1d0
      zstat(1) = 0
      zstat(2) = 0
      zstat(3) = 0
      rpe(1) = 28.3d0
      rpe(2) = 10.0d0
      rpe(3) = 2.7d0
      t_mean = (reentry(1)+reentry(2)+reentry(3))/3.d0
      t_std  = 59.3d0

      call orem_report(
     &   'output/test_report2.txt', 42928, t_obs,
     &   80.d0, 160.d0, 1,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, 3,
     &   rpe, t_mean, t_std, mxz, ierr_rep)

      call assert_i('R5 synthetic report ierr_rep=0', ierr_rep, 0,
     &              npass, nfail)

      ihdr  = 0
      izone = 0
      open(unit=45, file='output/test_report2.txt', status='old',
     &     err=810)
 30   read(45,'(A)',end=40) line
      if (index(line, 'PRIMARY estimate (latest zone, Z 3') .gt. 0)
     &   ihdr = 1
      if (index(line, 'latest-zone RPE') .gt. 0) izone = 1
      goto 30
 40   close(45)
      call assert_i('R6 PRIMARY = latest zone Z3', ihdr, 1,
     &              npass, nfail)
      call assert_i('R7 latest-zone RPE line present', izone, 1,
     &              npass, nfail)
      return

 800  call assert_i('R2 report file readable', 0, 1, npass, nfail)
      return
 810  call assert_i('R5 synthetic report readable', 0, 1,
     &              npass, nfail)
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

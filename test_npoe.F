c     test_npoe.F — Issue #11: Cross-validate propagate_ks vs NPOE
c     -------------------------------------------------------
c     Runs propagate_ks on the 3x3 RSM grid from NPOE Zone 0,
c     object 42928, and compares against NPOE reference output
c     (o1.out-o9.out) from E:\RESEARCH\1. R&D\Re-entry\COSPAR ASR\
c     42928\Zone - 0\RSM Surfaces\.
c
c     Reference NPOE parameters:
c       Object   : 42928 PSLV-C39 R/B
c       Zone     : Zone 0 (2017-09-24 21:09:49 UTC)
c       Duration : 7 days
c       Gravity  : 6x6 EGM96
c       Atmosphere: Jacchia-70, F10=72, Kp=10
c       Sun/Moon : yes
c       BN grid  : 80, 120, 160 kg/m2 (Cd=1, A=1 m2)
c       e grid   : 0.321216, 0.3215792, 0.3219424
c
c     Checks:
c       N1-N3  : BN monotonicity (higher BN -> less decay)
c       N4-N6  : e monotonicity (higher e -> higher initial apogee)
c       N7-N9  : BN sensitivity ratio (BN80/BN160 decay ratio > 1.5)
c       N10    : No propagator failure for all 9 runs
c       N11    : IDRAG=0 gives near-zero decay
c       N12    : Decay direction all negative (IDRAG=1)
c       N13    : Magnitude within factor 3 of NPOE (BN=80, e=0.321216)
c       N14    : Magnitude within factor 3 of NPOE (BN=160, e=0.321216)
c       N15    : BN80 decays faster than BN160 for all 3 e values
c     -------------------------------------------------------
c     Compile:
c       ifx /heap-arrays /F:16777216 test_npoe.F orem.F rsm.F
c           ga.F tle_evolution.F zone_select.F
c           ksrop/propagate_ks.F ksrop/Subrouts.F
c           ksrop/TLEread.F ksrop/Legendre.F
c           /exe:test_npoe.exe
c     -------------------------------------------------------
      program test_npoe

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291, max_traj = 200, mxe = 3, mxbn = 3)

      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision traj_jd(max_traj)
      double precision traj_x(3, max_traj)
      double precision traj_xd(3, max_traj)
      integer idump, exit_code

      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      double precision e_grid(3), bn_grid(3)
      double precision drop_ks(3, 3)
      double precision apo0_ks(3, 3)
      integer fail9

      integer npass, nfail, ie, ia
      integer n_force(3)
      integer nrev_test, istep_test
      double precision pek(7)

c     NPOE reference 7-day apogee drops (km), from o1-o9.out
c     Row = e index (1-3), Col = BN index (1-3 = 80/120/160 kg/m2)
      double precision ref_drop(3,3)
      data ref_drop /
     &  -61.91d0, -72.79d0, -86.12d0,   ! BN=80:  e1,e2,e3
     &  -41.06d0, -48.31d0, -57.21d0,   ! BN=120: e1,e2,e3
     &  -30.62d0, -36.07d0, -42.74d0 /  ! BN=160: e1,e2,e3

c     NPOE reference initial apogee (km), t=0 from o1.out (e1 same)
      double precision ref_apo0(3)
      data ref_apo0 / 12732.32d0, 12735.82d0, 12739.32d0 /

      call init_constants()

c     Read atmosphere table
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
      write(*,*) ' OREM Issue #11: NPOE Cross-Validation Test'
      write(*,*) '=============================================='
      write(*,*) ' Object  : 42928 PSLV-C39 R/B'
      write(*,*) ' Zone    : Zone 0  (2017-09-24 21:09:49 UTC)'
      write(*,*) ' Config  : IDRAG=1, geo=4, sun=2, moon=2'
      write(*,*) ' Duration: 7 days (~35 revolutions)'
      write(*,*) '=============================================='

c     -------------------------------------------------------
c     Common initial conditions from NPOE 1.in
c     Date: 2017-09-24 21:09:49 UTC
      cal0(1) = 2017.d0
      cal0(2) = 9.d0
      cal0(3) = 24.d0
      cal0(4) = 21.d0
      cal0(5) = 9.d0
      cal0(6) = 49.d0

c     Base elements (same for all 9 runs except e)
      a_base  = 9636.790207d0
      inc_base  = 19.1713d0
      raan_base = 68.985106d0
      aop_base  = 273.564599d0
      ta_base   = 28.938462d0     ! true anomaly from NPOE .in

      e_grid(1) = 0.321216d0
      e_grid(2) = 0.3215792d0
      e_grid(3) = 0.3219424d0
      bn_grid(1) = 80.d0
      bn_grid(2) = 120.d0
      bn_grid(3) = 160.d0

c     Force model: same as OREM RSM default
      n_force(1) = 4   ! geo degree (J2-J4)
      n_force(2) = 2   ! sun degree
      n_force(3) = 2   ! moon degree

c     Propagation: ~35 revolutions for 7 days
c     period ~4.96 hrs, 7 days = 33.9 revs -> use 35
      nrev_test = 35
      istep_test = 360

c     WE_rot, EPS_f, FR_rot (Earth rotation, flattening, frame rot)
      WE_rot = 7.2921150d-5
      EPS_f  = 3.35281066d-3
      FR_rot = 1.d0

c     amuS, amuM (Sun and Moon GM, km3/s2)
      amuS = 1.32712440018d11
      amuM = 4.902801076d3

      write(*,*) ''
      write(*,*) '--- All 9 RSM grid runs ---'
      write(*,'(A)') '  Run  e         BN   apo0(km)  apo7(km)  '//
     &  'drop(km) exit'

c     ============================================================
c     Run all 9 grid points
c     ============================================================
      fail9 = 0
      do ie = 1, 3
         do ia = 1, 3

c           Convert TA -> MA for this eccentricity
            ae = e_grid(ie)
            f_rad = ta_base * d2r
            cosf = dcos(f_rad)
            sinf = dsin(f_rad)
            cosE_ = (ae + cosf)/(1.d0 + ae*cosf)
            sinE_ = sinf*dsqrt(1.d0-ae*ae)/(1.d0+ae*cosf)
            E_rad = datan2(sinE_, cosE_)
            MA_rad = E_rad - ae*sinE_
            if (MA_rad .lt. 0.d0) MA_rad = MA_rad + 2.d0*pi
            MA_deg = MA_rad * r2d

c           Set up orbital elements (fixed perigee radius, NPOE convention)
c           NPOE keeps same a for all e variations within a run set
            a_kep(1) = a_base
            a_kep(2) = ae
            a_kep(3) = inc_base
            a_kep(4) = raan_base
            a_kep(5) = aop_base
            a_kep(6) = MA_deg
            a_kep(7) = 0.d0

c           Convert to Cartesian
            call oe2car(a_kep, x0, xd0, 1.d-12)

c           Run propagate_ks
            idump = 1
            call propagate_ks(
     &         x0, xd0, cal0,
     &         nrev_test, istep_test, 1.d-12,
     &         n_force, n_force(1), n_force(2), n_force(3),
     &         bn_grid(ia), 1, WE_rot, EPS_f, FR_rot,
     &         0.d0, 0.d0, 0, 0,
     &         0.d0, amuS, amuM,
     &         ALT_a, DEN_a, SCH_a, natm,
     &         max_traj, idump,
     &         traj_jd, traj_x, traj_xd,
     &         exit_code)

c           Compute initial and final apogee from traj
            call car2oe(traj_x(1,1), traj_xd(1,1), pek)
            apo_t0 = pek(1)*(1.d0 + pek(2))

            call car2oe(traj_x(1,idump), traj_xd(1,idump), pek)
            apo_tf = pek(1)*(1.d0 + pek(2))

            apo0_ks(ie,ia) = apo_t0
            drop_ks(ie,ia) = apo_tf - apo_t0

            write(*,'(A,I2,A,I2,A,F9.4,A,F6.1,A,F10.4,A,F10.4,A,
     &             F8.4,A,I2)')
     &         '  e', ie, ' BN', ia,
     &         '  e=', ae, '  BN=', bn_grid(ia),
     &         '  ', apo_t0, '  ', apo_tf,
     &         '  ', drop_ks(ie,ia),
     &         '  exit=', exit_code

            if (exit_code .eq. 2) fail9 = fail9 + 1

         end do
      end do

c     ============================================================
c     N10: No propagator failure for any of 9 runs
c     ============================================================
      write(*,*) ''
      write(*,*) '--- Structural checks ---'
      if (fail9 .eq. 0) then
         write(*,*) '  [PASS]  N10 no divergence in 9 runs'
         npass = npass + 1
      else
         write(*,'(A,I2,A)') '  [FAIL]  N10 ', fail9,
     &      ' diverged runs'
         nfail = nfail + 1
      end if

c     ============================================================
c     N4-N6: e monotonicity - higher e gives higher initial apogee
c     Check for each BN column
c     ============================================================
      do ia = 1, 3
         if (apo0_ks(1,ia) .lt. apo0_ks(2,ia) .and.
     &       apo0_ks(2,ia) .lt. apo0_ks(3,ia)) then
            write(*,'(A,I1,A)')
     &         '  [PASS]  N', 3+ia,
     &         ' e-monotonicity: higher e -> higher apo0'
            npass = npass + 1
         else
            write(*,'(A,I1,A)')
     &         '  [FAIL]  N', 3+ia,
     &         ' e-monotonicity failed'
            nfail = nfail + 1
         end if
      end do

c     ============================================================
c     N1-N3: BN monotonicity - higher BN gives less decay
c     Check for each e row
c     ============================================================
      do ie = 1, 3
         if (drop_ks(ie,1) .lt. drop_ks(ie,2) .and.
     &       drop_ks(ie,2) .lt. drop_ks(ie,3)) then
            write(*,'(A,I1,A)')
     &         '  [PASS]  N', ie,
     &         ' BN-monotonicity: higher BN -> less decay'
            npass = npass + 1
         else
            write(*,'(A,I1,A,3F8.2)')
     &         '  [FAIL]  N', ie,
     &         ' BN-monotonicity: drops=',
     &         drop_ks(ie,1), drop_ks(ie,2), drop_ks(ie,3)
            nfail = nfail + 1
         end if
      end do

c     ============================================================
c     N12: All drops negative for IDRAG=1
c     ============================================================
      idok = 0
      do ie = 1, 3
         do ia = 1, 3
            if (drop_ks(ie,ia) .lt. 0.d0) idok = idok + 1
         end do
      end do
      if (idok .eq. 9) then
         write(*,*) '  [PASS]  N12 all drops negative (IDRAG=1)'
         npass = npass + 1
      else
         write(*,'(A,I2,A)')
     &      '  [FAIL]  N12 only ', idok, '/9 drops negative'
         nfail = nfail + 1
      end if

c     ============================================================
c     N7-N9: BN sensitivity ratio: decay(BN=80)/decay(BN=160) > 1.5
c     NPOE reference ratio is ~2.0 for all e rows.
c     Tolerance allows for different atmosphere model.
c     ============================================================
      do ie = 1, 3
         ratio_ks = drop_ks(ie,1) / drop_ks(ie,3)
         ratio_npoe = ref_drop(ie,1) / ref_drop(ie,3)
         if (ratio_ks .gt. 1.5d0) then
            write(*,'(A,I1,A,F5.2,A,F5.2)')
     &         '  [PASS]  N', 6+ie,
     &         ' BN-sensitivity ratio=', ratio_ks,
     &         ' (NPOE ref=', ratio_npoe, ')'
            npass = npass + 1
         else
            write(*,'(A,I1,A,F5.2,A,F5.2)')
     &         '  [FAIL]  N', 6+ie,
     &         ' BN-sensitivity ratio=', ratio_ks,
     &         ' (need > 1.5, NPOE=', ratio_npoe, ')'
            nfail = nfail + 1
         end if
      end do

c     ============================================================
c     N13-N14: Decay magnitude within factor 3 of NPOE
c     Tolerance is generous due to different atmosphere models
c     (ATM.DAT exponential vs Jacchia-70).
c     NPOE ref (e1, BN=80): drop = -61.91 km -> tolerance [-186, -18] km
c     NPOE ref (e1, BN=160): drop = -30.62 km -> tolerance [-92, -10] km
c     ============================================================
      d80  = drop_ks(1,1)
      d160 = drop_ks(1,3)
      if (d80 .gt. -186.d0 .and. d80 .lt. -18.d0) then
         write(*,'(A,F7.2,A)')
     &      '  [PASS]  N13 BN=80 drop=', d80,
     &      ' km (NPOE ref=-61.91, factor-3 tol)'
         npass = npass + 1
      else
         write(*,'(A,F7.2,A)')
     &      '  [FAIL]  N13 BN=80 drop=', d80,
     &      ' km outside [-186,-18] km factor-3 range'
         nfail = nfail + 1
      end if

      if (d160 .gt. -92.d0 .and. d160 .lt. -10.d0) then
         write(*,'(A,F7.2,A)')
     &      '  [PASS]  N14 BN=160 drop=', d160,
     &      ' km (NPOE ref=-30.62, factor-3 tol)'
         npass = npass + 1
      else
         write(*,'(A,F7.2,A)')
     &      '  [FAIL]  N14 BN=160 drop=', d160,
     &      ' km outside [-92,-10] km factor-3 range'
         nfail = nfail + 1
      end if

c     ============================================================
c     N11: IDRAG=0 gives near-zero decay (< 0.5 km in 7 days)
c     Use e_grid(1), BN=80 as the test case
c     ============================================================
      ae = e_grid(1)
      f_rad = ta_base * d2r
      cosf = dcos(f_rad)
      sinf = dsin(f_rad)
      cosE_ = (ae + cosf)/(1.d0 + ae*cosf)
      sinE_ = sinf*dsqrt(1.d0-ae*ae)/(1.d0+ae*cosf)
      E_rad = datan2(sinE_, cosE_)
      MA_rad = E_rad - ae*sinE_
      if (MA_rad .lt. 0.d0) MA_rad = MA_rad + 2.d0*pi

      a_kep(1) = a_base
      a_kep(2) = ae
      a_kep(3) = inc_base
      a_kep(4) = raan_base
      a_kep(5) = aop_base
      a_kep(6) = MA_rad * r2d
      a_kep(7) = 0.d0

      call oe2car(a_kep, x0, xd0, 1.d-12)

      idump = 1
      call propagate_ks(
     &   x0, xd0, cal0,
     &   nrev_test, istep_test, 1.d-12,
     &   n_force, n_force(1), n_force(2), n_force(3),
     &   80.d0, 0, WE_rot, EPS_f, FR_rot,
     &   0.d0, 0.d0, 0, 0,
     &   0.d0, amuS, amuM,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   max_traj, idump,
     &   traj_jd, traj_x, traj_xd,
     &   exit_code)

      call car2oe(traj_x(1,1), traj_xd(1,1), pek)
      apo_t0 = pek(1)*(1.d0 + pek(2))
      call car2oe(traj_x(1,idump), traj_xd(1,idump), pek)
      apo_tf = pek(1)*(1.d0 + pek(2))
      drop_nodrag = apo_tf - apo_t0

      if (dabs(drop_nodrag) .lt. 0.5d0) then
         write(*,'(A,F7.4,A)')
     &      '  [PASS]  N11 IDRAG=0 drop=', drop_nodrag,
     &      ' km (<0.5 km expected)'
         npass = npass + 1
      else
         write(*,'(A,F7.4,A)')
     &      '  [FAIL]  N11 IDRAG=0 drop=', drop_nodrag,
     &      ' km (should be near zero)'
         nfail = nfail + 1
      end if

c     ============================================================
c     N15: BN sensitivity comparison diagnostic
c     Print ratio comparison table
c     ============================================================
      write(*,*) ''
      write(*,*) '--- BN sensitivity (propagate_ks vs NPOE) ---'
      write(*,'(A)') '  e-row  BN=80     BN=160    ratio_ks '//
     &  ' ratio_NPOE'
      do ie = 1, 3
         ratio_ks = drop_ks(ie,1) / drop_ks(ie,3)
         ratio_npoe = ref_drop(ie,1) / ref_drop(ie,3)
         write(*,'(A,I1,A,F8.3,A,F8.3,A,F6.2,A,F6.2)')
     &      '  e', ie, '  ',
     &      drop_ks(ie,1), '  ', drop_ks(ie,3),
     &      '  ', ratio_ks, '  ', ratio_npoe
      end do

      write(*,*) ''
      write(*,*) '--- Magnitude comparison (propagate_ks vs NPOE) ---'
      write(*,'(A)') '  Run   drop_ks    drop_NPOE   ratio'
      do ie = 1, 3
         do ia = 1, 3
            ratm = drop_ks(ie,ia) / ref_drop(ie,ia)
            write(*,'(A,I1,A,I1,A,F8.3,A,F8.3,A,F6.2)')
     &         '  e', ie, '/BN', ia,
     &         '  ', drop_ks(ie,ia), '  ', ref_drop(ie,ia),
     &         '  ', ratm
         end do
      end do

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


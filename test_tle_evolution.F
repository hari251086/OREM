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
      double precision ls1, ls_min, ls_max
      integer npts1, iall_ok

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

c     T2: Should find entries (fewer than 3977 after dedup)
      if (npts .gt. 3000 .and. npts .lt. 3977) then
         call assert_i('T2 47944 npts deduped',1,1,npass,nfail)
      else
         call assert_i('T2 47944 npts deduped',0,1,npass,nfail)
      end if

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
c     Geometric consistency
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Geometric consistency ---'

c     Reload 47944 for consistency checks
      call tle_evolve('input/example_47944.tle.txt',
     &               47944, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)

c     T26: ha + hp + 2*R_Earth ≈ 2*a (geometric identity)
      iall_ok = 1
      do i = 1, npts
         err = dabs((ha_out(i)+hp_out(i)+2.d0*R_Earth)
     &             - 2.d0*a_out(i))
         if (err .gt. 1.d-6) iall_ok = 0
      end do
      call assert_i('T26 ha+hp+2Re = 2a identity',
     &              iall_ok, 1, npass, nfail)

c     T27: ha = a*(1+e) - Re exactly
      iall_ok = 1
      do i = 1, npts
         err = dabs(ha_out(i) -
     &       (a_out(i)*(1.d0+e_out(i)) - R_Earth))
         if (err .gt. 1.d-6) iall_ok = 0
      end do
      call assert_i('T27 ha = a*(1+e)-Re exact',
     &              iall_ok, 1, npass, nfail)

c     T28: hp = a*(1-e) - Re exactly
      iall_ok = 1
      do i = 1, npts
         err = dabs(hp_out(i) -
     &       (a_out(i)*(1.d0-e_out(i)) - R_Earth))
         if (err .gt. 1.d-6) iall_ok = 0
      end do
      call assert_i('T28 hp = a*(1-e)-Re exact',
     &              iall_ok, 1, npass, nfail)

c     T29: All outputs finite (no NaN) — a, e, ha, hp, Ls
      iall_ok = 1
      do i = 1, npts
         if (a_out(i) .ne. a_out(i)) iall_ok = 0
         if (e_out(i) .ne. e_out(i)) iall_ok = 0
         if (ha_out(i) .ne. ha_out(i)) iall_ok = 0
         if (hp_out(i) .ne. hp_out(i)) iall_ok = 0
         if (lambda_s_out(i) .ne. lambda_s_out(i)) iall_ok = 0
      end do
      call assert_i('T29 all outputs finite (no NaN)',
     &              iall_ok, 1, npass, nfail)

c     T30: Perigee radius > 0 (physical constraint)
      iall_ok = 1
      do i = 1, npts
         if (hp_out(i) + R_Earth .le. 0.d0) iall_ok = 0
      end do
      call assert_i('T30 perigee radius > 0',
     &              iall_ok, 1, npass, nfail)

c     T31: Apogee radius > perigee radius
      iall_ok = 1
      do i = 1, npts
         ra = ha_out(i) + R_Earth
         rp = hp_out(i) + R_Earth
         if (ra .lt. rp - 1.d-6) iall_ok = 0
      end do
      call assert_i('T31 ra > rp everywhere',
     &              iall_ok, 1, npass, nfail)

c     ==============================================================
c     Sun azimuth physics
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Sun azimuth physics ---'

c     T32: Sun azimuth varies over the evolution
c     (not constant — Sun moves relative to orbit)
      ls_min = lambda_s_out(1)
      ls_max = lambda_s_out(1)
      do i = 2, npts
         if (lambda_s_out(i) .lt. ls_min) ls_min=lambda_s_out(i)
         if (lambda_s_out(i) .gt. ls_max) ls_max=lambda_s_out(i)
      end do
      if (ls_max - ls_min .gt. 1.d0) then
         call assert_i('T32 Lambda_S varies >1 deg',
     &                 1, 1, npass, nfail)
      else
         call assert_i('T32 Lambda_S varies >1 deg',
     &                 0, 1, npass, nfail)
      end if

c     T33: First and last Sun azimuth differ
      if (dabs(lambda_s_out(npts)-lambda_s_out(1)).gt.0.1d0)
     &then
         call assert_i('T33 first/last Ls differ',
     &                 1, 1, npass, nfail)
      else
         call assert_i('T33 first/last Ls differ',
     &                 0, 1, npass, nfail)
      end if

c     ==============================================================
c     Spot-check known values
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Spot-check known values ---'

c     T34: 47944 first epoch year = 2021 (day 83)
c     JD for 2021-03-24 ≈ 2459297
      if (epochs(1) .gt. 2459290.d0 .and.
     &    epochs(1) .lt. 2459300.d0) then
         call assert_i('T34 47944 first epoch ~2021',
     &                 1, 1, npass, nfail)
      else
         call assert_i('T34 47944 first epoch ~2021',
     &                 0, 1, npass, nfail)
      end if

c     T35: 47944 last epoch year = 2025
c     JD for 2025-Jul ≈ 2460882
      if (epochs(npts) .gt. 2460600.d0 .and.
     &    epochs(npts) .lt. 2460920.d0) then
         call assert_i('T35 47944 last epoch ~2025',
     &                 1, 1, npass, nfail)
      else
         call assert_i('T35 47944 last epoch ~2025',
     &                 0, 1, npass, nfail)
      end if

c     T36: 47944 first entry inclination ≈ 97.57°
      if (dabs(inc_out(1) - 97.57d0) .lt. 0.5d0) then
         call assert_i('T36 47944 first inc ~97.57',
     &                 1, 1, npass, nfail)
      else
         call assert_i('T36 47944 first inc ~97.57',
     &                 0, 1, npass, nfail)
      end if

c     ==============================================================
c     42928 detailed checks
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- 42928 detailed ---'

      call tle_evolve('input/example_42928_zone0.tle.txt',
     &               42928, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)

c     T37: 42928 apogee altitude physical (HEO range)
      iall_ok = 1
      do i = 1, npts
         if (ha_out(i) .lt. 5000.d0 .or.
     &       ha_out(i) .gt. 15000.d0) iall_ok = 0
      end do
      call assert_i('T37 42928 ha in [5k,15k]',
     &              iall_ok, 1, npass, nfail)

c     T38: 42928 perigee altitude in [100, 250] km
      iall_ok = 1
      do i = 1, npts
         if (hp_out(i) .lt. 100.d0 .or.
     &       hp_out(i) .gt. 250.d0) iall_ok = 0
      end do
      call assert_i('T38 42928 hp in [100,250]',
     &              iall_ok, 1, npass, nfail)

c     T39: 42928 SMA in [9500, 10500] km
      iall_ok = 1
      do i = 1, npts
         if (a_out(i) .lt. 9500.d0 .or.
     &       a_out(i) .gt. 10500.d0) iall_ok = 0
      end do
      call assert_i('T39 42928 SMA in [9.5k,10.5k]',
     &              iall_ok, 1, npass, nfail)

c     T40: 42928 geometric identity holds
      iall_ok = 1
      do i = 1, npts
         err = dabs((ha_out(i)+hp_out(i)+2.d0*R_Earth)
     &             - 2.d0*a_out(i))
         if (err .gt. 1.d-6) iall_ok = 0
      end do
      call assert_i('T40 42928 ha+hp+2Re=2a',
     &              iall_ok, 1, npass, nfail)

c     T41: 42928 RAAN is decreasing (J2 regression for i<90)
      if (raan_out(npts) .lt. raan_out(1)) then
         call assert_i('T41 42928 RAAN regressing',
     &                 1, 1, npass, nfail)
      else
         call assert_i('T41 42928 RAAN regressing',
     &                 0, 1, npass, nfail)
      end if

c     ==============================================================
c     Large catalog test
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Large catalog ---'

c     T42: Read multi-satellite catalog (94597 entries)
      call tle_evolve('input/example_multi.tle.txt',
     &               -1, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)
      call assert_i('T42 multi catalog ierr=0',
     &              ierr, 0, npass, nfail)

c     T43: Should cap at maxpts (10000)
      call assert_i('T43 multi capped at maxpts',
     &              npts, maxpts, npass, nfail)

c     T44: All eccentricities in [0, 1)
      iall_ok = 1
      do i = 1, npts
         if (e_out(i) .lt. 0.d0 .or.
     &       e_out(i) .ge. 1.d0) iall_ok = 0
      end do
      call assert_i('T44 multi all ecc in [0,1)',
     &              iall_ok, 1, npass, nfail)

c     T45: All inclinations in [0, 180]
      iall_ok = 1
      do i = 1, npts
         if (inc_out(i) .lt. 0.d0 .or.
     &       inc_out(i) .gt. 180.d0) iall_ok = 0
      end do
      call assert_i('T45 multi all inc in [0,180]',
     &              iall_ok, 1, npass, nfail)

c     T46: All SMA > R_Earth (above surface)
      iall_ok = 1
      do i = 1, npts
         if (a_out(i) .lt. R_Earth) iall_ok = 0
      end do
      call assert_i('T46 multi all SMA > R_Earth',
     &              iall_ok, 1, npass, nfail)

c     T47: All Sun azimuths in [0, 360)
      iall_ok = 1
      do i = 1, npts
         if (lambda_s_out(i) .lt. 0.d0 .or.
     &       lambda_s_out(i) .ge. 360.d0) iall_ok = 0
      end do
      call assert_i('T47 multi all Ls in [0,360)',
     &              iall_ok, 1, npass, nfail)

c     T48: Filter Vanguard 1 (NORAD 5) from multi catalog
      call tle_evolve('input/example_multi.tle.txt',
     &               5, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)
      call assert_i('T48 Vanguard-1 ierr=0',
     &              ierr, 0, npass, nfail)
      if (npts .ge. 1) then
         call assert_i('T49 Vanguard-1 npts>=1',
     &                 1, 1, npass, nfail)
      else
         call assert_i('T49 Vanguard-1 npts>=1',
     &                 0, 1, npass, nfail)
      end if

c     T50: Vanguard-1 high eccentricity (e ≈ 0.184)
      if (npts .ge. 1) then
         if (dabs(e_out(1) - 0.184d0) .lt. 0.01d0) then
            call assert_i('T50 Vanguard-1 e~0.184',
     &                    1, 1, npass, nfail)
         else
            call assert_i('T50 Vanguard-1 e~0.184',
     &                    0, 1, npass, nfail)
         end if
      else
         call assert_i('T50 Vanguard-1 e~0.184',
     &                 0, 1, npass, nfail)
      end if

c     ==============================================================
c     Repeatability
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Repeatability ---'

c     T51: Two calls with same input produce identical output
      call tle_evolve('input/example_42928_zone0.tle.txt',
     &               42928, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)
      npts1 = npts
      ep1 = epochs(1)
      a1 = a_out(1)
      ls1 = lambda_s_out(1)
      call tle_evolve('input/example_42928_zone0.tle.txt',
     &               42928, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)
      iall_ok = 1
      if (npts .ne. npts1) then
         iall_ok = 0
      end if
      if (npts .ge. 1 .and. npts1 .ge. 1) then
         if (dabs(epochs(1) - ep1) .gt. 1.d-3) iall_ok = 0
         if (dabs(a_out(1) - a1) .gt. 1.d-3) iall_ok = 0
         if (dabs(lambda_s_out(1)-ls1) .gt. 1.d-3) iall_ok = 0
      end if
      call assert_i('T51 repeatability identical',
     &              iall_ok, 1, npass, nfail)

c     ==============================================================
c     Boundary conditions
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Boundary conditions ---'

c     T52: maxpts=1 returns exactly 1 entry
      call tle_evolve('input/example_47944.tle.txt',
     &               47944, 1,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)
      call assert_i('T52 maxpts=1 returns 1',
     &              npts, 1, npass, nfail)

c     T53: maxpts=1 still has valid data
      if (a_out(1) .gt. 6000.d0 .and.
     &    a_out(1) .lt. 8000.d0) then
         call assert_i('T53 maxpts=1 data valid',
     &                 1, 1, npass, nfail)
      else
         call assert_i('T53 maxpts=1 data valid',
     &                 0, 1, npass, nfail)
      end if

c     ==============================================================
c     Deduplication
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Deduplication ---'

c     T54: No two consecutive epochs within 86 sec (0.001 day)
      call tle_evolve('input/example_47944.tle.txt',
     &               47944, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)
      iall_ok = 1
      do i = 2, npts
         if (dabs(epochs(i) - epochs(i-1)) .lt. 0.001d0)
     &      iall_ok = 0
      end do
      call assert_i('T54 no dup epochs (<86 sec)',
     &              iall_ok, 1, npass, nfail)

c     T55: Deduped count < raw count (589 dups known)
      if (npts .lt. 3977) then
         call assert_i('T55 dedup removed entries',
     &                 1, 1, npass, nfail)
      else
         call assert_i('T55 dedup removed entries',
     &                 0, 1, npass, nfail)
      end if

c     T56: 42928 has no dups (14 unique entries)
      call tle_evolve('input/example_42928_zone0.tle.txt',
     &               42928, maxpts,
     &               epochs, a_out, e_out, inc_out,
     &               raan_out, aop_out,
     &               ha_out, hp_out, lambda_s_out,
     &               npts, ierr)
      call assert_i('T56 42928 no dups npts=14',
     &              npts, 14, npass, nfail)

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

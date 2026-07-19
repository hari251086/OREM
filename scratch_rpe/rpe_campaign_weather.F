c     rpe_campaign_weather.F -- issue #26 regression check: identical
c     to rpe_campaign.F (7-object drag-enabled RPE campaign) except it
c     calls sw_load/atm2d_load before the object loop, so every zone's
c     propagate_ks call uses epoch-resolved space weather instead of
c     the static ATM.DAT table -- the same auto-detect main_orem.F
c     does. Purpose: measure whether enabling weather regresses the
c     7-object validation set (issue #26's "no regression" DoD),
c     compared against rpe_campaign.exe rebuilt fresh from the same
c     source tree (both must be regenerated together -- the old
c     scratch_rpe/rpe_campaign_8zone_gated.csv on disk predates both
c     the v1.23 weather commit and the 2026-07-18 array-bound fixes,
c     so it is not a valid comparison baseline).
c
c     Compile:
c       ifx /heap-arrays /F:16777216 scratch_rpe\rpe_campaign_weather.F
c           swx.F orem.F rsm.F ga.F tle_evolution.F zone_select.F
c           ksrop\propagate_ks.F ksrop\Subrouts.F
c           ksrop\TLEread.F ksrop\Legendre.F
c           /exe:scratch_rpe\rpe_campaign_weather.exe
c     Run from the OREM repo root:
c       scratch_rpe\rpe_campaign_weather.exe
c     Output: scratch_rpe\rpe_campaign_weather.csv
c     -------------------------------------------------------
      program rpe_campaign_weather

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

c     mxz=8 (v1.21 experiment; was 4): zone_select returns the top-R2
c     candidates, so a higher cap admits LATER zones -- and v1.20
c     showed the latest zone is the primary estimator. This run
c     measures how much later the latest zone gets and what it does
c     to the latest-zone RPE.
      parameter(ndim_a = 291, mxz = 8, nobj = 7)

      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr, nzv
      integer zstat(mxz)

      character*40 tle_file(nobj)
      character*24 obj_name(nobj)
      integer norad(nobj)
      double precision oyr(nobj), omo(nobj), ody(nobj)
      double precision t_obs_jd, horiz, rpe_ens
      integer io, iz, nre, ierr_sw, ierr_a2

      data norad / 42928, 35497, 37151, 39615, 27526, 32007, 37819 /
      data oyr / 2019.d0, 2016.d0, 2015.d0, 2017.d0, 2012.d0,
     &           2010.d0, 2013.d0 /
      data omo / 3.d0, 10.d0, 12.d0, 9.d0, 5.d0, 6.d0, 9.d0 /
      data ody / 3.d0, 31.d0, 3.d0, 15.d0, 9.d0, 6.d0, 12.d0 /

      tle_file(1) = 'input/example_42928.tle.txt'
      tle_file(2) = 'input/example_35497.tle.txt'
      tle_file(3) = 'input/example_37151.tle.txt'
      tle_file(4) = 'input/example_39615.tle.txt'
      tle_file(5) = 'input/example_27526.tle.txt'
      tle_file(6) = 'input/example_32007.tle.txt'
      tle_file(7) = 'input/example_37819.tle.txt'
      obj_name(1) = '42928 PSLV-C39'
      obj_name(2) = '35497 Ariane 5 ESC-A'
      obj_name(3) = '37151 Long March 3B'
      obj_name(4) = '39615 Proton-M Briz-M'
      obj_name(5) = '27526 Ariane 5 R/B'
      obj_name(6) = '32007 GSLV R/B'
      obj_name(7) = '37819 Proton-M R/B'

      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

c     Issue #26: enable epoch-resolved weather (same auto-detect
c     pattern as main_orem.F). ALT_a/DEN_a/SCH_a above stay as the
c     static fallback table -- unused per-rev once sw_tinf/atm2d_interp
c     report loaded data, per propagate_ks.F's existing auto-detect.
      call sw_load('input/SW-All.csv', ierr_sw)
      call atm2d_load('input/ATM2D.DAT', ierr_a2)
      if (ierr_sw .eq. 0 .and. ierr_a2 .eq. 0) then
         write(*,'(A)') '[SW] epoch-resolved density: ENABLED'
      else
         write(*,'(A,I2,A,I2)') '[SW] epoch-resolved density:'//
     &      ' DISABLED -- sw_load ierr=', ierr_sw,
     &      ' atm2d_load ierr=', ierr_a2
         stop 1
      end if

      open(unit=33, file='scratch_rpe/rpe_campaign_weather.csv',
     &     status='unknown')
      write(33,'(A)') 'norad,zone,e_opt,bn_opt,reentry_jd,rpe_pct,'//
     &   'zstat,t_mean,t_std,ens_rpe_pct'

      do io = 1, nobj
         write(*,'(/,A,A)') ' ===== ', obj_name(io)

         t_obs(1) = oyr(io)
         t_obs(2) = omo(io)
         t_obs(3) = ody(io)
         t_obs(4) = 0.d0
         t_obs(5) = 0.d0
         t_obs(6) = 0.d0

         call orem_run(
     &      tle_file(io), norad(io),
     &      t_obs, mxz, 8, 10.d0, 0.90d0, -1.d0,
     &      80.d0, 160.d0, 1,
     &      20, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &      20, 2, 3,
     &      7.2921150d-5, 3.35281066d-3, 1.d0,
     &      1.2d0, 0.01d0, 1, 2,
     &      4.56d-6, 1.32712440018d11, 4.902801076d3,
     &      ALT_a, DEN_a, SCH_a, natm,
     &      reentry, e_opt, bn_opt, rms_o,
     &      zepoch, nzu,
     &      zstat, nzv,
     &      rpe, t_mean, t_std, ierr)

         if (ierr .ne. 0 .and. ierr .ne. 3) then
            write(*,'(A,I3)') '  pipeline error ierr=', ierr
            write(33,'(I5,A,I3)') norad(io), ',ERR,,,,,,,,', ierr
            goto 90
         end if

c        Ensemble RPE: t_mean vs observed, horizon from zone 1 epoch
         call cal2jd(t_obs, t_obs_jd)
         rpe_ens = 0.d0
         nre = 0
         do iz = 1, nzu
            if (reentry(iz) .gt. 0.d0) nre = nre + 1
         end do
         if (nre .ge. 1 .and. t_obs_jd .gt. 0.d0) then
            horiz = t_obs_jd - zepoch(1)
            if (horiz .gt. 0.d0)
     &         rpe_ens = (t_mean - t_obs_jd)/horiz*100.d0
         end if

         do iz = 1, nzu
            if (reentry(iz) .gt. 0.d0) then
               write(*,'(A,I1,A,F8.4,A,F7.2,A,F8.2,A,I2)')
     &            '  Z', iz, ': e=', e_opt(iz), ' BN=', bn_opt(iz),
     &            ' RPE=', rpe(iz), '%  zstat=', zstat(iz)
            else
               write(*,'(A,I1,A,F8.4,A,F7.2,A,I2)')
     &            '  Z', iz, ': e=', e_opt(iz), ' BN=', bn_opt(iz),
     &            ' (no re-entry)  zstat=', zstat(iz)
            end if
            write(33,'(I5,A,I1,A,F9.6,A,F8.3,A,F13.4,A,F9.3,A,I2,
     &         A,F13.4,A,F8.2,A,F9.3)')
     &         norad(io), ',', iz, ',', e_opt(iz), ',', bn_opt(iz),
     &         ',', reentry(iz), ',', rpe(iz), ',', zstat(iz),
     &         ',', t_mean, ',', t_std, ',', rpe_ens
         end do

         if (nre .ge. 1) then
            write(*,'(A,I1,A,I1,A,F13.3,A,F7.1,A)')
     &         '  ensemble (', nre, ' of ', nzu,
     &         ' zones): mean JD=', t_mean, ' std=', t_std, ' d'
            write(*,'(A,F8.2,A)')
     &         '  ENSEMBLE RPE = ', rpe_ens, '%'
         else
            write(*,'(A)') '  no re-entry in any zone (weak-signal'//
     &         ' fits; lifetime beyond 5-yr cap)'
         end if

 90      continue
      end do

      close(33)
      write(*,'(/,A)') ' campaign complete -> scratch_rpe/'//
     &   'rpe_campaign_weather.csv'

      stop
      end

c     ============================================================
c     ATM reader (same as test_reentry.F)
c     ============================================================
      subroutine read_atm_c(alt, den, sch, maxn, natm)
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

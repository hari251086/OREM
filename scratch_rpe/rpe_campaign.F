c     rpe_campaign.F — 10-object drag-enabled RPE campaign (issue #16)
c     -------------------------------------------------------
c     Runs orem_run with IDRAG=1 and the full force model (the
c     test_e2e v1.11 config: geo=20, sun=2, moon=3, SRP on) on all
c     ten validation objects, and reports per-zone RPE plus the
c     per-object ensemble prediction (t_mean/t_std -> ensemble RPE).
c     Purpose: quantify post-v1.18 accuracy across the whole set --
c     test_e2e only exercises 42928 (+z0, 39615 z1, 35497 z2) and
c     test_reentry runs IDRAG=0. Results feed #16 (<10% target) and
c     the #13 ensemble-reporting design.
c
c     Extended from 7 to 10 objects (2026-07-19) with the three real-
c     decay cases from the same day's earlier individual validation
c     runs (11550/59347/40943, input/orem_<norad>.cfg) -- deliberately
c     run through the SAME nzones_max=8 (mxz below) as the original 7,
c     not the nzones_max=50 used for OREM-Watchlist's operational
c     37398 case (issue #9 there): the point of this campaign is to
c     test whether OREM's shipped production default generalizes
c     beyond the curated original set, not to tune per object.
c
c     Compile:
c       ifx /heap-arrays /F:16777216 scratch_rpe\rpe_campaign.F
c           orem.F rsm.F ga.F tle_evolution.F tle_filter.F
c           zone_select.F ksrop\propagate_ks.F ksrop\Subrouts.F
c           ksrop\TLEread.F ksrop\Legendre.F
c           /exe:scratch_rpe\rpe_campaign.exe
c     Run from the OREM repo root:
c       scratch_rpe\rpe_campaign.exe
c     Output: scratch_rpe\rpe_campaign.csv
c     -------------------------------------------------------
      program rpe_campaign

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

c     mxz=8 (v1.21 experiment; was 4): zone_select returns the top-R2
c     candidates, so a higher cap admits LATER zones -- and v1.20
c     showed the latest zone is the primary estimator. This run
c     measures how much later the latest zone gets and what it does
c     to the latest-zone RPE.
      parameter(ndim_a = 291, mxz = 8, nobj = 10)

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
      integer io, iz, nre

c     11550's date (2025 4 3) is INFERRED, not an official decay-date
c     file -- its last TLE shows mean perigee altitude collapsed to
c     96 km, the terminal-signature basis used throughout this
c     session's 11550 work. 59347/40943 dates are real, confirmed
c     (matching their individual input/orem_<norad>.cfg files).
      data norad / 42928, 35497, 37151, 39615, 27526, 32007, 37819,
     &             11550, 59347, 40943 /
      data oyr / 2019.d0, 2016.d0, 2015.d0, 2017.d0, 2012.d0,
     &           2010.d0, 2013.d0, 2025.d0, 2026.d0, 2020.d0 /
      data omo / 3.d0, 10.d0, 12.d0, 9.d0, 5.d0, 6.d0, 9.d0,
     &           4.d0, 6.d0, 2.d0 /
      data ody / 3.d0, 31.d0, 3.d0, 15.d0, 9.d0, 6.d0, 12.d0,
     &           3.d0, 8.d0, 20.d0 /

      tle_file(1) = 'input/example_42928.tle.txt'
      tle_file(2) = 'input/example_35497.tle.txt'
      tle_file(3) = 'input/example_37151.tle.txt'
      tle_file(4) = 'input/example_39615.tle.txt'
      tle_file(5) = 'input/example_27526.tle.txt'
      tle_file(6) = 'input/example_32007.tle.txt'
      tle_file(7) = 'input/example_37819.tle.txt'
      tle_file(8) = 'input/example_11550.tle.txt'
      tle_file(9) = 'input/example_59347.tle.txt'
      tle_file(10) = 'input/example_40943.tle.txt'
      obj_name(1) = '42928 PSLV-C39'
      obj_name(2) = '35497 Ariane 5 ESC-A'
      obj_name(3) = '37151 Long March 3B'
      obj_name(4) = '39615 Proton-M Briz-M'
      obj_name(5) = '27526 Ariane 5 R/B'
      obj_name(6) = '32007 GSLV R/B'
      obj_name(7) = '37819 Proton-M R/B'
      obj_name(8) = '11550 Molniya debris'
      obj_name(9) = '59347 GTO R/B'
      obj_name(10) = '40943 GTO debris'

      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      open(unit=33, file='scratch_rpe/rpe_campaign.csv',
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
     &   'rpe_campaign.csv'

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

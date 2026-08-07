c     ============================================================
c     two_pass_check.F -- issue #36 Finding 3 follow-up: measures
c     whether the two-pass RSM/GA refinement (orem_run_2pass.F)
c     actually improves real ensemble RPE, not just the interp-vs-
c     direct gap.
c
c     Runs the curated-7 under the same production parameters as
c     scratch_rpe/rpe_campaign.F (nzones_max=8, IDRAG=1, full force
c     model geo=20/sun=2/moon=3, maxgen=500, max_zone_days=20) --
c     current shipped baseline for the same objects/params is in
c     scratch_rpe/rpe_campaign.csv (mean|ens_rpe|=1.94%, median=1.19%,
c     2026-08-01 regeneration).
c
c     Writes two files:
c       scratch_rpe/bench/twopass/rms_gap_results.csv -- per zone,
c         pass-1 vs pass-2 rms_fit and interp-vs-direct RMS
c       scratch_rpe/bench/twopass/object_results.csv -- per-object
c         t_mean/t_std/ens_rpe_pct (compare directly against
c         rpe_campaign.csv's baseline rows)
c
c     Not wired into fpm.toml (ad-hoc). To reproduce, temporarily add:
c       [[executable]]
c       name = "two_pass_check"
c       source-dir = "scratch_rpe/bench/twopass"
c       main = "two_pass_check.F"
c     then: fpm run two_pass_check --compiler ifx --flag "-heap-arrays"
c
c     RESULT (2026-08-01, ifx, this machine): REJECTED. mean|ens_rpe|
c     1.94%->9.30% (4.8x worse), median 1.19%->8.12% (6.8x worse); 4/7
c     objects regress sharply, only 2 improve marginally. 5th confirmed
c     instance of the BN-carryover chain amplifying a locally-motivated
c     mechanism into an aggregate regression (see issue #35). NOT
c     shipped -- production orem.F is unchanged.
c
c     KNOWN BUG, not chased further since it doesn't affect the RPE
c     conclusion above: rms_gap_results.csv's rms_direct_pass2 column
c     is bogus (thousands of km, roughly matching raw apogee
c     magnitude -- signature of the diagnostic's own extra
c     propagate_ks call silently failing, unrelated to the real
c     downstream pipeline which is copied verbatim from production and
c     produces the sane t_mean/t_std/ens_rpe_pct values this file's
c     real conclusion is based on). Root cause not isolated -- if this
c     two-pass idea is ever revisited, fix that diagnostic block first
c     rather than trusting rms_gap_results.csv's direct-RMS column.
c     ============================================================
      program two_pass_check
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      character*30 f3_obj_name
      common /f3_label/ f3_obj_name

      parameter(ndim_a = 291, mxz = 8)
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      open(unit=93, file=
     &   'scratch_rpe/bench/twopass/rms_gap_results.csv',
     &   status='replace')
      write(93,'(A)') 'object,zone,nzone,bn_span_pass1,bn_lo2,'//
     &   'bn_hi2,rms_pass1,rms_pass2,rms_interp_pass2,'//
     &   'rms_direct_pass2'

      open(unit=94, file=
     &   'scratch_rpe/bench/twopass/object_results.csv',
     &   status='replace')
      write(94,'(A)') 'object,ierr,nzones,t_mean,t_std,ens_rpe_pct'

      write(*,*) '=============================================='
      write(*,*) ' Two-pass RSM/GA refinement -- real RPE check'
      write(*,*) ' Curated-7, production params (IDRAG=1, full'
      write(*,*) ' force model, maxgen=500, nzones_max=8)'
      write(*,*) '=============================================='

      call tp_case('input/example_42928.tle.txt', 42928,
     &   '42928', 2019.d0, 3.d0, 3.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call tp_case('input/example_35497.tle.txt', 35497,
     &   '35497', 2016.d0, 10.d0, 31.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call tp_case('input/example_37151.tle.txt', 37151,
     &   '37151', 2015.d0, 12.d0, 3.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call tp_case('input/example_39615.tle.txt', 39615,
     &   '39615', 2017.d0, 9.d0, 15.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call tp_case('input/example_27526.tle.txt', 27526,
     &   '27526', 2012.d0, 5.d0, 9.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call tp_case('input/example_32007.tle.txt', 32007,
     &   '32007', 2010.d0, 6.d0, 6.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call tp_case('input/example_37819.tle.txt', 37819,
     &   '37819', 2013.d0, 9.d0, 12.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      close(93)
      close(94)

      write(*,*) ''
      write(*,*) '=============================================='
      write(*,*) ' Done -- see rms_gap_results.csv, '//
     &   'object_results.csv'
      write(*,*) '=============================================='

      stop
      end

c     ============================================================
      subroutine tp_case(tle_file, norad, obj_name,
     &   obs_yr, obs_mo, obs_dy,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      character*30 f3_obj_name
      common /f3_label/ f3_obj_name

      character*(*) tle_file, obj_name
      integer norad, natm, ndim_a, mxz
      double precision obs_yr, obs_mo, obs_dy
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)

      double precision reentry(8), e_opt(8), bn_opt(8)
      double precision rms_o(8), zepoch(8), rpe(8)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr, nzv
      integer zstat(8)
      double precision ens_rpe_pct, t_obs_jd_loc, horiz
      integer nre, iz

      f3_obj_name = obj_name

      t_obs(1) = obs_yr
      t_obs(2) = obs_mo
      t_obs(3) = obs_dy
      t_obs(4) = 0.d0
      t_obs(5) = 0.d0
      t_obs(6) = 0.d0

      write(*,*) ''
      write(*,'(A,A)') '  --- ', obj_name

      call orem_run_2pass(
     &   tle_file, norad,
     &   t_obs, mxz, 8, 20.d0, 0.90d0, -1.d0,
     &   80.d0, 160.d0, 1,
     &   20, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   20, 2, 3,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   reentry, e_opt, bn_opt, rms_o,
     &   zepoch, nzu,
     &   zstat, nzv,
     &   rpe, t_mean, t_std, ierr)

c     Ensemble RPE: same formula as rpe_campaign.F -- t_mean vs
c     observed, normalized by horizon from zone 1 epoch (not by
c     t_obs_jd itself).
      ens_rpe_pct = 0.d0
      call cal2jd(t_obs, t_obs_jd_loc)
      nre = 0
      do iz = 1, nzu
         if (reentry(iz) .gt. 0.d0) nre = nre + 1
      end do
      if (nre .ge. 1 .and. t_obs_jd_loc .gt. 0.d0) then
         horiz = t_obs_jd_loc - zepoch(1)
         if (horiz .gt. 0.d0)
     &      ens_rpe_pct = (t_mean - t_obs_jd_loc) / horiz * 100.d0
      end if

      write(*,'(A,I2,A,I2,A,F14.4,A,F8.2,A,F8.3)')
     &   '      ierr=', ierr, ' zones=', nzu,
     &   ' t_mean=', t_mean, ' t_std=', t_std,
     &   ' ens_rpe_pct=', ens_rpe_pct

      write(94,930) obj_name, ierr, nzu, t_mean, t_std, ens_rpe_pct
 930  format(A,',',I2,',',I2,',',F14.4,',',F8.2,',',F8.3)

      return
      end

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

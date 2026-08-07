c     ============================================================
c     finding3_check.F -- issue #36 Finding 3 generalization check.
c
c     Runs orem_run_finding3 (see that file's header) on the curated-7
c     validation objects under the same production parameters
c     rpe_campaign.F uses (nzones_max=8, IDRAG=1, full force model
c     geo=20/sun=2/moon=3, maxgen=500, max_zone_days=20). Appends one
c     CSV row per zone to scratch_rpe/bench/finding3_results.csv:
c     object,zone,nzone,zone_dur_days,e_lo,e_hi,bn_lo,bn_hi,e_opt,
c     bn_opt,rms_fit_reported,rms_interp,rms_direct,ntraj_pts
c
c     Not wired into fpm.toml (ad-hoc, matching benchmark_stages.F's
c     pattern). To reproduce, temporarily add to fpm.toml:
c       [[executable]]
c       name = "finding3_check"
c       source-dir = "scratch_rpe/bench/finding3"
c       main = "finding3_check.F"
c     then: fpm run finding3_check --compiler ifx --flag "-heap-arrays"
c
c     Results (2026-08-01, ifx, this machine): 56 real zones, median
c     ratio 1.79x, mean 5.6x (right-skewed, up to 55.6x), dominant
c     correlate is BN grid span (corr=0.71 on log-ratio) -- see
c     README.md Version History / issue #32 for the full writeup.
c     Follow-up (two-pass BN/e refinement, scratch_rpe/bench/twopass/)
c     was tested and REJECTED -- regresses real ensemble RPE 4.8x via
c     the same BN-carryover-chain fragility as every other locally-
c     motivated fix in this investigation. Don't re-propose without a
c     genuinely different mechanism for feeding a corrected fit back
c     into the pipeline.
c     ============================================================
      program finding3_check
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      character*30 f3_obj_name
      common /f3_label/ f3_obj_name

      parameter(ndim_a = 291, mxz = 8)
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr, nzv
      integer zstat(mxz)

      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      open(unit=92, file='scratch_rpe/bench/finding3_results.csv',
     &     status='replace')
      write(92,'(A)') 'object,zone,nzone,zone_dur_days,e_lo,e_hi,'//
     &   'bn_lo,bn_hi,e_opt,bn_opt,rms_fit_reported,rms_interp,'//
     &   'rms_direct,ntraj_pts'

      write(*,*) '=============================================='
      write(*,*) ' Issue #36 Finding 3 generalization check'
      write(*,*) ' Curated-7, production params (IDRAG=1, full'
      write(*,*) ' force model, maxgen=500, nzones_max=8)'
      write(*,*) '=============================================='

      call f3_case('input/example_42928.tle.txt', 42928,
     &   '42928', 2019.d0, 3.d0, 3.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call f3_case('input/example_35497.tle.txt', 35497,
     &   '35497', 2016.d0, 10.d0, 31.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call f3_case('input/example_37151.tle.txt', 37151,
     &   '37151', 2015.d0, 12.d0, 3.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call f3_case('input/example_39615.tle.txt', 39615,
     &   '39615', 2017.d0, 9.d0, 15.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call f3_case('input/example_27526.tle.txt', 27526,
     &   '27526', 2012.d0, 5.d0, 9.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call f3_case('input/example_32007.tle.txt', 32007,
     &   '32007', 2010.d0, 6.d0, 6.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call f3_case('input/example_37819.tle.txt', 37819,
     &   '37819', 2013.d0, 9.d0, 12.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      close(92)

      write(*,*) ''
      write(*,*) '=============================================='
      write(*,*) ' Done -- see finding3_results.csv'
      write(*,*) '=============================================='

      stop
      end

c     ============================================================
      subroutine f3_case(tle_file, norad, obj_name,
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

      f3_obj_name = obj_name

      t_obs(1) = obs_yr
      t_obs(2) = obs_mo
      t_obs(3) = obs_dy
      t_obs(4) = 0.d0
      t_obs(5) = 0.d0
      t_obs(6) = 0.d0

      write(*,*) ''
      write(*,'(A,A)') '  --- ', obj_name

      call orem_run_finding3(
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

      write(*,'(A,I2,A,I2)') '      ierr=', ierr, ' zones=', nzu

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

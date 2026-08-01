c     ============================================================
c     benchmark_stages.F -- issue #17 (performance benchmarking).
c
c     Runs orem_run_timed (see orem_run_timed.F) on the 4 original
c     curated objects (42928, 35497, 37151, 39615) with the same
c     production parameters rpe_campaign.F uses (mxz=8, IDRAG=1, full
c     force model geo=20/sun=2/moon=3, maxgen=500) -- not
c     test_reentry.F's stripped IDRAG=0/maxgen=200 config, since that
c     config never exercises re-entry propagation at all and #17 asks
c     specifically about it. Reports wall-clock time per stage
c     (TLE load, zone select, RSM, GA, re-entry propagation) and per
c     object, plus a total-per-zone figure. Results (2026-08-01,
c     ifx, this machine): see README.md Version History / #17.
c
c     Not wired into fpm.toml (ad-hoc like the rest of scratch_rpe/,
c     not a permanent build target). To reproduce, temporarily add:
c       [[executable]]
c       name = "benchmark_stages"
c       source-dir = "scratch_rpe/bench"
c       main = "benchmark_stages.F"
c     to fpm.toml, then: fpm run benchmark_stages --compiler ifx
c     --flag "-heap-arrays"
c     ============================================================
      program benchmark_stages
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision t_tle, t_zone, t_rsm, t_ga, t_reentry
      integer nz_rsm
      common /bench_time/ t_tle, t_zone, t_rsm, t_ga, t_reentry, nz_rsm

      parameter(ndim_a = 291, mxz = 8)
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr, nzv
      integer zstat(mxz)

      character*60 tle_file
      character*30 obj_name
      integer norad
      integer(8) wc0, wc1, wcrate
      double precision t_wall

      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)
      call system_clock(count_rate=wcrate)

      write(*,*) '=============================================='
      write(*,*) ' OREM Stage Benchmark (issue #17)'
      write(*,*) ' 4 curated objects, production params (IDRAG=1,'
      write(*,*) ' full force model, maxgen=500, nzones_max=8)'
      write(*,*) '=============================================='

      call bench_case('input/example_42928.tle.txt', 42928,
     &   '42928 PSLV-C39',
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call bench_case('input/example_35497.tle.txt', 35497,
     &   '35497 Ariane 5 ESC-A',
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call bench_case('input/example_37151.tle.txt', 37151,
     &   '37151 Long March 3B',
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call bench_case('input/example_39615.tle.txt', 39615,
     &   '39615 Proton-M Briz-M',
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      write(*,*) ''
      write(*,*) '=============================================='
      write(*,*) ' Done'
      write(*,*) '=============================================='

      stop
      end

c     ============================================================
      subroutine bench_case(tle_file, norad, obj_name,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision t_tle, t_zone, t_rsm, t_ga, t_reentry
      integer nz_rsm
      common /bench_time/ t_tle, t_zone, t_rsm, t_ga, t_reentry, nz_rsm

      character*(*) tle_file, obj_name
      integer norad, natm, ndim_a, mxz
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)

      double precision reentry(8), e_opt(8), bn_opt(8)
      double precision rms_o(8), zepoch(8), rpe(8)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr, nzv
      integer zstat(8)
      integer(8) wc0, wc1, wcrate
      double precision t_wall, t_sum

      t_tle = 0.d0
      t_zone = 0.d0
      t_rsm = 0.d0
      t_ga = 0.d0
      t_reentry = 0.d0
      nz_rsm = 0

      do i = 1, 6
         t_obs(i) = 0.d0
      end do

      call system_clock(count_rate=wcrate)
      call system_clock(wc0)

      call orem_run_timed(
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

      call system_clock(wc1)
      t_wall = dble(wc1 - wc0) / dble(wcrate)
      t_sum = t_tle + t_zone + t_rsm + t_ga + t_reentry

      write(*,*) ''
      write(*,'(A,A,A,I2,A,I2)') '  --- ', obj_name,
     &   ' (ierr=', ierr, ', zones=', nzu
      write(*,'(A,F8.2,A)') '  TLE load+filter+read : ', t_tle, ' s'
      write(*,'(A,F8.2,A)') '  Zone selection        : ', t_zone, ' s'
      write(*,'(A,F8.2,A,I2,A,F8.2,A)') '  RSM (9x propagate_ks) : ',
     &   t_rsm, ' s  (', nz_rsm, ' zones, ',
     &   t_rsm/dmax1(1.d0,dble(nz_rsm)), ' s/zone)'
      write(*,'(A,F8.2,A,I2,A,F8.2,A)') '  GA (500gen x 20pop)   : ',
     &   t_ga, ' s  (', nz_rsm, ' zones, ',
     &   t_ga/dmax1(1.d0,dble(nz_rsm)), ' s/zone)'
      write(*,'(A,F8.2,A)') '  Re-entry propagation  : ',
     &   t_reentry, ' s'
      write(*,'(A,F8.2,A)') '  Sum of stages         : ', t_sum, ' s'
      write(*,'(A,F8.2,A)') '  Wall clock (full call): ', t_wall, ' s'
      write(*,'(A,F8.2,A)') '  Per zone (wall/zones) : ',
     &   t_wall/dmax1(1.d0,dble(nzu)), ' s'

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

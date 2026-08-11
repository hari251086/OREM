c     ============================================================
c     dumpsurf_check.F -- issue #40 option 3: dump real per-zone RSM
c     landscapes for offline Bayesian-optimization validation.
c
c     Runs orem_run_dumpsurf (see that file's header) on the curated-7
c     validation objects under the same production parameters
c     finding3_check.F/gasweep_check.F use. Writes
c     scratch_rpe/bench/surfaces_dump.csv:
c     object,zone,it,tobs,apobs,e1,e2,e3,bn1,bn2,bn3,rms_baseline,
c     s11,s21,s31,s12,s22,s32,s13,s23,s33
c
c     Not wired into fpm.toml (ad-hoc, matching finding3_check.F's
c     pattern). To reproduce, temporarily add to fpm.toml:
c       [[executable]]
c       name = "dumpsurf_check"
c       source-dir = "scratch_rpe/bench/gasweep"
c       main = "dumpsurf_check.F"
c     then: fpm run dumpsurf_check --compiler ifx --flag "-heap-arrays"
c     ============================================================
      program dumpsurf_check
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      character*30 ds_obj_name
      common /ds_label/ ds_obj_name

      parameter(ndim_a = 291, mxz = 8)
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      open(unit=92, file='scratch_rpe/bench/surfaces_dump.csv',
     &     status='replace')
      write(92,'(A)') 'object,zone,it,tobs,apobs,'//
     &   'e1,e2,e3,bn1,bn2,bn3,rms_baseline,'//
     &   's11,s21,s31,s12,s22,s32,s13,s23,s33'

      write(*,*) '=============================================='
      write(*,*) ' Issue #40 option 3: RSM landscape dump'
      write(*,*) ' Curated-7, production params'
      write(*,*) '=============================================='

      call ds_case('input/example_42928.tle.txt', 42928,
     &   '42928', 2019.d0, 3.d0, 3.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call ds_case('input/example_35497.tle.txt', 35497,
     &   '35497', 2016.d0, 10.d0, 31.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call ds_case('input/example_37151.tle.txt', 37151,
     &   '37151', 2015.d0, 12.d0, 3.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call ds_case('input/example_39615.tle.txt', 39615,
     &   '39615', 2017.d0, 9.d0, 15.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call ds_case('input/example_27526.tle.txt', 27526,
     &   '27526', 2012.d0, 5.d0, 9.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call ds_case('input/example_32007.tle.txt', 32007,
     &   '32007', 2010.d0, 6.d0, 6.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      call ds_case('input/example_37819.tle.txt', 37819,
     &   '37819', 2013.d0, 9.d0, 12.d0,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)

      close(92)

      write(*,*) ''
      write(*,*) '=============================================='
      write(*,*) ' Done -- see surfaces_dump.csv'
      write(*,*) '=============================================='

      stop
      end

c     ============================================================
      subroutine ds_case(tle_file, norad, obj_name,
     &   obs_yr, obs_mo, obs_dy,
     &   ALT_a, DEN_a, SCH_a, natm, ndim_a, mxz)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      character*30 ds_obj_name
      common /ds_label/ ds_obj_name

      character*(*) tle_file, obj_name
      integer norad, natm, ndim_a, mxz
      double precision obs_yr, obs_mo, obs_dy
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)

      double precision reentry(8), e_opt(8), bn_opt(8)
      double precision rms_o(8), zepoch(8), rpe(8)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr, nzv
      integer zstat(8)

      ds_obj_name = obj_name

      t_obs(1) = obs_yr
      t_obs(2) = obs_mo
      t_obs(3) = obs_dy
      t_obs(4) = 0.d0
      t_obs(5) = 0.d0
      t_obs(6) = 0.d0

      write(*,*) ''
      write(*,'(A,A)') '  --- ', obj_name

      call orem_run_dumpsurf(
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

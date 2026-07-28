c     diag_61734_zone1.F -- issue #29 follow-up: root-cause 61734's own
c     zone-1 defect (BN=16.286, boundary-pinned, RPE=759.7%), flagged
c     repeatedly this session (v1.41/v1.42/v1.43) as a DIFFERENT,
c     upstream problem from the BN-carryover chain -- zone 1 has no
c     prior zone to inherit a bad range from, so whatever is wrong here
c     is intrinsic to zone 1's own RSM/GA fit, not carryover.
c
c     Single-object re-run of the real 50-object campaign's exact call
c     (scratch_rpe/rpe_campaign.F, object index 13, same config values)
c     so all of orem.F's own console diagnostics (G2 BN floor, G3
c     BSTAR-informed prior, per-zone RSM/GA reports) print for this one
c     object without the other 49 objects' noise in between.
c
c     Compile:
c       ifx /heap-arrays /F:16777216 scratch_rpe\diag_61734_zone1.F
c           src\orem.F src\rsm.F src\ga.F src\tle_evolution.F
c           src\zone_select.F src\propagate_ks.F src\swx.F
c           build\dependencies\ksrop\src\Subrouts.F
c           build\dependencies\ksrop\src\Legendre.F
c           build\dependencies\ksrop\src\TLEread.F
c           /exe:scratch_rpe\diag_61734_zone1.exe
c
c     Analysis only -- not part of the build/test suite.
      program diag_61734_zone1
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291, mxz = 8)

      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr, nzv
      integer zstat(mxz)
      double precision t_obs_jd, horiz, rpe_ens
      integer iz

      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

c     61734 H-3 R/B -- real observed re-entry 2025-03-22
      t_obs(1) = 2025.d0
      t_obs(2) = 3.d0
      t_obs(3) = 22.d0
      t_obs(4) = 0.d0
      t_obs(5) = 0.d0
      t_obs(6) = 0.d0

      call orem_run(
     &   'input/example_61734.tle.txt', 61734,
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

      write(*,'(A,I3)') ' ierr=', ierr
      call cal2jd(t_obs, t_obs_jd)
      write(*,'(A,F14.4)') ' t_obs_jd=', t_obs_jd

      do iz = 1, nzu
         write(*,'(A,I2,A,F14.4,A,F9.6,A,F8.3,A,F14.4,A,F8.2,A,I2)')
     &      ' Z', iz, ': zepoch=', zepoch(iz),
     &      '  e=', e_opt(iz), '  BN=', bn_opt(iz),
     &      '  reentry_jd=', reentry(iz),
     &      '  RPE=', rpe(iz), '%  zstat=', zstat(iz)
      end do

      write(*,'(A,F14.4,A,F8.2)') ' ensemble t_mean=', t_mean,
     &   '  t_std=', t_std

      stop
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

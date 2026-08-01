c     ============================================================
c     bn_scale_test.F -- sanity check: does decay rate scale as
c     1/BN as the drag formula predicts? Same IC as ../npoe_xval.F,
c     short 10-day window, three BN values.
c
c     Not wired into fpm.toml. To reproduce, temporarily add:
c       [[executable]]
c       name = "bn_scale_test"
c       source-dir = "scratch_rpe/bench/npoe_xval/bnscale"
c       main = "bn_scale_test.F"
c     then: fpm run bn_scale_test --compiler ifx --flag "-heap-arrays"
c
c     RESULT (2026-08-01): BN=112.11 -> 19.46km/10d, BN=11.21 ->
c     192.97km/10d (ratio 9.92, expected ~10.0) -- confirms the drag
c     mechanism is correctly wired, not the source of the much larger
c     gap found in ../npoe_xval.F (see that file's header / issue #38).
c     ============================================================
      program bn_scale_test
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291)
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      integer n_force(3)

      parameter(max_pts = 200)
      double precision traj_jd(max_pts)
      double precision traj_x(3, max_pts)
      double precision traj_xd(3, max_pts)
      integer exit_code, idump

      double precision bn_test(3)
      double precision sma0, sma_end, r_mag, v_mag
      integer ibn

      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      a_kep(1) = 9594.412450d0
      a_kep(2) = 0.3185937d0
      a_kep(3) = 19.1779d0
      a_kep(4) = 36.853002d0
      a_kep(5) = 332.445197d0
      a_kep(6) = 96.766670d0
      a_kep(7) = 0.d0

      cal0(1) = 2017.d0
      cal0(2) = 10.d0
      cal0(3) = 3.d0
      cal0(4) = 14.d0
      cal0(5) = 35.d0
      cal0(6) = 56.d0

      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

      bn_test(1) = 112.1098d0
      bn_test(2) = 11.21098d0
      bn_test(3) = 8.0d0

      do ibn = 1, 3
         call oe2car(a_kep, x0, xd0, 1.d-12)
         sma0 = a_kep(1)

         idump = 1
         call propagate_ks(
     &      x0, xd0, cal0,
     &      92, 360, 1.d-10,
     &      n_force, 6, 2, 3,
     &      bn_test(ibn), 1, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &      1.2d0, 0.01d0, 0, 0,
     &      4.56d-6, 1.32712440018d11, 4.902801076d3,
     &      ALT_a, DEN_a, SCH_a, natm,
     &      max_pts, idump,
     &      traj_jd, traj_x, traj_xd,
     &      exit_code)

         r_mag = dsqrt(traj_x(1,idump)**2 + traj_x(2,idump)**2
     &               + traj_x(3,idump)**2)
         v_mag = dsqrt(traj_xd(1,idump)**2 + traj_xd(2,idump)**2
     &               + traj_xd(3,idump)**2)
         sma_end = 1.d0 / (2.d0/r_mag - v_mag**2/amue)

         write(*,'(A,F10.4,A,I2,A,F10.4,A,F10.4,A,F12.6,A,F12.6)')
     &      'BN=', bn_test(ibn), ' exit=', exit_code,
     &      ' idump=', dble(idump), ' sma0=', sma0,
     &      ' sma_end(10d)=', sma_end, ' decay=', sma0-sma_end
      end do

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

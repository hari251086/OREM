c     test_ga.F â€” Tests for GA optimizer
c     -------------------------------------------------------
c     Compile:
c       ifx test_ga.F ga.F /exe:test_ga.exe
c     -------------------------------------------------------
      program test_ga

      implicit double precision (a-h, o-z)

      integer npass, nfail, iall_ok

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' GA Optimizer Test Suite'
      write(*,*) '=============================================='

c     ==============================================================
c     TWOINT bilinear interpolation tests
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- TWOINT interpolation ---'

      call test_twoint(npass, nfail)

c     ==============================================================
c     GA decode tests
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Chromosome decode ---'

      call test_decode(npass, nfail)

c     ==============================================================
c     GA RNG tests
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Random number generator ---'

      call test_rng(npass, nfail)

c     ==============================================================
c     GA optimization â€” synthetic known-optimum
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Synthetic optimization ---'

      call test_synthetic(npass, nfail)

c     ==============================================================
c     GA robustness
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Robustness ---'

      call test_robustness(npass, nfail)

c     ==============================================================
c     Piecewise: ga_tinterp
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- ga_tinterp ---'

      call test_tinterp(npass, nfail)

c     ==============================================================
c     Extended TWOINT edge cases
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- TWOINT edge cases ---'

      call test_twoint_edge(npass, nfail)

c     ==============================================================
c     Extended decode edge cases
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Decode edge cases ---'

      call test_decode_edge(npass, nfail)

c     ==============================================================
c     GA convergence properties
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Convergence properties ---'

      call test_convergence(npass, nfail)

c     ==============================================================
c     GA edge-case surfaces
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Surface edge cases ---'

      call test_surface_edge(npass, nfail)

c     ==============================================================
c     High eccentricity orbits
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- High eccentricity orbits ---'

      call test_high_ecc(npass, nfail)

c     ==============================================================
c     Piecewise GA internals
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- GA internals ---'

      call test_ga_internals(npass, nfail)

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
c     TWOINT tests
c     ============================================================
      subroutine test_twoint(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

      double precision x3(3), y3(3), f3(3,3), fout

c     G1: Constant surface â†’ interpolation returns constant
      x3(1) = 0.3d0
      x3(2) = 0.32d0
      x3(3) = 0.34d0
      y3(1) = 80.d0
      y3(2) = 120.d0
      y3(3) = 160.d0
      do i = 1, 3
         do j = 1, 3
            f3(i,j) = 42.d0
         end do
      end do
      call ga_twoint(x3, y3, f3, 0.31d0, 100.d0, fout, 3, 3)
      call assert_d('G1 twoint const=42',
     &              fout, 42.d0, 1.d-10, npass, nfail)

c     G2: Linear surface f = 10*x + y
      do i = 1, 3
         do j = 1, 3
            f3(i,j) = 10.d0*x3(i) + y3(j)
         end do
      end do
      call ga_twoint(x3, y3, f3, 0.31d0, 100.d0,fout, 3, 3)
      call assert_d('G2 twoint linear',
     &              fout, 10.d0*0.31d0+100.d0,
     &              1.d-10, npass, nfail)

c     G3: Corner values exact
      call ga_twoint(x3, y3, f3, 0.3d0, 80.d0, fout, 3, 3)
      call assert_d('G3 twoint corner(1,1)',
     &              fout, 10.d0*0.3d0+80.d0,
     &              1.d-10, npass, nfail)

c     G4: Center point
      call ga_twoint(x3, y3, f3, 0.32d0, 120.d0,fout, 3, 3)
      call assert_d('G4 twoint center',
     &              fout, 10.d0*0.32d0+120.d0,
     &              1.d-10, npass, nfail)

c     G5: Edge value
      call ga_twoint(x3, y3, f3, 0.34d0, 160.d0,fout, 3, 3)
      call assert_d('G5 twoint edge(3,3)',
     &              fout, 10.d0*0.34d0+160.d0,
     &              1.d-10, npass, nfail)

      return
      end

c     ============================================================
c     Decode tests
c     ============================================================
      subroutine test_decode(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

      integer chrom(80)
      double precision xvec(2), alow(2), ahigh(2), factor(2)
      integer lsubstr(2)

      lsubstr(1) = 40
      lsubstr(2) = 40
      alow(1) = 0.d0
      ahigh(1) = 1.d0
      alow(2) = 0.d0
      ahigh(2) = 1.d0
      factor(1) = 2.d0**40 - 1.d0
      factor(2) = 2.d0**40 - 1.d0

c     G6: All zeros â†’ both variables = lower bound
      do i = 1, 80
         chrom(i) = 0
      end do
      call ga_decode(chrom, 80, xvec, 2,
     &               lsubstr, alow, ahigh, factor)
      call assert_d('G6 all-0 chr â†’ x=0',
     &              xvec(1), 0.d0, 1.d-10, npass, nfail)
      call assert_d('G7 all-0 chr â†’ y=0',
     &              xvec(2), 0.d0, 1.d-10, npass, nfail)

c     G8: All ones â†’ both variables = upper bound
      do i = 1, 80
         chrom(i) = 1
      end do
      call ga_decode(chrom, 80, xvec, 2,
     &               lsubstr, alow, ahigh, factor)
      call assert_d('G8 all-1 chr â†’ x=1',
     &              xvec(1), 1.d0, 1.d-10, npass, nfail)
      call assert_d('G9 all-1 chr â†’ y=1',
     &              xvec(2), 1.d0, 1.d-10, npass, nfail)

c     G10: First bit only â†’ x = 1/(2^40 - 1) â‰ˆ 0
      do i = 1, 80
         chrom(i) = 0
      end do
      chrom(1) = 1
      call ga_decode(chrom, 80, xvec, 2,
     &               lsubstr, alow, ahigh, factor)
      if (xvec(1) .gt. 0.d0 .and. xvec(1) .lt. 1.d-6) then
         call assert_i('G10 single bit â†’ small x',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G10 single bit â†’ small x',
     &                 0, 1, npass, nfail)
      end if

      return
      end

c     ============================================================
c     RNG tests
c     ============================================================
      subroutine test_rng(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

      double precision oldrand(55), oldrand2(55), rval
      integer jrand, jrand2, i, iall_ok

c     G11: All random numbers in [0, 1)
      call ga_warmup(0.123d0, oldrand, jrand)
      iall_ok = 1
      do i = 1, 1000
         rval = ga_random(oldrand, jrand)
         if (rval .lt. 0.d0 .or. rval .ge. 1.d0) iall_ok = 0
      end do
      call assert_i('G11 1000 randoms in [0,1)',
     &              iall_ok, 1, npass, nfail)

c     G12: Different seeds produce different sequences
      call ga_warmup(0.123d0, oldrand, jrand)
      r1 = ga_random(oldrand, jrand)
      r2 = ga_random(oldrand, jrand)

      call ga_warmup(0.987d0, oldrand2, jrand2)
      r3 = ga_random(oldrand2, jrand2)
      r4 = ga_random(oldrand2, jrand2)

      if (dabs(r1 - r3) .gt. 1.d-10) then
         call assert_i('G12 diff seeds â†’ diff seq',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G12 diff seeds â†’ diff seq',
     &                 0, 1, npass, nfail)
      end if

c     G13: Same seed reproduces sequence
      call ga_warmup(0.123d0, oldrand, jrand)
      r5 = ga_random(oldrand, jrand)
      call assert_d('G13 same seed reproducible',
     &              r1, r5, 1.d-15, npass, nfail)

      return
      end

c     ============================================================
c     Synthetic optimization tests
c     ============================================================
      subroutine test_synthetic(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

c     Create a synthetic 3Ã—3 surface grid where the optimal
c     (e, A) = (0.32, 120) produces apogee matching "observations"

      parameter(ne = 3, na = 3, npts = 10, nob = 10)
      double precision surf(npts, ne, na)
      double precision eg(ne), ag(na)
      double precision tob(nob), apob(nob)
      double precision e_opt, a_opt, rms_opt

c     Grid: e = [0.30, 0.32, 0.34], A = [80, 120, 160]
      eg(1) = 0.30d0
      eg(2) = 0.32d0
      eg(3) = 0.34d0
      ag(1) = 80.d0
      ag(2) = 120.d0
      ag(3) = 160.d0

c     Each surface: apogee = 12000 - 100*t - 50*(e-0.32) - 0.5*(A-120)
c     So (e=0.32, A=120) gives the "true" apogee
      do k = 1, na
         do j = 1, ne
            do i = 1, npts
               t = dble(i-1) / dble(npts-1)
               surf(i, j, k) = 12000.d0 - 100.d0*t
     &              - 50.d0*(eg(j) - 0.32d0)
     &              - 0.5d0*(ag(k) - 120.d0)
            end do
         end do
      end do

c     "Observations" = exact match at (0.32, 120)
      do i = 1, nob
         tob(i) = dble(i-1) / dble(npts-1)
         apob(i) = 12000.d0 - 100.d0*tob(i)
      end do

c     G14: GA finds optimal e near 0.32
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)

      if (dabs(e_opt - 0.32d0) .lt. 0.005d0) then
         call assert_i('G14 e_opt near 0.32',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G14 e_opt near 0.32',
     &                 0, 1, npass, nfail)
      end if

c     G15: GA finds optimal A near 120
      if (dabs(a_opt - 120.d0) .lt. 25.d0) then
         call assert_i('G15 A_opt near 120',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G15 A_opt near 120',
     &                 0, 1, npass, nfail)
      end if

c     G16: RMS at optimum is near zero
      if (rms_opt .lt. 0.1d0) then
         call assert_i('G16 rms_opt near 0',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G16 rms_opt near 0',
     &                 0, 1, npass, nfail)
      end if

c     G17: Fitness > 0.9 (high quality match)
      fit = 1.d0 / (1.d0 + rms_opt)
      if (fit .gt. 0.9d0) then
         call assert_i('G17 fitness > 0.9',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G17 fitness > 0.9',
     &                 0, 1, npass, nfail)
      end if

c     G18: Different seed still converges
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.987d0,
     &   e_opt, a_opt, rms_opt)
      if (dabs(e_opt - 0.32d0) .lt. 0.01d0 .and.
     &    rms_opt .lt. 0.5d0) then
         call assert_i('G18 diff seed same optimum',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G18 diff seed same optimum',
     &                 0, 1, npass, nfail)
      end if

c     G19: Wider bounds still finds optimum
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.28d0, 0.36d0, 60.d0, 180.d0,
     &   4, 500, 40, 20, 0.8d0, 0.01d0, 0.456d0,
     &   e_opt, a_opt, rms_opt)
      if (dabs(e_opt - 0.32d0) .lt. 0.01d0 .and.
     &    dabs(a_opt - 120.d0) .lt. 20.d0) then
         call assert_i('G19 wide bounds converges',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G19 wide bounds converges',
     &                 0, 1, npass, nfail)
      end if

c     G20: pop=4 with more generations (1000)
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 1000, 40, 16, 0.8d0, 0.01d0, 0.555d0,
     &   e_opt, a_opt, rms_opt)
      if (dabs(e_opt - 0.32d0) .lt. 0.01d0 .and.
     &    rms_opt .lt. 0.5d0) then
         call assert_i('G20 1000gen converges',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G20 pop=8 converges',
     &                 0, 1, npass, nfail)
      end if

      return
      end

c     ============================================================
c     Robustness tests
c     ============================================================
      subroutine test_robustness(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

      parameter(ne = 3, na = 3, npts = 5, nob = 5)
      double precision surf(npts, ne, na)
      double precision eg(ne), ag(na)
      double precision tob(nob), apob(nob)
      double precision e_opt, a_opt, rms_opt

      eg(1) = 0.30d0
      eg(2) = 0.32d0
      eg(3) = 0.34d0
      ag(1) = 80.d0
      ag(2) = 120.d0
      ag(3) = 160.d0

c     Flat surface
      do k = 1, na
         do j = 1, ne
            do i = 1, npts
               surf(i, j, k) = 10000.d0
            end do
         end do
      end do
      do i = 1, nob
         tob(i) = dble(i-1) / dble(npts-1)
         apob(i) = 10000.d0
      end do

c     G21: Flat surface â†’ RMS â‰ˆ 0 (any e,A matches)
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 100, 20, 20, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      if (rms_opt .lt. 0.01d0) then
         call assert_i('G21 flat surface rms~0',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G21 flat surface rms~0',
     &                 0, 1, npass, nfail)
      end if

c     G22: Output e within bounds
      if (e_opt .ge. 0.30d0 .and. e_opt .le. 0.34d0) then
         call assert_i('G22 e_opt in bounds',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G22 e_opt in bounds',
     &                 0, 1, npass, nfail)
      end if

c     G23: Output A within bounds
      if (a_opt .ge. 80.d0 .and. a_opt .le. 160.d0) then
         call assert_i('G23 A_opt in bounds',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G23 A_opt in bounds',
     &                 0, 1, npass, nfail)
      end if

c     G24: Fewer generations (50) still produces valid output
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 50, 20, 20, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      if (e_opt .ge. 0.30d0 .and. e_opt .le. 0.34d0 .and.
     &    a_opt .ge. 80.d0 .and. a_opt .le. 160.d0) then
         call assert_i('G24 50 gen valid output',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G24 50 gen valid output',
     &                 0, 1, npass, nfail)
      end if

c     G25: rms_opt is non-negative
      if (rms_opt .ge. 0.d0) then
         call assert_i('G25 rms_opt >= 0',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G25 rms_opt >= 0',
     &                 0, 1, npass, nfail)
      end if

      return
      end

c     ============================================================
c     ga_tinterp piecewise tests
c     ============================================================
      subroutine test_tinterp(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

      double precision s5(5), s1(1), fout

c     Surface: s = [100, 200, 300, 400, 500] at t = [0, 0.25, 0.5, 0.75, 1.0]
      s5(1) = 100.d0
      s5(2) = 200.d0
      s5(3) = 300.d0
      s5(4) = 400.d0
      s5(5) = 500.d0

c     G26: t=0 â†’ first value
      call ga_tinterp(s5, 5, 0.d0, fout)
      call assert_d('G26 tinterp t=0',
     &              fout, 100.d0, 1.d-10, npass, nfail)

c     G27: t=1 â†’ last value
      call ga_tinterp(s5, 5, 1.d0, fout)
      call assert_d('G27 tinterp t=1',
     &              fout, 500.d0, 1.d-10, npass, nfail)

c     G28: t=0.5 â†’ midpoint (300)
      call ga_tinterp(s5, 5, 0.5d0, fout)
      call assert_d('G28 tinterp t=0.5',
     &              fout, 300.d0, 1.d-10, npass, nfail)

c     G29: t=0.125 â†’ interpolated (150)
      call ga_tinterp(s5, 5, 0.125d0, fout)
      call assert_d('G29 tinterp t=0.125',
     &              fout, 150.d0, 1.d-10, npass, nfail)

c     G30: t>1 â†’ clamps to last value
      call ga_tinterp(s5, 5, 1.5d0, fout)
      call assert_d('G30 tinterp t>1 clamp',
     &              fout, 500.d0, 1.d-10, npass, nfail)

c     G31: Single point surface â†’ always returns that value
      s1(1) = 42.d0
      call ga_tinterp(s1, 1, 0.d0, fout)
      call assert_d('G31 tinterp 1pt',
     &              fout, 42.d0, 1.d-10, npass, nfail)

      return
      end

c     ============================================================
c     Extended TWOINT edge cases
c     ============================================================
      subroutine test_twoint_edge(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

      double precision x3(3), y3(3), f3(3,3), fout

      x3(1) = 0.30d0
      x3(2) = 0.32d0
      x3(3) = 0.34d0
      y3(1) = 80.d0
      y3(2) = 120.d0
      y3(3) = 160.d0

c     G32: Quadratic surface f = xÂ² + yÂ²/1000
      do i = 1, 3
         do j = 1, 3
            f3(i,j) = x3(i)**2 + y3(j)**2/1000.d0
         end do
      end do

c     At grid point (0.32, 120): exact
      call ga_twoint(x3, y3, f3, 0.32d0, 120.d0, fout, 3, 3)
      call assert_d('G32 twoint quad gridpt',
     &  fout, 0.32d0**2 + 120.d0**2/1000.d0,
     &  1.d-10, npass, nfail)

c     G33: Interpolation at x_min boundary
      call ga_twoint(x3, y3, f3, 0.30d0, 100.d0, fout, 3, 3)
      if (fout .ne. fout) then
         call assert_i('G33 twoint x_min finite',
     &                 0, 1, npass, nfail)
      else
         call assert_i('G33 twoint x_min finite',
     &                 1, 1, npass, nfail)
      end if

c     G34: Midpoint interpolation symmetric
c     For linear surface f = x + y/100
      do i = 1, 3
         do j = 1, 3
            f3(i,j) = x3(i) + y3(j)/100.d0
         end do
      end do
      call ga_twoint(x3, y3, f3, 0.31d0, 100.d0, fout, 3, 3)
      call assert_d('G34 twoint midpt',
     &  fout, 0.31d0 + 100.d0/100.d0,
     &  1.d-10, npass, nfail)

      return
      end

c     ============================================================
c     Extended decode edge cases
c     ============================================================
      subroutine test_decode_edge(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

      integer chrom(80)
      double precision xvec(2), alow(2), ahigh(2), factor(2)
      integer lsubstr(2)

c     G35: Asymmetric bit allocation (60 + 20)
      lsubstr(1) = 60
      lsubstr(2) = 20
      alow(1) = 0.30d0
      ahigh(1) = 0.34d0
      alow(2) = 80.d0
      ahigh(2) = 160.d0
      factor(1) = 2.d0**60 - 1.d0
      factor(2) = 2.d0**20 - 1.d0

c     All ones â†’ upper bounds
      do i = 1, 80
         chrom(i) = 1
      end do
      call ga_decode(chrom, 80, xvec, 2,
     &               lsubstr, alow, ahigh, factor)
      call assert_d('G35 asymm bits x=hi',
     &              xvec(1), 0.34d0, 1.d-10, npass, nfail)
      call assert_d('G36 asymm bits y=hi',
     &              xvec(2), 160.d0, 1.d-6, npass, nfail)

c     G37: Non-zero lower bound mapping
      lsubstr(1) = 40
      lsubstr(2) = 40
      alow(1) = 0.30d0
      ahigh(1) = 0.34d0
      alow(2) = 100.d0
      ahigh(2) = 200.d0
      factor(1) = 2.d0**40 - 1.d0
      factor(2) = 2.d0**40 - 1.d0

      do i = 1, 80
         chrom(i) = 0
      end do
      call ga_decode(chrom, 80, xvec, 2,
     &               lsubstr, alow, ahigh, factor)
      call assert_d('G37 nonzero lo: x=0.30',
     &              xvec(1), 0.30d0, 1.d-10, npass, nfail)
      call assert_d('G38 nonzero lo: y=100',
     &              xvec(2), 100.d0, 1.d-10, npass, nfail)

      return
      end

c     ============================================================
c     Convergence property tests
c     ============================================================
      subroutine test_convergence(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

      parameter(ne = 3, na = 3, npts = 10, nob = 10)
      double precision surf(npts, ne, na)
      double precision eg(ne), ag(na)
      double precision tob(nob), apob(nob)
      double precision e_opt, a_opt, rms_opt
      double precision e_opt2, a_opt2, rms_opt2

      eg(1) = 0.30d0
      eg(2) = 0.32d0
      eg(3) = 0.34d0
      ag(1) = 80.d0
      ag(2) = 120.d0
      ag(3) = 160.d0

c     Standard synthetic surface
      do k = 1, na
         do j = 1, ne
            do i = 1, npts
               t = dble(i-1) / dble(npts-1)
               surf(i, j, k) = 12000.d0 - 100.d0*t
     &              - 50.d0*(eg(j) - 0.32d0)
     &              - 0.5d0*(ag(k) - 120.d0)
            end do
         end do
      end do
      do i = 1, nob
         tob(i) = dble(i-1) / dble(npts-1)
         apob(i) = 12000.d0 - 100.d0*tob(i)
      end do

c     G39: Same seed â†’ identical results (deterministic)
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 200, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 200, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt2, a_opt2, rms_opt2)
      iall_ok = 1
      if (dabs(e_opt - e_opt2) .gt. 1.d-12) iall_ok = 0
      if (dabs(a_opt - a_opt2) .gt. 1.d-12) iall_ok = 0
      if (dabs(rms_opt - rms_opt2) .gt. 1.d-12) iall_ok = 0
      call assert_i('G39 deterministic (same seed)',
     &              iall_ok, 1, npass, nfail)

c     G40: More generations â†’ equal or better RMS
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 50, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt2, a_opt2, rms_opt2)
      if (rms_opt2 .le. rms_opt + 1.d-6) then
         call assert_i('G40 more gen <= rms',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G40 more gen <= rms',
     &                 0, 1, npass, nfail)
      end if

c     G41: Fitness is always positive
      fit = 1.d0 / (1.d0 + rms_opt2)
      if (fit .gt. 0.d0) then
         call assert_i('G41 fitness > 0',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G41 fitness > 0',
     &                 0, 1, npass, nfail)
      end if

c     G42: RMS is finite (no NaN)
      if (rms_opt2 .eq. rms_opt2) then
         call assert_i('G42 rms finite (no NaN)',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G42 rms finite (no NaN)',
     &                 0, 1, npass, nfail)
      end if

c     G43: e_opt and a_opt are finite
      iall_ok = 1
      if (e_opt2 .ne. e_opt2) iall_ok = 0
      if (a_opt2 .ne. a_opt2) iall_ok = 0
      call assert_i('G43 outputs finite',
     &              iall_ok, 1, npass, nfail)

c     G44: Very few generations (5) â†’ still valid bounds
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 5, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      iall_ok = 1
      if (e_opt .lt. 0.30d0 .or. e_opt .gt. 0.34d0) iall_ok=0
      if (a_opt .lt. 80.d0 .or. a_opt .gt. 160.d0) iall_ok=0
      call assert_i('G44 5gen in bounds',
     &              iall_ok, 1, npass, nfail)

c     G45: 1 generation â†’ still produces output
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 1, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      if (rms_opt .ge. 0.d0) then
         call assert_i('G45 1gen valid output',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G45 1gen valid output',
     &                 0, 1, npass, nfail)
      end if

      return
      end

c     ============================================================
c     Surface edge case tests
c     ============================================================
      subroutine test_surface_edge(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

      parameter(ne = 3, na = 3, npts = 10, nob = 10)
      double precision surf(npts, ne, na)
      double precision eg(ne), ag(na)
      double precision tob(nob), apob(nob)
      double precision e_opt, a_opt, rms_opt

      eg(1) = 0.30d0
      eg(2) = 0.32d0
      eg(3) = 0.34d0
      ag(1) = 80.d0
      ag(2) = 120.d0
      ag(3) = 160.d0

c     G46: Optimum at grid corner (e=0.30, A=80)
      do k = 1, na
         do j = 1, ne
            do i = 1, npts
               t = dble(i-1) / dble(npts-1)
               surf(i, j, k) = 12000.d0 - 100.d0*t
     &              - 50.d0*(eg(j) - 0.30d0)
     &              - 0.5d0*(ag(k) - 80.d0)
            end do
         end do
      end do
      do i = 1, nob
         tob(i) = dble(i-1) / dble(npts-1)
         apob(i) = 12000.d0 - 100.d0*tob(i)
      end do

      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      if (dabs(e_opt - 0.30d0) .lt. 0.01d0) then
         call assert_i('G46 corner opt e~0.30',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G46 corner opt e~0.30',
     &                 0, 1, npass, nfail)
      end if

c     G47: Optimum at grid corner (e=0.34, A=160)
      do k = 1, na
         do j = 1, ne
            do i = 1, npts
               t = dble(i-1) / dble(npts-1)
               surf(i, j, k) = 12000.d0 - 100.d0*t
     &              - 50.d0*(eg(j) - 0.34d0)
     &              - 0.5d0*(ag(k) - 160.d0)
            end do
         end do
      end do
      do i = 1, nob
         tob(i) = dble(i-1) / dble(npts-1)
         apob(i) = 12000.d0 - 100.d0*tob(i)
      end do

      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      if (rms_opt .lt. 1.d0 .and.
     &    e_opt .ge. 0.30d0 .and. e_opt .le. 0.34d0) then
         call assert_i('G47 upper corner converges',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G47 upper corner converges',
     &                 0, 1, npass, nfail)
      end if

c     G48: Steep gradient in e, flat in A
c     Surface: ha = 12000 - 1000*(e-0.32) (no A dependence)
      do k = 1, na
         do j = 1, ne
            do i = 1, npts
               t = dble(i-1) / dble(npts-1)
               surf(i, j, k) = 12000.d0 - 100.d0*t
     &              - 1000.d0*(eg(j) - 0.32d0)
            end do
         end do
      end do
      do i = 1, nob
         tob(i) = dble(i-1) / dble(npts-1)
         apob(i) = 12000.d0 - 100.d0*tob(i)
      end do

      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      if (dabs(e_opt - 0.32d0) .lt. 0.005d0) then
         call assert_i('G48 steep e-gradient found',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G48 steep e-gradient found',
     &                 0, 1, npass, nfail)
      end if

c     G49: Narrow bounds (e range = 0.001)
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.3195d0, 0.3205d0, 80.d0, 160.d0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      if (e_opt .ge. 0.3195d0 .and.
     &    e_opt .le. 0.3205d0) then
         call assert_i('G49 narrow bounds valid',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G49 narrow bounds valid',
     &                 0, 1, npass, nfail)
      end if

c     G50: Fewer obs (3 points) still works
      do i = 1, 3
         tob(i) = dble(i-1) * 0.5d0
         apob(i) = 12000.d0 - 100.d0*tob(i)
      end do
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, 3,
     &   0.30d0, 0.34d0, 80.d0, 160.d0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      if (rms_opt .ge. 0.d0 .and. rms_opt .eq. rms_opt) then
         call assert_i('G50 3obs valid output',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G50 3obs valid output',
     &                 0, 1, npass, nfail)
      end if

      return
      end

c     ============================================================
c     High eccentricity orbit tests
c     ============================================================
      subroutine test_high_ecc(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

      parameter(ne = 3, na = 3, npts = 14, nob = 14)
      double precision surf(npts, ne, na)
      double precision eg(ne), ag(na)
      double precision tob(nob), apob(nob)
      double precision e_opt, a_opt, rms_opt

c     --- 39615-like: e ~ 0.68, ha ~ 35000 km ---
      eg(1) = 0.6806d0
      eg(2) = 0.6808d0
      eg(3) = 0.6810d0
      ag(1) = 1.0d0
      ag(2) = 3.0d0
      ag(3) = 5.0d0

c     SMA ~ 20000 km at this ecc, so ha = a*(1+e)-Re ~ 27000
c     Sensitivity: dha/de = a ~ 20000 km per unit e
      do k = 1, na
         do j = 1, ne
            do i = 1, npts
               t = dble(i-1) / dble(npts-1)
               surf(i,j,k) = 27000.d0 - 200.d0*t
     &           - 20000.d0*(eg(j) - 0.6808d0)
     &           - 10.d0*(ag(k) - 3.d0)
            end do
         end do
      end do
      do i = 1, nob
         tob(i) = dble(i-1) / dble(npts-1)
         apob(i) = 27000.d0 - 200.d0*tob(i)
      end do

c     G51: High-e (0.68) finds optimum
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.6806d0, 0.6810d0, 1.d0, 5.d0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      if (dabs(e_opt - 0.6808d0) .lt. 0.0002d0) then
         call assert_i('G51 e~0.68 converges',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G51 e~0.68 converges',
     &                 0, 1, npass, nfail)
      end if

c     G52: RMS near zero at high-e optimum
      if (rms_opt .lt. 0.5d0) then
         call assert_i('G52 e~0.68 rms < 0.5',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G52 e~0.68 rms < 0.5',
     &                 0, 1, npass, nfail)
      end if

c     --- 35497-like: e ~ 0.63 ---
      eg(1) = 0.6290d0
      eg(2) = 0.6295d0
      eg(3) = 0.6300d0
      ag(1) = 1.0d0
      ag(2) = 2.0d0
      ag(3) = 3.0d0

      do k = 1, na
         do j = 1, ne
            do i = 1, npts
               t = dble(i-1) / dble(npts-1)
               surf(i,j,k) = 35000.d0 - 300.d0*t
     &           - 25000.d0*(eg(j) - 0.6295d0)
     &           - 5.d0*(ag(k) - 2.d0)
            end do
         end do
      end do
      do i = 1, nob
         tob(i) = dble(i-1) / dble(npts-1)
         apob(i) = 35000.d0 - 300.d0*tob(i)
      end do

c     G53: Ariane 5-like e~0.63
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.6290d0, 0.6300d0, 1.d0, 3.d0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      if (dabs(e_opt - 0.6295d0) .lt. 0.0005d0) then
         call assert_i('G53 e~0.63 converges',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G53 e~0.63 converges',
     &                 0, 1, npass, nfail)
      end if

c     --- 37151-like: e ~ 0.56 ---
      eg(1) = 0.5645d0
      eg(2) = 0.5650d0
      eg(3) = 0.5655d0
      ag(1) = 1.0d0
      ag(2) = 3.0d0
      ag(3) = 5.0d0

      do k = 1, na
         do j = 1, ne
            do i = 1, npts
               t = dble(i-1) / dble(npts-1)
               surf(i,j,k) = 15000.d0 - 150.d0*t
     &           - 15000.d0*(eg(j) - 0.5650d0)
     &           - 8.d0*(ag(k) - 3.d0)
            end do
         end do
      end do
      do i = 1, nob
         tob(i) = dble(i-1) / dble(npts-1)
         apob(i) = 15000.d0 - 150.d0*tob(i)
      end do

c     G54: Long March 3B-like e~0.56
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.5645d0, 0.5655d0, 1.d0, 5.d0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      if (dabs(e_opt - 0.5650d0) .lt. 0.0005d0) then
         call assert_i('G54 e~0.56 converges',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G54 e~0.56 converges',
     &                 0, 1, npass, nfail)
      end if

c     --- Extreme sensitivity: e 1000Ã— more sensitive than A ---
      eg(1) = 0.680d0
      eg(2) = 0.681d0
      eg(3) = 0.682d0
      ag(1) = 1.d0
      ag(2) = 3.d0
      ag(3) = 5.d0

      do k = 1, na
         do j = 1, ne
            do i = 1, npts
               t = dble(i-1) / dble(npts-1)
               surf(i,j,k) = 30000.d0 - 100.d0*t
     &           - 50000.d0*(eg(j) - 0.681d0)
     &           - 0.01d0*(ag(k) - 3.d0)
            end do
         end do
      end do
      do i = 1, nob
         tob(i) = dble(i-1) / dble(npts-1)
         apob(i) = 30000.d0 - 100.d0*tob(i)
      end do

c     G55: Extreme asymmetric sensitivity
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.680d0, 0.682d0, 1.d0, 5.d0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      if (dabs(e_opt - 0.681d0) .lt. 0.001d0) then
         call assert_i('G55 extreme asym sensitivity',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G55 extreme asym sensitivity',
     &                 0, 1, npass, nfail)
      end if

c     --- Large apogee altitude scale (GEO-like ~35000 km) ---
      eg(1) = 0.70d0
      eg(2) = 0.72d0
      eg(3) = 0.74d0
      ag(1) = 0.5d0
      ag(2) = 1.5d0
      ag(3) = 2.5d0

      do k = 1, na
         do j = 1, ne
            do i = 1, npts
               t = dble(i-1) / dble(npts-1)
               surf(i,j,k) = 36000.d0 - 500.d0*t
     &           - 30000.d0*(eg(j) - 0.72d0)
     &           - 20.d0*(ag(k) - 1.5d0)
            end do
         end do
      end do
      do i = 1, nob
         tob(i) = dble(i-1) / dble(npts-1)
         apob(i) = 36000.d0 - 500.d0*tob(i)
      end do

c     G56: GEO-scale apogee
      call ga_optimize(surf, npts, eg, ag, ne, na, npts,
     &   tob, apob, nob,
     &   0.70d0, 0.74d0, 0.5d0, 2.5d0,
     &   4, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, a_opt, rms_opt)
      if (dabs(e_opt - 0.72d0) .lt. 0.01d0 .and.
     &    rms_opt .lt. 2.d0) then
         call assert_i('G56 GEO-scale converges',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G56 GEO-scale converges',
     &                 0, 1, npass, nfail)
      end if

c     G57: All outputs finite at high-e
      iall_ok = 1
      if (e_opt .ne. e_opt) iall_ok = 0
      if (a_opt .ne. a_opt) iall_ok = 0
      if (rms_opt .ne. rms_opt) iall_ok = 0
      call assert_i('G57 high-e outputs finite',
     &              iall_ok, 1, npass, nfail)

c     G58: Outputs within bounds at high-e
      if (e_opt .ge. 0.70d0 .and. e_opt .le. 0.74d0 .and.
     &    a_opt .ge. 0.5d0 .and. a_opt .le. 2.5d0) then
         call assert_i('G58 high-e in bounds',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G58 high-e in bounds',
     &                 0, 1, npass, nfail)
      end if

      return
      end

c     ============================================================
c     GA internal piecewise tests
c     ============================================================
      subroutine test_ga_internals(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

      double precision fit(4), sumf, avg, amax, amin
      integer ibest
      double precision oldrand(55)
      integer jrand
      integer chrom(80), child1(20,80), child2(20,80)
      integer parent1(20,80), mate1, mate2
      integer ncross, nmute, iall_ok
      integer ga_mutate, ga_iflip, ga_irnd

c     G59: ga_stats with known fitness [0.1, 0.9, 0.5, 0.3]
      fit(1) = 0.1d0
      fit(2) = 0.9d0
      fit(3) = 0.5d0
      fit(4) = 0.3d0
      call ga_stats(fit, 4, sumf, avg, amax, amin, ibest)
      call assert_d('G59 stats sum=1.8',
     &              sumf, 1.8d0, 1.d-10, npass, nfail)
      call assert_d('G60 stats avg=0.45',
     &              avg, 0.45d0, 1.d-10, npass, nfail)
      call assert_d('G61 stats max=0.9',
     &              amax, 0.9d0, 1.d-10, npass, nfail)
      call assert_d('G62 stats min=0.1',
     &              amin, 0.1d0, 1.d-10, npass, nfail)
      call assert_i('G63 stats ibest=2',
     &              ibest, 2, npass, nfail)

c     G64: ga_stats with equal fitness
      fit(1) = 0.5d0
      fit(2) = 0.5d0
      fit(3) = 0.5d0
      fit(4) = 0.5d0
      call ga_stats(fit, 4, sumf, avg, amax, amin, ibest)
      call assert_d('G64 equal fitness avg=0.5',
     &              avg, 0.5d0, 1.d-10, npass, nfail)

c     G65: ga_mutate preserves value when no flip
      call ga_warmup(0.555d0, oldrand, jrand)
      nmute = 0
      ival = ga_mutate(0, 0.d0, nmute, oldrand, jrand)
c     With prob=0, should return original (0) most of the time
c     Just verify it returns 0 or 1 (valid bit)
      if (ival .eq. 0 .or. ival .eq. 1) then
         call assert_i('G65 mutate returns valid bit',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G65 mutate returns valid bit',
     &                 0, 1, npass, nfail)
      end if

c     G66: ga_mutate flips 0â†’1 when mutation occurs
c     Use pmute=1 to guarantee flip
      nmute = 0
      ival = ga_mutate(0, 1.d0, nmute, oldrand, jrand)
      call assert_i('G66 mutate 0â†’1 w/ pm=1',
     &              ival, 1, npass, nfail)

c     G67: ga_mutate flips 1â†’0 when mutation occurs
      ival = ga_mutate(1, 1.d0, nmute, oldrand, jrand)
      call assert_i('G67 mutate 1â†’0 w/ pm=1',
     &              ival, 0, npass, nfail)

c     G68: ga_iflip with prob=0 â†’ always 0
      call ga_warmup(0.123d0, oldrand, jrand)
      iall_ok = 1
      do i = 1, 100
         if (ga_iflip(0.d0, oldrand, jrand) .ne. 0)
     &      iall_ok = 0
      end do
      call assert_i('G68 iflip prob=0 â†’ 0',
     &              iall_ok, 1, npass, nfail)

c     G69: ga_irnd returns value in [lo, hi]
      call ga_warmup(0.123d0, oldrand, jrand)
      iall_ok = 1
      do i = 1, 500
         k = ga_irnd(1, 80, oldrand, jrand)
         if (k .lt. 1 .or. k .gt. 80) iall_ok = 0
      end do
      call assert_i('G69 irnd in [1,80]',
     &              iall_ok, 1, npass, nfail)

c     G70: ga_irnd produces varied output (not constant)
      k1 = ga_irnd(1, 80, oldrand, jrand)
      ivaried = 0
      do i = 1, 20
         k = ga_irnd(1, 80, oldrand, jrand)
         if (k .ne. k1) ivaried = 1
      end do
      call assert_i('G70 irnd not constant',
     &              ivaried, 1, npass, nfail)

c     G71: ga_fitness with perfect match â†’ fitness â‰ˆ 1
      call test_fitness_perfect(npass, nfail)

      return
      end

c     ============================================================
      subroutine test_fitness_perfect(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail

      parameter(ne = 3, na = 3, npts = 5, nob = 5)
      double precision surf(npts, ne, na)
      double precision eg(ne), ag(na)
      double precision tob(nob), apob(nob)
      double precision xvec(2)

      eg(1) = 0.30d0
      eg(2) = 0.32d0
      eg(3) = 0.34d0
      ag(1) = 80.d0
      ag(2) = 120.d0
      ag(3) = 160.d0

c     All surfaces = constant 10000
      do k = 1, na
         do j = 1, ne
            do i = 1, npts
               surf(i, j, k) = 10000.d0
            end do
         end do
      end do
      do i = 1, nob
         tob(i) = dble(i-1) / dble(npts-1)
         apob(i) = 10000.d0
      end do

c     Any (e, A) should give perfect match
      xvec(1) = 0.32d0
      xvec(2) = 120.d0
      fit = ga_fitness(xvec, surf, npts, eg, ag,
     &                 ne, na, npts, tob, apob, nob)

c     Fitness should be very close to 1.0
      if (fit .gt. 0.99d0) then
         call assert_i('G71 perfect match fit~1',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G71 perfect match fit~1',
     &                 0, 1, npass, nfail)
      end if

      return
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

c     ============================================================
      subroutine assert_d(label, got, expect, tol, npass, nfail)
      implicit double precision (a-h, o-z)
      character*(*) label
      double precision got, expect, tol
      integer npass, nfail
      if (dabs(got - expect) .le. tol) then
         write(*,'(A,A)') '  [PASS]  ', label
         npass = npass + 1
      else
         write(*,'(A,A)') '  [FAIL]  ', label
         write(*,'(A,E20.12,A,E20.12)') '           got=',
     &        got, '  expect=', expect
         nfail = nfail + 1
      end if
      return
      end

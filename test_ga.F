c     test_ga.F — Tests for GA optimizer
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
c     GA optimization — synthetic known-optimum
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

c     G1: Constant surface → interpolation returns constant
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

c     G6: All zeros → both variables = lower bound
      do i = 1, 80
         chrom(i) = 0
      end do
      call ga_decode(chrom, 80, xvec, 2,
     &               lsubstr, alow, ahigh, factor)
      call assert_d('G6 all-0 chr → x=0',
     &              xvec(1), 0.d0, 1.d-10, npass, nfail)
      call assert_d('G7 all-0 chr → y=0',
     &              xvec(2), 0.d0, 1.d-10, npass, nfail)

c     G8: All ones → both variables = upper bound
      do i = 1, 80
         chrom(i) = 1
      end do
      call ga_decode(chrom, 80, xvec, 2,
     &               lsubstr, alow, ahigh, factor)
      call assert_d('G8 all-1 chr → x=1',
     &              xvec(1), 1.d0, 1.d-10, npass, nfail)
      call assert_d('G9 all-1 chr → y=1',
     &              xvec(2), 1.d0, 1.d-10, npass, nfail)

c     G10: First bit only → x = 1/(2^40 - 1) ≈ 0
      do i = 1, 80
         chrom(i) = 0
      end do
      chrom(1) = 1
      call ga_decode(chrom, 80, xvec, 2,
     &               lsubstr, alow, ahigh, factor)
      if (xvec(1) .gt. 0.d0 .and. xvec(1) .lt. 1.d-6) then
         call assert_i('G10 single bit → small x',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G10 single bit → small x',
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
         call assert_i('G12 diff seeds → diff seq',
     &                 1, 1, npass, nfail)
      else
         call assert_i('G12 diff seeds → diff seq',
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

c     Create a synthetic 3×3 surface grid where the optimal
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
      call ga_optimize(surf, eg, ag, ne, na, npts,
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
      call ga_optimize(surf, eg, ag, ne, na, npts,
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
      call ga_optimize(surf, eg, ag, ne, na, npts,
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
      call ga_optimize(surf, eg, ag, ne, na, npts,
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

c     G21: Flat surface → RMS ≈ 0 (any e,A matches)
      call ga_optimize(surf, eg, ag, ne, na, npts,
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
      call ga_optimize(surf, eg, ag, ne, na, npts,
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

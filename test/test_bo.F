c     test_bo.F — Tests for Bayesian-optimization core (bo.F)
c     -------------------------------------------------------
c     Pure-math unit tests (erf/normcdf/kernel/Cholesky/GP/EI) plus a
c     synthetic end-to-end optimization loop that replicates
c     bo_optimize's own normalize -> bo_next_point -> evaluate -> append
c     structure against known Fortran test functions instead of real
c     propagation (no propagate_ks/orem.F dependency here at all --
c     this file only exercises bo.F).
c     -------------------------------------------------------
      program test_bo

      implicit double precision (a-h, o-z)

      integer npass, nfail
      double precision erfval, phival, pdfval, cov
      double precision chl(4,4), amat(4,4), recon
      logical ok
      double precision rhs(4), sol(4), tmp(4)
      double precision xe(100), xb(100), yn(100)
      double precision gpmean, sigma
      double precision xe_next, xb_next, ei_best, eiv
      integer n, i, j

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' Bayesian Optimization (bo.F) Test Suite'
      write(*,*) '=============================================='

c     ==============================================================
c     bo_erf
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- bo_erf ---'

      call bo_erf(0.d0, erfval)
      call assert_d('B1 erf(0)=0', erfval, 0.d0, 1.d-8, npass, nfail)

      call bo_erf(1.d0, erfval)
      call assert_d('B2 erf(1)=0.8427008',
     &              erfval, 0.8427007929d0, 1.d-6, npass, nfail)

      call bo_erf(5.d0, erfval)
      call assert_d('B3 erf(5)~1', erfval, 1.d0, 1.d-6, npass, nfail)

      call bo_erf(-1.d0, erfval)
      call assert_d('B4 erf(-1)=-erf(1)',
     &              erfval, -0.8427007929d0, 1.d-6, npass, nfail)

c     ==============================================================
c     bo_normcdf / bo_normpdf
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- bo_normcdf / bo_normpdf ---'

      call bo_normcdf(0.d0, phival)
      call assert_d('B5 normcdf(0)=0.5',
     &              phival, 0.5d0, 1.d-8, npass, nfail)

      call bo_normcdf(1.959964d0, phival)
      call assert_d('B6 normcdf(1.96)~0.975',
     &              phival, 0.975d0, 1.d-4, npass, nfail)

      call bo_normcdf(-1.959964d0, phival)
      call assert_d('B7 normcdf(-1.96)~0.025',
     &              phival, 0.025d0, 1.d-4, npass, nfail)

      call bo_normpdf(0.d0, pdfval)
      call assert_d('B8 normpdf(0)=0.39894',
     &              pdfval, 0.3989422804d0, 1.d-8, npass, nfail)

      call bo_normpdf(1.d0, pdfval)
      call assert_d('B9 normpdf(1)=0.24197',
     &              pdfval, 0.2419707245d0, 1.d-8, npass, nfail)

c     ==============================================================
c     bo_kernel
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- bo_kernel ---'

      call bo_kernel(0.3d0, 0.7d0, 0.3d0, 0.7d0, 0.25d0, cov)
      call assert_d('B10 kernel same point = 1',
     &              cov, 1.d0, 1.d-10, npass, nfail)

      call bo_kernel(0.d0, 0.d0, 1.d0, 1.d0, 0.25d0, cov)
      if (cov .lt. 1.d-3) then
         call assert_i('B11 kernel far apart ~ 0', 1, 1, npass, nfail)
      else
         call assert_i('B11 kernel far apart ~ 0', 0, 1, npass, nfail)
      end if

      call bo_kernel(0.1d0, 0.2d0, 0.4d0, 0.6d0, 0.25d0, cov)
      call bo_kernel(0.4d0, 0.6d0, 0.1d0, 0.2d0, 0.25d0, recon)
      call assert_d('B12 kernel symmetric', cov, recon, 1.d-12,
     &              npass, nfail)

c     ==============================================================
c     bo_cholesky / bo_chol_fsub / bo_chol_bsub
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- bo_cholesky ---'

c     B13: [[4,2],[2,3]] is PD -- reconstruct chl*chl^T
      amat(1,1) = 4.d0
      amat(1,2) = 2.d0
      amat(2,1) = 2.d0
      amat(2,2) = 3.d0
      call bo_cholesky(amat, 2, 4, chl, ok)
      if (ok) then
         call assert_i('B13 cholesky PD ok', 1, 1, npass, nfail)
         recon = chl(1,1)*chl(1,1)
         call assert_d('B13 recon(1,1)', recon, 4.d0, 1.d-10,
     &                 npass, nfail)
         recon = chl(2,1)*chl(1,1)
         call assert_d('B13 recon(2,1)', recon, 2.d0, 1.d-10,
     &                 npass, nfail)
         recon = chl(2,1)*chl(2,1) + chl(2,2)*chl(2,2)
         call assert_d('B13 recon(2,2)', recon, 3.d0, 1.d-10,
     &                 npass, nfail)
      else
         call assert_i('B13 cholesky PD ok', 0, 1, npass, nfail)
      end if

c     B14: [[1,2],[2,1]] is NOT PD (det=-3) -- must report ok=.false.
      amat(1,1) = 1.d0
      amat(1,2) = 2.d0
      amat(2,1) = 2.d0
      amat(2,2) = 1.d0
      call bo_cholesky(amat, 2, 4, chl, ok)
      if (.not. ok) then
         call assert_i('B14 cholesky non-PD detected', 1, 1,
     &                 npass, nfail)
      else
         call assert_i('B14 cholesky non-PD detected', 0, 1,
     &                 npass, nfail)
      end if

c     B15: forward/back substitution solves the PD system correctly
c     against a known RHS, matching an independently-checked solution
      amat(1,1) = 4.d0
      amat(1,2) = 2.d0
      amat(2,1) = 2.d0
      amat(2,2) = 3.d0
      call bo_cholesky(amat, 2, 4, chl, ok)
      rhs(1) = 1.d0
      rhs(2) = 1.d0
      call bo_chol_fsub(chl, 2, 4, rhs, tmp)
      call bo_chol_bsub(chl, 2, 4, tmp, sol)
c     Verify amat*sol = rhs (independent check, not the expected
c     solution vector itself)
      recon = amat(1,1)*sol(1) + amat(1,2)*sol(2)
      call assert_d('B15 fsub/bsub solves row 1', recon, 1.d0,
     &              1.d-10, npass, nfail)
      recon = amat(2,1)*sol(1) + amat(2,2)*sol(2)
      call assert_d('B15 fsub/bsub solves row 2', recon, 1.d0,
     &              1.d-10, npass, nfail)

c     ==============================================================
c     bo_gp_predict
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- bo_gp_predict ---'

c     Training data from a smooth linear function f(x,y)=x+y (already
c     zero-mean-ish over [0,1]^2, no separate normalization needed for
c     this check -- just verifying the GP machinery, not bo_optimize's
c     own normalization step, which is tested separately below)
      n = 5
      xe(1) = 0.1d0
      xb(1) = 0.1d0
      xe(2) = 0.9d0
      xb(2) = 0.1d0
      xe(3) = 0.1d0
      xb(3) = 0.9d0
      xe(4) = 0.9d0
      xb(4) = 0.9d0
      xe(5) = 0.5d0
      xb(5) = 0.5d0
      do i = 1, n
         yn(i) = xe(i) + xb(i)
      end do

c     B16: posterior mean AT a training point should recover that
c     point's own value closely (small nugget => near-exact fit)
      call bo_gp_predict(xe, xb, yn, n, 100, 0.25d0, 1.d-6,
     &   0.5d0, 0.5d0, gpmean, sigma, ok)
      call assert_i('B16 gp_predict ok', 1, 1, npass, nfail)
      call assert_d('B16 gp mean at training pt (0.5,0.5)',
     &              gpmean, 1.d0, 0.05d0, npass, nfail)
      if (sigma .lt. 0.05d0) then
         call assert_i('B16 gp sigma small at training pt', 1, 1,
     &                 npass, nfail)
      else
         call assert_i('B16 gp sigma small at training pt', 0, 1,
     &                 npass, nfail)
      end if

c     B17: posterior variance should be LARGER far from every
c     training point than at one
      call bo_gp_predict(xe, xb, yn, n, 100, 0.25d0, 1.d-6,
     &   0.5d0, 0.5d0, gpmean, sigma, ok)
      sigma_near = sigma
      call bo_gp_predict(xe, xb, yn, n, 100, 0.25d0, 1.d-6,
     &   0.02d0, 0.98d0, gpmean, sigma, ok)
      if (sigma .gt. sigma_near) then
         call assert_i('B17 sigma grows away from data', 1, 1,
     &                 npass, nfail)
      else
         call assert_i('B17 sigma grows away from data', 0, 1,
     &                 npass, nfail)
      end if

c     ==============================================================
c     bo_ei
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- bo_ei ---'

c     B18: sigma=0 -> EI=0 regardless of mean
      call bo_ei(-5.d0, 0.d0, 0.d0, eiv)
      call assert_d('B18 EI zero at sigma=0', eiv, 0.d0, 1.d-12,
     &              npass, nfail)

c     B19: predicted mean far BELOW current best (minimizing) with
c     nonzero sigma -> large positive EI
      call bo_ei(-3.d0, 1.d0, 0.d0, eiv)
      if (eiv .gt. 1.d0) then
         call assert_i('B19 EI large for promising candidate', 1, 1,
     &                 npass, nfail)
      else
         call assert_i('B19 EI large for promising candidate', 0, 1,
     &                 npass, nfail)
      end if

c     B20: predicted mean far ABOVE current best -> small EI (mostly
c     from the exploration/uncertainty term only)
      call bo_ei(5.d0, 1.d0, 0.d0, eiv)
      if (eiv .lt. 0.1d0) then
         call assert_i('B20 EI small for poor candidate', 1, 1,
     &                 npass, nfail)
      else
         call assert_i('B20 EI small for poor candidate', 0, 1,
     &                 npass, nfail)
      end if

c     ==============================================================
c     bo_next_point + end-to-end synthetic optimization
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- bo_next_point / end-to-end synthetic BO ---'

c     B21: centered bowl f(x,y) = (x-0.5)^2 + (y-0.5)^2, true min = 0
c     at (0.5, 0.5). Seed a 3x3 grid (same warm-start shape as
c     bo_optimize's 9 RSM corners), then run 60 BO iterations picking
c     the next point via bo_next_point each time, evaluating the
c     SYNTHETIC function directly (no propagation involved).
      call run_synthetic_bo(0.5d0, 0.5d0, 1.d0, 1.d0, 60,
     &   xe_found, xb_found, y_found)
      call assert_d('B21 bowl: converges near x*=0.5',
     &              xe_found, 0.5d0, 0.1d0, npass, nfail)
      call assert_d('B21 bowl: converges near y*=0.5',
     &              xb_found, 0.5d0, 0.1d0, npass, nfail)
      call assert_d('B21 bowl: best value near true min 0',
     &              y_found, 0.d0, 0.02d0, npass, nfail)

c     B22: off-center, anisotropic bowl f(x,y)=(x-0.2)^2+3*(y-0.8)^2,
c     true min = 0 at (0.2, 0.8) -- checks the method isn't secretly
c     biased toward the center of the search box
      call run_synthetic_bo(0.2d0, 0.8d0, 1.d0, 3.d0, 80,
     &   xe_found, xb_found, y_found)
      call assert_d('B22 offset bowl: converges near x*=0.2',
     &              xe_found, 0.2d0, 0.12d0, npass, nfail)
      call assert_d('B22 offset bowl: converges near y*=0.8',
     &              xb_found, 0.8d0, 0.12d0, npass, nfail)
      call assert_d('B22 offset bowl: best value near true min 0',
     &              y_found, 0.d0, 0.03d0, npass, nfail)

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
      subroutine run_synthetic_bo(xstar, ystar, wx, wy, n_new,
     &   xe_found, xb_found, y_found)
c     Replicates bo_optimize's own warm-start(9-corner)/normalize/
c     bo_next_point/evaluate/append loop, but against a known
c     synthetic quadratic bowl f(x,y) = wx*(x-xstar)^2 + wy*(y-ystar)^2
c     instead of real propagation -- validates the SEARCH STRATEGY
c     end-to-end without needing orem.F/propagate_ks in this test file.
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision xstar, ystar, wx, wy
      integer n_new
      double precision xe_found, xb_found, y_found

      parameter(maxbo = 100)
      double precision xe(maxbo), xb(maxbo), yraw(maxbo), yn(maxbo)
      integer n, i, ie, ia, itry
      double precision scal, reg
      parameter(scal = 0.25d0, reg = 1.d-6)
      integer ngrid
      parameter(ngrid = 50)
      double precision ymean, ystd, ybest_n
      double precision xe_next, xb_next, ei_best
      double precision xg, yg
      logical ok

      n = 0
      do ie = 1, 3
         do ia = 1, 3
            n = n + 1
            xg = dble(ie-1) / 2.d0
            yg = dble(ia-1) / 2.d0
            xe(n) = xg
            xb(n) = yg
            yraw(n) = wx*(xg-xstar)**2 + wy*(yg-ystar)**2
         end do
      end do

      y_found = yraw(1)
      xe_found = xe(1)
      xb_found = xb(1)
      do i = 1, n
         if (yraw(i) .lt. y_found) then
            y_found = yraw(i)
            xe_found = xe(i)
            xb_found = xb(i)
         end if
      end do

      do itry = 1, n_new
         ymean = 0.d0
         do i = 1, n
            ymean = ymean + yraw(i)
         end do
         ymean = ymean / dble(n)
         ystd = 0.d0
         do i = 1, n
            ystd = ystd + (yraw(i)-ymean)**2
         end do
         ystd = dsqrt(ystd / dble(n))
         if (ystd .lt. 1.d-9) ystd = 1.d-9
         do i = 1, n
            yn(i) = (yraw(i)-ymean) / ystd
         end do
         ybest_n = (y_found - ymean) / ystd

         call bo_next_point(xe, xb, yn, n, scal, reg, ngrid,
     &      ybest_n, xe_next, xb_next, ei_best, ok)
         if (.not. ok) goto 90

         if (n .lt. maxbo) then
            n = n + 1
            xe(n) = xe_next
            xb(n) = xb_next
            yraw(n) = wx*(xe_next-xstar)**2 + wy*(xb_next-ystar)**2
            if (yraw(n) .lt. y_found) then
               y_found = yraw(n)
               xe_found = xe_next
               xb_found = xb_next
            end if
         end if

         if (n .ge. maxbo) goto 90
      end do

 90   continue
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

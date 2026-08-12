c     ============================================================
c     bo.F — Bayesian optimization (GP + Expected Improvement) for
c     zone-level (e, BN) search against DIRECT-PROPAGATION fitness.
c
c     Issue #40 option 3 follow-up: option 3's own offline Python study
c     (scratch_rpe/bench/gasweep/bo_vs_ga.py) found a from-scratch GP-EI
c     search reaches ~72% of the GA's full-budget reliability at 1% of
c     the evaluation count -- but only against the CHEAP interpolated
c     RSM surface, explicitly flagged as untested against real direct-
c     propagation fitness (issue #40's caveat 2). This file is that
c     untested piece: the same GP-EI method, in Fortran, driving
c     rsm_direct_eval (real propagate_ks calls) instead of the
c     interpolated surface, as a full alternative search strategy to
c     ga_optimize (not a post-hoc refinement of it, unlike irefine).
c
c     Kernel: isotropic RBF, FIXED lengthscale (scal=0.25 in normalized
c     [0,1]x[0,1] (e,BN) space) -- same untuned choice as the Python
c     study, for direct comparability; a properly-tuned GP would likely
c     do better still (documented there as a floor, not a ceiling).
c     Observations are normalized (zero mean, unit variance) before GP
c     fitting so the fixed lengthscale/nugget choices are robust across
c     zones whose raw RMS values span very different physical scales.
c
c     Author: Harishkumar Sellamuthu
c     ============================================================

c     ============================================================
      subroutine bo_erf(x, erfout)
c     Abramowitz & Stegun 7.1.26 rational approximation, max error
c     1.5e-7 -- plenty for an acquisition function that only needs to
c     rank candidates, not integrate to high precision.
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision x, erfout
      double precision p, a1, a2, a3, a4, a5
      double precision t, ax, poly
      parameter(p = 0.3275911d0)
      parameter(a1 = 0.254829592d0, a2 = -0.284496736d0)
      parameter(a3 = 1.421413741d0, a4 = -1.453152027d0)
      parameter(a5 = 1.061405429d0)

      ax = dabs(x)
      t = 1.d0 / (1.d0 + p*ax)
      poly = t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))))
      erfout = 1.d0 - poly*dexp(-ax*ax)
      if (x .lt. 0.d0) erfout = -erfout

      return
      end

c     ============================================================
      subroutine bo_normcdf(z, phiout)
      implicit double precision (a-h, o-z)
      double precision z, phiout, erfval
      double precision sqrt2
      parameter(sqrt2 = 1.4142135623730951d0)

      call bo_erf(z/sqrt2, erfval)
      phiout = 0.5d0*(1.d0 + erfval)

      return
      end

c     ============================================================
      subroutine bo_normpdf(z, pdfout)
      implicit double precision (a-h, o-z)
      double precision z, pdfout
      double precision inv_sqrt2pi
      parameter(inv_sqrt2pi = 0.3989422804014327d0)

      pdfout = inv_sqrt2pi * dexp(-0.5d0*z*z)

      return
      end

c     ============================================================
      subroutine bo_kernel(xe1, xb1, xe2, xb2, scal, cov)
c     Isotropic RBF kernel over normalized 2D (e, BN) coordinates.
c     k(x,x) = 1 always (unit signal variance) -- keeps the posterior
c     variance formula below (var = 1 - v.v) simple.
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision xe1, xb1, xe2, xb2, scal, cov
      double precision d2

      d2 = (xe1-xe2)**2 + (xb1-xb2)**2
      cov = dexp(-d2 / (2.d0*scal*scal))

      return
      end

c     ============================================================
      subroutine bo_cholesky(amat, n, maxn, chl, ok)
c     Cholesky decomposition of an n x n symmetric positive-definite
c     matrix amat (stored in a maxn x maxn array): amat = chl * chl^T,
c     chl lower-triangular. ok=.false. if amat is not (numerically)
c     positive-definite -- caller's amat should always include the
c     regularizing nugget on its diagonal to make this exceedingly
c     unlikely in practice.
c     ============================================================
      implicit double precision (a-h, o-z)
      integer n, maxn
      double precision amat(maxn, maxn), chl(maxn, maxn)
      logical ok

      double precision ssum
      integer i, j, k

      ok = .true.

      do j = 1, n
         do i = 1, n
            chl(i, j) = 0.d0
         end do
      end do

      do j = 1, n
         ssum = amat(j, j)
         do k = 1, j-1
            ssum = ssum - chl(j, k)**2
         end do
         if (ssum .le. 0.d0) then
            ok = .false.
            return
         end if
         chl(j, j) = dsqrt(ssum)

         do i = j+1, n
            ssum = amat(i, j)
            do k = 1, j-1
               ssum = ssum - chl(i, k)*chl(j, k)
            end do
            chl(i, j) = ssum / chl(j, j)
         end do
      end do

      return
      end

c     ============================================================
      subroutine bo_chol_fsub(chl, n, maxn, rhs, sol)
c     Forward substitution: solves chl * sol = rhs (chl lower-tri).
c     ============================================================
      implicit double precision (a-h, o-z)
      integer n, maxn
      double precision chl(maxn, maxn), rhs(maxn), sol(maxn)
      double precision ssum
      integer i, k

      do i = 1, n
         ssum = rhs(i)
         do k = 1, i-1
            ssum = ssum - chl(i, k)*sol(k)
         end do
         sol(i) = ssum / chl(i, i)
      end do

      return
      end

c     ============================================================
      subroutine bo_chol_bsub(chl, n, maxn, rhs, sol)
c     Back substitution: solves chl^T * sol = rhs (chl lower-tri).
c     ============================================================
      implicit double precision (a-h, o-z)
      integer n, maxn
      double precision chl(maxn, maxn), rhs(maxn), sol(maxn)
      double precision ssum
      integer i, k

      do i = n, 1, -1
         ssum = rhs(i)
         do k = i+1, n
            ssum = ssum - chl(k, i)*sol(k)
         end do
         sol(i) = ssum / chl(i, i)
      end do

      return
      end

c     ============================================================
      subroutine bo_gp_predict(xe, xb, yn, n, maxn, scal, reg,
     &   xeq, xbq, gpmean, sigma, ok)
c     GP posterior mean/std at a single query point (xeq, xbq), fit
c     fresh from the n observed (xe, xb, yn) triples every call --
c     simple and correct, and fine for standalone use/testing. The
c     production search path (bo_next_point below) does NOT call this:
c     it needs the same Cholesky/alpha across ~2600 grid points per
c     iteration and would be needlessly O(grid^2 * n^3) refitting from
c     scratch at every grid point, so it inlines an equivalent fit-once,
c     query-many version instead.
c
c     Input
c       xe(n), xb(n) : observed points, normalized to [0,1]x[0,1]
c       yn(n)        : observed values, normalized (zero mean, unit
c                      variance) by the caller
c       scal         : RBF lengthscale
c       reg          : diagonal nugget (numerical jitter)
c       xeq, xbq     : query point
c     Output
c       gpmean, sigma : posterior mean/std (normalized units)
c       ok            : .false. if the Cholesky factorization failed
c     ============================================================
      implicit double precision (a-h, o-z)
      integer n, maxn
      double precision xe(maxn), xb(maxn), yn(maxn)
      double precision scal, reg, xeq, xbq
      double precision gpmean, sigma
      logical ok

      parameter(maxbo = 100)
      double precision amat(maxbo, maxbo), chl(maxbo, maxbo)
      double precision ztmp(maxbo), alpha(maxbo)
      double precision kstar(maxbo), vvec(maxbo)
      double precision cov, var
      integer i, j

      if (n .gt. maxbo) then
         ok = .false.
         return
      end if

      do i = 1, n
         do j = 1, n
            call bo_kernel(xe(i), xb(i), xe(j), xb(j), scal, cov)
            amat(i, j) = cov
         end do
         amat(i, i) = amat(i, i) + reg
      end do

      call bo_cholesky(amat, n, maxbo, chl, ok)
      if (.not. ok) return

      call bo_chol_fsub(chl, n, maxbo, yn, ztmp)
      call bo_chol_bsub(chl, n, maxbo, ztmp, alpha)

      do i = 1, n
         call bo_kernel(xe(i), xb(i), xeq, xbq, scal, cov)
         kstar(i) = cov
      end do

      gpmean = 0.d0
      do i = 1, n
         gpmean = gpmean + kstar(i)*alpha(i)
      end do

      call bo_chol_fsub(chl, n, maxbo, kstar, vvec)
      var = 1.d0
      do i = 1, n
         var = var - vvec(i)*vvec(i)
      end do
      if (var .lt. 0.d0) var = 0.d0
      sigma = dsqrt(var)

      return
      end

c     ============================================================
      subroutine bo_ei(gpmean, sigma, ybest, eiv)
c     Expected improvement for MINIMIZATION (lower fitness is better,
c     matching every fitness metric elsewhere in this repo -- RMS, not
c     GA's 1/(1+RMS)). ybest, gpmean, sigma are all in the GP's own
c     normalized units -- fine, since EI's ranking is invariant to any
c     shared affine rescaling of y.
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision gpmean, sigma, ybest, eiv
      double precision eimp, z, phival, pdfval

      if (sigma .lt. 1.d-9) then
         eiv = 0.d0
         return
      end if

      eimp = ybest - gpmean
      z = eimp / sigma
      call bo_normcdf(z, phival)
      call bo_normpdf(z, pdfval)
      eiv = eimp*phival + sigma*pdfval
      if (eiv .lt. 0.d0) eiv = 0.d0

      return
      end

c     ============================================================
      subroutine bo_next_point(xe, xb, yn, n, scal, reg, ngrid,
     &   ybest, xe_next, xb_next, ei_best, ok)
c     Picks the next point to evaluate: the (xe,xb) in [0,1]x[0,1]
c     maximizing Expected Improvement, by dense grid search
c     ((ngrid+1) x (ngrid+1) points) -- cheap (pure arithmetic, no
c     propagation) and matches this project's existing brute-force-
c     over-a-grid precedent (option 3's own 250x250 ground-truth scan).
c     Fits the GP ONCE (Cholesky + alpha), then re-uses that
c     factorization across the whole grid -- see bo_gp_predict's header
c     for why this can't just call bo_gp_predict per grid point.
c     ============================================================
      implicit double precision (a-h, o-z)
      integer n, ngrid
      double precision xe(n), xb(n), yn(n)
      double precision scal, reg, ybest
      double precision xe_next, xb_next, ei_best
      logical ok

      parameter(maxbo = 100)
      double precision amat(maxbo, maxbo), chl(maxbo, maxbo)
      double precision ztmp(maxbo), alpha(maxbo)
      double precision kstar(maxbo), vvec(maxbo)
      double precision cov, var, gpmean, sigma, eiv
      double precision xeq, xbq
      integer i, j, ig, jg

      ei_best = -1.d0
      xe_next = 0.5d0
      xb_next = 0.5d0
      ok = .false.

      if (n .gt. maxbo .or. n .lt. 1) return

      do i = 1, n
         do j = 1, n
            call bo_kernel(xe(i), xb(i), xe(j), xb(j), scal, cov)
            amat(i, j) = cov
         end do
         amat(i, i) = amat(i, i) + reg
      end do

      call bo_cholesky(amat, n, maxbo, chl, ok)
      if (.not. ok) return

      call bo_chol_fsub(chl, n, maxbo, yn, ztmp)
      call bo_chol_bsub(chl, n, maxbo, ztmp, alpha)

      do ig = 0, ngrid
         xeq = dble(ig) / dble(ngrid)
         do jg = 0, ngrid
            xbq = dble(jg) / dble(ngrid)

            do i = 1, n
               call bo_kernel(xe(i), xb(i), xeq, xbq, scal, cov)
               kstar(i) = cov
            end do

            gpmean = 0.d0
            do i = 1, n
               gpmean = gpmean + kstar(i)*alpha(i)
            end do

            call bo_chol_fsub(chl, n, maxbo, kstar, vvec)
            var = 1.d0
            do i = 1, n
               var = var - vvec(i)*vvec(i)
            end do
            if (var .lt. 0.d0) var = 0.d0
            sigma = dsqrt(var)

            call bo_ei(gpmean, sigma, ybest, eiv)
            if (eiv .gt. ei_best) then
               ei_best = eiv
               xe_next = xeq
               xb_next = xbq
            end if
         end do
      end do

      ok = .true.

      return
      end

c     ============================================================
      subroutine bo_optimize(
     &   a_zone1, e_zone1, inc_zone1, raan_zone1, aop_zone1, ma_zone1,
     &   epoch0_jd, zone_dur,
     &   surfaces, ld_surf, e_grid, bn_grid,
     &   tobs, apobs, nobs,
     &   istep_prop, tole_prop,
     &   n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &   WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD,
     &   PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   n_new_eval,
     &   e_opt, bn_opt, rms_opt, neval_total, neval_ok)
c
c     Issue #40: full alternative to ga_optimize -- searches (e, BN)
c     directly against rsm_direct_eval's REAL propagated fitness using
c     Bayesian optimization (GP + Expected Improvement), instead of the
c     GA searching the cheap 9-corner bilinear-interpolated surface.
c
c     Warm start: the 9 RSM grid corners (surfaces/e_grid/bn_grid,
c     already computed by rsm_generate for this zone) ARE exact real-
c     propagation evaluations -- a grid corner's own bilinear weight is
c     1, so reading surfaces(i,ie,ia) directly is a zero-error, zero-
c     extra-cost real fitness value, not an interpolation. Re-using
c     them as the GP's initial design means the n_new_eval budget below
c     is spent entirely on genuinely NEW information.
c
c     Input
c       a/e/inc/raan/aop/ma_zone1, epoch0_jd, zone_dur : as
c                                    rsm_direct_eval (zone's base state)
c       surfaces, ld_surf, e_grid, bn_grid : rsm_generate's own outputs
c                                    for this zone (9-corner grid)
c       tobs, apobs, nobs           : as rsm_generate/rsm_direct_eval
c       istep_prop, tole_prop, force model, atmosphere : as
c                                    rsm_direct_eval
c       n_new_eval                  : direct-propagation budget beyond
c                                    the 9 free corners (e.g. 91, for
c                                    100 total -- issue #40's option 3
c                                    validated figure)
c     Output
c       e_opt, bn_opt, rms_opt : best (e, BN, direct RMS) found across
c                                all evaluated points (warm-start +
c                                BO-selected)
c       neval_total : 9 + number of new-point attempts (successful or
c                     not)
c       neval_ok    : number of those new-point attempts that actually
c                     propagated successfully
c     ============================================================
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision a_zone1, e_zone1, inc_zone1
      double precision raan_zone1, aop_zone1, ma_zone1
      double precision epoch0_jd, zone_dur
      integer ld_surf
      double precision surfaces(ld_surf, 3, 3)
      double precision e_grid(3), bn_grid(3)
      integer nobs
      double precision tobs(nobs), apobs(nobs)
      integer istep_prop
      double precision tole_prop
      integer n_force(3), ngeo_deg, nsun_deg, nmoon_deg
      double precision WE_rot, EPS_f, FR_rot
      double precision CR_srp, AM_srp
      integer IPSR, ISHAD
      double precision PSR_srp, amuS, amuM
      integer ndim_atm
      double precision ALT_atm(ndim_atm)
      double precision DEN_atm(ndim_atm)
      double precision SCH_atm(ndim_atm)
      integer n_new_eval
      double precision e_opt, bn_opt, rms_opt
      integer neval_total, neval_ok

      parameter(maxbo = 100)
      double precision xe(maxbo), xb(maxbo), yraw(maxbo), yn(maxbo)
      integer n

      double precision scal, reg
      parameter(scal = 0.25d0, reg = 1.d-6)
      integer ngrid
      parameter(ngrid = 50)

      double precision ymean, ystd, ybest_n
      double precision xe_next, xb_next, ei_best
      double precision e_cand, bn_cand, rms_cand
      double precision e_lo, e_hi, bn_lo, bn_hi
      double precision ssq
      integer ie, ia, i, itry, ierr_de
      logical ok

      e_lo = e_grid(1)
      e_hi = e_grid(3)
      bn_lo = bn_grid(1)
      bn_hi = bn_grid(3)

      n = 0
      do ie = 1, 3
         do ia = 1, 3
            n = n + 1
            xe(n) = (e_grid(ie) - e_lo) / (e_hi - e_lo)
            xb(n) = (bn_grid(ia) - bn_lo) / (bn_hi - bn_lo)
            ssq = 0.d0
            do i = 1, nobs
               ssq = ssq + (surfaces(i, ie, ia) - apobs(i))**2
            end do
            yraw(n) = dsqrt(ssq / dble(nobs))
         end do
      end do

      e_opt = e_grid(2)
      bn_opt = bn_grid(2)
      rms_opt = yraw(1)
      do i = 1, n
         if (yraw(i) .lt. rms_opt) then
            rms_opt = yraw(i)
            e_opt = e_lo + xe(i)*(e_hi-e_lo)
            bn_opt = bn_lo + xb(i)*(bn_hi-bn_lo)
         end if
      end do

      neval_ok = 0

      do itry = 1, n_new_eval
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
         ybest_n = (rms_opt - ymean) / ystd

         call bo_next_point(xe, xb, yn, n, scal, reg, ngrid,
     &      ybest_n, xe_next, xb_next, ei_best, ok)
         if (.not. ok) goto 90

         e_cand = e_lo + xe_next*(e_hi-e_lo)
         bn_cand = bn_lo + xb_next*(bn_hi-bn_lo)

         call rsm_direct_eval(
     &      a_zone1, e_zone1, inc_zone1, raan_zone1, aop_zone1,
     &      ma_zone1, epoch0_jd, zone_dur,
     &      e_cand, bn_cand,
     &      tobs, apobs, nobs, istep_prop, tole_prop,
     &      n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &      WE_rot, EPS_f, FR_rot,
     &      CR_srp, AM_srp, IPSR, ISHAD,
     &      PSR_srp, amuS, amuM,
     &      ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &      rms_cand, ierr_de)

         if (ierr_de .eq. 0) then
            neval_ok = neval_ok + 1
            if (n .lt. maxbo) then
               n = n + 1
               xe(n) = xe_next
               xb(n) = xb_next
               yraw(n) = rms_cand
            end if
            if (rms_cand .lt. rms_opt) then
               rms_opt = rms_cand
               e_opt = e_cand
               bn_opt = bn_cand
            end if
         end if

         if (n .ge. maxbo) goto 90
      end do

 90   continue
      neval_total = 9 + n_new_eval

      return
      end

c     test_coverage.F -- tests for routines that had no direct test,
c     plus regression checks for the 2026-10-12 code-review findings.
c     -------------------------------------------------------
c     Covers: ga_campaign, ga_camp_fitness, ga_slope_optimize,
c     ga_slope_fitness, ga_preselect, ga_select2, ga_crossover,
c     ga_advance, rsm_tinterp, rsm_direct_eval, rsm_refine,
c     zone_select_next, linfit, compute_rpe, tle_find_ma,
c     tle_find_osc, tle_find_bstar, estimate_bn_bstar_prior,
c     estimate_bn_floor, propagate_ks force switches.
c     Run from the repo root (reads input/ATM.DAT).
c     -------------------------------------------------------
      program test_coverage
      implicit double precision (a-h, o-z)
      integer npass, nfail

      npass = 0
      nfail = 0
      call init_constants()

      write(*,*) '=============================================='
      write(*,*) ' Coverage / code-review regression tests'
      write(*,*) '=============================================='

      write(*,*) ''
      write(*,*) '  --- GA operators ---'
      call t_ga_ops(npass, nfail)

      write(*,*) ''
      write(*,*) '  --- GA population layout (pcross=pmute=0) ---'
      call t_ga_layout(npass, nfail)

      write(*,*) ''
      write(*,*) '  --- GA campaign / slope optimizers ---'
      call t_ga_campaign(npass, nfail)
      call t_ga_slope(npass, nfail)

      write(*,*) ''
      write(*,*) '  --- rsm_tinterp / compute_rpe / linfit ---'
      call t_tinterp(npass, nfail)
      call t_rpe(npass, nfail)
      call t_linfit(npass, nfail)

      write(*,*) ''
      write(*,*) '  --- zone_select_next ---'
      call t_zone_next(npass, nfail)

      write(*,*) ''
      write(*,*) '  --- TLE lookup helpers / BN priors ---'
      call t_tle_helpers(npass, nfail)

      write(*,*) ''
      write(*,*) '  --- propagate_ks force switches ---'
      call t_force_switch(npass, nfail)

      write(*,*) ''
      write(*,*) '  --- rsm_direct_eval / rsm_refine ---'
      call t_direct_refine(npass, nfail)

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
c     GA operators: ga_advance range, ga_crossover bit conservation,
c     ga_preselect termination on a flat population, ga_select2
c     range, NaN guards in the campaign/slope fitness functions.
c     ============================================================
      subroutine t_ga_ops(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail
      parameter(ld = 20, lc = 30)
      double precision oldrand(55), fit(8), frac(8), sl(3,3)
      double precision eg(3), bg(3), xv(2)
      integer oldchr(ld, lc), newchr(ld, lc), ich(8)
      integer jrand, ncross, nmute, ok, i, j, k, nseen
      integer ga_iflip
      double precision ga_slope_fitness

c     C1: ga_advance keeps every RNG slot in [0,1)
      call ga_warmup(0.5d0, oldrand, jrand)
      ok = 1
      do k = 1, 200
         call ga_advance(oldrand)
         do i = 1, 55
            if (oldrand(i) .lt. 0.d0 .or. oldrand(i) .ge. 1.d0) ok = 0
         end do
      end do
      call ck_i('C1 ga_advance keeps values in [0,1)', ok, 1,
     &          npass, nfail)

c     C2: single-point crossover (pmute=0) conserves bits: at every
c     locus the two children carry exactly the two parents' bits
      call ga_warmup(0.321d0, oldrand, jrand)
      do i = 1, 4
         do k = 1, lc
            oldchr(i, k) = ga_iflip(0.5d0, oldrand, jrand)
         end do
      end do
      ncross = 0
      nmute = 0
      call ga_crossover(1, 2, 4, oldchr, newchr, ld, lc,
     &   1.d0, 0.d0, ncross, nmute, oldrand, jrand)
      ok = 1
      do k = 1, lc
         if (newchr(1,k) + newchr(2,k) .ne.
     &       oldchr(2,k) + oldchr(4,k)) ok = 0
      end do
      call ck_i('C2 crossover conserves parent bits', ok, 1,
     &          npass, nfail)
      call ck_i('C3 crossover counted (pcross=1)', ncross, 1,
     &          npass, nfail)
      call ck_i('C4 no mutation at pmute=0', nmute, 0, npass, nfail)

c     C5: ga_preselect terminates and fills the pool with valid
c     indices when every fitness is 0 (all fractions 0)
      do i = 1, 8
         fit(i) = 0.d0
         ich(i) = -1
      end do
      call ga_preselect(fit, 8, 1.d-12, ich, frac, oldrand, jrand)
      ok = 1
      do i = 1, 8
         if (ich(i) .lt. 1 .or. ich(i) .gt. 8) ok = 0
      end do
      call ck_i('C5 preselect flat population: valid pool', ok, 1,
     &          npass, nfail)

c     C6: ga_select2 always returns an index from the pool. (INFO:
c     the pool range is never shrunk by the callers, so draws are
c     with replacement despite the "without replacement" comment.)
      do i = 1, 8
         ich(i) = i
      end do
      ok = 1
      nseen = 0
      do j = 1, 8
         call ga_select2(ich, 8, oldrand, jrand, k)
         if (k .lt. 1 .or. k .gt. 8) ok = 0
         fit(j) = dble(k)
      end do
      call ck_i('C6 select2 returns pool member', ok, 1, npass, nfail)
      do i = 1, 8
         do j = 1, 8
            if (nint(fit(j)) .eq. i) then
               nseen = nseen + 1
               goto 61
            end if
         end do
 61      continue
      end do
      write(*,'(A,I2,A)') '  [INFO]  select2: 8 draws hit ', nseen,
     &   ' distinct of 8 (a permutation would hit 8)'

c     C7: ga_slope_fitness exact match = 1, offset d gives 1/(1+d)
      do i = 1, 3
         eg(i) = 0.1d0*dble(i)
         bg(i) = 10.d0*dble(i)
         do j = 1, 3
            sl(i, j) = -2.d0
         end do
      end do
      xv(1) = 0.2d0
      xv(2) = 20.d0
      call ck_d('C7 slope fitness at match = 1',
     &   ga_slope_fitness(xv, sl, eg, bg, -2.d0), 1.d0, 1.d-14,
     &   npass, nfail)
      call ck_d('C8 slope fitness offset 3 = 1/4',
     &   ga_slope_fitness(xv, sl, eg, bg, 1.d0), 0.25d0, 1.d-14,
     &   npass, nfail)

c     C9: a NaN slope surface must give fitness 0, not NaN (the
c     same issue #30 guard ga_fitness already has)
      z1 = 0.d0
      z2 = 0.d0
      sl(2, 2) = z1 / z2
      f = ga_slope_fitness(xv, sl, eg, bg, -2.d0)
      ok = 0
      if (f .eq. f .and. f .ge. 0.d0 .and. f .le. 1.d0) ok = 1
      call ck_i('C9 slope fitness NaN-guarded', ok, 1, npass, nfail)

      return
      end

c     ============================================================
c     ga_optimize with no variation (pcross=pmute=0): children are
c     exact copies of parents, so the reported optimum must be one
c     of the initial individuals. Fails if ga_crossover addresses
c     the population arrays with the wrong leading dimension.
c     ============================================================
      subroutine t_ga_layout(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail
      parameter(np = 6)
      double precision surf(np, 3, 3), eg(3), ag(3), tob(np), apo(np)
      integer ipops(2), ip
      character*60 lab

      ipops(1) = 20
      ipops(2) = 4
      do i = 1, 3
         eg(i) = 0.30d0 + 0.02d0*dble(i-1)
         ag(i) = 80.d0 + 40.d0*dble(i-1)
      end do
      do k = 1, 3
         do j = 1, 3
            do i = 1, np
               t = dble(i-1)/dble(np-1)
               surf(i, j, k) = 12000.d0 - 100.d0*t
     &            - 50.d0*(eg(j)-0.32d0) - 0.5d0*(ag(k)-120.d0)
            end do
         end do
      end do
      do i = 1, np
         tob(i) = dble(i-1)/dble(np-1)
         apo(i) = 12000.d0 - 100.d0*tob(i)
      end do

      do ip = 1, 2
         call ga_optimize(surf, np, eg, ag, 3, 3, np,
     &      tob, apo, np, 0.30d0, 0.34d0, 80.d0, 160.d0,
     &      ipops(ip), 30, 20, 12, 0.d0, 0.d0, 0.777d0,
     &      e_opt, a_opt, rms_opt)
         call in_init_pop(ipops(ip), 0.777d0, 20, 12,
     &      0.30d0, 0.34d0, 80.d0, 160.d0, e_opt, a_opt, ihit)
         write(lab, '(A,I2,A)') 'C10 popsize ', ipops(ip),
     &      ': optimum from initial pop'
         call ck_i(lab, ihit, 1, npass, nfail)
      end do

c     C10b: the GA must actually search at any population size: after
c     300 generations with crossover and mutation the best point must
c     no longer be one of the random initial individuals.
c     Unique optimum near the upper corner (e=0.338, a=155): the
c     residual is -0.5(a-155) at t=0 and -500(e-0.338) at t=1.
      do k = 1, 3
         do j = 1, 3
            do i = 1, np
               t = dble(i-1)/dble(np-1)
               surf(i, j, k) = 12000.d0 - 100.d0*t
     &            - 500.d0*(eg(j)-0.338d0)*t
     &            - 0.5d0*(ag(k)-155.d0)*(1.d0-t)
            end do
         end do
      end do
      do ip = 1, 2
         call ga_optimize(surf, np, eg, ag, 3, 3, np,
     &      tob, apo, np, 0.30d0, 0.34d0, 80.d0, 160.d0,
     &      ipops(ip), 300, 20, 12, 0.8d0, 0.02d0, 0.555d0,
     &      e_opt, a_opt, rms_opt)
         write(*,'(A,I3,A,F9.5,A,F8.2,A,F8.3)') '  [INFO]  popsize',
     &      ipops(ip), ': e=', e_opt, ' a=', a_opt, ' rms=', rms_opt
         write(lab, '(A,I2,A)') 'C10b popsize ', ipops(ip),
     &      ': searches beyond initial pop'
         call in_init_pop(ipops(ip), 0.555d0, 20, 12,
     &      0.30d0, 0.34d0, 80.d0, 160.d0, e_opt, a_opt, ihit)
         call ck_i(lab, ihit, 0, npass, nfail)
      end do

      return
      end

c     Rebuild ga_optimize's initial population and report whether
c     (e,a) equals one of its decoded individuals.
      subroutine in_init_pop(ipop, seed, nbe, nba, elo, ehi,
     &                       alo, ahi, e, a, ihit)
      implicit double precision (a-h, o-z)
      integer ipop, nbe, nba, ihit
      double precision oldrand(55), xv(2), alow(2), ahigh(2), fac(2)
      integer chr(200), lsub(2), jrand, ga_iflip

      alow(1) = elo
      ahigh(1) = ehi
      alow(2) = alo
      ahigh(2) = ahi
      lsub(1) = nbe
      lsub(2) = nba
      fac(1) = 2.d0**dble(nbe) - 1.d0
      fac(2) = 2.d0**dble(nba) - 1.d0
      call ga_warmup(seed, oldrand, jrand)
      ihit = 0
      do j = 1, ipop
         do k = 1, nbe + nba
            chr(k) = ga_iflip(0.5d0, oldrand, jrand)
         end do
         call ga_decode(chr, nbe+nba, xv, 2, lsub, alow, ahigh, fac)
         if (dabs(xv(1)-e) .le. 1.d-14 .and.
     &       dabs(xv(2)-a) .le. 1.d-10) ihit = 1
      end do
      return
      end

c     ============================================================
c     ga_campaign / ga_camp_fitness on a synthetic two-zone case with
c     a unique BN optimum (surfaces linear in BN, so the bilinear
c     interpolation is exact).
c     ============================================================
      subroutine t_ga_campaign(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail
      parameter(lds = 20, lz = 3)
      double precision as(lds, 3, 3, lz), aeg(3, lz), abg(3, lz)
      double precision aem(lz), atob(500, lz), aapo(500, lz)
      integer ansurf(lz), anobs(lz), ok
      double precision ga_camp_fitness

      bn_true = 60.d0
      do iz = 1, 2
         ansurf(iz) = 10
         anobs(iz) = 10
         aem(iz) = 0.5d0
         do i = 1, 3
            aeg(i, iz) = 0.4d0 + 0.1d0*dble(i-1)
            abg(i, iz) = 20.d0 + 40.d0*dble(i-1)
         end do
         do it = 1, 10
            t = dble(it-1)
            atob(it, iz) = t
            aapo(it, iz) = 20000.d0 - t*(200.d0 - bn_true)
     &                   - 1000.d0*dble(iz)
            do ie = 1, 3
               do ib = 1, 3
                  as(it, ie, ib, iz) = 20000.d0
     &               - t*(200.d0 - abg(ib, iz)) - 1000.d0*dble(iz)
               end do
            end do
         end do
      end do

      f = ga_camp_fitness(bn_true, as, lds, lz, aeg, abg, aem,
     &   ansurf, atob, aapo, anobs, 2)
      call ck_d('C11 camp fitness at true BN = 1', f, 1.d0, 1.d-12,
     &          npass, nfail)

c     zone 3 has no surface: skipped, does not change the mean
      ansurf(3) = 0
      anobs(3) = 5
      f3 = ga_camp_fitness(bn_true + 10.d0, as, lds, lz, aeg, abg,
     &   aem, ansurf, atob, aapo, anobs, 3)
      f2 = ga_camp_fitness(bn_true + 10.d0, as, lds, lz, aeg, abg,
     &   aem, ansurf, atob, aapo, anobs, 2)
      call ck_d('C12 camp fitness skips empty zone', f3, f2, 1.d-14,
     &          npass, nfail)

      call ga_campaign(as, lds, lz, aeg, abg, aem, ansurf,
     &   atob, aapo, anobs, 2, 20.d0, 100.d0,
     &   20, 200, 16, 0.8d0, 0.01d0, 0.4242d0, bn_opt, rms_opt)
      call ck_d('C13 ga_campaign recovers BN=60', bn_opt, bn_true,
     &          1.d0, npass, nfail)
      call ck_d('C14 ga_campaign rms ~ 0 km', rms_opt, 0.d0, 3.d0,
     &          npass, nfail)

c     C15: NaN in one zone's surface must not give a NaN fitness
      z1 = 0.d0
      z2 = 0.d0
      as(3, 2, 2, 1) = z1 / z2
      f = ga_camp_fitness(bn_true, as, lds, lz, aeg, abg, aem,
     &   ansurf, atob, aapo, anobs, 2)
      ok = 0
      if (f .eq. f .and. f .ge. 0.d0 .and. f .le. 1.d0) ok = 1
      call ck_i('C15 camp fitness NaN-guarded', ok, 1, npass, nfail)

      return
      end

c     ============================================================
      subroutine t_ga_slope(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail
      double precision sl(3, 3), eg(3), bg(3), xv(2)
      double precision ga_slope_fitness

      do i = 1, 3
         eg(i) = 0.60d0 + 0.01d0*dble(i-1)
         bg(i) = 20.d0 + 40.d0*dble(i-1)
      end do
c     slope independent of e, linear in BN: unique optimum BN=80
      do i = 1, 3
         do j = 1, 3
            sl(i, j) = -10.d0 + 0.05d0*bg(j)
         end do
      end do
      obs = -10.d0 + 0.05d0*80.d0
      call ga_slope_optimize(sl, eg, bg, obs,
     &   eg(1), eg(3), bg(1), bg(3),
     &   20, 200, 12, 16, 0.8d0, 0.01d0, 0.31d0,
     &   e_opt, bn_opt, err_opt)
      call ck_d('C16 slope GA recovers BN=80', bn_opt, 80.d0, 1.d0,
     &          npass, nfail)
      call ck_d('C17 slope GA error ~ 0', err_opt, 0.d0, 0.05d0,
     &          npass, nfail)
      xv(1) = e_opt
      xv(2) = bn_opt
      f = ga_slope_fitness(xv, sl, eg, bg, obs)
      call ck_d('C18 reported error matches fitness',
     &          (1.d0 - f)/f, err_opt, 1.d-12, npass, nfail)
      return
      end

c     ============================================================
      subroutine t_tinterp(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail
      double precision ta(4), ya(4)

      t0 = 2460000.d0
      do i = 1, 4
         ta(i) = t0 + dble(i-1)
         ya(i) = 100.d0 - 10.d0*dble(i-1)
      end do
      call rsm_tinterp(ta, ya, 4, t0, 1.5d0, y)
      call ck_d('C19 tinterp midpoint', y, 85.d0, 1.d-9, npass, nfail)
      call rsm_tinterp(ta, ya, 4, t0, 2.d0, y)
      call ck_d('C20 tinterp at node', y, 80.d0, 1.d-9, npass, nfail)
      call rsm_tinterp(ta, ya, 4, t0, -3.d0, y)
      call ck_d('C21 tinterp clamps low', y, 100.d0, 0.d0,
     &          npass, nfail)
      call rsm_tinterp(ta, ya, 4, t0, 9.d0, y)
      call ck_d('C22 tinterp clamps high', y, 70.d0, 0.d0,
     &          npass, nfail)
      call rsm_tinterp(ta, ya, 1, t0, 9.d0, y)
      call ck_d('C23 tinterp n=1', y, 100.d0, 0.d0, npass, nfail)
      return
      end

c     ============================================================
      subroutine t_rpe(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail
      double precision rj(4), ze(4), rpe(4)

      do i = 1, 4
         ze(i) = 1000.d0 + 10.d0*dble(i-1)
      end do
      rj(1) = 1100.d0
      rj(2) = 0.d0
      rj(3) = 1120.d0
      rj(4) = 1090.d0
c     Mode 1: vs observed 1100
      call compute_rpe(rj, ze, 4, 1100.d0, rpe, tm, ts, 4)
      call ck_d('C24 median of 3 valid', tm, 1100.d0, 0.d0,
     &          npass, nfail)
      call ck_d('C25 Mode-1 RPE zone 3', rpe(3),
     &          20.d0/80.d0*100.d0, 1.d-12, npass, nfail)
      call ck_d('C26 invalid zone RPE = 0', rpe(2), 0.d0, 0.d0,
     &          npass, nfail)
      call ck_d('C27 std about median', ts,
     &          dsqrt((0.d0 + 400.d0 + 100.d0)/2.d0), 1.d-12,
     &          npass, nfail)
c     even count median, Mode 2
      rj(2) = 1110.d0
      call compute_rpe(rj, ze, 4, 0.d0, rpe, tm, ts, 4)
      call ck_d('C28 median of 4 valid', tm, 1105.d0, 0.d0,
     &          npass, nfail)
      call ck_d('C29 Mode-2 RPE zone 1', rpe(1),
     &          -5.d0/105.d0*100.d0, 1.d-12, npass, nfail)
      return
      end

c     ============================================================
c     linfit on JD-scale abscissae: exact straight line, short dense
c     window. Uncentred sums lose ~2-4% of the slope here.
c     ============================================================
      subroutine t_linfit(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail
      double precision x(5), y(5)

      do i = 1, 5
         x(i) = 2460000.3d0 + 0.125d0*dble(i-1)
         y(i) = 30000.d0 - 10.d0*(x(i) - x(1))
      end do
      call linfit(x, y, 5, s, c, r2)
      write(*,'(A,F12.6,A)') '  [INFO]  linfit slope ', s,
     &   ' (exact -10) on a 0.5-day JD window'
      call ck_known('C30 linfit slope, 0.5-day JD window', s, -10.d0,
     &          1.d-6, npass, nfail)
      call ck_known('C31 linfit R2 of exact line', r2, 1.d0, 1.d-9,
     &          npass, nfail)
      return
      end

c     ============================================================
      subroutine t_zone_next(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail
      parameter(n = 60)
      double precision ep(n), ha(n)
      integer zs, ze, zs2, ze2, ok
      logical found

      do i = 1, n
         ep(i) = 2460000.5d0 + dble(i-1)
         if (i .le. 20) then
            ha(i) = 30000.d0
         else if (i .le. 50) then
            ha(i) = 30000.d0 - 5.d0*dble(i-20)
         else
            ha(i) = 30000.d0 - 150.d0
         end if
      end do

      call zone_select_next(ep, ha, n, 1, 5, 20.d0, 0.9d0, -0.5d0,
     &                      zs, ze, found)
      ok = 0
      if (found .and. zs .ge. 15 .and. zs .le. 21 .and.
     &    ze - zs + 1 .ge. 5 .and. ep(ze) - ep(zs) .le. 20.d0)
     &   ok = 1
      call ck_i('C32 finds the decay segment', ok, 1, npass, nfail)
      write(*,'(A,I3,A,I3)') '  [INFO]  zone ', zs, ' -', ze

      call zone_select_next(ep, ha, n, ze + 1, 5, 20.d0, 0.9d0,
     &                      -0.5d0, zs2, ze2, found)
      ok = 0
      if (found .and. zs2 .gt. ze) ok = 1
      if (.not. found) ok = 1
      call ck_i('C33 next zone does not overlap', ok, 1, npass, nfail)

      call zone_select_next(ep, ha, n, 52, 5, 20.d0, 0.9d0, -0.5d0,
     &                      zs, ze, found)
      ok = 1
      if (found) ok = 0
      call ck_i('C34 flat tail: no zone', ok, 1, npass, nfail)

      call zone_select_next(ep, ha, 4, 1, 5, 20.d0, 0.9d0, -0.5d0,
     &                      zs, ze, found)
      ok = 1
      if (found) ok = 0
      call ck_i('C35 npts < min_zone_pts: no zone', ok, 1,
     &          npass, nfail)
      return
      end

c     ============================================================
c     TLE lookup helpers on hand-built raw arrays. Two TLEs 5 min
c     apart are kept as separate points by tle_evolve (>=86 s), so
c     each lookup must return the record for ITS epoch, not just the
c     first record within the 0.01-day tolerance.
c     ============================================================
      subroutine t_tle_helpers(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail
      parameter(mt = 4)
      integer nr(mt), iyr(mt), ok
      double precision eday(mt), ama(mt), bst(mt)
      double precision ainc(mt), raan(mt), ecc(mt), aop(mt), amm(mt)
      double precision epz(3), akep(7)
      double precision sm(3), ez(3), az(3)
      logical lok
      integer nf(3)
      double precision alt(1), den(1), sch(1)

      do i = 1, mt
         nr(i) = 11111
         iyr(i) = 2020
         ainc(i) = 63.4d0
         raan(i) = 100.d0
         ecc(i) = 0.7d0
         aop(i) = 270.d0
         amm(i) = 2.006d0
      end do
      nr(4) = 22222
      eday(1) = 100.50d0
      eday(2) = 100.50d0 + 5.d0/1440.d0
      eday(3) = 101.50d0
      eday(4) = 100.50d0 + 5.d0/1440.d0
      ama(1) = 10.d0
      ama(2) = 20.d0
      ama(3) = 30.d0
      ama(4) = 99.d0
      bst(1) = 1.d-4
      bst(2) = 2.d-4
      bst(3) = 1.d-4
      bst(4) = 9.d-4

      call tle_epoch2jd(2020, eday(2), ep2)
      call tle_find_ma(11111, ep2, nr, iyr, eday, ama, mt, mt, xma)
      call ck_known('C36 find_ma: exact epoch, not first in tolerance',
     &          xma, 20.d0, 0.d0, npass, nfail)
      call tle_find_bstar(11111, ep2, nr, iyr, eday, bst, mt, mt, b)
      call ck_known('C37 find_bstar: exact epoch', b, 2.d-4, 0.d0,
     &          npass, nfail)
      call tle_epoch2jd(2020, eday(3), ep3)
      call tle_find_ma(11111, ep3, nr, iyr, eday, ama, mt, mt, xma)
      call ck_d('C38 find_ma: isolated epoch', xma, 30.d0, 0.d0,
     &          npass, nfail)
      call tle_find_ma(33333, ep3, nr, iyr, eday, ama, mt, mt, xma)
      call ck_d('C39 find_ma: unknown object -> 0', xma, 0.d0, 0.d0,
     &          npass, nfail)

c     tle_find_osc: SGP4 state for record 3 (Molniya-like)
      call tle_find_osc(11111, ep3, nr, iyr, eday, ainc, raan, ecc,
     &   aop, ama, amm, bst, mt, mt, akep, lok)
      ok = 0
      if (lok .and. akep(2) .gt. 0.6d0 .and. akep(2) .lt. 0.8d0 .and.
     &    dabs(akep(3) - 63.4d0) .lt. 1.d0) ok = 1
      call ck_i('C40 find_osc: osculating e, i near mean', ok, 1,
     &          npass, nfail)
      call tle_find_osc(11111, ep3 + 0.5d0, nr, iyr, eday, ainc,
     &   raan, ecc, aop, ama, amm, bst, mt, mt, akep, lok)
      ok = 1
      if (lok) ok = 0
      call ck_i('C41 find_osc: no record -> .false.', ok, 1,
     &          npass, nfail)

c     estimate_bn_bstar_prior: BSTAR 1e-4 -> centre
c     10**(-0.4105*(-4)+0.5782) = 166.0, window x/ 10**(2*0.2689)
      call tle_epoch2jd(2020, eday(1), epz(1))
      epz(2) = ep3
      call estimate_bn_bstar_prior(11111, epz, 2, nr, iyr, eday,
     &   bst, mt, mt, 10.d0, 160.d0, blo, bhi, lok)
      ok = 0
      if (lok .and. dabs(blo - 10.d0**(2.2202d0 - 0.5378d0)) .lt.
     &    1.d-6*blo .and. bhi .eq. 160.d0) ok = 1
      call ck_i('C42 BSTAR prior narrows the low side', ok, 1,
     &          npass, nfail)
      call estimate_bn_bstar_prior(11111, epz, 2, nr, iyr, eday,
     &   bst, mt, mt, 600.d0, 900.d0, blo, bhi, lok)
      ok = 0
      if (.not. lok .and. blo .eq. 600.d0 .and. bhi .eq. 900.d0)
     &   ok = 1
      call ck_i('C43 no overlap -> caller range kept', ok, 1,
     &          npass, nfail)
      call estimate_bn_bstar_prior(44444, epz, 2, nr, iyr, eday,
     &   bst, mt, mt, 10.d0, 160.d0, blo, bhi, lok)
      ok = 0
      if (.not. lok .and. blo .eq. 10.d0 .and. bhi .eq. 160.d0)
     &   ok = 1
      call ck_i('C44 no BSTAR -> caller range kept', ok, 1,
     &          npass, nfail)

c     estimate_bn_floor early exits (no propagation reached)
      nf(1) = 1
      nf(2) = 0
      nf(3) = 0
      do i = 1, 3
         epz(i) = 2460000.d0 + dble(i)
         sm(i) = 26000.d0
         ez(i) = 0.7d0
         az(i) = 0.d0
      end do
      call estimate_bn_floor(epz, sm, ez, az, az, az, az, 2,
     &   nf, 2, 0, 0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0, 0, 0.d0,
     &   0.d0, 0.d0, alt, den, sch, 1, bfl, lok)
      ok = 1
      if (lok) ok = 0
      call ck_i('C45 bn_floor: <3 TLEs -> not ok', ok, 1,
     &          npass, nfail)
      call estimate_bn_floor(epz, sm, ez, az, az, az, az, 3,
     &   nf, 2, 0, 0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0, 0, 0.d0,
     &   0.d0, 0.d0, alt, den, sch, 1, bfl, lok)
      ok = 1
      if (lok) ok = 0
      call ck_i('C46 bn_floor: flat a(t) -> not ok', ok, 1,
     &          npass, nfail)
      return
      end

c     ============================================================
c     n_force(1)=0 must switch the geopotential off: with no drag,
c     no third body, the run must equal one with ngeo_deg=0.
c     ============================================================
      subroutine t_force_switch(npass, nfail)
      implicit double precision (a-h, o-z)
      integer npass, nfail
      parameter(ndim_a = 400, mp = 10)
      double precision alt(ndim_a), den(ndim_a), sch(ndim_a)
      double precision x0(3), xd0(3), cal0(6), ak(7)
      double precision tj1(mp), tx1(3,mp), tv1(3,mp)
      double precision tj2(mp), tx2(3,mp), tv2(3,mp)
      integer nf(3), natm, id1, id2, ie1, ie2, ok

      call read_atm2(alt, den, sch, ndim_a, natm)
      ak(1) = 7500.d0
      ak(2) = 0.02d0
      ak(3) = 51.6d0
      ak(4) = 30.d0
      ak(5) = 40.d0
      ak(6) = 50.d0
      ak(7) = 0.d0
      call oe2car(ak, x0, xd0, 1.d-12)
      cal0(1) = 2024.d0
      cal0(2) = 3.d0
      cal0(3) = 1.d0
      cal0(4) = 0.d0
      cal0(5) = 0.d0
      cal0(6) = 0.d0

      nf(1) = 0
      nf(2) = 0
      nf(3) = 0
      id1 = 1
      call propagate_ks(x0, xd0, cal0, 3, 360, 1.d-10,
     &   nf, 4, 0, 0, 100.d0, 0, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0, 0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm, mp, id1, tj1, tx1, tv1, ie1)
      nf(1) = 1
      id2 = 1
      call propagate_ks(x0, xd0, cal0, 3, 360, 1.d-10,
     &   nf, 0, 0, 0, 100.d0, 0, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0, 0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm, mp, id2, tj2, tx2, tv2, ie2)
      d = 0.d0
      do i = 1, 3
         d = max(d, dabs(tx1(i,id1) - tx2(i,id2)))
      end do
      write(*,'(A,ES10.3,A,4I4)') '  [INFO]  n_force(1)=0,ngeo=4 vs '//
     &   'ngeo=0: max |dr| = ', d, ' km; idump/exit ', id1, ie1,
     &   id2, ie2
      call ck_i('C47a both runs complete 3 rev', id1 + id2, 8,
     &          npass, nfail)
      call ck_d('C47 n_force(1)=0 disables geopotential', d, 0.d0,
     &          1.d-6, npass, nfail)
c     positive control: geopotential on (n_force(1)=1, ngeo=4)
      id1 = 1
      call propagate_ks(x0, xd0, cal0, 3, 360, 1.d-10,
     &   nf, 4, 0, 0, 100.d0, 0, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0, 0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm, mp, id1, tj1, tx1, tv1, ie1)
      d = 0.d0
      do i = 1, 3
         d = max(d, dabs(tx1(i,id1) - tx2(i,id2)))
      end do
      write(*,'(A,ES10.3,A)') '  [INFO]  ngeo=4 on vs off: max |dr| = ',
     &   d, ' km'
c     KNOWN-BUG while input/EGM2008_to2190_TideFree is the 0-byte
c     file committed in c820803: geo_coeff then reads no coefficients
c     and every zonal term (J2 included) is zero. C47 is only a real
c     check once the file has data.
      ok = 0
      if (d .gt. 1.d-2) ok = 1
      call ck_known('C47b geopotential on changes the orbit',
     &          dble(ok), 1.d0, 0.d0, npass, nfail)
      return
      end

c     ============================================================
c     rsm_direct_eval / rsm_refine on a drag-free orbit (perigee
c     ~970 km, above the 500 km table): apogee stays within a few km
c     of its initial value over one day.
c     ============================================================
      subroutine t_direct_refine(npass, nfail)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth
      integer npass, nfail
      parameter(ndim_a = 400, nob = 5)
      double precision alt(ndim_a), den(ndim_a), sch(ndim_a)
      double precision tob(nob), apo(nob), apo2(nob)
      integer nf(3), natm, ierr, nev, ok

      call read_atm2(alt, den, sch, ndim_a, natm)
      nf(1) = 1
      nf(2) = 0
      nf(3) = 0
      a1 = 7500.d0
      e1 = 0.02d0
      ha0 = a1*(1.d0 + e1) - R_Earth
      ep0 = 2460370.5d0
      do i = 1, nob
         tob(i) = 0.25d0*dble(i-1)
         apo(i) = ha0
         apo2(i) = ha0 + 100.d0
      end do

      call rsm_direct_eval(a1, e1, 51.6d0, 30.d0, 40.d0, 50.d0,
     &   ep0, 1.d0, e1, 100.d0, tob, apo, nob, 360, 1.d-10,
     &   nf, 2, 0, 0, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0, 0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm, rms1, ierr)
      call ck_i('C48 direct_eval ierr=0', ierr, 0, npass, nfail)
      ok = 0
      if (rms1 .ge. 0.d0 .and. rms1 .lt. 30.d0) ok = 1
      call ck_i('C49 direct_eval rms vs own apogee < 30 km', ok, 1,
     &          npass, nfail)
      call rsm_direct_eval(a1, e1, 51.6d0, 30.d0, 40.d0, 50.d0,
     &   ep0, 1.d0, e1, 100.d0, tob, apo2, nob, 360, 1.d-10,
     &   nf, 2, 0, 0, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0, 0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm, rms2, ierr)
      call ck_d('C50 direct_eval 100 km offset shows up', rms2,
     &          100.d0, 30.d0, npass, nfail)

c     rsm_refine: result inside bounds, never worse than the start,
c     evaluation budget respected
      call rsm_refine(a1, e1, 51.6d0, 30.d0, 40.d0, 50.d0,
     &   ep0, 1.d0, e1 + 0.004d0, 100.d0,
     &   e1 - 0.005d0, e1 + 0.005d0, 50.d0, 150.d0,
     &   tob, apo, nob, 360, 1.d-10,
     &   nf, 2, 0, 0, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0, 0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm, 9, er, br, rr, nev)
      call rsm_direct_eval(a1, e1, 51.6d0, 30.d0, 40.d0, 50.d0,
     &   ep0, 1.d0, e1 + 0.004d0, 100.d0, tob, apo, nob, 360, 1.d-10,
     &   nf, 2, 0, 0, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0, 0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm, rms_st, ierr)
      ok = 0
      if (nev .ge. 1 .and. nev .le. 9 .and. rr .le. rms_st .and.
     &    er .ge. e1 - 0.005d0 .and. er .le. e1 + 0.005d0 .and.
     &    br .ge. 50.d0 .and. br .le. 150.d0) ok = 1
      call ck_i('C51 refine: in bounds, no worse, budget kept', ok,
     &          1, npass, nfail)
      ok = 0
      if (rr .lt. rms_st) ok = 1
      call ck_i('C52 refine improves an off-optimum start', ok, 1,
     &          npass, nfail)

c     Degenerate ranges (e_lo=e_hi, bn_lo=bn_hi): must terminate
      call rsm_refine(a1, e1, 51.6d0, 30.d0, 40.d0, 50.d0,
     &   ep0, 1.d0, e1, 100.d0, e1, e1, 100.d0, 100.d0,
     &   tob, apo, nob, 360, 1.d-10,
     &   nf, 2, 0, 0, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   0.d0, 0.d0, 0, 0, 0.d0, 0.d0, 0.d0,
     &   alt, den, sch, natm, 9, er, br, rr, nev)
      call ck_i('C53 refine terminates on a zero-width range', nev,
     &          1, npass, nfail)
      return
      end

c     ============================================================
      subroutine read_atm2(alt, den, sch, maxn, natm)
      implicit double precision (a-h, o-z)
      integer maxn, natm
      double precision alt(maxn), den(maxn), sch(maxn)
      natm = 291
      alt(1) = 60.d0
      do i = 2, 141
         alt(i) = alt(i-1) + 1.d0
      end do
      do i = 142, natm
         alt(i) = alt(i-1) + 2.d0
      end do
      open(unit=88, file='input/ATM.DAT', status='old', err=900)
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
 900  write(*,*) 'read_atm2: cannot open input/ATM.DAT'
      stop 1
      end

c     ============================================================
      subroutine ck_i(label, got, expect, npass, nfail)
      character*(*) label
      integer got, expect, npass, nfail
      if (got .eq. expect) then
         write(*,'(A,A)') '  [PASS]  ', label
         npass = npass + 1
      else
         write(*,'(A,A)') '  [FAIL]  ', label
         write(*,'(A,I10,A,I10)') '           got=', got,
     &        '  expect=', expect
         nfail = nfail + 1
      end if
      return
      end

      subroutine ck_d(label, got, expect, tol, npass, nfail)
      character*(*) label
      double precision got, expect, tol
      integer npass, nfail
      if (dabs(got - expect) .le. tol) then
         write(*,'(A,A)') '  [PASS]  ', label
         npass = npass + 1
      else
         write(*,'(A,A)') '  [FAIL]  ', label
         write(*,'(A,E20.12,A,E20.12)') '           got=', got,
     &        '  expect=', expect
         nfail = nfail + 1
      end if
      return
      end

c     Check for a reported-but-unfixed bug: a failure is printed as
c     [KNOWN-BUG] and does not fail the suite; a pass means it was
c     fixed and the call can become ck_d.
      subroutine ck_known(label, got, expect, tol, npass, nfail)
      character*(*) label
      double precision got, expect, tol
      integer npass, nfail
      if (dabs(got - expect) .le. tol) then
         write(*,'(A,A)') '  [PASS]  ', label
         npass = npass + 1
      else
         write(*,'(A,A)') '  [KNOWN-BUG]  ', label
         write(*,'(A,E20.12,A,E20.12)') '           got=', got,
     &        '  expect=', expect
      end if
      return
      end

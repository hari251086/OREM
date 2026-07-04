c     ============================================================
c     ga.F — Binary-coded Genetic Algorithm for OREM
c
c     Refactored from Sharma's GENESIS GA (genen1.f).
c     Changes: real*4 → double precision, COMMON blocks removed,
c     file I/O removed, callable subroutine interface.
c
c     Algorithm: stochastic remainder selection, single-point
c     crossover, bit-flip mutation, 55-element Fibonacci RNG.
c
c     Author: Harishkumar Sellamuthu
c     ============================================================

c     ============================================================
      subroutine ga_optimize(
     &   surfaces, ld_surf, e_grid, a_grid,
     &   nsurf_e, nsurf_a, nsurf_pts,
     &   tobs, apobs, nobs,
     &   e_lo, e_hi, a_lo, a_hi,
     &   ipopsize, maxgen, nbits_e, nbits_a,
     &   pcross, pmute, seed,
     &   e_opt, a_opt, rms_opt)
c
c     Binary-coded GA searching pre-computed RSM surfaces.
c
c     Input
c       surfaces(ld_surf, nsurf_e, nsurf_a) : mean apogee
c       ld_surf            : leading dimension of surfaces array
c       e_grid(nsurf_e)    : eccentricity grid values
c       a_grid(nsurf_a)    : area grid values
c       nsurf_e, nsurf_a   : grid dimensions (typically 3)
c       nsurf_pts          : time points per surface (<=ld_surf)
c       tobs(nobs)         : observation times (days from zone start)
c       apobs(nobs)        : observed TLE mean apogee (km)
c       nobs               : number of observations
c       e_lo, e_hi         : eccentricity search bounds
c       a_lo, a_hi         : area search bounds (m²)
c       ipopsize           : population size (e.g. 4)
c       maxgen             : max generations (e.g. 500)
c       nbits_e            : bits for eccentricity (e.g. 40)
c       nbits_a            : bits for area (e.g. 16-28)
c       pcross             : crossover probability (e.g. 0.8)
c       pmute              : mutation probability (e.g. 0.01)
c       seed               : random seed (e.g. 0.123)
c
c     Output
c       e_opt              : optimal eccentricity
c       a_opt              : optimal area
c       rms_opt            : RMS fitness at optimum
c     ============================================================
      implicit double precision (a-h, o-z)
      integer nsurf_e, nsurf_a, nsurf_pts, ld_surf, nobs
      integer ipopsize, maxgen, nbits_e, nbits_a
      double precision surfaces(ld_surf, nsurf_e, nsurf_a)
      double precision e_grid(nsurf_e), a_grid(nsurf_a)
      double precision tobs(nobs), apobs(nobs)
      double precision e_lo, e_hi, a_lo, a_hi
      double precision pcross, pmute, seed
      double precision e_opt, a_opt, rms_opt

c     Internal parameters
      parameter(maxpop = 20, maxbits = 200)
      integer nparam
      parameter(nparam = 2)

      integer lchrom
      integer oldchr(maxpop, maxbits)
      integer newchr(maxpop, maxbits)
      integer chr(maxbits)
      integer lsubstr(nparam), choices(maxpop)
      double precision oldx(maxpop, nparam), newx(maxpop, nparam)
      double precision oldfit(maxpop), newfit(maxpop)
      double precision fraction(maxpop)
      double precision factor(nparam), xvec(nparam)
      double precision alow(nparam), ahigh(nparam)
      double precision bestsol(nparam)
      double precision bestfit
      integer ibest, igenbst

c     RNG state
      double precision oldrand(55)
      integer jrand

c     Working
      integer i, j, k, igen, ncross, nmute
      double precision sumfitness, avg, amax, amin

      lchrom = nbits_e + nbits_a
      alow(1) = e_lo
      ahigh(1) = e_hi
      alow(2) = a_lo
      ahigh(2) = a_hi
      lsubstr(1) = nbits_e
      lsubstr(2) = nbits_a

      do i = 1, nparam
         factor(i) = 2.d0**dble(lsubstr(i)) - 1.d0
      end do

      ncross = 0
      nmute = 0

c     Initialize RNG
      call ga_warmup(seed, oldrand, jrand)

c     Create initial population
      do j = 1, ipopsize
         do k = 1, lchrom
            chr(k) = ga_iflip(0.5d0, oldrand, jrand)
            oldchr(j, k) = chr(k)
         end do
         call ga_decode(chr, lchrom, xvec, nparam,
     &                  lsubstr, alow, ahigh, factor)
         do k = 1, nparam
            oldx(j, k) = xvec(k)
         end do
         oldfit(j) = ga_fitness(xvec,
     &      surfaces, ld_surf, e_grid, a_grid,
     &      nsurf_e, nsurf_a, nsurf_pts,
     &      tobs, apobs, nobs)
      end do

c     Initial statistics
      call ga_stats(oldfit, ipopsize,
     &              sumfitness, avg, amax, amin, ibest)
      bestfit = oldfit(ibest)
      do j = 1, nparam
         bestsol(j) = oldx(ibest, j)
      end do
      igenbst = 0

c     === Main GA loop ===
      do igen = 1, maxgen

c        Preselect mating pool
         call ga_preselect(oldfit, ipopsize, avg,
     &                     choices, fraction,
     &                     oldrand, jrand)

c        Create new population via crossover + mutation
         j = 1
 10      continue
         if (j .gt. ipopsize) goto 20

c        Select two mates
         call ga_select2(choices, ipopsize,
     &                   oldrand, jrand, mate1)
         call ga_select2(choices, ipopsize,
     &                   oldrand, jrand, mate2)

c        Crossover + mutation
         call ga_crossover(j, mate1, mate2,
     &      oldchr, newchr, ipopsize, lchrom,
     &      pcross, pmute, ncross, nmute,
     &      oldrand, jrand)

c        Decode and evaluate child 1
         do k = 1, lchrom
            chr(k) = newchr(j, k)
         end do
         call ga_decode(chr, lchrom, xvec, nparam,
     &                  lsubstr, alow, ahigh, factor)
         do k = 1, nparam
            newx(j, k) = xvec(k)
         end do
         newfit(j) = ga_fitness(xvec,
     &      surfaces, ld_surf, e_grid, a_grid,
     &      nsurf_e, nsurf_a, nsurf_pts,
     &      tobs, apobs, nobs)

c        Decode and evaluate child 2
         do k = 1, lchrom
            chr(k) = newchr(j+1, k)
         end do
         call ga_decode(chr, lchrom, xvec, nparam,
     &                  lsubstr, alow, ahigh, factor)
         do k = 1, nparam
            newx(j+1, k) = xvec(k)
         end do
         newfit(j+1) = ga_fitness(xvec,
     &      surfaces, ld_surf, e_grid, a_grid,
     &      nsurf_e, nsurf_a, nsurf_pts,
     &      tobs, apobs, nobs)

         j = j + 2
         goto 10
 20      continue

c        Statistics on new population
         call ga_stats(newfit, ipopsize,
     &                 sumfitness, avg, amax, amin, ibest)

c        Update best solution
         if (newfit(ibest) .gt. bestfit) then
            bestfit = newfit(ibest)
            igenbst = igen
            do k = 1, nparam
               bestsol(k) = newx(ibest, k)
            end do
         end if

c        Copy new → old
         do i = 1, ipopsize
            do k = 1, lchrom
               oldchr(i, k) = newchr(i, k)
            end do
            do k = 1, nparam
               oldx(i, k) = newx(i, k)
            end do
            oldfit(i) = newfit(i)
         end do

      end do

c     Output best solution
      e_opt = bestsol(1)
      a_opt = bestsol(2)
      rms_opt = (1.d0 - bestfit) / bestfit

      return
      end

c     ============================================================
      double precision function ga_fitness(xvec,
     &   surfaces, ld_surf, e_grid, a_grid,
     &   nsurf_e, nsurf_a, nsurf_pts,
     &   tobs, apobs, nobs)
c
c     Fitness = 1/(1 + RMS) where RMS is computed by
c     interpolating the 9 RSM surfaces at (e, A) and
c     comparing to observed TLE apogee.
c     ============================================================
      implicit double precision (a-h, o-z)
      integer nsurf_e, nsurf_a, nsurf_pts, ld_surf, nobs
      double precision xvec(2)
      double precision surfaces(ld_surf, nsurf_e, nsurf_a)
      double precision e_grid(nsurf_e), a_grid(nsurf_a)
      double precision tobs(nobs), apobs(nobs)

      double precision sout, fout
      double precision dum(500), tsurf(500)
      double precision f2d(3, 3)
      integer ih, ia, it

      e_val = xvec(1)
      a_val = xvec(2)

c     Build time array from first surface column
c     (surfaces are sampled at regular intervals; we interpolate
c      at tobs times)

      sout = 0.d0
      do it = 1, nobs
c        Surfaces are pre-interpolated at observation times
c        surfaces(it, ih, ia) = apogee at obs time it
         do ih = 1, nsurf_e
            do ia = 1, nsurf_a
               f2d(ih, ia) = surfaces(it, ih, ia)
            end do
         end do

c        2D interpolation at (e_val, a_val)
         call ga_twoint(e_grid, a_grid, f2d,
     &                  e_val, a_val, fout,
     &                  nsurf_e, nsurf_a)

c        Accumulate squared error
         sout = sout + ((fout - apobs(it)) / 100.d0)**2
      end do

      sout = dsqrt(sout / dble(nobs))
      ga_fitness = 1.d0 / (1.d0 + sout)

      return
      end

c     ============================================================
      subroutine ga_tinterp(surf_col, npts, t, fout)
c     Linear interpolation in time within a single surface.
c     Assumes surface is sampled at uniform dt from t=0.
c     surf_col(1:npts) = apogee values at uniform time steps.
c     ============================================================
      implicit double precision (a-h, o-z)
      integer npts
      double precision surf_col(npts), t, fout

      if (npts .le. 1) then
         fout = surf_col(1)
         return
      end if

c     Compute dt from total span (assume uniform spacing)
      dt = 1.d0 / dble(npts - 1)

c     Find bracketing index
      ridx = t / dt + 1.d0
      idx = int(ridx)
      if (idx .lt. 1) idx = 1
      if (idx .ge. npts) then
         fout = surf_col(npts)
         return
      end if

      frac = ridx - dble(idx)
      fout = surf_col(idx) + frac *
     &       (surf_col(idx+1) - surf_col(idx))

      return
      end

c     ============================================================
      subroutine ga_twoint(x, y, f, xinp, yinp, fout,
     &                     imax, jmax)
c     2D bilinear interpolation (from genpoen1.f TWOINT).
c     ============================================================
      implicit double precision (a-h, o-z)
      integer imax, jmax
      double precision x(imax), y(jmax), f(imax, jmax)
      double precision xinp, yinp, fout

      integer i, j

c     Find bracketing i-index
      do i = 2, imax
         if (xinp .le. x(i)) goto 10
      end do
      i = imax
 10   continue

c     Find bracketing j-index
      do j = 2, jmax
         if (yinp .le. y(j)) goto 30
      end do
      j = jmax
 30   continue

c     Bilinear interpolation
      ratio1 = (xinp - x(i-1)) / (x(i) - x(i-1))
      ratio2 = (yinp - y(j-1)) / (y(j) - y(j-1))
      fout1 = f(i-1, j-1) + ratio1*(f(i, j-1) - f(i-1, j-1))
      fout2 = f(i-1, j) + ratio1*(f(i, j) - f(i-1, j))
      fout = fout1 + ratio2*(fout2 - fout1)

      return
      end

c     ============================================================
c     GA internal subroutines
c     ============================================================

      subroutine ga_decode(chrom, lchrom, xvec, nparam,
     &                     lsubstr, alow, ahigh, factor)
c     Decode binary chromosome to real-valued parameters.
      implicit double precision (a-h, o-z)
      integer lchrom, nparam
      integer chrom(lchrom), lsubstr(nparam)
      double precision xvec(nparam), alow(nparam)
      double precision ahigh(nparam), factor(nparam)
      integer jpos, j, k

      jpos = 0
      do j = 1, nparam
c        Decode substring
         acc = 0.d0
         powof2 = 1.d0
         do k = 1, lsubstr(j)
            if (chrom(jpos + k) .eq. 1) acc = acc + powof2
            powof2 = powof2 * 2.d0
         end do
c        Map to [alow, ahigh]
         xvec(j) = alow(j) + (ahigh(j)-alow(j)) * acc/factor(j)
         jpos = jpos + lsubstr(j)
      end do

      return
      end

c     ============================================================
      subroutine ga_stats(fitness, ipopsize,
     &                    sumfit, avg, amax, amin, ibest)
c     Compute population statistics.
      implicit double precision (a-h, o-z)
      integer ipopsize, ibest
      double precision fitness(ipopsize)
      double precision sumfit, avg, amax, amin
      integer j

      sumfit = fitness(1)
      amin = fitness(1)
      amax = fitness(1)
      ibest = 1
      do j = 2, ipopsize
         sumfit = sumfit + fitness(j)
         if (fitness(j) .gt. amax) then
            amax = fitness(j)
            ibest = j
         end if
         if (fitness(j) .lt. amin) amin = fitness(j)
      end do
      avg = sumfit / dble(ipopsize)

      return
      end

c     ============================================================
      subroutine ga_preselect(oldfit, ipopsize, avg,
     &                        choices, fraction,
     &                        oldrand, jrand)
c     Stochastic remainder roulette wheel preselection.
      implicit double precision (a-h, o-z)
      integer ipopsize
      double precision oldfit(ipopsize), avg
      integer choices(ipopsize)
      double precision fraction(ipopsize)
      double precision oldrand(55)
      integer jrand

      integer j, k, jassign, winner

      k = 0
      do j = 1, ipopsize
         expected = oldfit(j) / avg
         jassign = int(expected)
         fraction(j) = expected - dble(jassign)
 42      if (jassign .le. 0) goto 43
         k = k + 1
         jassign = jassign - 1
         choices(k) = j
         goto 42
 43      continue
      end do

      j = 0
 44   if (k .ge. ipopsize) goto 45
      j = j + 1
      if (j .gt. ipopsize) j = 1
      if (fraction(j) .gt. 0.d0) then
         winner = ga_iflip(fraction(j), oldrand, jrand)
         if (winner .eq. 1) then
            k = k + 1
            choices(k) = j
            fraction(j) = fraction(j) - 1.d0
         end if
      end if
      goto 44
 45   continue

      return
      end

c     ============================================================
      subroutine ga_select2(choices, ipopsize,
     &                      oldrand, jrand, iselected)
c     Pick from mating pool without replacement.
      implicit double precision (a-h, o-z)
      integer ipopsize, iselected
      integer choices(ipopsize)
      double precision oldrand(55)
      integer jrand

      integer jpick

      jpick = ga_irnd(1, ipopsize, oldrand, jrand)
      iselected = choices(jpick)
      choices(jpick) = choices(ipopsize)

      return
      end

c     ============================================================
      subroutine ga_crossover(ipop, mate1, mate2,
     &   oldchr, newchr, ipopsize, lchrom,
     &   pcross, pmute, ncross, nmute,
     &   oldrand, jrand)
c     Single-point crossover with integrated mutation.
      implicit double precision (a-h, o-z)
      integer ipop, mate1, mate2, ipopsize, lchrom
      integer oldchr(ipopsize, lchrom)
      integer newchr(ipopsize, lchrom)
      double precision pcross, pmute
      integer ncross, nmute
      double precision oldrand(55)
      integer jrand

      integer j, jcross, iflipchk

      iflipchk = ga_iflip(pcross, oldrand, jrand)
      if (iflipchk .eq. 1) then
         jcross = ga_irnd(1, lchrom-1, oldrand, jrand)
         ncross = ncross + 1
      else
         jcross = lchrom
      end if

c     Copy + swap at crossover point, with mutation
      do j = 1, jcross
         newchr(ipop,   j) = ga_mutate(oldchr(mate1, j),
     &                                  pmute, nmute,
     &                                  oldrand, jrand)
         newchr(ipop+1, j) = ga_mutate(oldchr(mate2, j),
     &                                  pmute, nmute,
     &                                  oldrand, jrand)
      end do

      if (jcross .lt. lchrom) then
         do j = jcross+1, lchrom
            newchr(ipop,   j) = ga_mutate(oldchr(mate2, j),
     &                                     pmute, nmute,
     &                                     oldrand, jrand)
            newchr(ipop+1, j) = ga_mutate(oldchr(mate1, j),
     &                                     pmute, nmute,
     &                                     oldrand, jrand)
         end do
      end if

      return
      end

c     ============================================================
      integer function ga_mutate(ival, pmute, nmute,
     &                           oldrand, jrand)
c     Bit-flip mutation with probability pmute.
      implicit double precision (a-h, o-z)
      integer ival, nmute
      double precision pmute
      double precision oldrand(55)
      integer jrand

      ga_mutate = ival
      if (ga_iflip(pmute, oldrand, jrand) .eq. 1) then
         if (ival .eq. 0) then
            ga_mutate = 1
         else
            ga_mutate = 0
         end if
         nmute = nmute + 1
      end if

      return
      end

c     ============================================================
c     Random number generator (55-element Fibonacci)
c     ============================================================

      subroutine ga_warmup(seed, oldrand, jrand)
      implicit double precision (a-h, o-z)
      double precision seed, oldrand(55)
      integer jrand, j1, ii

      oldrand(55) = seed
      rand_new = 1.d-9
      prev_rand = seed
      do j1 = 1, 54
         ii = mod(21*j1, 55)
         if (ii .eq. 0) ii = 55
         oldrand(ii) = rand_new
         rand_new = prev_rand - rand_new
         if (rand_new .lt. 0.d0) rand_new = rand_new + 1.d0
         prev_rand = oldrand(ii)
      end do
      call ga_advance(oldrand)
      call ga_advance(oldrand)
      call ga_advance(oldrand)
      jrand = 0

      return
      end

c     ============================================================
      subroutine ga_advance(oldrand)
      implicit double precision (a-h, o-z)
      double precision oldrand(55), rand_new
      integer j1

      do j1 = 1, 24
         rand_new = oldrand(j1) - oldrand(j1+31)
         if (rand_new .lt. 0.d0) rand_new = rand_new + 1.d0
         oldrand(j1) = rand_new
      end do
      do j1 = 25, 55
         rand_new = oldrand(j1) - oldrand(j1-24)
         if (rand_new .lt. 0.d0) rand_new = rand_new + 1.d0
         oldrand(j1) = rand_new
      end do

      return
      end

c     ============================================================
      double precision function ga_random(oldrand, jrand)
      implicit double precision (a-h, o-z)
      double precision oldrand(55)
      integer jrand

      jrand = jrand + 1
      if (jrand .gt. 55) then
         jrand = 1
         call ga_advance(oldrand)
      end if
      ga_random = oldrand(jrand)

      return
      end

c     ============================================================
      integer function ga_iflip(prob, oldrand, jrand)
c     Biased coin flip: returns 1 with probability prob.
      implicit double precision (a-h, o-z)
      double precision prob, oldrand(55)
      integer jrand

      ga_iflip = 0
      if (ga_random(oldrand, jrand) .le. prob) ga_iflip = 1

      return
      end

c     ============================================================
      integer function ga_irnd(ilow, ihigh, oldrand, jrand)
c     Random integer in [ilow, ihigh].
      implicit double precision (a-h, o-z)
      integer ilow, ihigh
      double precision oldrand(55)
      integer jrand

      if (ilow .ge. ihigh) then
         ga_irnd = ilow
         return
      end if
      rchk = ga_random(oldrand, jrand)
      ga_irnd = int(rchk * dble(ihigh - ilow + 1)) + ilow
      if (ga_irnd .gt. ihigh) ga_irnd = ihigh

      return
      end

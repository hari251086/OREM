c     ============================================================
c     zone_select.F — Zone selection for re-entry prediction
c
c     Identifies contiguous intervals in the TLE orbital evolution
c     where the mean apogee altitude decays quasi-linearly. These
c     zones are used by RSM to generate response surfaces.
c
c     Algorithm:
c       1. Scan for decay segments (negative apogee slope)
c       2. Within each segment, test linearity (R² > threshold)
c       3. Split long segments at max_zone_days
c       4. Select best zones (highest R², prefer earlier epochs)
c
c     Author: Harishkumar Sellamuthu
c     ============================================================

c     ============================================================
      subroutine linfit(x, y, nn, slope, intercept, r2)
c     Least-squares linear fit: y = slope*x + intercept
c     Also computes R² (coefficient of determination).
c
c     Input
c       x(nn), y(nn) : data arrays
c       nn           : number of points (must be >= 2)
c     Output
c       slope, intercept : fit parameters
c       r2               : coefficient of determination [0,1]
c     ============================================================
      implicit double precision (a-h, o-z)
      integer nn
      double precision x(nn), y(nn)
      double precision slope, intercept, r2

      double precision sx, sy, sxx, sxy, syy, xmean, ymean
      double precision ss_res, ss_tot, denom
      integer i

      slope = 0.d0
      intercept = 0.d0
      r2 = 0.d0

      if (nn .lt. 2) return

      sx  = 0.d0
      sy  = 0.d0
      sxx = 0.d0
      sxy = 0.d0
      syy = 0.d0

      do i = 1, nn
         sx  = sx  + x(i)
         sy  = sy  + y(i)
         sxx = sxx + x(i)*x(i)
         sxy = sxy + x(i)*y(i)
         syy = syy + y(i)*y(i)
      end do

      denom = dble(nn)*sxx - sx*sx
      if (dabs(denom) .lt. 1.d-30) return

      slope     = (dble(nn)*sxy - sx*sy) / denom
      intercept = (sy - slope*sx) / dble(nn)

c     R² = 1 - SS_res / SS_tot
      xmean = sx / dble(nn)
      ymean = sy / dble(nn)

      ss_res = 0.d0
      ss_tot = 0.d0
      do i = 1, nn
         ss_res = ss_res + (y(i) - slope*x(i) - intercept)**2
         ss_tot = ss_tot + (y(i) - ymean)**2
      end do

      if (ss_tot .gt. 1.d-30) then
         r2 = 1.d0 - ss_res / ss_tot
      else
         r2 = 1.d0
      end if

      return
      end

c     ============================================================
      subroutine zone_select(epochs, ha, npts,
     &                       min_zone_pts, max_zone_days,
     &                       r2_threshold, slope_threshold,
     &                       nzones_max,
     &                       zone_start, zone_end, nzones)
c
c     Select zones of quasi-linear apogee decay from the TLE
c     orbital evolution time-series.
c
c     Input
c       epochs(npts)     : Julian dates (from tle_evolve)
c       ha(npts)         : mean apogee altitude in km
c       npts             : number of points
c       min_zone_pts     : minimum TLE observations per zone
c       max_zone_days    : maximum zone duration (days)
c       r2_threshold     : minimum R² for linearity (e.g. 0.90)
c       slope_threshold  : min negative slope km/day (e.g. -1.0)
c       nzones_max       : max zones to return
c
c     Output
c       zone_start(nzones_max) : start indices into epochs/ha
c       zone_end(nzones_max)   : end indices
c       nzones                 : actual number of zones found
c     ============================================================
      implicit double precision (a-h, o-z)
      integer npts, min_zone_pts, nzones_max, nzones
      double precision epochs(npts), ha(npts)
      double precision max_zone_days, r2_threshold
      double precision slope_threshold
      integer zone_start(nzones_max), zone_end(nzones_max)

c     Candidate zones (internal buffer)
      parameter(max_cand = 100)
      integer cand_s(max_cand), cand_e(max_cand)
      double precision cand_r2(max_cand), cand_sl(max_cand)
      integer ncand

c     Working variables
      double precision slope, intercept, r2
      double precision span_days
      integer i, j, k, istart, iend, best_idx
      double precision best_r2

      nzones = 0
      ncand = 0

      if (npts .lt. min_zone_pts) return

c     ---------------------------------------------------------
c     Pass 1: Find contiguous decay segments
c     A segment starts where ha(i+1) < ha(i) and continues
c     as long as the overall trend is downward.
c     We use a sliding window approach: for each possible
c     starting point, extend forward while R² stays high.
c     ---------------------------------------------------------
      i = 1
 10   continue
      if (i .gt. npts - min_zone_pts + 1) goto 50

c     Try starting a zone at index i
c     Extend forward to find the longest linear segment
      istart = i
      iend = i + min_zone_pts - 1

c     Check if minimum window fits within max_zone_days
      if (epochs(iend) - epochs(istart) .gt.
     &    max_zone_days) then
         i = i + 1
         goto 10
      end if

c     Check if this minimum window has negative slope
      call linfit(epochs(istart), ha(istart),
     &            iend - istart + 1, slope, intercept, r2)

      if (slope .ge. slope_threshold) then
c        Not decaying enough — skip forward
         i = i + 1
         goto 10
      end if

c     Extend the window as long as linearity holds
 20   continue
      if (iend .ge. npts) goto 30
      if (epochs(iend+1) - epochs(istart) .gt.
     &    max_zone_days) goto 30

      call linfit(epochs(istart), ha(istart),
     &            iend - istart + 2, slope, intercept, r2)

      if (r2 .ge. r2_threshold .and.
     &    slope .lt. slope_threshold) then
         iend = iend + 1
         goto 20
      end if

 30   continue

c     Verify the final window
      call linfit(epochs(istart), ha(istart),
     &            iend - istart + 1, slope, intercept, r2)

      if (r2 .ge. r2_threshold .and.
     &    slope .lt. slope_threshold .and.
     &    iend - istart + 1 .ge. min_zone_pts) then

c        Store as candidate zone
         if (ncand .lt. max_cand) then
            ncand = ncand + 1
            cand_s(ncand) = istart
            cand_e(ncand) = iend
            cand_r2(ncand) = r2
            cand_sl(ncand) = slope
         end if

c        Jump past this zone to find next one
         i = iend + 1
      else
         i = i + 1
      end if

      goto 10

 50   continue

c     ---------------------------------------------------------
c     Pass 2: Select best zones from candidates
c     Prefer higher R², earlier epochs. No overlaps.
c     ---------------------------------------------------------
      do k = 1, nzones_max
         if (ncand .eq. 0) goto 60

c        Find candidate with best R²
         best_idx = 0
         best_r2 = -1.d0
         do j = 1, ncand
            if (cand_r2(j) .gt. best_r2) then
               best_r2 = cand_r2(j)
               best_idx = j
            end if
         end do

         if (best_idx .eq. 0) goto 60

c        Add to output
         nzones = nzones + 1
         zone_start(nzones) = cand_s(best_idx)
         zone_end(nzones)   = cand_e(best_idx)

c        Remove this candidate and any overlapping ones
         istart = cand_s(best_idx)
         iend   = cand_e(best_idx)

c        Mark used candidate by setting R² to -1
         cand_r2(best_idx) = -1.d0

c        Remove overlapping candidates
         do j = 1, ncand
            if (cand_r2(j) .gt. 0.d0) then
               if (cand_s(j) .le. iend .and.
     &             cand_e(j) .ge. istart) then
                  cand_r2(j) = -1.d0
               end if
            end if
         end do
      end do

 60   continue

c     ---------------------------------------------------------
c     Pass 3: Sort zones by epoch (ascending)
c     ---------------------------------------------------------
      do i = 1, nzones - 1
         do j = i + 1, nzones
            if (zone_start(j) .lt. zone_start(i)) then
               k = zone_start(i)
               zone_start(i) = zone_start(j)
               zone_start(j) = k
               k = zone_end(i)
               zone_end(i) = zone_end(j)
               zone_end(j) = k
            end if
         end do
      end do

      return
      end

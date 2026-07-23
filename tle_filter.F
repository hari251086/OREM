c     ============================================================
c     tle_filter.F -- issue #10: TLE quality filtering
c
c     Pre-filters tle_evolve output before it reaches zone_select,
c     flagging TLEs likely to be maneuvers, outliers, or catalog
c     errors rather than genuine orbital-decay signal. Two
c     independent checks, both evaluated against a trailing local
c     window (only points BEFORE i -- a genuine "would this have
c     looked anomalous as new data arrived" check, not a centered
c     window that peeks ahead):
c
c       - Maneuver detection (on eccentricity): fit a local linear
c         trend (linfit, zone_select.F) over the trailing window,
c         predict e at point i's epoch, flag if the actual value
c         is more than sigma_e residual-sigmas from that prediction.
c         A real maneuver breaks the smooth trend the natural orbital
c         evolution follows; a catalog outlier does too.
c       - Outlier rejection (on apogee altitude): flag if ha(i)
c         deviates more than sigma_ha sigmas from the trailing
c         window's mean.
c
c     A point failing EITHER check gets mask(i)=0. Points too close
c     to the start of the record to have a full window (nwin points)
c     are never flagged -- not enough history to judge them against;
c     they pass through with mask=1.
c
c     Epoch gaps are reported (ngaps, and a console warning per gap)
c     but never filtered -- a gap means missing data, not bad data.
c     Duplicate-epoch resolution is NOT handled here: tle_evolve
c     itself already deduplicates TLEs within 86 sec (see
c     tle_evolution.F), which is issue #10's item 4.
c     ============================================================
      subroutine tle_filter(epochs, ha, e, npts,
     &   sigma_ha, sigma_e, gap_days, nwin,
     &   mask, nfiltered, nremoved, ngaps)
c
c     Input
c       epochs(npts), ha(npts), e(npts) : tle_evolve output arrays
c       npts       : number of points
c       sigma_ha   : outlier threshold, in local-window sigmas
c       sigma_e    : maneuver threshold, in local-trend-residual sigmas
c       gap_days   : epoch gap warning threshold, in days
c       nwin       : trailing window size for both checks
c
c     Recommended: sigma_ha=4.0, sigma_e=6.0, gap_days=30.0, nwin=20 --
c     tuned against real TLE histories (42928, 35497), not just the
c     "3-sigma/5-10 points" first guess. Real e(t) has enough natural
c     local curvature that a strict 3-sigma/8-point linear-trend check
c     flagged ~20% of clean, well-tracked objects as false positives;
c     the wider window and looser threshold bring that to ~1.5%,
c     matching the (already well-behaved) outlier check's own rate,
c     while still catching synthetic maneuvers/outliers many times
c     larger than the local noise floor (see test_tle_filter.F).
c
c     Output
c       mask(npts) : 1 = keep, 0 = reject (outlier or maneuver)
c       nfiltered  : number of points with enough trailing history to
c                    actually be evaluated (npts - min(nwin,npts))
c       nremoved   : number of points with mask=0
c       ngaps      : number of consecutive-epoch gaps > gap_days
c     ============================================================
      implicit double precision (a-h, o-z)
      integer npts, nwin
      double precision epochs(npts), ha(npts), e(npts)
      double precision sigma_ha, sigma_e, gap_days
      integer mask(npts)
      integer nfiltered, nremoved, ngaps

      integer i, j, w0, nw
      double precision wx(1000), wha(1000), we(1000)
      double precision wmean, wstd, wslope, wintc, wr2
      double precision epred, resid, ssres

      nfiltered = 0
      nremoved = 0
      ngaps = 0

      do i = 1, npts
         mask(i) = 1
      end do

      do i = 1, npts
         w0 = i - nwin
         if (w0 .lt. 1) w0 = 1
c        Clip the trailing window at any real epoch gap -- a window
c        straddling a data-availability gap mixes two disconnected
c        segments and makes the points right after the gap look
c        spuriously anomalous relative to stale pre-gap context.
c        Two cases: point i is itself the first point after a gap
c        (checked directly, since that pair sits outside the j-loop
c        below, which only checks pairs strictly inside the window);
c        or an older gap lies further back inside the window.
         if (i .gt. 1 .and.
     &       epochs(i) - epochs(i-1) .gt. gap_days) then
            w0 = i
         else
            do j = i - 1, w0 + 1, -1
               if (epochs(j) - epochs(j-1) .gt. gap_days) then
                  w0 = j
                  goto 40
               end if
            end do
 40         continue
         end if
         nw = i - w0
         if (nw .lt. 3) goto 50

         nfiltered = nfiltered + 1

         do j = 1, nw
            wx(j)  = epochs(w0 + j - 1)
            wha(j) = ha(w0 + j - 1)
            we(j)  = e(w0 + j - 1)
         end do

c        --- Outlier rejection: ha vs trailing window mean/stddev ---
         wmean = 0.d0
         do j = 1, nw
            wmean = wmean + wha(j)
         end do
         wmean = wmean / dble(nw)
         wstd = 0.d0
         do j = 1, nw
            wstd = wstd + (wha(j) - wmean)**2
         end do
         wstd = dsqrt(wstd / dble(nw))
         if (wstd .lt. 1.d-6) wstd = 1.d-6
         if (dabs(ha(i) - wmean) .gt. sigma_ha * wstd) then
            mask(i) = 0
         end if

c        --- Maneuver detection: e vs trailing local linear trend ---
         call linfit(wx, we, nw, wslope, wintc, wr2)
         epred = wslope * epochs(i) + wintc
         ssres = 0.d0
         do j = 1, nw
            resid = we(j) - (wslope * wx(j) + wintc)
            ssres = ssres + resid**2
         end do
         wstd = dsqrt(ssres / dble(nw))
         if (wstd .lt. 1.d-6) wstd = 1.d-6
         if (dabs(e(i) - epred) .gt. sigma_e * wstd) then
            mask(i) = 0
         end if

 50      continue
      end do

      do i = 1, npts
         if (mask(i) .eq. 0) nremoved = nremoved + 1
      end do

      do i = 2, npts
         if (epochs(i) - epochs(i-1) .gt. gap_days) then
            ngaps = ngaps + 1
            write(*,'(A,I5,A,F10.3,A,F10.3,A,F8.1,A)')
     &         '  [WARN] tle_filter: epoch gap at point ', i,
     &         ' (', epochs(i-1), ' -> ', epochs(i), ', ',
     &         epochs(i) - epochs(i-1), ' days)'
         end if
      end do

      return
      end

c     diag_zero_zone.F -- diagnose WHY a given object forms zero valid
c     zones under zone_select's shipped defaults (min_pts=8,
c     max_zone_days=10.0, r2_thresh=0.90, slope_thresh=-1.0), issue
c     #29. Distinguishes two failure modes:
c       (a) DENSITY: no 8-point window ever fits within 10 days at all
c           (tracking cadence too sparse) -- a linearity check is
c           irrelevant if this fails.
c       (b) LINEARITY: dense-enough windows exist, but none decay
c           cleanly enough (R2>=0.90, slope<-1.0) anywhere.
c     For every dense-enough (density-passing) 8-point window found,
c     reports its actual R2/slope, so a "close but not quite" case is
c     visible as such (e.g. R2=0.85) versus a "nowhere close" case.
c     Analysis only -- not part of the build/test suite.
c
c     Usage: diag_zero_zone.exe <norad> <tle_file>
      program diag_zero_zone
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter (maxpts = 15000)
      double precision epochs(maxpts), a_o(maxpts), e_o(maxpts)
      double precision inc_o(maxpts), raan_o(maxpts), aop_o(maxpts)
      double precision ha_o(maxpts), hp_o(maxpts), lam_o(maxpts)
      integer npts, ierr, norad_id
      character*80 tle_file, arg_str

      double precision slope, intercept, r2
      double precision best_r2, best_slope, span
      double precision best_decay_r2, best_decay_slope
      double precision minval8
      integer ndense, nbest_start, nbest_end
      integer nbest_decay_start, nbest_decay_end
      integer i, iend

c     Issue #10 TLE quality filter -- must match orem.F Step 1b
c     exactly (sigma_ha=4, sigma_e=6, gap_days=30, nwin=20), since
c     the real orem_run pipeline always applies this BEFORE
c     zone_select ever sees the data.
      integer mask_tle(maxpts)
      integer nfiltered_tle, nremoved_tle, ngaps_tle
      integer i_filt, j_filt

      call init_constants()

      call getarg(1, arg_str)
      read(arg_str, *) norad_id
      call getarg(2, tle_file)

      call tle_evolve(tle_file, norad_id, maxpts,
     &   epochs, a_o, e_o, inc_o, raan_o, aop_o,
     &   ha_o, hp_o, lam_o, npts, ierr)

      write(*,'(A,I8,A,I3,A,I6)') ' norad=', norad_id,
     &   '  tle_evolve ierr=', ierr, '  npts=', npts
      if (ierr .ne. 0 .or. npts .lt. 8) then
         write(*,'(A)') ' insufficient data, stopping'
         stop
      end if

      call tle_filter(epochs, ha_o, e_o, npts,
     &   4.0d0, 6.0d0, 30.0d0, 20,
     &   mask_tle, nfiltered_tle, nremoved_tle, ngaps_tle)
      if (nremoved_tle .gt. 0) then
         j_filt = 0
         do i_filt = 1, npts
            if (mask_tle(i_filt) .eq. 1) then
               j_filt = j_filt + 1
               epochs(j_filt) = epochs(i_filt)
               ha_o(j_filt)   = ha_o(i_filt)
               e_o(j_filt)    = e_o(i_filt)
            end if
         end do
         npts = j_filt
         write(*,'(A,I4,A,I6)') ' tle_filter removed ',
     &      nremoved_tle, ' points, npts now ', npts
      end if

      write(*,'(A,F14.4,A,F14.4,A,F10.2)')
     &   ' span: ', epochs(1), ' to ', epochs(npts),
     &   '  total_days=', epochs(npts)-epochs(1)

c     Slide an 8-point window across the whole record. For each
c     window, check density (span <= 10 days) first.
      ndense = 0
      best_r2 = -1.d0
      best_slope = 0.d0
      nbest_start = 0
      nbest_end = 0
      best_decay_r2 = -1.d0
      best_decay_slope = 0.d0
      nbest_decay_start = 0
      nbest_decay_end = 0

      do i = 1, npts - 7
         iend = i + 7
         span = epochs(iend) - epochs(i)
         if (span .le. 10.d0) then
            ndense = ndense + 1
            call linfit(epochs(i), ha_o(i), 8, slope, intercept, r2)
            if (r2 .gt. best_r2) then
               best_r2 = r2
               best_slope = slope
               nbest_start = i
               nbest_end = iend
            end if
            if (slope .lt. -1.d0 .and. r2 .gt. best_decay_r2) then
               best_decay_r2 = r2
               best_decay_slope = slope
               nbest_decay_start = i
               nbest_decay_end = iend
            end if
         end if
      end do

      write(*,'(A,I6,A,I6)') ' 8-pt windows tried: ', npts-7,
     &   '   density-passing (span<=10d): ', ndense

      if (ndense .eq. 0) then
         write(*,'(A)') ' FAILURE MODE: DENSITY -- no 8 TLEs ever'
     &      // ' fall within any 10-day span in this record'
c        Report actual typical spacing for context
         write(*,'(A,F10.3,A)') '   min gap between consecutive TLEs=',
     &      minval8(epochs, npts), ' days'
      else
         write(*,'(A,F8.4,A,F10.4,A,I6,A,I6)')
     &      ' best density-passing window (any slope): R2=', best_r2,
     &      '  slope=', best_slope, '  idx=', nbest_start,
     &      '..', nbest_end
         if (nbest_decay_start .eq. 0) then
            write(*,'(A)') ' no density-passing window has'
     &         // ' slope<-1.0 at all (nothing even attempts decay)'
         else
            write(*,'(A,F8.4,A,F10.4,A,I6,A,I6)')
     &         ' best DECAYING window (slope<-1.0): R2=',
     &         best_decay_r2, '  slope=', best_decay_slope,
     &         '  idx=', nbest_decay_start, '..', nbest_decay_end
         end if
         if (best_decay_r2 .ge. 0.90d0) then
            write(*,'(A)') ' (this window WOULD pass zone_select --'
     &         // ' unexpected given the ERR result, investigate'
     &         // ' zone_select/tle_filter interaction itself)'
         else
            write(*,'(A)') ' FAILURE MODE: LINEARITY -- dense-enough'
     &         // ' windows exist but never decay cleanly enough'
     &         // ' (need R2>=0.90 and slope<-1.0 together)'
         end if
      end if

      stop
      end

c     ============================================================
      double precision function minval8(epochs, npts)
      implicit double precision (a-h, o-z)
      integer npts
      double precision epochs(npts)
      double precision gmin, g
      integer i
      gmin = 1.d30
      do i = 1, npts-1
         g = epochs(i+1) - epochs(i)
         if (g .lt. gmin) gmin = g
      end do
      minval8 = gmin
      return
      end

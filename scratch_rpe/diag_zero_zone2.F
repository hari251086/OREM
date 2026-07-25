c     diag_zero_zone2.F -- issue #29 follow-up: for each of the 12
c     objects with zero valid zones at the shipped defaults
c     (min_pts=8, max_zone_days=10), sweep max_zone_days to see how
c     much widening would be needed to admit an 8-point window at
c     all, and whether that window (once dense enough) also clears
c     the linearity bar (R2>=0.90, slope<-1.0). Applies tle_filter
c     first, matching the real orem_run pipeline. Analysis only.
c
c     Usage: diag_zero_zone2.exe <norad> <tle_file>
      program diag_zero_zone2
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter (maxpts = 15000, nspan = 6)
      double precision epochs(maxpts), a_o(maxpts), e_o(maxpts)
      double precision inc_o(maxpts), raan_o(maxpts), aop_o(maxpts)
      double precision ha_o(maxpts), hp_o(maxpts), lam_o(maxpts)
      integer npts, ierr, norad_id
      character*80 tle_file, arg_str

      double precision slope, intercept, r2, span
      double precision spans(nspan)
      double precision best_decay_r2(nspan), best_decay_slope(nspan)
      integer ndense(nspan)
      integer i, iend, k

      integer mask_tle(maxpts)
      integer nfiltered_tle, nremoved_tle, ngaps_tle
      integer i_filt, j_filt

      data spans /10.d0, 15.d0, 20.d0, 30.d0, 45.d0, 60.d0/

      call init_constants()
      call getarg(1, arg_str)
      read(arg_str, *) norad_id
      call getarg(2, tle_file)

      call tle_evolve(tle_file, norad_id, maxpts,
     &   epochs, a_o, e_o, inc_o, raan_o, aop_o,
     &   ha_o, hp_o, lam_o, npts, ierr)
      if (ierr .ne. 0 .or. npts .lt. 8) then
         write(*,'(A,I8,A)') ' norad=', norad_id,
     &      ' insufficient data'
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
            end if
         end do
         npts = j_filt
      end if

      write(*,'(A,I8,A,I6)') ' norad=', norad_id,
     &   '  npts(post-filter)=', npts

      do k = 1, nspan
         ndense(k) = 0
         best_decay_r2(k) = -1.d0
         best_decay_slope(k) = 0.d0
      end do

      do i = 1, npts - 7
         iend = i + 7
         span = epochs(iend) - epochs(i)
         do k = 1, nspan
            if (span .le. spans(k)) then
               ndense(k) = ndense(k) + 1
               call linfit(epochs(i), ha_o(i), 8, slope,
     &            intercept, r2)
               if (slope .lt. -1.d0 .and.
     &             r2 .gt. best_decay_r2(k)) then
                  best_decay_r2(k) = r2
                  best_decay_slope(k) = slope
               end if
            end if
         end do
      end do

      write(*,'(A)') '  span  ndense  best_decay_r2  best_decay_slope'
      do k = 1, nspan
         write(*,'(F6.1,I8,F13.4,F16.4)') spans(k), ndense(k),
     &      best_decay_r2(k), best_decay_slope(k)
      end do

      stop
      end

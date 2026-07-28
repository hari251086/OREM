c     diag_zone_data.F -- issue #29 generic follow-up to
c     diag_61734_zone1_data.F: same inspection (zone_select's actual
c     chosen windows, their R2/slope/span, and per-point raw data for
c     a chosen zone), generalized to take <norad> <tle_file> <zone_idx>
c     as arguments instead of being hardcoded to one object. Reuses
c     zone_select's exact production settings (min_pts=8,
c     max_zone_days=20, r2=0.90, slope=-1.0, nzones_max=8).
c
c     Usage: diag_zone_data.exe <norad> <tle_file> <zone_idx>
c
c     Compile:
c       ifx /heap-arrays /F:16777216 scratch_rpe\diag_zone_data.F
c           src\tle_evolution.F src\tle_filter.F src\zone_select.F
c           build\dependencies\ksrop\src\Subrouts.F
c           build\dependencies\ksrop\src\Legendre.F
c           build\dependencies\ksrop\src\TLEread.F
c           /exe:scratch_rpe\diag_zone_data.exe
c
c     Analysis only -- not part of the build/test suite.
      program diag_zone_data
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter (maxpts = 15000, mxz = 8)
      double precision epochs(maxpts), a_o(maxpts), e_o(maxpts)
      double precision inc_o(maxpts), raan_o(maxpts), aop_o(maxpts)
      double precision ha_o(maxpts), hp_o(maxpts), lam_o(maxpts)
      integer npts, ierr, norad_id, zshow

      integer mask_tle(maxpts)
      integer nfiltered_tle, nremoved_tle, ngaps_tle
      integer i_filt, j_filt

      integer zone_start(mxz), zone_end(mxz), nzones
      double precision slope, intercept, r2
      integer iz, i, n1
      character*80 tle_file, arg_str

      call init_constants()

      call getarg(1, arg_str)
      read(arg_str, *) norad_id
      call getarg(2, tle_file)
      call getarg(3, arg_str)
      if (len_trim(arg_str) .eq. 0) then
         zshow = 1
      else
         read(arg_str, *) zshow
      end if

      call tle_evolve(tle_file, norad_id, maxpts,
     &   epochs, a_o, e_o, inc_o, raan_o, aop_o,
     &   ha_o, hp_o, lam_o, npts, ierr)
      write(*,'(A,I3,A,I6)') ' tle_evolve ierr=', ierr, '  npts=', npts
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
               a_o(j_filt)    = a_o(i_filt)
               inc_o(j_filt)  = inc_o(i_filt)
               raan_o(j_filt) = raan_o(i_filt)
               aop_o(j_filt)  = aop_o(i_filt)
            end if
         end do
         npts = j_filt
      end if
      write(*,'(A,I4,A,I6)') ' tle_filter removed ', nremoved_tle,
     &   ' points, npts now ', npts

      call zone_select(epochs, ha_o, npts,
     &   8, 20.d0, 0.90d0, -1.d0,
     &   mxz, zone_start, zone_end, nzones)
      write(*,'(A,I3)') ' nzones=', nzones

      do iz = 1, nzones
         n1 = zone_end(iz) - zone_start(iz) + 1
         call linfit(epochs(zone_start(iz)), ha_o(zone_start(iz)),
     &      n1, slope, intercept, r2)
         write(*,'(A,I2,A,I4,A,I4,A,I3,A,F9.3,A,F9.5,A,F9.4)')
     &      ' Z', iz, ': idx ', zone_start(iz), '..', zone_end(iz),
     &      '  npts=', n1,
     &      '  span_days=', epochs(zone_end(iz))-epochs(zone_start(iz)),
     &      '  R2=', r2, '  slope=', slope
      end do

      if (zshow .ge. 1 .and. zshow .le. nzones) then
         write(*,'(/,A,I2,A)') ' --- Zone ', zshow, ' raw '//
     &      '(epoch_day_since_zone_start, ha_km, e, a_km) ---'
         n1 = zone_end(zshow) - zone_start(zshow) + 1
         do i = zone_start(zshow), zone_end(zshow)
            write(*,'(A,F9.4,A,F10.3,A,F9.6,A,F12.3)')
     &         '   t=', epochs(i)-epochs(zone_start(zshow)),
     &         '  ha=', ha_o(i), '  e=', e_o(i), '  a=', a_o(i)
         end do
      end if

      write(*,'(/,A,F14.4,A,F14.4,A,F10.2)')
     &   ' full record span: ', epochs(1), ' to ', epochs(npts),
     &   '  total_days=', epochs(npts)-epochs(1)

      stop
      end

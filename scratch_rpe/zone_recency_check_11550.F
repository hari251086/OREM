c     zone_recency_check_11550.F -- before/after check for the Pass 1.5
c     recency-guarantee fix in zone_select.F. Recreates the ORIGINAL
c     crowding scenario (nzones_max=10, the pre-mxz=50-workaround cap)
c     against 11550's real 46-year TLE history, and reports how close
c     the LAST selected zone's end epoch is to the last real TLE epoch
c     in the whole dataset -- the terminal-decay zone was previously
c     crowded out entirely at small nzones_max (main_orem.F's mxz
c     comment). Analysis only -- not part of the build/test suite.
c
c     Compile (fpm-relative paths):
c       ifx /heap-arrays scratch_rpe\zone_recency_check_11550.F
c           src\tle_evolution.F src\zone_select.F
c           build\dependencies\ksrop\src\TLEread.F
c           build\dependencies\ksrop\src\Subrouts.F
c           build\dependencies\ksrop\src\Legendre.F
c           /exe:scratch_rpe\zone_recency_check_11550.exe
      program zone_recency_check_11550
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter (maxpts = 15000, nzmax = 10)
      double precision epochs(maxpts), a_o(maxpts), e_o(maxpts)
      double precision inc_o(maxpts), raan_o(maxpts), aop_o(maxpts)
      double precision ha_o(maxpts), hp_o(maxpts), lam_o(maxpts)
      integer npts, ierr
      integer zone_start(nzmax), zone_end(nzmax), nzones
      double precision last_epoch, gap_days

      call init_constants()

      call tle_evolve('input/example_11550.tle.txt',
     &   11550, maxpts,
     &   epochs, a_o, e_o, inc_o, raan_o, aop_o,
     &   ha_o, hp_o, lam_o, npts, ierr)

      write(*,'(A,I3,A,I6)') ' tle_evolve ierr=', ierr,
     &   '  npts=', npts
      last_epoch = epochs(npts)
      write(*,'(A,F14.4,A,F10.2)') ' last TLE epoch (JD)=',
     &   last_epoch, '  hp(last) km=', hp_o(npts)

      write(*,'(A,I3,A)') ' --- zone_select at nzones_max=', nzmax,
     &   ' (original pre-workaround cap) ---'
      call zone_select(epochs, ha_o, npts,
     &   8, 10.0d0, 0.90d0, -1.0d0,
     &   nzmax, zone_start, zone_end, nzones)

      write(*,'(A,I5)') ' nzones selected=', nzones
      write(*,'(A)') '  Z  npts   start_JD        end_JD' //
     &   '       gap_to_last_TLE(days)'
      do iz = 1, nzones
         n = zone_end(iz) - zone_start(iz) + 1
         gap_days = last_epoch - epochs(zone_end(iz))
         write(*,'(I4,I6,F14.4,F14.4,F14.2)') iz, n,
     &      epochs(zone_start(iz)), epochs(zone_end(iz)), gap_days
      end do

      gap_days = last_epoch - epochs(zone_end(nzones))
      write(*,'(A,F10.2,A)') ' Most-recent-selected-zone gap to '
     &   // 'last real TLE =', gap_days, ' days'

      stop
      end

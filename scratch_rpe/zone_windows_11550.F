c     zone_windows_11550.F -- dump each zone's actual TLE epoch window
c     (start, end, duration in days, npts) for the 11550 case, at the
c     default cfg params (min_pts=8, max_days=10.0, R2=0.90) across
c     its full 1979-2025 (46-year) TLE history. Cheap zone-select-only
c     pass before committing to a full RSM/GA/propagate_ks run.
c     Analysis only -- not part of the build/test suite.
c
c     Compile:
c       ifx scratch_rpe\zone_windows_11550.F tle_evolution.F
c           zone_select.F ksrop\TLEread.F ksrop\Subrouts.F
c           ksrop\Legendre.F /exe:scratch_rpe\zone_windows_11550.exe
      program zone_windows_11550
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter (maxpts = 15000, nzones_max = 2000)
      double precision epochs(maxpts), a_o(maxpts), e_o(maxpts)
      double precision inc_o(maxpts), raan_o(maxpts), aop_o(maxpts)
      double precision ha_o(maxpts), hp_o(maxpts), lam_o(maxpts)
      integer npts, ierr
      integer zone_start(nzones_max), zone_end(nzones_max), nzones

      call init_constants()

      call tle_evolve('input/example_40943.tle.txt',
     &   40943, maxpts,
     &   epochs, a_o, e_o, inc_o, raan_o, aop_o,
     &   ha_o, hp_o, lam_o, npts, ierr)

      write(*,'(A,I3,A,I6)') ' tle_evolve ierr=', ierr,
     &   '  npts=', npts

      write(*,'(A)') ' --- full candidate pool (nzones_max=2000) ---'
      call zone_select(epochs, ha_o, npts,
     &   8, 10.0d0, 0.90d0, -1.0d0,
     &   nzones_max, zone_start, zone_end, nzones)

      write(*,'(A,I5)') ' nzones=', nzones
      write(*,'(A)') '  Z  npts   start_JD        end_JD' //
     &   '       dur(days)'
      do iz = 1, nzones
         n = zone_end(iz) - zone_start(iz) + 1
         dur = epochs(zone_end(iz)) - epochs(zone_start(iz))
         write(*,'(I4,I6,F14.4,F14.4,F10.3)') iz, n,
     &      epochs(zone_start(iz)), epochs(zone_end(iz)), dur
      end do

      write(*,'(A)') ' --- what orem.F actually gets: nzones_max=50'
     &   // ' (mxz cap) ---'
      call zone_select(epochs, ha_o, npts,
     &   8, 10.0d0, 0.90d0, -1.0d0,
     &   50, zone_start, zone_end, nzones)
      write(*,'(A,I5)') ' nzones=', nzones
      do iz = 1, nzones
         n = zone_end(iz) - zone_start(iz) + 1
         dur = epochs(zone_end(iz)) - epochs(zone_start(iz))
         write(*,'(I4,I6,F14.4,F14.4,F10.3)') iz, n,
     &      epochs(zone_start(iz)), epochs(zone_end(iz)), dur
      end do

      write(*,'(A,F14.4,A,F14.4)') ' last TLE epoch=',
     &   epochs(npts), '  hp(last)=', hp_o(npts)
      write(*,'(A,F10.2,A,F10.2)') ' ha(last)=', ha_o(npts),
     &   '  e(last)=', e_o(npts)

      stop
      end

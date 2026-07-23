c     zone_windows_33587.F -- dump each zone's actual TLE epoch window
c     (start, end, duration in days, npts) for the 33587 case.
c     nzones_max raised to 15 (from cfg's 8) and min_zone_pts swept
c     8,6,5,4,3 to check whether the cfg's 8-point/10-day requirement
c     -- not nzones_max -- is why only 5 zones are ever found (user
c     request: "take more zones past zone 5").
c     Analysis only -- not part of the build/test suite.
c
c     Compile:
c       ifx scratch_rpe\zone_windows_33587.F tle_evolution.F
c           zone_select.F ksrop\TLEread.F ksrop\Subrouts.F
c           ksrop\Legendre.F /exe:scratch_rpe\zone_windows_33587.exe
      program zone_windows_33587
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter (maxpts = 400, nzones_max = 15)
      double precision epochs(maxpts), a_o(maxpts), e_o(maxpts)
      double precision inc_o(maxpts), raan_o(maxpts), aop_o(maxpts)
      double precision ha_o(maxpts), hp_o(maxpts), lam_o(maxpts)
      integer npts, ierr
      integer zone_start(nzones_max), zone_end(nzones_max), nzones
      integer minpts_try(5), itry

      call init_constants()

      call tle_evolve('input/example_33587.tle.txt', 33587, maxpts,
     &   epochs, a_o, e_o, inc_o, raan_o, aop_o,
     &   ha_o, hp_o, lam_o, npts, ierr)

      write(*,'(A,I3,A,I5)') ' tle_evolve ierr=', ierr,
     &   '  npts=', npts

      minpts_try(1) = 8
      minpts_try(2) = 6
      minpts_try(3) = 5
      minpts_try(4) = 4
      minpts_try(5) = 3

      do itry = 1, 5
         write(*,'(A,I2,A)') ' === min_zone_pts=',
     &      minpts_try(itry), ', max_days=10.0, R2=0.90 ==='

         call zone_select(epochs, ha_o, npts,
     &      minpts_try(itry), 10.0d0, 0.90d0, -1.0d0,
     &      nzones_max, zone_start, zone_end, nzones)

         write(*,'(A,I2)') ' nzones=', nzones
         write(*,'(A)') '  Z  npts   start_JD        end_JD' //
     &      '       dur(days)'
         do iz = 1, nzones
            n = zone_end(iz) - zone_start(iz) + 1
            dur = epochs(zone_end(iz)) - epochs(zone_start(iz))
            write(*,'(I3,I6,F14.4,F14.4,F10.3)') iz, n,
     &         epochs(zone_start(iz)), epochs(zone_end(iz)), dur
         end do
      end do

      stop
      end

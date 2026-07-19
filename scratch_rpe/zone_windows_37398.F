c     zone_windows_37398.F -- OREM-Watchlist issue #9: does 37398's
c     candidate zone pool contain anything recent (2024-2026), or
c     are the top-R2 candidates genuinely all a decade old regardless
c     of nzones_max? Sweeps nzones_max at the shipped default params
c     (min_pts=8, max_days=10, R2=0.90).
c
c     Compile:
c       ifx scratch_rpe\zone_windows_37398.F tle_evolution.F
c           zone_select.F ksrop\TLEread.F ksrop\Subrouts.F
c           ksrop\Legendre.F /exe:scratch_rpe\zone_windows_37398.exe
      program zone_windows_37398
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter (maxpts = 15000, nzones_max = 2000)
      double precision epochs(maxpts), a_o(maxpts), e_o(maxpts)
      double precision inc_o(maxpts), raan_o(maxpts), aop_o(maxpts)
      double precision ha_o(maxpts), hp_o(maxpts), lam_o(maxpts)
      integer npts, ierr
      integer zone_start(nzones_max), zone_end(nzones_max), nzones
      integer try(6), itry, iz, n
      double precision dur, cal(6), r2t

      call init_constants()

      call tle_evolve(
     &   'E:\Research\1. R&D\Re-entry\2026\Data\objects\37398'//
     &   '\37398.tle.txt',
     &   37398, maxpts,
     &   epochs, a_o, e_o, inc_o, raan_o, aop_o,
     &   ha_o, hp_o, lam_o, npts, ierr)

      write(*,'(A,I3,A,I6)') ' tle_evolve ierr=', ierr,
     &   '  npts=', npts
      call jd2cal(epochs(1), cal)
      write(*,'(A,I5,I3,I3)') ' first epoch: ',
     &   int(cal(1)), int(cal(2)), int(cal(3))
      call jd2cal(epochs(npts), cal)
      write(*,'(A,I5,I3,I3)') ' last epoch:  ',
     &   int(cal(1)), int(cal(2)), int(cal(3))

      try(1) = 8
      try(2) = 20
      try(3) = 50
      try(4) = 200
      try(5) = 500
      try(6) = 2000

      do itry = 1, 6
         call zone_select(epochs, ha_o, npts,
     &      8, 10.0d0, 0.90d0, -1.0d0,
     &      try(itry), zone_start, zone_end, nzones)
         write(*,'(A,I5,A,I5)') ' --- nzones_max=', try(itry),
     &      '  -> nzones=', nzones
c        Report the LATEST (most recent-epoch) zone found, and
c        whether ANY zone starts in 2023 or later.
         iz = 0
         do n = 1, nzones
            if (iz .eq. 0) then
               iz = n
            else if (epochs(zone_start(n)) .gt.
     &               epochs(zone_start(iz))) then
               iz = n
            end if
         end do
         if (iz .gt. 0) then
            call jd2cal(epochs(zone_start(iz)), cal)
            write(*,'(A,I5,I3,I3,A,I5)') '   latest zone starts ',
     &         int(cal(1)), int(cal(2)), int(cal(3)),
     &         '  (zone #', iz
         end if
      end do

      write(*,*) ''
      write(*,*) '=== R2 threshold sweep, nzones_max=2000 ==='
      do itry = 1, 5
         if (itry .eq. 1) r2t = 0.90d0
         if (itry .eq. 2) r2t = 0.70d0
         if (itry .eq. 3) r2t = 0.50d0
         if (itry .eq. 4) r2t = 0.30d0
         if (itry .eq. 5) r2t = 0.0d0
         call zone_select(epochs, ha_o, npts,
     &      8, 10.0d0, r2t, -1.0d0,
     &      nzones_max, zone_start, zone_end, nzones)
         iz = 0
         do n = 1, nzones
            if (iz .eq. 0) then
               iz = n
            else if (epochs(zone_start(n)) .gt.
     &               epochs(zone_start(iz))) then
               iz = n
            end if
         end do
         write(*,'(A,F5.2,A,I5)') ' R2>=', r2t, '  nzones=', nzones
         if (iz .gt. 0) then
            call jd2cal(epochs(zone_start(iz)), cal)
            write(*,'(A,I5,I3,I3)') '   latest zone starts ',
     &         int(cal(1)), int(cal(2)), int(cal(3))
         end if
      end do

c     Point density near the end of the record: how many TLEs land
c     in the last 60 days?
      n = 0
      do iz = 1, npts
         if (epochs(iz) .gt. epochs(npts) - 60.d0) n = n + 1
      end do
      write(*,'(A,I5,A)') ' TLE count in last 60 days: ', n

c     Compare against the 2015-2016 era that DID qualify.
      n = 0
      do iz = 1, npts
         if (epochs(iz) .gt. 2457300.d0 .and.
     &       epochs(iz) .lt. 2457360.d0) n = n + 1
      end do
      write(*,'(A,I5,A)') ' TLE count, 60-day window in 2015-2016: ',
     &   n

c     What does a wider max_zone_days (30 days) find recently?
      write(*,*) ''
      write(*,*) '=== max_zone_days=30 sweep: how big must'//
     &   ' nzones_max be to ADMIT the 2026 zone? ==='
      do itry = 1, 6
         call zone_select(epochs, ha_o, npts,
     &      8, 30.0d0, 0.90d0, -1.0d0,
     &      try(itry), zone_start, zone_end, nzones)
         iz = 0
         do n = 1, nzones
            if (iz .eq. 0) then
               iz = n
            else if (epochs(zone_start(n)) .gt.
     &               epochs(zone_start(iz))) then
               iz = n
            end if
         end do
         write(*,'(A,I5,A,I5)') ' nzones_max=', try(itry),
     &      '  -> nzones=', nzones
         if (iz .gt. 0) then
            call jd2cal(epochs(zone_start(iz)), cal)
            write(*,'(A,I5,I3,I3)') '   latest zone starts ',
     &         int(cal(1)), int(cal(2)), int(cal(3))
         end if
      end do

      stop
      end

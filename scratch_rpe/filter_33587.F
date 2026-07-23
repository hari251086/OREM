c     filter_33587.F -- issue #27 follow-up: does tle_filter (#10)
c     flag 33587's anomalous terminal TLEs as quality issues rather
c     than real physics? The record's last TLE has n-dot=0.03130027
c     vs the first TLE's 0.00021863 (~143x) -- a concrete lead from
c     the #27 BN-refit comment (2026-07-19) that neither third-body
c     nor BN misidentification explain the observed 616->341 km
c     perigee collapse.
c
c     Compile:
c       ifx scratch_rpe\filter_33587.F tle_filter.F zone_select.F
c           tle_evolution.F ksrop\TLEread.F ksrop\Subrouts.F
c           ksrop\Legendre.F /exe:scratch_rpe\filter_33587.exe
c     -------------------------------------------------------
      program filter_33587
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(maxpts = 400)
      double precision epochs(maxpts), a_o(maxpts), e_o(maxpts)
      double precision inc_o(maxpts), raan_o(maxpts), aop_o(maxpts)
      double precision ha_o(maxpts), hp_o(maxpts), lam_o(maxpts)
      integer mask(maxpts), npts, ierr
      integer nfiltered, nremoved, ngaps
      double precision cal(6)
      integer i

      call init_constants()

      call tle_evolve('input/example_33587.tle.txt', 33587, maxpts,
     &   epochs, a_o, e_o, inc_o, raan_o, aop_o,
     &   ha_o, hp_o, lam_o, npts, ierr)
      write(*,'(A,I3,A,I5)') 'tle_evolve ierr=', ierr,
     &   '  npts=', npts

      call tle_filter(epochs, ha_o, e_o, npts,
     &   4.0d0, 6.0d0, 30.0d0, 20,
     &   mask, nfiltered, nremoved, ngaps)

      write(*,'(A,I5,A,I5,A,I5)') 'nfiltered=', nfiltered,
     &   '  nremoved=', nremoved, '  ngaps=', ngaps
      write(*,*) ''
      write(*,*) 'Last 15 points in the record (mask=0 -> flagged):'
      write(*,*) '  idx  epoch(JD)     date        ha(km)',
     &   '     e        mask'
      do i = max(1,npts-14), npts
         call jd2cal(epochs(i), cal)
         write(*,'(I5,F14.4,I6,I3,I3,F12.2,F10.5,I6)')
     &      i, epochs(i), int(cal(1)), int(cal(2)), int(cal(3)),
     &      ha_o(i), e_o(i), mask(i)
      end do

      write(*,*) ''
      write(*,*) 'All flagged points in the full record:'
      write(*,*) '  idx  epoch(JD)     date        ha(km)     e'
      do i = 1, npts
         if (mask(i) .eq. 0) then
            call jd2cal(epochs(i), cal)
            write(*,'(I5,F14.4,I6,I3,I3,F12.2,F10.5)')
     &         i, epochs(i), int(cal(1)), int(cal(2)), int(cal(3)),
     &         ha_o(i), e_o(i)
         end if
      end do

      stop
      end

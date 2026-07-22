      program diag_42985_hang
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(maxpts = 15000, max_surf = 5000)
      double precision epochs(maxpts), a_tle(maxpts)
      double precision e_tle(maxpts), inc_tle(maxpts)
      double precision raan_tle(maxpts), aop_tle(maxpts)
      double precision ha_tle(maxpts), hp_tle(maxpts)
      double precision ls_tle(maxpts)
      integer npts, ierr_tle

      integer zone_start(20), zone_end(20), nzones

      double precision ALT_a(291), DEN_a(291), SCH_a(291)
      integer natm

      double precision epz(500), haz(500), ez(500)
      double precision incz(500), raanz(500), aopz(500)
      double precision maz(500), smaz(500)

      double precision surfaces(max_surf, 3, 3)
      double precision e_grid(3), bn_grid(3)
      double precision tobs(500), apobs(500)
      integer nsurf_pts, nobs, ierr_rsm

      integer n_force(3)
      integer ie, ia, ipt, nnan, ninf

      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, 291, natm)

      call tle_evolve('input/example_42985.tle.txt', 42985, maxpts,
     &   epochs, a_tle, e_tle, inc_tle, raan_tle, aop_tle,
     &   ha_tle, hp_tle, ls_tle, npts, ierr_tle)

      write(*,*) 'npts=', npts, ' ierr_tle=', ierr_tle

      call zone_select(epochs, ha_tle, npts,
     &   8, 10.d0, 0.90d0, -1.d0, 8,
     &   zone_start, zone_end, nzones)

      write(*,*) 'nzones=', nzones
      do i = 1, nzones
         write(*,'(A,I2,A,I5,A,I5,A,I3)') '  zone ', i,
     &      ' start=', zone_start(i), ' end=', zone_end(i),
     &      ' npts=', zone_end(i)-zone_start(i)+1
      end do

      iz = nzones
      write(*,*) 'Inspecting zone ', iz

      nzone = zone_end(iz) - zone_start(iz) + 1
      do i = 1, nzone
         epz(i)   = epochs(zone_start(iz) + i - 1)
         haz(i)   = ha_tle(zone_start(iz) + i - 1)
         ez(i)    = e_tle(zone_start(iz) + i - 1)
         incz(i)  = inc_tle(zone_start(iz) + i - 1)
         raanz(i) = raan_tle(zone_start(iz) + i - 1)
         aopz(i)  = aop_tle(zone_start(iz) + i - 1)
         smaz(i)  = a_tle(zone_start(iz) + i - 1)
         maz(i)   = 0.d0
      end do

      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

      e_mid = ez(1)
      de = 0.d0
      do i = 1, nzone
         if (ez(i) .gt. e_mid + de) de = ez(i) - e_mid
         if (ez(i) .lt. e_mid - de) de = e_mid - ez(i)
      end do
      if (de .lt. 1.d-5) de = 5.d-4

      write(*,*) 'calling rsm_generate, nzone=', nzone,
     &   ' e_mid=', e_mid, ' de=', de

      call rsm_generate(
     &   epz, haz, ez, incz, raanz, aopz, maz, smaz,
     &   nzone,
     &   e_mid - de, e_mid + de,
     &   80.d0, 160.d0,
     &   20, 500, 40, 16,
     &   n_force, 20, 2, 3,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   surfaces, nsurf_pts,
     &   tobs, apobs, nobs,
     &   e_grid, bn_grid,
     &   ierr_rsm)

      write(*,*) 'ierr_rsm=', ierr_rsm, ' nsurf_pts=', nsurf_pts,
     &   ' nobs=', nobs

      nnan = 0
      ninf = 0
      do ie = 1, 3
         do ia = 1, 3
            do ipt = 1, nsurf_pts
               val = surfaces(ipt, ie, ia)
               if (val .ne. val) nnan = nnan + 1
               if (dabs(val) .gt. 1.d30) ninf = ninf + 1
            end do
         end do
      end do
      write(*,*) 'NaN count=', nnan, '  huge/Inf count=', ninf,
     &   '  total pts=', 9*nsurf_pts

      do ie = 1, 3
         do ia = 1, 3
            write(*,'(A,I1,A,I1,A)') ' surface(e=', ie, ',bn=', ia,
     &         '):'
            do ipt = 1, min(nsurf_pts, nobs)
               write(*,'(A,I3,A,F16.4,A,F16.4)') '   pt', ipt,
     &            ' model=', surfaces(ipt,ie,ia),
     &            '  obs=', apobs(ipt)
            end do
         end do
      end do

      stop
      end

c     ============================================================
      subroutine read_atm_c(alt, den, sch, maxn, natm)
      implicit double precision (a-h, o-z)
      integer maxn, natm
      double precision alt(maxn), den(maxn), sch(maxn)

      natm = 291
      if (natm .gt. maxn) natm = maxn
      alt(1) = 60.d0
      do i = 2, 141
         alt(i) = alt(i-1) + 1.d0
      end do
      do i = 142, natm
         alt(i) = alt(i-1) + 2.d0
      end do

      open(unit=88, file='input/ATM.DAT', status='old',
     &     err=900)
      read(88,110) (sch(i), i=1,61)
      read(88,110) (sch(i), i=62,natm)
      read(88,111) (den(i), i=1,61)
      read(88,111) (den(i), i=62,natm)
      close(88)
 110  format(8f10.3)
 111  format(8e10.3)

      do i = 1, natm
         alt(i) = alt(i) * 1000.d0
         sch(i) = sch(i) * 1000.d0
         den(i) = den(i) * 1.0d10
      end do
      return
 900  natm = 0
      return
      end

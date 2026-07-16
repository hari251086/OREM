c     ============================================================
c     diag_rsm.F — scratch diagnostic (issue #12 follow-up)
c     -------------------------------------------------------
c     Dumps rsm_generate's raw 9-surface envelope against the
c     observed TLE apogee for a hardcoded list of known "nobound"
c     zones (from scratch_legacy_validation/orem_results.csv), to
c     see whether the envelope failure is a magnitude/shape
c     mismatch or caused by propagate_ks failures (zero-filled
c     surface columns). Not part of the build -- scratch only.
c
c     Compile:
c       ifx /heap-arrays /F:16777216 scratch_legacy_validation\diag_rsm.F
c           rsm.F tle_evolution.F
c           ksrop\propagate_ks.F ksrop\Subrouts.F
c           ksrop\TLEread.F ksrop\Legendre.F
c           /exe:scratch_legacy_validation\diag_rsm.exe
c     -------------------------------------------------------
      program diag_rsm

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291, maxpts = 10000, max_surf = 5000)
      parameter(maxtle = 100000, ncase = 3)

      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision epochs(maxpts), a_tle(maxpts)
      double precision e_tle(maxpts), inc_tle(maxpts)
      double precision raan_tle(maxpts), aop_tle(maxpts)
      double precision ha_tle(maxpts), hp_tle(maxpts)
      double precision ls_tle(maxpts)
      integer npts, ierr_tle

      character*24 names_raw(maxtle)
      character*5 catid_raw(maxtle)
      integer norad_raw(maxtle)
      integer epoch_yr_raw(maxtle), revnum_raw(maxtle)
      double precision epoch_day_raw(maxtle)
      double precision ndot_raw(maxtle), nddot_raw(maxtle)
      double precision bstar_raw(maxtle)
      double precision ainc_raw(maxtle), raan_raw(maxtle)
      double precision ecc_raw(maxtle), aop_raw(maxtle)
      double precision ama_raw(maxtle), amm_raw(maxtle)
      integer ntle_raw

      double precision epz(500), haz(500), ez(500)
      double precision incz(500), raanz(500), aopz(500)
      double precision maz(500), smaz(500)

      double precision surfaces(max_surf, 3, 3)
      double precision e_grid(3), bn_grid(3)
      double precision tobs(500), apobs(500)
      integer nsurf_pts, nobs, ierr_rsm

      integer n_force(3)

      character*100 tle_file(ncase)
      integer norad_c(ncase)
      character*10 label(ncase)
      character*3 merge_str

      integer ic, i, ie, ia

      norad_c(1) = 42928
      tle_file(1) = 'scratch_legacy_validation/tle/42928_Z0.tle.txt'
      label(1) = '42928-Z0 (isolated, run FIRST this time)'

      norad_c(2) = 37151
      tle_file(2) = 'scratch_legacy_validation/tle/37151_Z3.tle.txt'
      label(2) = '37151-Z3'

      norad_c(3) = 32007
      tle_file(3) = 'scratch_legacy_validation/tle/32007_Z0.tle.txt'
      label(3) = '32007-Z0'

      call init_constants()
      call read_atm_v(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

      do ic = 1, ncase
         write(*,'(A)') '============================================'
         write(*,'(A,A)') 'Case: ', trim(label(ic))

         call tle_evolve(trim(tle_file(ic)), norad_c(ic), maxpts,
     &      epochs, a_tle, e_tle, inc_tle,
     &      raan_tle, aop_tle,
     &      ha_tle, hp_tle, ls_tle,
     &      npts, ierr_tle)

         if (ierr_tle .ne. 0) then
            write(*,'(A,I3)') '  tle_evolve error ', ierr_tle
            goto 200
         end if

         open(unit=77, file=trim(tle_file(ic)), status='old')
         call read_tle(77, maxtle, names_raw, catid_raw,
     &      norad_raw, epoch_yr_raw, epoch_day_raw,
     &      ndot_raw, nddot_raw, bstar_raw,
     &      ainc_raw, raan_raw, ecc_raw,
     &      aop_raw, ama_raw, amm_raw,
     &      revnum_raw, ntle_raw)
         close(77)

         write(*,'(A,I4)') '  npts = ', npts

         do i = 1, npts
            epz(i)   = epochs(i)
            haz(i)   = ha_tle(i)
            ez(i)    = e_tle(i)
            incz(i)  = inc_tle(i)
            raanz(i) = raan_tle(i)
            aopz(i)  = aop_tle(i)
            smaz(i)  = a_tle(i)
         end do
         do i = 1, npts
            maz(i) = 0.d0
            call tle_find_ma(norad_c(ic), epz(i),
     &           norad_raw, epoch_yr_raw, epoch_day_raw,
     &           ama_raw, ntle_raw, maxtle, maz(i))
         end do

         e_mid = ez(1)
         de = 0.d0
         do i = 1, npts
            if (ez(i) .gt. e_mid + de) de = ez(i) - e_mid
            if (ez(i) .lt. e_mid - de) de = e_mid - ez(i)
         end do
         if (de .lt. 1.d-5) de = 5.d-4

         write(*,'(A,F10.6,A,F10.6)')
     &      '  e_mid=', e_mid, ' de=', de

         call rsm_generate(
     &      epz, haz, ez, incz, raanz, aopz, maz, smaz,
     &      npts,
     &      e_mid - de, e_mid + de,
     &      10.d0, 300.d0,
     &      200, 360, 1.d-10,
     &      n_force, 2, 0, 0,
     &      7.2921150d-5, 3.35281066d-3, 1.d0,
     &      0.d0, 0.d0, 0, 0,
     &      0.d0, 0.d0, 0.d0,
     &      ALT_a, DEN_a, SCH_a, natm,
     &      surfaces, nsurf_pts,
     &      tobs, apobs, nobs,
     &      e_grid, bn_grid,
     &      ierr_rsm)

         write(*,'(A,I3,A,I5,A,I4)')
     &      '  ierr_rsm=', ierr_rsm, ' nsurf_pts=', nsurf_pts,
     &      ' nobs=', nobs
         write(*,'(A,3F10.4)') '  e_grid  = ', (e_grid(ie),ie=1,3)
         write(*,'(A,3F10.4)') '  bn_grid = ', (bn_grid(ia),ia=1,3)

         if (ierr_rsm .ne. 0) goto 200

c        Flag which of the 9 grid points look zero-filled
c        (propagate_ks failure proxy: first surface point is
c        exactly 0.d0 -- real apogee altitudes are never exactly
c        zero).
         write(*,'(A)') '  grid-point status (F=zero-filled/failed):'
         do ie = 1, 3
            write(*,'(A,3(A,L1,A))')
     &         '    ',
     &         ' (', surfaces(1,ie,1).eq.0.d0, ')',
     &         ' (', surfaces(1,ie,2).eq.0.d0, ')',
     &         ' (', surfaces(1,ie,3).eq.0.d0, ')'
         end do

         write(*,'(A)')
     &      '  obs#   t(day)     apobs      smin       smax'//
     &      '   bounded?'
         do i = 1, nobs
            smin = surfaces(i,1,1)
            smax = surfaces(i,1,1)
            do ie = 1, 3
               do ia = 1, 3
                  if (surfaces(i,ie,ia) .lt. smin)
     &               smin = surfaces(i,ie,ia)
                  if (surfaces(i,ie,ia) .gt. smax)
     &               smax = surfaces(i,ie,ia)
               end do
            end do
            write(*,'(I5,F10.3,F11.3,F11.3,F11.3,A)')
     &         i, tobs(i), apobs(i), smin, smax,
     &         merge_str(apobs(i).ge.smin .and. apobs(i).le.smax)
         end do

 200     continue
      end do

      stop
      end

c     ============================================================
      character*3 function merge_str(flag)
      logical flag
      if (flag) then
         merge_str = 'YES'
      else
         merge_str = 'NO'
      end if
      return
      end

c     ============================================================
      subroutine read_atm_v(alt, den, sch, maxn, natm)
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

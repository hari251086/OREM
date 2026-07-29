c     ============================================================
c     diag_interp_vs_direct.F -- issue #36, Finding 3 isolation
c
c     Reproduces zone 1 of 42928 exactly as production does, then
c     calls rsm_generate directly to get the SAME 9-corner surfaces
c     the real GA search saw, and ga_twoint to interpolate at the
c     already-fitted (e_opt, bn_opt) -- reproducing EXACTLY what
c     ga_fitness computed internally during the real GA run (not an
c     approximation of it). Compares that interpolated value against
c     the zone's real TLE apogee (should closely match production's
c     own rms_fit, now correctly reported in km post-Finding-2-fix)
c     and against the direct-propagation result already measured in
c     diag_osc_vs_mean_apogee_v2.F (23.1 km RMS vs TLE).
c
c     Grid bounds: e_min/e_max reproduced exactly (deterministic from
c     zone TLE data, same formula orem.F uses). bn_lo/bn_hi use the
c     config default [80,160] unmodified -- zone 1 has no trust-gate
c     carryover history (it's the first zone), so the only possible
c     production modification is G2's floor/G3's BSTAR narrowing,
c     neither reconstructed here; a reasonable, disclosed
c     approximation for this isolation test, not an exact replay.
c
c     Read-only diagnostic. No production code touched.
c     ============================================================
      program diag_interp_vs_direct
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(maxpts = 15000, maxtle = 100000)
      character*80 tle_file
      integer norad_id

      double precision epochs(maxpts), a_tle(maxpts)
      double precision e_tle(maxpts), inc_tle(maxpts)
      double precision raan_tle(maxpts), aop_tle(maxpts)
      double precision ha_tle(maxpts), hp_tle(maxpts)
      double precision ls_tle(maxpts)
      integer npts, ierr_tle

      integer mask_tle(maxpts)
      integer nfiltered_tle, nremoved_tle, ngaps_tle
      integer i_filt, j_filt

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

      integer nzones_max
      parameter(nzones_max = 8)
      integer zone_start(nzones_max), zone_end(nzones_max)
      integer nzones

      double precision epz(500), haz(500), ez(500)
      double precision incz(500), raanz(500), aopz(500)
      double precision maz(500), smaz(500)

      parameter(max_surf = 5000)
      double precision surfaces(max_surf, 3, 3)
      double precision tobs(500), apobs(500)
      double precision e_grid(3), bn_grid(3)
      integer nsurf_pts, nobs, ierr_rsm

      double precision alt_a(291), den_a(291), sch_a(291)
      integer natm, n_force(3)

      double precision e_opt_zone1, bn_opt_zone1
      double precision f2d(3,3), fout
      integer iz, i, nzone1, ierr_sw, ierr_a2, ih, ia

      call init_constants()

      tle_file = 'input/example_42928.tle.txt'
      norad_id = 42928
      e_opt_zone1  = 0.321707d0
      bn_opt_zone1 = 81.195d0

      call tle_evolve(tle_file, norad_id, maxpts,
     &   epochs, a_tle, e_tle, inc_tle,
     &   raan_tle, aop_tle,
     &   ha_tle, hp_tle, ls_tle,
     &   npts, ierr_tle)
      if (ierr_tle .ne. 0) stop 1

      call tle_filter(epochs, ha_tle, e_tle, npts,
     &   4.0d0, 6.0d0, 30.0d0, 20,
     &   mask_tle, nfiltered_tle, nremoved_tle, ngaps_tle)
      if (nremoved_tle .gt. 0) then
         j_filt = 0
         do i_filt = 1, npts
            if (mask_tle(i_filt) .eq. 1) then
               j_filt = j_filt + 1
               epochs(j_filt)   = epochs(i_filt)
               a_tle(j_filt)    = a_tle(i_filt)
               e_tle(j_filt)    = e_tle(i_filt)
               inc_tle(j_filt)  = inc_tle(i_filt)
               raan_tle(j_filt) = raan_tle(i_filt)
               aop_tle(j_filt)  = aop_tle(i_filt)
               ha_tle(j_filt)   = ha_tle(i_filt)
               hp_tle(j_filt)   = hp_tle(i_filt)
               ls_tle(j_filt)   = ls_tle(i_filt)
            end if
         end do
         npts = j_filt
      end if

      open(unit=77, file=tle_file, status='old')
      call read_tle(77, maxtle, names_raw, catid_raw,
     &   norad_raw, epoch_yr_raw, epoch_day_raw,
     &   ndot_raw, nddot_raw, bstar_raw,
     &   ainc_raw, raan_raw, ecc_raw,
     &   aop_raw, ama_raw, amm_raw,
     &   revnum_raw, ntle_raw)
      close(77)

      call zone_select(epochs, ha_tle, npts,
     &   8, 20.0d0, 0.90d0, -1.0d0,
     &   nzones_max,
     &   zone_start, zone_end, nzones)
      if (nzones .lt. 1) stop 1

      nzone1 = zone_end(1) - zone_start(1) + 1
      do i = 1, nzone1
         epz(i)   = epochs(zone_start(1) + i - 1)
         haz(i)   = ha_tle(zone_start(1) + i - 1)
         ez(i)    = e_tle(zone_start(1) + i - 1)
         incz(i)  = inc_tle(zone_start(1) + i - 1)
         raanz(i) = raan_tle(zone_start(1) + i - 1)
         aopz(i)  = aop_tle(zone_start(1) + i - 1)
         smaz(i)  = a_tle(zone_start(1) + i - 1)
         call tle_find_ma(norad_id, epz(i),
     &        norad_raw, epoch_yr_raw, epoch_day_raw,
     &        ama_raw, ntle_raw, maxtle, maz(i))
      end do

c     e_min/e_max exactly as orem.F's own zone loop computes them
      e_mid = ez(1)
      de = 0.d0
      do i = 1, nzone1
         if (ez(i) .gt. e_mid + de) de = ez(i) - e_mid
         if (ez(i) .lt. e_mid - de) de = e_mid - ez(i)
      end do
      if (de .lt. 1.d-5) de = 5.d-4
      write(*,'(A,F10.6,A,F10.6,A,F10.6)') 'e_mid=', e_mid,
     &   ' e_min=', e_mid-de, ' e_max=', e_mid+de

      call read_atm_local(alt_a, den_a, sch_a, natm)
      call sw_load('input/SW-All.csv', ierr_sw)
      call atm2d_load('input/ATM2D.DAT', ierr_a2)
      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

      write(*,'(A)') 'Calling rsm_generate (9-corner grid, bn=[80,160])'
      call rsm_generate(
     &   epz, haz, ez, incz, raanz, aopz, maz, smaz,
     &   nzone1,
     &   e_mid - de, e_mid + de,
     &   80.d0, 160.d0,
     &   200, 360, 1.d-10,
     &   n_force, 20, 2, 3,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   alt_a, den_a, sch_a, natm,
     &   surfaces, nsurf_pts,
     &   tobs, apobs, nobs,
     &   e_grid, bn_grid,
     &   ierr_rsm)
      write(*,'(A,I3,A,I5,A,I5)') '  ierr_rsm=', ierr_rsm,
     &   ' nsurf_pts=', nsurf_pts, ' nobs=', nobs
      write(*,'(A,3F10.4)') '  e_grid=', e_grid
      write(*,'(A,3F10.4)') '  bn_grid=', bn_grid
      if (ierr_rsm .ne. 0) stop 1

      write(*,'(A)') ''
      write(*,'(A)')
     &  'epoch_day  ha_TLE(truth)  interp_pred   err_interp'
      ssq = 0.d0
      amax = 0.d0
      do i = 1, nobs
         do ih = 1, 3
            do ia = 1, 3
               f2d(ih, ia) = surfaces(i, ih, ia)
            end do
         end do
         call ga_twoint(e_grid, bn_grid, f2d,
     &                  e_opt_zone1, bn_opt_zone1, fout, 3, 3)
         err = fout - apobs(i)
         ssq = ssq + err**2
         if (dabs(err) .gt. amax) amax = dabs(err)
         write(*,'(F9.4,F14.4,F13.4,F12.4)')
     &      tobs(i), apobs(i), fout, err
      end do

      rms = dsqrt(ssq/dble(nobs))
      write(*,'(A)') ''
      write(*,'(A,F10.4,A)') 'RMS interpolated-vs-TLE = ', rms, ' km'
      write(*,'(A,F10.4,A)') 'MAX interpolated-vs-TLE = ', amax, ' km'
      write(*,'(A)') ''
      write(*,'(A)') 'Compare: production rms_fit (post-Finding-2-fix,'
      write(*,'(A)') '  same zone, real historical bn bounds) = 4.42 km'
      write(*,'(A)') 'Compare: direct-propagation RMS vs TLE'
      write(*,'(A)') '  (diag_osc_vs_mean_apogee_v2.F) = 23.08-23.12 km'

      stop
      end

c     ============================================================
      subroutine read_atm_local(alt, den, sch, natm)
      implicit double precision (a-h, o-z)
      parameter(ndim_a = 291)
      double precision alt(ndim_a), den(ndim_a), sch(ndim_a)
      integer natm, i
      natm = ndim_a
      alt(1) = 60.d0
      do i = 2, 141
         alt(i) = alt(i-1) + 1.d0
      end do
      do i = 142, natm
         alt(i) = alt(i-1) + 2.d0
      end do
      open(unit=60, file='input/ATM.DAT', status='old')
      read(60,110) (sch(i), i=1,61)
      read(60,110) (sch(i), i=62,natm)
      read(60,111) (den(i), i=1,61)
      read(60,111) (den(i), i=62,natm)
      close(60)
 110  format(8f10.3)
 111  format(8e10.3)
      do i = 1, natm
         alt(i) = alt(i) * 1000.d0
         sch(i) = sch(i) * 1000.d0
         den(i) = den(i) * 1.0d10
      end do
      return
      end

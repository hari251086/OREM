c     ============================================================
c     diag_zone1_faithful.F -- issue #36, faithful zone-1 replay
c
c     Corrects two confounds found in the earlier diagnostics
c     (diag_osc_vs_mean_apogee[_v2].F, diag_interp_vs_direct.F):
c       1. They seeded the propagation IC from RAW MEAN TLE elements
c          at zone-1 index-1. Current production (src/orem.F) applies
c          tle_find_osc's SGP4-osculating substitution at that same
c          point, for BOTH the G2 floor calibration AND the main
c          RSM/GA feed -- this diagnostic replicates that exactly.
c       2. They used the raw config BN bounds [80,160] for the RSM
c          grid. Production narrows/floors those via G2
c          (estimate_bn_floor) then G3 (estimate_bn_bstar_prior)
c          before rsm_generate ever sees them -- this diagnostic
c          reproduces that exact sequence by calling the real
c          production subroutines directly (not reimplementing them).
c
c     Produces, all apples-to-apples on the SAME faithfully-
c     reconstructed IC and grid bounds: (a) RSM-interpolated apogee
c     vs TLE truth [should now closely match production's real
c     rms_fit=4.42km if the reconstruction is faithful], (b) direct
c     (non-interpolated) propagation at (e_opt,bn_opt) vs TLE truth
c     [the real test of whether the fitted optimum itself tracks
c     reality], (c) within-revolution mean-averaged (via
c     propagate_ks_dense) vs TLE truth [re-test of Finding 1 with the
c     corrected IC].
c
c     Read-only diagnostic. No production code touched.
c     ============================================================
      program diag_zone1_faithful
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
      integer nzone1, i, j, ih, ia

      double precision a_kep_osc1(7)
      logical osc_ok, bn_floor_ok, bn_prior_ok

      double precision bn_lo, bn_hi, bn_lo_zone, bn_hi_zone
      double precision bn_floor_est

      integer n_force(3)
      double precision alt_a(291), den_a(291), sch_a(291)
      integer natm, ierr_sw, ierr_a2

      parameter(max_surf = 5000)
      double precision surfaces(max_surf, 3, 3)
      double precision tobs(500), apobs(500)
      double precision e_grid(3), bn_grid(3)
      integer nsurf_pts, nobs, ierr_rsm

      double precision e_opt_zone1, bn_opt_zone1
      double precision f2d(3,3), fout

      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      parameter(max_traj = 150000)
      double precision traj_jd(max_traj)
      double precision traj_x(3, max_traj)
      double precision traj_xd(3, max_traj)
      double precision osc_ha(max_traj), avg_ha(max_traj)
      integer idump_loc, exit_code_re, ntraj, istep_val

      call init_constants()

      tle_file = 'input/example_42928.tle.txt'
      norad_id = 42928
c     Real production-fitted values, current campaign CSV
      e_opt_zone1  = 0.321707d0
      bn_opt_zone1 = 81.195d0

      call tle_evolve(tle_file, norad_id, maxpts,
     &   epochs, a_tle, e_tle, inc_tle, raan_tle, aop_tle,
     &   ha_tle, hp_tle, ls_tle, npts, ierr_tle)
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
     &   nzones_max, zone_start, zone_end, nzones)
      if (nzones .lt. 1) stop 1
      nzone1 = zone_end(1) - zone_start(1) + 1

      call read_atm_local(alt_a, den_a, sch_a, natm)
      call sw_load('input/SW-All.csv', ierr_sw)
      call atm2d_load('input/ATM2D.DAT', ierr_a2)
      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

c     ============================================================
c     G2: physics-based BN floor, EXACTLY orem.F's own sequence
c     ============================================================
      do i = 1, nzone1
         epz(i)   = epochs(zone_start(1) + i - 1)
         ez(i)    = e_tle(zone_start(1) + i - 1)
         incz(i)  = inc_tle(zone_start(1) + i - 1)
         raanz(i) = raan_tle(zone_start(1) + i - 1)
         aopz(i)  = aop_tle(zone_start(1) + i - 1)
         smaz(i)  = a_tle(zone_start(1) + i - 1)
      end do
      do j = 1, nzone1
         maz(j) = 0.d0
         call tle_find_ma(norad_id, epz(j),
     &        norad_raw, epoch_yr_raw, epoch_day_raw,
     &        ama_raw, ntle_raw, maxtle, maz(j))
      end do
      call tle_find_osc(norad_id, epz(1),
     &   norad_raw, epoch_yr_raw, epoch_day_raw,
     &   ainc_raw, raan_raw, ecc_raw, aop_raw, ama_raw,
     &   amm_raw, bstar_raw, ntle_raw, maxtle,
     &   a_kep_osc1, osc_ok)
      if (osc_ok) then
         smaz(1)  = a_kep_osc1(1)
         ez(1)    = a_kep_osc1(2)
         incz(1)  = a_kep_osc1(3)
         raanz(1) = a_kep_osc1(4)
         aopz(1)  = a_kep_osc1(5)
         maz(1)   = a_kep_osc1(6)
      end if

      bn_lo = 80.d0
      bn_hi = 160.d0
      call estimate_bn_floor(epz, smaz, ez, incz, raanz, aopz, maz,
     &   nzone1, n_force, 20, 2, 3,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2, 4.56d-6,
     &   1.32712440018d11, 4.902801076d3,
     &   alt_a, den_a, sch_a, natm,
     &   bn_floor_est, bn_floor_ok)
      if (bn_floor_ok) then
         write(*,'(A,F8.2)') 'G2 bn_floor_est=', bn_floor_est
         if (0.5d0*bn_floor_est .lt. bn_lo) then
            bn_lo = 0.5d0*bn_floor_est
            write(*,'(A,F8.2)') '  -> bn_lo=', bn_lo
         end if
      else
         write(*,'(A)') 'G2: not ok, bn_lo unchanged'
      end if

c     ============================================================
c     Main per-zone loop, zone 1: REBUILD fresh (mean elements),
c     then the SAME osculating substitution at index 1 -- exactly
c     orem.F's own main loop sequence.
c     ============================================================
      do i = 1, nzone1
         epz(i)   = epochs(zone_start(1) + i - 1)
         haz(i)   = ha_tle(zone_start(1) + i - 1)
         ez(i)    = e_tle(zone_start(1) + i - 1)
         incz(i)  = inc_tle(zone_start(1) + i - 1)
         raanz(i) = raan_tle(zone_start(1) + i - 1)
         aopz(i)  = aop_tle(zone_start(1) + i - 1)
         smaz(i)  = a_tle(zone_start(1) + i - 1)
      end do
      do j = 1, nzone1
         maz(j) = 0.d0
         call tle_find_ma(norad_id, epz(j),
     &        norad_raw, epoch_yr_raw, epoch_day_raw,
     &        ama_raw, ntle_raw, maxtle, maz(j))
      end do
      call tle_find_osc(norad_id, epz(1),
     &   norad_raw, epoch_yr_raw, epoch_day_raw,
     &   ainc_raw, raan_raw, ecc_raw, aop_raw, ama_raw,
     &   amm_raw, bstar_raw, ntle_raw, maxtle,
     &   a_kep_osc1, osc_ok)
      if (osc_ok) then
         smaz(1)  = a_kep_osc1(1)
         ez(1)    = a_kep_osc1(2)
         incz(1)  = a_kep_osc1(3)
         raanz(1) = a_kep_osc1(4)
         aopz(1)  = a_kep_osc1(5)
         maz(1)   = a_kep_osc1(6)
      end if
      write(*,'(A,L1,A,F10.6,A,F12.4)') 'osc_ok=', osc_ok,
     &   ' osculating ez(1)=', ez(1), ' smaz(1)=', smaz(1)

      e_mid = ez(1)
      de = 0.d0
      do i = 1, nzone1
         if (ez(i) .gt. e_mid + de) de = ez(i) - e_mid
         if (ez(i) .lt. e_mid - de) de = e_mid - ez(i)
      end do
      if (de .lt. 1.d-5) de = 5.d-4

c     ============================================================
c     G3: BSTAR-informed narrowing, EXACTLY orem.F's own sequence
c     DEBUG: manually extract BSTAR per zone-1 epoch first, to see
c     WHY estimate_bn_bstar_prior returns ok or not-ok.
c     ============================================================
      nb_dbg = 0
      do i = 1, nzone1
         call tle_find_bstar(norad_id, epz(i),
     &      norad_raw, epoch_yr_raw, epoch_day_raw, bstar_raw,
     &      ntle_raw, maxtle, bstar_dbg)
         if (bstar_dbg .gt. 0.d0) then
            nb_dbg = nb_dbg + 1
            if (nb_dbg .eq. 1) write(*,'(A,F14.8)')
     &         'DEBUG first valid bstar=', bstar_dbg
         end if
      end do
      write(*,'(A,I4,A,I4)') 'DEBUG nb_valid_bstar=', nb_dbg,
     &   ' of nzone1=', nzone1

      call estimate_bn_bstar_prior(norad_id, epz, nzone1,
     &   norad_raw, epoch_yr_raw, epoch_day_raw, bstar_raw,
     &   ntle_raw, maxtle,
     &   bn_lo, bn_hi, bn_lo_zone, bn_hi_zone, bn_prior_ok)
      if (bn_prior_ok) then
         write(*,'(A,F8.2,A,F8.2)') 'G3 narrows BN to ',
     &      bn_lo_zone, ' .. ', bn_hi_zone
      else
         bn_lo_zone = bn_lo
         bn_hi_zone = bn_hi
         write(*,'(A,F8.2,A,F8.2)') 'G3: not ok, BN stays ',
     &      bn_lo_zone, ' .. ', bn_hi_zone
      end if

c     ============================================================
c     RSM grid + GA-equivalent interpolation at the real fitted
c     optimum, using the FAITHFULLY reconstructed bounds/IC
c     ============================================================
      write(*,'(A)') ''
      write(*,'(A)') 'Calling rsm_generate (faithful bounds)'
      call rsm_generate(
     &   epz, haz, ez, incz, raanz, aopz, maz, smaz,
     &   nzone1,
     &   e_mid - de, e_mid + de,
     &   bn_lo_zone, bn_hi_zone,
     &   200, 360, 1.d-10,
     &   n_force, 20, 2, 3,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   alt_a, den_a, sch_a, natm,
     &   surfaces, nsurf_pts, tobs, apobs, nobs,
     &   e_grid, bn_grid, ierr_rsm)
      write(*,'(A,I3,A,I5)') '  ierr_rsm=', ierr_rsm, ' nobs=', nobs
      write(*,'(A,3F10.4)') '  e_grid=', e_grid
      write(*,'(A,3F10.4)') '  bn_grid=', bn_grid
      if (ierr_rsm .ne. 0) stop 1

      ssq_i = 0.d0
      amax_i = 0.d0
      do i = 1, nobs
         do ih = 1, 3
            do ia = 1, 3
               f2d(ih, ia) = surfaces(i, ih, ia)
            end do
         end do
         call ga_twoint(e_grid, bn_grid, f2d,
     &                  e_opt_zone1, bn_opt_zone1, fout, 3, 3)
         erri = fout - apobs(i)
         ssq_i = ssq_i + erri**2
         if (dabs(erri) .gt. amax_i) amax_i = dabs(erri)
      end do
      rms_interp = dsqrt(ssq_i/dble(nobs))
      write(*,'(A,F10.4,A)') 'RMS interpolated-vs-TLE (faithful) = ',
     &   rms_interp, ' km'
      write(*,'(A,F10.4,A)') 'MAX interpolated-vs-TLE (faithful) = ',
     &   amax_i, ' km'

c     ============================================================
c     Direct propagation at the exact fitted optimum, faithful IC
c     ============================================================
      a_kep(1) = smaz(1)
      a_kep(2) = e_opt_zone1
      a_kep(3) = incz(1)
      a_kep(4) = raanz(1)
      a_kep(5) = aopz(1)
      a_kep(6) = maz(1)
      a_kep(7) = 0.d0
      rp = smaz(1) * (1.d0 - ez(1))
      a_kep(1) = rp / (1.d0 - e_opt_zone1)

      call oe2car(a_kep, x0, xd0, 1.d-12)
      call jd2cal(epz(1), cal0)

      period_s = 2.d0*pi*dsqrt(a_kep(1)**3/amue)
      period_d = period_s / 86400.d0
      zone_dur = epz(nzone1) - epz(1)
      nrev_zone = int(zone_dur / period_d) + 2
      istep_val = 360

      idump_loc = 1
      do i = 1, max_traj
         traj_jd(i) = 0.d0
      end do
      call propagate_ks_dense(
     &   x0, xd0, cal0,
     &   nrev_zone, istep_val, 1.d-10,
     &   n_force, 20, 2, 3,
     &   bn_opt_zone1, 1, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   alt_a, den_a, sch_a, natm,
     &   max_traj, idump_loc,
     &   traj_jd, traj_x, traj_xd,
     &   exit_code_re)
      write(*,'(A)') ''
      write(*,'(A,I3,A,I8)') 'propagate_ks_dense exit_code=',
     &   exit_code_re, ' idump=', idump_loc

      ntraj = 0
      do i = 1, max_traj
         if (traj_jd(i) .le. 0.d0) goto 16
         if (traj_jd(i) - epz(1) .gt. zone_dur + 0.5d0) goto 16
         ntraj = ntraj + 1
         r_mag = dsqrt(traj_x(1,i)**2+traj_x(2,i)**2+traj_x(3,i)**2)
         v_mag = dsqrt(traj_xd(1,i)**2+traj_xd(2,i)**2
     &                +traj_xd(3,i)**2)
         sma_pt = 1.d0/(2.d0/r_mag - v_mag**2/amue)
         hx = traj_x(2,i)*traj_xd(3,i) - traj_x(3,i)*traj_xd(2,i)
         hy = traj_x(3,i)*traj_xd(1,i) - traj_x(1,i)*traj_xd(3,i)
         hz = traj_x(1,i)*traj_xd(2,i) - traj_x(2,i)*traj_xd(1,i)
         h_mag = dsqrt(hx**2+hy**2+hz**2)
         p_sl = h_mag**2/amue
         ecc_pt = dsqrt(dabs(1.d0 - p_sl/sma_pt))
         osc_ha(ntraj) = sma_pt*(1.d0+ecc_pt) - R_Earth
      end do
 16   continue

      nhalf = istep_val / 2
      do i = 1, ntraj
         ilo = i - nhalf
         ihi = i + nhalf
         if (ilo .lt. 1) ilo = 1
         if (ihi .gt. ntraj) ihi = ntraj
         ssum = 0.d0
         ncnt = 0
         do k = ilo, ihi
            ssum = ssum + osc_ha(k)
            ncnt = ncnt + 1
         end do
         avg_ha(i) = ssum / dble(ncnt)
      end do

      ssq_d = 0.d0
      ssq_a = 0.d0
      amax_d = 0.d0
      amax_a = 0.d0
      do i = 1, nzone1
         t_want = epz(i) - epz(1)
         call rsm_tinterp(traj_jd, osc_ha, ntraj, epz(1), t_want,
     &                     y_osc)
         call rsm_tinterp(traj_jd, avg_ha, ntraj, epz(1), t_want,
     &                     y_avg)
         err_d = y_osc - haz(i)
         err_a = y_avg - haz(i)
         ssq_d = ssq_d + err_d**2
         ssq_a = ssq_a + err_a**2
         if (dabs(err_d) .gt. amax_d) amax_d = dabs(err_d)
         if (dabs(err_a) .gt. amax_a) amax_a = dabs(err_a)
      end do
      rms_direct = dsqrt(ssq_d/dble(nzone1))
      rms_avg    = dsqrt(ssq_a/dble(nzone1))

      write(*,'(A)') ''
      write(*,'(A,F10.4,A)') 'RMS direct-propagation-vs-TLE = ',
     &   rms_direct, ' km'
      write(*,'(A,F10.4,A)') 'RMS within-rev-averaged-vs-TLE = ',
     &   rms_avg, ' km'
      write(*,'(A)') ''
      write(*,'(A,F10.4,A,F10.4,A,F10.4,A)')
     &   'SUMMARY (faithful IC/bounds): interp=', rms_interp,
     &   'km direct=', rms_direct, 'km averaged=', rms_avg, 'km'

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

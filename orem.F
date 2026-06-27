c     ============================================================
c     orem.F — OREM Driver: Optimal Regularized re-Entry Method
c
c     Full pipeline: TLE → zone selection → RSM → GA → re-entry
c     prediction with iterative A-bound narrowing across zones.
c
c     Author: Harishkumar Sellamuthu
c     ============================================================

      subroutine orem_run(
     &   tle_file, norad_id,
     &   t_obs_cal,
     &   nzones_max, min_zone_pts, max_zone_days,
     &   r2_thresh, slope_thresh,
     &   a_min_init, a_max_init,
     &   fixed_mass, cd_drag, idrag_flag,
     &   ipopsize, maxgen, nbits_e, nbits_a,
     &   pcross, pmute, ga_seed,
     &   ngeo_deg, nsun_deg, nmoon_deg,
     &   WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD,
     &   PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   reentry_jd, e_opt_out, a_opt_out, rms_out,
     &   zone_epoch, nzones_used,
     &   rpe_out, t_mean, t_std,
     &   ierr)
c
c     Input
c       tle_file         : TLE file path
c       norad_id         : target NORAD catalog number
c       t_obs_cal(6)     : observed re-entry [yr,mo,dy,hr,mn,sc]
c                          (all zeros = operational mode)
c       zone selection params, area bounds, mass, Cd, GA params,
c       force model params, atmosphere tables
c
c     Output
c       reentry_jd(nzones_max) : predicted re-entry JD per zone
c       e_opt_out(nzones_max)  : optimal eccentricity per zone
c       a_opt_out(nzones_max)  : optimal area (m²) per zone
c       rms_out(nzones_max)    : GA RMS per zone
c       zone_epoch(nzones_max) : zone start epoch (JD)
c       nzones_used            : number of zones processed
c       rpe_out(nzones_max)    : RPE per zone (%)
c       t_mean, t_std          : mean/std re-entry prediction
c       ierr                   : 0=ok, 1=TLE error, 2=no zones,
c                                3=all zones failed
c     ============================================================
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

c     --- Input arguments ---
      character*(*) tle_file
      integer norad_id, nzones_max
      integer min_zone_pts
      double precision max_zone_days, r2_thresh, slope_thresh
      double precision t_obs_cal(6)
      double precision a_min_init, a_max_init
      double precision fixed_mass, cd_drag
      integer idrag_flag
      integer ipopsize, maxgen, nbits_e, nbits_a
      double precision pcross, pmute, ga_seed
      integer ngeo_deg, nsun_deg, nmoon_deg
      double precision WE_rot, EPS_f, FR_rot
      double precision CR_srp, AM_srp
      integer IPSR, ISHAD
      double precision PSR_srp, amuS, amuM
      integer ndim_atm
      double precision ALT_atm(ndim_atm)
      double precision DEN_atm(ndim_atm)
      double precision SCH_atm(ndim_atm)

c     --- Output arguments ---
      double precision reentry_jd(nzones_max)
      double precision e_opt_out(nzones_max)
      double precision a_opt_out(nzones_max)
      double precision rms_out(nzones_max)
      double precision zone_epoch(nzones_max)
      integer nzones_used
      double precision rpe_out(nzones_max)
      double precision t_mean, t_std
      integer ierr

c     --- Internal arrays ---
      parameter(maxpts = 10000, max_surf = 5000)
      parameter(maxtle = 100000)

c     TLE evolution
      double precision epochs(maxpts), a_tle(maxpts)
      double precision e_tle(maxpts), inc_tle(maxpts)
      double precision raan_tle(maxpts), aop_tle(maxpts)
      double precision ha_tle(maxpts), hp_tle(maxpts)
      double precision ls_tle(maxpts)
      integer npts, ierr_tle

c     Raw TLE for mean anomaly extraction
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
      save names_raw, catid_raw, norad_raw
      save epoch_yr_raw, epoch_day_raw, revnum_raw
      save ndot_raw, nddot_raw, bstar_raw
      save ainc_raw, raan_raw, ecc_raw
      save aop_raw, ama_raw, amm_raw

c     Zone selection
      integer zone_start(nzones_max), zone_end(nzones_max)
      integer nzones

c     RSM
      double precision surfaces(max_surf, 3, 3)
      double precision e_grid(3), a_grid(3)
      double precision tobs(500), apobs(500)
      integer nsurf_pts, nobs, ierr_rsm
      save surfaces

c     Zone data
      double precision epz(500), haz(500), ez(500)
      double precision incz(500), raanz(500), aopz(500)
      double precision maz(500), smaz(500)

c     Propagator output for re-entry
      parameter(max_reentry = 50000)
      double precision traj_jd(max_reentry)
      double precision traj_x(3, max_reentry)
      double precision traj_xd(3, max_reentry)
      integer exit_code_re, idump_re
      save traj_jd, traj_x, traj_xd

c     Force model
      integer n_force(3)

c     Working
      double precision a_lo, a_hi, e_mid, de
      double precision bn_val, t_obs_jd
      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      integer iz, i, nzone, nvalid

      ierr = 0
      nzones_used = 0
      t_mean = 0.d0
      t_std = 0.d0

c     ============================================================
c     Step 1: Load TLE evolution
c     ============================================================
      call tle_evolve(tle_file, norad_id, maxpts,
     &   epochs, a_tle, e_tle, inc_tle,
     &   raan_tle, aop_tle,
     &   ha_tle, hp_tle, ls_tle,
     &   npts, ierr_tle)

      if (ierr_tle .ne. 0) then
         ierr = 1
         return
      end if

c     ============================================================
c     Step 2: Read raw TLE for mean anomaly
c     ============================================================
      open(unit=77, file=tle_file, status='old', err=900)
      call read_tle(77, maxtle, names_raw, catid_raw,
     &   norad_raw, epoch_yr_raw, epoch_day_raw,
     &   ndot_raw, nddot_raw, bstar_raw,
     &   ainc_raw, raan_raw, ecc_raw,
     &   aop_raw, ama_raw, amm_raw,
     &   revnum_raw, ntle_raw)
      close(77)

c     ============================================================
c     Step 3: Zone selection
c     ============================================================
      call zone_select(epochs, ha_tle, npts,
     &   min_zone_pts, max_zone_days,
     &   r2_thresh, slope_thresh,
     &   nzones_max,
     &   zone_start, zone_end, nzones)

      if (nzones .eq. 0) then
         ierr = 2
         return
      end if

c     ============================================================
c     Step 4: Convert observed re-entry to JD (if provided)
c     ============================================================
      t_obs_jd = 0.d0
      if (t_obs_cal(1) .gt. 1900.d0) then
         call cal2jd(t_obs_cal, t_obs_jd)
      end if

c     ============================================================
c     Step 5: Zone-by-zone RSM → GA → re-entry prediction
c     ============================================================
      a_lo = a_min_init
      a_hi = a_max_init
      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1
      nvalid = 0

      do iz = 1, nzones

c        Extract TLE data from zone start to end of evolution
c        (RSM surfaces span full remaining lifetime, not just zone)
         nzone = npts - zone_start(iz) + 1
         if (nzone .gt. 500) nzone = 500
         do i = 1, nzone
            epz(i) = epochs(zone_start(iz) + i - 1)
            haz(i) = ha_tle(zone_start(iz) + i - 1)
            ez(i) = e_tle(zone_start(iz) + i - 1)
            incz(i) = inc_tle(zone_start(iz) + i - 1)
            raanz(i) = raan_tle(zone_start(iz) + i - 1)
            aopz(i) = aop_tle(zone_start(iz) + i - 1)
            smaz(i) = a_tle(zone_start(iz) + i - 1)
         end do

c        Find mean anomaly from raw TLE for each zone point
         do j = 1, nzone
            maz(j) = 0.d0
            call tle_find_ma(norad_id, epz(j),
     &           norad_raw, epoch_yr_raw, epoch_day_raw,
     &           ama_raw, ntle_raw, maxtle, maz(j))
         end do

c        Compute e-grid from TLE scatter
         e_mid = ez(1)
         de = 0.d0
         do i = 1, nzone
            if (ez(i) .gt. e_mid + de) de = ez(i) - e_mid
            if (ez(i) .lt. e_mid - de) de = e_mid - ez(i)
         end do
         if (de .lt. 1.d-5) de = 5.d-4

         write(*,'(A,I2,A,F10.4,A,F8.4,A,F8.2,A,F8.2,
     &      A,I4)')
     &      '  Zone ', iz, ': MA=', maz(1),
     &      ' e=', e_mid, ' A=[', a_lo, ',', a_hi,
     &      '] nz=', nzone

c        Generate RSM surfaces
         call rsm_generate(
     &      epz, haz, ez, incz, raanz, aopz, maz, smaz,
     &      nzone,
     &      e_mid - de, e_mid + de,
     &      a_lo, a_hi,
     &      fixed_mass, cd_drag, idrag_flag,
     &      10000, 1, 1.d-10,
     &      n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &      WE_rot, EPS_f, FR_rot,
     &      CR_srp, AM_srp, IPSR, ISHAD,
     &      PSR_srp, amuS, amuM,
     &      ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &      surfaces, nsurf_pts,
     &      tobs, apobs, nobs,
     &      e_grid, a_grid,
     &      ierr_rsm)

         write(*,'(A,I2,A,I2,A,I4)')
     &      '  RSM: ierr=', ierr_rsm,
     &      ' nsurf_pts=', nsurf_pts,
     &      ' nobs=', nobs

         if (ierr_rsm .ne. 0 .or. nsurf_pts .le. 0) then
            reentry_jd(iz) = 0.d0
            e_opt_out(iz) = 0.d0
            a_opt_out(iz) = 0.d0
            rms_out(iz) = -1.d0
            zone_epoch(iz) = epz(1)
            goto 50
         end if

c        GA optimization
         call ga_optimize(
     &      surfaces, e_grid, a_grid,
     &      3, 3, nsurf_pts,
     &      tobs, apobs, nobs,
     &      e_grid(1), e_grid(3),
     &      a_grid(1), a_grid(3),
     &      ipopsize, maxgen, nbits_e, nbits_a,
     &      pcross, pmute, ga_seed,
     &      e_opt_out(iz), a_opt_out(iz), rms_out(iz))

c        Re-entry propagation with optimal (e, A)
         a_kep(1) = smaz(1)
         a_kep(2) = e_opt_out(iz)
         a_kep(3) = incz(1)
         a_kep(4) = raanz(1)
         a_kep(5) = aopz(1)
         a_kep(6) = maz(1)
         a_kep(7) = 0.d0

c        Preserve perigee radius
         rp = smaz(1) * (1.d0 - ez(1))
         a_kep(1) = rp / (1.d0 - e_opt_out(iz))

         call oe2car(a_kep, x0, xd0, 1.d-12)
         call jd2cal(epz(1), cal0)

         bn_val = fixed_mass / (cd_drag * a_opt_out(iz))

c        Long propagation until re-entry
         period_s = 2.d0*pi*dsqrt(a_kep(1)**3/amue)
         period_d = period_s / 86400.d0
c        Propagate up to 5 years
         nrev_re = int(5.d0*365.25d0 / period_d) + 1
         if (nrev_re .gt. max_reentry) nrev_re = max_reentry

         idump_re = 1
         call propagate_ks(
     &      x0, xd0, cal0,
     &      nrev_re, 1, 1.d-10,
     &      n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &      bn_val, idrag_flag, WE_rot, EPS_f, FR_rot,
     &      CR_srp, AM_srp, IPSR, ISHAD,
     &      PSR_srp, amuS, amuM,
     &      ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &      max_reentry, idump_re,
     &      traj_jd, traj_x, traj_xd,
     &      exit_code_re)

         zone_epoch(iz) = epz(1)

         if (exit_code_re .eq. 1) then
c           Re-entry detected
            reentry_jd(iz) = traj_jd(idump_re)
            nvalid = nvalid + 1
         else
c           No re-entry within 5 years
            reentry_jd(iz) = 0.d0
         end if

c        Narrow A bounds for next zone
         a_range = (a_hi - a_lo) * 0.5d0
         a_lo = a_opt_out(iz) - a_range / 2.d0
         a_hi = a_opt_out(iz) + a_range / 2.d0
         if (a_lo .lt. 0.01d0) a_lo = 0.01d0

 50      continue
      end do

      nzones_used = nzones

c     ============================================================
c     Step 6: Compute RPE
c     ============================================================
      call compute_rpe(reentry_jd, zone_epoch, nzones,
     &   t_obs_jd, rpe_out, t_mean, t_std, nzones_max)

      if (nvalid .eq. 0) ierr = 3

      return

 900  continue
      ierr = 1
      return
      end

c     ============================================================
      subroutine compute_rpe(reentry_jd, zone_epoch, nzones,
     &   t_obs_jd, rpe, t_mean, t_std, maxz)
c
c     RPE metric — two modes:
c     Mode 1 (t_obs > 0): RPE vs observed re-entry
c     Mode 2 (t_obs = 0): RPE vs ensemble mean
c     ============================================================
      implicit double precision (a-h, o-z)
      integer nzones, maxz
      double precision reentry_jd(maxz), zone_epoch(maxz)
      double precision t_obs_jd
      double precision rpe(maxz), t_mean, t_std

      integer i, nvalid

c     Compute mean of valid predictions
      t_mean = 0.d0
      nvalid = 0
      do i = 1, nzones
         if (reentry_jd(i) .gt. 0.d0) then
            t_mean = t_mean + reentry_jd(i)
            nvalid = nvalid + 1
         end if
      end do

      if (nvalid .eq. 0) then
         t_mean = 0.d0
         t_std = 0.d0
         do i = 1, nzones
            rpe(i) = 0.d0
         end do
         return
      end if

      t_mean = t_mean / dble(nvalid)

c     Compute std deviation
      t_std = 0.d0
      do i = 1, nzones
         if (reentry_jd(i) .gt. 0.d0) then
            t_std = t_std + (reentry_jd(i) - t_mean)**2
         end if
      end do
      if (nvalid .gt. 1) then
         t_std = dsqrt(t_std / dble(nvalid - 1))
      else
         t_std = 0.d0
      end if

c     Compute RPE per zone
      do i = 1, nzones
         if (reentry_jd(i) .le. 0.d0) then
            rpe(i) = 0.d0
         else if (t_obs_jd .gt. 0.d0) then
c           Mode 1: vs observed
            denom = t_obs_jd - zone_epoch(i)
            if (dabs(denom) .gt. 0.1d0) then
               rpe(i) = (reentry_jd(i) - t_obs_jd)
     &                / denom * 100.d0
            else
               rpe(i) = 0.d0
            end if
         else
c           Mode 2: vs ensemble mean
            denom = t_mean - zone_epoch(i)
            if (dabs(denom) .gt. 0.1d0) then
               rpe(i) = (reentry_jd(i) - t_mean)
     &                / denom * 100.d0
            else
               rpe(i) = 0.d0
            end if
         end if
      end do

      return
      end

c     ============================================================
      subroutine tle_find_ma(norad_id, epoch_jd,
     &   norad_raw, epoch_yr_raw, epoch_day_raw,
     &   ama_raw, ntle_raw, maxtle, ma_out)
c     Find mean anomaly from raw TLE array matching epoch.
c     ============================================================
      implicit double precision (a-h, o-z)
      integer norad_id, ntle_raw, maxtle
      integer norad_raw(maxtle), epoch_yr_raw(maxtle)
      double precision epoch_day_raw(maxtle)
      double precision ama_raw(maxtle)
      double precision epoch_jd, ma_out

      ma_out = 0.d0
      do i = 1, ntle_raw
         if (norad_raw(i) .eq. norad_id) then
            call tle_epoch2jd(epoch_yr_raw(i),
     &           epoch_day_raw(i), djd_tmp)
            if (dabs(djd_tmp - epoch_jd) .lt. 0.01d0) then
               ma_out = ama_raw(i)
               return
            end if
         end if
      end do

      return
      end

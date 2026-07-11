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
     &   bn_min_init, bn_max_init,
     &   idrag_flag,
     &   ipopsize, maxgen, nbits_e, nbits_a,
     &   pcross, pmute, ga_seed,
     &   ngeo_deg, nsun_deg, nmoon_deg,
     &   WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD,
     &   PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   reentry_jd, e_opt_out, bn_opt_out, rms_out,
     &   zone_epoch, nzones_used,
     &   zone_status, nzones_valid,
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
c       bn_opt_out(nzones_max) : optimal ballistic number per zone
c       rms_out(nzones_max)    : GA RMS per zone
c       zone_epoch(nzones_max) : zone start epoch (JD)
c       nzones_used            : number of zones processed
c       zone_status(nzones_max): per-zone diagnostic (issue #12):
c                                0=ok, 1=skip_propfail, 2=boundary,
c                                3=nobound, 4=skip_toofewpts. Purely
c                                diagnostic -- does not change
c                                reentry_jd/e_opt_out/bn_opt_out/rpe_out,
c                                except propfail/toofewpts already
c                                skipped a zone before this existed.
c                                If more than one condition fires for a
c                                zone, whichever check runs last wins
c                                (toofewpts, then propfail, then
c                                nobound, then boundary).
c       nzones_valid            : count of zones with zone_status==0
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
      double precision bn_min_init, bn_max_init
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
      double precision bn_opt_out(nzones_max)
      double precision rms_out(nzones_max)
      double precision zone_epoch(nzones_max)
      integer nzones_used
      integer zone_status(nzones_max), nzones_valid
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
      double precision e_grid(3), bn_grid(3)
      double precision tobs(500), apobs(500)
      integer nsurf_pts, nobs, ierr_rsm
      save surfaces


c     Zone data
      double precision epz(500), haz(500), ez(500)
      double precision incz(500), raanz(500), aopz(500)
      double precision maz(500), smaz(500)

c     Propagator output for re-entry (1 point/rev, up to 20k rev)
      parameter(max_reentry = 20000)
      double precision traj_jd(max_reentry)
      double precision traj_x(3, max_reentry)
      double precision traj_xd(3, max_reentry)
      integer exit_code_re, idump_re
      save traj_jd, traj_x, traj_xd

c     Force model
      integer n_force(3)

c     Working
      double precision bn_lo, bn_hi, e_mid, de
      double precision bn_val, t_obs_jd
      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      integer iz, i, nzone, nvalid
      double precision smin, smax
      logical bn_at_bound, e_at_bound

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
c     Fitness = RMS(predicted_apogee - observed_apogee) over all TLE
c     points in the zone. Matches original genpoen1.f algorithm.
c     ============================================================
      bn_lo = bn_min_init
      bn_hi = bn_max_init
      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1
      nvalid = 0

      do iz = 1, nzones
         zone_status(iz) = 0
         bn_at_bound = .false.

c        Extract TLE data within zone
         nzone = zone_end(iz) - zone_start(iz) + 1

c        Issue #12: too few TLEs for reliable RSM. This is a fixed
c        RSM-reliability floor (3 TLEs is sufficient), independent of
c        min_zone_pts (a zone_select tuning parameter for detecting a
c        significant decay trend at selection time, not a statement
c        about RSM's own minimum viable input). Normally unreachable
c        since zone_select already enforces min_zone_pts >= 3 in
c        practice, but guarded here for when TLE quality filtering
c        (#10) can shrink a zone after selection.
         if (nzone .lt. 3) then
            zone_status(iz) = 4
            reentry_jd(iz) = 0.d0
            e_opt_out(iz)  = 0.d0
            bn_opt_out(iz) = 0.d0
            rms_out(iz)    = -1.d0
            zone_epoch(iz) = epochs(zone_start(iz))
            write(*,'(A,I2,A,I3,A)')
     &         '  [WARN] Z', iz, ': skip, too few TLEs for '//
     &         'reliable RSM: nzone=', nzone, ' < 3'
            goto 50
         end if

         do i = 1, nzone
            epz(i)   = epochs(zone_start(iz) + i - 1)
            haz(i)   = ha_tle(zone_start(iz) + i - 1)
            ez(i)    = e_tle(zone_start(iz) + i - 1)
            incz(i)  = inc_tle(zone_start(iz) + i - 1)
            raanz(i) = raan_tle(zone_start(iz) + i - 1)
            aopz(i)  = aop_tle(zone_start(iz) + i - 1)
            smaz(i)  = a_tle(zone_start(iz) + i - 1)
         end do

         do j = 1, nzone
            maz(j) = 0.d0
            call tle_find_ma(norad_id, epz(j),
     &           norad_raw, epoch_yr_raw, epoch_day_raw,
     &           ama_raw, ntle_raw, maxtle, maz(j))
         end do

         e_mid = ez(1)
         de = 0.d0
         do i = 1, nzone
            if (ez(i) .gt. e_mid + de) de = ez(i) - e_mid
            if (ez(i) .lt. e_mid - de) de = e_mid - ez(i)
         end do
         if (de .lt. 1.d-5) de = 5.d-4

         call rsm_generate(
     &      epz, haz, ez, incz, raanz, aopz, maz, smaz,
     &      nzone,
     &      e_mid - de, e_mid + de,
     &      bn_lo, bn_hi,
     &      200, 360, 1.d-10,
     &      n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &      WE_rot, EPS_f, FR_rot,
     &      CR_srp, AM_srp, IPSR, ISHAD,
     &      PSR_srp, amuS, amuM,
     &      ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &      surfaces, nsurf_pts,
     &      tobs, apobs, nobs,
     &      e_grid, bn_grid,
     &      ierr_rsm)

         if (ierr_rsm .ne. 0 .or. nsurf_pts .le. 0) then
            zone_status(iz) = 1
            reentry_jd(iz) = 0.d0
            e_opt_out(iz)  = 0.d0
            bn_opt_out(iz) = 0.d0
            rms_out(iz)    = -1.d0
            zone_epoch(iz) = epz(1)
            write(*,'(A,I2,A)')
     &         '  [WARN] Z', iz, ': skip, RSM/propagator failure '//
     &         '(propagate_ks diverged on more than 6 of 9 grid '//
     &         'points)'
            goto 50
         end if

c        Issue #12: does the RSM envelope bound the observations?
c        surfaces(i,ie,ia) is already interpolated at the same times
c        as apobs(i) (rsm_generate builds it that way), so this is a
c        direct per-observation min/max check, no interpolation needed.
c        Warn once per zone (with a count) rather than per point --
c        real zones can have dozens of observations, and one line per
c        violation would flood the console. Warn only -- still runs
c        GA on the existing surfaces.
         nobound_cnt = 0
         do i = 1, nobs
            smin = surfaces(i,1,1)
            smax = surfaces(i,1,1)
            do ie2 = 1, 3
               do ia2 = 1, 3
                  if (surfaces(i,ie2,ia2) .lt. smin)
     &               smin = surfaces(i,ie2,ia2)
                  if (surfaces(i,ie2,ia2) .gt. smax)
     &               smax = surfaces(i,ie2,ia2)
               end do
            end do
            if (apobs(i) .lt. smin .or. apobs(i) .gt. smax)
     &         nobound_cnt = nobound_cnt + 1
         end do
         if (nobound_cnt .gt. 0) then
            zone_status(iz) = 3
            write(*,'(A,I2,A,I3,A,I3,A)')
     &         '  [WARN] Z', iz, ': RSM envelope does not bound ',
     &         nobound_cnt, ' of ', nobs,
     &         ' observations -- proceeding anyway'
         end if

c        GA optimization: match predicted apogee trajectory to TLE observations
         call ga_optimize(
     &      surfaces, max_surf, e_grid, bn_grid,
     &      3, 3, nsurf_pts,
     &      tobs, apobs, nobs,
     &      e_grid(1), e_grid(3),
     &      bn_grid(1), bn_grid(3),
     &      ipopsize, maxgen, nbits_e, nbits_a,
     &      pcross, pmute, ga_seed,
     &      e_opt_out(iz), bn_opt_out(iz), rms_out(iz))

c        Issue #12: GA boundary solution (optimum within 15% of a
c        search bound suggests the true value may lie outside the
c        search range). 15% (not a tighter %) because the
c        trajectory-matching RMS fitness (RMS over all TLE points in
c        the zone, not just endpoint magnitude) does not produce a
c        purely monotonic response to BN within a narrow, badly-off
c        window -- empirically, a GA run whose true optimum is well
c        outside a 5 kg/m2 window still lands ~10-12% of the window
c        width from the edge, not exactly on it. Doesn't change this
c        zone's outputs -- only flags it and, for BN specifically,
c        widens (instead of narrows) the next zone's search range
c        below.
         e_at_bound =
     &      dabs(e_opt_out(iz)-e_grid(1)) .le.
     &         0.15d0*(e_grid(3)-e_grid(1)) .or.
     &      dabs(e_opt_out(iz)-e_grid(3)) .le.
     &         0.15d0*(e_grid(3)-e_grid(1))
         bn_at_bound =
     &      dabs(bn_opt_out(iz)-bn_grid(1)) .le.
     &         0.15d0*(bn_grid(3)-bn_grid(1)) .or.
     &      dabs(bn_opt_out(iz)-bn_grid(3)) .le.
     &         0.15d0*(bn_grid(3)-bn_grid(1))
         if (e_at_bound .or. bn_at_bound) then
            zone_status(iz) = 2
            write(*,'(A,I2,A)')
     &         '  [WARN] Z', iz, ': GA optimum at search boundary '//
     &         '-- true value may lie outside the search range'
         end if

c        Re-entry propagation with optimal (e, BN)
         a_kep(1) = smaz(1)
         a_kep(2) = e_opt_out(iz)
         a_kep(3) = incz(1)
         a_kep(4) = raanz(1)
         a_kep(5) = aopz(1)
         a_kep(6) = maz(1)
         a_kep(7) = 0.d0

         rp = smaz(1) * (1.d0 - ez(1))
         a_kep(1) = rp / (1.d0 - e_opt_out(iz))

         call oe2car(a_kep, x0, xd0, 1.d-12)
         call jd2cal(epz(1), cal0)

         bn_val = bn_opt_out(iz)
         zone_epoch(iz) = epz(1)

         if (idrag_flag .eq. 0) then
            reentry_jd(iz) = 0.d0
            exit_code_re   = 0
            idump_re       = 0
            write(*,'(A,I2,A,F8.2)')
     &         '  Z', iz, ': skip re-entry (IDRAG=0) BN=', bn_val
         else
            period_s = 2.d0*pi*dsqrt(a_kep(1)**3/amue)
            period_d = period_s / 86400.d0
            nrev_re = int(5.d0*365.25d0 / period_d) + 1
            if (nrev_re .gt. max_reentry) nrev_re = max_reentry

            idump_re = 1
            call propagate_ks(
     &         x0, xd0, cal0,
     &         nrev_re, 360, 1.d-10,
     &         n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &         bn_val, idrag_flag, WE_rot, EPS_f, FR_rot,
     &         CR_srp, AM_srp, IPSR, ISHAD,
     &         PSR_srp, amuS, amuM,
     &         ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &         max_reentry, idump_re,
     &         traj_jd, traj_x, traj_xd,
     &         exit_code_re)

            write(*,'(A,I2,A,I2,A,I6,A,F8.2)')
     &         '  Z', iz, ': exit=', exit_code_re,
     &         ' pts=', idump_re, ' BN=', bn_val

            if (exit_code_re .eq. 1) then
               reentry_jd(iz) = traj_jd(idump_re)
               nvalid = nvalid + 1
            else
               reentry_jd(iz) = 0.d0
            end if
         end if

c        Narrow BN bounds for next zone -- unless this zone's GA
c        optimum saturated the BN search bound (issue #12), in which
c        case widen instead: the true BN likely lies outside what was
c        just searched, and narrowing further would only entrench it.
         if (bn_at_bound) then
            bn_range = (bn_hi - bn_lo) * 1.5d0
         else
            bn_range = (bn_hi - bn_lo) * 0.5d0
         end if
         bn_lo = bn_opt_out(iz) - bn_range / 2.d0
         bn_hi = bn_opt_out(iz) + bn_range / 2.d0
         if (bn_lo .lt. 10.d0) bn_lo = 10.d0

 50      continue
      end do

      nzones_used = nzones

      nzones_valid = 0
      do iz = 1, nzones
         if (zone_status(iz) .eq. 0) nzones_valid = nzones_valid + 1
      end do

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

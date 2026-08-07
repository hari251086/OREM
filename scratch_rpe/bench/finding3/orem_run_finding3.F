c     ============================================================
c     orem_run_finding3.F -- issue #36 Finding 3 generalization check.
c
c     Verbatim copy of orem_run (like orem_run_timed.F -- same reason
c     for existing, same reason helper subprograms aren't duplicated
c     here since they resolve against liborem.a), with one diagnostic
c     block inserted right after ga_optimize returns e_opt_out(iz)/
c     bn_opt_out(iz), before anything else touches zone state.
c
c     At that point, using ONLY the real production variables already
c     in scope (surfaces/e_grid/bn_grid/tobs/apobs from the real
c     rsm_generate call, e_opt_out/bn_opt_out from the real GA), two
c     RMS figures are computed and appended to unit 92 (opened by the
c     caller):
c       rms_interp -- what the GA's own fitness function saw: bilinear
c         interpolation (ga_twoint, the real production routine) of
c         the 9 grid-corner surfaces at (e_opt,bn_opt), vs apobs.
c         Should match rms_out(iz)/rms_fit almost exactly -- serves as
c         a validity check on this diagnostic itself.
c       rms_direct -- a genuine propagate_ks run AT (e_opt,bn_opt)
c         (same IC-construction and rsm_tinterp time-sampling
c         rsm_generate itself uses for a grid corner, just pointed at
c         the GA's converged point instead of a grid corner), vs
c         apobs. This is what Finding 3 (issue #36) asks about: is
c         the GA's optimum actually as good as it looks, or does the
c         coarse 3x3 interpolation understate the true fit error?
c
c     No reconstruction of IC/grid bounds happens anywhere here --
c     that was the confound that inflated Finding 1's original n=1
c     numbers 3x (see project memory). Everything below is the exact
c     same variables the real GA/RSM already computed, read in place.
c     ============================================================
      subroutine orem_run_finding3(
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

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

c     issue #36 Finding 3: object label for the output CSV, set by
c     the caller before each orem_run_finding3 call.
      character*30 f3_obj_name
      common /f3_label/ f3_obj_name

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
      parameter(maxpts = 15000, max_surf = 5000)
      parameter(maxtle = 100000)

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
      save names_raw, catid_raw, norad_raw
      save epoch_yr_raw, epoch_day_raw, revnum_raw
      save ndot_raw, nddot_raw, bstar_raw
      save ainc_raw, raan_raw, ecc_raw
      save aop_raw, ama_raw, amm_raw

      integer zone_start(nzones_max), zone_end(nzones_max)
      integer nzones

      double precision surfaces(max_surf, 3, 3)
      double precision e_grid(3), bn_grid(3)
      double precision tobs(500), apobs(500)
      integer nsurf_pts, nobs, ierr_rsm
      save surfaces

      double precision epz(500), haz(500), ez(500)
      double precision incz(500), raanz(500), aopz(500)
      double precision maz(500), smaz(500)

      parameter(max_reentry = 20000)
      double precision traj_jd(max_reentry)
      double precision traj_x(3, max_reentry)
      double precision traj_xd(3, max_reentry)
      integer exit_code_re, idump_re
      save traj_jd, traj_x, traj_xd

      integer n_force(3)

      double precision bn_lo, bn_hi, e_mid, de
      double precision bn_val, t_obs_jd, bn_range, bn_center
      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      integer iz, i, nzone, nvalid, nzone1
      double precision smin, smax
      logical bn_at_bound, e_at_bound
      double precision bn_floor_est
      logical bn_floor_ok

      double precision a_kep_osc1(7)
      logical osc_ok

      double precision bn_lo_zone, bn_hi_zone
      logical bn_prior_ok

c     --- Finding 3 diagnostic working vars ---
      double precision f3_fout, f3_direct, f3_f2d(3,3)
      double precision f3_interp_apobs(500), f3_direct_apobs(500)
      double precision f3_sum_i, f3_sum_d, f3_rms_i, f3_rms_d
      double precision f3_zone_dur, f3_rp, f3_period_s, f3_period_d
      integer f3_nrev, f3_ntraj, f3_ipt, f3_it
      double precision f3_traj_ha(2000), f3_traj_t(2000)
      double precision f3_r_mag, f3_v_mag, f3_sma_pt
      double precision f3_hx, f3_hy, f3_hz, f3_h_mag, f3_p_sl
      double precision f3_ecc_pt
      integer f3_exit_code, f3_idump, f3_idrag_loc

      ierr = 0
      nzones_used = 0
      nzones_valid = 0
      t_mean = 0.d0
      t_std = 0.d0
      do i = 1, nzones_max
         zone_status(i) = 0
      end do

      call tle_evolve(tle_file, norad_id, maxpts,
     &   epochs, a_tle, e_tle, inc_tle,
     &   raan_tle, aop_tle,
     &   ha_tle, hp_tle, ls_tle,
     &   npts, ierr_tle)

      if (ierr_tle .ne. 0) then
         ierr = 1
         return
      end if

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

      open(unit=77, file=tle_file, status='old', err=900)
      call read_tle(77, maxtle, names_raw, catid_raw,
     &   norad_raw, epoch_yr_raw, epoch_day_raw,
     &   ndot_raw, nddot_raw, bstar_raw,
     &   ainc_raw, raan_raw, ecc_raw,
     &   aop_raw, ama_raw, amm_raw,
     &   revnum_raw, ntle_raw)
      close(77)

      call zone_select(epochs, ha_tle, npts,
     &   min_zone_pts, max_zone_days,
     &   r2_thresh, slope_thresh,
     &   nzones_max,
     &   zone_start, zone_end, nzones)

      if (nzones .eq. 0) then
         ierr = 2
         return
      end if

      t_obs_jd = 0.d0
      if (t_obs_cal(1) .gt. 1900.d0) then
         call cal2jd(t_obs_cal, t_obs_jd)
      end if

      bn_lo = bn_min_init
      bn_hi = bn_max_init
      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1
      nvalid = 0

      nzone1 = zone_end(1) - zone_start(1) + 1
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

      call estimate_bn_floor(epz, smaz, ez, incz, raanz, aopz, maz,
     &   nzone1, n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &   WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD, PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   bn_floor_est, bn_floor_ok)

      if (bn_floor_ok) then
         if (0.5d0*bn_floor_est .lt. bn_lo) then
            bn_lo = 0.5d0*bn_floor_est
         end if
      end if

      do iz = 1, nzones
         zone_status(iz) = 0
         bn_at_bound = .false.

         nzone = zone_end(iz) - zone_start(iz) + 1

         if (nzone .lt. 3) then
            zone_status(iz) = 4
            reentry_jd(iz) = 0.d0
            e_opt_out(iz)  = 0.d0
            bn_opt_out(iz) = 0.d0
            rms_out(iz)    = -1.d0
            zone_epoch(iz) = epochs(zone_start(iz))
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

         call tle_find_osc(norad_id, epz(1),
     &      norad_raw, epoch_yr_raw, epoch_day_raw,
     &      ainc_raw, raan_raw, ecc_raw, aop_raw, ama_raw,
     &      amm_raw, bstar_raw, ntle_raw, maxtle,
     &      a_kep_osc1, osc_ok)
         if (osc_ok) then
            smaz(1)  = a_kep_osc1(1)
            ez(1)    = a_kep_osc1(2)
            incz(1)  = a_kep_osc1(3)
            raanz(1) = a_kep_osc1(4)
            aopz(1)  = a_kep_osc1(5)
            maz(1)   = a_kep_osc1(6)
         end if

         e_mid = ez(1)
         de = 0.d0
         do i = 1, nzone
            if (ez(i) .gt. e_mid + de) de = ez(i) - e_mid
            if (ez(i) .lt. e_mid - de) de = e_mid - ez(i)
         end do
         if (de .lt. 1.d-5) de = 5.d-4

         call estimate_bn_bstar_prior(norad_id, epz, nzone,
     &      norad_raw, epoch_yr_raw, epoch_day_raw, bstar_raw,
     &      ntle_raw, maxtle,
     &      bn_lo, bn_hi, bn_lo_zone, bn_hi_zone, bn_prior_ok)

         call rsm_generate(
     &      epz, haz, ez, incz, raanz, aopz, maz, smaz,
     &      nzone,
     &      e_mid - de, e_mid + de,
     &      bn_lo_zone, bn_hi_zone,
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
            goto 50
         end if

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
         end if

         call ga_optimize(
     &      surfaces, max_surf, e_grid, bn_grid,
     &      3, 3, nsurf_pts,
     &      tobs, apobs, nobs,
     &      e_grid(1), e_grid(3),
     &      bn_grid(1), bn_grid(3),
     &      ipopsize, maxgen, nbits_e, nbits_a,
     &      pcross, pmute, ga_seed,
     &      e_opt_out(iz), bn_opt_out(iz), rms_out(iz))

c        ======================================================
c        issue #36 Finding 3: generalization check. Uses ONLY
c        the real production variables already in scope above
c        (surfaces/e_grid/bn_grid/tobs/apobs/e_opt_out/bn_opt_out)
c        -- nothing reconstructed.
c        ======================================================

c        --- rms_interp: what the GA's fitness fn actually saw ---
         f3_sum_i = 0.d0
         do f3_it = 1, nobs
            do ie2 = 1, 3
               do ia2 = 1, 3
                  f3_f2d(ie2,ia2) = surfaces(f3_it,ie2,ia2)
               end do
            end do
            call ga_twoint(e_grid, bn_grid, f3_f2d,
     &         e_opt_out(iz), bn_opt_out(iz), f3_fout, 3, 3)
            f3_interp_apobs(f3_it) = f3_fout
            f3_sum_i = f3_sum_i + (f3_fout - apobs(f3_it))**2
         end do
         f3_rms_i = dsqrt(f3_sum_i / dble(nobs))

c        --- rms_direct: real propagate_ks at (e_opt,bn_opt), same
c            IC-construction and rsm_tinterp sampling rsm_generate
c            itself uses for a grid corner ---
         a_kep(1) = smaz(1)
         a_kep(2) = e_opt_out(iz)
         a_kep(3) = incz(1)
         a_kep(4) = raanz(1)
         a_kep(5) = aopz(1)
         a_kep(6) = maz(1)
         a_kep(7) = 0.d0
         f3_rp = smaz(1) * (1.d0 - ez(1))
         a_kep(1) = f3_rp / (1.d0 - e_opt_out(iz))

         call oe2car(a_kep, x0, xd0, 1.d-12)
         call jd2cal(epz(1), cal0)

         f3_zone_dur = epz(nzone) - epz(1)
         f3_period_s = 2.d0*pi*dsqrt(a_kep(1)**3/amue)
         f3_period_d = f3_period_s / 86400.d0
         f3_nrev = int(f3_zone_dur / f3_period_d) + 2

         f3_idrag_loc = 1
         f3_idump = 1
         do f3_ipt = 1, 2000
            traj_jd(f3_ipt) = 0.d0
         end do
         call propagate_ks(
     &      x0, xd0, cal0,
     &      f3_nrev, 360, 1.d-10,
     &      n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &      bn_opt_out(iz), f3_idrag_loc, WE_rot, EPS_f, FR_rot,
     &      CR_srp, AM_srp, IPSR, ISHAD,
     &      PSR_srp, amuS, amuM,
     &      ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &      2000, f3_idump,
     &      traj_jd, traj_x, traj_xd,
     &      f3_exit_code)

         f3_ntraj = 0
         if (f3_exit_code .ne. 2) then
            do f3_ipt = 1, 2000
               if (traj_jd(f3_ipt) .le. 0.d0) goto 916
               if (traj_jd(f3_ipt) - epz(1) .gt.
     &             f3_zone_dur + 0.5d0) goto 916
               f3_ntraj = f3_ntraj + 1

               f3_r_mag = dsqrt(traj_x(1,f3_ipt)**2 +
     &                          traj_x(2,f3_ipt)**2 +
     &                          traj_x(3,f3_ipt)**2)
               f3_v_mag = dsqrt(traj_xd(1,f3_ipt)**2 +
     &                          traj_xd(2,f3_ipt)**2 +
     &                          traj_xd(3,f3_ipt)**2)
               f3_sma_pt = 1.d0 /
     &            (2.d0/f3_r_mag - f3_v_mag**2/amue)
               f3_hx = traj_x(2,f3_ipt)*traj_xd(3,f3_ipt)
     &               - traj_x(3,f3_ipt)*traj_xd(2,f3_ipt)
               f3_hy = traj_x(3,f3_ipt)*traj_xd(1,f3_ipt)
     &               - traj_x(1,f3_ipt)*traj_xd(3,f3_ipt)
               f3_hz = traj_x(1,f3_ipt)*traj_xd(2,f3_ipt)
     &               - traj_x(2,f3_ipt)*traj_xd(1,f3_ipt)
               f3_h_mag = dsqrt(f3_hx**2 + f3_hy**2 + f3_hz**2)
               f3_p_sl = f3_h_mag**2 / amue
               f3_ecc_pt = dsqrt(dabs(1.d0 - f3_p_sl/f3_sma_pt))
               f3_traj_ha(f3_ntraj) =
     &            f3_sma_pt*(1.d0+f3_ecc_pt) - R_Earth
               f3_traj_t(f3_ntraj) = traj_jd(f3_ipt)
            end do
         end if
 916     continue

         f3_sum_d = 0.d0
         do f3_it = 1, nobs
            call rsm_tinterp(f3_traj_t, f3_traj_ha,
     &         f3_ntraj, epz(1), epz(f3_it) - epz(1), f3_direct)
            f3_direct_apobs(f3_it) = f3_direct
            f3_sum_d = f3_sum_d + (f3_direct - apobs(f3_it))**2
         end do
         f3_rms_d = dsqrt(f3_sum_d / dble(nobs))

         write(92,910) f3_obj_name, iz, nzone, f3_zone_dur,
     &      e_grid(1), e_grid(3), bn_grid(1), bn_grid(3),
     &      e_opt_out(iz), bn_opt_out(iz), rms_out(iz),
     &      f3_rms_i, f3_rms_d, f3_ntraj
 910     format(A,',',I2,',',I3,',',F8.3,',',
     &      F9.6,',',F9.6,',',F9.4,',',F9.4,',',
     &      F9.6,',',F9.4,',',F10.4,',',F10.4,',',F10.4,',',I5)

c        ======================================================
c        end Finding 3 diagnostic block -- rest is unmodified
c        ======================================================

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
         end if

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

            if (exit_code_re .eq. 1) then
               reentry_jd(iz) = traj_jd(idump_re)
               nvalid = nvalid + 1
            else
               reentry_jd(iz) = 0.d0
            end if
         end if

         if ((idrag_flag .eq. 1 .and. reentry_jd(iz) .gt. 0.d0)
     &       .or.
     &       (idrag_flag .eq. 0 .and. zone_status(iz) .eq. 0)) then
            if (bn_at_bound) then
               bn_center = (bn_hi + bn_lo) / 2.d0
               bn_range = (bn_hi - bn_lo) * 1.5d0
               bn_lo = bn_center - bn_range / 2.d0
               bn_hi = bn_center + bn_range / 2.d0
            else
               bn_range = (bn_hi - bn_lo) * 0.5d0
               bn_lo = bn_opt_out(iz) - bn_range / 2.d0
               bn_hi = bn_opt_out(iz) + bn_range / 2.d0
            end if
            if (bn_lo .lt. 10.d0) bn_lo = 10.d0
         end if

 50      continue
      end do

      nzones_used = nzones

      nzones_valid = 0
      do iz = 1, nzones
         if (zone_status(iz) .eq. 0) nzones_valid = nzones_valid + 1
      end do

      call compute_rpe(reentry_jd, zone_epoch, nzones,
     &   t_obs_jd, rpe_out, t_mean, t_std, nzones_max)

      if (nvalid .eq. 0) ierr = 3

      return

 900  continue
      ierr = 1
      return
      end

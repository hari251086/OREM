c     ============================================================
c     orem_run_gasweep.F -- issue #40 option 1: GA evaluation-budget
c     floor, tested against REAL production landscapes.
c
c     Verbatim copy of orem_run (same reason orem_run_finding3.F and
c     orem_run_timed.F exist -- helper subprograms resolve against
c     liborem.a, only this file's own diagnostic block is new), with
c     one block inserted right after the PRODUCTION ga_optimize call
c     returns e_opt_out(iz)/bn_opt_out(iz)/rms_out(iz), before
c     anything else touches zone state.
c
c     Question: option 1 of issue #40 asks what the minimum GA
c     population/generation budget is that still converges reliably,
c     so a future expensive (direct-propagation) fitness function
c     knows its real affordable budget. Testing this against SYNTHETIC
c     linear surfaces (test_ga_sensitivity.F, pre-existing) cannot
c     answer it -- issue #12's pop=4 "range-invariant seed artifact"
c     only appears on real weak-gradient landscapes, and that test
c     never varies the RNG seed at all (single seed=0.123 throughout),
c     so it structurally cannot detect a seed-sensitivity failure mode.
c
c     This diagnostic instead sweeps (ipopsize, maxgen) x 10 distinct
c     seeds through ga_optimize on the REAL production surfaces/
c     e_grid/bn_grid/tobs/apobs already built above by this zone's own
c     rsm_generate call (no reconstruction, same "read production
c     variables in place" discipline as orem_run_finding3.F), and
c     measures, for each (pop,gen) combo, what fraction of the 10
c     seeds reach within 5% of this zone's own production baseline
c     rms_out(iz) (pop=20/gen=up to 500, the config actually shipped).
c
c     Downstream zone-chain state (e_opt_out/bn_opt_out/rms_out,
c     zone_status, BN carryover) is untouched by the sweep -- it reads
c     the same arrays the production call already populated and never
c     writes back to them, so this file's own zone-to-zone behavior
c     for reentry_jd/rpe stays bit-identical to plain orem_run.
c     ============================================================
      subroutine orem_run_gasweep(
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

c     issue #40: object label for the output CSV, set by the caller
c     before each orem_run_gasweep call.
      character*30 gs_obj_name
      common /gs_label/ gs_obj_name

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

c     --- GA sweep diagnostic working vars ---
      integer npoplist, ngenlist, nseedlist
      parameter(npoplist = 6, ngenlist = 7, nseedlist = 10)
      integer pop_list(npoplist), gen_list(ngenlist)
      double precision seed_list(nseedlist)
      double precision gs_e_s, gs_bn_s, gs_rms_s
      double precision gs_sw_bn(nseedlist), gs_sw_rms(nseedlist)
      double precision gs_mean_bn, gs_std_bn
      double precision gs_mean_rms, gs_min_rms, gs_max_rms
      double precision gs_tol, gs_bn_span
      integer gs_ip, gs_ig, gs_is, gs_nconv, gs_nevals
      integer gs_zone_dur_i

      data pop_list  /4, 6, 8, 10, 14, 20/
      data gen_list  /10, 20, 40, 75, 125, 200, 500/
      data seed_list /0.101d0, 0.213d0, 0.327d0, 0.441d0, 0.558d0,
     &                0.669d0, 0.777d0, 0.888d0, 0.912d0, 0.055d0/

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

c        --- PRODUCTION GA call: this is the real baseline this zone
c            actually ships (pop=ipopsize/gen=maxgen from the caller's
c            own config) -- untouched, feeds everything downstream. ---
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
c        issue #40 option 1: GA evaluation-budget sweep. Uses
c        ONLY the real production variables already in scope
c        above (surfaces/e_grid/bn_grid/tobs/apobs/nobs) -- same
c        landscape the production call just searched. Writes
c        nothing back to e_opt_out/bn_opt_out/rms_out/zone_status,
c        so production zone-chain behavior below is unaffected.
c        ======================================================
         gs_bn_span = bn_grid(3) - bn_grid(1)

         do gs_ip = 1, npoplist
            do gs_ig = 1, ngenlist

               gs_min_rms = 1.d30
               gs_max_rms = -1.d30
               gs_mean_rms = 0.d0
               gs_mean_bn = 0.d0

               do gs_is = 1, nseedlist
                  call ga_optimize(
     &               surfaces, max_surf, e_grid, bn_grid,
     &               3, 3, nsurf_pts,
     &               tobs, apobs, nobs,
     &               e_grid(1), e_grid(3),
     &               bn_grid(1), bn_grid(3),
     &               pop_list(gs_ip), gen_list(gs_ig),
     &               nbits_e, nbits_a,
     &               pcross, pmute, seed_list(gs_is),
     &               gs_e_s, gs_bn_s, gs_rms_s)

                  gs_sw_bn(gs_is)  = gs_bn_s
                  gs_sw_rms(gs_is) = gs_rms_s
                  gs_mean_rms = gs_mean_rms + gs_rms_s
                  gs_mean_bn  = gs_mean_bn + gs_bn_s
                  if (gs_rms_s .lt. gs_min_rms) gs_min_rms = gs_rms_s
                  if (gs_rms_s .gt. gs_max_rms) gs_max_rms = gs_rms_s
               end do

               gs_mean_rms = gs_mean_rms / dble(nseedlist)
               gs_mean_bn  = gs_mean_bn / dble(nseedlist)

               gs_std_bn = 0.d0
               do gs_is = 1, nseedlist
                  gs_std_bn = gs_std_bn +
     &               (gs_sw_bn(gs_is) - gs_mean_bn)**2
               end do
               gs_std_bn = dsqrt(gs_std_bn / dble(nseedlist))

c              Converged: within 5% (relative) of this zone's own
c              production baseline rms_out(iz), plus a small additive
c              floor so a near-zero baseline isn't an impossible bar.
               gs_tol = rms_out(iz)*1.05d0 + 1.d-6
               gs_nconv = 0
               do gs_is = 1, nseedlist
                  if (gs_sw_rms(gs_is) .le. gs_tol) then
                     gs_nconv = gs_nconv + 1
                  end if
               end do

               gs_nevals = pop_list(gs_ip) * (gen_list(gs_ig) + 1)

               write(92,920) gs_obj_name, iz, nzone,
     &            e_grid(3)-e_grid(1), gs_bn_span,
     &            rms_out(iz),
     &            pop_list(gs_ip), gen_list(gs_ig), gs_nevals,
     &            gs_nconv, gs_mean_rms, gs_min_rms, gs_max_rms,
     &            gs_mean_bn, gs_std_bn
 920           format(A,',',I2,',',I3,',',
     &            F10.6,',',F10.4,',',F12.6,',',
     &            I3,',',I4,',',I7,',',
     &            I3,',',F12.6,',',F12.6,',',F12.6,',',
     &            F10.4,',',F10.4)

            end do
         end do

c        ======================================================
c        end GA sweep diagnostic block -- rest is unmodified
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

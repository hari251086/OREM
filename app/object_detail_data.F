c     ============================================================
c     object_detail_data.F -- OREM-Watchlist issue #33 (object detail
c     page). Production sibling of main_orem.F: reads the SAME
c     orem.cfg format orem_wrapper.py already writes for every
c     watchlist object (write_orem_cfg), so no new cfg-generation code
c     is needed on the Python side.
c
c     Reuses OREM's own tle_evolve/tle_filter/zone_select/rsm_generate/
c     ga_optimize/propagate_ks subroutines directly -- the grid
c     construction and re-entry-propagation setup are the same recipe
c     orem.F's own zone loop uses (copied from orem.F ~line 572-577
c     and ~736-777), not re-derived. One documented simplification:
c     every zone uses the cfg's own initial bn_min/bn_max bracket
c     (the same per-object estimate orem_wrapper.bn_range_from_params
c     already computes from real object mass/area/Cd data) instead of
c     orem.F's adaptive per-zone G2/G3-narrowed bracket -- narrower
c     than a universal constant, but not identical to production's own
c     zone-by-zone narrowing. See OREM-Watchlist's
c     docs/mockups/object_detail/writeup.md open question #2.
c
c     The resonance/U-turn angle (psi = mod(RAAN+AOP-L_sun, 360)) is
c     ported directly from OREM's own
c     scratch_rpe/resonance_plots/prepare_plot_data.py (issue #56
c     lineage) -- same formula, moved into Fortran so this tool has no
c     Python dependency, rather than duplicating orbital-mechanics
c     logic in OREM-Watchlist's own Python code.
c
c     Zones are capped at MAX_ZONES_OUT (5) -- rsm_generate+ga_optimize
c     +propagate_ks is a real propagation per zone (not a DB read),
c     ~5-10s each; OREM-Watchlist's dashboard is expected to cache the
c     resulting CSVs per object (see issue #33's acceptance criteria),
c     not call this on every page view.
c
c     Usage:
c       object_detail_data.exe <orem.cfg> <output_prefix>
c     Writes:
c       <output_prefix>_evolution.csv
c       <output_prefix>_surfaces.csv
c       <output_prefix>_observed.csv
c       <output_prefix>_trajectory.csv
c     ============================================================
      program object_detail_data
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter (maxpts = 15000, maxtle = 100000, max_surf = 5000)
      parameter (max_traj = 20000, ndim_a = 291, max_zones_out = 5)
      integer norad_id
      character*200 cfg_file, tle_file, out_prefix
      character*260 path_evo, path_surf, path_obs, path_traj

c     --- cfg fields (same format as main_orem.F reads) ---
      double precision t_obs_cal(6)
      integer nzones_max, min_zone_pts
      double precision max_zone_days, r2_thresh, slope_thresh
      double precision bn_min, bn_max
      integer n_bn_sweep
      integer ipopsize, maxgen, nbits_e, nbits_a
      double precision pcross, pmute, ga_seed
      integer ngeo_deg, nsun_deg, nmoon_deg
      integer idrag_flag
      double precision WE_rot, EPS_f, FR_rot
      integer IPSR, ISHAD
      double precision CR_srp, AM_srp

c     tle_evolve output (full history)
      double precision epochs(maxpts), a_tle(maxpts), e_tle(maxpts)
      double precision inc_tle(maxpts), raan_tle(maxpts)
      double precision aop_tle(maxpts)
      double precision ha_tle(maxpts), hp_tle(maxpts), lam_tle(maxpts)
      integer npts, ierr_tle

c     tle_filter
      integer mask_tle(maxpts)
      integer nfiltered, nremoved, ngaps

c     filtered arrays
      double precision epochs_f(maxpts), a_f(maxpts), e_f(maxpts)
      double precision inc_f(maxpts), raan_f(maxpts), aop_f(maxpts)
      double precision ha_f(maxpts), hp_f(maxpts)
      integer npts_f

c     zone_select
      integer zone_start(50), zone_end(50), nzones, nzones_do

c     raw TLE (mean anomaly source, orem.F Step 2 pattern)
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

c     atmosphere
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

c     per-zone working arrays
      double precision epz(2000), haz(2000), ez(2000), incz(2000)
      double precision raanz(2000), aopz(2000), maz(2000), smaz(2000)
      integer nzone

c     rsm_generate output
      integer n_force(3)
      double precision surfaces(max_surf, 3, 3)
      double precision tobs(2000), apobs(2000)
      double precision e_grid(3), bn_grid(3)
      integer nsurf_pts, nobs, ierr_rsm

c     ga_optimize output
      double precision e_opt, bn_opt, rms_opt

c     re-entry propagation
      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      double precision traj_jd(max_traj)
      double precision traj_x(3, max_traj)
      double precision traj_xd(3, max_traj)
      integer idump_re, exit_code_re

      integer iz, i, j, iarg
      double precision e_mid, de, rp, period_d, lsun, psi
      integer nrev_re

      data n_force /1, 1, 1/

      iarg = iargc()
      if (iarg .lt. 2) then
         write(*,*) 'Usage: object_detail_data.exe <orem.cfg> '//
     &      '<output_prefix>'
         stop 1
      end if
      call getarg(1, cfg_file)
      call getarg(2, out_prefix)

      call init_constants()

c     ============================================================
c     Read config (same format as main_orem.F)
c     ============================================================
      open(unit=10, file=cfg_file, status='old', err=900)
      read(10,'(A)') tle_file
      read(10,*) norad_id
      read(10,*) t_obs_cal(1), t_obs_cal(2), t_obs_cal(3),
     &           t_obs_cal(4), t_obs_cal(5), t_obs_cal(6)
      read(10,*) nzones_max
      read(10,*) min_zone_pts, max_zone_days, r2_thresh, slope_thresh
      read(10,*) bn_min, bn_max, n_bn_sweep
      read(10,*) ipopsize, maxgen, nbits_e, nbits_a,
     &           pcross, pmute, ga_seed
      read(10,*) ngeo_deg, nsun_deg, nmoon_deg
      read(10,*) idrag_flag, WE_rot, EPS_f, FR_rot
      read(10,*) IPSR, CR_srp, AM_srp, ISHAD
      close(10)
      if (nzones_max .gt. 50) nzones_max = 50

      PSR_srp = 4.56d-6
      amuS = 1.32712440018d11
      amuM = 4.902801076d3

c     ============================================================
c     Atmosphere table (verbatim from app/main_orem.F)
c     ============================================================
      natm = ndim_a
      ALT_a(1) = 60.d0
      do i = 2, 141
         ALT_a(i) = ALT_a(i-1) + 1.d0
      end do
      do i = 142, natm
         ALT_a(i) = ALT_a(i-1) + 2.d0
      end do
      open(unit=60, file='input/ATM.DAT', status='old', err=910)
      read(60,110) (SCH_a(i), i=1,61)
      read(60,110) (SCH_a(i), i=62,natm)
      read(60,111) (DEN_a(i), i=1,61)
      read(60,111) (DEN_a(i), i=62,natm)
      close(60)
 110  format(8f10.3)
 111  format(8e10.3)
      do i = 1, natm
         ALT_a(i) = ALT_a(i) * 1000.d0
         SCH_a(i) = SCH_a(i) * 1000.d0
         DEN_a(i) = DEN_a(i) * 1.0d10
      end do

c     ============================================================
c     Full TLE evolution -> Fig 5.30-style CSV. Resonance angle
c     ported from resonance_plots/prepare_plot_data.py's L_sun_deg/
c     azimuth formula (issue #56 lineage), not tle_evolve's own
c     lambda_s_out -- that output was checked against this exact
c     object during the exemplar mockup and does not show a U-turn.
c     ============================================================
      call tle_evolve(tle_file, norad_id, maxpts,
     &   epochs, a_tle, e_tle, inc_tle, raan_tle, aop_tle,
     &   ha_tle, hp_tle, lam_tle, npts, ierr_tle)
      if (ierr_tle .ne. 0 .or. npts .lt. 2) then
         write(*,*) 'ERROR: tle_evolve failed, ierr=', ierr_tle,
     &      ' npts=', npts
         stop 1
      end if

      path_evo = trim(out_prefix)//'_evolution.csv'
      open(unit=21, file=path_evo)
      write(21,'(A)')
     &   'epoch_jd,years_since_first_tle,resonance_angle_deg,'//
     &   'apogee_km,perigee_km'
      do i = 1, npts
         lsun = dmod(280.46d0 + (360.d0/365.25d0)*
     &      (epochs(i)-2451545.d0), 360.d0)
         if (lsun .lt. 0.d0) lsun = lsun + 360.d0
         psi = dmod(raan_tle(i) + aop_tle(i) - lsun, 360.d0)
         if (psi .lt. 0.d0) psi = psi + 360.d0
         write(21,'(F14.5,A,F10.5,A,F10.4,A,F10.3,A,F10.3)')
     &      epochs(i), ',', (epochs(i)-epochs(1))/365.25d0, ',',
     &      psi, ',', ha_tle(i), ',', hp_tle(i)
      end do
      close(21)

c     ============================================================
c     Filter + zone selection (cfg's own parameters)
c     ============================================================
      call tle_filter(epochs, ha_tle, e_tle, npts,
     &   4.0d0, 6.0d0, 30.0d0, 20,
     &   mask_tle, nfiltered, nremoved, ngaps)
      j = 0
      do i = 1, npts
         if (mask_tle(i) .eq. 1) then
            j = j + 1
            epochs_f(j) = epochs(i)
            a_f(j)      = a_tle(i)
            e_f(j)      = e_tle(i)
            inc_f(j)    = inc_tle(i)
            raan_f(j)   = raan_tle(i)
            aop_f(j)    = aop_tle(i)
            ha_f(j)     = ha_tle(i)
            hp_f(j)     = hp_tle(i)
         end if
      end do
      npts_f = j

      call zone_select(epochs_f, ha_f, npts_f,
     &   min_zone_pts, max_zone_days, r2_thresh, slope_thresh,
     &   nzones_max, zone_start, zone_end, nzones)
      nzones_do = nzones
      if (nzones_do .gt. max_zones_out) nzones_do = max_zones_out
      if (nzones_do .le. 0) then
         write(*,*) 'No zones found for NORAD', norad_id
         stop 0
      end if

c     ============================================================
c     Raw TLE re-read for mean anomaly (orem.F Step 2 pattern)
c     ============================================================
      open(unit=77, file=tle_file, status='old', err=900)
      call read_tle(77, maxtle, names_raw, catid_raw,
     &   norad_raw, epoch_yr_raw, epoch_day_raw,
     &   ndot_raw, nddot_raw, bstar_raw,
     &   ainc_raw, raan_raw, ecc_raw,
     &   aop_raw, ama_raw, amm_raw,
     &   revnum_raw, ntle_raw)
      close(77)

      path_surf = trim(out_prefix)//'_surfaces.csv'
      path_obs  = trim(out_prefix)//'_observed.csv'
      path_traj = trim(out_prefix)//'_trajectory.csv'
      open(unit=22, file=path_surf)
      write(22,'(A)') 'zone,ecc_level,bn_level,t_days,apogee_km'
      open(unit=23, file=path_obs)
      write(23,'(A)') 'zone,t_days,apogee_km'
      open(unit=24, file=path_traj)
      write(24,'(A)') 'zone,epoch_jd,apogee_km,e_opt,bn_opt'

      do iz = 1, nzones_do
         nzone = zone_end(iz) - zone_start(iz) + 1
         if (nzone .gt. 2000) nzone = 2000
         do i = 1, nzone
            j = zone_start(iz) + i - 1
            epz(i)   = epochs_f(j)
            haz(i)   = ha_f(j)
            ez(i)    = e_f(j)
            incz(i)  = inc_f(j)
            raanz(i) = raan_f(j)
            aopz(i)  = aop_f(j)
            smaz(i)  = a_f(j)
            call tle_find_ma(norad_id, epz(i),
     &         norad_raw, epoch_yr_raw, epoch_day_raw,
     &         ama_raw, ntle_raw, maxtle, maz(i))
         end do

c        e_mid +- de : verbatim from orem.F's own zone-grid construction
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
     &      bn_min, bn_max,
     &      200, 360, 1.d-10,
     &      n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &      WE_rot, EPS_f, FR_rot,
     &      CR_srp, AM_srp, IPSR, ISHAD,
     &      PSR_srp, amuS, amuM,
     &      ALT_a, DEN_a, SCH_a, natm,
     &      surfaces, nsurf_pts,
     &      tobs, apobs, nobs,
     &      e_grid, bn_grid,
     &      ierr_rsm)
         if (ierr_rsm .ne. 0 .or. nsurf_pts .le. 0) then
            write(*,*) 'zone', iz, ': RSM/propagator failure, skipped'
            goto 50
         end if

         do i = 1, nsurf_pts
            do j1 = 1, 3
               do j2 = 1, 3
                  write(22,'(I2,A,I1,A,I1,A,F10.4,A,F10.3)')
     &               iz, ',', j1, ',', j2, ',',
     &               dble(i-1)*(epz(nzone)-epz(1))/dble(nsurf_pts-1),
     &               ',', surfaces(i,j1,j2)
               end do
            end do
         end do
         do i = 1, nobs
            write(23,'(I2,A,F10.4,A,F10.3)')
     &         iz, ',', tobs(i), ',', apobs(i)
         end do

         call ga_optimize(
     &      surfaces, max_surf, e_grid, bn_grid,
     &      3, 3, nsurf_pts,
     &      tobs, apobs, nobs,
     &      e_grid(1), e_grid(3), bn_grid(1), bn_grid(3),
     &      ipopsize, maxgen, nbits_e, nbits_a,
     &      pcross, pmute, ga_seed,
     &      e_opt, bn_opt, rms_opt)

c        Re-entry propagation with optimal (e, BN) -- verbatim from
c        orem.F's own a_kep/x0/xd0 construction (lines ~736-751).
         a_kep(1) = smaz(1)
         a_kep(2) = e_opt
         a_kep(3) = incz(1)
         a_kep(4) = raanz(1)
         a_kep(5) = aopz(1)
         a_kep(6) = maz(1)
         a_kep(7) = 0.d0
         rp = smaz(1) * (1.d0 - ez(1))
         a_kep(1) = rp / (1.d0 - e_opt)
         call oe2car(a_kep, x0, xd0, 1.d-12)
         call jd2cal(epz(1), cal0)

         period_d = 2.d0*pi*dsqrt(a_kep(1)**3/amue) / 86400.d0
         nrev_re = int(5.d0*365.25d0 / period_d) + 1
         if (nrev_re .gt. max_traj) nrev_re = max_traj

         idump_re = 1
         do i = 1, max_traj
            traj_jd(i) = 0.d0
         end do
         call propagate_ks(
     &      x0, xd0, cal0,
     &      nrev_re, 360, 1.d-10,
     &      n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &      bn_opt, idrag_flag, WE_rot, EPS_f, FR_rot,
     &      CR_srp, AM_srp, IPSR, ISHAD,
     &      PSR_srp, amuS, amuM,
     &      ALT_a, DEN_a, SCH_a, natm,
     &      max_traj, idump_re,
     &      traj_jd, traj_x, traj_xd,
     &      exit_code_re)

c        Osculating apogee from Cartesian state -- same method as
c        scratch_rpe/diag_osc_vs_mean_apogee.F.
         do i = 1, idump_re
            if (traj_jd(i) .le. 0.d0) goto 40
            r_mag = dsqrt(traj_x(1,i)**2+traj_x(2,i)**2
     &                   +traj_x(3,i)**2)
            v_mag = dsqrt(traj_xd(1,i)**2+traj_xd(2,i)**2
     &                   +traj_xd(3,i)**2)
            sma_pt = 1.d0/(2.d0/r_mag - v_mag**2/amue)
            hx = traj_x(2,i)*traj_xd(3,i) - traj_x(3,i)*traj_xd(2,i)
            hy = traj_x(3,i)*traj_xd(1,i) - traj_x(1,i)*traj_xd(3,i)
            hz = traj_x(1,i)*traj_xd(2,i) - traj_x(2,i)*traj_xd(1,i)
            h_mag = dsqrt(hx**2+hy**2+hz**2)
            p_sl = h_mag**2/amue
            ecc_pt = dsqrt(dabs(1.d0 - p_sl/sma_pt))
            ha_pt = sma_pt*(1.d0+ecc_pt) - R_Earth
            write(24,'(I2,A,F14.5,A,F10.3,A,F9.6,A,F8.2)')
     &         iz, ',', traj_jd(i), ',', ha_pt, ',',
     &         e_opt, ',', bn_opt
         end do
 40      continue

 50      continue
      end do

      close(22)
      close(23)
      close(24)

      write(*,'(A,I3,A)') 'DONE (', nzones_do, ' zones)'
      stop

 900  write(*,*) 'ERROR opening TLE/config file'
      stop 1
 910  write(*,*) 'ERROR opening input/ATM.DAT'
      stop 1
      end

c     ============================================================
c     tle_find_ma -- verbatim copy of orem.F's own helper (src/orem.F
c     ~line 1069). Not linked from orem.F directly to avoid pulling in
c     orem_run's whole dependency graph for one small lookup helper.
c     ============================================================
      subroutine tle_find_ma(norad_id, epoch_jd,
     &   norad_raw, epoch_yr_raw, epoch_day_raw,
     &   ama_raw, ntle_raw, maxtle, ma_out)
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

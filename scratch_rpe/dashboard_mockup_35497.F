c     ============================================================
c     dashboard_mockup_35497.F -- OREM-Watchlist issue (object detail
c     page mockup), thesis Section 5.4.2 / Figures 5.30-5.38 style.
c
c     Read-only diagnostic. No production code touched -- reuses
c     orem.F/rsm.F/ga.F/propagate_ks.F's own subroutines exactly as
c     orem_run itself calls them (same recipe, copied from orem.F's
c     own zone loop: rsm_generate's e_mid+-de grid construction at
c     ~line 572-577, the a_kep/x0/xd0 re-entry-propagation setup at
c     ~line 736-751, and diag_osc_vs_mean_apogee.F's osculating-apogee-
c     from-Cartesian-state method) rather than re-deriving any of the
c     physics. The only thing NOT reused from orem.F is its adaptive
c     G2/G3 BN-prior narrowing across zones (orem.F lines ~380-924) --
c     this uses a single fixed bn_lo/bn_hi bracket for every zone
c     instead, which is a real simplification, documented in the
c     dashboard mockup writeup, not a claim of matching production
c     exactly.
c
c     Dumps three CSVs for NORAD 35497 (thesis case (ii), Ariane 5
c     ESC-A upper stage -- chosen because it's one of the thesis's own
c     three worked examples AND already a first-class OREM test/
c     resonance-campaign object, decayed 2016-10-31 so it has real
c     ground truth):
c       scratch_rpe/dashboard_mockup_35497_evolution.csv   -- Fig 5.30
c       scratch_rpe/dashboard_mockup_35497_surfaces.csv    -- Fig 5.33
c       scratch_rpe/dashboard_mockup_35497_observed.csv    -- Fig 5.33 dots
c       scratch_rpe/dashboard_mockup_35497_trajectory.csv  -- Fig 5.36
c
c     Compile:
c       (vcvars64 + Intel vars.bat, see GitHub memory
c        reference_ksrop_compile.md)
c       ifx /heap-arrays /F:16777216
c           scratch_rpe\dashboard_mockup_35497.F
c           src\tle_evolution.F src\tle_filter.F src\zone_select.F
c           src\rsm.F src\ga.F src\propagate_ks.F
c           build\dependencies\ksrop\src\Subrouts.F
c           build\dependencies\ksrop\src\Legendre.F
c           build\dependencies\ksrop\src\LegendreTess.F
c           build\dependencies\ksrop\src\TLEread.F
c           build\dependencies\ksrop\src\swx.F
c           build\dependencies\ksrop\src\jr71_profile.F
c           build\dependencies\ksrop\src\DragOblateCorotating.F
c           /exe:scratch_rpe\dashboard_mockup_35497.exe
c     ============================================================
      program dashboard_mockup_35497
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter (maxpts = 15000, maxtle = 100000, max_surf = 5000)
      parameter (max_traj = 20000, ndim_a = 291)
      integer norad_id
      character*80 tle_file

c     tle_evolve output (full history)
      double precision epochs(maxpts), a_tle(maxpts), e_tle(maxpts)
      double precision inc_tle(maxpts), raan_tle(maxpts)
      double precision aop_tle(maxpts)
      double precision ha_tle(maxpts), hp_tle(maxpts), lam_tle(maxpts)
      integer npts, ierr_tle

c     tle_filter
      integer mask_tle(maxpts)
      integer nfiltered, nremoved, ngaps

c     filtered arrays (compacted in place, same pattern as the
c     existing scratch_rpe/diag_*.F scripts)
      double precision epochs_f(maxpts), a_f(maxpts), e_f(maxpts)
      double precision inc_f(maxpts), raan_f(maxpts), aop_f(maxpts)
      double precision ha_f(maxpts), hp_f(maxpts)
      integer npts_f

c     zone_select
      integer zone_start(8), zone_end(8), nzones

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
      double precision epz(500), haz(500), ez(500), incz(500)
      double precision raanz(500), aopz(500), maz(500), smaz(500)
      integer nzone

c     rsm_generate output
      integer n_force(3)
      double precision surfaces(max_surf, 3, 3)
      double precision tobs(500), apobs(500)
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

      integer iz, i, j
      double precision e_mid, de, rp, period_d
      integer nrev_re

c     force model / SRP defaults -- matched to input/orem_33587.cfg
c     and input/orem_40943.cfg (real, validated production configs),
c     not invented for this mockup.
      data n_force /1, 1, 1/
      parameter (ngeo_deg = 20, nsun_deg = 2, nmoon_deg = 3)
      parameter (idrag_flag = 1)
      parameter (WE_rot = 7.2921150d-5, EPS_f = 3.35281066d-3)
      parameter (FR_rot = 1.0d0)
      parameter (IPSR = 1, CR_srp = 1.2d0, AM_srp = 0.01d0, ISHAD = 2)
      parameter (PSR_srp = 4.56d-6)
      parameter (amuS = 1.32712440018d11, amuM = 4.902801076d3)

c     Simplified, fixed grid bounds -- see header comment. Table 5.4
c     of the thesis brackets case (ii)'s B around 80-120 kg/m^2; 60-140
c     is a slightly wider fixed bracket used for every zone here
c     instead of orem.F's per-zone G2/G3-narrowed bn_lo_zone/bn_hi_zone.
      parameter (bn_lo_fixed = 60.0d0, bn_hi_fixed = 140.0d0)

      norad_id = 35497
      tle_file = 'input/example_35497.tle.txt'

      call init_constants()

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
      open(unit=60, file='input/ATM.DAT', status='old', err=900)
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
c     Full TLE evolution -> Fig 5.30-style CSV (Sun-azimuth/apogee/
c     perigee vs. years since first TLE)
c     ============================================================
      call tle_evolve(tle_file, norad_id, maxpts,
     &   epochs, a_tle, e_tle, inc_tle, raan_tle, aop_tle,
     &   ha_tle, hp_tle, lam_tle, npts, ierr_tle)
      write(*,'(A,I3,A,I6)') ' tle_evolve ierr=', ierr_tle,
     &   '  npts=', npts
      if (ierr_tle .ne. 0) stop 1

      open(unit=21,
     &   file='scratch_rpe/dashboard_mockup_35497_evolution.csv')
      write(21,'(A)')
     &   'epoch_jd,years_since_first_tle,sun_azimuth_deg,'//
     &   'apogee_km,perigee_km'
      do i = 1, npts
         write(21,'(F14.5,A,F10.5,A,F10.4,A,F10.3,A,F10.3)')
     &      epochs(i), ',', (epochs(i)-epochs(1))/365.25d0, ',',
     &      lam_tle(i), ',', ha_tle(i), ',', hp_tle(i)
      end do
      close(21)

c     ============================================================
c     Filter (same params as production / test_tle_filter.F)
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
      write(*,'(A,I6,A,I6)') ' filtered npts=', npts_f,
     &   '  removed=', nremoved

c     ============================================================
c     Zone selection (same params as production write_orem_cfg)
c     ============================================================
      call zone_select(epochs_f, ha_f, npts_f,
     &   8, 30.0d0, 0.90d0, -1.0d0, 8,
     &   zone_start, zone_end, nzones)
      write(*,'(A,I3)') ' nzones=', nzones
      if (nzones .gt. 3) nzones = 3

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

      open(unit=22,
     &   file='scratch_rpe/dashboard_mockup_35497_surfaces.csv')
      write(22,'(A)') 'zone,series,t_days,apogee_km'
      open(unit=23,
     &   file='scratch_rpe/dashboard_mockup_35497_observed.csv')
      write(23,'(A)') 'zone,t_days,apogee_km'
      open(unit=24,
     &   file='scratch_rpe/dashboard_mockup_35497_trajectory.csv')
      write(24,'(A)') 'zone,epoch_jd,apogee_km,e_opt,bn_opt'

      do iz = 1, nzones
         nzone = zone_end(iz) - zone_start(iz) + 1
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
     &      bn_lo_fixed, bn_hi_fixed,
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
         write(*,'(A,I2,A,I3,A,I6,A,I4)') ' zone', iz,
     &      ': ierr_rsm=', ierr_rsm, ' nsurf_pts=', nsurf_pts,
     &      ' nobs=', nobs
         if (ierr_rsm .ne. 0 .or. nsurf_pts .le. 0) goto 50

         do i = 1, nsurf_pts
            do j1 = 1, 3
               do j2 = 1, 3
                  write(22,'(I2,A,A2,I1,A2,I1,A,F10.4,A,F10.3)')
     &               iz, ',', 'e', j1, 'b', j2, ',',
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
     &      20, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &      e_opt, bn_opt, rms_opt)
         write(*,'(A,I2,A,F9.6,A,F8.2,A,F10.4)') ' zone', iz,
     &      ': e_opt=', e_opt, ' bn_opt=', bn_opt,
     &      ' rms=', rms_opt

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
         write(*,'(A,I2,A,I2,A,I6)') ' zone', iz,
     &      ': exit=', exit_code_re, ' pts=', idump_re

c        Osculating apogee from Cartesian state -- verbatim method
c        from scratch_rpe/diag_osc_vs_mean_apogee.F.
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

      write(*,*) 'DONE'
      stop

 900  write(*,*) 'ERROR opening file'
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

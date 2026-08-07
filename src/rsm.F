c     ============================================================
c     rsm.F — Response Surface Methodology for OREM
c
c     Generates 9 mean-apogee surfaces per zone by calling
c     propagate_ks with a 3×3 grid of (eccentricity, area).
c     Extracts observed TLE mean apogee at matching time epochs
c     for GA fitness evaluation.
c
c     Author: Harishkumar Sellamuthu
c     ============================================================

c     ============================================================
      subroutine rsm_generate(
     &   epochs_zone, ha_zone, e_zone, inc_zone,
     &   raan_zone, aop_zone, ma_zone, a_zone,
     &   nzone,
     &   e_min, e_max, bn_min, bn_max,
     &   nrev_prop, istep_prop, tole_prop,
     &   n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &   WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD,
     &   PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   surfaces, nsurf_pts,
     &   tobs, apobs, nobs,
     &   e_grid, bn_grid,
     &   ierr)
c
c     Input
c       epochs_zone(nzone)  : Julian dates of zone TLEs
c       ha_zone(nzone)      : mean apogee altitude (km)
c       e_zone(nzone)       : eccentricities from TLE
c       inc/raan/aop/ma/a_zone : orbital elements (deg/km)
c       nzone               : number of TLEs in zone
c       e_min, e_max        : eccentricity grid bounds
c       a_min, a_max        : cross-sectional area bounds (m²)
c       fixed_mass          : spacecraft mass (kg)
c       cd_drag             : drag coefficient
c       nrev_prop, etc.     : propagator configuration
c       ALT/DEN/SCH_atm     : atmosphere tables
c
c     Output
c       surfaces(nsurf_pts, 3, 3) : 9 mean-apogee time-series
c       nsurf_pts           : time points per surface
c       tobs(nobs)          : observation times (days from zone start)
c       apobs(nobs)         : observed TLE mean apogee (km)
c       nobs                : number of observations
c       e_grid(3)           : eccentricity grid values
c       a_grid(3)           : area grid values
c       ierr                : 0=ok, 1=insufficient zone data,
c                             2=propagator failure
c     ============================================================
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      integer nzone
      double precision epochs_zone(nzone), ha_zone(nzone)
      double precision e_zone(nzone), inc_zone(nzone)
      double precision raan_zone(nzone), aop_zone(nzone)
      double precision ma_zone(nzone), a_zone(nzone)
      double precision e_min, e_max, bn_min, bn_max
      integer nrev_prop, istep_prop
      double precision tole_prop
      integer n_force(3), ngeo_deg, nsun_deg, nmoon_deg
      double precision WE_rot, EPS_f, FR_rot
      double precision CR_srp, AM_srp
      integer IPSR, ISHAD
      double precision PSR_srp, amuS, amuM
      integer ndim_atm
      double precision ALT_atm(ndim_atm)
      double precision DEN_atm(ndim_atm)
      double precision SCH_atm(ndim_atm)

c     Output
      parameter(max_surf = 5000)
      integer nsurf_pts, nobs, ierr
      double precision surfaces(max_surf, 3, 3)
      double precision tobs(nzone), apobs(nzone)
      double precision e_grid(3), bn_grid(3)

c     Propagator output buffers
      parameter(max_traj = 2000)
      double precision traj_jd(max_traj)
      double precision traj_x(3, max_traj)
      double precision traj_xd(3, max_traj)
      save traj_jd, traj_x, traj_xd
      integer exit_code

c     Working
      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      double precision bn_val, zone_dur
      double precision epoch0_jd, t_want
      integer ie, ia, ipt, nfail, ntraj
      integer idump_loc, idrag_loc
      double precision traj_ha(max_traj), traj_t(max_traj)
      save traj_ha, traj_t

      ierr = 0
      nsurf_pts = 0
      nobs = 0
      nfail = 0

      if (nzone .lt. 2) then
         ierr = 1
         return
      end if

c     --- Build 3×3 grid ---
      e_grid(1) = e_min
      e_grid(2) = (e_min + e_max) / 2.d0
      e_grid(3) = e_max
      bn_grid(1) = bn_min
      bn_grid(2) = (bn_min + bn_max) / 2.d0
      bn_grid(3) = bn_max

c     --- Extract observed TLE apogee (tobs relative to zone start) ---
      epoch0_jd = epochs_zone(1)
      nobs = nzone
      do i = 1, nzone
         tobs(i) = epochs_zone(i) - epoch0_jd
         apobs(i) = ha_zone(i)
      end do

c     --- Zone duration for propagation ---
      zone_dur = epochs_zone(nzone) - epoch0_jd

c     --- Convert zone start TLE epoch to calendar ---
      call jd2cal(epoch0_jd, cal0)

c     --- Base orbital elements from first TLE in zone ---
c     a_kep = [a, e, i, RAAN, AOP, MA] (angles in degrees)

c     --- Generate 9 surfaces ---
      do ie = 1, 3
         do ia = 1, 3

c           Set up orbital elements with perturbed eccentricity
            a_kep(1) = a_zone(1)
            a_kep(2) = e_grid(ie)
            a_kep(3) = inc_zone(1)
            a_kep(4) = raan_zone(1)
            a_kep(5) = aop_zone(1)
            a_kep(6) = ma_zone(1)
            a_kep(7) = 0.d0

c           Recalculate SMA for modified eccentricity
c           keeping perigee radius constant: rp = a*(1-e)
c           a_new = rp / (1 - e_new)
            rp = a_zone(1) * (1.d0 - e_zone(1))
            a_kep(1) = rp / (1.d0 - e_grid(ie))

c           Convert to Cartesian
            call oe2car(a_kep, x0, xd0, 1.d-12)

c           Ballistic number directly from grid
            bn_val = bn_grid(ia)

c           Estimate number of revolutions from zone duration
c           period = 2*pi*sqrt(a^3/mu)
            period_s = 2.d0*pi*dsqrt(a_kep(1)**3/amue)
            period_d = period_s / 86400.d0
            nrev_zone = int(zone_dur / period_d) + 2

c           Call propagator
c           Issue #12 follow-up: traj_jd/traj_x/traj_xd are SAVEd
c           across grid points and across separate rsm_generate
c           calls (different zones/objects). propagate_ks does not
c           clear indices beyond its own output, so a shorter run
c           here can otherwise inherit a longer prior run's leftover
c           trajectory tail -- and the end-of-data scan below only
c           excludes stale points *later* than the current zone, not
c           earlier ones. Zeroing traj_jd first makes every stale
c           slot fail the ".le. 0.d0" sentinel check regardless of
c           its old epoch.
            idrag_loc = 1
            idump_loc = 1
            do izero = 1, max_traj
               traj_jd(izero) = 0.d0
            end do
            call propagate_ks(
     &         x0, xd0, cal0,
     &         nrev_zone, istep_prop, tole_prop,
     &         n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &         bn_val, idrag_loc, WE_rot, EPS_f, FR_rot,
     &         CR_srp, AM_srp, IPSR, ISHAD,
     &         PSR_srp, amuS, amuM,
     &         ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &         max_traj, idump_loc,
     &         traj_jd, traj_x, traj_xd,
     &         exit_code)

c           Check for propagator failure
            if (exit_code .eq. 2) then
               nfail = nfail + 1
               if (nfail .gt. 6) then
                  ierr = 2
                  return
               end if
c              Fill surface with zeros for failed run
               do ipt = 1, max_surf
                  surfaces(ipt, ie, ia) = 0.d0
               end do
               goto 20
            end if

c           Compute apogee at each trajectory point within zone
c           Store in temporary arrays for interpolation
            ntraj = 0
            do ipt = 1, max_traj
               if (traj_jd(ipt) .le. 0.d0) goto 16
               if (traj_jd(ipt) - epoch0_jd .gt.
     &             zone_dur + 0.5d0) goto 16
               ntraj = ntraj + 1

               r_mag = dsqrt(traj_x(1,ipt)**2 +
     &                       traj_x(2,ipt)**2 +
     &                       traj_x(3,ipt)**2)
               v_mag = dsqrt(traj_xd(1,ipt)**2 +
     &                       traj_xd(2,ipt)**2 +
     &                       traj_xd(3,ipt)**2)

               sma_pt = 1.d0 /
     &           (2.d0/r_mag - v_mag**2/amue)

               hx = traj_x(2,ipt)*traj_xd(3,ipt)
     &            - traj_x(3,ipt)*traj_xd(2,ipt)
               hy = traj_x(3,ipt)*traj_xd(1,ipt)
     &            - traj_x(1,ipt)*traj_xd(3,ipt)
               hz = traj_x(1,ipt)*traj_xd(2,ipt)
     &            - traj_x(2,ipt)*traj_xd(1,ipt)
               h_mag = dsqrt(hx**2 + hy**2 + hz**2)
               p_sl = h_mag**2 / amue
               ecc_pt = dsqrt(dabs(1.d0 - p_sl/sma_pt))

               traj_ha(ntraj) = sma_pt*(1.d0+ecc_pt)-R_Earth
               traj_t(ntraj) = traj_jd(ipt)
            end do
 16         continue

c           Interpolate surface at observation times (zone only)
            nsurf_pts = nobs
            do ipt = 1, nobs
               t_want = epochs_zone(ipt) - epoch0_jd
               call rsm_tinterp(traj_t, traj_ha,
     &              ntraj, epoch0_jd, t_want,
     &              surfaces(ipt, ie, ia))
            end do

 20         continue
         end do
      end do

c     Normalize nsurf_pts to consistent value
c     (use the count from last successful run)

      return
      end

c     ============================================================
      subroutine jd2cal(jd, cal)
c     Convert Julian Date to calendar [yr, mo, day, hr, min, sec]
c     Algorithm from Meeus (1991)
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision jd, cal(6)
      integer ia, ib, ic, id, ie

      z = dint(jd + 0.5d0)
      f = jd + 0.5d0 - z

      if (z .lt. 2299161.d0) then
         ia = int(z)
      else
         alpha = dint((z - 1867216.25d0) / 36524.25d0)
         ia = int(z + 1.d0 + alpha - dint(alpha/4.d0))
      end if

      ib = ia + 1524
      ic = int((dble(ib) - 122.1d0) / 365.25d0)
      id = int(365.25d0 * dble(ic))
      ie = int(dble(ib - id) / 30.6001d0)

      day_frac = dble(ib - id) - dint(30.6001d0*dble(ie)) + f

      if (ie .lt. 14) then
         imo = ie - 1
      else
         imo = ie - 13
      end if

      if (imo .gt. 2) then
         iyr = ic - 4716
      else
         iyr = ic - 4715
      end if

      cal(1) = dble(iyr)
      cal(2) = dble(imo)
      cal(3) = dint(day_frac)
      frac_day = day_frac - cal(3)
      cal(4) = dint(frac_day * 24.d0)
      cal(5) = dint((frac_day*24.d0 - cal(4)) * 60.d0)
      cal(6) = (frac_day*86400.d0 - cal(4)*3600.d0
     &         - cal(5)*60.d0)

      return
      end

c     ============================================================
      subroutine rsm_tinterp(t_arr, y_arr, n,
     &                        t0, t_want, y_out)
c     Linear interpolation in JD-time arrays.
c     t_arr(n): JD times, y_arr(n): values.
c     t_want: days from t0. Interpolates at t0 + t_want.
c     ============================================================
      implicit double precision (a-h, o-z)
      integer n
      double precision t_arr(n), y_arr(n)
      double precision t0, t_want, y_out

      double precision t_target
      integer i

      t_target = t0 + t_want
      y_out = 0.d0

      if (n .le. 0) return
      if (n .eq. 1) then
         y_out = y_arr(1)
         return
      end if

c     Clamp to range
      if (t_target .le. t_arr(1)) then
         y_out = y_arr(1)
         return
      end if
      if (t_target .ge. t_arr(n)) then
         y_out = y_arr(n)
         return
      end if

c     Find bracketing interval
      do i = 2, n
         if (t_target .le. t_arr(i)) then
            frac = (t_target - t_arr(i-1)) /
     &             (t_arr(i) - t_arr(i-1))
            y_out = y_arr(i-1) + frac *
     &              (y_arr(i) - y_arr(i-1))
            return
         end if
      end do

      y_out = y_arr(n)
      return
      end

c     ============================================================
      subroutine rsm_direct_eval(
     &   a_zone1, e_zone1, inc_zone1, raan_zone1, aop_zone1, ma_zone1,
     &   epoch0_jd, zone_dur,
     &   e_cand, bn_cand,
     &   tobs, apobs, nobs,
     &   istep_prop, tole_prop,
     &   n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &   WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD,
     &   PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   rms_out, ierr)
c
c     Issue #40: direct-propagation fitness for a single (e, BN)
c     candidate -- same physics as one rsm_generate grid-corner call
c     (propagate over the zone's own span, extract apogee, interpolate
c     to the zone's real TLE observation times), but evaluated exactly
c     at the candidate point instead of interpolated between 9
c     pre-computed grid corners. Used by rsm_refine to locally search
c     near the GA's (interpolated-fitness) optimum using the true,
c     un-interpolated physics -- issue #36 found the interpolated
c     fitness is systematically ~1.8x more optimistic than this.
c
c     Input
c       a/e/inc/raan/aop/ma_zone1 : zone's base orbital elements
c                                    (index-1 TLE, same convention as
c                                    rsm_generate)
c       epoch0_jd, zone_dur       : zone start JD, span in days
c       e_cand, bn_cand           : candidate to evaluate
c       tobs(nobs), apobs(nobs)   : real TLE observations (days from
c                                    epoch0_jd, mean apogee km)
c       istep_prop, tole_prop     : propagator step/tolerance
c       force model / atmosphere params : as rsm_generate
c     Output
c       rms_out : RMS(predicted apogee - apobs) at this candidate
c       ierr    : 0=ok, 2=propagator failure
c     ============================================================
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision a_zone1, e_zone1, inc_zone1
      double precision raan_zone1, aop_zone1, ma_zone1
      double precision epoch0_jd, zone_dur
      double precision e_cand, bn_cand
      integer nobs
      double precision tobs(nobs), apobs(nobs)
      integer istep_prop
      double precision tole_prop
      integer n_force(3), ngeo_deg, nsun_deg, nmoon_deg
      double precision WE_rot, EPS_f, FR_rot
      double precision CR_srp, AM_srp
      integer IPSR, ISHAD
      double precision PSR_srp, amuS, amuM
      integer ndim_atm
      double precision ALT_atm(ndim_atm)
      double precision DEN_atm(ndim_atm)
      double precision SCH_atm(ndim_atm)
      double precision rms_out
      integer ierr

      parameter(max_traj = 2000)
      double precision traj_jd(max_traj)
      double precision traj_x(3, max_traj)
      double precision traj_xd(3, max_traj)
      save traj_jd, traj_x, traj_xd
      integer exit_code

      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      double precision rp, period_s, period_d
      integer nrev_zone, idrag_loc, idump_loc, izero, ipt, ntraj
      double precision traj_ha(max_traj), traj_t(max_traj)
      save traj_ha, traj_t
      double precision r_mag, v_mag, sma_pt
      double precision hx, hy, hz, h_mag, p_sl, ecc_pt
      double precision pred_ap, ssq
      integer i

      ierr = 0
      rms_out = 0.d0

      a_kep(1) = a_zone1
      a_kep(2) = e_cand
      a_kep(3) = inc_zone1
      a_kep(4) = raan_zone1
      a_kep(5) = aop_zone1
      a_kep(6) = ma_zone1
      a_kep(7) = 0.d0

      rp = a_zone1 * (1.d0 - e_zone1)
      a_kep(1) = rp / (1.d0 - e_cand)

      call oe2car(a_kep, x0, xd0, 1.d-12)
      call jd2cal(epoch0_jd, cal0)

      period_s = 2.d0*pi*dsqrt(a_kep(1)**3/amue)
      period_d = period_s / 86400.d0
      nrev_zone = int(zone_dur / period_d) + 2

      idrag_loc = 1
      idump_loc = 1
      do izero = 1, max_traj
         traj_jd(izero) = 0.d0
      end do
      call propagate_ks(
     &   x0, xd0, cal0,
     &   nrev_zone, istep_prop, tole_prop,
     &   n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &   bn_cand, idrag_loc, WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD,
     &   PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   max_traj, idump_loc,
     &   traj_jd, traj_x, traj_xd,
     &   exit_code)

      if (exit_code .eq. 2) then
         ierr = 2
         return
      end if

      ntraj = 0
      do ipt = 1, max_traj
         if (traj_jd(ipt) .le. 0.d0) goto 16
         if (traj_jd(ipt) - epoch0_jd .gt.
     &       zone_dur + 0.5d0) goto 16
         ntraj = ntraj + 1

         r_mag = dsqrt(traj_x(1,ipt)**2 +
     &                 traj_x(2,ipt)**2 +
     &                 traj_x(3,ipt)**2)
         v_mag = dsqrt(traj_xd(1,ipt)**2 +
     &                 traj_xd(2,ipt)**2 +
     &                 traj_xd(3,ipt)**2)

         sma_pt = 1.d0 / (2.d0/r_mag - v_mag**2/amue)

         hx = traj_x(2,ipt)*traj_xd(3,ipt)
     &      - traj_x(3,ipt)*traj_xd(2,ipt)
         hy = traj_x(3,ipt)*traj_xd(1,ipt)
     &      - traj_x(1,ipt)*traj_xd(3,ipt)
         hz = traj_x(1,ipt)*traj_xd(2,ipt)
     &      - traj_x(2,ipt)*traj_xd(1,ipt)
         h_mag = dsqrt(hx**2 + hy**2 + hz**2)
         p_sl = h_mag**2 / amue
         ecc_pt = dsqrt(dabs(1.d0 - p_sl/sma_pt))

         traj_ha(ntraj) = sma_pt*(1.d0+ecc_pt) - R_Earth
         traj_t(ntraj) = traj_jd(ipt)
      end do
 16   continue

      if (ntraj .lt. 1) then
         ierr = 2
         return
      end if

      ssq = 0.d0
      do i = 1, nobs
         call rsm_tinterp(traj_t, traj_ha, ntraj, epoch0_jd,
     &      tobs(i), pred_ap)
         ssq = ssq + (pred_ap - apobs(i))**2
      end do
      rms_out = dsqrt(ssq / dble(nobs))

      return
      end

c     ============================================================
      subroutine rsm_refine(
     &   a_zone1, e_zone1, inc_zone1, raan_zone1, aop_zone1, ma_zone1,
     &   epoch0_jd, zone_dur,
     &   e_start, bn_start,
     &   e_lo, e_hi, bn_lo, bn_hi,
     &   tobs, apobs, nobs,
     &   istep_prop, tole_prop,
     &   n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &   WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD,
     &   PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   max_eval,
     &   e_refined, bn_refined, rms_refined, neval)
c
c     Issue #40: bounded local refinement of the GA's (interpolated-
c     fitness) optimum using direct-propagation fitness
c     (rsm_direct_eval) instead of the RSM bilinear interpolation.
c     Coordinate pattern search (Hooke-Jeeves style): starts at
c     (e_start, bn_start), tries the 4 axis-aligned neighbors at the
c     current step size, moves to the best improving one and repeats
c     at the same step, or halves the step if none improve. Stops at
c     max_eval direct-propagation calls (~20-50 recommended -- bounds
c     the added cost to a fixed, small multiple of one zone's own RSM
c     grid cost, not the GA's full search). Uses only in-zone TLE data
c     (tobs/apobs) -- no future/oracle information, unlike issue #41's
c     rejected mechanism, so (unlike #41) this is usable in a real,
c     as-yet-undecayed prediction, not just an offline ceiling study.
c
c     Candidates are clamped to [e_lo,e_hi] x [bn_lo,bn_hi] -- the same
c     bounds the GA itself searched -- so the refined point never
c     leaves the physically-vetted search range.
c
c     Output
c       e_refined, bn_refined, rms_refined : best point found (direct
c                                             RMS), including the
c                                             starting point if no
c                                             neighbor ever improved on
c                                             it
c       neval : number of direct-propagation calls actually used
c                (0 if the starting point itself failed to propagate --
c                caller should keep the GA's original result in that
c                case)
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision a_zone1, e_zone1, inc_zone1
      double precision raan_zone1, aop_zone1, ma_zone1
      double precision epoch0_jd, zone_dur
      double precision e_start, bn_start
      double precision e_lo, e_hi, bn_lo, bn_hi
      integer nobs
      double precision tobs(nobs), apobs(nobs)
      integer istep_prop
      double precision tole_prop
      integer n_force(3), ngeo_deg, nsun_deg, nmoon_deg
      double precision WE_rot, EPS_f, FR_rot
      double precision CR_srp, AM_srp
      integer IPSR, ISHAD
      double precision PSR_srp, amuS, amuM
      integer ndim_atm
      double precision ALT_atm(ndim_atm)
      double precision DEN_atm(ndim_atm)
      double precision SCH_atm(ndim_atm)
      integer max_eval
      double precision e_refined, bn_refined, rms_refined
      integer neval

      double precision e_best, bn_best, rms_best
      double precision de, dbn, rms_try
      double precision e_cand(4), bn_cand(4)
      integer ierr_de, j, jbest
      logical improved

      neval = 0
      e_refined = e_start
      bn_refined = bn_start
      rms_refined = 0.d0

      call rsm_direct_eval(
     &   a_zone1, e_zone1, inc_zone1, raan_zone1, aop_zone1, ma_zone1,
     &   epoch0_jd, zone_dur, e_start, bn_start,
     &   tobs, apobs, nobs, istep_prop, tole_prop,
     &   n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &   WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD,
     &   PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   rms_best, ierr_de)
      neval = neval + 1
      if (ierr_de .ne. 0) then
         neval = 0
         return
      end if

      e_best = e_start
      bn_best = bn_start
      de = 0.05d0 * (e_hi - e_lo)
      dbn = 0.05d0 * (bn_hi - bn_lo)

 10   continue
      if (neval .ge. max_eval) goto 90
      if (de .lt. 1.d-6*(e_hi-e_lo) .and.
     &    dbn .lt. 1.d-6*(bn_hi-bn_lo)) goto 90

      e_cand(1) = e_best + de
      bn_cand(1) = bn_best
      e_cand(2) = e_best - de
      bn_cand(2) = bn_best
      e_cand(3) = e_best
      bn_cand(3) = bn_best + dbn
      e_cand(4) = e_best
      bn_cand(4) = bn_best - dbn

      improved = .false.
      jbest = 0
      do j = 1, 4
         if (e_cand(j) .lt. e_lo) e_cand(j) = e_lo
         if (e_cand(j) .gt. e_hi) e_cand(j) = e_hi
         if (bn_cand(j) .lt. bn_lo) bn_cand(j) = bn_lo
         if (bn_cand(j) .gt. bn_hi) bn_cand(j) = bn_hi
         if (e_cand(j) .eq. e_best .and.
     &       bn_cand(j) .eq. bn_best) goto 20

         if (neval .ge. max_eval) goto 20

         call rsm_direct_eval(
     &      a_zone1, e_zone1, inc_zone1, raan_zone1, aop_zone1,
     &      ma_zone1, epoch0_jd, zone_dur,
     &      e_cand(j), bn_cand(j),
     &      tobs, apobs, nobs, istep_prop, tole_prop,
     &      n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &      WE_rot, EPS_f, FR_rot,
     &      CR_srp, AM_srp, IPSR, ISHAD,
     &      PSR_srp, amuS, amuM,
     &      ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &      rms_try, ierr_de)
         neval = neval + 1

         if (ierr_de .eq. 0 .and. rms_try .lt. rms_best) then
            if (jbest .eq. 0) then
               rms_best = rms_try
               jbest = j
               improved = .true.
            else if (rms_try .lt. rms_best) then
               rms_best = rms_try
               jbest = j
            end if
         end if
 20      continue
      end do

      if (improved) then
         e_best = e_cand(jbest)
         bn_best = bn_cand(jbest)
      else
         de = de * 0.5d0
         dbn = dbn * 0.5d0
      end if

      goto 10

 90   continue
      e_refined = e_best
      bn_refined = bn_best
      rms_refined = rms_best

      return
      end

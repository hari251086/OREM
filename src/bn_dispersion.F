c     ============================================================
c     bn_dispersion.F — BN-sweep re-entry dispersion band (issue #43)
c
c     For an object with no observed re-entry epoch, compute_rpe's
c     Mode 1 (RPE vs. observed) is unavailable and Mode 2 (RPE vs.
c     ensemble median) is self-referential -- every zone can share
c     the same systematic bias and still "agree." Instead of a point
c     estimate, sweep the ballistic number BN across a zone's own
c     already-fitted search window (the range already treated as
c     "consistent with the observed TLE" by that zone's RSM/GA fit)
c     and re-propagate, producing a dispersion band on the predicted
c     re-entry date.
c
c     Directly re-implements the precedent documented in the ISRO
c     monograph deep-read (Fig 8.6, GSLV-D1 spent-stage lifetime):
c     real historical/forecast space weather is used UNCHANGED (via
c     propagate_ks's own epoch-resolved swx.F consumption); only BN
c     is varied between sweep points. No Monte Carlo, no synthesized
c     forecast-error model, no RNG.
c
c     Author: Harishkumar Sellamuthu
c     ============================================================
      subroutine bn_dispersion_band(
     &   x0, xd0, cal0, bn_lo, bn_hi, n_sweep,
     &   n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &   idrag_flag, WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD,
     &   PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   nrev_re, max_reentry,
     &   bn_band, t_band, n_band_valid)
c
c     Input
c       x0(3), xd0(3), cal0(6) : the zone's already-fitted initial
c                                 state (same state the deterministic
c                                 propagate_ks call already used)
c       bn_lo, bn_hi            : BN sweep range (kg/m2) -- pass the
c                                 zone's own bn_lo_zone/bn_hi_zone
c       n_sweep                 : number of sweep points, 2..NSWMAX
c       n_force, ngeo_deg, nsun_deg, nmoon_deg, idrag_flag, WE_rot,
c       EPS_f, FR_rot, CR_srp, AM_srp, IPSR, ISHAD, PSR_srp, amuS,
c       amuM, ALT_atm, DEN_atm, SCH_atm, ndim_atm : identical
c                                 force-model args passed straight
c                                 through to propagate_ks, unchanged
c                                 between sweep points
c       nrev_re, max_reentry    : same revolution cap / trajectory
c                                 buffer size the zone's own
c                                 deterministic call used
c
c     Output
c       bn_band(NSWMAX), t_band(NSWMAX) : BN value and resulting
c                                 re-entry JD for each sweep point
c                                 that actually re-entered within
c                                 the cap (unfilled slots are 0)
c       n_band_valid             : how many of n_sweep points
c                                 produced a re-entry within cap
c     ============================================================
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      integer NSWMAX
      parameter(NSWMAX = 10)

c     --- Input arguments ---
      double precision x0(3), xd0(3), cal0(6)
      double precision bn_lo, bn_hi
      integer n_sweep
      integer n_force(3)
      integer ngeo_deg, nsun_deg, nmoon_deg
      integer idrag_flag
      double precision WE_rot, EPS_f, FR_rot
      double precision CR_srp, AM_srp
      integer IPSR, ISHAD
      double precision PSR_srp, amuS, amuM
      integer ndim_atm
      double precision ALT_atm(ndim_atm)
      double precision DEN_atm(ndim_atm)
      double precision SCH_atm(ndim_atm)
      integer nrev_re, max_reentry

c     --- Output arguments ---
      double precision bn_band(NSWMAX), t_band(NSWMAX)
      integer n_band_valid

c     --- Local variables ---
c     max_reentry is a dummy arg here (unlike orem_run's own
c     compile-time parameter of the same name), so these are
c     automatic arrays -- SAVE is illegal on them (ifx #6754); rely
c     on -heap-arrays for placement instead, same requirement as the
c     rest of this codebase's large local arrays.
      double precision traj_jd(max_reentry)
      double precision traj_x(3, max_reentry)
      double precision traj_xd(3, max_reentry)
      double precision bn_i
      integer i, idump_i, exit_code_i, ns

      n_band_valid = 0
      do i = 1, NSWMAX
         bn_band(i) = 0.d0
         t_band(i) = 0.d0
      end do

      if (n_sweep .lt. 2) return
      ns = n_sweep
      if (ns .gt. NSWMAX) ns = NSWMAX
      if (bn_hi .le. bn_lo) return

      do i = 1, ns
         bn_i = bn_lo + (bn_hi - bn_lo) * dble(i-1) / dble(ns-1)

         idump_i = 1
         call propagate_ks(
     &      x0, xd0, cal0,
     &      nrev_re, 360, 1.d-10,
     &      n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &      bn_i, idrag_flag, WE_rot, EPS_f, FR_rot,
     &      CR_srp, AM_srp, IPSR, ISHAD,
     &      PSR_srp, amuS, amuM,
     &      ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &      max_reentry, idump_i,
     &      traj_jd, traj_x, traj_xd,
     &      exit_code_i)

         if (exit_code_i .eq. 1) then
            n_band_valid = n_band_valid + 1
            bn_band(n_band_valid) = bn_i
            t_band(n_band_valid) = traj_jd(idump_i)
         end if
      end do

      return
      end

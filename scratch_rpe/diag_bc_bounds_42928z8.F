c     diag_bc_bounds_42928z8.F -- issue #32 follow-up: does Sharma et al.
c     (2009) KSGEN's +/-10% ballistic-coefficient dispersion-bounds
c     technique give a useful uncertainty band, compared to OREM's own
c     existing ensemble spread (cross-zone t_std/ens_rpe_pct)?
c
c     Sharma's paper perturbs the FITTED BC by +/-10% and re-propagates,
c     reporting the resulting re-entry-time spread as the prediction's
c     uncertainty bounds -- validated (their own words) on SROSS-C2 and
c     SL-12 R/B, where the true re-entry time fell inside the band both
c     times. This is a DIFFERENT uncertainty mechanism than OREM's: OREM
c     never perturbs a single fit -- it looks at how much bn_opt/reentry_jd
c     already disagree ACROSS independently-fitted zones (the ensemble).
c
c     Reuses 42928 Z8's exact SGP4-osculating IC from the original GMAT
c     cross-check (scratch_gmat/gmat_xval_42928z8_2026.script header:
c     sma=7524.293124, e_osc=0.133643, inc=19.094467, raan=210.267993,
c     aop=73.147372, ma=301.305632, epoch JD=2458455.471721) with the
c     perigee-preserving SMA reshape already used by
c     scratch_gmat/bulge_xval_42928z8.F, applied to the POST-BULGE
c     e_opt/bn_opt this exact IC was actually fitted to at the time
c     (e_opt=0.051814, bn_opt=128.532, per that script's own header) --
c     NOT the current rpe_campaign_weather.csv zone-8 values, which now
c     belong to a different TLE window after the zone_select recency fix
c     and max_zone_days 10->20 change shifted zone boundaries. A first
c     attempt paired this fixed IC with the current-generation
c     e_opt/bn_opt and got a reentry_jd 37 days off the real pipeline's
c     own current answer -- confirming IC and fit generation must match.
c     What's tested here is the +/-10% MECHANISM's behavior on a genuine
c     (IC, e_opt, bn_opt) triple, not reproduction of today's exact
c     campaign figure for this object.
c
c     Runs propagate_ks 3x from the SAME IC with BN = 0.9x/1.0x/1.1x the
c     fitted value, using OREM's OWN re-entry criterion (h_alt<80km,
c     propagate_ks.F's own exit_code=1 -- not GMAT's 102km convention),
c     matching exactly what orem.F itself calls to produce reentry_jd.
c     Reports the resulting 3-point band vs the real observed re-entry
c     date (2019-03-03) and vs 42928's CURRENT already-computed ensemble
c     numbers (t_mean=2458593.3515, t_std=74.07d, ens_rpe_pct=9.100%,
c     from the present-day rpe_campaign_weather.csv) -- the ensemble
c     numbers are current-generation while this script's own 3-point
c     band is necessarily the older generation the fixed IC belongs to;
c     the comparison is about BAND WIDTH (a roughly generation-invariant
c     property of sensitivity to a 10% BN change), not about matching
c     absolute reentry dates across generations.
c
c     RESULT (see issue #32 comment): band width 5.4 days vs the current
c     ensemble's t_std=74.07 days -- the +/-10% BC band is ~14x
c     NARROWER than OREM's actual empirical cross-zone disagreement for
c     this object. Not implemented as a reported uncertainty metric.
c
c     Compile (fpm-layout paths, KSROP dependency under
c     build/dependencies/ksrop/src):
c       ifx /heap-arrays /F:16777216 scratch_rpe\diag_bc_bounds_42928z8.F
c           src\propagate_ks.F src\swx.F
c           build\dependencies\ksrop\src\Subrouts.F
c           build\dependencies\ksrop\src\Legendre.F
c           /exe:scratch_rpe\diag_bc_bounds_42928z8.exe
c
c     Analysis only -- not part of the build/test suite.
      program diag_bc_bounds_42928z8
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291, maxp = 20000)
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm, ierr_sw, ierr_a2
      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      double precision tj(maxp), tx(3,maxp), td(3,maxp)
      integer n_force(3), ec, idump
      double precision rp, sma_osc, e_osc, jd0
      double precision e_opt, bn_opt, bn_try(3), reentry_jd(3)
      double precision t_obs_cal(6), t_obs_jd, rpe_pct(3)
      integer i, k

      call init_constants()

      call sw_load('input/SW-All.csv', ierr_sw)
      call atm2d_load('input/ATM2D.DAT', ierr_a2)
      write(*,'(A,I2,A,I2)') ' sw_load ierr=', ierr_sw,
     &   '  atm2d_load ierr=', ierr_a2

      natm = ndim_a
      ALT_a(1) = 60.d0
      do i = 2, 141
         ALT_a(i) = ALT_a(i-1) + 1.d0
      end do
      do i = 142, natm
         ALT_a(i) = ALT_a(i-1) + 2.d0
      end do
      open(unit=60, file='input/ATM.DAT', status='old')
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

      sma_osc = 7524.293124d0
      e_osc   = 0.133643d0
      rp = sma_osc * (1.d0 - e_osc)
      jd0 = 2458455.471721d0
      call jd2cal(jd0, cal0)

      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

c     NOTE: zone boundaries have shifted since this IC/epoch was derived
c     (zone_select recency fix + max_zone_days 10->20, both after the
c     original bulge cross-check) -- "zone 8" under the CURRENT pipeline
c     is a different TLE window than it was when this IC was captured, so
c     pairing this fixed IC with the CURRENT rpe_campaign_weather.csv
c     zone-8 e_opt/bn_opt (0.081305/116.967) is inconsistent: first
c     attempt did that and got reentry_jd 37 days off the real pipeline's
c     own zone-8 answer. Using the POST-BULGE, SAME-GENERATION values
c     this exact IC was actually fitted to instead (from
c     scratch_gmat/bulge_xval_42928z8.F's own header) keeps IC and fit
c     internally consistent, even though it's no longer "the latest"
c     production number -- what's being tested here is the +/-10%
c     MECHANISM's behavior on a genuine GA-fitted (IC, e_opt, bn_opt)
c     triple, not reproduction of today's exact campaign figure.
      e_opt  = 0.051814d0
      bn_opt = 128.532d0

      bn_try(1) = bn_opt * 0.9d0
      bn_try(2) = bn_opt
      bn_try(3) = bn_opt * 1.1d0

c     Real observed re-entry date for 42928 (2019-03-03, per
c     scratch_rpe/rpe_campaign_7obj_issue31.F's own t_obs array)
      t_obs_cal(1) = 2019.d0
      t_obs_cal(2) = 3.d0
      t_obs_cal(3) = 3.d0
      t_obs_cal(4) = 0.d0
      t_obs_cal(5) = 0.d0
      t_obs_cal(6) = 0.d0
      call cal2jd(t_obs_cal, t_obs_jd)

      do k = 1, 3
         a_kep(1) = rp / (1.d0 - e_opt)
         a_kep(2) = e_opt
         a_kep(3) = 19.094467d0
         a_kep(4) = 210.267993d0
         a_kep(5) = 73.147372d0
         a_kep(6) = 301.305632d0
         a_kep(7) = 0.d0

         call oe2car(a_kep, x0, xd0, 1.d-12)

         idump = 1
         call propagate_ks(
     &      x0, xd0, cal0,
     &      3000, 360, 1.d-10,
     &      n_force, 20, 2, 3,
     &      bn_try(k), 1, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &      1.2d0, 0.01d0, 1, 2,
     &      4.56d-6, 1.32712440018d11, 4.902801076d3,
     &      ALT_a, DEN_a, SCH_a, natm,
     &      maxp, idump,
     &      tj, tx, td,
     &      ec)

         if (ec .eq. 1) then
            reentry_jd(k) = tj(idump)
            rpe_pct(k) = (reentry_jd(k) - t_obs_jd)
     &         / (t_obs_jd - jd0) * 100.d0
            write(*,'(A,F8.3,A,F14.4,A,F8.2,A)')
     &         ' BN=', bn_try(k), '  reentry_jd=', reentry_jd(k),
     &         '  RPE=', rpe_pct(k), '%'
         else
            reentry_jd(k) = 0.d0
            write(*,'(A,F8.3,A,I2)')
     &         ' BN=', bn_try(k), '  did not reach 80km, exit_code=',
     &         ec
         end if
      end do

      write(*,'(/,A)') ' --- +/-10% BC dispersion-bounds summary ---'
      if (reentry_jd(1) .gt. 0.d0 .and. reentry_jd(3) .gt. 0.d0) then
         write(*,'(A,F10.4,A,F10.4,A)') ' band: [',
     &      min(reentry_jd(1),reentry_jd(3)), ', ',
     &      max(reentry_jd(1),reentry_jd(3)), '] JD'
         write(*,'(A,F8.3,A)') ' band width = ',
     &      dabs(reentry_jd(3)-reentry_jd(1)), ' days'
         write(*,'(A,L1)') ' true re-entry date inside band? ',
     &      (t_obs_jd .ge. min(reentry_jd(1),reentry_jd(3)) .and.
     &       t_obs_jd .le. max(reentry_jd(1),reentry_jd(3)))
      end if
      write(*,'(A)') ' (compare to existing ensemble for 42928:'//
     &   ' t_std=74.07d, ens_rpe_pct=9.100%, from'//
     &   ' rpe_campaign_weather.csv)'

      stop
      end

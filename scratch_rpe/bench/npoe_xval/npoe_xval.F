c     ============================================================
c     npoe_xval.F -- ground-truth cross-validation against NPOE
c     (the independent third-party reference propagator the original
c     Sharma-lineage papers validated against -- see the literature
c     review in project memory / issue #36).
c
c     Reference data: E:\Research\1. R&D\Re-entry\COSPAR ASR\42928\
c     Zone - 2\ -- NPOE's own input file (RSM Surfaces\Z2.in) gives
c     exact IC + force model for a real propagation of object 42928
c     from 2017-10-03 14:35:56 UTC, and NPOE's own output (Z2.out
c     mean elements, Z2osc.out osculating elements) at 10-minute
c     cadence for ~444 days through to near re-entry.
c
c     This program propagates propagate_ks from the SAME IC under the
c     closest achievable match to NPOE's force model, dumps osculating
c     apogee altitude vs elapsed days to a CSV, for comparison against
c     Z2osc.out (parsed separately in Python -- see npoe_compare.py).
c
c     IC (from Z2.in, true anomaly converted to mean anomaly to match
c     a_kep's convention -- see file header comment at the call site):
c       epoch   2017-10-03 14:35:56 UTC
c       a       9594.412450 km
c       e       0.3185937
c       i       19.1779 deg
c       AOP     332.445197 deg
c       RAAN    36.853002 deg
c       MA      96.766670 deg (from true anomaly 129.53521 deg --
c               matches the TLE's own MA=96.773503 at this epoch to
c               0.007 deg, cross-validating the conversion itself)
c
c     Force model (from Z2.in, mapped to OREM's flags):
c       drag: BN = mass/(Cd*A) = 112.1098/(1*1) = 112.1098 kg/m2
c       gravity: NPOE degree/order 6/6 (zonal+tesseral) -- OREM is
c         zonal-only, so ngeo_deg=6 is the closest achievable match,
c         NOT identical (OREM structurally cannot include tesserals)
c       SRP: off (IPSR=0), matching NPOE's "no"
c       solar+lunar third-body: on, nsun_deg=2/nmoon_deg=3 (OREM's
c         standard choice -- NPOE doesn't expose a degree parameter)
c       atmosphere: OREM's static Jacchia-71 table (ATM.DAT) --
c         NPOE used density model "5" (a Jacchia 1970 variant),
c         same lineage but not bit-identical; NPOE's solar/geomagnetic
c         inputs (F10=72, F10B=72, Kp~10) are constant, matching
c         OREM's static (non-epoch-resolved) mode
c
c     Not wired into fpm.toml. To reproduce, temporarily add:
c       [[executable]]
c       name = "npoe_xval"
c       source-dir = "scratch_rpe/bench/npoe_xval"
c       main = "npoe_xval.F"
c     then: fpm run npoe_xval --compiler ifx --flag "-heap-arrays"
c     (Python comparison: scratch_rpe/bench/npoe_xval/npoe_compare.py)
c
c     RESULTS (2026-08-01, ifx, this machine) -- see issue #38 and
c     README.md Version History for the full writeup:
c       - bn_scale_test.F (companion, same dir) confirms decay rate
c         scales as ~1/BN as the drag formula requires -- mechanism
c         is not broken.
c       - Early-window (10-day) SMA decay rate: NPOE 3.013 km/day vs
c         OREM 1.946 km/day -- NPOE ~55% faster, independently
c         confirming the Phase 5 GMAT cross-check's "36-68% thin"
c         density finding via a different tool and object/zone.
c       - This modest per-rev bias COMPOUNDS over ~4260 revolutions
c         (460 days): by day 430, apogee altitude differs by ~3500km
c         (331% relative). New insight, not previously characterized:
c         this specifically threatens the RE-ENTRY PROPAGATION step
c         (single constant BN run to decay), not the per-zone RSM/GA
c         fit (recalibrated against fresh TLE every ~20 days).
c       - Caveat: this test starts from Zone 2, ~516 days before real
c         re-entry -- much longer than production's typical
c         latest-zone-to-re-entry horizon. Magnitude here is closer
c         to a stress test than a typical-case error bar.
c     ============================================================
      program npoe_xval
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291)
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      integer n_force(3)

      parameter(max_pts = 5000)
      double precision traj_jd(max_pts)
      double precision traj_x(3, max_pts)
      double precision traj_xd(3, max_pts)
      integer exit_code, idump

      double precision epoch0_jd
      double precision r_mag, v_mag, sma_pt, ecc_pt
      double precision hx, hy, hz, h_mag, p_sl, ha_pt
      integer ipt

      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

c     --- Initial conditions (Z2.in, MA converted from true anomaly) ---
      a_kep(1) = 9594.412450d0
      a_kep(2) = 0.3185937d0
      a_kep(3) = 19.1779d0
      a_kep(4) = 36.853002d0
      a_kep(5) = 332.445197d0
      a_kep(6) = 96.766670d0
      a_kep(7) = 0.d0

      call oe2car(a_kep, x0, xd0, 1.d-12)

      cal0(1) = 2017.d0
      cal0(2) = 10.d0
      cal0(3) = 3.d0
      cal0(4) = 14.d0
      cal0(5) = 35.d0
      cal0(6) = 56.d0
      call cal2jd(cal0, epoch0_jd)

      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

      open(unit=95, file=
     &   'scratch_rpe/bench/npoe_xval/orem_trajectory.csv',
     &   status='replace')
      write(95,'(A)') 'day,apogee_alt_km,sma_km,ecc'

c     Propagate ~460 days (~4260 revs at this ~2.6hr period) to cover
c     NPOE's full ~444-day run through to near re-entry, one dump per
c     revolution (propagate_ks's own convention).
      idump = 1
      call propagate_ks(
     &   x0, xd0, cal0,
     &   4260, 360, 1.d-10,
     &   n_force, 6, 2, 3,
     &   112.1098d0, 1, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 0, 0,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   max_pts, idump,
     &   traj_jd, traj_x, traj_xd,
     &   exit_code)

      write(*,*) 'exit_code=', exit_code, ' idump=', idump

      do ipt = 1, idump
         r_mag = dsqrt(traj_x(1,ipt)**2 + traj_x(2,ipt)**2
     &                + traj_x(3,ipt)**2)
         v_mag = dsqrt(traj_xd(1,ipt)**2 + traj_xd(2,ipt)**2
     &                + traj_xd(3,ipt)**2)
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
         ha_pt = sma_pt*(1.d0+ecc_pt) - R_Earth

         write(95,'(F10.4,A,F12.4,A,F12.4,A,F10.6)')
     &      traj_jd(ipt)-epoch0_jd, ',', ha_pt, ',', sma_pt, ',',
     &      ecc_pt
      end do

      close(95)

      write(*,*) 'wrote orem_trajectory.csv, ', idump, ' points'

      stop
      end

c     ============================================================
      subroutine read_atm_c(alt, den, sch, maxn, natm)
      implicit double precision (a-h, o-z)
      integer maxn, natm
      double precision alt(maxn), den(maxn), sch(maxn)

      natm = 291
      if (natm .gt. maxn) natm = maxn
      alt(1) = 60.d0
      do i = 2, 141
         alt(i) = alt(i-1) + 1.d0
      end do
      do i = 142, natm
         alt(i) = alt(i-1) + 2.d0
      end do

      open(unit=88, file='input/ATM.DAT', status='old',
     &     err=900)
      read(88,110) (sch(i), i=1,61)
      read(88,110) (sch(i), i=62,natm)
      read(88,111) (den(i), i=1,61)
      read(88,111) (den(i), i=62,natm)
      close(88)
 110  format(8f10.3)
 111  format(8e10.3)

      do i = 1, natm
         alt(i) = alt(i) * 1000.d0
         sch(i) = sch(i) * 1000.d0
         den(i) = den(i) * 1.0d10
      end do
      return
 900  natm = 0
      return
      end

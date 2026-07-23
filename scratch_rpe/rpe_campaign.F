c     rpe_campaign.F — 30-object drag-enabled RPE campaign (issue #16)
c     -------------------------------------------------------
c     Runs orem_run with IDRAG=1 and the full force model (the
c     test_e2e v1.11 config: geo=20, sun=2, moon=3, SRP on) on all
c     ten validation objects, and reports per-zone RPE plus the
c     per-object ensemble prediction (t_mean/t_std -> ensemble RPE).
c     Purpose: quantify post-v1.18 accuracy across the whole set --
c     test_e2e only exercises 42928 (+z0, 39615 z1, 35497 z2) and
c     test_reentry runs IDRAG=0. Results feed #16 (<10% target) and
c     the #13 ensemble-reporting design.
c
c     Extended from 7 to 10 objects (2026-07-19) with the three real-
c     decay cases from the same day's earlier individual validation
c     runs (11550/59347/40943, input/orem_<norad>.cfg) -- deliberately
c     run through the SAME nzones_max=8 (mxz below) as the original 7,
c     not the nzones_max=50 used for OREM-Watchlist's operational
c     37398 case (issue #9 there): the point of this campaign is to
c     test whether OREM's shipped production default generalizes
c     beyond the curated original set, not to tune per object. Result
c     (issue #29): 2/3 of the new objects fit dramatically worse at
c     nzones_max=8 than they had individually at nzones_max=50.
c
c     Extended again 10 -> 20 objects (2026-07-19), same session: 10
c     more real, confirmed-decay GTO/HEO rocket bodies, found by
c     filtering the already-on-disk satcat.csv (catalog/latest/
c     satcat.csv from that day's live catalog pull -- no new
c     Space-Track queries for the search itself) for
c     APOGEE>8000 km, PERIGEE<3000 km, OBJECT_TYPE=ROCKET BODY,
c     DECAY not null, decay after 2015 -- then fetching each
c     candidate's full TLE history via OREM-Watchlist's
c     fetch_object.py (one first-ever backfill call per object, per
c     SPACE_TRACK_USAGE_POLICY.md). 15 candidates found; excluded 2
c     anomalously high-apogee CZ-5 R/Bs (287434/537024 km -- almost
c     certainly trans-lunar-injection stages, a fundamentally
c     different regime than GTO) and 2 with too few TLEs (64785: 5,
c     60607: 21 -- not enough for min_pts=8 in any 10-day window),
c     leaving exactly 10 with 72-3203 TLEs each.
c
c     Extended again 20 -> 30 objects (2026-07-22): the ROCKET BODY-
c     only filter against the same on-disk satcat.csv snapshot (still
c     dated 2026-07-18 -- no fresher pull) returned ZERO new
c     candidates -- the 13-14 objects it matches were already all used
c     or excluded by the 20-object round. Broadened OBJECT_TYPE to
c     also include DEBRIS (same precedent as object 20, 48259, a
c     Cosmos 1124 DEB) -- 45 candidates. Fetched 16 via
c     fetch_object.py; excluded 2 for too few TLEs (40598: 2, 40777:
c     14), 1 for a saturated terminal ndot (36518: .99999999, the
c     exact issue #30 hang pattern -- not risked), and 2 for a
c     multi-decade tracking gap making the catalogued DECAY date
c     unreachable from any tracked TLE (38587: last TLE 2013 vs decay
c     2025; 27900: last TLE 2004 vs decay 2020) -- leaving exactly 10.
c     Several of the 10 kept still have real, non-trivial gaps between
c     their last tracked TLE and catalogued decay (28572: ~6 yr;
c     23647: ~2.8 yr; 44591: ~1.8 yr; 46429/39802: ~1 yr) -- kept
c     deliberately rather than cherry-picked away, since sparse
c     terminal tracking is exactly the real-world failure mode issue
c     #29 is about, not a reason to curate a rosier set.
c
c     Compile:
c       ifx /heap-arrays /F:16777216 scratch_rpe\rpe_campaign.F
c           orem.F rsm.F ga.F tle_evolution.F tle_filter.F
c           zone_select.F ksrop\propagate_ks.F ksrop\Subrouts.F
c           ksrop\TLEread.F ksrop\Legendre.F
c           /exe:scratch_rpe\rpe_campaign.exe
c     Run from the OREM repo root:
c       scratch_rpe\rpe_campaign.exe
c     Output: scratch_rpe\rpe_campaign.csv
c     -------------------------------------------------------
      program rpe_campaign

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

c     mxz=8 (v1.21 experiment; was 4): zone_select returns the top-R2
c     candidates, so a higher cap admits LATER zones -- and v1.20
c     showed the latest zone is the primary estimator. This run
c     measures how much later the latest zone gets and what it does
c     to the latest-zone RPE.
      parameter(ndim_a = 291, mxz = 8, nobj = 30)

      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr, nzv
      integer zstat(mxz)

      character*40 tle_file(nobj)
      character*24 obj_name(nobj)
      integer norad(nobj)
      double precision oyr(nobj), omo(nobj), ody(nobj)
      double precision t_obs_jd, horiz, rpe_ens
      integer io, iz, nre

c     11550's date corrected 2025-04-03 (inferred from its terminal
c     TLE signature) -> 2025-04-10, the real confirmed SATCAT DECAY
c     date found while searching for the next 10 objects. 59347/40943
c     dates are real, confirmed (matching their individual
c     input/orem_<norad>.cfg files). Objects 11-20 are real, confirmed
c     decay dates straight from satcat.csv's DECAY column.
c
c     Object 20 was originally 42985 (Falcon 9 R/B, decay 2018-10-27):
c     dropped after its zone 8 hung indefinitely in ga_optimize (not
c     Kesolve itself -- that has a hard 100-iteration cap, see
c     ksrop/Subrouts.F -- the hang was downstream, in the GA searching
c     a precomputed surface, no propagate_ks/Kesolve calls involved).
c     42985's terminal TLEs have a saturated ndot (.99999999, the
c     field's max representable value -- a degenerate catalog-fit
c     artifact, not real physics) likely poisoning the zone-8 RSM
c     surface with NaN/Inf that ga_optimize isn't robust to. Same
c     pattern independently seen on 33587's terminal TLE (issue #27).
c     Root cause not fixed here -- replaced with 48259 (Cosmos 1124
c     debris, decay 2026-06-18, 1473 TLEs, clean non-saturated
c     terminal ndot=0.157) instead of chasing the ga_optimize bug in
c     the middle of this campaign; that hang is real and worth its
c     own issue.
      data norad / 42928, 35497, 37151, 39615, 27526, 32007, 37819,
     &             11550, 59347, 40943,
     &             66587, 60328, 61734, 57804, 56758,
     &             30799, 44187, 35009, 41553, 48259,
     &             27906, 27882, 39802, 28572, 46429,
     &             44591, 41695, 52205, 45349, 23647 /
      data oyr / 2019.d0, 2016.d0, 2015.d0, 2017.d0, 2012.d0,
     &           2010.d0, 2013.d0, 2025.d0, 2026.d0, 2020.d0,
     &           2026.d0, 2025.d0, 2025.d0, 2025.d0, 2024.d0,
     &           2024.d0, 2021.d0, 2020.d0, 2019.d0, 2026.d0,
     &           2018.d0, 2019.d0, 2019.d0, 2020.d0, 2021.d0,
     &           2022.d0, 2022.d0, 2022.d0, 2022.d0, 2024.d0 /
      data omo / 3.d0, 10.d0, 12.d0, 9.d0, 5.d0, 6.d0, 9.d0,
     &           4.d0, 6.d0, 2.d0,
     &           2.d0, 4.d0, 3.d0, 1.d0, 8.d0,
     &           5.d0, 6.d0, 8.d0, 1.d0, 6.d0,
     &           11.d0, 4.d0, 11.d0, 10.d0, 10.d0,
     &           3.d0, 5.d0, 6.d0, 10.d0, 1.d0 /
      data ody / 3.d0, 31.d0, 3.d0, 15.d0, 9.d0, 6.d0, 12.d0,
     &           10.d0, 8.d0, 20.d0,
     &           10.d0, 14.d0, 22.d0, 3.d0, 25.d0,
     &           12.d0, 15.d0, 31.d0, 6.d0, 18.d0,
     &           19.d0, 7.d0, 5.d0, 4.d0, 18.d0,
     &           3.d0, 15.d0, 26.d0, 9.d0, 28.d0 /

      tle_file(1) = 'input/example_42928.tle.txt'
      tle_file(2) = 'input/example_35497.tle.txt'
      tle_file(3) = 'input/example_37151.tle.txt'
      tle_file(4) = 'input/example_39615.tle.txt'
      tle_file(5) = 'input/example_27526.tle.txt'
      tle_file(6) = 'input/example_32007.tle.txt'
      tle_file(7) = 'input/example_37819.tle.txt'
      tle_file(8) = 'input/example_11550.tle.txt'
      tle_file(9) = 'input/example_59347.tle.txt'
      tle_file(10) = 'input/example_40943.tle.txt'
      tle_file(11) = 'input/example_66587.tle.txt'
      tle_file(12) = 'input/example_60328.tle.txt'
      tle_file(13) = 'input/example_61734.tle.txt'
      tle_file(14) = 'input/example_57804.tle.txt'
      tle_file(15) = 'input/example_56758.tle.txt'
      tle_file(16) = 'input/example_30799.tle.txt'
      tle_file(17) = 'input/example_44187.tle.txt'
      tle_file(18) = 'input/example_35009.tle.txt'
      tle_file(19) = 'input/example_41553.tle.txt'
      tle_file(20) = 'input/example_48259.tle.txt'
      tle_file(21) = 'input/example_27906.tle.txt'
      tle_file(22) = 'input/example_27882.tle.txt'
      tle_file(23) = 'input/example_39802.tle.txt'
      tle_file(24) = 'input/example_28572.tle.txt'
      tle_file(25) = 'input/example_46429.tle.txt'
      tle_file(26) = 'input/example_44591.tle.txt'
      tle_file(27) = 'input/example_41695.tle.txt'
      tle_file(28) = 'input/example_52205.tle.txt'
      tle_file(29) = 'input/example_45349.tle.txt'
      tle_file(30) = 'input/example_23647.tle.txt'
      obj_name(1) = '42928 PSLV-C39'
      obj_name(2) = '35497 Ariane 5 ESC-A'
      obj_name(3) = '37151 Long March 3B'
      obj_name(4) = '39615 Proton-M Briz-M'
      obj_name(5) = '27526 Ariane 5 R/B'
      obj_name(6) = '32007 GSLV R/B'
      obj_name(7) = '37819 Proton-M R/B'
      obj_name(8) = '11550 Molniya debris'
      obj_name(9) = '59347 GTO R/B'
      obj_name(10) = '40943 GTO debris'
      obj_name(11) = '66587 CZ-3B R/B'
      obj_name(12) = '60328 CZ-3B R/B'
      obj_name(13) = '61734 H-3 R/B'
      obj_name(14) = '57804 H-2A R/B'
      obj_name(15) = '56758 Falcon9 R/B'
      obj_name(16) = '30799 Delta2 R/B'
      obj_name(17) = '44187 Falcon9 R/B'
      obj_name(18) = '35009 Fregat R/B'
      obj_name(19) = '41553 Falcon9 R/B'
      obj_name(20) = '48259 Cosmos1124 deb'
      obj_name(21) = '27906 Cosmos1481 deb'
      obj_name(22) = '27882 Cosmos903 deb'
      obj_name(23) = '39802 SL-12 deb'
      obj_name(24) = '28572 SL-12 deb'
      obj_name(25) = '46429 SL-12 deb'
      obj_name(26) = '44591 CZ-3A deb'
      obj_name(27) = '41695 SL-12 deb'
      obj_name(28) = '52205 Cosmos917 deb'
      obj_name(29) = '45349 SL-12 deb'
      obj_name(30) = '23647 Titan34B deb'

      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      open(unit=33, file='scratch_rpe/rpe_campaign.csv',
     &     status='unknown')
      write(33,'(A)') 'norad,zone,e_opt,bn_opt,reentry_jd,rpe_pct,'//
     &   'zstat,t_mean,t_std,ens_rpe_pct'

      do io = 1, nobj
         write(*,'(/,A,A)') ' ===== ', obj_name(io)

         t_obs(1) = oyr(io)
         t_obs(2) = omo(io)
         t_obs(3) = ody(io)
         t_obs(4) = 0.d0
         t_obs(5) = 0.d0
         t_obs(6) = 0.d0

         call orem_run(
     &      tle_file(io), norad(io),
     &      t_obs, mxz, 8, 10.d0, 0.90d0, -1.d0,
     &      80.d0, 160.d0, 1,
     &      20, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &      20, 2, 3,
     &      7.2921150d-5, 3.35281066d-3, 1.d0,
     &      1.2d0, 0.01d0, 1, 2,
     &      4.56d-6, 1.32712440018d11, 4.902801076d3,
     &      ALT_a, DEN_a, SCH_a, natm,
     &      reentry, e_opt, bn_opt, rms_o,
     &      zepoch, nzu,
     &      zstat, nzv,
     &      rpe, t_mean, t_std, ierr)

         if (ierr .ne. 0 .and. ierr .ne. 3) then
            write(*,'(A,I3)') '  pipeline error ierr=', ierr
            write(33,'(I5,A,I3)') norad(io), ',ERR,,,,,,,,', ierr
            goto 90
         end if

c        Ensemble RPE: t_mean vs observed, horizon from zone 1 epoch
         call cal2jd(t_obs, t_obs_jd)
         rpe_ens = 0.d0
         nre = 0
         do iz = 1, nzu
            if (reentry(iz) .gt. 0.d0) nre = nre + 1
         end do
         if (nre .ge. 1 .and. t_obs_jd .gt. 0.d0) then
            horiz = t_obs_jd - zepoch(1)
            if (horiz .gt. 0.d0)
     &         rpe_ens = (t_mean - t_obs_jd)/horiz*100.d0
         end if

         do iz = 1, nzu
            if (reentry(iz) .gt. 0.d0) then
               write(*,'(A,I1,A,F8.4,A,F7.2,A,F8.2,A,I2)')
     &            '  Z', iz, ': e=', e_opt(iz), ' BN=', bn_opt(iz),
     &            ' RPE=', rpe(iz), '%  zstat=', zstat(iz)
            else
               write(*,'(A,I1,A,F8.4,A,F7.2,A,I2)')
     &            '  Z', iz, ': e=', e_opt(iz), ' BN=', bn_opt(iz),
     &            ' (no re-entry)  zstat=', zstat(iz)
            end if
            write(33,'(I5,A,I1,A,F9.6,A,F8.3,A,F13.4,A,F9.3,A,I2,
     &         A,F13.4,A,F8.2,A,F9.3)')
     &         norad(io), ',', iz, ',', e_opt(iz), ',', bn_opt(iz),
     &         ',', reentry(iz), ',', rpe(iz), ',', zstat(iz),
     &         ',', t_mean, ',', t_std, ',', rpe_ens
         end do

         if (nre .ge. 1) then
            write(*,'(A,I1,A,I1,A,F13.3,A,F7.1,A)')
     &         '  ensemble (', nre, ' of ', nzu,
     &         ' zones): mean JD=', t_mean, ' std=', t_std, ' d'
            write(*,'(A,F8.2,A)')
     &         '  ENSEMBLE RPE = ', rpe_ens, '%'
         else
            write(*,'(A)') '  no re-entry in any zone (weak-signal'//
     &         ' fits; lifetime beyond 5-yr cap)'
         end if

 90      continue
c        Progress sentinel: a Monitor watching this log can grep
c        "^ \[PROGRESS\]" for a clean, unambiguous percent-complete
c        marker without polling -- covers both the normal-completion
c        and early-error (goto 90) exits per object.
         write(*,'(A,I3,A,I2,A,I2,A)') ' [PROGRESS] ',
     &      (io*100)/nobj, '% (', io, ' of ', nobj, ' objects done)'
      end do

      close(33)
      write(*,'(/,A)') ' campaign complete -> scratch_rpe/'//
     &   'rpe_campaign.csv'

      stop
      end

c     ============================================================
c     ATM reader (same as test_reentry.F)
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

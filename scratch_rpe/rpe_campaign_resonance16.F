c     rpe_campaign_resonance16.F -- Wang & Gurfil (2017) solar-apsidal
c     "U-turn" resonance restricted campaign.
c     -----------------------------------------------------------
c     Runs OREM's standard orem_run pipeline (identical zone/GA/force-
c     model parameters to the main scratch_rpe/rpe_campaign.F 97-object
c     campaign, for direct comparability) on only the 16 objects (of
c     that same 97) whose real TLE history shows a genuine crossing of
c     the solar-apsidal-resonance condition (Wang & Gurfil Eq. 12:
c     a(1-e^2)^(4/7) vs a computable critical value from J2/i/solar
c     mean motion) -- see scratch_rpe/resonance_campaign_sweep.py and
c     scratch_rpe/resonance_sweep_97obj.csv for the detection sweep and
c     full per-object results.
c
c     HARD RULE (explicit user requirement): each object's TLE input
c     here is NOT the object's full tracked history -- it is
c     input/example_<norad>_postresonance.tle.txt, truncated by
c     scratch_rpe/make_postresonance_tles.py to only the TLE records at
c     or after that object's own detected resonance epoch. zone_select/
c     the GA never sees any pre-resonance data for these 16 objects.
c     This follows directly from the paper's own conclusion: dynamics
c     are only reliably predictable BEFORE the resonance (a different,
c     not-directly-relevant regime for terminal re-entry timing) and
c     highly sensitive/chaotic AFTER it, which per the paper is also
c     the regime that is "the prelude to the final re-entry of GTO" --
c     i.e. the only regime this campaign is trying to characterize.
c
c     4 of the 16 have ZERO post-resonance TLEs at all (resonance was
c     detected at/coincident with the object's very last tracked TLE:
c     60328, 28141, 4570, 44911) -- these are included anyway so the
c     pipeline's own ierr/zero-zone response to "no data survives the
c     hard rule" is a real, reproducible result, not silently dropped.
c     3 more have 1-7 post-resonance TLEs (19881, 5980, 28140) -- below
c     zone_select's min_pts=8, so no zone can form; also run anyway for
c     the same reason. The remaining 9 have >=8 post-resonance points
c     and can form at least one real post-resonance zone.
      program rpe_campaign_resonance16

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      integer io_start, io_end, iarg
      character*80 csv_name
      character*16 arg_str
      integer irefine_cli

      parameter(ndim_a = 291, mxz = 8, nobj = 16)

      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr, nzv
      integer zstat(mxz)
      double precision bnb_t(10), tb_t(10)
      integer nbv_t

      character*54 tle_file(nobj)
      character*24 obj_name(nobj)
      integer norad(nobj)
      double precision oyr(nobj), omo(nobj), ody(nobj)
      double precision t_obs_jd, horiz, rpe_ens
      integer io, iz, nre

      data norad / 35497, 37151, 27526, 37819, 59347, 40943, 60328,
     &             18983, 13112, 28141, 19881, 5980, 4570, 10724,
     &             28140, 44911 /
      data oyr / 2016.d0, 2015.d0, 2012.d0, 2013.d0, 2026.d0, 2020.d0,
     &           2025.d0, 2002.d0, 1992.d0, 2007.d0, 1993.d0, 1976.d0,
     &           1972.d0, 1982.d0, 2007.d0, 2022.d0 /
      data omo / 10.d0, 12.d0, 5.d0, 9.d0, 6.d0, 2.d0,
     &           4.d0, 1.d0, 11.d0, 10.d0, 3.d0, 4.d0,
     &           5.d0, 8.d0, 10.d0, 5.d0 /
      data ody / 31.d0, 3.d0, 9.d0, 12.d0, 8.d0, 20.d0,
     &           14.d0, 23.d0, 14.d0, 14.d0, 1.d0, 12.d0,
     &           31.d0, 31.d0, 14.d0, 4.d0 /

      tle_file(1) = 'input/example_35497_postresonance.tle.txt'
      tle_file(2) = 'input/example_37151_postresonance.tle.txt'
      tle_file(3) = 'input/example_27526_postresonance.tle.txt'
      tle_file(4) = 'input/example_37819_postresonance.tle.txt'
      tle_file(5) = 'input/example_59347_postresonance.tle.txt'
      tle_file(6) = 'input/example_40943_postresonance.tle.txt'
      tle_file(7) = 'input/example_60328_postresonance.tle.txt'
      tle_file(8) = 'input/example_18983_postresonance.tle.txt'
      tle_file(9) = 'input/example_13112_postresonance.tle.txt'
      tle_file(10) = 'input/example_28141_postresonance.tle.txt'
      tle_file(11) = 'input/example_19881_postresonance.tle.txt'
      tle_file(12) = 'input/example_5980_postresonance.tle.txt'
      tle_file(13) = 'input/example_4570_postresonance.tle.txt'
      tle_file(14) = 'input/example_10724_postresonance.tle.txt'
      tle_file(15) = 'input/example_28140_postresonance.tle.txt'
      tle_file(16) = 'input/example_44911_postresonance.tle.txt'

      obj_name(1) = '35497 Ariane 5 ESC-A'
      obj_name(2) = '37151 Long March 3B'
      obj_name(3) = '27526 Ariane 5 R/B'
      obj_name(4) = '37819 Proton-M R/B'
      obj_name(5) = '59347 GTO R/B'
      obj_name(6) = '40943 GTO debris'
      obj_name(7) = '60328 CZ-3B R/B'
      obj_name(8) = '18983 SL-6 R/B(2)'
      obj_name(9) = '13112 SL-6 R/B(2)'
      obj_name(10) = '28141 CZ-2C R/B'
      obj_name(11) = '19881 Cosmos931 deb'
      obj_name(12) = '5980 Titan3C TS deb'
      obj_name(13) = '4570 SL-6 R/B(2)'
      obj_name(14) = '10724 AtlasCentaur deb'
      obj_name(15) = '28140 Doublestar TC-1'
      obj_name(16) = '44911 CZ-5 R/B'

      io_start = 1
      io_end = nobj
      irefine_cli = 0
      csv_name = 'scratch_rpe/rpe_campaign_resonance16.csv'
      iarg = iargc()
      if (iarg .ge. 2) then
         call getarg(1, arg_str)
         read(arg_str, *) io_start
         call getarg(2, arg_str)
         read(arg_str, *) io_end
         write(csv_name, '(A,I2.2,A,I2.2,A)')
     &      'scratch_rpe/rpe_campaign_resonance16_part_', io_start,
     &      '_', io_end, '.csv'
      end if
      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      open(unit=33, file=csv_name, status='unknown')
      write(33,'(A)') 'norad,zone,e_opt,bn_opt,reentry_jd,rpe_pct,'//
     &   'zstat,t_mean,t_std,ens_rpe_pct,zepoch,rms_fit'

      do io = io_start, io_end
         write(*,'(/,A,A)') ' ===== ', obj_name(io)

         t_obs(1) = oyr(io)
         t_obs(2) = omo(io)
         t_obs(3) = ody(io)
         t_obs(4) = 0.d0
         t_obs(5) = 0.d0
         t_obs(6) = 0.d0

         call orem_run(
     &      tle_file(io), norad(io),
     &      t_obs, mxz, 8, 20.d0, 0.90d0, -1.d0,
     &      80.d0, 160.d0, 0, 1, 0, irefine_cli,
     &      20, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &      20, 2, 3,
     &      7.2921150d-5, 3.35281066d-3, 1.d0,
     &      1.2d0, 0.01d0, 1, 2,
     &      4.56d-6, 1.32712440018d11, 4.902801076d3,
     &      ALT_a, DEN_a, SCH_a, natm,
     &      reentry, e_opt, bn_opt, rms_o,
     &      zepoch, nzu,
     &      zstat, nzv,
     &      rpe, t_mean, t_std,
     &      bnb_t, tb_t, nbv_t, ierr)

         if (ierr .ne. 0 .and. ierr .ne. 3) then
            write(*,'(A,I3)') '  pipeline error ierr=', ierr
            write(33,'(I5,A,I3)') norad(io), ',ERR,,,,,,,,,,', ierr
            goto 90
         end if

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
     &         A,F13.4,A,F8.2,A,F9.3,A,F13.4,A,F14.4)')
     &         norad(io), ',', iz, ',', e_opt(iz), ',', bn_opt(iz),
     &         ',', reentry(iz), ',', rpe(iz), ',', zstat(iz),
     &         ',', t_mean, ',', t_std, ',', rpe_ens,
     &         ',', zepoch(iz), ',', rms_o(iz)
         end do

         if (nre .ge. 1) then
            write(*,'(A,I1,A,I1,A,F13.3,A,F7.1,A)')
     &         '  ensemble (', nre, ' of ', nzu,
     &         ' zones): mean JD=', t_mean, ' std=', t_std, ' d'
            write(*,'(A,F8.2,A)')
     &         '  ENSEMBLE RPE = ', rpe_ens, '%'
         else
            write(*,'(A)') '  no re-entry in any zone (weak-signal'//
     &         ' fits, or too few post-resonance TLEs to form a'//
     &         ' zone at all -- see header comment)'
         end if

 90      continue
         write(*,'(A,I3,A,I2,A,I2,A)') ' [PROGRESS] ',
     &      ((io-io_start+1)*100)/(io_end-io_start+1), '% (',
     &      io, ' of ', io_end, ' objects done, this chunk)'
      end do

      close(33)
      write(*,'(/,A,A)') ' campaign complete -> ', trim(csv_name)

      stop
      end

c     ============================================================
c     ATM reader (same as rpe_campaign.F / test_reentry.F)
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

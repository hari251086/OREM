c     rpe_campaign_resonance_new7.F -- same hard-rule post-resonance-only
c     campaign as rpe_campaign_resonance16.F, applied to 7 new resonant
c     objects found by expanding the candidate pool beyond the original
c     97 (screen_new_candidates.py + fetch_new_candidates.py +
c     resonance_sweep_new_candidates.py, GitHub issue #56 follow-up).
c
c     3 other new resonant hits (39501/44187/42985) are Falcon 9 R/Bs --
c     excluded per this project's established policy (issue #29
c     follow-up, 2026-07-26 -- OREM has no thrust/maneuver modeling and
c     Falcon 9 second stages perform active post-separation passivation/
c     deorbit burns, unlike the purely passive expendable stages this
c     campaign otherwise uses). 1 more (13913) had an internally
c     inconsistent SATCAT decay date (last tracked TLE ~4.7yr AFTER its
c     catalogued decay, same class of exclusion as 9892 in the original
c     97-object round) -- also excluded.
c
c     TLE inputs are input/example_<norad>_postresonance.tle.txt,
c     truncated to only records at/after each object's own detected
c     resonance epoch -- same hard rule, same orem_run parameters as
c     rpe_campaign.F/rpe_campaign_resonance16.F for comparability.
      program rpe_campaign_resonance_new7

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      integer io_start, io_end, iarg
      character*80 csv_name
      character*16 arg_str
      integer irefine_cli

      parameter(ndim_a = 291, mxz = 8, nobj = 7)

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
      character*30 obj_name(nobj)
      integer norad(nobj)
      double precision oyr(nobj), omo(nobj), ody(nobj)
      double precision t_obs_jd, horiz, rpe_ens
      integer io, iz, nre

      data norad / 7096, 68538, 6294, 18954, 29462, 5941, 6231 /
      data oyr / 1974.d0, 2026.d0, 1976.d0, 1989.d0, 2012.d0, 1980.d0,
     &           1977.d0 /
      data omo / 1.d0, 4.d0, 2.d0, 2.d0, 10.d0, 1.d0, 11.d0 /
      data ody / 25.d0, 10.d0, 11.d0, 1.d0, 3.d0, 19.d0, 1.d0 /

      tle_file(1) = 'input/example_7096_postresonance.tle.txt'
      tle_file(2) = 'input/example_68538_postresonance.tle.txt'
      tle_file(3) = 'input/example_6294_postresonance.tle.txt'
      tle_file(4) = 'input/example_18954_postresonance.tle.txt'
      tle_file(5) = 'input/example_29462_postresonance.tle.txt'
      tle_file(6) = 'input/example_5941_postresonance.tle.txt'
      tle_file(7) = 'input/example_6231_postresonance.tle.txt'

      obj_name(1) = '7096 SKYNET 2A'
      obj_name(2) = '68538 ARTEMIS 2 (INTEGRITY)'
      obj_name(3) = '6294 MOLNIYA 1-22'
      obj_name(4) = '18954 ARIANE 3 DEB (SYLDA)'
      obj_name(5) = '29462 SL-12 DEB'
      obj_name(6) = '5941 PROGNOZ 1'
      obj_name(7) = '6231 MOLNIYA 1-21'

      io_start = 1
      io_end = nobj
      irefine_cli = 0
      csv_name = 'scratch_rpe/rpe_campaign_resonance_new7.csv'
      iarg = iargc()
      if (iarg .ge. 2) then
         call getarg(1, arg_str)
         read(arg_str, *) io_start
         call getarg(2, arg_str)
         read(arg_str, *) io_end
         write(csv_name, '(A,I2.2,A,I2.2,A)')
     &      'scratch_rpe/rpe_campaign_resonance_new7_part_', io_start,
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
     &         ' zone at all)'
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

c     validation_rpe_campaign_v2.F -- corrected-config re-run of the
c     7-object validation set (issue #56 manual review), diagnosing the
c     RPE gap vs COSPAR ASR found in rpe_accuracy_review.md.
c
c     Unlike scratch_rpe/validation_rpe/validation_rpe_campaign.F (the
c     340-zone full-tiling campaign), this driver:
c       - uses each object's FULL tracked TLE history (input/example_
c         <norad>.tle.txt), not a pre-split single-zone file per call
c       - restores PRODUCTION zone_select defaults (min_zone_pts=8,
c         max_zone_days=20, r2_thresh=0.90, slope_thresh=-1.0) instead
c         of the permissive override
c       - uses nzones_max=8 (mxz), so orem_run finds multiple zones
c         itself AND runs its own BN trust-gated carryover + ensemble
c         averaging across them -- the mechanism the 340-zone campaign
c         structurally could not exercise (mxz=1 per call there)
c     Same GA/force-model parameters as the main scratch_rpe/
c     rpe_campaign.F 97-object campaign, for direct comparability.
      program validation_rpe_campaign_v2

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      character*80 csv_name
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
      character*24 obj_name(nobj)
      integer norad(nobj)
      double precision oyr(nobj), omo(nobj), ody(nobj)
      double precision t_obs_jd, horiz, rpe_ens
      integer io, iz, nre

      data norad / 35497, 37151, 27526, 59347, 40943, 10724, 18954 /
      data oyr / 2016.d0, 2015.d0, 2012.d0, 2026.d0, 2020.d0, 1982.d0,
     &           1989.d0 /
      data omo / 10.d0, 12.d0, 5.d0, 6.d0, 2.d0, 8.d0, 2.d0 /
      data ody / 31.d0, 3.d0, 9.d0, 8.d0, 20.d0, 31.d0, 1.d0 /

      tle_file(1) = 'input/example_35497.tle.txt'
      tle_file(2) = 'input/example_37151.tle.txt'
      tle_file(3) = 'input/example_27526.tle.txt'
      tle_file(4) = 'input/example_59347.tle.txt'
      tle_file(5) = 'input/example_40943.tle.txt'
      tle_file(6) = 'input/example_10724.tle.txt'
      tle_file(7) = 'input/example_18954.tle.txt'

      obj_name(1) = '35497 Ariane 5 ESC-A'
      obj_name(2) = '37151 Long March 3B'
      obj_name(3) = '27526 Ariane 5 R/B'
      obj_name(4) = '59347 GTO R/B'
      obj_name(5) = '40943 GTO debris'
      obj_name(6) = '10724 AtlasCentaur deb'
      obj_name(7) = '18954 ARIANE 3 DEB'

      csv_name = 'scratch_rpe/validation_RPE/valrpe_v2_campaign.csv'
      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      open(unit=33, file=csv_name, status='unknown')
      write(33,'(A)') 'norad,zone,e_opt,bn_opt,reentry_jd,rpe_pct,'//
     &   'zstat,t_mean,t_std,ens_rpe_pct,zepoch,rms_fit'

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
     &      t_obs, mxz, 8, 20.d0, 0.90d0, -1.d0,
     &      80.d0, 160.d0, 0, 1, 0, 0,
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
            write(*,'(A)') '  no re-entry in any zone'
         end if

 90      continue
      end do

      close(33)
      write(*,'(/,A,A)') ' campaign complete -> ', trim(csv_name)

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

c     rpe_campaign_7obj_issue41.F -- issue #41 ceiling study: runs the
c     curated 7-object set twice per object, iadapt_zone=0 (baseline,
c     unchanged pipeline) vs iadapt_zone=1 (RPE-driven sequential zone
c     selection). Mode-1 RPE only, so this is meaningful ONLY because
c     these 7 objects have a real known re-entry date -- not a
c     deployable production path. Same objects/config as
c     rpe_campaign_7obj_issue31.F (nzones_max=8, min_zone_pts=8,
c     max_zone_days=20, r2=0.90).
c
c     Build: registered as the "campaign41" executable in fpm.toml.
c       fpm run campaign41 --compiler ifx
c     -------------------------------------------------------
      program rpe_campaign_7obj_issue41

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291, mxz = 8, nobj = 7)

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
      double precision t_obs_jd, horiz
      double precision rpe_ens0
      integer io, iadapt, nre, iz_last

      data norad / 42928, 35497, 37151, 39615, 27526, 32007, 37819 /
      data oyr / 2019.d0, 2016.d0, 2015.d0, 2017.d0, 2012.d0,
     &           2010.d0, 2013.d0 /
      data omo / 3.d0, 10.d0, 12.d0, 9.d0, 5.d0, 6.d0, 9.d0 /
      data ody / 3.d0, 31.d0, 3.d0, 15.d0, 9.d0, 6.d0, 12.d0 /

      tle_file(1) = 'input/example_42928.tle.txt'
      tle_file(2) = 'input/example_35497.tle.txt'
      tle_file(3) = 'input/example_37151.tle.txt'
      tle_file(4) = 'input/example_39615.tle.txt'
      tle_file(5) = 'input/example_27526.tle.txt'
      tle_file(6) = 'input/example_32007.tle.txt'
      tle_file(7) = 'input/example_37819.tle.txt'
      obj_name(1) = '42928 PSLV-C39'
      obj_name(2) = '35497 Ariane 5 ESC-A'
      obj_name(3) = '37151 Long March 3B'
      obj_name(4) = '39615 Proton-M Briz-M'
      obj_name(5) = '27526 Ariane 5 R/B'
      obj_name(6) = '32007 GSLV R/B'
      obj_name(7) = '37819 Proton-M R/B'

      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      open(unit=34, file='scratch_rpe/rpe_campaign_7obj_issue41.csv',
     &     status='unknown')
      write(34,'(A)') 'norad,iadapt_zone,nzu,latest_zone_rpe_pct,'//
     &   'ensemble_rpe_pct'

      do io = 1, nobj
         write(*,'(/,A,A)') ' ===== ', obj_name(io)

         t_obs(1) = oyr(io)
         t_obs(2) = omo(io)
         t_obs(3) = ody(io)
         t_obs(4) = 0.d0
         t_obs(5) = 0.d0
         t_obs(6) = 0.d0
         call cal2jd(t_obs, t_obs_jd)

         do iadapt = 0, 1
            write(*,'(A,I1)') '  --- iadapt_zone=', iadapt

            call orem_run(
     &         tle_file(io), norad(io),
     &         t_obs, mxz, 8, 20.d0, 0.90d0, -1.d0,
     &         80.d0, 160.d0, 1, iadapt, 0,
     &         20, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &         20, 2, 3,
     &         7.2921150d-5, 3.35281066d-3, 1.d0,
     &         1.2d0, 0.01d0, 1, 2,
     &         4.56d-6, 1.32712440018d11, 4.902801076d3,
     &         ALT_a, DEN_a, SCH_a, natm,
     &         reentry, e_opt, bn_opt, rms_o,
     &         zepoch, nzu,
     &         zstat, nzv,
     &         rpe, t_mean, t_std, ierr)

            if (ierr .ne. 0 .and. ierr .ne. 3) then
               write(*,'(A,I3)') '  pipeline error ierr=', ierr
               write(34,'(I5,A,I1,A,I3)') norad(io), ',', iadapt,
     &            ',ERR,,', ierr
               goto 80
            end if

            rpe_ens0 = 0.d0
            nre = 0
            iz_last = 0
            do iz = 1, nzu
               if (reentry(iz) .gt. 0.d0) then
                  nre = nre + 1
                  iz_last = iz
               end if
            end do
            if (nre .ge. 1 .and. t_obs_jd .gt. 0.d0) then
               horiz = t_obs_jd - zepoch(1)
               if (horiz .gt. 0.d0)
     &            rpe_ens0 = (t_mean - t_obs_jd)/horiz*100.d0
            end if

            if (iz_last .gt. 0) then
               write(*,'(A,I2,A,F8.2,A)') '  latest zone=', iz_last,
     &            ' RPE=', rpe(iz_last), '%'
            else
               write(*,'(A)') '  no zone predicted re-entry'
            end if
            write(*,'(A,F8.2,A)') '  ensemble RPE=', rpe_ens0, '%'

            if (iz_last .gt. 0) then
               write(34,'(I5,A,I1,A,I2,A,F9.3,A,F9.3)')
     &            norad(io), ',', iadapt, ',', nzu, ',',
     &            rpe(iz_last), ',', rpe_ens0
            else
               write(34,'(I5,A,I1,A,I2,A,A)')
     &            norad(io), ',', iadapt, ',', nzu,
     &            ',,'
            end if

 80         continue
         end do

         write(*,'(A,I3,A,I2,A,I2,A)') ' [PROGRESS] ',
     &      (io*100)/nobj, '% (', io, ' of ', nobj, ' objects done)'
      end do

      close(34)
      write(*,'(/,A)') ' campaign complete -> scratch_rpe/'//
     &   'rpe_campaign_7obj_issue41.csv'

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

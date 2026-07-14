c     ============================================================
c     report.F — Structured output reporting for OREM (issue #13)
c
c     Writes a human-readable prediction report from orem_run's
c     outputs: per-zone fit/prediction table, ensemble statistics
c     (mean/std/relative spread across zones with a re-entry), and
c     the zone_status diagnostic legend. The ensemble row is the
c     headline: as of v1.18 the ensemble mean beats every individual
c     zone on the validation set (42928: ensemble +11% vs best
c     single zone +3.2%, worst +52%).
c
c     Author: Harishkumar Sellamuthu
c     ============================================================
      subroutine orem_report(
     &   rep_file, norad_id, t_obs_cal,
     &   bn_min_init, bn_max_init, idrag_flag,
     &   reentry_jd, e_opt, bn_opt, rms_o,
     &   zone_epoch, nzones_used,
     &   zone_status, nzones_valid,
     &   rpe, t_mean, t_std, nzmax, ierr_rep)
c
c     Input
c       rep_file        : output file path
c       norad_id        : NORAD catalog number
c       t_obs_cal(6)    : observed re-entry [yr,mo,dy,hr,mn,sc]
c                         (all zeros = operational mode, no RPE)
c       bn_min_init, bn_max_init, idrag_flag : config echo
c       reentry_jd(nzmax), e_opt(nzmax), bn_opt(nzmax),
c       rms_o(nzmax), zone_epoch(nzmax), nzones_used,
c       zone_status(nzmax), nzones_valid,
c       rpe(nzmax), t_mean, t_std : orem_run outputs
c     Output
c       ierr_rep        : 0=ok, 1=cannot open rep_file
c     ============================================================
      implicit double precision (a-h, o-z)

      character*(*) rep_file
      integer norad_id
      double precision t_obs_cal(6)
      double precision bn_min_init, bn_max_init
      integer idrag_flag
      integer nzmax
      double precision reentry_jd(nzmax), e_opt(nzmax)
      double precision bn_opt(nzmax), rms_o(nzmax)
      double precision zone_epoch(nzmax)
      integer nzones_used
      integer zone_status(nzmax), nzones_valid
      double precision rpe(nzmax), t_mean, t_std
      integer ierr_rep

      character*10 zslab(0:4)
      character*20 datestr
      double precision cal(6), t_obs_jd, horiz, rpe_ens, rel_spread
      integer iz, nre, ivals(8)

      zslab(0) = 'ok        '
      zslab(1) = 'propfail  '
      zslab(2) = 'boundary  '
      zslab(3) = 'nobound   '
      zslab(4) = 'toofewpts '

      ierr_rep = 0
      open(unit=44, file=rep_file, status='unknown', err=900)

      call date_and_time(values=ivals)
      write(datestr,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')
     &   ivals(1),'-',ivals(2),'-',ivals(3),' ',ivals(5),':',ivals(6)

      write(44,'(A)') 'OREM Re-entry Prediction Report'
      write(44,'(A)') '================================'
      write(44,'(A,I6)')      'Object (NORAD): ', norad_id
      write(44,'(A,A)')       'Generated:      ', datestr
      write(44,'(A,F7.1,A,F7.1,A)')
     &   'BN search:      [', bn_min_init, ',', bn_max_init,
     &   '] kg/m2 (initial; G2 floor may extend zone 1 downward)'
      write(44,'(A,I2)')      'IDRAG:          ', idrag_flag
      write(44,'(A,I2,A,I2,A)') 'Zones:          ', nzones_used,
     &   ' processed, ', nzones_valid, ' clean (zone_status=0)'

      t_obs_jd = 0.d0
      if (t_obs_cal(1) .gt. 1900.d0) then
         call cal2jd(t_obs_cal, t_obs_jd)
         call jd2cal(t_obs_jd, cal)
         write(44,'(A,F13.4,A,I5,I3,I3,A)')
     &      'Observed:       JD ', t_obs_jd, '  (',
     &      int(cal(1)), int(cal(2)), int(cal(3)), ' )'
      else
         write(44,'(A)') 'Observed:       none (operational mode)'
      end if

      write(44,'(A)') ''
      write(44,'(A)') 'Zone  Epoch (JD)      e_opt     BN(kg/m2)'//
     &   '  Re-entry (JD)   Re-entry (UTC)   RPE(%)    status'
      do iz = 1, nzones_used
         if (reentry_jd(iz) .gt. 0.d0) then
            call jd2cal(reentry_jd(iz), cal)
            write(44,'(I4,F15.4,F10.5,F10.2,F15.4,I7,I3,I3,F11.2,
     &         3X,A)')
     &         iz, zone_epoch(iz), e_opt(iz), bn_opt(iz),
     &         reentry_jd(iz),
     &         int(cal(1)), int(cal(2)), int(cal(3)),
     &         rpe(iz), zslab(zone_status(iz))
         else
            write(44,'(I4,F15.4,F10.5,F10.2,A,3X,A)')
     &         iz, zone_epoch(iz), e_opt(iz), bn_opt(iz),
     &         '      no re-entry within cap        --',
     &         zslab(zone_status(iz))
         end if
      end do

c     Ensemble block
      nre = 0
      do iz = 1, nzones_used
         if (reentry_jd(iz) .gt. 0.d0) nre = nre + 1
      end do

      write(44,'(A)') ''
      if (nre .ge. 1) then
         call jd2cal(t_mean, cal)
         write(44,'(A,I2,A,I2,A)')
     &      'Ensemble (', nre, ' of ', nzones_used,
     &      ' zones with a predicted re-entry):'
         write(44,'(A,F13.4,A,I5,I3,I3,A,F8.2,A)')
     &      '  mean re-entry JD ', t_mean, '  (',
     &      int(cal(1)), int(cal(2)), int(cal(3)),
     &      ' )  +- ', t_std, ' days'
         rel_spread = 0.d0
         if (t_mean .gt. zone_epoch(1) .and. nre .ge. 2) then
            rel_spread = t_std/(t_mean - zone_epoch(1))*100.d0
            write(44,'(A,F7.2,A)')
     &         '  relative spread: ', rel_spread,
     &         ' % of the zone-1 -> mean horizon'
         end if
         if (t_obs_jd .gt. 0.d0) then
            horiz = t_obs_jd - zone_epoch(1)
            if (horiz .gt. 0.d0) then
               rpe_ens = (t_mean - t_obs_jd)/horiz*100.d0
               write(44,'(A,F8.2,A)')
     &            '  ENSEMBLE RPE:    ', rpe_ens, ' %'
            end if
         end if
      else
         write(44,'(A)') 'Ensemble: no zone predicted a re-entry'//
     &      ' within the propagation cap'
      end if

      write(44,'(A)') ''
      write(44,'(A)') 'zone_status legend: 0=ok, 1=skip_propfail,'//
     &   ' 2=GA_boundary, 3=RSM_nobound, 4=skip_toofewpts'
      close(44)
      return

 900  ierr_rep = 1
      return
      end

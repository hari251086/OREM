c     ============================================================
c     diag_osc_vs_mean_apogee_v2.F -- issue #36, Finding 1, proper
c     within-revolution phase averaging
c
c     v1 (diag_osc_vs_mean_apogee.F) mistakenly assumed propagate_ks
c     dumps 360 pts/rev and averaged across nearby REVOLUTIONS at the
c     SAME fixed anomaly phase each time -- which changed almost
c     nothing, because same-phase samples share nearly the same
c     phase-dependent short-period bias rather than averaging it out.
c     This version uses propagate_ks_dense (scratch_rpe/
c     propagate_ks_densedump.F, a diagnostic-only fork that dumps
c     EVERY integration step instead of once per revolution) to get
c     genuine multi-phase sampling within a single revolution, then
c     boxcar-averages across a full period -- an actual numerical
c     mean-element extraction, not same-phase resampling.
c
c     Read-only diagnostic (propagate_ks_dense is diagnostic-only,
c     never linked into any production executable). No production
c     code touched.
c     ============================================================
      program diag_osc_vs_mean_apogee_v2
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(maxpts = 15000, maxtle = 100000)
      character*80 tle_file
      integer norad_id

      double precision epochs(maxpts), a_tle(maxpts)
      double precision e_tle(maxpts), inc_tle(maxpts)
      double precision raan_tle(maxpts), aop_tle(maxpts)
      double precision ha_tle(maxpts), hp_tle(maxpts)
      double precision ls_tle(maxpts)
      integer npts, ierr_tle

      integer mask_tle(maxpts)
      integer nfiltered_tle, nremoved_tle, ngaps_tle
      integer i_filt, j_filt

      character*24 names_raw(maxtle)
      character*5 catid_raw(maxtle)
      integer norad_raw(maxtle)
      integer epoch_yr_raw(maxtle), revnum_raw(maxtle)
      double precision epoch_day_raw(maxtle)
      double precision ndot_raw(maxtle), nddot_raw(maxtle)
      double precision bstar_raw(maxtle)
      double precision ainc_raw(maxtle), raan_raw(maxtle)
      double precision ecc_raw(maxtle), aop_raw(maxtle)
      double precision ama_raw(maxtle), amm_raw(maxtle)
      integer ntle_raw

      integer nzones_max
      parameter(nzones_max = 8)
      integer zone_start(nzones_max), zone_end(nzones_max)
      integer nzones

      double precision epz(500), haz(500), ez(500)
      double precision incz(500), raanz(500), aopz(500)
      double precision maz(500), smaz(500)

c     Dense dump: up to istep(360) pts/rev * ~200 revs -- generous
c     headroom over the ~65,520 points actually expected for 42928
c     zone 1 (182 revs).
      parameter(max_traj = 150000)
      double precision traj_jd(max_traj)
      double precision traj_x(3, max_traj)
      double precision traj_xd(3, max_traj)
      double precision osc_ha(max_traj), avg_ha(max_traj)
      integer exit_code_re

      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      integer n_force(3)
      double precision alt_a(291), den_a(291), sch_a(291)
      integer natm

      double precision e_opt_zone1, bn_opt_zone1
      integer iz, i, nzone1, ntraj, ierr_sw, ierr_a2, idump_loc
      integer istep_val

      call init_constants()

      tle_file = 'input/example_42928.tle.txt'
      norad_id = 42928
      e_opt_zone1  = 0.321707d0
      bn_opt_zone1 = 81.195d0
      istep_val = 360

      write(*,'(A)') 'Step 1: tle_evolve'
      call tle_evolve(tle_file, norad_id, maxpts,
     &   epochs, a_tle, e_tle, inc_tle,
     &   raan_tle, aop_tle,
     &   ha_tle, hp_tle, ls_tle,
     &   npts, ierr_tle)
      if (ierr_tle .ne. 0) then
         write(*,*) 'tle_evolve failed, ierr=', ierr_tle
         stop 1
      end if

      write(*,'(A)') 'Step 2: tle_filter (production params)'
      call tle_filter(epochs, ha_tle, e_tle, npts,
     &   4.0d0, 6.0d0, 30.0d0, 20,
     &   mask_tle, nfiltered_tle, nremoved_tle, ngaps_tle)
      if (nremoved_tle .gt. 0) then
         j_filt = 0
         do i_filt = 1, npts
            if (mask_tle(i_filt) .eq. 1) then
               j_filt = j_filt + 1
               epochs(j_filt)   = epochs(i_filt)
               a_tle(j_filt)    = a_tle(i_filt)
               e_tle(j_filt)    = e_tle(i_filt)
               inc_tle(j_filt)  = inc_tle(i_filt)
               raan_tle(j_filt) = raan_tle(i_filt)
               aop_tle(j_filt)  = aop_tle(i_filt)
               ha_tle(j_filt)   = ha_tle(i_filt)
               hp_tle(j_filt)   = hp_tle(i_filt)
               ls_tle(j_filt)   = ls_tle(i_filt)
            end if
         end do
         npts = j_filt
      end if
      write(*,'(A,I6)') '  npts after filter=', npts

      write(*,'(A)') 'Step 3: read raw TLE for mean anomaly'
      open(unit=77, file=tle_file, status='old')
      call read_tle(77, maxtle, names_raw, catid_raw,
     &   norad_raw, epoch_yr_raw, epoch_day_raw,
     &   ndot_raw, nddot_raw, bstar_raw,
     &   ainc_raw, raan_raw, ecc_raw,
     &   aop_raw, ama_raw, amm_raw,
     &   revnum_raw, ntle_raw)
      close(77)

      write(*,'(A)') 'Step 4: zone_select (production params)'
      call zone_select(epochs, ha_tle, npts,
     &   8, 20.0d0, 0.90d0, -1.0d0,
     &   nzones_max,
     &   zone_start, zone_end, nzones)
      write(*,'(A,I3)') '  nzones=', nzones
      if (nzones .lt. 1) stop 1

      nzone1 = zone_end(1) - zone_start(1) + 1
      write(*,'(A,I4,A,F13.4,A,F13.4)') '  zone1 npts=', nzone1,
     &   ' epoch(1)=', epochs(zone_start(1)),
     &   ' epoch(end)=', epochs(zone_end(1))
      do i = 1, nzone1
         epz(i)   = epochs(zone_start(1) + i - 1)
         haz(i)   = ha_tle(zone_start(1) + i - 1)
         ez(i)    = e_tle(zone_start(1) + i - 1)
         incz(i)  = inc_tle(zone_start(1) + i - 1)
         raanz(i) = raan_tle(zone_start(1) + i - 1)
         aopz(i)  = aop_tle(zone_start(1) + i - 1)
         smaz(i)  = a_tle(zone_start(1) + i - 1)
         call tle_find_ma(norad_id, epz(i),
     &        norad_raw, epoch_yr_raw, epoch_day_raw,
     &        ama_raw, ntle_raw, maxtle, maz(i))
      end do

      a_kep(1) = smaz(1)
      a_kep(2) = e_opt_zone1
      a_kep(3) = incz(1)
      a_kep(4) = raanz(1)
      a_kep(5) = aopz(1)
      a_kep(6) = maz(1)
      a_kep(7) = 0.d0
      rp = smaz(1) * (1.d0 - ez(1))
      a_kep(1) = rp / (1.d0 - e_opt_zone1)

      call oe2car(a_kep, x0, xd0, 1.d-12)
      call jd2cal(epz(1), cal0)

      period_s = 2.d0*pi*dsqrt(a_kep(1)**3/amue)
      period_d = period_s / 86400.d0
      zone_dur = epz(nzone1) - epz(1)
      nrev_zone = int(zone_dur / period_d) + 2
      write(*,'(A,F10.5,A,F10.2,A,I5)') '  period_d=', period_d,
     &   ' zone_dur=', zone_dur, ' nrev=', nrev_zone
      write(*,'(A,I8)') '  expected dense dumps ~',
     &   nrev_zone*istep_val

      call read_atm_local(alt_a, den_a, sch_a, natm)
      call sw_load('input/SW-All.csv', ierr_sw)
      call atm2d_load('input/ATM2D.DAT', ierr_a2)
      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

      write(*,'(A)') 'Step 5: propagate_ks_dense (every-step dump)'
      idump_loc = 1
      do i = 1, max_traj
         traj_jd(i) = 0.d0
      end do
      call propagate_ks_dense(
     &   x0, xd0, cal0,
     &   nrev_zone, istep_val, 1.d-10,
     &   n_force, 20, 2, 3,
     &   bn_opt_zone1, 1, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   alt_a, den_a, sch_a, natm,
     &   max_traj, idump_loc,
     &   traj_jd, traj_x, traj_xd,
     &   exit_code_re)
      write(*,'(A,I3,A,I8)') '  exit_code=', exit_code_re,
     &   ' idump=', idump_loc

      ntraj = 0
      do i = 1, max_traj
         if (traj_jd(i) .le. 0.d0) goto 16
         if (traj_jd(i) - epz(1) .gt. zone_dur + 0.5d0) goto 16
         ntraj = ntraj + 1
         r_mag = dsqrt(traj_x(1,i)**2+traj_x(2,i)**2+traj_x(3,i)**2)
         v_mag = dsqrt(traj_xd(1,i)**2+traj_xd(2,i)**2
     &                +traj_xd(3,i)**2)
         sma_pt = 1.d0/(2.d0/r_mag - v_mag**2/amue)
         hx = traj_x(2,i)*traj_xd(3,i) - traj_x(3,i)*traj_xd(2,i)
         hy = traj_x(3,i)*traj_xd(1,i) - traj_x(1,i)*traj_xd(3,i)
         hz = traj_x(1,i)*traj_xd(2,i) - traj_x(2,i)*traj_xd(1,i)
         h_mag = dsqrt(hx**2+hy**2+hz**2)
         p_sl = h_mag**2/amue
         ecc_pt = dsqrt(dabs(1.d0 - p_sl/sma_pt))
         osc_ha(ntraj) = sma_pt*(1.d0+ecc_pt) - R_Earth
      end do
 16   continue
      write(*,'(A,I8)') '  ntraj=', ntraj

c     --- Genuine within-revolution boxcar average: +/-(istep/2)
c         INDEX positions now span +/- half a period across DIFFERENT
c         anomaly phases (dense dump), not the same phase across
c         different revolutions as v1 mistakenly did. This is what
c         actually tests "does removing short-period phase dependence
c         close the gap." ---
      nhalf = istep_val / 2
      do i = 1, ntraj
         ilo = i - nhalf
         ihi = i + nhalf
         if (ilo .lt. 1) ilo = 1
         if (ihi .gt. ntraj) ihi = ntraj
         ssum = 0.d0
         ncnt = 0
         do k = ilo, ihi
            ssum = ssum + osc_ha(k)
            ncnt = ncnt + 1
         end do
         avg_ha(i) = ssum / dble(ncnt)
      end do

      write(*,'(A)') ''
      write(*,'(A)')
     &  'epoch_day  ha_TLE(truth)  osc_pred   avg_pred'//
     &  '   err_osc   err_avg'
      ssq_osc = 0.d0
      ssq_avg = 0.d0
      amax_osc = 0.d0
      amax_avg = 0.d0
      do i = 1, nzone1
         t_want = epz(i) - epz(1)
         call rsm_tinterp(traj_jd, osc_ha, ntraj, epz(1), t_want,
     &                      y_osc)
         call rsm_tinterp(traj_jd, avg_ha, ntraj, epz(1), t_want,
     &                      y_avg)
         err_osc = y_osc - haz(i)
         err_avg = y_avg - haz(i)
         ssq_osc = ssq_osc + err_osc**2
         ssq_avg = ssq_avg + err_avg**2
         if (dabs(err_osc) .gt. amax_osc) amax_osc = dabs(err_osc)
         if (dabs(err_avg) .gt. amax_avg) amax_avg = dabs(err_avg)
         write(*,'(F9.4,F14.4,F11.4,F11.4,F10.4,F10.4)')
     &      t_want, haz(i), y_osc, y_avg, err_osc, err_avg
      end do

      rms_osc = dsqrt(ssq_osc/dble(nzone1))
      rms_avg = dsqrt(ssq_avg/dble(nzone1))
      write(*,'(A)') ''
      write(*,'(A,F10.4,A,F10.4,A)') 'RMS vs TLE ha: osculating=',
     &   rms_osc, ' km   within-rev-averaged=', rms_avg, ' km'
      write(*,'(A,F10.4,A,F10.4,A)') 'MAX vs TLE ha: osculating=',
     &   amax_osc, ' km   within-rev-averaged=', amax_avg, ' km'

      stop
      end

c     ============================================================
      subroutine read_atm_local(alt, den, sch, natm)
      implicit double precision (a-h, o-z)
      parameter(ndim_a = 291)
      double precision alt(ndim_a), den(ndim_a), sch(ndim_a)
      integer natm, i
      natm = ndim_a
      alt(1) = 60.d0
      do i = 2, 141
         alt(i) = alt(i-1) + 1.d0
      end do
      do i = 142, natm
         alt(i) = alt(i-1) + 2.d0
      end do
      open(unit=60, file='input/ATM.DAT', status='old')
      read(60,110) (sch(i), i=1,61)
      read(60,110) (sch(i), i=62,natm)
      read(60,111) (den(i), i=1,61)
      read(60,111) (den(i), i=62,natm)
      close(60)
 110  format(8f10.3)
 111  format(8e10.3)
      do i = 1, natm
         alt(i) = alt(i) * 1000.d0
         sch(i) = sch(i) * 1000.d0
         den(i) = den(i) * 1.0d10
      end do
      return
      end

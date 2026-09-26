c     manuscript23_leadtime.F -- manuscript (issues #45/#48) lead-time
c     hindcast campaign for the 23 solar-apsidal-resonance objects
c     (issue #56). Each run is one (object, lead time) pair whose TLE
c     file was truncated by prepare_leadtime_runs.py to epochs before
c     (SATCAT decay - lead days) -- no data after the cutoff is used.
c
c     orem_run parameters are IDENTICAL to scratch_rpe/rpe_campaign.F
c     (OREM's validated configuration: nzones_max=8, min_zone_pts=8,
c     max_zone_days=20, R2>=0.90, slope<-1 km/day, BN [80,160] with the
c     G2 floor / G3 BSTAR prior, drag on, full force model, GA 20/500)
c     -- one documented configuration for every result in the paper.
c
c     Runs are read at RUN time from scratch_rpe/manuscript23/
c     runs_manifest.txt (first line = number of runs; then per run one
c     line "norad yr mo dy lead" and one line with the TLE path), so the
c     driver never needs regenerating when the run list changes.
c
c     Usage, from the OREM repo root (input/ATM.DAT is opened relative
c     to it):
c       fpm run manuscript23 --compiler ifx --flag "/heap-arrays"
c     or split across at most 4 separate processes (GitHub\CLAUDE.md
c     #1 -- orem.F uses SAVE'd arrays, so processes, never threads):
c       manuscript23.exe 1 21
c       manuscript23.exe 22 42  ... etc.
c     Output: scratch_rpe/manuscript23/results_summary[_part_s_e].csv
c             scratch_rpe/manuscript23/results_zones[_part_s_e].csv
      program manuscript23_leadtime

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291, mxz = 8, maxrun = 200)

      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr, nzv
      integer zstat(mxz)
      double precision bnb_t(10), tb_t(10)
      integer nbv_t

      integer nrun, io, iz, nre, izp, io_start, io_end, iarg
      integer norad(maxrun), lead(maxrun)
      integer iyr(maxrun), imo(maxrun), idy(maxrun)
      character*120 tle_file(maxrun)
      character*100 sum_name, zon_name
      character*16 arg_str
      double precision t_obs_jd, horiz, rpe_ens, rpe_pri, t_pri

      open(unit=21, file='scratch_rpe/manuscript23/runs_manifest.txt',
     &     status='old')
      read(21, *) nrun
      if (nrun .gt. maxrun) stop 'raise maxrun'
      do io = 1, nrun
         read(21, *) norad(io), iyr(io), imo(io), idy(io), lead(io)
         read(21, '(A)') tle_file(io)
      end do
      close(21)

      io_start = 1
      io_end = nrun
      sum_name = 'scratch_rpe/manuscript23/results_summary.csv'
      zon_name = 'scratch_rpe/manuscript23/results_zones.csv'
      iarg = iargc()
      if (iarg .ge. 2) then
         call getarg(1, arg_str)
         read(arg_str, *) io_start
         call getarg(2, arg_str)
         read(arg_str, *) io_end
         if (io_end .gt. nrun) io_end = nrun
         write(sum_name, '(A,I3.3,A,I3.3,A)')
     &      'scratch_rpe/manuscript23/results_summary_part_',
     &      io_start, '_', io_end, '.csv'
         write(zon_name, '(A,I3.3,A,I3.3,A)')
     &      'scratch_rpe/manuscript23/results_zones_part_',
     &      io_start, '_', io_end, '.csv'
      end if

      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      open(unit=33, file=sum_name, status='unknown')
      write(33,'(A)') 'run_idx,norad,lead_days,ierr,nzones,nvalid,'//
     &   'npred,primary_zone,primary_reentry_jd,primary_rpe_pct,'//
     &   't_mean,t_std,ens_rpe_pct,zepoch1,decay_jd'
      open(unit=34, file=zon_name, status='unknown')
      write(34,'(A)') 'run_idx,norad,lead_days,zone,e_opt,bn_opt,'//
     &   'reentry_jd,rpe_pct,zstat,zepoch,rms_fit'

      do io = io_start, io_end
         write(*,'(/,A,I3,A,I6,A,I4,A)') ' ===== run ', io, ': NORAD ',
     &      norad(io), '  lead ', lead(io), ' d'

         t_obs(1) = dble(iyr(io))
         t_obs(2) = dble(imo(io))
         t_obs(3) = dble(idy(io))
         t_obs(4) = 0.d0
         t_obs(5) = 0.d0
         t_obs(6) = 0.d0
         call cal2jd(t_obs, t_obs_jd)

         call orem_run(
     &      trim(tle_file(io)), norad(io),
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
            write(33,'(I4,A,I6,A,I4,A,I3,A,F13.4)') io, ',', norad(io),
     &         ',', lead(io), ',', ierr, ',,,,,,,,,,,', t_obs_jd
            goto 90
         end if

c        PRIMARY = latest zone with a predicted re-entry (report.F's
c        own rule); ensemble = orem_run's own t_mean over those zones.
         nre = 0
         izp = 0
         do iz = 1, nzu
            if (reentry(iz) .gt. 0.d0) then
               nre = nre + 1
               izp = iz
            end if
         end do
         rpe_ens = 0.d0
         rpe_pri = 0.d0
         t_pri = 0.d0
         horiz = t_obs_jd - zepoch(1)
         if (nre .ge. 1 .and. horiz .gt. 0.d0) then
            rpe_ens = (t_mean - t_obs_jd) / horiz * 100.d0
            t_pri = reentry(izp)
            rpe_pri = rpe(izp)
         end if

         do iz = 1, nzu
            write(34,'(I4,A,I6,A,I4,A,I2,A,F9.6,A,F8.3,A,F13.4,A,
     &         F10.3,A,I2,A,F13.4,A,F14.4)')
     &         io, ',', norad(io), ',', lead(io), ',', iz, ',',
     &         e_opt(iz), ',', bn_opt(iz), ',', reentry(iz), ',',
     &         rpe(iz), ',', zstat(iz), ',', zepoch(iz), ',', rms_o(iz)
         end do

         write(33,'(I4,A,I6,A,I4,A,I3,A,I2,A,I2,A,I2,A,I2,A,F13.4,
     &      A,F10.3,A,F13.4,A,F9.2,A,F10.3,A,F13.4,A,F13.4)')
     &      io, ',', norad(io), ',', lead(io), ',', ierr, ',', nzu,
     &      ',', nzv, ',', nre, ',', izp, ',', t_pri, ',', rpe_pri,
     &      ',', t_mean, ',', t_std, ',', rpe_ens, ',', zepoch(1),
     &      ',', t_obs_jd

         if (nre .ge. 1) then
            write(*,'(A,I2,A,F9.2,A,F9.2,A)') '  PRIMARY Z', izp,
     &         ' RPE=', rpe_pri, '%   ENSEMBLE RPE=', rpe_ens, '%'
         else
            write(*,'(A)') '  no zone predicted a re-entry'
         end if

 90      continue
         write(*,'(A,I3,A,I3,A,I3)') ' [PROGRESS] run ', io,
     &      ' done (chunk ', io_start, '-', io_end
      end do

      close(33)
      close(34)
      write(*,'(/,A,A)') ' campaign complete -> ', trim(sum_name)
      stop
      end

c     ============================================================
c     ATM reader (same as rpe_campaign.F)
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

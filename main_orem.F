c     main_orem.F — Standalone OREM runner
c     Reads orem.cfg and runs the full pipeline.
c
c     Usage: orem.exe input/orem.cfg
c     ============================================================

      program main_orem
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291, mxz = 10)

c     Config
      character*200 cfg_file, tle_file
      integer norad_id
      double precision t_obs_cal(6)
      integer nzones_max, min_zone_pts
      double precision max_zone_days, r2_thresh, slope_thresh
      double precision a_min, a_max, fixed_mass, cd_drag
      integer idrag_flag
      integer ipopsize, maxgen, nbits_e, nbits_a
      double precision pcross, pmute, ga_seed
      integer ngeo_deg, nsun_deg, nmoon_deg
      double precision WE_rot, EPS_f, FR_rot
      integer IPSR, ISHAD
      double precision CR_srp, AM_srp, PSR_srp
      double precision amuS, amuM

c     Atmosphere
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

c     Output
      double precision reentry(mxz), e_opt(mxz), a_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_mean, t_std
      integer nzu, ierr

      double precision cal_out(6)
      integer iarg

c     ============================================================
c     Parse command line argument
c     ============================================================
      iarg = iargc()
      if (iarg .lt. 1) then
         write(*,*) 'Usage: orem.exe <config_file>'
         write(*,*) 'Example: orem.exe input/orem.cfg'
         stop 1
      end if
      call getarg(1, cfg_file)

c     ============================================================
c     Initialize constants
c     ============================================================
      call init_constants()

c     ============================================================
c     Read config file
c     ============================================================
      open(unit=10, file=cfg_file, status='old', err=900)
      read(10,'(A)') tle_file
      read(10,*) norad_id
      read(10,*) t_obs_cal(1), t_obs_cal(2), t_obs_cal(3),
     &           t_obs_cal(4), t_obs_cal(5), t_obs_cal(6)
      read(10,*) nzones_max
      read(10,*) min_zone_pts, max_zone_days,
     &           r2_thresh, slope_thresh
      read(10,*) a_min, a_max
      read(10,*) fixed_mass
      read(10,*) cd_drag
      read(10,*) ipopsize, maxgen, nbits_e, nbits_a,
     &           pcross, pmute, ga_seed
      read(10,*) ngeo_deg, nsun_deg, nmoon_deg
      read(10,*) idrag_flag, WE_rot, EPS_f, FR_rot
      read(10,*) IPSR, CR_srp, AM_srp, ISHAD
      close(10)

      if (nzones_max .gt. mxz) nzones_max = mxz

c     Sun/Moon GM
      amuS = 1.32712440018d11
      amuM = 4.902801076d3
      PSR_srp = 4.56d-6

c     ============================================================
c     Read atmosphere table
c     ============================================================
      natm = ndim_a
      ALT_a(1) = 60.d0
      do i = 2, 141
         ALT_a(i) = ALT_a(i-1) + 1.d0
      end do
      do i = 142, natm
         ALT_a(i) = ALT_a(i-1) + 2.d0
      end do
      open(unit=60, file='input/ATM.DAT', status='old',
     &     err=910)
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

c     ============================================================
c     Print config
c     ============================================================
      write(*,*) '=============================================='
      write(*,*) ' OREM — Optimal Regularized re-Entry Method'
      write(*,*) '=============================================='
      write(*,*) ''
      write(*,'(A,A)') '  TLE file:  ', trim(tle_file)
      write(*,'(A,I8)') '  NORAD ID:  ', norad_id
      if (t_obs_cal(1) .gt. 1900.d0) then
         write(*,'(A,I4,A,I2,A,I2)')
     &      '  Observed:  ',
     &      int(t_obs_cal(1)), '-',
     &      int(t_obs_cal(2)), '-',
     &      int(t_obs_cal(3))
      else
         write(*,*) '  Observed:  (operational mode)'
      end if
      write(*,'(A,I3)') '  Max zones: ', nzones_max
      write(*,'(A,F6.2,A,F6.2,A)')
     &   '  Area:      [', a_min, ', ', a_max, '] m^2'
      write(*,'(A,F8.1,A)') '  Mass:      ', fixed_mass, ' kg'
      write(*,'(A,F4.1)') '  Cd:        ', cd_drag
      write(*,'(A,I1)') '  Drag:      ', idrag_flag
      write(*,*) ''

c     ============================================================
c     Run OREM
c     ============================================================
      write(*,*) '  Running pipeline...'
      write(*,*) ''

      call orem_run(
     &   tle_file, norad_id,
     &   t_obs_cal,
     &   nzones_max, min_zone_pts, max_zone_days,
     &   r2_thresh, slope_thresh,
     &   a_min, a_max,
     &   fixed_mass, cd_drag, idrag_flag,
     &   ipopsize, maxgen, nbits_e, nbits_a,
     &   pcross, pmute, ga_seed,
     &   ngeo_deg, nsun_deg, nmoon_deg,
     &   WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD,
     &   PSR_srp, amuS, amuM,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   reentry, e_opt, a_opt, rms_o,
     &   zepoch, nzu,
     &   rpe, t_mean, t_std,
     &   ierr)

c     ============================================================
c     Print results
c     ============================================================
      write(*,*) '=============================================='
      write(*,*) ' Results'
      write(*,*) '=============================================='
      write(*,*) ''

      if (ierr .eq. 1) then
         write(*,*) '  ERROR: TLE file not found or no match'
         stop 1
      else if (ierr .eq. 2) then
         write(*,*) '  ERROR: No valid zones found'
         stop 1
      end if

      write(*,'(A,I2)') '  Zones processed: ', nzu
      write(*,*) ''

      write(*,'(A4,A14,A10,A8,A14,A18,A10)')
     &   'Zone', 'Epoch (JD)', 'e_opt', 'A (m2)',
     &   'Re-entry (JD)', 'Re-entry (UTC)', 'RPE(%)'
      write(*,'(A4,A14,A10,A8,A14,A18,A10)')
     &   '----', '--------------', '----------',
     &   '--------', '--------------',
     &   '------------------', '----------'

      do i = 1, nzu
         if (reentry(i) .gt. 0.d0) then
            call jd2cal(reentry(i), cal_out)
            write(*,'(I4,F14.2,F10.6,F8.3,F14.2,
     &         I5,A,I2.2,A,I2.2,F10.2)')
     &         i, zepoch(i), e_opt(i), a_opt(i),
     &         reentry(i),
     &         int(cal_out(1)), '-',
     &         int(cal_out(2)), '-',
     &         int(cal_out(3)),
     &         rpe(i)
         else
            write(*,'(I4,F14.2,F10.6,F8.3,A14,A18,A10)')
     &         i, zepoch(i), e_opt(i), a_opt(i),
     &         '  (no re-entry)', '', ''
         end if
      end do

      write(*,*) ''
      if (t_mean .gt. 0.d0) then
         call jd2cal(t_mean, cal_out)
         write(*,'(A,I4,A,I2.2,A,I2.2)')
     &      '  Mean re-entry: ',
     &      int(cal_out(1)), '-',
     &      int(cal_out(2)), '-',
     &      int(cal_out(3))
         write(*,'(A,F6.2,A)') '  Std deviation: ',
     &      t_std, ' days'
      else
         write(*,*) '  No re-entry predicted (drag off?)'
      end if

      if (ierr .eq. 3) then
         write(*,*) ''
         write(*,*) '  NOTE: No re-entry detected.'
         write(*,*) '  Enable drag (IDRAG=1) for prediction.'
      end if

      write(*,*) ''
      write(*,*) '=============================================='

      stop

 900  write(*,*) 'ERROR: cannot open config file: ', trim(cfg_file)
      stop 1
 910  write(*,*) 'ERROR: cannot open input/ATM.DAT'
      stop 1
      end

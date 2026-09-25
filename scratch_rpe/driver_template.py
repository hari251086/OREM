"""Fortran driver template for the validation RPE campaign -- kept as a
plain format string (not a Fortran INCLUDE) so the generated .F file is
fully self-contained and immune to relative-path breakage if the output
directory ever moves. Used by prepare_validation_zones.py."""

DRIVER_TEMPLATE = """\
c     validation_rpe_campaign.F -- Validation RPE campaign for the 7
c     confirmed/plausible-genuine solar-apsidal-resonance objects from
c     issue #56's manual review (35497, 37151, 27526, 59347, 40943,
c     10724, 18954). Runs OREM's standard orem_run pipeline once per
c     PRE-SELECTED zone ({nobj} total, from
c     scratch_rpe/validation_rpe_zones.py + prepare_validation_zones.py
c     -- 3-29 TLEs/zone, first zone >=6mo before each object's own
c     U-turn epoch, spread through to the last tracked TLE, matching
c     the "E:\\Research\\1. R&D\\Re-entry\\COSPAR ASR" reference's own
c     per-zone TLE-split style), rather than letting zone_select.F
c     pick zones itself.
c
c     GENERATED FILE -- do not hand-edit. Regenerate via
c     `python scratch_rpe/prepare_validation_zones.py` after any change
c     to scratch_rpe/validation_rpe_zones.py's OBJECTS/zone logic.
c
c     Each "object" below is really one (norad, zone) pair -- its own
c     TLE input file already contains ONLY that zone's TLE records
c     (input/example_<norad>_valzoneN.tle.txt). min_zone_pts=3 (down
c     from production's 8, since 3 is this campaign's own floor) and
c     r2_thresh/slope_thresh are maximally permissive (0.0 / +1.0d6)
c     so zone_select.F's decay-linearity filter -- built for the
c     terminal monotonic-decay regime -- doesn't reject the earlier,
c     still-oscillating pre-resonance windows this campaign
c     deliberately includes; the externally-chosen window is trusted
c     as-is, not re-filtered. nzones_max=1 since each input file is
c     already exactly one zone. max_zone_days=900 covers the widest
c     zone (10724's sparse-tracking zone 1, ~850 real days for just 29
c     points). All other GA/force-model parameters are identical to
c     the main scratch_rpe/rpe_campaign.F 97-object campaign for
c     direct comparability.
      program validation_rpe_campaign

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      integer io_start, io_end, iarg
      character*80 csv_name
      character*16 arg_str

      parameter(ndim_a = 291, mxz = 1, nobj = {nobj})

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

c     norad/oyr/omo/ody DATA statements
{data_block}
c     tle_file/obj_name string assignments (executable statements)
{assign_block}
      io_start = 1
      io_end = nobj
      csv_name = 'scratch_rpe/validation_RPE/valrpe_campaign.csv'
      iarg = iargc()
      if (iarg .ge. 2) then
         call getarg(1, arg_str)
         read(arg_str, *) io_start
         call getarg(2, arg_str)
         read(arg_str, *) io_end
         write(csv_name, '(A,I2.2,A,I2.2,A)')
     &      'scratch_rpe/validation_RPE/campaign_part_', io_start,
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
     &      t_obs, mxz, 3, 900.d0, 0.d0, 1.d6,
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
     &         '  RPE = ', rpe_ens, '%'
         else
            write(*,'(A)') '  no re-entry (weak-signal fit, or '//
     &         'zone_select still rejected this window)'
         end if

 90      continue
         write(*,'(A,I3,A,I2,A,I2,A)') ' [PROGRESS] ',
     &      ((io-io_start+1)*100)/(io_end-io_start+1), '% (',
     &      io, ' of ', io_end, ' zones done, this chunk)'
      end do

      close(33)
      write(*,'(/,A,A)') ' campaign complete -> ', trim(csv_name)

      stop
      end

c     ============================================================
c     ATM reader (same as rpe_campaign.F / rpe_campaign_resonance16.F)
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
"""

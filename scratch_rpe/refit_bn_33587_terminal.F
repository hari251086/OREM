c     refit_bn_33587_terminal.F -- issue #27's declined-twice next
c     step: refit (e, BN) directly against the in-record TLE window,
c     instead of propagating forward from an assumed BN=55.
c     -------------------------------------------------------
c     The #27 decisive experiment (scratch_gmat/gmat_hp_33587_
c     issue27_{full,grav}.script + prop_33587_hp.F) propagated the
c     2022-07-25 TLE state 131 days at a FIXED BN=55 and found flat
c     hp(t) in both OREM and near-exact-gravity GMAT, while the real
c     TLE record shows hp collapsing 616->341 km over the same
c     window -- ruling out third-body truncation (GMAT agrees with
c     OREM) and pointing at drag (raw TLE mean elements: sma
c     11128->9211 km with simultaneous circularization e 0.37->0.27,
c     the textbook drag signature). This driver tests the resulting
c     hypothesis directly: treat the whole 131-day span as ONE zone
c     (bypassing zone_select's 10-day cap -- that cap is why the
c     standard pipeline only ever finds ~9-10 day sub-windows here,
c     see zone_windows_33587.F) and let the standard RSM/GA machinery
c     fit (e, BN) against the REAL observed apogee trajectory across
c     it, with a wide BN search range (the old 80-160 assumption is
c     exactly what's in question).
c
c     If a low, physically-plausible BN fits this window well (RMS
c     low, not pinned at a search bound), that confirms BN
c     identifiability -- the object's effective BN changed in this
c     terminal regime (attitude/tumbling change is the most likely
c     physical cause for a fragment) -- rather than missing gravity
c     physics, consistent with the GMAT finding. Also re-propagates
c     from the same start state at the FITTED BN (weather on) and
c     dumps hp(t) for a direct before/after comparison against
c     hp_model_33587.txt (the old BN=55 result).
c
c     Compile:
c       ifx /heap-arrays /F:16777216 scratch_rpe\refit_bn_33587_terminal.F
c           swx.F rsm.F ga.F tle_evolution.F
c           ksrop\propagate_ks.F ksrop\Subrouts.F
c           ksrop\TLEread.F ksrop\Legendre.F
c           /exe:scratch_rpe\refit_bn_33587_terminal.exe
c     Output: scratch_rpe\refit_bn_33587_terminal.log (console),
c             scratch_rpe\hp_model_33587_refit.txt (re-propagation)
c     -------------------------------------------------------
      program refit_bn_33587_terminal
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(maxpts = 400, maxtle = 400, ndim_a = 291)
      parameter(max_surf = 5000)

      double precision epochs(maxpts), a_tle(maxpts), e_tle(maxpts)
      double precision inc_tle(maxpts), raan_tle(maxpts)
      double precision aop_tle(maxpts)
      double precision ha_tle(maxpts), hp_tle(maxpts), ls_tle(maxpts)
      integer npts, ierr_tle

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

      double precision epz(maxpts), haz(maxpts), ez(maxpts)
      double precision incz(maxpts), raanz(maxpts), aopz(maxpts)
      double precision maz(maxpts), smaz(maxpts)

      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision surfaces(max_surf, 3, 3)
      double precision e_grid(3), bn_grid(3)
      double precision tobs(500), apobs(500)
      integer nsurf_pts, nobs, ierr_rsm

      integer n_force(3)
      double precision e_opt, bn_opt, rms_opt

      integer i, i0, i1, nzone, ierr_sw, ierr_a2

c     Re-propagation arrays (same pattern as prop_33587_hp.F)
      parameter(maxp = 1300)
      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      double precision tj(maxp), tx(3,maxp), td(3,maxp)
      integer ec, idump
      double precision r, v2, a, hx, hy, hz, p, e, hp

      call init_constants()
      call read_atm_r(ALT_a, DEN_a, SCH_a, ndim_a, natm)
      call sw_load('input/SW-All.csv', ierr_sw)
      call atm2d_load('input/ATM2D.DAT', ierr_a2)
      if (ierr_sw .eq. 0 .and. ierr_a2 .eq. 0) then
         write(*,'(A)') '[SW] epoch-resolved density: ENABLED'
      else
         write(*,'(A)') '[SW] epoch-resolved density: DISABLED'//
     &      ' -- aborting, this run requires weather mode'
         stop 1
      end if

      call tle_evolve('input/example_33587.tle.txt', 33587, maxpts,
     &   epochs, a_tle, e_tle, inc_tle, raan_tle, aop_tle,
     &   ha_tle, hp_tle, ls_tle, npts, ierr_tle)
      if (ierr_tle .ne. 0) then
         write(*,'(A,I3)') 'tle_evolve error ierr=', ierr_tle
         stop 1
      end if
      write(*,'(A,I5)') 'tle_evolve npts=', npts

      open(unit=77, file='input/example_33587.tle.txt',
     &     status='old', err=900)
      call read_tle(77, maxtle, names_raw, catid_raw,
     &   norad_raw, epoch_yr_raw, epoch_day_raw,
     &   ndot_raw, nddot_raw, bstar_raw,
     &   ainc_raw, raan_raw, ecc_raw,
     &   aop_raw, ama_raw, amm_raw,
     &   revnum_raw, ntle_raw)
      close(77)

c     Window = the SAME 131-day span as the #27 GMAT decisive
c     experiment: 2022-07-25 (JD 2459785.93615) through the record's
c     last TLE (2022-12-03), i.e. i0 = first point at/after that JD,
c     i1 = npts (the record's end).
      i0 = 0
      do i = 1, npts
         if (i0 .eq. 0 .and. epochs(i) .ge. 2459785.9d0) i0 = i
      end do
      i1 = npts
      if (i0 .eq. 0) then
         write(*,'(A)') 'window start not found'
         stop 1
      end if
      nzone = i1 - i0 + 1
      write(*,'(A,I4,A,I4,A,I4,A,F10.3,A,F10.3)')
     &   'window i0=', i0, ' i1=', i1, ' nzone=', nzone,
     &   ' epoch0=', epochs(i0), ' epoch1=', epochs(i1)

      do i = 1, nzone
         epz(i)   = epochs(i0 + i - 1)
         haz(i)   = ha_tle(i0 + i - 1)
         ez(i)    = e_tle(i0 + i - 1)
         incz(i)  = inc_tle(i0 + i - 1)
         raanz(i) = raan_tle(i0 + i - 1)
         aopz(i)  = aop_tle(i0 + i - 1)
         smaz(i)  = a_tle(i0 + i - 1)
         call tle_find_ma(33587, epz(i),
     &      norad_raw, epoch_yr_raw, epoch_day_raw,
     &      ama_raw, ntle_raw, maxtle, maz(i))
      end do

      write(*,'(A,F8.4,A,F8.4)') 'e range in window: ',
     &   ez(1), ' -> ', ez(nzone)
      write(*,'(A,F10.2,A,F10.2)') 'apogee alt range: ',
     &   haz(1), ' -> ', haz(nzone)

      e_mid = 0.5d0*(ez(1) + ez(nzone))
      de = 0.5d0*dabs(ez(1) - ez(nzone))
      if (de .lt. 0.02d0) de = 0.02d0

      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

c     Wide BN search: the old 80-160 (or 55 point estimate) is
c     exactly what's in question. [1,150] brackets generously below
c     and keeps the upper end consistent with the rest of the
c     validation set.
      write(*,'(A)') 'calling rsm_generate over the full 131-day'//
     &   ' window as one zone...'
      call rsm_generate(
     &   epz, haz, ez, incz, raanz, aopz, maz, smaz,
     &   nzone,
     &   e_mid - de, e_mid + de,
     &   1.d0, 150.d0,
     &   200, 360, 1.d-10,
     &   n_force, 20, 2, 3,
     &   7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   surfaces, nsurf_pts,
     &   tobs, apobs, nobs,
     &   e_grid, bn_grid,
     &   ierr_rsm)

      write(*,'(A,I3,A,I5,A,I4)') 'ierr_rsm=', ierr_rsm,
     &   ' nsurf_pts=', nsurf_pts, ' nobs=', nobs
      if (ierr_rsm .ne. 0 .or. nsurf_pts .le. 0) then
         write(*,'(A)') 'RSM failed -- cannot fit this window'
         stop 1
      end if

      call ga_optimize(
     &   surfaces, max_surf, e_grid, bn_grid,
     &   3, 3, nsurf_pts,
     &   tobs, apobs, nobs,
     &   e_grid(1), e_grid(3),
     &   bn_grid(1), bn_grid(3),
     &   20, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &   e_opt, bn_opt, rms_opt)

      write(*,'(A)') '==================================='
      write(*,'(A,F10.5)') 'FITTED e_opt  = ', e_opt
      write(*,'(A,F10.3)') 'FITTED bn_opt = ', bn_opt
      write(*,'(A,F10.4)') 'RMS           = ', rms_opt
      write(*,'(A,F8.2,A,F8.2)') 'BN search range was [',
     &   bn_grid(1), ', ', bn_grid(3), ']'
      write(*,'(A)') '==================================='

c     Re-propagate from the SAME 2022-07-25 state at the FITTED BN,
c     weather on, for a direct before/after comparison against
c     hp_model_33587.txt (the old fixed-BN=55 result).
      a_kep(1) = 11128.043d0
      a_kep(2) = 0.3714527d0
      a_kep(3) = 65.0291d0
      a_kep(4) = 250.9049d0
      a_kep(5) = 105.5241d0
      a_kep(6) = 297.8486d0
      a_kep(7) = 0.d0

      call oe2car(a_kep, x0, xd0, 1.d-12)
      call jd2cal(2459785.93615d0, cal0)

      idump = 1
      call propagate_ks(
     &   x0, xd0, cal0,
     &   maxp-1, 360, 1.d-10,
     &   n_force, 20, 2, 3,
     &   bn_opt, 1, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   maxp, idump,
     &   tj, tx, td,
     &   ec)

      write(*,'(A,I2,A,I6)') 're-propagation exit=', ec,
     &   ' revs=', idump

      open(unit=44, file='scratch_rpe/hp_model_33587_refit.txt',
     &     status='unknown')
      do i = 1, idump
         r  = dsqrt(tx(1,i)**2 + tx(2,i)**2 + tx(3,i)**2)
         v2 = td(1,i)**2 + td(2,i)**2 + td(3,i)**2
         a  = 1.d0/(2.d0/r - v2/amue)
         hx = tx(2,i)*td(3,i) - tx(3,i)*td(2,i)
         hy = tx(3,i)*td(1,i) - tx(1,i)*td(3,i)
         hz = tx(1,i)*td(2,i) - tx(2,i)*td(1,i)
         p  = (hx*hx + hy*hy + hz*hz)/amue
         e  = dsqrt(dabs(1.d0 - p/a))
         hp = a*(1.d0 - e) - R_Earth
         write(44,'(F14.5,F10.2)') tj(i), hp
      end do
      close(44)
      write(*,'(A)') 'wrote scratch_rpe/hp_model_33587_refit.txt'

      stop
 900  write(*,'(A)') 'TLE file open error'
      stop 1
      end

c     ============================================================
      subroutine read_atm_r(alt, den, sch, maxn, natm)
      implicit double precision (a-h, o-z)
      integer maxn, natm
      double precision alt(maxn), den(maxn), sch(maxn)
      natm = 291
      alt(1) = 60.d0
      do i = 2, 141
         alt(i) = alt(i-1) + 1.d0
      end do
      do i = 142, natm
         alt(i) = alt(i-1) + 2.d0
      end do
      open(unit=88, file='input/ATM.DAT', status='old', err=900)
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

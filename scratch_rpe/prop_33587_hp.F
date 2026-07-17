c     prop_33587_hp.F — in-record lunisolar hp-descent diagnostic
c     (issue #26 follow-up / #9 attribution, 2026-07-17)
c     -------------------------------------------------------
c     The real 33587 TLE record shows perigee altitude falling
c     616 -> 341 km between 2022-07-25 and 2022-12-03 -- a descent
c     driven by lunisolar perigee cycling at the critical
c     inclination (i=65.0), not drag. This driver propagates the
c     2022-07-25 TLE state across exactly that window (weather
c     mode on, BN=55) and dumps per-revolution osculating perigee
c     altitude, so the modeled descent can be compared against the
c     observed one with data that is entirely inside the TLE
c     record. If the model cannot reproduce THIS descent, the
c     33587 lifetime error is attributed to the lunisolar/critical-
c     inclination dynamics, not to the atmosphere.
c
c     Output: scratch_rpe/hp_model_33587.txt  (JD, hp_km per rev)
c
c     Compile:
c       ifx /heap-arrays /F:16777216 scratch_rpe\prop_33587_hp.F
c           swx.F rsm.F ksrop\propagate_ks.F ksrop\Subrouts.F
c           ksrop\Legendre.F /exe:scratch_rpe\prop_33587_hp.exe
c     -------------------------------------------------------
      program prop_33587_hp
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291, maxp = 1300)
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm
      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      double precision tj(maxp), tx(3,maxp), td(3,maxp)
      integer n_force(3), ec, idump
      double precision r, v2, a, hx, hy, hz, p, e, hp
      integer i, ierr_sw, ierr_a2

      call init_constants()
      call read_atm_h(ALT_a, DEN_a, SCH_a, ndim_a, natm)
      call sw_load('input/SW-All.csv', ierr_sw)
      call atm2d_load('input/ATM2D.DAT', ierr_a2)
      if (ierr_sw .eq. 0 .and. ierr_a2 .eq. 0) then
         write(*,'(A)') ' [SW] epoch-resolved density ENABLED'
      end if

c     TLE 22206.43614790: a=11128.043 e=0.3714527 i=65.0291
c     raan=250.9049 aop=105.5241 ma=297.8486  (hp = 616.4 km)
      a_kep(1) = 11128.043d0
      a_kep(2) = 0.3714527d0
      a_kep(3) = 65.0291d0
      a_kep(4) = 250.9049d0
      a_kep(5) = 105.5241d0
      a_kep(6) = 297.8486d0
      a_kep(7) = 0.d0

      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

      call oe2car(a_kep, x0, xd0, 1.d-12)
      call jd2cal(2459785.93615d0, cal0)

      idump = 1
      call propagate_ks(
     &   x0, xd0, cal0,
     &   maxp-1, 360, 1.d-10,
     &   n_force, 20, 2, 3,
     &   55.d0, 1, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   maxp, idump,
     &   tj, tx, td,
     &   ec)

      write(*,'(A,I2,A,I6)') ' exit=', ec, ' revs=', idump

      open(unit=44, file='scratch_rpe/hp_model_33587.txt',
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
      write(*,'(A)') ' wrote scratch_rpe/hp_model_33587.txt'

      stop
      end

c     ============================================================
      subroutine read_atm_h(alt, den, sch, maxn, natm)
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

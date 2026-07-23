c     test_sw.F — Space weather + 2-D atmosphere tests (issue #26)
c     -------------------------------------------------------
c     Compile:
c       ifx /heap-arrays /F:16777216 test_sw.F swx.F orem.F report.F
c           rsm.F ga.F tle_evolution.F zone_select.F
c           ksrop/propagate_ks.F ksrop/Subrouts.F
c           ksrop/TLEread.F ksrop/Legendre.F
c           /exe:test_sw.exe
c     Needs input/SW-All.csv and input/ATM2D.DAT (committed).
c     -------------------------------------------------------
      program test_sw

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      integer npass, nfail, ierr
      double precision cal(6), jd, tinf, tinf2, rho, hs
      double precision rho72, hs72, rho_hi
      integer ist

      call init_constants()
      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' Space weather / ATM2D test suite (issue #26)'
      write(*,*) '=============================================='

c     ---------------------------------------------------------
c     W1-W2: defaults — nothing loaded, legacy contract holds
c     ---------------------------------------------------------
      call sw_tinf(2458000.d0, tinf)
      call asrt_d('W1 no data: sw_tinf=0', tinf, 0.d0, 1.d-12,
     &            npass, nfail)
      call atm2d_interp(170.d3, 700.d0, rho, hs, ist)
      call asrt_i('W2 no table: istat=2', ist, 2, npass, nfail)

c     ---------------------------------------------------------
c     W3-W6: sw_load + known-row lookups
c     ---------------------------------------------------------
      call sw_load('input/SW-All.csv', ierr)
      call asrt_i('W3 sw_load ierr=0', ierr, 0, npass, nfail)

c     2024-05-11 (G5 storm): F10.7_ADJ=218.0, ADJ_C81=180.5,
c     KP_SUM=670 -> mean Kp=8.375
c     Tc   = 379 + 3.24*180.5 + 1.3*(218.0-180.5) = 1012.57
c     Tinf = Tc + 14*8.375 + 0.02*exp(8.375)
c          = 1012.57 + 117.25 + 86.74 = 1216.56
      cal(1) = 2024.d0
      cal(2) = 5.d0
      cal(3) = 11.d0
      cal(4) = 12.d0
      cal(5) = 0.d0
      cal(6) = 0.d0
      call cal2jd(cal, jd)
      call sw_tinf(jd, tinf)
      call asrt_d('W4 2024-05-11 storm Tinf', tinf, 1216.56d0,
     &            2.d0, npass, nfail)

c     quiet solar minimum 2019-12-15: F10.7 ~ 69-72, Kp ~ 0-1:
c     Tinf must land in the 590-650 K band
      cal(1) = 2019.d0
      cal(2) = 12.d0
      cal(3) = 15.d0
      call cal2jd(cal, jd)
      call sw_tinf(jd, tinf)
      if (tinf .gt. 590.d0 .and. tinf .lt. 650.d0) then
         call asrt_i('W5 2019 solar-min Tinf in [590,650]',
     &               1, 1, npass, nfail)
      else
         write(*,'(A,F8.1)') '   got Tinf=', tinf
         call asrt_i('W5 2019 solar-min Tinf in [590,650]',
     &               0, 1, npass, nfail)
      end if

c     predicted era (monthly rows): 2030-06-15 must return
c     something physical, not zero
      cal(1) = 2030.d0
      cal(2) = 6.d0
      cal(3) = 15.d0
      call cal2jd(cal, jd)
      call sw_tinf(jd, tinf)
      if (tinf .gt. 500.d0 .and. tinf .lt. 1500.d0) then
         call asrt_i('W6 predicted era Tinf physical',
     &               1, 1, npass, nfail)
      else
         write(*,'(A,F8.1)') '   got Tinf=', tinf
         call asrt_i('W6 predicted era Tinf physical',
     &               0, 1, npass, nfail)
      end if

c     ---------------------------------------------------------
c     W7-W11: atm2d_load + interpolation consistency
c     ---------------------------------------------------------
      call atm2d_load('input/ATM2D.DAT', ierr)
      call asrt_i('W7 atm2d_load ierr=0', ierr, 0, npass, nfail)

c     at Tinf=625 (the static table's generation temperature is
c     626.33), rho(170 km) must match ATM.DAT's 4.866e-10 kg/m3
c     within the 25-K grid-cell interpolation error (~2%)
      call atm2d_interp(170.d3, 626.33d0, rho, hs, ist)
      call asrt_i('W8 interp istat=0', ist, 0, npass, nfail)
      rho72 = 4.866d-10 * 1.d10
      if (dabs(rho - rho72)/rho72 .lt. 0.02d0) then
         call asrt_i('W9 rho(170km,626K) ~ static table',
     &               1, 1, npass, nfail)
      else
         write(*,'(A,E12.4,A,E12.4)') '   got ', rho,
     &      ' want ~', rho72
         call asrt_i('W9 rho(170km,626K) ~ static table',
     &               0, 1, npass, nfail)
      end if

c     monotonicity in Tinf at 300 km (scaled units)
      call atm2d_interp(300.d3, 600.d0, rho, hs, ist)
      call atm2d_interp(300.d3, 1200.d0, rho_hi, hs72, ist)
      if (rho_hi .gt. 3.d0*rho) then
         call asrt_i('W10 rho(300km) grows strongly with Tinf',
     &               1, 1, npass, nfail)
      else
         call asrt_i('W10 rho(300km) grows strongly with Tinf',
     &               0, 1, npass, nfail)
      end if

c     altitude out of range -> istat=1 (drag-inactive contract)
      call atm2d_interp(550.d3, 700.d0, rho, hs, ist)
      call asrt_i('W11 h>500km istat=1', ist, 1, npass, nfail)

c     ---------------------------------------------------------
c     W12: weather-mode propagation smoke test — same state, one
c     epoch at solar minimum (2019) vs one at the 2024 maximum:
c     7-day apogee decay must be substantially larger at the max.
c     ---------------------------------------------------------
      call smoke_prop(npass, nfail)

      write(*,*) ''
      write(*,*) '=============================================='
      write(*,100) npass + nfail, npass, nfail
 100  format(' Total: ',i3,'   Passed: ',i3,'   Failed: ',i3)
      write(*,*) '=============================================='
      if (nfail .gt. 0) then
         write(*,*) 'SOME TESTS FAILED'
         stop 1
      else
         write(*,*) 'ALL TESTS PASSED'
      end if
      stop
      end

c     ============================================================
      subroutine smoke_prop(npass, nfail)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth
      integer npass, nfail

      parameter(ndim_a = 291, maxp = 200)
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm
      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      double precision tj(maxp), tx(3,maxp), td(3,maxp)
      integer n_force(3), ec, idump
      double precision apo0, apo1, drop_min, drop_max

      call read_atm_t(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

c     N13-like orbit, hp ~ 165 km
      a_kep(1) = 9636.790207d0
      a_kep(2) = 0.321216d0
      a_kep(3) = 19.1713d0
      a_kep(4) = 68.985106d0
      a_kep(5) = 273.564599d0
      a_kep(6) = 0.d0
      a_kep(7) = 0.d0

      drop_min = 0.d0
      drop_max = 0.d0
      do k = 1, 2
         call oe2car(a_kep, x0, xd0, 1.d-12)
         if (k .eq. 1) then
            cal0(1) = 2019.d0
         else
            cal0(1) = 2024.d0
         end if
         cal0(2) = 5.d0
         cal0(3) = 10.d0
         cal0(4) = 0.d0
         cal0(5) = 0.d0
         cal0(6) = 0.d0
         idump = 1
         call propagate_ks(
     &      x0, xd0, cal0,
     &      64, 360, 1.d-12,
     &      n_force, 4, 0, 0,
     &      80.d0, 1, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &      0.d0, 0.d0, 0, 0,
     &      0.d0, 0.d0, 0.d0,
     &      ALT_a, DEN_a, SCH_a, natm,
     &      maxp, idump,
     &      tj, tx, td,
     &      ec)
         call apoalt(tx(1,1), td(1,1), apo0)
         call apoalt(tx(1,idump), td(1,idump), apo1)
         if (k .eq. 1) drop_min = apo0 - apo1
         if (k .eq. 2) drop_max = apo0 - apo1
      end do

      write(*,'(A,F8.2,A,F8.2,A)') '   7-day drop: 2019 min=',
     &   drop_min, ' km, 2024 max=', drop_max, ' km'
      if (drop_max .gt. 1.5d0*drop_min .and. drop_min .gt. 0.d0)
     &then
         call asrt_i('W12 solar-max decay >> solar-min decay',
     &               1, 1, npass, nfail)
      else
         call asrt_i('W12 solar-max decay >> solar-min decay',
     &               0, 1, npass, nfail)
      end if
      return
      end

c     ============================================================
      subroutine apoalt(x, xd, ha)
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth
      double precision x(3), xd(3), ha
      r = dsqrt(x(1)**2 + x(2)**2 + x(3)**2)
      v2 = xd(1)**2 + xd(2)**2 + xd(3)**2
      a = 1.d0/(2.d0/r - v2/amue)
      hx = x(2)*xd(3) - x(3)*xd(2)
      hy = x(3)*xd(1) - x(1)*xd(3)
      hz = x(1)*xd(2) - x(2)*xd(1)
      p = (hx*hx + hy*hy + hz*hz)/amue
      e = dsqrt(dabs(1.d0 - p/a))
      ha = a*(1.d0 + e) - R_Earth
      return
      end

c     ============================================================
      subroutine read_atm_t(alt, den, sch, maxn, natm)
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

c     ============================================================
      subroutine asrt_i(label, got, expect, npass, nfail)
      character*(*) label
      integer got, expect, npass, nfail
      if (got .eq. expect) then
         write(*,'(A,A)') '  [PASS]  ', label
         npass = npass + 1
      else
         write(*,'(A,A)') '  [FAIL]  ', label
         write(*,'(A,I10,A,I10)') '           got=',
     &        got, '  expect=', expect
         nfail = nfail + 1
      end if
      return
      end

c     ============================================================
      subroutine asrt_d(label, got, expect, tol, npass, nfail)
      implicit double precision (a-h, o-z)
      character*(*) label
      double precision got, expect, tol
      integer npass, nfail
      if (dabs(got - expect) .le. tol) then
         write(*,'(A,A)') '  [PASS]  ', label
         npass = npass + 1
      else
         write(*,'(A,A)') '  [FAIL]  ', label
         write(*,'(A,E20.12,A,E20.12)') '           got=',
     &        got, '  expect=', expect
         nfail = nfail + 1
      end if
      return
      end

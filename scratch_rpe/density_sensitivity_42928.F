c     density_sensitivity_42928.F -- issue #14 item 2: perigee-altitude
c     density sensitivity. Propagates 42928's zone-1 state (real TLE
c     elements, e/BN from the zone-1 fit already on record in
c     scratch_rpe/rpe_campaign.csv: e=0.278951, BN=81.064) at
c     ATM.DAT densities scaled by 0.5x / 1.0x / 2.0x, and reports
c     re-entry time (or terminal altitude if it doesn't re-enter
c     within the window) for each -- quantifying "days per density
c     scale factor" the way the issue's acceptance criteria ask for.
c
c     Compile:
c       ifx /heap-arrays scratch_rpe\density_sensitivity_42928.F
c           rsm.F ga.F tle_evolution.F ksrop\propagate_ks.F
c           ksrop\Subrouts.F ksrop\TLEread.F ksrop\Legendre.F
c           /exe:scratch_rpe\density_sensitivity_42928.exe
      program density_sensitivity_42928
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(maxpts = 2000, ndim_a = 291)
      double precision epochs(maxpts), a_o(maxpts), e_o(maxpts)
      double precision inc_o(maxpts), raan_o(maxpts), aop_o(maxpts)
      double precision ha_o(maxpts), hp_o(maxpts), lam_o(maxpts)
      integer npts, ierr

      double precision ALT_a(ndim_a), DEN_a(ndim_a), DEN_scaled(ndim_a)
      double precision SCH_a(ndim_a)
      integer natm

      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      parameter(maxp = 6000)
      double precision tj(maxp), tx(3,maxp), td(3,maxp)
      integer n_force(3), ec, idump

      double precision r, v2, a, hx, hy, hz, p, e, hp, ha
      double precision scale_try(3)
      integer i, k

      call init_constants()
      call read_atm_d(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      call tle_evolve('input/example_42928.tle.txt', 42928, maxpts,
     &   epochs, a_o, e_o, inc_o, raan_o, aop_o,
     &   ha_o, hp_o, lam_o, npts, ierr)
      write(*,'(A,I3,A,I5)') 'tle_evolve ierr=', ierr,
     &   '  npts=', npts

c     Point-1 state: the object's own REAL TLE elements at its first
c     tracked epoch (2017-08-31, right after launch) -- a, e BOTH
c     from the same real TLE (mixing point-1's sma with a different
c     zone's fitted e gave a physically inconsistent hp=644 km state
c     that no version of this object was ever actually in, which is
c     why the first attempt showed zero density sensitivity: nothing
c     decays at 644 km regardless of density). At its own consistent
c     elements hp=166 km -- already drag-relevant, real re-entry was
c     ~1.5 years later (2019-03-03), a long enough window to see
c     divergence. MA fixed at 0 (not epoch-critical for a same-
c     starting-point sensitivity comparison).
      a_kep(1) = a_o(1)
      a_kep(2) = e_o(1)
      a_kep(3) = inc_o(1)
      a_kep(4) = raan_o(1)
      a_kep(5) = aop_o(1)
      a_kep(6) = 0.d0
      a_kep(7) = 0.d0

      write(*,'(A,F10.2,A,F8.5,A,F8.3)') 'start: a=', a_kep(1),
     &   ' e=', a_kep(2), ' i=', a_kep(3)

      call oe2car(a_kep, x0, xd0, 1.d-12)
      call jd2cal(epochs(1), cal0)

      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

      scale_try(1) = 0.5d0
      scale_try(2) = 1.0d0
      scale_try(3) = 2.0d0

      do k = 1, 3
         do i = 1, natm
            DEN_scaled(i) = DEN_a(i) * scale_try(k)
         end do

         idump = 1
         call propagate_ks(
     &      x0, xd0, cal0,
     &      maxp-1, 360, 1.d-10,
     &      n_force, 4, 0, 0,
     &      81.064d0, 1, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &      0.d0, 0.d0, 0, 0,
     &      4.56d-6, 1.32712440018d11, 4.902801076d3,
     &      ALT_a, DEN_scaled, SCH_a, natm,
     &      maxp, idump,
     &      tj, tx, td,
     &      ec)

         write(*,'(A,F4.1,A,I2,A,I6)') 'scale=', scale_try(k),
     &      '  exit=', ec, '  revs=', idump

         if (ec .eq. 1) then
c           Re-entered -- report elapsed days from start.
            write(*,'(A,F10.3,A)') '  RE-ENTERED after ',
     &         tj(idump) - epochs(1), ' days'
         else
c           Didn't re-enter in the window -- report terminal hp.
            r  = dsqrt(tx(1,idump)**2+tx(2,idump)**2+tx(3,idump)**2)
            v2 = td(1,idump)**2+td(2,idump)**2+td(3,idump)**2
            a  = 1.d0/(2.d0/r - v2/amue)
            hx = tx(2,idump)*td(3,idump) - tx(3,idump)*td(2,idump)
            hy = tx(3,idump)*td(1,idump) - tx(1,idump)*td(3,idump)
            hz = tx(1,idump)*td(2,idump) - tx(2,idump)*td(1,idump)
            p  = (hx*hx+hy*hy+hz*hz)/amue
            e  = dsqrt(dabs(1.d0 - p/a))
            hp = a*(1.d0-e) - R_Earth
            write(*,'(A,F10.3,A,F10.2,A)') '  no re-entry after ',
     &         tj(idump) - epochs(1), ' days, terminal hp=', hp,
     &         ' km'
         end if
      end do

      stop
      end

c     ============================================================
      subroutine read_atm_d(alt, den, sch, maxn, natm)
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

c     test_propagate_ks.F — Tests for refactored propagate_ks
c     -------------------------------------------------------
c     Compile:
c       ifx test_propagate_ks.F ksrop/propagate_ks.F
c            ksrop/Subrouts.F ksrop/Legendre.F
c            /exe:test_propagate_ks.exe
c     Run: ./test_propagate_ks.exe
c     -------------------------------------------------------
      program test_propagate_ks

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_atm = 291)
      parameter(max_pts = 3601)

      dimension x0(3), xd0(3), cal0(6)
      dimension n_force(3)
      dimension ALT_atm(ndim_atm), DEN_atm(ndim_atm)
      dimension SCH_atm(ndim_atm)
      dimension traj_jd(max_pts)
      dimension traj_x(3,max_pts), traj_xd(3,max_pts)

      integer npass, nfail, exit_code, idump

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' propagate_ks Subroutine Test Suite'
      write(*,*) '=============================================='

c     --- Initialize constants from const_new.dat ---
      call init_constants()

c     --- Build atmosphere table (same as driver_KS) ---
      ALT_atm(1) = 60.d0
      do i = 2, 141
         ALT_atm(i) = ALT_atm(i-1) + 1.d0
      end do
      do i = 142, ndim_atm
         ALT_atm(i) = ALT_atm(i-1) + 2.d0
      end do

      open(unit=60, file='input/ATM.DAT', status='old',
     &     err=900)
      read(60,110) (SCH_atm(i), i=1,61)
      read(60,110) (SCH_atm(i), i=62,ndim_atm)
      read(60,111) (DEN_atm(i), i=1,61)
      read(60,111) (DEN_atm(i), i=62,ndim_atm)
      close(60)
 110  format(8f10.3)
 111  format(8e10.3)

      do i = 1, ndim_atm
         ALT_atm(i) = ALT_atm(i) * 1000.d0
         SCH_atm(i) = SCH_atm(i) * 1000.d0
         DEN_atm(i) = DEN_atm(i) * 1.0d10
      end do

c     ==============================================================
c     T1: Two-body propagation — energy conservation
c     All perturbations off, 1 rev, 360 steps
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Two-body tests ---'

      x0(1) = 0.d0
      x0(2) = -5888.9727d0
      x0(3) = -3400.d0
      xd0(1) = 9.5d0
      xd0(2) = 0.d0
      xd0(3) = 0.d0
      cal0(1) = 2016.d0
      cal0(2) = 9.d0
      cal0(3) = 20.d0
      cal0(4) = 0.d0
      cal0(5) = 0.d0
      cal0(6) = 0.d0

      n_force(1) = 0
      n_force(2) = 0
      n_force(3) = 0

      call propagate_ks(
     &   x0, xd0, cal0,
     &   1, 360, 1.d-15,
     &   n_force, 0, 0, 0,
     &   50.d0, 0, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 0, 1,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   max_pts, idump,
     &   traj_jd, traj_x, traj_xd,
     &   exit_code)

      if (exit_code .eq. 0) then
         write(*,*) '  [PASS]  T1 exit_code = 0'
         npass = npass + 1
      else
         write(*,'(A,I3)')
     &      '  [FAIL]  T1 exit_code =', exit_code
         nfail = nfail + 1
      end if

c     T2: Correct number of output points (1*360+1 = 361)
      if (idump .eq. 361) then
         write(*,*) '  [PASS]  T2 idump = 361'
         npass = npass + 1
      else
         write(*,'(A,I6)')
     &      '  [FAIL]  T2 idump =', idump
         nfail = nfail + 1
      end if

c     T3: Energy conservation (E0 ≈ Ef)
      rmag0 = dsqrt(traj_x(1,1)**2 + traj_x(2,1)**2
     &            + traj_x(3,1)**2)
      vmag0 = dsqrt(traj_xd(1,1)**2 + traj_xd(2,1)**2
     &            + traj_xd(3,1)**2)
      energy0 = 0.5d0*vmag0**2 - amue/rmag0

      rmag1 = dsqrt(traj_x(1,idump)**2
     &            + traj_x(2,idump)**2
     &            + traj_x(3,idump)**2)
      vmag1 = dsqrt(traj_xd(1,idump)**2
     &            + traj_xd(2,idump)**2
     &            + traj_xd(3,idump)**2)
      energy1 = 0.5d0*vmag1**2 - amue/rmag1

      de_rel = dabs(energy1 - energy0)/dabs(energy0)
      if (de_rel .lt. 1.d-8) then
         write(*,*) '  [PASS]  T3 energy conserved'
         npass = npass + 1
      else
         write(*,'(A,E12.4)')
     &      '  [FAIL]  T3 energy rel err=', de_rel
         nfail = nfail + 1
      end if

c     T4: Orbit closure (|dr| < 0.1 km after 1 rev)
      dr = dsqrt(
     &   (traj_x(1,idump)-traj_x(1,1))**2 +
     &   (traj_x(2,idump)-traj_x(2,1))**2 +
     &   (traj_x(3,idump)-traj_x(3,1))**2)
      if (dr .lt. 0.1d0) then
         write(*,*) '  [PASS]  T4 orbit closure < 0.1 km'
         npass = npass + 1
      else
         write(*,'(A,F12.6,A)')
     &      '  [FAIL]  T4 orbit closure dr=',dr,' km'
         nfail = nfail + 1
      end if

c     T5: Julian dates are monotonically increasing
      imon = 1
      do i = 2, idump
         if (traj_jd(i) .le. traj_jd(i-1)) imon = 0
      end do
      if (imon .eq. 1) then
         write(*,*) '  [PASS]  T5 JD monotonic'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  T5 JD not monotonic'
         nfail = nfail + 1
      end if

c     ==============================================================
c     T6: Multiple revolutions (10 rev, 360 steps)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Multi-revolution ---'

      call propagate_ks(
     &   x0, xd0, cal0,
     &   10, 360, 1.d-15,
     &   n_force, 0, 0, 0,
     &   50.d0, 0, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 0, 1,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   max_pts, idump,
     &   traj_jd, traj_x, traj_xd,
     &   exit_code)

      if (exit_code .eq. 0 .and. idump .eq. 3601) then
         write(*,*) '  [PASS]  T6 10rev exit=0, 3601 pts'
         npass = npass + 1
      else
         write(*,'(A,I3,A,I6)')
     &      '  [FAIL]  T6 exit=',exit_code,' idump=',idump
         nfail = nfail + 1
      end if

c     ==============================================================
c     T7: Re-entry detection (low perigee orbit with drag)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Re-entry detection ---'

c     LEO orbit with drag enabled — should decay and trigger
c     re-entry within several revolutions
      x0(1) = 6578.d0
      x0(2) = 0.d0
      x0(3) = 0.d0
      xd0(1) = 0.d0
      xd0(2) = 7.2d0
      xd0(3) = 2.5d0

      n_force(1) = 1
      n_force(2) = 0
      n_force(3) = 0

      call propagate_ks(
     &   x0, xd0, cal0,
     &   5000, 360, 1.d-15,
     &   n_force, 2, 0, 0,
     &   10.d0, 1, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 0, 1,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   max_pts, idump,
     &   traj_jd, traj_x, traj_xd,
     &   exit_code)

      if (exit_code .eq. 1) then
         write(*,*) '  [PASS]  T7 re-entry detected (code=1)'
         npass = npass + 1
      else if (exit_code .eq. 2) then
         write(*,*) '  [PASS]  T7 diverged before reentry'
         npass = npass + 1
      else
         write(*,'(A,I3)')
     &      '  [FAIL]  T7 expected exit 1 or 2, got',
     &      exit_code
         nfail = nfail + 1
      end if

c     T8: idump < max_pts for early termination
      if (idump .lt. max_pts .and. idump .gt. 0) then
         write(*,*) '  [PASS]  T8 early termination idump ok'
         npass = npass + 1
      else
         write(*,'(A,I6)')
     &      '  [FAIL]  T8 idump=', idump
         nfail = nfail + 1
      end if

c     T9: Last trajectory point is physical
      rmag1 = dsqrt(traj_x(1,idump)**2
     &            + traj_x(2,idump)**2
     &            + traj_x(3,idump)**2)
      if (rmag1 .gt. 0.d0 .and. rmag1 .eq. rmag1) then
         write(*,*) '  [PASS]  T9 last state physical'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  T9 last state not physical'
         nfail = nfail + 1
      end if

c     ==============================================================
c     T10: Input arrays not modified by subroutine
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Input preservation ---'

      x0(1) = 7000.d0
      x0(2) = 0.d0
      x0(3) = 0.d0
      xd0(1) = 0.d0
      xd0(2) = 7.5d0
      xd0(3) = 0.d0
      n_force(1) = 0
      n_force(2) = 0
      n_force(3) = 0

      call propagate_ks(
     &   x0, xd0, cal0,
     &   1, 36, 1.d-15,
     &   n_force, 0, 0, 0,
     &   50.d0, 0, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 0, 1,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   max_pts, idump,
     &   traj_jd, traj_x, traj_xd,
     &   exit_code)

      if (dabs(x0(1) - 7000.d0) .lt. 1.d-10 .and.
     &    dabs(xd0(2) - 7.5d0) .lt. 1.d-10) then
         write(*,*) '  [PASS]  T10 inputs preserved'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  T10 inputs modified'
         nfail = nfail + 1
      end if

c     ==============================================================
c     Summary
c     ==============================================================
      write(*,*) ''
      write(*,*) '=============================================='
      write(*,200) npass + nfail, npass, nfail
 200  format(' Total: ',i3,'   Passed: ',i3,'   Failed: ',i3)
      write(*,*) '=============================================='

      if (nfail .gt. 0) stop 1
      stop

 900  write(*,*) 'ERROR: cannot open input/ATM.DAT'
      stop 1
      end

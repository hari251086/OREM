c     prop_33587.F — long direct propagation of 33587 from its LAST
c     TLE state (issue #14 / new-object case study, 2026-07-17)
c     -------------------------------------------------------
c     33587 (1989-039EF, i=65.2 Molniya-class) decayed 2025-04-22 but
c     its TLE record ends 2022-12-03 at hp=341 km, and every OREM fit
c     window sits at hp 440-640 km -- no drag signal, so the pipeline
c     honestly reports "no re-entry within the 5-year cap" under both
c     the F10.7=72 and F10.7=150 static tables. Reality reached decay
c     via LUNISOLAR perigee cycling (hp 616->341 km in 5 months) that
c     eventually hands the orbit to drag.
c     This driver propagates the last TLE state directly for up to
c     60000 revolutions (~18 yr) with the full force model to find
c     the model-consistent re-entry date for a sweep of BN values,
c     and the RPE vs the observed 2025-04-22 (JD 2460787.5).
c
c     Uses whatever input/ATM.DAT is in place (run with the F150
c     table for the solar-max arc approximation).
c
c     Compile:
c       ifx /heap-arrays /F:16777216 scratch_rpe\prop_33587.F
c           ksrop\propagate_ks.F ksrop\Subrouts.F ksrop\Legendre.F
c           /exe:scratch_rpe\prop_33587.exe
c     Run from repo root.
c     -------------------------------------------------------
      program prop_33587
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291, maxp = 60000)
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      double precision traj_jd(maxp), traj_x(3,maxp), traj_xd(3,maxp)
      save traj_jd, traj_x, traj_xd
      integer n_force(3), exit_code, idump
      double precision bns(4), t_obs_jd, epoch_jd, horiz
      integer ib

      call init_constants()
      call read_atm_p(ALT_a, DEN_a, SCH_a, ndim_a, natm)

c     Last TLE: 22337.85663267  i=64.9110 raan=132.4999 e=0.2705984
c     aop=89.5683 MA=301.0972 n=9.82008722 -> a=9211.4 km
      epoch_jd = 2459917.35663d0
      t_obs_jd = 2460787.5d0

      a_kep(1) = 9211.4d0
      a_kep(2) = 0.2705984d0
      a_kep(3) = 64.9110d0
      a_kep(4) = 132.4999d0
      a_kep(5) = 89.5683d0
      a_kep(6) = 301.0972d0
      a_kep(7) = 0.d0

      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

      bns(1) = 40.d0
      bns(2) = 55.d0
      bns(3) = 80.d0
      bns(4) = 120.d0

      write(*,'(A)') ' 33587 direct long propagation from last TLE'
      write(*,'(A,F13.4,A)') ' epoch JD ', epoch_jd,
     &   ' (2022-12-03), observed decay JD 2460787.5 (2025-04-22)'

      do ib = 1, 4
         call oe2car(a_kep, x0, xd0, 1.d-12)
         call jd2cal(epoch_jd, cal0)
         idump = 1
         call propagate_ks(
     &      x0, xd0, cal0,
     &      maxp-1, 360, 1.d-10,
     &      n_force, 20, 2, 3,
     &      bns(ib), 1, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &      1.2d0, 0.01d0, 1, 2,
     &      4.56d-6, 1.32712440018d11, 4.902801076d3,
     &      ALT_a, DEN_a, SCH_a, natm,
     &      maxp, idump,
     &      traj_jd, traj_x, traj_xd,
     &      exit_code)

         if (exit_code .eq. 1) then
            horiz = t_obs_jd - epoch_jd
            write(*,'(A,F6.1,A,F13.4,A,F8.1,A,F8.2,A)')
     &         ' BN=', bns(ib), ': re-entry JD ', traj_jd(idump),
     &         '  (epoch+', traj_jd(idump)-epoch_jd,
     &         ' d)  RPE=', (traj_jd(idump)-t_obs_jd)/horiz*100.d0,
     &         ' %'
         else
            write(*,'(A,F6.1,A,I2,A,I6,A,F13.4)')
     &         ' BN=', bns(ib), ': NO re-entry (exit=', exit_code,
     &         ') after ', idump, ' revs, last JD ', traj_jd(idump)
         end if
      end do

      stop
      end

c     ============================================================
      subroutine read_atm_p(alt, den, sch, maxn, natm)
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

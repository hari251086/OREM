c     bulge_xval_42928z8.F -- isolates the diurnal bulge term's effect
c     on 42928 Z8's own re-entry propagation (issue #32 Phase 5
c     diurnal-bulge follow-up), independent of the trust-gated
c     BN-carryover chain the full campaign RPE metric goes through.
c
c     Reuses the exact SGP4-osculating IC from the original GMAT
c     cross-check (scratch_gmat/gmat_xval_42928z8_2026.script header:
c     sma=7524.293124, e_osc=0.133643, inc=19.094467, raan=210.267993,
c     aop=73.147372, ma=301.305632, epoch JD=2458455.471721), with the
c     perigee-preserving SMA reshape (rp = sma*(1-e_osc);
c     a_new = rp/(1-e_opt)) using the WEATHER campaign's own pre/post
c     -bulge Z8 fits (from scratch_rpe/rpe_campaign_weather.csv):
c       pre-bulge:  e_opt=0.051538, bn_opt=128.040
c       post-bulge: e_opt=0.051814, bn_opt=128.532
c     Propagates each forward (same force model as the campaign) and
c     reports elapsed days to 102 km altitude (GMAT's own stop
c     threshold, per the original cross-check) -- comparable directly
c     to the original finding: OREM(no bulge, static)=86.4d,
c     GMAT(JacchiaRoberts)=48.3d.
c
c     Analysis only -- not part of the build/test suite. Link against
c     EITHER the pre-bulge or the post-bulge propagate_ks.F to get
c     each side of the comparison (same object file cannot contain
c     both).
      program bulge_xval_42928z8
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(ndim_a = 291, maxp = 20000)
      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm, ierr_sw, ierr_a2
      double precision a_kep(7), x0(3), xd0(3), cal0(6)
      double precision tj(maxp), tx(3,maxp), td(3,maxp)
      integer n_force(3), ec, idump
      double precision rp, sma_osc, e_osc, jd0
      double precision e_opt, bn_opt
      double precision r_mag, alt_km, alt_prev, jd_hit
      integer i

      call init_constants()

      call sw_load('input/SW-All.csv', ierr_sw)
      call atm2d_load('input/ATM2D.DAT', ierr_a2)
      write(*,'(A,I2,A,I2)') ' sw_load ierr=', ierr_sw,
     &   '  atm2d_load ierr=', ierr_a2

      natm = ndim_a
      ALT_a(1) = 60.d0
      do i = 2, 141
         ALT_a(i) = ALT_a(i-1) + 1.d0
      end do
      do i = 142, natm
         ALT_a(i) = ALT_a(i-1) + 2.d0
      end do
      open(unit=60, file='input/ATM.DAT', status='old')
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

      sma_osc = 7524.293124d0
      e_osc   = 0.133643d0
      rp = sma_osc * (1.d0 - e_osc)
      jd0 = 2458455.471721d0
      call jd2cal(jd0, cal0)

      n_force(1) = 1
      n_force(2) = 1
      n_force(3) = 1

c     Change ONLY these two lines between the pre-bulge and
c     post-bulge runs (values from rpe_campaign_weather.csv, norad
c     42928, zone 8):
      e_opt  = 0.051538d0
      bn_opt = 128.040d0

      a_kep(1) = rp / (1.d0 - e_opt)
      a_kep(2) = e_opt
      a_kep(3) = 19.094467d0
      a_kep(4) = 210.267993d0
      a_kep(5) = 73.147372d0
      a_kep(6) = 301.305632d0
      a_kep(7) = 0.d0

      call oe2car(a_kep, x0, xd0, 1.d-12)

      idump = 1
      call propagate_ks(
     &   x0, xd0, cal0,
     &   3000, 360, 1.d-10,
     &   n_force, 20, 2, 3,
     &   bn_opt, 1, 7.2921150d-5, 3.35281066d-3, 1.d0,
     &   1.2d0, 0.01d0, 1, 2,
     &   4.56d-6, 1.32712440018d11, 4.902801076d3,
     &   ALT_a, DEN_a, SCH_a, natm,
     &   maxp, idump,
     &   tj, tx, td,
     &   ec)

      write(*,'(A,I2,A,I8)') ' exit_code=', ec, '  pts=', idump

      jd_hit = 0.d0
      alt_prev = 1.d9
      do i = 1, idump
         r_mag = dsqrt(tx(1,i)**2 + tx(2,i)**2 + tx(3,i)**2)
         alt_km = r_mag - R_Earth
         if (alt_km .le. 102.d0 .and. jd_hit .eq. 0.d0) then
            jd_hit = tj(i)
            write(*,'(A,F10.3,A,F14.4)')
     &         '  crossed 102km at alt=', alt_km,
     &         '  JD=', jd_hit
         end if
         alt_prev = alt_km
      end do

      if (jd_hit .gt. 0.d0) then
         write(*,'(A,F10.3)') ' ELAPSED DAYS TO 102km = ',
     &      jd_hit - jd0
      else
         write(*,'(A,F10.3)') ' never reached 102km; final alt=',
     &      alt_prev
      end if

      stop
      end

c     Standalone check: TLE mean apogee vs true SGP4-osculating apogee
c     at the exact same TLE epoch, for one real 42928 TLE.
      program mean_vs_osc_check
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision r_j2k(3), v_j2k(3)
      double precision ecc0, ainc_deg, raan_deg, aop_deg, ama_deg
      double precision amm_revday, bstar
      integer ierr

      double precision sma_mean, ha_mean
      double precision r_mag, v_mag, hx, hy, hz, h_mag, p_sl
      double precision sma_osc, ecc_osc, ha_osc
      double precision ne_val

      call init_constants()

c     First TLE of input/example_42928.tle.txt (epoch 2017 day 243.632...)
      ainc_deg   = 19.1573d0
      raan_deg   = 127.9089d0
      ecc0       = 0.3280519d0
      aop_deg    = 165.9808d0
      ama_deg    = 200.4889d0
      amm_revday = 9.03267896d0
      bstar      = 0.0088244d0

c     --- "Mean apogee" the way tle_evolution.F computes it ---
      ne_val = amm_revday * 2.d0*pi / 86400.d0
      sma_mean = (amue*1.d9 / (ne_val*ne_val))**(1.d0/3.d0)
      sma_mean = sma_mean * 1.d-3
      ha_mean = sma_mean*(1.d0+ecc0) - R_Earth

c     --- True osculating apogee via SGP4 at tsince=0 (same epoch) ---
      call tle2sv(2017, 243.63239903d0, ainc_deg, raan_deg, ecc0,
     &            aop_deg, ama_deg, amm_revday, bstar,
     &            r_j2k, v_j2k, ierr)

      if (ierr .ne. 0) then
         write(*,*) 'tle2sv error, ierr=', ierr
         stop
      end if

      r_mag = dsqrt(r_j2k(1)**2 + r_j2k(2)**2 + r_j2k(3)**2)
      v_mag = dsqrt(v_j2k(1)**2 + v_j2k(2)**2 + v_j2k(3)**2)
      sma_osc = 1.d0 / (2.d0/r_mag - v_mag**2/amue)

      hx = r_j2k(2)*v_j2k(3) - r_j2k(3)*v_j2k(2)
      hy = r_j2k(3)*v_j2k(1) - r_j2k(1)*v_j2k(3)
      hz = r_j2k(1)*v_j2k(2) - r_j2k(2)*v_j2k(1)
      h_mag = dsqrt(hx**2 + hy**2 + hz**2)
      p_sl = h_mag**2 / amue
      ecc_osc = dsqrt(dabs(1.d0 - p_sl/sma_osc))

      ha_osc = sma_osc*(1.d0+ecc_osc) - R_Earth

      write(*,'(A,F14.6)') 'sma_mean (km)      = ', sma_mean
      write(*,'(A,F14.6)') 'sma_osc  (km)      = ', sma_osc
      write(*,'(A,F14.6)') 'ecc_mean           = ', ecc0
      write(*,'(A,F14.6)') 'ecc_osc            = ', ecc_osc
      write(*,'(A,F14.6)') 'ha_mean (km, alt)  = ', ha_mean
      write(*,'(A,F14.6)') 'ha_osc  (km, alt)  = ', ha_osc
      write(*,'(A,F14.6)') 'DIFFERENCE (km)    = ', ha_osc - ha_mean
      write(*,'(A,F10.4)') 'r_mag at epoch(km) = ', r_mag

      stop
      end

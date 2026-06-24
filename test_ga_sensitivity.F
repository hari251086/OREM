c     test_ga_sensitivity.F — GA parameter sensitivity study
c     -------------------------------------------------------
c     Runs the GA with varying parameters on synthetic surfaces
c     matching all 4 manuscript orbit regimes. Reports convergence
c     quality and optimal parameter recommendations.
c
c     Compile:
c       ifx test_ga_sensitivity.F ga.F /exe:test_ga_sensitivity.exe
c     -------------------------------------------------------
      program test_ga_sensitivity

      implicit double precision (a-h, o-z)

      parameter(ne = 3, na = 3, npts = 14, nob = 14)
      double precision surf(npts, ne, na)
      double precision eg(ne), ag(na)
      double precision tob(nob), apob(nob)
      double precision e_opt, a_opt, rms_opt

c     Test case parameters
      double precision e_true(4), ha_base(4), ha_slope(4)
      double precision e_delta(4), e_sens(4), a_sens(4)
      character*20 obj_name(4)
      integer ipops(3), igens(3)
      double precision pmuts(3)

c     Results storage
      double precision rms_pop(3, 4), rms_gen(3, 4), rms_pm(3, 4)

c     Orbit regimes
      obj_name(1) = '42928 PSLV-C39'
      e_true(1) = 0.3216d0
      ha_base(1) = 12000.d0
      ha_slope(1) = -100.d0
      e_delta(1) = 0.0004d0
      e_sens(1) = 10000.d0
      a_sens(1) = 0.5d0

      obj_name(2) = '37151 Long March 3B'
      e_true(2) = 0.5650d0
      ha_base(2) = 15000.d0
      ha_slope(2) = -150.d0
      e_delta(2) = 0.0005d0
      e_sens(2) = 15000.d0
      a_sens(2) = 8.d0

      obj_name(3) = '35497 Ariane 5'
      e_true(3) = 0.6295d0
      ha_base(3) = 35000.d0
      ha_slope(3) = -300.d0
      e_delta(3) = 0.0005d0
      e_sens(3) = 25000.d0
      a_sens(3) = 5.d0

      obj_name(4) = '39615 Proton-M'
      e_true(4) = 0.6808d0
      ha_base(4) = 27000.d0
      ha_slope(4) = -200.d0
      e_delta(4) = 0.0002d0
      e_sens(4) = 20000.d0
      a_sens(4) = 10.d0

c     Sweep values
      ipops(1) = 4
      ipops(2) = 8
      ipops(3) = 16

      igens(1) = 100
      igens(2) = 200
      igens(3) = 500

      pmuts(1) = 0.005d0
      pmuts(2) = 0.01d0
      pmuts(3) = 0.02d0

      ag(1) = 1.d0
      ag(2) = 3.d0
      ag(3) = 5.d0

      write(*,*) '=============================================='
      write(*,*) ' GA Parameter Sensitivity Study'
      write(*,*) '=============================================='

c     ==============================================================
c     Study 1: Population size sweep (gen=500, Pm=0.01)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  === Population Size Sweep ==='
      write(*,*) '  (gen=500, Pc=0.8, Pm=0.01)'
      write(*,*) ''
      write(*,'(A20,3A12)')
     &  'Object', 'pop=4', 'pop=8', 'pop=16'
      write(*,'(A20,3A12)')
     &  '--------------------','------------',
     &  '------------','------------'

      do iobj = 1, 4
         eg(1) = e_true(iobj) - e_delta(iobj)
         eg(2) = e_true(iobj)
         eg(3) = e_true(iobj) + e_delta(iobj)

         do k = 1, na
            do j = 1, ne
               do i = 1, npts
                  t = dble(i-1) / dble(npts-1)
                  surf(i,j,k) = ha_base(iobj) +
     &              ha_slope(iobj)*t
     &              - e_sens(iobj)*(eg(j) - e_true(iobj))
     &              - a_sens(iobj)*(ag(k) - 3.d0)
               end do
            end do
         end do
         do i = 1, nob
            tob(i) = dble(i-1) / dble(npts-1)
            apob(i) = ha_base(iobj) + ha_slope(iobj)*tob(i)
         end do

         do ip = 1, 3
            call ga_optimize(surf, eg, ag, ne, na, npts,
     &         tob, apob, nob,
     &         eg(1), eg(3), ag(1), ag(3),
     &         ipops(ip), 500, 40, 16,
     &         0.8d0, 0.01d0, 0.123d0,
     &         e_opt, a_opt, rms_opt)
            rms_pop(ip, iobj) = rms_opt
         end do
         write(*,'(A20,3F12.6)')
     &     obj_name(iobj),
     &     (rms_pop(ip, iobj), ip=1,3)
      end do

c     ==============================================================
c     Study 2: Generations sweep (pop=4, Pm=0.01)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  === Generations Sweep ==='
      write(*,*) '  (pop=4, Pc=0.8, Pm=0.01)'
      write(*,*) ''
      write(*,'(A20,3A12)')
     &  'Object', 'gen=100', 'gen=200', 'gen=500'
      write(*,'(A20,3A12)')
     &  '--------------------','------------',
     &  '------------','------------'

      do iobj = 1, 4
         eg(1) = e_true(iobj) - e_delta(iobj)
         eg(2) = e_true(iobj)
         eg(3) = e_true(iobj) + e_delta(iobj)

         do k = 1, na
            do j = 1, ne
               do i = 1, npts
                  t = dble(i-1) / dble(npts-1)
                  surf(i,j,k) = ha_base(iobj) +
     &              ha_slope(iobj)*t
     &              - e_sens(iobj)*(eg(j) - e_true(iobj))
     &              - a_sens(iobj)*(ag(k) - 3.d0)
               end do
            end do
         end do
         do i = 1, nob
            tob(i) = dble(i-1) / dble(npts-1)
            apob(i) = ha_base(iobj) + ha_slope(iobj)*tob(i)
         end do

         do ig = 1, 3
            call ga_optimize(surf, eg, ag, ne, na, npts,
     &         tob, apob, nob,
     &         eg(1), eg(3), ag(1), ag(3),
     &         4, igens(ig), 40, 16,
     &         0.8d0, 0.01d0, 0.123d0,
     &         e_opt, a_opt, rms_opt)
            rms_gen(ig, iobj) = rms_opt
         end do
         write(*,'(A20,3F12.6)')
     &     obj_name(iobj),
     &     (rms_gen(ig, iobj), ig=1,3)
      end do

c     ==============================================================
c     Study 3: Mutation rate sweep (pop=4, gen=500)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  === Mutation Rate Sweep ==='
      write(*,*) '  (pop=4, gen=500, Pc=0.8)'
      write(*,*) ''
      write(*,'(A20,3A12)')
     &  'Object', 'Pm=0.005', 'Pm=0.01', 'Pm=0.02'
      write(*,'(A20,3A12)')
     &  '--------------------','------------',
     &  '------------','------------'

      do iobj = 1, 4
         eg(1) = e_true(iobj) - e_delta(iobj)
         eg(2) = e_true(iobj)
         eg(3) = e_true(iobj) + e_delta(iobj)

         do k = 1, na
            do j = 1, ne
               do i = 1, npts
                  t = dble(i-1) / dble(npts-1)
                  surf(i,j,k) = ha_base(iobj) +
     &              ha_slope(iobj)*t
     &              - e_sens(iobj)*(eg(j) - e_true(iobj))
     &              - a_sens(iobj)*(ag(k) - 3.d0)
               end do
            end do
         end do
         do i = 1, nob
            tob(i) = dble(i-1) / dble(npts-1)
            apob(i) = ha_base(iobj) + ha_slope(iobj)*tob(i)
         end do

         do im = 1, 3
            call ga_optimize(surf, eg, ag, ne, na, npts,
     &         tob, apob, nob,
     &         eg(1), eg(3), ag(1), ag(3),
     &         4, 500, 40, 16,
     &         0.8d0, pmuts(im), 0.123d0,
     &         e_opt, a_opt, rms_opt)
            rms_pm(im, iobj) = rms_opt
         end do
         write(*,'(A20,3F12.6)')
     &     obj_name(iobj),
     &     (rms_pm(im, iobj), im=1,3)
      end do

c     ==============================================================
c     Summary and recommendation
c     ==============================================================
      write(*,*) ''
      write(*,*) '=============================================='
      write(*,*) ' Summary'
      write(*,*) '=============================================='

c     Check: default (pop=4, gen=500, Pm=0.01) achieves RMS < 0.1
c     for all objects
      iall_ok = 1
      do iobj = 1, 4
         if (rms_pop(1, iobj) .gt. 0.1d0) iall_ok = 0
      end do
      if (iall_ok .eq. 1) then
         write(*,*) '  Default params (4/500/0.01) achieve'
         write(*,*) '  RMS < 0.1 on all 4 objects: VALIDATED'
      else
         write(*,*) '  WARNING: Default params do not achieve'
         write(*,*) '  RMS < 0.1 on all objects'
      end if

c     Check: gen=200 sufficient?
      iall_ok = 1
      do iobj = 1, 4
         if (rms_gen(2, iobj) .gt. 0.1d0) iall_ok = 0
      end do
      if (iall_ok .eq. 1) then
         write(*,*) '  gen=200 sufficient: YES'
      else
         write(*,*) '  gen=200 sufficient: NO (need 500)'
      end if

c     Check: gen=100 sufficient?
      iall_ok = 1
      do iobj = 1, 4
         if (rms_gen(1, iobj) .gt. 0.1d0) iall_ok = 0
      end do
      if (iall_ok .eq. 1) then
         write(*,*) '  gen=100 sufficient: YES'
      else
         write(*,*) '  gen=100 sufficient: NO'
      end if

      write(*,*) ''

      stop
      end

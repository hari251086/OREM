c     ============================================================
c     screen_candidates.F -- issue: extend rpe_campaign to 100
c     objects. Runs the real tle_evolve -> tle_filter -> zone_select
c     chain (production params: min_zone_pts=8, max_zone_days=20,
c     r2_thresh=0.90, slope_thresh=-1.0, nzones_max=8) against each
c     candidate's fetched TLE history, to get the actual pass/fail
c     and zone count -- not a raw-TLE-count guess.
c
c     Not wired into fpm.toml. To reproduce, temporarily add:
c       [[executable]]
c       name = "screen_candidates"
c       source-dir = "scratch_rpe/bench/campaign100"
c       main = "screen_candidates.F"
c     then: fpm run screen_candidates --compiler ifx --flag "-heap-arrays"
c     (edit cand_norad's DATA statement + the output filename per round
c     -- this file was reused for 3 screening rounds, results in
c     screen_results.csv / screen_results2.csv / screen_results3.csv;
c     see rpe_campaign.F's own header for the full sourcing writeup
c     and which candidates each round's results fed into the final
c     50->97 object extension, 2026-08-01)
c     ============================================================
      program screen_candidates
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(maxpts = 15000, mxz = 8)
      double precision epochs(maxpts), a_tle(maxpts)
      double precision e_tle(maxpts), inc_tle(maxpts)
      double precision raan_tle(maxpts), aop_tle(maxpts)
      double precision ha_tle(maxpts), hp_tle(maxpts)
      double precision ls_tle(maxpts)
      integer npts, ierr_tle

      integer mask_tle(maxpts)
      integer nfiltered_tle, nremoved_tle, ngaps_tle

      integer zone_start(mxz), zone_end(mxz), nzones

      integer ncand, ic, norad
      character*80 tle_file

      parameter(ncand = 9)
      integer cand_norad(ncand)
      data cand_norad /
     &   41686, 42977, 42979, 43685, 44555, 45494, 45495, 45497,
     &   45499 /

      call init_constants()

      open(unit=97,
     &   file='scratch_rpe/bench/campaign100/screen_results3.csv',
     &   status='replace')
      write(97,'(A)') 'norad,npts,nfiltered,nremoved,nzones,ierr'

      do ic = 1, ncand
         norad = cand_norad(ic)
         write(tle_file, '(A,I0,A,I0,A)')
     &      'E:\Research\1. R&D\Re-entry\2026\Data\objects\',
     &      norad, '\', norad, '.tle.txt'

         call tle_evolve(tle_file, norad, maxpts,
     &      epochs, a_tle, e_tle, inc_tle,
     &      raan_tle, aop_tle,
     &      ha_tle, hp_tle, ls_tle,
     &      npts, ierr_tle)

         if (ierr_tle .ne. 0) then
            write(97,'(I6,A,I2)') norad, ',0,0,0,0,', ierr_tle
            write(*,'(I6,A,I2)') norad, ': tle_evolve ierr=', ierr_tle
            goto 50
         end if

         call tle_filter(epochs, ha_tle, e_tle, npts,
     &      4.0d0, 6.0d0, 30.0d0, 20,
     &      mask_tle, nfiltered_tle, nremoved_tle, ngaps_tle)

         if (nremoved_tle .gt. 0) then
            j_filt = 0
            do i_filt = 1, npts
               if (mask_tle(i_filt) .eq. 1) then
                  j_filt = j_filt + 1
                  epochs(j_filt) = epochs(i_filt)
                  ha_tle(j_filt) = ha_tle(i_filt)
               end if
            end do
            npts = j_filt
         end if

         call zone_select(epochs, ha_tle, npts,
     &      8, 20.d0, 0.90d0, -1.d0,
     &      mxz,
     &      zone_start, zone_end, nzones)

         write(97,'(I6,A,I6,A,I6,A,I6,A,I2,A,I2)')
     &      norad, ',', npts, ',', nfiltered_tle, ',',
     &      nremoved_tle, ',', nzones, ',0'
         write(*,'(I6,A,I6,A,I2)')
     &      norad, ': npts=', npts, ' nzones=', nzones

 50      continue
      end do

      close(97)
      write(*,*) 'Done -- see screen_results.csv'
      stop
      end

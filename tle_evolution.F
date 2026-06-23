c     ============================================================
c     tle_evolution.F — Batch TLE processing for orbital evolution
c
c     Reads a TLE file, filters by NORAD ID, and extracts the
c     orbital evolution time-series: epoch (JD), semi-major axis,
c     eccentricity, inclination, RAAN, argument of perigee,
c     mean apogee/perigee altitudes, and Sun azimuth angle.
c
c     Uses TLE mean elements directly (no osculating conversion
c     needed for evolution tracking). Sun azimuth computed via
c     solarnpv + sun_azimuth from Subrouts.F.
c
c     Compile:
c       ifx tle_evolution.F ksrop/Subrouts.F ksrop/TLEread.F
c           ksrop/Legendre.F
c
c     Author: Harishkumar Sellamuthu
c     ============================================================

c     ============================================================
      subroutine tle_evolve(tle_file, target_norad, maxpts,
     &                      epochs, a_out, e_out, inc_out,
     &                      raan_out, aop_out,
     &                      ha_out, hp_out, lambda_s_out,
     &                      npts, ierr)
c
c     Process TLE history for a single NORAD ID and extract
c     the orbital evolution time-series.
c
c     Input
c       tle_file     : character*(*) — path to TLE file
c       target_norad : integer — NORAD catalog number to filter
c                      (-1 to accept all entries)
c       maxpts       : integer — dimension of output arrays
c
c     Output (parallel arrays, one element per TLE entry)
c       epochs(maxpts)    : Julian dates (UTC)
c       a_out(maxpts)     : semi-major axis (km)
c       e_out(maxpts)     : eccentricity
c       inc_out(maxpts)   : inclination (deg)
c       raan_out(maxpts)  : RAAN (deg)
c       aop_out(maxpts)   : argument of perigee (deg)
c       ha_out(maxpts)    : mean apogee altitude (km)
c       hp_out(maxpts)    : mean perigee altitude (km)
c       lambda_s_out(maxpts): Sun azimuth angle (deg)
c       npts              : number of entries actually stored
c       ierr              : 0=ok, 1=file error, 2=no matching
c     ============================================================
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      character*(*) tle_file
      integer target_norad, maxpts, npts, ierr

      double precision epochs(maxpts)
      double precision a_out(maxpts), e_out(maxpts)
      double precision inc_out(maxpts)
      double precision raan_out(maxpts), aop_out(maxpts)
      double precision ha_out(maxpts), hp_out(maxpts)
      double precision lambda_s_out(maxpts)

c     TLE read_tle arrays — sized for large catalogs
      parameter(maxtle = 100000)
      integer norad(maxtle), epoch_yr(maxtle), revnum(maxtle)
      double precision epoch_day(maxtle)
      double precision ndot_a(maxtle), nddot_a(maxtle)
      double precision bstar_a(maxtle)
      double precision ainc_a(maxtle), raan_a(maxtle)
      double precision ecc_a(maxtle)
      double precision aop_a(maxtle), ama_a(maxtle)
      double precision amm_a(maxtle)
      character*24 names(maxtle)
      character*5  catid(maxtle)
      integer ntle

c     Local working variables
      double precision djd, sma, ha, hp
      double precision ts(3), sun_ra, sun_dec, sun_r
      double precision lambda_s
      double precision twopi_val, ne_val
      integer i, iunit

      ierr = 0
      npts = 0
      twopi_val = 2.d0 * pi

c     --- Read all TLEs from file ---
      iunit = 77
      open(unit=iunit, file=tle_file, status='old', err=900)
      call read_tle(iunit, maxtle, names, catid, norad,
     &             epoch_yr, epoch_day,
     &             ndot_a, nddot_a, bstar_a,
     &             ainc_a, raan_a, ecc_a,
     &             aop_a, ama_a, amm_a,
     &             revnum, ntle)
      close(iunit)

      if (ntle .eq. 0) then
         ierr = 2
         return
      end if

c     --- Filter by NORAD ID and extract evolution ---
      do i = 1, ntle
         if (target_norad .gt. 0 .and.
     &       norad(i) .ne. target_norad) goto 20

         if (npts .ge. maxpts) goto 30

         npts = npts + 1

c        Epoch → Julian Date
         call tle_epoch2jd(epoch_yr(i), epoch_day(i), djd)
         epochs(npts) = djd

c        Mean motion (rev/day) → semi-major axis (km)
c        a = (mu / n²)^(1/3) where n is in rad/s
c        n_rad_s = amm * 2*pi / 86400
         ne_val = amm_a(i) * twopi_val / 86400.d0
         sma = (amue*1.d9 / (ne_val*ne_val))**(1.d0/3.d0)
         sma = sma * 1.d-3

         a_out(npts)    = sma
         e_out(npts)    = ecc_a(i)
         inc_out(npts)  = ainc_a(i)
         raan_out(npts) = raan_a(i)
         aop_out(npts)  = aop_a(i)

c        Mean apogee/perigee altitude (km)
         ha = sma * (1.d0 + ecc_a(i)) - R_Earth
         hp = sma * (1.d0 - ecc_a(i)) - R_Earth
         ha_out(npts) = ha
         hp_out(npts) = hp

c        Sun azimuth angle
c        Get Sun position at this epoch
         call solarnpv(djd, ts)

c        Sun RA and Dec from position vector
         sun_r = dsqrt(ts(1)**2 + ts(2)**2 + ts(3)**2)
         sun_ra  = datan2(ts(2), ts(1))
         if (sun_ra .lt. 0.d0) sun_ra = sun_ra + twopi_val
         sun_dec = dasin(ts(3) / sun_r)

c        Compute Sun azimuth w.r.t. orbital plane
         call sun_azimuth(
     &        ainc_a(i)*d2r, aop_a(i)*d2r, raan_a(i)*d2r,
     &        sun_ra, sun_dec, lambda_s)

         lambda_s_out(npts) = lambda_s * r2d

 20      continue
      end do

 30   continue

      if (npts .eq. 0) then
         ierr = 2
      end if

      return

c     File open error
 900  ierr = 1
      return
      end

c     validation_rpe_campaign.F -- Validation RPE campaign for the 7
c     confirmed/plausible-genuine solar-apsidal-resonance objects from
c     issue #56's manual review (35497, 37151, 27526, 59347, 40943,
c     10724, 18954). Runs OREM's standard orem_run pipeline once per
c     PRE-SELECTED zone (340 total, from
c     scratch_rpe/validation_rpe_zones.py + prepare_validation_zones.py
c     -- 3-29 TLEs/zone, first zone >=6mo before each object's own
c     U-turn epoch, spread through to the last tracked TLE, matching
c     the "E:\Research\1. R&D\Re-entry\COSPAR ASR" reference's own
c     per-zone TLE-split style), rather than letting zone_select.F
c     pick zones itself.
c
c     GENERATED FILE -- do not hand-edit. Regenerate via
c     `python scratch_rpe/prepare_validation_zones.py` after any change
c     to scratch_rpe/validation_rpe_zones.py's OBJECTS/zone logic.
c
c     Each "object" below is really one (norad, zone) pair -- its own
c     TLE input file already contains ONLY that zone's TLE records
c     (input/example_<norad>_valzoneN.tle.txt). min_zone_pts=3 (down
c     from production's 8, since 3 is this campaign's own floor) and
c     r2_thresh/slope_thresh are maximally permissive (0.0 / +1.0d6)
c     so zone_select.F's decay-linearity filter -- built for the
c     terminal monotonic-decay regime -- doesn't reject the earlier,
c     still-oscillating pre-resonance windows this campaign
c     deliberately includes; the externally-chosen window is trusted
c     as-is, not re-filtered. nzones_max=1 since each input file is
c     already exactly one zone. max_zone_days=900 covers the widest
c     zone (10724's sparse-tracking zone 1, ~850 real days for just 29
c     points). All other GA/force-model parameters are identical to
c     the main scratch_rpe/rpe_campaign.F 97-object campaign for
c     direct comparability.
      program validation_rpe_campaign

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      integer io_start, io_end, iarg
      character*80 csv_name
      character*16 arg_str

      parameter(ndim_a = 291, mxz = 1, nobj = 340)

      double precision ALT_a(ndim_a), DEN_a(ndim_a), SCH_a(ndim_a)
      integer natm

      double precision reentry(mxz), e_opt(mxz), bn_opt(mxz)
      double precision rms_o(mxz), zepoch(mxz), rpe(mxz)
      double precision t_obs(6), t_mean, t_std
      integer nzu, ierr, nzv
      integer zstat(mxz)
      double precision bnb_t(10), tb_t(10)
      integer nbv_t

      character*54 tle_file(nobj)
      character*24 obj_name(nobj)
      integer norad(nobj)
      double precision oyr(nobj), omo(nobj), ody(nobj)
      double precision t_obs_jd, horiz, rpe_ens
      integer io, iz, nre

c     norad/oyr/omo/ody DATA statements
      data norad /
     &   35497, 35497, 35497, 35497, 35497, 35497,
     &   35497, 35497, 35497, 35497, 35497, 35497,
     &   35497, 35497, 35497, 35497, 35497, 35497,
     &   35497, 35497, 35497, 35497, 35497, 35497,
     &   35497, 35497, 35497, 35497, 35497, 35497,
     &   35497, 35497, 35497, 35497, 35497, 35497,
     &   35497, 35497, 35497, 35497, 35497, 35497,
     &   35497, 35497, 35497, 35497, 35497, 35497,
     &   35497, 35497, 35497, 35497, 35497, 35497,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 37151, 37151, 37151, 37151, 37151,
     &   37151, 27526, 27526, 27526, 27526, 27526,
     &   27526, 27526, 27526, 27526, 27526, 27526,
     &   27526, 27526, 27526, 27526, 27526, 27526,
     &   27526, 27526, 27526, 27526, 27526, 27526,
     &   27526, 27526, 27526, 27526, 27526, 27526,
     &   27526, 27526, 27526, 27526, 27526, 27526,
     &   27526, 27526, 27526, 27526, 27526, 27526,
     &   27526, 27526, 27526, 27526, 27526, 27526,
     &   59347, 59347, 59347, 59347, 59347, 59347,
     &   59347, 59347, 59347, 59347, 59347, 59347,
     &   59347, 59347, 59347, 59347, 59347, 59347,
     &   59347, 59347, 59347, 59347, 59347, 59347,
     &   59347, 59347, 59347, 59347, 59347, 59347,
     &   59347, 59347, 59347, 59347, 59347, 59347,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 40943,
     &   40943, 40943, 40943, 40943, 40943, 10724,
     &   10724, 10724, 10724, 18954 /
      data oyr /
     &   2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0,
     &   2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0,
     &   2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0,
     &   2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0,
     &   2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0,
     &   2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0,
     &   2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0,
     &   2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0,
     &   2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0, 2016.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0, 2015.d0,
     &   2015.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0,
     &   2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0,
     &   2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0,
     &   2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0,
     &   2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0,
     &   2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0,
     &   2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0,
     &   2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0, 2012.d0,
     &   2026.d0, 2026.d0, 2026.d0, 2026.d0, 2026.d0, 2026.d0,
     &   2026.d0, 2026.d0, 2026.d0, 2026.d0, 2026.d0, 2026.d0,
     &   2026.d0, 2026.d0, 2026.d0, 2026.d0, 2026.d0, 2026.d0,
     &   2026.d0, 2026.d0, 2026.d0, 2026.d0, 2026.d0, 2026.d0,
     &   2026.d0, 2026.d0, 2026.d0, 2026.d0, 2026.d0, 2026.d0,
     &   2026.d0, 2026.d0, 2026.d0, 2026.d0, 2026.d0, 2026.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0,
     &   2020.d0, 2020.d0, 2020.d0, 2020.d0, 2020.d0, 1982.d0,
     &   1982.d0, 1982.d0, 1982.d0, 1989.d0 /
      data omo /
     &   10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0,
     &   10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0,
     &   10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0,
     &   10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0,
     &   10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0,
     &   10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0,
     &   10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0,
     &   10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0,
     &   10.d0, 10.d0, 10.d0, 10.d0, 10.d0, 10.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 12.d0, 12.d0, 12.d0, 12.d0, 12.d0,
     &   12.d0, 5.d0, 5.d0, 5.d0, 5.d0, 5.d0,
     &   5.d0, 5.d0, 5.d0, 5.d0, 5.d0, 5.d0,
     &   5.d0, 5.d0, 5.d0, 5.d0, 5.d0, 5.d0,
     &   5.d0, 5.d0, 5.d0, 5.d0, 5.d0, 5.d0,
     &   5.d0, 5.d0, 5.d0, 5.d0, 5.d0, 5.d0,
     &   5.d0, 5.d0, 5.d0, 5.d0, 5.d0, 5.d0,
     &   5.d0, 5.d0, 5.d0, 5.d0, 5.d0, 5.d0,
     &   5.d0, 5.d0, 5.d0, 5.d0, 5.d0, 5.d0,
     &   6.d0, 6.d0, 6.d0, 6.d0, 6.d0, 6.d0,
     &   6.d0, 6.d0, 6.d0, 6.d0, 6.d0, 6.d0,
     &   6.d0, 6.d0, 6.d0, 6.d0, 6.d0, 6.d0,
     &   6.d0, 6.d0, 6.d0, 6.d0, 6.d0, 6.d0,
     &   6.d0, 6.d0, 6.d0, 6.d0, 6.d0, 6.d0,
     &   6.d0, 6.d0, 6.d0, 6.d0, 6.d0, 6.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0,
     &   2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 8.d0,
     &   8.d0, 8.d0, 8.d0, 2.d0 /
      data ody /
     &   31.d0, 31.d0, 31.d0, 31.d0, 31.d0, 31.d0,
     &   31.d0, 31.d0, 31.d0, 31.d0, 31.d0, 31.d0,
     &   31.d0, 31.d0, 31.d0, 31.d0, 31.d0, 31.d0,
     &   31.d0, 31.d0, 31.d0, 31.d0, 31.d0, 31.d0,
     &   31.d0, 31.d0, 31.d0, 31.d0, 31.d0, 31.d0,
     &   31.d0, 31.d0, 31.d0, 31.d0, 31.d0, 31.d0,
     &   31.d0, 31.d0, 31.d0, 31.d0, 31.d0, 31.d0,
     &   31.d0, 31.d0, 31.d0, 31.d0, 31.d0, 31.d0,
     &   31.d0, 31.d0, 31.d0, 31.d0, 31.d0, 31.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 3.d0, 3.d0, 3.d0, 3.d0, 3.d0,
     &   3.d0, 9.d0, 9.d0, 9.d0, 9.d0, 9.d0,
     &   9.d0, 9.d0, 9.d0, 9.d0, 9.d0, 9.d0,
     &   9.d0, 9.d0, 9.d0, 9.d0, 9.d0, 9.d0,
     &   9.d0, 9.d0, 9.d0, 9.d0, 9.d0, 9.d0,
     &   9.d0, 9.d0, 9.d0, 9.d0, 9.d0, 9.d0,
     &   9.d0, 9.d0, 9.d0, 9.d0, 9.d0, 9.d0,
     &   9.d0, 9.d0, 9.d0, 9.d0, 9.d0, 9.d0,
     &   9.d0, 9.d0, 9.d0, 9.d0, 9.d0, 9.d0,
     &   8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 8.d0,
     &   8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 8.d0,
     &   8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 8.d0,
     &   8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 8.d0,
     &   8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 8.d0,
     &   8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 8.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 20.d0,
     &   20.d0, 20.d0, 20.d0, 20.d0, 20.d0, 31.d0,
     &   31.d0, 31.d0, 31.d0, 1.d0 /

c     tle_file/obj_name string assignments (executable statements)
      tle_file(1) = 'input/example_35497_valzone1.tle.txt'
      tle_file(2) = 'input/example_35497_valzone2.tle.txt'
      tle_file(3) = 'input/example_35497_valzone3.tle.txt'
      tle_file(4) = 'input/example_35497_valzone4.tle.txt'
      tle_file(5) = 'input/example_35497_valzone5.tle.txt'
      tle_file(6) = 'input/example_35497_valzone6.tle.txt'
      tle_file(7) = 'input/example_35497_valzone7.tle.txt'
      tle_file(8) = 'input/example_35497_valzone8.tle.txt'
      tle_file(9) = 'input/example_35497_valzone9.tle.txt'
      tle_file(10) = 'input/example_35497_valzone10.tle.txt'
      tle_file(11) = 'input/example_35497_valzone11.tle.txt'
      tle_file(12) = 'input/example_35497_valzone12.tle.txt'
      tle_file(13) = 'input/example_35497_valzone13.tle.txt'
      tle_file(14) = 'input/example_35497_valzone14.tle.txt'
      tle_file(15) = 'input/example_35497_valzone15.tle.txt'
      tle_file(16) = 'input/example_35497_valzone16.tle.txt'
      tle_file(17) = 'input/example_35497_valzone17.tle.txt'
      tle_file(18) = 'input/example_35497_valzone18.tle.txt'
      tle_file(19) = 'input/example_35497_valzone19.tle.txt'
      tle_file(20) = 'input/example_35497_valzone20.tle.txt'
      tle_file(21) = 'input/example_35497_valzone21.tle.txt'
      tle_file(22) = 'input/example_35497_valzone22.tle.txt'
      tle_file(23) = 'input/example_35497_valzone23.tle.txt'
      tle_file(24) = 'input/example_35497_valzone24.tle.txt'
      tle_file(25) = 'input/example_35497_valzone25.tle.txt'
      tle_file(26) = 'input/example_35497_valzone26.tle.txt'
      tle_file(27) = 'input/example_35497_valzone27.tle.txt'
      tle_file(28) = 'input/example_35497_valzone28.tle.txt'
      tle_file(29) = 'input/example_35497_valzone29.tle.txt'
      tle_file(30) = 'input/example_35497_valzone30.tle.txt'
      tle_file(31) = 'input/example_35497_valzone31.tle.txt'
      tle_file(32) = 'input/example_35497_valzone32.tle.txt'
      tle_file(33) = 'input/example_35497_valzone33.tle.txt'
      tle_file(34) = 'input/example_35497_valzone34.tle.txt'
      tle_file(35) = 'input/example_35497_valzone35.tle.txt'
      tle_file(36) = 'input/example_35497_valzone36.tle.txt'
      tle_file(37) = 'input/example_35497_valzone37.tle.txt'
      tle_file(38) = 'input/example_35497_valzone38.tle.txt'
      tle_file(39) = 'input/example_35497_valzone39.tle.txt'
      tle_file(40) = 'input/example_35497_valzone40.tle.txt'
      tle_file(41) = 'input/example_35497_valzone41.tle.txt'
      tle_file(42) = 'input/example_35497_valzone42.tle.txt'
      tle_file(43) = 'input/example_35497_valzone43.tle.txt'
      tle_file(44) = 'input/example_35497_valzone44.tle.txt'
      tle_file(45) = 'input/example_35497_valzone45.tle.txt'
      tle_file(46) = 'input/example_35497_valzone46.tle.txt'
      tle_file(47) = 'input/example_35497_valzone47.tle.txt'
      tle_file(48) = 'input/example_35497_valzone48.tle.txt'
      tle_file(49) = 'input/example_35497_valzone49.tle.txt'
      tle_file(50) = 'input/example_35497_valzone50.tle.txt'
      tle_file(51) = 'input/example_35497_valzone51.tle.txt'
      tle_file(52) = 'input/example_35497_valzone52.tle.txt'
      tle_file(53) = 'input/example_35497_valzone53.tle.txt'
      tle_file(54) = 'input/example_35497_valzone54.tle.txt'
      tle_file(55) = 'input/example_37151_valzone1.tle.txt'
      tle_file(56) = 'input/example_37151_valzone2.tle.txt'
      tle_file(57) = 'input/example_37151_valzone3.tle.txt'
      tle_file(58) = 'input/example_37151_valzone4.tle.txt'
      tle_file(59) = 'input/example_37151_valzone5.tle.txt'
      tle_file(60) = 'input/example_37151_valzone6.tle.txt'
      tle_file(61) = 'input/example_37151_valzone7.tle.txt'
      tle_file(62) = 'input/example_37151_valzone8.tle.txt'
      tle_file(63) = 'input/example_37151_valzone9.tle.txt'
      tle_file(64) = 'input/example_37151_valzone10.tle.txt'
      tle_file(65) = 'input/example_37151_valzone11.tle.txt'
      tle_file(66) = 'input/example_37151_valzone12.tle.txt'
      tle_file(67) = 'input/example_37151_valzone13.tle.txt'
      tle_file(68) = 'input/example_37151_valzone14.tle.txt'
      tle_file(69) = 'input/example_37151_valzone15.tle.txt'
      tle_file(70) = 'input/example_37151_valzone16.tle.txt'
      tle_file(71) = 'input/example_37151_valzone17.tle.txt'
      tle_file(72) = 'input/example_37151_valzone18.tle.txt'
      tle_file(73) = 'input/example_37151_valzone19.tle.txt'
      tle_file(74) = 'input/example_37151_valzone20.tle.txt'
      tle_file(75) = 'input/example_37151_valzone21.tle.txt'
      tle_file(76) = 'input/example_37151_valzone22.tle.txt'
      tle_file(77) = 'input/example_37151_valzone23.tle.txt'
      tle_file(78) = 'input/example_37151_valzone24.tle.txt'
      tle_file(79) = 'input/example_37151_valzone25.tle.txt'
      tle_file(80) = 'input/example_37151_valzone26.tle.txt'
      tle_file(81) = 'input/example_37151_valzone27.tle.txt'
      tle_file(82) = 'input/example_37151_valzone28.tle.txt'
      tle_file(83) = 'input/example_37151_valzone29.tle.txt'
      tle_file(84) = 'input/example_37151_valzone30.tle.txt'
      tle_file(85) = 'input/example_37151_valzone31.tle.txt'
      tle_file(86) = 'input/example_37151_valzone32.tle.txt'
      tle_file(87) = 'input/example_37151_valzone33.tle.txt'
      tle_file(88) = 'input/example_37151_valzone34.tle.txt'
      tle_file(89) = 'input/example_37151_valzone35.tle.txt'
      tle_file(90) = 'input/example_37151_valzone36.tle.txt'
      tle_file(91) = 'input/example_37151_valzone37.tle.txt'
      tle_file(92) = 'input/example_37151_valzone38.tle.txt'
      tle_file(93) = 'input/example_37151_valzone39.tle.txt'
      tle_file(94) = 'input/example_37151_valzone40.tle.txt'
      tle_file(95) = 'input/example_37151_valzone41.tle.txt'
      tle_file(96) = 'input/example_37151_valzone42.tle.txt'
      tle_file(97) = 'input/example_37151_valzone43.tle.txt'
      tle_file(98) = 'input/example_37151_valzone44.tle.txt'
      tle_file(99) = 'input/example_37151_valzone45.tle.txt'
      tle_file(100) = 'input/example_37151_valzone46.tle.txt'
      tle_file(101) = 'input/example_37151_valzone47.tle.txt'
      tle_file(102) = 'input/example_37151_valzone48.tle.txt'
      tle_file(103) = 'input/example_37151_valzone49.tle.txt'
      tle_file(104) = 'input/example_37151_valzone50.tle.txt'
      tle_file(105) = 'input/example_37151_valzone51.tle.txt'
      tle_file(106) = 'input/example_37151_valzone52.tle.txt'
      tle_file(107) = 'input/example_37151_valzone53.tle.txt'
      tle_file(108) = 'input/example_37151_valzone54.tle.txt'
      tle_file(109) = 'input/example_37151_valzone55.tle.txt'
      tle_file(110) = 'input/example_37151_valzone56.tle.txt'
      tle_file(111) = 'input/example_37151_valzone57.tle.txt'
      tle_file(112) = 'input/example_37151_valzone58.tle.txt'
      tle_file(113) = 'input/example_37151_valzone59.tle.txt'
      tle_file(114) = 'input/example_37151_valzone60.tle.txt'
      tle_file(115) = 'input/example_37151_valzone61.tle.txt'
      tle_file(116) = 'input/example_37151_valzone62.tle.txt'
      tle_file(117) = 'input/example_37151_valzone63.tle.txt'
      tle_file(118) = 'input/example_37151_valzone64.tle.txt'
      tle_file(119) = 'input/example_37151_valzone65.tle.txt'
      tle_file(120) = 'input/example_37151_valzone66.tle.txt'
      tle_file(121) = 'input/example_37151_valzone67.tle.txt'
      tle_file(122) = 'input/example_37151_valzone68.tle.txt'
      tle_file(123) = 'input/example_37151_valzone69.tle.txt'
      tle_file(124) = 'input/example_37151_valzone70.tle.txt'
      tle_file(125) = 'input/example_37151_valzone71.tle.txt'
      tle_file(126) = 'input/example_37151_valzone72.tle.txt'
      tle_file(127) = 'input/example_37151_valzone73.tle.txt'
      tle_file(128) = 'input/example_37151_valzone74.tle.txt'
      tle_file(129) = 'input/example_37151_valzone75.tle.txt'
      tle_file(130) = 'input/example_37151_valzone76.tle.txt'
      tle_file(131) = 'input/example_37151_valzone77.tle.txt'
      tle_file(132) = 'input/example_37151_valzone78.tle.txt'
      tle_file(133) = 'input/example_37151_valzone79.tle.txt'
      tle_file(134) = 'input/example_37151_valzone80.tle.txt'
      tle_file(135) = 'input/example_37151_valzone81.tle.txt'
      tle_file(136) = 'input/example_37151_valzone82.tle.txt'
      tle_file(137) = 'input/example_37151_valzone83.tle.txt'
      tle_file(138) = 'input/example_37151_valzone84.tle.txt'
      tle_file(139) = 'input/example_37151_valzone85.tle.txt'
      tle_file(140) = 'input/example_37151_valzone86.tle.txt'
      tle_file(141) = 'input/example_37151_valzone87.tle.txt'
      tle_file(142) = 'input/example_37151_valzone88.tle.txt'
      tle_file(143) = 'input/example_37151_valzone89.tle.txt'
      tle_file(144) = 'input/example_37151_valzone90.tle.txt'
      tle_file(145) = 'input/example_37151_valzone91.tle.txt'
      tle_file(146) = 'input/example_37151_valzone92.tle.txt'
      tle_file(147) = 'input/example_37151_valzone93.tle.txt'
      tle_file(148) = 'input/example_37151_valzone94.tle.txt'
      tle_file(149) = 'input/example_37151_valzone95.tle.txt'
      tle_file(150) = 'input/example_37151_valzone96.tle.txt'
      tle_file(151) = 'input/example_37151_valzone97.tle.txt'
      tle_file(152) = 'input/example_37151_valzone98.tle.txt'
      tle_file(153) = 'input/example_37151_valzone99.tle.txt'
      tle_file(154) = 'input/example_37151_valzone100.tle.txt'
      tle_file(155) = 'input/example_37151_valzone101.tle.txt'
      tle_file(156) = 'input/example_37151_valzone102.tle.txt'
      tle_file(157) = 'input/example_37151_valzone103.tle.txt'
      tle_file(158) = 'input/example_27526_valzone1.tle.txt'
      tle_file(159) = 'input/example_27526_valzone2.tle.txt'
      tle_file(160) = 'input/example_27526_valzone3.tle.txt'
      tle_file(161) = 'input/example_27526_valzone4.tle.txt'
      tle_file(162) = 'input/example_27526_valzone5.tle.txt'
      tle_file(163) = 'input/example_27526_valzone6.tle.txt'
      tle_file(164) = 'input/example_27526_valzone7.tle.txt'
      tle_file(165) = 'input/example_27526_valzone8.tle.txt'
      tle_file(166) = 'input/example_27526_valzone9.tle.txt'
      tle_file(167) = 'input/example_27526_valzone10.tle.txt'
      tle_file(168) = 'input/example_27526_valzone11.tle.txt'
      tle_file(169) = 'input/example_27526_valzone12.tle.txt'
      tle_file(170) = 'input/example_27526_valzone13.tle.txt'
      tle_file(171) = 'input/example_27526_valzone14.tle.txt'
      tle_file(172) = 'input/example_27526_valzone15.tle.txt'
      tle_file(173) = 'input/example_27526_valzone16.tle.txt'
      tle_file(174) = 'input/example_27526_valzone17.tle.txt'
      tle_file(175) = 'input/example_27526_valzone18.tle.txt'
      tle_file(176) = 'input/example_27526_valzone19.tle.txt'
      tle_file(177) = 'input/example_27526_valzone20.tle.txt'
      tle_file(178) = 'input/example_27526_valzone21.tle.txt'
      tle_file(179) = 'input/example_27526_valzone22.tle.txt'
      tle_file(180) = 'input/example_27526_valzone23.tle.txt'
      tle_file(181) = 'input/example_27526_valzone24.tle.txt'
      tle_file(182) = 'input/example_27526_valzone25.tle.txt'
      tle_file(183) = 'input/example_27526_valzone26.tle.txt'
      tle_file(184) = 'input/example_27526_valzone27.tle.txt'
      tle_file(185) = 'input/example_27526_valzone28.tle.txt'
      tle_file(186) = 'input/example_27526_valzone29.tle.txt'
      tle_file(187) = 'input/example_27526_valzone30.tle.txt'
      tle_file(188) = 'input/example_27526_valzone31.tle.txt'
      tle_file(189) = 'input/example_27526_valzone32.tle.txt'
      tle_file(190) = 'input/example_27526_valzone33.tle.txt'
      tle_file(191) = 'input/example_27526_valzone34.tle.txt'
      tle_file(192) = 'input/example_27526_valzone35.tle.txt'
      tle_file(193) = 'input/example_27526_valzone36.tle.txt'
      tle_file(194) = 'input/example_27526_valzone37.tle.txt'
      tle_file(195) = 'input/example_27526_valzone38.tle.txt'
      tle_file(196) = 'input/example_27526_valzone39.tle.txt'
      tle_file(197) = 'input/example_27526_valzone40.tle.txt'
      tle_file(198) = 'input/example_27526_valzone41.tle.txt'
      tle_file(199) = 'input/example_27526_valzone42.tle.txt'
      tle_file(200) = 'input/example_27526_valzone43.tle.txt'
      tle_file(201) = 'input/example_27526_valzone44.tle.txt'
      tle_file(202) = 'input/example_27526_valzone45.tle.txt'
      tle_file(203) = 'input/example_27526_valzone46.tle.txt'
      tle_file(204) = 'input/example_27526_valzone47.tle.txt'
      tle_file(205) = 'input/example_59347_valzone1.tle.txt'
      tle_file(206) = 'input/example_59347_valzone2.tle.txt'
      tle_file(207) = 'input/example_59347_valzone3.tle.txt'
      tle_file(208) = 'input/example_59347_valzone4.tle.txt'
      tle_file(209) = 'input/example_59347_valzone5.tle.txt'
      tle_file(210) = 'input/example_59347_valzone6.tle.txt'
      tle_file(211) = 'input/example_59347_valzone7.tle.txt'
      tle_file(212) = 'input/example_59347_valzone8.tle.txt'
      tle_file(213) = 'input/example_59347_valzone9.tle.txt'
      tle_file(214) = 'input/example_59347_valzone10.tle.txt'
      tle_file(215) = 'input/example_59347_valzone11.tle.txt'
      tle_file(216) = 'input/example_59347_valzone12.tle.txt'
      tle_file(217) = 'input/example_59347_valzone13.tle.txt'
      tle_file(218) = 'input/example_59347_valzone14.tle.txt'
      tle_file(219) = 'input/example_59347_valzone15.tle.txt'
      tle_file(220) = 'input/example_59347_valzone16.tle.txt'
      tle_file(221) = 'input/example_59347_valzone17.tle.txt'
      tle_file(222) = 'input/example_59347_valzone18.tle.txt'
      tle_file(223) = 'input/example_59347_valzone19.tle.txt'
      tle_file(224) = 'input/example_59347_valzone20.tle.txt'
      tle_file(225) = 'input/example_59347_valzone21.tle.txt'
      tle_file(226) = 'input/example_59347_valzone22.tle.txt'
      tle_file(227) = 'input/example_59347_valzone23.tle.txt'
      tle_file(228) = 'input/example_59347_valzone24.tle.txt'
      tle_file(229) = 'input/example_59347_valzone25.tle.txt'
      tle_file(230) = 'input/example_59347_valzone26.tle.txt'
      tle_file(231) = 'input/example_59347_valzone27.tle.txt'
      tle_file(232) = 'input/example_59347_valzone28.tle.txt'
      tle_file(233) = 'input/example_59347_valzone29.tle.txt'
      tle_file(234) = 'input/example_59347_valzone30.tle.txt'
      tle_file(235) = 'input/example_59347_valzone31.tle.txt'
      tle_file(236) = 'input/example_59347_valzone32.tle.txt'
      tle_file(237) = 'input/example_59347_valzone33.tle.txt'
      tle_file(238) = 'input/example_59347_valzone34.tle.txt'
      tle_file(239) = 'input/example_59347_valzone35.tle.txt'
      tle_file(240) = 'input/example_59347_valzone36.tle.txt'
      tle_file(241) = 'input/example_40943_valzone1.tle.txt'
      tle_file(242) = 'input/example_40943_valzone2.tle.txt'
      tle_file(243) = 'input/example_40943_valzone3.tle.txt'
      tle_file(244) = 'input/example_40943_valzone4.tle.txt'
      tle_file(245) = 'input/example_40943_valzone5.tle.txt'
      tle_file(246) = 'input/example_40943_valzone6.tle.txt'
      tle_file(247) = 'input/example_40943_valzone7.tle.txt'
      tle_file(248) = 'input/example_40943_valzone8.tle.txt'
      tle_file(249) = 'input/example_40943_valzone9.tle.txt'
      tle_file(250) = 'input/example_40943_valzone10.tle.txt'
      tle_file(251) = 'input/example_40943_valzone11.tle.txt'
      tle_file(252) = 'input/example_40943_valzone12.tle.txt'
      tle_file(253) = 'input/example_40943_valzone13.tle.txt'
      tle_file(254) = 'input/example_40943_valzone14.tle.txt'
      tle_file(255) = 'input/example_40943_valzone15.tle.txt'
      tle_file(256) = 'input/example_40943_valzone16.tle.txt'
      tle_file(257) = 'input/example_40943_valzone17.tle.txt'
      tle_file(258) = 'input/example_40943_valzone18.tle.txt'
      tle_file(259) = 'input/example_40943_valzone19.tle.txt'
      tle_file(260) = 'input/example_40943_valzone20.tle.txt'
      tle_file(261) = 'input/example_40943_valzone21.tle.txt'
      tle_file(262) = 'input/example_40943_valzone22.tle.txt'
      tle_file(263) = 'input/example_40943_valzone23.tle.txt'
      tle_file(264) = 'input/example_40943_valzone24.tle.txt'
      tle_file(265) = 'input/example_40943_valzone25.tle.txt'
      tle_file(266) = 'input/example_40943_valzone26.tle.txt'
      tle_file(267) = 'input/example_40943_valzone27.tle.txt'
      tle_file(268) = 'input/example_40943_valzone28.tle.txt'
      tle_file(269) = 'input/example_40943_valzone29.tle.txt'
      tle_file(270) = 'input/example_40943_valzone30.tle.txt'
      tle_file(271) = 'input/example_40943_valzone31.tle.txt'
      tle_file(272) = 'input/example_40943_valzone32.tle.txt'
      tle_file(273) = 'input/example_40943_valzone33.tle.txt'
      tle_file(274) = 'input/example_40943_valzone34.tle.txt'
      tle_file(275) = 'input/example_40943_valzone35.tle.txt'
      tle_file(276) = 'input/example_40943_valzone36.tle.txt'
      tle_file(277) = 'input/example_40943_valzone37.tle.txt'
      tle_file(278) = 'input/example_40943_valzone38.tle.txt'
      tle_file(279) = 'input/example_40943_valzone39.tle.txt'
      tle_file(280) = 'input/example_40943_valzone40.tle.txt'
      tle_file(281) = 'input/example_40943_valzone41.tle.txt'
      tle_file(282) = 'input/example_40943_valzone42.tle.txt'
      tle_file(283) = 'input/example_40943_valzone43.tle.txt'
      tle_file(284) = 'input/example_40943_valzone44.tle.txt'
      tle_file(285) = 'input/example_40943_valzone45.tle.txt'
      tle_file(286) = 'input/example_40943_valzone46.tle.txt'
      tle_file(287) = 'input/example_40943_valzone47.tle.txt'
      tle_file(288) = 'input/example_40943_valzone48.tle.txt'
      tle_file(289) = 'input/example_40943_valzone49.tle.txt'
      tle_file(290) = 'input/example_40943_valzone50.tle.txt'
      tle_file(291) = 'input/example_40943_valzone51.tle.txt'
      tle_file(292) = 'input/example_40943_valzone52.tle.txt'
      tle_file(293) = 'input/example_40943_valzone53.tle.txt'
      tle_file(294) = 'input/example_40943_valzone54.tle.txt'
      tle_file(295) = 'input/example_40943_valzone55.tle.txt'
      tle_file(296) = 'input/example_40943_valzone56.tle.txt'
      tle_file(297) = 'input/example_40943_valzone57.tle.txt'
      tle_file(298) = 'input/example_40943_valzone58.tle.txt'
      tle_file(299) = 'input/example_40943_valzone59.tle.txt'
      tle_file(300) = 'input/example_40943_valzone60.tle.txt'
      tle_file(301) = 'input/example_40943_valzone61.tle.txt'
      tle_file(302) = 'input/example_40943_valzone62.tle.txt'
      tle_file(303) = 'input/example_40943_valzone63.tle.txt'
      tle_file(304) = 'input/example_40943_valzone64.tle.txt'
      tle_file(305) = 'input/example_40943_valzone65.tle.txt'
      tle_file(306) = 'input/example_40943_valzone66.tle.txt'
      tle_file(307) = 'input/example_40943_valzone67.tle.txt'
      tle_file(308) = 'input/example_40943_valzone68.tle.txt'
      tle_file(309) = 'input/example_40943_valzone69.tle.txt'
      tle_file(310) = 'input/example_40943_valzone70.tle.txt'
      tle_file(311) = 'input/example_40943_valzone71.tle.txt'
      tle_file(312) = 'input/example_40943_valzone72.tle.txt'
      tle_file(313) = 'input/example_40943_valzone73.tle.txt'
      tle_file(314) = 'input/example_40943_valzone74.tle.txt'
      tle_file(315) = 'input/example_40943_valzone75.tle.txt'
      tle_file(316) = 'input/example_40943_valzone76.tle.txt'
      tle_file(317) = 'input/example_40943_valzone77.tle.txt'
      tle_file(318) = 'input/example_40943_valzone78.tle.txt'
      tle_file(319) = 'input/example_40943_valzone79.tle.txt'
      tle_file(320) = 'input/example_40943_valzone80.tle.txt'
      tle_file(321) = 'input/example_40943_valzone81.tle.txt'
      tle_file(322) = 'input/example_40943_valzone82.tle.txt'
      tle_file(323) = 'input/example_40943_valzone83.tle.txt'
      tle_file(324) = 'input/example_40943_valzone84.tle.txt'
      tle_file(325) = 'input/example_40943_valzone85.tle.txt'
      tle_file(326) = 'input/example_40943_valzone86.tle.txt'
      tle_file(327) = 'input/example_40943_valzone87.tle.txt'
      tle_file(328) = 'input/example_40943_valzone88.tle.txt'
      tle_file(329) = 'input/example_40943_valzone89.tle.txt'
      tle_file(330) = 'input/example_40943_valzone90.tle.txt'
      tle_file(331) = 'input/example_40943_valzone91.tle.txt'
      tle_file(332) = 'input/example_40943_valzone92.tle.txt'
      tle_file(333) = 'input/example_40943_valzone93.tle.txt'
      tle_file(334) = 'input/example_40943_valzone94.tle.txt'
      tle_file(335) = 'input/example_40943_valzone95.tle.txt'
      tle_file(336) = 'input/example_10724_valzone1.tle.txt'
      tle_file(337) = 'input/example_10724_valzone2.tle.txt'
      tle_file(338) = 'input/example_10724_valzone3.tle.txt'
      tle_file(339) = 'input/example_10724_valzone4.tle.txt'
      tle_file(340) = 'input/example_18954_valzone1.tle.txt'
      obj_name(1) = '35497 z1 pre-resonance ('
      obj_name(2) = '35497 z2 pre-resonance ('
      obj_name(3) = '35497 z3 pre-resonance ('
      obj_name(4) = '35497 z4 pre-resonance ('
      obj_name(5) = '35497 z5 pre-resonance ('
      obj_name(6) = '35497 z6 pre-resonance ('
      obj_name(7) = '35497 z7 pre-resonance ('
      obj_name(8) = '35497 z8 pre-resonance ('
      obj_name(9) = '35497 z9 pre-resonance ('
      obj_name(10) = '35497 z10 pre-resonance '
      obj_name(11) = '35497 z11 pre-resonance '
      obj_name(12) = '35497 z12 pre-resonance '
      obj_name(13) = '35497 z13 pre-resonance '
      obj_name(14) = '35497 z14 pre-resonance '
      obj_name(15) = '35497 z15 pre-resonance '
      obj_name(16) = '35497 z16 pre-resonance '
      obj_name(17) = '35497 z17 pre-resonance '
      obj_name(18) = '35497 z18 pre-resonance '
      obj_name(19) = '35497 z19 pre-resonance '
      obj_name(20) = '35497 z20 pre-resonance '
      obj_name(21) = '35497 z21 pre-resonance '
      obj_name(22) = '35497 z22 pre-resonance '
      obj_name(23) = '35497 z23 pre-resonance '
      obj_name(24) = '35497 z24 pre-resonance '
      obj_name(25) = '35497 z25 pre-resonance '
      obj_name(26) = '35497 z26 pre-resonance '
      obj_name(27) = '35497 z27 pre-resonance '
      obj_name(28) = '35497 z28 pre-resonance '
      obj_name(29) = '35497 z29 pre-resonance '
      obj_name(30) = '35497 z30 pre-resonance '
      obj_name(31) = '35497 z31 pre-resonance '
      obj_name(32) = '35497 z32 pre-resonance '
      obj_name(33) = '35497 z33 pre-resonance '
      obj_name(34) = '35497 z34 pre-resonance '
      obj_name(35) = '35497 z35 pre-resonance '
      obj_name(36) = '35497 z36 pre-resonance '
      obj_name(37) = '35497 z37 pre-resonance '
      obj_name(38) = '35497 z38 pre-resonance '
      obj_name(39) = '35497 z39 pre-resonance '
      obj_name(40) = '35497 z40 pre-resonance '
      obj_name(41) = '35497 z41 pre-resonance '
      obj_name(42) = '35497 z42 post-resonance'
      obj_name(43) = '35497 z43 post-resonance'
      obj_name(44) = '35497 z44 post-resonance'
      obj_name(45) = '35497 z45 post-resonance'
      obj_name(46) = '35497 z46 post-resonance'
      obj_name(47) = '35497 z47 post-resonance'
      obj_name(48) = '35497 z48 post-resonance'
      obj_name(49) = '35497 z49 post-resonance'
      obj_name(50) = '35497 z50 post-resonance'
      obj_name(51) = '35497 z51 post-resonance'
      obj_name(52) = '35497 z52 post-resonance'
      obj_name(53) = '35497 z53 post-resonance'
      obj_name(54) = '35497 z54 post-resonance'
      obj_name(55) = '37151 z1 pre-resonance ('
      obj_name(56) = '37151 z2 pre-resonance ('
      obj_name(57) = '37151 z3 pre-resonance ('
      obj_name(58) = '37151 z4 pre-resonance ('
      obj_name(59) = '37151 z5 pre-resonance ('
      obj_name(60) = '37151 z6 pre-resonance ('
      obj_name(61) = '37151 z7 pre-resonance ('
      obj_name(62) = '37151 z8 pre-resonance ('
      obj_name(63) = '37151 z9 pre-resonance ('
      obj_name(64) = '37151 z10 pre-resonance '
      obj_name(65) = '37151 z11 pre-resonance '
      obj_name(66) = '37151 z12 pre-resonance '
      obj_name(67) = '37151 z13 pre-resonance '
      obj_name(68) = '37151 z14 pre-resonance '
      obj_name(69) = '37151 z15 pre-resonance '
      obj_name(70) = '37151 z16 pre-resonance '
      obj_name(71) = '37151 z17 pre-resonance '
      obj_name(72) = '37151 z18 pre-resonance '
      obj_name(73) = '37151 z19 pre-resonance '
      obj_name(74) = '37151 z20 pre-resonance '
      obj_name(75) = '37151 z21 pre-resonance '
      obj_name(76) = '37151 z22 pre-resonance '
      obj_name(77) = '37151 z23 pre-resonance '
      obj_name(78) = '37151 z24 pre-resonance '
      obj_name(79) = '37151 z25 pre-resonance '
      obj_name(80) = '37151 z26 pre-resonance '
      obj_name(81) = '37151 z27 pre-resonance '
      obj_name(82) = '37151 z28 pre-resonance '
      obj_name(83) = '37151 z29 pre-resonance '
      obj_name(84) = '37151 z30 pre-resonance '
      obj_name(85) = '37151 z31 pre-resonance '
      obj_name(86) = '37151 z32 pre-resonance '
      obj_name(87) = '37151 z33 pre-resonance '
      obj_name(88) = '37151 z34 pre-resonance '
      obj_name(89) = '37151 z35 pre-resonance '
      obj_name(90) = '37151 z36 pre-resonance '
      obj_name(91) = '37151 z37 pre-resonance '
      obj_name(92) = '37151 z38 pre-resonance '
      obj_name(93) = '37151 z39 pre-resonance '
      obj_name(94) = '37151 z40 pre-resonance '
      obj_name(95) = '37151 z41 pre-resonance '
      obj_name(96) = '37151 z42 pre-resonance '
      obj_name(97) = '37151 z43 pre-resonance '
      obj_name(98) = '37151 z44 pre-resonance '
      obj_name(99) = '37151 z45 pre-resonance '
      obj_name(100) = '37151 z46 pre-resonance '
      obj_name(101) = '37151 z47 pre-resonance '
      obj_name(102) = '37151 z48 pre-resonance '
      obj_name(103) = '37151 z49 pre-resonance '
      obj_name(104) = '37151 z50 pre-resonance '
      obj_name(105) = '37151 z51 pre-resonance '
      obj_name(106) = '37151 z52 pre-resonance '
      obj_name(107) = '37151 z53 pre-resonance '
      obj_name(108) = '37151 z54 pre-resonance '
      obj_name(109) = '37151 z55 pre-resonance '
      obj_name(110) = '37151 z56 pre-resonance '
      obj_name(111) = '37151 z57 pre-resonance '
      obj_name(112) = '37151 z58 pre-resonance '
      obj_name(113) = '37151 z59 pre-resonance '
      obj_name(114) = '37151 z60 pre-resonance '
      obj_name(115) = '37151 z61 pre-resonance '
      obj_name(116) = '37151 z62 pre-resonance '
      obj_name(117) = '37151 z63 pre-resonance '
      obj_name(118) = '37151 z64 pre-resonance '
      obj_name(119) = '37151 z65 pre-resonance '
      obj_name(120) = '37151 z66 pre-resonance '
      obj_name(121) = '37151 z67 pre-resonance '
      obj_name(122) = '37151 z68 post-resonance'
      obj_name(123) = '37151 z69 post-resonance'
      obj_name(124) = '37151 z70 post-resonance'
      obj_name(125) = '37151 z71 post-resonance'
      obj_name(126) = '37151 z72 post-resonance'
      obj_name(127) = '37151 z73 post-resonance'
      obj_name(128) = '37151 z74 post-resonance'
      obj_name(129) = '37151 z75 post-resonance'
      obj_name(130) = '37151 z76 post-resonance'
      obj_name(131) = '37151 z77 post-resonance'
      obj_name(132) = '37151 z78 post-resonance'
      obj_name(133) = '37151 z79 post-resonance'
      obj_name(134) = '37151 z80 post-resonance'
      obj_name(135) = '37151 z81 post-resonance'
      obj_name(136) = '37151 z82 post-resonance'
      obj_name(137) = '37151 z83 post-resonance'
      obj_name(138) = '37151 z84 post-resonance'
      obj_name(139) = '37151 z85 post-resonance'
      obj_name(140) = '37151 z86 post-resonance'
      obj_name(141) = '37151 z87 post-resonance'
      obj_name(142) = '37151 z88 post-resonance'
      obj_name(143) = '37151 z89 post-resonance'
      obj_name(144) = '37151 z90 post-resonance'
      obj_name(145) = '37151 z91 post-resonance'
      obj_name(146) = '37151 z92 post-resonance'
      obj_name(147) = '37151 z93 post-resonance'
      obj_name(148) = '37151 z94 post-resonance'
      obj_name(149) = '37151 z95 post-resonance'
      obj_name(150) = '37151 z96 post-resonance'
      obj_name(151) = '37151 z97 post-resonance'
      obj_name(152) = '37151 z98 post-resonance'
      obj_name(153) = '37151 z99 post-resonance'
      obj_name(154) = '37151 z100 post-resonanc'
      obj_name(155) = '37151 z101 post-resonanc'
      obj_name(156) = '37151 z102 post-resonanc'
      obj_name(157) = '37151 z103 post-resonanc'
      obj_name(158) = '27526 z1 pre-resonance ('
      obj_name(159) = '27526 z2 pre-resonance ('
      obj_name(160) = '27526 z3 pre-resonance ('
      obj_name(161) = '27526 z4 pre-resonance ('
      obj_name(162) = '27526 z5 pre-resonance ('
      obj_name(163) = '27526 z6 pre-resonance ('
      obj_name(164) = '27526 z7 pre-resonance ('
      obj_name(165) = '27526 z8 pre-resonance ('
      obj_name(166) = '27526 z9 pre-resonance ('
      obj_name(167) = '27526 z10 pre-resonance '
      obj_name(168) = '27526 z11 pre-resonance '
      obj_name(169) = '27526 z12 pre-resonance '
      obj_name(170) = '27526 z13 pre-resonance '
      obj_name(171) = '27526 z14 pre-resonance '
      obj_name(172) = '27526 z15 pre-resonance '
      obj_name(173) = '27526 z16 pre-resonance '
      obj_name(174) = '27526 z17 pre-resonance '
      obj_name(175) = '27526 z18 post-resonance'
      obj_name(176) = '27526 z19 post-resonance'
      obj_name(177) = '27526 z20 post-resonance'
      obj_name(178) = '27526 z21 post-resonance'
      obj_name(179) = '27526 z22 post-resonance'
      obj_name(180) = '27526 z23 post-resonance'
      obj_name(181) = '27526 z24 post-resonance'
      obj_name(182) = '27526 z25 post-resonance'
      obj_name(183) = '27526 z26 post-resonance'
      obj_name(184) = '27526 z27 post-resonance'
      obj_name(185) = '27526 z28 post-resonance'
      obj_name(186) = '27526 z29 post-resonance'
      obj_name(187) = '27526 z30 post-resonance'
      obj_name(188) = '27526 z31 post-resonance'
      obj_name(189) = '27526 z32 post-resonance'
      obj_name(190) = '27526 z33 post-resonance'
      obj_name(191) = '27526 z34 post-resonance'
      obj_name(192) = '27526 z35 post-resonance'
      obj_name(193) = '27526 z36 post-resonance'
      obj_name(194) = '27526 z37 post-resonance'
      obj_name(195) = '27526 z38 post-resonance'
      obj_name(196) = '27526 z39 post-resonance'
      obj_name(197) = '27526 z40 post-resonance'
      obj_name(198) = '27526 z41 post-resonance'
      obj_name(199) = '27526 z42 post-resonance'
      obj_name(200) = '27526 z43 post-resonance'
      obj_name(201) = '27526 z44 post-resonance'
      obj_name(202) = '27526 z45 post-resonance'
      obj_name(203) = '27526 z46 post-resonance'
      obj_name(204) = '27526 z47 post-resonance'
      obj_name(205) = '59347 z1 pre-resonance ('
      obj_name(206) = '59347 z2 pre-resonance ('
      obj_name(207) = '59347 z3 pre-resonance ('
      obj_name(208) = '59347 z4 pre-resonance ('
      obj_name(209) = '59347 z5 pre-resonance ('
      obj_name(210) = '59347 z6 pre-resonance ('
      obj_name(211) = '59347 z7 pre-resonance ('
      obj_name(212) = '59347 z8 pre-resonance ('
      obj_name(213) = '59347 z9 pre-resonance ('
      obj_name(214) = '59347 z10 pre-resonance '
      obj_name(215) = '59347 z11 pre-resonance '
      obj_name(216) = '59347 z12 pre-resonance '
      obj_name(217) = '59347 z13 pre-resonance '
      obj_name(218) = '59347 z14 pre-resonance '
      obj_name(219) = '59347 z15 pre-resonance '
      obj_name(220) = '59347 z16 pre-resonance '
      obj_name(221) = '59347 z17 pre-resonance '
      obj_name(222) = '59347 z18 pre-resonance '
      obj_name(223) = '59347 z19 pre-resonance '
      obj_name(224) = '59347 z20 pre-resonance '
      obj_name(225) = '59347 z21 pre-resonance '
      obj_name(226) = '59347 z22 pre-resonance '
      obj_name(227) = '59347 z23 pre-resonance '
      obj_name(228) = '59347 z24 pre-resonance '
      obj_name(229) = '59347 z25 pre-resonance '
      obj_name(230) = '59347 z26 pre-resonance '
      obj_name(231) = '59347 z27 pre-resonance '
      obj_name(232) = '59347 z28 pre-resonance '
      obj_name(233) = '59347 z29 pre-resonance '
      obj_name(234) = '59347 z30 pre-resonance '
      obj_name(235) = '59347 z31 pre-resonance '
      obj_name(236) = '59347 z32 pre-resonance '
      obj_name(237) = '59347 z33 pre-resonance '
      obj_name(238) = '59347 z34 pre-resonance '
      obj_name(239) = '59347 z35 pre-resonance '
      obj_name(240) = '59347 z36 post-resonance'
      obj_name(241) = '40943 z1 pre-resonance ('
      obj_name(242) = '40943 z2 pre-resonance ('
      obj_name(243) = '40943 z3 pre-resonance ('
      obj_name(244) = '40943 z4 pre-resonance ('
      obj_name(245) = '40943 z5 pre-resonance ('
      obj_name(246) = '40943 z6 pre-resonance ('
      obj_name(247) = '40943 z7 pre-resonance ('
      obj_name(248) = '40943 z8 pre-resonance ('
      obj_name(249) = '40943 z9 pre-resonance ('
      obj_name(250) = '40943 z10 pre-resonance '
      obj_name(251) = '40943 z11 pre-resonance '
      obj_name(252) = '40943 z12 pre-resonance '
      obj_name(253) = '40943 z13 pre-resonance '
      obj_name(254) = '40943 z14 pre-resonance '
      obj_name(255) = '40943 z15 pre-resonance '
      obj_name(256) = '40943 z16 pre-resonance '
      obj_name(257) = '40943 z17 pre-resonance '
      obj_name(258) = '40943 z18 pre-resonance '
      obj_name(259) = '40943 z19 pre-resonance '
      obj_name(260) = '40943 z20 pre-resonance '
      obj_name(261) = '40943 z21 pre-resonance '
      obj_name(262) = '40943 z22 pre-resonance '
      obj_name(263) = '40943 z23 pre-resonance '
      obj_name(264) = '40943 z24 pre-resonance '
      obj_name(265) = '40943 z25 pre-resonance '
      obj_name(266) = '40943 z26 pre-resonance '
      obj_name(267) = '40943 z27 pre-resonance '
      obj_name(268) = '40943 z28 pre-resonance '
      obj_name(269) = '40943 z29 pre-resonance '
      obj_name(270) = '40943 z30 pre-resonance '
      obj_name(271) = '40943 z31 pre-resonance '
      obj_name(272) = '40943 z32 pre-resonance '
      obj_name(273) = '40943 z33 pre-resonance '
      obj_name(274) = '40943 z34 pre-resonance '
      obj_name(275) = '40943 z35 pre-resonance '
      obj_name(276) = '40943 z36 pre-resonance '
      obj_name(277) = '40943 z37 pre-resonance '
      obj_name(278) = '40943 z38 pre-resonance '
      obj_name(279) = '40943 z39 pre-resonance '
      obj_name(280) = '40943 z40 pre-resonance '
      obj_name(281) = '40943 z41 pre-resonance '
      obj_name(282) = '40943 z42 pre-resonance '
      obj_name(283) = '40943 z43 pre-resonance '
      obj_name(284) = '40943 z44 pre-resonance '
      obj_name(285) = '40943 z45 pre-resonance '
      obj_name(286) = '40943 z46 pre-resonance '
      obj_name(287) = '40943 z47 pre-resonance '
      obj_name(288) = '40943 z48 pre-resonance '
      obj_name(289) = '40943 z49 pre-resonance '
      obj_name(290) = '40943 z50 pre-resonance '
      obj_name(291) = '40943 z51 pre-resonance '
      obj_name(292) = '40943 z52 pre-resonance '
      obj_name(293) = '40943 z53 pre-resonance '
      obj_name(294) = '40943 z54 pre-resonance '
      obj_name(295) = '40943 z55 pre-resonance '
      obj_name(296) = '40943 z56 pre-resonance '
      obj_name(297) = '40943 z57 pre-resonance '
      obj_name(298) = '40943 z58 pre-resonance '
      obj_name(299) = '40943 z59 pre-resonance '
      obj_name(300) = '40943 z60 pre-resonance '
      obj_name(301) = '40943 z61 pre-resonance '
      obj_name(302) = '40943 z62 pre-resonance '
      obj_name(303) = '40943 z63 pre-resonance '
      obj_name(304) = '40943 z64 pre-resonance '
      obj_name(305) = '40943 z65 pre-resonance '
      obj_name(306) = '40943 z66 pre-resonance '
      obj_name(307) = '40943 z67 pre-resonance '
      obj_name(308) = '40943 z68 pre-resonance '
      obj_name(309) = '40943 z69 pre-resonance '
      obj_name(310) = '40943 z70 pre-resonance '
      obj_name(311) = '40943 z71 pre-resonance '
      obj_name(312) = '40943 z72 pre-resonance '
      obj_name(313) = '40943 z73 pre-resonance '
      obj_name(314) = '40943 z74 pre-resonance '
      obj_name(315) = '40943 z75 pre-resonance '
      obj_name(316) = '40943 z76 pre-resonance '
      obj_name(317) = '40943 z77 pre-resonance '
      obj_name(318) = '40943 z78 pre-resonance '
      obj_name(319) = '40943 z79 pre-resonance '
      obj_name(320) = '40943 z80 pre-resonance '
      obj_name(321) = '40943 z81 pre-resonance '
      obj_name(322) = '40943 z82 pre-resonance '
      obj_name(323) = '40943 z83 pre-resonance '
      obj_name(324) = '40943 z84 post-resonance'
      obj_name(325) = '40943 z85 post-resonance'
      obj_name(326) = '40943 z86 post-resonance'
      obj_name(327) = '40943 z87 post-resonance'
      obj_name(328) = '40943 z88 post-resonance'
      obj_name(329) = '40943 z89 post-resonance'
      obj_name(330) = '40943 z90 post-resonance'
      obj_name(331) = '40943 z91 post-resonance'
      obj_name(332) = '40943 z92 post-resonance'
      obj_name(333) = '40943 z93 post-resonance'
      obj_name(334) = '40943 z94 post-resonance'
      obj_name(335) = '40943 z95 post-resonance'
      obj_name(336) = '10724 z1 pre-resonance ('
      obj_name(337) = '10724 z2 pre-resonance ('
      obj_name(338) = '10724 z3 post-resonance '
      obj_name(339) = '10724 z4 post-resonance '
      obj_name(340) = '18954 z1 pre-resonance ('

      io_start = 1
      io_end = nobj
      csv_name = 'scratch_rpe/validation_RPE/valrpe_campaign.csv'
      iarg = iargc()
      if (iarg .ge. 2) then
         call getarg(1, arg_str)
         read(arg_str, *) io_start
         call getarg(2, arg_str)
         read(arg_str, *) io_end
         write(csv_name, '(A,I2.2,A,I2.2,A)')
     &      'scratch_rpe/validation_RPE/campaign_part_', io_start,
     &      '_', io_end, '.csv'
      end if
      call init_constants()
      call read_atm_c(ALT_a, DEN_a, SCH_a, ndim_a, natm)

      open(unit=33, file=csv_name, status='unknown')
      write(33,'(A)') 'norad,zone,e_opt,bn_opt,reentry_jd,rpe_pct,'//
     &   'zstat,t_mean,t_std,ens_rpe_pct,zepoch,rms_fit'

      do io = io_start, io_end
         write(*,'(/,A,A)') ' ===== ', obj_name(io)

         t_obs(1) = oyr(io)
         t_obs(2) = omo(io)
         t_obs(3) = ody(io)
         t_obs(4) = 0.d0
         t_obs(5) = 0.d0
         t_obs(6) = 0.d0

         call orem_run(
     &      tle_file(io), norad(io),
     &      t_obs, mxz, 3, 900.d0, 0.d0, 1.d6,
     &      80.d0, 160.d0, 0, 1, 0, 0,
     &      20, 500, 40, 16, 0.8d0, 0.01d0, 0.123d0,
     &      20, 2, 3,
     &      7.2921150d-5, 3.35281066d-3, 1.d0,
     &      1.2d0, 0.01d0, 1, 2,
     &      4.56d-6, 1.32712440018d11, 4.902801076d3,
     &      ALT_a, DEN_a, SCH_a, natm,
     &      reentry, e_opt, bn_opt, rms_o,
     &      zepoch, nzu,
     &      zstat, nzv,
     &      rpe, t_mean, t_std,
     &      bnb_t, tb_t, nbv_t, ierr)

         if (ierr .ne. 0 .and. ierr .ne. 3) then
            write(*,'(A,I3)') '  pipeline error ierr=', ierr
            write(33,'(I5,A,I3)') norad(io), ',ERR,,,,,,,,,,', ierr
            goto 90
         end if

         call cal2jd(t_obs, t_obs_jd)
         rpe_ens = 0.d0
         nre = 0
         do iz = 1, nzu
            if (reentry(iz) .gt. 0.d0) nre = nre + 1
         end do
         if (nre .ge. 1 .and. t_obs_jd .gt. 0.d0) then
            horiz = t_obs_jd - zepoch(1)
            if (horiz .gt. 0.d0)
     &         rpe_ens = (t_mean - t_obs_jd)/horiz*100.d0
         end if

         do iz = 1, nzu
            if (reentry(iz) .gt. 0.d0) then
               write(*,'(A,I1,A,F8.4,A,F7.2,A,F8.2,A,I2)')
     &            '  Z', iz, ': e=', e_opt(iz), ' BN=', bn_opt(iz),
     &            ' RPE=', rpe(iz), '%  zstat=', zstat(iz)
            else
               write(*,'(A,I1,A,F8.4,A,F7.2,A,I2)')
     &            '  Z', iz, ': e=', e_opt(iz), ' BN=', bn_opt(iz),
     &            ' (no re-entry)  zstat=', zstat(iz)
            end if
            write(33,'(I5,A,I1,A,F9.6,A,F8.3,A,F13.4,A,F9.3,A,I2,
     &         A,F13.4,A,F8.2,A,F9.3,A,F13.4,A,F14.4)')
     &         norad(io), ',', iz, ',', e_opt(iz), ',', bn_opt(iz),
     &         ',', reentry(iz), ',', rpe(iz), ',', zstat(iz),
     &         ',', t_mean, ',', t_std, ',', rpe_ens,
     &         ',', zepoch(iz), ',', rms_o(iz)
         end do

         if (nre .ge. 1) then
            write(*,'(A,I1,A,I1,A,F13.3,A,F7.1,A)')
     &         '  ensemble (', nre, ' of ', nzu,
     &         ' zones): mean JD=', t_mean, ' std=', t_std, ' d'
            write(*,'(A,F8.2,A)')
     &         '  RPE = ', rpe_ens, '%'
         else
            write(*,'(A)') '  no re-entry (weak-signal fit, or '//
     &         'zone_select still rejected this window)'
         end if

 90      continue
         write(*,'(A,I3,A,I2,A,I2,A)') ' [PROGRESS] ',
     &      ((io-io_start+1)*100)/(io_end-io_start+1), '% (',
     &      io, ' of ', io_end, ' zones done, this chunk)'
      end do

      close(33)
      write(*,'(/,A,A)') ' campaign complete -> ', trim(csv_name)

      stop
      end

c     ============================================================
c     ATM reader (same as rpe_campaign.F / rpe_campaign_resonance16.F)
c     ============================================================
      subroutine read_atm_c(alt, den, sch, maxn, natm)
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

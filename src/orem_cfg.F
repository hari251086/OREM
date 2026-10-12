c     ============================================================
c     orem_cfg.F -- reader for the orem.cfg format shared by
c     app/main_orem.F and app/object_detail_data.F (and written by
c     OREM-Watchlist's orem_wrapper.write_orem_cfg).
c
c     Each line is read on its own and parsed from that line only,
c     so a line with too few values is an error. Until 2026-10-12
c     both apps read the file list-directed straight from the unit,
c     where a short line silently takes its missing values from the
c     NEXT line: every config shipped in input/ was one value short on
c     line 5 after primary_max_days was added (2026-09-27), which
c     shifted every later parameter (slope threshold 80, population 4,
c     2 generations) without any message.
c
c     ierr: 0 = ok, 1 = cannot open, 2 = line n_bad missing or short
c     ============================================================
      subroutine orem_read_cfg(cfg_file,
     &   tle_file, norad_id, t_obs_cal, nzones_max,
     &   min_zone_pts, max_zone_days, primary_max_days,
     &   r2_thresh, slope_thresh,
     &   bn_min, bn_max, n_bn_sweep,
     &   ipopsize, maxgen, nbits_e, nbits_a, pcross, pmute, ga_seed,
     &   ngeo_deg, nsun_deg, nmoon_deg,
     &   idrag_flag, WE_rot, EPS_f, FR_rot,
     &   IPSR, CR_srp, AM_srp, ISHAD,
     &   ierr, n_bad)
      implicit double precision (a-h, o-z)
      character*(*) cfg_file, tle_file
      integer norad_id, nzones_max, min_zone_pts, n_bn_sweep
      integer ipopsize, maxgen, nbits_e, nbits_a
      integer ngeo_deg, nsun_deg, nmoon_deg, idrag_flag
      integer IPSR, ISHAD, ierr, n_bad
      double precision t_obs_cal(6)
c     explicit: max_zone_days would otherwise be INTEGER (i-n rule)
      double precision max_zone_days, primary_max_days
      double precision r2_thresh, slope_thresh, bn_min, bn_max
      double precision pcross, pmute, ga_seed
      double precision WE_rot, EPS_f, FR_rot, CR_srp, AM_srp

      character*400 cline

      ierr = 0
      n_bad = 0
      open(unit=10, file=cfg_file, status='old', err=900)

      n_bad = 1
      read(10,'(A)',err=910,end=910) tle_file
      n_bad = 2
      read(10,'(A)',err=910,end=910) cline
      read(cline,*,err=910,end=910) norad_id
      n_bad = 3
      read(10,'(A)',err=910,end=910) cline
      read(cline,*,err=910,end=910) (t_obs_cal(i), i = 1, 6)
      n_bad = 4
      read(10,'(A)',err=910,end=910) cline
      read(cline,*,err=910,end=910) nzones_max
      n_bad = 5
      read(10,'(A)',err=910,end=910) cline
      read(cline,*,err=910,end=910) min_zone_pts, max_zone_days,
     &   primary_max_days, r2_thresh, slope_thresh
      n_bad = 6
      read(10,'(A)',err=910,end=910) cline
      read(cline,*,err=910,end=910) bn_min, bn_max, n_bn_sweep
      n_bad = 7
      read(10,'(A)',err=910,end=910) cline
      read(cline,*,err=910,end=910) ipopsize, maxgen, nbits_e,
     &   nbits_a, pcross, pmute, ga_seed
      n_bad = 8
      read(10,'(A)',err=910,end=910) cline
      read(cline,*,err=910,end=910) ngeo_deg, nsun_deg, nmoon_deg
      n_bad = 9
      read(10,'(A)',err=910,end=910) cline
      read(cline,*,err=910,end=910) idrag_flag, WE_rot, EPS_f, FR_rot
      n_bad = 10
      read(10,'(A)',err=910,end=910) cline
      read(cline,*,err=910,end=910) IPSR, CR_srp, AM_srp, ISHAD
      close(10)
      n_bad = 0
      return

 900  ierr = 1
      return
 910  ierr = 2
      close(10)
      return
      end

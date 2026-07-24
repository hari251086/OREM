c     ============================================================
c     swx.F — Space weather + 2-D atmosphere loaders (issue #26)
c
c     Fills two common blocks consumed by propagate_ks's per-
c     revolution drag setup. Neither loader is required: when
c     nsw=0 or n2d_t=0 (the compile-time default), propagate_ks
c     uses the legacy single-table ALT/DEN/SCH arguments unchanged,
c     so every existing call site and test is unaffected until a
c     caller explicitly loads both files.
c
c       call sw_load('input/SW-All.csv', ierr)      ! CelesTrak CSV
c       call atm2d_load('input/ATM2D.DAT', ierr)    ! gen_atm2d_jr71
c
c     /swdat/: per-row (JD, F10.7_ADJ, F10.7_ADJ_CENTER81, mean Kp)
c       from the CelesTrak file. Observed rows are daily; predicted
c       rows (PRM/PRD types) are monthly — lookup is by JD binary
c       search (latest row with jd <= t), never by index arithmetic.
c       Missing Kp fields (predicted rows) default to Kp=2.0.
c     /atm2dc/: rho(alt,Tinf) and H(alt,Tinf) tables, stored in the
c       SAME internal scaled units as the legacy path (alt in m,
c       sch in m, den x 1e10) so the drag formulas are untouched.
c
c     Author: Harishkumar Sellamuthu
c     ============================================================
      subroutine sw_load(fname, ierr)
      implicit none
      character*(*) fname
      integer ierr

      integer NSWMAX
      parameter (NSWMAX = 31000)
      double precision sw_jd(NSWMAX), sw_f107(NSWMAX)
      double precision sw_f81(NSWMAX), sw_kp(NSWMAX)
      integer nsw
      common /swdat/ sw_jd, sw_f107, sw_f81, sw_kp, nsw

      character*400 line
      character*40 tok(31)
      double precision cal(6), jd, f107, f81, akp, ksum
      integer i, j, ntok, ip, lt, ios

      ierr = 0
      nsw = 0
      open(unit=91, file=fname, status='old', err=900)
      read(91,'(A)',iostat=ios) line          ! header
      if (ios .ne. 0) goto 900

 10   read(91,'(A)',end=50,err=50) line
      if (len_trim(line) .lt. 10) goto 10

c     --- split on commas into up to 31 tokens ---
      ntok = 1
      tok(1) = ' '
      ip = 0
      do i = 1, len_trim(line)
         if (line(i:i) .eq. ',') then
            if (ntok .ge. 31) goto 12
            ntok = ntok + 1
            tok(ntok) = ' '
            ip = 0
         else
            ip = ip + 1
            tok(ntok)(ip:ip) = line(i:i)
         end if
      end do
 12   continue
      if (ntok .lt. 30) goto 10

c     --- DATE yyyy-mm-dd -> JD (0h UTC) ---
      read(tok(1)(1:4),'(F4.0)',err=10) cal(1)
      read(tok(1)(6:7),'(F2.0)',err=10) cal(2)
      read(tok(1)(9:10),'(F2.0)',err=10) cal(3)
      cal(4) = 0.d0
      cal(5) = 0.d0
      cal(6) = 0.d0
      call cal2jd(cal, jd)

c     --- KP_SUM (token 12, x10 sum of 8) -> daily mean Kp ---
      akp = 2.0d0
      lt = len_trim(tok(12))
      if (lt .gt. 0) then
         read(tok(12)(1:lt),*,err=14) ksum
         akp = ksum / 80.d0
      end if
 14   continue

c     --- F10.7_ADJ (26) and F10.7_ADJ_CENTER81 (30) ---
      lt = len_trim(tok(26))
      if (lt .le. 0) goto 10
      read(tok(26)(1:lt),*,err=10) f107
      lt = len_trim(tok(30))
      if (lt .le. 0) goto 10
      read(tok(30)(1:lt),*,err=10) f81

      if (nsw .lt. NSWMAX) then
         nsw = nsw + 1
         sw_jd(nsw)   = jd
         sw_f107(nsw) = f107
         sw_f81(nsw)  = f81
         sw_kp(nsw)   = akp
      end if
      goto 10

 50   close(91)
      if (nsw .lt. 100) then
         nsw = 0
         ierr = 2
      end if
      return

 900  ierr = 1
      nsw = 0
      return
      end

c     ============================================================
c     NOTE: the runtime consumers of these commons — sw_tinf()
c     (JD -> exospheric temperature) and atm2d_interp() (bilinear
c     rho/H at (h, Tinf)) — live in ksrop/propagate_ks.F so that
c     executables which never load space weather link without this
c     file and run the legacy single-table path untouched.
c     ============================================================

c     ============================================================
c     atm2d_load — read ATM2D.DAT (gen_atm2d_jr71.F) and store in
c     the legacy path's internal scaled units.
c     ============================================================
      subroutine atm2d_load(fname, ierr)
      implicit none
      character*(*) fname
      integer ierr

      integer NA2D, NT2D
      parameter (NA2D = 291, NT2D = 39)
      double precision t2d(NT2D), alt2d(NA2D)
      double precision den2d(NA2D,NT2D), sch2d(NA2D,NT2D)
      integer n2d_alt, n2d_t
      common /atm2dc/ t2d, alt2d, den2d, sch2d, n2d_alt, n2d_t

      character*10 tag
      integer i, j, na, nt

      ierr = 0
      n2d_t = 0
      open(unit=92, file=fname, status='old', err=900)
      read(92,'(A10,2I6)',err=910) tag, na, nt
      if (tag .ne. 'ATM2D_JR71' .or. na .ne. NA2D .or.
     &    nt .ne. NT2D) goto 910
      read(92,*,err=910) (t2d(j), j=1,nt)
      do j = 1, nt
         read(92,*,err=910) (sch2d(i,j), i=1,na)
         read(92,*,err=910) (den2d(i,j), i=1,na)
      end do
      close(92)

      alt2d(1) = 60.d0
      do i = 2, 141
         alt2d(i) = alt2d(i-1) + 1.d0
      end do
      do i = 142, NA2D
         alt2d(i) = alt2d(i-1) + 2.d0
      end do

c     scale to the legacy internal units (m, m, x1e10)
      do i = 1, NA2D
         alt2d(i) = alt2d(i) * 1.d3
         do j = 1, nt
            sch2d(i,j) = sch2d(i,j) * 1.d3
            den2d(i,j) = den2d(i,j) * 1.d10
         end do
      end do

      n2d_alt = na
      n2d_t   = nt
      return

 900  ierr = 1
      return
 910  ierr = 2
      n2d_t = 0
      close(92)
      return
      end


c     COPYRIGHT 2018, HARISHKUMAR SELLAMUTHU, ALL RIGHTS RESERVED
C     -------------------------------------------------------------
c     propagate_ks — KS orbital propagation subroutine
c       Description:
c     Orbit propagation using Kustaanheimo-Stiefel regular elements
c     with numerical integrator RKG4. (Earth)
c     Perturbations included are:
c     - Earth's oblateness ( up to 2190th degree )
c     - Luni-solar perturbations
c     - Atmospheric drag ( Oblate, co-rotating exponential
c       atmosphere referenced to perigee conditions, tabulated
c       density model via caller-provided tables )
c     - Solar radiation pressure (cannon-ball, KSJLSDNP2)
c     -------------------------------------------------------------
c     Author: Harishkumar Sellamuthu
c     hari251086@gmail.com
c     -------------------------------------------------------------
c     Refactored from driver_KS.F (standalone main program) into a
c     callable subroutine. All file I/O removed; data passed via
c     arguments. Caller must call init_constants() before this
c     subroutine to populate common block /xy/.
c     -------------------------------------------------------------
c     Subroutines used
c     geo_coeff,car2ks,ks2car,aLegP,car2oe,ks2ksr,INTPOL
c     -------------------------------------------------------------c

      subroutine propagate_ks(
     &   x0, xd0, cal0,
     &   nrev, istep, tole,
     &   n_force, ngeo_deg, nsun_deg, nmoon_deg,
     &   BN, IDRAG, WE_rot, EPS_f, FR_rot,
     &   CR_srp, AM_srp, IPSR, ISHAD,
     &   PSR_srp, amuS, amuM,
     &   ALT_atm, DEN_atm, SCH_atm, ndim_atm,
     &   max_pts, idump,
     &   traj_jd, traj_x, traj_xd,
     &   exit_code)

      implicit double precision (a-h, o-z)
      parameter(n = 10, cn0 = 2190)
c     Value of pi, degree to radians and radians to degree
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth

c     --- Input arguments ---
      double precision x0(3), xd0(3), cal0(6)
      integer nrev, istep
      double precision tole
      integer n_force(3)
      integer ngeo_deg, nsun_deg, nmoon_deg
      double precision BN
      integer IDRAG
      double precision WE_rot, EPS_f, FR_rot
      double precision CR_srp, AM_srp
      integer IPSR, ISHAD
      double precision PSR_srp, amuS, amuM
      integer ndim_atm
      double precision ALT_atm(ndim_atm)
      double precision DEN_atm(ndim_atm)
      double precision SCH_atm(ndim_atm)
      integer max_pts
c     --- Output arguments ---
      integer idump, exit_code
      double precision traj_jd(max_pts)
      double precision traj_x(3,max_pts)
      double precision traj_xd(3,max_pts)

c     --- Local variables ---
c     Large arrays use SAVE to avoid stack overflow
      dimension pek(7),x(3),xd(3),q(4),u(4),us(4),y(10),z(10)
      dimension qu(3),q_til(4),u_til(4),uu(4)
      dimension ZbyR(cn0*2),ObyR(cn0*2),Re(cn0)
      dimension cn_j(cn0),c_j(cn0),V_p(cn0),aleg(cn0)
      dimension sleg(cn0),oleg(cn0),R(n*2)
      dimension cal(6),ts(3),tm(3),opt(3)
      dimension RbyRs(cn0),RbyRm(cn0),cn_m(cn0),cn_s(cn0)
      dimension Rm(n*2)
      dimension qsun(4),qmoon(4),slambda(cn0)
      dimension olambda(cn0),Rs(n*2)
      save ZbyR, ObyR, Re, cn_j, c_j, V_p, aleg
      save sleg, oleg, RbyRs, RbyRm, cn_m, cn_s
      save slambda, olambda
c     Drag: KS drag perturbation vectors
      dimension qdrag(4),p_drag(3)
c     Drag model state (oblate, co-rotating exponential atmosphere
c     referenced to perigee conditions)
      double precision A_dg,EC_dg,AI_dg,SOMGA_dg,EA_dg
      double precision SI_dg,CSI_dg,SIOMGA_dg,VIPP_dg,DE_dg
      double precision RREQ_dg,SQ_dg,HSQ_dg
      double precision RHOPI_dg,H_dg,BETA_dg,RPO_dg,VPO_dg
      double precision F_dg,XX_dg,XXX_dg,C_dg
      double precision SEA_dg,CEA_dg,TER_dg,TER2_dg
      double precision CSEF_dg,SIEF_dg,FA_dg,exp_arg
c     SRP: cannonball SRP perturbation
      dimension p_srp(3), qsrp(4), r_srp(3)
      double precision PSR_loc

c     --- Initialize exit code ---
      exit_code = 0
      idump = 0

c     --- Copy inputs to local working arrays ---
      do i = 1,3
         x(i)  = x0(i)
         xd(i) = xd0(i)
      end do
      do i = 1,6
         cal(i) = cal0(i)
      end do

c     --- Local copy of PSR so we can zero it if IPSR=0 ---
      PSR_loc = PSR_srp

c     Local copies of constants that force_models may zero
      amuS_loc = amuS
      amuM_loc = amuM
      ngeo_loc = ngeo_deg

      call car2oe(x,xd,pek)
      call oe2car(pek,x,xd,tole)
      call cal2jd(cal,djulian)
      call force_models(n_force,ngeo_loc,amuS_loc,amuM_loc)
      if (IPSR .eq. 0) PSR_loc = 0.d0
      SRPC_srp = CR_srp * AM_srp

      if(ngeo_deg.eq.0.or.ngeo_deg.eq.1) then
         do i = 1,cn0
            c_j(i) = 0.d0
         end do
      else
!     compute j_n values
         call geo_coeff(ngeo_deg,c_j)
      end if

!     Initialization
      E_0 = 0.d0                !Gen. Eccentric Anomaly in radians
      E = E_0
      T = 0.d0                  !Physical time in seconds
      nt = 0                    !Counter in the integrator

      dE_0 = 2.d0*pi/istep      !Step-size in radians

      num_steps = nrev*istep

      call solarnpv(djulian,ts)
      call lunarpv(djulian,tm)

      R(1) = vmn(x)
      V = vmn(xd)

      R(2) = R(1)*R(1)
      V2 = V*V

      do i = 2, 20
         R(i) = R(i-1) * R(1)
      end do

!     Perturbing potential, V_pot

!     Oblateness terms
      Re(1) = R_Earth
      ZbyR(1) = x(3)/R(1)
      ObyR(1) = 1.d0/R(1)

!     Third-body terms
      call third_body_aux(nsun_deg,  x, ts, R(1), amuS_loc,
     &                    Rs(1), RbyRs, first_s, cphi_s, cphi_s2)
      call third_body_aux(nmoon_deg, x, tm, R(1), amuM_loc,
     &                    Rm(1), RbyRm, first_m, cphi_m, cphi_m2)

      do i = 1,ngeo_deg
         ObyR(i+1) = ObyR(i)*ObyR(1)
         Re(i+1) = Re(i)*Re(1)
      end do

      V_pot = 0.d0

      call aLegP(ngeo_deg,ZbyR(1),aleg)
      call aLegP(nsun_deg,cphi_s,sleg)
      call aLegP(nmoon_deg,cphi_m,oleg)

      do i = 2,ngeo_deg
         cn_j(i) = amue*Re(i)*c_j(i)*ObyR(i+1)
         V_p(i) = cn_j(i)*aleg(i+1)
         V_pot = V_pot+V_p(i)
      end do

      do i = 2,nsun_deg
         cn_s(i) = first_s * RbyRs(i)
         V_p(i) = cn_s(i) * sleg(i+1)
         V_pot = V_pot + V_p(i)
      end do

      do i = 2,nmoon_deg
         cn_m(i) = first_m * RbyRm(i)
         V_p(i) = cn_m(i) * oleg(i+1)
         V_pot = V_pot + V_p(i)
      end do

!     Angular frequency
      w = dsqrt(0.5d0*(amue/R(1)-V2/2.d0-V_pot))
      w_kep = dsqrt(0.5d0*(amue/R(1)-V2/2.d0))
      Gam = w/w_kep

      dE = dE_0*Gam

      call car2ks(x,xd,u,us,w)
      call ks2car(u,us,x,xd,w)

      Tow = T+dotp4(u,us)/w     !Time-element

!     Initial KS regular elements
      y(1) = Tow
      y(10) = w

      do i = 1,4
         y(i+1) = u(i)
         y(i+5) = 2.d0*us(i)
      end do

      sie = dsin(E/2.d0)
      cse = dcos(E/2.d0)
      iE = dint(E)

      call car2oe(x,xd,pek)

c     --- Buffer initial state point ---
      idump = idump + 1
      traj_jd(idump) = djulian + T/86400.d0
      do i = 1,3
         traj_x(i,idump)  = x(i)
         traj_xd(i,idump) = xd(i)
      end do

c     Integration
      do ik = 1,nrev            !Revolutions

c     -----------------------------------------------------------
c     Drag: per-revolution reference setup
c     Oblate, co-rotating exponential atmosphere referenced to
c     perigee conditions (same physical model as KSJLSDNP2.F).
c     -----------------------------------------------------------
         call car2oe(x,xd,pek)
         A_dg     = pek(1)*1.d3       ! semi-major axis: km -> m
         EC_dg    = pek(2)
         AI_dg    = pek(3)*d2r        ! inclination: deg -> rad
         SOMGA_dg = pek(5)*d2r        ! argument of perigee: deg -> rad
         EA_dg    = pek(7)*d2r        ! eccentric anomaly: deg -> rad

c     Eccentric-anomaly step that advances EA_dg to the next
c     perigee/apogee crossing over the revolution's istep steps
         VIPP_dg = 4.d0
         if (EA_dg .lt. pi) VIPP_dg = 2.d0
         DE_dg = (VIPP_dg*pi - EA_dg)/istep

         SI_dg     = dsin(AI_dg)
         CSI_dg    = dcos(AI_dg)
         SIOMGA_dg = dsin(SOMGA_dg)

c     Perigee altitude above the oblate Earth (drag is active for
c     the revolution only when this is below 500 km)
         RREQ_dg = R_Earth*1.d3*(1.d0
     &        - SI_dg*SI_dg*SIOMGA_dg*SIOMGA_dg
     &        *EPS_f)
         SQ_dg  = A_dg*(1.d0 - EC_dg)
         HSQ_dg = SQ_dg - RREQ_dg

         if (HSQ_dg .ge. ALT_atm(1) .and.
     &       HSQ_dg .le. ALT_atm(ndim_atm)) then
c     Reference density/scale-height at perigee
            call INTPOL(ALT_atm, DEN_atm, ndim_atm,
     &                  HSQ_dg, RHOPI_dg)
            call INTPOL(ALT_atm, SCH_atm, ndim_atm,
     &                  HSQ_dg, H_dg)
            if (H_dg .le. 0.d0) H_dg = 1.d0
            BETA_dg = 1.d0/H_dg
            RPO_dg  = SQ_dg
            VPO_dg  = dsqrt(amue*1.d9
     &           *(2.d0/RPO_dg - 1.d0/A_dg))
            F_dg    = (1.d0 - RPO_dg*WE_rot
     &           *FR_rot*CSI_dg/VPO_dg)**2
            XX_dg   = A_dg*EC_dg
            XXX_dg  = BETA_dg*XX_dg
            C_dg    = EPS_f*BETA_dg*RPO_dg
     &           *SI_dg*SI_dg*0.5d0
         end if

         do ki = 1,istep        !Steps
 60         continue
            call u2uu(u,uu)

c     -----------------------------------------------------------
c     Atmospheric drag perturbation
c     -----------------------------------------------------------
            if (IDRAG .eq. 1 .and.
     &          HSQ_dg .ge. ALT_atm(1) .and.
     &          HSQ_dg .le. ALT_atm(ndim_atm)) then
               SEA_dg  = dsin(EA_dg)
               CEA_dg  = dcos(EA_dg)
               TER_dg  = 1.d0 - EC_dg*CEA_dg
               CSEF_dg = (CEA_dg - EC_dg)/TER_dg
               SIEF_dg = SEA_dg
     &              *dsqrt(1.d0 - EC_dg*EC_dg)/TER_dg
               FA_dg   = atan3(SIEF_dg,CSEF_dg)
               TER2_dg = dcos(SOMGA_dg)**2
     &                 - dcos(SOMGA_dg+FA_dg)**2
               exp_arg = -2.d0*C_dg*TER2_dg
     &                  - XXX_dg*(1.d0-CEA_dg)
               if (exp_arg .gt. 500.d0) then
                  TER_dg = 0.d0
               else if (exp_arg .lt. -500.d0) then
                  TER_dg = 0.d0
               else
                  TER_dg = dexp(exp_arg)
               end if
               DENS = -F_dg*RHOPI_dg*TER_dg
     &              *V/(BN*2.d0)*1.d-7
               do i = 1,3
                  p_drag(i) = DENS * xd(i)
               end do
            else
               do i = 1,3
                  p_drag(i) = 0.d0
               end do
            end if

c     L(u) transformation: qdrag = L(u) * p_drag
            qdrag(1) =  u(1)*p_drag(1)+u(2)*p_drag(2)
     &           +u(3)*p_drag(3)
            qdrag(2) = -u(2)*p_drag(1)+u(1)*p_drag(2)
     &           +u(4)*p_drag(3)
            qdrag(3) = -u(3)*p_drag(1)-u(4)*p_drag(2)
     &           +u(1)*p_drag(3)
            qdrag(4) =  u(4)*p_drag(1)-u(3)*p_drag(2)
     &           +u(2)*p_drag(3)

c     -----------------------------------------------------------
c     SRP: cannonball solar radiation pressure perturbation.
c     -----------------------------------------------------------
            if (PSR_loc .eq. 0.d0) then
               do i = 1,4
                  qsrp(i) = 0.d0
               end do
            else
               if (ISHAD .eq. 1) then
                  ROHP = shadfncyl(x(1),x(2),x(3),
     &                 ts(1),ts(2),ts(3))
               else if (ISHAD .eq. 2) then
                  ROHP = shadfncone(x(1),x(2),x(3),
     &                 ts(1),ts(2),ts(3))
               else
                  ROHP = 1.d0
               end if
               r_srp(1) = ts(1) - x(1)
               r_srp(2) = ts(2) - x(2)
               r_srp(3) = ts(3) - x(3)
               xsrr = vmn(r_srp)
               do i = 1,3
                  p_srp(i) = -ROHP*PSR_loc
     &                 *SRPC_srp*r_srp(i)/xsrr*1.d-3
               end do
               qsrp(1) =  u(1)*p_srp(1)
     &              +u(2)*p_srp(2)+u(3)*p_srp(3)
               qsrp(2) = -u(2)*p_srp(1)
     &              +u(1)*p_srp(2)+u(4)*p_srp(3)
               qsrp(3) = -u(3)*p_srp(1)
     &              -u(4)*p_srp(2)+u(1)*p_srp(3)
               qsrp(4) =  u(4)*p_srp(1)
     &              -u(3)*p_srp(2)+u(2)*p_srp(3)
            end if

            Tau_geo = 0.d0
            do i = 2,ngeo_deg
               Tau_geo = Tau_geo
     &              -(i-1)*c_j(i)*Re(i)*aleg(i+1)
     &              *ObyR(i)
            end do

            Tau_3body1 = 0.d0
            do i = 2,nsun_deg
               Tau_3body1=Tau_3body1
     &              +(i+2)*sleg(i+1)*RbyRs(i)
     &              *RbyRs(1)*amuS_loc
            end do
            Tau_3body2 = 0.d0
            do i = 2,nmoon_deg
               Tau_3body2=Tau_3body2
     &              +(i+2)*oleg(i+1)*RbyRm(i)
     &              *RbyRm(1)*amuM_loc
            end do

c     Energy element: non-conservative drag + SRP
            z(10) = -0.5d0*(us(1)*(qdrag(1)+qsrp(1))
     &                    + us(2)*(qdrag(2)+qsrp(2))
     &                    + us(3)*(qdrag(3)+qsrp(3))
     &                    + us(4)*(qdrag(4)+qsrp(4)))
     &                    /w

c     Time element
            z(1) = (amue + Tau_geo
     &           + Tau_3body1 + Tau_3body2
     &           + R(1)*(dotp4(u,qdrag)
     &           +dotp4(u,qsrp)))/(8.d0*w*w*w)
     &           - z(10)*(2.d0/(w*w))
     &           *dotp4(u,us)

            if (ngeo_deg  .gt. 1)
     &         call aLegP(ngeo_deg, ZbyR(1), aleg)
            if (nsun_deg  .gt. 1)
     &         call aLegP(nsun_deg, cphi_s,  sleg)
            if (nmoon_deg .gt. 1)
     &         call aLegP(nmoon_deg,cphi_m,  oleg)

            do j = 1,4
               q(j) = 0.d0
               qsun(j) = 0.d0
               qmoon(j) = 0.d0
            end do

            do j = 1,4

!     oblateness terms
               z2mR2=x(3)**2-R(1)**2
               qj=0.d0
               do i = 2,ngeo_deg
                  qj = dble(ObyR(i+1)*Re(i)*c_j(i)*
     &              (R(1)*(i+1)*aleg(i+2)
     &              *(R(1)*uu(j)-x(3)*u(j))+
     &              (aleg(i+1)*(u(j)*(i*R(2)
     &              +x(3)**2)
     &              -uu(j)*(i+1)*R(1)*x(3))))
     &              /z2mR2)
                  q(j)=qj+q(j)
               end do

!     third-body terms
               call u2qu(u,qu,j)

               den_s = cphi_s2 - 1.d0
               den_m = cphi_m2 - 1.d0

               sq = ts(1)*qu(1)+ts(2)*qu(2)
     &              +ts(3)*qu(3)

               do i = 2,nsun_deg
                  slambda(i)=2.d0*R(i-1)*(i+1)
     &                 *(sleg(i+2)
     &                 -cphi_s*sleg(i+1))/den_s

                  qs = (RbyRs(i)*(sq*slambda(i)
     &                 /Rs(1)+
     &                 u(j)*(1-i)*sleg(i+1)
     &                 -slambda(i)*cphi_s))
                  qsun(j) = qs + qsun(j)
               end do

               oq = tm(1)*qu(1)+tm(2)*qu(2)
     &              +tm(3)*qu(3)

               do i = 2,nmoon_deg
                  olambda(i)=2.d0*R(i-1)*(i+1)
     &                 *(oleg(i+2)
     &                 -cphi_m*oleg(i+1))/den_m
                  qm = (RbyRm(i)*(oq*olambda(i)
     &                 /Rm(1)+
     &                 u(j)*(1-i)*oleg(i+1)
     &                 -olambda(i)*cphi_m))
                  qmoon(j) = qm + qmoon(j)
               end do

c     State: conservative + drag + SRP, energy corr
               aqt = sie * 0.25d0
               bqt = cse * 0.25d0

               z(j+1) = (q(j)*amue
     &              + first_s*qsun(j)
     &              + first_m*qmoon(j)
     &              - R(1)*(qdrag(j)+qsrp(j)))
     &              * aqt/(w*w)
     &              + 2.d0*z(10)*us(j)*sie/w

               z(j+5) = -(q(j)*amue
     &              + first_s*qsun(j)
     &              + first_m*qmoon(j)
     &              - R(1)*(qdrag(j)+qsrp(j)))
     &              * bqt/(w*w)
     &              - 2.d0*z(10)*us(j)*cse/w

            end do

            I_rkg = rkgil(n,y,z,E,dE,nt)
            call ks2ksr(y,u,us,E,cse,sie)
            call ks2car(u,us,x,xd,w)
            call car2oe(x,xd,pek)

            w = y(10)
            T = y(1) - dotp4(u,us)/w

            R(1) = dotp4(u,u)

            do i = 2, 20
               R(i) = R(i-1) * R(1)
            end do

            V = vmn(xd)
            ZbyR(1) = x(3)/R(1)
            ObyR(1) = 1.d0/R(1)
            do i = 1,ngeo_deg
               ZbyR(i+1) = ZbyR(i)*ZbyR(1)
               ObyR(i+1) = ObyR(i)*ObyR(1)
            end do
            FracJd = T/86400.d0
            djulian1 = djulian + FracJd
            call solarnpv(djulian1, ts)
            call lunarpv(djulian1, tm)

            call third_body_aux(nsun_deg,
     &           x, ts, R(1), amuS_loc,
     &           Rs(1), RbyRs, first_s,
     &           cphi_s, cphi_s2)
            call third_body_aux(nmoon_deg,
     &           x, tm, R(1), amuM_loc,
     &           Rm(1), RbyRm, first_m,
     &           cphi_m, cphi_m2)

            if (I_rkg .le. 1) goto 60
c     Refresh drag-model eccentricity and advance EA
            EC_dg = pek(2)
            EA_dg = EA_dg + DE_dg

            E_0 = E_0+dE_0
            iE = dint(E_0*r2d)

!     Buffer state at every completed integration step
            idump = idump + 1
            traj_jd(idump) = djulian + T/86400.d0
            do i = 1,3
               traj_x(i,idump)  = x(i)
               traj_xd(i,idump) = xd(i)
            end do

c     -----------------------------------------------------------
c     Re-entry / divergence check
c     -----------------------------------------------------------
            h_alt = R(1) - R_Earth
            if (h_alt .ne. h_alt) then
               idump = idump - 1
               exit_code = 2
               return
            else if (h_alt .lt. 80.d0) then
               exit_code = 1
               return
            end if

         end do

      end do

c     Normal completion
      exit_code = 0
      return
      end

!!$subroutine nsa_initial_conditions(outvar)
!!$  !-----------------------------------------------------------------------
!!$  !****f* Nastal/nsa_inimet
!!$  ! NAME 
!!$  !    nsa_inimet
!!$  ! DESCRIPTION
!!$  !    This routine sets up the initial condition for the meteo benchmarks
!!$  ! USED BY
!!$  !    nsa_iniunk
!!$  !***
!!$  !-----------------------------------------------------------------------
!!$  use      def_master
!!$  use      def_domain
!!$  use      mod_postpr
!!$  use      def_nastal
!!$  implicit none
!!$
!!$  integer(ip) :: idime,icomp,idofn,ipoin,itime,ibubb,nbubb,kshbu,ix,iy,iz
!!$  real(rp)    :: velmi,xfact,xfac2,rgacp,xinve,xradi,xdifi,xrano(3),xtemp,rauxi,exner
!!$  real(rp)    :: u0,v0,w0,u_k,v_k,w_k,thetac,theta0,press0,rho0,dtheta,dtracer,tracerc,tracer0,tracer_k,p00
!!$  real(rp)    :: dqv,dqc,dqr
!!$  real(rp)    :: rc,rc0,rc1,rc2,rc3,xc,xc2,yc,zc,r,xr,yr,zr,zl,xr_tr,yr_tr,zr_tr,sigma
!!$  real(rp)    :: x,y,z, xmin,xmax,ymin,ymax,zmin,zmax
!!$  real(rp)    :: theta_k,pi_k,rho_k,press_k, theta_ref,theta_z0,dtdz
!!$  real(rp)    :: theta_hyd,pi_hyd,rho_hyd,press_hyd, rho_baro
!!$  real(rp)    :: a,b,c,d,e,l,lz,hc,ac,l2,n,n2,c2
!!$  real(rp)    :: rho,theta,pressure, qv,bv,bv2,g2,static_stability,es,gam,pii
!!$  real(rp)    :: xlv, ep2,svp1,svp2,svp3,svpt0,rhowater,qvs,temp,RH_k,qv_k
!!$
!!$  character   :: fname*72
!!$
!!$  !
!!$  ! outvar = 1 --> densi and theta
!!$  ! outvar = 2 --> press and theta
!!$  !
!!$  integer(ip) :: outvar
!!$  real(rp)    :: q(ndofn_nsa+3,npoin), q_ref(ndofn_nsa+3,npoin), q_exact(ndofn_nsa,npoin)
!!$  !  real(rp)    :: q(7,npoin), q_ref(7,npoin), q_exact(4,npoin)
!!$  real(rp)    :: tmpt(nelz_nsa+1)
!!$
!!$  !
!!$  ! Sponge parameters:
!!$  !
!!$  real(rp)   :: aa(npoin),bb(npoin)
!!$  real(rp)   :: xsr,xsl,ysr,ysl,zs,dxs,dys,dzs
!!$  real(rp)   :: amp,ampx,ampy,ampz,ctop,cside
!!$
!!$  !
!!$  !Exact solution for mountains
!!$  !
!!$  real(rp)   ::  rho00, p_k, rhofac, afac, hafac, delta_k, temp0
!!$  real(rp)   ::  ddelta_dx_k, ddelta_dy_k, ddelta_dz_k, xappa
!!$  real(rp)   ::  zd,alpha,ct,cs,xbound,xl,dsx,dsy,zdcorrect,ztop
!!$  real(rp)   ::  dbl,beta,xid,abstaud,dt,dbt
!!$  integer(ip)::  nelx,nely,nelz
!!$  real(rp)   ::  rho_exact, pi_exact, u_exact, v_exact, w_exact, theta_exact
!!$  real(rp)   ::  rho_pert, p_pert, t_pert, theta_pert
!!$  real(rp)   ::  m, intercept
!!$
!!$  !
!!$  !  input values for squall-line test
!!$  !
!!$  xlv      = 2500000.0_rp   !Latent heat of vaporization
!!$  ep2      = 0.6217504_rp   !Rd/Rv
!!$  svp1     = 0.6112000_rp   
!!$  svp2     = 17.67000_rp
!!$  svp3     = 29.65000_rp
!!$  svpt0    = 273.1500_rp
!!$  rhowater = 1000.000_rp    !density of rainwater
!!$
!!$
!!$  select case(kfl_benme_nsa)
!!$
!!$  case (15) !Inertia-Gravity Wave
!!$     !Case Constants
!!$
!!$     temp0  = 300.0_rp
!!$     theta0 = tempe_nsa
!!$
!!$     xmin = 0.0_rp
!!$     xmax = 300000.0_rp
!!$     zmin = 0.0_rp
!!$     zmax = 10000.0_rp
!!$
!!$     c = cvcoe_nsa/rgasc_nsa
!!$     c2 = rgasc_nsa/cpcoe_nsa
!!$
!!$     u0 =20.0_rp
!!$     w0 = 0.0_rp
!!$
!!$     bv=0.01_rp
!!$     bv2=bv*bv
!!$
!!$     g2=grnor_nsa*grnor_nsa
!!$
!!$     ac = 5000.0_rp
!!$     xl = 300000.0_rp
!!$     zl = 10000.0_rp
!!$     xc = xl/3.0_rp
!!$     zc = zl/2.0_rp
!!$     thetac = 0.01_rp
!!$
!!$     do ipoin = 1,npoin
!!$        x = coord(1,ipoin)
!!$        z = coord(ndime,ipoin)
!!$
!!$        u_k = u0
!!$        w_k = w0
!!$
!!$        pi_hyd =g2/(cpcoe_nsa*temp0*bv2)*(exp(- bv2/grnor_nsa*z ) - 1.0_rp) + 1.0_rp
!!$        theta_hyd = temp0*exp( bv2/grnor_nsa*z )
!!$        dtheta    = thetac*sin( pi*z/zl )/( 1.0_rp + ( (x-xc)/ac )**2 )
!!$
!!$        theta_k   = theta_hyd + dtheta
!!$        rho_k     = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_hyd)**c
!!$        rho_hyd   = pbaro_nsa/(rgasc_nsa*theta_hyd)*(pi_hyd)**c
!!$        call nsa_stalaw(2,rho_k,press_k,theta_k,z,0.0_rp)
!!$        call nsa_stalaw(2,rho_hyd,press_hyd,theta_hyd,z,0.0_rp)
!!$
!!$        !Press_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho_hyd
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta_hyd
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press_hyd
!!$
!!$        !Perturbation fields: (total - reference)
!!$        bvess_nsa(1,ipoin,1)           = u0 !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo
!!$
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$     end do !ip
!!$
!!$     !Return density and theta
!!$     outvar = 1
!!$
!!$  case(1)
!!$     !
!!$     !Rising Thermal Bubble
!!$     !
!!$     !Case Constants
!!$     c=cvcoe_nsa/rgasc_nsa
!!$     theta0=tempe_nsa
!!$
!!$     xmin = 0.0_rp
!!$     xmax = 1000.0_rp
!!$     zmin = 0.0_rp
!!$     zmax = 1000.0_rp
!!$     !zmax = 1500.0_rp
!!$     if(ndime > 2) then
!!$        ymin = ymin_nsa
!!$        ymax = ymax_nsa
!!$     end if
!!$
!!$     u0=0.0_rp
!!$     w0=0.0_rp
!!$     xc=0.5_rp*(xmin + xmax)
!!$     zc=350.0_rp
!!$     if(ndime > 2) &
!!$          yc = 0.5_rp*(ymin + ymax)
!!$     rc=250.0_rp
!!$
!!$     !
!!$     ! Perturbation intensity:
!!$     !
!!$     thetac=0.5_rp   ! With bubble
!!$
!!$     do ipoin=1,npoin
!!$        x=coord(1,ipoin)
!!$        z=coord(ndime,ipoin)
!!$        if(ndime > 2) &
!!$             y=coord(ndime-1,ipoin)
!!$
!!$
!!$        dtheta=0.0_rp
!!$        if(ndime < 3) then
!!$           !
!!$           ! 2D
!!$           !
!!$           r=sqrt( (x-xc)**2.0_rp + (z-zc)**2.0_rp )
!!$           if (r < rc) then
!!$              dtheta = 0.5_rp*thetac*(1.0_rp + cos(pi*r/rc) )
!!$           end if
!!$        else
!!$           !
!!$           ! 3D
!!$           !
!!$           if( kfl_cylind_nsa > 0 ) then
!!$              !Cylinder
!!$              r=sqrt( (x-xc)**2.0_rp + (z-zc)**2.0_rp )
!!$           else
!!$              !Sphere
!!$              r=sqrt( (x-xc)**2.0_rp + (y-yc)**2_rp + (z-zc)**2.0_rp ) !Sphere
!!$           end if
!!$
!!$           if (r < rc) then
!!$              dtheta = 0.5_rp*thetac*(1.0_rp + cos(pi*r/rc) )
!!$           end if
!!$        end if
!!$
!!$
!!$        !Total fields:
!!$        theta_k = theta0 + dtheta
!!$        pi_k    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta0)*z
!!$        rho_k   = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_k)**c
!!$        call nsa_stalaw(2,rho_k,press_k,theta_k,z,0.0_rp)
!!$
!!$        !Reference Fields (hydrostatic equilibrium):
!!$        theta_hyd = tempe_nsa
!!$        pi_hyd    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta_hyd)*z
!!$        rho_hyd   = pbaro_nsa/(rgasc_nsa*theta_hyd)*(pi_hyd)**c
!!$        call nsa_stalaw(2,rho_hyd,press_hyd,theta_hyd,z,0.0_rp)
!!$
!!$        !Press_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho_hyd
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta_hyd
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press_hyd
!!$
!!$        !Perturbation fields: (total - reference)
!!$        bvess_nsa(1,ipoin,1)           = u0 !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo
!!$
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$     end do !ipoin
!!$
!!$     !Return density and theta
!!$     outvar = 1
!!$
!!$  case(11)
!!$     !
!!$     !Rising Thermal Bubble of Robert 1993
!!$     !
!!$     !Case Constants
!!$     c=cvcoe_nsa/rgasc_nsa
!!$     theta0=tempe_nsa
!!$
!!$     xmin = 0.0_rp
!!$     xmax = 1000.0_rp
!!$     zmin = 0.0_rp
!!$     zmax = 1000.0_rp
!!$     !zmax = 1500.0_rp
!!$
!!$     u0 = 0.0_rp
!!$     w0 = 0.0_rp
!!$     xc = 500.0_rp
!!$     zc = 260.0_rp
!!$     rc = 50.0_rp
!!$     sigma = 100.0_rp
!!$     sigma = sigma*sigma
!!$
!!$     !Perturbation intensity:
!!$     thetac=0.5_rp   ! With bubble
!!$
!!$     do ipoin=1,npoin
!!$        x=coord(1,ipoin)
!!$        y=coord(ndime-1,ipoin)
!!$        z=coord(ndime,ipoin)
!!$
!!$        if(ndime < 3) then
!!$           r=sqrt( (x-xc)**2.0_rp + (z-zc)**2.0_rp )
!!$           if (r > rc) then
!!$              dtheta = thetac*exp(-(r-rc)*(r-rc)/sigma)
!!$           else
!!$              dtheta = thetac
!!$           end if
!!$        else
!!$           !Cosine function as Kelly and Giraldo 2012 JCP 3D
!!$           xc = 500.0_rp
!!$           yc = 500.0_rp
!!$           zc = 260.0_rp
!!$           xr = 250.0_rp
!!$           yr = 250.0_rp
!!$           zr = 250.0_rp
!!$           rc =   1.0_rp
!!$
!!$           dtheta = 0.0_rp
!!$           r=sqrt( ((x-xc)/xr)**2_rp  + ((y-yc)/yr)**2_rp + ((z-zc)/zr)**2_rp )
!!$           if (r < rc) then
!!$              dtheta=thetac*(1.0_rp + cos(pi*r/rc))
!!$           end if
!!$        end if
!!$
!!$        !Total fields:
!!$        theta_k = theta0 + dtheta
!!$        pi_k    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta0)*z
!!$        rho_k   = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_k)**c
!!$        call nsa_stalaw(2,rho_k,press_k,theta_k,z,0.0_rp)
!!$
!!$        !Reference Fields (hydrostatic equilibrium):
!!$        theta_hyd = tempe_nsa
!!$        pi_hyd    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta_hyd)*z
!!$        rho_hyd   = pbaro_nsa/(rgasc_nsa*theta_hyd)*(pi_hyd)**c
!!$        call nsa_stalaw(2,rho_hyd,press_hyd,theta_hyd,z,0.0_rp)
!!$
!!$        !Press_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho_hyd
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta_hyd
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press_hyd
!!$
!!$        !Perturbation fields: (total - reference)
!!$        bvess_nsa(1,ipoin,1)           = u0 !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo
!!$
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$     end do !ipoin
!!$
!!$     !Return density and theta
!!$     outvar = 1
!!$
!!$  case(12)
!!$     !
!!$     !Rising Thermal Bubble as in Ahmad
!!$     !
!!$     !Case Constants
!!$     c=cvcoe_nsa/rgasc_nsa
!!$     theta0 = tempe_nsa
!!$
!!$     xmin = xmin_nsa
!!$     xmax = xmax_nsa
!!$     zmin = 0.0_rp
!!$     zmax = zmax_nsa
!!$     if(ndime > 2) then
!!$        ymin = ymin_nsa
!!$        ymax = ymax_nsa
!!$     end if
!!$
!!$     u0=0.0_rp
!!$     w0=0.0_rp
!!$
!!$     xc = 0.5_rp*(xmin + xmax)
!!$     zc = zc_nsa
!!$     if(ndime > 2) &
!!$          yc = 0.5_rp*(ymin + ymax)
!!$
!!$     rc = 2000.0_rp
!!$
!!$     xr = xradi_nsa
!!$     yr = yradi_nsa
!!$     zr = zradi_nsa
!!$     if(xradi_nsa < 0) xr = 2000.0_rp
!!$     if(yradi_nsa < 0) yr = 2000.0_rp
!!$     if(zradi_nsa < 0) zr = 2000.0_rp
!!$
!!$     !Perturbation intensity:
!!$     if(thetac_nsa < -999.0_rp) then
!!$        thetac=0.0_rp       ! With bubble
!!$     else
!!$        thetac=thetac_nsa   ! With bubble
!!$     end if
!!$
!!$     do ipoin=1,npoin
!!$        x=coord(1,ipoin)
!!$        z=coord(ndime,ipoin)
!!$        if(ndime > 2) &
!!$             y=coord(ndime-1,ipoin)
!!$
!!$        dtheta=0.0_rp
!!$        if(ndime < 3) then
!!$           r=sqrt( (x-xc)**2.0_rp + (z-zc)**2.0_rp )   
!!$           if (r <= rc) then
!!$              dtheta=thetac*(1.0_rp - r)
!!$           end if
!!$        else
!!$           !r=sqrt( ((x-xc)/xr)**2_rp  + ((y-yc)/yr)**2_rp + ((z-zc)/zr)**2_rp )
!!$           r=sqrt( ((x-xc)/xr)**2.0_rp + ((z-zc)/zr)**2.0_rp )   
!!$           if (r <= 1.0_rp) then
!!$              dtheta=thetac*(1.0_rp - r)
!!$           end if
!!$        end if
!!$
!!$        !Total fields:
!!$        theta_k = theta0 + dtheta
!!$        pi_k    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta0)*z
!!$        rho_k   = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_k)**c
!!$        call nsa_stalaw(2,rho_k,press_k,theta_k,z,0.0_rp)
!!$
!!$        !Reference Fields (hydrostatic equilibrium):
!!$        theta_hyd = tempe_nsa
!!$        pi_hyd    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta_hyd)*z
!!$        rho_hyd   = pbaro_nsa/(rgasc_nsa*theta_hyd)*(pi_hyd)**c
!!$        call nsa_stalaw(2,rho_hyd,press_hyd,theta_hyd,z,0.0_rp)
!!$
!!$        !Press_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho_hyd
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta_hyd
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press_hyd
!!$
!!$        !Perturbation fields: (total - reference)
!!$        bvess_nsa(1,ipoin,1)           = u0 !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0 !wvelo
!!$        if(ndime > 2) &
!!$             bvess_nsa(ndime-1,ipoin,1)       = v0 !vvelo
!!$
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$     end do !ipoin
!!$
!!$     !Return density and theta
!!$     outvar = 1
!!$
!!$  case(2)
!!$     !
!!$     ! Density current
!!$     !
!!$     c=cvcoe_nsa/rgasc_nsa
!!$     theta0=tempe_nsa
!!$     u0=0.0_rp
!!$     w0=0.0_rp
!!$
!!$     xmin = -25600.0_rp
!!$     xmax =  25600.0_rp
!!$     zmin =  0.0_rp
!!$     zmax =  6400.0_rp
!!$
!!$     xc=0.5_rp*(xmin + xmax)
!!$     zc=3000.0_rp
!!$     xr=4000.0_rp
!!$     zr=2000.0_rp
!!$     rc=1.0_rp
!!$
!!$     !Perturbation intensity:
!!$     thetac = -15.0_rp
!!$
!!$     do ipoin=1,npoin
!!$        x=coord(1,ipoin)
!!$        z=coord(ndime,ipoin)
!!$
!!$        r=sqrt( ((x-xc)/xr)**2 + ((z-zc)/zr)**2 )
!!$        dtheta=0
!!$        if (r <= rc) then
!!$           dtheta=thetac/2*(1 + cos(pi*r) )
!!$        end if
!!$        !Total fields:
!!$        theta_k = theta0 + dtheta
!!$        pi_k    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta0)*z
!!$        rho_k   = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_k)**c
!!$        call nsa_stalaw(2,rho_k,press_k,theta_k,z,0.0_rp)
!!$
!!$        !Reference Fields (hydrostatic equilibrium):
!!$        theta_hyd = tempe_nsa
!!$        pi_hyd    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta_hyd)*z
!!$        rho_hyd   = pbaro_nsa/(rgasc_nsa*theta_hyd)*(pi_hyd)**c
!!$        call nsa_stalaw(2,rho_hyd,press_hyd,theta_hyd,z,0.0_rp)
!!$
!!$        !Press_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho_hyd
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta_hyd
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press_hyd
!!$
!!$        !Perturbation fields: (total - reference)
!!$        bvess_nsa(1,ipoin,1)           = u0 !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo
!!$
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$     end do !ipoin
!!$
!!$     !Return density and theta
!!$     outvar = 1
!!$
!!$  case(3)
!!$     !
!!$     ! HS Linear mountain
!!$     !
!!$     !Case Constants
!!$     theta0    = tempe_nsa
!!$     temp0     = tempe_nsa
!!$     rho_baro  = 1.0_rp/(rgasc_nsa*theta0)*pbaro_nsa; ! surface density
!!$     rho00     = 1.0_rp/(rgasc_nsa*theta0)*pbaro_nsa; ! surface density
!!$
!!$     u0=20.0_rp
!!$     w0=0.0_rp
!!$
!!$     ac = 10000.0_rp
!!$     hc = 1.0_rp !Small Mountain
!!$     xc = 240000.0_rp/2.0_rp
!!$
!!$     if(kfl_botta_nsa > 0) then
!!$        !
!!$        !Botta and Klein
!!$        !
!!$        u0 = 20.0_rp
!!$        ac = 800.0_rp
!!$        hc = 2000.0_rp
!!$        xc = 0.0_rp
!!$     end if
!!$
!!$     n2 = grnor_nsa**2/(cpcoe_nsa*temp0); ! 1/s^2
!!$     l2 = n2/u0**2;              ! 1/m^2
!!$     l = sqrt(l2);               ! 1/m
!!$     xappa=0.286_rp
!!$
!!$     do ipoin=1,npoin
!!$        x=coord(1,ipoin)
!!$        z=coord(ndime,ipoin)
!!$
!!$        !Total variables
!!$        pi_k    = exp(-grnor_nsa/(cpcoe_nsa*theta0)*z);
!!$        theta_k = theta0/pi_k;
!!$        press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);
!!$        call nsa_stalaw(1,rho_k,press_k,theta_k,z,0.0_rp)
!!$
!!$        !Reference Solution1
!!$        rho_hyd   = rho_k
!!$        press_hyd = press_k
!!$        theta_hyd = theta_k
!!$
!!$        !Press_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho_hyd
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta_hyd
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press_hyd
!!$
!!$        if(kfl_ncons_nsa == 0) then
!!$           !
!!$           ! Momentum in Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0*rho_hyd
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0*rho_hyd
!!$        else
!!$           !
!!$           ! Momentum in NON-Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0
!!$
!!$        end if
!!$
!!$        !Perturbation fields: (total - reference)
!!$        bvess_nsa(1,ipoin,1)           = u0 !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo
!!$
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$     end do !ipoin
!!$
!!$     !Return press and theta
!!$     outvar = 1
!!$
!!$  case (4) !Linear Hydrostatic: Large Mountain
!!$
!!$     !Case Constants
!!$     theta0    = tempe_nsa
!!$     rho_baro  = 1.0_rp/(rgasc_nsa*theta0)*pbaro_nsa; ! surface density
!!$     u0 = 20.0_rp
!!$     w0 =  0.0_rp
!!$
!!$     do ipoin=1,npoin
!!$        x=coord(1,ipoin)
!!$        z=coord(ndime,ipoin)
!!$
!!$        !Total variables
!!$        pi_k    = exp(-grnor_nsa/(cpcoe_nsa*theta0)*z);
!!$        theta_k = theta0/pi_k;
!!$        press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);
!!$        call nsa_stalaw(1,rho_k,press_k,theta_k,z,0.0_rp)
!!$
!!$        !Reference Solution1
!!$        rho_hyd   = rho_k
!!$        press_hyd = press_k
!!$        theta_hyd = theta_k
!!$
!!$        !Press_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho_hyd
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta_hyd
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press_hyd
!!$
!!$        if(kfl_ncons_nsa == 0) then
!!$           !
!!$           ! Momentum in Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0*rho_hyd
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0*rho_hyd
!!$        else
!!$           !
!!$           ! Momentum in NON-Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0
!!$
!!$        end if
!!$
!!$        !Perturbation fields: (total - reference)
!!$        bvess_nsa(1,ipoin,1)           = u0 !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo
!!$
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$     end do
!!$
!!$     !Return press and theta
!!$     outvar = 1
!!$
!!$  case(5)
!!$     !
!!$     ! NH Linear mountain
!!$     !
!!$     !Case Constants
!!$     temp0 = tempe_nsa
!!$     theta0    = tempe_nsa
!!$     c     = cvcoe_nsa/rgasc_nsa
!!$     c2    = rgasc_nsa/cpcoe_nsa
!!$
!!$     u0    =  10.0_rp
!!$     w0    =  0.0_rp
!!$     bv    =  0.01_rp
!!$     bv2   = bv*bv
!!$     g2    = grnor_nsa*grnor_nsa
!!$
!!$     ac = 1000.0_rp
!!$     hc = 1.0_rp !Small Mountain
!!$     xc = 0.0_rp
!!$
!!$     do ipoin=1,npoin
!!$        x=coord(1,ipoin)
!!$        z=coord(ndime,ipoin)
!!$
!!$        !Total variables
!!$        pi_k    = g2/(cpcoe_nsa*theta0*bv2)*(exp(- bv2/grnor_nsa*z ) - 1.0_rp) + 1.0_rp
!!$        theta_k = theta0*exp( bv2/grnor_nsa*z )
!!$        rho_k   = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_k)**c
!!$        call nsa_stalaw(2,rho_k,press_k,theta_k,z,0.0_rp)
!!$        !press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);
!!$
!!$        !Reference Solution1
!!$        rho_hyd   = rho_k
!!$        press_hyd = press_k
!!$        theta_hyd = theta_k
!!$
!!$        !Press_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho_hyd
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta_hyd
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press_hyd
!!$
!!$        if(kfl_ncons_nsa == 0) then
!!$           !
!!$           ! Momentum in Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0*rho_hyd
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0*rho_hyd
!!$        else
!!$           !
!!$           ! Momentum in NON-Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0
!!$
!!$        end if
!!$
!!$        !Perturbation fields: (total - reference)
!!$        bvess_nsa(1,ipoin,1)           = u0 !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo
!!$
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$     end do !ipoin
!!$
!!$     !Return press and theta
!!$     outvar = 1
!!$
!!$  case(200)
!!$     !
!!$     ! Moisture: this case is similar to case 200
!!$     ! except that it already has a supersaturaded cloud 
!!$     ! in the initial field. qv initial is the same as for case 200.
!!$     !
!!$     if(ndime < 3) then
!!$        !
!!$        ! 2D:
!!$        !
!!$        nvar_nsa = 7_ip
!!$        c=cvcoe_nsa/rgasc_nsa
!!$        c2=rgasc_nsa/cpcoe_nsa
!!$        
!!$        !Center of thermal:
!!$        xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
!!$        zc = 2000.0_rp
!!$        if(ndime > 2) &
!!$             yc = 0.5_rp*(maxval(coord(2,:)) + minval(coord(2,:)))
!!$        
!!$        
!!$        !Reference theta and perturbation theta
!!$        theta0 = 300.0_rp
!!$        if(thetac_nsa > -999.0_rp) then
!!$           thetac = thetac_nsa
!!$        else
!!$           thetac = 3.0_rp
!!$        end if
!!$        if(tracerc_nsa > -999.0_rp) then
!!$           tracerc = tracerc_nsa
!!$        else
!!$           tracerc = 0.0_rp
!!$        end if
!!$        
!!$        !Radii of thermal and tracer
!!$        rc     = 1.0_rp
!!$        
!!$        !Thermal
!!$        xr     = 10000.0_rp
!!$        zr     = 1500.0_rp    
!!$        if(ndime > 2) &
!!$             yr = 10000.0_rp
!!$        
!!$        !Tracer
!!$        xr_tr  = 2500.0_rp
!!$        zr_tr  = 1000.0_rp
!!$        if(ndime > 2) &
!!$             yr_tr = 2500.0_rp
!!$        
!!$        !Load Reference Values
!!$        call nsa_loadref(q,q_ref)
!!$        
!!$        do ipoin=1,npoin
!!$           x=coord(1,ipoin)
!!$           z=coord(ndime,ipoin)
!!$           if(ndime > 2) &
!!$                y = coord(ndime-1,ipoin)
!!$           
!!$           !Background density and theta (from sounding):
!!$           rho0   = q_ref(1,ipoin)
!!$           theta0 = q_ref(4,ipoin)
!!$           press0 = pbaro_nsa*(rho0*rgasc_nsa*theta0/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)
!!$
!!$
!!$           if(ndime < 3) then
!!$              r=sqrt( ((x-xc)/xr)**2_rp  + ((z-zc)/zr)**2_rp )
!!$           else
!!$              r=sqrt( ((x-xc)/xr)**2_rp  + ((y-yc)/yr)**2_rp + ((z-zc)/zr)**2_rp )
!!$           end if
!!$
!!$           dtheta = 0.0_rp
!!$           qv     = 0.0_rp
!!$           if (r < rc) then
!!$              dtheta = thetac*(cos(pi*r/2.0_rp))**2_rp
!!$           end if
!!$
!!$           dtracer = 0.0_rp
!!$           if(ndime < 3) then
!!$              r=sqrt( ((x-xc)/xr_tr)**2  + ((z-zc)/zr_tr)**2 )
!!$           else
!!$              r=sqrt( ((x-xc)/xr_tr)**2_rp  + ((y-yc)/yr_tr)**2_rp + ((z-zc)/zr_tr)**2_rp )
!!$           end if
!!$
!!$           if (r < rc) then
!!$              dtracer = tracerc!*(cos(pi*r/2.0))**2
!!$           end if
!!$
!!$           theta_k = theta0 + dtheta
!!$           press_k = pbaro_nsa*(rho0*rgasc_nsa*theta_k/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)
!!$           rho_k   = pbaro_nsa*(press_k/pbaro_nsa)**(cvcoe_nsa/cpcoe_nsa)/(rgasc_nsa*theta_k)
!!$           !call nsa_stalaw(1,rho_k,press_k,theta_k,z,0.0_rp)
!!$
!!$           !
!!$           ! Total and perturbation Fields (for plotting only):
!!$           !
!!$           q(1,ipoin) = rho_k
!!$           if( kfl_uniformvelocity_nsa >= 0) then
!!$              q(2,ipoin) =  uvelo_nsa !if commented, horizontal velocity comes from sounding
!!$           else
!!$              !ELSE q takes the value read in the file.
!!$           end if
!!$
!!$           q(3,ipoin) = 0.0_rp
!!$           q(4,ipoin) = theta_k
!!$           q(5,ipoin) = q_ref(5,ipoin)
!!$           q(6,ipoin) = 0.0_rp
!!$           q(7,ipoin) = 0.0_rp
!!$
!!$           !Densi_hyd
!!$           rekee_nsa(ndime+1, ipoin) = rho0
!!$           !Theta_hyd
!!$           rekee_nsa(ndime+2, ipoin) = theta0
!!$           !Press_hyd
!!$           rekee_nsa(ndime+3, ipoin) = press0
!!$
!!$           !Water tracers
!!$           rekee_nsa(ndime+3 + 1, ipoin) = q_ref(5,ipoin)
!!$           rekee_nsa(ndime+3 + 2, ipoin) = 0.0_rp
!!$           rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp
!!$
!!$           if(kfl_ncons_nsa == 0) then
!!$              !
!!$              ! Momentum in Conservative form
!!$              ! 
!!$              !u-velo
!!$              rekee_nsa(1, ipoin)       = q(2,ipoin)*rho0
!!$              !w-velo
!!$              rekee_nsa(ndime, ipoin)   = w0*rho0
!!$           else
!!$              !
!!$              ! Momentum in NON-Conservative form
!!$              ! 
!!$              !u-velo
!!$              rekee_nsa(1, ipoin)       = q(2,ipoin)
!!$              !w-velo
!!$              rekee_nsa(ndime, ipoin)   = w0
!!$
!!$           end if
!!$
!!$           !Velocity:
!!$           bvess_nsa(1,ipoin,1)           = q(2,ipoin) !uvelo
!!$           bvess_nsa(ndime,ipoin,1)       = w0         !vvelo
!!$
!!$           !Total fields:
!!$           bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$           bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$           !Water tracers:
!!$           bvess_nsa(ndofn_nsa+1,ipoin,1)   = q(5,ipoin)
!!$           bvess_nsa(ndofn_nsa+2,ipoin,1)   = q(6,ipoin)
!!$           bvess_nsa(ndofn_nsa+3,ipoin,1)   = q(7,ipoin)
!!$
!!$        end do !ipoin
!!$
!!$        !Return press and theta
!!$        outvar = 1
!!$
!!$        !
!!$        ! Build the mapping from 2D grid to column.
!!$        ! This is done here because it's only used by the Kessler module:
!!$        !
!!$        call nsa_build_column
!!$
!!$     else
!!$        !
!!$        ! 3D:
!!$        !
!!$        nvar_nsa = 8_ip
!!$        c=cvcoe_nsa/rgasc_nsa
!!$        c2=rgasc_nsa/cpcoe_nsa
!!$
!!$        !Center of thermal:
!!$        xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
!!$        zc = 2000.0_rp
!!$        if(ndime > 2) &
!!$             yc = 0.5_rp*(maxval(coord(2,:)) + minval(coord(2,:)))
!!$
!!$
!!$        !Reference theta and perturbation theta
!!$        theta0 = 300.0_rp
!!$        if(thetac_nsa > -999.0_rp) then
!!$           thetac = thetac_nsa
!!$        else
!!$           thetac = 3.0_rp
!!$        end if
!!$        if(tracerc_nsa > -999.0_rp) then
!!$           tracerc = tracerc_nsa
!!$        else
!!$           tracerc = 0.0_rp
!!$        end if
!!$
!!$        !Radii of thermal and tracer
!!$        rc     = 1.0_rp
!!$
!!$        !Thermal
!!$        xr     = 10000.0_rp
!!$        zr     = 1500.0_rp    
!!$        if(ndime > 2) &
!!$             yr = 10000.0_rp
!!$
!!$        !Tracer
!!$        xr_tr  = 2500.0_rp
!!$        zr_tr  = 1000.0_rp
!!$        if(ndime > 2) &
!!$             yr_tr = 2500.0_rp
!!$
!!$        !Load Reference Values
!!$        call nsa_loadref(q,q_ref)
!!$
!!$        do ipoin=1,npoin
!!$           x=coord(1,ipoin)
!!$           z=coord(ndime,ipoin)
!!$           if(ndime > 2) &
!!$                y = coord(ndime-1,ipoin)
!!$
!!$           !Background density and theta (from sounding):
!!$           rho0   = q_ref(1,ipoin)
!!$           theta0 = q_ref(5,ipoin)
!!$           press0 = pbaro_nsa*(rho0*rgasc_nsa*theta0/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)
!!$           
!!$           if( kfl_cylind_nsa > 0) then
!!$              !Cylinder:
!!$              r=sqrt( ((x-xc)/xr)**2_rp + ((z-zc)/zr)**2_rp )
!!$           else
!!$              !Sphere:
!!$              r=sqrt( ((x-xc)/xr)**2_rp  + ((y-yc)/yr)**2_rp + ((z-zc)/zr)**2_rp ) !sphere
!!$           end if
!!$
!!$           dtheta = 0.0_rp
!!$           qv     = 0.0_rp
!!$           if (r < rc) then
!!$              dtheta = thetac*(cos(pi*r/2.0_rp))**2_rp
!!$           end if
!!$
!!$           dtracer = 0.0_rp
!!$           if(ndime < 3) then
!!$              r=sqrt( ((x-xc)/xr_tr)**2  + ((z-zc)/zr_tr)**2 )
!!$           else
!!$              !r=sqrt( ((x-xc)/xr_tr)**2_rp  + ((y-yc)/yr_tr)**2_rp + ((z-zc)/zr_tr)**2_rp ) !sphere
!!$              r=sqrt( ((x-xc)/xr_tr)**2_rp  + ((z-zc)/zr_tr)**2_rp ) !Cylinder
!!$           end if
!!$
!!$           if (r < rc) then
!!$              dtracer = tracerc!*(cos(pi*r/2.0))**2
!!$           end if
!!$
!!$           theta_k = theta0 + dtheta
!!$           press_k = pbaro_nsa*(rho0*rgasc_nsa*theta_k/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)
!!$           rho_k   = pbaro_nsa*(press_k/pbaro_nsa)**(cvcoe_nsa/cpcoe_nsa)/(rgasc_nsa*theta_k)
!!$           !call nsa_stalaw(1,rho_k,press_k,theta_k,z,0.0_rp)
!!$
!!$           !
!!$           ! Total and perturbation Fields (for plotting only):
!!$           !
!!$           q(1,ipoin) = rho_k
!!$           if( kfl_uniformvelocity_nsa >= 0) then
!!$              q(2,ipoin) =  uvelo_nsa !if commented, horizontal velocity comes from sounding
!!$           else
!!$              !ELSE q takes the value read in the file.
!!$           end if
!!$
!!$           q(3,ipoin) = 0.0_rp  !vvelo
!!$           q(4,ipoin) = 0.0_rp  !wvelo
!!$           q(5,ipoin) = theta_k
!!$           q(6,ipoin) = q_ref(6,ipoin)
!!$           q(7,ipoin) = 0.0_rp
!!$           q(8,ipoin) = 0.0_rp
!!$
!!$           !Densi_hyd
!!$           rekee_nsa(ndime+1, ipoin) = rho0
!!$           !Theta_hyd
!!$           rekee_nsa(ndime+2, ipoin) = theta0
!!$           !Press_hyd
!!$           rekee_nsa(ndime+3, ipoin) = press0
!!$
!!$           !Water tracers
!!$           rekee_nsa(ndime+3 + 1, ipoin) = q_ref(6,ipoin)
!!$           rekee_nsa(ndime+3 + 2, ipoin) = 0.0_rp
!!$           rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp
!!$
!!$           if(kfl_ncons_nsa == 0) then
!!$              !
!!$              ! Momentum in Conservative form
!!$              ! 
!!$              !u-velo
!!$              rekee_nsa(1, ipoin)       = q(2,ipoin)*rho0
!!$              !w-velo
!!$              rekee_nsa(ndime, ipoin)   = w0*rho0
!!$           else
!!$              !
!!$              ! Momentum in NON-Conservative form
!!$              ! 
!!$              !u-velo
!!$              rekee_nsa(1, ipoin)       = q(2,ipoin)
!!$              !w-velo
!!$              rekee_nsa(ndime, ipoin)   = w0
!!$
!!$           end if
!!$
!!$           !Velocity:
!!$           bvess_nsa(1,ipoin,1)           = q(2,ipoin) !uvelo
!!$           bvess_nsa(ndime,ipoin,1)       = w0         !vvelo
!!$
!!$           !Total fields:
!!$           bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$           bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$           !Water tracers:
!!$           bvess_nsa(ndofn_nsa+1,ipoin,1)   = q(6,ipoin)
!!$           bvess_nsa(ndofn_nsa+2,ipoin,1)   = q(7,ipoin)
!!$           bvess_nsa(ndofn_nsa+3,ipoin,1)   = q(8,ipoin)
!!$
!!$        end do !ipoin
!!$
!!$        !Return press and theta
!!$        outvar = 1
!!$
!!$        !
!!$        ! Build the mapping from 2D grid to column.
!!$        ! This is done here because it's only used by the Kessler module:
!!$        !
!!$        call nsa_build_column
!!$
!!$     end if
!!$
!!$  case (204) !Dynamics as in LH mountain (large mountain), but also with moisture
!!$
!!$     !
!!$     ! HS Linear mountain with tracer (nvar_nsa = 7)
!!$     !
!!$     !Case Constants
!!$     nvar_nsa  = 7_ip
!!$     theta0    = tempe_nsa
!!$     temp0     = tempe_nsa
!!$     rho_baro  = 1.0_rp/(rgasc_nsa*theta0)*pbaro_nsa; ! surface density
!!$     rho00     = 1.0_rp/(rgasc_nsa*theta0)*pbaro_nsa; ! surface density
!!$
!!$     u0=20.0_rp
!!$     w0=0.0_rp
!!$
!!$     ac = 10000.0_rp
!!$     hc = 1.0_rp !Small Mountain
!!$     xc = 240000.0_rp/2.0_rp
!!$
!!$     n2 = grnor_nsa**2/(cpcoe_nsa*temp0); ! 1/s^2
!!$     l2 = n2/u0**2;                       ! 1/m^2
!!$     l = sqrt(l2);                        ! 1/m
!!$     xappa=0.286_rp
!!$
!!$     !Load Reference Values
!!$     call nsa_loadref(q,q_ref)
!!$
!!$     do ipoin=1,npoin
!!$        x=coord(1,ipoin)
!!$        z=coord(ndime,ipoin)
!!$
!!$        !Total variables
!!$        pi_k    = exp(-grnor_nsa/(cpcoe_nsa*theta0)*z);
!!$        theta_k = theta0/pi_k;
!!$        press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);
!!$        call nsa_stalaw(1,rho_k,press_k,theta_k,z,0.0_rp)
!!$
!!$        !Reference Solution1
!!$        rho_hyd   = rho_k
!!$        press_hyd = press_k
!!$        theta_hyd = theta_k
!!$
!!$        !Press_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho_hyd
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta_hyd
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press_hyd
!!$
!!$        !Water tracers
!!$        q(5,ipoin) = q_ref(5,ipoin)
!!$        q(6,ipoin) = 0.0_rp
!!$        q(7,ipoin) = 0.0_rp
!!$
!!$        if(kfl_ncons_nsa == 0) then
!!$           !
!!$           ! Momentum in Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0*rho_hyd
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0*rho_hyd
!!$        else
!!$           !
!!$           ! Momentum in NON-Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0
!!$
!!$        end if
!!$
!!$        !Perturbation fields: (total - reference)
!!$        bvess_nsa(1,ipoin,1)           = u0 !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo
!!$
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$        !Water tracers:
!!$        bvess_nsa(ndofn_nsa+1,ipoin,1)   = q(5,ipoin)
!!$        bvess_nsa(ndofn_nsa+2,ipoin,1)   = q(6,ipoin)
!!$        bvess_nsa(ndofn_nsa+3,ipoin,1)   = q(7,ipoin)
!!$     end do !ipoin
!!$
!!$     !Return press and theta
!!$     outvar = 1
!!$
!!$  case(205)
!!$     !
!!$     ! NH Linear mountain with tracer (nvar_nsa = 7)
!!$     !
!!$     !Case Constants
!!$     nvar_nsa = ndofn_nsa + 3_ip
!!$     temp0    = tempe_nsa
!!$     theta0   = tempe_nsa
!!$     c        = cvcoe_nsa/rgasc_nsa
!!$     c2       = rgasc_nsa/cpcoe_nsa
!!$
!!$     u0       = 10.0_rp
!!$     w0       =  0.0_rp
!!$     bv       =  0.01_rp
!!$     bv2      = bv*bv
!!$     g2       = grnor_nsa*grnor_nsa
!!$
!!$     ac = 1000.0_rp
!!$     hc = 2500.0_rp !Small Mountain
!!$     xc = 0.0
!!$
!!$     rc     = 1.0_rp
!!$     xr     = 10000.0_rp
!!$     zr     = 1500.0_rp
!!$     xr_tr  = 2500.0_rp
!!$     zr_tr  = 1000.0_rp
!!$
!!$     theta0 = 300.0_rp
!!$     if(thetac_nsa > -999.0_rp) then
!!$        thetac = thetac_nsa
!!$     else
!!$        thetac = 3.0_rp
!!$     end if
!!$     if(tracerc_nsa > -999.0_rp) then
!!$        tracerc = tracerc_nsa
!!$     else
!!$        tracerc = 0.0_rp
!!$     end if
!!$
!!$     !Load Reference Values
!!$     call nsa_loadref(q,q_ref)
!!$
!!$     do ipoin=1,npoin
!!$        x=coord(1,ipoin)
!!$        z=coord(ndime,ipoin)
!!$
!!$        !Total variables
!!$        pi_k    = g2/(cpcoe_nsa*theta0*bv2)*(exp(- bv2/grnor_nsa*z ) - 1.0_rp) + 1.0_rp
!!$        theta_k = theta0*exp( bv2/grnor_nsa*z )
!!$        rho_k   = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_k)**c
!!$        call nsa_stalaw(2,rho_k,press_k,theta_k,z,0.0_rp)
!!$        !press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);
!!$
!!$        !Reference Solution1
!!$        rho_hyd   = rho_k
!!$        press_hyd = press_k
!!$        theta_hyd = theta_k
!!$
!!$        !
!!$        ! Total and perturbation Fields (for plotting only):
!!$        !
!!$        q(1,ipoin) = rho_k
!!$        if( kfl_uniformvelocity_nsa > 0) then
!!$           q(2,ipoin) =  uvelo_nsa !if commented, horizontal velocity comes from sounding
!!$        else
!!$           !ELSE q takes the value read in the file.
!!$        end if
!!$        q(3,ipoin) = 0.0_rp
!!$        q(4,ipoin) = theta_k
!!$        q(5,ipoin) = q_ref(5,ipoin)
!!$        q(6,ipoin) = 0.0_rp
!!$        q(7,ipoin) = 0.0_rp
!!$
!!$
!!$        !Press_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho_hyd
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta_hyd
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press_hyd
!!$
!!$        !Water tracers
!!$        rekee_nsa(ndime+3 + 1, ipoin) = q_ref(5,ipoin)
!!$        rekee_nsa(ndime+3 + 2, ipoin) = 0.0_rp
!!$        rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp
!!$
!!$        if(kfl_ncons_nsa == 0) then
!!$           !
!!$           ! Momentum in Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0*rho_hyd
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0*rho_hyd
!!$        else
!!$           !
!!$           ! Momentum in NON-Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0
!!$
!!$        end if
!!$
!!$        !Perturbation fields: (total - reference)
!!$        bvess_nsa(1,ipoin,1)           = u0 !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo
!!$
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$        !Water tracers:
!!$        bvess_nsa(ndofn_nsa+1,ipoin,1)   = q(5,ipoin)
!!$        bvess_nsa(ndofn_nsa+2,ipoin,1)   = q(6,ipoin)
!!$        bvess_nsa(ndofn_nsa+3,ipoin,1)   = q(7,ipoin)
!!$
!!$     end do !ipoin
!!$
!!$     !Return press and theta
!!$     outvar = 1
!!$
!!$     !
!!$     ! Build the mapping from 2D grid to column.
!!$     ! This is done here because it's only used by the Kessler module:
!!$     !
!!$     call nsa_build_column
!!$     !
!!$
!!$  case(202)
!!$     !
!!$     ! Moisture: this case is similar to case 200
!!$     ! except that computes the HS state from analytic equations,
!!$     ! and computes theta_k as 
!!$     ! theta_k = theta_HYD + (theta0 - theta_HYD) + dtheta
!!$     ! The same applies to p_k and rho_k
!!$     !
!!$     nvar_nsa = ndofn_nsa + 3_ip
!!$     c=cvcoe_nsa/rgasc_nsa
!!$     c2=rgasc_nsa/cpcoe_nsa
!!$
!!$     xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
!!$     zc = 2000.0_rp
!!$
!!$     if(thetac_nsa > -999.0_rp) then
!!$        thetac = thetac_nsa
!!$     else
!!$        thetac = 3.0_rp
!!$     end if
!!$     if(tracerc_nsa > -999.0_rp) then
!!$        tracerc = tracerc_nsa
!!$     else
!!$        tracerc = 0.0_rp
!!$     end if
!!$
!!$     rc     = 1.0_rp
!!$     xr     = 10000.0_rp
!!$     zr     = 1500.0_rp
!!$     xr_tr  = 2500.0_rp
!!$     zr_tr  = 1000.0_rp
!!$
!!$     !Load Reference Values
!!$     call nsa_loadref(q,q_ref)
!!$
!!$     do ipoin=1,npoin
!!$        x=coord(1,ipoin)
!!$        z=coord(2,ipoin)
!!$
!!$        !Reference Fields (hydrostatic equilibrium):
!!$        theta_hyd = tempe_nsa
!!$        pi_hyd    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta_hyd)*z
!!$        rho_hyd   = pbaro_nsa/(rgasc_nsa*theta_hyd)*(pi_hyd)**c
!!$        call nsa_stalaw(2,rho_hyd,press_hyd,theta_hyd,z,0.0_rp)
!!$
!!$        !Background density and theta (from sounding):
!!$        rho0   = q_ref(1,ipoin)
!!$        theta0 = q_ref(4,ipoin)
!!$        press0 = pbaro_nsa*(rho0*rgasc_nsa*theta0/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)
!!$
!!$        r=sqrt( ((x-xc)/xr)**2_rp  + ((z-zc)/zr)**2_rp )
!!$        dtheta = 0.0_rp
!!$        qv     = 0.0_rp
!!$        if (r < rc) then
!!$           dtheta = thetac*(cos(pi*r/2.0_rp))**2_rp
!!$        end if
!!$
!!$        dtracer = 0.0_rp
!!$        r=sqrt( ((x-xc)/xr_tr)**2  + ((z-zc)/zr_tr)**2 )
!!$        if (r < rc) then
!!$           dtracer = tracerc!*(cos(pi*r/2.0))**2
!!$        end if
!!$
!!$        theta_k = theta_hyd + (theta0 - theta_hyd) + dtheta
!!$        press_k = pbaro_nsa*(rho0*rgasc_nsa*theta_k/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)
!!$        rho_k   = pbaro_nsa*(press_k/pbaro_nsa)**(cvcoe_nsa/cpcoe_nsa)/(rgasc_nsa*theta_k)
!!$        !
!!$        ! Total and perturbation Fields (for plotting only):
!!$        !
!!$        q(1,ipoin) = rho_k
!!$        !q(2,ipoin) = 0.0_rp !if commented, horizontal velocity comes from sounding
!!$        q(3,ipoin) = 0.0_rp
!!$        q(4,ipoin) = theta_k
!!$        q(5,ipoin) = q_ref(5,ipoin)
!!$        q(6,ipoin) = 0.0_rp
!!$        q(7,ipoin) = 0.0_rp
!!$
!!$        !Initial values of tracers:
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho_hyd
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta_hyd
!!$        !Press_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press_hyd
!!$
!!$        !Water tracers
!!$        rekee_nsa(ndime+3 + 1, ipoin) = q_ref(5,ipoin)
!!$        rekee_nsa(ndime+3 + 2, ipoin) = 0.0_rp
!!$        rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp
!!$
!!$        if(kfl_ncons_nsa == 0) then
!!$           !
!!$           ! Momentum in Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = q(2,ipoin)*rho0
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0*rho0
!!$        else
!!$           !
!!$           ! Momentum in NON-Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = q(2,ipoin)
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0
!!$
!!$        end if
!!$
!!$        !Velocity:
!!$        bvess_nsa(1,ipoin,1)           = q(2,ipoin) !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0         !vvelo
!!$
!!$        !Total fields:
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$        !Water tracers:
!!$        bvess_nsa(ndofn_nsa+1,ipoin,1)   = q(5,ipoin)
!!$        bvess_nsa(ndofn_nsa+2,ipoin,1)   = q(6,ipoin)
!!$        bvess_nsa(ndofn_nsa+3,ipoin,1)   = q(7,ipoin)
!!$
!!$     end do !ipoin
!!$
!!$     !Return press and theta
!!$     outvar = 1
!!$
!!$     !
!!$     ! Build the mapping from 2D grid to column.
!!$     ! This is done here because it's only used by the Kessler module:
!!$     !
!!$     call nsa_build_column
!!$
!!$  case(201)
!!$     !
!!$     ! KESSLER, SIMPLE
!!$     ! Moisture: this case is similar to case 200
!!$     ! except that it already has a supersaturaded cloud 
!!$     ! in the initial field. qv initial is the same as for case 200.
!!$     !
!!$     nvar_nsa = ndofn_nsa + 3_ip
!!$     c=cvcoe_nsa/rgasc_nsa
!!$     c2=rgasc_nsa/cpcoe_nsa
!!$
!!$     xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
!!$     !zc = 1200.0_rp
!!$     zc = 2000.0_rp     
!!$     theta0 = 250.0_rp
!!$
!!$     if(thetac_nsa > -999.0_rp) then
!!$        thetac = thetac_nsa
!!$     else
!!$        thetac = 3.0_rp
!!$     end if
!!$     if(tracerc_nsa > -999.0_rp) then
!!$        tracerc = tracerc_nsa
!!$     else
!!$        tracerc = 0.0_rp
!!$     end if
!!$
!!$     !     rc     = 1.0_rp
!!$     !     xr     = 10000.0_rp
!!$     !     zr     = 1500.0_rp
!!$     !     xr_tr  = 2500.0_rp
!!$     !     zr_tr  = 1000.0_rp
!!$
!!$     rc     = 1.0_rp
!!$     xr     = 500.0_rp
!!$     zr     = 250.0_rp
!!$     xr_tr  = 100.0_rp
!!$     zr_tr  = 25.0_rp
!!$
!!$     !Load Reference Values
!!$     call nsa_loadref(q,q_ref)
!!$
!!$     do ipoin=1,npoin
!!$        x=coord(1,ipoin)
!!$        z=coord(2,ipoin)
!!$
!!$        !Background density and theta (from sounding):
!!$        rho0   = q_ref(1,ipoin)
!!$        theta0 = q_ref(4,ipoin)
!!$        press0 = pbaro_nsa*(rho0*rgasc_nsa*theta0/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)
!!$
!!$        r=sqrt( ((x-xc)/xr)**2_rp  + ((z-zc)/zr)**2_rp )
!!$        dtheta = 0.0_rp
!!$        qv     = 0.0_rp
!!$        if (r < rc) then
!!$           dtheta = thetac*(cos(pi*r/2.0_rp))**2_rp
!!$        end if
!!$
!!$        dtracer = 0.0_rp
!!$        r=sqrt( ((x-xc)/xr_tr)**2  + ((z-zc)/zr_tr)**2 )
!!$        if (r < rc) then
!!$           dtracer = tracerc!*(cos(pi*r/2.0))**2
!!$        end if
!!$
!!$        theta_k = theta0 + dtheta
!!$        press_k = pbaro_nsa*(rho0*rgasc_nsa*theta_k/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)
!!$        rho_k   = pbaro_nsa*(press_k/pbaro_nsa)**(cvcoe_nsa/cpcoe_nsa)/(rgasc_nsa*theta_k)
!!$        !call nsa_stalaw(1,rho_k,press_k,theta_k,z,0.0_rp)
!!$
!!$        !
!!$        ! Total and perturbation Fields (for plotting only):
!!$        !
!!$        q(1,ipoin) = rho_k
!!$        q(2,ipoin) = 0.0 !if commented, horizontal velocity comes from sounding
!!$        q(3,ipoin) = 0.0_rp
!!$        q(4,ipoin) = theta_k
!!$        q(5,ipoin) = q_ref(5,ipoin)
!!$        q(6,ipoin) = 0.0_rp!dtracer
!!$        q(7,ipoin) = 0.0_rp
!!$
!!$        !Initial values of tracers:
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho0
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta0
!!$        !Press_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press0
!!$
!!$        !Water tracers
!!$        rekee_nsa(ndime+3 + 1, ipoin) = q_ref(5,ipoin)
!!$        rekee_nsa(ndime+3 + 2, ipoin) = dtracer
!!$        rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp
!!$
!!$        if(kfl_ncons_nsa == 0) then
!!$           !
!!$           ! Momentum in Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = q(2,ipoin)*rho0
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0*rho0
!!$        else
!!$           !
!!$           ! Momentum in NON-Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = q(2,ipoin)
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0
!!$
!!$        end if
!!$
!!$
!!$        !Velocity:
!!$        bvess_nsa(1,ipoin,1)           = q(2,ipoin) !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0         !vvelo
!!$
!!$        !Total fields:
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$        !Water tracers:
!!$        bvess_nsa(ndofn_nsa+1,ipoin,1)   = q(5,ipoin)
!!$        bvess_nsa(ndofn_nsa+2,ipoin,1)   = q(6,ipoin)
!!$        bvess_nsa(ndofn_nsa+3,ipoin,1)   = q(7,ipoin)
!!$
!!$     end do !ipoin
!!$
!!$     !Return press and theta
!!$     outvar = 1
!!$
!!$     !
!!$     ! Build the mapping from 2D grid to column.
!!$     ! This is done here because it's only used by the Kessler module:
!!$     !
!!$     call nsa_build_column
!!$
!!$  case(203)
!!$     !
!!$     ! Klaassen and Klark 1985:
!!$     !
!!$     nvar_nsa = ndofn_nsa + 3_ip
!!$     c=cvcoe_nsa/rgasc_nsa
!!$     c2=rgasc_nsa/cpcoe_nsa
!!$
!!$     u0 = 0.0_rp
!!$     w0 = 0.0_rp
!!$
!!$     theta_z0 = 250.16_rp
!!$
!!$     if(thetac_nsa > -999.0_rp) then
!!$        thetac = thetac_nsa
!!$     else
!!$        thetac = 3.0_rp
!!$     end if
!!$     if(tracerc_nsa > -999.0_rp) then
!!$        tracerc = tracerc_nsa
!!$     else
!!$        tracerc = 0.0_rp
!!$     end if
!!$
!!$
!!$     xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
!!$     zc = 500.0_rp
!!$
!!$     rc     = 1.0_rp
!!$     xr     = 1500.0_rp
!!$     zr     = 750.0_rp
!!$
!!$     !Load Reference Values
!!$     !call nsa_loadref(q,q_ref)
!!$
!!$     do ipoin=1,npoin
!!$        x=coord(1,ipoin)
!!$        z=coord(2,ipoin)
!!$
!!$        !Background theta0:
!!$        !computed from the points of the graph in Klaassen and Klark 1985:
!!$        !Background density and theta (from sounding):
!!$        if( z <= 1400.0_rp) then
!!$           dtdz = 0.36e-03        
!!$           theta0 = dtdz*(z - zmin) + theta_z0
!!$
!!$        else if( z > 1400.0_rp .and. z <= 1800.0_rp) then
!!$           dtdz = 4.0e-03
!!$           zmin = 1400.0_rp
!!$           theta_ref = 0.36e-03*1400.0_rp + theta_z0
!!$           theta0 = dtdz*(z - zmin) + theta_ref
!!$
!!$        else if( z > 1800.0_rp .and. z <= 5000.0_rp) then
!!$           dtdz = 4.5e-03
!!$           zmin = 1800.0_rp
!!$           theta_ref = 4.0e-03*(1800.0_rp - 1400.0_rp) + 0.36e-03*1400.0_rp + theta_z0
!!$           theta0 = dtdz*(z - zmin) + theta_ref
!!$        end if
!!$
!!$        !Background water vapor:
!!$        !Background density and theta (from sounding):
!!$        if( z <= 1400.0_rp) then
!!$
!!$           call nsa_compline(0.0_rp,1400.0_rp,0.0085_rp,0.008_rp, m, intercept)
!!$           qv_k = m*z + intercept
!!$
!!$        else if( z > 1400.0_rp .and. z <= 1800.0_rp) then
!!$
!!$           call nsa_compline(1400.0_rp,1800.0_rp, 0.008_rp,0.0028_rp, m, intercept)
!!$           qv_k = m*z + intercept
!!$
!!$        else if( z > 1800.0_rp .and. z <= 3000.0_rp) then
!!$
!!$           call nsa_compline(1800.0_rp,3000.0_rp, 0.0028_rp,0.0015_rp, m, intercept)
!!$           qv_k = m*z + intercept
!!$
!!$        else if( z > 3000.0_rp .and. z <= 5000.0_rp) then
!!$
!!$           call nsa_compline(3000.0_rp,5000.0_rp, 0.0015_rp,0.0005_rp, m, intercept)
!!$           qv_k = m*z + intercept
!!$
!!$        end if
!!$
!!$        !PErturbation of theta_k
!!$        r=sqrt( ((x-xc)/xr)**2_rp  + ((z-zc)/zr)**2_rp )
!!$        dtheta = 0.0_rp
!!$        if (r < rc) then
!!$           dtheta = thetac*(cos(pi*r/2.0_rp))**2_rp
!!$        end if
!!$
!!$        !Total variables:
!!$        theta_k = theta0
!!$        pi_k    = theta_z0/theta_k
!!$        press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);
!!$        call nsa_stalaw(1,rho_k,press_k,theta_k,z,0.0_rp)
!!$
!!$        !Reference Solution1
!!$        rho_hyd   = rho_k
!!$        press_hyd = press_k
!!$        theta_hyd = theta_k
!!$
!!$        !Press_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho_hyd
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta_hyd
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press_hyd
!!$
!!$        !Water tracers:
!!$        q(5,ipoin) = qv_k
!!$        q(6,ipoin) = 0.0_rp
!!$        q(7,ipoin) = 0.0_rp
!!$
!!$        !
!!$        ! Add perturbation to theta and correct the other variables accordingly:
!!$        ! 
!!$        theta_k = theta0 + dtheta
!!$        pi_k    = theta_z0/theta_k
!!$        press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);
!!$        call nsa_stalaw(1,rho_k,press_k,theta_k,z,0.0_rp)
!!$
!!$        !Water tracers
!!$        rekee_nsa(ndime+3 + 1, ipoin) = qv_k
!!$        rekee_nsa(ndime+3 + 2, ipoin) = 0.0_rp
!!$        rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp
!!$
!!$        if(kfl_ncons_nsa == 0) then
!!$           !
!!$           ! Momentum in Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0*rho_hyd
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0*rho_hyd
!!$        else
!!$           !
!!$           ! Momentum in NON-Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0
!!$
!!$        end if
!!$
!!$        !Perturbation fields: (total - reference)
!!$        bvess_nsa(1,ipoin,1)           = u0 !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo
!!$
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
!!$
!!$        !Water tracers:
!!$        bvess_nsa(ndofn_nsa+1,ipoin,1)   = q(5,ipoin)
!!$        bvess_nsa(ndofn_nsa+2,ipoin,1)   = q(6,ipoin)
!!$        bvess_nsa(ndofn_nsa+3,ipoin,1)   = q(7,ipoin)
!!$     end do !ipoin
!!$
!!$     !Return press and theta
!!$     outvar = 1
!!$
!!$     !
!!$     ! Build the mapping from 2D grid to column.
!!$     ! This is done here because it's only used by the Kessler module:
!!$     !
!!$     call nsa_build_column
!!$
!!$  case(210)
!!$     !
!!$     ! Moist atmosphere in a stably stratified environment: 
!!$     ! Grabowski JAS 2007 and Grbowski & Clark 1991
!!$     ! 
!!$
!!$     ! Case Constants
!!$     nvar_nsa = ndofn_nsa + 3_ip
!!$
!!$     temp0    = 283.0_rp
!!$     theta0   = 283.0_rp  
!!$     tracer0  = 0.1_rp    !background value
!!$     if(thetac_nsa > -999.0_rp) then
!!$        thetac = thetac_nsa
!!$     else
!!$        thetac = 0.05_rp
!!$     end if
!!$     if(tracerc_nsa > -999.0_rp) then
!!$        tracerc = tracerc_nsa
!!$     else
!!$        tracerc = 0.95_rp
!!$     end if
!!$
!!$     !p00   = 85000.0_rp
!!$     p00   = pbaro_nsa
!!$
!!$     xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
!!$     zc = 1200.0_rp
!!$
!!$     rc    = 500.0_rp
!!$     rc0   = 100.0_rp
!!$     rc1   = 200.0_rp
!!$     rc2   = 300.0_rp
!!$     c     = cvcoe_nsa/rgasc_nsa
!!$     c2    = rgasc_nsa/cpcoe_nsa
!!$
!!$     u0    = 0.0_rp
!!$     w0    = 0.0_rp
!!$
!!$     static_stability = 1.0e-5_rp
!!$     bv    = sqrt(static_stability*grnor_nsa)
!!$     bv2   = bv*bv
!!$     g2    = grnor_nsa*grnor_nsa
!!$
!!$     do ipoin=1,npoin
!!$        x=coord(1,ipoin)
!!$        z=coord(ndime,ipoin)
!!$
!!$        !Total variables
!!$        pi_k    = g2/(cpcoe_nsa*theta0*bv2)*(exp(- bv2/grnor_nsa*z ) - 1.0_rp) + 1.0_rp
!!$        theta_k = theta0*exp( bv2/grnor_nsa*z )
!!$        rho_k   = p00/(rgasc_nsa*theta_k)*(pi_k)**c
!!$        call nsa_stalaw(2,rho_k,press_k,theta_k,z,0.0_rp)
!!$
!!$        !Reference Solution1
!!$        rho_hyd   = rho_k
!!$        press_hyd = press_k
!!$        theta_hyd = theta_k
!!$
!!$        !Perturbation of RH:
!!$        dtracer = 0.0_rp
!!$        r=sqrt( ((x-xc))**2  + ((z-zc))**2 )
!!$        if(r <= 400.0_rp) then
!!$           dtracer = 0.9_rp
!!$        else if (r > 400.0_rp .and. r <= 500.0_rp) then
!!$           dtracer = tracerc*(1.0_rp - r/500.0_rp)
!!$        end if
!if(r <= rc1) then
!   dtracer = tracerc         
!else if (r > rc1 .and. r <= rc2) then
!   dtracer = tracerc*(cos(0.5_rp*pi * (r - rc1)/rc0))**2
!end if
!!$        tracer_k = tracer0 + dtracer
!!$
!!$        !Relative humidity:
!!$        RH_k = tracer_k
!!$
!!$        temp     = theta_k*pi_k
!!$        pressure = press_k
!!$
!!$        es  = 1000.0_rp*svp1*exp(svp2*(temp-svpt0)/(temp-svp3))
!!$        qvs = ep2*es/(pressure-es)
!!$
!!$        !Derive qv from relative humidity:
!!$        qv_k = RH_k*qvs
!!$
!!$        dtheta = 0.0_rp
!!$        if (r < rc) then
!!$           dtheta = thetac*(cos(pi*r/2.0_rp))**2_rp
!!$        end if
!!$
!!$        !Press_hyd
!!$        rekee_nsa(ndime+1, ipoin) = rho_hyd
!!$        !Theta_hyd
!!$        rekee_nsa(ndime+2, ipoin) = theta_hyd
!!$        !Densi_hyd
!!$        rekee_nsa(ndime+3, ipoin) = press_hyd
!!$
!!$        !Water tracers
!!$        rekee_nsa(ndime+3 + 1, ipoin) = qv_k
!!$        rekee_nsa(ndime+3 + 2, ipoin) = 0.0_rp
!!$        rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp
!!$
!!$        if(kfl_ncons_nsa == 0) then
!!$           !
!!$           ! Momentum in Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0*rho_hyd
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0*rho_hyd
!!$        else
!!$           !
!!$           ! Momentum in NON-Conservative form
!!$           ! 
!!$           !u-velo
!!$           rekee_nsa(1, ipoin)       = u0
!!$           !w-velo
!!$           rekee_nsa(ndime, ipoin)   = w0
!!$
!!$        end if
!!$
!!$        !Perturbation fields: (total - reference)
!!$        bvess_nsa(1,ipoin,1)           = u0 !uvelo
!!$        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo
!!$
!!$        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
!!$        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k + dtheta
!!$
!!$        !Water tracers:
!!$        bvess_nsa(ndofn_nsa+1,ipoin,1)   = qv_k
!!$        bvess_nsa(ndofn_nsa+2,ipoin,1)   = 0.0_rp
!!$        bvess_nsa(ndofn_nsa+3,ipoin,1)   = 0.0_rp
!!$     end do !ipoin
!!$
!!$     !Return press and theta
!!$     outvar = 1
!!$
!!$     !
!!$     ! Build the mapping from 2D grid to column.
!!$     ! This is done here because it's only used by the Kessler module:
!!$     !
!!$     call nsa_build_column
!!$
!!$  end select
!!$
!!$  !
!!$  ! If RESTART_FILE ON (see *.nsa.dat), read the restart file to initialize
!!$  ! ONLY bvess_nsa(:,:). rekee_nsa are computed and stored above as usual.
!!$  !
!!$  if(kfl_rearst_nsa > 0) &
!!$       call nsa_rearst(outvar)
!!$  !
!!$  ! Compute sponge prameters:
!!$  !   
!!$  sponge: if( kfl_benme_nsa == 3 .or. kfl_benme_nsa == 4 .or. kfl_benme_nsa == 5 .or. &
!!$       kfl_benme_nsa == 5 .or. kfl_benme_nsa == 6 .or. kfl_benme_nsa == 7 .or. &
!!$       kfl_benme_nsa >= 200) then
!!$
!!$     !if(kfl_benme_nsa >= 200 .and. kfl_sponge_nsa < 0) & 
!!$     !     return
!!$
!!$     !HSlinear with the domain and sponge of the paper
!!$     xmin = xmin_nsa
!!$     xmax = xmax_nsa
!!$     ymin = ymin_nsa
!!$     ymax = ymax_nsa
!!$     zmin = zmin_nsa
!!$     zmax = zmax_nsa
!!$
!!$     !Sponge
!!$     dxs=dxs_nsa
!!$     dys=dys_nsa
!!$     dzs=dzs_nsa
!!$     ampx=ampx_nsa
!!$     ampy=ampy_nsa
!!$     ampz=ampz_nsa
!!$
!!$     if(kfl_botta_nsa > 0) then
!!$        xmin = -8000.0_rp
!!$        xmax =  8000.0_rp
!!$        zmin =     0.0_rp
!!$        zmax =  8000.0_rp
!!$
!!$        dxs=2000.0_rp
!!$        dzs=1000.0_rp
!!$        ampx=0.5_rp
!!$        ampz=1.0_rp
!!$     end if
!!$
!!$     if( kfl_benme_nsa == 205 .and. kfl_sptyp_nsa == 2) then     !NHLINEAR-Kessler     
!!$        xmin = -72000.0_rp
!!$        xmax =  72000.0_rp
!!$        zmin =      0.0_rp
!!$        zmax =  30000.0_rp
!!$
!!$        dxs=40000.0_rp
!!$        dzs=15000.0_rp
!!$        ampx=1.0_rp
!!$        ampz=0.5_rp
!!$
!!$     end if
!!$
!!$     if(kfl_sptyp_nsa == 1) then
!!$        !
!!$        ! Simple sponge
!!$        !
!!$        zs  = zmax - dzs
!!$        xsl = xmin + dxs
!!$        xsr = xmax - dxs
!!$
!!$        do ipoin= 1,npoin
!!$
!!$           cside=0.0_rp; ctop=0.0_rp
!!$           x=coord(1,ipoin)
!!$           z=coord(ndime,ipoin)
!!$
!!$           !Top Sponge
!!$           if (z >= zs) then
!!$              ctop=ampz*( (z - zs)/(zmax - zs) )**4
!!$           end if
!!$
!!$           !X-Lateral Sponge
!!$           if (x <= xsl) then
!!$              cside = ampx*( (x - xsl)/(xmin - xsl) )**4
!!$           else if (x >= xsr) then
!!$              cside = ampx*( (x - xsr)/(xmax - xsr) )**4
!!$           end if
!!$
!!$           ! boundary damping
!!$           !Store Values
!!$           bb(ipoin) = ctop + cside
!!$           bb(ipoin) = min(bb(ipoin), 1.0_rp)
!!$           aa(ipoin) = 1.0_rp - bb(ipoin)
!!$
!!$           bspon_nsa(1,ipoin) = aa(ipoin)
!!$           bspon_nsa(2,ipoin) = bb(ipoin)
!!$
!!$        end do !ipoin        
!!$
!!$     else
!!$        !
!!$        ! Lilly and Klemp 1978
!!$        !
!!$        nelx = nelx_nsa
!!$        nely = nelz_nsa
!!$
!!$        !zd    = 17700.0_rp   ! m
!!$        !zd    = 11500.0_rp   ! m
!!$        zd    = 13500.0_rp   ! m
!!$        alpha = 0.02_rp      ! s^-1
!!$
!!$        !Constants
!!$        ct = 0.5_rp
!!$        cs = 1.0_rp
!!$
!!$        !Constants
!!$        ztop   = zmax
!!$        xbound = dxs
!!$
!!$        xl     = xmin + xbound
!!$        xr     = xmax - xbound
!!$        !      zbound=0.5*ztop
!!$        !      dz=ztop - zbound
!!$        dsx    = 8000.0_rp
!!$        dsy    =  660.0_rp
!!$        zdcorrect = zd
!!$
!!$        !Initialize
!!$        aa=0; bb=0
!!$        do ipoin=1,npoin
!!$           ctop=0; cside=0
!!$           x=coord(1,ipoin)
!!$           z=coord(ndime,ipoin)
!!$
!!$           !
!!$           ! boundary damping
!!$           !
!!$           dbl   = min(x - xmin, xmax - x) ! distance from the boundary. xs in Restelli's thesis
!!$           beta  = (1.0_rp - tanh(dbl/dsx) )/tanh(dbl/dsx)
!!$           cside = cs*beta
!!$
!!$           !
!!$           ! top damping
!!$           ! first layer: damp lee waves
!!$           !
!!$           if(z .le. zdcorrect) then
!!$              ctop = 0
!!$           else
!!$              xid = (z-zd)/(ztop-zd) ! normalized coordinate
!!$              if(xid .lt. 0.5_rp) then
!!$                 abstaud = 0.5_rp*alpha*(1.0_rp - cos(xid*pi))
!!$
!!$              else
!!$                 abstaud = 0.5_rp*alpha*(1.0_rp + (xid - 0.5_rp)*pi)
!!$
!!$              endif
!!$
!!$              !ctop = ct * (dt*abstaud)/(1.0_rp + 0.5_rp*dt*abstaud)
!!$              ctop = ct*abstaud
!!$           endif
!!$
!!$           !
!!$           ! second layer: damp short waves
!!$           !
!!$           dbt = ztop - z ! distance from the boundary
!!$           beta = (1.0_rp - tanh(dbt/dsy))/tanh(dbt/dsy)
!!$
!!$           ctop = ctop + ct*beta
!!$
!!$           !Store Values
!!$           bb(ipoin) = min(ctop + cside,1.0_rp)
!!$           aa(ipoin) = 1.0_rp - bb(ipoin)
!!$
!!$           bspon_nsa(1,ipoin) = aa(ipoin)
!!$           bspon_nsa(2,ipoin) = bb(ipoin)
!!$
!!$        end do !i
!!$     end if
!!$
!!$  end if sponge
!!$  
!!$end subroutine nsa_initial_conditions
!!$
!!$subroutine nsa_loadref(q,q_ref)
!!$
!!$  use def_master
!!$  use def_domain
!!$  use def_nastal
!!$  use def_parame
!!$
!!$  implicit none
!!$  integer(ip) :: ipoin, i,j,k
!!$  real(rp)    :: q(ndofn_nsa+3,npoin), q_ref(ndofn_nsa+3,npoin)
!!$  real(rp)    :: qaux(ndofn_nsa+3,nelz_nsa), qaux_ref(ndofn_nsa+3,nelz_nsa)
!!$
!!$  integer(ip) :: icase_npoin
!!$
!!$  !flag to define if the case200ref.out has npoin lines, or
!!$  !nlevels lines.
!!$  icase_npoin = 1!
!!$
!!$  q_ref = 0.0
!!$  open(1,file='case200ref.out')
!!$
!!$  !
!!$  ! The case200 file is stored in npoin rows (obtained from the matlab routines):
!!$  ! 
!!$  if(ndime < 3) then
!!$     !
!!$     ! 2D
!!$     !
!!$     do ipoin=1,npoin
!!$
!!$        !do ipoin = 1,npoin
!!$        !          densi             theta               qv         u-velo (almost sure)
!!$        read(1,*) q_ref(1,ipoin), q_ref(4,ipoin), q_ref(5,ipoin), q(2,ipoin) 
!!$     end do
!!$  else
!!$     !
!!$     ! 3D
!!$     !
!!$     do ipoin=1,npoin
!!$
!!$        !do ipoin = 1,npoin
!!$        !          densi             theta            qv          u-velo (almost sure)
!!$        read(1,*) q_ref(1,ipoin), q_ref(5,ipoin), q_ref(6,ipoin), q(2,ipoin) 
!!$        !print*,'', q_ref(1,ipoin), q_ref(5,ipoin), q_ref(6,ipoin), q(2,ipoin) 
!!$     end do
!!$  end if
!!$
!!$  close(1)
!!$
!!$end subroutine nsa_loadref
!!$
!!$subroutine nsa_loadref_field(q,q_ref)
!!$
!!$  use def_master
!!$  use def_domain
!!$  use def_nastal
!!$  use def_parame
!!$
!!$  implicit none
!!$  integer(ip) :: ipoin, i,j,k
!!$  real(rp)    :: q(ndofn_nsa+3,npoin), q_ref(ndofn_nsa+3,npoin)
!!$  real(rp)    :: qaux(ndofn_nsa+3,nelz_nsa), qaux_ref(ndofn_nsa+3,nelz_nsa)
!!$
!!$  integer(ip) :: icase_npoin
!!$
!!$  !flag to define if the case200ref.out has npoin lines, or
!!$  !nlevels lines.
!!$  icase_npoin = 1!
!!$
!!$  q_ref = 0.0
!!$  
!!$  ! Read the field from file passed through *.dom.dat as
!!$  ! FIELDA, END_FIELDS
!!$  !
!!$  do ipoin=1,npoin
!!$     densi_ref_nsa => xfiel(-nfiel_nsa(1))%a
!!$     tempe_ref_nsa => xfiel(-nfiel_nsa(2))%a
!!$     qvapo_ref_nsa => xfiel(-nfiel_nsa(3))%a
!!$     uvelo_ref_nsa => xfiel(-nfiel_nsa(4))%a
!!$  end do
!!$  
!!$end subroutine nsa_loadref_field

subroutine nsa_initial_conditions(outvar)
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_inimet
  ! NAME 
  !    nsa_inimet
  ! DESCRIPTION
  !    This routine sets up the initial condition for the meteo benchmarks
  ! USED BY
  !    nsa_iniunk
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_parame
  use      def_nastal
  implicit none

  integer(ip) :: idime,icomp,idofn,ipoin,itime,ibubb,nbubb,kshbu,ix,iy,iz
  real(rp)    :: velmi,xfact,xfac2,rgacp,xinve,xradi,xdifi,xrano(3),xtemp,rauxi,exner
  real(rp)    :: u0,v0,w0,u_k,v_k,w_k,thetac,theta0,press0,rho0,dtheta,dtracer,tracerc,tracer0,tracer_k,p00
  real(rp)    :: dqv,dqc,dqr
  real(rp)    :: rc,rc0,rc1,rc2,rc3,xc,xc2,yc,zc,r,xr,yr,zr,zl,xr_tr,yr_tr,zr_tr,sigma
  real(rp)    :: x,y,z, xmin,xmax,ymin,ymax,zmin,zmax
  real(rp)    :: theta_k,pi_k,rho_k,press_k, theta_ref,theta_z0,dtdz
  real(rp)    :: theta_hyd,pi_hyd,rho_hyd,press_hyd, rho_baro
  real(rp)    :: a,b,c,d,e,l,lz,hc,ac,l2,n,n2,c2
  real(rp)    :: rho,theta,pressure, qv,bv,bv2,g2,static_stability,es,gam,pii
  real(rp)    :: xlv, ep2,svp1,svp2,svp3,svpt0,rhowater,qvs,temp,RH_k,qv_k

  character   :: fname*72

  !
  ! outvar = 1 --> densi and theta
  ! outvar = 2 --> press and theta
  !
  integer(ip) :: outvar
  real(rp)    :: q(ndofn_nsa+3,npoin), q_ref(ndofn_nsa+3,npoin), q_exact(ndofn_nsa,npoin)
  !  real(rp)    :: q(7,npoin), q_ref(7,npoin), q_exact(4,npoin)
  real(rp)    :: tmpt(nelz_nsa+1)

  !
  ! Sponge parameters:
  !
!!  real(rp)   :: aa(npoin),bb(npoin)
  real(rp)   :: xsr,xsl,ysr,ysl,zs,dxs,dys,dzs
  real(rp)   :: amp,ampx,ampy,ampz,ctop,cside

  !
  !Exact solution for mountains
  !
  real(rp)   ::  rho00, p_k, rhofac, afac, hafac, delta_k, temp0
  real(rp)   ::  ddelta_dx_k, ddelta_dy_k, ddelta_dz_k, xappa
  real(rp)   ::  zd,alpha,ct,cs,xbound,xl,dsx,dsy,zdcorrect,ztop
  real(rp)   ::  dbl,beta,xid,abstaud,dt,dbt
  integer(ip)::  nelx,nely,nelz
  real(rp)   ::  rho_exact, pi_exact, u_exact, v_exact, w_exact, theta_exact
  real(rp)   ::  rho_pert, p_pert, t_pert, theta_pert
  real(rp)   ::  m, intercept

  !
  !  input values for squall-line test
  !
  xlv      = 2500000.0_rp   !Latent heat of vaporization
  ep2      = 0.6217504_rp   !Rd/Rv
  svp1     = 0.6112000_rp   
  svp2     = 17.67000_rp
  svp3     = 29.65000_rp
  svpt0    = 273.1500_rp
  rhowater = 1000.000_rp    !density of rainwater
  
  !Initialize what will be used later (it is necessary for parallel)
  u0 = 0.0_rp
  v0 = 0.0_rp
  w0 = 0.0_rp

  if (nzone > 1) call runend("NSA_INITIAL_CONDITIONS: THIS SUB IS NOT PREPARED TO RUN WITH ZONES.")

  select case(kfl_benme_nsa)

  case (15) !Inertia-Gravity Wave
     !Case Constants

     temp0  = 300.0_rp
     theta0 = tempe_nsa

     xmin = 0.0_rp
     xmax = 300000.0_rp
     zmin = 0.0_rp
     zmax = 10000.0_rp

     c = cvcoe_nsa/rgasc_nsa
     c2 = rgasc_nsa/cpcoe_nsa

     u0 =20.0_rp
     w0 = 0.0_rp

     bv=0.01_rp
     bv2=bv*bv

     g2=grnor_nsa*grnor_nsa

     ac = 5000.0_rp
     xl = 300000.0_rp
     zl = 10000.0_rp
     xc = xl/3.0_rp
     zc = zl/2.0_rp
     thetac = 0.01_rp

     do ipoin = 1,npoin
        x = coord(1,ipoin)
        z = coord(ndime,ipoin)

        u_k = u0
        w_k = w0

        pi_hyd =g2/(cpcoe_nsa*temp0*bv2)*(exp(- bv2/grnor_nsa*z ) - 1.0_rp) + 1.0_rp
        theta_hyd = temp0*exp( bv2/grnor_nsa*z )
        dtheta    = thetac*sin( pi*z/zl )/( 1.0_rp + ( (x-xc)/ac )**2 )

        theta_k   = theta_hyd + dtheta
        rho_k     = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_hyd)**c
        rho_hyd   = pbaro_nsa/(rgasc_nsa*theta_hyd)*(pi_hyd)**c
        call nsa_stalaw(2,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)
        call nsa_stalaw(2,0,rho_hyd,press_hyd,theta_hyd,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Press_hyd
        rekee_nsa(ndime+1, ipoin) = rho_hyd
        !Theta_hyd
        rekee_nsa(ndime+2, ipoin) = theta_hyd
        !Densi_hyd
        rekee_nsa(ndime+3, ipoin) = press_hyd

        !Perturbation fields: (total - reference)
        bvess_nsa(1,ipoin,1)           = u0 !uvelo
        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo

        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

     end do !ip

     !Return density and theta
     outvar = 1

     case(1)
     !
     !Rising Thermal Bubble
     !
     !Case Constants
     c=cvcoe_nsa/rgasc_nsa
     theta0=tempe_nsa

     xmin = 0.0_rp
     xmax = 1000.0_rp
     zmin = 0.0_rp
     zmax = 1000.0_rp
     !zmax = 1500.0_rp
     if(ndime > 2) then
        ymin = ymin_nsa
        ymax = ymax_nsa
     end if

     u0=0.0_rp
     w0=0.0_rp
     xc=0.5_rp*(xmin + xmax)
     zc=350.0_rp
     if(ndime > 2) &
          yc = 0.5_rp*(ymin + ymax)
     rc=250.0_rp

     !
     ! Perturbation intensity:
     !
     thetac=0.5_rp   ! With bubble

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(ndime,ipoin)
        if(ndime > 2) &
             y=coord(ndime-1,ipoin)


        dtheta=0.0_rp
        if(ndime < 3) then
           !
           ! 2D
           !
           r=sqrt( (x-xc)**2.0_rp + (z-zc)**2.0_rp )
           if (r < rc) then
              dtheta = 0.5_rp*thetac*(1.0_rp + cos(pi*r/rc) )
           end if
        else
           !
           ! 3D
           !
           if( kfl_cylind_nsa > 0 ) then
              !Cylinder
              r=sqrt( (x-xc)**2.0_rp + (z-zc)**2.0_rp )
           else
              !Sphere
              r=sqrt( (x-xc)**2.0_rp + (y-yc)**2_rp + (z-zc)**2.0_rp ) !Sphere
           end if

           if (r < rc) then
              dtheta = 0.5_rp*thetac*(1.0_rp + cos(pi*r/rc) )
           end if
        end if


        !Total fields:
        theta_k = theta0 + dtheta
        pi_k    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta0)*z
        rho_k   = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_k)**c
        call nsa_stalaw(2,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Reference Fields (hydrostatic equilibrium):
        theta_hyd = tempe_nsa
        pi_hyd    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta_hyd)*z
        rho_hyd   = pbaro_nsa/(rgasc_nsa*theta_hyd)*(pi_hyd)**c
        call nsa_stalaw(2,0,rho_hyd,press_hyd,theta_hyd,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Press_hyd
        rekee_nsa(ndime+1, ipoin) = rho_hyd
        !Theta_hyd
        rekee_nsa(ndime+2, ipoin) = theta_hyd
        !Densi_hyd
        rekee_nsa(ndime+3, ipoin) = press_hyd

        !Perturbation fields: (total - reference)
        bvess_nsa(1,ipoin,1)           = u0 !uvelo
        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo

        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

     end do !ipoin

     !Return density and theta
     outvar = 1

     case(11)
     !
     !Rising Thermal Bubble of Robert 1993
     !
     !Case Constants
     c=cvcoe_nsa/rgasc_nsa
     theta0=tempe_nsa

     xmin = 0.0_rp
     xmax = 1000.0_rp
     zmin = 0.0_rp
     zmax = 1000.0_rp
     !zmax = 1500.0_rp

     u0 = 0.0_rp
     w0 = 0.0_rp
     xc = 500.0_rp
     zc = 260.0_rp
     rc = 50.0_rp
     sigma = 100.0_rp
     sigma = sigma*sigma

     !Perturbation intensity:
     thetac=0.5_rp   ! With bubble

     do ipoin=1,npoin
        x=coord(1,ipoin)
        y=coord(ndime-1,ipoin)
        z=coord(ndime,ipoin)

        if(ndime < 3) then
           r=sqrt( (x-xc)**2.0_rp + (z-zc)**2.0_rp )
           if (r > rc) then
              dtheta = thetac*exp(-(r-rc)*(r-rc)/sigma)
           else
              dtheta = thetac
           end if
        else
           !Cosine function as Kelly and Giraldo 2012 JCP 3D
           xc = 500.0_rp
           yc = 500.0_rp
           zc = 260.0_rp
           xr = 250.0_rp
           yr = 250.0_rp
           zr = 250.0_rp
           rc =   1.0_rp

           dtheta = 0.0_rp
           r=sqrt( ((x-xc)/xr)**2_rp  + ((y-yc)/yr)**2_rp + ((z-zc)/zr)**2_rp )
           if (r < rc) then
              dtheta=thetac*(1.0_rp + cos(pi*r/rc))
           end if
        end if

        !Total fields:
        theta_k = theta0 + dtheta
        pi_k    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta0)*z
        rho_k   = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_k)**c
        call nsa_stalaw(2,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Reference Fields (hydrostatic equilibrium):
        theta_hyd = tempe_nsa
        pi_hyd    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta_hyd)*z
        rho_hyd   = pbaro_nsa/(rgasc_nsa*theta_hyd)*(pi_hyd)**c
        call nsa_stalaw(2,0,rho_hyd,press_hyd,theta_hyd,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Press_hyd
        rekee_nsa(ndime+1, ipoin) = rho_hyd
        !Theta_hyd
        rekee_nsa(ndime+2, ipoin) = theta_hyd
        !Densi_hyd
        rekee_nsa(ndime+3, ipoin) = press_hyd

        !Perturbation fields: (total - reference)
        bvess_nsa(1,ipoin,1)           = u0 !uvelo
        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo

        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

     end do !ipoin

     !Return density and theta
     outvar = 1

     case(12)
     !
     !Rising Thermal Bubble as in Ahmad
     !
     !Case Constants
     c=cvcoe_nsa/rgasc_nsa
     theta0 = tempe_nsa

     xmin = xmin_nsa
     xmax = xmax_nsa
     zmin = 0.0_rp
     zmax = zmax_nsa
     if(ndime > 2) then
        ymin = ymin_nsa
        ymax = ymax_nsa
     end if

     u0=0.0_rp
     w0=0.0_rp

     xc = 0.5_rp*(xmin + xmax)
     zc = zc_nsa
     if(ndime > 2) &
          yc = 0.5_rp*(ymin + ymax)

     rc = 2000.0_rp

     xr = xradi_nsa
     yr = yradi_nsa
     zr = zradi_nsa
     if(xradi_nsa < 0) xr = 2000.0_rp
     if(yradi_nsa < 0) yr = 2000.0_rp
     if(zradi_nsa < 0) zr = 2000.0_rp

     !Perturbation intensity:
     if(thetac_nsa < -999.0_rp) then
        thetac=0.0_rp       ! With bubble
     else
        thetac=thetac_nsa   ! With bubble
     end if

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(ndime,ipoin)
        if(ndime > 2) &
             y=coord(ndime-1,ipoin)

        dtheta=0.0_rp
        if(ndime < 3) then
           r=sqrt( (x-xc)**2.0_rp + (z-zc)**2.0_rp )   
           if (r <= rc) then
              dtheta=thetac*(1.0_rp - r)
           end if
        else
           !r=sqrt( ((x-xc)/xr)**2_rp  + ((y-yc)/yr)**2_rp + ((z-zc)/zr)**2_rp )
           r=sqrt( ((x-xc)/xr)**2.0_rp + ((z-zc)/zr)**2.0_rp )   
           if (r <= 1.0_rp) then
              dtheta=thetac*(1.0_rp - r)
           end if
        end if

        !Total fields:
        theta_k = theta0 + dtheta
        pi_k    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta0)*z
        rho_k   = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_k)**c
        call nsa_stalaw(2,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Reference Fields (hydrostatic equilibrium):
        theta_hyd = tempe_nsa
        pi_hyd    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta_hyd)*z
        rho_hyd   = pbaro_nsa/(rgasc_nsa*theta_hyd)*(pi_hyd)**c
        call nsa_stalaw(2,0,rho_hyd,press_hyd,theta_hyd,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Press_hyd
        rekee_nsa(ndime+1, ipoin) = rho_hyd
        !Theta_hyd
        rekee_nsa(ndime+2, ipoin) = theta_hyd
        !Densi_hyd
        rekee_nsa(ndime+3, ipoin) = press_hyd

        !Perturbation fields: (total - reference)
        bvess_nsa(1,ipoin,1)           = u0 !uvelo
        bvess_nsa(ndime,ipoin,1)       = w0 !wvelo
        if(ndime > 2) &
             bvess_nsa(ndime-1,ipoin,1)       = v0 !vvelo

        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

     end do !ipoin

     !Return density and theta
     outvar = 1

     case(2)
     !
     ! Density current
     !
     c=cvcoe_nsa/rgasc_nsa
     theta0=tempe_nsa
     u0=0.0_rp
     w0=0.0_rp

     xmin = -25600.0_rp
     xmax =  25600.0_rp
     zmin =  0.0_rp
     zmax =  6400.0_rp

     xc=0.5_rp*(xmin + xmax)
     zc=3000.0_rp
     xr=4000.0_rp
     zr=2000.0_rp
     rc=1.0_rp

     !Perturbation intensity:
     thetac = -15.0_rp

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(ndime,ipoin)

        r=sqrt( ((x-xc)/xr)**2 + ((z-zc)/zr)**2 )
        dtheta=0
        if (r <= rc) then
           dtheta=thetac/2*(1 + cos(pi*r) )
        end if
        !Total fields:
        theta_k = theta0 + dtheta
        pi_k    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta0)*z
        rho_k   = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_k)**c
        call nsa_stalaw(2,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Reference Fields (hydrostatic equilibrium):
        theta_hyd = tempe_nsa
        pi_hyd    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta_hyd)*z
        rho_hyd   = pbaro_nsa/(rgasc_nsa*theta_hyd)*(pi_hyd)**c
        call nsa_stalaw(2,0,rho_hyd,press_hyd,theta_hyd,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Press_hyd
        rekee_nsa(ndime+1, ipoin) = rho_hyd
        !Theta_hyd
        rekee_nsa(ndime+2, ipoin) = theta_hyd
        !Densi_hyd
        rekee_nsa(ndime+3, ipoin) = press_hyd

        !Perturbation fields: (total - reference)
        bvess_nsa(1,ipoin,1)           = u0 !uvelo
        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo

        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

     end do !ipoin

     !Return density and theta
     outvar = 1

     case(3)
     !
     ! HS Linear mountain
     !
     !Case Constants
     theta0    = tempe_nsa
     temp0     = tempe_nsa
     rho_baro  = 1.0_rp/(rgasc_nsa*theta0)*pbaro_nsa; ! surface density
     rho00     = 1.0_rp/(rgasc_nsa*theta0)*pbaro_nsa; ! surface density

     u0=20.0_rp
     w0=0.0_rp

     ac = 10000.0_rp
     hc = 1.0_rp !Small Mountain
     xc = 240000.0_rp/2.0_rp

     if(kfl_botta_nsa > 0) then
        !
        !Botta and Klein
        !
        u0 = 20.0_rp
        ac = 800.0_rp
        hc = 2000.0_rp
        xc = 0.0_rp
     end if

     n2 = grnor_nsa**2/(cpcoe_nsa*temp0); ! 1/s^2
     l2 = n2/u0**2;              ! 1/m^2
     l = sqrt(l2);               ! 1/m
     xappa=0.286_rp

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(ndime,ipoin)

        !Total variables
        pi_k    = exp(-grnor_nsa/(cpcoe_nsa*theta0)*z);
        theta_k = theta0/pi_k;
        press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);
        call nsa_stalaw(1,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Reference Solution1
        rho_hyd   = rho_k
        press_hyd = press_k
        theta_hyd = theta_k

        !Press_hyd
        rekee_nsa(ndime+1, ipoin) = rho_hyd
        !Theta_hyd
        rekee_nsa(ndime+2, ipoin) = theta_hyd
        !Densi_hyd
        rekee_nsa(ndime+3, ipoin) = press_hyd

        if(kfl_ncons_nsa == 0) then
           !
           ! Momentum in Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = u0*rho_hyd
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0*rho_hyd
        else
           !
           ! Momentum in NON-Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = u0
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0

        end if

        !Perturbation fields: (total - reference)
        bvess_nsa(1,ipoin,1)           = u0 !uvelo
        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo

        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

     end do !ipoin

     !Return press and theta
     outvar = 1

     case (4) !Linear Hydrostatic: Large Mountain

     if(ndime < 3) then
        !
        ! 2D:
        !
        
        !Case Constants
        theta0    = tempe_nsa
        rho_baro  = 1.0_rp/(rgasc_nsa*theta0)*pbaro_nsa; ! surface density
        u0 = 20.0_rp
        w0 =  0.0_rp

        do ipoin=1,npoin
           x=coord(1,ipoin)
           z=coord(ndime,ipoin)

           !Total variables
           pi_k    = exp(-grnor_nsa/(cpcoe_nsa*theta0)*z);
           theta_k = theta0/pi_k;
           press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);
           call nsa_stalaw(1,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

           !Reference Solution1
           rho_hyd   = rho_k
           press_hyd = press_k
           theta_hyd = theta_k

           !Press_hyd
           rekee_nsa(ndime+1, ipoin) = rho_hyd
           !Theta_hyd
           rekee_nsa(ndime+2, ipoin) = theta_hyd
           !Densi_hyd
           rekee_nsa(ndime+3, ipoin) = press_hyd

           if(kfl_ncons_nsa == 0) then
              !
              ! Momentum in Conservative form
              ! 
              !u-velo
              rekee_nsa(1, ipoin)       = u0*rho_hyd
              !w-velo
              rekee_nsa(ndime, ipoin)   = w0*rho_hyd
           else
              !
              ! Momentum in NON-Conservative form
              ! 
              !u-velo
              rekee_nsa(1, ipoin)       = u0
              !w-velo
              rekee_nsa(ndime, ipoin)   = w0

           end if

           !Perturbation fields: (total - reference)
           bvess_nsa(1,ipoin,1)           = u0 !uvelo
           bvess_nsa(ndime,ipoin,1)       = w0 !vvelo

           bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
           bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
        end do

        !Return press and theta
        outvar = 1

     else
        !
        ! 3D
        !

        !Case Constants
        theta0    = tempe_nsa
        rho_baro  = 1.0_rp/(rgasc_nsa*theta0)*pbaro_nsa; ! surface density
        u0 = 20.0_rp
        v0 =  0.0_rp
        w0 =  0.0_rp

        do ipoin=1,npoin
           x=coord(1,ipoin)
           y=coord(ndime-1,ipoin)
           z=coord(ndime,ipoin)

           !Total variables
           pi_k    = exp(-grnor_nsa/(cpcoe_nsa*theta0)*z);
           theta_k = theta0/pi_k;
           press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);
           call nsa_stalaw(1,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

           !Reference Solution1
           rho_hyd   = rho_k
           press_hyd = press_k
           theta_hyd = theta_k

           !Press_hyd
           rekee_nsa(ndime+1, ipoin) = rho_hyd
           !Theta_hyd
           rekee_nsa(ndime+2, ipoin) = theta_hyd
           !Densi_hyd
           rekee_nsa(ndime+3, ipoin) = press_hyd

           if(kfl_ncons_nsa == 0) then
              !
              ! Momentum in Conservative form
              ! 
              !u-velo
              rekee_nsa(1, ipoin)       = u0*rho_hyd
              !v-velo
              rekee_nsa(ndime-1, ipoin) = v0*rho_hyd
              !w-velo
              rekee_nsa(ndime, ipoin)   = w0*rho_hyd
           else
              !
              ! Momentum in NON-Conservative form
              ! 
              !u-velo
              rekee_nsa(1, ipoin)       = u0
              !v-velo
              !rekee_nsa(ndime-1, ipoin) = v0
              !w-velo
              rekee_nsa(ndime, ipoin)   = w0

           end if

           !Perturbation fields: (total - reference)
           bvess_nsa(1,ipoin,1)           = u0 !uvelo
           bvess_nsa(ndime-1,ipoin,1)     = v0 !vvelo
           bvess_nsa(ndime,ipoin,1)       = w0 !wvelo

           bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
           bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k
        end do

        !Return press and theta
        outvar = 1

     end if

     case(5)
     !
     ! NH Linear mountain
     !
     !Case Constants
     temp0 = tempe_nsa
     theta0    = tempe_nsa
     c     = cvcoe_nsa/rgasc_nsa
     c2    = rgasc_nsa/cpcoe_nsa

     u0    =  10.0_rp
     w0    =  0.0_rp
     bv    =  0.01_rp
     bv2   = bv*bv
     g2    = grnor_nsa*grnor_nsa

     ac = 1000.0_rp
     hc = 1.0_rp !Small Mountain
     xc = 0.0_rp

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(ndime,ipoin)

        !Total variables
        pi_k    = g2/(cpcoe_nsa*theta0*bv2)*(exp(- bv2/grnor_nsa*z ) - 1.0_rp) + 1.0_rp
        theta_k = theta0*exp( bv2/grnor_nsa*z )
        rho_k   = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_k)**c
        call nsa_stalaw(2,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)
        !press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);

        !Reference Solution1
        rho_hyd   = rho_k
        press_hyd = press_k
        theta_hyd = theta_k

        !Press_hyd
        rekee_nsa(ndime+1, ipoin) = rho_hyd
        !Theta_hyd
        rekee_nsa(ndime+2, ipoin) = theta_hyd
        !Densi_hyd
        rekee_nsa(ndime+3, ipoin) = press_hyd

        if(kfl_ncons_nsa == 0) then
           !
           ! Momentum in Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = u0*rho_hyd
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0*rho_hyd
        else
           !
           ! Momentum in NON-Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = u0
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0

        end if

        !Perturbation fields: (total - reference)
        bvess_nsa(1,ipoin,1)           = u0 !uvelo
        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo

        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

     end do !ipoin

     !Return press and theta
     outvar = 1

     case(200)
     !
     ! Moisture: this case is similar to case 200
     ! except that it already has a supersaturaded cloud 
     ! in the initial field. qv initial is the same as for case 200.
     !
     !
     u0 = 0.0_rp
     v0 = 0.0_rp
     w0 = 0.0_rp
     
     if(ndime < 3) then
        !
        ! 2D:
        !
        nvar_nsa = 7_ip
        c=cvcoe_nsa/rgasc_nsa
        c2=rgasc_nsa/cpcoe_nsa

        !Center of thermal:
        xc = 0.5_rp*(xmin_nsa + xmax_nsa)
        zc = 2000.0_rp

        !Reference theta and perturbation theta
        theta0 = 300.0_rp
        if(thetac_nsa > -999.0_rp) then
           thetac = thetac_nsa
        else
           thetac = 3.0_rp
        end if

        !Radii of thermal and tracer
        rc     = 1.0_rp

        !Thermal
        xr     = 10000.0_rp
        zr     =  1500.0_rp    

        !Tracer
        xr_tr  = 2500.0_rp
        zr_tr  = 1000.0_rp

        !Load Reference Values
        !call nsa_loadref_field

        do ipoin=1,npoin
           x=coord(1,ipoin)
           z=coord(ndime,ipoin)

           !Background density and theta (from sounding):
           rho0   = xfiel(1)%a(1,ipoin,1)!densiref_nsa(ipoin)
           theta0 = xfiel(2)%a(1,ipoin,1) !temperef_nsa(ipoin)
           press0 = pbaro_nsa*(rho0*rgasc_nsa*theta0/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)

           !Cylinder:
           r=sqrt( ((x-xc)/xr)**2_rp + ((z-zc)/zr)**2_rp )

           dtheta = 0.0_rp
           qv     = 0.0_rp
           if (r < rc) then
              dtheta = thetac*(cos(pi*r/2.0_rp))**2_rp
           end if

           theta_k = theta0 + dtheta
           press_k = pbaro_nsa*(rho0*rgasc_nsa*theta_k/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)
           rho_k   = pbaro_nsa*(press_k/pbaro_nsa)**(cvcoe_nsa/cpcoe_nsa)/(rgasc_nsa*theta_k)

           !
           ! Total and perturbation Fields (for plotting only):
           !
           q_nsa(1,ipoin) = rho_k
           if( kfl_uniformvelocity_nsa >= 0) then
              q_nsa(2,ipoin) =  xfiel(4)%a(1,ipoin,1) !if commented, horizontal velocity comes from sounding
           else
              !Read from file:
              q_nsa(2,ipoin) = xfiel(4)%a(1,ipoin,1)
           end if
           q_nsa(3,ipoin) = 0.0_rp  !wvelo
           q_nsa(4,ipoin) = theta_k
           q_nsa(5,ipoin) = xfiel(3)%a(1,ipoin,1)
           q_nsa(6,ipoin) = 0.0_rp
           q_nsa(7,ipoin) = 0.0_rp

           !Densi_hyd
           rekee_nsa(ndime+1, ipoin) = rho0
           !Theta_hyd
           rekee_nsa(ndime+2, ipoin) = theta0
           !Press_hyd
           rekee_nsa(ndime+3, ipoin) = press0

           !Water tracers
           rekee_nsa(ndime+3 + 1, ipoin) = xfiel(3)%a(1,ipoin,1)
           rekee_nsa(ndime+3 + 2, ipoin) = 0.0_rp
           rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp

           if(kfl_ncons_nsa == 0) then
              !
              ! Momentum in Conservative form
              ! 
              !u-velo
              rekee_nsa(1, ipoin)       = q_nsa(2,ipoin)*rho0
              !w-velo
              rekee_nsa(ndime, ipoin)   = w0*rho0
           else
              !
              ! Momentum in NON-Conservative form
              ! 
              !u-velo
              rekee_nsa(1, ipoin)       = q_nsa(2,ipoin)
              !w-velo
              rekee_nsa(ndime, ipoin)   = w0

           end if

           !Velocity:
           bvess_nsa(1,ipoin,1)           = q_nsa(2,ipoin) !uvelo
           bvess_nsa(ndime,ipoin,1)       = w0         !vvelo

           !Total fields:
           bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
           bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

           !Water tracers:
           bvess_nsa(ndofn_nsa+1,ipoin,1)   = q_nsa(5,ipoin)
           bvess_nsa(ndofn_nsa+2,ipoin,1)   = q_nsa(6,ipoin)
           bvess_nsa(ndofn_nsa+3,ipoin,1)   = q_nsa(7,ipoin)

        end do !ipoin

        !Return press and theta
        outvar = 1

        !
        ! Build the mapping from 2D grid to column.
        ! This is done here because it's only used by the Kessler module:
        !
        call nsa_build_column

     else
        !
        ! 3D:
        !
        nvar_nsa = 8_ip
        c=cvcoe_nsa/rgasc_nsa
        c2=rgasc_nsa/cpcoe_nsa

        !Center of thermal:
        xc = 0.5_rp*(xmin_nsa + xmax_nsa)
        yc = 0.5_rp*(ymin_nsa + ymax_nsa)
        zc = 2000.0_rp

        !Reference theta and perturbation theta
        theta0 = 300.0_rp
        if(thetac_nsa > -999.0_rp) then
           thetac = thetac_nsa
        else
           thetac = 3.0_rp
        end if
        
        !Radii of thermal and tracer
        rc     = 1.0_rp

        !Thermal
        xr     = 10000.0_rp
        yr     = 10000.0_rp
        !zr     =  1500.0_rp
        zr     =  2000.0_rp

        !Tracer
        xr_tr  = 2500.0_rp
        yr_tr  = 2500.0_rp
        zr_tr  = 1000.0_rp
        
        do ipoin=1,npoin
           x=coord(1,ipoin)
           y=coord(ndime-1,ipoin)
           z=coord(ndime,ipoin)
           
           !Background density and theta (from sounding):
           rho0   = xfiel(1)%a(1,ipoin,1) !densiref_nsa(ipoin)
           theta0 = xfiel(2)%a(1,ipoin,1) !temperef_nsa(ipoin)
           press0 = pbaro_nsa*(rho0*rgasc_nsa*theta0/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)
           if( kfl_cylind_nsa > 0) then
              !Cylinder:
              r=sqrt( ((x-xc)/xr)**2_rp + ((z-zc)/zr)**2_rp )
           else
              !Sphere:
              r=sqrt( ((x-xc)/xr)**2_rp  + ((y-yc)/yr)**2_rp + ((z-zc)/zr)**2_rp ) !sphere
           end if

           dtheta = 0.0_rp
           qv     = 0.0_rp
           if( kfl_cylind_nsa > 0) then
              if (r < rc .and. y <= yc + yr .and. y>= yc - yr) then
                 dtheta = thetac*(cos(pi*r/2.0_rp))**2_rp
              end if
           else
              if (r < rc) then
                 dtheta = thetac*(cos(pi*r/2.0_rp))**2_rp
              end if
           end if

           theta_k = theta0 + dtheta
           press_k = pbaro_nsa*(rho0*rgasc_nsa*theta_k/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)

           rho_k   = pbaro_nsa*(press_k/pbaro_nsa)**(cvcoe_nsa/cpcoe_nsa)/(rgasc_nsa*theta_k)
           !
           ! Total and perturbation Fields (for plotting only):
           !
           q_nsa(1,ipoin) = rho_k
           if( kfl_uniformvelocity_nsa >= 0) then
              q_nsa(2,ipoin) =  uvelo_nsa !if commented, horizontal velocity comes from sounding
           else
              !Read from file:
              q_nsa(2,ipoin) = xfiel(4)%a(1,ipoin,1) !uveloref_nsa(ipoin)
           end if

           q_nsa(3,ipoin) = 0.0_rp  !vvelo
           q_nsa(4,ipoin) = 0.0_rp  !wvelo
           q_nsa(5,ipoin) = theta_k
           q_nsa(6,ipoin) = xfiel(3)%a(1,ipoin,1) !qvaporef_nsa(ipoin)
           q_nsa(7,ipoin) = 0.0_rp
           q_nsa(8,ipoin) = 0.0_rp

           !Densi_hyd
           rekee_nsa(ndime+1, ipoin) = rho0
           !Theta_hyd
           rekee_nsa(ndime+2, ipoin) = theta0
           !Press_hyd
           rekee_nsa(ndime+3, ipoin) = press0

           !Water tracers
           rekee_nsa(ndime+3 + 1, ipoin) = xfiel(3)%a(1,ipoin,1) !qvaporef_nsa(ipoin)
           rekee_nsa(ndime+3 + 2, ipoin) = 0.0_rp
           rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp

           if(kfl_ncons_nsa == 0) then
              !
              ! Momentum in Conservative form
              ! 
              !u-velo
              rekee_nsa(1, ipoin)       = q_nsa(2,ipoin)*rho0
              !w-velo
              rekee_nsa(ndime, ipoin)   = w0*rho0
           else
              !
              ! Momentum in NON-Conservative form
              ! 
              !u-velo
              rekee_nsa(1, ipoin)       = q_nsa(2,ipoin)
              !w-velo
              rekee_nsa(ndime, ipoin)   = w0

           end if

           !Velocity:
           bvess_nsa(1,ipoin,1)           = q_nsa(2,ipoin) !uvelo
           bvess_nsa(ndime,ipoin,1)       = w0         !vvelo

           !Total fields:
           bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
           bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

           !Water tracers:
           bvess_nsa(ndofn_nsa+1,ipoin,1)   = q_nsa(6,ipoin)
           bvess_nsa(ndofn_nsa+2,ipoin,1)   = q_nsa(7,ipoin)
           bvess_nsa(ndofn_nsa+3,ipoin,1)   = q_nsa(8,ipoin)

        end do !ipoin

        !Return press and theta
        outvar = 1

        !
        ! Build the mapping from 2D grid to column.
        ! This is done here because it's only used by the Kessler module:
        !
        call nsa_build_column

     end if

     case (204) !Dynamics as in LH mountain (large mountain), but also with moisture

     !
     ! HS Linear mountain with tracer (nvar_nsa = 7)
     !
     !Case Constants
     nvar_nsa  = 7_ip
     theta0    = tempe_nsa
     temp0     = tempe_nsa
     rho_baro  = 1.0_rp/(rgasc_nsa*theta0)*pbaro_nsa; ! surface density
     rho00     = 1.0_rp/(rgasc_nsa*theta0)*pbaro_nsa; ! surface density

     u0=20.0_rp
     w0=0.0_rp

     ac = 10000.0_rp
     hc = 1.0_rp !Small Mountain
     xc = 240000.0_rp/2.0_rp

     n2 = grnor_nsa**2/(cpcoe_nsa*temp0); ! 1/s^2
     l2 = n2/u0**2;                       ! 1/m^2
     l = sqrt(l2);                        ! 1/m
     xappa=0.286_rp

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(ndime,ipoin)

        !Total variables
        pi_k    = exp(-grnor_nsa/(cpcoe_nsa*theta0)*z);
        theta_k = theta0/pi_k;
        press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);
        call nsa_stalaw(1,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Reference Solution1
        rho_hyd   = rho_k
        press_hyd = press_k
        theta_hyd = theta_k

        !Press_hyd
        rekee_nsa(ndime+1, ipoin) = rho_hyd
        !Theta_hyd
        rekee_nsa(ndime+2, ipoin) = theta_hyd
        !Densi_hyd
        rekee_nsa(ndime+3, ipoin) = press_hyd

        !Water tracers
        q_nsa(5,ipoin) = q_ref(5,ipoin)
        q_nsa(6,ipoin) = 0.0_rp
        q_nsa(7,ipoin) = 0.0_rp

        if(kfl_ncons_nsa == 0) then
           !
           ! Momentum in Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = u0*rho_hyd
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0*rho_hyd
        else
           !
           ! Momentum in NON-Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = u0
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0

        end if

        !Perturbation fields: (total - reference)
        bvess_nsa(1,ipoin,1)           = u0 !uvelo
        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo

        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

        !Water tracers:
        bvess_nsa(ndofn_nsa+1,ipoin,1)   = q_nsa(5,ipoin)
        bvess_nsa(ndofn_nsa+2,ipoin,1)   = q_nsa(6,ipoin)
        bvess_nsa(ndofn_nsa+3,ipoin,1)   = q_nsa(7,ipoin)
     end do !ipoin

     !Return press and theta
     outvar = 1

     case(205)
     !
     ! NH Linear mountain with tracer (nvar_nsa = 7)
     !
     !Case Constants
     nvar_nsa = ndofn_nsa + 3_ip
     temp0    = tempe_nsa
     theta0   = tempe_nsa
     c        = cvcoe_nsa/rgasc_nsa
     c2       = rgasc_nsa/cpcoe_nsa

     u0       = 10.0_rp
     w0       =  0.0_rp
     bv       =  0.01_rp
     bv2      = bv*bv
     g2       = grnor_nsa*grnor_nsa

     ac = 1000.0_rp
     hc = 2500.0_rp !Small Mountain
     xc = 0.0

     rc     = 1.0_rp
     xr     = 10000.0_rp
     zr     = 1500.0_rp
     xr_tr  = 2500.0_rp
     zr_tr  = 1000.0_rp

     theta0 = 300.0_rp
     if(thetac_nsa > -999.0_rp) then
        thetac = thetac_nsa
     else
        thetac = 3.0_rp
     end if
     if(tracerc_nsa > -999.0_rp) then
        tracerc = tracerc_nsa
     else
        tracerc = 0.0_rp
     end if

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(ndime,ipoin)

        !Total variables
        pi_k    = g2/(cpcoe_nsa*theta0*bv2)*(exp(- bv2/grnor_nsa*z ) - 1.0_rp) + 1.0_rp
        theta_k = theta0*exp( bv2/grnor_nsa*z )
        rho_k   = pbaro_nsa/(rgasc_nsa*theta_k)*(pi_k)**c
        call nsa_stalaw(2,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)
        !press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);

        !Reference Solution1
        rho_hyd   = rho_k
        press_hyd = press_k
        theta_hyd = theta_k

        !
        ! Total and perturbation Fields (for plotting only):
        !
        q_nsa(1,ipoin) = rho_k
        if( kfl_uniformvelocity_nsa > 0) then
           q_nsa(2,ipoin) =  uvelo_nsa !if commented, horizontal velocity comes from sounding
        else
           !ELSE q takes the value read in the file.
        end if
        q_nsa(3,ipoin) = 0.0_rp
        q_nsa(4,ipoin) = theta_k
        q_nsa(5,ipoin) = q_ref(5,ipoin)
        q_nsa(6,ipoin) = 0.0_rp
        q_nsa(7,ipoin) = 0.0_rp


        !Press_hyd
        rekee_nsa(ndime+1, ipoin) = rho_hyd
        !Theta_hyd
        rekee_nsa(ndime+2, ipoin) = theta_hyd
        !Densi_hyd
        rekee_nsa(ndime+3, ipoin) = press_hyd

        !Water tracers
        rekee_nsa(ndime+3 + 1, ipoin) = q_ref(5,ipoin)
        rekee_nsa(ndime+3 + 2, ipoin) = 0.0_rp
        rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp

        if(kfl_ncons_nsa == 0) then
           !
           ! Momentum in Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = u0*rho_hyd
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0*rho_hyd
        else
           !
           ! Momentum in NON-Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = u0
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0

        end if

        !Perturbation fields: (total - reference)
        bvess_nsa(1,ipoin,1)           = u0 !uvelo
        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo

        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

        !Water tracers:
        bvess_nsa(ndofn_nsa+1,ipoin,1)   = q_nsa(5,ipoin)
        bvess_nsa(ndofn_nsa+2,ipoin,1)   = q_nsa(6,ipoin)
        bvess_nsa(ndofn_nsa+3,ipoin,1)   = q_nsa(7,ipoin)

     end do !ipoin

     !Return press and theta
     outvar = 1

     !
     ! Build the mapping from 2D grid to column.
     ! This is done here because it's only used by the Kessler module:
     !
     call nsa_build_column
     !

     case(202)
     !
     ! Moisture: this case is similar to case 200
     ! except that computes the HS state from analytic equations,
     ! and computes theta_k as 
     ! theta_k = theta_HYD + (theta0 - theta_HYD) + dtheta
     ! The same applies to p_k and rho_k
     !
     nvar_nsa = ndofn_nsa + 3_ip
     c=cvcoe_nsa/rgasc_nsa
     c2=rgasc_nsa/cpcoe_nsa

     xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
     zc = 2000.0_rp

     if(thetac_nsa > -999.0_rp) then
        thetac = thetac_nsa
     else
        thetac = 3.0_rp
     end if
     if(tracerc_nsa > -999.0_rp) then
        tracerc = tracerc_nsa
     else
        tracerc = 0.0_rp
     end if

     rc     = 1.0_rp
     xr     = 10000.0_rp
     zr     = 1500.0_rp
     xr_tr  = 2500.0_rp
     zr_tr  = 1000.0_rp

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(2,ipoin)

        !Reference Fields (hydrostatic equilibrium):
        theta_hyd = tempe_nsa
        pi_hyd    = 1.0_rp - grnor_nsa/(cpcoe_nsa*theta_hyd)*z
        rho_hyd   = pbaro_nsa/(rgasc_nsa*theta_hyd)*(pi_hyd)**c
        call nsa_stalaw(2,0,rho_hyd,press_hyd,theta_hyd,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Background density and theta (from sounding):
        rho0   = q_ref(1,ipoin)
        theta0 = q_ref(4,ipoin)
        press0 = pbaro_nsa*(rho0*rgasc_nsa*theta0/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)

        r=sqrt( ((x-xc)/xr)**2_rp  + ((z-zc)/zr)**2_rp )
        dtheta = 0.0_rp
        qv     = 0.0_rp
        if (r < rc) then
           dtheta = thetac*(cos(pi*r/2.0_rp))**2_rp
        end if

        dtracer = 0.0_rp
        r=sqrt( ((x-xc)/xr_tr)**2  + ((z-zc)/zr_tr)**2 )
        if (r < rc) then
           dtracer = tracerc!*(cos(pi*r/2.0))**2
        end if

        theta_k = theta_hyd + (theta0 - theta_hyd) + dtheta
        press_k = pbaro_nsa*(rho0*rgasc_nsa*theta_k/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)
        rho_k   = pbaro_nsa*(press_k/pbaro_nsa)**(cvcoe_nsa/cpcoe_nsa)/(rgasc_nsa*theta_k)
        !
        ! Total and perturbation Fields (for plotting only):
        !
        q_nsa(1,ipoin) = rho_k
        !q_nsa(2,ipoin) = 0.0_rp !if commented, horizontal velocity comes from sounding
        q_nsa(3,ipoin) = 0.0_rp
        q_nsa(4,ipoin) = theta_k
        q_nsa(5,ipoin) = q_ref(5,ipoin)
        q_nsa(6,ipoin) = 0.0_rp
        q_nsa(7,ipoin) = 0.0_rp

        !Initial values of tracers:
        !Densi_hyd
        rekee_nsa(ndime+1, ipoin) = rho_hyd
        !Theta_hyd
        rekee_nsa(ndime+2, ipoin) = theta_hyd
        !Press_hyd
        rekee_nsa(ndime+3, ipoin) = press_hyd

        !Water tracers
        rekee_nsa(ndime+3 + 1, ipoin) = q_ref(5,ipoin)
        rekee_nsa(ndime+3 + 2, ipoin) = 0.0_rp
        rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp

        if(kfl_ncons_nsa == 0) then
           !
           ! Momentum in Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = q_nsa(2,ipoin)*rho0
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0*rho0
        else
           !
           ! Momentum in NON-Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = q_nsa(2,ipoin)
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0

        end if

        !Velocity:
        bvess_nsa(1,ipoin,1)           = q_nsa(2,ipoin) !uvelo
        bvess_nsa(ndime,ipoin,1)       = w0         !vvelo

        !Total fields:
        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

        !Water tracers:
        bvess_nsa(ndofn_nsa+1,ipoin,1)   = q_nsa(5,ipoin)
        bvess_nsa(ndofn_nsa+2,ipoin,1)   = q_nsa(6,ipoin)
        bvess_nsa(ndofn_nsa+3,ipoin,1)   = q_nsa(7,ipoin)

     end do !ipoin

     !Return press and theta
     outvar = 1

     !
     ! Build the mapping from 2D grid to column.
     ! This is done here because it's only used by the Kessler module:
     !
     call nsa_build_column

     case(201)
     !
     ! KESSLER, SIMPLE
     ! Moisture: this case is similar to case 200
     ! except that it already has a supersaturaded cloud 
     ! in the initial field. qv initial is the same as for case 200.
     !
     nvar_nsa = ndofn_nsa + 3_ip
     c=cvcoe_nsa/rgasc_nsa
     c2=rgasc_nsa/cpcoe_nsa

     xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
     !zc = 1200.0_rp
     zc = 2000.0_rp     
     theta0 = 250.0_rp

     if(thetac_nsa > -999.0_rp) then
        thetac = thetac_nsa
     else
        thetac = 3.0_rp
     end if
     if(tracerc_nsa > -999.0_rp) then
        tracerc = tracerc_nsa
     else
        tracerc = 0.0_rp
     end if

     !     rc     = 1.0_rp
     !     xr     = 10000.0_rp
     !     zr     = 1500.0_rp
     !     xr_tr  = 2500.0_rp
     !     zr_tr  = 1000.0_rp

     rc     = 1.0_rp
     xr     = 500.0_rp
     zr     = 250.0_rp
     xr_tr  = 100.0_rp
     zr_tr  = 25.0_rp

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(2,ipoin)

        !Background density and theta (from sounding):
        rho0   = q_ref(1,ipoin)
        theta0 = q_ref(4,ipoin)
        press0 = pbaro_nsa*(rho0*rgasc_nsa*theta0/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)

        r=sqrt( ((x-xc)/xr)**2_rp  + ((z-zc)/zr)**2_rp )
        dtheta = 0.0_rp
        qv     = 0.0_rp
        if (r < rc) then
           dtheta = thetac*(cos(pi*r/2.0_rp))**2_rp
        end if

        dtracer = 0.0_rp
        r=sqrt( ((x-xc)/xr_tr)**2  + ((z-zc)/zr_tr)**2 )
        if (r < rc) then
           dtracer = tracerc!*(cos(pi*r/2.0))**2
        end if

        theta_k = theta0 + dtheta
        press_k = pbaro_nsa*(rho0*rgasc_nsa*theta_k/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa)
        rho_k   = pbaro_nsa*(press_k/pbaro_nsa)**(cvcoe_nsa/cpcoe_nsa)/(rgasc_nsa*theta_k)
        !call nsa_stalaw(1,rho_k,press_k,theta_k,z,0.0_rp)

        !
        ! Total and perturbation Fields (for plotting only):
        !
        q_nsa(1,ipoin) = rho_k
        q_nsa(2,ipoin) = 0.0 !if commented, horizontal velocity comes from sounding
        q_nsa(3,ipoin) = 0.0_rp
        q_nsa(4,ipoin) = theta_k
        q_nsa(5,ipoin) = q_ref(5,ipoin)
        q_nsa(6,ipoin) = 0.0_rp!dtracer
        q_nsa(7,ipoin) = 0.0_rp

        !Initial values of tracers:
        !Densi_hyd
        rekee_nsa(ndime+1, ipoin) = rho0
        !Theta_hyd
        rekee_nsa(ndime+2, ipoin) = theta0
        !Press_hyd
        rekee_nsa(ndime+3, ipoin) = press0

        !Water tracers
        rekee_nsa(ndime+3 + 1, ipoin) = q_ref(5,ipoin)
        rekee_nsa(ndime+3 + 2, ipoin) = dtracer
        rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp

        if(kfl_ncons_nsa == 0) then
           !
           ! Momentum in Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = q_nsa(2,ipoin)*rho0
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0*rho0
        else
           !
           ! Momentum in NON-Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = q_nsa(2,ipoin)
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0

        end if


        !Velocity:
        bvess_nsa(1,ipoin,1)           = q_nsa(2,ipoin) !uvelo
        bvess_nsa(ndime,ipoin,1)       = w0         !vvelo

        !Total fields:
        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

        !Water tracers:
        bvess_nsa(ndofn_nsa+1,ipoin,1)   = q_nsa(5,ipoin)
        bvess_nsa(ndofn_nsa+2,ipoin,1)   = q_nsa(6,ipoin)
        bvess_nsa(ndofn_nsa+3,ipoin,1)   = q_nsa(7,ipoin)

     end do !ipoin

     !Return press and theta
     outvar = 1

     !
     ! Build the mapping from 2D grid to column.
     ! This is done here because it's only used by the Kessler module:
     !
     call nsa_build_column

     case(203)
     !
     ! Klaassen and Klark 1985:
     !
     nvar_nsa = ndofn_nsa + 3_ip
     c=cvcoe_nsa/rgasc_nsa
     c2=rgasc_nsa/cpcoe_nsa

     u0 = 0.0_rp
     w0 = 0.0_rp

     theta_z0 = 250.16_rp

     if(thetac_nsa > -999.0_rp) then
        thetac = thetac_nsa
     else
        thetac = 3.0_rp
     end if
     if(tracerc_nsa > -999.0_rp) then
        tracerc = tracerc_nsa
     else
        tracerc = 0.0_rp
     end if


     xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
     zc = 500.0_rp

     rc     = 1.0_rp
     xr     = 1500.0_rp
     zr     = 750.0_rp

     !Load Reference Values
     !call nsa_loadref(q,q_ref)

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(2,ipoin)

        !Background theta0:
        !computed from the points of the graph in Klaassen and Klark 1985:
        !Background density and theta (from sounding):
        if( z <= 1400.0_rp) then
           dtdz = 0.36e-03        
           theta0 = dtdz*(z - zmin) + theta_z0

        else if( z > 1400.0_rp .and. z <= 1800.0_rp) then
           dtdz = 4.0e-03
           zmin = 1400.0_rp
           theta_ref = 0.36e-03*1400.0_rp + theta_z0
           theta0 = dtdz*(z - zmin) + theta_ref

        else if( z > 1800.0_rp .and. z <= 5000.0_rp) then
           dtdz = 4.5e-03
           zmin = 1800.0_rp
           theta_ref = 4.0e-03*(1800.0_rp - 1400.0_rp) + 0.36e-03*1400.0_rp + theta_z0
           theta0 = dtdz*(z - zmin) + theta_ref
        end if

        !Background water vapor:
        !Background density and theta (from sounding):
        if( z <= 1400.0_rp) then

           call nsa_compline(0.0_rp,1400.0_rp,0.0085_rp,0.008_rp, m, intercept)
           qv_k = m*z + intercept

        else if( z > 1400.0_rp .and. z <= 1800.0_rp) then

           call nsa_compline(1400.0_rp,1800.0_rp, 0.008_rp,0.0028_rp, m, intercept)
           qv_k = m*z + intercept

        else if( z > 1800.0_rp .and. z <= 3000.0_rp) then

           call nsa_compline(1800.0_rp,3000.0_rp, 0.0028_rp,0.0015_rp, m, intercept)
           qv_k = m*z + intercept

        else if( z > 3000.0_rp .and. z <= 5000.0_rp) then

           call nsa_compline(3000.0_rp,5000.0_rp, 0.0015_rp,0.0005_rp, m, intercept)
           qv_k = m*z + intercept

        end if

        !PErturbation of theta_k
        r=sqrt( ((x-xc)/xr)**2_rp  + ((z-zc)/zr)**2_rp )
        dtheta = 0.0_rp
        if (r < rc) then
           dtheta = thetac*(cos(pi*r/2.0_rp))**2_rp
        end if

        !Total variables:
        theta_k = theta0
        pi_k    = theta_z0/theta_k
        press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);
        call nsa_stalaw(1,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Reference Solution1
        rho_hyd   = rho_k
        press_hyd = press_k
        theta_hyd = theta_k

        !Press_hyd
        rekee_nsa(ndime+1, ipoin) = rho_hyd
        !Theta_hyd
        rekee_nsa(ndime+2, ipoin) = theta_hyd
        !Densi_hyd
        rekee_nsa(ndime+3, ipoin) = press_hyd

        !Water tracers:
        q_nsa(5,ipoin) = qv_k
        q_nsa(6,ipoin) = 0.0_rp
        q_nsa(7,ipoin) = 0.0_rp

        !
        ! Add perturbation to theta and correct the other variables accordingly:
        ! 
        theta_k = theta0 + dtheta
        pi_k    = theta_z0/theta_k
        press_k = pbaro_nsa*pi_k**(cpcoe_nsa/rgasc_nsa);
        call nsa_stalaw(1,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Water tracers
        rekee_nsa(ndime+3 + 1, ipoin) = qv_k
        rekee_nsa(ndime+3 + 2, ipoin) = 0.0_rp
        rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp

        if(kfl_ncons_nsa == 0) then
           !
           ! Momentum in Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = u0*rho_hyd
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0*rho_hyd
        else
           !
           ! Momentum in NON-Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = u0
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0

        end if

        !Perturbation fields: (total - reference)
        bvess_nsa(1,ipoin,1)           = u0 !uvelo
        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo

        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k

        !Water tracers:
        bvess_nsa(ndofn_nsa+1,ipoin,1)   = q_nsa(5,ipoin)
        bvess_nsa(ndofn_nsa+2,ipoin,1)   = q_nsa(6,ipoin)
        bvess_nsa(ndofn_nsa+3,ipoin,1)   = q_nsa(7,ipoin)
     end do !ipoin

     !Return press and theta
     outvar = 1

     !
     ! Build the mapping from 2D grid to column.
     ! This is done here because it's only used by the Kessler module:
     !
     call nsa_build_column

     case(210)
     !
     ! Moist atmosphere in a stably stratified environment: 
     ! Grabowski JAS 2007 and Grbowski & Clark 1991
     ! 

     ! Case Constants
     nvar_nsa = ndofn_nsa + 3_ip

     temp0    = 283.0_rp
     theta0   = 283.0_rp  
     tracer0  = 0.1_rp    !background value
     if(thetac_nsa > -999.0_rp) then
        thetac = thetac_nsa
     else
        thetac = 0.05_rp
     end if
     if(tracerc_nsa > -999.0_rp) then
        tracerc = tracerc_nsa
     else
        tracerc = 0.95_rp
     end if

     !p00   = 85000.0_rp
     p00   = pbaro_nsa

     xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
     zc = 1200.0_rp

     rc    = 500.0_rp
     rc0   = 100.0_rp
     rc1   = 200.0_rp
     rc2   = 300.0_rp
     c     = cvcoe_nsa/rgasc_nsa
     c2    = rgasc_nsa/cpcoe_nsa

     u0    = 0.0_rp
     w0    = 0.0_rp

     static_stability = 1.0e-5_rp
     bv    = sqrt(static_stability*grnor_nsa)
     bv2   = bv*bv
     g2    = grnor_nsa*grnor_nsa

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(ndime,ipoin)

        !Total variables
        pi_k    = g2/(cpcoe_nsa*theta0*bv2)*(exp(- bv2/grnor_nsa*z ) - 1.0_rp) + 1.0_rp
        theta_k = theta0*exp( bv2/grnor_nsa*z )
        rho_k   = p00/(rgasc_nsa*theta_k)*(pi_k)**c
        call nsa_stalaw(2,0,rho_k,press_k,theta_k,z,0.0_rp,0.0_rp,0.0_rp,0.0_rp)

        !Reference Solution1
        rho_hyd   = rho_k
        press_hyd = press_k
        theta_hyd = theta_k

        !Perturbation of RH:
        dtracer = 0.0_rp
        r=sqrt( ((x-xc))**2  + ((z-zc))**2 )
        if(r <= 400.0_rp) then
           dtracer = 0.9_rp
        else if (r > 400.0_rp .and. r <= 500.0_rp) then
           dtracer = tracerc*(1.0_rp - r/500.0_rp)
        end if
!!$        if(r <= rc1) then
!!$           dtracer = tracerc         
!!$        else if (r > rc1 .and. r <= rc2) then
!!$           dtracer = tracerc*(cos(0.5_rp*pi * (r - rc1)/rc0))**2
!!$        end if
        tracer_k = tracer0 + dtracer

        !Relative humidity:
        RH_k = tracer_k

        temp     = theta_k*pi_k
        pressure = press_k

        es  = 1000.0_rp*svp1*exp(svp2*(temp-svpt0)/(temp-svp3))
        qvs = ep2*es/(pressure-es)

        !Derive qv from relative humidity:
        qv_k = RH_k*qvs

        dtheta = 0.0_rp
        if (r < rc) then
           dtheta = thetac*(cos(pi*r/2.0_rp))**2_rp
        end if

        !Press_hyd
        rekee_nsa(ndime+1, ipoin) = rho_hyd
        !Theta_hyd
        rekee_nsa(ndime+2, ipoin) = theta_hyd
        !Densi_hyd
        rekee_nsa(ndime+3, ipoin) = press_hyd

        !Water tracers
        rekee_nsa(ndime+3 + 1, ipoin) = qv_k
        rekee_nsa(ndime+3 + 2, ipoin) = 0.0_rp
        rekee_nsa(ndime+3 + 3, ipoin) = 0.0_rp

        if(kfl_ncons_nsa == 0) then
           !
           ! Momentum in Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = u0*rho_hyd
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0*rho_hyd
        else
           !
           ! Momentum in NON-Conservative form
           ! 
           !u-velo
           rekee_nsa(1, ipoin)       = u0
           !w-velo
           rekee_nsa(ndime, ipoin)   = w0

        end if

        !Perturbation fields: (total - reference)
        bvess_nsa(1,ipoin,1)           = u0 !uvelo
        bvess_nsa(ndime,ipoin,1)       = w0 !vvelo

        bvess_nsa(ndofn_nsa-1,ipoin,1) = rho_k
        bvess_nsa(ndofn_nsa,ipoin,1)   = theta_k + dtheta

        !Water tracers:
        bvess_nsa(ndofn_nsa+1,ipoin,1)   = qv_k
        bvess_nsa(ndofn_nsa+2,ipoin,1)   = 0.0_rp
        bvess_nsa(ndofn_nsa+3,ipoin,1)   = 0.0_rp
     end do !ipoin

     !Return press and theta
     outvar = 1

     !
     ! Build the mapping from 2D grid to column.
     ! This is done here because it's only used by the Kessler module:
     !
     call nsa_build_column

     end select

!
! If RESTART_FILE ON (see *.nsa.dat), read the restart file to initialize
! ONLY bvess_nsa(:,:). rekee_nsa are computed and stored above as usual.
!
     if(kfl_rearst_nsa > 0) &
       call nsa_rearst(outvar)
!
! Compute sponge prameters:
!
     sponge: if( kfl_sponge_nsa > 0) then

     if( kfl_benme_nsa == 3 .or. kfl_benme_nsa == 4 .or. kfl_benme_nsa == 5 .or. &
          kfl_benme_nsa == 5 .or. kfl_benme_nsa == 6 .or. kfl_benme_nsa == 7 .or. &
          kfl_benme_nsa >= 200) then

        !HSlinear with the domain and sponge of the paper
        xmin = xmin_nsa
        xmax = xmax_nsa
        ymin = ymin_nsa
        ymax = ymax_nsa
        zmin = zmin_nsa
        zmax = zmax_nsa

        !Sponge
        dxs=dxs_nsa
        dys=dys_nsa
        dzs=dzs_nsa
        ampx=ampx_nsa
        ampy=ampy_nsa
        ampz=ampz_nsa

        if(kfl_botta_nsa > 0) then
           xmin = -8000.0_rp
           xmax =  8000.0_rp
           zmin =     0.0_rp
           zmax =  8000.0_rp

           dxs=2000.0_rp
           dzs=1000.0_rp
           ampx=0.5_rp
           ampz=1.0_rp
        end if

        if( kfl_benme_nsa == 205 .and. kfl_sptyp_nsa == 2) then     !NHLINEAR-Kessler     
           xmin = -72000.0_rp
           xmax =  72000.0_rp
           zmin =      0.0_rp
           zmax =  30000.0_rp

           dxs=40000.0_rp
           dzs=15000.0_rp
           ampx=1.0_rp
           ampz=0.5_rp

        end if

        if(kfl_sptyp_nsa == 1) then
           !
           ! Simple sponge
           !
           zs  = zmax - dzs
           xsl = xmin + dxs
           xsr = xmax - dxs

           do ipoin= 1,npoin

              cside=0.0_rp; ctop=0.0_rp
              x=coord(1,ipoin)
              z=coord(ndime,ipoin)
              if(ndime > 2) &
                   y=coord(ndime-1,ipoin)

              !Top Sponge
              if (z >= zs) then
                 ctop=ampz*( (z - zs)/(zmax - zs) )**4
              end if

              !X-Lateral Sponge
              if (x <= xsl) then
                 cside = ampx*( (x - xsl)/(xmin - xsl) )**4
              else if (x >= xsr) then
                 cside = ampx*( (x - xsr)/(xmax - xsr) )**4
              end if

              ! boundary damping
              !Store Values
              bb_nsa(ipoin) = ctop + cside
              bb_nsa(ipoin) = min(bb_nsa(ipoin), 1.0_rp)
              aa_nsa(ipoin) = 1.0_rp - bb_nsa(ipoin)

              bspon_nsa(1,ipoin) = aa_nsa(ipoin)
              bspon_nsa(2,ipoin) = bb_nsa(ipoin)

           end do !ipoin        

        else
           !
           ! Lilly and Klemp 1978
           !
           nelx = nelx_nsa
           nely = nelz_nsa

           !zd    = 17700.0_rp   ! m
           !zd    = 11500.0_rp   ! m
           zd    = 13500.0_rp   ! m
           alpha = 0.02_rp      ! s^-1

           !Constants
           ct = 0.5_rp
           cs = 1.0_rp

           !Constants
           ztop   = zmax
           xbound = dxs

           xl     = xmin + xbound
           xr     = xmax - xbound
           !      zbound=0.5*ztop
           !      dz=ztop - zbound
           dsx    = 8000.0_rp
           dsy    =  660.0_rp
           zdcorrect = zd

           !Initialize
           aa_nsa=0; bb_nsa=0
           do ipoin=1,npoin
              ctop=0; cside=0
              x=coord(1,ipoin)
              z=coord(ndime,ipoin)

              !
              ! boundary damping
              !
              dbl   = min(x - xmin, xmax - x) ! distance from the boundary. xs in Restelli's thesis
              if(dbl == 0.0_rp) then
                 dbl   = min(x - xmin + 0.1, xmax - x + 0.1)
              end if
              beta  = (1.0_rp - tanh(dbl/dsx) )/(tanh(dbl/dsx) + 0.1_rp)
              cside = cs*beta

              !
              ! top damping
              ! first layer: damp lee waves
              !
              if(z .le. zdcorrect) then
                 ctop = 0
              else
                 xid = (z-zd)/(ztop-zd) ! normalized coordinate
                 if(xid .lt. 0.5_rp) then
                    abstaud = 0.5_rp*alpha*(1.0_rp - cos(xid*pi))

                 else
                    abstaud = 0.5_rp*alpha*(1.0_rp + (xid - 0.5_rp)*pi)

                 endif

                 !ctop = ct * (dt*abstaud)/(1.0_rp + 0.5_rp*dt*abstaud)
                 ctop = ct*abstaud
              endif

              !
              ! second layer: damp short waves
              !
              dbt = ztop - z + 0.1_rp ! distance from the boundary
              !write(*,*) 'dbt:',dbt,'dsy',dsy
              !write(*,*) tanh(dbt/dsy)
              !!if (dbt == 0.0_rp) then
              !!   dbt = ztop
              !!endif
              beta = (1.0_rp - tanh(dbt/dsy))/(tanh(dbt/dsy) + 0.1_rp)

              ctop = ctop + ct*beta

              !Store Values
              bb_nsa(ipoin) = min(ctop + cside,1.0_rp)
              aa_nsa(ipoin) = 1.0_rp - bb_nsa(ipoin)

              bspon_nsa(1,ipoin) = aa_nsa(ipoin)
              bspon_nsa(2,ipoin) = bb_nsa(ipoin)

           end do !i
        end if

     end if
end if sponge

end subroutine nsa_initial_conditions

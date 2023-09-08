
subroutine chm_initial_conditions(outvar)
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
  use      def_chemic
  use      def_parame
  implicit none

  integer(ip) :: idime,icomp,idofn,ipoin,itime,ibubb,nbubb,kshbu,dummi,i,j,k
  real(rp)    :: velmi,xfact,xfac2,rgacp,xinve,xradi,xdifi,xrano(3),xtemp,rauxi,exner

  real(rp)    :: u0,w0,thetac,tracerc,dtracer,tracer_k,tracer0
  real(rp)    :: temp0,theta0,dtheta,press0,rho0,p00
  real(rp)    :: rc,rc0,rc1,rc2,rc3,xc,yc,zc,r,xr,yr,zr,xr_tr,yr_tr,zr_tr
  real(rp)    :: x,y,z, xmin, xmax, ymin, ymax, zmin, zmax
  real(rp)    :: theta_k,pi_k,rho_k,press_k
  real(rp)    :: theta_hyd,pi_hyd,rho_hyd,press_hyd, rho_baro
  real(rp)    :: rho,theta,pressure,qv_k
  real(rp)    :: a,b,c,d,e,c2,bv,bv2,g2,es,gam,pii,m,intercept
  real(rp)    :: xlv, ep2,svp1,svp2,svp3,svpt0,rhowater,qvs,temp,RH_k,static_stability

  real(rp)    :: q(ndime+5,npoin), q_ref(ndime+5,npoin)
  character   :: fname*72

  !
  ! outvar = 1 --> densi and theta
  ! outvar = 2 --> press and theta
  !
  integer(ip) :: outvar, ntracers
  dummi = 0_ip

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
  
  select case(kfl_benme_chm)

  case(100)
     !
     ! Ahmad warm bubble with tracer:
     !
     !Domain characteristics:
     xmin = -10000.0_rp
     xmax =  10000.0_rp
     zmin = 0.0_rp
     zmax = 10000.0_rp
     if(ndime > 2) then
        ymin = -5000.0_rp
        ymax =  5000.0_rp
     end if

     xc =    0.5_rp*(xmin + xmax)
     zc = 2000.0_rp
     if(ndime > 2) &
          yc = 0.5_rp*(ymin + ymax)

     rc = 2000.0_rp
     xr = 2000.0_rp
     yr = 1000.0_rp
     zr = 2000.0_rp

     !Perturbation intensity:
     tracerc=2.0_rp   ! With bubble

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(ndime,ipoin)
        if(ndime > 2) &
             y=coord(ndime-1,ipoin)

        dtracer=0.0_rp
        if(ndime < 3) then
           !
           ! 2D
           !
           r=sqrt( (x-xc)**2.0_rp + (z-zc)**2.0_rp )   
           if (r <= rc) then
              dtracer=tracerc*(1.0_rp - r)
           end if
        else
           !
           ! 3D
           !
           r=sqrt( ((x-xc)/xr)**2_rp  + ((y-yc)/yr)**2_rp + ((z-zc)/zr)**2_rp )
           if (r <= 1.0_rp) then
              dtracer=tracerc*(1.0_rp - r)

           end if
        end if

        !Tracers:
        !qv
        bvess_chm(1,ipoin) = dtracer
        !qc
        bvess_chm(2,ipoin) = 0.0_rp
        !qr
        bvess_chm(3,ipoin) = 0.0_rp

        q(5,ipoin) = dtracer
        !qc
        q(6,ipoin) = 0.0_rp
        !qr
        q(7,ipoin) = 0.0_rp

     end do !ipoin
     
  case(200)
     !
     ! Moisture (from Jim's):
     !
     ! Load Reference Values
     !
     if(ndime < 3) then
        !
        ! 2D:
        !
        !call chm_loadref_field

        do ipoin=1,npoin
           x=coord(1,ipoin)
           z=coord(ndime,ipoin)
           
           !Tracers:
           !qv
           bvess_chm(1,ipoin) = xfiel(3)%a(1,ipoin,1) !qvaporef_chm(ipoin)
           !qc
           bvess_chm(2,ipoin) = 0.0_rp
           !qr
           bvess_chm(3,ipoin) = 0.0_rp

           q(ndime+2+1,ipoin) = xfiel(3)%a(1,ipoin,1) !qvaporef_chm(ipoin) 
           !qc
           q(ndime+2+1,ipoin) = 0.0_rp
           !qr
           q(ndime+2+1,ipoin) = 0.0_rp
           
        end do !ipoin

     else
        !
        ! 3D:
        !
        !call chm_loadref_field

        do ipoin=1,npoin
           x=coord(1,ipoin)
           z=coord(ndime,ipoin)
           if(ndime > 2) &
                y=coord(ndime-1,ipoin)

           !Tracers:
           !qv
           bvess_chm(1,ipoin) = xfiel(3)%a(1,ipoin,1) !qvaporef_chm(ipoin)

           !qc
           bvess_chm(2,ipoin) = 0.0_rp
           !qr
           bvess_chm(3,ipoin) = 0.0_rp

           q(ndime+2+1,ipoin) = xfiel(3)%a(1,ipoin,1) !qvaporef_chm(ipoin)
           !qc
           q(ndime+2+1,ipoin) = 0.0_rp
           !qr
           q(ndime+2+1,ipoin) = 0.0_rp

        end do !ipoin
     end if

  case(200000)
     !
     ! Moisture (from Jim's):
     !
     ! Load Reference Values
     !
     call chm_loadref(q,q_ref)

     !Tracer perturbation amplitude
     tracerc = 0.002_rp

     !Radii of tracer perturbation
     rc = 1.0_rp
     xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
     zc = 2000.0_rp
     if(ndime > 2) &
          yc = 0.5_rp*(maxval(coord(ndime-1,:)) + minval(coord(ndime-1,:)))

     xr_tr  = 2500.0_rp
     zr_tr  = 1000.0_rp
     if(ndime > 2) &
          yr_tr = 1500.0_rp

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(ndime,ipoin)
        if(ndime > 2) &
             y=coord(2,ipoin)

        dtracer = 0.0_rp
        if(ndime < 3) then
           !
           ! 2D
           !
           r=sqrt( ((x-xc)/xr_tr)**2  + ((z-zc)/zr_tr)**2 )
        else
           !
           ! 3D
           !
           r=sqrt( ((x-xc)/xr_tr)**2 + ((y-yc)/yr_tr)**2 + ((z-zc)/zr_tr)**2 )
        end if

        if (r < rc) then
           dtracer = tracerc!*(cos(pi*r/2.0))**2
        end if

        !Tracers:
        !qv
        bvess_chm(1,ipoin) = q_ref(ndime+3,ipoin)

        !qc
        bvess_chm(2,ipoin) = 0.0_rp
        !qr
        bvess_chm(3,ipoin) = 0.0_rp

        q(ndime+3,ipoin) = q_ref(ndime+3,ipoin)
        !qc
        q(ndime+4,ipoin) = 0.0_rp
        !qr
        q(ndime+5,ipoin) = 0.0_rp

     end do !ipoin

  case(201)
     !
     ! KESSLER, SIMPLE
     !
     ! Moisture: this case is similar to case 200
     ! except that it already has a supersaturaded cloud 
     ! in the initial field. qv initial is the same as for case 200.
     !
     xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
     !zc = 1200.0_rp
     zc = 2000.0_rp

     tracerc= 0.05_rp

     rc     = 1.0_rp
     !xr     = 10000.0_rp
     !zr     = 1500.0_rp

     xr     = 500.0_rp
     zr     = 250.0_rp

     ! Load Reference Values
     call chm_loadref(q,q_ref)

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(2,ipoin)

        dtracer = 0.0_rp
        r=sqrt( ((x-xc)/xr)**2  + ((z-zc)/zr)**2 )
        if (r < rc) then
           dtracer = tracerc*(cos(pi*r/2.0))**2
        end if

        !Tracers:
        !qv
        bvess_chm(1,ipoin) = q_ref(5,ipoin) + dtracer
        !qc
        bvess_chm(2,ipoin) = 0.0_rp!dtracer
        !qr
        bvess_chm(3,ipoin) = 0.0_rp

        q(5,ipoin) = q_ref(5,ipoin) + dtracer
        !qc
        q(6,ipoin) = 0.0_rp!dtracer
        !qr
        q(7,ipoin) = 0.0_rp

     end do !ipoin

  case(203)
     !
     ! Klaassen and Klark 1985:
     !    
     xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
     zc = 800.0_rp
     tracerc = 0.0001_rp

     rc     = 1.0_rp
     xr     = 1000.0_rp
     zr     = 500.0_rp

     ! Load Reference Values
     !    call chm_loadref(q,q_ref)

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(ndime,ipoin)

        !Background density and theta (from sounding):
        if( z <= 1400.0_rp) then

           call chm_compline(0.0_rp,1400.0_rp,0.0085_rp,0.008_rp, m, intercept)
           qv_k = m*z + intercept

        else if( z > 1400.0_rp .and. z <= 1800.0_rp) then

           call chm_compline(1400.0_rp,1800.0_rp, 0.008_rp,0.0028_rp, m, intercept)
           qv_k = m*z + intercept

        else if( z > 1800.0_rp .and. z <= 3000.0_rp) then

           call chm_compline(1800.0_rp,3000.0_rp, 0.0028_rp,0.0015_rp, m, intercept)
           qv_k = m*z + intercept

        else if( z > 3000.0_rp .and. z <= 5000.0_rp) then

           call chm_compline(3000.0_rp,5000.0_rp, 0.0015_rp,0.0005_rp, m, intercept)
           qv_k = m*z + intercept

        end if

        !Vapor perturbation (see Grabowski and Clark 1991):
        dtracer = 0.0_rp
        r=sqrt( ((x-xc)/xr)**2  + ((z-zc)/zr)**2 )
        if (r < rc) then
           dtracer = tracerc*(cos(pi*r/2.0))**2
        end if

        !Tracers:
        !qv
        bvess_chm(1,ipoin) = qv_k + dtracer
        !qc
        bvess_chm(2,ipoin) = 0.0_rp
        !qr
        bvess_chm(3,ipoin) = 0.0_rp

        q(5,ipoin) = qv_k
        !qc
        q(6,ipoin) = 0.0_rp
        !qr
        q(7,ipoin) = 0.0_rp

     end do !ipoin

  case(210)
     !
     ! KESSLER, G2007
     ! Grabowski JAS 2007 and Grbowski & Clark 1991
     ! 
     !
     temp0   = 283.0_rp
     theta0  = 283.0_rp

     !p00   = 85000.0_rp
     p00   = pbaro_chm

     tracer0 = 0.1_rp  !background value
     tracerc = 1.0_rp  !amplitude of the tracer bubble     

     xc = 0.5_rp*(maxval(coord(1,:)) + minval(coord(1,:)))
     zc = 1200.0_rp

     rc0     = 100.0_rp
     rc1     = 200.0_rp
     rc2     = 300.0_rp

     c     = cvcoe_chm/rgasc_chm
     c2    = rgasc_chm/cpcoe_chm

     static_stability = 1.0e-5_rp
     !bv    = 0.01_rp
     bv    = sqrt(static_stability*grnor_chm)
     bv2   = bv*bv
     g2    = grnor_chm*grnor_chm

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(2,ipoin)

        !Perturbation of RH:
        dtracer = 0.0_rp
        r=sqrt( ((x-xc))**2  + ((z-zc))**2 )
        if(r <= 400.0_rp) then
           dtracer = 0.9_rp
        else if (r > 400.0_rp .and. r <= 500.0_rp) then
           dtracer = tracerc*(1.0_rp - r/500.0_rp)
        end if
!!$        
!!$        if(r <= rc1) then
!!$           dtracer = tracerc
!!$           
!!$        else if (r > rc1 .and. r <= rc2) then
!!$           dtracer = tracerc*(cos(0.5_rp*pi * (r - rc1)/rc0))**2
!!$
!!$        end if
        tracer_k = tracer0 + dtracer

        !Total relative humidity:
        RH_k = tracer_k

        !
        ! Convert RH to qv: qv = RH*qvs
        !
        !Total variables
        pi_k    = g2/(cpcoe_chm*theta0*bv2)*(exp(- bv2/grnor_chm*z ) - 1.0_rp) + 1.0_rp
        theta_k = theta0*exp( bv2/grnor_chm*z )
        rho_k   = p00/(rgasc_chm*theta_k)*(pi_k)**c
        call chm_stalaw(2,rho_k,press_k,theta_k,z,0.0_rp)

        temp     = theta_k*pi_k
        pressure = press_k

        es  = 1000.0_rp*svp1*exp(svp2*(temp-svpt0)/(temp-svp3))
        qvs = ep2*es/(pressure-es)

        !Derive qv from relative humidity:
        qv_k = RH_k*qvs

        !Tracers:
        !qv
        bvess_chm(1,ipoin) = qv_k
        !qc
        bvess_chm(2,ipoin) = 0.0_rp
        !qr
        bvess_chm(3,ipoin) = 0.0_rp

        q(5,ipoin) = qv_k
        !qc
        q(6,ipoin) = 0.0_rp
        !qr
        q(7,ipoin) = 0.0_rp

     end do !ipoin

  case(110)
     !
     !Rising Thermal Bubble with passive tracer
     !
     !Case Constants
     ntracers = 1_ip
     c=cvcoe_chm/rgasc_chm
     theta0=300.0_rp
     u0=0.0_rp
     w0=0.0_rp
     xc=500.0_rp
     zc=350.0_rp
     rc=250.0_rp

     !Perturbation intensity:
     thetac=10.0_rp

     do ipoin=1,npoin
        x=coord(1,ipoin)
        z=coord(ndime,ipoin)

        r = sqrt( (x-xc)**2.0_rp + (z-zc)**2.0_rp )
        dtheta = 0.0_rp
        if (r < rc) then
           dtheta=thetac!/2.0_rp*(1.0_rp + cos(pi*r/rc) )
        end if

        !Tracer initial field:
        bvess_chm(1,ipoin) = dtheta

     end do !ipoin

  end select
  
end subroutine chm_initial_conditions


subroutine chm_loadref(q,q_ref)

  use def_master
  use def_domain
  use def_chemic
  use def_parame

  implicit none
  integer(ip) :: ipoin
  real(rp)    :: q(ndime+5,npoin), q_ref(ndime+5,npoin)

  q_ref = 0.0
  open(1,file='case200ref.out')
  if(ndime < 3) then
     !
     ! 2D
     !
     do ipoin = 1,npoin
        !          densi        theta         qv          u-velo (almost sure)
        read(1,*) q_ref(1,ipoin), q_ref(4,ipoin), q_ref(5,ipoin), q(2,ipoin) 
        !print*,'case200: ', q_ref(1,ipoin), q_ref(4,ipoin), q_ref(5,ipoin), q(2,ipoin) 

     end do
  else
     !
     ! 3D
     !
     do ipoin = 1,npoin
        !          densi        theta                   qv          u-velo (almost sure)
        read(1,*) q_ref(1,ipoin), q_ref(5,ipoin), q_ref(6,ipoin), q(2,ipoin) 
        !print*,'case200: ', q_ref(1,ipoin), q_ref(4,ipoin), q_ref(5,ipoin), q(2,ipoin) 

     end do

  end if
  close(1)

end subroutine chm_loadref

subroutine chm_loadref_field

  use def_master
  use def_domain
  use def_chemic
  use def_parame

  implicit none
  integer(ip) :: ipoin, i,j,k
  
  !
  ! Read the field from file passed through *.dom.dat as
  ! FIELDA, END_FIELDS
  !
  qvaporef_chm => xfiel(3)%a(1,:,1)

end subroutine chm_loadref_field

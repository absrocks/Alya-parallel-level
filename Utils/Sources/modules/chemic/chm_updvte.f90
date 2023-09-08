subroutine chm_updvte()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_updvte
  ! NAME 
  !    chm_updvte
  ! DESCRIPTION
  !    This routine updates terminal velocity
  ! USES
  ! USED BY
  !    Chemic
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  implicit none
  integer(ip) :: ipoin,maxit,iclas,it
  real(rp)    :: T,visc,visc0,gi,eps,vold,cd100
  real(rp)    :: rhop,diam,rhoa,Tair0,vset,cd,rey,rk1,rk2,psi
  real(rp)    :: a,b,gama,dnd,M,B0,w
  logical(lg) :: conve

  if( INOTMASTER .and. lawvt_chm >= 1 ) then

     Tair0 = 291.15_rp         ! Reference temperature
     visc0 = 1.827e-5_rp       ! Viscosity at Tair0
     gi    = 9.81_rp           ! Gravity acceleration
     eps   = 1e-3_rp           ! Tolerance
     maxit = 100               ! Maximum number of iterations

     if( lawvt_chm == 1 ) then

        !----------------------------------------------------------------
        !
        ! GANSER model: ut = ut(T,rho)
        !
        !----------------------------------------------------------------

        do iclas=1,nclas_chm
           rhop = rhopa_chm(iclas)
           diam = diame_chm(iclas)
           psi  = shape_chm(iclas)
           call get_gama(diam,spher_chm(iclas),gama)                ! Get gama=a/c
           dnd=0.5_rp*(1.0_rp+gama)*gama**(-2.0_rp/3.0_rp)
!         
           do ipoin = 1,npoin
              T     = tempe_chm(ipoin)
              rhoa  = densi_chm(ipoin)
              visc  = visc0*((Tair0+120.0_rp)/(T+120.0_rp))*((T/Tair0)**1.5_rp)  ! Sutherland's law
              vold  = 1e5_rp 
              cd    = 1.0_rp
              it    = 0
              conve = .false.
              do while( .not. conve )
                 it   = it+1
                 vset = sqrt(4.0_rp*gi*diam*rhop/(3.0_rp*cd*rhoa))
                 rey  = rhoa*vset*diam/visc
                 rk1  = 1.0_rp/(dnd/3.0_rp+2.0_rp/(3.0_rp*sqrt(psi)))
                 rk2  = 10.0_rp**(1.8148_rp*(-log10(psi))**0.5743_rp)
                 cd   = 24.0_rp/(rey*rk1)*(1.0_rp+0.1118_rp*           &
                        (rey*rk1*rk2)**0.6567_rp)+0.4305_rp*rk2/       &
                        (1.0_rp+3305.0_rp/(rey*rk1*rk2))
                 if( it >= maxit .or. abs(vset-vold) <= eps ) conve = .true.
                 vold = vset
              end do
              vterm_chm(ipoin,iclas) = sqrt(4.0_rp*gi*diam*rhop/(3.0_rp*cd*rhoa))
           end do
        end do

     else if( lawvt_chm == 2 ) then

        !----------------------------------------------------------------
        !
        ! WILSON model: ut = ut(T,rho)
        !
        !----------------------------------------------------------------

        do iclas=1,nclas_chm
           rhop = rhopa_chm(iclas)
           diam = diame_chm(iclas)
           psi  = shape_chm(iclas)
           do ipoin=1,npoin
              T     = tempe_chm(ipoin)
              rhoa  = densi_chm(ipoin)
              visc  = visc0*((Tair0+120.0_rp)/(T+120.0_rp))*((T/Tair0)**1.5_rp)  ! Sutherland's law
              vold  = 1e5_rp 
              cd    = 1.0_rp
              it    = 0
              conve = .false.
              do while( .not. conve )
                 it   = it+1
                 vset = sqrt(4.0_rp*gi*diam*rhop/(3.0_rp*cd*rhoa))
                 rey  = rhoa*vset*diam/visc
                 if(rey.le.100.0_rp) then
                     cd=24.0_rp/rey*psi**(-0.828_rp)+2.0_rp*sqrt(1.0_rp-psi)
                 else if(rey.gt.100.0_rp.and.rey.lt.1000.0_rp) then
                     cd100=0.24_rp*psi**(-0.828_rp)+2.0_rp*sqrt(1.0_rp-psi)
                     a=(1.0_rp-cd100)/900.0_rp
                     b=1.0_rp-1000.0_rp*a
                     cd=a*rey+b
                else
                     cd=1.0_rp
                endif
                if( it >= maxit .or. abs(vset-vold) <= eps ) conve = .true.
                vold = vset
              end do
              vterm_chm(ipoin,iclas) = sqrt(4.0_rp*gi*diam*rhop/(3.0_rp*cd*rhoa))
           end do
        end do                

     else if( lawvt_chm == 3 ) then

        !----------------------------------------------------------------
        !
        ! DELLINO model: ut = ut(rho)
        !
        !----------------------------------------------------------------

        do iclas = 1,nclas_chm
           rhop = rhopa_chm(iclas)
           diam = diame_chm(iclas)
           psi  = shape_chm(iclas)
           do ipoin = 1,npoin
              rhoa = densi_chm(ipoin)
              visc = visc0*((Tair0+120.0_rp)/(T+120.0_rp))*((T/Tair0)**1.5_rp)  ! Sutherland's law
              vset = ((diam*diam*diam*gi*(rhop-rhoa)*rhoa*(psi**1.6_rp))/(visc*visc))**0.5206_rp
              vset = (1.2065_rp*visc*vset)/(diam*rhoa)
              vterm_chm(ipoin,iclas) = vset 
          end do
        end do   

     else if( lawvt_chm == 4 ) then

        !----------------------------------------------------------------
        !
        ! W2Plastics: Modify GANSER model: ut = ut(rho)
        !
        !----------------------------------------------------------------

        do iclas=1,nclas_chm
           rhop = rhopa_chm(iclas)
           diam = diame_chm(iclas)
           psi  = shape_chm(iclas)
           call get_gama(diam,spher_chm(iclas),gama)                ! Get gama=a/c
           dnd=0.5_rp*(1.0_rp+gama)*gama**(-2.0_rp/3.0_rp)

           do ipoin = 1,npoin
              rhoa  = denma_chm
              visc  = 1.0e-3_rp
              vold  = 1e5_rp 
              cd    = 1.0_rp
              it    = 0
              !M     = 470.0_rp
              M     = 300.0_rp
              B0    = 0.6_rp
              w     = 0.24_rp

              !rhoa  = rhoa - 2.0_rp*pi*M*B0/(gi*w)*exp(-2.0_rp*pi*(0.1_rp-coord(ndime,ipoin))/w)
              rhoa  = rhoa + 2.0_rp*pi*M*B0/(gi*w)*exp(-2.0_rp*pi*(coord(ndime,ipoin)+0.04_rp)/w)

              conve = .false.
              do while( .not. conve )
                 it   = it+1
                 vset = sqrt(4.0_rp*gi*diam*abs(rhop-rhoa)/(3.0_rp*cd*abs(rhoa)))
                 rey  = abs(rhoa)*vset*diam/visc
                 rk1  = 1.0_rp/(dnd/3.0_rp+2.0_rp/(3.0_rp*sqrt(psi)))
                 rk2  = 10.0_rp**(1.8148_rp*(-log10(psi))**0.5743_rp)
                 cd   = 24.0_rp/(rey*rk1)*(1.0_rp+0.1118_rp*           &
                        (rey*rk1*rk2)**0.6567_rp)+0.4305_rp*rk2/       &
                        (1.0_rp+3305.0_rp/(rey*rk1*rk2))
                 if( it >= maxit .or. abs(vset-vold) <= eps ) conve = .true.
                 vold = vset
              end do
!print*,rhop,rhoa
              if( rhop > rhoa ) then
                 vterm_chm(ipoin,iclas) =  sqrt(4.0_rp*gi*diam*abs(rhop-rhoa)/(3.0_rp*cd*abs(rhoa)))
              else
                 vterm_chm(ipoin,iclas) = -sqrt(4.0_rp*gi*diam*abs(rhop-rhoa)/(3.0_rp*cd*abs(rhoa))) 
              end if

           end do
        end do
       
     end if

  end if


end subroutine chm_updvte

subroutine elmshc(&
     pnode,pgaus,ptopo,pelty,nddif,kfl_shock,shock,dtinv,gpden,&
     gpvel,gprhs,gpstp,gpdif,gprea,gpgrd,elunk,gpsha,gpcar,&
     gplap,gphes,chale,CD,CF)
  !-----------------------------------------------------------------------
  !****f* mathru/shocel
  ! NAME
  !   shocel
  ! DESCRIPTION
  !   This routine computes the contribution from the shock capturing
  ! OUTPUT
  ! USES
  ! USED BY
  !   elmadr
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mgaus,mnode,ntens
  use def_elmtyp, only       :  TRI03,TET04
  use def_master, only       :  zeror
  implicit none  
  integer(ip), intent(in)    :: pnode,pgaus,ptopo,pelty,nddif,kfl_shock
  real(rp),    intent(in)    :: shock,dtinv
  real(rp),    intent(in)    :: gpvel(ndime,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gplap(mnode,pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)
  real(rp),    intent(in)    :: elunk(pnode)
  real(rp),    intent(in)    :: gpden(pgaus),gprhs(pgaus)
  real(rp),    intent(in)    :: gpdif(nddif,pgaus)
  real(rp),    intent(in)    :: gprea(pgaus)
  real(rp),    intent(in)    :: gpgrd(ndime,pgaus)
  real(rp),    intent(in)    :: chale,gpstp(pgaus) 
  real(rp),    intent(out)   :: CD(pgaus),CF(pgaus)                 
  real(rp)                   :: gpve2,grnor,vepar,facta
  real(rp)                   :: gpgru(3),factt,gpres
  real(rp)                   :: F1,F2,F3,F4,cdv2,umbra
  integer(ip)                :: inode,igaus
  !
  ! Factor according to element topology
  !
  if( ptopo == 0 .or. ptopo == -1 ) then
     factt = 1.5_rp                ! Quadrilateral/Hexahedra or Bars
  else if( ptopo == 1 ) then
     factt = 0.7_rp                ! Triangles/Tetrahedra
  else
     factt = 1.0_rp                ! Other elements
  end if
  factt = 0.5_rp*factt
  umbra = 1.0e-6_rp
  !
  ! Isotropic/anisotropic shock capturing
  !
  facta = real(kfl_shock-1_ip,rp)
  !
  ! Coarse grid residual |R(u)|
  ! |R(u)| = | f - rho*u/(theta*dt) - rho*a.grad(u) + grad(k).grad(u) + k*Lap(u) - s*u|
  ! f  = Q + rho*u^n/(theta*dt)
  !

  if( nddif == 1 ) then

     if( pgaus == 1 .and. pelty == TRI03 ) then

        !-------------------------------------------------------------------
        !
        ! P1 in 2D and 1 Gauss point
        !
        !-------------------------------------------------------------------

        gpres =  gprhs(1)                                    ! Residual
        F1    = -gpden(1) * dtinv - gprea(1)
        F2    = -gpden(1) * gpvel(1,1) + gpgrd(1,1)
        F3    = -gpden(1) * gpvel(2,1) + gpgrd(2,1)     
        do inode=1,pnode
           gpres = gpres+elunk(inode)*(&
                &  +F1*gpsha(inode,1)&
                &  +F2*gpcar(1,inode,1)&
                &  +F3*gpcar(2,inode,1)&
                &  +gpdif(1,1)*gplap(inode,1))
        end do
        gpres    = abs(gpres)
        gpve2    =   gpvel(1,1)*gpvel(1,1)&              ! a^2
             &     + gpvel(2,1)*gpvel(2,1)
        gpgru(1) = 0.0_rp
        gpgru(2) = 0.0_rp
        do inode=1,pnode                                         ! grad(u)
           gpgru(1) = gpgru(1)+gpcar(1,inode,1)*elunk(inode)
           gpgru(2) = gpgru(2)+gpcar(2,inode,1)*elunk(inode)
        end do
        grnor=sqrt( gpgru(1)*gpgru(1) + gpgru(2)*gpgru(2) )      ! | grad(u) | 

        if( grnor>umbra ) then
           vepar = gpres/grnor
           CD(1) = factt*shock*chale*vepar-gpdif(1,1)
           if( CD(1)>0.0_rp .and. gpve2>umbra ) then
              cdv2  = CD(1)/gpve2
              CF(1) = max(0.0_rp,cdv2-gpden(1)*gpstp(1)*gpden(1)) - cdv2
              CF(1) = facta*CF(1)
           end if
        end if

     else if( ndime == 1 ) then

        !-------------------------------------------------------------------
        !
        ! 1D
        !
        !-------------------------------------------------------------------

        do igaus=1,pgaus
           gpres =  gprhs(igaus)                                    ! Residual
           F1    = -gpden(igaus) * dtinv - gprea(igaus)
           F2    = -gpden(igaus) * gpvel(1,igaus) + gpgrd(1,igaus)  
           do inode=1,pnode
              gpres = gpres+elunk(inode)*(&
                   &  +F1*gpsha(inode,igaus)&
                   &  +F2*gpcar(1,inode,igaus)&
                   &  +gpdif(1,igaus)*gplap(inode,igaus))
           end do
           gpres    = abs(gpres)
           gpve2    =   gpvel(1,igaus)*gpvel(1,igaus)               ! a^2
           gpgru(1) = 0.0_rp
           do inode=1,pnode                                         ! grad(u)
              gpgru(1) = gpgru(1)+gpcar(1,inode,igaus)*elunk(inode)
           end do
           grnor=abs(gpgru(1))                                      ! | grad(u) | 

           if( grnor>umbra ) then
              vepar     = gpres/grnor
              CD(igaus) = factt*shock*chale*vepar-gpdif(1,igaus)
              if( CD(igaus)>0.0_rp .and. gpve2>umbra ) then
                 cdv2      = CD(igaus)/gpve2
                 CF(igaus) = max(0.0_rp,cdv2-gpden(igaus)*gpstp(igaus)*gpden(igaus)) - cdv2
                 CF(igaus) = facta*CF(igaus)
              end if
           end if

        end do

     else if( ndime == 2 ) then

        !-------------------------------------------------------------------
        !
        ! 2D
        !
        !-------------------------------------------------------------------

        do igaus=1,pgaus
           gpres    =  gprhs(igaus)                                 ! Residual
           F1       = -gpden(igaus) * dtinv - gprea(igaus)
           F2       = -gpden(igaus) * gpvel(1,igaus) + gpgrd(1,igaus)
           F3       = -gpden(igaus) * gpvel(2,igaus) + gpgrd(2,igaus)   
           gpgru(1) =  0.0_rp
           gpgru(2) =  0.0_rp
           do inode = 1,pnode
              gpres   = gpres + elunk(inode)*(     &
                   &   + F1 * gpsha(inode,igaus)   &
                   &   + F2 * gpcar(1,inode,igaus) &
                   &   + F3 * gpcar(2,inode,igaus) &
                   &   + gpdif(1,igaus) * gplap(inode,igaus))
              gpgru(1) = gpgru(1) + gpcar(1,inode,igaus) * elunk(inode)
              gpgru(2) = gpgru(2) + gpcar(2,inode,igaus) * elunk(inode)
           end do
           gpres = abs(gpres)
           gpve2 =   gpvel(1,igaus)*gpvel(1,igaus) &              ! a^2
                &  + gpvel(2,igaus)*gpvel(2,igaus)
           grnor = sqrt( gpgru(1)*gpgru(1) + gpgru(2)*gpgru(2) + zeror )  ! | grad(u) | 

           if( grnor > umbra ) then
              vepar     = gpres / grnor
              CD(igaus) = factt * shock * chale * vepar - gpdif(1,igaus)
              if( CD(igaus) > 0.0_rp .and. gpve2 > umbra ) then
                 cdv2      = CD(igaus) / gpve2
                 CF(igaus) = max(0.0_rp,cdv2-gpden(igaus)*gpstp(igaus)*gpden(igaus)) - cdv2
                 CF(igaus) = facta * CF(igaus)
              end if
           end if

        end do

     else

        !-------------------------------------------------------------------
        !
        ! 3D
        !
        !-------------------------------------------------------------------

        do igaus=1,pgaus
           gpres    =  gprhs(igaus)                                    ! Residual
           F1       = -gpden(igaus) * dtinv - gprea(igaus)
           F2       = -gpden(igaus) * gpvel(1,igaus) + gpgrd(1,igaus)
           F3       = -gpden(igaus) * gpvel(2,igaus) + gpgrd(2,igaus)
           F4       = -gpden(igaus) * gpvel(3,igaus) + gpgrd(3,igaus)  
           gpgru(1) = 0.0_rp
           gpgru(2) = 0.0_rp
           gpgru(3) = 0.0_rp
           do inode=1,pnode
              gpres    = gpres + elunk(inode)*(&
                   &     + F1 * gpsha(inode,igaus)&
                   &     + F2 * gpcar(1,inode,igaus)&
                   &     + F3 * gpcar(2,inode,igaus)&
                   &     + F4 * gpcar(3,inode,igaus)&
                   &     + gpdif(1,igaus) * gplap(inode,igaus))
              gpgru(1) = gpgru(1) + gpcar(1,inode,igaus) * elunk(inode) ! grad(u)
              gpgru(2) = gpgru(2) + gpcar(2,inode,igaus) * elunk(inode)
              gpgru(3) = gpgru(3) + gpcar(3,inode,igaus) * elunk(inode)
           end do
           gpres    = abs(gpres)
           gpve2    =   gpvel(1,igaus)*gpvel(1,igaus)&                  ! a^2
                &     + gpvel(2,igaus)*gpvel(2,igaus)& 
                &     + gpvel(3,igaus)*gpvel(3,igaus)
           grnor=sqrt(   gpgru(1)*gpgru(1) &                            ! | grad(u) | 
                &      + gpgru(2)*gpgru(2) &
                &      + gpgru(3)*gpgru(3) + zeror ) 

           if( grnor>umbra ) then
              vepar     = gpres/grnor
              CD(igaus) = factt*0.5_rp*shock*chale*vepar-gpdif(1,igaus)
              if( CD(igaus)>0.0_rp .and. gpve2>umbra ) then
                 cdv2      = CD(igaus)/gpve2
                 CF(igaus) = max(0.0_rp,cdv2-gpden(igaus)*gpstp(igaus)*gpden(igaus)) - cdv2
                 CF(igaus) = facta*CF(igaus)
              end if
           end if

        end do

     end if

  else if( nddif == ndime ) then

     if( pgaus == 1 .and. pelty == TRI03 ) then

        !-------------------------------------------------------------------
        !
        ! P1 in 2D and 1 Gauss point
        !
        !-------------------------------------------------------------------

        gpres =  gprhs(1)                                    ! Residual
        F1    = -gpden(1) * dtinv - gprea(1)
        F2    = -gpden(1) * gpvel(1,1) + gpgrd(1,1)
        F3    = -gpden(1) * gpvel(2,1) + gpgrd(2,1)     
        do inode=1,pnode
           gpres = gpres+elunk(inode)*(&
                &  +F1*gpsha(inode,1)&
                &  +F2*gpcar(1,inode,1)&
                &  +F3*gpcar(2,inode,1)&
                &  +gpdif(1,1)*gphes(1,inode,1)&           ! kx*d^2N/dx^2
                &  +gpdif(2,1)*gphes(2,inode,1)&           ! ky*d^2N/dy^2
                )
        end do
        gpres    = abs(gpres)
        gpve2    =   gpvel(1,1)*gpvel(1,1)&              ! a^2
             &     + gpvel(2,1)*gpvel(2,1)
        gpgru(1) = 0.0_rp
        gpgru(2) = 0.0_rp
        do inode=1,pnode                                         ! grad(u)
           gpgru(1) = gpgru(1)+gpcar(1,inode,1)*elunk(inode)
           gpgru(2) = gpgru(2)+gpcar(2,inode,1)*elunk(inode)
        end do
        grnor=sqrt( gpgru(1)*gpgru(1) + gpgru(2)*gpgru(2) + zeror )      ! | grad(u) | 

        if( grnor>umbra ) then
           vepar = gpres/grnor
           CD(1) = factt*shock*chale*vepar-sqrt(gpdif(1,1)**2+gpdif(2,1)**2+ zeror )
           if( CD(1)>0.0_rp .and. gpve2>umbra ) then
              cdv2  = CD(1)/gpve2
              CF(1) = max(0.0_rp,cdv2-gpden(1)*gpstp(1)*gpden(1)) - cdv2
              CF(1) = facta*CF(1)
           end if
        end if

     else if( ndime == 1 ) then

        !-------------------------------------------------------------------
        !
        ! 1D
        !
        !-------------------------------------------------------------------

        do igaus=1,pgaus
           gpres =  gprhs(igaus)                                    ! Residual
           F1    = -gpden(igaus) * dtinv - gprea(igaus)
           F2    = -gpden(igaus) * gpvel(1,igaus) + gpgrd(1,igaus)  
           do inode=1,pnode
              gpres = gpres+elunk(inode)*(&
                   &  +F1*gpsha(inode,igaus)&
                   &  +F2*gpcar(1,inode,igaus)&
                   &  +gpdif(1,igaus)*gphes(1,inode,igaus))
           end do
           gpres    = abs(gpres)
           gpve2    =   gpvel(1,igaus)*gpvel(1,igaus)               ! a^2
           gpgru(1) = 0.0_rp
           do inode=1,pnode                                         ! grad(u)
              gpgru(1) = gpgru(1)+gpcar(1,inode,igaus)*elunk(inode)
           end do
           grnor=abs(gpgru(1))                                      ! | grad(u) | 

           if( grnor>umbra ) then
              vepar     = gpres/grnor
              CD(igaus) = factt*shock*chale*vepar-gpdif(1,igaus)
              if( CD(igaus)>0.0_rp .and. gpve2>umbra ) then
                 cdv2      = CD(igaus)/gpve2
                 CF(igaus) = max(0.0_rp,cdv2-gpden(igaus)*gpstp(igaus)*gpden(igaus)) - cdv2
                 CF(igaus) = facta*CF(igaus)
              end if
           end if

        end do

     else if( ndime == 2 ) then

        !-------------------------------------------------------------------
        !
        ! 2D
        !
        !-------------------------------------------------------------------

        do igaus=1,pgaus
           gpres    =  gprhs(igaus)                                 ! Residual
           F1       = -gpden(igaus) * dtinv - gprea(igaus)
           F2       = -gpden(igaus) * gpvel(1,igaus) + gpgrd(1,igaus)
           F3       = -gpden(igaus) * gpvel(2,igaus) + gpgrd(2,igaus)   
           gpgru(1) =  0.0_rp
           gpgru(2) =  0.0_rp
           do inode = 1,pnode
              gpres   = gpres + elunk(inode)*(     &
                   &   + F1 * gpsha(inode,igaus)   &
                   &   + F2 * gpcar(1,inode,igaus) &
                   &   + F3 * gpcar(2,inode,igaus) &
                   &   + gpdif(1,igaus) * gphes(1,inode,igaus) &
                   &   + gpdif(2,igaus) * gphes(2,inode,igaus) &
                   )
              gpgru(1) = gpgru(1) + gpcar(1,inode,igaus) * elunk(inode)
              gpgru(2) = gpgru(2) + gpcar(2,inode,igaus) * elunk(inode)
           end do
           gpres = abs(gpres)
           gpve2 =   gpvel(1,igaus)*gpvel(1,igaus) &              ! a^2
                &  + gpvel(2,igaus)*gpvel(2,igaus)
           grnor = sqrt( gpgru(1)*gpgru(1) + gpgru(2)*gpgru(2) + zeror )  ! | grad(u) | 

           if( grnor > umbra ) then
              vepar     = gpres / grnor
              CD(igaus) = factt * shock * chale * vepar - sqrt(gpdif(1,igaus)**2+gpdif(2,igaus)**2+ zeror )
              if( CD(igaus) > 0.0_rp .and. gpve2 > umbra ) then
                 cdv2      = CD(igaus) / gpve2
                 CF(igaus) = max(0.0_rp,cdv2-gpden(igaus)*gpstp(igaus)*gpden(igaus)) - cdv2
                 CF(igaus) = facta * CF(igaus)
              end if
           end if

        end do

     else

        !-------------------------------------------------------------------
        !
        ! 3D
        !
        !-------------------------------------------------------------------

        do igaus=1,pgaus
           gpres    =  gprhs(igaus)                                    ! Residual
           F1       = -gpden(igaus) * dtinv - gprea(igaus)
           F2       = -gpden(igaus) * gpvel(1,igaus) + gpgrd(1,igaus)
           F3       = -gpden(igaus) * gpvel(2,igaus) + gpgrd(2,igaus)
           F4       = -gpden(igaus) * gpvel(3,igaus) + gpgrd(3,igaus)  
           gpgru(1) = 0.0_rp
           gpgru(2) = 0.0_rp
           gpgru(3) = 0.0_rp
           do inode=1,pnode
              gpres    = gpres + elunk(inode)*(&
                   &     + F1 * gpsha(inode,igaus)&
                   &     + F2 * gpcar(1,inode,igaus)&
                   &     + F3 * gpcar(2,inode,igaus)&
                   &     + F4 * gpcar(3,inode,igaus)&
                   &     + gpdif(1,igaus) * gphes(1,inode,igaus)&
                   &     + gpdif(2,igaus) * gphes(2,inode,igaus)&
                   &     + gpdif(3,igaus) * gphes(3,inode,igaus)&
                   )
              gpgru(1) = gpgru(1) + gpcar(1,inode,igaus) * elunk(inode) ! grad(u)
              gpgru(2) = gpgru(2) + gpcar(2,inode,igaus) * elunk(inode)
              gpgru(3) = gpgru(3) + gpcar(3,inode,igaus) * elunk(inode)
           end do
           gpres    = abs(gpres)
           gpve2    =   gpvel(1,igaus)*gpvel(1,igaus)&                  ! a^2
                &     + gpvel(2,igaus)*gpvel(2,igaus)& 
                &     + gpvel(3,igaus)*gpvel(3,igaus)
           grnor=sqrt(   gpgru(1)*gpgru(1) &                            ! | grad(u) | 
                &      + gpgru(2)*gpgru(2) &
                &      + gpgru(3)*gpgru(3) ) 

           if( grnor>umbra ) then
              vepar     = gpres/grnor
              CD(igaus) = factt*0.5_rp*shock*chale*vepar&
                   &      -sqrt(gpdif(1,igaus)**2+gpdif(2,igaus)**2+gpdif(3,igaus)**2+ zeror )
              if( CD(igaus)>0.0_rp .and. gpve2>umbra ) then
                 cdv2      = CD(igaus)/gpve2
                 CF(igaus) = max(0.0_rp,cdv2-gpden(igaus)*gpstp(igaus)*gpden(igaus)) - cdv2
                 CF(igaus) = facta*CF(igaus)
              end if
           end if

        end do

     end if

  end if

  do igaus = 1,pgaus
     CD(igaus) = max( CD(igaus) , 0.0_rp )
     CF(igaus) = min( CF(igaus) , 0.0_rp )
  end do

end subroutine elmshc
!------------------------------------------------------------------------
! NOTES
!
! Shock capturing for the ADR equation
! 
!     du
! rho -- + rho*a.grad(u) - div[k*grad(u)] + s*u = f
!     dt
!
! k    = Diffusion                  [M/(L*T)]    
! R    = Residual of the equation   [M*U/(L^3*T)]
! C    = Shock capturing constant 
! tau  = Stabilization parameter: 
!        Its units are h/(rho*a)=   [L^3*T/M]
!        so that rho*tau is in [T]
! kiso = Isotropic SC diffusion     [M/(L*T)]
! k'   = Anisotropic SC diffusion   [M/(L*T)]
!
!        1    |R|    h  
! Pe   = - --------- - 
!        2 |grad(u)| k 
!
!            +          2k   +           R
! ac   = max | 0 , C - ----- | , a*= ----------- grad(u), therefore
!            +         |a*|h +       |grad(u)|^2
!
!            +           1   +
! ac   = max | 0 , C - ----- |
!            +          Pe   +
!        1          |R|
! kiso = - ac*h  --------- , k'=(rho*tau)*rho*a^2
!        2       |grad(u)|
!
! Isotropic:     kiso*grad(u).grad(v)  
!                                         a x a
! Anisotropic:   (<kiso-k'>-kiso)*grad(u).-----.grad(v)
!                                          a^2
!***
!------------------------------------------------------------------------

subroutine chm_inibcs()
  !-----------------------------------------------------------------------
  !****f* partis/chm_inibcs
  ! NAME
  !    chm_inibcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions
  ! OUTPUT 
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_chemic
  use mod_opebcs
  use def_kermod , only : thicl
  use mod_ker_proper

  implicit none
  integer(ip)  :: ipoin,iclas,iboun,inode,ielem,ieset
  real(rp)     :: T,value_xxx,f
  real(rp)     :: Heaviside,x,y,rc,xc,yc
  integer(ip)  :: dummi=1

  if( INOTMASTER ) then

     !-------------------------------------------------------------
     !
     ! Allocate memory
     !
     !-------------------------------------------------------------

     call chm_membcs(1_ip)

     !-------------------------------------------------------------
     !
     ! Node codes
     !
     !-------------------------------------------------------------

     if( kfl_icodn > 0 ) then

        iffun     =  0
        ifloc     =  0
        ifbop     =  0

        if( kfl_allcl_chm == 1 ) then

           iclas     =  1
           kfl_fixno => kfl_fixno_chm(iclas:iclas,:)
           bvess     => bvess_chm(iclas:iclas,:)
           tncod     => tncod_chm(iclas:)
           call reacod(10_ip)
           do ipoin = 1,npoin
              do iclas = 2,nclas_chm
                 kfl_fixno_chm(iclas,ipoin) = kfl_fixno_chm(1,ipoin) 
                 bvess_chm(iclas,ipoin)     = bvess_chm(1,ipoin)
              end do
           end do

        else

           do iclas = 1,nclas_chm
              kfl_fixno => kfl_fixno_chm(iclas:iclas,:)
              bvess     => bvess_chm(iclas:iclas,:)
              tncod     => tncod_chm(iclas:)
              call reacod(10_ip)
           end do
        end if

     end if

     !-------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------

     if( kfl_icodb > 0 ) then

        do iclas = 1,nclas_chm
           kfl_fixbo => kfl_fixbo_chm(iclas,:)
           tbcod     => tbcod_chm(iclas:)
           call reacod(20_ip)           
        end do
      
     end if

     !-------------------------------------------------------------
     !
     ! Initial solution
     !  
     !-------------------------------------------------------------
     
     do iclas = 1,nclas_chm
        if( kfl_initi_chm(iclas) == 1 ) then
           do ipoin=1,npoin
              if(kfl_fixno_chm(iclas,ipoin)==0) then
                 call chm_usrtem(T)
                 bvess_chm(iclas,ipoin)= equil_chm(1,iclas)*exp(-equil_chm(2,iclas)/(boltz_chm*T))
              end if
           end do
        else if( kfl_initi_chm(iclas) == 2 ) then
           do ipoin=1,npoin
              if(kfl_fixno_chm(iclas,ipoin)==0) then
                 bvess_chm(iclas,ipoin)=xinit_chm(iclas,1)
              end if
           end do
        else if( kfl_initi_chm(iclas) == 3 ) then ! Levels
           do ipoin=1,npoin
              if(kfl_fixno_chm(iclas,ipoin)==0) then
                 if ( fleve(ipoin,1)>=0) then
                    Heaviside=1.0_rp
                 else
                    Heaviside=0.0_rp
                 endif
                 bvess_chm(iclas,ipoin)=xinit_chm(iclas,1) + (xinit_chm(iclas,2)-xinit_chm(iclas,1))*Heaviside
              endif
           end do

        else if( kfl_initi_chm(iclas) == 4 ) then
           !
           !  circle of radius rc arround xc,yc
           !
           xc=0.005_rp
           yc=0.005_rp
           rc=0.002_rp
           
           do ipoin=1,npoin
              if(kfl_fixno_chm(iclas,ipoin)==0) then
                 x=coord(1,ipoin)
                 y=coord(2,ipoin)
                 if ( ((x-xc)**2+(y-yc)**2) - rc**2 >= 0.0_rp) then
                    Heaviside=1.0_rp
                 else
                    Heaviside=0.0_rp
                 endif
                 bvess_chm(iclas,ipoin)=xinit_chm(iclas,1) + (xinit_chm(iclas,2)-xinit_chm(iclas,1))*Heaviside
              endif
           end do


        else if( kfl_initi_chm(iclas) < 0 ) then
           ieset = -kfl_initi_chm(iclas)
           if( neset == 0 ) call runend('CHM_INIBCS: ELEMENT SETS NOT DEFINED')
           do ielem = 1,nelem
              if( leset(ielem) == ieset ) then
                 do inode = 1,lnnod(ielem)
                    ipoin = lnods(inode,ielem)
                    if(kfl_fixno_chm(iclas,ipoin)==0) then
                       bvess_chm(iclas,ipoin)=xinit_chm(iclas,1)
                    end if
                 end do
              end if
           end do
        end if
     end do
!     do iclas = nclas_chm+1,nspec_chm
     do iclas = 1,nspec_chm
        if( kfl_initi_chm(iclas) == 2 ) then
           do ipoin=1,npoin
              conce(ipoin,iclas,1)=xinit_chm(iclas,1)
              conce(ipoin,iclas,2)=xinit_chm(iclas,1)
              conce(ipoin,iclas,3)=xinit_chm(iclas,1)
           end do
        else if (kfl_initi_chm(iclas) == 3) then
           
           do ipoin=1,npoin
              if ( fleve(ipoin,1)>=0) then
                 Heaviside=1.0_rp
              else
                 Heaviside=0.0_rp
              endif
              value_xxx= xinit_chm(iclas,1)  + (xinit_chm(iclas,2)-xinit_chm(iclas,1)) * Heaviside
              conce(ipoin,iclas,1)=value_xxx
              conce(ipoin,iclas,2)=value_xxx
              conce(ipoin,iclas,3)=value_xxx
           end do

        else if (kfl_initi_chm(iclas) == 4) then
           
           xc= 0.005
           yc= 0.005
           rc= 0.002
           do ipoin=1,npoin
              x=coord(1,ipoin)
              y=coord(2,ipoin)
              if ( ((x-xc)**2+(y-yc)**2) - rc**2 >= 0.0_rp) then
                 Heaviside=1.0_rp
              else
                 Heaviside=0.0_rp
              endif
              value_xxx= xinit_chm(iclas,1)  + (xinit_chm(iclas,2)-xinit_chm(iclas,1)) * Heaviside
              conce(ipoin,iclas,1)=value_xxx
              conce(ipoin,iclas,2)=value_xxx
              conce(ipoin,iclas,3)=value_xxx
           end do

        end if
     end do

     !-------------------------------------------------------------
     !
     ! User bc
     !
     !-------------------------------------------------------------

     do iclas = 1,nclas_chm
        if( kfl_usrbc_chm(iclas) /= 0 ) then
           call chm_usrbcs(iclas,kfl_usrbc_chm(iclas))
        end if
     end do
 
     !-------------------------------------------------------------
     !
     ! KFL_ROBIN_CHM: Count number of Robin boundaries
     !
     !-------------------------------------------------------------

     kfl_robin_chm = 0
     iboun         = 0
     do while( iboun < nboun )
        iboun = iboun+1
        iclas = 0
        do while( iclas < nclas_chm )
           iclas = iclas+1
           if( kfl_fixbo_chm(iclas,iboun) == 2 ) then
              kfl_robin_chm = kfl_robin_chm+ 1 
              iclas         = nclas_chm
           end if
        end do
     end do

  end if
!

end subroutine chm_inibcs

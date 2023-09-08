subroutine rad_inibcs()
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_inibcs
  ! NAME
  !    rad_inibcs
  ! DESCRIPTION
  !    This routine applied boundary conditions
  ! OUTPUT 
  ! USES
  ! USED BY
  !    rad_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_radiat
  use mod_opebcs
  implicit none
  integer(ip)  :: ipoin,pnodb,iboun,inodb,ifunc,ipara,ibsta,knodb(mnodb)
  integer(ip)  :: pblty,ncodf,nbcod

  if( INOTMASTER ) then

     !-------------------------------------------------------------
     !
     ! Allocate memory
     !
     !-------------------------------------------------------------

     call rad_membcs(1_ip)
     if( kfl_conbc_rad == 0 ) then
        call rad_membcs(2_ip)
     end if

     !-------------------------------------------------------------
     !
     ! Node codes
     !
     !-------------------------------------------------------------

     if( kfl_icodn > 0 ) then

        if(kfl_conbc_rad==0) then
           iffun     =  1
           kfl_funno => kfl_funno_rad
        else
           iffun      =  0
        end if
        ifloc     =  0
        ifbop     =  0
        kfl_fixno => kfl_fixno_rad
        bvess     => bvess_rad
        tncod     => tncod_rad
        call reacod(10_ip)
        
     end if

     !-------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------

     if( kfl_icodb > 0 ) then

        kfl_fixbo => kfl_fixbo_rad
        bvnat     => bvnat_rad(:,:,1)
        tbcod     => tbcod_rad
        call reacod(20_ip)
        
     end if

     call rad_bcntoe() 

     !-------------------------------------------------------------
     !
     ! Put -1 value to 0
     !
     !-------------------------------------------------------------

     do ipoin=1,npoin
        if(kfl_fixno_rad(1,ipoin)==-1) kfl_fixno_rad(1,ipoin)=0
     end do

  end if

end subroutine rad_inibcs

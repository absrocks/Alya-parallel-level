subroutine ale_inibcs()
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_inibcs
  ! NAME
  !    ale_inibcs
  ! DESCRIPTION
  !    This routine applied boundary conditions
  ! OUTPUT 
  ! USES
  ! USED BY
  !    ale_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_alefor
  use mod_opebcs
  use mod_memory,         only :  memory_alloca_min

  implicit none
  integer(ip)  :: ipoin,iboun,inodb,idime

  if( INOTMASTER ) then

     !-------------------------------------------------------------
     !
     ! Allocate memory
     !
     !-------------------------------------------------------------

     call ale_membcs(1_ip)
     if( kfl_conbc_ale == 0 ) then
        call ale_membcs(45_ip)
        call ale_membcs(48_ip)
     end if

     !-------------------------------------------------------------
     !
     ! Default: All boundary nodes should be prescribed
     !
     !-------------------------------------------------------------

     !do ipoin = 1,npoin
     !   if( lpoty(ipoin) /= 0 ) then
     !      do idime = 1,ndime
     !         kfl_fixno_ale(idime,ipoin) = 1
     !      end do
     !   end if
     !end do

     !-------------------------------------------------------------
     !
     ! Node codes
     !
     !-------------------------------------------------------------

     if( kfl_icodn > 0 ) then

        if( kfl_conbc_ale == 0 ) then
           iffun     =  1
           kfl_funno => kfl_funno_ale
        else
           iffun      =  0
        end if
        ifloc     =  0
        ifbop     =  0
        kfl_fixno => kfl_fixno_ale
        bvess     => bvess_ale
        tncod     => tncod_ale
        call reacod(10_ip)
        
     end if

     !-------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------

     if( kfl_icodb > 0 ) then

        !kfl_fixbo => kfl_fixbo_ale
        !tbcod     => tbcod_ale
        !call reacod(20_ip)
        
     end if
     
     !-------------------------------------------------------------
     !
     ! Boundary to node fixity and function
     !
     !-------------------------------------------------------------

     if( kfl_conbc_ale == 0 ) then
        call memgen(1_ip,npoin,0_ip)
        do iboun = 1,nboun
           if( kfl_fixbo_ale(iboun) == 1 ) then
              do inodb = 1,nnode(ltypb(iboun))
                 ipoin = lnodb(inodb,iboun)
                 kfl_fixno_ale(1,ipoin) = 1
                 kfl_funno_ale(ipoin)   = kfl_funbo_ale(iboun)
              end do
           end if
        end do
        call parari('SLX',NPOIN_TYPE,npoin,gisca)
        do ipoin = 1,npoin
           if( gisca(ipoin) >= 1 ) then
              do idime = 1,ndime
                 kfl_fixno_ale(idime,ipoin) = 1
              end do
           end if
        end do
        call memgen(3_ip,npoin,0_ip)
     end if 

     !-------------------------------------------------------------
     !
     ! If one dof is prescribed presribed the other
     !
     !-------------------------------------------------------------

     !do ipoin = 1,npoin
     !   do idime = 2,ndime
     !      kfl_fixno_ale(idime,ipoin) = kfl_fixno_ale(1,ipoin) 
     !   end do
     !end do

     !-------------------------------------------------------------
     !
     ! Support geometry
     !
     !-------------------------------------------------------------

     if( kfl_suppo_ale == 1 ) then
        call ale_suppor()
     end if 

  else

     call memory_alloca_min(bvess_ale) 

  end if

end subroutine ale_inibcs

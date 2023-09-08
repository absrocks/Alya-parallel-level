subroutine Conblk()
  !-----------------------------------------------------------------------
  !****f* master/Conblk
  ! NAME
  !    Alya
  ! DESCRIPTION
  !    Increase block number
  ! USES
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use mod_ker_proper
  use def_domain
  use def_coupli, only : mcoup
  use def_coupli, only : kfl_gozon
  use def_coupli, only : coupling_driver_iteration
  use def_coupli, only : coupling_driver_number_couplings
  use mod_commdom_alya,     only: INONE
  use mod_messages, only : livinf
  implicit none
  !
  ! Write message if block is in a coupling loop
  !
  if( mcoup > 0 ) then
     if( coupling_driver_number_couplings(iblok) /= 0 ) then
        call livinf(-13_ip,'END ZONAL COUPLING: ',coupling_driver_iteration(iblok))
     end if
  end if
  !
#ifdef COMMDOM
  call moduls(ITASK_CONBLK)
#endif 
  ! 
  ! Initialize
  !   
  coupling_driver_iteration(iblok) = 0
  kfl_gozon = 1
  iblok     = iblok+1  
  kfl_gocou = 1
  itcou     = 1
  if( iblok > nblok ) then
     kfl_goblk = 0
  end if
  if( nblok > 1 ) then
     call livinf(eight,' ',zero)
  end if

  call ker_updpro(ITASK_CONBLK) ! update property

end subroutine Conblk

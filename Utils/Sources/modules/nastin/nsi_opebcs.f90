subroutine nsi_opebcs()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_opebcs
  ! NAME 
  !    nsi_opebcs
  ! DESCRIPTION
  !    This routine imposes boiundary conditions according to codes
  ! USES
  ! USED BY
  !    nsi_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_nastin 
  use mod_memchk
  implicit none

  !if( INOTMASTER .and. kfl_wwork == 1 ) then

  !if( kfl_conbc_nsi == 0 ) then
  !   iffun      =  1
  !   kfl_funno  => kfl_funno_nsi
  !else
  !   iffun      =  0
  !end if
  !ifbop     =  0
  !ifloc     =  1
  !skcos_nsi => skcos
  !kfl_fixrs => kfl_fixrs_nsi
  !kfl_fixno => kfl_fixno_nsi
  !bvess     => bvess_nsi(:,:,1)
  !call reacod(1_ip)
  
end subroutine nsi_opebcs

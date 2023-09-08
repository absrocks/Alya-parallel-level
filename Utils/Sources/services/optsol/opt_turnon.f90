subroutine opt_turnon()
  !-----------------------------------------------------------------------
  !****f* Optsol/opt_turnon
  ! NAME 
  !    opt_turnon
  ! DESCRIPTION
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_optsol
  implicit none
  !
  ! Exchange variables between Master and Slaves
  !
  call opt_parall(-1_ip)
  !
  ! Allocate memory (it needs opt_parall(-1) before) 
  !
  call opt_memunk()
  !
  ! Exchange variables between Master and Slaves
  !
  call opt_parall(1_ip)


end subroutine opt_turnon

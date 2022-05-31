subroutine ibm_timste()
  !-----------------------------------------------------------------------
  !****f* Immbou/ibm_timste
  ! NAME 
  !    ibm_timste
  ! DESCRIPTION
  !    This routine computes the time step
  ! USES
  !    ibm_iniunk
  !    ibm_updtss
  !    ibm_updbcs
  !    ibm_updunk
  ! USED BY
  !    Nastin
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  implicit none
  !
  ! Time step size
  !
  call ibm_updtss()

end subroutine ibm_timste

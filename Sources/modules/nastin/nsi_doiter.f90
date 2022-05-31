subroutine nsi_doiter()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_doiter
  ! NAME 
  !    nsi_doiter
  ! DESCRIPTION
  !    This routine solves an iteration of the linearized incompressible NS
  !    equations.
  ! USES
  !    nsi_begite
  !    nsi_updbcs
  !    nsi_solite
  ! USED BY
  !    Nastin
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_domain
  use def_nastin
  implicit none
  
  if( kfl_stead_nsi == 0 ) then
     call nsi_begite()

     do while( kfl_goite_nsi == 1 )
        call nsi_solite()
        call nsi_endite(ITASK_ENDINN)
     end do

     call nsi_endite(ITASK_ENDITE)

  end if

end subroutine nsi_doiter


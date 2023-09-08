subroutine Optsol(order)
  !-----------------------------------------------------------------------
  !****f* optsol/Optsol
  ! NAME 
  !    Optsol
  ! DESCRIPTION
  !    This routine is the bridge for Optsol service
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_optsol
  use def_domain
  use def_elmtyp
  use def_inpout
  !use mod_memchk
  implicit none
  integer(ip), intent(in) :: order
 
  if(order==ITASK_REAPRO) then
    ! Read service data 
    call opt_reapro()          
  else
     if( kfl_servi(ID_OPTSOL) == 1 ) then 

        if(order==ITASK_TURNON) then
           ! Allocate memory for several arrays and
           ! send parameters and variables between slaves and master 
           call opt_turnon()
        elseif(order==ITASK_INIUNK) then
           ! Nothing
        elseif(order==ITASK_BEGSTE) then 
           ! Nothing
        elseif(order==ITASK_ENDSTE) then 
           ! Nothing
        elseif(order==ITASK_DOOPTI) then
           ! Calculate descent direction and update candidate
           call opt_doopti()
        elseif(order==ITASK_ENDOPT) then
           ! Accept/reject update candidate, update step length, update control flags and variables
           call opt_endopt() 
        elseif(order==ITASK_TURNOF) then
           ! Nothing
        else
           ! Not implemented
        end if
     end if
  end if

end subroutine Optsol


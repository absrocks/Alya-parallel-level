subroutine pts_lagdef(itask,numpr)
  !-----------------------------------------------------------------------
  !****f* Partis/pts_lagdef
  ! NAME
  !    pts_lagdef
  ! DESCRIPTION
  !    This subroutine deos the following:
  !    ITASK = 0 ... Initialize the solver type
  !    ITASK = 1 ... Bridge between modules and parall service
  !    ITASK = 4 ... Define matrices, RHS and preconditioner sizes
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master 
  use def_kermod
  use def_partis
  use def_solver
  use def_inpout
  use mod_memchk
  implicit none  
  integer(ip), intent(in)    :: itask
  integer(ip), intent(inout) :: numpr
  integer(ip)                :: ilagr,ii,mlagr_tmp,iprop
  integer(ip)                :: ilag1,ilag2
  integer(4)                 :: istat

  if( itask == 3_ip ) then

     !-------------------------------------------------------------------
     !
     ! Modules request property slots, return initial slot assigned
     !
     !-------------------------------------------------------------------

     if( numpr > 0 ) then
        if (nlapr+numpr <= mlapr ) then
           nlapr = nlapr + numpr ! Increase number of assigned slots          
           numpr = nlapr - numpr + 1 ! Return the first assigned slot in the properties vector
        else
           call runend("LAGDEF: Requested too many properties for lagrangian particles")
        endif
     else
        call runend("LAGDEF called without number of memory slots requested")
     endif

  end if

end subroutine pts_lagdef

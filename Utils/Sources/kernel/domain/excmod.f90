subroutine excmod(itask)
  !-----------------------------------------------------------------------
  !****f* domain/excmod
  ! NAME
  !    excmod
  ! DESCRIPTION
  !    This routine does some exceptional things with modules
  ! USED BY
  !    Turnon
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  implicit none
  integer(ip), intent(in) :: itask
  
  select case ( itask )
     
  case ( 1_ip ) 
     !
     ! Nothing to do
     !

  case ( 2_ip ) 
     !
     ! Immbou module data should be read here
     !
     modul = ID_IMMBOU
     call moduls(-ITASK_TURNON) 
   
  end select
  
end subroutine excmod

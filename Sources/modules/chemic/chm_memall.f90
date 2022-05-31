subroutine chm_memall()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_memall
  ! NAME 
  !    chm_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use mod_chm_arrays, only : chm_arrays

  implicit none
  !
  ! Allocate 
  !
  call chm_arrays('ALLOCATE')

end subroutine chm_memall


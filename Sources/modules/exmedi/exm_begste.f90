subroutine exm_begste
!-----------------------------------------------------------------------
!****f* Exmedi/exm_begste
! NAME 
!    exm_begste
! DESCRIPTION
!    This routine prepares for a new time step 
! USES
!    exm_iniunk
!    exm_updtss
!    exm_updbcs
!    exm_updunk
! USED BY
!    Exmedi
!***
!-----------------------------------------------------------------------

  use def_parame
  use def_master
  use def_domain

  use def_exmedi

  implicit none
  
  !
  ! Initial guess fo the unknowns: U(n,0,*) <-- U(n-1,*,*).
  !
  call exm_updunk(ITASK_BEGSTE)     ! u(,ITER_AUX) <-- u(,TIME_N) 
  call exm_upcell(ITASK_BEGSTE)     !  (,ITER_AUX) <--  (,TIME_N) 

end subroutine exm_begste

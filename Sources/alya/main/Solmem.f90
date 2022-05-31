!-----------------------------------------------------------------------
!> @addtogroup Solmem
!> @{
!> @file    Solmem.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Solver memory
!> @details Compute solver requirements
!> @} 
!-----------------------------------------------------------------------

subroutine Solmem()

  use def_kintyp,             only : ip
  use def_kintyp,             only : lg
  use def_master,             only : iblok,kfl_paral,ittim
  use def_master,             only : nblok
  use def_master,             only : modul
  use def_master,             only : mmodu
  use def_master,             only : momod
  use def_master,             only : kfl_modul
  use def_master,             only : ITASK_SOLMEM
  use mod_output_postprocess, only : output_postprocess_allocate_sets_and_witness
  use mod_moduls,             only : moduls

  implicit none

  call Kermod(ITASK_SOLMEM)
  do iblok = 1,nblok
    call moduls(ITASK_SOLMEM)
  end do
  !
  ! Postprocess NOT GOOD, CHECK ORDER OF THINGS... put open file in turnon after moddef...
  !
  do modul = 1,mmodu
     if( kfl_modul(modul) /= 0 ) then
        call output_postprocess_allocate_sets_and_witness(momod(modul) % postp(1))
     end if
  end do
  modul = 0
  !
  ! Allocate memory for all unknowns of the problem and coefficients
  !
  call memunk(1_ip)

end subroutine Solmem

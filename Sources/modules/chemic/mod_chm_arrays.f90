!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    mod_chm_arrays.f90
!> @author  houzeaux
!> @date    2019-11-16
!> @brief   Nastin arrays
!> @details Nastin arrays
!-----------------------------------------------------------------------

module mod_chm_arrays

  use def_master
  use def_domain 
  use def_chemic
  use def_kermod
  use def_solver
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  use mod_arrays,              only : arrays
  use mod_arrays,              only : arrays_number
  use mod_memory,              only : memory_alloca
  use mod_ADR,                 only : ADR_initialize_type
  use mod_ADR,                 only : ADR_check_and_compute_data
  use mod_ADR,                 only : ADR_allocate_projections_bubble_sgs
  use mod_ADR,                 only : ADR_arrays

  implicit none

  private

  public :: chm_arrays
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-16
  !> @brief   Nastin arrays
  !> @details Do what you have to do with nastin arrays
  !> 
  !-----------------------------------------------------------------------
  
  subroutine chm_arrays(wtask)

    character(len=*), intent(in) :: wtask
    integer(ip)                  :: iclas
    ! 
    ! CONCE
    !
    call arrays(arrays_number('CONCE'),wtask,conce,npoin,nclas_chm,ncomp_chm)
    if( wtask == 'ALLOCATE' .and. npoin == 0 ) then
       call memory_alloca(mem_modul(1:2,modul),'CONCE','chm_memall',conce,1_ip,nclas_chm,ncomp_chm)      
    end if
    !
    ! ADR type
    !
    do iclas=1,nclas_chm
       call ADR_check_and_compute_data(ADR_chm(iclas))
       !call ADR_arrays(wtask,ADR_chm(iclas),'BUBBT','PROJ1','PROJ2','TESG2',TAG1=iclas)       
    end do
    
  end subroutine chm_arrays
   
end module mod_chm_arrays
!> @}

!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    mod_exm_arrays.f90
!> @author  houzeaux
!> @date    2019-11-16
!> @brief   Nastin arrays
!> @details Nastin arrays
!-----------------------------------------------------------------------

module mod_exm_arrays

  use def_master
  use def_domain 
  use def_exmedi
  use def_kermod
  use mod_arrays, only : arrays
  use mod_arrays, only : arrays_number
  use mod_memory, only : memory_alloca

  
  implicit none

  private

  public :: exm_arrays
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-16
  !> @brief   Nastin arrays
  !> @details Do what you have to do with nastin arrays
  !> 
  !-----------------------------------------------------------------------
  
  subroutine exm_arrays(wtask)

    character(len=*), intent(in) :: wtask

    call arrays(arrays_number('TAULO'),trim(wtask),taulo,npoin)
    call arrays(arrays_number('FISOC'),trim(wtask),fisoc,npoin)
    call arrays(arrays_number('KWAVE'),trim(wtask),kwave_exm,npoin)  
    call arrays(arrays_number('ISOCC'),trim(wtask),isoch_modified,npoin) !ISOCH
    call arrays(arrays_number('ELMAG'),trim(wtask),elmag,npoin,ncomp_exm)
    call arrays(arrays_number('VCONC'),trim(wtask),vconc,nconc_exm,npoin,ncomp_exm)
    call arrays(arrays_number('QNETT'),trim(wtask),qneto_exm,npoin)      !QNETO
    call arrays(arrays_number('REFHN'),trim(wtask),refhn_exm,npoin,ncomp_exm)
    call arrays(arrays_number('VAUXI'),trim(wtask),vauxi_exm,nauxi_exm,npoin,ncomp_exm)
    call arrays(arrays_number('VICEL'),trim(wtask),vicel_exm,nicel_exm,npoin,1_ip)    
    
  end subroutine exm_arrays

end module mod_exm_arrays

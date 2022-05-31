!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_ker_arrays.f90
!> @author  houzeaux
!> @date    2019-11-16
!> @brief   Nastin arrays
!> @details Nastin arrays
!-----------------------------------------------------------------------

module mod_ker_arrays

  use def_master
  use def_domain 
  use def_kermod
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  use mod_arrays,              only : arrays
  use mod_arrays,              only : arrays_number
  use mod_memory,              only : memory_alloca
  use mod_communications,      only : PAR_MAX
  implicit none

  private

  public :: ker_arrays
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-16
  !> @brief   Nastin arrays
  !> @details Do what you have to do with nastin arrays
  !> 
  !-----------------------------------------------------------------------
  
  subroutine ker_arrays(wtask)

    character(len=*), intent(in) :: wtask
    integer(ip)                  :: nclas,ifunc,kfl_value
    !
    ! Time-averaged velocity for wall law
    !
    if( kfl_wlaav_ker /= 0 ) then
       call arrays(75_ip,wtask,velav_ker,ndime,mgaub,nboun)
       call arrays(76_ip,wtask,avupo_ker,ndime,npoin)
    end if
    !
    ! Wall law
    !
    if( kfl_noslw_ker /= 0 ) then
       call arrays(77_ip,wtask,avta1_nsw_ker,ndime,nelem)
       call arrays(78_ip,wtask,fact_nsw_ker,npoin)
    end if
    !
    ! Advection vector
    ! KFL_VEFUN = 0 ... ADVEC points to VELOC 
    !           < 0 ... Copy values from a field
    !           > 0 ... User-defined function
    !
    if( kfl_vefun == 0 ) then
       advec => veloc
    else
       call arrays(arrays_number('ADVEC'),wtask,advec,ndime,npoin,3_ip)
    end if
    !
    ! Temperature
    ! KFL_TEFUN = 0 ... Do not do anything
    !           < 0 ... Copy values from a field
    !           > 0 ... User-defined function
    !
    if( kfl_tefun == 0 ) then
       continue
    else
       call arrays(arrays_number('TEMPE'),wtask,tempe,npoin,3_ip)
       therm => tempe
    end if
    !
    ! Ddisplacement and mesh velocity
    ! KFL_DIFUN = 0 ... Do not do anything
    !           < 0 ... Copy values from a field
    !           > 0 ... User-defined function
    ! 
    if( kfl_difun == 0 ) then
       continue
    else
       call arrays(arrays_number('DISPM'),wtask,dispm,ndime,npoin,3_ip)
    end if
    !
    ! Advection vector
    ! KFL_COFUN = 0 ... Do not do anything
    !           < 0 ... Copy values from a field
    !           > 0 ... User-defined function
    !
    if( kfl_cofun == 0 ) then
       continue
    else
       !
       ! Guess number of classes
       !
       nclas = 4_ip 
       if( kfl_cofun < 0 ) then
          kfl_value = -kfl_cofun
          if( associated(xfiel(kfl_value) % a) ) nclas = size(xfiel(kfl_value) % a(:,1,1), KIND=ip)
          call PAR_MAX(nclas)
       else
           if( kfl_cofun > 1000 ) then 
               !
               ! Use function dimension to determine size
               !
               ifunc = kfl_cofun - 1000
               nclas = size(space_time_function(ifunc) % expression, KIND=ip)
           endif
           if( kfl_cofun == 666_ip) then
               nclas = 1_ip 
           end if          
       end if
       call arrays(arrays_number('CONCE'),wtask,conce,npoin,nclas,3_ip)
    end if
    
  end subroutine ker_arrays
   
end module mod_ker_arrays
!> @}

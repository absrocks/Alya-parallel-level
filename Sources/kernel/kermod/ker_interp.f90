!-----------------------------------------------------------------------
!> @addtogroup kermod
!> @{
!> @file    ker_interp.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Array redistribution
!> @details Array redistribution and reallocation in case values are
!>          not needed
!> @} 
!-----------------------------------------------------------------------

subroutine ker_interp()
  
  use def_master
  use def_domain
  use def_kermod
  use mod_parall,        only : par_memor
  use mod_ker_proper,    only : NUMBER_OF_PROPERTIES
  use mod_ker_proper,    only : ker_proper_property_array_exists
  use mod_ker_proper,    only : ker_proper_point_to_property
  use mod_ker_arrays,    only : ker_arrays
  use def_AMR
  use mod_AMR_interpolate
  implicit none
  integer(ip)                   :: iprop
  type(typ_valpr_ker), pointer  :: prope_ker
  real(rp),            pointer  :: xx(:,:)

  nullify(xx)
  !
  ! Boundary conditions
  !  
  call AMR_interpolate(walld,              'NPOIN',MEMOR=par_memor,VARIABLE_NAME='WALLD')
  call AMR_interpolate(kfl_fixno_walld_ker,'NPOIN',MEMOR=par_memor,VARIABLE_NAME='KFL_FIXNO_WALLD_KER')
  call AMR_interpolate(walln,              'NPOIN',MEMOR=par_memor,VARIABLE_NAME='WALLN')
  call AMR_interpolate(kfl_fixno_walln_ker,'NPOIN',MEMOR=par_memor,VARIABLE_NAME='KFL_FIXNO_WALLN_KER')
  !
  ! Arrays
  !
  call ker_arrays('INTERPOLATE')
  !
  ! Properties
  !  
  do iprop = 1,NUMBER_OF_PROPERTIES
     call ker_proper_point_to_property(iprop,prope_ker)
     
!!$     if( prope_ker % kfl_exist == 1 ) then
!!$        if( ker_proper_property_array_exists('VALUE_IPOIN',prope_ker) ) &
!!$             call AMR_interpolate(prope_ker % value_ipoin,    'NPOIN',MEMOR=par_memor,VARIABLE_NAME='PROPE_KER % VALUE_IPOIN')
!!$        if( ker_proper_property_array_exists('GRVAL_IPOIN',prope_ker) ) &
!!$             call AMR_interpolate(prope_ker % grval_ipoin,    'NELEM',MEMOR=par_memor,VARIABLE_NAME='PROPE_KER % GRVAL_IPOIN')
!!$        if( ker_proper_property_array_exists('VALUE_IELEM',prope_ker) ) &
!!$             call AMR_interpolate(prope_ker % value_ielem,    'NELEM',MEMOR=par_memor,VARIABLE_NAME='PROPE_KER % VALUE_IELEM')
!!$        if( ker_proper_property_array_exists('VALUE_IBOUN',prope_ker) ) &
!!$             call AMR_interpolate(prope_ker % value_iboun,    'NBOUN',MEMOR=par_memor,VARIABLE_NAME='PROPE_KER % VALUE_IBOUN')
!!$        if( ker_proper_property_array_exists('GRVAL_IELEM',prope_ker) ) &
!!$             call AMR_interpolate(prope_ker % grval_ielem,    'NELEM',MEMOR=par_memor,VARIABLE_NAME='PROPE_KER % GRVAL_IELEM')
!!$        if( ker_proper_property_array_exists('DRVAL_IELEM',prope_ker) ) &
!!$             call AMR_interpolate(prope_ker % drval_ielem,    'NELEM',MEMOR=par_memor,VARIABLE_NAME='PROPE_KER % DRVAL_IELEM')
!!$        if( ker_proper_property_array_exists('DRVAL_VEL_IELEM',prope_ker) ) &
!!$             call AMR_interpolate(prope_ker % drval_vel_ielem,'NELEM',MEMOR=par_memor,VARIABLE_NAME='PROPE_KER % DRVAL_VEL_IELEM')
!!$        if( ker_proper_property_array_exists('DRVAL_TUR_IELEM',prope_ker) ) &
!!$             call AMR_interpolate(prope_ker % drval_tur_ielem,'NELEM',MEMOR=par_memor,VARIABLE_NAME='PROPE_KER % DRVAL_TUR_IELEM')
!!$        if( ker_proper_property_array_exists('GDVAL_IELEM',prope_ker) ) &
!!$             call AMR_interpolate(prope_ker % gdval_ielem,    'NELEM',MEMOR=par_memor,VARIABLE_NAME='PROPE_KER % GDVAL_IELEM')
!!$     end if
  end do

end subroutine ker_interp


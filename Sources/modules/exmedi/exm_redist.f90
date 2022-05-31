!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_redist.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Redistribution of arrays
!> @details Redistribution of arrays. mainly allocated exm_memall
!> @} 
!-----------------------------------------------------------------------

subroutine exm_redist()

  use def_master
  use def_domain
  use def_kermod
  use mod_redistribution
  use mod_parall
  use def_parall
  use mod_memory
  use mod_mesh_type
  use def_exmedi
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_exm_arrays,         only : exm_arrays
  implicit none
  !
  ! Variables in exm_membcs: boundary conditions
  !
  call redistribution_array(kfl_fixno_exm,'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXNO_EXM')
  call redistribution_array(kfl_fixbo_exm,'NBOUN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXBO_EXM')
  call redistribution_array(bvess_exm,    'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVESS_EXM')
  call redistribution_array(bvnat_exm,    'NBOUN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVNAT_EXM')
  !
  ! Variables in memall
  !
  call exm_arrays('REDISTRIBUTE')

end subroutine exm_redist

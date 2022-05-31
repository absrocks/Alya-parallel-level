!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    nsi_begste.f90
!> @author  Guillaume Houzeaux
!> @date    05/06/2013
!> @brief   Begin a time step
!> @details Begin a time step
!> @} 
!-----------------------------------------------------------------------
subroutine ker_begste()
  
  use def_master,        only : ITASK_BEGSTE
  use mod_mass_matrix,   only : mass_matrix_consistent
  use mod_ker_subdomain, only : ker_subdomain_update_coordinates

  implicit none
  !
  ! Transient fields
  !
  call calc_kx_tran_fiel()
  !
  ! Update fields: Velocity, temperature, concentration and displacement functions
  !
  call ker_velfun(ITASK_BEGSTE)
  call ker_temfun(ITASK_BEGSTE)
  call ker_confun(ITASK_BEGSTE)
  call ker_disfun(ITASK_BEGSTE)
  !
  ! Update subdomains
  !
  call ker_subdomain_update_coordinates()  
  !
  ! Weighted consistent mass... Kermod should have been read
  !
  call mass_matrix_consistent(CONSISTENT_WEIGHTED_MASS=.true.)

end subroutine ker_begste

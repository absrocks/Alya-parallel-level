  !----------------------------------------------------------------------
  !> @addtogroup Solidz
  !> @{
  !> @file    ale_plugin.f90
  !> @author  J.C. Cajas
  !> @date    04/11/2014
  !> @brief   Plugin for zonal coupling with other modules
  !> @details Plugin for zonal coupling with other modules. 
  !> @        Used in the FSI coupling with Nastin and Solidz.
  !> @} 
  !----------------------------------------------------------------------

subroutine ale_plugin(icoup)
  !
  ! Obligatory variables 
  !
  use def_coupli,        only :  coupling_type
  use def_domain,        only :  npoin
  use def_domain,        only :  ndime
  use def_master,        only :  solve_sol,modul
  use def_kintyp,        only :  ip,rp
  use def_master,        only :  current_code
  use def_master,        only :  INOTMASTER
  use mod_couplings,     only :  COU_INTERPOLATE_NODAL_VALUES
  use mod_couplings,     only :  COU_SET_FIXITY_ON_TARGET
  use mod_memory,        only :  memory_deallo
  use mod_memory,        only :  memory_alloca
  use mod_matrix,        only :  matrix_initialize
  !
  ! Possible variables 
  !
  use def_master,        only :  displ,kfl_fixno_ale
  use def_master,        only :  bvess_ale
  use mod_parall,        only :  PAR_GLOBAL_TO_LOCAL_NODE

  implicit none

  real(rp),    pointer    :: svalu(:,:)

  integer(ip)             :: idime,ipoin,kpoin
  integer(ip), intent(in) :: icoup
  character(5)            :: variable
  integer(ip), save       :: ipass = 0

  nullify(svalu)
  variable = coupling_type(icoup) % variable
  
   if( variable == 'ALEFO' ) then
     !
     ! Coupling with solidz
     !
     if ( INOTMASTER ) then
        
        allocate(svalu(ndime, npoin))
        
     else
        
        allocate(svalu(1_ip,1_ip))
        
     end if
     if( ipass == 0 ) then
        ipass = 1
        call COU_SET_FIXITY_ON_TARGET('ALEFO',modul,kfl_fixno_ale)
     end if     
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,bvess_ale,svalu,kfl_fixno_ale)

  end if

  if( associated(svalu) ) deallocate( svalu )
  
end subroutine ale_plugin
!> @} 
!-----------------------------------------------------------------------

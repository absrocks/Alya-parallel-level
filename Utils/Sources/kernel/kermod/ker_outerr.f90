!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    ker_outerr.f90
!> @author  Guillaume Houzeaux
!> @date    19/02/2016
!> @brief   Check errors
!> @details Check errors
!> @} 
!-----------------------------------------------------------------------

subroutine ker_outerr()

  use def_master
  use def_kermod
  use def_domain
  use def_solver
  use mod_communications, only : PAR_SUM
  use mod_outfor, only : outfor
  implicit none
  integer(ip) :: ierro,iwarn

  ierro = 0
  iwarn = 0
  !
  ! Direct solver
  ! 
  if( kfl_direct_solver == SOL_DIRECT_SOLVER_PASTIX ) then
#ifndef PASTIX 
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'COMPILE ALYA WITH -DPASTIX AND LINK IT WITH PASTIX')     
#endif
  end if
  if( kfl_direct_solver == SOL_DIRECT_SOLVER_MUMPS ) then
#ifndef MUMPS
     ierro = ierro + 1 
     call outfor(1_ip,momod(modul)%lun_outpu,'COMPILE ALYA WITH -DMUMPS AND LINK IT WITH MUMPS')     
#endif
  end if
  if( kfl_direct_solver == SOL_DIRECT_SOLVER_WSMP ) then
#ifndef WSMP
     ierro = ierro + 1
     print*,'ierro',ierro
     call outfor(1_ip,momod(modul)%lun_outpu,'COMPILE ALYA WITH -DWSMP AND LINK IT WITH WSMP')
#endif
  end if
  if( kfl_direct_solver == SOL_DIRECT_SOLVER_PWSMP) then
#ifndef PWSMP
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'COMPILE ALYA WITH -DPWSMP AND LINK IT WITH PWSMP')
#endif
  end if

!  if( kfl_noslw_ker /= 0_ip .and. abs(delmu_dom-2.0_rp) > 1.0e-8_rp) then
!     ierro = ierro + 1
!     call outfor(1_ip,momod(modul)%lun_outpu,'When no slip wall is used MULTI must be 2 for wall law')  ! else ywalb and thus ywale are incorrect
!  end if

  

#ifdef VECTOR_SIZE_VARIABLE
  if( VECTOR_SIZE <= 0 .and. VECTOR_SIZE /= -1 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'VECTOR SIZE SHOULD BE DEFINED > 0 or =-1')
  end if
#endif

  !
  ! Periodicity does not work with SFC 
  !
  if( nperi/=0_ip .and. kfl_ngrou /= -1_ip  ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'FOR PERIODIC CASES YOU NEED TO USE GROUPS = int, SEQUENTIAL_FRONTAL IN THE DOM.DAT')
  end if

  !
  ! Exchange location needs elsest 
  !
  if( kfl_waexl_ker/=0_ip .and. kfl_elses == 0_ip  ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'EXCHANGE LOCATION NEEDS ELSEST')
  end if

  

  call PAR_SUM(ierro,'IN MY CODE','INCLUDE MASTER')
  if( ierro /= 0 ) call runend('AN ERROR HAS BEEN FOUND IN KERMOD')

end subroutine ker_outerr

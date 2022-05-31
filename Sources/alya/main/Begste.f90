!-----------------------------------------------------------------------
!> @addtogroup Begste
!> @{
!> @file    Begste.f90
!> @author  Guillaume Houzeaux
!> @brief   Begin a time step
!> @details Begin a time step:
!>          - Modules update boundary conditions, unknowns, etc.
!> @} 
!-----------------------------------------------------------------------
subroutine Begste()
  use def_kintyp,        only : ip
  use def_master
  use def_parame
  use mod_ker_proper
  use def_domain,        only : kfl_domar_world
  use def_coupli,        only : mcoup
  use def_coupli,        only : coupling_driver_iteration
  use mod_couplings,     only : COU_TEMPORAL_PREDICTOR
  use mod_messages,      only : livinf
  use mod_moduls,        only : moduls 
  use mod_ker_subdomain, only : ker_subdomain_update_coordinates
  implicit none 

  integer(ip)            :: icoup
  !
  ! Turn back reset flag to previous step
  !
  if (kfl_reset == 1) then
     call iniste(2_ip)
     cutim  = cutim - dtime
     call setgts(ITASK_TIMSTE)     
     call livinf(201_ip, ' ',1_ip)
  endif
  !
  ! Coupling
  !
  call cou_begste()
  !
  ! Begin a time step for each module
  !
  call Kermod(ITASK_BEGSTE)
  do iblok = 1,nblok
     call moduls(ITASK_BEGSTE)
  end do
  iblok = 1
  !
  ! Turn back properties if reset
  !
  if (kfl_reset == 1) then
     call ker_updpro()
     ! Deactivate reset request
     kfl_reset = 0
  end if
  !
  ! Coupling: Put counters to zero
  ! 
  if( mcoup > 0 ) then
     coupling_driver_iteration(1:nblok) = 0
  end if
  !
  ! Moving meshes, recompute some things
  !
  if( kfl_domar_world == 1 ) call domarr(3_ip)

  contains
!-------------------------------------------------------------------------||---!
!                                                                              !
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
    subroutine coupling_set_dt()
      use def_master, only : ITASK_TIMSTE
      use mod_communications, only : par_min, par_max  
      use mod_couplings,      only : THERE_EXISTS_A_ZONE_COUPLING
      use def_coupli,         only : mcoup, coudt
      implicit none 
      real(rp) :: dt_inv(1_ip)
      !
      !
      !print *, coudt, "<-----"
      !
      if( (mcoup > 0).and.(coudt==1_ip).and.THERE_EXISTS_A_ZONE_COUPLING() ) then 
         dt_inv(1_ip) = dtinv  
         call PAR_MAX(1_ip, dt_inv, 'IN CURRENT COUPLING') !< 2016Mar23 
         !call PAR_MIN(1_ip, dt_inv, 'IN CURRENT COUPLING') 
         !
         !    if(imaster) print *, 1.0/dtinv, 1.0/dt_inv(1_ip), "<----"
         !
         dtinv = dt_inv(1_ip)      
         !
         cutim  = cutim - dtime          ! 
         !call iniste(2_ip)               ! 
         call setgts(ITASK_TIMSTE) 
         call livinf(201_ip, ' ',1_ip)   ! 
         !
      endif
      ! 

    end subroutine coupling_set_dt
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!

  end subroutine Begste

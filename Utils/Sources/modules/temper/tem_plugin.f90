!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_plugin.f90
!> @author  Guillaume Houzeaux
!> @date    13/04/2014
!> @brief   Plugin with coupling
!> @details Plugin for the zonal coupling
!> @} 
!-----------------------------------------------------------------------

subroutine tem_plugin(icoup)
  !
  ! Obligatory variables 
  !
  use def_kintyp,         only :  ip,rp
  use def_coupli,         only :  coupling_type
  use def_master,         only :  solve_sol,modul,lninv_loc
  use mod_couplings,      only :  COU_INTERPOLATE_NODAL_VALUES
  use mod_couplings,      only :  COU_SET_FIXITY_ON_TARGET
  use mod_memory,         only :  memory_deallo
  use mod_memory,         only :  memory_alloca
  use mod_matrix,         only :  matrix_initialize
  !
  ! Possible variables 
  !
  use def_master,         only :  tempe, therm
  use def_temper,         only :  bvess_tem
  use mod_parall,         only :  PAR_GLOBAL_TO_LOCAL_NODE

  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use def_master,         only : current_code,inotmaster,kfl_paral
  use mod_parall
  use def_temper,         only : kfl_fixno_tem
  use def_domain
  implicit none 
  integer(ip), intent(in) :: icoup
  character(5)            :: variable
  integer(ip)             :: ipoin
  integer(ip), save       :: ipass=0

  variable = coupling_type(icoup) % variable 

  if( variable == 'TEMPE' ) then   
     !
     ! Temperature 
     !
     if( ipass == 0 ) then
        ipass = 1
        call COU_SET_FIXITY_ON_TARGET('TEMPE',modul,kfl_fixno_tem)
     end if
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,bvess_tem,therm,kfl_fixno_tem)      
     
  else if( variable == 'RESID' ) then
     !
     ! Residual
     !
     call matrix_initialize(solve_sol(1) % bvnat)
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,solve_sol(1) % bvnat,solve_sol(1) % reaction,kfl_fixno_tem)
    
  else if( variable == 'HFLUX' ) then
     ! 
     !     call tem_plugin_hflux( icoup ) 
     ! 
  else if( variable == 'FAKE ' ) then
     !! This variable is defined to allow for the use of the coupling tools

  end if

  call tem_plugin_hflux( icoup ) !< 2016Mar12 


  !=============================================================| contains |===!
contains
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine tem_plugin_hflux( icoup )
    use mod_tem_commdom, only : commdom_alya_cht_node_flux, commdom_alya_cht_nodes2bvnat 
    use def_master,      only : title, momod, modul, ID_TEMPER  
    use def_domain,      only : ndime, npoin, nboun 
    use def_coupli,      only : coupling_type, coudt 
    use def_temper,      only : bvnat_tem !, kfl_regim_tem, bvess_tem, kfl_fixno_tem
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    implicit none
    integer(ip), intent(in) :: icoup
    integer(ip)             :: idime, n_wets=-1        
    logical(ip)             :: codei, codej, vari, varj  
    integer(ip), pointer    ::    wets(:) => null()
    real(rp),    pointer    ::   varij(:) => null()
    real(rp),    pointer    ::   varji(:) => null()
    !
    codei = .False. 
    codej = .False.
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    if(inotmaster) then
       allocate( varij(npoin) ) 
       allocate( varji(npoin) )    
    else
       allocate( varij(1) )    
       allocate( varji(1) )
    endif
    !-----------------------------------------------------------------------||---!
    vari     =      'TEMPE' == coupling_type(icoup) % variable  
    codei    = current_code == coupling_type(icoup) % code_source 
    varj     =      'HFLUX' == coupling_type(icoup) % variable 
    codej    = current_code == coupling_type(icoup) % code_target 
    !-----------------------------------------------------------------------||---!
    !-----------------------------------------------------------------------||---!
    if(varj) then
       !---------------------------------------------------------------------||---!
       if( (codej.and.trim(title)/='NEUMANN').or.(codei.and.trim(title)/='DIRICHL')) then 
          print *, "ERROR: codei/='DIRICHL'.or.codej/='NEUMANN' !!"
          stop
       endif
       !
       !< 2017SEP22
       ! 
       !    if(coudt /= 1) then 
       !      print *, "ERROR: coudt/=1. set *cou.dat -> NUMERICAL_TREATMENT -> 'TIMES' !!"
       !      stop
       !    endif
       !
       !-----------------------------------------------------------| HEATF -> |---!
       if(codei.and.inotmaster) then     
          !print *, "["//trim(title)//"]."//trim( coupling_type(icoup)%variable ), "->" !, coupling_type(icoup) % relax 
          !
          varij(1:npoin) = 0.0_rp  
          call commdom_alya_cht_node_flux( varij(1:npoin) )
          varij(1:npoin) = -varij(1:npoin)  
       endif
       !---------------------------------------------------------------------||---!
       call COU_INTERPOLATE_NODAL_VALUES( icoup, 1_ip, varji(:),  varij(:) ) 
       !-----------------------------------------------------------| HEATF <- |---!
       if(codej.and.inotmaster) then
          !print *, "["//trim(title)//"]."//trim( coupling_type(icoup)%variable ), "<-"
          ! 
          call commdom_alya_cht_nodes2bvnat( varji(1:npoin), bvnat_tem(3,1:nboun,1) )
          varji(1:npoin) = 0.0_rp                                                       ! <--- reset values for the next step...
       endif
       !---------------------------------------------------------------------||---!
    endif
    !-----------------------------------------------------------------------||---!
    !-----------------------------------------------------------------------||---!
    deallocate( varij )
    deallocate( varji )
    !deallocate( aux )
    !deallocate( touched )
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
  end subroutine tem_plugin_hflux
  !-----------------------------------------------------------------------||---!


end subroutine tem_plugin


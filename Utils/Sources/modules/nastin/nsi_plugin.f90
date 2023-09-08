!------------------------------------------------------------------------
!> @addtogroup Nastin 
!> @{
!> @file    nsi_plugin.f90
!> @date    14/10/2014
!> @author  Guillaume Houzeaux
!> @brief   Plugin for coupling
!> @details Plugin for coupling
!>          1. Allocate a minimum memory so that all ranks can enter 
!>             COU_INTERPOLATE_NODAL_VALUES without blowing up
!>             (see nsi_membcs.f90 as an example)
!> @}
!------------------------------------------------------------------------
subroutine nsi_plugin(icoup)

  use def_kintyp,    only :  ip,rp
  use def_master,    only :  momod
  use def_master,    only :  modul
  use def_domain,    only :  ndime
  use def_coupli,    only :  coupling_type
  use def_coupli,    only :  UNKNOWN
  use def_coupli,    only :  RESIDUAL
  use mod_couplings, only :  COU_INTERPOLATE_NODAL_VALUES
  use def_solver,    only :  solve_sol
  use mod_matrix,    only :  matrix_initialize
  !
  ! Possible variables => 
  ! 
  use def_master,    only :  veloc
  use def_master,    only :  press,INOTMASTER
  use def_nastin,    only :  bvess_nsi
  use def_nastin,    only :  bpess_nsi
  use def_kermod,    only :  kfl_twola_ker 
  use def_nastin,    only :  btrac_nsi, tracr_nsi, tluav_nsi

  use def_master,    only : current_code,kfl_paral,inotmaster,lninv_loc
  use def_domain,    only : npoin,nboun,lnodb,lnnob,coord
  use def_kintyp,    only : lg
  use def_nastin,    only : kfl_fixno_nsi,kfl_fixpr_nsi,kfl_fixpp_nsi
  !use mod_projec
  implicit none
  integer(ip) :: ipoin,iboun,kboun,inodb
  !
  ! <= end coupling variables
  !
  integer(ip), intent(in) :: icoup    !< Coupling number
  character(5)            :: variable
  logical(lg), pointer    :: gboun(:)
  logical(lg), pointer    :: gpoin(:)

  variable = coupling_type(icoup) % variable 
  !
  ! Velocity 
  ! 
  if( variable == 'VELOC' .or. variable == 'UNKNO' ) then  
     if (kfl_twola_ker == 0_ip) then
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,bvess_nsi,veloc,kfl_fixno_nsi)
     else  ! Two-layer wall model
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,bvess_nsi,tluav_nsi,kfl_fixno_nsi)
     end if
  end if
  !
  ! Pressure
  ! 
  if( variable == 'PRESS' .or. variable == 'UNKNO' ) then   
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,bpess_nsi,press) 
     !if( current_code == 2 .and. INOTMASTER ) then
     !   !do ipoin = 1,npoin
     !      !kfl_fixpr_nsi(1,ipoin) = 0
     !      !kfl_fixpp_nsi(1,ipoin) = 0
     !    !  if(lninv_loc(ipoin)==19 ) then
     !    !     kfl_fixpr_nsi(1,ipoin) = 1
     !    !     kfl_fixpp_nsi(1,ipoin) = 1
     !    !     bpess_nsi(1,ipoin) = 1.0_rp     ! Ojo bpess_nsi cambio a (:,:,:)
     !       end if
     !   end do
     !end if
  end if
  !
  ! Momentum residual
  !
  if( variable == 'MOMEN' .or. variable == 'RESID' ) then
     call matrix_initialize(momod(modul) % solve(1) % block_array(1) % bvnat)
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,momod(modul) % solve(1) % block_array(1) % bvnat,momod(modul) % solve(1) % reaction)
  end if
  !
  ! Continuity residual
  !
  if( variable == 'CONTI' .or. variable == 'RESID' ) then
     call matrix_initialize(momod(modul) % solve(1) % block_array(2) % bvnat)
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,momod(modul) % solve(1) % block_array(2) % bvnat,momod(modul) % solve(2) % reaction)
  end if

  !
  ! Traction for two-layer wall modelling (RANS/LES coupling)
  !

  if( variable == 'TRACT' .or. variable == 'UNKNO' ) then
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,btrac_nsi,tracr_nsi,kfl_fixno_nsi)
  end if

end subroutine nsi_plugin


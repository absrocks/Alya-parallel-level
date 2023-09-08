!------------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_plugin.f90
!> @date    14/10/2014
!> @author  Guillaume Houzeaux
!> @brief   Plugin for coupling
!> @details Plugin for coupling
!> @}
!------------------------------------------------------------------------
subroutine nsa_plugin(icoup)

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
  use def_master,    only :  veloc, densi, tempe, press 
  !
  use mod_couplings,   only: COU_LIST_SOURCE_NODES
  use mod_nsa_commdom, only: nsa_commdom_phys2cons 
  use def_master,      only: TIME_N, ITER_AUX, ITER_K, inotmaster, current_code 
  use def_nastal,      only: bvess_nsa  
  use def_domain,      only: ndime, npoin
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  integer(ip), intent(in) :: icoup
  character(5)            :: variable
  real(rp)                :: dummr(2,2)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  logical(ip)          :: code_i, code_j
  integer(ip)          :: n_wets, idime
  logical(ip), pointer :: touched(:)
  integer(ip), pointer :: wets(:) => null()
  real(rp),    pointer :: aux(:,:)   
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(inotmaster) then 
    allocate( aux(ndime,npoin) ) 
    allocate( touched(npoin) )
  else 
    allocate( aux(ndime,1) )  
    allocate( touched(1) )
  endif 
  aux     = huge(1.0_rp)  
  touched = .false.
  n_wets  = -1
  !
  code_j = current_code == coupling_type(icoup) % code_target 
  if(code_j) then
    n_wets =  coupling_type(icoup) % wet % npoin_wet
    wets   => coupling_type(icoup) % wet % lpoin_wet(:)
    if(inotmaster) touched( wets(1:n_wets) ) = .true.
  endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  variable = coupling_type(icoup) % variable 
  if( variable == 'VELOC' ) then  
     call COU_INTERPOLATE_NODAL_VALUES(icoup, ndime, aux(:,:),  veloc(:,1:npoin,1) )
    do idime=1,ndime
      where(touched)
    bvess_nsa(idime,:,       1) = aux(idime,:)
        veloc(idime,:,ITER_K  ) = aux(idime,:) 
        veloc(idime,:,ITER_AUX) = aux(idime,:)  
        veloc(idime,:,TIME_N  ) = aux(idime,:)
      endwhere
    enddo 
    ! 
    call nsa_commdom_phys2cons(ITER_K)
    call nsa_commdom_phys2cons(ITER_AUX)
    call nsa_commdom_phys2cons(TIME_N)
    !
  end if
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( variable == 'DENSI' ) then
    call COU_INTERPOLATE_NODAL_VALUES(icoup,  1_ip,  aux(1,:),  densi(  1:npoin,1) )  
    where(touched) 
      bvess_nsa(ndime+1,:,       1) = aux(1,:)   
          densi(        :,ITER_K  ) = aux(1,:)
          densi(        :,ITER_AUX) = aux(1,:)
          densi(        :,TIME_N  ) = aux(1,:)
    endwhere 
    ! 
    call nsa_commdom_phys2cons(ITER_K)
    call nsa_commdom_phys2cons(ITER_AUX)
    call nsa_commdom_phys2cons(TIME_N)
    !
  end if
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( variable == 'TEMPE' ) then
    call COU_INTERPOLATE_NODAL_VALUES(icoup,  1_ip,  aux(1,:),  press(  1:npoin,1) )
    where(touched)
          press(        :,ITER_K  ) = aux(1,:)
          press(        :,ITER_AUX) = aux(1,:)
          press(        :,TIME_N  ) = aux(1,:)
          ! 
      aux(1,:) = aux(1,:)/densi(1:npoin,1)/ (8.314621_rp/0.0289531_rp) !<  T = p / rho / Raire =  101253.0 / 1.176 / (8.314621 / 0.0289531) 
      bvess_nsa(ndime+2,:,       1) = aux(1,:)
          tempe(        :,ITER_K  ) = aux(1,:)
          tempe(        :,ITER_AUX) = aux(1,:)
          tempe(        :,TIME_N  ) = aux(1,:)
    endwhere
    ! 
    call nsa_commdom_phys2cons(ITER_K)
    call nsa_commdom_phys2cons(ITER_AUX)
    call nsa_commdom_phys2cons(TIME_N)
    !
  endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( variable == 'RESID' ) then
     call matrix_initialize(solve_sol(1) % bvnat)
     call COU_INTERPOLATE_NODAL_VALUES(icoup, ndime+2_ip, solve_sol(1)%bvnat, solve_sol(1)%reaction)
  end if
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  deallocate( aux )
  deallocate( touched )
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
end subroutine nsa_plugin

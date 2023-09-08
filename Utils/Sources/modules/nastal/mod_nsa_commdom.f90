!==============================================================================!
!  I am your father...
!
!< 2015Feb09, debug: Caught signal 8 (Floating point exception) 
!
!==============================================================================!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| INIT |---!
  ! +XXXXX
  !  ! |_xxx_concou
  !  !   |_call xxx_coupli(ITASK_CONCOU)
  !  case(ITASK_CONCOU)
  !    call xxx_commdom_lmc_code_i(-1, -1)
  !    call xxx_commdom_lmc_code_j(-1, -1)
  !  !
  !  ! |_xxx_begste
  !  !   |_call nsa_coupli(ITASK_BEGSTE)
  !  case(ITASK_BEGSTE)
  !    call xxx_commdom_lmc_code_i(-1, -1)
  !    call xxx_commdom_lmc_code_j(-1, -1)
  !  !
  !  ! |_xxx_iniunk
  !  !   |_call xxx_coupli(ITASK_INIUNK)
  !  case ( ITASK_INIUNK )
  !    call xxx_commdom_lmc_code_i(-1, -1)
  !    call xxx_commdom_lmc_code_j(-1, -1)
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
module mod_nsa_commdom
  use def_parame,           only: ip, rp
  use def_domain,           only: ndime, npoin  
  use def_master,           only: IMASTER, INOTMASTER, ISEQUEN
  use def_master,           only: tempe, veloc, wmean, umome 
  use def_master,           only: veloc, densi, press
  use mod_commdom_alya,     only: ISEND, IRECV, ISENDRECV
  use mod_commdom_alya,     only: commdom_alya_calculate_driver
  use def_master,           only: ID_NASTAL
  use def_nastal,           only: adgam_nsa   ! gamma = Cp/(Cp-R)
  use def_nastal,           only: shecp_nsa   ! Specific heat Cp [J/K Kg]
  use def_nastal,           only: bvess_nsa  
  !
  use def_kintyp,           only: soltyp
  use def_master,           only: momod, modul
  !
  use mod_commdom_alya,     only: CPLNG, CC
#ifdef COMMDOM
  use mod_commdom_plepp,    only: commdom_plepp_set_vals
#endif
  implicit none
  type(soltyp), pointer :: solve(:) 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  private
  public:: nsa_commdom_lmc_code_i 
  public:: nsa_commdom_lmc_code_j 
  public:: nsa_commdom_phys2cons
  public:: nsa_commdom_cons2phys
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  !=============================================================| contains |===!
contains
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| INIT |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine nsa_commdom_lmc_code_i(sendrecv_code, when)
  use def_domain,           only: npoin, nboun
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  integer(ip), intent(in) :: sendrecv_code, when
  logical(ip) :: regim, restriction(2_ip) = .false.  
  solve => momod(modul) % solve(1:)
!print *, "solve_i", solve(1)%ndofn
#ifdef COMMDOM 
#if COMMDOM==-4
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  restriction(CPLNG(CC)%code_i) = .true.
  call commdom_alya_calculate_driver(CPLNG(CC), restriction) 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(INOTMASTER) then                                                  
    if( CPLNG(CC)%setgetvar(CPLNG(CC)%code_i,1) ) then !< ITASK_INIUNK
      !-------------------------------------------------------------------||---!
      CPLNG(CC)%var_ij(1:ndime,1:npoin) = veloc(1:ndime,1:npoin,1)
      CPLNG(CC)%var_ij(ndime+1,1:npoin) = densi(        1:npoin,1)
      CPLNG(CC)%var_ij(ndime+2,1:npoin) = press(        1:npoin,1)
      !-------------------------------------------------------------------||---!
    else&
    if( CPLNG(CC)%setgetvar(CPLNG(CC)%code_i,3) ) then !< ITASK_BEGSTE
!
!print *, "[commdom_plepp_inivar]", " 'RESIDUALi'", count( solve(1) % lpoin_reaction(1:npoin) ,KIND=ip ), solve(1) % kfl_bvnat, solve(1) % kfl_react
!
      !-------------------------------------------------------------------||---!
!      if( CPLNG(CC)%current_what(2_ip) ) then            !< WHATj == RESIDUAL  
      if( solve(1)%kfl_bvnat==1 ) then 
        solve(1)%bvnat(1:ndime+2,1:npoin) = 0.0_rp
        solve(1)%bvnat(1:ndime+2,1:npoin) = CPLNG(CC)%var_ji(1:ndime+2,1:npoin)
      endif  
      !-------------------------------------------------------------------||---!
    else&
    if( CPLNG(CC)%setgetvar(CPLNG(CC)%code_i,2) ) then !< ITASK_CONCOU
      !-------------------------------------------------------------------||---! 
      CPLNG(CC)%var_ij(1:ndime,1:npoin) = veloc(1:ndime,1:npoin,1)
      CPLNG(CC)%var_ij(ndime+1,1:npoin) = densi(        1:npoin,1)
      CPLNG(CC)%var_ij(ndime+2,1:npoin) = press(        1:npoin,1)
      !-------------------------------------------------------------------||---!
    endif
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#endif
#endif 
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine nsa_commdom_lmc_code_j(sendrecv_code, when)
  use def_domain, only: npoin, nboun
  use def_master, only: TIME_N, ITER_AUX, ITER_K
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  integer(ip), intent(in) :: sendrecv_code, when
  integer(ip) :: idime
  logical(ip) :: regim,  restriction(2_ip) = .false. 
  solve => momod(modul) % solve(1:)
#ifdef COMMDOM 
#if COMMDOM==-4
!
! print *, "REACTION", count( solve(1) % lpoin_reaction(1:npoin) ,KIND=ip ), solve(1) % kfl_react == 1, solve(1) % kfl_bvnat == 1
!      call matrix_initialize(solve_sol(1) % bvnat)
!
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  restriction(CPLNG(CC)%code_j) = .true.
  call commdom_alya_calculate_driver(CPLNG(CC), restriction)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(INOTMASTER) then
    if( CPLNG(CC)%setgetvar(CPLNG(CC)%code_j,1) ) then !< ITASK_INIUNK
      !-------------------------------------------------------------------||---!
!
!print *, "[commdom_plepp_inivar]", " 'RESIDUALj'", count( solve(1) % lpoin_reaction(1:npoin) , KIND=ip ), solve(1) % kfl_bvnat, solve(1) % kfl_react
!
      CPLNG(CC)%var_ij(1:ndime+2,1:npoin) = 0.0
!      if( CPLNG(CC)%current_what(2_ip) ) then ! what_j == RESIDUAL
      if( solve(1) % kfl_react==1 ) then 
        CPLNG(CC)%var_ij(1:ndime+2,1:npoin) = solve(1)%reaction(1:ndime+2,1:npoin)
      endif  
      !-------------------------------------------------------------------||---!
    else&
    if( CPLNG(CC)%setgetvar(CPLNG(CC)%code_j,3) ) then !< ITASK_BEGSTE
      !-------------------------------------------------------------------||---!
      do idime = 1,ndime
        call commdom_plepp_set_vals( CPLNG(CC)%var_ji(idime,1:npoin),     veloc(idime,1:npoin,ITER_K  ) )
        call commdom_plepp_set_vals( CPLNG(CC)%var_ji(idime,1:npoin),     veloc(idime,1:npoin,ITER_AUX) )
        call commdom_plepp_set_vals( CPLNG(CC)%var_ji(idime,1:npoin),     veloc(idime,1:npoin,TIME_N  ) )
        call commdom_plepp_set_vals( CPLNG(CC)%var_ji(idime,1:npoin), bvess_nsa(idime,1:npoin,       1) )
      enddo

      call   commdom_plepp_set_vals( CPLNG(CC)%var_ji(ndime+1,1:npoin),     densi(        1:npoin,ITER_K  ) )
      call   commdom_plepp_set_vals( CPLNG(CC)%var_ji(ndime+1,1:npoin),     densi(        1:npoin,ITER_AUX) )
      call   commdom_plepp_set_vals( CPLNG(CC)%var_ji(ndime+1,1:npoin),     densi(        1:npoin,TIME_N  ) )
      call   commdom_plepp_set_vals( CPLNG(CC)%var_ji(ndime+1,1:npoin), bvess_nsa(ndime+1,1:npoin,       1) )

      call   commdom_plepp_set_vals( CPLNG(CC)%var_ji(ndime+2,1:npoin),     press(        1:npoin,ITER_K  ) )
      call   commdom_plepp_set_vals( CPLNG(CC)%var_ji(ndime+2,1:npoin),     press(        1:npoin,ITER_AUX) )
      call   commdom_plepp_set_vals( CPLNG(CC)%var_ji(ndime+2,1:npoin),     press(        1:npoin,TIME_N  ) )

      !> T = p / rho / Raire =  101253.0 / 1.176 / (8.314621 / 0.0289531) = 101253.0 Pa
      where( .not.(CPLNG(CC)%var_ji(ndime+1,1:npoin)==0.0) ) CPLNG(CC)%var_ji(ndime+2,1:npoin) = CPLNG(CC)%var_ji(ndime+2,1:npoin)/CPLNG(CC)%var_ji(ndime+1,1:npoin) / (8.314621 / 0.0289531) !< 2015Feb09
      call   commdom_plepp_set_vals( CPLNG(CC)%var_ji(ndime+2,1:npoin),     tempe(        1:npoin,ITER_K  ) )
      call   commdom_plepp_set_vals( CPLNG(CC)%var_ji(ndime+2,1:npoin),     tempe(        1:npoin,ITER_AUX) )
      call   commdom_plepp_set_vals( CPLNG(CC)%var_ji(ndime+2,1:npoin),     tempe(        1:npoin,TIME_N  ) )
      call   commdom_plepp_set_vals( CPLNG(CC)%var_ji(ndime+2,1:npoin), bvess_nsa(ndime+2,1:npoin,       1) )

      call nsa_commdom_phys2cons(ITER_K)
      call nsa_commdom_phys2cons(ITER_AUX)
      call nsa_commdom_phys2cons(TIME_N)
      call nsa_setvar(0_ip,1_ip)
      !-------------------------------------------------------------------||---!
    else&
    if( CPLNG(CC)%setgetvar(CPLNG(CC)%code_j,2) ) then !< ITASK_CONCOU 
      !-------------------------------------------------------------------||---!
      CPLNG(CC)%var_ij(1:ndime+2,1:npoin) = 0.0
      if( solve(1) % kfl_react==1 ) then
!      if( CPLNG(CC)%current_what(2_ip) ) then            !< WHATj == RESIDUAL
!print *, "RESIDUALj", count( solve(1) % lpoin_reaction(1:npoin) ,KIND=ip  ), sum( solve(1) % reaction(1:ndime+2,1:npoin) ), "<--------- J"
        CPLNG(CC)%var_ij(1:ndime+2,1:npoin) = solve(1)%reaction(1:ndime+2,1:npoin)
      endif
      !-------------------------------------------------------------------||---!
    endif
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#endif
#endif
  end subroutine
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!----------------------------------------------------------------| PRIVATE |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !---------------------------------------------------------------| upcons |---!
  subroutine nsa_commdom_phys2cons( itime )
  use def_master, only: veloc, densi, press
  use def_master, only: umome, energ
  use def_domain, only: ndime, npoin 
  use def_nastal, only: gamma_nsa
  implicit none
  integer(ip), intent(in) :: itime
  integer(ip) :: ipoin
  real(rp)    :: gamme, phit(ndime+2)
  !---------------------------------------------------------------------||---!
  !
  !---------------------------------------------------------------------||---!
  !> [rho,vel,p]->[mom,rho,et]
  if(INOTMASTER) then
    do ipoin = 1,npoin
      phit(1        ) = densi(         ipoin, itime)
      phit(2:ndime+1) = veloc(1:ndime, ipoin, itime)
      phit(  ndime+2) = press(         ipoin, itime)

      gamme = gamma_nsa(ipoin)
      call phys2cons( phit(1:ndime+2), gamme, phit(1:ndime+2) )

      umome(1:ndime,ipoin, itime) = phit(1:ndime) 
      densi(        ipoin, itime) = phit(ndime+1) 
      energ(        ipoin, itime) = phit(ndime+2)
    enddo
  endif
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !---------------------------------------------------------------| upcons |---!
  subroutine nsa_commdom_cons2phys( itime )
  use def_master, only: veloc, densi, press
  use def_master, only: umome, energ
  use def_domain, only: ndime, npoin
  use def_nastal, only: gamma_nsa
  implicit none
  integer(ip), intent(in) :: itime
  integer(ip) :: ipoin
  real(rp)    :: gamme, phit(ndime+2)
  !---------------------------------------------------------------------||---!
  !
  !---------------------------------------------------------------------||---!
  !> [mom,rho,et] -> [rho,vel,p]
  if(INOTMASTER) then
    do ipoin = 1,npoin
      phit(1:ndime) = umome(1:ndime,ipoin, itime)
      phit(ndime+1) = densi(        ipoin, itime)
      phit(ndime+2) = energ(        ipoin, itime)

      gamme = gamma_nsa(ipoin)
      call cons2phys( phit(1:ndime+2), gamme, phit(1:ndime+2) )

      densi(         ipoin, itime) = phit(1        )
      veloc(1:ndime, ipoin, itime) = phit(2:ndime+1) 
      press(         ipoin, itime) = phit(  ndime+2) 
    enddo
  endif
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine cons2phys(Uk, gamme, Wk)
    implicit none
    real(rp), intent(in)  :: Uk(ndime+2), gamme
    real(rp), intent(out) :: Wk(ndime+2)
    real(rp) :: aux(ndime+2)
    !> P = (gamma-1.0) * rho * (Et - 1/2 |V|) 
    !>   = (gamme-1.0) * (ENE - 0.5*sum(MOM(1:ndime)*MOM(1:ndime))/RHO)
    !Wk(1:ndime) = Uk(1:ndime)/Uk(ndime+1) 
    !Wk(ndime+1) = Uk(ndime+1) 
    !Wk(ndime+2) = (gamme-1.0)*(Uk(ndime+2) - 0.5/Uk(ndime+1)*sum(Uk(1:ndime)*Uk(1:ndime)))  
    !> [mom, rho, et] -> [rho, vel, p] 
    aux(      1)   = Uk(ndime+1)
    aux(2:ndime+1) = Uk(1:ndime)/Uk(ndime+1)
    aux(ndime+2)   = (gamme-1.0_rp)*(Uk(ndime+2) - 0.5_rp/Uk(ndime+1)*sum(Uk(1:ndime)*Uk(1:ndime)))
    Wk = aux
  end subroutine cons2phys

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---! 
  subroutine phys2cons(Wk, gamme, Uk)
    implicit none
    real(rp), intent(in)  :: Wk(ndime+2), gamme
    real(rp), intent(out) :: Uk(ndime+2)
    real(rp) :: aux(ndime+2)
    !> rho*Et = P/(gamma-1) + 1/2*rho*|V| 
    !>        =   rho*Cv*T  + 1/2*rho*|V|
    !Uk(1:ndime) = Wk(1:ndime)*Wk(ndime+1) 
    !Uk(ndime+1) = Wk(ndime+1)
    !Uk(ndime+2) = Wk(ndime+2)/(gamme-1.0) + 0.5*Wk(ndime+1)*sum(Wk(1:ndime)*Wk(1:ndime))
    !> [rho, vel, p] -> [mom, rho, et]
    aux(1:ndime) = Wk(2:ndime+1)*Wk(1)
    aux(ndime+1) = Wk(1)
    aux(ndime+2) = Wk(ndime+2)/(gamme-1.0_rp) + 0.5_rp*Wk(1)*sum(Wk(2:ndime+1)*Wk(2:ndime+1))
    Uk = aux
  end subroutine phys2cons

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---! 
  subroutine csound(Uk, gamme, c, vmod)
    implicit none
    real(rp), intent(in)  :: Uk(ndime+2), gamme
    real(rp), intent(out) :: c, vmod
    !> c2 = gamma*P/rho = gamma(gamma-1)(Et-1/2|V|) 
    vmod = sum(Uk(1:ndime)/Uk(ndime+1)*Uk(1:ndime)/Uk(ndime+1))
    c    = Uk(ndime+2)/Uk(ndime+1) - 0.5_rp*vmod
    vmod = sqrt(vmod)
    c    = sqrt(gamme*(gamme-1.0_rp)*c)
  end subroutine csound

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine psound(Wk, gamme, c, vmod)
    implicit none
    real(rp), intent(in)  :: Wk(ndime+2), gamme
    real(rp), intent(out) :: c, vmod
    !> c^2 = gamma*p/rho = gamma*(gamma-1)*cv*T = gamma*(gamma-1)*e
    !> W = [rho, vel, p] 
    c    = Wk(ndime+2)/Wk(1)
    c    = sqrt( gamme*c )
    vmod = sqrt( dot_product(Wk(2:ndime+1), Wk(2:ndime+1)) )
  end subroutine psound


end module mod_nsa_commdom 
!==============================================================================!
!==============================================================================!

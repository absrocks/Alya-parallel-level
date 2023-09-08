!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| INIT |---!
!-------------------------------------------------------------------------||---!
module mod_nsi_commdom
  use def_parame,           only: ip, rp
  use def_domain,           only: ndime, npoin  
  use def_master,           only: IMASTER, INOTMASTER, ISEQUEN
  use def_master,           only: tempe, veloc, wmean
  use mod_commdom_alya,     only: ISEND, IRECV, ISENDRECV
  use mod_commdom_alya,     only: commdom_alya_calculate_driver 
  use mod_commdom_lm2,      only: LM2_CPLNG
  use mod_commdom_lmc,      only: LMC_CPLNG
  use def_domain,           only: walld
  use def_kermod,           only: kfl_walld
  use def_master, only: ID_NASTIN
  use def_nastin, only: kfl_regim_nsi
  use def_nastin, only: gamth_nsi   ! gamma = Cp/(Cp-R)
  use def_nastin, only: sphea_nsi   ! Specific heat Cp [J/K Kg]
  use def_nastin, only: prthe_nsi   ! Thermodynamic pressure = 101325.0_rp
  use def_nastin, only: bvess_nsi
  implicit none
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  private
    !   
    ! +Nastin  
    ! |_nsi_outvar
    !   |_ case(63_ip)
    public:: nsi_commdom_lm2_outvar 
    !
    ! +Nastin
    ! |_nsi_concou                 _
    !   |_nsi_coupli(ITASK_CONCOU)  |
    ! |_nsi_begste                 _|  
    !   |_nsi_coupli(ITASK_BEGSTE)  |
    ! |_nsi_iniunk                 _|   
    !   |_nsi_coupli(ITASK_INIUNK)  |
    !  _____________________________| 
    ! |_nsi_coupling 
    !   |_case(ITASK_CONCOU|ITASK_BEGSTE|ITASK_INIUNK)  
    public:: nsi_commdom_lm2_code_i 
    public:: nsi_commdom_lm2_code_j
    !
    public:: nsi_commdom_lmc_code_i
    public:: nsi_commdom_lmc_code_j 
    !
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  !=============================================================| contains |===!
contains
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| INIT |---!
!-------------------------------------------------------------------------||---!

  subroutine commdom_nsi_plugin()
#ifdef COMMDOM
  use mod_commdom_dynamic,  only: commdom_dynamic_check_fixno
  use mod_commdom_dynamic,  only: commdom_dynamic_set_vals
  use mod_commdom_dynamic,  only: commdom_dynamic_reduce_sum
#endif
  use mod_commdom_alya,    only: COMMDOM_COUPLING
  use mod_commdom_driver,  only: CNT_SENDRECV, CNT_SMS
  use mod_commdom_driver,  only: CNT_CPLNG, commdom_driver_exchange02
  use def_master,          only: displ, unkno
  implicit none
  integer(ip) :: idime, ipoin, idof, icomp
  real(rp) :: d_relax = 1.0_rp
  real(rp) :: n_relax = 1.0_rp
  real(rp) :: residual2(3,2) = 0.0_rp
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( any(CNT_SENDRECV) ) then
#ifdef COMMDOM
#if   COMMDOM==2 
    !---------------------------------------------------- code_i==1==FLUID |---!
    code_i: &
    if(CNT_CPLNG%current_code==CNT_CPLNG%code_i) then
      if( CNT_SENDRECV(7) ) then
        !-----------------------------------------------------------| U--> |---!
        to_send01: &
        if(inotmaster) then
          CNT_CPLNG%var_ij(1:ndime,1:npoin) = 0.0
          CNT_CPLNG%var_ij(1:ndime,1:npoin) = veloc(1:ndime,1:npoin,1)
        endif to_send01
        !-----------------------------------------------------------------||---!
        !
        call commdom_driver_exchange02( CNT_CPLNG, debug=.true.)
        !
        !--------------------------------------------------------| dUdn<-- |---!
        to_recv01: &
        if(inotmaster) then
!
!          call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(3,1:npoin), solve(1)%bvnat(3,1:npoin), relax_op=n_relax, res2=residual2(3,1:2) )
!
        endif to_recv01
        !-----------------------------------------------------------------||---!
      endif
    endif code_i
    !---------------------------------------------------------------------||---!
    !---------------------------------------------------| code_j==2==SOLID |---!
    code_j: &
    if(CNT_CPLNG%current_code==CNT_CPLNG%code_j) then
      if( CNT_SENDRECV(7) ) then
        !-----------------------------------------------------------| U--> |---!
        to_send02: &
        if(inotmaster) then
!
!          CNT_CPLNG%var_ij(1:ndime,1:npoin) = 0.0
!
        endif to_send02
        !-----------------------------------------------------------------||---!
        !
        call commdom_driver_exchange02( CNT_CPLNG, debug=.true.)
        !
        !--------------------------------------------------------| dUdn<-- |---!
        to_recv02: &
        if(inotmaster) then
          do idime = 1,ndime
            call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(idime,1:npoin), bvess_nsi(idime,1:npoin,1), relax_op=n_relax, res2=residual2(3,1:2), debug=.false. )
            call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(idime,1:npoin),     veloc(idime,1:npoin,1), relax_op=n_relax, res2=residual2(3,1:2), debug=.false. )
            CNT_CPLNG%var_ji(idime,1:npoin) = 0.0
          enddo
        endif to_recv02
        !-----------------------------------------------------------------||---!
      endif
    endif code_j
    !---------------------------------------------------------------------||---!
#endif
#endif 
    ! 
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine nsi_commdom_lm2_code_i(sendrecv_code, when)
  use def_domain,           only: npoin, nboun
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  integer(ip), intent(in) :: sendrecv_code, when
  logical(ip) :: set_var_ij, get_var_ji, ini_var_ij 
  logical(ip) :: regim
  integer(ip) :: idx_exch = 2_ip  
#ifdef COMMDOM 
#if COMMDOM==2
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  call commdom_alya_calculate_driver(LM2_CPLNG) 
  regim = .true. !kfl_regim_nsi==0 
  ini_var_ij = LM2_CPLNG%setgetvar(LM2_CPLNG%code_i,1).and.regim
  set_var_ij = LM2_CPLNG%setgetvar(LM2_CPLNG%code_i,2).and.regim
  get_var_ji = LM2_CPLNG%setgetvar(LM2_CPLNG%code_i,3).and.regim
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(INOTMASTER) then                                                  
    if(ini_var_ij) then !< ITASK_INIUNK
      print *, "[nsi_commdom_lm2_i]", " ini_var_ij"
      call nsi_commdom_lm2_set_vel( veloc(1,1:npoin,1), LM2_CPLNG%var_ij(1,1:LM2_CPLNG%n_pts) )
      call nsi_commdom_lm2_set_vel( veloc(2,1:npoin,1), LM2_CPLNG%var_ij(2,1:LM2_CPLNG%n_pts) )
      call nsi_commdom_lm2_set_vel( veloc(3,1:npoin,1), LM2_CPLNG%var_ij(3,1:LM2_CPLNG%n_pts) )
    else&
    if(get_var_ji) then !< ITASK_BEGSTE
      print *, "[nsi_commdom_lm2_i]", " set_var_ij"
    else&
    if(set_var_ij) then !< ITASK_CONCOU 
      print *, "[nsi_commdom_lm2_i]", " get_var_ji"
      call nsi_commdom_lm2_set_vel( veloc(1,1:npoin,1), LM2_CPLNG%var_ij(1,1:LM2_CPLNG%n_pts) )
      call nsi_commdom_lm2_set_vel( veloc(2,1:npoin,1), LM2_CPLNG%var_ij(2,1:LM2_CPLNG%n_pts) )
      call nsi_commdom_lm2_set_vel( veloc(3,1:npoin,1), LM2_CPLNG%var_ij(3,1:LM2_CPLNG%n_pts) )
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
  subroutine nsi_commdom_lm2_code_j(sendrecv_code, when)
  use def_domain,           only: npoin, nboun
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  integer(ip), intent(in) :: sendrecv_code, when
  logical(ip) :: set_var_ij, get_var_ji, ini_var_ij
  logical(ip) :: regim
  integer(ip) :: idx_exch = 2_ip
#ifdef COMMDOM 
#if COMMDOM==2
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  call commdom_alya_calculate_driver(LM2_CPLNG)
  regim = .true. !kfl_regim_nsi==0 
  ini_var_ij = LM2_CPLNG%setgetvar(LM2_CPLNG%code_j,1).and.regim
  set_var_ij = LM2_CPLNG%setgetvar(LM2_CPLNG%code_j,2).and.regim
  get_var_ji = LM2_CPLNG%setgetvar(LM2_CPLNG%code_j,3).and.regim
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(INOTMASTER) then
    if(ini_var_ij) then !< ITASK_INIUNK
      print *, "[nsi_commdom_lm2_j]", " ini_var_ij"
    else&
    if(get_var_ji) then !< ITASK_BEGSTE
      print *, "[nsi_commdom_lm2_j]", " set_var_ij"
      call nsi_commdom_lm2_get_vel(LM2_CPLNG%var_ji(1,1:LM2_CPLNG%n_pts),     veloc(1,1:npoin,1))
      call nsi_commdom_lm2_get_vel(LM2_CPLNG%var_ji(2,1:LM2_CPLNG%n_pts),     veloc(2,1:npoin,1))
      call nsi_commdom_lm2_get_vel(LM2_CPLNG%var_ji(3,1:LM2_CPLNG%n_pts),     veloc(3,1:npoin,1))
      call nsi_commdom_lm2_get_vel(LM2_CPLNG%var_ji(1,1:LM2_CPLNG%n_pts), bvess_nsi(1,1:npoin,1))
      call nsi_commdom_lm2_get_vel(LM2_CPLNG%var_ji(2,1:LM2_CPLNG%n_pts), bvess_nsi(2,1:npoin,1))
      call nsi_commdom_lm2_get_vel(LM2_CPLNG%var_ji(3,1:LM2_CPLNG%n_pts), bvess_nsi(3,1:npoin,1))
    else&
    if(set_var_ij) then !< ITASK_CONCOU 
      print *, "[nsi_commdom_lm2_j]", " get_var_ji"
!      LM2_CPLNG%var_ij(idx_exch,1_ip:LM2_CPLNG%n_pts) = -4.0
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
!--------------------------------------------------------------------| LMC |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine nsi_commdom_lmc_code_i(sendrecv_code, when)
  use def_domain,           only: npoin, nboun
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  integer(ip), intent(in) :: sendrecv_code, when
  logical(ip) :: set_var_ij, get_var_ji, ini_var_ij
  logical(ip) :: regim
  integer(ip) :: idx_exch = 2_ip
#ifdef COMMDOM 
#if COMMDOM==-3
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  call commdom_alya_calculate_driver(LMC_CPLNG)
  regim = .true. !kfl_regim_nsi==0 
  ini_var_ij = LM2_CPLNG%setgetvar(LMC_CPLNG%code_i,1).and.regim
  set_var_ij = LM2_CPLNG%setgetvar(LMC_CPLNG%code_i,2).and.regim
  get_var_ji = LM2_CPLNG%setgetvar(LMC_CPLNG%code_i,3).and.regim
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(INOTMASTER) then
    if(ini_var_ij) then !< ITASK_INIUNK
      print *, "[nsi_commdom_lmc_i]", " ini_var_ij"
      call nsi_commdom_lm2_set_vel( veloc(1,1:npoin,1), LMC_CPLNG%var_ij(1,1:LMC_CPLNG%n_pts) )
      call nsi_commdom_lm2_set_vel( veloc(2,1:npoin,1), LMC_CPLNG%var_ij(2,1:LMC_CPLNG%n_pts) )
      call nsi_commdom_lm2_set_vel( veloc(3,1:npoin,1), LMC_CPLNG%var_ij(3,1:LMC_CPLNG%n_pts) )
    else&
    if(get_var_ji) then !< ITASK_BEGSTE
      print *, "[nsi_commdom_lmc_i]", " set_var_ij"
    else&
    if(set_var_ij) then !< ITASK_CONCOU 
      print *, "[nsi_commdom_lmc_i]", " get_var_ji"
      call nsi_commdom_lm2_set_vel( veloc(1,1:npoin,1), LMC_CPLNG%var_ij(1,1:LMC_CPLNG%n_pts) )
      call nsi_commdom_lm2_set_vel( veloc(2,1:npoin,1), LMC_CPLNG%var_ij(2,1:LMC_CPLNG%n_pts) )
      call nsi_commdom_lm2_set_vel( veloc(3,1:npoin,1), LMC_CPLNG%var_ij(3,1:LMC_CPLNG%n_pts) )
    endif
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#endif
#endif
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine nsi_commdom_lmc_code_j(sendrecv_code, when)
  use def_domain,           only: npoin, nboun
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  integer(ip), intent(in) :: sendrecv_code, when
  logical(ip) :: set_var_ij, get_var_ji, ini_var_ij
  logical(ip) :: regim
  integer(ip) :: idx_exch = 2_ip
#ifdef COMMDOM 
#if COMMDOM==-3
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  call commdom_alya_calculate_driver(LMC_CPLNG)
  regim = .true. !kfl_regim_nsi==0 
  ini_var_ij = LM2_CPLNG%setgetvar(LMC_CPLNG%code_j,1).and.regim
  set_var_ij = LM2_CPLNG%setgetvar(LMC_CPLNG%code_j,2).and.regim
  get_var_ji = LM2_CPLNG%setgetvar(LMC_CPLNG%code_j,3).and.regim
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(INOTMASTER) then
    if(ini_var_ij) then !< ITASK_INIUNK
      print *, "[nsi_commdom_lmc_j]", " ini_var_ij"
    else&
    if(get_var_ji) then !< ITASK_BEGSTE
      print *, "[nsi_commdom_lmc_j]", " set_var_ij"
      call nsi_commdom_lm2_get_vel(LMC_CPLNG%var_ji(1,1:LMC_CPLNG%n_pts),     veloc(1,1:npoin,1))
      call nsi_commdom_lm2_get_vel(LMC_CPLNG%var_ji(2,1:LMC_CPLNG%n_pts),     veloc(2,1:npoin,1))
      call nsi_commdom_lm2_get_vel(LMC_CPLNG%var_ji(3,1:LMC_CPLNG%n_pts),     veloc(3,1:npoin,1))
      call nsi_commdom_lm2_get_vel(LMC_CPLNG%var_ji(1,1:LMC_CPLNG%n_pts), bvess_nsi(1,1:npoin,1))
      call nsi_commdom_lm2_get_vel(LMC_CPLNG%var_ji(2,1:LMC_CPLNG%n_pts), bvess_nsi(2,1:npoin,1))
      call nsi_commdom_lm2_get_vel(LMC_CPLNG%var_ji(3,1:LMC_CPLNG%n_pts), bvess_nsi(3,1:npoin,1))
    else&
    if(set_var_ij) then !< ITASK_CONCOU 
      print *, "[nsi_commdom_lmc_j]", " get_var_ji"
!      LM2_CPLNG%var_ij(idx_exch,1_ip:LM2_CPLNG%n_pts) = -4.0
    endif
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#endif
#endif 
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| INIT |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine nsi_commdom_lm2_set_vel( prop, var_ij )
  implicit none
  real(rp), intent(in ) ::   prop(npoin)
  real(rp), intent(out) :: var_ij(npoin)
  integer(ip) :: ipoin
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( INOTMASTER ) then
    var_ij(1:npoin) = prop(1:npoin)
  end if
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine nsi_commdom_lm2_get_vel( var_ji, prop )
!use mod_commdom_plepp, only: PLEPP_CPLNG

  implicit none
  real(rp), intent(in ) :: var_ji(npoin)
  real(rp), intent(out) ::   prop(npoin)
  integer(ip) :: ipoin
  integer(ip), pointer :: in_list(:) => null()
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( INOTMASTER ) then
    prop(1:npoin) = var_ji(1:npoin)
!
!in_list => PLEPP_CPLNG%interior_list_j
!print *, in_list
!print *, bvess_nsi(1:ndime,in_list,1) 
!print *, "", sum( bvess_nsi(1:ndime,in_list,1), dim=2 )/PLEPP_CPLNG%n_ji
!print *, "", sum( var_ji(1:ndime,in_list,1) )/PLEPP_CPLNG%n_ji
!
  end if
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| LODI |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine nsi_commdom_lm2_outvar( param ) 
  use def_parame, only: zero
  use def_master, only: gevec
  implicit none
  real(rp), intent(in) :: param(npoin,ndime) 
  integer(ip) :: ipoin 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( INOTMASTER ) then
    call memgen(zero,ndime,npoin)
    gevec = 0.0_rp
    do ipoin = 1,npoin
      gevec(1:ndime,ipoin) = param(ipoin,1:ndime) 
    enddo
  end if
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine nsi_commdom_lm2_lodi( param )
  use def_parame, only: zero
  use def_master, only: gevec, tempe, veloc, wmean, kfl_coupl 
  use def_master, only: ID_CHEMIC
  use def_kermod, only: gasco       ! Gas constant R [J/K Kg] 
  implicit none
  real(rp), intent(in) :: param(npoin,ndime)
  integer(ip) :: ipoin
  real(rp)    :: sofac, xsoun, vmodu
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( INOTMASTER ) then
        sofac=sqrt(gamth_nsi*gasco)
        do ipoin=1,npoin
           xsoun = sofac*sqrt(abs(tempe(ipoin,1)))
           if(kfl_coupl(ID_NASTIN,ID_CHEMIC) > 0 ) xsoun = xsoun / sqrt(abs(wmean(ipoin,1)))
           vmodu = sqrt( dot_product( veloc(1:ndime,ipoin,1), veloc(1:ndime,ipoin,1) ) )
        end do
    !allocate( gradv_nsi(ntens,npoin), stat=istat)
    !call memchk(zero,istat,mem_modul(1:2,modul),'GRADV_NSI','nsi_wyplus',gradv_nsi)
    !call graten( unkno(1:npoin),         grunk(1:ndime,1:npoin) )
    !call graten( unkno(1:ndime,1:npoin), grunk(1:ntens,1:npoin) )
    !call gravec( unkno(1:ndime,1:npoin), grunk(1:ndime,1:ndime,1:npoin) )
  end if
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !--------------------------------------------------------------| PRIVATE |---!
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine SMTRX(Ak, Wk, kk, cvel)
    implicit none
    real(rp), intent(inout) :: Ak(ndime+2, ndime+2)
    real(rp), intent(in)    :: Wk(ndime+2), cvel
    integer(ip), intent(in) :: kk
    real(rp)    :: inv_c, inv_c2, inv_c_rho
    integer(ip) :: ii
    !> W  = [rho, u1, u2, u3, p]^T 
    inv_c     = 1.0_rp/cvel
    inv_c2    = inv_c*inv_c
    inv_c_rho = 1.0_rp/Wk(1)*inv_c
    Ak(1:ndime+2,1:ndime+2)   = 0.0_rp
    !> |
    Ak(1,(/1_ip,ndime+2/)) = inv_c2*0.5_rp
    Ak(1,         kk+1) = inv_c2
    Ak(ndime+2, (/1_ip,ndime+2/)) = 0.5_rp
    !> _
    Ak(kk+1,      1) = -0.5_rp*inv_c_rho
    Ak(kk+1,ndime+2) =  0.5_rp*inv_c_rho
    !> |_ 
    do ii = 2,ndime+1
      Ak(ii,ii) = 1
    enddo
    Ak(kk+1,kk+1) = 0.0_rp
  end subroutine
 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine SINVMTRX(Ak, Wk, kk, cvel)
    implicit none
    real(rp), intent(out)   :: Ak(ndime+2, ndime+2)
    real(rp), intent(in)    :: Wk(ndime+2), cvel
    integer(ip), intent(in) :: kk
    real(rp)    :: c2, c_rho, rho, vk
    real(rp)    :: lambda(ndime+2)
    integer(ip) :: ii
    !> W  = [rho, u1, u2, u3, p]^T 
    !> S^{-1}_k * F^{k} S_k = \Lambda_k
    !> F^k = [rho*u_i*u_k+p, rho*u_k, (rho*e+p)*u_k]^T  
    rho   = Wk(1)
    vk    = Wk(kk+1)
    c2    = cvel*cvel
    c_rho = cvel*rho
    Ak(:,:) = 0.0_rp
    !> |
    Ak(1,ndime+2) =  1.0_rp
    Ak(1,   kk+1) = -c_rho
    Ak(ndime+2,   kk+1) = c_rho
    Ak(ndime+2,ndime+2) = 1.0_rp
    !> _
    Ak(kk+1,      1) =  c2
    Ak(kk+1,ndime+2) = -1.0_rp
    !> |_ 
    do ii = 2,ndime+2
      Ak(ii,ii) = 1
    enddo
    Ak(kk+1,kk+1) = 0.0_rp

    lambda(2:ndime+1) = vk
    lambda(        1) = vk - cvel
    lambda(ndime+2) = vk + cvel

    do ii = 1,ndime+2
      Ak(1:ndime+2,ii) = Ak(1:ndime+2,ii)*lambda(1:ndime+2)
    enddo
  end subroutine SINVMTRX

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

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---! 
  !=============================================================| contains |===!
end module mod_nsi_commdom 
!==============================================================================!
!==============================================================================!
!
!< 'nsi_inivar.f90'
!   postp(1) % wopos ( 1,63) = 'CHARZ' !< commdom 
!   postp(1) % wopos ( 1,63) = 'VECTO'
!
!==============================================================================!
!==============================================================================!

!==============================================================================!
  !
  !< 2015Mar20 -> created (from 'mod_commdom_plepp') 
  !
!==============================================================================!
module mod_commdom_aitken 
  use mod_commdom_alya,  only: INONE 
#ifdef COMMDOM 
  use def_parame,        only: ip, rp
  use def_master,        only: inotmaster 
  use def_domain,        only: coord, mnode, nelem, ndime, npoin
  use def_domain,        only: LESET 
  use def_domain,        only: ltype, nnode, ngaus, lnods, coord
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none 
  integer(ip),   parameter :: cp   = 64  
  character(cp), parameter :: frmt = '(E11.4)' 
  !
  logical(ip) :: initiated = .false.  
  logical(ip) :: debug     = .true.  
  !-----------------------------------------------------------------------||---!
  type COMMDOM_DYNAMIC_RELAXATION
    real(rp),    pointer ::      PHI(:,:,:,:) => null()   !< PHI_{m}^{n} -> PHI(ninter,ndof,niter,ntime); niter:iteration loop; ntime:time loop 
    real(rp),    pointer ::      RES(:,:,:)   => null()   !<      RES(ninter, ndof, niter-1)   
    real(rp),    pointer ::    omega(:,:)     => null()   !<    omega(        ndof, niter-1)  
    real(rp),    pointer :: mod_res2(:,:)     => null()   !< mod_res2(        ndof, niter-1)  
    real(rp)             :: omega_max    =  1.0_rp 
    real(rp)             :: error        =  1e-3_rp 
    real(rp)             :: convergence  =  0_ip
    integer(ip)          :: ninter       =  0_ip 
    integer(ip)          :: ndof         =  0_ip 
    integer(ip)          :: ntime        =  0_ip 
    integer(ip)          :: niter        =  0_ip
    integer(ip)          :: now(2_ip)    =  0_ip
    logical(rp)          :: converged    = .False.
    logical(rp)          :: ready        = .False.
  end type COMMDOM_DYNAMIC_RELAXATION
  !
  type(COMMDOM_DYNAMIC_RELAXATION) RELAXATION
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  private
    public :: RELAXATION
    public :: commdom_aitken_init
    public :: commdom_aitken_deallocate
    public :: commdom_aitken_allocate
    public :: commdom_aitken_set_prop
    public :: commdom_aitken_set_vals
    !public :: 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  contains
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| INIT |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_aitken_init( RELAX, ninter, ndof, niter, omega_max) 
  implicit none
  type(COMMDOM_DYNAMIC_RELAXATION), intent(inout) :: RELAX 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  integer(ip), intent(in) :: ninter 
  integer(ip), intent(in) :: ndof
  integer(ip), intent(in) :: niter
  real(rp),    intent(in) :: omega_max
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  call commdom_aitken_deallocate( RELAX ) 
  !
  if(INOTMASTER) then 
    RELAX%ninter     = ninter
    RELAX%ndof       = ndof
    RELAX%ntime      = 1
    RELAX%niter      = niter
    RELAX%omega_max  = omega_max 
  endif 
  !
  if(debug) print*, "[commdom_aitken_init]"
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_aitken_deallocate( RELAX ) 
  implicit none
  type(COMMDOM_DYNAMIC_RELAXATION), intent(inout) :: RELAX 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( initiated ) then 
    if( associated(RELAX%PHI     ) ) deallocate(  RELAX%PHI      )  
    if( associated(RELAX%RES     ) ) deallocate(  RELAX%RES      ) 
    if( associated(RELAX%mod_res2) ) deallocate(  RELAX%mod_res2 ) 
    if( associated(RELAX%omega   ) ) deallocate(  RELAX%omega    ) 
    !
    initiated = .false.
    !
    if(debug) print*, "[commdom_aitken_deallocate]"
  endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_aitken_allocate( RELAX ) 
  implicit none
  type(COMMDOM_DYNAMIC_RELAXATION), intent(inout) :: RELAX 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( .not.initiated ) then 
    !
    if( .not.associated( RELAX%PHI      ) ) allocate( RELAX%PHI     ( RELAX%ndof,RELAX%ninter,RELAX%niter,RELAX%ntime) )
    if( .not.associated( RELAX%RES      ) ) allocate( RELAX%RES     ( RELAX%ndof,RELAX%ninter,RELAX%niter            ) )
    if( .not.associated( RELAX%mod_res2 ) ) allocate( RELAX%mod_res2( RELAX%ndof,             RELAX%niter            ) )
    if( .not.associated( RELAX%omega    ) ) allocate( RELAX%omega   ( RELAX%ndof,             RELAX%niter            ) )
    !
    RELAX%PHI      = 0.0 
    RELAX%RES      = 0.0 
    RELAX%mod_res2 = 0.0
    RELAX%omega    = 0.0
    !
    initiated = .true.
    !
    if(debug) print*, "[commdom_aitken_allocate]"
  endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  


!-------------------------------------------------------------------------||---!
!---------------------------------------------------------------| EXCHANGE |---!
!-------------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_aitken_set_prop(RELAX, PHI, iiter, itime_opt)
  implicit none
  type(COMMDOM_DYNAMIC_RELAXATION), intent(inout) :: RELAX 
  real(rp),                         intent(in   ) :: PHI(RELAX%ndof,RELAX%ninter) 
  integer(ip),                      intent(in   ) :: iiter 
  integer(ip), optional,            intent(in   ) :: itime_opt
  integer(ip) :: itime  = 1
  integer(ip) :: ninter = 0, ndof = 0, k = 0, idof = 0 
  logical(ip) :: ready  =.false. 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  if( present(itime_opt) ) itime = itime_opt
  !
  if( (itime<=0).or.(itime>RELAX%ntime) ) then 
    print*, "[commdom_aitken_set_prop]", " ERROR: itime>RELAX%ntime!!", itime, RELAX%ntime
    call runend("EXIT!!")
  endif 
  !
  if( .not.all( shape(PHI)==(/RELAX%ndof,RELAX%ninter/) ) ) then 
    print*, "[commdom_aitken_set_prop]", " ERROR: array shapes!!", shape(PHI), "/=", (/RELAX%ndof, RELAX%ninter/) 
    call runend("EXIT!!")
  endif 
  !
  if(.not.initiated) then 
    print*, "[commdom_aitken_set_prop]", " ERROR: initiated='false'" 
    call runend("EXIT!!")
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  !
  ! omega(v) = -omega(v-1) * R(v-1) * [R(v) - R(v-1)] / [R(v) - R(v-1)]**2
  RELAX%now = (/ itime, iiter /)
  !
  if(INOTMASTER) then
    !
    ninter = RELAX%ninter
    ndof   = RELAX%ndof
    k      = RELAX%niter !iiter
    !
    if(iiter-2>0) then !< PHI_{k-2}^{n} = PHI_{k-1}^{n} 
      RELAX%PHI(1:ndof,1:ninter,k-2,itime) = RELAX%PHI(1:ndof,1:ninter,k-1,itime) 
    endif 
    !
    if(iiter-1>0) then !< PHI_{k-1}^{n} = PHI_{k-0}^{n} 
      RELAX%PHI(1:ndof,1:ninter,k-1,itime) = RELAX%PHI(1:ndof,1:ninter,k-0,itime) 
      !
      ! R_{k-1}^{n} = PHI_{k-1}^{n} - PHI_{k-2}^{n}
      RELAX%RES(1:ndof,1:ninter,k-1      ) = RELAX%PHI(1:ndof,1:ninter,k-1,itime) - &
                                             RELAX%PHI(1:ndof,1:ninter,k-2,itime)
      !
      do idof = 1,ndof
        RELAX%mod_res2(idof,k-1) = dot_product( RELAX%RES(idof,1:ninter,k-1), RELAX%RES(idof,1:ninter,k-1) ) 
      enddo 
      !print *, "mod_res2(k-1)", RELAX%mod_res2(idof,k-1)
      !
    endif
    ! 
    if(iiter-0>0) then !< PHI_{k-0}^{n} = U^{n} 
      RELAX%PHI(1:ndof,1:ninter,k-0,itime) = PHI(1:ndof,1:ninter) 
      !
      ! R_{k-0}^{n} = PHI_{k-0}^{n} - PHI_{k-1}^{n}
      RELAX%RES(1:ndof,1:ninter,k-0      ) = RELAX%PHI(1:ndof,1:ninter,k-0,itime) - &
                                             RELAX%PHI(1:ndof,1:ninter,k-1,itime)
      !
      do idof = 1,ndof
        RELAX%mod_res2(idof,k-0) = dot_product( RELAX%RES(idof,1:ninter,k-0), RELAX%RES(idof,1:ninter,k-0) ) 
      enddo 
      print *, "mod_res2(k-0)", RELAX%mod_res2(1:ndof,k-0)
      !
    endif
    !
  endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---! 
  ready = all(iiter + (/0,-1,-2/) > 0)   
  if( ready ) then 
    RELAX%ready = .true.
    if(debug) print*, "[commdom_aitken_set_prop]" 
  endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_aitken_get_relax( RELAX )
  use mod_communications, only :  PAR_SUM
  implicit none
  type(COMMDOM_DYNAMIC_RELAXATION), intent(inout) :: RELAX 
  integer(ip), parameter :: v = 2
  integer(ip)            :: ninter = 0_ip
  integer(ip)            :: iiter  = 1_ip 
  integer(ip)            :: idof   = 1_ip 
  integer(ip)            :: itime  = 1_ip
  ! 
  real(rp) :: omega_num = 0.0
  real(rp) :: omega_den = 0.0
  real(rp) :: mod_res2  = 0.0
  !
  ninter = RELAX%ninter
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  !
  ! omega(v+1) = -omega(v) * R(v) * [R(v+1) - R(v)] / [R(v+1) - R(v)]**2
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( RELAX%ready ) then 
    if(INOTMASTER) then  
      ! R(v+1) = PHI(v+1) - PHI(v  ) 
      RELAX%RES(:,idof,v+1) = RELAX%PHI(:,idof,v+1,itime) - RELAX%PHI(:,idof,v+0,itime)
      ! R(v  ) = PHI(v  ) - PHI(v-1)
      RELAX%RES(:,idof,v+0) = RELAX%PHI(:,idof,v+0,itime) - RELAX%PHI(:,idof,v-1,itime)
      !
      mod_res2  = dot_product( RELAX%RES(:,idof,v+1), RELAX%RES(:,idof,v+1) ) 
      !
      ! omega(v+1) = -omega(v) * R(v) * [R(v+1) - R(v)] / [R(v+1) - R(v)]**2
      omega_num = dot_product( RELAX%RES(:,idof,v+0)                        , RELAX%RES(:,idof,v+1) - RELAX%RES(:,idof,v+0) )
      omega_den = dot_product( RELAX%RES(:,idof,v+1) - RELAX%RES(:,idof,v+0), RELAX%RES(:,idof,v+1) - RELAX%RES(:,idof,v+0) ) 
      !
    endif 
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    call PAR_SUM(mod_res2,  'IN MY CODE' )
    call PAR_SUM(omega_num, 'IN MY CODE' )
    call PAR_SUM(omega_den, 'IN MY CODE' )
!    RELAX%mod_res2 = mod_res2
!    RELAX%omega    = omega_num/omega_den
    !
  endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  



  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_aitken_set_vals(PHI_1, PHI_2, relax, modR2, debug)
  !use mod_commdom_plepp, only: PLEPP_CPLNG
  implicit none
  real(rp),              intent(in ) :: PHI_1(npoin)
  real(rp),              intent(out) :: PHI_2(npoin)
  real(rp),    optional, intent(in ) :: relax
  real(rp),    optional, intent(out) :: modR2
  logical(ip), optional, intent(in ) :: debug
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
!  call commdom_dynamic_set_vals(PLEPP_CPLNG, PHI_1, PHI_2, relax, modR2, debug) 
PHI_2 = 0.0
modR2 = 0.0 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_aitken_relaxation( RELAX ) 
  implicit none
  type(COMMDOM_DYNAMIC_RELAXATION), intent(inout) :: RELAX 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  ! 
  ! i: Stage
  ! v: fixed point iteration 
  ! 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  integer(ip), parameter :: v = 2
  real(rp), pointer      :: PHI(:,:) !   PHI  = <PHI(v+2), PHI(v+1), PHI(v+0)>
  real(rp), pointer      ::   R(:,:) ! R(v+1) = PHI(v+1) - PHI(v+0)
  !
  real(rp)               :: modR  
  real(rp)               :: tolerance 
  logical(ip)            :: convergence 
  integer(ip)            :: n_inter
  !
  real(rp)               :: omega(3)
  real(rp)               :: omega_0   = -1, omega_n = -1
  real(rp)               :: omega_max = 0.8 !< philipp, tobias, 2014 
  !
  R(:,v+1) = PHI(:,v+1) - PHI(:,v+0) ! R(v+1) = PHI(v+1) - PHI(v  ) 
  R(:,v+0) = PHI(:,v+0) - PHI(:,v-1) ! R(v  ) = PHI(v  ) - PHI(v-1)
  !
  modR        = dot_product( R(:,v+1), R(:,v+1) ) 
  convergence = sqrt( modR/n_inter ) < tolerance 
  ! 

  ! 
  omega(v+1) = & 
               dot_product( R(:,v+0)           , R(:,v+1) - R(:,v+0) ) / & 
               dot_product( R(:,v+1) - R(:,v+0), R(:,v+1) - R(:,v+0) ) 
  !
  omega(v+1) = -omega(v+1) * omega(v)
  !
  !< ulrich wolfgang 2007 
  omega_0 = max( omega_n, omega_max)
  !
  !< Dregroote, Souto, 2010 
  omega_0 =  min( abs(omega_n), omega_max)
  omega_0 = sign(      omega_0, omega_n  )
  !  
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| ++++ |---!
!-------------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| ++++ |---!
!-------------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| ++++ |---!
!-------------------------------------------------------------------------||---!
#endif 
end module mod_commdom_aitken 
!==============================================================================!
  !-----------------------------------------------------| OUTER_ITERATIONS |---!
  !
  !   |_Alya                                       
  !     |_call Turnon()                            
  !     |_call Iniunk()                             
  !     |_time: do while
  !       |_call Timste()                          
  !       |_reset: do 
  !         |_call Begste()                                    coupling_driver_iteration(1:max_block_cou)  = 0
  !           |_block: do while                          
  !             |_coupling: do while                    
  !               |_call Begzon()                              coupling_driver_iteration( iblok ) += 1
  !               |_modules: do while                               
  !                 |_call Doiter()                
  !                 |_call Concou()                
  !               |_call Endzon()                              call COU_CHECK_CONVERGENCE
  !                                                            call cou_cvgunk                
  !             |_call Conblk()                                coupling_driver_iteration( iblok )  = 0 
  !       |_call Endste()                                     
  !
  !-----------------------------------------------------------------------||---!
!==============================================================================!

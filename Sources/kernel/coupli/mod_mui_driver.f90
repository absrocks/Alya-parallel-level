!==============================================================================!
!  I am your father...
! 
!==============================================================================!
module mod_mui_driver
  use def_parame,           only: ip, rp
  use def_master,           only: inotmaster, imaster, isequen, title ,ITASK_TIMSTE
  use def_domain,           only: coord, mnode, ndime, npoin 
  use def_domain,           only: ltype, lnods
  use def_kintyp,           only: soltyp
  use def_master,           only: momod, modul
  use mod_mui,              only: MUI_COUPLING
  use mod_messages,         only : livinf
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  implicit none 
  !
  logical(ip)  :: debug = .false.
  !
  character(6) :: name_task(20) = (/ 'REAPRO', 'TURNON', 'INIUNK', &
                                     'TIMSTE', 'BEGSTE', 'DOITER', &
                                     'CONCOU', 'CONBLK', 'NEWMSH', &
                                     'ENDSTE', 'FILTER', 'OUTPUT', &
                                     'TURNOF', 'BEGITE', 'ENDITE', &
                                     'MATRIX', 'DOOPTI', 'ENDOPT', &
                                     'BEGZON', 'ENDZON' /) 

  character(6) :: name_when(2) = (/ 'BEFORE', 'AFTERE'/) 
  !
  integer(ip), parameter   :: n_times = 9                                      !< 2016MAR26. 7 -> 8. 2016JUN06 8 -> 9 
  logical(ip)              :: MUI_SENDRECV(n_times) = .false. 
  character(64)            :: MUI_SMS               = ' '
  !
  real(rp)                 :: n_max_its
  logical(ip)              :: dynamic
  !
  logical(ip) :: residual   = .false. 
  !
  logical(ip) :: exchanged  = .false.
  !
  type(soltyp), pointer    :: solve(:)
  !
  integer(ip) :: modules_i(3) = -1  !< 2016JUN06, 2016MAR29
  integer(ip) :: modules_j(3) = -1
  !
  private
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  !
  public :: MUI_SENDRECV
  public :: MUI_SMS
  public :: mui_driver_init
  public :: mui_driver_sendrecv
  public :: mui_driver_exchange02 
  public :: mui_driver_set_mesh        !< 2016MAR30 
  !
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  contains
!-------------------------------------------------------------------------||---!
!-----------------------------------------------------------------| PUBLIC |---!
!-------------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  subroutine mui_driver_init_fsi(CPLNG)                                    !< 2016MAR29 
  use def_coupli,       only: UNKNOWN, RESIDUAL
  use def_coupli,       only: mcoup
  use def_master,       only: ID_NASTIN, ID_TEMPER, ID_SOLIDZ, ID_ALEFOR
  use def_master,       only: current_code, modul
  implicit none
  type(MUI_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  ! DIRIC
  CPLNG%code_i       =  2_ip        !< CODEi                                   !< 2016MAR26 1 -> 2 
  CPLNG%module_i     =  ID_TEMPER   !< MODULEi --> 
  modules_i(2)       = -ID_SOLIDZ   !< MODULEi <--
  modules_i(3)       = -ID_TEMPER                                              !< 2016JUN06  
  CPLNG%fixbo_i      = -1_ip
  !
  ! NEUMA  
  CPLNG%code_j       =  1_ip        !< CODEj                                   !< 2016MAR26 2 -> 1 
  CPLNG%module_j     =  ID_TEMPER   !< MODULEj <-- !< 2016MAR30 
  modules_j(2)       = -ID_ALEFOR   !< MODULEj -->
  modules_j(3)       = -ID_TEMPER                                              !< 2016JUN06    
  CPLNG%fixbo_j      = -1_ip
  !
  ! 'ALGEBRAIC' DIRIC-NEUMA system activation...   !< 2016MAR30
  CPLNG%what_i       =  RESIDUAL    !<   
  CPLNG%what_j       = -RESIDUAL    !< 'physical' (<0) or 'numerical' (>0) coupling  
  ! 
  CPLNG%tolerance    = -1e-4
  n_max_its          =  1_ip       !< 2015JUL17  
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  debug              = .false.
  dynamic            = .false.
  CPLNG%n_dof        =  ndime      !< D.O.F.                                   !< 2016ABR04   
 !CPLNG%n_dof        =  ndime + 1  !< D.O.F.                                   !< 2016JUN06
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef MUI
  if( (CPLNG%app_name(1:5) == 'NEUMA').and.(CPLNG%app_name(1:5) == title(1:5)).and.(current_code==CPLNG%code_j) ) then !< 2016MAR26 i -> j 
  else&
  if( (CPLNG%app_name(1:5) == 'DIRIC').and.(CPLNG%app_name(1:5) == title(1:5)).and.(current_code==CPLNG%code_i) ) then !< 2016MAR26 j -> i 
  else
    print *, "[", trim(title),"] ERROR: 'DIRIC' or 'NEUMA' not found!!", " app_name:'", trim(CPLNG%app_name),"' title:'", trim(title) , "' "
    print *, "or DIRIC==2, NEUMA==1 CODE=", current_code
    stop
  endif
#endif 
  !
  ! 
  if( .not.( (CPLNG%what_i==RESIDUAL).and.(CPLNG%what_j/=RESIDUAL) ) ) then
    print *, "[", trim(title),"] ERROR: ", "[mui_driver_init_fsi] SET 'what_i==RESIDUAL and what_j/=RESIDUAL' "
    print *, " 'what_i==RESIDUAL and what_j/=RESIDUAL' <-> DIRIC: kfl_bvnat=0, kfl_react=1; NEUMA: kfl_bvnat=1, kfl_react=0"
    stop
  endif
  !
  !-----------------------------------------------------------------------||---!
  if(CPLNG%code_i==current_code) CPLNG%current_fixbo = CPLNG%fixbo_i
  if(CPLNG%code_j==current_code) CPLNG%current_fixbo = CPLNG%fixbo_j
  !
  CPLNG%current_code = current_code !< *.dat CODE: ID_CODE 
  current_code       = 1_ip         !< trick!! -> PAR_COLOR_COMMUNICATORS: CODE NUMBER EXCEED NUMBER OF CODES 
  mcoup              = 0_ip         !< evoid cou_turnon 
  !
  if(IMASTER) print*, "[mui_driver_init_fsi]"
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine mui_driver_sendrecv(CPLNG, current_when, current_task)
  use def_master,           only: iblok, ittim, itcou 
  use def_master,           only: modul, current_code
  use def_master,           only: nblok
  use def_master,           only: namod, mmodu
  use def_master,           only: ITASK_INIUNK, ITASK_TURNOF
  use def_master,           only: ITASK_TIMSTE, ITASK_ENDSTE 
  use def_master,           only: ITASK_BEGZON, ITASK_ENDZON 
  use def_master,           only: ITASK_AFTER,  ITASK_BEFORE 
  use def_master,           only: ITASK_BEGSTE, ITASK_CONBLK 
  use def_master,           only: ID_KERMOD, ITASK_DOITER, ITASK_CONCOU 
  use def_master,           only: ID_NASTIN, ID_ALEFOR, ID_SOLIDZ, ID_TEMPER
  !
  use def_master,           only: dtinv
  use def_master,           only: displ
  use def_domain,           only: coord
  !
  !< ITERATIONs
  use def_coupli,           only: coupling_driver_iteration 
  !
  implicit none
  type(MUI_COUPLING), intent(inout) :: CPLNG
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  integer(ip),  intent(in)  :: current_when
  integer(ip),  intent(in)  :: current_task
  !
 !integer(ip), parameter   :: n_times = 6
  integer(ip)   :: n_dirichlet(3), n_neumann(3)
  integer(ip)              :: itime 
  integer(ip)              ::       now(8_ip) 
  integer(ip)              ::  the_time(8_ip,n_times) = -1_ip 
  logical(ip)              :: sendrecv(n_times) = .false. 
  !
  character(16)            :: saux(8_ip) = ' '
  character(64)            :: sms        = ' '
  character( 4), parameter :: frmt = '(I2)'
  !
  MUI_SMS = '+-+-+-+-+-+-+-+-'
!print*, "mui_driver_sendrecv"
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  !  ittim: Current time step
  !  itcou: Current global iteration
  !  iblok: Current block
  !  is_when = itcou==micou(iblok)
  !
  now = (/ CPLNG%current_code, iblok, modul, current_task, current_when, ittim, -1_ip, -1_ip /) 
  !
  !-----------------------------------------------------------------------||---!
  write(saux(1), frmt) CPLNG%current_code
  write(saux(2), frmt) iblok
  write(saux(3), frmt) ittim
  write(saux(4), frmt) coupling_driver_iteration(iblok)
  !
  if(current_when==ITASK_BEFORE) then     
    saux(5) = "+"
    saux(6) = "-"
  else& 
  if(current_when==ITASK_AFTER ) then 
    saux(5) = "-"
    saux(6) = "+"
  endif 
  !
  sms = "'"//trim(saux(1))//&
        "."//namod(modul)//&
        "."//trim(saux(5))//name_task(current_task)//trim(saux(6))// & 
       !"."//name_when(current_when)//&
        ".B"//trim(saux(2))// & 
        ".T"//trim(saux(3))// &
        ".I"//trim(saux(4))//"'"
  !-----------------------------------------------------------------------||---!
  if(nblok>1) then 
    print *, "[mui_driver_sendrecv] ", sms 
    print *, "[mui_driver_sendrecv] ", "nblok==1 !!" 
    call runend("EXIT!!")
  endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------| ITERATIONS |---!
  !   +
  !   |_Alya                                       
  !     |_call Turnon()                            
  !     |_call Iniunk()                             
  !     |_time: do while
  !       |_call Timste()                (-1+)
  !       |_reset: do                    
  !         |_call Begste()              (+2-)  
  !           |_block: do while          
  !             |_coupling: do while     
  !               |_call Begzon()        (+4-)  AFTER<-never INTO the module: into 'coupli/mod_coupling_driver.f90' 
  !                                                                           add 'if(current_when==ITASK_AFTER) call Temper(-1_ip)'
  !               |_modules: do while                           
  !                 |_call Doiter()      
  !                 |_call Concou()      (+7-)
  !               |_call Endzon()        (-5+) 
  !             |_call Conblk()          (-6+) 
  !       |_call Endste()                                     
  !   |_call Turnof()                    (-3+)
  !
  !-----------------------------------------------------------------------||---!
  if(CPLNG%current_code==CPLNG%code_i) then 
    the_time(:,1_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_TIMSTE, ITASK_AFTER,  ittim, -1_ip, -1_ip /)   
    the_time(:,2_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_BEGSTE, ITASK_BEFORE, ittim, -1_ip, -1_ip /)   
    the_time(:,3_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_TURNOF, ITASK_AFTER,  ittim, -1_ip, -1_ip /) 
    the_time(:,4_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_BEGZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /) 
    the_time(:,5_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /) 
    the_time(:,6_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_CONBLK, ITASK_BEFORE, ittim, -1_ip, -1_ip /) 
    !
    the_time(:,7_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /) 
    the_time(:,8_ip)  = (/ CPLNG%current_code, iblok,   modules_i(2), ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)  !< 2016MAR26  
    the_time(:,9_ip)  = (/ CPLNG%current_code, iblok,   modules_i(3), ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)  !< 2016JUN06  
  endif 
  if(CPLNG%current_code==CPLNG%code_j) then 
    the_time(:,1_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_TIMSTE, ITASK_AFTER,  ittim, -1_ip, -1_ip /)   
    the_time(:,2_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_BEGSTE, ITASK_BEFORE, ittim, -1_ip, -1_ip /)   
    the_time(:,3_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_TURNOF, ITASK_AFTER,  ittim, -1_ip, -1_ip /) 
    the_time(:,4_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_BEGZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /) 
    the_time(:,5_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /) 
    the_time(:,6_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_CONBLK, ITASK_BEFORE, ittim, -1_ip, -1_ip /) 
    !
    the_time(:,7_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_BEGZON, ITASK_BEFORE, ittim, -1_ip, -1_ip /) 
    the_time(:,8_ip)  = (/ CPLNG%current_code, iblok,   modules_j(2), ITASK_BEGZON, ITASK_BEFORE, ittim, -1_ip, -1_ip /)  !< 2016MAR26 
    the_time(:,9_ip)  = (/ CPLNG%current_code, iblok,   modules_j(3), ITASK_BEGZON, ITASK_BEFORE, ittim, -1_ip, -1_ip /)  !< 2016JUN06   
  endif 
  !
  sendrecv = (/ (all(the_time(:,itime)==now), itime=1,n_times) /)
  !
  do itime = 1,n_times
    if( sendrecv(itime) ) then 
      if(debug) print *, " [", trim(title),"] ",  trim(sms)   
    endif 
  enddo
  !
  MUI_SENDRECV = sendrecv                                                      !< 2016MAR29 
  if(any(sendrecv)) MUI_SMS = trim(sms)
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-------------------------------------------------------------| -TIMSTE+ |---!
#ifdef MUI
  if( sendrecv(1_ip) ) then 
!    call mui_localize(CPLNG)                                               
  endif 
  !-------------------------------------------------------------| +BEGSTE- |---!
  if( sendrecv(2_ip) ) then 
    call mui_driver_begste()
  endif 
  !-------------------------------------------------------------| -TURNOF+ |---!
  if( sendrecv(3_ip) ) then 
!    call mui_plepp_compare_dtinv(dtinv) 
  endif 
  !-------------------------------------------------------------| -BEGZON+ |---! see ../master/Begzon.f90
  if( sendrecv(4_ip) ) then 
    call mui_driver_begzon() 
  endif 
  !-------------------------------------------------------------| -ENDZON+ |---! see ../master/Endzon.f90
  if( sendrecv(5_ip) ) then !.and.exchanged ) then 
    call mui_driver_endzon() 
  endif 
  !-------------------------------------------------------------| +CONBLK- |---! see ../master/Conblk.f90
  if( sendrecv(6_ip) ) then 
    call mui_driver_conblk() 
  endif
  !-------------------------------------------------------------| EXCHANGE |---!
  if( any(sendrecv(7:n_times)) ) then 
    !< AFTER<-never INTO the module !! 
    if(current_when==ITASK_AFTER) then                                            
      if( modul==ID_SOLIDZ ) call Solidz( -1_ip )
      if( modul==ID_NASTIN ) call Nastin( -1_ip )
      if( modul==ID_ALEFOR ) call Alefor( -1_ip )
      if( modul==ID_TEMPER ) call Temper( -1_ip )
    endif
    !
  endif
  !-----------------------------------------------------------------------||---!  
#endif 
  !-----------------------------------------------------------------------||---!  
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine mui_localize(CPLNG) 
  use mod_mui,              only: mui_exchange 
  implicit none
  type(MUI_COUPLING), intent(inout) :: CPLNG
  !
  real(rp), pointer       :: aux(:,:) => null()
  integer(ip)  :: n_recv  
  ! 
  !-----------------------------------------------------------------------||---!
  n_recv = 0 
  if( INOTMASTER ) n_recv = npoin  

  allocate( aux(ndime,n_recv) )
  aux = huge(0_ip)

  !-----------------------------------------------------------------------||---!
  call mui_exchange( CPLNG, &  
                     nameij="XXXX", tij=0.0_rp, varij=aux(1_ip,1:n_recv), & 
                     nameji="YYYY", tji=0.0_rp, varji=aux(2_ip,1:n_recv)  & 
                   )

  !-----------------------------------------------------------------------||---!
  deallocate( aux )
  !-----------------------------------------------------------------------||---!
  if(debug) & 
    print *, "[mui_localize]"
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine mui_driver_exchange02( CPLNG, debug) 
  implicit none 
  type(MUI_COUPLING), intent(inout) :: CPLNG
  logical(ip), optional,  intent(in   ) :: debug
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine mui_driver_init( )
#ifdef MUI
  use mod_mui,              only: MUI_CPLNG
  implicit none
  !-----------------------------------------------------------------------||---!
  call mui_driver_init_fsi( MUI_CPLNG )
  !-----------------------------------------------------------------------||---!
#endif
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine mui_driver_set_mesh( CPLNG )
  implicit none
  type(MUI_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



!-------------------------------------------------------------------------||---!
!----------------------------------------------------------------| PRIVATE |---!
!-------------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine mui_driver_begste() 
  use def_master,  only: iblok
  use def_coupli,  only: coupling_driver_iteration 
  use def_coupli,  only: coupling_driver_number_couplings
  use def_coupli,  only: coupling_driver_max_iteration
  use def_coupli,  only: max_block_cou
  use def_coupli,  only: kfl_gozon
  use def_master,  only: kfl_gocou 
  !
  use def_master,        only: dtinv, cutim, dtime  
  !
  implicit none 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  coupling_driver_iteration(1:max_block_cou) = 0 
  coupling_driver_number_couplings(iblok)    = 1 
  coupling_driver_max_iteration(iblok)       = n_max_its  
  !-----------------------------------------------------------------------||---!
#ifdef MUI
!  call mui_plepp_compare_dtinv(dtinv) 
#endif 
  cutim  = cutim - dtime       
  call setgts(ITASK_TIMSTE)             
  call livinf(201_ip, ' ',1_ip)
  !
  if(debug) print*, "[mui_driver_begste]"
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine mui_driver_begzon() 
  use def_master,    only : iblok
  use def_master,    only : mmodu
  use def_master,    only : lmord
  use def_master,    only : itinn
  use def_coupli,    only : coupling_driver_iteration
  use def_coupli,    only : coupling_driver_number_couplings
  implicit none 
  integer(ip) :: iorde,imodu
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
     if( coupling_driver_number_couplings(iblok) /= 0 .and. coupling_driver_iteration(iblok) == 0 ) then
        call livinf(-6_ip,'ZONAL COUPLING FOR BLOCK ', iblok)
     end if
     !
     ! Put inner iterations to zero
     !
     do iorde = 1,mmodu
        imodu = lmord(iorde,iblok)
        itinn(imodu) = 0
     end do
     !
     ! Iteration counter
     !
     coupling_driver_iteration(iblok) = coupling_driver_iteration(iblok) + 1
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine mui_driver_endzon() 
  use def_master,  only: iblok
  use def_master,  only: kfl_gocou 
  use def_coupli,  only: coupling_driver_iteration 
  use def_coupli,  only: coupling_driver_number_couplings
  use def_coupli,  only: coupling_driver_max_iteration
  use def_coupli,  only: kfl_gozon
  implicit none 
  integer(ip) :: iorde,imodu
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
     if( coupling_driver_number_couplings(iblok) /= 0 ) then
        !
        kfl_gozon = 0
        if( coupling_driver_iteration(iblok) >= coupling_driver_max_iteration(iblok) ) then
 kfl_gozon = 0 !< kernel/coupli/mod_couplings.f90 
!          return
        else
! if( resid_cou(1,icoup) > coupling_driver_tolerance(iblok) ) kfl_gozon = 1 !< kernel/coupli/mod_couplings.f90 
          kfl_gozon = 1
        endif  
        !
        if( kfl_gozon == 1 ) kfl_gocou = 1
     else
        kfl_gozon = 0
     end if
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine mui_driver_conblk() 
  use def_master,  only: iblok
  use def_coupli,  only: coupling_driver_number_couplings
  use def_coupli,  only : coupling_driver_iteration
  implicit none 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( coupling_driver_number_couplings(iblok) /= 0 ) then
    call livinf(-13_ip,'END ZONAL COUPLING: ', coupling_driver_iteration(iblok) )
  end if
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
!------------------------------------------------------------------| TOOLs |---!
!-------------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
end module mod_mui_driver 
!==============================================================================!
!==============================================================================!

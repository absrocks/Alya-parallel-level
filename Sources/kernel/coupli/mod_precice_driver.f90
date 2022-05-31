!==============================================================================!
!  I am your father...
!
!< 2015Jun02 -> created    
!
!==============================================================================!
  !-----------------------------------------------------------------------||---!
  !   + current_code                                      ___________current_task 
  !   |_Alya                                       ______|_____
  !     |_call Turnon()                            ITASK_TURNON  02  
  !     |_call Iniunk()                            ITASK_INIUNK  03 
  !     |_time: do while
  !       |_call Timste()                          ITASK_TIMSTE  04 
  !       |_reset: do 
  !       | |_call Begste()                        ITASK_BEGSTE  05 
  !       |    |_block: do while                                     
  !       |       |_coupling: do while                              
  !       |         |_call Begzon()                ITASK_BEGZON  19   _
  !       |         |_modules: do while                              / TASK_BEGITE  14 
  !       |           |_call Doiter()              ITASK_DOITER  06-|  
  !       |           |_call Concou()              ITASK_CONCOU  07  \_ITASK_ENDITE 15 
  !       |         |_call Endzon()                ITASK_ENDZON  20
  !       |       |_call Conblk()                  ITASK_CONBLK  08 
  !       |_call Endste()                          ITASK_ENDSTE  10 
  !   |_call Turnof()                              ITASK_TURNOF  13   
  !          __
  ! BLOCK 3_   | 
  !   1 X   |  |--current_block  -> CPLNG%blocks_list 
  !   2 Y Z |  |     
  !   3 W  _|-----current_module -> CPLNG%moduls_list 
  ! END_BLOCK__| 
  !
  !-----------------------------------------------------------------------||---!
  !
  ! <code, block, modul, task, when, send|recv>
  !
  !-----------------------------------------------------------------------||---!
module  mod_precice_driver
  use def_parame,           only: ip, rp
  use def_master,           only: inotmaster, imaster, ITASK_TIMSTE
  use def_domain,           only: coord, mnode, ndime, npoin 
  use def_domain,           only: ltype, lnods
  use mod_precice,          only: PRECICE_COUPLING 
  use def_kintyp,           only: soltyp
  use def_master,           only: momod, modul
  use mod_messages, only : livinf
#ifdef PRECICE  
#endif
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  implicit none 
  !
  logical(ip), parameter  :: debug = .true.
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
  integer(ip), parameter   :: n_times = 7
  logical(ip)              :: CNT_SENDRECV(n_times) = .false. 
  character(64)            :: CNT_SMS               = ' '
  type(PRECICE_COUPLING), SAVE   :: PRCC  
  !
  real(rp)                 :: n_max_its
  logical(ip)              :: dynamic
  !
  logical(ip) :: residual = .false. 
  !
  type(soltyp), pointer    :: solve(:)
  !
  private
    public :: PRCC 
    public :: precice_driver_sendrecv
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  !
  !
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  contains
!-------------------------------------------------------------------------||---!
!-----------------------------------------------------------------| PUBLIC |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine precice_driver_sendrecv(CPLNG, current_when, current_task)
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
  !
  use def_master,           only: dtinv
  use def_master,           only: displ
  use def_domain,           only: coord
  !
  !< ITERATIONs
  use def_coupli,           only: coupling_driver_iteration 
  !
  implicit none
  type(PRECICE_COUPLING), intent(inout) :: CPLNG
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  integer(ip),  intent(in)  :: current_when
  integer(ip),  intent(in)  :: current_task
  !
  integer(ip)              :: n_dirichlet(3), n_neumann(3)
  integer(ip)              :: itime 
  integer(ip)              ::       now(8_ip) 
  integer(ip)              ::  the_time(8_ip,n_times) = -1_ip 
  logical(ip)              :: sendrecv(n_times) = .false. 
  !
  character(16)            :: saux(8_ip) = ' '
  character(64)            :: sms        = ' '
  character( 4), parameter :: frmt = '(I2)'
  !
  CNT_SMS = '+-+-+-+-+-+-+-+-'
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
    print *, "[precice_driver_sendrecv] ", sms 
    print *, "[precice_driver_sendrecv] ", "nblok==1 !!" 
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
    the_time(:,6_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_CONBLK, ITASK_BEFORE,  ittim, -1_ip, -1_ip /) 
    !
    the_time(:,7_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_ENDZON, ITASK_AFTER, ittim, -1_ip, -1_ip /) 
  endif 
  if(CPLNG%current_code==CPLNG%code_j) then 
    the_time(:,1_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_TIMSTE, ITASK_AFTER,  ittim, -1_ip, -1_ip /)   
    the_time(:,2_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_BEGSTE, ITASK_BEFORE, ittim, -1_ip, -1_ip /)   
    the_time(:,3_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_TURNOF, ITASK_AFTER,  ittim, -1_ip, -1_ip /) 
    the_time(:,4_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_BEGZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /) 
    the_time(:,5_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /) 
    the_time(:,6_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_CONBLK, ITASK_BEFORE,  ittim, -1_ip, -1_ip /) 
    !
    the_time(:,7_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_BEGZON, ITASK_BEFORE, ittim, -1_ip, -1_ip /) 
  endif 
  !
  sendrecv = (/ (all(the_time(:,itime)==now), itime=1,n_times) /)
  !
  do itime = 1,n_times
    if( sendrecv(itime) ) then 
      if(debug) print *,"[precice_driver_sendrecv] ",  trim(sms)   
    endif 
  enddo
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-------------------------------------------------------------| -TIMSTE+ |---!
#ifdef PRECICE  
  if( sendrecv(1_ip) ) then 
    !--------------------------------------------------------------| here? |---!
!    call precicef_advance(CPLNG%tsl)
    !---------------------------------------------------------------------||---!
  endif 
  !-------------------------------------------------------------| +BEGSTE- |---!
  if( sendrecv(2_ip) ) then 
    call precice_driver_begste()
  endif 
  !-------------------------------------------------------------| -TURNOF+ |---!
  if( sendrecv(3_ip) ) then 
!    call precice_compare_dtinv(dtinv) 
  endif 
  !-------------------------------------------------------------| -BEGZON+ |---!
  if( sendrecv(4_ip) ) then 
    call precice_driver_begzon() 
  endif 
  !-------------------------------------------------------------| -ENDZON+ |---!
  if( sendrecv(5_ip) ) then 
    call precice_driver_endzon() 
  endif 
  !-------------------------------------------------------------| +CONBLK- |---!
  if( sendrecv(6_ip) ) then 
    call precice_driver_conblk() 
    !-----------------------------------------------------------| or here? |---!
    !---------------------------------------------------------------------||---!
  endif 
#endif 
  !-----------------------------------------------------------------------||---!  
  CNT_SENDRECV = sendrecv
  if(any(sendrecv)) CNT_SMS = trim(sms)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
!----------------------------------------------------------------| PRIVATE |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine precice_driver_begste() 
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
#ifdef PRECICE  
!
!  call precicef_advance(CPLNG%tsl)
!  dtinv = 1.0/CPLNG%tsl  
!   
#endif 
  cutim  = cutim - dtime       
  call setgts(ITASK_TIMSTE)             
  call livinf(201_ip, ' ',1_ip)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine precice_driver_begzon() 
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
  subroutine precice_driver_endzon() 
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
          return
        else
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
  subroutine precice_driver_conblk() 
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
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
end module mod_precice_driver   
!==============================================================================!
!==============================================================================!

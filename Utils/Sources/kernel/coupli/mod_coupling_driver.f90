!------------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @name    Coupling functions
!> @file    mod_coupling_driver.f90
!> @author  Guillaume Houzeaux
!> @date    11/06/2014
!> @brief   Driver for coupling
!> @details Driver for coupling
!> @{
!------------------------------------------------------------------------

module mod_coupling_driver
  use def_kintyp,         only : ip,rp,lg
  use def_master,         only : intost
  use def_master,         only : modul
  use def_master,         only : ittim
  use def_master,         only : ID_NASTIN
  use def_master,         only : ID_TEMPER
  use def_master,         only : ID_NASTAL,IMASTER,ITASK_INIUNK
  use def_master,         only : ID_SOLIDZ
  use def_master,         only : ID_ALEFOR
  use def_master,         only : ID_PARTIS
  use def_master,         only : ID_KERMOD
  use def_master,         only : ID_TURBUL
  use def_master,         only : ID_INSITU
  use def_master,         only : ID_EXMEDI
  use def_master,         only : ID_SOLFE2
  use def_master,         only : ITASK_AFTER
  use def_master,         only : namod
  use mod_parall,         only : I_AM_IN_COLOR
  use def_coupli,         only : coupling_type
  use def_coupli,         only : BETWEEN_ZONES
  use def_coupli,         only : mcoup
  use def_kermod,         only : kfl_timeline
  use mod_ker_timeline,   only : ker_timeline
  use mod_commdom_driver, only: CNT_CPLNG, commdom_driver_sendrecv   !< 2016MAR30  
  use mod_measurements,   only : measurements_driver                 !< 2017JAN07  
  use mod_coupling_timer, only : coupling_timer_driver               !< 2017ABR06  
  use mod_cantera,        only : cantera_driver                      !< 2018NOV12   
  use mod_messages,       only : livinf
  use mod_alya2dlb,       only : alya2dlb_DLB_Barrier
  
  implicit none

contains

  subroutine COU_DRIVER(current_when,current_task)
    integer(ip),  intent(in)  :: current_when 
    integer(ip),  intent(in)  :: current_task
    integer(ip)               :: icoup
    integer(ip)               :: ITASK_COUPL
    integer(ip)               :: module_source
    integer(ip)               :: module_target
    integer(ip)               :: color_source
    integer(ip)               :: color_target
    logical(lg)               :: i_compute_and_send
    logical(lg)               :: i_recv_and_assemble

    !
    ! Loop over couplings
    !
    !
#ifdef COMMDOM 
    !
    call   measurements_driver( current_when, current_task ) !< 2017JAN07  
    call coupling_timer_driver( current_when, current_task ) !< 2017ABR06  
    call commdom_driver_sendrecv( CNT_CPLNG, current_when, current_task )  
    ! 
#elif MUI 
    !< 2016Ago23  
#else            
    ! 
    call cantera_driver( current_when, current_task )  
    !
    do icoup = 1,mcoup

       if( coupling_type(icoup) % kind == BETWEEN_ZONES ) then

          module_source = coupling_type(icoup) % module_source
          module_target = coupling_type(icoup) % module_target
          color_source  = coupling_type(icoup) % color_source
          color_target  = coupling_type(icoup) % color_target
          ITASK_COUPL   = icoup + 1000
          !
          ! Should I stay or should I go
          !
          ! Only activate coupling if the current time step is a multiple of the coupling frequency 
          !
          if(    current_task == coupling_type(icoup) % task_compute_and_send .and. &
               & current_when == coupling_type(icoup) % when_compute_and_send .and. &
               & I_AM_IN_COLOR(color_source).and. modul == module_source      .and. &
               & mod( ittim,coupling_type(icoup) % frequ_send ) == 0_ip )   then
             i_compute_and_send  = .true.
          else
             i_compute_and_send  = .false.
          end if

          if(    current_task == coupling_type(icoup) % task_recv_and_assemble .and. &
               & current_when == coupling_type(icoup) % when_recv_and_assemble .and. &
               & I_AM_IN_COLOR(color_target).and. modul == module_target       .and. &
               & mod( ittim,coupling_type(icoup) % frequ_recv ) == 0_ip)     then
             i_recv_and_assemble = .true.
          else
             i_recv_and_assemble = .false.
          end if
          !
          ! Call corresponding module (should call the plugin of the module_source/target and that's it)
          ! 
          if( module_source == modul .or. module_target == modul ) then
             if(  ( i_compute_and_send  .and. .not. i_recv_and_assemble ) .or. &
                & ( i_recv_and_assemble .and. .not. i_compute_and_send  ) ) then
                !
                ! I am source or target
                ! 
                if( i_compute_and_send ) then
                   call livinf(-14_ip,'-> SEND '//trim(coupling_type(icoup) % variable)//' AS A SOURCE FOR COUPLING ',icoup)
                   call ker_timeline(module_source,'INI_COUPLING',icoup,module_source)
                else if( i_recv_and_assemble ) then
                   call livinf(-14_ip,'<- RECV '//trim(coupling_type(icoup) % variable)//' AS A TARGET FOR COUPLING ',icoup)
                end if
                !
                ! DLB Barrier, to calm down MPI...
                !
#ifdef ALYA_DLB_BARRIER   
                if( modul > 0 ) then
                   if( IMASTER ) call messages_live(namod(modul)//': BEFORE DLB BARRIER')
                   ierror = alya2dlb_DLB_Barrier()
                   if( IMASTER ) call messages_live(namod(modul)//': BEFORE DLB BARRIER')
                end if
#endif

                select case ( modul )
                   
                case ( ID_NASTIN )
                   call Nastin(ITASK_COUPL) 
                case ( ID_TEMPER )
                   call Temper(ITASK_COUPL)
                case ( ID_NASTAL )
                   call Nastal(ITASK_COUPL)
                case ( ID_SOLFE2 )
                   call Solfe2(ITASK_COUPL)
                case ( ID_SOLIDZ ) 
                   call Solidz(ITASK_COUPL)
                case ( ID_TURBUL ) 
                   call Turbul(ITASK_COUPL)
                case ( ID_ALEFOR ) 
                   call Alefor(ITASK_COUPL)
                case ( ID_PARTIS )
                   call Partis(ITASK_COUPL)
                case ( ID_KERMOD ) 
                   call Kermod(ITASK_COUPL)
                case ( ID_EXMEDI ) 
                   call Exmedi(ITASK_COUPL)
                case ( ID_INSITU ) 
                   call Insitu(ITASK_COUPL)
                case default
                   call runend('COU_DRIVER: MODULE NOT CODED')
                   
                end select

             else if( i_compute_and_send .and. i_recv_and_assemble ) then

                call runend('COU_DRIVER: WE ARE IN TROUBLE')

             end if

          end if

          if( i_compute_and_send ) then
             call ker_timeline(module_source,'END_COUPLING',icoup,module_target)
          end if

       end if

    end do
    ! 
#endif 

  end subroutine COU_DRIVER


end module mod_coupling_driver
!> @} 
!-----------------------------------------------------------------------

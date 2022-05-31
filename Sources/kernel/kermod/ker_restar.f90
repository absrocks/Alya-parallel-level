subroutine ker_restar(itask)
  !------------------------------------------------------------------------
  !****f* Nastin/ker_restar
  ! NAME 
  !    ker_restar
  ! DESCRIPTION
  !    This routine:
  !    ITASK = 1 ... Reads the initial values from the restart file
  !            2 ... Writes restart file
  ! USES
  ! USED BY
  !    ker_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_kermod
  use mod_postpr
  use mod_memchk
  use mod_ker_arrays,     only : ker_arrays
  use mod_communications, only : PAR_BROADCAST
  use mod_messages, only : messages_live

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: kfl_gores
  integer(ip)             :: ifunc, ndxs, max_windk_systems_local
  integer(ip)             :: read_status
  character(len = 5)      :: var_id

  !----------------------------------------------------------------------
  !
  ! Primary arrays
  !
  !----------------------------------------------------------------------
  
  if( itask == ITASK_READ_RESTART ) then
     call ker_arrays('READ RESTART')
  else if( itask == ITASK_WRITE_RESTART ) then
     call ker_arrays('WRITE RESTART')
  end if


  !----------------------------------------------------------------------
  !
  ! Variables
  !
  !---------------------------------------------------------------------- 

  if( itask == ITASK_READ_RESTART ) then
     !
     ! Read restart
     !
     if( INOTSLAVE ) then

        !----------------------------------------------------------------------
        !       
        ! read windkessel params
        !
        ! there are file reading checks insode to be compatible with old restarts. The checks can be removed after some time, when people will move on to the new restarts.
        read(momod(modul) % lun_rstar, IOSTAT=read_status) var_id !windkessel id
        if (read_status/=0) then
            call messages_live('ker_restart: Error reading kernel windkessel restart information. Possibly the restart files are old. If you are NOT using kernel windkessel, you can ignore this warning','WARNING')
        else
            if ( var_id == "WINDK" ) then
                read(momod(modul) % lun_rstar, IOSTAT=read_status) max_windk_systems_local
                if ( max_windk_systems_local /= max_windk_systems ) then
                    call messages_live('ker_restart: max_windk_systems stored in restart != max_windk_systems in compiled alya. Possibly the restart files are old. Using the smallest of the two to read restart.','WARNING')
                end if
            

                do ifunc = 1,min(max_windk_systems_local, max_windk_systems)
                    read(momod(modul) % lun_rstar, IOSTAT=read_status) ndxs
                    if (windk_systems(ifunc) % ndxs  > 0) then

                        if ( windk_systems(ifunc) % ndxs /= ndxs ) then
                            call runend( "Windkessel "//trim(intost(ifunc))//" was saved with "//trim(intost(ndxs))//" derivatives in restart but expected to have "//trim(intost(windk_systems(ifunc) % ndxs))//". Probably dat files changed between restarts." )
                        end if

                        read(momod(modul) % lun_rstar, IOSTAT=read_status) windk_systems(ifunc) % stored_time_step
                        read(momod(modul) % lun_rstar, IOSTAT=read_status) windk_systems(ifunc) % yprev
                        read(momod(modul) % lun_rstar, IOSTAT=read_status) windk_systems(ifunc) % xprev
                        read(momod(modul) % lun_rstar, IOSTAT=read_status) windk_systems(ifunc) % y_out
                        read(momod(modul) % lun_rstar, IOSTAT=read_status) windk_systems(ifunc) % x_in
                    end if

                    if (read_status/=0) then
                        call runend( "ker_restart: Something failed while reading restart for Windkessel "//trim(intost(ifunc)) )
                    end if

                end do

            else
                call runend('ker_restart: Error reading kernel windkessel restart information. Instead of windkessel parameters found '//var_id)                
            end if
        end if
        !       
        ! read windkessel params
        !
        !----------------------------------------------------------------------


     end if
     if( IPARALL ) then
        do ifunc = 1,max_windk_systems
            if (windk_systems(ifunc) % ndxs  > 0) then
                call PAR_BROADCAST ( windk_systems(ifunc) % stored_time_step )
                call PAR_BROADCAST ( windk_systems(ifunc) % yprev )
                call PAR_BROADCAST ( windk_systems(ifunc) % xprev )
                call PAR_BROADCAST ( windk_systems(ifunc) % y_out )
                call PAR_BROADCAST ( windk_systems(ifunc) % x_in )
            end if
        end do
     end if

  else if( itask == ITASK_WRITE_RESTART ) then 
     !
     ! Write restart file
     !
     if( INOTSLAVE ) then

        !----------------------------------------------------------------------
        !  
        ! write windkessel params
        !
        write(momod(modul) % lun_rstar) "WINDK" !windkessel id
        write(momod(modul) % lun_rstar) max_windk_systems

        do ifunc = 1,max_windk_systems
            write(momod(modul) % lun_rstar) windk_systems(ifunc) % ndxs
            if ( windk_systems(ifunc) % ndxs > 0 ) then
                write(momod(modul) % lun_rstar) windk_systems(ifunc) % stored_time_step
                write(momod(modul) % lun_rstar) windk_systems(ifunc) % yprev
                write(momod(modul) % lun_rstar) windk_systems(ifunc) % xprev
                write(momod(modul) % lun_rstar) windk_systems(ifunc) % y_out
                write(momod(modul) % lun_rstar) windk_systems(ifunc) % x_in
            end if
        end do
        !  
        ! write windkessel params
        !
        !----------------------------------------------------------------------

     end if
     
  end if


end subroutine ker_restar
 

!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_par_configure.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   MPI-IO Wrapper
!> @details This module is a wrapper of MPI-IO functions.
!>          It adds error management and makes type management transparent to the user.
!> @}
!-----------------------------------------------------------------------


module mod_mpio_par_configure

    use def_kintyp,                     only : ip,rp,lg,r1p
    use def_master,                     only : comm_data_par
    use def_master,                     only : intost, retost
    use def_master,                     only : IPARSLAVE, IMASTER, INOTMASTER, IPARALL, ISEQUEN
    use def_master,                     only : lun_pos00
    use def_master,                     only : kfl_ptask
    use mod_memory,                     only : memory_alloca, memory_deallo
    use def_domain,                     only : ndime, npoin, nboun, nelem, npoin_2, nelem_2, nboun_2
    use def_parall,                     only : kfl_parseq_par
    use mod_parall,                     only : PAR_PARALLEL_PARTITION
    use mod_messages,                   only : livinf
    use def_mpio
    use mod_mpio_files

#ifndef MPI_OFF
    use mod_communications,             only : PAR_DEFINE_COMMUNICATOR
#endif

    implicit none

    private
#ifndef MPI_OFF
    include 'mpif.h'
#endif

    character(150)                          ::  wherein_world="IN THE WORLD"        !
    character(150)                          ::  wherein_code="IN MY CODE"           !
    character(150)                          ::  wherein_p="IN MPIO WITH MASTER"     !
    character(150)                          ::  wherein_pwm="IN MPIO"               !

    integer(ip)                             ::  dboxc_p(3)      ! #bin boxes per direction in the coarse bin

    integer(ip)                             ::  i



    public                                  ::  mpio_init,&
                                                par_mpio_finalize,&
                                                dboxc_p,compute_dimensions

    contains

      subroutine mpio_init()
        use mod_communications,     only : PAR_BROADCAST

#ifndef MPI_OFF
        call PAR_BROADCAST(kfl_ptask)
#endif
        mpio_memor=0_8

        if (kfl_ptask == 2) then
           kfl_mpio_input=IO_DISABLED
           if (mpio_flag_geometry_export == PAR_MPIO_ON) then
              call runend("MPI-IO MESH EXPORTING IS NOT COMPATIBLE WITH MESH RESTART")
           end if
        else if (mpio_flag_geometry_export == PAR_MPIO_ON) then
           mpio_flag_post = PAR_MPIO_ON
           kfl_mpio_input=IO_CLASSIC
           kfl_mpio_export=1
        else if (mpio_flag_geometry == PAR_MPIO_ON) then
           if (IPARALL) then
              kfl_mpio_input=IO_MPIO_PAR
           else
              kfl_mpio_input=IO_MPIO_SEQ
           end if
        else if (mpio_flag_geometry == PAR_MPIO_FORCE_SEQ) then
           kfl_mpio_input=IO_MPIO_SEQ
        else
           kfl_mpio_input=IO_CLASSIC
        end if

        if (mpio_flag_post == PAR_MPIO_ON) then
           if (IPARALL) then
              kfl_mpio_post=IO_MPIO_PAR
           else
              kfl_mpio_post=IO_MPIO_SEQ
           end if
        else if (mpio_flag_post == PAR_MPIO_FORCE_SEQ) then
           kfl_mpio_post=IO_MPIO_SEQ
        else
           kfl_mpio_post=IO_CLASSIC
        end if

        if (mpio_flag_rst == PAR_MPIO_ON) then
           if (IPARALL) then
              kfl_mpio_rst=IO_MPIO_PAR
           else
              kfl_mpio_rst=IO_MPIO_SEQ
           end if
        else if (mpio_flag_rst == PAR_MPIO_FORCE_SEQ) then
           kfl_mpio_rst=IO_MPIO_SEQ
        else
           kfl_mpio_rst=IO_CLASSIC
        end if

        if (IPARALL .and. mpio_flag_enabled == PAR_MPIO_ON) then
           call par_mpio_init()
        end if
        if (mpio_flag_geometry /= PAR_MPIO_OFF) then
           if (mpio_flag_geometry_read_post == PAR_MPIO_ON) then
              call find_mpio_restart_mesh_name()
              call compute_dimensions()
           else
              call find_mpio_mesh_name()
           end if
        end if

      end subroutine mpio_init

    subroutine par_mpio_init()
#ifndef MPI_OFF
        if (mpio_flag_enabled==PAR_MPIO_ON) then
            call livinf(-4_ip,'INIT MPI-IO COMMUNICATORS',0_ip)
            if (mpio_flag_communicator==PAR_MPIO_COMM_ALL) then
                call par_create_communicators_all()
            else
                call runend('PARALLEL IO COMMUNICATOR STRATEGY UNKNOWN')
            end if
            if (kfl_parseq_par == PAR_PARALLEL_PARTITION) then
                mpio_flag_all_par=PAR_MPIO_ON
            end if
            call livinf(-5_ip,'END INIT MPI-IO COMMUNICATORS',0_ip)
        end if
#endif
    end subroutine

    subroutine compute_dimensions()
      use mod_mpio_seq_io, only : SEQ_FILE_READ_HEADER
      type(mpio_header) :: header
      if (IMASTER .or. ISEQUEN) then
        call SEQ_FILE_READ_HEADER(fil_coord, header)
        npoin=header%lines
        call SEQ_FILE_READ_HEADER(fil_lnods, header)
        nelem=header%lines
        call SEQ_FILE_READ_HEADER(fil_lnodb, header)
        nboun=header%lines
        npoin_2 = npoin
        nelem_2 = nelem
        nboun_2 = nboun
      end if
    end subroutine

    subroutine par_mpio_finalize()
#ifndef MPI_OFF
        use mod_mpio_par_async_io
        implicit none
        call PAR_ASYNCHRONOUS_WRITE_POP_ALL(0_ip)
#endif
    end subroutine

    subroutine par_create_communicators_all()
#ifndef MPI_OFF
        use mod_communications, only : PAR_COMM_SPLIT, PAR_BROADCAST
        use mod_parall,         only : PAR_CODE_SIZE, PAR_MY_CODE_RANK
        use mod_parall,         only : PAR_COMM_MPIO, PAR_COMM_MPIO_WM, PAR_COMM_MPIO_RANK_WM, PAR_COMM_MPIO_WM_SIZE
        use def_master,         only : IIOSLAVE

        implicit none
        integer(ip)         :: icolor, iproc, ndims
        character(100), PARAMETER :: vacal = "par_create_communicators_all"

        PAR_COMM_MPIO_WM_SIZE = int((PAR_CODE_SIZE-1_ip), ip)

        if (PAR_COMM_MPIO_WM_SIZE<2) then
            mpio_flag_enabled=PAR_MPIO_OFF
            mpio_flag_communicator=PAR_MPIO_OFF
            if (mpio_flag_geometry==PAR_MPIO_ON) then
                kfl_mpio_input=IO_MPIO_SEQ
            end if
            if (mpio_flag_post==PAR_MPIO_ON) then
                kfl_mpio_post=IO_MPIO_SEQ
            end if
            if (mpio_flag_rst==PAR_MPIO_ON) then
                kfl_mpio_rst=IO_MPIO_SEQ
            end if
            call livinf(0_ip,"ONLY ONE MPI-IO WORKER: SWITCHING TO SEQUENTIAL VERSION",0_ip)
            return
        end if

        call livinf(0_ip,"NUMBER OF MPI-IO WORKERS (GEOMETRY READING): "//intost(PAR_COMM_MPIO_WM_SIZE),0_ip)
        call livinf(0_ip,"NUMBER OF MPI-IO WORKERS (POST/RST): "//intost(PAR_CODE_SIZE-1),0_ip)
        !
        ! Create the slaves communicator
        !
        IIOSLAVE = .FALSE.
        if (INOTMASTER .and. (PAR_MY_CODE_RANK-1 < PAR_COMM_MPIO_WM_SIZE)) then
            IIOSLAVE = .TRUE.
            !print*, "process: ", PAR_MY_CODE_RANK
        end if

        icolor = 0
        PAR_COMM_MPIO_WM = 0
        PAR_COMM_MPIO    = 0
        
        !
        ! Creates partitioners and master communicator
        !
        icolor = 0
        if( IIOSLAVE .or. IMASTER ) then
           icolor = 1_4
        else
           icolor = 0_4
        end if
        call PAR_COMM_SPLIT(icolor,PAR_COMM_MPIO,iproc,wherein_code)
        
        if( IIOSLAVE ) then
           icolor = 1_4
        else
           icolor = 0_4
        end if
        call PAR_COMM_SPLIT(icolor,PAR_COMM_MPIO_WM,iproc,wherein_code)
        PAR_COMM_MPIO_RANK_WM = int(iproc,ip)
        
#endif
    end subroutine

end module

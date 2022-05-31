!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    reampio.f90
!> @author  houzeaux
!> @date    2018-11-16
!> @brief   Read
!> @details Read MPI IO strategy
!> @}
!-----------------------------------------------------------------------

subroutine reampio()

  use def_parame
  use def_mpio
  use def_master
  use def_inpout
  use def_parall
  use mod_parall
  use mod_ecoute,   only : ecoute
  use mod_messages, only : messages_live
  use mod_messages, only : livinf

  implicit none

  if( IMASTER .or. ISEQUEN ) then
     !
     ! These 3 variables are initialized at the beginning
     ! in order to define them with the --export option
     !
     !mpio_flag_geometry             = PAR_MPIO_OFF              ! Parallel mesh reading (disabled by default)
     !mpio_flag_geometry_export      = PAR_MPIO_OFF              ! Parallel mesh export (disabled by default)
     !mpio_flag_post_merge           = PAR_MPIO_OFF
     mpio_flag_enabled              = PAR_MPIO_OFF              ! Parallel IO
     mpio_flag_geometry_read_post   = PAR_MPIO_OFF              ! Read post process files (disabled by default)
     mpio_flag_post                 = PAR_MPIO_OFF              ! Parallel IO for post-process (disabled by default)
     mpio_flag_post_light           = PAR_MPIO_OFF
     mpio_flag_rst                  = PAR_MPIO_OFF              ! Parallel IO for restarts (disabled by default)
     mpio_flag_autoromio            = PAR_MPIO_OFF              ! Configure romio automatically (disabled by default)
     mpio_flag_synchro              = PAR_MPIO_OFF              ! Synchronous IO calls (disabled by default)
     mpio_flag_collective           = PAR_MPIO_ON               ! Collective IO calls (enabled by default)
     mpio_flag_all_par              = PAR_MPIO_OFF
     mpio_val_asyncbuffer           = -1                        ! Asynchronous calls buffer size:
     mpio_val_merge_block           = -1
     ! -1: unlimited, synchronize only at next step
     ! 0: synchronize at each write, 1+: synchronize when reach the limit
     mpio_val_hybrid_threshold      = 0
     mpio_flag_communicator         = PAR_MPIO_COMM_ALL
     !.md<module>kernel
     !.md<input>case.dat
     !.md<pos>4
     !.md<sec>
     !.md<0># MPI-IO section
     !.md<code>
     !.mdMPI_IO: ON | OFF
     do while( words(1) /= 'MPIIO' )
        call ecoute('READ_MPI_IO',STOP_END_OF_FILE=.false.)
        if( words(1) == 'ENDFI' ) return
     end do
     if( option_not_off('MPIIO') ) then
        mpio_flag_enabled = PAR_MPIO_ON
        mpio_flag_geometry = PAR_MPIO_ON
        mpio_flag_post = PAR_MPIO_ON
        mpio_flag_rst =PAR_MPIO_ON
        call ecoute('READ_MPI_IO')
        do while( words(1) /= 'ENDMP' )
           !.md<1>GEOMETRY: ON | SEQUENTIAL | EXPORT | OFF
           !.md<field>GEOMETRY
           !.md<com>Geometry reading method.
           !.md<com>    -  `ON`:         Read the `mpio` format, in parallel if possible
           !.md<com>    -  `SEQUENTIAL`: Read the `mpio` format, in sequential
           !.md<com>    -  `EXPORT`:     Read the `ASCII` format, and export the mesh to the `mpio` format without running any step
           !.md<com>    -  `OFF`:        Read the `ASCII` format
           if ( words(1)=='GEOME') then
              if( option_not_off('GEOME') ) then
                 mpio_flag_geometry = PAR_MPIO_ON
                 if( words(2)=='SEQUE' ) then
                    mpio_flag_geometry = PAR_MPIO_FORCE_SEQ
                 else if(words(2)=='POST ') then
                    mpio_flag_geometry_read_post = PAR_MPIO_ON
                 else if(words(2)=='EXPOR') then
                    mpio_flag_geometry = PAR_MPIO_OFF
                    mpio_flag_geometry_export = PAR_MPIO_ON
                    if (IPARALL) then
                       mpio_flag_post_merge=PAR_MPIO_ON
                    end if
                 end if
              else
                 mpio_flag_geometry = PAR_MPIO_OFF
              end if
              !.md<1>MERGE: ON | OFF
              !.md<field>MERGE
              !.md<com>Merge the post-process subdomains
           else if ( words(1)=='MERGE') then
              if( option_not_off('MERGE') ) then
                 mpio_flag_post_merge=PAR_MPIO_ON
              end if
              call ecoute('READ_MPI_IO_MERGE')
              do while( words(1) /= 'ENDME' )
                 !.md<2>BLOCK: (int)
                 !.md<com>    -  **BLOCK** (int): Minimum redistribution block size by process, 0 by default
                 if( words(1) == 'BLOCK' ) then
                    mpio_val_merge_block=int(param(1), ip)
                 end if
                 call ecoute('READ_MPI_IO_MERGE')
              end do
              !.md<1>END_MERGE
              !.md<1>AUTOROMIO: ON | OFF
              !.md<field>AUTOROMIO
              !.md<com>Romio auto configuration (experimental)
           else if ( words(1)=='AUTOR') then
              if( option_not_off('AUTOR') ) then
                 mpio_flag_autoromio = PAR_MPIO_ON
              end if
              !.md<1>SYNCHRONOUS: ON | OFF
              !.md<field>SYNCHRONOUS
              !.md<com>Force synchronous file writing, OFF by default. Use this option if the asynchronous writing fails.
           else if ( words(1)=='SYNCH') then
              if( option('SYNCH') ) then
                 mpio_flag_synchro = PAR_MPIO_ON
              end if
              !.md<1>MINIMUM: (int)
              !.md<field>MINIMUM
              !.md<com>Minimal block size en MB to perform MPI-IO calls (sequential I/O if < this value)
           else if( words(1) == 'MINIM' ) then
              mpio_val_hybrid_threshold=param(1)
              !.md<1>BUFFER | ASYNC: (int)
              !.md<com>-  **ASYNC**: Flush asynchronous writes every [int] files.
           else if( words(1) == 'BUFFE' .or. words(1) == 'ASYNC') then
              mpio_val_asyncbuffer=int(param(1), ip)
              !.md<1>POSTPROCESS: ON | SEQUENTIAL | OFF, LIGHT
              !.md<field>POSTPROCESS
              !.md<com>Post-process writing method
              !.md<com>    -  `ON`:         Write to the `mpio` format, in parallel if possible
              !.md<com>    -  `SEQUENTIAL`: Write to the `mpio` format, in sequential
              !.md<com>    -  `OFF`:        Write to the `alyabin` format
              !.md<com>    -  `LIGHT`:      (optional, if `ON` or `SEQUENTIAL`) Only export the mesh main files necessary for post-process. Incompatible with mesh restart!

           else if( words(1) == 'POSTP' ) then
              if( option_not_off('POSTP') ) then
                 mpio_flag_post = PAR_MPIO_ON
                 if( words(2)=='SEQUE' ) then
                    mpio_flag_post = PAR_MPIO_FORCE_SEQ
                 end if
                 if (exists('LIGHT')) then
                    mpio_flag_post_light = PAR_MPIO_ON
                 end if
              else
                 mpio_flag_post = PAR_MPIO_OFF
              end if
              !.md<1>RESTART: ON | SEQUENTIAL | OFF
              !.md<field>RESTART
              !.md<com>Restart writing method
              !.md<com>    -  `ON`:         Write to the `mpio` format, in parallel if possible
              !.md<com>    -  `SEQUENTIAL`: Write to the `mpio` format, in sequential
              !.md<com>    -  `OFF`:        Write to the `alyabin` format
           else if( words(1) == 'RESTA' ) then
              if( option_not_off('RESTA') ) then
                 mpio_flag_rst = PAR_MPIO_ON
                 if( words(2)=='SEQUE' ) then
                    mpio_flag_rst= PAR_MPIO_FORCE_SEQ
                 end if
              else
                 mpio_flag_rst = PAR_MPIO_OFF
              end if
           end if
           call ecoute('par_reapro')
        end do
        !.mdEND_MPI_IO
        !.md</code>
        !
        ! Check errors and warnings
        !
        if (mpio_flag_enabled == PAR_MPIO_ON) then
           call livinf(0_ip,'MPI-IO ENABLED',0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_geometry_export == PAR_MPIO_ON) then
           call livinf(0_ip,'MPI-IO FORMAT MESH WILL BE EXPORTED AND EXECUTION WILL END JUST AFTER',0_ip)
           return
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_geometry == PAR_MPIO_ON) then
           call livinf(0_ip,'PARALLEL MPI-IO MESH READING ENABLED',0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_geometry == PAR_MPIO_FORCE_SEQ) then
           call livinf(0_ip,'SEQUENTIAL MESH READING ENABLED (FORCED)',0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_geometry == PAR_MPIO_ON.and.mpio_flag_geometry_read_post == PAR_MPIO_ON) then
           call livinf(0_ip,'READING FROM POST-PROCESSED MESH',0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_post == PAR_MPIO_ON.and.mpio_flag_post_light == PAR_MPIO_OFF) then
           call livinf(0_ip,'PARALLEL POST-PROCESS IO ENABLED',0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_post_light == PAR_MPIO_ON) then
           call livinf(0_ip,'PARALLEL POST-PROCESS IO ENABLED (LIGHT)',0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_post_merge == PAR_MPIO_ON) then
           call livinf(0_ip,'PARALLEL POST-PROCESS MERGING ENABLED',0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_post == PAR_MPIO_FORCE_SEQ) then
           call livinf(0_ip,'SEQUENTIAL POST-PROCESS IO ENABLED (FORCED)',0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_rst == PAR_MPIO_ON) then
           call livinf(0_ip,'PARALLEL RESTART IO ENABLED',0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_rst == PAR_MPIO_FORCE_SEQ) then
           call livinf(0_ip,'SEQUENTIAL RESTART IO (FORCED)',0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_autoromio == PAR_MPIO_ON) then
           call livinf(0_ip,'ROMIO AUTO CONFIGURATION ENABLED',0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_synchro == PAR_MPIO_ON) then
           call livinf(0_ip,'SYNCHRONOUS MPI-IO WRITE CALLS ENABLED',0_ip)
        else if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_synchro == PAR_MPIO_OFF) then
           call livinf(0_ip,'ASYNCHRONOUS MPI-IO WRITE CALLS ENABLED',0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_communicator == PAR_MPIO_COMM_ALL) then
           call livinf(0_ip,'PARALLEL IO COMMUNICATOR STRATEGY: ALL',0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_communicator /= PAR_MPIO_COMM_ALL) then
           call runend('UNKNOWN PARALLEL IO COMMUNICATOR STRATEGY!')
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_val_hybrid_threshold/=0.0_rp) then
           call livinf(0_ip,'PARALLEL IO MINIMUM SIZE IS '//trim(retost(mpio_val_hybrid_threshold)),0_ip)
        end if
        if (mpio_flag_enabled == PAR_MPIO_ON.and.mpio_flag_synchro == PAR_MPIO_OFF) then
           call livinf(0_ip,'ASYNCHRONOUS IO FILE BUFFER SIZE IS '//trim(intost(mpio_val_asyncbuffer)),0_ip)
        end if
        if (mpio_flag_post_merge == PAR_MPIO_ON) then
           call livinf(0_ip,'MINIMUM MERGING BLOCK SIZE IS '//trim(intost(mpio_val_merge_block)),0_ip)
        end if
     end if

  end if

end subroutine reampio


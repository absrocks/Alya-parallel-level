!> @file    Alya.f90
!! @author  Guillaume Houzeaux
!! @brief   Ayla main
!! @details Ayla main  \n
!!          A L Y A \n
!!          COMPUTATIONAL MECHANICS AND DESIGN  \n
!!          \n
!!          Contact and general info: \n
!!           \n
!!          guillaume.houzeaux@bsc.es   \n
!!          mariano.vazquez@bcs.es  \n
program Alya
  use def_kintyp, only : ip
  use def_master, only : ITASK_INITIA
  use def_master, only : ITASK_ENDTIM
  use def_master, only : ITASK_ENDRUN
  use def_master, only : INOTMASTER
  use def_master, only : kfl_reset
  use def_master, only : kfl_goopt
  use def_master, only : kfl_gotim
  use def_master, only : kfl_goblk
  use def_master, only : kfl_gocou
  use def_master, only : IPARALL
  use def_coupli, only : kfl_gozon
  use mod_mpio_par_async_io, only : PAR_ASYNCHRONOUS_WRITE_POP_ALL
  use def_master

  implicit none
#ifdef ALYA_OMPSS
  integer :: i
  !$omp do schedule(static)
   do i=1, 2
   enddo
  !$omp do schedule(dynamic)
   do i=1, 2
   enddo
#endif
  call Turnon()

  if (IPARALL) then
    call PAR_ASYNCHRONOUS_WRITE_POP_ALL(0_ip)
  end if

  optimization: do while ( kfl_goopt == 1 )

     call Iniunk()

     call Filter(ITASK_INITIA)
     call Output(ITASK_INITIA)

     time: do while ( kfl_gotim == 1 )

        call Timste()

        reset: do
           call Begste()

           block: do while ( kfl_goblk == 1 )

              zone_coupling: do while ( kfl_gozon == 1 )

                 call Begzon()

                 coupling_modules: do while ( kfl_gocou == 1 )
                    call Doiter()
                    call Concou()
                 end do coupling_modules

                 call Endzon()

              end do zone_coupling

              call Conblk()

           end do block
           if( kfl_reset /= 1 ) exit reset

        enddo reset

        call Endste()

        call Filter(ITASK_ENDTIM)
        call Output(ITASK_ENDTIM)

     end do time

     call Doopti()

     call Endopt()

     call Filter(ITASK_ENDRUN)
     call Output(ITASK_ENDRUN)

  end do optimization
#ifdef NINJA
  call lasterrorninja()
#endif
  call Turnof()

end program Alya

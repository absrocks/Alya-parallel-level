!----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @author  J.C. Cajas
!> @date    17/02/2014
!> @brief   Split MPI_COMM_WORLD in code communicators
!> @details Split MPI_COMM_WORLD in code communicators, based on the 
!>          name of the executable the code communicator is stored 
!>          to PAR_COMM_MY_CODE_ARRAY(:). This is done before knowing 
!>          the number of zones and subdomains.
!>
!>                        MPI_COMM_WORLD = PAR_COMM_UNIVERSE              |
!>                                 (THE UNIVERSE)                         |
!>                                       ||                               |
!>                                       \/                               |
!>                                                                        |
!>                                 MPI_COMM_SPLIT                         |
!>                                                                        | 
!>                                //            \\                        |        PERFORMED IN 
!>                               //              \\                       |  par_code_split_universe
!>                                                                        | 
!>                   PAR_COMM_WORLD             REST OF THE WORLD         |
!>              (ALYA's COMMUNICATOR)        (OTHER CODES COMMUNICATOR)   |
!>                        ||                                              | 
!>                        \/                                              |
!> ----------------------------------------------------------------------------------------
!>                  MPI_COMM_SPLIT                                        |
!>                                                                        |
!>                        ||                                              |
!>                        \/                                              |
!>                                                                        |        PERFORMED IN 
!> PAR_COMM_MY_CODE_ARRAY(1) % PAR_COMM_WORLD = PAR_COMM_MY_CODE          |   par_code_split_world
!>                 (par_sencom.f90)                                       |
!>                                                                        |
!>                        ||  MPI_COMM_SPLIT                              |
!>                        \/                                              |
!>                                                                        |
!> PAR_COMM_COLOR_ARRAY(:) % PAR_COMM_WORLD (par_color_communicators.f90) |
!>                                                                        |
!>                                                                        |       
!> @} 
!
!----------------------------------------------------------------------

subroutine par_code_split_world()
  use def_kintyp,         only : ip
  use def_master,         only : IPARALL,kfl_paral
  use mod_parall,         only : PAR_COMM_MY_CODE
  use mod_parall,         only : PAR_COMM_MY_CODE4
  use mod_parall,         only : PAR_COMM_MY_CODE_WM4
  use mod_parall,         only : PAR_COMM_WORLD
  use mod_parall,         only : PAR_WORLD_SIZE
  use mod_parall,         only : PAR_MY_WORLD_RANK
  use mod_parall,         only : PAR_MY_WORLD_RANK_WM
  use mod_parall,         only : PAR_COMM_CURRENT
  use mod_parall,         only : mcode
  use mod_communications, only : PAR_INIT
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
  use mod_communications, only : PAR_GATHER
  use mod_communications, only : PAR_COMM_SPLIT
  use mod_communications, only : PAR_SCATTER
  use mod_communications, only : PAR_BROADCAST
  use mod_communications, only : PAR_LENGTH_INTEGER
  use mod_communications, only : PAR_DEFINE_COMMUNICATOR
  implicit none

  integer(ip)            :: i
  integer(ip),   pointer :: app_id_arra(:)
  integer(ip)            :: app_id,my_new_rank
  character(32)          :: app_type,app_name
  character(32), pointer :: app_arra(:)
  integer(4)             :: PAR_COMM_WORLD4
  integer(4)             :: my_new_rank4,my_new_rank_wm4
  !
  ! Initialize variables and pointers
  !
!   call par_initialize_nullify()
  !
  ! Initialize MPI
  !
!   call PAR_INIT()
  ! 
  ! Nullify pointers
  !
  nullify(app_arra)
  nullify(app_id_arra)
  !
  ! Get the names of the different Alya codes
  !
  call get_command_argument(1_4,app_name)
  !
  ! PAR_INTEGER: Length of integers for MPI communciations (4/8)
  !
  call PAR_LENGTH_INTEGER()

#ifndef MPI_OFF
  !
  ! Get your rank, and the world size. PAR_COMM is MPI_COMM_WORLD
  !
  call PAR_DEFINE_COMMUNICATOR('IN THE WORLD',PAR_COMM_WORLD4)
  PAR_COMM_WORLD = int(PAR_COMM_WORLD4,ip)
  call PAR_COMM_RANK_AND_SIZE(PAR_COMM_WORLD,PAR_MY_WORLD_RANK,PAR_WORLD_SIZE)

  if( PAR_WORLD_SIZE > 1 ) then
     !
     ! Rank 0 gathers the names, compares them and assigns app_id4 to each process
     !
     IPARALL = .true.
     if( PAR_MY_WORLD_RANK == 0 ) then
        allocate(app_id_arra(PAR_WORLD_SIZE))
        allocate(app_arra(PAR_WORLD_SIZE))
     end if

     call PAR_GATHER(app_name,app_arra,'IN THE WORLD')

     if( PAR_MY_WORLD_RANK == 0 ) then
        app_id      =  1
        app_id_arra = -1
        do i = 1_4,PAR_WORLD_SIZE-1_4
           if( trim(app_arra(i)) /= trim(app_arra(i+1)) ) then 
              app_id_arra(i)   = app_id
              app_id_arra(i+1) = app_id + 1
              app_id           = app_id + 1
           else
              app_id_arra(i)   = app_id
              app_id_arra(i+1) = app_id
           endif
        end do
        mcode = app_id
     end if
     !
     ! Maximum number of codes
     !
     call PAR_BROADCAST(mcode,'IN THE WORLD')
     ! 
     ! Rank 0 scatters the app_id_arra4, and the split of COMM_WORLD is performed
     !
     call PAR_SCATTER(app_id_arra,app_id,'IN THE WORLD')
     !
     ! Split communicator
     !
     
     call PAR_COMM_SPLIT(app_id,PAR_COMM_MY_CODE,my_new_rank,'IN THE WORLD')
     
     PAR_COMM_MY_CODE4 = int(PAR_COMM_MY_CODE,4_ip)
     PAR_COMM_CURRENT  = PAR_COMM_MY_CODE 
     !
     ! My code communicator without master
     !
     my_new_rank4 = min(1_4,int(my_new_rank,4))
     call PAR_COMM_SPLIT(my_new_rank4,PAR_COMM_MY_CODE_WM4,my_new_rank_wm4,'IN MY CODE')
     PAR_MY_WORLD_RANK_WM = int(my_new_rank_wm4,ip)
     !
     ! Deallocate arrays 
     !
     if( PAR_MY_WORLD_RANK == 0 ) then
        deallocate(app_id_arra)
        deallocate(app_arra)
     end if

  else

     call PAR_DEFINE_COMMUNICATOR('IN THE WORLD',PAR_COMM_MY_CODE4)
     PAR_COMM_MY_CODE = int(PAR_COMM_MY_CODE4,ip)
     PAR_COMM_CURRENT = PAR_COMM_MY_CODE

  end if

#else
   
  PAR_WORLD_SIZE = 1

#endif
  
end subroutine par_code_split_world

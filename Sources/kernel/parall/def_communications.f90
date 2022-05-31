!-----------------------------------------------------------------------
!> @addtogroup Communication_Toolbox
!> @{
!> @file    def_communications.f90
!> @author  houzeaux
!> @date    2020-05-13
!> @brief   Variables
!> @details Variables needed for communications
!-----------------------------------------------------------------------

module def_communications

  use def_kintyp, only : ip,rp,lg,r1p,i1p
  use def_domain, only : npoin
  use def_domain, only : npoin_2
  use def_domain, only : nelem_2
  use def_domain, only : nboun_2
  use def_master, only : npoi1,kfl_paral
  use def_master, only : IPARALL,IMASTER,ISLAVE,INOTSLAVE,ISEQUEN
  use def_master, only : INOTMASTER
  use def_master, only : comm_data_par
  use def_master, only : current_code,current_zone,current_subd
  use mod_maths,  only : maths_heap_sort
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use mod_memory, only : memory_size
  use mod_memory, only : memory_resize
  use mod_parall, only : par_code_zone_subd_to_color
  use mod_parall, only : PAR_COMM_COLOR
  use mod_parall, only : PAR_COMM_MY_CODE
  use mod_parall, only : PAR_COMM_MY_CODE_WM4
  use mod_parall, only : PAR_COMM_COLOR_ARRAY
  use mod_parall, only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall, only : PAR_MY_CODE_RANK, PAR_MY_WORLD_RANK
  use mod_parall, only : PAR_COMM_WORLD
  use mod_parall, only : PAR_COMM_CURRENT
  use mod_parall, only : PAR_INTEGER
  use mod_parall, only : PAR_COMM_SFC
  use mod_parall, only : PAR_COMM_SFC_WM
  use mod_parall, only : PAR_COMM_MPIO
  use mod_parall, only : PAR_COMM_MPIO_WM
  use mod_parall, only : commd
  use mod_parall, only : color_target
  use mod_parall, only : color_source
  use mod_parall, only : par_memor
  use mod_parall, only : PAR_CODE_SIZE
  use def_parall, only : kfl_order_exchange_par   
  use mod_std

  implicit none
#ifndef MPI_OFF
  include 'mpif.h'
  integer(4)  :: status(MPI_STATUS_SIZE)
#endif

  integer(ip) :: PAR_COMM_NULL      ! NULL MPI communicator
  integer(4)  :: PAR_COMM_NULL4     ! NULL MPI communicator

  public :: PAR_COMM_NULL
  public :: PAR_COMM_NULL4
  
end module def_communications

!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_prepro.f90
!> @author  houzeaux
!> @date    2018-10-18
!> @brief   Parallel preprocess
!> @details Partitioning and redistribution to the slaves
!> @} 
!-----------------------------------------------------------------------

subroutine par_prepro()

  use def_kintyp,                      only : ip
  use def_master,                      only : routp
  use mod_parall,                      only : PAR_PARALLEL_PARTITION
  use mod_parall,                      only : PAR_SEQUENTIAL_PARTITION
  use mod_parall,                      only : PAR_COMM_MY_CODE4
  use def_parall,                      only : kfl_parseq_par 
  use def_parall,                      only : kfl_partition_par
  use def_parall,                      only : boxes_fine_par 
  use def_parall,                      only : boxes_coarse_par 
  use mod_parall,                      only : PAR_SFC
  use mod_par_parallel_partitioning,   only : par_parallel_partitioning
  use mod_par_sequential_partitioning, only : par_sequential_partitioning
  use mod_partition_sfc,               only : partition_sfc_statistics
  use mod_outfor,                      only : outfor

  implicit none 

  call outfor(77_ip)
  
  select case ( kfl_parseq_par )

  case ( PAR_PARALLEL_PARTITION ) 
     !
     ! Parallel partitioning
     !
     call par_parallel_partitioning()

  case ( PAR_SEQUENTIAL_PARTITION ) 
     !
     ! Sequential partitioning
     !
     call par_sequential_partitioning()

  end select
  !
  ! Output statistics when using SFC 
  !  
  if( kfl_partition_par ==  PAR_SFC ) then
     call partition_sfc_statistics(&
          routp(1),routp(2),routp(3),routp(4),routp(5),boxes_fine_par,&
          boxes_coarse_par,COMM4=PAR_COMM_MY_CODE4)
  end if
  call outfor(95_ip)

end subroutine par_prepro

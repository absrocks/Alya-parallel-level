!-----------------------------------------------------------------------
!> @addtogroup Partit
!> @{
!> @file    Partit.f90
!> @author  houzeaux
!> @date    2019-06-11
!> @brief   Partition domain
!> @details Partition domain
!> @} 
!-----------------------------------------------------------------------

subroutine Partit()

  use def_master
  use def_domain
  use def_kermod,           only : kfl_timeline
  use mod_ker_timeline,     only : ker_timeline
  use def_parall,           only : kfl_partition_par
  use def_parall,           only : kfl_parseq_par
  use def_parall,           only : kfl_virfi_par
  use mod_parall,           only : PAR_USING_RANK
  use mod_parall,           only : PAR_PARALLEL_PARTITION
  use def_mpio,             only : kfl_mpio_export
  use mod_domain,           only : domain_memory_deallocate 
  use mod_par_partitioning, only : par_partitioning
  use mod_par_virfil,       only : par_dumbuf
  use mod_performance,      only : performance_outcpu
  use mod_performance,      only : performance_outmem

  implicit none
  real(rp) :: time1,time2

  if( IPARALL ) then
     !
     ! When exporting mesh in parallel, do not partition neither redistribute
     !
     if( kfl_mpio_export == 1 ) then
        kfl_partition_par = PAR_USING_RANK
        kfl_parseq_par    = PAR_PARALLEL_PARTITION
     end if
     call cputim(time1)
     if( kfl_timeline /= 0 ) call ker_timeline(0_ip,'INI_PARTITION_MESH')
     !
     ! Partition mesh
     !
     call par_partitioning()

     call par_errors(2_ip)
     call par_openfi(2_ip)
     !
     ! Output partition info
     !
     if( kfl_virfi_par == 1 ) then                       ! Virtual file
        call par_dumbuf(-1_ip)
        kfl_virfi_par = 0
     end if
     !
     ! Info and possibly end the run
     !
     if( PART_AND_WRITE() ) then                         ! Partition only: end of the run
        call performance_outmem()
        call performance_outcpu() 
        call outlat(2_ip)
        call outlat(3_ip)
        call par_turnof()
        call runend('O.K.!')
     end if
     if( IMASTER .and. .not. READ_AND_RUN() ) then
        call domain_memory_deallocate('ALL MESH')        ! Deallocate Master geometry memory
        call domain_memory_deallocate('LESET')           ! Sets
        call domain_memory_deallocate('LBSET')           ! Sets
        call domain_memory_deallocate('LNSET')           ! Sets
     end if
     if( IMASTER ) then
        call par_memory(4_ip)                            ! Deallocate memory of partition arrays
     end if
     kfl_ptask = 1                                       ! Switch to normal execution
     call vocabu(-1_ip,0_ip,0_ip)
     if( kfl_timeline /= 0 ) call ker_timeline(0_ip,'END_PARTITION_MESH')
     call cputim(time2)
     cpu_start(CPU_MESH_PARTITION) = time2 - time1

  end if

end subroutine Partit

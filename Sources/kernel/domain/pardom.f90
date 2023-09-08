subroutine pardom()
  !-----------------------------------------------------------------------
  !****f* domain/domain
  ! NAME
  !    domain
  ! DESCRIPTION
  !    This is the main routine to partition the domain.
  ! USED BY
  !    Turnon
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_kermod,       only : kfl_timeline
  use mod_ker_timeline, only : ker_timeline
  use def_parall,       only : kfl_partition_par
  use def_parall,       only : kfl_parseq_par
  use mod_parall,       only : PAR_USING_RANK
  use mod_parall,       only : PAR_PARALLEL_PARTITION
  use def_mpio,         only : kfl_mpio_export
  implicit none
  real(rp) :: time1,time2
  !
  ! When exporting mesh in parallel, do not partition neither redistribute
  !
  if( kfl_mpio_export == 1 ) then
     kfl_partition_par = PAR_USING_RANK
     kfl_parseq_par    = PAR_PARALLEL_PARTITION
  end if
  !
  ! Services
  !
  call cputim(time1)
  if( kfl_timeline /= 0 ) call ker_timeline(0_ip,'INI_PARTITION_MESH')
  call Parall(1_ip)                        ! Partition graph
  if( kfl_timeline /= 0 ) call ker_timeline(0_ip,'END_PARTITION_MESH')
  call cputim(time2)
  cpu_start(CPU_MESH_PARTITION) = time2 - time1
  !
  ! temporary writing of 2d mesh
  !
  !call wriaux()

end subroutine pardom

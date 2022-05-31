!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_parallel_efficiency.f90
!> @author  bsc21943
!> @date    2020-03-05
!> @brief   Check parallel efficiency
!> @details Check parallel efficiency
!> @} 
!-----------------------------------------------------------------------

subroutine par_parallel_efficiency()

  use def_master
  use mod_communications
  use mod_alya2talp, only : alya2talp_parallel_efficiency
  use mod_alya2talp, only : alya2talp_MonitoringRegionreset
  implicit none    
  real(rp) :: elapsed_time
  real(rp) :: accumulated_MPI_time
  real(rp) :: accumulated_computation_time
  real(rp) :: sum_accumulated_MPI_time
  real(rp) :: sum_accumulated_computation_time
  real(rp) :: sum_elapsed_time
  real(rp) :: max_accumulated_MPI_time
  real(rp) :: max_accumulated_computation_time
  real(rp) :: max_elapsed_time
  real(rp) :: ave_accumulated_MPI_time
  real(rp) :: ave_accumulated_computation_time
  real(rp) :: ave_elapsed_time
  real(rp) :: PE,PE_rms,LB,LB_rms

  if(mod(ittim,10_ip)==0) then

     call alya2talp_parallel_efficiency(PE,PE_rms,LB,LB_rms,DYNAMIC_ALLOCATION_REGION=.true.)

     if( INOTSLAVE ) then
        write(*,10) &
             npart,PE,PE_rms,LB,LB_rms
     end if

     call alya2talp_MonitoringRegionreset(DYNAMIC_ALLOCATION_REGION=.true.)

  end if

10 format(i6,1x,20(1x,e12.6))

end subroutine par_parallel_efficiency

!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_turnon.f90
!> @author  Guillaume Houzeaux
!> @brief   Turn on Nastin module
!> @details Read data and allocate memory
!> @}
!------------------------------------------------------------------------
subroutine nsi_turnon()

  use def_kintyp
  use def_master
  use def_domain
  use def_nastin
#ifdef OPENACCHHH
  use def_parall,                only : kfl_cores_per_gpu
  use openacc
#endif

  implicit none
#ifdef OPENACCHHH
  integer(ip)  :: gpunum, ngpus
#endif
  !
  ! Initial variables
  !
  call nsi_inivar(0_ip)
  !
  ! Read the physical problem
  !
  call nsi_reaphy()
  !
  ! Read the numerical treatment
  !
  call nsi_reanut()
  !
  ! Read the output strategy
  !
  call nsi_reaous()
  !
  ! Read the boundary conditions
  !
  call nsi_reabcs()
  !
  ! Service: Parall
  !
  call nsi_parall(1_ip)
  !
  ! Modify boundary conditions
  !
  call nsi_inibcs()
  call nsi_updbcs(0_ip)
  !
  ! Initial variables
  !
  call nsi_inivar(1_ip)
  !
  ! Allocate memory
  !
  call nsi_memall()
  !
  ! Warnings and errors
  !
  call nsi_outerr()
  !
  ! Open additional files
  !
  call nsi_openfi(2_ip)
  !
  ! For openacc data exchange cpu-gpu
  ! Perhaps add contiguous to pointers
  !

#ifdef OPENACCHHH
  ngpus = acc_get_num_devices(acc_device_nvidia)
  if(ngpus/=0) then
     if (kfl_cores_per_gpu==0) then
        gpunum = mod(kfl_paral, ngpus)
     else
        gpunum = mod(kfl_paral/kfl_cores_per_gpu, ngpus)
     end if
     call acc_set_device_num(gpunum,acc_device_nvidia)
  else
     ! call acc_set_device_type(acc_device_host)
  end if
  
!!  ngpus = acc_get_num_devices(acc_device_nvidia)
!!  if(ngpus/=0) then
!!     if (kfl_cores_per_gpu==0) then
!!        gpunum = mod(kfl_paral, ngpus)
!!     else
!!        gpunum = mod(kfl_paral/kfl_cores_per_gpu, ngpus)
!!     end if
!!     call acc_set_device_num(gpunum,acc_device_nvidia)
!!  else
!!     call acc_set_device_type(acc_device_host)
!!  end if
!!

  if (kfl_paral /= 0 ) then
    
     if( kfl_savda == 0 ) then
        !$acc enter data copyin(coord,ltype,lnods,lnodb,gravi_nsi)
     else
        !$acc enter data copyin(coord,ltype,lnods,lnodb,gravi_nsi,   &
        !$acc                   elmda_gpvol, elmda_gpcar )
     end if
  end if  
#endif
  
  porfo_nsi = 0.0_rp
  
end subroutine nsi_turnon


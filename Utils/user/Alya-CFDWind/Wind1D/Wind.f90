program Wind
  !-----------------------------------------------------------------------
  !****f* master/zephyr
  ! NAME
  !    Zephyr
  ! DESCRIPTION
  !    This is the main routine of the program. It controls the main flow 
  !    of the program by calling routines that drive specific tasks 
  ! USES
  !    Reapro
  !    Turnon
  !    Begste
  !    Doiter
  !    Endste
  !    Turnof
  !***
  !-----------------------------------------------------------------------
  use def_master
  implicit none
  integer :: kfl_gotim

  call Reapro  ! Reads program data and loads some variables 
    
  ! Gets the mesoscale input variables if GABLS3 type
  ! This allows to compile without read_meso subroutine
  ! and avoid netCDF library
#ifdef _NETCDF_
  ! GABLS3
  if (kfl_case.eq.3) call read_meso
  
#endif
  ! call read_qrad
!  call read_johans
  ! call read_horna_tower_2
!  call read_radiat ! read radiative heat
!  call read_vegeo  ! read geostrophic velocity
  ! Turns on problem
  istep = 0
  ctime = 0.0d0

  call Turnon  
!  if (kfl_trtem) call postpr ! prints first step
  kfl_gotim = 1

  ! Advances in time
  time: do while (kfl_gotim==1)  
     call Begste
     iiter = 0  ! Global iteration
     do while(iiter.lt.miite)   
        iiter = iiter +1
        call Doiter ! solves all variables
     end do
     call Endste     
     ! if (istep.eq.nstep) kfl_gotim =0
     if (istep.ne.1.and.mod(istep, nstep)==0) kfl_gotim =0
  end do time

  ! Writes restart and postprocess info
  call Turnof  

end program Wind


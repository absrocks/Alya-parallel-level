subroutine endste
  !--------------------------------------------------
  !This routine updates unknowns and calls postporcess
  !--------------------------------------------------
  use      def_master
  implicit none


  ! Updates unknowns
  do ipoin = 1, npoin
     veloc(ipoin, 1,2)=veloc(ipoin, 1,1)
     veloc(ipoin, 2,2)=veloc(ipoin, 2,1)
     keyva(ipoin,2)=keyva(ipoin,1)
     epsil(ipoin,2)=epsil(ipoin,1)
     tempe(ipoin,2)=tempe(ipoin,1)
  end do
  
  ! postprocess time step results only for transient problems
  if (kfl_trtem) call postpr
  !.and.mod(istep, stepr)==0) call postpr

  ! Filling netCDF output variables every freq_rst seconds
#ifdef _NETCDF_      ! This allows to compile without read_meso subroutine and avoid netCDF library
     if (kfl_netCDF.and.mod(istep, freq_netcdf)==0) call netCDF_Postpr
#endif 
  
     ! Writes Restart file every freq_rst seconds
!     if (restart_out.and.mod(ctime, freq_rst)==0) call write_restart (Jordi)
     if (restart_out.and.mod(istep, freq_rst)==0) call write_restart

End subroutine endste

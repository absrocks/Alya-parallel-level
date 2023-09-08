subroutine Iniunk()
  !-----------------------------------------------------------------------
  !****f* master/Iniunk
  ! NAME
  !    Iniunk
  ! DESCRIPTION
  !    This routine ask the modules to compute their initial condition
  ! USES
  !    moduls
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  
  use def_master
  use mod_coupling_driver, only : COU_DRIVER
  use mod_communications,  only : PAR_BARRIER
  use mod_ker_noslwa,      only : ker_noslwa
  use mod_outfor,          only : outfor
  use mod_messages,        only : livinf
  use mod_mass_matrix,     only : mass_matrix_consistent

  implicit none
  !
  ! Call modules
  !
  modul = ID_KERMOD
  call COU_DRIVER(ITASK_BEFORE,ITASK_INIUNK)
  call Kermod(-ITASK_INIUNK)
  modul = ID_KERMOD
  call COU_DRIVER(ITASK_AFTER,ITASK_INIUNK)
  !
  ! Weighted consistent mass... Kermod should have been read
  !
  call mass_matrix_consistent(CONSISTENT_MASS=.false.,CONSISTENT_WEIGHTED_MASS=.true.)
  !
  ! Wall exchange strategy
  !
  call ker_waexlo()
!  call ker_extrap_boundary() 
  !
  ! NO Slip Wall law - wall law adding extra viscosity 
  !
  call ker_noslwa() 
  
  modul = 0 

  call livinf(75_ip,' ',0_ip)
  do iblok = 1,nblok
     call moduls(ITASK_INIUNK)
  end do
  !
  ! No time step is performed
  !
  if( mitim == 0 )  then
     kfl_gotim = 0
     cutim     = 0.0_rp
  end if
  call livinf(76_ip,' ',0_ip)
  !
  ! Start iterating!
  !
  call outfor(91_ip,lun_outpu,' ')

end subroutine Iniunk

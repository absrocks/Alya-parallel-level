subroutine ibm_outerr()
  !------------------------------------------------------------------------
  !****f* Immbou/ibm_outerr
  ! NAME 
  !    ibm_outerr
  ! DESCRIPTION
  !    This routine checks if there are errros and warnings
  ! USES
  ! USED BY
  !    ibm_turnon
  !***
  !------------------------------------------------------------------------
  use def_master
  use def_kermod
  use def_domain
  use def_immbou
  use mod_outfor, only : outfor
  implicit none
  integer(ip)    :: ierro=0,iwarn=0
  character(20)  :: messa
  character(200) :: wmess
  !
  ! Boussinesq without temperature
  !
  if(  kfl_cofor == -1 .and. ( denme == 0.0_rp .or. visme == 0.0_rp ) ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'SPHERE FORCE: DENSITY AND VISCOSITY SHOULD NOT BE ZERO')     
  end if
  !
  ! ERROR MESSAGE
  !
  call errors(3_ip,ierro,iwarn,' ')

end subroutine ibm_outerr

subroutine ibm_reaous()
  !------------------------------------------------------------------------
  !****f* Immbou/ibm_reaous
  ! NAME 
  !    ibm_reaphy
  ! DESCRIPTION
  !    This routine reads the postprocess
  ! USES
  ! USED BY
  !    ibm_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_postpr, only: npp_stepo
  use mod_ecoute, only :  ecoute
  implicit none
  !
  ! Guillaume: Le puse acá porque immbou se ejecuta antes que kermod. Sorry!
  !
  npp_stepo     = -2                                     ! Do not step over the defined postprocesing steps
  call reaous()

end subroutine ibm_reaous

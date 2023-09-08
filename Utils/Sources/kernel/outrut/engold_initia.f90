subroutine engold_initia()
  !-----------------------------------------------------------------------
  !****f* outrut/engold_initia
  ! NAME 
  !    engold_initia
  ! DESCRIPTION
  !    This routine initializes Ensight
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use mod_iofile
  use def_postpr
  implicit none
  !
  ! Default values
  !     
  nppva_ens            = 0
  npart_pos            = 1
  parts_pos(1)%name    = 'Volume Mesh'
  parts_pos(1)%numepart = 1

end subroutine engold_initia

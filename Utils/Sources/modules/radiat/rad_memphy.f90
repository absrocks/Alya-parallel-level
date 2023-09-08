subroutine rad_memphy(itask)
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_memphy
  ! NAME
  !    rad_memphy
  ! DESCRIPTION
  !    Allocate memory for the physical problem
  ! OUTPUT 
  ! USES
  ! USED BY
  !    rad_reaphy
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_radiat
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)

  case(1)
     !
     ! Properties
     !
     allocate(scatt_rad(nspec_rad),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'scatt_rad','rad_memphy',scatt_rad)
     allocate(absor_rad(nspec_rad),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'absor_rad','rad_memphy',absor_rad)
     allocate(aniso_rad(nspec_rad),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'aniso_rad','rad_memphy',aniso_rad)

  case(2)
     !
     ! Interpolation of properties
     !
!!$     allocate(coefk_rad(2,mkint_rad,nmate),stat=istat)
!!$     call memchk(zero,istat,mem_modul(1:2,modul),'COEFK_TEM','rad_memphy',coefk_rad)
!!$     allocate(coefc_rad(2,mcint_rad,nmate),stat=istat)
!!$     call memchk(zero,istat,mem_modul(1:2,modul),'COEFC_TEM','rad_memphy',coefc_rad)
 
  end select

end subroutine rad_memphy

subroutine tem_chemic(&
     ielem,pgaus,gprhs)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_chemic
  ! NAME
  !   tem_radiat
  ! DESCRIPTION
  !    Couple to the heat source from chemical processes
  ! USES
  ! USED BY
  !    tem_elmop2 
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  div_enthalpy_transport,chemical_heat,kfl_coupl, &
                                ID_TEMPER,ID_CHEMIC,radiative_heat
  use def_domain, only       :  mgaus
  implicit none 
  integer(ip), intent(in)    :: ielem,pgaus
  real(rp),    intent(inout) :: gprhs(pgaus)
  integer(ip)                :: igaus

  !
  ! Coupling with CHEMIC (it comes in Gauss points !!!!)
  !
  if (kfl_coupl(ID_TEMPER,ID_CHEMIC) >= 1 ) then
     ! Add to RHS
     do igaus=1,pgaus
        gprhs(igaus)= gprhs(igaus)+div_enthalpy_transport(ielem)%a(igaus)+ &
                      chemical_heat(ielem)%a(igaus)+radiative_heat(ielem)%a(igaus)
     end do
  endif

end subroutine tem_chemic

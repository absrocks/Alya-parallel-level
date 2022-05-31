subroutine rubclo(ngaus,posgp,weigp,ierro)
  !-----------------------------------------------------------------------
  !****f* Domain/rubclo
  ! NAME
  !    rubclo
  ! DESCRIPTION
  !     This routine sets up the integration constants of closed rules in
  !     one dimension.
  ! OUTPUT
  !    
  ! USED BY
  !    rulepw
  !***
  !-----------------------------------------------------------------------
  use  def_kintyp, only    :  ip,rp
  implicit none
  integer(ip), intent(in)  :: ngaus
  integer(ip), intent(out) :: ierro
  real(rp),    intent(out) :: posgp(ngaus),weigp(ngaus)

  ierro=0

  if(ngaus==2) then
     posgp(1)=-1.0_rp
     posgp(2)= 1.0_rp
     weigp(1)= 1.0_rp
     weigp(2)= 1.0_rp

  else if(ngaus==3) then
     posgp(1)=-1.0_rp
     posgp(2)= 0.0_rp
     posgp(3)= 1.0_rp
     weigp(1)= 1.0_rp/3.0_rp
     weigp(2)= 4.0_rp/3.0_rp
     weigp(3)= 1.0_rp/3.0_rp

  else if(ngaus==4) then
     posgp(1)=-1.0_rp
     posgp(2)=-1.0_rp/3.0_rp
     posgp(3)= 1.0_rp/3.0_rp
     posgp(4)= 1.0_rp
     weigp(1)= 1.0_rp/4.0_rp
     weigp(2)= 3.0_rp/4.0_rp
     weigp(3)= 3.0_rp/4.0_rp
     weigp(4)= 1.0_rp/4.0_rp

  else

     ierro=1

  end if

end subroutine rubclo

subroutine linint(ncoef,value,xcoef,varia,dvari)

  !-----------------------------------------------------------------------
  !****f* mathru/linint
  ! NAME 
  !    linint
  ! DESCRIPTION
  !    Compute a variable and its derivatives using linear interpolation
  ! USES
  ! USED BY
  !    *
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: ncoef
  real(rp),    intent(in)  :: value,xcoef(2,ncoef)
  real(rp),    intent(out) :: varia,dvari
  integer(ip)              :: icoef,jcoef
  real(rp)                 :: xvalu

  if(value<xcoef(1,1)) then
     xvalu=xcoef(1,1)
  else if(value>xcoef(1,ncoef)) then
     xvalu=xcoef(1,ncoef)
  else
     xvalu=value
  end if

  icoef=1
  do while(icoef<ncoef)
     icoef=icoef+1
     if(xcoef(1,icoef)>=xvalu) then
        jcoef=icoef-1
        icoef=ncoef
     end if
  end do

  dvari=(xcoef(2,jcoef+1)-xcoef(2,jcoef))/(xcoef(1,jcoef+1)-xcoef(1,jcoef))
  varia=dvari*(xvalu-xcoef(1,jcoef))+xcoef(2,jcoef)
  
end subroutine linint

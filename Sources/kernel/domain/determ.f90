subroutine determ(ndime,pnode,elcod,deriv,xjacm,gpcar,gpdet)

  !-----------------------------------------------------------------------
  !
  ! This routine evaluates the Determan
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime,pnode
  real(rp),    intent(in)  :: elcod(ndime,pnode),deriv(ndime,pnode)
  real(rp),    intent(out) :: xjacm(ndime,ndime)
  real(rp),    intent(out) :: gpcar(ndime,pnode),gpdet
  integer(ip)              :: k
  real(rp)                 :: t1,t2,t3

  if( ndime == 2 .and. pnode == 3 ) then
     !
     ! 2D P1 element
     !
     gpdet = (-elcod(1,1)+elcod(1,2))*(-elcod(2,1)+elcod(2,3)) &
          & -(-elcod(2,1)+elcod(2,2))*(-elcod(1,1)+elcod(1,3))

  else if( ndime == 3 .and. pnode == 4 ) then
     !
     ! 3D P1 element
     !
     gpcar(1,1) =  elcod(1,2) - elcod(1,1)
     gpcar(1,2) =  elcod(1,3) - elcod(1,1)
     gpcar(1,3) =  elcod(1,4) - elcod(1,1)
     gpcar(2,1) =  elcod(2,2) - elcod(2,1)
     gpcar(2,2) =  elcod(2,3) - elcod(2,1)
     gpcar(2,3) =  elcod(2,4) - elcod(2,1)
     gpcar(3,1) =  elcod(3,2) - elcod(3,1)
     gpcar(3,2) =  elcod(3,3) - elcod(3,1)
     gpcar(3,3) =  elcod(3,4) - elcod(3,1)
     t1         =  gpcar(2,2) * gpcar(3,3) - gpcar(3,2) * gpcar(2,3)
     t2         = -gpcar(2,1) * gpcar(3,3) + gpcar(3,1) * gpcar(2,3)
     t3         =  gpcar(2,1) * gpcar(3,2) - gpcar(3,1) * gpcar(2,2)
     gpdet      =  gpcar(1,1) * t1 + gpcar(1,2) * t2 + gpcar(1,3) * t3

  else if ( ndime == 1 ) then
     !
     ! 1D
     !
     xjacm(1,1) = 0.0_rp
     do k = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
     end do
     gpdet = xjacm(1,1)

  else if ( ndime == 2 ) then
     !
     ! 2D
     !
     xjacm(1,1) = 0.0_rp
     xjacm(1,2) = 0.0_rp
     xjacm(2,1) = 0.0_rp
     xjacm(2,2) = 0.0_rp
     do k = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
        xjacm(1,2) = xjacm(1,2) + elcod(1,k) * deriv(2,k)
        xjacm(2,1) = xjacm(2,1) + elcod(2,k) * deriv(1,k)
        xjacm(2,2) = xjacm(2,2) + elcod(2,k) * deriv(2,k)
     end do

     gpdet = xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)

  else if ( ndime == 3 ) then
     !
     ! 3D
     !
     xjacm(1,1) = 0.0_rp ! xjacm = elcod * deriv^t
     xjacm(1,2) = 0.0_rp ! xjaci = xjacm^-1
     xjacm(1,3) = 0.0_rp ! gpcar = xjaci^t * deriv 
     xjacm(2,1) = 0.0_rp
     xjacm(2,2) = 0.0_rp
     xjacm(2,3) = 0.0_rp
     xjacm(3,1) = 0.0_rp
     xjacm(3,2) = 0.0_rp
     xjacm(3,3) = 0.0_rp
     do k = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
        xjacm(1,2) = xjacm(1,2) + elcod(1,k) * deriv(2,k)
        xjacm(1,3) = xjacm(1,3) + elcod(1,k) * deriv(3,k)
        xjacm(2,1) = xjacm(2,1) + elcod(2,k) * deriv(1,k)
        xjacm(2,2) = xjacm(2,2) + elcod(2,k) * deriv(2,k)
        xjacm(2,3) = xjacm(2,3) + elcod(2,k) * deriv(3,k)
        xjacm(3,1) = xjacm(3,1) + elcod(3,k) * deriv(1,k)
        xjacm(3,2) = xjacm(3,2) + elcod(3,k) * deriv(2,k)
        xjacm(3,3) = xjacm(3,3) + elcod(3,k) * deriv(3,k)
     end do

     t1    =  xjacm(2,2) * xjacm(3,3) - xjacm(3,2) * xjacm(2,3)
     t2    = -xjacm(2,1) * xjacm(3,3) + xjacm(3,1) * xjacm(2,3)
     t3    =  xjacm(2,1) * xjacm(3,2) - xjacm(3,1) * xjacm(2,2)
     gpdet =  xjacm(1,1) * t1 + xjacm(1,2) * t2 + xjacm(1,3) * t3

  end if

end subroutine determ

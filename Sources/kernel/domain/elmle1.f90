subroutine elmle1(&
     ndime,nnode,dercg,tragl,elcod,hnatu,hleng)
  !-----------------------------------------------------------------------
  !****f* Domain/elmlen
  ! NAME
  !    elmlen
  ! DESCRIPTION
  !    Compute element lengths
  ! OUTPUT
  !    HLENG ... Element length but not sorted as in elmlen
  ! USED BY
  !    
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     : ip,rp
  implicit none 
  integer(ip), intent(in)  :: ndime,nnode
  real(rp),    intent(in)  :: tragl(ndime,ndime)
  real(rp),    intent(in)  :: hnatu,elcod(ndime,nnode)
  real(rp),    intent(in)  :: dercg(ndime,nnode)
  real(rp),    intent(out) :: hleng(ndime)
  integer(ip)              :: idime,jdime
  real(rp)                 :: enor0,detjm,xjacm(ndime,ndime)
  !
  ! Jacobian
  !
  call mbmabt(xjacm,elcod,dercg,ndime,ndime,nnode) ! J^t
  call invmtx(xjacm,tragl,detjm,ndime)             ! J^(-t)
  !
  ! Element length
  !
  do idime=1,ndime
     enor0=0.0_rp
     do jdime=1,ndime
        enor0=enor0+tragl(idime,jdime)*tragl(idime,jdime)
     end do
     enor0=sqrt(enor0)
     hleng(idime)=hnatu/enor0
  end do

end subroutine elmle1

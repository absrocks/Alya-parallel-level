subroutine d2sdx2(igaus,pnode,pgaus,xjaci,elcod,deriv,wmat2,d2sdx)
  !-----------------------------------------------------------------------
  !****f* domain/d2sdx2
  ! NAME
  !    d2sdx2
  ! DESCRIPTION
  !    This routine calculates D2SDX and WMAT2
  ! USES
  ! USED BY
  !    elmcar
  ! SOURCE
  !-----------------------------------------------------------------------
  use def_domain, only     :  ndime
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: igaus,pnode,pgaus
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(in)  :: xjaci(ndime,ndime)
  real(rp),    intent(in)  :: deriv(ndime,pnode,pgaus)
  real(rp),    intent(out) :: wmat2(ndime,ndime,pnode)
  real(rp),    intent(out) :: d2sdx(ndime,ndime,ndime)
  integer(ip)              :: kdime,idime,jdime,ldime,inode
  !
  ! Obtains (d^2 s_k / d x_i d x_j) as the solution of the system
  ! (d x_l / d s_k) (d^2 s_k / d x_i d x_j) 
  !     = - (d^2 x_l / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j), 
  ! for l,i,j = 1,...,NDIME
  !
  do kdime=1,ndime
     do idime=1,ndime
        do jdime=1,ndime
           d2sdx(kdime,idime,jdime)=0.0_rp
           do ldime=1,ndime
              do inode=1,pnode
                 d2sdx(kdime,idime,jdime)=d2sdx(kdime,idime,jdime)&
                      -xjaci(kdime,ldime)*wmat2(idime,jdime,inode)&
                      *elcod(ldime,inode)
              end do
           end do
        end do
     end do
  end do
  !
  ! Computes the second Cartesian derivatives of the shape functions
  !
  do inode=1,pnode
     do idime=1,ndime
        do jdime=1,ndime
           do kdime=1,ndime
              wmat2(idime,jdime,inode)=wmat2(idime,jdime,inode)&
                   +deriv(kdime,inode,igaus)*d2sdx(kdime,idime,jdime)
           end do
        end do
     end do
  end do

end subroutine d2sdx2

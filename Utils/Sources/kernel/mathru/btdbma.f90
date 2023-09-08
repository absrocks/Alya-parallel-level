subroutine btdbma(aglob,aloca,bmatr,n1,n2)
!------------------------------------------------------------------------
!                                      
!    This routine computes Ag = Bt Al B  when Ag and Al are stored 
!    as full matrices (Ag := aglob, Al := aloca, B := bmatr). The di-
!    mensions are Al -> Mat(n1,n1), Ag -> Mat(n2,n2), B -> Mat(n2,n1) 
!
!------------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: n1,n2
  real(rp),    intent(in)  :: aloca(n1,n1), bmatr(n1,n2)
  real(rp),    intent(out) :: aglob(n2,n2)
  integer(ip)              :: i,j,k,l

  do i=1,n2
     do j=1,n2
        aglob(i,j)=0.0_rp
        do k=1,n1
           do l=1,n1
              aglob(i,j)=aglob(i,j)+bmatr(k,i)*aloca(k,l)*bmatr(l,j)
           end do
        end do
     end do
  end do
  
end subroutine btdbma

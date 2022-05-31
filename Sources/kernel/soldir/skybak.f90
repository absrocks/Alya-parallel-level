subroutine skybak(gstdi,gstlo,gstup,lpont,neqns,rhsid,lures,outso)
!-----------------------------------------------------------------------
!
! This routine performs the back-solution of symmetric equations
! stored in profile form .
! coefficient matrix must be decomposed into its triangular
! factors using skylin before using skybak.
! 
! Input parameters
! 
!     gstdi(neqns)      : reciprocal of diagonal of triangular 
!                           factor
!     gstlo(lpont(neqns): lower triangular factor of matrix
!     gstup(lpont(neqns): upper triangular factor of matrix
!                           (gstup and gstlo have same calling 
!                           address for symmetric matrices)
!     rhsid(neqns)      : right hand side vector in equations
!     lpont(neqns)      : pointer array to bottom of columns 
!                         of gstlo and gstup
!     neqns             : number of equations to be solved
! 
! Output parameter
! 
!     rhsid(neqns)      : solution of equations
!
!-----------------------------------------------------------------------
  use def_parame
  implicit none
  integer(ip) :: lures,neqns
  integer(ip) :: lpont(neqns)
  real(rp)    :: gstlo(lpont(neqns)), gstup(lpont(neqns))
  real(rp)    :: gstdi(neqns), rhsid(neqns)
  integer(ip) :: outso
  integer(ip) :: keqns,ieqns,jeqns,jrows,jheig,kterm,iterm
  real(rp)    :: zeros
  real(rp)    :: fvecdo

  zeros=0.0_rp

!
! Find the first nonzero entry in the right hand side
!
  do keqns=1,neqns
     ieqns=keqns
     if(rhsid(ieqns)/=zeros) go to 200
  end do
  if((neqns>0).and.(outso==1)) write(lures,2000)
  return
200 if(ieqns<neqns) then
!
! Reduce the right hand side
!
     do jeqns=ieqns+1,neqns
        jrows=lpont(jeqns-1)
        jheig=lpont(jeqns)-jrows
        if(jheig>0) then
           rhsid(jeqns)=rhsid(jeqns)-fvecdo(jheig,gstlo(jrows+1),rhsid(jeqns-jheig))
           !rhsid(jeqns)=rhsid(jeqns)-dot_product(&
           !     gstlo(jrows+1:jrows+1+jheig-1),&
           !     rhsid(jeqns-jheig:jeqns-jheig+jheig-1))
        end if
     end do
  end if
!
! Multiply by inverse of diagonal elements
!
  do jeqns=ieqns,neqns
     rhsid(jeqns)=rhsid(jeqns)*gstdi(jeqns)
  end do
!
! Backsubstitution
!
  if(neqns>1) then
     do jeqns=neqns,2,-1
        jrows=lpont(jeqns-1)
        jheig=lpont(jeqns)-jrows
        if(jheig>0) then
           kterm=jeqns-jheig
           do iterm=1,jheig
              rhsid(kterm)=rhsid(kterm)-gstup(jrows+iterm)*rhsid(jeqns)
              kterm=kterm+1
           end do
        end if
     end do
  end if
 
2000 format(11x,'*** WARNING: ZERO RIGHT-HAND-SIDE VECTOR')

end subroutine skybak

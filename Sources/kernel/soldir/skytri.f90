subroutine skytri(gstdi,gstlo,gstup,lpont,neqns,lures,outso)
!-----------------------------------------------------------------------
!
! This routine performs the triangular decomposition 
! of a matrix stored in profile form
! 
! Input parameters:
! 
! gstlo(lpont(neqns)):  lower triangular part of matrix
! gstup(lpont(neqns)):  upper part of triangular matrix
! gstdi(neqns)       :  diagonals of triangular matrix
! lpont(neqns)       :  pointers to bottom of colums of 
!                       gstlo and gstup arrays
! neqns              :  number of equations to be solved
! 
! Output parameters:
! 
! gstlo(lpont(neqns)):  lower triangular factor of matrix
! gstup(lpont(neqns)):  upper triangular factor of matrix
! gstdi(neqns)       :  inverse of diagonal matrix in triangular 
!                       factor
! 
! Local  parameters:
! 
! tolle              :  tollerance used to check for null pivot
!                       it should be set to approximate half-word
!                       precision 
! 
!-----------------------------------------------------------------------
  use def_kintyp
  implicit none
  integer(ip) :: neqns,lures
  integer(ip) :: lpont(neqns)
  real(rp)    :: gstlo(lpont(neqns))
  real(rp)    :: gstup(lpont(neqns))
  real(rp)    :: gstdi(neqns) 
  integer(ip) :: outso
  integer(ip) :: nillc,kpivo,idiag,ieqns,irows,iheig
  integer(ip) :: istar,iends,jdiag,jheig,irhei,jeqns
  integer(ip) :: jdhei,ihei1
  real(rp)    :: zeros,unity,tolle,pivot,saval
  real(rp)    :: fvecdo

  nillc=0                                               ! Number of ILL-C
  zeros=0.0_rp
  unity=1.0_rp
  tolle=0.5d-07
!
! Loop through the columns to perform the triangular decomposition
!
  kpivo=0
  idiag=1
  do ieqns=1,neqns
     irows=idiag+1
     idiag=lpont(ieqns)
     iheig=idiag-irows
     if(iheig>0) then
        istar=ieqns-iheig
        iends=ieqns-1
!
! If diagonal is zero compute a norm for singularity test
!
        if(gstdi(ieqns)==zeros) call skycek(gstup(irows:),iheig,saval)
        do jeqns=istar,iends
           irows=irows+1
           jdiag=lpont(jeqns)
           jheig=min(jdiag-lpont(jeqns-1),jeqns-istar+1_ip)
           if(jheig>0) then
              irhei=irows-jheig
              jdhei=jdiag-jheig+1
              gstup(irows)=gstup(irows)-fvecdo(jheig,gstup(irhei),gstlo(jdhei)) 
              gstlo(irows)=gstlo(irows)-fvecdo(jheig,gstlo(irhei),gstup(jdhei)) 
              !gstup(irows)=gstup(irows)-dot_product(&
              !     gstup(irhei:irhei+jheig-1),&
              !     gstlo(jdhei:jdhei+jheig-1)) 
              !gstlo(irows)=gstlo(irows)-dot_product(&
              !     gstlo(irhei:irhei+jheig-1),&
              !     gstup(jdhei:jdhei+jheig-1)) 
           end if
        end do
     end if
!
! Reduce the diagonal
!
     if(iheig>=0) then
        pivot=gstdi(ieqns)
        irows=idiag-iheig
        irhei=ieqns-iheig-1
        ihei1=iheig+1
        call skydia(gstlo(irows:irows+iheig),gstup(irows:irows+iheig),&
             gstdi(irhei:irhei+iheig),ihei1,gstdi(ieqns))
!
! Check for possible errors and print warnings
!
        if(pivot<zeros) kpivo=kpivo+1
        if(abs(gstdi(ieqns))<tolle*abs(pivot)) nillc=nillc+1
        if((abs(gstdi(ieqns))<1.0d-15).and.(outso==1)) write(lures,2001) ieqns
!
! Complete rank test for a zero diagonal case
!
        if(pivot==zeros.and.iheig>0) then
           if(abs(gstdi(ieqns))<tolle*saval) then
              write(lures,2003) ieqns
              stop
           end if
        end if
     end if
!
! Store reciprocal of diagonal
!
     if(gstdi(ieqns)/=zeros) gstdi(ieqns)=unity/gstdi(ieqns)
  end do
  
  if((nillc/=0).and.(outso==1)) write(lures,2000) nillc
!
! Formats
!
2000 format(11x,'*** ILL-CONDITIONING.',&
          ' LOSS OF AT LEAST 7 DIGITS FOR ',i12,' EQS.')
2001 format(11x,'*** SINGULAR MATRIX.',&
          ' REDUCED DIAGONAL IS ZERO FOR EQ.',I12)
2003 format(11x,'*** SINGULAR MATRIX.',&
          ' RANK FAILURE FOR ZERO UNREDUCED DIAGONAL IN EQ. ',I12)

end subroutine skytri
